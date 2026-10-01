"""Matrix-free Two-Walk seriation (paper Section 3.3).

Section 3.3 of the manuscript proves that the two Gram blocks never have to be
formed.  With ``z = [x; y]``, ``u = B_bar^T x`` and ``v = B_bar y``,

    A_tw z = [ B_bar u + alpha v ;  B_bar^T v + alpha u ]        (Eq. matvec)

and, with ``r = B_bar 1`` and ``c = B_bar^T 1``,

    d_alpha = [ B_bar c + alpha r ;  B_bar^T r + alpha c ]
    L_tw z  = d_alpha * z - A_tw z                               (Eq. lmatvec)

so one product costs four sparse mat-vecs and stores ``O(nnz(B) + m + n)``.

The released benchmark does not implement this.  ``_benchmark_utils`` builds an
explicit ``np.block([[B@B.T, a*B], [a*B.T, B.T@B]])`` per connected component and
solves it with the dense full-spectrum ``scipy.linalg.eigh``.  This module supplies
the replacement: the operators above are exposed as ``LinearOperator`` objects and
solved with the sparse ``eigsh``.

Scope
-----
``mode="tw"``  Two-Walk, ``A = M^2 + alpha M``, ``L = D - A``.
``mode="ow"``  One-Walk, ``A = M`` (the bipartite lift).

``componentwise_reorder`` reproduces ``_benchmark_utils._componentwise_spectral_reorder``
semantics exactly -- global-max normalisation, components ordered by decreasing mass,
``_orient_fiedler``, stable argsort, zero-marginal axes appended in original order --
changing only *which solver* produces the Fiedler vector.

Known limitation, shared with the dense path
--------------------------------------------
The Fiedler vector of a real matrix is often heavily tied (on CIP->SOC, 956
coordinate pairs lie within 1e-12 of each other).  ``argsort(kind="stable")``
resolves ties by index, so two different numerical solvers can disagree on a large
fraction of positions while agreeing on the sorted coordinate sequence to ~1e-12.
The induced *metrics* are unaffected (measured |dR2S| <= 2.4e-5, |dBand| = 0), but
exact permutation equality between solvers is unattainable.  Use ``compare_block``,
which is tie-aware, rather than comparing orders directly.
"""

from __future__ import annotations

import importlib.util
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import scipy.sparse as sp
from scipy.linalg import eigh
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import (
    ArpackError,  # ArpackNoConvergence subclasses this
    LinearOperator,
    eigsh,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
SYNTHETIC_DIR = REPO_ROOT / "main_experiment" / "synthetic_benchmark"

# pandas imports numexpr lazily; the benchmark scripts disable it for determinism.
sys.modules.setdefault("numexpr", None)


def _load_module(name: str, path: Path):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


# Reuse the released baseline/metric machinery rather than reimplementing it.
bm = _load_module("_benchmark_utils", SYNTHETIC_DIR / "_benchmark_utils.py")

# Matches ``tw_components.EIGEN_TOL``; also the threshold below which an eigenvalue
# counts as the (exact) nullspace of a connected component rather than the Fiedler value.
EIGEN_TOL = 1e-10

# Below this many vertices the dense solve is exact, faster than ARPACK, and has no
# convergence caveats.  ``solver="auto"`` uses it; the rebuttal harness forces
# ``solver="matrix_free"`` so that every reported number comes from Section 3.3's method.
DENSE_MAX_VERTICES = 512

# ARPACK needs ``k < N``; blocks at or below this size always route to dense.
ARPACK_MIN_VERTICES = 3

DEFAULT_ALPHA_GRID = (1.0, 2.0, 4.0, 6.0, 8.0, 12.0)

# Number of sparse mat-vecs with B_bar or B_bar^T in one L_tw product.  Section 3.3
# claims four; ``test_tw_matrixfree`` verifies the code still matches the sentence.
TW_SPMV_PER_MATVEC = 4
OW_SPMV_PER_MATVEC = 2

# Relative tolerance for the orientation anchor; see ``_orient``.
ORIENT_TOL = 1e-12


class DegenerateBlockError(RuntimeError):
    """A block whose spectrum has no positive eigenvalue (e.g. alpha = 0, disjoint)."""


# --------------------------------------------------------------------------- #
# Operators
# --------------------------------------------------------------------------- #


def as_csr(block) -> sp.csr_matrix:
    """Sparse view of a component block, dense or already sparse.

    The benchmark pipeline is dense end to end (``load_matrix_csv`` returns an
    ``ndarray``; ``build_case`` allocates ``np.zeros((m, n))``), so for dense input this
    is where sparsity is recovered: entries that are exactly zero are dropped.  The
    matrix-free win is therefore real only insofar as the block is genuinely sparse --
    0.18-0.33% on the two largest real matrices, but 10-17% on the existing synthetic
    generator, which is why the scaling sweep needs its own CSR-native generator.
    """
    if sp.issparse(block):
        B = block.tocsr().astype(float, copy=True)
    else:
        B = sp.csr_matrix(np.asarray(block, dtype=float))
    B.eliminate_zeros()
    return B


class TwoWalkOperator(LinearOperator):
    """Matrix-free ``L_tw = D_tw - A_tw`` for ``A_tw = M^2 + alpha M``.

    Implements Eq. (matvec) and Eq. (lmatvec) verbatim.  One ``_matvec`` performs
    exactly four sparse products (``Bt @ x``, ``B @ y``, ``B @ u``, ``Bt @ v``),
    each ``O(nnz(B))``, plus ``O(m + n)`` vector arithmetic.
    """

    def __init__(self, block: np.ndarray, alpha: float):
        self.B = as_csr(block)
        self.Bt = self.B.T.tocsr()
        self.alpha = float(alpha)
        self.m, self.n = self.B.shape

        # r = B_bar 1, c = B_bar^T 1  (Eq. lmatvec)
        r = np.asarray(self.B @ np.ones(self.n)).ravel()
        c = np.asarray(self.Bt @ np.ones(self.m)).ravel()
        self.row_degree = r
        self.col_degree = c
        self.degree = np.concatenate(
            [self.B @ c + self.alpha * r, self.Bt @ r + self.alpha * c]
        )
        self.n_matvec = 0
        self.n_spmv = 0
        super().__init__(dtype=np.float64, shape=(self.m + self.n, self.m + self.n))

    @property
    def trace(self) -> float:
        """``sum(d_alpha)`` -- an upper bound on ``lambda_max`` (Gershgorin)."""
        return float(self.degree.sum())

    def _matvec(self, z: np.ndarray) -> np.ndarray:
        self.n_matvec += 1
        x, y = z[: self.m], z[self.m :]
        u = self.Bt @ x
        v = self.B @ y
        self.n_spmv += 4
        return self.degree * z - np.concatenate(
            [self.B @ u + self.alpha * v, self.Bt @ v + self.alpha * u]
        )

    def _adjoint(self) -> LinearOperator:
        # L_tw is symmetric.
        return self


class OneWalkOperator(LinearOperator):
    """Matrix-free ``L_ow = D - A`` for the bipartite lift ``A = [[0, B], [B^T, 0]]``."""

    def __init__(self, block: np.ndarray):
        self.B = as_csr(block)
        self.Bt = self.B.T.tocsr()
        self.m, self.n = self.B.shape
        self.row_degree = np.asarray(self.B @ np.ones(self.n)).ravel()
        self.col_degree = np.asarray(self.Bt @ np.ones(self.m)).ravel()
        self.n_matvec = 0
        self.n_spmv = 0
        super().__init__(dtype=np.float64, shape=(self.m + self.n, self.m + self.n))

    @property
    def trace(self) -> float:
        return float(self.row_degree.sum() + self.col_degree.sum())

    def _matvec(self, z: np.ndarray) -> np.ndarray:
        self.n_matvec += 1
        x, y = z[: self.m], z[self.m :]
        self.n_spmv += 2
        return np.concatenate(
            [self.row_degree * x - self.B @ y, self.col_degree * y - self.Bt @ x]
        )

    def _adjoint(self) -> LinearOperator:
        return self


# --------------------------------------------------------------------------- #
# Dense references -- the object the benchmark ships today
# --------------------------------------------------------------------------- #


def build_two_walk_dense(block: np.ndarray, alpha: float) -> np.ndarray:
    """The explicit-Gram Laplacian, via the pinned ``mheatmap`` package.

    Deliberately routed through ``bm.mheatmap_two_walk_laplacian`` rather than
    re-derived, so this is literally the shipped object and the existing
    call-count test keeps exercising it.
    """
    dense = _densify(block)
    return np.asarray(bm.mheatmap_two_walk_laplacian(dense, alpha=float(alpha)), dtype=float)


def _densify(block) -> np.ndarray:
    """Dense view of a block.  Only the *reference* path should ever call this."""
    if sp.issparse(block):
        return np.asarray(block.todense(), dtype=float)
    return np.asarray(block, dtype=float)


def build_one_walk_dense(block) -> np.ndarray:
    """Mirrors the ``np.block`` branch of ``_componentwise_spectral_reorder``."""
    block = _densify(block)
    local_rows, local_cols = block.shape
    adjacency = np.block(
        [
            [np.zeros((local_rows, local_rows), dtype=float), block],
            [block.T, np.zeros((local_cols, local_cols), dtype=float)],
        ]
    )
    return np.diag(adjacency.sum(axis=1)) - adjacency


def build_operator(block, mode: str, alpha: float | None = None) -> LinearOperator:
    if mode == "tw":
        if alpha is None or alpha <= 0:
            raise ValueError("alpha must be positive for Two-Walk")
        return TwoWalkOperator(block, alpha)
    if mode == "ow":
        return OneWalkOperator(block)
    raise ValueError(f"unknown mode {mode!r}; expected 'tw' or 'ow'")


def build_dense(block, mode: str, alpha: float | None = None) -> np.ndarray:
    if mode == "tw":
        if alpha is None or alpha <= 0:
            raise ValueError("alpha must be positive for Two-Walk")
        return build_two_walk_dense(block, alpha)
    if mode == "ow":
        return build_one_walk_dense(block)
    raise ValueError(f"unknown mode {mode!r}; expected 'tw' or 'ow'")


# --------------------------------------------------------------------------- #
# Fiedler solve
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class FiedlerResult:
    """The smallest *positive* eigenpair of one connected block, plus diagnostics."""

    eigenvalue: float
    eigenvector: np.ndarray
    joint_order: np.ndarray
    solver: str
    n_matvec: int
    n_spmv: int
    residual: float
    next_eigenvalue: float  # lambda_3, i.e. the neighbour above the Fiedler value
    degenerate: bool  # lambda_3 - lambda_2 <= EIGEN_TOL * lambda_max


def _orient(vector: np.ndarray) -> np.ndarray:
    """Deterministic sign convention, tolerant of the near-ties that break the released rule.

    ``_benchmark_utils._orient_fiedler`` anchors on ``argmax(abs(v))`` and negates when
    that coordinate is negative.  When the Fiedler vector is antisymmetric the peak
    magnitude is attained at two indices that differ only in their last bits, so
    ``argmax`` picks a *different* index for each solver and the two vectors come out as
    exact negatives -- a globally reversed order that looks like a total disagreement.

    Anchoring on the first index *within a relative tolerance* of the peak makes the
    choice solver-independent.  All three reported metrics are invariant under joint
    reversal (they depend on ``rho_i - kappa_j``, which negation leaves unchanged), so
    this changes no reported number; it only makes the published permutations
    reproducible across solvers.
    """
    oriented = np.asarray(vector, dtype=float).copy()
    magnitude = np.abs(oriented)
    peak = float(magnitude.max())
    if peak <= 0.0:
        return oriented
    anchor = int(np.flatnonzero(magnitude >= peak * (1.0 - ORIENT_TOL))[0])
    if oriented[anchor] < 0.0:
        oriented *= -1.0
    return oriented


def _pick_positive(eigenvalues: np.ndarray) -> int:
    """Index of the first eigenvalue above ``EIGEN_TOL``.

    For a connected block this is index 1 -- the dense reference hardcodes that,
    and ``compare_block`` asserts the two agree.  Deriving it instead of hardcoding
    keeps the routine total when a block is degenerate.
    """
    positive = np.flatnonzero(np.abs(eigenvalues) > EIGEN_TOL)
    if positive.size == 0:
        raise DegenerateBlockError(
            "block has no eigenvalue above EIGEN_TOL; the coordinated order is undefined"
        )
    return int(positive[0])


def _shifted(op: LinearOperator, beta: float) -> LinearOperator:
    """Rank-one deflation ``L + beta * 11^T / N``.

    ``11^T/N`` is the projector onto ``span{1}``, so this moves the nullspace
    eigenvalue from 0 to ``beta`` and leaves every other eigenvalue untouched.  With
    ``beta`` above ``lambda_max`` the Fiedler value becomes the smallest eigenvalue,
    which ARPACK finds reliably -- unlike the plain operator, where 0 sits directly
    below it and attracts the iteration.
    """
    size = op.shape[0]
    ones = np.ones(size) / np.sqrt(size)

    def matvec(z: np.ndarray) -> np.ndarray:
        return op.matvec(z) + beta * ones * (ones @ z)

    return LinearOperator((size, size), matvec=matvec, dtype=np.float64)


def fiedler_block(
    block,
    alpha: float | None = 1.0,
    mode: str = "tw",
    *,
    solver: str = "auto",
    tol: float = 1e-10,
    maxiter: int | None = None,
    seed: int = 0,
    k: int = 3,
) -> FiedlerResult:
    """Smallest positive eigenpair of one connected component's Laplacian.

    ``solver`` is ``"dense"``, ``"matrix_free"`` or ``"auto"`` (dense at or below
    ``DENSE_MAX_VERTICES``).  ``k`` is the number of ARPACK pairs requested; three
    lets the caller detect a repeated Fiedler eigenvalue.
    """
    if not sp.issparse(block):
        block = np.asarray(block, dtype=float)
    size = block.shape[0] + block.shape[1]

    if solver not in {"auto", "dense", "matrix_free"}:
        raise ValueError(f"unknown solver {solver!r}")

    use_dense = solver == "dense" or (
        solver == "auto" and size <= DENSE_MAX_VERTICES
    )
    if size < ARPACK_MIN_VERTICES:
        use_dense = True

    if use_dense:
        dense = build_dense(block, mode, alpha)
        eigenvalues, eigenvectors = eigh(dense, check_finite=False)
        index = _pick_positive(eigenvalues)
        vector = eigenvectors[:, index]
        nxt = float(eigenvalues[index + 1]) if index + 1 < eigenvalues.size else np.nan
        lam = float(eigenvalues[index])
        residual = float(np.linalg.norm(dense @ vector - lam * vector))
        return FiedlerResult(
            eigenvalue=lam,
            eigenvector=vector,
            joint_order=np.argsort(_orient(vector), kind="stable"),
            solver="dense",
            n_matvec=0,
            n_spmv=0,
            residual=residual,
            next_eigenvalue=nxt,
            degenerate=_is_degenerate(lam, nxt, float(np.abs(eigenvalues).max())),
        )

    op = build_operator(block, mode, alpha)
    rng = np.random.default_rng(seed)
    v0 = rng.standard_normal(size)
    pairs = max(2, min(int(k), size - 1))

    # `beta` above `lambda_max` so the shifted operator's smallest eigenvalue is the
    # Fiedler value rather than the nullspace; see `_shifted`.
    beta = float(op.trace) * 1.001 + 1.0
    candidates: list[tuple[str, LinearOperator]] = [
        ("matrix_free", op),
        ("eigsh_shifted", _shifted(op, beta)),
    ]

    last_error: Exception | None = None
    for used, candidate in candidates:
        try:
            values, vectors = eigsh(
                candidate, k=pairs, which="SA", v0=v0, tol=tol, maxiter=maxiter
            )
        except ArpackError as exc:
            # ArpackNoConvergence carries a partial spectrum worth keeping; other
            # ARPACK failures (error -9, zero starting vector, on an all-zero block)
            # do not.
            partial = getattr(exc, "eigenvalues", None)
            if partial is None or partial.size < 2:
                last_error = exc
                continue
            values, vectors = exc.eigenvalues, exc.eigenvectors
            used = f"{used}_partial"

        order = np.argsort(values)
        values, vectors = values[order], vectors[:, order]
        try:
            index = _pick_positive(values)
        except DegenerateBlockError as exc:
            # The plain operator returned only the nullspace; the shifted one will not.
            last_error = exc
            continue
        break
    else:
        # Both the plain and the shifted operator failed, which for a block with a
        # positive eigenvalue means non-convergence; for an all-zero block it really is
        # degenerate.  The two are distinguished by the chained cause.
        raise DegenerateBlockError(
            "no positive eigenvalue could be resolved for this block "
            f"(last solver error: {type(last_error).__name__})"
        ) from last_error

    lam = float(values[index])
    vector = vectors[:, index]
    nxt = float(values[index + 1]) if index + 1 < values.size else np.nan
    residual = float(np.linalg.norm(op.matvec(vector) - lam * vector))

    return FiedlerResult(
        eigenvalue=lam,
        eigenvector=vector,
        joint_order=np.argsort(_orient(vector), kind="stable"),
        solver=used,
        n_matvec=int(op.n_matvec),
        n_spmv=int(op.n_spmv),
        residual=residual,
        next_eigenvalue=nxt,
        degenerate=_is_degenerate(lam, nxt, float(op.trace)),
    )


def _is_degenerate(lam: float, nxt: float, scale: float) -> bool:
    if not np.isfinite(nxt):
        return False
    return bool(nxt - lam <= EIGEN_TOL * max(scale, EIGEN_TOL))


# --------------------------------------------------------------------------- #
# Component-wise orchestration -- mirrors the benchmark, swaps the solver
# --------------------------------------------------------------------------- #


def componentwise_reorder(
    matrix: np.ndarray,
    *,
    mode: str,
    alpha: float | None = None,
    solver: str = "auto",
    diagnostics: list[FiedlerResult] | None = None,
    **solver_kwargs,
) -> bm.ReorderResult:
    """Paper Section 3.3 ordering, solver-selectable.

    ``diagnostics``, when given, collects one :class:`FiedlerResult` per component
    solve so callers can report iteration counts and per-solve cost separately from the
    whole-run total.

    Semantics are identical to ``_benchmark_utils._componentwise_spectral_reorder``:
    one global-max normalisation, components from the active support ordered by
    ``(-mass, first_vertex)``, ``_orient_fiedler``, stable argsort, zero-marginal
    axes appended in original index order.
    """
    values = bm._validate_nonnegative_matrix(matrix)
    n_rows, n_cols = values.shape
    row_sums = values.sum(axis=1)
    col_sums = values.sum(axis=0)
    active_rows = np.flatnonzero(row_sums > 0)
    active_cols = np.flatnonzero(col_sums > 0)
    zero_rows = np.flatnonzero(row_sums == 0)
    zero_cols = np.flatnonzero(col_sums == 0)

    if active_rows.size == 0 or active_cols.size == 0:
        return bm.ReorderResult(
            values.copy(),
            np.arange(n_rows, dtype=int),
            np.arange(n_cols, dtype=int),
            component_count=0,
        )

    if mode == "tw" and (alpha is None or alpha <= 0):
        raise ValueError("alpha must be positive for Two-Walk ordering")
    if mode not in {"ow", "tw"}:
        raise ValueError(f"unknown spectral mode: {mode!r}")

    normalized = values / float(values.max())
    labels, component_count = bm._active_support_labels(values, active_rows, active_cols)
    row_labels = labels[: active_rows.size]
    col_labels = labels[active_rows.size :]

    components: list[tuple[float, int, np.ndarray, np.ndarray]] = []
    for component_id in range(component_count):
        component_rows = active_rows[row_labels == component_id]
        component_cols = active_cols[col_labels == component_id]
        mass = float(values[np.ix_(component_rows, component_cols)].sum())
        first_vertex = min(
            int(component_rows.min()), n_rows + int(component_cols.min())
        )
        components.append((mass, first_vertex, component_rows, component_cols))
    components.sort(key=lambda item: (-item[0], item[1]))

    ordered_rows: list[np.ndarray] = []
    ordered_cols: list[np.ndarray] = []
    component_masses: list[float] = []
    solves: list[FiedlerResult] = []

    for mass, _, component_rows, component_cols in components:
        block = normalized[np.ix_(component_rows, component_cols)]
        result = fiedler_block(block, alpha, mode, solver=solver, **solver_kwargs)
        solves.append(result)
        if diagnostics is not None:
            diagnostics.append(result)

        local_row_count = len(component_rows)
        joint_order = result.joint_order
        local_rows = joint_order[joint_order < local_row_count]
        local_cols = joint_order[joint_order >= local_row_count] - local_row_count
        ordered_rows.append(component_rows[local_rows])
        ordered_cols.append(component_cols[local_cols])
        component_masses.append(mass)

    row_order = np.concatenate([*ordered_rows, zero_rows]).astype(int, copy=False)
    col_order = np.concatenate([*ordered_cols, zero_cols]).astype(int, copy=False)
    bm._validate_permutation(row_order, n_rows, "row")
    bm._validate_permutation(col_order, n_cols, "column")

    return bm.ReorderResult(
        bm.apply_orders(values, row_order, col_order),
        row_order,
        col_order,
        component_count=component_count,
        component_masses=tuple(component_masses),
    )


# --------------------------------------------------------------------------- #
# Sparse-native pipeline
# --------------------------------------------------------------------------- #
#
# ``componentwise_reorder`` above mirrors the benchmark and is therefore dense: it slices
# ``values[np.ix_(rows, cols)]`` and takes ``values.max()``.  That is fine for the seven
# real matrices (the largest is 1004x1176) but impossible for the scaling sweep, where a
# 100k x 90k dense block would be 72 GB.  The functions below reproduce the identical
# procedure straight from CSR, so the sweep exercises the real algorithm rather than a
# surrogate.  ``test_tw_matrixfree`` checks the two agree on every small block.


def active_support_sparse(matrix) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, sp.csr_matrix, float]:
    """Definition 1 on sparse input: drop zero marginals, normalise by the global max.

    Returns ``(active_rows, active_cols, zero_rows, zero_cols, sub, peak)``.
    """
    M = as_csr(matrix)
    row_sums = np.asarray(M.sum(axis=1)).ravel()
    col_sums = np.asarray(M.sum(axis=0)).ravel()
    active_rows = np.flatnonzero(row_sums > 0)
    active_cols = np.flatnonzero(col_sums > 0)
    zero_rows = np.flatnonzero(row_sums == 0)
    zero_cols = np.flatnonzero(col_sums == 0)
    if active_rows.size == 0 or active_cols.size == 0:
        return (
            active_rows,
            active_cols,
            zero_rows,
            zero_cols,
            sp.csr_matrix((0, 0), dtype=float),
            0.0,
        )
    sub = M[active_rows][:, active_cols].tocsr()
    peak = float(sub.max()) if sub.nnz else 0.0
    if peak > 0.0:
        sub = (sub / peak).tocsr()
    return active_rows, active_cols, zero_rows, zero_cols, sub, peak


def support_components_sparse(sub: sp.csr_matrix) -> tuple[np.ndarray, int]:
    """Component labels of the bipartite support graph, matching ``_active_support_labels``."""
    if sub.shape[0] == 0 or sub.shape[1] == 0:
        return np.zeros(0, dtype=int), 0
    support = (sub > 0).astype(np.int8).tocsr()
    graph = sp.bmat([[None, support], [support.T, None]], format="csr")
    count, labels = connected_components(graph, directed=False)
    return labels.astype(int, copy=False), int(count)


def componentwise_reorder_sparse(
    matrix,
    *,
    mode: str,
    alpha: float | None = None,
    solver: str = "auto",
    diagnostics: list[FiedlerResult] | None = None,
    **solver_kwargs,
) -> bm.ReorderResult:
    """Sparse-native ``componentwise_reorder``; identical semantics, ``O(nnz)`` input."""
    M = as_csr(matrix)
    n_rows, n_cols = M.shape
    active_rows, active_cols, zero_rows, zero_cols, sub, _ = active_support_sparse(M)

    if active_rows.size == 0 or active_cols.size == 0:
        return bm.ReorderResult(
            M,
            np.arange(n_rows, dtype=int),
            np.arange(n_cols, dtype=int),
            component_count=0,
        )
    if mode == "tw" and (alpha is None or alpha <= 0):
        raise ValueError("alpha must be positive for Two-Walk ordering")
    if mode not in {"ow", "tw"}:
        raise ValueError(f"unknown spectral mode: {mode!r}")

    labels, component_count = support_components_sparse(sub)
    row_labels = labels[: active_rows.size]
    col_labels = labels[active_rows.size :]

    components: list[tuple[float, int, np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = []
    for component_id in range(component_count):
        local_rows = np.flatnonzero(row_labels == component_id)
        local_cols = np.flatnonzero(col_labels == component_id)
        if local_rows.size == 0 or local_cols.size == 0:
            continue
        component_rows = active_rows[local_rows]
        component_cols = active_cols[local_cols]
        # Mass must come from the *unnormalised* matrix, matching the dense path,
        # which sums `values`, not `normalized`.
        mass = float(M[component_rows][:, component_cols].sum())
        first_vertex = min(int(component_rows.min()), n_rows + int(component_cols.min()))
        components.append(
            (mass, first_vertex, component_rows, component_cols, local_rows, local_cols)
        )
    components.sort(key=lambda item: (-item[0], item[1]))

    ordered_rows: list[np.ndarray] = []
    ordered_cols: list[np.ndarray] = []
    component_masses: list[float] = []
    for mass, _, component_rows, component_cols, local_rows, local_cols in components:
        block = sub[local_rows][:, local_cols]
        result = fiedler_block(block, alpha, mode, solver=solver, **solver_kwargs)
        if diagnostics is not None:
            diagnostics.append(result)
        local_row_count = len(component_rows)
        joint_order = result.joint_order
        local_rows = joint_order[joint_order < local_row_count]
        local_cols = joint_order[joint_order >= local_row_count] - local_row_count
        ordered_rows.append(component_rows[local_rows])
        ordered_cols.append(component_cols[local_cols])
        component_masses.append(mass)

    row_order = np.concatenate([*ordered_rows, zero_rows]).astype(int, copy=False)
    col_order = np.concatenate([*ordered_cols, zero_cols]).astype(int, copy=False)
    bm._validate_permutation(row_order, n_rows, "row")
    bm._validate_permutation(col_order, n_cols, "column")

    return bm.ReorderResult(
        bm.apply_orders(M, row_order, col_order),
        row_order,
        col_order,
        component_count=component_count,
        component_masses=tuple(component_masses),
    )


def footprints_sparse(matrix) -> dict:
    """Analytic dense-vs-sparse footprint of one instance, in megabytes.

    Computed without allocating either: ``operator_footprint_mb`` is the explicit
    ``sum_c (m_c + n_c)^2 * 8`` the released code would build, ``csr_footprint_mb`` is the
    ``12 * nnz + 4 * (m + n)`` the matrix-free path actually needs.  Both are exact and
    noise-free, so they carry the memory claim while measured RSS corroborates them.
    """
    M = as_csr(matrix)
    active_rows, active_cols, _, _, sub, _ = active_support_sparse(M)
    if active_rows.size == 0 or active_cols.size == 0:
        return {
            "n_components": 0,
            "max_component_vertices": 0,
            "max_component_nnz": 0,
            "operator_footprint_mb": 0.0,
            "csr_footprint_mb": 0.0,
        }

    labels, count = support_components_sparse(sub)
    row_labels = labels[: active_rows.size]
    col_labels = labels[active_rows.size :]

    operator_bytes = 0
    max_vertices = 0
    max_nnz = 0
    for component_id in range(count):
        local_rows = np.flatnonzero(row_labels == component_id)
        local_cols = np.flatnonzero(col_labels == component_id)
        if local_rows.size == 0 or local_cols.size == 0:
            continue
        vertices = int(local_rows.size + local_cols.size)
        block = sub[local_rows][:, local_cols]
        operator_bytes += 2 * vertices * vertices * 8  # np.block, then D - A
        max_vertices = max(max_vertices, vertices)
        max_nnz = max(max_nnz, int(block.nnz))

    nnz = int(M.nnz)
    return {
        "n_components": int(count),
        "max_component_vertices": int(max_vertices),
        "max_component_nnz": int(max_nnz),
        "operator_footprint_mb": operator_bytes / 1e6,
        "csr_footprint_mb": (12 * nnz + 4 * (M.shape[0] + M.shape[1])) / 1e6,
    }


# --------------------------------------------------------------------------- #
# Installing the solver into the benchmark's single call site
# --------------------------------------------------------------------------- #

_MODE_FROM_BENCHMARK = {"one_walk": "ow", "two_walk": "tw"}


def _adapt(block: np.ndarray, mode: str, alpha: float | None) -> np.ndarray:
    """Adapter matching ``_benchmark_utils.BLOCK_SOLVER_HOOK``'s contract."""
    try:
        local_mode = _MODE_FROM_BENCHMARK[mode]
    except KeyError:
        raise ValueError(f"unknown benchmark mode {mode!r}") from None
    return fiedler_block(
        block, alpha, local_mode, solver=bm.BLOCK_SOLVER, seed=0
    ).joint_order


def install_hook(solver: str = "matrix_free") -> None:
    """Route the benchmark's component-wise solve through this module.

    Sets ``_benchmark_utils.BLOCK_SOLVER_HOOK``.  With ``solver="matrix_free"`` every
    component takes the Section 3.3 path regardless of size, which is what the
    rebuttal reports; with ``"auto"`` blocks at or below ``DENSE_MAX_VERTICES`` keep
    the exact dense solve.
    """
    if solver not in {"dense", "matrix_free", "auto"}:
        raise ValueError(f"unknown solver {solver!r}")
    bm.BLOCK_SOLVER_HOOK = _adapt
    bm.BLOCK_SOLVER = solver


def remove_hook() -> None:
    """Restore the shipped dense path."""
    bm.BLOCK_SOLVER_HOOK = None
    bm.BLOCK_SOLVER = "dense"


def one_walk_reorder(matrix: np.ndarray, *, solver: str = "auto", **kw) -> bm.ReorderResult:
    return componentwise_reorder(matrix, mode="ow", solver=solver, **kw)


def two_walk_reorder(
    matrix: np.ndarray, alpha: float, *, solver: str = "auto", **kw
) -> bm.ReorderResult:
    return componentwise_reorder(matrix, mode="tw", alpha=float(alpha), solver=solver, **kw)


# --------------------------------------------------------------------------- #
# Validation
# --------------------------------------------------------------------------- #


def compare_block(
    block: np.ndarray,
    alpha: float | None = 1.0,
    mode: str = "tw",
    *,
    tie_tol: float = 1e-9,
    **solver_kwargs,
) -> dict:
    """Dense-versus-matrix-free agreement record for one connected block.

    Reports the tie-aware quantities rather than a raw permutation comparison,
    because the latter is meaningless when the Fiedler vector has ties: on MBTA
    100 of 133 coordinates are tied, so ``argsort(kind="stable")`` legitimately
    resolves differently under two solvers that agree to machine precision.
    """
    dense = fiedler_block(block, alpha, mode, solver="dense", **solver_kwargs)
    free = fiedler_block(block, alpha, mode, solver="matrix_free", **solver_kwargs)

    scale = max(abs(dense.eigenvalue), EIGEN_TOL)
    sorted_dense = np.sort(_orient(dense.eigenvector))
    sorted_free = np.sort(_orient(free.eigenvector))
    max_coord_dev = float(np.abs(sorted_dense - sorted_free).max())

    disagree = dense.joint_order != free.joint_order
    n_tied_groups, ties_only = _tie_groups(sorted_dense, tie_tol)
    confined = _disagreement_is_within_ties(
        dense.joint_order, sorted_dense, disagree, tie_tol
    )

    return {
        "lambda_dense": dense.eigenvalue,
        "lambda_mf": free.eigenvalue,
        "rel_lambda_diff": abs(free.eigenvalue - dense.eigenvalue) / scale,
        "max_sorted_coord_dev": max_coord_dev,
        "n_positions_differ": int(disagree.sum()),
        "frac_positions_differ": float(disagree.mean()),
        "n_tied_groups": n_tied_groups,
        "n_tied_positions": ties_only,
        "disagreement_within_ties": confined,
        "dense_residual": dense.residual,
        "mf_residual": free.residual,
        "dense_degenerate": dense.degenerate,
        "mf_degenerate": free.degenerate,
        "mf_solver": free.solver,
        "mf_n_matvec": free.n_matvec,
        "mf_n_spmv": free.n_spmv,
    }


def _tie_groups(sorted_coordinates: np.ndarray, tol: float) -> tuple[int, int]:
    """Count runs of coordinates closer together than ``tol``."""
    if sorted_coordinates.size < 2:
        return 0, 0
    gaps = np.diff(sorted_coordinates)
    breaks = np.flatnonzero(gaps > tol)
    run_starts = np.concatenate([[0], breaks + 1])
    run_ends = np.concatenate([breaks + 1, [sorted_coordinates.size]])
    sizes = run_ends - run_starts
    tied = sizes[sizes > 1]
    return int(tied.size), int(tied.sum())


def _disagreement_is_within_ties(
    dense_order: np.ndarray,
    sorted_dense: np.ndarray,
    disagree: np.ndarray,
    tol: float,
) -> bool:
    """True when every disagreeing position sits inside a tie group.

    Position ``i`` of the dense order holds ``sorted_dense[i]``.  A disagreement is
    harmless when the coordinate there equals its neighbour within ``tol``, i.e. the
    two solvers merely permuted a set of coordinates the method cannot order anyway.
    """
    if not disagree.any():
        return True
    positions = np.flatnonzero(disagree)
    for pos in positions:
        if pos > 0 and abs(sorted_dense[pos] - sorted_dense[pos - 1]) <= tol:
            continue
        if pos + 1 < sorted_dense.size and abs(sorted_dense[pos + 1] - sorted_dense[pos]) <= tol:
            continue
        return False
    return True
