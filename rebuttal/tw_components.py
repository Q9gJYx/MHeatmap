"""Single-vector baseline and component census for MHeatMap (paper Section 3.3).

The benchmark originally ordered a *single* Fiedler vector of the whole
active-support Laplacian. Section 3.3 of the paper specifies something different:

    "If the active support has several connected components, we order each
     component independently, arrange components by decreasing total mass, and
     append zero-marginal rows and columns in their original order."

`_benchmark_utils.py` in the repository now **implements Section 3.3 natively**
(`_componentwise_spectral_reorder`, used by `one_walk_reorder` and
`tw_alpha_reorder`), so that is the canonical procedure and
`order_componentwise` below simply delegates to it. This module is the audit
layer around it, and keeps two things the benchmark does not provide:

``order_single``
    The pre-fix behaviour -- one Fiedler vector over the whole active support,
    with the released `orient_orders_for_diagonal` step. This is the "as
    submitted" baseline of record, and is the only remaining implementation of
    it, so Table A depends on it.

``support_components`` / ``component_table``
    The per-component census (count, shape, mass) used by the audit table.

Modes
-----
``"tw"``  Two-Walk: ``A_alpha = M^2 + alpha * M``, ``L = D - A``.
``"ow"``  One-Walk:  ``A = M`` (the bipartite lift).
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

REPO_ROOT = Path(__file__).resolve().parents[1]
SYNTHETIC_DIR = REPO_ROOT / "main_experiment" / "synthetic_benchmark"

# pandas imports numexpr lazily; the benchmark scripts disable it for determinism.
sys.modules.setdefault("numexpr", None)


def _load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


# Reuse the released baselines/metrics rather than reimplementing them.
bm = _load_module("_benchmark_utils", SYNTHETIC_DIR / "_benchmark_utils.py")

from mheatmap.graph import (  # noqa: E402
    copermute_from_bipermute,
    spectral_permute,
    two_walk_laplacian,
)

EIGEN_TOL = 1e-10


@dataclass(frozen=True)
class Component:
    """One connected component of the active bipartite support graph."""

    row_indices: np.ndarray  # global row indices
    col_indices: np.ndarray  # global column indices
    mass: float

    @property
    def n_rows(self) -> int:
        return int(self.row_indices.size)

    @property
    def n_cols(self) -> int:
        return int(self.col_indices.size)


def active_submatrix(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Definition 1: drop rows/columns with zero marginal, normalise by max entry."""
    row_sums = np.sum(matrix, axis=1)
    col_sums = np.sum(matrix, axis=0)
    active_rows = np.where(row_sums > 0)[0]
    active_cols = np.where(col_sums > 0)[0]
    if active_rows.size == 0 or active_cols.size == 0:
        return active_rows, active_cols, np.zeros((0, 0), dtype=float)
    sub = matrix[np.ix_(active_rows, active_cols)].astype(float)
    peak = sub.max()
    if peak > 0:
        sub = sub / peak
    return active_rows, active_cols, sub


def support_components(matrix: np.ndarray) -> list[Component]:
    """Connected components of the active bipartite support graph, by decreasing mass.

    Components are computed on the *active* support (zero-marginal rows/columns
    removed). Counting components before that removal inflates the count by one
    per empty row/column, since each is an isolated vertex.
    """
    active_rows, active_cols, sub = active_submatrix(matrix)
    if sub.size == 0:
        return []

    support = sp.csr_matrix((sub > 0).astype(np.int8))
    graph = sp.bmat([[None, support], [support.T, None]], format="csr")
    _, labels = connected_components(graph, directed=False)
    row_labels = labels[: active_rows.size]
    col_labels = labels[active_rows.size :]

    components: list[Component] = []
    for label in range(int(labels.max()) + 1):
        rows = np.where(row_labels == label)[0]
        cols = np.where(col_labels == label)[0]
        if rows.size == 0 or cols.size == 0:
            continue
        mass = float(sub[np.ix_(rows, cols)].sum())
        components.append(
            Component(active_rows[rows], active_cols[cols], mass)
        )
    components.sort(key=lambda c: -c.mass)
    return components


def _fiedler_bipermutation(sub: np.ndarray, alpha: float, mode: str) -> np.ndarray:
    """Fiedler-vector ordering of a single connected block, in released conventions."""
    if mode == "ow":
        n_rows, n_cols = sub.shape
        adjacency = np.block(
            [
                [np.zeros((n_rows, n_rows)), sub],
                [sub.T, np.zeros((n_cols, n_cols))],
            ]
        )
        laplacian = np.diag(adjacency.sum(axis=1)) - adjacency
    elif mode == "tw":
        laplacian = two_walk_laplacian(sub, alpha=alpha)
    else:
        raise ValueError(f"unknown mode {mode!r}; expected 'tw' or 'ow'")

    eigenvalues, eigenvectors = eigh(laplacian)
    positive = np.where(np.abs(eigenvalues) > EIGEN_TOL)[0]
    if positive.size == 0:
        # Degenerate: every eigenvalue is ~0 (e.g. alpha = 0 on a block-diagonal
        # adjacency). The coordinated order is undefined; fall back to identity.
        return np.arange(laplacian.shape[0])
    return np.argsort(eigenvectors[:, positive[0]])


def _recover_column_order(
    original_matrix: np.ndarray,
    row_order: np.ndarray,
    reordered_matrix: np.ndarray,
) -> np.ndarray:
    """Submitted helper: match permuted columns back to their originals by value."""
    row_reordered_original = original_matrix[row_order, :]

    column_lookup: dict[bytes, list[int]] = {}
    for col_index in range(row_reordered_original.shape[1]):
        key = np.ascontiguousarray(row_reordered_original[:, col_index]).tobytes()
        column_lookup.setdefault(key, []).append(col_index)

    recovered = []
    for col_index in range(reordered_matrix.shape[1]):
        key = np.ascontiguousarray(reordered_matrix[:, col_index]).tobytes()
        matches = column_lookup.get(key)
        if not matches:
            raise ValueError("could not recover reordered column order from matrix")
        recovered.append(matches.pop(0))
    return np.array(recovered, dtype=int)


def _submitted_tw_alpha1(matrix: np.ndarray):
    """The submitted Two-Walk path at alpha = 1.

    The submitted `tw_reorder_for_alpha` short-circuited alpha == 1 to the external
    `mheatmap.graph.spectral_permute`, whose preprocessing prunes rows below
    ``1e-3`` of the total, and then recovered the column order by matching permuted
    column vectors. Every alpha > 1 took a different, local code path. A baseline
    claiming to be "as submitted" has to reproduce that split rather than
    approximate it, otherwise its alpha selection is not the submitted one.
    """
    row_labels = np.arange(matrix.shape[0], dtype=int)
    tw_matrix, tw_row_labels = spectral_permute(matrix, row_labels, mode="tw")
    row_order = tw_row_labels.astype(int)
    col_order = _recover_column_order(matrix, row_order, tw_matrix)
    return bm.orient_orders_for_diagonal(matrix, row_order, col_order)


def order_single(matrix: np.ndarray, alpha: float = 1.0, mode: str = "tw"):
    """Released behaviour: one Fiedler vector over the whole active support.

    Kept as the baseline of record for Table 1; see module docstring.
    """
    if mode == "tw" and abs(alpha - 1.0) <= 1e-12:
        return _submitted_tw_alpha1(matrix)

    n_rows, n_cols = matrix.shape
    active_rows, active_cols, sub = active_submatrix(matrix)
    if sub.size == 0:
        return np.arange(n_rows), np.arange(n_cols)

    bipermutation = _fiedler_bipermutation(sub, alpha, mode)
    row_order, col_order = copermute_from_bipermute(
        [n_rows, n_cols], active_rows, active_cols, bipermutation
    )
    return bm.orient_orders_for_diagonal(matrix, row_order, col_order)


def order_componentwise(matrix: np.ndarray, alpha: float = 1.0, mode: str = "tw"):
    """Paper Section 3.3, via the benchmark's canonical implementation.

    Delegates to `_benchmark_utils`, which owns the procedure. Two conventions
    differ from the local implementation this module used to carry, and the
    benchmark's are authoritative:

    * components are normalised by the *global* matrix maximum rather than
      renormalised per block (this moves Two-Walk only -- the One-Walk
      combinatorial Laplacian is scale-homogeneous);
    * the assembled order is *not* re-oriented by `orient_orders_for_diagonal`,
      which searched four row/column reversals and kept the lowest 2-SUM. That
      search minimised the metric being reported, so its removal is deliberate.
    """
    if mode == "ow":
        result = bm.one_walk_reorder(matrix)
    elif mode == "tw":
        result = bm.tw_alpha_reorder(matrix, alpha=float(alpha))
    else:
        raise ValueError(f"unknown mode {mode!r}; expected 'tw' or 'ow'")
    return result.row_order, result.col_order


def component_table(matrix: np.ndarray) -> list[dict]:
    """Per-component diagnostics for the audit table."""
    _, _, sub = active_submatrix(matrix)
    total = float(sub.sum()) if sub.size else 0.0
    rows = []
    for rank, component in enumerate(support_components(matrix), start=1):
        rows.append(
            {
                "rank": rank,
                "n_rows": component.n_rows,
                "n_cols": component.n_cols,
                "mass": component.mass,
                "mass_fraction": (component.mass / total) if total > 0 else 0.0,
            }
        )
    return rows
