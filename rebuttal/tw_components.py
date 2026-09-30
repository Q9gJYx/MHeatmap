"""Component-wise spectral ordering for MHeatMap (paper Section 3.3).

The released benchmark (`main_experiment/synthetic_benchmark/_benchmark_utils.py`)
orders a *single* Fiedler vector of the whole active-support Laplacian. Section 3.3
of the paper specifies something different:

    "If the active support has several connected components, we order each
     component independently, arrange components by decreasing total mass, and
     append zero-marginal rows and columns in their original order."

This module implements that procedure. It reuses the released conventions
(`two_walk_laplacian`, `copermute_from_bipermute`, `orient_orders_for_diagonal`)
so that when the active support has exactly one connected component the output is
identical to the released code -- see `run_components_audit.py` for the self-check.

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

from mheatmap.graph import copermute_from_bipermute, two_walk_laplacian  # noqa: E402

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


def order_single(matrix: np.ndarray, alpha: float = 1.0, mode: str = "tw"):
    """Released behaviour: one Fiedler vector over the whole active support.

    Kept as the baseline of record for Table 1; see module docstring.
    """
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
    """Paper Section 3.3: order each component, arrange by decreasing mass."""
    n_rows, n_cols = matrix.shape
    active_rows, active_cols, sub = active_submatrix(matrix)
    if sub.size == 0:
        return np.arange(n_rows), np.arange(n_cols)

    ordered_rows: list[np.ndarray] = []
    ordered_cols: list[np.ndarray] = []
    for component in support_components(matrix):
        local = np.ix_(component.row_indices, component.col_indices)
        block = matrix[local].astype(float)
        peak = block.max()
        if peak > 0:
            block = block / peak
        bipermutation = _fiedler_bipermutation(block, alpha, mode)
        local_rows, local_cols = copermute_from_bipermute(
            [component.n_rows, component.n_cols],
            np.arange(component.n_rows),
            np.arange(component.n_cols),
            bipermutation,
        )
        ordered_rows.append(component.row_indices[local_rows])
        ordered_cols.append(component.col_indices[local_cols])

    zero_rows = np.setdiff1d(np.arange(n_rows), active_rows, assume_unique=True)
    zero_cols = np.setdiff1d(np.arange(n_cols), active_cols, assume_unique=True)
    row_order = np.concatenate(ordered_rows + [zero_rows])
    col_order = np.concatenate(ordered_cols + [zero_cols])

    _assert_permutation(row_order, n_rows)
    _assert_permutation(col_order, n_cols)
    return bm.orient_orders_for_diagonal(matrix, row_order, col_order)


def _assert_permutation(order: np.ndarray, size: int) -> None:
    """Guards the active-index -> global-index mapping: a silent mismatch here
    yields a malformed 'order' that still looks plausible in aggregate metrics."""
    if order.size != size or np.unique(order).size != size:
        raise AssertionError(
            f"order is not a permutation of range({size}): size={order.size}, "
            f"unique={np.unique(order).size}"
        )


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
