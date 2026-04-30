from __future__ import annotations

import csv
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.cluster.hierarchy import leaves_list, linkage, optimal_leaf_ordering
from scipy.linalg import eigh
from scipy.sparse.linalg import svds
from scipy.spatial.distance import pdist

sys.modules.setdefault("numexpr", None)

from mheatmap.graph import (
    copermute_from_bipermute,
    spectral_permute,
    two_walk_laplacian,
)

DEFAULT_WIDTH_GRID = np.linspace(0.02, 0.50, 25)
EPS = 1e-12


def diagonal_band_mass(matrix: np.ndarray, frac: float) -> float:
    total = float(matrix.sum())
    if total <= 0:
        return 0.0

    n_rows, n_cols = matrix.shape
    band = max(1, int(frac * min(n_rows, n_cols)))
    score = 0.0

    for row_index in range(n_rows):
        center = row_index * n_cols / n_rows
        start = max(0, int(np.floor(center - band)))
        stop = min(n_cols, int(np.ceil(center + band + 1)))
        score += float(matrix[row_index, start:stop].sum())

    return score / total


@dataclass(frozen=True)
class ReorderResult:
    matrix: np.ndarray
    row_order: np.ndarray
    col_order: np.ndarray


@dataclass(frozen=True)
class AutoTWResult:
    reorder: ReorderResult
    alpha: float
    selection_score: float


def load_matrix_csv(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with path.open(encoding="utf-8") as handle:
        rows = list(csv.reader(handle))
    col_labels = np.array(rows[0][1:])
    row_labels = np.array([row[0] for row in rows[1:]])
    matrix = np.array([[float(value) for value in row[1:]] for row in rows[1:]])
    return matrix, row_labels, col_labels


def apply_orders(
    matrix: np.ndarray,
    row_order: np.ndarray,
    col_order: np.ndarray,
) -> np.ndarray:
    return matrix[row_order][:, col_order]


def recover_column_order(
    original_matrix: np.ndarray,
    row_order: np.ndarray,
    reordered_matrix: np.ndarray,
) -> np.ndarray:
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
            raise ValueError("Could not recover reordered column order from matrix.")
        recovered.append(matches.pop(0))
    return np.array(recovered, dtype=int)


def normalized_two_sum(matrix: np.ndarray) -> float:
    total = float(np.sum(matrix))
    if total <= 0:
        return 0.0

    n_rows, n_cols = matrix.shape
    row_positions = np.linspace(0.0, 1.0, n_rows)[:, None]
    if n_cols == 1:
        col_positions = np.zeros((1, 1))
    else:
        col_positions = np.linspace(0.0, 1.0, n_cols)[None, :]
    squared_distance = (row_positions - col_positions) ** 2
    return float(np.sum(matrix * squared_distance) / total)


def mwb_auc(matrix: np.ndarray, widths: np.ndarray) -> float:
    scores = np.array([diagonal_band_mass(matrix, width) for width in widths])
    return float(np.trapezoid(scores, widths) / (widths[-1] - widths[0]))


def orient_orders_for_diagonal(
    matrix: np.ndarray,
    row_order: np.ndarray,
    col_order: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    candidates = [
        (row_order, col_order),
        (row_order[::-1], col_order),
        (row_order, col_order[::-1]),
        (row_order[::-1], col_order[::-1]),
    ]

    best_score = None
    best_pair = candidates[0]
    for candidate_rows, candidate_cols in candidates:
        candidate_matrix = apply_orders(matrix, candidate_rows, candidate_cols)
        score = normalized_two_sum(candidate_matrix)
        if best_score is None or score < best_score:
            best_score = score
            best_pair = (candidate_rows, candidate_cols)
    return best_pair


def orient_orders_for_mwb_auc(
    matrix: np.ndarray,
    row_order: np.ndarray,
    col_order: np.ndarray,
    widths: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    widths = DEFAULT_WIDTH_GRID if widths is None else widths
    candidates = [
        (row_order, col_order),
        (row_order[::-1], col_order[::-1]),
    ]

    best_score = None
    best_two_sum = None
    best_pair = candidates[0]
    for candidate_rows, candidate_cols in candidates:
        candidate_matrix = apply_orders(matrix, candidate_rows, candidate_cols)
        score = mwb_auc(candidate_matrix, widths)
        two_sum = normalized_two_sum(candidate_matrix)
        if (
            best_score is None
            or score > best_score + 1e-12
            or (abs(score - best_score) <= 1e-12 and two_sum < best_two_sum - 1e-12)
        ):
            best_score = score
            best_two_sum = two_sum
            best_pair = (candidate_rows, candidate_cols)
    return best_pair


def nonzero_axis_indices(
    matrix: np.ndarray,
    axis: int,
) -> tuple[np.ndarray, np.ndarray]:
    sums = np.sum(matrix, axis=axis)
    nonzero = np.where(sums > 0)[0]
    zero = np.where(sums <= 0)[0]
    return nonzero, zero


def _append_zero_axes(
    row_nonzero: np.ndarray,
    row_local: np.ndarray,
    row_zero: np.ndarray,
    col_nonzero: np.ndarray,
    col_local: np.ndarray,
    col_zero: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    row_order = np.concatenate([row_nonzero[row_local], row_zero])
    col_order = np.concatenate([col_nonzero[col_local], col_zero])
    return row_order, col_order


def _top_svd_components(
    matrix: np.ndarray,
    max_axes: int = 2,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    min_dim = min(matrix.shape)
    if min_dim <= 1:
        raise ValueError("SVD requires both matrix dimensions to exceed one.")

    if min_dim <= 64:
        u, s, vt = np.linalg.svd(matrix, full_matrices=False)
        k = min(max_axes, u.shape[1], vt.shape[0])
        return u[:, :k], s[:k], vt[:k, :]

    k = min(max_axes, min_dim - 1)
    if k <= 0:
        u, s, vt = np.linalg.svd(matrix, full_matrices=False)
        return u[:, :1], s[:1], vt[:1, :]

    try:
        u, s, vt = svds(matrix, k=k, which="LM")
        order = np.argsort(s)[::-1]
        return u[:, order], s[order], vt[order, :]
    except Exception:
        u, s, vt = np.linalg.svd(matrix, full_matrices=False)
        k = min(max_axes, u.shape[1], vt.shape[0])
        return u[:, :k], s[:k], vt[:k, :]


def _best_svd_axis_reorder(
    matrix: np.ndarray,
    scores_left: np.ndarray,
    scores_right: np.ndarray,
    row_nonzero: np.ndarray,
    row_zero: np.ndarray,
    col_nonzero: np.ndarray,
    col_zero: np.ndarray,
    widths: np.ndarray | None = None,
) -> ReorderResult:
    widths = DEFAULT_WIDTH_GRID if widths is None else widths

    best_result: ReorderResult | None = None
    best_score: float | None = None
    best_two_sum: float | None = None

    n_axes = min(scores_left.shape[1], scores_right.shape[0])
    for axis in range(n_axes):
        row_local = np.argsort(scores_left[:, axis], kind="mergesort")
        col_local = np.argsort(scores_right[axis, :], kind="mergesort")
        row_order, col_order = _append_zero_axes(
            row_nonzero,
            row_local,
            row_zero,
            col_nonzero,
            col_local,
            col_zero,
        )
        row_order, col_order = orient_orders_for_mwb_auc(
            matrix,
            row_order,
            col_order,
            widths=widths,
        )
        candidate_matrix = apply_orders(matrix, row_order, col_order)
        candidate_score = mwb_auc(candidate_matrix, widths)
        candidate_two_sum = normalized_two_sum(candidate_matrix)
        if (
            best_result is None
            or candidate_score > best_score + 1e-12
            or (
                abs(candidate_score - best_score) <= 1e-12
                and candidate_two_sum < best_two_sum - 1e-12
            )
        ):
            best_result = ReorderResult(candidate_matrix, row_order, col_order)
            best_score = candidate_score
            best_two_sum = candidate_two_sum

    if best_result is None:
        raise ValueError("Could not derive an SVD-based ordering.")
    return best_result


def _center_of_mass(values: np.ndarray, axis_length: int) -> np.ndarray:
    positions = np.arange(axis_length, dtype=float)
    totals = values.sum(axis=1)
    weighted = values @ positions
    return np.divide(
        weighted,
        totals,
        out=np.full(values.shape[0], axis_length, dtype=float),
        where=totals > 0,
    )


def marginal_sort_reorder(matrix: np.ndarray) -> ReorderResult:
    row_sums = np.sum(matrix, axis=1)
    col_sums = np.sum(matrix, axis=0)
    row_centers = _center_of_mass(matrix, matrix.shape[1])
    col_centers = _center_of_mass(matrix.T, matrix.shape[0])

    row_order = np.lexsort((row_centers, -row_sums))
    col_order = np.lexsort((col_centers, -col_sums))
    row_order, col_order = orient_orders_for_diagonal(matrix, row_order, col_order)
    return ReorderResult(
        apply_orders(matrix, row_order, col_order),
        row_order,
        col_order,
    )


def ca_svd_reorder(
    matrix: np.ndarray,
    widths: np.ndarray | None = None,
    eps: float = EPS,
) -> ReorderResult:
    row_nonzero, row_zero = nonzero_axis_indices(matrix, axis=1)
    col_nonzero, col_zero = nonzero_axis_indices(matrix, axis=0)

    if len(row_nonzero) == 0 or len(col_nonzero) == 0:
        row_order = np.arange(matrix.shape[0], dtype=int)
        col_order = np.arange(matrix.shape[1], dtype=int)
        return ReorderResult(matrix.copy(), row_order, col_order)

    B_sub = matrix[np.ix_(row_nonzero, col_nonzero)].astype(float, copy=False)
    total = float(np.sum(B_sub))
    if total <= 0:
        row_order = np.arange(matrix.shape[0], dtype=int)
        col_order = np.arange(matrix.shape[1], dtype=int)
        return ReorderResult(matrix.copy(), row_order, col_order)

    P = B_sub / total
    r = P.sum(axis=1)
    c = P.sum(axis=0)
    expected = np.outer(r, c)
    residual = (P - expected) / np.sqrt(expected + eps)
    u, _, vt = _top_svd_components(residual, max_axes=2)
    return _best_svd_axis_reorder(
        matrix,
        u,
        vt,
        row_nonzero,
        row_zero,
        col_nonzero,
        col_zero,
        widths=widths,
    )


def normsvd_reorder(
    matrix: np.ndarray,
    widths: np.ndarray | None = None,
    eps: float = EPS,
) -> ReorderResult:
    row_nonzero, row_zero = nonzero_axis_indices(matrix, axis=1)
    col_nonzero, col_zero = nonzero_axis_indices(matrix, axis=0)

    if len(row_nonzero) == 0 or len(col_nonzero) == 0:
        row_order = np.arange(matrix.shape[0], dtype=int)
        col_order = np.arange(matrix.shape[1], dtype=int)
        return ReorderResult(matrix.copy(), row_order, col_order)

    B_sub = matrix[np.ix_(row_nonzero, col_nonzero)].astype(float, copy=False)
    row_sums = np.sum(B_sub, axis=1)
    col_sums = np.sum(B_sub, axis=0)
    normalized = (
        (1.0 / np.sqrt(row_sums + eps))[:, None]
        * B_sub
        * (1.0 / np.sqrt(col_sums + eps))[None, :]
    )
    u, _, vt = _top_svd_components(normalized, max_axes=2)
    return _best_svd_axis_reorder(
        matrix,
        u,
        vt,
        row_nonzero,
        row_zero,
        col_nonzero,
        col_zero,
        widths=widths,
    )


def _olo_order(vectors: np.ndarray) -> np.ndarray:
    count = vectors.shape[0]
    if count <= 1:
        return np.arange(count, dtype=int)
    if count == 2:
        totals = np.sum(vectors, axis=1)
        return np.argsort(-totals, kind="mergesort")

    distances = pdist(vectors, metric="cosine")
    if not np.isfinite(distances).all() or np.allclose(distances, 0):
        totals = np.sum(vectors, axis=1)
        centers = _center_of_mass(vectors, vectors.shape[1])
        return np.lexsort((centers, -totals))

    tree = linkage(distances, method="average")
    ordered_tree = optimal_leaf_ordering(tree, distances)
    return leaves_list(ordered_tree).astype(int)


def hierarchical_olo_reorder(matrix: np.ndarray) -> ReorderResult:
    row_nonzero, row_zero = nonzero_axis_indices(matrix, axis=1)
    col_nonzero, col_zero = nonzero_axis_indices(matrix, axis=0)

    if len(row_nonzero) == 0 or len(col_nonzero) == 0:
        row_order = np.arange(matrix.shape[0], dtype=int)
        col_order = np.arange(matrix.shape[1], dtype=int)
        return ReorderResult(matrix.copy(), row_order, col_order)

    row_local = _olo_order(matrix[row_nonzero][:, col_nonzero])
    col_local = _olo_order(matrix[row_nonzero][:, col_nonzero].T)
    row_order = np.concatenate([row_nonzero[row_local], row_zero])
    col_order = np.concatenate([col_nonzero[col_local], col_zero])
    row_order, col_order = orient_orders_for_diagonal(matrix, row_order, col_order)
    return ReorderResult(
        apply_orders(matrix, row_order, col_order),
        row_order,
        col_order,
    )


def one_walk_reorder(matrix: np.ndarray) -> ReorderResult:
    rows, cols = matrix.shape
    row_sums = np.sum(matrix, axis=1)
    col_sums = np.sum(matrix, axis=0)
    nonzero_rows = np.where(row_sums > 0)[0]
    nonzero_cols = np.where(col_sums > 0)[0]

    if len(nonzero_rows) == 0 or len(nonzero_cols) == 0:
        row_order = np.arange(rows, dtype=int)
        col_order = np.arange(cols, dtype=int)
        return ReorderResult(matrix.copy(), row_order, col_order)

    B_sub = matrix[np.ix_(nonzero_rows, nonzero_cols)].astype(float)
    if np.max(B_sub) > 0:
        B_sub = B_sub / np.max(B_sub)

    zeros_rr = np.zeros((B_sub.shape[0], B_sub.shape[0]), dtype=float)
    zeros_cc = np.zeros((B_sub.shape[1], B_sub.shape[1]), dtype=float)
    adjacency = np.block([[zeros_rr, B_sub], [B_sub.T, zeros_cc]])
    degree = np.diag(np.sum(adjacency, axis=1))
    laplacian = degree - adjacency
    eigenvalues, eigenvectors = eigh(laplacian)
    fiedler_index = np.where(np.abs(eigenvalues) > 1e-10)[0][0]
    bipartite_order = np.argsort(eigenvectors[:, fiedler_index])
    row_order, col_order = copermute_from_bipermute(
        [rows, cols],
        nonzero_rows,
        nonzero_cols,
        bipartite_order,
    )
    row_order, col_order = orient_orders_for_diagonal(matrix, row_order, col_order)
    return ReorderResult(
        apply_orders(matrix, row_order, col_order),
        row_order,
        col_order,
    )


def _weighted_median(positions: np.ndarray, weights: np.ndarray) -> float:
    if len(positions) == 0:
        return 0.0
    order = np.argsort(positions, kind="mergesort")
    sorted_positions = positions[order]
    sorted_weights = weights[order]
    cumulative = np.cumsum(sorted_weights)
    threshold = 0.5 * float(np.sum(sorted_weights))
    index = int(np.searchsorted(cumulative, threshold, side="left"))
    index = min(index, len(sorted_positions) - 1)
    return float(sorted_positions[index])


def _alternating_layer_reorder(
    matrix: np.ndarray,
    use_median: bool,
    max_iter: int = 20,
    init: str = "marginal",
    widths: np.ndarray | None = None,
    eps: float = EPS,
) -> ReorderResult:
    row_nonzero, row_zero = nonzero_axis_indices(matrix, axis=1)
    col_nonzero, col_zero = nonzero_axis_indices(matrix, axis=0)

    if len(row_nonzero) == 0 or len(col_nonzero) == 0:
        row_order = np.arange(matrix.shape[0], dtype=int)
        col_order = np.arange(matrix.shape[1], dtype=int)
        return ReorderResult(matrix.copy(), row_order, col_order)

    B_sub = matrix[np.ix_(row_nonzero, col_nonzero)].astype(float, copy=False)
    n_rows, n_cols = B_sub.shape

    if init == "marginal":
        init_result = marginal_sort_reorder(B_sub)
        row_order = init_result.row_order.copy()
        col_order = init_result.col_order.copy()
    else:
        row_order = np.arange(n_rows, dtype=int)
        col_order = np.arange(n_cols, dtype=int)

    row_neighbors = [np.flatnonzero(B_sub[i] > 0) for i in range(n_rows)]
    col_neighbors = [np.flatnonzero(B_sub[:, j] > 0) for j in range(n_cols)]
    row_weights = [B_sub[i, nbrs] for i, nbrs in enumerate(row_neighbors)]
    col_weights = [B_sub[nbrs, j] for j, nbrs in enumerate(col_neighbors)]

    for _ in range(max_iter):
        old_row_order = row_order.copy()
        old_col_order = col_order.copy()

        col_pos = np.empty(n_cols, dtype=float)
        col_pos[col_order] = np.arange(n_cols, dtype=float)
        old_row_pos = np.empty(n_rows, dtype=float)
        old_row_pos[row_order] = np.arange(n_rows, dtype=float)
        row_score = np.empty(n_rows, dtype=float)
        row_mass = np.sum(B_sub, axis=1)
        for row_idx in range(n_rows):
            if row_mass[row_idx] <= eps or len(row_neighbors[row_idx]) == 0:
                row_score[row_idx] = old_row_pos[row_idx]
            elif use_median:
                row_score[row_idx] = _weighted_median(
                    col_pos[row_neighbors[row_idx]],
                    row_weights[row_idx],
                )
            else:
                row_score[row_idx] = float(
                    np.dot(row_weights[row_idx], col_pos[row_neighbors[row_idx]])
                    / (row_mass[row_idx] + eps)
                )
        row_order = np.lexsort((old_row_pos, row_score))

        row_pos = np.empty(n_rows, dtype=float)
        row_pos[row_order] = np.arange(n_rows, dtype=float)
        old_col_pos = np.empty(n_cols, dtype=float)
        old_col_pos[old_col_order] = np.arange(n_cols, dtype=float)
        col_score = np.empty(n_cols, dtype=float)
        col_mass = np.sum(B_sub, axis=0)
        for col_idx in range(n_cols):
            if col_mass[col_idx] <= eps or len(col_neighbors[col_idx]) == 0:
                col_score[col_idx] = old_col_pos[col_idx]
            elif use_median:
                col_score[col_idx] = _weighted_median(
                    row_pos[col_neighbors[col_idx]],
                    col_weights[col_idx],
                )
            else:
                col_score[col_idx] = float(
                    np.dot(col_weights[col_idx], row_pos[col_neighbors[col_idx]])
                    / (col_mass[col_idx] + eps)
                )
        col_order = np.lexsort((old_col_pos, col_score))

        if np.array_equal(row_order, old_row_order) and np.array_equal(
            col_order,
            old_col_order,
        ):
            break

    row_order, col_order = _append_zero_axes(
        row_nonzero,
        row_order,
        row_zero,
        col_nonzero,
        col_order,
        col_zero,
    )
    row_order, col_order = orient_orders_for_mwb_auc(
        matrix,
        row_order,
        col_order,
        widths=widths,
    )
    return ReorderResult(
        apply_orders(matrix, row_order, col_order),
        row_order,
        col_order,
    )


def barycenter_reorder(
    matrix: np.ndarray,
    max_iter: int = 20,
    init: str = "marginal",
    widths: np.ndarray | None = None,
    eps: float = EPS,
) -> ReorderResult:
    return _alternating_layer_reorder(
        matrix,
        use_median=False,
        max_iter=max_iter,
        init=init,
        widths=widths,
        eps=eps,
    )


def median_reorder(
    matrix: np.ndarray,
    max_iter: int = 20,
    init: str = "marginal",
    widths: np.ndarray | None = None,
    eps: float = EPS,
) -> ReorderResult:
    return _alternating_layer_reorder(
        matrix,
        use_median=True,
        max_iter=max_iter,
        init=init,
        widths=widths,
        eps=eps,
    )


def tw_reorder(matrix: np.ndarray) -> ReorderResult:
    row_labels = np.arange(matrix.shape[0], dtype=int)
    tw_matrix, tw_row_labels = spectral_permute(matrix, row_labels, mode="tw")
    row_order = tw_row_labels.astype(int)
    col_order = recover_column_order(matrix, row_order, tw_matrix)
    row_order, col_order = orient_orders_for_diagonal(matrix, row_order, col_order)
    return ReorderResult(
        apply_orders(matrix, row_order, col_order),
        row_order,
        col_order,
    )


def tw_reorder_for_alpha(matrix: np.ndarray, alpha: float) -> ReorderResult:
    if abs(alpha - 1.0) <= 1e-12:
        return tw_reorder(matrix)
    return tw_alpha_reorder(matrix, alpha=alpha)


def tw_auto_reorder(
    matrix: np.ndarray,
    alphas: tuple[float, ...],
    widths: np.ndarray,
) -> AutoTWResult:
    best_result: ReorderResult | None = None
    best_alpha: float | None = None
    best_score: float | None = None

    for alpha in alphas:
        result = tw_reorder_for_alpha(matrix, alpha=alpha)
        score = mwb_auc(result.matrix, widths)
        if (
            best_result is None
            or score > best_score + 1e-12
            or (abs(score - best_score) <= 1e-12 and alpha < best_alpha)
        ):
            best_result = result
            best_alpha = alpha
            best_score = score

    if best_result is None or best_alpha is None or best_score is None:
        raise ValueError("TW-Auto could not select an alpha.")

    return AutoTWResult(
        reorder=best_result,
        alpha=float(best_alpha),
        selection_score=float(best_score),
    )


def tw_alpha_reorder(matrix: np.ndarray, alpha: float) -> ReorderResult:
    rows, cols = matrix.shape
    row_sums = np.sum(matrix, axis=1)
    col_sums = np.sum(matrix, axis=0)
    nonzero_rows = np.where(row_sums > 0)[0]
    nonzero_cols = np.where(col_sums > 0)[0]

    if len(nonzero_rows) == 0 or len(nonzero_cols) == 0:
        row_order = np.arange(rows, dtype=int)
        col_order = np.arange(cols, dtype=int)
        return ReorderResult(matrix.copy(), row_order, col_order)

    B_sub = matrix[np.ix_(nonzero_rows, nonzero_cols)].astype(float)
    if np.max(B_sub) > 0:
        B_sub = B_sub / np.max(B_sub)

    laplacian = two_walk_laplacian(B_sub, alpha=alpha)
    eigenvalues, eigenvectors = eigh(laplacian)
    fiedler_index = np.where(np.abs(eigenvalues) > 1e-10)[0][0]
    bipartite_order = np.argsort(eigenvectors[:, fiedler_index])
    row_order, col_order = copermute_from_bipermute(
        [rows, cols],
        nonzero_rows,
        nonzero_cols,
        bipartite_order,
    )
    row_order, col_order = orient_orders_for_diagonal(matrix, row_order, col_order)
    return ReorderResult(
        apply_orders(matrix, row_order, col_order),
        row_order,
        col_order,
    )

