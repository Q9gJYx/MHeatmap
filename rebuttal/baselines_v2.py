"""Extra baselines for the MHeatMap rebuttal.

``gram_one_walk_order`` -- reviewer XFNF W5 (workstream A3 / D10)
----------------------------------------------------------------

    "An important baseline appears to be missing. Specifically, the paper does not
     evaluate the approach of applying one-walk Laplacian reordering separately to
     B_bar @ B_bar.T and B_bar.T @ B_bar. Including this baseline would help clarify
     the benefit provided specifically by the proposed two-walk formulation."

This is the natural ablation of the paper's central claim. Two-Walk builds
``A_alpha = M^2 + alpha * M`` on the bipartite lift ``M = [[0, B], [B.T, 0]]``, whose
``M^2`` block is exactly ``blockdiag(B B^T, B^T B)``. Seriating the two Gram blocks
*separately* therefore keeps the shared-neighbour information while dropping both the
joint solve and the ``alpha * M`` term that ties the two axes together. Whatever the
two formulations differ by is what this baseline isolates.

The procedure, per connected component of the active support (Section 3.3), is:

* the row order is the Fiedler order of ``L = D - B B^T``;
* the column order is the Fiedler order of ``L = D - B.T B``;
* components are concatenated by decreasing mass and zero-marginal rows and columns
  are appended in their original order, exactly as the paper-spec pipeline does.

Self-loops need no special handling: in a combinatorial Laplacian the diagonal of
``D`` absorbs ``W_ii`` and then cancels against it, so the diagonal of the Gram --
which is not a similarity between distinct rows -- drops out on its own.

Why the orientation is resolved by ``orient_orders_for_diagonal``
-----------------------------------------------------------------

Solving the two Grams independently gives two eigenproblems with **independent** sign
ambiguity, so relative to the joint case there is a genuinely free choice: flipping
the row order, or the column order, or both. A joint eigenvector ``[v_rows; v_cols]``
cannot express "rows reversed, columns not"; two separate ones can, and nothing in
the mathematics prefers any of the four combinations. Getting it wrong is fatal --
the mass lands on the anti-diagonal -- so it has to be resolved by a stated rule.

This module uses the released rule, ``orient_orders_for_diagonal``: try all four
reversals and keep the lowest 2-SUM. That is deliberately *more* generous than the
treatment One-Walk and Two-Walk get in the corrected tables, which is precisely the
point: commit ``abea5bc`` removed the 2-SUM orientation search for OW/TW because
their joint eigenvector makes it unnecessary, and because it fits a reported metric.
Gram-OW needs it, so it keeps it. If the baseline still loses while being allowed to
orient itself against the metric, the loss is not an artifact of how it was set up.
State that asymmetry wherever the column is reported.

Callers
-------

``make_table1_aligned.baseline_orderings`` uses this as the sixth non-spectral
baseline, so it appears in the Table 1 real-data output. ``make_table1_synthetic``
inherits it through the same function. Access the benchmark module as ``tc.bm``
rather than re-implementing any of its helpers.
"""

from __future__ import annotations

import numpy as np
from scipy.linalg import eigh

import tw_components as tc

bm = tc.bm


def _gram_fiedler_order(gram: np.ndarray) -> np.ndarray:
    """Fiedler order of the combinatorial Laplacian of one Gram block.

    ``gram`` is ``B B^T`` (rows) or ``B.T B`` (columns) for a single connected
    component. The orientation and sort conventions are the benchmark's own:
    ``bm._orient_fiedler`` then a stable argsort, as in ``_block_fiedler``.
    """
    size = gram.shape[0]
    if size < 2:
        # A one-row or one-column component has no second eigenvector to sort by.
        return np.arange(size)

    laplacian = np.diag(gram.sum(axis=1)) - gram
    eigenvalues, eigenvectors = eigh(laplacian, check_finite=False)
    positive = np.where(np.abs(eigenvalues) > tc.EIGEN_TOL)[0]
    if positive.size == 0:
        # Every eigenvalue is ~0: the Gram has no off-diagonal structure at all, so
        # there is no graph to seriate. Leave the component in its original order.
        return np.arange(size)
    return np.argsort(bm._orient_fiedler(eigenvectors[:, positive[0]]), kind="stable")


def gram_one_walk_order(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """One-walk Laplacian reordering applied separately to ``B B^T`` and ``B.T B``.

    Returns ``(row_order, col_order)``, both validated permutations of the input axes.
    """
    n_rows, n_cols = matrix.shape
    active_rows, active_cols, sub = tc.active_submatrix(matrix)
    zero_rows = np.setdiff1d(np.arange(n_rows), active_rows)
    zero_cols = np.setdiff1d(np.arange(n_cols), active_cols)

    if sub.size == 0:
        return (
            np.arange(n_rows, dtype=int),
            np.arange(n_cols, dtype=int),
        )

    ordered_rows: list[np.ndarray] = []
    ordered_cols: list[np.ndarray] = []
    for component in tc.support_components(matrix):
        # Component indices are global; the Gram is built on the compact block.
        local_rows = np.searchsorted(active_rows, component.row_indices)
        local_cols = np.searchsorted(active_cols, component.col_indices)
        block = sub[np.ix_(local_rows, local_cols)]
        ordered_rows.append(component.row_indices[_gram_fiedler_order(block @ block.T)])
        ordered_cols.append(component.col_indices[_gram_fiedler_order(block.T @ block)])

    row_order = np.concatenate([*ordered_rows, zero_rows]).astype(int, copy=False)
    col_order = np.concatenate([*ordered_cols, zero_cols]).astype(int, copy=False)
    bm._validate_permutation(row_order, n_rows, "row")
    bm._validate_permutation(col_order, n_cols, "column")

    # The free relative flip between the two independent solves. See the module
    # docstring: this is intentionally the generous, released convention.
    return bm.orient_orders_for_diagonal(matrix, row_order, col_order)
