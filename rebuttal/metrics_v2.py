"""Corrected and extended metrics for the MHeatMap rebuttal.

Two definitions of each diagonal metric are provided, deliberately:

``*_released``
    What the released implementation computes, and therefore what produced
    Table 1. Kept so corrected numbers can be compared against the published
    ones and so reviewer kFbK W1 can be answered with both figures side by side.

``band_mass_eq16`` / ``mwb_auc_block``
    Paper-faithful definitions from Section 4.5 / Eq. (16). These are what the
    manuscript actually describes.

Known divergences (reviewer kFbK W1):

* Band@10% -- Section 4.5 / Eq. (16) is the mass fraction with
  ``|rho_i - kappa_j| <= tau`` on normalised axes. The released code instead uses
  an integer index window ``band = max(1, int(frac * min(m, n)))``, which agrees
  only when the aspect ratio is near 1 (MBTA 8x125 uses +-1 column, Eq. (16)
  implies +-12.4).
* MWB-AUC -- Section 4.5 describes matched contiguous diagonal blocks over a cut
  sweep. The released code integrates diagonal-band mass over a width grid.

``r2s`` matches Eq. (12) and the released implementation; it is not affected.
"""

from __future__ import annotations

import numpy as np

DEFAULT_WIDTH_GRID = np.linspace(0.02, 0.50, 25)
DEFAULT_CUTS = 20


def _positions(n_rows: int, n_cols: int, row_order, col_order):
    """Normalised axis coordinates rho_i and kappa_j (Eq. 15)."""
    row_pos = np.empty(n_rows)
    row_pos[row_order] = np.arange(n_rows) / max(n_rows - 1, 1)
    col_pos = np.empty(n_cols)
    col_pos[col_order] = np.arange(n_cols) / max(n_cols - 1, 1)
    return row_pos, col_pos


def r2s(matrix: np.ndarray, row_order, col_order) -> float:
    """Normalised 2-SUM x100 (Eq. 12). Lower is better."""
    n_rows, n_cols = matrix.shape
    row_pos, col_pos = _positions(n_rows, n_cols, row_order, col_order)
    total = float(matrix.sum())
    if total <= 0:
        return 0.0
    squared = (row_pos[:, None] - col_pos[None, :]) ** 2
    return float((matrix * squared).sum() / total) * 100.0


def band_mass_eq16(matrix: np.ndarray, row_order, col_order, tau: float = 0.10) -> float:
    """Eq. (16) Band@tau. Higher is better."""
    n_rows, n_cols = matrix.shape
    row_pos, col_pos = _positions(n_rows, n_cols, row_order, col_order)
    total = float(matrix.sum())
    if total <= 0:
        return 0.0
    mask = np.abs(row_pos[:, None] - col_pos[None, :]) <= tau
    return float((matrix * mask).sum() / total)


def mwb_auc_block(
    matrix: np.ndarray, row_order, col_order, cuts: int = DEFAULT_CUTS
) -> float:
    """Section 4.5 MWB-AUC: mass retained within matched contiguous diagonal blocks.

    At cut scale ``k`` each axis is split into ``k`` equal bins and a cell counts
    as retained when both coordinates land in the same bin index. The score is the
    mean retained mass fraction over ``k = 2 .. cuts``. Higher is better.
    """
    n_rows, n_cols = matrix.shape
    row_pos, col_pos = _positions(n_rows, n_cols, row_order, col_order)
    total = float(matrix.sum())
    if total <= 0 or cuts < 2:
        return 0.0
    scores = []
    for k in range(2, cuts + 1):
        row_bins = np.minimum((row_pos * k).astype(int), k - 1)
        col_bins = np.minimum((col_pos * k).astype(int), k - 1)
        mask = row_bins[:, None] == col_bins[None, :]
        scores.append(float((matrix * mask).sum() / total))
    return float(np.mean(scores))


def released_band_mass(bm, matrix: np.ndarray, frac: float = 0.10) -> float:
    """The released integer-window Band@10%, via `_benchmark_utils.diagonal_band_mass`."""
    return float(bm.diagonal_band_mass(matrix, frac))


def released_mwb_auc(bm, matrix: np.ndarray, widths=DEFAULT_WIDTH_GRID) -> float:
    """The released MWB-AUC (band-mass integral), via `_benchmark_utils.mwb_auc`."""
    return float(bm.mwb_auc(matrix, np.asarray(widths)))
