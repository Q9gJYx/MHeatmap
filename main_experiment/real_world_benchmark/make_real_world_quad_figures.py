# ruff: noqa: E402, I001

from __future__ import annotations

import sys
from pathlib import Path

sys.modules.setdefault("numexpr", None)

REPO_ROOT = Path(__file__).resolve().parents[2]
BENCHMARK_DIR = Path(__file__).resolve().parent
SYNTHETIC_DIR = BENCHMARK_DIR.parent / "synthetic_benchmark"
for path in (BENCHMARK_DIR, SYNTHETIC_DIR):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap, LogNorm, Normalize

from _benchmark_utils import DEFAULT_WIDTH_GRID, mwb_auc, tw_reorder_for_alpha
from _real_world_data import FIGURES_DIR, prepare_real_world_datasets
from mheatmap import mosaic_heatmap

ALPHAS = (1.0, 2.0, 4.0, 6.0, 8.0, 12.0)
CIP_BLUE = "#2b8cbe"
CIP_MONO_CMAP = ListedColormap([CIP_BLUE])
DISPLAY_THRESHOLDS = {
    "lodes_clean": 2000.0,
}


def _select_tw_for_figure(matrix: np.ndarray):
    best_result = None
    best_score = -np.inf
    best_alpha = -np.inf
    for alpha in ALPHAS:
        result = tw_reorder_for_alpha(matrix, alpha=alpha)
        score = mwb_auc(result.matrix, DEFAULT_WIDTH_GRID)
        if score > best_score + 1e-12 or (
            abs(score - best_score) <= 1e-12 and alpha > best_alpha
        ):
            best_result = result
            best_score = score
            best_alpha = alpha
    if best_result is None:
        raise ValueError("Could not select a TW ordering for plotting.")
    return best_result


def _make_norm(values: np.ndarray, threshold: float = 0.0):
    if threshold > 0:
        values = values[values > threshold]
    vmin = float(values.min())
    vmax = float(values.max())
    if np.isclose(vmin, vmax):
        spread = max(abs(vmin) * 0.1, 1e-3)
        return Normalize(vmin=vmin - spread, vmax=vmax + spread)
    return LogNorm(vmin=vmin, vmax=vmax)


def _style_panel(ax: plt.Axes, title: str) -> None:
    ax.set_title(title, fontsize=10, pad=3)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.45)
        spine.set_color("black")


def _imshow_panel(
    ax: plt.Axes,
    matrix: np.ndarray,
    norm,
    threshold: float = 0.0,
    cmap: str = "YlGnBu",
) -> None:
    masked = np.ma.masked_where(matrix <= threshold, matrix)
    ax.imshow(
        masked,
        aspect="auto",
        interpolation="nearest",
        cmap=cmap,
        norm=norm,
    )


def _draw_support_panel(ax: plt.Axes, matrix: np.ndarray) -> None:
    rows, cols = np.nonzero(matrix > 0)
    if rows.size:
        ax.scatter(
            cols,
            rows,
            s=2.0,
            marker="s",
            c=CIP_BLUE,
            linewidths=0.0,
            rasterized=True,
        )
    ax.set_xlim(-0.5, matrix.shape[1] - 0.5)
    ax.set_ylim(matrix.shape[0] - 0.5, -0.5)
    ax.set_aspect("auto")


def _draw_mosaic(
    ax: plt.Axes,
    matrix: np.ndarray,
    norm,
    cmap: str | ListedColormap,
    threshold: float = 0.0,
) -> None:
    mosaic_heatmap(
        matrix,
        ax=ax,
        cmap=cmap,
        mask=matrix <= threshold,
        norm=norm,
        cbar=False,
        square=False,
        xticklabels=False,
        yticklabels=False,
        rasterized=True,
        linewidths=0.0,
    )


def _plot_dataset_quad(dataset_key: str, matrix: np.ndarray, output_path: Path) -> None:
    tw_matrix = _select_tw_for_figure(matrix).matrix
    positive = np.concatenate([matrix[matrix > 0], tw_matrix[tw_matrix > 0]])
    threshold = DISPLAY_THRESHOLDS.get(dataset_key, 0.0)
    support_style = dataset_key in {"cip_soc", "naics_sic"}

    if support_style:
        norm = Normalize(vmin=0.0, vmax=1.0)
        mosaic_cmap: str | ListedColormap = CIP_MONO_CMAP
    else:
        norm = _make_norm(positive, threshold)
        mosaic_cmap = "Blues" if dataset_key == "openalex" else "YlGnBu"

    fig, axes = plt.subplots(
        1,
        4,
        figsize=(9.8, 2.9),
        gridspec_kw={"wspace": 0.08, "width_ratios": [1.0, 1.0, 1.15, 1.15]},
        constrained_layout=False,
    )

    if support_style:
        _draw_support_panel(axes[0], matrix)
        _draw_support_panel(axes[1], tw_matrix)
    else:
        imshow_cmap = "Blues" if dataset_key == "openalex" else "YlGnBu"
        _imshow_panel(axes[0], matrix, norm, threshold, imshow_cmap)
        _imshow_panel(axes[1], tw_matrix, norm, threshold, imshow_cmap)
    _draw_mosaic(axes[2], matrix, norm, mosaic_cmap, threshold)
    _draw_mosaic(axes[3], tw_matrix, norm, mosaic_cmap, threshold)

    titles = (
        "Original",
        "TW(Two-Walk Laplacian)" if dataset_key == "openalex" else "TW",
        "Mosaic",
        "Mheatmap(TW+Mosaic)",
    )
    for ax, title in zip(axes, titles, strict=True):
        _style_panel(ax, title)

    fig.subplots_adjust(left=0.015, right=0.995, top=0.88, bottom=0.06, wspace=0.06)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    datasets = prepare_real_world_datasets()
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    for dataset in datasets:
        output_path = FIGURES_DIR / f"real_world_quad_{dataset.key}.png"
        _plot_dataset_quad(dataset.key, dataset.matrix.astype(float), output_path)
        print(output_path)


if __name__ == "__main__":
    main()
