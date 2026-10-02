# ruff: noqa: E402, I001

from __future__ import annotations

import sys
import time
from pathlib import Path

sys.modules.setdefault("numexpr", None)

REPO_ROOT = Path(__file__).resolve().parents[2]
BENCHMARK_DIR = Path(__file__).resolve().parent
SYNTHETIC_DIR = BENCHMARK_DIR.parent / "synthetic_benchmark"
for path in (BENCHMARK_DIR, SYNTHETIC_DIR):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import numpy as np
import pandas as pd

from _benchmark_utils import (
    DEFAULT_WIDTH_GRID,
    ReorderResult,
    ca_svd_reorder,
    hierarchical_olo_reorder,
    marginal_sort_reorder,
    median_reorder,
    mwb_auc,
    normalized_band_mass,
    normalized_two_sum,
    one_walk_reorder,
    tw_auto_reorder,
)
from _real_world_data import OUTPUT_DIR, prepare_real_world_datasets

ALPHAS = (1.0, 2.0, 4.0, 6.0, 8.0, 12.0)
METHOD_ORDER = (
    "Original",
    "Marginal",
    "HC+OLO",
    "One-walk",
    "CA-SVD",
    "Median",
    "TW",
)
METHOD_SHORT = {
    "Original": "O",
    "Marginal": "M",
    "HC+OLO": "HC",
    "One-walk": "OW",
    "CA-SVD": "CA",
    "Median": "MD",
    "TW": "TW",
}


def save_paper_table(records: list[dict[str, object]]) -> None:
    df = pd.DataFrame.from_records(records)
    selected = df.loc[
        df["method"] == "TW",
        ["dataset_index", "dataset", "shape", "selected_alpha"],
    ].copy()
    selected.to_csv(
        OUTPUT_DIR / "real_world_selected_alpha.csv",
        index=False,
        lineterminator="\n",
    )
    pivot = df.pivot_table(
        index=["dataset_index", "dataset", "shape"],
        columns="method",
        values=["two_sum", "band_mass_10", "mwb_auc"],
    )
    pivot.columns = [f"{metric}|{method}" for metric, method in pivot.columns]
    pivot = pivot.reset_index().sort_values("dataset_index").reset_index(drop=True)

    metric_specs = (
        ("two_sum", "$10^2 \\times$ 2-SUM $\\downarrow$"),
        ("band_mass_10", "Band@10\\% $\\uparrow$"),
        ("mwb_auc", "MWB-AUC $\\uparrow$"),
    )

    paper = pivot.loc[:, ["dataset", "shape"]].copy()
    alpha_by_dataset = selected.set_index("dataset")["selected_alpha"]
    paper["best_alpha"] = paper["dataset"].map(alpha_by_dataset).astype(float)
    for metric, _ in metric_specs:
        for method in METHOD_ORDER:
            values = pivot[f"{metric}|{method}"].astype(float)
            if metric == "two_sum":
                values = values * 100.0
            paper[f"{metric}|{method}"] = values.round(2)

    csv_columns = ["dataset", "shape", "best_alpha"]
    for metric, _ in metric_specs:
        for method in METHOD_ORDER:
            csv_columns.append(f"{metric}|{method}")
    paper.loc[:, csv_columns].to_csv(
        OUTPUT_DIR / "real_world_benchmark_paper.csv",
        index=False,
        lineterminator="\n",
    )

    md_lines = [
        "# Real-World Rectangular Benchmark (Paper Table)",
        "",
        (
            "Best $\\alpha$ maximizes MWB-AUC on the submitted grid; exact "
            "ties use the smaller value."
        ),
        "",
        "| Dataset | Shape | Best $\\alpha$ | "
        + " | ".join(
            f"{title} {METHOD_SHORT[method]}"
            for metric, title in metric_specs
            for method in METHOD_ORDER
        )
        + " |",
        "|---|---|---:|" + "---:|" * (len(metric_specs) * len(METHOD_ORDER)),
    ]
    for _, row in paper.iterrows():
        cells = [row["dataset"], row["shape"], f"{float(row['best_alpha']):g}"]
        for metric, _ in metric_specs:
            for method in METHOD_ORDER:
                cells.append(f"{float(row[f'{metric}|{method}']):.2f}")
        md_lines.append("| " + " | ".join(cells) + " |")

    (OUTPUT_DIR / "real_world_benchmark_paper.md").write_text(
        "\n".join(md_lines) + "\n",
        encoding="utf-8",
        newline="\n",
    )


def main() -> None:
    datasets = prepare_real_world_datasets()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    reorderers = {
        "Marginal": marginal_sort_reorder,
        "HC+OLO": hierarchical_olo_reorder,
        "One-walk": one_walk_reorder,
        "CA-SVD": lambda matrix: ca_svd_reorder(matrix, widths=DEFAULT_WIDTH_GRID),
        "Median": lambda matrix: median_reorder(matrix, widths=DEFAULT_WIDTH_GRID),
    }

    records: list[dict[str, object]] = []
    for dataset_index, dataset in enumerate(datasets):
        dataset_started = time.perf_counter()
        original_matrix = dataset.matrix
        results: dict[str, ReorderResult] = {
            "Original": ReorderResult(
                original_matrix.copy(),
                np.arange(original_matrix.shape[0], dtype=int),
                np.arange(original_matrix.shape[1], dtype=int),
            )
        }

        for method, reorder in reorderers.items():
            results[method] = reorder(original_matrix)

        tw_auto = tw_auto_reorder(original_matrix, alphas=ALPHAS, widths=DEFAULT_WIDTH_GRID)
        results["TW"] = tw_auto.reorder

        for method in METHOD_ORDER:
            matrix = results[method].matrix
            records.append(
                {
                    "dataset_index": dataset_index,
                    "dataset": dataset.name,
                    "method": method,
                    "shape": f"{matrix.shape[0]}x{matrix.shape[1]}",
                    "two_sum": normalized_two_sum(matrix),
                    "band_mass_10": normalized_band_mass(matrix, 0.10),
                    "mwb_auc": mwb_auc(matrix, DEFAULT_WIDTH_GRID),
                    "selected_alpha": tw_auto.alpha if method == "TW" else np.nan,
                }
            )
        print(
            f"{dataset.name}: alpha={tw_auto.alpha:g}, "
            f"components={tw_auto.reorder.component_count}, "
            f"{time.perf_counter() - dataset_started:.1f}s",
            flush=True,
        )

    save_paper_table(records)


if __name__ == "__main__":
    main()
