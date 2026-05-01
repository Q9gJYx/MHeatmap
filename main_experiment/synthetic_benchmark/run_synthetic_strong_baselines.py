# ruff: noqa: E402, I001

from __future__ import annotations

import sys
from pathlib import Path

sys.modules.setdefault("numexpr", None)

REPO_ROOT = Path(__file__).resolve().parents[2]
BENCHMARK_DIR = Path(__file__).resolve().parent
for path in (BENCHMARK_DIR,):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import numpy as np
import pandas as pd

from _benchmark_utils import (
    DEFAULT_WIDTH_GRID,
    ReorderResult,
    ca_svd_reorder,
    diagonal_band_mass,
    hierarchical_olo_reorder,
    marginal_sort_reorder,
    median_reorder,
    mwb_auc,
    normalized_two_sum,
    one_walk_reorder,
    tw_auto_reorder,
)
from run_synthetic_evaluation import FAMILIES, SIZES, build_case

OUTPUT_DIR = REPO_ROOT / "output" / "main_experiment" / "synthetic_benchmark"
PROCESSED_DIR = OUTPUT_DIR

NUM_SEEDS = 20
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
PAPER_METHOD_ORDER = METHOD_ORDER
PAPER_SIZE_ORDER = ("Large", "Medium", "Small")
PAPER_FAMILY_ORDER = (
    "Family A: Clean one-to-one",
    "Family B: Paired subgroup overlap",
    "Family C: Shared super-prototype",
    "Family D: Shared prototype with noise",
    "Family E: Cross-block leakage",
)


def ensure_directories() -> None:
    for directory in (OUTPUT_DIR,):
        directory.mkdir(parents=True, exist_ok=True)


def save_paper_table(pivot: pd.DataFrame) -> None:
    methods = PAPER_METHOD_ORDER
    metric_specs = (
        ("two_sum_mean", "$10^2 \\times$ 2-SUM $\\downarrow$", False),
        ("band_mass_10_mean", "Band@10\\% $\\uparrow$", True),
        ("mwb_auc_mean", "MWB-AUC $\\uparrow$", True),
    )

    family_rank = {family: rank for rank, family in enumerate(PAPER_FAMILY_ORDER)}
    size_rank = {size: rank for rank, size in enumerate(PAPER_SIZE_ORDER)}
    paper = pivot.loc[:, ["family_key", "family", "size_key", "size"]].copy()
    paper["family_rank"] = paper["family"].map(family_rank)
    paper["size_rank"] = paper["size"].map(size_rank)
    for metric, _, _ in metric_specs:
        for method in methods:
            values = pivot[f"{metric}|{method}"].astype(float)
            if metric == "two_sum_mean":
                values = values * 100.0
            paper[f"{metric}|{method}"] = values
    paper = paper.sort_values(["family_rank", "size_rank"]).reset_index(drop=True)

    csv_columns = ["family", "size"]
    for metric, _, _ in metric_specs:
        for method in methods:
            csv_columns.append(f"{metric}|{method}")
            paper[f"{metric}|{method}"] = paper[f"{metric}|{method}"].round(2)
    paper.loc[:, csv_columns].to_csv(
        PROCESSED_DIR / "synthetic_strong_baselines_paper.csv",
        index=False,
    )

    md_lines = [
        "# Synthetic Rectangular Benchmark (Paper Table)",
        "",
        "| Family | Size | "
        + " | ".join(
            f"{title} {METHOD_SHORT[method]}"
            for _, title, _ in metric_specs
            for method in methods
        )
        + " |",
        "|---|---|" + "---:|" * (len(metric_specs) * len(methods)),
    ]
    for _, row in paper.iterrows():
        cells = [row["family"], row["size"]]
        for metric, _, _ in metric_specs:
            for method in methods:
                cells.append(f"{float(row[f'{metric}|{method}']):.2f}")
        md_lines.append("| " + " | ".join(cells) + " |")
    (PROCESSED_DIR / "synthetic_strong_baselines_paper.md").write_text(
        "\n".join(md_lines) + "\n",
        encoding="utf-8",
    )


def save_outputs(records: list[dict[str, object]]) -> None:
    df = pd.DataFrame.from_records(records)

    summary = (
        df.groupby(
            ["family", "family_key", "size", "size_key", "method"],
            as_index=False,
        )
        .agg(
            two_sum_mean=("two_sum", "mean"),
            band_mass_10_mean=("band_mass_10", "mean"),
            mwb_auc_mean=("mwb_auc", "mean"),
        )
    )

    pivot = summary.pivot_table(
        index=["family_key", "family", "size_key", "size"],
        columns="method",
        values=["two_sum_mean", "band_mass_10_mean", "mwb_auc_mean"],
    )
    pivot.columns = [f"{metric}|{method}" for metric, method in pivot.columns]
    pivot = pivot.reset_index()
    pivot = pivot.sort_values(["family_key", "size_key"]).reset_index(drop=True)

    save_paper_table(pivot)


def main() -> None:
    ensure_directories()

    reorderers = {
        "Marginal": marginal_sort_reorder,
        "HC+OLO": hierarchical_olo_reorder,
        "One-walk": one_walk_reorder,
        "CA-SVD": lambda matrix: ca_svd_reorder(matrix, widths=DEFAULT_WIDTH_GRID),
        "Median": lambda matrix: median_reorder(matrix, widths=DEFAULT_WIDTH_GRID),
    }

    records: list[dict[str, object]] = []

    for family in FAMILIES:
        for size in SIZES:
            for seed in range(NUM_SEEDS):
                payload = build_case(family, size, seed)
                observed = payload["observed"]

                results: dict[str, ReorderResult] = {
                    "Original": ReorderResult(
                        observed.copy(),
                        np.arange(observed.shape[0], dtype=int),
                        np.arange(observed.shape[1], dtype=int),
                    )
                }

                for method, reorder in reorderers.items():
                    results[method] = reorder(observed)

                tw_auto = tw_auto_reorder(observed, alphas=ALPHAS, widths=DEFAULT_WIDTH_GRID)
                results["TW"] = tw_auto.reorder

                for method in METHOD_ORDER:
                    matrix = results[method].matrix
                    records.append(
                        {
                            "family": family.name,
                            "family_key": family.key,
                            "size": size.name,
                            "size_key": size.key,
                            "seed": seed,
                            "method": method,
                            "two_sum": normalized_two_sum(matrix),
                            "band_mass_10": diagonal_band_mass(matrix, 0.10),
                            "mwb_auc": mwb_auc(matrix, DEFAULT_WIDTH_GRID),
                        }
                    )

    save_outputs(records)


if __name__ == "__main__":
    main()
