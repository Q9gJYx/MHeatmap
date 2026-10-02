# ruff: noqa: E402, I001

"""Run the rebuttal alpha ablation with the paper-spec TW implementation.

The submitted-grid Table 1 is kept separate from this expanded sensitivity
analysis.  Alpha zero is diagnosed spectrally rather than assigned an ordering:
Eq. (5) is block diagonal at zero, so its joint row/column Fiedler coordinate is
not identifiable.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

sys.modules.setdefault("numexpr", None)

REPO_ROOT = Path(__file__).resolve().parents[2]
SYNTHETIC_DIR = REPO_ROOT / "main_experiment" / "synthetic_benchmark"
REAL_DIR = REPO_ROOT / "main_experiment" / "real_world_benchmark"
for path in (SYNTHETIC_DIR, REAL_DIR):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.linalg import eigvalsh
from scipy.sparse import bmat, csr_matrix
from scipy.sparse.csgraph import connected_components

from _benchmark_utils import (
    DEFAULT_WIDTH_GRID,
    mheatmap_two_walk_laplacian,
    mwb_auc,
    normalized_band_mass,
    normalized_two_sum,
    one_walk_reorder,
    tw_reorder_for_alpha,
)
from _real_world_data import prepare_real_world_datasets
from run_synthetic_evaluation import FAMILIES, SIZES, build_case


OUTPUT_DIR = REPO_ROOT / "output" / "main_experiment" / "alpha_sensitivity"
FIGURE_DIR = OUTPUT_DIR / "figures"
NUM_SEEDS = 20
SUBMITTED_ALPHAS = (1.0, 2.0, 4.0, 6.0, 8.0, 12.0)
ALPHAS = (0.01, 0.05, 0.1, 0.25, 0.5, 0.75, *SUBMITTED_ALPHAS, 24.0)
TIE_TOLERANCE = 1e-12


def choose_global_alpha(records: pd.DataFrame) -> float:
    """Choose one alpha by mean synthetic MWB-AUC; ties favor smaller alpha."""
    means = records.groupby("alpha", as_index=False)["tw_mwb_auc"].mean()
    best = float(means["tw_mwb_auc"].max())
    tied = means.loc[
        np.abs(means["tw_mwb_auc"] - best) <= TIE_TOLERANCE,
        "alpha",
    ]
    return float(tied.min())


def select_best_rows(
    records: pd.DataFrame,
    alphas: tuple[float, ...],
) -> pd.DataFrame:
    """Select maximum MWB-AUC per case; ties favor smaller alpha."""
    subset = records.loc[records["alpha"].isin(alphas)].copy()
    subset = subset.sort_values(
        ["case_id", "tw_mwb_auc", "alpha"],
        ascending=[True, False, True],
        kind="stable",
    )
    return subset.drop_duplicates("case_id", keep="first").reset_index(drop=True)


def _metric_values(matrix: np.ndarray) -> dict[str, float]:
    return {
        "two_sum": normalized_two_sum(matrix),
        "band_mass_10": normalized_band_mass(matrix, 0.10),
        "mwb_auc": mwb_auc(matrix, DEFAULT_WIDTH_GRID),
    }


def _evaluate_case(
    matrix: np.ndarray,
    metadata: dict[str, object],
) -> list[dict[str, object]]:
    ow = one_walk_reorder(matrix)
    ow_metrics = _metric_values(ow.matrix)
    rows: list[dict[str, object]] = []
    for alpha in ALPHAS:
        started = time.perf_counter()
        tw = tw_reorder_for_alpha(matrix, alpha=alpha)
        tw_metrics = _metric_values(tw.matrix)
        rows.append(
            {
                **metadata,
                "alpha": alpha,
                "component_count": tw.component_count,
                "runtime_seconds": time.perf_counter() - started,
                **{f"tw_{key}": value for key, value in tw_metrics.items()},
                **{f"ow_{key}": value for key, value in ow_metrics.items()},
            }
        )
    return rows


def run_synthetic() -> pd.DataFrame:
    records: list[dict[str, object]] = []
    for family in FAMILIES:
        for size in SIZES:
            started = time.perf_counter()
            for seed in range(NUM_SEEDS):
                matrix = build_case(family, size, seed)["observed"]
                case_id = f"synthetic/{family.key}/{size.key}/seed-{seed:02d}"
                records.extend(
                    _evaluate_case(
                        matrix,
                        {
                            "suite": "synthetic",
                            "case_id": case_id,
                            "family_key": family.key,
                            "family": family.name,
                            "size_key": size.key,
                            "size": size.name,
                            "dataset": "",
                            "seed": seed,
                            "n_rows": matrix.shape[0],
                            "n_cols": matrix.shape[1],
                        },
                    )
                )
            print(
                f"alpha sweep {family.name} / {size.name}: "
                f"{time.perf_counter() - started:.1f}s",
                flush=True,
            )
    return pd.DataFrame.from_records(records)


def run_real() -> pd.DataFrame:
    records: list[dict[str, object]] = []
    for dataset in prepare_real_world_datasets(rebuild=False):
        started = time.perf_counter()
        matrix = dataset.matrix
        records.extend(
            _evaluate_case(
                matrix,
                {
                    "suite": "real",
                    "case_id": f"real/{dataset.key}",
                    "family_key": "",
                    "family": "",
                    "size_key": "",
                    "size": "",
                    "dataset": dataset.name,
                    "seed": np.nan,
                    "n_rows": matrix.shape[0],
                    "n_cols": matrix.shape[1],
                },
            )
        )
        print(
            f"alpha sweep {dataset.name}: {time.perf_counter() - started:.1f}s",
            flush=True,
        )
    return pd.DataFrame.from_records(records)


def _support_components(matrix: np.ndarray) -> list[tuple[np.ndarray, np.ndarray]]:
    active_rows = np.flatnonzero(matrix.sum(axis=1) > 0)
    active_cols = np.flatnonzero(matrix.sum(axis=0) > 0)
    support = csr_matrix(matrix[np.ix_(active_rows, active_cols)] > 0)
    graph = bmat([[None, support], [support.T, None]], format="csr")
    count, labels = connected_components(graph, directed=False, return_labels=True)
    row_labels = labels[: len(active_rows)]
    col_labels = labels[len(active_rows) :]
    return [
        (
            active_rows[row_labels == component_id],
            active_cols[col_labels == component_id],
        )
        for component_id in range(int(count))
    ]


def alpha_zero_diagnostic(
    matrix: np.ndarray,
    metadata: dict[str, object],
) -> dict[str, object]:
    """Numerically verify the two zero modes per support component at alpha=0."""
    components = _support_components(np.asarray(matrix, dtype=float))
    normalized = matrix / float(matrix.max())
    observed_nullity = 0
    component_nullities: list[int] = []
    for rows, cols in components:
        block = normalized[np.ix_(rows, cols)]
        laplacian = np.asarray(
            mheatmap_two_walk_laplacian(block, alpha=0.0),
            dtype=float,
        )
        eigenvalues = eigvalsh(laplacian, check_finite=False)
        tolerance = (
            100.0
            * np.finfo(float).eps
            * max(1, laplacian.shape[0])
            * max(1.0, float(np.max(np.abs(eigenvalues))))
        )
        nullity = int(np.count_nonzero(np.abs(eigenvalues) <= tolerance))
        component_nullities.append(nullity)
        observed_nullity += nullity
    component_count = len(components)
    return {
        **metadata,
        "alpha": 0.0,
        "component_count": component_count,
        "expected_nullity": 2 * component_count,
        "observed_nullity": observed_nullity,
        "component_nullities": json.dumps(component_nullities),
        "joint_order_identifiable": False,
    }


def run_alpha_zero_diagnostics() -> pd.DataFrame:
    records: list[dict[str, object]] = []
    # One deterministic representative from every synthetic regime is enough to
    # verify the algebraic claim; the claim itself applies to every instance.
    for family in FAMILIES:
        for size in SIZES:
            matrix = build_case(family, size, 0)["observed"]
            records.append(
                alpha_zero_diagnostic(
                    matrix,
                    {
                        "suite": "synthetic",
                        "case_id": f"synthetic/{family.key}/{size.key}/seed-00",
                        "family": family.name,
                        "size": size.name,
                        "dataset": "",
                        "n_rows": matrix.shape[0],
                        "n_cols": matrix.shape[1],
                    },
                )
            )
    for dataset in prepare_real_world_datasets(rebuild=False):
        matrix = dataset.matrix
        records.append(
            alpha_zero_diagnostic(
                matrix,
                {
                    "suite": "real",
                    "case_id": f"real/{dataset.key}",
                    "family": "",
                    "size": "",
                    "dataset": dataset.name,
                    "n_rows": matrix.shape[0],
                    "n_cols": matrix.shape[1],
                },
            )
        )
    return pd.DataFrame.from_records(records)


def summarize_sensitivity(records: pd.DataFrame) -> pd.DataFrame:
    data = records.copy()
    data["tw_wins_two_sum"] = data["tw_two_sum"] < data["ow_two_sum"] - TIE_TOLERANCE
    data["tw_wins_band"] = (
        data["tw_band_mass_10"] > data["ow_band_mass_10"] + TIE_TOLERANCE
    )
    data["tw_wins_mwb"] = data["tw_mwb_auc"] > data["ow_mwb_auc"] + TIE_TOLERANCE
    data["delta_two_sum"] = data["tw_two_sum"] - data["ow_two_sum"]
    data["delta_band_mass_10"] = data["tw_band_mass_10"] - data["ow_band_mass_10"]
    data["delta_mwb_auc"] = data["tw_mwb_auc"] - data["ow_mwb_auc"]

    pieces: list[pd.DataFrame] = []
    for suite, keys in (
        ("synthetic", ["suite", "family", "size", "alpha"]),
        ("real", ["suite", "dataset", "alpha"]),
    ):
        subset = data.loc[data["suite"] == suite]
        summary = subset.groupby(keys, as_index=False).agg(
            n_cases=("case_id", "nunique"),
            tw_two_sum_mean=("tw_two_sum", "mean"),
            tw_two_sum_std=("tw_two_sum", "std"),
            tw_band_mass_10_mean=("tw_band_mass_10", "mean"),
            tw_band_mass_10_std=("tw_band_mass_10", "std"),
            tw_mwb_auc_mean=("tw_mwb_auc", "mean"),
            tw_mwb_auc_std=("tw_mwb_auc", "std"),
            ow_two_sum_mean=("ow_two_sum", "mean"),
            ow_band_mass_10_mean=("ow_band_mass_10", "mean"),
            ow_mwb_auc_mean=("ow_mwb_auc", "mean"),
            delta_two_sum_mean=("delta_two_sum", "mean"),
            delta_band_mass_10_mean=("delta_band_mass_10", "mean"),
            delta_mwb_auc_mean=("delta_mwb_auc", "mean"),
            tw_win_rate_two_sum=("tw_wins_two_sum", "mean"),
            tw_win_rate_band=("tw_wins_band", "mean"),
            tw_win_rate_mwb=("tw_wins_mwb", "mean"),
        )
        pieces.append(summary)
    return pd.concat(pieces, ignore_index=True, sort=False)


def overall_sensitivity(records: pd.DataFrame) -> pd.DataFrame:
    data = records.copy()
    data["tw_wins_two_sum"] = data["tw_two_sum"] < data["ow_two_sum"] - TIE_TOLERANCE
    data["tw_wins_band"] = (
        data["tw_band_mass_10"] > data["ow_band_mass_10"] + TIE_TOLERANCE
    )
    data["tw_wins_mwb"] = data["tw_mwb_auc"] > data["ow_mwb_auc"] + TIE_TOLERANCE
    data["delta_two_sum"] = data["tw_two_sum"] - data["ow_two_sum"]
    data["delta_band_mass_10"] = data["tw_band_mass_10"] - data["ow_band_mass_10"]
    data["delta_mwb_auc"] = data["tw_mwb_auc"] - data["ow_mwb_auc"]
    return data.groupby(["suite", "alpha"], as_index=False).agg(
        n_cases=("case_id", "nunique"),
        tw_two_sum_mean=("tw_two_sum", "mean"),
        tw_two_sum_std=("tw_two_sum", "std"),
        tw_band_mass_10_mean=("tw_band_mass_10", "mean"),
        tw_band_mass_10_std=("tw_band_mass_10", "std"),
        tw_mwb_auc_mean=("tw_mwb_auc", "mean"),
        tw_mwb_auc_std=("tw_mwb_auc", "std"),
        ow_two_sum_mean=("ow_two_sum", "mean"),
        ow_band_mass_10_mean=("ow_band_mass_10", "mean"),
        ow_mwb_auc_mean=("ow_mwb_auc", "mean"),
        delta_two_sum_mean=("delta_two_sum", "mean"),
        delta_band_mass_10_mean=("delta_band_mass_10", "mean"),
        delta_mwb_auc_mean=("delta_mwb_auc", "mean"),
        tw_win_rate_two_sum=("tw_wins_two_sum", "mean"),
        tw_win_rate_band=("tw_wins_band", "mean"),
        tw_win_rate_mwb=("tw_wins_mwb", "mean"),
    )


def best_alpha_outputs(records: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    submitted = select_best_rows(records, SUBMITTED_ALPHAS).assign(
        search_grid="submitted"
    )
    expanded = select_best_rows(records, ALPHAS).assign(search_grid="expanded")
    selected = pd.concat([submitted, expanded], ignore_index=True)

    summaries: list[dict[str, object]] = []
    for (grid, suite), subset in selected.groupby(["search_grid", "suite"]):
        if suite == "synthetic":
            groups = subset.groupby(["family", "size"], sort=False)
        else:
            groups = subset.groupby(["dataset"], sort=False)
        for group_key, group in groups:
            label = group_key if isinstance(group_key, str) else " / ".join(group_key)
            counts = group["alpha"].value_counts().sort_index()
            summaries.append(
                {
                    "search_grid": grid,
                    "suite": suite,
                    "group": label,
                    "n_cases": len(group),
                    "mean_best_alpha": float(group["alpha"].mean()),
                    "median_best_alpha": float(group["alpha"].median()),
                    "mode_best_alpha": float(counts.index[np.argmax(counts.to_numpy())]),
                    "best_alpha_counts": json.dumps(
                        {f"{alpha:g}": int(count) for alpha, count in counts.items()},
                        sort_keys=True,
                    ),
                }
            )
    return selected, pd.DataFrame.from_records(summaries)


def method_comparison(records: pd.DataFrame, global_alpha: float) -> pd.DataFrame:
    selections = {
        "One-Walk": records.drop_duplicates("case_id", keep="first"),
        "TW fixed alpha=1": records.loc[records["alpha"] == 1.0],
        f"TW global alpha={global_alpha:g}": records.loc[
            records["alpha"] == global_alpha
        ],
        "TW adaptive submitted grid": select_best_rows(records, SUBMITTED_ALPHAS),
        "TW adaptive expanded grid": select_best_rows(records, ALPHAS),
    }
    rows: list[dict[str, object]] = []
    for method, data in selections.items():
        for suite, subset in data.groupby("suite"):
            if method == "One-Walk":
                metric_prefix = "ow"
                selected_alpha_mean = np.nan
            else:
                metric_prefix = "tw"
                selected_alpha_mean = float(subset["alpha"].mean())
            rows.append(
                {
                    "suite": suite,
                    "method": method,
                    "n_cases": subset["case_id"].nunique(),
                    "selected_alpha_mean": selected_alpha_mean,
                    "two_sum_mean": float(subset[f"{metric_prefix}_two_sum"].mean()),
                    "band_mass_10_mean": float(
                        subset[f"{metric_prefix}_band_mass_10"].mean()
                    ),
                    "mwb_auc_mean": float(subset[f"{metric_prefix}_mwb_auc"].mean()),
                }
            )
    return pd.DataFrame.from_records(rows)


def save_figure(overall: pd.DataFrame, global_alpha: float) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    specs = (
        ("tw_two_sum_mean", "ow_two_sum_mean", "Mean normalized 2-SUM"),
        ("tw_band_mass_10_mean", "ow_band_mass_10_mean", "Mean Band@10%"),
        ("tw_mwb_auc_mean", "ow_mwb_auc_mean", "Mean MWB-AUC"),
    )
    suites = (("synthetic", "Synthetic (300 instances)"), ("real", "Real (7 datasets)"))
    fig, axes = plt.subplots(2, 3, figsize=(12.0, 6.2), constrained_layout=True)
    for row, (suite, suite_label) in enumerate(suites):
        data = overall.loc[overall["suite"] == suite].sort_values("alpha")
        for column_index, (tw_column, ow_column, label) in enumerate(specs):
            axis = axes[row, column_index]
            axis.plot(
                data["alpha"],
                data[tw_column],
                color="tab:blue",
                marker="o",
                label="Two-Walk",
            )
            axis.axhline(
                float(data[ow_column].iloc[0]),
                color="black",
                linestyle="--",
                linewidth=1.2,
                label="One-Walk",
            )
            axis.axvline(
                1.0,
                color="0.5",
                linestyle=":",
                linewidth=1.0,
                label=r"fixed $\alpha=1$" if column_index == 0 else None,
            )
            axis.axvline(
                global_alpha,
                color="tab:red",
                linestyle=":",
                linewidth=1.2,
                label=(
                    rf"global $\alpha={global_alpha:g}$"
                    if column_index == 0
                    else None
                ),
            )
            axis.set_xscale("log")
            axis.set_xlabel(r"$\alpha$")
            axis.set_ylabel(label)
            axis.grid(alpha=0.25)
            if column_index == 0:
                axis.set_title(suite_label, loc="left", fontweight="bold")
                axis.legend(frameon=False, fontsize=8)
    fig.savefig(FIGURE_DIR / "alpha_sensitivity.png", dpi=200)
    fig.savefig(FIGURE_DIR / "alpha_sensitivity.pdf")
    plt.close(fig)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    synthetic = run_synthetic()
    real = run_real()
    records = pd.concat([synthetic, real], ignore_index=True)
    records.to_csv(
        OUTPUT_DIR / "alpha_sensitivity_per_instance.csv",
        index=False,
        lineterminator="\n",
    )

    summary = summarize_sensitivity(records)
    summary.to_csv(
        OUTPUT_DIR / "alpha_sensitivity_by_group.csv",
        index=False,
        lineterminator="\n",
    )
    overall = overall_sensitivity(records)
    overall.to_csv(
        OUTPUT_DIR / "alpha_sensitivity_overall.csv",
        index=False,
        lineterminator="\n",
    )

    global_alpha = choose_global_alpha(synthetic)
    selected, selected_summary = best_alpha_outputs(records)
    selected.to_csv(
        OUTPUT_DIR / "best_alpha_per_instance.csv",
        index=False,
        lineterminator="\n",
    )
    selected_summary.to_csv(
        OUTPUT_DIR / "best_alpha_summary.csv",
        index=False,
        lineterminator="\n",
    )
    comparison = method_comparison(records, global_alpha)
    comparison.to_csv(
        OUTPUT_DIR / "fixed_and_adaptive_comparison.csv",
        index=False,
        lineterminator="\n",
    )

    zero = run_alpha_zero_diagnostics()
    zero.to_csv(
        OUTPUT_DIR / "alpha_zero_diagnostics.csv",
        index=False,
        lineterminator="\n",
    )
    save_figure(overall, global_alpha)

    metadata = {
        "positive_alpha_grid": list(ALPHAS),
        "submitted_alpha_grid": list(SUBMITTED_ALPHAS),
        "global_alpha_selection_suite": "synthetic only (300 instances)",
        "global_alpha_selection_objective": "mean retained MWB-AUC",
        "global_alpha": global_alpha,
        "tie_break": "smaller alpha within 1e-12",
        "alpha_zero_policy": (
            "spectral degeneracy diagnostic only; no ordering metric is reported"
        ),
        "elapsed_seconds": time.perf_counter() - started,
    }
    (OUTPUT_DIR / "experiment_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    print(f"global alpha (synthetic-only): {global_alpha:g}", flush=True)
    print(f"outputs: {OUTPUT_DIR}", flush=True)


if __name__ == "__main__":
    main()
