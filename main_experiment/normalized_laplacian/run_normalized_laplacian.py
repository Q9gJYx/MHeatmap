"""Paired final-Laplacian normalization ablation for reviewer JfQK Q2.

The Two-Walk adjacency, component handling, metrics and alpha grid are held
fixed. Both controls are recomputed; the existing Table 1 is not rewritten.
"""

from __future__ import annotations

import argparse
import hashlib
from importlib.metadata import version
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

REPO_ROOT = Path(__file__).resolve().parents[2]
SYNTHETIC_DIR = REPO_ROOT / "main_experiment" / "synthetic_benchmark"
REAL_DIR = REPO_ROOT / "main_experiment" / "real_world_benchmark"
for directory in (SYNTHETIC_DIR, REAL_DIR):
    sys.path.insert(0, str(directory))
sys.modules.setdefault("numexpr", None)
os.environ["MPLCONFIGDIR"] = str(REPO_ROOT / ".cache" / "matplotlib")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_info

from _benchmark_utils import (
    DEFAULT_WIDTH_GRID,
    mwb_auc,
    normalized_band_mass,
    normalized_laplacian_reorder,
    normalized_two_sum,
    one_walk_reorder,
    tw_reorder_for_alpha,
)
from _real_world_data import prepare_real_world_datasets
from run_synthetic_evaluation import FAMILIES, SIZES, build_case

OUTPUT_DIR = REPO_ROOT / "output" / "main_experiment" / "normalized_laplacian"
ALPHAS = (1.0, 2.0, 4.0, 6.0, 8.0, 12.0)
NUM_SEEDS = 20
GLOBAL_ALPHA = 12.0
METRICS = ("two_sum", "band_mass_10", "mwb_auc")
TOLERANCE = 1e-12
BOOTSTRAP_SAMPLES = 10000
BOOTSTRAP_SEED = 20261004
VARIANTS = (
    "c_ow", "n_ow", "c_tw_fixed_1", "n_tw_fixed_1",
    "c_tw_fixed_12", "n_tw_fixed_12", "c_tw_adaptive", "n_tw_adaptive",
)
COMPARISONS = (
    ("OW", "c_ow", "n_ow"),
    ("TW fixed alpha=1", "c_tw_fixed_1", "n_tw_fixed_1"),
    ("TW fixed alpha=12", "c_tw_fixed_12", "n_tw_fixed_12"),
    ("TW adaptive", "c_tw_adaptive", "n_tw_adaptive"),
)
META_COLUMNS = (
    "suite", "case_id", "family_key", "family", "size_key", "size",
    "dataset", "seed", "n_rows", "n_cols", "matrix_sha256",
)


def matrix_hash(matrix: np.ndarray) -> str:
    canonical = np.ascontiguousarray(matrix, dtype="<f8")
    digest = hashlib.sha256()
    digest.update(np.asarray(canonical.shape, dtype="<i8").tobytes())
    digest.update(canonical.tobytes())
    return digest.hexdigest()


def save_csv(frame: pd.DataFrame, filename: str) -> None:
    frame.to_csv(OUTPUT_DIR / filename, index=False, lineterminator="\n")


def evaluate_case(matrix: np.ndarray, metadata: dict) -> tuple[list, list]:
    metadata = {**metadata, "matrix_sha256": matrix_hash(matrix)}
    rows, diagnostics = [], []
    for laplacian in ("combinatorial", "normalized"):
        for method, alpha in (("OW", np.nan), *(("TW", a) for a in ALPHAS)):
            started = time.perf_counter()
            if laplacian == "normalized":
                result = normalized_laplacian_reorder(
                    matrix,
                    mode="one_walk" if method == "OW" else "two_walk",
                    alpha=1.0 if method == "OW" else alpha,
                )
            elif method == "OW":
                result = one_walk_reorder(matrix)
            else:
                result = tw_reorder_for_alpha(matrix, alpha=alpha)
            rows.append({
                **metadata, "laplacian": laplacian, "method": method,
                "alpha": alpha, "component_count": result.component_count,
                "runtime_seconds": time.perf_counter() - started,
                "two_sum": normalized_two_sum(result.matrix),
                "band_mass_10": normalized_band_mass(result.matrix, 0.10),
                "mwb_auc": mwb_auc(result.matrix, DEFAULT_WIDTH_GRID),
                "row_order": json.dumps(result.row_order.tolist()),
                "col_order": json.dumps(result.col_order.tolist()),
            })
            for index, diagnostic in enumerate(result.spectral_diagnostics):
                diagnostics.append({
                    "suite": metadata["suite"], "case_id": metadata["case_id"],
                    "method": method, "alpha": alpha, "component_index": index,
                    **diagnostic,
                })
    return rows, diagnostics


def run_cases() -> tuple[pd.DataFrame, pd.DataFrame]:
    rows, diagnostics = [], []
    for family in FAMILIES:
        for size in SIZES:
            started = time.perf_counter()
            for seed in range(NUM_SEEDS):
                matrix = build_case(family, size, seed)["observed"]
                records, audit = evaluate_case(matrix, {
                    "suite": "synthetic",
                    "case_id": f"synthetic/{family.key}/{size.key}/seed-{seed:02d}",
                    "family_key": family.key, "family": family.name,
                    "size_key": size.key, "size": size.name, "dataset": "",
                    "seed": seed, "n_rows": matrix.shape[0], "n_cols": matrix.shape[1],
                })
                rows.extend(records)
                diagnostics.extend(audit)
                if (seed + 1) % 5 == 0:
                    print(f"{family.key}/{size.key}: {seed + 1}/20 seeds, "
                          f"{time.perf_counter() - started:.1f}s", flush=True)
    for dataset in prepare_real_world_datasets(rebuild=False):
        started = time.perf_counter()
        matrix = dataset.matrix
        records, audit = evaluate_case(matrix, {
            "suite": "real", "case_id": f"real/{dataset.key}",
            "family_key": "", "family": "", "size_key": "", "size": "",
            "dataset": dataset.name, "seed": np.nan,
            "n_rows": matrix.shape[0], "n_cols": matrix.shape[1],
        })
        rows.extend(records)
        diagnostics.extend(audit)
        print(f"{dataset.name}: {time.perf_counter() - started:.1f}s", flush=True)
    return pd.DataFrame(rows), pd.DataFrame(diagnostics)


def select_adaptive(group: pd.DataFrame) -> pd.Series:
    """Exact Table 1 selector policy: MWB and a 1e-12 smaller-alpha tie rule."""
    best = None
    for _, row in group.sort_values("alpha").iterrows():
        if best is None or row["mwb_auc"] > best["mwb_auc"] + TOLERANCE:
            best = row
    if best is None:
        raise ValueError("cannot select alpha from an empty group")
    return best


def protocol_rows(records: pd.DataFrame) -> pd.DataFrame:
    selected = []
    for (case_id, laplacian), group in records.groupby(["case_id", "laplacian"], sort=False):
        prefix = "c" if laplacian == "combinatorial" else "n"
        tw = group.loc[group["method"] == "TW"]
        choices = (
            (f"{prefix}_ow", group.loc[group["method"] == "OW"].iloc[0]),
            (f"{prefix}_tw_fixed_1", tw.loc[tw["alpha"] == 1.0].iloc[0]),
            (f"{prefix}_tw_fixed_12", tw.loc[tw["alpha"] == GLOBAL_ALPHA].iloc[0]),
            (f"{prefix}_tw_adaptive", select_adaptive(tw)),
        )
        for variant, row in choices:
            selected.append({
                **{key: row[key] for key in META_COLUMNS},
                "variant": variant, "selected_alpha": row["alpha"],
                **{key: row[key] for key in METRICS},
            })
    return pd.DataFrame(selected)


def method_summary(protocols: pd.DataFrame) -> pd.DataFrame:
    return protocols.groupby(["suite", "variant"], as_index=False).agg(
        n_cases=("case_id", "nunique"),
        selected_alpha_mean=("selected_alpha", "mean"),
        **{f"{m}_{stat}": (m, stat) for m in METRICS for stat in ("mean", "std")},
    )


def alpha_summary(records: pd.DataFrame) -> pd.DataFrame:
    return records.loc[records["method"] == "TW"].groupby(
        ["suite", "alpha", "laplacian"], as_index=False
    ).agg(n_cases=("case_id", "nunique"), **{
        f"{m}_mean": (m, "mean") for m in METRICS
    })


def paired_deltas(protocols: pd.DataFrame) -> pd.DataFrame:
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    results = []
    for label, control, normalized in COMPARISONS:
        left = protocols.loc[protocols["variant"] == control]
        right = protocols.loc[protocols["variant"] == normalized, ["case_id", *METRICS]]
        pairs = left.merge(right, on="case_id", suffixes=("_c", "_n"), validate="one_to_one")
        for suite, suite_pairs in pairs.groupby("suite", sort=False):
            if suite == "synthetic":
                groups = suite_pairs.groupby(["family", "size"], sort=False)
            else:
                groups = suite_pairs.groupby("dataset", sort=False)
            for key, group in groups:
                differences = np.column_stack([
                    group[f"{m}_n"].to_numpy() - group[f"{m}_c"].to_numpy()
                    for m in METRICS
                ])
                if suite == "synthetic":
                    indexes = rng.integers(0, len(group), (BOOTSTRAP_SAMPLES, len(group)))
                    boot = differences[indexes].mean(axis=1)
                    intervals = np.quantile(boot, [0.025, 0.975], axis=0)
                else:
                    intervals = np.full((2, 3), np.nan)
                common = {
                    "suite": suite, "group_scope": "group",
                    "family": key[0] if suite == "synthetic" else "",
                    "size": key[1] if suite == "synthetic" else "",
                    "dataset": key if suite == "real" else "", "comparison": label,
                }
                results.extend(_paired_metric_rows(group, differences, intervals, common))
            differences = np.column_stack([
                suite_pairs[f"{m}_n"].to_numpy() - suite_pairs[f"{m}_c"].to_numpy()
                for m in METRICS
            ])
            if suite == "synthetic":
                # The generators restart default_rng(seed) in every regime.
                # Resample seed blocks jointly to preserve cross-regime
                # dependence, instead of assuming 15 independent streams.
                macro_differences = []
                for metric in METRICS:
                    panel = suite_pairs.assign(
                        difference=suite_pairs[f"{metric}_n"] - suite_pairs[f"{metric}_c"]
                    ).pivot(index="seed", columns=["family", "size"], values="difference")
                    if panel.isna().any().any():
                        raise ValueError("synthetic seed blocks must form a complete regime panel")
                    macro_differences.append(panel.mean(axis=1).to_numpy())
                macro = np.column_stack(macro_differences)
                indexes = rng.integers(0, len(macro), (BOOTSTRAP_SAMPLES, len(macro)))
                intervals = np.quantile(macro[indexes].mean(axis=1), [0.025, 0.975], axis=0)
            else:
                intervals = np.full((2, 3), np.nan)
            results.extend(_paired_metric_rows(suite_pairs, differences, intervals, {
                "suite": suite, "group_scope": "overall", "family": "", "size": "",
                "dataset": "", "comparison": label,
            }))
    return pd.DataFrame(results)


def _paired_metric_rows(group, differences, intervals, common):
    rows = []
    for index, metric in enumerate(METRICS):
        delta = differences[:, index]
        win = delta < -TOLERANCE if metric == "two_sum" else delta > TOLERANCE
        rows.append({
            **common, "metric": metric, "n_cases": len(group),
            "combinatorial_mean": float(group[f"{metric}_c"].mean()),
            "normalized_mean": float(group[f"{metric}_n"].mean()),
            "delta_mean": float(delta.mean()),
            "delta_std": float(delta.std(ddof=1)) if len(delta) > 1 else np.nan,
            "delta_ci_low": intervals[0, index], "delta_ci_high": intervals[1, index],
            "normalized_win_rate": float(win.mean()),
            "tie_rate": float((np.abs(delta) <= TOLERANCE).mean()),
        })
    return rows


def validate_control_replay(records: pd.DataFrame) -> dict:
    prior_path = REPO_ROOT / "output/main_experiment/alpha_sensitivity/alpha_sensitivity_per_instance.csv"
    prior = pd.read_csv(prior_path)
    controls = records.loc[(records["laplacian"] == "combinatorial") & (records["method"] == "TW")]
    joined = controls.merge(prior, on=["case_id", "alpha"], validate="one_to_one", suffixes=("", "_prior"))
    regression_rows = []
    for metric in METRICS:
        delta = joined[metric] - joined[f"tw_{metric}"]
        regression_rows.extend({
            "case_id": row.case_id, "alpha": row.alpha, "metric": metric,
            "absolute_delta": float(abs(value)),
        } for row, value in zip(joined.itertuples(), delta, strict=True) if abs(value) > 1e-10)
    save_csv(pd.DataFrame(regression_rows, columns=["case_id", "alpha", "metric", "absolute_delta"]),
             "control_regression_differences.csv")
    maximum = max(float((joined[m] - joined[f"tw_{m}"]).abs().max()) for m in METRICS)
    return {
        "prior_alpha_csv_sha256": hashlib.sha256(prior_path.read_bytes()).hexdigest(),
        "replayed_case_alpha_count": len(joined),
        "maximum_metric_absolute_delta": maximum,
        "metric_cells_above_1e_minus_10": len(regression_rows),
    }


def save_figure(summary: pd.DataFrame) -> None:
    target = OUTPUT_DIR / "figures"
    target.mkdir(exist_ok=True)
    fig, axes = plt.subplots(2, 3, figsize=(11.5, 6.0), constrained_layout=True)
    for row, suite in enumerate(("synthetic", "real")):
        data = summary.loc[summary["suite"] == suite].set_index("variant")
        for column, metric in enumerate(METRICS):
            ax = axes[row, column]
            multiplier = 100.0 if metric == "two_sum" else 1.0
            x = np.arange(4)
            c_values = [data.loc[pair[1], f"{metric}_mean"] * multiplier for pair in COMPARISONS]
            n_values = [data.loc[pair[2], f"{metric}_mean"] * multiplier for pair in COMPARISONS]
            ax.bar(x - 0.19, c_values, 0.38, color="#355c87", label="Combinatorial")
            ax.bar(x + 0.19, n_values, 0.38, color="#db8a39", label="Normalized")
            ax.set_xticks(x, ["OW", "TW1", "TW12", "TW adaptive"], rotation=12)
            ax.set_ylabel({"two_sum": "100 x R2S (lower is better)",
                           "band_mass_10": "Band@10% (higher is better)",
                           "mwb_auc": "MWB-AUC (higher is better)"}[metric])
            ax.set_title("Synthetic: 300 instances" if suite == "synthetic" else "Real: 7 datasets")
            ax.grid(axis="y", alpha=0.2)
            ax.set_axisbelow(True)
            if row == 0 and column == 0:
                ax.legend(frameon=False, fontsize=8)
    fig.savefig(target / "normalized_laplacian_comparison.png", dpi=200)
    fig.savefig(target / "normalized_laplacian_comparison.pdf")
    plt.close(fig)


def provenance(records: pd.DataFrame, diagnostics: pd.DataFrame, regression: dict, elapsed: float) -> dict:
    import mheatmap.graph._two_walk_laplacian as package_module
    sources = {
        "runner": Path(__file__), "benchmark_helper": SYNTHETIC_DIR / "_benchmark_utils.py",
        "synthetic_generator": SYNTHETIC_DIR / "run_synthetic_evaluation.py",
        "real_loader": REAL_DIR / "_real_world_data.py",
        "mheatmap_laplacian": Path(package_module.__file__),
    }
    return {
        "case_count": int(records["case_id"].nunique()), "record_count": len(records),
        "synthetic_instances": 300, "real_datasets": 7, "alpha_grid": list(ALPHAS),
        "fixed_alphas": [1.0, GLOBAL_ALPHA], "global_alpha": GLOBAL_ALPHA,
        "global_alpha_origin": "previously selected on synthetic data; frozen, not retuned",
        "coordinate": "f = D^(-1/2) u; L f = lambda D f",
        "self_loops": "excluded from normalization degree; d = diag(L)",
        "projection_normalization": "none; original A_alpha is unchanged",
        "adaptive_selector": "MWB-AUC; same Table 1 grid and 1e-12 smaller-alpha tie policy",
        "bootstrap": {"samples": BOOTSTRAP_SAMPLES, "seed": BOOTSTRAP_SEED,
                      "confidence": 0.95,
                      "overall_synthetic": "joint seed-block bootstrap; equal-weight 15 regimes",
                      "preserves_shared_rng_seeds": True,
                      "interpretation": "conditional seed variability; no real-data CIs"},
        "python": platform.python_version(),
        "packages": {name: version(name) for name in ("numpy", "scipy", "pandas", "mheatmap", "matplotlib")},
        "source_sha256": {name: hashlib.sha256(path.read_bytes()).hexdigest() for name, path in sources.items()},
        "git_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, text=True).strip(),
        "git_status": subprocess.check_output(["git", "status", "--short"], cwd=REPO_ROOT, text=True).splitlines(),
        "threadpools": [{key: pool.get(key) for key in ("internal_api", "num_threads", "version")}
                        for pool in threadpool_info()],
        "control_replay": regression,
        "diagnostics": {
            "component_solves": len(diagnostics),
            "max_generalized_relative_residual": float(diagnostics["generalized_relative_residual"].max()),
            "max_null_mode_residual": float(diagnostics["null_mode_residual"].max()),
            "near_repeated_fiedler_components": int((diagnostics["fiedler_gap"] <= 1e-10).sum()),
            "flat_nontrivial_row_components": int(((diagnostics["n_rows"] > 1) &
                (diagnostics["row_relative_spread"] <= 1e-10)).sum()),
            "flat_nontrivial_col_components": int(((diagnostics["n_cols"] > 1) &
                (diagnostics["col_relative_spread"] <= 1e-10)).sum()),
        }, "elapsed_seconds": elapsed,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summarize-only", action="store_true")
    args = parser.parse_args()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    if args.summarize_only:
        records = pd.read_csv(OUTPUT_DIR / "normalized_laplacian_per_instance.csv", keep_default_na=False)
        records["alpha"] = pd.to_numeric(records["alpha"], errors="coerce")
        records["seed"] = pd.to_numeric(records["seed"], errors="coerce")
        diagnostics = pd.read_csv(OUTPUT_DIR / "spectral_diagnostics.csv")
    else:
        records, diagnostics = run_cases()
        save_csv(records, "normalized_laplacian_per_instance.csv")
        save_csv(diagnostics, "spectral_diagnostics.csv")
    assert len(records) == 307 * 14 and records["case_id"].nunique() == 307
    protocols = protocol_rows(records)
    save_csv(protocols, "protocol_comparison.csv")
    summary = method_summary(protocols)
    save_csv(summary, "method_summary.csv")
    save_csv(alpha_summary(records), "alpha_summary.csv")
    save_csv(paired_deltas(protocols), "paired_deltas.csv")
    regression = validate_control_replay(records)
    save_figure(summary)
    if not args.summarize_only:
        metadata = provenance(records, diagnostics, regression, time.perf_counter() - started)
    else:
        metadata = json.loads((OUTPUT_DIR / "experiment_metadata.json").read_text(encoding="utf-8"))
        metadata["bootstrap"]["overall_synthetic"] = (
            "joint seed-block bootstrap; equal-weight 15 regimes"
        )
        metadata["bootstrap"]["preserves_shared_rng_seeds"] = True
        metadata["postprocessing_source_sha256"] = {
            "runner": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        }
        metadata["postprocessing_elapsed_seconds"] = time.perf_counter() - started
        metadata["control_replay"] = regression
    (OUTPUT_DIR / "experiment_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    print(summary.to_string(index=False), flush=True)
    print("Prior control replay:", regression, flush=True)
    print(f"Outputs: {OUTPUT_DIR}", flush=True)


if __name__ == "__main__":
    main()
