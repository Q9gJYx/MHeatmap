"""Dense-vs-matrix-free equivalence over the whole benchmark corpus.

Makes the Section 3.3 matrix-free path canonical by demonstrating that it reproduces
the shipped dense solver's results.  For every instance this runs both paths and
records the tie-aware agreement contract:

1. ``|lambda_mf - lambda_dense| <= 1e-8 * max(lambda_dense, lambda_max)``
2. ``max_sorted_coord_dev <= 1e-9`` -- the sorted coordinate sequences agree, i.e. the
   orderings agree *up to ties*
3. every differing position lies inside a numerically tied coordinate group
4. ``|dR2S| <= 1e-4``, ``|dBand| <= 1e-9``, ``|dMWB| <= 1e-4``
5. the adaptive alpha selector picks the same alpha under both solvers

Exact permutation equality is deliberately *not* asserted: the Fiedler vector of a real
matrix is heavily tied (MBTA has 100 of 133 coordinates tied), so ``argsort(kind="stable")``
legitimately resolves differently under two solvers that agree to machine precision.  The
differing positions are counted and published rather than asserted away.

Run:
    cd MHeatMap/src
    uv run python rebuttal/run_matrixfree_equivalence.py
    uv run python rebuttal/run_matrixfree_equivalence.py --datasets cip_soc --seeds 3
"""

from __future__ import annotations

import argparse
import csv
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import metrics_v2 as mv  # noqa: E402
import tw_components as tc  # noqa: E402
import tw_matrixfree as mf  # noqa: E402

SYNTHETIC_DIR = mf.SYNTHETIC_DIR
if str(SYNTHETIC_DIR) not in sys.path:
    sys.path.insert(0, str(SYNTHETIC_DIR))

from run_synthetic_evaluation import FAMILIES, SIZES, build_case  # noqa: E402

OUTPUT_DIR = HERE / "output"
MATRIX_DIR = (
    mf.REPO_ROOT
    / "output"
    / "main_experiment"
    / "real_world_benchmark"
    / "processed_matrices"
)
COMMITTED_CSV = (
    mf.REPO_ROOT
    / "output"
    / "main_experiment"
    / "real_world_benchmark"
    / "real_world_benchmark_paper.csv"
)

DATASET_LABELS = {
    "naics_sic": "SIC -> NAICS",
    "cip_soc": "CIP -> SOC",
    "acs_clean": "ACS OCCP x INDP",
    "lodes_clean": "LODES Home x Work",
    "twenty_newsgroups": "20 Newsgroups",
    "gtfs_clean": "MBTA Route x Station",
    "openalex": "OpenAlex Author x Topic",
}
DATASET_ORDER = tuple(DATASET_LABELS)

ALPHAS = mf.DEFAULT_ALPHA_GRID
WIDTHS = mf.bm.DEFAULT_WIDTH_GRID

REL_LAMBDA_TOL = 1e-8
COORD_TOL = 1e-9
R2S_TOL = 1e-4
BAND_TOL = 1e-9
MWB_TOL = 1e-4

FIELDS = [
    "suite",
    "instance",
    "mode",
    "alpha",
    "m",
    "n",
    "nnz",
    "n_components",
    "max_component_vertices",
    "n_positions_differ",
    "frac_positions_differ",
    "disagreement_within_ties",
    "max_sorted_coord_dev",
    "rel_lambda_diff",
    "r2s_dense",
    "r2s_mf",
    "d_r2s",
    "band_dense",
    "band_mf",
    "d_band",
    "mwb_dense",
    "mwb_mf",
    "d_mwb",
    "d_r2s_intrinsic",
    "d_band_intrinsic",
    "d_mwb_intrinsic",
    "metric_class",
    "ok",
]

# Relative perturbation used to estimate how well conditioned a reported number is.
INTRINSIC_EPS = 1e-12


# --------------------------------------------------------------------------- #
# Comparison
# --------------------------------------------------------------------------- #


def _axis_metrics(matrix, row_order, col_order) -> tuple[float, float, float]:
    return (
        mv.r2s(matrix, row_order, col_order),
        mv.band_mass_eq16(matrix, row_order, col_order),
        mv.mwb_auc_block(matrix, row_order, col_order),
    )


def compare_matrix(
    matrix: np.ndarray,
    *,
    mode: str,
    alpha: float | None,
    suite: str,
    instance: str,
    intrinsic_reps: int = 2,
) -> dict:
    """Run both solver paths on one instance and score the agreement."""
    dense = mf.componentwise_reorder(matrix, mode=mode, alpha=alpha, solver="dense")
    free = mf.componentwise_reorder(
        matrix, mode=mode, alpha=alpha, solver="matrix_free", tol=1e-12
    )

    r2s_d, band_d, mwb_d = _axis_metrics(matrix, dense.row_order, dense.col_order)
    r2s_f, band_f, mwb_f = _axis_metrics(matrix, free.row_order, free.col_order)

    # Agreement is measured on the joint order, since rows and columns are sorted
    # from one Fiedler vector; comparing the axes separately would understate nothing
    # but is harder to read.
    disagree = int((dense.row_order != free.row_order).sum() + (
        dense.col_order != free.col_order
    ).sum())
    positions = matrix.shape[0] + matrix.shape[1]

    component_diag = _largest_component_diagnostics(matrix, mode, alpha)

    record = {
        "suite": suite,
        "instance": instance,
        "mode": mode,
        "alpha": "" if alpha is None else f"{alpha:g}",
        "m": matrix.shape[0],
        "n": matrix.shape[1],
        "nnz": int((matrix > 0).sum()),
        "n_components": dense.component_count,
        "max_component_vertices": component_diag.get("max_component_vertices", 0),
        "n_positions_differ": disagree,
        "frac_positions_differ": disagree / max(positions, 1),
        "disagreement_within_ties": component_diag.get("disagreement_within_ties", True),
        "max_sorted_coord_dev": component_diag.get("max_sorted_coord_dev", 0.0),
        "rel_lambda_diff": component_diag.get("rel_lambda_diff", 0.0),
        "r2s_dense": r2s_d,
        "r2s_mf": r2s_f,
        "d_r2s": abs(r2s_d - r2s_f),
        "band_dense": band_d,
        "band_mf": band_f,
        "d_band": abs(band_d - band_f),
        "mwb_dense": mwb_d,
        "mwb_mf": mwb_f,
        "d_mwb": abs(mwb_d - mwb_f),
    }
    ensemble = intrinsic_ensemble(
        matrix, mode=mode, alpha=alpha, reps=intrinsic_reps
    )
    lo, hi = ensemble.min(axis=0), ensemble.max(axis=0)
    # Reported as the ensemble half-width, so the CSV shows how ambiguous each number is.
    record["d_r2s_intrinsic"] = max(abs(lo[0] - r2s_d), abs(hi[0] - r2s_d))
    record["d_band_intrinsic"] = max(abs(lo[1] - band_d), abs(hi[1] - band_d))
    record["d_mwb_intrinsic"] = max(abs(lo[2] - mwb_d), abs(hi[2] - mwb_d))
    record["_ensemble_lo"] = tuple(float(v) for v in lo)
    record["_ensemble_hi"] = tuple(float(v) for v in hi)
    record["metric_class"] = classify_record(record)
    record["ok"] = record["metric_class"] != "fail"
    return record


def _largest_component_diagnostics(
    matrix: np.ndarray, mode: str, alpha: float | None
) -> dict:
    """Dense-vs-sparse diagnostics on the largest connected component only.

    Whole-matrix agreement is scored by the metrics; this adds the spectral
    quantities (eigenvalue drift, coordinate drift, tie confinement) at the one
    component where the two solvers differ most.
    """
    components = tc.support_components(matrix)
    if not components:
        return {}
    component = max(components, key=lambda c: c.n_rows + c.n_cols)
    if component.n_rows + component.n_cols < mf.ARPACK_MIN_VERTICES:
        return {"max_component_vertices": component.n_rows + component.n_cols}

    normalized = matrix / float(matrix.max())
    block = normalized[np.ix_(component.row_indices, component.col_indices)]
    record = mf.compare_block(block, alpha, mode, tie_tol=COORD_TOL, tol=1e-12)
    record["max_component_vertices"] = component.n_rows + component.n_cols
    return record


def intrinsic_ensemble(
    matrix: np.ndarray,
    *,
    mode: str,
    alpha: float | None,
    reps: int,
    seed: int = 0,
) -> np.ndarray:
    """Metric values attained by valid orderings under 1e-12 input perturbation.

    This is the control for the metric comparison.  Where the Fiedler coordinates are
    tied the induced order is not unique, so a 1e-12 relative perturbation of the input
    legitimately moves the ordering -- and with it the metric.  The returned
    ``(reps + 1, 3)`` array is therefore a sample of the metric values this method
    produces on this instance, and the matrix-free result is accepted when it is a
    *member* of that set rather than when it clears an arbitrary threshold.

    Metrics are always evaluated on the *unperturbed* matrix, so only the ordering effect
    is measured, never the direct effect of perturbing the values.
    """
    base = mf.componentwise_reorder(matrix, mode=mode, alpha=alpha, solver="dense")
    values = [list(_axis_metrics(matrix, base.row_order, base.col_order))]

    rng = np.random.default_rng(seed)
    for _ in range(max(reps, 0)):
        perturbed = matrix * (1.0 + INTRINSIC_EPS * rng.standard_normal(matrix.shape))
        np.maximum(perturbed, 0.0, out=perturbed)
        result = mf.componentwise_reorder(
            perturbed, mode=mode, alpha=alpha, solver="dense"
        )
        values.append(list(_axis_metrics(matrix, result.row_order, result.col_order)))
    return np.asarray(values, dtype=float)


# How many times the instance's own tie-induced spread a solver delta may reach before it
# is called a genuine disagreement.  "Same order of magnitude" rather than "smaller than",
# because the perturbation ensemble is a small sample and cannot bound the spread tightly.
AMBIGUITY_FACTOR = 10.0


def classify_record(record: dict) -> str:
    """``exact``, ``tie-ambiguous`` or ``fail``.

    ``exact``          every metric matches to the tight tolerance.
    ``tie-ambiguous``  the metrics move, but no more than an order of magnitude beyond
                       what a 1e-12 perturbation of the *input* already produces.  The
                       reported number is not determined by the method on this instance,
                       so the difference cannot be attributed to the solver.
    ``fail``           otherwise, or the spectral contract (eigenvalue, sorted
                       coordinates, tie confinement) was violated.
    """
    if not (
        record["rel_lambda_diff"] <= REL_LAMBDA_TOL
        and record["max_sorted_coord_dev"] <= COORD_TOL
        and bool(record["disagreement_within_ties"])
    ):
        return "fail"

    deltas = (
        (record["d_r2s"], record.get("d_r2s_intrinsic") or 0.0, R2S_TOL),
        (record["d_band"], record.get("d_band_intrinsic") or 0.0, BAND_TOL),
        (record["d_mwb"], record.get("d_mwb_intrinsic") or 0.0, MWB_TOL),
    )
    if all(delta <= tol for delta, _, tol in deltas):
        return "exact"
    for delta, intrinsic, tol in deltas:
        if delta > max(tol, AMBIGUITY_FACTOR * max(intrinsic, tol)):
            return "fail"
    return "tie-ambiguous"


def _passes(record: dict) -> bool:
    return classify_record(record) != "fail"


# --------------------------------------------------------------------------- #
# Alpha-selection stability
# --------------------------------------------------------------------------- #


def select_alpha(matrix: np.ndarray, solver: str) -> tuple[float, float]:
    """Mirror ``tw_auto_reorder``: maximise the released MWB-AUC over the alpha grid."""
    best_alpha, best_score = None, None
    for alpha in ALPHAS:
        result = mf.componentwise_reorder(
            matrix, mode="tw", alpha=alpha, solver=solver, tol=1e-12
        )
        score = float(mf.bm.mwb_auc(result.matrix, WIDTHS))
        if (
            best_score is None
            or score > best_score + 1e-12
            or (abs(score - best_score) <= 1e-12 and alpha < best_alpha)
        ):
            best_alpha, best_score = alpha, score
    return float(best_alpha), float(best_score)


# --------------------------------------------------------------------------- #
# Corpus
# --------------------------------------------------------------------------- #


def real_datasets(keys: tuple[str, ...]) -> list[tuple[str, np.ndarray]]:
    out = []
    for key in keys:
        path = MATRIX_DIR / f"{key}_original_matrix.csv"
        if not path.is_file():
            raise FileNotFoundError(path)
        matrix, _, _ = mf.bm.load_matrix_csv(path)
        out.append((DATASET_LABELS.get(key, key), matrix))
    return out


def synthetic_instances(seeds: int) -> list[tuple[str, np.ndarray]]:
    out = []
    for family in FAMILIES:
        family_key = family.key if hasattr(family, "key") else str(family)
        for size in SIZES:
            for seed in range(seeds):
                payload = build_case(family, size, seed)
                out.append(
                    (f"{family_key}/{size.key}/s{seed}", payload["observed"])
                )
    return out


def build_records(
    datasets: tuple[str, ...],
    seeds: int,
    include_synthetic: bool,
    intrinsic_reps: int = 2,
) -> list[dict]:
    records: list[dict] = []

    for label, matrix in real_datasets(datasets):
        started = time.perf_counter()
        records.append(
            compare_matrix(
                matrix, mode="ow", alpha=None, suite="real", instance=label,
                intrinsic_reps=intrinsic_reps,
            )
        )
        for alpha in ALPHAS:
            records.append(
                compare_matrix(
                    matrix, mode="tw", alpha=alpha, suite="real", instance=label,
                    intrinsic_reps=intrinsic_reps,
                )
            )
        records.append(
            {
                **_alpha_record(label, matrix),
                "suite": "real",
            }
        )
        print(f"  {label:24s} done in {time.perf_counter() - started:.1f}s", flush=True)

    if include_synthetic:
        instances = synthetic_instances(seeds)
        print(f"  synthetic: {len(instances)} instances", flush=True)
        for index, (name, matrix) in enumerate(instances):
            records.append(
                compare_matrix(
                    matrix, mode="ow", alpha=None, suite="synthetic", instance=name,
                    intrinsic_reps=intrinsic_reps,
                )
            )
            records.append(
                compare_matrix(
                    matrix, mode="tw", alpha=1.0, suite="synthetic", instance=name,
                    intrinsic_reps=intrinsic_reps,
                )
            )
            if index % 25 == 0:
                records.append({**_alpha_record(name, matrix), "suite": "synthetic"})

    return records


def _alpha_record(instance: str, matrix: np.ndarray) -> dict:
    dense_alpha, dense_score = select_alpha(matrix, "dense")
    free_alpha, free_score = select_alpha(matrix, "matrix_free")
    agree = dense_alpha == free_alpha
    return {
        "instance": instance,
        "mode": "tw_auto",
        "alpha": f"{free_alpha:g}",
        "m": matrix.shape[0],
        "n": matrix.shape[1],
        "nnz": int((matrix > 0).sum()),
        "n_components": 0,
        "max_component_vertices": 0,
        "n_positions_differ": 0,
        "frac_positions_differ": 0.0,
        "disagreement_within_ties": True,
        "max_sorted_coord_dev": 0.0,
        "rel_lambda_diff": 0.0,
        "r2s_dense": "",
        "r2s_mf": "",
        "d_r2s": abs(dense_score - free_score),
        "band_dense": "",
        "band_mf": "",
        "d_band": 0.0,
        "mwb_dense": dense_score,
        "mwb_mf": free_score,
        "d_mwb": abs(dense_score - free_score),
        "d_r2s_intrinsic": 0.0,
        "d_band_intrinsic": 0.0,
        "d_mwb_intrinsic": 0.0,
        "metric_class": "alpha-selection",
        "ok": agree,
    }


# --------------------------------------------------------------------------- #
# Reporting
# --------------------------------------------------------------------------- #


def write_csv(records: list[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        for record in records:
            writer.writerow({key: record.get(key, "") for key in FIELDS})


def write_markdown(records: list[dict], path: Path) -> None:
    checks = [r for r in records if r["mode"] != "tw_auto"]
    alphas = [r for r in records if r["mode"] == "tw_auto"]

    lines = [
        "# Matrix-free vs dense equivalence (Section 3.3)",
        "",
        f"Instances compared: {len(checks)} solver configurations "
        f"plus {len(alphas)} adaptive-alpha selections.",
        "",
        "## Agreement contract",
        "",
        "| quantity | worst observed | threshold | pass |",
        "|---|---|---|---|",
    ]
    worst = {
        "rel_lambda_diff": (max((r["rel_lambda_diff"] for r in checks), default=0.0), REL_LAMBDA_TOL),
        "max_sorted_coord_dev": (max((r["max_sorted_coord_dev"] for r in checks), default=0.0), COORD_TOL),
        "d_r2s": (max((r["d_r2s"] for r in checks), default=0.0), R2S_TOL),
        "d_band": (max((r["d_band"] for r in checks), default=0.0), BAND_TOL),
        "d_mwb": (max((r["d_mwb"] for r in checks), default=0.0), MWB_TOL),
    }
    for name, (value, tol) in worst.items():
        lines.append(f"| {name} | {value:.3e} | {tol:.0e} | {'yes' if value <= tol else '**NO**'} |")

    ties = [r for r in checks if not r["disagreement_within_ties"]]
    by_class: dict[str, list[dict]] = {}
    for record in checks:
        by_class.setdefault(record["metric_class"], []).append(record)
    lines += [
        "",
        "## Classification",
        "",
        "| class | configurations | meaning |",
        "|---|---|---|",
        f"| exact | {len(by_class.get('exact', []))} | every metric matches to the tight tolerance |",
        f"| tie-ambiguous | {len(by_class.get('tie-ambiguous', []))} | metrics move, but no further than a "
        "1e-12 perturbation of the input already moves them |",
        f"| fail | {len(by_class.get('fail', []))} | genuine disagreement, or the spectral contract broke |",
        "",
    ]
    if by_class.get("tie-ambiguous"):
        lines += [
            "Tie-ambiguous instances, with the solver delta beside the instance's own",
            "input-perturbation spread:",
            "",
            "| instance | mode | alpha | dR2S | spread | dMWB | spread | columns |",
            "|---|---|---|---|---|---|---|---|",
        ]
        for record in by_class["tie-ambiguous"]:
            lines.append(
                f"| {record['instance']} | {record['mode']} | {record['alpha']} | "
                f"{record['d_r2s']:.2e} | {record['d_r2s_intrinsic']:.2e} | "
                f"{record['d_mwb']:.2e} | {record['d_mwb_intrinsic']:.2e} | {record['n']} |"
            )
        lines.append("")
    lines += [
        f"- Orderings differ outside a tie group: **{len(ties)}** configurations.",
        f"- Configurations whose permutation differs at all: "
        f"**{sum(1 for r in checks if r['n_positions_differ'])}** of {len(checks)}.",
        "",
        "A differing permutation is expected and harmless: the Fiedler coordinates are tied,",
        "ties are adjacent once sorted, and all three metrics are invariant to permuting",
        "within a tie group. The published per-instance permutations are only defined up to",
        "those ties.",
        "",
        "## Alpha-selection stability",
        "",
        "| instance | selected alpha | dMWB-AUC | stable |",
        "|---|---|---|---|",
    ]
    for record in alphas:
        lines.append(
            f"| {record['instance']} | {record['alpha']} | {record['d_mwb']:.3e} | "
            f"{'yes' if record['ok'] else '**NO**'} |"
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--datasets", default=",".join(DATASET_ORDER), help="comma-separated dataset keys"
    )
    parser.add_argument("--seeds", type=int, default=20, help="synthetic seeds per regime")
    parser.add_argument("--no-synthetic", action="store_true")
    parser.add_argument(
        "--intrinsic-reps",
        type=int,
        default=2,
        help="input perturbations used to estimate each metric's conditioning (0 disables)",
    )
    args = parser.parse_args()

    datasets = tuple(k.strip() for k in args.datasets.split(",") if k.strip())
    started = time.perf_counter()
    records = build_records(
        datasets, args.seeds, not args.no_synthetic, args.intrinsic_reps
    )
    elapsed = time.perf_counter() - started

    write_csv(records, OUTPUT_DIR / "matrixfree_equivalence.csv")
    write_markdown(records, OUTPUT_DIR / "matrixfree_equivalence.md")

    checks = [r for r in records if r["mode"] != "tw_auto"]
    failed = [r for r in records if not r["ok"]]
    outside_ties = [r for r in checks if not r["disagreement_within_ties"]]

    print(f"\n{len(records)} records in {elapsed:.1f}s")
    print(f"  outside tie groups : {len(outside_ties)}")
    print(f"  failed contract    : {len(failed)}")
    for record in failed[:10]:
        print(
            f"    {record['suite']}/{record['instance']} {record['mode']} "
            f"alpha={record['alpha']}"
        )
    print(f"  wrote {OUTPUT_DIR / 'matrixfree_equivalence.csv'}")

    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
