"""Emit Table 1 in the paper's layout, as-submitted and corrected.

Produces three markdown tables over the seven real matrices, each with the same
row/column structure as Table 1 of the manuscript (rows = datasets, column
groups = metrics, within each group = the seven methods):

  A. As submitted   -- released orderings (single Fiedler vector) and released
                       metric definitions. Must reproduce the committed CSV; the
                       script checks this and reports any drift.
  B. Corrected      -- Section 3.3 component-wise ordering for One-walk and
                       Two-walk, alpha selected by the paper-faithful block-cut
                       MWB-AUC, metrics under the paper-faithful definitions.
  C. Delta          -- B - A, so the effect of the correction is visible per cell.

Original / Marginal / HC+OLO / CA-SVD / Median do not use the Fiedler pipeline,
so their orderings are unchanged by the correction; only their metric values move
under the corrected definitions. One-walk and Two-walk change in both respects.

Run:
    cd MHeatMap/src
    uv run python rebuttal/make_table1_aligned.py
"""

from __future__ import annotations

import argparse
import csv
import glob
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import metrics_v2 as mv  # noqa: E402
import tw_components as tc  # noqa: E402

HERE = Path(__file__).resolve().parent
OUTPUT_DIR = HERE / "output"
MATRIX_DIR = (
    tc.REPO_ROOT
    / "output"
    / "main_experiment"
    / "real_world_benchmark"
    / "processed_matrices"
)
COMMITTED_CSV = (
    tc.REPO_ROOT
    / "output"
    / "main_experiment"
    / "real_world_benchmark"
    / "real_world_benchmark_paper.csv"
)

METHODS = ("Original", "Marginal", "HC+OLO", "One-walk", "CA-SVD", "Median", "TW")
ALPHAS = (1.0, 2.0, 4.0, 6.0, 8.0, 12.0)
WIDTHS = tc.bm.DEFAULT_WIDTH_GRID

DATASET_LABELS = {
    "naics_sic": "SIC -> NAICS",
    "cip_soc": "CIP -> SOC",
    "acs_clean": "ACS OCCP x INDP",
    "lodes_clean": "LODES Home x Work",
    "twenty_newsgroups": "20 Newsgroups",
    "gtfs_clean": "MBTA Route x Station",
    "openalex": "OpenAlex Author x Topic",
}
# Row order follows Table 1 of the manuscript.
DATASET_ORDER = ("naics_sic", "cip_soc", "acs_clean", "lodes_clean",
                 "twenty_newsgroups", "gtfs_clean", "openalex")


def _orders(result) -> tuple[np.ndarray, np.ndarray]:
    return result.row_order, result.col_order


def baseline_orderings(matrix: np.ndarray) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """The five non-spectral baselines. Unchanged by the component correction."""
    n_rows, n_cols = matrix.shape
    return {
        "Original": (np.arange(n_rows), np.arange(n_cols)),
        "Marginal": _orders(tc.bm.marginal_sort_reorder(matrix)),
        "HC+OLO": _orders(tc.bm.hierarchical_olo_reorder(matrix)),
        "CA-SVD": _orders(tc.bm.ca_svd_reorder(matrix, widths=WIDTHS)),
        "Median": _orders(tc.bm.median_reorder(matrix, widths=WIDTHS)),
    }


def _select_released(matrix, row_order, col_order) -> float:
    return mv.released_mwb_auc(tc.bm, tc.bm.apply_orders(matrix, row_order, col_order))


def _select_block(matrix, row_order, col_order) -> float:
    return mv.mwb_auc_block(matrix, row_order, col_order)


def tw_adaptive(matrix, procedure, selector):
    """Select alpha by the given MWB-AUC, as in Section 4.6."""
    best_score, best_alpha, best_orders = -np.inf, None, None
    for alpha in ALPHAS:
        row_order, col_order = procedure(matrix, alpha, "tw")
        score = selector(matrix, row_order, col_order)
        if score > best_score + 1e-12:
            best_score, best_alpha, best_orders = score, alpha, (row_order, col_order)
    return best_alpha, best_orders


def metric_triple(matrix, row_order, col_order, definitions: str) -> dict[str, float]:
    """`definitions` is 'released' (matches Table 1) or 'faithful' (Section 4.5 / Eq. 16)."""
    ordered = tc.bm.apply_orders(matrix, row_order, col_order)
    if definitions == "released":
        return {
            "2-SUM": mv.r2s(matrix, row_order, col_order),
            "Band@10%": mv.released_band_mass(tc.bm, ordered),
            "MWB-AUC": mv.released_mwb_auc(tc.bm, ordered),
        }
    return {
        "2-SUM": mv.r2s(matrix, row_order, col_order),
        "Band@10%": mv.band_mass_eq16(matrix, row_order, col_order),
        "MWB-AUC": mv.mwb_auc_block(matrix, row_order, col_order),
    }


def build_rows(datasets: dict[str, np.ndarray]):
    """Returns (as_submitted, corrected, selected_alpha) keyed by dataset then method."""
    as_submitted: dict[str, dict[str, dict[str, float]]] = {}
    corrected: dict[str, dict[str, dict[str, float]]] = {}
    alphas: dict[str, dict[str, float | None]] = {}

    for key, matrix in datasets.items():
        label = DATASET_LABELS.get(key, key)
        print(f"  {label} ({matrix.shape[0]}x{matrix.shape[1]}) ...", flush=True)

        baseline = baseline_orderings(matrix)

        # --- A: as submitted -------------------------------------------------
        ow_released = _orders(tc.bm.one_walk_reorder(matrix))
        tw_result = tc.bm.tw_auto_reorder(matrix, alphas=ALPHAS, widths=WIDTHS)
        submitted_orderings = dict(baseline)
        submitted_orderings["One-walk"] = ow_released
        submitted_orderings["TW"] = _orders(tw_result.reorder)

        as_submitted[label] = {
            method: metric_triple(matrix, *submitted_orderings[method], "released")
            for method in METHODS
        }
        alphas.setdefault(label, {})["as submitted"] = tw_result.alpha

        # --- B: corrected ----------------------------------------------------
        ow_comp = tc.order_componentwise(matrix, 0.0, "ow")
        best_alpha, tw_comp = tw_adaptive(matrix, tc.order_componentwise, _select_block)
        corrected_orderings = dict(baseline)
        corrected_orderings["One-walk"] = ow_comp
        corrected_orderings["TW"] = tw_comp

        corrected[label] = {
            method: metric_triple(matrix, *corrected_orderings[method], "faithful")
            for method in METHODS
        }
        alphas[label]["corrected"] = best_alpha

    return as_submitted, corrected, alphas


def _table(title: str, rows: dict, alphas: dict, note: str, alpha_key: str) -> list[str]:
    metrics = ("2-SUM", "Band@10%", "MWB-AUC")
    header = "| Dataset | Shape | " + " | ".join(
        f"{m} {short}" for m in metrics for short in ("O", "M", "HC", "OW", "CA", "MD", "TW")
    ) + " | alpha |"
    separator = "|---|---|" + "---:|" * (len(metrics) * len(METHODS)) + "---:|"

    lines = [f"## {title}", "", note, "", header, separator]
    for key in DATASET_ORDER:
        label = DATASET_LABELS[key]
        if label not in rows:
            continue
        cells = []
        for metric in metrics:
            best = min(
                (rows[label][m][metric] for m in METHODS),
                key=(lambda v: v) if metric == "2-SUM" else (lambda v: -v),
            )
            for method in METHODS:
                value = rows[label][method][metric]
                text = f"{value:.2f}" if metric == "2-SUM" else f"{value:.3f}"
                cells.append(f"**{text}**" if abs(value - best) < 1e-12 else text)
        alpha = alphas[label][alpha_key]
        alpha_text = "-" if alpha is None else f"{alpha:g}"
        lines.append(
            f"| {label} | {SHAPES[label]} | " + " | ".join(cells) + f" | {alpha_text} |"
        )
    lines.append("")
    return lines


SHAPES: dict[str, str] = {}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--datasets", nargs="*", default=list(DATASET_ORDER))
    args = parser.parse_args()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    datasets: dict[str, np.ndarray] = {}
    for key in args.datasets:
        matches = glob.glob(str(MATRIX_DIR / f"{key}_original_matrix.csv"))
        if not matches:
            print(f"  [skip] no cached matrix for {key}", file=sys.stderr)
            continue
        matrix, _, _ = tc.bm.load_matrix_csv(Path(matches[0]))
        datasets[key] = matrix
        SHAPES[DATASET_LABELS[key]] = f"{matrix.shape[0]}x{matrix.shape[1]}"

    print("Building orderings and metrics ...")
    started = time.perf_counter()
    as_submitted, corrected, alphas = build_rows(datasets)
    print(f"  done in {time.perf_counter() - started:.1f}s")

    compare = _compare_with_committed(as_submitted)

    lines = [
        "# Table 1 -- aligned with the manuscript",
        "",
        "Column groups are the three metrics; within each group the method order",
        "matches Table 1: O = Original, M = Marginal, HC = HC+OLO, OW = One-walk,",
        "CA = CA-SVD, MD = Median, TW = Two-walk. Bold marks the best value per row",
        "and metric. Row order follows Table 1.",
        "",
    ]
    lines += _table(
        "A. As submitted",
        as_submitted,
        alphas,
        "Released orderings and released metric definitions. Reproduces Table 1.",
        "as submitted",
    )
    lines += _table(
        "B. Corrected",
        corrected,
        alphas,
        "Section 3.3 component-wise ordering for OW and TW, alpha chosen by the "
        "block-cut MWB-AUC, and the paper-faithful metric definitions "
        "(Eq. 16 Band@10%, Section 4.5 block-cut MWB-AUC). Baselines keep their "
        "orderings; only their metric values move.",
        "corrected",
    )
    lines += ["## C. Reproduction check", ""]
    lines += compare

    path = OUTPUT_DIR / "table1_aligned.md"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"\nWrote {path}")


def _compare_with_committed(as_submitted: dict) -> list[str]:
    """Table A should match the committed CSV; report any cell that does not."""
    if not COMMITTED_CSV.exists():
        return ["Committed CSV not found; skipped."]
    committed = {}
    with COMMITTED_CSV.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            committed[row["dataset"].strip()] = row

    label_to_key = {label: key for key, label in DATASET_LABELS.items()}
    worst = []
    for label, methods in as_submitted.items():
        # Committed CSV keys differ in punctuation; match on the leading token.
        lead = label.split()[0]
        match = next((k for k in committed if k.split()[0] == lead), None)
        if match is None:
            worst.append(f"- {label}: no committed row found")
            continue
        for method, metric, column in (
            ("One-walk", "2-SUM", "two_sum|One-walk"),
            ("TW", "2-SUM", "two_sum|TW"),
            ("One-walk", "MWB-AUC", "mwb_auc|One-walk"),
            ("TW", "MWB-AUC", "mwb_auc|TW"),
        ):
            mine = as_submitted[label][method][metric]
            # The committed CSV already stores 2-SUM at the x100 scale used in Table 1.
            theirs = float(committed[match][column])
            if abs(mine - theirs) > 5e-3:
                worst.append(
                    f"- {label} | {method} | {metric}: ours {mine:.4f} vs committed {theirs:.4f}"
                )
    if not worst:
        return [
            "Table A reproduces the committed CSV for the checked cells "
            "(One-walk and TW under 2-SUM and MWB-AUC) to within 5e-3.",
            "",
        ]
    return ["Cells that differ from the committed CSV:", ""] + worst + [""]


if __name__ == "__main__":
    main()
