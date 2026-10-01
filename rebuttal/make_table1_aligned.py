"""Emit Table 1 in the paper's layout, as-submitted and corrected.

Produces markdown tables over the seven real matrices, each with the same
row/column structure as Table 1 of the manuscript (rows = datasets, column
groups = metrics, within each group = the eight methods -- the manuscript's seven
plus Gram-OW, the reviewer XFNF W5 baseline):

  A. As submitted   -- single Fiedler vector for One-walk and Two-walk, released
                       metric definitions. Must reproduce the committed CSV; the
                       script checks this and reports any drift.
  B. Corrected      -- Section 3.3 component-wise ordering for One-walk and
                       Two-walk, metrics under the paper-faithful definitions
                       (Eq. 16 Band@10%, Section 4.5 block-cut MWB-AUC).
  C. Reproduction   -- A against the *submitted* CSV, B against the *current* one.
  D. Alpha          -- the alpha each MWB-AUC selector picks: the legacy band
                       integral the code actually uses, versus the block-cut
                       quantity Section 4.5 describes.

The two tables are checked against different baselines on purpose. A is the code
as submitted, so it reproduces the CSV at the pre-fix revision (SUBMITTED_REV); B
is the corrected code, so it reproduces the working-tree CSV.

Both tables select alpha with the *released* MWB-AUC, because that is what the
repository's adaptive Two-walk actually uses and what produced the committed CSV.

Table A cannot be built from `bm.one_walk_reorder` / `bm.tw_auto_reorder` any
more: the benchmark's paper-spec commit turned those into the Section 3.3
component-wise procedure, so using them would make A identical to B and destroy
the contrast this document exists to show. A comes from
`tw_components.order_single`, the surviving implementation of the pre-fix
behaviour.

Original / Marginal / HC+OLO / CA-SVD / Median do not use the Fiedler pipeline,
so their orderings are unchanged by the correction; only their metric values move
under the corrected definitions. One-walk and Two-walk change in both respects.
Gram-OW is new here rather than corrected: it was never in the submitted table, and
it keeps one ordering across A and B like the other non-adaptive methods.

Run:
    cd MHeatMap/src
    uv run python rebuttal/make_table1_aligned.py
"""

from __future__ import annotations

import argparse
import csv
import glob
import io
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import baselines_v2 as bv  # noqa: E402
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
# The benchmark as submitted, before the paper-spec metric and ordering fixes.
# Table A is that code, so that is the CSV it should reproduce.
SUBMITTED_REV = "445334d"

# `Gram-OW` is the reviewer XFNF W5 baseline (see `baselines_v2`). It sits next to
# One-walk and Two-walk because that is the comparison it exists to serve.
METHODS = ("Original", "Marginal", "HC+OLO", "One-walk", "Gram-OW", "CA-SVD", "Median", "TW")
METHOD_SHORT = ("O", "M", "HC", "OW", "GO", "CA", "MD", "TW")
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
    """The orderings shared by tables A and B.

    The five non-spectral baselines plus Gram-OW. None of them goes through the
    adaptive Two-walk path, so the single-vector -> component-wise correction leaves
    their orderings alone and only their metric values move between A and B.
    """
    n_rows, n_cols = matrix.shape
    return {
        "Original": (np.arange(n_rows), np.arange(n_cols)),
        "Marginal": _orders(tc.bm.marginal_sort_reorder(matrix)),
        "HC+OLO": _orders(tc.bm.hierarchical_olo_reorder(matrix)),
        "Gram-OW": bv.gram_one_walk_order(matrix),
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
        # `tc.order_single` is the single-Fiedler baseline of record. It must not
        # be `bm.one_walk_reorder` / `bm.tw_auto_reorder`: those are now the
        # Section 3.3 component-wise procedure, i.e. Table B.
        submitted_alpha, tw_submitted = tw_adaptive(matrix, tc.order_single, _select_released)
        submitted_orderings = dict(baseline)
        submitted_orderings["One-walk"] = tc.order_single(matrix, 1.0, "ow")
        submitted_orderings["TW"] = tw_submitted

        as_submitted[label] = {
            method: metric_triple(matrix, *submitted_orderings[method], "released")
            for method in METHODS
        }
        alphas.setdefault(label, {})["as submitted"] = submitted_alpha

        # --- B: corrected ----------------------------------------------------
        # Alpha is selected by the released MWB-AUC in both tables, matching the
        # adaptive procedure the repository actually runs.
        best_alpha, tw_comp = tw_adaptive(matrix, tc.order_componentwise, _select_released)
        corrected_orderings = dict(baseline)
        corrected_orderings["One-walk"] = tc.order_componentwise(matrix, 0.0, "ow")
        corrected_orderings["TW"] = tw_comp

        corrected[label] = {
            method: metric_triple(matrix, *corrected_orderings[method], "faithful")
            for method in METHODS
        }
        alphas[label]["corrected"] = best_alpha

        # --- D: what the paper-as-written selector would have chosen ---------
        block_alpha, _ = tw_adaptive(matrix, tc.order_componentwise, _select_block)
        alphas[label]["block-cut"] = block_alpha

    return as_submitted, corrected, alphas


def _table(title: str, rows: dict, alphas: dict, note: str, alpha_key: str) -> list[str]:
    metrics = ("2-SUM", "Band@10%", "MWB-AUC")
    header = "| Dataset | Shape | " + " | ".join(
        f"{m} {short}" for m in metrics for short in METHOD_SHORT
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


def _alpha_section(alphas: dict) -> list[str]:
    """Section D: alpha under the legacy selector versus the block-cut one.

    Kept compact on purpose -- the point is whether the choice moves, not to be a
    third full metric table.
    """
    def fmt(value) -> str:
        return "-" if value is None else f"{value:g}"

    lines = [
        "## D. Alpha selection under each MWB-AUC definition",
        "",
        "Section 4.5 describes MWB-AUC as a matched contiguous block-cut quantity, "
        "but the repository's adaptive Two-walk selects alpha with the legacy "
        "integer-window band integral. This shows whether adopting the "
        "paper-as-written selector would have changed the choice.",
        "",
        "| Dataset | As submitted | Corrected, legacy selector | Corrected, block-cut selector |",
        "|---|---:|---:|---:|",
    ]
    for key in DATASET_ORDER:
        label = DATASET_LABELS[key]
        if label not in alphas:
            continue
        row = alphas[label]
        lines.append(
            f"| {label} | {fmt(row.get('as submitted'))} "
            f"| {fmt(row.get('corrected'))} | {fmt(row.get('block-cut'))} |"
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

    compare = _reproduction_check(as_submitted, corrected)

    lines = [
        "# Table 1 -- aligned with the manuscript",
        "",
        "Column groups are the three metrics; within each group the method order",
        "matches Table 1 with one addition: O = Original, M = Marginal, HC = HC+OLO,",
        "OW = One-walk, GO = Gram-OW, CA = CA-SVD, MD = Median, TW = Two-walk. Bold",
        "marks the best value per row and metric. Row order follows Table 1.",
        "",
        "`Gram-OW` is the baseline reviewer XFNF W5 asked for: one-walk Laplacian",
        "reordering applied separately to `B B^T` (row order) and `B^T B` (column",
        "order). Because the two axes are solved as independent eigenproblems, their",
        "relative reversal is a free choice that the joint One-walk/Two-walk",
        "eigenvector cannot express. It is resolved here with the released",
        "`orient_orders_for_diagonal` rule -- try four reversals, keep the lowest",
        "2-SUM -- which is deliberately *more* generous than the treatment OW and TW",
        "receive in table B, where that metric-fitted step was removed. Gram-OW is",
        "therefore, if anything, flattered by this table. See `baselines_v2.py`.",
        "",
    ]
    lines += _table(
        "A. As submitted",
        as_submitted,
        alphas,
        "The submitted code: one Fiedler vector over the whole active support, "
        "finished with `orient_orders_for_diagonal` and an unstable sort, under the "
        f"released metric definitions. Should reproduce the `{SUBMITTED_REV}` table.",
        "as submitted",
    )
    lines += _table(
        "B. Corrected",
        corrected,
        alphas,
        "Section 3.3 component-wise ordering for OW and TW, alpha chosen by the "
        "released MWB-AUC (the selector the repository actually runs), and the "
        "paper-faithful metric definitions (Eq. 16 Band@10%, Section 4.5 "
        "block-cut MWB-AUC). Baselines keep their orderings; only their metric "
        "values move.",
        "corrected",
    )
    lines += ["## C. Reproduction check", ""]
    lines += compare
    lines += _alpha_section(alphas)

    path = OUTPUT_DIR / "table1_aligned.md"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"\nWrote {path}")


def _read_csv_rows(text: str) -> dict[str, dict[str, str]]:
    return {row["dataset"].strip(): row for row in csv.DictReader(io.StringIO(text))}


def _submitted_rows() -> dict[str, dict[str, str]] | None:
    """The submitted CSV, read from git. None when that revision is unavailable."""
    rel = COMMITTED_CSV.relative_to(tc.REPO_ROOT).as_posix()
    try:
        proc = subprocess.run(
            ["git", "show", f"{SUBMITTED_REV}:{rel}"],
            cwd=tc.REPO_ROOT,
            capture_output=True,
            text=True,
            check=True,
        )
    except (subprocess.CalledProcessError, FileNotFoundError, OSError):
        return None
    return _read_csv_rows(proc.stdout)


def _compare(rows: dict, committed: dict, columns, tolerance: float = 5e-3) -> list[str]:
    """Cells of `rows` that differ from `committed` beyond `tolerance`."""
    drift = []
    for label, methods in rows.items():
        # CSV keys differ in punctuation from the table labels; match the lead token.
        lead = label.split()[0]
        match = next((k for k in committed if k.split()[0] == lead), None)
        if match is None:
            drift.append(f"- {label}: no matching row in the reference CSV")
            continue
        for method, metric, column in columns:
            mine = methods[method][metric]
            theirs = float(committed[match][column])
            if abs(mine - theirs) > tolerance:
                drift.append(
                    f"- {label} | {method} | {metric}: ours {mine:.4f} "
                    f"vs reference {theirs:.4f}"
                )
    return drift


# A is the submitted code under the released metric definitions, so all three
# columns are comparable with the submitted CSV. B is the corrected code under the
# paper-faithful definitions, so its MWB-AUC column is the block-cut quantity and
# is deliberately not compared.
_SUBMITTED_COLUMNS = tuple(
    (method, metric, column)
    for method, key in (("One-walk", "One-walk"), ("TW", "TW"))
    for metric, column in (
        ("2-SUM", f"two_sum|{key}"),
        ("Band@10%", f"band_mass_10|{key}"),
        ("MWB-AUC", f"mwb_auc|{key}"),
    )
)
_CORRECTED_COLUMNS = tuple(
    (method, metric, column)
    for method, key in (("One-walk", "One-walk"), ("TW", "TW"))
    for metric, column in (
        ("2-SUM", f"two_sum|{key}"),
        ("Band@10%", f"band_mass_10|{key}"),
    )
)


def _reproduction_check(as_submitted: dict, corrected: dict) -> list[str]:
    """A against the submitted CSV, B against the working-tree CSV."""
    lines: list[str] = []

    lines += [
        f"### A. As submitted vs the `{SUBMITTED_REV}` CSV",
        "",
        "The submitted code should reproduce its own table. All seven datasets are",
        "checked; a mismatch on One-walk or Two-walk means the single-vector",
        "baseline here has drifted from the released behaviour.",
        "",
    ]
    submitted = _submitted_rows()
    if submitted is None:
        lines += [f"Revision `{SUBMITTED_REV}` not available; skipped.", ""]
    else:
        drift = _compare(as_submitted, submitted, _SUBMITTED_COLUMNS)
        lines += (
            ["All checked cells reproduce the submitted CSV to within 5e-3.", ""]
            if not drift
            else ["Cells that differ:", ""] + drift + [""]
        )

    lines += [
        "### B. Corrected vs the current CSV",
        "",
        "The corrected ordering pipeline should reproduce the repository's current",
        "table exactly on the shared columns. MWB-AUC is excluded because Table B",
        "reports the Section 4.5 block-cut quantity while the CSV carries the",
        "selector's band integral.",
        "",
    ]
    if not COMMITTED_CSV.exists():
        lines += ["Working-tree CSV not found; skipped.", ""]
    else:
        current = _read_csv_rows(COMMITTED_CSV.read_text(encoding="utf-8"))
        drift = _compare(corrected, current, _CORRECTED_COLUMNS)
        lines += (
            ["All checked cells reproduce the current CSV to within 5e-3.", ""]
            if not drift
            else ["Cells that differ:", ""] + drift + [""]
        )
    return lines


if __name__ == "__main__":
    main()
