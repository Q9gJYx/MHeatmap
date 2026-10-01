"""D4 -- component-wise ordering audit (reviewer kFbK W2 / XFNF W2).

Compares the pre-fix single-Fiedler-vector procedure against the Section 3.3
component-wise procedure on the real benchmarks, under both One-Walk and Two-Walk,
over an alpha grid, reporting every metric under both the released and the
paper-faithful definition.

The two procedures differ by *three* independent changes, so a gap here does not
isolate the component correction:

1. component-wise ordering -- the subject of A7, and the only change reviewers
   asked for;
2. the final orientation step. The pre-fix path ended with
   `orient_orders_for_diagonal`, which searched four row/column reversals and kept
   the lowest 2-SUM -- the metric being reported. The paper-spec commit dropped
   that from One-Walk and Two-Walk;
3. the `argsort` tie-break, quicksort before and `kind="stable"` now. On MBTA 100
   of the 133 Fiedler entries are tied, so this alone reorders the result.

Changes 2 and 3 are **paper-conformance fixes, not arbitrary edits**. Section 3.3
already specifies the new behaviour verbatim: "we orient the vector so that its
first coordinate of maximum absolute value is positive, then apply a stable
sort". The submitted code did neither -- it searched four reversals for the lowest
2-SUM and used an unstable sort -- so the submitted Table 1 numbers were produced
by code that did not implement the procedure the paper describes. That is the same
class of defect as the Band@10% and MWB-AUC mismatches, and unlike those it is not
confined to disconnected support.

Because changes 2 and 3 apply to connected matrices, they move the published
One-Walk and Two-Walk numbers on ACS, LODES, 20 Newsgroups and MBTA as well as on
the three disconnected ones. The self-check below therefore validates the
component plumbing directly rather than asserting the two procedures agree, which
is no longer true by design.

Run:
    cd MHeatMap/src
    uv run python rebuttal/run_components_audit.py                # all datasets, full grid
    uv run python rebuttal/run_components_audit.py --datasets openalex --alphas 2 8

Outputs (under rebuttal/output/):
    component_audit.csv   machine-readable, one row per dataset x method x alpha
    component_audit.md    paper-style tables, one per metric family
    components.csv        per-component structure per dataset
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

DATASET_LABELS = {
    "naics_sic": "SIC -> NAICS",
    "cip_soc": "CIP -> SOC",
    "acs_clean": "ACS OCCP x INDP",
    "lodes_clean": "LODES Home x Work",
    "twenty_newsgroups": "20 Newsgroups",
    "gtfs_clean": "MBTA Route x Station",
    "openalex": "OpenAlex Author x Topic",
}
DEFAULT_ALPHAS = (1.0, 2.0, 4.0, 6.0, 8.0, 12.0)

def load_datasets(keys: list[str]) -> dict[str, np.ndarray]:
    out: dict[str, np.ndarray] = {}
    for key in keys:
        matches = glob.glob(str(MATRIX_DIR / f"{key}_original_matrix.csv"))
        if not matches:
            print(f"  [skip] no cached matrix for {key}", file=sys.stderr)
            continue
        matrix, _, _ = tc.bm.load_matrix_csv(Path(matches[0]))
        out[key] = matrix
    return out


def measure(matrix: np.ndarray, row_order, col_order) -> dict[str, float]:
    """All five metrics for one ordering."""
    ordered = tc.bm.apply_orders(matrix, row_order, col_order)
    return {
        "r2s": mv.r2s(matrix, row_order, col_order),
        "band_released": mv.released_band_mass(tc.bm, ordered),
        "band_eq16": mv.band_mass_eq16(matrix, row_order, col_order),
        "mwb_band": mv.released_mwb_auc(tc.bm, ordered),
        "mwb_block": mv.mwb_auc_block(matrix, row_order, col_order),
    }


def evaluate(matrix: np.ndarray, alpha: float | None, mode: str) -> dict[str, float]:
    """Both procedures under every metric, for one (alpha, mode) configuration."""
    effective_alpha = 0.0 if alpha is None else alpha
    single_rows, single_cols = tc.order_single(matrix, effective_alpha, mode)
    comp_rows, comp_cols = tc.order_componentwise(matrix, effective_alpha, mode)
    row = {}
    for key, value in measure(matrix, single_rows, single_cols).items():
        row[f"single_{key}"] = value
    for key, value in measure(matrix, comp_rows, comp_cols).items():
        row[f"comp_{key}"] = value
    return row


def _component_plumbing_ok(matrix: np.ndarray, components: list) -> str | None:
    """Return a description of the first fault, or None when the machinery is sound.

    On a single-component matrix the component path has nothing to assemble, so it
    must reduce to one spectral solve over exactly the full active support. This
    test used to be written as "the two procedures return identical orderings"; that
    stopped holding when the paper-spec commit also changed the orientation step and
    the argsort tie-break, so it is now checked directly.
    """
    if len(components) != 1:
        return f"expected a single component in the active support, found {len(components)}"

    values = np.asarray(matrix, dtype=float)
    active_rows = np.flatnonzero(values.sum(axis=1) > 0)
    active_cols = np.flatnonzero(values.sum(axis=0) > 0)
    component = components[0]
    if not np.array_equal(component.row_indices, active_rows):
        return "the single component does not cover the active rows"
    if not np.array_equal(component.col_indices, active_cols):
        return "the single component does not cover the active columns"

    n_rows, n_cols = values.shape
    for rows, cols in (
        tc.order_single(matrix, 1.0, "ow"),
        tc.order_componentwise(matrix, 0.0, "ow"),
    ):
        if not np.array_equal(np.sort(rows), np.arange(n_rows)):
            return "One-walk returned a malformed row order"
        if not np.array_equal(np.sort(cols), np.arange(n_cols)):
            return "One-walk returned a malformed column order"
    return None


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--datasets", nargs="*", default=list(DATASET_LABELS))
    parser.add_argument("--alphas", nargs="*", type=float, default=list(DEFAULT_ALPHAS))
    args = parser.parse_args()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    datasets = load_datasets(args.datasets)
    if not datasets:
        raise SystemExit("no datasets loaded; check --datasets")

    records: list[dict] = []
    component_rows: list[dict] = []
    failures: list[str] = []
    checks_run = 0

    for key, matrix in datasets.items():
        label = DATASET_LABELS.get(key, key)
        n_rows, n_cols = matrix.shape
        components = tc.support_components(matrix)
        for entry in tc.component_table(matrix):
            component_rows.append({"dataset": label, **entry})

        single_component = len(components) == 1
        print(f"\n=== {label} ({n_rows}x{n_cols}) -- {len(components)} component(s)")

        if single_component:
            checks_run += 1
            problem = _component_plumbing_ok(matrix, components)
            if problem is not None:
                failures.append(f"{label} | {problem}")

        configurations: list[tuple[str, float | None]] = [("One-walk", None)]
        configurations += [("Two-walk", alpha) for alpha in args.alphas]

        for method, alpha in configurations:
            mode = "ow" if method == "One-walk" else "tw"
            started = time.perf_counter()
            row = evaluate(matrix, alpha, mode)
            seconds = time.perf_counter() - started
            records.append(
                {
                    "dataset": label,
                    "method": method,
                    "alpha": "" if alpha is None else alpha,
                    **row,
                    "seconds": seconds,
                }
            )


            alpha_text = "   -" if alpha is None else f"{alpha:4g}"
            print(
                f"  {method:9s} a={alpha_text}  R2S {row['single_r2s']:6.3f} -> "
                f"{row['comp_r2s']:6.3f}  | Band16 {row['single_band_eq16']:.3f} -> "
                f"{row['comp_band_eq16']:.3f}  | MWBblk {row['single_mwb_block']:.3f} -> "
                f"{row['comp_mwb_block']:.3f}"
            )

    _write_csv(OUTPUT_DIR / "component_audit.csv", records)
    _write_components_csv(OUTPUT_DIR / "components.csv", component_rows)
    _write_markdown(
        OUTPUT_DIR / "component_audit.md", records, component_rows, checks_run, failures
    )

    print(f"\nWrote results to {OUTPUT_DIR}")
    for name in ("component_audit.csv", "component_audit.md", "components.csv"):
        print(f"  {name}")
    print(f"\nSelf-check: {checks_run} single-component datasets checked.")
    if failures:
        print("  FAILURES:")
        for line in failures:
            print(f"    {line}")
    else:
        print("  Component machinery sound on all of them.")


def _write_csv(path: Path, records: list[dict]) -> None:
    if not records:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def _write_components_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _metric_table(
    records: list[dict], metric: str, better: str, scale: float = 1.0
) -> list[str]:
    lines = [
        f"*{better} is better.*",
        "",
        "| Dataset | Method | alpha | single | comp | delta |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for row in records:
        alpha = row["alpha"]
        alpha_text = "-" if alpha == "" else f"{alpha:g}"
        single = row[f"single_{metric}"] * scale
        comp = row[f"comp_{metric}"] * scale
        lines.append(
            f"| {row['dataset']} | {row['method']} | {alpha_text} "
            f"| {single:.3f} | {comp:.3f} | {comp - single:+.3f} |"
        )
    return lines


def _write_markdown(
    path: Path,
    records: list[dict],
    component_rows: list[dict],
    checks_run: int,
    failures: list[str],
) -> None:
    lines = [
        "# Component-wise ordering audit (Section 3.3)",
        "",
        "- `single` = pre-fix procedure: one Fiedler vector over the whole active",
        "  support, finished with the released `orient_orders_for_diagonal` step.",
        "- `comp` = Section 3.3 procedure, via the benchmark's canonical",
        "  implementation: order each component, arrange by decreasing mass.",
        "- `delta` = comp - single. Negative is an improvement for R2S only.",
        "",
        "**Three changes, not one.** The two procedures differ by the component",
        "correction *and* by two fixes the paper-spec commit bundled with it:",
        "",
        "| # | Change | Datasets affected |",
        "|---|---|---|",
        "| 1 | Component-wise ordering (Section 3.3) | SIC -> NAICS, CIP -> SOC, OpenAlex |",
        "| 2 | Orientation: the 2-SUM search in `orient_orders_for_diagonal` dropped | all |",
        "| 3 | `argsort` quicksort -> `kind=\"stable\"` | 20 Newsgroups, MBTA |",
        "",
        "Changes 2 and 3 are paper-conformance fixes -- Section 3.3 specifies",
        '"orient the vector so that its first coordinate of maximum absolute value',
        'is positive, then apply a stable sort", and the submitted code did neither.',
        "So a nonzero delta on a **single-component** dataset is expected, and on",
        "those four datasets it is entirely changes 2 and 3. On MBTA 100 of the 133",
        "Fiedler entries are tied, so change 3 alone reorders the result.",
        "",
        "**Metric definitions.** `band_released` and `mwb_band` reproduce the released",
        "implementation and therefore Table 1. `band_eq16` (Eq. 16) and `mwb_block`",
        "(Section 4.5 matched block cuts) are the paper-faithful definitions. They",
        "differ from the released ones in general -- see reviewer kFbK W1.",
        "",
        "## Component structure",
        "",
        "| Dataset | Rank | Rows | Cols | Mass | Mass fraction |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in component_rows:
        if row["rank"] > 8:
            continue
        lines.append(
            f"| {row['dataset']} | {row['rank']} | {row['n_rows']} | {row['n_cols']} "
            f"| {row['mass']:.1f} | {row['mass_fraction']:.1%} |"
        )

    lines += ["", "## R2S (x100, lower is better)", ""]
    lines += _metric_table(records, "r2s", "Lower")

    lines += ["", "## Band@10% (higher is better)", ""]
    lines += ["**Released integer-window definition (matches Table 1):**", ""]
    lines += _metric_table(records, "band_released", "Higher")
    lines += ["", "**Eq. (16) definition (paper-faithful):**", ""]
    lines += _metric_table(records, "band_eq16", "Higher")

    lines += ["", "## MWB-AUC (higher is better)", ""]
    lines += ["**Released band-mass integral (matches Table 1):**", ""]
    lines += _metric_table(records, "mwb_band", "Higher")
    lines += ["", "**Section 4.5 matched block cuts (paper-faithful):**", ""]
    lines += _metric_table(records, "mwb_block", "Higher")

    lines += ["", "## Self-check", ""]
    lines.append(
        f"{checks_run} single-component datasets were checked. With one component "
        "there is nothing to assemble, so the component path must reduce to a "
        "single solve over exactly the full active support and return a valid "
        "permutation.",
    )
    lines.append(
        "",
    )
    lines.append(
        "This does **not** assert that the two procedures agree on those datasets, "
        "because they no longer do: the paper-spec commit also dropped the "
        "2-SUM orientation step and switched `argsort` to a stable sort. See the "
        "module docstring.",
    )
    if failures:
        lines.append("")
        lines.append("**Failures:**")
        lines += [f"- {line}" for line in failures]
    else:
        lines.append("")
        lines.append("Component machinery is sound on all of them.")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
