"""Emit the synthetic component of Table 1 in the paper's layout.

Same column structure as `table1_aligned.md` -- three metric groups, one column per
method, plus the selected alpha -- but the rows are the 15 family-size regimes of the
synthetic benchmark instead of the seven real matrices. Each row is a mean over the
20 seeds, as in the manuscript. `mt.METHODS` and `mt.METHOD_SHORT` are imported
rather than restated, so the two tables stay in step.

    A. As submitted   -- single Fiedler vector, released metric definitions.
                         Reproduces `synthetic_strong_baselines_paper.csv` at the
                         pre-fix revision (SUBMITTED_REV).
    B. Corrected      -- Section 3.3 component-wise ordering, paper-faithful
                         metrics (Eq. 16 Band@10%, Section 4.5 block-cut MWB-AUC).
                         Reproduces the working-tree CSV on the shared columns.
    C. Reproduction   -- A against the submitted CSV, B against the current one.
    D. Dispersion     -- mean +- sd for One-walk / Gram-OW / Two-walk, and for every
                         method in `synthetic_summary.csv`.
    E. Paired         -- per-seed TW - One-walk, which is the comparative claim.

Every table cell is `mean+-sd` over the 20 seeds; that is the answer to JfQK Q4.
The deviation is the *sample* sd (ddof=1), and for the two adaptive methods it folds
in the alpha the selector picked on each draw as well as the draw itself, so it is a
spread and not a confidence interval on the mean.

The synthetic matrices are all connected on their active support, so the component
correction itself does nothing here. A and B still differ, because the paper-spec
commit also dropped the 2-SUM orientation step and switched `argsort` to a stable
sort -- both of which Section 3.3 already specified. That is why the synthetic rows
move even though A7 does not touch them.

Alpha is selected per instance, so a regime has 20 selected values rather than one.
The alpha column reports the modal value and flags regimes where the seeds did not
agree; `synthetic_per_seed.csv` carries the full per-seed records.

Run:
    cd MHeatMap/src
    uv run python rebuttal/make_table1_synthetic.py
    uv run python rebuttal/make_table1_synthetic.py --families clean_block --sizes small
"""

from __future__ import annotations

import argparse
import csv
import io
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import make_table1_aligned as mt  # noqa: E402  (reuses tw_adaptive / metric_triple)
import tw_components as tc  # noqa: E402

sys.path.insert(0, str(tc.SYNTHETIC_DIR))
from run_synthetic_evaluation import FAMILIES, SIZES, build_case  # noqa: E402

HERE = Path(__file__).resolve().parent
OUTPUT_DIR = HERE / "output"
# The benchmark as submitted, matching make_table1_aligned.SUBMITTED_REV.
SUBMITTED_REV = mt.SUBMITTED_REV
COMMITTED_CSV = (
    tc.REPO_ROOT
    / "output"
    / "main_experiment"
    / "synthetic_benchmark"
    / "synthetic_strong_baselines_paper.csv"
)

NUM_SEEDS = 20
ALPHAS = mt.ALPHAS
METHODS = mt.METHODS
# Shared with the real-data table so the two cannot drift apart.
METHOD_SHORT = mt.METHOD_SHORT
METRICS = ("2-SUM", "Band@10%", "MWB-AUC")
# The manuscript lists large first within each family.
PAPER_FAMILY_ORDER = (
    "clean_block",
    "paired_overlap",
    "shared_super",
    "shared_super_noisy",
    "cross_block_leakage",
)
PAPER_SIZE_ORDER = ("large", "medium", "small")


def regime_label(family) -> str:
    """`Family A: Clean one-to-one` -> `A: Clean one-to-one`, as in Table 1."""
    return family.name.removeprefix("Family ")


def _empty_samples() -> dict[str, dict[str, list[float]]]:
    """Per-method, per-metric lists of the 20 seed values."""
    return {method: {metric: [] for metric in METRICS} for method in METHODS}


def _mean_table(samples: dict[str, dict[str, list[float]]]) -> dict[str, dict[str, float]]:
    return {
        method: {metric: float(np.mean(values)) for metric, values in samples[method].items()}
        for method in METHODS
    }


def _sd_table(samples: dict[str, dict[str, list[float]]]) -> dict[str, dict[str, float]]:
    """Sample standard deviation (ddof=1) over the seeds.

    Sample rather than population: the 20 seeds are draws from each family's
    generator, not the whole population, so ddof=1 is the right estimator. With
    n = 1 there is nothing to estimate and the deviation is reported as 0.
    """
    return {
        method: {
            metric: float(np.std(values, ddof=1)) if len(values) > 1 else 0.0
            for metric, values in samples[method].items()
        }
        for method in METHODS
    }


def _grouped(per_seed: list[dict]) -> dict[tuple[str, str], list[dict]]:
    """Per-seed records keyed by (family name, size name)."""
    grouped: dict[tuple[str, str], list[dict]] = {}
    for record in per_seed:
        grouped.setdefault((record["family"], record["size"]), []).append(record)
    return grouped


def _paired_deltas(records: list[dict], side: str, metric: str) -> np.ndarray:
    """Per-seed TW - One-walk on one metric, for one regime.

    Paired on the seed, which is the point: both methods see the same instance, so
    the difference removes the between-instance variation that dominates the raw
    spread and leaves only the systematic gap.
    """
    return np.array(
        [
            record[f"{side}_{metric}|TW"] - record[f"{side}_{metric}|One-walk"]
            for record in records
        ]
    )


def _attribute(row: dict, prefix: str) -> dict[str, float]:
    return {
        f"{prefix}_{metric}|{method}": row[method][metric]
        for method in METHODS
        for metric in METRICS
    }


def run_case(observed: np.ndarray):
    """Both configurations for one instance. Returns (a_row, b_row, a_alpha, b_alpha).

    The five baselines keep one ordering across both tables; only One-walk and
    Two-walk move, so they are computed once and shared.
    """
    shared = mt.baseline_orderings(observed)

    # --- A: as submitted -----------------------------------------------------
    a_alpha, tw_submitted = mt.tw_adaptive(observed, tc.order_single, mt._select_released)
    a_orders = dict(shared)
    a_orders["One-walk"] = tc.order_single(observed, 1.0, "ow")
    a_orders["TW"] = tw_submitted
    a_row = {
        method: mt.metric_triple(observed, *a_orders[method], "released")
        for method in METHODS
    }

    # --- B: corrected --------------------------------------------------------
    b_alpha, tw_corrected = mt.tw_adaptive(
        observed, tc.order_componentwise, mt._select_released
    )
    b_orders = dict(shared)
    b_orders["One-walk"] = tc.order_componentwise(observed, 0.0, "ow")
    b_orders["TW"] = tw_corrected
    b_row = {
        method: mt.metric_triple(observed, *b_orders[method], "faithful")
        for method in METHODS
    }
    return a_row, b_row, a_alpha, b_alpha


def collect(families, sizes, seeds: int):
    """Per-regime metric means, standard deviations and alpha counts, plus records."""
    regimes: dict[tuple[str, str], dict] = {}
    per_seed: list[dict] = []

    for family in families:
        for size in sizes:
            label = regime_label(family)
            a_samples, b_samples = _empty_samples(), _empty_samples()
            a_alphas: Counter = Counter()
            b_alphas: Counter = Counter()

            started = time.perf_counter()
            for seed in range(seeds):
                observed = build_case(family, size, seed)["observed"]
                a_row, b_row, a_alpha, b_alpha = run_case(observed)
                for method in METHODS:
                    for metric in METRICS:
                        a_samples[method][metric].append(a_row[method][metric])
                        b_samples[method][metric].append(b_row[method][metric])
                a_alphas[a_alpha] += 1
                b_alphas[b_alpha] += 1
                per_seed.append(
                    {
                        "family": family.name,
                        "size": size.name,
                        "seed": seed,
                        **_attribute(a_row, "a"),
                        "a_alpha": a_alpha,
                        **_attribute(b_row, "b"),
                        "b_alpha": b_alpha,
                    }
                )

            regimes[(family.key, size.key)] = {
                "label": label,
                "family_name": family.name,
                "size_name": size.name,
                "shape": f"{size.n_rows}x{size.n_cols}",
                "n_seeds": seeds,
                "a": _mean_table(a_samples),
                "a_sd": _sd_table(a_samples),
                "b": _mean_table(b_samples),
                "b_sd": _sd_table(b_samples),
                "a_alpha": a_alphas,
                "b_alpha": b_alphas,
            }
            print(
                f"  {label} / {size.name}: {time.perf_counter() - started:.1f}s",
                flush=True,
            )
    return regimes, per_seed


def _ordered(regimes: dict) -> list[dict]:
    """Regimes in the manuscript's family-then-size order."""
    return [
        regimes[(family_key, size_key)]
        for family_key in PAPER_FAMILY_ORDER
        for size_key in PAPER_SIZE_ORDER
        if (family_key, size_key) in regimes
    ]


def _modal(counter: Counter):
    """Most common alpha, smallest on a tie. Returns (value, agreeing, total)."""
    total = sum(counter.values())
    value, count = min(counter.items(), key=lambda kv: (-kv[1], kv[0]))
    return value, count, total


def _table(title: str, note: str, regimes: dict, side: str) -> list[str]:
    header = (
        "| Dataset | Shape | "
        + " | ".join(
            f"{metric} {short}" for metric in METRICS for short in METHOD_SHORT
        )
        + " | alpha |"
    )
    separator = "|---|---|" + "---:|" * (len(METRICS) * len(METHODS)) + "---:|"
    lines = [f"## {title}", "", note, "", header, separator]

    for entry in _ordered(regimes):
        rows, sds = entry[side], entry[f"{side}_sd"]
        cells = []
        for metric in METRICS:
            places = 2 if metric == "2-SUM" else 3
            values = [rows[method][metric] for method in METHODS]
            best = min(values) if metric == "2-SUM" else max(values)
            best_text = f"{best:.{places}f}"
            for method in METHODS:
                # Bold on the mean only -- the deviation must not decide the winner.
                mean_text = f"{rows[method][metric]:.{places}f}"
                text = f"{mean_text}±{sds[method][metric]:.{places}f}"
                cells.append(f"**{text}**" if mean_text == best_text else text)
        alpha, agree, total = _modal(entry[f"{side}_alpha"])
        alpha_text = f"{alpha:g}" + ("" if agree == total else f" ({agree}/{total})")
        lines.append(
            f"| {entry['label']} | {entry['shape']} | "
            + " | ".join(cells)
            + f" | {alpha_text} |"
        )
    lines.append("")
    return lines


# The methods the two-walk claim actually rests on.
DISPERSION_METHODS = ("One-walk", "Gram-OW", "TW")
DISPERSION_SHORT = ("OW", "GO", "TW")
# Which table the paired summary is computed on. B is the corrected pipeline and
# the one the paper reports, so that is the one a reviewer's question is about.
PAIRED_SIDE = "b"
PAIRED_LABEL = "the corrected pipeline (table B)"


def _places(metric: str) -> int:
    return 2 if metric == "2-SUM" else 3


def _dispersion_section(regimes: dict, seeds: int) -> list[str]:
    """Section D: mean +- sd for the three methods the comparison turns on."""
    lines = [
        "## D. Seed dispersion, key methods",
        "",
        f"Mean $\\pm$ sample standard deviation (ddof=1) over the {seeds} seeds, for "
        "the three methods the two-walk claim rests on. The main tables above carry "
        "`mean±sd` in every cell; `synthetic_summary.csv` has the same numbers for "
        "all eight methods, including the deviations the tables above do not show.",
        "",
        "| Regime | Shape | "
        + " | ".join(
            f"{metric} {short}" for metric in METRICS for short in DISPERSION_SHORT
        )
        + " |",
        "|---|---|" + "---:|" * (len(METRICS) * len(DISPERSION_METHODS)),
    ]
    for entry in _ordered(regimes):
        cells = []
        for metric in METRICS:
            places = _places(metric)
            for method in DISPERSION_METHODS:
                mean = entry["b"][method][metric]
                sd = entry["b_sd"][method][metric]
                cells.append(f"{mean:.{places}f}±{sd:.{places}f}")
        lines.append(
            f"| {entry['label']} | {entry['shape']} | " + " | ".join(cells) + " |"
        )
    lines.append("")
    return lines


def _paired_section(regimes: dict, per_seed: list[dict], seeds: int) -> list[str]:
    """Section E: per-seed paired TW - One-walk, which is the comparative claim."""
    grouped = _grouped(per_seed)
    lines = [
        "## E. Paired difference against One-walk",
        "",
        f"TW minus One-walk on the same seed, under {PAIRED_LABEL}. Every method sees "
        "the identical instance, so pairing removes the between-instance variation "
        "that dominates the raw spread in section D and leaves only the systematic "
        "gap. `wins` counts seeds where TW is the better of the two "
        f"(lower for 2-SUM, higher for the other two), out of {seeds}. A negative "
        "2-SUM gap and a positive Band/MWB gap both favour Two-walk.",
        "",
        "| Regime | "
        + " | ".join(
            f"{metric} gap | {metric} sd | {metric} wins" for metric in METRICS
        )
        + " |",
        "|---|" + "---:|" * (len(METRICS) * 3),
    ]

    totals = {metric: [0, 0] for metric in METRICS}
    for entry in _ordered(regimes):
        records = grouped[(entry["family_name"], entry["size_name"])]
        cells = []
        for metric in METRICS:
            places = _places(metric)
            deltas = _paired_deltas(records, PAIRED_SIDE, metric)
            wins = int(
                (deltas < -1e-12).sum() if metric == "2-SUM" else (deltas > 1e-12).sum()
            )
            totals[metric][0] += wins
            totals[metric][1] += len(deltas)
            cells += [
                f"{deltas.mean():+.{places}f}",
                f"{deltas.std(ddof=1):.{places}f}",
                f"{wins}/{len(deltas)}",
            ]
        lines.append(f"| {entry['label']} | " + " | ".join(cells) + " |")

    lines.append("")
    lines.append(
        "Across all regimes: "
        + "; ".join(
            f"{metric} {totals[metric][0]}/{totals[metric][1]} seeds" for metric in METRICS
        )
        + "."
    )
    lines.append("")
    return lines


def _write_rows(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _summary_rows(regimes: dict) -> list[dict]:
    """Long-format mean and sd, one row per regime x config x method x metric."""
    rows = []
    for entry in _ordered(regimes):
        for side, config in (("a", "as-submitted"), ("b", "corrected")):
            for method in METHODS:
                for metric in METRICS:
                    rows.append(
                        {
                            "family": entry["family_name"],
                            "size": entry["size_name"],
                            "config": config,
                            "method": method,
                            "metric": metric,
                            "mean": entry[side][method][metric],
                            "sd": entry[f"{side}_sd"][method][metric],
                            "n_seeds": entry["n_seeds"],
                        }
                    )
    return rows


def _paired_rows(regimes: dict, per_seed: list[dict]) -> list[dict]:
    """Per-seed paired TW - One-walk, one row per regime x metric."""
    grouped = _grouped(per_seed)
    rows = []
    for entry in _ordered(regimes):
        records = grouped[(entry["family_name"], entry["size_name"])]
        for metric in METRICS:
            deltas = _paired_deltas(records, PAIRED_SIDE, metric)
            wins = int(
                (deltas < -1e-12).sum() if metric == "2-SUM" else (deltas > 1e-12).sum()
            )
            rows.append(
                {
                    "family": entry["family_name"],
                    "size": entry["size_name"],
                    "config": PAIRED_SIDE,
                    "metric": metric,
                    "tw_minus_ow_mean": float(deltas.mean()),
                    "tw_minus_ow_sd": float(deltas.std(ddof=1)),
                    "tw_wins": wins,
                    "n_seeds": len(deltas),
                }
            )
    return rows


def _read_reference(text: str) -> dict[tuple[str, str], dict[str, str]]:
    return {
        (row["family"], row["size"]): row for row in csv.DictReader(io.StringIO(text))
    }


def _submitted_reference():
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
    return _read_reference(proc.stdout)


def _compare(regimes, side: str, reference, columns) -> list[str]:
    """Cells whose mean differs from the reference CSV beyond its rounding."""
    drift = []
    for entry in _ordered(regimes):
        match = reference.get((entry["family_name"], entry["size_name"]))
        if match is None:
            drift.append(f"- {entry['label']} / {entry['size_name']}: no reference row")
            continue
        for method, metric, column in columns:
            mine = entry[side][method][metric]
            theirs = float(match[column])
            # The CSV stores two decimals, so half a unit in the last place is noise.
            if abs(mine - theirs) > 5e-3:
                drift.append(
                    f"- {entry['label']} / {entry['size_name']} | {method} | "
                    f"{metric}: ours {mine:.4f} vs reference {theirs:.4f}"
                )
    return drift


def _reproduction_check(regimes) -> list[str]:
    submitted_columns = tuple(
        (method, metric, f"{csv_metric}|{method}")
        for method in ("One-walk", "TW")
        for metric, csv_metric in (
            ("2-SUM", "two_sum_mean"),
            ("Band@10%", "band_mass_10_mean"),
            ("MWB-AUC", "mwb_auc_mean"),
        )
    )
    # Table B's MWB-AUC column is the block-cut quantity, so the CSV's band
    # integral is not comparable there.
    corrected_columns = tuple(
        (method, metric, f"{csv_metric}|{method}")
        for method in ("One-walk", "TW")
        for metric, csv_metric in (
            ("2-SUM", "two_sum_mean"),
            ("Band@10%", "band_mass_10_mean"),
        )
    )

    lines = [
        f"### A. As submitted vs the `{SUBMITTED_REV}` CSV",
        "",
        "Means over the seeds, compared against the submitted table's two-decimal",
        "values.",
        "",
    ]
    submitted = _submitted_reference()
    if submitted is None:
        lines += [f"Revision `{SUBMITTED_REV}` not available; skipped.", ""]
    else:
        drift = _compare(regimes, "a", submitted, submitted_columns)
        lines += (
            ["All checked cells reproduce the submitted CSV to within 5e-3.", ""]
            if not drift
            else ["Cells that differ:", ""] + drift + [""]
        )

    lines += [
        "### B. Corrected vs the current CSV",
        "",
        "MWB-AUC is excluded: Table B reports the Section 4.5 block-cut quantity",
        "while the CSV carries the selector's band integral.",
        "",
    ]
    if not COMMITTED_CSV.exists():
        lines += ["Working-tree CSV not found; skipped.", ""]
    else:
        current = _read_reference(COMMITTED_CSV.read_text(encoding="utf-8"))
        drift = _compare(regimes, "b", current, corrected_columns)
        lines += (
            ["All checked cells reproduce the current CSV to within 5e-3.", ""]
            if not drift
            else ["Cells that differ:", ""] + drift + [""]
        )
    return lines


def _split_sections(text: str) -> tuple[str, dict[str, list[str]]]:
    """Split a generated file into its preamble and its `## `-headed sections."""
    preamble: list[str] = []
    sections: dict[str, list[str]] = {}
    current: str | None = None
    for line in text.splitlines():
        if line.startswith("## "):
            current = line[3:].strip()
            sections[current] = [line]
        elif current is None:
            preamble.append(line)
        else:
            sections[current].append(line)
    return "\n".join(preamble).strip(), sections


def combine_with_real(output_dir: Path) -> Path | None:
    """Splice the synthetic and real-data tables into one Table 1.

    The manuscript's Table 1 is a single table with a synthetic component on top and
    a real-world component below, so the two generated files are interleaved here
    rather than left for the reader to join by hand. Section order follows the
    real-data file, which carries the extra alpha section.
    """
    real_path = output_dir / "table1_aligned.md"
    if not real_path.exists():
        return None

    _, real = _split_sections(real_path.read_text(encoding="utf-8"))
    _, synth = _split_sections(
        (output_dir / "table1_synthetic.md").read_text(encoding="utf-8")
    )

    def body(sections: dict[str, list[str]], title: str) -> list[str]:
        return sections.get(title, [])[1:]  # drop the section's own `## ` line

    lines = [
        "# Table 1 -- synthetic and real-world components, aligned with the manuscript",
        "",
        "The two components of Table 1 in one file. Column groups are the three",
        "metrics; within each group the methods are O = Original, M = Marginal,",
        "HC = HC+OLO, OW = One-walk, GO = Gram-OW, CA = CA-SVD, MD = Median,",
        "TW = Two-walk. Bold marks the best value per row and metric. Synthetic rows",
        "are `mean±sd` over 20 seeds (sample sd, ddof=1); real-world rows are single",
        "full-size matrices and carry no deviation.",
        "",
        "Sources: `table1_synthetic.md` and `table1_aligned.md`.",
        "",
    ]
    # The real-data file sets the section skeleton (it carries the extra alpha
    # section), but the synthetic file has sections of its own -- the dispersion and
    # paired summaries. Union the titles rather than iterating one file's, or those
    # get dropped from the combined output.
    titles = list(real) + [title for title in synth if title not in real]
    for title in titles or list(synth):
        lines += [f"## {title}", ""]
        for heading, sections in (
            ("Synthetic component", synth),
            ("Real-world component", real),
        ):
            if title not in sections:
                continue
            lines += [f"### {heading}", ""]
            lines += body(sections, title)
            lines.append("")

    path = output_dir / "table1_full.md"
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")
    return path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--families", nargs="*", default=None)
    parser.add_argument("--sizes", nargs="*", default=None)
    parser.add_argument("--seeds", type=int, default=NUM_SEEDS)
    parser.add_argument(
        "--combine-only",
        action="store_true",
        help="re-splice table1_full.md from the two existing outputs, no recomputation",
    )
    args = parser.parse_args()

    if args.combine_only:
        combined = combine_with_real(OUTPUT_DIR)
        if combined is None:
            raise SystemExit(
                "need both table1_synthetic.md and table1_aligned.md to combine"
            )
        print(f"Wrote {combined}")
        return

    families = [f for f in FAMILIES if args.families is None or f.key in args.families]
    sizes = [s for s in SIZES if args.sizes is None or s.key in args.sizes]
    if not families or not sizes:
        raise SystemExit("no families or sizes selected")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    print(
        f"Building {len(families)} families x {len(sizes)} sizes x {args.seeds} seeds ..."
    )
    started = time.perf_counter()
    regimes, per_seed = collect(families, sizes, args.seeds)
    print(f"  done in {time.perf_counter() - started:.1f}s")

    a_note = (
        "The submitted code: one Fiedler vector over the whole active support, "
        "finished with `orient_orders_for_diagonal` and an unstable sort, under the "
        f"released metric definitions. Should reproduce the `{SUBMITTED_REV}` table. "
        "Each row is a mean over the seeds."
    )
    b_note = (
        "Section 3.3 component-wise ordering, alpha chosen by the released MWB-AUC, "
        "and the paper-faithful metric definitions (Eq. 16 Band@10%, Section 4.5 "
        "block-cut MWB-AUC). The synthetic active supports are all connected, so the "
        "component rule itself changes nothing here; the movement comes from the "
        "orientation and stable-sort fixes the paper-spec commit bundled with it. "
        "Each row is a mean over the seeds."
    )

    lines = [
        "# Table 1, synthetic component -- aligned with the manuscript",
        "",
        "Same column structure as the real-data tables in `table1_aligned.md`: three",
        "metric groups, one column per method, plus the selected alpha. O = Original,",
        "M = Marginal, HC = HC+OLO, OW = One-walk, GO = Gram-OW, CA = CA-SVD,",
        "MD = Median, TW = Two-walk. Bold marks the best displayed value per row and",
        "metric. Rows follow the manuscript's family order, large to small within each",
        "family.",
        "",
        "GO is the reviewer XFNF W5 baseline, described in `baselines_v2.py`. As in",
        "the real-data table it is the sixth ordering shared by both configurations,",
        "and it is allowed the released metric-fitted orientation that OW and TW no",
        "longer receive -- so it is flattered here, if anything.",
        "",
        "Alpha is chosen per instance, so a regime has one selection per seed. The",
        "alpha column gives the modal value; a parenthesised `k/N` means the seeds did",
        "not agree. `synthetic_per_seed.csv` carries every seed's record.",
        "",
        f"Every cell is `mean±sd` over the {args.seeds} seeds, with the sample standard",
        "deviation (ddof=1). The deviation expresses seed-to-seed variation in the",
        "generator draw *and*, for One-walk and Two-walk, in the alpha the selector",
        "picks on that draw; it is not a confidence interval on the mean. Rows are",
        "bolded on the mean alone. Two caveats worth reading before quoting a spread:",
        "a method with a large sd is not necessarily worse on the mean, and the paired",
        "comparison in section E is the tighter statement because it removes the",
        "between-instance variation entirely.",
        "",
    ]
    lines += _table("A. As submitted", a_note, regimes, "a")
    lines += _table("B. Corrected", b_note, regimes, "b")
    lines += ["## C. Reproduction check", ""]
    lines += _reproduction_check(regimes)
    lines += _dispersion_section(regimes, args.seeds)
    lines += _paired_section(regimes, per_seed, args.seeds)

    (OUTPUT_DIR / "table1_synthetic.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )
    _write_rows(OUTPUT_DIR / "synthetic_per_seed.csv", per_seed)
    summary = _summary_rows(regimes)
    _write_rows(OUTPUT_DIR / "synthetic_summary.csv", summary)
    paired = _paired_rows(regimes, per_seed)
    _write_rows(OUTPUT_DIR / "synthetic_paired.csv", paired)
    print(f"\nWrote {OUTPUT_DIR / 'table1_synthetic.md'}")
    print(f"Wrote {OUTPUT_DIR / 'synthetic_per_seed.csv'} ({len(per_seed)} records)")
    print(f"Wrote {OUTPUT_DIR / 'synthetic_summary.csv'} ({len(summary)} rows)")
    print(f"Wrote {OUTPUT_DIR / 'synthetic_paired.csv'} ({len(paired)} rows)")
    combined = combine_with_real(OUTPUT_DIR)
    if combined is not None:
        print(f"Wrote {combined}")
    else:
        print("Skipped the combined table: run make_table1_aligned.py first.")


if __name__ == "__main__":
    main()
