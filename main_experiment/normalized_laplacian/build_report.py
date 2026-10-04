"""Build the canonical English normalized-Laplacian ablation reports.

The runner writes measurements and paired uncertainty estimates. This script
only formats those measurements; it does not rerun, select, or alter orderings.
The LaTeX output is a single-section fragment, requiring amsmath, booktabs,
graphicx, longtable, and pdflscape in the enclosing document.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Sequence

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = REPO_ROOT / "output" / "main_experiment" / "normalized_laplacian"
METRICS = ("two_sum", "band_mass_10", "mwb_auc")
COMPARISONS = ("OW", "TW fixed alpha=1", "TW fixed alpha=12", "TW adaptive")
VARIANTS = (
    "c_ow", "n_ow", "c_tw_fixed_1", "n_tw_fixed_1",
    "c_tw_fixed_12", "n_tw_fixed_12", "c_tw_adaptive", "n_tw_adaptive",
)
VARIANT_LABELS = {
    "c_ow": "OW combinatorial",
    "n_ow": "OW normalized",
    "c_tw_fixed_1": "TW combinatorial, fixed alpha=1",
    "n_tw_fixed_1": "TW normalized, fixed alpha=1",
    "c_tw_fixed_12": "TW combinatorial, fixed alpha=12",
    "n_tw_fixed_12": "TW normalized, fixed alpha=12",
    "c_tw_adaptive": "TW combinatorial, adaptive",
    "n_tw_adaptive": "TW normalized, adaptive",
}
FAMILY_ORDER = (
    "Family A: Clean one-to-one",
    "Family B: Paired subgroup overlap",
    "Family C: Shared super-prototype",
    "Family D: Shared prototype with noise",
    "Family E: Cross-block leakage",
)
FAMILY_SHORT = {
    FAMILY_ORDER[0]: "A: Clean one-to-one",
    FAMILY_ORDER[1]: "B: Paired overlap",
    FAMILY_ORDER[2]: "C: Shared super-prototype",
    FAMILY_ORDER[3]: "D: Shared + noise",
    FAMILY_ORDER[4]: "E: Cross-block leakage",
}
SIZE_ORDER = ("Large", "Medium", "Small")
DATASET_ORDER = (
    "SIC -> NAICS", "CIP -> SOC", "ACS OCCP x INDP",
    "LODES Home x Work County", "20 Newsgroups Term x Document",
    "MBTA GTFS Route x Station", "OpenAlex Author x Topic",
)
DATASET_SHORT = {
    DATASET_ORDER[0]: "SIC -> NAICS",
    DATASET_ORDER[1]: "CIP -> SOC",
    DATASET_ORDER[2]: "ACS OCCP x INDP",
    DATASET_ORDER[3]: "LODES Home x Work",
    DATASET_ORDER[4]: "20NG Term x Document",
    DATASET_ORDER[5]: "MBTA Route x Station",
    DATASET_ORDER[6]: "OpenAlex Author x Topic",
}


def tex_escape(value: object) -> str:
    replacements = {
        "\\": r"\textbackslash{}", "&": r"\&", "%": r"\%",
        "$": r"\$", "#": r"\#", "_": r"\_", "{": r"\{",
        "}": r"\}", "~": r"\textasciitilde{}", "^": r"\textasciicircum{}",
    }
    return "".join(replacements.get(char, char) for char in str(value))


def fmt(value: float, *, signed: bool = False, scale: float = 1.0) -> str:
    if pd.isna(value):
        return "--"
    return format(float(value) * scale, "+.4f" if signed else ".4f")


def alpha_fmt(value: float) -> str:
    return "--" if pd.isna(value) else f"{float(value):g}"


def md_table(headers: Sequence[str], rows: Sequence[Sequence[object]]) -> str:
    lines = ["| " + " | ".join(headers) + " |",
             "| " + " | ".join("---" if i == 0 else "---:" for i in range(len(headers))) + " |"]
    lines.extend("| " + " | ".join(str(value) for value in row) + " |" for row in rows)
    return "\n".join(lines)


def tex_table(headers: Sequence[str], rows: Sequence[Sequence[object]],
              caption: str, label: str, *, landscape: bool = False) -> str:
    n_columns = len(headers)
    header = " & ".join(tex_escape(value) for value in headers) + r" \\"
    lines = []
    if landscape:
        lines.append(r"\begin{landscape}")
    lines.extend([
        r"\begingroup", r"\scriptsize", r"\setlength{\tabcolsep}{3pt}",
        "\\begin{longtable}{l" + "r" * (n_columns - 1) + "}",
        "\\caption{" + tex_escape(caption) + "}\\label{" + label + r"}\\",
        r"\toprule", header, r"\midrule", r"\endfirsthead",
        "\\multicolumn{" + str(n_columns) + r"}{c}{\tablename\ \thetable\ (continued)}\\",
        r"\toprule", header, r"\midrule", r"\endhead", r"\midrule",
        "\\multicolumn{" + str(n_columns) + r"}{r}{Continued on next page}\\",
        r"\endfoot", r"\bottomrule", r"\endlastfoot",
    ])
    lines.extend(" & ".join(tex_escape(value) for value in row) + r" \\" for row in rows)
    lines.extend([r"\end{longtable}", r"\endgroup"])
    if landscape:
        lines.append(r"\end{landscape}")
    return "\n".join(lines)


def group_label(row: pd.Series) -> str:
    if row["suite"] == "synthetic":
        return FAMILY_SHORT.get(str(row["family"]), str(row["family"])) + " / " + str(row["size"])
    return DATASET_SHORT.get(str(row["dataset"]), str(row["dataset"]))


def group_key(row: pd.Series) -> tuple[object, ...]:
    if row["suite"] == "synthetic":
        return (0, FAMILY_ORDER.index(str(row["family"])), SIZE_ORDER.index(str(row["size"])))
    return (1, DATASET_ORDER.index(str(row["dataset"])))


def paired_blocks(data: pd.DataFrame, *, comparison: str, scope: str,
                  suite: str | None = None) -> list[pd.DataFrame]:
    subset = data.loc[(data["comparison"] == comparison) & (data["group_scope"] == scope)]
    if suite is not None:
        subset = subset.loc[subset["suite"] == suite]
    columns = ["suite", "family", "size", "dataset"]
    blocks = [group.set_index("metric") for _, group in subset.groupby(columns, dropna=False, sort=False)]
    if scope == "overall":
        return sorted(blocks, key=lambda block: 0 if block.iloc[0]["suite"] == "synthetic" else 1)
    return sorted(blocks, key=lambda block: group_key(block.iloc[0]))


SCORE_HEADERS = (
    "Group", "n", "C R2S x100", "N R2S x100", "Delta R2S x100",
    "C Band", "N Band", "Delta Band", "C MWB", "N MWB", "Delta MWB",
)
CI_HEADERS = (
    "Group", "Delta R2S x100 [95% CI]", "R2S wins",
    "Delta Band [95% CI]", "Band wins", "Delta MWB [95% CI]", "MWB wins",
)


def score_row(block: pd.DataFrame, label: str) -> list[object]:
    values: list[object] = [label, int(block.iloc[0]["n_cases"])]
    for metric in METRICS:
        scale = 100.0 if metric == "two_sum" else 1.0
        row = block.loc[metric]
        values.extend((fmt(row["combinatorial_mean"], scale=scale),
                       fmt(row["normalized_mean"], scale=scale),
                       fmt(row["delta_mean"], signed=True, scale=scale)))
    return values


def ci_row(block: pd.DataFrame, label: str) -> list[object]:
    values: list[object] = [label]
    for metric in METRICS:
        row = block.loc[metric]
        scale = 100.0 if metric == "two_sum" else 1.0
        interval = (fmt(row["delta_mean"], signed=True, scale=scale) + " ["
                    + fmt(row["delta_ci_low"], signed=True, scale=scale) + ", "
                    + fmt(row["delta_ci_high"], signed=True, scale=scale) + "]")
        values.extend((interval, f"{100.0 * float(row['normalized_win_rate']):.1f}%"))
    return values


def mean_row(summary: pd.DataFrame, suite: str, variant: str) -> pd.Series:
    subset = summary.loc[(summary["suite"] == suite) & (summary["variant"] == variant)]
    if len(subset) != 1:
        raise ValueError(f"Expected exactly one summary for {suite}/{variant}")
    return subset.iloc[0]


def relative_r2s(summary: pd.DataFrame, suite: str, c_variant: str, n_variant: str) -> float:
    c_value = float(mean_row(summary, suite, c_variant)["two_sum_mean"])
    n_value = float(mean_row(summary, suite, n_variant)["two_sum_mean"])
    return 100.0 * (n_value - c_value) / c_value


def interpret(summary: pd.DataFrame, deltas: pd.DataFrame) -> list[str]:
    paragraphs = [
        "The OW control changes mean R2S by "
        + f"{relative_r2s(summary, 'synthetic', 'c_ow', 'n_ow'):+.1f}% on synthetic data and "
        + f"{relative_r2s(summary, 'real', 'c_ow', 'n_ow'):+.1f}% on real data "
        + "under degree normalization. This control measures the corresponding effect on one-step spectral ordering."
    ]
    for suffix, label in (("fixed_1", "Fixed alpha=1"), ("fixed_12", "Frozen global alpha=12"),
                          ("adaptive", "Adaptive six-alpha search")):
        c_variant, n_variant = "c_tw_" + suffix, "n_tw_" + suffix
        parts = []
        for suite in ("synthetic", "real"):
            c_row = mean_row(summary, suite, c_variant)
            n_row = mean_row(summary, suite, n_variant)
            parts.append(
                suite + " mean R2S " + f"{relative_r2s(summary, suite, c_variant, n_variant):+.1f}%"
                + ", Band " + fmt(n_row["band_mass_10_mean"] - c_row["band_mass_10_mean"], signed=True)
                + ", MWB-AUC " + fmt(n_row["mwb_auc_mean"] - c_row["mwb_auc_mean"], signed=True)
            )
        paragraphs.append(label + ": normalized minus combinatorial changes are " + "; ".join(parts) + ".")
    for comparison in COMPARISONS[1:]:
        counts = []
        for suite in ("synthetic", "real"):
            blocks = paired_blocks(deltas, comparison=comparison, scope="group", suite=suite)
            wins = sum(all(float(block.loc[metric, "delta_mean"]) < -1e-12
                           if metric == "two_sum" else float(block.loc[metric, "delta_mean"]) > 1e-12
                           for metric in METRICS) for block in blocks)
            counts.append(f"{wins}/{len(blocks)} {suite} groups")
        paragraphs.append(comparison + ": normalization improves all three group-mean metrics in "
                          + " and ".join(counts) + ". A group-mean improvement is not an instance-level win rate.")
    paragraphs.append(
        "Fixed-alpha comparisons isolate the choice of Laplacian and coordinate. Adaptive comparisons also allow "
        "the selected alpha to change, while keeping the candidate grid and selector identical. These experiments "
        "evaluate degree normalization of the final Two-Walk graph; they do not test normalized Gram projections. "
        "MWB-AUC is the adaptive selection objective, so its adaptive score is selection-aligned evidence. "
        "The seven real matrices are descriptive cases, not independent repeated trials or evidence of universal superiority."
    )
    return paragraphs


def reviewer_answer(summary: pd.DataFrame) -> str:
    rates = {
        (suite, suffix): relative_r2s(summary, suite, "c_tw_" + suffix, "n_tw_" + suffix)
        for suite in ("synthetic", "real") for suffix in ("fixed_1", "fixed_12", "adaptive")
    }
    text = (
        "We additionally evaluated a degree-normalized Two-Walk variant on the same 300 synthetic instances "
        "and 7 real matrices. We retain A_alpha=M^2+alpha M, remove self-loops, and use the symmetric normalized "
        "Laplacian D^(-1/2)(D-W)D^(-1/2), sorting the generalized Fiedler coordinate f=D^(-1/2)u component-wise. "
        "We compare fixed alpha=1, the previously frozen global alpha=12, and the same six-alpha adaptive search; "
        "all metrics use the original matrix weights. Relative to combinatorial TW, normalized TW changes mean "
        f"R2S by {rates['synthetic', 'fixed_1']:+.1f}%/{rates['real', 'fixed_1']:+.1f}% at alpha=1, "
        f"{rates['synthetic', 'fixed_12']:+.1f}%/{rates['real', 'fixed_12']:+.1f}% at alpha=12, and "
        f"{rates['synthetic', 'adaptive']:+.1f}%/{rates['real', 'adaptive']:+.1f}% adaptively "
        "(synthetic/real; lower is better). We report all three metrics and paired synthetic confidence intervals. "
        "This directly tests Laplacian normalization; normalized projections are a separate variant."
    )
    if len(text) >= 1200:
        raise ValueError("Reviewer answer exceeded the report's character budget")
    return text


def validate(raw: pd.DataFrame, protocols: pd.DataFrame, summary: pd.DataFrame,
             alpha_summary: pd.DataFrame, deltas: pd.DataFrame) -> None:
    if set(raw["laplacian"]) != {"combinatorial", "normalized"} or set(raw["method"]) != {"OW", "TW"}:
        raise ValueError("Both Laplacians and OW/TW controls are required")
    if set(protocols["variant"]) != set(VARIANTS):
        raise ValueError("All eight protocol variants are required")
    counts = raw.groupby("suite")["case_id"].nunique().to_dict()
    if counts != {"real": 7, "synthetic": 300} or len(raw) != 4298 or len(protocols) != 2456:
        raise ValueError(f"Incomplete benchmark: case counts {counts}, {len(raw)} raw, {len(protocols)} protocol rows")
    if raw.duplicated(["case_id", "laplacian", "method", "alpha"]).any():
        raise ValueError("Duplicate raw case/operator/method/alpha rows")
    if protocols.duplicated(["case_id", "variant"]).any():
        raise ValueError("Duplicate case/protocol rows")
    if len(summary) != 16 or len(alpha_summary) != 24:
        raise ValueError("Expected 16 method summaries and 24 alpha summaries")
    for comparison in COMPARISONS:
        if len(paired_blocks(deltas, comparison=comparison, scope="group")) != 22:
            raise ValueError(f"Expected 22 paired groups for {comparison}")
    if not all(math.isfinite(float(value)) for metric in METRICS for value in raw[metric]):
        raise ValueError("Non-finite benchmark metric")


def build(output: Path) -> None:
    raw = pd.read_csv(output / "normalized_laplacian_per_instance.csv")
    protocols = pd.read_csv(output / "protocol_comparison.csv")
    summary = pd.read_csv(output / "method_summary.csv")
    alpha_summary = pd.read_csv(output / "alpha_summary.csv")
    deltas = pd.read_csv(output / "paired_deltas.csv")
    metadata = json.loads((output / "experiment_metadata.json").read_text(encoding="utf-8"))
    validate(raw, protocols, summary, alpha_summary, deltas)

    md: list[str] = ["# Normalized-Laplacian ablation", ""]
    tex: list[str] = [
        "% Requires amsmath, booktabs, graphicx, longtable, pdflscape.",
        r"\section{Normalized-Laplacian Ablation}", r"\label{sec:normalized-laplacian-ablation}",
    ]

    def text(md_text: str, tex_text: str | None = None) -> None:
        md.extend([md_text, ""])
        tex.extend([tex_escape(md_text) if tex_text is None else tex_text, ""])

    def heading(title: str) -> None:
        md.extend(["## " + title, ""])
        tex.extend(["\\subsection{" + tex_escape(title) + "}", ""])

    def table(title: str, headers: Sequence[str], rows: Sequence[Sequence[object]], label: str,
              *, landscape: bool = False) -> None:
        md.extend(["### " + title, "", md_table(headers, rows), ""])
        tex.extend([tex_table(headers, rows, title, label, landscape=landscape), ""])

    text(
        "This controlled ablation addresses reviewer JfQK Q2: whether degree-normalized Laplacians were considered. "
        "It compares the existing combinatorial method with a normalized operator on identical inputs. "
        "The canonical Table 1 remains the research baseline; this report is a separate normalization ablation.")
    heading("Definition and correspondence to the reviewer question")
    text(
        "Both variants first apply the existing global scale normalization, then use the same Two-Walk adjacency:",
        "Both variants first apply the existing global scale normalization, then use the same Two-Walk adjacency:")
    definition = r"""\begin{aligned}
\bar B&=\frac{B}{\max_{i,j}B_{ij}},\qquad
M=\begin{bmatrix}0&\bar B\\\bar B^\top&0\end{bmatrix},\\
A_\alpha&=M^2+\alpha M=
\begin{bmatrix}\bar B\bar B^\top&\alpha\bar B\\
\alpha\bar B^\top&\bar B^\top\bar B\end{bmatrix}.
\end{aligned}"""
    text("$$\n" + definition + "\n$$", "\\[\n" + definition + "\n\\]")
    text(
        "To make the degree convention explicit, self-loops are excluded from the normalization. This does not change "
        "the combinatorial Laplacian, because each self-loop cancels between its diagonal degree and adjacency term:")
    normalized = r"""\begin{aligned}
W_\alpha&=A_\alpha-\operatorname{diag}(\operatorname{diag}(A_\alpha)),\qquad
D_\alpha=\operatorname{diag}(W_\alpha\mathbf1),\\
L_\alpha&=D_\alpha-W_\alpha,\qquad
L_{\mathrm{sym},\alpha}=D_\alpha^{-1/2}L_\alpha D_\alpha^{-1/2}.
\end{aligned}"""
    text("$$\n" + normalized + "\n$$", "\\[\n" + normalized + "\n\\]")
    coordinate = r"""L_{\mathrm{sym},\alpha}u=\lambda u,\qquad
f=D_\alpha^{-1/2}u,\qquad L_\alpha f=\lambda D_\alpha f."""
    text(
        "The normalized variant sorts the generalized Fiedler coordinate f, obtained from the second eigenvector u "
        "of the symmetric normalized operator. It therefore uses degree-weighted centering and scale constraints "
        "rather than the Euclidean constraints of combinatorial spectral ordering:")
    text("$$\n" + coordinate + "\n$$", "\\[\n" + coordinate + "\n\\]")
    text(
        "For OW, the same comparison uses W=M. Active support components, component placement, sign and orientation "
        "conventions, and placement of zero rows and columns follow the existing benchmark policy. Degrees and "
        "normalization are computed within each active component. Positive alpha preserves the support components "
        "of the bipartite graph. This directly tests the reviewer's normalized-Laplacian alternative. It does not "
        "change the shared-neighbor projections to degree- or cosine-normalized Gram matrices.")

    heading("Protocol and evaluation")
    text(
        "The data are the same 300 synthetic matrices (five families, three sizes, 20 seeds per regime) and seven "
        "real matrices. Each matrix is evaluated with both Laplacians for OW and for TW at alpha in {1, 2, 4, 6, 8, 12}. "
        "This produces 4,298 raw ordering/metric records and 2,456 protocol records. The three TW protocols are "
        "fixed alpha=1, the previously selected global alpha=12, and adaptive per-matrix selection from this same "
        "six-value grid. Global alpha=12 is frozen from the previous synthetic-only study and is not retuned here. "
        "Both adaptive methods maximize the same retained MWB-AUC selector; ties within 1e-12 favor smaller alpha.")
    text(
        "Every score uses the permuted original matrix B, not the normalized graph weights. For normalized row/column "
        "positions x_i and y_j in [0,1], the two coordinate metrics are:")
    metrics = r"""\begin{aligned}
\mathrm{R2S}(B')&=\frac{\sum_{i,j}B'_{ij}(x_i-y_j)^2}{\sum_{i,j}B'_{ij}},\\
\mathrm{Band}_{0.10}(B')&=\frac{\sum_{i,j}B'_{ij}\mathbf1\{|x_i-y_j|\leq0.10\}}{\sum_{i,j}B'_{ij}}.
\end{aligned}"""
    text("$$\n" + metrics + "\n$$", "\\[\n" + metrics + "\n\\]")
    text(
        "R2S is minimized; Band@10% and MWB-AUC are maximized. MWB-AUC retains the existing integer-band score "
        "Q(w), integrating across 25 uniformly spaced widths from 0.02 to 0.50 by the trapezoidal rule:")
    mwb = r"""\mathrm{MWB\!\!\!-AUC}(B')\approx
\frac{1}{0.50-0.02}\sum_{k=1}^{24}\frac{Q(w_k)+Q(w_{k+1})}{2}(w_{k+1}-w_k),
\qquad w_k=0.02+0.02(k-1)."""
    text("$$\n" + mwb + "\n$$", "\\[\n" + mwb + "\n\\]")
    text(
        "For a zero-indexed n-by-m matrix, Q(w) retains entries j from max(0, floor(i*m/n-b)) inclusive to "
        "min(m, ceil(i*m/n+b+1)) exclusive, where b=max(1, floor(w*min(n,m))), and divides retained mass by total mass. "
        "Thus this retained selection metric is an integer-index band integral, distinct from Band@10%.")
    text(
        "All deltas are normalized minus combinatorial, paired by the same input matrix. Negative R2S deltas and "
        "positive Band/MWB deltas favor normalization. Tables display R2S multiplied by 100; Band and MWB remain "
        "in [0,1]. Synthetic group means aggregate 20 seeds; their 95% paired bootstrap intervals use 10,000 "
        "resamples with seed 20261004. Overall synthetic intervals use a seed-block bootstrap over the 20 "
        "seed-level macro deltas, each combining the 15 regimes with equal weight. The same seed is resampled "
        "jointly across regimes to preserve the common-random-number dependence induced by the data generator. "
        "These intervals describe variability within the synthetic "
        "benchmark; no real-data confidence interval is reported from a single matrix per dataset. Win rates "
        "are paired per-instance wins with tolerance 1e-12, excluding ties from the win numerator.")

    heading("Overall results and interpretation")
    headers = ("Suite / method", "n", "Mean alpha", "R2S x100", "Band@10%", "MWB-AUC")
    rows = []
    for suite in ("synthetic", "real"):
        for variant in VARIANTS:
            row = mean_row(summary, suite, variant)
            rows.append([suite.title() + " / " + VARIANT_LABELS[variant], int(row["n_cases"]),
                         alpha_fmt(row["selected_alpha_mean"]), fmt(row["two_sum_mean"], scale=100.0),
                         fmt(row["band_mass_10_mean"]), fmt(row["mwb_auc_mean"])])
    table("All OW and TW protocol means", headers, rows, "tab:normalized-method-means")
    rows = []
    ci_rows = []
    for comparison in COMPARISONS:
        for block in paired_blocks(deltas, comparison=comparison, scope="overall"):
            label = str(block.iloc[0]["suite"]).title() + " / " + comparison
            rows.append(score_row(block, label))
            if block.iloc[0]["suite"] == "synthetic":
                ci_rows.append(ci_row(block, comparison))
    table("Overall paired normalized-minus-combinatorial changes", SCORE_HEADERS, rows,
          "tab:normalized-overall-pairs", landscape=True)
    table("Overall synthetic paired uncertainty and instance win rates", CI_HEADERS, ci_rows,
          "tab:normalized-overall-ci", landscape=True)
    for paragraph in interpret(summary, deltas):
        text(paragraph)
    md.extend(["![Normalized-Laplacian comparison](figures/normalized_laplacian_comparison.png)", ""])
    tex.extend([
        r"\begin{figure}[htbp]", r"\centering",
        r"\includegraphics[width=\linewidth]{figures/normalized_laplacian_comparison.pdf}",
        r"\caption{Combinatorial and normalized Laplacian scores under matched alpha protocols. Lower R2S is better; higher Band@10\% and MWB-AUC are better.}",
        r"\label{fig:normalized-laplacian-comparison}", r"\end{figure}", "",
    ])

    heading("Complete six-alpha aggregate comparison")
    rows = []
    for suite in ("synthetic", "real"):
        alphas = sorted(alpha_summary.loc[alpha_summary["suite"] == suite, "alpha"].unique())
        for alpha in alphas:
            subset = alpha_summary.loc[(alpha_summary["suite"] == suite) & (alpha_summary["alpha"] == alpha)]
            c_row = subset.loc[subset["laplacian"] == "combinatorial"].iloc[0]
            n_row = subset.loc[subset["laplacian"] == "normalized"].iloc[0]
            values: list[object] = [suite.title() + " / alpha=" + alpha_fmt(alpha), int(c_row["n_cases"])]
            for metric in METRICS:
                scale = 100.0 if metric == "two_sum" else 1.0
                c_value, n_value = c_row[metric + "_mean"], n_row[metric + "_mean"]
                values.extend([fmt(c_value, scale=scale), fmt(n_value, scale=scale),
                               fmt(n_value - c_value, signed=True, scale=scale)])
            rows.append(values)
    table("Every alpha on both suites", SCORE_HEADERS, rows, "tab:normalized-every-alpha", landscape=True)

    heading("Complete results by synthetic regime and real dataset")
    text("C denotes the combinatorial method, N the normalized method. Every group is included for each of the three TW protocols.")
    for index, comparison in enumerate(COMPARISONS[1:], start=1):
        blocks = paired_blocks(deltas, comparison=comparison, scope="group")
        rows = [score_row(block, group_label(block.iloc[0])) for block in blocks]
        table(comparison + ": all 22 groups", SCORE_HEADERS, rows,
              f"tab:normalized-group-{index}", landscape=True)
        ci_rows = [ci_row(block, group_label(block.iloc[0])) for block in blocks if block.iloc[0]["suite"] == "synthetic"]
        table(comparison + ": paired synthetic intervals and instance win rates", CI_HEADERS, ci_rows,
              f"tab:normalized-group-ci-{index}", landscape=True)

    heading("Adaptive alpha selections on real matrices")
    rows = []
    for dataset in DATASET_ORDER:
        subset = protocols.loc[(protocols["suite"] == "real") & (protocols["dataset"] == dataset)]
        c_row = subset.loc[subset["variant"] == "c_tw_adaptive"].iloc[0]
        n_row = subset.loc[subset["variant"] == "n_tw_adaptive"].iloc[0]
        rows.append([dataset, alpha_fmt(c_row["selected_alpha"]), alpha_fmt(n_row["selected_alpha"])])
    table("Selected alpha under the same adaptive grid", ("Dataset", "Combinatorial alpha", "Normalized alpha"),
          rows, "tab:normalized-selected-alpha")
    text("Selections use the stated MWB-AUC tie policy. A different selected alpha is part of the adaptive variant and does not by itself identify an effect of degree normalization at a common alpha.")

    heading("Numerical audit and ordering identifiability")
    audit_labels = {
        "component_solves": "Normalized component eigensolves",
        "max_generalized_relative_residual": "Maximum generalized eigenpair relative residual",
        "max_null_mode_residual": "Maximum normalized null-mode residual",
        "near_repeated_fiedler_components": "Component solves with Fiedler gap <=1e-10",
        "flat_nontrivial_row_components": "Component solves with nontrivial flat row coordinates",
        "flat_nontrivial_col_components": "Component solves with nontrivial flat column coordinates",
    }
    audit_rows = []
    for key, value in metadata.get("diagnostics", {}).items():
        audit_rows.append([audit_labels.get(key, key), f"{value:.6e}" if isinstance(value, float) else str(value)])
    control_labels = {
        "replayed_case_alpha_count": "Replayed combinatorial TW matrix-alpha controls",
        "maximum_metric_absolute_delta": "Maximum metric delta from previous alpha sweep",
        "metric_cells_above_1e_minus_10": "Replayed metric cells with delta >1e-10",
    }
    for key, label in control_labels.items():
        if key in metadata.get("control_replay", {}):
            value = metadata["control_replay"][key]
            audit_rows.append([label, f"{value:.6e}" if isinstance(value, float) else str(value)])
    table("Solver diagnostics and combinatorial replay", ("Audit quantity", "Value"), audit_rows,
          "tab:normalized-numerical-audit")
    text(
        "Small eigenpair residuals verify the implemented spectral equations, but do not establish a unique or "
        "stable ordering. A repeated Fiedler eigenvalue gives an eigenspace rather than a unique direction; a "
        "constant or nearly constant within-row/within-column coordinate does not identify an ordering on that "
        "side. Sign conventions and stable sorting do not resolve eigenspace or near-tie ambiguity. The diagnostic "
        "CSV records the eigenvalue gaps and relative within-side coordinate spreads, and the raw CSV preserves "
        "the resulting permutations. The reported metrics describe these recorded deterministic solver outputs; "
        "they should not be interpreted as proving robustness to perturbations in such ambiguous cases. "
        "The paired confidence intervals describe synthetic seed variation conditional on this implementation, "
        "and a favorable point estimate alone is not a claim of a statistically resolved gain.")

    heading("Draft reviewer answer")
    answer = reviewer_answer(summary)
    text(f"Character count: {len(answer)} (including spaces; answer text only).")
    md.extend(["> " + answer, ""])
    tex.extend([r"\begin{quote}", tex_escape(answer), r"\end{quote}", ""])

    heading("Reproducibility and complete artifacts")
    text(
        "The raw CSV stores each input identifier, matrix hash, Laplacian, alpha, all metrics, and row/column "
        "permutations. Component spectral diagnostics and experiment metadata record normalization conventions, "
        "source hashes, package versions, and numerical/control checks. The complete results are generated by "
        "the runner, and this document is regenerated by build_report.py from those outputs.")
    artifacts = [
        ("All 4,298 raw orderings and metrics", "normalized_laplacian_per_instance.csv"),
        ("All 2,456 protocol records", "protocol_comparison.csv"),
        ("Method summary", "method_summary.csv"),
        ("All six-alpha suite summaries", "alpha_summary.csv"),
        ("Paired means, intervals, and win rates", "paired_deltas.csv"),
        ("Component spectral diagnostics", "spectral_diagnostics.csv"),
        ("Control replay differences", "control_regression_differences.csv"),
        ("Experiment metadata", "experiment_metadata.json"),
    ]
    for description, filename in artifacts:
        md.append(f"- [{description}]({filename})")
        tex.append(tex_escape(description + ": " + filename) + r"\par")
    md.append("")
    text("The experiment metadata reports the coordinate convention as " + str(metadata.get("coordinate", "f = D^(-1/2) u"))
         + " and self-loops as " + str(metadata.get("self_loops", "excluded")) + ".")
    (output / "normalized_laplacian_summary.md").write_text("\n".join(md).rstrip() + "\n", encoding="utf-8", newline="\n")
    (output / "normalized_laplacian_summary.tex").write_text("\n".join(tex).rstrip() + "\n", encoding="utf-8", newline="\n")
    print(f"Wrote English MD and single-section LaTeX reports; reviewer answer {len(answer)} characters.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    build(args.output_dir.resolve())


if __name__ == "__main__":
    main()
