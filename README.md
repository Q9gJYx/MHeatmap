# MHeatmap — Reproducibility Companion

The MHeatmap method itself is published as the
[`mheatmap`](https://github.com/qqgjyx/mheatmap) Python package
(PyPI v1.2.5). This repository contains only the benchmark, preprocessing,
and figure code specific to the paper.

## Repository layout

```
.
├── main_experiment/
│   ├── synthetic_benchmark/    # Families A–E, 3 sizes, 20 seeds (Table 1 top)
│   ├── real_world_benchmark/   # 7 real datasets (Teaser, Table 1 bottom)
│   ├── alpha_sensitivity/      # Fixed alpha and adaptive-selection controls
│   └── normalized_laplacian/   # Degree-normalized final-Laplacian ablation
├── examples/
│   └── tabula_sapiens/         # Fig 3 case study (4 tissues × 3 panels)
├── output/
│   └── main_experiment/        # Versioned paper-table artifacts (CSV / Markdown)
└── data/                       # Reserved for locally cached raw data (not committed)
```

## Install

Requires [uv](https://docs.astral.sh/uv/).

```bash
uv sync # Execute under project root to sync all dependencies.
```

The version pin matches the package version used for every result in the paper.

## Experiments, benchmarks, and case studies

### Reproduce Table 1

**Synthetic component (top)** — 5 families × 3 sizes × 20 seeds = 300 instances:

```bash
uv run python main_experiment/synthetic_benchmark/run_synthetic_strong_baselines.py
```

Writes
`output/main_experiment/synthetic_benchmark/synthetic_strong_baselines_paper.{csv,md}`.

**Real-world component (bottom)** — 7 rectangular matrices:

```bash
uv run python main_experiment/real_world_benchmark/run_real_world_benchmark.py
```

Writes
`output/main_experiment/real_world_benchmark/real_world_benchmark_paper.{csv,md}`
plus the seven cleaned matrix CSVs under `processed_matrices/`.

Both runs evaluate the same six baselines plus adaptive Two-Walk:
Original (`O`), Marginal-total (`M`), Hierarchical clustering with optimal leaf
ordering (`HC`), One-walk bipartite spectral (`OW`), CA-SVD (`CA`),
Alternating weighted median (`MD`), and adaptive Two-Walk (`TW`). Adaptive TW
selects α from `{1, 2, 4, 6, 8, 12}` per matrix by maximizing internal MWB-AUC.

#### Corrected metric and ordering definitions

- **Band@10%** implements Equation (16): row and column ranks are mapped to
  `i / (m - 1)` and `j / (n - 1)`, and the score is the fraction of total mass
  with normalized coordinate distance at most `0.10`.
- **MWB-AUC** retains the submitted alpha selector. It is the normalized
  trapezoidal integral of the legacy integer-index diagonal-band mass over 25
  widths from `0.02` through `0.50`; it is not a matched-block-cut metric.
- **OW/TW component handling** follows the manuscript specification. The
  active matrix is normalized once by its global maximum, each positive-support
  component is ordered independently, components are concatenated by decreasing
  raw mass, and zero-marginal axes are appended in original order. Every
  positive alpha, including `alpha=1`, follows the same TW path.
- Each TW component delegates its Laplacian construction to
  `mheatmap.graph.two_walk_laplacian` from the pinned package. Component
  discovery, global normalization, deterministic Fiedler orientation, stable
  sorting, and concatenation are benchmark-level responsibilities.

Run the paper-specification regression tests with:

```bash
uv run python -m unittest discover -s tests -v
```

### Normalized-Laplacian comparison (JfQK Q2)

This paired ablation keeps the Two-Walk adjacency and the six-point alpha grid
unchanged. It compares the combinatorial Laplacian with loop-free symmetric
normalization, sorting the generalized coordinate `f = D^(-1/2) u`.
Fixed `alpha=1`, the previously selected global `alpha=12`, and equal-grid
adaptive selection are evaluated on all 300 synthetic and seven real matrices.

```bash
uv run python main_experiment/normalized_laplacian/run_normalized_laplacian.py
uv run python main_experiment/normalized_laplacian/build_report.py
```

The [experiment README](main_experiment/normalized_laplacian/README.md) specifies
the formulas, coordinate and self-loop conventions, paired intervals, and
permutation/provenance artifacts. Results supplement Table 1 rather than
creating a second Table 1 code path.

### Reproduce qualitative figures

**Figure 1 (teaser) and the real-world four-panel figures:**

```bash
uv run python main_experiment/real_world_benchmark/make_real_world_quad_figures.py
```

Writes seven PNGs under
`output/main_experiment/real_world_benchmark/figures/`.

**Figure 3 (Tabula Sapiens case study):**

```bash
uv run python examples/tabula_sapiens/tabula_sapiens_v2.py
```

Per-tissue PNGs (raw / spectral / mosaic / combined / pareto / umap) are
checked into `examples/tabula_sapiens/<tissue>/` for browsing convenience.

### Data sources

The **synthetic suite** is generated procedurally by code; no raw synthetic data is
required.

The seven **real-world matrices** are built from public sources:

| Dataset           | Source                                                                                                                                                                   |
| ----------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| **SIC → NAICS**   | [1987 SIC → 2002 NAICS concordance](https://www2.census.gov/library/reference/naics/technical-documentation/concordance/1987_sic_to_2002_naics.xls) (U.S. Census Bureau) |
| **CIP → SOC**     | [CIP 2020 → SOC 2018 crosswalk](https://nces.ed.gov/ipeds/cipcode/Files/CIP2020_SOC2018_Crosswalk.xlsx) (NCES / BLS)                                                     |
| **ACS PUMS**      | [ACS 2023 1-Year PUMS, Massachusetts](https://www2.census.gov/programs-surveys/acs/data/pums/2023/1-Year/csv_pma.zip) (U.S. Census Bureau)                               |
| **LODES**         | [LODES8 2022 Tennessee OD + crosswalk](https://lehd.ces.census.gov/data/lodes/LODES8/tn/) (U.S. Census LEHD)                                                             |
| **20 Newsgroups** | scikit-learn's `fetch_20newsgroups`; the fixed document and term lists for the paper's matrix are hard-coded in the preprocessing script.                                |
| **MBTA**          | [Static MBTA GTFS feed](https://cdn.mbta.com/MBTA_GTFS.zip)                                                                                                              |
| **OpenAlex**      | [OpenAlex `/authors` API](https://api.openalex.org/authors); the exact author-topic query is hard-coded in the preprocessing script.                                     |

The full preprocessing pipeline lives in
`main_experiment/real_world_benchmark/_real_world_data.py`.

The **Tabula Sapiens v2 single-cell atlas** used for Figure 3 is downloaded from
[figshare 27921984](https://figshare.com/articles/dataset/Tabula_Sapiens_v2/27921984).
Update `DATA_DIR` at the top of `examples/tabula_sapiens/tabula_sapiens_v2.py`
to point at your local copy of the four `.h5ad` files
(`bone_marrow`, `spleen`, `stomach`, `vasculature`).

## License

MIT. See `LICENSE`.
