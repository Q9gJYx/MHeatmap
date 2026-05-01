# Synthetic rectangular benchmark

This benchmark procedurally generates five synthetic rectangular matrix
families at three sizes and averages results over 20 pre-specified random
seeds per family-size regime.

The evaluated methods are:

- Original order (`O`)
- Marginal-total ordering (`M`)
- Hierarchical clustering with optimal leaf ordering (`HC`)
- One-walk bipartite spectral ordering (`OW`)
- CA-SVD co-ordering (`CA`)
- Alternating weighted median ordering (`MD`)
- Adaptive Two-Walk (`TW`)

Adaptive TW selects `alpha` from `{1, 2, 4, 6, 8, 12}` independently for each
matrix instance by maximizing internal MWB-AUC.

## Requirements

Install the paper-benchmark dependencies and the public `mheatmap` package:

```bash
uv pip install numpy pandas scipy scikit-learn matplotlib
uv pip install git+https://github.com/qqgjyx/mheatmap.git
```

## Reproduce

From the repository root:

```bash
uv run python main_experiment/synthetic_benchmark/run_synthetic_strong_baselines.py
```

This regenerates:

- `output/main_experiment/synthetic_benchmark/synthetic_strong_baselines_paper.csv`
- `output/main_experiment/synthetic_benchmark/synthetic_strong_baselines_paper.md`

These two files contain the machine-readable and Markdown versions of the
synthetic component of the main paper table.

