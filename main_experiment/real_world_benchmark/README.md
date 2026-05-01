# Real-world rectangular benchmark

This benchmark rebuilds the seven real-world rectangular matrices used in the
main paper table, evaluates the same representative reordering methods, and
generates qualitative four-panel figures.

The evaluated methods are:

- Original order (`O`)
- Marginal-total ordering (`M`)
- Hierarchical clustering with optimal leaf ordering (`HC`)
- One-walk bipartite spectral ordering (`OW`)
- CA-SVD co-ordering (`CA`)
- Alternating weighted median ordering (`MD`)
- Adaptive Two-Walk (`TW`)

Adaptive TW selects `alpha` from `{1, 2, 4, 6, 8, 12}` independently for each
input matrix by maximizing internal MWB-AUC.

## Inputs

The seven cleaned matrix CSVs used by the paper table are versioned under
`output/main_experiment/real_world_benchmark/processed_matrices/`. Raw public
sources are not committed; source locations are listed in
`data/main_experiment/real_world_benchmark/README.md`.

The preprocessing code is in `_real_world_data.py`, not in the data directory.

## Reproduce the benchmark table

From the repository root:

```bash
uv run python main_experiment/real_world_benchmark/run_real_world_benchmark.py
```

This regenerates:

- `output/main_experiment/real_world_benchmark/real_world_benchmark_paper.csv`
- `output/main_experiment/real_world_benchmark/real_world_benchmark_paper.md`
- seven cleaned matrix CSVs under
  `output/main_experiment/real_world_benchmark/processed_matrices/`

## Reproduce the four-panel figures

From the repository root:

```bash
uv run python main_experiment/real_world_benchmark/make_real_world_quad_figures.py
```

This regenerates seven PNG figures under:

```text
output/main_experiment/real_world_benchmark/figures/
```
