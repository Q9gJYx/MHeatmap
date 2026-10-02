# Alpha sensitivity experiment

This experiment uses the same paper-spec, component-wise Two-Walk implementation
as Table 1. Each positive-alpha solve delegates the component Laplacian to
`mheatmap.graph.two_walk_laplacian` from `mheatmap==1.2.5`.

Run from the repository root:

```powershell
.\.venv\Scripts\python.exe main_experiment\alpha_sensitivity\run_alpha_sensitivity.py
```

The positive grid is
`{0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1, 2, 4, 6, 8, 12, 24}`.
The submitted Table 1 grid remains `{1, 2, 4, 6, 8, 12}`. Alpha zero is not
assigned an ordering: the script verifies its two zero modes per connected
support component and records the result in `alpha_zero_diagnostics.csv`.

The single global alpha is selected once by maximum mean MWB-AUC over the 300
synthetic instances, with smaller-alpha tie breaking, and then frozen for the
seven real matrices. Outputs are written to
`output/main_experiment/alpha_sensitivity/`.

Key files:

- `alpha_sensitivity_per_instance.csv`: 307 cases x 13 positive alpha values.
- `alpha_sensitivity_by_group.csv`: per-regime and per-dataset summaries.
- `alpha_sensitivity_overall.csv`: synthetic/real sensitivity curves.
- `best_alpha_per_instance.csv`: submitted-grid and expanded-grid selections.
- `best_alpha_summary.csv`: selection distributions and means.
- `fixed_and_adaptive_comparison.csv`: OW, fixed-alpha, and adaptive comparison.
- `alpha_zero_diagnostics.csv`: numerical zero-mode audit.
- `figures/alpha_sensitivity.{png,pdf}`: sensitivity figure.
