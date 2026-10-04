# Normalized-Laplacian ablation (JfQK Q2)

This experiment tests the reviewer's **normalized Laplacian** alternative.
The submitted Two-Walk adjacency is unchanged; this is not a normalized-Gram
or cosine-projection experiment.

For the globally maximum-scaled matrix, construct the same
`A_alpha = M^2 + alpha M` using `mheatmap==1.2.5` inside each support component.
Let `L = diag(A_alpha @ 1) - A_alpha`. Gram self loops cancel in `L`.
Normalize using the **loop-free** degrees `d = diag(L)`:

$$
L_{\mathrm{sym}}=D^{-1/2}LD^{-1/2},\qquad D=\operatorname{diag}(d).
$$

Compute the second eigenvector `u` of `L_sym` and sort the generalized / random
walk coordinate `f = D^(-1/2) u`, so `L f = lambda D f`. This single coordinate
convention is fixed before examining results. Self-loop degrees and sorting
`u` directly would define different variants.

All other handling follows the existing Table 1 implementation: one global
maximum scaling, independent support components, decreasing component mass,
deterministic sign, stable sorting, and trailing zero-marginal axes. Evaluate
R2S, Equation (16) Band@10%, and retained MWB-AUC on the reordered **original
matrix**, not the normalized adjacency.

The paired experiment recomputes combinatorial and normalized OW, and both TW
versions for `alpha in {1,2,4,6,8,12}`, on the same 300 synthetic instances and
seven cached real matrices. Report fixed `alpha=1`, the previously
synthetic-selected and frozen `alpha=12`, and per-matrix adaptive selection
using the same six-point grid and Table 1 MWB-AUC tie policy. MWB-AUC is
selection-aligned for the adaptive comparison. The frozen 12 is not retuned
for the normalized method.

Run from the repository root:

```powershell
.\.venv\Scripts\python.exe main_experiment\normalized_laplacian\run_normalized_laplacian.py
.\.venv\Scripts\python.exe main_experiment\normalized_laplacian\build_report.py
```

Alternatively use `uv run python` with the same script paths. Rebuild
aggregates without repeating the eigensolves using the runner's
`--summarize-only` flag. Outputs and Matplotlib cache remain within this
repository, under `output/main_experiment/normalized_laplacian/` and `.cache/`.

Artifacts:

- `normalized_laplacian_per_instance.csv`: 4,298 method/case/alpha records,
  including matrix hashes and row/column permutations.
- `protocol_comparison.csv`: 2,456 case/protocol records (eight protocols).
- `method_summary.csv`: aggregate OW, fixed-alpha, and adaptive results.
- `alpha_summary.csv`: all six alpha values for both suites and Laplacians.
- `paired_deltas.csv`: normalized-minus-combinatorial differences. Synthetic
  95% percentile bootstrap intervals use 10,000 paired resamples over seeds;
  the overall synthetic interval resamples seed blocks jointly across all
  15 family/size regimes, preserving their shared random-number streams.
  They describe conditional seed variation, not arbitrary matrix populations.
  Single real matrices have no resampling CIs; the real mean is a seven-dataset
  macro average.
- `spectral_diagnostics.csv`: eigenvalues, residuals, degree ranges and
  within-side coordinate spread for every normalized component solve.
- `control_regression_differences.csv`: control metric changes above `1e-10`
  relative to the previous alpha experiment, if any.
- `experiment_metadata.json`: protocol, software/source hashes, diagnostics,
  and regression audit.
- `figures/normalized_laplacian_comparison.{png,pdf}`: overall comparison.
- `normalized_laplacian_summary.{md,tex}`: full English analysis and draft
  rebuttal; the LaTeX fragment lives under one section.

These ablation artifacts supplement the canonical Table 1; its existing files
and runners remain the sole Table 1 generation path.
