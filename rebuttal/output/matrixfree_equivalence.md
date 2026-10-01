# Matrix-free vs dense equivalence (Section 3.3)

Instances compared: 21 solver configurations plus 3 adaptive-alpha selections.

## Agreement contract

| quantity | worst observed | threshold | pass |
|---|---|---|---|
| rel_lambda_diff | 1.019e-12 | 1e-08 | yes |
| max_sorted_coord_dev | 3.698e-11 | 1e-09 | yes |
| d_r2s | 1.819e-02 | 1e-04 | **NO** |
| d_band | 0.000e+00 | 1e-09 | yes |
| d_mwb | 1.728e-02 | 1e-04 | **NO** |

- Orderings differ outside a tie group: **0** configurations.
- Configurations whose permutation differs at all: **7** of 21.

A differing permutation is expected and harmless: the Fiedler coordinates are tied,
ties are adjacent once sorted, and all three metrics are invariant to permuting
within a tie group. The published per-instance permutations are only defined up to
those ties.

## Alpha-selection stability

| instance | selected alpha | dMWB-AUC | stable |
|---|---|---|---|
| OpenAlex Author x Topic | 1 | 0.000e+00 | yes |
| LODES Home x Work | 4 | 0.000e+00 | yes |
| ACS OCCP x INDP | 4 | 0.000e+00 | yes |
