# Table 1, synthetic component -- aligned with the manuscript

Same column structure as the real-data tables in `table1_aligned.md`: three
metric groups, seven methods each, plus the selected alpha. O = Original,
M = Marginal, HC = HC+OLO, OW = One-walk, CA = CA-SVD, MD = Median,
TW = Two-walk. Bold marks the best displayed value per row and metric. Rows
follow the manuscript's family order, large to small within each family.

Alpha is chosen per instance, so a regime has one selection per seed. The
alpha column gives the modal value; a parenthesised `k/N` means the seeds did
not agree. `synthetic_per_seed.csv` carries every seed's record.

## A. As submitted

The submitted code: one Fiedler vector over the whole active support, finished with `orient_orders_for_diagonal` and an unstable sort, under the released metric definitions. Should reproduce the `445334d` table. Each row is a mean over the seeds.

| Dataset | Shape | 2-SUM O | 2-SUM M | 2-SUM HC | 2-SUM OW | 2-SUM CA | 2-SUM MD | 2-SUM TW | Band@10% O | Band@10% M | Band@10% HC | Band@10% OW | Band@10% CA | Band@10% MD | Band@10% TW | MWB-AUC O | MWB-AUC M | MWB-AUC HC | MWB-AUC OW | MWB-AUC CA | MWB-AUC MD | MWB-AUC TW | alpha |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| A: Clean one-to-one | 400x360 | 16.70 | 13.84 | 10.36 | 0.50 | 0.59 | 0.79 | **0.43** | 0.195 | 0.214 | 0.237 | 0.866 | 0.849 | 0.765 | **0.902** | 0.436 | 0.473 | 0.530 | 0.924 | 0.918 | 0.890 | **0.933** | 2 (5/20) |
| A: Clean one-to-one | 200x180 | 16.88 | 13.72 | 10.67 | 0.57 | 0.63 | 0.89 | **0.43** | 0.201 | 0.226 | 0.260 | 0.857 | 0.860 | 0.775 | **0.914** | 0.438 | 0.484 | 0.542 | 0.921 | 0.919 | 0.893 | **0.938** | 1 (5/20) |
| A: Clean one-to-one | 100x90 | 17.09 | 13.25 | 9.14 | 0.63 | 0.79 | 0.96 | **0.50** | 0.211 | 0.240 | 0.342 | 0.878 | 0.846 | 0.883 | **0.919** | 0.442 | 0.494 | 0.603 | 0.922 | 0.912 | 0.914 | **0.937** | 12 (6/20) |
| B: Paired subgroup overlap | 400x360 | 16.72 | 14.04 | 8.10 | 1.03 | 1.11 | 1.08 | **0.86** | 0.194 | 0.213 | 0.314 | 0.684 | 0.675 | 0.693 | **0.733** | 0.435 | 0.472 | 0.595 | 0.872 | 0.868 | 0.871 | **0.887** | 12 (8/20) |
| B: Paired subgroup overlap | 200x180 | 16.76 | 14.15 | 8.53 | 1.06 | 1.20 | 1.34 | **0.81** | 0.200 | 0.216 | 0.271 | 0.697 | 0.698 | 0.717 | **0.777** | 0.438 | 0.472 | 0.558 | 0.874 | 0.869 | 0.867 | **0.897** | 12 (10/20) |
| B: Paired subgroup overlap | 100x90 | 16.90 | 13.68 | 10.39 | 1.03 | 1.32 | 1.93 | **0.93** | 0.216 | 0.238 | 0.317 | 0.759 | 0.715 | 0.726 | **0.796** | 0.445 | 0.491 | 0.554 | 0.886 | 0.870 | 0.857 | **0.897** | 1 (6/20) |
| C: Shared super-prototype | 400x360 | 16.92 | 14.07 | 8.71 | 1.15 | 1.35 | 1.78 | **0.92** | 0.195 | 0.216 | 0.267 | 0.669 | 0.643 | 0.687 | **0.726** | 0.434 | 0.475 | 0.575 | 0.865 | 0.853 | 0.850 | **0.884** | 12 (9/20) |
| C: Shared super-prototype | 200x180 | 16.88 | 13.62 | 9.15 | 1.01 | 1.24 | 2.07 | **0.83** | 0.200 | 0.229 | 0.306 | 0.721 | 0.686 | 0.723 | **0.778** | 0.438 | 0.485 | 0.585 | 0.880 | 0.866 | 0.854 | **0.897** | 4 (5/20) |
| C: Shared super-prototype | 100x90 | 17.24 | 13.12 | 7.44 | 1.01 | 1.39 | 3.23 | **0.85** | 0.209 | 0.250 | 0.371 | 0.770 | 0.719 | 0.732 | **0.826** | 0.438 | 0.502 | 0.641 | 0.890 | 0.868 | 0.830 | **0.906** | 4 (9/20) |
| D: Shared prototype with noise | 400x360 | 16.63 | 15.12 | 9.85 | 1.28 | 1.54 | 1.18 | **1.00** | 0.196 | 0.218 | 0.238 | 0.628 | 0.609 | 0.644 | **0.702** | 0.437 | 0.466 | 0.521 | 0.852 | 0.841 | 0.861 | **0.876** | 1 (5/20) |
| D: Shared prototype with noise | 200x180 | 16.74 | 14.35 | 8.94 | 1.24 | 1.67 | 1.51 | **1.04** | 0.200 | 0.234 | 0.294 | 0.689 | 0.623 | 0.649 | **0.724** | 0.439 | 0.480 | 0.568 | 0.866 | 0.839 | 0.853 | **0.880** | 12 (5/20) |
| D: Shared prototype with noise | 100x90 | 17.07 | 13.51 | 7.70 | 1.13 | 1.65 | 2.67 | **0.96** | 0.207 | 0.254 | 0.402 | 0.739 | 0.670 | 0.723 | **0.789** | 0.439 | 0.502 | 0.641 | 0.880 | 0.851 | 0.838 | **0.895** | 4 (4/20) |
| E: Cross-block leakage | 400x360 | 16.65 | 14.70 | 9.46 | 2.04 | 2.09 | 2.65 | **1.91** | 0.195 | 0.216 | 0.260 | 0.551 | 0.568 | **0.605** | 0.554 | 0.436 | 0.469 | 0.560 | 0.812 | 0.815 | 0.816 | **0.817** | 4 (7/20) |
| E: Cross-block leakage | 200x180 | 16.88 | 13.88 | 7.59 | 1.82 | 2.16 | 3.28 | **1.63** | 0.201 | 0.232 | 0.382 | 0.602 | 0.586 | 0.607 | **0.635** | 0.438 | 0.488 | 0.628 | 0.831 | 0.817 | 0.799 | **0.843** | 12 (7/20) |
| E: Cross-block leakage | 100x90 | 16.87 | 13.15 | 8.36 | **1.58** | 2.12 | 4.25 | **1.58** | 0.216 | 0.251 | 0.331 | 0.685 | 0.643 | 0.671 | **0.700** | 0.443 | 0.503 | 0.610 | 0.855 | 0.833 | 0.794 | **0.857** | 12 (8/20) |

## B. Corrected

Section 3.3 component-wise ordering, alpha chosen by the released MWB-AUC, and the paper-faithful metric definitions (Eq. 16 Band@10%, Section 4.5 block-cut MWB-AUC). The synthetic active supports are all connected, so the component rule itself changes nothing here; the movement comes from the orientation and stable-sort fixes the paper-spec commit bundled with it. Each row is a mean over the seeds.

| Dataset | Shape | 2-SUM O | 2-SUM M | 2-SUM HC | 2-SUM OW | 2-SUM CA | 2-SUM MD | 2-SUM TW | Band@10% O | Band@10% M | Band@10% HC | Band@10% OW | Band@10% CA | Band@10% MD | Band@10% TW | MWB-AUC O | MWB-AUC M | MWB-AUC HC | MWB-AUC OW | MWB-AUC CA | MWB-AUC MD | MWB-AUC TW | alpha |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| A: Clean one-to-one | 400x360 | 16.70 | 13.84 | 10.36 | 0.56 | 0.59 | 0.79 | **0.45** | 0.190 | 0.209 | 0.230 | 0.830 | 0.837 | 0.743 | **0.886** | 0.137 | 0.149 | 0.166 | 0.494 | 0.529 | 0.419 | **0.543** | 1 (7/20) |
| A: Clean one-to-one | 200x180 | 16.88 | 13.72 | 10.67 | 0.57 | 0.63 | 0.89 | **0.45** | 0.190 | 0.215 | 0.247 | 0.835 | 0.838 | 0.736 | **0.897** | 0.136 | 0.155 | 0.175 | 0.520 | 0.529 | 0.437 | **0.568** | 4 (6/20) |
| A: Clean one-to-one | 100x90 | 17.09 | 13.25 | 9.14 | 0.54 | 0.79 | 0.96 | **0.53** | 0.191 | 0.218 | 0.312 | 0.870 | 0.806 | 0.838 | **0.875** | 0.137 | 0.162 | 0.214 | **0.552** | 0.521 | 0.513 | **0.552** | 12 (7/20) |
| B: Paired subgroup overlap | 400x360 | 16.72 | 14.04 | 8.10 | 1.00 | 1.11 | 1.08 | **0.86** | 0.189 | 0.207 | 0.306 | 0.679 | 0.660 | 0.675 | **0.720** | 0.137 | 0.149 | 0.208 | 0.419 | 0.422 | 0.404 | **0.439** | 12 (9/20) |
| B: Paired subgroup overlap | 200x180 | 16.76 | 14.15 | 8.53 | **0.94** | 1.20 | 1.34 | 0.95 | 0.189 | 0.204 | 0.256 | **0.711** | 0.668 | 0.681 | 0.706 | 0.137 | 0.149 | 0.178 | **0.440** | 0.428 | 0.405 | 0.434 | 12 (8/20) |
| B: Paired subgroup overlap | 100x90 | 16.90 | 13.68 | 10.39 | 1.13 | 1.32 | 1.93 | **1.07** | 0.193 | 0.215 | 0.290 | 0.670 | 0.663 | 0.669 | **0.712** | 0.141 | 0.158 | 0.202 | 0.414 | **0.436** | 0.415 | **0.436** | 1 (5/20) |
| C: Shared super-prototype | 400x360 | 16.92 | 14.07 | 8.71 | 1.12 | 1.35 | 1.78 | **0.91** | 0.189 | 0.210 | 0.258 | 0.654 | 0.630 | 0.670 | **0.703** | 0.137 | 0.152 | 0.178 | 0.402 | 0.407 | 0.404 | **0.433** | 8 (6/20) |
| C: Shared super-prototype | 200x180 | 16.88 | 13.62 | 9.15 | 1.00 | 1.24 | 2.07 | **0.85** | 0.190 | 0.218 | 0.291 | 0.698 | 0.658 | 0.688 | **0.744** | 0.137 | 0.158 | 0.207 | 0.433 | 0.424 | 0.410 | **0.460** | 12 (7/20) |
| C: Shared super-prototype | 100x90 | 17.24 | 13.12 | 7.44 | 1.00 | 1.39 | 3.23 | **0.90** | 0.190 | 0.225 | 0.341 | 0.712 | 0.663 | 0.675 | **0.763** | 0.137 | 0.166 | 0.233 | 0.441 | 0.440 | 0.397 | **0.471** | 1 (6/20) |
| D: Shared prototype with noise | 400x360 | 16.63 | 15.12 | 9.85 | 1.23 | 1.54 | 1.18 | **0.93** | 0.191 | 0.213 | 0.231 | 0.623 | 0.595 | 0.626 | **0.707** | 0.137 | 0.157 | 0.164 | 0.384 | 0.394 | 0.380 | **0.431** | 1 (6/20) |
| D: Shared prototype with noise | 200x180 | 16.74 | 14.35 | 8.94 | 1.17 | 1.67 | 1.51 | **1.06** | 0.190 | 0.223 | 0.279 | 0.651 | 0.593 | 0.611 | **0.683** | 0.137 | 0.165 | 0.190 | 0.400 | 0.396 | 0.368 | **0.422** | 12 (9/20) |
| D: Shared prototype with noise | 100x90 | 17.07 | 13.51 | 7.70 | 1.05 | 1.65 | 2.67 | **0.95** | 0.188 | 0.230 | 0.367 | 0.701 | 0.619 | 0.663 | **0.737** | 0.137 | 0.175 | 0.245 | 0.430 | 0.418 | 0.410 | **0.450** | 12 (5/20) |
| E: Cross-block leakage | 400x360 | 16.65 | 14.70 | 9.46 | **1.77** | 2.09 | 2.65 | 1.85 | 0.190 | 0.211 | 0.253 | 0.587 | 0.554 | **0.589** | 0.550 | 0.136 | 0.154 | 0.179 | **0.369** | 0.367 | 0.368 | 0.342 | 12 (15/20) |
| E: Cross-block leakage | 200x180 | 16.88 | 13.88 | 7.59 | 1.76 | 2.16 | 3.28 | **1.62** | 0.190 | 0.222 | 0.366 | 0.579 | 0.557 | 0.575 | **0.605** | 0.137 | 0.164 | 0.244 | 0.361 | **0.373** | 0.356 | 0.371 | 12 (12/20) |
| E: Cross-block leakage | 100x90 | 16.87 | 13.15 | 8.36 | **1.52** | 2.12 | 4.25 | 1.61 | 0.194 | 0.224 | 0.305 | 0.632 | 0.588 | 0.618 | **0.644** | 0.139 | 0.171 | 0.213 | 0.394 | **0.397** | 0.390 | 0.393 | 12 (7/20) |

## C. Reproduction check

### A. As submitted vs the `445334d` CSV

Means over the seeds, compared against the submitted table's two-decimal
values.

Cells that differ:

- A: Clean one-to-one / Large | One-walk | 2-SUM: ours 0.5046 vs reference 0.5200
- A: Clean one-to-one / Large | One-walk | Band@10%: ours 0.8658 vs reference 0.8600
- A: Clean one-to-one / Medium | One-walk | 2-SUM: ours 0.5706 vs reference 0.5500
- A: Clean one-to-one / Medium | One-walk | Band@10%: ours 0.8572 vs reference 0.8700
- A: Clean one-to-one / Medium | TW | 2-SUM: ours 0.4348 vs reference 0.4200
- A: Clean one-to-one / Medium | TW | Band@10%: ours 0.9139 vs reference 0.9200
- A: Clean one-to-one / Small | One-walk | 2-SUM: ours 0.6305 vs reference 0.6000
- A: Clean one-to-one / Small | One-walk | Band@10%: ours 0.8783 vs reference 0.8900
- A: Clean one-to-one / Small | One-walk | MWB-AUC: ours 0.9219 vs reference 0.9300
- A: Clean one-to-one / Small | TW | 2-SUM: ours 0.4990 vs reference 0.5100
- B: Paired subgroup overlap / Large | One-walk | 2-SUM: ours 1.0331 vs reference 1.0400
- B: Paired subgroup overlap / Large | TW | 2-SUM: ours 0.8608 vs reference 0.8800
- B: Paired subgroup overlap / Large | TW | Band@10%: ours 0.7326 vs reference 0.7400
- B: Paired subgroup overlap / Medium | One-walk | 2-SUM: ours 1.0615 vs reference 1.1500
- B: Paired subgroup overlap / Medium | One-walk | Band@10%: ours 0.6966 vs reference 0.6700
- B: Paired subgroup overlap / Medium | TW | 2-SUM: ours 0.8116 vs reference 0.8200
- B: Paired subgroup overlap / Small | TW | 2-SUM: ours 0.9308 vs reference 0.9400
- B: Paired subgroup overlap / Small | TW | Band@10%: ours 0.7962 vs reference 0.7900
- C: Shared super-prototype / Large | One-walk | 2-SUM: ours 1.1458 vs reference 1.1300
- C: Shared super-prototype / Large | One-walk | MWB-AUC: ours 0.8646 vs reference 0.8700
- C: Shared super-prototype / Large | TW | 2-SUM: ours 0.9226 vs reference 0.9300
- C: Shared super-prototype / Medium | TW | 2-SUM: ours 0.8343 vs reference 0.8600
- C: Shared super-prototype / Medium | TW | Band@10%: ours 0.7778 vs reference 0.7700
- C: Shared super-prototype / Small | TW | 2-SUM: ours 0.8527 vs reference 0.8600
- C: Shared super-prototype / Small | TW | MWB-AUC: ours 0.9056 vs reference 0.9000
- D: Shared prototype with noise / Medium | One-walk | 2-SUM: ours 1.2441 vs reference 1.2300
- D: Shared prototype with noise / Small | One-walk | 2-SUM: ours 1.1333 vs reference 1.1400
- D: Shared prototype with noise / Small | TW | 2-SUM: ours 0.9571 vs reference 0.9400
- D: Shared prototype with noise / Small | TW | MWB-AUC: ours 0.8947 vs reference 0.9000
- E: Cross-block leakage / Large | One-walk | 2-SUM: ours 2.0428 vs reference 2.0900
- E: Cross-block leakage / Large | TW | 2-SUM: ours 1.9054 vs reference 1.8700
- E: Cross-block leakage / Large | TW | Band@10%: ours 0.5539 vs reference 0.5600
- E: Cross-block leakage / Medium | One-walk | 2-SUM: ours 1.8241 vs reference 1.8600
- E: Cross-block leakage / Medium | TW | 2-SUM: ours 1.6336 vs reference 1.6500
- E: Cross-block leakage / Small | One-walk | 2-SUM: ours 1.5835 vs reference 1.6100
- E: Cross-block leakage / Small | One-walk | Band@10%: ours 0.6852 vs reference 0.6700
- E: Cross-block leakage / Small | TW | 2-SUM: ours 1.5842 vs reference 1.6000

### B. Corrected vs the current CSV

MWB-AUC is excluded: Table B reports the Section 4.5 block-cut quantity
while the CSV carries the selector's band integral.

All checked cells reproduce the current CSV to within 5e-3.

