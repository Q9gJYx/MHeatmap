# Table 1 -- aligned with the manuscript

Column groups are the three metrics; within each group the method order
matches Table 1: O = Original, M = Marginal, HC = HC+OLO, OW = One-walk,
CA = CA-SVD, MD = Median, TW = Two-walk. Bold marks the best value per row
and metric. Row order follows Table 1.

## A. As submitted

The submitted code: one Fiedler vector over the whole active support, finished with `orient_orders_for_diagonal` and an unstable sort, under the released metric definitions. Should reproduce the `445334d` table.

| Dataset | Shape | 2-SUM O | 2-SUM M | 2-SUM HC | 2-SUM OW | 2-SUM CA | 2-SUM MD | 2-SUM TW | Band@10% O | Band@10% M | Band@10% HC | Band@10% OW | Band@10% CA | Band@10% MD | Band@10% TW | MWB-AUC O | MWB-AUC M | MWB-AUC HC | MWB-AUC OW | MWB-AUC CA | MWB-AUC MD | MWB-AUC TW | alpha |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| SIC -> NAICS | 1004x1176 | 2.61 | 10.78 | 2.11 | 0.66 | 2.88 | 0.76 | **0.21** | 0.660 | 0.302 | 0.515 | **0.976** | 0.509 | 0.901 | 0.954 | 0.822 | 0.564 | 0.779 | 0.944 | 0.734 | 0.917 | **0.963** | 4 |
| CIP -> SOC | 2143x868 | 15.90 | 16.04 | 12.18 | 2.52 | 9.77 | 12.37 | **1.68** | 0.190 | 0.268 | 0.205 | 0.435 | 0.305 | 0.165 | **0.619** | 0.444 | 0.500 | 0.533 | 0.768 | 0.573 | 0.401 | **0.830** | 8 |
| ACS OCCP x INDP | 24x10 | 20.67 | 8.87 | 15.94 | 3.39 | **2.23** | 3.05 | 3.61 | 0.262 | 0.518 | 0.319 | 0.875 | **0.910** | 0.850 | 0.881 | 0.518 | 0.735 | 0.534 | 0.921 | **0.946** | 0.910 | 0.922 | 4 |
| LODES Home x Work | 95x95 | 10.98 | 3.98 | 2.78 | **1.89** | 5.42 | 4.54 | 1.90 | 0.601 | 0.723 | 0.693 | 0.775 | **0.792** | 0.766 | 0.776 | 0.690 | 0.824 | 0.815 | 0.871 | 0.858 | 0.831 | **0.871** | 4 |
| 20 Newsgroups | 60x48 | 6.98 | 0.99 | 2.74 | **0.43** | 0.56 | 0.47 | 0.43 | 0.272 | 0.916 | 0.067 | 0.944 | 0.929 | 0.937 | **0.946** | 0.681 | 0.932 | 0.715 | 0.945 | **0.968** | 0.947 | 0.946 | 6 |
| MBTA Route x Station | 8x125 | 17.96 | 9.39 | 3.33 | **1.25** | 3.09 | 2.19 | 1.25 | 0.024 | 0.014 | 0.089 | 0.101 | 0.140 | **0.189** | 0.110 | 0.035 | 0.020 | 0.134 | 0.139 | 0.193 | **0.265** | 0.151 | 1 |
| OpenAlex Author x Topic | 60x8 | 19.72 | 12.88 | 11.79 | **0.74** | 2.53 | 3.36 | 0.79 | 0.394 | 0.556 | 0.512 | 0.971 | 0.884 | 0.978 | **0.996** | 0.555 | 0.675 | 0.675 | 0.984 | 0.940 | 0.986 | **0.996** | 1 |

## B. Corrected

Section 3.3 component-wise ordering for OW and TW, alpha chosen by the released MWB-AUC (the selector the repository actually runs), and the paper-faithful metric definitions (Eq. 16 Band@10%, Section 4.5 block-cut MWB-AUC). Baselines keep their orderings; only their metric values move.

| Dataset | Shape | 2-SUM O | 2-SUM M | 2-SUM HC | 2-SUM OW | 2-SUM CA | 2-SUM MD | 2-SUM TW | Band@10% O | Band@10% M | Band@10% HC | Band@10% OW | Band@10% CA | Band@10% MD | Band@10% TW | MWB-AUC O | MWB-AUC M | MWB-AUC HC | MWB-AUC OW | MWB-AUC CA | MWB-AUC MD | MWB-AUC TW | alpha |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| SIC -> NAICS | 1004x1176 | 2.61 | 10.78 | 2.11 | 0.22 | 2.88 | 0.76 | **0.16** | 0.715 | 0.364 | 0.558 | 0.982 | 0.557 | 0.906 | **0.983** | 0.444 | 0.278 | 0.408 | 0.579 | 0.399 | 0.643 | **0.684** | 1 |
| CIP -> SOC | 2143x868 | 15.90 | 16.04 | 12.18 | 2.81 | 9.77 | 12.37 | **1.95** | 0.190 | 0.267 | 0.205 | 0.377 | 0.305 | 0.164 | **0.558** | 0.134 | 0.215 | 0.189 | 0.251 | 0.268 | 0.139 | **0.323** | 6 |
| ACS OCCP x INDP | 24x10 | 20.67 | 8.87 | 15.94 | 3.39 | **2.23** | 3.05 | 3.61 | 0.138 | 0.291 | 0.150 | 0.431 | **0.768** | 0.765 | 0.339 | 0.077 | 0.277 | 0.135 | 0.283 | **0.581** | 0.563 | 0.251 | 4 |
| LODES Home x Work | 95x95 | 10.98 | 3.98 | 2.78 | **1.89** | 5.42 | 4.54 | 1.90 | 0.601 | 0.723 | 0.693 | 0.775 | **0.792** | 0.766 | 0.776 | 0.579 | 0.617 | 0.393 | 0.546 | 0.692 | **0.693** | 0.528 | 4 |
| 20 Newsgroups | 60x48 | 6.98 | 0.99 | 2.74 | **0.32** | 0.56 | 0.47 | 0.42 | 0.272 | 0.912 | 0.063 | 0.943 | 0.926 | 0.929 | **0.943** | 0.265 | 0.848 | 0.248 | 0.678 | 0.722 | **0.861** | 0.533 | 12 |
| MBTA Route x Station | 8x125 | 17.96 | 9.39 | 3.33 | **1.25** | 3.09 | 2.19 | 1.25 | 0.227 | 0.252 | 0.577 | 0.874 | 0.517 | 0.606 | **0.875** | 0.169 | 0.219 | 0.361 | **0.551** | 0.389 | 0.388 | **0.551** | 2 |
| OpenAlex Author x Topic | 60x8 | 19.72 | 12.88 | 11.79 | **0.44** | 2.53 | 3.36 | 0.45 | 0.075 | 0.302 | 0.193 | **0.847** | 0.555 | 0.191 | **0.847** | 0.078 | 0.227 | 0.154 | **0.560** | 0.486 | 0.252 | 0.543 | 1 |

## C. Reproduction check

### A. As submitted vs the `445334d` CSV

The submitted code should reproduce its own table. All seven datasets are
checked; a mismatch on One-walk or Two-walk means the single-vector
baseline here has drifted from the released behaviour.

Cells that differ:

- SIC -> NAICS | One-walk | 2-SUM: ours 0.6641 vs reference 0.6800
- SIC -> NAICS | TW | 2-SUM: ours 0.2117 vs reference 0.1800
- SIC -> NAICS | TW | Band@10%: ours 0.9542 vs reference 0.9600
- SIC -> NAICS | TW | MWB-AUC: ours 0.9628 vs reference 0.9700
- CIP -> SOC | One-walk | Band@10%: ours 0.4353 vs reference 0.4300
- CIP -> SOC | TW | 2-SUM: ours 1.6775 vs reference 1.6700
- CIP -> SOC | TW | Band@10%: ours 0.6192 vs reference 0.6100
- OpenAlex Author x Topic | One-walk | 2-SUM: ours 0.7444 vs reference 0.9800
- OpenAlex Author x Topic | One-walk | Band@10%: ours 0.9714 vs reference 0.9900
- OpenAlex Author x Topic | One-walk | MWB-AUC: ours 0.9840 vs reference 0.9900
- OpenAlex Author x Topic | TW | 2-SUM: ours 0.7886 vs reference 0.9200

### B. Corrected vs the current CSV

The corrected ordering pipeline should reproduce the repository's current
table exactly on the shared columns. MWB-AUC is excluded because Table B
reports the Section 4.5 block-cut quantity while the CSV carries the
selector's band integral.

All checked cells reproduce the current CSV to within 5e-3.

## D. Alpha selection under each MWB-AUC definition

Section 4.5 describes MWB-AUC as a matched contiguous block-cut quantity, but the repository's adaptive Two-walk selects alpha with the legacy integer-window band integral. This shows whether adopting the paper-as-written selector would have changed the choice.

| Dataset | As submitted | Corrected, legacy selector | Corrected, block-cut selector |
|---|---:|---:|---:|
| SIC -> NAICS | 4 | 1 | 1 |
| CIP -> SOC | 8 | 6 | 12 |
| ACS OCCP x INDP | 4 | 4 | 2 |
| LODES Home x Work | 4 | 4 | 6 |
| 20 Newsgroups | 6 | 12 | 12 |
| MBTA Route x Station | 1 | 2 | 1 |
| OpenAlex Author x Topic | 1 | 1 | 8 |

