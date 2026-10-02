# Synthetic Rectangular Benchmark (Paper Table)

Best $\alpha$ maximizes MWB-AUC on the submitted grid; exact ties use the smaller value. Synthetic entries are means over 20 per-instance selections.

| Family | Size | Mean best $\alpha$ | $10^2 \times$ 2-SUM $\downarrow$ O | $10^2 \times$ 2-SUM $\downarrow$ M | $10^2 \times$ 2-SUM $\downarrow$ HC | $10^2 \times$ 2-SUM $\downarrow$ OW | $10^2 \times$ 2-SUM $\downarrow$ CA | $10^2 \times$ 2-SUM $\downarrow$ MD | $10^2 \times$ 2-SUM $\downarrow$ TW | Band@10\% $\uparrow$ O | Band@10\% $\uparrow$ M | Band@10\% $\uparrow$ HC | Band@10\% $\uparrow$ OW | Band@10\% $\uparrow$ CA | Band@10\% $\uparrow$ MD | Band@10\% $\uparrow$ TW | MWB-AUC $\uparrow$ O | MWB-AUC $\uparrow$ M | MWB-AUC $\uparrow$ HC | MWB-AUC $\uparrow$ OW | MWB-AUC $\uparrow$ CA | MWB-AUC $\uparrow$ MD | MWB-AUC $\uparrow$ TW |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Family A: Clean one-to-one | Large | 3.45 | 16.70 | 13.84 | 10.36 | 0.56 | 0.59 | 0.79 | 0.45 | 0.19 | 0.21 | 0.23 | 0.83 | 0.84 | 0.74 | 0.89 | 0.44 | 0.47 | 0.53 | 0.91 | 0.92 | 0.89 | 0.93 |
| Family A: Clean one-to-one | Medium | 4.35 | 16.88 | 13.72 | 10.67 | 0.57 | 0.63 | 0.89 | 0.45 | 0.19 | 0.22 | 0.25 | 0.84 | 0.84 | 0.74 | 0.90 | 0.44 | 0.48 | 0.54 | 0.92 | 0.92 | 0.89 | 0.94 |
| Family A: Clean one-to-one | Small | 6.50 | 17.09 | 13.25 | 9.14 | 0.54 | 0.79 | 0.96 | 0.53 | 0.19 | 0.22 | 0.31 | 0.87 | 0.81 | 0.84 | 0.87 | 0.44 | 0.49 | 0.60 | 0.93 | 0.91 | 0.91 | 0.93 |
| Family B: Paired subgroup overlap | Large | 8.30 | 16.72 | 14.04 | 8.10 | 1.00 | 1.11 | 1.08 | 0.86 | 0.19 | 0.21 | 0.31 | 0.68 | 0.66 | 0.67 | 0.72 | 0.44 | 0.47 | 0.60 | 0.87 | 0.87 | 0.87 | 0.89 |
| Family B: Paired subgroup overlap | Medium | 6.15 | 16.76 | 14.15 | 8.53 | 0.94 | 1.20 | 1.34 | 0.95 | 0.19 | 0.20 | 0.26 | 0.71 | 0.67 | 0.68 | 0.71 | 0.44 | 0.47 | 0.56 | 0.89 | 0.87 | 0.87 | 0.88 |
| Family B: Paired subgroup overlap | Small | 4.65 | 16.90 | 13.68 | 10.39 | 1.13 | 1.32 | 1.93 | 1.07 | 0.19 | 0.21 | 0.29 | 0.67 | 0.66 | 0.67 | 0.71 | 0.44 | 0.49 | 0.55 | 0.88 | 0.87 | 0.86 | 0.88 |
| Family C: Shared super-prototype | Large | 6.60 | 16.92 | 14.07 | 8.71 | 1.12 | 1.35 | 1.78 | 0.91 | 0.19 | 0.21 | 0.26 | 0.65 | 0.63 | 0.67 | 0.70 | 0.43 | 0.47 | 0.58 | 0.87 | 0.85 | 0.85 | 0.88 |
| Family C: Shared super-prototype | Medium | 6.45 | 16.88 | 13.62 | 9.15 | 1.00 | 1.24 | 2.07 | 0.85 | 0.19 | 0.22 | 0.29 | 0.70 | 0.66 | 0.69 | 0.74 | 0.44 | 0.49 | 0.58 | 0.88 | 0.87 | 0.85 | 0.89 |
| Family C: Shared super-prototype | Small | 4.70 | 17.24 | 13.12 | 7.44 | 1.00 | 1.39 | 3.23 | 0.90 | 0.19 | 0.22 | 0.34 | 0.71 | 0.66 | 0.68 | 0.76 | 0.44 | 0.50 | 0.64 | 0.89 | 0.87 | 0.83 | 0.90 |
| Family D: Shared prototype with noise | Large | 5.50 | 16.63 | 15.12 | 9.85 | 1.23 | 1.54 | 1.18 | 0.93 | 0.19 | 0.21 | 0.23 | 0.62 | 0.60 | 0.63 | 0.71 | 0.44 | 0.47 | 0.52 | 0.86 | 0.84 | 0.86 | 0.88 |
| Family D: Shared prototype with noise | Medium | 7.15 | 16.74 | 14.35 | 8.94 | 1.17 | 1.67 | 1.51 | 1.06 | 0.19 | 0.22 | 0.28 | 0.65 | 0.59 | 0.61 | 0.68 | 0.44 | 0.48 | 0.57 | 0.87 | 0.84 | 0.85 | 0.88 |
| Family D: Shared prototype with noise | Small | 6.30 | 17.07 | 13.51 | 7.70 | 1.05 | 1.65 | 2.67 | 0.95 | 0.19 | 0.23 | 0.37 | 0.70 | 0.62 | 0.66 | 0.74 | 0.44 | 0.50 | 0.64 | 0.88 | 0.85 | 0.84 | 0.89 |
| Family E: Cross-block leakage | Large | 9.70 | 16.65 | 14.70 | 9.46 | 1.77 | 2.09 | 2.65 | 1.85 | 0.19 | 0.21 | 0.25 | 0.59 | 0.55 | 0.59 | 0.55 | 0.44 | 0.47 | 0.56 | 0.83 | 0.81 | 0.82 | 0.82 |
| Family E: Cross-block leakage | Medium | 8.55 | 16.88 | 13.88 | 7.59 | 1.76 | 2.16 | 3.28 | 1.62 | 0.19 | 0.22 | 0.37 | 0.58 | 0.56 | 0.57 | 0.61 | 0.44 | 0.49 | 0.63 | 0.83 | 0.82 | 0.80 | 0.84 |
| Family E: Cross-block leakage | Small | 7.15 | 16.87 | 13.15 | 8.36 | 1.52 | 2.12 | 4.25 | 1.61 | 0.19 | 0.22 | 0.30 | 0.63 | 0.59 | 0.62 | 0.64 | 0.44 | 0.50 | 0.61 | 0.86 | 0.83 | 0.79 | 0.85 |
