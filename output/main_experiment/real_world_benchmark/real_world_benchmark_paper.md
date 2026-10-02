# Real-World Rectangular Benchmark (Paper Table)

Best $\alpha$ maximizes MWB-AUC on the submitted grid; exact ties use the smaller value.

| Dataset | Shape | Best $\alpha$ | $10^2 \times$ 2-SUM $\downarrow$ O | $10^2 \times$ 2-SUM $\downarrow$ M | $10^2 \times$ 2-SUM $\downarrow$ HC | $10^2 \times$ 2-SUM $\downarrow$ OW | $10^2 \times$ 2-SUM $\downarrow$ CA | $10^2 \times$ 2-SUM $\downarrow$ MD | $10^2 \times$ 2-SUM $\downarrow$ TW | Band@10\% $\uparrow$ O | Band@10\% $\uparrow$ M | Band@10\% $\uparrow$ HC | Band@10\% $\uparrow$ OW | Band@10\% $\uparrow$ CA | Band@10\% $\uparrow$ MD | Band@10\% $\uparrow$ TW | MWB-AUC $\uparrow$ O | MWB-AUC $\uparrow$ M | MWB-AUC $\uparrow$ HC | MWB-AUC $\uparrow$ OW | MWB-AUC $\uparrow$ CA | MWB-AUC $\uparrow$ MD | MWB-AUC $\uparrow$ TW |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| SIC -> NAICS | 1004x1176 | 1 | 2.61 | 10.78 | 2.11 | 0.22 | 2.88 | 0.76 | 0.16 | 0.71 | 0.36 | 0.56 | 0.98 | 0.56 | 0.91 | 0.98 | 0.82 | 0.56 | 0.78 | 0.94 | 0.73 | 0.92 | 0.96 |
| CIP -> SOC | 2143x868 | 6 | 15.90 | 16.04 | 12.18 | 2.81 | 9.77 | 12.37 | 1.95 | 0.19 | 0.27 | 0.20 | 0.38 | 0.30 | 0.16 | 0.56 | 0.44 | 0.50 | 0.53 | 0.74 | 0.57 | 0.40 | 0.81 |
| ACS OCCP x INDP | 24x10 | 4 | 20.67 | 8.87 | 15.94 | 3.39 | 2.23 | 3.05 | 3.61 | 0.14 | 0.29 | 0.15 | 0.43 | 0.77 | 0.77 | 0.34 | 0.52 | 0.73 | 0.53 | 0.92 | 0.95 | 0.91 | 0.92 |
| LODES Home x Work County | 95x95 | 4 | 10.98 | 3.98 | 2.78 | 1.89 | 5.42 | 4.54 | 1.90 | 0.60 | 0.72 | 0.69 | 0.78 | 0.79 | 0.77 | 0.78 | 0.69 | 0.82 | 0.82 | 0.87 | 0.86 | 0.83 | 0.87 |
| 20 Newsgroups Term x Document | 60x48 | 12 | 6.98 | 0.99 | 2.74 | 0.32 | 0.56 | 0.47 | 0.42 | 0.27 | 0.91 | 0.06 | 0.94 | 0.93 | 0.93 | 0.94 | 0.68 | 0.93 | 0.72 | 0.95 | 0.97 | 0.95 | 0.95 |
| MBTA GTFS Route x Station | 8x125 | 2 | 17.96 | 9.39 | 3.33 | 1.25 | 3.09 | 2.19 | 1.25 | 0.23 | 0.25 | 0.58 | 0.87 | 0.52 | 0.61 | 0.88 | 0.04 | 0.02 | 0.13 | 0.20 | 0.19 | 0.26 | 0.20 |
| OpenAlex Author x Topic | 60x8 | 1 | 19.72 | 12.88 | 11.79 | 0.44 | 2.53 | 3.36 | 0.45 | 0.08 | 0.30 | 0.19 | 0.85 | 0.55 | 0.19 | 0.85 | 0.55 | 0.67 | 0.68 | 1.00 | 0.94 | 0.99 | 1.00 |
