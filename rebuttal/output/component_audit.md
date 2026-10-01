# Component-wise ordering audit (Section 3.3)

- `single` = pre-fix procedure: one Fiedler vector over the whole active
  support, finished with the released `orient_orders_for_diagonal` step.
- `comp` = Section 3.3 procedure, via the benchmark's canonical
  implementation: order each component, arrange by decreasing mass.
- `delta` = comp - single. Negative is an improvement for R2S only.

**A nonzero delta on a single-component dataset is expected.** The two
procedures differ by the component correction *and* by the final
orientation step, which the paper-spec commit removed from One-walk and
Two-walk. Only the four datasets with one component isolate the second
change; on those, any delta is entirely the orientation effect.

**Metric definitions.** `band_released` and `mwb_band` reproduce the released
implementation and therefore Table 1. `band_eq16` (Eq. 16) and `mwb_block`
(Section 4.5 matched block cuts) are the paper-faithful definitions. They
differ from the released ones in general -- see reviewer kFbK W1.

## Component structure

| Dataset | Rank | Rows | Cols | Mass | Mass fraction |
|---|---:|---:|---:|---:|---:|
| SIC -> NAICS | 1 | 516 | 644 | 1527.0 | 70.6% |
| SIC -> NAICS | 2 | 8 | 8 | 26.0 | 1.2% |
| SIC -> NAICS | 3 | 10 | 10 | 20.0 | 0.9% |
| SIC -> NAICS | 4 | 3 | 9 | 13.0 | 0.6% |
| SIC -> NAICS | 5 | 4 | 6 | 13.0 | 0.6% |
| SIC -> NAICS | 6 | 6 | 6 | 11.0 | 0.5% |
| SIC -> NAICS | 7 | 5 | 4 | 11.0 | 0.5% |
| SIC -> NAICS | 8 | 4 | 5 | 9.0 | 0.4% |
| CIP -> SOC | 1 | 1872 | 591 | 5578.0 | 91.5% |
| CIP -> SOC | 2 | 194 | 1 | 194.0 | 3.2% |
| CIP -> SOC | 3 | 1 | 180 | 180.0 | 3.0% |
| CIP -> SOC | 4 | 6 | 20 | 31.0 | 0.5% |
| CIP -> SOC | 5 | 7 | 2 | 14.0 | 0.2% |
| CIP -> SOC | 6 | 5 | 7 | 13.0 | 0.2% |
| CIP -> SOC | 7 | 8 | 3 | 10.0 | 0.2% |
| CIP -> SOC | 8 | 1 | 8 | 8.0 | 0.1% |
| ACS OCCP x INDP | 1 | 24 | 10 | 8.7 | 100.0% |
| LODES Home x Work | 1 | 95 | 95 | 8.5 | 100.0% |
| 20 Newsgroups | 1 | 59 | 48 | 3.3 | 100.0% |
| MBTA Route x Station | 1 | 8 | 125 | 100.2 | 100.0% |
| OpenAlex Author x Topic | 1 | 27 | 4 | 9.2 | 39.6% |
| OpenAlex Author x Topic | 2 | 15 | 1 | 6.6 | 28.6% |
| OpenAlex Author x Topic | 3 | 8 | 1 | 4.4 | 19.0% |
| OpenAlex Author x Topic | 4 | 9 | 1 | 2.7 | 11.8% |
| OpenAlex Author x Topic | 5 | 1 | 1 | 0.2 | 1.0% |

## R2S (x100, lower is better)

*Lower is better.*

| Dataset | Method | alpha | single | comp | delta |
|---|---|---:|---:|---:|---:|
| SIC -> NAICS | One-walk | - | 0.664 | 0.219 | -0.445 |
| SIC -> NAICS | Two-walk | 1 | 0.316 | 0.164 | -0.152 |
| SIC -> NAICS | Two-walk | 2 | 0.343 | 0.186 | -0.157 |
| SIC -> NAICS | Two-walk | 4 | 0.212 | 0.188 | -0.023 |
| SIC -> NAICS | Two-walk | 6 | 0.215 | 0.192 | -0.023 |
| SIC -> NAICS | Two-walk | 8 | 0.215 | 0.197 | -0.017 |
| SIC -> NAICS | Two-walk | 12 | 0.216 | 0.202 | -0.014 |
| CIP -> SOC | One-walk | - | 2.518 | 2.806 | +0.289 |
| CIP -> SOC | Two-walk | 1 | 2.155 | 2.141 | -0.014 |
| CIP -> SOC | Two-walk | 2 | 2.011 | 2.038 | +0.027 |
| CIP -> SOC | Two-walk | 4 | 1.687 | 1.973 | +0.285 |
| CIP -> SOC | Two-walk | 6 | 1.868 | 1.954 | +0.087 |
| CIP -> SOC | Two-walk | 8 | 1.678 | 1.965 | +0.288 |
| CIP -> SOC | Two-walk | 12 | 2.029 | 2.050 | +0.021 |
| ACS OCCP x INDP | One-walk | - | 3.385 | 3.385 | +0.000 |
| ACS OCCP x INDP | Two-walk | 1 | 4.188 | 4.188 | +0.000 |
| ACS OCCP x INDP | Two-walk | 2 | 4.353 | 4.353 | +0.000 |
| ACS OCCP x INDP | Two-walk | 4 | 3.606 | 3.606 | +0.000 |
| ACS OCCP x INDP | Two-walk | 6 | 3.606 | 3.606 | +0.000 |
| ACS OCCP x INDP | Two-walk | 8 | 3.606 | 3.606 | -0.000 |
| ACS OCCP x INDP | Two-walk | 12 | 3.606 | 3.606 | +0.000 |
| LODES Home x Work | One-walk | - | 1.895 | 1.895 | +0.000 |
| LODES Home x Work | Two-walk | 1 | 2.447 | 2.447 | +0.000 |
| LODES Home x Work | Two-walk | 2 | 1.990 | 1.990 | +0.000 |
| LODES Home x Work | Two-walk | 4 | 1.899 | 1.899 | +0.000 |
| LODES Home x Work | Two-walk | 6 | 1.899 | 1.899 | +0.000 |
| LODES Home x Work | Two-walk | 8 | 1.892 | 1.892 | -0.000 |
| LODES Home x Work | Two-walk | 12 | 1.881 | 1.881 | +0.000 |
| 20 Newsgroups | One-walk | - | 0.426 | 0.324 | -0.103 |
| 20 Newsgroups | Two-walk | 1 | 0.864 | 0.631 | -0.233 |
| 20 Newsgroups | Two-walk | 2 | 0.630 | 0.630 | +0.000 |
| 20 Newsgroups | Two-walk | 4 | 0.429 | 0.429 | +0.000 |
| 20 Newsgroups | Two-walk | 6 | 0.427 | 0.427 | +0.000 |
| 20 Newsgroups | Two-walk | 8 | 0.579 | 0.427 | -0.152 |
| 20 Newsgroups | Two-walk | 12 | 0.571 | 0.420 | -0.151 |
| MBTA Route x Station | One-walk | - | 1.246 | 1.246 | +0.000 |
| MBTA Route x Station | Two-walk | 1 | 1.250 | 1.250 | +0.000 |
| MBTA Route x Station | Two-walk | 2 | 1.246 | 1.246 | +0.000 |
| MBTA Route x Station | Two-walk | 4 | 1.246 | 1.246 | +0.000 |
| MBTA Route x Station | Two-walk | 6 | 1.246 | 1.246 | +0.000 |
| MBTA Route x Station | Two-walk | 8 | 1.246 | 1.246 | +0.000 |
| MBTA Route x Station | Two-walk | 12 | 1.246 | 1.246 | +0.000 |
| OpenAlex Author x Topic | One-walk | - | 0.744 | 0.439 | -0.305 |
| OpenAlex Author x Topic | Two-walk | 1 | 0.789 | 0.450 | -0.338 |
| OpenAlex Author x Topic | Two-walk | 2 | 1.079 | 0.453 | -0.626 |
| OpenAlex Author x Topic | Two-walk | 4 | 0.987 | 0.453 | -0.534 |
| OpenAlex Author x Topic | Two-walk | 6 | 1.013 | 0.453 | -0.560 |
| OpenAlex Author x Topic | Two-walk | 8 | 0.750 | 0.449 | -0.300 |
| OpenAlex Author x Topic | Two-walk | 12 | 1.945 | 0.472 | -1.473 |

## Band@10% (higher is better)

**Released integer-window definition (matches Table 1):**

*Higher is better.*

| Dataset | Method | alpha | single | comp | delta |
|---|---|---:|---:|---:|---:|
| SIC -> NAICS | One-walk | - | 0.976 | 0.973 | -0.003 |
| SIC -> NAICS | Two-walk | 1 | 0.906 | 0.963 | +0.058 |
| SIC -> NAICS | Two-walk | 2 | 0.890 | 0.957 | +0.067 |
| SIC -> NAICS | Two-walk | 4 | 0.954 | 0.956 | +0.002 |
| SIC -> NAICS | Two-walk | 6 | 0.956 | 0.957 | +0.002 |
| SIC -> NAICS | Two-walk | 8 | 0.954 | 0.955 | +0.001 |
| SIC -> NAICS | Two-walk | 12 | 0.957 | 0.959 | +0.001 |
| CIP -> SOC | One-walk | - | 0.435 | 0.377 | -0.059 |
| CIP -> SOC | Two-walk | 1 | 0.459 | 0.457 | -0.002 |
| CIP -> SOC | Two-walk | 2 | 0.493 | 0.480 | -0.013 |
| CIP -> SOC | Two-walk | 4 | 0.578 | 0.525 | -0.053 |
| CIP -> SOC | Two-walk | 6 | 0.567 | 0.558 | -0.008 |
| CIP -> SOC | Two-walk | 8 | 0.619 | 0.564 | -0.055 |
| CIP -> SOC | Two-walk | 12 | 0.578 | 0.564 | -0.013 |
| ACS OCCP x INDP | One-walk | - | 0.875 | 0.884 | +0.009 |
| ACS OCCP x INDP | Two-walk | 1 | 0.753 | 0.800 | +0.047 |
| ACS OCCP x INDP | Two-walk | 2 | 0.783 | 0.783 | +0.000 |
| ACS OCCP x INDP | Two-walk | 4 | 0.881 | 0.881 | +0.000 |
| ACS OCCP x INDP | Two-walk | 6 | 0.881 | 0.881 | +0.000 |
| ACS OCCP x INDP | Two-walk | 8 | 0.875 | 0.881 | +0.005 |
| ACS OCCP x INDP | Two-walk | 12 | 0.881 | 0.881 | +0.000 |
| LODES Home x Work | One-walk | - | 0.775 | 0.775 | +0.000 |
| LODES Home x Work | Two-walk | 1 | 0.764 | 0.764 | +0.000 |
| LODES Home x Work | Two-walk | 2 | 0.777 | 0.777 | +0.000 |
| LODES Home x Work | Two-walk | 4 | 0.776 | 0.776 | +0.000 |
| LODES Home x Work | Two-walk | 6 | 0.777 | 0.777 | +0.000 |
| LODES Home x Work | Two-walk | 8 | 0.774 | 0.774 | +0.000 |
| LODES Home x Work | Two-walk | 12 | 0.773 | 0.773 | +0.000 |
| 20 Newsgroups | One-walk | - | 0.944 | 0.945 | +0.001 |
| 20 Newsgroups | Two-walk | 1 | 0.789 | 0.944 | +0.155 |
| 20 Newsgroups | Two-walk | 2 | 0.947 | 0.947 | -0.001 |
| 20 Newsgroups | Two-walk | 4 | 0.944 | 0.944 | -0.000 |
| 20 Newsgroups | Two-walk | 6 | 0.946 | 0.945 | -0.000 |
| 20 Newsgroups | Two-walk | 8 | 0.943 | 0.945 | +0.002 |
| 20 Newsgroups | Two-walk | 12 | 0.944 | 0.945 | +0.001 |
| MBTA Route x Station | One-walk | - | 0.101 | 0.146 | +0.045 |
| MBTA Route x Station | Two-walk | 1 | 0.110 | 0.146 | +0.036 |
| MBTA Route x Station | Two-walk | 2 | 0.101 | 0.146 | +0.045 |
| MBTA Route x Station | Two-walk | 4 | 0.101 | 0.146 | +0.045 |
| MBTA Route x Station | Two-walk | 6 | 0.101 | 0.146 | +0.045 |
| MBTA Route x Station | Two-walk | 8 | 0.101 | 0.146 | +0.045 |
| MBTA Route x Station | Two-walk | 12 | 0.101 | 0.146 | +0.045 |
| OpenAlex Author x Topic | One-walk | - | 0.971 | 1.000 | +0.029 |
| OpenAlex Author x Topic | Two-walk | 1 | 0.996 | 1.000 | +0.004 |
| OpenAlex Author x Topic | Two-walk | 2 | 0.996 | 1.000 | +0.004 |
| OpenAlex Author x Topic | Two-walk | 4 | 0.996 | 1.000 | +0.004 |
| OpenAlex Author x Topic | Two-walk | 6 | 0.996 | 1.000 | +0.004 |
| OpenAlex Author x Topic | Two-walk | 8 | 0.996 | 1.000 | +0.004 |
| OpenAlex Author x Topic | Two-walk | 12 | 0.986 | 1.000 | +0.014 |

**Eq. (16) definition (paper-faithful):**

*Higher is better.*

| Dataset | Method | alpha | single | comp | delta |
|---|---|---:|---:|---:|---:|
| SIC -> NAICS | One-walk | - | 0.978 | 0.982 | +0.004 |
| SIC -> NAICS | Two-walk | 1 | 0.940 | 0.983 | +0.043 |
| SIC -> NAICS | Two-walk | 2 | 0.924 | 0.981 | +0.056 |
| SIC -> NAICS | Two-walk | 4 | 0.975 | 0.976 | +0.002 |
| SIC -> NAICS | Two-walk | 6 | 0.977 | 0.979 | +0.002 |
| SIC -> NAICS | Two-walk | 8 | 0.974 | 0.975 | +0.001 |
| SIC -> NAICS | Two-walk | 12 | 0.975 | 0.977 | +0.001 |
| CIP -> SOC | One-walk | - | 0.433 | 0.377 | -0.057 |
| CIP -> SOC | Two-walk | 1 | 0.459 | 0.455 | -0.003 |
| CIP -> SOC | Two-walk | 2 | 0.492 | 0.479 | -0.013 |
| CIP -> SOC | Two-walk | 4 | 0.576 | 0.523 | -0.053 |
| CIP -> SOC | Two-walk | 6 | 0.566 | 0.558 | -0.009 |
| CIP -> SOC | Two-walk | 8 | 0.618 | 0.563 | -0.055 |
| CIP -> SOC | Two-walk | 12 | 0.576 | 0.564 | -0.013 |
| ACS OCCP x INDP | One-walk | - | 0.431 | 0.431 | +0.000 |
| ACS OCCP x INDP | Two-walk | 1 | 0.394 | 0.394 | +0.000 |
| ACS OCCP x INDP | Two-walk | 2 | 0.304 | 0.304 | +0.000 |
| ACS OCCP x INDP | Two-walk | 4 | 0.339 | 0.339 | +0.000 |
| ACS OCCP x INDP | Two-walk | 6 | 0.339 | 0.339 | +0.000 |
| ACS OCCP x INDP | Two-walk | 8 | 0.339 | 0.339 | +0.000 |
| ACS OCCP x INDP | Two-walk | 12 | 0.339 | 0.339 | +0.000 |
| LODES Home x Work | One-walk | - | 0.775 | 0.775 | +0.000 |
| LODES Home x Work | Two-walk | 1 | 0.764 | 0.764 | +0.000 |
| LODES Home x Work | Two-walk | 2 | 0.777 | 0.777 | +0.000 |
| LODES Home x Work | Two-walk | 4 | 0.776 | 0.776 | +0.000 |
| LODES Home x Work | Two-walk | 6 | 0.777 | 0.777 | +0.000 |
| LODES Home x Work | Two-walk | 8 | 0.774 | 0.774 | +0.000 |
| LODES Home x Work | Two-walk | 12 | 0.773 | 0.773 | +0.000 |
| 20 Newsgroups | One-walk | - | 0.941 | 0.943 | +0.002 |
| 20 Newsgroups | Two-walk | 1 | 0.785 | 0.943 | +0.158 |
| 20 Newsgroups | Two-walk | 2 | 0.946 | 0.946 | +0.000 |
| 20 Newsgroups | Two-walk | 4 | 0.942 | 0.942 | +0.000 |
| 20 Newsgroups | Two-walk | 6 | 0.944 | 0.944 | +0.000 |
| 20 Newsgroups | Two-walk | 8 | 0.938 | 0.944 | +0.006 |
| 20 Newsgroups | Two-walk | 12 | 0.941 | 0.943 | +0.002 |
| MBTA Route x Station | One-walk | - | 0.874 | 0.874 | +0.000 |
| MBTA Route x Station | Two-walk | 1 | 0.875 | 0.875 | +0.000 |
| MBTA Route x Station | Two-walk | 2 | 0.875 | 0.875 | +0.000 |
| MBTA Route x Station | Two-walk | 4 | 0.875 | 0.875 | +0.000 |
| MBTA Route x Station | Two-walk | 6 | 0.874 | 0.874 | +0.000 |
| MBTA Route x Station | Two-walk | 8 | 0.874 | 0.874 | +0.000 |
| MBTA Route x Station | Two-walk | 12 | 0.874 | 0.874 | +0.000 |
| OpenAlex Author x Topic | One-walk | - | 0.810 | 0.847 | +0.037 |
| OpenAlex Author x Topic | Two-walk | 1 | 0.777 | 0.847 | +0.070 |
| OpenAlex Author x Topic | Two-walk | 2 | 0.663 | 0.847 | +0.184 |
| OpenAlex Author x Topic | Two-walk | 4 | 0.686 | 0.847 | +0.161 |
| OpenAlex Author x Topic | Two-walk | 6 | 0.697 | 0.847 | +0.151 |
| OpenAlex Author x Topic | Two-walk | 8 | 0.818 | 0.847 | +0.029 |
| OpenAlex Author x Topic | Two-walk | 12 | 0.661 | 0.847 | +0.186 |

## MWB-AUC (higher is better)

**Released band-mass integral (matches Table 1):**

*Higher is better.*

| Dataset | Method | alpha | single | comp | delta |
|---|---|---:|---:|---:|---:|
| SIC -> NAICS | One-walk | - | 0.944 | 0.942 | -0.002 |
| SIC -> NAICS | Two-walk | 1 | 0.947 | 0.961 | +0.014 |
| SIC -> NAICS | Two-walk | 2 | 0.941 | 0.955 | +0.013 |
| SIC -> NAICS | Two-walk | 4 | 0.963 | 0.954 | -0.009 |
| SIC -> NAICS | Two-walk | 6 | 0.962 | 0.952 | -0.010 |
| SIC -> NAICS | Two-walk | 8 | 0.958 | 0.950 | -0.008 |
| SIC -> NAICS | Two-walk | 12 | 0.957 | 0.947 | -0.010 |
| CIP -> SOC | One-walk | - | 0.768 | 0.744 | -0.024 |
| CIP -> SOC | Two-walk | 1 | 0.790 | 0.789 | -0.001 |
| CIP -> SOC | Two-walk | 2 | 0.801 | 0.797 | -0.004 |
| CIP -> SOC | Two-walk | 4 | 0.827 | 0.804 | -0.023 |
| CIP -> SOC | Two-walk | 6 | 0.813 | 0.807 | -0.006 |
| CIP -> SOC | Two-walk | 8 | 0.830 | 0.807 | -0.023 |
| CIP -> SOC | Two-walk | 12 | 0.808 | 0.804 | -0.004 |
| ACS OCCP x INDP | One-walk | - | 0.921 | 0.924 | +0.002 |
| ACS OCCP x INDP | Two-walk | 1 | 0.871 | 0.885 | +0.015 |
| ACS OCCP x INDP | Two-walk | 2 | 0.880 | 0.880 | +0.000 |
| ACS OCCP x INDP | Two-walk | 4 | 0.922 | 0.922 | +0.000 |
| ACS OCCP x INDP | Two-walk | 6 | 0.922 | 0.922 | +0.000 |
| ACS OCCP x INDP | Two-walk | 8 | 0.921 | 0.922 | +0.001 |
| ACS OCCP x INDP | Two-walk | 12 | 0.922 | 0.922 | +0.000 |
| LODES Home x Work | One-walk | - | 0.871 | 0.871 | +0.000 |
| LODES Home x Work | Two-walk | 1 | 0.861 | 0.861 | +0.000 |
| LODES Home x Work | Two-walk | 2 | 0.864 | 0.864 | +0.000 |
| LODES Home x Work | Two-walk | 4 | 0.871 | 0.871 | +0.000 |
| LODES Home x Work | Two-walk | 6 | 0.866 | 0.866 | +0.000 |
| LODES Home x Work | Two-walk | 8 | 0.871 | 0.871 | +0.000 |
| LODES Home x Work | Two-walk | 12 | 0.870 | 0.870 | +0.000 |
| 20 Newsgroups | One-walk | - | 0.945 | 0.946 | +0.001 |
| 20 Newsgroups | Two-walk | 1 | 0.871 | 0.911 | +0.040 |
| 20 Newsgroups | Two-walk | 2 | 0.911 | 0.911 | +0.000 |
| 20 Newsgroups | Two-walk | 4 | 0.945 | 0.946 | +0.000 |
| 20 Newsgroups | Two-walk | 6 | 0.946 | 0.946 | +0.000 |
| 20 Newsgroups | Two-walk | 8 | 0.919 | 0.946 | +0.027 |
| 20 Newsgroups | Two-walk | 12 | 0.920 | 0.946 | +0.026 |
| MBTA Route x Station | One-walk | - | 0.139 | 0.204 | +0.065 |
| MBTA Route x Station | Two-walk | 1 | 0.151 | 0.204 | +0.054 |
| MBTA Route x Station | Two-walk | 2 | 0.139 | 0.204 | +0.065 |
| MBTA Route x Station | Two-walk | 4 | 0.139 | 0.204 | +0.065 |
| MBTA Route x Station | Two-walk | 6 | 0.139 | 0.204 | +0.065 |
| MBTA Route x Station | Two-walk | 8 | 0.139 | 0.204 | +0.065 |
| MBTA Route x Station | Two-walk | 12 | 0.139 | 0.204 | +0.065 |
| OpenAlex Author x Topic | One-walk | - | 0.984 | 1.000 | +0.016 |
| OpenAlex Author x Topic | Two-walk | 1 | 0.996 | 1.000 | +0.004 |
| OpenAlex Author x Topic | Two-walk | 2 | 0.996 | 1.000 | +0.004 |
| OpenAlex Author x Topic | Two-walk | 4 | 0.996 | 1.000 | +0.004 |
| OpenAlex Author x Topic | Two-walk | 6 | 0.996 | 1.000 | +0.004 |
| OpenAlex Author x Topic | Two-walk | 8 | 0.996 | 1.000 | +0.004 |
| OpenAlex Author x Topic | Two-walk | 12 | 0.988 | 1.000 | +0.012 |

**Section 4.5 matched block cuts (paper-faithful):**

*Higher is better.*

| Dataset | Method | alpha | single | comp | delta |
|---|---|---:|---:|---:|---:|
| SIC -> NAICS | One-walk | - | 0.683 | 0.579 | -0.104 |
| SIC -> NAICS | Two-walk | 1 | 0.687 | 0.684 | -0.004 |
| SIC -> NAICS | Two-walk | 2 | 0.665 | 0.651 | -0.014 |
| SIC -> NAICS | Two-walk | 4 | 0.741 | 0.663 | -0.079 |
| SIC -> NAICS | Two-walk | 6 | 0.726 | 0.644 | -0.082 |
| SIC -> NAICS | Two-walk | 8 | 0.701 | 0.630 | -0.071 |
| SIC -> NAICS | Two-walk | 12 | 0.706 | 0.616 | -0.090 |
| CIP -> SOC | One-walk | - | 0.283 | 0.251 | -0.032 |
| CIP -> SOC | Two-walk | 1 | 0.297 | 0.291 | -0.006 |
| CIP -> SOC | Two-walk | 2 | 0.320 | 0.312 | -0.009 |
| CIP -> SOC | Two-walk | 4 | 0.360 | 0.323 | -0.037 |
| CIP -> SOC | Two-walk | 6 | 0.334 | 0.323 | -0.010 |
| CIP -> SOC | Two-walk | 8 | 0.361 | 0.330 | -0.031 |
| CIP -> SOC | Two-walk | 12 | 0.351 | 0.341 | -0.010 |
| ACS OCCP x INDP | One-walk | - | 0.282 | 0.283 | +0.001 |
| ACS OCCP x INDP | Two-walk | 1 | 0.222 | 0.184 | -0.037 |
| ACS OCCP x INDP | Two-walk | 2 | 0.266 | 0.266 | +0.000 |
| ACS OCCP x INDP | Two-walk | 4 | 0.251 | 0.251 | +0.000 |
| ACS OCCP x INDP | Two-walk | 6 | 0.251 | 0.251 | +0.000 |
| ACS OCCP x INDP | Two-walk | 8 | 0.260 | 0.251 | -0.009 |
| ACS OCCP x INDP | Two-walk | 12 | 0.251 | 0.251 | +0.000 |
| LODES Home x Work | One-walk | - | 0.549 | 0.546 | -0.003 |
| LODES Home x Work | Two-walk | 1 | 0.524 | 0.525 | +0.001 |
| LODES Home x Work | Two-walk | 2 | 0.512 | 0.511 | -0.001 |
| LODES Home x Work | Two-walk | 4 | 0.528 | 0.528 | +0.000 |
| LODES Home x Work | Two-walk | 6 | 0.530 | 0.530 | +0.000 |
| LODES Home x Work | Two-walk | 8 | 0.509 | 0.512 | +0.003 |
| LODES Home x Work | Two-walk | 12 | 0.510 | 0.510 | +0.000 |
| 20 Newsgroups | One-walk | - | 0.533 | 0.678 | +0.144 |
| 20 Newsgroups | Two-walk | 1 | 0.378 | 0.389 | +0.011 |
| 20 Newsgroups | Two-walk | 2 | 0.388 | 0.388 | +0.000 |
| 20 Newsgroups | Two-walk | 4 | 0.533 | 0.533 | +0.000 |
| 20 Newsgroups | Two-walk | 6 | 0.532 | 0.532 | +0.000 |
| 20 Newsgroups | Two-walk | 8 | 0.508 | 0.532 | +0.024 |
| 20 Newsgroups | Two-walk | 12 | 0.510 | 0.533 | +0.023 |
| MBTA Route x Station | One-walk | - | 0.565 | 0.551 | -0.014 |
| MBTA Route x Station | Two-walk | 1 | 0.558 | 0.552 | -0.006 |
| MBTA Route x Station | Two-walk | 2 | 0.565 | 0.551 | -0.014 |
| MBTA Route x Station | Two-walk | 4 | 0.565 | 0.551 | -0.014 |
| MBTA Route x Station | Two-walk | 6 | 0.565 | 0.551 | -0.014 |
| MBTA Route x Station | Two-walk | 8 | 0.565 | 0.551 | -0.014 |
| MBTA Route x Station | Two-walk | 12 | 0.565 | 0.551 | -0.014 |
| OpenAlex Author x Topic | One-walk | - | 0.559 | 0.560 | +0.001 |
| OpenAlex Author x Topic | Two-walk | 1 | 0.506 | 0.543 | +0.037 |
| OpenAlex Author x Topic | Two-walk | 2 | 0.433 | 0.541 | +0.107 |
| OpenAlex Author x Topic | Two-walk | 4 | 0.437 | 0.541 | +0.103 |
| OpenAlex Author x Topic | Two-walk | 6 | 0.447 | 0.541 | +0.094 |
| OpenAlex Author x Topic | Two-walk | 8 | 0.554 | 0.545 | -0.009 |
| OpenAlex Author x Topic | Two-walk | 12 | 0.438 | 0.529 | +0.091 |

## Self-check

4 single-component datasets were checked. With one component there is nothing to assemble, so the component path must reduce to a single solve over exactly the full active support and return a valid permutation.

This does **not** assert that the two procedures agree on those datasets, because they no longer do: the paper-spec commit also dropped the 2-SUM orientation step and switched `argsort` to a stable sort. See the module docstring.

Component machinery is sound on all of them.
