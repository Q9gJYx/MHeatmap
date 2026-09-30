# Component-wise ordering audit (Section 3.3)

- `single` = released procedure: one Fiedler vector over the whole active support.
- `comp` = Section 3.3 procedure: order each component, arrange by decreasing mass.
- `delta` = comp - single. Negative is an improvement for R2S only.

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
| SIC -> NAICS | One-walk | - | 0.664 | 0.095 | -0.569 |
| SIC -> NAICS | Two-walk | 1 | 0.180 | 0.153 | -0.027 |
| SIC -> NAICS | Two-walk | 2 | 0.343 | 0.139 | -0.204 |
| SIC -> NAICS | Two-walk | 4 | 0.212 | 0.122 | -0.089 |
| SIC -> NAICS | Two-walk | 6 | 0.215 | 0.119 | -0.097 |
| SIC -> NAICS | Two-walk | 8 | 0.215 | 0.115 | -0.100 |
| SIC -> NAICS | Two-walk | 12 | 0.216 | 0.107 | -0.110 |
| CIP -> SOC | One-walk | - | 2.518 | 2.806 | +0.289 |
| CIP -> SOC | Two-walk | 1 | 2.145 | 4.276 | +2.131 |
| CIP -> SOC | Two-walk | 2 | 2.011 | 2.038 | +0.027 |
| CIP -> SOC | Two-walk | 4 | 1.687 | 1.973 | +0.285 |
| CIP -> SOC | Two-walk | 6 | 1.868 | 1.955 | +0.087 |
| CIP -> SOC | Two-walk | 8 | 1.678 | 1.965 | +0.288 |
| CIP -> SOC | Two-walk | 12 | 2.029 | 3.448 | +1.420 |
| ACS OCCP x INDP | One-walk | - | 3.385 | 3.385 | +0.000 |
| ACS OCCP x INDP | Two-walk | 1 | 4.188 | 4.188 | +0.000 |
| ACS OCCP x INDP | Two-walk | 2 | 4.353 | 4.353 | +0.000 |
| ACS OCCP x INDP | Two-walk | 4 | 3.606 | 3.606 | +0.000 |
| ACS OCCP x INDP | Two-walk | 6 | 3.606 | 3.606 | +0.000 |
| ACS OCCP x INDP | Two-walk | 8 | 3.606 | 3.606 | +0.000 |
| ACS OCCP x INDP | Two-walk | 12 | 3.606 | 3.606 | +0.000 |
| LODES Home x Work | One-walk | - | 1.895 | 1.895 | +0.000 |
| LODES Home x Work | Two-walk | 1 | 2.447 | 2.447 | +0.000 |
| LODES Home x Work | Two-walk | 2 | 1.990 | 1.990 | +0.000 |
| LODES Home x Work | Two-walk | 4 | 1.899 | 1.899 | +0.000 |
| LODES Home x Work | Two-walk | 6 | 1.899 | 1.899 | +0.000 |
| LODES Home x Work | Two-walk | 8 | 1.892 | 1.892 | +0.000 |
| LODES Home x Work | Two-walk | 12 | 1.881 | 1.881 | +0.000 |
| 20 Newsgroups | One-walk | - | 0.426 | 0.426 | +0.000 |
| 20 Newsgroups | Two-walk | 1 | 0.864 | 0.864 | +0.000 |
| 20 Newsgroups | Two-walk | 2 | 0.630 | 0.630 | +0.000 |
| 20 Newsgroups | Two-walk | 4 | 0.429 | 0.429 | +0.000 |
| 20 Newsgroups | Two-walk | 6 | 0.427 | 0.427 | +0.000 |
| 20 Newsgroups | Two-walk | 8 | 0.579 | 0.579 | +0.000 |
| 20 Newsgroups | Two-walk | 12 | 0.571 | 0.571 | +0.000 |
| MBTA Route x Station | One-walk | - | 1.246 | 1.246 | +0.000 |
| MBTA Route x Station | Two-walk | 1 | 1.250 | 1.250 | +0.000 |
| MBTA Route x Station | Two-walk | 2 | 1.246 | 1.246 | +0.000 |
| MBTA Route x Station | Two-walk | 4 | 1.246 | 1.246 | +0.000 |
| MBTA Route x Station | Two-walk | 6 | 1.246 | 1.246 | +0.000 |
| MBTA Route x Station | Two-walk | 8 | 1.246 | 1.246 | +0.000 |
| MBTA Route x Station | Two-walk | 12 | 1.246 | 1.246 | +0.000 |
| OpenAlex Author x Topic | One-walk | - | 0.744 | 0.430 | -0.315 |
| OpenAlex Author x Topic | Two-walk | 1 | 1.062 | 0.437 | -0.625 |
| OpenAlex Author x Topic | Two-walk | 2 | 1.079 | 0.437 | -0.642 |
| OpenAlex Author x Topic | Two-walk | 4 | 0.987 | 0.407 | -0.580 |
| OpenAlex Author x Topic | Two-walk | 6 | 1.013 | 0.411 | -0.601 |
| OpenAlex Author x Topic | Two-walk | 8 | 0.750 | 0.442 | -0.307 |
| OpenAlex Author x Topic | Two-walk | 12 | 1.945 | 0.433 | -1.512 |

## Band@10% (higher is better)

**Released integer-window definition (matches Table 1):**

*Higher is better.*

| Dataset | Method | alpha | single | comp | delta |
|---|---|---:|---:|---:|---:|
| SIC -> NAICS | One-walk | - | 0.976 | 0.990 | +0.013 |
| SIC -> NAICS | Two-walk | 1 | 0.961 | 0.981 | +0.019 |
| SIC -> NAICS | Two-walk | 2 | 0.890 | 0.975 | +0.085 |
| SIC -> NAICS | Two-walk | 4 | 0.954 | 0.977 | +0.023 |
| SIC -> NAICS | Two-walk | 6 | 0.956 | 0.984 | +0.029 |
| SIC -> NAICS | Two-walk | 8 | 0.954 | 0.986 | +0.032 |
| SIC -> NAICS | Two-walk | 12 | 0.957 | 0.988 | +0.031 |
| CIP -> SOC | One-walk | - | 0.435 | 0.379 | -0.057 |
| CIP -> SOC | Two-walk | 1 | 0.463 | 0.358 | -0.105 |
| CIP -> SOC | Two-walk | 2 | 0.493 | 0.480 | -0.013 |
| CIP -> SOC | Two-walk | 4 | 0.578 | 0.525 | -0.053 |
| CIP -> SOC | Two-walk | 6 | 0.567 | 0.558 | -0.008 |
| CIP -> SOC | Two-walk | 8 | 0.619 | 0.564 | -0.055 |
| CIP -> SOC | Two-walk | 12 | 0.578 | 0.366 | -0.212 |
| ACS OCCP x INDP | One-walk | - | 0.875 | 0.875 | +0.000 |
| ACS OCCP x INDP | Two-walk | 1 | 0.753 | 0.753 | +0.000 |
| ACS OCCP x INDP | Two-walk | 2 | 0.783 | 0.783 | +0.000 |
| ACS OCCP x INDP | Two-walk | 4 | 0.881 | 0.881 | +0.000 |
| ACS OCCP x INDP | Two-walk | 6 | 0.881 | 0.881 | +0.000 |
| ACS OCCP x INDP | Two-walk | 8 | 0.875 | 0.875 | +0.000 |
| ACS OCCP x INDP | Two-walk | 12 | 0.881 | 0.881 | +0.000 |
| LODES Home x Work | One-walk | - | 0.775 | 0.775 | +0.000 |
| LODES Home x Work | Two-walk | 1 | 0.764 | 0.764 | +0.000 |
| LODES Home x Work | Two-walk | 2 | 0.777 | 0.777 | +0.000 |
| LODES Home x Work | Two-walk | 4 | 0.776 | 0.776 | +0.000 |
| LODES Home x Work | Two-walk | 6 | 0.777 | 0.777 | +0.000 |
| LODES Home x Work | Two-walk | 8 | 0.774 | 0.774 | +0.000 |
| LODES Home x Work | Two-walk | 12 | 0.773 | 0.773 | +0.000 |
| 20 Newsgroups | One-walk | - | 0.944 | 0.944 | +0.000 |
| 20 Newsgroups | Two-walk | 1 | 0.789 | 0.789 | +0.000 |
| 20 Newsgroups | Two-walk | 2 | 0.947 | 0.947 | +0.000 |
| 20 Newsgroups | Two-walk | 4 | 0.944 | 0.944 | +0.000 |
| 20 Newsgroups | Two-walk | 6 | 0.946 | 0.946 | +0.000 |
| 20 Newsgroups | Two-walk | 8 | 0.943 | 0.943 | +0.000 |
| 20 Newsgroups | Two-walk | 12 | 0.944 | 0.944 | +0.000 |
| MBTA Route x Station | One-walk | - | 0.101 | 0.101 | +0.000 |
| MBTA Route x Station | Two-walk | 1 | 0.110 | 0.110 | +0.000 |
| MBTA Route x Station | Two-walk | 2 | 0.101 | 0.101 | +0.000 |
| MBTA Route x Station | Two-walk | 4 | 0.101 | 0.101 | +0.000 |
| MBTA Route x Station | Two-walk | 6 | 0.101 | 0.101 | +0.000 |
| MBTA Route x Station | Two-walk | 8 | 0.101 | 0.101 | +0.000 |
| MBTA Route x Station | Two-walk | 12 | 0.101 | 0.101 | +0.000 |
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
| SIC -> NAICS | One-walk | - | 0.978 | 0.995 | +0.017 |
| SIC -> NAICS | Two-walk | 1 | 0.981 | 0.985 | +0.005 |
| SIC -> NAICS | Two-walk | 2 | 0.924 | 0.987 | +0.063 |
| SIC -> NAICS | Two-walk | 4 | 0.975 | 0.988 | +0.013 |
| SIC -> NAICS | Two-walk | 6 | 0.977 | 0.989 | +0.012 |
| SIC -> NAICS | Two-walk | 8 | 0.974 | 0.990 | +0.017 |
| SIC -> NAICS | Two-walk | 12 | 0.975 | 0.992 | +0.016 |
| CIP -> SOC | One-walk | - | 0.433 | 0.377 | -0.057 |
| CIP -> SOC | Two-walk | 1 | 0.462 | 0.357 | -0.105 |
| CIP -> SOC | Two-walk | 2 | 0.492 | 0.479 | -0.013 |
| CIP -> SOC | Two-walk | 4 | 0.576 | 0.523 | -0.053 |
| CIP -> SOC | Two-walk | 6 | 0.566 | 0.558 | -0.009 |
| CIP -> SOC | Two-walk | 8 | 0.618 | 0.563 | -0.055 |
| CIP -> SOC | Two-walk | 12 | 0.576 | 0.365 | -0.211 |
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
| 20 Newsgroups | One-walk | - | 0.941 | 0.941 | +0.000 |
| 20 Newsgroups | Two-walk | 1 | 0.785 | 0.785 | +0.000 |
| 20 Newsgroups | Two-walk | 2 | 0.946 | 0.946 | +0.000 |
| 20 Newsgroups | Two-walk | 4 | 0.942 | 0.942 | +0.000 |
| 20 Newsgroups | Two-walk | 6 | 0.944 | 0.944 | +0.000 |
| 20 Newsgroups | Two-walk | 8 | 0.938 | 0.938 | +0.000 |
| 20 Newsgroups | Two-walk | 12 | 0.941 | 0.941 | +0.000 |
| MBTA Route x Station | One-walk | - | 0.874 | 0.874 | +0.000 |
| MBTA Route x Station | Two-walk | 1 | 0.875 | 0.875 | +0.000 |
| MBTA Route x Station | Two-walk | 2 | 0.875 | 0.875 | +0.000 |
| MBTA Route x Station | Two-walk | 4 | 0.875 | 0.875 | +0.000 |
| MBTA Route x Station | Two-walk | 6 | 0.874 | 0.874 | +0.000 |
| MBTA Route x Station | Two-walk | 8 | 0.874 | 0.874 | +0.000 |
| MBTA Route x Station | Two-walk | 12 | 0.874 | 0.874 | +0.000 |
| OpenAlex Author x Topic | One-walk | - | 0.810 | 0.898 | +0.088 |
| OpenAlex Author x Topic | Two-walk | 1 | 0.661 | 0.874 | +0.213 |
| OpenAlex Author x Topic | Two-walk | 2 | 0.663 | 0.874 | +0.212 |
| OpenAlex Author x Topic | Two-walk | 4 | 0.686 | 0.898 | +0.212 |
| OpenAlex Author x Topic | Two-walk | 6 | 0.697 | 0.895 | +0.199 |
| OpenAlex Author x Topic | Two-walk | 8 | 0.818 | 0.874 | +0.056 |
| OpenAlex Author x Topic | Two-walk | 12 | 0.661 | 0.895 | +0.234 |

## MWB-AUC (higher is better)

**Released band-mass integral (matches Table 1):**

*Higher is better.*

| Dataset | Method | alpha | single | comp | delta |
|---|---|---:|---:|---:|---:|
| SIC -> NAICS | One-walk | - | 0.944 | 0.975 | +0.031 |
| SIC -> NAICS | Two-walk | 1 | 0.970 | 0.963 | -0.006 |
| SIC -> NAICS | Two-walk | 2 | 0.941 | 0.966 | +0.025 |
| SIC -> NAICS | Two-walk | 4 | 0.963 | 0.969 | +0.006 |
| SIC -> NAICS | Two-walk | 6 | 0.962 | 0.970 | +0.009 |
| SIC -> NAICS | Two-walk | 8 | 0.958 | 0.971 | +0.012 |
| SIC -> NAICS | Two-walk | 12 | 0.957 | 0.972 | +0.014 |
| CIP -> SOC | One-walk | - | 0.768 | 0.745 | -0.023 |
| CIP -> SOC | Two-walk | 1 | 0.791 | 0.686 | -0.105 |
| CIP -> SOC | Two-walk | 2 | 0.801 | 0.797 | -0.004 |
| CIP -> SOC | Two-walk | 4 | 0.827 | 0.804 | -0.023 |
| CIP -> SOC | Two-walk | 6 | 0.813 | 0.807 | -0.006 |
| CIP -> SOC | Two-walk | 8 | 0.830 | 0.807 | -0.023 |
| CIP -> SOC | Two-walk | 12 | 0.808 | 0.719 | -0.090 |
| ACS OCCP x INDP | One-walk | - | 0.921 | 0.921 | +0.000 |
| ACS OCCP x INDP | Two-walk | 1 | 0.871 | 0.871 | +0.000 |
| ACS OCCP x INDP | Two-walk | 2 | 0.880 | 0.880 | +0.000 |
| ACS OCCP x INDP | Two-walk | 4 | 0.922 | 0.922 | +0.000 |
| ACS OCCP x INDP | Two-walk | 6 | 0.922 | 0.922 | +0.000 |
| ACS OCCP x INDP | Two-walk | 8 | 0.921 | 0.921 | +0.000 |
| ACS OCCP x INDP | Two-walk | 12 | 0.922 | 0.922 | +0.000 |
| LODES Home x Work | One-walk | - | 0.871 | 0.871 | +0.000 |
| LODES Home x Work | Two-walk | 1 | 0.861 | 0.861 | +0.000 |
| LODES Home x Work | Two-walk | 2 | 0.864 | 0.864 | +0.000 |
| LODES Home x Work | Two-walk | 4 | 0.871 | 0.871 | +0.000 |
| LODES Home x Work | Two-walk | 6 | 0.866 | 0.866 | +0.000 |
| LODES Home x Work | Two-walk | 8 | 0.871 | 0.871 | +0.000 |
| LODES Home x Work | Two-walk | 12 | 0.870 | 0.870 | +0.000 |
| 20 Newsgroups | One-walk | - | 0.945 | 0.945 | +0.000 |
| 20 Newsgroups | Two-walk | 1 | 0.871 | 0.871 | +0.000 |
| 20 Newsgroups | Two-walk | 2 | 0.911 | 0.911 | +0.000 |
| 20 Newsgroups | Two-walk | 4 | 0.945 | 0.945 | +0.000 |
| 20 Newsgroups | Two-walk | 6 | 0.946 | 0.946 | +0.000 |
| 20 Newsgroups | Two-walk | 8 | 0.919 | 0.919 | +0.000 |
| 20 Newsgroups | Two-walk | 12 | 0.920 | 0.920 | +0.000 |
| MBTA Route x Station | One-walk | - | 0.139 | 0.139 | +0.000 |
| MBTA Route x Station | Two-walk | 1 | 0.151 | 0.151 | +0.000 |
| MBTA Route x Station | Two-walk | 2 | 0.139 | 0.139 | +0.000 |
| MBTA Route x Station | Two-walk | 4 | 0.139 | 0.139 | +0.000 |
| MBTA Route x Station | Two-walk | 6 | 0.139 | 0.139 | +0.000 |
| MBTA Route x Station | Two-walk | 8 | 0.139 | 0.139 | +0.000 |
| MBTA Route x Station | Two-walk | 12 | 0.139 | 0.139 | +0.000 |
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
| SIC -> NAICS | One-walk | - | 0.683 | 0.752 | +0.069 |
| SIC -> NAICS | Two-walk | 1 | 0.758 | 0.689 | -0.069 |
| SIC -> NAICS | Two-walk | 2 | 0.665 | 0.705 | +0.040 |
| SIC -> NAICS | Two-walk | 4 | 0.741 | 0.717 | -0.024 |
| SIC -> NAICS | Two-walk | 6 | 0.726 | 0.730 | +0.004 |
| SIC -> NAICS | Two-walk | 8 | 0.701 | 0.737 | +0.036 |
| SIC -> NAICS | Two-walk | 12 | 0.706 | 0.731 | +0.025 |
| CIP -> SOC | One-walk | - | 0.283 | 0.251 | -0.032 |
| CIP -> SOC | Two-walk | 1 | 0.298 | 0.255 | -0.043 |
| CIP -> SOC | Two-walk | 2 | 0.320 | 0.312 | -0.009 |
| CIP -> SOC | Two-walk | 4 | 0.360 | 0.323 | -0.037 |
| CIP -> SOC | Two-walk | 6 | 0.334 | 0.323 | -0.010 |
| CIP -> SOC | Two-walk | 8 | 0.361 | 0.330 | -0.031 |
| CIP -> SOC | Two-walk | 12 | 0.351 | 0.272 | -0.079 |
| ACS OCCP x INDP | One-walk | - | 0.282 | 0.282 | +0.000 |
| ACS OCCP x INDP | Two-walk | 1 | 0.222 | 0.222 | +0.000 |
| ACS OCCP x INDP | Two-walk | 2 | 0.266 | 0.266 | +0.000 |
| ACS OCCP x INDP | Two-walk | 4 | 0.251 | 0.251 | +0.000 |
| ACS OCCP x INDP | Two-walk | 6 | 0.251 | 0.251 | +0.000 |
| ACS OCCP x INDP | Two-walk | 8 | 0.260 | 0.260 | +0.000 |
| ACS OCCP x INDP | Two-walk | 12 | 0.251 | 0.251 | +0.000 |
| LODES Home x Work | One-walk | - | 0.549 | 0.549 | +0.000 |
| LODES Home x Work | Two-walk | 1 | 0.525 | 0.525 | +0.000 |
| LODES Home x Work | Two-walk | 2 | 0.512 | 0.512 | +0.000 |
| LODES Home x Work | Two-walk | 4 | 0.528 | 0.528 | +0.000 |
| LODES Home x Work | Two-walk | 6 | 0.530 | 0.530 | +0.000 |
| LODES Home x Work | Two-walk | 8 | 0.509 | 0.509 | +0.000 |
| LODES Home x Work | Two-walk | 12 | 0.510 | 0.510 | +0.000 |
| 20 Newsgroups | One-walk | - | 0.533 | 0.533 | +0.000 |
| 20 Newsgroups | Two-walk | 1 | 0.378 | 0.378 | +0.000 |
| 20 Newsgroups | Two-walk | 2 | 0.388 | 0.388 | +0.000 |
| 20 Newsgroups | Two-walk | 4 | 0.533 | 0.533 | +0.000 |
| 20 Newsgroups | Two-walk | 6 | 0.532 | 0.532 | +0.000 |
| 20 Newsgroups | Two-walk | 8 | 0.508 | 0.508 | +0.000 |
| 20 Newsgroups | Two-walk | 12 | 0.510 | 0.510 | +0.000 |
| MBTA Route x Station | One-walk | - | 0.565 | 0.565 | +0.000 |
| MBTA Route x Station | Two-walk | 1 | 0.558 | 0.558 | +0.000 |
| MBTA Route x Station | Two-walk | 2 | 0.565 | 0.565 | +0.000 |
| MBTA Route x Station | Two-walk | 4 | 0.565 | 0.565 | +0.000 |
| MBTA Route x Station | Two-walk | 6 | 0.565 | 0.565 | +0.000 |
| MBTA Route x Station | Two-walk | 8 | 0.565 | 0.565 | +0.000 |
| MBTA Route x Station | Two-walk | 12 | 0.565 | 0.565 | +0.000 |
| OpenAlex Author x Topic | One-walk | - | 0.559 | 0.503 | -0.055 |
| OpenAlex Author x Topic | Two-walk | 1 | 0.442 | 0.515 | +0.074 |
| OpenAlex Author x Topic | Two-walk | 2 | 0.433 | 0.516 | +0.083 |
| OpenAlex Author x Topic | Two-walk | 4 | 0.437 | 0.509 | +0.072 |
| OpenAlex Author x Topic | Two-walk | 6 | 0.447 | 0.507 | +0.059 |
| OpenAlex Author x Topic | Two-walk | 8 | 0.554 | 0.514 | -0.040 |
| OpenAlex Author x Topic | Two-walk | 12 | 0.438 | 0.495 | +0.056 |

## Self-check

28 single-component configurations were run through both procedures; these must be identical.

All identical, as required.
