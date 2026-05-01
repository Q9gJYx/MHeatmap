# Real-world benchmark data sources

This directory records the public data sources used by
`main_experiment/real_world_benchmark`. Raw data files are not committed to the
repository. The cleaned matrix CSVs used by the paper table are versioned under
`output/main_experiment/real_world_benchmark/processed_matrices/`.

Sources:

- 1987 SIC to 2002 NAICS concordance, downloaded from the U.S. Census Bureau:
  `https://www2.census.gov/library/reference/naics/technical-documentation/concordance/1987_sic_to_2002_naics.xls`
- CIP 2020 to SOC 2018 crosswalk, downloaded from NCES:
  `https://nces.ed.gov/ipeds/cipcode/Files/CIP2020_SOC2018_Crosswalk.xlsx`
- ACS 2023 Massachusetts PUMS, downloaded from the U.S. Census Bureau:
  `https://www2.census.gov/programs-surveys/acs/data/pums/2023/1-Year/csv_pma.zip`
- Tennessee LODES8 2022 OD data and crosswalk, downloaded from LEHD:
  `https://lehd.ces.census.gov/data/lodes/LODES8/tn/od/tn_od_main_JT00_2022.csv.gz`
  and `https://lehd.ces.census.gov/data/lodes/LODES8/tn/tn_xwalk.csv.gz`
- 20 Newsgroups training data, downloaded through scikit-learn's
  `fetch_20newsgroups` dataset helper. The fixed document and term lists used
  for the paper's four-topic matrix are encoded in the preprocessing script.
- MBTA static GTFS feed, downloaded from MBTA:
  `https://cdn.mbta.com/MBTA_GTFS.zip`
- OpenAlex author-topic payload, queried from the OpenAlex API:
  `https://api.openalex.org/authors`

The preprocessing code is intentionally kept with the experiment scripts:

```text
main_experiment/real_world_benchmark/_real_world_data.py
```
