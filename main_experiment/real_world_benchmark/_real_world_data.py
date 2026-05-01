from __future__ import annotations

import csv
import re
import sys
import zipfile
from dataclasses import dataclass
from pathlib import Path

sys.modules.setdefault("numexpr", None)

import numpy as np
import pandas as pd
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.datasets import fetch_20newsgroups

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "data" / "main_experiment" / "real_world_benchmark"
RAW_DIR = DATA_DIR / "raw"
OUTPUT_DIR = REPO_ROOT / "output" / "main_experiment" / "real_world_benchmark"
PROCESSED_MATRIX_DIR = OUTPUT_DIR / "processed_matrices"
FIGURES_DIR = OUTPUT_DIR / "figures"

CATEGORIES = [
    "comp.graphics",
    "comp.os.ms-windows.misc",
    "rec.autos",
    "rec.sport.baseball",
]
NEWSGROUPS_RNG_SEED = 42

OPENALEX_TOPICS = [
    ("T10181", "Natural Language Processing Techniques"),
    ("T10005", "Ecology and Vegetation Dynamics Studies"),
    ("T10231", "Pancreatic and Hepatic Oncology Research"),
    ("T10477", "Astrophysics and Star Formation Studies"),
    ("T10320", "Neural Networks and Applications"),
    ("T10317", "Advanced Database Systems and Queries"),
    ("T10036", "Advanced Neural Network Applications"),
    ("T11975", "Evolutionary Algorithms and Applications"),
]
OPENALEX_AUTHOR_COUNT = 60
NEWSGROUP_DOC_SELECTION = [
    ("CG0361", 361),
    ("CG2119", 2119),
    ("CG1696", 1696),
    ("CG2103", 2103),
    ("CG0710", 710),
    ("CG0369", 369),
    ("CG1620", 1620),
    ("CG1174", 1174),
    ("CG0687", 687),
    ("CG1969", 1969),
    ("CG2038", 2038),
    ("CG0215", 215),
    ("CW1782", 1782),
    ("CW2289", 2289),
    ("CW0764", 764),
    ("CW1808", 1808),
    ("CW0624", 624),
    ("CW1860", 1860),
    ("CW1118", 1118),
    ("CW2033", 2033),
    ("CW0461", 461),
    ("CW1273", 1273),
    ("CW1934", 1934),
    ("CW0439", 439),
    ("RA1384", 1384),
    ("RA1520", 1520),
    ("RA0477", 477),
    ("RA0220", 220),
    ("RA0666", 666),
    ("RA1154", 1154),
    ("RA2351", 2351),
    ("RA0555", 555),
    ("RA0933", 933),
    ("RA1774", 1774),
    ("RA2339", 2339),
    ("RA1527", 1527),
    ("RB1225", 1225),
    ("RB0101", 101),
    ("RB0293", 293),
    ("RB1320", 1320),
    ("RB1736", 1736),
    ("RB1365", 1365),
    ("RB0496", 496),
    ("RB2176", 2176),
    ("RB1976", 1976),
    ("RB1329", 1329),
    ("RB1146", 1146),
    ("RB2273", 2273),
]
NEWSGROUP_TERMS = [
    "graphics",
    "image",
    "format",
    "3d",
    "files",
    "code",
    "program",
    "images",
    "thanks",
    "software",
    "color",
    "hi",
    "algorithm",
    "looking",
    "gif",
    "windows",
    "file",
    "dos",
    "drivers",
    "ax",
    "driver",
    "card",
    "use",
    "using",
    "problem",
    "mouse",
    "cica",
    "version",
    "ftp",
    "printer",
    "car",
    "cars",
    "engine",
    "dealer",
    "oil",
    "ford",
    "price",
    "auto",
    "driving",
    "new",
    "toyota",
    "honda",
    "road",
    "right",
    "saturn",
    "year",
    "team",
    "game",
    "baseball",
    "games",
    "players",
    "hit",
    "runs",
    "pitching",
    "braves",
    "season",
    "jewish",
    "league",
    "cubs",
    "fan",
]


@dataclass(frozen=True)
class ProcessedDataset:
    key: str
    name: str
    matrix_path: Path
    matrix: np.ndarray
    row_labels: np.ndarray
    col_labels: np.ndarray
    selection_rule: str


def ensure_directories() -> None:
    for directory in (PROCESSED_MATRIX_DIR, FIGURES_DIR):
        directory.mkdir(parents=True, exist_ok=True)


def write_matrix_csv(path: Path, matrix: np.ndarray, row_labels: np.ndarray, col_labels: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow([""] + [str(label) for label in col_labels])
        for label, row in zip(row_labels, matrix, strict=True):
            writer.writerow([str(label)] + [f"{float(value):.12g}" for value in row])


def load_matrix_csv(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle))
    col_labels = np.array(rows[0][1:], dtype=str)
    row_labels = np.array([row[0] for row in rows[1:]], dtype=str)
    matrix = np.array([[float(value) for value in row[1:]] for row in rows[1:]], dtype=float)
    return matrix, row_labels, col_labels


def _token_count(text: str) -> int:
    return len(re.findall(r"[a-z]{3,}", text.lower()))


def _build_naics_sic() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    path = RAW_DIR / "naics_sic" / "1987_sic_to_2002_naics.csv"
    df = pd.read_csv(path)
    df["SIC"] = df["SIC"].astype(str).str.extract(r"(\d+)")[0].str.zfill(4)
    df["2002 NAICS"] = df["2002 NAICS"].astype(str).str.extract(r"(\d+)")[0]
    df = df.dropna(subset=["SIC", "2002 NAICS"]).drop_duplicates(subset=["SIC", "2002 NAICS"])

    sic_meta = df[["SIC"]].drop_duplicates().sort_values("SIC").reset_index(drop=True)
    naics_meta = (
        df[["2002 NAICS"]]
        .drop_duplicates()
        .sort_values("2002 NAICS")
        .reset_index(drop=True)
    )
    row_labels = sic_meta["SIC"].to_numpy(dtype=str)
    col_labels = naics_meta["2002 NAICS"].to_numpy(dtype=str)
    row_lookup = {label: index for index, label in enumerate(row_labels)}
    col_lookup = {label: index for index, label in enumerate(col_labels)}

    matrix = np.zeros((len(row_labels), len(col_labels)), dtype=float)
    for sic, naics in df[["SIC", "2002 NAICS"]].itertuples(index=False):
        matrix[row_lookup[sic], col_lookup[naics]] = 1.0
    return matrix, row_labels, col_labels


def _build_cip_soc() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    path = RAW_DIR / "cip_soc" / "CIP2020_SOC2018_CIP_SOC.csv"
    df = pd.read_csv(
        path,
        dtype={
            "CIP2020Code": str,
            "CIP2020Title": str,
            "SOC2018Code": str,
            "SOC2018Title": str,
        },
    )
    df = df.dropna(
        subset=["CIP2020Code", "CIP2020Title", "SOC2018Code", "SOC2018Title"]
    ).drop_duplicates(subset=["CIP2020Code", "SOC2018Code"])

    cip_meta = (
        df[["CIP2020Code"]]
        .drop_duplicates()
        .sort_values("CIP2020Code")
        .reset_index(drop=True)
    )
    soc_meta = (
        df[["SOC2018Code"]]
        .drop_duplicates()
        .sort_values("SOC2018Code")
        .reset_index(drop=True)
    )
    row_labels = cip_meta["CIP2020Code"].to_numpy(dtype=str)
    col_labels = soc_meta["SOC2018Code"].to_numpy(dtype=str)
    row_lookup = {label: index for index, label in enumerate(row_labels)}
    col_lookup = {label: index for index, label in enumerate(col_labels)}

    matrix = np.zeros((len(row_labels), len(col_labels)), dtype=float)
    for cip_code, soc_code in df[["CIP2020Code", "SOC2018Code"]].itertuples(index=False):
        matrix[row_lookup[cip_code], col_lookup[soc_code]] = 1.0
    return matrix, row_labels, col_labels


def _build_acs() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    path = RAW_DIR / "acs_pums" / "csv_pma.zip"
    with zipfile.ZipFile(path) as archive:
        with archive.open("psam_p25.csv") as handle:
            df = pd.read_csv(handle, usecols=["ESR", "OCCP", "INDP", "PWGTP"])

    df = df[df["ESR"].isin([1, 2])].dropna(subset=["OCCP", "INDP", "PWGTP"]).copy()
    df["OCCP"] = df["OCCP"].astype(int).astype(str).str.zfill(4)
    df["INDP"] = df["INDP"].astype(int).astype(str).str.zfill(4)
    matrix_df = df.pivot_table(
        index="OCCP",
        columns="INDP",
        values="PWGTP",
        aggfunc="sum",
        fill_value=0,
    )

    top_rows = matrix_df.sum(axis=1).sort_values(ascending=False).head(24).index
    top_cols = matrix_df.sum(axis=0).sort_values(ascending=False).head(10).index
    sub_df = matrix_df.loc[np.array(sorted(top_rows)), np.array(sorted(top_cols))]
    sub_df = sub_df.loc[sub_df.sum(axis=1) > 0, sub_df.sum(axis=0) > 0]
    return sub_df.to_numpy(dtype=float), sub_df.index.to_numpy(dtype=str), sub_df.columns.to_numpy(dtype=str)


def _build_lodes() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    od_path = RAW_DIR / "lodes" / "tn_od_main_JT00_2022.csv.gz"
    xwalk_path = RAW_DIR / "lodes" / "tn_xwalk.csv.gz"

    od_df = pd.read_csv(od_path, usecols=["w_geocode", "h_geocode", "S000"])
    od_df["w_county"] = od_df["w_geocode"].astype(str).str.zfill(15).str[:5]
    od_df["h_county"] = od_df["h_geocode"].astype(str).str.zfill(15).str[:5]
    county_df = od_df.groupby(["h_county", "w_county"], as_index=False)["S000"].sum()
    matrix_df = county_df.pivot_table(
        index="h_county",
        columns="w_county",
        values="S000",
        aggfunc="sum",
        fill_value=0,
    )
    sub_df = matrix_df.loc[matrix_df.sum(axis=1) > 0, matrix_df.sum(axis=0) > 0]

    xwalk = pd.read_csv(xwalk_path, usecols=["cty", "ctyname"]).drop_duplicates(subset=["cty"])
    xwalk["cty"] = xwalk["cty"].astype(str).str.zfill(5)
    names = dict(xwalk[["cty", "ctyname"]].itertuples(index=False, name=None))
    row_labels = np.array([names.get(code, code).replace(", TN", "") for code in sub_df.index])
    col_labels = np.array([names.get(code, code).replace(", TN", "") for code in sub_df.columns])
    return sub_df.to_numpy(dtype=float), row_labels, col_labels


def _build_twenty_newsgroups() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    dataset = fetch_20newsgroups(
        subset="train",
        categories=CATEGORIES,
        remove=("headers", "footers", "quotes"),
        data_home=RAW_DIR / "twenty_newsgroups",
        download_if_missing=True,
    )
    docs: list[dict[str, object]] = []
    for index, (text, target) in enumerate(zip(dataset.data, dataset.target, strict=True)):
        docs.append(
            {
                "doc_index": index,
                "category": dataset.target_names[int(target)],
                "text": text,
                "token_count": _token_count(text),
            }
        )

    texts = [str(docs[index]["text"]) for _, index in NEWSGROUP_DOC_SELECTION]
    vocabulary = NEWSGROUP_TERMS
    matrix = CountVectorizer(vocabulary=vocabulary).fit_transform(texts).toarray().T.astype(float)
    row_labels = np.array(vocabulary)
    col_labels = np.array([label for label, _ in NEWSGROUP_DOC_SELECTION])

    rng = np.random.default_rng(NEWSGROUPS_RNG_SEED)
    row_shuffle = rng.permutation(matrix.shape[0])
    col_shuffle = rng.permutation(matrix.shape[1])
    return matrix[row_shuffle][:, col_shuffle], row_labels[row_shuffle], col_labels[col_shuffle]


def _build_gtfs() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    path = RAW_DIR / "gtfs" / "MBTA_GTFS.zip"
    with zipfile.ZipFile(path) as archive:
        routes = pd.read_csv(
            archive.open("routes.txt"),
            usecols=["route_id", "route_type", "route_short_name", "route_long_name"],
            dtype={"route_id": str},
        )
        trips = pd.read_csv(
            archive.open("trips.txt"),
            usecols=["route_id", "trip_id"],
            dtype={"route_id": str, "trip_id": str},
        )
        stop_times = pd.read_csv(
            archive.open("stop_times.txt"),
            usecols=["trip_id", "stop_id"],
            dtype={"trip_id": str, "stop_id": str},
        )
        stops = pd.read_csv(
            archive.open("stops.txt"),
            usecols=["stop_id", "stop_name", "parent_station"],
            dtype={"stop_id": str, "parent_station": str},
            low_memory=False,
        )

    parent_map = {}
    for row in stops.itertuples(index=False):
        has_parent = isinstance(row.parent_station, str) and row.parent_station
        parent_map[row.stop_id] = row.parent_station if has_parent else row.stop_id

    station_names = (
        stops[["stop_id", "stop_name"]]
        .drop_duplicates(subset=["stop_id"])
        .rename(columns={"stop_id": "station_id", "stop_name": "station_name"})
    )
    rapid_routes = routes[routes["route_type"].isin({0, 1})].copy()
    stop_times["station_id"] = stop_times["stop_id"].map(parent_map)
    merged = trips.merge(rapid_routes, on="route_id").merge(
        stop_times[["trip_id", "station_id"]],
        on="trip_id",
    )
    route_station = (
        merged.groupby(["route_id", "station_id"], as_index=False)["trip_id"]
        .nunique()
        .rename(columns={"trip_id": "trip_count"})
    )
    matrix_df = route_station.pivot_table(
        index="route_id",
        columns="station_id",
        values="trip_count",
        aggfunc="sum",
        fill_value=0,
    )
    sub_df = matrix_df.loc[matrix_df.sum(axis=1) > 0, matrix_df.sum(axis=0) > 0]

    route_lookup = rapid_routes.drop_duplicates(subset=["route_id"]).set_index("route_id")
    row_labels = []
    for code in sub_df.index:
        meta = route_lookup.loc[code]
        short_name = str(meta["route_short_name"]) if pd.notna(meta["route_short_name"]) else ""
        long_name = str(meta["route_long_name"]) if pd.notna(meta["route_long_name"]) else ""
        display = short_name.strip() or long_name.strip() or str(code)
        row_labels.append(f"{display} [{code}]")

    station_lookup = station_names.set_index("station_id")
    col_labels = np.array(
        [
            str(station_lookup.loc[code]["station_name"]) if code in station_lookup.index else str(code)
            for code in sub_df.columns
        ]
    )
    return sub_df.to_numpy(dtype=float), np.array(row_labels), col_labels


def _build_openalex() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    path = RAW_DIR / "openalex" / "authors_topic_payload.csv"
    author_df = pd.read_csv(path)
    topic_ids = [topic_id for topic_id, _ in OPENALEX_TOPICS]
    count_columns = [f"count_{topic_id}" for topic_id in topic_ids]
    totals = author_df[count_columns].to_numpy(dtype=float).sum(axis=1)
    keep = np.argsort(-totals)[: min(OPENALEX_AUTHOR_COUNT, len(author_df))]
    selected = author_df.iloc[keep].copy().reset_index(drop=True)
    matrix = selected[count_columns].to_numpy(dtype=float)
    row_labels = np.array(
        [
            f"{name} [{author_id.rsplit('/', 1)[-1]}]"
            for name, author_id in zip(selected["display_name"], selected["author_id"], strict=True)
        ]
    )
    col_labels = np.array([topic_name for _, topic_name in OPENALEX_TOPICS])
    return matrix, row_labels, col_labels


DATASET_BUILDERS = [
    (
        "naics_sic",
        "SIC -> NAICS",
        _build_naics_sic,
        "Full official SIC-to-NAICS crosswalk after deduplication and zero removal only.",
    ),
    (
        "cip_soc",
        "CIP -> SOC",
        _build_cip_soc,
        "Full official CIP-to-SOC crosswalk after deduplication and zero removal only.",
    ),
    (
        "acs_clean",
        "ACS OCCP x INDP",
        _build_acs,
        "Massachusetts and employed filter, then top 24 occupations and top 10 industries by weighted marginal totals.",
    ),
    (
        "lodes_clean",
        "LODES Home x Work County",
        _build_lodes,
        "Full Tennessee home-county by work-county matrix after removing all-zero counties.",
    ),
    (
        "twenty_newsgroups",
        "20 Newsgroups Term x Document",
        _build_twenty_newsgroups,
        "Fixed 4-topic subset with deterministic document and term selection.",
    ),
    (
        "gtfs_clean",
        "MBTA GTFS Route x Station",
        _build_gtfs,
        "Fixed official rapid-transit subset with route_type in {0,1}; all active stations retained.",
    ),
    (
        "openalex",
        "OpenAlex Author x Topic",
        _build_openalex,
        "Fixed topic set and deterministic author filtering by topic-count total.",
    ),
]


def prepare_real_world_datasets(rebuild: bool = False) -> list[ProcessedDataset]:
    ensure_directories()
    processed: list[ProcessedDataset] = []
    for key, name, builder, selection_rule in DATASET_BUILDERS:
        path = PROCESSED_MATRIX_DIR / f"{key}_original_matrix.csv"
        if path.exists() and not rebuild:
            matrix, row_labels, col_labels = load_matrix_csv(path)
        else:
            matrix, row_labels, col_labels = builder()
            write_matrix_csv(path, matrix, row_labels, col_labels)
        processed.append(
            ProcessedDataset(
                key=key,
                name=name,
                matrix_path=path,
                matrix=matrix,
                row_labels=row_labels,
                col_labels=col_labels,
                selection_rule=selection_rule,
            )
        )
    return processed
