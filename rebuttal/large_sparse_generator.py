"""CSR-native synthetic generator for the rebuttal scaling sweep.

The dense generator in ``main_experiment/synthetic_benchmark/run_synthetic_evaluation.py``
reproduces the five benchmark families faithfully, but it cannot be scaled up:
``build_case`` materialises an ``n_rows x n_cols`` dense buffer and sizes every
prototype as a *fraction of its subgroup*, so ``nnz`` grows like ``O(N^2 / 25)``
(about 80 GB of dense buffer at 1e5 vertices).

This module reproduces the same *family semantics* -- five super-blocks, two row
subgroups per super-block, one unique prototype per row subgroup, a paired
prototype drawn from the sibling row subgroup, a shared super-prototype, a cross
prototype aimed at the next super-block, Poisson-weighted hits and integer
off-target noise -- while emitting sparse triples directly. The one deliberate
change is the sizing rule:

    prototype sizes and hit counts derive from ``avg_row_degree``, not from the
    subgroup size. A row is attached to ``n_proto_active`` active prototypes and
    receives ``Poisson(avg_row_degree / n_proto_active)`` hits in each, so
    ``E[degree] ~= avg_row_degree`` and ``E[nnz] ~= n_rows * avg_row_degree``
    regardless of ``N``.

That makes the sweep linear in ``N`` for both time and memory, and it is exactly
what forces a dense baseline to run out of memory rather than trivially
decomposing into small blocks.

Connectivity
------------
The random structure alone leaves families A and B as disconnected super-blocks
(and family A as two components per super-block). A small deterministic spanning
backbone -- about twenty extra entries, negligible against ``n_rows *
avg_row_degree`` -- links every prototype set inside a super-block and chains
consecutive super-blocks, so the active support is a single connected component
for every family. Every row also receives at least one hit in its primary
prototype, so no row is left isolated. ``active_components`` reports the realised
count and the largest component's vertex count so callers can assert this.

Matrix conventions match ``build_case``: the returned ``"matrix"`` is the
*observed* matrix, i.e. the ground truth shuffled by ``row_perm`` / ``col_perm``.
The unshuffled ground truth is not returned (it would double peak memory at the
largest sizes); recover it with ``observed[row_perm][:, col_perm]``.

The family parameters below are transcribed from
``main_experiment/synthetic_benchmark/run_synthetic_evaluation.py`` (see
``DENSE_SEMANTICS_SOURCE``); keep the two in step.

Usage
-----
    python rebuttal/large_sparse_generator.py                  # full sweep
    python rebuttal/large_sparse_generator.py --specs x1e4     # one size
"""

from __future__ import annotations

import argparse
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components

REPO_ROOT = Path(__file__).resolve().parents[1]
DENSE_SEMANTICS_SOURCE = (
    REPO_ROOT / "main_experiment" / "synthetic_benchmark" / "run_synthetic_evaluation.py"
)

SUPER_BLOCKS = 5
ROW_SUBGROUPS_PER_SUPER = 2

# A prototype holds this many columns per unit of ``avg_row_degree`` (modulated by
# the family's keep fraction). The duty cycle is at most one hit per column, so
# anything comfortably above 1 keeps the without-replacement draw uncensored.
PROTOTYPE_DEGREE_MULTIPLE = 3.0

# Lower bounds for the multinomial group splits; they only guard against an empty
# sub-group at the smallest size, not against the ordinary near-uniform split.
ROW_SUBGROUP_MIN = 4
COL_SUBGROUP_MIN = 4


@dataclass(frozen=True)
class ScaleSpec:
    """One point of the scaling sweep."""

    key: str
    name: str
    n_rows: int
    n_cols: int
    avg_row_degree: float = 12.0


SIZE_SPECS: tuple[ScaleSpec, ...] = (
    ScaleSpec("x1e3", "1.0k x 0.9k", 1_000, 900),
    ScaleSpec("x1e4", "10k x 9k", 10_000, 9_000),
    ScaleSpec("x3e4", "30k x 27k", 30_000, 27_000),
    ScaleSpec("x1e5", "100k x 90k", 100_000, 90_000),
    ScaleSpec("x3e5", "300k x 270k", 300_000, 270_000),
)


@dataclass(frozen=True)
class _FamilySpec:
    """Sparse-side family parameters (subset of the dense ``FamilySpec``)."""

    key: str
    col_subgroups_per_super: int
    unique_keep: float
    unique_prob: float
    unique_lambda: float
    paired_keep: float = 0.0
    paired_prob: float = 0.0
    paired_lambda: float = 0.0
    shared_keep: float = 0.0
    shared_prob: float = 0.0
    shared_lambda: float = 0.0
    cross_keep: float = 0.0
    cross_prob: float = 0.0
    cross_lambda: float = 0.0
    off_target_choices: tuple[int, int] = (0, 0)
    off_target_lambda: float = 1.0


_FAMILIES: tuple[_FamilySpec, ...] = (
    _FamilySpec(
        key="clean_block",
        col_subgroups_per_super=2,
        unique_keep=0.82,
        unique_prob=0.72,
        unique_lambda=5.0,
        off_target_choices=(0, 1),
        off_target_lambda=1.0,
    ),
    _FamilySpec(
        key="paired_overlap",
        col_subgroups_per_super=2,
        unique_keep=0.78,
        unique_prob=0.60,
        unique_lambda=4.8,
        paired_keep=0.72,
        paired_prob=0.34,
        paired_lambda=3.0,
        off_target_choices=(0, 1),
        off_target_lambda=1.0,
    ),
    _FamilySpec(
        key="shared_super",
        col_subgroups_per_super=3,
        unique_keep=0.74,
        unique_prob=0.54,
        unique_lambda=4.4,
        shared_keep=0.78,
        shared_prob=0.42,
        shared_lambda=4.0,
        off_target_choices=(0, 1),
        off_target_lambda=1.0,
    ),
    _FamilySpec(
        key="shared_super_noisy",
        col_subgroups_per_super=3,
        unique_keep=0.72,
        unique_prob=0.40,
        unique_lambda=4.2,
        paired_keep=0.60,
        paired_prob=0.0,
        paired_lambda=2.0,
        shared_keep=0.78,
        shared_prob=0.48,
        shared_lambda=4.8,
        off_target_choices=(0, 1),
        off_target_lambda=1.0,
    ),
    _FamilySpec(
        key="cross_block_leakage",
        col_subgroups_per_super=3,
        unique_keep=0.72,
        unique_prob=0.48,
        unique_lambda=4.0,
        shared_keep=0.74,
        shared_prob=0.40,
        shared_lambda=4.0,
        cross_keep=0.68,
        cross_prob=0.12,
        cross_lambda=2.2,
        off_target_choices=(0, 1),
        off_target_lambda=1.1,
    ),
)

FAMILY_KEYS: tuple[str, ...] = tuple(family.key for family in _FAMILIES)

_FAMILY_BY_KEY: dict[str, _FamilySpec] = {family.key: family for family in _FAMILIES}

SIZE_SPEC_BY_KEY: dict[str, ScaleSpec] = {spec.key: spec for spec in SIZE_SPECS}


# --------------------------------------------------------------------------- #
# Structure helpers (mirrors of the dense generator's helpers)
# --------------------------------------------------------------------------- #


def _split_sizes(
    total: int,
    n_groups: int,
    minimum: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Split ``total`` into ``n_groups`` sizes, each at least ``minimum``."""
    base = np.full(n_groups, minimum, dtype=int)
    remaining = total - int(base.sum())
    if remaining < 0:
        msg = f"Invalid group-size configuration: total={total}, minimum={minimum}"
        raise ValueError(msg)
    if remaining == 0:
        return base
    allocation = rng.multinomial(remaining, np.full(n_groups, 1.0 / n_groups))
    return base + allocation


def _subgroup_intervals(sizes: np.ndarray) -> list[tuple[int, int]]:
    starts = np.concatenate([[0], np.cumsum(sizes[:-1])])
    return [(int(start), int(start + size)) for start, size in zip(starts, sizes, strict=True)]


def _prototype_size(avg_row_degree: float, keep_fraction: float) -> int:
    """Prototype width in columns, pinned to the row degree rather than to N."""
    return max(2, int(np.ceil(PROTOTYPE_DEGREE_MULTIPLE * avg_row_degree * keep_fraction)))


def _choose_prototype(
    cols: np.ndarray,
    size: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Uniform without-replacement subset of ``cols`` (sorted), capped at len(cols)."""
    if cols.size == 0 or size <= 0:
        return np.empty(0, dtype=cols.dtype)
    size = int(min(cols.size, size))
    return np.sort(rng.choice(cols, size=size, replace=False))


def _sample_without_replacement(
    proto: np.ndarray,
    counts: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    """Per-row uniform without-replacement subsets of ``proto``.

    ``counts[r]`` columns are drawn for row ``r``; the result is the flat,
    row-major concatenation of the selected columns. Uses the random-key trick so
    the whole row subgroup is drawn with numpy rather than one row at a time.
    """
    n_rows = counts.size
    width = proto.size
    counts = np.minimum(counts, width)
    keys = rng.random((n_rows, width))
    order = np.argsort(keys, axis=1)
    take = np.arange(width, dtype=np.int64)[None, :] < counts[:, None]
    return proto[order[take]]


# --------------------------------------------------------------------------- #
# Generator
# --------------------------------------------------------------------------- #


def build_sparse_case(family_key: str, spec: ScaleSpec, seed: int) -> dict:
    """Build one sparse case; ``"matrix"`` is the observed (shuffled) matrix.

    Returns the contract keys ``matrix``, ``row_perm``, ``col_perm``,
    ``row_labels``, ``col_labels``, ``family_key``, ``spec``, ``seed``, ``nnz``,
    plus ``row_super_labels`` / ``col_super_labels`` / ``n_components`` /
    ``largest_component`` for parity with ``build_case`` and the connectivity
    assertion.
    """
    try:
        family = _FAMILY_BY_KEY[family_key]
    except KeyError:
        raise KeyError(
            f"unknown family {family_key!r}; expected one of {FAMILY_KEYS}"
        ) from None

    rng = np.random.default_rng(seed)
    n_rows = int(spec.n_rows)
    n_cols = int(spec.n_cols)
    avg_row_degree = float(spec.avg_row_degree)

    # --- group structure ---------------------------------------------------
    row_super_sizes = _split_sizes(
        n_rows, SUPER_BLOCKS, ROW_SUBGROUPS_PER_SUPER * ROW_SUBGROUP_MIN, rng
    )
    col_super_sizes = _split_sizes(
        n_cols, SUPER_BLOCKS, family.col_subgroups_per_super * COL_SUBGROUP_MIN, rng
    )

    row_sub_sizes: list[int] = []
    for super_size in row_super_sizes:
        row_sub_sizes.extend(
            int(size)
            for size in _split_sizes(
                int(super_size), ROW_SUBGROUPS_PER_SUPER, ROW_SUBGROUP_MIN, rng
            )
        )
    col_sub_sizes: list[int] = []
    for super_size in col_super_sizes:
        col_sub_sizes.extend(
            int(size)
            for size in _split_sizes(
                int(super_size), family.col_subgroups_per_super, COL_SUBGROUP_MIN, rng
            )
        )

    row_sub_ranges = _subgroup_intervals(np.array(row_sub_sizes, dtype=int))
    col_sub_ranges = _subgroup_intervals(np.array(col_sub_sizes, dtype=int))

    row_labels = np.empty(n_rows, dtype=np.int64)
    row_super_labels = np.empty(n_rows, dtype=np.int64)
    for global_sub, (start, stop) in enumerate(row_sub_ranges):
        row_labels[start:stop] = global_sub
        row_super_labels[start:stop] = global_sub // ROW_SUBGROUPS_PER_SUPER

    col_labels = np.empty(n_cols, dtype=np.int64)
    col_super_labels = np.empty(n_cols, dtype=np.int64)
    for global_sub, (start, stop) in enumerate(col_sub_ranges):
        col_labels[start:stop] = global_sub
        col_super_labels[start:stop] = global_sub // family.col_subgroups_per_super

    # --- prototypes --------------------------------------------------------
    unique_prototypes: dict[tuple[int, int], np.ndarray] = {}
    shared_prototypes: dict[int, np.ndarray] = {}
    for super_id in range(SUPER_BLOCKS):
        col_base = super_id * family.col_subgroups_per_super
        for local_sub in range(ROW_SUBGROUPS_PER_SUPER):
            cols = np.arange(*col_sub_ranges[col_base + local_sub], dtype=np.int64)
            unique_prototypes[(super_id, local_sub)] = _choose_prototype(
                cols, _prototype_size(avg_row_degree, family.unique_keep), rng
            )
        if family.col_subgroups_per_super >= 3:
            shared_cols = np.arange(*col_sub_ranges[col_base + 2], dtype=np.int64)
            shared_prototypes[super_id] = _choose_prototype(
                shared_cols, _prototype_size(avg_row_degree, family.shared_keep), rng
            )
        else:
            shared_prototypes[super_id] = np.empty(0, dtype=np.int64)

    if family.cross_keep > 0:
        cross_prototypes = {
            super_id: shared_prototypes[(super_id + 1) % SUPER_BLOCKS]
            for super_id in range(SUPER_BLOCKS)
        }
    else:
        cross_prototypes = {super_id: np.empty(0, dtype=np.int64) for super_id in range(SUPER_BLOCKS)}

    # --- hits --------------------------------------------------------------
    row_chunks: list[np.ndarray] = []
    col_chunks: list[np.ndarray] = []
    val_chunks: list[np.ndarray] = []

    def emit(rows: np.ndarray, cols: np.ndarray, vals: np.ndarray) -> None:
        if rows.size:
            row_chunks.append(rows)
            col_chunks.append(cols)
            val_chunks.append(vals)

    backbone_rows: list[int] = []
    backbone_cols: list[int] = []

    for global_row_sub, (row_start, row_stop) in enumerate(row_sub_ranges):
        super_id = global_row_sub // ROW_SUBGROUPS_PER_SUPER
        local_sub = global_row_sub % ROW_SUBGROUPS_PER_SUPER
        pair_sub = 1 - local_sub
        n_subgroup_rows = row_stop - row_start

        primary = unique_prototypes[(super_id, local_sub)]
        paired = _choose_prototype(
            unique_prototypes[(super_id, pair_sub)],
            _prototype_size(avg_row_degree, family.paired_keep),
            rng,
        )
        shared = shared_prototypes[super_id]
        cross = cross_prototypes[super_id]

        # The four prototype roles, in the dense generator's order. Only those
        # with a positive hit probability carry mass.
        roles = (
            (primary, family.unique_prob, family.unique_lambda),
            (paired, family.paired_prob, family.paired_lambda),
            (shared, family.shared_prob, family.shared_lambda),
            (cross, family.cross_prob, family.cross_lambda),
        )
        active = [role for role in roles if role[1] > 0 and role[0].size > 0]
        if not active:
            continue
        hits_per_role = avg_row_degree / len(active)

        for role_index, (proto, _prob, lam) in enumerate(active):
            counts = rng.poisson(hits_per_role, size=n_subgroup_rows).astype(np.int64)
            if role_index == 0:
                # Guarantee every row touches its primary prototype, so no row is
                # left isolated (and the active support stays one component).
                np.maximum(counts, 1, out=counts)
            selected = _sample_without_replacement(proto, counts, rng)
            row_ids = np.repeat(
                np.arange(row_start, row_stop, dtype=np.int64), np.minimum(counts, proto.size)
            )
            weights = 1.0 + rng.poisson(lam, size=selected.size)
            emit(row_ids, selected, weights.astype(np.float64))

        # Integer off-target noise: exclude every prototype column (all four
        # roles, active or not) exactly as the dense generator does.
        noise_low, noise_high = family.off_target_choices
        if noise_high > 0:
            targeted = np.unique(np.concatenate([role[0] for role in roles if role[0].size]))
            keep_mask = np.ones(n_cols, dtype=bool)
            keep_mask[targeted] = False
            non_target = np.flatnonzero(keep_mask)
            if non_target.size:
                counts = rng.integers(noise_low, noise_high + 1, size=n_subgroup_rows)
                total = int(counts.sum())
                if total > 0:
                    _emit_noise(
                        rng,
                        row_start,
                        counts,
                        non_target,
                        family.off_target_lambda,
                        emit,
                    )

        # Deterministic spanning backbone (see module docstring).
        anchor0 = row_sub_ranges[super_id * ROW_SUBGROUPS_PER_SUPER][0]
        anchor1 = row_sub_ranges[super_id * ROW_SUBGROUPS_PER_SUPER + 1][0]
        for proto in (unique_prototypes[(super_id, 0)], unique_prototypes[(super_id, 1)], shared):
            if proto.size:
                backbone_rows.append(anchor0)
                backbone_cols.append(int(proto[0]))
        if unique_prototypes[(super_id, 1)].size:
            backbone_rows.append(anchor1)
            backbone_cols.append(int(unique_prototypes[(super_id, 1)][0]))
        next_super = (super_id + 1) % SUPER_BLOCKS
        next_primary = unique_prototypes[(next_super, 0)]
        if next_primary.size:
            backbone_rows.append(anchor0)
            backbone_cols.append(int(next_primary[0]))

    if backbone_rows:
        emit(
            np.array(backbone_rows, dtype=np.int64),
            np.array(backbone_cols, dtype=np.int64),
            np.ones(len(backbone_rows), dtype=np.float64),
        )

    rows_gt = np.concatenate(row_chunks)
    cols_gt = np.concatenate(col_chunks)
    vals_gt = np.concatenate(val_chunks)

    # --- shuffle to the observed matrix (same contract as build_case) -------
    row_perm = rng.permutation(n_rows)
    col_perm = rng.permutation(n_cols)
    inv_row = np.empty(n_rows, dtype=np.int64)
    inv_row[row_perm] = np.arange(n_rows, dtype=np.int64)
    inv_col = np.empty(n_cols, dtype=np.int64)
    inv_col[col_perm] = np.arange(n_cols, dtype=np.int64)

    observed = sp.coo_matrix(
        (vals_gt, (inv_row[rows_gt], inv_col[cols_gt])),
        shape=(n_rows, n_cols),
        dtype=np.float64,
    ).tocsr()
    observed.sum_duplicates()
    observed.eliminate_zeros()

    n_components, largest_component = active_components(observed)

    return {
        "matrix": observed,
        "row_perm": row_perm,
        "col_perm": col_perm,
        "row_labels": row_labels,
        "col_labels": col_labels,
        "row_super_labels": row_super_labels,
        "col_super_labels": col_super_labels,
        "row_sub_labels": row_labels,
        "col_sub_labels": col_labels,
        "family_key": family_key,
        "spec": spec,
        "seed": int(seed),
        "nnz": int(observed.nnz),
        "n_components": n_components,
        "largest_component": largest_component,
    }


def _emit_noise(
    rng: np.random.Generator,
    row_start: int,
    counts: np.ndarray,
    non_target: np.ndarray,
    lam: float,
    emit,
) -> None:
    """Append off-target noise hits for one row subgroup.

    The configured families draw 0 or 1 noise targets per row, which vectorises
    exactly (a global without-replacement draw is per-row without-replacement
    when every row takes at most one). A general per-row path covers denser
    off-target settings; it is never reached by the five benchmark families.
    """
    if counts.max() <= 1:
        selected_rows = np.flatnonzero(counts == 1)
        draw = int(min(selected_rows.size, non_target.size))
        if draw == 0:
            return
        cols = rng.choice(non_target, size=draw, replace=False)
        weights = 1.0 + rng.poisson(lam, size=draw)
        emit(selected_rows.astype(np.int64) + row_start, cols, weights.astype(np.float64))
        return

    for offset, count in enumerate(counts):
        count = int(min(int(count), non_target.size))
        if count <= 0:
            continue
        cols = rng.choice(non_target, size=count, replace=False)
        weights = 1.0 + rng.poisson(lam, size=count)
        emit(
            np.full(count, row_start + offset, dtype=np.int64),
            cols,
            weights.astype(np.float64),
        )


# --------------------------------------------------------------------------- #
# Diagnostics
# --------------------------------------------------------------------------- #


def active_components(matrix) -> tuple[int, int]:
    """Connected components of the *active* bipartite support.

    Zero-marginal rows and columns are dropped first (as in the density
    conventions of ``tw_components.support_components``), then the bipartite
    support graph is built and its weak components counted. Returns
    ``(count, largest vertex count)`` where a vertex is one active row or active
    column.
    """
    csr = matrix.tocsr(copy=True) if sp.issparse(matrix) else sp.csr_matrix(matrix)
    csr.sum_duplicates()
    csr.eliminate_zeros()

    if csr.nnz == 0:
        return 0, 0

    row_sums = np.asarray(csr.sum(axis=1)).ravel()
    col_sums = np.asarray(csr.sum(axis=0)).ravel()
    active_rows = np.flatnonzero(row_sums > 0)
    active_cols = np.flatnonzero(col_sums > 0)
    if active_rows.size == 0 or active_cols.size == 0:
        return 0, 0

    support = (csr[active_rows][:, active_cols] != 0).astype(np.int8)
    graph = sp.bmat([[None, support], [support.T, None]], format="csr", dtype=np.int8)
    n_components, labels = connected_components(graph, directed=False)

    sizes = np.bincount(labels)
    return int(n_components), int(sizes.max())


def _select_specs(raw: str | None) -> tuple[ScaleSpec, ...]:
    if not raw:
        return SIZE_SPECS
    wanted = [token.strip() for token in raw.split(",") if token.strip()]
    unknown = [key for key in wanted if key not in SIZE_SPEC_BY_KEY]
    if unknown:
        raise SystemExit(
            f"unknown spec key(s): {', '.join(unknown)}; "
            f"available: {', '.join(spec.key for spec in SIZE_SPECS)}"
        )
    return tuple(SIZE_SPEC_BY_KEY[key] for key in wanted)


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Sparse CSR sweep over the five synthetic families; prints nnz, "
            "density, component count and build time per case."
        ),
        epilog=f"family semantics mirror {DENSE_SEMANTICS_SOURCE.relative_to(REPO_ROOT)}",
    )
    parser.add_argument(
        "--specs",
        default=None,
        help="comma-separated size keys to run (e.g. x1e3,x1e4); default: all",
    )
    parser.add_argument("--seed", type=int, default=0, help="generator seed (default: 0)")
    args = parser.parse_args()

    specs = _select_specs(args.specs)

    print(f"seed={args.seed}  families={len(FAMILY_KEYS)}  specs={len(specs)}")
    header = (
        f"{'spec':<6} {'family':<20} {'shape':>13} {'nnz':>10} "
        f"{'density':>10} {'comps':>6} {'largest':>9} {'secs':>7}"
    )
    print(header)
    print("-" * len(header))

    for spec in specs:
        for family_key in FAMILY_KEYS:
            started = time.perf_counter()
            case = build_sparse_case(family_key, spec, args.seed)
            elapsed = time.perf_counter() - started

            nnz = case["nnz"]
            density = nnz / (spec.n_rows * spec.n_cols)
            shape = f"{spec.n_rows}x{spec.n_cols}"
            print(
                f"{spec.key:<6} {family_key:<20} {shape:>13} {nnz:>10} "
                f"{density:>10.3e} {case['n_components']:>6} "
                f"{case['largest_component']:>9} {elapsed:>7.3f}"
            )
        print("-" * len(header))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
