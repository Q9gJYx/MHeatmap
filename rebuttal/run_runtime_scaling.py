"""Runtime and memory measurement harness for workstream A2.

Answers wPzm W3/Q2/Q4, JfQK W5/Q3 and XFNF W6: wall-clock and peak memory for the
spectral methods, on the real matrices, the controlled synthetic benchmark, and a new
sparse scaling sweep that reaches 300k x 270k.

Measurement methodology
-----------------------
Three mechanisms, each used for what it is actually good at.

* **Wall clock** -- ``time.perf_counter`` per phase, with ``time.process_time`` alongside
  so the cpu/wall ratio exposes BLAS threading.
* **Peak memory per cell** -- ``resource.getrusage(RUSAGE_SELF).ru_maxrss`` is a monotonic,
  process-wide high-water mark, so it is only attributable if the cell owns a fresh
  process.  Hence one ``spawn`` subprocess per (instance, method, solver, rep).
* **Per-phase memory** -- ``tracemalloc`` cannot see numpy buffers, so it is useless here.
  ``psutil`` samples RSS from a daemon thread every 5 ms instead, and ``memory_full_info``
  supplies USS, which excludes pages shared with the BLAS/numpy libraries and therefore
  attributes the algorithm's own footprint more honestly than RSS.

The analytic footprints (``operator_footprint_mb`` for the explicit Gram the released code
builds, ``csr_footprint_mb`` for the matrix-free path) are exact and noise-free, and carry
the headline memory claim; measured RSS corroborates them.

Noise control
-------------
``OMP_NUM_THREADS`` and friends are pinned to 1 **before numpy is imported**, in this
process and therefore in every spawned child.  One untimed warmup plus ``--reps`` timed
repetitions are recorded per cell, every repetition persisted, so medians and spread can
be computed downstream rather than a single number being reported.

Failure is data
---------------
The dense baseline is expected to die at large N.  A cell whose analytic dense footprint
exceeds ``--mem-cap-mb`` is still attempted once under an ``RLIMIT_AS`` cap so the failure
is *measured* (``status="oom"``), and only skipped as ``predicted_oom`` when the footprint
exceeds a hard ceiling.  The dense curve therefore terminates at a real N while the
matrix-free curve continues -- that pair is the scalability evidence.

Run:
    cd MHeatMap/src
    uv run python rebuttal/run_runtime_scaling.py --suites real
    uv run python rebuttal/run_runtime_scaling.py --suites scale --specs x1e3,x1e4,x1e5
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import multiprocessing as mp
import os
import platform
import resource
import sys
import threading
import time
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from pathlib import Path

# --------------------------------------------------------------------------- #
# Pin BLAS threads BEFORE numpy is imported, so the child inherits the setting.
# --------------------------------------------------------------------------- #

BLAS_VARS = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
)
PRIOR_THREAD_ENV = {name: os.environ.get(name) for name in BLAS_VARS}
for _name in BLAS_VARS:
    os.environ[_name] = "1"

import numpy as np  # noqa: E402
import psutil  # noqa: E402

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import tw_matrixfree as mf  # noqa: E402

SYNTHETIC_DIR = mf.SYNTHETIC_DIR
if str(SYNTHETIC_DIR) not in sys.path:
    sys.path.insert(0, str(SYNTHETIC_DIR))

from run_synthetic_evaluation import FAMILIES, SIZES, build_case  # noqa: E402

OUTPUT_DIR = HERE / "output"
MATRIX_DIR = (
    mf.REPO_ROOT
    / "output"
    / "main_experiment"
    / "real_world_benchmark"
    / "processed_matrices"
)

DATASET_ORDER = (
    "naics_sic",
    "cip_soc",
    "acs_clean",
    "lodes_clean",
    "twenty_newsgroups",
    "gtfs_clean",
    "openalex",
)
DATASET_LABELS = {
    "naics_sic": "SIC->NAICS",
    "cip_soc": "CIP->SOC",
    "acs_clean": "ACS",
    "lodes_clean": "LODES",
    "twenty_newsgroups": "20Newsgroups",
    "gtfs_clean": "MBTA",
    "openalex": "OpenAlex",
}

METHODS = ("One-walk", "Two-walk", "Two-walk-Auto")
SOLVERS = ("dense", "matrix_free")
ALPHAS = mf.DEFAULT_ALPHA_GRID
WIDTHS = mf.bm.DEFAULT_WIDTH_GRID

# Synthetic sizes are named small/medium/large in the dense generator.
SIZE_LABELS = {size.key: size.name for size in SIZES}

RUN_FIELDS = [
    "instance_id", "suite", "family_key", "size_key", "seed",
    "m", "n", "nnz", "density",
    "n_components", "max_component_vertices", "max_component_nnz",
    "method", "alpha", "alpha_sweep", "solver", "repeat_index", "is_warmup", "isolation",
    "wall_s", "cpu_s", "peak_rss_mb", "peak_uss_mb", "ru_maxrss_mb",
    "operator_footprint_mb", "csr_footprint_mb",
    "n_solves", "total_matvecs", "max_iters", "mean_iters",
    "fiedler_lambda2", "fiedler_cond",
    "status", "error_msg", "host_id", "git_rev", "timestamp",
]

PHASE_FIELDS = [
    "instance_id", "suite", "method", "solver", "repeat_index",
    "phase", "wall_s", "cpu_s", "rss_peak_mb", "uss_peak_mb", "n_matvec",
]

PHASES = ("load", "reorder", "score", "total")


# --------------------------------------------------------------------------- #
# Cells
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class Cell:
    suite: str  # "real" | "synthetic" | "scale"
    instance_id: str
    family_key: str
    size_key: str
    seed: int
    method: str
    solver: str

    def as_dict(self) -> dict:
        return asdict(self)


def method_alpha(method: str) -> float | None:
    return 1.0 if method == "Two-walk" else None


# --------------------------------------------------------------------------- #
# Memory sampling
# --------------------------------------------------------------------------- #


class MemorySampler:
    """Peak RSS/USS per phase, sampled from a daemon thread.

    ``memory_info().rss`` is cheap enough to poll every few milliseconds;
    ``memory_full_info().uss`` reads smaps and is expensive, so USS is taken once at the
    end of the measured region rather than on every tick.
    """

    def __init__(self, interval: float = 0.005):
        self.interval = interval
        self.process = psutil.Process()
        self.peak_rss = 0
        self.peak_uss = 0
        self._stop = threading.Event()
        self._thread: threading.Thread | None = None
        self._current: list[int] | None = None

    def _sample(self) -> None:
        try:
            rss = self.process.memory_info().rss
        except psutil.Error:
            return
        self.peak_rss = max(self.peak_rss, rss)
        if self._current is not None:
            self._current[0] = max(self._current[0], rss)

    def _run(self) -> None:
        while not self._stop.wait(self.interval):
            self._sample()

    def start(self) -> None:
        self._sample()
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._thread.start()

    def stop(self) -> None:
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=1.0)
        self._sample()
        try:
            self.peak_uss = max(self.peak_uss, self.process.memory_full_info().uss)
        except psutil.Error:
            self.peak_uss = self.peak_rss

    @contextmanager
    def phase(self, name: str):
        slot = [0, 0]
        previous = self._current
        self._current = slot
        wall = time.perf_counter()
        cpu = time.process_time()
        try:
            yield slot
        finally:
            self._current = previous
            PHASE_ACCUMULATOR.append(
                {
                    "phase": name,
                    "wall_s": time.perf_counter() - wall,
                    "cpu_s": time.process_time() - cpu,
                    "rss_peak_mb": slot[0] / 1e6,
                    "uss_peak_mb": slot[1] / 1e6,
                }
            )


PHASE_ACCUMULATOR: list[dict] = []


# --------------------------------------------------------------------------- #
# Workers
# --------------------------------------------------------------------------- #


def _materialize(cell: Cell, sampler: MemorySampler):
    """Load or generate the instance; returns a dense ndarray or a CSR matrix."""
    if cell.suite == "real":
        # ``instance_id`` is the display label; ``family_key`` carries the file stem.
        matrix, _, _ = mf.bm.load_matrix_csv(
            MATRIX_DIR / f"{cell.family_key}_original_matrix.csv"
        )
        return matrix

    if cell.suite == "synthetic":
        family = next(f for f in FAMILIES if f.key == cell.family_key)
        size = next(s for s in SIZES if s.key == cell.size_key)
        return build_case(family, size, cell.seed)["observed"]

    from large_sparse_generator import SIZE_SPECS, build_sparse_case

    spec = next(s for s in SIZE_SPECS if s.key == cell.size_key)
    return build_sparse_case(cell.family_key, spec, cell.seed)["matrix"]


def _run_reorder(matrix, cell: Cell, sampler: MemorySampler, diagnostics: list):
    """Run the requested method, attributing reorder and scoring to separate phases."""
    if cell.method == "One-walk":
        with sampler.phase("reorder"):
            return mf.componentwise_reorder_sparse(
                matrix, mode="ow", solver=cell.solver, diagnostics=diagnostics,
                tol=1e-12,
            )

    if cell.method == "Two-walk":
        with sampler.phase("reorder"):
            return mf.componentwise_reorder_sparse(
                matrix, mode="tw", alpha=1.0, solver=cell.solver,
                diagnostics=diagnostics, tol=1e-12,
            )

    # Two-walk-Auto: replicate `tw_auto_reorder` so the sweep is attributed honestly.
    best = None
    for alpha in ALPHAS:
        with sampler.phase("reorder"):
            result = mf.componentwise_reorder_sparse(
                matrix, mode="tw", alpha=alpha, solver=cell.solver,
                diagnostics=diagnostics, tol=1e-12,
            )
        with sampler.phase("score"):
            score = float(mf.bm.mwb_auc(_scoring_view(result.matrix), WIDTHS))
        if (
            best is None
            or score > best[1] + 1e-12
            or (abs(score - best[1]) <= 1e-12 and alpha < best[0])
        ):
            best = (alpha, score, result)
    if best is None:
        raise RuntimeError("TW-Auto selected no alpha")
    return best[2]


def _scoring_view(matrix):
    """``mwb_auc`` needs a dense array; the scale suite never calls this."""
    return matrix.toarray() if hasattr(matrix, "toarray") else matrix


def _worker(cell_dict: dict, queue) -> None:
    global PHASE_ACCUMULATOR
    payload = dict(cell_dict)
    hard_limit = payload.pop("_hard_limit_bytes", 0)
    cell = Cell(**payload)
    PHASE_ACCUMULATOR = []

    baseline_rss = psutil.Process().memory_info().rss
    if hard_limit:
        # `hard_limit` is headroom *above the post-import address space*, not an absolute
        # cap: numpy and scipy map well over a gigabyte of virtual memory at import, so an
        # absolute RLIMIT_AS below that would kill the interpreter before it ran anything.
        current = psutil.Process().memory_info().vms
        resource.setrlimit(resource.RLIMIT_AS, (current + hard_limit, current + hard_limit))

    sampler = MemorySampler()
    diagnostics: list = []
    record = {
        **{key: "" for key in RUN_FIELDS},
        **cell.as_dict(),
        "isolation": "process",
        "host_id": _host_id(),
        "timestamp": time.time(),
        "status": "ok",
    }
    started = time.perf_counter()
    cpu_started = time.process_time()
    sampler.start()
    try:
        with sampler.phase("load"):
            matrix = _materialize(cell, sampler)
        # Record shape and the analytic footprints *before* solving: an OOM cell must
        # still carry the number that explains why it died.
        record.update(
            {
                "m": matrix.shape[0],
                "n": matrix.shape[1],
                "nnz": (
                    int(matrix.nnz)
                    if hasattr(matrix, "nnz")
                    else int((matrix > 0).sum())
                ),
            }
        )
        record.update(mf.footprints_sparse(matrix))
        with sampler.phase("total"):
            result = _run_reorder(matrix, cell, sampler, diagnostics)
        record.update(
            {
                "n_components": result.component_count,
                "fiedler_lambda2": _stat(diagnostics, "eigenvalue"),
                "fiedler_cond": _conditioning(diagnostics),
                "n_solves": len(diagnostics),
                "total_matvecs": sum(d.n_matvec for d in diagnostics),
                "max_iters": max((d.n_matvec for d in diagnostics), default=0),
                "mean_iters": _mean(d.n_matvec for d in diagnostics),
                "alpha": "" if method_alpha(cell.method) is None else method_alpha(cell.method),
                "alpha_sweep": ",".join(f"{a:g}" for a in ALPHAS) if cell.method == "Two-walk-Auto" else "",
            }
        )
    except MemoryError:
        record["status"] = "oom"
        record["error_msg"] = "MemoryError"
    except Exception as exc:  # noqa: BLE001 - the harness must record, not crash
        record["status"] = "error"
        record["error_msg"] = f"{type(exc).__name__}: {exc}"
    finally:
        sampler.stop()
        record["wall_s"] = time.perf_counter() - started
        record["cpu_s"] = time.process_time() - cpu_started
        record["peak_rss_mb"] = sampler.peak_rss / 1e6
        record["peak_uss_mb"] = sampler.peak_uss / 1e6
        record["ru_maxrss_mb"] = (
            resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0
        )
        record["baseline_rss_mb"] = baseline_rss / 1e6
        if record["nnz"]:
            record["density"] = record["nnz"] / (record["m"] * record["n"])

    phases = [
        {**{key: "" for key in PHASE_FIELDS}, **entry,
         "instance_id": cell.instance_id, "suite": cell.suite,
         "method": cell.method, "solver": cell.solver}
        for entry in PHASE_ACCUMULATOR
    ]
    queue.put((record, phases))


def _mean(values) -> float:
    values = list(values)
    return float(np.mean(values)) if values else 0.0


def _stat(diagnostics: list, attribute: str) -> float:
    if not diagnostics:
        return 0.0
    return float(np.median([getattr(d, attribute) for d in diagnostics]))


def _conditioning(diagnostics: list) -> float:
    """Median lambda_max / lambda_2, the ratio that governs Lanczos convergence."""
    ratios = [
        (d.next_eigenvalue / d.eigenvalue)
        for d in diagnostics
        if d.eigenvalue > 0 and np.isfinite(d.next_eigenvalue)
    ]
    return float(np.median(ratios)) if ratios else 0.0


def _run_cell(cell: Cell, timeout: float, hard_limit: int) -> tuple[dict, list[dict]]:
    """Run one cell in a fresh spawn process so its peak RSS is attributable."""
    context = mp.get_context("spawn")
    queue = context.Queue()
    payload = cell.as_dict()
    payload["_hard_limit_bytes"] = hard_limit
    process = context.Process(target=_worker, args=(payload, queue))
    process.start()
    process.join(timeout)

    if process.is_alive():
        process.terminate()
        process.join(5)
        if process.is_alive():
            process.kill()
        return (
            {
                **{key: "" for key in RUN_FIELDS},
                **cell.as_dict(),
                "isolation": "process",
                "status": "timeout",
                "error_msg": f"exceeded {timeout:g}s",
                "host_id": _host_id(),
                "timestamp": time.time(),
            },
            [],
        )

    if queue.empty():
        return (
            {
                **{key: "" for key in RUN_FIELDS},
                **cell.as_dict(),
                "isolation": "process",
                "status": "error",
                "error_msg": f"worker exited with code {process.exitcode}",
                "host_id": _host_id(),
                "timestamp": time.time(),
            },
            [],
        )
    return queue.get()


# --------------------------------------------------------------------------- #
# Cell enumeration
# --------------------------------------------------------------------------- #


def enumerate_cells(args) -> list[Cell]:
    suites = tuple(s.strip() for s in args.suites.split(",") if s.strip())
    methods = tuple(m.strip() for m in args.methods.split(",") if m.strip())
    cells: list[Cell] = []

    if "real" in suites:
        wanted_datasets = tuple(
            k.strip() for k in getattr(args, "datasets", "").split(",") if k.strip()
        )
        for key in DATASET_ORDER:
            if wanted_datasets and key not in wanted_datasets:
                continue
            for method in methods:
                for solver in SOLVERS:
                    cells.append(
                        Cell("real", DATASET_LABELS[key], key, "", 0, method, solver)
                    )

    if "synthetic" in suites:
        for family in FAMILIES:
            for size in SIZES:
                for seed in range(args.seeds):
                    for method in methods:
                        for solver in SOLVERS:
                            cells.append(
                                Cell(
                                    "synthetic",
                                    f"{family.key}/{size.key}/s{seed}",
                                    family.key,
                                    size.key,
                                    seed,
                                    method,
                                    solver,
                                )
                            )

    if "scale" in suites:
        from large_sparse_generator import FAMILY_KEYS, SIZE_SPECS

        wanted = tuple(s.strip() for s in args.specs.split(",") if s.strip())
        for family_key in FAMILY_KEYS:
            for spec in SIZE_SPECS:
                if wanted and spec.key not in wanted:
                    continue
                for method in methods:
                    for solver in SOLVERS:
                        cells.append(
                            Cell("scale", f"{family_key}/{spec.key}", family_key, spec.key, 0, method, solver)
                        )

    return cells


def predicted_dense_mb(cell: Cell) -> float:
    """Analytic explicit-Gram footprint, so OOM cells can be flagged before running."""
    if cell.suite == "scale":
        from large_sparse_generator import SIZE_SPECS

        spec = next(s for s in SIZE_SPECS if s.key == cell.size_key)
        return 2 * (spec.n_rows + spec.n_cols) ** 2 * 8 / 1e6
    return 0.0


# --------------------------------------------------------------------------- #
# Hardware block
# --------------------------------------------------------------------------- #


def _host_id(prefix: str = "machine") -> str:
    """Stable per-machine label that does not disclose the hostname.

    ``platform.node()`` typically embeds the account and hardware model, which would
    break the submission's anonymity the moment the hardware block is pasted into the
    paper.  A truncated digest still distinguishes two machines in a results file.
    """
    digest = hashlib.sha256(platform.node().encode("utf-8")).hexdigest()[:8]
    return f"{prefix}-{digest}"


def _scipy_version() -> str:
    try:
        import scipy

        return scipy.__version__
    except Exception:  # noqa: BLE001
        return "unknown"


def _cpu_model() -> str:
    try:
        with open("/proc/cpuinfo") as handle:
            for line in handle:
                if line.lower().startswith("model name"):
                    return line.split(":", 1)[1].strip()
    except OSError:
        pass
    return platform.processor() or "unknown"


def _total_ram_gb() -> float:
    try:
        with open("/proc/meminfo") as handle:
            for line in handle:
                if line.startswith("MemTotal"):
                    return int(line.split()[1]) / 1024 / 1024
    except OSError:
        pass
    return psutil.virtual_memory().total / 1e9


def _blas_info() -> tuple[str, str]:
    try:
        config = np.show_config(mode="dicts")
        blas = config.get("Build Dependencies", {}).get("blas", {})
        return blas.get("name", "unknown"), blas.get("version", "unknown")
    except Exception:  # noqa: BLE001 - show_config's schema has changed across versions
        return "unknown", "unknown"


def _lock_digest() -> str:
    lock = mf.REPO_ROOT / "uv.lock"
    if not lock.is_file():
        return "absent"
    return hashlib.sha256(lock.read_bytes()).hexdigest()[:16]


def hardware_block() -> str:
    blas_name, blas_version = _blas_info()
    lines = [
        "# Hardware and software",
        "",
        "| field | value |",
        "|---|---|",
        f"| CPU | {_cpu_model()} |",
        f"| physical cores | {psutil.cpu_count(logical=False)} |",
        f"| logical cores | {psutil.cpu_count(logical=True)} |",
        f"| RAM (GB) | {_total_ram_gb():.1f} |",
        f"| kernel | {platform.release()} |",
        f"| host | {_host_id()} (hostname digested, not disclosed) |",
        f"| Python | {platform.python_version()} |",
        f"| numpy | {np.__version__} |",
        f"| scipy | {_scipy_version()} |",
        f"| psutil | {psutil.__version__} |",
        f"| BLAS | {blas_name} {blas_version} |",
        f"| uv.lock sha256[:16] | {_lock_digest()} |",
        "",
        "Thread environment, pinned before numpy is imported:",
        "",
        "| variable | set to | prior value |",
        "|---|---|---|",
    ]
    for name in BLAS_VARS:
        lines.append(f"| {name} | 1 | {PRIOR_THREAD_ENV.get(name) or '(unset)'} |")
    lines += [
        "",
        "Every measurement ran single-threaded. Each cell ran in a fresh `spawn` process,",
        "so `ru_maxrss` is attributable to that cell rather than inherited from earlier work.",
    ]
    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #


def run(args) -> int:
    cells = enumerate_cells(args)
    if not cells:
        print("no cells selected")
        return 1

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    suffix = f"_{args.suffix}" if args.suffix else ""
    runs_path = OUTPUT_DIR / f"runtime_scaling{suffix}.csv"
    phases_path = OUTPUT_DIR / f"runtime_scaling_phases{suffix}.csv"

    total = len(cells) * (args.reps + 1)
    print(f"{len(cells)} cells x {args.reps + 1} repetitions ({total} runs), "
          f"timeout {args.timeout:g}s, soft cap {args.mem_cap_mb:g} MB")
    print(f"predicting dense OOM above the cap for scale cells\n")

    done = 0
    with runs_path.open("w", newline="") as runs_handle, phases_path.open(
        "w", newline=""
    ) as phases_handle:
        runs_writer = csv.DictWriter(runs_handle, fieldnames=RUN_FIELDS, extrasaction="ignore")
        phases_writer = csv.DictWriter(
            phases_handle, fieldnames=PHASE_FIELDS, extrasaction="ignore"
        )
        runs_writer.writeheader()
        phases_writer.writeheader()

        for cell in cells:
            footprint = predicted_dense_mb(cell)
            if (
                cell.solver == "dense"
                and footprint > args.mem_cap_mb
                and footprint > 0.5 * _total_ram_gb() * 1024
            ):
                print(f"  skip {cell.instance_id:32s} {cell.method:14s} dense "
                      f"(predicted {footprint:.0f} MB)")
                runs_writer.writerow(
                    {
                        **{key: "" for key in RUN_FIELDS},
                        **cell.as_dict(),
                        "isolation": "none",
                        "status": "predicted_oom",
                        "operator_footprint_mb": footprint,
                        "error_msg": f"analytic dense footprint {footprint:.0f} MB",
                    }
                )
                done += 1
                continue

            hard_limit = 0
            if cell.solver == "dense" and footprint > args.mem_cap_mb:
                hard_limit = int(args.mem_cap_mb * 1024 * 1024)

            for rep in range(-1, args.reps):  # rep == -1 is the untimed warmup
                record, phase_rows = _run_cell(cell, args.timeout, hard_limit)
                record["repeat_index"] = rep
                record["is_warmup"] = rep < 0
                runs_writer.writerow(record)
                for row in phase_rows:
                    row["repeat_index"] = rep
                    phases_writer.writerow(row)
                runs_handle.flush()
                phases_handle.flush()
                done += 1
                if rep < 0:
                    continue

            print(
                f"  {cell.instance_id[:32]:32s} {cell.method:14s} {cell.solver:12s} "
                f"{record['status']:13s} {record.get('wall_s') or 0:.3f}s "
                f"rss={record.get('peak_rss_mb') or 0:.1f}MB "
                f"gram={record.get('operator_footprint_mb') or 0:.1f}MB"
            )

    (OUTPUT_DIR / f"hardware{suffix}.md").write_text(hardware_block())
    print(f"\n{done} runs -> {runs_path}")
    print(f"           -> {phases_path}")
    print(f"           -> {OUTPUT_DIR / f'hardware{suffix}.md'}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--suites", default="real", help="real,synthetic,scale")
    parser.add_argument("--methods", default=",".join(METHODS))
    parser.add_argument("--specs", default="", help="scale specs, e.g. x1e3,x1e4")
    parser.add_argument("--datasets", default="", help="real datasets, e.g. cip_soc,lodes")
    parser.add_argument("--suffix", default="", help="suffix for the output filenames")
    parser.add_argument("--seeds", type=int, default=2, help="synthetic seeds per regime")
    parser.add_argument("--reps", type=int, default=3)
    parser.add_argument("--timeout", type=float, default=900.0, help="seconds per run")
    parser.add_argument(
        "--mem-cap-mb",
        type=float,
        default=8192.0,
        help="dense cells above this analytic footprint run under RLIMIT_AS and are "
        "expected to record status=oom",
    )
    args = parser.parse_args()
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())
