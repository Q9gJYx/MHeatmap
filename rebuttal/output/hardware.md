# Hardware and software

| field | value |
|---|---|
| CPU | AMD Ryzen 9 7940H w/ Radeon 780M Graphics |
| physical cores | 8 |
| logical cores | 16 |
| RAM (GB) | 30.1 |
| kernel | 7.0.0-34-generic |
| host | machine-3a00f8af (hostname digested, not disclosed) |
| Python | 3.13.13 |
| numpy | 2.4.4 |
| scipy | 1.17.1 |
| psutil | 7.2.2 |
| BLAS | scipy-openblas 0.3.31.188.0 |
| uv.lock sha256[:16] | 79321e4d42418884 |

Thread environment, pinned before numpy is imported:

| variable | set to | prior value |
|---|---|---|
| OMP_NUM_THREADS | 1 | (unset) |
| OPENBLAS_NUM_THREADS | 1 | (unset) |
| MKL_NUM_THREADS | 1 | (unset) |
| NUMEXPR_NUM_THREADS | 1 | (unset) |

Every measurement ran single-threaded. Each cell ran in a fresh `spawn` process,
so `ru_maxrss` is attributable to that cell rather than inherited from earlier work.
