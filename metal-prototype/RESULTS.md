# Metal prototype results

Tested on an Apple M3 MacBook Air with the real AutoDock Vina 1IEP affinity maps and 40 imatinib atoms. The kernel evaluates both trilinear grid energy and the analytic translation gradient. Timings exclude one-time runtime shader compilation and include command-buffer execution and synchronization.

| Poses | CPU | Metal | Speedup | Maximum energy error | P99 gradient error | Maximum gradient error |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 0.020 ms | 0.365 ms | 0.06× | 0.000002 | 0.000007 | 0.000007 |
| 64 | 0.802 ms | 0.356 ms | 2.25× | 0.000061 | 0.000019 | 0.000019 |
| 1,024 | 9.909 ms | 0.329 ms | 30.12× | 0.000061 | 0.000023 | 0.000038 |
| 16,384 | 152.422 ms | 0.846 ms | 180.16× | 0.000084 | 0.000023 | 1.043417 |
| 65,536 | 604.346 ms | 1.904 ms | 317.35× | 0.000111 | 0.000023 | 1.589083 |

Metal is counterproductive for individual evaluations but clearly useful once independent poses are batched. The rare maximum-gradient outliers occur at trilinear grid-cell boundaries, where energy is continuous but the analytic derivative has different valid one-sided values. The 99th-percentile derivative agreement remains approximately `2.5e-5` at scale.

These numbers isolate the translated Vina grid primitive. Complete docking measurements are reported below.

## End-to-end docking

All generated poses were rescored with the official Vina 1.2.7 executable.

| Experiment | Search time | Metal search energy | Official Vina intermolecular | Official estimated binding |
|---|---:|---:|---:|---:|
| Metal rigid, 8,192 × 512 | 0.353 s | -16.793 | -16.794 | -11.918 |
| Metal flexible, 8,192 × 512, hybrid refinement | 3.936 s | -16.643 | -16.057 | -11.394 |
| Metal flexible, 8,192 × 512, inverse-BFGS | 11.354 s | -17.551 | -17.036 (grid/no-refine) | -12.089 |
| Metal BFGS + Vina final local refinement | 12.344 s total | n/a | -18.223 | -12.931 |
| CPU Vina, 8 threads, exhaustiveness 8 | 19.248 s | n/a | approximately -18 | -13.210 |

The rigid result validates the translated grid path to 0.001 kcal/mol. The flexible search energy includes intramolecular energy; for its output Vina reported -16.057 intermolecular and -0.587 internal, matching the Metal combined energy to approximately 0.001.

The on-device inverse-BFGS result is internally consistent: its -17.551 Metal combined energy matches Vina's no-refine intermolecular plus internal energy (-17.036 - 0.512 = -17.548) within 0.003 kcal/mol. Vina's explicit-receptor local refinement then takes 0.99 seconds and improves the binding estimate to -12.931.

That hybrid production path is 1.56× faster than the earlier 8-thread CPU run while ending within 0.279 kcal/mol of its best reported binding estimate. This is a single-complex engineering benchmark, not evidence of general docking accuracy or throughput; broader receptor/ligand validation remains necessary.

## Broader validation

After porting Vina compatibility and optimizing per-lane storage and scheduling, the five-case three-repetition benchmark took 27.30 seconds for stock eight-thread Vina, 10.40 seconds for Metal search, and 12.92 seconds for Metal plus official Vina refinement. The end-to-end workflow is 2.11× faster. Four cases stayed within 0.008–0.157 kcal/mol of Vina; 1FPU was a 0.557 kcal/mol median stochastic outlier, although one seed reached -11.849 versus approximately -11.9 for Vina. See `../broader-benchmark/results/REPORT.md` for methods, per-case results, raw artifacts, and feature exclusions.
