# Broader Vina versus Metal benchmark

This benchmark uses prepared complexes distributed with the official AutoDock Vina 1.2.7 source tree. It compares the stock Apple ARM64 Vina 1.2.7 executable with the experimental Metal rewrite using fixed seeds and one reported mode.

Run from the repository root:

```sh
bun broader-benchmark/benchmark.ts
```

The runner generates version-matched Vina affinity maps, performs three repetitions per engine and case by default, rescores and locally refines Metal poses with official Vina, and writes raw logs plus `results.json`, `results.csv`, and `REPORT.md` under `broader-benchmark/results/`.

The comparison is deliberately limited to the Vina scoring function and ligands supported by the current Metal kernel: at most 64 heavy atoms, at most eight ordinary torsions, and no macrocycle glue atoms, hydrated-docking water pseudo-atoms, flexible receptor, or AD4-specific zinc pseudo-atoms.
