# AutoDock Vina Metal

Vina-compatible Apple Metal docking engine plus reproducible benchmark harnesses.

This repository contains:

- `metal-prototype/`: the Metal backend that mirrors Vina scoring, search, clustering, and optional final Vina refinement
- `benchmark-100/`: a 100-docking benchmark built from five DUD-E Diverse targets
- `broader-benchmark/`: the earlier smaller compatibility benchmark
- `bin/vina`: the official AutoDock Vina 1.2.7 executable used for baseline runs and final refinement

The raw DUD-E Diverse archive is treated as an external input cache, not a versioned source artifact. The prepared benchmark inputs and generated reports remain in the repository.

The working name I would use for GitHub is `autodock-vina-metal`. It is explicit, searchable, and matches the upstream project naming better than a generic benchmark title.

## What This Experiment Measures

The benchmark compares:

- stock AutoDock Vina 1.2.7 on CPU
- Metal search plus authoritative Vina 1.2.7 explicit-receptor refinement

The final scores are always reported by official Vina. The Metal engine is responsible for the search phase and candidate generation.

## Build The Metal Engine

```sh
cd metal-prototype
swift build -c release
```

The release binary is written to:

```sh
metal-prototype/.build/release/VinaMetal
```

## Run A Single Dock

The Metal engine accepts the same receptor and map inputs used by Vina-compatible workflows:

```sh
metal-prototype/.build/release/VinaMetal \
  --dock \
  --flexible \
  --maps benchmark-100/prepared/hivrt/maps \
  --ligand benchmark-100/prepared/hivrt/ligands/01_CHEMBL262184.pdbqt \
  --output /tmp/metal-out.pdbqt \
  --vina-binary bin/vina \
  --vina-receptor benchmark-100/prepared/hivrt/receptor.pdbqt \
  --vina-config benchmark-100/prepared/hivrt/receptor.box.txt \
  --vina-output /tmp/metal-final.pdbqt
```

For a direct comparison case, use the benchmark runner below instead of hand-assembling inputs.

## Run The 100-Docking Benchmark

The benchmark uses five proteins and 20 ligands per target for a total of 100 cases.

```sh
bun benchmark-100/benchmark.ts
```

Useful flags:

```sh
bun benchmark-100/benchmark.ts --limit 10
bun benchmark-100/benchmark.ts --force
```

The runner is resumable. It writes each case result immediately under `benchmark-100/results/<case-id>/result.json`, then regenerates the aggregate artifacts:

- `benchmark-100/results/summary.json`
- `benchmark-100/results/results.csv`
- `benchmark-100/results/REPORT.md`

## Benchmark Protocol

The 100-case benchmark uses:

- DUD-E Diverse subset
- targets: `ampc`, `cxcr4`, `gcr`, `hivpr`, `hivrt`
- 20 ligands per target
- AutoDock Vina exhaustiveness `8`
- `8` CPU threads for the stock Vina baseline
- one output mode per run
- matching per-case seeds
- map preparation excluded from timing

The Metal engine adapts its lane count to ligand complexity before search and then runs official Vina refinement on the retained modes.

## Results

Benchmark version: `dude-diverse-100-v1`

Completed on July 17, 2026.

| Protein | Dockings | CPU Vina s | Metal + Vina s | Speedup | Mean abs score delta | Max abs score delta |
|---|---:|---:|---:|---:|---:|---:|
| AMPC | 20 | 62.47 | 34.32 | 1.82x | 0.102 | 0.577 |
| CXCR4 | 20 | 162.49 | 60.75 | 2.67x | 0.149 | 0.513 |
| GCR | 20 | 120.47 | 55.57 | 2.17x | 0.172 | 1.662 |
| HIVPR | 20 | 231.27 | 92.81 | 2.49x | 0.177 | 1.091 |
| HIVRT | 20 | 169.04 | 55.71 | 3.03x | 0.289 | 1.382 |

Total:

- 100 dockings
- CPU Vina: `745.74 s`
- Metal + Vina: `299.15 s`
- end-to-end speedup: `2.49x`
- Metal search only: `174.31 s`
- mean absolute final-score delta: `0.178 kcal/mol`
- maximum absolute final-score delta: `1.662 kcal/mol`

The aggregate report is in `benchmark-100/results/REPORT.md`.

Interpretation:

- 57/100 cases stayed within `0.1 kcal/mol`
- 78/100 cases stayed within `0.2 kcal/mol`
- 22 cases scored better than stock Vina, 77 scored worse, and 1 matched exactly
- the worst outliers also had large pose RMSD, which means they are driven by search-path divergence, not just float noise

## Release Status

The engine release notes are in `RELEASE.md`.

Publishing the GitHub release still requires a valid authenticated GitHub session. The local CLI session on this machine is currently invalid, so the repository can be committed locally now but cannot be pushed or released until authentication is repaired.
