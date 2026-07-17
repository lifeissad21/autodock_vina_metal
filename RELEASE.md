# AutoDock Vina Metal v0.1.0

Release tag: `v0.1.0`

Release name: `autodock-vina-metal`

This is the first published release of the Metal-backed Vina-compatible docking engine and its benchmark suite.

## Included Artifacts

- Source code archive from the `v0.1.0` tag
- macOS arm64 binary bundle for `VinaMetal`
- benchmark inputs and reports under `benchmark-100/`

The raw DUD-E download archive is treated as a local input cache and is not part of the published release payload.

## Downloadable Assets

- `VinaMetal-macos-arm64.zip`
- GitHub-generated source archive for `v0.1.0`

## What Changed

- Added an Apple Metal docking engine that keeps Vina-compatible search, scoring, clustering, and final refinement semantics
- Added a resumable 100-docking benchmark spanning five DUD-E Diverse targets
- Added benchmark reports and per-case artifacts for reproducibility
- Added release packaging guidance for the macOS arm64 binary

## Validation

The 100-docking benchmark completed with:

- CPU Vina: `745.74 s`
- Metal + Vina: `299.15 s`
- end-to-end speedup: `2.49x`
- mean absolute final-score delta: `0.178 kcal/mol`
- maximum absolute final-score delta: `1.662 kcal/mol`

See `benchmark-100/results/REPORT.md` for the full table and per-target breakdown. The benchmark result files remain in the repository under `benchmark-100/results/`.

## Build

```sh
cd metal-prototype
swift build -c release
```

## Run

```sh
metal-prototype/.build/arm64-apple-macosx/release/VinaMetal \
  --dock \
  --flexible \
  --maps benchmark-100/prepared/hivrt/maps \
  --ligand benchmark-100/prepared/hivrt/ligands/01_CHEMBL262184.pdbqt \
  --vina-binary bin/vina \
  --vina-receptor benchmark-100/prepared/hivrt/receptor.pdbqt \
  --vina-config benchmark-100/prepared/hivrt/receptor.box.txt
```

## Benchmark

```sh
bun benchmark-100/benchmark.ts
```

The benchmark runner is resumable and writes results to `benchmark-100/results/`.
