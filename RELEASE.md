# Vina Metal Release Notes

Proposed release name: `autodock-vina-metal`.

Proposed release tag: `v0.1.0`.

This release packages the Metal-backed Vina-compatible docking engine and the benchmark artifacts that validate it:

- `metal-prototype/`: Metal search engine and explicit Vina refinement path
- `benchmark-100/`: 100-case DUD-E Diverse benchmark
- `broader-benchmark/`: smaller compatibility benchmark used during development

The raw DUD-E Diverse download archive is treated as an external input cache. The prepared benchmark inputs and the generated result reports are the versioned experiment artifacts.

## Release Highlights

- Vina-compatible scoring and refinement through the official `bin/vina` executable
- adaptive GPU lane scheduling based on ligand complexity
- clustered mode retention with one-mode benchmark output
- resumable benchmark runner that writes per-case JSON immediately
- documented 100-docking benchmark on five proteins

## Validation Summary

The 100-docking benchmark completed successfully with:

- CPU Vina: `745.74 s`
- Metal + Vina: `299.15 s`
- end-to-end speedup: `2.49x`
- mean absolute final-score delta: `0.178 kcal/mol`
- maximum absolute final-score delta: `1.662 kcal/mol`

The aggregate report is available at `benchmark-100/results/REPORT.md`.

## Release Notes For GitHub

Use this file as the body for the GitHub release once the repository is pushed and the GitHub session is authenticated again.

Suggested artifact list:

- source tree at the release tag
- `benchmark-100/results/summary.json`
- `benchmark-100/results/results.csv`
- `benchmark-100/results/REPORT.md`

## Build

```sh
cd metal-prototype
swift build -c release
```

## Docking

```sh
metal-prototype/.build/release/VinaMetal \
  --dock \
  --flexible \
  --maps benchmark-100/prepared/hivrt/maps \
  --ligand benchmark-100/prepared/hivrt/ligands/01_CHEMBL262184.pdbqt \
  --vina-binary bin/vina \
  --vina-receptor benchmark-100/prepared/hivrt/receptor.pdbqt \
  --vina-config benchmark-100/prepared/hivrt/receptor.box.txt
```

The benchmark runner in `benchmark-100/benchmark.ts` is the preferred way to reproduce the full experiment.
