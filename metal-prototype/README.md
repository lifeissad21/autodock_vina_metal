# Vina Metal prototype

Experimental Apple Metal docking engine derived from AutoDock Vina's scoring and search design. It supports exact Vina atom typing, real affinity maps, rigid motion, up to eight ligand torsions, intramolecular Vina potentials, parallel Monte Carlo lanes, local refinement, and PDBQT output. Read [DESIGN.md](DESIGN.md) and [RESULTS.md](RESULTS.md) for scope and measured results.

## Build and run

```sh
cd /Users/gpm/Desktop/Repos/vina-multicore-benchmark/metal-prototype
./tools/build-extractor.sh
swift build -c release
.build/release/VinaMetal
```

Run docking:

```sh
# Rigid-body search with Vina-derived adaptive effort
.build/release/VinaMetal --dock --exhaustiveness 8

# Flexible search, clustering, and authoritative Vina finalization
.build/release/VinaMetal --dock --flexible --exhaustiveness 8 --num-modes 9 \
  --vina-receptor ../data/1iep_receptor.pdbqt \
  --vina-config ../data/1iep_receptor.box.txt
```

The Metal shader is compiled from source once at process startup because only the Command Line Tools are installed on this machine; the standalone `metal` compiler from full Xcode is unavailable. Runtime compilation is excluded from the kernel timings.

The helper under `tools/` links against the version-matched Vina 1.2.7 source to export the exact ligand atom types and internal interaction pairs. This avoids subtly different topology reconstruction in the Metal host.

The flexible optimizer performs inverse-BFGS entirely inside each Metal lane. It differentiates grid and intramolecular energies, projects atomic forces onto translation, orientation, and torsion coordinates, and follows Vina's Armijo line search and initial Hessian scaling. Global/local effort, mutation selection, torsion replacement, rotation amplitude, Metropolis temperature, and RMSD clustering follow Vina 1.2.7 semantics. Total work is distributed across shorter GPU lanes, so stochastic trajectories are not identical to Vina's eight CPU tasks.

When `--vina-receptor` and `--vina-config` are supplied, all retained modes are refined concurrently by the version-matched official Vina executable. The best authoritative double-precision score and pose are written separately from the raw Metal search and clustered-mode outputs.

For maximum single-result throughput use `--num-modes 1`. The adaptive scheduler chooses the lane count automatically. Explicit `--lanes`, `--steps`, and `--local-steps` remain available for controlled experiments, but overriding them can change search quality.
