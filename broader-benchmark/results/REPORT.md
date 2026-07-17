# AutoDock Vina 1.2.7 versus Metal rewrite

Apple M3, 3 repetitions per case, median reported. CPU Vina used 8 threads, exhaustiveness 8, one mode, and fixed seeds. Metal used Vina's ligand-complexity formulas for total Monte Carlo work and local BFGS steps, distributed across an adaptive power-of-two lane count that preserves at least about 32 sequential mutations per lane. Hybrid time includes official Vina explicit-receptor local refinement of the best Metal pose. Map generation is excluded.

| Case | CPU Vina s | Metal search s | Metal + refine s | Hybrid speedup | CPU score | Refined Metal score | Score delta | Pose RMSD Å |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 1iep | 8.68 | 2.89 | 3.40 | 2.55× | -13.230 | -13.222 | 0.008 | 0.86 |
| 1fpu_rigid | 9.99 | 2.94 | 3.44 | 2.90× | -11.910 | -11.353 | 0.557 | 1.00 |
| 1s63 | 4.83 | 1.89 | 2.42 | 1.99× | -9.011 | -8.929 | 0.082 | 0.17 |
| 5x72_p59 | 1.86 | 1.34 | 1.83 | 1.02× | -10.980 | -10.904 | 0.076 | 0.09 |
| 5x72_p69 | 1.94 | 1.33 | 1.82 | 1.06× | -10.470 | -10.313 | 0.157 | 0.07 |

Aggregate median-case time: CPU 27.30 s, Metal search 10.40 s, Metal plus refinement 12.92 s. Aggregate hybrid speedup: 2.11×.

Score delta is refined Metal minus CPU Vina; positive values are worse. RMSD is a direct atom-order heavy-atom RMSD between independently selected best poses and is not symmetry corrected, so it diagnoses different minima but is not a formal pose-accuracy metric.

## Why results differ

The engines now share Vina's ligand-complexity effort formula, one-entity mutation distribution, torsion randomization, gyration-scaled rotation, Armijo line search, inverse-BFGS initialization, RMSD clustering rule, and authoritative final score/refinement. They still do not execute identical stochastic trajectories: Metal distributes the same total mutation count across an adaptive set of shorter parallel lanes, uses a GPU random-number stream, and evaluates in Float32. Stock Vina uses eight longer CPU trajectories, Boost MT19937/distributions, and double precision. Those differences can alter Metropolis and line-search decisions. Raw Metal energy intentionally remains a fast search objective; reported final scores come from official Vina 1.2.7 and include torsional, unbound-state, and explicit-receptor corrections.

## Interpretation

Triangular Hessian storage, one atom-gradient buffer, and occupancy-aware adaptive lanes reduced aggregate Metal search time to about 38% of stock Vina time. Including authoritative refinement, the workflow is 2.11× faster overall. Four cases retain score differences of 0.008–0.157 kcal/mol and direct pose RMSDs of 0.07–0.86 Å; 1FPU remains a stochastic outlier at 0.557 kcal/mol and 1.00 Å for the median seed, although one of its three Metal runs reached -11.849 versus Vina's approximately -11.9.

Unsupported official examples were excluded: BACE-1 macrocycle (22 branches and glue-atom treatment), hydrated 1UW6 (AD4 water pseudo-atoms), and flexible receptor sidechains.
