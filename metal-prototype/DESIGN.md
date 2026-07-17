# Metal acceleration design note

## What was studied

- AutoDock Vina 1.2.7 `parallel_mc`, `monte_carlo`, `quasi_newton`, `cache`, and `grid` paths.
- Vina-GPU 2.1's OpenCL `kernel1` (map construction) and `kernel2` (independent Monte Carlo lanes containing mutation and BFGS).
- Apple's Metal compute dispatch, runtime shader compilation, thread-grid sizing, and unified-memory model.

## Chosen boundary

Vina's CPU parallelism assigns one independent Monte Carlo trajectory to each CPU worker. Within every trajectory, mutation is followed by BFGS local optimization; BFGS repeatedly evaluates affinity grids and their derivatives. Existing GPU Vina work moves the complete trajectory onto the device so thousands of independent lanes can run concurrently.

Calling the GPU once for each scalar BFGS evaluation would be slower because command-buffer dispatch and synchronization would dominate the small amount of arithmetic. The useful first Metal primitive is therefore a **batched affinity-grid evaluator**: one GPU thread evaluates one pose, looping over its movable atoms. This establishes:

1. correct Vina map parsing and trilinear interpolation;
2. Float32 numerical tolerances against a CPU reference;
3. the batch size at which Metal amortizes dispatch overhead;
4. the buffer layouts needed by a later on-device BFGS/Monte Carlo kernel.

## Prototype scope

The prototype uses real 1IEP Vina affinity maps and real imatinib coordinates. It evaluates batches of translated query poses on both CPU and Metal and rejects any result outside the configured tolerance.

The derivative comparison reports a 99th-percentile error as well as a maximum. Trilinear interpolation energy is continuous at grid-cell boundaries, but its derivative is not; CPU and GPU Float32 rounding can select different valid one-sided derivatives at an exact boundary. Energy remains the primary strict invariant, while the percentile check catches systematic derivative errors without misclassifying isolated boundary choices.

The implementation now includes exact AD-to-XS typing, exact parser-generated internal pairs, Vina intramolecular potentials, affinity-grid derivatives, nested torsion transforms, Vina-style mutation and Metropolis acceptance, multiple-minimum retention, RMSD clustering, and authoritative Vina finalization. Parsing, map construction, clustering, and PDBQT writing remain on the CPU.

Local optimization uses inverse-BFGS over translation, orientation, and all torsions. Each lane accumulates grid and exact Vina intramolecular Cartesian derivatives, projects them into the ligand's degrees of freedom, updates a 14-by-14 maximum inverse Hessian, and follows Vina's ten-trial Armijo backtracking and first-step Hessian scaling without returning to the CPU.

The optimized kernel stores the symmetric inverse Hessian in Vina-style triangular form (105 rather than 196 floats at the maximum dimension) and reuses one 64-atom force buffer for local and world gradients. The host preserves Vina's total heuristic mutation count but selects a power-of-two lane count that leaves roughly 32 or more sequential mutations per trajectory. On Apple M3, 4,096 lanes is the measured quality/throughput point for the seven-torsion imatinib cases; 8,192 is only slightly faster but begins losing minima, while 16,384 degrades clearly.

## Next porting order

1. Recover performance through kernel specialization and occupancy improvements without changing the parity contract.
2. Generalize fixed limits (`64` scored atoms and `8` torsions) through generated kernel specializations.
3. Add macrocycle glue atoms, hydrated-docking pseudo-atoms, and flexible receptor support.
4. Precompile the Metal library when full Xcode is available to remove startup compilation time.
