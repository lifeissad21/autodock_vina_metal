# dude-diverse-100-v1 results

Completed 2026-07-17T17:07:55.956Z on Apple M3. Stock AutoDock Vina 1.2.7 used 8 CPU threads; the hybrid used Metal search followed by authoritative Vina 1.2.7 explicit-receptor refinement. Both used exhaustiveness 8, one output mode, and matching per-case seeds. Map preparation is excluded.

| Protein | Dockings | CPU Vina s | Metal + Vina s | Speedup | Mean abs score delta | Max abs score delta |
|---|---:|---:|---:|---:|---:|---:|
| AMPC | 20 | 62.47 | 34.32 | 1.82x | 0.102 | 0.577 |
| CXCR4 | 20 | 162.49 | 60.75 | 2.67x | 0.149 | 0.513 |
| GCR | 20 | 120.47 | 55.57 | 2.17x | 0.172 | 1.662 |
| HIVPR | 20 | 231.27 | 92.81 | 2.49x | 0.177 | 1.091 |
| HIVRT | 20 | 169.04 | 55.71 | 3.03x | 0.289 | 1.382 |

**Total:** 100 dockings; CPU Vina 745.74 s; Metal + Vina 299.15 s; 2.49x end-to-end speedup. Metal kernels totaled 174.31 s. Mean absolute final-score delta was 0.178 kcal/mol; maximum was 1.662 kcal/mol.

Score delta is hybrid minus stock Vina. Both final scores are calculated by official Vina 1.2.7, but independent stochastic searches can select different minima. Pose RMSD in the CSV is direct atom-order heavy-atom RMSD and is not symmetry corrected.
