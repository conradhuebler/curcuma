# EEQ solver control (GFN-FF)

Claude Generated (Oct 2026) — 🤖 AI-generated, ⚙️ machine-tested (`cli_gfnff_06_eeq_solver_control`, polymer_2x, the 1410-atom chain and caffeine on CPU and one RTX 5080); human production testing pending.

## Rule

`-eeq_solver.solve_method NAME` (alias `-gfnff.solve NAME`; the `eeq_solver` spelling wins when both are set to a non-auto name) chooses the
linear solve of the EEQ charge system. **The same names run on the CPU and, with `-gpu cuda`, on the GPU.** An explicit name is never
overridden by a size threshold.

| name | algorithm | CPU | GPU (CUDA) |
|---|---|---|---|
| `auto` (default) | `ppcg` if N >= `eeq_ppcg_min_atoms` (500) and nfrag >= `eeq_ppcg_min_nfrag` (1), else `cholesky`. Never picks an inexact method | yes | yes |
| `cholesky` | dense Cholesky of A_nn + Schur complement, at any size | yes | WP5-A (one fragment) / WP7-A |
| `ppcg` | projected PCG, one iterative solve for any nfrag | yes | WP7-E |
| `pcg` | per-fragment PCG with block-Jacobi, nfrag+1 solves; stops at a relative residual of 1e-6 (about 4e-7 Eh on a 1410-atom chain); cost grows with nfrag, see below | yes | WP7-C |
| `batched` | per-fragment solves, **drops cross-fragment Coulomb**: inexact; warns when the estimated error exceeds `eeq_batched_max_error_eh` (below) | yes | WP7-B (one fragment: same as `cholesky`) |
| `lu` | dense LU of the augmented (N+nfrag) system | yes | augmented LU, cuSOLVER getrf/getrs (Oct 2026; ROCm/Vulkan: error) |
| `ldlt` | dense Bunch-Kaufman LDL^T of A_nn + Schur (`dsytrf` + `dsytrs2`) | yes | cuSOLVER `sytrf` + `Xsytrs` (in blocks of 1024 right-hand sides) + Schur (Oct 2026; ROCm/Vulkan: error) |

Aliases `schur_cholesky` = `cholesky`, `projected_pcg` = `ppcg`. **An unknown name is an error** that lists the valid ones (it used to mean `auto`).
The thresholds are arguments: `-eeq_solver.eeq_ppcg_min_atoms`, `-eeq_solver.eeq_ppcg_min_nfrag` (0 = auto never picks `ppcg`).

## One implementation of the choice

`EEQSolver::resolveMethod()` (`eeq_solver.cpp`) is the only place that turns a requested name into the method that runs; the CPU dispatch and
the GPU wrapper (`gfnff_gpu_method_impl.h`) both call it. `fixEEQSolverChoice()` (`method_factory.cpp`) reads the CLI scopes once, rejects typos and
GPU gaps, and writes the canonical name back so CPU and GPU see the same value.

## Output (verbosity 1)

One line per distinct choice and run: `EEQ solver (<phase>): <name>[ (auto)] [CPU | GPU CUDA, <path>], N=.., nfrag=..`. The CPU prints
`topology charges` and `final charges`; the GPU prints `final charges` with its path (`WP7-E projected PCG`, `WP5-A Cholesky`, ...) at the first
device solve. A single point with `-gpu cuda` still solves its charges on the CPU (setup), so its line says `[CPU]`; the GPU path appears in MD/opt.

## What changed against the old behaviour

- Default was `cholesky` and a hidden switch turned it into `ppcg` at 500 atoms; an explicit `cholesky` was switched too. Now `auto` is the default and `cholesky` is dense.
- The CPU `auto` benchmark (timing Cholesky against PCG at the first call) and its Gershgorin bypass are removed; `auto` never selects `batched` any more.
- `-gfnff.solve` was ignored on the GPU; the GPU had its own `Auto` branch; unknown names fell back to `auto` silently.
- GPU stalls of `ppcg`/`pcg`/`batched` (fall-through to the dense solve) are counted GPU fallbacks: `-gpu_strict true` ends the run with exit code 3.

## batched: energy criterion (`eeq_batched_max_error_eh`, default 1e-6 Eh)

The batched charges q_b leave the residual r = A q_b - x of the full EEQ system. Its energy excess is
1/2 r^T (P A P)^+ r; with A replaced by its fragment blocks this becomes
dE ~ 1/2 sum_f r_f^T d_f, [A_ff 1; 1^T 0][d_f; mu] = [r_f; 0] - one O(N^2) product plus the fragment blocks, made once per
topology and phase (CPU) or once per topology (GPU path, on the host), only at verbosity >= 1. Above the threshold a warning gives
the estimate. It replaces `eeq_batched_min_distance` (15 Bohr) and the unused `eeq_contact_prefer_exact` (both removed; the
per-step O(N^2) contact scan that fed them is gone). Checked against single points (dE measured = E(batched) - E(cholesky)):

| system | points | measured / estimate |
|---|---:|---|
| water dimer, dmin 2-15 A | 13 (20/25 A below print precision) | 0.95-1.11, at 2.0 A 3.2 |
| acetic-acid dimer, dmin 2-20 A | 14 | 0.97-1.12 (20 A: 0.86, dE = 2e-8 has one digit) |
| H3O+(H2O)2, charge +1, dmin 2-25 A | 15 | 1.42-1.99 (the estimate is low) |
| polymer_2x (7320 atoms, 1502 fragments) | 1 | 0.09502 / 0.0956 Eh = 0.99 |

Cost on polymer_2x: 82 ms on the CPU, 336 ms on the GPU path (host loop), once (n = 1 each).

## Measurements (Oct 2026, polymer_2x = 7320 atoms / 1502 fragments, 1 thread, RTX 5080)

🤖 Claude Generated, ⚙️ machine-tested; one run per configuration, MD per step median over 49 steps (pcg: 2 steps).

| solver | CPU ms/step | CUDA ms/step | Epot after 50 fs / Eh |
|---|---:|---:|---:|
| auto = ppcg | 2325 | 66 | -902.161126 |
| cholesky | 5417 | 610 | -902.161126 |
| lu | 7290 | 641 | -902.161126 |
| ldlt | 5499 | 1039 | -902.161126 |
| batched (inexact) | 1609 | 150 | -901.887405 |
| pcg (3 fs) | 48672 | 36387 | -900.545812 / -900.545813 |

- **batched on CUDA was wrong from the second device solve on** (fixed Oct 9): `dpotrs` overwrote the all-ones constraint column,
  so every later solve used the previous step's `Z2` as right-hand side (+840 Eh at step 2).
- **lu on CUDA**: augmented (N+nfrag) LU, charges = CPU lu to 7e-12 e; device matrix (N+nfrag)^2, 623 MB here.
- **ldlt on CUDA**: `cusolverDnDsytrf` + `cusolverDnXsytrs` on A_nn, then the WP7-A Schur step; charges = CPU ldlt to 3e-14 e.
  Memory as `cholesky` (N^2 matrix, 428 MB, + N x (nfrag+1) right-hand sides, 88 MB). **cuSOLVER 12.2 `Xsytrs` returns wrong
  columns with info = 0 for more than 1024 right-hand sides** (standalone check, n = 500-7320), hence the blocks of 1024.
- **pcg** is slow by construction: nfrag+1 = 1503 right-hand sides. CPU: 21 block-PCG iterations, each a 7320 x 7320 x 1503 GEMM
  (2 s); the block-Jacobi preconditioner now solves all columns per fragment at once (13.0 -> 4.1 s per solve, charges move by
  up to 4e-9 e). GPU: 1503 sequential PCG solves with diagonal Jacobi, 74-90 iterations each, 0.32 ms per iteration.
- **ldlt on the CPU** uses `dsytrs2` (level 3) instead of `dsytrs`: 7.3 -> 5.5 s per step, charges move by 2e-13 e.
  **lu on the CPU** is at the cost of the LU itself (Eigen and MKL `dgetrf` both 5.9 s for the 8822-row matrix).

## Local-inverse preconditioner of ppcg (`eeq_ppcg_precond localinv`, opt-in)

Claude Generated (Oct 10, 2026) — 🤖 AI-generated, ⚙️ machine-tested (`cli_gfnff_07_eeq_localinv`, measurements below); human production testing pending.

`-eeq_solver.eeq_ppcg_precond localinv` (default `diag`) replaces the diagonal Jacobi of the projected PCG, CPU and `-gpu cuda` (WP7-E).
Per atom the dense block of A over all atoms within `eeq_localinv_ext_bohr` (default 10) is factorised, the column of the atom is solved and
its entries within `eeq_localinv_radius_bohr` (default 8) are kept; G = (Gl + Gl^T)/2. It enters as a constraint preconditioner,
z = G r - G C^T S^-1 C G r with S = C G C^T, so the iterates stay on the fragment-charge constraints; derivation and references in
`eeq_solver.h` (`LocalInversePreconditioner`). Only the Phase-2 (final) charges use it; Phase 1 has no geometry and stays on `diag`.

- **Definiteness**: G is accepted only if a sparse Cholesky of G and of S succeeds; otherwise, or when the CG shows rz <= 0 / pAp <= 0
  (GPU: also a stall), the solve runs on `diag` with a counted warning. No rejection occurred on the systems below (0 of all builds).
  `CURCUMA_LOCALINV_TEST_REJECT=1` forces a rejection (test hook).
- **Rebuild**: `eeq_localinv_rebuild_bohr` (default 1.0; 0 = every solve) rebuilds when any atom moved further since the last build. A stale
  G changes only the iteration count, never the converged charges (same `eeq_ppcg_tol`).
- **Reproducibility**: the build is independent of the thread count (one thread per atom row, fixed summation order); MD charges at a fixed
  thread count are bit-identical between runs. GPU: CSR kernels + dense S^-1 on the device (nfrag <= 4096), built on the host; not bit-reproducible, like the GPU path as a whole.
- **CUDA only** on the device side; ROCm keeps `diag`.

Measurements (Oct 10, 2026; MD 50 fs, dt 1 fs, 300 K; per-step means over steps 1-49, n = 3 runs unless noted, spread = range of the run means;
RTX 5080, 32 cores). Epot after 50 fs -902.161126 Eh in every run; per-step charges within 3.5e-11 e (CPU) / 3.4e-11 e (GPU) of `diag`.

| polymer_2x (7320 atoms, 1502 fragments) | iterations/step | EEQ ms/step | step ms |
|---|---:|---:|---:|
| CPU 1 thread, diag | 68 | 1208 [1208-1208] | 2339 [2338-2340] |
| CPU 1 thread, localinv rebuild 0 | 15 | 599 [598-599] | 1726 [1723-1730] |
| CPU 32 threads, diag | 68 | 1004 [999-1008] | 1252 [1244-1259] |
| CPU 32 threads, localinv rebuild 0 | 15 | 411 [407-417] | 657 [651-665] |
| CUDA, diag | 68 | 32.2 [32.2-32.2] | 73.6 [73.5-73.7] |
| CUDA, localinv rebuild 0 | 15 | 51.0 [50.9-51.1] | 92.4 [92.2-92.7] |
| CUDA, localinv rebuild 0.5 (13 builds) | 16 | 25.4 [25.3-25.6] | 66.9 [66.7-67.2] |
| CUDA, localinv rebuild 1.0 (7 builds) | 16 | 21.3 [21.2-21.3] | 62.7 [62.5-62.9] |

- Build on polymer_2x (R 8 / R_ext 10): 48 ms at 1 thread, 22 ms at 32 threads, of which the serial definiteness check (sparse Cholesky of G)
  is 17 ms; GPU path 31 ms (adds the dense S^-1, 1502^2). Peak RSS +6 MB over `diag` (3277 MB, SP, 32 threads); R 12/16: +20 MB, R 16/20: +54 MB.
- Other systems, 16 threads, EEQ ms/step diag -> localinv (rebuild 0): peo201-water500 104 -> 54; NaCl in water (peo201-water500 with 40 waters
  replaced by 20 Na + 20 Cl, neutral) 96 -> 51; the same with 24 Na / 20 Cl, `-charge 4 -gfnff.frag_charge_model reference`, n = 1: 95 -> 49;
  urea400-water1000 440 -> 362 (build 166 ms, check-dominated); peo201-chain (nfrag 1) 13.4 -> 13.0.
- Radius scan (SP, 1 thread, Phase-2 iterations; diag 72 / 62 / 63 / 37 for polymer_2x / chain / water500 / urea): R 6: 21-22 / 19 / 19 / 11;
  R 8: 16 / 14 / 14-15 / 9; R 10: 13 / 12 / 12 / 8; R 12: 12 / 10 / 11 / 7. R_ext above R + 2 Bohr lowers no count but costs the cube of the block.
- The Python reconstruction that motivated this (vault note "curcuma Schwarz-Attention") predicted 81 -> 18 on polymer_2x; the C++ measures 68 -> 15.

## Threaded matvec of the CPU ppcg (Oct 10, 2026)

The dense products `A_nn * p` (once per iteration) and `A_nn * q` (residual start) run over fixed blocks of `kPpcgMatvecBlockRows` = 64 rows
(`solveWithProjectedPCG`, `eeq_solver.cpp`), distributed over `-threads` like `buildLocalInverse` (`m_ppcg_pool`, `m_ppcg_threads`). Every output row is
computed by one thread, the blocks do not depend on the thread count, there is no reduction: the product is bitwise identical for every `-threads` and
equal to the former unblocked `A_nn * x` (checked once with a temporary debug hook on polymer_2x, 4 and 32 threads, block sizes 64 and 256; the hook is removed).
`CURCUMA_PPCG_MATVEC_BLOCK=<rows>` overrides the block size (measurement only). The product is memory-bandwidth bound (428 MB per product on polymer_2x):
`-sp`, Phase-2 solve 902 -> 715 (4 threads), 710 (16), 754 ms (32 threads = SMT, slower than 16); 1 thread 904 -> 901 ms (n = 3 each).
Machine ceiling (OpenMP read of the same 428 MB, 1/2/16 threads): 21 / 43 / 60 GB/s; the matvec reaches 35 GB/s on one thread and 45 GB/s threaded, so
2x over serial is not reachable here. Open: exploiting the symmetry of `A_nn` (half the bytes). The matrix itself differs by ~1e-13 between thread counts
(threaded assembly, also before this change), so charges are not bitwise equal across `-threads`; energy and iteration counts are (test 07).

## Not tested / open

- ROCm and Vulkan: same code path, not built or run here. `lu` and `ldlt` are errors there.
- Multi-GPU split and the batch-worker device pool were not run after this change.
- localinv: the GPU-side breakdown fallback (rz <= 0 in the device CG) was never triggered (only the host rejection, via the test hook); not run
  with ALPB solvation, `eeq_distance_cutoff > 0`, transition metals, MD longer than 50 fs, optimisations, or nfrag > 4096 on the GPU (falls back to diag).

## Test note (Oct 10, 2026)

All tests that use the GPU carry `RESOURCE_LOCK "gpu"`, so `ctest -j8` runs them one after another. Without it `gfnff_gpu_vs_cpu_*` failed with
"Insufficient GPU memory" when several GPU tests overlapped on one 16 GB card: 4 of 7 `-j8 -L gfnff` runs on the state before the Oct 9 solver changes
(worktree at `a2b7c5f0`, built with the same flags) and 3 of 7 after them; with the lock 0 of 7 (87 tests each). Serial runs were green throughout.
