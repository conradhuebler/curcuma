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

## Not tested / open

- ROCm and Vulkan: same code path, not built or run here. `lu` and `ldlt` are errors there.
- Multi-GPU split and the batch-worker device pool were not run after this change.

## Test note (Oct 10, 2026)

All tests that use the GPU carry `RESOURCE_LOCK "gpu"`, so `ctest -j8` runs them one after another. Without it `gfnff_gpu_vs_cpu_*` failed with
"Insufficient GPU memory" when several GPU tests overlapped on one 16 GB card: 4 of 7 `-j8 -L gfnff` runs on the state before the Oct 9 solver changes
(worktree at `a2b7c5f0`, built with the same flags) and 3 of 7 after them; with the lock 0 of 7 (87 tests each). Serial runs were green throughout.
