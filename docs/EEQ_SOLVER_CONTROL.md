# EEQ solver control (GFN-FF)

Claude Generated (Oct 2026) — 🤖 AI-generated, ⚙️ machine-tested (`cli_gfnff_06_eeq_solver_control`, polymer_2x and caffeine on CPU and one RTX 5080); human production testing pending.

## Rule

`-eeq_solver.solve_method NAME` (alias `-gfnff.solve NAME`; the `eeq_solver` spelling wins when both are set to a non-auto name) chooses the
linear solve of the EEQ charge system. **The same names run on the CPU and, with `-gpu cuda`, on the GPU.** An explicit name is never
overridden by a size threshold.

| name | algorithm | CPU | GPU (CUDA) |
|---|---|---|---|
| `auto` (default) | `ppcg` if N >= `eeq_ppcg_min_atoms` (500) and nfrag >= `eeq_ppcg_min_nfrag` (1), else `cholesky`. Never picks an inexact method | yes | yes |
| `cholesky` | dense Cholesky of A_nn + Schur complement, at any size | yes | WP5-A (one fragment) / WP7-A |
| `ppcg` | projected PCG, one iterative solve for any nfrag | yes | WP7-E |
| `pcg` | per-fragment PCG with block-Jacobi, nfrag+1 solves; stops at a relative residual of 1e-6 (about 4e-7 Eh on a 1410-atom chain) | yes | WP7-C |
| `batched` | per-fragment solves, **drops cross-fragment Coulomb**: inexact unless the fragments are well separated | yes | WP7-B (one fragment: same as `cholesky`) |
| `lu` | dense LU of the augmented system | yes | **not implemented: error with `-gpu`** |
| `ldlt` | dense Bunch-Kaufman LDL^T of A_nn + Schur | yes | **not implemented: error with `-gpu`** |

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

## Not tested / open

- ROCm and Vulkan: same code path, not built or run here. `lu` and `ldlt` on the GPU are errors until implemented.
- `batched` on CUDA reports "batched failed, falling through to WP7-A" on polymer_2x (now a counted fallback); cause not investigated.
- Multi-GPU split and the batch-worker device pool were not run after this change.
