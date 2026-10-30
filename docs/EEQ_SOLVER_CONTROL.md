# EEQ solver control (GFN-FF)

Claude Generated (Oct 2026) — 🤖 AI-generated, ⚙️ machine-tested (`cli_gfnff_06_eeq_solver_control`, polymer_2x, the 1410-atom chain and caffeine on CPU and one RTX 5080); human production testing pending.

## Rule

`-gfnff.solve_method NAME` (flat `-solve_method NAME`; the older spellings `-eeq_solver.solve_method NAME` and `-gfnff.solve NAME` are still
accepted, see "Parameters" below for which wins) chooses the linear solve of the EEQ charge system. **The same names run on the CPU and, with `-gpu cuda`, on the GPU; `-gpu rocm` accepts `auto`, `cholesky` and `ppcg` (Oct 2026), `-gpu vulkan` accepts all names (GFN-FF runs on the CPU there).** An explicit name is never
overridden by a size threshold.

| name | algorithm | CPU | GPU (CUDA) |
|---|---|---|---|
| `auto` (default) | `ppcg` if N >= `eeq_ppcg_min_atoms` (500) and nfrag >= `eeq_ppcg_min_nfrag` (1), else `cholesky`. Never picks an inexact method | yes | yes |
| `cholesky` | dense Cholesky of A_nn + Schur complement, at any size | yes | WP5-A (one fragment) / WP7-A |
| `ppcg` | projected PCG, one iterative solve for any nfrag | yes | WP7-E |
| `pcg` | per-fragment PCG with block-Jacobi, nfrag+1 solves; stops at a relative residual of 1e-6 (about 4e-7 Eh on a 1410-atom chain); cost grows with nfrag, see below | yes | WP7-C |
| `batched` | per-fragment solves, **drops cross-fragment Coulomb**: inexact; warns when the estimated error exceeds `eeq_batched_max_error_eh` (below) | yes | WP7-B (one fragment: same as `cholesky`) |
| `lu` | dense LU of the augmented (N+nfrag) system | yes | augmented LU, cuSOLVER getrf/getrs (Oct 2026; ROCm: error; Vulkan: CPU solver) |
| `ldlt` | dense Bunch-Kaufman LDL^T of A_nn + Schur (`dsytrf` + `dsytrs2`) | yes | cuSOLVER `sytrf` + `Xsytrs` (in blocks of 1024 right-hand sides) + Schur (Oct 2026; ROCm: error; Vulkan: CPU solver) |

ROCm has a device `ppcg` (Oct 2026, `has_device_ppcg`; diagonal Jacobi or, with the default `eeq_ppcg_precond localinv`, the local inverse built on the device (default, `has_localinv_device_build`) or, with `CURCUMA_LOCALINV_HOST_BUILD=1`, built on the host and uploaded as CSR (`has_localinv_upload`)) but no device `pcg`/`batched` (`has_device_schur=false`): `pcg`, `batched`, `lu`, `ldlt` with `-gpu rocm` are an error naming the valid names (`auto, cholesky, ppcg`). `cholesky` runs the dense device Cholesky of the WP2 path plus the host Schur step (verbosity 1: `cholesky ... [GPU ROCm, WP2 dense Cholesky + CPU Schur]`). `auto` and `ppcg` run the device projected PCG (`... WP7-E projected PCG, diagonal Jacobi`); with 16 or more fragments a stalled device solve uses the CPU solver for that step, and giving `-gfnff.eeq_rocm_cpu_fragment_threshold` explicitly forces the CPU solver from that fragment count on. Names checked by `eeq_gpu_supports` (no device); the device solve measured in [ROCM_PERFORMANCE_ANALYSIS_2026-10-13.md](ROCM_PERFORMANCE_ANALYSIS_2026-10-13.md).

| | CPU | CUDA | ROCm / Vulkan |
|---|---|---|---|
| `ppcg` stop | relative: `eeq_ppcg_tol` 1e-12, `|P r| <= tol (|b|+1)` | absolute on `|r|`: `pcg_tolerance` 1e-10 | ROCm: as CUDA (absolute `|r| <= tol`); Vulkan: no device `ppcg` |
| far-field shortcut (`farfield_shortcut`) | gamma_ij*r >= 6.5 (`COULOMB_FARFIELD_X2` = 42.25) | gamma_ij*r >= 6.0 (`GFNFF_COUL_FARFIELD_X2` = 36) | none |
| `localinv` preconditioner | yes | yes | ROCm: device build (default) or host build + device apply (Oct 13, 2026; median CG iterations 56-69 -> 12-16 at 1410-7320 atoms); Vulkan runs the CPU |

Aliases `schur_cholesky` = `cholesky`, `projected_pcg` = `ppcg`. **An unknown name is an error** that lists the valid ones (it used to mean `auto`).
The thresholds are arguments: `-gfnff.eeq_ppcg_min_atoms`, `-gfnff.eeq_ppcg_min_nfrag` (0 = auto never picks `ppcg`).

`ldlt` equals `cholesky` for an SPD A_nn (every gas-phase system). It is not a fallback: for an indefinite A_nn (implicit-solvation reaction
field) the robust path is the augmented-LU fallback of `cholesky`, which `ldlt` does not reproduce.

## Parameters (module `gfnff`, category EEQ)

The EEQ solver PARAMs are defined in the `gfnff` block of `src/core/energy_calculators/ff_methods/gfnff.h` (category EEQ, GPU-only ones under
GPU); `curcuma -help_module gfnff` lists them. Until Oct 2026 they formed a registry module `eeq_solver`, which was merged into `gfnff` because
only GFN-FF uses the solver. Names, types, defaults and aliases did not change. Every spelling acts:

| spelling | lands in | precedence |
|---|---|---|
| `-eeq_solver.<name>` | `controller["eeq_solver"]`, the JSON scope `EEQSolver` reads | wins over the other spellings |
| `-gfnff.<name>`, flat `-<name>` | `controller["gfnff"]`, reaches the solver as a flat key (CPU), through `eeq_value()` (GPU) and, for `solve_method`, through `fixEEQSolverChoice()` | used when the `eeq_solver` scope does not set the name |
| `-gfnff.solve`, `eeq_max_iterations`, `eeq_tolerance`, `eeq_accuracy` | older gfnff names, copied into the `eeq_solver` scope by `GFNFF::forwardEEQSolverParams()` | CPU: the last three beat `max_pcg_iterations` / `pcg_tolerance` / `pcg_large_system_tol_factor` of the gfnff scope; the GPU reads only those canonical names. `solve` acts only while every `solve_method` is auto |

`ParameterRegistry::canonicalModule()` answers lookups by the name `eeq_solver` (defaults, aliases such as `ppcg_precond`, `-help_module
eeq_solver`) from `gfnff`. The GPU wrapper reads `skip_phase2` and `allow_unconverged_charges` from the gfnff scope only, so the `-eeq_solver.`
spelling of these two does not reach the GPU path. `accuracy` has no reader and no effect.

Measurements behind some defaults (one machine each, not repeated):

- `eeq_refactor_eps_bohr` (cached Cholesky factor, default 0 = off): 450-atom water cluster 16.05 ms per energy call cached vs 15.53 ms uncached,
  so below `eeq_ppcg_min_atoms` it saves nothing; in MD at eps 0.05 Bohr one refinement step (`eeq_refine_iters` 1) left a charge error of 1e-5 to
  8e-4 e, three steps about 1e-7 e.
- `eeq_ppcg_tol`: 1e-10 instead of 1e-12 saves about 15 % of the iterations.
- `gpu_block_jacobi_max_frag_atoms` / `gpu_block_jacobi_max_nfrag` (CUDA `pcg`, EEQ solver benchmark, Sep 2026): at nfrag 500, largest fragment
  3 / 250 / 500 / 1410 atoms took 5.7 / 7.6 / 11.4 / 53.0 s per 2 steps; 3-atom fragments at nfrag 500 / 1000 / ~2440 took 5.8 / 21.1 / > 150 s;
  a polymer_2x-like system (nfrag ~1500) exceeded 120 s. The preconditioner build is a serial host loop with two device synchronisations per fragment.

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

## Local-inverse preconditioner of ppcg (`eeq_ppcg_precond localinv`, default since Oct 10, 2026)

Claude Generated (Oct 10, 2026) — 🤖 AI-generated, ⚙️ machine-tested (`cli_gfnff_07_eeq_localinv`, measurements below); human production testing pending.

`-gfnff.eeq_ppcg_precond localinv` (or `-eeq_solver.eeq_ppcg_precond`; the default since Oct 10, 2026; `diag` restores the diagonal Jacobi) replaces the diagonal Jacobi of the projected PCG, CPU and `-gpu cuda` (WP7-E).
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
  thread count are bit-identical between runs. The GPU path as a whole is not bit-reproducible (two diag runs: charges within 9e-12 e).
- **CUDA only** on the device side: `pcg`, `batched`, `lu`, `ldlt`; ROCm has `ppcg` with Jacobi or `localinv` (device build, host build as comparison). ROCm device build: 3.6 ms (1410 atoms), 3.9 (2910), 11.0 (6200, urea400-water1000), 7.9 (7320) per build against 11.8 / 15.1 / 247 / 37 ms for the host build (median of 8-15 builds in 100 MD steps, subagent measurement); MD step urea400-water1000 124.4 ms with the host build, 88.7 ms with the device build (Jacobi 104.7), peo201x2-water1500 58.7 -> 55.2 ms (CPU 144.9), author's measurement, 3 alternating runs. Test hooks `CURCUMA_LOCALINV_TEST_INDEFINITE`/`_REJECT` exist on ROCm as on CUDA; no ctest runs the ROCm device build.

### GPU: build on the device (Oct 10, 2026)

`buildLocalInverseOnDevice` (`cuda/eeq_solver_gpu.cu`), requested by `prepareLocalInverseGPU` (`gfnff_gpu_method_impl.h`) when the rebuild rule fires,
runs inside the WP7-E solve after the EEQ matrix is assembled: `k_li_count` (block size and kept entries per atom, distance scan with a float
prefilter and the exact double test), `k_li_blocks` (one 128-thread block per atom, launched in buckets of block size <= 32/64/128/larger so the
shared memory fits the bucket: neighbours in ascending order, B = A[E,E] gathered from the device matrix `d_A`, i.e. the values of
`k_eeq_build_matrix`, Cholesky + one column of B^-1 in shared memory, global scratch if a block does not fit), `k_li_symmetrize` (G = (Gl + Gl^T)/2).
No floating-point atomics. The host build stays: `CURCUMA_LOCALINV_HOST_BUILD=1` (comparison) and the multi-GPU split (no full matrix on one device).

- **No S on the device.** The device applies z = P G r (P = orthogonal fragment projection), not G - H S^-1 H^T: dense S = C G C^T costs
  4.2 ms (potrf) + 9.4 ms (potri) at nfrag 1502 on the RTX 5080 (FP64 1/64), a sparse factor needs the host. P G P is symmetric and positive
  definite on null(C) whenever G is; it removes the nfrag <= 4096 limit. Cost in iterations (MD, median per step): polymer_2x 15 vs 15 (host build),
  peo201-water500 13 vs 13, peo201-chain 12 vs 12, urea400-water1000 11 vs 8.
- **Definiteness** without the host's sparse Cholesky of G: (1) every local pivot must be > 0, else the build is rejected (Jacobi until the next
  rebuild, a warning, no GPU fallback); (2) rz <= 0, pAp <= 0 or a stall stop the CG; (3) after convergence one more matvec checks the TRUE residual,
  |P(b - A q)| <= 2 tol (observed ratio to tol: at most 0.999 in all runs below). (2) and (3) redo the solve with Jacobi from the same start
  and are counted GPU fallbacks (`-gpu_strict true`: exit code 3). Test hooks: `CURCUMA_LOCALINV_TEST_INDEFINITE=1` (G -> -G), `=2` (every 10th
  diagonal negated), `=3` (true-residual check forced to fail); `CURCUMA_LOCALINV_TEST_REJECT=1` rejects every build (test 07).
- **Timing hooks**: `CURCUMA_PPCG_ITERS=1` prints `LOCALINV GPU device build ok= build_ms= ...` and the true-residual ratio;
  `CURCUMA_PPCG_GPU_TIMING=1` synchronises after each phase of the WP7-E call and prints matrix / build / CG / residual times.

MD 50 fs, dt 1 fs, 300 K, `-gpu cuda -gpu_strict true -threads 4`, per-step means over steps 1-49, RTX 5080; polymer_2x n = 3 runs (alternating
configurations), the others n = 2; [range of the run means]. "host" = `CURCUMA_LOCALINV_HOST_BUILD=1` (the build before this change).
Foreign load during the series: another session's CPU job at ~500 % for the first run, later none above 100 %; load average 5.6-15.5.

| polymer_2x | iterations | build ms | EEQ ms/step | step ms |
|---|---:|---:|---:|---:|
| diag | 68 | - | 32.4 [32.2-32.7] | 76.4 [74.7-79.5] |
| localinv rebuild 1.0, host build (7 builds) | 15 | 54.5 [51.4-60.5] | 23.8 [23.3-24.7] | 67.2 [65.8-69.7] |
| localinv rebuild 1.0, device build (7 builds) | 16 | 2.3 [2.2-2.4] | 16.5 [16.4-16.6] | 60.9 [59.0-64.4] |
| localinv rebuild 0, host build | 15 | 59.3 [50.5-76.5] | 76.8 [67.8-94.5] | 122.8 [110.8-146.7] |
| localinv rebuild 0, device build | 15 | 2.1 [2.1-2.2] | 18.1 [17.9-18.3] | 63.4 [60.3-69.0] |

Phase split of one device-build solve (polymer_2x, `CURCUMA_PPCG_GPU_TIMING`, n = 8 steps, median): matrix 6.6 ms, build 2.0, CG 5.2 (15
iterations, 0.34 ms each; diag 0.33), true residual 0.3; the rest of the EEQ time (~3.8 ms, also in diag) is outside the solver call.
Epot after 50 fs -902.161126 Eh in all 15 runs; per-step charges vs diag: device rebuild 0 3.5e-11 e, rebuild 1.0 3.6e-11, host 3.4e-11.

| other systems, EEQ ms/step (iterations) | diag | host build, rebuild 0 | device, rebuild 0 | device, rebuild 1.0 |
|---|---:|---:|---:|---:|
| peo201-water500 (2910 atoms) | 7.0 (58) | 23.9 (13) | 4.7 (13), build 1.0 ms | 3.9 (14) |
| urea400-water1000 | 22.9 (37) | 279 (8) | 20.0 (11), build 3.2 ms | 17.5 (12) |
| peo201-chain (nfrag 1) | 3.7 (56) | 15.6 (12) | 2.4 (12), build 0.9 ms | 1.6 (13) |

- **Rebuild threshold**: rebuild 1.0 stays faster than 0 on the device on all four systems (a build costs 1-3 ms, one more iteration 0.1-0.35 ms), so
  the default 1.0 is kept for the GPU as well.
- urea400-water1000: charges vs diag up to 1.1e-10 e per step for the device AND the host build (1.10e-10); two diag runs agree to 6e-14. This is the
  spread of two solves converged to the absolute GPU tolerance 1e-10 on |r|, not an error of the build.

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

## Threaded symmetric matvec of the CPU ppcg (Oct 10, 2026)

The dense products `A_nn * p` (once per iteration) and `A_nn * q` (residual start) in `solveWithProjectedPCG` (`eeq_solver.cpp`) read only the blocks (I,J), J >= I,
of fixed `kPpcgMatvecBlockRows` = 64 tiles (half the bytes; both triangles of `A_nn` are stored and identical, the assembly evaluates one erf per pair i > j).
Each tile gives `y_I += A_IJ x_J` and `y_J += A_IJ^T x_I`. Block row I is handled by one thread (atomic work counter, `m_ppcg_pool`/`m_ppcg_threads`, std::thread fallback);
the transposed parts go to a private slice per block row (nblocks x N doubles), `y_J = own_J + sum_{I<J} slice_I[J]` is summed with I ascending. No atomics, no
thread-dependent order: the hash of every product (213 per SP, both phases) is identical at 1, 4, 32 threads (checked once with a temporary hook, removed); it differs
from the former row-wise product by rounding only (charges 1e-14 e, same iteration counts, SP -898.61727307 Eh). `CURCUMA_PPCG_MATVEC_BLOCK=<rows>` overrides the block size.
polymer_2x `-sp`, Phase-2 solve (median of 3): 1 thread 1138 -> 937 ms, 4: 868 -> 467, 16: 855 -> 350, 32: 827 -> 357. Machine-tested only; MD numbers in `AIChangelog.md`.
The matrix itself differs by ~1e-13 between thread counts (threaded assembly), so charges are not bitwise equal across `-threads`; energy and iteration counts are (test 07).

## Related

GPU MD step breakdown and the CN-list / overlap / staged-upload switches: [GFNFF-PERFORMANCE-ROADMAP.md](GFNFF-PERFORMANCE-ROADMAP.md), [GPU_TUNING.md](GPU_TUNING.md).

## Not tested / open

- ROCm: the rejection of `pcg`/`batched`/`lu`/`ldlt` is tested by name only (`eeq_gpu_supports`); the device `ppcg` was run on one RX 9060 XT (gfx1200), see the analysis document. Vulkan runs GFN-FF on the CPU, so the CPU solvers apply.
- Multi-GPU split and the batch-worker device pool were not run after this change.
- localinv: the GPU fallbacks were triggered only through the test hooks (no natural rejection or breakdown in any run); not run with ALPB
  solvation, `eeq_distance_cutoff > 0`, transition metals, MD longer than 50 fs, optimisations, nfrag > 4096, or `eeq_localinv_ext_bohr` large
  enough for the global-scratch path of `k_li_blocks` (blocks above about 150 atoms at 99 KB shared memory per block, estimated; polymer_2x has at most 86).

## Test note (Oct 10, 2026)

All tests that use the GPU carry `RESOURCE_LOCK "gpu"`, so `ctest -j8` runs them one after another. Without it `gfnff_gpu_vs_cpu_*` failed with
"Insufficient GPU memory" when several GPU tests overlapped on one 16 GB card: 4 of 7 `-j8 -L gfnff` runs on the state before the Oct 9 solver changes
(worktree at `a2b7c5f0`, built with the same flags) and 3 of 7 after them; with the lock 0 of 7 (87 tests each). Serial runs were green throughout.

## Validation of `localinv` (Oct 10, 2026)

`scripts/localinv_validation.py` (sub-commands run, grad, opt, md, analyze; memory capped through systemd-run) compares ppcg+`localinv` with ppcg+`diag` and the exact dense
solver. Machine-measured: 703 single points (588 GMTKN55 structures incl. all 327 charged ones, 95 MOR41 complexes incl. 72 with a metal or Se, 10 ALPB systems), |dE| <= 1e-10 Eh
against `cholesky`; 13 gradients, |dg| <= 9e-13 Eh/Angstrom; NVE 2 ps on polymer_2x (GPU) and peo201-water500 (CPU), drift and spread of Etot equal to `diag`; no Jacobi
fallback with the default radii. With truncated radii (4/6 Bohr) 6 of 588 builds were rejected (small ions) and caught by the Jacobi fallback. Not measured: open shell, `eeq_distance_cutoff` > 0,
MD beyond 2 ps, ROCm. A -opt of peo201-chain ends in different basins from run to run (thread count, rounding at 1e-13); at both end geometries cholesky, diag and localinv give the same energy.

## Phase 2 without per-step copies (Oct 11, 2026)

`EEQSolver::calculateFinalCharges` no longer copies or zeroes per step: (1) the plausibility test after the solve reads A_nn through a `MatrixCRef` view (was a 428 MB copy at N = 7320); (2) the
packed distance triangle `m_external_srab` (GFNFF::m_shared_srab, valid for the whole energy call, read only) is used in place and its atoms-too-close scan runs on the pool (was `setZero` + 214 MB copy + serial
scan); (3) the constraint rows are filled from `fraglist` in O(N), and the same atom -> fragment map is handed to the projected PCG (`m_ppcg_atom_frag`), which otherwise reads it from C. Machine-tested: md5 of the
whole `-dump_gradient` file identical old/new in 58 comparisons (6 bulk/caffeine incl. polymer_2x at 1 and 16 threads, 31 GMTKN55 incl. 8 charged and 3 open-shell, 15 MOR41, 6 solver-path runs caffeine/peo201-chain with
`diag`, `cholesky`, `lu`); 50 fs MD polymer_2x (16 threads) `md_diagnostics` file identical in all 50 steps; with `-eeq_solver.dump_charges true` even `e_coulomb_phase1/2` are bit-identical (5 systems). polymer_2x MD 20 fs,
alternating, n=3: `eeq_solve` 239 -> 164 ms, step 469 -> 392 ms at 16 threads; 1 thread (noisy, other load) 390-527 -> 312-402 ms. Not done: one pass over A_nn for both energy forms (gemm instead of two gemv, rounding not guaranteed the same, at most 13 ms).

## Distance buffer, plausibility test and CN without serial passes (Oct 11, 2026)

Three exact changes in the CPU step (machine-tested, not human-tested):
- `GFNFF::computeSharedDistances`: the squared-distance twin `m_shared_sqrab` had no reader and is gone (with its accessors); no `setZero` (2 x 214 MB at N = 7320), because only the
  off-diagonal slots j < i are read (`dist(i, j)` in `EEQSolver::calculateFinalCharges`); the atoms-too-close test (r < 1e-10) is computed while writing and handed over with
  `setExternalDistances(&srab, &flag)`, so the EEQ solver's second scan over the buffer is skipped (same message, same fallback charges). The distance is written as
  `std::fma(dz, dz, std::fma(dy, dy, dx * dx))`: this is how GCC (-O3 -march=native, default `-ffp-contract=fast`) contracted the old expression in its vectorised loop; the
  rewritten loop has a scalar remainder that GCC contracts differently (dy first), which changed ~280 of 26.8 M distances in the last bit and the MD charges at 1e-15 (found by the MD
  md5, invisible in single points, which reuse the topology-stage charges).
- Plausibility test after the Phase-2 solve: both q^T A_nn q in one pass over the rows of A_nn on the pool, fixed thread-independent order (4 partial sums per row, blocks of 64
  rows summed in order) instead of two MKL dgemv. `e_coulomb_phase1/2` differ from before by 1-3 ulp (MKL's order cannot be rebuilt; MKL gemv per row block was tried and also differs);
  they only feed the ratio test (> 5) and the `-eeq_solver.dump_charges` file, never charges, energy or gradient.
- `CNCalculator::calculateGFNFFCNWithNeighbors`: pair scan i < j and CN sums parallel over rows when given the pool (`threads > 1`, N > 64); the neighbour lists are merged in the old
  order (partners k < i ascending, then j > i), so each CN sum runs in the old order. One thread keeps the old serial code.

Checks (md5 of the whole `-dump_gradient` file, old = HEAD e3eac812 build): 60 of 60 identical (8 bulk incl. polymer_2x at 1 and 16 threads and GFN2 caffeine, 6 solver paths
diag/cholesky/lu, 31 GMTKN55 `[20::80]` incl. 4 charged and 6 open shell, 15 MOR41); `md_diagnostics` files identical: polymer_2x 50 fs 16 threads, 4 fs 1 thread, urea400-water1000 and
peo201-water500 20 fs 16 threads. GPU (`-gpu cuda`, the changed code runs there only at setup): SP and 100 fs MD old/new inside the GPU run-to-run spread (n=3), step 21.26 -> 21.17 ms.
Timing polymer_2x, MD 20 fs, alternating, n=3 (median of medians): 16 threads step 294.2 -> 227.8 ms (`eeq_solve` 157.9 -> 136.5, `cn` 19.4 -> 5.5, the distance buffer is not
timed separately), wall clock (250-50 fs) 287.8 -> 224.1 ms per step; 1 thread step 674.9 -> 670.8 ms, where `cn` reads 19.0 -> 23.8 ms although the one-thread CN code is unchanged
(not explained; process page faults and user time equal). Not tested: ROCm, periodic systems, ALPB, `eeq_matrix_rebuild_eps_bohr` > 0, other compilers (the FMA spelling fixes GCC's choice).
The `GFNFF` class layout changed: rebuild the CUDA plugin and the test executables with it (`make` without target); stale `test_gfnff_gpu` aborts in malloc.

## Phase-2 matrix built in tiles (Oct 11, 2026)

`EEQSolver::calculateFinalCharges` builds the N x N Coulomb block in 64 x 64 tiles of the lower triangle instead of one row per thread with a stride partition (machine-tested, not
human-tested). Each element is computed once with the former expression and inputs (cutoff 0, far field 1/r, erf(r/sqrt(alpha_i+alpha_j))/r, diagonal gam_i + sqrt(2/pi)/sqrt(alpha_i))
into a thread-local tile buffer, then stored row by row to A(i, J-block) and transposed to A(j, I-block). Threads take contiguous ranges of the tile list. For N > 2048 the stores are
non-temporal (MOVNTI, `stream_store_double`): the 428 MB block at N = 7320 fits no cache, and a normal store reads each line first. The expression contains no a*b+c, so FMA
contraction cannot differ; no element depends on the order or the thread. Unchanged: the matrix-cache reuse path (`eeq_matrix_rebuild_eps_bohr` > 0), constraint rows, Born matrix
(ALPB), Phase 1, GPU code. No full `setZero` was or is needed: every cell of the N x N block is written; only the constraint rows/columns are zeroed.

Checks: hash of the whole (N+nfrag)^2 matrix, old/new, each with the buffer pre-filled with 0 and with 777 (shows that every cell is written), 21 cases / 60 matrices identical
(polymer_2x, peo201-chain, caffeine at 1 and 16 (4) threads, charged +1 peo201-chain and hydronium-water2, ALPB water caffeine/peo201-chain incl. Born matrix, MOR41 PR16 (Pd) and ED08
(Ir), peo201-water500, urea400-water1000, `eeq_distance_cutoff 20` at 1 and 16 threads, MD 4 steps peo201-chain (also with `eeq_matrix_rebuild_eps_bohr 0.05`) and polymer_2x);
md5 of the whole `-dump_gradient` file identical in 60 of 60 (8 bulk incl. GFN2 caffeine, 6 solver paths diag/cholesky/lu, 31 GMTKN55 `[10::80]` with 6 charged and 5 open shell,
15 MOR41 `[2::6]`); `md_diagnostics` polymer_2x 50 fs 16 threads and 4 fs 1 thread identical; `e_coulomb_phase1/2` (`dump_charges`) identical (3 systems x 2 thread counts).
GPU (`-gpu cuda`): the CPU build runs there only during setup (2 builds per run, SP and MD alike, counted with a temporary hash print).
Timing polymer_2x (16 threads unless noted): matrix build alone (SP, temporary timer, n=3 x 2 calls, median) 37.5 -> 19.6 ms, 1 thread 55.0 -> 40.8 ms; MD 20 fs, alternating,
n=3: `eeq_solve` 135.8 -> 120.0 ms, step 231.1 -> 212.2 ms; 1 thread step 670.5 -> 655.4 ms; wall clock (250-50 fs) 224.6 -> 207.2 ms per step.
Not tested: ROCm, periodic systems, CPUs without SSE2 (plain stores), other compilers.
