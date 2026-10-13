# ROCm GFN-FF and GFN2: where the time goes and what to change (Oct 13, 2026)

Status: AI analysis, machine-measured, nothing implemented. Follows [ROCM_EVALUATION_2026-10-13.md](ROCM_EVALUATION_2026-10-13.md).
Machine: RX 9060 XT (gfx1200, 32 CU, 16 GB), Core Ultra 5 225, 10 threads. All gains marked "est." are extrapolations, not measurements.
The first version of this analysis used micro benchmarks and the built-in timers (`-verbosity 2|3`) because no profiler was installed;
`rocprofv3` was installed during the session and the GFN-FF numbers below it were re-measured with it. Micro benchmarks: [scripts/rocm_microbench/](../scripts/rocm_microbench/README.md).

## Hardware ceiling (measured, rocBLAS)

| n | dgemm | sgemm | ratio |
|---|---|---|---|
| 512 | 397 GF/s | 3.6 TF/s | 9 |
| 1024 | 402 GF/s | 6.5 TF/s | 16 |
| 2048 | 427 GF/s | 7.4 TF/s | 17 |

FP64 is about 1/17 of FP32 on this card. Every FP64-bound kernel pays this; FP32 or fewer FP64 operations are the levers. The 890M
numbers in `SQM_ROCM.md` come from a different FP64 balance and are not transferable.

## 1. GFN-FF (7320 atoms, 1502 fragments, MD step ~225 ms on ROCm, ~150 ms on the CPU)

Per-step phases from `-verbosity 2` (step with HB refresh, n=1 run; steps vary 200-243 ms):

| Phase | ms | Note |
|---|---|---|
| EEQ solve on the CPU (ppcg, `nfrag >= 16` routing) | 100-120 | sequential, ~50 % of the step |
| Coulomb + postprocess + DMA (GPU) | 88 | ~40 % |
| dlogdcn + HB/XB detection | 13 + 22 | HB refresh steps only |
| everything else on the GPU | see correction below | the first version of this table said ~5 ms; that was wrong |

**Correction (rocprofv3, same day).** The per-term "GPU ms" column of `-verbosity 2` shares one value between several terms (for
example Dispersion 5.4 ms and Repulsion 5.4 ms are one phase), so it cannot be read as per-kernel time. Kernel trace of 41 MD steps on
7320 atoms, baseline binary (`rocprofv3 --kernel-trace`, n = 1 run, avg per launch): `k_coulomb_dense` 49.9 ms, **`k_dispersion_gather`
37.9 ms**, `k_dc6dcn_per_pair` 8.1 ms, `k_generate_cn_pairs_count`+`_write` 6.4 + 6.4 ms (37 launches), `k_repulsion` 2.5 ms (82 launches),
`k_hbonds` 3.1 ms; total kernel time 117 ms per step. Dispersion was the second-largest kernel and the plan above missed it.

### F1. Coulomb kernel: no far-field shortcut, every pair twice (measured, large)
`k_coulomb_dense` (`ff_methods/rocm/gfnff_rocm.hip:581`) evaluates `sqrt`, `1/sqrt`, `erf`, `exp` and three FP64 divisions for all N^2
ordered pairs. CUDA skips `erf`/`exp` for `(gamma r)^2 >= 36` (`GFNFF_COUL_FARFIELD_X2`, `gfnff_kernels.cuh:152`) and evaluates each
pair once (`k_coulomb_tiles`); CUDA measured `k_coulomb_tiles` 11.2 -> 2.4 ms on polymer_2x, 99.9 % of pairs far field
([GPU_TUNING.md](GPU_TUNING.md)).
- Standalone copy of the kernel on this GPU, 7320 atoms from `peo201x2-water1500.gfnff-min.xyz`, random alpha in [0.3, 1.0] and charges in
  [-0.4, 0.4], no PBC, best of 5: **49 ms -> 8.5 ms with the far-field branch** (block 64/128/256 identical). Energy and
  |gradient| equal to all 12 printed digits (E = 6.223329600489, |g| = 1.498901515311e+00).
- HIP `erf` is exactly 1.0 for x >= 6: 0 of 3e8 sampled points in [6, 9] differ; between 5.9007 and 6.0 21 % are still below 1.0, so the
  threshold must stay at 6 (36 in squared form), not CUDA's saturation point 5.9007.
- Not checked: realistic alpha from a real EEQ, the periodic path (`applyMIC`), the difference between the 49 ms kernel and the 88 ms
  phase (the other ~40 ms are postprocess/DMA and not profiled).
- Est. step at 7320 atoms: -40 ms (225 -> ~185 ms) from the kernel alone. Halving the pair count (tiles, each pair once, as CUDA) is a second
  step; it needs shared-memory gradient accumulation as in `k_coulomb_tiles`.

### F2. EEQ on the CPU (largest single block)
For `nfrag >= 16` (`-gfnff.eeq_rocm_cpu_fragment_threshold`) ROCm routes to the CPU PPCG (log: "EEQ GPU Phase 2: nfrag=1502 >= 16 - routing
to exact CPU PCG"), 100-120 ms per call, run in parallel with charge-independent GPU work but still 84 % of the calculation in sequence
("Sequential work (EEQ + finish): 190.0 ms"). CUDA has the device PPCG (WP7-D, `k_ppcg_*`, local-inverse preconditioner `k_li_*`):
EEQ 16.1 -> 12.3 ms on polymer_2x (RTX 5080).
- Est. with a ported device PPCG and far-field matrix build: dense FP64 matrix 7320^2 x 8 B = 429 MB per matvec; at the card's memory
  bandwidth (not measured here) a few ms per matvec; 16 iterations gives ~20-40 ms, against 100-120 ms now. Combined with F1, step
  ~225 -> ~100-120 ms, i.e. about 1.3-1.5x faster than the CPU instead of 1.4x slower. This is an estimate with no prototype.
- Cost: a port of ~1000 lines of CUDA (`eeq_solver_gpu.cu` WP7). The ROCm solver names are limited to `auto|cholesky` today.

### F3. Single-fragment EEQ (peo201-chain, N=1410)
Normal steps 13.8 ms, but a topology/HB refresh step took 101 ms with 94 ms CPU EEQ (`-verbosity 2`, 2 observed steps, no statistics).
The dense Cholesky for N<=~2000 is acceptable; the refresh cost is the CPU side (matrix build + Schur). Not analysed further.

### F4. FP32 for Coulomb (not recommended before F1)
FP64 erf/exp dominate only until F1; after F1 the near field (a few 0.1 % of pairs) is small. An FP32 near field would break
"100 % with respect to the reference" and would need an error budget; defer.

### F1/F5 result: implemented on branch `fix/rocm-coulomb-farfield` (uncommitted, 2026-10-13)
Changes (all in `ff_methods/rocm/`, CUDA and CPU files untouched):
- Far-field branch (gamma r >= 6, i.e. r^2 >= 36 (alp_i + alp_j)) in `k_coulomb`, `k_coulomb_gather`, `k_coulomb_dense` and the two EEQ matrix
  builds (`gfnff_rocm.hip`, `eeq_solver_hip.hiph`); constant `GFNFF_COUL_FARFIELD_GAMMA_R` in `gfnff_kernels_hip.h`. sqrt and division are
  kept (no `rsqrt`), so the energy path stays bit-identical; the dropped Gaussian term is below 2.3e-16 relative.
- `k_dispersion_gather`: one workgroup (64 lanes) per atom with a shared-memory tree reduction instead of one thread per atom (the CSR rows of
  neighbouring threads were ~2600 entries apart, so no load was coalesced). Per-atom sums are added in a different order.

Measured, 7320 atoms (`peo201x2-water1500.gfnff-min.xyz`), `rocprofv3` kernel trace, 41 steps, n = 1 per binary:

| Kernel | before | after |
|---|---|---|
| `k_coulomb_dense` | 49.9 ms | 20.4 - 23.3 ms (two runs) |
| `k_dispersion_gather` | 37.9 ms | 11.2 ms |
| total kernel time per step | 117 ms | 63 ms |

The standalone prototype predicted 8.5 ms for the Coulomb kernel; the real kernel is 2.4x slower than the prototype (real alpha values make
fewer pairs far field, plus `applyMIC` and NaN checks; not decomposed).
Whole-process MD time, 3 alternating runs per binary, increment of 100 `-MaxTime` units (about 100 steps), spread < 0.02 s:

| System | atoms | before | after | CPU (`-gpu none`) |
|---|---|---|---|---|
| peo201x2-water1500.gfnff-min | 7320 | 20.5 s | 15.0 s | 14.5 s |
| urea400-water1000 | 6200 | 25.0 s | 18.0 s | 18.0 s |

So ROCm went from 1.4x slower than the CPU to equal (ratio 1.03 and 1.00). The GPU kernels (63 ms) are now well below the step time
(~150 ms): the CPU EEQ (F2) is the critical path.

Checks against the unchanged baseline binary (copy of the pre-change `release/curcuma` and plugin):
- `-sp` energies and the MD Epot/Ekin/Etot columns after 200 steps equal to 6 printed decimals on 1410, 2910, 6200, 7320 atoms; last-frame
  coordinates equal to the 6 printed decimals (max |dx| = 0 on 1410, 6200, 7320 atoms).
- `-verbosity 2` at step 1: total energy equal to 10 decimals and gradient norm equal to 7 digits between baseline, new ROCm and CPU on
  1410, 6200, 7320 atoms (for 7320: 9.290235e-02 in all three).
- `ctest -R "gpu|rocm"`: 34 of 34 pass (includes `cli_gfnff_gpu_02_rocm_singlepoint`, 32 xtb GPU validation tests that do not use these kernels).
- Not checked: gradient components beyond the norm and a printed-coordinate trajectory (a bit-exact comparison would need more digits);
  periodic cells (`applyMIC` path); charged or multi-fragment systems with very different alpha; other GPUs; the CPU-vs-ROCm trajectory
  divergence over long runs.

### F2 result: device projected PCG on ROCm (branch `feature/rocm-eeq-ppcg`, Oct 13, 2026)
Port of CUDA WP7-E into `ff_methods/rocm/eeq_solver_hip.hiph` (written by a subagent, checked by the author of this document):
projected PCG with the diagonal Jacobi preconditioner, `rocblas_dsymv` for the matrix-vector product, fragment sums as one block per
fragment (the first version with `atomicAdd` was slower than the CPU route, 172 ms/step). Not ported: `localinv` (`k_li_*`), multi-GPU,
input events. `-eeq_solver.solve_method ppcg` is accepted with `-gpu rocm`; many-fragment systems whose solver resolves to ppcg use it
by default, a stalled device solve uses the CPU EEQ for that step, and giving `-gfnff.eeq_rocm_cpu_fragment_threshold` explicitly forces
the old CPU route. CUDA code paths are unchanged in the source; the CUDA build was not compiled or run here.

MD time per step (t(MaxTime 120) - t(MaxTime 20)) / 100, 3 alternating runs per binary, spread < 2 ms (author's own measurement):

| System | atoms / fragments | ROCm before | ROCm now | CPU |
|---|---|---|---|---|
| urea400-water1000 | 6200 / 1400 | 182 ms | 125 ms | 180 ms |
| peo201x2-water1500.gfnff-min | 7320 / 1502 | 150 ms | 123 ms | 146 ms |

ROCm is now 1.4x (6200 atoms) and 1.2x (7320 atoms) faster than the CPU. The subagent measured 157 -> 130 ms and 185 -> 128 ms for the
same systems. Charges against the exact Cholesky solution: max |dq| 3.7e-11 to 1.1e-10 over 20 MD steps on 1410, 2910, 6200 and 7320
atoms (subagent measurement, 4 systems; baseline CPU PPCG vs exact: up to 6.0e-11); 35-79 CG iterations. Epot/Ekin/Etot of a 100-step
MD equal to 6 printed decimals between new ROCm, previous ROCm and CPU on the same four systems (author's own check, 2 printed rows each).
`ctest -R "gpu|rocm|eeq|gfnff"`: 104 of 104 pass; `eeq_gpu_supports` was changed to expect `auto, cholesky, ppcg` for ROCm.
Remaining cost: 0.82 ms per CG iteration (dense `dsymv`), EEQ about 55-65 ms per step at 7320 atoms; `localinv` would cut the
iterations by about half on the CPU test (14 vs 32 at 2910 atoms) and is the next lever, together with the dispersion/Coulomb kernels
(about 20 + 11 ms).

### F6 result: CN pair skin, lower-triangle EEQ matrix, Coulomb pair once (branch `feature/rocm-gfnff-step-speed`, Oct 13, 2026)
Profile before (rocprofv3, 7320 atoms, per step, GPU 91 % busy): EEQ `hemv` 57 ms, Coulomb 19.6, dispersion 10.0, `dc6dcn` 7.9, EEQ matrix build 7.1, CN pair generation 12.3 (list rebuilt every step), repulsion 4.7, H bonds 2.8. Three changes by a subagent, checked by the author:
- CN pair-list Verlet skin on ROCm (CUDA K1; `-gfnff.gpu_cn_pair_skin_bohr`, trait `has_cn_pair_skin`): 12.1 -> 0.95 ms per step, list rebuilt 4 times in 51 steps.
- EEQ matrix written as the lower triangle only, coalesced (`k_eeq_build_matrix_hip`): 7.13 -> 4.95 ms. The LU fallback (dgetrf) still gets the full matrix; that path did not run in any test.
- `k_coulomb_tiles_hip`: each pair once, 128-atom tiles (replaces `k_coulomb_dense`): 19.7 -> 9.75 ms; Coulomb energy equal to 10 decimals.
- Not kept: extrapolated warm start of the PPCG: median CG iterations 37 -> 36 (6200 atoms) and 69 -> 67 (7320), 100 solves each; below the 10 % threshold, and the wrapper cannot tell MD from optimisation.

MD time per step ((t120 - t20)/100, 3 alternating runs, spread < 2 ms, author's own measurement): urea400-water1000 124.6 -> 106.1 ms (CPU 179.8), peo201x2-water1500 122.7 -> 101.2 ms (CPU 145.2). Subagent: 127.2 -> 107.5 and 130.7 -> 108.8 ms, profile step period 130.1 -> 107.8 ms. Energy, Epot/Ekin/Etot and last-frame coordinates equal to the printed 6-8 decimals against the previous binary on 7 systems (subagent) and on 4 bulk systems (author); `ctest -R "gpu|rocm|eeq|gfnff"` 104 of 104. Not checked: periodic cells, the LU fallback, other GPUs.
Remaining per step (7320 atoms): EEQ `hemv` 57 ms (69 CG iterations, memory-bandwidth bound), dispersion 10.0, `dc6dcn` 8.1, Coulomb 9.75, EEQ matrix build 5.0.

### F7 result: host-built local-inverse preconditioner on ROCm (branch `feature/rocm-eeq-localinv`, Oct 13, 2026)
Port of the CUDA host-build mode (`EEQSolver::buildLocalInverse` on the CPU, CSR upload, device apply z = G r - H S^-1 C G r, fallback to Jacobi on rz <= 0, pAp <= 0 or stall; the S^-1 product with one block per row). Not ported: the device build, the true-residual check, multi-GPU. Median CG iterations per step (100 solves): peo201-chain 56 -> 13, peo201-water500 58 -> 14, urea400-water1000 37 -> 10, peo201x2-water1500 69 -> 15. Max |dq| against the exact Cholesky solution over 20 steps: 1.7e-11, 2.6e-11, 1.5e-10, 1.4e-11 (subagent; Jacobi: 4.0e-11, 3.7e-11, 1.1e-10, 3.7e-11). Epot/Ekin/Etot of 100 steps equal to 6 decimals against the previous binary and the CPU on the four systems (author's check).
MD time per step (3 alternating runs, spread < 1 ms, author's measurement): peo201x2-water1500 (7320 atoms) 101.9 -> 58.7 ms (CPU 146.0); urea400-water1000 (6200 atoms) 104.7 -> 124.1 ms (CPU 179.5), slower than Jacobi. EEQ span per step at 7320 atoms 73.0 -> 22.6 ms (subagent, rocprofv3, 29 steps).
Cause of the urea result (subagent measurement): each host build takes a median 247 ms (7320 atoms: 37 ms), 15 builds in 100 steps, 204 ms of each in the serial sparse Cholesky check of G; the cause of the difference between the two systems was not examined. With `-eeq_solver.eeq_localinv_rebuild_bohr 2.0` (n = 1) 7 builds and 106.6 ms/step. Open: device build, a cheaper or asynchronous definiteness check, or a ROCm default of `diag`.

### F8 result: local inverse built on the device (branch `feature/rocm-eeq-localinv`, Oct 13, 2026)
Port of the CUDA device build (kernels `k_li_pack_f_hip`, `k_li_count_hip`, `k_li_blocks_hip`, `k_li_symmetrize_hip`, `k_li_test_perturb_hip`; z = P G r; true-residual check |P(b - A q)| <= 2 tol with Jacobi re-solve; blocks gathered from the lower triangle of the EEQ matrix). Default on ROCm; `CURCUMA_LOCALINV_HOST_BUILD=1` selects the host build. A build takes 3.6 / 3.9 / 11.0 / 7.9 ms at 1410 / 2910 / 6200 / 7320 atoms (host build 11.8 / 15.1 / 247 / 37 ms; medians of 8-15 builds per 100 steps; subagent). Latency bound, not FP64 bound; no block used the global-memory path at the default radii (block median/max 20-57 / 87-91 atoms; with `eeq_localinv_ext_bohr 14` at 7320 atoms 1913 atoms did, build 52 ms).
Max |dq| against the exact Cholesky solution over 20 steps (device / host / diag): 2.0e-11 / 1.7e-11 / 4.0e-11 (1410), 3.3e-11 / 2.6e-11 / 3.7e-11 (2910), 1.0e-10 / 1.5e-10 / 1.1e-10 (6200), 3.1e-11 / 1.4e-11 / 3.7e-11 (7320); median CG iterations 13, 14, 12, 16. Epot/Ekin/Etot of 100 steps equal at 6 decimals against the host-build binary and the CPU (author's check, 4 systems).
MD time per step (3 alternating runs, author's measurement): urea400-water1000 124.4 (host build) -> 88.7 ms (CPU 180.1; Jacobi 104.7); peo201x2-water1500 58.7 -> 55.2 ms (CPU 144.9). Subagent: 127.2 -> 92.4 and 66.3 -> 62.8 ms. `ctest -R "gpu|rocm|eeq|gfnff"` 104 of 104; no ctest exercises the ROCm device build; the CUDA build was not compiled after the trait changes.

Cumulative at 7320 atoms: 157 ms (ROCm after the first PPCG port, subagent, n=3) -> 55.2 ms, CPU 144.9 ms, i.e. 2.6x faster than the CPU; before this work ROCm was 1.4x slower than the CPU.

### F9 result: dispersion kernels (branch `feature/rocm-dispersion-speed`, Oct 13, 2026)
Changes in `gfnff_rocm.hip`, `gfnff_kernels_hip.h` (subagent, checked by the author): (1) per-atom vectors u/v (partner element type x 7 references) so that C6 and dC6/dCN of a pair are 7-term dot products instead of the 7x7 double loop; (2) in the dispersion pair math no `sqrt`, no division by r (the factor cancels), the cutoff test on r^2, one reciprocal for t6 and t8; (3) fused gather `k_dispersion_gather_uv` that forms C6 and dC6/dCN per edge in gradient mode and skips `k_dc6dcn_per_pair`. The fused path requires the per-step C6 refresh and a bitwise symmetric C6 reference table (checked at upload); `CURCUMA_ROCM_DISP_FUSED=0` selects the unfused path with (1)+(2). Energy-only calls still run the per-pair kernel. Dropped: CSR-ordered copies of the static pair constants (7.44 -> 6.67 ms for +607 MB at 7320 atoms).
Kernel time per step at 7320 atoms (rocprofv3, 29 steps, n = 1): `k_dispersion_gather` 10.15 + `k_dc6dcn_per_pair` 8.27 ms -> `k_dispersion_gather_uv` 7.44 + `k_disp_uv` 0.08 ms; step period 59.1 -> 51.5 ms.
MD time per step ((t120 - t20)/100, 3 alternating runs, author's measurement): urea400-water1000 89.3 -> 78.3 ms (CPU 180.6), peo201x2-water1500 55.4 -> 47.4 ms (CPU 145.1). Results against the previous binary (subagent, 10 systems): `-sp` and energy-only paths identical, gradient (fused) max |dg| 4.4e-16 Eh/A (relative to max|g| <= 3.4e-15), 100-step MD averages <= 4e-13 Eh, last frames equal at 6 decimals; author: Epot/Ekin/Etot of 100 steps equal at 6 decimals on 4 systems; `ctest -R "gpu|rocm|eeq|gfnff"` 104 of 104 (the energy-only ctest tests run with `-gpu none`; the subagent ran `test_gfnff_energy_only_geometry` by hand with `rocm`, diff 0.0 Eh). Not run at that point: the frozen-CN energy-only refill after fused steps and the asymmetric-table fallback; periodic cells; other GPUs. Follow-up below. The same restructuring applies to the CUDA `k_dc6dcn_per_pair` (not edited).
Follow-up on the untested C6 paths (Oct 13, 2026, author): `test_gfnff_energy_only_geometry` now also runs on ROCm (ctest `gfnff_rocm_energy_only_geometry_{triose,caffeine}` fused, `..._unfused_...` with `CURCUMA_ROCM_DISP_FUSED=0`; 6 ROCm and CPU variants pass) and has a third scenario with `-gfnff.static_cn` (energy-only call after a gradient call, frozen CN). A mutation test (refill before `k_dispersion` switched off) did NOT make any of them fail: the refill is needed only when an energy-only call skips the weight refresh (`frozenCN()`) while the per-pair C6 is older than the frozen CN, and the workspace is created at the first call, which also freezes the CN, so both belong to the same geometry; no public switch combination reached the stale state. The refill is therefore a guard that no test exercises. The asymmetric-table fallback selects the unfused path, which `CURCUMA_ROCM_DISP_FUSED=0` covers.

A discrepancy that predates all these changes: `cl2-anion.xyz` with `-charge -1`, 100 MD steps: ROCm and CPU agree at step 0 and at 20 fs but not at 100 steps (Epot -0.9983 against -1.1458 Eh, CPU independent of the thread count, the binary from the start of the session behaves like the current one); cause not examined (TODO.md).

### F10 result: HB and repulsion kernels (branch `feature/rocm-hb-repulsion-kernels`, Oct 13, 2026)
Findings by the subagent (disassembly plus before/after timing; `rocprofv3 --pmc` SQ_INSTS_* return 0 on this card): on gfx1200 every FP64 `atomicAdd` compiles to a compare-and-swap retry loop (`global_atomic_cmpswap_b64`; 14 in `k_repulsion`, 26 in `k_hbonds`), and the pair lists are grouped by atom, so the lanes of a wavefront hit the same address and serialise. `k_hbonds` was additionally limited by branch divergence: the case type (1, 2, 3) changes about every 4.4 triples in list order (132,272 changes over 582,105 triples on urea), so most wavefronts ran every branch.
Kept: (1) non-bonded repulsion as one workgroup per atom over a per-atom neighbour list, rebuilt with the Verlet list (skin unchanged; bonded repulsion still uses `k_repulsion`); (2) HB gradient pieces written to per-triple scratch slots and summed per atom by `k_slot_gather`; (3) the three charge-scaling exponentials per triple computed once per HB list upload (frozen charges) and exp/pow/sqrt values of the energy part reused in the gradient; (4) triples grouped by case type with a counting sort on the host. Dropped: reciprocals instead of divisions in the HB gradient (3.17 -> 3.04 ms, rounding change), a far-field cut for repulsion (no exactly negligible region below 20 Bohr found), CSR copies of the per-pair constants (the kernel was atomic-bound). Switches `CURCUMA_ROCM_REP_GATHER=0` and `CURCUMA_ROCM_HB_GATHER=0` restore the old paths (not run as separate tests).
Kernel time per step (rocprofv3, 29 steps, n = 1; subagent): urea400-water1000 `k_hbonds` 15.5 -> 3.18 + 0.41 ms gather, repulsion 8.28 -> 1.66 ms; 7320 atoms HB 1.07 -> 0.28 + 0.09, repulsion 3.50 -> 1.33 ms (an earlier table quoted 2.9 / 4.9 ms for the 7320 system; the subagent's own baseline runs are used here).
MD time per step ((t120 - t20)/100, 3 alternating runs, author's measurement): urea400-water1000 78.1 -> 71.5 ms (CPU 180.3), peo201x2-water1500 47.7 -> 46.7 ms (CPU 145.5). The step gains less than the kernels because each HB list refresh now costs about 22 ms more on the host for 582k triples (sort 3 ms, permuting the arrays about 10 ms, slot tables about 9 ms; 6 refreshes in 41 steps): the GPU is busy only about 41 of 70 ms per step on urea, the CPU HB re-detection and other host work are the limit. Memory +115 MB on urea, +24 MB at 7320 atoms.
Correctness against the previous binary (subagent, 13 structures, caffeine `-opt` 38 identical steps): step-0 term energies equal at 10 printed decimals, |g| at 7 digits, MD rows and last frames equal at 6 decimals, full-precision averages within 2.3e-13 Eh; against `-gpu none` terms within 1e-10 and |g| within 1.1e-9. Author: Epot/Ekin/Etot of 100 steps equal at 6 decimals on 6 systems (incl. fi-pyridine, acetic-acid-dimer); `ctest -R "gpu|rocm|eeq|gfnff"` 109 of 109 (no ctest covers the new kernels on ROCm). The lone-pair case 4 was exercised only with a hand-built pyridine plus two waters (218 HBs); that case-4 entries are present is inferred from the detection rule.
Charged systems: besides cl2-anion, `clusters/hydronium-water2.xyz` (+1) also differs between ROCm and CPU after 100 steps (Epot row 0.17 Eh, max |dx| 1.21 A), identical for the binary before this change.

### F11 result: cheaper HB/XB list refresh (branch `feature/rocm-hb-detection`, Oct 13, 2026)
Host cost per refresh on urea400-water1000 (subagent, medians, n = 11 refreshes each): 226.9 -> 52.2 ms; peo201x2-water1500 23.2 -> 12.6 ms. Split before (5 refreshes): detection ~100 ms (21.0M hydrogen visits, 3.1M pass geometry, 582k kept), copy into the CPU workspace 12, `rebuildBondHBData` 16, freeing the old list 11, HIP fill 52 (pack 21, permute 12, upload 3.5, slot tables 9), debug mirror 13-30.
Kept, all exact: (1) detection (shared CPU code) caches the hydrogens per acceptor and decides the case-1 distance and strength tests on plain doubles first, borderline cases (margins 1e-9 and 1e-6) go to the original expression; old list entries are reused as slots and filled in parallel; (2) `rebuildBondHBData` uses a per-hydrogen key table instead of a hash map (same order); (3) ROCm only, trait `lean_hb_refresh`: `setSkipHostHBWorkspaceUpdate` (the wrapper never runs the CPU HB workspace), parallel pack in device order, staging buffers kept, parallel stable counting sort for the slot tables, debug mirror only at verbosity >= 3.
Correctness: HB/XB list dumps (hexfloat) identical at 3 refreshes on urea, peo201-water500, peo201x2 for `-gpu none` and `rocm` (subagent, 18/18); 11 systems `-sp` energies identical to the previous binary, ROCm 100-step averages differ by <= 4.5e-13 Eh (run-to-run noise of the previous binary 9.1e-13). Author: last frame of a 60-step MD on urea and peo201x2 identical to the previous binary (max |dx| = 0, n = 2). The CPU path is bit-identical; `refset_regression.py` was not run (MOR41/GMTKN55 not on this machine; all are below 800 atoms and do not reach the changed ROCm path, but the shared detection code runs for them).
MD time per step ((t120 - t20)/100, 3 alternating runs, author): urea400-water1000 71.8 -> 54.7 ms (CPU 172.9), peo201x2-water1500 46.5 -> 45.7 ms (CPU 148.6). Remaining refresh cost on urea: ~30 ms detection, ~15 fill, 6.6 `rebuildBondHBData`; idle gaps of ~12.6 ms/step remain, cause not analysed.

### F12 result: FP32 matvec in the device EEQ PCG (branch `feature/rocm-eeq-fp32-matvec`, Oct 13, 2026)
Measured before (subagent, rocprofv3, 7320 atoms, 40 steps, n = 1): 17.5 `rocblas_dsymv` per step at 0.84 ms (14.6 ms, bandwidth-bound at 240-264 GB/s, `ssymv` takes half), matrix build 5.2 ms (compute-bound by the FP64 far-field sqrt and division).
Kept: the matrix kernel writes the FP64 lower triangle plus an FP32 copy (`d_A_f32`, +N^2*4 bytes, 214 MB at 7320) and uses `rsqrt` for the far field in that mode (max relative difference 4e-16 in the subagent's check); the CG product A*p uses `rocblas_ssymv`, vectors and dot products stay FP64; reliable update replaces r by the true FP64 residual when |r| fell by 1e-4 (`CURCUMA_ROCM_EEQ_FP32_DELTA`) and before convergence is accepted; convergence tolerance unchanged but tested on the true residual. `CURCUMA_ROCM_EEQ_FP32=0` restores the FP64 path (bit-identical to the previous binary in the subagent's runs, n = 13 systems). FP32 switches itself off if its allocation fails or N^2 exceeds the int range.
Effect at 7320 atoms: 4 FP64 + 16 FP32 matvecs per step instead of 17-18 FP64; EEQ GPU time 20.3 -> 13.3 ms, build 5.2 -> 2.9 ms (subagent).
MD time per step ((t120 - t20)/100, 3 alternating runs, author): urea400-water1000 54.6 -> 50.6 ms (CPU 173.5), peo201x2-water1500 45.5 -> 38.3 ms (CPU 145.5). Subagent table (3 runs each): 4600 atoms 31.4 -> 29.3, 3000 atoms 14.0 -> 13.1, 2910 atoms 14.2 -> 12.9, 1410 atoms 6.05 -> 5.93, 999 atoms 3.74 -> 3.71 ms.
Accuracy: author, 100-step MD vs previous binary on urea and peo201x2: last frame identical (max |dx| = 0), average Epot differs 1.5e-12 / 1.9e-12 Eh, Etot 3.1e-12 / 1.5e-12 Eh (n = 2). Subagent: per-step energy max |dE| 1.0e-10 Eh (n = 51 steps, 5 systems, one last-digit flip), charges max |dq| 7.1e-11 (n = 100 calls per system). Not exercised: cl2-anion, hydronium-water2, trisaccharide, caffeine never reach the device PPCG; systems above 7320 atoms, memory limit, `-opt`, long MD, other GPUs.

### F13 result: small-system step latency (branch `feature/rocm-launch-latency`, Oct 13, 2026)
Timeline of the 66-atom step (subagent, rocprofv3, 250 steps): 1.55 ms period, 69 GPU operations, GPU busy 0.78 ms, host waiting 0.99 ms in 5.2 stream syncs. Items: rocsolver potf2 110-130 us and potrs 65-90 us plus sync and host Schur step; `k_coulomb_tiles_hip` 147 us on a single workgroup; two displacement checks with pageable copies and a sync (~150 us); a dead dc6dcn skip-flag chain (5 ops, `shouldSkipDc6dcn()` has no caller), 3 unused CN downloads, 4 memsets; ~20 us per hardware-queue switch; forced HB refresh every 10 steps with 24-29 blocking pageable uploads (~0.9 ms).
Kept (switch back to the old path in brackets): (1) device EEQ for N <= 128 and <= 8 fragments, one workgroup (Cholesky in LDS plus Schur step) on the auxiliary stream (`CURCUMA_ROCM_EEQ_SMALL=0`); (2) Coulomb tile shape 32x4 / 32x8 / 32x16 below 2000 atoms, 128x128 above unchanged (`CURCUMA_ROCM_COUL_TILE=128x128`; 66 atoms 141 -> 27 us); (3) lean step bookkeeping: one launch for both displacement checks, CN downloads only when the host EEQ needs them, dead skip check removed, one zeroing kernel, pinned charge upload (`CURCUMA_ROCM_LEAN_STEP=0`); (4) single-stream mode below 500 atoms (`CURCUMA_ROCM_SINGLE_STREAM=0/1`; slower at 600 and 1002 atoms); (5) unchanged HB lists not re-uploaded, rest through a pinned async buffer (single-stream only); (6) dispersion launched early in single-stream mode (`has_stream_overlap` is true on ROCm now). No net gain: EEQ on its own stream overlapping phase 1 (66 atoms 0.67 -> 0.83 ms), hipGraph (micro-benchmark about 0.5 us saved per kernel, not pursued), EEQ on the null stream (fragile).
MD time per step (ms, median of 3 alternating runs, author; base -> new, CPU): caffeine 24 atoms 1.176 -> 0.439 (0.075); 66 atoms 1.526 -> 0.571 (0.271); 114 atoms 1.965 -> 0.863 (0.434); water300 2.076 -> 1.150 (1.244); water600 2.747 -> 2.438 (2.087); water1002 3.412 -> 3.003 (4.115). Subagent: urea400-water1000 50.74 -> 50.45, peo201x2-water1500 38.80 -> 38.29 ms.
Break-even moves from ~900 to between 300 and 600 atoms for water cuts (water300 GPU 8 % faster, water600 CPU 14 % faster); the `-gpu auto` threshold (`gpu_auto_min_atoms`, 900) is not changed here; other system types were not measured.
Accuracy (author, 7 systems incl. +1 and -1, 100-step MD vs previous binary): last frame identical (max |dx| = 0), average Epot differs <= 4.3e-14 Eh, Etot <= 5.7e-14 Eh. Subagent: `-sp` at 10 decimals identical on 12 systems. Not verified: CUDA build, gfx1150, threaded/multi-instance runs (stream handles are per thread), LU fallback of the small device EEQ for an indefinite matrix, optimisation and energy-only calls beyond ctest.

## 2. GFN2 (native xTB, `-gpu rocm`)

114 atoms (frame 1 of c33h45n9-agal-frames44), `-sp`, `-verbosity 3`, ROCm vs CPU:

| Phase | CPU | ROCm |
|---|---|---|
| setup | 11 ms | 96 ms |
| SCF, 14 iterations | 92 ms (eigen 74, 5.3/it) | 420 ms (eigen 392, 28/it) |
| post-SCF energy | 9 ms | 31 ms (D4 ATM 27.6 ms on the device, 207 triples) |
| TOTAL | 113 ms | 572 ms |

Iteration 0 took 113 ms and iterations 1-2 about 21-23 ms: ~90 ms one-time warm-up (library/kernel loading), plus ~85 ms extra setup.
Neither matters in MD/optimisation; both matter for many short jobs (`-confscan`, `-sp` batches). Cause not measured.

### G1. The per-iteration solve repeats work that is constant per geometry (measured, medium)
`XtbHipContext::residentSolve*` copies S0 and calls `rocsolver_dsygvd` every iteration (`xtb_hip_context.hip:1065-1069`, `1354-1358`). That
is potrf + sygst + syevd + trsm on a fixed S. CUDA caches L and reduces with triangular solves (`xtb_gpu_context.cu:3363`).
Micro benchmark, same matrices, per iteration:

| n (AOs) | dsygvd now | 2 dtrsm + dsyevd + dtrsm | ssygvd now | same in FP32 | dsygst alone |
|---|---|---|---|---|---|
| 256 | 23.0 ms | 11.6 | 17.2 | 7.7 | 9.3 |
| 512 | 48.5 | 29.5 | 34.2 | 16.4 | 19.2 |
| 1024 | 145 | 107 | 77 | 38 | 42 |
| 2048 | 452 | 361 | 194 | 106 | 110 |

So the cached-L reduction is 1.25x (n=2048) to 2x (n=256) faster in FP64 and 1.8-2.2x in FP32, in the matrices of the benchmark. The same
frame has nao ~300 (inferred from 28 ms/it, not printed), so est. 28 -> ~14 ms per iteration: still 2.6x slower than the CPU's 5.3 ms.
Not checked: that the reduced problem reproduces the `dsygvd` eigenvectors to the SCF tolerance in a real SCF (CUDA does it, 12 digits
claimed in its tests).

### G2. The eigensolver itself is the floor, and rocSOLVER is slow at small n (measured)
`rocsolver_dsyevd` alone: 4.9 / 10.7 / 17.5 / 25.6 ms at n = 128 / 256 / 384 / 512, i.e. time grows slowly with n: it is
latency-bound (kernel launches, host synchronisation), not FLOP-bound. The CPU does the same job in ~5 ms (nao~300, 10 threads, in-process
timer). `rocsolver_dsyevj` is 2-6x slower, `rocsolver_dsyevdx` for the lowest n/4 pairs is 2-11x slower than the full `dsyevd`: neither
helps, so CUDA's partial diagonalisation (`syevdx`) must not be copied.
Consequence: for nao below roughly 1000 the GPU cannot beat the CPU eigensolve with rocSOLVER, whatever G1 does.

### G3. Avoid the eigensolve: GEMM-only steps (est., large lever, risk medium)
Four GEMMs of size n cost 0.41 ms (dgemm) / 0.12 ms (sgemm) at n=256, 2.8 / 0.23 at 512, 20.5 / 1.25 at 1024, 160 / 9.3 at 2048: 30-100x
below one `dsyevd`. CUDA's `-scf_pseudo_diag` (first-order Jacobi rotation in the AO basis plus Cholesky re-orthonormalisation,
[GFN2_GPU_COST_PLAN.md](GFN2_GPU_COST_PLAN.md) stage 2/3) replaces the middle iterations that way. On CUDA it gained only 11.8 -> 10.6 s SCF
(polymer, 1 GPU) because cuSOLVER's eigensolve is already good; on ROCm the eigensolve/GEMM ratio is far larger, so the same idea should pay
more. It is opt-in and not bit-identical on CUDA (1.3e-4 bias fixed by a W correction, 2459/2460 GMTKN55 identical); a ROCm port
inherits that validation burden. No prototype; the SCF iteration count rising (CUDA saw 19 -> 24 before the correction) would eat part of the gain.

### G4. Hybrid: host eigensolver below a crossover (est., cheapest, risk low)
`XtbGpuAdapter` already has an `ExternalEigensolver` hook and the CPU solve is 2-5x faster than rocSOLVER for n < ~1000. Doing the
integrals, gradient and D4 on the GPU and the solve on the host (transfer of n^2 doubles: 0.5 MB at n=256, 8 MB at n=1024) would remove
G2 for small and medium systems. Break-even is not known: at 1410 atoms (nao ~4000, inferred) the whole `-sp` is 15.0 s ROCm vs 18.3 s CPU
(eigen 977 vs 1190 ms/it), so the GPU eigensolve starts to win somewhere between nao ~300 and ~4000. Needs the crossover measured with
the real hook before choosing a default. Also relevant for `-gpu auto`: today it selects the GPU even where it is 5x slower than the CPU.

### G5. Remaining GFN2 items (unanalysed)
- Device D4 ATM 27.6 ms for 207 triples is likely first-launch cost, not kernel cost (CPU post-SCF total 9 ms); measure a second call.
- Frequent blocking `hipMemcpy` plus `hipDeviceSynchronize` per solver call; at ~28 ms per iteration it is not the bottleneck now, but
  it will be after G1/G3.
- The GPU reports "low-power state" in `rocm-smi`; clock ramp-up on short jobs not measured.

## 3. Priority

| # | Change | Evidence | Est. effect | Effort | Risk |
|---|---|---|---|---|---|
| 1 | DONE (branch, uncommitted): far-field branch in `k_coulomb_dense`, `k_coulomb_gather` and EEQ matrix builds (threshold x >= 6) | measured 49 -> 8.5 ms, equal to 12 digits | gfnff 7320 atoms: -40 ms/step | small (one condition per kernel) | low; keep the pair kernels bit-comparable by testing against CPU as in `ROCM_EVALUATION` |
| 1b | DONE (branch, uncommitted): workgroup-per-atom `k_dispersion_gather` | profiler: 37.9 -> 11.2 ms | step -27 ms | small | low |
| 2 | DONE (branch `feature/rocm-eeq-ppcg`): device PPCG for the EEQ on ROCm (port of CUDA WP7-D) | CPU EEQ 100-120 ms is half the step | gfnff many-fragment systems ROCm faster than CPU (est.) | large | medium |
| 3 | G4: choose host eigensolver below a measured crossover; `-gpu auto` must not pick a 5x slower path | rocSOLVER syevd 10.7 ms vs CPU 4.5 ms at n~256-300 | gfn2 small/medium: ROCm time ~ CPU time instead of 5x | small | low |
| 4 | G1: cached-L reduction (2 trsm + syevd + trsm) instead of `dsygvd` | measured 1.25-2x on the solve | gfn2 large: modest | medium | low-medium (SCF equivalence) |
| 5 | G3: GEMM-only pseudo-diagonalisation | CUDA precedent, ROCm GEMM/eig ratio 30-100x | large for nao > 1000 | large | medium-high (validation) |
| - | Not worth it: `syevdx`, `syevj`, FP32 Coulomb, EEQ FP32 factor on GFN-FF before 1 and 2 | measured / reasoning above | - | - | - |

## 4. What to measure next (before coding)
1. Profiler (`rocprofv3`, installed Oct 13, 2026): done for GFN-FF at 7320 atoms; still to do for GFN2 and for the remaining kernels
   (`k_dc6dcn_per_pair` 8.1 ms, CN pair generation 12.8 ms, `k_coulomb_dense` 20 ms).
2. Re-measure with a statistic: current numbers are single runs on a machine where other jobs may run; use n >= 5 and report the median.
3. Time `-opt`/`-md` for gfn2 per SCF iteration on 3 systems (66, 114, 231 atoms) with the host-eigensolver hook to place the G4 crossover.
4. Run the same MD step breakdown on the CUDA hosts to confirm that the CUDA kernels listed here are the cause of its speed, not the GPU.

## Not covered
Multi-GPU, solvation, gradients of GFN2 on the device (only the `-sp` energy path was profiled), periodic systems, other AMD GPUs,
power/clock behaviour, memory footprint.
