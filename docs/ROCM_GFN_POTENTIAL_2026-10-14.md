# ROCm GFN1/GFN2: where the time goes and what could be gained (Oct 14, 2026)

Status: measurement and analysis of Oct 14 (AI-generated), followed on Oct 15 by the solver changes in
[Solver switches and results (Oct 15, 2026)](#solver-switches-and-results-oct-15-2026). Machine: RX 9060 XT (gfx1200, FP64
about 1/17 of FP32), Core Ultra 5 225, branch `feature/rocm-step-pauses-mirror`. The Oct 14 tables use `-threads 10`; all
wall times are single runs (n = 1) unless noted. The Oct 14 sections below describe the code before the Oct 15 changes.

## Solver switches and results (Oct 15, 2026)

AI-generated, machine-tested (ctest `gpu|rocm|sqm|xtb|grad`: 205 of 207 pass, `xtb_cpscf` and `cli_sqm_12` fail as before).
peo201-chain is gfn2 `-sp` on 1410 atoms (nao 3222), energy -2088.25340678 Eh in every run listed here.

| Switch | Default | What it does |
|---|---|---|
| `-xtb.rocm_cached_cholesky` | false (true until Oct 17) | S = L L^T once per geometry, each SCF step sygst, syevd, trsm instead of sygvd (FP64 and FP32); every new factor is probed (L L^T v = S v) and recomputed when it fails |
| `-xtb.rocm_density_syrk` | false | density P = C diag(occ) C^T by `rocblas_dsyrk` on the occupied columns, lower triangle mirrored |
| `-xtb.rocm_explicit_linv` | false (true Oct 15 to Oct 17) | with the cached factor (probed L (Linv v) = v): Linv = L^-1 once per geometry (`rocsolver_dtrtri`, FP64, FP32 copy), each step F' = Linv F Linv^T by two trmm, syevd, C = Linv^T Y by trmm. The speedups quoted below (-9.6 percent) belong to this switch and to the cached factor, both off by default since Oct 17 |
| `-xtb.rocm_device_multipole` | false | gfn2 only: multipole potential and energy of the host SCF loop on the device (FP64), see [Device multipole potential](#device-multipole-potential-oct-17-2026) |

Defaults live in the PARAM block of `native_xtb_method.h`, in `cfgBool(...)` in `xtb_hip_method.cpp` and in the member
initialisers of `XtbHipContext::Impl` (`xtb_hip_context.hip`). Other backends ignore these switches.

Measured (each line: structure, sample size, what was not checked):

- **Density GEMM over the occupied columns** (always on): peo201-chain 15.1 -> 14.3 s wall (n = 2), energy equal to 8
  digits; water1002 and caffeine energies unchanged.
- **Cached Cholesky / syrk**: peo201-chain median wall 14.35 s before, 14.12 cached factor, 13.81 syrk, 13.58 both (n = 4
  each); gfn1 17.41 -> 16.36 s (n = 3). The cached factor saves only potrf (about 15 ms per FP32 and 30 ms per FP64 step),
  syrk cuts the density from 94 to 49 ms per step. Max |dE| against the old path 1.5e-11 Eh, against the CPU 8.1e-11 Eh
  (6 cases).
- **Explicit L^-1**: peo201-chain median 13.71 -> 12.39 s (n = 4, alternating on/off); rocprof FP32 solve 48,281 -> 22,818
  kernels and 520 -> 322 ms span, FP64 span 1422 -> 1279 ms. |dE| against sygst/trsm <= 1.7e-11 Eh over 14 cases (7
  structures, gfn1 and gfn2: caffeine, caffeine+ doublet, Cl2- doublet, MOR41 PR22, water1002, peo201-chain, caffeine
  compressed to 0.6). cond2(S) <= 56 over 234 screened structures (MOR41, C60ISO, HEAVY28, DC13, UPU23, AL2X6, HEAVYSB11
  and local test structures, gfn2); caffeine compressed to 0.5 reaches 1.1e3 and 1.8e3 (the two methods' bases), |dE| <= 3.4e-12 Eh.
  Stress test: caffeine, an 8-water cluster and a Rh complex scaled about the centroid from 1.0 to 0.10 reach cond2(S)
  4.5e6; over 40 valid geometries x gfn1/gfn2 x 7 solver paths max |dE| against the CPU is 6.0e-11 Eh (L^-1), 4.9e-11
  (sygst), 5.7e-11 (sygvd). At the CPU, the gfn2 D4 gradient of these compressed geometries disagreed with finite
  differences (up to 6.1 Eh/A, caffeine scaled 0.40 to 0.55); that was a host-side D4 weight-derivative underflow, fixed
  for every backend (`AIChangelog.md`). At n = 6222 the FP64
  trmm reduction is slower than sygst in a standalone benchmark (1139 against 993 ms); peo201-water500 (2910 atoms) off
  72.82 s / on 71.87 s (n = 1), peak VRAM 8.0 of 17.1 GB.
- **All 8 switch combinations** on peo201-chain (n = 2 each): max |dE| 1.5e-11 Eh, median wall 14.45 s (all off) to
  12.31 s (all on); the new default (cached factor + L^-1, dgemm density) 12.84 s. Re-measured after the default change:
  12.85 s median (12.94 / 12.85 / 12.83 s, n = 3), `-xtb.rocm_explicit_linv false` 14.18 s (n = 3, run after the default
  runs, not alternating).
- **Dynamics and optimisation**: 100-step MD of caffeine and a 114-atom molecule, max |dEpot| against the CPU 1.9e-7 Eh,
  the same as the sygst path; trtri runs once per geometry (101 calls in 100 MD steps). gfn1 NVE MD, 100 steps, against
  the CPU: caffeine max |dEpot| 1.2e-7 Eh, the 114-atom molecule 4.8e-7 Eh (n = 1 each); gfn2 peo201-chain (1410 atoms)
  30 steps max |dEpot| 1.8e-6 Eh, peak VRAM 2.63 GB (n = 1). Caffeine `-opt` 33 steps, same
  final energy with and without L^-1. 114-atom `-opt` (n = 1 per path): final energies of all paths and the CPU within
  1.4e-6 Eh; the 1e-7 Eh agreement criterion fails for every path including the sygst path (different optimiser
  trajectories); wall gfn2 38.0 s (L^-1), 58.0 s (cached factor), 60.3 s (sygvd), 22.2 s (CPU), cause of the L^-1 gap not
  examined.
- **GFN2 multipole potential and energy on the host** (`xtb_multipole.cpp`, all backends and the CPU): reading contiguous
  columns cut the host part of one SCF iteration on peo201-chain (1 host thread, 12 iterations, n = 1) from about 112 to
  44 ms; wall peo201-chain 15.21 -> 14.23 s (n = 5), water1002 ROCm 4.53 -> 4.28 s (n = 5), CPU 4.56 -> 4.32 s (n = 1).
  The results differ from the old code in the last bit of `v_qp`: the old strided loop was auto-vectorised without FMA
  (GCC 16, `-O3 -mavx2 -mfma`), the new loops contract every product into FMA; with `-ffp-contract=off` old and new agree
  bit for bit. Converged energies of caffeine, water1002 and peo201 agree to 1e-10 Eh. Gradients and CPSCF only through
  ctest.
- **GPU gradient with open shells** (`xtb_gradient.cpp`, `calculateGradientGpu`): the device energy-weighted density was
  built with occupation 2 for every column up to round(nocc/2), so a doublet SOMO and fractional occupations were
  weighted wrongly. Now the host passes eps_w = 0.5 f eps and the last occupied column. GPU against CPU max |dG| for
  caffeine+ doublet gfn2 6.4e-2 -> 7.4e-7 Eh/A (6 open-shell or charged cases, gfn1 and gfn2, n = 1 each); finite
  differences (11 coordinates) agree with the CPU; closed shell unchanged. Only ROCm was run: CUDA and Vulkan receive
  the same host input but were not compiled or run here, and the CUDA `w_from_goo` branch (`xtb_gpu_context.cu`) still
  assumes closed shell.

### Fallbacks of the solver and their notice

The resident solver keeps running on the device when a switch cannot be honoured: L^-1 allocation fails (FP64 or FP32
buffer, then sygst/trsm for the rest of the run), trtri fails or reports a singular L (sygst/trsm for that geometry, the
next geometry tries again), or the cached factor cannot be formed (allocation, or potrf reports S not positive definite;
sygvd for that geometry, an allocation failure for the rest of the run). The context records each kind as a bit
(`XtbHipContext::solverFallbackEvents()`); `XtbHipComputationalMethod` prints one `[WARN]` line per kind and method
object, not per geometry. While an optimiser or MD driver runs the energy calls at verbosity 0 the notice stays pending
and is printed at level 1 when the method object is destroyed (as the GPU fallback summary does). These are GPU-internal
path changes, not CPU fallbacks, so they are not counted in `gpu_fallback` and `-gpu_strict` does not stop on them.

Diagnostic: `CURCUMA_ROCM_FAIL_LINV=alloc|trtri|chol` (read once per process) forces the corresponding failure; the notice
then ends with "[forced by CURCUMA_ROCM_FAIL_LINV]". Check on caffeine (n = 1 run per value): gfn2 and gfn1 `-sp`, each
value prints exactly one notice, |dE| against the unset run <= 2.1e-13 Eh; gfn2 `-opt` (33 steps, 41 energy calls) prints
one notice per value and reaches the same final energy. Real allocation or trtri failures were not provoked; the
injected failures skip the call (alloc, chol) or overwrite its result (trtri).

### Device multipole potential (Oct 17, 2026)

AI-generated, machine-tested. With `-xtb.rocm_device_multipole true` the GFN2 multipole potential (`v_at` shift, `v_dp`,
`v_qp`) and the multipole energy of each SCF iteration run on the device (`k_mp_potential_stored`, `k_mp_energy_stored` in
`xtb_hip_context.hip`, one block per atom, shared-memory reduction). The host still builds the interaction matrices in
`setupMultipole()`; the 15 distinct ones (`amat_dd` is symmetric in its two indices) are uploaded once per
`setupMultipole()` call, 15 nat^2 doubles (239 MB at 1410 atoms). Per iteration the host sends q_at, dp_at, qp_at and
receives v_at, v_dp, v_qp (10 and 10 nat doubles); `v_dp`/`v_qp` stay on the device for the Fock build. Above 4 GB of
matrices (about 5770 atoms), or when an allocation, a kernel or the solve without upload fails, the run continues on the
host path and one `[WARN]` line is printed (same reporting rule as the solver fallbacks above). The multipole gradient
stays on the host. With the switch off the host code path is unchanged.

Measured against the switch off (gfn2 `-sp`, `-gpu rocm`, energy and `-dump_gradient`):

| Structures | n | max dE (Eh) | max dG (Eh/A) |
|---|---|---|---|
| caffeine, caffeine +1 doublet, MOR41 AcCl AcI Bz CO ED01 ED07 ED16a ED24 PR07 PR16 PR25 PR41 | 14 (1 run each) | 4.3e-14 | 8.0e-14 |
| peo201-chain (1410 atoms), sygvd and cached factor + L^-1 | 4 pairs | 4.5e-13 | 1.3e-14 |

Against the CPU (`-gpu none`, same 15 structures, n = 1 each) the switch on and off show the same differences (max dE
5.5e-12 Eh, max dG 2.7e-6 Eh/A, both MOR41 ED01; these come from the rest of the ROCm path, not from this switch).

peo201-chain `-sp` at the default 1 host thread, 12 SCF iterations, runs alternating off/on (n = 3 each, median):
wall 14.92 -> 14.38 s with the default sygvd solver, 13.35 -> 12.86 s with `-xtb.rocm_cached_cholesky true
-xtb.rocm_explicit_linv true`. Per iteration (verbosity 3 timers): potential build 20.6 -> 6.9 ms (of which D4 5.8 and
5.0 ms), multipole energy 14.75 -> 0.91 ms. The eigensolve (780 to 900 ms per iteration) stays the dominant part.
Not run: optimisation or MD (re-upload per geometry), implicit solvent, open-shell MOR41-type metals beyond caffeine +1,
the fallback paths (not provoked), peak VRAM.

### Open levers (Oct 15)

- **Host threads**: the ROCm adapter passes `-threads` as given, so the default host work runs on 1 thread
  (`xtb_gpu_adapter.h`); the CUDA wrapper defaults to min(hardware threads, 8) (`xtb_gpu_method.cpp`). peo201-chain with
  `-threads 8`: 14.2 -> 13.78 s (n = 4, before the L^-1 default); not made the default.
- **Multipole on the device**: opt-in since Oct 17 (`-xtb.rocm_device_multipole`, see above); by default potential and
  energy run on the host (about 29 of the 44 ms host part per iteration at 1 thread; 4.0 + 6.7 ms at 8 threads).
- **Remaining host gap** per SCF iteration (D4 potential 5.9 ms serial, populations and moments 7 ms, Coulomb about 2 ms
  at 1 thread, peo201-chain, n = 1 run).
- **D4 ATM 3-body** about 1.2 s per peo201-chain run (one FP64 kernel, `k_d4_atm`), unchanged.
- **Pseudo-diagonalisation** on ROCm (lever 3 below), not started.
- Later FP64 switch (lever 2 below) and the `-gpu auto` threshold for xTB (lever 6) are unchanged.

## CPU against ROCm, `-sp`, same energies to 8 decimals

| System | Atoms | gfn2 CPU / ROCm | gfn1 CPU / ROCm |
|---|---|---|---|
| caffeine | 24 | 0.02 / 0.33 s | 0.02 / 0.31 s |
| urea-water-cluster197 | 197 | 0.44 / 1.55 s | 0.31 / 1.02 s |
| pep20 | 309 | 0.80 / 1.81 s | 0.75 / 1.54 s |
| water300 | 300 | 0.35 / 0.93 s | 0.40 / 1.02 s |
| water1002 | 1002 | 4.51 / 4.42 s | 5.90 / 7.17 s |
| peo201-chain | 1410 | 18.98 / 15.28 s | 24.78 / 19.06 s |

`-gpu auto` selects ROCm for GFN1/GFN2 at every size (checked on caffeine: "using native xTB on GPU (rocm)"); the
`gfnff.gpu_auto_min_atoms` threshold covers GFN-FF only. Break-even is near 1000 atoms for gfn2, above 1000 for gfn1.
ROCm `-opt` of a 114-atom and a 66-atom molecule is 3x slower than the CPU (n = 3, `ROCM_EVALUATION_2026-10-13.md`).

## Where the ROCm time goes (gfn2, peo201-chain, 1410 atoms, nao 3222, 15.0 s)

SCF 12.4 s (82 %), of which the "solve eigen" timer shows 11.7 s: 7 FP32 iterations at 0.64 s and 5 FP64 iterations at
1.55 s (switch at `scf_fp32_threshold` 1e-3); the solver itself is 9.7 s of that (7 x 428 + 5 x 1337 ms), the remaining
2.0 s is the FP64 density GEMM (see below). Setup 1.2 s (`k_overlap_h0` 0.46 s), D4 ATM 3-body 1.2 s (one FP64 kernel,
`k_d4_atm`). peo201-water500 (2910 atoms, nao 6222, 79 s): "solve eigen" timer 86 % (FP32 3.4 s and FP64 9.4 s per iteration, density GEMM included; at standalone dgemm speed 1.1 s per iteration it would be about 13 s of the 68 s, estimate).
rocprofv3 on the 1410-atom run: GPU busy 13.0 of about 15 s; FP64 GEMM kernels (inside sygst/syevd) 5.5 s,
`latrd` (tridiagonalisation, gemv) 1.8 s.

## Standalone rocSOLVER timing (subagent, median of 3, ms, FP32 / FP64)

| n | sygvd | syevd | sygst | potrf + trsm | gemm TFLOPS |
|---|---|---|---|---|---|
| 3222 | 410 / 1327 | 240 / 966 | 148 / 247 | 23 / 149 | 6.8 / 0.43 |
| 6222 | 1992 / 7998 | 1534 / 6212 | 365 / 994 | 88 / 885 | 7.2 / 0.43 |

In curcuma the solver call is as fast as standalone: temporary timers in `resSolveMP` (not committed) gave 416 to 435 ms
(ssygvd, hipEvent 417 to 426 ms, n = 7 calls) and 1330 to 1396 ms (dsygvd, n = 5 calls) at n = 3222, so the earlier
impression of "1.2 to 1.7 times slower" came from comparing the SCF iteration time (640 / 1545 ms) and the
"solve eigen" timer of `-verbosity 3` (11.7 s = 7 x 428 + 5 x 1337 ms + about 2.0 s) with the standalone number.
The difference is the density matrix: `resDensity` computes P = (C diag(occ)) C^T with a full n x n x n `rocblas_dgemm`
(`xtb_hip_context.hip`, `resDensity`), although columns beyond `ncol` are zero. rocprofv3 shows one FP64 GEMM kernel of
165 ms per iteration (n = 3222; 0.43 TFLOPS peak FP64 gives 2 n^3 / 0.43e12 = 160 ms), and thread count does not matter
(`-threads` 1 / 4 / 10: 710 / 655 / 660 ms per FP32 iteration, n = 1). Jacobi, QR and bisection variants are 10 to 50 times slower than
syevd (`syevdj` 80 s at n = 6222). FP32 purification (2 GEMM + geam per iteration, 17.3 ms at n = 3222, 120 ms at 6222)
only beats sygvd below about 10 to 16 iterations (arithmetic from the table, not a converging purification).

## Density GEMM (Oct 14 finding; done Oct 15, see above)

P is built with a full n^3 FP64 GEMM. Using only the `ncol` occupied columns and a symmetric rank-k update
(`rocblas_dsyrk` on C diag(sqrt(occ)), valid because occupations are non-negative) would cut the flops by about a factor 4
(E: 165 -> about 40 ms per iteration at n = 3222, 12 iterations: about 1.5 s of the 15 s run; at n = 6222 about 10 s of
79 s). Implemented on Oct 15 (occupied columns always, syrk as opt-in `-xtb.rocm_density_syrk`); results change at
rounding level only. In FP32 iterations the same product could still run in FP32 (not done).

## Levers (estimates marked E are arithmetic on the numbers above, not measured in curcuma)

0. **Density GEMM** (previous section). Done Oct 15.
1. **Precomputed L instead of potrf + sygst per iteration.** Done Oct 15: the cached factor alone saves only potrf
   (-1.6 % wall); the explicit L^-1 with trmm (default since Oct 15) gives -9.6 % (13.71 -> 12.39 s, n = 4). The Oct 14
   estimate was: F' = L^-1 F L^-T needs two trsm (7 ms each FP32, 85 ms FP64 at n = 3222). E: FP32 solve 410 -> 265 ms
   (-35 %), FP64 1327 -> 1190 ms (-10 %).
2. **Later FP64 switch.** Measured on peo201-chain: `-scf_fp32_threshold` 1e-3 / 1e-4 / 1e-5 / 1e-6 gives 15.1 / 13.3 / 12.2 /
   15.0 s wall (12 / 12 / 13 / 17 iterations), energy -2088.25340678 Eh in all four (n = 1 structure). 1e-6 stalls in
   FP32 and loses. `TODO.md` asks not to change `scf_fp32_stall_patience` until the fixed-point question is settled.
3. **Pseudo-diagonalisation on ROCm.** CUDA has `-scf_pseudo_diag` (113 -> 89 s on polymer_2x, 4 GPUs, MOR41 95/95
   identical, GMTKN55 2459/2460); ROCm has none. A pseudo step costs three GEMMs; E: FP32 about 40 ms against 640 ms,
   FP64 about 470 ms against 1550 ms at n = 3222. Largest single lever, largest effort.
4. **FP32 eigenvectors plus FP64 correction instead of the FP64 solves** (CUDA plan item 0.2, `GFN2_GPU_COST_PLAN.md`).
5. **D4 ATM 3-body (1.2 s) and integrals (1.2 s)**: FP64-bound single kernels; cutoff and tiling untested here.
6. **Small systems**: plugin start 0.2 s, per-call overhead, host-side potential and gradient parts
   (`ROCM_EVALUATION_2026-10-13.md`). Cheapest fix: `-gpu auto` uses the CPU for GFN1/GFN2 below about 1000 atoms (the
   GFN-FF rule applied to xTB); the threshold needs a measurement on more system types than water and peo.

## Not measured

MD or `-opt` per-step cost of GFN2 on ROCm (only whole optimisations, n = 3); gradient phases on the device; gfn1
eigensolve split; memory beyond the two peak-VRAM readings above; gfx1150; transition metals beyond MOR41 PR22 and the
Rh complex of the stress test (energies only); CPU large-system modes (`-large_system_mode`) as the
comparison at 1410 atoms (the CPU time above is the default mode).

## Wrong rocSOLVER results with several processes on one GPU (Oct 16/17, 2026)

AI-generated, machine-tested; measured on RX 9060 XT (gfx1200), ROCm 7.2.4, Linux 6.18.49, neutral caffeine, gfn2 `-sp`,
`-scf_threshold 1e-10`, 145 geometries (the FD set of the gradient check) times repeats, 4 processes at a time unless noted.

- **Symptom.** A run reports "SCF converged" and an energy that is off from the median of its geometry by 1e-8 to 9e-3 Eh,
  without a warning. Share of such runs: 14 of 20300 and 7 of 8700 (cached factor on, the default until Oct 17); 0 of 13050 with
  one process at a time (4350 + 8700), 1 of 8700 with two processes. CPU: 0 of 8700.
- **Bisection.** 0 of 8700 at `9ef9b98c`, `948d7c12` and `e2c72cc1`; 7 of 8700 at `ca2ce18e` (introduces the cached
  factor); FP32 steps without cache 10 of 8700, FP64 steps without cache 0 of 8700.
- **Cause as far as found.** An in-process comparison of the cached solve with `rocsolver_dsygvd` on the same F and S
  (4350 runs) found 146 runs where the eigenvalues differ by 8e-6 to 2.4e-3 Eh, constant over all SCF steps of the run.
  In 83 recorded cases the cached factor was wrong (max |L L^T - S| up to 3.3e-3, normal 4e-16) or the factor was correct
  and the `sygvd` result was wrong (generalized residual up to 3.5e-3, the cached result 1e-15). A program without curcuma
  (`scripts/rocsolver_concurrency_check.cpp`, n = 66) shows the same: four processes at a time, `rocsolver_dpotrf` wrong in 6
  and `rocsolver_dsygvd` in 8 of 12000 calls each; one process 0 of 12000; `rocblas_dgemm` 0 of 48000 (short run). Whether
  the fault is in rocSOLVER, the driver or the hardware was not tested; no report was filed.
- **Not the cause.** Host threads (`-threads 1`: 7 of 8700), shifted allocation addresses (0 of 200, one process),
  the FP32 phase (`-scf_fp32_threshold 1e3`: 7 of 8700), reuse of the factor from the device integral build (8 of 8700 with
  it disabled). Serialised kernels (`AMD_SERIALIZE_KERNEL=3 AMD_SERIALIZE_COPY=3`) left 2 of 8700.
- **What was changed (Oct 17).** Both cache switches default to false. A new factor, and a new Linv, is probed once per
  geometry with a fixed vector (`factorCheck` in `xtb_hip_context.hip`: |L L^T v - S v| and |L (Linv v) - v| below 1e-10
  relative); a failing factor is recomputed twice at most, then `sygvd` is used; a notice is printed once per run
  (`kFallbackFactorRetried`). `CURCUMA_ROCM_FAIL_LINV=factor|factorall` forces a failing probe for tests.
- **Limits.** The probe sees a wrong factor, not a wrong `sygvd` result (the default path), a wrong `syevd`, or a wrong
  `dgemm`. The default path and parallel jobs on one GPU therefore stay at risk; a GPU job per process is the safe use.
  CUDA was not tested.

