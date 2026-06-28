# SQM/FF Precision & BLAS-Threading Work Packages

Status: 🤖 AI-generated, machine-tested. Captures the June 2026 performance session
(BLAS thread starvation, mixed-precision SCF, GFN-FF EEQ, per-card levers).

## Background: the FP64 wall + the OMP=1 pin

Two orthogonal performance facts drive everything here:

1. **System OpenBLAS is the OpenMP build**, and curcuma pins OMP/MKL to 1 thread globally
   (`CxxThreadPool` `omp_set_num_threads(1)`, `main.cpp` `MKL_Set_Num_Threads(1)`) to avoid
   oversubscription in molecule-parallel batches. Side effect: a single large dense
   BLAS/LAPACK solve (EEQ Cholesky, xTB eigensolve) was starved to ONE thread. On the
   OpenMP OpenBLAS build `openblas_set_num_threads()` is a no-op — only `omp_set_num_threads()`
   takes effect. See `src/core/blas_threads.h`.
2. **GFN2's hot path is the FP64 generalized eigensolve** (O(nao^3) per SCF iteration).
   Consumer/workstation NVIDIA cards cripple FP64 (1/32 to 1/64 of FP32), so the eigensolve
   is FP64-bound and never saturates the card (100% "util", low power).

## Benchmark reference — polymer (1410 atoms), GFN2 single point

CPU scaling (release/, with the WP-1 eigensolve threading; energy −2088.25340678 throughout):

| threads | 1 | 2 | 4 | 8 | 16 |
|---------|-----|-----|-----|-----|-----|
| time    | 123 s | 76 s | 54 s | **48 s** | 70 s |

8 threads is the sweet spot; 16 regresses (SMT oversubscription — the eigensolve saturates ~8).

GPU (GTX 1660, working build, `-gpu cuda` → mixed precision ON by default; energy identical):
- **30 s total, SCF 23 s / 11 iters.** Per-iteration `t/ms`: iters 0–6 ≈ **1200 ms (FP32)**,
  iters 7–10 ≈ **3620 ms (FP64 polish)**. The 4 FP64 iterations = 14.5 s = **63 % of the SCF**.
- So even though the GPU wins here (30 s < 48 s CPU), the FP64 polish iterations are the
  bottleneck and the source of the low power draw. They are the tuning target.

(Caveat: an earlier benchmark used a stale `release_cuda` binary that ran ~50 s/iter with the
GPU idle — that was a broken build, not representative. The 30 s figure is the real one.)

## Hardware FP64 matrix + GFN2 verdict

| Card | arch | FP64 ratio | FP64 ~TFLOPS | GFN2 verdict |
|------|------|-----------|--------------|--------------|
| GT 1030 | Pascal | 1/32 | ~0.03 | useless for compute |
| GTX 1660 | Turing | 1/32 | ~0.16 | GPU ~1.6× CPU here, FP64-bound |
| RTX 5080 | Blackwell | 1/64 | ~0.9 | FP32-mixed only; high abs FP32 |
| RTX 4500 Ada | Ada | 1/64 | ~0.6 | FP32-mixed; 24 GB for large systems |
| RTX 5000 Ada | Ada | 1/64 | ~1.0 | FP32-mixed; 32 GB |
| A100 / H100 | DC | **1/2** | 9.7 / 34 | FP64 path flies — right tool |
| AMD MI250 / MI300 | CDNA | 1/2+ | 24 / 80+ | best FP64; ROCm path exists |

## Work Packages

### WP-1 ✅ DONE (committed) — BLAS thread starvation
- EEQ Schur-Cholesky: `ScopedBlasThreads` (RAII, driven by the EnergyCalculator `num_threads`)
  around the Phase-2 solve. `src/core/blas_threads.h`, `eeq_solver.cpp`. ~4.5× (mixture2 8 thr).
- xTB SCF eigensolve: `MklThreadScope` (`xtb_native.h`) extended to also drive per-thread
  `omp_set_num_threads`, so the GFN1/GFN2 eigensolve threads on OpenBLAS too. complex/231 8 thr:
  gfn1 2.69×, gfn2 2.41×. Budget = `effectiveIntraThreads` (1 in batch/suppressed).
- Also fixed: GFN-FF HB-coordination-number data race (separate, threaded-energy drift).
- Audited NOT beneficial: rf_solver/lbfgs (Lanczos for large → dsyevd only tiny), hessian/ancopt
  (Eigen's own non-BLAS solver, would need EIGEN_USE_LAPACKE first).

### WP-2 🟡 OPEN — per-card `scf_fp32_threshold` tuning
- Mixed precision is ON by default for `-gpu`; threshold `-scf_fp32_threshold` default 1e-3.
- Logic: `eig_fp32 = mixed && (dq_prev > threshold)` → **SMALLER threshold = stays FP32 longer
  = fewer expensive FP64 polish iterations** (faster on FP64-weak GPUs). Must stay above
  `scf_threshold` so ≥1 FP64 step converges, and above the FP32 accuracy floor (~1e-6) or the
  FP32 fixed point never crosses the threshold. (Help text was backwards — fixed `xtbinterface.h`.)
- TODO: sweep 1e-4 / 3e-4 / 1e-3 per card on a representative system; record the optimum.
  Expect lower-is-better on 1/32–1/64 FP64 cards, less critical on A100/H100/MI300.

### WP-3 🟡 OPEN — auto-default precision/threshold by device FP64 ratio
- Query the GPU FP64:FP32 ratio (compute capability / `cudaDeviceProp`) at init:
  1/32–1/64 (consumer/workstation) → mixed precision on + lower threshold; 1/2 (DC/MI300) →
  FP64 default. Removes manual per-card tuning. CPU path unchanged.

### WP-4 🟢 IMPLEMENTED (CUDA, Jun 2026) — GFN-FF mixed-precision iterative refinement (EEQ solve)
- GFN-FF has NO SCF/eigensolve, so the "FP32 early-iter" pattern does not map. The structural
  analog is the **EEQ linear solve** (dpotrf/dpotrs). Applied classic mixed-precision
  **iterative refinement** (LAPACK dsposv pattern): factor A in FP32 (cusolverDnSpotrf on an
  FP32 copy, d_A kept as the FP64 matrix), FP32 spotrs, then refine with the FP64 residual
  `b − A·x` (cublasDsymm) + an FP32 correction, 1–2 steps → full FP64 accuracy. Opt-in
  `-gfnff.eeq_mixed_precision` (+ `eeq_mixed_precision_iters`), default OFF on CUDA.
  `EEQSolverGPU::setMixedPrecision` / `mixedFactor` / `mixedSolveRefine` in
  `eeq_solver_gpu.cu`; wired into the factor-dominated paths (solve / solveWithDeviceRHS /
  GPU-Schur nfrag=1) with auto-fallback to FP64 dpotrf/LU when the FP32 factor is not SPD; the
  nfrag>1 general path stays FP64 (nrhs≈N makes the FP64 residual GEMM as costly as the solve).
  **Validated (GTX 1660, 1/32 FP64): triose nfrag=1 SPD 41-step MD bit-identical to FP64.**
  Honest: the factor-dominated wall-time win needs a LARGE nfrag=1 SPD system (the available GPU
  GFN-FF set lacks one: triose tiny, polymer indefinite→LU, mixtures multi-frag); win expected
  but not yet measured. **ROCm mirror DONE (Jun 2026)**: `mixedFactorHip`/`mixedSolveRefineHip`
  + `k_cast_d2f_hip`/`k_cast_f2d_hip`/`k_axpy_f2d_hip` in `rocm/eeq_solver_hip.hiph` (rocSOLVER
  `spotrf`/`spotrs` + rocBLAS `dsymm`, same hand-rolled SPD iterative refinement), gated to
  nfrag<=1 in `eeqBuildFactorSolve`, opt-in via `-gfnff.eeq_mixed_precision` (default OFF on ROCm
  too); Radeon 890M caffeine + 231-atom `complex` energy and caffeine opt bit-identical to the
  FP64 ROCm path. Pairwise FF energy/gradient terms remain a poorer candidate (FP32 risks the
  1e-6 Eh goal).

### WP-5 🟢 IDEAS — bigger GPU levers (newer cards / HPC)
- **Tensor-core GEMM (TF32/FP16)** for the non-eigensolve BLAS (Fock build, density C·Cᵀ,
  Löwdin) — 10–50× FP32 rate on Ada/Blackwell/Hopper; the eigensolve itself does not use them.
- **GPU O(N) density purification** (`large_system_mode sparse|dc` is CPU-only today) — avoids
  the O(nao^3) FP64 eigensolve for 1000+ atoms entirely.
- **Multi-GPU fragment SCF** (`large_system_mode fragments`) — one fragment per card on the HPC.
- **Reduce per-iteration host syncs** (documented dominant cost beyond the eigensolve).
- **Investigate the GPU-build single-threaded CPU setup** (observed on the stale build; confirm
  on a current build — the GPU path should still thread / offload the CPU-side setup).

## Practical routing today
- GFN2 daily: CPU `-threads 8` (8 is optimal; avoid 16). On a working GPU build the GTX 1660 is
  ~1.6× faster but FP64-bound.
- 5080 / Ada: only with FP32-mixed (+ future tensor-GEMM); mainly for VRAM on large systems.
- A100/H100 or MI300: the real FP64 GPU path. GT 1030: never for compute.
