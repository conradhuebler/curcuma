# GFN-FF GPU — cross-backend work packages

> Claude Generated (June 2026). Spun out of the ROCm GFN-FF port + a code-side audit comparing
> the CPU / CUDA / ROCm GFN-FF paths. These are the items that **cannot be implemented or validated
> on the available ROCm box** (no NVIDIA / no CDNA / Vulkan GFN-FF absent) or that are optional
> enhancements deferred from that session. Everything ROCm-testable from the audit was already done
> (true-D4 verification + doc fix, dense-Coulomb cutoff_active gate).

## Background — what the ROCm port revealed

The GFN-FF MD step is dominated by two O(N²) pairwise gradient kernels (Coulomb ~993k pairs no
cutoff; dispersion ~943k pairs). The CUDA path computes their gradient with per-pair `atomicAdd`
into `grad`/`dEdcn`. On RDNA this contention (~6–8M FP64 atomics into N addresses) was the bottleneck;
rewriting both as **per-atom GATHER** kernels (`k_coulomb_gather` / `k_dispersion_gather`, CSR
adjacency; dynamic dc6dcn read via a pair-index + is_i flag) took polymer/1410 MD from slower-than-CPU
to ~2.3× faster. FP64 atomics themselves are fine on RDNA (the gfn2 ROCm path uses them) — the issue
was *contention at ~1M pairs*, which the gather removes.

## WP-A — Port the gradient gather to CUDA (large-N MD)  [needs NVIDIA to validate]

CUDA optimized the **energy** atomics (Phase-5 `blockReduceAddEnergy`) but the **gradient** still
does per-pair `add_grad` + `atomicAdd(&dEdcn[...])`. At ~1M pairs even NVIDIA's fast FP64 atomics pay
a contention cost the CUDA GFN-FF was apparently never profiled for. Port `k_coulomb_gather` /
`k_dispersion_gather` (+ the CSR build in the SoA upload) into `cuda/gfnff_kernels.cu` +
`ff_workspace_gpu.cu`, gated on a size threshold (keep the pair kernels for small N where the gather's
per-atom serialism doesn't pay). The ROCm kernels in `rocm/gfnff_rocm.hip` are a 1:1 blueprint.
*Cannot be measured here (no NVIDIA GPU).*

## WP-B — Wave-size-agnostic warp reduction (CDNA / wave64)  [needs MI hardware to validate]

`blockReduceAddEnergy` / `warpReduceSum` hard-code 32-lane warps (`& 31`, `>> 5`, `warp_sums[32]`,
`offset = 16`). Correct on gfx11xx/RDNA (HIP compute = wave32), **wrong on CDNA (MI200/MI300, wave64)**.
Make the reduction use `warpSize` / `__AMDGCN_WAVEFRONT_SIZE` (loop start `warpSize/2`, lane mask
`& (warpSize-1)`, `warp_sums[64]`, `num_warps` from `warpSize`). Affects both the hipified energy
kernels and the new gather kernels. *Cannot be validated here (gfx1150 is wave32 only).*

## WP-C — Device-resident EEQ Schur for ROCm  [ROCm-testable; deferred perf]

ROCm currently does `rocsolver_dpotrf/dpotrs` (+ LU fallback) → z1/Z2, then the host applies the exact
CPU Schur complement (a D2H of z1/Z2 + CPU work + H2D per step). Exact and cheap for single molecules,
but for multi-fragment / very large MD the round-trip adds up. CUDA keeps the charges on-device via
`k_eeq_reduce_sums` / `k_eeq_schur_*` / PCG (cuSOLVER). Port those kernels to HIP (rocSOLVER/rocBLAS)
and wire `EEQSolverHip::solveWithDeviceRHSAndGPUSchur*` to return `true` + `getDeviceChargesPtr()`.
Optional — the current path is correct, this is residency/perf only.

## WP-D — Vulkan GFN-FF  [Vulkan GFN-FF does not exist yet]

Only native xTB runs on Vulkan today; there is no Vulkan GFN-FF. Vulkan has **no FP64 atomics**, so
the per-atom gather is mandatory — the ROCm `k_*_gather` kernels are the direct SPIR-V blueprint (the
xTB Vulkan gradient already uses the same per-atom-gather pattern). A port would mirror the ROCm TU:
SPIR-V kernels for the bonded/nonbonded/Coulomb/dispersion terms + a host eigensolve-free EEQ (or the
existing Vulkan FP64 solver) + CPU Schur.

## WP-E — Re-enable FP32 mixed precision  [ROCm-testable; small]

`use_mixed_precision` defaults to `false` in both CUDA and ROCm (`// TODO: restore to true after GPU
gradient debugging`). The `k_repulsion_mixed` / `k_batm_mixed` / `k_xbonds_mixed` kernels exist but are
unused. Re-enabling (with validation vs the FP64 path) would help RDNA (weak FP64) more than NVIDIA.
Gate it behind `-scf_mixed_precision` for GFN-FF as the xTB ROCm path does.

## Also worth a cleanup (CPU-side, not GPU)

The standalone external D4 (`dftd4interface` / `curcuma_d4`, which `USE_D4` builds and which hard-links
LAPACKE) is **not used by GFN-FF** at all — GFN-FF has its own self-contained `D4ParameterGenerator`.
The `USE_D4` CMake block / the standalone `-d4` method are the only consumers. If the standalone `-d4`
path is itself legacy, the LAPACKE dependency could be dropped entirely.

## EEQ for many fragments (solvent) — Schur / batched / PCG analysis

Fragments = connected components of the bond graph (BFS, `gfnff_method.cpp:8583`), so each solvent
molecule is one fragment → nfrag large, N large. The charge constraint is one bordered row/column
per fragment: `A_aug = [[Coulomb (N×N), Cᵀ];[C, 0]]`, C = nfrag×N fragment indicator. Current state:

| Path | CPU | CUDA | ROCm |
|---|---|---|---|
| Full direct (exact) | Cholesky/LU on (N+nfrag) | cuSOLVER potrf/getrf | rocSOLVER potrf/getrf |
| Constraint | bordered solve | GPU-Schur (WP5-A/WP7-A) | **host CPU-Schur** |
| Batched per-fragment (drops cross-fragment Coulomb, **approximate**) | ✅ auto (nfrag≥16 & nfrag/N>20%) | ✅ WP7-B | ❌ |
| PCG iterative (**exact**, full matvec) | ✅ **block-Jacobi PC + multi-step MD warm-start** | ✅ but **plain diagonal Jacobi PC** | ❌ |

**Cost for many small fragments (e.g. N=6200, nfrag=1400):** full N×N Cholesky is O(N³) ≈ prohibitive;
batched is O(Σ N_f³) ≈ free but drops inter-molecular Coulomb (no inter-fragment charge polarization
— an approximation vs the Fortran-exact full solve); PCG+block-Jacobi is O(k·N²) with small k for
well-separated fragments and, warm-started, ≈ 1–2 iters/MD step — **exact AND fast** (the best option).

### Findings / WPs
- **ROCm — biggest gap (no advanced path).** Every nfrag>1 step does the full N×N rocSOLVER Cholesky
  (the batched/PCG/general device variants are stubbed). Two fixes: (a) route highly-fragmented EEQ to
  the CPU's `EEQSolveMethod::PCG` (block-Jacobi + warm-start, exact) or `Batched` path and upload the
  charges — reuses validated code, near-free, and **matches whatever the CPU auto-selects** (today CPU
  silently auto-batches at nfrag≥16 while ROCm does the exact full solve → the two DIVERGE for solvent;
  this also means the "ROCm==CPU" validation does not cover many-fragment systems); (b) port PCG +
  block-Jacobi to HIP (rocBLAS `dsymv` matvec + per-fragment Cholesky preconditioner). (a) first.
- **CUDA — under-preconditioned PCG.** `solveWithDeviceRHSAndGPUPCG` uses a plain diagonal Jacobi PC
  (`d_pcg_M_inv = 1/A[i,i]`), which converges slowly for many fragments. Port the CPU's **block-Jacobi**
  (per-fragment) preconditioner → far fewer iters, exact. Also the CUDA batched (WP7-B) drops
  cross-fragment Coulomb (approximate) — PCG+block-Jacobi is the exact replacement.
- **CPU — already advanced** (batched + PCG/block-Jacobi/warm-start). Remaining "more advanced": an
  FMM/treecode matvec to make PCG's cross-fragment Coulomb O(N log N) instead of O(N²) for very large
  solvent boxes (big effort, niche); and preferring the **exact** PCG+block-Jacobi over the auto-batched
  **approximation** by default when fragments are in contact (accuracy), since batched silently omits
  inter-molecular charge coupling.
- **Cannot validate here** — no many-fragment solvent test molecule in-tree (water8_cluster is nfrag=8,
  below the batched/PCG thresholds). Needs a solvent box (e.g. the referenced mixture.xyz N=6200,
  nfrag=1400) to benchmark + verify exact-vs-batched accuracy.
