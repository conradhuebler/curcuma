# Native GFN1/GFN2/GFN-FF on AMD GPUs (ROCm/HIP)

> Status: 🤖 AI-generated, ⚙️ machine-tested. NOT human production tested.
> **GFN1 = Stage 4 (fully device-resident):** the integral build (CN, overlap S, bare
> Hamiltonian H0, Cholesky L, Coulomb γ), the SCF (Fock/density/populations/eigensolve) and
> the nuclear gradient (repulsion + on-site CN + H0/Pulay + Coulomb) all run on the GPU via
> HIP `__global__` kernels + rocBLAS + rocSOLVER; only the dispersion gradient + CN
> chain-rule stay on the host. `-opt`/`-md` need no host integral or gradient build.
> **GFN2 = device-resident multipole SCF (Stage 2b/R-AP2):** the Fock build (incl. the
> anisotropic dp/qp term via `k_add_fock_multipole`), the rocSOLVER eigensolve, density and
> the atomic multipole moments (`k_multipole_moments`) run on the device; only `v_dp`/`v_qp`
> up + eps/moments down per iteration. The potential stays on the host. **The nuclear gradient
> incl. the multipole-integral Pulay term now runs on the device too (Stage 4 / R-AP3); only
> the multipole SD/DD/SQ interaction gradient + dispersion + CN chain-rule stay on the host.**
> Energies + gradient + the `-opt` trajectory match the CPU (AMD 890M / gfx1150).
> **FP32 mixed precision (X-AP3) is ON by default for `-gpu rocm`** and is a real win:
> far-from-convergence iterations use `rocsolver_ssygvd` (FP32), reverting to FP64 near
> convergence. On `complex` (231 atoms) the resident path is already ~6-8× the CPU
> (rocSOLVER is mature), and FP32 adds **GFN1 1.44× (1844→1281 ms), GFN2 1.27× (1924→1510
> ms)** on top. Energies bit-identical; the gradient at the loose default `scf_threshold`
> (1e-5) differs ~1e-7 from CPU (use `-scf_threshold 1e-8` for bit-identical gradients, or
> `-scf_mixed_precision false` to disable). See SQM_GPU_ROADMAP.md X-AP3.

## What it is

`-gpu rocm` selects the ROCm/HIP backend for the native methods, mirroring the existing
CUDA path (`-gpu cuda`). It is the AMD sibling of [docs/SQM_GPU.md](SQM_GPU.md): the
per-iteration eigensolve uses **rocSOLVER** (Stage 1, live); later stages port the CUDA
`.cu` integral/SCF/gradient kernels into a separate `rocm/` source tree with `hipcc`
(separate HIP files, not shared-via-hipify, so the CUDA and ROCm builds evolve
independently).

```
./curcuma -sp mol.xyz -method gfn2 -gpu rocm    # explicit ROCm
./curcuma -sp mol.xyz -method gfn1 -gpu auto    # GPU if a backend is compiled
./curcuma -sp mol.xyz -method gfn2              # CPU (default)
```

`-gpu rocm` on a build without ROCm warns and falls back to the CPU.

## Dependencies

ROCm is **not** pulled in by the default build; it is only needed for `-DUSE_ROCM=ON`.
All ROCm pieces live under `/opt/rocm`, so pass `-DCMAKE_PREFIX_PATH=/opt/rocm`.

| Component | Purpose | Arch/Manjaro package | CMake / flag |
|-----------|---------|----------------------|--------------|
| HIP runtime + `hipcc` + headers | compile/run HIP, `libamdhip64` | `hip-runtime-amd` | `find_package(hip)` → `hip::host`, `hip::device` |
| HIP/Clang compiler (`amdclang++`) | the HIP language compiler | `rocm-llvm` | auto-detected by `enable_language(HIP)` |
| ROCm core + device bitcode | gfx device libs, base | `rocm-core`, `rocm-device-libs` | — |
| Device query | resolve your `gfx` arch | `rocminfo` | `rocminfo \| grep gfx` |
| **rocBLAS** | GEMM / TRSM (Stage 1+) | `rocblas` | `find_package(rocblas)` → `HAVE_ROCBLAS` |
| **rocSOLVER** | symmetric eigensolver (Stage 1) | `rocsolver` | `find_package(rocsolver)` → `HAVE_ROCSOLVER` |
| hipBLAS (optional) | portable BLAS wrapper over rocBLAS | `hipblas` | not required (rocBLAS used directly) |
| hipify (optional) | seed HIP from the CUDA `.cu` kernels | `hipify-clang` | dev-time only |
| GPU + kernel driver | a ROCm-capable AMD GPU | `amdgpu` (kernel) | set `CMAKE_HIP_ARCHITECTURES` to your `gfx` |

- **Stage 0** (current — device handshake + CPU fallback) calls only the *host-side* HIP
  runtime API (device query + stream), so it is compiled as **plain C++** against
  `<hip/hip_runtime_api.h>` and linked against **`libamdhip64`** — no HIP-language
  compilation, no `amdclang`, no `--offload-arch`. It needs only `hip-runtime-amd`
  (the runtime + headers). rocBLAS/rocSOLVER are detected (`find_package(... QUIET)`,
  `HAVE_ROCBLAS`/`HAVE_ROCSOLVER`) but not yet linked.
- **Stage 1** (GPU eigensolver) **requires `rocsolver`** (and `rocblas`), and introduces
  real HIP-language device kernels. Those build in a dedicated unit with `hipcc`/`amdclang`;
  the HIP-only `--offload-arch`/`--hip-link` flags must be kept off the GNU/`g++` link of
  `curcuma_core` and the executables (otherwise the link flips to `ld.lld` and drops the
  GNU OpenMP runtime, `libgomp`).

## Build

```
cmake -S . -B release_rocm -DCMAKE_BUILD_TYPE=Release \
      -DUSE_ROCM=ON \
      -DCMAKE_PREFIX_PATH=/opt/rocm \
      -DCMAKE_HIP_ARCHITECTURES=gfx1150     # match your GPU: rocminfo | grep gfx
cmake --build release_rocm -j4
```

`gfx1150` = AMD Radeon 880M/890M (RDNA 3.5); use e.g. `gfx1100` (RX 7900), `gfx1030`
(RX 6800/6900), `gfx90a` (MI200). One GPU backend per build dir (like `release_cuda/`);
the default `release/` is unchanged.

## Architecture

The core `XTB` SCF loop is backend-neutral — it offloads through the abstract
`GpuScfBackend` seam (`xtb_native.h`), where any hook returning `false` falls back to the
CPU for the whole calculation. The ROCm engine is therefore brought up the same staged way
as CUDA, correct at every step:

| File | Role |
|------|------|
| `qm_methods/rocm/xtb_hip_context.{h,cpp}` | HIP context: device handshake + `solveGeneralized` (rocSOLVER `dsygvd`). Host-callable — device memory via `hipMalloc`/`hipMemcpy`, no device kernels yet, so it compiles as plain C++ |
| `qm_methods/xtb_hip_method.{h,cpp}` | `ComputationalMethod` wrapper; owns the context + the CPU `NativeXtbMethod`; installs the `ExternalEigensolver` hook |
| `ff_methods/rocm/` | (later stage) hipified GFN-FF kernels + workspace |

Dispatch: `method_factory.cpp` `resolveNativeXtbGpuMode()` returns `"rocm"` when
`-gpu rocm` and `USE_ROCM` are set, then constructs `XtbHipComputationalMethod`.

Unlike the Vulkan path (which hand-writes a Jacobi eigensolver + Löwdin reduction),
rocSOLVER provides the **generalized** symmetric-definite solver directly, so the host
reduction/back-transform is unnecessary and the same hook serves GFN1 and GFN2.

The HIP context (`xtb_hip_context.hip`) is compiled by **`hipcc` into a plain relocatable
object** (CMake `add_custom_command`) and linked into `curcuma_core` with `g++` as an
`EXTERNAL_OBJECT`. This builds the device kernels while keeping the HIP-only
`--offload-arch`/`--hip-link` flags off the GNU link (which would otherwise flip the link
to `ld.lld` and drop GNU OpenMP / `libgomp`). No `enable_language(HIP)` is used.

## Stages

0. **Build + dispatch + device handshake** — done.
1. **GPU eigensolver via rocSOLVER `dsygvd`** — done (used by GFN2). The per-iteration
   `F C = S C ε` is solved on the device (`ExternalEigensolver` hook: host builds
   `S = L·Lᵀ`, rocSOLVER returns `C`/`ε`).
2. **Device-resident GFN1 SCF via `GpuScfBackend`** — done. `solve` builds the Fock
   (`k_fock`) and calls rocSOLVER `dsygvd` (C kept resident, ascending — no host sort);
   `density` forms P = C·diag(occ)·Cᵀ (`k_scale_cols` + rocBLAS `dgemm`) + Mulliken
   populations / band (`k_popband`). Only `v_ao`/`occ`/`eps`/`pop`/`band` cross per iteration.
3. **On-device integral build via `beginBasis`/`beginComputed`** — done. `beginBasis` uploads
   the flattened basis + element tables once; `beginComputed` runs `k_cn` → `k_self_energy`
   → `k_overlap_h0` (s/p contracted Gaussian overlap + the GFN1/GFN2 H0 scaling) →
   `rocsolver_dpotrf` (Cholesky) → `k_gamma`, writing S/H0 into the resident buffers, so the
   host skips its integral build (`downloadOverlap`/`H0`/`Cholesky`/`Gamma` for its
   bookkeeping). The device math is a verbatim port of the proven CUDA kernels
   (`rocm/xtb_hip_integrals.hiph`). Used by both GFN1 and GFN2.
3m. **On-device GFN2 multipole integrals (R-AP1)** — done. `beginMultipoleComputed` runs
   `k_multipole_ints` (one thread per AO pair: global-origin `d_cgto_multipole` then the
   per-column origin shift + traceless transform with the resident overlap `dS0`);
   `downloadMultipoleInts` fetches dp_int(3·nao²)/qp_int(6·nao²) so the host GFN2 SCF skips
   its O(nao²) `setupMultipole` integral loop. `d_moment1d`/`d_type_to_cart`/
   `d_primitive_multipole`/`d_cgto_multipole` are verbatim ports of the CUDA `.cuh` into
   `xtb_hip_integrals.hiph`. GFN2 runs the host SCF (rocSOLVER eigensolver), so the device
   integrals feed the host Fock — bit-identical GFN2 energy proves them. (The GFN2 resident
   multipole SCF later landed as Stage 2b/R-AP2 and the multipole-integral gradient as R-AP3.)
4. **On-device nuclear gradient via `GpuScfBackend::gradient`** — done (GFN1). Kernels
   `k_grad_repulsion` (section 1), `k_grad_cn_onsite` (2a, dEdcn diagonal), `k_grad_h0_pulay`
   (2b, the overlap-derivative Pulay term using the energy-weighted density W = C·diag(2ε)·Cᵀ,
   built via `k_scale_cols` + rocBLAS), `k_grad_coulomb` (3, γ-derivative). The overlap
   derivative `d_cgto_overlap_grad` (Obara-Saika, s/p) is in `xtb_hip_integrals.hiph`. The
   dispersion gradient + CN chain-rule stay on the host. **GFN2 multipole-integral Pulay
   gradient done (R-AP3):** `k_grad_h0_pulay` adds the dp/qp-integral derivative term
   (`d_cgto_multipole_grad_transformed`) contracted with the converged `v_dp`/`v_qp`, so the
   GFN2 nuclear gradient is device-resident too; only the multipole SD/DD/SQ interaction
   gradient (§5) stays on the host. GFN-FF ROCm still pending.

The remaining work (GFN-FF, device solvation) is broken into work packages in
[SQM_GPU_ROADMAP.md](SQM_GPU_ROADMAP.md) (ROCm = `R-AP*`).

## What was tested

On an **AMD Radeon 890M (gfx1150)**, build `release_rocm/` (`-DUSE_ROCM=ON`, rocSOLVER):

- **gfn1 / gfn2 single point** `-gpu rocm` vs `-gpu none` (CN/S/H0/L/γ built on the GPU,
  SCF log "CN/S/H0/L built on GPU; no nao^2 upload"): H2O and benzoic acid (15 atoms)
  energies bit-identical at 8 decimals (|dE| = 0) — confirms the device overlap/H0/γ are
  correct (any integral error would shift the energy).
- **gfn1 `-opt` with the device gradient**: the *entire* trajectory matches the CPU
  step-by-step — identical energies AND gradient norms at every step (H2O 8 steps to
  -5.768775; benzoic acid to -27.417443). A wrong gradient would diverge immediately, so
  this confirms the device repulsion/Pulay/Coulomb gradient is correct.
- **gfn2 device multipole integrals (R-AP1)**: with `k_multipole_ints` building dp_int/qp_int
  on the device (SCF log "GFN2 multipole integrals built on GPU device"), gfn2 `-sp` energy
  bit-identical to CPU (8 dp) over the full 12-molecule `sqm_reference` set incl. the
  231-atom `complex`; gfn2 `-opt` (NH3) identical trajectory (integrals rebuilt each geometry,
  feeding the host Fock AND the host multipole gradient).
- **gfn2 device-resident multipole SCF (R-AP2)**: `k_add_fock_multipole` + `k_multipole_moments`
  (log "device-resident GFN2 path"), gfn2 `-sp` energy bit-identical + gradient norm matches
  CPU, gfn2 `-opt` (NH3) step-by-step identical. NOT faster than CPU on the iGPU (FP64
  eigensolve-bound).
- **gfn2 device nuclear gradient incl. multipole-integral Pulay (R-AP3)**: `k_grad_h0_pulay`
  with the dp/qp-integral derivative term (`d_cgto_multipole_grad_transformed` ← verbatim CUDA
  port), `v_dp`/`v_qp` uploaded once per gradient. With `-scf_mixed_precision false` (FP64),
  gfn2 `-sp` gradient norm bit-identical to CPU over the `sqm_reference` set incl. the 231-atom
  `complex` (reldiff = 0). The ROCm-routed finite-difference gradient check
  (`test_xtb_gradient`) is byte-identical to the CPU path (same MaxErr/RMSErr per molecule),
  proving the device analytic gradient equals the FD-validated CPU AP5 gradient. gfn2 `-opt`
  (triose) tracks the CPU step-by-step then diverges in the LBFGS tail exactly like the
  known-good GFN1 device gradient (rocSOLVER vs Eigen eigenvectors ~1e-13, amplified by LBFGS
  — not an R-AP3 error). At the default mixed-precision path the gradient differs ~1e-7 (the
  documented FP32 lever; use `-scf_threshold 1e-8` or `-scf_mixed_precision false`).
- Stage 0 (no rocSOLVER): device handshake + CPU fallback, energies bit-identical.
- Default non-ROCm `release/` build stays green (cli_curcumaopt_*/cli_rmsd_* 11/11).

d shells (main-group S/P/Cl/Si/…) are now device-resident on ROCm (X-I1 B6 port,
2026-06-28): the HIP integral/SCF/gradient kernels carry the cartesian->spherical dtrafo, so
`-gpu rocm` energy is bit-identical to CPU at 8 dp on H2S/PH3/SiH4/HCl (gfn1+gfn2) and `-opt`
tracks the CPU trajectory. See [SQM_DSHELL_WP.md](SQM_DSHELL_WP.md) B6.

What was **NOT** tested (native xTB): large systems (iGPU FP64 is slow — correctness
milestone, not performance), discrete/CDNA GPUs, MD; transition-metal d (enabled but
unvalidated — only main-group d checked, same as the CPU/CUDA paths).

## GFN-FF on ROCm (`-gpu rocm`, `USE_ROCM`, June 2026) — 🤖 machine-tested

Separate from the native xTB ROCm path above; mirrors the CUDA GFN-FF GPU pipeline.

- **Build**: `cmake -DUSE_ROCM=ON -DUSE_D4=OFF -DUSE_D3=OFF
  -DROCM_GPU_ARCH=gfx1150 -DCMAKE_PREFIX_PATH=/opt/rocm` (a dedicated `release_rocm/` dir, not
  the canonical CPU `release/`). Needs `rocsolver` + `rocblas`.
- **Device code** (`ff_methods/rocm/`): `gfnff_rocm.hip` is a **single** hipcc TU (the
  single-TU rule avoids cross-TU device linking so the final g++ link stays libgomp-safe, same
  reason as `xtb_hip_context.hip`). It bundles the hipified kernels (`gfnff_kernels_hip.h`),
  SoA helpers (`gfnff_soa_hip.h`), the workspace impl (port of `cuda/ff_workspace_gpu.cu` —
  3 streams + hipGraph), gpu-utils, and the EEQ solver (`eeq_solver_hip.hiph`). Distinct class
  names `FFWorkspaceHip`/`EEQSolverHip` ⇒ **zero change to the CUDA path**, both backends can
  coexist. Host wrapper `qm_methods/gfnff_hip_method.cpp` (g++) is a rename-copy of the CUDA
  wrapper; factory routes `-gpu rocm` → `GFNFFHipComputationalMethod`.
- **Warp reduction**: HIP 7.x requires a 64-bit `__shfl_down_sync` mask (AMD wavefronts can be
  64 lanes). gfx11xx/RDNA HIP compute is wave32, so the 32-lane indexing stays correct;
  `warpReduceSum` now passes `0xFFFFFFFFFFFFFFFFull`.
- **EEQ**: `solveWithDeviceRHS`/`solve` build the N×N Coulomb matrix on device and factor it
  with rocSOLVER `dpotrf` (LU `dgetrf` fallback for indefinite matrices), solving `[b | Cᵀ]`
  with `dpotrs`/`dgetrs` → z1/Z2. The host wrapper's **exact CPU Schur complement** applies the
  charge constraint (works for nfrag=1 and nfrag>1). The device-resident GPU-Schur/PCG/batched
  variants are intentionally not ported (they return false → CPU-Schur path).
- **Validated** (Radeon 890M/gfx1150): single-point **energy and gradient norm** match CPU
  ≤1e-7 Eh / Eh·Bohr⁻¹ on water, CH4, caffeine, 231-atom `complex` (FP reduction-order only);
  `-opt` trajectories track CPU (small FP accumulation over steps). `ctest -R
  cli_gfnff_gpu_02_rocm_singlepoint` (caffeine, label `rocm`).
- **Performance — per-atom gather (the real MD win)**: profiling the phase-2 kernels (polymer/1410,
  gradient mode) showed the pairwise **Coulomb** (~993k pairs, no cutoff) and the deferred
  **dispersion gradient** (~943k pairs, grad + dEdcn atomics) dominated — ~6M+ FP64 `atomicAdd` on
  RDNA. Both were rewritten as **per-atom GATHER** kernels (`k_coulomb_gather` / `k_dispersion_gather`,
  built from a per-atom CSR adjacency; the dynamic dc6dcn stays pair-indexed and is read via a
  CSR pair-index + is_i flag → no per-step rebuild; Coulomb additionally has a dense, coalesced,
  shared-mem-tiled `k_coulomb_dense` that recomputes gamma from the per-atom EEQ alpha). Each atom
  accumulates locally and writes grad with 3 (+1 for dEdcn) atomicAdd, bit-identical to the pair
  kernels. Result on polymer/1410: phase-2 (Coulomb+Disp+DMA) **122 → 15 ms**, full energy+gradient
  eval **~150 → ~26 ms**, and the **MD went from slower-than-CPU to ~2.3× faster** (30 fs:
  10.1 s CPU-powersave vs 4.4 s ROCm, back-to-back). Note: the gfn2 ROCm path is efficient because
  its O(N²)/O(N³) work is rocSOLVER/rocBLAS on small systems; GFN-FF's O(N²) FP64-erf Coulomb has no
  BLAS analog, so the gather (not BLAS) is the lever. FP64 atomics themselves are fine on RDNA
  (gfn2 uses them) — the issue was *contention* at ~1M pairs, which the gather removes.
- **GPU-wedge hardening**: `FFWorkspaceHip` installs a SIGINT/SIGTERM/SIGABRT handler that
  best-effort `hipDeviceReset()`s, restores the default disposition and re-raises. Killing a process
  with a HIP kernel in flight *without* device teardown wedges the RDNA compute queue (every later
  submission then hangs at init until a GPU reset/reboot); the handler makes Ctrl-C / `kill` /
  `timeout` clean up instead. Only `kill -9` (SIGKILL) can't be caught — avoid it on in-flight GPU
  runs, and do not run multiple GPU jobs concurrently on one iGPU. Validated: a polymer MD killed
  with SIGTERM mid-run no longer wedges the GPU.
- **D4 is real, not a fallback** (corrected): the build is `USE_D4=OFF`, but GFN-FF's dispersion
  does NOT depend on that flag. `generateGFNFFDispersionPairs()` (gfnff_method.cpp:9726) always
  runs the **self-contained** `D4ParameterGenerator` (Casimir-Polder C6, `dispersion/d4param_generator`,
  compiled unconditionally), which has no `#ifdef USE_D4` gate and never calls the external
  `dftd4interface`/`curcuma_d4`/LAPACKE lib (that lib backs only the standalone `-d4` method).
  Verified at runtime: the "D4 C6 Reference Matrix Pre-computation … Casimir-Polder" path executes
  on the `USE_D4=OFF` build. So the ROCm GFN-FF runs **true-D4 GFN-FF** and matches the CPU (same
  generator). The earlier "free-atom fallback / not-true-D4" note was wrong — `generateFreeAtomDispersion()`
  is only reached if D4 construction throws.
- **Opt-in GPU flags mirrored from CUDA (Jun 2026, default OFF)**:
  - `-gfnff.eeq_mixed_precision` (**WP-B**): FP32 `rocsolver_spotrf` factor + FP64 iterative
    refinement (FP64 residual via `rocblas_dsymm`, FP32 `spotrs` correction, dsposv pattern) for
    full FP64 accuracy at a fraction of the FP64-factor cost on FP64-weak RDNA. `mixedFactorHip`/
    `mixedSolveRefineHip` + `k_cast_d2f_hip`/`k_cast_f2d_hip`/`k_axpy_f2d_hip` in
    `eeq_solver_hip.hiph`; gated to nfrag<=1 in `eeqBuildFactorSolve`, auto-falls back to FP64
    dpotrf/LU when the FP32 factor is not SPD. `eeq_mixed_precision_iters` (default 2).
  - `-gfnff.gpu_disp_pairs_on_device` (**WP-A**): two-pass on-device D4 dispersion pair build
    (`k_disp_pairs_count`/`k_disp_pairs_build` + `generateDispersionPairListOnGPU`) replacing the
    host O(N²) loop + upload; ROCm additionally **rebuilds the per-atom CSR adjacency** for
    `k_dispersion_gather` from the device list (no CUDA precedent). Both bit-identical to the host/
    FP64 ROCm path on caffeine + 231-atom `complex` (energy + opt). Neither is a measured SP
    speedup (residency/correctness milestones — same as CUDA).
- **Many-fragment EEQ → exact CPU PCG** (`-gfnff.eeq_rocm_cpu_fragment_threshold`, default 16):
  for `nfrag >= threshold` (solvent boxes) the dispatch routes to the exact CPU PCG+block-Jacobi
  +warm-start solver (`finalizeCNForCPU` then `prepareCNAndEEQ`) instead of the dense N×N device
  Cholesky (O(N³), intractable; the device PCG/batched/general-Schur variants are stubbed). This
  also fixes the device-vs-CPU charge divergence for multi-fragment systems: water8_cluster
  (threshold 4) opt now tracks the pure-CPU `-gpu none` reference (the device path was the
  divergent one). Set 0 to force the device solve.
- **NOT done / honest**: untested — MD long-run stability, very large solvent boxes (no such
  molecule in the test tree; only correctness CPU==ROCm is CI-checkable), discrete/CDNA GPUs (the
  warp reduction is wave32-only — see roadmap), and the **device-resident HIP EEQ PCG/block-Jacobi**
  (CUDA has WP7-D; ROCm uses the exact CPU PCG routing above — the device-resident port is the open
  follow-up).
