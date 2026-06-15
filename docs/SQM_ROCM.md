# Native GFN1/GFN2/GFN-FF on AMD GPUs (ROCm/HIP)

> Status: 🤖 AI-generated, ⚙️ machine-tested. NOT human production tested.
> **GFN1 = Stage 4 (fully device-resident):** the integral build (CN, overlap S, bare
> Hamiltonian H0, Cholesky L, Coulomb γ), the SCF (Fock/density/populations/eigensolve) and
> the nuclear gradient (repulsion + on-site CN + H0/Pulay + Coulomb) all run on the GPU via
> HIP `__global__` kernels + rocBLAS + rocSOLVER; only the dispersion gradient + CN
> chain-rule stay on the host. `-opt`/`-md` need no host integral or gradient build.
> **GFN2:** uses the same device integral build (S/H0/γ **and now the dp/qp AO multipole
> integrals, Stage 3m/R-AP1**), then the Stage-1 rocSOLVER eigensolver per iteration; the
> multipole Fock + gradient stay on the host. Energies (and the full `-opt` trajectory)
> match the CPU bit-for-bit (AMD 890M / gfx1150).

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
      -DUSE_ROCM=ON -DUSE_ROCM_XTB=ON \
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
`-gpu rocm` and `USE_ROCM_XTB` are set, then constructs `XtbHipComputationalMethod`.

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
   integrals feed the host Fock — bit-identical GFN2 energy proves them. The GFN2 resident
   multipole SCF (Stage 2b) + multipole gradient stay pending.
4. **On-device nuclear gradient via `GpuScfBackend::gradient`** — done (GFN1). Kernels
   `k_grad_repulsion` (section 1), `k_grad_cn_onsite` (2a, dEdcn diagonal), `k_grad_h0_pulay`
   (2b, the overlap-derivative Pulay term using the energy-weighted density W = C·diag(2ε)·Cᵀ,
   built via `k_scale_cols` + rocBLAS), `k_grad_coulomb` (3, γ-derivative). The overlap
   derivative `d_cgto_overlap_grad` (Obara-Saika, s/p) is in `xtb_hip_integrals.hiph`. The
   dispersion gradient + CN chain-rule stay on the host. GFN2 multipole gradient + GFN-FF
   ROCm still pending.

The remaining work (GFN2 resident multipole SCF + multipole gradient, GFN-FF, device solvation) is
broken into work packages in [SQM_GPU_ROADMAP.md](SQM_GPU_ROADMAP.md) (ROCm = `R-AP*`).

## What was tested

On an **AMD Radeon 890M (gfx1150)**, build `release_rocm/` (`-DUSE_ROCM_XTB=ON`, rocSOLVER):

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
- **gfn2 `-opt`**: device integrals + Stage-1 eigensolver; bit-identical to the CPU.
- Stage 0 (no rocSOLVER): device handshake + CPU fallback, energies bit-identical.
- Default non-ROCm `release/` build stays green (cli_curcumaopt_*/cli_rmsd_* 11/11).

What was **NOT** tested: large systems (iGPU FP64 is slow — correctness milestone, not
performance), discrete/CDNA GPUs, MD, the GFN2 multipole gradient (GFN2 gradient on host),
GFN-FF; the overlap/gradient cover s/p only (H/C/N/O...), no d shells (as the CPU native path).
