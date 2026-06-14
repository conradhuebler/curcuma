# Native GFN1/GFN2/GFN-FF on AMD GPUs (ROCm/HIP)

> Status: 🤖 AI-generated, ⚙️ machine-tested. NOT human production tested.
> **GFN1 = Stage 2 (device-resident SCF):** the Fock build and Mulliken populations are
> HIP `__global__` kernels (compiled by `hipcc`), the density is a rocBLAS GEMM, and the
> generalized eigensolve is rocSOLVER `dsygvd`; H0/S stay resident, only `v_ao`/`occ` up
> and `eps`/`pop`/`band` down per iteration. **GFN2 = Stage 1 (GPU eigensolver only):**
> rocSOLVER `dsygvd` per iteration, the rest on the CPU. Integrals and the nuclear gradient
> are on the CPU for both. Energies match the CPU path bit-for-bit (AMD 890M / gfx1150).

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
2. **Device-resident GFN1 SCF via `GpuScfBackend`** — done. `begin` uploads H0/S; `solve`
   builds the Fock (`k_fock` HIP kernel) and calls rocSOLVER `dsygvd` (C kept resident,
   ascending — no host sort); `density` forms P = C·diag(occ)·Cᵀ (`k_scale_cols` + rocBLAS
   `dgemm`) + Mulliken populations / band (`k_popband`). Only `v_ao`/`occ` up and
   `eps`/`pop`/`band` down per iteration. GFN2 (no multipole here) falls back to stage 1.
3-4. On-device integral build + nuclear gradient — pending (still CPU). GFN-FF ROCm pending.

## What was tested

On an **AMD Radeon 890M (gfx1150)**, build `release_rocm/` (`-DUSE_ROCM_XTB=ON`, rocSOLVER):

- **gfn1 (Stage-2 resident) / gfn2 (Stage-1) single point** `-gpu rocm` vs `-gpu none`:
  H2O and benzoic acid (15 atoms) energies bit-identical at 8 decimals (|dE| = 0).
- **gfn1 / gfn2 `-opt`**: converge to the same minima as the CPU (H2O -5.768775 / -5.070544),
  exercising the resident HIP kernels across SCF iterations and opt steps.
- Stage 0 (no rocSOLVER): device handshake + CPU fallback, energies bit-identical.
- Default non-ROCm `release/` build stays green (cli_curcumaopt_*/cli_rmsd_* 11/11).

What was **NOT** tested: large systems (iGPU FP64 is slow — correctness milestone, not
performance), discrete/CDNA GPUs, MD, on-device integrals/gradient (still CPU), GFN-FF.
