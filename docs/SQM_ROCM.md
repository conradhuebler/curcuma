# Native GFN1/GFN2/GFN-FF on AMD GPUs (ROCm/HIP)

> Status: 🤖 AI-generated, **Stage 0 (build + dispatch scaffold)**. NOT yet functional
> on-device, NOT human production tested. The ROCm path currently performs the device
> handshake and then runs the validated CPU pipeline; the HIP kernels (hipified from the
> CUDA stack) land in later stages.

## What it is

`-gpu rocm` selects the ROCm/HIP backend for the native methods, mirroring the existing
CUDA path (`-gpu cuda`). It is the AMD sibling of [docs/SQM_GPU.md](SQM_GPU.md):
cuSOLVER/cuBLAS map to rocSOLVER/hipBLAS, and the CUDA `.cu` kernels are hipified into a
separate `rocm/` source tree (the operator chose separate HIP files over shared-via-hipify
so the CUDA and ROCm builds evolve independently).

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
| `qm_methods/rocm/xtb_hip_context.{h,hip}` | HIP device engine (hipBLAS/rocSOLVER + stream); mirrors `XtbGpuContext` |
| `qm_methods/xtb_hip_method.{h,cpp}` | `ComputationalMethod` wrapper; owns the context + the CPU `NativeXtbMethod` |
| `ff_methods/rocm/` | (later stage) hipified GFN-FF kernels + workspace |

Dispatch: `method_factory.cpp` `resolveNativeXtbGpuMode()` returns `"rocm"` when
`-gpu rocm` and `USE_ROCM_XTB` are set, then constructs `XtbHipComputationalMethod`.

## Stages (mirroring the CUDA port)

0. **Build + dispatch + device handshake** — done (this commit). `-gpu rocm` builds,
   selects the backend, probes the device, runs CPU.
1. GPU eigensolver (`rocsolver_dsygst` + `rocsolver_dsyevd`) via `ExternalEigensolver`.
2. Device-resident SCF (GFN1 isotropic; GFN2 multipole).
3. On-device integral build (CN/S/H0/L/γ/multipole).
4. On-device nuclear gradient (`-opt`/`-md` device-resident).
   GFN-FF HIP kernels are ported alongside stages 2-4.

## Validation plan

Reuse the CUDA harness relabelled `rocm_*`: elementwise match vs the CPU path
(integrals ~1e-15, gradient ~1e-15 Eh/A) and the 12-molecule SQM set ≤1e-8 Eh vs tblite.

## What was tested

- Stage 0 only: factory dispatch (`-gpu rocm|auto`, fallback warnings) verified on the CPU
  build; the `xtb_hip_method.cpp` host wrapper compiles against the project flags. The
  `.hip` device TU is compile-pending a full `hipcc` toolchain.
- NOT tested: any on-device numerics (no functional ROCm kernels yet).
