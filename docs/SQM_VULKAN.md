# Native GFN1/GFN2/GFN-FF on Vulkan compute

> Status: 🤖 AI-generated, ⚙️ machine-tested, **Stage 1 (GPU eigensolver)**. NOT human
> production tested. `-gpu vulkan` runs the native GFN1/GFN2 SCF with the per-iteration
> dense symmetric eigensolve on the GPU (FP64 two-sided Jacobi); integrals, Fock build,
> density and the nuclear gradient still run on the CPU. Energies match the CPU path
> bit-for-bit on the validation runs so far (see "What was tested").

## What it is

`-gpu vulkan` selects a hand-written Vulkan compute backend for the native methods —
vendor-neutral GPU acceleration (AMD/NVIDIA/Intel) with no CUDA/ROCm dependency. Unlike
CUDA/ROCm there is **no BLAS/LAPACK on Vulkan**, so every kernel (GEMM, triangular solve,
the integral/Fock/density/gradient kernels, and the symmetric eigensolver) is a
hand-written SPIR-V compute shader.

```
./curcuma -sp mol.xyz -method gfn2 -gpu vulkan   # explicit Vulkan
./curcuma -sp mol.xyz -method gfn1 -gpu auto     # GPU if a backend is compiled
```

`-gpu vulkan` on a build without Vulkan (or with no FP64-capable device) warns and falls
back to the CPU. **FP64 (`shaderFloat64`) is required** — the SCF numerics are double
precision, and `VkContext` rejects devices without it.

## Build

```
cmake -S . -B release_vulkan -DCMAKE_BUILD_TYPE=Release \
      -DUSE_VULKAN=ON -DUSE_VULKAN_XTB=ON
cmake --build release_vulkan -j4
```

Requires the Vulkan SDK (loader + headers + `glslc`/`glslangValidator`). One GPU backend
per build dir; the default `release/` is unchanged.

## Architecture

The core `XTB` SCF loop offloads through the abstract `GpuScfBackend` seam
(`xtb_native.h`); any hook returning `false` falls back to the CPU, so the Vulkan engine is
brought up stage-by-stage and is correct at every step.

| File | Role |
|------|------|
| `qm_methods/vulkan/vk_context.{h,cpp}` | Generic Vulkan compute context: instance / FP64 device / compute queue / command pool |
| `qm_methods/vulkan/xtb_vulkan_context.{h,cpp}` | xTB device engine over a `VkContext`; mirrors `XtbGpuContext` |
| `qm_methods/xtb_vulkan_method.{h,cpp}` | `ComputationalMethod` wrapper; owns the context + the CPU `NativeXtbMethod` |
| `qm_methods/vulkan/shaders/*.comp` | Hand-written FP64 GLSL compute shaders (compiled to SPIR-V) |

Dispatch: `method_factory.cpp` `resolveNativeXtbGpuMode()` returns `"vulkan"` when
`-gpu vulkan` and `USE_VULKAN_XTB` are set, then constructs `XtbVulkanComputationalMethod`.

## Eigensolver choice

The generalized problem `F C = S C ε` reduces via the Cholesky factor `L` (two triangular
solves + back-transform) to a **standard** symmetric eigenproblem, exactly like the CUDA
path — so Vulkan needs only GEMM + TRSM + a standard symmetric eigensolver. That solver is
**cyclic one-sided Jacobi** (`shaders/jacobi_apply_f64.comp` + host pairing): the simplest
symmetric eigensolver to parallelize correctly in compute shaders, returning the full
spectrum (required for the Fermi occupation). It costs more flops than divide-and-conquer
but is chosen for correctness first.

## Stages

0. **Build + dispatch + device handshake** — done. `-gpu vulkan` builds, selects the
   backend, probes an FP64 device, falls back to CPU when absent.
1. **GPU symmetric eigensolver wired via `ExternalEigensolver`** — done. The cyclic
   two-sided Jacobi (`shaders/{angles,col,row,vec}.comp`, embedded SPIR-V) solves the
   per-iteration SCF eigenproblem on the device; the host reduces the generalized
   problem to standard form (Cholesky) and back-transforms. Validated vs CPU (below).
2-4. Resident SCF, on-device integral build, nuclear gradient (as CUDA/ROCm) — pending.
   Until then the integrals / Fock / density / gradient run on the CPU.

## What was tested

On an **AMD Radeon 890M (RADV, integrated, shaderFloat64)**, build `release_vulkan/`
(`-DUSE_VULKAN=ON -DUSE_VULKAN_XTB=ON`):

- **Eigensolver vs Eigen** (standalone, `prototype/`): random symmetric n=4..128, eigenvalues
  to ~1e-13, reconstruction / orthogonality to ~1e-14.
- **gfn2 / gfn1 single point** `-gpu vulkan` vs `-gpu none`: H2O and benzoic acid
  (C6H5COOH, 15 atoms) — energies bit-identical at 8 decimals (|dE| = 0).
- **gfn2 -opt** `-gpu vulkan`: H2O converges in 7 steps to the same minimum as the CPU
  (-5.070544 Eh), exercising repeated GPU eigensolves across SCF iterations and opt steps.
- Default non-Vulkan `release/` build stays green.

What was **NOT** tested: large systems (>~100 nao; the iGPU FP64 Jacobi is slow — this is a
correctness milestone, not a performance one), discrete GPUs, NVIDIA/Intel devices, MD,
and the on-device integral/Fock/gradient stages (still CPU). Human production testing pending.
