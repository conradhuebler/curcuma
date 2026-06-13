# Native GFN1/GFN2/GFN-FF on Vulkan compute

> Status: 🤖 AI-generated, ⚙️ machine-tested. NOT human production tested.
> **GFN1 = Stage 2 (device-resident SCF):** Fock build, Löwdin S⁻¹ᐟ² reduction, the
> symmetric eigensolve, density, Mulliken populations and the band energy all run on the
> GPU; H0/S stay resident and only `v_ao` (up) + eigenvalues/populations (down) cross the
> bus per iteration. **GFN2 = Stage 1 (GPU eigensolver only):** the per-iteration
> eigensolve is on the GPU, the rest on the CPU. Integrals and the nuclear gradient are
> on the CPU for both. Energies match the CPU path bit-for-bit on the validation runs.

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
1. **GPU symmetric eigensolver via `ExternalEigensolver`** — done (used by GFN2). The
   cyclic two-sided Jacobi (`shaders/{angles,col,row,vec}.comp`) solves the per-iteration
   eigenproblem; the host reduces (Cholesky) and back-transforms.
2. **Device-resident GFN1 SCF via `GpuScfBackend`** — done. `begin` uploads H0/S and
   builds the resident X = S⁻¹ᐟ² (Jacobi(S) + `scale_cols` + `gemm`, no triangular
   solve); `solve` does the Fock build (`fock`), Ã = X·F·X (`gemm`), Jacobi, C = X·C̃
   (`gemm`); `density` forms P = C·diag(occ)·Cᵀ + Mulliken populations + band
   (`scale_cols`/`gemm`/`popband`). Only `v_ao`/`occ` up and `eps`/`pop`/`band` down
   per iteration. GFN2 (no multipole here) falls back to stage 1.
3. On-device integral build (CN/S/H0/γ + GFN2 multipole) — pending (built on CPU,
   uploaded once per geometry in `begin`).
4. On-device nuclear gradient — pending (CPU); GFN2 multipole resident SCF (Stage 2b)
   — pending.

## What was tested

On an **AMD Radeon 890M (RADV, integrated, shaderFloat64)**, build `release_vulkan/`
(`-DUSE_VULKAN=ON -DUSE_VULKAN_XTB=ON`):

- **Eigensolver vs Eigen** (standalone, `prototype/`): random symmetric n=4..128, eigenvalues
  to ~1e-13, reconstruction / orthogonality to ~1e-14. Generalized solve (Löwdin S⁻¹ᐟ²
  chain) vs Eigen GeneralizedSelfAdjointEigenSolver: residual ~1e-14.
- **gfn1 single point + opt** (Stage-2 device-resident path) `-gpu vulkan` vs `-gpu none`:
  H2O and benzoic acid (C6H5COOH, 15 atoms) bit-identical at 8 decimals (|dE| = 0); H2O
  opt converges to the same minimum (-5.768775 Eh).
- **gfn2 single point + opt** (Stage-1 eigensolver path): H2O / benzoic acid bit-identical;
  opt to -5.070544 Eh.
- Default non-Vulkan `release/` build stays green (cli_curcumaopt_*/cli_rmsd_* 11/11).

What was **NOT** tested: large systems (>~100 nao; the iGPU FP64 Jacobi is slow — this is a
correctness milestone, not a performance one), discrete GPUs, NVIDIA/Intel devices, MD,
and the on-device integral/Fock/gradient stages (still CPU). Human production testing pending.
