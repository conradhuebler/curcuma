# Native GFN1/GFN2/GFN-FF on Vulkan compute

> Status: 🤖 AI-generated, ⚙️ machine-tested. NOT human production tested.
> **Stage 3 (on-device integral build):** the CN, overlap S, bare Hamiltonian H0, the
> Löwdin X = S⁻¹ᐟ² and the Coulomb γ are built by SPIR-V kernels on the device and consumed
> in place (no per-geometry nao² upload). **GFN1 = device-resident SCF (Stage 2):** Fock
> build, the symmetric eigensolve, density, Mulliken populations and the band energy all run
> on the GPU; only `v_ao` (up) + eigenvalues/populations (down) cross the bus per iteration.
> **GFN2 = device integrals + GPU eigensolver (Stage 1):** the per-iteration eigensolve is
> on the GPU; the multipole integrals and the rest of the potential are on the CPU. The
> nuclear gradient is on the CPU for both (Stage 4 pending). Energies match the CPU path
> bit-for-bit on the validation runs.

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

## Dependencies

Vulkan is **not** pulled in by the default build; it is only needed for `-DUSE_VULKAN=ON`.
It is much lighter than ROCm — just the loader, headers and an FP64-capable driver. No
vendor SDK is required.

| Component | Purpose | Arch/Manjaro package | CMake / flag |
|-----------|---------|----------------------|--------------|
| Vulkan loader | `libvulkan.so` (runtime) | `vulkan-icd-loader` | `find_package(Vulkan)` → `Vulkan::Vulkan` |
| Vulkan headers | `vulkan/vulkan.h` (build) | `vulkan-headers` | `find_package(Vulkan)` |
| Driver / ICD (FP64) | the actual GPU backend | AMD `vulkan-radeon` (RADV, **no ROCm needed**) or `amdvlk`; NVIDIA `nvidia-utils`; Intel `vulkan-intel` | — |
| `glslc` / `glslangValidator` | recompile shaders → SPIR-V | `shaderc` / `glslang` | **dev-time only** — the `.spv.inc` are committed |
| `vulkaninfo` (optional) | check the device + `shaderFloat64` | `vulkan-tools` | — |

- The build itself needs only `vulkan-icd-loader` + `vulkan-headers` (the SPIR-V is
  committed as `shaders/*.spv.inc`, embedded by `spirv_kernels.h`). `glslc`/`glslangValidator`
  are needed only to regenerate them (`shaders/compile_shaders.sh`).
- **Runtime requirement: a device with `shaderFloat64`** (FP64 in compute shaders). Check
  with `vulkaninfo | grep shaderFloat64`. `VkContext` rejects devices without it (→ CPU).
- AMD users get Vulkan via the open Mesa **RADV** driver — `-gpu vulkan` therefore works
  on AMD **without** installing the ROCm stack.

## Build

```
cmake -S . -B release_vulkan -DCMAKE_BUILD_TYPE=Release \
      -DUSE_VULKAN=ON -DUSE_VULKAN_XTB=ON
cmake --build release_vulkan -j4
```

One GPU backend per build dir; the default `release/` is unchanged.

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
3. **On-device integral build (CN/S/H0/L/γ)** — done. `beginBasis` uploads the
   molecule-constant flattened basis + element tables once; `beginComputed` (per geometry)
   runs `cn` → `self_energy` → `overlap_h0` → `gamma` SPIR-V kernels, then builds the Löwdin
   X from the device S. `xtb_native.cpp` downloads S/H0/γ (and derives L = chol(S) host-side
   via Eigen LLT — no device triangular solve) in place of the host build. Active for both
   GFN1 and GFN2 (the s/p subset); GFN2 multipole integrals (dp/qp) are still CPU.
4. On-device nuclear gradient — pending (CPU); GFN2 multipole resident SCF (Stage 2b)
   — pending.

The remaining stages (gradient, GFN2 multipole stack, GFN-FF, device solvation) are broken
into work packages in [SQM_GPU_ROADMAP.md](SQM_GPU_ROADMAP.md) (Vulkan = `V-AP*`).

## What was tested

On an **AMD Radeon 890M (RADV, integrated, shaderFloat64)**, build `release_vulkan/`
(`-DUSE_VULKAN=ON -DUSE_VULKAN_XTB=ON`):

- **Eigensolver vs Eigen** (standalone, `prototype/`): random symmetric n=4..128, eigenvalues
  to ~1e-13, reconstruction / orthogonality to ~1e-14. Generalized solve (Löwdin S⁻¹ᐟ²
  chain) vs Eigen GeneralizedSelfAdjointEigenSolver: residual ~1e-14.
- **Stage-3 device integral build active** (confirmed at `-verbosity 2`): both methods log
  "SCF: integrals built on GPU device (S/H0/Lowdin/gamma downloaded; host build skipped)";
  GFN1 additionally logs "device-resident GFN1 path (CN/S/H0/L built on GPU; no nao^2 upload)".
- **gfn1 + gfn2 single point** `-gpu vulkan` vs `-gpu none` over the full 12-molecule
  `test_cases/sqm_reference` set (H2, He2, LiH, H2O, NH3, CH4, HCN, C6H6, triose, caffeine,
  acetic_acid_dimer, and the 231-atom `complex`): all 24 bit-identical at 8 decimals (|dE| = 0).
- **gfn1 + gfn2 opt** (H2O): converges to the CPU minimum (gfn1 -5.768775, gfn2 -5.070544 Eh).
- Default non-Vulkan `release/` build stays green (cli_curcumaopt_*/cli_rmsd_* 11/11).

What was **NOT** tested: large systems (>~100 nao; the iGPU FP64 Jacobi is slow — this is a
correctness milestone, not a performance one), discrete GPUs, NVIDIA/Intel devices, MD, the
GFN2 multipole integrals + multipole resident SCF, and the on-device nuclear gradient (still
CPU, Stage 4). Human production testing pending.
