# GPU backends as runtime plugins — fast CPU startup, one binary for every machine

> 🤖 AI-generated (Jul 2026), ⚙️ machine-tested (200/200 GPU ctests, 472/472 non-GPU
> ctests, CPU==GPU gradient 5e-10). Not human production-tested.

## Problem

The CUDA math libraries `libcublasLt` (~30 ms) and `libcusolver` (~35 ms) run heavy
`DT_INIT` constructors (building kernel dispatch tables) when the dynamic loader maps
them. Linking them into the main `curcuma` binary taxed **every** run — including
CPU-only ones that never touch the GPU — with ~40 ms of startup. Measured: an empty
`main(){}` linked against them starts in 40 ms; `/bin/true` in 0.5 ms; gxtb in 0.9 ms.
For fast methods this dominated: GFN-FF on a 20-atom molecule was 62 ms wall, of which
~40 ms was this tax and only ~3 ms was actual compute.

## Fix

The CUDA backend is now a **runtime-`dlopen`'d plugin** (`libcurcuma_cuda.so`) instead of
being linked into `curcuma_core`. The core stays CUDA-symbol-free — it already drove the
GPU only through the CUDA-free `GpuScfBackend` seam (`xtb_native.h`) — and loads the
plugin lazily the first time `-gpu cuda` is requested.

- `src/core/energy_calculators/gpu_plugin.{h,cpp}` — the loader: `dlopen`s
  `libcurcuma_cuda.so` (next to the executable via `/proc/self/exe`, then the default
  search path), caches the handle, `dlsym`s the C entry points. Returns `nullptr` on any
  failure so `method_factory` falls back to the CPU method.
- `qm_methods/cuda/gpu_plugin_entry_cuda.cpp` — the plugin's `extern "C"` entry points
  (`curcuma_cuda_create_native_xtb`, `curcuma_cuda_create_gfnff`) construct the concrete
  GPU `ComputationalMethod`. Config crosses the ABI as a JSON string (no C++ container
  layout on the boundary).
- `method_factory.cpp` — the `-gpu` construction sites call the loader (since Sep 2026 for every backend, see below).
- CMake: the CUDA sources + `xtb_gpu_method.cpp` build into `add_library(curcuma_cuda
  SHARED …)` linking `CUDA::*`; `curcuma_core` no longer links CUDA. The plugin's
  undefined core symbols resolve against the `curcuma` executable at load time (the build
  already links everything `-rdynamic`, no `--gc-sections`). Plugin compile definitions
  are matched to the core's via `$<TARGET_PROPERTY:curcuma_core,COMPILE_DEFINITIONS>` so
  every shared class has one ABI on both sides.

## Result

| | before | after |
|---|---|---|
| `curcuma -version` startup | ~50 ms | **~9 ms** |
| main binary CUDA deps (`ldd`) | 6 libs | **0** |
| GFN-FF `-sp` on 20 atoms | 62 ms | **21 ms** |
| GFN-FF `-sp` on 231 atoms | 96 ms | **49 ms** |
| `-gpu cuda` (RTX 5080) | works | works (via plugin), gradient == CPU to 5e-10 |

The residual ~9 ms is MKL library loading (BLAS/eigensolve, used by every SQM run — not
deferrable). GPU is unchanged functionally: the plugin is kept loaded for the process
lifetime, so the returned object's vtable/destructor stay valid.

See [MOR41_CPU_GPU_GXTB_EVAL.md](MOR41_CPU_GPU_GXTB_EVAL.md) for the CPU-vs-gxtb benchmark
this work came from.

## Sep 2026: plugin symmetry — ROCm and Vulkan are plugins too

> 🤖 AI-generated, ⚙️ machine-tested: CPU build (plugin absent → warning + CPU energy
> unchanged), CUDA plugin (`release_cuda`) and Vulkan plugin (`release_vulkan`, RTX 5080)
> reproduce the pre-plugin energies to the printed digit; ROCm could not be compiled here
> (no SDK) — its CMake block mirrors the CUDA one line by line and is **unverified**.

ROCm and Vulkan used to be compiled *into* `curcuma_core` (`target_sources(curcuma_core …)`,
`target_compile_definitions(curcuma_core PUBLIC USE_ROCM/USE_VULKAN)`), so every backend
choice was a different core binary and `method_factory.cpp` / `energycalculator.cpp` carried
a `#if USE_CUDA … #elif USE_ROCM …` ladder. Now all three backends follow the same recipe:

| backend | plugin | entry points (`extern "C"`) | GFN-FF |
|---|---|---|---|
| CUDA | `libcurcuma_cuda.so` | `curcuma_cuda_create_native_xtb`, `curcuma_cuda_create_gfnff` | yes |
| ROCm | `libcurcuma_rocm.so` | `curcuma_rocm_create_native_xtb`, `curcuma_rocm_create_gfnff` | yes |
| Vulkan | `libcurcuma_vulkan.so` | `curcuma_vulkan_create_native_xtb`, `curcuma_vulkan_create_gfnff` (returns `nullptr` + warning: shaders not ported → CPU) | no |

What changed:

- **No backend `#ifdef` in the core.** `USE_CUDA/USE_ROCM/USE_VULKAN` are defined only on
  the plugin targets; `grep USE_CUDA src/core/energy_calculators/*.cpp` is empty. The CPU
  objects are byte-identical whether or not any plugin is configured, and any subset of
  plugins can sit next to the same executable.
- **Runtime dispatch** (`method_factory.cpp` `resolveGpuMode()`): `-gpu cuda|rocm|vulkan`
  probes `gpu_plugin::available(backend)` (a quiet `dlopen`); missing plugin → warning
  *"the plugin libcurcuma_<b>.so is not present … cmake -DUSE_<B>=ON"* and CPU fallback.
  `-gpu auto` takes the first plugin present in the order cuda, rocm, vulkan. Unknown
  values warn and use CPU. `curcuma -methods` lists the plugins found next to the binary.
- **`EnergyCalculator::m_gpu_fallback`** (the "requested GPU but running on CPU" flag the
  capabilities print) uses the same runtime probe.
- **CMake**: `add_library(curcuma_rocm SHARED …)` holds the two `hipcc`-compiled
  `EXTERNAL_OBJECT`s + the g++ host wrappers + the entry file and is the *only* target that
  links `libamdhip64`/rocSOLVER/rocBLAS (rocSOLVER is now unconditionally required — the
  GFN-FF HIP EEQ solve always needed it, the old "Stage 0 without rocSOLVER" branch was
  dead). `add_library(curcuma_vulkan SHARED …)` holds the Vulkan context + wrapper + entry
  and is the only target linking `Vulkan::Vulkan`. Both mirror the core's compile
  definitions via `$<TARGET_PROPERTY:curcuma_core,COMPILE_DEFINITIONS>` for one ABI.
- The two CUDA-only test executables (`test_gfnff_gpu`, `test_gpu_numgrad`) get
  `USE_CUDA` explicitly and `test_gfnff_gpu` links the plugin (it `dynamic_cast`s to
  `GFNFFGPUComputationalMethod`).

Adding a fourth backend now means: one `gpu_plugin_entry_<b>.cpp` with the two C entry
points, one `add_library(curcuma_<b> SHARED …)` block, and appending `"<b>"` to
`gpu_plugin::knownBackends()`. Nothing else in the core changes.
