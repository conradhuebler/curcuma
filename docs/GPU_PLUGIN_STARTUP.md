# CUDA GPU as a runtime plugin — fast CPU startup

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
- `method_factory.cpp` — the 3 `-gpu cuda` construction sites now call the loader.
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

**ROCm/Vulkan** still link at compile time (they are OFF in the default build); the same
plugin pattern generalizes to them when their startup cost matters. See
[MOR41_CPU_GPU_GXTB_EVAL.md](MOR41_CPU_GPU_GXTB_EVAL.md) for the CPU-vs-gxtb benchmark
this work came from.
