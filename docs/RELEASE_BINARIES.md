# Prebuilt binaries: download, GPU plugins, CI

Built by `.github/workflows/ccpp.yml` (AutomaticBuild), which calls `linux-gpu.yml` for Linux. Published as pre-releases on the [Releases page](https://github.com/conradhuebler/curcuma/releases): `Curcuma 0.0.<n>` from `master`, `Curcuma CI (<branch>)` for every `feature/**` and `fix/**` branch (one release per branch, assets named with the run number; older assets of a branch release are not removed automatically).

## Files

| File | Content | In which releases |
|---|---|---|
| `curcuma-<v>-avx2-x86_64.AppImage`, `curcuma-<v>-avx512-x86_64.AppImage` | CPU-only program, libraries bundled | all |
| `curcuma-<v>-linux-<level>.tar.gz` | the same program as an archive, top directory `curcuma-<level>/` | all |
| `curcuma-<v>-plugin-cuda-<level>.tar.gz` | `libcurcuma_cuda.so` and the CUDA libraries (about 1.1 GB) | manual runs only |
| `curcuma-<v>-plugin-rocm-<level>.tar.gz` | `libcurcuma_rocm.so`, HIP, rocBLAS, rocSOLVER and rocBLAS kernels (about 1.0 GB) | manual runs only |
| `curcuma-<v>-x86_64-Windows.zip` | Windows build (baseline, no AVX2/AVX512, no GPU) | all |

`<level>` is `avx512` if `grep -c avx512f /proc/cpuinfo` is above 0, otherwise `avx2`. Linux builds need glibc 2.35 or newer. Vulkan is not shipped.

## Using the GPU plugins

1. Download the AppImage and the plugin archive of the **same level**. The loader compares the Eigen alignment (`curcuma_<backend>_eigen_align()`) and refuses a plugin of the other level with a warning, then runs on the CPU.
2. Unpack the plugin archive next to the AppImage: `tar -xzf curcuma-<v>-plugin-cuda-<level>.tar.gz --strip-components=1` (the archives have the top directory `curcuma-<level>/`; with the executable archive, unpack everything into one directory without the option).
3. Run `./curcuma-<v>-<level>-x86_64.AppImage -sp mol.xyz -method gfn2 -gpu cuda` (or `-gpu rocm`).

Plugin search order: `$CURCUMA_PLUGIN_DIR`, the executable's directory, the directory of the `.AppImage` file, `~/.local/share/curcuma/plugins`, then the system library path. Details: [GPU_PLUGIN_STARTUP.md](GPU_PLUGIN_STARTUP.md).

Requirements: CUDA needs the NVIDIA driver (`libcuda` is taken from the system). ROCm needs the `amdgpu` kernel driver (`/dev/kfd`); the plugin is built for `gfx1100`, `gfx1150` and `gfx1200` and falls back to the CPU on other cards. `-gpu_strict true` stops at the first CPU fallback.

## Getting the plugins from a branch

Plugin archives are not built on pushes, neither on branches nor on `master` (two GPU builds take about two hours). Start Actions > AutomaticBuild or CiVariants with "Run workflow" on the branch; the plugins are in the artifacts `linux-avx2` and `linux-avx512` (login required, kept 90 days), or `gh run download <run id> -n linux-avx512`.

## Status

AI-generated CI. Checked by CI: the packages build, resolve all libraries in a clean Ubuntu 22.04 (`ldd`), the CPU-only program starts, and the AVX2 AppImage loads the CUDA plugin placed next to it on a runner without GPU. Operator report (Oct 30, 2026): both CPU levels (AVX2 and AVX512) with the CUDA and the ROCm plugin placed next to the program work, each on a different machine; which machine ran which combination, and the GPU models, are not recorded. Not checked: other distributions, the rejection of a plugin of the wrong level, ROCm cards outside `gfx1100`, `gfx1150`, `gfx1200`.
