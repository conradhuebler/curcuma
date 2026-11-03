/*
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 */
// Stand-in GPU plugin for the ABI version check (ctest gpu_plugin_abi) - Claude Generated (Nov 2026, S0).
// Built several times with different definitions (test_cases/CMakeLists.txt):
//   STUB_BACKEND   name in the exported symbols (curcuma_<STUB_BACKEND>_*)
//   STUB_ABI       exported plugin ABI version; undefined = no curcuma_<b>_abi_version symbol (pre-S0 plugin)
//   STUB_ALIGN     exported Eigen alignment; undefined = the value of gpu_plugin_abi.h's user (EIGEN_MAX_ALIGN_BYTES)
// It has no GPU code and no create_* entry points: the loader must decide from the version symbol alone.
#include "src/core/energy_calculators/gpu_plugin_abi.h"
#include "src/core/curcuma_eigen_config.h"
#include <Eigen/Core>

#define STUB_CAT2(a, b, c) a##b##c
#define STUB_CAT(a, b, c) STUB_CAT2(a, b, c)
#define STUB_SYM(suffix) STUB_CAT(curcuma_, STUB_BACKEND, suffix)

extern "C" {
#ifdef STUB_ABI
int STUB_SYM(_abi_version)() { return STUB_ABI; }
#endif
#ifdef STUB_ALIGN
int STUB_SYM(_eigen_align)() { return STUB_ALIGN; }
#else
int STUB_SYM(_eigen_align)() { return EIGEN_MAX_ALIGN_BYTES; }
#endif
int STUB_SYM(_device_count)() { return 0; }
}
