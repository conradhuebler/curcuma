/*
 * < GPU backend plugin loader >
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 */

// Claude Generated (Jul 2026): runtime loader for the GPU backend plugins.
//
// The CUDA/ROCm/Vulkan backends are built into separate shared libraries
// (libcurcuma_cuda.so, ...) that are dlopen'd ON DEMAND — only when a GPU run is
// actually requested. This keeps the heavy GPU runtime libraries (cuBLAS/cuSOLVER,
// which run ~40 ms of DT_INIT table-building at load) OUT of the main binary's
// startup path, so CPU-only runs start in ~1 ms instead of ~50 ms.
//
// The plugin exposes C entry points (curcuma_<backend>_create_native_xtb /
// curcuma_<backend>_create_gfnff) that construct the concrete GPU ComputationalMethod
// and hand back a raw pointer the core wraps in a unique_ptr. The plugin .so is kept
// loaded for the process lifetime, so the returned object's vtable/destructor stay valid.

#pragma once

#include <memory>
#include <string>

#include "computational_method.h"   // ComputationalMethod, json (global namespace)

namespace gpu_plugin {

/**
 * Construct a native-xTB (GFN1/GFN2) GPU method from the given backend plugin.
 * @param backend    "cuda" | "rocm" | "vulkan"
 * @param method_type curcuma::xtb::MethodType value (GFN1=1, GFN2=2)
 * Returns nullptr if the plugin is unavailable or construction fails (caller must
 * fall back to the CPU method). Claude Generated.
 */
std::unique_ptr<ComputationalMethod> createNativeXtb(const std::string& backend,
                                                     int method_type, const json& config);

/**
 * Construct a GFN-FF GPU method from the given backend plugin.
 * Returns nullptr on failure (caller falls back to CPU). Claude Generated.
 */
std::unique_ptr<ComputationalMethod> createGfnff(const std::string& backend, const json& config);

} // namespace gpu_plugin
