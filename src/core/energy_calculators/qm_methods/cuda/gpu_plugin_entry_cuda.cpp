/*
 * < CUDA GPU backend plugin — C entry points >
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

// Claude Generated (Jul 2026): the C ABI the gpu_plugin loader (core) dlsym's from
// libcurcuma_cuda.so. These functions construct the concrete CUDA ComputationalMethods
// (which live in this plugin .so together with the cuBLAS/cuSOLVER dependencies) and hand
// back a raw pointer the core wraps in a unique_ptr. The config crosses the ABI as a JSON
// string so the boundary carries no C++ container layout. See gpu_plugin.h.

#include "src/core/curcuma_logger.h"
#include "../gfnff_gpu_method.h"      // GFNFFGPUComputationalMethod (global namespace)
#include "../xtb_gpu_method.h"        // XtbGpuComputationalMethod (global namespace)
#include "../xtb_native.h"            // curcuma::xtb::MethodType

#include "json.hpp"

extern "C" {

// method_type: curcuma::xtb::MethodType (GFN1=1, GFN2=2). Returns a new heap object the
// caller (core) owns, or nullptr on failure so the core falls back to the CPU method.
ComputationalMethod* curcuma_cuda_create_native_xtb(int method_type, const char* config_json)
{
    try {
        const json config = json::parse(config_json ? config_json : "{}");
        return new XtbGpuComputationalMethod(
            static_cast<curcuma::xtb::MethodType>(method_type), config);
    } catch (const std::exception& e) {
        CurcumaLogger::error(std::string("CUDA plugin: native-xTB construction failed: ") + e.what());
        return nullptr;
    }
}

ComputationalMethod* curcuma_cuda_create_gfnff(const char* config_json)
{
    try {
        const json config = json::parse(config_json ? config_json : "{}");
        return new GFNFFGPUComputationalMethod("gfnff", config);
    } catch (const std::exception& e) {
        CurcumaLogger::error(std::string("CUDA plugin: GFN-FF construction failed: ") + e.what());
        return nullptr;
    }
}

} // extern "C"
