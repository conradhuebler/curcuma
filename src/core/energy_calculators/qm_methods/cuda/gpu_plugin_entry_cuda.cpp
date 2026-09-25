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

#include <cuda_runtime.h>

#include <cstring>
#include "xtb_distributed_eigensolver.h"   // availableBackends() for the -methods query

namespace {
// Copy a JSON dump into the caller's buffer. Returns the length needed (excluding the NUL),
// so a caller with a too-small buffer can retry. Claude Generated (Sep 2026, multi-GPU).
int writeJson(const json& j, char* buf, int len)
{
    const std::string s = j.dump();
    if (buf && len > 0) {
        const size_t n = std::min(s.size(), static_cast<size_t>(len - 1));
        std::memcpy(buf, s.data(), n);
        buf[n] = '\0';
    }
    return static_cast<int>(s.size());
}
} // namespace

extern "C" {

// ---- Multi-GPU device discovery (Claude Generated, Sep 2026) ----------------------------
// Optional symbols: the core treats a plugin without them as "1 device, no details".
// Indices are the CUDA runtime's, i.e. AFTER CUDA_VISIBLE_DEVICES (which SLURM sets per job).

int curcuma_cuda_device_count()
{
    int n = 0;
    return cudaGetDeviceCount(&n) == cudaSuccess ? n : 0;
}

// Claude Generated (Sep 2026): which distributed-eigensolver backends libcurcuma_cuda_mgpu.so
// was built with - "mp" (cuSOLVERMp + NCCL), "mg" (the deprecated cusolverMg fallback, measured
// 15x slower than one GPU on polymer_2x), "mp,mg", or "none (<why>)" when the library is
// missing. Queried by `curcuma -methods`, so a cluster build can be checked without running a
// calculation. Loading the library here is the same dlopen the SCF would do.
const char* curcuma_cuda_mgpu_backends()
{
    static const std::string s = curcuma::xtb::gpu::DistributedEigensolver::availableBackends();
    return s.c_str();
}

// JSON: {"index","name","memory_total_bytes","compute_capability","pci_bus_id","pci_device_id"}.
// Reads device properties only - no context is created, so this is cheap and has no side effect
// on the device memory of the calling process.
int curcuma_cuda_device_info(int index, char* buf, int len)
{
    cudaDeviceProp prop{};
    if (cudaGetDeviceProperties(&prop, index) != cudaSuccess)
        return -1;
    json j;
    j["index"] = index;
    j["name"] = std::string(prop.name);
    j["memory_total_bytes"] = static_cast<std::uint64_t>(prop.totalGlobalMem);
    j["compute_capability"] = std::to_string(prop.major) + "." + std::to_string(prop.minor);
    j["pci_bus_id"] = prop.pciBusID;
    j["pci_device_id"] = prop.pciDeviceID;
    return writeJson(j, buf, len);
}

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
