// Claude Generated (Sep 2026): C ABI of the ROCm/HIP plugin (libcurcuma_rocm.so), the
// twin of cuda/gpu_plugin_entry_cuda.cpp. The core dlopen's this library on demand
// (gpu_plugin.cpp), so the curcuma executable stays free of HIP/rocBLAS/rocSOLVER
// symbols and their DT_INIT cost. Each function constructs the concrete HIP
// ComputationalMethod and hands back a raw pointer the core wraps in a unique_ptr;
// the config crosses the ABI as a JSON string. See gpu_plugin.h.
#include "src/core/curcuma_logger.h"
#include "../gfnff_hip_method.h"      // GFNFFHipComputationalMethod
#include "../xtb_hip_method.h"        // XtbHipComputationalMethod
#include "../xtb_native.h"            // curcuma::xtb::MethodType
#include "json.hpp"

#include <hip/hip_runtime.h>

#include <cstring>
#include <string>

#ifndef CURCUMA_ROCM_BUILT_ARCHS
#define CURCUMA_ROCM_BUILT_ARCHS ""
#endif

namespace {
/**
 * Claude Generated: true if every visible device has a code object in this build.
 * A device whose architecture was not compiled in makes the HIP runtime abort in
 * hip::StatCO::getStatFunc at the first kernel launch, so the plugin must refuse before.
 * gcnArchName looks like "gfx1150" or "gfx90a:sramecc+:xnack-"; only the part before the
 * first ':' names the code object. All visible devices are checked because the device
 * index is chosen further down (multi-GPU).
 */
bool allDevicesHaveCodeObject()
{
    static const std::string built = std::string(",") + CURCUMA_ROCM_BUILT_ARCHS + ",";
    int n = 0;
    if (hipGetDeviceCount(&n) != hipSuccess || n <= 0)
        return false;
    for (int i = 0; i < n; ++i) {
        hipDeviceProp_t prop;
        if (hipGetDeviceProperties(&prop, i) != hipSuccess)
            return false;
        const std::string arch = std::string(prop.gcnArchName).substr(0, std::string(prop.gcnArchName).find(':'));
        if (built.find("," + arch + ",") == std::string::npos) {
            CurcumaLogger::warn("ROCm plugin: device " + std::to_string(i) + " (" + prop.name + ") is " + arch
                + ", this build contains code for " + std::string(CURCUMA_ROCM_BUILT_ARCHS)
                + ". Using the CPU. Rebuild with -DROCM_GPU_ARCH=" + arch + ".");
            return false;
        }
    }
    return true;
}

// Copy a JSON dump into the caller's buffer; returns the needed length. Claude Generated.
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

// ---- Multi-GPU device discovery (Claude Generated, Sep 2026; unverified - no ROCm SDK here) --
// Indices are after ROCR_VISIBLE_DEVICES / HIP_VISIBLE_DEVICES.
int curcuma_rocm_device_count()
{
    int n = 0;
    return hipGetDeviceCount(&n) == hipSuccess ? n : 0;
}

int curcuma_rocm_device_info(int index, char* buf, int len)
{
    hipDeviceProp_t prop;
    if (hipGetDeviceProperties(&prop, index) != hipSuccess)
        return -1;
    json j;
    j["index"] = index;
    j["name"] = std::string(prop.name);
    // hipDeviceProp_t mirrors cudaDeviceProp: the member is totalGlobalMem, there is no
    // totalMemory (which is what the CUDA twin in cuda/gpu_plugin_entry_cuda.cpp uses).
    j["memory_total_bytes"] = static_cast<std::uint64_t>(prop.totalGlobalMem);
    j["compute_capability"] = std::string(prop.gcnArchName);
    return writeJson(j, buf, len);
}

// method_type: curcuma::xtb::MethodType (GFN1=1, GFN2=2). nullptr on failure -> CPU fallback.
ComputationalMethod* curcuma_rocm_create_native_xtb(int method_type, const char* config_json)
{
    if (!allDevicesHaveCodeObject())
        return nullptr;
    try {
        const json config = json::parse(config_json ? config_json : "{}");
        return new XtbHipComputationalMethod(
            static_cast<curcuma::xtb::MethodType>(method_type), config);
    } catch (const std::exception& e) {
        CurcumaLogger::error(std::string("ROCm plugin: native-xTB construction failed: ") + e.what());
        return nullptr;
    }
}

ComputationalMethod* curcuma_rocm_create_gfnff(const char* config_json)
{
    if (!allDevicesHaveCodeObject())
        return nullptr;
    try {
        const json config = json::parse(config_json ? config_json : "{}");
        return new GFNFFHipComputationalMethod("gfnff", config);
    } catch (const std::exception& e) {
        CurcumaLogger::error(std::string("ROCm plugin: GFN-FF construction failed: ") + e.what());
        return nullptr;
    }
}

} // extern "C"
