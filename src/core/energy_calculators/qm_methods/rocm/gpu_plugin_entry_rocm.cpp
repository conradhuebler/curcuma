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

namespace {
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
    j["memory_total_bytes"] = static_cast<std::uint64_t>(prop.totalMemory);
    j["compute_capability"] = std::string(prop.gcnArchName);
    return writeJson(j, buf, len);
}

// method_type: curcuma::xtb::MethodType (GFN1=1, GFN2=2). nullptr on failure -> CPU fallback.
ComputationalMethod* curcuma_rocm_create_native_xtb(int method_type, const char* config_json)
{
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
    try {
        const json config = json::parse(config_json ? config_json : "{}");
        return new GFNFFHipComputationalMethod("gfnff", config);
    } catch (const std::exception& e) {
        CurcumaLogger::error(std::string("ROCm plugin: GFN-FF construction failed: ") + e.what());
        return nullptr;
    }
}

} // extern "C"
