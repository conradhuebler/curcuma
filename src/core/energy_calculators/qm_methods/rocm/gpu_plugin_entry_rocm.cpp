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

extern "C" {

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
