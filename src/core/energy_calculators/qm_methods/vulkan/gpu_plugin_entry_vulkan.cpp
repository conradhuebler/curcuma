// Claude Generated (Sep 2026): C ABI of the Vulkan plugin (libcurcuma_vulkan.so), the
// twin of cuda/gpu_plugin_entry_cuda.cpp. The core dlopen's this library on demand
// (gpu_plugin.cpp), so the curcuma executable stays free of the Vulkan loader.
// Vulkan implements the native-xTB backend only; GFN-FF compute shaders are not ported,
// so the GFN-FF entry point reports that and returns nullptr (the core falls back to CPU).
#include "src/core/curcuma_logger.h"
#include "../xtb_vulkan_method.h"     // XtbVulkanComputationalMethod
#include "../xtb_native.h"            // curcuma::xtb::MethodType
#include "json.hpp"

extern "C" {

ComputationalMethod* curcuma_vulkan_create_native_xtb(int method_type, const char* config_json)
{
    try {
        const json config = json::parse(config_json ? config_json : "{}");
        return new XtbVulkanComputationalMethod(
            static_cast<curcuma::xtb::MethodType>(method_type), config);
    } catch (const std::exception& e) {
        CurcumaLogger::error(std::string("Vulkan plugin: native-xTB construction failed: ") + e.what());
        return nullptr;
    }
}

ComputationalMethod* curcuma_vulkan_create_gfnff(const char* /*config_json*/)
{
    CurcumaLogger::warn("GFN-FF -gpu vulkan: compute shaders not yet ported; using CPU.");
    return nullptr;
}

} // extern "C"
