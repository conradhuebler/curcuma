// Claude Generated (Sep 2026): C ABI of the Vulkan plugin (libcurcuma_vulkan.so), the
// twin of cuda/gpu_plugin_entry_cuda.cpp. The core dlopen's this library on demand
// (gpu_plugin.cpp), so the curcuma executable stays free of the Vulkan loader.
// Vulkan implements the native-xTB backend only; GFN-FF compute shaders are not ported,
// so the GFN-FF entry point reports that and returns nullptr (the core falls back to CPU).
#include "src/core/curcuma_logger.h"
#include "src/core/energy_calculators/gpu_plugin_abi.h"
#include "../xtb_vulkan_method.h"     // XtbVulkanComputationalMethod
#include "../xtb_native.h"            // curcuma::xtb::MethodType
#include "json.hpp"
#include "vk_context.h"

#include <cstring>

// Plugin ABI version (gpu_plugin_abi.h); the loader refuses this plugin if it differs from the executable's
// (Claude Generated, Nov 2026, S0).
CURCUMA_GPU_PLUGIN_DEFINE_ABI_VERSION(vulkan)

extern "C" {

// ---- Multi-GPU device discovery (Claude Generated, Sep 2026) -------------------------------
// Indices are vkEnumeratePhysicalDevices positions, the same numbers `gpu_device` selects.
int curcuma_vulkan_device_count()
{
    return static_cast<int>(curcuma::vk::VkContext::enumerateDevices().size());
}

int curcuma_vulkan_device_info(int index, char* buf, int len)
{
    const auto devs = curcuma::vk::VkContext::enumerateDevices();
    if (index < 0 || index >= static_cast<int>(devs.size()))
        return -1;
    json j;
    j["index"] = index;
    j["name"] = devs[index].name;
    j["memory_total_bytes"] = devs[index].memory_bytes;
    j["usable"] = devs[index].usable;
    j["discrete"] = devs[index].discrete;
    const std::string s = j.dump();
    if (buf && len > 0) {
        const size_t n = std::min(s.size(), static_cast<size_t>(len - 1));
        std::memcpy(buf, s.data(), n);
        buf[n] = '\0';
    }
    return static_cast<int>(s.size());
}

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
