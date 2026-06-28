/*
 * <Vulkan Compute Context — implementation>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0.
 *
 * Claude Generated (2026-06): Stage-0 device handshake for the Vulkan compute path —
 * VkInstance + VkPhysicalDevice (discrete GPU preferred, shaderFloat64 required) +
 * one compute queue + a resettable command pool. No validation layers (production).
 * Every failure path leaves ok() == false so the caller falls back to the CPU.
 */

#ifdef USE_VULKAN

#include "vk_context.h"

#include <vulkan/vulkan.h>

#include <cstring>
#include <vector>

namespace curcuma {
namespace vk {

struct VkContext::Impl {
    VkInstance       instance     = VK_NULL_HANDLE;
    VkPhysicalDevice phys         = VK_NULL_HANDLE;
    VkDevice         device       = VK_NULL_HANDLE;
    VkQueue          queue        = VK_NULL_HANDLE;
    VkCommandPool    cmdPool      = VK_NULL_HANDLE;
    uint32_t         queueFamily  = 0;
    int              deviceIndex  = -1;
    bool             fp64         = false;
    bool             ok           = false;
    std::string      name;
};

namespace {

/// Create a bare compute-capable instance (no layers, no extensions). Returns
/// VK_NULL_HANDLE on failure.
VkInstance createInstance()
{
    VkApplicationInfo app{};
    app.sType              = VK_STRUCTURE_TYPE_APPLICATION_INFO;
    app.pApplicationName   = "curcuma";
    app.applicationVersion = VK_MAKE_API_VERSION(0, 1, 0, 0);
    app.pEngineName        = "curcuma";
    app.engineVersion      = VK_MAKE_API_VERSION(0, 1, 0, 0);
    app.apiVersion         = VK_API_VERSION_1_1;

    VkInstanceCreateInfo ci{};
    ci.sType            = VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO;
    ci.pApplicationInfo = &app;

    VkInstance inst = VK_NULL_HANDLE;
    if (vkCreateInstance(&ci, nullptr, &inst) != VK_SUCCESS)
        return VK_NULL_HANDLE;
    return inst;
}

/// Find a queue family supporting compute on `phys`. Returns true + the family index
/// on success. Prefers a compute-only family (no graphics) for a dedicated queue.
bool findComputeQueueFamily(VkPhysicalDevice phys, uint32_t& family_out)
{
    uint32_t n = 0;
    vkGetPhysicalDeviceQueueFamilyProperties(phys, &n, nullptr);
    if (n == 0) return false;
    std::vector<VkQueueFamilyProperties> props(n);
    vkGetPhysicalDeviceQueueFamilyProperties(phys, &n, props.data());

    int any_compute = -1;
    for (uint32_t i = 0; i < n; ++i) {
        if (!(props[i].queueFlags & VK_QUEUE_COMPUTE_BIT)) continue;
        if (any_compute < 0) any_compute = static_cast<int>(i);
        // Prefer a dedicated compute family (compute without graphics).
        if (!(props[i].queueFlags & VK_QUEUE_GRAPHICS_BIT)) {
            family_out = i;
            return true;
        }
    }
    if (any_compute >= 0) {
        family_out = static_cast<uint32_t>(any_compute);
        return true;
    }
    return false;
}

/// Score a candidate device: must have compute + shaderFloat64. Discrete GPUs rank
/// highest. Returns a score (<0 = unusable).
int scoreDevice(VkPhysicalDevice phys, uint32_t& family_out)
{
    VkPhysicalDeviceFeatures feats{};
    vkGetPhysicalDeviceFeatures(phys, &feats);
    if (!feats.shaderFloat64) return -1;
    if (!findComputeQueueFamily(phys, family_out)) return -1;

    VkPhysicalDeviceProperties prop{};
    vkGetPhysicalDeviceProperties(phys, &prop);
    int score = 1;
    if (prop.deviceType == VK_PHYSICAL_DEVICE_TYPE_DISCRETE_GPU) score = 1000;
    else if (prop.deviceType == VK_PHYSICAL_DEVICE_TYPE_INTEGRATED_GPU) score = 100;
    return score;
}

/// Pick the best usable physical device on `inst`. Returns true + handle/index/family.
bool pickPhysicalDevice(VkInstance inst, VkPhysicalDevice& phys_out, int& index_out,
                        uint32_t& family_out)
{
    uint32_t n = 0;
    if (vkEnumeratePhysicalDevices(inst, &n, nullptr) != VK_SUCCESS || n == 0)
        return false;
    std::vector<VkPhysicalDevice> devs(n);
    if (vkEnumeratePhysicalDevices(inst, &n, devs.data()) != VK_SUCCESS)
        return false;

    int best_score = -1;
    for (uint32_t i = 0; i < n; ++i) {
        uint32_t fam = 0;
        const int s = scoreDevice(devs[i], fam);
        if (s > best_score) {
            best_score = s;
            phys_out   = devs[i];
            index_out  = static_cast<int>(i);
            family_out = fam;
        }
    }
    return best_score >= 0;
}

} // namespace

VkContext::VkContext()
    : m_impl(std::make_unique<Impl>())
{
    m_impl->instance = createInstance();
    if (m_impl->instance == VK_NULL_HANDLE) return;

    if (!pickPhysicalDevice(m_impl->instance, m_impl->phys, m_impl->deviceIndex,
                            m_impl->queueFamily))
        return;

    VkPhysicalDeviceProperties prop{};
    vkGetPhysicalDeviceProperties(m_impl->phys, &prop);
    m_impl->name = prop.deviceName;
    m_impl->fp64 = true;  // pickPhysicalDevice only accepts shaderFloat64 devices

    // Logical device with the one compute queue and shaderFloat64 enabled.
    const float priority = 1.0f;
    VkDeviceQueueCreateInfo qci{};
    qci.sType            = VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO;
    qci.queueFamilyIndex = m_impl->queueFamily;
    qci.queueCount       = 1;
    qci.pQueuePriorities = &priority;

    VkPhysicalDeviceFeatures enable{};
    enable.shaderFloat64 = VK_TRUE;

    VkDeviceCreateInfo dci{};
    dci.sType                = VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO;
    dci.queueCreateInfoCount = 1;
    dci.pQueueCreateInfos    = &qci;
    dci.pEnabledFeatures     = &enable;

    if (vkCreateDevice(m_impl->phys, &dci, nullptr, &m_impl->device) != VK_SUCCESS)
        return;

    vkGetDeviceQueue(m_impl->device, m_impl->queueFamily, 0, &m_impl->queue);

    VkCommandPoolCreateInfo pci{};
    pci.sType            = VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO;
    pci.flags            = VK_COMMAND_POOL_CREATE_RESET_COMMAND_BUFFER_BIT;
    pci.queueFamilyIndex = m_impl->queueFamily;
    if (vkCreateCommandPool(m_impl->device, &pci, nullptr, &m_impl->cmdPool) != VK_SUCCESS)
        return;

    m_impl->ok = true;
}

VkContext::~VkContext()
{
    if (!m_impl) return;
    if (m_impl->cmdPool != VK_NULL_HANDLE)
        vkDestroyCommandPool(m_impl->device, m_impl->cmdPool, nullptr);
    if (m_impl->device != VK_NULL_HANDLE)
        vkDestroyDevice(m_impl->device, nullptr);
    if (m_impl->instance != VK_NULL_HANDLE)
        vkDestroyInstance(m_impl->instance, nullptr);
}

bool VkContext::ok() const { return m_impl && m_impl->ok; }
bool VkContext::hasFloat64() const { return m_impl && m_impl->fp64; }
std::string VkContext::deviceName() const { return m_impl ? m_impl->name : std::string(); }
int VkContext::deviceId() const { return m_impl ? m_impl->deviceIndex : -1; }

bool VkContext::deviceAvailable()
{
    VkInstance inst = createInstance();
    if (inst == VK_NULL_HANDLE) return false;
    VkPhysicalDevice phys = VK_NULL_HANDLE;
    int idx = -1;
    uint32_t fam = 0;
    const bool found = pickPhysicalDevice(inst, phys, idx, fam);
    vkDestroyInstance(inst, nullptr);
    return found;
}

void*    VkContext::instance() const { return m_impl ? static_cast<void*>(m_impl->instance) : nullptr; }
void*    VkContext::physicalDevice() const { return m_impl ? static_cast<void*>(m_impl->phys) : nullptr; }
void*    VkContext::device() const { return m_impl ? static_cast<void*>(m_impl->device) : nullptr; }
void*    VkContext::computeQueue() const { return m_impl ? static_cast<void*>(m_impl->queue) : nullptr; }
uint32_t VkContext::computeQueueFamily() const { return m_impl ? m_impl->queueFamily : 0; }
void*    VkContext::commandPool() const { return m_impl ? static_cast<void*>(m_impl->cmdPool) : nullptr; }

} // namespace vk
} // namespace curcuma

#endif // USE_VULKAN
