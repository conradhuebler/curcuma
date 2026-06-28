/*
 * <Vulkan Compute Context — instance / device / compute queue owner>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0.
 *
 * Claude Generated (2026-06): Generic Vulkan compute context shared by the GFN-FF and
 * native-xTB Vulkan backends. Owns a VkInstance, the selected VkPhysicalDevice, a
 * logical VkDevice with one compute queue, and a command pool. FP64 (shaderFloat64)
 * is REQUIRED for the SCF numerics — ok() is false (→ CPU fallback) if no suitable
 * device with shaderFloat64 is found.
 *
 * A pimpl keeps <vulkan/vulkan.h> out of the host translation units that merely hold a
 * context. Device handles are exposed as opaque void* (cast back to the Vulkan types
 * inside the .cpp engines) so this header stays Vulkan-header-free.
 *
 * Stage 0 (current): instance + device + queue + command-pool handshake only. The
 * compute pipelines (GEMM, TRSM, integral/Fock/density/gradient kernels, the cyclic
 * Jacobi symmetric eigensolver) are added in later stages.
 */

#pragma once

#ifdef USE_VULKAN

#include <cstdint>
#include <memory>
#include <string>

namespace curcuma {
namespace vk {

class VkContext {
public:
    VkContext();
    ~VkContext();

    VkContext(const VkContext&) = delete;
    VkContext& operator=(const VkContext&) = delete;

    /// True once instance + device + compute queue + command pool are created and the
    /// device advertises shaderFloat64.
    bool ok() const;

    /// True if the selected device supports FP64 in shaders (required; ok() implies it).
    bool hasFloat64() const;

    /// Selected device name (e.g. "AMD Radeon RX 7900 XTX" / "NVIDIA ..."); empty if none.
    std::string deviceName() const;

    /// Index of the selected physical device in vkEnumeratePhysicalDevices, or -1.
    int deviceId() const;

    /// True if a Vulkan instance can be created and at least one device with a compute
    /// queue + shaderFloat64 is visible (static probe; creates and tears down a probe
    /// instance). Claude Generated.
    static bool deviceAvailable();

    // ---- Opaque device handles (Stage 1+ engines cast these back) -----------
    // Returned as void* so this header needs no <vulkan/vulkan.h>. nullptr when !ok().
    void*    instance() const;             ///< VkInstance
    void*    physicalDevice() const;       ///< VkPhysicalDevice
    void*    device() const;               ///< VkDevice
    void*    computeQueue() const;         ///< VkQueue
    uint32_t computeQueueFamily() const;   ///< queue-family index of computeQueue()
    void*    commandPool() const;          ///< VkCommandPool (compute, resettable)

private:
    struct Impl;
    std::unique_ptr<Impl> m_impl;
};

} // namespace vk
} // namespace curcuma

#endif // USE_VULKAN
