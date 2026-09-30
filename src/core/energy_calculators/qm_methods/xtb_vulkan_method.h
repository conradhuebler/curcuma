/*
 * <Native xTB Vulkan Method Wrapper — GFN1 / GFN2 on Vulkan>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0.
 *
 * Claude Generated (2026-06): ComputationalMethod adapter for the Vulkan path of
 * native GFN1/GFN2 — the Vulkan sibling of XtbGpuComputationalMethod (CUDA) /
 * XtbHipComputationalMethod (ROCm). Selected by the factory when `-gpu vulkan` (or
 * `-gpu auto` on a Vulkan-only build) is given and the build has USE_VULKAN.
 *
 * Design: OWNS a NativeXtbMethod (validated CPU pipeline) + an XtbVulkanContext (device
 * handles); every ComputationalMethod call forwards to the CPU pipeline, and the GPU is
 * reached through the hooks the constructor installs on the owned XTB. That shared
 * wrapper logic lives ONCE in XtbGpuAdapter (xtb_gpu_adapter.h, Claude Generated Sep
 * 2026); only the Vulkan-specific eigensolver hook, the resident backend and the
 * mixed-precision policy are in the .cpp.
 *
 * Usage:
 *   ./curcuma -sp mol.xyz -method gfn2 -gpu vulkan   # explicit Vulkan
 *   ./curcuma -sp mol.xyz -method gfn1 -gpu auto     # GPU if available
 *   ./curcuma -sp mol.xyz -method gfn2               # CPU (default)
 */

#pragma once

#ifdef USE_VULKAN

#include "xtb_gpu_adapter.h"         // XtbGpuAdapter (shared wrapper logic)
#include "vulkan/xtb_vulkan_context.h"  // XtbVulkanContext (pimpl — no Vulkan headers)

/**
 * @brief GPU-accelerated native GFN1/GFN2 (method "gfn1"/"gfn2" with -gpu vulkan).
 *
 * A thin instantiation: the whole ComputationalMethod surface (forwarding, device
 * handshake + logging, CPU fallback, gpuActive()) comes from the shared adapter; the
 * constructor adds the Vulkan Löwdin/Jacobi eigensolver hook and the device-resident
 * SCF backend. Claude Generated.
 */
class XtbVulkanComputationalMethod
    : public curcuma::xtb::gpu::XtbGpuAdapter<curcuma::xtb::gpu::XtbVulkanContext> {
public:
    explicit XtbVulkanComputationalMethod(curcuma::xtb::MethodType method,
                                          const json& config = json{});
    ~XtbVulkanComputationalMethod() override;
};

#endif // USE_VULKAN
