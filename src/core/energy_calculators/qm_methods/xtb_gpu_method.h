/*
 * <Native xTB GPU Method Wrapper — GFN1 / GFN2 on CUDA>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0.
 *
 * Claude Generated (2026-06): ComputationalMethod adapter for the GPU path of
 * native GFN1/GFN2. Selected by the factory when `-gpu cuda` (or `-gpu auto`)
 * is given and the build has USE_CUDA.
 *
 * Design: OWNS a NativeXtbMethod (the validated CPU pipeline: config, large-system
 * modes, error handling, property accessors) + an XtbGpuContext (cuSOLVER/cuBLAS handles
 * + stream). Not a subclass of XTB. Every ComputationalMethod call forwards to the CPU
 * pipeline; the GPU is reached through the hooks the constructor installs on the owned
 * XTB. That shared wrapper logic lives ONCE in XtbGpuAdapter (xtb_gpu_adapter.h, Claude
 * Generated Sep 2026), together with the ROCm and Vulkan wrappers. The CUDA
 * GpuScfBackend, however, stays hand-written in the .cpp: it drives device stages the
 * other two backends do not have (the fused resident SCF loop, the device potential /
 * solvation build, atomic + shell charges, the device SCC energy).
 *
 * Usage:
 *   ./curcuma -sp mol.xyz -method gfn2 -gpu cuda    # explicit GPU
 *   ./curcuma -sp mol.xyz -method gfn1 -gpu auto    # GPU if available
 *   ./curcuma -sp mol.xyz -method gfn2              # CPU (default)
 */

#pragma once

#ifdef USE_CUDA

#include "xtb_gpu_adapter.h"        // XtbGpuAdapter (shared wrapper logic)
#include "cuda/xtb_gpu_context.h"   // XtbGpuContext (pimpl — no CUDA headers)

#include <memory>

/**
 * @brief GPU-accelerated native GFN1/GFN2 (method "gfn1"/"gfn2" with -gpu cuda).
 *
 * A thin instantiation: the whole ComputationalMethod surface (forwarding, device
 * handshake + logging, CPU fallback, gpuActive()) comes from the shared adapter. CUDA
 * adds the cuSOLVER eigensolver hook, the device-resident SCF backend and one behavioural
 * override — setThreadCount, which defaults the remaining HOST work to several cores.
 * Claude Generated.
 */
class XtbGpuComputationalMethod
    : public curcuma::xtb::gpu::XtbGpuAdapter<curcuma::xtb::gpu::XtbGpuContext> {
public:
    explicit XtbGpuComputationalMethod(curcuma::xtb::MethodType method,
                                       const json& config = json{});
    ~XtbGpuComputationalMethod() override;

    /// CUDA-specific: default the host-side work to several cores (see the .cpp).
    void setThreadCount(int threads) override;
};

/// Build a device-resident SCF+gradient backend over an existing context
/// (Stage 4 validation, Claude Generated). Lets a bare XTB be driven on the GPU
/// path via XTB::setGpuScfBackend without the full XtbGpuComputationalMethod.
std::unique_ptr<curcuma::xtb::GpuScfBackend>
createXtbGpuScfBackend(curcuma::xtb::gpu::XtbGpuContext* ctx);

#endif // USE_CUDA
