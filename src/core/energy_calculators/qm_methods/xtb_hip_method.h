/*
 * <Native xTB ROCm/HIP Method Wrapper — GFN1 / GFN2 on ROCm>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0.
 *
 * Claude Generated (2026-06): ComputationalMethod adapter for the ROCm path of
 * native GFN1/GFN2 — the AMD/HIP sibling of XtbGpuComputationalMethod (CUDA).
 * Selected by the factory when `-gpu rocm` (or `-gpu auto` on a ROCm-only build) is
 * given and the build has USE_ROCM.
 *
 * Design: OWNS a NativeXtbMethod (the validated CPU pipeline: config, large-system
 * modes, error handling, property accessors) + an XtbHipContext (hipBLAS/rocSOLVER
 * handles + stream). Every ComputationalMethod call forwards to the CPU pipeline; the
 * GPU is reached through the hooks the constructor installs on the owned XTB. That
 * shared wrapper logic lives ONCE in XtbGpuAdapter (xtb_gpu_adapter.h, Claude Generated
 * Sep 2026); only the ROCm-specific eigensolver hook, the resident backend and the
 * mixed-precision policy are in the .cpp.
 *
 * Usage:
 *   ./curcuma -sp mol.xyz -method gfn2 -gpu rocm   # explicit ROCm
 *   ./curcuma -sp mol.xyz -method gfn1 -gpu auto   # GPU if available
 *   ./curcuma -sp mol.xyz -method gfn2             # CPU (default)
 */

#pragma once

#ifdef USE_ROCM

#include "xtb_gpu_adapter.h"        // XtbGpuAdapter (shared wrapper logic)
#include "rocm/xtb_hip_context.h"   // XtbHipContext (pimpl — no HIP headers)

/**
 * @brief GPU-accelerated native GFN1/GFN2 (method "gfn1"/"gfn2" with -gpu rocm).
 *
 * A thin instantiation: the whole ComputationalMethod surface (forwarding, device
 * handshake + logging, CPU fallback, gpuActive()) comes from the shared adapter; the
 * constructor adds the rocSOLVER eigensolver hook and the device-resident SCF backend.
 * Claude Generated.
 */
class XtbHipComputationalMethod
    : public curcuma::xtb::gpu::XtbGpuAdapter<curcuma::xtb::gpu::XtbHipContext> {
public:
    explicit XtbHipComputationalMethod(curcuma::xtb::MethodType method,
                                       const json& config = json{});
    ~XtbHipComputationalMethod() override;
};

#endif // USE_ROCM
