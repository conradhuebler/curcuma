/*
 * <Native xTB ROCm/HIP Context — rocSOLVER / hipBLAS handle + stream owner>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0.
 *
 * Claude Generated (2026-06): Device-side context for the native GFN1/GFN2 ROCm
 * path — the AMD/HIP sibling of curcuma::xtb::gpu::XtbGpuContext (CUDA). Owns one
 * HIP stream and (in later stages) a hipBLAS handle + a rocSOLVER handle, all bound
 * to the same stream. A pimpl keeps every HIP header out of the host translation
 * units that merely hold an XtbHipContext (the method wrapper, the factory), so only
 * the .hip is compiled by hipcc.
 *
 * Stage 0 (current): the context only proves the build/link and performs the device
 * handshake (probe + stream). The numerical kernels (eigensolve, SCF, integrals,
 * gradient) are hipified from cuda/xtb_gpu_context.cu on this object in later stages,
 * exposing the SAME public method signatures so the GpuScfBackend adapter is shared.
 */

#pragma once

#ifdef USE_ROCM_XTB

#include <memory>
#include <string>

namespace curcuma {
namespace xtb {
namespace gpu {

/**
 * @brief Opaque HIP context (hipBLAS + rocSOLVER + stream) for native xTB.
 *
 * Non-copyable; one instance per XtbHipComputationalMethod. Construction is
 * non-throwing: if no usable HIP device is present, ok() returns false and the
 * caller falls back to the CPU path. Claude Generated.
 */
class XtbHipContext {
public:
    XtbHipContext();
    ~XtbHipContext();

    XtbHipContext(const XtbHipContext&) = delete;
    XtbHipContext& operator=(const XtbHipContext&) = delete;

    /// True once the stream (and, later, hipBLAS/rocSOLVER handles) are created on a
    /// real device.
    bool ok() const;

    /// Selected device name (e.g. "AMD Radeon RX 7900 XTX"); empty if none.
    std::string deviceName() const;

    /// Selected HIP device id, or -1 if none.
    int deviceId() const;

    /// True if at least one HIP device is visible (static probe, no allocation).
    static bool deviceAvailable();

    // ---- Stage 1+ (hipified from XtbGpuContext) -----------------------------
    // The resident SCF / integral / gradient API mirrors XtbGpuContext exactly
    // (residentBegin/Solve/Density/Finalize, beginBasis/computeIntegrals,
    // computeGradient, the Stage-5/6 hooks). Added with the kernel port so the
    // shared GpuScfBackend adapter forwards to either context unchanged.

private:
    struct Impl;
    std::unique_ptr<Impl> m_impl;
};

} // namespace gpu
} // namespace xtb
} // namespace curcuma

#endif // USE_ROCM_XTB
