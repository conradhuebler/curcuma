/*
 * <Native xTB ROCm/HIP Context — implementation>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0.
 *
 * Claude Generated (2026-06): Stage-0 device handshake for the native GFN1/GFN2
 * ROCm path. Stage 0 calls ONLY the host-side HIP runtime API (device query + stream),
 * so it is compiled as plain C++ against <hip/hip_runtime_api.h> and linked against
 * libamdhip64 — no HIP-language compilation / amdclang / --offload-arch is needed yet.
 * The numerical kernels (hipified from cuda/xtb_gpu_context.cu) arrive in later stages;
 * the file is then split so the device kernels build with hipcc.
 */

#ifdef USE_ROCM_XTB

#include "xtb_hip_context.h"

#include <hip/hip_runtime_api.h>

#include <memory>
#include <string>

namespace curcuma {
namespace xtb {
namespace gpu {

struct XtbHipContext::Impl {
    bool        ok = false;
    int         device = -1;
    std::string name;
    hipStream_t stream = nullptr;
    // Stage 1+: hipblasHandle_t blas; rocblas_handle solver; (rocSOLVER reuses rocBLAS)
};

XtbHipContext::XtbHipContext()
    : m_impl(std::make_unique<Impl>())
{
    int count = 0;
    if (hipGetDeviceCount(&count) != hipSuccess || count <= 0)
        return;
    if (hipSetDevice(0) != hipSuccess)
        return;

    hipDeviceProp_t prop;
    if (hipGetDeviceProperties(&prop, 0) != hipSuccess)
        return;
    if (hipStreamCreate(&m_impl->stream) != hipSuccess)
        return;

    m_impl->device = 0;
    m_impl->name   = prop.name;
    m_impl->ok     = true;
}

XtbHipContext::~XtbHipContext()
{
    if (m_impl && m_impl->stream)
        hipStreamDestroy(m_impl->stream);
}

bool XtbHipContext::ok() const { return m_impl && m_impl->ok; }

std::string XtbHipContext::deviceName() const
{
    return m_impl ? m_impl->name : std::string();
}

int XtbHipContext::deviceId() const { return m_impl ? m_impl->device : -1; }

bool XtbHipContext::deviceAvailable()
{
    int count = 0;
    return hipGetDeviceCount(&count) == hipSuccess && count > 0;
}

} // namespace gpu
} // namespace xtb
} // namespace curcuma

#endif // USE_ROCM_XTB
