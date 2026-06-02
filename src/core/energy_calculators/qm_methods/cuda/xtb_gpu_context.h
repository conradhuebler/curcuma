/*
 * <Native xTB GPU Context — cuSOLVER / cuBLAS handle + stream owner>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0.
 *
 * Claude Generated (2026-06): Device-side context for the native GFN1/GFN2 GPU
 * path. Owns one CUDA stream, a cuBLAS handle, and a cuSOLVER dense handle, all
 * bound to the same stream. A pimpl keeps every CUDA header out of the host
 * translation units that merely hold an XtbGpuContext (the method wrapper, the
 * factory), so only the .cu is compiled by nvcc.
 *
 * Stage 0 (current): the context only proves the build/link and performs the
 * device handshake. The numerical GPU kernels (eigensolve, SCF, integrals,
 * gradient) are added on this object in later stages.
 */

#pragma once

#ifdef USE_CUDA_XTB

#include <memory>
#include <string>

namespace curcuma {
namespace xtb {
namespace gpu {

/**
 * @brief Opaque CUDA context (cuSOLVER + cuBLAS + stream) for native xTB.
 *
 * Non-copyable; one instance per XtbGpuComputationalMethod. Construction is
 * non-throwing: if no usable CUDA device is present, ok() returns false and the
 * caller falls back to the CPU path. Claude Generated.
 */
class XtbGpuContext {
public:
    XtbGpuContext();
    ~XtbGpuContext();

    XtbGpuContext(const XtbGpuContext&) = delete;
    XtbGpuContext& operator=(const XtbGpuContext&) = delete;

    /// True once stream + cuBLAS + cuSOLVER handles are created on a real device.
    bool ok() const;

    /// Selected device name (e.g. "NVIDIA GeForce RTX 5080"); empty if none.
    std::string deviceName() const;

    /// Selected CUDA device id, or -1 if none.
    int deviceId() const;

    /// True if at least one CUDA device is visible (static probe, no allocation).
    static bool deviceAvailable();

    /**
     * @brief Solve the generalized symmetric eigenproblem F C = S C ε on the GPU,
     * with S = L·Lᵀ, reusing the host-supplied lower Cholesky factor L.
     *
     * Reduction Ã = L⁻¹·F·L⁻ᵀ via two triangular solves (cublasDtrsm), standard
     * eigensolve cusolverDnDsyevd, back-transform C = L⁻ᵀ·C̃ via one more trsm —
     * the device analogue of the CPU "native" path (cuSOLVER has no standalone
     * sygst). All matrices are column-major, size n×n (F symmetric, L lower).
     * On success C holds the generalized eigenvectors (Cᵀ S C = I) and eps the
     * ascending eigenvalues. FP64. Returns false on any CUDA/cuSOLVER error.
     * Claude Generated (Stage 1).
     */
    bool solveGeneralizedEigenF64(const double* F_colmajor,
                                  const double* L_colmajor,
                                  int           n,
                                  double*       C_colmajor,
                                  double*       eps);

private:
    struct Impl;
    std::unique_ptr<Impl> m_impl;
};

} // namespace gpu
} // namespace xtb
} // namespace curcuma

#endif // USE_CUDA_XTB
