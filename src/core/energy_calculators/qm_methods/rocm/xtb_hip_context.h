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

    /// Solve the generalized symmetric-definite eigenproblem F C = S C ε on the GPU via
    /// rocSOLVER (rocsolver_dsygvd). F_colmajor / S_colmajor are n×n symmetric
    /// (column-major; symmetric so the storage order is irrelevant). On success eps_out
    /// (length n) holds the ASCENDING eigenvalues and C_colmajor the generalized
    /// eigenvectors (CᵀSC = I, column j ↔ eps_out[j]). Returns false on any HIP/rocSOLVER
    /// error, or always when built without rocSOLVER (HAVE_ROCSOLVER) — the caller then
    /// falls back to the CPU eigensolver. Device buffers + the rocBLAS handle are cached
    /// across calls. Host-callable (no device kernels). Claude Generated (Stage 1).
    bool solveGeneralized(const double* F_colmajor, const double* S_colmajor, int n,
                          double* eps_out, double* C_colmajor);

    // ---- Device-resident GFN1 SCF (Stage 2) ---------------------------------
    // Keeps H0/S (and the per-iteration density + MO coefficients) RESIDENT on the
    // device for one geometry; only length-n vectors cross the bus per iteration. The
    // Fock build + Mulliken populations are HIP __global__ kernels (hipcc), the density
    // is a rocBLAS GEMM, and the generalized eigensolve is rocSOLVER dsygvd (which
    // returns ascending eigenpairs directly — no host reduction / sort). All matrices
    // n×n column-major; H0/S/P symmetric. Any false → caller falls back to the CPU.
    // Claude Generated (Stage 2, GFN1 isotropic).

    /// Upload the geometry-constant H0, overlap S; allocate the resident work buffers.
    bool residentBegin(const double* H0_colmajor, const double* S_colmajor, int n);

    /// One SCF step: F = H0 − ½·S·(v_ao⊕v_ao) (kernel), solve F C = S C ε on the device
    /// (rocSOLVER), keep C resident. Writes the ascending eigenvalues to eps_out (n).
    bool residentSolve(const double* v_ao, double* eps_out);

    /// Density P = C·diag(occ)·Cᵀ over the leading ncol (ascending-eps) columns from the
    /// resident C (rocBLAS); return pop_ao(μ)=Σ_ν P_μν·S_μν and band = Σ_μν P_μν·H0_μν.
    bool residentDensity(const double* occ, int ncol, double* pop_ao, double* band);

    /// Download the converged density P and MO coefficients C (n×n column-major).
    bool residentFinalize(double* P_colmajor, double* C_colmajor);

private:
    struct Impl;
    std::unique_ptr<Impl> m_impl;
};

} // namespace gpu
} // namespace xtb
} // namespace curcuma

#endif // USE_ROCM_XTB
