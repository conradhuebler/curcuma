/*
 * <Native xTB Vulkan Context — compute engine over a VkContext>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0.
 *
 * Claude Generated (2026-06): Device-side engine for the native GFN1/GFN2 Vulkan path
 * — the Vulkan sibling of XtbGpuContext (CUDA) / XtbHipContext (ROCm). Owns a
 * curcuma::vk::VkContext (instance/device/compute-queue/command-pool) and, in later
 * stages, the compute pipelines + device buffers for the resident SCF. A pimpl keeps
 * <vulkan/vulkan.h> out of the host translation units that merely hold the context.
 *
 * Stage 0 (current): device handshake only. Stage 1 adds the GEMM/TRSM + cyclic
 * Jacobi symmetric eigensolver (the generalized problem F C = S C ε reduces via the
 * Cholesky factor L exactly like the CUDA path); later stages add the integral build,
 * resident SCF and nuclear gradient — exposing the SAME public API as XtbGpuContext so
 * the shared GpuScfBackend adapter forwards unchanged.
 */

#pragma once

#ifdef USE_VULKAN_XTB

#include <memory>
#include <string>

namespace curcuma {
namespace xtb {
namespace gpu {

class XtbVulkanContext {
public:
    XtbVulkanContext();
    ~XtbVulkanContext();

    XtbVulkanContext(const XtbVulkanContext&) = delete;
    XtbVulkanContext& operator=(const XtbVulkanContext&) = delete;

    /// True once the Vulkan device + compute queue are live (and shaderFloat64 present).
    bool ok() const;

    /// Selected device name; empty if none.
    std::string deviceName() const;

    /// Index of the selected physical device, or -1.
    int deviceId() const;

    /// True if a suitable compute+FP64 Vulkan device is visible (static probe).
    static bool deviceAvailable();

    /// Standard symmetric eigensolve on the GPU (FP64 two-sided cyclic Jacobi).
    /// A_colmajor is an n*n symmetric matrix (column-major; symmetric so the storage
    /// order is irrelevant). On success eps_out (length n) holds the ASCENDING
    /// eigenvalues and V_colmajor (n*n, column-major) the eigenvectors — column j is
    /// the eigenvector for eps_out[j], with VᵀV = I. Returns false on any Vulkan
    /// error so the caller can fall back to the CPU eigensolver. The generalized
    /// problem F C = S C ε is reduced to this standard form on the host (Cholesky L).
    /// Claude Generated (Stage 1, validated on AMD 890M / RADV vs Eigen).
    bool solveSymmetric(const double* A_colmajor, int n, double* eps_out, double* V_colmajor);

    // ---- Device-resident GFN1 SCF (Stage 2) ---------------------------------
    // Keeps H0/S (and the per-iteration density + MO coefficients) RESIDENT on the
    // device for one geometry, so only length-n vectors cross the bus per iteration.
    // The generalized problem F C = S C ε is reduced via the Löwdin transform
    // X = S⁻¹ᐟ² (built once on the device by residentBegin), not the host Cholesky L
    // — XᵀFX is two GEMMs, avoiding a GPU triangular solve. All matrices n×n
    // column-major; H0/S/P symmetric. Any false → caller falls back to the CPU.
    // Claude Generated (Stage 2a, GFN1 isotropic).

    /// Upload the geometry-constant H0, overlap S; build the resident X = S⁻¹ᐟ²;
    /// allocate the resident F/Ã/C/P work buffers. Once per geometry, before the loop.
    bool residentBegin(const double* H0_colmajor, const double* S_colmajor, int n);

    /// One SCF step: F = H0 − ½·S·(v_ao⊕v_ao), Ã = X·F·X, eigensolve, C = X·C̃
    /// (resident). Writes the ascending eigenvalues to eps_out (length n).
    bool residentSolve(const double* v_ao, double* eps_out);

    /// Build the density P = C·diag(occ)·Cᵀ over the leading ncol (ascending-eps)
    /// columns from the resident C; return Mulliken AO populations pop_ao(μ) =
    /// Σ_ν P_μν·S_μν (length n) and band = Σ_μν P_μν·H0_μν. P stays resident.
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

#endif // USE_VULKAN_XTB
