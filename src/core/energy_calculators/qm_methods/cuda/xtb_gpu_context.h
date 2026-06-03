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

    /* ----- Device-resident GFN1 SCF (Stage 2) ---------------------------- *
     * These four calls keep H0, S, the Cholesky factor L and the per-iteration
     * density P and MO coefficients C RESIDENT on the device for one geometry,
     * so only length-n vectors cross the bus inside the SCF loop. All matrices
     * are column-major n×n; H0 and S are symmetric, so a row-major host buffer
     * may be passed unchanged. The driving SCF loop (mixing, energies,
     * occupation, convergence) stays on the host. Claude Generated. */

    /// Upload the geometry-constant H0, overlap S and lower Cholesky L (S=L·Lᵀ);
    /// allocate the resident F/C/P work buffers + the cuSOLVER workspace. Once
    /// per geometry, before the SCF loop. Returns false on any CUDA error.
    bool residentBegin(const double* H0_colmajor, const double* S_colmajor,
                       const double* L_colmajor, int n);

    /// One SCF step: build F = H0 − ½·S·(v_ao⊕v_ao) on the device, solve
    /// F C = S C ε with the cached L, write the ascending eigenvalues to eps_out
    /// (length n). The eigenvectors C stay resident. Returns false on error.
    bool residentSolve(const double* v_ao, int n, double* eps_out);

    /// Build the density P = C·diag(occ)·Cᵀ over the leading ncol columns, then
    /// return the Mulliken AO populations pop_ao_out(μ)=Σ_ν P_μν·S_μν (length n)
    /// and the band energy band_out = Σ_μν P_μν·H0_μν. P stays resident.
    bool residentDensity(const double* occ, int ncol, int n,
                         double* pop_ao_out, double* band_out);

    /// Download the converged density and MO coefficients to the host (column
    /// major; symmetric P may be read as row-major). Called once after the loop.
    bool residentFinalize(double* P_colmajor, double* C_colmajor, int n);

    /* ----- Device-resident GFN2 multipole (Stage 2b) --------------------- *
     * Layered on residentBegin: the geometry-constant dipole (3) and quadrupole
     * (6) AO integrals and the AO→atom map are uploaded once, then each iteration
     * adds the anisotropic Fock term and reads off the atomic moments. Matrices
     * column-major n×n; dp_int/qp_int are NOT symmetric (both triangles used).
     * Claude Generated. */

    /// Upload the geometry-constant multipole integrals (3 dipole + 6 quadrupole,
    /// each n×n column-major) and the AO→atom map (length n). After residentBegin.
    bool residentBeginMultipole(const double* dp_int3, const double* qp_int6,
                                const int* ao2at, int n, int nat);

    /// Like residentSolve, but the device Fock additionally gets the GFN2 multipole
    /// contribution −½·Σ(dp_int·v_dp + qp_int·v_qp) before the eigensolve.
    /// v_dp is 3×nat, v_qp is 6×nat (column-major). Eigenvalues → eps_out.
    bool residentSolveMultipole(const double* v_ao, const double* v_dp,
                                const double* v_qp, int n, double* eps_out);

    /// Atom-resolved multipole moments from the resident density:
    ///   dp_at(k,iat) = −Σ_{μ∈iat} Σ_ν P_νμ·dp_int[k]_νμ   (3×nat)
    ///   qp_at(k,iat) = −Σ_{μ∈iat} Σ_ν P_νμ·qp_int[k]_νμ   (6×nat)
    /// Call after residentDensity (reads the resident P). Column-major outputs.
    bool residentMultipoleMoments(double* dp_at3, double* qp_at6, int n, int nat);

private:
    /// Reduce the resident Fock in dC to standard form with the cached L, solve
    /// it (cusolverDnDsyevd) and back-transform → generalized eigenvectors in dC,
    /// ascending eigenvalues downloaded to eps_out. Shared by residentSolve and
    /// residentSolveMultipole (which differ only in how dC's Fock is built).
    bool eigensolveResidentFock(double* eps_out);

    struct Impl;
    std::unique_ptr<Impl> m_impl;
};

} // namespace gpu
} // namespace xtb
} // namespace curcuma

#endif // USE_CUDA_XTB
