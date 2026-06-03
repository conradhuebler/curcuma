/*
 * <Native xTB GPU Context — implementation>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0.
 *
 * Claude Generated (2026-06): cuSOLVER/cuBLAS handle + stream lifecycle. Kept
 * free of CurcumaLogger so nvcc never sees host-only logging headers; the host
 * wrapper does all user-facing logging based on ok()/deviceName().
 */

#ifdef USE_CUDA_XTB

#include "xtb_gpu_context.h"

#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cusolverDn.h>

// CudaBuffer<T> (RAII cudaMalloc/cudaFree) — the project's established device
// allocation helper; routing every allocation through it is the mitigation for
// the GFN-FF GPU heap-corruption class of bug.
#include "../../ff_methods/cuda/gfnff_soa.h"

namespace curcuma {
namespace xtb {
namespace gpu {

struct XtbGpuContext::Impl {
    cudaStream_t       stream   = nullptr;
    cublasHandle_t     cublas   = nullptr;
    cusolverDnHandle_t cusolver = nullptr;
    int                device   = -1;
    std::string        name;
    bool               ok       = false;

    // Device-resident GFN1 SCF state (Stage 2). Allocated once per geometry in
    // residentBegin and reused across SCF iterations: H0/S/L are uploaded once,
    // C/P are produced on the device and never leave until residentFinalize.
    int                resident_n = 0;
    CudaBuffer<double> dH0, dS, dL;      // geometry-constant (uploaded once)
    CudaBuffer<double> dC, dP, dCw;      // eigenvectors / density / scaled cols
    CudaBuffer<double> dEps, dVao, dOcc, dPop;  // length-n working vectors
    CudaBuffer<double> dWork;            // cuSOLVER dsyevd workspace
    CudaBuffer<int>    dInfo;            // cuSOLVER devInfo
    int                lwork = 0;
};

// ---- Device kernels (Stage 2 resident SCF) -------------------------------
// All matrices column-major n×n: element (i,j) at index i + j*n.

// F = H0 − ½·S·(v_ao(i) + v_ao(j))   (GFN1 isotropic Fock; H0, S symmetric).
__global__ void k_build_fock_iso(double* F, const double* H0, const double* S,
                                 const double* vao, int n)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i < n && j < n) {
        const size_t idx = static_cast<size_t>(i) + static_cast<size_t>(j) * n;
        F[idx] = H0[idx] - 0.5 * S[idx] * (vao[i] + vao[j]);
    }
}

// Scale the leading ncol columns of C by the occupation: Cw(:,k) = occ[k]·C(:,k).
__global__ void k_scale_cols(double* Cw, const double* C, const double* occ,
                             int n, int ncol)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int k = blockIdx.y * blockDim.y + threadIdx.y;
    if (i < n && k < ncol) {
        const size_t idx = static_cast<size_t>(i) + static_cast<size_t>(k) * n;
        Cw[idx] = C[idx] * occ[k];
    }
}

// Mulliken AO populations pop(μ) = Σ_ν P(μ,ν)·S(μ,ν)  (S symmetric → = Σ P_μν S_νμ).
__global__ void k_pop_ao(double* pop, const double* P, const double* S, int n)
{
    const int mu = blockIdx.x * blockDim.x + threadIdx.x;
    if (mu < n) {
        double s = 0.0;
        for (int nu = 0; nu < n; ++nu) {
            const size_t idx = static_cast<size_t>(mu) + static_cast<size_t>(nu) * n;
            s += P[idx] * S[idx];
        }
        pop[mu] = s;
    }
}

XtbGpuContext::XtbGpuContext()
    : m_impl(std::make_unique<Impl>())
{
    int count = 0;
    if (cudaGetDeviceCount(&count) != cudaSuccess || count == 0)
        return; // no device — caller falls back to CPU

    if (cudaGetDevice(&m_impl->device) != cudaSuccess)
        return;

    cudaDeviceProp prop{};
    if (cudaGetDeviceProperties(&prop, m_impl->device) == cudaSuccess)
        m_impl->name = prop.name;

    if (cudaStreamCreate(&m_impl->stream) != cudaSuccess)
        return;

    if (cublasCreate(&m_impl->cublas) != CUBLAS_STATUS_SUCCESS)
        return;
    cublasSetStream(m_impl->cublas, m_impl->stream);

    if (cusolverDnCreate(&m_impl->cusolver) != CUSOLVER_STATUS_SUCCESS)
        return;
    cusolverDnSetStream(m_impl->cusolver, m_impl->stream);

    m_impl->ok = true;
}

XtbGpuContext::~XtbGpuContext()
{
    if (!m_impl)
        return;
    if (m_impl->cusolver) cusolverDnDestroy(m_impl->cusolver);
    if (m_impl->cublas)   cublasDestroy(m_impl->cublas);
    if (m_impl->stream)   cudaStreamDestroy(m_impl->stream);
}

bool XtbGpuContext::ok() const { return m_impl && m_impl->ok; }

std::string XtbGpuContext::deviceName() const
{
    return m_impl ? m_impl->name : std::string();
}

int XtbGpuContext::deviceId() const { return m_impl ? m_impl->device : -1; }

bool XtbGpuContext::deviceAvailable()
{
    int count = 0;
    return cudaGetDeviceCount(&count) == cudaSuccess && count > 0;
}

bool XtbGpuContext::solveGeneralizedEigenF64(const double* F, const double* L, int n,
                                             double* C, double* eps)
{
    if (!ok() || n <= 0 || !F || !L || !C || !eps)
        return false;

    cudaStream_t stream = m_impl->stream;
    const size_t nn = static_cast<size_t>(n) * static_cast<size_t>(n);

    CudaBuffer<double> dA, dL, dW;
    try {
        dA.alloc(static_cast<int>(nn)); // Fock → Ã → C̃ → C (in place)
        dL.alloc(static_cast<int>(nn)); // lower Cholesky factor L (constant per geometry)
        dW.alloc(n);                    // eigenvalues
    } catch (...) {
        return false;
    }

    if (cudaMemcpyAsync(dA.ptr, F, sizeof(double) * nn, cudaMemcpyHostToDevice, stream) != cudaSuccess)
        return false;
    if (cudaMemcpyAsync(dL.ptr, L, sizeof(double) * nn, cudaMemcpyHostToDevice, stream) != cudaSuccess)
        return false;

    const double one = 1.0;

    // Step 1: dA <- F · L⁻ᵀ      (solve X·Lᵀ = F)
    if (cublasDtrsm(m_impl->cublas, CUBLAS_SIDE_RIGHT, CUBLAS_FILL_MODE_LOWER,
                    CUBLAS_OP_T, CUBLAS_DIAG_NON_UNIT, n, n, &one, dL.ptr, n, dA.ptr, n)
        != CUBLAS_STATUS_SUCCESS)
        return false;
    // Step 2: dA <- L⁻¹ · dA     (solve L·X = dA)  ⇒  Ã = L⁻¹·F·L⁻ᵀ
    if (cublasDtrsm(m_impl->cublas, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_LOWER,
                    CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, n, n, &one, dL.ptr, n, dA.ptr, n)
        != CUBLAS_STATUS_SUCCESS)
        return false;

    // Standard symmetric eigensolve of Ã (lower): eigenvectors overwrite dA
    // (column-major), ascending eigenvalues into dW.
    int lwork = 0;
    if (cusolverDnDsyevd_bufferSize(m_impl->cusolver, CUSOLVER_EIG_MODE_VECTOR,
                                    CUBLAS_FILL_MODE_LOWER, n, dA.ptr, n, dW.ptr, &lwork)
        != CUSOLVER_STATUS_SUCCESS)
        return false;

    CudaBuffer<double> dWork;
    CudaBuffer<int>    dInfo;
    try {
        dWork.alloc(lwork > 0 ? lwork : 1);
        dInfo.alloc(1);
    } catch (...) {
        return false;
    }

    if (cusolverDnDsyevd(m_impl->cusolver, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_LOWER,
                         n, dA.ptr, n, dW.ptr, dWork.ptr, lwork, dInfo.ptr)
        != CUSOLVER_STATUS_SUCCESS)
        return false;

    // Step 3: dA <- L⁻ᵀ · C̃      (solve Lᵀ·C = C̃)  ⇒  generalized eigenvectors C
    if (cublasDtrsm(m_impl->cublas, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_LOWER,
                    CUBLAS_OP_T, CUBLAS_DIAG_NON_UNIT, n, n, &one, dL.ptr, n, dA.ptr, n)
        != CUBLAS_STATUS_SUCCESS)
        return false;

    int info = 1;
    if (cudaMemcpyAsync(&info, dInfo.ptr, sizeof(int), cudaMemcpyDeviceToHost, stream) != cudaSuccess)
        return false;
    if (cudaMemcpyAsync(C, dA.ptr, sizeof(double) * nn, cudaMemcpyDeviceToHost, stream) != cudaSuccess)
        return false;
    if (cudaMemcpyAsync(eps, dW.ptr, sizeof(double) * n, cudaMemcpyDeviceToHost, stream) != cudaSuccess)
        return false;
    if (cudaStreamSynchronize(stream) != cudaSuccess)
        return false;

    return info == 0; // cusolver devInfo: 0 = success
}

/* ====================================================================== *
 *  Device-resident GFN1 SCF (Stage 2). H0/S/L stay on the GPU for the
 *  whole loop; only length-n vectors cross the bus per iteration.
 * ====================================================================== */

bool XtbGpuContext::residentBegin(const double* H0, const double* S,
                                  const double* L, int n)
{
    if (!ok() || n <= 0 || !H0 || !S || !L)
        return false;
    const size_t nn = static_cast<size_t>(n) * static_cast<size_t>(n);
    try {
        // Geometry-constant matrices, uploaded exactly once for this geometry.
        m_impl->dH0.upload(H0, static_cast<int>(nn), m_impl->stream);
        m_impl->dS.upload(S, static_cast<int>(nn), m_impl->stream);
        m_impl->dL.upload(L, static_cast<int>(nn), m_impl->stream);
        // Resident work buffers, reused across iterations.
        m_impl->dC.alloc(static_cast<int>(nn));
        m_impl->dP.alloc(static_cast<int>(nn));
        m_impl->dCw.alloc(static_cast<int>(nn));
        m_impl->dEps.alloc(n);
        m_impl->dVao.alloc(n);
        m_impl->dOcc.alloc(n);
        m_impl->dPop.alloc(n);
        // cuSOLVER dsyevd workspace (size is geometry-constant for fixed n).
        int lwork = 0;
        if (cusolverDnDsyevd_bufferSize(m_impl->cusolver, CUSOLVER_EIG_MODE_VECTOR,
                                        CUBLAS_FILL_MODE_LOWER, n, m_impl->dC.ptr, n,
                                        m_impl->dEps.ptr, &lwork)
            != CUSOLVER_STATUS_SUCCESS)
            return false;
        m_impl->lwork = lwork;
        m_impl->dWork.alloc(lwork > 0 ? lwork : 1);
        m_impl->dInfo.alloc(1);
    } catch (...) {
        return false;
    }
    m_impl->resident_n = n;
    return cudaStreamSynchronize(m_impl->stream) == cudaSuccess;
}

bool XtbGpuContext::residentSolve(const double* v_ao, int n, double* eps_out)
{
    if (!ok() || n <= 0 || n != m_impl->resident_n || !v_ao || !eps_out)
        return false;
    cudaStream_t stream = m_impl->stream;

    // Upload the AO potential (the only matrix-sized input is already resident).
    m_impl->dVao.upload(v_ao, n, stream);

    // F = H0 − ½·S·(v_ao⊕v_ao), built straight into the eigenvector buffer dC.
    const dim3 block(16, 16);
    const dim3 grid((n + block.x - 1) / block.x, (n + block.y - 1) / block.y);
    k_build_fock_iso<<<grid, block, 0, stream>>>(m_impl->dC.ptr, m_impl->dH0.ptr,
                                                 m_impl->dS.ptr, m_impl->dVao.ptr, n);
    if (cudaGetLastError() != cudaSuccess)
        return false;

    // Reduce to standard form Ã = L⁻¹·F·L⁻ᵀ via two triangular solves, then
    // dsyevd, then back-transform C = L⁻ᵀ·C̃ — identical to the Stage-1 path but
    // reusing the resident L (no per-iteration L upload). dC holds F → Ã → C̃ → C.
    const double one = 1.0;
    if (cublasDtrsm(m_impl->cublas, CUBLAS_SIDE_RIGHT, CUBLAS_FILL_MODE_LOWER,
                    CUBLAS_OP_T, CUBLAS_DIAG_NON_UNIT, n, n, &one,
                    m_impl->dL.ptr, n, m_impl->dC.ptr, n) != CUBLAS_STATUS_SUCCESS)
        return false;
    if (cublasDtrsm(m_impl->cublas, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_LOWER,
                    CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, n, n, &one,
                    m_impl->dL.ptr, n, m_impl->dC.ptr, n) != CUBLAS_STATUS_SUCCESS)
        return false;
    if (cusolverDnDsyevd(m_impl->cusolver, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_LOWER,
                         n, m_impl->dC.ptr, n, m_impl->dEps.ptr, m_impl->dWork.ptr,
                         m_impl->lwork, m_impl->dInfo.ptr) != CUSOLVER_STATUS_SUCCESS)
        return false;
    if (cublasDtrsm(m_impl->cublas, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_LOWER,
                    CUBLAS_OP_T, CUBLAS_DIAG_NON_UNIT, n, n, &one,
                    m_impl->dL.ptr, n, m_impl->dC.ptr, n) != CUBLAS_STATUS_SUCCESS)
        return false;

    int info = 1;
    if (cudaMemcpyAsync(&info, m_impl->dInfo.ptr, sizeof(int),
                        cudaMemcpyDeviceToHost, stream) != cudaSuccess)
        return false;
    m_impl->dEps.download(eps_out, n, stream);
    if (cudaStreamSynchronize(stream) != cudaSuccess)
        return false;
    return info == 0;
}

bool XtbGpuContext::residentDensity(const double* occ, int ncol, int n,
                                    double* pop_ao_out, double* band_out)
{
    if (!ok() || n <= 0 || n != m_impl->resident_n || !pop_ao_out || !band_out)
        return false;
    cudaStream_t stream = m_impl->stream;
    const double one = 1.0, zero = 0.0;

    if (ncol > 0) {
        m_impl->dOcc.upload(occ, ncol, stream);
        // Cw(:,k) = occ[k]·C(:,k)  for k < ncol.
        const dim3 block(16, 16);
        const dim3 grid((n + block.x - 1) / block.x, (ncol + block.y - 1) / block.y);
        k_scale_cols<<<grid, block, 0, stream>>>(m_impl->dCw.ptr, m_impl->dC.ptr,
                                                 m_impl->dOcc.ptr, n, ncol);
        if (cudaGetLastError() != cudaSuccess)
            return false;
        // P = Cw[:, :ncol] · C[:, :ncol]ᵀ  ⇒  P_ij = Σ_k occ_k C_ik C_jk.
        if (cublasDgemm(m_impl->cublas, CUBLAS_OP_N, CUBLAS_OP_T, n, n, ncol,
                        &one, m_impl->dCw.ptr, n, m_impl->dC.ptr, n,
                        &zero, m_impl->dP.ptr, n) != CUBLAS_STATUS_SUCCESS)
            return false;
    } else {
        // No occupied orbitals — zero density.
        if (cudaMemsetAsync(m_impl->dP.ptr, 0, sizeof(double) * static_cast<size_t>(n) * n,
                            stream) != cudaSuccess)
            return false;
    }

    // pop_ao(μ) = Σ_ν P(μ,ν)·S(μ,ν).
    const int b1 = 128;
    k_pop_ao<<<(n + b1 - 1) / b1, b1, 0, stream>>>(m_impl->dPop.ptr, m_impl->dP.ptr,
                                                   m_impl->dS.ptr, n);
    if (cudaGetLastError() != cudaSuccess)
        return false;

    // Band energy = Σ_μν P_μν·H0_μν (host-pointer dot is blocking; sync to be safe).
    const size_t nn = static_cast<size_t>(n) * static_cast<size_t>(n);
    if (cublasDdot(m_impl->cublas, static_cast<int>(nn), m_impl->dP.ptr, 1,
                   m_impl->dH0.ptr, 1, band_out) != CUBLAS_STATUS_SUCCESS)
        return false;

    m_impl->dPop.download(pop_ao_out, n, stream);
    return cudaStreamSynchronize(stream) == cudaSuccess;
}

bool XtbGpuContext::residentFinalize(double* P_colmajor, double* C_colmajor, int n)
{
    if (!ok() || n <= 0 || n != m_impl->resident_n || !P_colmajor || !C_colmajor)
        return false;
    cudaStream_t stream = m_impl->stream;
    const size_t nn = static_cast<size_t>(n) * static_cast<size_t>(n);
    m_impl->dP.download(P_colmajor, static_cast<int>(nn), stream);
    m_impl->dC.download(C_colmajor, static_cast<int>(nn), stream);
    return cudaStreamSynchronize(stream) == cudaSuccess;
}

} // namespace gpu
} // namespace xtb
} // namespace curcuma

#endif // USE_CUDA_XTB
