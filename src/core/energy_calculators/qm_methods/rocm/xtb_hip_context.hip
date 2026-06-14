/*
 * <Native xTB ROCm/HIP Context — implementation>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0.
 *
 * Claude Generated (2026-06): native GFN1/GFN2 ROCm device engine.
 *   Stage 0: device handshake (query + stream).
 *   Stage 1: solveGeneralized — the per-iteration F C = S C ε via rocSOLVER dsygvd.
 *   Stage 2: device-resident GFN1 SCF (residentBegin/Solve/Density/Finalize) — Fock
 *            build + populations via the __global__ kernels below, density via rocBLAS,
 *            eigensolve via rocSOLVER; H0/S stay resident.
 * Compiled by hipcc into a plain object (CMake add_custom_command), then linked into
 * curcuma_core with g++ — so the device kernels build, but the HIP-only --offload-arch /
 * --hip-link flags never reach the GNU link (which would drop libgomp). Links libamdhip64
 * + librocsolver + librocblas.
 */

#ifdef USE_ROCM_XTB

#include "xtb_hip_context.h"

#include <hip/hip_runtime.h>
#ifdef HAVE_ROCSOLVER
#include <rocblas/rocblas.h>
#include <rocsolver/rocsolver.h>
#endif

#include <memory>
#include <string>
#include <vector>

namespace curcuma {
namespace xtb {
namespace gpu {

#ifdef HAVE_ROCSOLVER
namespace {
// GFN1 isotropic Fock: F[i,j] = H0[i,j] − ½·S[i,j]·(v[i] + v[j]). Column-major n×n.
__global__ void k_fock(const double* H0, const double* S, const double* v, double* F, int n)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;  // row
    const int j = blockIdx.y * blockDim.y + threadIdx.y;  // col
    if (i >= n || j >= n) return;
    const int idx = i + j * n;
    F[idx] = H0[idx] - 0.5 * S[idx] * (v[i] + v[j]);
}
// Column scaling: B[i,j] = A[i,j]·d[j] (forms Cw = C·diag(occ)).
__global__ void k_scale_cols(const double* A, const double* d, double* B, int n)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i >= n || j >= n) return;
    const int idx = i + j * n;
    B[idx] = A[idx] * d[j];
}
// Mulliken populations + band-energy partials from the resident density P (symmetric):
//   out[i]   = Σ_j P[i,j]·S[i,j],   out[n+i] = Σ_j P[i,j]·H0[i,j].
__global__ void k_popband(const double* P, const double* S, const double* H0, double* out, int n)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    double pop = 0.0, bnd = 0.0;
    for (int j = 0; j < n; ++j) {
        const double p = P[i + j * n];
        pop += p * S[i + j * n];
        bnd += p * H0[i + j * n];
    }
    out[i] = pop; out[n + i] = bnd;
}
} // namespace
#endif // HAVE_ROCSOLVER

struct XtbHipContext::Impl {
    bool        ok = false;
    int         device = -1;
    std::string name;
    hipStream_t stream = nullptr;

#ifdef HAVE_ROCSOLVER
    rocblas_handle handle = nullptr;   // rocSOLVER reuses the rocBLAS handle
    bool    handle_ok = false;
    int     solver_n  = 0;             // cached device-buffer size
    double* dF   = nullptr;            // F → eigenvectors (n×n)
    double* dS   = nullptr;            // S (overwritten by its Cholesky factor)
    double* dEps = nullptr;            // eigenvalues (n)
    double* dE   = nullptr;            // rocSOLVER workspace (n)
    rocblas_int* dInfo = nullptr;      // convergence flag

    bool ensureSolver(int n) {
        if (!handle_ok) {
            if (rocblas_create_handle(&handle) != rocblas_status_success) return false;
            rocblas_set_stream(handle, stream);
            handle_ok = true;
        }
        if (n == solver_n && dF) return true;
        freeSolver();
        const size_t nn = static_cast<size_t>(n) * n;
        if (hipMalloc(&dF, nn * sizeof(double)) != hipSuccess) return false;
        if (hipMalloc(&dS, nn * sizeof(double)) != hipSuccess) return false;
        if (hipMalloc(&dEps, n * sizeof(double)) != hipSuccess) return false;
        if (hipMalloc(&dE, n * sizeof(double)) != hipSuccess) return false;
        if (hipMalloc(&dInfo, sizeof(rocblas_int)) != hipSuccess) return false;
        solver_n = n;
        return true;
    }
    void freeSolver() {
        for (void* p : { (void*)dF, (void*)dS, (void*)dEps, (void*)dE, (void*)dInfo })
            if (p) hipFree(p);
        dF = dS = dEps = dE = nullptr; dInfo = nullptr; solver_n = 0;
    }

    // ----- resident GFN1 SCF (Stage 2) ------------------------------------
    // Reuses dF (eigenvectors C), dS (Cholesky work copy), dEps/dE/dInfo from the
    // solver cache; adds the geometry-resident H0/S + density work buffers.
    int     res_n = 0;
    double* dH0 = nullptr;   // resident H0 (n×n)
    double* dS0 = nullptr;   // resident S  (n×n, untouched; dS is the work copy)
    double* dvao = nullptr;  // v_ao (n)
    double* dP  = nullptr;   // density (n×n)
    double* dCw = nullptr;   // C·diag(occ) (n×n)
    double* dOcc = nullptr;  // occupations (n)
    double* dPB = nullptr;   // pop + bandpart (2n)

    bool ensureResident(int n) {
        if (!ensureSolver(n)) return false;
        if (n == res_n && dH0) return true;
        freeResident();
        const size_t nn = static_cast<size_t>(n) * n;
        if (hipMalloc(&dH0,  nn * sizeof(double)) != hipSuccess) return false;
        if (hipMalloc(&dS0,  nn * sizeof(double)) != hipSuccess) return false;
        if (hipMalloc(&dvao, n  * sizeof(double)) != hipSuccess) return false;
        if (hipMalloc(&dP,   nn * sizeof(double)) != hipSuccess) return false;
        if (hipMalloc(&dCw,  nn * sizeof(double)) != hipSuccess) return false;
        if (hipMalloc(&dOcc, n  * sizeof(double)) != hipSuccess) return false;
        if (hipMalloc(&dPB,  2 * n * sizeof(double)) != hipSuccess) return false;
        res_n = n;
        return true;
    }
    void freeResident() {
        for (void* p : { (void*)dH0, (void*)dS0, (void*)dvao, (void*)dP, (void*)dCw, (void*)dOcc, (void*)dPB })
            if (p) hipFree(p);
        dH0 = dS0 = dvao = dP = dCw = dOcc = dPB = nullptr; res_n = 0;
    }

    bool resBegin(const double* H0, const double* S, int n) {
        if (!ok || n <= 0 || !ensureResident(n)) return false;
        const size_t nn = static_cast<size_t>(n) * n;
        if (hipMemcpy(dH0, H0, nn * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
        if (hipMemcpy(dS0, S,  nn * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
        return true;
    }
    bool resSolve(const double* v_ao, double* eps_out) {
        const int n = res_n; if (n <= 0) return false;
        const size_t nn = static_cast<size_t>(n) * n;
        if (hipMemcpy(dvao, v_ao, n * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
        dim3 blk(16, 16), grd((n + 15) / 16, (n + 15) / 16);
        hipLaunchKernelGGL(k_fock, grd, blk, 0, stream, dH0, dS0, dvao, dF, n);   // F → dF
        if (hipMemcpyAsync(dS, dS0, nn * sizeof(double), hipMemcpyDeviceToDevice, stream) != hipSuccess) return false;
        if (hipStreamSynchronize(stream) != hipSuccess) return false;
        // F C = S C ε on the device; dF ← eigenvectors (ascending), dS destroyed.
        if (rocsolver_dsygvd(handle, rocblas_eform_ax, rocblas_evect_original, rocblas_fill_lower,
                             n, dF, n, dS, n, dEps, dE, dInfo) != rocblas_status_success) return false;
        if (hipDeviceSynchronize() != hipSuccess) return false;
        rocblas_int info = 0;
        if (hipMemcpy(&info, dInfo, sizeof(rocblas_int), hipMemcpyDeviceToHost) != hipSuccess) return false;
        if (info != 0) return false;
        if (hipMemcpy(eps_out, dEps, n * sizeof(double), hipMemcpyDeviceToHost) != hipSuccess) return false;
        return true;   // dF holds C resident, ascending order
    }
    bool resDensity(const double* occ, int ncol, double* pop_ao, double* band) {
        const int n = res_n; if (n <= 0 || ncol < 0 || ncol > n) return false;
        std::vector<double> occ_full(n, 0.0);          // ascending order matches C columns
        for (int r = 0; r < ncol; ++r) occ_full[r] = occ[r];
        if (hipMemcpy(dOcc, occ_full.data(), n * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
        dim3 blk(16, 16), grd((n + 15) / 16, (n + 15) / 16);
        hipLaunchKernelGGL(k_scale_cols, grd, blk, 0, stream, dF, dOcc, dCw, n);   // Cw = C·diag(occ)
        const double one = 1.0, zero = 0.0;
        // P = Cw · Cᵀ
        if (rocblas_dgemm(handle, rocblas_operation_none, rocblas_operation_transpose,
                          n, n, n, &one, dCw, n, dF, n, &zero, dP, n) != rocblas_status_success) return false;
        dim3 blk1(64), grd1((n + 63) / 64);
        hipLaunchKernelGGL(k_popband, grd1, blk1, 0, stream, dP, dS0, dH0, dPB, n);
        if (hipStreamSynchronize(stream) != hipSuccess) return false;
        std::vector<double> pb(static_cast<size_t>(2) * n);
        if (hipMemcpy(pb.data(), dPB, 2 * n * sizeof(double), hipMemcpyDeviceToHost) != hipSuccess) return false;
        double b = 0.0;
        for (int i = 0; i < n; ++i) { pop_ao[i] = pb[i]; b += pb[n + i]; }
        *band = b;
        return true;
    }
    bool resFinalize(double* P, double* C) {
        const int n = res_n; if (n <= 0) return false;
        const size_t nn = static_cast<size_t>(n) * n;
        if (hipMemcpy(P, dP, nn * sizeof(double), hipMemcpyDeviceToHost) != hipSuccess) return false;
        if (hipMemcpy(C, dF, nn * sizeof(double), hipMemcpyDeviceToHost) != hipSuccess) return false;
        return true;
    }
#endif // HAVE_ROCSOLVER
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
    if (!m_impl) return;
#ifdef HAVE_ROCSOLVER
    m_impl->freeResident();
    m_impl->freeSolver();
    if (m_impl->handle_ok) rocblas_destroy_handle(m_impl->handle);
#endif
    if (m_impl->stream)
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

bool XtbHipContext::solveGeneralized(const double* F, const double* S, int n,
                                     double* eps_out, double* C_colmajor)
{
#ifndef HAVE_ROCSOLVER
    (void)F; (void)S; (void)n; (void)eps_out; (void)C_colmajor;
    return false;   // built without rocSOLVER → caller uses the CPU eigensolver
#else
    if (!m_impl || !m_impl->ok || n <= 0 || !F || !S || !eps_out || !C_colmajor)
        return false;
    Impl& d = *m_impl;
    if (!d.ensureSolver(n)) return false;
    const size_t nn = static_cast<size_t>(n) * n;

    // F C = S C ε  (itype = AX), eigenvectors of the ORIGINAL problem (CᵀSC = I).
    // dF is overwritten with the eigenvectors, dS with its Cholesky factor.
    if (hipMemcpy(d.dF, F, nn * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;
    if (hipMemcpy(d.dS, S, nn * sizeof(double), hipMemcpyHostToDevice) != hipSuccess) return false;

    const rocblas_status st = rocsolver_dsygvd(
        d.handle, rocblas_eform_ax, rocblas_evect_original, rocblas_fill_lower,
        n, d.dF, n, d.dS, n, d.dEps, d.dE, d.dInfo);
    if (st != rocblas_status_success) return false;
    if (hipDeviceSynchronize() != hipSuccess) return false;

    rocblas_int info = 0;
    if (hipMemcpy(&info, d.dInfo, sizeof(rocblas_int), hipMemcpyDeviceToHost) != hipSuccess) return false;
    if (info != 0) return false;   // S not positive-definite or no convergence

    if (hipMemcpy(C_colmajor, d.dF, nn * sizeof(double), hipMemcpyDeviceToHost) != hipSuccess) return false;
    if (hipMemcpy(eps_out, d.dEps, n * sizeof(double), hipMemcpyDeviceToHost) != hipSuccess) return false;
    return true;
#endif
}

bool XtbHipContext::residentBegin(const double* H0, const double* S, int n)
{
#ifndef HAVE_ROCSOLVER
    (void)H0; (void)S; (void)n; return false;
#else
    return m_impl && m_impl->ok && H0 && S && n > 0 && m_impl->resBegin(H0, S, n);
#endif
}
bool XtbHipContext::residentSolve(const double* v_ao, double* eps_out)
{
#ifndef HAVE_ROCSOLVER
    (void)v_ao; (void)eps_out; return false;
#else
    return m_impl && m_impl->ok && v_ao && eps_out && m_impl->resSolve(v_ao, eps_out);
#endif
}
bool XtbHipContext::residentDensity(const double* occ, int ncol, double* pop_ao, double* band)
{
#ifndef HAVE_ROCSOLVER
    (void)occ; (void)ncol; (void)pop_ao; (void)band; return false;
#else
    return m_impl && m_impl->ok && occ && pop_ao && band && m_impl->resDensity(occ, ncol, pop_ao, band);
#endif
}
bool XtbHipContext::residentFinalize(double* P_colmajor, double* C_colmajor)
{
#ifndef HAVE_ROCSOLVER
    (void)P_colmajor; (void)C_colmajor; return false;
#else
    return m_impl && m_impl->ok && P_colmajor && C_colmajor && m_impl->resFinalize(P_colmajor, C_colmajor);
#endif
}

} // namespace gpu
} // namespace xtb
} // namespace curcuma

#endif // USE_ROCM_XTB
