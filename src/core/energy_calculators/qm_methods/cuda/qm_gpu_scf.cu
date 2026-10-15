/*
 * < CUDA device-resident RHF SCF state of the native QM engine (Q-G4) >
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0
 */

// Claude Generated (Oct 2026, Q-G4 of docs/QM_GPU_ROADMAP.md).
//
// Device side of the resident closed-shell RHF SCF loop. The loop itself (iteration control,
// convergence test, DIIS/ADIIS coefficient solve, incremental-Fock cadence) stays in
// QMEngine::runSCF (qm_scf.cpp), so host and device loop share one driver; this class owns
// the n x n matrices on the device and performs the steps of one iteration on them:
//
//   Fock      F = H + J - 1/2 K              J, K from the device J/K builder (Q-G1..Q-G3),
//                                            full or incremental on P - P_ref
//   DIIS      e = X^T (F P S - S P F) X     commutator in the orthonormal basis (Pulay 1982);
//                                            Gram matrix B_ij = <e_i, e_j> and the ADIIS model
//                                            (Hu, Yang, J. Chem. Phys. 132, 054109 (2010)) are
//                                            reduced on the device, solved on the host (m <= 8)
//   eigen     F' = X^T F X,  F' C' = C' eps  cusolverDnDsyevd (as the xTB resident loop)
//             C = X C'
//   density   P_new = 2 C_occ C_occ^T       cublasDsyrk + mirror (exactly symmetric)
//   scalars   ||P_new - P||_F, E = 1/2 Tr(P_new (H + F)), max|e|
//
// cuBLAS works column-major. H, S, P, F, J, K are symmetric, so row- and column-major storage
// coincide; X is uploaded transposed (column-major) and C is returned row-major. All
// reductions use a fixed two-stage order (deterministic for a given n). FP64 throughout.

#include "../qm_gpu_scf_backend.h"

#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cusolverDn.h>

#include <chrono>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

namespace {

constexpr int kRedBlocks = 120;   // stage-1 blocks per reduction
constexpr int kRedThreads = 256;

/// F = (H + J) - 0.5 K and HF = H + F, element-wise; explicit rounding (no FMA contraction),
/// the expression of the host loop (m_H + J_acc - 0.5 * K_acc).
__global__ void fockKernel(const double* __restrict__ H, const double* __restrict__ J, const double* __restrict__ K,
                           double* __restrict__ F, double* __restrict__ HF, size_t nn)
{
    const size_t k = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= nn) return;
    const double f = __dsub_rn(__dadd_rn(H[k], J[k]), __dmul_rn(0.5, K[k]));
    F[k] = f;
    HF[k] = __dadd_rn(H[k], f);
}

/// C = A - B
__global__ void subKernel(const double* __restrict__ A, const double* __restrict__ B, double* __restrict__ C, size_t nn)
{
    const size_t k = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (k < nn) C[k] = __dsub_rn(A[k], B[k]);
}

/// A += B
__global__ void addKernel(double* __restrict__ A, const double* __restrict__ B, size_t nn)
{
    const size_t k = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (k < nn) A[k] = __dadd_rn(A[k], B[k]);
}

/// D = T - T^T (column-major n x n)
__global__ void antisymKernel(const double* __restrict__ T, double* __restrict__ D, int n)
{
    const size_t k = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= (size_t)n * n) return;
    const size_t i = k % n, j = k / n;
    D[k] = __dsub_rn(T[k], T[i * n + j]);
}

/// Copy the lower triangle (column-major, as written by syrk) to the upper one.
__global__ void mirrorLowerKernel(double* __restrict__ A, int n)
{
    const size_t k = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= (size_t)n * n) return;
    const size_t i = k % n, j = k / n;
    if (i < j) A[k] = A[j + i * n];
}

/// out = sum_i c_i hist[slot_i], summed in the host order (0 + c_0 F_0 + c_1 F_1 + ...).
__global__ void combineKernel(const double* __restrict__ hist, const int* __restrict__ slots,
                              const double* __restrict__ c, int m, double* __restrict__ out, size_t nn)
{
    const size_t k = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= nn) return;
    double acc = 0.0;
    for (int i = 0; i < m; ++i) acc = __dadd_rn(acc, __dmul_rn(c[i], hist[(size_t)slots[i] * nn + k]));
    out[k] = acc;
}

/// P = d P + (1 - d) Pn
__global__ void dampKernel(double* __restrict__ P, const double* __restrict__ Pn, double d, size_t nn)
{
    const size_t k = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (k < nn) P[k] = __dadd_rn(__dmul_rn(d, P[k]), __dmul_rn(1.0 - d, Pn[k]));
}

/// One reduction term: sum_k (A - Aref)_k (B - Bref)_k; a null ref counts as zero.
struct PairSpec {
    const double* A;
    const double* Aref;
    const double* B;
    const double* Bref;
};

__device__ inline void blockSum(double* sh, double v)
{
    sh[threadIdx.x] = v;
    __syncthreads();
    for (int o = blockDim.x / 2; o > 0; o >>= 1) {
        if (threadIdx.x < o) sh[threadIdx.x] += sh[threadIdx.x + o];
        __syncthreads();
    }
}

/// Stage 1: partial[p * gridDim.x + b] for pair p = blockIdx.y.
__global__ void pairDotStage1(const PairSpec* __restrict__ specs, size_t nn, double* __restrict__ partial)
{
    __shared__ double sh[kRedThreads];
    const PairSpec s = specs[blockIdx.y];
    double acc = 0.0;
    for (size_t k = (size_t)blockIdx.x * blockDim.x + threadIdx.x; k < nn; k += (size_t)gridDim.x * blockDim.x) {
        const double a = s.Aref ? __dsub_rn(s.A[k], s.Aref[k]) : s.A[k];
        const double b = s.Bref ? __dsub_rn(s.B[k], s.Bref[k]) : s.B[k];
        acc = __dadd_rn(acc, __dmul_rn(a, b));
    }
    blockSum(sh, acc);
    if (threadIdx.x == 0) partial[(size_t)blockIdx.y * gridDim.x + blockIdx.x] = sh[0];
}

/// Stage 2: one block per pair sums its nb partials (fixed order).
__global__ void pairDotStage2(const double* __restrict__ partial, int nb, double* __restrict__ out)
{
    __shared__ double sh[kRedThreads];
    double acc = 0.0;
    for (int b = threadIdx.x; b < nb; b += blockDim.x) acc += partial[(size_t)blockIdx.x * nb + b];
    blockSum(sh, acc);
    if (threadIdx.x == 0) out[blockIdx.x] = sh[0];
}

/// max |A| (two stages; max is order-independent).
__global__ void maxAbsStage1(const double* __restrict__ A, size_t nn, double* __restrict__ partial)
{
    __shared__ double sh[kRedThreads];
    double m = 0.0;
    for (size_t k = (size_t)blockIdx.x * blockDim.x + threadIdx.x; k < nn; k += (size_t)gridDim.x * blockDim.x)
        m = fmax(m, fabs(A[k]));
    sh[threadIdx.x] = m;
    __syncthreads();
    for (int o = blockDim.x / 2; o > 0; o >>= 1) {
        if (threadIdx.x < o) sh[threadIdx.x] = fmax(sh[threadIdx.x], sh[threadIdx.x + o]);
        __syncthreads();
    }
    if (threadIdx.x == 0) partial[blockIdx.x] = sh[0];
}

__global__ void maxStage2(const double* __restrict__ partial, int nb, double* __restrict__ out)
{
    __shared__ double sh[kRedThreads];
    double m = 0.0;
    for (int b = threadIdx.x; b < nb; b += blockDim.x) m = fmax(m, partial[b]);
    sh[threadIdx.x] = m;
    __syncthreads();
    for (int o = blockDim.x / 2; o > 0; o >>= 1) {
        if (threadIdx.x < o) sh[threadIdx.x] = fmax(sh[threadIdx.x], sh[threadIdx.x + o]);
        __syncthreads();
    }
    if (threadIdx.x == 0) out[0] = sh[0];
}

/// A = identity (column-major n x n); used once to warm up cuSOLVER/cuBLAS.
__global__ void identityKernel(double* __restrict__ A, int n)
{
    const size_t k = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (k < (size_t)n * n) A[k] = (k % n == k / n) ? 1.0 : 0.0;
}

inline unsigned grid1(size_t nn) { return (unsigned)((nn + 255) / 256); }

class QmGpuScfCuda final : public qmint::QMGpuSCF {
public:
    ~QmGpuScfCuda() override
    {
        release();
        if (m_cusolver) cusolverDnDestroy(m_cusolver);
        if (m_cublas) cublasDestroy(m_cublas);
        for (auto& e : m_ev)
            if (e) cudaEventDestroy(e);
    }

    bool init(int n, int device, std::string& why)
    {
        m_n = n;
        if (cublasCreate(&m_cublas) != CUBLAS_STATUS_SUCCESS) { why = "cublasCreate failed"; return false; }
        if (cusolverDnCreate(&m_cusolver) != CUSOLVER_STATUS_SUCCESS) { why = "cusolverDnCreate failed"; return false; }
        for (auto& e : m_ev)
            if (cudaEventCreate(&e) != cudaSuccess) { why = "cudaEventCreate failed"; return false; }
        const size_t nn = (size_t)n * n;
        bool ok = alloc(&m_H, nn) && alloc(&m_S, nn) && alloc(&m_X, nn) && alloc(&m_P, nn) && alloc(&m_Pn, nn)
            && alloc(&m_Pref, nn) && alloc(&m_dP, nn) && alloc(&m_J, nn) && alloc(&m_K, nn) && alloc(&m_dJ, nn)
            && alloc(&m_dK, nn) && alloc(&m_F, nn) && alloc(&m_HF, nn) && alloc(&m_Fuse, nn) && alloc(&m_T1, nn)
            && alloc(&m_T2, nn) && alloc(&m_C, nn) && alloc(&m_eps, (size_t)n) && alloc(&m_partial, (size_t)kRedBlocks * 1024)
            && alloc(&m_out, 1024) && alloc(&m_coef, 64) && alloc(&m_info, 1) && alloc(&m_slots, 64)
            && alloc(&m_specs, 1024);
        if (!ok) { why = "device allocation for the resident SCF failed"; return false; }
        int lwork = 0;
        if (cusolverDnDsyevd_bufferSize(m_cusolver, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_LOWER, n, m_T1, n,
                                        m_eps, &lwork) != CUSOLVER_STATUS_SUCCESS) {
            why = "cusolverDnDsyevd_bufferSize failed";
            return false;
        }
        if (!alloc(&m_work, (size_t)(lwork > 0 ? lwork : 1))) { why = "syevd workspace allocation failed"; return false; }
        m_lwork = lwork;
        // Warm-up: the first cuSOLVER/cuBLAS calls load their kernels lazily (measured ~65 ms for
        // n = 114 on an RTX 5090); do it here once, so it is accounted as setup, not as SCF time.
        const size_t nn0 = (size_t)n * n;
        identityKernel<<<grid1(nn0), 256>>>(m_T1, n);
        identityKernel<<<grid1(nn0), 256>>>(m_X, n);
        if (!gemm(CUBLAS_OP_T, CUBLAS_OP_N, m_X, m_T1, m_T2) || !gemm(CUBLAS_OP_N, CUBLAS_OP_N, m_T2, m_X, m_C)
            || cusolverDnDsyevd(m_cusolver, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_LOWER, n, m_T1, n, m_eps, m_work,
                                m_lwork, m_info) != CUSOLVER_STATUS_SUCCESS) {
            why = "cuSOLVER/cuBLAS warm-up failed";
            return false;
        }
        const double two = 2.0, zero = 0.0;
        cublasDsyrk(m_cublas, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_N, n, n, &two, m_C, n, &zero, m_Pn, n);
        if (cudaDeviceSynchronize() != cudaSuccess) { why = "warm-up synchronisation failed"; return false; }
        cudaDeviceProp prop{};
        cudaGetDeviceProperties(&prop, device >= 0 ? device : 0);
        m_name = std::string("cuda:") + prop.name;
        return cudaGetLastError() == cudaSuccess;
    }

    std::string name() const override { return m_name; }
    int n() const override { return m_n; }

    bool setup(const double* H, const double* S, const double* X, int n_occ, int nslots) override
    {
        const int n = m_n;
        const size_t nn = (size_t)n * n;
        if (n_occ < 0 || n_occ > n || nslots < 1 || nslots > 32) return false;
        m_nocc = n_occ;
        if (nslots != m_nslots) {
            cudaFree(m_histF); cudaFree(m_histP); cudaFree(m_histE);
            m_histF = m_histP = m_histE = nullptr;
            if (cudaMalloc(&m_histF, nslots * nn * sizeof(double)) != cudaSuccess
                || cudaMalloc(&m_histP, nslots * nn * sizeof(double)) != cudaSuccess
                || cudaMalloc(&m_histE, nslots * nn * sizeof(double)) != cudaSuccess) {
                m_nslots = 0;
                return false;
            }
            m_nslots = nslots;
        }
        std::vector<double> Xc(nn);  // row-major host X -> column-major device X
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) Xc[(size_t)j * n + i] = X[(size_t)i * n + j];
        m_t_fock = m_t_diis = m_t_eig = m_t_dens = 0.0;
        m_iters = 0;
        return cudaMemcpy(m_H, H, nn * sizeof(double), cudaMemcpyHostToDevice) == cudaSuccess
            && cudaMemcpy(m_S, S, nn * sizeof(double), cudaMemcpyHostToDevice) == cudaSuccess
            && cudaMemcpy(m_X, Xc.data(), nn * sizeof(double), cudaMemcpyHostToDevice) == cudaSuccess;
    }

    bool setDensity(const double* P) override
    {
        return cudaMemcpy(m_P, P, (size_t)m_n * m_n * sizeof(double), cudaMemcpyHostToDevice) == cudaSuccess;
    }

    bool buildFock(qmint::QMGpuJK* jk, bool full) override
    {
        if (!jk || jk->n() != m_n) return false;
        const size_t nn = (size_t)m_n * m_n;
        const auto t0 = std::chrono::steady_clock::now();
        if (full) {
            if (!jk->buildDevice(m_P, m_J, m_K)) return false;
        } else {
            subKernel<<<grid1(nn), 256>>>(m_P, m_Pref, m_dP, nn);
            if (cudaGetLastError() != cudaSuccess || !jk->buildDevice(m_dP, m_dJ, m_dK)) return false;
            addKernel<<<grid1(nn), 256>>>(m_J, m_dJ, nn);
            addKernel<<<grid1(nn), 256>>>(m_K, m_dK, nn);
        }
        if (cudaMemcpyAsync(m_Pref, m_P, nn * sizeof(double), cudaMemcpyDeviceToDevice, 0) != cudaSuccess) return false;
        fockKernel<<<grid1(nn), 256>>>(m_H, m_J, m_K, m_F, m_HF, nn);
        // J/K builders are synchronous; the few element-wise kernels after them are included
        // in the next synchronising call's time instead of forcing a sync here.
        m_t_fock += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
        return cudaGetLastError() == cudaSuccess;
    }

    bool diisPush(int slot, double& max_err) override
    {
        if (slot < 0 || slot >= m_nslots) return false;
        const int n = m_n;
        const size_t nn = (size_t)n * n;
        cudaEventRecord(m_ev[0]);
        double* Fs = m_histF + slot * nn;
        double* Ps = m_histP + slot * nn;
        double* Es = m_histE + slot * nn;
        if (cudaMemcpyAsync(Fs, m_F, nn * sizeof(double), cudaMemcpyDeviceToDevice, 0) != cudaSuccess
            || cudaMemcpyAsync(Ps, m_P, nn * sizeof(double), cudaMemcpyDeviceToDevice, 0) != cudaSuccess)
            return false;
        // (F P) S, as the host's F * P * S
        if (!gemm(CUBLAS_OP_N, CUBLAS_OP_N, m_F, m_P, m_T1) || !gemm(CUBLAS_OP_N, CUBLAS_OP_N, m_T1, m_S, m_T2))
            return false;
        antisymKernel<<<grid1(nn), 256>>>(m_T2, m_T1, n);  // T1 = FPS - (FPS)^T
        // (X^T D) X
        if (!gemm(CUBLAS_OP_T, CUBLAS_OP_N, m_X, m_T1, m_T2) || !gemm(CUBLAS_OP_N, CUBLAS_OP_N, m_T2, m_X, Es))
            return false;
        maxAbsStage1<<<kRedBlocks, kRedThreads>>>(Es, nn, m_partial);
        maxStage2<<<1, kRedThreads>>>(m_partial, kRedBlocks, m_out);
        cudaEventRecord(m_ev[1]);
        if (cudaMemcpy(&max_err, m_out, sizeof(double), cudaMemcpyDeviceToHost) != cudaSuccess) return false;
        accumulate(m_t_diis, 0, 1);
        return cudaGetLastError() == cudaSuccess;
    }

    bool diisErrorGram(const int* slots, int m, double* B) override
    {
        const size_t nn = (size_t)m_n * m_n;
        std::vector<PairSpec> specs;
        for (int i = 0; i < m; ++i)
            for (int j = 0; j < m; ++j)
                specs.push_back({ m_histE + slots[i] * nn, nullptr, m_histE + slots[j] * nn, nullptr });
        std::vector<double> out;
        if (!reducePairs(specs, out, m_t_diis)) return false;
        for (int k = 0; k < m * m; ++k) B[k] = out[k];
        return true;
    }

    bool adiisData(const int* slots, int m, double* a, double* M) override
    {
        const size_t nn = (size_t)m_n * m_n;
        const double* Pn = m_histP + slots[m - 1] * nn;
        const double* Fn = m_histF + slots[m - 1] * nn;
        std::vector<PairSpec> specs;
        for (int i = 0; i < m; ++i) specs.push_back({ m_histP + slots[i] * nn, Pn, Fn, nullptr });
        for (int i = 0; i < m; ++i)
            for (int j = 0; j < m; ++j) specs.push_back({ m_histP + slots[i] * nn, Pn, m_histF + slots[j] * nn, Fn });
        std::vector<double> out;
        if (!reducePairs(specs, out, m_t_diis)) return false;
        for (int i = 0; i < m; ++i) a[i] = out[i];
        for (int k = 0; k < m * m; ++k) M[k] = out[m + k];
        return true;
    }

    bool selectFock(const int* slots, const double* c, int m) override
    {
        const size_t nn = (size_t)m_n * m_n;
        if (m <= 0) return cudaMemcpyAsync(m_Fuse, m_F, nn * sizeof(double), cudaMemcpyDeviceToDevice, 0) == cudaSuccess;
        if (m > 64) return false;
        cudaEventRecord(m_ev[0]);
        if (cudaMemcpy(m_slots, slots, m * sizeof(int), cudaMemcpyHostToDevice) != cudaSuccess
            || cudaMemcpy(m_coef, c, m * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess)
            return false;
        combineKernel<<<grid1(nn), 256>>>(m_histF, m_slots, m_coef, m, m_Fuse, nn);
        cudaEventRecord(m_ev[1]);
        m_pending_diis = true;
        return cudaGetLastError() == cudaSuccess;
    }

    bool solveAndDensity(double& delta_P, double& e_elec) override
    {
        const int n = m_n;
        const size_t nn = (size_t)n * n;
        cudaEventRecord(m_ev[2]);
        // F' = (X^T F) X  -> T1 (overwritten by the eigenvectors C')
        if (!gemm(CUBLAS_OP_T, CUBLAS_OP_N, m_X, m_Fuse, m_T2) || !gemm(CUBLAS_OP_N, CUBLAS_OP_N, m_T2, m_X, m_T1))
            return false;
        cudaEventRecord(m_ev[3]);
        if (cusolverDnDsyevd(m_cusolver, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_LOWER, n, m_T1, n, m_eps, m_work,
                             m_lwork, m_info) != CUSOLVER_STATUS_SUCCESS)
            return false;
        cudaEventRecord(m_ev[4]);
        // C = X C'
        if (!gemm(CUBLAS_OP_N, CUBLAS_OP_N, m_X, m_T1, m_C)) return false;
        // P_new = 2 C_occ C_occ^T (lower) + mirror
        cudaMemsetAsync(m_Pn, 0, nn * sizeof(double), 0);
        if (m_nocc > 0) {
            const double two = 2.0, zero = 0.0;
            if (cublasDsyrk(m_cublas, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_N, n, m_nocc, &two, m_C, n, &zero, m_Pn, n)
                != CUBLAS_STATUS_SUCCESS)
                return false;
        }
        mirrorLowerKernel<<<grid1(nn), 256>>>(m_Pn, n);
        cudaEventRecord(m_ev[5]);
        // ||P_new - P||^2 and Tr(P_new (H + F))
        std::vector<PairSpec> specs = { { m_Pn, m_P, m_Pn, m_P }, { m_Pn, nullptr, m_HF, nullptr } };
        std::vector<double> out;
        double dummy = 0.0;
        if (!reducePairs(specs, out, dummy, /*info=*/true)) return false;
        if (m_last_info != 0) return false;  // syevd did not converge / illegal value
        delta_P = std::sqrt(out[0]);
        e_elec = 0.5 * out[1];
        // all events of this iteration are complete now
        if (m_pending_diis) { accumulate(m_t_diis, 0, 1); m_pending_diis = false; }
        accumulate(m_t_eig, 2, 4);
        accumulate(m_t_dens, 4, 5);
        ++m_iters;
        return true;
    }

    bool acceptDensity(double damping) override
    {
        const size_t nn = (size_t)m_n * m_n;
        dampKernel<<<grid1(nn), 256>>>(m_P, m_Pn, damping, nn);
        return cudaGetLastError() == cudaSuccess;
    }

    bool downloadDensity(bool newest, double* P) const override
    {
        return cudaMemcpy(P, newest ? m_Pn : m_P, (size_t)m_n * m_n * sizeof(double), cudaMemcpyDeviceToHost)
            == cudaSuccess;
    }

    bool downloadOrbitals(double* C, double* eps) const override
    {
        const int n = m_n;
        std::vector<double> Cc((size_t)n * n);
        if (cudaMemcpy(Cc.data(), m_C, Cc.size() * sizeof(double), cudaMemcpyDeviceToHost) != cudaSuccess
            || cudaMemcpy(eps, m_eps, (size_t)n * sizeof(double), cudaMemcpyDeviceToHost) != cudaSuccess)
            return false;
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) C[(size_t)i * n + j] = Cc[(size_t)j * n + i];
        return true;
    }

    std::string stats() const override
    {
        char buf[320];
        std::snprintf(buf, sizeof buf,
                      "%d resident iterations: Fock (J/K) %.1f ms, DIIS %.1f ms, eigensolve (F'=X^T F X + syevd) %.1f ms, "
                      "C = X C' + density %.1f ms",
                      m_iters, m_t_fock, m_t_diis, m_t_eig, m_t_dens);
        return buf;
    }

private:
    template <typename T>
    bool alloc(T** p, size_t count)
    {
        if (cudaMalloc(p, count * sizeof(T)) != cudaSuccess) { *p = nullptr; return false; }
        m_owned.push_back(*p);
        return true;
    }

    void release()
    {
        for (void* p : m_owned) cudaFree(p);
        m_owned.clear();
        cudaFree(m_histF); cudaFree(m_histP); cudaFree(m_histE);
        m_histF = m_histP = m_histE = nullptr;
    }

    /// C = op(A) op(B), n x n column-major.
    bool gemm(cublasOperation_t ta, cublasOperation_t tb, const double* A, const double* B, double* C)
    {
        const double one = 1.0, zero = 0.0;
        return cublasDgemm(m_cublas, ta, tb, m_n, m_n, m_n, &one, A, m_n, B, m_n, &zero, C, m_n)
            == CUBLAS_STATUS_SUCCESS;
    }

    void accumulate(double& bucket, int e0, int e1)
    {
        float ms = 0.0f;
        if (cudaEventElapsedTime(&ms, m_ev[e0], m_ev[e1]) == cudaSuccess) bucket += ms;
    }

    /// Deterministic two-stage reductions for every spec; synchronous download of the results
    /// (and of the syevd info word when info = true).
    bool reducePairs(const std::vector<PairSpec>& specs, std::vector<double>& out, double& bucket, bool info = false)
    {
        const int np = (int)specs.size();
        if (np == 0 || np > 1024) return false;
        const auto t0 = std::chrono::steady_clock::now();
        if (cudaMemcpy(m_specs, specs.data(), np * sizeof(PairSpec), cudaMemcpyHostToDevice) != cudaSuccess) return false;
        const size_t nn = (size_t)m_n * m_n;
        pairDotStage1<<<dim3(kRedBlocks, np), kRedThreads>>>(m_specs, nn, m_partial);
        pairDotStage2<<<np, kRedThreads>>>(m_partial, kRedBlocks, m_out);
        if (cudaGetLastError() != cudaSuccess) return false;
        out.resize(np);
        if (cudaMemcpy(out.data(), m_out, np * sizeof(double), cudaMemcpyDeviceToHost) != cudaSuccess) return false;
        if (info && cudaMemcpy(&m_last_info, m_info, sizeof(int), cudaMemcpyDeviceToHost) != cudaSuccess) return false;
        bucket += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
        return true;
    }

    int m_n = 0, m_nocc = 0, m_nslots = 0, m_lwork = 0, m_iters = 0, m_last_info = 0;
    cublasHandle_t m_cublas = nullptr;
    cusolverDnHandle_t m_cusolver = nullptr;
    cudaEvent_t m_ev[6] = {};
    std::vector<void*> m_owned;
    double *m_H = nullptr, *m_S = nullptr, *m_X = nullptr, *m_P = nullptr, *m_Pn = nullptr, *m_Pref = nullptr,
           *m_dP = nullptr, *m_J = nullptr, *m_K = nullptr, *m_dJ = nullptr, *m_dK = nullptr, *m_F = nullptr,
           *m_HF = nullptr, *m_Fuse = nullptr, *m_T1 = nullptr, *m_T2 = nullptr, *m_C = nullptr, *m_eps = nullptr,
           *m_partial = nullptr, *m_out = nullptr, *m_coef = nullptr, *m_work = nullptr;
    double *m_histF = nullptr, *m_histP = nullptr, *m_histE = nullptr;
    int *m_info = nullptr, *m_slots = nullptr;
    PairSpec* m_specs = nullptr;
    bool m_pending_diis = false;
    double m_t_fock = 0.0, m_t_diis = 0.0, m_t_eig = 0.0, m_t_dens = 0.0;
    std::string m_name = "cuda";
};

void setWhyScf(char* why, int whylen, const std::string& s)
{
    if (why && whylen > 0) std::snprintf(why, (size_t)whylen, "%s", s.c_str());
}

} // namespace

// Resident SCF state for n basis functions on `device` (QMGpuSCF*), or nullptr with a reason.
// Claude Generated (Oct 2026, Q-G4).
extern "C" void* curcuma_cuda_create_qm_scf(int n, int device, char* why, int whylen)
{
    if (n <= 0) { setWhyScf(why, whylen, "empty basis"); return nullptr; }
    if (device >= 0 && cudaSetDevice(device) != cudaSuccess) { setWhyScf(why, whylen, "cudaSetDevice failed"); return nullptr; }
    auto s = std::make_unique<QmGpuScfCuda>();
    std::string err;
    if (!s->init(n, device, err)) {
        cudaGetLastError();
        setWhyScf(why, whylen, err);
        return nullptr;
    }
    return static_cast<qmint::QMGpuSCF*>(s.release());
}
