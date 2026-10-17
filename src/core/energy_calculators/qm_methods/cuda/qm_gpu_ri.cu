/*
 * < CUDA density-fitting (RI-J / RI-JK) contractions of the native QM engine (WP10) >
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0
 */

// Claude Generated (Oct 2026, WP10). AI-generated, machine-tested.
//
// The fitted 3-index tensor T (nfit x np, np = n(n+1)/2, built on the host by qmri::RIJK::setup) is
// uploaded once per geometry. Row Q of the row-major host array is the PACKED lower triangle of the
// symmetric n x n block T_Q (column m(m+1)/2 + k for m >= k; packed since Oct 2026, half the memory and
// PCIe traffic of the former n^2 layout). As a column-major array the same memory is
// Mp = [T_0, T_1, ...] (np x nfit). The contractions of one J/K build:
//
//   g = Mp^T w                             cublasDgemv (T)   w = packed weights of P (P_mm, P_mk + P_km)
//   J_p = Mp c                             cublasDgemv (N)   (c = g for B, J^-1 g in mode j; host unpacks J)
//   K = sum_Q (T_Q V+)(T_Q V+)^T - (T_Q V-)(T_Q V-)^T, per chunk of nq rows:
//       U = full blocks [T_Q0 .. T_Q0+nq-1] (n x nq n)  unpackPairs kernel
//       Y = U^T V (row block q = T_q V)    cublasDgemm (T)
//       K += X X^T                         cublasDsyrk (lower), +1 for V+, -1 for V-
//
// V (n x r, column-major) are the eigenvectors of P scaled by sqrt|lambda| or sqrt(2) C_occ. FP64 only.

#include "../qm_gpu_ri_backend.h"

#include <cublas_v2.h>
#include <cuda_runtime.h>

#include <algorithm>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

namespace {

/// U[q][i][j] = T_{Q0+q}[max(i,j)][min(i,j)] from the packed rows (blocks are symmetric, so row- and
/// column-major are the same). One thread per output element, grid-stride. Claude Generated (Oct 2026).
__global__ void unpackPairs(const double* __restrict__ Tp, double* __restrict__ U, int n, long long np, int nq)
{
    const long long n2 = (long long)n * n, total = (long long)nq * n2;
    for (long long idx = (long long)blockIdx.x * blockDim.x + threadIdx.x; idx < total; idx += (long long)gridDim.x * blockDim.x) {
        const long long q = idx / n2, rem = idx - q * n2;
        const int i = (int)(rem / n), j = (int)(rem - (long long)i * n);
        const int a = i > j ? i : j, b = i > j ? j : i;
        U[idx] = Tp[q * np + (long long)a * (a + 1) / 2 + b];
    }
}

class QmGpuRiCuda : public qmint::QMGpuRI {
public:
    ~QmGpuRiCuda() override
    {
        for (void* p : m_owned) cudaFree(p);
        if (m_cublas) cublasDestroy(m_cublas);
    }

    bool init(const double* T, int nfit, int n, int device, std::string& err)
    {
        m_nfit = nfit;
        m_n = n;
        m_np = (long long)n * (n + 1) / 2;
        m_device = device;
        const size_t n2 = (size_t)n * n, np = (size_t)m_np;
        size_t free_b = 0, total_b = 0;
        if (cudaMemGetInfo(&free_b, &total_b) != cudaSuccess) { err = "cudaMemGetInfo failed"; return false; }
        const size_t tensor_b = np * (size_t)nfit * sizeof(double);
        // chunk buffers for X and the unpacked blocks U: each at most 512 MB or a fifth of what is left, at least
        // one block (n x n, r <= n); plus 64 MB slack for cuBLAS and the small arrays
        const size_t blk_b = n2 * sizeof(double);
        if (tensor_b + 2 * blk_b + (size_t)64 * 1024 * 1024 > free_b) {
            char b[256];
            std::snprintf(b, sizeof b, "RI tensor %.1f MB (packed) does not fit into %.1f MB free device memory", tensor_b / 1e6,
                          free_b / 1e6);
            err = b;
            return false;
        }
        const size_t left = free_b - tensor_b - (size_t)64 * 1024 * 1024;
        m_xbuf_b = std::max(std::min<size_t>((size_t)512 * 1024 * 1024, left / 5), blk_b);
        m_ubuf_b = m_xbuf_b;
        if (cublasCreate(&m_cublas) != CUBLAS_STATUS_SUCCESS) { err = "cublasCreate failed"; return false; }
        if (!alloc(&m_T, np * nfit) || !alloc(&m_P, np) || !alloc(&m_J, np) || !alloc(&m_K, n2) || !alloc(&m_g, nfit)
            || !alloc(&m_V, n2) || !alloc(&m_X, m_xbuf_b / sizeof(double)) || !alloc(&m_U, m_ubuf_b / sizeof(double))) {
            err = "device allocation failed";
            return false;
        }
        if (cudaMemcpy(m_T, T, tensor_b, cudaMemcpyHostToDevice) != cudaSuccess) { err = "tensor upload failed"; return false; }
        char b[128];
        std::snprintf(b, sizeof b, "cuda:%d (cuBLAS FP64)", device);
        m_name = b;
        return true;
    }

    std::string name() const override { return m_name; }

    bool contract(const double* w, double* g) override
    {
        const int np = (int)m_np;
        const double one = 1.0, zero = 0.0;
        if (cudaMemcpy(m_P, w, sizeof(double) * np, cudaMemcpyHostToDevice) != cudaSuccess) return false;
        if (cublasDgemv(m_cublas, CUBLAS_OP_T, np, m_nfit, &one, m_T, np, m_P, 1, &zero, m_g, 1) != CUBLAS_STATUS_SUCCESS) return false;
        return cudaMemcpy(g, m_g, sizeof(double) * m_nfit, cudaMemcpyDeviceToHost) == cudaSuccess;
    }

    bool coulomb(const double* c, double* J) override
    {
        const int np = (int)m_np;
        const double one = 1.0, zero = 0.0;
        if (cudaMemcpy(m_g, c, sizeof(double) * m_nfit, cudaMemcpyHostToDevice) != cudaSuccess) return false;
        if (cublasDgemv(m_cublas, CUBLAS_OP_N, np, m_nfit, &one, m_T, np, m_g, 1, &zero, m_J, 1) != CUBLAS_STATUS_SUCCESS) return false;
        return cudaMemcpy(J, m_J, sizeof(double) * np, cudaMemcpyDeviceToHost) == cudaSuccess;
    }

    bool exchange(const double* V, int r, int npos, double* K) override
    {
        const int n = m_n;
        const size_t n2 = (size_t)n * n;
        if (r > n || r < 0 || npos < 0 || npos > r) return false;
        if (cudaMemcpy(m_V, V, sizeof(double) * (size_t)n * r, cudaMemcpyHostToDevice) != cudaSuccess) return false;
        if (cudaMemset(m_K, 0, sizeof(double) * n2) != cudaSuccess) return false;
        const double one = 1.0, zero = 0.0;
        // per chunk: unpack the rows once, then the positive columns with +1 and the negative ones with -1
        const int maxcols = std::max(1, std::max(npos, r - npos));
        const int chunk = (int)std::max<size_t>(
            1, std::min({ (size_t)m_nfit, m_xbuf_b / sizeof(double) / ((size_t)n * maxcols), m_ubuf_b / sizeof(double) / n2 }));
        for (int Q0 = 0; Q0 < m_nfit; Q0 += chunk) {
            const int nq = std::min(chunk, m_nfit - Q0);
            const long long total = (long long)nq * (long long)n2;
            const int threads = 256;
            const int blocks = (int)std::min<long long>((total + threads - 1) / threads, 65535LL * 8);
            unpackPairs<<<blocks, threads>>>(m_T + (size_t)Q0 * (size_t)m_np, m_U, n, m_np, nq);
            if (cudaGetLastError() != cudaSuccess) return false;
            for (int pass = 0; pass < 2; ++pass) {
                const int cols = pass == 0 ? npos : r - npos;
                if (cols == 0) continue;
                const double* Vp = m_V + (size_t)(pass == 0 ? 0 : npos) * n;
                const double alpha = pass == 0 ? 1.0 : -1.0;
                // Y = U^T V (nq n x cols), U = [T_Q0, ...] (n x nq n, column-major, blocks symmetric): row block q
                // of Y is X_q = T_q V. Read as column-major n x (nq cols) the same memory is X = [.. (X_q)_{:,k} ..].
                if (cublasDgemm(m_cublas, CUBLAS_OP_T, CUBLAS_OP_N, nq * n, cols, n, &one, m_U, n, Vp, n, &zero, m_X, nq * n)
                    != CUBLAS_STATUS_SUCCESS)
                    return false;
                // K += alpha X X^T, X (n x nq*cols)
                if (cublasDsyrk(m_cublas, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_N, n, nq * cols, &alpha, m_X, n, &one, m_K, n)
                    != CUBLAS_STATUS_SUCCESS)
                    return false;
            }
        }
        std::vector<double> Kc(n2);
        if (cudaMemcpy(Kc.data(), m_K, sizeof(double) * n2, cudaMemcpyDeviceToHost) != cudaSuccess) return false;
        // column-major lower triangle -> full symmetric row-major
        for (int j = 0; j < n; ++j)
            for (int i = j; i < n; ++i) {
                const double v = Kc[(size_t)j * n + i];
                K[(size_t)i * n + j] = v;
                K[(size_t)j * n + i] = v;
            }
        return true;
    }

private:
    bool alloc(double** p, size_t count)
    {
        void* q = nullptr;
        if (cudaMalloc(&q, sizeof(double) * std::max<size_t>(1, count)) != cudaSuccess) return false;
        m_owned.push_back(q);
        *p = static_cast<double*>(q);
        return true;
    }

    int m_nfit = 0, m_n = 0, m_device = 0;
    long long m_np = 0;
    size_t m_xbuf_b = 0, m_ubuf_b = 0;
    cublasHandle_t m_cublas = nullptr;
    std::vector<void*> m_owned;
    double *m_T = nullptr, *m_P = nullptr, *m_J = nullptr, *m_K = nullptr, *m_g = nullptr, *m_V = nullptr, *m_X = nullptr, *m_U = nullptr;
    std::string m_name = "cuda";
};

void setWhyRi(char* why, int whylen, const std::string& s)
{
    if (why && whylen > 0) std::snprintf(why, (size_t)whylen, "%s", s.c_str());
}

} // namespace

// Device RI contractions for the fitted tensor T (nfit x n(n+1)/2 packed, row-major host array), or nullptr
// with a reason. Claude Generated (Oct 2026, WP10; packed layout and entry name since Oct 2026).
extern "C" void* curcuma_cuda_create_qm_ri_packed(const double* T, int nfit, int n, int device, char* why, int whylen)
{
    if (!T || nfit <= 0 || n <= 0) { setWhyRi(why, whylen, "empty RI tensor"); return nullptr; }
    if (device >= 0 && cudaSetDevice(device) != cudaSuccess) { setWhyRi(why, whylen, "cudaSetDevice failed"); return nullptr; }
    auto s = std::make_unique<QmGpuRiCuda>();
    std::string err;
    if (!s->init(T, nfit, n, device, err)) {
        cudaGetLastError();
        setWhyRi(why, whylen, err);
        return nullptr;
    }
    return static_cast<qmint::QMGpuRI*>(s.release());
}
