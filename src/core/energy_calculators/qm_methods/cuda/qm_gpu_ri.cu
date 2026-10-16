/*
 * < CUDA density-fitting (RI-J / RI-JK) contractions of the native QM engine (WP10) >
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0
 */

// Claude Generated (Oct 2026, WP10). AI-generated, machine-tested.
//
// The fitted 3-index tensor T (nfit x n^2, built on the host by qmri::RIJK::setup) is uploaded once
// per geometry. Row Q of the row-major host array is the symmetric n x n block T_Q; as a
// column-major array the same memory is M = [vec T_0, vec T_1, ...] (n^2 x nfit), and every block is
// also a column-major n x n matrix (symmetric). The three contractions of one J/K build:
//
//   g = M^T vec(P)                         cublasDgemv (T)          (P symmetric: layout-free)
//   vec(J) = M c                           cublasDgemv (N)          (c = g for B, J^-1 g in mode j)
//   K = sum_Q (T_Q V+)(T_Q V+)^T - (T_Q V-)(T_Q V-)^T
//       Y = M_chunk^T V (row block q = T_q V)  cublasDgemm (T)
//       K += X X^T                         cublasDsyrk (lower), +1 for V+, -1 for V-
//
// V (n x r, column-major) are the eigenvectors of P scaled by sqrt|lambda| (host, n^3). FP64 only.

#include "../qm_gpu_ri_backend.h"

#include <cublas_v2.h>
#include <cuda_runtime.h>

#include <algorithm>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

namespace {

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
        m_device = device;
        const size_t n2 = (size_t)n * n;
        size_t free_b = 0, total_b = 0;
        if (cudaMemGetInfo(&free_b, &total_b) != cudaSuccess) { err = "cudaMemGetInfo failed"; return false; }
        const size_t tensor_b = n2 * (size_t)nfit * sizeof(double);
        // chunk buffer for X: at most 512 MB or a quarter of what is left
        const size_t left = free_b > tensor_b ? free_b - tensor_b : 0;
        if (tensor_b + (size_t)64 * 1024 * 1024 > free_b) {
            char b[256];
            std::snprintf(b, sizeof b, "RI tensor %.1f MB does not fit into %.1f MB free device memory", tensor_b / 1e6, free_b / 1e6);
            err = b;
            return false;
        }
        // at least one block X_Q (n x r, r <= n) must fit
        m_xbuf_b = std::max(std::min<size_t>((size_t)512 * 1024 * 1024, left / 4), n2 * sizeof(double));
        if (cublasCreate(&m_cublas) != CUBLAS_STATUS_SUCCESS) { err = "cublasCreate failed"; return false; }
        if (!alloc(&m_T, n2 * nfit) || !alloc(&m_P, n2) || !alloc(&m_J, n2) || !alloc(&m_K, n2) || !alloc(&m_g, nfit)
            || !alloc(&m_V, n2) || !alloc(&m_X, m_xbuf_b / sizeof(double))) {
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

    bool contract(const double* P, double* g) override
    {
        const int n2 = m_n * m_n;
        const double one = 1.0, zero = 0.0;
        if (cudaMemcpy(m_P, P, sizeof(double) * n2, cudaMemcpyHostToDevice) != cudaSuccess) return false;
        if (cublasDgemv(m_cublas, CUBLAS_OP_T, n2, m_nfit, &one, m_T, n2, m_P, 1, &zero, m_g, 1) != CUBLAS_STATUS_SUCCESS) return false;
        return cudaMemcpy(g, m_g, sizeof(double) * m_nfit, cudaMemcpyDeviceToHost) == cudaSuccess;
    }

    bool coulomb(const double* c, double* J) override
    {
        const int n2 = m_n * m_n;
        const double one = 1.0, zero = 0.0;
        if (cudaMemcpy(m_g, c, sizeof(double) * m_nfit, cudaMemcpyHostToDevice) != cudaSuccess) return false;
        if (cublasDgemv(m_cublas, CUBLAS_OP_N, n2, m_nfit, &one, m_T, n2, m_g, 1, &zero, m_J, 1) != CUBLAS_STATUS_SUCCESS) return false;
        return cudaMemcpy(J, m_J, sizeof(double) * n2, cudaMemcpyDeviceToHost) == cudaSuccess;
    }

    bool exchange(const double* V, int r, int npos, double* K) override
    {
        const int n = m_n;
        const size_t n2 = (size_t)n * n;
        if (r > n || r < 0 || npos < 0 || npos > r) return false;
        if (cudaMemcpy(m_V, V, sizeof(double) * (size_t)n * r, cudaMemcpyHostToDevice) != cudaSuccess) return false;
        if (cudaMemset(m_K, 0, sizeof(double) * n2) != cudaSuccess) return false;
        const double one = 1.0, zero = 0.0;
        // two passes: positive columns with +1, negative columns with -1
        for (int pass = 0; pass < 2; ++pass) {
            const int cols = pass == 0 ? npos : r - npos;
            if (cols == 0) continue;
            const double* Vp = m_V + (size_t)(pass == 0 ? 0 : npos) * n;
            const double alpha = pass == 0 ? 1.0 : -1.0;
            const int chunk = (int)std::max<size_t>(1, std::min<size_t>((size_t)m_nfit, m_xbuf_b / sizeof(double) / ((size_t)n * cols)));
            for (int Q0 = 0; Q0 < m_nfit; Q0 += chunk) {
                const int nq = std::min(chunk, m_nfit - Q0);
                // Y = M^T V (nq n x cols), M = [T_Q0, ...] (n x nq n, column-major, blocks symmetric): row block q
                // of Y is X_q = T_q V. Read as column-major n x (nq cols) the same memory is X = [.. (X_q)_{:,k} ..].
                if (cublasDgemm(m_cublas, CUBLAS_OP_T, CUBLAS_OP_N, nq * n, cols, n, &one, m_T + (size_t)Q0 * n2, n, Vp, n, &zero,
                                m_X, nq * n)
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
    size_t m_xbuf_b = 0;
    cublasHandle_t m_cublas = nullptr;
    std::vector<void*> m_owned;
    double *m_T = nullptr, *m_P = nullptr, *m_J = nullptr, *m_K = nullptr, *m_g = nullptr, *m_V = nullptr, *m_X = nullptr;
    std::string m_name = "cuda";
};

void setWhyRi(char* why, int whylen, const std::string& s)
{
    if (why && whylen > 0) std::snprintf(why, (size_t)whylen, "%s", s.c_str());
}

} // namespace

// Device RI contractions for the fitted tensor T (nfit x n^2, row-major host array), or nullptr with a
// reason. Claude Generated (Oct 2026, WP10).
extern "C" void* curcuma_cuda_create_qm_ri(const double* T, int nfit, int n, int device, char* why, int whylen)
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
