/*
 * < CUDA J/K build from a device-resident ERI tensor (Q-G1) >
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0
 */

// Claude Generated (Oct 2026): the stored (n^2 x n^2) ERI tensor is uploaded once per
// geometry; every Fock build then only moves P (n^2) in and J, K (2 n^2) out.
// J_mn = sum_t E[(mn), t] P[t]            (t = lam*n + sig, contiguous)
// K_mn = sum_{l,s} E[m,l,n,s] P[l,s]      (s contiguous, l strided by n^2)
// One thread block per (m,n); FP64 accumulation, fixed reduction order (deterministic).

#include "../qm_gpu_backend.h"

#include <cuda_runtime.h>

#include <cstdio>
#include <memory>
#include <string>

namespace {

__global__ void jkKernel(const double* __restrict__ E, const double* __restrict__ P,
                         double* __restrict__ J, double* __restrict__ K, int n)
{
    const int mn = blockIdx.x;
    const int m = mn / n, nu = mn % n;
    const size_t n2 = (size_t)n * n;
    const double* rowJ = E + (size_t)mn * n2;
    const double* baseK = E + (size_t)m * n2 * n + (size_t)nu * n;

    double sj = 0.0, sk = 0.0;
    for (size_t t = threadIdx.x; t < n2; t += blockDim.x) {
        const int l = (int)(t / n), s = (int)(t % n);
        const double p = P[t];
        sj += rowJ[t] * p;
        sk += baseK[(size_t)l * n2 + s] * p;
    }
    __shared__ double rj[256], rk[256];
    rj[threadIdx.x] = sj;
    rk[threadIdx.x] = sk;
    __syncthreads();
    for (int o = blockDim.x / 2; o > 0; o >>= 1) {
        if (threadIdx.x < o) {
            rj[threadIdx.x] += rj[threadIdx.x + o];
            rk[threadIdx.x] += rk[threadIdx.x + o];
        }
        __syncthreads();
    }
    if (threadIdx.x == 0) {
        J[mn] = rj[0];
        K[mn] = rk[0];
    }
}

class QmGpuJkCuda final : public qmint::QMGpuJK {
public:
    QmGpuJkCuda(const double* eri, int n, int device)
        : m_n(n)
    {
        if (device >= 0 && cudaSetDevice(device) != cudaSuccess) return;
        const size_t n2 = (size_t)n * n;
        if (cudaMalloc(&m_E, n2 * n2 * sizeof(double)) != cudaSuccess) { m_E = nullptr; return; }
        if (cudaMalloc(&m_P, n2 * sizeof(double)) != cudaSuccess
            || cudaMalloc(&m_J, n2 * sizeof(double)) != cudaSuccess
            || cudaMalloc(&m_K, n2 * sizeof(double)) != cudaSuccess) { release(); return; }
        if (cudaMemcpy(m_E, eri, n2 * n2 * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess)
            release();
        cudaDeviceProp prop{};
        cudaGetDeviceProperties(&prop, device >= 0 ? device : 0);
        m_name = std::string("cuda:") + prop.name;
    }
    ~QmGpuJkCuda() override { release(); }
    bool ok() const { return m_E != nullptr; }

    bool build(const double* P, double* J, double* K) override
    {
        if (!m_E) return false;
        const size_t n2 = (size_t)m_n * m_n;
        if (cudaMemcpy(m_P, P, n2 * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess) return false;
        jkKernel<<<(unsigned)n2, 256>>>(m_E, m_P, m_J, m_K, m_n);
        if (cudaGetLastError() != cudaSuccess) return false;
        if (cudaMemcpy(J, m_J, n2 * sizeof(double), cudaMemcpyDeviceToHost) != cudaSuccess) return false;
        return cudaMemcpy(K, m_K, n2 * sizeof(double), cudaMemcpyDeviceToHost) == cudaSuccess;
    }
    std::string name() const override { return m_name; }
    int n() const override { return m_n; }

private:
    void release()
    {
        cudaFree(m_E); cudaFree(m_P); cudaFree(m_J); cudaFree(m_K);
        m_E = m_P = m_J = m_K = nullptr;
    }
    int m_n;
    double *m_E = nullptr, *m_P = nullptr, *m_J = nullptr, *m_K = nullptr;
    std::string m_name = "cuda";
};

} // namespace

extern "C" void* curcuma_cuda_create_qm_jk(const double* eri, int n, int device)
{
    auto* b = new QmGpuJkCuda(eri, n, device);
    if (!b->ok()) { delete b; return nullptr; }
    return static_cast<qmint::QMGpuJK*>(b);
}
