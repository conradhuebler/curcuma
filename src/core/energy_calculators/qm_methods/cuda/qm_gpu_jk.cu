/*
 * < CUDA J/K build from a device-resident ERI tensor (Q-G1) >
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0
 */

// Claude Generated (Oct 2026): the stored (n^2 x n^2) ERI tensor is uploaded once per
// geometry (Q-G1) or built on the device (Q-G2, qm_gpu_eri.cu); every Fock build then only
// moves P (n^2) in and J, K (2 n^2) out.
// J_mn = sum_t E[(mn), t] P[t]            (t = lam*n + sig, contiguous)
// K_mn = sum_{l,s} E[m,l,n,s] P[l,s]      (s contiguous, l strided by n^2)
// One thread block per (m,n); FP64 accumulation, fixed reduction order (deterministic).

#include "qm_gpu_jk_cuda.h"

#include <cuda_runtime.h>

#include <cstdio>
#include <memory>
#include <string>
#include <vector>

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

/// Stage 6a (open shell, Claude Generated Nov 2026): one pass over the stored tensor for two spin densities,
///   J_mn = sum_t E[(mn), t] Pt[t] (Pt = Pa + Pb),  Ka_mn = sum_ls E[m,l,n,s] Pa[l,s],  Kb_mn = (same with Pb).
/// Same block/thread layout and reduction order as jkKernel; each tensor element is read once for three contractions.
__global__ void jkSpinKernel(const double* __restrict__ E, const double* __restrict__ Pt, const double* __restrict__ Pa,
                             const double* __restrict__ Pb, double* __restrict__ J, double* __restrict__ Ka,
                             double* __restrict__ Kb, int n)
{
    const int mn = blockIdx.x;
    const int m = mn / n, nu = mn % n;
    const size_t n2 = (size_t)n * n;
    const double* rowJ = E + (size_t)mn * n2;
    const double* baseK = E + (size_t)m * n2 * n + (size_t)nu * n;

    double sj = 0.0, ska = 0.0, skb = 0.0;
    for (size_t t = threadIdx.x; t < n2; t += blockDim.x) {
        const int l = (int)(t / n), s = (int)(t % n);
        const double ek = baseK[(size_t)l * n2 + s];
        sj += rowJ[t] * Pt[t];
        ska += ek * Pa[t];
        skb += ek * Pb[t];
    }
    __shared__ double rj[256], ra[256], rb[256];
    rj[threadIdx.x] = sj;
    ra[threadIdx.x] = ska;
    rb[threadIdx.x] = skb;
    __syncthreads();
    for (int o = blockDim.x / 2; o > 0; o >>= 1) {
        if (threadIdx.x < o) {
            rj[threadIdx.x] += rj[threadIdx.x + o];
            ra[threadIdx.x] += ra[threadIdx.x + o];
            rb[threadIdx.x] += rb[threadIdx.x + o];
        }
        __syncthreads();
    }
    if (threadIdx.x == 0) {
        J[mn] = rj[0];
        Ka[mn] = ra[0];
        Kb[mn] = rb[0];
    }
}

} // namespace

namespace qmgpu {

void warmJkKernel()
{
    cudaFuncAttributes a{};
    cudaFuncGetAttributes(&a, jkKernel);
}

bool QmGpuJkCuda::allocWork(int device)
{
    const size_t n2 = (size_t)m_n * m_n;
    if (cudaMalloc(&m_P, n2 * sizeof(double)) != cudaSuccess
        || cudaMalloc(&m_J, n2 * sizeof(double)) != cudaSuccess
        || cudaMalloc(&m_K, n2 * sizeof(double)) != cudaSuccess)
        return false;
    cudaDeviceProp prop{};
    cudaGetDeviceProperties(&prop, device >= 0 ? device : 0);
    m_name = std::string("cuda:") + prop.name;
    return true;
}

QmGpuJkCuda::QmGpuJkCuda(const double* eri, int n, int device)
    : m_n(n)
{
    if (device >= 0 && cudaSetDevice(device) != cudaSuccess) return;
    const size_t n2 = (size_t)n * n;
    if (cudaMalloc(&m_E, n2 * n2 * sizeof(double)) != cudaSuccess) { m_E = nullptr; return; }
    if (!allocWork(device)) { release(); return; }
    if (cudaMemcpy(m_E, eri, n2 * n2 * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess)
        release();
}

QmGpuJkCuda::QmGpuJkCuda(double* d_eri, int n, int device, bool /*adopt*/)
    : m_n(n)
{
    m_E = d_eri;
    if (!allocWork(device)) release();
}

QmGpuJkCuda::~QmGpuJkCuda() { release(); }

bool QmGpuJkCuda::build(const double* P, double* J, double* K)
{
    if (!m_E) return false;
    const size_t n2 = (size_t)m_n * m_n;
    if (cudaMemcpy(m_P, P, n2 * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess) return false;
    jkKernel<<<(unsigned)n2, 256>>>(m_E, m_P, m_J, m_K, m_n);
    if (cudaGetLastError() != cudaSuccess) return false;
    if (cudaMemcpy(J, m_J, n2 * sizeof(double), cudaMemcpyDeviceToHost) != cudaSuccess) return false;
    return cudaMemcpy(K, m_K, n2 * sizeof(double), cudaMemcpyDeviceToHost) == cudaSuccess;
}

// Claude Generated (Oct 2026, Q-G4): device-pointer build for the resident SCF loop; the same
// kernel as build(), so J/K are bit-identical to the host-transfer variant.
bool QmGpuJkCuda::buildDevice(const double* dP, double* dJ, double* dK)
{
    if (!m_E) return false;
    const size_t n2 = (size_t)m_n * m_n;
    jkKernel<<<(unsigned)n2, 256>>>(m_E, dP, dJ, dK, m_n);
    if (cudaGetLastError() != cudaSuccess) return false;
    return cudaDeviceSynchronize() == cudaSuccess;
}

// Stage 6a (Claude Generated, Nov 2026): J[Pa + Pb], K[Pa], K[Pb] in one pass over the tensor (jkSpinKernel). Pt is
// formed on the host (the same sum the CPU path hands to buildCoulomb). Spin buffers are allocated on first use.
bool QmGpuJkCuda::buildSpin(const double* Pa, const double* Pb, double* J, double* Ka, double* Kb)
{
    if (!m_E) return false;
    const size_t n2 = (size_t)m_n * m_n;
    if (!m_S) {
        if (cudaMalloc(&m_S, 4 * n2 * sizeof(double)) != cudaSuccess) { m_S = nullptr; return false; }
    }
    double *dPa = m_S, *dPb = m_S + n2, *dKa = m_S + 2 * n2, *dKb = m_S + 3 * n2;
    std::vector<double> Pt(n2);
    for (size_t k = 0; k < n2; ++k) Pt[k] = Pa[k] + Pb[k];
    if (cudaMemcpy(m_P, Pt.data(), n2 * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess
        || cudaMemcpy(dPa, Pa, n2 * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess
        || cudaMemcpy(dPb, Pb, n2 * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess)
        return false;
    jkSpinKernel<<<(unsigned)n2, 256>>>(m_E, m_P, dPa, dPb, m_J, dKa, dKb, m_n);
    if (cudaGetLastError() != cudaSuccess) return false;
    return cudaMemcpy(J, m_J, n2 * sizeof(double), cudaMemcpyDeviceToHost) == cudaSuccess
        && cudaMemcpy(Ka, dKa, n2 * sizeof(double), cudaMemcpyDeviceToHost) == cudaSuccess
        && cudaMemcpy(Kb, dKb, n2 * sizeof(double), cudaMemcpyDeviceToHost) == cudaSuccess;
}

bool QmGpuJkCuda::downloadERI(double* host) const
{
    if (!m_E) return false;
    const size_t n2 = (size_t)m_n * m_n;
    return cudaMemcpy(host, m_E, n2 * n2 * sizeof(double), cudaMemcpyDeviceToHost) == cudaSuccess;
}

void QmGpuJkCuda::release()
{
    cudaFree(m_E); cudaFree(m_P); cudaFree(m_J); cudaFree(m_K); cudaFree(m_S);
    m_E = m_P = m_J = m_K = m_S = nullptr;
}

} // namespace qmgpu

// _spin (Nov 2026, stage 6a): QMGpuJK gained buildSpin(); entry renamed with the vtable (gpu_plugin.cpp).
extern "C" void* curcuma_cuda_create_qm_jk_spin(const double* eri, int n, int device)
{
    auto* b = new qmgpu::QmGpuJkCuda(eri, n, device);
    if (!b->ok()) { delete b; return nullptr; }
    return static_cast<qmint::QMGpuJK*>(b);
}
