/*
 * < CUDA integral-direct J/K in FP32 for the early SCF phase (Q-G5) >
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0
 */

// Claude Generated (Oct 20 2026, Q-G5 of docs/QM_GPU_ROADMAP.md).
//
// FP32 instances of the thread-per-quartet McMurchie-Davidson kernels of qm_gpu_direct_md_impl.cuh
// (RT = float) for the early phase of a direct SCF with -qm.gpu_eri_precision mixed (default auto =
// mixed on cards whose FP64 rate is at most 1/16 of FP32, e.g. GeForce). Semantics as the GFN2
// -scf_mixed_precision and the RI -qm.ri_mixed_precision: FP32 far from convergence, FP64 near it.
// The SCF driver (qm_scf.cpp, DirectFockCadence) switches back to FP64 with a FULL rebuild once
// ||dP|| < -qm.gpu_fp32_switch or the FP32 builds stall, and accepts convergence only on FP64 builds,
// so no FP32 contribution survives into the converged Fock matrix or the energy.
//
// Accuracy measures in the kernel (qm_gpu_direct_md_impl.cuh): the product-centre difference
// W = P - Q is formed in FP64 before the cast; Boys from an FP32 copy of the same grid/Taylor table;
// Hermite tables F and the density as FP32 copies; per quartet FP32 sums; FP64 atomics into Jt, Kt.
// Screening is unchanged (FP64 Schwarz factors and Pmax), so the FP32 build computes the same quartets.

#include "qm_gpu_direct_md_impl.cuh"

#include <string>
#include <vector>

namespace {

__global__ void toFloatKernel(const double* __restrict__ in, float* __restrict__ out, long long n)
{
    const long long k = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (k < n) out[k] = (float)in[k];
}

void toFloat(const double* in, float* out, long long n, cudaStream_t st)
{
    if (n > 0) toFloatKernel<<<(unsigned)((n + 255) / 256), 256, 0, st>>>(in, out, n);
}

} // namespace

namespace qmgpu {

bool md32Init(const DevTables& d, MdTables& m, std::vector<void*>& owned, std::string& why, size_t& bytes)
{
    bytes = 0;
    float *F32 = nullptr, *B32 = nullptr, *P32 = nullptr;
    const size_t n2 = (size_t)d.n_out * d.n_out;
    if ((m.nF > 0 && cudaMalloc(&F32, (size_t)m.nF * sizeof(float)) != cudaSuccess)
        || cudaMalloc(&B32, (size_t)m.nBoys * sizeof(float)) != cudaSuccess
        || cudaMalloc(&P32, n2 * sizeof(float)) != cudaSuccess) {
        cudaGetLastError();
        if (F32) cudaFree(F32);
        if (B32) cudaFree(B32);
        why = "device memory for the FP32 tables";
        return false;
    }
    if (F32) owned.push_back(F32);
    owned.push_back(B32);
    owned.push_back(P32);
    toFloat(m.F, F32, m.nF, 0);
    toFloat(d.boys_f, B32, m.nBoys, 0);
    if (cudaDeviceSynchronize() != cudaSuccess || cudaGetLastError() != cudaSuccess) {
        why = "FP32 table conversion failed";
        return false;
    }
    m.F32 = F32;
    m.boys_f32 = B32;
    m.P32 = P32;
    bytes = ((size_t)m.nF + (size_t)m.nBoys + n2) * sizeof(float);
    return true;
}

void md32Density(const double* P, const MdTables& m, size_t n2, cudaStream_t st)
{
    toFloat(P, const_cast<float*>(m.P32), (long long)n2, st);
}

bool mdLaunch32(int LA, int LB, int LC, int LD, const DevTables& d, const MdTables& m, const QSrc& src, long long nq,
                cudaStream_t st)
{
    if (nq <= 0) return true;
#define QG5_CASE(a, b, c, e) \
    case ((a * 3 + b) * 3 + c) * 3 + e: return mdLaunchT<float, a, b, c, e>(d, m, src, nq, st);
    switch (((LA * 3 + LB) * 3 + LC) * 3 + LD) {
        // every class with bra pair class X <= ket pair class Y (the launches of qm_gpu_direct.cu)
        QG5_CASE(0, 0, 0, 0) QG5_CASE(0, 0, 0, 1) QG5_CASE(0, 0, 0, 2) QG5_CASE(0, 0, 1, 0) QG5_CASE(0, 0, 1, 1) QG5_CASE(0, 0, 1, 2) QG5_CASE(0, 0, 2, 0) QG5_CASE(0, 0, 2, 1) QG5_CASE(0, 0, 2, 2)
        QG5_CASE(0, 1, 0, 1) QG5_CASE(0, 1, 0, 2) QG5_CASE(0, 1, 1, 0) QG5_CASE(0, 1, 1, 1) QG5_CASE(0, 1, 1, 2) QG5_CASE(0, 1, 2, 0) QG5_CASE(0, 1, 2, 1) QG5_CASE(0, 1, 2, 2)
        QG5_CASE(0, 2, 0, 2) QG5_CASE(0, 2, 1, 0) QG5_CASE(0, 2, 1, 1) QG5_CASE(0, 2, 1, 2) QG5_CASE(0, 2, 2, 0) QG5_CASE(0, 2, 2, 1) QG5_CASE(0, 2, 2, 2)
        QG5_CASE(1, 0, 1, 0) QG5_CASE(1, 0, 1, 1) QG5_CASE(1, 0, 1, 2) QG5_CASE(1, 0, 2, 0) QG5_CASE(1, 0, 2, 1) QG5_CASE(1, 0, 2, 2)
        QG5_CASE(1, 1, 1, 1) QG5_CASE(1, 1, 1, 2) QG5_CASE(1, 1, 2, 0) QG5_CASE(1, 1, 2, 1) QG5_CASE(1, 1, 2, 2)
        QG5_CASE(1, 2, 1, 2) QG5_CASE(1, 2, 2, 0) QG5_CASE(1, 2, 2, 1) QG5_CASE(1, 2, 2, 2)
        QG5_CASE(2, 0, 2, 0) QG5_CASE(2, 0, 2, 1) QG5_CASE(2, 0, 2, 2)
        QG5_CASE(2, 1, 2, 1) QG5_CASE(2, 1, 2, 2)
        QG5_CASE(2, 2, 2, 2)
    default: return false;
    }
#undef QG5_CASE
}

} // namespace qmgpu
