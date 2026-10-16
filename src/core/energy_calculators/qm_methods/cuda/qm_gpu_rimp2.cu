/*
 * < CUDA RI-MP2 pair energies of the native QM engine >
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0
 */

// Claude Generated (Oct 2026, RI-MP2). AI-generated, machine-tested.
//
// The fitted factors B[i][a][P] (built on the host, qm_rimp2.cpp) are uploaded once; per batch of
// ib active occupied orbitals one cuBLAS DGEMM (FP64) forms all (ia|jb) of the batch. cuBLAS is
// column-major: the row-major host array B (o v x naux) is the column-major matrix Bc (naux x o v),
//
//   K (o v x ib v) = Bc^T (o v x naux) * Bc[:, i0 v : (i0+ib) v]   ->  K(jb, ia) at jb + ia * o v,
//
// i.e. the row-major layout [i][a][j][b] of the canonical code, so (ib|ja) = K[i][b][j][a] lies in
// the same batch and the pair-energy kernel is the one of cuda/qm_gpu_mp2.cu (copied: a shared
// header would pull that TU into every rebuild). Deterministic for given o, v, naux and batch size
// (block partials summed on the host in block order). No use of the i <-> j symmetry here (the CPU
// path uses it); the device GEMM is fast enough that the factor 2 was not worth the extra kernel.

#include "../qm_gpu_rimp2_backend.h"

#include <cublas_v2.h>
#include <cuda_runtime.h>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <string>
#include <vector>

namespace {

constexpr int kRedBlocks = 240;
constexpr int kRedThreads = 256;

/// Partial E_OS / E_SS sums of one batch (same expression and order as mp2PairEnergyKernel).
__global__ void riMp2PairEnergyKernel(const double* __restrict__ K, const double* __restrict__ eo,
                                      const double* __restrict__ ev, int i0, int ib, int o, int v,
                                      double* __restrict__ part_os, double* __restrict__ part_ss)
{
    __shared__ double s_os[kRedThreads];
    __shared__ double s_ss[kRedThreads];
    const long long total = (long long)ib * v * o * v;
    double os = 0.0, ss = 0.0;
    for (long long idx = (long long)blockIdx.x * blockDim.x + threadIdx.x; idx < total;
         idx += (long long)gridDim.x * blockDim.x) {
        const int b = (int)(idx % v);
        long long r = idx / v;
        const int j = (int)(r % o);
        r /= o;
        const int a = (int)(r % v);
        const int i = (int)(r / v);
        const double t = K[idx];
        const double tx = K[(((long long)i * v + b) * o + j) * v + a];  // (ib|ja)
        const double inv_d = 1.0 / (eo[i0 + i] + eo[j] - ev[a] - ev[b]);
        os += t * t * inv_d;
        ss += t * (t - tx) * inv_d;
    }
    s_os[threadIdx.x] = os;
    s_ss[threadIdx.x] = ss;
    __syncthreads();
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (threadIdx.x < s) {
            s_os[threadIdx.x] += s_os[threadIdx.x + s];
            s_ss[threadIdx.x] += s_ss[threadIdx.x + s];
        }
        __syncthreads();
    }
    if (threadIdx.x == 0) {
        part_os[blockIdx.x] = s_os[0];
        part_ss[blockIdx.x] = s_ss[0];
    }
}

void setWhy(char* why, int len, const std::string& s)
{
    if (why && len > 0) std::snprintf(why, (size_t)len, "%s", s.c_str());
}

/// Device buffers of one call; freed in the destructor (also on every error return).
struct RiMp2Work {
    double *B = nullptr, *eo = nullptr, *ev = nullptr, *K = nullptr, *part = nullptr;
    cublasHandle_t h = nullptr;
    ~RiMp2Work()
    {
        for (double* p : { B, eo, ev, K, part })
            if (p) cudaFree(p);
        if (h) cublasDestroy(h);
    }
};

double msSince(std::chrono::steady_clock::time_point t0)
{
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}

}  // namespace

// RI-MP2 pair energies on the device; returns 1 on success, 0 with a reason in `why`.
// Claude Generated (Oct 2026).
extern "C" int curcuma_cuda_qm_rimp2_energy(const qmint::QMGpuRIMP2Input* in, qmint::QMGpuRIMP2Output* out, char* why,
                                            int whylen)
{
    if (!in || !out) { setWhy(why, whylen, "null argument"); return 0; }
    const int o = in->o, v = in->v, naux = in->naux;
    *out = qmint::QMGpuRIMP2Output();
    if (o <= 0 || v <= 0) return 1;
    if (!in->B || naux <= 0) { setWhy(why, whylen, "no B tensor"); return 0; }
    const size_t ov = (size_t)o * v;
    if (ov > (size_t)2147483647) {
        setWhy(why, whylen, "dimensions too large for 32-bit cuBLAS arguments");
        return 0;
    }
    if (in->device >= 0 && cudaSetDevice(in->device) != cudaSuccess) { setWhy(why, whylen, "cudaSetDevice failed"); return 0; }

    RiMp2Work w;
    auto fail = [&](const std::string& what) {
        cudaGetLastError();  // clear a sticky non-fatal error for the CPU fallback
        setWhy(why, whylen, what);
        return 0;
    };
#define RIMP2_CUDA(call, what) do { if ((call) != cudaSuccess) return fail(std::string(what) + ": " + cudaGetErrorString(cudaGetLastError())); } while (0)
#define RIMP2_CUBLAS(call, what) do { if ((call) != CUBLAS_STATUS_SUCCESS) return fail(what); } while (0)

    const auto t0 = std::chrono::steady_clock::now();
    if (cublasCreate(&w.h) != CUBLAS_STATUS_SUCCESS) return fail("cublasCreate failed");
    cublasSetMathMode(w.h, CUBLAS_DEFAULT_MATH);  // plain FP64

    // Memory: B (o v naux) + one batch of ib v o v doubles.
    size_t free_b = 0, total_b = 0;
    RIMP2_CUDA(cudaMemGetInfo(&free_b, &total_b), "cudaMemGetInfo");
    double avail = (double)free_b;
    if (in->max_device_mb > 0.0) avail = std::min(avail, in->max_device_mb * 1.0e6);
    avail -= (double)(ov * naux + o + v + 2 * kRedBlocks) * 8.0;
    avail *= 0.9;  // cuBLAS workspace and fragmentation
    const double per_i = (double)v * ov * 8.0;
    if (in->batch_mb > 0.0) avail = std::min(avail, std::max(per_i, in->batch_mb * 1.0e6));  // at least one i
    const int batch = (int)std::min<double>(o, std::floor(avail / per_i));
    if (batch < 1)
        return fail("device memory too small (" + std::to_string((long long)(free_b / 1000000)) + " MB free, need "
                    + std::to_string((long long)(((double)ov * naux * 8.0 + per_i) / 1.0e6)) + " MB)");
    if ((size_t)batch * v > (size_t)2147483647) return fail("batch too large for 32-bit cuBLAS arguments");

    RIMP2_CUDA(cudaMalloc(&w.B, ov * naux * sizeof(double)), "cudaMalloc B");
    RIMP2_CUDA(cudaMalloc(&w.eo, (size_t)o * sizeof(double)), "cudaMalloc eo");
    RIMP2_CUDA(cudaMalloc(&w.ev, (size_t)v * sizeof(double)), "cudaMalloc ev");
    RIMP2_CUDA(cudaMalloc(&w.part, 2 * kRedBlocks * sizeof(double)), "cudaMalloc partials");
    RIMP2_CUDA(cudaMalloc(&w.K, (size_t)batch * v * ov * sizeof(double)), "cudaMalloc batch buffer");
    RIMP2_CUDA(cudaMemcpy(w.B, in->B, ov * naux * sizeof(double), cudaMemcpyHostToDevice), "upload B");
    RIMP2_CUDA(cudaMemcpy(w.eo, in->eo, (size_t)o * sizeof(double), cudaMemcpyHostToDevice), "upload eo");
    RIMP2_CUDA(cudaMemcpy(w.ev, in->ev, (size_t)v * sizeof(double), cudaMemcpyHostToDevice), "upload ev");
    RIMP2_CUDA(cudaDeviceSynchronize(), "upload sync");
    out->t_upload_ms = msSince(t0);

    const double one = 1.0, zero = 0.0;
    std::vector<double> part(2 * kRedBlocks);
    double e_os = 0.0, e_ss = 0.0;
    int nbatches = 0;
    for (int i0 = 0; i0 < o; i0 += batch) {
        const int ib = std::min(batch, o - i0);
        ++nbatches;
        const auto tg = std::chrono::steady_clock::now();
        RIMP2_CUBLAS(cublasDgemm(w.h, CUBLAS_OP_T, CUBLAS_OP_N, (int)ov, ib * v, naux, &one, w.B, naux,
                                 w.B + (size_t)i0 * v * naux, naux, &zero, w.K, (int)ov),
                     "(ia|jb) DGEMM");
        RIMP2_CUDA(cudaDeviceSynchronize(), "GEMM");
        out->t_gemm_ms += msSince(tg);

        const auto te = std::chrono::steady_clock::now();
        riMp2PairEnergyKernel<<<kRedBlocks, kRedThreads>>>(w.K, w.eo, w.ev, i0, ib, o, v, w.part, w.part + kRedBlocks);
        RIMP2_CUDA(cudaGetLastError(), "pair-energy kernel launch");
        RIMP2_CUDA(cudaMemcpy(part.data(), w.part, part.size() * sizeof(double), cudaMemcpyDeviceToHost), "partials download");
        for (int k = 0; k < kRedBlocks; ++k) {
            e_os += part[k];
            e_ss += part[kRedBlocks + k];
        }
        out->t_energy_ms += msSince(te);
    }
#undef RIMP2_CUDA
#undef RIMP2_CUBLAS
    out->e_os = e_os;
    out->e_ss = e_ss;
    out->batch = batch;
    out->nbatches = nbatches;
    return 1;
}
