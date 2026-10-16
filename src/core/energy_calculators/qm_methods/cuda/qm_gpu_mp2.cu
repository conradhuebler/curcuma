/*
 * < CUDA canonical closed-shell MP2 of the native QM engine (WP11) >
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0
 */

// Claude Generated (Oct 2026, WP11 step 1, GPU variant). AI-generated, machine-tested.
//
// The same four one-index transformations as qmint::mp2EnergyHost (qm_mp2.cpp), as cuBLAS DGEMMs
// in FP64, batched over the active occupied index i. cuBLAS is column-major; with n the basis
// size, o the active occupied, v the virtual count:
//
//   AO tensor  E   (n^3 x n)  E(nu lam sig, mu)      = (mu nu|lam sig)  (row-major ERITensor)
//   step 1     X1  (n^3 x ib) X1(rest, i)            = sum_mu E(rest, mu) Co(mu, i)
//   step 2     W_i (n^2 x v)  W_i(lam sig, a)        = sum_nu X1_i(lam sig, nu) Cv(nu, a)   batched over i
//   step 3     Z_ia (n x o)   Z_ia(sig, j)           = sum_lam W_ia(sig, lam) Co(lam, j)    batched over (i,a)
//   step 4     K_ia (v x o)   K_ia(b, j) = (ia|jb)   = sum_sig Cv(sig, b) Z_ia(sig, j)       batched over (i,a)
//
// Z reuses the X1 buffer (v o <= n^2), K reuses the W buffer (o v <= n^2), so a batch of ib
// occupied orbitals needs ib (n^3 + n^2 v) doubles besides the tensor. K is stored as
// [i][a][j][b], the layout of the CPU code, so (ib|ja) = K[i][b][j][a] lies in the same batch.
// The pair energies E_OS = sum T^2/D, E_SS = sum T (T - T')/D are reduced with a fixed
// two-stage order (block partials summed on the host in block order): deterministic for a
// given n, o, v and batch size.

#include "../qm_gpu_mp2_backend.h"

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

/// Partial E_OS / E_SS sums of one batch: one pair of doubles per block.
__global__ void mp2PairEnergyKernel(const double* __restrict__ K, const double* __restrict__ eo,
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
struct Mp2Work {
    double *E = nullptr, *Co = nullptr, *Cv = nullptr, *eo = nullptr, *ev = nullptr;
    double *buf1 = nullptr, *buf2 = nullptr, *part = nullptr;
    bool own_E = false;
    cublasHandle_t h = nullptr;
    ~Mp2Work()
    {
        if (own_E && E) cudaFree(E);
        for (double* p : { Co, Cv, eo, ev, buf1, buf2, part })
            if (p) cudaFree(p);
        if (h) cublasDestroy(h);
    }
};

double msSince(std::chrono::steady_clock::time_point t0)
{
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}

}  // namespace

// Canonical MP2 pair energies on the device; returns 1 on success, 0 with a reason in `why`.
// Claude Generated (Oct 2026, WP11).
extern "C" int curcuma_cuda_qm_mp2_energy(const qmint::QMGpuMP2Input* in, qmint::QMGpuMP2Output* out, char* why,
                                          int whylen)
{
    if (!in || !out) { setWhy(why, whylen, "null argument"); return 0; }
    const int n = in->n, nfc = in->n_frozen, nocc = in->n_occ;
    const int o = nocc - nfc, v = n - nocc;
    *out = qmint::QMGpuMP2Output();
    if (o <= 0 || v <= 0) return 1;  // nothing to correlate
    if (!in->eri_device && !in->eri_host) { setWhy(why, whylen, "no AO tensor"); return 0; }
    const size_t n2 = (size_t)n * n, n3 = n2 * n, n4 = n3 * n;
    if (n3 > (size_t)2147483647) { setWhy(why, whylen, "basis too large for 32-bit cuBLAS dimensions"); return 0; }
    if (in->device >= 0 && cudaSetDevice(in->device) != cudaSuccess) { setWhy(why, whylen, "cudaSetDevice failed"); return 0; }

    Mp2Work w;
    auto fail = [&](const std::string& what) {
        cudaGetLastError();  // clear a sticky non-fatal error for the CPU fallback
        setWhy(why, whylen, what);
        return 0;
    };
#define MP2_CUDA(call, what) do { if ((call) != cudaSuccess) return fail(std::string(what) + ": " + cudaGetErrorString(cudaGetLastError())); } while (0)
#define MP2_CUBLAS(call, what) do { if ((call) != CUBLAS_STATUS_SUCCESS) return fail(what); } while (0)

    const auto t0 = std::chrono::steady_clock::now();
    if (cublasCreate(&w.h) != CUBLAS_STATUS_SUCCESS) return fail("cublasCreate failed");
    cublasSetMathMode(w.h, CUBLAS_DEFAULT_MATH);  // plain FP64, no emulation / reduced precision

    // Memory: (uploaded tensor) + C blocks + one batch of ib (n^3 + n^2 v) doubles.
    size_t free_b = 0, total_b = 0;
    MP2_CUDA(cudaMemGetInfo(&free_b, &total_b), "cudaMemGetInfo");
    double avail = (double)free_b;
    if (in->max_device_mb > 0.0) avail = std::min(avail, in->max_device_mb * 1.0e6);
    if (!in->eri_device) avail -= (double)n4 * 8.0;
    avail -= (double)(2 * n2 + 2 * n + 2 * kRedBlocks) * 8.0;
    avail *= 0.9;  // cuBLAS workspace and fragmentation
    const double per_i = (double)(n3 + n2 * v) * 8.0;
    if (in->batch_mb > 0.0) avail = std::min(avail, std::max(per_i, in->batch_mb * 1.0e6));  // at least one i
    const int batch = (int)std::min<double>(o, std::floor(avail / per_i));
    if (batch < 1)
        return fail("device memory too small (" + std::to_string((long long)(free_b / 1000000)) + " MB free, need "
                    + std::to_string((long long)(((in->eri_device ? 0.0 : (double)n4 * 8.0) + per_i) / 1.0e6)) + " MB)");

    // Uploads: tensor (if not resident), active occupied / virtual MO blocks column-major, energies.
    if (in->eri_device) {
        w.E = const_cast<double*>(in->eri_device);
    } else {
        MP2_CUDA(cudaMalloc(&w.E, n4 * sizeof(double)), "cudaMalloc tensor");
        w.own_E = true;
        MP2_CUDA(cudaMemcpy(w.E, in->eri_host, n4 * sizeof(double), cudaMemcpyHostToDevice), "tensor upload");
    }
    std::vector<double> co((size_t)n * o), cv((size_t)n * v), eo(o), ev(v);
    for (int i = 0; i < o; ++i) {
        eo[i] = in->eps[nfc + i];
        for (int mu = 0; mu < n; ++mu) co[(size_t)i * n + mu] = in->C[(size_t)mu * n + nfc + i];
    }
    for (int a = 0; a < v; ++a) {
        ev[a] = in->eps[nocc + a];
        for (int mu = 0; mu < n; ++mu) cv[(size_t)a * n + mu] = in->C[(size_t)mu * n + nocc + a];
    }
    MP2_CUDA(cudaMalloc(&w.Co, co.size() * sizeof(double)), "cudaMalloc Co");
    MP2_CUDA(cudaMalloc(&w.Cv, cv.size() * sizeof(double)), "cudaMalloc Cv");
    MP2_CUDA(cudaMalloc(&w.eo, eo.size() * sizeof(double)), "cudaMalloc eo");
    MP2_CUDA(cudaMalloc(&w.ev, ev.size() * sizeof(double)), "cudaMalloc ev");
    MP2_CUDA(cudaMalloc(&w.part, 2 * kRedBlocks * sizeof(double)), "cudaMalloc partials");
    MP2_CUDA(cudaMemcpy(w.Co, co.data(), co.size() * sizeof(double), cudaMemcpyHostToDevice), "upload Co");
    MP2_CUDA(cudaMemcpy(w.Cv, cv.data(), cv.size() * sizeof(double), cudaMemcpyHostToDevice), "upload Cv");
    MP2_CUDA(cudaMemcpy(w.eo, eo.data(), eo.size() * sizeof(double), cudaMemcpyHostToDevice), "upload eo");
    MP2_CUDA(cudaMemcpy(w.ev, ev.data(), ev.size() * sizeof(double), cudaMemcpyHostToDevice), "upload ev");
    MP2_CUDA(cudaMalloc(&w.buf1, (size_t)batch * n3 * sizeof(double)), "cudaMalloc batch buffer 1");
    MP2_CUDA(cudaMalloc(&w.buf2, (size_t)batch * n2 * v * sizeof(double)), "cudaMalloc batch buffer 2");
    MP2_CUDA(cudaDeviceSynchronize(), "upload sync");
    out->t_upload_ms = msSince(t0);

    const double one = 1.0, zero = 0.0;
    std::vector<double> part(2 * kRedBlocks);
    double e_os = 0.0, e_ss = 0.0;
    int nbatches = 0;
    for (int i0 = 0; i0 < o; i0 += batch) {
        const int ib = std::min(batch, o - i0);
        ++nbatches;
        const auto tt = std::chrono::steady_clock::now();
        // 1. X1 (n^3 x ib) = E (n^3 x n) * Co[:, i0:i0+ib]
        MP2_CUBLAS(cublasDgemm(w.h, CUBLAS_OP_N, CUBLAS_OP_N, (int)n3, ib, n, &one, w.E, (int)n3, w.Co + (size_t)i0 * n, n,
                               &zero, w.buf1, (int)n3), "step 1 DGEMM");
        // 2. W_i (n^2 x v) = X1_i (n^2 x n) * Cv, i = 0..ib-1
        MP2_CUBLAS(cublasDgemmStridedBatched(w.h, CUBLAS_OP_N, CUBLAS_OP_N, (int)n2, v, n, &one, w.buf1, (int)n2,
                                             (long long)n3, w.Cv, n, 0LL, &zero, w.buf2, (int)n2, (long long)n2 * v, ib),
                   "step 2 batched DGEMM");
        // 3. Z_ia (n x o) = W_ia (n x n) * Co, (i,a) = 0..ib*v-1; into buf1
        MP2_CUBLAS(cublasDgemmStridedBatched(w.h, CUBLAS_OP_N, CUBLAS_OP_N, n, o, n, &one, w.buf2, n, (long long)n2, w.Co, n,
                                             0LL, &zero, w.buf1, n, (long long)n * o, ib * v),
                   "step 3 batched DGEMM");
        // 4. K_ia (v x o) = Cv^T (v x n) * Z_ia (n x o); into buf2
        MP2_CUBLAS(cublasDgemmStridedBatched(w.h, CUBLAS_OP_T, CUBLAS_OP_N, v, o, n, &one, w.Cv, n, 0LL, w.buf1, n,
                                             (long long)n * o, &zero, w.buf2, v, (long long)v * o, ib * v),
                   "step 4 batched DGEMM");
        MP2_CUDA(cudaDeviceSynchronize(), "transformation");
        out->t_transform_ms += msSince(tt);

        const auto te = std::chrono::steady_clock::now();
        mp2PairEnergyKernel<<<kRedBlocks, kRedThreads>>>(w.buf2, w.eo, w.ev, i0, ib, o, v, w.part, w.part + kRedBlocks);
        MP2_CUDA(cudaGetLastError(), "pair-energy kernel launch");
        MP2_CUDA(cudaMemcpy(part.data(), w.part, part.size() * sizeof(double), cudaMemcpyDeviceToHost), "partials download");
        for (int k = 0; k < kRedBlocks; ++k) {
            e_os += part[k];
            e_ss += part[kRedBlocks + k];
        }
        out->t_energy_ms += msSince(te);
    }
#undef MP2_CUDA
#undef MP2_CUBLAS
    out->e_os = e_os;
    out->e_ss = e_ss;
    out->batch = batch;
    out->nbatches = nbatches;
    return 1;
}
