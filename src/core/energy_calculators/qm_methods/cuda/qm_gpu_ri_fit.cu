/*
 * < CUDA metric solve B = L^-1 T of the native RI-JK setup (WP10) >
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0
 */

// Claude Generated (Oct 2026, WP10 scaling). AI-generated, machine-tested.
//
// qmri::RIJK::setup (mode jk, Cholesky metric) turns the Coulomb-normalised 3-index tensor T
// (N x ncol, row-major on the host, ncol = n^2) into B = L^-1 T with J_n = L L^T. That solve costs
// N^2 ncol flops and was half of the GPU RI-JK wall time when it ran on the host. Here it runs as
// cublasDtrsm on column chunks of T: the column-major view of columns c0..c0+w is X = T(:, c0:c0+w)^T
// (w x N), and T := L^-1 T is X := X L^-T (side right, lower, transposed). Two streams, each with a
// pinned staging buffer: while the device solves chunk i, the host packs chunk i+1 and unpacks
// chunk i-1 (OpenMP memcpy of the strided rows). T is overwritten with B; nothing stays on the device.
// FP64 only. The result differs from the host solve only in the summation order (last bits).

#include <cublas_v2.h>
#include <cuda_runtime.h>

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

void setWhyFit(char* why, int whylen, const std::string& s)
{
    if (why && whylen > 0) std::snprintf(why, (size_t)whylen, "%s", s.c_str());
}

struct FitResources {
    cublasHandle_t h = nullptr;
    cudaStream_t s[2] = { nullptr, nullptr };
    double* dL = nullptr;
    double* dX[2] = { nullptr, nullptr };
    double* hX[2] = { nullptr, nullptr };  // pinned staging
    ~FitResources()
    {
        for (int k = 0; k < 2; ++k) {
            if (s[k]) cudaStreamSynchronize(s[k]);
            if (dX[k]) cudaFree(dX[k]);
            if (hX[k]) cudaFreeHost(hX[k]);
            if (s[k]) cudaStreamDestroy(s[k]);
        }
        if (dL) cudaFree(dL);
        if (h) cublasDestroy(h);
    }
};

// X (w x N, column-major, ld w) <-> columns c0..c0+w of the row-major host T (N x ncol)
void packChunk(const double* T, long ncol, int N, long c0, long w, double* X, int threads)
{
#pragma omp parallel for schedule(static) num_threads(threads)
    for (int P = 0; P < N; ++P) std::memcpy(X + (size_t)P * w, T + (size_t)P * ncol + c0, sizeof(double) * (size_t)w);
}
void unpackChunk(double* T, long ncol, int N, long c0, long w, const double* X, int threads)
{
#pragma omp parallel for schedule(static) num_threads(threads)
    for (int P = 0; P < N; ++P) std::memcpy(T + (size_t)P * ncol + c0, X + (size_t)P * w, sizeof(double) * (size_t)w);
}

} // namespace

// In place T := L^-1 T (T: N x ncol row-major host array, L: N x N lower triangular, COLUMN-major).
// Returns 1 on success; 0 with a reason. On failure *touched says whether part of T was already
// overwritten (then T is lost and the caller must rebuild it; otherwise it can solve on the host).
// Claude Generated (Oct 2026).
extern "C" int curcuma_cuda_qm_ri_lower_solve(double* T, int N, long ncol, const double* L, int device, int host_threads,
                                              int* touched, char* why, int whylen)
{
    if (touched) *touched = 0;
    if (!T || !L || N <= 0 || ncol <= 0) { setWhyFit(why, whylen, "empty RI tensor"); return 0; }
    if (device >= 0 && cudaSetDevice(device) != cudaSuccess) { setWhyFit(why, whylen, "cudaSetDevice failed"); return 0; }
    const int threads = std::max(1, host_threads);
    FitResources r;
    size_t free_b = 0, total_b = 0;
    if (cudaMemGetInfo(&free_b, &total_b) != cudaSuccess) { setWhyFit(why, whylen, "cudaMemGetInfo failed"); return 0; }
    const size_t l_b = (size_t)N * N * sizeof(double);
    // chunk: about 256 MB per buffer (two device and two pinned buffers), at least 32 columns
    long w = (long)std::min<size_t>((size_t)ncol, std::max<size_t>(32, ((size_t)256 << 20) / ((size_t)N * sizeof(double))));
    size_t x_b = (size_t)w * N * sizeof(double);
    if (l_b + 2 * x_b + ((size_t)64 << 20) > free_b) {
        char b[256];
        std::snprintf(b, sizeof b, "metric solve needs %.1f MB, %.1f MB free on the device", (l_b + 2 * x_b) / 1e6, free_b / 1e6);
        setWhyFit(why, whylen, b);
        return 0;
    }
    if (cublasCreate(&r.h) != CUBLAS_STATUS_SUCCESS) { setWhyFit(why, whylen, "cublasCreate failed"); return 0; }
    if (cudaMalloc(&r.dL, l_b) != cudaSuccess) { setWhyFit(why, whylen, "device allocation failed"); return 0; }
    for (int k = 0; k < 2; ++k) {
        if (cudaStreamCreateWithFlags(&r.s[k], cudaStreamNonBlocking) != cudaSuccess || cudaMalloc(&r.dX[k], x_b) != cudaSuccess
            || cudaMallocHost(&r.hX[k], x_b) != cudaSuccess) {
            setWhyFit(why, whylen, "device/pinned allocation failed");
            return 0;
        }
    }
    if (cudaMemcpy(r.dL, L, l_b, cudaMemcpyHostToDevice) != cudaSuccess) { setWhyFit(why, whylen, "metric upload failed"); return 0; }
    const double one = 1.0;
    const long nchunks = (ncol + w - 1) / w;
    long pending_c0[2] = { -1, -1 }, pending_w[2] = { 0, 0 };
    auto drain = [&](int k) -> bool {
        if (pending_c0[k] < 0) return true;
        if (cudaStreamSynchronize(r.s[k]) != cudaSuccess) return false;
        unpackChunk(T, ncol, N, pending_c0[k], pending_w[k], r.hX[k], threads);
        if (touched) *touched = 1;
        pending_c0[k] = -1;
        return true;
    };
    for (long b = 0; b < nchunks; ++b) {
        const int k = (int)(b & 1);
        if (!drain(k)) { setWhyFit(why, whylen, "device metric solve failed (stream)"); return 0; }
        const long c0 = b * w, wb = std::min(w, ncol - c0);
        packChunk(T, ncol, N, c0, wb, r.hX[k], threads);
        const size_t bytes = (size_t)wb * N * sizeof(double);
        if (cudaMemcpyAsync(r.dX[k], r.hX[k], bytes, cudaMemcpyHostToDevice, r.s[k]) != cudaSuccess) {
            setWhyFit(why, whylen, "chunk upload failed");
            return 0;
        }
        if (cublasSetStream(r.h, r.s[k]) != CUBLAS_STATUS_SUCCESS
            || cublasDtrsm(r.h, CUBLAS_SIDE_RIGHT, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_T, CUBLAS_DIAG_NON_UNIT, (int)wb, N, &one, r.dL, N,
                           r.dX[k], (int)wb)
                != CUBLAS_STATUS_SUCCESS) {
            setWhyFit(why, whylen, "cublasDtrsm failed");
            return 0;
        }
        if (cudaMemcpyAsync(r.hX[k], r.dX[k], bytes, cudaMemcpyDeviceToHost, r.s[k]) != cudaSuccess) {
            setWhyFit(why, whylen, "chunk download failed");
            return 0;
        }
        pending_c0[k] = c0;
        pending_w[k] = wb;
    }
    if (!drain(0) || !drain(1)) { setWhyFit(why, whylen, "device metric solve failed (stream)"); return 0; }
    if (cudaGetLastError() != cudaSuccess) { setWhyFit(why, whylen, "CUDA error after the metric solve"); return 0; }
    return 1;
}
