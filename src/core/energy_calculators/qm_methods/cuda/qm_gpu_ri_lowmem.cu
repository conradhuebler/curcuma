/*
 * < CUDA resource-saving density fitting of the native QM engine (WP10): streamed / FP32 tensor, batched exchange >
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0
 */

// Claude Generated (Oct 2026, WP10 "RI-JK ressourcensparend"). AI-generated, machine-tested.
//
// 1. QmGpuRiLowMem (entry curcuma_cuda_create_qm_ri_lowmem): the fitted tensor T (nfit x np, np = n(n+1)/2,
//    packed rows exactly as qm_gpu_ri.cu) in FP64 or FP32, either RESIDENT on the device or STREAMED from host
//    memory (page-locked in place with cudaHostRegister) in row chunks. Streaming uses one copy stream and one
//    compute stream with two device buffers: the copy of chunk i+1 runs while chunk i is contracted (events order
//    the reuse of a buffer). Per chunk of rows (column-major view Mp = [T_q0 .. T_q0+cnt-1], np x cnt):
//      g_q = T_q . w             FP64: cublasDgemv(T)      FP32 tensor: kernel, FP64 accumulation
//      J  += sum_q g_q T_q       FP64: cublasDgemv(N)      FP32 tensor: kernel, FP64 accumulation
//      K  += sum_q s (T_q V)(T_q V)^T: unpack nq rows (n x nq n), Y = U^T V (GEMM), SYRK
//                                FP64: DGEMM + DSYRK into K  FP32: SGEMM + SSYRK into an FP32 partial K per
//                                sub-batch, added into the FP64 K (kernel), so the FP32 rounding does not grow with nfit
//    mode jk: buildJK does g, J and K in ONE pass (c = g). contract/coulomb/exchange (QMGpuRI) are one pass each.
// 2. QmGpuRiYSolver (entry curcuma_cuda_create_qm_ri_ysolver): batched storage (no tensor). L (Cholesky factor of
//    the normalised metric) stays on the device; per call the host's Y (naux x n rb, row P = column-major n x rb
//    block T^P V) is uploaded in orbital column chunks, X = L^-1 Y (cublasDtrsm on the aux index, columns are
//    independent) and K += sign X X^T (one DSYRK: the chunk read as n x (kc naux) column-major holds the vectors X^Q_k).

#include "../qm_gpu_ri_lowmem_backend.h"

#include <cublas_v2.h>
#include <cuda_runtime.h>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

namespace {

double msNow(std::chrono::steady_clock::time_point t0)
{
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}

template <typename T, typename TO>
__global__ void unpackPairsT(const T* __restrict__ Tp, TO* __restrict__ U, int n, long long np, int nq)
{
    const long long n2 = (long long)n * n, total = (long long)nq * n2;
    for (long long idx = (long long)blockIdx.x * blockDim.x + threadIdx.x; idx < total; idx += (long long)gridDim.x * blockDim.x) {
        const long long q = idx / n2, rem = idx - q * n2;
        const int i = (int)(rem / n), j = (int)(rem - (long long)i * n);
        const int a = i > j ? i : j, b = i > j ? j : i;
        U[idx] = (TO)Tp[q * np + (long long)a * (a + 1) / 2 + b];
    }
}

/// g[q] = sum_p T[q][p] w[p] for FP32 rows, FP64 accumulation; one block per row.
__global__ void gRowsF32(const float* __restrict__ T, const double* __restrict__ w, double* __restrict__ g, long long np)
{
    __shared__ double red[256];
    const float* t = T + (long long)blockIdx.x * np;
    double s = 0.0;
    for (long long p = threadIdx.x; p < np; p += blockDim.x) s += (double)t[p] * w[p];
    red[threadIdx.x] = s;
    __syncthreads();
    for (int h = blockDim.x / 2; h > 0; h >>= 1) {
        if ((int)threadIdx.x < h) red[threadIdx.x] += red[threadIdx.x + h];
        __syncthreads();
    }
    if (threadIdx.x == 0) g[blockIdx.x] = red[0];
}

/// J[p] += sum_q c[q] T[q][p] for FP32 rows, FP64 accumulation (q in order per column).
__global__ void jColsF32(const float* __restrict__ T, const double* __restrict__ c, double* __restrict__ J, long long np, int cnt)
{
    for (long long p = (long long)blockIdx.x * blockDim.x + threadIdx.x; p < np; p += (long long)gridDim.x * blockDim.x) {
        double s = J[p];
        for (int q = 0; q < cnt; ++q) s += c[q] * (double)T[(long long)q * np + p];
        J[p] = s;
    }
}

__global__ void toFloat(const double* __restrict__ a, float* __restrict__ b, long long m)
{
    for (long long i = (long long)blockIdx.x * blockDim.x + threadIdx.x; i < m; i += (long long)gridDim.x * blockDim.x) b[i] = (float)a[i];
}

/// Kd(lower) += Kf(lower), column-major n x n
__global__ void addLowerF32(const float* __restrict__ Kf, double* __restrict__ Kd, int n)
{
    const long long n2 = (long long)n * n;
    for (long long idx = (long long)blockIdx.x * blockDim.x + threadIdx.x; idx < n2; idx += (long long)gridDim.x * blockDim.x) {
        const int j = (int)(idx / n), i = (int)(idx - (long long)j * n);
        if (i >= j) Kd[idx] += (double)Kf[idx];
    }
}

int gridFor(long long total, int threads = 256) { return (int)std::max<long long>(1, std::min<long long>((total + threads - 1) / threads, 65535LL * 8)); }

class QmGpuRiLowMem : public qmint::QMGpuRILowMem {
public:
    ~QmGpuRiLowMem() override
    {
        if (m_sc) cudaStreamSynchronize(m_sc);
        if (m_sx) cudaStreamSynchronize(m_sx);
        for (cudaEvent_t e : m_events)
            if (e) cudaEventDestroy(e);
        for (void* p : m_owned) cudaFree(p);
        if (m_registered) cudaHostUnregister(const_cast<void*>(m_host));
        if (m_cublas) cublasDestroy(m_cublas);
        if (m_sc) cudaStreamDestroy(m_sc);
        if (m_sx) cudaStreamDestroy(m_sx);
    }

    bool init(const void* T, bool fp32, int nfit, int n, int device, int stream, double max_dev_mb, double chunk_mb, std::string& err)
    {
        m_host = T;
        m_fp32 = fp32;
        m_es = fp32 ? sizeof(float) : sizeof(double);
        m_nfit = nfit;
        m_n = n;
        m_np = (long long)n * (n + 1) / 2;
        m_device = device;
        const size_t n2 = (size_t)n * n, np = (size_t)m_np;
        size_t free_b = 0, total_b = 0;
        if (cudaMemGetInfo(&free_b, &total_b) != cudaSuccess) { err = "cudaMemGetInfo failed"; return false; }
        if (max_dev_mb > 0.0) free_b = std::min(free_b, (size_t)(max_dev_mb * 1.0e6));
        const size_t slack = (size_t)64 << 20;
        // small arrays: w, J, g (FP64), K (FP64), V (FP64 + FP32), FP32 partial K
        const size_t small_b = (np * 2 + (size_t)nfit) * sizeof(double) + n2 * (2 * sizeof(double) + 2 * sizeof(float));
        const size_t blk_b = n2 * sizeof(double);  // one unpacked block; FP64 also for an FP32 tensor (second phase)
        const size_t tensor_b = np * (size_t)nfit * m_es;
        // K work buffers (unpacked rows U and Y): each 512 MB or less, at least one block
        auto workFor = [&](size_t left) { return std::max(std::min<size_t>((size_t)512 << 20, left / 5), blk_b); };
        bool resident = false;
        if (stream != 1 && tensor_b + small_b + slack + 2 * blk_b <= free_b) {
            m_work_b = workFor(free_b - tensor_b - small_b - slack);
            resident = tensor_b + small_b + slack + 2 * m_work_b <= free_b;
        }
        if (!resident && stream == 2) {
            char b[256];
            std::snprintf(b, sizeof b, "RI tensor %.1f MB (%s) does not fit into %.1f MB device memory and streaming is off", tensor_b / 1e6,
                          fp32 ? "FP32" : "FP64", free_b / 1e6);
            err = b;
            return false;
        }
        if (cublasCreate(&m_cublas) != CUBLAS_STATUS_SUCCESS) { err = "cublasCreate failed"; return false; }
        cublasSetMathMode(m_cublas, CUBLAS_DEFAULT_MATH);  // no TF32: FP32 inputs stay FP32
        if (cudaStreamCreateWithFlags(&m_sc, cudaStreamNonBlocking) != cudaSuccess
            || cudaStreamCreateWithFlags(&m_sx, cudaStreamNonBlocking) != cudaSuccess) {
            err = "cudaStreamCreate failed";
            return false;
        }
        cublasSetStream(m_cublas, m_sc);
        if (resident) {
            m_rows_chunk = nfit;
            if (!alloc(&m_T, tensor_b)) { err = "device allocation of the tensor failed"; return false; }
            if (cudaMemcpy(m_T, T, tensor_b, cudaMemcpyHostToDevice) != cudaSuccess || cudaDeviceSynchronize() != cudaSuccess) { err = "tensor upload failed"; return false; }
        } else {
            m_streamed = true;
            if (small_b + slack + 2 * blk_b + 2 * np * m_es > free_b) {
                char b[256];
                std::snprintf(b, sizeof b, "streamed RI needs at least %.1f MB device memory (two rows, one block), %.1f MB free",
                              (small_b + slack + 2 * blk_b + 2 * np * m_es) / 1e6, free_b / 1e6);
                err = b;
                return false;
            }
            const size_t left = free_b - small_b - slack;
            size_t chunk_b = chunk_mb > 0.0 ? (size_t)(chunk_mb * 1.0e6) : std::min<size_t>((size_t)512 << 20, left / 4);
            chunk_b = std::max(chunk_b, np * m_es);
            m_work_b = std::max(std::min<size_t>((size_t)512 << 20, (left - 2 * std::min(chunk_b, left / 3)) / 3), blk_b);
            if (2 * chunk_b + 2 * m_work_b > left) chunk_b = std::max(np * m_es, (left - 2 * m_work_b) / 2);
            m_rows_chunk = (int)std::max<size_t>(1, std::min<size_t>((size_t)nfit, chunk_b / (np * m_es)));
            for (int k = 0; k < 2; ++k)
                if (!alloc(&m_buf[k], (size_t)m_rows_chunk * np * m_es)) { err = "device allocation of the stream buffers failed"; return false; }
            // page-lock the host tensor in place (async copies need it); without it the copies still work, only slower
            m_registered = cudaHostRegister(const_cast<void*>(T), tensor_b, cudaHostRegisterDefault) == cudaSuccess;
            if (!m_registered) cudaGetLastError();
            for (int k = 0; k < 4; ++k) {
                cudaEvent_t e = nullptr;
                if (cudaEventCreateWithFlags(&e, cudaEventDisableTiming) != cudaSuccess) { err = "cudaEventCreate failed"; return false; }
                m_events.push_back(e);
            }
        }
        if (!alloc(&m_U, m_work_b) || !alloc(&m_X, m_work_b) || !alloc(&m_w, np * sizeof(double)) || !alloc(&m_J, np * sizeof(double))
            || !alloc(&m_g, (size_t)nfit * sizeof(double)) || !alloc(&m_K, n2 * sizeof(double)) || !alloc(&m_V, n2 * sizeof(double))
            || !alloc(&m_Vf, n2 * sizeof(float)) || !alloc(&m_Kf, n2 * sizeof(float))) {
            err = "device allocation failed";
            return false;
        }
        char b[192];
        std::snprintf(b, sizeof b, "cuda:%d (cuBLAS %s, %s)", device, fp32 ? "FP32 tensor, FP64 accumulation" : "FP64",
                      m_streamed ? (m_registered ? "streamed from pinned host memory" : "streamed from pageable host memory") : "resident");
        m_name = b;
        return true;
    }

    std::string name() const override { return m_name; }
    bool streamed() const override { return m_streamed; }
    bool fp32() const override { return m_fp32; }
    void setFP64Arithmetic(bool on) override { m_fp64arith = on; }
    std::string stats() const override
    {
        char b[256];
        std::snprintf(b, sizeof b, "%ld device passes, %.2f GB streamed, %.1f ms in device J/K calls, %d rows per chunk", m_passes,
                      m_streamed_b / 1e9, m_ms, m_rows_chunk);
        return b;
    }

    bool contract(const double* w, double* g) override
    {
        if (cudaMemcpyAsync(m_w, w, sizeof(double) * m_np, cudaMemcpyHostToDevice, m_sc) != cudaSuccess) return false;
        if (!pass(true, false, false, 0, 0)) return false;
        return download(g, m_g, sizeof(double) * m_nfit);
    }

    bool coulomb(const double* c, double* J) override
    {
        if (cudaMemcpyAsync(m_g, c, sizeof(double) * m_nfit, cudaMemcpyHostToDevice, m_sc) != cudaSuccess) return false;
        if (cudaMemsetAsync(m_J, 0, sizeof(double) * m_np, m_sc) != cudaSuccess) return false;
        if (!pass(false, true, false, 0, 0)) return false;
        return download(J, m_J, sizeof(double) * m_np);
    }

    bool exchange(const double* V, int r, int npos, double* K) override
    {
        if (!uploadV(V, r, npos)) return false;
        if (cudaMemsetAsync(m_K, 0, sizeof(double) * (size_t)m_n * m_n, m_sc) != cudaSuccess) return false;
        if (!pass(false, false, true, r, npos)) return false;
        return downloadK(K);
    }

    bool buildJK(const double* w, const double* V, int r, int npos, double* g, double* Jp, double* K) override
    {
        if (cudaMemcpyAsync(m_w, w, sizeof(double) * m_np, cudaMemcpyHostToDevice, m_sc) != cudaSuccess) return false;
        if (!uploadV(V, r, npos)) return false;
        if (cudaMemsetAsync(m_J, 0, sizeof(double) * m_np, m_sc) != cudaSuccess) return false;
        if (cudaMemsetAsync(m_K, 0, sizeof(double) * (size_t)m_n * m_n, m_sc) != cudaSuccess) return false;
        if (!pass(true, true, true, r, npos)) return false;
        if (!download(g, m_g, sizeof(double) * m_nfit)) return false;
        if (!download(Jp, m_J, sizeof(double) * m_np)) return false;
        return downloadK(K);
    }

private:
    bool alloc(void* pp, size_t bytes)
    {
        void* q = nullptr;
        if (cudaMalloc(&q, std::max<size_t>(1, bytes)) != cudaSuccess) return false;
        m_owned.push_back(q);
        *static_cast<void**>(pp) = q;
        return true;
    }

    /// device -> host on the compute stream, then wait (all uploads/kernels run on m_sc, so the order holds)
    bool download(void* dst, const void* src, size_t bytes)
    {
        return cudaMemcpyAsync(dst, src, bytes, cudaMemcpyDeviceToHost, m_sc) == cudaSuccess && cudaStreamSynchronize(m_sc) == cudaSuccess;
    }

    bool uploadV(const double* V, int r, int npos)
    {
        if (r > m_n || r < 0 || npos < 0 || npos > r) return false;
        if (cudaMemcpyAsync(m_V, V, sizeof(double) * (size_t)m_n * r, cudaMemcpyHostToDevice, m_sc) != cudaSuccess) return false;
        if (m_fp32) {
            const long long m = (long long)m_n * r;
            toFloat<<<gridFor(m), 256, 0, m_sc>>>(static_cast<const double*>(m_V), static_cast<float*>(m_Vf), m);
            if (cudaGetLastError() != cudaSuccess) return false;
        }
        return true;
    }

    bool downloadK(double* K)
    {
        const int n = m_n;
        const size_t n2 = (size_t)n * n;
        std::vector<double> Kc(n2);
        if (!download(Kc.data(), m_K, sizeof(double) * n2)) return false;
        for (int j = 0; j < n; ++j)
            for (int i = j; i < n; ++i) {
                const double v = Kc[(size_t)j * n + i];
                K[(size_t)i * n + j] = v;
                K[(size_t)j * n + i] = v;
            }
        return true;
    }

    /// Contractions of one chunk of rows [q0, q0 + cnt) at device pointer Tc, on the compute stream.
    bool chunk(const void* Tc, int q0, int cnt, bool doG, bool doJ, bool doK, int r, int npos)
    {
        const int n = m_n;
        const long long np = m_np;
        const double one = 1.0, zero = 0.0;
        double* g = static_cast<double*>(m_g) + q0;
        if (doG) {
            if (m_fp32)
                gRowsF32<<<cnt, 256, 0, m_sc>>>(static_cast<const float*>(Tc), static_cast<const double*>(m_w), g, np);
            else if (cublasDgemv(m_cublas, CUBLAS_OP_T, (int)np, cnt, &one, static_cast<const double*>(Tc), (int)np, static_cast<const double*>(m_w),
                                 1, &zero, g, 1)
                     != CUBLAS_STATUS_SUCCESS)
                return false;
        }
        if (doJ) {
            if (m_fp32)
                jColsF32<<<gridFor(np), 256, 0, m_sc>>>(static_cast<const float*>(Tc), g, static_cast<double*>(m_J), np, cnt);
            else if (cublasDgemv(m_cublas, CUBLAS_OP_N, (int)np, cnt, &one, static_cast<const double*>(Tc), (int)np, g, 1, &one,
                                 static_cast<double*>(m_J), 1)
                     != CUBLAS_STATUS_SUCCESS)
                return false;
        }
        if (cudaGetLastError() != cudaSuccess) return false;
        if (!doK) return true;
        const size_t n2 = (size_t)n * n;
        const int maxcols = std::max(1, std::max(npos, r - npos));
        const bool f32k = m_fp32 && !m_fp64arith;  // FP32 GEMM/SYRK
        const size_t ea = f32k ? sizeof(float) : sizeof(double);
        const int sub = (int)std::max<size_t>(1, std::min({ (size_t)cnt, m_work_b / ea / ((size_t)n * maxcols), m_work_b / ea / n2 }));
        for (int s0 = 0; s0 < cnt; s0 += sub) {
            const int nq = std::min(sub, cnt - s0);
            const long long total = (long long)nq * (long long)n2;
            if (f32k)
                unpackPairsT<float, float><<<gridFor(total), 256, 0, m_sc>>>(static_cast<const float*>(Tc) + (size_t)s0 * np, static_cast<float*>(m_U), n, np, nq);
            else if (m_fp32)
                unpackPairsT<float, double><<<gridFor(total), 256, 0, m_sc>>>(static_cast<const float*>(Tc) + (size_t)s0 * np, static_cast<double*>(m_U), n, np, nq);
            else
                unpackPairsT<double, double><<<gridFor(total), 256, 0, m_sc>>>(static_cast<const double*>(Tc) + (size_t)s0 * np, static_cast<double*>(m_U), n, np, nq);
            if (cudaGetLastError() != cudaSuccess) return false;
            for (int ps = 0; ps < 2; ++ps) {
                const int cols = ps == 0 ? npos : r - npos;
                if (cols == 0) continue;
                const int off = ps == 0 ? 0 : npos;
                if (f32k) {
                    const float onef = 1.0f, zerof = 0.0f, alpha = ps == 0 ? 1.0f : -1.0f;
                    if (cublasSgemm(m_cublas, CUBLAS_OP_T, CUBLAS_OP_N, nq * n, cols, n, &onef, static_cast<const float*>(m_U), n,
                                    static_cast<const float*>(m_Vf) + (size_t)off * n, n, &zerof, static_cast<float*>(m_X), nq * n)
                            != CUBLAS_STATUS_SUCCESS
                        || cublasSsyrk(m_cublas, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_N, n, nq * cols, &alpha, static_cast<const float*>(m_X), n,
                                       &zerof, static_cast<float*>(m_Kf), n)
                            != CUBLAS_STATUS_SUCCESS)
                        return false;
                    addLowerF32<<<gridFor((long long)n2), 256, 0, m_sc>>>(static_cast<const float*>(m_Kf), static_cast<double*>(m_K), n);
                    if (cudaGetLastError() != cudaSuccess) return false;
                } else {
                    const double alpha = ps == 0 ? 1.0 : -1.0;
                    if (cublasDgemm(m_cublas, CUBLAS_OP_T, CUBLAS_OP_N, nq * n, cols, n, &one, static_cast<const double*>(m_U), n,
                                    static_cast<const double*>(m_V) + (size_t)off * n, n, &zero, static_cast<double*>(m_X), nq * n)
                            != CUBLAS_STATUS_SUCCESS
                        || cublasDsyrk(m_cublas, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_N, n, nq * cols, &alpha, static_cast<const double*>(m_X), n,
                                       &one, static_cast<double*>(m_K), n)
                            != CUBLAS_STATUS_SUCCESS)
                        return false;
                }
            }
        }
        return true;
    }

    /// One pass over the tensor (resident: chunk = whole tensor; streamed: copy stream + compute stream, two buffers).
    bool pass(bool doG, bool doJ, bool doK, int r, int npos)
    {
        const auto t0 = std::chrono::steady_clock::now();
        ++m_passes;
        bool ok = true;
        if (!m_streamed) {
            ok = chunk(m_T, 0, m_nfit, doG, doJ, doK, r, npos);
        } else {
            const size_t row_b = (size_t)m_np * m_es;
            const int nchunks = (m_nfit + m_rows_chunk - 1) / m_rows_chunk;
            cudaEvent_t copied[2] = { m_events[0], m_events[1] }, used[2] = { m_events[2], m_events[3] };
            auto issueCopy = [&](int b) -> bool {
                const int k = b & 1, q0 = b * m_rows_chunk, cnt = std::min(m_rows_chunk, m_nfit - q0);
                if (b >= 2 && cudaStreamWaitEvent(m_sx, used[k], 0) != cudaSuccess) return false;  // buffer k free again
                const char* src = static_cast<const char*>(m_host) + (size_t)q0 * row_b;
                if (cudaMemcpyAsync(m_buf[k], src, (size_t)cnt * row_b, cudaMemcpyHostToDevice, m_sx) != cudaSuccess) return false;
                m_streamed_b += (double)cnt * row_b;
                return cudaEventRecord(copied[k], m_sx) == cudaSuccess;
            };
            ok = issueCopy(0);
            for (int b = 0; ok && b < nchunks; ++b) {
                const int k = b & 1, q0 = b * m_rows_chunk, cnt = std::min(m_rows_chunk, m_nfit - q0);
                if (b + 1 < nchunks) ok = issueCopy(b + 1);
                ok = ok && cudaStreamWaitEvent(m_sc, copied[k], 0) == cudaSuccess;
                ok = ok && chunk(m_buf[k], q0, cnt, doG, doJ, doK, r, npos);
                ok = ok && cudaEventRecord(used[k], m_sc) == cudaSuccess;
            }
        }
        ok = ok && cudaStreamSynchronize(m_sc) == cudaSuccess && cudaStreamSynchronize(m_sx) == cudaSuccess && cudaGetLastError() == cudaSuccess;
        m_ms += msNow(t0);
        return ok;
    }

    const void* m_host = nullptr;
    bool m_fp32 = false, m_streamed = false, m_registered = false, m_fp64arith = false;
    size_t m_es = 8, m_work_b = 0;
    int m_nfit = 0, m_n = 0, m_device = 0, m_rows_chunk = 0;
    long long m_np = 0;
    long m_passes = 0;
    double m_streamed_b = 0.0, m_ms = 0.0;
    cublasHandle_t m_cublas = nullptr;
    cudaStream_t m_sc = nullptr, m_sx = nullptr;
    std::vector<cudaEvent_t> m_events;
    std::vector<void*> m_owned;
    void *m_T = nullptr, *m_buf[2] = { nullptr, nullptr }, *m_U = nullptr, *m_X = nullptr, *m_w = nullptr, *m_J = nullptr, *m_g = nullptr,
         *m_K = nullptr, *m_V = nullptr, *m_Vf = nullptr, *m_Kf = nullptr;
    std::string m_name = "cuda";
};

class QmGpuRiYSolverCuda : public qmint::QMGpuRIYSolver {
public:
    ~QmGpuRiYSolverCuda() override
    {
        if (m_dL) cudaFree(m_dL);
        if (m_dZ) cudaFree(m_dZ);
        if (m_dK) cudaFree(m_dK);
        if (m_h) cublasDestroy(m_h);
    }
    bool init(const double* L, int naux, int n, int device, double max_dev_mb, std::string& err)
    {
        m_naux = naux;
        m_n = n;
        size_t free_b = 0, total_b = 0;
        if (cudaMemGetInfo(&free_b, &total_b) != cudaSuccess) { err = "cudaMemGetInfo failed"; return false; }
        if (max_dev_mb > 0.0) free_b = std::min(free_b, (size_t)(max_dev_mb * 1.0e6));
        const size_t l_b = (size_t)naux * naux * sizeof(double), k_b = (size_t)n * n * sizeof(double), col_b = (size_t)n * naux * sizeof(double);
        const size_t slack = (size_t)64 << 20;
        if (l_b + k_b + col_b + slack > free_b) {
            char b[256];
            std::snprintf(b, sizeof b, "batched RI exchange needs at least %.1f MB device memory, %.1f MB free", (l_b + k_b + col_b + slack) / 1e6,
                          free_b / 1e6);
            err = b;
            return false;
        }
        m_kc_max = (int)std::max<size_t>(1, (free_b - l_b - k_b - slack) / col_b);
        m_kc_max = std::min(m_kc_max, std::max(1, (int)(((size_t)4 << 30) / col_b)));  // at most 4 GB per chunk
        if (cublasCreate(&m_h) != CUBLAS_STATUS_SUCCESS) { err = "cublasCreate failed"; return false; }
        if (cudaMalloc(&m_dL, l_b) != cudaSuccess || cudaMalloc(&m_dK, k_b) != cudaSuccess
            || cudaMalloc(&m_dZ, (size_t)m_kc_max * col_b) != cudaSuccess) {
            err = "device allocation failed";
            return false;
        }
        if (cudaMemcpy(m_dL, L, l_b, cudaMemcpyHostToDevice) != cudaSuccess) { err = "metric factor upload failed"; return false; }
        char b[128];
        std::snprintf(b, sizeof b, "cuda:%d (cuBLAS FP64, up to %d orbitals per chunk)", device, m_kc_max);
        m_name = b;
        return true;
    }
    std::string name() const override { return m_name; }
    std::string stats() const override
    {
        char b[160];
        std::snprintf(b, sizeof b, "%ld calls, %.2f GB uploaded, %.1f ms", m_calls, m_up_b / 1e9, m_ms);
        return b;
    }
    bool exchange(const double* Y, int rb, double sign, double* K) override
    {
        const auto t0 = std::chrono::steady_clock::now();
        ++m_calls;
        const int n = m_n, N = m_naux;
        const double one = 1.0;
        if (cudaMemset(m_dK, 0, sizeof(double) * (size_t)n * n) != cudaSuccess) return false;
        for (int k0 = 0; k0 < rb; k0 += m_kc_max) {
            const int kc = std::min(m_kc_max, rb - k0);
            const size_t w = (size_t)kc * n;  // columns of Y in this chunk, contiguous per row P
            // Z (w x N, column-major): column P = Y(P, k0 n .. k0 n + w)
            if (cudaMemcpy2D(m_dZ, w * sizeof(double), Y + (size_t)k0 * n, (size_t)n * rb * sizeof(double), w * sizeof(double), (size_t)N,
                             cudaMemcpyHostToDevice)
                != cudaSuccess)
                return false;
            m_up_b += (double)w * N * sizeof(double);
            // Z := Z L^-T  <=>  X = L^-1 Y on the aux index
            if (cublasDtrsm(m_h, CUBLAS_SIDE_RIGHT, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_T, CUBLAS_DIAG_NON_UNIT, (int)w, N, &one, m_dL, N, m_dZ,
                            (int)w)
                != CUBLAS_STATUS_SUCCESS)
                return false;
            // the chunk read as n x (kc N) column-major: column (Q kc + k) = X^Q_{:, k0 + k}
            if (cublasDsyrk(m_h, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_N, n, kc * N, &sign, m_dZ, n, &one, m_dK, n) != CUBLAS_STATUS_SUCCESS)
                return false;
        }
        std::vector<double> Kc((size_t)n * n);
        if (cudaMemcpy(Kc.data(), m_dK, sizeof(double) * Kc.size(), cudaMemcpyDeviceToHost) != cudaSuccess) return false;
        for (int j = 0; j < n; ++j)
            for (int i = j; i < n; ++i) {
                const double v = Kc[(size_t)j * n + i];
                K[(size_t)i * n + j] += v;
                if (i != j) K[(size_t)j * n + i] += v;
            }
        m_ms += msNow(t0);
        return true;
    }

private:
    int m_naux = 0, m_n = 0, m_kc_max = 1;
    long m_calls = 0;
    double m_up_b = 0.0, m_ms = 0.0;
    double *m_dL = nullptr, *m_dZ = nullptr, *m_dK = nullptr;
    cublasHandle_t m_h = nullptr;
    std::string m_name = "cuda";
};

void setWhyLm(char* why, int whylen, const std::string& s)
{
    if (why && whylen > 0) std::snprintf(why, (size_t)whylen, "%s", s.c_str());
}

} // namespace

// Resident or streamed RI tensor (FP64 or FP32), see the header comment. Claude Generated (Oct 2026).
extern "C" void* curcuma_cuda_create_qm_ri_lowmem(const void* T, int fp32, int nfit, int n, int device, int stream, double max_dev_mb,
                                                  double chunk_mb, int host_threads, char* why, int whylen)
{
    (void)host_threads;
    if (!T || nfit <= 0 || n <= 0) { setWhyLm(why, whylen, "empty RI tensor"); return nullptr; }
    if (device >= 0 && cudaSetDevice(device) != cudaSuccess) { setWhyLm(why, whylen, "cudaSetDevice failed"); return nullptr; }
    auto s = std::make_unique<QmGpuRiLowMem>();
    std::string err;
    if (!s->init(T, fp32 != 0, nfit, n, device, stream, max_dev_mb, chunk_mb, err)) {
        cudaGetLastError();
        setWhyLm(why, whylen, err);
        return nullptr;
    }
    return static_cast<qmint::QMGpuRILowMem*>(s.release());
}

// Batched RI-JK exchange X = L^-1 Y, K += sign X X^T. Claude Generated (Oct 2026).
extern "C" void* curcuma_cuda_create_qm_ri_ysolver(const double* L, int naux, int n, int device, double max_dev_mb, int host_threads,
                                                   char* why, int whylen)
{
    (void)host_threads;
    if (!L || naux <= 0 || n <= 0) { setWhyLm(why, whylen, "empty metric factor"); return nullptr; }
    if (device >= 0 && cudaSetDevice(device) != cudaSuccess) { setWhyLm(why, whylen, "cudaSetDevice failed"); return nullptr; }
    auto s = std::make_unique<QmGpuRiYSolverCuda>();
    std::string err;
    if (!s->init(L, naux, n, device, max_dev_mb, err)) {
        cudaGetLastError();
        setWhyLm(why, whylen, err);
        return nullptr;
    }
    return static_cast<qmint::QMGpuRIYSolver*>(s.release());
}
