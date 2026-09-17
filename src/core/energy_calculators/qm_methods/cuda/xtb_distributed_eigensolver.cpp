/*
 * < Multi-GPU dense symmetric eigensolver for the native xTB SCF - implementation >
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 */

// Claude Generated (Sep 2026, multi-GPU step 3). See xtb_distributed_eigensolver.h.
//
// Built into its own library, libcurcuma_cuda_mgpu.so, so that a missing cuSOLVERMp/NCCL on the
// run host only disables the multi-GPU eigensolve instead of the whole CUDA plugin.

#ifdef USE_CUDA

#include "xtb_distributed_eigensolver.h"

#include <cuda_runtime.h>

#ifdef CURCUMA_HAVE_CUSOLVERMG
#define DISABLE_CUSOLVERMG_DEPRECATED
#include <cusolverMg.h>
#endif
#ifdef CURCUMA_HAVE_CUSOLVERMP
#include <cublasmp.h>
#include <cusolverMp.h>
#include <nccl.h>
#endif

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstring>
#include <thread>

namespace curcuma {
namespace xtb {
namespace gpu {

namespace {

// One column block of the global matrix and where it lives.
struct ColumnBlock {
    int rank;          // index into the device list
    long col0;         // first global column
    long ncol;         // columns in the block
    long local_col0;   // first column inside the rank's local array
};

// Cyclic 1-D column distribution (block b -> rank b % ndev), the layout of a 1 x ndev
// ScaLAPACK/cusolverMg grid. local_cols[r] = columns held by rank r.
std::vector<ColumnBlock> columnLayout(int n, int block, int ndev, std::vector<long>& local_cols)
{
    std::vector<ColumnBlock> blocks;
    local_cols.assign(ndev, 0);
    const int nblocks = (n + block - 1) / block;
    for (int b = 0; b < nblocks; ++b) {
        const int r = b % ndev;
        const long c0 = static_cast<long>(b) * block;
        const long nc = std::min<long>(block, n - c0);
        blocks.push_back({ r, c0, nc, local_cols[r] });
        local_cols[r] += nc;
    }
    return blocks;
}

double msSince(std::chrono::steady_clock::time_point t0)
{
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}

void enablePeers(const std::vector<int>& devs)
{
    for (int a : devs) {
        cudaSetDevice(a);
        for (int b : devs) {
            int can = 0;
            if (a != b && cudaDeviceCanAccessPeer(&can, a, b) == cudaSuccess && can)
                cudaDeviceEnablePeerAccess(b, 0);   // "already enabled" is harmless
        }
    }
}

// Run fn(rank) on one host thread per rank and report whether every rank succeeded.
template <class Fn>
bool onAllRanks(int ndev, Fn&& fn)
{
    std::vector<char> ok(ndev, 0);
    std::vector<std::thread> threads;
    threads.reserve(ndev);
    for (int r = 0; r < ndev; ++r)
        threads.emplace_back([&, r]() { ok[r] = fn(r) ? 1 : 0; });
    for (auto& t : threads) t.join();
    return std::all_of(ok.begin(), ok.end(), [](char c) { return c != 0; });
}

#ifdef CURCUMA_HAVE_CUSOLVERMP
// ---------------------------------------------------------------------------------------------
// cuSOLVERMp: one rank per GPU, collective calls on one host thread each.
// ---------------------------------------------------------------------------------------------
class MpEigensolver final : public DistributedEigensolver {
public:
    MpEigensolver(const std::vector<int>& devs, int block) : m_devs(devs), m_block(block) {}

    ~MpEigensolver() override
    {
        releaseMatrices();
        // Grid destruction is collective (cuBLASMp documents it; a sequential loop deadlocks at
        // exit), so every rank tears down on its own thread.
        onAllRanks(static_cast<int>(m_ranks.size()), [&](int r) {
            Rank& R = m_ranks[r];
            cudaSetDevice(m_devs[r]);
            if (R.bgrid) cublasMpGridDestroy(R.bgrid);
            if (R.bhandle) cublasMpDestroy(R.bhandle);
            if (R.grid) cusolverMpDestroyGrid(R.grid);
            if (R.handle) cusolverMpDestroy(R.handle);
            if (R.stream) cudaStreamDestroy(R.stream);
            if (R.comm) ncclCommDestroy(R.comm);
            return true;
        });
    }

    bool init(std::string& why)
    {
        const int ndev = static_cast<int>(m_devs.size());
        enablePeers(m_devs);
        std::vector<ncclComm_t> comms(ndev);
        if (ncclCommInitAll(comms.data(), ndev, m_devs.data()) != ncclSuccess) {
            why = "ncclCommInitAll failed";
            return false;
        }
        m_ranks.resize(ndev);
        for (int r = 0; r < ndev; ++r) {
            Rank& R = m_ranks[r];
            R.comm = comms[r];
            cudaSetDevice(m_devs[r]);
            if (cudaStreamCreate(&R.stream) != cudaSuccess
                || cusolverMpCreate(&R.handle, m_devs[r], R.stream) != CUSOLVER_STATUS_SUCCESS
                || cusolverMpCreateDeviceGrid(R.handle, &R.grid, R.comm, 1, ndev,
                                              CUSOLVERMP_GRID_MAPPING_COL_MAJOR) != CUSOLVER_STATUS_SUCCESS) {
                why = "cuSOLVERMp handle/grid creation failed on device " + std::to_string(m_devs[r]);
                return false;
            }
        }
        // cuBLASMp (distributed trsm for the back-transform) on the same communicators. Optional:
        // without it only solve() is offered.
        m_have_blas = onAllRanks(ndev, [&](int r) {
            Rank& R = m_ranks[r];
            cudaSetDevice(m_devs[r]);
            return cublasMpCreate(&R.bhandle, R.stream) == CUBLASMP_STATUS_SUCCESS
                && cublasMpGridCreate(1, ndev, CUBLASMP_GRID_LAYOUT_COL_MAJOR, R.comm, &R.bgrid)
                       == CUBLASMP_STATUS_SUCCESS;
        });
        return true;
    }

    const char* name() const override { return "cuSOLVERMp"; }
    bool supportsGeneralized() const override { return m_have_blas; }
    void releaseBuffers() override { releaseMatrices(); }
    bool hasMetric(int n, bool fp32, long l_generation) const override
    {
        return m_metric_valid && n == m_n && fp32 == m_fp32 && l_generation == m_metric_gen;
    }
    int deviceCount() const override { return static_cast<int>(m_devs.size()); }

    bool solve(int n, void* A, void* eig, bool fp32, int src_device) override
    {
        const int ndev = static_cast<int>(m_devs.size());
        m_input_intact = true;
        if (!ensureMatrices(n, fp32)) return false;
        const size_t esize = fp32 ? sizeof(float) : sizeof(double);
        char jobz[] = "V";

        // Scatter: each rank pulls its column blocks from the source device (peer copies).
        auto t0 = std::chrono::steady_clock::now();
        bool ok = onAllRanks(ndev, [&](int r) {
            cudaSetDevice(m_devs[r]);
            for (const auto& b : m_layout) {
                if (b.rank != r) continue;
                const size_t bytes = esize * static_cast<size_t>(n) * b.ncol;
                char* dst = static_cast<char*>(m_ranks[r].dA) + esize * static_cast<size_t>(n) * b.local_col0;
                const char* src = static_cast<const char*>(A) + esize * static_cast<size_t>(n) * b.col0;
                if (cudaMemcpyPeer(dst, m_devs[r], src, src_device, bytes) != cudaSuccess) return false;
            }
            return true;
        });
        if (!ok) return false;
        m_scatter_ms = msSince(t0);

        // Collective eigensolve.
        t0 = std::chrono::steady_clock::now();
        ok = onAllRanks(ndev, [&](int r) {
            Rank& R = m_ranks[r];
            cudaSetDevice(m_devs[r]);
            const auto st = cusolverMpSyevd(R.handle, jobz, CUBLAS_FILL_MODE_LOWER, n, R.dA, 1, 1, R.descA,
                                            R.dD, R.dQ, 1, 1, R.descQ, fp32 ? CUDA_R_32F : CUDA_R_64F,
                                            R.dWork, R.work_dev, R.hWork.data(), R.hWork.size(),
                                            static_cast<int*>(R.dInfo));
            if (cudaStreamSynchronize(R.stream) != cudaSuccess) return false;
            int info = -1;
            if (cudaMemcpy(&info, R.dInfo, sizeof(int), cudaMemcpyDeviceToHost) != cudaSuccess) return false;
            return st == CUSOLVER_STATUS_SUCCESS && info == 0;
        });
        if (!ok) return false;
        m_solve_ms = msSince(t0);

        // Gather eigenvectors back into A, eigenvalues from rank 0.
        t0 = std::chrono::steady_clock::now();
        m_input_intact = false;
        ok = onAllRanks(ndev, [&](int r) {
            cudaSetDevice(m_devs[r]);
            for (const auto& b : m_layout) {
                if (b.rank != r) continue;
                const size_t bytes = esize * static_cast<size_t>(n) * b.ncol;
                const char* src = static_cast<const char*>(m_ranks[r].dQ) + esize * static_cast<size_t>(n) * b.local_col0;
                char* dst = static_cast<char*>(A) + esize * static_cast<size_t>(n) * b.col0;
                if (cudaMemcpyPeer(dst, src_device, src, m_devs[r], bytes) != cudaSuccess) return false;
            }
            return true;
        });
        if (!ok) return false;
        const bool eig_ok = cudaMemcpyPeer(eig, src_device, m_ranks[0].dD, m_devs[0], esize * n) == cudaSuccess;
        m_gather_ms = msSince(t0);
        return eig_ok;
    }

    bool solveGeneralized(int n, void* A, const void* L, long l_generation, void* eig, bool fp32,
                          int src_device) override
    {
        const int ndev = static_cast<int>(m_devs.size());
        m_input_intact = true;
        m_reduce_ms = m_back_ms = 0.0;
        if (!m_have_blas || !ensureMatrices(n, fp32) || !ensureGeneralizedBuffers(n, fp32)) return false;
        const size_t esize = fp32 ? sizeof(float) : sizeof(double);
        const cudaDataType dt = fp32 ? CUDA_R_32F : CUDA_R_64F;
        const cublasComputeType_t ct = fp32 ? CUBLAS_COMPUTE_32F : CUBLAS_COMPUTE_64F;
        char jobz[] = "V";
        const bool need_metric = !hasMetric(n, fp32, l_generation);
        if (need_metric && !L) return false;

        // Scatter F (and L when it changed) column blocks.
        auto t0 = std::chrono::steady_clock::now();
        bool ok = scatterColumns(n, esize, A, src_device, [](Rank& R) { return R.dA; })
               && (!need_metric || scatterColumns(n, esize, L, src_device, [](Rank& R) { return R.dL; }));
        if (!ok) return false;
        if (need_metric) {
            m_metric_valid = true;
            m_metric_gen = l_generation;
        }
        m_scatter_ms = msSince(t0);

        t0 = std::chrono::steady_clock::now();
        std::vector<double> t_red(ndev, 0.0), t_eig(ndev, 0.0), t_back(ndev, 0.0);
        ok = onAllRanks(ndev, [&](int r) {
            Rank& R = m_ranks[r];
            cudaSetDevice(m_devs[r]);
            auto tr = std::chrono::steady_clock::now();
            // A = L^-1 F L^-T (lower triangle).
            if (cusolverMpSygst(R.handle, CUSOLVER_EIG_TYPE_1, CUBLAS_FILL_MODE_LOWER, n, R.dA, 1, 1, R.descA,
                                R.dL, 1, 1, R.descA, dt, R.gstWork, R.gst_dev, R.gstHost.data(),
                                R.gstHost.size(), R.mInfo) != CUSOLVER_STATUS_SUCCESS)
                return false;
            if (cudaStreamSynchronize(R.stream) != cudaSuccess || *R.mInfo != 0) return false;
            t_red[r] = msSince(tr);
            tr = std::chrono::steady_clock::now();
            const auto st = cusolverMpSyevd(R.handle, jobz, CUBLAS_FILL_MODE_LOWER, n, R.dA, 1, 1, R.descA,
                                            R.dD, R.dQ, 1, 1, R.descQ, dt, R.dWork, R.work_dev,
                                            R.hWork.data(), R.hWork.size(), static_cast<int*>(R.dInfo));
            if (st != CUSOLVER_STATUS_SUCCESS || cudaStreamSynchronize(R.stream) != cudaSuccess) return false;
            int info = -1;
            if (cudaMemcpy(&info, R.dInfo, sizeof(int), cudaMemcpyDeviceToHost) != cudaSuccess || info != 0)
                return false;
            t_eig[r] = msSince(tr);
            tr = std::chrono::steady_clock::now();
            // C = L^-T Q.
            const double one64 = 1.0;
            const float one32 = 1.0f;
            const void* alpha = fp32 ? static_cast<const void*>(&one32) : static_cast<const void*>(&one64);
            if (cublasMpTrsm(R.bhandle, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_T,
                             CUBLAS_DIAG_NON_UNIT, n, n, alpha, R.dL, 1, 1, R.bdesc, R.dQ, 1, 1, R.bdesc, ct,
                             R.trsmWork, R.trsm_dev, R.trsmHost.data(), R.trsmHost.size())
                != CUBLASMP_STATUS_SUCCESS)
                return false;
            if (cudaStreamSynchronize(R.stream) != cudaSuccess) return false;
            t_back[r] = msSince(tr);
            return true;
        });
        if (!ok) return false;
        m_reduce_ms = *std::max_element(t_red.begin(), t_red.end());
        m_back_ms = *std::max_element(t_back.begin(), t_back.end());
        m_solve_ms = msSince(t0);

        t0 = std::chrono::steady_clock::now();
        m_input_intact = false;
        ok = gatherColumns(n, esize, A, src_device, [](Rank& R) { return R.dQ; });
        if (!ok) return false;
        const bool eig_ok = cudaMemcpyPeer(eig, src_device, m_ranks[0].dD, m_devs[0], esize * n) == cudaSuccess;
        m_gather_ms = msSince(t0);
        return eig_ok;
    }

private:
    struct Rank {
        ncclComm_t comm = nullptr;
        cudaStream_t stream = nullptr;
        cusolverMpHandle_t handle = nullptr;
        cusolverMpGrid_t grid = nullptr;
        cusolverMpMatrixDescriptor_t descA = nullptr, descQ = nullptr;
        void* dA = nullptr;
        void* dQ = nullptr;
        void* dD = nullptr;
        void* dWork = nullptr;
        void* dInfo = nullptr;
        size_t work_dev = 0;
        std::vector<char> hWork;
        // Generalized path (cuBLASMp trsm + cuSOLVERMp sygst).
        cublasMpHandle_t bhandle = nullptr;
        cublasMpGrid_t bgrid = nullptr;
        cublasMpMatrixDescriptor_t bdesc = nullptr;
        void* dL = nullptr;
        void* gstWork = nullptr;
        size_t gst_dev = 0;
        std::vector<char> gstHost;
        void* trsmWork = nullptr;
        size_t trsm_dev = 0;
        std::vector<char> trsmHost;
        int* mInfo = nullptr;   // managed: readable on host and device whichever sygst expects
    };

    template <class BufFn>
    bool scatterColumns(int n, size_t esize, const void* src_mat, int src_device, BufFn buf)
    {
        return onAllRanks(static_cast<int>(m_devs.size()), [&](int r) {
            cudaSetDevice(m_devs[r]);
            for (const auto& b : m_layout) {
                if (b.rank != r) continue;
                const size_t bytes = esize * static_cast<size_t>(n) * b.ncol;
                char* dst = static_cast<char*>(buf(m_ranks[r])) + esize * static_cast<size_t>(n) * b.local_col0;
                const char* src = static_cast<const char*>(src_mat) + esize * static_cast<size_t>(n) * b.col0;
                if (cudaMemcpyPeer(dst, m_devs[r], src, src_device, bytes) != cudaSuccess) return false;
            }
            return true;
        });
    }

    template <class BufFn>
    bool gatherColumns(int n, size_t esize, void* dst_mat, int dst_device, BufFn buf)
    {
        return onAllRanks(static_cast<int>(m_devs.size()), [&](int r) {
            cudaSetDevice(m_devs[r]);
            for (const auto& b : m_layout) {
                if (b.rank != r) continue;
                const size_t bytes = esize * static_cast<size_t>(n) * b.ncol;
                const char* src = static_cast<const char*>(buf(m_ranks[r])) + esize * static_cast<size_t>(n) * b.local_col0;
                char* dst = static_cast<char*>(dst_mat) + esize * static_cast<size_t>(n) * b.col0;
                if (cudaMemcpyPeer(dst, dst_device, src, m_devs[r], bytes) != cudaSuccess) return false;
            }
            return true;
        });
    }

    void releaseMatrices()
    {
        for (int r = 0; r < static_cast<int>(m_ranks.size()); ++r) {
            Rank& R = m_ranks[r];
            cudaSetDevice(m_devs[r]);
            for (void** p : { &R.dA, &R.dQ, &R.dD, &R.dWork, &R.dInfo, &R.dL, &R.gstWork, &R.trsmWork })
                if (*p) { cudaFree(*p); *p = nullptr; }
            if (R.mInfo) { cudaFree(R.mInfo); R.mInfo = nullptr; }
            if (R.bdesc) { cublasMpMatrixDescriptorDestroy(R.bdesc); R.bdesc = nullptr; }
            R.gstHost.clear();
            R.trsmHost.clear();
            R.gst_dev = R.trsm_dev = 0;
            if (R.descA) { cusolverMpDestroyMatrixDesc(R.descA); R.descA = nullptr; }
            if (R.descQ) { cusolverMpDestroyMatrixDesc(R.descQ); R.descQ = nullptr; }
            R.hWork.clear();
            R.work_dev = 0;
        }
        m_n = 0;
        m_metric_valid = false;
        m_gen_alloc = false;
    }

    /// Distributed L plus sygst/trsm workspaces, only for solveGeneralized() (the plain solve()
    /// path does not carry them). Queries are collective.
    bool ensureGeneralizedBuffers(int n, bool fp32)
    {
        if (m_gen_alloc) return true;
        const int ndev = static_cast<int>(m_devs.size());
        std::vector<long> local_cols;
        columnLayout(n, m_block, ndev, local_cols);
        const size_t esize = fp32 ? sizeof(float) : sizeof(double);
        const cudaDataType dt = fp32 ? CUDA_R_32F : CUDA_R_64F;
        bool ok = true;
        // Generalized path buffers: distributed L plus sygst and trsm workspaces (queries collective).
        if (ok) {
            const cublasComputeType_t ct = fp32 ? CUBLAS_COMPUTE_32F : CUBLAS_COMPUTE_64F;
            for (int r = 0; r < ndev && ok; ++r) {
                Rank& R = m_ranks[r];
                cudaSetDevice(m_devs[r]);
                const size_t lbytes = esize * static_cast<size_t>(n) * std::max<long>(1, local_cols[r]);
                void* info_mem = nullptr;
                ok = cudaMalloc(&R.dL, lbytes) == cudaSuccess
                  && cudaMallocManaged(&info_mem, sizeof(int), cudaMemAttachGlobal) == cudaSuccess
                  && cublasMpMatrixDescriptorCreate(n, n, m_block, m_block, 0, 0, n, dt, R.bgrid, &R.bdesc)
                         == CUBLASMP_STATUS_SUCCESS;
                R.mInfo = static_cast<int*>(info_mem);
            }
            ok = ok && onAllRanks(ndev, [&](int r) {
                Rank& R = m_ranks[r];
                cudaSetDevice(m_devs[r]);
                size_t gh = 0, th = 0;
                const double one64 = 1.0;
                const float one32 = 1.0f;
                const void* alpha = fp32 ? static_cast<const void*>(&one32) : static_cast<const void*>(&one64);
                if (cusolverMpSygst_bufferSize(R.handle, CUSOLVER_EIG_TYPE_1, CUBLAS_FILL_MODE_LOWER, n, 1, 1,
                                               R.descA, 1, 1, R.descA, dt, &R.gst_dev, &gh)
                    != CUSOLVER_STATUS_SUCCESS)
                    return false;
                if (cublasMpTrsm_bufferSize(R.bhandle, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_T,
                                            CUBLAS_DIAG_NON_UNIT, n, n, alpha, R.dL, 1, 1, R.bdesc, R.dQ, 1, 1,
                                            R.bdesc, ct, &R.trsm_dev, &th)
                    != CUBLASMP_STATUS_SUCCESS)
                    return false;
                R.gstHost.assign(std::max<size_t>(1, gh), 0);
                R.trsmHost.assign(std::max<size_t>(1, th), 0);
                return cudaMalloc(&R.gstWork, std::max<size_t>(1, R.gst_dev)) == cudaSuccess
                    && cudaMalloc(&R.trsmWork, std::max<size_t>(1, R.trsm_dev)) == cudaSuccess;
            });
        }
        m_gen_alloc = ok;
        return ok;
    }

    bool ensureMatrices(int n, bool fp32)
    {
        if (n == m_n && fp32 == m_fp32) return true;
        releaseMatrices();
        const int ndev = static_cast<int>(m_devs.size());
        std::vector<long> local_cols;
        m_layout = columnLayout(n, m_block, ndev, local_cols);
        const size_t esize = fp32 ? sizeof(float) : sizeof(double);
        const cudaDataType dt = fp32 ? CUDA_R_32F : CUDA_R_64F;
        char jobz[] = "V";
        bool ok = true;
        for (int r = 0; r < ndev && ok; ++r) {
            Rank& R = m_ranks[r];
            cudaSetDevice(m_devs[r]);
            const size_t lbytes = esize * static_cast<size_t>(n) * std::max<long>(1, local_cols[r]);
            ok = cusolverMpCreateMatrixDesc(&R.descA, R.grid, dt, n, n, m_block, m_block, 0, 0, n) == CUSOLVER_STATUS_SUCCESS
              && cusolverMpCreateMatrixDesc(&R.descQ, R.grid, dt, n, n, m_block, m_block, 0, 0, n) == CUSOLVER_STATUS_SUCCESS
              && cudaMalloc(&R.dA, lbytes) == cudaSuccess
              && cudaMalloc(&R.dQ, lbytes) == cudaSuccess
              && cudaMalloc(&R.dD, esize * n) == cudaSuccess
              && cudaMalloc(&R.dInfo, sizeof(int)) == cudaSuccess;
        }
        // Workspace query is itself collective in cuSOLVERMp.
        ok = ok && onAllRanks(ndev, [&](int r) {
            Rank& R = m_ranks[r];
            cudaSetDevice(m_devs[r]);
            size_t whost = 0;
            if (cusolverMpSyevd_bufferSize(R.handle, jobz, CUBLAS_FILL_MODE_LOWER, n, R.dA, 1, 1, R.descA,
                                           R.dD, R.dQ, 1, 1, R.descQ, dt, &R.work_dev, &whost)
                != CUSOLVER_STATUS_SUCCESS)
                return false;
            R.hWork.assign(std::max<size_t>(1, whost), 0);
            return cudaMalloc(&R.dWork, std::max<size_t>(1, R.work_dev)) == cudaSuccess;
        });
        if (!ok) {
            releaseMatrices();
            return false;
        }
        m_n = n;
        m_fp32 = fp32;
        return true;
    }

    std::vector<int> m_devs;
    int m_block;
    bool m_have_blas = false;
    bool m_metric_valid = false;
    bool m_gen_alloc = false;
    long m_metric_gen = -1;
    std::vector<Rank> m_ranks;
    std::vector<ColumnBlock> m_layout;
    int m_n = 0;
    bool m_fp32 = false;
};
#endif // CURCUMA_HAVE_CUSOLVERMP

#ifdef CURCUMA_HAVE_CUSOLVERMG
// ---------------------------------------------------------------------------------------------
// cusolverMg: single call over all devices (CUDA toolkit, deprecated in CUDA 13).
// ---------------------------------------------------------------------------------------------
class MgEigensolver final : public DistributedEigensolver {
public:
    MgEigensolver(const std::vector<int>& devs, int block) : m_devs(devs), m_block(block) {}

    ~MgEigensolver() override
    {
        releaseMatrices();
        if (m_grid) cusolverMgDestroyGrid(m_grid);
        if (m_handle) cusolverMgDestroy(m_handle);
    }

    bool init(std::string& why)
    {
        enablePeers(m_devs);
        if (cusolverMgCreate(&m_handle) != CUSOLVER_STATUS_SUCCESS
            || cusolverMgDeviceSelect(m_handle, static_cast<int>(m_devs.size()), m_devs.data()) != CUSOLVER_STATUS_SUCCESS
            || cusolverMgCreateDeviceGrid(&m_grid, 1, static_cast<int>(m_devs.size()), m_devs.data(),
                                          CUDALIBMG_GRID_MAPPING_COL_MAJOR) != CUSOLVER_STATUS_SUCCESS) {
            why = "cusolverMg handle/grid creation failed";
            return false;
        }
        return true;
    }

    const char* name() const override { return "cusolverMg"; }
    void releaseBuffers() override { releaseMatrices(); }
    int deviceCount() const override { return static_cast<int>(m_devs.size()); }

    bool solve(int n, void* A, void* eig, bool fp32, int src_device) override
    {
        const int ndev = static_cast<int>(m_devs.size());
        if (!ensureMatrices(n, fp32)) return false;
        const size_t esize = fp32 ? sizeof(float) : sizeof(double);
        const cudaDataType dt = fp32 ? CUDA_R_32F : CUDA_R_64F;
        auto t0 = std::chrono::steady_clock::now();
        for (const auto& b : m_layout) {
            const size_t bytes = esize * static_cast<size_t>(n) * b.ncol;
            char* dst = static_cast<char*>(m_dA[b.rank]) + esize * static_cast<size_t>(n) * b.local_col0;
            const char* src = static_cast<const char*>(A) + esize * static_cast<size_t>(n) * b.col0;
            if (cudaMemcpyPeer(dst, m_devs[b.rank], src, src_device, bytes) != cudaSuccess) return false;
        }
        int info = 0;
        m_input_intact = true;
        m_scatter_ms = msSince(t0);
        t0 = std::chrono::steady_clock::now();
        const auto st = cusolverMgSyevd(m_handle, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_LOWER, n,
                                        m_dA.data(), 1, 1, m_desc, m_hW.data(), dt, dt, m_dWork.data(),
                                        m_lwork, &info);
        for (int d : m_devs) { cudaSetDevice(d); cudaDeviceSynchronize(); }
        if (st != CUSOLVER_STATUS_SUCCESS || info != 0) return false;
        m_solve_ms = msSince(t0);
        t0 = std::chrono::steady_clock::now();
        m_input_intact = false;
        for (const auto& b : m_layout) {
            const size_t bytes = esize * static_cast<size_t>(n) * b.ncol;
            const char* src = static_cast<const char*>(m_dA[b.rank]) + esize * static_cast<size_t>(n) * b.local_col0;
            char* dst = static_cast<char*>(A) + esize * static_cast<size_t>(n) * b.col0;
            if (cudaMemcpyPeer(dst, src_device, src, m_devs[b.rank], bytes) != cudaSuccess) return false;
        }
        cudaSetDevice(src_device);
        (void)ndev;
        const bool eig_ok = cudaMemcpy(eig, m_hW.data(), esize * n, cudaMemcpyHostToDevice) == cudaSuccess;
        m_gather_ms = msSince(t0);
        return eig_ok;
    }

private:
    void releaseMatrices()
    {
        for (size_t r = 0; r < m_dA.size(); ++r) {
            cudaSetDevice(m_devs[r]);
            if (m_dA[r]) cudaFree(m_dA[r]);
            if (m_dWork[r]) cudaFree(m_dWork[r]);
        }
        m_dA.clear();
        m_dWork.clear();
        if (m_desc) { cusolverMgDestroyMatrixDesc(m_desc); m_desc = nullptr; }
        m_n = 0;
    }

    bool ensureMatrices(int n, bool fp32)
    {
        if (n == m_n && fp32 == m_fp32) return true;
        releaseMatrices();
        const int ndev = static_cast<int>(m_devs.size());
        std::vector<long> local_cols;
        m_layout = columnLayout(n, m_block, ndev, local_cols);
        const size_t esize = fp32 ? sizeof(float) : sizeof(double);
        const cudaDataType dt = fp32 ? CUDA_R_32F : CUDA_R_64F;
        if (cusolverMgCreateMatrixDesc(&m_desc, n, n, n, m_block, dt, m_grid) != CUSOLVER_STATUS_SUCCESS)
            return false;
        m_dA.assign(ndev, nullptr);
        m_dWork.assign(ndev, nullptr);
        for (int r = 0; r < ndev; ++r) {
            cudaSetDevice(m_devs[r]);
            if (cudaMalloc(&m_dA[r], esize * static_cast<size_t>(n) * std::max<long>(1, local_cols[r])) != cudaSuccess) {
                releaseMatrices();
                return false;
            }
        }
        m_hW.assign(esize * n, 0);
        if (cusolverMgSyevd_bufferSize(m_handle, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_LOWER, n, m_dA.data(),
                                       1, 1, m_desc, m_hW.data(), dt, dt, &m_lwork) != CUSOLVER_STATUS_SUCCESS) {
            releaseMatrices();
            return false;
        }
        for (int r = 0; r < ndev; ++r) {
            cudaSetDevice(m_devs[r]);
            if (cudaMalloc(&m_dWork[r], esize * static_cast<size_t>(std::max<int64_t>(1, m_lwork))) != cudaSuccess) {
                releaseMatrices();
                return false;
            }
        }
        m_n = n;
        m_fp32 = fp32;
        return true;
    }

    std::vector<int> m_devs;
    int m_block;
    cusolverMgHandle_t m_handle = nullptr;
    cudaLibMgGrid_t m_grid = nullptr;
    cudaLibMgMatrixDesc_t m_desc = nullptr;
    std::vector<ColumnBlock> m_layout;
    std::vector<void*> m_dA, m_dWork;
    std::vector<char> m_hW;
    int64_t m_lwork = 0;
    int m_n = 0;
    bool m_fp32 = false;
};
#endif // CURCUMA_HAVE_CUSOLVERMG

} // namespace

// ---------------------------------------------------------------------------------------------
// C entry points of libcurcuma_cuda_mgpu.so, resolved by xtb_distributed_eigensolver_loader.cpp.
// ---------------------------------------------------------------------------------------------
static std::string mgpuBackends()
{
    std::string s;
#ifdef CURCUMA_HAVE_CUSOLVERMP
    s += "mp";
#endif
#ifdef CURCUMA_HAVE_CUSOLVERMG
    s += s.empty() ? "mg" : ",mg";
#endif
    return s.empty() ? "none" : s;
}

static DistributedEigensolver* mgpuCreate(const std::string& backend, const std::vector<int>& devices,
                                          int block, std::string& why)
{
    const int blk = block > 0 ? block : 128;
    int saved = 0;
    cudaGetDevice(&saved);
    std::unique_ptr<DistributedEigensolver> out;
    const bool want_mp = (backend == "auto" || backend == "mp");
    const bool want_mg = (backend == "auto" || backend == "mg");
#ifdef CURCUMA_HAVE_CUSOLVERMP
    if (!out && want_mp) {
        auto mp = std::make_unique<MpEigensolver>(devices, blk);
        if (mp->init(why)) out = std::move(mp);
    }
#endif
#ifdef CURCUMA_HAVE_CUSOLVERMG
    if (!out && want_mg) {
        why.clear();
        auto mg = std::make_unique<MgEigensolver>(devices, blk);
        if (mg->init(why)) out = std::move(mg);
    }
#endif
    (void)want_mp; (void)want_mg; (void)blk;
    if (!out && why.empty())
        why = "backend '" + backend + "' not built (available: " + mgpuBackends() + ")";
    cudaSetDevice(saved);
    return out.release();
}

} // namespace gpu
} // namespace xtb
} // namespace curcuma

extern "C" {
// Returns "mp,mg" / "mp" / "mg" / "none"; the string lives for the process lifetime.
const char* curcuma_cuda_mgpu_backends()
{
    static const std::string s = curcuma::xtb::gpu::mgpuBackends();
    return s.c_str();
}

// Caller owns the returned object (deleted through the virtual destructor while this library
// stays loaded). On nullptr, `why` holds the reason.
curcuma::xtb::gpu::DistributedEigensolver* curcuma_cuda_mgpu_create(const std::string* backend,
                                                                    const std::vector<int>* devices,
                                                                    int block, std::string* why)
{
    return curcuma::xtb::gpu::mgpuCreate(*backend, *devices, block, *why);
}
}

#endif // USE_CUDA
