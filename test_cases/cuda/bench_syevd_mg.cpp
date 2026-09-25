/*
 * < Multi-GPU dense symmetric eigensolver microbenchmark >
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

// Claude Generated (Sep 2026): decision gate G0d of the multi-GPU plan (docs/MULTI_GPU.md).
//
// The dense SCF eigensolve (F C = S C e, after Cholesky reduction a standard symmetric
// problem) is O(n^3) and dominates native GFN1/GFN2 for large systems. This standalone tool
// measures, for a random symmetric n x n matrix:
//   - cusolverDnXsyevd on ONE device (FP64 / FP32), data already on the device
//   - cusolverMgSyevd over 1..k devices (FP64 / FP32), 1D column-block distribution,
//     plus the host<->device transfer cost of scattering the matrix and gathering the vectors
// and reports the eigenvalue agreement between the two solvers.
//
// Build (no CMake needed):
//   g++ -O2 -std=c++17 bench_syevd_mg.cpp -I/opt/cuda/include -L/opt/cuda/lib64 \
//       -lcusolverMg -lcusolver -lcublas -lcudart -o bench_syevd_mg
// Run:
//   ./bench_syevd_mg <n> <fp64|fp32> <devices e.g. 0,1,2,3> [block_size=1024]
// Use CUDA_VISIBLE_DEVICES to restrict what the process sees (SLURM does this for you).

#define DISABLE_CUSOLVERMG_DEPRECATED

#include <cuda.h>
#include <cudalibxt.h>
#include <cusolverDn.h>
#include <cusolverMg.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <sstream>
#include <string>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;
double secondsSince(Clock::time_point t0)
{
    return std::chrono::duration<double>(Clock::now() - t0).count();
}

void check(bool ok, const char* what)
{
    if (!ok) {
        std::fprintf(stderr, "FAILED: %s\n", what);
        std::exit(1);
    }
}

// Random symmetric matrix with entries ~ U(-1,1) plus a graded diagonal, so the spectrum is
// spread (closer to a Fock matrix than a pure Wigner matrix). Column-major, full storage.
template <typename T>
std::vector<T> makeSymmetric(int n)
{
    std::mt19937_64 rng(42);
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    std::vector<T> a(static_cast<size_t>(n) * n);
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i <= j; ++i) {
            const double v = u(rng) / std::sqrt(static_cast<double>(n));
            a[i + static_cast<size_t>(j) * n] = static_cast<T>(v);
            a[j + static_cast<size_t>(i) * n] = static_cast<T>(v);
        }
        a[j + static_cast<size_t>(j) * n] += static_cast<T>(-20.0 + 40.0 * j / n);
    }
    return a;
}

void syncAll(const std::vector<int>& devs)
{
    for (int d : devs) {
        cudaSetDevice(d);
        cudaDeviceSynchronize();
    }
}

template <typename T>
std::vector<T> runDn(int dev, int n, const std::vector<T>& host, cudaDataType dtype)
{
    cudaSetDevice(dev);
    cusolverDnHandle_t h = nullptr;
    check(cusolverDnCreate(&h) == CUSOLVER_STATUS_SUCCESS, "cusolverDnCreate");
    cusolverDnParams_t params = nullptr;
    check(cusolverDnCreateParams(&params) == CUSOLVER_STATUS_SUCCESS, "params");

    const size_t bytes = sizeof(T) * static_cast<size_t>(n) * n;
    void* dA = nullptr;
    void* dW = nullptr;
    check(cudaMalloc(&dA, bytes) == cudaSuccess, "Dn malloc A");
    check(cudaMalloc(&dW, sizeof(T) * n) == cudaSuccess, "Dn malloc W");

    auto t0 = Clock::now();
    cudaMemcpy(dA, host.data(), bytes, cudaMemcpyHostToDevice);
    cudaDeviceSynchronize();
    const double t_h2d = secondsSince(t0);

    size_t lwork_dev = 0, lwork_host = 0;
    check(cusolverDnXsyevd_bufferSize(h, params, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_LOWER,
                                      n, dtype, dA, n, dtype, dW, dtype, &lwork_dev, &lwork_host)
              == CUSOLVER_STATUS_SUCCESS,
          "Dn bufferSize");
    void* dWork = nullptr;
    check(cudaMalloc(&dWork, lwork_dev) == cudaSuccess, "Dn malloc work");
    std::vector<char> hWork(lwork_host);
    void* dInfo = nullptr;
    cudaMalloc(&dInfo, sizeof(int));

    size_t free_b = 0, total_b = 0;
    cudaMemGetInfo(&free_b, &total_b);

    t0 = Clock::now();
    check(cusolverDnXsyevd(h, params, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_LOWER, n, dtype, dA, n,
                           dtype, dW, dtype, dWork, lwork_dev, hWork.data(), lwork_host,
                           static_cast<int*>(dInfo))
              == CUSOLVER_STATUS_SUCCESS,
          "Dn syevd");
    cudaDeviceSynchronize();
    const double t_solve = secondsSince(t0);

    std::vector<T> w(n);
    t0 = Clock::now();
    std::vector<T> vecs(static_cast<size_t>(n) * n);
    cudaMemcpy(vecs.data(), dA, bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(w.data(), dW, sizeof(T) * n, cudaMemcpyDeviceToHost);
    const double t_d2h = secondsSince(t0);

    std::printf("DN   dev=%d n=%d %s solve=%.2f s  h2d=%.2f s d2h=%.2f s  used=%.2f GB\n", dev, n,
                sizeof(T) == 8 ? "fp64" : "fp32", t_solve, t_h2d, t_d2h,
                (total_b - free_b) / 1073741824.0);

    cudaFree(dA); cudaFree(dW); cudaFree(dWork); cudaFree(dInfo);
    cusolverDnDestroyParams(params);
    cusolverDnDestroy(h);
    return w;
}

template <typename T>
std::vector<T> runMg(const std::vector<int>& devs, int n, int block, const std::vector<T>& host,
                     cudaDataType dtype)
{
    const int ndev = static_cast<int>(devs.size());
    cusolverMgHandle_t h = nullptr;
    check(cusolverMgCreate(&h) == CUSOLVER_STATUS_SUCCESS, "cusolverMgCreate");
    check(cusolverMgDeviceSelect(h, ndev, const_cast<int*>(devs.data())) == CUSOLVER_STATUS_SUCCESS,
          "DeviceSelect");
    cudaLibMgGrid_t grid = nullptr;
    check(cusolverMgCreateDeviceGrid(&grid, 1, ndev, devs.data(), CUDALIBMG_GRID_MAPPING_COL_MAJOR)
              == CUSOLVER_STATUS_SUCCESS,
          "DeviceGrid");
    cudaLibMgMatrixDesc_t desc = nullptr;
    check(cusolverMgCreateMatrixDesc(&desc, n, n, n, block, dtype, grid) == CUSOLVER_STATUS_SUCCESS,
          "MatrixDesc");

    // 1D column distribution (verified empirically, Sep 2026): the n columns are cut into
    // ceil(n/block) blocks and the blocks are dealt out CYCLICALLY, block b -> device b % ndev.
    // Each device stores its own blocks contiguously, in increasing block order, lda = n.
    // (A contiguous "first half on device 0" layout returns wrong eigenvalues without error.)
    const int nblocks = (n + block - 1) / block;
    struct Piece { int dev; long col0; long ncol; long local0; };
    std::vector<Piece> pieces;
    std::vector<long> local_cols(ndev, 0);
    for (int b = 0; b < nblocks; ++b) {
        const int g = b % ndev;
        const long c0 = static_cast<long>(b) * block;
        const long nc = std::min<long>(block, n - c0);
        pieces.push_back({ g, c0, nc, local_cols[g] });
        local_cols[g] += nc;
    }
    std::vector<void*> dA(ndev, nullptr);
    auto copyPieces = [&](bool toDevice, T* hostbuf) {
        for (const auto& pc : pieces) {
            cudaSetDevice(devs[pc.dev]);
            const size_t bytes = sizeof(T) * static_cast<size_t>(n) * pc.ncol;
            T* h_ptr = hostbuf + static_cast<size_t>(pc.col0) * n;
            T* d_ptr = static_cast<T*>(dA[pc.dev]) + static_cast<size_t>(pc.local0) * n;
            if (toDevice)
                cudaMemcpy(d_ptr, h_ptr, bytes, cudaMemcpyHostToDevice);
            else
                cudaMemcpy(h_ptr, d_ptr, bytes, cudaMemcpyDeviceToHost);
        }
    };

    for (int g = 0; g < ndev; ++g) {
        cudaSetDevice(devs[g]);
        const size_t b = sizeof(T) * static_cast<size_t>(n) * std::max<long>(1, local_cols[g]);
        check(cudaMalloc(&dA[g], b) == cudaSuccess, "Mg malloc A");
    }

    auto t0 = Clock::now();
    std::vector<T> work_host(host);
    copyPieces(true, work_host.data());
    syncAll(devs);
    const double t_h2d = secondsSince(t0);

    std::vector<T> w(n);
    int64_t lwork = 0;
    check(cusolverMgSyevd_bufferSize(h, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_LOWER, n, dA.data(),
                                     1, 1, desc, w.data(), dtype, dtype, &lwork)
              == CUSOLVER_STATUS_SUCCESS,
          "Mg bufferSize");
    std::vector<void*> dWork(ndev, nullptr);
    for (int g = 0; g < ndev; ++g) {
        cudaSetDevice(devs[g]);
        check(cudaMalloc(&dWork[g], sizeof(T) * static_cast<size_t>(lwork)) == cudaSuccess,
              "Mg malloc work");
    }
    std::ostringstream used;
    for (int g = 0; g < ndev; ++g) {
        cudaSetDevice(devs[g]);
        size_t f = 0, t = 0;
        cudaMemGetInfo(&f, &t);
        used << (g ? "," : "") << static_cast<int>((t - f) / 1048576);
    }

    int info = 0;
    t0 = Clock::now();
    check(cusolverMgSyevd(h, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_LOWER, n, dA.data(), 1, 1, desc,
                          w.data(), dtype, dtype, dWork.data(), lwork, &info)
              == CUSOLVER_STATUS_SUCCESS,
          "Mg syevd");
    syncAll(devs);
    const double t_solve = secondsSince(t0);
    check(info == 0, "Mg syevd info != 0");

    std::vector<T> vecs(static_cast<size_t>(n) * n);
    t0 = Clock::now();
    copyPieces(false, vecs.data());
    const double t_d2h = secondsSince(t0);

    // Claude Generated (Sep 2026): is the result actually usable? A correct syevd returns
    // eigenvalues in ASCENDING order with column k belonging to w[k]. A permutation-invariant
    // check (trace, or A ~ Q L Q^T) does not catch a wrong order, but the SCF does: it fills the
    // leading columns. Report both the order and a per-column residual.
    {
        int desc_asc = 0;
        for (int i = 1; i < n; ++i)
            if (w[i] < w[i - 1]) ++desc_asc;
        // Residual ||A x_k - w_k x_k|| / |w_k| for a few k, recomputing A x_k on the host would be
        // O(n^2) per column: do 3 columns only.
        double worst = 0.0;
        for (int kk : { 0, n / 2, n - 1 }) {
            std::vector<double> y(n, 0.0);
            for (int j = 0; j < n; ++j) {
                const double xj = static_cast<double>(vecs[static_cast<size_t>(j) + static_cast<size_t>(kk) * n]);
                if (xj == 0.0) continue;
                for (int i = 0; i < n; ++i)
                    y[i] += static_cast<double>(host[static_cast<size_t>(i) + static_cast<size_t>(j) * n]) * xj;
            }
            double num = 0.0, den = 0.0;
            for (int i = 0; i < n; ++i) {
                const double xi = static_cast<double>(vecs[static_cast<size_t>(i) + static_cast<size_t>(kk) * n]);
                const double r = y[i] - static_cast<double>(w[kk]) * xi;
                num += r * r; den += xi * xi;
            }
            const double rel = std::sqrt(num) / (std::sqrt(den) * std::max(1e-30, std::fabs(static_cast<double>(w[kk]))));
            if (rel > worst) worst = rel;
        }
        std::printf("     check: %d of %d eigenvalues out of ascending order, worst column residual %.3e\n",
                    desc_asc, n - 1, worst);
    }

    std::printf("MG   ndev=%d n=%d %s block=%d lwork/dev=%.2f GB solve=%.2f s  h2d=%.2f s d2h=%.2f s"
                "  used_MiB=[%s]\n",
                ndev, n, sizeof(T) == 8 ? "fp64" : "fp32", block,
                sizeof(T) * static_cast<double>(lwork) / 1073741824.0, t_solve, t_h2d, t_d2h,
                used.str().c_str());

    for (int g = 0; g < ndev; ++g) {
        cudaSetDevice(devs[g]);
        cudaFree(dA[g]);
        cudaFree(dWork[g]);
    }
    cusolverMgDestroyMatrixDesc(desc);
    cusolverMgDestroyGrid(grid);
    cusolverMgDestroy(h);
    return w;
}

template <typename T>
void benchmark(int n, const std::vector<int>& devs, int block, cudaDataType dtype)
{
    const auto host = makeSymmetric<T>(n);
    const auto w_dn = runDn<T>(devs[0], n, host, dtype);
    const auto w_mg = runMg<T>(devs, n, block, host, dtype);
    double maxdiff = 0.0;
    for (int i = 0; i < n; ++i)
        maxdiff = std::max(maxdiff, std::abs(static_cast<double>(w_dn[i]) - w_mg[i]));
    double tr = 0.0, s_dn = 0.0, s_mg = 0.0;
    for (int i = 0; i < n; ++i) {
        tr += host[i + static_cast<size_t>(i) * n];
        s_dn += w_dn[i];
        s_mg += w_mg[i];
    }
    std::printf("max |eig_Dn - eig_Mg| = %.3e   trace=%.6f sum_Dn=%.6f sum_Mg=%.6f\n", maxdiff, tr, s_dn, s_mg);
}

} // namespace

int main(int argc, char** argv)
{
    if (argc < 4) {
        std::fprintf(stderr, "usage: %s <n> <fp64|fp32> <dev,dev,...> [block]\n", argv[0]);
        return 1;
    }
    const int n = std::atoi(argv[1]);
    const std::string prec = argv[2];
    std::vector<int> devs;
    std::stringstream ss(argv[3]);
    for (std::string tok; std::getline(ss, tok, ',');)
        devs.push_back(std::atoi(tok.c_str()));
    const int block = argc > 4 ? std::atoi(argv[4]) : 1024;

    // Peer access between all selected devices (cusolverMg uses P2P when available).
    for (int a : devs) {
        cudaSetDevice(a);
        for (int b : devs) {
            int can = 0;
            if (a != b && cudaDeviceCanAccessPeer(&can, a, b) == cudaSuccess && can)
                cudaDeviceEnablePeerAccess(b, 0);
        }
    }

    if (prec == "fp32")
        benchmark<float>(n, devs, block, CUDA_R_32F);
    else
        benchmark<double>(n, devs, block, CUDA_R_64F);
    return 0;
}
