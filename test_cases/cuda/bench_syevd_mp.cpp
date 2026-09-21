/*
 * < Multi-GPU dense symmetric eigensolver microbenchmark: cuSOLVERMp (NCCL) >
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

// Claude Generated (Sep 2026): companion of bench_syevd_mg.cpp for cuSOLVERMp, the NCCL-based
// successor of the deprecated cusolverMg. cuSOLVERMp is ScaLAPACK-style: every GPU is a "rank"
// with its own handle, grid and local matrix block, and every call is collective. curcuma is a
// single process, so this benchmark drives one host thread per GPU with communicators from
// ncclCommInitAll - the same model a curcuma integration would use.
//
// Build (paths for the pip wheels nvidia-cusolvermp-cu13 / nvidia-cublasmp-cu13):
//   g++ -O2 -std=c++17 bench_syevd_mp.cpp -I/opt/cuda/include -I$MP/include \
//       -L/opt/cuda/lib64 -L$MP/lib -L$BMP/lib -lcusolverMp -lcublasmp -lnccl -lcudart -lpthread \
//       -o bench_syevd_mp
// Run:
//   LD_LIBRARY_PATH=$MP/lib:$BMP/lib ./bench_syevd_mp <n> <fp64|fp32> <devices e.g. 0,1,2,3> [block=1024]

#include <cuda.h>
#include <cusolverMp.h>
#include <nccl.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;
double secondsSince(Clock::time_point t0) { return std::chrono::duration<double>(Clock::now() - t0).count(); }

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

template <typename T>
int run(int n, const std::vector<int>& devs, int nb, cudaDataType dtype)
{
    const int ndev = static_cast<int>(devs.size());
    const auto host = makeSymmetric<T>(n);
    std::vector<ncclComm_t> comms(ndev);
    if (ncclCommInitAll(comms.data(), ndev, devs.data()) != ncclSuccess) {
        std::fprintf(stderr, "ncclCommInitAll failed\n");
        return 1;
    }
    std::vector<double> t_scatter(ndev), t_solve(ndev), used_gb(ndev);
    std::vector<T> eig(n);
    std::vector<int> ok(ndev, 0);
    char jobz[] = "V";

    auto worker = [&](int r) {
        cudaSetDevice(devs[r]);
        cudaStream_t stream = nullptr;
        cudaStreamCreate(&stream);
        cusolverMpHandle_t h = nullptr;
        if (cusolverMpCreate(&h, devs[r], stream) != CUSOLVER_STATUS_SUCCESS) return;
        cusolverMpGrid_t grid = nullptr;
        if (cusolverMpCreateDeviceGrid(h, &grid, comms[r], 1, ndev, CUSOLVERMP_GRID_MAPPING_COL_MAJOR)
            != CUSOLVER_STATUS_SUCCESS) return;
        const int64_t lrows = cusolverMpNUMROC(n, nb, 0, 0, 1);
        const int64_t lcols = cusolverMpNUMROC(n, nb, static_cast<uint32_t>(r), 0, static_cast<uint32_t>(ndev));
        cusolverMpMatrixDescriptor_t descA = nullptr, descQ = nullptr;
        cusolverMpCreateMatrixDesc(&descA, grid, dtype, n, n, nb, nb, 0, 0, lrows);
        cusolverMpCreateMatrixDesc(&descQ, grid, dtype, n, n, nb, nb, 0, 0, lrows);
        void* dA = nullptr; void* dQ = nullptr; void* dD = nullptr; void* dInfo = nullptr;
        const size_t lbytes = sizeof(T) * static_cast<size_t>(lrows) * std::max<int64_t>(1, lcols);
        cudaMalloc(&dA, lbytes);
        cudaMalloc(&dQ, lbytes);
        cudaMalloc(&dD, sizeof(T) * n);
        cudaMalloc(&dInfo, sizeof(int));

        auto t0 = Clock::now();
        cusolverMpMatrixScatterH2D(h, n, n, dA, 1, 1, descA, 0, r == 0 ? host.data() : nullptr, n);
        cudaStreamSynchronize(stream);
        t_scatter[r] = secondsSince(t0);

        size_t wdev = 0, whost = 0;
        cusolverMpSyevd_bufferSize(h, jobz, CUBLAS_FILL_MODE_LOWER, n, dA, 1, 1, descA, dD, dQ, 1, 1,
                                   descQ, dtype, &wdev, &whost);
        void* dW = nullptr;
        cudaMalloc(&dW, std::max<size_t>(1, wdev));
        std::vector<char> hW(std::max<size_t>(1, whost));
        size_t f = 0, tot = 0;
        cudaMemGetInfo(&f, &tot);
        used_gb[r] = (tot - f) / 1073741824.0;

        t0 = Clock::now();
        const auto st = cusolverMpSyevd(h, jobz, CUBLAS_FILL_MODE_LOWER, n, dA, 1, 1, descA, dD, dQ, 1, 1,
                                        descQ, dtype, dW, wdev, hW.data(), whost, static_cast<int*>(dInfo));
        cudaStreamSynchronize(stream);
        t_solve[r] = secondsSince(t0);
        int info = -1;
        cudaMemcpy(&info, dInfo, sizeof(int), cudaMemcpyDeviceToHost);
        if (st == CUSOLVER_STATUS_SUCCESS && info == 0) ok[r] = 1;
        if (r == 0) cudaMemcpy(eig.data(), dD, sizeof(T) * n, cudaMemcpyDeviceToHost);

        cudaFree(dA); cudaFree(dQ); cudaFree(dD); cudaFree(dInfo); cudaFree(dW);
        cusolverMpDestroyMatrixDesc(descA);
        cusolverMpDestroyMatrixDesc(descQ);
        cusolverMpDestroyGrid(grid);
        cusolverMpDestroy(h);
        cudaStreamDestroy(stream);
    };
    std::vector<std::thread> threads;
    for (int r = 0; r < ndev; ++r) threads.emplace_back(worker, r);
    for (auto& t : threads) t.join();
    for (auto c : comms) ncclCommDestroy(c);

    double tr = 0.0, se = 0.0;
    for (int i = 0; i < n; ++i) { tr += host[i + static_cast<size_t>(i) * n]; se += eig[i]; }
    std::ostringstream mem;
    for (int r = 0; r < ndev; ++r) mem << (r ? "," : "") << static_cast<int>(used_gb[r] * 1024);
    std::printf("MP   ndev=%d n=%d %s nb=%d solve=%.2f s scatter=%.2f s ok=%d  used_MiB=[%s]  trace=%.6f sum_eig=%.6f\n",
                ndev, n, sizeof(T) == 8 ? "fp64" : "fp32", nb,
                *std::max_element(t_solve.begin(), t_solve.end()), *std::max_element(t_scatter.begin(), t_scatter.end()),
                static_cast<int>(std::count(ok.begin(), ok.end(), 1) == ndev), mem.str().c_str(), tr, se);
    return 0;
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
    for (std::string tok; std::getline(ss, tok, ',');) devs.push_back(std::atoi(tok.c_str()));
    const int nb = argc > 4 ? std::atoi(argv[4]) : 1024;
    return prec == "fp32" ? run<float>(n, devs, nb, CUDA_R_32F) : run<double>(n, devs, nb, CUDA_R_64F);
}
