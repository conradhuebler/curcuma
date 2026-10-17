// Copyright (C) 2026 Conrad Huebler <Conrad.Huebler@gmx.net>
//
// Claude Generated (Oct 2026): reproducer for sporadically wrong rocSOLVER results when several
// processes share one GPU. Not part of the curcuma build.
//
// Each iteration uploads a fixed symmetric matrix A and a fixed positive definite matrix B
// (n = 66, the size of the caffeine GFN2 problem), runs rocsolver_dpotrf on B and
// rocsolver_dsygvd on (A, B), and checks the result on the host:
//   potrf: max |L L^T - B| <= 1e-10
//   sygvd: max over eigenpairs of |A c - w B c| <= 1e-8
// and prints how many iterations failed. Measured on RX 9060 XT (gfx1200), ROCm 7.2.4,
// Linux 6.18.49: one process 0 of 12000 iterations wrong for both routines, four processes at
// the same time potrf 6 and sygvd 8 wrong of 12000 each. Details: docs/ROCM_GFN_POTENTIAL_2026-10-14.md.
//
// Build:  hipcc -O2 -o rocsolver_concurrency_check scripts/rocsolver_concurrency_check.cpp -lrocsolver -lrocblas
// Run:    ./rocsolver_concurrency_check 12000                      (one process, control)
//         for i in 1 2 3 4; do ./rocsolver_concurrency_check 3000 & done; wait
#include <hip/hip_runtime.h>
#include <rocblas/rocblas.h>
#include <rocsolver/rocsolver.h>
#include <vector>
#include <cstdio>
#include <cmath>
#include <cstdlib>
int main(int argc, char** argv) {
    const int n = 66, iters = argc > 1 ? atoi(argv[1]) : 2000; const size_t nn = (size_t)n * n;
    srand(7);
    std::vector<double> A(nn), B(nn), M(nn);
    for (auto& x : M) x = (rand() / (double)RAND_MAX - 0.5);
    for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) { double a = 0; for (int k = 0; k < n; ++k) a += M[i + k * n] * M[j + k * n]; B[i + j * n] = a / n + (i == j ? 1.0 : 0.0); }
    for (auto& x : M) x = (rand() / (double)RAND_MAX - 0.5);
    for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) A[i + j * n] = 0.5 * (M[i + j * n] + M[j + i * n]);
    rocblas_handle h; rocblas_create_handle(&h);
    double *dA, *dB, *dW, *dE; rocblas_int* dI; hipMalloc(&dA, nn * 8); hipMalloc(&dB, nn * 8); hipMalloc(&dW, n * 8); hipMalloc(&dE, n * 8); hipMalloc(&dI, 4);
    int badp = 0, bads = 0; std::vector<double> L(nn), C(nn), w(n);
    for (int it = 0; it < iters; ++it) {
        hipMemcpy(dB, B.data(), nn * 8, hipMemcpyHostToDevice);
        rocsolver_dpotrf(h, rocblas_fill_lower, n, dB, n, dI); hipDeviceSynchronize();
        hipMemcpy(L.data(), dB, nn * 8, hipMemcpyDeviceToHost);
        double e = 0; for (int i = 0; i < n; ++i) for (int j = 0; j <= i; ++j) { double a = 0; for (int k = 0; k <= j; ++k) a += L[i + k * n] * L[j + k * n]; e = fmax(e, fabs(a - B[i + j * n])); }
        if (e > 1e-10) ++badp;
        hipMemcpy(dA, A.data(), nn * 8, hipMemcpyHostToDevice); hipMemcpy(dB, B.data(), nn * 8, hipMemcpyHostToDevice);
        rocsolver_dsygvd(h, rocblas_eform_ax, rocblas_evect_original, rocblas_fill_lower, n, dA, n, dB, n, dW, dE, dI); hipDeviceSynchronize();
        hipMemcpy(C.data(), dA, nn * 8, hipMemcpyDeviceToHost); hipMemcpy(w.data(), dW, n * 8, hipMemcpyDeviceToHost);
        double r = 0; for (int k = 0; k < n; ++k) for (int i = 0; i < n; ++i) { double fx = 0, sx = 0; for (int j = 0; j < n; ++j) { double aij = i >= j ? A[i + j * n] : A[j + i * n]; double bij = i >= j ? B[i + j * n] : B[j + i * n]; fx += aij * C[j + k * n]; sx += bij * C[j + k * n]; } r = fmax(r, fabs(fx - w[k] * sx)); }
        if (r > 1e-8) ++bads;
    }
    printf("iters %d bad_potrf %d bad_sygvd %d\n", iters, badp, bads);
}
