// Agent N: bitwise run-to-run check of single GPU routines (no curcuma). Each iteration recomputes the same
// result from the same input; counts iterations whose output bits differ from iteration 0.
// usage: bitcheck <n> <iters> <routine>   routine: sygvd | potrf | syevd | ssyevd | dgemm | kernel
#include <hip/hip_runtime.h>
#include <rocblas/rocblas.h>
#include <rocsolver/rocsolver.h>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <cmath>
#include <algorithm>

__global__ void k_dot(int n, const double* a, const double* b, double* out)
{   // one thread per row, sequential, deterministic
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    double s = 0; for (int k = 0; k < n; ++k) s = fma(a[i + (size_t)k * n], b[k], s);
    out[i] = s;
}

int main(int argc, char** argv)
{
    const int n = atoi(argv[1]), iters = atoi(argv[2]); std::string r = argv[3];
    const size_t nn = (size_t)n * n;
    srand(7);
    std::vector<double> A(nn), B(nn), M(nn);
    for (auto& x : M) x = (rand() / (double)RAND_MAX - 0.5);
    for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) { double a = 0; for (int k = 0; k < n; ++k) a += M[i + k * n] * M[j + k * n]; B[i + j * n] = a / n + (i == j ? 1.0 : 0.0); }
    for (auto& x : M) x = (rand() / (double)RAND_MAX - 0.5);
    for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) A[i + j * n] = 0.5 * (M[i + j * n] + M[j + i * n]);
    std::vector<float> Af(nn); for (size_t i = 0; i < nn; ++i) Af[i] = (float)A[i];
    rocblas_handle h; rocblas_create_handle(&h);
    double *dA, *dB, *dC, *dW, *dE; float *dAf, *dWf, *dEf; rocblas_int* dI;
    hipMalloc(&dA, nn * 8); hipMalloc(&dB, nn * 8); hipMalloc(&dC, nn * 8); hipMalloc(&dW, n * 8); hipMalloc(&dE, n * 8);
    hipMalloc(&dAf, nn * 4); hipMalloc(&dWf, n * 4); hipMalloc(&dEf, n * 4); hipMalloc(&dI, 4);
    std::vector<unsigned char> ref, cur;
    int bad = 0, firstbad = -1; double maxdiff = 0;
    for (int it = 0; it < iters; ++it) {
        hipMemcpy(dA, A.data(), nn * 8, hipMemcpyHostToDevice);
        hipMemcpy(dB, B.data(), nn * 8, hipMemcpyHostToDevice);
        size_t bytes = 0; const void* src = nullptr;
        if (r == "sygvd") {
            rocsolver_dsygvd(h, rocblas_eform_ax, rocblas_evect_original, rocblas_fill_lower, n, dA, n, dB, n, dW, dE, dI);
            src = dA; bytes = nn * 8;
        } else if (r == "potrf") {
            rocsolver_dpotrf(h, rocblas_fill_lower, n, dB, n, dI); src = dB; bytes = nn * 8;
        } else if (r == "syevd") {
            rocsolver_dsyevd(h, rocblas_evect_original, rocblas_fill_lower, n, dA, n, dW, dE, dI); src = dA; bytes = nn * 8;
        } else if (r == "ssyevd") {
            hipMemcpy(dAf, Af.data(), nn * 4, hipMemcpyHostToDevice);
            rocsolver_ssyevd(h, rocblas_evect_original, rocblas_fill_lower, n, dAf, n, dWf, dEf, dI); src = dAf; bytes = nn * 4;
        } else if (r == "ssygvd") {
            std::vector<float> Bf(B.begin(), B.end());
            float* dBf = nullptr; hipMalloc(&dBf, nn * 4);
            hipMemcpy(dAf, Af.data(), nn * 4, hipMemcpyHostToDevice); hipMemcpy(dBf, Bf.data(), nn * 4, hipMemcpyHostToDevice);
            rocsolver_ssygvd(h, rocblas_eform_ax, rocblas_evect_original, rocblas_fill_lower, n, dAf, n, dBf, n, dWf, dEf, dI);
            hipDeviceSynchronize(); hipFree(dBf); src = dAf; bytes = nn * 4;
        } else if (r == "dgemm") {
            const double one = 1.0, zero = 0.0;
            rocblas_dgemm(h, rocblas_operation_none, rocblas_operation_transpose, n, n, n, &one, dA, n, dB, n, &zero, dC, n);
            src = dC; bytes = nn * 8;
        } else {
            k_dot<<<(n + 127) / 128, 128>>>(n, dA, dB, dC); src = dC; bytes = n * 8;
        }
        hipDeviceSynchronize();
        cur.resize(bytes); hipMemcpy(cur.data(), src, bytes, hipMemcpyDeviceToHost);
        if (it == 0) ref = cur;
        else if (std::memcmp(ref.data(), cur.data(), bytes) != 0) {
            ++bad; if (firstbad < 0) firstbad = it;
            double md = 0, mr = 0;
            const bool f32 = (r == "ssyevd" || r == "ssygvd");
            const size_t cnt = bytes / (f32 ? 4 : 8);
            for (size_t k = 0; k < cnt; ++k) {
                double a = f32 ? ((const float*)ref.data())[k] : ((const double*)ref.data())[k];
                double b = f32 ? ((const float*)cur.data())[k] : ((const double*)cur.data())[k];
                if (r == "sygvd" || r == "syevd" || r == "ssyevd" || r == "ssygvd") { a = std::fabs(a); b = std::fabs(b); } // eigenvector sign
                md = std::max(md, std::fabs(a - b)); mr = std::max(mr, std::fabs(a));
            }
            if (md > maxdiff) maxdiff = md;
        }
    }
    std::printf("%s n %d iters %d differ_from_first %d first_at %d max_abs_diff %.3e\n", r.c_str(), n, iters, bad, firstbad, maxdiff);
    return 0;
}
