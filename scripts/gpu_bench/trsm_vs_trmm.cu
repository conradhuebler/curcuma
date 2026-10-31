// Build and run: nvcc -O2 -arch=sm_86 -o trsm_vs_trmm scripts/gpu_bench/trsm_vs_trmm.cu -lcublas && ./trsm_vs_trmm
// Claude Generated: cublas trsm (two per Fock reduction, as the cached-L path) against trmm (explicit L^-1), FP64 and FP32.
#include <cstdio>
#include <vector>
#include <cublas_v2.h>
#include <cuda_runtime.h>
template <typename T> struct Ops;
template <> struct Ops<double> {
    static cublasStatus_t trsm(cublasHandle_t h, cublasSideMode_t s, cublasOperation_t o, int n, const double* a, double* b) { double one = 1; return cublasDtrsm(h, s, CUBLAS_FILL_MODE_LOWER, o, CUBLAS_DIAG_NON_UNIT, n, n, &one, a, n, b, n); }
    static cublasStatus_t trmm(cublasHandle_t h, cublasSideMode_t s, cublasOperation_t o, int n, const double* a, double* b, double* c) { double one = 1; return cublasDtrmm(h, s, CUBLAS_FILL_MODE_LOWER, o, CUBLAS_DIAG_NON_UNIT, n, n, &one, a, n, b, n, c, n); }
};
template <> struct Ops<float> {
    static cublasStatus_t trsm(cublasHandle_t h, cublasSideMode_t s, cublasOperation_t o, int n, const float* a, float* b) { float one = 1; return cublasStrsm(h, s, CUBLAS_FILL_MODE_LOWER, o, CUBLAS_DIAG_NON_UNIT, n, n, &one, a, n, b, n); }
    static cublasStatus_t trmm(cublasHandle_t h, cublasSideMode_t s, cublasOperation_t o, int n, const float* a, float* b, float* c) { float one = 1; return cublasStrmm(h, s, CUBLAS_FILL_MODE_LOWER, o, CUBLAS_DIAG_NON_UNIT, n, n, &one, a, n, b, n, c, n); }
};
template <typename T> void run(int n, const char* name) {
    T *L, *A, *C; size_t bytes = (size_t)n * n * sizeof(T);
    cudaMalloc(&L, bytes); cudaMalloc(&A, bytes); cudaMalloc(&C, bytes);
    std::vector<T> h((size_t)n * n, T(0)); for (int i = 0; i < n; ++i) for (int j = 0; j <= i; ++j) h[(size_t)j * n + i] = (i == j) ? T(2) : T(1e-3);
    cudaMemcpy(L, h.data(), bytes, cudaMemcpyHostToDevice); cudaMemcpy(A, h.data(), bytes, cudaMemcpyHostToDevice);
    cublasHandle_t cb; cublasCreate(&cb); cudaEvent_t e0, e1; cudaEventCreate(&e0); cudaEventCreate(&e1);
    auto t = [&](auto f) { f(); cudaDeviceSynchronize(); cudaEventRecord(e0); for (int i = 0; i < 5; ++i) f(); cudaEventRecord(e1); cudaEventSynchronize(e1); float ms; cudaEventElapsedTime(&ms, e0, e1); return ms / 5; };
    float ts = t([&] { Ops<T>::trsm(cb, CUBLAS_SIDE_RIGHT, CUBLAS_OP_T, n, L, A); Ops<T>::trsm(cb, CUBLAS_SIDE_LEFT, CUBLAS_OP_N, n, L, A); });
    float tm = t([&] { Ops<T>::trmm(cb, CUBLAS_SIDE_LEFT, CUBLAS_OP_N, n, L, A, C); Ops<T>::trmm(cb, CUBLAS_SIDE_RIGHT, CUBLAS_OP_T, n, L, C, A); });
    printf("%s n=%5d: 2 x trsm %8.2f ms, 2 x trmm %8.2f ms, ratio trsm/trmm %.2f\n", name, n, ts, tm, ts / tm);
    cudaFree(L); cudaFree(A); cudaFree(C); cublasDestroy(cb);
}
int main() {
    cudaDeviceProp p; cudaGetDeviceProperties(&p, 0); printf("%s\n", p.name);
    for (int n : {1500, 3222, 6000}) { run<double>(n, "FP64"); run<float>(n, "FP32"); }
}
