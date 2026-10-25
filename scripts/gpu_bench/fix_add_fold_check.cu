// Claude Generated (Oct 2026): checks that the folded fix_add of gfnff_kernels.cu / gfix_add of xtb_gpu_context.cu
// gives the accumulator bits of the former 64-chunk loop, and that NaN / Inf / out-of-range sets the flag without a hang.
// Build and run: nvcc -O2 -arch=sm_86 -o fix_add_check scripts/gpu_bench/fix_add_fold_check.cu && ./fix_add_check
// Expected: '0 differ' for the three scales, fault injection 'bad flag 1', accumulators 0 0 0 0 0 10000000000000.
// Claude Generated: old chunk loop against the folded fix_add, accumulator bits compared on the device.
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <random>
__device__ int bad;
__device__ void old_add(unsigned long long* a, double v, double scale) {
    while (fabs(v) >= 64.0) { const double c = copysign(64.0, v); atomicAdd(a, (unsigned long long)llrint(c * scale)); v -= c; }
    atomicAdd(a, (unsigned long long)llrint(v * scale));
}
__device__ void new_add(unsigned long long* a, double v, double scale) {
    const double av = fabs(v);
    if (!(av * scale < 9.2e18)) { atomicOr(&bad, 1); if (!(av < 5.7e17)) return; }
    if (av < 64.0) { atomicAdd(a, (unsigned long long)llrint(v * scale)); return; }
    const double n = trunc(av * (1.0 / 64.0));
    const double r = v - copysign(n * 64.0, v);
    const unsigned long long chunk = (unsigned long long)llrint(copysign(64.0, v) * scale);
    atomicAdd(a, (unsigned long long)n * chunk + (unsigned long long)llrint(r * scale));
}
__device__ bool OLD_OK(double v) { return isfinite(v) && fabs(v) < 2e6; }
__global__ void kn(int n, const double* v, double scale, unsigned long long* nw) {
    int i = blockIdx.x * blockDim.x + threadIdx.x; if (i < n) new_add(&nw[i], v[i], scale);
}
__global__ void k(int n, const double* v, double scale, unsigned long long* o, unsigned long long* nw) {
    int i = blockIdx.x * blockDim.x + threadIdx.x; if (i >= n) return;
    o[i] = 0; nw[i] = 0;               // global memory accumulators, one per value
    if (OLD_OK(v[i])) old_add(&o[i], v[i], scale);
    new_add(&nw[i], v[i], scale);
}
int main() {
    std::vector<double> v;
    const double sp[] = {0, 1e-300, 0.5, 63.999999999999993, 64.0, 64.00000000000001, 127.99999999999999, 128.0, 128.5, 1000.1, 2047.9, 2048.0, 3000.7, 1e5, 123456.789, 1e8, 1e12};
    for (double x : sp) { v.push_back(x); v.push_back(-x); }
    std::mt19937_64 g(7); std::uniform_real_distribution<double> u(-1, 1);
    for (int e = -10; e <= 40; ++e) for (int j = 0; j < 4000; ++j) v.push_back(u(g) * std::pow(2.0, e));
    const int n = v.size(); double* dv; unsigned long long *o, *w;
    cudaMalloc(&dv, n * 8); cudaMalloc(&o, n * 8); cudaMalloc(&w, n * 8); cudaMemcpy(dv, v.data(), n * 8, cudaMemcpyHostToDevice);
    const double scales[] = {4503599627370496.0 /*2^52*/, 72057594037927936.0 /*2^56*/, 281474976710656.0 /*2^48*/};
    for (double sc : scales) {
        int zero = 0; cudaMemcpyToSymbol(bad, &zero, 4);
        k<<<(n + 255) / 256, 256>>>(n, dv, sc, o, w); cudaError_t ke = cudaDeviceSynchronize(); if (ke != cudaSuccess) { printf("kernel error %s\n", cudaGetErrorString(ke)); return 1; }
        std::vector<unsigned long long> ho(n), hw(n); cudaMemcpy(ho.data(), o, n * 8, cudaMemcpyDeviceToHost); cudaMemcpy(hw.data(), w, n * 8, cudaMemcpyDeviceToHost);
        int hb; cudaMemcpyFromSymbol(&hb, bad, 4);
        long diff = 0, inrange = 0, nz = 0; for (int i = 0; i < n; ++i) { if (std::fabs(v[i]) * sc < 9.2e18 && std::fabs(v[i]) < 2e6) { ++inrange; if (ho[i] != hw[i]) ++diff; if (ho[i] != 0) ++nz; } }
        printf("scale 2^%d: %d values, %ld in range, %ld differ, %ld nonzero (old loop vs fold), bad flag %d\n", (int)std::lround(std::log2(sc)), n, inrange, diff, nz, hb);
    }
    // fault injection: new path only (the old loop never ends for Inf)
    double bad_in[] = {NAN, INFINITY, -INFINITY, 1e30, 3e17, 1.0};
    cudaMemcpy(dv, bad_in, sizeof(bad_in), cudaMemcpyHostToDevice);
    int zero = 0; cudaMemcpyToSymbol(bad, &zero, 4);
    cudaMemset(w, 0, 6 * 8);
    kn<<<1, 6>>>(6, dv, 4503599627370496.0, w); cudaError_t e = cudaDeviceSynchronize();
    int hb = -1; cudaMemcpyFromSymbol(&hb, bad, 4);
    std::vector<unsigned long long> hw(6); cudaMemcpy(hw.data(), w, 48, cudaMemcpyDeviceToHost);
    printf("fault injection: cuda %s, bad flag %d, accumulators:", cudaGetErrorString(e), hb);
    for (int i = 0; i < 6; ++i) printf(" %llx", hw[i]);
    printf("\n");
    return 0;
}
