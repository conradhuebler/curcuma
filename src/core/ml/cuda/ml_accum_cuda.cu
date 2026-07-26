/*
 * < CUDA accumulation of the Delta-ML normal equations >
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
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
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Claude Generated (July 2026): GPU path of the Delta-ML work package.
 *
 * Why this and not the descriptor: the fit spends its time in one place. Per training
 * structure the design block Phi is (1 + 3N) x M, and folding it into the normal equations
 * costs M^2 * 3N -- 2.2 GFlop at M = 2000 and N = 90, against roughly 0.16 GFlop for
 * building Phi itself. Moving that single symmetric rank-k update to the device therefore
 * captures most of the training cost, and it maps onto a stock cuBLAS call rather than a
 * hand-written kernel.
 *
 * Storage convention: the host hands over Phi in ROW-major (rows x M). Read with leading
 * dimension M, the very same buffer is the COLUMN-major matrix Phi^T of shape (M x rows),
 * which is exactly the operand DSYRK and DGEMV need. No transpose or repacking is required
 * anywhere -- the two conventions cancel.
 */

#include <cublas_v2.h>
#include <cuda_runtime.h>

#include <cstdio>
#include <cstdlib>
#include <vector>

namespace {

struct MLAccumContext {
    cublasHandle_t blas = nullptr;
    int M = 0;
    double* d_A = nullptr;   ///< M x M, lower triangle (column major)
    double* d_b = nullptr;   ///< M
    double* d_Phi = nullptr; ///< capacity_rows x M
    double* d_y = nullptr;   ///< capacity_rows
    int capacity_rows = 0;
};

/*! \brief Grow the staging buffers so a block of @p rows fits. */
bool ensureCapacity(MLAccumContext* ctx, int rows)
{
    if (rows <= ctx->capacity_rows)
        return true;
    if (ctx->d_Phi)
        cudaFree(ctx->d_Phi);
    if (ctx->d_y)
        cudaFree(ctx->d_y);
    ctx->d_Phi = nullptr;
    ctx->d_y = nullptr;

    const size_t phi_bytes = static_cast<size_t>(rows) * ctx->M * sizeof(double);
    if (cudaMalloc(&ctx->d_Phi, phi_bytes) != cudaSuccess)
        return false;
    if (cudaMalloc(&ctx->d_y, static_cast<size_t>(rows) * sizeof(double)) != cudaSuccess)
        return false;
    ctx->capacity_rows = rows;
    return true;
}

} // namespace

extern "C" {

/*! \brief Allocate the device-side accumulator for an M x M system. Returns nullptr on failure. */
void* curcuma_cuda_ml_accum_create(int M)
{
    if (M <= 0)
        return nullptr;

    MLAccumContext* ctx = new MLAccumContext();
    ctx->M = M;

    if (cublasCreate(&ctx->blas) != CUBLAS_STATUS_SUCCESS) {
        delete ctx;
        return nullptr;
    }

    const size_t a_bytes = static_cast<size_t>(M) * M * sizeof(double);
    if (cudaMalloc(&ctx->d_A, a_bytes) != cudaSuccess
        || cudaMalloc(&ctx->d_b, static_cast<size_t>(M) * sizeof(double)) != cudaSuccess) {
        cublasDestroy(ctx->blas);
        if (ctx->d_A)
            cudaFree(ctx->d_A);
        delete ctx;
        return nullptr;
    }

    cudaMemset(ctx->d_A, 0, a_bytes);
    cudaMemset(ctx->d_b, 0, static_cast<size_t>(M) * sizeof(double));
    return ctx;
}

/**
 * @brief Fold one block of observations into the normal equations.
 *
 *   A += alpha * Phi^T Phi      (symmetric rank-k update, lower triangle)
 *   b += alpha * Phi^T y
 *
 * @param Phi row-major, rows x M
 * @param y   rows
 * @return 0 on success
 */
int curcuma_cuda_ml_accum_add(void* handle, int rows, const double* Phi, const double* y, double alpha)
{
    MLAccumContext* ctx = static_cast<MLAccumContext*>(handle);
    if (!ctx || rows <= 0 || !Phi || !y)
        return 1;
    if (!ensureCapacity(ctx, rows))
        return 2;

    const size_t phi_bytes = static_cast<size_t>(rows) * ctx->M * sizeof(double);
    if (cudaMemcpy(ctx->d_Phi, Phi, phi_bytes, cudaMemcpyHostToDevice) != cudaSuccess)
        return 3;
    if (cudaMemcpy(ctx->d_y, y, static_cast<size_t>(rows) * sizeof(double), cudaMemcpyHostToDevice)
        != cudaSuccess)
        return 4;

    const double one = 1.0;

    // Row-major (rows x M) read as column-major with ld = M is Phi^T (M x rows), so the
    // update C := alpha * A * A^T + C with A = Phi^T gives exactly alpha * Phi^T Phi.
    if (cublasDsyrk(ctx->blas, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_N, ctx->M, rows, &alpha,
            ctx->d_Phi, ctx->M, &one, ctx->d_A, ctx->M)
        != CUBLAS_STATUS_SUCCESS)
        return 5;

    if (cublasDgemv(ctx->blas, CUBLAS_OP_N, ctx->M, rows, &alpha, ctx->d_Phi, ctx->M, ctx->d_y, 1,
            &one, ctx->d_b, 1)
        != CUBLAS_STATUS_SUCCESS)
        return 6;

    return 0;
}

/**
 * @brief Copy the accumulated system back to the host.
 * @param A out, M x M; only the lower triangle in column-major order is filled, which the
 *          caller mirrors -- symmetrising here would just move the same work.
 * @param b out, M
 */
int curcuma_cuda_ml_accum_fetch(void* handle, double* A, double* b)
{
    MLAccumContext* ctx = static_cast<MLAccumContext*>(handle);
    if (!ctx || !A || !b)
        return 1;
    if (cudaDeviceSynchronize() != cudaSuccess)
        return 2;
    const size_t a_bytes = static_cast<size_t>(ctx->M) * ctx->M * sizeof(double);
    if (cudaMemcpy(A, ctx->d_A, a_bytes, cudaMemcpyDeviceToHost) != cudaSuccess)
        return 3;
    if (cudaMemcpy(b, ctx->d_b, static_cast<size_t>(ctx->M) * sizeof(double), cudaMemcpyDeviceToHost)
        != cudaSuccess)
        return 4;
    return 0;
}

void curcuma_cuda_ml_accum_destroy(void* handle)
{
    MLAccumContext* ctx = static_cast<MLAccumContext*>(handle);
    if (!ctx)
        return;
    if (ctx->d_A)
        cudaFree(ctx->d_A);
    if (ctx->d_b)
        cudaFree(ctx->d_b);
    if (ctx->d_Phi)
        cudaFree(ctx->d_Phi);
    if (ctx->d_y)
        cudaFree(ctx->d_y);
    if (ctx->blas)
        cublasDestroy(ctx->blas);
    delete ctx;
}

} // extern "C"

// =================================================================================
// Inference
// =================================================================================
//
// The prediction is a different shape of problem from the fit: only ~22 MFlop for 90 atoms
// at M = 2000, but repeated every dynamics step. What decides whether the GPU helps is
// therefore traffic, not arithmetic.
//
// So the centres, weights and feature scaling stay RESIDENT on the device (about 1 MB, set
// up once when the model is loaded), and the compact Jacobian blocks -- roughly 5 MB per
// step, by far the largest object -- never cross the bus at all: the device returns the
// per-atom adjoint in feature space (90 x 60, some 43 KB) and the host contracts it with
// the Jacobian it already holds. Traffic per step is thus two small buffers instead of
// megabytes.

namespace {

struct MLInferContext {
    int M = 0;
    int F = 0;
    double sigma = 1.0;
    double* d_centers = nullptr;   ///< M x F, row major, standardised
    double* d_beta = nullptr;      ///< M
    double* d_inv_scale = nullptr; ///< F
    double* d_X = nullptr;         ///< natoms x F (capacity)
    double* d_k = nullptr;         ///< natoms x M (capacity)
    double* d_s = nullptr;         ///< natoms
    double* d_kmax = nullptr;      ///< natoms
    double* d_W = nullptr;         ///< natoms x F (capacity)
    int* d_start = nullptr;
    int* d_count = nullptr;
    int capacity_atoms = 0;
};

/*! \brief k[i,m] = exp(-|x_i - z_m|^2 / 2 sigma^2) inside atom i's element block, else 0. */
__global__ void kernelK(int natoms, int M, int F, double inv_two_sigma2, const double* __restrict__ X,
    const double* __restrict__ Z, const int* __restrict__ start, const int* __restrict__ count,
    double* __restrict__ k)
{
    const int m = blockIdx.x * blockDim.x + threadIdx.x;
    const int i = blockIdx.y;
    if (i >= natoms || m >= M)
        return;

    const int s0 = start[i];
    const int c0 = count[i];
    if (c0 <= 0 || m < s0 || m >= s0 + c0) {
        k[static_cast<size_t>(i) * M + m] = 0.0;
        return;
    }

    double d2 = 0.0;
    for (int f = 0; f < F; ++f) {
        const double diff = X[static_cast<size_t>(i) * F + f] - Z[static_cast<size_t>(m) * F + f];
        d2 += diff * diff;
    }
    k[static_cast<size_t>(i) * M + m] = exp(-d2 * inv_two_sigma2);
}

/*! \brief Per atom: s_i = sum_m beta_m k_im and the largest k, which is the indicator. */
__global__ void kernelReduce(int natoms, int M, const double* __restrict__ k,
    const double* __restrict__ beta, double* __restrict__ s, double* __restrict__ kmax)
{
    extern __shared__ double shared[];
    double* sh_s = shared;
    double* sh_m = shared + blockDim.x;

    const int i = blockIdx.x;
    double acc = 0.0;
    double best = 0.0;
    for (int m = threadIdx.x; m < M; m += blockDim.x) {
        const double kv = k[static_cast<size_t>(i) * M + m];
        acc += beta[m] * kv;
        best = fmax(best, kv);
    }
    sh_s[threadIdx.x] = acc;
    sh_m[threadIdx.x] = best;
    __syncthreads();

    for (int stride = blockDim.x / 2; stride > 0; stride >>= 1) {
        if (threadIdx.x < stride) {
            sh_s[threadIdx.x] += sh_s[threadIdx.x + stride];
            sh_m[threadIdx.x] = fmax(sh_m[threadIdx.x], sh_m[threadIdx.x + stride]);
        }
        __syncthreads();
    }
    if (threadIdx.x == 0) {
        s[i] = sh_s[0];
        kmax[i] = sh_m[0];
    }
}

/*! \brief Adjoint in raw descriptor space: W_if = (s_i x_if - sum_m beta_m k_im z_mf) * (-1/sigma^2) / scale_f. */
__global__ void kernelW(int natoms, int M, int F, double neg_inv_sigma2,
    const double* __restrict__ X, const double* __restrict__ Z, const double* __restrict__ k,
    const double* __restrict__ beta, const double* __restrict__ s,
    const double* __restrict__ inv_scale, const int* __restrict__ start,
    const int* __restrict__ count, double* __restrict__ W)
{
    const int f = blockIdx.x * blockDim.x + threadIdx.x;
    const int i = blockIdx.y;
    if (i >= natoms || f >= F)
        return;

    const int s0 = start[i];
    const int c0 = count[i];
    double dot = 0.0;
    for (int m = s0; m < s0 + c0; ++m)
        dot += beta[m] * k[static_cast<size_t>(i) * M + m] * Z[static_cast<size_t>(m) * F + f];

    const double v = s[i] * X[static_cast<size_t>(i) * F + f] - dot;
    W[static_cast<size_t>(i) * F + f] = v * neg_inv_sigma2 * inv_scale[f];
}

bool ensureInferCapacity(MLInferContext* c, int natoms)
{
    if (natoms <= c->capacity_atoms)
        return true;
    if (c->d_X) cudaFree(c->d_X);
    if (c->d_k) cudaFree(c->d_k);
    if (c->d_s) cudaFree(c->d_s);
    if (c->d_kmax) cudaFree(c->d_kmax);
    if (c->d_W) cudaFree(c->d_W);
    if (c->d_start) cudaFree(c->d_start);
    if (c->d_count) cudaFree(c->d_count);
    c->d_X = c->d_k = c->d_s = c->d_kmax = c->d_W = nullptr;
    c->d_start = c->d_count = nullptr;

    const size_t n = static_cast<size_t>(natoms);
    if (cudaMalloc(&c->d_X, n * c->F * sizeof(double)) != cudaSuccess) return false;
    if (cudaMalloc(&c->d_k, n * c->M * sizeof(double)) != cudaSuccess) return false;
    if (cudaMalloc(&c->d_s, n * sizeof(double)) != cudaSuccess) return false;
    if (cudaMalloc(&c->d_kmax, n * sizeof(double)) != cudaSuccess) return false;
    if (cudaMalloc(&c->d_W, n * c->F * sizeof(double)) != cudaSuccess) return false;
    if (cudaMalloc(&c->d_start, n * sizeof(int)) != cudaSuccess) return false;
    if (cudaMalloc(&c->d_count, n * sizeof(int)) != cudaSuccess) return false;
    c->capacity_atoms = natoms;
    return true;
}

} // namespace

extern "C" {

void curcuma_cuda_ml_infer_destroy(void* handle);

void* curcuma_cuda_ml_infer_create(int M, int F, const double* centers, const double* beta,
    const double* inv_scale, double sigma)
{
    if (M <= 0 || F <= 0 || !centers || !beta || !inv_scale)
        return nullptr;
    MLInferContext* c = new MLInferContext();
    c->M = M;
    c->F = F;
    c->sigma = sigma;

    const size_t zb = static_cast<size_t>(M) * F * sizeof(double);
    if (cudaMalloc(&c->d_centers, zb) != cudaSuccess
        || cudaMalloc(&c->d_beta, static_cast<size_t>(M) * sizeof(double)) != cudaSuccess
        || cudaMalloc(&c->d_inv_scale, static_cast<size_t>(F) * sizeof(double)) != cudaSuccess) {
        curcuma_cuda_ml_infer_destroy(c);
        return nullptr;
    }
    cudaMemcpy(c->d_centers, centers, zb, cudaMemcpyHostToDevice);
    cudaMemcpy(c->d_beta, beta, static_cast<size_t>(M) * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(c->d_inv_scale, inv_scale, static_cast<size_t>(F) * sizeof(double), cudaMemcpyHostToDevice);
    return c;
}

/**
 * @brief One prediction pass.
 * @param X       natoms x F, row major, already standardised
 * @param start   per atom, first centre of its element block
 * @param count   per atom, number of centres (0 for an element the model does not know)
 * @param energy  out, the summed kernel contribution (no element baseline)
 * @param W       out, natoms x F adjoint in raw descriptor space
 * @param uncertainty out, max over atoms of 1 - max_m k
 */
int curcuma_cuda_ml_infer_run(void* handle, int natoms, const double* X, const int* start,
    const int* count, double* energy, double* W, double* uncertainty)
{
    MLInferContext* c = static_cast<MLInferContext*>(handle);
    if (!c || natoms <= 0 || !X || !start || !count || !energy || !W)
        return 1;
    if (!ensureInferCapacity(c, natoms))
        return 2;

    const size_t n = static_cast<size_t>(natoms);
    cudaMemcpy(c->d_X, X, n * c->F * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(c->d_start, start, n * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(c->d_count, count, n * sizeof(int), cudaMemcpyHostToDevice);

    const double inv_two_sigma2 = 1.0 / (2.0 * c->sigma * c->sigma);
    dim3 blockK(128), gridK((c->M + 127) / 128, natoms);
    kernelK<<<gridK, blockK>>>(natoms, c->M, c->F, inv_two_sigma2, c->d_X, c->d_centers,
        c->d_start, c->d_count, c->d_k);

    const int red_threads = 256;
    kernelReduce<<<natoms, red_threads, 2 * red_threads * sizeof(double)>>>(
        natoms, c->M, c->d_k, c->d_beta, c->d_s, c->d_kmax);

    dim3 blockW(64), gridW((c->F + 63) / 64, natoms);
    kernelW<<<gridW, blockW>>>(natoms, c->M, c->F, -1.0 / (c->sigma * c->sigma), c->d_X,
        c->d_centers, c->d_k, c->d_beta, c->d_s, c->d_inv_scale, c->d_start, c->d_count, c->d_W);

    if (cudaDeviceSynchronize() != cudaSuccess)
        return 3;

    std::vector<double> s_host(n), kmax_host(n);
    cudaMemcpy(s_host.data(), c->d_s, n * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(kmax_host.data(), c->d_kmax, n * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(W, c->d_W, n * c->F * sizeof(double), cudaMemcpyDeviceToHost);

    double e = 0.0, worst = 0.0;
    for (int i = 0; i < natoms; ++i) {
        e += s_host[i];
        worst = fmax(worst, (count[i] > 0) ? 1.0 - kmax_host[i] : 1.0);
    }
    *energy = e;
    if (uncertainty)
        *uncertainty = worst;
    return 0;
}

void curcuma_cuda_ml_infer_destroy(void* handle)
{
    MLInferContext* c = static_cast<MLInferContext*>(handle);
    if (!c)
        return;
    if (c->d_centers) cudaFree(c->d_centers);
    if (c->d_beta) cudaFree(c->d_beta);
    if (c->d_inv_scale) cudaFree(c->d_inv_scale);
    if (c->d_X) cudaFree(c->d_X);
    if (c->d_k) cudaFree(c->d_k);
    if (c->d_s) cudaFree(c->d_s);
    if (c->d_kmax) cudaFree(c->d_kmax);
    if (c->d_W) cudaFree(c->d_W);
    if (c->d_start) cudaFree(c->d_start);
    if (c->d_count) cudaFree(c->d_count);
    delete c;
}

} // extern "C"
