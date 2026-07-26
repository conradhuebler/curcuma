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
