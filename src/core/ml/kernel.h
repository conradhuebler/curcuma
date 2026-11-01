/*
 * < Kernel functions for Delta-Machine-Learning >
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
 * Claude Generated (July 2026): kernel layer of the Delta-ML work package.
 */

#pragma once

#include "src/core/global.h"

#include <cmath>

namespace curcuma {
namespace ml {

/**
 * @brief Squared Euclidean distances between all rows of A and all rows of B.
 *
 * Uses the expansion |a-b|^2 = |a|^2 + |b|^2 - 2 a.b so the bulk of the work is a single
 * matrix product. That is not just an optimisation: it is the reason this routine maps
 * directly onto a GEMM call when the kernel evaluation is later moved to the GPU, where
 * curcuma already has cuBLAS/rocBLAS available through the plugin layer.
 *
 * @param A n x F matrix of descriptors
 * @param B m x F matrix of descriptors
 * @return n x m matrix of squared distances, clamped at zero
 */
inline Matrix squaredDistances(const Matrix& A, const Matrix& B)
{
    const Vector a2 = A.rowwise().squaredNorm();
    const Vector b2 = B.rowwise().squaredNorm();

    Matrix d2 = A * B.transpose();
    d2 *= -2.0;
    d2.colwise() += a2;
    d2.rowwise() += b2.transpose();

    // Cancellation can push near-zero distances slightly negative.
    return d2.cwiseMax(0.0);
}

/**
 * @brief Gaussian (squared-exponential) kernel from precomputed squared distances.
 *
 *   k(x,y) = exp(-|x-y|^2 / (2 sigma^2))
 *
 * The Gaussian choice matters for the failure mode of delta-learning: far outside the
 * training region k goes to zero, so both the correction and its gradient vanish smoothly
 * and the model degrades back to the pure base method instead of extrapolating wildly.
 */
inline Matrix gaussianKernel(const Matrix& d2, double sigma)
{
    const double inv = -1.0 / (2.0 * sigma * sigma);
    return (d2 * inv).array().exp().matrix();
}

/*! \brief Convenience overload computing distances and kernel in one step. */
inline Matrix gaussianKernel(const Matrix& A, const Matrix& B, double sigma)
{
    return gaussianKernel(squaredDistances(A, B), sigma);
}

} // namespace ml
} // namespace curcuma
