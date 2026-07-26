/*
 * < GPU offload for the Delta-ML normal equations >
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
 */

#pragma once

#include "src/core/global.h"

#include <string>
#include <vector>

namespace curcuma {
namespace ml {

/**
 * @brief Device-side accumulator for (A += alpha Phi^T Phi, b += alpha Phi^T y).
 *
 * Thin RAII wrapper around the entry points in the GPU backend plugin. The plugin is
 * dlopen'd on demand through the same cache the method factory uses, so a CPU-only build
 * or a machine without a working driver simply reports unavailable and the caller keeps
 * its CPU path -- no CUDA symbols are linked into the main binary.
 *
 * Numerical note: cuBLAS sums in a different order than the CPU BLAS, so the accumulated
 * system agrees to roughly rounding level rather than bit-for-bit. That difference
 * propagates into the fitted weights, which is why the validation compares predictions
 * within a tolerance instead of asserting identity.
 */
class GpuAccumulator {
public:
    /*! \brief Create an accumulator for an M x M system; check valid() afterwards. */
    GpuAccumulator(const std::string& backend, int M);
    ~GpuAccumulator();

    GpuAccumulator(const GpuAccumulator&) = delete;
    GpuAccumulator& operator=(const GpuAccumulator&) = delete;

    /*! \brief True if the plugin loaded, the symbols resolved and the device context exists. */
    bool valid() const { return m_ctx != nullptr; }

    /**
     * @brief Fold one block of observations in.
     * @param rows number of observation rows in this block
     * @param phi  row-major, rows x M
     * @param y    observations, length rows
     */
    bool add(int rows, const double* phi, const double* y, double alpha);

    /*! \brief Retrieve the system, symmetrised into row-major @p A and @p b. */
    bool fetch(Matrix& A, Vector& b);

private:
    void* m_ctx = nullptr;
    int m_M = 0;

    void* (*m_create)(int) = nullptr;
    int (*m_add)(void*, int, const double*, const double*, double) = nullptr;
    int (*m_fetch)(void*, double*, double*) = nullptr;
    void (*m_destroy)(void*) = nullptr;
};

/**
 * @brief Device-resident evaluation of the learned correction.
 *
 * Centres, weights and feature scaling are uploaded once; each call transfers only the
 * standardised descriptor in and the feature-space adjoint out, both a few tens of
 * kilobytes. The compact Jacobian blocks -- the largest object by an order of magnitude --
 * stay on the host, which is where they are contracted. That choice is what decides
 * whether the GPU helps at all here: the arithmetic is small, the traffic is not.
 */
class GpuInference {
public:
    GpuInference(const std::string& backend, const Matrix& centers, const Vector& beta,
        const Vector& inv_scale, double sigma);
    ~GpuInference();

    GpuInference(const GpuInference&) = delete;
    GpuInference& operator=(const GpuInference&) = delete;

    bool valid() const { return m_ctx != nullptr; }

    /**
     * @brief Evaluate the kernel part of the correction and its feature-space adjoint.
     * @param X      natoms x F, standardised
     * @param start  per atom, first centre of its element block
     * @param count  per atom, size of that block (0 for an unknown element)
     * @param W      out, natoms x F adjoint in raw descriptor space
     */
    bool run(const Matrix& X, const std::vector<int>& start, const std::vector<int>& count,
        double& energy, Matrix& W, double& uncertainty);

private:
    void* m_ctx = nullptr;
    int m_F = 0;

    void* (*m_create)(int, int, const double*, const double*, const double*, double) = nullptr;
    int (*m_run)(void*, int, const double*, const int*, const int*, double*, double*, double*)
        = nullptr;
    void (*m_destroy)(void*) = nullptr;
};

} // namespace ml
} // namespace curcuma
