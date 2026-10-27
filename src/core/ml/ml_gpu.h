/*
 * < GPU offload of the Delta-ML model >
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
 * Claude Generated (July 2026, rewritten October 2026): GPU path of the Delta-ML work package.
 */

#pragma once

#include "descriptor.h"

#include "src/core/global.h"

#include <mutex>
#include <string>
#include <vector>

namespace curcuma {
namespace ml {

/**
 * @brief Device-resident Delta-ML model: descriptor, kernel, forces and the fit's normal equations.
 *
 * Thin RAII wrapper around the plugin entry points curcuma_<backend>_ml_model_*. The plugin is
 * dlopen'd on demand through the cache the method factory uses, so a CPU-only build or a machine
 * without a device reports unavailable (valid() false, reason in error()) and the caller keeps
 * its CPU path -- no CUDA symbol is linked into the main binary.
 *
 * Per prediction only coordinates go up and energy, indicator and gradient come down; the
 * element ids are re-sent only when they change. The ACSF descriptor is evaluated on the device
 * with the definitions of ACSFDescriptor::compute (see ml_model_cuda.cu), so results agree with
 * the CPU to rounding, not bit for bit.
 *
 * Calls are serialised by a mutex: one model may be shared by threads, one device context not.
 */
class GpuModel {
public:
    /**
     * @param centers   M x F standardised centres, grouped by element
     * @param mean      F feature means, @param scale F feature scales (never zero)
     */
    GpuModel(const std::string& backend, const ACSFConfig& acsf, const Matrix& centers,
        const Vector& mean, const Vector& scale, double sigma);
    ~GpuModel();

    GpuModel(const GpuModel&) = delete;
    GpuModel& operator=(const GpuModel&) = delete;

    /*! \brief True if the plugin loaded, the symbols resolved and the model is on the device. */
    bool valid() const { return m_ctx != nullptr; }
    /*! \brief Why the model is not valid, or why the last call failed. */
    const std::string& error() const { return m_error; }

    /*! \brief Upload the regression weights (length M); required before predict(). */
    bool setWeights(const Vector& beta);

    /**
     * @brief Kernel part of the correction, indicator and optionally the gradient.
     * @param sidx  per atom, ACSF species index (-1 unknown); @param start/@param count its centre block
     * @param gradient nullptr for energy only, else natoms x 3 [Eh/Angstrom] (resized)
     */
    bool predict(const std::vector<int>& sidx, const std::vector<int>& start,
        const std::vector<int>& count, const Geometry& geometry, double& e_kernel,
        double& uncertainty, Matrix* gradient);

    /*! \brief Zero the device normal equations (M x M). */
    bool beginAccumulation();
    /*! \brief Fold one structure in; @p y_grad (3N, row major) nullptr for an energy-only fit. */
    bool accumulate(const std::vector<int>& sidx, const std::vector<int>& start,
        const std::vector<int>& count, const Geometry& geometry, double y_energy,
        const double* y_grad, double wE, double wF);
    /*! \brief Retrieve the system, symmetrised into row-major @p A and @p b; frees it on the device. */
    bool fetch(Matrix& A, Vector& b);

private:
    bool fail(int code, const char* what);

    void* m_ctx = nullptr;
    int m_M = 0;
    std::string m_error;
    std::mutex m_mutex;

    void (*m_destroy)(void*) = nullptr;
    int (*m_set_beta)(void*, const double*) = nullptr;
    int (*m_predict)(void*, int, const double*, const int*, const int*, const int*, int, double*) = nullptr;
    int (*m_accum_begin)(void*) = nullptr;
    int (*m_accum_add)(void*, int, const double*, const int*, const int*, const int*, double,
        const double*, double, double) = nullptr;
    int (*m_accum_fetch)(void*, double*, double*) = nullptr;

    std::vector<double> m_out;
};

} // namespace ml
} // namespace curcuma
