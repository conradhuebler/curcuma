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

#include "ml_gpu.h"

#include "src/core/curcuma_logger.h"
#include "src/core/energy_calculators/gpu_plugin.h"

#include <vector>

namespace curcuma {
namespace ml {

GpuAccumulator::GpuAccumulator(const std::string& backend, int M)
    : m_M(M)
{
    if (M <= 0)
        return;

    const std::string prefix = "curcuma_" + backend + "_ml_accum_";
    m_create = reinterpret_cast<void* (*)(int)>(
        gpu_plugin::resolveEntryPoint(backend, prefix + "create"));
    m_add = reinterpret_cast<int (*)(void*, int, const double*, const double*, double)>(
        gpu_plugin::resolveEntryPoint(backend, prefix + "add"));
    m_fetch = reinterpret_cast<int (*)(void*, double*, double*)>(
        gpu_plugin::resolveEntryPoint(backend, prefix + "fetch"));
    m_destroy = reinterpret_cast<void (*)(void*)>(
        gpu_plugin::resolveEntryPoint(backend, prefix + "destroy"));

    if (!m_create || !m_add || !m_fetch || !m_destroy) {
        CurcumaLogger::warn_fmt(
            "delta-ML: backend '{}' provides no accumulation entry points; fitting on the CPU",
            backend);
        return;
    }

    m_ctx = m_create(M);
    if (!m_ctx) {
        CurcumaLogger::warn_fmt(
            "delta-ML: backend '{}' could not allocate a {}x{} device system ({:.2f} GB); "
            "fitting on the CPU",
            backend, M, M, static_cast<double>(M) * M * 8.0 / 1e9);
    }
}

GpuAccumulator::~GpuAccumulator()
{
    if (m_ctx && m_destroy)
        m_destroy(m_ctx);
}

bool GpuAccumulator::add(int rows, const double* phi, const double* y, double alpha)
{
    if (!m_ctx)
        return false;
    return m_add(m_ctx, rows, phi, y, alpha) == 0;
}

bool GpuAccumulator::fetch(Matrix& A, Vector& b)
{
    if (!m_ctx)
        return false;

    std::vector<double> a_raw(static_cast<std::size_t>(m_M) * m_M, 0.0);
    std::vector<double> b_raw(static_cast<std::size_t>(m_M), 0.0);
    if (m_fetch(m_ctx, a_raw.data(), b_raw.data()) != 0)
        return false;

    // The device buffer holds the lower triangle in column-major order, i.e. entry (i,j)
    // with i >= j sits at j*M + i. Mirror it into the row-major result.
    A = Matrix::Zero(m_M, m_M);
    for (int j = 0; j < m_M; ++j) {
        for (int i = j; i < m_M; ++i) {
            const double v = a_raw[static_cast<std::size_t>(j) * m_M + i];
            A(i, j) = v;
            A(j, i) = v;
        }
    }
    b = Eigen::Map<const Vector>(b_raw.data(), m_M);
    return true;
}

} // namespace ml
} // namespace curcuma
