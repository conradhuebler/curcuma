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

#include "ml_gpu.h"
#include "ml_gpu_abi.h"

#include "src/core/energy_calculators/gpu_plugin.h"

#include <cmath>

namespace curcuma {
namespace ml {

namespace {
    const char* codeText(int code)
    {
        switch (code) {
        case CURCUMA_ML_ERR_ARGS:
            return "invalid arguments";
        case CURCUMA_ML_ERR_ALLOC:
            return "device memory allocation failed";
        case CURCUMA_ML_ERR_CUDA:
            return "device runtime error";
        case CURCUMA_ML_ERR_BLAS:
            return "device BLAS error";
        case CURCUMA_ML_ERR_NOT_READY:
            return "device model not initialised";
        case CURCUMA_ML_ERR_NEIGHBOURS:
            return "an atom has more neighbours inside the descriptor cutoff than the device "
                   "kernels hold (CURCUMA_ML_MAX_NEIGHBOURS = 512)";
        default:
            return "unknown device error";
        }
    }
} // namespace

GpuModel::GpuModel(const std::string& backend_in, const ACSFConfig& acsf, const Matrix& centers,
    const Vector& mean, const Vector& scale, double sigma)
    : m_M(static_cast<int>(centers.rows()))
{
    const std::string backend = (backend_in == "auto") ? std::string("cuda") : backend_in;
    const int F = static_cast<int>(centers.cols());
    if (m_M <= 0 || F <= 0 || mean.size() != F || scale.size() != F) {
        m_error = "model has no centres or inconsistent feature scaling";
        return;
    }

    // The device kernels hold the grid in fixed arrays; refuse a grid that does not fit
    // rather than run a truncated one.
    CurcumaMLDeviceConfig cfg {};
    cfg.n_species = acsf.speciesCount();
    cfg.n_radial = acsf.radialCount();
    cfg.n_zeta = static_cast<int>(acsf.angular_zeta.size());
    cfg.n_lambda = static_cast<int>(acsf.angular_lambda.size());
    if (cfg.n_radial > CURCUMA_ML_MAX_RADIAL || cfg.n_zeta > CURCUMA_ML_MAX_ZETA
        || cfg.n_lambda > CURCUMA_ML_MAX_LAMBDA || acsf.featureCount() != F) {
        m_error = "descriptor grid exceeds the device limits (radial <= 32, zeta <= 8, lambda <= 8)";
        return;
    }
    cfg.rcut_radial = acsf.rcut_radial;
    cfg.rcut_angular = acsf.rcut_angular;
    cfg.radial_eta = acsf.radial_eta;
    cfg.angular_eta = acsf.angular_eta;
    for (int p = 0; p < cfg.n_radial; ++p)
        cfg.radial_rs[p] = acsf.radial_rs[p];
    for (int z = 0; z < cfg.n_zeta; ++z) {
        cfg.zeta[z] = acsf.angular_zeta[z];
        cfg.zeta_pref[z] = std::pow(2.0, 1.0 - acsf.angular_zeta[z]); // as in descriptor.cpp
        const double zr = std::round(acsf.angular_zeta[z]);
        cfg.zeta_int[z] = (zr == acsf.angular_zeta[z] && zr >= 1.0 && zr <= 64.0) ? static_cast<int>(zr) : 0;
    }
    for (int l = 0; l < cfg.n_lambda; ++l)
        cfg.lambda[l] = static_cast<double>(acsf.angular_lambda[l]);

    const std::string prefix = "curcuma_" + backend + "_ml_model_";
    auto create = reinterpret_cast<void* (*)(const CurcumaMLDeviceConfig*, int, int, const double*,
        const double*, const double*, double)>(gpu_plugin::resolveEntryPoint(backend, prefix + "create"));
    m_destroy = reinterpret_cast<void (*)(void*)>(gpu_plugin::resolveEntryPoint(backend, prefix + "destroy"));
    m_set_beta = reinterpret_cast<int (*)(void*, const double*)>(
        gpu_plugin::resolveEntryPoint(backend, prefix + "set_beta"));
    m_predict = reinterpret_cast<int (*)(void*, int, const double*, const int*, const int*, const int*,
        int, double*)>(gpu_plugin::resolveEntryPoint(backend, prefix + "predict"));
    m_accum_begin = reinterpret_cast<int (*)(void*)>(
        gpu_plugin::resolveEntryPoint(backend, prefix + "accum_begin"));
    m_accum_add = reinterpret_cast<int (*)(void*, int, const double*, const int*, const int*, const int*,
        double, const double*, double, double)>(gpu_plugin::resolveEntryPoint(backend, prefix + "accum_add"));
    m_accum_fetch = reinterpret_cast<int (*)(void*, double*, double*)>(
        gpu_plugin::resolveEntryPoint(backend, prefix + "accum_fetch"));

    if (!create || !m_destroy || !m_set_beta || !m_predict || !m_accum_begin || !m_accum_add || !m_accum_fetch) {
        m_error = "backend '" + backend + "' provides no delta-ML entry points (plugin not built or not found)";
        return;
    }

    const Vector inv_scale = scale.cwiseInverse(); // same expression as DeltaMLModel::standardise
    m_ctx = create(&cfg, m_M, F, centers.data(), mean.data(), inv_scale.data(), sigma);
    if (!m_ctx)
        m_error = "backend '" + backend + "' found no usable device or could not upload the model";
}

GpuModel::~GpuModel()
{
    if (m_ctx && m_destroy)
        m_destroy(m_ctx);
}

bool GpuModel::fail(int code, const char* what)
{
    m_error = std::string(what) + ": " + codeText(code);
    return false;
}

bool GpuModel::setWeights(const Vector& beta)
{
    std::lock_guard<std::mutex> lock(m_mutex);
    if (!m_ctx || beta.size() != m_M)
        return false;
    const int rc = m_set_beta(m_ctx, beta.data());
    return rc == CURCUMA_ML_OK ? true : fail(rc, "uploading the weights");
}

bool GpuModel::predict(const std::vector<int>& sidx, const std::vector<int>& start,
    const std::vector<int>& count, const Geometry& geometry, double& e_kernel, double& uncertainty,
    Matrix* gradient)
{
    std::lock_guard<std::mutex> lock(m_mutex);
    const int n = static_cast<int>(sidx.size());
    if (!m_ctx || n == 0 || geometry.rows() != n || geometry.cols() != 3)
        return false;
    m_out.assign(2 + 3 * static_cast<std::size_t>(n), 0.0);
    // Geometry is a row-major N x 3 typedef, which is the layout the device reads.
    const int rc = m_predict(m_ctx, n, geometry.data(), sidx.data(), start.data(), count.data(),
        gradient ? 1 : 0, m_out.data());
    if (rc != CURCUMA_ML_OK)
        return fail(rc, "prediction");
    e_kernel = m_out[0];
    uncertainty = m_out[1];
    if (gradient)
        *gradient = Eigen::Map<const Matrix>(m_out.data() + 2, n, 3);
    return true;
}

bool GpuModel::beginAccumulation()
{
    std::lock_guard<std::mutex> lock(m_mutex);
    if (!m_ctx)
        return false;
    const int rc = m_accum_begin(m_ctx);
    return rc == CURCUMA_ML_OK ? true : fail(rc, "allocating the normal equations");
}

bool GpuModel::accumulate(const std::vector<int>& sidx, const std::vector<int>& start,
    const std::vector<int>& count, const Geometry& geometry, double y_energy, const double* y_grad,
    double wE, double wF)
{
    std::lock_guard<std::mutex> lock(m_mutex);
    const int n = static_cast<int>(sidx.size());
    if (!m_ctx || n == 0 || geometry.rows() != n || geometry.cols() != 3)
        return false;
    const int rc = m_accum_add(m_ctx, n, geometry.data(), sidx.data(), start.data(), count.data(),
        y_energy, y_grad, wE, wF);
    return rc == CURCUMA_ML_OK ? true : fail(rc, "building the design block");
}

bool GpuModel::fetch(Matrix& A, Vector& b)
{
    std::lock_guard<std::mutex> lock(m_mutex);
    if (!m_ctx)
        return false;
    std::vector<double> a_raw(static_cast<std::size_t>(m_M) * m_M, 0.0);
    std::vector<double> b_raw(static_cast<std::size_t>(m_M), 0.0);
    const int rc = m_accum_fetch(m_ctx, a_raw.data(), b_raw.data());
    if (rc != CURCUMA_ML_OK)
        return fail(rc, "reading the normal equations back");

    // Lower triangle in column-major order: entry (i,j), i >= j, sits at j*M + i.
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
