/*
 * <Path collective variables (Branduardi s,z) - implementation.>
 * Copyright (C) 2019 - 2026 Conrad Huebler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated 2026 (AI-generated, machine-tested pending - not TESTED/APPROVED).
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
 */

#include "src/capabilities/path_cv.h"

#include <cmath>
#include <limits>

#include "src/core/curcuma_logger.h"
#include "src/core/parameter_registry.h"

PathCV::PathCV(const std::vector<Geometry>& path, const Molecule& templ, double lambda)
    : m_path(path)
    , m_template(templ)
{
    if (m_path.size() < 2)
        CurcumaLogger::error("PathCV: need at least two reference frames.");
    m_lambda = (lambda > 0.0) ? lambda : computeLambdaFromPath(m_path);
}

double PathCV::computeLambdaFromPath(const std::vector<Geometry>& path)
{
    // lambda ~ 1/<d^2> with d the RMSD spacing of neighbouring frames, so that the
    // Gaussians of adjacent frames overlap appreciably. Using the plain (non-best-fit)
    // spacing here is sufficient for the scale estimate.
    if (path.size() < 2)
        return 1.0;
    const int nat = static_cast<int>(path[0].rows());
    double sum = 0.0;
    int cnt = 0;
    for (size_t i = 1; i < path.size(); ++i) {
        const double d2 = (path[i] - path[i - 1]).squaredNorm() / std::max(nat, 1);
        sum += d2;
        ++cnt;
    }
    const double mean_d2 = (cnt > 0) ? sum / cnt : 1.0;
    return (mean_d2 > 1e-12) ? 1.0 / mean_d2 : 1.0;
}

double PathCV::rmsdTo(int i, const Geometry& R, Geometry* grad) const
{
    // A fresh driver per call keeps this const and thread-safe; the RMSD of a few
    // dozen atoms is cheap compared with the energy evaluation it accompanies.
    // Start from the registry defaults - handing ConfigManager a partial rmsd block
    // makes it throw on the first parameter that is missing.
    json rmsd_mod = ParameterRegistry::getInstance().getDefaultJson("rmsd");
    rmsd_mod["method"] = "inertia";
    rmsd_mod["no_reorder"] = true;
    rmsd_mod["threads"] = 1;
    json cfg = json::object();
    cfg["rmsd"] = rmsd_mod;
    RMSDDriver driver(cfg, true);

    Molecule ref = m_template;
    ref.setGeometry(m_path[i]);
    Molecule tar = m_template;
    tar.setGeometry(R);
    driver.setReference(ref);
    driver.setTarget(tar);

    const double d = driver.BestFitRMSD();
    if (grad)
        *grad = driver.Gradient();
    return d;
}

void PathCV::evaluate(const Geometry& R, double* s_out, double* z_out, Geometry* ds_dR,
                      Geometry* dz_dR) const
{
    const int N = static_cast<int>(m_path.size());
    const int nat = static_cast<int>(R.rows());
    if (N < 2)
        return;

    std::vector<double> d(N, 0.0);
    std::vector<Geometry> gd(N);
    const bool need_grad = (ds_dR != nullptr) || (dz_dR != nullptr);
    for (int i = 0; i < N; ++i)
        d[i] = rmsdTo(i, R, need_grad ? &gd[i] : nullptr);

    // Softmax weights with the max subtracted for numerical stability.
    double emax = -std::numeric_limits<double>::max();
    for (int i = 0; i < N; ++i)
        emax = std::max(emax, -m_lambda * d[i] * d[i]);

    std::vector<double> w(N, 0.0);
    double wsum = 0.0;
    for (int i = 0; i < N; ++i) {
        w[i] = std::exp(-m_lambda * d[i] * d[i] - emax);
        wsum += w[i];
    }
    if (wsum < 1e-300)
        wsum = 1e-300;

    double s = 0.0;
    for (int i = 0; i < N; ++i)
        s += static_cast<double>(i) * w[i];
    s /= (wsum * (N - 1));

    if (s_out)
        *s_out = s;
    if (z_out)
        *z_out = -(1.0 / m_lambda) * (std::log(wsum) + emax);

    // dz/dR = sum_i (w_i/wsum) * 2 * d_i * dd_i/dR. The driver gradient is
    // -d(RMSD)/dR (see the sign note below), hence the leading minus here.
    if (dz_dR) {
        *dz_dR = Geometry::Zero(nat, 3);
        for (int i = 0; i < N; ++i) {
            if (gd[i].rows() != nat || gd[i].cols() != 3)
                continue;
            *dz_dR += (-(w[i] / wsum) * 2.0 * d[i]) * gd[i];
        }
    }

    if (!ds_dR)
        return;

    // ds/dR = sum_i [ i/(N-1) - s ] * (w_i/wsum) * (-2 lambda d_i) * dd_i/dR
    //
    // NOTE on the sign: RMSDDriver::Gradient() returns (ref - tar)/(RMSD*N), i.e.
    // MINUS d(RMSD)/dR_target (RMSD itself is sqrt(sum|delta|^2/N), see
    // RMSDFunctions::getRMSD). Folding that minus into the prefactor turns the
    // -2*lambda*d of the chain rule into +2*lambda*d. Getting this wrong shows up as
    // an exact factor -2 against finite differences, which is what the unit test
    // caught.
    *ds_dR = Geometry::Zero(nat, 3);
    for (int i = 0; i < N; ++i) {
        const double pref = (static_cast<double>(i) / (N - 1) - s) * (w[i] / wsum)
            * (2.0 * m_lambda * d[i]);
        if (gd[i].rows() == nat && gd[i].cols() == 3)
            *ds_dR += pref * gd[i];
    }
}

double PathCV::s(const Geometry& R) const
{
    double s_val = 0.0;
    evaluate(R, &s_val, nullptr, nullptr);
    return s_val;
}

double PathCV::z(const Geometry& R) const
{
    double z_val = 0.0;
    evaluate(R, nullptr, &z_val, nullptr);
    return z_val;
}
