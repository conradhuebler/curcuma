/*
 * < Atom-Centered Symmetry Functions (ACSF) for Delta-Machine-Learning >
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
 * Claude Generated (July 2026): descriptor layer of the Delta-ML work package.
 *
 * Literature:
 *   J. Behler, M. Parrinello, Phys. Rev. Lett. 98, 146401 (2007)
 *   J. Behler, J. Chem. Phys. 134, 074106 (2011)  -- the G2/G4 definitions used here
 */

#include "descriptor.h"

#include <algorithm>
#include <cmath>
#include <sstream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace curcuma {
namespace ml {

// =================================================================================
// ACSFConfig
// =================================================================================

ACSFConfig ACSFConfig::defaultFor(const std::vector<int>& species_list)
{
    ACSFConfig cfg;

    cfg.species = species_list;
    std::sort(cfg.species.begin(), cfg.species.end());
    cfg.species.erase(std::unique(cfg.species.begin(), cfg.species.end()), cfg.species.end());

    // Radial: eight Gaussian shells spanning the bonded region up to the cutoff tail.
    // A shifted grid (fixed width, varying centre) is better conditioned than the
    // alternative of a fixed centre with varying width, because the resulting features
    // are far less correlated.
    cfg.radial_rs.clear();
    const int n_shell = 8;
    const double r_first = 0.9;
    const double r_last = 5.2;
    for (int i = 0; i < n_shell; ++i)
        cfg.radial_rs.push_back(r_first + (r_last - r_first) * i / static_cast<double>(n_shell - 1));
    cfg.radial_eta = 8.0;

    // Angular: three sharpness levels times two orientations.
    cfg.angular_zeta = { 1.0, 4.0, 16.0 };
    cfg.angular_lambda = { 1, -1 };
    cfg.angular_eta = 0.05;

    cfg.rcut_radial = 6.0;
    cfg.rcut_angular = 4.0;

    return cfg;
}

int ACSFConfig::speciesIndex(int z) const
{
    const auto it = std::lower_bound(species.begin(), species.end(), z);
    if (it == species.end() || *it != z)
        return -1;
    return static_cast<int>(std::distance(species.begin(), it));
}

int ACSFConfig::pairIndex(int a, int b) const
{
    const int lo = std::min(a, b);
    const int hi = std::max(a, b);
    const int s = speciesCount();
    // Row-wise index of the upper triangle including the diagonal.
    return lo * s - (lo * (lo - 1)) / 2 + (hi - lo);
}

json ACSFConfig::toJson() const
{
    json j;
    j["rcut_radial"] = rcut_radial;
    j["rcut_angular"] = rcut_angular;
    j["radial_rs"] = radial_rs;
    j["radial_eta"] = radial_eta;
    j["angular_eta"] = angular_eta;
    j["angular_zeta"] = angular_zeta;
    j["angular_lambda"] = angular_lambda;
    j["species"] = species;
    return j;
}

ACSFConfig ACSFConfig::fromJson(const json& j)
{
    ACSFConfig cfg;
    cfg.rcut_radial = j.value("rcut_radial", 6.0);
    cfg.rcut_angular = j.value("rcut_angular", 4.0);
    cfg.radial_eta = j.value("radial_eta", 8.0);
    cfg.angular_eta = j.value("angular_eta", 0.05);
    if (j.contains("radial_rs"))
        cfg.radial_rs = j["radial_rs"].get<std::vector<double>>();
    if (j.contains("angular_zeta"))
        cfg.angular_zeta = j["angular_zeta"].get<std::vector<double>>();
    if (j.contains("angular_lambda"))
        cfg.angular_lambda = j["angular_lambda"].get<std::vector<int>>();
    if (j.contains("species"))
        cfg.species = j["species"].get<std::vector<int>>();
    std::sort(cfg.species.begin(), cfg.species.end());
    return cfg;
}

bool ACSFConfig::validate(std::string& error) const
{
    std::ostringstream msg;
    if (species.empty())
        msg << "no species defined; ";
    if (radial_rs.empty())
        msg << "no radial shifts defined; ";
    if (angular_zeta.empty() || angular_lambda.empty())
        msg << "no angular parameters defined; ";
    if (rcut_radial <= 0.0 || rcut_angular <= 0.0)
        msg << "cutoff must be positive; ";
    if (radial_eta <= 0.0 || angular_eta <= 0.0)
        msg << "eta must be positive; ";
    for (double rs : radial_rs) {
        if (rs < 0.0 || rs > rcut_radial) {
            msg << "radial shift outside [0, rcut_radial]; ";
            break;
        }
    }
    for (int lam : angular_lambda) {
        if (lam != 1 && lam != -1) {
            msg << "angular lambda must be +1 or -1; ";
            break;
        }
    }
    if (!std::is_sorted(species.begin(), species.end()))
        msg << "species list must be sorted; ";

    error = msg.str();
    return error.empty();
}

// =================================================================================
// ACSFDescriptor
// =================================================================================

ACSFDescriptor::ACSFDescriptor(const ACSFConfig& config)
    : m_config(config)
{
}

void ACSFDescriptor::compute(const std::vector<int>& elements, const Geometry& geometry,
    Matrix& X, std::vector<AtomJacobian>* jacobians) const
{
    const int natoms = static_cast<int>(elements.size());
    const int nfeat = m_config.featureCount();

    X = Matrix::Zero(natoms, nfeat);
    if (jacobians) {
        jacobians->clear();
        jacobians->resize(natoms);
    }
    if (natoms == 0 || nfeat == 0)
        return;

    const double rc_r = m_config.rcut_radial;
    const double rc_a = m_config.rcut_angular;
    const double rc_max = std::max(rc_r, rc_a);
    const int n_radial = m_config.radialCount();
    const int n_zeta = static_cast<int>(m_config.angular_zeta.size());
    const int n_lambda = static_cast<int>(m_config.angular_lambda.size());

    // Species index per atom; -1 marks an element the descriptor was not trained for.
    std::vector<int> sidx(natoms);
    for (int i = 0; i < natoms; ++i)
        sidx[i] = m_config.speciesIndex(elements[i]);

    // Reused scratch, allocated once instead of per centre atom.
    std::vector<int> nb;
    std::vector<Position> rij;
    std::vector<double> dij;
    nb.reserve(64);
    rij.reserve(64);
    dij.reserve(64);

    for (int i = 0; i < natoms; ++i) {
        AtomJacobian* J = jacobians ? &(*jacobians)[i] : nullptr;

        // An unknown element contributes no descriptor and no derivative.
        if (sidx[i] < 0) {
            if (J) {
                J->partners.assign(1, i);
                J->d = Matrix::Zero(3, nfeat);
            }
            continue;
        }

        // --- neighbour list of the centre atom -------------------------------------
        nb.clear();
        rij.clear();
        dij.clear();
        for (int j = 0; j < natoms; ++j) {
            if (j == i || sidx[j] < 0)
                continue;
            const Position v(geometry(j, 0) - geometry(i, 0), geometry(j, 1) - geometry(i, 1),
                geometry(j, 2) - geometry(i, 2));
            const double r = v.norm();
            if (r >= rc_max || r < 1e-8)
                continue;
            nb.push_back(j);
            rij.push_back(v);
            dij.push_back(r);
        }
        const int nnb = static_cast<int>(nb.size());

        if (J) {
            // Partner slot 0 is the centre atom itself, slots 1..nnb are its neighbours.
            J->partners.resize(nnb + 1);
            J->partners[0] = i;
            for (int t = 0; t < nnb; ++t)
                J->partners[t + 1] = nb[t];
            J->d = Matrix::Zero(3 * (nnb + 1), nfeat);
        }

        // --- radial part: Behler G2 -------------------------------------------------
        for (int t = 0; t < nnb; ++t) {
            const double r = dij[t];
            if (r >= rc_r)
                continue;

            const double fc = cutoff(r, rc_r);
            const double dfc = cutoffDerivative(r, rc_r);
            const int sj = sidx[nb[t]];
            const Position u = rij[t] / r; // dR_ij / dr_j

            for (int p = 0; p < n_radial; ++p) {
                const double dr = r - m_config.radial_rs[p];
                const double e = std::exp(-m_config.radial_eta * dr * dr);
                const int f = m_config.radialFeature(sj, p);

                X(i, f) += e * fc;
                if (!J)
                    continue;

                // d/dR_ij [ exp(-eta dr^2) fc(R) ]
                const double dval = e * (-2.0 * m_config.radial_eta * dr) * fc + e * dfc;
                for (int g = 0; g < 3; ++g) {
                    const double contrib = dval * u(g);
                    J->d(3 * (t + 1) + g, f) += contrib; // neighbour j
                    J->d(g, f) -= contrib;               // centre i
                }
            }
        }

        // --- angular part: Behler G4 ------------------------------------------------
        for (int t = 0; t < nnb; ++t) {
            const double a = dij[t];
            if (a >= rc_a)
                continue;

            for (int v = t + 1; v < nnb; ++v) {
                const double b = dij[v];
                if (b >= rc_a)
                    continue;

                const Position r_jk = rij[v] - rij[t]; // r_k - r_j
                const double c = r_jk.norm();
                // fc(R_jk) is part of G4, so a distant j-k pair contributes nothing.
                if (c >= rc_a || c < 1e-8)
                    continue;

                const double cos_theta = rij[t].dot(rij[v]) / (a * b);
                const int pi = m_config.pairIndex(sidx[nb[t]], sidx[nb[v]]);

                const double fa = cutoff(a, rc_a);
                const double fb = cutoff(b, rc_a);
                const double fcc = cutoff(c, rc_a);
                const double F = fa * fb * fcc;
                const double E = std::exp(-m_config.angular_eta * (a * a + b * b + c * c));

                double dfa = 0.0, dfb = 0.0, dfcc = 0.0;
                Position u_a, u_b, u_c, dcos_dj, dcos_dk;
                if (J) {
                    dfa = cutoffDerivative(a, rc_a);
                    dfb = cutoffDerivative(b, rc_a);
                    dfcc = cutoffDerivative(c, rc_a);
                    u_a = rij[t] / a;  // da / dr_j
                    u_b = rij[v] / b;  // db / dr_k
                    u_c = r_jk / c;    // dc / dr_k  (and -u_c for dc / dr_j)
                    dcos_dj = rij[v] / (a * b) - cos_theta * rij[t] / (a * a);
                    dcos_dk = rij[t] / (a * b) - cos_theta * rij[v] / (b * b);
                }

                int q = 0;
                for (int z = 0; z < n_zeta; ++z) {
                    const double zeta = m_config.angular_zeta[z];
                    const double pref = std::pow(2.0, 1.0 - zeta);

                    for (int l = 0; l < n_lambda; ++l) {
                        const double lambda = static_cast<double>(m_config.angular_lambda[l]);
                        // cos_theta is in [-1,1] and lambda is +-1, so base is in [0,2];
                        // the clamp only absorbs rounding at the endpoints.
                        double base = 1.0 + lambda * cos_theta;
                        if (base < 0.0)
                            base = 0.0;

                        const double P = std::pow(base, zeta);
                        const int f = m_config.angularFeature(pi, q);
                        X(i, f) += pref * P * E * F;

                        if (J) {
                            const double dP = zeta * lambda * std::pow(base, zeta - 1.0);
                            const double dt_dcos = pref * dP * E * F;
                            const double dt_da
                                = pref * P * ((-2.0 * m_config.angular_eta * a) * E * F + E * dfa * fb * fcc);
                            const double dt_db
                                = pref * P * ((-2.0 * m_config.angular_eta * b) * E * F + E * fa * dfb * fcc);
                            const double dt_dc
                                = pref * P * ((-2.0 * m_config.angular_eta * c) * E * F + E * fa * fb * dfcc);

                            for (int g = 0; g < 3; ++g) {
                                const double gj = dt_dcos * dcos_dj(g) + dt_da * u_a(g) - dt_dc * u_c(g);
                                const double gk = dt_dcos * dcos_dk(g) + dt_db * u_b(g) + dt_dc * u_c(g);
                                J->d(3 * (t + 1) + g, f) += gj;
                                J->d(3 * (v + 1) + g, f) += gk;
                                // Translation invariance fixes the centre contribution.
                                J->d(g, f) -= (gj + gk);
                            }
                        }
                        ++q;
                    }
                }
            }
        }
    }
}

} // namespace ml
} // namespace curcuma
