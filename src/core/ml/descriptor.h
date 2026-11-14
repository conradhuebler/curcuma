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
 */

#pragma once

#include "src/core/global.h"
#include "json.hpp"

#include <string>
#include <vector>

using json = nlohmann::json;

namespace curcuma {
namespace ml {

/**
 * @brief Compact Jacobian of one atom's descriptor with respect to atomic positions.
 *
 * The descriptor x_i of atom i depends only on atoms inside the cutoff sphere. Storing
 * a dense (3*N_atoms x n_features) block per atom would scale quadratically in memory
 * (288 MB for 500 atoms); this compact form keeps only the atoms that actually enter,
 * which is linear (~17 MB for 500 atoms).
 *
 * Row layout: derivative with respect to component gamma of partner p lives in row
 * 3*p + gamma, i.e. d(3*p + gamma, f) = d x_{i,f} / d R_{partners[p], gamma}.
 *
 * Units: descriptor values are dimensionless, positions are in Angstrom, so the entries
 * carry 1/Angstrom.
 */
struct AtomJacobian {
    std::vector<int> partners; ///< atom indices entering x_i, partners[0] is always i itself
    Matrix d;                  ///< (3*partners.size()) x n_features
};

/**
 * @brief User-facing knobs from which ACSFConfig::defaultFor builds a descriptor.
 *
 * Defaults (October 2026): six radial shells 0.9..5.2 Angstrom, eta_r 4, cutoffs 6/4 Angstrom,
 * eta_a 0.05, zeta {1,4,16}, lambda +-1. Until October 2026 the default grid had eight shells
 * with eta_r 8; on a 90-atom conformer set the 6/4 grid gave a held-out relative energy MAE
 * 0.215 +- 0.100 kcal/mol lower (10 split seeds, see docs/DELTA_ML.md). Models written with
 * the old grid keep it: the model file stores the full descriptor block, not these settings.
 * Radial shifts form an equidistant grid R_p = r_first + (r_last - r_first) p/(n_shell - 1).
 *
 * Claude Generated (October 2026): exposed for hyperparameter studies.
 */
struct ACSFSettings {
    int n_shell = 6;                ///< number of radial shells (was 8 before Oct 2026)
    double r_first = 0.9;           ///< first radial shift [Angstrom]
    double r_last = 5.2;            ///< last radial shift [Angstrom]
    double radial_eta = 4.0;        ///< radial width [1/Angstrom^2] (was 8.0 before Oct 2026)
    double rcut_radial = 6.0;       ///< radial cutoff [Angstrom]
    double rcut_angular = 4.0;      ///< angular cutoff [Angstrom]
    double angular_eta = 0.05;      ///< angular damping [1/Angstrom^2]
    std::vector<double> angular_zeta = { 1.0, 4.0, 16.0 }; ///< angular sharpness exponents

    /*! \brief Parse a comma list such as "1,4,16"; returns false on a malformed or non-positive entry. */
    static bool parseZetaList(const std::string& text, std::vector<double>& out);
};

/**
 * @brief Parameters of the ACSF descriptor (Behler & Parrinello).
 *
 * Radial part (Behler G2), for centre i and neighbour element s:
 *   G2_{i,s,p} = sum_{j: Z_j=s} exp(-eta_r (R_ij - R_p)^2) * fc(R_ij, rcut_radial)
 *
 * Angular part (Behler G4), for centre i and unordered element pair (s,t):
 *   G4_{i,st,q} = 2^(1-zeta) * sum_{j<k} (1 + lambda*cos(theta_ijk))^zeta
 *                 * exp(-eta_a (R_ij^2 + R_ik^2 + R_jk^2))
 *                 * fc(R_ij) * fc(R_ik) * fc(R_jk)
 *
 * The cutoff function fc(R) = 0.5*(cos(pi*R/Rc) + 1) vanishes together with its first
 * derivative at Rc. That C1 property is not cosmetic: a descriptor with a kink produces
 * a force discontinuity, which shows up as energy drift in NVE molecular dynamics.
 *
 * Feature layout (see index helpers below):
 *   [ radial: n_species * n_radial ] [ angular: n_pairs * n_angular ]
 */
struct ACSFConfig {
    double rcut_radial = 6.0;  ///< radial cutoff [Angstrom]
    double rcut_angular = 4.0; ///< angular cutoff [Angstrom], limits the O(n_nb^2) cost

    std::vector<double> radial_rs;   ///< radial shifts R_p [Angstrom]
    double radial_eta = 8.0;         ///< radial width parameter [1/Angstrom^2]

    double angular_eta = 0.05;         ///< angular radial damping [1/Angstrom^2]
    std::vector<double> angular_zeta;  ///< angular sharpness exponents
    std::vector<int> angular_lambda;   ///< +1 / -1, selects maxima at theta = 0 / pi

    std::vector<int> species; ///< sorted, unique atomic numbers the descriptor knows about

    /*! \brief Sensible starting point for organic molecules; `settings` overrides the grid
     *  (default-constructed settings give the current default grid, see ACSFSettings). */
    static ACSFConfig defaultFor(
        const std::vector<int>& species_list, const ACSFSettings& settings = ACSFSettings{});

    int radialCount() const { return static_cast<int>(radial_rs.size()); }
    int angularCount() const
    {
        return static_cast<int>(angular_zeta.size() * angular_lambda.size());
    }
    int speciesCount() const { return static_cast<int>(species.size()); }
    int pairCount() const
    {
        const int s = speciesCount();
        return s * (s + 1) / 2;
    }
    int featureCount() const
    {
        return speciesCount() * radialCount() + pairCount() * angularCount();
    }

    /*! \brief Position of element Z in the species list, or -1 if unknown. */
    int speciesIndex(int z) const;

    /*! \brief Triangular index of the unordered species pair (a,b), both species indices. */
    int pairIndex(int a, int b) const;

    /*! \brief Feature index of radial function p around a neighbour of species index si. */
    int radialFeature(int si, int p) const { return si * radialCount() + p; }

    /*! \brief Feature index of angular function q for species pair index pi. */
    int angularFeature(int pi, int q) const
    {
        return speciesCount() * radialCount() + pi * angularCount() + q;
    }

    json toJson() const;
    static ACSFConfig fromJson(const json& j);
    bool validate(std::string& error) const;
};

/**
 * @brief Evaluates ACSF descriptors and their analytic position derivatives.
 *
 * The value and derivative code share one traversal on purpose: duplicating the geometric
 * loops is the classic way to end up with an energy and a gradient that disagree.
 */
class ACSFDescriptor {
public:
    explicit ACSFDescriptor(const ACSFConfig& config);

    const ACSFConfig& config() const { return m_config; }
    int featureCount() const { return m_config.featureCount(); }

    /**
     * @brief Compute descriptors, optionally with derivatives.
     * @param elements  atomic numbers, length N
     * @param geometry  N x 3 coordinates [Angstrom]
     * @param X         output, N x n_features (resized by the call)
     * @param jacobians optional output, one compact block per atom; nullptr to skip
     *
     * Descriptors of elements outside ACSFConfig::species are left at zero, and the
     * corresponding atoms contribute nothing. That is deliberate: an unseen element must
     * degrade to "no correction" instead of silently reusing another element's parameters.
     */
    void compute(const std::vector<int>& elements, const Geometry& geometry, Matrix& X,
        std::vector<AtomJacobian>* jacobians = nullptr) const;

private:
    /*! \brief Cutoff function and its derivative. */
    static inline double cutoff(double r, double rc)
    {
        if (r >= rc)
            return 0.0;
        return 0.5 * (std::cos(M_PI * r / rc) + 1.0);
    }
    static inline double cutoffDerivative(double r, double rc)
    {
        if (r >= rc)
            return 0.0;
        return -0.5 * M_PI / rc * std::sin(M_PI * r / rc);
    }

    ACSFConfig m_config;
};

} // namespace ml
} // namespace curcuma
