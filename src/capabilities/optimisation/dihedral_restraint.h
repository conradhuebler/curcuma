/*
 * <Harmonic dihedral restraints for geometry optimisation>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
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
 * Claude Generated (Jul 2026)
 */

#pragma once

#include <cmath>
#include <vector>

#include "json.hpp"

#include "src/core/energy_calculators/ff_methods/gfnff_geometry.h"
#include "src/tools/geometry.h"

namespace Optimization {

/**
 * @brief One harmonic restraint on a dihedral angle i-j-k-l.
 *
 * Physics (the whole method, in one line):
 *
 *     E_restraint = 1/2 * k * (phi - phi_target)^2
 *
 * with the deviation wrapped into (-pi, pi] so that a restraint to +170 deg pulls a torsion sitting
 * at -175 deg the short way round (15 deg), not the long way (345 deg).
 *
 * The gradient follows the chain rule, dE/dr = k * dphi * dphi/dr, with dphi/dr from
 * GFNFF_Geometry::calculateDihedralAngle (the same analytic derivative the force field's torsion
 * term uses).
 *
 * Why it exists: recombining conformers by setting torsions on a rigid template creates clashes --
 * measured on a compact 90-atom molecule, 72 % of the proposals were unusable for that reason. With
 * a restraint the torsion can be DRIVEN to its target while the rest of the molecule relaxes out of
 * the way, and released afterwards. See docs/CONFSEARCH_PROPOSALS.md (P0).
 */
struct DihedralRestraint {
    int i = -1, j = -1, k = -1, l = -1; ///< atom indices (0-based) defining the dihedral
    double target = 0.0;                ///< target angle in RADIANS
    double force = 0.1;                 ///< force constant in Eh/rad^2
};

/**
 * @brief A set of dihedral restraints, applied to an energy/gradient pair.
 *
 * UNITS -- the single most error-prone part, so it is stated explicitly:
 *  - the geometry passed in is in ANGSTROM (curcuma's Molecule geometry),
 *  - the gradient is in Eh/BOHR (curcuma's optimiser convention),
 *  - dphi/dr therefore comes out in rad/Angstrom and is converted with
 *    d/dBohr = d/dAngstrom * (Angstrom per Bohr) = * 0.529177...
 *  - the returned energy is in Hartree and MUST be added to the method energy: the LBFGSpp line
 *    search only accepts a step that lowers the RETURNED energy, so a gradient-only bias would be
 *    fought by the line search instead of followed (the same lesson as the interactive-grab bias in
 *    LBFGSppObjectiveFunction).
 */
class DihedralRestraints {
public:
    DihedralRestraints() = default;

    /** Parse from an optimiser config: config["dihedral_restraints"] = [ {i,j,k,l,target,force}, ... ]
     *  "target" is read in DEGREES (that is what a user or a torsion-space state centre speaks),
     *  "force" in Eh/rad^2. Malformed entries are skipped. */
    static DihedralRestraints fromJson(const nlohmann::json& config)
    {
        DihedralRestraints out;
        if (!config.contains("dihedral_restraints") || !config["dihedral_restraints"].is_array())
            return out;
        for (const auto& entry : config["dihedral_restraints"]) {
            if (!entry.is_object())
                continue;
            DihedralRestraint r;
            r.i = entry.value("i", -1);
            r.j = entry.value("j", -1);
            r.k = entry.value("k", -1);
            r.l = entry.value("l", -1);
            r.target = entry.value("target", 0.0) * M_PI / 180.0;
            r.force = entry.value("force", 0.1);
            if (r.i < 0 || r.j < 0 || r.k < 0 || r.l < 0)
                continue;
            out.m_restraints.push_back(r);
        }
        return out;
    }

    /** Serialise back into the json form fromJson() reads (target in degrees). */
    static nlohmann::json toJson(const std::vector<DihedralRestraint>& restraints)
    {
        nlohmann::json arr = nlohmann::json::array();
        for (const auto& r : restraints) {
            nlohmann::json e;
            e["i"] = r.i;
            e["j"] = r.j;
            e["k"] = r.k;
            e["l"] = r.l;
            e["target"] = r.target * 180.0 / M_PI;
            e["force"] = r.force;
            arr.push_back(e);
        }
        return arr;
    }

    bool empty() const { return m_restraints.empty(); }
    std::size_t size() const { return m_restraints.size(); }
    const std::vector<DihedralRestraint>& restraints() const { return m_restraints; }
    void add(const DihedralRestraint& r) { m_restraints.push_back(r); }

    /**
     * @brief Add the restraint energy and gradient.
     * @param geometry  atom positions in Angstrom (natoms x 3)
     * @param gradient  flat gradient in Eh/Bohr, atom-major (3*natoms); restraint terms are ADDED
     * @return the restraint energy in Hartree (to be added to the method energy)
     */
    double apply(const Geometry& geometry, Vector& gradient) const
    {
        if (m_restraints.empty())
            return 0.0;
        const int natoms = static_cast<int>(geometry.rows());
        double energy = 0.0;
        Matrix dphi(4, 3);
        for (const DihedralRestraint& r : m_restraints) {
            if (r.i >= natoms || r.j >= natoms || r.k >= natoms || r.l >= natoms)
                continue;
            const Eigen::Vector3d ri = geometry.row(r.i);
            const Eigen::Vector3d rj = geometry.row(r.j);
            const Eigen::Vector3d rk = geometry.row(r.k);
            const Eigen::Vector3d rl = geometry.row(r.l);
            const double phi = GFNFF_Geometry::calculateDihedralAngle(ri, rj, rk, rl, dphi, true);

            // Shortest angular deviation: wrap into (-pi, pi].
            double delta = phi - r.target;
            while (delta > M_PI)
                delta -= 2.0 * M_PI;
            while (delta <= -M_PI)
                delta += 2.0 * M_PI;

            energy += 0.5 * r.force * delta * delta;

            // dE/dr [Eh/Angstrom] = k * delta * dphi/dr [rad/Angstrom]; -> Eh/Bohr
            const double prefactor = r.force * delta * kAngstromPerBohr;
            const int idx[4] = { r.i, r.j, r.k, r.l };
            for (int a = 0; a < 4; ++a)
                for (int c = 0; c < 3; ++c)
                    gradient[3 * idx[a] + c] += prefactor * dphi(a, c);
        }
        return energy;
    }

    /** Current deviation of the worst restraint, in degrees (diagnostics / convergence reporting). */
    double maxDeviationDegrees(const Geometry& geometry) const
    {
        double worst = 0.0;
        Matrix dphi(4, 3);
        const int natoms = static_cast<int>(geometry.rows());
        for (const DihedralRestraint& r : m_restraints) {
            if (r.i >= natoms || r.j >= natoms || r.k >= natoms || r.l >= natoms)
                continue;
            const double phi = GFNFF_Geometry::calculateDihedralAngle(
                Eigen::Vector3d(geometry.row(r.i)), Eigen::Vector3d(geometry.row(r.j)),
                Eigen::Vector3d(geometry.row(r.k)), Eigen::Vector3d(geometry.row(r.l)), dphi, false);
            double delta = phi - r.target;
            while (delta > M_PI)
                delta -= 2.0 * M_PI;
            while (delta <= -M_PI)
                delta += 2.0 * M_PI;
            worst = std::max(worst, std::abs(delta) * 180.0 / M_PI);
        }
        return worst;
    }

private:
    // CODATA-2018 Bohr radius in Angstrom (mirrors CurcumaUnit; kept local so this header stays
    // usable from the optimiser without pulling in the unit system).
    static constexpr double kAngstromPerBohr = 0.52917721067;

    std::vector<DihedralRestraint> m_restraints;
};

} // namespace Optimization
