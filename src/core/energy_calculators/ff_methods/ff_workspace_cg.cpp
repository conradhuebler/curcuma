/*
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
 */

// Claude Generated (Sep 2026): coarse-grained pair kernel for FFWorkspace.
//
// The CG method (`-method cg`, spheres/ellipsoids interacting through the LJ-type
// potentials in cg_potentials.h) used to live in the legacy ForceField thread engine,
// energy-only with numerical forces. It is now a workspace method like UFF/QMDFF:
// every `vdW` entry of type 3 is one CG pair; the energy is CGPotentials::
// calculateCGPairEnergy (unchanged physics, distances in the workspace geometry
// units = Bohr, as before). Gradient: analytic for spheres (the effective sigma is
// then a constant), central finite differences of the same energy function for
// ellipsoids (their orientation is fixed, so the pair energy depends on r_i - r_j
// only and the two gradients are antisymmetric).

#include "ff_workspace.h"
#include "cg_potentials.h"

#include <cmath>

void FFWorkspace::executeCG(int p)
{
    calcCGPairs(p);
}

void FFWorkspace::calcCGPairs(int p)
{
    auto& acc = m_accumulators[p];
    const auto [beg, end] = m_partitions[p].vdws;
    for (int idx = beg; idx < end; ++idx) {
        const vdW& pair = m_vdws[idx];
        if (pair.type != 3) continue;
        const Eigen::Vector3d ri = m_geometry.row(pair.i);
        const Eigen::Vector3d rj = m_geometry.row(pair.j);
        acc.energy.vdw += CGPotentials::calculateCGPairEnergy(pair, ri, rj);
        if (!m_do_gradient) continue;

        const Eigen::Vector3d d = ri - rj;
        const double r = d.norm();
        Eigen::Vector3d g = Eigen::Vector3d::Zero();
        const bool spheres = CGPotentials::isSpherical(pair.shape_i) && CGPotentials::isSpherical(pair.shape_j);
        if (spheres && r >= 0.1) {
            // E(r) with sigma_eff = sigma:  type 1: eps*(s^12 - s^6 + 1),  type 2: 4 eps (s^12 - s^6),  s = sigma/r
            const double s = pair.sigma / r;
            const double s6 = std::pow(s, 6);
            const double s12 = s6 * s6;
            double dEdr = 0.0;
            if (pair.cg_potential_type == 1)      dEdr = pair.epsilon * (-12.0 * s12 + 6.0 * s6) / r;
            else if (pair.cg_potential_type == 2) dEdr = 4.0 * pair.epsilon * (-12.0 * s12 + 6.0 * s6) / r;
            g = dEdr * d / r;
        } else if (r >= 0.1) {
            // Ellipsoids (or near-contact): central differences of the pair energy.
            const double h = 1e-5;
            for (int k = 0; k < 3; ++k) {
                Eigen::Vector3d rp = ri, rm = ri;
                rp[k] += h; rm[k] -= h;
                g[k] = (CGPotentials::calculateCGPairEnergy(pair, rp, rj)
                      - CGPotentials::calculateCGPairEnergy(pair, rm, rj)) / (2.0 * h);
            }
        }
        acc.gradient.row(pair.i) += g.transpose();
        acc.gradient.row(pair.j) -= g.transpose();
    }
}
