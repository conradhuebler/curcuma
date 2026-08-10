/*
 * <GFN-FF bonded terms at unit force constant, for Hessian parametrisation>
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
 * =============================================================================
 * WHY THIS FILE EXISTS
 * =============================================================================
 *
 * Every GFN-FF bonded term is LINEAR in exactly one constant:
 *
 *   bond      E = fc * exp(-alpha * dr^2)                     -> Bond::fc
 *   angle     E = fc * (cos(theta) - cos(theta0))^2 * damp    -> Angle::fc
 *   torsion   E = V  * (1 + cos(n(phi-phi0) + pi)) * damp     -> Dihedral::V
 *   inversion E = V  * (1 - cos(omega)) * damp                -> Inversion::fc
 *
 * which is what allows GFN-FF's force constants to be determined by the same linear
 * Hessian fit that QMDFF uses (qmdff_parametrisation.h). The fit needs each term's
 * Hessian evaluated with its own constant set to 1, i.e. the term in isolation — but
 * the production kernels live inside FFWorkspace::calc{Bonds,Angles,Dihedrals,
 * Inversions}, which always evaluate whole partitioned lists.
 *
 * These are therefore standalone, single-term copies of those kernels. They are
 * transcriptions, not reimplementations: energy and gradient must agree with the
 * workspace to machine precision or the fit is meaningless.
 *
 * THAT AGREEMENT IS CHECKED, not assumed: `ctest -R gfnff_unit_terms` compares the
 * assembled sum_p k_p U_p against a Hessian produced by the real FFWorkspace, and every
 * `-qmdfffit` run repeats the same comparison on its own result. A divergence here shows
 * up immediately as a large verification residual rather than as silently wrong
 * force constants.
 *
 * Reference: gfnff_engrad.F90 via ff_workspace_gfnff.cpp.
 *
 * UNITS: atomic units throughout (Bohr, Hartree), matching FFWorkspace's GFN-FF path
 * where m_au = 1.
 *
 * Claude Generated (August 2026)
 */

#pragma once

#include "forcefieldfunctions.h"
#include "gfnff_geometry.h"
#include "gfnff_par.h"
#include "src/core/global.h"

#include <Eigen/Dense>
#include <cmath>

namespace GFNFFUnitTerms {

/// Covalent radius used by the bend/torsion damping (ff_workspace_gfnff.cpp: rcov_scale 4/3).
inline double covalentRadius(int atomic_number)
{
    constexpr double rcov_scale = 4.0 / 3.0;
    if (atomic_number >= 1
        && atomic_number <= static_cast<int>(GFNFFParameters::covalent_rad_d3.size()))
        return GFNFFParameters::covalent_rad_d3[atomic_number - 1] * rcov_scale;
    return 1.0 * GFNFFParameters::gfnff_aatoau * rcov_scale;
}

/// damp = 1 / (1 + (r^2 / (atcut (rcov_a + rcov_b)^2))^2), and its derivative in the
/// Fortran convention (the factor 2 of d(r^2)/dx_i = 2 (x_i - x_j) is folded in).
inline void distanceDamping(double r_sq, double rcov_a, double rcov_b, double atcut,
    double& damp, double& ddamp)
{
    const double rcut = atcut * (rcov_a + rcov_b) * (rcov_a + rcov_b);
    double rr = r_sq / rcut;
    rr *= rr;
    const double denom = 1.0 + rr;
    damp = 1.0 / denom;
    ddamp = (r_sq > 1e-8) ? -4.0 * rr / (r_sq * denom * denom) : 0.0;
}

/**
 * @brief GFN-FF bond at unit force constant: E = exp(-alpha * dr^2).
 *
 * @param r0    equilibrium distance in Bohr — pass exactly what the workspace uses,
 *              i.e. the DYNAMIC r0 = (r0_base_i + cnfak_i cn_i + r0_base_j
 *              + cnfak_j cn_j + rabshift) * ff when coordination numbers are available,
 *              otherwise Bond::r0_ij
 * @param alpha bond exponent, already including the hydrogen-bridge modulation
 *              (1 - 0.1 hb_cn_H) * alpha when Bond::nr_hb >= 1
 */
inline double bondEnergy(const Eigen::Vector3d& ri, const Eigen::Vector3d& rj,
    double r0, double alpha, Eigen::Vector3d& gi, Eigen::Vector3d& gj, bool do_gradient)
{
    const Eigen::Vector3d d = ri - rj;
    const double r = d.norm();
    if (r < 1e-12)
        return 0.0;

    const double dr = r - r0;
    const double energy = std::exp(-alpha * dr * dr);

    if (do_gradient) {
        // dE/dr = -2 alpha dr E ; dr/dx_i = d/r
        const double dEdr = -2.0 * alpha * dr * energy;
        const Eigen::Vector3d unit = d / r;
        gi += dEdr * unit;
        gj -= dEdr * unit;
    }
    return energy;
}

/**
 * @brief GFN-FF angle at unit force constant, damped.
 *
 * E = (cos(theta) - cos(theta0))^2 * damp, or (theta - theta0)^2 * damp for a linear
 * reference angle. @p rb is the CENTRAL atom.
 */
inline double angleEnergy(const Eigen::Vector3d& ra, const Eigen::Vector3d& rb,
    const Eigen::Vector3d& rc, int za, int zb, int zc, double theta0,
    Eigen::Vector3d& ga, Eigen::Vector3d& gb, Eigen::Vector3d& gc, bool do_gradient)
{
    constexpr double linear_threshold = 1.0e-6;
    constexpr double atcuta = 0.595;
    const double pi = 3.14159265358979323846;

    Matrix derivate;
    Vector va = ra, vb = rb, vc = rc;
    double costheta = UFF::AngleBending(va, vb, vc, derivate, do_gradient);
    costheta = std::max(-1.0, std::min(1.0, costheta));
    const double theta = std::acos(costheta);

    double energy = 0.0, dedtheta = 0.0;
    if (std::abs(pi - theta0) < linear_threshold) {
        const double dtheta = theta - theta0;
        energy = dtheta * dtheta;
        dedtheta = 2.0 * dtheta;
    } else {
        const double costheta0 = std::cos(theta0);
        const double dcostheta = costheta - costheta0;
        energy = dcostheta * dcostheta;
        dedtheta = 2.0 * std::sin(theta) * (costheta0 - costheta);
    }

    const double r_ab_sq = (ra - rb).squaredNorm();
    const double r_cb_sq = (rc - rb).squaredNorm();
    double damp_ab = 0.0, ddamp_ab = 0.0, damp_cb = 0.0, ddamp_cb = 0.0;
    distanceDamping(r_ab_sq, covalentRadius(za), covalentRadius(zb), atcuta, damp_ab, ddamp_ab);
    distanceDamping(r_cb_sq, covalentRadius(zb), covalentRadius(zc), atcuta, damp_cb, ddamp_cb);
    // The workspace uses the -2*2*rr/... form for angles, identical to -4*rr/...
    const double damp = damp_ab * damp_cb;

    if (do_gradient) {
        const Eigen::Vector3d vab = ra - rb;
        const Eigen::Vector3d vcb = rc - rb;
        const Eigen::Vector3d term1 = energy * ddamp_ab * damp_cb * vab;
        const Eigen::Vector3d term2 = energy * ddamp_cb * damp_ab * vcb;

        const Eigen::Vector3d g_a = dedtheta * damp * Eigen::Vector3d(derivate.row(0));
        const Eigen::Vector3d g_b = dedtheta * damp * Eigen::Vector3d(derivate.row(1));
        const Eigen::Vector3d g_c = dedtheta * damp * Eigen::Vector3d(derivate.row(2));

        ga += g_a + term1;
        gb += g_b - term1 - term2;
        gc += g_c + term2;
    }
    return energy * damp;
}

/**
 * @brief GFN-FF proper torsion at unit barrier: E = (1 + cos(n(phi-phi0)+pi)) * damp.
 * @param g out: 4x3, rows are atoms i, j, k, l (accumulated)
 */
inline double torsionEnergy(const Eigen::Vector3d& ri, const Eigen::Vector3d& rj,
    const Eigen::Vector3d& rk, const Eigen::Vector3d& rl,
    int zi, int zj, int zk, int zl, double n, double phi0, bool is_nci,
    Eigen::Matrix<double, 4, 3>& g, bool do_gradient)
{
    constexpr double atcutt = 0.505;
    constexpr double atcutt_nci = 0.305;
    const double atcut = is_nci ? atcutt_nci : atcutt;

    Matrix derivate;
    const double phi = GFNFF_Geometry::calculateDihedralAngle(ri, rj, rk, rl, derivate, do_gradient);

    const double c1 = n * (phi - phi0) + M_PI;
    const double energy = 1.0 + std::cos(c1);

    const double r_ij_sq = (ri - rj).squaredNorm();
    const double r_jk_sq = (rj - rk).squaredNorm();
    const double r_kl_sq = (rk - rl).squaredNorm();

    double damp_ij = 0.0, ddamp_ij = 0.0, damp_jk = 0.0, ddamp_jk = 0.0, damp_kl = 0.0, ddamp_kl = 0.0;
    distanceDamping(r_ij_sq, covalentRadius(zi), covalentRadius(zj), atcut, damp_ij, ddamp_ij);
    distanceDamping(r_jk_sq, covalentRadius(zj), covalentRadius(zk), atcut, damp_jk, ddamp_jk);
    distanceDamping(r_kl_sq, covalentRadius(zk), covalentRadius(zl), atcut, damp_kl, ddamp_kl);
    const double damp = damp_ij * damp_jk * damp_kl;

    if (do_gradient) {
        const double dEdphi = -n * std::sin(c1) * damp;
        for (int r = 0; r < 4; ++r)
            g.row(r) += dEdphi * derivate.row(r);

        const Eigen::Vector3d vij = ri - rj;
        const Eigen::Vector3d vjk = rj - rk;
        const Eigen::Vector3d vkl = rk - rl;
        const Eigen::Vector3d t1 = energy * ddamp_ij * damp_jk * damp_kl * vij;
        const Eigen::Vector3d t2 = energy * damp_ij * ddamp_jk * damp_kl * vjk;
        const Eigen::Vector3d t3 = energy * damp_ij * damp_jk * ddamp_kl * vkl;

        g.row(0) += t1.transpose();
        g.row(1) += (-t1 + t2).transpose();
        g.row(2) += (-t2 + t3).transpose();
        g.row(3) += (-t3).transpose();
    }
    return energy * damp;
}

/**
 * @brief GFN-FF out-of-plane term at unit barrier.
 *
 * potential_type 0 : E = (1 - cos(omega)) * damp                 (planar sp2)
 * otherwise        : E = (cos(omega) - cos(omega0))^2 * damp     (saturated N)
 *
 * Atom layout follows Inversion: i = centre, j = nb1, k = nb2, l = nb3, and the damping
 * uses nb1 as the hub (NOT a star around the centre) — gfnff_engrad.F90:1356-1365.
 */
inline double outOfPlaneEnergy(const Eigen::Vector3d& r_center, const Eigen::Vector3d& r_nb1,
    const Eigen::Vector3d& r_nb2, const Eigen::Vector3d& r_nb3,
    int z_center, int z_nb1, int z_nb2, int z_nb3,
    int potential_type, double omega0,
    Eigen::Matrix<double, 4, 3>& g, bool do_gradient)
{
    Matrix derivate;
    const double omega = GFNFF_Geometry::calculateOutOfPlaneAngle(
        r_center, r_nb1, r_nb2, r_nb3, derivate, do_gradient);

    const double rcov_c = covalentRadius(z_center);
    const double rcov_1 = covalentRadius(z_nb1);
    const double rcov_2 = covalentRadius(z_nb2);
    const double rcov_3 = covalentRadius(z_nb3);

    const double rij_sq = (r_nb1 - r_center).squaredNorm();
    const double rjk_sq = (r_nb1 - r_nb2).squaredNorm();
    const double rjl_sq = (r_nb1 - r_nb3).squaredNorm();

    double damp_ij = 0.0, ddamp_ij = 0.0, damp_jk = 0.0, ddamp_jk = 0.0, damp_jl = 0.0, ddamp_jl = 0.0;
    distanceDamping(rij_sq, rcov_c, rcov_1, GFNFFParameters::atcutt, damp_ij, ddamp_ij);
    distanceDamping(rjk_sq, rcov_2, rcov_1, GFNFFParameters::atcutt, damp_jk, ddamp_jk);
    distanceDamping(rjl_sq, rcov_1, rcov_3, GFNFFParameters::atcutt, damp_jl, ddamp_jl);
    const double damp = damp_ij * damp_jk * damp_jl;

    double et = 0.0, dEdomega = 0.0;
    if (potential_type == 0) {
        et = 1.0 - std::cos(omega);
        dEdomega = std::sin(omega) * damp;
    } else {
        const double diff = std::cos(omega) - std::cos(omega0);
        et = diff * diff;
        dEdomega = -2.0 * std::sin(omega) * diff * damp;
    }

    if (do_gradient) {
        for (int r = 0; r < 4; ++r)
            g.row(r) += dEdomega * derivate.row(r);

        const Eigen::Vector3d vab = r_nb1 - r_center;
        const Eigen::Vector3d vcb = r_nb1 - r_nb2;
        const Eigen::Vector3d vdc = r_nb1 - r_nb3;
        const Eigen::Vector3d t1 = et * ddamp_ij * damp_jk * damp_jl * vab;
        const Eigen::Vector3d t2 = et * damp_ij * ddamp_jk * damp_jl * vcb;
        const Eigen::Vector3d t3 = et * damp_ij * damp_jk * ddamp_jl * vdc;

        // nb1 is the hub of all three damping distances
        g.row(0) += (-t1).transpose();
        g.row(1) += (t1 + t2 + t3).transpose();
        g.row(2) += (-t2).transpose();
        g.row(3) += (-t3).transpose();
    }
    return et * damp;
}

} // namespace GFNFFUnitTerms
