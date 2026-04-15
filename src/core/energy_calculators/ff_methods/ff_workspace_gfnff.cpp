/*
 * <FFWorkspace GFN-FF Energy Term Calculators>
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
 * Claude Generated (March 2026): GFN-FF energy term calculators ported from
 * ForceFieldThread to FFWorkspace. Physics formulas are identical — only the
 * data access pattern changes (ranges on shared master lists, per-partition
 * accumulators instead of per-thread copies).
 *
 * Reference: Spicher/Grimme J. Chem. Theory Comput. 2020 (GFN-FF)
 * Reference: Fortran gfnff_engrad.F90 (energy and gradient formulas)
 */

#include "ff_workspace.h"
#include "gfnff_par.h"
#include "forcefieldfunctions.h"
#include "forcefieldderivaties.h"
#include "gfnff_geometry.h"
#include "src/core/units.h"
#include "src/core/curcuma_logger.h"

#include <fmt/core.h>
#include <fmt/format.h>

#include <cmath>
#include <unordered_map>

// ============================================================================
// HB Coordination Numbers (must run before bonds, only on partition 0)
// Reference: Fortran gfnff_engrad.F90:361, gfnff_data_types.f90:88
// ============================================================================

void FFWorkspace::computeHBCoordinationNumbers(int p)
{
    // Only partition 0 runs this — it modifies the shared m_bonds array,
    // but each bond.hb_cn_H is written by only one entry, so no race.
    if (p != 0) return;
    if (m_bond_hb_data.empty()) return;

    static const std::vector<double>& rcov_base = GFNFFParameters::covalent_rad_d3;
    constexpr double rcov_43 = 4.0 / 3.0;
    constexpr double kn = 27.5;
    constexpr double rcov_scal = 1.78;
    constexpr double thr = 900.0;

    std::unordered_map<int, double> hb_cn_map;
    m_hb_grad_entries.clear();
    constexpr double inv_sqrt_pi = 0.5641895835477563;

    for (const auto& entry : m_bond_hb_data) {
        int H = entry.H;
        int ati = m_atom_types[H];

        for (int B : entry.B_atoms) {
            int atj = m_atom_types[B];
            double dx = m_geometry(B, 0) - m_geometry(H, 0);
            double dy = m_geometry(B, 1) - m_geometry(H, 1);
            double dz = m_geometry(B, 2) - m_geometry(H, 2);
            double r2 = dx * dx + dy * dy + dz * dz;
            if (r2 > thr) continue;
            double r = std::sqrt(r2);

            double rcovij = rcov_scal * rcov_43 * (rcov_base[ati - 1] + rcov_base[atj - 1]);
            double arg = -kn * (r - rcovij) / rcovij;
            double tmp = 0.5 * (1.0 + std::erf(arg));
            hb_cn_map[H] += tmp;

            double dCN_dr = inv_sqrt_pi * (-kn / rcovij) * std::exp(-arg * arg) / r;
            Eigen::Vector3d r_HB(dx, dy, dz);
            m_hb_grad_entries.push_back({H, B, -dCN_dr * r_HB, dCN_dr * r_HB});
        }
    }

    for (auto& bond : m_bonds) {
        if (bond.nr_hb < 1) continue;
        int H = -1;
        if (m_atom_types[bond.i] == 1) H = bond.i;
        else if (m_atom_types[bond.j] == 1) H = bond.j;
        if (H >= 0) {
            auto it = hb_cn_map.find(H);
            bond.hb_cn_H = (it != hb_cn_map.end()) ? it->second : 0.0;
        }
    }
}

// ============================================================================
// Bond Stretching (GFN-FF exponential potential)
// Reference: Fortran gfnff_engrad.F90:675-721
// ============================================================================

void FFWorkspace::calcBonds(int p)
{
    auto& acc = m_accumulators[p];
    auto [begin, end] = m_partitions[p].bonds;
    Matrix grad_before;
    if (acc.has_components && m_do_gradient) grad_before = acc.gradient;

    bool use_dynamic_r0 = (m_d3_cn.size() > 0);

    for (int idx = begin; idx < end; ++idx) {
        const auto& bond = m_bonds[idx];

        Eigen::VectorXd vi = m_geometry.row(bond.i);
        Eigen::VectorXd vj = m_geometry.row(bond.j);
        Matrix derivate;
        double rij = UFF::BondStretching(vi, vj, derivate, m_do_gradient);

        double r0_ij;
        if (use_dynamic_r0 && bond.z_i > 0 && bond.z_j > 0 &&
            bond.i < m_d3_cn.size() && bond.j < m_d3_cn.size()) {
            double cn_i = m_d3_cn(bond.i);
            double cn_j = m_d3_cn(bond.j);
            double ra = bond.r0_base_i + bond.cnfak_i * cn_i;
            double rb = bond.r0_base_j + bond.cnfak_j * cn_j;
            r0_ij = (ra + rb + bond.rabshift) * bond.ff;
        } else {
            r0_ij = bond.r0_ij;
        }

        double dr = rij - r0_ij;
        double alpha_orig = bond.exponent;
        double alpha = alpha_orig;
        double k_b = bond.fc;

        // egbond_hb: Modified exponent for HB X-H bonds
        // Reference: Fortran gfnff_engrad.F90:957-958
        if (bond.nr_hb >= 1) {
            constexpr double VBOND_SCALE = 0.9;
            double t1 = 1.0 - VBOND_SCALE;
            alpha = (-t1 * bond.hb_cn_H + 1.0) * alpha_orig;
        }

        double exp_term = std::exp(-alpha * dr * dr);
        double energy = k_b * exp_term;
        acc.energy.bond += energy;

        if (m_do_gradient) {
            double dEdr = -2.0 * alpha * dr * energy;
            acc.gradient.row(bond.i) += dEdr * derivate.row(0);
            acc.gradient.row(bond.j) += dEdr * derivate.row(1);

            // HB alpha-modulation chain-rule gradient
            if (bond.nr_hb >= 1) {
                constexpr double t1 = 0.1;
                double zz = t1 * alpha_orig * dr * dr * energy;
                int H = (m_atom_types[bond.i] == 1) ? bond.i : bond.j;
                for (const auto& hbg : m_hb_grad_entries) {
                    if (hbg.H_atom != H) continue;
                    acc.gradient.row(H) += zz * hbg.dCN_dH.transpose();
                    acc.gradient.row(hbg.B_atom) += zz * hbg.dCN_dB.transpose();
                }
            }

            // Bond dr0/dCN chain-rule contribution
            if (use_dynamic_r0 && bond.z_i > 0 && bond.z_j > 0) {
                double yy = -dEdr;
                acc.dEdcn(bond.i) += yy * bond.ff * bond.cnfak_i;
                acc.dEdcn(bond.j) += yy * bond.ff * bond.cnfak_j;
                acc.dEdcn_bond(bond.i) += yy * bond.ff * bond.cnfak_i;
                acc.dEdcn_bond(bond.j) += yy * bond.ff * bond.cnfak_j;
            }
        }
    }

    if (acc.has_components && m_do_gradient)
        acc.grad_bond += (acc.gradient - grad_before);
}

// ============================================================================
// Angle Bending (GFN-FF cosine + distance damping)
// Reference: Fortran gfnff_engrad.F90:857-916
// ============================================================================

void FFWorkspace::calcAngles(int p)
{
    auto& acc = m_accumulators[p];
    auto [begin, end] = m_partitions[p].angles;
    Matrix grad_before;
    if (acc.has_components && m_do_gradient) grad_before = acc.gradient;

    const double pi = 3.14159265358979323846;
    const double linear_threshold = 1.0e-6;
    const double atcuta = 0.595;
    constexpr double rcov_scale_angle = 4.0 / 3.0;

    auto get_rcov_bohr = [&](int atomic_number) -> double {
        if (atomic_number >= 1 && atomic_number <= static_cast<int>(GFNFFParameters::covalent_rad_d3.size()))
            return GFNFFParameters::covalent_rad_d3[atomic_number - 1] * rcov_scale_angle;
        return 1.0 * GFNFFParameters::gfnff_aatoau * rcov_scale_angle;
    };

    for (int idx = begin; idx < end; ++idx) {
        const auto& angle = m_angles[idx];
        auto i = m_geometry.row(angle.i);
        auto j = m_geometry.row(angle.j);
        auto k = m_geometry.row(angle.k);
        Matrix derivate;
        double costheta = UFF::AngleBending(i, j, k, derivate, m_do_gradient);
        costheta = std::max(-1.0, std::min(1.0, costheta));

        double theta = std::acos(costheta);
        double theta0 = angle.theta0_ijk;
        double k_ijk = angle.fc;
        double energy, dedtheta;

        if (std::abs(pi - theta0) < linear_threshold) {
            double dtheta = theta - theta0;
            energy = k_ijk * dtheta * dtheta;
            dedtheta = 2.0 * k_ijk * dtheta;
        } else {
            double costheta0 = std::cos(theta0);
            double dcostheta = costheta - costheta0;
            energy = k_ijk * dcostheta * dcostheta;
            double sintheta = std::sin(theta);
            dedtheta = 2.0 * k_ijk * sintheta * (costheta0 - costheta);
        }

        double r_ij_sq = (i - j).squaredNorm();
        double r_jk_sq = (k - j).squaredNorm();

        double rcov_i = get_rcov_bohr(m_atom_types[angle.i]);
        double rcov_j = get_rcov_bohr(m_atom_types[angle.j]);
        double rcov_k = get_rcov_bohr(m_atom_types[angle.k]);

        double rcut_ij_sq = atcuta * (rcov_i + rcov_j) * (rcov_i + rcov_j);
        double rcut_jk_sq = atcuta * (rcov_j + rcov_k) * (rcov_j + rcov_k);

        double rr_ij = (r_ij_sq / rcut_ij_sq); rr_ij = rr_ij * rr_ij;
        double rr_jk = (r_jk_sq / rcut_jk_sq); rr_jk = rr_jk * rr_jk;

        double damp_ij = 1.0 / (1.0 + rr_ij);
        double damp_jk = 1.0 / (1.0 + rr_jk);
        double damp = damp_ij * damp_jk;

        double damp2ij = (r_ij_sq > 1e-8) ? -2.0 * 2.0 * rr_ij / (r_ij_sq * (1.0 + rr_ij) * (1.0 + rr_ij)) : 0.0;
        double damp2jk = (r_jk_sq > 1e-8) ? -2.0 * 2.0 * rr_jk / (r_jk_sq * (1.0 + rr_jk) * (1.0 + rr_jk)) : 0.0;

        acc.energy.angle += energy * damp;

        if (m_do_gradient) {
            // Claude Generated (Mar 2026): Use intermediate Vector to force column orientation
            // Mixing derivate.row() (1×3 row) with VectorXd (3×1 column) in += is undefined
            // in Eigen. Force all to column vectors first, matching ForceFieldThread pattern.
            Vector vab = (i - j).transpose();
            Vector vcb = (k - j).transpose();

            Vector term1 = energy * damp2ij * damp_jk * vab;
            Vector term2 = energy * damp2jk * damp_ij * vcb;

            Vector grad_i = dedtheta * damp * Vector(derivate.row(0));
            Vector grad_j = dedtheta * damp * Vector(derivate.row(1));
            Vector grad_k = dedtheta * damp * Vector(derivate.row(2));

            acc.gradient.row(angle.i) += (grad_i + term1).transpose();
            acc.gradient.row(angle.j) += (grad_j - term1 - term2).transpose();
            acc.gradient.row(angle.k) += (grad_k + term2).transpose();
        }
    }

    if (acc.has_components && m_do_gradient)
        acc.grad_angle += (acc.gradient - grad_before);
}

// ============================================================================
// Dihedral Torsion (GFN-FF with distance damping)
// Reference: Fortran gfnff_engrad.F90:1041-1122
// ============================================================================

void FFWorkspace::calcDihedrals(int p)
{
    auto& acc = m_accumulators[p];
    auto [begin, end] = m_partitions[p].dihedrals;
    Matrix grad_before;
    if (acc.has_components && m_do_gradient) grad_before = acc.gradient;

    const double atcutt = 0.505;
    const double atcutt_nci = 0.305;
    constexpr double rcov_scale_tors = 4.0 / 3.0;

    auto get_rcov_bohr = [&](int atomic_number) -> double {
        if (atomic_number >= 1 && atomic_number <= static_cast<int>(GFNFFParameters::covalent_rad_d3.size()))
            return GFNFFParameters::covalent_rad_d3[atomic_number - 1] * rcov_scale_tors;
        return 1.0 * GFNFFParameters::gfnff_aatoau * rcov_scale_tors;
    };

    for (int idx = begin; idx < end; ++idx) {
        const auto& dih = m_dihedrals[idx];

        Matrix derivate;
        double phi = GFNFF_Geometry::calculateDihedralAngle(
            m_geometry.row(dih.i).transpose(),
            m_geometry.row(dih.j).transpose(),
            m_geometry.row(dih.k).transpose(),
            m_geometry.row(dih.l).transpose(),
            derivate, m_do_gradient);

        double V = dih.V;
        double n = dih.n;
        double phi0 = dih.phi0;
        // Primary torsion: E = V * (1 + cos(n*(phi - phi0) + pi)) * damp
        // Reference: gfnff_engrad.F90:1268 — c1 = n*(phi-phi0) + pi
        double c1 = n * (phi - phi0) + M_PI;
        double energy = V * (1.0 + std::cos(c1));

        // Three-bond distance damping
        Eigen::Vector3d ri = m_geometry.row(dih.i).transpose();
        Eigen::Vector3d rj = m_geometry.row(dih.j).transpose();
        Eigen::Vector3d rk = m_geometry.row(dih.k).transpose();
        Eigen::Vector3d rl = m_geometry.row(dih.l).transpose();

        double r_ij_sq = (ri - rj).squaredNorm();
        double r_jk_sq = (rj - rk).squaredNorm();
        double r_kl_sq = (rk - rl).squaredNorm();

        double atcut = dih.is_nci ? atcutt_nci : atcutt;

        double rcov_i = get_rcov_bohr(m_atom_types[dih.i]);
        double rcov_j = get_rcov_bohr(m_atom_types[dih.j]);
        double rcov_k = get_rcov_bohr(m_atom_types[dih.k]);
        double rcov_l = get_rcov_bohr(m_atom_types[dih.l]);

        double rcut_ij = atcut * (rcov_i + rcov_j) * (rcov_i + rcov_j);
        double rcut_jk = atcut * (rcov_j + rcov_k) * (rcov_j + rcov_k);
        double rcut_kl = atcut * (rcov_k + rcov_l) * (rcov_k + rcov_l);

        double rr_ij = (r_ij_sq / rcut_ij); rr_ij *= rr_ij;
        double rr_jk = (r_jk_sq / rcut_jk); rr_jk *= rr_jk;
        double rr_kl = (r_kl_sq / rcut_kl); rr_kl *= rr_kl;

        double damp_ij = 1.0 / (1.0 + rr_ij);
        double damp_jk = 1.0 / (1.0 + rr_jk);
        double damp_kl = 1.0 / (1.0 + rr_kl);
        double damp = damp_ij * damp_jk * damp_kl;

        acc.energy.dihedral += energy * damp;

        if (m_do_gradient) {
            double dEdphi = -V * n * std::sin(c1) * damp;

            // derivate.row() is 1×3: no .transpose() needed for row() +=
            acc.gradient.row(dih.i) += dEdphi * derivate.row(0);
            acc.gradient.row(dih.j) += dEdphi * derivate.row(1);
            acc.gradient.row(dih.k) += dEdphi * derivate.row(2);
            acc.gradient.row(dih.l) += dEdphi * derivate.row(3);

            // Damping gradient terms
            double damp2ij = (r_ij_sq > 1e-8) ? -4.0 * rr_ij / (r_ij_sq * (1.0 + rr_ij) * (1.0 + rr_ij)) : 0.0;
            double damp2jk = (r_jk_sq > 1e-8) ? -4.0 * rr_jk / (r_jk_sq * (1.0 + rr_jk) * (1.0 + rr_jk)) : 0.0;
            double damp2kl = (r_kl_sq > 1e-8) ? -4.0 * rr_kl / (r_kl_sq * (1.0 + rr_kl) * (1.0 + rr_kl)) : 0.0;

            Eigen::Vector3d vij = ri - rj;
            Eigen::Vector3d vjk = rj - rk;
            Eigen::Vector3d vkl = rk - rl;

            Eigen::Vector3d t1 = energy * damp2ij * damp_jk * damp_kl * vij;
            Eigen::Vector3d t2 = energy * damp_ij * damp2jk * damp_kl * vjk;
            Eigen::Vector3d t3 = energy * damp_ij * damp_jk * damp2kl * vkl;

            // .transpose() converts Vector3d (3×1) to RowVector (1×3) for row() +=
            acc.gradient.row(dih.i) += t1.transpose();
            acc.gradient.row(dih.j) += (-t1 + t2).transpose();
            acc.gradient.row(dih.k) += (-t2 + t3).transpose();
            acc.gradient.row(dih.l) += (-t3).transpose();
        }
    }

    if (acc.has_components && m_do_gradient)
        acc.grad_torsion += (acc.gradient - grad_before);
}

// ============================================================================
// Extra Torsions (sp3-sp3 gauche, same formula)
// ============================================================================

void FFWorkspace::calcExtraTorsions(int p)
{
    auto& acc = m_accumulators[p];
    auto [begin, end] = m_partitions[p].extra_dihedrals;
    if (begin == end) return;

    Matrix grad_before;
    if (acc.has_components && m_do_gradient) grad_before = acc.gradient;

    const double atcutt = 0.505;
    constexpr double rcov_scale_tors = 4.0 / 3.0;

    auto get_rcov_bohr = [&](int atomic_number) -> double {
        if (atomic_number >= 1 && atomic_number <= static_cast<int>(GFNFFParameters::covalent_rad_d3.size()))
            return GFNFFParameters::covalent_rad_d3[atomic_number - 1] * rcov_scale_tors;
        return 1.0 * GFNFFParameters::gfnff_aatoau * rcov_scale_tors;
    };

    for (int idx = begin; idx < end; ++idx) {
        const auto& dih = m_extra_dihedrals[idx];

        Matrix derivate;
        double phi = GFNFF_Geometry::calculateDihedralAngle(
            m_geometry.row(dih.i).transpose(),
            m_geometry.row(dih.j).transpose(),
            m_geometry.row(dih.k).transpose(),
            m_geometry.row(dih.l).transpose(),
            derivate, m_do_gradient);

        double V = dih.V;
        double n = dih.n;
        double phi0 = dih.phi0;
        // Extra torsion: same formula as primary — c1 = n*(phi-phi0) + π
        // Reference: gfnff_engrad.F90:1268 — ALL torsions with tlist(5)>0 use +π
        double c1 = n * (phi - phi0) + M_PI;
        double energy = V * (1.0 + std::cos(c1));

        Eigen::Vector3d ri = m_geometry.row(dih.i).transpose();
        Eigen::Vector3d rj = m_geometry.row(dih.j).transpose();
        Eigen::Vector3d rk = m_geometry.row(dih.k).transpose();
        Eigen::Vector3d rl = m_geometry.row(dih.l).transpose();

        double r_ij_sq = (ri - rj).squaredNorm();
        double r_jk_sq = (rj - rk).squaredNorm();
        double r_kl_sq = (rk - rl).squaredNorm();

        double rcov_i = get_rcov_bohr(m_atom_types[dih.i]);
        double rcov_j = get_rcov_bohr(m_atom_types[dih.j]);
        double rcov_k = get_rcov_bohr(m_atom_types[dih.k]);
        double rcov_l = get_rcov_bohr(m_atom_types[dih.l]);

        double rcut_ij = atcutt * (rcov_i + rcov_j) * (rcov_i + rcov_j);
        double rcut_jk = atcutt * (rcov_j + rcov_k) * (rcov_j + rcov_k);
        double rcut_kl = atcutt * (rcov_k + rcov_l) * (rcov_k + rcov_l);

        double rr_ij = (r_ij_sq / rcut_ij); rr_ij *= rr_ij;
        double rr_jk = (r_jk_sq / rcut_jk); rr_jk *= rr_jk;
        double rr_kl = (r_kl_sq / rcut_kl); rr_kl *= rr_kl;

        double damp = (1.0 / (1.0 + rr_ij)) * (1.0 / (1.0 + rr_jk)) * (1.0 / (1.0 + rr_kl));

        acc.energy.dihedral += energy * damp;

        if (m_do_gradient) {
            double damp_ij = 1.0 / (1.0 + rr_ij);
            double damp_jk = 1.0 / (1.0 + rr_jk);
            double damp_kl = 1.0 / (1.0 + rr_kl);

            double dEdphi = -V * n * std::sin(c1) * damp;
            // derivate.row() is 1×3: no .transpose() needed for row() +=
            acc.gradient.row(dih.i) += dEdphi * derivate.row(0);
            acc.gradient.row(dih.j) += dEdphi * derivate.row(1);
            acc.gradient.row(dih.k) += dEdphi * derivate.row(2);
            acc.gradient.row(dih.l) += dEdphi * derivate.row(3);

            double damp2ij = (r_ij_sq > 1e-8) ? -4.0 * rr_ij / (r_ij_sq * (1.0 + rr_ij) * (1.0 + rr_ij)) : 0.0;
            double damp2jk = (r_jk_sq > 1e-8) ? -4.0 * rr_jk / (r_jk_sq * (1.0 + rr_jk) * (1.0 + rr_jk)) : 0.0;
            double damp2kl = (r_kl_sq > 1e-8) ? -4.0 * rr_kl / (r_kl_sq * (1.0 + rr_kl) * (1.0 + rr_kl)) : 0.0;

            Eigen::Vector3d t1 = energy * damp2ij * damp_jk * damp_kl * (ri - rj);
            Eigen::Vector3d t2 = energy * damp_ij * damp2jk * damp_kl * (rj - rk);
            Eigen::Vector3d t3 = energy * damp_ij * damp_jk * damp2kl * (rk - rl);

            // .transpose() converts Vector3d (3×1) to RowVector (1×3) for row() +=
            acc.gradient.row(dih.i) += t1.transpose();
            acc.gradient.row(dih.j) += (-t1 + t2).transpose();
            acc.gradient.row(dih.k) += (-t2 + t3).transpose();
            acc.gradient.row(dih.l) += (-t3).transpose();
        }
    }

    if (acc.has_components && m_do_gradient)
        acc.grad_torsion += (acc.gradient - grad_before);
}

// ============================================================================
// Inversions (out-of-plane bending)
// Reference: ForceFieldThread::CalculateGFNFFInversionContribution
// Uses domegadr analytical derivatives from gfnff_inversions.cpp
// ============================================================================

void FFWorkspace::calcInversions(int p)
{
    // Claude Generated (Mar 2026): Complete inversion with gradient
    // Ported from ForceFieldThread::CalculateGFNFFInversionContribution
    // Reference: gfnff_engrad.F90:1355-1387
    auto& acc = m_accumulators[p];
    auto [begin, end] = m_partitions[p].inversions;
    if (begin == end) return;

    Matrix grad_before;
    if (acc.has_components && m_do_gradient) grad_before = acc.gradient;

    constexpr double rcov_scale = 4.0 / 3.0;

    auto get_rcov = [&](int at) -> double {
        if (at >= 1 && at <= static_cast<int>(GFNFFParameters::covalent_rad_d3.size()))
            return GFNFFParameters::covalent_rad_d3[at - 1] * rcov_scale;
        return 1.0 * GFNFFParameters::gfnff_aatoau * rcov_scale;
    };

    auto calc_damp = [](double rsq, double rcov_a, double rcov_b) -> double {
        double rcut = GFNFFParameters::atcutt * (rcov_a + rcov_b) * (rcov_a + rcov_b);
        double rr = (rsq / rcut) * (rsq / rcut);
        return 1.0 / (1.0 + rr);
    };

    auto calc_ddamp = [](double r2_val, double rcov_a, double rcov_b) -> double {
        if (r2_val < 1e-8) return 0.0;
        double rcut_val = GFNFFParameters::atcutt * (rcov_a + rcov_b) * (rcov_a + rcov_b);
        double rr_val = (r2_val / rcut_val) * (r2_val / rcut_val);
        double one_plus_rr = 1.0 + rr_val;
        return -4.0 * rr_val / (r2_val * one_plus_rr * one_plus_rr);
    };

    for (int idx = begin; idx < end; ++idx) {
        const auto& inv = m_inversions[idx];

        // Atom layout: i=center, j=nb1, k=nb2, l=nb3
        Eigen::Vector3d r_center = m_geometry.row(inv.i).transpose();
        Eigen::Vector3d r_nb1 = m_geometry.row(inv.j).transpose();
        Eigen::Vector3d r_nb2 = m_geometry.row(inv.k).transpose();
        Eigen::Vector3d r_nb3 = m_geometry.row(inv.l).transpose();

        // Out-of-plane angle via calculateOutOfPlaneAngle
        Matrix derivate;
        double omega = GFNFF_Geometry::calculateOutOfPlaneAngle(
            r_center, r_nb1, r_nb2, r_nb3, derivate, m_do_gradient);

        // Damping: nb1 as hub (NOT star topology)
        // Reference: gfnff_engrad.F90:1356-1365
        double rcov_c = get_rcov(m_atom_types[inv.i]);
        double rcov_1 = get_rcov(m_atom_types[inv.j]);
        double rcov_2 = get_rcov(m_atom_types[inv.k]);
        double rcov_3 = get_rcov(m_atom_types[inv.l]);

        double rij_sq = (r_nb1 - r_center).squaredNorm();  // nb1-center: bond
        double rjk_sq = (r_nb1 - r_nb2).squaredNorm();     // nb1-nb2: 1-3 distance
        double rjl_sq = (r_nb1 - r_nb3).squaredNorm();     // nb1-nb3: 1-3 distance

        double damp_ij = calc_damp(rij_sq, rcov_c, rcov_1);
        double damp_jk = calc_damp(rjk_sq, rcov_2, rcov_1);
        double damp_jl = calc_damp(rjl_sq, rcov_1, rcov_3);
        double damp = damp_ij * damp_jk * damp_jl;

        // Energy and dE/domega
        double V = inv.fc;
        double et = 0.0;
        double dEdomega = 0.0;

        if (inv.potential_type == 0) {
            // Planar sp2: E = V*(1 - cos(omega)) * damp
            et = V * (1.0 - std::cos(omega));
            dEdomega = V * std::sin(omega) * damp;
        } else {
            // Saturated N: E = V*(cos(omega) - cos(omega0))^2 * damp
            double diff = std::cos(omega) - std::cos(inv.omega0);
            et = V * diff * diff;
            dEdomega = -2.0 * V * std::sin(omega) * diff * damp;
        }

        acc.energy.inversion += et * damp;

        if (m_do_gradient) {
            // PART 1: Omega derivative
            // Note: derivate.row() is 1×3; no transpose needed for row() += row()
            acc.gradient.row(inv.i) += dEdomega * derivate.row(0);
            acc.gradient.row(inv.j) += dEdomega * derivate.row(1);
            acc.gradient.row(inv.k) += dEdomega * derivate.row(2);
            acc.gradient.row(inv.l) += dEdomega * derivate.row(3);

            // PART 2: Damping gradient (nb1 as hub)
            // Reference: gfnff_engrad.F90:1379-1385
            double ddamp_ij = calc_ddamp(rij_sq, rcov_c, rcov_1);
            double ddamp_jk = calc_ddamp(rjk_sq, rcov_2, rcov_1);
            double ddamp_jl = calc_ddamp(rjl_sq, rcov_1, rcov_3);

            Eigen::Vector3d vab = r_nb1 - r_center;  // j - i
            Eigen::Vector3d vcb = r_nb1 - r_nb2;     // j - k
            Eigen::Vector3d vdc = r_nb1 - r_nb3;     // j - l

            Eigen::Vector3d term1 = (et * ddamp_ij * damp_jk * damp_jl) * vab;
            Eigen::Vector3d term2 = (et * ddamp_jk * damp_ij * damp_jl) * vcb;
            Eigen::Vector3d term3 = (et * ddamp_jl * damp_ij * damp_jk) * vdc;

            // center: -term1, nb1(hub): +term1+term2+term3, nb2: -term2, nb3: -term3
            // Note: .transpose() converts Vector3d (3×1) to RowVector (1×3) for row() +=
            acc.gradient.row(inv.i) -= term1.transpose();
            acc.gradient.row(inv.j) += (term1 + term2 + term3).transpose();
            acc.gradient.row(inv.k) -= term2.transpose();
            acc.gradient.row(inv.l) -= term3.transpose();
        }
    }

    if (acc.has_components && m_do_gradient)
        acc.grad_torsion += (acc.gradient - grad_before);
}

// ============================================================================
// Triple Bond Torsions (sTors_eg)
// Reference: Fortran gfnff_engrad.F90:3454
// ============================================================================

void FFWorkspace::calcSTorsions(int p)
{
    auto& acc = m_accumulators[p];
    auto [begin, end] = m_partitions[p].storsions;
    if (begin == end) return;

    Matrix grad_before;
    if (acc.has_components && m_do_gradient) grad_before = acc.gradient;

    for (int idx = begin; idx < end; ++idx) {
        const auto& stor = m_storsions[idx];
        Matrix derivate;
        // GFN-FF uses Bohr coordinates internally (m_au = 1.0)
        double phi = GFNFF_Geometry::calculateDihedralAngle(
            m_geometry.row(stor.i).transpose(),
            m_geometry.row(stor.j).transpose(),
            m_geometry.row(stor.k).transpose(),
            m_geometry.row(stor.l).transpose(),
            derivate, m_do_gradient);

        double erefhalf = stor.erefhalf;
        double energy = -erefhalf * std::cos(2.0 * phi) + erefhalf;
        acc.energy.stors += energy;

        if (m_do_gradient) {
            double dEdphi = 2.0 * erefhalf * std::sin(2.0 * phi);
            // derivate.row() is 1×3: no .transpose() needed for row() +=
            acc.gradient.row(stor.i) += dEdphi * derivate.row(0);
            acc.gradient.row(stor.j) += dEdphi * derivate.row(1);
            acc.gradient.row(stor.k) += dEdphi * derivate.row(2);
            acc.gradient.row(stor.l) += dEdphi * derivate.row(3);
        }
    }

    // Track sTorsion gradient in torsion component (matches ForceFieldThread)
    if (acc.has_components && m_do_gradient)
        acc.grad_torsion += (acc.gradient - grad_before);
}

// ============================================================================
// Dispersion (GFN-FF modified BJ damping)
// Reference: gfnff_gdisp0.f90:365-377
// ============================================================================

void FFWorkspace::calcDispersion(int p)
{
    auto& acc = m_accumulators[p];
    auto [begin, end] = m_partitions[p].dispersions;
    if (begin == end) return;

    Matrix grad_before;
    if (acc.has_components && m_do_gradient) grad_before = acc.gradient;

    for (int idx = begin; idx < end; ++idx) {
        const auto& disp = m_dispersions[idx];

        Eigen::Vector3d ri = m_geometry.row(disp.i);
        Eigen::Vector3d rj = m_geometry.row(disp.j);
        Eigen::Vector3d rij_vec = ri - rj;
        double rij = rij_vec.norm();

        if (rij > disp.r_cut || rij < 1e-8) continue;

        double r2 = rij * rij;
        double r6 = r2 * r2 * r2;
        double r0_6 = disp.r0_squared * disp.r0_squared * disp.r0_squared;
        double t6 = 1.0 / (r6 + r0_6);
        double r8 = r6 * r2;
        double r0_8 = r0_6 * disp.r0_squared;
        double t8 = 1.0 / (r8 + r0_8);

        double disp_sum = t6 + 2.0 * disp.r4r2ij * t8;
        double energy = -disp.C6 * disp_sum * disp.zetac6;
        acc.energy.dispersion += energy;

        if (m_do_gradient) {
            double d6 = -6.0 * r2 * r2 * t6 * t6;
            double d8 = -8.0 * r2 * r2 * r2 * t8 * t8;
            double ddisp_dr2 = d6 + 2.0 * disp.r4r2ij * d8;
            double dEdr = -disp.C6 * disp.zetac6 * ddisp_dr2 * rij;

            Eigen::Vector3d grad = dEdr * rij_vec / rij;
            acc.gradient.row(disp.i) += grad.transpose();
            acc.gradient.row(disp.j) -= grad.transpose();

            // dc6dcn chain-rule
            if (m_dc6dcn_ptr && m_dc6dcn_ptr->size() > 0 &&
                disp.i < m_dc6dcn_ptr->rows() && disp.j < m_dc6dcn_ptr->cols()) {
                double disp_value = disp_sum * disp.zetac6;
                acc.dEdcn(disp.i) -= (*m_dc6dcn_ptr)(disp.i, disp.j) * disp_value;
                acc.dEdcn(disp.j) -= (*m_dc6dcn_ptr)(disp.j, disp.i) * disp_value;
            }
        }
    }

    if (acc.has_components && m_do_gradient)
        acc.grad_dispersion += (acc.gradient - grad_before);
}

// ============================================================================
// D4 Dispersion (same formula, separate parameter list)
// Reference: gfnff_gdisp0.f90:365-377
// ============================================================================

void FFWorkspace::calcD4Dispersion(int p)
{
    auto& acc = m_accumulators[p];
    auto [begin, end] = m_partitions[p].d4_dispersions;
    if (begin == end) return;

    Matrix grad_before;
    if (acc.has_components && m_do_gradient) grad_before = acc.gradient;

    for (int idx = begin; idx < end; ++idx) {
        const auto& disp = m_d4_dispersions[idx];

        Eigen::Vector3d ri = m_geometry.row(disp.i);
        Eigen::Vector3d rj = m_geometry.row(disp.j);
        Eigen::Vector3d rij_vec = ri - rj;
        double rij = rij_vec.norm();

        if (rij > disp.r_cut || rij < 1e-10) continue;

        double r2 = rij * rij;
        double r6 = r2 * r2 * r2;
        double r0_6 = disp.r0_squared * disp.r0_squared * disp.r0_squared;
        double t6 = 1.0 / (r6 + r0_6);
        double r8 = r6 * r2;
        double r0_8 = r0_6 * disp.r0_squared;
        double t8 = 1.0 / (r8 + r0_8);

        double disp_sum = t6 + 2.0 * disp.r4r2ij * t8;
        double pair_energy = -disp.C6 * disp_sum * disp.zetac6;
        acc.energy.dispersion += pair_energy;

        if (m_do_gradient) {
            double d6 = -6.0 * r2 * r2 * t6 * t6;
            double d8 = -8.0 * r2 * r2 * r2 * t8 * t8;
            double ddisp_dr2 = d6 + 2.0 * disp.r4r2ij * d8;
            double dEdr = -disp.C6 * disp.zetac6 * ddisp_dr2 * rij;

            Eigen::Vector3d grad = dEdr * rij_vec / rij;
            acc.gradient.row(disp.i) += grad.transpose();
            acc.gradient.row(disp.j) -= grad.transpose();

            if (m_dc6dcn_ptr && m_dc6dcn_ptr->size() > 0 &&
                disp.i < m_dc6dcn_ptr->rows() && disp.j < m_dc6dcn_ptr->cols()) {
                double disp_value = disp_sum * disp.zetac6;
                acc.dEdcn(disp.i) -= (*m_dc6dcn_ptr)(disp.i, disp.j) * disp_value;
                acc.dEdcn(disp.j) -= (*m_dc6dcn_ptr)(disp.j, disp.i) * disp_value;
            }
        }
    }

    if (acc.has_components && m_do_gradient)
        acc.grad_dispersion += (acc.gradient - grad_before);
}

// ============================================================================
// Bonded Repulsion — E = repab * exp(-α*r^1.5) / r
// Reference: Fortran gfnff_engrad.F90:467-495
// ============================================================================

void FFWorkspace::calcBondedRepulsion(int p)
{
    auto& acc = m_accumulators[p];
    auto [begin, end] = m_partitions[p].bonded_reps;
    if (begin == end) return;

    Matrix grad_before;
    if (acc.has_components && m_do_gradient) grad_before = acc.gradient;

    for (int idx = begin; idx < end; ++idx) {
        const auto& rep = m_bonded_reps[idx];

        Eigen::Vector3d ri = m_geometry.row(rep.i);
        Eigen::Vector3d rj = m_geometry.row(rep.j);
        Eigen::Vector3d rij_vec = ri - rj;
        double rij = rij_vec.norm();
        if (rij > rep.r_cut || rij < 1e-8) continue;

        double r_1_5 = rij * std::sqrt(rij);
        double exp_term = std::exp(-rep.alpha * r_1_5);
        double base_energy = rep.repab * exp_term / rij;
        acc.energy.bonded_rep += base_energy;

        if (m_do_gradient) {
            double dEdr = (-base_energy / rij - 1.5 * rep.alpha * std::sqrt(rij) * base_energy);
            Eigen::Vector3d grad = dEdr * rij_vec / rij;
            acc.gradient.row(rep.i) += grad.transpose();
            acc.gradient.row(rep.j) -= grad.transpose();
        }
    }

    if (acc.has_components && m_do_gradient)
        acc.grad_repulsion += (acc.gradient - grad_before);
}

// ============================================================================
// Non-bonded Repulsion — same formula, different parameter set
// Reference: Fortran gfnff_engrad.F90:255-276
// ============================================================================

void FFWorkspace::calcNonbondedRepulsion(int p)
{
    auto& acc = m_accumulators[p];
    auto [begin, end] = m_partitions[p].nonbonded_reps;
    if (begin == end) return;
    // Non-bonded repulsion goes to total gradient only, NOT to g_rep component

    for (int idx = begin; idx < end; ++idx) {
        const auto& rep = m_nonbonded_reps[idx];

        Eigen::Vector3d ri = m_geometry.row(rep.i);
        Eigen::Vector3d rj = m_geometry.row(rep.j);
        Eigen::Vector3d rij_vec = ri - rj;
        double rij = rij_vec.norm();
        if (rij > rep.r_cut || rij < 1e-8) continue;

        double r_1_5 = rij * std::sqrt(rij);
        double exp_term = std::exp(-rep.alpha * r_1_5);
        double base_energy = rep.repab * exp_term / rij;
        acc.energy.nonbonded_rep += base_energy;

        if (m_do_gradient) {
            double dEdr = (-base_energy / rij - 1.5 * rep.alpha * std::sqrt(rij) * base_energy);
            Eigen::Vector3d grad = dEdr * rij_vec / rij;
            acc.gradient.row(rep.i) += grad.transpose();
            acc.gradient.row(rep.j) -= grad.transpose();
        }
    }
}

// ============================================================================
// Coulomb (TERM 1: pairwise electrostatics)
// TERM 2+3 (self-energy) handled in postProcess()
// Reference: Fortran gfnff_engrad.F90:1378-1389
// ============================================================================

void FFWorkspace::calcCoulomb(int p)
{
    auto& acc = m_accumulators[p];
    auto [begin, end] = m_partitions[p].coulombs;
    if (begin == end) return;

    Matrix grad_before;
    if (acc.has_components && m_do_gradient) grad_before = acc.gradient;

    for (int idx = begin; idx < end; ++idx) {
        const auto& coul = m_coulombs[idx];

        Eigen::Vector3d ri = m_geometry.row(coul.i);
        Eigen::Vector3d rj = m_geometry.row(coul.j);
        Eigen::Vector3d rij_vec = ri - rj;
        double rij = rij_vec.norm();
        if (rij > coul.r_cut || rij < 1e-10) continue;

        // Dynamic EEQ charges (fall back to static if unavailable/NaN)
        double qi = coul.q_i, qj = coul.q_j;
        if (m_eeq_charges.size() > 0) {
            double qi_dyn = m_eeq_charges(coul.i);
            double qj_dyn = m_eeq_charges(coul.j);
            if (!std::isnan(qi_dyn) && !std::isnan(qj_dyn)) {
                qi = qi_dyn;
                qj = qj_dyn;
            }
        }

        double gamma_r = coul.gamma_ij * rij;
        double erf_term = std::erf(gamma_r);
        double energy_pair = qi * qj * erf_term / rij;
        acc.energy.coulomb += energy_pair;

        if (m_do_gradient) {
            const double sqrt_pi = 1.772453850905516;
            double exp_term = std::exp(-gamma_r * gamma_r);
            double derf_dr = coul.gamma_ij * exp_term * (2.0 / sqrt_pi);
            double dEdr_pair = qi * qj * (derf_dr / rij - erf_term / (rij * rij));

            Eigen::Vector3d grad = dEdr_pair * rij_vec / rij;
            acc.gradient.row(coul.i) += grad.transpose();
            acc.gradient.row(coul.j) -= grad.transpose();
        }
    }

    if (acc.has_components && m_do_gradient)
        acc.grad_coulomb += (acc.gradient - grad_before);
}

// ============================================================================
// HB/XB Damping Helper Functions (file-local)
// ============================================================================

namespace {

inline double ws_damping_out_of_line(double r_AH, double r_HB, double r_AB, double radab, double bacut)
{
    double ratio = (r_AH + r_HB) / r_AB;
    double exponent = (bacut / radab) * (ratio - 1.0);
    if (exponent > 15.0) return 0.0;
    return 2.0 / (1.0 + std::exp(exponent));
}

inline double ws_damping_short_range(double r, double r_vdw, double scut, double alp)
{
    double ratio = scut * r_vdw / (r * r);
    return 1.0 / (1.0 + std::pow(ratio, alp));
}

inline double ws_damping_long_range(double r, double longcut, double alp)
{
    return 1.0 / (1.0 + std::pow(r * r / longcut, alp));
}

inline double ws_charge_scaling(double q, double st, double sf)
{
    double exp_term = std::exp(st * q);
    return exp_term / (exp_term + sf);
}

} // anonymous namespace

// ============================================================================
// Hydrogen Bonds (three-body A-H...B)
// Reference: gfnff_engrad.F90 - abhgfnff_eg1, abhgfnff_eg2new, abhgfnff_eg3
// ============================================================================

void FFWorkspace::calcHydrogenBonds(int p)
{
    auto& acc = m_accumulators[p];
    auto [begin, end] = m_partitions[p].hbonds;
    if (begin == end) return;

    Matrix grad_before;
    if (acc.has_components && m_do_gradient) grad_before = acc.gradient;

    using namespace GFNFFParameters;

    for (int idx = begin; idx < end; ++idx) {
        const auto& hb = m_hbonds[idx];

        Eigen::Vector3d pos_A = m_geometry.row(hb.i).transpose();
        Eigen::Vector3d pos_H = m_geometry.row(hb.j).transpose();
        Eigen::Vector3d pos_B = m_geometry.row(hb.k).transpose();

        double r_AH = (pos_H - pos_A).norm();
        double r_HB = (pos_B - pos_H).norm();
        double r_AB = (pos_B - pos_A).norm();

        if (r_AB > hb.r_cut) continue;

        double r_AH_4 = r_AH * r_AH * r_AH * r_AH;
        double r_HB_4 = r_HB * r_HB * r_HB * r_HB;
        double denom_DA = 1.0 / (r_AH_4 + r_HB_4);

        // HBond uses Phase 1 topology charges from struct (same as ForceFieldThread)
        // NOT Phase 2 EEQ charges — Fortran gfnff_engrad uses nhb1/nhb2 list charges
        double Q_H = ws_charge_scaling(hb.q_H, HB_ST, HB_SF);
        double Q_A = ws_charge_scaling(-hb.q_A, HB_ST, HB_SF);
        double Q_B = ws_charge_scaling(-hb.q_B, HB_ST, HB_SF);

        double bas = (Q_A * hb.basicity_A * r_AH_4 + Q_B * hb.basicity_B * r_HB_4) * denom_DA;
        double aci = (hb.acidity_B * r_AH_4 + hb.acidity_A * r_HB_4) * denom_DA;

        int elem_A = m_atom_types[hb.i];
        int elem_B = m_atom_types[hb.k];
        double r_vdw_AB = covalent_radii[elem_A - 1] + covalent_radii[elem_B - 1];

        double damp_short = ws_damping_short_range(r_AB, r_vdw_AB, HB_SCUT, HB_ALP);
        double damp_long = ws_damping_long_range(r_AB, HB_LONGCUT, HB_ALP);
        double damp_env = damp_short * damp_long;

        double damp_outl = ws_damping_out_of_line(r_AH, r_HB, r_AB, r_vdw_AB, HB_BACUT);

        double outl_nb_tot = 1.0;
        if (hb.case_type >= 2) {
            double hbnbcut_save = (elem_B == 7 && hb.neighbors_B.size() == 1) ? 2.0 : HB_NBCUT;
            for (int nb : hb.neighbors_B) {
                Eigen::Vector3d pos_nb = m_geometry.row(nb).transpose();
                double r_Anb = (pos_nb - pos_A).norm();
                double r_Bnb = (pos_nb - pos_B).norm();
                double expo_nb = (hbnbcut_save / r_vdw_AB) * ((r_Anb + r_Bnb) / r_AB - 1.0);
                outl_nb_tot *= (2.0 / (1.0 + std::exp(-expo_nb)) - 1.0);
            }
        }

        double rdamp;
        if (hb.case_type >= 2) {
            rdamp = damp_env * (1.8 / (r_HB * r_HB * r_HB) - 0.8 / (r_AB * r_AB * r_AB));
        } else {
            rdamp = damp_env / (r_AB * r_AB * r_AB);
        }

        // Case 4: Virtual LP (with gradient variables)
        // Claude Generated (Mar 2026): Port from ForceFieldThread with LP gradient support
        double outl_lp = 1.0;
        Eigen::Vector3d lp_pos = pos_B;
        double lp_dist = 0.0;
        int nbb_lp = 0;
        Eigen::Vector3d lp_vector = Eigen::Vector3d::Zero();
        if (hb.case_type == 4) {
            int z_B = m_atom_types[hb.k];
            double repz_B = (z_B >= 1 && z_B <= static_cast<int>(GFNFFParameters::repz.size()))
                          ? GFNFFParameters::repz[z_B - 1] : 1.0;
            lp_dist = 0.50 - 0.018 * repz_B;
            static constexpr double HBLPCUT = 56.0;

            nbb_lp = static_cast<int>(hb.neighbors_B.size());
            for (int nb_idx : hb.neighbors_B) {
                lp_vector += (m_geometry.row(nb_idx).transpose() - pos_B);
            }
            double vnorm = lp_vector.norm();
            if (vnorm > 1e-10 && nbb_lp > 0) {
                lp_pos = pos_B + (-lp_dist) * (lp_vector / vnorm);
                double ralp = (pos_A - lp_pos).norm();
                double expo_lp_val = (HBLPCUT / r_vdw_AB) * ((ralp + lp_dist + 1e-12) / r_AB - 1.0);
                outl_lp = 2.0 / (1.0 + std::exp(expo_lp_val));
            } else {
                nbb_lp = 0;
                outl_lp = 1.0;
            }
        }

        double qhoutl = Q_H * damp_outl * outl_nb_tot * outl_lp;

        // Case 3: eangl + etors (with gradient data)
        // Claude Generated (Mar 2026): Full port from ForceFieldThread with gradient
        // Reference: gfnff_engrad.F90:2574-2583 (egbend_nci_mul), 2529-2557 (egtors_nci_mul)
        double eangl = 1.0, etors = 1.0;
        struct TorsGrad {
            double energy;
            Eigen::Matrix<double, 4, 3> grad;
            int D_idx;
        };
        std::vector<TorsGrad> tors_data;
        Eigen::Matrix<double, 3, 3> gangl_3body = Eigen::Matrix<double, 3, 3>::Zero();
        int jj_idx = -1, kk_idx = -1, ll_idx = -1;

        if (hb.case_type == 3 && hb.acceptor_parent_index != -1) {
            int B_idx = hb.k;
            int C_idx = hb.acceptor_parent_index;
            int H_idx = hb.j;
            jj_idx = B_idx;
            kk_idx = C_idx;
            ll_idx = H_idx;

            // eangl: angle bending H...B=C (matching ForceFieldThread exactly)
            {
                double c0 = 120.0 * M_PI / 180.0;
                double fc_bend = 1.0 - BEND_HB;
                double kijk = fc_bend / ((std::cos(0.0) - std::cos(c0)) * (std::cos(0.0) - std::cos(c0)));

                Eigen::Vector3d va = m_geometry.row(C_idx).transpose();
                Eigen::Vector3d vb = m_geometry.row(B_idx).transpose();
                Eigen::Vector3d vc = m_geometry.row(H_idx).transpose();

                Eigen::Vector3d vab = va - vb;
                Eigen::Vector3d vcb = vc - vb;
                double rab2_a = vab.squaredNorm();
                double rcb2_a = vcb.squaredNorm();
                Eigen::Vector3d vp = vcb.cross(vab);
                double rp = vp.norm() + 1e-14;

                double cosa = vab.dot(vcb) / (std::sqrt(rab2_a) * std::sqrt(rcb2_a) + 1e-14);
                cosa = std::clamp(cosa, -1.0, 1.0);
                double theta_a = std::acos(cosa);

                double ea, deddt;
                if (M_PI - c0 < 1e-6) {
                    double dt = theta_a - c0;
                    ea = kijk * dt * dt;
                    deddt = 2.0 * kijk * dt;
                } else {
                    ea = kijk * (cosa - std::cos(c0)) * (cosa - std::cos(c0));
                    deddt = 2.0 * kijk * std::sin(theta_a) * (std::cos(c0) - cosa);
                }
                eangl = 1.0 - ea;

                if (m_do_gradient) {
                    Eigen::Vector3d deda_v = vab.cross(vp) * (-deddt / (rab2_a * rp));
                    Eigen::Vector3d dedc_v = vcb.cross(vp) * (deddt / (rcb2_a * rp));
                    Eigen::Vector3d dedb_v = deda_v + dedc_v;
                    gangl_3body.row(0) = dedb_v.transpose();
                    gangl_3body.row(1) = -deda_v.transpose();
                    gangl_3body.row(2) = -dedc_v.transpose();
                }
            }

            // etors: product of torsion terms D-B-C-H (with gradient)
            for (int D_idx : hb.neighbors_C) {
                if (D_idx == B_idx) continue;
                Matrix dihedral_grad;
                double phi = GFNFF_Geometry::calculateDihedralAngle(
                    m_geometry.row(D_idx).transpose(),
                    m_geometry.row(B_idx).transpose(),
                    m_geometry.row(C_idx).transpose(),
                    m_geometry.row(H_idx).transpose(),
                    dihedral_grad, m_do_gradient);

                double tshift = TORS_HB;
                double fc_tors = (1.0 - tshift) / 2.0;
                double phi0_tors = M_PI / 2.0;
                int rn = 2;
                double dphi1 = phi - phi0_tors;
                double c1 = rn * dphi1 + M_PI;
                double et = (1.0 + std::cos(c1)) * fc_tors + tshift;
                double dij = -rn * std::sin(c1) * fc_tors;

                TorsGrad tg;
                tg.energy = et;
                tg.D_idx = D_idx;
                tg.grad = Eigen::Matrix<double, 4, 3>::Zero();
                if (m_do_gradient && dihedral_grad.rows() == 4) {
                    for (int a = 0; a < 4; ++a) {
                        tg.grad.row(a) = dij * dihedral_grad.row(a);
                    }
                }
                tors_data.push_back(tg);
            }

            etors = 1.0;
            for (const auto& tg : tors_data) {
                etors *= tg.energy;
            }
        }

        // Energy: case-specific formula (matching ForceFieldThread)
        double global_scale = 1.0;
        if (hb.case_type == 2 || hb.case_type == 4) global_scale = XHACI_GLOBABH;
        else if (hb.case_type == 3) global_scale = XHACI_COH;

        double E_HB;
        if (hb.case_type >= 2) {
            double const_val = hb.acidity_A * hb.basicity_B * Q_A * Q_B * global_scale;
            E_HB = -rdamp * qhoutl * const_val * eangl * etors;
        } else {
            E_HB = -bas * aci * rdamp * qhoutl;
        }
        acc.energy.hbond += E_HB;

        // ========== ANALYTICAL GRADIENT CALCULATION ==========
        // Claude Generated (Mar 2026): Complete port from ForceFieldThread
        // Reference: gfnff_engrad.F90 - abhgfnff_eg1/eg2new/eg2_rnr/eg3
        if (m_do_gradient) {
            // Fortran convention distance vectors (all in Bohr)
            Eigen::Vector3d drab = pos_A - pos_B;  // A - B
            Eigen::Vector3d drah = pos_A - pos_H;  // A - H
            Eigen::Vector3d drbh = pos_B - pos_H;  // B - H

            double rab2 = r_AB * r_AB;
            double rbh2 = r_HB * r_HB;
            double rah2 = r_AH * r_AH;
            double rahprbh = r_AH + r_HB + 1e-12;

            // Damping derivative intermediates
            double ratio1 = std::pow(rab2 / HB_LONGCUT, HB_ALP);
            double shortcut = HB_SCUT * r_vdw_AB;
            double ratio3 = std::pow(shortcut / rab2, HB_ALP);
            double ddamp = (-2.0 * HB_ALP * ratio1 / (1.0 + ratio1))
                         + ( 2.0 * HB_ALP * ratio3 / (1.0 + ratio3));

            // Out-of-line intermediates for gradient
            double expo = (HB_BACUT / r_vdw_AB) * (rahprbh / r_AB - 1.0);
            if (expo > 15.0) continue;  // Fortran early return (gfnff_engrad.F90:1781)
            double ratio2 = std::exp(expo);

            Eigen::Vector3d ga = Eigen::Vector3d::Zero();
            Eigen::Vector3d gb = Eigen::Vector3d::Zero();
            Eigen::Vector3d gh = Eigen::Vector3d::Zero();

            if (hb.case_type >= 2) {
                // ===== Case 2/3/4: abhgfnff_eg2new gradient =====
                double p_bh = 1.8;   // 1 + hbabmix
                double p_ab = -0.8;  // -hbabmix
                double rbhdamp = damp_env * p_bh / (rbh2 * r_HB);
                double rabdamp = damp_env * p_ab / (rab2 * r_AB);

                double const_val = hb.acidity_A * hb.basicity_B * Q_A * Q_B * global_scale;
                double dterm  = -qhoutl * eangl * etors * const_val;
                double aterm  = -rdamp * Q_H * outl_nb_tot * outl_lp * eangl * etors * const_val;
                double nbterm = -rdamp * Q_H * damp_outl * outl_lp * eangl * etors * const_val;

                // Damping part: rab
                double gi = ((rabdamp + rbhdamp) * ddamp - 3.0 * rabdamp) / rab2;
                gi *= dterm;
                Eigen::Vector3d dg = gi * drab;
                ga = dg;
                gb = -dg;

                // Damping part: rbh
                gi = -3.0 * rbhdamp / rbh2;
                gi *= dterm;
                dg = gi * drbh;
                gb += dg;
                gh = -dg;

                // Out-of-line: rab
                double tmp1 = -2.0 * aterm * ratio2 * expo
                            / ((1.0 + ratio2) * (1.0 + ratio2))
                            / (rahprbh - r_AB);
                gi = -tmp1 * rahprbh / rab2;
                dg = gi * drab;
                ga += dg;
                gb -= dg;

                // Out-of-line: rah, rbh
                gi = tmp1 / r_AH;
                Eigen::Vector3d dga_outl = gi * drah;
                ga += dga_outl;
                gi = tmp1 / r_HB;
                Eigen::Vector3d dgb_outl = gi * drbh;
                gb += dgb_outl;
                gh += -dga_outl - dgb_outl;

                // Neighbor out-of-line gradient (Case >= 2)
                double hbnbcut_g = (elem_B == 7 && hb.neighbors_B.size() == 1) ? 2.0 : HB_NBCUT;
                for (size_t nb_i = 0; nb_i < hb.neighbors_B.size(); ++nb_i) {
                    int nb = hb.neighbors_B[nb_i];
                    Eigen::Vector3d pos_nb = m_geometry.row(nb).transpose();
                    Eigen::Vector3d dranb = pos_A - pos_nb;
                    Eigen::Vector3d drbnb = pos_B - pos_nb;
                    double ranb = dranb.norm();
                    double rbnb = drbnb.norm();
                    double ranbprbnb = ranb + rbnb + 1e-12;

                    double expo_nb_i = (hbnbcut_g / r_vdw_AB) * (ranbprbnb / r_AB - 1.0);
                    double ratio2_nb_i = std::exp(-expo_nb_i);
                    double outl_nb_i = 2.0 / (1.0 + ratio2_nb_i) - 1.0;

                    double outl_nb_others = 1.0;
                    if (std::abs(outl_nb_i) > 1e-12) {
                        outl_nb_others = outl_nb_tot / outl_nb_i;
                    }

                    double tmp2 = 2.0 * nbterm * outl_nb_others * ratio2_nb_i * expo_nb_i
                                / ((1.0 + ratio2_nb_i) * (1.0 + ratio2_nb_i))
                                / (ranbprbnb - r_AB);

                    double gi_nb = -tmp2 * ranbprbnb / rab2;
                    dg = gi_nb * drab;
                    ga += dg;
                    gb -= dg;

                    gi_nb = tmp2 / ranb;
                    Eigen::Vector3d dga_nb = gi_nb * dranb;
                    ga += dga_nb;
                    gi_nb = tmp2 / rbnb;
                    Eigen::Vector3d dgb_nb = gi_nb * drbnb;
                    gb += dgb_nb;
                    acc.gradient.row(nb) += (-dga_nb - dgb_nb).transpose();
                }

                // Case 4: LP out-of-line gradient
                if (hb.case_type == 4 && nbb_lp > 0) {
                    double lpterm = -rdamp * Q_H * damp_outl * outl_nb_tot * const_val;
                    double ralp = (pos_A - lp_pos).norm();
                    double rblp = lp_dist;
                    double ralpprblp = ralp + rblp + 1e-12;
                    static constexpr double HBLPCUT = 56.0;
                    double expo_lp_val = (HBLPCUT / r_vdw_AB) * (ralpprblp / r_AB - 1.0);
                    double ratio2_lp_val = std::exp(expo_lp_val);

                    // LP out-of-line: rab
                    double tmp3 = -2.0 * lpterm * ratio2_lp_val * expo_lp_val
                                / ((1.0 + ratio2_lp_val) * (1.0 + ratio2_lp_val))
                                / (ralpprblp - r_AB);
                    double gi_lp = -tmp3 * ralpprblp / rab2;
                    Eigen::Vector3d dg_lp = gi_lp * drab;
                    ga += dg_lp;
                    gb -= dg_lp;

                    // LP out-of-line: ralp
                    Eigen::Vector3d dralp = pos_A - lp_pos;
                    gi_lp = tmp3 / (ralp + 1e-12);
                    Eigen::Vector3d dga_lp = gi_lp * dralp;
                    ga += dga_lp;

                    // Fortran: gb -= dga (uses dga, not dgb)
                    gb -= dga_lp;
                    Eigen::Vector3d glp = -dga_lp;

                    // LP neighbor chain rule
                    double vnorm = lp_vector.norm();
                    if (vnorm > 1e-10) {
                        Eigen::Matrix3d gii = Eigen::Matrix3d::Zero();
                        for (int col = 0; col < 3; ++col) {
                            Eigen::Vector3d unit_vec = Eigen::Vector3d::Zero();
                            unit_vec(col) = -1.0;
                            gii.col(col) = -lp_dist * static_cast<double>(nbb_lp)
                                         * (unit_vec / vnorm + lp_vector * lp_vector(col) / std::pow(vnorm, 3.0));
                        }
                        Eigen::Vector3d gnb_lp = gii * glp;
                        gb += gnb_lp;
                        Eigen::Vector3d gnb_lp_share = gnb_lp / static_cast<double>(nbb_lp);
                        for (int nb_idx : hb.neighbors_B) {
                            acc.gradient.row(nb_idx) -= gnb_lp_share.transpose();
                        }
                    }
                }

                // Case 3: angle bending and torsion gradient contributions
                if (hb.case_type == 3 && jj_idx >= 0) {
                    double bterm_c3 = -rdamp * qhoutl * etors * const_val;
                    double tterm_c3 = -rdamp * qhoutl * eangl * const_val;

                    acc.gradient.row(jj_idx) += bterm_c3 * gangl_3body.row(0);
                    acc.gradient.row(kk_idx) += bterm_c3 * gangl_3body.row(1);
                    acc.gradient.row(ll_idx) += bterm_c3 * gangl_3body.row(2);

                    for (size_t k = 0; k < tors_data.size(); ++k) {
                        double factor_k = (std::abs(tors_data[k].energy) > 1e-12)
                                        ? etors / tors_data[k].energy : 0.0;
                        double t_k = factor_k * tterm_c3;
                        acc.gradient.row(tors_data[k].D_idx) += t_k * tors_data[k].grad.row(0);
                        acc.gradient.row(jj_idx)             += t_k * tors_data[k].grad.row(1);
                        acc.gradient.row(kk_idx)             += t_k * tors_data[k].grad.row(2);
                        acc.gradient.row(ll_idx)             += t_k * tors_data[k].grad.row(3);
                    }
                }
            } else {
                // ===== Case 1: abhgfnff_eg1 gradient =====
                double caa = Q_A * hb.basicity_A;
                double cbb = Q_B * hb.basicity_B;

                double rterm = -aci * rdamp * qhoutl;
                double dterm = -aci * bas * qhoutl;
                double sterm = -rdamp * bas * qhoutl;
                double aterm = -aci * bas * rdamp * Q_H;

                double denom_val = 1.0 / (r_AH_4 + r_HB_4);
                double tmp = denom_val * denom_val * 4.0;
                double dd24a = rah2 * r_HB_4 * tmp;
                double dd24b = rbh2 * r_AH_4 * tmp;

                // Donor-acceptor: bas
                double gi = (caa - cbb) * dd24a * rterm;
                ga = gi * drah;
                gi = (cbb - caa) * dd24b * rterm;
                gb = gi * drbh;
                gh = -ga - gb;

                // Donor-acceptor: aci
                gi = (hb.acidity_B - hb.acidity_A) * dd24a;
                Eigen::Vector3d dga_aci = gi * drah * sterm;
                ga += dga_aci;
                gi = (hb.acidity_A - hb.acidity_B) * dd24b;
                Eigen::Vector3d dgb_aci = gi * drbh * sterm;
                gb += dgb_aci;
                gh += -dga_aci - dgb_aci;

                // Damping: rab
                gi = rdamp * (ddamp - 3.0) / rab2;
                Eigen::Vector3d dg = gi * drab * dterm;
                ga += dg;
                gb -= dg;

                // Out-of-line: rab
                gi = aterm * 2.0 * ratio2 * expo * rahprbh
                   / ((1.0 + ratio2) * (1.0 + ratio2))
                   / (rahprbh - r_AB) / rab2;
                dg = gi * drab;
                ga += dg;
                gb -= dg;

                // Out-of-line: rah, rbh
                double tmp_outl = -2.0 * aterm * ratio2 * expo
                                / ((1.0 + ratio2) * (1.0 + ratio2))
                                / (rahprbh - r_AB);
                Eigen::Vector3d dga_outl = drah * tmp_outl / r_AH;
                ga += dga_outl;
                Eigen::Vector3d dgb_outl = drbh * tmp_outl / r_HB;
                gb += dgb_outl;
                gh += -dga_outl - dgb_outl;
            }

            // Accumulate: A=hb.i, B=hb.k, H=hb.j
            acc.gradient.row(hb.i) += ga.transpose();
            acc.gradient.row(hb.k) += gb.transpose();
            acc.gradient.row(hb.j) += gh.transpose();
        }
    }

    if (acc.has_components && m_do_gradient)
        acc.grad_hb += (acc.gradient - grad_before);
}

// ============================================================================
// Halogen Bonds (three-body A-X...B)
// Reference: gfnff_engrad.F90 - rbxgfnff_eg
// ============================================================================

void FFWorkspace::calcHalogenBonds(int p)
{
    auto& acc = m_accumulators[p];
    auto [begin, end] = m_partitions[p].xbonds;
    if (begin == end) return;

    Matrix grad_before;
    if (acc.has_components && m_do_gradient) grad_before = acc.gradient;

    using namespace GFNFFParameters;

    for (int idx = begin; idx < end; ++idx) {
        const auto& xb = m_xbonds[idx];

        Eigen::Vector3d pos_A = m_geometry.row(xb.i).transpose();
        Eigen::Vector3d pos_X = m_geometry.row(xb.j).transpose();
        Eigen::Vector3d pos_B = m_geometry.row(xb.k).transpose();

        Eigen::Vector3d r_AX_vec = pos_X - pos_A;
        Eigen::Vector3d r_XB_vec = pos_B - pos_X;
        Eigen::Vector3d r_AB_vec = pos_B - pos_A;

        double r_AX = r_AX_vec.norm();
        double r_XB = r_XB_vec.norm();
        double r_AB = r_AB_vec.norm();

        if (r_XB > xb.r_cut) continue;

        // Claude Generated (Mar 2026): Fixed to match ForceFieldThread exactly
        int elem_A = m_atom_types[xb.i];
        int elem_B = m_atom_types[xb.k];
        double r_vdw_AB = covalent_radii[elem_A - 1] + covalent_radii[elem_B - 1];

        double damp_short = ws_damping_short_range(r_XB, r_vdw_AB, XB_SCUT, HB_ALP);
        double damp_long = ws_damping_long_range(r_XB, HB_LONGCUT_XB, HB_ALP);

        // XB out-of-line: uses xbacut directly, not divided by radab
        double ratio_outl = (r_AX + r_XB) / r_AB;
        double expo_outl = XB_BACUT * (ratio_outl - 1.0);
        if (expo_outl > 15.0) continue;  // Early exit matching Fortran
        double damp_outl = 2.0 / (1.0 + std::exp(expo_outl));

        double Q_X = ws_charge_scaling(xb.q_X, XB_ST, XB_SF);
        double Q_B = ws_charge_scaling(-xb.q_B, XB_ST, XB_SF);

        // Energy: R_damp = damp_short * damp_long * damp_outl / r_XB^3
        // (Fixed: was r_AB^3, now r_XB^3 matching ForceFieldThread/Fortran)
        double R_damp = damp_short * damp_long * damp_outl / (r_XB * r_XB * r_XB);
        double E_XB = -R_damp * Q_B * xb.acidity_X * Q_X;
        acc.energy.xbond += E_XB;

        // ========== ANALYTICAL GRADIENT CALCULATION ==========
        // Claude Generated (Mar 2026): Port from ForceFieldThread
        // Reference: gfnff_engrad.F90 - rbxgfnff_eg()
        if (m_do_gradient) {
            // Short-range damping derivative w.r.t. r_XB
            double ratio_short = XB_SCUT * r_vdw_AB / (r_XB * r_XB);
            double damp_short_term = std::pow(ratio_short, HB_ALP);
            double ddamp_short_dr = -2.0 * HB_ALP * damp_short * damp_short_term
                                  / (r_XB * (1.0 + damp_short_term));

            // Long-range damping derivative w.r.t. r_XB
            double ratio_long = (r_XB * r_XB) / HB_LONGCUT_XB;
            double damp_long_term = std::pow(ratio_long, HB_ALP);
            double ddamp_long_dr = -2.0 * HB_ALP * r_XB * damp_long * damp_long_term
                                 / (HB_LONGCUT_XB * (1.0 + damp_long_term));

            // Out-of-line damping derivatives
            double exp_term = std::exp(expo_outl);
            double denom_outl = 1.0 + exp_term;
            double ddamp_outl_drAX = -2.0 * exp_term * XB_BACUT
                                   / (r_AB * denom_outl * denom_outl);
            double ddamp_outl_drXB = ddamp_outl_drAX;
            double ddamp_outl_drAB = 2.0 * exp_term * XB_BACUT * (r_AX + r_XB)
                                   / (r_AB * r_AB * denom_outl * denom_outl);

            // R_damp chain rule derivatives
            double dRdamp_drXB = (ddamp_short_dr * damp_long * damp_outl
                                + damp_short * ddamp_long_dr * damp_outl
                                + damp_short * damp_long * ddamp_outl_drXB) / (r_XB * r_XB * r_XB)
                               - 3.0 * R_damp / r_XB;
            double dRdamp_drAX = damp_short * damp_long * ddamp_outl_drAX / (r_XB * r_XB * r_XB);
            double dRdamp_drAB = damp_short * damp_long * ddamp_outl_drAB / (r_XB * r_XB * r_XB);

            double E_prefactor = -Q_B * xb.acidity_X * Q_X;
            double dE_drXB = E_prefactor * dRdamp_drXB;
            double dE_drAX = E_prefactor * dRdamp_drAX;
            double dE_drAB = E_prefactor * dRdamp_drAB;

            Eigen::Vector3d grad_rAX_unit = r_AX_vec / (r_AX + 1e-14);
            Eigen::Vector3d grad_rXB_unit = r_XB_vec / (r_XB + 1e-14);
            Eigen::Vector3d grad_rAB_unit = r_AB_vec / (r_AB + 1e-14);

            // A: dr_AX/dA = -unit_AX, dr_AB/dA = -unit_AB
            Eigen::Vector3d grad_A = -dE_drAX * grad_rAX_unit - dE_drAB * grad_rAB_unit;
            // X: dr_AX/dX = +unit_AX, dr_XB/dX = -unit_XB
            Eigen::Vector3d grad_X = dE_drAX * grad_rAX_unit - dE_drXB * grad_rXB_unit;
            // B: dr_XB/dB = +unit_XB, dr_AB/dB = +unit_AB
            Eigen::Vector3d grad_B = dE_drXB * grad_rXB_unit + dE_drAB * grad_rAB_unit;

            acc.gradient.row(xb.i) += grad_A.transpose();
            acc.gradient.row(xb.j) += grad_X.transpose();
            acc.gradient.row(xb.k) += grad_B.transpose();
        }
    }

    if (acc.has_components && m_do_gradient)
        acc.grad_xb += (acc.gradient - grad_before);
}

// ============================================================================
// ATM Three-Body Dispersion (Axilrod-Teller-Muto)
// Reference: external/cpp-d4/src/damping/atm.cpp:70-138
// ============================================================================

void FFWorkspace::calcATM(int p)
{
    auto& acc = m_accumulators[p];
    auto [begin, end] = m_partitions[p].atm_triples;
    if (begin == end) return;

    Matrix grad_before;
    if (acc.has_components && m_do_gradient) grad_before = acc.gradient;

    using namespace GFNFFParameters;

    for (int idx = begin; idx < end; ++idx) {
        const auto& triple = m_atm_triples[idx];

        Eigen::Vector3d pos_i = m_geometry.row(triple.i).transpose();
        Eigen::Vector3d pos_j = m_geometry.row(triple.j).transpose();
        Eigen::Vector3d pos_k = m_geometry.row(triple.k).transpose();

        double rij = (pos_i - pos_j).norm();
        double rik = (pos_i - pos_k).norm();
        double rjk = (pos_j - pos_k).norm();

        double r2ij = rij * rij, r2ik = rik * rik, r2jk = rjk * rjk;

        double c9 = triple.s9 * std::sqrt(std::fabs(triple.C6_ij * triple.C6_ik * triple.C6_jk));

        int zi = m_atom_types[triple.i], zj = m_atom_types[triple.j], zk = m_atom_types[triple.k];
        double r_cov_i = (zi > 0 && zi <= static_cast<int>(rcov_bohr.size())) ? rcov_bohr[zi - 1] : 1.0;
        double r_cov_j = (zj > 0 && zj <= static_cast<int>(rcov_bohr.size())) ? rcov_bohr[zj - 1] : 1.0;
        double r_cov_k = (zk > 0 && zk <= static_cast<int>(rcov_bohr.size())) ? rcov_bohr[zk - 1] : 1.0;

        double r0ij = triple.a1 * std::sqrt(3.0 * r_cov_i * r_cov_j) + triple.a2;
        double r0ik = triple.a1 * std::sqrt(3.0 * r_cov_i * r_cov_k) + triple.a2;
        double r0jk = triple.a1 * std::sqrt(3.0 * r_cov_j * r_cov_k) + triple.a2;

        double rijk = rij * rik * rjk;
        double r2ijk = r2ij * r2ik * r2jk;
        double r3ijk = rijk * r2ijk;

        double fdmp = 1.0 / (1.0 + 6.0 * std::pow(r0ij * r0ik * r0jk / rijk, triple.alp / 3.0));

        double A = r2ij + r2jk - r2ik;
        double B = r2ij + r2ik - r2jk;
        double C = r2ik + r2jk - r2ij;
        double ang = (0.375 * A * B * C / r2ijk + 1.0) / r3ijk;

        acc.energy.atm += ang * fdmp * c9 / 3.0 * triple.triple_scale;
    }

    if (acc.has_components && m_do_gradient)
        acc.grad_atm += (acc.gradient - grad_before);
}

// ============================================================================
// ATM Gradient (analytical)
// Reference: external/cpp-d4/src/damping/atm.cpp:141-289
// ============================================================================

void FFWorkspace::calcATMGradient(int p)
{
    auto& acc = m_accumulators[p];
    auto [begin, end] = m_partitions[p].atm_triples;
    if (begin == end) return;

    using namespace GFNFFParameters;

    for (int idx = begin; idx < end; ++idx) {
        const auto& triple = m_atm_triples[idx];

        Eigen::Vector3d pos_i = m_geometry.row(triple.i).transpose();
        Eigen::Vector3d pos_j = m_geometry.row(triple.j).transpose();
        Eigen::Vector3d pos_k = m_geometry.row(triple.k).transpose();

        Eigen::Vector3d rij_vec = pos_j - pos_i;
        Eigen::Vector3d rik_vec = pos_k - pos_i;
        Eigen::Vector3d rjk_vec = pos_k - pos_j;

        double rij = rij_vec.norm(), rik = rik_vec.norm(), rjk = rjk_vec.norm();
        double r2ij = rij * rij, r2ik = rik * rik, r2jk = rjk * rjk;

        // Negative c9 for gradient (dispersion is attractive)
        double c9 = -triple.s9 * std::sqrt(std::fabs(triple.C6_ij * triple.C6_ik * triple.C6_jk));

        int zi = m_atom_types[triple.i], zj = m_atom_types[triple.j], zk = m_atom_types[triple.k];
        double r_cov_i = (zi > 0 && zi <= static_cast<int>(rcov_bohr.size())) ? rcov_bohr[zi - 1] : 1.0;
        double r_cov_j = (zj > 0 && zj <= static_cast<int>(rcov_bohr.size())) ? rcov_bohr[zj - 1] : 1.0;
        double r_cov_k = (zk > 0 && zk <= static_cast<int>(rcov_bohr.size())) ? rcov_bohr[zk - 1] : 1.0;

        double r0ij = triple.a1 * std::sqrt(3.0 * r_cov_i * r_cov_j) + triple.a2;
        double r0ik = triple.a1 * std::sqrt(3.0 * r_cov_i * r_cov_k) + triple.a2;
        double r0jk = triple.a1 * std::sqrt(3.0 * r_cov_j * r_cov_k) + triple.a2;
        double r0ijk = r0ij * r0ik * r0jk;

        double rijk = rij * rik * rjk;
        double r2ijk = r2ij * r2ik * r2jk;
        double r3ijk = rijk * r2ijk;
        double r5ijk = r2ijk * r3ijk;

        double tmp = std::pow(r0ijk / rijk, triple.alp / 3.0);
        double fdmp = 1.0 / (1.0 + 6.0 * tmp);
        double dfdmp = -2.0 * triple.alp * tmp * fdmp * fdmp;

        double A = r2ij + r2jk - r2ik;
        double B = r2ij + r2ik - r2jk;
        double C = r2ik + r2jk - r2ij;
        double ang = (0.375 * A * B * C / r2ijk + 1.0) / r3ijk;

        double dang_ij = -0.375 * (std::pow(r2ij, 3) + std::pow(r2ij, 2) * (r2jk + r2ik)
                          + r2ij * (3.0 * std::pow(r2jk, 2) + 2.0 * r2jk * r2ik + 3.0 * std::pow(r2ik, 2))
                          - 5.0 * std::pow(r2jk - r2ik, 2) * (r2jk + r2ik)) / r5ijk;

        double dang_ik = -0.375 * (std::pow(r2ik, 3) + std::pow(r2ik, 2) * (r2jk + r2ij)
                          + r2ik * (3.0 * std::pow(r2jk, 2) + 2.0 * r2jk * r2ij + 3.0 * std::pow(r2ij, 2))
                          - 5.0 * std::pow(r2jk - r2ij, 2) * (r2jk + r2ij)) / r5ijk;

        double dang_jk = -0.375 * (std::pow(r2jk, 3) + std::pow(r2jk, 2) * (r2ik + r2ij)
                          + r2jk * (3.0 * std::pow(r2ik, 2) + 2.0 * r2ik * r2ij + 3.0 * std::pow(r2ij, 2))
                          - 5.0 * std::pow(r2ik - r2ij, 2) * (r2ik + r2ij)) / r5ijk;

        double prefactor = c9 * triple.triple_scale / 3.0;
        Eigen::Vector3d dgij = prefactor * (-dang_ij * fdmp + ang * dfdmp) / r2ij * rij_vec;
        Eigen::Vector3d dgik = prefactor * (-dang_ik * fdmp + ang * dfdmp) / r2ik * rik_vec;
        Eigen::Vector3d dgjk = prefactor * (-dang_jk * fdmp + ang * dfdmp) / r2jk * rjk_vec;

        acc.gradient.row(triple.i) += -(dgij + dgik);
        acc.gradient.row(triple.j) += (dgij - dgjk);
        acc.gradient.row(triple.k) += (dgik + dgjk);
    }
}

// ============================================================================
// Bonded ATM (BATM) for 1,4-pairs
// Reference: gfnff_engrad.F90:3267-3334 (batmgfnff_eg)
// ============================================================================

void FFWorkspace::calcBATM(int p)
{
    auto& acc = m_accumulators[p];
    auto [begin, end] = m_partitions[p].batm_triples;
    if (begin == end) return;

    Matrix grad_before;
    if (acc.has_components && m_do_gradient) grad_before = acc.gradient;

    const double fqq = 3.0;

    for (int idx = begin; idx < end; ++idx) {
        const auto& batm = m_batm_triples[idx];

        Eigen::Vector3d i_pos = m_geometry.row(batm.i).transpose();
        Eigen::Vector3d j_pos = m_geometry.row(batm.j).transpose();
        Eigen::Vector3d k_pos = m_geometry.row(batm.k).transpose();

        Eigen::Vector3d rij_vec = j_pos - i_pos;
        Eigen::Vector3d rik_vec = k_pos - i_pos;
        Eigen::Vector3d rjk_vec = k_pos - j_pos;

        double r2ij = rij_vec.squaredNorm();
        double r2jk = rjk_vec.squaredNorm();
        double r2ik = rik_vec.squaredNorm();

        double rij = std::sqrt(r2ij), rjk = std::sqrt(r2jk), rik = std::sqrt(r2ik);

        double rijk3 = r2ij * r2jk * r2ik;
        double mijk = -r2ij + r2jk + r2ik;
        double imjk = r2ij - r2jk + r2ik;
        double ijmk = r2ij + r2jk - r2ik;

        double ang = 0.375 * ijmk * imjk * mijk / rijk3;
        double rav3 = std::pow(rijk3, 1.5);
        double angr9 = (ang + 1.0) / rav3;

        // Phase-1 topology charges (fixed)
        double fi = std::min(std::max(1.0 - fqq * m_topology_charges(batm.i), -4.0), 4.0);
        double fj = std::min(std::max(1.0 - fqq * m_topology_charges(batm.j), -4.0), 4.0);
        double fk = std::min(std::max(1.0 - fqq * m_topology_charges(batm.k), -4.0), 4.0);

        double c9 = fi * fj * fk * batm.zb3atm_i * batm.zb3atm_j * batm.zb3atm_k;
        double energy = c9 * angr9;
        acc.energy.batm += energy;

        if (m_do_gradient) {
            double dang_ij = -0.375 * (std::pow(r2ij, 3) + std::pow(r2ij, 2) * (r2jk + r2ik)
                                  + r2ij * (3.0 * std::pow(r2jk, 2) + 2.0 * r2jk * r2ik + 3.0 * std::pow(r2ik, 2))
                                  - 5.0 * std::pow(r2jk - r2ik, 2) * (r2jk + r2ik))
                                  / (rij * rijk3 * rav3);

            double dang_jk = -0.375 * (std::pow(r2jk, 3) + std::pow(r2jk, 2) * (r2ik + r2ij)
                                  + r2jk * (3.0 * std::pow(r2ik, 2) + 2.0 * r2ik * r2ij + 3.0 * std::pow(r2ij, 2))
                                  - 5.0 * std::pow(r2ik - r2ij, 2) * (r2ik + r2ij))
                                  / (rjk * rijk3 * rav3);

            double dang_ik = -0.375 * (std::pow(r2ik, 3) + std::pow(r2ik, 2) * (r2jk + r2ij)
                                  + r2ik * (3.0 * std::pow(r2jk, 2) + 2.0 * r2jk * r2ij + 3.0 * std::pow(r2ij, 2))
                                  - 5.0 * std::pow(r2jk - r2ij, 2) * (r2jk + r2ij))
                                  / (rik * rijk3 * rav3);

            Eigen::Vector3d dgij = -dang_ij * c9 * (rij_vec / rij);
            Eigen::Vector3d dgjk = -dang_jk * c9 * (rjk_vec / rjk);
            Eigen::Vector3d dgik = -dang_ik * c9 * (rik_vec / rik);

            acc.gradient.row(batm.j) += (-dgij + dgjk);
            acc.gradient.row(batm.k) += (-dgik - dgjk);
            acc.gradient.row(batm.i) += (dgij + dgik);
        }
    }

    if (acc.has_components && m_do_gradient)
        acc.grad_batm += (acc.gradient - grad_before);
}
