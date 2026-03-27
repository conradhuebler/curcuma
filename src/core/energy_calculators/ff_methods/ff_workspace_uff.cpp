/*
 * <FFWorkspace UFF/QMDFF Energy Term Calculators>
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
 * Claude Generated (March 2026): UFF/QMDFF energy term calculators ported from
 * ForceFieldThread to FFWorkspace. Physics formulas are identical — only the
 * data access pattern changes (ranges on shared master lists, per-partition
 * accumulators instead of per-thread copies).
 *
 * Unit convention: UFF/QMDFF geometry is in Angström.
 * m_au = 1.889726125 converts Å → Bohr for vdW distances.
 */

#include "ff_workspace.h"
#include "forcefieldfunctions.h"
#include "uff_par.h"
#include "src/core/units.h"
#include "src/core/curcuma_logger.h"

#include <fmt/core.h>
#include <fmt/format.h>

#include <cmath>

// ============================================================================
// UFF execution entry point
// Calls: bonds, angles, dihedrals, inversions, vdW, dispersion (if uff-d3)
// ============================================================================

void FFWorkspace::executeUFF(int p)
{
    calcUFFBonds(p);
    calcUFFAngles(p);
    calcUFFDihedrals(p);
    calcUFFInversions(p);
    calcUFFvdW(p);

    // D3/D4 dispersion (for uff-d3: dispersions vector is populated)
    if (m_dispersion_enabled) {
        if (!m_dispersions.empty())
            calcDispersion(p);
        if (!m_d4_dispersions.empty())
            calcD4Dispersion(p);
    }
}

// ============================================================================
// QMDFF execution entry point
// Calls: QMDFF bonds+angles, then UFF dihedrals/inversions/vdW (same formulas)
// ============================================================================

void FFWorkspace::executeQMDFF(int p)
{
    calcQMDFFBonds(p);
    calcQMDFFAngles(p);
    calcUFFDihedrals(p);   // QMDFF reuses UFF torsion formula
    calcUFFInversions(p);  // QMDFF reuses UFF inversion formula
    calcUFFvdW(p);         // QMDFF reuses UFF LJ formula

    if (m_dispersion_enabled) {
        if (!m_dispersions.empty())
            calcDispersion(p);
        if (!m_d4_dispersions.empty())
            calcD4Dispersion(p);
    }
}

// ============================================================================
// UFF Bond Stretching  — E = 0.5 * fc * (r - r0)²
// Reference: ForceFieldThread::CalculateUFFBondContribution()
// ============================================================================

void FFWorkspace::calcUFFBonds(int p)
{
    auto& acc = m_accumulators[p];
    const auto [beg, end] = m_partitions[p].bonds;

    for (int j = beg; j < end; ++j) {
        const auto& bond = m_bonds[j];
        Matrix derivate;
        double rij = UFF::BondStretching(m_geometry.row(bond.i),
                                          m_geometry.row(bond.j),
                                          derivate, m_do_gradient);

        acc.energy.bond += 0.5 * bond.fc * (rij - bond.r0_ij) * (rij - bond.r0_ij);

        if (m_do_gradient) {
            double diff = bond.fc * (rij - bond.r0_ij);
            acc.gradient.row(bond.i) += diff * derivate.row(0);
            acc.gradient.row(bond.j) += diff * derivate.row(1);
        }
    }
}

// ============================================================================
// UFF Angle Bending  — E = fc * (C0 + C1*cos(θ) + C2*cos(2θ))
// Reference: ForceFieldThread::CalculateUFFAngleContribution()
// ============================================================================

void FFWorkspace::calcUFFAngles(int p)
{
    auto& acc = m_accumulators[p];
    const auto [beg, end] = m_partitions[p].angles;

    for (int j = beg; j < end; ++j) {
        const auto& angle = m_angles[j];
        Matrix derivate;
        double costheta = UFF::AngleBending(m_geometry.row(angle.i),
                                             m_geometry.row(angle.j),
                                             m_geometry.row(angle.k),
                                             derivate, m_do_gradient);

        double cos2theta = 2.0 * costheta * costheta - 1.0;
        acc.energy.angle += angle.fc * (angle.C0 + angle.C1 * costheta + angle.C2 * cos2theta);

        if (m_do_gradient) {
            double diff = angle.fc * (angle.C1 + 4.0 * angle.C2 * costheta);
            acc.gradient.row(angle.i) += diff * derivate.row(0);
            acc.gradient.row(angle.j) += diff * derivate.row(1);
            acc.gradient.row(angle.k) += diff * derivate.row(2);
        }
    }
}

// ============================================================================
// UFF Dihedral Torsion — E = 0.5*V*(1 - cos(n*phi0)*cos(n*phi))
// Reference: ForceFieldThread::CalculateUFFDihedralContribution()
// ============================================================================

void FFWorkspace::calcUFFDihedrals(int p)
{
    auto& acc = m_accumulators[p];
    const auto [beg, end] = m_partitions[p].dihedrals;

    for (int j = beg; j < end; ++j) {
        const auto& dihedral = m_dihedrals[j];
        Eigen::Vector3d i = m_geometry.row(dihedral.i);
        Eigen::Vector3d jj = m_geometry.row(dihedral.j);
        Eigen::Vector3d k = m_geometry.row(dihedral.k);
        Eigen::Vector3d l = m_geometry.row(dihedral.l);

        Eigen::Vector3d nijk = UFF::NormalVector(i, jj, k);
        Eigen::Vector3d njkl = UFF::NormalVector(jj, k, l);
        double n_ijk = nijk.norm();
        double n_jkl = njkl.norm();
        double dotpr = nijk.dot(njkl);
        Eigen::Vector3d ji = jj - i;

        double sign = ((-ji).dot(njkl) < 0) ? -1.0 : 1.0;
        double phi = M_PI + sign * std::acos(dotpr / (n_ijk * n_jkl));
        double tmp_energy = 0.5 * dihedral.V * (1.0 - std::cos(dihedral.n * dihedral.phi0) * std::cos(dihedral.n * phi));
        if (std::isnan(tmp_energy))
            continue;

        acc.energy.dihedral += tmp_energy;

        if (m_do_gradient) {
            double dEdphi = 0.5 * dihedral.V * dihedral.n
                            * std::cos(dihedral.n * dihedral.phi0)
                            * std::sin(dihedral.n * phi);
            if (std::isnan(dEdphi))
                continue;

            Eigen::Vector3d kj = k - jj;
            Eigen::Vector3d kl = k - l;

            Eigen::Vector3d dEdi = dEdphi * kj.norm() / (nijk.squaredNorm()) * nijk;
            Eigen::Vector3d dEdl = -dEdphi * kj.norm() / (njkl.squaredNorm()) * njkl;
            Eigen::Vector3d dEdj = -dEdi + ((-ji).dot(kj) / kj.squaredNorm() * dEdi)
                                          - (kl.dot(kj)  / kj.squaredNorm() * dEdl);
            Eigen::Vector3d dEdk = -(dEdi + dEdj + dEdl);

            if (std::isnan(dEdi.sum()) || std::isnan(dEdj.sum()) ||
                std::isnan(dEdk.sum()) || std::isnan(dEdl.sum()))
                continue;

            acc.gradient.row(dihedral.i) += dEdi.transpose();
            acc.gradient.row(dihedral.j) += dEdj.transpose();
            acc.gradient.row(dihedral.k) += dEdk.transpose();
            acc.gradient.row(dihedral.l) += dEdl.transpose();
        }
    }
}

// ============================================================================
// UFF Out-of-Plane Inversion — E = fc * (C0 + C1*sinY + C2*cos2Y)
// Reference: ForceFieldThread::CalculateUFFInversionContribution()
// ============================================================================

void FFWorkspace::calcUFFInversions(int p)
{
    auto& acc = m_accumulators[p];
    const auto [beg, end] = m_partitions[p].inversions;

    for (int j = beg; j < end; ++j) {
        const auto& inv = m_inversions[j];
        Eigen::Vector3d i = m_geometry.row(inv.i);
        Eigen::Vector3d jj = m_geometry.row(inv.j);
        Eigen::Vector3d k = m_geometry.row(inv.k);
        Eigen::Vector3d l = m_geometry.row(inv.l);

        Eigen::Vector3d ail = SubVector(i, l);
        Eigen::Vector3d nijk = UFF::NormalVector(i, jj, k);

        double cosY = nijk.dot(ail) / (nijk.norm() * ail.norm());
        double sinYSq = 1.0 - cosY * cosY;
        double sinY = (sinYSq > 0.0) ? std::sqrt(sinYSq) : 0.0;
        double cos2Y = sinY * sinY - 1.0;

        double tmp_energy = inv.fc * (inv.C0 + inv.C1 * sinY + inv.C2 * cos2Y);
        if (std::isnan(tmp_energy))
            continue;
        acc.energy.inversion += tmp_energy;

        if (m_do_gradient) {
            Eigen::Vector3d ji = jj - i;
            Eigen::Vector3d jk = k - i;
            Eigen::Vector3d jl = l - i;

            if (ji.norm() < 1e-5 || jk.norm() < 1e-5 || jl.norm() < 1e-5)
                continue;

            double dji = ji.norm(), djk = jk.norm(), djl = jl.norm();
            ji /= dji; jk /= djk; jl /= djl;

            Eigen::Vector3d nijk2 = ji.cross(jk);
            nijk2 /= nijk2.norm();

            double cosY2 = nijk2.dot(jl);
            double sinYSq2 = 1.0 - cosY2 * cosY2;
            double sinY2 = (sinYSq2 > 0.0) ? std::sqrt(sinYSq2) : 0.0;
            double cosTheta = ji.dot(jk);
            double sinThetaSq = std::max(1.0 - cosTheta * cosTheta, 1.0e-8);
            double sinTheta = std::sqrt(sinThetaSq);

            double dEdY = -inv.fc * (inv.C1 * cosY2 - 4.0 * inv.C2 * cosY2 * sinY2);

            Eigen::Vector3d p1 = ji.cross(jk);
            Eigen::Vector3d p2 = jk.cross(jl);
            Eigen::Vector3d p3 = jl.cross(ji);

            double sin_dl = p1.dot(jl) / sinTheta;
            double dll = std::asin(sin_dl);
            (void)dll;  // used only for cos_dl concept; kept for clarity

            Eigen::Vector3d dYdl = (p1 / sinTheta - jl * sin_dl) / djl;
            Eigen::Vector3d dYdi = ((p2 + ((-ji + jk * cosTheta) * sin_dl) / sinTheta) / dji) / sinTheta;
            Eigen::Vector3d dYdk = ((p3 + ((-jk + ji * cosTheta) * sin_dl) / sinTheta) / djk) / sinTheta;
            Eigen::Vector3d dYdj = -(dYdi + dYdk + dYdl);

            acc.gradient.row(inv.i) += (dEdY * dYdj).transpose();
            acc.gradient.row(inv.j) += (dEdY * dYdi).transpose();
            acc.gradient.row(inv.k) += (dEdY * dYdk).transpose();
            acc.gradient.row(inv.l) += (dEdY * dYdl).transpose();
        }
    }
}

// ============================================================================
// UFF/QMDFF van der Waals (12-6 LJ) — E = C_ij * [-(r0/r)^6 + (r0/r)^12]/100
// Reference: ForceFieldThread::CalculateUFFvdWContribution()
// m_au converts Angström geometry → Bohr for LJ distance
// ============================================================================

void FFWorkspace::calcUFFvdW(int p)
{
    auto& acc = m_accumulators[p];
    const auto [beg, end] = m_partitions[p].vdws;

    for (int j = beg; j < end; ++j) {
        const auto& vdw = m_vdws[j];
        double ij = (m_geometry.row(vdw.i) - m_geometry.row(vdw.j)).norm() * m_au;
        double pow6 = std::pow(vdw.r0_ij / ij, 6);

        acc.energy.vdw += vdw.C_ij * (-2.0 * pow6) / 100.0;
        acc.energy.rep += vdw.C_ij * (pow6 * pow6)  / 100.0;

        if (m_do_gradient) {
            // dE/dr = 12*C_ij*(pow6 - pow6^2) / r^2 / 100
            // But we need to convert back: r_geom = ij/m_au
            double diff = 12.0 * vdw.C_ij * (pow6 - pow6 * pow6) / (ij * ij) / 100.0 * m_au;
            Eigen::Vector3d dr = (m_geometry.row(vdw.i) - m_geometry.row(vdw.j)).transpose();
            acc.gradient.row(vdw.i) += (diff * dr).transpose();
            acc.gradient.row(vdw.j) -= (diff * dr).transpose();
        }
    }
}

// ============================================================================
// QMDFF Bond Stretching — E = fc * (1 + (r0/r)^n - 2*(r0/r)^(0.75n))
// Reference: ForceFieldThread::CalculateQMDFFBondContribution()
// ============================================================================

void FFWorkspace::calcQMDFFBonds(int p)
{
    auto& acc = m_accumulators[p];
    const auto [beg, end] = m_partitions[p].bonds;

    for (int j = beg; j < end; ++j) {
        const auto& bond = m_bonds[j];
        Eigen::Vector3d vi = m_geometry.row(bond.i);
        Eigen::Vector3d vj = m_geometry.row(bond.j);
        Eigen::Vector3d ij = vi - vj;
        double distance = ij.norm();
        double fc = bond.fc;              // Force constant from fc field (Bond struct)
        double ratio = bond.r0_ij / distance;  // Equilibrium distance from r0_ij

        acc.energy.bond += fc * (1.0 + std::pow(ratio, bond.exponent)
                                     - 2.0 * std::pow(ratio, bond.exponent * 0.75));

        if (m_do_gradient) {
            double diff = fc * (-bond.exponent * std::pow(ratio, bond.exponent - 1.0)
                                + 2.0 * bond.exponent * 0.75 * std::pow(ratio, bond.exponent * 0.75 - 1.0));
            acc.gradient.row(bond.i) += (diff * ij / distance).transpose();
            acc.gradient.row(bond.j) -= (diff * ij / distance).transpose();
        }
    }
}

// ============================================================================
// QMDFF Angle Bending — E = fc * (cos(θ) - cos(θ0))²
// Reference: ForceFieldThread::CalculateQMDFFAngleContribution()
// ============================================================================

void FFWorkspace::calcQMDFFAngles(int p)
{
    auto& acc = m_accumulators[p];
    const auto [beg, end] = m_partitions[p].angles;

    for (int j = beg; j < end; ++j) {
        const auto& angle = m_angles[j];
        Matrix derivate;
        double costheta0 = std::cos(angle.theta0_ijk * M_PI / 180.0);
        double costheta = UFF::AngleBending(m_geometry.row(angle.i),
                                             m_geometry.row(angle.j),
                                             m_geometry.row(angle.k),
                                             derivate, true);

        double dcos = costheta - costheta0;
        acc.energy.angle += angle.fc * dcos * dcos;

        if (m_do_gradient) {
            double dEdTheta = 2.0 * angle.fc * dcos;
            acc.gradient.row(angle.i) -= (dEdTheta * derivate.row(0)).transpose();
            acc.gradient.row(angle.j) -= (dEdTheta * derivate.row(1)).transpose();
            acc.gradient.row(angle.k) -= (dEdTheta * derivate.row(2)).transpose();
        }
    }
}
