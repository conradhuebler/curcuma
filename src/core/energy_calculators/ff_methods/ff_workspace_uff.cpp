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
#include "qmdff_terms.h"
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
//
// Claude Generated (Aug 2026): the full QMDFF energy expression, mirroring the
// three Fortran drivers of xtb's removed `src/qmdff.f90`:
//   ff_eg   -> bonds, angles, torsions, out-of-plane   (qmdff.f90:285-509)
//   ff_nonb -> dispersion + electrostatics + repulsion (qmdff.f90:514-590)
//   ff_hb   -> hydrogen and halogen bonds              (qmdff.f90:596-626)
//
// The UFF torsion/inversion/LJ fall-backs are only used when a parameter set
// carries no genuine QMDFF torsion or NCI lists — i.e. for legacy parameter
// files written before the port. That keeps old `qmdff_param.json` files usable.
// ============================================================================

void FFWorkspace::executeQMDFF(int p)
{
    calcQMDFFBonds(p);
    calcQMDFFAngles(p);

    if (!m_qmdff_torsions.empty()) {
        // Proper torsions AND out-of-plane terms — QMDFF keeps both in one list
        // (tors(6,·) == 2 selects the inversion branch), so no UFF improper term
        // is evaluated on this path.
        calcQMDFFTorsions(p);
    } else {
        // Legacy fall-back for parameter sets written before the port
        calcUFFDihedrals(p);
        calcUFFInversions(p);
    }

    if (!m_qmdff_ncis.empty())
        calcQMDFFNonCovalent(p);
    else
        calcUFFvdW(p);          // legacy fall-back

    if (m_hbond_enabled && !m_qmdff_hbonds.empty())
        calcQMDFFHBonds(p);

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
// QMDFF term calculators
// ----------------------------------------------------------------------------
// Claude Generated (Aug 2026). Reference: xtb b7dbd36^:src/qmdff.f90.
//
// UNITS: the UFF/QMDFF geometry stored in m_geometry is in Angstrom, while the
// QMDFF energy expression is defined in atomic units. Every kernel therefore
// scales positions by m_au (Angstrom -> Bohr) on the way in and the resulting
// gradient by m_au on the way out (dE/dx_Angstrom = dE/dx_Bohr * m_au).
// Bond/angle/torsion energies are unit-free in the distance (they depend on
// r0/r or on angles only), but their damping factors are not — hence the
// conversion is applied uniformly rather than term by term.
// ============================================================================

namespace {

/// Element number of atom `i`, or 1 (hydrogen) if the type table is short.
inline int atomZ(const std::vector<int>& types, int i)
{
    return (i >= 0 && i < static_cast<int>(types.size())) ? types[i] : 1;
}

} // namespace

// ============================================================================
// QMDFF Bond Stretching (qmdff.f90:316-354)
//   LJ    : E = k [1 + (r0/r)^a - 2 (r0/r)^(a/2)]
//   Morse : E = k [1 - exp(-a/(2 r0) (r - r0))]^2
// The potential type per bond comes from Bond::qmdff_potential (Fortran bond(3,·)).
// ============================================================================

void FFWorkspace::calcQMDFFBonds(int p)
{
    auto& acc = m_accumulators[p];
    const auto [beg, end] = m_partitions[p].bonds;

    for (int b = beg; b < end; ++b) {
        const auto& bond = m_bonds[b];
        const Eigen::Vector3d ri = m_geometry.row(bond.i).transpose() * m_au;
        const Eigen::Vector3d rj = m_geometry.row(bond.j).transpose() * m_au;

        Eigen::Vector3d gi = Eigen::Vector3d::Zero();
        Eigen::Vector3d gj = Eigen::Vector3d::Zero();

        const auto potential = (bond.qmdff_potential == 1)
            ? QMDFFTerms::BondPotential::Morse
            : QMDFFTerms::BondPotential::LennardJones;

        acc.energy.bond += QMDFFTerms::bondEnergy(ri, rj,
                                                  bond.r0_ij * m_au, bond.fc, bond.exponent,
                                                  potential, gi, gj, m_do_gradient);

        if (m_do_gradient) {
            acc.gradient.row(bond.i) += (gi * m_au).transpose();
            acc.gradient.row(bond.j) += (gj * m_au).transpose();
        }
    }
}

// ============================================================================
// QMDFF Angle Bending (qmdff.f90:357-402)
//   E = k (cos(theta) - cos(theta0))^2 * damp(r_ij) * damp(r_kj)
//   (harmonic in theta instead of cos for linear reference angles)
//
// Angle::j is the central atom in curcuma's convention, which matches the
// Fortran angl(1,·) entry.
// ============================================================================

void FFWorkspace::calcQMDFFAngles(int p)
{
    auto& acc = m_accumulators[p];
    const auto [beg, end] = m_partitions[p].angles;

    for (int a = beg; a < end; ++a) {
        const auto& angle = m_angles[a];
        const Eigen::Vector3d ra = m_geometry.row(angle.i).transpose() * m_au;
        const Eigen::Vector3d rb = m_geometry.row(angle.j).transpose() * m_au; // centre
        const Eigen::Vector3d rc = m_geometry.row(angle.k).transpose() * m_au;

        Eigen::Vector3d ga = Eigen::Vector3d::Zero();
        Eigen::Vector3d gb = Eigen::Vector3d::Zero();
        Eigen::Vector3d gc = Eigen::Vector3d::Zero();

        acc.energy.angle += QMDFFTerms::angleEnergy(
            ra, rb, rc,
            atomZ(m_atom_types, angle.i), atomZ(m_atom_types, angle.j), atomZ(m_atom_types, angle.k),
            angle.theta0_ijk, angle.fc,
            ga, gb, gc, m_do_gradient);

        if (m_do_gradient) {
            acc.gradient.row(angle.i) += (ga * m_au).transpose();
            acc.gradient.row(angle.j) += (gb * m_au).transpose();
            acc.gradient.row(angle.k) += (gc * m_au).transpose();
        }
    }
}

// ============================================================================
// QMDFF torsions and out-of-plane terms (qmdff.f90:405-507)
// Both live in the same Fortran list; tors(6,·) == 2 selects the inversion
// branch, which is why they share one loop here as well.
// ============================================================================

void FFWorkspace::calcQMDFFTorsions(int p)
{
    auto& acc = m_accumulators[p];
    const auto [beg, end] = m_partitions[p].qmdff_torsions;

    for (int t = beg; t < end; ++t) {
        const auto& tor = m_qmdff_torsions[t];
        const Eigen::Vector3d ri = m_geometry.row(tor.i).transpose() * m_au;
        const Eigen::Vector3d rj = m_geometry.row(tor.j).transpose() * m_au;
        const Eigen::Vector3d rk = m_geometry.row(tor.k).transpose() * m_au;
        const Eigen::Vector3d rl = m_geometry.row(tor.l).transpose() * m_au;

        const int zi = atomZ(m_atom_types, tor.i);
        const int zj = atomZ(m_atom_types, tor.j);
        const int zk = atomZ(m_atom_types, tor.k);
        const int zl = atomZ(m_atom_types, tor.l);

        Eigen::Matrix<double, 4, 3> g = Eigen::Matrix<double, 4, 3>::Zero();

        if (tor.out_of_plane) {
            acc.energy.inversion += QMDFFTerms::outOfPlaneEnergy(
                ri, rj, rk, rl, zi, zj, zk, zl, tor, g, m_do_gradient);
        } else {
            acc.energy.dihedral += QMDFFTerms::torsionEnergy(
                ri, rj, rk, rl, zi, zj, zk, zl, tor, g, m_do_gradient);
        }

        if (m_do_gradient) {
            acc.gradient.row(tor.i) += g.row(0) * m_au;
            acc.gradient.row(tor.j) += g.row(1) * m_au;
            acc.gradient.row(tor.k) += g.row(2) * m_au;
            acc.gradient.row(tor.l) += g.row(3) * m_au;
        }
    }
}

// ============================================================================
// QMDFF non-covalent pairs (qmdff.f90:514-590, ff_nonb)
// D3(BJ) dispersion + point-charge electrostatics + Born-Mayer repulsion,
// each scaled by the topological class of the pair (eps1/eps2).
// ============================================================================

void FFWorkspace::calcQMDFFNonCovalent(int p)
{
    auto& acc = m_accumulators[p];
    const auto [beg, end] = m_partitions[p].qmdff_ncis;

    for (int n = beg; n < end; ++n) {
        const auto& nci = m_qmdff_ncis[n];
        const Eigen::Vector3d ri = m_geometry.row(nci.i).transpose() * m_au;
        const Eigen::Vector3d rj = m_geometry.row(nci.j).transpose() * m_au;

        Eigen::Vector3d gi = Eigen::Vector3d::Zero();
        Eigen::Vector3d gj = Eigen::Vector3d::Zero();
        double e_disp = 0.0, e_es = 0.0, e_rep = 0.0;

        QMDFFTerms::nonCovalentPair(ri, rj, nci.c6, nci.r0_bj, nci.sr42,
                                    nci.zab, nci.alpha, nci.qq, nci.nk,
                                    e_disp, e_es, e_rep, gi, gj, m_do_gradient);

        acc.energy.dispersion += e_disp;
        acc.energy.coulomb += e_es;
        acc.energy.nonbonded_rep += e_rep;

        if (m_do_gradient) {
            acc.gradient.row(nci.i) += (gi * m_au).transpose();
            acc.gradient.row(nci.j) += (gj * m_au).transpose();
        }
    }
}

// ============================================================================
// QMDFF hydrogen and halogen bonds (qmdff.f90:596-626, ff_hb)
// ============================================================================

void FFWorkspace::calcQMDFFHBonds(int p)
{
    auto& acc = m_accumulators[p];
    const auto [beg, end] = m_partitions[p].qmdff_hbonds;

    for (int h = beg; h < end; ++h) {
        const auto& hb = m_qmdff_hbonds[h];
        const Eigen::Vector3d ra = m_geometry.row(hb.a).transpose() * m_au;
        const Eigen::Vector3d rb = m_geometry.row(hb.b).transpose() * m_au;
        const Eigen::Vector3d rh = m_geometry.row(hb.h).transpose() * m_au;

        // The reference skips widely separated A...B pairs (qmdff.f90:614)
        if ((ra - rb).norm() > 25.0)
            continue;

        Eigen::Vector3d ga = Eigen::Vector3d::Zero();
        Eigen::Vector3d gb = Eigen::Vector3d::Zero();
        Eigen::Vector3d gh = Eigen::Vector3d::Zero();

        if (hb.halogen) {
            acc.energy.xbond += QMDFFTerms::halogenBondEnergyGradient(
                ra, rb, rh, hb.c1, ga, gb, gh, m_do_gradient);
        } else {
            acc.energy.hbond += QMDFFTerms::hydrogenBondEnergy(
                ra, rb, rh, hb.c1, hb.c2, ga, gb, gh, m_do_gradient);
        }

        if (m_do_gradient) {
            acc.gradient.row(hb.a) += (ga * m_au).transpose();
            acc.gradient.row(hb.b) += (gb * m_au).transpose();
            acc.gradient.row(hb.h) += (gh * m_au).transpose();
        }
    }
}
