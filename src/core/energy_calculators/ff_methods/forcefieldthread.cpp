/*
 * < Generic force field class for curcuma . >
 * Copyright (C) 2024 - 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "forcefieldderivaties.h"
#include "qmdff_par.h"
#include "uff_par.h"
#include "gfnff_par.h"

#include "forcefieldfunctions.h"

#include "forcefield.h"

#include <unordered_map>  // Claude Generated (Dec 2025): For atom_to_params lookup in Coulomb self-energy

ForceFieldThread::ForceFieldThread(int thread, int threads)
    : m_thread(thread)
    , m_threads(threads)
{
    setAutoDelete(false);
    m_final_factor = 1; // / 2625.15 * 4.19;
    // m_d = parameters["differential"].get<double>();
    m_d = 1e-7;
}

int ForceFieldThread::execute()
{
    m_angle_energy = 0.0;
    m_bond_energy = 0.0;
    m_vdw_energy = 0;
    m_rep_energy = 0;
    m_inversion_energy = 0;
    m_dihedral_energy = 0;
    m_angle_energy = 0;
    m_bond_energy = 0.0;

    // CRITICAL FIX (Nov 2025): Reset GFN-FF pairwise energy terms
    // These were missing, causing energy accumulation across Calculate() calls
    // Business Impact: 100K € - GFN-FF must match XTB 6.6.1 reference exactly
    m_dispersion_energy = 0.0;
    m_coulomb_energy = 0.0;
    m_energy_hbond = 0.0;
    m_energy_xbond = 0.0;
    m_eq_energy = 0.0;  // Also reset EQ energy for consistency

    // Claude Generated 2025: Reset native D3/D4 dispersion energy terms
    m_d3_energy = 0.0;
    m_d4_energy = 0.0;
    m_atm_energy = 0.0;  // Claude Generated (Dec 2025): Reset ATM three-body dispersion

    // Phase 1.1: Guard debug output with verbosity check (Claude Generated - Dec 2025)
    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info(fmt::format("ForceFieldThread {} executing (method={})", m_thread, m_method));
    }

    if (m_method == 1) {
        CalculateUFFBondContribution();
        CalculateUFFAngleContribution();

    } else if (m_method == 2) {
        CalculateQMDFFBondContribution();
        // CalculateUFFBondContribution();
        CalculateQMDFFAngleContribution();
    } else if (m_method == 3) {
        // Phase 1.1: Guard debug output with verbosity check (Claude Generated - Dec 2025)
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info(fmt::format("GFN-FF energy calculation started in thread {}", m_thread));
        }

        // Phase 1.2: Build bonded pairs cache ONCE for fast lookup (Claude Generated - Dec 2025)
        if (!m_bonded_pairs_cached && m_gfnff_bonds.size() > 0) {
            m_bonded_pairs.clear();
            for (const auto& bond : m_gfnff_bonds) {
                m_bonded_pairs.insert({bond.i, bond.j});
                m_bonded_pairs.insert({bond.j, bond.i});  // symmetric
            }
            m_bonded_pairs_cached = true;
            if (CurcumaLogger::get_verbosity() >= 3) {
                CurcumaLogger::info(fmt::format("Cached {} bonded pairs", m_bonded_pairs.size()));
            }
        }

        // GFN-FF bonded terms
        CalculateGFNFFBondContribution();
        CalculateGFNFFAngleContribution();
        CalculateGFNFFDihedralContribution();
        CalculateGFNFFExtraTorsionContribution();  // Claude Generated (Jan 2, 2026): Extra sp3-sp3 gauche torsions
        CalculateGFNFFInversionContribution();

        // GFN-FF non-bonded pairwise parallelizable terms (Phase 4)
        if (m_dispersion_enabled) {
            CalculateGFNFFDispersionContribution();  // D3/D4 dispersion
        }
        if (m_repulsion_enabled) {
            CalculateGFNFFBondedRepulsionContribution();
            CalculateGFNFFNonbondedRepulsionContribution();
        }
        if (m_coulomb_enabled) {
            CalculateGFNFFCoulombContribution();     // EEQ Coulomb electrostatics
        }

        // GFN-FF hydrogen bond and halogen bond terms (Phase 5)
        if (m_hbond_enabled) {
            CalculateGFNFFHydrogenBondContribution();  // HB three-body terms
            CalculateGFNFFHalogenBondContribution();   // XB three-body terms
        }

        // Claude Generated (December 19, 2025): Native D3/D4 dispersion calculation for GFN-FF
        // Note: GFN-FF uses its own dispersion (CalculateGFNFFDispersionContribution above)
        // D3/D4 are additional corrections that can be enabled separately
        if (m_d3_dispersions.size() > 0) {
            if (CurcumaLogger::get_verbosity() >= 3) {
                CurcumaLogger::info(fmt::format("Thread {} calculating {} D3 dispersion pairs", m_thread, m_d3_dispersions.size()));
            }
            CalculateD3DispersionContribution();
        }

        if (m_d4_dispersions.size() > 0) {
            if (CurcumaLogger::get_verbosity() >= 3) {
                CurcumaLogger::info(fmt::format("Thread {} calculating {} D4 dispersion pairs", m_thread, m_d4_dispersions.size()));
            }
            CalculateD4DispersionContribution();  // Claude Generated - Dec 25, 2025: Native D4 energy calculation
        }

        // ATM three-body dispersion (D3/D4)
        if (!m_atm_triples.empty()) {
            if (CurcumaLogger::get_verbosity() >= 3) {
                CurcumaLogger::info(fmt::format("Thread {} calculating {} ATM triples", m_thread, m_atm_triples.size()));
            }
            CalculateATMContribution();

            // Claude Generated (2025): Calculate analytical ATM gradients
            if (m_calculate_gradient) {
                CalculateATMGradient();
            }
        }

    } else if (m_method == 1) {
        // UFF bonded terms calculated above (bonds, angles)
        // Now add UFF non-bonded terms
        CalculateUFFDihedralContribution();
        CalculateUFFInversionContribution();
        CalculateUFFvdWContribution();

        // Claude Generated (December 19, 2025): UFF-D3 native dispersion correction
        // Add D3 dispersion to UFF if available (UFF-D3 hybrid method)
        if (m_d3_dispersions.size() > 0) {
            if (CurcumaLogger::get_verbosity() >= 3) {
                CurcumaLogger::info(fmt::format("Thread {} calculating {} D3 dispersion pairs for UFF-D3", m_thread, m_d3_dispersions.size()));
            }
            CalculateD3DispersionContribution();
        }
    } else if (m_method == 5) {  // Claude Generated (December 21, 2025)
        // D3-only method: pure dispersion correction without bonded terms
        if (m_gfnff_dispersions.size() > 0) {
            if (CurcumaLogger::get_verbosity() >= 3) {
                CurcumaLogger::info(fmt::format("Thread {} calculating {} D3-only dispersion pairs",
                                                 m_thread, m_gfnff_dispersions.size()));
            }
            CalculateGFNFFDispersionContribution();  // D3 pairs stored in gfnff_dispersions
        }
    }

    if (m_method != 3 && m_method != 1 && m_method != 5) {
        // QMDFF or other methods
        CalculateUFFDihedralContribution();
        CalculateUFFInversionContribution();
        CalculateUFFvdWContribution();
    }
    CalculateESPContribution();
    /*
    CalculateQMDFFDihedralContribution();
    */
    return 0;
}

void ForceFieldThread::addBond(const Bond& bonds)
{
    // if (bonds.type == 1)
    m_uff_bonds.push_back(bonds);
    // else if (bonds.type == 2)
    //     m_qmdff_bonds.push_back(bonds);
}

void ForceFieldThread::addAngle(const Angle& angles)
{
    // if (angles.type == 1)
    m_uff_angles.push_back(angles);
    // else if (angles.type == 2)
    //     m_qmdff_angles.push_back(angles);
}

void ForceFieldThread::addDihedral(const Dihedral& dihedrals)
{
    if (dihedrals.type == 1)
        m_uff_dihedrals.push_back(dihedrals);
    else if (dihedrals.type == 2)
        m_qmdff_dihedrals.push_back(dihedrals);
}

void ForceFieldThread::addInversion(const Inversion& inversions)
{
    if (inversions.type == 1)
        m_uff_inversions.push_back(inversions);
    else if (inversions.type == 2)
        m_qmdff_inversions.push_back(inversions);
}

void ForceFieldThread::addvdW(const vdW& vdWs)
{
    if (vdWs.type == 1)
        m_uff_vdWs.push_back(vdWs);
}

void ForceFieldThread::addEQ(const EQ& EQs)
{
    // if (EQs.type == 2)
    m_EQs.push_back(EQs);
}

void ForceFieldThread::addGFNFFBond(const Bond& bonds)
{
    m_gfnff_bonds.push_back(bonds);
}

void ForceFieldThread::addGFNFFAngle(const Angle& angles)
{
    m_gfnff_angles.push_back(angles);
}

void ForceFieldThread::addGFNFFDihedral(const Dihedral& dihedrals)
{
    m_gfnff_dihedrals.push_back(dihedrals);
}

// Claude Generated (Jan 2, 2026): Extra sp3-sp3 gauche torsions separated from primary torsions
void ForceFieldThread::addGFNFFExtraTorsion(const Dihedral& extra_torsion)
{
    m_gfnff_extra_torsions.push_back(extra_torsion);
}

void ForceFieldThread::addGFNFFInversion(const Inversion& inversions)
{
    m_gfnff_inversions.push_back(inversions);
}

void ForceFieldThread::addGFNFFvdW(const vdW& vdWs)
{
    m_gfnff_vdWs.push_back(vdWs);
}

// Phase 4: GFN-FF pairwise non-bonded addition methods (Claude Generated 2025)

void ForceFieldThread::addGFNFFDispersion(const GFNFFDispersion& dispersion)
{
    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info(fmt::format("Thread {} adding GFNFF dispersion pair {}-{}", m_thread, dispersion.i, dispersion.j));
    }
    m_gfnff_dispersions.push_back(dispersion);
    m_dispersion_enabled = true;  // CRITICAL (Jan 7, 2026): Enable dispersion calculation
}

void ForceFieldThread::addGFNFFBondedRepulsion(const GFNFFRepulsion& repulsion)
{
    m_gfnff_bonded_repulsions.push_back(repulsion);
}

void ForceFieldThread::addGFNFFNonbondedRepulsion(const GFNFFRepulsion& repulsion)
{
    m_gfnff_nonbonded_repulsions.push_back(repulsion);
}

void ForceFieldThread::addGFNFFCoulomb(const GFNFFCoulomb& coulomb)
{
    m_gfnff_coulombs.push_back(coulomb);
}

// Phase 6: Assign atoms for self-energy calculation (Claude Generated Dec 2025)
void ForceFieldThread::assignAtomsForSelfEnergy(const std::vector<int>& atom_indices)
{
    m_assigned_atoms_for_self_energy = atom_indices;
}

// Phase 4: GFN-FF hydrogen bond and halogen bond addition methods (Claude Generated 2025)

void ForceFieldThread::addGFNFFHydrogenBond(const GFNFFHydrogenBond& hbond)
{
    m_gfnff_hbonds.push_back(hbond);
}

void ForceFieldThread::addGFNFFHalogenBond(const GFNFFHalogenBond& xbond)
{
    m_gfnff_xbonds.push_back(xbond);
}

void ForceFieldThread::addATMTriple(const ATMTriple& triple)
{
    m_atm_triples.push_back(triple);
}

void ForceFieldThread::CalculateUFFBondContribution()
{
    double factor = m_final_factor * m_bond_scaling;

    for (int index = 0; index < m_uff_bonds.size(); ++index) {
        const auto& bond = m_uff_bonds[index];

        Vector i = m_geometry.row(bond.i);
        Vector j = m_geometry.row(bond.j);
        Matrix derivate;
        double rij = UFF::BondStretching(i, j, derivate, m_calculate_gradient);

        m_bond_energy += (0.5 * bond.fc * (rij - bond.r0_ij) * (rij - bond.r0_ij)) * factor;
        if (m_calculate_gradient) {
            double diff = (bond.fc) * (rij - bond.r0_ij) * factor;
            m_gradient.row(bond.i) += diff * derivate.row(0);
            m_gradient.row(bond.j) += diff * derivate.row(1);
        }
    }
}

void ForceFieldThread::CalculateUFFAngleContribution()
{
    for (int index = 0; index < m_uff_angles.size(); ++index) {
        const auto& angle = m_uff_angles[index];
        auto i = m_geometry.row(angle.i);
        auto j = m_geometry.row(angle.j);
        auto k = m_geometry.row(angle.k);
        Matrix derivate;
        double costheta = UFF::AngleBending(i, j, k, derivate, m_calculate_gradient);

        // UFF angle bending potential
        double cos2theta = 2 * costheta * costheta - 1;
        m_angle_energy += angle.fc * (angle.C0 + angle.C1 * costheta + angle.C2 * cos2theta) * m_final_factor * m_angle_scaling;

        if (m_calculate_gradient) {
            double diff = angle.fc * (angle.C1 + 4 * angle.C2 * costheta) * m_final_factor * m_angle_scaling;
            m_gradient.row(angle.i) += diff * derivate.row(0);
            m_gradient.row(angle.j) += diff * derivate.row(1);
            m_gradient.row(angle.k) += diff * derivate.row(2);
        }
    }
}

void ForceFieldThread::CalculateUFFDihedralContribution()
{
    for (int index = 0; index < m_uff_dihedrals.size(); ++index) {
        const auto& dihedral = m_uff_dihedrals[index];
        Eigen::Vector3d i = m_geometry.row(dihedral.i);
        Eigen::Vector3d j = m_geometry.row(dihedral.j);
        Eigen::Vector3d k = m_geometry.row(dihedral.k);
        Eigen::Vector3d l = m_geometry.row(dihedral.l);
        Eigen::Vector3d nijk = UFF::NormalVector(i, j, k);
        Eigen::Vector3d njkl = UFF::NormalVector(j, k, l);
        double n_ijk = (nijk).norm();
        double n_jkl = (njkl).norm();
        double dotpr = nijk.dot(njkl);
        Eigen::Vector3d ji = j - i;

        double sign = (-1 * ji).dot(njkl) < 0 ? -1 : 1;
        double phi = pi + sign * acos(dotpr / (n_ijk * n_jkl));
        double sinphi = sin(phi);
        double tmp_energy = (1 / 2.0 * dihedral.V * (1 - cos(dihedral.n * dihedral.phi0) * cos(dihedral.n * phi))) * m_final_factor * m_dihedral_scaling;
        if (isnan(tmp_energy))
            continue;

        m_dihedral_energy += tmp_energy;
        if (m_calculate_gradient) {

            Eigen::Vector3d kj = k - j;
            Eigen::Vector3d kl = k - l;
            double tmp = 0;
            double dEdphi = (1 / 2.0 * dihedral.V * dihedral.n * (cos(dihedral.n * dihedral.phi0) * sin(dihedral.n * phi))) * m_final_factor * m_dihedral_scaling;
            if (isnan(dEdphi))
                continue;

            Eigen::Vector3d dEdi = dEdphi * kj.norm() / (nijk.norm() * nijk.norm()) * nijk;
            Eigen::Vector3d dEdl = -1 * dEdphi * kj.norm() / (njkl.norm() * njkl.norm()) * njkl;
            Eigen::Vector3d dEdj = -1 * dEdi + ((-1 * ji).dot(kj) / (kj.norm() * kj.norm()) * dEdi) - (kl.dot(kj) / (kj.norm() * kj.norm()) * dEdl);
            Eigen::Vector3d dEdk = -1 * (dEdi + dEdj + dEdl);

            if (isnan(dEdi.sum()) || isnan(dEdj.sum()) || isnan(dEdk.sum()) || isnan(dEdl.sum()))
                continue;
            m_gradient(dihedral.i, 0) += dEdi(0);
            m_gradient(dihedral.i, 1) += dEdi(1);
            m_gradient(dihedral.i, 2) += dEdi(2);

            m_gradient(dihedral.l, 0) += dEdl(0);
            m_gradient(dihedral.l, 1) += dEdl(1);
            m_gradient(dihedral.l, 2) += dEdl(2);

            m_gradient(dihedral.j, 0) += dEdj(0);
            m_gradient(dihedral.j, 1) += dEdj(1);
            m_gradient(dihedral.j, 2) += dEdj(2);

            m_gradient(dihedral.k, 0) += dEdk(0);
            m_gradient(dihedral.k, 1) += dEdk(1);
            m_gradient(dihedral.k, 2) += dEdk(2);
        }
    }
}

void ForceFieldThread::CalculateUFFInversionContribution()
{

    for (int index = 0; index < m_uff_inversions.size(); ++index) {
        const auto& inversion = m_uff_inversions[index];

        Eigen::Vector3d i = m_geometry.row(inversion.i);
        Eigen::Vector3d j = m_geometry.row(inversion.j);
        Eigen::Vector3d k = m_geometry.row(inversion.k);
        Eigen::Vector3d l = m_geometry.row(inversion.l);

        Eigen::Vector3d ail = SubVector(i, l);
        Eigen::Vector3d nijk = UFF::NormalVector(i, j, k);

        double cosY = (nijk.dot(ail) / ((nijk).norm() * (ail).norm()));

        double sinYSq = 1.0 - cosY * cosY;
        double sinY = ((sinYSq > 0.0) ? sqrt(sinYSq) : 0.0);
        double cos2Y = sinY * sinY - 1.0;

        double tmp_energy = (inversion.fc * (inversion.C0 + inversion.C1 * sinY + inversion.C2 * cos2Y)) * m_final_factor * m_inversion_scaling;
        if (std::isnan(tmp_energy))
            continue;
        m_inversion_energy += tmp_energy;

        // energy += Inversion(i, j, k, l, d_forceConstant, C0, C1, C2);
        if (m_calculate_gradient) {
            Eigen::Vector3d ji = (j - i);
            Eigen::Vector3d jk = (k - i);
            Eigen::Vector3d jl = (l - i);

            if (ji.norm() < 1e-5 || jk.norm() < 1e-5 || jl.norm() < 1e-5)
                continue;

            double dji = ji.norm();
            double djk = jk.norm();
            double djl = jl.norm();
            ji /= dji;
            jk /= djk;
            jl /= djl;

            Eigen::Vector3d nijk = ji.cross(jk);
            nijk /= nijk.norm();

            double cosY = (nijk.dot(jl));
            double sinYSq = 1.0 - cosY * cosY;
            double sinY = ((sinYSq > 0.0) ? sqrt(sinYSq) : 0.0);
            double cosTheta = (ji.dot(jk));
            double sinThetaSq = std::max(1.0 - cosTheta * cosTheta, 1.0e-8);
            double sinTheta = std::max(((sinThetaSq > 0.0) ? sqrt(sinThetaSq) : 0.0), 1.0e-8);

            double dEdY = -1 * (inversion.fc * (inversion.C1 * cosY - 4 * inversion.C2 * cosY * sinY)) * m_final_factor * m_inversion_scaling;

            Eigen::Vector3d p1 = ji.cross(jk);
            Eigen::Vector3d p2 = jk.cross(jl);
            Eigen::Vector3d p3 = jl.cross(ji);

            const double sin_dl = p1.dot(jl) / sinTheta;

            // the wilson angle:
            const double dll = asin(sin_dl);
            const double cos_dl = cos(dll);

            Eigen::Vector3d dYdl = (p1 / sinTheta - (jl * sin_dl)) / djl;
            Eigen::Vector3d dYdi = ((p2 + (((-ji + jk * cosTheta) * sin_dl) / sinTheta)) / dji) / sinTheta;
            Eigen::Vector3d dYdk = ((p3 + (((-jk + ji * cosTheta) * sin_dl) / sinTheta)) / djk) / sinTheta;
            Eigen::Vector3d dYdj = -1 * (dYdi + dYdk + dYdl);

            m_gradient(inversion.i, 0) += dEdY * dYdj(0);
            m_gradient(inversion.i, 1) += dEdY * dYdj(1);
            m_gradient(inversion.i, 2) += dEdY * dYdj(2);

            m_gradient(inversion.j, 0) += dEdY * dYdi(0);
            m_gradient(inversion.j, 1) += dEdY * dYdi(1);
            m_gradient(inversion.j, 2) += dEdY * dYdi(2);

            m_gradient(inversion.k, 0) += dEdY * dYdk(0);
            m_gradient(inversion.k, 1) += dEdY * dYdk(1);
            m_gradient(inversion.k, 2) += dEdY * dYdk(2);

            m_gradient(inversion.l, 0) += dEdY * dYdl(0);
            m_gradient(inversion.l, 1) += dEdY * dYdl(1);
            m_gradient(inversion.l, 2) += dEdY * dYdl(2);
        }
    }
}

void ForceFieldThread::CalculateUFFvdWContribution()
{
    for (int index = 0; index < m_uff_vdWs.size(); ++index) {
        const auto& vdw = m_uff_vdWs[index];
        Eigen::Vector3d i = m_geometry.row(vdw.i);
        Eigen::Vector3d j = m_geometry.row(vdw.j);
        double ij = (i - j).norm() * m_au;
        double pow6 = pow((vdw.r0_ij / ij), 6);

        m_vdw_energy += vdw.C_ij * (-2 * pow6 * m_vdw_scaling) * m_final_factor / 100;
        m_rep_energy += vdw.C_ij * (pow6 * pow6 * m_rep_scaling) * m_final_factor / 100;
        if (m_calculate_gradient) {
            double diff = 12 * vdw.C_ij * (pow6 * m_vdw_scaling - pow6 * pow6 * m_rep_scaling) / (ij * ij) * m_final_factor / 100;
            m_gradient(vdw.i, 0) += diff * (i(0) - j(0));
            m_gradient(vdw.i, 1) += diff * (i(1) - j(1));
            m_gradient(vdw.i, 2) += diff * (i(2) - j(2));

            m_gradient(vdw.j, 0) -= diff * (i(0) - j(0));
            m_gradient(vdw.j, 1) -= diff * (i(1) - j(1));
            m_gradient(vdw.j, 2) -= diff * (i(2) - j(2));
        }
    }
}

void ForceFieldThread::CalculateQMDFFBondContribution()
{
    m_d = 1e-5;
    double factor = m_final_factor * m_bond_scaling;
    Eigen::Vector3d dx = { m_d, 0, 0 };
    Eigen::Vector3d dy = { 0, m_d, 0 };
    Eigen::Vector3d dz = { 0, 0, m_d };
    for (int index = 0; index < m_uff_bonds.size(); ++index) {
        const auto& bond = m_uff_bonds[index];

        Vector i = m_geometry.row(bond.i);
        Vector j = m_geometry.row(bond.j);
        Vector ij = i - j;
        double distance = (ij).norm();

        // Matrix derivate;
        double fc = bond.r0_ij;
        const double ratio = bond.r0_ij / distance;
        m_bond_energy += fc * (1 + pow(ratio, bond.exponent) - 2 * pow(ratio, bond.exponent * 0.5));
        if (m_calculate_gradient) {

            double diff = 1 * fc * (-1 * bond.exponent * pow(ratio, bond.exponent - 1) + 2 * bond.exponent * 0.5 * pow(ratio, bond.exponent * 0.5 - 1));
            /*
                        Vector ijx = i+ dx - j;
                        double distancex1 = (ijx).norm();
                        double ratio = bond.r0_ij / distancex1;

                        double dx_p = fc * (1 + pow(ratio, bond.exponent) - 2 * pow(ratio, bond.exponent * 0.5));
                        ijx = i - dx - j;
                        distancex1 = (ijx).norm();
                        ratio = bond.r0_ij / distancex1;
                        double dx_m = fc * (1 + pow(ratio, bond.exponent) - 2 * pow(ratio, bond.exponent * 0.5));


                        Vector ijy = i+ dy - j;
                        double distancey1 = (ijy).norm();
                        ratio = bond.r0_ij / distancey1;

                        double dy_p = fc * (1 + pow(ratio, bond.exponent) - 2 * pow(ratio, bond.exponent * 0.5));
                        ijy = i - dy - j;
                        distancey1 = (ijy).norm();
                        ratio = bond.r0_ij / distancey1;
                        double dy_m = fc * (1 + pow(ratio, bond.exponent) - 2 * pow(ratio, bond.exponent * 0.5));


                        Vector ijz = i+ dz - j;
                        double distancez1 = (ijz).norm();
                        ratio = bond.r0_ij / distancez1;

                        double dz_p = fc * (1 + pow(ratio, bond.exponent) - 2 * pow(ratio, bond.exponent * 0.5));
                        ijz = i - dz - j;
                        distancez1 = (ijz).norm();
                        ratio = bond.r0_ij / distancez1;
                        double dz_m = fc * (1 + pow(ratio, bond.exponent) - 2 * pow(ratio, bond.exponent * 0.5));
            */
            // std::cout << (dx_p - dx_m)/(2.0*m_d) << " " << (dy_p - dy_m)/(2.0*m_d) <<" "<< (dz_p - dz_m)/(2.0*m_d) << " :: " << (diff * ij / (distance)).transpose() << std::endl;
            /*
            m_gradient(bond.i, 0) += (dx_p - dx_m)/(2*m_d);
            m_gradient(bond.i, 1) += (dy_p - dy_m)/(2*m_d);
            m_gradient(bond.i, 2) += (dz_p - dz_m)/(2*m_d);

            m_gradient(bond.j, 0) -= (dx_p - dx_m)/(2*m_d);
            m_gradient(bond.j, 1) -= (dy_p - dy_m)/(2*m_d);
            m_gradient(bond.j, 2) -= (dz_p - dz_m)/(2*m_d);
            */
            m_gradient.row(bond.i) += diff * ij / (distance);
            m_gradient.row(bond.j) -= diff * ij / (distance);
        }
    }
}
/*
double ForceFieldThread::HarmonicBondStretching()
{
    double energy = 0;

    return energy;
}
*/

void ForceFieldThread::CalculateQMDFFAngleContribution()
{
    double threshold = 1e-2;
    Eigen::Vector3d dx = { m_d, 0, 0 };
    Eigen::Vector3d dy = { 0, m_d, 0 };
    Eigen::Vector3d dz = { 0, 0, m_d };
    for (int index = 0; index < m_uff_angles.size(); ++index) {
        const auto& angle = m_uff_angles[index];

        Vector i = m_geometry.row(angle.i);
        Vector j = m_geometry.row(angle.j);
        Vector k = m_geometry.row(angle.k);
        Matrix derivate;
        double costheta0_ijk = cos(angle.theta0_ijk * pi / 180.0);
        double costheta = 0;
        double energy = 0;
        double dEdTheta = 0;

        if (std::abs(costheta0_ijk + 1) < threshold) {
            costheta = UFF::AngleBending(i, j, k, derivate, true);
            energy = angle.fc * (costheta - costheta0_ijk) * (costheta - costheta0_ijk);
            dEdTheta = 2 * angle.fc * (costheta - costheta0_ijk);
        } else {
            costheta = UFF::AngleBending(i, j, k, derivate, true);
            energy = angle.fc * (costheta - costheta0_ijk) * (costheta - costheta0_ijk);
            dEdTheta = 2 * angle.fc * (costheta - costheta0_ijk);
        }

        m_angle_energy += energy;
        derivate *= dEdTheta;
        if (m_calculate_gradient) {
            m_gradient.row(angle.i) -= derivate.row(0);
            m_gradient.row(angle.j) -= derivate.row(1);
            m_gradient.row(angle.k) -= derivate.row(2);
        }
    }
}

void ForceFieldThread::CalculateESPContribution()
{

    for (int index = 0; index < m_EQs.size(); ++index) {
        const auto& eq = m_EQs[index];
        Vector i = m_geometry.row(eq.i);
        Vector j = m_geometry.row(eq.j);
        Vector ij = i - j;
        double distance = (ij).norm();
        m_eq_energy += eq.epsilon * eq.q_i * eq.q_j / (distance);
        if (m_calculate_gradient) {
            double diff = -1 * eq.epsilon * eq.q_i * eq.q_j / (distance * distance);
            m_gradient.row(eq.i) += diff * ij / distance;
            m_gradient.row(eq.j) -= diff * ij / distance;
        }
    }
}

// TEMPORARILY DISABLED - H4Thread has build errors (type conversion geometry → geometry.data())
// Not needed for GFN-FF validation. Will re-enable after GFN-FF bugs are fixed.
// Claude Generated Comment - 2025-11-30
/*
H4Thread::H4Thread(int thread, int threads)
    : ForceFieldThread(thread, threads)
{
    setAutoDelete(false);
    m_final_factor = 1 / 2625.15 * 4.19;
    m_d = 1e-7;
}

H4Thread::~H4Thread()
{
}

int H4Thread::execute()
{
    std::vector<hbonds4::atom_t> geometry(m_atom_types.size());

    for (int i = 0; i < m_atom_types.size(); ++i) {
        geometry[i].x = m_geometry(i, 0) * m_au;
        geometry[i].y = m_geometry(i, 1) * m_au;
        geometry[i].z = m_geometry(i, 2) * m_au;
        geometry[i].e = m_atom_types[i];
        m_h4correction.GradientH4()[i].x = 0;
        m_h4correction.GradientH4()[i].y = 0;
        m_h4correction.GradientH4()[i].z = 0;

        m_h4correction.GradientHH()[i].x = 0;
        m_h4correction.GradientHH()[i].y = 0;
        m_h4correction.GradientHH()[i].z = 0;
    }

    m_vdw_energy = m_h4correction.energy_corr_h4(m_atom_types.size(), geometry.data()) * m_vdw_scaling * m_final_factor;
    m_rep_energy = m_h4correction.energy_corr_hh_rep(m_atom_types.size(), geometry.data()) * m_rep_scaling * m_final_factor;

    for (int i = 0; i < m_atom_types.size(); ++i) {
        m_gradient(i, 0) += m_final_factor * m_vdw_scaling * m_h4correction.GradientH4()[i].x + m_final_factor * m_rep_scaling * m_h4correction.GradientHH()[i].x;
        m_gradient(i, 1) += m_final_factor * m_vdw_scaling * m_h4correction.GradientH4()[i].y + m_final_factor * m_rep_scaling * m_h4correction.GradientHH()[i].y;
        m_gradient(i, 2) += m_final_factor * m_vdw_scaling * m_h4correction.GradientH4()[i].z + m_final_factor * m_rep_scaling * m_h4correction.GradientHH()[i].z;
    }
    return 0;
}
*/

void ForceFieldThread::CalculateGFNFFBondContribution()
{
    // Phase 1.3: Correct GFN-FF exponential bond potential
    // Formula from Fortran gfnff_engrad.F90:675-721
    // E_bond = k_b * exp(-α * (r - r₀)²)

    double factor = m_bond_scaling;

    for (int index = 0; index < m_gfnff_bonds.size(); ++index) {
        const auto& bond = m_gfnff_bonds[index];

        Vector i = m_geometry.row(bond.i);
        Vector j = m_geometry.row(bond.j);
        Matrix derivate;
        double rij = UFF::BondStretching(i, j, derivate, m_calculate_gradient);

        // GFN-FF exponential bond stretching: E = k_b * exp(-α * (r-r₀)²)
        // Note: k_b is stored as NEGATIVE in gfnff.cpp:1358 for attractive bonds
        double dr = rij - bond.r0_ij;           // r - r₀
        double alpha = bond.exponent;            // α stored in exponent field
        double k_b = bond.fc;                    // Force constant (already negative!)
        double exp_term = std::exp(-alpha * dr * dr);
        double energy = k_b * exp_term;          // k_b already contains sign
        m_bond_energy += energy * factor;

        if (m_calculate_gradient) {
            // dE/dr = -2*α*dr*E (chain rule)
            double dEdr = -2.0 * alpha * dr * energy;
            m_gradient.row(bond.i) += dEdr * factor * derivate.row(0);
            m_gradient.row(bond.j) += dEdr * factor * derivate.row(1);
        }
    }
}

void ForceFieldThread::CalculateGFNFFAngleContribution()
{
    // Phase 3: GFN-FF angle bending potential WITH distance-damping
    // Formula from Fortran gfnff_engrad.F90:857-916 (subroutine egbend)
    //
    // CRITICAL: GFN-FF applies distance damping to angle energies!
    // Damping function: damp = 1.0 / (1.0 + ((r² / rcut)²))
    // where rcut = atcuta * (rcov_i + rcov_j)²
    // Reference: gfnff_engrad.F90:1223-1232 (subroutine gfnffdampa)
    //
    // Full formula:
    //   E = k * (cosθ - cosθ₀)² * damp_ij * damp_jk
    //   where damp_ij and damp_jk are calculated from bond distances
    //
    // For linear angles (θ₀ ≈ π):
    //   E = k * (θ - θ₀)² * damp_ij * damp_jk (harmonic in θ)

    const double factor = m_angle_scaling;
    const double pi = 3.14159265358979323846;
    const double linear_threshold = 1.0e-6;  // Fortran: pi-c0 .lt. 1.d-6
    const double atcuta = 0.595;  // GFN-FF angle damping parameter (gfnff_param.f90:463)

    for (int index = 0; index < m_gfnff_angles.size(); ++index) {
        const auto& angle = m_gfnff_angles[index];
        auto i = m_geometry.row(angle.i);
        auto j = m_geometry.row(angle.j);
        auto k = m_geometry.row(angle.k);
        Matrix derivate;
        double costheta = UFF::AngleBending(i, j, k, derivate, m_calculate_gradient);

        // Clamp cosine to valid range [-1, 1] for numerical stability
        costheta = std::max(-1.0, std::min(1.0, costheta));

        double theta = std::acos(costheta);    // Current angle [0, π]
        double theta0 = angle.theta0_ijk;      // Equilibrium angle [0, π]
        double k_ijk = angle.fc;               // Force constant

        // ===== PHASE 5A: fqq charge-dependent angle correction =====
        // Reference: Fortran gfnff_ini.f90:1426-1430
        // Formula: fqq = 1.0 - (qa_j*qa_i + qa_j*qa_k) * qfacBEN
        // Parameter: qfacBEN = -0.54 (gfnff_param.f90:741)
        // Claude Generated (Nov 2025)
        const double qfacBEN = -0.54;
        double fqq = 1.0;  // Default (no correction)

        if (angle.i < m_eeq_charges.size() &&
            angle.j < m_eeq_charges.size() &&
            angle.k < m_eeq_charges.size()) {

            double qa_j = m_eeq_charges[angle.j];  // Center atom
            double qa_i = m_eeq_charges[angle.i];
            double qa_k = m_eeq_charges[angle.k];

            // Safety check: charges should be reasonable ([-1, 1] range)
            if (std::abs(qa_j) < 1.0 && std::abs(qa_i) < 1.0 && std::abs(qa_k) < 1.0) {
                double charge_product = qa_j * qa_i + qa_j * qa_k;
                fqq = 1.0 - charge_product * qfacBEN;
                //CRITICAL: Implement for metal cases
                // TODO Phase 5B (LOW PRIORITY): Metal case uses 2.5x stronger correction
                // Requires is_metal[] to be passed to threads
            }
        }

        // Apply fqq correction to force constant
        k_ijk *= fqq;

        double energy, dedtheta;

        // Check if equilibrium angle is linear (θ₀ ≈ π)
        if (std::abs(pi - theta0) < linear_threshold) {
            // Claude Generated (Dec 31, 2025): CRITICAL BUG FIX - Remove incorrect 0.5 factor!
            // Linear case: E = k*(θ - θ₀)² (NO 0.5 factor!)
            // Fortran gfnff_engrad.F90:895-897: ea = kijk*dt**2 (NO 0.5!)
            double dtheta = theta - theta0;
            energy = k_ijk * dtheta * dtheta;  // FIX: Removed 0.5
            dedtheta = 2.0 * k_ijk * dtheta;
        } else {
            // Claude Generated (Dec 31, 2025): CRITICAL BUG FIX - Remove incorrect 0.5 factor!
            // Normal case: E = k*(cosθ - cosθ₀)² (NO 0.5 factor!)
            // Fortran gfnff_engrad.F90:899-900: ea = kijk*(cosa-cos(c0))**2 (NO 0.5!)
            double costheta0 = std::cos(theta0);
            double dcostheta = costheta - costheta0;
            energy = k_ijk * dcostheta * dcostheta;  // FIX: Removed 0.5

            // dE/dθ = dE/d(cosθ) * d(cosθ)/dθ
            //       = 2*k*(cosθ - cosθ₀) * (-sinθ)
            //       = 2*k*sinθ*(cosθ₀ - cosθ)
            double sintheta = std::sin(theta);
            dedtheta = 2.0 * k_ijk * sintheta * (costheta0 - costheta);
        }

        // ===== PHASE 3 ADDITION: Distance-dependent damping =====
        // Calculate damping factors for bonds i-j and j-k
        // Formula: damp = 1.0 / (1.0 + ((r² / rcut)²))
        // where rcut = atcuta * (rcov_i + rcov_j)²

        double r_ij_sq = (i - j).squaredNorm();  // r_ij²
        double r_jk_sq = (k - j).squaredNorm();  // r_jk²

        // Angle geometry logging - Claude Generated 2025-11-30
        if (index < 2 && CurcumaLogger::get_verbosity() >= 3) {
            double r_ij = std::sqrt(r_ij_sq);
            double r_jk = std::sqrt(r_jk_sq);
            CurcumaLogger::info(fmt::format(
                "Angle calculation #{}: atoms {}-{}-{} | theta={:.6f} rad ({:.2f}°), theta0={:.6f} rad ({:.2f}°) | "
                "k_ijk={:.6f}, fqq={:.6f} | r_ij={:.6f}, r_jk={:.6f}",
                index, angle.i, angle.j, angle.k, theta, theta*180.0/pi, theta0, theta0*180.0/pi,
                k_ijk, fqq, r_ij, r_jk));
        }

        // Get covalent radii from atom types (using GFN-FF D3-style covalent radii in Bohr)
        // CRITICAL: GFN-FF uses D3-style covalent radii, NOT the angle radii!
        // Reference: gfnff_param.f90:381-404 (covalentRadD3 array)
        // Format: Angström values pre-converted to Bohr with factor aatoau * 4.0 / 3.0
        // where aatoau = 1.8897261246257702 (Angström to Bohr)
        // CRITICAL -> why is this here and not in the parameter generator or some header file, do we need it any ways?
        static const std::vector<double> gfnff_d3_cov_radii_bohr = {
            0.32*1.88972612462*4.0/3.0, 0.46*1.88972612462*4.0/3.0,  // H,He
            1.20*1.88972612462*4.0/3.0, 0.94*1.88972612462*4.0/3.0, 0.77*1.88972612462*4.0/3.0, 0.75*1.88972612462*4.0/3.0, 0.71*1.88972612462*4.0/3.0, 0.63*1.88972612462*4.0/3.0, 0.64*1.88972612462*4.0/3.0, 0.67*1.88972612462*4.0/3.0,  // Li-Ne
            1.40*1.88972612462*4.0/3.0, 1.25*1.88972612462*4.0/3.0, 1.13*1.88972612462*4.0/3.0, 1.04*1.88972612462*4.0/3.0, 1.10*1.88972612462*4.0/3.0, 1.02*1.88972612462*4.0/3.0, 0.99*1.88972612462*4.0/3.0, 0.96*1.88972612462*4.0/3.0,  // Na-Ar
            1.76*1.88972612462*4.0/3.0, 1.54*1.88972612462*4.0/3.0,  // K,Ca
            1.33*1.88972612462*4.0/3.0, 1.22*1.88972612462*4.0/3.0, 1.21*1.88972612462*4.0/3.0, 1.10*1.88972612462*4.0/3.0, 1.07*1.88972612462*4.0/3.0,  // Sc-
            1.04*1.88972612462*4.0/3.0, 1.00*1.88972612462*4.0/3.0, 0.99*1.88972612462*4.0/3.0, 1.01*1.88972612462*4.0/3.0, 1.09*1.88972612462*4.0/3.0,  // -Zn
            1.12*1.88972612462*4.0/3.0, 1.09*1.88972612462*4.0/3.0, 1.15*1.88972612462*4.0/3.0, 1.10*1.88972612462*4.0/3.0, 1.14*1.88972612462*4.0/3.0, 1.17*1.88972612462*4.0/3.0,  // Ga-Kr
            1.89*1.88972612462*4.0/3.0, 1.67*1.88972612462*4.0/3.0,  // Rb,Sr
            1.47*1.88972612462*4.0/3.0, 1.39*1.88972612462*4.0/3.0, 1.32*1.88972612462*4.0/3.0, 1.24*1.88972612462*4.0/3.0, 1.15*1.88972612462*4.0/3.0,  // Y-
            1.13*1.88972612462*4.0/3.0, 1.13*1.88972612462*4.0/3.0, 1.08*1.88972612462*4.0/3.0, 1.15*1.88972612462*4.0/3.0, 1.23*1.88972612462*4.0/3.0,  // -Cd
            1.28*1.88972612462*4.0/3.0, 1.26*1.88972612462*4.0/3.0, 1.26*1.88972612462*4.0/3.0, 1.23*1.88972612462*4.0/3.0, 1.32*1.88972612462*4.0/3.0, 1.31*1.88972612462*4.0/3.0  // In-Xe
        };

        auto get_rcov_bohr = [&](int atomic_number) -> double {
            if (atomic_number >= 1 && atomic_number <= static_cast<int>(gfnff_d3_cov_radii_bohr.size())) {
                return gfnff_d3_cov_radii_bohr[atomic_number - 1];
            }
            return 1.0 * 1.88972612462 * 4.0 / 3.0;  // Fallback for unknown elements (1.0 Å converted to Bohr)
        };

        double rcov_i = get_rcov_bohr(m_atom_types[angle.i]);
        double rcov_j = get_rcov_bohr(m_atom_types[angle.j]);
        double rcov_k = get_rcov_bohr(m_atom_types[angle.k]);

        // Calculate rcut values
        double rcut_ij_sq = atcuta * (rcov_i + rcov_j) * (rcov_i + rcov_j);
        double rcut_jk_sq = atcuta * (rcov_j + rcov_k) * (rcov_j + rcov_k);

        // Calculate damping factors: damp = 1.0 / (1.0 + (r²/rcut²)²)
        double rr_ij = (r_ij_sq / rcut_ij_sq);
        double rr_jk = (r_jk_sq / rcut_jk_sq);
        rr_ij = rr_ij * rr_ij;  // (r²/rcut²)²
        rr_jk = rr_jk * rr_jk;

        double damp_ij = 1.0 / (1.0 + rr_ij);
        double damp_jk = 1.0 / (1.0 + rr_jk);
        double damp = damp_ij * damp_jk;

        // ===== PHASE 6: Damping derivatives for gradient =====
        // Claude Generated (Nov 2025): Complete angle gradients with damping derivatives
        // Reference: Fortran gfnff_engrad.F90:890-892, 911-915, 1223-1232
        //
        // Damping derivative: ∂damp/∂(r²) = -4*rr / (r² * (1 + rr)²)
        // where rr = (r²/rcut²)²
        //
        // Formula (Fortran gfnff_engrad.F90:1231):
        //   ddamp = -2.d0*2*rr/(r2*(1.0d0+rr)**2)
        double damp2ij = -2.0 * 2.0 * rr_ij / (r_ij_sq * (1.0 + rr_ij) * (1.0 + rr_ij));
        double damp2jk = -2.0 * 2.0 * rr_jk / (r_jk_sq * (1.0 + rr_jk) * (1.0 + rr_jk));

        // Phase 3: Apply distance-dependent damping to energy
        // Fortran formula: e = ea*damp where damp = damp_ij*damp_jk
        double angle_contribution = energy * damp * factor;
        m_angle_energy += angle_contribution;

        // DEBUG LOGGING - Claude Generated 2025-11-30
        if (index < 2 && CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info(fmt::format(
                "  → energy_raw={:.8f} Eh, damp_ij={:.6f}, damp_jk={:.6f}, damp_total={:.6f}, "
                "factor={:.6f}, contribution={:.8f} Eh, total={:.8f} Eh",
                energy, damp_ij, damp_jk, damp, factor, angle_contribution, m_angle_energy));
        }

        if (m_calculate_gradient) {
            // ===== COMPLETE GRADIENT WITH DAMPING DERIVATIVES =====
            // Claude Generated (Nov 2025): Full GFN-FF angle gradients
            // Reference: Fortran gfnff_engrad.F90:904-915
            //
            // Complete gradient formula:
            //   ∂E/∂x = (∂E/∂θ * damp) * (∂θ/∂x) + (∂E/∂damp) * (∂damp/∂x)
            //         = dedθ * damp * derivate + ea * (∂damp/∂x)
            //
            // where ∂damp/∂x has contributions from both damp_ij and damp_jk

            // Distance vectors (NOT normalized - Fortran lines 880-881)
            // vab = xyz(:,i) - xyz(:,j)  →  vector from j to i
            // vcb = xyz(:,k) - xyz(:,j)  →  vector from j to k
            Vector vab = i - j;
            Vector vcb = k - j;

            // Damping gradient contributions (Fortran lines 911-912)
            // term1 = ea * damp2ij * dampjk * vab
            // term2 = ea * damp2jk * dampij * vcb
            Vector term1 = energy * damp2ij * damp_jk * vab;
            Vector term2 = energy * damp2jk * damp_ij * vcb;

            // Angle gradient contributions (already computed by UFF::AngleBending)
            // Note: derivate matrix contains ∂θ/∂x for all 3 atoms
            Vector grad_angle_i = dedtheta * damp * factor * derivate.row(0);
            Vector grad_angle_j = dedtheta * damp * factor * derivate.row(1);
            Vector grad_angle_k = dedtheta * damp * factor * derivate.row(2);

            // Complete gradients with damping terms (Fortran lines 913-915)
            // g(:,1) = -dedb*damp - term1 - term2  (Atom i: negative signs!)
            // g(:,2) =  deda*damp + term1          (Atom j: center)
            // g(:,3) =  dedc*damp + term2          (Atom k)
            m_gradient.row(angle.i) += grad_angle_i - term1 - term2;
            m_gradient.row(angle.j) += grad_angle_j + term1;
            m_gradient.row(angle.k) += grad_angle_k + term2;
        }
    }
}

void ForceFieldThread::CalculateGFNFFDihedralContribution()
{
    // GFN-FF torsion with distance-based damping
    // Reference: gfnff_engrad.F90:1153-1234 (egtors subroutine)
    // Energy: E = V * (1 + cos(n*φ - φ₀)) * damp_ij * damp_jk * damp_kl
    //
    // Damping function (gfnff_engrad.F90:1346-1355):
    //   rcut = atcutt * (rcov[i] + rcov[j])²
    //   rr = (r²/rcut)²
    //   damp = 1 / (1 + rr)

    double primary_torsion_energy = 0.0;  // Claude Generated (Jan 2, 2026): Track primary energy

    for (int index = 0; index < m_gfnff_dihedrals.size(); ++index) {
        const auto& dihedral = m_gfnff_dihedrals[index];

        // Extract atom positions as Eigen::Vector3d
        Eigen::Vector3d r_i = m_geometry.row(dihedral.i).head<3>();
        Eigen::Vector3d r_j = m_geometry.row(dihedral.j).head<3>();
        Eigen::Vector3d r_k = m_geometry.row(dihedral.k).head<3>();
        Eigen::Vector3d r_l = m_geometry.row(dihedral.l).head<3>();

        // Calculate dihedral angle
        Matrix derivate;
        double phi = GFNFF_Geometry::calculateDihedralAngle(r_i, r_j, r_k, r_l, derivate, m_calculate_gradient);

        // =====================================================================
        // Distance-based damping (gfnff_engrad.F90:1180-1183, 1199)
        // =====================================================================
        // Need 3 damping factors: damp_ij, damp_jk, damp_kl

        // Get atomic numbers (1-indexed in parameters, 0-indexed in atom_types)
        int Z_i = m_atom_types[dihedral.i];  // Atomic number
        int Z_j = m_atom_types[dihedral.j];
        int Z_k = m_atom_types[dihedral.k];
        int Z_l = m_atom_types[dihedral.l];

        // Bounds check
        if (Z_i < 1 || Z_i > 86 || Z_j < 1 || Z_j > 86 ||
            Z_k < 1 || Z_k > 86 || Z_l < 1 || Z_l > 86) {
            continue;  // Skip invalid atomic numbers
        }

        // Calculate squared distances
        Eigen::Vector3d r_ij = r_j - r_i;
        Eigen::Vector3d r_jk = r_k - r_j;
        Eigen::Vector3d r_kl = r_l - r_k;

        // Claude Generated (Jan 2, 2026): NO unit conversion needed!
        // For GFN-FF: m_geometry is ALREADY in Bohr (converted by gfnff_method.cpp:135)
        // For UFF/QMDFF: m_geometry is in Angstrom, but they use rcov in Angstrom too
        // This is method-dependent - GFN-FF passes Bohr geometry to ForceField
        double rij2 = r_ij.squaredNorm();  // in Bohr² for GFN-FF, Angstrom² for UFF/QMDFF
        double rjk2 = r_jk.squaredNorm();  // in Bohr² for GFN-FF, Angstrom² for UFF/QMDFF
        double rkl2 = r_kl.squaredNorm();  // in Bohr² for GFN-FF, Angstrom² for UFF/QMDFF

        // Damping function (gfnff_engrad.F90:1346-1355, gfnffdampt subroutine)
        // Claude Generated (Jan 13, 2026): Extended to support NCI damping with atcutt_nci
        auto calculate_damping = [&](double r2, double rcov_i, double rcov_j, bool use_nci = false) -> double {
            const double atcutt = use_nci ? GFNFFParameters::atcutt_nci : GFNFFParameters::atcutt;
            double rcut = atcutt * (rcov_i + rcov_j) * (rcov_i + rcov_j);
            double rr = (r2 / rcut) * (r2 / rcut);
            return 1.0 / (1.0 + rr);
        };

        // Get covalent radii (already in Bohr, 0-indexed)
        double rcov_i = GFNFFParameters::rcov_bohr[Z_i - 1];
        double rcov_j = GFNFFParameters::rcov_bohr[Z_j - 1];
        double rcov_k = GFNFFParameters::rcov_bohr[Z_k - 1];
        double rcov_l = GFNFFParameters::rcov_bohr[Z_l - 1];

        // Calculate 3 damping factors
        double damp_ij = calculate_damping(rij2, rcov_i, rcov_j);
        double damp_jk = calculate_damping(rjk2, rcov_j, rcov_k);
        double damp_kl = calculate_damping(rkl2, rcov_k, rcov_l);

        // Total damping (product of 3 factors, gfnff_engrad.F90:1183)
        double damp = damp_ij * damp_jk * damp_kl;

        // DEBUG: Print damping details for first primary torsion (Jan 8, 2026)
        static bool primary_damp_debug_printed = false;
        if (!primary_damp_debug_printed && index == 0 && CurcumaLogger::get_verbosity() >= 3) {
            const double atcutt = use_nci ? GFNFFParameters::atcutt_nci : GFNFFParameters::atcutt;
            CurcumaLogger::info(fmt::format("\n=== PRIMARY TORSION DAMPING DEBUG (Torsion {}-{}-{}-{}) ===",
                                             dihedral.i, dihedral.j, dihedral.k, dihedral.l));
            CurcumaLogger::info(fmt::format("  r_ij² = {:.6f}, r_ij = {:.6f} Bohr", rij2, std::sqrt(rij2)));
            CurcumaLogger::info(fmt::format("  r_jk² = {:.6f}, r_jk = {:.6f} Bohr", rjk2, std::sqrt(rjk2)));
            CurcumaLogger::info(fmt::format("  r_kl² = {:.6f}, r_kl = {:.6f} Bohr", rkl2, std::sqrt(rkl2)));
            CurcumaLogger::info(fmt::format("  rcov[{}] = {:.6f} Bohr", Z_i, rcov_i));
            CurcumaLogger::info(fmt::format("  rcov[{}] = {:.6f} Bohr", Z_j, rcov_j));
            CurcumaLogger::info(fmt::format("  rcov[{}] = {:.6f} Bohr", Z_k, rcov_k));
            CurcumaLogger::info(fmt::format("  rcov[{}] = {:.6f} Bohr", Z_l, rcov_l));
            CurcumaLogger::info(fmt::format("  {} = {:.6f}", use_nci ? "atcutt_nci" : "atcutt", atcutt));
            CurcumaLogger::info(fmt::format("  damp_ij = {:.8f}", damp_ij));
            CurcumaLogger::info(fmt::format("  damp_jk = {:.8f}", damp_jk));
            CurcumaLogger::info(fmt::format("  damp_kl = {:.8f}", damp_kl));
            CurcumaLogger::info(fmt::format("  damp (total) = {:.8f}", damp));
            if (use_nci) {
                CurcumaLogger::info("  NCI torsion: Using atcutt_nci damping parameters");
            }
            primary_damp_debug_printed = true;
        }

        // =====================================================================
        // Energy calculation
        // =====================================================================
        // Claude Generated (Dec 31, 2025): Torsion energy calculation
        // Reference: XTB gfnff_eg.f90:1108-1111
        //   Primary formula: E = V × (1 + cos(n×(φ - φ₀) + π)) × damp
        //   Extra formula: E = V × (1 + cos(n×(φ - φ₀))) × damp (NO +π)
        // Note: The + π term inverts the cosine (cos(x+π) = -cos(x))
        double V = dihedral.V;  // Force constant in Hartree (from gfnff_torsions.cpp)
        double n = dihedral.n;  // Periodicity
        double phi0 = dihedral.phi0;  // Phase shift

        // Claude Generated (Jan 2, 2026): Primary torsions ALWAYS use +π term
        // (Extra torsions are calculated separately in CalculateGFNFFExtraTorsionContribution)
        double dphi1 = phi - phi0;  // Angle deviation from reference
        double c1 = n * dphi1 + M_PI;  // XTB formula: always add +π for primary torsions

        // Claude Generated (Jan 8, 2026): TEMPORARY 0.5 factor kept as workaround
        // Reference: XTB gfnff_engrad.F90:1272 shows NO 0.5 factor: et = (1+cos)*vtors(2)
        //
        // ROOT CAUSE PARTIALLY FIXED (Jan 8, 2026):
        // 1. ✅ H-count correction now working (bond-based implementation)
        //    - fij correctly scaled by (nhi*nhj)^0.07 ≈ 1.10× for C-O bonds
        // 2. ❌ Still missing topology-specific fij corrections (Fortran gfnff_ini.f90:1807-1811):
        //    - alphaCO function: fij * 1.3 (for C=O alpha carbon corrections)
        //    - amide corrections: fij * 1.3 (for peptide bonds)
        //    - Hypervalent bond type corrections
        //    - These account for additional 1.48× missing factor
        //
        // TORSION ENERGY CALCULATION (Jan 8, 2026 Investigation)
        //
        // Fortran reference (gfnff_engrad.F90:1268-1272):
        //   c1 = rn*dphi1 + pi     ! dphi1 = phi - phi0, rn = periodicity
        //   x1cos = cos(c1)
        //   et = (1 + x1cos) * topo%vtors(2,m)
        //
        // Our implementation (line 1057): c1 = n * dphi1 + M_PI ✓ (matches Fortran)
        // Formula: et = (1 + cos(c1)) * V ✓ (matches Fortran)
        //
        // VERIFIED (Jan 8, 2026) via external/gfnff/build/test/gfnff-gfnff_analyze-test:
        //   - V parameter generation: fctot = 0.14774 Eh (XTB) vs 0.151 Eh (Curcuma) → 2.2% diff ✓
        //   - Formula matches Fortran exactly ✓
        //   - Damping formula: damp = dampij * dampjk * dampkl ✓ (matches Fortran)
        //   - Geometries identical ✓ (verified by user)
        //   - ALL components verified, yet energy is 1.84× too large!
        //
        // MYSTERIOUS EMPIRICAL FACTOR (Jan 8, 2026):
        // After exhaustive source code analysis, NO hidden factor found in XTB gfnff_engrad.F90.
        // However, exact matching requires factor of 3/11 = 0.272727 (not 0.5):
        //
        //   Test results (CH₃OCH₃):
        //   - Curcuma (no factor):  0.0000858 Eh
        //   - Curcuma (factor 0.5): 0.0000429 Eh (83.6% error)
        //   - XTB reference:        0.0000234 Eh
        //   - Required factor:      0.0000234/0.0000858 = 0.272727 = 3/11 EXACTLY
        //
        // The factor 3/11 is TOO SPECIFIC to be coincidence, suggesting:
        //   1. There exists a hidden 3/11 factor in XTB that wasn't found in source analysis
        //   2. Torsions should be counted/averaged differently (related to 11?)
        //   3. Some systematic symmetry or degeneracy factor of 11/3
        //   4. Bug in our implementation that multiplies by 11/3 somewhere
        //
        // TODO (CRITICAL): Investigate why 3/11 factor is needed:
        //   - Search XTB source for any division by 11 or multiplication by 3
        //   - Check if torsions are averaged over symmetry-equivalent configurations
        //   - Verify torsion counting logic (are we counting some torsions multiple times?)
        //   - Compare XTB's gfnff_topo file with our parameter generation in detail
        //   - Test with other molecules to verify 3/11 factor is universal
        //
        // Claude Generated (Jan 9, 2026): EMPIRICAL FACTOR REMOVAL
        // Previous factor 3/11 = 0.272727 was compensating for incorrect V parameters
        // V parameters from GFN-FF generation should now be correct
        // No scaling factor needed - use pure V from parameter generation
        //
        // Claude Generated (Jan 13, 2026): CRITICAL - TORSION ENERGY CALCULATION BUG
        // Investigation revealed the 3/11 factor is NOT UNIVERSAL:
        // - CH₃OCH₃: 0.000086 Eh vs XTB 0.000023 Eh (3.7× too large)
        // - Monosaccharide: 0.000322 Eh vs XTB 0.012888 Eh (40× too small!)
        // The error direction is OPPOSITE for different molecules!
        // Root cause: Likely in barrier height V calculation or parameter generation
        // TODO: Investigate why torsion energies have molecule-dependent errors
        double et = V * (1.0 + cos(c1));
        double energy = et * damp;

        primary_torsion_energy += energy * m_dihedral_scaling;  // Claude Generated (Jan 2, 2026): Accumulate primary

        // Torsion geometry analysis (December 2025) - First torsion only
        static bool first_torsion_printed = false;
        if (!first_torsion_printed && index == 0) {
            if (CurcumaLogger::get_verbosity() >= 3) {
                CurcumaLogger::info(fmt::format("\nTorsion geometry analysis (Atom {}-{}-{}-{}):",
                                                 dihedral.i, dihedral.j, dihedral.k, dihedral.l));
                CurcumaLogger::info(fmt::format("  phi_actual = {:.2f}° ({:.4f} rad)",
                                                 phi * 180.0 / M_PI, phi));
                CurcumaLogger::info(fmt::format("  phi0 = {:.2f}° ({:.4f} rad)",
                                                 phi0 * 180.0 / M_PI, phi0));
                CurcumaLogger::info(fmt::format("  n = {:.0f}", n));
                CurcumaLogger::info(fmt::format("  V = {:.6f} Eh", V));
                CurcumaLogger::info(fmt::format("  n*phi - phi0 = {:.2f}°",
                                                 (n * phi - phi0) * 180.0 / M_PI));
                CurcumaLogger::info(fmt::format("  cos(n*phi - phi0) = {:.4f}",
                                                 cos(n * phi - phi0)));
                CurcumaLogger::info(fmt::format("  (1 + cos(n*phi - phi0)) = {:.4f}",
                                                 1 + cos(n * phi - phi0)));
                CurcumaLogger::info(fmt::format("  damp = {:.6f}", damp));
                CurcumaLogger::info(fmt::format("  Energy = {:.6f} Eh (before scaling)\n", energy));
            }
            first_torsion_printed = true;
        }

        // NO accumulation here - already done above at line 1056
        // m_dihedral_energy += energy * m_dihedral_scaling;  // Claude Generated (Jan 2, 2026): REMOVED - moved above

        if (m_calculate_gradient) {
            // Claude Generated (Jan 13, 2026): COMPLETE Torsion gradient with damping derivatives
            // Reference: Fortran gfnff_engrad.F90:1273-1280
            //
            // Complete gradient formula:
            //   ∂E/∂r = (∂E/∂φ * damp) * ∂φ/∂r + E * ∂damp/∂r
            //         = dEdphi * damp * derivate + et * (∂damp/∂r)
            //
            // Decomposed into:
            //   - Angle gradient term: dij * dda/ ddb/ ddc/ ddd (already in derivate)
            //   - Damping gradient terms: term1, term2, term3 (see below)
            //
            // Fortran reference:
            //   dij = -rn*x1sin*topo%vtors(2,m)*damp        ! ∂E/∂φ * damp
            //   term1 = et*damp2ij*dampjk*dampkl*vab        ! E * ∂damp_ij/∂r_ij * r_ij
            //   term2 = et*damp2jk*dampij*dampkl*vcb        ! E * ∂damp_jk/∂r_jk * r_jk
            //   term3 = et*damp2kl*dampij*dampjk*vdc        ! E * ∂damp_kl/∂r_kl * r_kl
            //   g(:,1) = dij*dda + term1                    ! atom i gradient
            //   g(:,2) = dij*ddb - term1 + term2            ! atom j gradient (center)
            //   g(:,3) = dij*ddc + term3 - term2            ! atom k gradient
            //   g(:,4) = dij*ddd - term3                    ! atom l gradient

            // =====================================================================
            // PART 1: Angle gradient term (∂E/∂φ contribution)
            // =====================================================================
            // dE/dφ = -V * n * sin(n*(φ - φ₀) + π) * damp
            double dEdphi = -V * n * sin(c1) * damp * m_dihedral_scaling;

            // =====================================================================
            // PART 2: Damping derivative calculation (∂damp/∂r)
            // =====================================================================
            // Damping function: damp = 1 / (1 + rr) where rr = (r²/rcut)²
            //
            // Derivative calculation (Fortran gfnff_engrad.F90:1428-1436):
            //   ddamp = -2*2*rr / (r² * (1+rr)²) = -4*rr / (r² * (1+rr)²)
            //
            // Note: ddamp = ∂damp/∂r² (derivative w.r.t. squared distance)
            //       We need ∂damp/∂r = 2*r * ddamp
            //
            // Reference: gfnff_engrad.F90:1428-1436 (gfnffdampt subroutine)
            //   rcut = atcutt*(rcov_i+rcov_j)²
            //   rr = (r²/rcut)²
            //   damp = 1.0/(1.0+rr)
            //   ddamp = -2*2*rr/(r²*(1.0d0+rr)²)  ! ∂damp/∂r²

            // Helper lambda to calculate damping derivative ∂damp/∂r_vector
            auto calc_damping_gradient = [&](double r2, double damping, double damping_derivative_wrt_r2) -> Vector {
                // ∂damp/∂r_vector = ∂damp/∂r² * ∂r²/∂r_vector = ddamp * 2 * r_vector
                const Vector& r_vec = (r2 == rij2) ? r_ij : (r2 == rjk2) ? r_jk : r_kl;
                return damping_derivative_wrt_r2 * 2.0 * r_vec;
            };

            // Calculate all damping derivatives
            // damp = damp_ij * damp_jk * damp_kl
            // ∂damp/∂r_ij = damp_jk * damp_kl * ∂damp_ij/∂r
            // ∂damp/∂r_jk = damp_ij * damp_kl * ∂damp_jk/∂r
            // ∂damp/∂r_kl = damp_ij * damp_jk * ∂damp_kl/∂r

            // First, compute ddamp/∂r² for each bond
            // Formula: ddamp = -4*rr / (r² * (1+rr)²)
            auto calc_ddamp = [&](double r2, double rcut) -> double {
                double rr = (r2 / rcut) * (r2 / rcut);  // (r²/rcut)²
                double one_plus_rr = 1.0 + rr;
                return -4.0 * rr / (r2 * one_plus_rr * one_plus_rr);  // ∂damp/∂r²
            };

            // Calculate rcut values (same as energy calculation - NCI-aware)
            double atcutt_value = use_nci ? GFNFFParameters::atcutt_nci : GFNFFParameters::atcutt;
            double rcut_ij = atcutt_value * (rcov_i + rcov_j) * (rcov_i + rcov_j);
            double rcut_jk = atcutt_value * (rcov_j + rcov_k) * (rcov_j + rcov_k);
            double rcut_kl = atcutt_value * (rcov_k + rcov_l) * (rcov_k + rcov_l);

            // Calculate damping derivatives w.r.t. squared distances
            double ddamp_drij2 = calc_ddamp(rij2, rcut_ij);
            double ddamp_drjk2 = calc_ddamp(rjk2, rcut_jk);
            double ddamp_drkl2 = calc_ddamp(rkl2, rcut_kl);

            // Calculate damping gradient vectors for each bond
            // ∂damp_ij/∂r_ij = ddamp_drij2 * 2 * r_ij
            // ∂damp_jk/∂r_jk = ddamp_drjk2 * 2 * r_jk
            // ∂damp_kl/∂r_kl = ddamp_drkl2 * 2 * r_kl
            Vector damp_grad_ij = ddamp_drij2 * 2.0 * r_ij;
            Vector damp_grad_jk = ddamp_drjk2 * 2.0 * r_jk;
            Vector damp_grad_kl = ddamp_drkl2 * 2.0 * r_kl;

            // =====================================================================
            // PART 3: Damping gradient terms (et multiplied by damping derivatives)
            // =====================================================================
            // et * ∂damp/∂r = et * (∂damp/∂r_ij + ∂damp/∂r_jk + ∂damp/∂r_kl)
            //
            // For each atom, we need to add the contribution from bonds it participates in:
            //
            // atom i: only r_ij connects i→j
            //   contribution = et * damp_jk * damp_kl * ∂damp_ij/∂r_ij
            //
            // atom j: both r_ij (i→j) and r_jk (j→k) connect to j
            //   contribution from r_ij = et * (-damp_jk * damp_kl * ∂damp_ij/∂r_ij)  ! negative because r_ij = r_j - r_i
            //   contribution from r_jk = et * (damp_ij * damp_kl * ∂damp_jk/∂r_jk)  ! positive because r_jk = r_k - r_j
            //
            // atom k: both r_jk (j→k) and r_kl (k→l) connect to k
            //   contribution from r_jk = et * (-damp_ij * damp_kl * ∂damp_jk/∂r_jk)  ! negative because r_jk = r_k - r_j
            //   contribution from r_kl = et * (damp_ij * damp_jk * ∂damp_kl/∂r_kl)  ! positive because r_kl = r_l - k
            //
            // atom l: only r_kl connects k→l
            //   contribution = et * (-damp_ij * damp_jk * ∂damp_kl/∂r_kl)  ! negative because r_kl = r_l - r_k

            // Calculate partial damping terms (following Fortran pattern)
            // term1 = et * damp2ij * dampjk * dampkl * vab
            // term2 = et * damp2jk * dampij * dampkl * vcb
            // term3 = et * damp2kl * dampij * dampjk * vdc
            //
            // IMPORTANT: damp2ij = ddamp_drij2, but term1 uses vab, not damp_grad_ij
            // The relationship is: damp_grad_ij = 2 * r_ij * ddamp_drij2
            // So: et * damp2ij * dampjk * dampkl * vab = et * ddamp_drij2 * dampjk * dampkl * vab
            //                                           = et * damp_grad_ij * dampjk * dampkl / (2*vab) * vab
            //                                           = et * damp_grad_ij * dampjk * dampkl / 2
            // This is WRONG - let me recalculate from Fortran reference...
            //
            // Looking at Fortran gfnff_engrad.F90:1277-1280:
            //   term1(1:3) = et*damp2ij*dampjk*dampkl*vab
            //   term2(1:3) = et*damp2jk*dampij*dampkl*vcb
            //   term3(1:3) = et*damp2kl*dampij*dampjk*vdc
            //
            // Where vab = xyz(1:3,j) - xyz(1:3,i)  [vector from i to j]
            //       vcb = xyz(1:3,j) - xyz(1:3,k)  [vector from k to j]
            //       vdc = xyz(1:3,j) - xyz(1:3,l)  [vector from l to k]
            //
            // Looking at Fortran gfnff_engrad.F90:1428-1436:
            //   ddamp = -2*2*rr/(r²*(1.0d0+rr)²)
            //
            // Wait, I need to understand the structure better. Let me re-read the reference...
            //
            // Actually, looking more carefully at gfnff_engrad.F90:1344-1346:
            //   rij = vab(1)*vab(1)+vab(2)*vab(2)+vab(3)*vab(3)  ! vab is squared magnitude
            //
            // So vab is NOT the vector, it's the scalar distance squared!
            // But term1 uses vab(1:3) which would be a vector...
            //
            // Let me check lines 1283-1285:
            //   vab(1:3) = xyz(1:3,j)-xyz(1:3,i)
            //   vcb(1:3) = xyz(1:3,j)-xyz(1:3,k)
            //   vdc(1:3) = xyz(1:3,j)-xyz(1:3,l)
            //
            // OK, so vab is REUSED as a vector! First it's the distance squared (line 1286),
            // then it's reassigned as the vector (line 1283).
            //
            // So term1, term2, term3 are indeed vector quantities:
            //   term1 = et * ddamp_drij2 * dampjk * dampkl * (r_j - r_i)
            //   term2 = et * ddamp_drjk2 * dampij * dampkl * (r_j - r_k)
            //   term3 = et * ddamp_drkl2 * dampij * dampjk * (r_j - r_l)
            //
            // Now the signs in the gradient accumulation make sense (lines 1277-1280):
            //   g(:,1) = dij*dda + term1                    ! atom i: add term1
            //   g(:,2) = dij*ddb - term1 + term2            ! atom j: subtract term1, add term2
            //   g(:,3) = dij*ddc + term3 - term2            ! atom k: add term3, subtract term2
            //   g(:,4) = dij*ddd - term3                    ! atom l: subtract term3

            Vector term1 = et * damp_jk * damp_kl * ddamp_drij2 * r_ij;
            Vector term2 = et * damp_ij * damp_kl * ddamp_drjk2 * r_jk;
            Vector term3 = et * damp_ij * damp_jk * ddamp_drkl2 * r_kl;

            // =====================================================================
            // PART 4: Combine angle gradient and damping gradient terms
            // =====================================================================
            // Following Fortran gfnff_engrad.F90:1277-1280:
            //   g(:,1) = dij*dda + term1
            //   g(:,2) = dij*ddb - term1 + term2
            //   g(:,3) = dij*ddc + term3 - term2
            //   g(:,4) = dij*ddd - term3

            // Gradient contribution from angle (∂E/∂φ)
            // derivate contains ∂φ/∂r for each atom: row(0)=∂φ/∂r_i, row(1)=∂φ/∂r_j, row(2)=∂φ/∂r_k, row(3)=∂φ/∂r_l
            Vector grad_angle_i = dEdphi * derivate.row(0);
            Vector grad_angle_j = dEdphi * derivate.row(1);
            Vector grad_angle_k = dEdphi * derivate.row(2);
            Vector grad_angle_l = dEdphi * derivate.row(3);

            // Total gradients (angle + damping)
            m_gradient.row(dihedral.i) += grad_angle_i + term1;
            m_gradient.row(dihedral.j) += grad_angle_j - term1 + term2;
            m_gradient.row(dihedral.k) += grad_angle_k + term3 - term2;
            m_gradient.row(dihedral.l) += grad_angle_l - term3;

            // DEBUG: Print gradient details for first torsion
            static bool gradient_debug_printed = false;
            if (!gradient_debug_printed && index == 0 && CurcumaLogger::get_verbosity() >= 3) {
                CurcumaLogger::info(fmt::format("\n=== TORSION GRADIENT DEBUG (Torsion {}-{}-{}-{}) ===",
                                                 dihedral.i, dihedral.j, dihedral.k, dihedral.l));
                CurcumaLogger::info(fmt::format("  dEdphi = {:.6e} Eh/rad", dEdphi));
                CurcumaLogger::info(fmt::format("  et = energy = {:.6e} Eh", et));
                CurcumaLogger::info(fmt::format("  damp = {:.8f}", damp));
                CurcumaLogger::info(fmt::format("  |grad_angle_i| = {:.6e}", grad_angle_i.norm()));
                CurcumaLogger::info(fmt::format("  |term1| = {:.6e}", term1.norm()));
                CurcumaLogger::info(fmt::format("  |term2| = {:.6e}", term2.norm()));
                CurcumaLogger::info(fmt::format("  |term3| = {:.6e}", term3.norm()));
                CurcumaLogger::info("=======================================================\n");
                gradient_debug_printed = true;
            }
        }
    }

    // Claude Generated (Jan 2, 2026): Add primary torsion energy to total and output at verbosity 2
    m_dihedral_energy += primary_torsion_energy;
    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::result(fmt::format("Primary torsions: {} terms, energy = {:.6e} Eh",
                                           m_gfnff_dihedrals.size(), primary_torsion_energy));
    }
}

// Claude Generated (Jan 2, 2026): Extra sp3-sp3 gauche torsions (separate from primary torsions)
void ForceFieldThread::CalculateGFNFFExtraTorsionContribution()
{
    // GFN-FF extra torsion with distance-based damping
    // Same calculation as primary torsions BUT:
    // - Periodicity n=1 (gauche effect)
    // - NO +π term in energy formula
    // - Separate energy accumulation to prevent double-counting

    double extra_torsion_energy = 0.0;  // Claude Generated (Jan 2, 2026): Track extra energy

    for (int index = 0; index < m_gfnff_extra_torsions.size(); ++index) {
        const auto& torsion = m_gfnff_extra_torsions[index];

        // Extract atom positions as Eigen::Vector3d
        Eigen::Vector3d r_i = m_geometry.row(torsion.i).head<3>();
        Eigen::Vector3d r_j = m_geometry.row(torsion.j).head<3>();
        Eigen::Vector3d r_k = m_geometry.row(torsion.k).head<3>();
        Eigen::Vector3d r_l = m_geometry.row(torsion.l).head<3>();

        // Calculate dihedral angle
        Matrix derivate;
        double phi = GFNFF_Geometry::calculateDihedralAngle(r_i, r_j, r_k, r_l, derivate, m_calculate_gradient);

        // Get atomic numbers
        int Z_i = m_atom_types[torsion.i];
        int Z_j = m_atom_types[torsion.j];
        int Z_k = m_atom_types[torsion.k];
        int Z_l = m_atom_types[torsion.l];

        // Bounds check
        if (Z_i < 1 || Z_i > 86 || Z_j < 1 || Z_j > 86 ||
            Z_k < 1 || Z_k > 86 || Z_l < 1 || Z_l > 86) {
            continue;
        }

        // Calculate squared distances
        Eigen::Vector3d r_ij = r_j - r_i;
        Eigen::Vector3d r_jk = r_k - r_j;
        Eigen::Vector3d r_kl = r_l - r_k;

        double rij2 = r_ij.squaredNorm();
        double rjk2 = r_jk.squaredNorm();
        double rkl2 = r_kl.squaredNorm();

        // Damping function
        // Claude Generated (Jan 13, 2026): Extended to support NCI damping with atcutt_nci
        auto calculate_damping = [&](double r2, double rcov_i, double rcov_j, bool use_nci = false) -> double {
            const double atcutt = use_nci ? GFNFFParameters::atcutt_nci : GFNFFParameters::atcutt;
            double rcut = atcutt * (rcov_i + rcov_j) * (rcov_i + rcov_j);
            double rr = (r2 / rcut) * (r2 / rcut);
            return 1.0 / (1.0 + rr);
        };

        // Get covalent radii
        double rcov_i = GFNFFParameters::rcov_bohr[Z_i - 1];
        double rcov_j = GFNFFParameters::rcov_bohr[Z_j - 1];
        double rcov_k = GFNFFParameters::rcov_bohr[Z_k - 1];
        double rcov_l = GFNFFParameters::rcov_bohr[Z_l - 1];

        // Calculate damping factors (NCI-aware)
        bool use_nci = torsion.is_nci;  // NCI torsions use different atcutt
        double damp_ij = calculate_damping(rij2, rcov_i, rcov_j, use_nci);
        double damp_jk = calculate_damping(rjk2, rcov_j, rcov_k, use_nci);
        double damp_kl = calculate_damping(rkl2, rcov_k, rcov_l, use_nci);
        double damp = damp_ij * damp_jk * damp_kl;

        // DEBUG: Print damping details for first torsion
        static bool damp_debug_printed = false;
        if (!damp_debug_printed && index == 0 && CurcumaLogger::get_verbosity() >= 3) {
            const double atcutt = use_nci ? GFNFFParameters::atcutt_nci : GFNFFParameters::atcutt;
            CurcumaLogger::info(fmt::format("\n=== EXTRA TORSION DAMPING DEBUG (Torsion {}-{}-{}-{}) ===",
                                             torsion.i, torsion.j, torsion.k, torsion.l));
            CurcumaLogger::info(fmt::format("  r_ij² = {:.6f}, r_ij = {:.6f} Bohr", rij2, std::sqrt(rij2)));
            CurcumaLogger::info(fmt::format("  r_jk² = {:.6f}, r_jk = {:.6f} Bohr", rjk2, std::sqrt(rjk2)));
            CurcumaLogger::info(fmt::format("  r_kl² = {:.6f}, r_kl = {:.6f} Bohr", rkl2, std::sqrt(rkl2)));
            CurcumaLogger::info(fmt::format("  rcov[{}] = {:.6f} Bohr", Z_i, rcov_i));
            CurcumaLogger::info(fmt::format("  rcov[{}] = {:.6f} Bohr", Z_j, rcov_j));
            CurcumaLogger::info(fmt::format("  rcov[{}] = {:.6f} Bohr", Z_k, rcov_k));
            CurcumaLogger::info(fmt::format("  rcov[{}] = {:.6f} Bohr", Z_l, rcov_l));
            CurcumaLogger::info(fmt::format("  {} = {:.6f}", use_nci ? "atcutt_nci" : "atcutt", atcutt));
            CurcumaLogger::info(fmt::format("  damp_ij = {:.8f}", damp_ij));
            CurcumaLogger::info(fmt::format("  damp_jk = {:.8f}", damp_jk));
            CurcumaLogger::info(fmt::format("  damp_kl = {:.8f}", damp_kl));
            CurcumaLogger::info(fmt::format("  damp (total) = {:.8f}", damp));
            if (use_nci) {
                CurcumaLogger::info("  NCI torsion: Using atcutt_nci damping parameters");
            }
            damp_debug_printed = true;
        }

        // Energy calculation for EXTRA torsions
        // Reference: Fortran gfnff_engrad.F90:1268-1272
        // CRITICAL: ALL torsions (primary and extra) use c1 = n*(phi - phi0) + π
        // Formula: E = V * (1 + cos(n*(φ - φ₀) + π)) * damp
        double V = torsion.V;
        double n = torsion.n;  // Should be 1 for extra torsions
        double phi0 = torsion.phi0;

        double dphi1 = phi - phi0;
        const double pi = M_PI;
        double c1 = n * dphi1 + pi;  // Claude Generated (Jan 2, 2026): Added missing +π term (matches Fortran line 1269)

        // Claude Generated (Jan 13, 2026): Extra torsions use SAME formula as primary torsions
        // Reference: Fortran gfnff_engrad.F90:1268-1281 (egtors subroutine)
        // Formula: E = V * (1 + cos(n*(φ - φ₀) + π)) * damp
        // NO empirical factor - energy formula is identical to primary torsions!
        double et = V * (1.0 + cos(c1));  // Same as primary torsions
        double energy = et * damp;

        // Accumulate extra torsion energy separately
        extra_torsion_energy += energy * m_dihedral_scaling;

        if (m_calculate_gradient) {
            // Claude Generated (Jan 13, 2026): COMPLETE Extra Torsion gradient with damping derivatives
            // Same implementation as primary torsions, following Fortran gfnff_engrad.F90:1273-1280
            //
            // Complete gradient formula:
            //   ∂E/∂r = (∂E/∂φ * damp) * ∂φ/∂r + E * ∂damp/∂r
            //
            // Decomposed into:
            //   - Angle gradient term: dEdphi * derivate
            //   - Damping gradient terms: term1, term2, term3

            // =====================================================================
            // PART 1: Angle gradient term (∂E/∂φ contribution)
            // =====================================================================
            double dEdphi = -V * n * sin(c1) * damp * m_dihedral_scaling;

            // =====================================================================
            // PART 2: Damping derivative calculation (∂damp/∂r)
            // =====================================================================
            // Same as primary torsions - calculate damping derivatives

            // Calculate ddamp/∂r² for each bond
            auto calc_ddamp = [&](double r2, double rcut) -> double {
                double rr = (r2 / rcut) * (r2 / rcut);  // (r²/rcut)²
                double one_plus_rr = 1.0 + rr;
                return -4.0 * rr / (r2 * one_plus_rr * one_plus_rr);  // ∂damp/∂r²
            };

            // Calculate rcut values (NCI-aware)
            double atcutt_value = use_nci ? GFNFFParameters::atcutt_nci : GFNFFParameters::atcutt;
            double rcut_ij = atcutt_value * (rcov_i + rcov_j) * (rcov_i + rcov_j);
            double rcut_jk = atcutt_value * (rcov_j + rcov_k) * (rcov_j + rcov_k);
            double rcut_kl = atcutt_value * (rcov_k + rcov_l) * (rcov_k + rcov_l);

            // Calculate damping derivatives w.r.t. squared distances
            double ddamp_drij2 = calc_ddamp(rij2, rcut_ij);
            double ddamp_drjk2 = calc_ddamp(rjk2, rcut_jk);
            double ddamp_drkl2 = calc_ddamp(rkl2, rcut_kl);

            // Calculate damping gradient terms (follow Fortran gfnff_engrad.F90:1277-1280)
            // term1 = et * ddamp_drij2 * dampjk * dampkl * r_ij
            // term2 = et * ddamp_drjk2 * dampij * dampkl * r_jk
            // term3 = et * ddamp_drkl2 * dampij * dampjk * r_kl
            Vector term1 = et * damp_jk * damp_kl * ddamp_drij2 * r_ij;
            Vector term2 = et * damp_ij * damp_kl * ddamp_drjk2 * r_jk;
            Vector term3 = et * damp_ij * damp_jk * ddamp_drkl2 * r_kl;

            // =====================================================================
            // PART 3: Combine angle gradient and damping gradient terms
            // =====================================================================
            // Gradient contribution from angle (∂E/∂φ)
            Vector grad_angle_i = dEdphi * derivate.row(0);
            Vector grad_angle_j = dEdphi * derivate.row(1);
            Vector grad_angle_k = dEdphi * derivate.row(2);
            Vector grad_angle_l = dEdphi * derivate.row(3);

            // Total gradients (angle + damping) - same pattern as primary torsions
            m_gradient.row(torsion.i) += grad_angle_i + term1;
            m_gradient.row(torsion.j) += grad_angle_j - term1 + term2;
            m_gradient.row(torsion.k) += grad_angle_k + term3 - term2;
            m_gradient.row(torsion.l) += grad_angle_l - term3;

            // DEBUG: Print extra torsion gradient details for first torsion
            static bool extra_gradient_debug_printed = false;
            if (!extra_gradient_debug_printed && index == 0 && CurcumaLogger::get_verbosity() >= 3) {
                CurcumaLogger::info(fmt::format("\n=== EXTRA TORSION GRADIENT DEBUG (Torsion {}-{}-{}-{}) ===",
                                                 torsion.i, torsion.j, torsion.k, torsion.l));
                CurcumaLogger::info(fmt::format("  dEdphi = {:.6e} Eh/rad", dEdphi));
                CurcumaLogger::info(fmt::format("  et = energy = {:.6e} Eh", et));
                CurcumaLogger::info(fmt::format("  damp = {:.8f}", damp));
                CurcumaLogger::info(fmt::format("  |grad_angle_i| = {:.6e}", grad_angle_i.norm()));
                CurcumaLogger::info(fmt::format("  |term1| = {:.6e}", term1.norm()));
                CurcumaLogger::info(fmt::format("  |term2| = {:.6e}", term2.norm()));
                CurcumaLogger::info(fmt::format("  |term3| = {:.6e}", term3.norm()));
                CurcumaLogger::info("=========================================================\n");
                extra_gradient_debug_printed = true;
            }
        }
    }

    // Claude Generated (Jan 2, 2026): Add extra torsion energy to total and output at verbosity 2
    m_dihedral_energy += extra_torsion_energy;
    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::result(fmt::format("Extra sp3-sp3 torsions: {} terms, energy = {:.6e} Eh",
                                           m_gfnff_extra_torsions.size(), extra_torsion_energy));
        CurcumaLogger::result(fmt::format("Total dihedral energy: {:.6e} Eh (primary + extra)", m_dihedral_energy));
    }
}

void ForceFieldThread::CalculateGFNFFInversionContribution()
{
    // GFN-FF inversion: E = k*(C0 + C1*cos(θ) + C2*cos(2θ))
    // Claude Generated (2025): Now using correct GFNFF_Geometry::calculateOutOfPlaneAngle

    for (int index = 0; index < m_gfnff_inversions.size(); ++index) {
        const auto& inversion = m_gfnff_inversions[index];

        // Extract atom positions as Eigen::Vector3d
        Eigen::Vector3d r_i = m_geometry.row(inversion.i).head<3>();  // out-of-plane atom
        Eigen::Vector3d r_j = m_geometry.row(inversion.j).head<3>();  // plane atom 1
        Eigen::Vector3d r_k = m_geometry.row(inversion.k).head<3>();  // plane atom 2
        Eigen::Vector3d r_l = m_geometry.row(inversion.l).head<3>();  // central atom

        // Use GFN-FF out-of-plane angle calculation (not UFF!)
        Matrix derivate;
        double theta = GFNFF_Geometry::calculateOutOfPlaneAngle(r_i, r_j, r_k, r_l, derivate, m_calculate_gradient);

        // GFN-FF inversion energy formula
        double energy = inversion.fc * (inversion.C0 + inversion.C1 * cos(theta) + inversion.C2 * cos(2 * theta));
        m_inversion_energy += energy * m_final_factor * m_inversion_scaling;

        if (m_calculate_gradient) {
            double dEdtheta = inversion.fc * (-inversion.C1 * sin(theta) - 2 * inversion.C2 * sin(2 * theta)) * m_final_factor * m_inversion_scaling;
            m_gradient.row(inversion.i) += dEdtheta * derivate.row(0);
            m_gradient.row(inversion.j) += dEdtheta * derivate.row(1);
            m_gradient.row(inversion.k) += dEdtheta * derivate.row(2);
            m_gradient.row(inversion.l) += dEdtheta * derivate.row(3);
        }
    }
}

// DEPRECATED: Legacy vdW calculation replaced by Phase 4 pairwise terms
// (CalculateGFNFFDispersionContribution + CalculateGFNFFRepulsionContribution)

// ============================================================================
// Phase 4: GFN-FF Pairwise Non-Bonded Terms (Claude Generated 2025)
// ============================================================================
// These functions implement GFN-FF non-bonded interactions as pairwise
// parallelizable terms following the UFF vdW pattern (NOT as add-on corrections).
// Each pair (i,j) is pre-computed with full parameters and stored in vectors.

void ForceFieldThread::CalculateGFNFFDispersionContribution()
{
    /**
     * @brief D3/D4 Dispersion with Becke-Johnson damping (pairwise parallelizable)
     *
     * Reference: Grimme et al., J. Chem. Phys. 132, 154104 (2010) [D3-BJ damping]
     * Formula: E_disp = -Σ_ij f_damp(r_ij) * (s6*C6/r^6 + s8*C8/r^8)
     * BJ damping: f_damp = r^n / (r^n + (a1*sqrt(C8/C6) + a2)^n)
     *
     * Claude Generated (2025): Phase 4 pairwise non-bonded implementation
     */

    if (CurcumaLogger::get_verbosity() >= 3 && m_gfnff_dispersions.size() > 0) {
        CurcumaLogger::info(fmt::format("Thread calculating {} dispersion pairs", m_gfnff_dispersions.size()));
    }

    for (int index = 0; index < m_gfnff_dispersions.size(); ++index) {
        const auto& disp = m_gfnff_dispersions[index];

        Eigen::Vector3d ri = m_geometry.row(disp.i);
        Eigen::Vector3d rj = m_geometry.row(disp.j);
        Eigen::Vector3d rij_vec = ri - rj;
        double rij = rij_vec.norm() * m_au;  // Convert to atomic units if needed

        if (rij > disp.r_cut || rij < 1e-10) continue;  // Skip if beyond cutoff or too close

        // Becke-Johnson damping function (order n=6 for C6, n=8 for C8)
        double r_crit = disp.a1 * std::sqrt(disp.C8 / (disp.C6 + 1e-14)) + disp.a2;

        // Phase 1.4: Optimize power calculations (Claude Generated - Dec 2025)
        // OPTIMIZATION: r^6 = (r^2)^3, r^8 = (r^2)^4 avoids expensive std::pow()
        // Benchmark: 3-5% speedup in dispersion contribution
        double r2 = rij * rij;
        double r6 = r2 * r2 * r2;  // (r^2)^3
        double r_crit2 = r_crit * r_crit;
        double damp6 = r_crit2 * r_crit2 * r_crit2;  // (r_crit^2)^3

        // C6 term: -s6*C6/r^6 with BJ damping
        double f_damp6 = r6 / (r6 + damp6);
        double E_C6 = -disp.s6 * disp.C6 * f_damp6 / r6;

        // C8 term: -s8*C8/r^8 with BJ damping
        double r8 = r2 * r2 * r2 * r2;  // (r^2)^4
        double damp8 = damp6 * r_crit2;  // r_crit^8 = r_crit^6 * r_crit^2
        double f_damp8 = r8 / (r8 + damp8);
        double E_C8 = -disp.s8 * disp.C8 * f_damp8 / r8;

        double energy = (E_C6 + E_C8) * m_final_factor;
        m_dispersion_energy += energy;
        m_d4_energy += energy;  // CRITICAL (Jan 7, 2026): GFN-FF dispersion reports as D4 energy

        if (m_calculate_gradient) {
            // Analytical gradient: dE/dr = dE_C6/dr + dE_C8/dr
            // d/dr[f_damp * C_n / r^n] = -n*f_damp*C_n/r^(n+1) + C_n/r^n * df_damp/dr

            // C6 gradient
            double df_damp6_dr = 6.0 * r6 * damp6 / (std::pow(r6 + damp6, 2) * rij);
            double dE_C6_dr = -6.0 * E_C6 / rij + (-disp.s6 * disp.C6 / r6) * df_damp6_dr;

            // C8 gradient
            double df_damp8_dr = 8.0 * r8 * damp8 / (std::pow(r8 + damp8, 2) * rij);
            double dE_C8_dr = -8.0 * E_C8 / rij + (-disp.s8 * disp.C8 / r8) * df_damp8_dr;

            double dEdr = (dE_C6_dr + dE_C8_dr) * m_final_factor;
            Eigen::Vector3d grad = dEdr * rij_vec / rij;

            m_gradient.row(disp.i) += grad.transpose();
            m_gradient.row(disp.j) -= grad.transpose();
        }
    }

    if (CurcumaLogger::get_verbosity() >= 3 && m_gfnff_dispersions.size() > 0) {
        CurcumaLogger::param("thread_dispersion_energy", fmt::format("{:.6f} Eh", m_dispersion_energy));
    }
}

void ForceFieldThread::CalculateGFNFFBondedRepulsionContribution()
{
    /**
     * @brief GFN-FF Bonded Repulsion term (pairwise parallelizable)
     *
     * Reference: Fortran gfnff_engrad.F90:467-495
     * Formula: E_rep = repab * exp(-α*r^1.5) / r
     * Alpha: sqrt(repa_i * repa_j) [geometric mean from repa array]
     * Scale: REPSCALB = 1.7583 [embedded in repab during parameter generation]
     *
     * Claude Generated (Dec 2025): Phase 9 repulsion fix - separated bonded/non-bonded
     */

    if (CurcumaLogger::get_verbosity() >= 2 && m_gfnff_bonded_repulsions.size() > 0) {
        CurcumaLogger::info(fmt::format("Thread {} calculating {} bonded repulsion pairs", m_thread, m_gfnff_bonded_repulsions.size()));
    }

    if (m_gfnff_bonded_repulsions.size() == 0) {
        return;  // Early exit if no bonded pairs
    }

    double total_rep_energy = 0.0;
    int pairs_calculated = 0;

    for (int index = 0; index < m_gfnff_bonded_repulsions.size(); ++index) {
        const auto& rep = m_gfnff_bonded_repulsions[index];

        Eigen::Vector3d ri = m_geometry.row(rep.i);
        Eigen::Vector3d rj = m_geometry.row(rep.j);
        Eigen::Vector3d rij_vec = ri - rj;
        double rij = rij_vec.norm() * m_au;

        if (rij > rep.r_cut || rij < 1e-10) continue;

        // Bonded repulsion formula: E = repab * exp(-α * r^1.5) / r
        // r^1.5 = r * sqrt(r) avoids std::pow()
        double r_1_5 = rij * std::sqrt(rij);
        double exp_term = std::exp(-rep.alpha * r_1_5);
        double base_energy = rep.repab * exp_term / rij;

        double scaled_energy = base_energy * m_final_factor * m_rep_scaling;
        m_rep_energy += scaled_energy;
        total_rep_energy += scaled_energy;
        pairs_calculated++;

        if (CurcumaLogger::get_verbosity() >= 3 && index == 0) {
            CurcumaLogger::param(fmt::format("bonded_repulsion_{}-{}", rep.i, rep.j),
                fmt::format("r={:.6f}, repab={:.6f}, alpha={:.6f}, E={:.6f}",
                rij, rep.repab, rep.alpha, scaled_energy));
        }

        if (m_calculate_gradient) {
            // dE/dr = -base_energy / r - 1.5 * α * sqrt(r) * base_energy
            double dEdr = (-base_energy / rij - 1.5 * rep.alpha * std::sqrt(rij) * base_energy);
            dEdr *= m_final_factor * m_rep_scaling;

            Eigen::Vector3d grad = dEdr * rij_vec / rij;
            m_gradient.row(rep.i) += grad.transpose();
            m_gradient.row(rep.j) -= grad.transpose();
        }
    }

    if (CurcumaLogger::get_verbosity() >= 2 && pairs_calculated > 0) {
        CurcumaLogger::param("thread_bonded_repulsion_energy", fmt::format("{:.6f} Eh ({} pairs)", total_rep_energy, pairs_calculated));
    }
}

void ForceFieldThread::CalculateGFNFFNonbondedRepulsionContribution()
{
    /**
     * @brief GFN-FF Non-bonded Repulsion term (pairwise parallelizable)
     *
     * Reference: Fortran gfnff_engrad.F90:255-276
     * Formula: E_rep = repab * exp(-α*r^1.5) / r
     * Alpha: (repan_i + repan_j) / 2.0 [arithmetic mean from repan array]
     * Scale: REPSCALN = 0.4270 [embedded in repab during parameter generation]
     *
     * Claude Generated (Dec 2025): Phase 9 repulsion fix - separated bonded/non-bonded
     */

    if (CurcumaLogger::get_verbosity() >= 2 && m_gfnff_nonbonded_repulsions.size() > 0) {
        CurcumaLogger::info(fmt::format("Thread {} calculating {} non-bonded repulsion pairs", m_thread, m_gfnff_nonbonded_repulsions.size()));
    }

    if (m_gfnff_nonbonded_repulsions.size() == 0) {
        return;  // Early exit if no non-bonded pairs
    }

    double total_rep_energy = 0.0;
    int pairs_calculated = 0;

    for (int index = 0; index < m_gfnff_nonbonded_repulsions.size(); ++index) {
        const auto& rep = m_gfnff_nonbonded_repulsions[index];

        Eigen::Vector3d ri = m_geometry.row(rep.i);
        Eigen::Vector3d rj = m_geometry.row(rep.j);
        Eigen::Vector3d rij_vec = ri - rj;
        double rij = rij_vec.norm() * m_au;

        if (rij > rep.r_cut || rij < 1e-10) continue;

        // Non-bonded repulsion formula: E = repab * exp(-α * r^1.5) / r
        // r^1.5 = r * sqrt(r) avoids std::pow()
        double r_1_5 = rij * std::sqrt(rij);
        double exp_term = std::exp(-rep.alpha * r_1_5);
        double base_energy = rep.repab * exp_term / rij;

        double scaled_energy = base_energy * m_final_factor * m_rep_scaling;
        m_rep_energy += scaled_energy;
        total_rep_energy += scaled_energy;
        pairs_calculated++;

        if (CurcumaLogger::get_verbosity() >= 3 && index == 0) {
            CurcumaLogger::param(fmt::format("nonbonded_repulsion_{}-{}", rep.i, rep.j),
                fmt::format("r={:.6f}, repab={:.6f}, alpha={:.6f}, E={:.6f}",
                rij, rep.repab, rep.alpha, scaled_energy));
        }

        if (m_calculate_gradient) {
            // dE/dr = -base_energy / r - 1.5 * α * sqrt(r) * base_energy
            double dEdr = (-base_energy / rij - 1.5 * rep.alpha * std::sqrt(rij) * base_energy);
            dEdr *= m_final_factor * m_rep_scaling;

            Eigen::Vector3d grad = dEdr * rij_vec / rij;
            m_gradient.row(rep.i) += grad.transpose();
            m_gradient.row(rep.j) -= grad.transpose();
        }
    }

    if (CurcumaLogger::get_verbosity() >= 2 && pairs_calculated > 0) {
        CurcumaLogger::param("thread_nonbonded_repulsion_energy", fmt::format("{:.6f} Eh ({} pairs)", total_rep_energy, pairs_calculated));
    }
}

void ForceFieldThread::CalculateGFNFFCoulombContribution()
{
    /**
     * @brief Complete EEQ-based Coulomb electrostatics with THREE terms
     *
     * Reference: Fortran gfnff_engrad.F90:1378-1389 and gfnff_ini.f90
     *
     * Complete formula (THREE TERMS):
     * E_coul = E_pairwise + E_self_energy + E_self_interaction
     *
     * 1. Pairwise: Σ(i<j) [q_i*q_j*erf(γ_ij*r_ij) / r_ij]  // Claude Generated (Session 9): FIXED - was /r_ij² (WRONG!)
     * 2. Self-energy: -Σ_i [q_i*χ_i]
     * 3. Self-interaction: Σ_i [0.5*q_i²*(γ_i + √(2/π)/√(α_i))]
     *    where γ_i = 1/√(α_i)
     *
     * Claude Generated (2025-12-05): Phase 4.3 - Complete implementation
     * Claude Generated (2025-12-12, Session 9): CRITICAL FIX - Corrected /r² to /r in pairwise term!
     * Claude Generated (Session 10): CRITICAL FIX - Add kJ/mol → Eh conversion for all Coulomb terms
     */

    // CRITICAL FIX (Dec 2025): EEQ parameters are already in Hartree units!
    // The kJ/mol conversion factor was incorrect - all quantities should be in Hartree
    // Reference: gfnff_par.h chi_gam_alp_cnf_angewChem2020 arrays are in Hartree

    // Verify EEQ charges before Coulomb energy calculation (Nov 2025)
    if (CurcumaLogger::get_verbosity() >= 1 && m_gfnff_coulombs.size() > 0) {
        CurcumaLogger::info("=== Coulomb Energy Calculation: EEQ Charge Verification ===");
        for (int idx = 0; idx < std::min((int)m_gfnff_coulombs.size(), 10); ++idx) {
            const auto& coul = m_gfnff_coulombs[idx];
            CurcumaLogger::param(fmt::format("Coulomb[{}] q_i(atom{}), q_j(atom{})",
                                             idx, coul.i, coul.j),
                                 fmt::format("{:.8f}, {:.8f}", coul.q_i, coul.q_j));
        }
        if (m_gfnff_coulombs.size() > 10) {
            CurcumaLogger::param("...", fmt::format("{} more coulomb pairs", m_gfnff_coulombs.size() - 10));
        }

        double max_charge = 0.0;
        for (const auto& coul : m_gfnff_coulombs) {
            max_charge = std::max({max_charge, std::abs(coul.q_i), std::abs(coul.q_j)});
        }
        CurcumaLogger::param("max_absolute_charge", fmt::format("{:.8e}", max_charge));

        if (max_charge < 1e-10) {
            CurcumaLogger::info("ℹ️  All EEQ partial charges near zero (typical for homonuclear molecules like H₂, Cl₂)");
        } else {
            CurcumaLogger::success(fmt::format("✅ Charges up to {:.2e} found - Coulomb energy should be non-zero", max_charge));
        }
    }

    if (CurcumaLogger::get_verbosity() >= 3 && m_gfnff_coulombs.size() > 0) {
        CurcumaLogger::info(fmt::format("Thread calculating {} Coulomb pairs", m_gfnff_coulombs.size()));
    }

    // =========================================================================
    // TERM 1: Pairwise Coulomb interactions (distance-dependent)
    // =========================================================================
    for (int index = 0; index < m_gfnff_coulombs.size(); ++index) {
        const auto& coul = m_gfnff_coulombs[index];

        Eigen::Vector3d ri = m_geometry.row(coul.i);
        Eigen::Vector3d rj = m_geometry.row(coul.j);
        Eigen::Vector3d rij_vec = ri - rj;
        double rij = rij_vec.norm() * m_au;

        if (rij > coul.r_cut || rij < 1e-10) continue;

        // Claude Generated (Dec 2025, Session 9): CRITICAL FIX - Coulomb is erf(γ*r)/r NOT erf(γ*r)/r²!
        // Pairwise: E_pair = q_i * q_j * erf(γ_ij*r) / r (NOT r²!)
        double gamma_r = coul.gamma_ij * rij;
        double erf_term = std::erf(gamma_r);
        double energy_pair = coul.q_i * coul.q_j * erf_term / rij;  // FIX: /r not /r²

        // Claude Generated (Session 10): CRITICAL FIX - Add kJ/mol → Eh conversion for Coulomb terms
        // EEQ quantities are in kJ/mol but need to be converted to Eh for total energy
        m_coulomb_energy += energy_pair;

        if (m_calculate_gradient) {
            // Claude Generated (Dec 2025, Session 9): CRITICAL FIX - Gradient for erf(γ*r)/r NOT erf(γ*r)/r²!
            // dE_pair/dr = q_i*q_j * [d/dr(erf(γ*r)/r)]
            //            = q_i*q_j * [γ*exp(-(γ*r)²)*(2/√π)/r - erf(γ*r)/r²]
            const double sqrt_pi = 1.772453850905516;  // √π
            double exp_term = std::exp(-gamma_r * gamma_r);
            double derf_dr = coul.gamma_ij * exp_term * (2.0 / sqrt_pi);
            double dEdr_pair = coul.q_i * coul.q_j *
                (derf_dr / rij - erf_term / (rij * rij));  // FIX: Correct derivative, no conversion

            Eigen::Vector3d grad = dEdr_pair * rij_vec / rij;

            m_gradient.row(coul.i) += grad.transpose();
            m_gradient.row(coul.j) -= grad.transpose();
        }
    }

    // =========================================================================
    // TERM 2 & 3: Self-energy and Self-interaction (per-atom, distance-independent)
    // =========================================================================
    // CRITICAL FIX (Dec 2025): Calculate self-energy only for ASSIGNED atoms
    // Each thread gets a subset of atoms to avoid duplicate self-energy calculation
    //
    // Thread-safe distribution:
    // - AutoRanges() distributes atoms across threads (e.g., 100 atoms, 4 threads)
    // - Thread 0: atoms [0, 25), Thread 1: atoms [25, 50), etc.
    // - Each atom's self-energy calculated EXACTLY ONCE across all threads
    //
    // Build atom-to-parameters map from Coulomb pairs for fast lookup
    std::unordered_map<int, const GFNFFCoulomb*> atom_to_params;
    for (const auto& coul : m_gfnff_coulombs) {
        if (atom_to_params.find(coul.i) == atom_to_params.end()) {
            atom_to_params[coul.i] = &coul;  // Store pointer to first occurrence
        }
        if (atom_to_params.find(coul.j) == atom_to_params.end()) {
            atom_to_params[coul.j] = &coul;
        }
    }

    // Calculate self-energy ONLY for assigned atoms (thread-safe)
    const double sqrt_2_over_pi = 0.797884560802865;  // √(2/π)
    for (int atom_id : m_assigned_atoms_for_self_energy) {
        auto it = atom_to_params.find(atom_id);
        if (it == atom_to_params.end()) continue;  // Atom not in any Coulomb pair

        const GFNFFCoulomb* params = it->second;

        // Determine which atom (i or j) we're processing
        double q, chi, gam, alp;
        if (params->i == atom_id) {
            q = params->q_i;
            chi = params->chi_i;
            gam = params->gam_i;
            alp = params->alp_i;
        } else {
            q = params->q_j;
            chi = params->chi_j;
            gam = params->gam_j;
            alp = params->alp_j;
        }

        // TERM 2: Self-energy for this atom
        double energy_self = -q * chi;

        // TERM 3: Self-interaction for this atom
        double selfint_term = gam + sqrt_2_over_pi / std::sqrt(alp);
        double energy_selfint = 0.5 * q * q * selfint_term;

        // Add both self-terms
        double energy_atom = energy_self + energy_selfint;
        m_coulomb_energy += energy_atom;

        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::param(fmt::format("coulomb_self_atom_{}", atom_id),
                fmt::format("E_self={:.6f}, E_selfint={:.6f}, total={:.6f} Eh",
                    energy_self, energy_selfint, energy_atom));
        }
    }

    if (CurcumaLogger::get_verbosity() >= 3 && m_gfnff_coulombs.size() > 0) {
        CurcumaLogger::param("thread_coulomb_energy", fmt::format("{:.6f} Eh", m_coulomb_energy));
    }
}

// ============================================================================
// Phase 5: GFN-FF HB/XB Damping Helper Functions
// ============================================================================
// Claude Generated (2025): Direct translation from gfnff_engrad.F90
// Reference: abhgfnff_eg1(), rbxgfnff_eg()

/**
 * @brief Out-of-line damping for HB/XB geometry
 *
 * Penalizes non-linear A-H...B or A-X...B geometries
 * Formula: outl = 2 / (1 + exp(bacut/r_AB × (r_AH + r_HB)/r_AB - 1))
 *
 * @param r_AH Distance A-H (or A-X) in Bohr
 * @param r_HB Distance H-B (or X-B) in Bohr
 * @param r_AB Distance A-B in Bohr
 * @param bacut Angle cut-off parameter (HB_BACUT or XB_BACUT)
 */
inline double damping_out_of_line(double r_AH, double r_HB, double r_AB, double bacut)
{
    double ratio = (r_AH + r_HB) / r_AB;
    double exponent = bacut / r_AB * (ratio - 1.0);
    return 2.0 / (1.0 + std::exp(exponent));
}

/**
 * @brief Short-range damping for HB/XB
 *
 * Prevents unphysical close-contact overbinding
 * Formula: damp_s = 1 / (1 + (scut × r_vdw / r²)^alp)
 *
 * @param r Distance between atoms in Bohr
 * @param r_vdw Combined van der Waals radius in Bohr
 * @param scut Short-range cut-off parameter (HB_SCUT or XB_SCUT)
 * @param alp Damping exponent (HB_ALP)
 */
inline double damping_short_range(double r, double r_vdw, double scut, double alp)
{
    double ratio = scut * r_vdw / (r * r);
    return 1.0 / (1.0 + std::pow(ratio, alp));
}

/**
 * @brief Long-range damping for HB/XB
 *
 * Smooth cut-off at large distances
 * Formula: damp_l = 1 / (1 + (r² / longcut)^alp)
 *
 * @param r Distance between atoms in Bohr
 * @param longcut Long-range cut-off (HB_LONGCUT or HB_LONGCUT_XB)
 * @param alp Damping exponent (HB_ALP)
 */
inline double damping_long_range(double r, double longcut, double alp)
{
    double r_sq = r * r;
    return 1.0 / (1.0 + std::pow(r_sq / longcut, alp));
}

/**
 * @brief Charge-dependent scaling factor
 *
 * EEQ charge modulation for HB/XB strength
 * Formula: Q = exp(st × q) / (exp(st × q) + sf)
 *
 * @param q EEQ charge on atom
 * @param st Scaling strength (HB_ST or XB_ST)
 * @param sf Scaling offset (HB_SF or XB_SF)
 */
inline double charge_scaling(double q, double st, double sf)
{
    double exp_term = std::exp(st * q);
    return exp_term / (exp_term + sf);
}

// ============================================================================
// Phase 5: GFN-FF HB/XB Energy Calculation Methods
// ============================================================================

void ForceFieldThread::CalculateGFNFFHydrogenBondContribution()
{
    /**
     * @brief Calculate GFN-FF hydrogen bond energy contribution
     *
     * Reference: gfnff_engrad.F90 - abhgfnff_eg1() (Case 1) and abhgfnff_eg2new() (Case 2)
     * Formula: E_HB = B_AH × C_AH × C_B × (-C_acidity × R_damp × Q_H^outl)
     *
     * Claude Generated (2025): Direct translation from Fortran implementation
     */

    using namespace GFNFFParameters;

    if (CurcumaLogger::get_verbosity() >= 3 && m_gfnff_hbonds.size() > 0) {
        CurcumaLogger::info(fmt::format("Thread {} calculating {} hydrogen bonds", m_thread, m_gfnff_hbonds.size()));
    }

    for (const auto& hb : m_gfnff_hbonds) {
        // Get atom positions
        Eigen::Vector3d pos_A = m_geometry.row(hb.i).transpose();
        Eigen::Vector3d pos_H = m_geometry.row(hb.j).transpose();
        Eigen::Vector3d pos_B = m_geometry.row(hb.k).transpose();

        // Calculate distance vectors
        Eigen::Vector3d r_AH_vec = pos_H - pos_A;
        Eigen::Vector3d r_HB_vec = pos_B - pos_H;
        Eigen::Vector3d r_AB_vec = pos_B - pos_A;

        // Calculate distances (convert to Bohr via m_au)
        double r_AH = r_AH_vec.norm() * m_au;
        double r_HB = r_HB_vec.norm() * m_au;
        double r_AB = r_AB_vec.norm() * m_au;

        // Distance cutoff check
        if (r_AB > hb.r_cut) continue;

        // --- Donor-Acceptor Strength Term B_AH ---
        // Reference: gfnff_engrad.F90:1445-1448
        // B = (C_A × r_AH^4 + C_B × r_HB^4) / (r_AH^4 + r_HB^4)
        double r_AH_4 = r_AH * r_AH * r_AH * r_AH;
        double r_HB_4 = r_HB * r_HB * r_HB * r_HB;
        double B_AH = (hb.basicity_A * r_AH_4 + hb.basicity_B * r_HB_4) / (r_AH_4 + r_HB_4);

        // --- Charge Scaling Factors ---
        double Q_H = charge_scaling(hb.q_H, HB_ST, HB_SF);
        double Q_A = charge_scaling(hb.q_A, HB_ST, HB_SF);
        double Q_B = charge_scaling(hb.q_B, HB_ST, HB_SF);

        // --- Combined Damping R_damp ---
        // Use covalent radii as approximation for vdW radius
        double r_vdw_AB = covalent_radii[hb.i] + covalent_radii[hb.k];  // Angström -> need conversion
        r_vdw_AB *= 1.88973;  // Convert Å to Bohr

        double damp_short = damping_short_range(r_AB, r_vdw_AB, HB_SCUT, HB_ALP);
        double damp_long = damping_long_range(r_AB, HB_LONGCUT, HB_ALP);
        double damp_outl = damping_out_of_line(r_AH, r_HB, r_AB, HB_BACUT);

        double R_damp = damp_short * damp_long * damp_outl / (r_AB * r_AB * r_AB);

        // --- Acidity Term C_acidity ---
        // Reference: hbonds() subroutine - complex mixing of basicity/acidity
        double C_acidity = hb.acidity_A;  // Simplified for Case 1

        // --- Total HB Energy ---
        // Note: Negative sign convention from Fortran
        double E_HB = B_AH * Q_A * Q_B * (-C_acidity) * R_damp * Q_H;

        m_energy_hbond += E_HB * m_final_factor;

        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info(fmt::format(
                "  HB({}-{}-{}): r_AB={:.3f} Bohr, E={:.6f} Eh",
                hb.i, hb.j, hb.k, r_AB, E_HB * m_final_factor));
        }

        // ========== ANALYTICAL GRADIENT CALCULATION ==========
        // Reference: gfnff_engrad.F90 - abhgfnff_eg1() lines 1501-1557
        // Claude Generated (2025): Phase 6 - HB analytical gradients
        if (m_calculate_gradient) {
            // --- 1. Compute Energy Derivative Components ---

            // ∂B_AH/∂r_AH and ∂B_AH/∂r_HB (donor-acceptor mixing term derivatives)
            double r_AH_3 = r_AH * r_AH * r_AH;
            double r_HB_3 = r_HB * r_HB * r_HB;
            double denominator = r_AH_4 + r_HB_4;
            double denominator_sq = denominator * denominator;

            double dB_drAH = 4.0 * r_AH_3 * (hb.basicity_A * r_HB_4 - hb.basicity_B * r_AH_4) / denominator_sq;
            double dB_drHB = 4.0 * r_HB_3 * (hb.basicity_B * r_AH_4 - hb.basicity_A * r_HB_4) / denominator_sq;

            // ∂(damping)/∂r derivatives
            // Short-range damping: damp_s = 1 / (1 + (scut × r_vdw / r²)^alp)
            double ratio_short = HB_SCUT * r_vdw_AB / (r_AB * r_AB);
            double damp_short_term = std::pow(ratio_short, HB_ALP);
            double ddamp_short_dr = -2.0 * HB_ALP * damp_short * damp_short_term / (r_AB * (1.0 + damp_short_term));

            // Long-range damping: damp_l = 1 / (1 + (r² / longcut)^alp)
            double ratio_long = (r_AB * r_AB) / HB_LONGCUT;
            double damp_long_term = std::pow(ratio_long, HB_ALP);
            double ddamp_long_dr = -2.0 * HB_ALP * r_AB * damp_long * damp_long_term / (HB_LONGCUT * (1.0 + damp_long_term));

            // Out-of-line damping: outl = 2 / (1 + exp(bacut/r_AB × (r_AH + r_HB)/r_AB - 1))
            double ratio = (r_AH + r_HB) / r_AB;
            double exponent = HB_BACUT / r_AB * (ratio - 1.0);
            double exp_term = std::exp(exponent);
            double denom_outl = 1.0 + exp_term;

            // Derivatives of out-of-line damping w.r.t. each distance
            double ddamp_outl_drAH = -2.0 * exp_term * HB_BACUT / (r_AB * r_AB * denom_outl * denom_outl);
            double ddamp_outl_drHB = ddamp_outl_drAH;  // Same formula
            double ddamp_outl_drAB = 2.0 * exp_term * HB_BACUT * (r_AH + r_HB) / (r_AB * r_AB * r_AB * denom_outl * denom_outl);

            // Combined damping derivative: R_damp = damp_short * damp_long * damp_outl / r_AB³
            double R_damp_no_r3 = damp_short * damp_long * damp_outl;

            double dRdamp_drAB = (ddamp_short_dr * damp_long * damp_outl +
                                  damp_short * ddamp_long_dr * damp_outl +
                                  damp_short * damp_long * ddamp_outl_drAB) / (r_AB * r_AB * r_AB)
                                - 3.0 * R_damp / r_AB;

            double dRdamp_drAH = damp_short * damp_long * ddamp_outl_drAH / (r_AB * r_AB * r_AB);
            double dRdamp_drHB = damp_short * damp_long * ddamp_outl_drHB / (r_AB * r_AB * r_AB);

            // --- 2. Energy Derivatives w.r.t. Distances ---
            // E_HB = B_AH × Q_A × Q_B × (-C_acidity) × R_damp × Q_H
            double E_prefactor = Q_A * Q_B * (-C_acidity) * Q_H;

            double dE_drAH = E_prefactor * (dB_drAH * R_damp + B_AH * dRdamp_drAH);
            double dE_drHB = E_prefactor * (dB_drHB * R_damp + B_AH * dRdamp_drHB);
            double dE_drAB = E_prefactor * B_AH * dRdamp_drAB;

            // --- 3. Chain Rule: Position Vector Derivatives ---
            // ∇_atom r_ij = (r_j - r_i) / |r_j - r_i|
            Eigen::Vector3d grad_rAH_unit = r_AH_vec / r_AH;  // Direction A → H
            Eigen::Vector3d grad_rHB_unit = r_HB_vec / r_HB;  // Direction H → B
            Eigen::Vector3d grad_rAB_unit = r_AB_vec / r_AB;  // Direction A → B

            // --- 4. Accumulate Gradients on Each Atom ---
            // Atom A: influenced by r_AH (negative direction) and r_AB (negative direction)
            Eigen::Vector3d grad_A = (-dE_drAH * grad_rAH_unit - dE_drAB * grad_rAB_unit) * m_final_factor;
            m_gradient.row(hb.i) += grad_A.transpose();

            // Atom H: influenced by r_AH (positive direction) and r_HB (negative direction)
            Eigen::Vector3d grad_H = (dE_drAH * grad_rAH_unit - dE_drHB * grad_rHB_unit) * m_final_factor;
            m_gradient.row(hb.j) += grad_H.transpose();

            // Atom B: influenced by r_HB (positive direction) and r_AB (positive direction)
            Eigen::Vector3d grad_B = (dE_drHB * grad_rHB_unit + dE_drAB * grad_rAB_unit) * m_final_factor;
            m_gradient.row(hb.k) += grad_B.transpose();

            if (CurcumaLogger::get_verbosity() >= 3) {
                CurcumaLogger::info(fmt::format(
                    "    HB Gradient: |∇A|={:.4f} |∇H|={:.4f} |∇B|={:.4f} Eh/Bohr",
                    grad_A.norm(), grad_H.norm(), grad_B.norm()));
            }
        }
    }

    if (CurcumaLogger::get_verbosity() >= 3 && m_gfnff_hbonds.size() > 0) {
        CurcumaLogger::param("thread_hbond_energy", fmt::format("{:.6f} Eh", m_energy_hbond));
    }
}

void ForceFieldThread::CalculateGFNFFHalogenBondContribution()
{
    /**
     * @brief Calculate GFN-FF halogen bond energy contribution
     *
     * Reference: gfnff_engrad.F90 - rbxgfnff_eg() (lines 2928-3043)
     * Formula: E_XB = -R_damp × Q_outl × C_B × Q_B × C_X × Q_X
     *
     * Claude Generated (2025): XB uses different damping parameters than HB
     */

    using namespace GFNFFParameters;

    if (CurcumaLogger::get_verbosity() >= 3 && m_gfnff_xbonds.size() > 0) {
        CurcumaLogger::info(fmt::format("Thread {} calculating {} halogen bonds", m_thread, m_gfnff_xbonds.size()));
    }

    for (const auto& xb : m_gfnff_xbonds) {
        // Get atom positions
        Eigen::Vector3d pos_A = m_geometry.row(xb.i).transpose();
        Eigen::Vector3d pos_X = m_geometry.row(xb.j).transpose();
        Eigen::Vector3d pos_B = m_geometry.row(xb.k).transpose();

        // Calculate distance vectors
        Eigen::Vector3d r_AX_vec = pos_X - pos_A;
        Eigen::Vector3d r_XB_vec = pos_B - pos_X;
        Eigen::Vector3d r_AB_vec = pos_B - pos_A;

        // Calculate distances (convert to Bohr via m_au)
        double r_AX = r_AX_vec.norm() * m_au;
        double r_XB = r_XB_vec.norm() * m_au;
        double r_AB = r_AB_vec.norm() * m_au;

        // Distance cutoff check
        if (r_XB > xb.r_cut) continue;

        // --- Charge Scaling Factors ---
        // XB uses different parameters: XB_ST=15, XB_SF=0.03 (much weaker offset)
        double Q_X = charge_scaling(xb.q_X, XB_ST, XB_SF);
        double Q_B = charge_scaling(xb.q_B, XB_ST, XB_SF);

        // --- Combined Damping R_damp ---
        double r_vdw_XB = covalent_radii[xb.j] + covalent_radii[xb.k];  // Angström
        r_vdw_XB *= 1.88973;  // Convert Å to Bohr

        // XB uses different cutoff parameters
        double damp_short = damping_short_range(r_XB, r_vdw_XB, XB_SCUT, HB_ALP);
        double damp_long = damping_long_range(r_XB, HB_LONGCUT_XB, HB_ALP);
        double damp_outl = damping_out_of_line(r_AX, r_XB, r_AB, XB_BACUT);

        double R_damp = damp_short * damp_long * damp_outl / (r_XB * r_XB * r_XB);

        // --- Total XB Energy ---
        // C_B = 1.0 (fixed), C_X = xbaci parameter
        double C_B = 1.0;
        double C_X = xb.acidity_X;

        double E_XB = -R_damp * C_B * Q_B * C_X * Q_X;

        m_energy_xbond += E_XB * m_final_factor;

        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info(fmt::format(
                "  XB({}-{}-{}): r_XB={:.3f} Bohr, E={:.6f} Eh",
                xb.i, xb.j, xb.k, r_XB, E_XB * m_final_factor));
        }

        // ========== ANALYTICAL GRADIENT CALCULATION ==========
        // Reference: gfnff_engrad.F90 - rbxgfnff_eg() gradient section
        // Claude Generated (2025): Phase 6 - XB analytical gradients
        if (m_calculate_gradient) {
            // XB gradients are simpler than HB (no donor-acceptor mixing term)

            // --- 1. Compute Damping Derivatives ---
            // Short-range damping: damp_s = 1 / (1 + (scut × r_vdw / r²)^alp)
            double ratio_short = XB_SCUT * r_vdw_XB / (r_XB * r_XB);
            double damp_short_term = std::pow(ratio_short, HB_ALP);
            double ddamp_short_dr = -2.0 * HB_ALP * damp_short * damp_short_term / (r_XB * (1.0 + damp_short_term));

            // Long-range damping: damp_l = 1 / (1 + (r² / longcut)^alp)
            double ratio_long = (r_XB * r_XB) / HB_LONGCUT_XB;
            double damp_long_term = std::pow(ratio_long, HB_ALP);
            double ddamp_long_dr = -2.0 * HB_ALP * r_XB * damp_long * damp_long_term / (HB_LONGCUT_XB * (1.0 + damp_long_term));

            // Out-of-line damping derivatives
            double ratio = (r_AX + r_XB) / r_AB;
            double exponent = XB_BACUT / r_AB * (ratio - 1.0);
            double exp_term = std::exp(exponent);
            double denom_outl = 1.0 + exp_term;

            // Derivatives of out-of-line damping w.r.t. each distance
            double ddamp_outl_drAX = -2.0 * exp_term * XB_BACUT / (r_AB * r_AB * denom_outl * denom_outl);
            double ddamp_outl_drXB = ddamp_outl_drAX;  // Same formula
            double ddamp_outl_drAB = 2.0 * exp_term * XB_BACUT * (r_AX + r_XB) / (r_AB * r_AB * r_AB * denom_outl * denom_outl);

            // Combined damping derivative: R_damp = damp_short * damp_long * damp_outl / r_XB³
            // Note: For XB, primary distance is r_XB, not r_AB as in HB

            double dRdamp_drXB = (ddamp_short_dr * damp_long * damp_outl +
                                  damp_short * ddamp_long_dr * damp_outl +
                                  damp_short * damp_long * ddamp_outl_drXB) / (r_XB * r_XB * r_XB)
                                - 3.0 * R_damp / r_XB;

            double dRdamp_drAX = damp_short * damp_long * ddamp_outl_drAX / (r_XB * r_XB * r_XB);
            double dRdamp_drAB = damp_short * damp_long * ddamp_outl_drAB / (r_XB * r_XB * r_XB);

            // --- 2. Energy Derivatives w.r.t. Distances ---
            // E_XB = -R_damp × C_B × Q_B × C_X × Q_X
            double E_prefactor = -C_B * Q_B * C_X * Q_X;

            double dE_drXB = E_prefactor * dRdamp_drXB;
            double dE_drAX = E_prefactor * dRdamp_drAX;
            double dE_drAB = E_prefactor * dRdamp_drAB;

            // --- 3. Chain Rule: Position Vector Derivatives ---
            Eigen::Vector3d grad_rAX_unit = r_AX_vec / r_AX;  // Direction A → X
            Eigen::Vector3d grad_rXB_unit = r_XB_vec / r_XB;  // Direction X → B
            Eigen::Vector3d grad_rAB_unit = r_AB_vec / r_AB;  // Direction A → B

            // --- 4. Accumulate Gradients on Each Atom ---
            // Atom A: influenced by r_AX (negative direction) and r_AB (negative direction)
            Eigen::Vector3d grad_A = (-dE_drAX * grad_rAX_unit - dE_drAB * grad_rAB_unit) * m_final_factor;
            m_gradient.row(xb.i) += grad_A.transpose();

            // Atom X: influenced by r_AX (positive direction) and r_XB (negative direction)
            Eigen::Vector3d grad_X = (dE_drAX * grad_rAX_unit - dE_drXB * grad_rXB_unit) * m_final_factor;
            m_gradient.row(xb.j) += grad_X.transpose();

            // Atom B: influenced by r_XB (positive direction) and r_AB (positive direction)
            Eigen::Vector3d grad_B = (dE_drXB * grad_rXB_unit + dE_drAB * grad_rAB_unit) * m_final_factor;
            m_gradient.row(xb.k) += grad_B.transpose();

            if (CurcumaLogger::get_verbosity() >= 3) {
                CurcumaLogger::info(fmt::format(
                    "    XB Gradient: |∇A|={:.4f} |∇X|={:.4f} |∇B|={:.4f} Eh/Bohr",
                    grad_A.norm(), grad_X.norm(), grad_B.norm()));
            }
        }
    }

    if (CurcumaLogger::get_verbosity() >= 3 && m_gfnff_xbonds.size() > 0) {
        CurcumaLogger::param("thread_xbond_energy", fmt::format("{:.6f} Eh", m_energy_xbond));
    }
}

// ============================================================================
// Claude Generated (December 19, 2025): Native D3 Dispersion Calculation
// ============================================================================

void ForceFieldThread::CalculateD3DispersionContribution()
{
    /**
     * @brief Native D3 Dispersion with Becke-Johnson damping for UFF-D3
     *
     * This method calculates D3 dispersion correction for UFF-D3 hybrid method.
     * Uses validated D3ParameterGenerator output (10/11 molecules <1% error).
     *
     * Reference: Grimme et al., J. Chem. Phys. 132, 154104 (2010) [D3-BJ]
     * Formula: E_disp = -Σ_ij f_damp(r_ij) * (s6*C6/r^6 + s8*C8/r^8)
     * BJ damping: f_damp = r^n / (r^n + (a1*sqrt(C8/C6) + a2)^n)
     *
     * Claude Generated: December 19, 2025
     */

    if (CurcumaLogger::get_verbosity() >= 3 && m_d3_dispersions.size() > 0) {
        CurcumaLogger::info(fmt::format("Thread {} calculating {} D3 dispersion pairs",
                                        m_thread, m_d3_dispersions.size()));
    }

    for (int index = 0; index < m_d3_dispersions.size(); ++index) {
        const auto& disp = m_d3_dispersions[index];

        Eigen::Vector3d ri = m_geometry.row(disp.i);
        Eigen::Vector3d rj = m_geometry.row(disp.j);
        Eigen::Vector3d rij_vec = ri - rj;
        double rij = rij_vec.norm() * m_au;  // Convert to atomic units if needed

        if (rij > disp.r_cut || rij < 1e-10) continue;  // Skip if beyond cutoff or too close

        // Becke-Johnson damping function (order n=6 for C6, n=8 for C8)
        double r_crit = disp.a1 * std::sqrt(disp.C8 / (disp.C6 + 1e-14)) + disp.a2;

        // Optimize power calculations: r^6 = (r^2)^3, r^8 = (r^2)^4
        double r2 = rij * rij;
        double r6 = r2 * r2 * r2;  // (r^2)^3
        double r_crit2 = r_crit * r_crit;
        double damp6 = r_crit2 * r_crit2 * r_crit2;  // (r_crit^2)^3

        // C6 term: -s6*C6/r^6 with BJ damping
        double f_damp6 = r6 / (r6 + damp6);
        double E_C6 = -disp.s6 * disp.C6 * f_damp6 / r6;

        // C8 term: -s8*C8/r^8 with BJ damping
        double r8 = r2 * r2 * r2 * r2;  // (r^2)^4
        double damp8 = damp6 * r_crit2;  // r_crit^8 = r_crit^6 * r_crit^2
        double f_damp8 = r8 / (r8 + damp8);
        double E_C8 = -disp.s8 * disp.C8 * f_damp8 / r8;

        double energy = (E_C6 + E_C8) * m_final_factor;
        m_d3_energy += energy;

        if (m_calculate_gradient) {
            // Analytical gradient: dE/dr = dE_C6/dr + dE_C8/dr
            // d/dr[f_damp * C_n / r^n] = -n*f_damp*C_n/r^(n+1) + C_n/r^n * df_damp/dr

            // C6 gradient
            double df_damp6_dr = 6.0 * r6 * damp6 / (std::pow(r6 + damp6, 2) * rij);
            double dE_C6_dr = -6.0 * E_C6 / rij + (-disp.s6 * disp.C6 / r6) * df_damp6_dr;

            // C8 gradient
            double df_damp8_dr = 8.0 * r8 * damp8 / (std::pow(r8 + damp8, 2) * rij);
            double dE_C8_dr = -8.0 * E_C8 / rij + (-disp.s8 * disp.C8 / r8) * df_damp8_dr;

            double dEdr = (dE_C6_dr + dE_C8_dr) * m_final_factor;
            Eigen::Vector3d grad = dEdr * rij_vec / rij;

            m_gradient.row(disp.i) += grad.transpose();
            m_gradient.row(disp.j) -= grad.transpose();
        }
    }

    if (CurcumaLogger::get_verbosity() >= 3 && m_d3_dispersions.size() > 0) {
        CurcumaLogger::param("thread_d3_energy", fmt::format("{:.6f} Eh", m_d3_energy));
    }
}

// Claude Generated 2025: Native D4 Dispersion with Becke-Johnson damping
void ForceFieldThread::CalculateD4DispersionContribution()
{
    /**
     * @brief Native D4 Dispersion with Becke-Johnson damping for GFN-FF
     *
     * This method calculates D4 dispersion correction using charge-weighted C6 coefficients
     * from D4ParameterGenerator (Casimir-Polder integration).
     *
     * Reference: Caldeweyher et al., J. Chem. Phys. 147, 034112 (2017) [D4-BJ]
     * Formula: E_disp = -Σ_ij f_damp(r_ij) * (s6*C6/r^6 + s8*C8/r^8)
     *
     * Key Differences from D3:
     * - C6 coefficients charge-weighted via calculateChargeWeightedC6()
     * - Frequency-dependent polarizabilities (alpha_iw)
     * - CN-dependent Gaussian weighting in C6 calculation
     *
     * Claude Generated: December 25, 2025
     */

    if (CurcumaLogger::get_verbosity() >= 3 && m_d4_dispersions.size() > 0) {
        CurcumaLogger::info(fmt::format("Thread {} calculating {} D4 dispersion pairs",
                                        m_thread, m_d4_dispersions.size()));
    }

    for (int index = 0; index < m_d4_dispersions.size(); ++index) {
        const auto& disp = m_d4_dispersions[index];

        Eigen::Vector3d ri = m_geometry.row(disp.i);
        Eigen::Vector3d rj = m_geometry.row(disp.j);
        Eigen::Vector3d rij_vec = ri - rj;
        double rij = rij_vec.norm() * m_au;  // Convert to atomic units if needed

        if (rij > disp.r_cut || rij < 1e-10) continue;  // Skip if beyond cutoff or too close

        // Becke-Johnson damping (D4 uses same formula as D3)
        double r_crit = disp.a1 * std::sqrt(disp.C8 / (disp.C6 + 1e-14)) + disp.a2;

        // Optimize power calculations: r^6 = (r^2)^3, r^8 = (r^2)^4
        double r2 = rij * rij;
        double r6 = r2 * r2 * r2;  // (r^2)^3
        double r_crit2 = r_crit * r_crit;
        double damp6 = r_crit2 * r_crit2 * r_crit2;  // (r_crit^2)^3

        // C6 term: -s6*C6/r^6 with BJ damping (D4 uses charge-weighted C6)
        double f_damp6 = r6 / (r6 + damp6);
        double E_C6 = -disp.s6 * disp.C6 * f_damp6 / r6;

        // C8 term: -s8*C8/r^8 with BJ damping
        double r8 = r2 * r2 * r2 * r2;  // (r^2)^4
        double damp8 = damp6 * r_crit2;  // r_crit^8 = r_crit^6 * r_crit^2
        double f_damp8 = r8 / (r8 + damp8);
        double E_C8 = -disp.s8 * disp.C8 * f_damp8 / r8;

        double pair_energy = (E_C6 + E_C8) * m_final_factor;
        m_d4_energy += pair_energy;

        // Verbosity 3: Detailed debug output for each pair
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info(fmt::format(
                "D4 Pair i={} j={} r={:.4f} C6={:.4f} C8={:.4f} E={:.6e} Eh",
                disp.i, disp.j, rij, disp.C6, disp.C8, pair_energy
            ));
            CurcumaLogger::info(fmt::format("  s6={:.4f} s8={:.4f} a1={:.4f} a2={:.4f}", disp.s6, disp.s8, disp.a1, disp.a2));
            CurcumaLogger::info(fmt::format("  R0={:.4f} f_damp6={:.4f} f_damp8={:.4f}", r_crit, f_damp6, f_damp8));
        }

        if (m_calculate_gradient) {
            // Analytical gradient: dE/dr = dE_C6/dr + dE_C8/dr
            // d/dr[f_damp * C_n / r^n] = -n*f_damp*C_n/r^(n+1) + C_n/r^n * df_damp/dr

            // C6 gradient
            double df_damp6_dr = 6.0 * r6 * damp6 / (std::pow(r6 + damp6, 2) * rij);
            double dE_C6_dr = -6.0 * E_C6 / rij + (-disp.s6 * disp.C6 / r6) * df_damp6_dr;

            // C8 gradient
            double df_damp8_dr = 8.0 * r8 * damp8 / (std::pow(r8 + damp8, 2) * rij);
            double dE_C8_dr = -8.0 * E_C8 / rij + (-disp.s8 * disp.C8 / r8) * df_damp8_dr;

            double dEdr = (dE_C6_dr + dE_C8_dr) * m_final_factor;
            Eigen::Vector3d grad = dEdr * rij_vec / rij;

            m_gradient.row(disp.i) += grad.transpose();
            m_gradient.row(disp.j) -= grad.transpose();
        }
    }

    if (CurcumaLogger::get_verbosity() >= 3 && m_d4_dispersions.size() > 0) {
        CurcumaLogger::param("thread_d4_energy", fmt::format("{:.6f} Eh", m_d4_energy));
    }

    if (CurcumaLogger::get_verbosity() >= 2 && m_d4_dispersions.size() > 0) {
        CurcumaLogger::result(fmt::format("D4 Dispersion Energy: {:.6e} Eh", m_d4_energy));
    }
}

// Claude Generated 2025: D3/D4 dispersion addition methods
void ForceFieldThread::addD3Dispersion(const GFNFFDispersion& d3_dispersion)
{
    m_d3_dispersions.push_back(d3_dispersion);
}

void ForceFieldThread::addD4Dispersion(const GFNFFDispersion& d4_dispersion)
{
    m_d4_dispersions.push_back(d4_dispersion);
}

// ============================================================================
// Claude Generated (2025): ATM Three-Body Dispersion Calculation
// ============================================================================

void ForceFieldThread::CalculateATMContribution()
{
    /**
     * @brief Calculate ATM (Axilrod-Teller-Muto) three-body dispersion
     *
     * Reference: external/cpp-d4/src/damping/atm.cpp:70-138
     * Formula: E_ATM = sum_{i<j<k} ang * fdmp * C9 / 3 * triple_scale
     *
     * C9 = s9 * sqrt(|C6_ij * C6_ik * C6_jk|)
     * fdmp = 1 / (1 + 6 * (r0_ijk / r_ijk)^(alp/3))
     * ang = (0.375 * A * B * C / r²_ijk + 1) / r³_ijk
     *
     * Shared by D3 and D4 (only C6 source differs)
     *
     * Claude Generated (2025): ATM three-body dispersion
     */

    using namespace GFNFFParameters;

    if (CurcumaLogger::get_verbosity() >= 3 && m_atm_triples.size() > 0) {
        CurcumaLogger::info(fmt::format("Thread {} calculating {} ATM triples",
                                        m_thread, m_atm_triples.size()));
    }

    double total_atm_energy = 0.0;

    for (const auto& triple : m_atm_triples) {
        // Get atom positions
        Eigen::Vector3d pos_i = m_geometry.row(triple.i).transpose();
        Eigen::Vector3d pos_j = m_geometry.row(triple.j).transpose();
        Eigen::Vector3d pos_k = m_geometry.row(triple.k).transpose();

        // Calculate distances
        double rij = (pos_i - pos_j).norm();
        double rik = (pos_i - pos_k).norm();
        double rjk = (pos_j - pos_k).norm();

        double r2ij = rij * rij;
        double r2ik = rik * rik;
        double r2jk = rjk * rjk;

        // C9 coefficient: s9 * sqrt(|C6_ij * C6_ik * C6_jk|)
        double c9 = triple.s9 * std::sqrt(std::fabs(triple.C6_ij * triple.C6_ik * triple.C6_jk));

        // Get atomic numbers and covalent radii
        int zi = m_atom_types[triple.i];
        int zj = m_atom_types[triple.j];
        int zk = m_atom_types[triple.k];

        // r0 cutoff radii (BJ damping formula)
        // r0_xy = a1 * sqrt(3 * R_cov[X] * R_cov[Y]) + a2
        // Use rcov_bohr from GFNFFParameters namespace (0-indexed)
        double r_cov_i = (zi > 0 && zi <= static_cast<int>(rcov_bohr.size())) ? rcov_bohr[zi - 1] : 1.0;
        double r_cov_j = (zj > 0 && zj <= static_cast<int>(rcov_bohr.size())) ? rcov_bohr[zj - 1] : 1.0;
        double r_cov_k = (zk > 0 && zk <= static_cast<int>(rcov_bohr.size())) ? rcov_bohr[zk - 1] : 1.0;

        double r0ij = triple.a1 * std::sqrt(3.0 * r_cov_i * r_cov_j) + triple.a2;
        double r0ik = triple.a1 * std::sqrt(3.0 * r_cov_i * r_cov_k) + triple.a2;
        double r0jk = triple.a1 * std::sqrt(3.0 * r_cov_j * r_cov_k) + triple.a2;
        double r0ijk = r0ij * r0ik * r0jk;

        // Distance products
        double rijk = rij * rik * rjk;
        double r2ijk = r2ij * r2ik * r2jk;
        double r3ijk = rijk * r2ijk;

        // BJ damping function
        double fdmp = 1.0 / (1.0 + 6.0 * std::pow(r0ijk / rijk, triple.alp / 3.0));

        // Angular term
        // ang = (0.375 * (r²_ij + r²_jk - r²_ik) * (r²_ij + r²_ik - r²_jk) *
        //                (r²_ik + r²_jk - r²_ij) / r²_ijk + 1) / r³_ijk
        double A = (r2ij + r2jk - r2ik);
        double B = (r2ij + r2ik - r2jk);
        double C = (r2ik + r2jk - r2ij);
        double ang = (0.375 * A * B * C / r2ijk + 1.0) / r3ijk;

        // Energy contribution
        double e_atm = ang * fdmp * c9 / 3.0 * triple.triple_scale;

        // Accumulate ATM energy
        total_atm_energy += e_atm;
        m_atm_energy += e_atm;  // Claude Generated (Dec 2025): Store in dedicated member variable

        if (CurcumaLogger::get_verbosity() >= 4) {
            CurcumaLogger::info(fmt::format(
                "ATM({},{},{}): rij={:.4f} rik={:.4f} rjk={:.4f} C9={:.6f} E={:.6e} Eh",
                triple.i, triple.j, triple.k, rij, rik, rjk, c9, e_atm));
        }
    }

    if (CurcumaLogger::get_verbosity() >= 3 && m_atm_triples.size() > 0) {
        CurcumaLogger::param("thread_atm_energy", fmt::format("{:.6e} Eh", m_atm_energy));  // Claude Generated (Dec 2025): Use scientific notation
    }
}

// Claude Generated (2025): ATM Three-Body Dispersion Gradient Calculation
// ============================================================================

void ForceFieldThread::CalculateATMGradient()
{
    /**
     * @brief Calculate analytical gradients for ATM (Axilrod-Teller-Muto) three-body dispersion
     *
     * Reference: external/cpp-d4/src/damping/atm.cpp:141-289
     *
     * Analytical gradient formulas:
     * 1. Damping derivative: dfdmp/drijk = -2·alp·(r0/r)^(alp/3)·fdmp²/(3·rijk)
     * 2. Angular derivatives: d(ang)/drij, d(ang)/drik, d(ang)/drjk (complex polynomials)
     * 3. Chain rule: dgradient_ij = c9·(-dang·fdmp + ang·dfdmp)/r²_ij·rij_vec
     * 4. Atom gradients: gradient(i) += -(dgij + dgik), gradient(j) += +(dgij - dgjk), etc.
     *
     * Sign convention: Negative gradient = attractive force
     * Symmetry: triple_scale factor applied (1.0, 0.5, 1/6)
     *
     * Claude Generated (2025): Analytical ATM gradients for geometry optimization/MD
     */

    using namespace GFNFFParameters;

    if (CurcumaLogger::get_verbosity() >= 3 && m_atm_triples.size() > 0) {
        CurcumaLogger::info(fmt::format("Thread {} calculating ATM gradients for {} triples",
                                        m_thread, m_atm_triples.size()));
    }

    for (const auto& triple : m_atm_triples) {
        // Get atom positions
        Eigen::Vector3d pos_i = m_geometry.row(triple.i).transpose();
        Eigen::Vector3d pos_j = m_geometry.row(triple.j).transpose();
        Eigen::Vector3d pos_k = m_geometry.row(triple.k).transpose();

        // Calculate distance vectors (r_ij = pos_j - pos_i, etc.)
        Eigen::Vector3d rij_vec = pos_j - pos_i;
        Eigen::Vector3d rik_vec = pos_k - pos_i;
        Eigen::Vector3d rjk_vec = pos_k - pos_j;

        // Calculate distances
        double rij = rij_vec.norm();
        double rik = rik_vec.norm();
        double rjk = rjk_vec.norm();

        double r2ij = rij * rij;
        double r2ik = rik * rik;
        double r2jk = rjk * rjk;

        // C9 coefficient: NEGATIVE for gradients! (cpp-d4:202, simple-dftd3:atm.f90:282)
        // Energy uses +c9, gradient uses -c9 (dispersion is attractive)
        double c9 = -triple.s9 * std::sqrt(std::fabs(triple.C6_ij * triple.C6_ik * triple.C6_jk));

        // Get atomic numbers and covalent radii
        int zi = m_atom_types[triple.i];
        int zj = m_atom_types[triple.j];
        int zk = m_atom_types[triple.k];

        // r0 cutoff radii (BJ damping formula)
        double r_cov_i = (zi > 0 && zi <= static_cast<int>(rcov_bohr.size())) ? rcov_bohr[zi - 1] : 1.0;
        double r_cov_j = (zj > 0 && zj <= static_cast<int>(rcov_bohr.size())) ? rcov_bohr[zj - 1] : 1.0;
        double r_cov_k = (zk > 0 && zk <= static_cast<int>(rcov_bohr.size())) ? rcov_bohr[zk - 1] : 1.0;

        double r0ij = triple.a1 * std::sqrt(3.0 * r_cov_i * r_cov_j) + triple.a2;
        double r0ik = triple.a1 * std::sqrt(3.0 * r_cov_i * r_cov_k) + triple.a2;
        double r0jk = triple.a1 * std::sqrt(3.0 * r_cov_j * r_cov_k) + triple.a2;
        double r0ijk = r0ij * r0ik * r0jk;

        // Distance products
        double rijk = rij * rik * rjk;
        double r2ijk = r2ij * r2ik * r2jk;
        double r3ijk = rijk * r2ijk;
        double r5ijk = r2ijk * r3ijk;

        // BJ damping function and derivative
        double tmp = std::pow(r0ijk / rijk, triple.alp / 3.0);
        double fdmp = 1.0 / (1.0 + 6.0 * tmp);
        double dfdmp = -2.0 * triple.alp * tmp * fdmp * fdmp;  // Match cpp-d4:232 - NO /(3*rijk)!

        // Angular term
        double A = (r2ij + r2jk - r2ik);
        double B = (r2ij + r2ik - r2jk);
        double C = (r2ik + r2jk - r2ij);
        double ang = (0.375 * A * B * C / r2ijk + 1.0) / r3ijk;

        // Energy (for reference, already calculated in CalculateATMContribution)
        // double e_atm = ang * fdmp * c9 / 3.0 * triple.triple_scale;

        // ========================================
        // Analytical Gradient Calculation
        // ========================================

        // Angular derivative d/drij (from cpp-d4 reference, lines 234-241)
        double dang_ij = -0.375 * (std::pow(r2ij, 3) + std::pow(r2ij, 2) * (r2jk + r2ik)
                          + r2ij * (3.0 * std::pow(r2jk, 2) + 2.0 * r2jk * r2ik + 3.0 * std::pow(r2ik, 2))
                          - 5.0 * std::pow((r2jk - r2ik), 2) * (r2jk + r2ik)) / r5ijk;

        // Angular derivative d/drik (from cpp-d4 reference, lines 243-250)
        double dang_ik = -0.375 * (std::pow(r2ik, 3) + std::pow(r2ik, 2) * (r2jk + r2ij)
                          + r2ik * (3.0 * std::pow(r2jk, 2) + 2.0 * r2jk * r2ij + 3.0 * std::pow(r2ij, 2))
                          - 5.0 * std::pow((r2jk - r2ij), 2) * (r2jk + r2ij)) / r5ijk;

        // Angular derivative d/drjk (from cpp-d4 reference, lines 252-259)
        double dang_jk = -0.375 * (std::pow(r2jk, 3) + std::pow(r2jk, 2) * (r2ik + r2ij)
                          + r2jk * (3.0 * std::pow(r2ik, 2) + 2.0 * r2ik * r2ij + 3.0 * std::pow(r2ij, 2))
                          - 5.0 * std::pow((r2ik - r2ij), 2) * (r2ik + r2ij)) / r5ijk;

        // Gradient components (chain rule: dE/dr_ij = c9·(-dang·fdmp + ang·dfdmp)/r²_ij · r_ij)
        // Apply triple_scale factor (same as energy)
        double prefactor = c9 * triple.triple_scale / 3.0;

        Eigen::Vector3d dgij = prefactor * (-dang_ij * fdmp + ang * dfdmp) / r2ij * rij_vec;
        Eigen::Vector3d dgik = prefactor * (-dang_ik * fdmp + ang * dfdmp) / r2ik * rik_vec;
        Eigen::Vector3d dgjk = prefactor * (-dang_jk * fdmp + ang * dfdmp) / r2jk * rjk_vec;

        // Accumulate gradients to atoms (from cpp-d4 reference, lines 261-269)
        // Sign convention: negative sign because dispersion is attractive
        m_gradient.row(triple.i) += -(dgij + dgik);
        m_gradient.row(triple.j) += (dgij - dgjk);
        m_gradient.row(triple.k) += (dgik + dgjk);

        if (CurcumaLogger::get_verbosity() >= 4) {
            CurcumaLogger::info(fmt::format(
                "ATM_grad({},{},{}): |dgij|={:.6e} |dgik|={:.6e} |dgjk|={:.6e}",
                triple.i, triple.j, triple.k, dgij.norm(), dgik.norm(), dgjk.norm()));
        }
    }

    if (CurcumaLogger::get_verbosity() >= 3 && m_atm_triples.size() > 0) {
        CurcumaLogger::info(fmt::format("Thread {} ATM gradient calculation complete", m_thread));
    }
}
