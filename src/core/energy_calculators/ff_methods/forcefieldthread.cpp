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
#include "src/core/units.h"

// Claude Generated (Jan 17, 2026): Required for diagnostic logging in batm calculation
#include "src/core/curcuma_logger.h"
#include <fmt/core.h>
#include <fmt/format.h>

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
    // Set method-specific unit conversion factor (m_au)
    // CRITICAL: GFNFF wrapper already passes coordinates in Bohr.
    // Thus m_au (Distance multiplier) must be 1.0 for method 3.
    // For methods using Angstrom coordinates (UFF/QMDFF), m_au is ANGSTROM_TO_BOHR.
    if (m_method == 3 || m_method == 5) {
        m_au = 1.0;
    } else {
        m_au = 1.889726125; // Angstrom to Bohr for UFF-based non-bonded terms
    }

    m_angle_energy = 0.0;
    m_bond_energy = 0.0;
    m_vdw_energy = 0;
    m_rep_energy = 0;
    m_bonded_rep_energy = 0.0;
    m_nonbonded_rep_energy = 0.0;
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
    m_stors_energy = 0.0; // Claude Generated (March 2026): Reset triple bond torsion energy
    m_eq_energy = 0.0;  // Also reset EQ energy for consistency

    // Claude Generated 2025: Reset native D3/D4 dispersion energy terms
    m_d3_energy = 0.0;
    m_d4_energy = 0.0;
    m_atm_energy = 0.0;  // Claude Generated (Dec 2025): Reset ATM three-body dispersion
    m_batm_energy = 0.0;  // Claude Generated (Jan 17, 2026): Reset batm three-body energy - CRITICAL FIX

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

        // Claude Generated (February 2026): Wrap all energy term calculations with timing
        // When m_store_gradient_components is active, capture gradient delta per term

        // Helper lambda: run calculation and capture gradient delta into component matrix
        // This avoids modifying every inner Calculate*() method individually
        auto runWithGradCapture = [this](const std::string& name, auto calc_fn, Matrix& comp_grad) {
            if (m_store_gradient_components && m_calculate_gradient) {
                Matrix before = m_gradient;
                timeEnergyTerm(name, calc_fn);
                comp_grad += (m_gradient - before);
            } else {
                timeEnergyTerm(name, calc_fn);
            }
        };

        // Claude Generated (Feb 21, 2026): Compute HB coordination numbers before bond energy
        // Reference: Fortran gfnff_engrad.F90:361 - called before bond loop if sum(nr_hb)>0
        computeHBCoordinationNumbers();

        // GFN-FF bonded terms
        runWithGradCapture("bonds", [this]() { CalculateGFNFFBondContribution(); }, m_gradient_bond);
        runWithGradCapture("angles", [this]() { CalculateGFNFFAngleContribution(); }, m_gradient_angle);
        runWithGradCapture("torsions", [this]() { CalculateGFNFFDihedralContribution(); }, m_gradient_torsion);
        runWithGradCapture("extra_torsions", [this]() { CalculateGFNFFExtraTorsionContribution(); }, m_gradient_torsion);
        runWithGradCapture("inversions", [this]() { CalculateGFNFFInversionContribution(); }, m_gradient_torsion);

        // GFN-FF non-bonded pairwise parallelizable terms (Phase 4)
        if (m_dispersion_enabled) {
            runWithGradCapture("dispersion", [this]() { CalculateGFNFFDispersionContribution(); }, m_gradient_dispersion);
        }
        if (m_repulsion_enabled) {
            runWithGradCapture("bonded_repulsion", [this]() { CalculateGFNFFBondedRepulsionContribution(); }, m_gradient_repulsion);
            // Claude Generated (Mar 2026): Non-bonded repulsion goes to total gradient only, NOT g_rep
            // Reference: Fortran gfnff_engrad.F90:280-304 uses reduction(+:erep, g) without g_rep
            timeEnergyTerm("nonbonded_repulsion", [this]() { CalculateGFNFFNonbondedRepulsionContribution(); });
        }
        if (m_coulomb_enabled) {
            runWithGradCapture("coulomb", [this]() { CalculateGFNFFCoulombContribution(); }, m_gradient_coulomb);
        }

        // GFN-FF hydrogen bond and halogen bond terms (Phase 5)
        if (m_hbond_enabled) {
            runWithGradCapture("hydrogen_bonds", [this]() { CalculateGFNFFHydrogenBondContribution(); }, m_gradient_hb);
            runWithGradCapture("halogen_bonds", [this]() { CalculateGFNFFHalogenBondContribution(); }, m_gradient_xb);
        }

        // Claude Generated (December 19, 2025): Native D3/D4 dispersion calculation for GFN-FF
        if (m_d3_dispersions.size() > 0) {
            if (CurcumaLogger::get_verbosity() >= 3) {
                CurcumaLogger::info(fmt::format("Thread {} calculating {} D3 dispersion pairs", m_thread, m_d3_dispersions.size()));
            }
            runWithGradCapture("d3_dispersion", [this]() { CalculateD3DispersionContribution(); }, m_gradient_dispersion);
        }

        if (m_d4_dispersions.size() > 0) {
            if (CurcumaLogger::get_verbosity() >= 3) {
                CurcumaLogger::info(fmt::format("Thread {} calculating {} D4 dispersion pairs", m_thread, m_d4_dispersions.size()));
            }
            runWithGradCapture("d4_dispersion", [this]() { CalculateD4DispersionContribution(); }, m_gradient_dispersion);
        }

        // ATM three-body dispersion (D3/D4)
        if (!m_atm_triples.empty()) {
            if (CurcumaLogger::get_verbosity() >= 3) {
                CurcumaLogger::info(fmt::format("Thread {} calculating {} ATM triples", m_thread, m_atm_triples.size()));
            }
            // Claude Generated (Mar 2026): ATM captured in own component for structural
            // correctness — Fortran puts ATM outside g_disp (gfnff_gdisp0.f90:308-400)
            runWithGradCapture("atm_dispersion", [this]() {
                CalculateATMContribution();
                if (m_calculate_gradient) {
                    CalculateATMGradient();
                }
            }, m_gradient_atm);
        }

        // BF (Bonded ATM/GFN-FF) - Claude Generated (January 17, 2026)
        if (!m_gfnff_batms.empty()) {
            if (CurcumaLogger::get_verbosity() >= 3) {
                CurcumaLogger::info(fmt::format("Thread {} calculating {} batm triples", m_thread, m_gfnff_batms.size()));
            }
            runWithGradCapture("batm", [this]() { CalculateGFNFFBatmContribution(); }, m_gradient_batm);
        }

        if (!m_gfnff_storsions.empty()) {
            runWithGradCapture("storsions", [this]() { CalculateGFNFFSTorsionContribution(); }, m_gradient_torsion);
        }

        if (CurcumaLogger::get_verbosity() >= 1) {
            CurcumaLogger::info(fmt::format("Thread {}: Bond={:.6f}, Angle={:.6f}, Torsion={:.6f}, Rep={:.6f}, Coul={:.6f}, Disp={:.6f}, Batm={:.6f}",
                m_thread, m_bond_energy, m_angle_energy, m_dihedral_energy,
                m_rep_energy, m_coulomb_energy, m_dispersion_energy, m_batm_energy));
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

void ForceFieldThread::addGFNFFSTorsion(const GFNFFSTorsion& storsion)
{
    m_gfnff_storsions.push_back(storsion);
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

/**
 * @brief Set HB list (bulk update) for dynamic MD updates
 *
 * Claude Generated (Feb 15, 2026): Bulk setter for HB/XB dynamic updates
 * Reference: Fortran gfnff_engrad.F90:246-260 - rebuild HB/XB lists at each energy call
 *
 * @param hbonds Vector of hydrogen bond parameters (copy is efficient)
 */
void ForceFieldThread::setGFNFFHBonds(const std::vector<GFNFFHydrogenBond>& hbonds)
{
    m_gfnff_hbonds = hbonds;  // Efficient vector copy
}

/**
 * @brief Set XB list (bulk update) for dynamic MD updates
 *
 * Claude Generated (Feb 15, 2026): Bulk setter for XB dynamic updates
 *
 * @param xbonds Vector of halogen bond parameters
 */
void ForceFieldThread::setGFNFFHalogenBonds(const std::vector<GFNFFHalogenBond>& xbonds)
{
    m_gfnff_xbonds = xbonds;  // Efficient vector copy
}

void ForceFieldThread::addATMTriple(const ATMTriple& triple)
{
    m_atm_triples.push_back(triple);
}

// BF (Bonded ATM/GFN-FF) - Claude Generated (January 17, 2026)
void ForceFieldThread::addGFNFFBatmTriple(const GFNFFBatmTriple& batm_triple)
{
    m_gfnff_batms.push_back(batm_triple);
}

void ForceFieldThread::CalculateUFFBondContribution()
{
    double factor = m_final_factor * m_bond_scaling;

    for (int index = 0; index < m_uff_bonds.size(); ++index) {
        const auto& bond = m_uff_bonds[index];

        Vector i = geom().row(bond.i);
        Vector j = geom().row(bond.j);
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
        auto i = geom().row(angle.i);
        auto j = geom().row(angle.j);
        auto k = geom().row(angle.k);
        Matrix derivate;
        double costheta = UFF::AngleBending(i, j, k, derivate, m_calculate_gradient);

        // UFF angle bending potential
        double cos2theta = 2 * costheta * costheta - 1;
        m_angle_energy += angle.fc * (angle.C0 + angle.C1 * costheta + angle.C2 * cos2theta) * m_final_factor * m_angle_scaling;

        if (m_calculate_gradient) {
            // UFF gradient fix (Claude Generated Mar 2026):
            // AngleBending returns derivate = dθ/dr_i (gradient of angle θ, NOT cosθ).
            // diff must be dE/dθ, not dE/d(cosθ):
            //   dE/dθ = dE/d(cosθ) * d(cosθ)/dθ = fc*(C1+4*C2*cosθ) * (-sinθ)
            double sintheta = std::sin(std::acos(costheta));
            double dEdTheta = -sintheta * angle.fc * (angle.C1 + 4 * angle.C2 * costheta) * m_final_factor * m_angle_scaling;
            m_gradient.row(angle.i) += dEdTheta * derivate.row(0);
            m_gradient.row(angle.j) += dEdTheta * derivate.row(1);
            m_gradient.row(angle.k) += dEdTheta * derivate.row(2);
        }
    }
}

void ForceFieldThread::CalculateUFFDihedralContribution()
{
    for (int index = 0; index < m_uff_dihedrals.size(); ++index) {
        const auto& dihedral = m_uff_dihedrals[index];
        Eigen::Vector3d i = geom().row(dihedral.i);
        Eigen::Vector3d j = geom().row(dihedral.j);
        Eigen::Vector3d k = geom().row(dihedral.k);
        Eigen::Vector3d l = geom().row(dihedral.l);
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

        Eigen::Vector3d i = geom().row(inversion.i);
        Eigen::Vector3d j = geom().row(inversion.j);
        Eigen::Vector3d k = geom().row(inversion.k);
        Eigen::Vector3d l = geom().row(inversion.l);

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
        Eigen::Vector3d i = geom().row(vdw.i);
        Eigen::Vector3d j = geom().row(vdw.j);
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

        Vector i = geom().row(bond.i);
        Vector j = geom().row(bond.j);
        Vector ij = i - j;
        double distance = (ij).norm();

        // Matrix derivate;
        double fc = bond.r0_ij;
        const double ratio = bond.r0_ij / distance;
        // Claude Generated Fix (Jan 18, 2026): Correct bond order scaling exponent from 2.0 to 1.5
        // Issue: bond energy discrepancy in caffeine molecule (~0.6 Eh difference vs xtb reference)
        // Root cause: Incorrect exponent used in bond stretching energy formula
        // Reference: GFN-FF implementation should use 1.5 exponent for proper bond dissociation energy scaling
        m_bond_energy += fc * (1 + pow(ratio, bond.exponent) - 2 * pow(ratio, bond.exponent * 0.75));
        if (m_calculate_gradient) {

            double diff = 1 * fc * (-1 * bond.exponent * pow(ratio, bond.exponent - 1) + 2 * bond.exponent * 0.75 * pow(ratio, bond.exponent * 0.75 - 1));
            /*
                        Vector ijx = i+ dx - j;
                        double distancex1 = (ijx).norm();
                        double ratio = bond.r0_ij / distancex1;

                        // Claude Generated Fix (Jan 18, 2026): Updated commented code to match corrected exponent
                        double dx_p = fc * (1 + pow(ratio, bond.exponent) - 2 * pow(ratio, bond.exponent * 0.75));
                        ijx = i - dx - j;
                        distancex1 = (ijx).norm();
                        ratio = bond.r0_ij / distancex1;
                        double dx_m = fc * (1 + pow(ratio, bond.exponent) - 2 * pow(ratio, bond.exponent * 0.75));


                        Vector ijy = i+ dy - j;
                        double distancey1 = (ijy).norm();
                        ratio = bond.r0_ij / distancey1;

                        // Claude Generated Fix (Jan 18, 2026): Updated commented code to match corrected exponent
                        double dy_p = fc * (1 + pow(ratio, bond.exponent) - 2 * pow(ratio, bond.exponent * 0.75));
                        ijy = i - dy - j;
                        distancey1 = (ijy).norm();
                        ratio = bond.r0_ij / distancey1;
                        double dy_m = fc * (1 + pow(ratio, bond.exponent) - 2 * pow(ratio, bond.exponent * 0.75));


                        Vector ijz = i+ dz - j;
                        double distancez1 = (ijz).norm();
                        ratio = bond.r0_ij / distancez1;

                        // Claude Generated Fix (Jan 18, 2026): Updated commented code to match corrected exponent
                        double dz_p = fc * (1 + pow(ratio, bond.exponent) - 2 * pow(ratio, bond.exponent * 0.75));
                        ijz = i - dz - j;
                        distancez1 = (ijz).norm();
                        ratio = bond.r0_ij / distancez1;
                        double dz_m = fc * (1 + pow(ratio, bond.exponent) - 2 * pow(ratio, bond.exponent * 0.75));
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

        Vector i = geom().row(angle.i);
        Vector j = geom().row(angle.j);
        Vector k = geom().row(angle.k);
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
        Vector i = geom().row(eq.i);
        Vector j = geom().row(eq.j);
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
        geometry[i].x = geom()(i, 0) * m_au;
        geometry[i].y = geom()(i, 1) * m_au;
        geometry[i].z = geom()(i, 2) * m_au;
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

// Claude Generated (Feb 21, 2026): Compute HB coordination numbers for bond-HB coupling
// Reference: Fortran gfnff_engrad.F90:1069-1120 (dncoord_erf subroutine)
//
// For each A-H bond participating in hydrogen bridges, computes hb_cn_H by counting
// B atoms within erf-damped covalent radius distance from H. This CN is then used
// in egbond_hb to modulate the bond exponent: alpha_mod = (1 - 0.1*hb_cn_H) * alpha
//
// Parameters (from Fortran):
//   kn = 27.5 (error function steepness)
//   rcov_scal = 1.78 (covalent radius scaling)
//   thr = 900.0 Bohr² (distance threshold)
void ForceFieldThread::computeHBCoordinationNumbers()
{
    if (m_bond_hb_data.empty())
        return;

    // Use covalent radii from D3 (Bohr) with 4/3 scaling - matching Fortran param%rcov
    // Reference: gfnff_param.f90:381-405 — covalentRadD3 = raw_Å * aatoau * 4/3
    // Our covalent_rad_d3 in gfnff_par.h stores only raw_Å * aatoau (WITHOUT 4/3),
    // so we must apply the 4/3 factor here to match Fortran's param%rcov.
    // Claude Generated (Mar 4, 2026): Fix missing 4/3 scaling — caused hb_cn_H=0 for all bonds
    static const std::vector<double>& rcov_base = GFNFFParameters::covalent_rad_d3;
    constexpr double rcov_43 = 4.0 / 3.0;  // Fortran 4/3 scaling factor

    constexpr double kn = 27.5;
    constexpr double rcov_scal = 1.78;
    constexpr double thr = 900.0;  // Bohr² distance threshold

    // Build map from H atom index to accumulated hb_cn
    std::unordered_map<int, double> hb_cn_map;

    // Claude Generated (Feb 22, 2026): Clear and rebuild HB gradient entries
    // Reference: Fortran gfnff_engrad.F90:1054-1063 (hb_dcn chain-rule)
    m_hb_grad_entries.clear();

    constexpr double inv_sqrt_pi = 0.5641895835477563;  // 1/sqrt(pi)

    for (const auto& entry : m_bond_hb_data) {
        int H = entry.H;
        int ati = m_atom_types[H];  // Atomic number of H (should be 1)

        for (int B : entry.B_atoms) {
            int atj = m_atom_types[B];  // Atomic number of B

            // Distance H-B (both in Bohr for GFN-FF)
            double dx = geom()(B, 0) - geom()(H, 0);
            double dy = geom()(B, 1) - geom()(H, 1);
            double dz = geom()(B, 2) - geom()(H, 2);
            double r2 = dx * dx + dy * dy + dz * dz;

            if (r2 > thr) continue;
            double r = std::sqrt(r2);

            // rcov indices are 0-based (ati-1 for 1-based atomic number)
            // Apply 4/3 scaling to match Fortran param%rcov = covalentRadD3 * aatoau * 4/3
            double rcovij = rcov_scal * rcov_43 * (rcov_base[ati - 1] + rcov_base[atj - 1]);

            // erf-based coordination number contribution
            // Fortran: tmp = 0.5*(1 + erf(-kn*(r-rcovij)/rcovij))
            double arg = -kn * (r - rcovij) / rcovij;
            double tmp = 0.5 * (1.0 + std::erf(arg));

            hb_cn_map[H] += tmp;

            // Claude Generated (Feb 22, 2026): Gradient of erf-CN w.r.t. positions
            // d/dr [0.5*(1 + erf(-kn*(r-rcovij)/rcovij))]
            //   = inv_sqrt_pi * (-kn/rcovij) * exp(-arg²) / r * (r_B - r_H)
            // Reference: Fortran dncoord_erf derivative
            double dCN_dr = inv_sqrt_pi * (-kn / rcovij) * std::exp(-arg * arg) / r;
            Eigen::Vector3d r_HB(dx, dy, dz);  // B - H vector (in Bohr)
            Eigen::Vector3d dCN_dH = -dCN_dr * r_HB;
            Eigen::Vector3d dCN_dB = dCN_dr * r_HB;

            m_hb_grad_entries.push_back({H, B, dCN_dH, dCN_dB});
        }
    }

    // Update hb_cn_H on all bonds with nr_hb >= 1
    for (auto& bond : m_gfnff_bonds) {
        if (bond.nr_hb < 1) continue;

        // Identify which atom is H
        int H = -1;
        if (m_atom_types[bond.i] == 1) H = bond.i;
        else if (m_atom_types[bond.j] == 1) H = bond.j;

        if (H >= 0) {
            auto it = hb_cn_map.find(H);
            bond.hb_cn_H = (it != hb_cn_map.end()) ? it->second : 0.0;
        }
    }
}

void ForceFieldThread::CalculateGFNFFBondContribution()
{
    // Phase 1.3: Correct GFN-FF exponential bond potential
    // Formula from Fortran gfnff_engrad.F90:675-721
    // E_bond = k_b * exp(-α * (r - r₀)²)
    //
    // Claude Generated (Jan 18, 2026): DYNAMIC r0 CALCULATION
    // Reference: Fortran gfnff_engrad.F90:432-433, gfnff_rab.f90:147-153
    // r0 is recalculated at each energy evaluation using current D3 CN
    // Formula: r0 = (r0_base_i + cnfak_i*cn_i + r0_base_j + cnfak_j*cn_j + rabshift) * ff

    double factor = m_bond_scaling;

    // Check if we have D3 CN data for dynamic r0 calculation
    const bool have_cn_ptr = (m_d3_cn_ptr && m_d3_cn_ptr->size() > 0);
    const bool have_cn_local = (m_d3_cn.size() > 0);
    bool use_dynamic_r0 = have_cn_ptr || have_cn_local;
    const int cn_size = have_cn_ptr ? static_cast<int>(m_d3_cn_ptr->size())
                                    : static_cast<int>(m_d3_cn.size());

    for (int index = 0; index < m_gfnff_bonds.size(); ++index) {
        const auto& bond = m_gfnff_bonds[index];

        Vector i = geom().row(bond.i);
        Vector j = geom().row(bond.j);
        Matrix derivate;
        double rij = UFF::BondStretching(i, j, derivate, m_calculate_gradient);

        // Calculate r0 - either dynamic (using current CN) or static (from initialization)
        double r0_ij;
        if (use_dynamic_r0 && bond.z_i > 0 && bond.z_j > 0 &&
            bond.i < cn_size && bond.j < cn_size) {
            // Dynamic r0 calculation using current D3 coordination numbers
            // Reference: gfnff_method.cpp:1381 (initialization formula)
            // Formula: r0 = (ra + rb + rabshift) * ff
            double cn_i = d3cn(bond.i);
            double cn_j = d3cn(bond.j);

            // ra = r0(ati) + cnfak(ati) * cn(i)
            // rb = r0(atj) + cnfak(atj) * cn(j)
            double ra = bond.r0_base_i + bond.cnfak_i * cn_i;
            double rb = bond.r0_base_j + bond.cnfak_j * cn_j;

            // r0 = (ra + rb + rabshift) * ff (EXACT Fortran formula)
            // Reference: gfnff_rab.f90:153: rab(k) = (ra + rb + rab(k)) * ff
            // CRITICAL FIX (Jan 31, 2026): Reverting to INSIDE formula which matches XTB gradient source!
            r0_ij = (ra + rb + bond.rabshift) * bond.ff;
        } else {
            // Static r0 from initialization (fallback)
            r0_ij = bond.r0_ij;
        }

        // GFN-FF exponential bond stretching: E = k_b * exp(-α * (r-r₀)²)
        // Note: k_b is stored as NEGATIVE in gfnff.cpp:1358 for attractive bonds
        double dr = rij - r0_ij;                 // r - r₀ (dynamic or static)
        double alpha_orig = bond.exponent;       // α_orig before HB modification
        double alpha = alpha_orig;               // α stored in exponent field
        double k_b = bond.fc;                    // Force constant (already negative!)

        // Claude Generated (Jan 24, 2026): egbond_hb - Modified exponent for HB X-H bonds
        // Reference: Fortran gfnff_engrad.F90:957-958
        // For X-H bonds participating in hydrogen bridges (nr_hb >= 1):
        //   t1 = 1.0 - VBOND_SCALE = 0.1
        //   alpha_modified = (1.0 - t1 * hb_cn_H) * alpha = (1.0 - 0.1 * hb_cn_H) * alpha
        // This reduces the exponent by ~10% for typical H with CN=1, weakening the bond well
        if (bond.nr_hb >= 1) {
            constexpr double VBOND_SCALE = 0.9;
            double t1 = 1.0 - VBOND_SCALE;  // = 0.1
            alpha = (-t1 * bond.hb_cn_H + 1.0) * alpha_orig;  // Exact Fortran formula
        }

        double exp_term = std::exp(-alpha * dr * dr);
        double energy = k_b * exp_term;          // k_b already contains sign
        m_bond_energy += energy * factor;

        // Claude Generated (Feb 14, 2026): Per-bond CSV diagnostic for parameter comparison
        // Format matches Fortran analyzer "b" mode output for direct comparison
        if (CurcumaLogger::get_verbosity() >= 3) {
            if (index == 0) {
                CurcumaLogger::info("BOND_CSV: idx, atom_i, atom_j, Z_i, Z_j, rij, r0, fc, alpha, fqq, energy");
            }
            CurcumaLogger::info(fmt::format("BOND_CSV: {:3d}, {:3d}, {:3d}, {:2d}, {:2d}, {:.6f}, {:.6f}, {:.9f}, {:.9f}, {:.6f}, {:.12f}",
                                             index, bond.i, bond.j, bond.z_i, bond.z_j,
                                             rij, r0_ij, k_b, alpha, bond.fqq, energy * factor));
        }

        if (m_calculate_gradient) {
            // dE/dr = -2α * dr * E (exact chain rule formula)
            // E = k_b * exp(-α*dr²) where k_b < 0 (attractive)
            // dE/dr = k_b * (-2α*dr) * exp(-α*dr²) = -2α * dr * E
            double dEdr = -2.0 * alpha * dr * energy;  // Correct sign
            m_gradient.row(bond.i) += dEdr * factor * derivate.row(0);
            m_gradient.row(bond.j) += dEdr * factor * derivate.row(1);

            // Claude Generated (Feb 22, 2026): HB alpha-modulation chain-rule gradient
            // Reference: Fortran gfnff_engrad.F90:1054-1063 (egbond_hb, hb_dcn term)
            //
            // When alpha_mod = (1 - t1*hb_cn_H)*alpha_orig, the gradient has an extra term:
            //   dE/d(hb_cn_H) = dE/d(alpha_mod) * d(alpha_mod)/d(hb_cn_H)
            //                 = (-dr²*energy) * (-t1*alpha_orig) = t1*alpha_orig*dr²*energy
            //   zz = t1 * alpha_orig * dr² * energy * factor
            //
            // For each (H, B) pair: grad_H += zz * dCN_dH, grad_B += zz * dCN_dB
            if (bond.nr_hb >= 1) {
                constexpr double t1 = 0.1;  // = 1 - VBOND_SCALE = 1 - 0.9
                // alpha_orig already saved above (before HB modification)
                double zz = t1 * alpha_orig * dr * dr * energy * factor;

                int H = (m_atom_types[bond.i] == 1) ? bond.i : bond.j;
                for (const auto& hbg : m_hb_grad_entries) {
                    if (hbg.H_atom != H) continue;
                    m_gradient.row(H) += zz * hbg.dCN_dH.transpose();
                    m_gradient.row(hbg.B_atom) += zz * hbg.dCN_dB.transpose();
                }
            }

            // Claude Generated (Feb 15, 2026): Bond dr0/dCN chain-rule contribution
            // Reference: Fortran gfnff_engrad.F90:973-974, gfnff_rab.f90:147-156
            //
            // r0 = (r0_base_i + cnfak_i*cn_i + r0_base_j + cnfak_j*cn_j + rabshift) * ff
            // dr0/dCN_i = ff * cnfak_i, dr0/dCN_j = ff * cnfak_j
            // dE/dr0 = -dEdr (opposite sign: stretching r increases E, stretching r0 decreases E)
            // yy = -dEdr = 2*α*dr*E (Fortran convention: positive for stretched bond)
            //
            // Accumulate: dEdcn(atom) += yy * ff * cnfak
            if (use_dynamic_r0 && bond.z_i > 0 && bond.z_j > 0) {
                double yy = -dEdr;  // = 2*alpha*dr*energy (Fortran convention)
                m_dEdcn(bond.i) += yy * factor * bond.ff * bond.cnfak_i;
                m_dEdcn(bond.j) += yy * factor * bond.ff * bond.cnfak_j;
                // Claude Generated (Mar 2026): Bond-specific dEdcn for component gradient attribution
                m_dEdcn_bond(bond.i) += yy * factor * bond.ff * bond.cnfak_i;
                m_dEdcn_bond(bond.j) += yy * factor * bond.ff * bond.cnfak_j;
            }
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
        auto i = geom().row(angle.i);
        auto j = geom().row(angle.j);
        auto k = geom().row(angle.k);
        Matrix derivate;
        double costheta = UFF::AngleBending(i, j, k, derivate, m_calculate_gradient);

        // Clamp cosine to valid range [-1, 1] for numerical stability
        costheta = std::max(-1.0, std::min(1.0, costheta));

        double theta = std::acos(costheta);    // Current angle [0, π]
        double theta0 = angle.theta0_ijk;      // Equilibrium angle [0, π]
        double k_ijk = angle.fc;               // Force constant

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

        // Angle geometry logging - Claude Generated 2025-11-30, extended Feb 2026
        if (CurcumaLogger::get_verbosity() >= 3) {
            double r_ij = std::sqrt(r_ij_sq);
            double r_jk = std::sqrt(r_jk_sq);
            CurcumaLogger::info(fmt::format(
                "Angle #{}: atoms {}-{}-{} | theta={:.6f} rad ({:.2f}°), theta0={:.6f} rad ({:.2f}°) | "
                "k_ijk={:.8f} | r_ij={:.6f}, r_jk={:.6f}",
                index, angle.i, angle.j, angle.k, theta, theta*180.0/pi, theta0, theta0*180.0/pi,
                k_ijk, r_ij, r_jk));
        }

        // Get covalent radii: D3-style (gfnff_param.f90:381-404), 4/3 scaling applied at runtime
        // Reference: GFNFFParameters::covalent_rad_d3 stores raw_Å * aatoau (WITHOUT 4/3 factor)
        constexpr double rcov_scale_angle = 4.0 / 3.0;
        auto get_rcov_bohr = [&](int atomic_number) -> double {
            if (atomic_number >= 1 && atomic_number <= static_cast<int>(GFNFFParameters::covalent_rad_d3.size())) {
                return GFNFFParameters::covalent_rad_d3[atomic_number - 1] * rcov_scale_angle;
            }
            return 1.0 * GFNFFParameters::gfnff_aatoau * rcov_scale_angle;  // Fallback
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
        //
        // CRITICAL FIX (Feb 2026): Guard against division by near-zero squared distances
        double damp2ij = (r_ij_sq > 1e-8) ? -2.0 * 2.0 * rr_ij / (r_ij_sq * (1.0 + rr_ij) * (1.0 + rr_ij)) : 0.0;
        double damp2jk = (r_jk_sq > 1e-8) ? -2.0 * 2.0 * rr_jk / (r_jk_sq * (1.0 + rr_jk) * (1.0 + rr_jk)) : 0.0;

        // Phase 3: Apply distance-dependent damping to energy
        // Fortran formula: e = ea*damp where damp = damp_ij*damp_jk
        double angle_contribution = energy * damp * factor;
        m_angle_energy += angle_contribution;

        // DEBUG LOGGING - Claude Generated 2025-11-30, extended Feb 2026
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info(fmt::format(
                "  → energy_raw={:.10f} Eh, damp_ij={:.8f}, damp_jk={:.8f}, damp_total={:.8f}, "
                "factor={:.6f}, contribution={:.10f} Eh, cumulative={:.10f} Eh",
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
            // energy is raw energy, we need to apply damping from the other bond
            // term1 = E_raw * (d_damp_ij/d_xi) * damp_jk * factor
            // term2 = E_raw * (d_damp_jk/d_xk) * damp_ij * factor
            Vector term1 = energy * damp2ij * damp_jk * factor * vab;
            Vector term2 = energy * damp2jk * damp_ij * factor * vcb;

            // Angle gradient contributions (already computed by UFF::AngleBending)
            // Note: derivate matrix contains ∂θ/∂x for all 3 atoms
            Vector grad_angle_i = dedtheta * damp * factor * derivate.row(0);
            Vector grad_angle_j = dedtheta * damp * factor * derivate.row(1);
            Vector grad_angle_k = dedtheta * damp * factor * derivate.row(2);

            // Complete gradients with damping terms (Fortran lines 913-915)
            // g(:,1) =  deda*damp + term1          (Atom i: terminal)
            // g(:,2) = -dedb*damp - term1 - term2  (Atom j: center)
            // g(:,3) =  dedc*damp + term2          (Atom k: terminal)
            m_gradient.row(angle.i) += grad_angle_i + term1;
            m_gradient.row(angle.j) += grad_angle_j - term1 - term2;
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
        Eigen::Vector3d r_i = geom().row(dihedral.i).head<3>();
        Eigen::Vector3d r_j = geom().row(dihedral.j).head<3>();
        Eigen::Vector3d r_k = geom().row(dihedral.k).head<3>();
        Eigen::Vector3d r_l = geom().row(dihedral.l).head<3>();

        // Standard dihedral angle for bonded chain i-j-k-l
        // Reference: gfnff_engrad.F90:1260 — phi = valijklff(n,xyz,i,j,k,l)
        Matrix derivate;
        double phi = GFNFF_Geometry::calculateDihedralAngle(r_i, r_j, r_k, r_l, derivate, m_calculate_gradient);

        // Mapping for i-j-k-l call: derivate row 0=i, 1=j, 2=k, 3=l

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

        // Damping distances in Fortran ll-ii-jj-kk convention (Feb 4, 2026)
        // Reference: gfnff_engrad.F90:1250-1259 (primary torsions, tlist(5)>0)
        // After atom-ordering fix: i=ll, j=ii, k=jj, l=kk
        //   i-j = ll-ii: 1-3 distance (ll is neighbor of jj, NOT of ii)
        //   j-k = ii-jj: central bond
        //   k-l = jj-kk: 1-3 distance (kk is neighbor of ii, NOT of jj)
        //   damp = dampij * dampjk * dampkl
        Eigen::Vector3d vab = r_i - r_j;   // ll-ii: 1-3 distance
        Eigen::Vector3d vcb = r_j - r_k;   // ii-jj: central bond
        Eigen::Vector3d vdc = r_k - r_l;   // jj-kk: 1-3 distance

        double rij2 = vab.squaredNorm();
        double rjk2 = vcb.squaredNorm();
        double rkl2 = vdc.squaredNorm();

        // Damping function (gfnff_engrad.F90:1422-1431, gfnffdampt subroutine)
        auto calculate_damping = [&](double r2, double rcov_a, double rcov_b, bool use_nci = false) -> double {
            const double atcutt = use_nci ? GFNFFParameters::atcutt_nci : GFNFFParameters::atcutt;
            double rcut = atcutt * (rcov_a + rcov_b) * (rcov_a + rcov_b);
            double rr = (r2 / rcut) * (r2 / rcut);
            return 1.0 / (1.0 + rr);
        };

        // Covalent radii for damping (D3 radii scaled by 4/3, as in gfnff_param.f90:381-405)
        constexpr double rcov_scale = 4.0 / 3.0;
        double rcov_i = GFNFFParameters::covalent_rad_d3[Z_i - 1] * rcov_scale;
        double rcov_j = GFNFFParameters::covalent_rad_d3[Z_j - 1] * rcov_scale;
        double rcov_k = GFNFFParameters::covalent_rad_d3[Z_k - 1] * rcov_scale;
        double rcov_l = GFNFFParameters::covalent_rad_d3[Z_l - 1] * rcov_scale;

        // Three damping factors on consecutive bonds (gfnff_engrad.F90:1256-1259)
        double damp_ij = calculate_damping(rij2, rcov_i, rcov_j);
        double damp_jk = calculate_damping(rjk2, rcov_j, rcov_k);
        double damp_kl = calculate_damping(rkl2, rcov_k, rcov_l);
        double damp = damp_ij * damp_jk * damp_kl;

        // DEBUG: Print damping details for first primary torsion
        static bool primary_damp_debug_printed = false;
        if (!primary_damp_debug_printed && index == 0 && CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info(fmt::format("\n=== GFN-FF PRIMARY TORSION DAMPING (Torsion {}-{}-{}-{}) ===",
                                             dihedral.i, dihedral.j, dihedral.k, dihedral.l));
            CurcumaLogger::info(fmt::format("  r_ij = {:.6f} Bohr (bond i-j)", std::sqrt(rij2)));
            CurcumaLogger::info(fmt::format("  r_jk = {:.6f} Bohr (bond j-k)", std::sqrt(rjk2)));
            CurcumaLogger::info(fmt::format("  r_kl = {:.6f} Bohr (bond k-l)", std::sqrt(rkl2)));
            CurcumaLogger::info(fmt::format("  damp_ij = {:.8f}, damp_jk = {:.8f}, damp_kl = {:.8f}",
                                             damp_ij, damp_jk, damp_kl));
            CurcumaLogger::info(fmt::format("  damp (total) = {:.8f}", damp));
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

        // Claude Generated (Feb 4, 2026): ROOT CAUSE RESOLVED - atom ordering
        // All previous "mystery factors" (3/11, 0.5, molecule-dependent errors) were caused by
        // incorrect atom ordering in torsion storage. Torsions were stored as i-j-k-l (bonded chain)
        // but Fortran uses ll-ii-jj-kk where pos0-pos1 and pos2-pos3 are 1-3 distances.
        // This made damping ~40× too large (bond distances vs 1-3 distances).
        // Fix: Swap terminal atoms in gfnff_torsions.cpp storage (Feb 4, 2026).
        // CH₃OCH₃ torsion: 0.001 Eh → 0.000023 Eh (exact match with XTB reference).
        // Reference: gfnff_engrad.F90:1268-1272: et = (1+cos(c1)) * vtors(2,m)
        double et = V * (1.0 + cos(c1));
        double energy = et * damp;

        primary_torsion_energy += energy * m_dihedral_scaling;  // Claude Generated (Jan 2, 2026): Accumulate primary

        // Torsion geometry analysis (December 2025) - Significant energy torsions
        if (std::abs(energy) > 1e-4) {
            if (CurcumaLogger::get_verbosity() >= 3) {
                CurcumaLogger::info(fmt::format("\nSignificant Torsion analysis (Atom {}-{}-{}-{}):",
                                                 dihedral.i, dihedral.j, dihedral.k, dihedral.l));
                CurcumaLogger::info(fmt::format("  phi_actual = {:.2f}°", phi * 180.0 / M_PI));
                CurcumaLogger::info(fmt::format("  phi0 = {:.2f}°", phi0 * 180.0 / M_PI));
                CurcumaLogger::info(fmt::format("  n = {:.0f}", n));
                CurcumaLogger::info(fmt::format("  V = {:.6f} Eh", V));
                CurcumaLogger::info(fmt::format("  damp = {:.6f}", damp));
                CurcumaLogger::info(fmt::format("  Energy = {:.6f} Eh\n", energy));
            }
        }

        // NO accumulation here - already done above at line 1056
        // m_dihedral_energy += energy * m_dihedral_scaling;  // Claude Generated (Jan 2, 2026): REMOVED - moved above

        if (m_calculate_gradient) {
            // CRITICAL FIX (Feb 2026): Stricter threshold for aromatic systems
            // Small V parameters combined with numerical noise can amplify errors
            // Raised from 1e-10 to 1e-6 for aromatic stability
            if (std::abs(V) < 1e-6) {
                // V too small → gradient negligible, skip calculation
                // This prevents numerical noise amplification in benzene-like systems
                continue;
            }

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

            // First, compute ddamp/∂r² for each bond
            // Formula: ddamp = -4*rr / (r² * (1+rr)²)
            // CRITICAL FIX (Feb 2026): Guard against division by near-zero r2_val
            auto calc_ddamp = [&](double r2_val, double rcut_val) -> double {
                // Prevent division by zero when atoms are extremely close
                if (r2_val < 1e-8) {
                    return 0.0;
                }
                double rr_val = (r2_val / rcut_val) * (r2_val / rcut_val);  // (r²/rcut)²
                double one_plus_rr_val = 1.0 + rr_val;
                return -4.0 * rr_val / (r2_val * one_plus_rr_val * one_plus_rr_val);  // ∂damp/∂r²
            };

            // rcut values for consecutive bonds (same pairs as energy calculation)
            bool use_nci = dihedral.is_nci;
            double atcutt_value = use_nci ? GFNFFParameters::atcutt_nci : GFNFFParameters::atcutt;
            double rcut_ij = atcutt_value * (rcov_i + rcov_j) * (rcov_i + rcov_j);
            double rcut_jk = atcutt_value * (rcov_j + rcov_k) * (rcov_j + rcov_k);
            double rcut_kl = atcutt_value * (rcov_k + rcov_l) * (rcov_k + rcov_l);

            // Damping derivatives w.r.t. squared bond distances
            double ddamp_drij2 = calc_ddamp(rij2, rcut_ij);
            double ddamp_drjk2 = calc_ddamp(rjk2, rcut_jk);
            double ddamp_drkl2 = calc_ddamp(rkl2, rcut_kl);

            // Damping gradient terms following Fortran gfnff_engrad.F90:1268-1270
            //   term1 = et*damp2ij*dampjk*dampkl*vab   (vab = i-j vector)
            //   term2 = et*damp2jk*dampij*dampkl*vcb   (vcb = j-k vector)
            //   term3 = et*damp2kl*dampij*dampjk*vdc   (vdc = k-l vector)
            Vector term_ab = (et * damp_jk * damp_kl * ddamp_drij2) * vab;
            Vector term_bc = (et * damp_ij * damp_kl * ddamp_drjk2) * vcb;
            Vector term_cd = (et * damp_ij * damp_jk * ddamp_drkl2) * vdc;

            // =====================================================================
            // Combine angle and damping gradient (Fortran gfnff_engrad.F90:1271-1274)
            // =====================================================================
            // g(1) = dij*dda + term1          ← atom i: +term_ab
            // g(2) = dij*ddb - term1 + term2  ← atom j: -term_ab + term_bc
            // g(3) = dij*ddc + term3 - term2  ← atom k: +term_cd - term_bc
            // g(4) = dij*ddd - term3          ← atom l: -term_cd
            // derivate mapping (i-j-k-l call): row 0=i, 1=j, 2=k, 3=l

            double t_scaling = m_final_factor * m_dihedral_scaling;

            m_gradient.row(dihedral.i) += (dEdphi * derivate.row(0)).transpose() + term_ab * t_scaling;
            m_gradient.row(dihedral.j) += (dEdphi * derivate.row(1)).transpose() + (-term_ab + term_bc) * t_scaling;
            m_gradient.row(dihedral.k) += (dEdphi * derivate.row(2)).transpose() + (term_cd - term_bc) * t_scaling;
            m_gradient.row(dihedral.l) += (dEdphi * derivate.row(3)).transpose() - term_cd * t_scaling;

            // DEBUG
            static bool gradient_debug_printed = false;
            if (!gradient_debug_printed && index == 0 && CurcumaLogger::get_verbosity() >= 3) {
                CurcumaLogger::info(fmt::format("\n=== PRIMARY TORSION GRADIENT (Torsion {}-{}-{}-{}) ===",
                                                 dihedral.i, dihedral.j, dihedral.k, dihedral.l));
                CurcumaLogger::info(fmt::format("  damp = {:.8f}, et = {:.6e} Eh", damp, et));
                CurcumaLogger::info(fmt::format("  |term_ab| = {:.6e}, |term_bc| = {:.6e}, |term_cd| = {:.6e}",
                                                 term_ab.norm(), term_bc.norm(), term_cd.norm()));
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
    // Debug for torsion investigation (verbosity >= 3 only)
    if (CurcumaLogger::get_verbosity() >= 3 && m_thread == 0) {
        CurcumaLogger::info(fmt::format("DEBUG_THREAD0 primary_torsion: {} terms, E={:.8e}",
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
        Eigen::Vector3d r_i = geom().row(torsion.i).head<3>();
        Eigen::Vector3d r_j = geom().row(torsion.j).head<3>();
        Eigen::Vector3d r_k = geom().row(torsion.k).head<3>();
        Eigen::Vector3d r_l = geom().row(torsion.l).head<3>();

        // Standard dihedral angle — extra torsions have tlist(5)=1 > 0, so they
        // enter the same branch as primary torsions in gfnff_engrad.F90:1249
        // using valijklff(n,xyz,i,j,k,l), NOT omega() (which is for inversions)
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

        // Consecutive damping — identical to primary torsions (Feb 4, 2026)
        // Extra torsions have tlist(5)=1 > 0, entering the same branch as primary
        // in gfnff_engrad.F90:1249-1260: vab=xyz(i)-xyz(j), vcb=xyz(j)-xyz(k), vdc=xyz(k)-xyz(l)
        Eigen::Vector3d vab = r_i - r_j;   // ll-ii: 1-3 distance
        Eigen::Vector3d vcb = r_j - r_k;   // ii-jj: central bond
        Eigen::Vector3d vdc = r_k - r_l;   // jj-kk: 1-3 distance

        double rij2 = vab.squaredNorm();
        double rjk2 = vcb.squaredNorm();
        double rkl2 = vdc.squaredNorm();

        // Damping function
        auto calculate_damping = [&](double r2, double rcov_a, double rcov_b, bool use_nci_flag = false) -> double {
            const double atcutt_val = use_nci_flag ? GFNFFParameters::atcutt_nci : GFNFFParameters::atcutt;
            double rcut_val = atcutt_val * (rcov_a + rcov_b) * (rcov_a + rcov_b);
            double rr_val = (r2 / rcut_val) * (r2 / rcut_val);
            return 1.0 / (1.0 + rr_val);
        };

        // Covalent radii for damping (D3 radii scaled by 4/3)
        constexpr double rcov_scale = 4.0 / 3.0;
        double rcov_i = GFNFFParameters::covalent_rad_d3[Z_i - 1] * rcov_scale;
        double rcov_j = GFNFFParameters::covalent_rad_d3[Z_j - 1] * rcov_scale;
        double rcov_k = GFNFFParameters::covalent_rad_d3[Z_k - 1] * rcov_scale;
        double rcov_l = GFNFFParameters::covalent_rad_d3[Z_l - 1] * rcov_scale;

        // Consecutive damping factors — same as primary (gfnff_engrad.F90:1256-1259)
        bool use_nci = torsion.is_nci;
        double damp_ij = calculate_damping(rij2, rcov_i, rcov_j, use_nci);
        double damp_jk = calculate_damping(rjk2, rcov_j, rcov_k, use_nci);
        double damp_kl = calculate_damping(rkl2, rcov_k, rcov_l, use_nci);
        double damp = damp_ij * damp_jk * damp_kl;

        // DEBUG
        static bool damp_debug_printed = false;
        if (!damp_debug_printed && index == 0 && CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info(fmt::format("\n=== EXTRA TORSION DAMPING (Torsion {}-{}-{}-{}) ===",
                                             torsion.i, torsion.j, torsion.k, torsion.l));
            CurcumaLogger::info(fmt::format("  r_ij = {:.6f}, r_jk = {:.6f}, r_kl = {:.6f} Bohr",
                                             std::sqrt(rij2), std::sqrt(rjk2), std::sqrt(rkl2)));
            CurcumaLogger::info(fmt::format("  damp = {:.8f}", damp));
            damp_debug_printed = true;
        }

        // Energy calculation for EXTRA torsions
        // Reference: Fortran gfnff_engrad.F90:1268-1272
        // CRITICAL: ALL torsions (primary and extra) use c1 = n*(phi - phi0) + π
        double V = torsion.V;
        double n = torsion.n;
        double phi0 = torsion.phi0;
        double c1 = n * (phi - phi0) + M_PI;

        double et = V * (1.0 + cos(c1));
        double energy = et * damp;

        // Accumulate extra torsion energy separately
        extra_torsion_energy += energy * m_dihedral_scaling;

        if (m_calculate_gradient) {
            // Analytical gradient — same formula as primary torsions (Feb 4, 2026)
            double dEdphi = -V * n * sin(c1) * damp * m_dihedral_scaling;

            auto calc_ddamp = [&](double r2_v, double rcut_v) -> double {
                if (r2_v < 1e-8) return 0.0;
                double rr_v = (r2_v / rcut_v) * (r2_v / rcut_v);
                double one_plus_rr_v = 1.0 + rr_v;
                return -4.0 * rr_v / (r2_v * one_plus_rr_v * one_plus_rr_v);
            };

            // Consecutive damping rcut values (gfnff_engrad.F90:1256-1259)
            double atcutt_val = use_nci ? GFNFFParameters::atcutt_nci : GFNFFParameters::atcutt;
            double rcut_ij = atcutt_val * (rcov_i + rcov_j) * (rcov_i + rcov_j);
            double rcut_jk = atcutt_val * (rcov_j + rcov_k) * (rcov_j + rcov_k);
            double rcut_kl = atcutt_val * (rcov_k + rcov_l) * (rcov_k + rcov_l);

            double ddamp_drij2 = calc_ddamp(rij2, rcut_ij);
            double ddamp_drjk2 = calc_ddamp(rjk2, rcut_jk);
            double ddamp_drkl2 = calc_ddamp(rkl2, rcut_kl);

            // Gradient terms following Fortran gfnff_engrad.F90:1271-1274
            //   term1 = et*damp2ij*dampjk*dampkl*vab
            //   term2 = et*damp2jk*dampij*dampkl*vcb
            //   term3 = et*damp2kl*dampij*dampjk*vdc
            Vector term_ab = (et * damp_jk * damp_kl * ddamp_drij2) * vab;
            Vector term_bc = (et * damp_ij * damp_kl * ddamp_drjk2) * vcb;
            Vector term_cd = (et * damp_ij * damp_jk * ddamp_drkl2) * vdc;

            double t_scal = m_final_factor * m_dihedral_scaling;

            // Consecutive damping gradient signs (gfnff_engrad.F90:1271-1274)
            // g(1) = dij*dda + term1   ← atom i
            // g(2) = dij*ddb - term1 + term2  ← atom j
            // g(3) = dij*ddc + term3 - term2  ← atom k
            // g(4) = dij*ddd - term3   ← atom l
            // derivate mapping (i-j-k-l call): row 0=i, 1=j, 2=k, 3=l
            m_gradient.row(torsion.i) += (dEdphi * derivate.row(0)).transpose() + term_ab * t_scal;
            m_gradient.row(torsion.j) += (dEdphi * derivate.row(1)).transpose() + (-term_ab + term_bc) * t_scal;
            m_gradient.row(torsion.k) += (dEdphi * derivate.row(2)).transpose() + (term_cd - term_bc) * t_scal;
            m_gradient.row(torsion.l) += (dEdphi * derivate.row(3)).transpose() - term_cd * t_scal;
        }
    }

    // Claude Generated (Jan 2, 2026): Add extra torsion energy to total and output at verbosity 2
    m_dihedral_energy += extra_torsion_energy;
    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::result(fmt::format("Extra sp3-sp3 torsions: {} terms, energy = {:.6e} Eh",
                                           m_gfnff_extra_torsions.size(), extra_torsion_energy));
        CurcumaLogger::result(fmt::format("Total dihedral energy: {:.6e} Eh (primary + extra)", m_dihedral_energy));
    }
    // Debug for torsion investigation (verbosity >= 3 only)
    if (CurcumaLogger::get_verbosity() >= 3 && m_thread == 0) {
        CurcumaLogger::info(fmt::format("DEBUG_THREAD0 extra_torsion: {} terms, E={:.8e}",
                                       m_gfnff_extra_torsions.size(), extra_torsion_energy));
    }
}

void ForceFieldThread::CalculateGFNFFInversionContribution()
{
    // Claude Generated (Feb 2026): Rewritten to match Fortran egtors inversion branch
    // Reference: gfnff_engrad.F90:1355-1387
    //
    // GFN-FF inversions use the omega (out-of-plane angle) function with star-topology damping.
    // Atom convention: i=center, j/k/l=three neighbors (sorted by distance in generation)
    //
    // Two potential types:
    //   potential_type=0  (tlist(5)=0): E = V*(1-cos(omega)) * damp   [planar sp2 centers]
    //   potential_type=-1 (tlist(5)<0): E = V*(cos(omega)-cos(omega0))^2 * damp  [saturated N]

    for (int index = 0; index < m_gfnff_inversions.size(); ++index) {
        const auto& inv = m_gfnff_inversions[index];

        // For UFF/QMDFF inversions (type != 3), use old formula
        if (inv.type != 3) {
            Eigen::Vector3d r_i = geom().row(inv.i).head<3>();
            Eigen::Vector3d r_j = geom().row(inv.j).head<3>();
            Eigen::Vector3d r_k = geom().row(inv.k).head<3>();
            Eigen::Vector3d r_l = geom().row(inv.l).head<3>();
            Matrix derivate;
            double theta = GFNFF_Geometry::calculateOutOfPlaneAngle(r_i, r_j, r_k, r_l, derivate, m_calculate_gradient);
            double energy = inv.fc * (inv.C0 + inv.C1 * cos(theta) + inv.C2 * cos(2 * theta));
            m_inversion_energy += energy * m_final_factor * m_inversion_scaling;
            if (m_calculate_gradient) {
                double dEdtheta = inv.fc * (-inv.C1 * sin(theta) - 2 * inv.C2 * sin(2 * theta)) * m_final_factor * m_inversion_scaling;
                m_gradient.row(inv.i) += dEdtheta * derivate.row(0);
                m_gradient.row(inv.j) += dEdtheta * derivate.row(1);
                m_gradient.row(inv.k) += dEdtheta * derivate.row(2);
                m_gradient.row(inv.l) += dEdtheta * derivate.row(3);
            }
            continue;
        }

        // =====================================================================
        // GFN-FF inversion (type == 3)
        // =====================================================================
        // Atom layout: i=center (3 neighbors), j=nb1, k=nb2, l=nb3
        Eigen::Vector3d r_center = geom().row(inv.i).head<3>();
        Eigen::Vector3d r_nb1 = geom().row(inv.j).head<3>();
        Eigen::Vector3d r_nb2 = geom().row(inv.k).head<3>();
        Eigen::Vector3d r_nb3 = geom().row(inv.l).head<3>();

        // Out-of-plane angle (Fortran omega function, gfnff_helpers.f90:427-448)
        // Center atom is i, calculateOutOfPlaneAngle expects i=center
        Matrix derivate;
        double omega = GFNFF_Geometry::calculateOutOfPlaneAngle(r_center, r_nb1, r_nb2, r_nb3, derivate, m_calculate_gradient);

        // =====================================================================
        // Star-topology damping (3 distances from center to each neighbor)
        // Reference: gfnff_engrad.F90:1372-1381
        // =====================================================================
        int Z_center = m_atom_types[inv.i];
        int Z_nb1 = m_atom_types[inv.j];
        int Z_nb2 = m_atom_types[inv.k];
        int Z_nb3 = m_atom_types[inv.l];

        if (Z_center < 1 || Z_center > 86 || Z_nb1 < 1 || Z_nb1 > 86 ||
            Z_nb2 < 1 || Z_nb2 > 86 || Z_nb3 < 1 || Z_nb3 > 86) {
            continue;
        }

        constexpr double rcov_scale = 4.0 / 3.0;
        double rcov_c = GFNFFParameters::covalent_rad_d3[Z_center - 1] * rcov_scale;
        double rcov_1 = GFNFFParameters::covalent_rad_d3[Z_nb1 - 1] * rcov_scale;
        double rcov_2 = GFNFFParameters::covalent_rad_d3[Z_nb2 - 1] * rcov_scale;
        double rcov_3 = GFNFFParameters::covalent_rad_d3[Z_nb3 - 1] * rcov_scale;

        // Fortran inversion damping (gfnff_engrad.F90:1356-1365):
        // NOT star topology! Uses nb1 as hub for two of three distances:
        //   vab = xyz(j) - xyz(i)  → nb1 minus center    (bond distance)
        //   vcb = xyz(j) - xyz(k)  → nb1 minus nb2       (1-3 distance across angle)
        //   vdc = xyz(j) - xyz(l)  → nb1 minus nb3       (1-3 distance across angle)
        // This gives much smaller damping than star topology because 1-3 distances
        // (~4.5 Bohr) damp much more than bond distances (~2.65 Bohr).
        double rij_sq = (r_nb1 - r_center).squaredNorm();  // nb1-center: bond
        double rjk_sq = (r_nb1 - r_nb2).squaredNorm();     // nb1-nb2: 1-3 distance
        double rjl_sq = (r_nb1 - r_nb3).squaredNorm();     // nb1-nb3: 1-3 distance

        auto calculate_damping = [](double r2, double rcov_a, double rcov_b) -> double {
            double rcut = GFNFFParameters::atcutt * (rcov_a + rcov_b) * (rcov_a + rcov_b);
            double rr = (r2 / rcut) * (r2 / rcut);
            return 1.0 / (1.0 + rr);
        };

        // Fortran: damp = dampij * dampjk * dampjl
        // gfnffdampt(at(i),at(j),rij) → center, nb1 pair
        // gfnffdampt(at(k),at(j),rjk) → nb2, nb1 pair
        // gfnffdampt(at(j),at(l),rjl) → nb1, nb3 pair
        double damp_ij = calculate_damping(rij_sq, rcov_c, rcov_1);  // center-nb1
        double damp_jk = calculate_damping(rjk_sq, rcov_2, rcov_1);  // nb2-nb1
        double damp_jl = calculate_damping(rjl_sq, rcov_1, rcov_3);  // nb1-nb3
        double damp = damp_ij * damp_jk * damp_jl;

        // =====================================================================
        // Energy calculation
        // =====================================================================
        double V = inv.fc;
        double et = 0.0;       // Energy before damping (needed for damping gradient)
        double dEdomega = 0.0; // dE/domega including damp factor

        if (inv.potential_type == 0) {
            // Planar sp2: E = V*(1 - cos(omega)) * damp
            // Reference: gfnff_engrad.F90:1373-1374
            et = V * (1.0 - cos(omega));
            dEdomega = V * sin(omega) * damp;
        } else {
            // Saturated N (double minima at ±omega0):
            // E = V*(cos(omega) - cos(omega0))^2 * damp
            // Reference: gfnff_engrad.F90:1376-1377
            double diff = cos(omega) - cos(inv.omega0);
            et = V * diff * diff;
            dEdomega = -2.0 * V * sin(omega) * diff * damp;
        }

        double energy = et * damp;
        m_inversion_energy += energy * m_final_factor * m_inversion_scaling;

        if (m_calculate_gradient) {
            // =====================================================================
            // PART 1: Omega derivative term (dij * dda/ddb/ddc/ddd)
            // =====================================================================
            double grad_scale = dEdomega * m_final_factor * m_inversion_scaling;

            // =====================================================================
            // PART 2: Damping derivative terms (Fortran gfnff_engrad.F90:1379-1385)
            // =====================================================================
            // ddamp = d(damp)/d(r²) with factor of 2 built in:
            //   ddamp = -4*rr / (r² * (1+rr)²)
            // Matches Fortran gfnffdampt subroutine (gfnff_engrad.F90:1501-1510)
            auto calc_ddamp = [](double r2_val, double rcov_a, double rcov_b) -> double {
                if (r2_val < 1e-8) return 0.0;
                double rcut_val = GFNFFParameters::atcutt * (rcov_a + rcov_b) * (rcov_a + rcov_b);
                double rr_val = (r2_val / rcut_val) * (r2_val / rcut_val);
                double one_plus_rr = 1.0 + rr_val;
                return -4.0 * rr_val / (r2_val * one_plus_rr * one_plus_rr);
            };

            double ddamp_ij = calc_ddamp(rij_sq, rcov_c, rcov_1);
            double ddamp_jk = calc_ddamp(rjk_sq, rcov_2, rcov_1);
            double ddamp_jl = calc_ddamp(rjl_sq, rcov_1, rcov_3);

            // Damping vectors (same as energy: nb1 is hub)
            Eigen::Vector3d vab = r_nb1 - r_center;  // j - i
            Eigen::Vector3d vcb = r_nb1 - r_nb2;     // j - k
            Eigen::Vector3d vdc = r_nb1 - r_nb3;     // j - l

            // Fortran gfnff_engrad.F90:1379-1381
            Eigen::Vector3d term1 = (et * ddamp_ij * damp_jk * damp_jl) * vab;
            Eigen::Vector3d term2 = (et * ddamp_jk * damp_ij * damp_jl) * vcb;
            Eigen::Vector3d term3 = (et * ddamp_jl * damp_ij * damp_jk) * vdc;

            double t_scal = m_final_factor * m_inversion_scaling;

            // Fortran gfnff_engrad.F90:1382-1385 (inversion gradient distribution)
            // g(:,1) = dij*dda - term1               ← center
            // g(:,2) = dij*ddb + term1 + term2 + term3  ← nb1 (hub)
            // g(:,3) = dij*ddc - term2               ← nb2
            // g(:,4) = dij*ddd - term3               ← nb3
            m_gradient.row(inv.i) += (grad_scale * derivate.row(0)).transpose() - term1 * t_scal;
            m_gradient.row(inv.j) += (grad_scale * derivate.row(1)).transpose() + (term1 + term2 + term3) * t_scal;
            m_gradient.row(inv.k) += (grad_scale * derivate.row(2)).transpose() - term2 * t_scal;
            m_gradient.row(inv.l) += (grad_scale * derivate.row(3)).transpose() - term3 * t_scal;
        }
    }
}

void ForceFieldThread::CalculateGFNFFSTorsionContribution()
{
    /**
     * @brief GFN-FF Triple Bond Torsion (sTors_eg)
     *
     * Specialized torsion for rotation around sp-sp systems (alkynes).
     *
     * Reference: gfnff_engrad.F90:3454 - sTors_eg()
     * Potential: E = erefhalf * (1 - cos(2*phi))
     */

    for (const auto& stor : m_gfnff_storsions) {
        Matrix derivate;
        // phi in radians
        double phi = GFNFF_Geometry::calculateDihedralAngle(
            geom().row(stor.i).transpose() * m_au,
            geom().row(stor.j).transpose() * m_au,
            geom().row(stor.k).transpose() * m_au,
            geom().row(stor.l).transpose() * m_au,
            derivate, m_calculate_gradient);

        double erefhalf = stor.erefhalf;
        double energy = -erefhalf * std::cos(2.0 * phi) + erefhalf;

        m_stors_energy += energy * m_final_factor;

        if (m_calculate_gradient) {
            double dEdphi = 2.0 * erefhalf * std::sin(2.0 * phi) * m_final_factor;
            m_gradient.row(stor.i) += (dEdphi * derivate.row(0)).transpose();
            m_gradient.row(stor.j) += (dEdphi * derivate.row(1)).transpose();
            m_gradient.row(stor.k) += (dEdphi * derivate.row(2)).transpose();
            m_gradient.row(stor.l) += (dEdphi * derivate.row(3)).transpose();
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
     * @brief GFN-FF Dispersion with MODIFIED Becke-Johnson damping (pairwise parallelizable)
     *
     * CLAUDE GENERATED (January 25, 2026): Fix to match XTB 6.6.1 reference
     *
     * Reference: gfnff_gdisp0.f90:365-377
     *
     * GFN-FF uses a MODIFIED BJ damping formula (NOT standard D3/D4):
     *   E = -0.5 * C6 * (t6 + 2*r4r2ij*t8)
     * where:
     *   t6 = 1/(r^6 + R0^6)
     *   t8 = 1/(r^8 + R0^8)
     *   r4r2ij = 3 * sqrtZr4r2_i * sqrtZr4r2_j  (implicit C8/C6 ratio)
     *   R0^2 = (a1*sqrt(r4r2ij) + a2)^2 with a1=0.58, a2=4.80
     *
     * Key differences from standard BJ damping:
     * 1. R0 computed from sqrtZr4r2 product (NOT from C8/C6 ratio)
     * 2. C8 is implicit: factor 2*r4r2ij*t8 instead of separate C8*t8
     * 3. 0.5 factor for pair counting (each pair counted once)
     */

    if (CurcumaLogger::get_verbosity() >= 3 && m_gfnff_dispersions.size() > 0) {
        CurcumaLogger::info(fmt::format("Thread calculating {} GFN-FF dispersion pairs", m_gfnff_dispersions.size()));
    }
    // Claude Generated (Feb 8, 2026): Per-pair diagnostic for dispersion accuracy investigation (verbosity >= 3 only)
    // Enable for small molecules only (< 50 pairs) to avoid excessive output
    bool disp_diag = (CurcumaLogger::get_verbosity() >= 3 && m_gfnff_dispersions.size() > 0 && m_gfnff_dispersions.size() <= 50);
    if (disp_diag) {
        CurcumaLogger::info("DISP_CSV: idx,i,j,rij_bohr,C6,r4r2ij,r0_sq,zetac6,t6,t8,energy");
    }

    for (int index = 0; index < m_gfnff_dispersions.size(); ++index) {
        const auto& disp = m_gfnff_dispersions[index];

        Eigen::Vector3d ri = geom().row(disp.i);
        Eigen::Vector3d rj = geom().row(disp.j);
        Eigen::Vector3d rij_vec = ri - rj;
        double rij = rij_vec.norm() * m_au;  // Convert to atomic units if needed

        // HIGH PRIORITY FIX (Feb 2026): Reduce epsilon threshold for gradient robustness
        // Gradient has division by rij → strengthen guard from 1e-10 to 1e-8
        if (rij > disp.r_cut || rij < 1e-8) continue;  // Skip if beyond cutoff or too close

        // GFN-FF modified BJ damping (NOT standard D3/D4!)
        // Reference: gfnff_gdisp0.f90:365-377
        //
        // r0 = radii(lin(ati, atj))  // This is R0^2 (pre-computed)
        // t6 = 1._wp/(r2**3+r0**3)   // 1/(r^6 + R0^6)
        // t8 = 1._wp/(r2**4+r0**4)   // 1/(r^8 + R0^8)

        // Optimized power calculations: r^6 = (r^2)^3, r^8 = (r^2)^4
        double r2 = rij * rij;
        double r6 = r2 * r2 * r2;     // (r^2)^3 = r^6

        // R0^6 = (R0^2)^3 where r0_squared is pre-computed as (a1*sqrt(r4r2ij)+a2)^2
        double r0_6 = disp.r0_squared * disp.r0_squared * disp.r0_squared;

        // t6 = 1/(r^6 + R0^6)
        double t6 = 1.0 / (r6 + r0_6);

        // t8 = 1/(r^8 + R0^8)
        double r8 = r6 * r2;          // r^8 = r^6 * r^2
        double r0_8 = r0_6 * disp.r0_squared;  // R0^8 = R0^6 * R0^2
        double t8 = 1.0 / (r8 + r0_8);

        // GFN-FF dispersion formula: disp = (t6 + 2*r4r2ij*t8) * zetac6
        // Energy: dE = -c6 * disp (no 0.5 - Fortran accumulates to both atoms, we count once per pair)
        // Reference: gfnff_gdisp0.f90:374,377
        // Claude Generated (Jan 31, 2026): Added zetac6 charge scaling
        // Claude Generated (Feb 8, 2026): Removed 0.5 factor - Fortran's 0.5 is cancelled by dual accumulation
        double disp_sum = t6 + 2.0 * disp.r4r2ij * t8;
        double energy = -disp.C6 * disp_sum * disp.zetac6 * m_final_factor;

        m_dispersion_energy += energy;
        m_d4_energy += energy;  // GFN-FF dispersion reports as D4 energy

        // Debug: Per-pair dispersion breakdown for accuracy investigation (verbosity >= 3 only)
        if (disp_diag) {
            CurcumaLogger::info(fmt::format("DISP_CSV: {},{},{},{:.6f},{:.6e},{:.6f},{:.6f},{:.6e},{:.6e},{:.6e},{:.10e}",
                                           index, disp.i, disp.j,
                                           rij, disp.C6, disp.r4r2ij, disp.r0_squared, disp.zetac6,
                                           t6, t8, energy));
        }

        if (m_calculate_gradient) {
            // Analytical gradient for GFN-FF dispersion formula
            // Reference: gfnff_gdisp0.f90:371-372

            double d6 = -6.0 * r2 * r2 * t6 * t6;
            double d8 = -8.0 * r2 * r2 * r2 * t8 * t8;

            double ddisp_dr2 = d6 + 2.0 * disp.r4r2ij * d8;
            double dEdr = -disp.C6 * disp.zetac6 * ddisp_dr2 * rij * m_final_factor;

            Eigen::Vector3d grad = dEdr * rij_vec / rij;

            m_gradient.row(disp.i) += grad.transpose();
            m_gradient.row(disp.j) -= grad.transpose();

            // Claude Generated (Feb 15, 2026): Dispersion dC6/dCN chain-rule contribution
            // Reference: Fortran gfnff_gdisp0.f90:382-395
            //
            // dEdcn(iat) -= dc6dcn(iat,jat) * disp_value
            // dEdcn(jat) -= dc6dcn(jat,iat) * disp_value
            // where disp_value = (t6 + 2*r4r2ij*t8) * zetac6 (the dispersion "strength")
            if (m_dc6dcn_ptr && m_dc6dcn_ptr->size() > 0 &&
                disp.i < m_dc6dcn_ptr->rows() && disp.j < m_dc6dcn_ptr->cols()) {
                double disp_value = disp_sum * disp.zetac6 * m_final_factor;
                m_dEdcn(disp.i) -= (*m_dc6dcn_ptr)(disp.i, disp.j) * disp_value;
                m_dEdcn(disp.j) -= (*m_dc6dcn_ptr)(disp.j, disp.i) * disp_value;
            }
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

        Eigen::Vector3d ri = geom().row(rep.i);
        Eigen::Vector3d rj = geom().row(rep.j);
        Eigen::Vector3d rij_vec = ri - rj;
        double rij = rij_vec.norm() * m_au;

        // HIGH PRIORITY FIX (Feb 2026): Strengthen distance check to prevent gradient Inf/NaN
        // Previous threshold 1e-10 too small, gradient division by rij needs robustness
        if (rij > rep.r_cut || rij < 1e-8) continue;

        // Bonded repulsion formula: E = repab * exp(-α * r^1.5) / r
        // r^1.5 = r * sqrt(r) avoids std::pow()
        double r_1_5 = rij * std::sqrt(rij);
        double exp_term = std::exp(-rep.alpha * r_1_5);
        double base_energy = rep.repab * exp_term / rij;

        double scaled_energy = base_energy * m_final_factor * m_rep_scaling;
        m_rep_energy += scaled_energy;
        m_bonded_rep_energy += scaled_energy;
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

        Eigen::Vector3d ri = geom().row(rep.i);
        Eigen::Vector3d rj = geom().row(rep.j);
        Eigen::Vector3d rij_vec = ri - rj;
        double rij = rij_vec.norm() * m_au;

        // HIGH PRIORITY FIX (Feb 2026): Strengthen distance check to prevent gradient Inf/NaN
        // Previous threshold 1e-10 too small, gradient division by rij needs robustness
        if (rij > rep.r_cut || rij < 1e-8) continue;

        // Non-bonded repulsion formula: E = repab * exp(-α * r^1.5) / r
        // r^1.5 = r * sqrt(r) avoids std::pow()
        double r_1_5 = rij * std::sqrt(rij);
        double exp_term = std::exp(-rep.alpha * r_1_5);
        double base_energy = rep.repab * exp_term / rij;

        double scaled_energy = base_energy * m_final_factor * m_rep_scaling;
        m_rep_energy += scaled_energy;
        m_nonbonded_rep_energy += scaled_energy;
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
    // TERM 1: Pairwise Coulomb interactions (distance-dependent, parallelizable)
    // =========================================================================
    // Claude Generated (Feb 23, 2026): TERM 2+3 (self-energy) and Term 1b (CN chain-rule)
    // moved to parent reduction in ForceField::Calculate() for thread-count independence.
    // Reference: Fortran gfnff_engrad.F90:1670-1680 (single sequential loop)
    double E_interaction = 0.0;

    for (int index = 0; index < m_gfnff_coulombs.size(); ++index) {
        const auto& coul = m_gfnff_coulombs[index];

        Eigen::Vector3d ri = geom().row(coul.i);
        Eigen::Vector3d rj = geom().row(coul.j);
        Eigen::Vector3d rij_vec = ri - rj;
        double rij = rij_vec.norm() * m_au;

        if (rij > coul.r_cut || rij < 1e-10) continue;

        // Use dynamic EEQ charges via pointer or local copy, fall back to static if unavailable or NaN
        double qi = coul.q_i;
        double qj = coul.q_j;
        const bool have_eeq = m_eeq_charges_ptr ? (m_eeq_charges_ptr->size() > 0) : (m_eeq_charges.size() > 0);
        if (have_eeq) {
            qi = eeq_q(coul.i);
            qj = eeq_q(coul.j);
            if (std::isnan(qi) || std::isnan(qj)) {
                qi = coul.q_i;
                qj = coul.q_j;
            }
        }

        // Pairwise: E_pair = q_i * q_j * erf(γ_ij*r) / r
        double gamma_r = coul.gamma_ij * rij;
        double erf_term = std::erf(gamma_r);
        double energy_pair = qi * qj * erf_term / rij;

        E_interaction += energy_pair;
        m_coulomb_energy += energy_pair;

        if (index < 5 && CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info(fmt::format("Coulomb pair {}: i={}, j={}, rij={:.6f} Bohr, gamma_ij={:.6f}, "
                                             "q_i={:.6f}, q_j={:.6f}, erf={:.6f}, E_pair={:.9f} Eh",
                                             index, coul.i, coul.j, rij, coul.gamma_ij,
                                             qi, qj, erf_term, energy_pair));
        }

        if (m_calculate_gradient) {
            const double sqrt_pi = 1.772453850905516;  // √π
            double exp_term = std::exp(-gamma_r * gamma_r);
            double derf_dr = coul.gamma_ij * exp_term * (2.0 / sqrt_pi);
            double dEdr_pair = qi * qj *
                (derf_dr / rij - erf_term / (rij * rij));

            Eigen::Vector3d grad = dEdr_pair * rij_vec / rij;

            m_gradient.row(coul.i) += grad.transpose();
            m_gradient.row(coul.j) -= grad.transpose();
        }
    }

    // Debug summary: pairwise interaction only (TERM 2+3 computed in parent)
    if (CurcumaLogger::get_verbosity() >= 1 && m_gfnff_coulombs.size() > 0) {
        CurcumaLogger::info(fmt::format("  Coulomb pairwise (thread): {:+.12f} Eh", E_interaction));
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
 * FIX (January 18, 2026): Corrected formula to match Fortran reference
 * Formula: outl = 2 / (1 + exp((bacut/radab) × ((r_AH + r_HB)/r_AB - 1)))
 *
 * Reference: gfnff_engrad.F90:1671
 *   expo = (param%hbacut/radab)*(rahprbh/rab-1.d0)
 *   outl = 2.d0/(1.d0+ratio2)
 *
 * @param r_AH Distance A-H (or A-X) in Bohr
 * @param r_HB Distance H-B (or X-B) in Bohr
 * @param r_AB Distance A-B in Bohr
 * @param radab Combined covalent radii (A+B) in Ångström
 * @param bacut Angle cut-off parameter (HB_BACUT or XB_BACUT)
 */
inline double damping_out_of_line(double r_AH, double r_HB, double r_AB, double radab, double bacut)
{
    // Fortran: expo = (param%hbacut/radab)*(rahprbh/rab-1.d0)
    // radab is in Ångström (covalent radii sum)
    double ratio = (r_AH + r_HB) / r_AB;
    double exponent = (bacut / radab) * (ratio - 1.0);

    // Fortran line 1672: early return if exponent too large (avoid overflow)
    if (exponent > 15.0) return 0.0;

    return 2.0 / (1.0 + std::exp(exponent));
}

/**
 * @brief GFN-FF out-of-line damping function for halogen bonds (XB)
 *
 * Reference: gfnff_engrad.F90:3309 - expo = param%xbacut*((rax+rbx)/rab-1.d0)
 * Note: XB does NOT divide xbacut by radab.
 */
inline double damping_out_of_line_xb(double r_AX, double r_XB, double r_AB, double bacut)
{
    double ratio = (r_AX + r_XB) / r_AB;
    double exponent = bacut * (ratio - 1.0);

    // Fortran line 3310: early return if exponent too large (avoid overflow)
    if (exponent > 15.0) return 0.0;

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

    //CurcumaLogger::error(fmt::format("Thread {} calculating {} hydrogen bonds", m_thread, m_gfnff_hbonds.size()));

    for (const auto& hb : m_gfnff_hbonds) {
        // Get atom positions
        Eigen::Vector3d pos_A = geom().row(hb.i).transpose();
        Eigen::Vector3d pos_H = geom().row(hb.j).transpose();
        Eigen::Vector3d pos_B = geom().row(hb.k).transpose();

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

        // --- Donor-Acceptor Strength Term B_AH and acidity mix ---
        // Reference: gfnff_engrad.F90:1704-1710 (eg1) / 1850-1855 (eg2new)
        double r_AH_4 = r_AH * r_AH * r_AH * r_AH;
        double r_HB_4 = r_HB * r_HB * r_HB * r_HB;
        double denom_DA = 1.0 / (r_AH_4 + r_HB_4);

        // --- Charge Scaling Factors ---
        // CRITICAL FIX (January 17, 2026): A and B use NEGATIVE charge (Fortran: exp(-hbst*q))
        // Reference: gfnff_engrad.F90:1693-1699
        double Q_H = charge_scaling(hb.q_H, HB_ST, HB_SF);    // qh = exp(st*q_H) / (exp(st*q_H) + sf)
        double Q_A = charge_scaling(-hb.q_A, HB_ST, HB_SF);   // qa = exp(-st*q_A) / ...
        double Q_B = charge_scaling(-hb.q_B, HB_ST, HB_SF);   // qb = exp(-st*q_B) / ...

        // Weighted Mix for Basicity (bas) and Acidity (aci)
        // eg1: bas = (qa*bas_A*r_AH^4 + qb*bas_B*r_HB^4) * denom
        // eg1: aci = (aci_B*r_AH^4 + aci_A*r_HB^4) * denom
        double bas = (Q_A * hb.basicity_A * r_AH_4 + Q_B * hb.basicity_B * r_HB_4) * denom_DA;
        double aci = (hb.acidity_B * r_AH_4 + hb.acidity_A * r_HB_4) * denom_DA;

        // --- Combined Damping R_damp ---
        int elem_A = m_atom_types[hb.i];
        int elem_B = m_atom_types[hb.k];
        double r_vdw_AB = covalent_radii[elem_A - 1] + covalent_radii[elem_B - 1]; // [Å] (Matches Fortran param%rad)

        // Important: damping_short_range expects scut*r_vdw in unit [Bohr^2] internally
        // In GFN-FF, HB_SCUT=22.0 is calibrated for Angstrom radii vs Bohr distance
        double damp_short = damping_short_range(r_AB, r_vdw_AB, HB_SCUT, HB_ALP);
        double damp_long = damping_long_range(r_AB, HB_LONGCUT, HB_ALP);
        double damp_env = damp_short * damp_long;

        // Out-of-line damping (A-H...B)
        // Reference: gfnff_engrad.F90:1671 - expo = (hbacut/radab)*(rahprbh/rab-1)
        // radab is sum of radii in Ångström!
        double damp_outl = damping_out_of_line(r_AH, r_HB, r_AB, r_vdw_AB, HB_BACUT);

        // Neighbor environmental damping (outl_nb_tot) - Case 2/3
        double outl_nb_tot = 1.0;
        if (hb.case_type >= 2) {
            // Reference: gfnff_engrad.F90:1875-1881
            // expo_nb = (hbnbcut/radab) * ((r_Anb + r_Bnb)/r_AB - 1)
            // outl_nb = 2/(1+exp(-expo_nb)) - 1.0

            double hbnbcut_save = (elem_B == 7 && hb.neighbors_B.size() == 1) ? 2.0 : HB_NBCUT;

            for (int nb : hb.neighbors_B) {
                Eigen::Vector3d pos_nb = geom().row(nb).transpose();
                double r_Anb = (pos_nb - pos_A).norm() * m_au;
                double r_Bnb = (pos_nb - pos_B).norm() * m_au;

                double expo_nb = (hbnbcut_save / r_vdw_AB) * ((r_Anb + r_Bnb) / r_AB - 1.0);
                double ratio_nb = std::exp(-expo_nb);
                outl_nb_tot *= (2.0 / (1.0 + ratio_nb) - 1.0);
            }
        }

        // Distance Damping Factor (rdamp)
        double rdamp;
        if (hb.case_type >= 2) {
            // Reference: gfnff_engrad.F90:1894-1896
            // rdamp = damp_env * (1.8/rbh³ - 0.8/rab³)
            rdamp = damp_env * (1.8 / (r_HB * r_HB * r_HB) - 0.8 / (r_AB * r_AB * r_AB));
        } else {
            // Reference: gfnff_engrad.F90:1686
            // rdamp = damp_env / (rab*rab2) = damp_env / rab³
            rdamp = damp_env / (r_AB * r_AB * r_AB);
        }

        // --- Case 4: N heteroaromatic virtual lone pair (eg2_rnr) ---
        // Claude Generated (Feb 25, 2026): Virtual LP out-of-line damping
        // Reference: gfnff_engrad.F90:2153-2219 (abhgfnff_eg2_rnr)
        double outl_lp = 1.0;
        Eigen::Vector3d lp_pos = pos_B;  // fallback if no LP computed
        double lp_dist = 0.0;
        int nbb_lp = 0;
        Eigen::Vector3d lp_vector = Eigen::Vector3d::Zero();

        if (hb.case_type == 4) {
            // Compute virtual lone pair position for sp2 N with 2 neighbors
            // lp_dist = 0.50 - 0.018 * repz(at(B))
            int z_B = m_atom_types[hb.k];
            double repz_B = (z_B >= 1 && z_B <= static_cast<int>(GFNFFParameters::repz.size()))
                          ? GFNFFParameters::repz[z_B - 1] : 1.0;
            lp_dist = 0.50 - 0.018 * repz_B;
            static constexpr double HBLPCUT = 56.0;  // Fortran: hblpcut = 56

            // Sum bond vectors from B to neighbors (in Bohr)
            nbb_lp = static_cast<int>(hb.neighbors_B.size());
            for (int nb_idx : hb.neighbors_B) {
                Eigen::Vector3d pos_nb = geom().row(nb_idx).transpose();
                Eigen::Vector3d bond_vec = (pos_nb - pos_B) * m_au;  // B→nb in Bohr
                lp_vector += bond_vec;
            }
            double vnorm = lp_vector.norm();

            if (vnorm > 1e-10 && nbb_lp > 0) {
                // LP sits opposite to the average bond direction
                lp_pos = pos_B * m_au + (-lp_dist) * (lp_vector / vnorm);
                // Note: lp_pos is now in Bohr, but we need it relative to pos_B in Bohr

                // Compute LP out-of-line damping
                double rblp = lp_dist;  // |B - LP| = lp_dist by construction
                Eigen::Vector3d r_A_bohr = pos_A * m_au;
                double ralp = (r_A_bohr - lp_pos).norm();
                double ralpprblp = ralp + rblp + 1e-12;

                double expo_lp_val = (HBLPCUT / r_vdw_AB) * (ralpprblp / r_AB - 1.0);
                double ratio2_lp = std::exp(expo_lp_val);
                outl_lp = 2.0 / (1.0 + ratio2_lp);
            } else {
                // Degenerate case (vnorm~0): LP at B, no damping
                nbb_lp = 0;
                outl_lp = 1.0;
            }
        }

        // Charge scaling and out-of-line scaling (qhoutl)
        double qhoutl = Q_H * damp_outl * outl_nb_tot * outl_lp;

        // --- Case 3 eangl/etors (Carbonyl/Nitro: H...B=C<D) ---
        // Claude Generated (Feb 24, 2026): Full Fortran abhgfnff_eg3 port
        // Reference: gfnff_engrad.F90:2545-2583
        // eangl = angle bending multiplier for H...B=C angle
        // etors = product of torsion multipliers for D-B-C-H dihedral angles
        double eangl = 1.0;
        double etors = 1.0;
        // gangl[atom_idx] = gradient of eangl at each atom (3-body: H=ll, B=jj, C=kk)
        // Per-torsion energy+gradient for product rule application
        struct TorsGrad {
            double energy;
            Eigen::Matrix<double, 4, 3> grad;  // rows: D,B,C,H
            int D_idx;  // atom index of D
        };
        std::vector<TorsGrad> tors_data;
        Eigen::Matrix<double, 3, 3> gangl_3body = Eigen::Matrix<double, 3, 3>::Zero();  // rows: jj(B), kk(C), ll(H)
        int jj_idx = -1, kk_idx = -1, ll_idx = -1;  // B, C, H atom indices for gangl

        if (hb.case_type == 3 && hb.acceptor_parent_index != -1) {
            int B_idx = hb.k;
            int C_idx = hb.acceptor_parent_index;
            int H_idx = hb.j;
            jj_idx = B_idx;
            kk_idx = C_idx;
            ll_idx = H_idx;

            // --- eangl: angle bending H...B=C ---
            // Reference: gfnff_engrad.F90:2574-2583 (egbend_nci_mul)
            // phi0 = 120°, fc = 1 - bend_hb, kijk = fc / (1 - cos(120°))²
            {
                double c0 = 120.0 * M_PI / 180.0;
                double fc_bend = 1.0 - BEND_HB;
                double kijk = fc_bend / ((std::cos(0.0) - std::cos(c0)) * (std::cos(0.0) - std::cos(c0)));

                // Atoms: j=B (center), i=C, k=H (Fortran convention: egbend_nci_mul(jj,kk,ll,...))
                // jj=B, kk=C, ll=H → center=B(jj), ends=C(kk) and H(ll)
                // In Fortran egbend_nci_mul(j,i,k,...): j=center, i,k=ends
                // So: center=jj=B, i=kk=C, k=ll=H
                Eigen::Vector3d va = geom().row(C_idx).transpose();  // atom i = C
                Eigen::Vector3d vb = geom().row(B_idx).transpose();  // atom j = B (center)
                Eigen::Vector3d vc = geom().row(H_idx).transpose();  // atom k = H

                Eigen::Vector3d vab = (va - vb) * m_au;
                Eigen::Vector3d vcb = (vc - vb) * m_au;
                double rab2_a = vab.squaredNorm();
                double rcb2_a = vcb.squaredNorm();
                Eigen::Vector3d vp = vcb.cross(vab);
                double rp = vp.norm() + 1e-14;

                double cosa = vab.dot(vcb) / (std::sqrt(rab2_a) * std::sqrt(rcb2_a) + 1e-14);
                cosa = std::clamp(cosa, -1.0, 1.0);
                double theta_a = std::acos(cosa);

                double ea, deddt;
                if (M_PI - c0 < 1e-6) {  // linear
                    double dt = theta_a - c0;
                    ea = kijk * dt * dt;
                    deddt = 2.0 * kijk * dt;
                } else {
                    ea = kijk * (cosa - std::cos(c0)) * (cosa - std::cos(c0));
                    deddt = 2.0 * kijk * std::sin(theta_a) * (std::cos(c0) - cosa);
                }
                eangl = 1.0 - ea;

                // Gradient of eangl: returns g(3,3) = [center(B), end1(C), end2(H)]
                Eigen::Vector3d deda_v = vab.cross(vp) * (-deddt / (rab2_a * rp));
                Eigen::Vector3d dedc_v = vcb.cross(vp) * (deddt / (rcb2_a * rp));
                Eigen::Vector3d dedb_v = deda_v + dedc_v;

                // gangl: row 0 = B (center/jj), row 1 = C (kk), row 2 = H (ll)
                // Fortran: g(1:3,1) = dedb, g(1:3,2) = -deda, g(1:3,3) = -dedc
                gangl_3body.row(0) = dedb_v.transpose();   // B (center)
                gangl_3body.row(1) = -deda_v.transpose();  // C (end1)
                gangl_3body.row(2) = -dedc_v.transpose();  // H (end2)
            }

            // --- etors: product of torsion terms D-B-C-H ---
            // Reference: gfnff_engrad.F90:2529-2557 (egtors_nci_mul)
            // For each D = neighbor of C (excluding B): torsion D-B-C-H
            // Parameters: rn=2, phi0=pi/2, tshift=tors_hb, fc=(1-tshift)/2
            for (int D_idx : hb.neighbors_C) {
                int ii = D_idx;     // tlist(1) = D
                int jj_t = B_idx;   // tlist(2) = B
                int kk_t = C_idx;   // tlist(3) = C
                int ll_t = H_idx;   // tlist(4) = H

                // Compute dihedral angle and gradient
                Matrix dihedral_grad;
                bool calc_grad = m_calculate_gradient;
                double phi = GFNFF_Geometry::calculateDihedralAngle(
                    geom().row(ii).transpose() * m_au,
                    geom().row(jj_t).transpose() * m_au,
                    geom().row(kk_t).transpose() * m_au,
                    geom().row(ll_t).transpose() * m_au,
                    dihedral_grad, calc_grad);

                // Fortran: fc=(1-tshift)/2, phi0=pi/2, rn=2
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
                if (calc_grad && dihedral_grad.rows() == 4) {
                    // g(1:3,atom) = dij * dphi/dr(atom)
                    for (int a = 0; a < 4; ++a) {
                        tg.grad.row(a) = dij * dihedral_grad.row(a);
                    }
                }
                tors_data.push_back(tg);
            }

            // etors = product of all torsion energies
            etors = 1.0;
            for (const auto& tg : tors_data) {
                etors *= tg.energy;
            }
        }

        // Global Scaling and Final Energy
        double global_scale = 1.0;
        if (hb.case_type == 2 || hb.case_type == 4) global_scale = XHACI_GLOBABH;
        else if (hb.case_type == 3) global_scale = XHACI_COH;

        double E_HB;
        if (hb.case_type >= 2) {
            // Const part: acidity_A * basicity_B * qa * qb * global_scale
            // Reference: gfnff_engrad.F90:2603 (eg3) / 1913 (eg2new)
            double const_val = hb.acidity_A * hb.basicity_B * Q_A * Q_B * global_scale;
            // Case 3: energy = -rdamp * qhoutl * eangl * etors * const
            // Case 2: eangl = etors = 1 (no angle/torsion terms)
            E_HB = -rdamp * qhoutl * const_val * eangl * etors;
        } else {
            // Case 1: energy = bas * (-aci * rdamp * qhoutl)
            // Reference: gfnff_engrad.F90:1798-1799
            E_HB = bas * (-aci) * rdamp * qhoutl;
        }

        m_energy_hbond += E_HB * m_final_factor;

        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info(fmt::format(
                "  HB({}-{}-{}): r_AB={:.3f} Bohr, E={:.6e} Eh",
                hb.i, hb.j, hb.k, r_AB, E_HB * m_final_factor));
            CurcumaLogger::info(fmt::format(
                "    DEBUG: bas={:.6f}, aci={:.6f}, rdamp={:.6e}, qhoutl={:.6f}, case={}",
                bas, aci, rdamp, qhoutl, hb.case_type));
        }

        // ========== ANALYTICAL GRADIENT CALCULATION ==========
        // REWRITTEN (Feb 19, 2026): Direct translation from Fortran gfnff_engrad.F90
        //   Case 1: abhgfnff_eg1() lines 1692-1853
        //   Case 2: abhgfnff_eg2new() lines 1856-2087
        // Previous implementation had 3 bugs:
        //   1. Short damping derivative sign error (positive → negative)
        //   2. Long damping derivative magnitude error (factor rab²/longcut)
        //   3. Missing neighbor out-of-line gradient terms (case >= 2)
        if (m_calculate_gradient) {
            // Fortran distance vectors: drab = xyz(A)-xyz(B), drah = xyz(A)-xyz(H), drbh = xyz(B)-xyz(H)
            // Curcuma vectors: r_AH_vec = H-A, r_HB_vec = B-H, r_AB_vec = B-A
            Eigen::Vector3d drab = -r_AB_vec;  // A - B (Fortran convention)
            Eigen::Vector3d drah = -r_AH_vec;  // A - H (Fortran convention)
            Eigen::Vector3d drbh =  r_HB_vec;  // B - H (same as Fortran)

            double rab2 = r_AB * r_AB;
            double rbh2 = r_HB * r_HB;
            double rah2 = r_AH * r_AH;
            double rahprbh = r_AH + r_HB + 1e-12;

            // Recompute damping intermediates for gradient
            // Fortran: ratio1 = (rab2/hblongcut)^hbalp, ratio3 = (shortcut/rab2)^hbalp
            double ratio1 = std::pow(rab2 / HB_LONGCUT, HB_ALP);
            double shortcut = HB_SCUT * r_vdw_AB;
            double ratio3 = std::pow(shortcut / rab2, HB_ALP);
            // ddamp = rab * d(ln(damp))/d(rab), exactly Fortran gfnff_engrad.F90:1969
            double ddamp = (-2.0 * HB_ALP * ratio1 / (1.0 + ratio1))
                         + ( 2.0 * HB_ALP * ratio3 / (1.0 + ratio3));

            // Recompute out-of-line intermediates for gradient
            // Fortran: expo = (hbacut/radab)*(rahprbh/rab - 1)
            double expo = (HB_BACUT / r_vdw_AB) * (rahprbh / r_AB - 1.0);

            // Claude Generated (Feb 25, 2026): Match Fortran early return (gfnff_engrad.F90:1781)
            // When expo > 15, outl → 0 and energy → 0; Fortran returns before gradient.
            // Without this guard, aterm (which lacks outl factor) produces huge gradients.
            if (expo > 15.0) continue;

            double ratio2 = std::exp(expo);  // Fortran ratio2 = exp(expo) for outl

            Eigen::Vector3d ga = Eigen::Vector3d::Zero();
            Eigen::Vector3d gb = Eigen::Vector3d::Zero();
            Eigen::Vector3d gh = Eigen::Vector3d::Zero();

            if (hb.case_type >= 2) {
                // ===== Case 2/3: abhgfnff_eg2new gradient (lines 1997-2086) =====
                double p_bh = 1.8;   // 1 + hbabmix (hbabmix = 0.8)
                double p_ab = -0.8;  // -hbabmix
                double rbhdamp = damp_env * p_bh / (rbh2 * r_HB);
                double rabdamp = damp_env * p_ab / (rab2 * r_AB);

                double const_val = hb.acidity_A * hb.basicity_B * Q_A * Q_B * global_scale;
                // Claude Generated (Feb 25, 2026): Include eangl*etors*outl_lp in gradient coefficients
                // Reference: Fortran gfnff_engrad.F90:2275-2278 (eg2_rnr), 2613-2617 (eg3)
                // For Case 2: eangl=etors=outl_lp=1, so these reduce to the previous form
                // For Case 4: eangl=etors=1, outl_lp is the LP damping factor
                double dterm  = -qhoutl * eangl * etors * const_val;
                double aterm  = -rdamp * Q_H * outl_nb_tot * outl_lp * eangl * etors * const_val;
                double nbterm = -rdamp * Q_H * damp_outl * outl_lp * eangl * etors * const_val;

                // --- Damping part: rab (Fortran lines 2008-2012) ---
                double gi = ((rabdamp + rbhdamp) * ddamp - 3.0 * rabdamp) / rab2;
                gi *= dterm;
                Eigen::Vector3d dg = gi * drab;
                ga = dg;
                gb = -dg;

                // --- Damping part: rbh (Fortran lines 2016-2020) ---
                gi = -3.0 * rbhdamp / rbh2;
                gi *= dterm;
                dg = gi * drbh;
                gb += dg;
                gh = -dg;

                // --- Out-of-line: rab (Fortran lines 2026-2030) ---
                double tmp1 = -2.0 * aterm * ratio2 * expo
                            / ((1.0 + ratio2) * (1.0 + ratio2))
                            / (rahprbh - r_AB);
                gi = -tmp1 * rahprbh / rab2;
                dg = gi * drab;
                ga += dg;
                gb -= dg;

                // --- Out-of-line: rah, rbh (Fortran lines 2033-2040) ---
                gi = tmp1 / r_AH;
                Eigen::Vector3d dga_outl = gi * drah;
                ga += dga_outl;
                gi = tmp1 / r_HB;
                Eigen::Vector3d dgb_outl = gi * drbh;
                gb += dgb_outl;
                Eigen::Vector3d dgh_outl = -dga_outl - dgb_outl;
                gh += dgh_outl;

                // --- Neighbor out-of-line: rab + ranb/rbnb (Fortran lines 2047-2068) ---
                // Previously MISSING - causes significant gradient error
                double hbnbcut_save = (elem_B == 7 && hb.neighbors_B.size() == 1) ? 2.0 : HB_NBCUT;
                for (size_t idx = 0; idx < hb.neighbors_B.size(); ++idx) {
                    int nb = hb.neighbors_B[idx];
                    Eigen::Vector3d pos_nb = geom().row(nb).transpose();
                    Eigen::Vector3d dranb = pos_A - pos_nb;  // A - nb (Fortran convention)
                    Eigen::Vector3d drbnb = pos_B - pos_nb;  // B - nb (Fortran convention)
                    double ranb = dranb.norm();
                    double rbnb = drbnb.norm();
                    double ranbprbnb = ranb + rbnb + 1e-12;

                    // Recompute this neighbor's outl_nb and expo_nb
                    double expo_nb_i = (hbnbcut_save / r_vdw_AB) * (ranbprbnb / r_AB - 1.0);
                    double ratio2_nb_i = std::exp(-expo_nb_i);
                    double outl_nb_i = 2.0 / (1.0 + ratio2_nb_i) - 1.0;

                    // Product of all OTHER outl_nb values (Fortran: product(outl_nb, mask))
                    double outl_nb_others = 1.0;
                    if (std::abs(outl_nb_i) > 1e-12) {
                        outl_nb_others = outl_nb_tot / outl_nb_i;
                    }

                    // Fortran: tmp2 = 2*nbterm*product(outl_nb,mask)*ratio2_nb*expo_nb /
                    //                  (1+ratio2_nb)^2 / (ranbprbnb - rab)
                    double tmp2 = 2.0 * nbterm * outl_nb_others * ratio2_nb_i * expo_nb_i
                                / ((1.0 + ratio2_nb_i) * (1.0 + ratio2_nb_i))
                                / (ranbprbnb - r_AB);

                    // rab contribution (Fortran lines 2051-2054)
                    double gi_nb = -tmp2 * ranbprbnb / rab2;
                    dg = gi_nb * drab;
                    ga += dg;
                    gb -= dg;

                    // ranb, rbnb contributions (Fortran lines 2060-2067)
                    gi_nb = tmp2 / ranb;
                    Eigen::Vector3d dga_nb = gi_nb * dranb;
                    ga += dga_nb;
                    gi_nb = tmp2 / rbnb;
                    Eigen::Vector3d dgb_nb = gi_nb * drbnb;
                    gb += dgb_nb;
                    Eigen::Vector3d dgnb = -dga_nb - dgb_nb;
                    m_gradient.row(nb) += dgnb.transpose() * m_final_factor;
                }

                // --- Case 4 only: LP out-of-line gradient (eg2_rnr) ---
                // Claude Generated (Feb 25, 2026): Fortran abhgfnff_eg2_rnr lines 2316-2387
                if (hb.case_type == 4 && nbb_lp > 0) {
                    // lpterm = -rdamp * qh * outl * outl_nb_tot * const (d/d(outl_lp))
                    double lpterm = -rdamp * Q_H * damp_outl * outl_nb_tot * const_val;

                    // LP coordinates relative to atom positions (all in Bohr)
                    Eigen::Vector3d r_A_bohr = pos_A * m_au;
                    Eigen::Vector3d r_B_bohr = pos_B * m_au;
                    double ralp = (r_A_bohr - lp_pos).norm();
                    double rblp = lp_dist;
                    double ralpprblp = ralp + rblp + 1e-12;
                    static constexpr double HBLPCUT = 56.0;
                    double expo_lp_val = (HBLPCUT / r_vdw_AB) * (ralpprblp / r_AB - 1.0);
                    double ratio2_lp_val = std::exp(expo_lp_val);

                    // --- LP out-of-line: rab (Fortran lines 2319-2324) ---
                    double tmp3 = -2.0 * lpterm * ratio2_lp_val * expo_lp_val
                                / ((1.0 + ratio2_lp_val) * (1.0 + ratio2_lp_val))
                                / (ralpprblp - r_AB);
                    double gi_lp = -tmp3 * ralpprblp / rab2;
                    Eigen::Vector3d dg_lp = gi_lp * drab;
                    ga += dg_lp;
                    gb -= dg_lp;

                    // --- LP out-of-line: ralp (Fortran lines 2326-2329) ---
                    Eigen::Vector3d dralp = (r_A_bohr - lp_pos);  // A - LP
                    gi_lp = tmp3 / (ralp + 1e-12);
                    Eigen::Vector3d dga_lp = gi_lp * dralp;
                    ga += dga_lp;

                    // --- LP out-of-line: rblp (Fortran lines 2330-2333) ---
                    // Fortran: gb += -dga_lp (not dgb_lp!), glp = -dga_lp
                    // This is because LP is a virtual point derived from B's neighbors
                    Eigen::Vector3d drblp = (r_B_bohr - lp_pos);  // B - LP
                    // Fortran line 2332: gb = gb - dga (note: uses dga, not dgb!)
                    gb -= dga_lp;
                    Eigen::Vector3d glp = -dga_lp;  // Gradient on LP

                    // --- LP neighbor chain rule (Fortran lines 2336-2342, 2383-2387) ---
                    // LP = B - lp_dist * (vector / |vector|)
                    // d(LP)/d(nb_pos) propagated to B's neighbors
                    double vnorm = lp_vector.norm();
                    if (vnorm > 1e-10) {
                        // Jacobian: gii(i,j) = d(LP_i)/d(nb_j) = -lp_dist*nbb*(delta_ij/vnorm + v_i*v_j/vnorm³)
                        // gnb_lp = gii * glp (chain rule)
                        Eigen::Matrix3d gii = Eigen::Matrix3d::Zero();
                        for (int col = 0; col < 3; ++col) {
                            Eigen::Vector3d unit_vec = Eigen::Vector3d::Zero();
                            unit_vec(col) = -1.0;
                            gii.col(col) = -lp_dist * static_cast<double>(nbb_lp)
                                         * (unit_vec / vnorm + lp_vector * lp_vector(col) / std::pow(vnorm, 3.0));
                        }
                        Eigen::Vector3d gnb_lp = gii * glp;

                        // gdr(B) += gnb_lp; gdr(nb_i) += gnb(:,i) - gnb_lp/nbb
                        gb += gnb_lp;
                        Eigen::Vector3d gnb_lp_share = gnb_lp / static_cast<double>(nbb_lp);
                        for (int nb_idx : hb.neighbors_B) {
                            // Note: gnb(nb) gradient from neighbor out-of-line already applied above
                            // This adds the LP chain rule: -gnb_lp / nbb per neighbor
                            m_gradient.row(nb_idx) -= gnb_lp_share.transpose() * m_final_factor;
                        }
                    }
                }

                // --- Case 3 only: angle bending and torsion gradient contributions ---
                // Claude Generated (Feb 24, 2026): Fortran abhgfnff_eg3 lines 2694-2717
                // bterm = -rdamp * qhoutl * etors * const  (d(E)/d(eangl) chain rule)
                // tterm = -rdamp * qhoutl * eangl * const  (d(E)/d(etors) chain rule)
                if (hb.case_type == 3 && jj_idx >= 0) {
                    double bterm_c3 = -rdamp * qhoutl * etors * const_val;
                    double tterm_c3 = -rdamp * qhoutl * eangl * const_val;

                    // Angle bending gradient: gangl_3body stores d(eangl)/d(atom)
                    // rows: 0=B(jj_idx), 1=C(kk_idx), 2=H(ll_idx)
                    m_gradient.row(jj_idx) += bterm_c3 * gangl_3body.row(0) * m_final_factor;
                    m_gradient.row(kk_idx) += bterm_c3 * gangl_3body.row(1) * m_final_factor;
                    m_gradient.row(ll_idx) += bterm_c3 * gangl_3body.row(2) * m_final_factor;

                    // Torsion gradient: product rule d(etors)/dr_atom = sum_k [prod_other_k * dtors_k/dr]
                    for (size_t k = 0; k < tors_data.size(); ++k) {
                        // product of all OTHER torsion energies (excluding torsion k)
                        double factor_k = (std::abs(tors_data[k].energy) > 1e-12)
                                        ? etors / tors_data[k].energy : 0.0;
                        double t_k = factor_k * tterm_c3;

                        // Gradient rows: 0=D, 1=B, 2=C, 3=H
                        m_gradient.row(tors_data[k].D_idx) += t_k * tors_data[k].grad.row(0) * m_final_factor;
                        m_gradient.row(jj_idx)             += t_k * tors_data[k].grad.row(1) * m_final_factor;
                        m_gradient.row(kk_idx)             += t_k * tors_data[k].grad.row(2) * m_final_factor;
                        m_gradient.row(ll_idx)             += t_k * tors_data[k].grad.row(3) * m_final_factor;
                    }
                }
            } else {
                // ===== Case 1: abhgfnff_eg1 gradient (lines 1795-1851) =====
                double caa = Q_A * hb.basicity_A;   // Fortran: qa*ca(1)
                double cbb = Q_B * hb.basicity_B;   // Fortran: qb*cb(1)

                // Fortran terms
                double rterm = -aci * rdamp * qhoutl;         // -aci*rdamp*qhoutl
                double dterm = -aci * bas * qhoutl;           // -aci*bas*qhoutl
                double sterm = -rdamp * bas * qhoutl;         // -rdamp*bas*qhoutl
                double aterm = -aci * bas * rdamp * Q_H;      // -aci*bas*rdamp*qh

                double denom_val = 1.0 / (r_AH_4 + r_HB_4);
                double tmp = denom_val * denom_val * 4.0;
                double dd24a = rah2 * r_HB_4 * tmp;  // rah²*rbh⁴*4/denom²
                double dd24b = rbh2 * r_AH_4 * tmp;  // rbh²*rah⁴*4/denom²

                // --- Donor-acceptor part: bas (Fortran lines 1809-1813) ---
                double gi = (caa - cbb) * dd24a * rterm;
                ga = gi * drah;
                gi = (cbb - caa) * dd24b * rterm;
                gb = gi * drbh;
                gh = -ga - gb;

                // --- Donor-acceptor part: aci (Fortran lines 1816-1825) ---
                gi = (hb.acidity_B - hb.acidity_A) * dd24a;
                Eigen::Vector3d dga_aci = gi * drah * sterm;
                ga += dga_aci;
                gi = (hb.acidity_A - hb.acidity_B) * dd24b;
                Eigen::Vector3d dgb_aci = gi * drbh * sterm;
                gb += dgb_aci;
                Eigen::Vector3d dgh_aci = -dga_aci - dgb_aci;
                gh += dgh_aci;

                // --- Damping part: rab (Fortran line 1828) ---
                // eg1: rdamp = damp/rab³, so d(rdamp)/d(rab) = damp*(ddamp-3)/rab⁴
                gi = rdamp * (ddamp - 3.0) / rab2;
                Eigen::Vector3d dg = gi * drab * dterm;
                ga += dg;
                gb -= dg;

                // --- Out-of-line: rab (Fortran line 1834) ---
                gi = aterm * 2.0 * ratio2 * expo * rahprbh
                   / ((1.0 + ratio2) * (1.0 + ratio2))
                   / (rahprbh - r_AB) / rab2;
                dg = gi * drab;
                ga += dg;
                gb -= dg;

                // --- Out-of-line: rah, rbh (Fortran lines 1840-1846) ---
                double tmp_outl = -2.0 * aterm * ratio2 * expo
                                / ((1.0 + ratio2) * (1.0 + ratio2))
                                / (rahprbh - r_AB);
                Eigen::Vector3d dga_outl = drah * tmp_outl / r_AH;
                ga += dga_outl;
                Eigen::Vector3d dgb_outl = drbh * tmp_outl / r_HB;
                gb += dgb_outl;
                Eigen::Vector3d dgh_outl = -dga_outl - dgb_outl;
                gh += dgh_outl;
            }

            // Accumulate gradients: A=hb.i, B=hb.k, H=hb.j
            m_gradient.row(hb.i) += ga.transpose() * m_final_factor;
            m_gradient.row(hb.k) += gb.transpose() * m_final_factor;
            m_gradient.row(hb.j) += gh.transpose() * m_final_factor;

            if (CurcumaLogger::get_verbosity() >= 3) {
                CurcumaLogger::info(fmt::format(
                    "    HB Gradient: |∇A|={:.4f} |∇H|={:.4f} |∇B|={:.4f} Eh/Bohr",
                    ga.norm(), gb.norm(), gh.norm()));
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
        Eigen::Vector3d pos_A = geom().row(xb.i).transpose();
        Eigen::Vector3d pos_X = geom().row(xb.j).transpose();
        Eigen::Vector3d pos_B = geom().row(xb.k).transpose();

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
        // CRITICAL FIX (January 17, 2026): B uses NEGATIVE charge (Fortran: exp(-xbst*q))
        // Reference: gfnff_engrad.F90:3210-3218
        // XB uses different parameters: XB_ST=15, XB_SF=0.03 (much weaker offset)
        double Q_X = charge_scaling(xb.q_X, XB_ST, XB_SF);    // exp(+st*q_X)
        double Q_B = charge_scaling(-xb.q_B, XB_ST, XB_SF);   // exp(-st*q_B) ← NEGATIVE!

        // --- Combined Damping R_damp ---
        // CRITICAL FIX (March 10, 2026): XB uses radii of donor A and acceptor B (not halogen X and B)
        // Reference: gfnff_engrad.F90:3319 - shortcut = xbscut * (rad(A) + rad(B))
        int elem_A = m_atom_types[xb.i];
        int elem_B_xb = m_atom_types[xb.k];
        double r_vdw_AB_xb = covalent_radii[elem_A - 1] + covalent_radii[elem_B_xb - 1];  // Ångström (Matches Fortran)

        // XB uses different cutoff parameters
        double damp_short = damping_short_range(r_XB, r_vdw_AB_xb, XB_SCUT, HB_ALP);
        double damp_long = damping_long_range(r_XB, HB_LONGCUT_XB, HB_ALP);
        // FIX (March 10, 2026): XB does NOT divide xbacut by radab
        // Reference: gfnff_engrad.F90:3309 - expo = xbacut * ((rax+rbx)/rab - 1)
        double damp_outl = damping_out_of_line_xb(r_AX, r_XB, r_AB, XB_BACUT);

        double R_damp = damp_short * damp_long * damp_outl / (r_XB * r_XB * r_XB);

        // --- Total XB Energy ---
        // C_B = 1.0 (fixed), C_X = xbaci parameter
        double C_B = 1.0;
        double C_X = xb.acidity_X;

        double E_XB = -R_damp * C_B * Q_B * C_X * Q_X;

        // Reference: gfnff_engrad.F90:rbxgfnff_eg - No sigmoid in Fortran
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
            double ratio_short = XB_SCUT * r_vdw_AB_xb / (r_XB * r_XB);
            double damp_short_term = std::pow(ratio_short, HB_ALP);
            double ddamp_short_dr = -2.0 * HB_ALP * damp_short * damp_short_term / (r_XB * (1.0 + damp_short_term));

            // Long-range damping: damp_l = 1 / (1 + (r² / longcut)^alp)
            double ratio_long = (r_XB * r_XB) / HB_LONGCUT_XB;
            double damp_long_term = std::pow(ratio_long, HB_ALP);
            double ddamp_long_dr = -2.0 * HB_ALP * r_XB * damp_long * damp_long_term / (HB_LONGCUT_XB * (1.0 + damp_long_term));

            // Out-of-line damping derivatives
            // FIX (March 10, 2026): XB does NOT divide xbacut by radab
            double ratio = (r_AX + r_XB) / r_AB;
            double scale_outl_xb = XB_BACUT;
            double exponent = scale_outl_xb * (ratio - 1.0);
            double exp_term = std::exp(exponent);
            double denom_outl = 1.0 + exp_term;

            // Derivatives of out-of-line damping w.r.t. each distance
            double ddamp_outl_drAX = -2.0 * exp_term * scale_outl_xb / (r_AB * denom_outl * denom_outl);
            double ddamp_outl_drXB = ddamp_outl_drAX;  // Same formula
            double ddamp_outl_drAB = 2.0 * exp_term * scale_outl_xb * (r_AX + r_XB) / (r_AB * r_AB * denom_outl * denom_outl);

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

            // Claude Generated (Feb 15, 2026): Sigmoid strength gradient contribution for XB
            // REMOVED: The sigmoid derivative adds an unphysical force that destabilizes MD
            // The XB energy already has natural distance dependence through rdamp
            // The sigmoid is only for smooth energy transition at cutoff, NOT for forces

            // --- 3. Chain Rule: Position Vector Derivatives ---
            Eigen::Vector3d grad_rAX_unit = r_AX_vec / r_AX;  // Direction A → X
            Eigen::Vector3d grad_rXB_unit = r_XB_vec / r_XB;  // Direction X → B
            Eigen::Vector3d grad_rAB_unit = r_AB_vec / r_AB;  // Direction A → B

            // --- 4. Accumulate Gradients on Each Atom ---
            // Reference: gfnff_engrad.F90 - No sigmoid in Fortran, direct gradient scaling

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

        Eigen::Vector3d ri = geom().row(disp.i);
        Eigen::Vector3d rj = geom().row(disp.j);
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
        m_dispersion_energy += energy; // Claude Generated: Ensure dispersion energy is tracked for D3 fallback

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

// Claude Generated 2025: Native D4/GFN-FF Dispersion with MODIFIED Becke-Johnson damping
// UPDATED (January 25, 2026): Fix to match XTB 6.6.1 reference
void ForceFieldThread::CalculateD4DispersionContribution()
{
    /**
     * @brief GFN-FF Dispersion with MODIFIED Becke-Johnson damping
     *
     * CLAUDE GENERATED (January 25, 2026): Fix to match XTB 6.6.1 reference
     *
     * Reference: gfnff_gdisp0.f90:365-377
     *
     * GFN-FF uses a MODIFIED BJ damping formula (NOT standard D3/D4):
     *   E = -0.5 * C6 * (t6 + 2*r4r2ij*t8)
     * where:
     *   t6 = 1/(r^6 + R0^6)
     *   t8 = 1/(r^8 + R0^8)
     *   r4r2ij = 3 * sqrtZr4r2_i * sqrtZr4r2_j  (implicit C8/C6 ratio)
     *   R0^2 = (a1*sqrt(r4r2ij) + a2)^2 with a1=0.58, a2=4.80
     *
     * Key differences from standard BJ damping:
     * 1. R0 computed from sqrtZr4r2 product (NOT from C8/C6 ratio)
     * 2. C8 is implicit: factor 2*r4r2ij*t8 instead of separate C8*t8
     * 3. 0.5 factor for pair counting (each pair counted once)
     */

    if (CurcumaLogger::get_verbosity() >= 3 && m_d4_dispersions.size() > 0) {
        CurcumaLogger::info(fmt::format("Thread {} calculating {} GFN-FF/D4 dispersion pairs",
                                        m_thread, m_d4_dispersions.size()));
    }
    // Claude Generated (Feb 8, 2026): Per-pair diagnostic for dispersion accuracy investigation (verbosity >= 3 only)
    bool d4_disp_diag = (CurcumaLogger::get_verbosity() >= 3 && m_d4_dispersions.size() > 0 && m_d4_dispersions.size() <= 50);
    if (d4_disp_diag) {
        CurcumaLogger::info("D4_DISP_CSV: idx,i,j,rij_bohr,C6,r4r2ij,r0_sq,zetac6,t6,t8,energy");
    }

    for (int index = 0; index < m_d4_dispersions.size(); ++index) {
        const auto& disp = m_d4_dispersions[index];

        Eigen::Vector3d ri = geom().row(disp.i);
        Eigen::Vector3d rj = geom().row(disp.j);
        Eigen::Vector3d rij_vec = ri - rj;
        double rij = rij_vec.norm() * m_au;  // Convert to atomic units if needed

        if (rij > disp.r_cut || rij < 1e-10) continue;  // Skip if beyond cutoff or too close

        // GFN-FF modified BJ damping (NOT standard D3/D4!)
        // Reference: gfnff_gdisp0.f90:365-377
        //
        // Optimized power calculations: r^6 = (r^2)^3, r^8 = (r^2)^4
        double r2 = rij * rij;
        double r6 = r2 * r2 * r2;     // (r^2)^3 = r^6

        // R0^6 = (R0^2)^3 where r0_squared is pre-computed as (a1*sqrt(r4r2ij)+a2)^2
        double r0_6 = disp.r0_squared * disp.r0_squared * disp.r0_squared;

        // t6 = 1/(r^6 + R0^6)
        double t6 = 1.0 / (r6 + r0_6);

        // t8 = 1/(r^8 + R0^8)
        double r8 = r6 * r2;          // r^8 = r^6 * r^2
        double r0_8 = r0_6 * disp.r0_squared;  // R0^8 = R0^6 * R0^2
        double t8 = 1.0 / (r8 + r0_8);

        // GFN-FF dispersion formula: disp = (t6 + 2*r4r2ij*t8) * zetac6
        // Energy: dE = -c6 * disp (no 0.5 - see line 1586 comment)
        // Reference: gfnff_gdisp0.f90:374,377
        // Claude Generated (Jan 31, 2026): Added zetac6 charge scaling
        // Claude Generated (Feb 8, 2026): Removed 0.5 factor - consistent with energy fix
        double disp_sum = t6 + 2.0 * disp.r4r2ij * t8;
        double pair_energy = -disp.C6 * disp_sum * disp.zetac6 * m_final_factor;

        // Claude Generated (Jan 31, 2026): Fix dispersion energy reporting
        // Previously only m_d4_energy was updated, causing total_gfnff_dispersion to show 0
        // Now both are updated, matching CalculateGFNFFDispersionContribution() behavior
        m_dispersion_energy += pair_energy;
        m_d4_energy += pair_energy;

        // Debug: Per-pair diagnostic CSV (verbosity >= 3 only)
        if (d4_disp_diag) {
            CurcumaLogger::info(fmt::format("D4_DISP_CSV: {},{},{},{:.6f},{:.6e},{:.6f},{:.6f},{:.6e},{:.6e},{:.6e},{:.10e}",
                                           index, disp.i, disp.j,
                                           rij, disp.C6, disp.r4r2ij, disp.r0_squared, disp.zetac6,
                                           t6, t8, pair_energy));
        }

        if (m_calculate_gradient) {
            // Analytical gradient for GFN-FF dispersion formula
            // d/dr[E] = d/dr[-0.5 * C6 * (t6 + 2*r4r2ij*t8) * zetac6]
            //
            // Reference: gfnff_gdisp0.f90:371-372, 375
            // d6 = -6*r2**2*t6**2  -> derivative of t6 w.r.t. r
            // d8 = -8*r2**3*t8**2  -> derivative of t8 w.r.t. r
            // Note: zetac6 is constant during gradient evaluation (fixed charges)

            double d6 = -6.0 * r2 * r2 * t6 * t6;  // d(t6)/d(r^2) scaled appropriately
            double d8 = -8.0 * r2 * r2 * r2 * t8 * t8;  // d(t8)/d(r^2) scaled appropriately

            // ddisp/d(r^2) = d6 + 2*r4r2ij*d8
            double ddisp_dr2 = d6 + 2.0 * disp.r4r2ij * d8;

            // dE/dr = -0.5 * C6 * zetac6 * ddisp/d(r^2) * d(r^2)/dr = -0.5 * C6 * zetac6 * ddisp_dr2 * 2*r
            //       = -C6 * zetac6 * ddisp_dr2 * r
            // Claude Generated (Jan 31, 2026): Added zetac6 to gradient
            double dEdr = -disp.C6 * disp.zetac6 * ddisp_dr2 * rij * m_final_factor;

            Eigen::Vector3d grad = dEdr * rij_vec / rij;

            m_gradient.row(disp.i) += grad.transpose();
            m_gradient.row(disp.j) -= grad.transpose();

            // Claude Generated (Mar 2026): Dispersion dC6/dCN chain-rule contribution
            // Reference: Fortran gfnff_gdisp0.f90:382-395
            // CRITICAL FIX: Was missing in D4 path, causing dispersion gradient errors
            // (D3 path in CalculateGFNFFDispersionContribution already had this)
            if (m_dc6dcn_ptr && m_dc6dcn_ptr->size() > 0 &&
                disp.i < m_dc6dcn_ptr->rows() && disp.j < m_dc6dcn_ptr->cols()) {
                double disp_value = disp_sum * disp.zetac6 * m_final_factor;
                m_dEdcn(disp.i) -= (*m_dc6dcn_ptr)(disp.i, disp.j) * disp_value;
                m_dEdcn(disp.j) -= (*m_dc6dcn_ptr)(disp.j, disp.i) * disp_value;
            }
        }
    }

    if (CurcumaLogger::get_verbosity() >= 3 && m_d4_dispersions.size() > 0) {
        CurcumaLogger::param("thread_d4_energy", fmt::format("{:.6f} Eh", m_d4_energy));
    }

    if (CurcumaLogger::get_verbosity() >= 2 && m_d4_dispersions.size() > 0) {
        CurcumaLogger::result(fmt::format("GFN-FF D4 Dispersion Energy: {:.6e} Eh", m_d4_energy));
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
        Eigen::Vector3d pos_i = geom().row(triple.i).transpose();
        Eigen::Vector3d pos_j = geom().row(triple.j).transpose();
        Eigen::Vector3d pos_k = geom().row(triple.k).transpose();

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
        Eigen::Vector3d pos_i = geom().row(triple.i).transpose();
        Eigen::Vector3d pos_j = geom().row(triple.j).transpose();
        Eigen::Vector3d pos_k = geom().row(triple.k).transpose();

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

// BF (Bonded ATM/GFN-FF) - Claude Generated (January 17, 2026)
// Calculate bonded ATM (batm) energy and gradients for 1,4-pairs
// Reference: external/gfnff/src/gfnff_engrad.F90:3267-3334 (batmgfnff_eg subroutine)
// Formula: E_batm = c9 * (ang + 1.0) / rav3
// where c9 = ff * zb3atm_i * zb3atm_j * zb3atm_k and ff = (1 - 3*q_i)*(1 - 3*q_j)*(1 - 3*q_k)
void ForceFieldThread::CalculateGFNFFBatmContribution()
{
    m_batm_energy = 0.0;

    if (m_gfnff_batms.empty()) {
        return;
    }

    // Clamp to [-4, 4] to prevent runaway (matches Fortran gfnff_engrad.F90:3282-3287)
    const double fqq = 3.0;

    for (const auto& batm : m_gfnff_batms) {
        // Get atom positions and convert to Bohr
        Eigen::Vector3d i_pos = geom().row(batm.i).transpose() * m_au;
        Eigen::Vector3d j_pos = geom().row(batm.j).transpose() * m_au;
        Eigen::Vector3d k_pos = geom().row(batm.k).transpose() * m_au;

        // Calculate distance vectors (in Bohr)
        Eigen::Vector3d rij_vec = j_pos - i_pos;
        Eigen::Vector3d rik_vec = k_pos - i_pos;
        Eigen::Vector3d rjk_vec = k_pos - j_pos;

        // Squared distances
        double r2ij = rij_vec.squaredNorm();
        double r2jk = rjk_vec.squaredNorm();
        double r2ik = rik_vec.squaredNorm();

        // Distances
        double rij = std::sqrt(r2ij);
        double rjk = std::sqrt(r2jk);
        double rik = std::sqrt(r2ik);

        // Angular term: ang = 0.375 * (r_ik*r_jk - r_ij^2) * (r_ij*r_ik - r_jk^2) * (r_ij*r_jk - r_ik^2) / r_ijk^3
        // where r_ijk^3 = (r_ij*r_jk*r_ik)^2 = r2ij * r2jk * r2ik
        double rijk3 = r2ij * r2jk * r2ik;

        // mijk = -r_ij^2 + r_jk^2 + r_ik^2
        double mijk = -r2ij + r2jk + r2ik;
        // imjk = r_ij^2 - r_jk^2 + r_ik^2
        double imjk = r2ij - r2jk + r2ik;
        // ijmk = r_ij^2 + r_jk^2 - r_ik^2
        double ijmk = r2ij + r2jk - r2ik;

        // Angular term
        double ang = 0.375 * ijmk * imjk * mijk / rijk3;

        // rav3 = (r_ij*r_jk*r_ik)^1.5
        double rav3 = std::pow(rijk3, 1.5);  // This is (r_ij*r_jk*r_ik)^1.5 = rijk^1.5

        // Combined angular term
        double angr9 = (ang + 1.0) / rav3;

        // Charge factor: ff = (1 - 3*q_i) * (1 - 3*q_j) * (1 - 3*q_k)
        // Reference: Fortran gfnff_engrad.F90:620 passes topo%qa (Phase-1) to batmgfnff_eg
        // CRITICAL FIX (Feb 21, 2026): Use Phase-1 topology charges (fixed), NOT Phase-2 EEQ charges (dynamic)
        // Phase-2 charges change with geometry, but BATM gradient has no dq/dx term → energy drift in MD
        // Claude Generated (Mar 2026): Use pointer-based accessor for zero-copy access
        double fi = (1.0 - fqq * topo_q(batm.i));
        fi = std::min(std::max(fi, -4.0), 4.0);

        double fj = (1.0 - fqq * topo_q(batm.j));
        fj = std::min(std::max(fj, -4.0), 4.0);

        double fk = (1.0 - fqq * topo_q(batm.k));
        fk = std::min(std::max(fk, -4.0), 4.0);

        double ff_charge = fi * fj * fk;

        // Strength: c9 = ff * zb3atm_i * zb3atm_j * zb3atm_k
        // zb3atm parameters are already in unit Bohr^3
        double c9 = ff_charge * batm.zb3atm_i * batm.zb3atm_j * batm.zb3atm_k;

        // Energy (matches Fortran: energy = c9 * (ang + 1.0) / rijk^1.5)
        double energy = c9 * angr9;
        m_batm_energy += energy * m_final_factor;


        // ============================
        // Analytical Gradient Calculation
        // ============================

        if (m_calculate_gradient) {
            // Angular derivatives (from Fortran gfnff_engrad.F90:3308-3322)
            double dang_ij = -0.375 * (std::pow(r2ij, 3) + std::pow(r2ij, 2) * (r2jk + r2ik)
                                  + r2ij * (3.0 * std::pow(r2jk, 2) + 2.0 * r2jk * r2ik + 3.0 * std::pow(r2ik, 2))
                                  - 5.0 * std::pow((r2jk - r2ik), 2) * (r2jk + r2ik))
                                  / (rij * rijk3 * rav3);

            double dang_jk = -0.375 * (std::pow(r2jk, 3) + std::pow(r2jk, 2) * (r2ik + r2ij)
                                  + r2jk * (3.0 * std::pow(r2ik, 2) + 2.0 * r2ik * r2ij + 3.0 * std::pow(r2ij, 2))
                                  - 5.0 * std::pow((r2ik - r2ij), 2) * (r2ik + r2ij))
                                  / (rjk * rijk3 * rav3);

            double dang_ik = -0.375 * (std::pow(r2ik, 3) + std::pow(r2ik, 2) * (r2jk + r2ij)
                                  + r2ik * (3.0 * std::pow(r2jk, 2) + 2.0 * r2jk * r2ij + 3.0 * std::pow(r2ij, 2))
                                  - 5.0 * std::pow((r2jk - r2ij), 2) * (r2jk + r2ij))
                                  / (rik * rijk3 * rav3);

            // Gradient components
            Eigen::Vector3d dgij = -dang_ij * c9 * (rij_vec / rij);
            Eigen::Vector3d dgjk = -dang_jk * c9 * (rjk_vec / rjk);
            Eigen::Vector3d dgik = -dang_ik * c9 * (rik_vec / rik);

            // Accumulate gradients (from Fortran gfnff_engrad.F90:3324-3332)
            // Atom j: -dg_ij + dg_jk
            m_gradient.row(batm.j) += (-dgij + dgjk) * m_final_factor;
            // Atom k: -dg_ik - dg_jk
            m_gradient.row(batm.k) += (-dgik - dgjk) * m_final_factor;
            // Atom i: dg_ij + dg_ik
            m_gradient.row(batm.i) += (dgij + dgik) * m_final_factor;
        }
    }

}
