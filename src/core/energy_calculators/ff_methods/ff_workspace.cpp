/*
 * <FFWorkspace - Unified Force Field Workspace for Curcuma>
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
 * Claude Generated (March 2026): Core workspace logic — init, partition, calculate, reduce.
 */

#include "ff_workspace.h"
#include "src/core/curcuma_logger.h"
#include "src/core/units.h"

#include <fmt/core.h>
#include <fmt/format.h>

#include <algorithm>
#include <cmath>

FFWorkspace::FFWorkspace(int num_threads)
    : m_num_threads(std::max(1, num_threads))
{
}

void FFWorkspace::setInteractionLists(GFNFFParameterSet&& params)
{
    // Move all interaction lists from parameter set (zero-copy)
    m_bonds = std::move(params.bonds);
    m_angles = std::move(params.angles);
    m_dihedrals = std::move(params.dihedrals);
    m_extra_dihedrals = std::move(params.extra_dihedrals);
    m_inversions = std::move(params.inversions);
    m_storsions = std::move(params.storsions);

    // Dispersion routing: D4 vs D3 (same logic as ForceField::setGFNFFParameters)
    m_dispersions.clear();
    m_d4_dispersions.clear();
    if (params.dispersion_method == "d4") {
        m_d4_dispersions = std::move(params.dispersions);
    } else {
        m_dispersions = std::move(params.dispersions);
    }

    m_bonded_reps = std::move(params.bonded_repulsions);
    m_nonbonded_reps = std::move(params.nonbonded_repulsions);
    m_coulombs = std::move(params.coulombs);

    m_hbonds = std::move(params.hbonds);
    m_xbonds = std::move(params.xbonds);
    m_atm_triples = std::move(params.atm_triples);
    m_batm_triples = std::move(params.batm_triples);

    m_bond_hb_data = std::move(params.bond_hb_data);

    m_eeq_charges = std::move(params.eeq_charges);
    m_topology_charges = std::move(params.topology_charges);

    m_e0 = params.e0;

    m_dispersion_enabled = params.dispersion_enabled;
    m_hbond_enabled = params.hbond_enabled;
    m_repulsion_enabled = params.repulsion_enabled;
    m_coulomb_enabled = params.coulomb_enabled;

    // Build bonded pairs cache for fast repulsion lookup
    m_bonded_pairs.clear();
    for (const auto& bond : m_bonds) {
        m_bonded_pairs.insert({bond.i, bond.j});
        m_bonded_pairs.insert({bond.j, bond.i});
    }

    // Extract per-atom Coulomb self-energy parameters from pairs (for TERM 2+3 in postProcess)
    // Same logic as ForceField::setGFNFFParameters lines 418-444
    if (!m_coulombs.empty() && m_natoms > 0) {
        m_coul_chi_base = Vector::Zero(m_natoms);
        m_coul_gam = Vector::Zero(m_natoms);
        m_coul_alp = Vector::Zero(m_natoms);
        m_coul_cnf = Vector::Zero(m_natoms);
        m_coul_chi_static = Vector::Zero(m_natoms);

        std::vector<bool> atom_seen(m_natoms, false);
        for (const auto& coul : m_coulombs) {
            if (!atom_seen[coul.i]) {
                m_coul_chi_base(coul.i) = coul.chi_base_i;
                m_coul_gam(coul.i) = coul.gam_i;
                m_coul_alp(coul.i) = coul.alp_i;
                m_coul_cnf(coul.i) = coul.cnf_i;
                m_coul_chi_static(coul.i) = coul.chi_i;
                atom_seen[coul.i] = true;
            }
            if (!atom_seen[coul.j]) {
                m_coul_chi_base(coul.j) = coul.chi_base_j;
                m_coul_gam(coul.j) = coul.gam_j;
                m_coul_alp(coul.j) = coul.alp_j;
                m_coul_cnf(coul.j) = coul.cnf_j;
                m_coul_chi_static(coul.j) = coul.chi_j;
                atom_seen[coul.j] = true;
            }
        }
    }
}

void FFWorkspace::setAtomTypes(const std::vector<int>& atoms)
{
    m_atom_types = atoms;
    m_natoms = static_cast<int>(atoms.size());
}

void FFWorkspace::partition()
{
    int T = m_num_threads;
    m_partitions.resize(T);
    m_accumulators.resize(T);

    for (int t = 0; t < T; ++t) {
        auto& pr = m_partitions[t];
        pr.bonds = linearRange(m_bonds.size(), t, T);
        pr.angles = linearRange(m_angles.size(), t, T);
        pr.dihedrals = linearRange(m_dihedrals.size(), t, T);
        pr.extra_dihedrals = linearRange(m_extra_dihedrals.size(), t, T);
        pr.inversions = linearRange(m_inversions.size(), t, T);
        pr.storsions = linearRange(m_storsions.size(), t, T);
        pr.dispersions = linearRange(m_dispersions.size(), t, T);
        pr.d4_dispersions = linearRange(m_d4_dispersions.size(), t, T);
        pr.bonded_reps = linearRange(m_bonded_reps.size(), t, T);
        pr.nonbonded_reps = linearRange(m_nonbonded_reps.size(), t, T);
        pr.coulombs = linearRange(m_coulombs.size(), t, T);
        pr.hbonds = linearRange(m_hbonds.size(), t, T);
        pr.xbonds = linearRange(m_xbonds.size(), t, T);
        pr.atm_triples = linearRange(m_atm_triples.size(), t, T);
        pr.batm_triples = linearRange(m_batm_triples.size(), t, T);
    }

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info(fmt::format("FFWorkspace: partitioned for {} threads, {} atoms",
            T, m_natoms));
        CurcumaLogger::param("bonds", std::to_string(m_bonds.size()));
        CurcumaLogger::param("angles", std::to_string(m_angles.size()));
        CurcumaLogger::param("dihedrals", std::to_string(m_dihedrals.size()));
        CurcumaLogger::param("dispersions", std::to_string(m_dispersions.size()));
        CurcumaLogger::param("coulombs", std::to_string(m_coulombs.size()));
    }
}

void FFWorkspace::setCNDerivatives(const Vector& cn, const Vector& cnf,
                                    const std::vector<SpMatrix>& dcn)
{
    m_cn = cn;
    m_cnf = cnf;
    m_dcn = dcn;
}

void FFWorkspace::updateHBonds(const std::vector<GFNFFHydrogenBond>& hbonds)
{
    m_hbonds = hbonds;
    // Re-partition HB ranges only
    int T = m_num_threads;
    for (int t = 0; t < T; ++t) {
        m_partitions[t].hbonds = linearRange(m_hbonds.size(), t, T);
    }
}

void FFWorkspace::updateXBonds(const std::vector<GFNFFHalogenBond>& xbonds)
{
    m_xbonds = xbonds;
    // Re-partition XB ranges only
    int T = m_num_threads;
    for (int t = 0; t < T; ++t) {
        m_partitions[t].xbonds = linearRange(m_xbonds.size(), t, T);
    }
}

void FFWorkspace::setCoulombSelfEnergyParams(const Vector& chi_base, const Vector& gam,
                                               const Vector& alp, const Vector& cnf,
                                               const Vector& chi_static)
{
    m_coul_chi_base = chi_base;
    m_coul_gam = gam;
    m_coul_alp = alp;
    m_coul_cnf = cnf;
    m_coul_chi_static = chi_static;
}

double FFWorkspace::calculate(bool gradient)
{
    m_do_gradient = gradient;

    if (m_num_threads == 1) {
        // T=1: Direct call, zero pool overhead
        m_accumulators[0].reset(m_natoms, gradient, m_store_components);
        executeGFNFF(0);

        // acc[0] IS the result — zero-copy swap
        m_result_energy = m_accumulators[0].energy;
        if (gradient) {
            m_result_gradient.swap(m_accumulators[0].gradient);
            m_dEdcn_total = m_accumulators[0].dEdcn;
            m_dEdcn_bond_total = m_accumulators[0].dEdcn_bond;
        }
        if (m_store_components && gradient) {
            m_result_grad_bond.swap(m_accumulators[0].grad_bond);
            m_result_grad_angle.swap(m_accumulators[0].grad_angle);
            m_result_grad_torsion.swap(m_accumulators[0].grad_torsion);
            m_result_grad_repulsion.swap(m_accumulators[0].grad_repulsion);
            m_result_grad_coulomb.swap(m_accumulators[0].grad_coulomb);
            m_result_grad_dispersion.swap(m_accumulators[0].grad_dispersion);
            m_result_grad_hb.swap(m_accumulators[0].grad_hb);
            m_result_grad_xb.swap(m_accumulators[0].grad_xb);
            m_result_grad_batm.swap(m_accumulators[0].grad_batm);
            m_result_grad_atm.swap(m_accumulators[0].grad_atm);
        }
    } else {
        // T>1: Pool + barrier
        for (int t = 0; t < m_num_threads; ++t)
            m_accumulators[t].reset(m_natoms, gradient, m_store_components);

        std::vector<std::future<void>> futures;
        futures.reserve(m_num_threads - 1);
        for (int t = 1; t < m_num_threads; ++t)
            futures.push_back(m_pool->enqueue([this, t]() { executeGFNFF(t); }));
        executeGFNFF(0);  // Main thread works on partition 0
        for (auto& f : futures)
            f.get();

        reduce();
    }

    postProcess(gradient);
    return m_e0 + m_result_energy.total();
}

void FFWorkspace::executeGFNFF(int p)
{
    // HB coordination numbers must be computed before bond energy
    computeHBCoordinationNumbers(p);

    // Bonded terms
    calcBonds(p);
    calcAngles(p);
    calcDihedrals(p);
    calcExtraTorsions(p);
    calcInversions(p);

    // Non-bonded pairwise terms
    if (m_dispersion_enabled) {
        calcDispersion(p);
        calcD4Dispersion(p);
    }
    if (m_repulsion_enabled) {
        calcBondedRepulsion(p);
        calcNonbondedRepulsion(p);
    }
    if (m_coulomb_enabled) {
        calcCoulomb(p);
    }

    // HB/XB three-body terms
    if (m_hbond_enabled) {
        calcHydrogenBonds(p);
        calcHalogenBonds(p);
    }

    // Three-body dispersion
    if (!m_atm_triples.empty()) {
        calcATM(p);
        if (m_do_gradient)
            calcATMGradient(p);
    }

    // Bonded ATM
    if (!m_batm_triples.empty()) {
        calcBATM(p);
    }

    // Triple bond torsions
    if (!m_storsions.empty()) {
        calcSTorsions(p);
    }
}

void FFWorkspace::reduce()
{
    // Sum all accumulators into result
    m_result_energy = m_accumulators[0].energy;
    if (m_do_gradient) {
        m_result_gradient = m_accumulators[0].gradient;
        m_dEdcn_total = m_accumulators[0].dEdcn;
        m_dEdcn_bond_total = m_accumulators[0].dEdcn_bond;
    }
    if (m_store_components && m_do_gradient) {
        m_result_grad_bond = m_accumulators[0].grad_bond;
        m_result_grad_angle = m_accumulators[0].grad_angle;
        m_result_grad_torsion = m_accumulators[0].grad_torsion;
        m_result_grad_repulsion = m_accumulators[0].grad_repulsion;
        m_result_grad_coulomb = m_accumulators[0].grad_coulomb;
        m_result_grad_dispersion = m_accumulators[0].grad_dispersion;
        m_result_grad_hb = m_accumulators[0].grad_hb;
        m_result_grad_xb = m_accumulators[0].grad_xb;
        m_result_grad_batm = m_accumulators[0].grad_batm;
        m_result_grad_atm = m_accumulators[0].grad_atm;
    }

    for (int t = 1; t < m_num_threads; ++t) {
        m_result_energy += m_accumulators[t].energy;
        if (m_do_gradient) {
            m_result_gradient += m_accumulators[t].gradient;
            m_dEdcn_total += m_accumulators[t].dEdcn;
            m_dEdcn_bond_total += m_accumulators[t].dEdcn_bond;
        }
        if (m_store_components && m_do_gradient) {
            m_result_grad_bond += m_accumulators[t].grad_bond;
            m_result_grad_angle += m_accumulators[t].grad_angle;
            m_result_grad_torsion += m_accumulators[t].grad_torsion;
            m_result_grad_repulsion += m_accumulators[t].grad_repulsion;
            m_result_grad_coulomb += m_accumulators[t].grad_coulomb;
            m_result_grad_dispersion += m_accumulators[t].grad_dispersion;
            m_result_grad_hb += m_accumulators[t].grad_hb;
            m_result_grad_xb += m_accumulators[t].grad_xb;
            m_result_grad_batm += m_accumulators[t].grad_batm;
            m_result_grad_atm += m_accumulators[t].grad_atm;
        }
    }
}

void FFWorkspace::postProcess(bool gradient)
{
    // =========================================================================
    // Coulomb TERM 2+3: Self-energy (sequential, thread-count-independent)
    // Reference: Fortran gfnff_engrad.F90:1678-1679
    // =========================================================================
    if (m_coul_gam.size() == m_natoms && m_eeq_charges.size() == m_natoms) {
        const double sqrt_2_over_pi = 0.797884560802865;
        const bool has_cn = (m_cn.size() == m_natoms);
        double E_en = 0.0, E_self = 0.0;

        for (int i = 0; i < m_natoms; ++i) {
            if (m_coul_alp(i) <= 0.0) continue;
            double q = m_eeq_charges(i);
            if (std::isnan(q)) continue;
            double chi;
            if (m_coul_cnf(i) != 0.0 && has_cn) {
                chi = m_coul_chi_base(i) + m_coul_cnf(i) * std::sqrt(std::max(m_cn(i), 0.0));
            } else {
                chi = m_coul_chi_static(i);
            }
            E_en -= q * chi;
            E_self += 0.5 * q * q * (m_coul_gam(i) + sqrt_2_over_pi / std::sqrt(m_coul_alp(i)));
        }
        m_result_energy.coulomb += E_en + E_self;

        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info(fmt::format("  Coulomb self-energy (workspace): EN={:+.12f}, Self={:+.12f} Eh", E_en, E_self));
        }
    }

    // =========================================================================
    // dEdcn chain-rule gradient + Coulomb TERM 1b
    // Reference: Fortran gfnff_engrad.F90:418-422 (bond/disp), 449-454 (coulomb)
    // =========================================================================
    if (gradient && !m_dcn.empty() && m_dcn.size() == 3) {
        // Compute TERM 1b qtmp
        Vector qtmp = Vector::Zero(m_natoms);
        bool has_term1b = (m_eeq_charges.size() == m_natoms &&
                          m_cnf.size() == m_natoms &&
                          m_cn.size() == m_natoms);
        if (has_term1b) {
            for (int i = 0; i < m_natoms; ++i) {
                double cn_i = std::max(m_cn(i), 0.0);
                qtmp(i) = m_eeq_charges(i) * m_cnf(i) / (2.0 * std::sqrt(cn_i) + 1e-16);
            }
        }

        // Combined matvec: gradient += dcn * (dEdcn_total - qtmp)
        Vector dEdcn_combined = has_term1b ? (m_dEdcn_total - qtmp).eval() : m_dEdcn_total;
        for (int dim = 0; dim < 3; ++dim) {
            if (m_dcn[dim].rows() == m_natoms && m_dcn[dim].cols() == m_natoms) {
                m_result_gradient.col(dim) += m_dcn[dim] * dEdcn_combined;
            }
        }

        // Per-component CN corrections
        if (m_store_components) {
            Vector dEdcn_disp = m_dEdcn_total - m_dEdcn_bond_total;
            for (int dim = 0; dim < 3; ++dim) {
                if (m_dcn[dim].rows() == m_natoms && m_dcn[dim].cols() == m_natoms) {
                    m_result_grad_bond.col(dim) += m_dcn[dim] * m_dEdcn_bond_total;
                    m_result_grad_dispersion.col(dim) += m_dcn[dim] * dEdcn_disp;
                    if (has_term1b)
                        m_result_grad_coulomb.col(dim) -= m_dcn[dim] * qtmp;
                }
            }
        }
    }
}
