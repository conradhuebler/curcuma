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
    if (m_record_hb_energies)
        m_hb_triple_energy.assign(m_hbonds.size(), 0.0);
    m_xbonds = std::move(params.xbonds);
    m_atm_triples = std::move(params.atm_triples);
    m_batm_triples = std::move(params.batm_triples);

    m_bond_hb_data = std::move(params.bond_hb_data);

    // UFF/QMDFF non-bonded pairs
    m_vdws = std::move(params.vdws);

    m_eeq_charges = std::move(params.eeq_charges);
    m_topology_charges = std::move(params.topology_charges);

    m_e0 = params.e0;

    // Method type and distance unit factor
    m_method_type = params.method_type;
    m_au = (m_method_type != FFMethodType::GFN_FF) ? CurcumaUnit::Length::angstrom_to_bohr_or_legacy(1.889726125) : 1.0;

    m_dispersion_enabled = params.dispersion_enabled;
    m_hbond_enabled = params.hbond_enabled;
    m_repulsion_enabled = params.repulsion_enabled;
    m_coulomb_enabled = params.coulomb_enabled;
    m_coulomb_implicit = params.coulomb_implicit;
    m_coulomb_implicit_rcut = params.coulomb_implicit_rcut;

    // Build bonded pairs cache for fast repulsion lookup
    m_bonded_pairs.clear();
    for (const auto& bond : m_bonds) {
        m_bonded_pairs.insert({bond.i, bond.j});
        m_bonded_pairs.insert({bond.j, bond.i});
    }

    // Per-atom Coulomb self-energy parameters (for TERM 2+3 in postProcess).
    // Claude Generated (Sep 2026): prefer the dedicated per-atom fields
    // (GFNFF::generateCoulombSelfEnergyNative()), which are populated regardless
    // of pair count. Falling back to scanning m_coulombs (the previous, buggy
    // behaviour) only when a producer of GFNFFParameterSet doesn't fill the new
    // fields (e.g. UFF/QMDFF, which share this struct but have no EEQ self-energy)
    // — that fallback is empty for an isolated atom (0 pairs for N=1), which was
    // exactly the bug: native GFN-FF returned 0.0 Eh for any single charged atom
    // because the pair-derived vectors stayed empty and postProcess()'s
    // `m_coul_gam.size() == m_natoms` guard silently skipped the self-energy.
    if (params.coul_self_gam.size() == m_natoms && m_natoms > 0) {
        m_coul_chi_base = std::move(params.coul_self_chi_base);
        m_coul_gam = std::move(params.coul_self_gam);
        m_coul_alp = std::move(params.coul_self_alp);
        m_coul_cnf = std::move(params.coul_self_cnf);
        m_coul_chi_static = std::move(params.coul_self_chi_static);
    } else if (!m_coulombs.empty() && m_natoms > 0) {
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
        // Implicit Coulomb: split the OUTER atom index so that every thread gets roughly the
        // same number of (i<j) pairs - atom i carries natoms-1-i of them, so equal atom counts
        // would leave thread 0 with most of the work. Claude Generated (Sep 2026).
        if (m_coulomb_implicit && m_natoms > 1) {
            const double total = 0.5 * static_cast<double>(m_natoms) * (m_natoms - 1);
            auto atom_at = [&](int k) {                 // first atom whose prefix >= k/T of total
                if (k <= 0) return 0;
                if (k >= T) return m_natoms;
                const double want = total * k / T;
                // pairs(i) = i*natoms - i*(i+1)/2 solved for i
                const double n = m_natoms;
                const double disc = (n - 0.5) * (n - 0.5) - 2.0 * want;
                const int i = static_cast<int>(std::ceil((n - 0.5) - std::sqrt(std::max(0.0, disc))));
                return std::min(m_natoms, std::max(0, i));
            };
            pr.coulomb_atoms = { atom_at(t), atom_at(t + 1) };
        } else {
            pr.coulomb_atoms = { 0, 0 };
        }
        pr.hbonds = linearRange(m_hbonds.size(), t, T);
        pr.xbonds = linearRange(m_xbonds.size(), t, T);
        pr.atm_triples = linearRange(m_atm_triples.size(), t, T);
        pr.batm_triples = linearRange(m_batm_triples.size(), t, T);
        pr.vdws = linearRange(m_vdws.size(), t, T);
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
                                    const CNDerivStore& dcn)
{
    m_cn = cn;
    m_cnf = cnf;
    m_dcn = dcn;
}

// Claude Generated (Aug 2026): react topology mode — full list swap on a live workspace.
void FFWorkspace::rebuildInteractionLists(GFNFFParameterSet&& params)
{
    setInteractionLists(std::move(params));
    partition();
}

void FFWorkspace::setRecordHBondEnergies(bool on)
{
    m_record_hb_energies = on;
    m_hb_triple_energy.clear();
    if (on)
        m_hb_triple_energy.assign(m_hbonds.size(), 0.0);
}

std::vector<GFNFFHBondEnergy> FFWorkspace::hydrogenBondEnergies() const
{
    std::vector<GFNFFHBondEnergy> out;
    if (!m_record_hb_energies || m_hb_triple_energy.size() != m_hbonds.size())
        return out;
    out.reserve(m_hbonds.size());
    for (std::size_t i = 0; i < m_hbonds.size(); ++i)
        out.push_back({ m_hbonds[i].i, m_hbonds[i].j, m_hbonds[i].k, m_hbonds[i].case_type, m_hb_triple_energy[i] });
    return out;
}

void FFWorkspace::updateHBonds(const std::vector<GFNFFHydrogenBond>& hbonds)
{
    m_hbonds = hbonds;
    if (m_record_hb_energies)
        m_hb_triple_energy.assign(m_hbonds.size(), 0.0);
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

void FFWorkspace::updateRepulsion(const std::vector<GFNFFRepulsion>& bonded_reps,
                                   const std::vector<GFNFFRepulsion>& nonbonded_reps)
{
    m_bonded_reps = bonded_reps;
    m_nonbonded_reps = nonbonded_reps;
    // Re-partition repulsion ranges only (same pattern as updateHBonds/updateXBonds)
    int T = m_num_threads;
    for (int t = 0; t < T; ++t) {
        m_partitions[t].bonded_reps = linearRange(m_bonded_reps.size(), t, T);
        m_partitions[t].nonbonded_reps = linearRange(m_nonbonded_reps.size(), t, T);
    }
}

void FFWorkspace::updateD4Dispersions(std::vector<GFNFFDispersion>&& pairs)
{
    // Claude Generated (Sep 2026): same swap + re-partition pattern as updateRepulsion().
    m_d4_dispersions = std::move(pairs);
    int T = m_num_threads;
    for (int t = 0; t < T; ++t)
        m_partitions[t].d4_dispersions = linearRange(m_d4_dispersions.size(), t, T);
}

void FFWorkspace::updateCoulombPairs(std::vector<GFNFFCoulomb>&& pairs)
{
    // Claude Generated (Sep 2026): same swap + re-partition pattern as updateRepulsion().
    // The per-atom self-energy parameters (m_coul_*) are topology-only and stay untouched.
    m_coulombs = std::move(pairs);
    int T = m_num_threads;
    for (int t = 0; t < T; ++t)
        m_partitions[t].coulombs = linearRange(m_coulombs.size(), t, T);
}

void FFWorkspace::updateBondHBData(const std::vector<int>& nr_hb, std::vector<BondHBEntry>&& data)
{
    // Claude Generated (Sep 2026): see the declaration. nr_hb is indexed like m_bonds.
    const size_t nb = std::min(nr_hb.size(), m_bonds.size());
    for (size_t b = 0; b < nb; ++b)
        m_bonds[b].nr_hb = nr_hb[b];
    m_bond_hb_data = std::move(data);
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

// Claude Generated (Nov 2026, WP-PERIODIC-NONBONDED B2)
void FFWorkspace::setEwald(const curcuma::ewald::EwaldParameters& p)
{
    if (!m_cell_bohr.active)
        throw std::runtime_error("FFWorkspace::setEwald: no periodic cell set");
    m_ewald = std::make_unique<curcuma::ewald::EwaldSum>(m_cell_bohr, p);  // throws for r_cut > L_min/2
    m_coulomb_ewald = true;
}

// Claude Generated (Nov 2026, S2): the Ewald pair part alone, for the GPU hybrid (see ff_workspace.h)
double FFWorkspace::evaluateEwaldCoulomb(const Matrix& geometry_bohr, const Vector& charges, bool gradient,
                                         Matrix* gradient_out, curcuma::ewald::EwaldEnergyTerms* terms_out)
{
    if (!ewaldActive() || m_method_type != FFMethodType::GFN_FF)
        throw std::runtime_error("FFWorkspace::evaluateEwaldCoulomb: Ewald sum not set up (setEwald)");
    if (geometry_bohr.rows() != m_natoms || m_coul_alp.size() != m_natoms || charges.size() != m_natoms)
        throw std::runtime_error("FFWorkspace::evaluateEwaldCoulomb: geometry, alp and charges must have natoms rows");
    m_geometry = geometry_bohr;
    m_eeq_charges = charges;
    m_do_gradient = gradient;
    m_ewald_terms = curcuma::ewald::EwaldEnergyTerms();
    m_ewald_setup_ms = 0.0;
    {
        // the same build (or lookup, when the periodic EEQ solve of this step built it) as in calculate()
        auto t0 = std::chrono::high_resolution_clock::now();
        if (!m_ewald_cache)
            m_ewald_cache = std::make_shared<curcuma::ewald::EwaldGeometryCache>();
        const Eigen::MatrixXd pos = m_geometry;
        auto run = [this](int T, const std::function<void(int)>& f) {
            std::vector<std::future<void>> futures;
            futures.reserve(T - 1);
            for (int t = 1; t < T; ++t)
                futures.push_back(m_pool->enqueue([&f, t]() { f(t); }));
            f(0);
            for (auto& fu : futures)
                fu.get();
        };
        const int T = (m_pool && m_num_threads > 1) ? m_num_threads : 1;
        m_ewald->buildGeometryCache(pos, &m_coul_alp, *m_ewald_cache, m_ewald_kernel_cap_bytes, T,
            T > 1 ? curcuma::ewald::EwaldParallelFor(run) : curcuma::ewald::EwaldParallelFor());
        m_ewald_setup_ms = std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - t0).count();
    }
    for (int t = 0; t < m_num_threads; ++t)
        m_accumulators[t].reset(m_natoms, gradient, /*store_components=*/false, /*want_potentials=*/false);
    m_virial_valid = false;  // N1: the hybrid result holds no atomic virial of the local terms
    if (m_num_threads == 1 || !m_pool) {
        calcCoulombEwald(0);
    } else {
        std::vector<std::future<void>> futures;
        futures.reserve(m_num_threads - 1);
        for (int t = 1; t < m_num_threads; ++t)
            futures.push_back(m_pool->enqueue([this, t]() { calcCoulombEwald(t); }));
        calcCoulombEwald(0);
        for (auto& f : futures)
            f.get();
    }
    // partition order, as reduce()
    double e_pair = 0.0;
    m_ewald_real_ms = m_ewald_recip_ms = 0.0;
    if (gradient && gradient_out) {
        if (gradient_out->rows() != m_natoms || gradient_out->cols() != 3)
            gradient_out->resize(m_natoms, 3);
        gradient_out->setZero();
    }
    for (int t = 0; t < m_num_threads; ++t) {
        const auto& acc = m_accumulators[t];
        e_pair += acc.energy.coulomb;
        m_ewald_terms.real += acc.ewald_real;
        m_ewald_terms.reciprocal += acc.ewald_recip;
        m_ewald_real_ms += acc.ewald_real_ms;
        m_ewald_recip_ms += acc.ewald_recip_ms;
        if (gradient && gradient_out)
            *gradient_out += acc.gradient;
    }
    curcuma::ewald::EwaldResult tail;
    m_ewald->selfAndBackground(m_eeq_charges, &m_coul_alp, /*want_potentials=*/false, tail);
    m_ewald_terms.self = tail.terms.self;
    m_ewald_terms.background = tail.terms.background;
    m_ewald_terms.gaussian_background = tail.terms.gaussian_background;
    e_pair += tail.terms.self + tail.terms.background + tail.terms.gaussian_background;
    if (terms_out)
        *terms_out = m_ewald_terms;
    return e_pair;
}

double FFWorkspace::calculate(bool gradient)
{
    const bool do_timing = (CurcumaLogger::get_verbosity() >= 2);
    auto t_calc_start = do_timing ? std::chrono::high_resolution_clock::now() : std::chrono::time_point<std::chrono::high_resolution_clock>{};
    double t_execute = 0.0, t_reduce = 0.0, t_post = 0.0;

    m_do_gradient = gradient;

    // Periodic Coulomb (Claude Generated, Nov 2026, B2): the real-space pair list (minimum image, r < r_cut,
    // linked cells) and the phase table of this geometry, built once here before the partitions read them.
    // Not reached without setEwald(): the non-periodic path is bitwise what it was.
    m_ewald_terms = curcuma::ewald::EwaldEnergyTerms();
    m_ewald_setup_ms = 0.0;
    const bool ewald_now = ewaldActive() && m_coulomb_enabled && m_method_type == FFMethodType::GFN_FF;
    if (ewald_now) {
        if (m_coul_alp.size() != m_natoms || m_eeq_charges.size() != m_natoms)
            throw std::runtime_error("FFWorkspace: Ewald Coulomb needs per-atom alp and EEQ charges of size natoms");
        auto t0 = std::chrono::high_resolution_clock::now();
        // P1 (Nov 2026): built once per geometry; when the periodic EEQ solve of this step already built it for the
        // same positions, widths and parameters (shared object, bitwise key), this is a lookup.
        if (!m_ewald_cache)
            m_ewald_cache = std::make_shared<curcuma::ewald::EwaldGeometryCache>();
        const Eigen::MatrixXd pos = m_geometry;  // column-major copy (RowMajor workspace geometry)
        auto run = [this](int T, const std::function<void(int)>& f) {
            std::vector<std::future<void>> futures;
            futures.reserve(T - 1);
            for (int t = 1; t < T; ++t)
                futures.push_back(m_pool->enqueue([&f, t]() { f(t); }));
            f(0);
            for (auto& fu : futures)
                fu.get();
        };
        const int T = (m_pool && m_num_threads > 1) ? m_num_threads : 1;
        m_ewald->buildGeometryCache(pos, &m_coul_alp, *m_ewald_cache, m_ewald_kernel_cap_bytes, T,
            T > 1 ? curcuma::ewald::EwaldParallelFor(run) : curcuma::ewald::EwaldParallelFor());
        m_ewald_setup_ms = std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - t0).count();
    }
    const bool want_pot = ewald_now && m_want_coulomb_potentials;
    // N1 (Nov 2026): the atomic virial is a by-product of the gradient kernels (virial instantiations)
    const bool want_virial = m_want_virial && gradient && m_method_type == FFMethodType::GFN_FF;
    m_virial_valid = false;
    if (want_virial)
        m_result_virial.reset();

    // Select execute function based on method type
    auto executeMethod = [this](int t) {
        if (m_method_type == FFMethodType::UFF)
            executeUFF(t);
        else if (m_method_type == FFMethodType::QMDFF)
            executeQMDFF(t);
        else if (m_method_type == FFMethodType::CG)
            executeCG(t);
        else
            executeGFNFF(t);
    };

    // Thread-safety fix (Jun 2026): GFN-FF HB coordination numbers write the SHARED
    // m_bonds[].hb_cn_H and m_hb_grad_entries, which calcBonds() then reads. Previously
    // computeHBCoordinationNumbers(0) ran inside partition 0's executeGFNFF with NO barrier
    // before partitions 1..N ran calcBonds — so the other threads read a half-written
    // hb_cn_H (wrong bond exponent -> non-deterministic bond energy, ~0.2 Eh drift on
    // many-fragment systems) and raced on the m_hb_grad_entries push_back. Compute it ONCE
    // here on the main thread, before the parallel dispatch; partitions now only read it.
    if (m_method_type == FFMethodType::GFN_FF)
        computeHBCoordinationNumbers(0);

    auto t0 = do_timing ? std::chrono::high_resolution_clock::now() : std::chrono::time_point<std::chrono::high_resolution_clock>{};
    if (m_num_threads == 1) {
        // T=1: Direct call, zero pool overhead
        m_accumulators[0].reset(m_natoms, gradient, m_store_components, want_pot, want_virial);
        executeMethod(0);

        // acc[0] IS the result — zero-copy swap
        m_result_energy = m_accumulators[0].energy;
        m_result_timings = m_accumulators[0].timings;  // Claude Generated (May 2026): per-term timing
        m_ewald_terms.real = m_accumulators[0].ewald_real;       // B2
        m_ewald_terms.reciprocal = m_accumulators[0].ewald_recip;
        m_ewald_real_ms = m_accumulators[0].ewald_real_ms;
        m_ewald_recip_ms = m_accumulators[0].ewald_recip_ms;
        if (want_virial)
            m_result_virial = m_accumulators[0].virial;  // N1
        if (want_pot)
            m_result_coulomb_potentials.swap(m_accumulators[0].coulomb_potentials);
        else
            m_result_coulomb_potentials.resize(0);
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
            m_accumulators[t].reset(m_natoms, gradient, m_store_components, want_pot, want_virial);

        std::vector<std::future<void>> futures;
        futures.reserve(m_num_threads - 1);
        for (int t = 1; t < m_num_threads; ++t)
            futures.push_back(m_pool->enqueue([this, t, &executeMethod]() { executeMethod(t); }));
        executeMethod(0);  // Main thread works on partition 0
        for (auto& f : futures)
            f.get();

        reduce();
        if (do_timing) {
            t_reduce = std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - t0).count() - t_execute;
        }
    }
    if (do_timing) {
        t_execute = std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - t0).count();
    }

    t0 = do_timing ? std::chrono::high_resolution_clock::now() : std::chrono::time_point<std::chrono::high_resolution_clock>{};
    postProcess(gradient);
    if (do_timing) {
        t_post = std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - t0).count();
    }
    m_virial_valid = want_virial;  // N1: components summed in reduce()/the T=1 copy, remainder in postProcess()

    // CPU ENERGY TERMS (verbosity >= 3)
    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("=== CPU ENERGY TERMS ===");
        fmt::print("  bond      = {:+.15e}\n", m_result_energy.bond);
        fmt::print("  angle     = {:+.15e}\n", m_result_energy.angle);
        fmt::print("  dihedral  = {:+.15e}\n", m_result_energy.dihedral);
        fmt::print("  inversion = {:+.15e}\n", m_result_energy.inversion);
        fmt::print("  stors     = {:+.15e}\n", m_result_energy.stors);
        fmt::print("  batm      = {:+.15e}\n", m_result_energy.batm);
        fmt::print("  atm       = {:+.15e}\n", m_result_energy.atm);
        fmt::print("  disp      = {:+.15e}\n", m_result_energy.dispersion);
        fmt::print("  brep      = {:+.15e}\n", m_result_energy.bonded_rep);
        fmt::print("  nbrep     = {:+.15e}\n", m_result_energy.nonbonded_rep);
        fmt::print("  coulomb   = {:+.15e}\n", m_result_energy.coulomb);
        fmt::print("  hbond     = {:+.15e}\n", m_result_energy.hbond);
        fmt::print("  xbond     = {:+.15e}\n", m_result_energy.xbond);
        CurcumaLogger::info("=== CPU ENERGY END ===");
    }

    return m_e0 + m_result_energy.total();
}

void FFWorkspace::executeGFNFF(int p)
{
    // Claude Generated (May 2026): Per-term timing — fills m_accumulators[p].timings.
    // Aggregated across partitions in reduce() for the GFNFFEnergyReport CPU-sum column.
    auto& timings = m_accumulators[p].timings;
    auto tic = []() { return std::chrono::high_resolution_clock::now(); };
    auto toc = [](auto t0) {
        return std::chrono::duration<double, std::milli>(
            std::chrono::high_resolution_clock::now() - t0).count();
    };

    // HB coordination numbers are computed once on the main thread in calculate()
    // BEFORE the parallel dispatch (they write shared m_bonds[].hb_cn_H / m_hb_grad_entries
    // that every partition's calcBonds reads) — see the thread-safety note there.

    // Bonded terms
    auto t = tic(); calcBonds(p);          timings.bonds      = toc(t);
    t = tic();      calcAngles(p);         timings.angles     = toc(t);
    t = tic();      calcDihedrals(p);
                    calcExtraTorsions(p);  timings.dihedrals  = toc(t);
    t = tic();      calcInversions(p);     timings.inversions = toc(t);

    // Non-bonded pairwise terms
    if (m_dispersion_enabled) {
        t = tic();
        calcDispersion(p);
        calcD4Dispersion(p);
        timings.dispersion = toc(t);
    }
    if (m_repulsion_enabled) {
        t = tic(); calcBondedRepulsion(p);    timings.bonded_rep    = toc(t);
        t = tic(); calcNonbondedRepulsion(p); timings.nonbonded_rep = toc(t);
    }
    if (m_coulomb_enabled) {
        t = tic();
        if (ewaldActive())
            calcCoulombEwald(p);  // Nov 2026 (B2): periodic cell with coulomb_ewald
        else
            calcCoulomb(p);
        timings.coulomb = toc(t);
    }

    // HB/XB three-body terms
    if (m_hbond_enabled) {
        t = tic(); calcHydrogenBonds(p); timings.hbond = toc(t);
        t = tic(); calcHalogenBonds(p);  timings.xbond = toc(t);
    }

    // Three-body dispersion
    if (!m_atm_triples.empty()) {
        t = tic();
        calcATM(p);
        if (m_do_gradient) calcATMGradient(p);
        timings.atm = toc(t);
    }

    // Bonded ATM
    if (!m_batm_triples.empty()) {
        t = tic(); calcBATM(p); timings.batm = toc(t);
    }

    // Triple bond torsions
    if (!m_storsions.empty()) {
        t = tic(); calcSTorsions(p); timings.stors = toc(t);
    }
}

void FFWorkspace::reduce()
{
    // Sum all accumulators into result
    m_result_energy = m_accumulators[0].energy;
    m_result_timings = m_accumulators[0].timings;  // Claude Generated (May 2026)
    m_ewald_terms.real = m_accumulators[0].ewald_real;  // B2 (Nov 2026): fixed partition order, reproducible per thread count
    m_ewald_terms.reciprocal = m_accumulators[0].ewald_recip;
    m_ewald_real_ms = m_accumulators[0].ewald_real_ms;
    m_ewald_recip_ms = m_accumulators[0].ewald_recip_ms;
    m_result_coulomb_potentials = m_accumulators[0].coulomb_potentials;
    const bool virial = m_want_virial && m_do_gradient && m_method_type == FFMethodType::GFN_FF;  // N1
    if (virial)
        m_result_virial = m_accumulators[0].virial;
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
        m_result_timings += m_accumulators[t].timings;  // Claude Generated (May 2026): CPU-sum
        m_ewald_terms.real += m_accumulators[t].ewald_real;
        m_ewald_terms.reciprocal += m_accumulators[t].ewald_recip;
        m_ewald_real_ms += m_accumulators[t].ewald_real_ms;
        m_ewald_recip_ms += m_accumulators[t].ewald_recip_ms;
        if (m_result_coulomb_potentials.size() == m_natoms && m_accumulators[t].coulomb_potentials.size() == m_natoms)
            m_result_coulomb_potentials += m_accumulators[t].coulomb_potentials;
        if (virial)
            m_result_virial += m_accumulators[t].virial;  // N1, partition order
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
    const bool do_timing = (CurcumaLogger::get_verbosity() >= 2);
    auto t_post_start = do_timing ? std::chrono::high_resolution_clock::now() : std::chrono::time_point<std::chrono::high_resolution_clock>{};
    double t_self_energy = 0.0, t_chainrule = 0.0;

    // =========================================================================
    // Coulomb TERM 2+3: Self-energy (sequential, thread-count-independent)
    // Reference: Fortran gfnff_engrad.F90:1678-1679
    // =========================================================================
    auto t0 = do_timing ? std::chrono::high_resolution_clock::now() : std::chrono::time_point<std::chrono::high_resolution_clock>{};
    // Nov 2026 (S1): `coulomb false` switches the whole Coulomb term off, the self term and TERM 1b included, as the
    // GPU workspace does (before, the CPU kept +0.71 Eh of self energy on water8 and the two sides disagreed).
    if (m_coulomb_enabled && m_coul_gam.size() == m_natoms && m_eeq_charges.size() == m_natoms) {
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
    // Periodic Coulomb (Claude Generated, Nov 2026, B2): Ewald self term, neutralising background and the
    // Gaussian background (pi Q/V) sum q_i alp_i of the point-charge-plus-correction sum (ewald.h). Sequential,
    // O(N), thread-count independent. The GFN-FF self term above (Gaussian self-energy, chi) is untouched.
    if (ewaldActive() && m_coulomb_enabled && m_method_type == FFMethodType::GFN_FF
        && m_coul_alp.size() == m_natoms && m_eeq_charges.size() == m_natoms) {
        curcuma::ewald::EwaldResult tail;
        const bool want_pot = m_result_coulomb_potentials.size() == m_natoms;
        if (want_pot)
            tail.potentials = Eigen::VectorXd::Zero(m_natoms);
        m_ewald->selfAndBackground(m_eeq_charges, &m_coul_alp, want_pot, tail);
        m_ewald_terms.self = tail.terms.self;
        m_ewald_terms.background = tail.terms.background;
        m_ewald_terms.gaussian_background = tail.terms.gaussian_background;
        m_result_energy.coulomb += tail.terms.self + tail.terms.background + tail.terms.gaussian_background;
        if (want_pot)
            m_result_coulomb_potentials += tail.potentials;
        if (CurcumaLogger::get_verbosity() >= 3)
            CurcumaLogger::info(fmt::format("  Coulomb Ewald (workspace): real={:+.12f} recip={:+.12f} self={:+.12f} bg={:+.12f} gauss_bg={:+.12f} Eh",
                m_ewald_terms.real, m_ewald_terms.reciprocal, m_ewald_terms.self, m_ewald_terms.background,
                m_ewald_terms.gaussian_background));
    }
    if (do_timing) {
        t_self_energy = std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - t0).count();
    }

    // Homogeneous long-range D4 tail (Claude Generated, Nov 2026, A3c): energy computed by GFNFF from the
    // current CN weights, see setDispersionTail(); its dE/dCN joins the chain rule below.
    if (m_disp_tail_active)
        m_result_energy.disp_tail = m_disp_tail_energy;

    // =========================================================================
    // dEdcn chain-rule gradient + Coulomb TERM 1b
    // Reference: Fortran gfnff_engrad.F90:418-422 (bond/disp), 449-454 (coulomb)
    // =========================================================================
    t0 = do_timing ? std::chrono::high_resolution_clock::now() : std::chrono::time_point<std::chrono::high_resolution_clock>{};
    if (gradient && !m_dcn.empty() && m_dcn.natoms == m_natoms) {
        // Snapshot gradient before CN chain-rule (diagnostic)
        m_grad_before_cn = m_result_gradient;

        // CPU PRE-CN GRADIENT (verbosity >= 3)
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info(fmt::format("=== CPU PRE-CN GRADIENT (natoms={}) ===", m_natoms));
            for (int i = 0; i < m_natoms; ++i) {
                fmt::print("  atom {:3d}: {:+.15e} {:+.15e} {:+.15e}\n",
                           i, m_grad_before_cn(i,0), m_grad_before_cn(i,1), m_grad_before_cn(i,2));
            }
            CurcumaLogger::info("=== CPU PRE-CN GRADIENT END ===");
        }

        // Compute TERM 1b qtmp
        Vector qtmp = Vector::Zero(m_natoms);
        bool has_term1b = (m_coulomb_enabled &&   // S1: off with `coulomb false`, as the GPU postprocess
                          m_eeq_charges.size() == m_natoms &&
                          m_cnf.size() == m_natoms &&
                          m_cn.size() == m_natoms);
        if (has_term1b) {
            for (int i = 0; i < m_natoms; ++i) {
                double cn_i = std::max(m_cn(i), 0.0);
                // Claude Generated (Oct 2026): gfnff.cn_zero_guard (opt-in with the reference charge model, always on
                // with a fragment-charge override, GFNFF::cnZeroGuardEffective, Oct 31, 2026). At CN_i = 0 exactly (isolated
                // atom: 0.5*(1+erf(x)) rounded to 0, exp(-x^2) in dcn not) the reference formula gives
                // q*cnf*1e16 times a tiny dcn row, a spurious force. The energy term -q*cnf*sqrt(max(CN,0))
                // is flat there, so its derivative is 0. Exact test: for CN > 0 the formula is the true
                // derivative of the evaluated energy (smallest nonzero log-compressed CN is ~1e-15).
                if (m_cn_zero_guard && cn_i == 0.0) {
                    qtmp(i) = 0.0;
                    continue;
                }
                qtmp(i) = m_eeq_charges(i) * m_cnf(i) / (2.0 * std::sqrt(cn_i) + 1e-16);
            }
        }

        // Combined matvec: gradient += dcn * (dEdcn_total - qtmp)
        Vector dEdcn_combined = has_term1b ? (m_dEdcn_total - qtmp).eval() : m_dEdcn_total;
        // A3c tail (Nov 2026): only with a cell; without one this line is not reached, bitwise as before
        const bool tail_cn = m_disp_tail_active && m_disp_tail_dEdcn.size() == m_natoms;
        if (tail_cn)
            dEdcn_combined += m_disp_tail_dEdcn;
        // CPU dEdcn_combined (verbosity >= 3)
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info(fmt::format("=== CPU dEdcn_combined (natoms={}) ===", m_natoms));
            for (int i = 0; i < m_natoms; ++i) {
                fmt::print("  atom {:3d}: dEdcn={:+.15e} qtmp={:+.15e} comb={:+.15e} cn={:+.15e} cnf={:+.15e}\n",
                           i, m_dEdcn_total(i), (has_term1b ? qtmp(i) : 0.0), dEdcn_combined(i),
                           (m_cn.size() > i ? m_cn(i) : 0.0), (m_cnf.size() > i ? m_cnf(i) : 0.0));
            }
            CurcumaLogger::info("=== CPU dEdcn_combined END ===");
        }
        // Claude Generated (WP4, May 2026): CNDerivStore::applyAdd replaces 3× SpMatrix*v
        m_dcn.applyAdd(dEdcn_combined, m_result_gradient);
        // Per-component CN corrections
        if (m_store_components) {
            Vector dEdcn_disp = m_dEdcn_total - m_dEdcn_bond_total;
            if (tail_cn)
                dEdcn_disp += m_disp_tail_dEdcn;  // tail force counted as dispersion
            m_dcn.applyAdd(m_dEdcn_bond_total, m_result_grad_bond);
            m_dcn.applyAdd(dEdcn_disp, m_result_grad_dispersion);
            if (has_term1b)
                m_dcn.applyAdd(qtmp, m_result_grad_coulomb, -1.0);
        }
        // N1 (Nov 2026): the CN chain rule, the Ewald background and the D4 tail parts of the atomic virial
        if (m_want_virial && m_method_type == FFMethodType::GFN_FF)
            postProcessVirial(has_term1b, qtmp, tail_cn);
    } else if (gradient && m_want_virial && m_method_type == FFMethodType::GFN_FF) {
        postProcessVirial(false, Vector(), false);  // no CN derivative store: background and tail only
    }
    // CPU POST-CN GRADIENT (verbosity >= 3)
    if (gradient) {
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info(fmt::format("=== CPU POST-CN GRADIENT (natoms={}) ===", m_natoms));
            for (int i = 0; i < m_natoms; ++i) {
                fmt::print("  atom {:3d}: {:+.15e} {:+.15e} {:+.15e}\n",
                           i, m_result_gradient(i,0), m_result_gradient(i,1), m_result_gradient(i,2));
            }
            CurcumaLogger::info("=== CPU POST-CN GRADIENT END ===");

            // Per-component CPU gradient decomposition (verbosity >= 3)
            if (m_store_components) {
                const char* names[] = {"REPULSION", "BONDS", "ANGLES", "DISPERSION", "COULOMB", "HB"};
                const GeoGradMatrix* comps[] = {&m_result_grad_repulsion, &m_result_grad_bond, &m_result_grad_angle,
                                                 &m_result_grad_dispersion, &m_result_grad_coulomb, &m_result_grad_hb};
                for (int c = 0; c < 6; ++c) {
                    if (comps[c]->rows() != m_natoms) continue;
                    CurcumaLogger::info(fmt::format("=== CPU {} GRADIENT ===", names[c]));
                    for (int i = 0; i < m_natoms; ++i) {
                        fmt::print("  atom {:3d}: {:+.15e} {:+.15e} {:+.15e}\n",
                                   i, (*comps[c])(i,0), (*comps[c])(i,1), (*comps[c])(i,2));
                    }
                }
                CurcumaLogger::info("=== CPU PER-COMPONENT END ===");
            }
        }
    }

}

// ============================================================================
// Atomic virial: the parts that are not a kernel sum (Claude Generated, Nov 2026, WP-PERIODIC-NONBONDED N1)
//
// CN chain rule: the gradient of every CN-dependent term is G_i += sum_j v_j dCN_j/dx_i (CNDerivStore::applyAdd,
// entries (i, j, c) with c = dCN_j/dx_i times dlogCN/dcn, plus the diagonal dCN_j/dx_j = -sum_i c). Under a strain
// the pair vector d = x_i - x_j (minimum image, the one the derivative was built with) scales affinely, so
//   Xi_ab -= sum_entries v_j d_a c_b
// (the diagonal entry is the negative row sum and drops out, as the frame origin of every term does). v is the
// same dEdcn_combined that the gradient uses, split as the gradient components are: bond (dr0/dCN), dispersion
// (dC6/dCN and the D4 tail), Coulomb (TERM 1b, -qtmp).
// Ewald background terms ~ 1/V: Xi_ab = (E_bg + E_gbg) delta_ab; the self term is strain-independent.
// D4 tail at fixed r_c (A3c): E_tail ~ 1/V, so Xi_ab = E_tail delta_ab; its dE/dCN is in the chain rule above.
//
// V-LRC long-range pressure correction (Claude Generated, Nov 2026; the classical tail correction of the virial
// pressure for a truncated pair potential, e.g. Allen and Tildesley, Computer Simulation of Liquids, section on
// long-range corrections). For one pair class u(r) = -C_n r^-n, N atoms, rho = N/V, g(r) = 1 beyond r_c:
//   energy tail    E_n = (N rho/2) int_rc^inf 4 pi r^2 u dr       = -2 pi N rho C_n / ((n-3) r_c^(n-3))
//   virial tail    W_n = -(N rho/2) int_rc^inf 4 pi r^2 r u'(r) dr = n E_n   (r u' = -n u)
//   pressure       P_n = W_n / (3V) = (n/3) E_n / V.
// The model derivative above (fixed r_c, fixed pair set: E_n ~ 1/V) gives E_n/V, factor 1 instead of n/3. The
// missing (n-3)/3 E_n/V is the impulsive term of the pairs that cross the sharp radius r_c under the strain
// (energy jump u(r_c) per crossing pair; it equals -(2 pi/3) rho^2 r_c^3 u(r_c) for g = 1). A finite-difference
// strain at constant pair count does not see it. Trace of the correction: (n-3) E_n, i.e. with the two parts of
// the D4 tail (C6: n = 6, C8: n = 8)  tr Xi_lrc = 3 E6_tail + 5 E8_tail, isotropic Xi_lrc = tr/3 * I.
// This is a statistical estimator in the homogeneous approximation (g = 1 beyond r_c, large-r form of the damping
// as the tail energy), NOT a derivative of the model energy; kept in its own component tail_lrc.
// ============================================================================
void FFWorkspace::postProcessVirial(bool has_term1b, const Vector& qtmp, bool tail_cn)
{
    if (ewaldActive() && m_coulomb_enabled)
        m_result_virial.coulomb_bg += curcuma::ewald::EwaldSum::backgroundVirial(m_ewald_terms);
    if (m_disp_tail_active)
        m_result_virial.disp_tail += m_disp_tail_energy * Eigen::Matrix3d::Identity();
    if (m_disp_tail_active && m_disp_tail_lrc)
        m_result_virial.tail_lrc += ((3.0 * m_disp_tail_e6 + 5.0 * m_disp_tail_e8) / 3.0) * Eigen::Matrix3d::Identity();
    if (m_dcn.empty() || m_dcn.natoms != m_natoms || m_dEdcn_total.size() != m_natoms)
        return;
    const bool pbc = m_cell_bohr.active;
    Vector v_bond = m_dEdcn_bond_total;
    Vector v_disp = m_dEdcn_total - m_dEdcn_bond_total;
    if (tail_cn)
        v_disp += m_disp_tail_dEdcn;
    Vector v_coul = has_term1b ? Vector(-qtmp) : Vector(Vector::Zero(m_natoms));
    Eigen::Matrix3d xb = Eigen::Matrix3d::Zero(), xd = Eigen::Matrix3d::Zero(), xc = Eigen::Matrix3d::Zero();
    for (const auto& e : m_dcn.pairs) {
        double d0 = m_geometry(e.i, 0) - m_geometry(e.j, 0);
        double d1 = m_geometry(e.i, 1) - m_geometry(e.j, 1);
        double d2 = m_geometry(e.i, 2) - m_geometry(e.j, 2);
        if (pbc)
            PBCUtils::applyMinimumImageOrtho(m_cell_bohr, d0, d1, d2);
        const Eigen::Vector3d d(d0, d1, d2);
        const Eigen::Vector3d c(e.cx, e.cy, e.cz);
        const Eigen::Matrix3d outer = d * c.transpose();
        xb -= v_bond(e.j) * outer;
        xd -= v_disp(e.j) * outer;
        xc -= v_coul(e.j) * outer;
    }
    m_result_virial.cn_bond += xb;
    m_result_virial.cn_dispersion += xd;
    m_result_virial.cn_coulomb += xc;
}
