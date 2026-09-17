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
#include "cn_calculator.h"
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

    // UFF/QMDFF non-bonded pairs
    m_vdws = std::move(params.vdws);

    m_eeq_charges = std::move(params.eeq_charges);
    m_topology_charges = std::move(params.topology_charges);

    m_e0 = params.e0;

    // Method type and distance unit factor
    m_method_type = params.method_type;
    m_au = (m_method_type != FFMethodType::GFN_FF) ? 1.889726125 : 1.0;

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
    // rev-gfnff stage 3a(i) (Claude Generated, Sep 2026): the CN radii the pair correction in
    // calcBonds() needs. Built here (not lazily in a kernel) so that no worker thread writes it.
    m_rev_cn_rcov.resize(atoms.size());
    for (size_t i = 0; i < atoms.size(); ++i)
        m_rev_cn_rcov[i] = CNCalculator::gfnffCNRadiusBohr(atoms[i]);
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

double FFWorkspace::calculateSingle(bool gradient)
{
    const bool do_timing = (CurcumaLogger::get_verbosity() >= 2);
    auto t_calc_start = do_timing ? std::chrono::high_resolution_clock::now() : std::chrono::time_point<std::chrono::high_resolution_clock>{};
    double t_execute = 0.0, t_reduce = 0.0, t_post = 0.0;

    m_do_gradient = gradient;

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

    // rev-gfnff stage 3a(ii) (Sep 2026): the valence share of a bond needs the atom's complete
    // bond-order sum, so it is built once per corner on the main thread, before the partitions.
    if (m_rev.enabled && m_rev.valence_share && m_method_type == FFMethodType::GFN_FF)
        prepareValenceShare();
    // rev-gfnff stage 3a(iii) (Sep 2026): the per-bond well-form parameters of this corner. Same
    // place and the same reason - the erf-Morse form needs a bisection for its offset, which must
    // not run inside the per-partition energy loop.
    if (m_rev.enabled && m_rev.well_form != 0 && m_method_type == FFMethodType::GFN_FF)
        prepareWellForms();

    auto t0 = do_timing ? std::chrono::high_resolution_clock::now() : std::chrono::time_point<std::chrono::high_resolution_clock>{};
    if (m_num_threads == 1) {
        // T=1: Direct call, zero pool overhead
        m_accumulators[0].reset(m_natoms, gradient, m_store_components);
        executeMethod(0);

        // acc[0] IS the result — zero-copy swap
        m_result_energy = m_accumulators[0].energy;
        m_result_timings = m_accumulators[0].timings;  // Claude Generated (May 2026): per-term timing
        if (gradient) {
            m_result_gradient.swap(m_accumulators[0].gradient);
            m_dEdcn_total = m_accumulators[0].dEdcn;
            m_dEdcn_bond_total = m_accumulators[0].dEdcn_bond;
            m_dEdshare_total = m_accumulators[0].dEdshare;
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

    // rev-gfnff stage 1 (Sep 2026): the over-coordination energy needs the complete
    // bond-order sums, so it runs once on the main thread after the partitions.
    if (m_rev.enabled && m_rev.over_coord && m_method_type == FFMethodType::GFN_FF)
        calcOverCoordination(gradient);
    // rev-gfnff stage 3a(ii) (Sep 2026): the valence share's own chain rule - the share factor of
    // every bond depends on the sum over both atoms' OTHER pairs, so each pair also carries the
    // derivative of those. Same place as E_over: after the complete sums exist.
    if (gradient && m_rev.enabled && m_rev.valence_share && m_method_type == FFMethodType::GFN_FF)
        applyValenceShareGradient();
    // rev-gfnff stage 2 (Sep 2026): the bond-hardness term of the split-charge model. Same
    // place as E_over: it is a per-pair term over the whole corner, not a partitioned list.
    if (m_rev.enabled && !m_sqe_pairs.empty() && m_method_type == FFMethodType::GFN_FF)
        calcSqeHardness(gradient);

    t0 = do_timing ? std::chrono::high_resolution_clock::now() : std::chrono::time_point<std::chrono::high_resolution_clock>{};
    postProcess(gradient);
    if (do_timing) {
        t_post = std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - t0).count();
    }

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
        fmt::print("  overcoord = {:+.15e}\n", m_result_energy.over_coord);
        fmt::print("  sqe_hard  = {:+.15e}\n", m_result_energy.sqe_hardness);
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

    // Bonded terms. rev-gfnff stage 1b (Sep 2026): during a topology transition both
    // topologies are evaluated and blended, E = (1 - s) E_alt + s E_primary.
    auto t = tic();
    runBonded(m_accumulators[p], m_partitions[p], &timings);
    (void)t;

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
        t = tic(); calcCoulomb(p); timings.coulomb = toc(t);
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

    // Triple bond torsions: evaluated inside runBonded() with the other bonded terms.
}

// rev-gfnff stage 1b (Claude Generated, Sep 2026): the six bonded kernels of one topology
void FFWorkspace::runBonded(FFAccumulator& acc, const PartitionRanges& pr, FFTermTimings* timings)
{
    auto tic = []() { return std::chrono::high_resolution_clock::now(); };
    auto toc = [](auto t0) { return std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - t0).count(); };
    const auto& bonds = m_bonds;
    const auto& angles = m_angles;
    const auto& dihedrals = m_dihedrals;
    const auto& extra = m_extra_dihedrals;
    const auto& inversions = m_inversions;
    const auto& storsions = m_storsions;
    auto t = tic(); calcBonds(acc, bonds, pr.bonds);            if (timings) timings->bonds = toc(t);
    t = tic();      calcAngles(acc, angles, pr.angles);         if (timings) timings->angles = toc(t);
    t = tic();      calcDihedrals(acc, dihedrals, pr.dihedrals);
                    calcExtraTorsions(acc, extra, pr.extra_dihedrals); if (timings) timings->dihedrals = toc(t);
    t = tic();      calcInversions(acc, inversions, pr.inversions); if (timings) timings->inversions = toc(t);
    if (!storsions.empty()) {
        t = tic();  calcSTorsions(acc, storsions, pr.storsions); if (timings) timings->stors = toc(t);
    }
}



void FFWorkspace::reduce()
{
    // Sum all accumulators into result
    m_result_energy = m_accumulators[0].energy;
    m_result_timings = m_accumulators[0].timings;  // Claude Generated (May 2026)
    if (m_do_gradient) {
        m_result_gradient = m_accumulators[0].gradient;
        m_dEdcn_total = m_accumulators[0].dEdcn;
        m_dEdcn_bond_total = m_accumulators[0].dEdcn_bond;
        // Claude Generated (Sep 14, 2026): dEdshare was never reduced, so with more than one
        // thread the valence share's chain rule (applyValenceShareGradient) read a stale or empty
        // m_dEdshare_total - the share gradient was silently thread-count dependent. The T=1 path
        // always took the accumulator directly, which is why only the threaded path was wrong.
        m_dEdshare_total = m_accumulators[0].dEdshare;
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
        if (m_do_gradient) {
            m_result_gradient += m_accumulators[t].gradient;
            m_dEdcn_total += m_accumulators[t].dEdcn;
            m_dEdcn_bond_total += m_accumulators[t].dEdcn_bond;
            if (m_dEdshare_total.size() == m_accumulators[t].dEdshare.size())
                m_dEdshare_total += m_accumulators[t].dEdshare;
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
    if (do_timing) {
        t_self_energy = std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - t0).count();
    }


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
            m_dcn.applyAdd(m_dEdcn_bond_total, m_result_grad_bond);
            m_dcn.applyAdd(dEdcn_disp, m_result_grad_dispersion);
            if (has_term1b)
                m_dcn.applyAdd(qtmp, m_result_grad_coulomb, -1.0);
        }
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

// ============================================================================================
// rev-gfnff stage 1b (Claude Generated, Sep 12, 2026): multi-transition blending over 2^k corners
// ============================================================================================

void FFWorkspace::swapState(TopologyState& st)
{
    std::swap(m_bonds, st.bonds);
    std::swap(m_angles, st.angles);
    std::swap(m_dihedrals, st.dihedrals);
    std::swap(m_extra_dihedrals, st.extra_dihedrals);
    std::swap(m_inversions, st.inversions);
    std::swap(m_storsions, st.storsions);
    std::swap(m_dispersions, st.dispersions);
    std::swap(m_d4_dispersions, st.d4_dispersions);
    std::swap(m_bonded_reps, st.bonded_reps);
    std::swap(m_nonbonded_reps, st.nonbonded_reps);
    std::swap(m_coulombs, st.coulombs);
    std::swap(m_hbonds, st.hbonds);
    std::swap(m_xbonds, st.xbonds);
    std::swap(m_atm_triples, st.atm_triples);
    std::swap(m_batm_triples, st.batm_triples);
    std::swap(m_bond_hb_data, st.bond_hb_data);
    std::swap(m_hb_grad_entries, st.hb_grad_entries);
    std::swap(m_hb_grad_offsets, st.hb_grad_offsets);
    std::swap(m_hb_grad_list, st.hb_grad_list);
    std::swap(m_bonded_pairs, st.bonded_pairs);
    std::swap(m_partitions, st.partitions);
    std::swap(m_eeq_charges, st.eeq_charges);
    std::swap(m_topology_charges, st.topology_charges);
    std::swap(m_coul_chi_base, st.coul_chi_base);
    std::swap(m_coul_gam, st.coul_gam);
    std::swap(m_coul_alp, st.coul_alp);
    std::swap(m_coul_cnf, st.coul_cnf);
    std::swap(m_coul_chi_static, st.coul_chi_static);
    std::swap(m_sqe_pairs, st.sqe_pairs); // rev-gfnff stage 2: the corner's split charges travel with it
    std::swap(m_e0, st.e0);
}

void FFWorkspace::updateTransitions()
{
    for (RevTransition& tr : m_transitions) {
        if (!tr.active || m_geometry.rows() <= std::max(tr.i, tr.j))
            continue;
        Eigen::Vector3d ri = m_geometry.row(tr.i), rj = m_geometry.row(tr.j);
        const double r = (ri - rj).norm();
        double dw = 0.0, dc = 0.0;
        const double w = revWeightRaw(tr.i, tr.j, r, &dw); // the scan's thresholds are on the raw switch
        const double c = tr.tight ? revOrder(tr.i, tr.j, r, &dc) : revCoord(tr.i, tr.j, r, &dc);
        // smoothstep in the window: C1 at both ends (the linear clamp had a kink at s = 0 and
        // s = 1, i.e. a force that switches on abruptly at the detected coordinate)
        const double span = tr.w_b - tr.w_a;
        double u = (c - tr.w_a) / span;
        double s, dsdw;
        if (u <= 0.0) { s = 0.0; dsdw = 0.0; }
        else if (u >= 1.0) { s = 1.0; dsdw = 0.0; }
        else { s = u * u * (3.0 - 2.0 * u); dsdw = 6.0 * u * (1.0 - u) / span; }
        tr.r = r;
        tr.w = w;
        tr.c = c;
        tr.dwdr = dc;
        tr.s = s;
        tr.dsdw = dsdw;
    }
}

void FFWorkspace::beginTransition(RevTransition tr, std::vector<GFNFFParameterSet>&& params)
{
    const int k = static_cast<int>(m_transitions.size());
    const int n_old = 1 << k;
    if (static_cast<int>(m_corners.size()) != n_old) { // k == 0: the slot is the only corner
        m_corners.assign(n_old, TopologyState {});
        m_slot_mask = n_old - 1;
    }
    swapState(m_corners[m_slot_mask]); // the slot now holds the placeholder
    std::vector<TopologyState> corners(2 * n_old);
    for (int m = 0; m < n_old; ++m)
        corners[m] = std::move(m_corners[m]);
    for (int m = 0; m < n_old; ++m) {
        rebuildInteractionLists(std::move(params[m])); // into the empty slot
        swapState(corners[m | (1 << k)]);              // and out again
    }
    if (tr.forming && !tr.well_blend) {
        // the well of the forming pair belongs to every corner, so the blend never touches it
        // (only correct while the join sits where the term weight is ~0; with rev_form_switch =
        // order the join is tight, the well is NOT copied and the blend ramps it in over s)
        const int a = std::min(tr.i, tr.j), b = std::max(tr.i, tr.j);
        auto find_pair = [&](const std::vector<Bond>& list) -> const Bond* {
            for (const auto& bd : list)
                if (std::min(bd.i, bd.j) == a && std::max(bd.i, bd.j) == b)
                    return &bd;
            return nullptr;
        };
        for (int m = 0; m < n_old; ++m) {
            if (find_pair(corners[m].bonds))
                continue;
            if (const Bond* bd = find_pair(corners[m | (1 << k)].bonds)) {
                corners[m].bonds.push_back(*bd);
                swapState(corners[m]);
                partition();
                swapState(corners[m]);
            }
        }
    }
    m_corners = std::move(corners);
    tr.active = true;
    m_transitions.push_back(tr);
    m_slot_mask = (1 << (k + 1)) - 1;
    swapState(m_corners[m_slot_mask]);
    updateTransitions();
}

void FFWorkspace::endTransition(int t, bool keep_new)
{
    const int k = static_cast<int>(m_transitions.size());
    if (t < 0 || t >= k)
        return;
    const int n = 1 << k;
    swapState(m_corners[m_slot_mask]); // store the slot
    std::vector<TopologyState> kept(n / 2);
    for (int mask = 0; mask < n; ++mask) {
        if (((mask >> t) & 1) != (keep_new ? 1 : 0))
            continue;
        const int nm = (mask & ((1 << t) - 1)) | ((mask >> (t + 1)) << t);
        kept[nm] = std::move(m_corners[mask]);
    }
    m_corners = std::move(kept);
    m_transitions.erase(m_transitions.begin() + t);
    m_slot_mask = (1 << (k - 1)) - 1;
    swapState(m_corners[m_slot_mask]);
    if (m_transitions.empty()) {
        m_corners.clear();
        m_slot_mask = 0;
    }
    updateTransitions();
}

static void addScaledComponents(FFEnergyComponents& a, const FFEnergyComponents& b, double w)
{
    a.bond += w * b.bond; a.angle += w * b.angle; a.dihedral += w * b.dihedral; a.inversion += w * b.inversion;
    a.dispersion += w * b.dispersion; a.vdw += w * b.vdw; a.rep += w * b.rep;
    a.bonded_rep += w * b.bonded_rep; a.nonbonded_rep += w * b.nonbonded_rep;
    a.coulomb += w * b.coulomb; a.hbond += w * b.hbond; a.xbond += w * b.xbond;
    a.atm += w * b.atm; a.batm += w * b.batm; a.stors += w * b.stors; a.over_coord += w * b.over_coord;
    a.sqe_hardness += w * b.sqe_hardness;
    a.hbond_case1 += w * b.hbond_case1; a.hbond_case2 += w * b.hbond_case2; a.hbond_case3 += w * b.hbond_case3; a.hbond_case4 += w * b.hbond_case4;
}

double FFWorkspace::calculate(bool gradient)
{
    if (m_transitions.empty())
        return calculateSingle(gradient);
    updateTransitions();
    const int k = static_cast<int>(m_transitions.size());
    const int n = 1 << k;
    m_corner_energy.assign(n, 0.0);
    m_corner_gradient.resize(n);
    m_corner_components.assign(n, FFEnergyComponents {});
    std::vector<Vector> dcn(n), dcnb(n);
    for (int mask = 0; mask < n; ++mask) {
        const bool swap_in = (mask != m_slot_mask);
        if (swap_in)
            swapState(m_corners[mask]);
        if (m_corner_prepare)
            m_corner_prepare(mask);
        m_corner_energy[mask] = calculateSingle(gradient);
        m_corner_components[mask] = m_result_energy;
        if (gradient) {
            m_corner_gradient[mask] = m_result_gradient;
            dcn[mask] = m_dEdcn_total;
            dcnb[mask] = m_dEdcn_bond_total;
        }
        if (swap_in)
            swapState(m_corners[mask]);
    }
    auto weight = [&](int mask, int skip) {
        double w = 1.0;
        for (int t = 0; t < k; ++t) {
            if (t == skip)
                continue;
            w *= ((mask >> t) & 1) ? m_transitions[t].s : (1.0 - m_transitions[t].s);
        }
        return w;
    };
    double E = 0.0;
    FFEnergyComponents comp {};
    GeoGradMatrix G;
    Vector dcn_tot, dcnb_tot;
    if (gradient) {
        G = GeoGradMatrix::Zero(m_natoms, 3);
        if (dcn[0].size() > 0) dcn_tot = Vector::Zero(dcn[0].size());
        if (dcnb[0].size() > 0) dcnb_tot = Vector::Zero(dcnb[0].size());
    }
    for (int mask = 0; mask < n; ++mask) {
        const double W = weight(mask, -1);
        E += W * m_corner_energy[mask];
        addScaledComponents(comp, m_corner_components[mask], W);
        if (gradient) {
            G += W * m_corner_gradient[mask];
            if (dcn_tot.size() == dcn[mask].size()) dcn_tot += W * dcn[mask];
            if (dcnb_tot.size() == dcnb[mask].size()) dcnb_tot += W * dcnb[mask];
        }
    }
    if (gradient) {
        // d E / d s_t  (d s_t / d r_t) on the transition pair
        for (int t = 0; t < k; ++t) {
            const RevTransition& tr = m_transitions[t];
            if (tr.dsdw == 0.0 || tr.r <= 1e-8)
                continue;
            double D = 0.0;
            for (int mask = 0; mask < n; ++mask)
                D += (((mask >> t) & 1) ? 1.0 : -1.0) * weight(mask, t) * m_corner_energy[mask];
            const double f = D * tr.dsdw * tr.dwdr / tr.r;
            Eigen::Vector3d ri = m_geometry.row(tr.i), rj = m_geometry.row(tr.j);
            Eigen::Vector3d g = f * (ri - rj);
            G.row(tr.i) += g.transpose();
            G.row(tr.j) -= g.transpose();
        }
        m_result_gradient = G;
        m_dEdcn_total = dcn_tot;
        m_dEdcn_bond_total = dcnb_tot;
    }
    m_result_energy = comp;
    // the corners' e0 differ; callers read m_e0 + total(), so the blended offset goes into the Coulomb field
    m_result_energy.coulomb += E - (m_e0 + m_result_energy.total());
    return E;
}
