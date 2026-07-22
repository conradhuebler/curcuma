/*
 * <Conformational Search based on Molecular Dynamics>
 * Copyright (C) 2022 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "src/global_config.h"

#include "src/capabilities/confscan.h"
#include "src/capabilities/optimizer_factory.h"
#include "src/capabilities/simplemd.h"

#include "src/core/fileiterator.h"
#include "src/core/molecule.h"

#include "src/tools/general.h"
#include "src/core/parameter_registry.h"
#include "src/capabilities/rmsd/rmsd_functions.h"  // Claude Generated (Jun 2026): Kabsch helpers for bias calibration

#include <algorithm>
#include <cmath>

#include "external/CxxThreadPool/include/CxxThreadPool.hpp"

#include <fstream>
#include <iostream>
#include <stdio.h>

#include "confsearch.h"
using curcuma::Molecule;

ConfSearch::ConfSearch(const json& controller, bool silent)
    : CurcumaMethod(ParameterRegistry::getInstance().getDefaultJson("confsearch"), controller, silent)
    , m_config("confsearch", controller)
{
    UpdateController(controller);
}

ConfSearch::~ConfSearch()
{
    delete m_bias_pool;
}

void ConfSearch::setFile(const std::string& filename)
{
    CurcumaMethod::setFile(filename);

    FileIterator file(Filename());
    while (!file.AtEnd()) {
        Molecule* mol = new Molecule(file.Next());
        // Claude Generated (Jul 2026): seed the requested charge/spin onto the molecule. This is
        // the ONLY channel that reaches a QM method (EnergyCalculator reads Mol::m_charge, never
        // the controller), and an XYZ file carries no charge unless curcuma itself wrote it.
        // Matches curcumaopt.cpp / simplemd.cpp / main.cpp -sp. LoadControlJson has already run
        // (ctor -> UpdateController), so m_charge/m_spin are valid here.
        mol->setCharge(m_charge);
        mol->setSpin(m_spin);
        m_in_stack.push_back(mol);
        m_topo_matrix = mol->DistanceMatrix().second;
    }
}

bool ConfSearch::Initialise()
{
    return true;
}

void ConfSearch::start()
{
    const std::string p = Basename();

    // Claude Generated (Jul 2026): the MD configuration is layered, most general first.
    // ConfSearch now owns its parameter names in the registry, so the flat CLI flags land in
    // controller["confsearch"] and the old hand-picked rescue from controller["simplemd"] is gone.
    //
    // Layer 1: SimpleMD's own registry defaults.
    nlohmann::json md = ParameterRegistry::getInstance().getDefaultJson("simplemd");

    // Layer 2: everything the user aimed explicitly at the MD engine (-md.x / -simplemd.x).
    // Merged wholesale now -- previously only 8 hand-listed keys were rescued from here, so
    // -thermostat / -coupling / -rattle / -wall_* were silently dropped.
    if (m_controller.contains("simplemd") && m_controller["simplemd"].is_object())
        md.merge_patch(m_controller["simplemd"]);

    // Layer 3: global section (gpu, threads, charge, spin, verbosity, method).
    if (m_controller.contains("global") && m_controller["global"].is_object()) {
        for (auto& [key, value] : m_controller["global"].items()) {
            if (!value.is_object() && key != "confsearch" && key != "confscan")
                md[key] = value;
        }
    }

    // Layer 4: values the user set on the ConfSearch command itself. Most specific to the active
    // command, so it wins over Layer 2 ("-md.thermostat berendsen -thermostat csvr" -> csvr).
    // Nested objects (method sub-scopes) are handled by ChildConfig in Layer 5.
    if (m_controller.contains("confsearch") && m_controller["confsearch"].is_object()) {
        for (auto& [key, value] : m_controller["confsearch"].items()) {
            if (value.is_object())
                continue;
            md[key] = value;
        }
    }

    // Layer 5: system identity + runtime + method sub-scopes, shared with every other child.
    // When ConfSearch parallelises cycles (m_threads > 1) each MD must stay single-threaded to
    // avoid nested CxxThreadPools.
    md.merge_patch(ChildConfig(m_md_method, (m_threads > 1) ? 1 : m_threads));

    // GPU + Multi-Threading Safety: Deactivate GPU when threads > 1 to prevent
    // GPU contention. Multiple MD instances cannot share the GPU simultaneously.
    if (m_threads > 1 && md.contains("gpu") && !md["gpu"].is_null() && md["gpu"] != "none") {
        CurcumaLogger::warn("GPU cannot be used with multiple threads simultaneously. Disabling GPU for this run.");
        md["gpu"] = "none";
        m_gpu = "none";
    }

    // Layer 6: ConfSearch's non-negotiable overrides.
    md["unique"] = true;
    md["rmsd"] = m_rmsd;
    md["time_step"] = m_dT;
    md["max_time"] = m_time;
    md["restart"] = false;       // ConfSearch manages its own state, no MD restart
    md["norestart"] = true;

    // Robustness controls: these self-reference inside SimpleMD (start fragment count / start
    // energy / target T / per-run inherited pool), so only the enable flags and windows are
    // forwarded. temp_abort and rmsd_mtd_freeze_inherited default to ON in ConfSearch but OFF in
    // SimpleMD -- pinned here so the divergence survives every layer above.
    md["topo_check"] = m_topo_check;
    md["topo_check_interval"] = m_topo_check_interval;
    md["epot_abort"] = m_epot_abort;
    md["epot_abort_window"] = m_epot_abort_window;
    md["temp_abort"] = m_temp_abort;
    md["temp_abort_factor"] = m_temp_abort_factor;
    md["temp_abort_delta"] = m_temp_abort_delta;
    md["rmsd_mtd_max_height"] = m_rmsd_mtd_max_height;
    md["rmsd_mtd_freeze_inherited"] = m_freeze_inherited;

    // RMSD metadynamics is the default driver for conformational exploration.
    // The SimpleMD default is false, but ConfSearch enables it by default.
    // Only disable if the user explicitly passed -rmsd_mtd false.
    if (!(m_controller.contains("confsearch") && m_controller["confsearch"].contains("rmsd_mtd"))
        && !(m_controller.contains("simplemd") && m_controller["simplemd"].contains("rmsd_mtd"))
        && !m_controller.contains("rmsd_mtd"))
        md["rmsd_mtd"] = true;

    // Log ConfSearch configuration (visible at verbosity >= 1)
    CurcumaLogger::result_fmt("ConfSearch: MD method={}, Opt method={}, Thermostat={}, Threads={}", m_md_method, m_opt_method, m_thermostat, m_threads);
    // Claude Generated (Jun 2026): explicit phase-method mapping so the user can verify their choices took effect
    if (m_opt_method != m_md_method) {
        CurcumaLogger::result_fmt("ConfSearch: Dual-method mode (explore={}, refine={})", m_md_method, m_opt_method);
        CurcumaLogger::result_fmt("ConfSearch: Phase methods: explore={}, pre-opt={}, refine={}, final rank={}",
            m_md_method, m_md_method, m_opt_method, m_opt_method);
    } else {
        CurcumaLogger::result("ConfSearch: Single-method mode (Phase 3b skipped)");
    }
    CurcumaLogger::result_fmt("ConfSearch: Temperature={}K -> {}K, step={}K", m_startT, m_endT, m_deltaT);
    CurcumaLogger::result_fmt("ConfSearch: Repetitions={}, RMSD threshold={} A, Energy window={} kJ/mol, Seed rank={}", m_repeat, m_rmsd, m_energy_window, m_seed_rank);
    // Debug: log full MD parameter set for diagnosing dynamics issues
    CurcumaLogger::result_fmt("ConfSearch MD config: temperature={}, T={}, impuls={}, time_step={}, max_time={}, rmsd_mtd={}",
        md.value("temperature", -1.0), md.value("T", -1.0), md.value("impuls", -1.0),
        md.value("time_step", -1.0), md.value("max_time", -1.0), md.value("rmsd_mtd", false));
    CurcumaLogger::result_fmt("ConfSearch MD config: method={}, thermostat={}, coupling={}, remove_com={}, remove_com_mode={}, no_center={}",
        md.value("method", "?"), md.value("thermostat", "?"), md.value("coupling", -1.0),
        md.value("remove_com_motion", -1.0), md.value("remove_com_mode", -1), md.value("no_center", false));
    // Claude Generated (Jun 2026): adaptive MTD hill width (experimental, opt-in). The default
    // alpha=10 (half-max ~0.26 A) is far narrower than the dedup RMSD scale (m_rmsd ~1.25 A), so
    // filling a basin needs many hills and exploration is slow. "couple" sets the hill half-max at
    // bias_couple_factor*m_rmsd -> alpha = ln2/(factor*rmsd)^2, so one hill ~ one basin. RMSD and
    // m_rmsd are both in Angstrom (m_eigen_geometry is set without unit conversion) -> alpha in A^-2.
    // Applied before the dump/log below so they reflect the calibrated value.
    // "couple" sets alpha once from the dedup threshold. "cluster" uses the same formula only as
    // the cycle-1 bootstrap (no cluster data yet) and then relearns alpha each cycle in CalibrateBias.
    if (m_bias_calibration == "couple" || m_bias_calibration == "cluster") {
        double r = m_bias_couple_factor * m_rmsd;
        if (r > 1e-6) {
            double alpha = std::log(2.0) / (r * r);
            md["rmsd_mtd_alpha"] = alpha;
            CurcumaLogger::result_fmt("ConfSearch: bias_calibration={} -> initial RMSD-MTD alpha={:.4f} A^-2 (hill half-max at {:.3f} A)",
                m_bias_calibration, alpha, r);
        }
    } else if (m_bias_calibration != "off") {
        CurcumaLogger::warn_fmt("ConfSearch: bias_calibration='{}' unknown -- treating as off", m_bias_calibration);
    }

    // Dump all MD parameters to file for debugging
    {
        std::ofstream debug_file(outputPath(p + "_md_params.json"));
        debug_file << md.dump(2) << std::endl;
        debug_file.close();
        CurcumaLogger::result_fmt("ConfSearch: Full MD parameters written to {}_md_params.json", p);

        // Claude Generated (Jul 2026): the configs every non-MD child gets. Makes it verifiable
        // that charge/spin/gpu and the method sub-scopes actually reach the optimisations and
        // the ConfScan filters, which is exactly what the old ad-hoc {method,threads,gpu} JSONs
        // silently dropped.
        nlohmann::json children;
        children["opt_md"] = ChildConfig(m_md_method, m_threads);
        children["opt"] = ChildConfig(m_opt_method, m_threads);
        children["filter"] = FilterConfig(m_opt_method, m_threads);
        std::ofstream child_file(outputPath(p + "_child_params.json"));
        child_file << children.dump(2) << std::endl;
        child_file.close();
        CurcumaLogger::result_fmt("ConfSearch: Child computation parameters written to {}_child_params.json", p);
    }

    if (md.value("rmsd_mtd", false)) {
        CurcumaLogger::result("ConfSearch: RMSD-MTD Enabled");
        double k = md.value("rmsd_mtd_k", 0.1);
        double alpha = md.value("rmsd_mtd_alpha", 10.0);
        int pace = md.value("rmsd_mtd_pace", 1);
        bool wtmtd = md.value("wtmtd", false);
        CurcumaLogger::result_fmt("ConfSearch: RMSD-MTD k={} Eh, alpha={} A^-2, pace={} steps", k, alpha, pace);
        if (wtmtd)
            CurcumaLogger::result_fmt("ConfSearch: RMSD-MTD Well-tempered (dT={})", md.value("rmsd_mtd_dt", 1000000.0));
        else
            CurcumaLogger::result("ConfSearch: RMSD-MTD Well-tempered Off");
    } else {
        CurcumaLogger::result("ConfSearch: RMSD-MTD Disabled");
    }

    // Claude Generated (Jun 2026): restart. Try to resume from a checkpoint in the start directory.
    // On a successful resume we restore the search state and SKIP the pre-optimisation below.
    bool resumed = false;
    int pending_entry = 0; // entry phase for the FIRST resumed cycle (0=md, 1=post_md); reset to 0 after
    if (m_do_restart && loadCheckpoint()) {
        resumed = true;
        pending_entry = (m_restart.entry_phase >= 1) ? 1 : 0; // v1: only md / post_md resume granularity
        m_elements = m_restart.elements;
        m_topo_ref = m_restart.topo_ref;
        m_topo_matrix = m_topo_ref.DistanceMatrix().second;
        m_global_min = m_restart.global_min;
        m_initial_energy = m_restart.initial_energy;
        m_best_energy = m_restart.best_energy;
        m_permutation_cache = m_restart.permutations;
        for (auto* mol : m_in_stack) delete mol;
        m_in_stack.clear();
        for (const auto& mm : m_restart.seeds)
            m_in_stack.push_back(new Molecule(mm));
        CurcumaLogger::header("=== ConfSearch: RESUMING from checkpoint ===");
        CurcumaLogger::result_fmt("ConfSearch: resume T={}K, {} cycles done, {} bias structures, {} seeds, {} cumulative conformers",
            static_cast<int>(m_restart.next_T), m_restart.temperature_cycle,
            static_cast<int>(m_restart.bias.size()), static_cast<int>(m_in_stack.size()),
            static_cast<int>(m_restart.cumulative.size()));
    }

    // Optimise all input structures before any MD run (skipped on resume).
    // m_topo_matrix is updated from the first optimised structure so Phase 4
    // topology checks compare against the relaxed geometry, not the raw input.
    if (!resumed) {
        CurcumaLogger::header("=== ConfSearch: Initial Geometry Optimisation ===");
        bool first = true;
        for (auto* mol : m_in_stack) {
            if (first) { mol->writeXYZFile(outputPath(p + ".input.xyz")); first = false; }
            else          mol->appendXYZFile(outputPath(p + ".input.xyz"));
        }
        // pre-optimization at md_method ("die md-methode macht die voroptimierung")
        nlohmann::json opt_init = ChildConfig(m_md_method, m_threads);
        PerformOptimisation(p + ".input", opt_init);

        for (auto* mol : m_in_stack) delete mol;
        m_in_stack.clear();
        FileIterator opt_in(outputPath(p + ".input.opt.xyz"));
        while (!opt_in.AtEnd()) {
            Molecule mol = opt_in.Next();
            mol.setCharge(m_charge);
            mol.setSpin(m_spin);
            if (mol.AtomCount() > 0)
                m_in_stack.push_back(new Molecule(mol));
        }
        if (!m_in_stack.empty()) {
            m_topo_matrix = m_in_stack[0]->DistanceMatrix().second;
            m_topo_ref = *m_in_stack[0];           // reference structure for restart topology
            m_elements = m_in_stack[0]->Atoms();   // shared atomic-number list for checkpoint frames
        }
        CurcumaLogger::result_fmt("ConfSearch: {} input structures optimised", m_in_stack.size());

        // Claude Generated (Jun 2026): dual initial optimization -- in dual-method mode,
        // also optimise the md_method-minimised structures at opt_method so we can track
        // the energy gain on both PES from the very start. The opt_method structures are
        // for reporting only; m_in_stack keeps the md_method structures (they feed the MD loop).
        if (m_opt_method != m_md_method) {
            CurcumaLogger::header("=== ConfSearch: Initial Geometry Optimisation (" + m_opt_method + ") ===");
            bool first_opt = true;
            for (auto* mol : m_in_stack) {
                if (first_opt) { mol->writeXYZFile(outputPath(p + ".input_mdopt.xyz")); first_opt = false; }
                else              mol->appendXYZFile(outputPath(p + ".input_mdopt.xyz"));
            }
            nlohmann::json opt_hi = ChildConfig(m_opt_method, m_threads);
            PerformOptimisation(p + ".input_mdopt", opt_hi);

            // Read back opt_method-optimized structures for energy reporting
            std::vector<Molecule*> opt_init_stack;
            FileIterator opt_hi_in(outputPath(p + ".input_mdopt.opt.xyz"));
            while (!opt_hi_in.AtEnd()) {
                Molecule mol = opt_hi_in.Next();
                mol.setCharge(m_charge);
                mol.setSpin(m_spin);
                if (mol.AtomCount() > 0)
                    opt_init_stack.push_back(new Molecule(mol));
            }
            // Report both initial energies
            if (!m_in_stack.empty() && !opt_init_stack.empty()) {
                double md_e = m_in_stack[0]->Energy();
                double opt_e = opt_init_stack[0]->Energy();
                CurcumaLogger::result_fmt("ConfSearch: Initial energies: {}={:.6f} Eh, {}={:.6f} Eh (delta={:.2f} kJ/mol)",
                    m_md_method, md_e, m_opt_method, opt_e, (md_e - opt_e) * 2625.5);
                m_initial_energy_opt = opt_e;
            }
            for (auto* mol : opt_init_stack) delete mol;
        }
    }

    // Create shared bias pool for parallel ConfSearch.
    // When rmsd_mtd is enabled, workers share bias structures for better exploration.
    bool use_shared_pool = md.value("rmsd_mtd", false);
    if (use_shared_pool) {
        m_bias_pool = new SharedBiasPool();
        if (resumed && !m_restart.bias.empty()) {
            m_bias_pool->restoreStructures(m_restart.bias); // preserves index/counter exactly
            CurcumaLogger::result_fmt("ConfSearch: restored {} bias structures into the shared pool", m_bias_pool->biasStructureCount());
        }
        CurcumaLogger::success("Shared bias pool: Active (cross-worker bias sharing enabled)");
    }

    // Cumulative output: all accepted conformers across all temperature cycles.
    // The structures are already optimised, so the file is named ".cumulative.opt.xyz"
    // to match PerformFilter's "<f>.opt.xyz" convention for the final ConfScan below.
    m_cumulative_file = outputPath(p + ".cumulative.opt.xyz");
    const std::string& cumulative_file = m_cumulative_file; // alias for the existing in-loop appends
    if (resumed)
        writeMolVectorToFile(m_restart.cumulative, m_cumulative_file); // rebuild in the new BMT dir
    else
        std::ofstream(m_cumulative_file).close();

    // Energy reference: cycle 1's best sets the baseline; best_energy tracks the running minimum.
    // Bound to members so writeCheckpoint() can persist them (restored above on resume).
    double& initial_energy = m_initial_energy;
    double& best_energy = m_best_energy;
    // Claude Generated (Jun 2026): dual-method opt_method energy tracking (local, not checkpointed)
    double best_energy_opt = m_initial_energy_opt;

    // Claude Generated (Jun 2026): baseline RATTLE setting (whatever the user/registry chose).
    // Hot cycles override it with rattle_hot_mode; cooler cycles restore this baseline.
    nlohmann::json rattle_base = md.contains("rattle") ? md["rattle"] : nlohmann::json(0);

    // Save a copy of the initial optimised input structures as fallback seeds.
    // When a temperature cycle leaves m_in_stack empty (all structures topo/energy-rejected),
    // subsequent cycles would run 0 MD steps and re-process the same stale bias pool forever.
    // The fallback restores the initial seeds so lower-T cycles still get a fresh start.
    std::vector<Molecule*> initial_seeds;
    for (auto* mol : m_in_stack)
        initial_seeds.push_back(new Molecule(*mol));

    // On resume, continue the temperature schedule from the checkpoint; otherwise start at startT.
    int temperature_cycle = resumed ? m_restart.temperature_cycle : 0;
    const double loop_start_T = resumed ? m_restart.next_T : m_startT;
    bool stop_requested = false; // set when a 'stop' file is seen at a checkpoint boundary
    for (m_currentT = loop_start_T; m_currentT >= m_endT; m_currentT -= m_deltaT) {
        temperature_cycle++;
        // Claude Generated (Jun 2026): per-cycle wall-clock timing
        RunTimer cycle_timer;
        const int entry = pending_entry; // 0 = run MD; 1 = MD already done (resume), skip it
        pending_entry = 0;               // only the first resumed cycle re-enters mid-way
        // Verbosity is scoped by the CurcumaMethod base (ctor captures, dtor restores), so each
        // sub-phase below (Opt/MD/ConfScan) restores this level on its own — no re-assert needed.
        CurcumaLogger::header("=== ConfSearch Temperature Cycle " + std::to_string(temperature_cycle)
            + " / " + std::to_string(static_cast<int>((m_startT - m_endT) / m_deltaT) + 1)
            + " : T = " + std::to_string(static_cast<int>(m_currentT)) + " K ===");
        CurcumaLogger::result_fmt("ConfSearch: T={}K -- {} independent MD runs per structure, {} input structures, {} total runs",
            m_currentT, m_repeat, m_in_stack.size(), m_repeat * m_in_stack.size());
        CurcumaLogger::info("Each repetition starts from the same input geometry with fresh velocities (exploration via shared bias pool).");
        md["T"] = m_currentT;
        md["temperature"] = m_currentT;
        // Claude Generated (Jun 2026): auto-enable RATTLE for hot cycles. A 1 fs step at high T
        // under-samples X-H stretches (period ~10 fs) -> energy drift / spurious bond breaking;
        // constraining them (mode 2 = H-only) stabilises the dynamics. Cooler cycles keep the
        // user's/registry baseline so flexible low-T sampling is unaffected.
        if (m_currentT >= m_rattle_threshold_temp) {
            md["rattle"] = m_rattle_hot_mode;
            CurcumaLogger::result_fmt("ConfSearch: T={}K >= {}K -> RATTLE auto-enabled (mode {})",
                m_currentT, m_rattle_threshold_temp, m_rattle_hot_mode);
        } else {
            md["rattle"] = rattle_base;
        }
        // Impuls: optional initial velocity boost (reinitializes velocities each
        // step while m_impuls > m_T). ConfSearch does NOT set impuls by default
        // because SimpleMD already initializes velocities at m_T0 during Initialise().
        // Setting impuls = m_currentT would cause re-initialization EVERY step,
        // destroying dynamics. User can explicitly enable with -impuls <value>.
        if (m_defaults.contains("impuls") || m_controller.contains("impuls")) {
            md["impuls"] = m_defaults.value("impuls", m_controller.value("impuls", 0));
        }

        // Cross-temperature: log pool statistics before MD phase
        const bool in_stack_empty_before_md = m_in_stack.empty();
        const std::size_t bias_pool_size_before_md = m_bias_pool ? m_bias_pool->biasStructureCount() : 0;
        bool no_new_bias_structures = false;
        if (entry <= 0) {
            if (m_bias_pool) {
                CurcumaLogger::result_fmt("ConfSearch: Bias pool has {} structures before T={}K cycle",
                    bias_pool_size_before_md, m_currentT);
            }

#ifdef CURCUMA_DEBUG
            if (CurcumaLogger::get_verbosity() >= 3) {
                CurcumaLogger::info("MD Parameters:");
                CurcumaLogger::param_table(md, "Molecular Dynamics Settings");
            }
#endif
            // Claude Generated (Jun 2026): expose the accumulated symmetry permutations to the MTD
            // bias for this cycle (empty until the first ConfScan -> identity-only, unchanged).
            if (m_bias_pool && m_mtd_permutation && !m_permutation_cache.empty()) {
                m_bias_pool->setPermutations(m_permutation_cache);
                CurcumaLogger::result_fmt("ConfSearch: {} symmetry permutation(s) active in RMSD-MTD bias (smooth sum-over-images)",
                    static_cast<int>(m_permutation_cache.size()));
            }
            PerformMolecularDynamics(m_in_stack, md);

            // Cross-temperature: log pool statistics after MD phase and prune
            if (m_bias_pool) {
                const std::size_t post_md = m_bias_pool->biasStructureCount();
                // Prune structures with very low counter (rarely visited regions)
                // Keep at least 2 structures to maintain bias coverage
                if (post_md > 2) {
                    m_bias_pool->pruneByCounter(1);
                }
                const std::size_t post_prune = m_bias_pool->biasStructureCount();
                // Claude Generated (Jun 2026): compact bias pool delta instead of separate before/after lines
                CurcumaLogger::result_fmt("ConfSearch: Bias pool: {} -> {} after MD, {} after prune (+{} new deposits)",
                    bias_pool_size_before_md, post_md, post_prune, post_prune > bias_pool_size_before_md ? post_prune - bias_pool_size_before_md : 0);
            }
            const std::size_t bias_pool_size_after_md = m_bias_pool ? m_bias_pool->biasStructureCount() : 0;
            no_new_bias_structures = in_stack_empty_before_md
                && (bias_pool_size_after_md == bias_pool_size_before_md);

            // Within-cycle checkpoint (Claude Generated, Jun 2026): the MD exploration + grown bias
            // pool are now persisted, so an interrupt during the optimisation phases below does not
            // lose them. next_T = current T, so a resume re-enters THIS cycle at the post-MD phase.
            writeCheckpoint("post_md", m_currentT, temperature_cycle - 1);
            // Claude Generated (Jun 2026): MD phase timing
            CurcumaLogger::result_fmt("ConfSearch: MD phase took {:.1f} s", cycle_timer.Elapsed() / 1000.0);
            if (CheckStop()) {
                CurcumaLogger::warn("ConfSearch: 'stop' file detected after MD -- checkpoint written, halting.");
                stop_requested = true;
                break;
            }
        } else {
            // RESUME (post_md): the MD for this cycle already ran in the interrupted session and the
            // grown bias pool was restored above. Re-export its raw MD snapshots so Phase 2 finds its
            // input file in the new BMT dir, then fall through to the optimisation phases.
            CurcumaLogger::result_fmt("ConfSearch: resume -- skipping MD for T={}K (bias pool restored: {} structures)",
                m_currentT, m_bias_pool ? m_bias_pool->biasStructureCount() : 0);
            if (m_bias_pool && !m_in_stack.empty()) {
                auto snapshot = m_bias_pool->snapshot();
                const Molecule& ref_mol = *m_in_stack[0];
                const std::string bias_path = outputPath(p + ".bias.xyz");
                bool first = true;
                for (const auto& bs : snapshot) {
                    if (bs.persistent)
                        continue; // only raw MD snapshots are optimised in Phase 2 (matches PerformMolecularDynamics)
                    Molecule mol(ref_mol);
                    mol.setGeometry(bs.geometry);
                    mol.setName("bias_" + std::to_string(bs.index));
                    if (first) { mol.writeXYZFile(bias_path); first = false; }
                    else          mol.appendXYZFile(bias_path);
                }
                if (first)
                    std::ofstream(bias_path).close(); // no raw snapshots -> empty input
            }
        }

        CurcumaLogger::result("ConfSearch: === Phase 2: Geometry Optimisation of Bias Structures ===");
        // Skip Phase 2 when no new MD ran (empty in_stack going in) and the bias pool did not grow.
        if (no_new_bias_structures) {
            CurcumaLogger::warn_fmt("ConfSearch: T={}K -- no new MD runs and bias pool unchanged -- skipping Phase 2/3.",
                m_currentT);
        } else {
            // Fast per-cycle optimization at md_method.
            // Single-threaded per optimization when ConfSearch parallelizes externally.
            nlohmann::json opt = ChildConfig(m_md_method, (m_threads > 1) ? 1 : m_threads);
            // Bias structures are the primary conformers discovered by RMSD-MTD.
            PerformOptimisation(p + ".bias", opt);
            int opt_count = 0;
            {
                FileIterator opt_file(outputPath(p + ".bias.opt.xyz"));
                while (!opt_file.AtEnd()) { opt_file.Next(); opt_count++; }
            }
            CurcumaLogger::result_fmt("ConfSearch: Optimisation complete. {} bias structures optimised.", opt_count);
        }
        // Claude Generated (Jun 2026): Phase 2 timing
        CurcumaLogger::result_fmt("ConfSearch: Opt phase took {:.1f} s", cycle_timer.Elapsed() / 1000.0);

        CurcumaLogger::result("ConfSearch: === Phase 3: RMSD-Based Conformer Filtering ===");
        int rmsd_count = 0;
        if (no_new_bias_structures) {
            CurcumaLogger::warn_fmt("ConfSearch: T={}K -- skipping Phase 3 (no new structures).", m_currentT);
        } else {
            // "filter between": dedup at md level before the accurate re-opt.
            // Single-threaded per ConfScan when ConfSearch parallelizes externally.
            nlohmann::json scan = FilterConfig(m_md_method, (m_threads > 1) ? 1 : m_threads);
            PerformFilter(p + ".bias", scan);
            {
                FileIterator rmsd_file(outputPath(p + ".bias.opt.accepted.xyz"));
                while (!rmsd_file.AtEnd()) { rmsd_file.Next(); rmsd_count++; }
            }
            CurcumaLogger::result_fmt("ConfSearch: RMSD filtering complete. {} structures accepted.", rmsd_count);
        }

        // Claude Generated (Jun 2026): experimental adaptive calibration (Phase C). Learns the MTD
        // hill width (cluster) and/or per-atom RMSF weights (weighted) from this cycle's opt+filter
        // clustering and applies them to the NEXT cycle's MD. No-op for the default off/global.
        if (m_bias_calibration == "cluster" || m_bias_scale_mode == "weighted") {
            CalibrateBias(p, md);
        }

        // Phase 3b (Claude Generated, Jun 2026): high-level re-optimization at opt_method.
        // The md-level optimize+filter above ("filter between both optimisations") has reduced
        // the per-cycle set, so the accurate method only runs on the deduplicated survivors.
        // Skipped when opt_method == md_method, so a single-method run is unchanged and pays no
        // extra optimization. The output (".bias.opt.accepted.opt.xyz") is consumed by the
        // REFINEMENT side of Phase 4 (cumulative pool + bias); the EXPLORATION side keeps using
        // the md_method file -- the two PES are never mixed (see Phase 4).
        if (!no_new_bias_structures && m_opt_method != m_md_method) {
            CurcumaLogger::result("ConfSearch: === Phase 3b: High-Level Re-Optimisation (" + m_opt_method + ") ===");
            nlohmann::json opt_hi = ChildConfig(m_opt_method, (m_threads > 1) ? 1 : m_threads);
            // PerformOptimisation reads "<f>.xyz" and writes "<f>.opt.xyz".
            PerformOptimisation(p + ".bias.opt.accepted", opt_hi);
            int hi_count = 0;
            {
                FileIterator hi_file(outputPath(p + ".bias.opt.accepted.opt.xyz"));
                while (!hi_file.AtEnd()) { hi_file.Next(); hi_count++; }
            }
            CurcumaLogger::result_fmt("ConfSearch: High-level re-optimisation complete. {} structures re-optimised at {}.",
                hi_count, m_opt_method);
        } else if (m_opt_method == m_md_method) {
            // Claude Generated (Jun 2026): explicit skip notice at result level
            CurcumaLogger::result("ConfSearch: Phase 3b skipped (single-method mode)");
        } else if (no_new_bias_structures) {
            CurcumaLogger::result("ConfSearch: Phase 3b skipped (no new bias structures this cycle)");
        }

        CurcumaLogger::result("ConfSearch: === Phase 4: Energy Window and Topology Filter ===");
        for (auto* m : m_in_stack) delete m;
        m_in_stack.clear();
        // EXPLORATION side stays strictly on the md_method (gfnff) PES: seed selection, the
        // exploration global minimum and the bias all read the md_method minima. A region
        // explored by gfnff must never be discarded because opt_method (a different PES, e.g.
        // r2scan) ranks it higher. The opt_method structures are handled in the refinement step
        // below and only feed the FINAL ranking + an extra bias geometry -- their energies are
        // never compared to md_method energies.
        const bool dual_method = (m_opt_method != m_md_method);
        const std::string md_accepted = outputPath(p + ".bias.opt.accepted.xyz");
        double lowest_energy = std::numeric_limits<double>::infinity(); // md_method (exploration)
        int accepted = 0, rejected_topo = 0, rejected_energy = 0;
        std::vector<Molecule*> candidates;
        if (!no_new_bias_structures) {
            FileIterator file(md_accepted);
            while (!file.AtEnd()) {
                Molecule* mol = new Molecule(file.Next());
                mol->setCharge(m_charge);
                mol->setSpin(m_spin);
                // Topology check: compare bond connectivity (0/1 matrix) against reference.
                // A broken or formed bond changes >=2 entries by 1.0 -> sum >> 1e-4.
                // Log the first mismatched pair to help distinguish GFN-FF artefacts from
                // genuine chemical reactions (proton transfer, ring opening, etc.).
                auto topo_cur = mol->DistanceMatrix().second;
                double topo_diff_sum = (m_topo_matrix - topo_cur).cwiseAbs().sum();
                if (topo_diff_sum > 1e-4) {
                    if (rejected_topo == 0) {
                        // Find first differing bond for diagnostic output
                        int natoms = mol->AtomCount();
                        for (int ii = 0; ii < natoms; ++ii) {
                            bool found = false;
                            for (int jj = ii + 1; jj < natoms; ++jj) {
                                if (std::abs(m_topo_matrix(ii, jj) - topo_cur(ii, jj)) > 0.5) {
                                    CurcumaLogger::warn_fmt(
                                        "ConfSearch: topo reject (first diff): atoms {}-{} ref_bond={} cur_bond={} "
                                        "(energy {:.0f} kJ/mol above ref; total bond changes: {:.0f})",
                                        ii, jj,
                                        static_cast<int>(std::round(m_topo_matrix(ii, jj))),
                                        static_cast<int>(std::round(topo_cur(ii, jj))),
                                        (mol->Energy() - initial_energy) * 2625.5,
                                        topo_diff_sum / 2.0);
                                    found = true;
                                    break;
                                }
                            }
                            if (found) break;
                        }
                    }
                    rejected_topo++;
                    delete mol;
                    continue;
                }
                candidates.push_back(mol);
                lowest_energy = std::min(lowest_energy, mol->Energy());
            }
        }
        // Update the running global minimum across all cycles (anchor for seed selection).
        if (lowest_energy < m_global_min)
            m_global_min = lowest_energy;

        // Funnel: the seed window may shrink each cycle so later cycles only refine the
        // deepest basins. "static" keeps it constant; "exp" multiplies by decay^(cycle-1).
        double eff_seed_window = m_seed_energy_window;
        if (m_seed_window_schedule == "exp")
            eff_seed_window *= std::pow(m_seed_window_decay, temperature_cycle - 1);

        // Two independent windows (Claude Generated, Jun 2026):
        //  - cumulative OUTPUT keeps the wide energy_window relative to this cycle's lowest,
        //    so the final pool stays rich for the closing dedup;
        //  - next-cycle SEEDS use the tighter seed_energy_window relative to the GLOBAL
        //    minimum, so no MD time is spent exploring from irrelevant high-energy basins.
        // seed_rank (default 1 = only most stable): when > 0, only the N lowest-energy seeds that pass
        //  the energy window are kept for the next cycle. This focuses MD time on the
        //  currently most stable structures rather than spreading it across many basins.
        // Every topology-valid optimised minimum is also fed back into the shared bias pool
        // (opt_feedback_bias) so the next MTD cycle is biased away from what we already found.

        // Collect all candidates that pass the energy window, then optionally trim by count.
        std::vector<Molecule*> window_seeds;
        std::vector<BiasStructure> feedback;
        for (auto* mol : candidates) {
            // Cumulative output only when single-method: then md_method IS the ranking level.
            // In dual mode the cumulative pool is filled with the opt_method structures below,
            // so the final ranking never mixes the two PES.
            if (!dual_method && (mol->Energy() - lowest_energy) * 2625.5 < m_energy_window)
                mol->appendXYZFile(cumulative_file);

            // md_method minimum -> bias pool (drives the gfnff MD next cycle).
            if (m_opt_feedback_bias && m_bias_pool) {
                BiasStructure bs;
                bs.geometry = mol->getGeometry();  // full-atom, Angstrom (same units as the pool)
                bs.energy = mol->Energy();          // md_method energy (metadata only, never in the force)
                bs.counter = m_opt_feedback_height;  // hill height W = k*counter
                bs.temperature = m_currentT;
                bs.persistent = true;                // never pruned: represents a real basin
                feedback.push_back(std::move(bs));
            }

            if ((mol->Energy() - m_global_min) * 2625.5 < eff_seed_window) {
                window_seeds.push_back(mol);
            } else {
                rejected_energy++;
                delete mol;
            }
        }

        // seed_rank trimming: keep only the N lowest-energy seeds.
        if (m_seed_rank > 0 && static_cast<int>(window_seeds.size()) > m_seed_rank) {
            int rejected_by_rank = static_cast<int>(window_seeds.size()) - m_seed_rank;
            std::partial_sort(window_seeds.begin(),
                               window_seeds.begin() + m_seed_rank,
                               window_seeds.end(),
                               [](const Molecule* a, const Molecule* b) {
                                   return a->Energy() < b->Energy();
                               });
            // Delete the excess seeds.
            for (int i = m_seed_rank; i < static_cast<int>(window_seeds.size()); ++i) {
                delete window_seeds[i];
            }
            window_seeds.resize(m_seed_rank);
            rejected_energy += rejected_by_rank;
            CurcumaLogger::result_fmt("ConfSearch: seed_rank={}: keeping {} lowest-energy seeds ({} rejected by rank)",
                m_seed_rank, m_seed_rank, rejected_by_rank);
        }
        for (auto* mol : window_seeds) {
            m_in_stack.push_back(mol);
            accepted++;
        }

        // REFINEMENT side (dual-method only): the opt_method (e.g. gfn2/r2scan) re-optimised
        // structures fill the cumulative pool for the FINAL ranking and add their geometry to the
        // bias pool as a second, independent minimum. Their energies live on the opt_method PES
        // and are NEVER compared to md_method energies -- the cumulative window here is relative
        // to THIS cycle's lowest opt_method energy (same PES), and the next-cycle seeds were
        // already chosen above purely from md_method energies. So a gfnff-explored basin is kept
        // even if opt_method finds it less stable.
        if (dual_method && !no_new_bias_structures) {
            std::vector<Molecule*> opt_candidates;
            double opt_lowest = std::numeric_limits<double>::infinity();
            int opt_rejected_topo = 0;
            FileIterator ofile(outputPath(p + ".bias.opt.accepted.opt.xyz"));
            while (!ofile.AtEnd()) {
                Molecule* mol = new Molecule(ofile.Next());
                mol->setCharge(m_charge);
                mol->setSpin(m_spin);
                auto topo_cur = mol->DistanceMatrix().second;
                if ((m_topo_matrix - topo_cur).cwiseAbs().sum() > 1e-4) {
                    opt_rejected_topo++;
                    delete mol;
                    continue;
                }
                opt_candidates.push_back(mol);
                opt_lowest = std::min(opt_lowest, mol->Energy());
            }
            for (auto* mol : opt_candidates) {
                if ((mol->Energy() - opt_lowest) * 2625.5 < m_energy_window)
                    mol->appendXYZFile(cumulative_file);
                if (m_opt_feedback_bias && m_bias_pool) {
                    BiasStructure bs;
                    bs.geometry = mol->getGeometry();
                    bs.energy = mol->Energy();   // opt_method energy (metadata only, never in the force)
                    bs.counter = m_opt_feedback_height;
                    bs.temperature = m_currentT;
                    bs.persistent = true;
                    feedback.push_back(std::move(bs));
                }
                delete mol;
            }
            CurcumaLogger::result_fmt("ConfSearch: opt_method ({}) refinement: {} structures -> cumulative + bias, {} rejected (topo). Energies on the {} PES (not compared to {}).",
                m_opt_method, static_cast<int>(opt_candidates.size()), opt_rejected_topo, m_opt_method, m_md_method);
            if (opt_lowest < std::numeric_limits<double>::infinity()) {
                CurcumaLogger::result_fmt("ConfSearch: cycle lowest {} energy = {:.6f} Eh", m_opt_method, opt_lowest);
                // Claude Generated (Jun 2026): track opt_method best across cycles
                best_energy_opt = std::min(best_energy_opt, opt_lowest);
            }
        }

        if (!feedback.empty()) {
            m_bias_pool->depositBatch(feedback);
            if (m_opt_feedback_prune_snapshots) {
                // Default ON: remove raw MD snapshots now that their basins are
                // represented by the optimised persistent minima. Unoptimised
                // snapshot geometries would otherwise re-enter the bias next cycle
                // and cause marginal re-optimisation artifacts.
                // Disable with -opt_feedback_prune_snapshots false to keep all snapshots.
                m_bias_pool->pruneNonPersistent();
                CurcumaLogger::result_fmt("ConfSearch: {} optimised minima fed back, raw snapshots removed -- pool now {} structures",
                    static_cast<int>(feedback.size()), m_bias_pool->biasStructureCount());
            } else {
                CurcumaLogger::result_fmt("ConfSearch: {} optimised minima fed back (snapshots kept) -- pool now {} structures",
                    static_cast<int>(feedback.size()), m_bias_pool->biasStructureCount());
            }
        }
        if (m_seed_window_schedule == "exp")
            CurcumaLogger::result_fmt("ConfSearch: seed window (funnel) = {:.1f} kJ/mol vs. global min {:.6f} Eh",
                eff_seed_window, m_global_min);

        // Energy tracking: cycle 1 sets the initial reference; subsequent cycles compare against both.
        // This is the EXPLORATION (md_method) energy progression -- it narrates the gfnff search and
        // is intentionally NOT the opt_method ranking (logged separately above / in the final stats).
        if (temperature_cycle == 1) {
            initial_energy = lowest_energy;
            best_energy = lowest_energy;
            CurcumaLogger::result_fmt("ConfSearch: Initial best ({}) energy: {:.6f} Eh", m_md_method, initial_energy);
            // Claude Generated (Jun 2026): report initial opt_method energy in dual mode
            if (m_opt_method != m_md_method && m_initial_energy_opt < std::numeric_limits<double>::infinity())
                CurcumaLogger::result_fmt("ConfSearch: Initial best ({}) energy: {:.6f} Eh", m_opt_method, m_initial_energy_opt);
        } else if (lowest_energy < std::numeric_limits<double>::infinity()) {
            double delta_best = (best_energy - lowest_energy) * 2625.5;    // >0 = improvement vs. last best
            double delta_initial = (initial_energy - lowest_energy) * 2625.5; // >0 = improvement vs. start
            if (lowest_energy < best_energy) {
                CurcumaLogger::success_fmt("ConfSearch: New best ({})! {:.6f} Eh (+{:.2f} kJ/mol vs. prev best {:.6f} Eh, +{:.2f} kJ/mol vs. initial {:.6f} Eh)",
                    m_md_method, lowest_energy, delta_best, best_energy, delta_initial, initial_energy);
                best_energy = lowest_energy;
            } else {
                CurcumaLogger::result_fmt("ConfSearch: No new best this cycle ({}): lowest {:.6f} Eh (best still {:.6f} Eh, {:.2f} kJ/mol vs. initial {:.6f} Eh)",
                    m_md_method, lowest_energy, best_energy, delta_initial, initial_energy);
            }
            // Claude Generated (Jun 2026): report opt_method best alongside md_method in dual mode
            if (m_opt_method != m_md_method && best_energy_opt < std::numeric_limits<double>::infinity()) {
                double opt_delta_initial = (m_initial_energy_opt - best_energy_opt) * 2625.5;
                CurcumaLogger::result_fmt("ConfSearch: {} best: {:.6f} Eh ({:+.2f} kJ/mol vs. initial {:.6f} Eh)",
                    m_opt_method, best_energy_opt, opt_delta_initial, m_initial_energy_opt);
            }
        }

        CurcumaLogger::result_fmt("ConfSearch: T={}K cycle complete -- {} accepted, {} rejected (topo), {} rejected (energy), {} in next cycle",
            m_currentT, accepted, rejected_topo, rejected_energy, static_cast<int>(m_in_stack.size()));
        // Claude Generated (Jun 2026): report cumulative conformer count so the user can track progress
        {
            int cumulative_count = 0;
            FileIterator cf(cumulative_file);
            while (!cf.AtEnd()) { cf.Next(); cumulative_count++; }
            CurcumaLogger::result_fmt("ConfSearch: Cumulative conformers: {} (after T={}K cycle)", cumulative_count, static_cast<int>(m_currentT));
        }
        // Claude Generated (Jun 2026): cycle timing summary
        CurcumaLogger::result_fmt("ConfSearch: T={}K cycle took {:.1f} s", static_cast<int>(m_currentT), cycle_timer.Elapsed() / 1000.0);

        // Fallback: if Phase 4 left m_in_stack empty, restore the initial optimised seeds.
        // Without this, all remaining temperature cycles run 0 MD steps and repeatedly
        // re-optimise the same stale bias pool, wasting time and producing the same rejections.
        if (m_in_stack.empty() && !initial_seeds.empty()) {
            CurcumaLogger::warn_fmt(
                "ConfSearch: T={}K produced no valid seeds -- falling back to initial {} input structure(s) for next cycle.",
                m_currentT, static_cast<int>(initial_seeds.size()));
            for (auto* mol : initial_seeds)
                m_in_stack.push_back(new Molecule(*mol));
        }

        // End-of-cycle checkpoint (Claude Generated, Jun 2026): the cycle is complete -- cumulative
        // pool, seeds, energies and the bias pool are all final for this T. next_T points at the
        // next (lower) temperature, so a resume starts the following cycle fresh.
        writeCheckpoint("post_cycle", m_currentT - m_deltaT, temperature_cycle);
        if (CheckStop()) {
            CurcumaLogger::warn("ConfSearch: 'stop' file detected after cycle -- checkpoint written, halting.");
            stop_requested = true;
            break;
        }

        CurcumaLogger::header("=== End Temperature Cycle T = " + std::to_string(static_cast<int>(m_currentT)) + " K ===");
    }  // end temperature loop

    for (auto* mol : initial_seeds) delete mol;
    initial_seeds.clear();

    // Graceful stop: a 'stop' file was seen at a checkpoint boundary. The full state is in the
    // restart file (start dir + BMT dir); resume with -restart. Skip the final dedup (the search
    // is incomplete) and return cleanly.
    if (stop_requested) {
        CurcumaLogger::success_fmt("ConfSearch: halted by 'stop' file -- resume with: curcuma -confsearch <input> -restart  (state in {})",
            restartFileName());
        for (auto* mol : m_in_stack) delete mol;
        m_in_stack.clear();
        return;
    }

    // Final deduplication pass over all conformers collected across all temperature cycles.
    CurcumaLogger::header("=== ConfSearch: Final Deduplication Pass ===");
    int total_cumulative = 0;
    {
        FileIterator cf(cumulative_file);
        while (!cf.AtEnd()) { cf.Next(); total_cumulative++; }
        CurcumaLogger::result_fmt("ConfSearch: {} structures in cumulative pool before final filter", total_cumulative);
    }
    if (total_cumulative == 0) {
        CurcumaLogger::warn("ConfSearch: cumulative pool is empty -- all structures were rejected by topology or energy filters.");
        CurcumaLogger::warn("ConfSearch: Check the topo-reject diagnostics above. Common causes: reactive conditions (high T + strong bias),");
        CurcumaLogger::warn("ConfSearch: GFN-FF bond-length sensitivity near covalent-radii cutoffs, or mismatched input topology.");
        // Write an empty accepted file so downstream tools don't crash on a missing path.
        std::ofstream(outputPath(p + ".cumulative.opt.accepted.xyz")).close();
        CurcumaLogger::warn_fmt("ConfSearch: Empty result written to {}.cumulative.opt.accepted.xyz", p);
    } else {
        // Final ranking at the accurate level (opt_method).
        nlohmann::json final_scan = FilterConfig(m_opt_method, m_threads);
        PerformFilter(p + ".cumulative", final_scan);
        CurcumaLogger::success_fmt("ConfSearch: Final result in {}.cumulative.opt.accepted.xyz", p);

        // Claude Generated (Jun 2026): Final energy statistics over the deduplicated conformer set.
        // Read the accepted conformers back, collect their energies, and report the spread relative
        // to the lowest-energy conformer (= the deepest minimum found across all temperature cycles).
        std::vector<double> energies;
        {
            FileIterator af(outputPath(p + ".cumulative.opt.accepted.xyz"));
            while (!af.AtEnd()) {
                Molecule m = af.Next();
                energies.push_back(m.Energy());
            }
        }
        if (!energies.empty()) {
            std::sort(energies.begin(), energies.end());
            const double e_min = energies.front();
            const double e_max = energies.back();
            const double span_kj = (e_max - e_min) * 2625.5;
            CurcumaLogger::header("=== ConfSearch: Final Energy Statistics ===");
            CurcumaLogger::result_fmt("ConfSearch: {} unique conformer(s); global minimum {:.6f} Eh",
                static_cast<int>(energies.size()), e_min);
            CurcumaLogger::result_fmt("ConfSearch: energy span {:.2f} kJ/mol (lowest {:.6f} Eh, highest kept {:.6f} Eh)",
                span_kj, e_min, e_max);
            if (initial_energy < std::numeric_limits<double>::infinity()) {
                const double gain_kj = (initial_energy - e_min) * 2625.5;
                if (gain_kj > 1e-3)
                    CurcumaLogger::success_fmt("ConfSearch: search lowered the energy by {:.2f} kJ/mol vs. the initial structure ({}: {:.6f} -> {:.6f} Eh)",
                        gain_kj, m_md_method, initial_energy, e_min);
                else
                    CurcumaLogger::result_fmt("ConfSearch: initial structure remains the global minimum ({:.6f} Eh)", e_min);
            }
            // Claude Generated (Jun 2026): dual-method -- also report opt_method energy gain
            if (m_opt_method != m_md_method && m_initial_energy_opt < std::numeric_limits<double>::infinity()) {
                const double opt_gain_kj = (m_initial_energy_opt - e_min) * 2625.5;
                if (opt_gain_kj > 1e-3)
                    CurcumaLogger::success_fmt("ConfSearch: search lowered the energy by {:.2f} kJ/mol vs. the initial structure ({}: {:.6f} -> {:.6f} Eh)",
                        opt_gain_kj, m_opt_method, m_initial_energy_opt, e_min);
                else
                    CurcumaLogger::result_fmt("ConfSearch: initial structure remains the global minimum ({}: {:.6f} Eh)",
                        m_opt_method, e_min);
            }
            // Relative energies of the lowest few conformers, for a quick conformer-landscape readout.
            const int n_show = std::min(static_cast<int>(energies.size()), 10);
            for (int i = 0; i < n_show; ++i)
                CurcumaLogger::result_fmt("ConfSearch:   conformer {:>3}: {:.6f} Eh  (+{:.2f} kJ/mol)",
                    i + 1, energies[i], (energies[i] - e_min) * 2625.5);
            if (static_cast<int>(energies.size()) > n_show)
                CurcumaLogger::result_fmt("ConfSearch:   ... and {} more within the energy window",
                    static_cast<int>(energies.size()) - n_show);
        }
    }

    // Claude Generated (Apr 2026): Clean up shared bias pool
    delete m_bias_pool;
    m_bias_pool = nullptr;
}

void ConfSearch::PerformMolecularDynamics(const std::vector<Molecule*>& molecules, const nlohmann::json& parameter)
{
    CxxThreadPool* pool = new CxxThreadPool;
    int index = 0;
    CurcumaLogger::result_fmt("ConfSearch MD: Starting {} independent runs ({} repeats x {} structures), {} threads in parallel",
        m_repeat * static_cast<int>(molecules.size()), m_repeat, static_cast<int>(molecules.size()), m_threads);
    for (int repeat = 0; repeat < m_repeat; ++repeat) {
        for (size_t i = 0; i < molecules.size(); ++i) {
            MDThread* thread = new MDThread(parameter);
            thread->setThreadId(index++);
            thread->setBasename(outputPath(Basename() + ".r" + std::to_string(repeat)));
            thread->setMolecule(molecules[i]);
            thread->setSharedBiasPool(m_bias_pool);
            pool->addThread(thread);
        }
    }
    pool->setActiveThreadCount(m_threads);
    pool->StartAndWait();

    CurcumaLogger::result_fmt("ConfSearch: {} MD runs finished. Bias pool: {} structures.",
        index, m_bias_pool ? m_bias_pool->biasStructureCount() : 0);

    // Claude Generated (Jul 2026): removed the empty "confsearch.unique.xyz" stub that was
    // created here. It hardcoded the basename instead of Basename(), was always written empty
    // and was never read back -- the bias pool is the only conformer source (see Phase 2).

    // Export bias pool structures to confsearch.bias.xyz (primary conformer source).
    // Only raw MD snapshots (persistent=false) are exported for optimization — persistent
    // structures are already-optimized fed-back minima and must not be re-optimized every
    // cycle (would cause false "New best!" triggers via numerical noise).
    // confsearch.mtd.xyz gets the full unsampled pool (including persistent) for inspection.
    if (m_bias_pool && m_bias_pool->biasStructureCount() > 0 && !m_in_stack.empty()) {
        auto snapshot = m_bias_pool->snapshot();
        const Molecule& ref_mol = *m_in_stack[0];

        // .mtd.xyz: full pool for inspection
        bool first = true;
        for (const auto& bs : snapshot) {
            Molecule mol(ref_mol);
            mol.setGeometry(bs.geometry);
            mol.setName("bias_" + std::to_string(bs.index) + " t=" + std::to_string(static_cast<int>(bs.time)));
            if (first) { mol.writeXYZFile(outputPath(Basename() + ".mtd.xyz")); first = false; }
            else          mol.appendXYZFile(outputPath(Basename() + ".mtd.xyz"));
        }

        // .bias.xyz: only new MD snapshots, not already-optimized persistent minima
        std::vector<BiasStructure> new_snapshots;
        std::copy_if(snapshot.begin(), snapshot.end(), std::back_inserter(new_snapshots),
            [](const BiasStructure& bs) { return !bs.persistent; });

        int bias_count = static_cast<int>(new_snapshots.size());
        int stride = (bias_count <= m_max_bias_export) ? 1
                     : static_cast<int>(std::ceil(static_cast<double>(bias_count) / m_max_bias_export));
        int exported = 0;
        first = true;
        for (int i = 0; i < bias_count; i += stride) {
            Molecule mol(ref_mol);
            mol.setGeometry(new_snapshots[i].geometry);
            mol.setName("bias_" + std::to_string(new_snapshots[i].index));
            if (first) { mol.writeXYZFile(outputPath(Basename() + ".bias.xyz")); first = false; }
            else          mol.appendXYZFile(outputPath(Basename() + ".bias.xyz"));
            exported++;
        }
        CurcumaLogger::result_fmt("ConfSearch: {} new MD snapshots (of {} total pool, {} persistent skipped, stride={}) written to {}.bias.xyz",
            exported, static_cast<int>(snapshot.size()), static_cast<int>(snapshot.size()) - bias_count, stride, Basename());
    }

    delete pool;

    // Thread-pool boundary (Claude Generated, Jun 2026): the MD workers run and are destroyed on
    // CxxThreadPool worker threads, where the global CurcumaLogger verbosity (a shared static) is
    // left unreliable (e.g. clamped to 0 by an energy-eval on a worker) — the CurcumaMethod base
    // RAII only restores the level on the thread that owns the object, not this main thread. So the
    // pool-owning helper restores the orchestrator's level here, before returning to start().
    CurcumaLogger::set_verbosity(m_verbosity);
}

nlohmann::json ConfSearch::ChildConfig(const std::string& method, int threads) const
{
    nlohmann::json cfg;
    cfg["method"] = method;
    cfg["threads"] = threads;
    cfg["charge"] = m_charge;
    cfg["spin"] = m_spin;
    cfg["verbosity"] = m_verbosity;
    if (!m_gpu.empty() && m_gpu != "none")
        cfg["gpu"] = m_gpu;

    // Method sub-scopes that EnergyCalculator re-merges before building the method. Mirrors
    // kEnergyCalcMethodScopes in src/core/energycalculator.cpp -- without this,
    // "curcuma -confsearch mol.xyz -method gfn2 -xtb.solvent water" never reached ANY child
    // calculation: the sub-scope sits at controller["xtb"] and was simply never forwarded.
    static const char* const scopes[] = { "gfnff", "eeq_solver", "xtb", "tblite", "ulysses",
        "d3", "d4", "uff", "qmdff", "eht", "orca" };
    for (const char* s : scopes) {
        if (m_controller.contains(s) && m_controller[s].is_object())
            cfg[s] = m_controller[s];
    }
    return cfg;
}

nlohmann::json ConfSearch::FilterConfig(const std::string& energy_method, int threads) const
{
    nlohmann::json scan = ParameterRegistry::getInstance().getDefaultJson("confscan");
    scan.merge_patch(ChildConfig(energy_method, threads)); // deep merge keeps the sub-scopes

    // ConfScan's "method" is the RMSD ALIGNMENT method (default "subspace"), NOT an energy
    // method -- ChildConfig just wrote the energy method into it, so take it back.
    scan["method"] = "subspace";
    scan["energy_method"] = energy_method; // the real energy channel (ConfScan::LoadControlJson)
    scan["rmsdmethod"] = "inertia";        // cross-module alias -> rmsd.method (see rmsd.h)
    scan["fewer_file"] = true;
    scan["rmsd"] = m_rmsd;                 // the user's dedup threshold, not confscan's own 0.9
    scan["max_energy"] = m_energy_window;
    // ConfScan's own "restart" default is TRUE. A nested filter must never try to resume a
    // stale scan; the old code only got this right by accident (it inherited ConfSearch's false).
    scan["restart"] = false;
    // The input was just optimised at exactly this energy_method, so recomputing every energy is
    // pure waste. Safe because the pools handed to PerformFilter are homogeneous in level of
    // theory -- see the cumulative-append sites in start().
    scan["reuse_energies"] = true;
    return scan;
}

std::string ConfSearch::PerformOptimisation(const std::string& f, const nlohmann::json& parameter)
{
    std::string basename = f;
    std::string input_file = outputPath(basename + ".xyz");
    std::string output_file = outputPath(basename + ".opt.xyz");

    // Suppress per-step output and trajectory for batch intermediate optimizations
    json local_param = parameter;
    local_param["verbosity"] = 0;
    local_param["write_trajectory"] = false;

    // Clear output file
    std::ofstream(output_file).close();

    // Load molecules from file
    std::vector<Molecule> molecules;
    FileIterator file(input_file);
    while (!file.AtEnd()) {
        Molecule mol = file.Next();
        // Claude Generated (Jul 2026): the single choke point covering ALL four optimisation
        // sites. OptimizationDispatcher never applies a config charge (its context.charge is
        // dead storage), so without this every optimisation ran the ion as neutral.
        mol.setCharge(m_charge);
        mol.setSpin(m_spin);
        if (mol.AtomCount() > 0)
            molecules.push_back(std::move(mol));
    }

    if (molecules.empty())
        return basename;

    int total = static_cast<int>(molecules.size());

    // Claude Generated (Jun 2026): include method name in optimization output
    const std::string opt_method_name = parameter.value("method", std::string("?"));

    // Write criterion: accept the final geometry whenever it has atoms,
    // regardless of convergence. For conformational search, a partially
    // optimised structure is still useful as input for the next MD cycle.
    int written = 0;
    auto write_result = [&](const Optimization::OptimizationResult& res, const Molecule& fallback, int idx) {
        double e_start = res.energy_trajectory.empty() ? 0.0 : res.energy_trajectory.front();
        double e_end   = res.final_energy;
        double dE_kjmol = (e_end - e_start) * 2625.5;  // Eh -> kJ/mol
        if (res.final_molecule.AtomCount() > 0) {
            res.final_molecule.appendXYZFile(output_file);
            CurcumaLogger::result_fmt("  Struct {:2d} [{}]: {:4d} steps, E = {:+.6f} Eh, dE = {:+.2f} kJ/mol{}",
                idx + 1, opt_method_name, res.iterations_performed, e_end, dE_kjmol,
                res.success ? "" : "  (not converged)");
            ++written;
        } else if (fallback.AtomCount() > 0) {
            fallback.appendXYZFile(output_file);
            CurcumaLogger::result_fmt("  Struct {:2d}: optimizer failed, using input geometry", idx + 1);
            ++written;
        }
    };

    // Unified path: CxxThreadPool runs the optimizations inline on the calling
    // thread when m_threads == 1 (no worker spawned) and in parallel otherwise.
    // OptThread uses autoDelete=false, so the pool never touches these pointers
    // in its destructor — we own and delete them here.
    CxxThreadPool pool;
    pool.setActiveThreadCount(m_threads);

    std::vector<OptThread*> threads;
    threads.reserve(total);
    for (int i = 0; i < total; ++i) {
        OptThread* t = new OptThread(molecules[i], local_param);
        pool.addThread(t);
        threads.push_back(t);
    }

    CurcumaLogger::result_fmt("Optimizing {} structures using {} thread(s) [{}]", total, m_threads, opt_method_name);
    pool.StartAndWait();

    int failed = 0;
    for (int i = 0; i < static_cast<int>(threads.size()); ++i) {
        if (!threads[i]->result().success) ++failed;
        write_result(threads[i]->result(), molecules[i], i);
        delete threads[i];
    }

    if (failed > 0)
        CurcumaLogger::result_fmt("Optimization batch complete [{}]: {}/{} structures written ({} failed: zero step / gradient failure)",
            opt_method_name, written, total, failed);
    else
        CurcumaLogger::result_fmt("Optimization batch complete [{}]: {}/{} structures written to {}", opt_method_name, written, total, output_file);
    // Thread-pool boundary: restore the orchestrator's verbosity (workers leave the shared-static
    // CurcumaLogger level unreliable). See the matching note in PerformMolecularDynamics.
    CurcumaLogger::set_verbosity(m_verbosity);
    return basename;
}

std::string ConfSearch::PerformFilter(const std::string& f, const nlohmann::json& parameter)
{
    ConfScan* scan = new ConfScan(parameter, false);
    // Claude Generated (Jul 2026): propagate the BMT output directory to the nested ConfScan
    // instance. Without this, ConfScan's own outputPath() resolves without the BMT prefix while
    // ConfSearch reads the accepted file back through its (correctly prefixed) outputPath(),
    // so the two paths never agree -- ConfScan writes "<basename>.accepted.xyz" into the CWD,
    // ConfSearch looks for "<bmt_dir>/<basename>.accepted.xyz" and crashes with "File not found"
    // (uncaught std::runtime_error -> terminate). Only manifests in default BMT mode (-no_bmt
    // sidesteps it, both sides resolve to the bare filename).
    scan->setOutputDir(OutputDir());
    scan->setFileName(outputPath(f + ".opt.xyz"));
    scan->start();
    // Claude Generated (Jun 2026): harvest the symmetry/atom-permutation rules ConfScan found
    // (Hungarian reorder after inertia prealignment) and accumulate the distinct, non-identity
    // ones across cycles. They are fed into the RMSD-MTD bias as extra smooth Gaussian images.
    if (m_mtd_permutation) {
        for (const auto& rule : scan->getReorderRules()) {
            if (rule.empty())
                continue;
            bool identity = true;
            for (int i = 0; i < static_cast<int>(rule.size()); ++i)
                if (rule[i] != i) { identity = false; break; }
            if (identity)
                continue;
            if (std::find(m_permutation_cache.begin(), m_permutation_cache.end(), rule) == m_permutation_cache.end())
                m_permutation_cache.push_back(rule);
        }
    }
    delete scan;
    return f;
}

// Claude Generated (Jun 2026): permutation-aware best-fit (Kabsch) RMSD in Angstrom. Centres both
// geometries by their centroid, computes the optimal rotation, and returns the minimum RMSD over
// the identity plus every cached symmetry permutation (applied to the target). Both geometries must
// share the canonical atom order. Built directly on RMSDFunctions to avoid RMSDDriver state.
double ConfSearch::PermRMSD(const Geometry& reference, const Geometry& target) const
{
    auto kabsch = [](const Geometry& ref, const Geometry& tar) -> double {
        const int n = ref.rows();
        if (n == 0 || tar.rows() != n)
            return std::numeric_limits<double>::infinity();
        Geometry r = ref, t = tar;
        Eigen::Vector3d cr = Eigen::Vector3d::Zero(), ct = Eigen::Vector3d::Zero();
        for (int i = 0; i < n; ++i) { cr += r.row(i).transpose(); ct += t.row(i).transpose(); }
        cr /= n; ct /= n;
        for (int i = 0; i < n; ++i) { r.row(i) -= cr.transpose(); t.row(i) -= ct.transpose(); }
        Eigen::Matrix3d R = RMSDFunctions::BestFitRotation(r, t, 1);
        return RMSDFunctions::getRMSD(r, RMSDFunctions::applyRotation(t, R));
    };
    double best = kabsch(reference, target);
    const int natoms = reference.rows();
    for (const auto& rule : m_permutation_cache) {
        if (static_cast<int>(rule.size()) != natoms)
            continue;
        Geometry tp(natoms, 3);
        bool ok = true;
        for (int j = 0; j < natoms; ++j) {
            int s = rule[j];
            if (s < 0 || s >= natoms) { ok = false; break; }
            tp.row(j) = target.row(s);
        }
        if (ok)
            best = std::min(best, kabsch(reference, tp));
    }
    return best;
}

void ConfSearch::CalibrateBias(const std::string& p, nlohmann::json& md)
{
    auto load = [](const std::string& file) {
        std::vector<Molecule> v;
        FileIterator f(file);
        while (!f.AtEnd()) { Molecule m = f.Next(); if (m.AtomCount() > 0) v.push_back(m); }
        return v;
    };
    std::vector<Molecule> pre = load(outputPath(p + ".bias.xyz"));           // pre-optimisation MD snapshots
    std::vector<Molecule> post = load(outputPath(p + ".bias.opt.xyz"));      // optimised (index-aligned with pre)
    std::vector<Molecule> minima = load(outputPath(p + ".bias.opt.accepted.xyz")); // distinct minima (deduped)
    if (minima.empty() || post.empty()) {
        CurcumaLogger::warn("ConfSearch: bias calibration skipped (no clusters available this cycle)");
        return;
    }
    const int natoms = minima[0].AtomCount();
    const int nmin = static_cast<int>(minima.size());
    const bool index_aligned = (pre.size() == post.size());

    // Assign each optimised structure to its nearest distinct minimum (permutation-aware RMSD AND
    // an energy tolerance -> "same minimum" requires geometric AND energetic match).
    std::vector<int> assign(post.size(), -1);
    for (size_t k = 0; k < post.size(); ++k) {
        double bestr = std::numeric_limits<double>::infinity();
        int bestj = -1;
        for (int j = 0; j < nmin; ++j) {
            if (std::abs(post[k].Energy() - minima[j].Energy()) * 2625.5 > m_bias_energy_tol)
                continue;
            double r = PermRMSD(minima[j].getGeometry(), post[k].getGeometry());
            if (r < bestr) { bestr = r; bestj = j; }
        }
        assign[k] = bestj;
    }

    // Calibration needs >=2 distinct minima: with a single accepted minimum every structure is
    // lumped into one "cluster", which makes both the basin radius and the inter-minimum distance
    // meaningless (this produced absurd alpha values in testing). Degrade gracefully -> keep the
    // bootstrap alpha / uniform weights. Claude Generated (Jun 2026).
    if (nmin < 2 || !index_aligned) {
        CurcumaLogger::warn_fmt("ConfSearch: bias calibration skipped this cycle ({} distinct minima, index_aligned={}) -- keeping current alpha/uniform weights",
            nmin, index_aligned);
        return;
    }

    // Inter-minimum separation: the closest pair of distinct minima. A hill wider than this would
    // merge two genuinely different conformers, so it is the hard upper bound on the learned width.
    double d_inter = std::numeric_limits<double>::infinity();
    for (int a = 0; a < nmin; ++a)
        for (int b = a + 1; b < nmin; ++b)
            d_inter = std::min(d_inter, PermRMSD(minima[a].getGeometry(), minima[b].getGeometry()));
    if (!std::isfinite(d_inter) || d_inter < 1e-3) {
        CurcumaLogger::warn("ConfSearch: bias calibration skipped (degenerate inter-minimum distance)");
        return;
    }

    // Intra-cluster spread (basin capture radius) and per-atom RMSF, from the PRE-optimisation
    // geometries of structures that optimised to the same minimum, Kabsch-aligned to the cluster's
    // lowest-energy representative. ROBUSTNESS: only pairs closer than d_inter count as genuinely
    // same-basin -- a larger pre-opt RMSD signals a mis-assignment (energy-degenerate but distinct)
    // and would otherwise blow up the scale, so it is dropped from both the radius and the RMSF.
    std::vector<double> sigma2(natoms, 0.0);
    std::vector<double> basin_rmsds; // for a robust median
    long pair_count = 0;
    for (int j = 0; j < nmin; ++j) {
        int rep = -1;
        double repE = std::numeric_limits<double>::infinity();
        std::vector<int> members;
        for (size_t k = 0; k < post.size(); ++k)
            if (assign[k] == j) {
                members.push_back(static_cast<int>(k));
                if (post[k].Energy() < repE) { repE = post[k].Energy(); rep = static_cast<int>(k); }
            }
        if (rep < 0 || members.size() < 2)
            continue;
        Geometry refc = pre[rep].getGeometry();
        Eigen::Vector3d cr = Eigen::Vector3d::Zero();
        for (int i = 0; i < natoms; ++i) cr += refc.row(i).transpose();
        cr /= natoms;
        for (int i = 0; i < natoms; ++i) refc.row(i) -= cr.transpose();
        for (int k : members) {
            if (k == rep) continue;
            Geometry t = pre[k].getGeometry();
            Eigen::Vector3d ct = Eigen::Vector3d::Zero();
            for (int i = 0; i < natoms; ++i) ct += t.row(i).transpose();
            ct /= natoms;
            for (int i = 0; i < natoms; ++i) t.row(i) -= ct.transpose();
            Eigen::Matrix3d R = RMSDFunctions::BestFitRotation(refc, t, 1);
            Geometry ta = RMSDFunctions::applyRotation(t, R);
            double rmsd = RMSDFunctions::getRMSD(refc, ta);
            if (!std::isfinite(rmsd) || rmsd >= d_inter)
                continue; // drop cross-basin outliers
            basin_rmsds.push_back(rmsd);
            for (int i = 0; i < natoms; ++i)
                sigma2[i] += (refc.row(i) - ta.row(i)).squaredNorm();
            pair_count++;
        }
    }

    // --- cluster mode: set the MTD hill width from the learned, bounded basin radius ---
    if (m_bias_calibration == "cluster") {
        double r_eff;
        if (!basin_rmsds.empty()) {
            std::sort(basin_rmsds.begin(), basin_rmsds.end());
            r_eff = basin_rmsds[basin_rmsds.size() / 2]; // median capture radius (robust to outliers)
        } else {
            r_eff = 0.5 * d_inter; // no multi-member basins -> half the nearest-minima distance
        }
        // Keep the hill inside its basin: lower floor + hard cap below the nearest distinct minimum.
        r_eff = std::min(std::max(r_eff, 0.05), 0.8 * d_inter);
        double alpha = std::log(2.0) / (r_eff * r_eff);
        md["rmsd_mtd_alpha"] = alpha;
        CurcumaLogger::result_fmt("ConfSearch: bias_calibration=cluster -> alpha={:.4f} A^-2 (r_basin={:.3f} A, d_inter={:.3f} A, {} minima, {} same-basin pairs)",
            alpha, r_eff, d_inter, nmin, pair_count);
    }

    // --- weighted mode: per-atom flexibility weights w_i = 1/(sigma_i^2 + floor), clamped + normalised ---
    if (m_bias_scale_mode == "weighted") {
        if (pair_count > 0) {
            const double floor2 = 1e-2 * 1e-2; // 0.01 A floor: a near-rigid atom must not dominate
            std::vector<double> w(natoms, 1.0);
            double wmean = 0.0;
            for (int i = 0; i < natoms; ++i) {
                double s2 = sigma2[i] / static_cast<double>(pair_count);
                w[i] = 1.0 / (s2 + floor2);
                wmean += w[i];
            }
            wmean /= natoms;
            // Clamp the dynamic range to 10x around the mean (cap rigid/floppy ratio ~100) so a few
            // noisy atoms cannot dominate the metric, then normalise the mean to 1.
            if (wmean > 0)
                for (int i = 0; i < natoms; ++i)
                    w[i] = std::min(std::max(w[i] / wmean, 0.1), 10.0);
            if (m_bias_pool)
                m_bias_pool->setWeights(w);
            CurcumaLogger::result_fmt("ConfSearch: bias_scale_mode=weighted -> RMSF weights from {} same-basin pairs (rigid/floppy ratio {:.1f}, clamped <=100)",
                pair_count, *std::max_element(w.begin(), w.end()) / std::max(1e-9, *std::min_element(w.begin(), w.end())));
        } else {
            CurcumaLogger::warn("ConfSearch: weighted scale found no same-basin pairs this cycle (uniform weights kept)");
        }
    }
}

// ===================================================================================
// Claude Generated (Jun 2026): ConfSearch restart / checkpoint serialization layer.
// The whole search state is stored self-contained in "<basename>.confsearch.restart.json"
// (written into the BMT dir and copied back to the start directory). All frames share one
// atomic-number list, so each structure is just a flat "x|y|z|..." geometry string + energy.
// ===================================================================================

std::string ConfSearch::restartFileName() const
{
    return Basename() + ".confsearch.restart.json";
}

nlohmann::json ConfSearch::molToJson(const Molecule& mol) const
{
    nlohmann::json j;
    j["geometry"] = Tools::Geometry2String(mol.getGeometry());
    j["energy"] = mol.Energy();
    j["name"] = mol.Name();
    return j;
}

nlohmann::json ConfSearch::molVectorToJson(const std::vector<Molecule>& mols) const
{
    nlohmann::json arr = nlohmann::json::array();
    for (const auto& m : mols)
        arr.push_back(molToJson(m));
    return arr;
}

nlohmann::json ConfSearch::molPtrVectorToJson(const std::vector<Molecule*>& mols) const
{
    nlohmann::json arr = nlohmann::json::array();
    for (const auto* m : mols)
        if (m)
            arr.push_back(molToJson(*m));
    return arr;
}

nlohmann::json ConfSearch::fileFramesToJson(const std::string& path) const
{
    nlohmann::json arr = nlohmann::json::array();
    std::ifstream test(path);
    if (!test.good())
        return arr; // missing/empty -> empty set (e.g. a phase that never produced this file)
    test.close();
    FileIterator it(path, true);
    while (!it.AtEnd()) {
        Molecule m = it.Next();
        if (m.AtomCount() > 0)
            arr.push_back(molToJson(m));
    }
    return arr;
}

Molecule ConfSearch::jsonToMol(const std::vector<int>& elements, const nlohmann::json& entry) const
{
    const int natoms = static_cast<int>(elements.size());
    Geometry g(natoms, 3);
    Tools::String2Geometry(g, entry.value("geometry", std::string("")));
    Molecule m; // build atom-by-atom so element types are correct (no coordinate duplication)
    for (int i = 0; i < natoms; ++i)
        m.addPair({ elements[i], Position(g(i, 0), g(i, 1), g(i, 2)) });
    m.setEnergy(entry.value("energy", 0.0));
    m.setName(entry.value("name", std::string("")));
    // Claude Generated (Jul 2026): the checkpoint stores geometry + energy only, so a resumed
    // run would otherwise fall back to charge 0 / spin 0 for every restored structure.
    m.setCharge(m_charge);
    m.setSpin(m_spin);
    return m;
}

std::vector<Molecule> ConfSearch::jsonToMolVector(const std::vector<int>& elements, const nlohmann::json& arr) const
{
    std::vector<Molecule> out;
    if (!arr.is_array())
        return out;
    out.reserve(arr.size());
    for (const auto& entry : arr)
        out.push_back(jsonToMol(elements, entry));
    return out;
}

void ConfSearch::writeMolVectorToFile(const std::vector<Molecule>& mols, const std::string& path) const
{
    std::ofstream(path).close(); // truncate so an empty set yields an empty file
    bool first = true;
    for (const auto& m : mols) {
        if (first) { m.writeXYZFile(path); first = false; }
        else          m.appendXYZFile(path);
    }
}

void ConfSearch::writeCheckpoint(const std::string& phase, double next_T, int temperature_cycle)
{
    if (!m_do_restart)
        return;

    nlohmann::json state;
    state["format"] = "confsearch-restart-1";
    state["phase"] = phase; // post_md | post_filter | post_refine | post_cycle
    state["next_T"] = next_T;
    state["temperature_cycle"] = temperature_cycle;
    state["md_method"] = m_md_method;
    state["opt_method"] = m_opt_method;
    state["natoms"] = static_cast<int>(m_elements.size());
    state["elements"] = m_elements;
    state["global_min"] = m_global_min;
    state["best_energy"] = m_best_energy;
    state["initial_energy"] = m_initial_energy;

    // Full bias pool: geometry + all metadata (counter/index preserved exactly).
    nlohmann::json bias = nlohmann::json::array();
    if (m_bias_pool) {
        for (const auto& bs : m_bias_pool->snapshot()) {
            nlohmann::json b;
            b["geometry"] = Tools::Geometry2String(bs.geometry);
            b["time"] = bs.time;
            b["rmsd_reference"] = bs.rmsd_reference;
            b["energy"] = bs.energy;
            b["factor"] = bs.factor;
            b["index"] = bs.index;
            b["counter"] = bs.counter;
            b["temperature"] = bs.temperature;
            b["persistent"] = bs.persistent;
            bias.push_back(std::move(b));
        }
    }
    state["bias"] = std::move(bias);

    state["seeds"] = molPtrVectorToJson(m_in_stack);
    state["cumulative"] = fileFramesToJson(m_cumulative_file);
    state["topo_ref"] = molToJson(m_topo_ref);
    state["permutations"] = m_permutation_cache;

    // Intermediate accepted sets, only needed to resume mid-cycle without redoing the opts.
    if (phase == "post_filter" || phase == "post_refine")
        state["accepted_md"] = fileFramesToJson(outputPath(Basename() + ".bias.opt.accepted.xyz"));
    if (phase == "post_refine")
        state["accepted_opt"] = fileFramesToJson(outputPath(Basename() + ".bias.opt.accepted.opt.xyz"));

    nlohmann::json root;
    root[MethodName()[0]] = state;
    const std::string dump = root.dump();

    try {
        std::ofstream bmt(outputPath(restartFileName()));
        bmt << dump << std::endl;
        bmt.close();
        std::ofstream cwd(restartFileName()); // copy back to the start directory (stable name)
        cwd << dump << std::endl;
        cwd.close();
    } catch (...) {
        CurcumaLogger::warn("ConfSearch: failed to write restart checkpoint.");
        return;
    }
    CurcumaLogger::result_fmt("ConfSearch: checkpoint (phase={}, next_T={}K, cycles_done={}, bias={}) -> {}",
        phase, static_cast<int>(next_T), temperature_cycle,
        m_bias_pool ? m_bias_pool->biasStructureCount() : 0, restartFileName());
}

bool ConfSearch::loadCheckpoint()
{
    std::ifstream f(restartFileName());
    if (!f.good())
        return false;
    nlohmann::json root;
    try {
        f >> root;
    } catch (...) {
        CurcumaLogger::warn_fmt("ConfSearch: restart file {} is not valid JSON; starting fresh.", restartFileName());
        return false;
    }
    if (!root.contains(MethodName()[0]))
        return false;
    const nlohmann::json s = root[MethodName()[0]];

    RestartState st;
    st.md_method = s.value("md_method", std::string(""));
    st.opt_method = s.value("opt_method", std::string(""));
    if (st.md_method != m_md_method || st.opt_method != m_opt_method) {
        CurcumaLogger::warn_fmt("ConfSearch: restart method mismatch (file {}/{} vs run {}/{}); starting fresh.",
            st.md_method, st.opt_method, m_md_method, m_opt_method);
        return false;
    }
    st.natoms = s.value("natoms", 0);
    st.elements = s.value("elements", std::vector<int>{});

    const std::string phase = s.value("phase", std::string("post_cycle"));
    st.entry_phase = (phase == "post_md") ? 1 : (phase == "post_filter") ? 2
        : (phase == "post_refine")                                       ? 3
                                                                         : 0; // post_cycle -> next cycle from MD
    st.next_T = s.value("next_T", m_startT);
    st.temperature_cycle = s.value("temperature_cycle", 0);
    st.global_min = s.value("global_min", std::numeric_limits<double>::infinity());
    st.best_energy = s.value("best_energy", std::numeric_limits<double>::infinity());
    st.initial_energy = s.value("initial_energy", std::numeric_limits<double>::infinity());

    if (s.contains("bias") && s["bias"].is_array()) {
        for (const auto& b : s["bias"]) {
            BiasStructure bs;
            Geometry g(st.natoms, 3);
            Tools::String2Geometry(g, b.value("geometry", std::string("")));
            bs.geometry = g;
            bs.time = b.value("time", 0.0);
            bs.rmsd_reference = b.value("rmsd_reference", 0.0);
            bs.energy = b.value("energy", 0.0);
            bs.factor = b.value("factor", 1.0);
            bs.index = b.value("index", 0);
            bs.counter = b.value("counter", 0);
            bs.temperature = b.value("temperature", 0.0);
            bs.persistent = b.value("persistent", false);
            st.bias.push_back(std::move(bs));
        }
    }
    st.seeds = jsonToMolVector(st.elements, s.value("seeds", nlohmann::json::array()));
    st.cumulative = jsonToMolVector(st.elements, s.value("cumulative", nlohmann::json::array()));
    st.accepted_md = jsonToMolVector(st.elements, s.value("accepted_md", nlohmann::json::array()));
    st.accepted_opt = jsonToMolVector(st.elements, s.value("accepted_opt", nlohmann::json::array()));
    if (s.contains("topo_ref"))
        st.topo_ref = jsonToMol(st.elements, s["topo_ref"]);
    st.permutations = s.value("permutations", std::vector<std::vector<int>>{});

    st.valid = true;
    m_restart = std::move(st);
    return true;
}

nlohmann::json ConfSearch::WriteRestartInformation()
{
    // ConfSearch manages its own self-contained checkpoint via writeCheckpoint() (BMT + start dir).
    // The base TriggerWriteRestart() path is unused here; return an empty object.
    return nlohmann::json::object();
}

bool ConfSearch::LoadRestartInformation()
{
    return loadCheckpoint();
}

void ConfSearch::ReadControlFile()
{
}

void ConfSearch::LoadControlJson()
{
    // Claude Generated (Jul 2026): registry-backed reads via ConfigManager. It resolves the
    // intra-module aliases declared in the PARAM block (dt -> time_step, velo ->
    // initial_velocity_scale, Spin -> spin, ...), which Json2KeyWord(m_defaults, ...) cannot do.
    m_method = m_config.get<std::string>("method");
    // Claude Generated (Jun 2026): dual-method exploration/refinement. Empty values
    // fall back to "method", so a single -method run is unchanged.
    m_md_method = m_config.get<std::string>("md_method");
    if (m_md_method.empty())
        m_md_method = m_method;
    m_opt_method = m_config.get<std::string>("opt_method");
    if (m_opt_method.empty())
        m_opt_method = m_method;
    m_thermostat = m_config.get<std::string>("thermostat");
    m_spin = m_config.get<int>("spin");
    m_charge = m_config.get<int>("charge");
    // Claude Generated (Jul 2026): GPU backend, forwarded to every child computation via
    // ChildConfig(). "gpu" is a global CLI parameter, so ConfigManager picks it up from the
    // global section even when it was not given as -confsearch.gpu.
    m_gpu = m_config.get<std::string>("gpu", std::string("none"));
    m_time = m_config.get<double>("time");
    m_startT = m_config.get<double>("startT");
    m_endT = m_config.get<double>("endT");
    m_deltaT = m_config.get<double>("deltaT");
    m_repeat = m_config.get<int>("repeat");
    m_rmsd = m_config.get<double>("rmsd");
    m_threads = m_config.get<int>("threads");
    m_energy_window = m_config.get<double>("energy_window");
    m_dT = m_config.get<double>("time_step");
    m_max_bias_export = m_config.get<int>("max_bias_export");
    // Claude Generated (Jun 2026): efficiency/robustness controls
    m_rattle_threshold_temp = m_config.get<double>("rattle_threshold_temp");
    m_rattle_hot_mode = m_config.get<int>("rattle_hot_mode");
    m_topo_check = m_config.get<bool>("topo_check");
    m_topo_check_interval = m_config.get<int>("topo_check_interval");
    m_seed_energy_window = m_config.get<double>("seed_energy_window");
    m_seed_rank = m_config.get<int>("seed_rank");
    m_seed_window_schedule = m_config.get<std::string>("seed_window_schedule");
    m_seed_window_decay = m_config.get<double>("seed_window_decay");
    m_epot_abort = m_config.get<bool>("epot_abort");
    m_epot_abort_window = m_config.get<double>("epot_abort_window");
    m_temp_abort = m_config.get<bool>("temp_abort");
    m_temp_abort_factor = m_config.get<double>("temp_abort_factor");
    m_temp_abort_delta = m_config.get<double>("temp_abort_delta");
    m_rmsd_mtd_max_height = m_config.get<int>("rmsd_mtd_max_height");
    m_freeze_inherited = m_config.get<bool>("rmsd_mtd_freeze_inherited");
    m_opt_feedback_bias = m_config.get<bool>("opt_feedback_bias");
    m_opt_feedback_height = m_config.get<int>("opt_feedback_height");
    m_opt_feedback_prune_snapshots = m_config.get<bool>("opt_feedback_prune_snapshots");
    m_mtd_permutation = m_config.get<bool>("mtd_permutation");
    m_bias_calibration = m_config.get<std::string>("bias_calibration");
    m_bias_couple_factor = m_config.get<double>("bias_couple_factor");
    m_bias_scale_mode = m_config.get<std::string>("bias_scale_mode");
    m_bias_energy_tol = m_config.get<double>("bias_energy_tol");
    m_do_restart = m_config.get<bool>("restart");
}
