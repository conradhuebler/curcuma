/*
 * <Conformational Search based on Molecular Dynamics>
 * Copyright (C) 2022 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "external/CxxThreadPool/include/CxxThreadPool.hpp"

#include <fstream>
#include <iostream>
#include <stdio.h>

#include "confsearch.h"
using curcuma::Molecule;

ConfSearch::ConfSearch(const json& controller, bool silent)
    : CurcumaMethod(ConfSearchJson, controller, silent)
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
    nlohmann::json md = ParameterRegistry::getInstance().getDefaultJson("simplemd");

    // Forward user-specified parameters from the ConfSearch controller to the MD config.
    // Only forward parameters that the user explicitly set (present in m_controller),
    // not ConfSearch-specific defaults from ConfSearchJson (startT, endT, deltaT, etc.)
    // which would corrupt SimpleMD settings (e.g., rmrottrans:0 disabling COM removal).
    // Also forward from the global section (e.g., gpu settings for EnergyCalculator).
    if (m_controller.contains("confsearch") && m_controller["confsearch"].is_object()) {
        for (auto& [key, value] : m_controller["confsearch"].items()) {
            if (value.is_object())
                continue;
            md[key] = value;
        }
    }
    // Forward from global section (GPU, threads, etc.) — these are always user-specified
    if (m_controller.contains("global") && m_controller["global"].is_object()) {
        for (auto& [key, value] : m_controller["global"].items()) {
            if (!value.is_object() && key != "confsearch" && key != "confscan")
                md[key] = value;
        }
    }

    // GPU + Multi-Threading Safety: Deactivate GPU when threads > 1 to prevent
    // GPU contention. Multiple MD instances cannot share the GPU simultaneously.
    if (m_threads > 1 && md.contains("gpu") && !md["gpu"].is_null() && md["gpu"] != "none") {
        CurcumaLogger::warn("GPU cannot be used with multiple threads simultaneously. Disabling GPU for this run.");
        md["gpu"] = "none";
    }

    md["unique"] = true;
    md["method"] = m_method;
    md["rmsd"] = m_rmsd;
    // When ConfSearch itself runs multiple cycles in parallel (m_threads > 1),
    // each individual MD simulation should run single-threaded to avoid nested
    // thread pools (CxxThreadPool inside CxxThreadPool), which causes crashes.
    md["threads"] = (m_threads > 1) ? 1 : m_threads;
    md["time_step"] = m_dT;
    md["max_time"] = m_time;
    md["restart"] = false;       // ConfSearch manages its own state, no MD restart
    md["norestart"] = true;

    // Claude Generated (Jun 2026): forward the robustness controls into every MD run.
    // Both checks self-reference inside SimpleMD (start fragment count / start energy),
    // so only the enable flags and windows need to be passed; no per-seed data required.
    md["topo_check"] = m_topo_check;
    md["topo_check_interval"] = m_topo_check_interval;
    md["epot_abort"] = m_epot_abort;
    md["epot_abort_window"] = m_epot_abort_window;

    // RMSD metadynamics is the default driver for conformational exploration.
    // The SimpleMD default is false, but ConfSearch enables it by default.
    // Only disable if the user explicitly passed -rmsd_mtd false.
    if (!m_defaults.contains("rmsd_mtd") && !(m_controller.contains("rmsd_mtd")))
        md["rmsd_mtd"] = true;

    // Log ConfSearch configuration (visible at verbosity >= 1)
    CurcumaLogger::result_fmt("ConfSearch: Method={}, Thermostat={}, Threads={}", m_method, m_thermostat, m_threads);
    CurcumaLogger::result_fmt("ConfSearch: Temperature={}K -> {}K, step={}K", m_startT, m_endT, m_deltaT);
    CurcumaLogger::result_fmt("ConfSearch: Repetitions={}, RMSD threshold={} A, Energy window={} kJ/mol", m_repeat, m_rmsd, m_energy_window);
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
    if (m_bias_calibration == "couple") {
        double r = m_bias_couple_factor * m_rmsd;
        if (r > 1e-6) {
            double alpha = std::log(2.0) / (r * r);
            md["rmsd_mtd_alpha"] = alpha;
            CurcumaLogger::result_fmt("ConfSearch: bias_calibration=couple -> RMSD-MTD alpha={:.4f} A^-2 (hill half-max at {:.3f} A)", alpha, r);
        }
    } else if (m_bias_calibration != "off") {
        CurcumaLogger::warn_fmt("ConfSearch: bias_calibration='{}' not yet implemented (cluster/weighted pending) -- treating as off", m_bias_calibration);
    }

    // Dump all MD parameters to file for debugging
    {
        std::ofstream debug_file(p + "_md_params.json");
        debug_file << md.dump(2) << std::endl;
        debug_file.close();
        CurcumaLogger::result_fmt("ConfSearch: Full MD parameters written to {}_md_params.json", p);
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

    // Optimise all input structures before any MD run.
    // m_topo_matrix is updated from the first optimised structure so Phase 4
    // topology checks compare against the relaxed geometry, not the raw input.
    CurcumaLogger::header("=== ConfSearch: Initial Geometry Optimisation ===");
    {
        bool first = true;
        for (auto* mol : m_in_stack) {
            if (first) { mol->writeXYZFile(p + ".input.xyz"); first = false; }
            else          mol->appendXYZFile(p + ".input.xyz");
        }
        nlohmann::json opt_init;
        opt_init["method"] = m_method;
        opt_init["threads"] = m_threads;
        if (md.contains("gpu") && !md["gpu"].is_null())
            opt_init["gpu"] = md["gpu"];
        PerformOptimisation(p + ".input", opt_init);
        CurcumaLogger::set_verbosity(m_verbosity);  // re-assert after initial optimisation (see temperature loop)

        for (auto* mol : m_in_stack) delete mol;
        m_in_stack.clear();
        FileIterator opt_in(p + ".input.opt.xyz");
        while (!opt_in.AtEnd()) {
            Molecule mol = opt_in.Next();
            if (mol.AtomCount() > 0)
                m_in_stack.push_back(new Molecule(mol));
        }
        if (!m_in_stack.empty())
            m_topo_matrix = m_in_stack[0]->DistanceMatrix().second;
        CurcumaLogger::result_fmt("ConfSearch: {} input structures optimised", m_in_stack.size());
    }

    // Create shared bias pool for parallel ConfSearch.
    // When rmsd_mtd is enabled, workers share bias structures for better exploration.
    bool use_shared_pool = md.value("rmsd_mtd", false);
    if (use_shared_pool) {
        m_bias_pool = new SharedBiasPool();
        CurcumaLogger::success("Shared bias pool: Active (cross-worker bias sharing enabled)");
    }

    // Cumulative output: all accepted conformers across all temperature cycles.
    // The structures are already optimised, so the file is named ".cumulative.opt.xyz"
    // to match PerformFilter's "<f>.opt.xyz" convention for the final ConfScan below.
    const std::string cumulative_file = p + ".cumulative.opt.xyz";
    std::ofstream(cumulative_file).close();

    // Energy reference: cycle 1's best sets the baseline; best_energy tracks the running minimum.
    double initial_energy = std::numeric_limits<double>::infinity();
    double best_energy = std::numeric_limits<double>::infinity();

    // Claude Generated (Jun 2026): baseline RATTLE setting (whatever the user/registry chose).
    // Hot cycles override it with rattle_hot_mode; cooler cycles restore this baseline.
    nlohmann::json rattle_base = md.contains("rattle") ? md["rattle"] : nlohmann::json(0);

    int temperature_cycle = 0;
    for (m_currentT = m_startT; m_currentT >= m_endT; m_currentT -= m_deltaT) {
        temperature_cycle++;
        // Claude Generated (Jun 2026): sub-methods (initial Opt, per-cycle Opt/ConfScan/MD)
        // drive the global CurcumaLogger verbosity down to their own level and do not restore
        // it, so ConfSearch's per-cycle progress would be invisible. Re-assert our level here
        // and after each sub-phase below so the cycle logs are actually shown.
        CurcumaLogger::set_verbosity(m_verbosity);
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
        if (m_bias_pool) {
            CurcumaLogger::result_fmt("ConfSearch: Bias pool has {} structures before T={}K cycle",
                m_bias_pool->biasStructureCount(), m_currentT);
        }

#ifdef CURCUMA_DEBUG
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info("MD Parameters:");
            CurcumaLogger::param_table(md, "Molecular Dynamics Settings");
        }
#endif
        // std::vector<Molecule*> uniques;
        // Claude Generated (Jun 2026): expose the accumulated symmetry permutations to the MTD
        // bias for this cycle (empty until the first ConfScan -> identity-only, unchanged).
        if (m_bias_pool && m_mtd_permutation && !m_permutation_cache.empty()) {
            m_bias_pool->setPermutations(m_permutation_cache);
            CurcumaLogger::result_fmt("ConfSearch: {} symmetry permutation(s) active in RMSD-MTD bias (smooth sum-over-images)",
                static_cast<int>(m_permutation_cache.size()));
        }
        PerformMolecularDynamics(m_in_stack, md);
        CurcumaLogger::set_verbosity(m_verbosity);  // re-assert after MD sub-instances (see cycle top)

        // Cross-temperature: log pool statistics after MD phase and prune
        if (m_bias_pool) {
            CurcumaLogger::result_fmt("ConfSearch: Bias pool has {} structures after T={}K MD",
                m_bias_pool->biasStructureCount(), m_currentT);
            // Prune structures with very low counter (rarely visited regions)
            // Keep at least 2 structures to maintain bias coverage
            if (m_bias_pool->biasStructureCount() > 2) {
                m_bias_pool->pruneByCounter(1);
                CurcumaLogger::result_fmt("ConfSearch: Bias pool pruned to {} structures",
                    m_bias_pool->biasStructureCount());
            }
        }

        CurcumaLogger::result("ConfSearch: === Phase 2: Geometry Optimisation of Bias Structures ===");
        nlohmann::json opt;
        opt["method"] = m_method;
        // Single-threaded per optimization when ConfSearch parallelizes externally
        opt["threads"] = (m_threads > 1) ? 1 : m_threads;
        if (md.contains("gpu") && !md["gpu"].is_null())
            opt["gpu"] = md["gpu"];
        // Bias structures are the primary conformers discovered by RMSD-MTD.
        // MD unique snapshots (confsearch.unique.xyz) are secondary and not used here.
        PerformOptimisation(p + ".bias", opt);
        CurcumaLogger::set_verbosity(m_verbosity);  // re-assert after optimisation (see cycle top)
        int opt_count = 0;
        {
            FileIterator opt_file(p + ".bias.opt.xyz");
            while (!opt_file.AtEnd()) { opt_file.Next(); opt_count++; }
        }
        CurcumaLogger::result_fmt("ConfSearch: Optimisation complete. {} bias structures optimised.", opt_count);

        CurcumaLogger::result("ConfSearch: === Phase 3: RMSD-Based Conformer Filtering ===");
        nlohmann::json scan = ConfSearchJson;
        scan["rmsdmethod"] = "inertia";
        scan["fewerFile"] = true;
        // Single-threaded per ConfScan when ConfSearch parallelizes externally
        scan["threads"] = (m_threads > 1) ? 1 : m_threads;
        scan["energy_method"] = m_method;
        scan["max_energy"] = m_energy_window;
        if (md.contains("gpu") && !md["gpu"].is_null())
            scan["gpu"] = md["gpu"];
        PerformFilter(p + ".bias", scan);
        CurcumaLogger::set_verbosity(m_verbosity);  // re-assert after ConfScan filter (see cycle top)
        int rmsd_count = 0;
        {
            FileIterator rmsd_file(p + ".bias.opt.accepted.xyz");
            while (!rmsd_file.AtEnd()) { rmsd_file.Next(); rmsd_count++; }
        }
        CurcumaLogger::result_fmt("ConfSearch: RMSD filtering complete. {} structures accepted.", rmsd_count);

        CurcumaLogger::result("ConfSearch: === Phase 4: Energy Window and Topology Filter ===");
        for (auto* m : m_in_stack) delete m;
        m_in_stack.clear();
        double lowest_energy = std::numeric_limits<double>::infinity();
        int accepted = 0, rejected_topo = 0, rejected_energy = 0;
        std::vector<Molecule*> candidates;
        FileIterator file(p + ".bias.opt.accepted.xyz");
        while (!file.AtEnd()) {
            Molecule* mol = new Molecule(file.Next());
            // Float-safe topology check: any single broken/formed bond changes distances by >> 1e-4 A
            if ((m_topo_matrix - mol->DistanceMatrix().second).cwiseAbs().sum() > 1e-4) {
                rejected_topo++;
                delete mol;
                continue;
            }
            candidates.push_back(mol);
            lowest_energy = std::min(lowest_energy, mol->Energy());
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
        // Every topology-valid optimised minimum is also fed back into the shared bias pool
        // (opt_feedback_bias) so the next MTD cycle is biased away from what we already found.
        std::vector<BiasStructure> feedback;
        for (auto* mol : candidates) {
            if ((mol->Energy() - lowest_energy) * 2625.5 < m_energy_window)
                mol->appendXYZFile(cumulative_file);

            if (m_opt_feedback_bias && m_bias_pool) {
                BiasStructure bs;
                bs.geometry = mol->getGeometry();  // full-atom, Angstrom (same units as the pool)
                bs.energy = mol->Energy();
                bs.counter = m_opt_feedback_height;  // hill height W = k*counter
                bs.temperature = m_currentT;
                bs.persistent = true;                // never pruned: represents a real basin
                feedback.push_back(std::move(bs));
            }

            if ((mol->Energy() - m_global_min) * 2625.5 < eff_seed_window) {
                m_in_stack.push_back(mol);
                accepted++;
            } else {
                rejected_energy++;
                delete mol;
            }
        }
        if (!feedback.empty()) {
            m_bias_pool->depositBatch(feedback);
            CurcumaLogger::result_fmt("ConfSearch: {} optimised minima fed back into bias pool (height={}), pool now {} structures",
                static_cast<int>(feedback.size()), m_opt_feedback_height, m_bias_pool->biasStructureCount());
        }
        if (m_seed_window_schedule == "exp")
            CurcumaLogger::result_fmt("ConfSearch: seed window (funnel) = {:.1f} kJ/mol vs. global min {:.6f} Eh",
                eff_seed_window, m_global_min);

        // Energy tracking: cycle 1 sets the initial reference; subsequent cycles compare against both.
        if (temperature_cycle == 1) {
            initial_energy = lowest_energy;
            best_energy = lowest_energy;
            CurcumaLogger::result_fmt("ConfSearch: Initial best energy: {:.6f} Eh", initial_energy);
        } else if (lowest_energy < std::numeric_limits<double>::infinity()) {
            double vs_initial = (lowest_energy - initial_energy) * 2625.5;
            if (lowest_energy < best_energy) {
                double vs_best = (lowest_energy - best_energy) * 2625.5;
                CurcumaLogger::success_fmt("ConfSearch: New best! {:.6f} Eh | vs. initial: {:.2f} kJ/mol | improvement: {:.2f} kJ/mol",
                    lowest_energy, vs_initial, vs_best);
                best_energy = lowest_energy;
            } else {
                CurcumaLogger::result_fmt("ConfSearch: Lowest this cycle: {:.6f} Eh ({:.2f} kJ/mol vs. initial, best still {:.6f} Eh)",
                    lowest_energy, vs_initial, best_energy);
            }
        }

        CurcumaLogger::result_fmt("ConfSearch: T={}K cycle complete -- {} accepted, {} rejected (topo), {} rejected (energy), {} in next cycle",
            m_currentT, accepted, rejected_topo, rejected_energy, static_cast<int>(m_in_stack.size()));
        CurcumaLogger::header("=== End Temperature Cycle T = " + std::to_string(static_cast<int>(m_currentT)) + " K ===");
    }  // end temperature loop

    // Final deduplication pass over all conformers collected across all temperature cycles.
    CurcumaLogger::header("=== ConfSearch: Final Deduplication Pass ===");
    {
        int total_cumulative = 0;
        FileIterator cf(cumulative_file);
        while (!cf.AtEnd()) { cf.Next(); total_cumulative++; }
        CurcumaLogger::result_fmt("ConfSearch: {} structures in cumulative pool before final filter", total_cumulative);
    }
    nlohmann::json final_scan = ConfSearchJson;
    final_scan["rmsdmethod"] = "inertia";
    final_scan["fewerFile"] = true;
    final_scan["threads"] = m_threads;
    final_scan["energy_method"] = m_method;
    final_scan["max_energy"] = m_energy_window;
    PerformFilter(p + ".cumulative", final_scan);
    CurcumaLogger::set_verbosity(m_verbosity);  // re-assert after final ConfScan filter (see temperature loop)
    CurcumaLogger::success_fmt("ConfSearch: Final result in {}.cumulative.opt.accepted.xyz", p);

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
            thread->setBasename(Basename() + ".r" + std::to_string(repeat));
            thread->setMolecule(molecules[i]);
            thread->setSharedBiasPool(m_bias_pool);
            pool->addThread(thread);
        }
    }
    pool->setActiveThreadCount(m_threads);
    pool->StartAndWait();

    CurcumaLogger::result_fmt("ConfSearch: {} MD runs finished. Bias pool: {} structures.",
        index, m_bias_pool ? m_bias_pool->biasStructureCount() : 0);

    // Export bias pool structures to confsearch.bias.xyz (primary conformer source).
    // confsearch.mtd.xyz gets the full unsampled pool for inspection.
    if (m_bias_pool && m_bias_pool->biasStructureCount() > 0 && !m_in_stack.empty()) {
        auto snapshot = m_bias_pool->snapshot();
        int bias_count = static_cast<int>(snapshot.size());
        const Molecule& ref_mol = *m_in_stack[0];

        bool first = true;
        for (const auto& bs : snapshot) {
            Molecule mol(ref_mol);
            mol.setGeometry(bs.geometry);
            mol.setName("bias_" + std::to_string(bs.index) + " t=" + std::to_string(static_cast<int>(bs.time)));
            if (first) { mol.writeXYZFile(Basename() + ".mtd.xyz"); first = false; }
            else          mol.appendXYZFile(Basename() + ".mtd.xyz");
        }

        int stride = (bias_count <= m_max_bias_export) ? 1
                     : static_cast<int>(std::ceil(static_cast<double>(bias_count) / m_max_bias_export));
        int exported = 0;
        first = true;
        for (int i = 0; i < bias_count; i += stride) {
            Molecule mol(ref_mol);
            mol.setGeometry(snapshot[i].geometry);
            mol.setName("bias_" + std::to_string(snapshot[i].index));
            if (first) { mol.writeXYZFile(Basename() + ".bias.xyz"); first = false; }
            else          mol.appendXYZFile(Basename() + ".bias.xyz");
            exported++;
        }
        CurcumaLogger::result_fmt("ConfSearch: {} bias structures (of {}, stride={}) written to {}.bias.xyz",
            exported, bias_count, stride, Basename());
    }

    delete pool;
}

std::string ConfSearch::PerformOptimisation(const std::string& f, const nlohmann::json& parameter)
{
    std::string basename = f;
    std::string input_file = basename + ".xyz";
    std::string output_file = basename + ".opt.xyz";

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
        if (mol.AtomCount() > 0)
            molecules.push_back(std::move(mol));
    }

    if (molecules.empty())
        return basename;

    int total = static_cast<int>(molecules.size());

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
            CurcumaLogger::result_fmt("  Struct {:2d}: {:4d} steps, E = {:+.6f} Eh, dE = {:+.2f} kJ/mol{}",
                idx + 1, res.iterations_performed, e_end, dE_kjmol,
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

    CurcumaLogger::result_fmt("Optimizing {} structures using {} thread(s)", total, m_threads);
    pool.StartAndWait();

    for (int i = 0; i < static_cast<int>(threads.size()); ++i) {
        write_result(threads[i]->result(), molecules[i], i);
        delete threads[i];
    }

    CurcumaLogger::result_fmt("Optimization batch complete: {}/{} structures written to {}", written, total, output_file);
    return basename;
}

std::string ConfSearch::PerformFilter(const std::string& f, const nlohmann::json& parameter)
{
    ConfScan* scan = new ConfScan(parameter, false);
    scan->setFileName(f + ".opt.xyz");
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

nlohmann::json ConfSearch::WriteRestartInformation()
{
    nlohmann::json restart;
    return restart;
}

bool ConfSearch::LoadRestartInformation()
{
    return true;
}

void ConfSearch::ReadControlFile()
{
}

void ConfSearch::LoadControlJson()
{
    m_method = Json2KeyWord<std::string>(m_defaults, "method");
    m_thermostat = Json2KeyWord<std::string>(m_defaults, "thermostat");
    m_rattle = Json2KeyWord<bool>(m_defaults, "rattle");
    m_spin = Json2KeyWord<int>(m_defaults, "spin");
    m_charge = Json2KeyWord<int>(m_defaults, "charge");
    //    m_single_step = Json2KeyWord<double>(m_defaults, "dT"); // * fs2amu;
    m_time = Json2KeyWord<double>(m_defaults, "time");
    m_startT = Json2KeyWord<double>(m_defaults, "startT");
    m_endT = Json2KeyWord<double>(m_defaults, "endT");
    m_deltaT = Json2KeyWord<double>(m_defaults, "deltaT");
    m_repeat = Json2KeyWord<int>(m_defaults, "repeat");
    m_rmsd = Json2KeyWord<double>(m_defaults, "rmsd");
    m_threads = Json2KeyWord<int>(m_defaults, "threads");
    m_energy_window = Json2KeyWord<double>(m_defaults, "energy_window");
    m_dT = Json2KeyWord<double>(m_defaults, "dT");
    m_max_bias_export = Json2KeyWord<int>(m_defaults, "max_bias_export");
    // Claude Generated (Jun 2026): efficiency/robustness controls
    m_rattle_threshold_temp = Json2KeyWord<double>(m_defaults, "rattle_threshold_temp");
    m_rattle_hot_mode = Json2KeyWord<int>(m_defaults, "rattle_hot_mode");
    m_topo_check = Json2KeyWord<bool>(m_defaults, "topo_check");
    m_topo_check_interval = Json2KeyWord<int>(m_defaults, "topo_check_interval");
    m_seed_energy_window = Json2KeyWord<double>(m_defaults, "seed_energy_window");
    m_seed_window_schedule = Json2KeyWord<std::string>(m_defaults, "seed_window_schedule");
    m_seed_window_decay = Json2KeyWord<double>(m_defaults, "seed_window_decay");
    m_epot_abort = Json2KeyWord<bool>(m_defaults, "epot_abort");
    m_epot_abort_window = Json2KeyWord<double>(m_defaults, "epot_abort_window");
    m_opt_feedback_bias = Json2KeyWord<bool>(m_defaults, "opt_feedback_bias");
    m_opt_feedback_height = Json2KeyWord<int>(m_defaults, "opt_feedback_height");
    m_mtd_permutation = Json2KeyWord<bool>(m_defaults, "mtd_permutation");
    m_bias_calibration = Json2KeyWord<std::string>(m_defaults, "bias_calibration");
    m_bias_couple_factor = Json2KeyWord<double>(m_defaults, "bias_couple_factor");
}
