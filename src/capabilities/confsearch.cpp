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
    // Dump all MD parameters to file for debugging
    {
        std::ofstream debug_file("confsearch_md_params.json");
        debug_file << md.dump(2) << std::endl;
        debug_file.close();
        CurcumaLogger::result("ConfSearch: Full MD parameters written to confsearch_md_params.json");
    }

    if (md.value("rmsd_mtd", false)) {
        CurcumaLogger::result("ConfSearch: RMSD-MTD Enabled");
        double k = md.value("rmsd_mtd_k", 0.1);
        double alpha = md.value("rmsd_mtd_alpha", 10.0);
        int pace = md.value("rmsd_mtd_pace", 1);
        bool wtmtd = md.value("wtmtd", false);
        CurcumaLogger::result_fmt("ConfSearch: RMSD-MTD k={} Eh, alpha={} Bohr^-2, pace={} steps", k, alpha, pace);
        if (wtmtd)
            CurcumaLogger::result_fmt("ConfSearch: RMSD-MTD Well-tempered (dT={})", md.value("rmsd_mtd_dt", 1000000.0));
        else
            CurcumaLogger::result("ConfSearch: RMSD-MTD Well-tempered Off");
    } else {
        CurcumaLogger::result("ConfSearch: RMSD-MTD Disabled");
    }

    // Claude Generated (Apr 2026): Create shared bias pool for parallel ConfSearch.
    // When rmsd_mtd is enabled, workers share bias structures for better exploration.
    bool use_shared_pool = md.value("rmsd_mtd", false);
    if (use_shared_pool) {
        m_bias_pool = new SharedBiasPool();
        CurcumaLogger::success("Shared bias pool: Active (cross-worker bias sharing enabled)");
    }

    // Cumulative output: all accepted conformers across all temperature cycles.
    // confsearch.cumulative.xyz grows each cycle; a final ConfScan at the end deduplicates.
    const std::string cumulative_file = "confsearch.cumulative.xyz";
    std::ofstream(cumulative_file).close();

    // Energy reference: cycle 1's best sets the baseline; best_energy tracks the running minimum.
    double initial_energy = std::numeric_limits<double>::infinity();
    double best_energy = std::numeric_limits<double>::infinity();

    int temperature_cycle = 0;
    for (m_currentT = m_startT; m_currentT >= m_endT; m_currentT -= m_deltaT) {
        temperature_cycle++;
        CurcumaLogger::header("=== ConfSearch Temperature Cycle " + std::to_string(temperature_cycle)
            + " / " + std::to_string(static_cast<int>((m_startT - m_endT) / m_deltaT) + 1)
            + " : T = " + std::to_string(static_cast<int>(m_currentT)) + " K ===");
        CurcumaLogger::result_fmt("ConfSearch: T={}K -- {} independent MD runs per structure, {} input structures, {} total runs",
            m_currentT, m_repeat, m_in_stack.size(), m_repeat * m_in_stack.size());
        CurcumaLogger::info("Each repetition starts from the same input geometry with fresh velocities (exploration via shared bias pool).");
        md["T"] = m_currentT;
        md["temperature"] = m_currentT;
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
        PerformMolecularDynamics(m_in_stack, md);

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
        PerformOptimisation("confsearch.bias", opt);
        int opt_count = 0;
        {
            FileIterator opt_file("confsearch.bias.opt.xyz");
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
        PerformFilter("confsearch.bias", scan);
        int rmsd_count = 0;
        {
            FileIterator rmsd_file("confsearch.bias.opt.accepted.xyz");
            while (!rmsd_file.AtEnd()) { rmsd_file.Next(); rmsd_count++; }
        }
        CurcumaLogger::result_fmt("ConfSearch: RMSD filtering complete. {} structures accepted.", rmsd_count);

        CurcumaLogger::result("ConfSearch: === Phase 4: Energy Window and Topology Filter ===");
        for (auto* m : m_in_stack) delete m;
        m_in_stack.clear();
        double lowest_energy = std::numeric_limits<double>::infinity();
        int accepted = 0, rejected_topo = 0, rejected_energy = 0;
        std::vector<Molecule*> candidates;
        FileIterator file("confsearch.bias.opt.accepted.xyz");
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
        for (auto* mol : candidates) {
            if ((mol->Energy() - lowest_energy) * 2625.5 < m_energy_window) {
                m_in_stack.push_back(mol);
                accepted++;
            } else {
                rejected_energy++;
                delete mol;
            }
        }

        // Append accepted structures to the cumulative pool for the final output.
        for (auto* mol : m_in_stack)
            mol->appendXYZFile(cumulative_file);

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
        FileIterator cf(cumulative_file.c_str());
        while (!cf.AtEnd()) { cf.Next(); total_cumulative++; }
        CurcumaLogger::result_fmt("ConfSearch: {} structures in cumulative pool before final filter", total_cumulative);
    }
    nlohmann::json final_scan = ConfSearchJson;
    final_scan["rmsdmethod"] = "inertia";
    final_scan["fewerFile"] = true;
    final_scan["threads"] = m_threads;
    final_scan["energy_method"] = m_method;
    final_scan["max_energy"] = m_energy_window;
    PerformFilter("confsearch.cumulative", final_scan);
    CurcumaLogger::success("ConfSearch: Final result in confsearch.cumulative.accepted.xyz");

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
            thread->setBasename("confsearch.r" + std::to_string(repeat));
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
            if (first) { mol.writeXYZFile("confsearch.mtd.xyz"); first = false; }
            else          mol.appendXYZFile("confsearch.mtd.xyz");
        }

        int stride = (bias_count <= m_max_bias_export) ? 1
                     : static_cast<int>(std::ceil(static_cast<double>(bias_count) / m_max_bias_export));
        int exported = 0;
        first = true;
        for (int i = 0; i < bias_count; i += stride) {
            Molecule mol(ref_mol);
            mol.setGeometry(snapshot[i].geometry);
            mol.setName("bias_" + std::to_string(snapshot[i].index));
            if (first) { mol.writeXYZFile("confsearch.bias.xyz"); first = false; }
            else          mol.appendXYZFile("confsearch.bias.xyz");
            exported++;
        }
        CurcumaLogger::result_fmt("ConfSearch: {} bias structures (of {}, stride={}) written to confsearch.bias.xyz",
            exported, bias_count, stride);
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

    // Claude Generated (May 2026): Unified optimization path using CxxThreadPool.
    // Worker pool mode (default) is used instead of legacy mode to avoid
    // conflicts with GFN-FF's internal ForceFieldThread/CxxThreadPool.
    int total = static_cast<int>(molecules.size());
    if (total == 0)
        return basename;

    CxxThreadPool pool;
    pool.setActiveThreadCount(m_threads);

    std::vector<OptThread*> threads;
    for (const auto& mol : molecules) {
        OptThread* t = new OptThread(mol, local_param);
        pool.addThread(t);
        threads.push_back(t);
    }

    CurcumaLogger::result_fmt("Optimizing {} structures using {} threads", total, m_threads);
    pool.StartAndWait();
    CurcumaLogger::result("Optimization batch complete.");

    for (auto* t : threads) {
        if (t->result().success) {
            t->result().final_molecule.appendXYZFile(output_file);
        }
        delete t;
    }

    return basename;
}

std::string ConfSearch::PerformFilter(const std::string& f, const nlohmann::json& parameter)
{
    ConfScan* scan = new ConfScan(parameter, false);
    scan->setFileName(f + ".opt.xyz");
    scan->start();
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
}
