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
    // These checks self-reference inside SimpleMD (start fragment count / start energy /
    // target T / per-run inherited pool), so only the enable flags and windows need to be
    // forwarded; no per-seed data required.
    //
    // IMPORTANT: these are SimpleMD-owned PARAMs, so a flat CLI flag (e.g. -topo_check,
    // -temp_abort, -rmsd_mtd_freeze_inherited) is routed by the auto-router to
    // controller["simplemd"], NOT controller["confsearch"]. ConfSearch does not otherwise
    // consume controller["simplemd"], so we must read the user value from there (falling back
    // to the ConfSearch default) -- otherwise the flag is silently ignored. The dotted form
    // -confsearch.<key> still works via the ConfSearch member default.
    json sd = (m_controller.contains("simplemd") && m_controller["simplemd"].is_object())
        ? m_controller["simplemd"] : json::object();
    md["topo_check"] = sd.value("topo_check", m_topo_check);
    md["topo_check_interval"] = sd.value("topo_check_interval", m_topo_check_interval);
    md["epot_abort"] = sd.value("epot_abort", m_epot_abort);
    md["epot_abort_window"] = sd.value("epot_abort_window", m_epot_abort_window);
    md["temp_abort"] = sd.value("temp_abort", m_temp_abort);
    md["temp_abort_factor"] = sd.value("temp_abort_factor", m_temp_abort_factor);
    md["temp_abort_delta"] = sd.value("temp_abort_delta", m_temp_abort_delta);
    md["rmsd_mtd_max_height"] = sd.value("rmsd_mtd_max_height", m_rmsd_mtd_max_height);
    md["rmsd_mtd_freeze_inherited"] = sd.value("rmsd_mtd_freeze_inherited", m_freeze_inherited);

    // RMSD metadynamics is the default driver for conformational exploration.
    // The SimpleMD default is false, but ConfSearch enables it by default.
    // Only disable if the user explicitly passed -rmsd_mtd false.
    if (!m_defaults.contains("rmsd_mtd") && !(m_controller.contains("rmsd_mtd")))
        md["rmsd_mtd"] = true;

    // Log ConfSearch configuration (visible at verbosity >= 1)
    CurcumaLogger::result_fmt("ConfSearch: Method={}, Thermostat={}, Threads={}", m_method, m_thermostat, m_threads);
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
            if (first) { mol->writeXYZFile(outputPath(p + ".input.xyz")); first = false; }
            else          mol->appendXYZFile(outputPath(p + ".input.xyz"));
        }
        nlohmann::json opt_init;
        opt_init["method"] = m_method;
        opt_init["threads"] = m_threads;
        if (md.contains("gpu") && !md["gpu"].is_null())
            opt_init["gpu"] = md["gpu"];
        PerformOptimisation(p + ".input", opt_init);

        for (auto* mol : m_in_stack) delete mol;
        m_in_stack.clear();
        FileIterator opt_in(outputPath(p + ".input.opt.xyz"));
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
    const std::string cumulative_file = outputPath(p + ".cumulative.opt.xyz");
    std::ofstream(cumulative_file).close();

    // Energy reference: cycle 1's best sets the baseline; best_energy tracks the running minimum.
    double initial_energy = std::numeric_limits<double>::infinity();
    double best_energy = std::numeric_limits<double>::infinity();

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

    int temperature_cycle = 0;
    for (m_currentT = m_startT; m_currentT >= m_endT; m_currentT -= m_deltaT) {
        temperature_cycle++;
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
        // std::vector<Molecule*> uniques;
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
        // Skip Phase 2 when no new MD ran (empty in_stack going in) and the bias pool
        // did not grow. Re-optimising the same stale bias structures would produce the
        // same topo/energy rejections and waste the entire remaining temperature schedule.
        const std::size_t bias_pool_size_after_md = m_bias_pool ? m_bias_pool->biasStructureCount() : 0;
        const bool no_new_bias_structures = in_stack_empty_before_md
            && (bias_pool_size_after_md == bias_pool_size_before_md);
        if (no_new_bias_structures) {
            CurcumaLogger::warn_fmt("ConfSearch: T={}K -- no new MD runs and bias pool unchanged -- skipping Phase 2/3.",
                m_currentT);
        } else {
            nlohmann::json opt;
            opt["method"] = m_method;
            // Single-threaded per optimization when ConfSearch parallelizes externally
            opt["threads"] = (m_threads > 1) ? 1 : m_threads;
            if (md.contains("gpu") && !md["gpu"].is_null())
                opt["gpu"] = md["gpu"];
            // Bias structures are the primary conformers discovered by RMSD-MTD.
            // MD unique snapshots (confsearch.unique.xyz) are secondary and not used here.
            PerformOptimisation(p + ".bias", opt);
            int opt_count = 0;
            {
                FileIterator opt_file(outputPath(p + ".bias.opt.xyz"));
                while (!opt_file.AtEnd()) { opt_file.Next(); opt_count++; }
            }
            CurcumaLogger::result_fmt("ConfSearch: Optimisation complete. {} bias structures optimised.", opt_count);
        }

        CurcumaLogger::result("ConfSearch: === Phase 3: RMSD-Based Conformer Filtering ===");
        int rmsd_count = 0;
        if (no_new_bias_structures) {
            CurcumaLogger::warn_fmt("ConfSearch: T={}K -- skipping Phase 3 (no new structures).", m_currentT);
        } else {
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

        CurcumaLogger::result("ConfSearch: === Phase 4: Energy Window and Topology Filter ===");
        for (auto* m : m_in_stack) delete m;
        m_in_stack.clear();
        double lowest_energy = std::numeric_limits<double>::infinity();
        int accepted = 0, rejected_topo = 0, rejected_energy = 0;
        std::vector<Molecule*> candidates;
        if (!no_new_bias_structures) {
            FileIterator file(outputPath(p + ".bias.opt.accepted.xyz"));
            while (!file.AtEnd()) {
                Molecule* mol = new Molecule(file.Next());
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
        if (temperature_cycle == 1) {
            initial_energy = lowest_energy;
            best_energy = lowest_energy;
            CurcumaLogger::result_fmt("ConfSearch: Initial best energy: {:.6f} Eh", initial_energy);
        } else if (lowest_energy < std::numeric_limits<double>::infinity()) {
            double delta_best = (best_energy - lowest_energy) * 2625.5;    // >0 = improvement vs. last best
            double delta_initial = (initial_energy - lowest_energy) * 2625.5; // >0 = improvement vs. start
            if (lowest_energy < best_energy) {
                CurcumaLogger::success_fmt("ConfSearch: New best! {:.6f} Eh (+{:.2f} kJ/mol vs. prev best {:.6f} Eh, +{:.2f} kJ/mol vs. initial {:.6f} Eh)",
                    lowest_energy, delta_best, best_energy, delta_initial, initial_energy);
                best_energy = lowest_energy;
            } else {
                CurcumaLogger::result_fmt("ConfSearch: No new best this cycle: lowest {:.6f} Eh (best still {:.6f} Eh, {:.2f} kJ/mol vs. initial {:.6f} Eh)",
                    lowest_energy, best_energy, delta_initial, initial_energy);
            }
        }

        CurcumaLogger::result_fmt("ConfSearch: T={}K cycle complete -- {} accepted, {} rejected (topo), {} rejected (energy), {} in next cycle",
            m_currentT, accepted, rejected_topo, rejected_energy, static_cast<int>(m_in_stack.size()));

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

        CurcumaLogger::header("=== End Temperature Cycle T = " + std::to_string(static_cast<int>(m_currentT)) + " K ===");
    }  // end temperature loop

    for (auto* mol : initial_seeds) delete mol;
    initial_seeds.clear();

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
        nlohmann::json final_scan = ConfSearchJson;
        final_scan["rmsdmethod"] = "inertia";
        final_scan["fewerFile"] = true;
        final_scan["threads"] = m_threads;
        final_scan["energy_method"] = m_method;
        final_scan["max_energy"] = m_energy_window;
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
                    CurcumaLogger::success_fmt("ConfSearch: search lowered the energy by {:.2f} kJ/mol vs. the initial structure ({:.6f} -> {:.6f} Eh)",
                        gain_kj, initial_energy, e_min);
                else
                    CurcumaLogger::result_fmt("ConfSearch: initial structure remains the global minimum ({:.6f} Eh)", e_min);
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

    std::string file = outputPath("confsearch.unique.xyz");
    std::ofstream result_file;
    result_file.open(file);
    result_file.close();

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

    int failed = 0;
    for (int i = 0; i < static_cast<int>(threads.size()); ++i) {
        if (!threads[i]->result().success) ++failed;
        write_result(threads[i]->result(), molecules[i], i);
        delete threads[i];
    }

    if (failed > 0)
        CurcumaLogger::result_fmt("Optimization batch complete: {}/{} structures written ({} failed: zero step / gradient failure)",
            written, total, failed);
    else
        CurcumaLogger::result_fmt("Optimization batch complete: {}/{} structures written to {}", written, total, output_file);
    // Thread-pool boundary: restore the orchestrator's verbosity (workers leave the shared-static
    // CurcumaLogger level unreliable). See the matching note in PerformMolecularDynamics.
    CurcumaLogger::set_verbosity(m_verbosity);
    return basename;
}

std::string ConfSearch::PerformFilter(const std::string& f, const nlohmann::json& parameter)
{
    ConfScan* scan = new ConfScan(parameter, false);
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
    m_seed_rank = Json2KeyWord<int>(m_defaults, "seed_rank");
    m_seed_window_schedule = Json2KeyWord<std::string>(m_defaults, "seed_window_schedule");
    m_seed_window_decay = Json2KeyWord<double>(m_defaults, "seed_window_decay");
    m_epot_abort = Json2KeyWord<bool>(m_defaults, "epot_abort");
    m_epot_abort_window = Json2KeyWord<double>(m_defaults, "epot_abort_window");
    m_temp_abort = Json2KeyWord<bool>(m_defaults, "temp_abort");
    m_temp_abort_factor = Json2KeyWord<double>(m_defaults, "temp_abort_factor");
    m_temp_abort_delta = Json2KeyWord<double>(m_defaults, "temp_abort_delta");
    m_rmsd_mtd_max_height = Json2KeyWord<int>(m_defaults, "rmsd_mtd_max_height");
    m_freeze_inherited = Json2KeyWord<bool>(m_defaults, "rmsd_mtd_freeze_inherited");
    m_opt_feedback_bias = Json2KeyWord<bool>(m_defaults, "opt_feedback_bias");
    m_opt_feedback_height = Json2KeyWord<int>(m_defaults, "opt_feedback_height");
    m_opt_feedback_prune_snapshots = Json2KeyWord<bool>(m_defaults, "opt_feedback_prune_snapshots");
    m_mtd_permutation = Json2KeyWord<bool>(m_defaults, "mtd_permutation");
    m_bias_calibration = Json2KeyWord<std::string>(m_defaults, "bias_calibration");
    m_bias_couple_factor = Json2KeyWord<double>(m_defaults, "bias_couple_factor");
    m_bias_scale_mode = Json2KeyWord<std::string>(m_defaults, "bias_scale_mode");
    m_bias_energy_tol = Json2KeyWord<double>(m_defaults, "bias_energy_tol");
}
