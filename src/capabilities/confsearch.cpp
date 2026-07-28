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

#include "src/capabilities/confgen.h"
#include "src/capabilities/confscan.h"
#include "src/capabilities/optimizer_factory.h"
#include "src/capabilities/simplemd.h"

#include "src/core/fileiterator.h"
#include "src/core/molecule.h"

#include "src/tools/general.h"
#include "src/core/parameter_registry.h"
#include "src/core/elements.h"
#include "src/capabilities/rmsd/rmsd_functions.h"  // Claude Generated (Jun 2026): Kabsch helpers for bias calibration

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <map>
#include <mutex>
#include <set>

#include <fmt/core.h>

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
        m_topo_matrix = TopologyMatrix(*mol);
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

    // Claude Generated (Jul 2026): normalize legacy SimpleMD aliases to their canonical names.
    // The CLI stores a flat -hmass as "hmass", but SimpleMD reads the canonical "hydrogen_mass"
    // (simplemd.cpp:229). Layer 1 already inserted the SimpleMD default "hydrogen_mass":1, so
    // both keys sit in md; the ConfigManager merge resolves "hmass"->"hydrogen_mass" only when
    // the canonical key is present, but the iteration order lets the default win, so the user
    // value is silently dropped. (-dt worked only because Layer 6 overwrites md["time_step"].)
    // Rewriting every legacy alias (hmass, velo, dump, print, writeXYZ, nocenter, rm_COM,
    // rmrottrans, initfile, writeinit, writerestart, norestart, rattle_maxiter, ...) to its
    // canonical name here makes the user value reach SimpleMD regardless of merge order.
    // resolveAlias is module-scoped and returns "" for non-SimpleMD names (startT, opt_method,
    // ...), which are passed through unchanged. Edge case: setting both -hydrogen_mass and its
    // alias -hmass lets the alias win here; accepted, that is not a real-world invocation.
    {
        auto& reg = ParameterRegistry::getInstance();
        std::vector<std::pair<std::string, std::string>> remap; // (aliasKey, canonicalKey)
        for (auto& [key, value] : md.items()) {
            if (value.is_object())
                continue;
            std::string canon = reg.resolveAlias("simplemd", key);
            if (!canon.empty() && canon != key)
                remap.emplace_back(key, canon);
        }
        for (auto& [alias, canon] : remap)
            md[canon] = md[alias];
        for (auto& [alias, canon] : remap)
            md.erase(alias);
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
    // Claude Generated (Jul 2026): forward the bias-evaluation speedup controls to each MD run.
    md["rmsd_mtd_max_gaussians"] = m_rmsd_mtd_max_gaussians;
    md["rmsd_mtd_screen"] = m_rmsd_mtd_screen;
    md["rmsd_mtd_cutoff_tol"] = m_rmsd_mtd_cutoff_tol;
    md["rmsd_mtd_screen_margin"] = m_rmsd_mtd_screen_margin;
    // Claude Generated (Jul 2026): strided scheme is inherited from the SimpleMD defaults (Layer 1).
    // Suppress per-child provenance diagnostics in ConfSearch (many MD children would each dump files);
    // provenance is meant for pure -md -rmsd_mtd runs.
    md["rmsd_mtd_diag"] = false;

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
        if (m_seed_pes == "both")
            CurcumaLogger::result_fmt("ConfSearch: next-cycle seeds are picked on BOTH PES -- up to {} from {} plus "
                                      "{} from {} (-seed_pes both)",
                m_seed_rank, m_md_method, m_seed_rank, m_opt_method);
        else
            CurcumaLogger::result_fmt("ConfSearch: next-cycle seeds are picked on the {} PES (-seed_pes {})",
                m_seed_pes == "opt" ? m_opt_method : m_md_method, m_seed_pes);
        // Claude Generated (Jul 2026): state how Phase 3b is run -- same reasoning as the Phase 3c
        // line below, an absent stage message must never be ambiguous.
        if (m_phase3b_two_stage)
            CurcumaLogger::result_fmt("ConfSearch: Phase 3b TWO-STAGE -- crude pre-opt (preset '{}'{}), {}, accurate opt (preset '{}')",
                m_phase3b_preopt_preset,
                m_phase3b_preopt_max_iter > 0 ? fmt::format(", max {} steps", m_phase3b_preopt_max_iter) : "",
                m_phase3b_filter ? "dedup filter" : "no filter", m_phase3b_preset);
        else
            CurcumaLogger::result_fmt("ConfSearch: Phase 3b single-stage (preset '{}') -- enable the crude pre-optimisation with -phase3b_two_stage true",
                m_phase3b_preset);
    } else {
        CurcumaLogger::result("ConfSearch: Single-method mode (Phase 3b skipped)");
    }
    CurcumaLogger::result_fmt("ConfSearch: Temperature={}K -> {}K, step={}K", m_startT, m_endT, m_deltaT);
    CurcumaLogger::result_fmt("ConfSearch: Repetitions={}, RMSD threshold={} A, Energy window={} kJ/mol, Seed rank={}", m_repeat, m_rmsd, m_energy_window, m_seed_rank);
    // Claude Generated (Jul 2026): say whether the recombination phase is armed. It used to be
    // silent when off, so a run simply had no Phase 3c line and there was no way to tell whether the
    // phase had been skipped, had found nothing, or was never enabled.
    if (m_confgen_phase)
        CurcumaLogger::result_fmt("ConfSearch: Phase 3c (torsion recombination) ON -- up to {} proposals "
                                  "per cycle from {} template(s), depth {}",
            m_confgen_max_proposals, m_confgen_templates, m_confgen_depth);
    else
        CurcumaLogger::result("ConfSearch: Phase 3c (torsion recombination) OFF -- enable with -confgen_phase true");
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
        m_topo_matrix = TopologyMatrix(m_topo_ref);
        if (m_restart.topo_ref_opt.AtomCount() > 0) {
            m_topo_ref_opt = m_restart.topo_ref_opt;
            m_topo_matrix_opt = TopologyMatrix(m_topo_ref_opt);
        }
        m_global_min = m_restart.global_min;
        m_initial_energy = m_restart.initial_energy;
        m_initial_energy_opt = m_restart.initial_energy_opt; // opt-PES anchor for the final statistics
        m_best_energy = m_restart.best_energy;
        m_permutation_cache = m_restart.permutations;
        for (auto* mol : m_in_stack) delete mol;
        m_in_stack.clear();
        for (const auto& mm : m_restart.seeds)
            m_in_stack.push_back(new Molecule(mm));
        CurcumaLogger::section("ConfSearch: RESUMING from checkpoint", true);
        CurcumaLogger::result_fmt("ConfSearch: resume T={}K, {} cycles done, {} bias structures, {} seeds, {} cumulative conformers",
            static_cast<int>(m_restart.next_T), m_restart.temperature_cycle,
            static_cast<int>(m_restart.bias.size()), static_cast<int>(m_in_stack.size()),
            static_cast<int>(m_restart.cumulative.size()));
    }

    // Optimise all input structures before any MD run (skipped on resume).
    // m_topo_matrix is updated from the first optimised structure so Phase 4
    // topology checks compare against the relaxed geometry, not the raw input.
    if (!resumed) {
        CurcumaLogger::section("ConfSearch: Initial Geometry Optimisation (" + m_md_method + ")", true);
        // Claude Generated (Jul 2026): stage-named files -- "<base>.initial.<method>[.opt].xyz"
        // instead of the old "<base>.input[.opt].xyz", which said neither what the file is for nor
        // which method produced it (and "input.input.opt.xyz" read like a typo).
        const std::string initial_md = stageBase("initial", m_md_method);
        bool first = true;
        for (auto* mol : m_in_stack) {
            if (first) { mol->writeXYZFile(outputPath(initial_md + ".xyz")); first = false; }
            else          mol->appendXYZFile(outputPath(initial_md + ".xyz"));
        }
        // pre-optimization at md_method ("die md-methode macht die voroptimierung")
        nlohmann::json opt_init = ChildConfig(m_md_method, m_threads);
        PerformOptimisation(initial_md, opt_init);

        for (auto* mol : m_in_stack) delete mol;
        m_in_stack.clear();
        FileIterator opt_in(outputPath(initial_md + ".opt.xyz"));
        while (!opt_in.AtEnd()) {
            Molecule mol = opt_in.Next();
            mol.setCharge(m_charge);
            mol.setSpin(m_spin);
            if (mol.AtomCount() > 0)
                m_in_stack.push_back(new Molecule(mol));
        }
        if (!m_in_stack.empty()) {
            m_topo_matrix = TopologyMatrix(*m_in_stack[0]);
            m_topo_ref = *m_in_stack[0];           // reference structure for restart topology
            m_elements = m_in_stack[0]->Atoms();   // shared atomic-number list for checkpoint frames
        }
        CurcumaLogger::result_fmt("ConfSearch: {} input structures optimised", m_in_stack.size());

        // Claude Generated (Jul 2026): THE md-level reference energy is fixed HERE, at the optimised
        // input structure -- not at the lowest structure of the first metadynamics cycle (which is
        // what the temperature loop used to do). The reference has to answer "what did the search
        // gain over the structure the user handed in", so it must be measured before any MD runs.
        // Taking it after cycle 1 silently absorbed the whole first cycle's gain into the baseline.
        for (const auto* mol : m_in_stack)
            if (mol->AtomCount() > 0 && std::abs(mol->Energy()) > 1e-10)
                m_initial_energy = std::min(m_initial_energy, mol->Energy());
        if (m_initial_energy < std::numeric_limits<double>::infinity()) {
            m_best_energy = m_initial_energy;
            // The optimised input IS a minimum on the md PES, so it also anchors the seed window
            // (which is measured against the running global minimum). Consequence to be aware of:
            // if the input is far below everything a cycle finds, its candidates can fall outside
            // seed_energy_window -- the existing fallback then re-seeds from the input structures,
            // which is the intended behaviour (do not seed MD from clearly worse basins).
            m_global_min = std::min(m_global_min, m_initial_energy);
        }

        // Claude Generated (Jun 2026): dual initial optimization -- in dual-method mode,
        // also optimise the md_method-minimised structures at opt_method so we can track
        // the energy gain on both PES from the very start. The opt_method structures are
        // for reporting only; m_in_stack keeps the md_method structures (they feed the MD loop).
        if (m_opt_method != m_md_method) {
            CurcumaLogger::section("ConfSearch: Initial Geometry Optimisation (" + m_opt_method + ")", true);
            const std::string initial_opt = stageBase("initial", m_opt_method);
            bool first_opt = true;
            for (auto* mol : m_in_stack) {
                if (first_opt) { mol->writeXYZFile(outputPath(initial_opt + ".xyz")); first_opt = false; }
                else              mol->appendXYZFile(outputPath(initial_opt + ".xyz"));
            }
            nlohmann::json opt_hi = ChildConfig(m_opt_method, m_threads);
            PerformOptimisation(initial_opt, opt_hi);

            // Read back opt_method-optimized structures for energy reporting
            std::vector<Molecule*> opt_init_stack;
            FileIterator opt_hi_in(outputPath(initial_opt + ".opt.xyz"));
            while (!opt_hi_in.AtEnd()) {
                Molecule mol = opt_hi_in.Next();
                mol.setCharge(m_charge);
                mol.setSpin(m_spin);
                if (mol.AtomCount() > 0)
                    opt_init_stack.push_back(new Molecule(mol));
            }
            // Claude Generated (Jul 2026): the opt-PES reference, taken from the SAME point in the
            // run as the md-PES one above (the optimised input structure), and over the whole input
            // stack rather than its first frame.
            for (const auto* mol : opt_init_stack)
                if (mol->AtomCount() > 0 && std::abs(mol->Energy()) > 1e-10)
                    m_initial_energy_opt = std::min(m_initial_energy_opt, mol->Energy());

            // Claude Generated (Jul 2026): topology reference for the RANKING PES.
            //
            // The refinement side of Phase 4 checked the opt_method geometries against the
            // md_method reference. Two methods do not agree on bond lengths, so an atom pair
            // sitting near the covalent-radius cutoff can be bonded for one and not for the other.
            // The check then fires for EVERY structure the accurate method produces -- measured on
            // a real run: "opt_method (gfn2) refinement: 0 structures -> cumulative + bias, 38
            // rejected (topo)", i.e. the entire ranking side of that cycle was thrown away while
            // the exploration side kept 10 of the same structures. The reference is now taken from
            // the input structure optimised at opt_method, so like is compared with like.
            if (!opt_init_stack.empty() && opt_init_stack[0]->AtomCount() > 0) {
                m_topo_ref_opt = *opt_init_stack[0];
                m_topo_matrix_opt = TopologyMatrix(m_topo_ref_opt);
                // If the two methods disagree about the topology of the SAME molecule, say so --
                // it explains any later divergence in the reject counts and is worth knowing.
                if (m_topo_matrix.rows() == m_topo_matrix_opt.rows()) {
                    const double diff = (m_topo_matrix - m_topo_matrix_opt).cwiseAbs().sum();
                    if (diff > 1e-4)
                        CurcumaLogger::warn_fmt("ConfSearch: {} and {} disagree about the bond topology of the "
                                                "input structure ({:.0f} bond difference(s)). Each PES is now "
                                                "checked against its own reference.",
                            m_md_method, m_opt_method, diff / 2.0);
                }
            }
            for (auto* mol : opt_init_stack) delete mol;
        }

        // Claude Generated (Jul 2026): report the two reference energies WITHOUT a difference.
        // The old line printed "gfnff=-18.739260 Eh, gfn2=-161.600419 Eh (delta=375081.97 kJ/mol)":
        // the two numbers live on different potential-energy surfaces, so their difference is not a
        // physical quantity at all -- it is dominated by the different zero of energy of the two
        // methods. Only differences WITHIN one method are reported anywhere in this run.
        CurcumaLogger::result("ConfSearch: Reference energies (optimised input structure, one per level of theory):");
        if (m_initial_energy < std::numeric_limits<double>::infinity())
            CurcumaLogger::result_fmt("ConfSearch:   {:<10} (exploration) = {:.6f} Eh", m_md_method, m_initial_energy);
        if (m_opt_method != m_md_method) {
            if (m_initial_energy_opt < std::numeric_limits<double>::infinity())
                CurcumaLogger::result_fmt("ConfSearch:   {:<10} (ranking)     = {:.6f} Eh", m_opt_method, m_initial_energy_opt);
            CurcumaLogger::result("ConfSearch:   the two are different potential-energy surfaces -- no difference between them is reported");
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
        // Claude Generated (Jul 2026): every file this cycle writes carries this tag, so the cycles
        // no longer overwrite each other's intermediates and a listing groups them by temperature.
        m_cycle_tag = fmt::format("cycle{:02d}_T{}K", temperature_cycle, static_cast<int>(m_currentT));
        // Stem of the exploration files: MD snapshots (<stem>.xyz), their optimisation
        // (<stem>.opt.xyz) and the dedup result (<stem>.opt.accepted.xyz).
        const std::string explore = cycleStage("explore", m_md_method);
        // Claude Generated (Jun 2026): per-cycle wall-clock timing
        RunTimer cycle_timer;
        const int entry = pending_entry; // 0 = run MD; 1 = MD already done (resume), skip it
        pending_entry = 0;               // only the first resumed cycle re-enters mid-way
        // Verbosity is scoped by the CurcumaMethod base (ctor captures, dtor restores), so each
        // sub-phase below (Opt/MD/ConfScan) restores this level on its own — no re-assert needed.
        // Claude Generated (Jul 2026): section() instead of header() -- header() is gated at
        // verbosity >= 2, so at the default level the run had no visible block structure at all.
        CurcumaLogger::section(fmt::format("ConfSearch Temperature Cycle {} / {}   T = {} K",
                                   temperature_cycle,
                                   static_cast<int>((m_startT - m_endT) / m_deltaT) + 1,
                                   static_cast<int>(m_currentT)),
            true);
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
                // Claude Generated (Jul 2026): enforce the pool-size cap (rmsd_mtd_max_gaussians) so
                // the per-step bias cost stays bounded as structures accumulate across cycles.
                if (m_rmsd_mtd_max_gaussians > 0) {
                    int removed = m_bias_pool->capToSize(m_rmsd_mtd_max_gaussians);
                    if (removed > 0)
                        CurcumaLogger::result_fmt("ConfSearch: Bias pool capped to {} (rmsd_mtd_max_gaussians), dropped {} low-counter snapshot(s)",
                            m_rmsd_mtd_max_gaussians, removed);
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
                const std::string bias_path = outputPath(explore + ".xyz");
                bool first = true;
                for (const auto& bs : snapshot) {
                    if (bs.persistent)
                        continue; // only raw MD snapshots are optimised in Phase 2 (matches PerformMolecularDynamics)
                    Molecule mol(ref_mol);
                    mol.setGeometry(bs.geometry);
                    mol.setEnergy(bs.energy);
                    mol.setName("bias_" + std::to_string(bs.index));
                    if (first) { mol.writeXYZFile(bias_path); first = false; }
                    else          mol.appendXYZFile(bias_path);
                }
                if (first)
                    std::ofstream(bias_path).close(); // no raw snapshots -> empty input
            }
        }

        CurcumaLogger::section(fmt::format("Phase 2: Geometry Optimisation of Bias Structures ({})", m_md_method));
        // Skip Phase 2 when no new MD ran (empty in_stack going in) and the bias pool did not grow.
        if (no_new_bias_structures) {
            CurcumaLogger::warn_fmt("ConfSearch: T={}K -- no new MD runs and bias pool unchanged -- skipping Phase 2/3.",
                m_currentT);
        } else {
            // Claude Generated (Jul 2026): enforce the reference topology BEFORE optimising.
            // A snapshot that formed or broke a bond is not a conformer of this molecule and is
            // rejected by the Phase 4 filter anyway -- but only after a full optimisation was spent
            // on it, and those are exactly the geometries for which GFN-FF returns a finite energy
            // with a NaN gradient (overlapping atoms in the 1/r hb and three-body terms).
            if (m_snapshot_topology_gate && FilterSnapshotsByTopology(outputPath(explore + ".xyz")) == 0) {
                // Nothing survived: the MD produced only broken structures this cycle. Treat it like
                // a cycle without new bias structures instead of running the phases on an empty
                // file (ConfScan on an empty ensemble is not a path worth relying on).
                CurcumaLogger::warn_fmt("ConfSearch: T={}K -- every MD snapshot changed its topology; "
                                        "skipping Phase 2/3 for this cycle. The dynamics is destroying the "
                                        "molecule -- lower the temperature, or set -topo_check true to abort "
                                        "such MD runs early.",
                    m_currentT);
                no_new_bias_structures = true;
            }

            // Fast per-cycle optimization at md_method.
            // Single-threaded per optimization when ConfSearch parallelizes externally.
            nlohmann::json opt = ChildConfig(m_md_method, (m_threads > 1) ? 1 : m_threads);
            // Bias structures are the primary conformers discovered by RMSD-MTD.
            PerformOptimisation(explore, opt);
            int opt_count = 0;
            {
                FileIterator opt_file(outputPath(explore + ".opt.xyz"));
                while (!opt_file.AtEnd()) { opt_file.Next(); opt_count++; }
            }
            CurcumaLogger::result_fmt("ConfSearch: Optimisation complete. {} bias structures optimised.", opt_count);
            ReportEnsemble("Phase 2 optimised bias structures", m_md_method, outputPath(explore + ".opt.xyz"));
        }
        // Claude Generated (Jun 2026): Phase 2 timing
        CurcumaLogger::result_fmt("ConfSearch: Opt phase took {:.1f} s", cycle_timer.Elapsed() / 1000.0);

        CurcumaLogger::section(fmt::format("Phase 3: RMSD-Based Conformer Filtering ({})", m_md_method));
        int rmsd_count = 0;
        if (no_new_bias_structures) {
            CurcumaLogger::warn_fmt("ConfSearch: T={}K -- skipping Phase 3 (no new structures).", m_currentT);
        } else {
            // "filter between": dedup at md level before the accurate re-opt.
            // Single-threaded per ConfScan when ConfSearch parallelizes externally.
            nlohmann::json scan = FilterConfig(m_md_method, (m_threads > 1) ? 1 : m_threads);
            PerformFilter(explore, scan);
            {
                FileIterator rmsd_file(outputPath(explore + ".opt.accepted.xyz"));
                while (!rmsd_file.AtEnd()) { rmsd_file.Next(); rmsd_count++; }
            }
            CurcumaLogger::result_fmt("ConfSearch: RMSD filtering complete. {} structures accepted.", rmsd_count);
            // Claude Generated (Jul 2026): the accepted count alone says nothing about WHAT was
            // accepted. Report the energies of the surviving ensemble right where the count is.
            ReportEnsemble("Phase 3 accepted ensemble", m_md_method, outputPath(explore + ".opt.accepted.xyz"));
        }

        // Claude Generated (Jul 2026): Phase 3c -- recombine the torsion states of this cycle's
        // minima. Placed after the dedup (so it works on distinct minima) and before Phase 3b (so its
        // structures go through the accurate re-optimisation and the Phase 4 filters like any other).
        if (m_confgen_phase && !no_new_bias_structures && rmsd_count >= 2) {
            CurcumaLogger::section(fmt::format("Phase 3c: Torsion Recombination / ConfGen ({})", m_md_method));
            const int added = PerformConfGen(explore + ".opt.accepted", m_md_method);
            if (added > 0) {
                rmsd_count += added;
                CurcumaLogger::success_fmt("ConfSearch: Phase 3c added {} new conformer(s) that the "
                                           "metadynamics had not found ({} structures now in this cycle)",
                    added, rmsd_count);
                ReportEnsemble("Phase 3c ensemble (metadynamics + recombination)", m_md_method,
                    outputPath(explore + ".opt.accepted.xyz"));
            } else {
                CurcumaLogger::result("ConfSearch: Phase 3c found no new conformer this cycle");
            }
        } else if (m_confgen_phase && no_new_bias_structures) {
            CurcumaLogger::result("ConfSearch: Phase 3c skipped -- no new structures this cycle");
        } else if (m_confgen_phase) {
            CurcumaLogger::result_fmt("ConfSearch: Phase 3c skipped -- {} distinct minimum/minima this "
                                      "cycle, recombination needs at least 2",
                rmsd_count);
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
        // Claude Generated (Jul 2026): the file that holds the opt_method ensemble of this cycle.
        // Not a fixed name any more -- the two-stage mode produces it one optimisation later, so
        // Phase 4 below reads whatever Phase 3b actually wrote.
        std::string hi_level_file;
        if (!no_new_bias_structures && m_opt_method != m_md_method) {
            CurcumaLogger::section(fmt::format("Phase 3b: High-Level Re-Optimisation ({})", m_opt_method));
            hi_level_file = PerformHighLevelOptimisation(explore + ".opt.accepted");
        } else if (m_opt_method == m_md_method) {
            // Claude Generated (Jun 2026): explicit skip notice at result level
            CurcumaLogger::result("ConfSearch: Phase 3b skipped (single-method mode)");
        } else if (no_new_bias_structures) {
            CurcumaLogger::result("ConfSearch: Phase 3b skipped (no new bias structures this cycle)");
        }

        CurcumaLogger::section("Phase 4: Energy Window and Topology Filter");
        for (auto* m : m_in_stack) delete m;
        m_in_stack.clear();
        // EXPLORATION side stays strictly on the md_method (gfnff) PES: seed selection, the
        // exploration global minimum and the bias all read the md_method minima. A region
        // explored by gfnff must never be discarded because opt_method (a different PES, e.g.
        // r2scan) ranks it higher. The opt_method structures are handled in the refinement step
        // below and only feed the FINAL ranking + an extra bias geometry -- their energies are
        // never compared to md_method energies.
        const bool dual_method = (m_opt_method != m_md_method);
        // Claude Generated (Jul 2026): which PES picks the next cycle's seeds. Default "md": the
        // exploration stays on the cheap surface and a region it likes is never discarded because
        // the accurate method ranks it higher. "opt" hands that decision to opt_method -- motivated
        // by measurement: on a 107-atom peptide the two surfaces correlate at only r = 0.40 and the
        // gfn2 minimum sat at gfnff rank 59 of 141, so a gfnff-picked seed is close to a random
        // choice as far as the ranking method is concerned. The MD itself still runs at md_method;
        // only the starting geometries change.
        const bool seed_from_opt = dual_method && (m_seed_pes == "opt" || m_seed_pes == "both");
        const bool seed_from_md = !dual_method || (m_seed_pes != "opt"); // "md" and "both"
        const std::string& cycle_tag = m_cycle_tag;
        const std::string md_accepted = outputPath(explore + ".opt.accepted.xyz");
        double lowest_energy = std::numeric_limits<double>::infinity(); // md_method (exploration)
        int accepted = 0, rejected_topo = 0, rejected_energy = 0;
        std::vector<Molecule*> candidates;
        if (!no_new_bias_structures) {
            FileIterator file(md_accepted);
            int md_index = 0;
            while (!file.AtEnd()) {
                Molecule* mol = new Molecule(file.Next());
                mol->setCharge(m_charge);
                mol->setSpin(m_spin);
                // Claude Generated (Jul 2026): a traceable identity, so a reported seed can be found
                // again in the cycle's ensemble file instead of printing "(unnamed)".
                mol->setName(fmt::format("{}.{}#{:03d}", m_cycle_tag, m_md_method, ++md_index));
                // Topology check: compare bond connectivity (0/1 matrix) against reference.
                // A broken or formed bond changes >=2 entries by 1.0 -> sum >> 1e-4.
                // Log the first mismatched pair to help distinguish GFN-FF artefacts from
                // genuine chemical reactions (proton transfer, ring opening, etc.).
                auto topo_cur = TopologyMatrix(*mol);
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
        // Claude Generated (Jul 2026): per-cycle ensemble + most stable structure on the exploration
        // level. Written here, from the topology-valid candidates, i.e. exactly the set the cycle
        // contributes -- the working files (".bias.opt.accepted.xyz") are overwritten by the next
        // cycle, so without this the per-cycle result was not recoverable after the run.
        if (m_cycle_output && !candidates.empty())
            WriteCycleEnsemble(cycle_tag, m_md_method, candidates);

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

            // Claude Generated (Jul 2026): with seed_pes = opt the next cycle is seeded from the
            // opt_method structures instead (below), so the md_method candidates are only used for
            // the bias feedback and are released here.
            if (!seed_from_md) {
                delete mol;
            } else if ((mol->Energy() - m_global_min) * 2625.5 < eff_seed_window) {
                window_seeds.push_back(mol);
            } else {
                rejected_energy++;
                delete mol;
            }
        }

        // Seed selection: energy ranking (seed_rank) plus, in "diverse" mode, an RMSD spacing
        // requirement so the seeds do not all sit in the same basin. Deletes what it rejects.
        if (seed_from_md) {
            rejected_energy += SelectSeeds(window_seeds, m_md_method, m_global_min);
            for (auto* mol : window_seeds) {
                m_in_stack.push_back(mol);
                accepted++;
            }
        }

        // REFINEMENT side (dual-method only): the opt_method (e.g. gfn2/r2scan) re-optimised
        // structures fill the cumulative pool for the FINAL ranking and add their geometry to the
        // bias pool as a second, independent minimum. Their energies live on the opt_method PES
        // and are NEVER compared to md_method energies -- the cumulative window here is relative
        // to THIS cycle's lowest opt_method energy (same PES), and the next-cycle seeds were
        // already chosen above purely from md_method energies. So a gfnff-explored basin is kept
        // even if opt_method finds it less stable.
        if (dual_method && !no_new_bias_structures && !hi_level_file.empty()) {
            std::vector<Molecule*> opt_candidates;
            double opt_lowest = std::numeric_limits<double>::infinity();
            int opt_rejected_topo = 0, opt_total = 0;
            // Claude Generated (Jul 2026): check against the OPT_METHOD reference (see where it is
            // set). Falls back to the md reference only when there is none (e.g. a resume from a
            // checkpoint written before this existed).
            const Matrix& opt_reference = (m_topo_matrix_opt.rows() == m_topo_matrix.rows())
                ? m_topo_matrix_opt
                : m_topo_matrix;
            FileIterator ofile(hi_level_file);
            while (!ofile.AtEnd()) {
                Molecule* mol = new Molecule(ofile.Next());
                mol->setCharge(m_charge);
                mol->setSpin(m_spin);
                opt_total++;
                mol->setName(fmt::format("{}.{}#{:03d}", m_cycle_tag, m_opt_method, opt_total));
                auto topo_cur = TopologyMatrix(*mol);
                const double topo_diff_sum = (opt_reference - topo_cur).cwiseAbs().sum();
                if (topo_diff_sum > 1e-4) {
                    // Same diagnostic as the exploration side: name the first differing pair. A
                    // silent count is exactly what let "38 of 38 rejected" pass unnoticed.
                    if (opt_rejected_topo == 0) {
                        const int natoms = mol->AtomCount();
                        bool found = false;
                        for (int ii = 0; ii < natoms && !found; ++ii) {
                            for (int jj = ii + 1; jj < natoms; ++jj) {
                                if (std::abs(opt_reference(ii, jj) - topo_cur(ii, jj)) > 0.5) {
                                    CurcumaLogger::warn_fmt("ConfSearch: {} topo reject (first diff): atoms {}-{} "
                                                            "ref_bond={} cur_bond={} (total bond changes: {:.0f})",
                                        m_opt_method, ii, jj,
                                        static_cast<int>(std::round(opt_reference(ii, jj))),
                                        static_cast<int>(std::round(topo_cur(ii, jj))),
                                        topo_diff_sum / 2.0);
                                    found = true;
                                    break;
                                }
                            }
                        }
                    }
                    opt_rejected_topo++;
                    delete mol;
                    continue;
                }
                opt_candidates.push_back(mol);
                opt_lowest = std::min(opt_lowest, mol->Energy());
            }
            // Losing the WHOLE ranking side of a cycle is not a detail -- it means the cumulative
            // pool gets nothing from this temperature and the final ranking silently thins out.
            if (opt_total > 0 && opt_candidates.empty())
                CurcumaLogger::warn_fmt("ConfSearch: ALL {} {} structures of this cycle were topology-rejected -- "
                                        "the ranking side contributes nothing here. If the diff above is a single "
                                        "bond, {} and {} likely disagree about a contact near the covalent-radius "
                                        "cutoff rather than a real reaction.",
                    opt_total, m_opt_method, m_md_method, m_opt_method);
            // Per-cycle ensemble + most stable structure on the ranking level (see the md-level
            // counterpart above). Claude Generated (Jul 2026).
            if (m_cycle_output && !opt_candidates.empty())
                WriteCycleEnsemble(cycle_tag, m_opt_method, opt_candidates);
            // Claude Generated (Jul 2026): running opt-PES minimum, the anchor of the opt-side seed
            // window (mirrors m_global_min on the md side).
            if (opt_lowest < m_global_min_opt)
                m_global_min_opt = opt_lowest;

            std::vector<Molecule*> opt_window_seeds;
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
                if (seed_from_opt && (mol->Energy() - m_global_min_opt) * 2625.5 < eff_seed_window)
                    opt_window_seeds.push_back(mol); // ownership moves to the seed selection
                else
                    delete mol;
            }
            if (seed_from_opt) {
                const int md_seeds = static_cast<int>(m_in_stack.size()); // 0 unless seed_pes=both
                rejected_energy += SelectSeeds(opt_window_seeds, m_opt_method, m_global_min_opt);
                // Claude Generated (Jul 2026): in "both" mode the two lists overlap by construction
                // -- the opt_method structures ARE the re-optimised md_method ones, so the gfn2
                // favourite is usually the gfn2-relaxed version of a gfnff favourite. Seeding the
                // same basin twice wastes an MD run, so an opt seed closer than the seed spacing to
                // an already chosen md seed is dropped.
                const double r_min = (m_seed_min_rmsd > 0.0) ? m_seed_min_rmsd
                                                             : m_seed_diversity_factor * m_rmsd;
                int duplicates = 0, added = 0;
                for (auto* mol : opt_window_seeds) {
                    bool same_basin = false;
                    for (int k = 0; k < md_seeds; ++k)
                        if (PermRMSD(m_in_stack[k]->getGeometry(), mol->getGeometry()) < r_min) {
                            same_basin = true;
                            break;
                        }
                    if (same_basin) {
                        duplicates++;
                        delete mol;
                        continue;
                    }
                    m_in_stack.push_back(mol);
                    accepted++;
                    added++;
                }
                if (m_seed_pes == "both")
                    CurcumaLogger::result_fmt("ConfSearch: next cycle seeded from BOTH PES: {} from {} + {} from {} "
                                              "= {} seed(s){}",
                        md_seeds, m_md_method, added, m_opt_method, static_cast<int>(m_in_stack.size()),
                        duplicates > 0 ? fmt::format(" ({} {} seed(s) dropped -- same basin as a {} seed, "
                                                     "closer than {:.2f} A)",
                                             duplicates, m_opt_method, m_md_method, r_min)
                                       : "");
                else
                    CurcumaLogger::result_fmt("ConfSearch: next cycle seeded from the {} PES ({} seed(s)) -- "
                                              "-seed_pes opt",
                        m_opt_method, static_cast<int>(m_in_stack.size()));
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
        // Claude Generated (Jul 2026): the reference is the optimised INPUT structure, set before the
        // first MD ran (see the initial optimisation). The first cycle is therefore reported like
        // every other one -- it used to overwrite the reference with its own result, which made the
        // gain of cycle 1 invisible by construction. Only a pre-optimisation that produced no usable
        // energy at all falls back to the old behaviour.
        if (!std::isfinite(initial_energy) && std::isfinite(lowest_energy)) {
            initial_energy = lowest_energy;
            best_energy = lowest_energy;
            CurcumaLogger::warn_fmt("ConfSearch: no energy from the initial optimisation -- using this cycle's lowest {} structure as reference ({:.6f} Eh)",
                m_md_method, initial_energy);
        }
        if (lowest_energy < std::numeric_limits<double>::infinity()) {
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

        CurcumaLogger::result_fmt("ConfSearch: === End Temperature Cycle {} (T = {} K) ===",
            temperature_cycle, static_cast<int>(m_currentT));
        CurcumaLogger::blank_line(); // visual separation between cycles
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
    CurcumaLogger::section("ConfSearch: Final Deduplication Pass", true);
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
            CurcumaLogger::section("ConfSearch: Final Energy Statistics", true);
            CurcumaLogger::result_fmt("ConfSearch: {} unique conformer(s); global minimum {:.6f} Eh",
                static_cast<int>(energies.size()), e_min);
            CurcumaLogger::result_fmt("ConfSearch: energy span {:.2f} kJ/mol (lowest {:.6f} Eh, highest kept {:.6f} Eh)",
                span_kj, e_min, e_max);
            // Claude Generated (Jul 2026): e_min is the lowest energy of the ACCEPTED pool, and that
            // pool is filled on the opt_method PES in dual mode (Phase 4 refinement side) and on the
            // md_method PES otherwise. Absolute energies of two different methods are not comparable
            // -- the old code always subtracted the md_method initial energy from e_min, so a dual
            // run reported nonsense like "lowered the energy by 199567 kJ/mol (gfnff: -9.005 ->
            // -85.017 Eh)". Each PES is now only ever compared against itself.
            const bool dual_run = (m_opt_method != m_md_method);
            if (!dual_run) {
                if (initial_energy < std::numeric_limits<double>::infinity()) {
                    const double gain_kj = (initial_energy - e_min) * 2625.5;
                    if (gain_kj > 1e-3)
                        CurcumaLogger::success_fmt("ConfSearch: search lowered the energy by {:.2f} kJ/mol vs. the initial structure ({}: {:.6f} -> {:.6f} Eh)",
                            gain_kj, m_md_method, initial_energy, e_min);
                    else
                        CurcumaLogger::result_fmt("ConfSearch: initial structure remains the global minimum ({:.6f} Eh)", e_min);
                }
            } else {
                // Exploration side: md_method initial -> md_method running best (both on the md PES).
                // This narrates the gfnff search; it is NOT the ranking and must not touch e_min.
                if (initial_energy < std::numeric_limits<double>::infinity()
                    && best_energy < std::numeric_limits<double>::infinity()) {
                    const double explore_gain_kj = (initial_energy - best_energy) * 2625.5;
                    CurcumaLogger::result_fmt("ConfSearch: exploration ({}) lowered its own minimum by {:.2f} kJ/mol ({:.6f} -> {:.6f} Eh) -- separate PES, not comparable to the {} ranking",
                        m_md_method, explore_gain_kj, initial_energy, best_energy, m_opt_method);
                }
                // Ranking side: opt_method initial -> opt_method global minimum (both on the opt PES).
                if (m_initial_energy_opt < std::numeric_limits<double>::infinity()) {
                    const double opt_gain_kj = (m_initial_energy_opt - e_min) * 2625.5;
                    if (opt_gain_kj > 1e-3)
                        CurcumaLogger::success_fmt("ConfSearch: search lowered the energy by {:.2f} kJ/mol vs. the initial structure ({}: {:.6f} -> {:.6f} Eh)",
                            opt_gain_kj, m_opt_method, m_initial_energy_opt, e_min);
                    else
                        CurcumaLogger::result_fmt("ConfSearch: initial structure remains the global minimum ({}: {:.6f} Eh)",
                            m_opt_method, e_min);
                }
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
    // Claude Generated (Jul 2026): MD-module output is verbosity-2+. At verbosity 1 each parallel
    // MD run is silenced (its per-step tables would interleave across threads anyway) and replaced
    // by a compact "X/N runs done" counter with per-run + cumulative timing, emitted from a run-
    // completion hook. Uses fmt::print (not CurcumaLogger) because the global logger level is 0
    // inside the workers, which would swallow result()/result_fmt.
    const int total_runs = m_repeat * static_cast<int>(molecules.size());
    // Claude Generated (Jul 2026): at verbosity 1 the completion hook now drives a single live bar
    // (stdout, honours -noprogress/non-TTY) instead of one scrolling line per finished run. Without
    // a TTY the line form is kept, so a redirected log still records the progress.
    const bool live_bar = (m_verbosity == 1) && CurcumaLogger::progress_enabled();
    const bool report_counter = (m_verbosity == 1);
    const std::string bar_label = fmt::format("MD runs [{}]", m_md_method);
    nlohmann::json md_param = parameter;
    md_param["verbosity"] = (m_verbosity >= 2) ? m_verbosity : 0;

    std::atomic<int> done_runs{ 0 };
    std::mutex counter_mtx;
    const auto phase_start = std::chrono::steady_clock::now();

    CxxThreadPool* pool = new CxxThreadPool;
    int index = 0;
    CurcumaLogger::result_fmt("ConfSearch MD: Starting {} independent runs ({} repeats x {} structures), {} threads in parallel",
        total_runs, m_repeat, static_cast<int>(molecules.size()), m_threads);
    for (int repeat = 0; repeat < m_repeat; ++repeat) {
        for (size_t i = 0; i < molecules.size(); ++i) {
            MDThread* thread = new MDThread(md_param);
            thread->setThreadId(index++);
            thread->setBasename(outputPath(Basename() + ".r" + std::to_string(repeat)));
            thread->setMolecule(molecules[i]);
            thread->setSharedBiasPool(m_bias_pool);
            if (report_counter) {
                thread->setOnComplete([&](double run_seconds) {
                    const int d = ++done_runs;
                    const double elapsed = std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - phase_start)
                                               .count();
                    const int running = std::max(0, std::min(m_threads, total_runs - d));
                    std::lock_guard<std::mutex> lock(counter_mtx);
                    if (live_bar) {
                        CurcumaLogger::progress_bar(d, total_runs, bar_label);
                    } else {
                        fmt::print("  MD runs: {}/{} done | ~{} running | last {:.2f} s | elapsed {:.1f} s\n",
                            d, total_runs, running, run_seconds, elapsed);
                        fflush(stdout);
                    }
                });
            }
            pool->addThread(thread);
        }
    }
    pool->setActiveThreadCount(m_threads);
    // Claude Generated (Jul 2026): the CxxThreadPool stderr bar is per-run detail -- keep only at
    // verbosity 3. The phase lines here + the verbosity-1 run counter above give progress otherwise.
    pool->setProgressBar(m_verbosity >= 3 ? CxxThreadPool::ProgressBarType::Continously
                                          : CxxThreadPool::ProgressBarType::None);
    pool->StartAndWait();
    // Restore the orchestrator's level before anything is reported (the workers leave the shared
    // static logger level at their own), then close the bar line.
    CurcumaLogger::set_verbosity(m_verbosity);
    if (live_bar)
        CurcumaLogger::progress_done();

    const double total_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - phase_start)
                                     .count();
    CurcumaLogger::result_fmt("ConfSearch: {} MD runs finished in {:.1f} s. Bias pool: {} structures.",
        index, total_seconds, m_bias_pool ? m_bias_pool->biasStructureCount() : 0);

    // Claude Generated (Jul 2026): removed the empty "confsearch.unique.xyz" stub that was
    // created here. It hardcoded the basename instead of Basename(), was always written empty
    // and was never read back -- the bias pool is the only conformer source (see Phase 2).

    // Export the bias pool. Only raw MD snapshots (persistent=false) go into the file that Phase 2
    // optimises -- persistent structures are already-optimised fed-back minima and must not be
    // re-optimised every cycle (would cause false "New best!" triggers via numerical noise).
    // The bias-pool file gets the full unsampled pool (including persistent) for inspection.
    //
    // Claude Generated (Jul 2026): stage-named files. "<base>.explore.<method>.xyz" (the snapshots
    // Phase 2 consumes) and "<base>.bias_pool.<method>.xyz" (the full pool, was ".mtd.xyz").
    const std::string snapshot_file = cycleStage("explore", m_md_method) + ".xyz";
    const std::string pool_file = cycleStage("bias_pool", m_md_method) + ".xyz";
    if (m_bias_pool && m_bias_pool->biasStructureCount() > 0 && !m_in_stack.empty()) {
        auto snapshot = m_bias_pool->snapshot();
        const Molecule& ref_mol = *m_in_stack[0];

        // full pool for inspection
        bool first = true;
        for (const auto& bs : snapshot) {
            Molecule mol(ref_mol);
            mol.setGeometry(bs.geometry);
            // Claude Generated (Jul 2026): the copy inherits ref_mol's ENERGY together with its atom
            // list. Without this the exported snapshots all carried the first seed's energy -- a
            // wrong but plausible number attached to a different geometry (measured: a snapshot
            // labelled -8.985207 Eh is really worth -8.840638 Eh).
            mol.setEnergy(bs.energy);
            mol.setName("bias_" + std::to_string(bs.index) + " t=" + std::to_string(static_cast<int>(bs.time)));
            if (first) { mol.writeXYZFile(outputPath(pool_file)); first = false; }
            else          mol.appendXYZFile(outputPath(pool_file));
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
            mol.setEnergy(new_snapshots[i].energy); // see the .mtd.xyz export above
            mol.setName("bias_" + std::to_string(new_snapshots[i].index));
            if (first) { mol.writeXYZFile(outputPath(snapshot_file)); first = false; }
            else          mol.appendXYZFile(outputPath(snapshot_file));
            exported++;
        }
        CurcumaLogger::result_fmt("ConfSearch: {} new MD snapshots (of {} total pool, {} persistent skipped, stride={}) written to {}",
            exported, static_cast<int>(snapshot.size()), static_cast<int>(snapshot.size()) - bias_count, stride, snapshot_file);
    }

    delete pool; // MDThread sets autoDelete=true, so the pool frees the threads here (no leak)

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

    // Claude Generated (Jul 2026): freeze the GFN-FF topology for every child computation.
    //
    // A conformational search explores ONE molecule: if the bonding topology changes, the structure
    // is a reaction product and is rejected downstream anyway. The default "auto" mode re-derives the
    // topology mid-run whenever an atom moved more than 0.26 Bohr, and every re-derivation shifts the
    // energy SCALE (new topology charges / coordination numbers). During a long optimisation of a hot
    // MD snapshot the optimiser then follows those jumps instead of the physics: measured on a
    // penta-alanine snapshot, an optimisation "converged" to -9.168083 Eh while its final geometry is
    // really worth -8.668213 Eh (1312 kJ/mol apart, identical bond topology). One such structure in
    // the pool becomes the energy reference and pushes every real conformer out of the energy window
    // -- the run then reports "1 unique conformer" out of 482. With a frozen topology the same
    // optimisation lands on -8.945487 Eh, matching an independent optimisation of the same basin to
    // 0.1 kJ/mol. An explicit -gfnff.topology_mode still wins.
    if (!(cfg.contains("gfnff") && cfg["gfnff"].is_object() && cfg["gfnff"].contains("topology_mode")))
        cfg["gfnff"]["topology_mode"] = "constant";

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

    // Claude Generated (Jul 2026): a cycle can produce no structures to optimise (e.g. -rmsd_mtd
    // false, or an MD phase that deposited nothing), so the input file may not exist. FileIterator
    // would throw an uncaught std::runtime_error ("File not found") -> terminate. Skip gracefully.
    {
        std::ifstream check(input_file);
        if (!check.good()) {
            CurcumaLogger::warn("PerformOptimisation: no input file (nothing to optimise), skipping: " + input_file);
            return basename;
        }
    }

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
    // Claude Generated (Jul 2026): structures the method could not handle are written out for
    // analysis instead of only being counted. What is saved is the INPUT geometry -- that is the
    // thing that makes the force field produce NaN, and it is what one wants to reload, inspect and
    // reduce to a bug report. The optimised geometry is useless here (it is the diverged one).
    // "Failed" deliberately does NOT include plain non-convergence: ConfSearch accepts unconverged
    // structures on purpose, so writing those would dump most of the batch.
    const std::string failed_file = outputPath(basename + ".failed.xyz");
    int failed_written = 0;
    auto save_failed = [&](const Molecule& input, int idx, const std::string& reason) {
        Molecule copy = input;
        copy.setName(fmt::format("failed_{:03d}_{} [{}]", idx + 1, reason, opt_method_name));
        if (failed_written == 0)
            copy.writeXYZFile(failed_file);
        else
            copy.appendXYZFile(failed_file);
        failed_written++;
    };

    int written = 0;
    auto write_result = [&](const Optimization::OptimizationResult& res, const Molecule& fallback, int idx) {
        double e_start = res.energy_trajectory.empty() ? 0.0 : res.energy_trajectory.front();
        double e_end   = res.final_energy;
        double dE_kjmol = (e_end - e_start) * 2625.5;  // Eh -> kJ/mol
        // Claude Generated (Jul 2026): never write a diverged structure into the pool. A geometry
        // with NaN/Inf coordinates passes every downstream numeric test (comparisons with NaN are
        // false), so it would travel through the dedup into the cumulative ensemble and poison the
        // energy reference. The input geometry is the honest fallback.
        const bool finite_result = res.final_molecule.AtomCount() > 0
            && res.final_molecule.getGeometry().allFinite() && std::isfinite(e_end);
        if (finite_result) {
            res.final_molecule.appendXYZFile(output_file);
            // Claude Generated (Jul 2026): the per-structure table is the detail record and moves to
            // verbosity >= 2. At verbosity 1 the live progress bar during the batch plus the summary
            // line below carry the information -- a 400-structure batch printed 400 lines here.
            if (m_verbosity >= 2)
                CurcumaLogger::result_fmt("  Struct {:2d} [{}]: {:4d} steps, E = {:+.6f} Eh, dE = {:+.2f} kJ/mol{}",
                    idx + 1, opt_method_name, res.iterations_performed, e_end, dE_kjmol,
                    res.success ? "" : "  (not converged)");
            ++written;
        } else if (fallback.AtomCount() > 0 && fallback.getGeometry().allFinite()) {
            fallback.appendXYZFile(output_file);
            CurcumaLogger::warn_fmt("  Struct {:2d}: optimisation unusable ({}), using the input geometry",
                idx + 1,
                res.final_molecule.AtomCount() == 0 ? "no geometry returned" : "non-finite energy/geometry");
            ++written;
        } else {
            CurcumaLogger::warn_fmt("  Struct {:2d}: dropped -- neither the optimised nor the input geometry is finite",
                idx + 1);
        }
    };

    // Unified path: CxxThreadPool runs the optimizations inline on the calling thread when
    // m_threads == 1 (no worker spawned) and in parallel otherwise.
    //
    // Claude Generated (Jul 2026): use-after-free fix. The CxxThreadPool destructor iterates its
    // internal thread queues and reads each thread's m_autodelete (CxxThreadPool.h ~185). OptThread
    // sets autoDelete=false so the pool must not free them -- but deleting the threads while the pool
    // still references them makes that destructor read freed memory: a garbage (non-zero) m_autodelete
    // makes it call the virtual destructor through a corrupted vtable -> SIGSEGV (observed as a crash
    // in ConfSearch::PerformOptimisation). So the pool must be destroyed FIRST, while the threads are
    // still alive (AutoDelete()==false is read correctly, nothing is freed), and only then are the
    // threads deleted here. Scope the pool to enforce that order.
    // Claude Generated (Jul 2026): live batch progress.
    //
    // The pool here is ONE pool for ONE batch (all structures are queued, one StartAndWait), unlike
    // ConfScan's reorder/check pools which are Reset()+StartAndWait() per candidate -- that is the
    // per-permutation flood the CxxThreadPool bar was suppressed for, and it stays suppressed there.
    // The pool's own bar is still not used here either: it writes percent-only to stderr, carries no
    // "structure i of N" and ignores -noprogress / non-TTY. The OptThread completion hook feeds
    // CurcumaLogger::progress_bar instead (stdout, labelled, i/N, honours -noprogress).
    std::atomic<int> done_structures{ 0 };
    std::mutex progress_mtx;
    const auto batch_start = std::chrono::steady_clock::now();
    const bool live_bar = (m_verbosity == 1) && CurcumaLogger::progress_enabled();
    const std::string bar_label = fmt::format("Optimising [{}]", opt_method_name);

    std::vector<OptThread*> threads;
    threads.reserve(total);
    int failed = 0;
    {
        CxxThreadPool pool;
        pool.setActiveThreadCount(m_threads);
        // The pool's own stderr bar stays off below verbosity 3 (see above); progress comes from the
        // completion hook.
        pool.setProgressBar(m_verbosity >= 3 ? CxxThreadPool::ProgressBarType::Continously
                                             : CxxThreadPool::ProgressBarType::None);
        for (int i = 0; i < total; ++i) {
            OptThread* t = new OptThread(molecules[i], local_param);
            t->setOnComplete([&, i](double seconds, const Optimization::OptimizationResult& res) {
                const int d = ++done_structures;
                std::lock_guard<std::mutex> lock(progress_mtx);
                if (live_bar) {
                    CurcumaLogger::progress_bar(d, total, bar_label);
                } else if (m_verbosity >= 2) {
                    // Raw fmt::print, not CurcumaLogger: this runs on a pool worker where the shared
                    // static logger level is 0 (see OptThread::execute).
                    const double elapsed = std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - batch_start)
                                               .count();
                    fmt::print("  Opt [{}]: {}/{} done | struct {}: {} steps, E = {:.6f} Eh{} | {:.2f} s (elapsed {:.1f} s)\n",
                        opt_method_name, d, total, i + 1, res.iterations_performed, res.final_energy,
                        res.success ? "" : ", not converged", seconds, elapsed);
                    std::fflush(stdout);
                }
            });
            pool.addThread(t);
            threads.push_back(t);
        }
        CurcumaLogger::result_fmt("Optimizing {} structures using {} thread(s) [{}]", total, m_threads, opt_method_name);
        pool.StartAndWait();
        // Claude Generated (Jul 2026): restore the level HERE, not at the end of the function.
        // OptThread::execute() sets the shared static logger level to the worker level (0), so
        // everything below -- the per-structure lines AND the batch summary -- was swallowed at
        // verbosity 1: a 43-structure batch printed "Optimizing 43 structures ..." and then nothing
        // at all, neither during nor after.
        CurcumaLogger::set_verbosity(m_verbosity);
        if (live_bar)
            CurcumaLogger::progress_done();
        for (int i = 0; i < static_cast<int>(threads.size()); ++i) {
            const Optimization::OptimizationResult& res = threads[i]->result();
            if (!res.success) ++failed;
            // Save the pathological ones (method NaN / no geometry / non-finite result), not the
            // merely unconverged ones.
            const bool no_geometry = res.final_molecule.AtomCount() == 0;
            const bool non_finite = !no_geometry
                && (!res.final_molecule.getGeometry().allFinite() || !std::isfinite(res.final_energy));
            if (threads[i]->methodReportedNaN() || no_geometry || non_finite)
                save_failed(molecules[i], i,
                    threads[i]->methodReportedNaN() ? "method_nan" : (no_geometry ? "no_geometry" : "non_finite"));
            write_result(res, molecules[i], i);
        }
    } // pool destroyed here while the threads are still alive -> no use-after-free
    for (auto* t : threads)
        delete t;

    const double batch_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - batch_start)
                                     .count();
    if (failed_written > 0)
        CurcumaLogger::warn_fmt("ConfSearch: {} structure(s) the method could not handle written to {} "
                                "(input geometries, named with the reason) -- reproduce with: "
                                "curcuma -sp <structure> -method {}",
            failed_written, failed_file, opt_method_name);
    if (failed > 0)
        CurcumaLogger::result_fmt("Optimization batch complete [{}]: {}/{} structures written in {:.1f} s ({} failed: zero step / gradient failure)",
            opt_method_name, written, total, batch_seconds, failed);
    else
        CurcumaLogger::result_fmt("Optimization batch complete [{}]: {}/{} structures written in {:.1f} s to {}",
            opt_method_name, written, total, batch_seconds, output_file);
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

// Claude Generated (Jul 2026): Phase 3b -- re-optimisation of the cycle's deduplicated minima at
// opt_method.
//
// Default (single stage): one accurate optimisation per structure, as before.
//
// Two-stage (-phase3b_two_stage true): the accurate optimisation is the most expensive step of a
// dual-method run, and it is paid for structures that are not distinct minima at that level -- the
// md_method dedup ran on a DIFFERENT potential-energy surface, so structures it kept apart can
// collapse onto one another as soon as the accurate method relaxes them. The two-stage mode
// therefore relaxes everything crudely first (loose convergence preset), deduplicates at that level
// and only then optimises the survivors accurately. The crude stage is not a shortcut in accuracy:
// its structures are thrown away, only the SELECTION it produces is used.
//
// Claude Generated (Jul 2026): stage-named files. Each stage COPIES its input to a file whose name
// states purpose and method instead of chaining another suffix onto the previous stage's name:
//   single stage: <base>.refine.<opt_method>.xyz -> .opt.xyz
//   two stage:    <base>.refine_crude.<opt_method>.xyz -> .opt.xyz -> .opt.accepted.xyz
//                 <base>.refine.<opt_method>.xyz       -> .opt.xyz
// The copy costs one file write per cycle and is what keeps a dual-method run readable -- the old
// chain ended in "<base>.bias.opt.accepted.opt.accepted.opt.xyz", which named neither the method
// nor the stage.
std::string ConfSearch::PerformHighLevelOptimisation(const std::string& f)
{
    const int child_threads = (m_threads > 1) ? 1 : m_threads;
    nlohmann::json opt_accurate = ChildConfig(m_opt_method, child_threads);
    opt_accurate["convergence_preset"] = m_phase3b_preset;

    auto count_frames = [](const std::string& path) {
        std::ifstream check(path);
        if (!check.good())
            return 0;
        int n = 0;
        FileIterator it(path);
        while (!it.AtEnd()) { it.Next(); n++; }
        return n;
    };

    const std::string refine = cycleStage("refine", m_opt_method);

    if (!m_phase3b_two_stage) {
        CopyFrames(outputPath(f + ".xyz"), outputPath(refine + ".xyz"));
        PerformOptimisation(refine, opt_accurate); // reads "<refine>.xyz", writes "<refine>.opt.xyz"
        const std::string out = outputPath(refine + ".opt.xyz");
        CurcumaLogger::result_fmt("ConfSearch: High-level re-optimisation complete. {} structures re-optimised at {} (preset '{}').",
            count_frames(out), m_opt_method, m_phase3b_preset);
        ReportEnsemble("Phase 3b re-optimised ensemble", m_opt_method, out);
        return out;
    }

    // --- Stage 1: crude optimisation of every structure ---
    const std::string crude = cycleStage("refine_crude", m_opt_method);
    nlohmann::json opt_crude = ChildConfig(m_opt_method, child_threads);
    opt_crude["convergence_preset"] = m_phase3b_preopt_preset;
    if (m_phase3b_preopt_max_iter > 0)
        opt_crude["max_iterations"] = m_phase3b_preopt_max_iter;
    CurcumaLogger::result_fmt("ConfSearch: Phase 3b stage 1/2 -- crude pre-optimisation at {} (preset '{}'{})",
        m_opt_method, m_phase3b_preopt_preset,
        m_phase3b_preopt_max_iter > 0 ? fmt::format(", max {} steps", m_phase3b_preopt_max_iter) : "");
    CopyFrames(outputPath(f + ".xyz"), outputPath(crude + ".xyz"));
    PerformOptimisation(crude, opt_crude);
    const std::string crude_file = outputPath(crude + ".opt.xyz");
    const int n_crude = count_frames(crude_file);
    if (n_crude == 0) {
        CurcumaLogger::warn("ConfSearch: Phase 3b stage 1 produced no structures -- skipping the accurate stage");
        return crude_file;
    }
    ReportEnsemble("Phase 3b stage 1 (crude)", m_opt_method, crude_file);

    // --- Stage 2: dedup at the crude level, so the accurate stage only pays for distinct minima ---
    std::string accurate_source = crude_file; // no filter: the crude output goes on directly
    if (m_phase3b_filter) {
        nlohmann::json scan = FilterConfig(m_opt_method, child_threads);
        // The intermediate filter deduplicates, it does NOT rank: its input carries CRUDE energies,
        // and a structure that sits high after a loose relaxation can still fall into the window
        // once it is optimised properly. Discarding it here on that number would lose a conformer
        // for good. -1 disables ConfScan's energy cutoff; the real energy window is applied in
        // Phase 4, on fully optimised energies.
        scan["max_energy"] = -1.0;
        PerformFilter(crude, scan); // reads "<crude>.opt.xyz", writes "<crude>.opt.accepted.xyz"
        const std::string filtered = outputPath(crude + ".opt.accepted.xyz");
        const int n_filtered = count_frames(filtered);
        if (n_filtered > 0) {
            accurate_source = filtered;
            CurcumaLogger::result_fmt("ConfSearch: Phase 3b filter at {}: {} of {} crude structures are distinct minima ({} accurate optimisation(s) saved)",
                m_opt_method, n_filtered, n_crude, n_crude - n_filtered);
        } else {
            CurcumaLogger::warn("ConfSearch: Phase 3b filter accepted nothing -- optimising all crude structures accurately");
        }
    }

    // --- Stage 3: accurate optimisation of the survivors ---
    CurcumaLogger::result_fmt("ConfSearch: Phase 3b stage 2/2 -- accurate optimisation at {} (preset '{}')",
        m_opt_method, m_phase3b_preset);
    CopyFrames(accurate_source, outputPath(refine + ".xyz"));
    PerformOptimisation(refine, opt_accurate);
    const std::string out = outputPath(refine + ".opt.xyz");
    CurcumaLogger::result_fmt("ConfSearch: High-level re-optimisation complete. {} structures re-optimised at {}.",
        count_frames(out), m_opt_method);
    ReportEnsemble("Phase 3b re-optimised ensemble", m_opt_method, out);
    return out;
}

// Claude Generated (Jul 2026): bond matrix with an explicit factor -- see the header for why 1.5 is
// unusable for this decision.
Matrix ConfSearch::TopologyMatrix(const Molecule& mol) const
{
    const int n = mol.AtomCount();
    Matrix topo = Matrix::Zero(n, n);
    const Geometry geom = mol.getGeometry();
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            const double cutoff = m_topology_factor
                * (Elements::CovalentRadius[mol.Atom(i).first] + Elements::CovalentRadius[mol.Atom(j).first]);
            if ((Eigen::Vector3d(geom.row(i)) - Eigen::Vector3d(geom.row(j))).norm() <= cutoff) {
                topo(i, j) = 1.0;
                topo(j, i) = 1.0;
            }
        }
    }
    return topo;
}

// Claude Generated (Jul 2026): restrained repair of a snapshot -- see the header for the rationale.
bool ConfSearch::RepairSnapshot(Molecule& mol, EnergyCalculator& calculator) const
{
    const Matrix topo = TopologyMatrix(mol);
    if (topo.rows() != m_topo_matrix.rows())
        return false;

    // One restraint per offending pair. A MISSING bond is pulled back to the length it has in the
    // reference structure; a SPURIOUS bond is pushed just outside bonding range (the covalent sum
    // times the same 1.5 factor the bond criterion uses, plus a margin) so the contact opens up.
    std::vector<Optimization::DistanceRestraint> restraints;
    const Geometry ref_geom = m_topo_ref.getGeometry();
    const bool have_reference_geometry = (ref_geom.rows() == mol.AtomCount());
    for (int i = 0; i < mol.AtomCount(); ++i) {
        for (int j = i + 1; j < mol.AtomCount(); ++j) {
            const double d = topo(i, j) - m_topo_matrix(i, j);
            if (std::abs(d) < 0.5)
                continue;
            Optimization::DistanceRestraint r;
            r.i = i;
            r.j = j;
            r.force = m_repair_force;
            const double covalent_sum = Elements::CovalentRadius[mol.Atom(i).first]
                + Elements::CovalentRadius[mol.Atom(j).first];
            if (d < 0) {
                // bond missing -> target = its length in the reference structure
                r.target = have_reference_geometry
                    ? (Eigen::Vector3d(ref_geom.row(i)) - Eigen::Vector3d(ref_geom.row(j))).norm()
                    : covalent_sum;
            } else {
                // bond formed -> target just beyond the bond criterion (1.5 * covalent sum)
                r.target = 1.5 * covalent_sum * 1.15;
            }
            restraints.push_back(r);
        }
    }
    if (restraints.empty() || static_cast<int>(restraints.size()) > m_repair_max_bonds)
        return false; // nothing to do, or too broken to be a conformer of this molecule

    nlohmann::json cfg = ChildConfig(m_md_method, 1);
    cfg["verbosity"] = 0;
    cfg["write_trajectory"] = false;

    // Stage 1: restrained -- pull the offending pairs back while the rest of the molecule follows.
    nlohmann::json restrained = cfg;
    restrained["max_iterations"] = m_repair_max_iterations;
    restrained["distance_restraints"] = Optimization::GeometryRestraints::toJson(restraints);
    Molecule work = mol;
    auto stage1 = Optimization::OptimizationDispatcher::optimizeStructure(
        &work, Optimization::OptimizerType::LBFGSPP, &calculator, restrained);
    if (stage1.final_molecule.AtomCount() == 0)
        return false;

    // Stage 2: RELEASED -- the reported structure must be a minimum of the plain force field.
    Molecule free_mol = stage1.final_molecule;
    auto stage2 = Optimization::OptimizationDispatcher::optimizeStructure(
        &free_mol, Optimization::OptimizerType::LBFGSPP, &calculator, cfg);
    if (stage2.final_molecule.AtomCount() == 0)
        return false;

    Molecule repaired = stage2.final_molecule;
    repaired.setCharge(m_charge);
    repaired.setSpin(m_spin);
    // Claude Generated (Jul 2026): validate BEFORE the topology comparison. A NaN coordinate makes
    // the difference NaN, and `NaN > 1e-4` is FALSE -- a diverged repair would be accepted as
    // "topology restored" and handed to the next optimisation, where every GFN-FF term returns NaN.
    if (!repaired.getGeometry().allFinite() || !std::isfinite(stage2.final_energy))
        return false;
    if ((TopologyMatrix(repaired) - m_topo_matrix).cwiseAbs().sum() > 1e-4)
        return false; // still not this molecule -- the break was real chemistry, not an artefact

    mol = repaired;
    return true;
}

// Claude Generated (Jul 2026): pre-optimisation topology gate -- see the header for the rationale.
int ConfSearch::FilterSnapshotsByTopology(const std::string& path) const
{
    std::ifstream check(path);
    if (!check.good() || m_topo_matrix.rows() == 0)
        return 0;

    std::vector<Molecule> keep;
    std::vector<Molecule> broken; // candidates for the restrained repair (m_repair_snapshots)
    int total = 0, rejected_formed = 0, rejected_broken = 0, rejected_clash = 0, rejected_nonfinite = 0;
    bool reported = false, reported_clash = false;
    {
        FileIterator it(path);
        while (!it.AtEnd()) {
            Molecule mol = it.Next();
            if (mol.AtomCount() == 0)
                continue;
            total++;
            // Claude Generated (Jul 2026): a geometry with NaN/Inf coordinates would SILENTLY PASS
            // the topology test -- every comparison with NaN is false, so the loop below finds "no
            // difference" and keeps the structure. It then reaches GFN-FF, whose energy and every
            // gradient term come back NaN. Screen it here.
            if (!mol.getGeometry().allFinite()) {
                rejected_nonfinite++;
                continue;
            }
            // Claude Generated (Jul 2026): a clash that does NOT change the topology still breaks
            // the force field. The 1/r factors in the hb and three-body terms blow up when any pair
            // gets close enough, and for a pair that is ALREADY BONDED in the reference no bond is
            // "formed", so the topology test above sees nothing wrong. Observed exactly that: the
            // gate passed structures whose optimisation then reported NaN in hb/batm/atm while the
            // energy stayed finite. Screen on the raw distance instead: anything below
            // clash_ratio x (sum of covalent radii) is a collapsed geometry, bonded or not.
            if (m_snapshot_clash_ratio > 0.0) {
                const Geometry geom = mol.getGeometry();
                bool collapsed = false;
                for (int i = 0; i < mol.AtomCount() && !collapsed; ++i) {
                    for (int j = i + 1; j < mol.AtomCount(); ++j) {
                        const double limit = m_snapshot_clash_ratio
                            * (Elements::CovalentRadius[mol.Atom(i).first]
                                + Elements::CovalentRadius[mol.Atom(j).first]);
                        if ((Eigen::Vector3d(geom.row(i)) - Eigen::Vector3d(geom.row(j))).norm() < limit) {
                            if (!reported_clash) {
                                CurcumaLogger::warn_fmt("ConfSearch: clash screen (first reject): atoms {}-{} are "
                                                        "{:.3f} A apart, below {:.0f}% of their covalent radii sum "
                                                        "({:.3f} A) -- the force field's 1/r terms break there",
                                    i, j,
                                    (Eigen::Vector3d(geom.row(i)) - Eigen::Vector3d(geom.row(j))).norm(),
                                    m_snapshot_clash_ratio * 100.0, limit);
                                reported_clash = true;
                            }
                            collapsed = true;
                            break;
                        }
                    }
                }
                if (collapsed) {
                    rejected_clash++;
                    continue;
                }
            }
            const Matrix topo = TopologyMatrix(mol);
            if (topo.rows() != m_topo_matrix.rows()) {
                rejected_broken++;
                continue;
            }
            // Count the two failure modes separately: an EXTRA bond means two atoms collided (this
            // is what produces the NaN in the 1/r terms), a MISSING one means the molecule came
            // apart. Both disqualify the structure, but they say different things about the run.
            int formed = 0, broken_count = 0;
            int first_i = -1, first_j = -1;
            for (int i = 0; i < mol.AtomCount(); ++i) {
                for (int j = i + 1; j < mol.AtomCount(); ++j) {
                    const double d = topo(i, j) - m_topo_matrix(i, j);
                    if (std::abs(d) < 0.5)
                        continue;
                    if (d > 0)
                        formed++;
                    else
                        broken_count++;
                    if (first_i < 0) { first_i = i; first_j = j; }
                }
            }
            if (formed == 0 && broken_count == 0) {
                keep.push_back(std::move(mol));
                continue;
            }
            if (formed > 0)
                rejected_formed++;
            else
                rejected_broken++;
            // Only near-misses are worth repairing; a structure with many changed bonds is a
            // different molecule, not a conformer with a force-field artefact.
            if (m_repair_snapshots && formed + broken_count <= m_repair_max_bonds)
                broken.push_back(mol);
            if (!reported) {
                CurcumaLogger::warn_fmt("ConfSearch: topology gate (first reject): atoms {}-{}, {} bond(s) formed, "
                                        "{} broken (structure energy {:.6f} Eh)",
                    first_i, first_j, formed, broken_count, mol.Energy());
                reported = true;
            }
        }
    }

    const int rejected_total = rejected_formed + rejected_broken + rejected_clash + rejected_nonfinite;
    if (rejected_total == 0) {
        CurcumaLogger::result_fmt("ConfSearch: topology gate: all {} MD snapshots keep the reference topology "
                                  "and are free of collapsed contacts", total);
        return total;
    }

    // Claude Generated (Jul 2026): restrained repair instead of discarding (opt-in). Candidates are
    // taken lowest-energy first -- a structure 500 kJ/mol up will not become the conformer that
    // matters, and each repair costs two optimisations.
    int repaired = 0, repair_attempts = 0;
    if (m_repair_snapshots && !broken.empty()) {
        std::multimap<double, const Molecule*> by_energy; // the container sorts
        for (const Molecule& mol : broken)
            by_energy.insert(std::pair<double, const Molecule*>(mol.Energy(), &mol));

        // ONE calculator for the whole repair batch: identical topology/parameters for every
        // structure, and a second GFN-FF instance for the same molecule in one process has crashed
        // before (see ConfGen).
        nlohmann::json calc_cfg = ChildConfig(m_md_method, 1);
        calc_cfg["verbosity"] = 0; // otherwise every repair logs "Molecule set: N atoms" + timings
        const int saved_verbosity = CurcumaLogger::get_verbosity();
        CurcumaLogger::set_verbosity(0);
        EnergyCalculator calculator(m_md_method, calc_cfg);
        for (const auto& entry : by_energy) {
            if (repair_attempts >= m_repair_max)
                break;
            repair_attempts++;
            Molecule candidate = *entry.second;
            if (RepairSnapshot(candidate, calculator)) {
                keep.push_back(candidate);
                repaired++;
            }
        }
        CurcumaLogger::set_verbosity(saved_verbosity);
        CurcumaLogger::result_fmt("ConfSearch: topology repair: {} of {} attempted snapshot(s) pulled back to the "
                                  "reference topology (restrained, then released; {} candidate(s) had at most {} "
                                  "changed bond(s))",
            repaired, repair_attempts, static_cast<int>(broken.size()), m_repair_max_bonds);
    }

    // Rewrite the file with the survivors only.
    bool first = true;
    for (const Molecule& mol : keep) {
        if (first) { mol.writeXYZFile(path); first = false; }
        else          mol.appendXYZFile(path);
    }
    if (first)
        std::ofstream(path).close(); // nothing survived -> empty input, Phase 2 skips it

    // Report the chain, not two numbers that do not add up: N failed the check, R of them were
    // repaired and kept, N-R are gone.
    CurcumaLogger::result_fmt("ConfSearch: topology gate: {} of {} MD snapshots failed the check ({} formed a bond "
                              "(collision), {} broke one, {} collapsed contact, {} non-finite geometry){} -- {} "
                              "snapshots go into the optimisation",
        rejected_total, total, rejected_formed, rejected_broken, rejected_clash, rejected_nonfinite,
        repaired > 0 ? fmt::format(", {} of them repaired and kept", repaired)
                     : std::string(", none repaired"),
        static_cast<int>(keep.size()));
    return static_cast<int>(keep.size());
}

// Claude Generated (Jul 2026): copy every frame of an XYZ file (see the header for why a stage
// starts from a copy rather than chaining another suffix onto the previous stage's output).
int ConfSearch::CopyFrames(const std::string& source, const std::string& destination) const
{
    std::ifstream check(source);
    if (!check.good())
        return 0;
    int frames = 0;
    FileIterator it(source);
    while (!it.AtEnd()) {
        Molecule mol = it.Next();
        if (mol.AtomCount() == 0)
            continue;
        if (frames == 0)
            mol.writeXYZFile(destination);
        else
            mol.appendXYZFile(destination);
        frames++;
    }
    return frames;
}

// Claude Generated (Jul 2026): energies of an XYZ ensemble on disk, ascending. std::multiset does
// the ordering (the container sorts, no explicit sort call).
std::vector<double> ConfSearch::EnsembleEnergies(const std::string& path) const
{
    std::vector<double> energies;
    std::ifstream check(path);
    if (!check.good())
        return energies;
    std::multiset<double> ordered;
    FileIterator it(path);
    while (!it.AtEnd()) {
        Molecule mol = it.Next();
        if (mol.AtomCount() > 0)
            ordered.insert(mol.Energy());
    }
    energies.assign(ordered.begin(), ordered.end());
    return energies;
}

// Claude Generated (Jul 2026): compact energy summary of an ensemble file. A bare structure COUNT
// ("36 structures accepted") says nothing about what the phase produced -- this reports the lowest
// energy, the spread and the first few conformers, read back from the file that was just written.
void ConfSearch::ReportEnsemble(const std::string& label, const std::string& method, const std::string& path) const
{
    const std::vector<double> e = EnsembleEnergies(path);
    if (e.empty()) {
        CurcumaLogger::result_fmt("ConfSearch: {} [{}]: no structures", label, method);
        return;
    }
    const double e_min = e.front(), e_max = e.back();
    CurcumaLogger::result_fmt("ConfSearch: {} [{}]: {} structure(s), lowest {:.6f} Eh, span {:.2f} kJ/mol",
        label, method, static_cast<int>(e.size()), e_min, (e_max - e_min) * 2625.5);
    const int n_show = std::min(m_ensemble_report, static_cast<int>(e.size()));
    for (int i = 0; i < n_show; ++i)
        CurcumaLogger::result_fmt("ConfSearch:     #{:<3} {:.6f} Eh  (+{:.2f} kJ/mol)", i + 1, e[i], (e[i] - e_min) * 2625.5);
    if (n_show > 0 && static_cast<int>(e.size()) > n_show)
        CurcumaLogger::result_fmt("ConfSearch:     ... and {} more, up to +{:.2f} kJ/mol",
            static_cast<int>(e.size()) - n_show, (e_max - e_min) * 2625.5);
}

// Claude Generated (Jul 2026): per-cycle result files. The working files of a cycle
// (".bias.opt.accepted.xyz" and its opt_method counterpart) are overwritten by the next cycle, so
// what a given temperature actually produced was not recoverable after the run. This writes it out
// per cycle and per level of theory, energy-sorted, and keeps the most stable structure of every
// cycle in one trajectory so the progress of the search is a single file to look at.
void ConfSearch::WriteCycleEnsemble(const std::string& cycle_tag, const std::string& method,
    const std::vector<Molecule*>& molecules) const
{
    std::multimap<double, Molecule*> ordered; // the container sorts -- no std::sort anywhere
    for (Molecule* mol : molecules)
        if (mol && mol->AtomCount() > 0)
            ordered.insert(std::pair<double, Molecule*>(mol->Energy(), mol));
    if (ordered.empty())
        return;

    // "<base>.cycleNN_TxxxK.ensemble.<method>.xyz" -- same purpose+method convention as every other
    // file of the cycle (Claude Generated, Jul 2026).
    const std::string ensemble_file = outputPath(Basename() + "." + cycle_tag + ".ensemble." + method + ".xyz");
    bool first = true;
    for (const auto& entry : ordered) {
        if (first) { entry.second->writeXYZFile(ensemble_file); first = false; }
        else          entry.second->appendXYZFile(ensemble_file);
    }

    const std::string best_file = outputPath(Basename() + ".best_per_cycle." + method + ".xyz");
    ordered.begin()->second->appendXYZFile(best_file);

    const double e_min = ordered.begin()->first;
    const double e_max = ordered.rbegin()->first;
    CurcumaLogger::result_fmt("ConfSearch: {} ensemble [{}]: {} structure(s), most stable {:.6f} Eh, span {:.2f} kJ/mol -> {}",
        cycle_tag, method, static_cast<int>(ordered.size()), e_min, (e_max - e_min) * 2625.5, ensemble_file);
    CurcumaLogger::result_fmt("ConfSearch:   most stable structure of this cycle appended to {}", best_file);
}

int ConfSearch::PerformConfGen(const std::string& f, const std::string& method)
{
    const std::string input_file = outputPath(f + ".xyz");
    {
        std::ifstream check(input_file);
        if (!check.good()) {
            CurcumaLogger::warn("ConfSearch: Phase 3c skipped, no input file: " + input_file);
            return 0;
        }
    }

    // ConfGen's own defaults plus this phase's overrides. The couplings/model comparison are the
    // analysis half of ConfGen and are not needed here -- Phase 3c only wants the structures, and the
    // per-cycle ensembles are far too small for that statistic anyway.
    nlohmann::json cfg = ParameterRegistry::getInstance().getDefaultJson("confgen");
    cfg.merge_patch(ChildConfig(method, (m_threads > 1) ? 1 : m_threads));
    cfg["generate"] = true;
    cfg["couplings"] = false;
    cfg["max_proposals"] = m_confgen_max_proposals;
    cfg["proposal_templates"] = m_confgen_templates;
    cfg["proposal_depth"] = m_confgen_depth;
    cfg["new_rmsd"] = m_rmsd;   // "new" means the same thing here as everywhere else in the search
    cfg["verbosity"] = (m_verbosity >= 2) ? m_verbosity : 0;
    // Claude Generated (Jul 2026): the standalone -confgen path sets this (main.cpp), and GFN-FF keys
    // its topology/parameter cache on it. Without it the key is empty, so a cache written for ANOTHER
    // basename can be picked up -- which is how this phase first crashed inside the GFN-FF parameter
    // generation. Give the child its own, unambiguous cache identity.
    cfg["geometry_file"] = input_file;

    int added = 0;
    {
        ConfGen gen(cfg, false);
        gen.setOutputDir(OutputDir()); // same BMT directory as the rest of the run
        gen.setFile(input_file);
        gen.start();
    }

    // Append what survived ConfGen's own topology and novelty checks. From here on these structures
    // are indistinguishable from metadynamics hits: Phase 3b re-optimises them, Phase 4 topology- and
    // energy-filters them, feeds them into the bias pool and lets them compete as seeds.
    const std::string new_file = outputPath(f + ".proposals.new.xyz");
    std::ifstream check_new(new_file);
    if (!check_new.good())
        return 0;
    FileIterator it(new_file);
    while (!it.AtEnd()) {
        Molecule mol = it.Next();
        if (mol.AtomCount() == 0)
            continue;
        mol.setCharge(m_charge);
        mol.setSpin(m_spin);
        mol.appendXYZFile(input_file);
        added++;
    }
    return added;
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

// Claude Generated (Jul 2026): seed selection for the next temperature cycle.
//
// The old rule was purely energetic: sort the survivors of the seed energy window by energy and
// keep the seed_rank lowest ones. That is degenerate whenever the deepest minima are variations of
// the SAME basin -- three seeds 0.4 A apart start three MD runs that all explore the same region,
// while a genuinely different conformer 4 A away 2 kJ/mol higher is thrown away.
//
// "diverse" walks the same energy ranking (so the global minimum is always seed #1) but only
// accepts a candidate whose permutation-aware best-fit RMSD to EVERY already chosen seed is at
// least r_min. If that leaves fewer seeds than requested, the spacing is relaxed (r_min/2, r_min/4)
// and finally dropped altogether, so the MD phase is never starved -- the fallback is exactly the
// old energy ranking. Cost: O(kept * candidates) Kabsch fits on a handful of structures, i.e.
// nothing next to one MD step.
int ConfSearch::SelectSeeds(std::vector<Molecule*>& window_seeds, const std::string& method,
    double global_min) const
{
    const int n = static_cast<int>(window_seeds.size());
    if (n <= 1)
        return 0;

    // Energy ranking is the backbone of both strategies.
    std::sort(window_seeds.begin(), window_seeds.end(),
        [](const Molecule* a, const Molecule* b) { return a->Energy() < b->Energy(); });
    const double e_ref = window_seeds[0]->Energy();
    const int limit = (m_seed_rank > 0) ? std::min(m_seed_rank, n) : n;

    std::vector<Molecule*> keep;
    if (m_seed_selection != "diverse") {
        keep.assign(window_seeds.begin(), window_seeds.begin() + limit);
        CurcumaLogger::result_fmt("ConfSearch: seed selection (energy) on the {} PES: {} of {} candidate(s) kept "
                                  "({} rejected by rank {})",
            method, limit, n, n - limit, m_seed_rank);
    } else {
        const double r_min = (m_seed_min_rmsd > 0.0) ? m_seed_min_rmsd
                                                     : m_seed_diversity_factor * m_rmsd;
        std::vector<char> taken(n, 0);
        keep.push_back(window_seeds[0]); // the most stable structure is always a seed
        taken[0] = 1;
        // Distance of each accepted seed to the seeds accepted before it (diagnostics).
        std::vector<double> spacing(1, 0.0);
        double r = r_min;
        bool relaxed = false;
        while (static_cast<int>(keep.size()) < limit && r > 1e-3) {
            for (int i = 1; i < n && static_cast<int>(keep.size()) < limit; ++i) {
                if (taken[i])
                    continue;
                double dmin = std::numeric_limits<double>::infinity();
                for (const Molecule* s : keep)
                    dmin = std::min(dmin, PermRMSD(s->getGeometry(), window_seeds[i]->getGeometry()));
                if (dmin >= r) {
                    keep.push_back(window_seeds[i]);
                    spacing.push_back(dmin);
                    taken[i] = 1;
                }
            }
            if (static_cast<int>(keep.size()) < limit) {
                r *= 0.5; // not enough distinct basins available -> relax the spacing
                relaxed = true;
            }
        }
        CurcumaLogger::result_fmt("ConfSearch: seed selection (diverse) on the {} PES: {} of {} candidate(s) kept, "
                                  "spacing >= {:.2f} A{}",
            method, static_cast<int>(keep.size()), n, r_min,
            relaxed ? fmt::format(" (relaxed to {:.2f} A -- not enough distinct basins)", r) : "");
        m_seed_spacing = spacing;
    }

    // Claude Generated (Jul 2026): say exactly WHICH structures the next cycle starts from -- name,
    // absolute energy on the PES that picked them, distance to this cycle's best and to the running
    // global minimum, and how far each seed sits from the ones chosen before it. "seeded from the
    // opt PES" is otherwise an unverifiable claim.
    CurcumaLogger::result_fmt("ConfSearch: the next cycle starts from these {} structure(s) [{}]:",
        static_cast<int>(keep.size()), method);
    for (int i = 0; i < static_cast<int>(keep.size()); ++i) {
        const double e = keep[i]->Energy();
        std::string spacing_str = "--";
        if (i > 0) {
            if (m_seed_selection == "diverse" && i < static_cast<int>(m_seed_spacing.size()))
                spacing_str = fmt::format("{:.2f} A", m_seed_spacing[i]);
            else
                spacing_str = fmt::format("{:.2f} A", PermRMSD(keep[0]->getGeometry(), keep[i]->getGeometry()));
        }
        CurcumaLogger::result_fmt("ConfSearch:   seed {:>2}: {:<16} E = {:.6f} Eh | {:+7.2f} kJ/mol vs this cycle's "
                                  "best | {:+7.2f} vs the running best | RMSD to previous seeds: {}",
            i + 1, keep[i]->Name().empty() ? std::string("(unnamed)") : keep[i]->Name(), e,
            (e - e_ref) * 2625.5,
            std::isfinite(global_min) ? (e - global_min) * 2625.5 : 0.0, spacing_str);
    }

    // Everything not kept is dropped (and owned here, so deleted).
    int rejected = 0;
    for (Molecule* mol : window_seeds) {
        if (std::find(keep.begin(), keep.end(), mol) == keep.end()) {
            delete mol;
            rejected++;
        }
    }
    window_seeds = std::move(keep);
    return rejected;
}

void ConfSearch::CalibrateBias(const std::string& p, nlohmann::json& md)
{
    auto load = [](const std::string& file) {
        std::vector<Molecule> v;
        FileIterator f(file);
        while (!f.AtEnd()) { Molecule m = f.Next(); if (m.AtomCount() > 0) v.push_back(m); }
        return v;
    };
    const std::string explore = cycleStage("explore", m_md_method);
    std::vector<Molecule> pre = load(outputPath(explore + ".xyz"));            // pre-optimisation MD snapshots
    std::vector<Molecule> post = load(outputPath(explore + ".opt.xyz"));       // optimised (index-aligned with pre)
    std::vector<Molecule> minima = load(outputPath(explore + ".opt.accepted.xyz")); // distinct minima (deduped)
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
    // Claude Generated (Jul 2026): the opt_method reference energy of the input structure. Without it
    // a resumed dual-method run had no opt-PES anchor left (m_initial_energy_opt = inf after the
    // resume skips the initial optimisation), so the final "search lowered the energy by X kJ/mol
    // (opt_method: ...)" line -- the only one on the ranking PES -- silently disappeared.
    state["initial_energy_opt"] = m_initial_energy_opt;

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
    // Claude Generated (Jul 2026): the opt-PES topology reference. Without it a resumed dual-method
    // run would fall back to comparing opt_method geometries against the md_method reference again.
    if (m_topo_ref_opt.AtomCount() > 0)
        state["topo_ref_opt"] = molToJson(m_topo_ref_opt);
    state["permutations"] = m_permutation_cache;

    // Intermediate accepted sets, only needed to resume mid-cycle without redoing the opts.
    if (phase == "post_filter" || phase == "post_refine")
        state["accepted_md"] = fileFramesToJson(outputPath(cycleStage("explore", m_md_method) + ".opt.accepted.xyz"));
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
    // Claude Generated (Jul 2026): energies that are still infinity when the checkpoint is written
    // (e.g. initial_energy_opt in a single-method run) are serialised by nlohmann as JSON *null*,
    // and `value<double>()` THROWS on null instead of returning the fallback. That threw out of
    // loadCheckpoint, uncaught -> abort: every `-restart` resume of a single-method run died with a
    // stack trace in get_arithmetic_value. Read them null-tolerantly.
    auto read_energy = [&s](const char* key) {
        if (!s.contains(key) || s[key].is_null() || !s[key].is_number())
            return std::numeric_limits<double>::infinity();
        return s[key].get<double>();
    };
    st.global_min = read_energy("global_min");
    st.best_energy = read_energy("best_energy");
    st.initial_energy = read_energy("initial_energy");
    // Missing in checkpoints written before Jul 2026 -> stays infinity, i.e. the old behaviour.
    st.initial_energy_opt = read_energy("initial_energy_opt");

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
            bs.counter = b.value("counter", 0.0);
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
    if (s.contains("topo_ref_opt"))
        st.topo_ref_opt = jsonToMol(st.elements, s["topo_ref_opt"]);
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
    m_snapshot_topology_gate = m_config.get<bool>("snapshot_topology_gate"); // Claude Generated (Jul 2026)
    m_snapshot_clash_ratio = m_config.get<double>("snapshot_clash_ratio");
    m_topology_factor = m_config.get<double>("topology_factor");
    m_repair_snapshots = m_config.get<bool>("repair_snapshots");
    m_repair_max = m_config.get<int>("repair_max");
    m_repair_max_bonds = m_config.get<int>("repair_max_bonds");
    m_repair_force = m_config.get<double>("repair_force");
    m_repair_max_iterations = m_config.get<int>("repair_max_iterations");
    m_topo_check = m_config.get<bool>("topo_check");
    m_topo_check_interval = m_config.get<int>("topo_check_interval");
    m_seed_energy_window = m_config.get<double>("seed_energy_window");
    m_seed_rank = m_config.get<int>("seed_rank");
    // Claude Generated (Jul 2026): RMSD-aware seed selection
    m_seed_selection = m_config.get<std::string>("seed_selection");
    m_seed_pes = m_config.get<std::string>("seed_pes");
    if (m_seed_pes != "md" && m_seed_pes != "opt" && m_seed_pes != "both") {
        CurcumaLogger::warn_fmt("ConfSearch: seed_pes='{}' unknown -- falling back to 'md'", m_seed_pes);
        m_seed_pes = "md";
    }
    m_seed_min_rmsd = m_config.get<double>("seed_min_rmsd");
    m_seed_diversity_factor = m_config.get<double>("seed_diversity_factor");
    if (m_seed_selection != "diverse" && m_seed_selection != "energy") {
        CurcumaLogger::warn_fmt("ConfSearch: seed_selection='{}' unknown -- falling back to 'energy'", m_seed_selection);
        m_seed_selection = "energy";
    }
    m_seed_window_schedule = m_config.get<std::string>("seed_window_schedule");
    m_seed_window_decay = m_config.get<double>("seed_window_decay");
    m_epot_abort = m_config.get<bool>("epot_abort");
    m_epot_abort_window = m_config.get<double>("epot_abort_window");
    m_temp_abort = m_config.get<bool>("temp_abort");
    m_temp_abort_factor = m_config.get<double>("temp_abort_factor");
    m_temp_abort_delta = m_config.get<double>("temp_abort_delta");
    m_rmsd_mtd_max_height = m_config.get<int>("rmsd_mtd_max_height");
    m_freeze_inherited = m_config.get<bool>("rmsd_mtd_freeze_inherited");
    m_rmsd_mtd_max_gaussians = m_config.get<int>("rmsd_mtd_max_gaussians");  // Claude Generated (Jul 2026)
    m_rmsd_mtd_screen = m_config.get<bool>("rmsd_mtd_screen");
    m_rmsd_mtd_cutoff_tol = m_config.get<double>("rmsd_mtd_cutoff_tol");
    m_rmsd_mtd_screen_margin = m_config.get<double>("rmsd_mtd_screen_margin");
    m_opt_feedback_bias = m_config.get<bool>("opt_feedback_bias");
    m_opt_feedback_height = m_config.get<int>("opt_feedback_height");
    m_opt_feedback_prune_snapshots = m_config.get<bool>("opt_feedback_prune_snapshots");
    m_mtd_permutation = m_config.get<bool>("mtd_permutation");
    m_bias_calibration = m_config.get<std::string>("bias_calibration");
    m_bias_couple_factor = m_config.get<double>("bias_couple_factor");
    m_bias_scale_mode = m_config.get<std::string>("bias_scale_mode");
    m_bias_energy_tol = m_config.get<double>("bias_energy_tol");
    m_do_restart = m_config.get<bool>("restart");
    // Claude Generated (Jul 2026): torsion-recombination phase
    m_confgen_phase = m_config.get<bool>("confgen_phase");
    m_confgen_max_proposals = m_config.get<int>("confgen_max_proposals");
    m_confgen_templates = m_config.get<int>("confgen_templates");
    m_confgen_depth = m_config.get<int>("confgen_depth");
    // Claude Generated (Jul 2026): two-stage Phase 3b + per-cycle reporting
    m_phase3b_two_stage = m_config.get<bool>("phase3b_two_stage");
    m_phase3b_preopt_preset = m_config.get<std::string>("phase3b_preopt_preset");
    m_phase3b_preopt_max_iter = m_config.get<int>("phase3b_preopt_max_iter");
    m_phase3b_preset = m_config.get<std::string>("phase3b_preset");
    m_phase3b_filter = m_config.get<bool>("phase3b_filter");
    m_cycle_output = m_config.get<bool>("cycle_output");
    m_ensemble_report = m_config.get<int>("ensemble_report");
    static const std::vector<std::string> presets = { "loose", "normal", "tight", "verytight" };
    for (std::string* p : { &m_phase3b_preopt_preset, &m_phase3b_preset }) {
        if (std::find(presets.begin(), presets.end(), *p) == presets.end()) {
            CurcumaLogger::warn_fmt("ConfSearch: convergence preset '{}' unknown -- falling back to 'normal'", *p);
            *p = "normal";
        }
    }
}
