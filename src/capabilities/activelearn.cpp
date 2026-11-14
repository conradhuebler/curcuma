/*
 * < Delta-ML training and active learning driver >
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
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
 * Claude Generated (July 2026): driver of the Delta-ML work package.
 */

#include "activelearn.h"

#include "src/core/citation_registry.h"
#include "src/core/curcuma_logger.h"
#include "src/core/energycalculator.h"
#include "src/core/fileiterator.h"
#include "src/core/intra_parallel_context.h"
#include "src/core/parameter_registry.h"
#include "src/core/units.h"

#include "simplemd.h"

#include "src/core/energy_calculators/ff_methods/gfnff.h"

#include "external/CxxThreadPool/include/CxxThreadPool.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>

using curcuma::ml::DeltaMLErrors;
using curcuma::ml::DeltaMLModel;
using curcuma::ml::DeltaMLTrainingOptions;
using curcuma::ml::DeltaMLTrainingReport;
using curcuma::ml::MLDataset;
using curcuma::ml::MLStructure;

namespace {

/**
 * @brief Computes base and reference energy plus gradient for a batch of structures.
 *
 * Every structure gets freshly constructed calculators, exactly as `-sp` and the deployed
 * model (DeltaMLMethod::initialiseBase, MethodFactory::create per setMolecule) do. Only the
 * anchor bond graph is carried over, via Mol::m_bonds, which DeltaMLMethod applies in the
 * same way; so training and deployment see the same base energy for the same geometry.
 *
 * Results are kept with their input index, so the caller can restore the input order.
 */
class LabelThread : public CxxThread {
public:
    LabelThread(const std::string& base_method, const std::string& ref_method, const json& controller,
        const std::vector<std::pair<int, int>>& anchor_bonds)
        : m_base_method(base_method)
        , m_ref_method(ref_method)
        , m_controller(controller)
        , m_anchor_bonds(anchor_bonds)
    {
        setAutoDelete(false);
    }

    void addStructure(const Molecule& m, std::size_t index)
    {
        m_molecules.push_back(m);
        m_indices.push_back(index);
    }

    int execute() override
    {
        if (m_molecules.empty())
            return 0;

        // Molecule-level batch worker: every method it calls stays serial (same guard as the
        // -sp/-opt batch workers), so a label does not depend on how many workers run.
        curcuma::SuppressIntraParallel intra_guard;

        json cfg = m_controller;
        cfg["threads"] = 1;
        // Thousands of per-structure SCF banners would bury the actual result.
        cfg["verbosity"] = 0;
        // ConfigManager resolves controller["global"] BEFORE flat top-level keys, and the CLI
        // puts "-threads N" there too, so the line above alone left every GFN-FF instance of a
        // -threads 4 run on 4 threads of its own (labels then differed from -threads 1 by
        // 1.8e-14 Eh, Oct 2026). Override both places.
        if (cfg.contains("global") && cfg["global"].is_object()) {
            cfg["global"]["threads"] = 1;
            if (cfg["global"].contains("verbosity"))
                cfg["global"]["verbosity"] = 0;
        }
        // Claude Generated (Oct 2026): a topology file cache (<basename>.topo.json) would hand
        // every structure the topology of whichever geometry wrote it first. The label must
        // come from the structure's own geometry (plus the pinned graph), as in `-sp`.
        json gfnff_scope = cfg.value("gfnff", json::object());
        gfnff_scope["cache_topology"] = false;
        cfg["gfnff"] = gfnff_scope;

        // External programs that write scratch files under a fixed name cannot share a
        // working directory. ORCA in particular uses one 'orca_calc' basename, so parallel
        // labelling would have every thread overwriting the same input and reading back
        // whichever output happened to survive. Give each thread its own.
        const std::string tag = "orca_calc_t" + std::to_string(getThreadId());
        cfg["orca_basename"] = tag;
        json orca = cfg.value("orca", json::object());
        orca["orca_basename"] = tag;
        cfg["orca"] = orca;

        try {
            for (std::size_t k = 0; k < m_molecules.size(); ++k) {
                const Molecule& mol = m_molecules[k];

                MLStructure s;
                s.elements = mol.Atoms();
                s.geometry = mol.getGeometry();
                s.charge = mol.Charge();
                s.spin = mol.Spin();

                Mol pinned = mol.getMolInfo();
                if (!m_anchor_bonds.empty())
                    pinned.m_bonds = m_anchor_bonds;

                // Claude Generated (Oct 2026): new calculators per structure. One calculator per
                // thread, re-set with setMolecule() for each structure, stored base energies that
                // differed from `-sp` by -5.1 +- 2.0 kcal/mol (trisaccharide MD set, n=401) and were
                // exact only for the first structure of each thread. Cause: a re-initialised GFNFF
                // keeps EEQSolver::m_external_srab pointing at the distance array of the previous
                // structure, so the Phase-2 EEQ solve during topology setup uses stale distances,
                // and Calculation() skips its own Phase-2 solve because the geometry "has not
                // changed" since setup (docs/DELTA_ML.md, "Open: the anchor freezes more than it
                // should"). The reference is treated the same way, so a warm start can never make
                // a label depend on the order or the thread chunking.
                EnergyCalculator base(m_base_method, cfg);
                base.setMolecule(pinned);
                s.energy_base = base.CalculateEnergy(true);
                if (base.Error())
                    throw std::runtime_error("base method failed on structure "
                        + std::to_string(m_indices[k] + 1) + ": " + base.ErrorMessage());
                // EnergyCalculator::Gradient() is Eh/Angstrom (getGradient() contract, Known
                // Issue #28), which is also what the descriptor and the whole model work in.
                // An earlier version assumed Eh/Bohr and multiplied by ANGSTROM_TO_BOHR; after the
                // merge with master that made every gradient label 1.89 times too large
                // (dataset gradient_base vs `-sp -dump_gradient`: ratio 1.8906, n=1 structure),
                // so the energy rows and the gradient rows of the fit contradicted each other.
                s.gradient_base = base.Gradient();

                EnergyCalculator ref(m_ref_method, cfg);
                ref.setMolecule(pinned);
                s.energy_ref = ref.CalculateEnergy(true);
                if (ref.Error())
                    throw std::runtime_error("reference method failed on structure "
                        + std::to_string(m_indices[k] + 1) + ": " + ref.ErrorMessage());
                s.gradient_ref = ref.Gradient();
                s.has_gradients = true;

                m_results.emplace_back(m_indices[k], s);
            }
        } catch (const std::exception& e) {
            m_error = e.what();
            return 1;
        }
        return 0;
    }

    /// Labelled structures with their index in the input list.
    const std::vector<std::pair<std::size_t, MLStructure>>& results() const { return m_results; }
    const std::string& error() const { return m_error; }

private:
    std::string m_base_method;
    std::string m_ref_method;
    json m_controller;
    std::vector<std::pair<int, int>> m_anchor_bonds;
    std::vector<Molecule> m_molecules;
    std::vector<std::size_t> m_indices;
    std::vector<std::pair<std::size_t, MLStructure>> m_results;
    std::string m_error;
};

std::vector<int> parseIntList(const std::string& text)
{
    std::vector<int> out;
    std::stringstream ss(text);
    std::string item;
    while (std::getline(ss, item, ',')) {
        try {
            const int v = std::stoi(item);
            if (v > 0)
                out.push_back(v);
        } catch (const std::exception&) {
            // ignore malformed entries
        }
    }
    return out;
}

} // namespace

ActiveLearn::ActiveLearn(const json& controller, bool silent)
    : CurcumaMethod(ParameterRegistry::getInstance().getDefaultJson("activelearn"), controller, silent)
{
    // UpdateController merges the 'activelearn' sub-document into the defaults and calls
    // LoadControlJson() itself.
    UpdateController(controller);
}

void ActiveLearn::LoadControlJson()
{
    m_mode = Json2KeyWord<std::string>(m_defaults, "al_mode");
    m_ref_method = Json2KeyWord<std::string>(m_defaults, "ref_method");
    m_model_path = Json2KeyWord<std::string>(m_defaults, "model");
    m_max_structures = Json2KeyWord<int>(m_defaults, "max_structures");
    m_test_fraction = Json2KeyWord<double>(m_defaults, "test_fraction");
    m_n_centers = Json2KeyWord<int>(m_defaults, "n_centers");
    m_kernel_sigma = Json2KeyWord<double>(m_defaults, "kernel_sigma");
    m_lambda = Json2KeyWord<double>(m_defaults, "lambda");
    m_jitter_relative = Json2KeyWord<double>(m_defaults, "jitter_relative");
    m_sigma_energy = Json2KeyWord<double>(m_defaults, "sigma_energy");
    m_sigma_force = Json2KeyWord<double>(m_defaults, "sigma_force");
    m_use_forces = Json2KeyWord<bool>(m_defaults, "use_forces");
    m_seed = static_cast<unsigned>(Json2KeyWord<int>(m_defaults, "seed"));
    m_write_dataset = Json2KeyWord<bool>(m_defaults, "write_dataset");
    m_al_cycles = Json2KeyWord<int>(m_defaults, "al_cycles");
    m_al_md_method = Json2KeyWord<std::string>(m_defaults, "al_md_method");
    m_al_temperature = Json2KeyWord<double>(m_defaults, "al_temperature");
    m_al_md_time = Json2KeyWord<double>(m_defaults, "al_md_time");
    m_al_md_runs = Json2KeyWord<int>(m_defaults, "al_md_runs");
    m_al_dump_frequency = Json2KeyWord<int>(m_defaults, "al_dump_frequency");
    m_al_candidates = Json2KeyWord<int>(m_defaults, "al_candidates");
    m_al_uncertainty = Json2KeyWord<double>(m_defaults, "al_uncertainty");
    m_al_initial = Json2KeyWord<int>(m_defaults, "al_initial");
    m_al_onthefly = Json2KeyWord<bool>(m_defaults, "al_onthefly");
    m_al_onthefly_threshold = Json2KeyWord<double>(m_defaults, "al_onthefly_threshold");
    // Prospective validation and convergence (Claude Generated, October 2026).
    m_al_validation_frames = std::max(0, Json2KeyWord<int>(m_defaults, "al_validation_frames"));
    m_al_criterion = Json2KeyWord<std::string>(m_defaults, "al_criterion");
    std::transform(m_al_criterion.begin(), m_al_criterion.end(), m_al_criterion.begin(),
        [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    m_al_target_energy = Json2KeyWord<double>(m_defaults, "al_target_energy");
    m_al_target_hole = Json2KeyWord<double>(m_defaults, "al_target_hole");
    m_al_patience = std::max(1, Json2KeyWord<int>(m_defaults, "al_patience"));
    m_al_window = std::max(1, Json2KeyWord<int>(m_defaults, "al_window"));
    m_al_hole_quantile = std::clamp(Json2KeyWord<double>(m_defaults, "al_hole_quantile"), 0.5, 1.0);
    m_center_sweep = parseIntList(Json2KeyWord<std::string>(m_defaults, "center_sweep"));
    m_dataset_path = Json2KeyWord<std::string>(m_defaults, "dataset");

    // Energy-only second model (Claude Generated, October 2026).
    m_energy_model = Json2KeyWord<bool>(m_defaults, "energy_model");
    m_energy_model_n_centers = Json2KeyWord<int>(m_defaults, "energy_model_n_centers");
    m_energy_model_lambda = Json2KeyWord<double>(m_defaults, "energy_model_lambda");
    m_energy_model_kernel_sigma = Json2KeyWord<double>(m_defaults, "energy_model_kernel_sigma");

    // Descriptor grid (Claude Generated, October 2026). Defaults equal the built-in grid.
    m_acsf.n_shell = Json2KeyWord<int>(m_defaults, "acsf_n_shell");
    m_acsf.r_first = Json2KeyWord<double>(m_defaults, "acsf_r_first");
    m_acsf.r_last = Json2KeyWord<double>(m_defaults, "acsf_r_last");
    m_acsf.radial_eta = Json2KeyWord<double>(m_defaults, "acsf_radial_eta");
    m_acsf.rcut_radial = Json2KeyWord<double>(m_defaults, "acsf_rcut_radial");
    m_acsf.rcut_angular = Json2KeyWord<double>(m_defaults, "acsf_rcut_angular");
    m_acsf.angular_eta = Json2KeyWord<double>(m_defaults, "acsf_angular_eta");
    // A single value ("-acsf_zeta 4") may arrive as a JSON number rather than a string.
    const json zeta_raw = Json2KeyWord<json>(m_defaults, "acsf_zeta");
    std::string zeta_text;
    if (zeta_raw.is_string())
        zeta_text = zeta_raw.get<std::string>();
    else if (zeta_raw.is_array())
        for (const auto& z : zeta_raw)
            zeta_text += (zeta_text.empty() ? "" : ",") + z.dump();
    else
        zeta_text = zeta_raw.dump();
    m_acsf_valid = curcuma::ml::ACSFSettings::parseZetaList(zeta_text, m_acsf.angular_zeta);
    if (!m_acsf_valid)
        CurcumaLogger::error_fmt(
            "activelearn: acsf_zeta '{}' is not a comma list of positive numbers", zeta_text);
    if (m_acsf.n_shell < 1) {
        CurcumaLogger::error_fmt("activelearn: acsf_n_shell must be >= 1 (got {})", m_acsf.n_shell);
        m_acsf_valid = false;
    }

    // The base method is the global -method setting; that keeps the CLI consistent with
    // every other capability.
    if (m_controller.contains("method") && m_controller["method"].is_string())
        m_base_method = m_controller["method"].get<std::string>();
    if (m_controller.contains("threads"))
        m_threads = m_controller["threads"].get<int>();
    if (m_threads < 1)
        m_threads = 1;
}

bool ActiveLearn::loadStructures(std::vector<Molecule>& molecules)
{
    molecules.clear();
    FileIterator iter(Filename(), true);
    while (!iter.AtEnd()) {
        molecules.push_back(iter.Next());
        if (m_max_structures > 0 && static_cast<int>(molecules.size()) >= m_max_structures)
            break;
    }

    if (molecules.empty()) {
        CurcumaLogger::error_fmt("activelearn: no structures read from '{}'", Filename());
        return false;
    }

    // Every structure must describe the same molecule; a delta model with a frozen
    // topology is defined for one composition only.
    const int n0 = molecules[0].AtomCount();
    for (std::size_t i = 1; i < molecules.size(); ++i) {
        if (molecules[i].AtomCount() != n0) {
            CurcumaLogger::error_fmt(
                "activelearn: structure {} has {} atoms but the first has {}. This capability "
                "trains a model for a single molecule at fixed topology.",
                i + 1, molecules[i].AtomCount(), n0);
            return false;
        }
    }

    m_anchor = molecules[0];
    captureAnchorBonds();
    return true;
}

void ActiveLearn::captureAnchorBonds()
{
    m_anchor_bonds.clear();
    // Only GFN-FF-style methods carry a frozen graph worth pinning; for anything else the
    // empty list simply means "detect as usual".
    if (m_base_method.rfind("gfnff", 0) != 0)
        return;

    const int saved = CurcumaLogger::get_verbosity();
    CurcumaLogger::set_verbosity(0);
    GFNFF probe(m_controller);
    const bool ok = probe.InitialiseMolecule(m_anchor.getMolInfo());
    if (ok)
        m_anchor_bonds = probe.getCachedBondList();
    CurcumaLogger::set_verbosity(saved);

    if (m_anchor_bonds.empty())
        CurcumaLogger::warn("activelearn: could not capture a bond graph from the anchor; "
                            "topology will be redetected for every structure");
    else
        CurcumaLogger::result_fmt(
            "activelearn: anchor bond graph captured ({} bonds); it is pinned for every "
            "structure while coordination numbers and parameters follow each geometry",
            m_anchor_bonds.size());
}

bool ActiveLearn::labelStructures(const std::vector<Molecule>& molecules, MLDataset& data)
{
    const auto t0 = std::chrono::high_resolution_clock::now();
    CurcumaLogger::result_fmt("activelearn: labelling {} structures  base = {}  reference = {}",
        molecules.size(), m_base_method, m_ref_method);

    // One worker per thread, each with a round-robin share of the structures. Every structure
    // still gets its own calculators (see LabelThread); the chunking only saves pool overhead.
    const int n_threads = std::max(1, m_threads);
    std::vector<LabelThread*> threads;
    threads.reserve(n_threads);
    for (int t = 0; t < n_threads; ++t) {
        LabelThread* lt
            = new LabelThread(m_base_method, m_ref_method, m_controller, m_anchor_bonds);
        lt->setThreadId(t);
        threads.push_back(lt);
    }
    for (std::size_t i = 0; i < molecules.size(); ++i)
        threads[i % n_threads]->addStructure(molecules[i], i);

    CxxThreadPool* pool = new CxxThreadPool;
    pool->setActiveThreadCount(n_threads);
    for (auto* t : threads)
        pool->addThread(t);

    // The logger verbosity is a shared static, so a per-thread setting is not available
    // (see the verbosity-scoping note in the top-level CLAUDE.md). Silence it around the
    // whole batch and restore afterwards.
    const int saved_verbosity = CurcumaLogger::get_verbosity();
    CurcumaLogger::set_verbosity(0);
    pool->StartAndWait();
    CurcumaLogger::set_verbosity(saved_verbosity);

    // Structures were dealt out round-robin (structure i to thread i % n_threads), so the
    // per-thread result lists interleave. Restore the input order: the k-th dataset entry is
    // the k-th input structure, independent of -threads (Claude Generated, Oct 2026).
    bool ok = true;
    std::vector<std::pair<std::size_t, const MLStructure*>> labelled;
    labelled.reserve(molecules.size());
    for (auto* t : threads) {
        if (!t->error().empty()) {
            CurcumaLogger::error_fmt("activelearn: labelling failed: {}", t->error());
            ok = false;
        }
        for (const auto& r : t->results())
            labelled.emplace_back(r.first, &r.second);
    }
    std::sort(labelled.begin(), labelled.end(),
        [](const auto& a, const auto& b) { return a.first < b.first; });
    for (const auto& r : labelled)
        data.add(*r.second);

    // The pool destructor still inspects every thread it holds, so it has to go first;
    // deleting the threads before it is a use-after-free.
    delete pool;
    for (auto* t : threads)
        delete t;

    if (!ok || data.empty())
        return false;

    const auto t1 = std::chrono::high_resolution_clock::now();
    const double seconds = std::chrono::duration<double>(t1 - t0).count();

    // The spread of the learning target, so the errors below have a scale to be judged against.
    double dmin = data[0].delta(), dmax = data[0].delta(), dsum = 0.0;
    for (std::size_t i = 0; i < data.size(); ++i) {
        dmin = std::min(dmin, data[i].delta());
        dmax = std::max(dmax, data[i].delta());
        dsum += data[i].delta();
    }
    CurcumaLogger::result_fmt("activelearn: labelled {} structures in {:.1f} s ({:.2f} s each)",
        data.size(), seconds, seconds / static_cast<double>(data.size()));
    CurcumaLogger::result_fmt(
        "activelearn: delta = E_ref - E_base spans {:.6f} Eh ({:.2f} kcal/mol), mean {:.6f} Eh",
        dmax - dmin, (dmax - dmin) * 627.5094740631, dsum / static_cast<double>(data.size()));

    return true;
}

void ActiveLearn::logErrors(const std::string& label, const DeltaMLErrors& e) const
{
    const double kcal = 627.5094740631;
    CurcumaLogger::result_fmt("  {}  ({} structures)", label, e.n_structures);
    CurcumaLogger::result_fmt("    relative energy MAE : model {:9.4f} kcal/mol   base {:9.4f} kcal/mol",
        e.energy_mae_centered * kcal, e.base_energy_mae_centered * kcal);
    if (e.force_rmse > 0.0)
        CurcumaLogger::result_fmt("    gradient RMSE       : model {:9.6f} Eh/A       base {:9.6f} Eh/A",
            e.force_rmse, e.base_force_rmse);
    CurcumaLogger::result_fmt("    rank correlation    : model {:9.4f}            base {:9.4f}",
        e.spearman, e.base_spearman);
}

DeltaMLTrainingOptions ActiveLearn::trainingOptions(int n_centers) const
{
    DeltaMLTrainingOptions opt;
    opt.n_centers = n_centers;
    opt.kernel_sigma = m_kernel_sigma;
    opt.lambda = m_lambda;
    opt.jitter_relative = m_jitter_relative;
    opt.sigma_energy = m_sigma_energy;
    opt.sigma_force = m_sigma_force;
    opt.use_forces = m_use_forces;
    opt.seed = m_seed;
    opt.verbosity = m_verbosity;
    opt.acsf = m_acsf; // same descriptor grid for every fit, in every mode and cycle
    // The global -gpu setting drives the fit as well, not just the energy methods.
    if (m_controller.contains("gpu") && m_controller["gpu"].is_string())
        opt.gpu = m_controller["gpu"].get<std::string>();
    return opt;
}

std::string ActiveLearn::energyModelPath(const std::string& main_model_path)
{
    // <x>.json -> <x>.energy.json, so the default <basename>.deltaml.json gives
    // <basename>.deltaml.energy.json in the same directory.
    const std::string ext = ".json";
    if (main_model_path.size() > ext.size()
        && main_model_path.compare(main_model_path.size() - ext.size(), ext.size(), ext) == 0)
        return main_model_path.substr(0, main_model_path.size() - ext.size()) + ".energy.json";
    return main_model_path + ".energy.json";
}

bool ActiveLearn::fitAndReport(const MLDataset& train, const MLDataset& test, int n_centers,
    DeltaMLModel& model, bool log_details)
{
    return fitAndReport(train, test, trainingOptions(n_centers), model, log_details, "", true);
}

bool ActiveLearn::fitAndReport(const MLDataset& train, const MLDataset& test,
    const DeltaMLTrainingOptions& opt, DeltaMLModel& model, bool log_details,
    const std::string& label, bool log_heldout)
{
    // Fitting a model is as much a use of these methods as evaluating one.
    CitationRegistry::cite("deltaml");
    CitationRegistry::cite("acsf", "deltaml");
    CitationRegistry::cite("gap_sparse", "deltaml");

    DeltaMLTrainingReport report;
    if (!model.train(train, opt, report))
        return false;

    model.setBaseMethod(m_base_method);
    model.setReferenceMethod(m_ref_method);
    model.setAnchor(m_anchor.Atoms(), m_anchor.getGeometry(), m_anchor_bonds);
    model.setRole(opt.use_forces ? "joint" : "energy");

    const std::string suffix = label.empty() ? std::string() : " (" + label + ")";
    if (log_details) {
        CurcumaLogger::result_fmt(
            "activelearn{}: fitted {} centres from {} environments, {} features, sigma = {:.4f}",
            suffix, report.n_centers, report.n_environments, report.n_features, report.kernel_sigma);
        CurcumaLogger::result_fmt("activelearn: {} energy rows, {} gradient rows, solved in {:.1f} s",
            report.n_energy_rows, report.n_force_rows, report.train_seconds);
        if (report.jitter_used > 0.0)
            CurcumaLogger::warn_fmt(
                "activelearn: the normal equations needed a jitter of {:.1e} to factorise; "
                "consider a larger lambda or fewer centres",
                report.jitter_used);
        logErrors("training set" + suffix, report.train_errors);
    }

    if (log_heldout && !test.empty())
        logErrors(log_details ? "held-out set" + suffix : fmt::format("M = {}", report.n_centers),
            model.evaluate(test));

    return true;
}

void ActiveLearn::runTrain()
{
    std::vector<Molecule> molecules;
    if (!loadStructures(molecules))
        return;

    MLDataset data;
    if (!m_dataset_path.empty()) {
        // Reuse labels from an earlier run (Claude Generated, October 2026). The input
        // structures only provide the basename and the anchor; they are not labelled.
        const auto t0 = std::chrono::high_resolution_clock::now();
        if (!std::filesystem::exists(m_dataset_path)) {
            CurcumaLogger::error_fmt("activelearn: dataset file '{}' does not exist", m_dataset_path);
            return;
        }
        if (!data.load(m_dataset_path))
            return;
        if (data.empty()) {
            CurcumaLogger::error_fmt("activelearn: dataset file '{}' contains no structures",
                m_dataset_path);
            return;
        }
        if (data[0].AtomCount() != m_anchor.AtomCount() || data[0].elements != m_anchor.Atoms())
            CurcumaLogger::warn_fmt(
                "activelearn: dataset structures ({} atoms) do not match the input molecule "
                "({} atoms); the stored anchor in the model comes from the input",
                data[0].AtomCount(), m_anchor.AtomCount());
        const double seconds = std::chrono::duration<double>(
            std::chrono::high_resolution_clock::now() - t0).count();
        CurcumaLogger::result_fmt(
            "activelearn: {} labelled structures read from {} in {:.2f} s (no new labels)",
            data.size(), m_dataset_path, seconds);
    } else {
        if (!labelStructures(molecules, data))
            return;
    }

    if (m_write_dataset && m_dataset_path.empty()) {
        const std::string ds = outputPath(Basename() + ".deltaml.dataset.json");
        data.save(ds);
        data.writeXYZ(outputPath(Basename() + ".deltaml.dataset.xyz"));
        CurcumaLogger::result_fmt("activelearn: dataset written to {}", ds);
    }

    MLDataset train, test;
    data.split(m_test_fraction, m_seed, train, test);
    CurcumaLogger::result_fmt("activelearn: {} training / {} held-out structures", train.size(),
        test.size());
    if (train.empty()) {
        CurcumaLogger::error("activelearn: training split is empty; lower test_fraction");
        return;
    }

    // Optional sweep over the number of centres. This is the measurement that decides how
    // sparse the model may be: accuracy against evaluation cost, for this system and this
    // pair of methods.
    if (!m_center_sweep.empty()) {
        CurcumaLogger::result("activelearn: centre sweep (accuracy over evaluation cost)");
        for (int m : m_center_sweep) {
            DeltaMLModel probe;
            if (!fitAndReport(train, test, m, probe, false))
                CurcumaLogger::warn_fmt("activelearn: fit with M = {} failed", m);
        }
    }

    // Main (joint energy + gradient) model. With the energy model enabled its held-out block
    // is printed after the energy model's, so it stays the last block in the log; scripts and
    // the CLI tests read the main model's numbers from the last "gradient RMSE" line.
    const std::string main_label = m_energy_model ? "joint model" : "";
    DeltaMLModel model;
    if (!fitAndReport(train, test, trainingOptions(m_n_centers), model, true, main_label,
            !m_energy_model))
        return;

    std::string path = m_model_path;
    if (path.empty())
        path = Basename() + ".deltaml.json";
    path = outputPath(path);
    if (model.save(path))
        CurcumaLogger::success_fmt("activelearn: model written to {}", path);

    // Energy-only second model (Claude Generated, October 2026). Same data, split and
    // descriptor; only the gradient rows are dropped and M / lambda / sigma have their own
    // settings. Motivation (90-atom conformer set, GFN-FF -> GFN2, 5 split seeds): the joint
    // fit trades energy accuracy for forces and ends above the base method in relative
    // energy MAE, the energy-only fit ends below it.
    if (m_energy_model) {
        DeltaMLTrainingOptions eopt = trainingOptions(
            m_energy_model_n_centers > 0 ? m_energy_model_n_centers : m_n_centers);
        eopt.use_forces = false;
        eopt.lambda = m_energy_model_lambda;
        eopt.kernel_sigma = m_energy_model_kernel_sigma;

        DeltaMLModel emodel;
        if (fitAndReport(train, test, eopt, emodel, true, "energy model", false)) {
            const std::string epath = energyModelPath(path);
            if (emodel.save(epath))
                CurcumaLogger::success_fmt("activelearn: energy model written to {}", epath);

            if (!test.empty()) {
                CurcumaLogger::result("activelearn: held-out comparison of the two models");
                logErrors("held-out set (energy model: energies only, for -sp energies)",
                    emodel.evaluate(test));
                logErrors("held-out set (joint model: energies + gradients, for -opt/-md)",
                    model.evaluate(test));
            }
            CurcumaLogger::result_fmt(
                "activelearn: relative energies with  -method deltaml:{}  (forces no better than the base)",
                epath);
        } else {
            CurcumaLogger::warn("activelearn: energy-only fit failed; only the main model was written");
            if (!test.empty())
                logErrors("held-out set (joint model)", model.evaluate(test));
        }
    }

    CurcumaLogger::result_fmt("activelearn: use it with  -method deltaml:{}", path);
}

void ActiveLearn::runTest()
{
    if (m_model_path.empty()) {
        CurcumaLogger::error("activelearn: -al_mode test needs -model <file.json>");
        return;
    }

    DeltaMLModel model;
    if (!model.load(m_model_path))
        return;

    m_base_method = model.baseMethod();
    if (!model.referenceMethod().empty())
        m_ref_method = model.referenceMethod();

    std::vector<Molecule> molecules;
    if (!loadStructures(molecules))
        return;

    MLDataset data;
    if (!labelStructures(molecules, data))
        return;

    CurcumaLogger::result_fmt("activelearn: evaluating model {} ({} centres)", m_model_path,
        model.centerCount());
    logErrors("test set", model.evaluate(data));
}

// =================================================================================
// Active learning
// =================================================================================

std::vector<Molecule> ActiveLearn::sampleStructures(
    const std::string& md_method, const std::vector<Molecule>& seeds, int cycle)
{
    std::vector<Molecule> frames;
    if (seeds.empty())
        return frames;

    json md = m_controller;

    // SimpleMD reads its settings through ConfigManager("simplemd", ...), which resolves the
    // module sub-document first. A bare top-level key is not reliably seen -- 'method' in
    // particular falls back to the module default ("uff") -- and from the command line these
    // flags only work because the registry routes them into controller["simplemd"] for us.
    // So write them where the module actually looks, and keep the flat copies for anything
    // downstream that reads the controller directly.
    json sim = md.value("simplemd", json::object());
    const auto set_both = [&md, &sim](const std::string& key, const json& value) {
        sim[key] = value;
        md[key] = value;
    };
    set_both("method", md_method);
    set_both("temperature", m_al_temperature);
    set_both("initial_temperature", m_al_temperature);
    set_both("max_time", m_al_md_time);
    set_both("dump_frequency", m_al_dump_frequency);
    set_both("print_frequency", 1000000);
    set_both("write_xyz", true);
    set_both("threads", 1);
    set_both("verbosity", 0);
    md["simplemd"] = sim;

    // On-the-fly: end each run where the model first loses its footing, and have the method
    // write those structures out. They are the candidates -- no separate selection needed,
    // because the trajectory itself picked them.
    const bool trigger_mode = m_al_onthefly && md_method.rfind("deltaml", 0) == 0;
    std::string record_path;
    if (trigger_mode) {
        record_path = outputPath(Basename() + ".al" + std::to_string(cycle) + ".triggered.xyz");
        std::error_code ec;
        std::filesystem::remove(record_path, ec);
        json dml = md.value("deltaml", json::object());
        dml["uncertainty_abort"] = m_al_onthefly_threshold;
        dml["uncertainty_record"] = record_path;
        // Record a little before the stop, so the approach to the excursion is captured too.
        dml["uncertainty_warn"] = 0.9 * m_al_onthefly_threshold;
        md["deltaml"] = dml;
        CurcumaLogger::result_fmt(
            "activelearn: on-the-fly sampling on '{}' -- runs end above indicator {:.2f}, "
            "trigger structures go to {}",
            md_method, m_al_onthefly_threshold, record_path);
    }

    // SimpleMD redirects its whole output into "<basename>.snapshots" and then builds file
    // paths as outputDir + "/" + basename. A basename containing a directory separator
    // therefore points into a nested directory that nobody creates, and the trajectory is
    // dropped without an error. So the run tag stays a plain name and the files are
    // collected from the working directory afterwards.
    std::string tag = Basename() + "_al" + std::to_string(cycle);
    std::replace(tag.begin(), tag.end(), '/', '_');

    std::vector<MDThread*> threads;
    CxxThreadPool* pool = new CxxThreadPool;
    for (int r = 0; r < m_al_md_runs; ++r) {
        json cfg = md;
        // Distinct seeds, otherwise every run reproduces the same trajectory.
        cfg["seed"] = static_cast<int>(m_seed) + 1000 * cycle + r;
        cfg["simplemd"]["seed"] = cfg["seed"];
        MDThread* t = new MDThread(cfg);
        t->setAutoDelete(false);
        t->setThreadId(r);
        t->setBasename(tag);
        t->setMolecule(seeds[r % seeds.size()]);
        threads.push_back(t);
        pool->addThread(t);
    }
    pool->setActiveThreadCount(std::max(1, m_threads));

    const int saved_verbosity = CurcumaLogger::get_verbosity();
    CurcumaLogger::set_verbosity(0);
    pool->StartAndWait();

    delete pool;
    for (auto* t : threads)
        delete t;
    // Restore only after the threads are gone: each MDThread owns a SimpleMD, and the
    // CurcumaMethod destructor restores the global logger level it captured at
    // construction -- which is the silenced one. Restoring earlier gets clobbered.
    CurcumaLogger::set_verbosity(saved_verbosity);

    // Collect the trajectories, then remove the scratch the dynamics left behind.
    // Claude Generated (Oct 2026): directory_iterator has no defined order, and the frame order
    // feeds the selection (start point and ties of the farthest-point pick) and the draw of the
    // validation frames. Read the runs and their files in sorted order, so a fixed seed gives a
    // fixed candidate list.
    std::vector<std::filesystem::path> scratch;
    for (const auto& entry : std::filesystem::directory_iterator(std::filesystem::current_path())) {
        const std::string name = entry.path().filename().string();
        if (name.rfind(tag + ".t", 0) != 0)
            continue;
        scratch.push_back(entry.path());
    }
    std::sort(scratch.begin(), scratch.end());
    for (const auto& dir : scratch) {
        if (!std::filesystem::is_directory(dir))
            continue;
        std::vector<std::string> trajectories;
        for (const auto& f : std::filesystem::directory_iterator(dir)) {
            const std::string p = f.path().string();
            if (p.size() < 8 || p.substr(p.size() - 8) != ".trj.xyz")
                continue;
            trajectories.push_back(p);
        }
        std::sort(trajectories.begin(), trajectories.end());
        for (const auto& p : trajectories) {
            FileIterator iter(p, true);
            while (!iter.AtEnd())
                frames.push_back(iter.Next());
        }
    }
    for (const auto& p : scratch) {
        std::error_code ec;
        std::filesystem::remove_all(p, ec);
    }

    m_triggered.clear();
    if (trigger_mode && std::filesystem::exists(record_path)) {
        FileIterator iter(record_path, true);
        while (!iter.AtEnd())
            m_triggered.push_back(iter.Next());
        CurcumaLogger::result_fmt(
            "activelearn: {} structures tripped the {:.2f} threshold and ended their run",
            m_triggered.size(), m_al_onthefly_threshold);
    }

    return frames;
}

std::vector<Molecule> ActiveLearn::selectCandidates(const std::vector<Molecule>& candidates,
    const DeltaMLModel& model, bool have_model, int wanted, double& max_uncertainty,
    std::vector<std::size_t>* picked_indices, bool ignore_threshold)
{
    max_uncertainty = 0.0;
    std::vector<Molecule> picked;
    if (picked_indices)
        picked_indices->clear();
    if (candidates.empty() || wanted <= 0)
        return picked;

    // Filter by unfamiliarity first.
    std::vector<std::size_t> pool;
    for (std::size_t i = 0; i < candidates.size(); ++i) {
        double u = 1.0;
        if (have_model)
            u = model.uncertainty(candidates[i].Atoms(), candidates[i].getGeometry());
        max_uncertainty = std::max(max_uncertainty, u);
        if (!have_model || ignore_threshold || u > m_al_uncertainty)
            pool.push_back(i);
    }
    if (pool.empty())
        return picked;

    // Structure-level fingerprint: the mean atomic descriptor. Cheap, and it is the same
    // space the model measures similarity in.
    const std::vector<int> species = m_anchor.Atoms();
    curcuma::ml::ACSFConfig cfg = have_model
        ? model.descriptorConfig()
        : curcuma::ml::ACSFConfig::defaultFor(species, m_acsf);
    const curcuma::ml::ACSFDescriptor desc(cfg);

    Matrix fingerprints(static_cast<int>(pool.size()), cfg.featureCount());
    for (std::size_t k = 0; k < pool.size(); ++k) {
        Matrix X;
        desc.compute(candidates[pool[k]].Atoms(), candidates[pool[k]].getGeometry(), X, nullptr);
        fingerprints.row(static_cast<int>(k)) = X.colwise().mean();
    }

    // Greedy farthest point, started from the most unfamiliar structure.
    const int n = static_cast<int>(pool.size());
    std::vector<char> taken(n, 0);
    int current = 0;
    if (have_model) {
        double worst = -1.0;
        for (int k = 0; k < n; ++k) {
            const double u
                = model.uncertainty(candidates[pool[k]].Atoms(), candidates[pool[k]].getGeometry());
            if (u > worst) {
                worst = u;
                current = k;
            }
        }
    }

    Vector mind = Vector::Constant(n, std::numeric_limits<double>::max());
    const int take = std::min(wanted, n);
    for (int p = 0; p < take; ++p) {
        taken[current] = 1;
        picked.push_back(candidates[pool[current]]);
        if (picked_indices)
            picked_indices->push_back(pool[current]);
        for (int k = 0; k < n; ++k) {
            if (taken[k])
                continue;
            const double d = (fingerprints.row(k) - fingerprints.row(current)).squaredNorm();
            if (d < mind(k))
                mind(k) = d;
        }
        double far = -1.0;
        int next = -1;
        for (int k = 0; k < n; ++k)
            if (!taken[k] && mind(k) > far) {
                far = mind(k);
                next = k;
            }
        if (next < 0)
            break;
        current = next;
    }

    return picked;
}

std::vector<Molecule> ActiveLearn::drawValidationFrames(const std::vector<Molecule>& frames,
    const std::vector<std::size_t>& excluded, int cycle, std::vector<std::size_t>* drawn) const
{
    if (drawn)
        drawn->clear();
    std::vector<Molecule> out;
    if (m_al_validation_frames <= 0 || frames.empty())
        return out;
    std::vector<char> skip(frames.size(), 0);
    for (std::size_t i : excluded)
        if (i < skip.size())
            skip[i] = 1;
    std::vector<std::size_t> free;
    for (std::size_t i = 0; i < frames.size(); ++i)
        if (!skip[i])
            free.push_back(i);

    // Partial Fisher-Yates shuffle with an engine of its own, seeded from (seed, cycle): the
    // draw neither consumes nor perturbs any other random stream of the run.
    std::seed_seq seq{ static_cast<unsigned>(m_seed), static_cast<unsigned>(cycle), 0x76616cu };
    std::mt19937 rng(seq);
    const std::size_t take = std::min(free.size(), static_cast<std::size_t>(m_al_validation_frames));
    for (std::size_t k = 0; k < take; ++k) {
        std::uniform_int_distribution<std::size_t> pick(k, free.size() - 1);
        std::swap(free[k], free[pick(rng)]);
    }
    free.resize(take);
    // Label in trajectory order, so the dataset order does not depend on the shuffle.
    std::sort(free.begin(), free.end());
    for (std::size_t i : free)
        out.push_back(frames[i]);
    if (drawn)
        *drawn = free;
    return out;
}

ActiveLearn::PooledValidation ActiveLearn::evaluatePool(
    const DeltaMLModel& model, const std::vector<const MLDataset*>& pool) const
{
    const double kcal = 627.5094740631;
    std::vector<double> err; // model minus reference [kcal/mol]
    for (const MLDataset* set : pool)
        for (std::size_t s = 0; s < set->size(); ++s) {
            const MLStructure& st = (*set)[s];
            const double corr = model.predictEnergy(st.elements, st.geometry);
            err.push_back((st.energy_base + corr - st.energy_ref) * kcal);
        }
    PooledValidation p;
    p.n = static_cast<int>(err.size());
    if (err.empty())
        return p;
    double mean = 0.0;
    for (double e : err)
        mean += e;
    mean /= static_cast<double>(err.size());
    p.offset = mean;
    std::vector<double> c(err.size());
    double mae = 0.0;
    for (std::size_t i = 0; i < err.size(); ++i) {
        c[i] = err[i] - mean;
        mae += std::abs(c[i]);
    }
    p.mae_centred = mae / static_cast<double>(c.size());
    // Lower (1 - q) quantile, type 7: position (1 - q)(n - 1) in the sorted list.
    std::sort(c.begin(), c.end());
    const double pos = (1.0 - m_al_hole_quantile) * static_cast<double>(c.size() - 1);
    const std::size_t lo = static_cast<std::size_t>(std::floor(pos));
    const std::size_t hi = std::min(lo + 1, c.size() - 1);
    p.hole = c[lo] + (pos - static_cast<double>(lo)) * (c[hi] - c[lo]);
    return p;
}

void ActiveLearn::runLearn()
{
    std::vector<Molecule> seeds;
    if (!loadStructures(seeds))
        return;

    CurcumaLogger::result_fmt(
        "activelearn: {} cycles, base = {}, reference = {}, sampling at {:.0f} K",
        m_al_cycles, m_base_method, m_ref_method, m_al_temperature);
    CurcumaLogger::warn("activelearn: sampling runs on the surface being learned, so the "
                        "first cycle only sees what the base method considers relevant. "
                        "Elevated temperature mitigates this bias but does not remove it.");

    MLDataset data;
    DeltaMLModel model;
    bool have_model = false;

    const std::string model_path = outputPath(
        m_model_path.empty() ? (Basename() + ".deltaml.json") : m_model_path);
    const std::string curve_path = outputPath(Basename() + ".deltaml.learning_curve.dat");
    std::ofstream curve(curve_path);
    // Columns 1-6 as before October 2026. The val_* columns (Claude Generated, October 2026) are
    // the prospective validation of the model ENTERING the cycle, on frames of this cycle's
    // sampling dynamics, measured before the refit; nan where there is none (cycle 0, or
    // al_validation_frames 0). val_hole_kcal is the most negative model-minus-reference energy.
    // The pool_* columns (Claude Generated, October 2026) are what the validation criterion
    // tests: the validation frames of the last al_window cycles re-predicted with the model
    // entering the cycle; centred MAE, offset (mean signed error) and the hole quantile of the
    // centred signed errors. A line whose training columns are nan is the cycle that converged
    // (no refit).
    curve << "# cycle  n_train  E_MAE_kcal  F_RMSE_Eh_A  spearman  max_uncertainty"
             "  val_n  val_E_MAE_kcal  val_mean_signed_kcal  val_hole_kcal  val_F_RMSE_Eh_A"
             "  pool_n  pool_E_MAE_centred_kcal  pool_offset_kcal  pool_hole_q_kcal\n";

    const double kcal = 627.5094740631;
    const bool criterion_known = m_al_criterion == "uncertainty" || m_al_criterion == "validation"
        || m_al_criterion == "both";
    if (!criterion_known) {
        CurcumaLogger::error_fmt(
            "activelearn: unknown al_criterion '{}' (uncertainty | validation | both)", m_al_criterion);
        return;
    }
    if (m_al_criterion != "uncertainty" && m_al_validation_frames <= 0)
        CurcumaLogger::warn_fmt("activelearn: al_criterion '{}' needs al_validation_frames > 0; "
                                "with 0 the run cannot converge and runs all {} cycles",
            m_al_criterion, m_al_cycles);
    // Validation pool (Claude Generated, October 2026). The validation frames of the last
    // `window` cycles, kept out of the training data while they are in the pool so that the
    // pooled error stays out of sample. With 'uncertainty' the frames join training at once, as
    // before October 2026 (the learning stays bit-identical), so the pool is one cycle there.
    // Cycle 0 has no model: its frames are no validation and go straight into training.
    const bool hold_out = m_al_criterion != "uncertainty";
    const std::size_t window = hold_out ? static_cast<std::size_t>(m_al_window) : 1;
    struct PoolSet {
        int cycle;
        MLDataset set;
    };
    std::vector<PoolSet> pool_sets; // oldest first
    // Move all but the newest `keep` pool sets into the training data, oldest first.
    const auto flush_pool = [&](std::size_t keep) {
        while (pool_sets.size() > keep) {
            for (std::size_t i = 0; i < pool_sets.front().set.size(); ++i) {
                MLStructure s = pool_sets.front().set[i];
                s.cycle = pool_sets.front().cycle;
                data.add(s);
            }
            pool_sets.erase(pool_sets.begin());
        }
    };
    // Columns 7-15 of the curve: single-cycle validation (val_*) and pool (pool_*).
    const auto write_validation_columns
        = [&](bool validated, const DeltaMLErrors& ve, const PooledValidation& pv) {
              if (validated)
                  curve << "  " << ve.n_structures << "  " << ve.energy_mae_centered * kcal << "  "
                        << ve.energy_mean_signed * kcal << "  " << ve.energy_min_signed * kcal
                        << "  " << ve.force_rmse << "  " << pv.n << "  " << pv.mae_centred << "  "
                        << pv.offset << "  " << pv.hole;
              else
                  curve << "  0  nan  nan  nan  nan  0  nan  nan  nan";
          };

    CurcumaLogger::result_fmt(
        "activelearn: convergence criterion '{}' in {} consecutive cycle(s); indicator <= {:.4f}, "
        "validation pool of {} cycle(s) x {} frames: centred E_MAE <= {:.3f} kcal/mol and hole "
        "quantile (lowest {:.0f} %) >= -{:.3f} kcal/mol",
        m_al_criterion, m_al_patience, m_al_uncertainty, window, m_al_validation_frames,
        m_al_target_energy, 100.0 * (1.0 - m_al_hole_quantile), m_al_target_hole);
    int streak = 0;
    bool converged = false;
    bool stopped_early = false; ///< loop left before al_cycles for a reason other than convergence

    for (int cycle = 0; cycle < m_al_cycles; ++cycle) {
        CurcumaLogger::result_fmt("activelearn: --- cycle {} of {} ---", cycle + 1, m_al_cycles);

        // Sampling surface: explicit override, else the base method until a model exists
        // and the corrected surface from then on.
        std::string md_method = m_al_md_method;
        if (md_method.empty())
            md_method = have_model ? ("deltaml:" + model_path) : m_base_method;
        std::vector<Molecule> pool_seeds = seeds;
        if (!data.empty()) {
            // Restart the dynamics from structures already in the set, so later cycles
            // explore outward from what has been visited rather than from the input again.
            pool_seeds.clear();
            for (int r = 0; r < m_al_md_runs; ++r) {
                const std::size_t idx = (static_cast<std::size_t>(r) * 7 + cycle) % data.size();
                // Reuse the input molecule for its element list and metadata, swapping in
                // the stored geometry.
                pool_seeds.push_back(seeds[0]);
                pool_seeds.back().setGeometry(data[idx].geometry);
            }
        }

        const std::vector<Molecule> frames = sampleStructures(md_method, pool_seeds, cycle);
        if (frames.empty()) {
            CurcumaLogger::warn("activelearn: sampling produced no frames; stopping");
            stopped_early = true;
            break;
        }

        double max_unc = 0.0;
        const int wanted = (cycle == 0) ? m_al_initial : m_al_candidates;
        // In on-the-fly mode the trajectories already nominated their candidates by stopping
        // on them; selection then only has to thin them down for diversity.
        const std::vector<Molecule>& pool = m_triggered.empty() ? frames : m_triggered;
        std::vector<std::size_t> picked_idx;
        std::vector<Molecule> picked
            = selectCandidates(pool, model, have_model, wanted, max_unc, &picked_idx);

        // Report the indicator together with the descriptor-space distance it corresponds
        // to, in units of the kernel width. The bounded indicator alone saturates near one
        // and then says nothing about how much further out a structure actually is.
        const double d_over_sigma
            = (max_unc >= 1.0) ? std::numeric_limits<double>::infinity()
                               : std::sqrt(std::max(0.0, -2.0 * std::log(1.0 - max_unc)));
        CurcumaLogger::result_fmt(
            "activelearn: {} frames sampled, largest extrapolation indicator {:.4f} "
            "(nearest centre at {:.2f} sigma), {} selected for reference calculation",
            frames.size(), max_unc, d_over_sigma, picked.size());

        // Prospective validation (Claude Generated, October 2026): frames drawn at random from
        // this cycle's dynamics, labelled and scored with the model that drove the dynamics,
        // before any refit. They come from the frame list, never from the candidates picked
        // above, so no structure is labelled twice. Cycle 0 has no model: its frames are only
        // added to the training data.
        MLDataset validation;
        std::vector<std::size_t> validation_idx; // indices in `frames`
        bool validated = false;
        DeltaMLErrors ve;
        if (m_al_validation_frames > 0) {
            const std::vector<std::size_t> excluded
                = (&pool == &frames) ? picked_idx : std::vector<std::size_t>{};
            const std::vector<Molecule> vframes
                = drawValidationFrames(frames, excluded, cycle, &validation_idx);
            if (!vframes.empty()) {
                if (!labelStructures(vframes, validation)) {
                    stopped_early = true;
                    break;
                }
                if (have_model) {
                    ve = model.evaluate(validation);
                    validated = true;
                    CurcumaLogger::result_fmt(
                        "activelearn: validation on {} new frames (model of the previous cycle, before "
                        "refit): E_MAE {:.4f} kcal/mol, mean signed {:+.4f} kcal/mol, deepest hole "
                        "{:+.4f} kcal/mol, gradient RMSE {:.6f} Eh/A",
                        ve.n_structures, ve.energy_mae_centered * kcal, ve.energy_mean_signed * kcal,
                        ve.energy_min_signed * kcal, ve.force_rmse);
                } else {
                    CurcumaLogger::result_fmt("activelearn: {} validation frames labelled; no model "
                                              "yet, they only enter the training data",
                        validation.size());
                }
            }
        }
        // Old single-cycle path, kept for 'uncertainty' (window 1) and for cycle 0.
        const auto add_validation = [&]() {
            for (std::size_t i = 0; i < validation.size(); ++i) {
                MLStructure s = validation[i];
                s.cycle = cycle;
                data.add(s);
            }
            validation = MLDataset{};
        };
        // Pooled validation (Claude Generated, October 2026): the current model on the frames
        // of the last `window` cycles, centred errors, hole quantile.
        PooledValidation pv;
        bool pool_full = false;
        if (validated && hold_out) {
            pool_sets.push_back(PoolSet{ cycle, validation });
            validation = MLDataset{}; // now owned by the pool, joins training via flush_pool
            std::vector<const MLDataset*> pool;
            for (const PoolSet& p : pool_sets)
                pool.push_back(&p.set);
            pv = evaluatePool(model, pool);
            pool_full = pool_sets.size() >= window;
        } else if (validated) {
            pv = evaluatePool(model, { &validation });
            pool_full = true;
        }
        if (validated)
            CurcumaLogger::result_fmt(
                "activelearn: validation pool {} frames from {} cycle(s){} (model entering this "
                "cycle): centred E_MAE {:.4f} kcal/mol, offset {:+.4f} kcal/mol, hole (lowest "
                "{:.0f} %) {:+.4f} kcal/mol",
                pv.n, hold_out ? pool_sets.size() : std::size_t{ 1 },
                pool_full ? "" : fmt::format(", window of {} not yet full", window),
                pv.mae_centred, pv.offset, 100.0 * (1.0 - m_al_hole_quantile), pv.hole);

        // Convergence: the chosen criterion in m_al_patience consecutive cycles.
        if (have_model) {
            const bool unc_ok = max_unc <= m_al_uncertainty;
            const bool val_ok = validated && pool_full && pv.mae_centred <= m_al_target_energy
                && pv.hole >= -m_al_target_hole;
            // The plain uncertainty criterion keeps its old meaning: stop at the first cycle in which
            // it is met. Patience applies to the validation part only (Claude Generated, Oct 2026):
            // 'both' stops when the validation held in al_patience consecutive cycles AND the
            // indicator holds in the current one. Requiring the indicator twice in a row as well
            // made 'both' unreachable: on peo8 (tuning seeds 4-6, 20 cycles) the indicator was
            // <= 0.3 in single cycles only, never in two consecutive ones.
            const int patience_needed = (m_al_criterion == "uncertainty") ? 1 : m_al_patience;
            const bool counted = (m_al_criterion == "uncertainty") ? unc_ok : val_ok;
            streak = counted ? streak + 1 : 0;
            const bool stop_now = streak >= patience_needed
                && (m_al_criterion != "both" || unc_ok);
            CurcumaLogger::result_fmt(
                "activelearn: criterion '{}': uncertainty {} (this cycle), validation {}; {} met in "
                "{} of {} consecutive cycle(s)",
                m_al_criterion, unc_ok ? "met" : "not met",
                validated ? (val_ok ? "met" : "not met") : "n/a",
                (m_al_criterion == "uncertainty") ? "uncertainty" : "validation", streak,
                patience_needed);
            if (stop_now) {
                // The validation labels are paid for: they stay in the dataset file, but the
                // model is not refitted, so the written model is the one that was validated.
                add_validation();
                flush_pool(0);
                const std::string what = (m_al_criterion == "uncertainty")
                    ? fmt::format("uncertainty (no sampled structure above the indicator threshold "
                                  "{:.4f}; certifies coverage of the sampled region, not accuracy)",
                          m_al_uncertainty)
                    : (m_al_criterion == "validation")
                    ? fmt::format("validation (pool of {} frames: centred E_MAE {:.4f} <= {:.4f} "
                                  "kcal/mol, hole {:+.4f} >= -{:.4f} kcal/mol)",
                          pv.n, pv.mae_centred, m_al_target_energy, pv.hole, m_al_target_hole)
                    : fmt::format("both (indicator {:.4f} <= {:.4f}; pool of {} frames: centred "
                                  "E_MAE {:.4f} <= {:.4f} kcal/mol, hole {:+.4f} >= -{:.4f} kcal/mol)",
                          max_unc, m_al_uncertainty, pv.n, pv.mae_centred, m_al_target_energy,
                          pv.hole, m_al_target_hole);
                CurcumaLogger::success_fmt("activelearn: converged in cycle {} by criterion {}; "
                                           "met in {} consecutive cycle(s)",
                    cycle + 1, what, streak);
                // Curve line of the converged cycle: no refit, so the training columns are nan;
                // n_train is the final number of labelled structures.
                curve << cycle << "  " << data.size() << "  nan  nan  nan  " << max_unc;
                write_validation_columns(validated, ve, pv);
                curve << "\n";
                curve.flush();
                converged = true;
                break;
            }
        }
        // Claude Generated (Oct 2026): no frame above the indicator threshold, yet the run has
        // not converged (a validation criterion still fails). Stopping here would end the run
        // because the indicator is satisfied although the error is not (peo8: indicator 0.26,
        // validation hole -7.9 kcal/mol). Pick among all frames instead, farthest point first.
        if (picked.empty() && have_model && m_al_criterion != "uncertainty") {
            // The validation frames are labelled already; leave them out so none is labelled twice.
            std::vector<Molecule> rest;
            if (&pool == &frames) {
                std::vector<char> used(frames.size(), 0);
                for (std::size_t i : validation_idx)
                    used[i] = 1;
                for (std::size_t i = 0; i < frames.size(); ++i)
                    if (!used[i])
                        rest.push_back(frames[i]);
            } else {
                rest = pool;
            }
            double unused = 0.0;
            picked = selectCandidates(rest, model, have_model, wanted, unused, nullptr, true);
            CurcumaLogger::warn_fmt(
                "activelearn: no frame above the indicator threshold {:.4f}, but criterion '{}' "
                "is not met; {} candidates selected among all frames instead",
                m_al_uncertainty, m_al_criterion, picked.size());
        }
        if (picked.empty()) {
            add_validation();
            flush_pool(0);
            CurcumaLogger::warn("activelearn: nothing selected; stopping");
            stopped_early = true;
            break;
        }

        MLDataset fresh;
        if (!labelStructures(picked, fresh)) {
            stopped_early = true;
            break;
        }
        for (std::size_t i = 0; i < fresh.size(); ++i) {
            MLStructure s = fresh[i];
            s.cycle = cycle;
            data.add(s);
        }
        add_validation();
        // Pool sets that will not be in the next cycle's pool join the training data now; in the
        // last cycle all of them, so the final model is fitted on every label.
        flush_pool(cycle + 1 >= m_al_cycles ? 0 : window - 1);

        MLDataset train, test;
        data.split(m_test_fraction, m_seed, train, test);
        if (train.empty()) {
            CurcumaLogger::warn("activelearn: training split empty, collecting more data");
            continue;
        }

        DeltaMLModel candidate_model;
        if (!fitAndReport(train, test, m_n_centers, candidate_model, true)) {
            stopped_early = true;
            break;
        }
        model = candidate_model;
        have_model = true;
        model.save(model_path);

        const DeltaMLErrors e = test.empty() ? DeltaMLErrors{} : model.evaluate(test);
        curve << cycle << "  " << data.size() << "  " << e.energy_mae_centered * 627.5094740631
              << "  " << e.force_rmse << "  " << e.spearman << "  " << max_unc;
        write_validation_columns(validated, ve, pv);
        curve << "\n";
        curve.flush();
    }
    flush_pool(0); // every paid label ends up in the dataset file, also after an early stop
    if (!converged && !stopped_early)
        CurcumaLogger::warn_fmt("activelearn: not converged -- criterion '{}' was not met in {} "
                                "consecutive cycle(s) within al_cycles = {}",
            m_al_criterion, m_al_patience, m_al_cycles);
    else if (!converged)
        CurcumaLogger::warn_fmt("activelearn: not converged -- the run stopped early (see above) "
                                "before criterion '{}' was met", m_al_criterion);

    if (!have_model) {
        CurcumaLogger::error("activelearn: no model was produced");
        return;
    }

    if (m_write_dataset) {
        data.save(outputPath(Basename() + ".deltaml.dataset.json"));
        data.writeXYZ(outputPath(Basename() + ".deltaml.dataset.xyz"));
    }
    CurcumaLogger::success_fmt("activelearn: model written to {} ({} labelled structures)",
        model_path, data.size());
    CurcumaLogger::result_fmt("activelearn: learning curve in {}", curve_path);
    CurcumaLogger::result_fmt("activelearn: use it with  -method deltaml:{}", model_path);
    CurcumaLogger::warn("activelearn: before trusting the model in dynamics, check energy "
                        "conservation in an NVE run -- a slightly rippled learned surface is "
                        "invisible in the energy error and still ruins a trajectory.");
}

void ActiveLearn::start()
{
    if (Filename().empty()) {
        CurcumaLogger::error("activelearn: no input structure file given");
        return;
    }
    if (!m_acsf_valid) {
        CurcumaLogger::error("activelearn: invalid descriptor settings (acsf_*), see above");
        return;
    }
    if (!m_dataset_path.empty() && m_mode != "train")
        CurcumaLogger::warn_fmt("activelearn: -dataset is only used in 'train' mode; ignored in '{}'",
            m_mode);

    if (m_mode == "train") {
        runTrain();
    } else if (m_mode == "test") {
        runTest();
    } else if (m_mode == "learn") {
        runLearn();
    } else {
        CurcumaLogger::error_fmt("activelearn: unknown mode '{}' (train | test | learn)", m_mode);
    }
}

void ActiveLearn::printHelp() const
{
    std::cout << "\n";
    std::cout << "curcuma -activelearn <structures.xyz> [options]\n\n";
    std::cout << "Fits a delta-ML model that lifts a cheap base method to a more expensive\n";
    std::cout << "reference, so conformer search and pre-optimisation run at base-method cost.\n\n";
    std::cout << "Modes:\n";
    std::cout << "  -al_mode train   fit a model to the given structures (default)\n";
    std::cout << "  -al_mode learn   generate structures by dynamics, label and refit in cycles\n";
    std::cout << "  -al_mode test    evaluate an existing -model against fresh references\n\n";
    std::cout << "Examples:\n";
    std::cout << "  curcuma -activelearn conf.xyz -method gfnff -ref_method gfn2 -threads 8\n";
    std::cout << "  curcuma -activelearn conf.xyz -center_sweep 200,500,1000,2000\n";
    std::cout << "  curcuma -activelearn conf.xyz -activelearn.dataset conf.deltaml.dataset.json \\\n";
    std::cout << "          -activelearn.acsf_n_shell 12 -activelearn.acsf_zeta 1,2,4,8,16\n";
    std::cout << "  curcuma -activelearn mol.xyz -al_mode learn -al_cycles 6 -threads 8\n";
    std::cout << "  curcuma -opt mol.xyz -method deltaml:conf.deltaml.json\n";
    std::cout << "  curcuma -sp  mol.xyz -method deltaml:conf.deltaml.energy.json\n\n";
    std::cout << "'train' writes two models: <basename>.deltaml.json (energies + gradients, for\n";
    std::cout << "-opt/-md) and <basename>.deltaml.energy.json (energies only, for -sp; its\n";
    std::cout << "forces are no better than the base method). -activelearn.energy_model false\n";
    std::cout << "skips the second one.\n\n";
    std::cout << "The reference must carry the same dispersion treatment as the base method,\n";
    std::cout << "otherwise the long-range dispersion term ends up in the learned difference\n";
    std::cout << "where a local descriptor cannot represent it.\n\n";
}
