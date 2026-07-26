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

#include "src/core/curcuma_logger.h"
#include "src/core/energycalculator.h"
#include "src/core/fileiterator.h"
#include "src/core/parameter_registry.h"
#include "src/core/units.h"

#include "external/CxxThreadPool/include/CxxThreadPool.hpp"

#include <chrono>
#include <fstream>
#include <sstream>

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
 * Each thread owns its own pair of calculators. Both are initialised from the anchor
 * geometry and then moved to the actual structure, so a method that freezes its topology
 * at setup produces the same bond list in every thread. Without this the learned
 * difference absorbs setup noise instead of physics.
 */
class LabelThread : public CxxThread {
public:
    LabelThread(const std::string& base_method, const std::string& ref_method, const json& controller,
        const Molecule& anchor)
        : m_base_method(base_method)
        , m_ref_method(ref_method)
        , m_controller(controller)
        , m_anchor(anchor)
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

        json cfg = m_controller;
        cfg["threads"] = 1;
        // Thousands of per-structure SCF banners would bury the actual result.
        cfg["verbosity"] = 0;

        try {
            EnergyCalculator base(m_base_method, cfg);
            EnergyCalculator ref(m_ref_method, cfg);
            base.setMolecule(m_anchor.getMolInfo());
            ref.setMolecule(m_anchor.getMolInfo());

            for (std::size_t k = 0; k < m_molecules.size(); ++k) {
                const Molecule& mol = m_molecules[k];

                MLStructure s;
                s.elements = mol.Atoms();
                s.geometry = mol.getGeometry();
                s.charge = mol.Charge();
                s.spin = mol.Spin();

                base.updateGeometry(s.geometry);
                s.energy_base = base.CalculateEnergy(true);
                s.gradient_base = base.Gradient();

                ref.updateGeometry(s.geometry);
                s.energy_ref = ref.CalculateEnergy(true);
                s.gradient_ref = ref.Gradient();
                s.has_gradients = true;

                m_results.push_back(s);
            }
        } catch (const std::exception& e) {
            m_error = e.what();
            return 1;
        }
        return 0;
    }

    const std::vector<MLStructure>& results() const { return m_results; }
    const std::string& error() const { return m_error; }

private:
    std::string m_base_method;
    std::string m_ref_method;
    json m_controller;
    Molecule m_anchor;
    std::vector<Molecule> m_molecules;
    std::vector<std::size_t> m_indices;
    std::vector<MLStructure> m_results;
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
    m_center_sweep = parseIntList(Json2KeyWord<std::string>(m_defaults, "center_sweep"));

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
    return true;
}

bool ActiveLearn::labelStructures(const std::vector<Molecule>& molecules, MLDataset& data)
{
    const auto t0 = std::chrono::high_resolution_clock::now();
    CurcumaLogger::result_fmt("activelearn: labelling {} structures  base = {}  reference = {}",
        molecules.size(), m_base_method, m_ref_method);

    // One thread per chunk rather than per structure: constructing the calculators
    // (topology setup, parameter generation) costs far more than a single point.
    const int n_threads = std::max(1, m_threads);
    std::vector<LabelThread*> threads;
    threads.reserve(n_threads);
    for (int t = 0; t < n_threads; ++t)
        threads.push_back(new LabelThread(m_base_method, m_ref_method, m_controller, m_anchor));
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

    bool ok = true;
    for (auto* t : threads) {
        if (!t->error().empty()) {
            CurcumaLogger::error_fmt("activelearn: labelling failed: {}", t->error());
            ok = false;
        }
        for (const auto& s : t->results())
            data.add(s);
    }

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

bool ActiveLearn::fitAndReport(const MLDataset& train, const MLDataset& test, int n_centers,
    DeltaMLModel& model, bool log_details)
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

    DeltaMLTrainingReport report;
    if (!model.train(train, opt, report))
        return false;

    model.setBaseMethod(m_base_method);
    model.setReferenceMethod(m_ref_method);
    model.setAnchor(m_anchor.Atoms(), m_anchor.getGeometry());

    if (log_details) {
        CurcumaLogger::result_fmt(
            "activelearn: fitted {} centres from {} environments, {} features, sigma = {:.4f}",
            report.n_centers, report.n_environments, report.n_features, report.kernel_sigma);
        CurcumaLogger::result_fmt("activelearn: {} energy rows, {} gradient rows, solved in {:.1f} s",
            report.n_energy_rows, report.n_force_rows, report.train_seconds);
        if (report.jitter_used > 0.0)
            CurcumaLogger::warn_fmt(
                "activelearn: the normal equations needed a jitter of {:.1e} to factorise; "
                "consider a larger lambda or fewer centres",
                report.jitter_used);
        logErrors("training set", report.train_errors);
    }

    if (!test.empty())
        logErrors(log_details ? "held-out set" : fmt::format("M = {}", report.n_centers),
            model.evaluate(test));

    return true;
}

void ActiveLearn::runTrain()
{
    std::vector<Molecule> molecules;
    if (!loadStructures(molecules))
        return;

    MLDataset data;
    if (!labelStructures(molecules, data))
        return;

    if (m_write_dataset) {
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

    DeltaMLModel model;
    if (!fitAndReport(train, test, m_n_centers, model, true))
        return;

    std::string path = m_model_path;
    if (path.empty())
        path = Basename() + ".deltaml.json";
    path = outputPath(path);
    if (model.save(path))
        CurcumaLogger::success_fmt("activelearn: model written to {}", path);

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

void ActiveLearn::start()
{
    if (Filename().empty()) {
        CurcumaLogger::error("activelearn: no input structure file given");
        return;
    }

    if (m_mode == "train") {
        runTrain();
    } else if (m_mode == "test") {
        runTest();
    } else if (m_mode == "learn") {
        CurcumaLogger::error("activelearn: -al_mode learn is not implemented yet; "
                             "use 'train' with an existing structure set");
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
    std::cout << "  -al_mode test    evaluate an existing -model against fresh references\n\n";
    std::cout << "Examples:\n";
    std::cout << "  curcuma -activelearn conf.xyz -method gfnff -ref_method gfn2 -threads 8\n";
    std::cout << "  curcuma -activelearn conf.xyz -center_sweep 200,500,1000,2000\n";
    std::cout << "  curcuma -opt mol.xyz -method deltaml:conf.deltaml.json\n\n";
    std::cout << "The reference must carry the same dispersion treatment as the base method,\n";
    std::cout << "otherwise the long-range dispersion term ends up in the learned difference\n";
    std::cout << "where a local descriptor cannot represent it.\n\n";
}
