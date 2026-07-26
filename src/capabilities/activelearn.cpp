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

#include "simplemd.h"

#include "external/CxxThreadPool/include/CxxThreadPool.hpp"

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <limits>
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
    m_al_cycles = Json2KeyWord<int>(m_defaults, "al_cycles");
    m_al_temperature = Json2KeyWord<double>(m_defaults, "al_temperature");
    m_al_md_time = Json2KeyWord<double>(m_defaults, "al_md_time");
    m_al_md_runs = Json2KeyWord<int>(m_defaults, "al_md_runs");
    m_al_dump_frequency = Json2KeyWord<int>(m_defaults, "al_dump_frequency");
    m_al_candidates = Json2KeyWord<int>(m_defaults, "al_candidates");
    m_al_uncertainty = Json2KeyWord<double>(m_defaults, "al_uncertainty");
    m_al_initial = Json2KeyWord<int>(m_defaults, "al_initial");
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
    md["method"] = md_method;
    md["temperature"] = m_al_temperature;
    md["initial_temperature"] = m_al_temperature;
    md["max_time"] = m_al_md_time;
    md["dump_frequency"] = m_al_dump_frequency;
    md["print_frequency"] = 1000000;
    md["write_xyz"] = true;
    md["threads"] = 1;
    md["verbosity"] = 0;

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
    std::vector<std::filesystem::path> scratch;
    for (const auto& entry : std::filesystem::directory_iterator(std::filesystem::current_path())) {
        const std::string name = entry.path().filename().string();
        if (name.rfind(tag + ".t", 0) != 0)
            continue;
        scratch.push_back(entry.path());
        if (!entry.is_directory())
            continue;
        for (const auto& f : std::filesystem::directory_iterator(entry.path())) {
            const std::string p = f.path().string();
            if (p.size() < 8 || p.substr(p.size() - 8) != ".trj.xyz")
                continue;
            FileIterator iter(p, true);
            while (!iter.AtEnd())
                frames.push_back(iter.Next());
        }
    }
    for (const auto& p : scratch) {
        std::error_code ec;
        std::filesystem::remove_all(p, ec);
    }

    return frames;
}

std::vector<Molecule> ActiveLearn::selectCandidates(const std::vector<Molecule>& candidates,
    const DeltaMLModel& model, bool have_model, int wanted, double& max_uncertainty)
{
    max_uncertainty = 0.0;
    std::vector<Molecule> picked;
    if (candidates.empty() || wanted <= 0)
        return picked;

    // Filter by unfamiliarity first.
    std::vector<std::size_t> pool;
    for (std::size_t i = 0; i < candidates.size(); ++i) {
        double u = 1.0;
        if (have_model)
            u = model.uncertainty(candidates[i].Atoms(), candidates[i].getGeometry());
        max_uncertainty = std::max(max_uncertainty, u);
        if (!have_model || u > m_al_uncertainty)
            pool.push_back(i);
    }
    if (pool.empty())
        return picked;

    // Structure-level fingerprint: the mean atomic descriptor. Cheap, and it is the same
    // space the model measures similarity in.
    const std::vector<int> species = m_anchor.Atoms();
    curcuma::ml::ACSFConfig cfg = have_model
        ? model.descriptorConfig()
        : curcuma::ml::ACSFConfig::defaultFor(species);
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
    curve << "# cycle  n_train  E_MAE_kcal  F_RMSE_Eh_A  spearman  max_uncertainty\n";

    for (int cycle = 0; cycle < m_al_cycles; ++cycle) {
        CurcumaLogger::result_fmt("activelearn: --- cycle {} of {} ---", cycle + 1, m_al_cycles);

        const std::string md_method = have_model ? ("deltaml:" + model_path) : m_base_method;
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
            break;
        }

        double max_unc = 0.0;
        const int wanted = (cycle == 0) ? m_al_initial : m_al_candidates;
        const std::vector<Molecule> picked
            = selectCandidates(frames, model, have_model, wanted, max_unc);

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

        if (have_model && max_unc <= m_al_uncertainty) {
            CurcumaLogger::success_fmt(
                "activelearn: converged -- no sampled structure exceeds the indicator "
                "threshold {:.4f}. Note this certifies coverage of the sampled region, "
                "not accuracy.",
                m_al_uncertainty);
            break;
        }
        if (picked.empty()) {
            CurcumaLogger::warn("activelearn: nothing selected; stopping");
            break;
        }

        MLDataset fresh;
        if (!labelStructures(picked, fresh))
            break;
        for (std::size_t i = 0; i < fresh.size(); ++i) {
            MLStructure s = fresh[i];
            s.cycle = cycle;
            data.add(s);
        }

        MLDataset train, test;
        data.split(m_test_fraction, m_seed, train, test);
        if (train.empty()) {
            CurcumaLogger::warn("activelearn: training split empty, collecting more data");
            continue;
        }

        DeltaMLModel candidate_model;
        if (!fitAndReport(train, test, m_n_centers, candidate_model, true))
            break;
        model = candidate_model;
        have_model = true;
        model.save(model_path);

        const DeltaMLErrors e = test.empty() ? DeltaMLErrors{} : model.evaluate(test);
        curve << cycle << "  " << data.size() << "  " << e.energy_mae_centered * 627.5094740631
              << "  " << e.force_rmse << "  " << e.spearman << "  " << max_unc << "\n";
        curve.flush();
    }

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
    std::cout << "  curcuma -activelearn mol.xyz -al_mode learn -al_cycles 6 -threads 8\n";
    std::cout << "  curcuma -opt mol.xyz -method deltaml:conf.deltaml.json\n\n";
    std::cout << "The reference must carry the same dispersion treatment as the base method,\n";
    std::cout << "otherwise the long-range dispersion term ends up in the learned difference\n";
    std::cout << "where a local descriptor cannot represent it.\n\n";
}
