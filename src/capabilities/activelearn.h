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

#pragma once

#include "curcumamethod.h"

#include "src/core/molecule.h"
#include "src/core/ml/krr_model.h"
#include "src/core/ml/ml_dataset.h"
#include "src/core/parameter_macros.h"

#include <string>
#include <vector>

BEGIN_PARAMETER_DEFINITION(activelearn)
PARAM(al_mode, String, "train",
    "What to do: 'train' fits a model to a given structure set, 'test' evaluates an "
    "existing model against reference calculations, 'learn' runs the active-learning "
    "cycle that generates its own structures.",
    "Basic", { "mode" })
PARAM(ref_method, String, "gfn2",
    "Reference method providing the labels, i.e. the level the base method is lifted to. "
    "IMPORTANT: it must carry the same dispersion treatment as the base method. Pairing a "
    "dispersion-free reference with a D4-corrected base pushes the whole long-range "
    "dispersion term into the learned difference, where a local descriptor cannot represent it.",
    "Basic", {})
PARAM(model, String, "",
    "Model file. Output path in 'train' and 'learn', input path in 'test'. Empty derives it "
    "from the input basename.",
    "Basic", {})
PARAM(max_structures, Int, 0,
    "Use at most this many structures from the input (0 = all). Useful to check how the "
    "error falls with training-set size before committing to a full run.",
    "Basic", {})
PARAM(test_fraction, Double, 0.2,
    "Share of structures held out for validation. The split is random rather than a "
    "contiguous tail, because neighbouring frames of a trajectory or conformer file are "
    "correlated and a tail split would flatter the model.",
    "Basic", {})

PARAM(n_centers, Int, 1000,
    "Number of kernel centres M. This is the accuracy/cost lever: a prediction costs O(M) "
    "regardless of dataset size. Whether a large M pays off depends on how expensive the "
    "reference is relative to the base method, so measure rather than guess.",
    "Algorithm", {})
PARAM(kernel_sigma, Double, 0.0,
    "Gaussian kernel width in standardised descriptor space (0 = median heuristic). Strong "
    "lever: too narrow and the model cannot interpolate, too wide and it oversmooths.",
    "Algorithm", {})
PARAM(lambda, Double, 1e-4,
    "RKHS regulariser strength relative to the data term: the solved system is "
    "(Phi^T W Phi + lambda*s*K_MM) beta = Phi^T W y with s the mean diagonal of the data "
    "block. lambda = 1 regularises as strongly as the data constrains, small values approach "
    "interpolation. Too small gives a rippled surface that is invisible in the energy error "
    "and does not conserve energy in dynamics -- raise it if an NVE run drifts. Use the "
    "dotted form -activelearn.lambda; the flat -lambda does not route. Measured on a "
    "90-atom conformer set: 1e-8 gives the best energies but an NVE drift 3.7x the base "
    "method, 1e-2 matches the base drift but loses most of the accuracy, 1e-4 sits between.",
    "Algorithm", {})
PARAM(jitter_relative, Double, 1e-8,
    "Ridge added to the normal equations relative to their diagonal scale. Accumulating "
    "normal equations squares the condition number, and the centre kernel is near-singular "
    "once many centres crowd together, so without this the factorisation can succeed and "
    "still return a useless solution. Symptom: accuracy gets worse as n_centers grows.",
    "Algorithm", {})
PARAM(sigma_energy, Double, 1e-4,
    "Assumed energy noise [Eh]; sets the weight of energy observations relative to gradients.",
    "Algorithm", {})
PARAM(sigma_force, Double, 1e-3,
    "Assumed gradient noise [Eh/Angstrom]; sets the weight of gradient observations.",
    "Algorithm", {})
PARAM(use_forces, Bool, true,
    "Train on gradients as well as energies. One reference calculation yields a single "
    "energy but 3N gradient components, which is what makes a few hundred points enough. "
    "Required for optimisation and dynamics on the corrected surface.",
    "Algorithm", {})
PARAM(seed, Int, 42,
    "Seed for the train/test split and sampling heuristics, fixed so runs are reproducible.",
    "Algorithm", {})

PARAM(al_cycles, Int, 5,
    "Active-learning cycles in 'learn' mode. Each cycle samples on the currently corrected "
    "surface, picks the least familiar structures, computes references for them and refits.",
    "ActiveLearning", {})
PARAM(al_temperature, Double, 500.0,
    "Sampling temperature [K]. Deliberately above the temperature of interest: sampling "
    "happens on the surface being learned, so cycle 0 only ever visits what the base method "
    "considers relevant. Heat and metadynamics mitigate that bias, they do not remove it.",
    "ActiveLearning", {})
PARAM(al_md_time, Double, 2000.0,
    "Simulated time per sampling run [fs].", "ActiveLearning", {})
PARAM(al_md_runs, Int, 4,
    "Independent sampling runs per cycle, started from different structures. Several short "
    "runs cover more of the surface than one long one.",
    "ActiveLearning", {})
PARAM(al_dump_frequency, Int, 50,
    "Keep every Nth molecular-dynamics frame as a candidate.", "ActiveLearning", {})
PARAM(al_candidates, Int, 15,
    "Reference calculations per cycle, i.e. how many candidates are labelled.",
    "ActiveLearning", {})
PARAM(al_uncertainty, Double, 0.5,
    "Convergence threshold on the extrapolation indicator 1 - max_m k(x,z_m). The cycle "
    "stops once no sampled structure exceeds it. Mind the scale: with a Gaussian kernel the "
    "indicator saturates quickly, since a nearest centre at distance d gives "
    "1 - exp(-d^2/2 sigma^2), so 0.5 corresponds to d = 1.18 sigma and 0.9 already to "
    "d = 2.1 sigma. Values below about 0.2 are effectively unreachable. And it measures "
    "descriptor-space distance, not error: it certifies coverage, never accuracy.",
    "ActiveLearning", {})
PARAM(al_initial, Int, 15,
    "Structures labelled in cycle 0, before any model exists.", "ActiveLearning", {})

PARAM(center_sweep, String, "",
    "Comma-separated list of M values to fit and report in one run, e.g. '200,500,1000,2000'. "
    "Produces the accuracy-over-cost curve that decides how sparse the model can be.",
    "Output", {})
PARAM(write_dataset, Bool, true,
    "Write the labelled structures as JSON (and XYZ for inspection) next to the model, so a "
    "refit with different hyperparameters needs no new reference calculations.",
    "Output", {})
END_PARAMETER_DEFINITION

/**
 * @brief Builds delta-ML models: offline fitting today, active learning on top.
 *
 * Offline mode ('train') takes an existing structure set -- a conformer file, a trajectory --
 * computes base and reference energies and gradients for every structure, fits the model and
 * reports how much it actually improved over the bare base method. That contrast is the point:
 * an error without its baseline says nothing.
 *
 * Claude Generated (July 2026)
 */
class ActiveLearn : public CurcumaMethod {
public:
    ActiveLearn(const json& controller, bool silent);
    ~ActiveLearn() = default;

    void start() override;
    void printHelp() const override;

    /*! \brief Set the input structure file; also derives the basename used for outputs. */
    inline void setFileName(const std::string& filename) { setFile(filename); }

private:
    /*! \brief Read structures from the input file, honouring max_structures. */
    bool loadStructures(std::vector<Molecule>& molecules);

    /**
     * @brief Compute base and reference energies and gradients for every structure.
     *
     * All calculators are anchored on the first structure so that methods which freeze a
     * topology at setup (GFN-FF) use an identical bond list everywhere. Without that the
     * difference picks up setup noise that no amount of training can remove.
     */
    bool labelStructures(const std::vector<Molecule>& molecules, curcuma::ml::MLDataset& data);

    void runTrain();
    void runTest();
    void runLearn();

    /**
     * @brief Sample candidate structures with molecular dynamics.
     * @param md_method method driving the dynamics: the base method in cycle 0, the
     *        corrected surface afterwards
     * @param seeds structures the runs start from
     * @param cycle current cycle, used for output naming and seeding
     */
    std::vector<Molecule> sampleStructures(
        const std::string& md_method, const std::vector<Molecule>& seeds, int cycle);

    /**
     * @brief Pick the structures worth spending reference calculations on.
     *
     * Least-familiar-first, but with a diversity term: maximising the uncertainty alone
     * piles every pick onto one pathological region, which buys almost no new information
     * for the same cost. Selection is greedy farthest-point in descriptor space among the
     * candidates that exceed the uncertainty threshold.
     *
     * @param max_uncertainty out, largest indicator seen; the convergence signal
     */
    std::vector<Molecule> selectCandidates(const std::vector<Molecule>& candidates,
        const curcuma::ml::DeltaMLModel& model, bool have_model, int wanted,
        double& max_uncertainty);

    /*! \brief Fit one model and log its held-out performance. */
    bool fitAndReport(const curcuma::ml::MLDataset& train, const curcuma::ml::MLDataset& test,
        int n_centers, curcuma::ml::DeltaMLModel& model, bool log_details);

    void logErrors(const std::string& label, const curcuma::ml::DeltaMLErrors& e) const;

    /* Reads controller values into members. */
    void LoadControlJson() override;

    /* CurcumaMethod plumbing; this capability keeps no restart state. */
    nlohmann::json WriteRestartInformation() override { return json{}; }
    bool LoadRestartInformation() override { return true; }
    StringList MethodName() const override { return { "activelearn" }; }
    void ReadControlFile() override {}

    std::string m_mode = "train";
    std::string m_base_method = "gfnff";
    std::string m_ref_method = "gfn2";
    std::string m_model_path;
    std::vector<int> m_center_sweep;

    int m_max_structures = 0;
    double m_test_fraction = 0.2;
    int m_n_centers = 1000;
    double m_kernel_sigma = 0.0;
    double m_lambda = 1e-8;
    double m_jitter_relative = 1e-8;
    double m_sigma_energy = 1e-4;
    double m_sigma_force = 1e-3;
    bool m_use_forces = true;
    unsigned m_seed = 42;
    bool m_write_dataset = true;

    int m_al_cycles = 5;
    double m_al_temperature = 500.0;
    double m_al_md_time = 2000.0;
    int m_al_md_runs = 4;
    int m_al_dump_frequency = 50;
    int m_al_candidates = 15;
    double m_al_uncertainty = 0.05;
    int m_al_initial = 15;

    /* Anchor structure that fixes the topology for every calculator in this run. */
    Molecule m_anchor;
};
