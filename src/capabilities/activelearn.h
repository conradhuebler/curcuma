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
PARAM(lambda, Double, 1e-8,
    "Weight of the RKHS regulariser. Too small interpolates label noise into a rippled "
    "surface, which is invisible in the energy error but destroys molecular dynamics.",
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
    double m_sigma_energy = 1e-4;
    double m_sigma_force = 1e-3;
    bool m_use_forces = true;
    unsigned m_seed = 42;
    bool m_write_dataset = true;

    /* Anchor structure that fixes the topology for every calculator in this run. */
    Molecule m_anchor;
};
