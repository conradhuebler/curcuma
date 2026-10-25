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
PARAM(lambda, Double, 1e-6,
    "RKHS regulariser strength relative to the data term: the solved system is "
    "(Phi^T W Phi + lambda*s*K_MM) beta = Phi^T W y with s the mean diagonal of the data "
    "block. lambda = 1 regularises as strongly as the data constrains, small values approach "
    "interpolation. Too small gives a rippled surface that is invisible in the energy error "
    "and does not conserve energy in dynamics -- raise it if an NVE run drifts. Use the "
    "dotted form -activelearn.lambda; the flat -lambda does not route. Measured on a 90-atom "
    "conformer set: 1e-6 is best for ranking and energy conservation, 1e-2 gives the closest "
    "optimised geometries but the worst ranking. There is no single optimum -- pick for the "
    "application.",
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
    "Seed for the train/test split, the sampling dynamics (initial velocities and thermostat "
    "of every run) and the draw of the validation frames. Same seed, same input, same binary: "
    "bit-identical 'learn' runs (checked at -threads 1, 4 and 8 on a 9-atom test).",
    "Algorithm", {})

PARAM(al_cycles, Int, 5,
    "Active-learning cycles in 'learn' mode. Each cycle samples on the currently corrected "
    "surface, picks the least familiar structures, computes references for them and refits.",
    "ActiveLearning", {})
PARAM(al_md_method, String, "",
    "Method driving the sampling dynamics. Empty (default) uses the base method in cycle 0 "
    "and the corrected surface afterwards, so reference evaluations are spent on labels -- "
    "where they buy exact values -- rather than on sampling, where approximate geometries "
    "suffice. Set it explicitly (e.g. gfn2) to sample on the reference instead: that removes "
    "the self-confirmation bias, at a steep price. Measured on 90 atoms, the default four "
    "runs of 2000 steps cost about 7 s on gfnff and 11 min on gfn2 -- per cycle.",
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
PARAM(al_onthefly, Bool, false,
    "Stop each sampling run at the first structure that leaves the trained region, label it "
    "immediately and refit before continuing. Batch sampling instead runs a fixed length on "
    "a model that may already be out of its depth, and only notices afterwards. The stronger "
    "scheme, at the price of many short trajectories: each run ends where the information is.",
    "ActiveLearning", {})
PARAM(al_onthefly_threshold, Double, 0.9,
    "Extrapolation indicator that ends a sampling run in on-the-fly mode. Same scale as "
    "al_uncertainty: 0.9 corresponds to a nearest kernel centre 2.1 sigma away. Lower values "
    "stop sooner and label more often.",
    "ActiveLearning", {})
PARAM(al_initial, Int, 15,
    "Structures labelled in cycle 0, before any model exists.", "ActiveLearning", {})
// Prospective validation and convergence (Claude Generated, October 2026).
PARAM(al_validation_frames, Int, 10,
    "Frames drawn uniformly at random from each cycle's sampling dynamics (excluding the "
    "selected candidates) and labelled in addition to them. The model entering the cycle is "
    "evaluated on them BEFORE it is refitted: an unbiased error on exactly the regions the "
    "dynamics on the learned surface visits, including energy holes the training error cannot "
    "show. Afterwards they join the training data. Cycle 0 has no model, its frames only go "
    "into training. 0 switches the validation off.",
    "ActiveLearning", {})
PARAM(al_criterion, String, "both",
    "Convergence criterion: 'uncertainty' (no sampled structure above al_uncertainty, the "
    "original rule), 'validation' (validation error within al_target_energy and no hole deeper "
    "than al_target_hole) or 'both' (both at once). Must hold in al_patience consecutive cycles. "
    "'validation' and 'both' need al_validation_frames > 0. Measured on peo8 (58 atoms): the "
    "uncertainty rule alone stopped at 87 labels with the torsion distribution still wrong.",
    "ActiveLearning", {})
PARAM(al_target_energy, Double, 1.0,
    "Validation criterion: largest allowed relative (centred) energy MAE of the model on the "
    "cycle's validation frames [kcal/mol].",
    "ActiveLearning", {})
PARAM(al_target_hole, Double, 2.0,
    "Validation criterion: deepest allowed energy hole [kcal/mol], i.e. no validation frame may "
    "have model minus reference energy below -al_target_hole. Absolute energies, so a constant "
    "offset of the model counts as well.",
    "ActiveLearning", {})
PARAM(al_patience, Int, 2,
    "Consecutive cycles in which al_criterion must hold before the run stops (criteria "
    "'validation' and 'both'). The plain 'uncertainty' criterion always stops at the first cycle "
    "that meets it (the behaviour before October 2026).",
    "ActiveLearning", {})

PARAM(dataset, String, "",
    "Labelled dataset (JSON as written by -activelearn, <basename>.deltaml.dataset.json) to "
    "train on in 'train' mode instead of labelling the input structures again. The input file "
    "is still required (basename of the outputs, anchor structure) but is not labelled. Use it "
    "for hyperparameter studies: a refit then costs seconds instead of reference calculations.",
    "Basic", {})

PARAM(acsf_n_shell, Int, 6,
    "Number of radial ACSF shells (Gaussian centres on an equidistant grid from acsf_r_first "
    "to acsf_r_last). Sets the radial feature count per element; more shells resolve distances "
    "more finely at proportionally higher descriptor cost. Default 6 since October 2026 (was 8): "
    "together with acsf_radial_eta 4 it lowered the held-out relative energy MAE by "
    "0.215 +- 0.100 kcal/mol on a 90-atom conformer set (10 split seeds). Models store their "
    "own grid, so older model files are unaffected.",
    "Descriptor", {})
PARAM(acsf_r_first, Double, 0.9,
    "Centre of the first radial shell [Angstrom].", "Descriptor", {})
PARAM(acsf_r_last, Double, 5.2,
    "Centre of the last radial shell [Angstrom]; must not exceed acsf_rcut_radial.",
    "Descriptor", {})
PARAM(acsf_radial_eta, Double, 4.0,
    "Width parameter eta of the radial Gaussians exp(-eta (R - R_s)^2) [1/Angstrom^2]. Larger "
    "values give narrower shells; with too few shells for the width, gaps in the grid appear. "
    "Default 4.0 since October 2026 (was 8.0), see acsf_n_shell.",
    "Descriptor", {})
PARAM(acsf_rcut_radial, Double, 6.0,
    "Cutoff radius of the radial functions [Angstrom].", "Descriptor", {})
PARAM(acsf_rcut_angular, Double, 4.0,
    "Cutoff radius of the angular functions [Angstrom]; the angular cost grows with the square "
    "of the neighbour count inside it.",
    "Descriptor", {})
PARAM(acsf_angular_eta, Double, 0.05,
    "Radial damping eta of the angular functions [1/Angstrom^2].", "Descriptor", {})
PARAM(acsf_zeta, String, "1,4,16",
    "Comma-separated angular sharpness exponents zeta, each used with lambda = +1 and -1. "
    "Example: '1,2,4,8,16'. Each entry adds two angular features per element pair.",
    "Descriptor", {})

// Energy model (Claude Generated, October 2026). Defaults from the October 2026 campaign
// (docs/DELTA_ML.md, "Held-out accuracy"); to change one, edit it here and in the matching
// member initialiser of ActiveLearn (m_energy_model*), nothing else reads them.
PARAM(energy_model, Bool, true,
    "After the main (joint energy+gradient) fit, also fit an energy-only model "
    "(use_forces = false, same data, split and descriptor) and write it as "
    "<basename>.deltaml.energy.json next to the main model, with \"role\": \"energy\". Measured "
    "on a 90-atom conformer set (GFN-FF -> GFN2, 5 split seeds): the joint model's held-out "
    "relative energy MAE is worse than the base method's (+0.26 kcal/mol) while its gradient "
    "is 4x better; the energy-only fit (lambda 1e-4) has a smaller relative energy error than "
    "the base (-0.19 +- 0.05 kcal/mol, 10 seeds) but a rank correlation no better than GFN-FF, "
    "and forces no better than the base. Use the energy model for relative energies (-sp, "
    "ConfScan energies), the main model for -opt/-md. Only in al_mode train (also with -dataset); not in learn mode "
    "and not in the center_sweep.",
    "Energy model", {})
PARAM(energy_model_n_centers, Int, 0,
    "Kernel centres of the energy-only model; 0 uses n_centers.",
    "Energy model", {})
PARAM(energy_model_lambda, Double, 1e-4,
    "Regulariser of the energy-only model, same meaning as lambda. 1e-4 measured best on a "
    "90-atom conformer set (held-out relative energy MAE 0.623 +- 0.028 kcal/mol against "
    "0.695 +- 0.052 at 1e-6, 10 split seeds); 1e-8 and below is clearly worse.",
    "Energy model", {})
PARAM(energy_model_kernel_sigma, Double, 0.0,
    "Kernel width of the energy-only model (0 = median heuristic), same meaning as kernel_sigma.",
    "Energy model", {})

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

    /*! \brief Structures a stopped sampling run flagged as outside the trained region. */
    std::vector<Molecule> m_triggered;

    /**
     * @brief Pick the structures worth spending reference calculations on.
     *
     * Least-familiar-first, but with a diversity term: maximising the uncertainty alone
     * piles every pick onto one pathological region, which buys almost no new information
     * for the same cost. Selection is greedy farthest-point in descriptor space among the
     * candidates that exceed the uncertainty threshold.
     *
     * @param max_uncertainty out, largest indicator seen; the convergence signal
     * @param picked_indices optional out, index of each pick in `candidates`
     * @param ignore_threshold true: every candidate is eligible, not only those above
     *        al_uncertainty (used when the indicator is satisfied but the validation is not)
     */
    std::vector<Molecule> selectCandidates(const std::vector<Molecule>& candidates,
        const curcuma::ml::DeltaMLModel& model, bool have_model, int wanted,
        double& max_uncertainty, std::vector<std::size_t>* picked_indices = nullptr,
        bool ignore_threshold = false);

    /**
     * @brief Draw up to al_validation_frames frames uniformly at random, without replacement.
     *
     * Frames listed in `excluded` (the selected candidates) are skipped, so no structure is
     * labelled twice. The draw is seeded from (seed, cycle) and therefore reproducible.
     * Claude Generated (October 2026)
     */
    std::vector<Molecule> drawValidationFrames(const std::vector<Molecule>& frames,
        const std::vector<std::size_t>& excluded, int cycle,
        std::vector<std::size_t>* drawn = nullptr) const;

    /*! \brief Training options of the main (joint) fit with M = n_centers centres. */
    curcuma::ml::DeltaMLTrainingOptions trainingOptions(int n_centers) const;

    /*! \brief Fit one model with the main-fit options and log its held-out performance. */
    bool fitAndReport(const curcuma::ml::MLDataset& train, const curcuma::ml::MLDataset& test,
        int n_centers, curcuma::ml::DeltaMLModel& model, bool log_details);

    /**
     * @brief Fit one model with explicit options.
     * @param label  name of the model in the log ("joint model", "energy model")
     * @param log_heldout false: skip the held-out block (the caller prints it later)
     * Claude Generated (October 2026)
     */
    bool fitAndReport(const curcuma::ml::MLDataset& train, const curcuma::ml::MLDataset& test,
        const curcuma::ml::DeltaMLTrainingOptions& opt, curcuma::ml::DeltaMLModel& model,
        bool log_details, const std::string& label, bool log_heldout);

    /*! \brief Path of the energy-only model next to the main model path. */
    static std::string energyModelPath(const std::string& main_model_path);

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
    double m_lambda = 1e-6;
    double m_jitter_relative = 1e-8;
    double m_sigma_energy = 1e-4;
    double m_sigma_force = 1e-3;
    bool m_use_forces = true;
    unsigned m_seed = 42;
    bool m_write_dataset = true;
    std::string m_dataset_path;             ///< labelled dataset to train on instead of labelling
    curcuma::ml::ACSFSettings m_acsf;       ///< descriptor grid for every fit and selection
    bool m_acsf_valid = true;               ///< false if acsf_zeta did not parse

    // Energy-only second model (Claude Generated, October 2026). Runtime values come from the
    // PARAM block above; keep these initialisers equal to it when the defaults are updated.
    bool m_energy_model = true;
    int m_energy_model_n_centers = 0;       ///< 0 = m_n_centers
    double m_energy_model_lambda = 1e-4;
    double m_energy_model_kernel_sigma = 0.0; ///< 0 = median heuristic

    int m_al_cycles = 5;
    std::string m_al_md_method;
    double m_al_temperature = 500.0;
    double m_al_md_time = 2000.0;
    int m_al_md_runs = 4;
    int m_al_dump_frequency = 50;
    int m_al_candidates = 15;
    double m_al_uncertainty = 0.05;
    int m_al_initial = 15;
    bool m_al_onthefly = false;
    double m_al_onthefly_threshold = 0.9;
    // Prospective validation and convergence (Claude Generated, October 2026); values from the PARAM block.
    int m_al_validation_frames = 10;
    std::string m_al_criterion = "both";
    double m_al_target_energy = 1.0; ///< kcal/mol
    double m_al_target_hole = 2.0;   ///< kcal/mol
    int m_al_patience = 2;

    /* Anchor structure and, more importantly, the bond graph taken from it. Only the graph
       is pinned for later structures; everything geometry dependent is recomputed where it
       is needed. */
    Molecule m_anchor;
    std::vector<std::pair<int, int>> m_anchor_bonds;

    /*! \brief Determine the bond graph of the anchor once, with the base method itself. */
    void captureAnchorBonds();
};
