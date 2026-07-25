/*
 * <Conformer proposals from an energy-decomposed ensemble>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
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
 * Claude Generated (Jul 2026)
 */

#pragma once

#include <map>
#include <string>
#include <vector>

#include "json.hpp"

#include "src/capabilities/curcumamethod.h"
#include "src/capabilities/torsion_space.h"
#include "src/core/config_manager.h"
#include "src/core/molecule.h"
#include "src/core/parameter_macros.h"
#include "src/core/parameter_registry.h"

using namespace curcuma;

/**
 * @brief Learn which structural features carry energy in a conformer ensemble, and propose new
 *        conformers by recombining the favourable ones.
 *
 * THE IDEA
 * A conformer search (ConfSearch/MTD) produces an ensemble of optimised structures. Every member is
 * a point in the discrete torsion-state space (see TorsionSpace) and carries a force-field energy
 * that GFN-FF can decompose per term (bond/angle/torsion/repulsion/dispersion/Coulomb/H-bond/
 * halogen-bond). That decomposition is thrown away today -- the search only ever uses total energies.
 *
 * MATCHED PAIRS (this stage)
 * The contribution of ONE torsion is measured, not fitted: take all pairs of ensemble members whose
 * state vectors differ in exactly one torsion (Hamming distance 1). Their per-term energy differences
 * are the contribution of that single state change in a real molecular environment:
 *
 *     torsion C4-C7, state 180 deg -> 60 deg:  dE_total -8.2 kJ/mol
 *                                              = Torsion +4.1, HBond -12.3, Dispersion +0.1, ...
 *
 * Several pairs per transition give the SPREAD, which measures how much the contribution depends on
 * the rest of the molecule -- i.e. how strongly this torsion is coupled to the others. No model fit,
 * no extra electronic-structure calculation beyond one single point per ensemble member.
 *
 * WHAT IT CANNOT DO (honest scope)
 * States the ensemble never visited stay invisible: this recombines, it does not extrapolate. A
 * torsion that appears in only one state in the whole ensemble yields no contrast and no information.
 * Relaxed torsion scans are the tool for those two cases and are a planned add-on, not the basis.
 *
 * See docs/CONFSEARCH_PROPOSALS.md.
 */
class ConfGen : public CurcumaMethod {
public:
    explicit ConfGen(const json& controller, bool silent);
    ~ConfGen() = default;

    void setFile(const std::string& filename) override;
    bool Initialise() override { return true; }
    void start() override;

    void printHelp() const override
    {
        std::cout << "Usage: curcuma -confgen ensemble.xyz [parameters]\n\n";
        ParameterRegistry::getInstance().printHelp("confgen");
        std::cout << "\nThe input is a multi-structure XYZ of ALREADY OPTIMISED conformers,\n"
                     "e.g. <basename>.cumulative.opt.accepted.xyz from a ConfSearch run.\n";
    }

private:
    /// Per-structure data collected in the analysis pass.
    struct Frame {
        Molecule molecule;
        std::vector<double> angles;   ///< dihedral of every torsion, degrees
        std::vector<int> states;      ///< discrete state per torsion
        double energy = 0.0;          ///< total energy, Hartree
        std::map<std::string, double> terms; ///< energy decomposition, Hartree
        bool valid = false;           ///< topology matches the reference structure
    };

    /// Aggregated result for one state transition of one torsion.
    struct Transition {
        int torsion = -1;
        int state_from = -1, state_to = -1; ///< canonical order: state_from < state_to
        int pairs = 0;
        int distinct_from = 0, distinct_to = 0; ///< distinct structures on each side of the change
        double d_total_mean = 0.0, d_total_min = 0.0, d_total_max = 0.0; ///< Hartree
        double d_total_sd = 0.0;                                        ///< Hartree, sample std dev
        std::map<std::string, double> d_terms_mean;                     ///< Hartree, per term
        /// Full range. For many pairs the standard deviation (d_total_sd) is the honest measure --
        /// max-min grows with the sample size and is driven by single outliers.
        double range() const { return d_total_max - d_total_min; }
    };

    /* Read the ensemble, run one single point per member (shared calculator -> identical topology
     * and parameters for every frame, so the term differences are comparable), and fill m_frames. */
    bool analyseEnsemble();

    /* All Hamming-1 pairs -> per-transition term statistics. */
    std::vector<Transition> matchedPairs() const;

    /**
     * @brief One measured coupling between two torsions (double-mutant cycle).
     *
     * Four ensemble members that form a RECTANGLE in state space -- (a,c), (b,c), (a,d), (b,d) with
     * every other torsion identical -- give the coupling without any model fit:
     *
     *     J = [E(b,d) - E(a,d)] - [E(b,c) - E(a,c)]
     *
     * J = 0 means the two state changes are independent (their effects simply add); J != 0 is exactly
     * the non-additivity. This is the double-mutant cycle of protein biochemistry (Carter/Fersht/
     * Horovitz 1984/1990), applied to torsions instead of side chains.
     */
    struct Coupling {
        int torsion_a = -1, torsion_b = -1;
        int a_from = -1, a_to = -1, b_from = -1, b_to = -1;
        int cycles = 0;             ///< number of independent rectangles found
        double j_mean = 0.0;        ///< Hartree
        double j_sd = 0.0;          ///< Hartree
        std::map<std::string, double> j_terms_mean; ///< per-term coupling, Hartree
    };

    /* Search the ensemble for double-mutant cycles and average their coupling. */
    std::vector<Coupling> doubleMutantCycles() const;

    /// Result of one cross-validated model fit.
    struct ModelFit {
        std::string name;
        int columns = 0;      ///< design-matrix columns offered
        int rank = 0;         ///< linearly independent ones (what the ensemble can actually resolve)
        double rmse_cv = 0.0; ///< kJ/mol, out-of-sample
        double mae_cv = 0.0;  ///< kJ/mol, out-of-sample median absolute error (outlier-robust)
        double rmse_in = 0.0; ///< kJ/mol, in-sample (for reference only -- always improves)
        double r2_cv = 0.0;   ///< fraction of the energy variance explained out of sample
    };

    /**
     * @brief Fit E(s) with indicator variables and score it OUT OF SAMPLE.
     *
     * level 0: constant only (the null model -- its error IS the energy spread)
     * level 1: + one coefficient per torsion state    E ~ c + sum_i h_i(s_i)
     * level 2: + one per state pair of two torsions   E ~ ... + sum_ij J_ij(s_i,s_j)
     *
     * k-fold cross-validation is not optional here: level 2 has many more parameters and would win
     * any in-sample comparison by construction. Only a better prediction on data the fit has not seen
     * shows that the couplings carry real information.
     */
    ModelFit fitModel(int level) const;

    /* Per torsion and state: population, lowest relative energy, Boltzmann-weighted mean. */
    void writeStateStatistics(const std::string& path) const;
    void writeFrameTable(const std::string& path) const;
    void writeTransitionTable(const std::string& path, const std::vector<Transition>& transitions) const;
    void reportTransitions(const std::vector<Transition>& transitions) const;
    void reportCouplings(const std::vector<Coupling>& couplings) const;
    void reportModelComparison(const std::vector<ModelFit>& fits) const;
    void writeCouplingTable(const std::string& path, const std::vector<Coupling>& couplings) const;

    /// Torsions with at least two populated states -- the only ones that carry information.
    std::vector<int> informativeTorsions() const;

    /// One generated candidate: a state vector that does not occur in the ensemble.
    struct Proposal {
        std::vector<int> states;   ///< target state vector
        int template_frame = -1;   ///< ensemble member the geometry was built from
        int distance = 0;          ///< Hamming distance to that template
        double predicted = 0.0;    ///< additive-model estimate, kJ/mol (ORDERING ONLY, see below)
        Molecule geometry;         ///< built structure (before optimisation)
        // filled after the optimisation
        bool optimised = false;
        double energy = 0.0;              ///< Hartree
        std::vector<int> states_after;    ///< state vector of the optimised structure
        double min_rmsd_to_ensemble = 0.0;///< Angstrom, best-fit RMSD to the closest input structure
        bool topology_ok = false;         ///< optimised structure still has the reference bond topology
        bool is_new = false;              ///< survived the topology AND novelty checks
    };

    /**
     * @brief Enumerate state vectors that the ensemble does NOT contain, build them, optimise them
     *        and check whether they are new conformers.
     *
     * The cross-validated model comparison showed that a torsion-state energy model explains only
     * ~15 % of the energy variation, so the model is used for ONE thing only: deciding which untried
     * combinations to build first. The force field decides everything else -- every proposal is
     * optimised and compared against the input ensemble, so a bad proposal costs one optimisation and
     * can never produce a wrong result. That asymmetry is what makes generation worthwhile even
     * though the model is weak.
     */
    std::vector<Proposal> generateProposals() const;
    void optimiseProposals(std::vector<Proposal>& proposals) const;

    /**
     * @brief Re-optimise the template structures to get a comparable energy reference.
     *
     * Proposals are optimised; the input ensemble may not be (it can come from another method, or
     * another optimiser setting). Comparing the two directly makes the optimisation gain look like a
     * discovery -- measured: a proposal appeared "106 kJ/mol below the ensemble minimum" purely
     * because the input structures had never been optimised with this method. Re-optimising the
     * templates costs a handful of optimisations and makes the comparison like-for-like; a large gain
     * is reported as the warning it is.
     */
    double referenceEnergyOptimised(double& worst_gain_kJ) const;
    void reportProposals(const std::vector<Proposal>& proposals, const std::string& base) const;

    /// Additive-model coefficients (intercept + per torsion state), kJ/mol. Ordering heuristic only.
    std::vector<std::vector<double>> additiveCoefficients() const;

    /// Cheapest possible sanity filter on a built geometry: any atom pair far below bonding distance.
    bool hasClash(const Molecule& mol, double factor) const;

    /**
     * @brief Sorted list of bonded atom pairs, with an EXPLICIT covalent-radius factor.
     *
     * Deliberately not Molecule::DistanceMatrix(): its default scaling of 1.5 is generous enough to
     * call a compressed 1-3 contact a bond. Measured on a 114-atom ensemble: at 1.5 it flagged 45 of
     * 46 optimised proposals as "topology changed" while an independent check at 1.3 found 44 of them
     * bit-identical to the reference -- i.e. the test, not the structures, was wrong. The factor is
     * exposed as -topology_factor so the criterion is visible rather than hidden in a default.
     */
    std::vector<std::pair<int, int>> topologyFingerprint(const Molecule& mol, double factor) const;

    StringList MethodName() const override { return { "ConfGen" }; }
    void ReadControlFile() override { }
    void LoadControlJson() override;
    nlohmann::json WriteRestartInformation() override { return nlohmann::json::object(); }
    bool LoadRestartInformation() override { return false; }

    ConfigManager m_config;
    std::string m_method = "gfnff";
    int m_charge = 0, m_spin = 0, m_threads = 1;
    double m_state_tolerance = 40.0;
    double m_temperature = 298.15;
    int m_min_pairs = 1;
    double m_report_threshold = 1.0;
    int m_cv_folds = 5;
    bool m_couplings = true;
    bool m_generate = false;
    int m_max_proposals = 50, m_proposal_templates = 5, m_proposal_depth = 2;
    double m_clash_factor = 1.2, m_new_rmsd = 1.0, m_topology_factor = 1.3;

    std::vector<TorsionSpace::Torsion> m_torsions;
    std::vector<std::vector<double>> m_state_centres; ///< per torsion
    std::vector<Frame> m_frames;
    std::vector<std::string> m_term_names; ///< decomposition keys actually present, stable order
    double m_reference_energy = 0.0;       ///< lowest total energy in the ensemble (Hartree)

    // vvvvvvvvvvvv PARAMETER DEFINITION BLOCK vvvvvvvvvvvv
    // Every Double MUST carry a decimal-point literal (see the note in confsearch.h).
    BEGIN_PARAMETER_DEFINITION(confgen)

    PARAM(method, String, "gfnff", "Energy method used for the single point behind every ensemble member. Only force fields provide a per-term decomposition; gfnff is the intended method.", "Methods", {})
    PARAM(charge, Int, 0, "Total charge of the system.", "System", {})
    PARAM(spin, Int, 0, "Total spin of the system.", "System", {})
    PARAM(threads, Int, 1, "Threads handed to the energy method.", "System", {})

    PARAM(state_tolerance, Double, 40.0, "Angular tolerance in degrees for grouping observed torsion values into rotamer states. Larger values merge neighbouring basins, smaller values split noisy ones.", "Analysis", {})
    PARAM(temperature, Double, 298.15, "Temperature in Kelvin for the Boltzmann-weighted state statistics.", "Analysis", {})
    PARAM(min_pairs, Int, 1, "Minimum number of matched pairs required before a state transition is reported.", "Analysis", {})
    PARAM(generate, Bool, false, "Generate new conformer proposals: enumerate state vectors that the ensemble does not contain, build them, optimise them and report which ones are genuinely new. Off by default because it runs one geometry optimisation per proposal.", "Generation", {})
    PARAM(max_proposals, Int, 50, "Maximum number of proposals built and optimised, ordered by the additive-model estimate.", "Generation", {})
    PARAM(proposal_templates, Int, 5, "Number of lowest-energy ensemble members used as geometric templates.", "Generation", {})
    PARAM(proposal_depth, Int, 2, "Maximum number of torsions changed simultaneously relative to a template (Hamming distance).", "Generation", {})
    PARAM(clash_factor, Double, 1.2, "A built structure is rejected when a non-bonded atom pair comes closer than this factor times the sum of their covalent radii. The default is deliberately close to the BOND-DETECTION criterion (~1.3): a built structure that puts two atoms inside bonding distance makes the force field derive a new bond, and the optimisation then relaxes into a different molecule, not a conformer.", "Generation", {})
    PARAM(topology_factor, Double, 1.3, "Covalent-radius factor for the topology check of optimised proposals. A proposal whose bond list differs from the reference is a reaction product, not a conformer, and is rejected. Lower than Molecule's default 1.5, which counts compressed 1-3 contacts as bonds.", "Generation", {})
    PARAM(new_rmsd, Double, 1.0, "Best-fit RMSD in Angstrom above which an optimised proposal counts as a new conformer.", "Generation", {})
    PARAM(cv_folds, Int, 5, "Number of cross-validation folds for the model comparison (additive vs. additive+couplings). Values below 2 disable the comparison.", "Analysis", {})
    PARAM(couplings, Bool, true, "Measure torsion-torsion couplings from double-mutant cycles (four ensemble members forming a rectangle in state space) and run the model comparison.", "Analysis", {})
    PARAM(report_threshold, Double, 1.0, "Only transitions whose mean total energy difference exceeds this many kJ/mol are printed in the summary. All of them are written to the CSV.", "Analysis", {})

    END_PARAMETER_DEFINITION
    // ^^^^^^^^^^^^ PARAMETER DEFINITION BLOCK ^^^^^^^^^^^^
};
