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
#include <memory>
#include <string>
#include <vector>

#include "json.hpp"

#include "src/capabilities/curcumamethod.h"
#include "src/capabilities/torsion_space.h"
#include "src/core/config_manager.h"
#include "src/core/energycalculator.h"
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
        /// Presence (1) or absence (0) of every non-covalent contact in m_nci_pairs. Second
        /// descriptor next to the torsion states -- see buildNCISpace().
        std::vector<int> nci;
        std::vector<double> charges;  ///< partial charges of the single point (for the NCI detection)
        double energy = 0.0;          ///< total energy, Hartree
        std::map<std::string, double> terms; ///< energy decomposition, Hartree
        bool valid = false;           ///< topology matches the reference structure
        /** Claude Generated (Aug 2026): this frame comes from the file the run was called with (the
         *  current cycle) rather than from -analysis_file. Only such frames may serve as geometric
         *  templates: building around a structure of an earlier cycle would re-explore ground the
         *  search has already covered. The DESCRIPTION uses every frame -- that is the point of the
         *  separation, since the state and contact statistics need as many structures as possible. */
        bool from_input = true;
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

    /**
     * @brief One non-covalent interaction of the NCI pattern -- the second descriptor of a conformer.
     *
     * MOTIVATION (Jul 2026). The torsion-state vector describes what the covalent skeleton does and
     * nothing else. Measured on a 142-structure peptide ensemble, that is the wrong half of the
     * physics: the per-term attribution of a torsion change is dominated by the H-bond term in 8 of
     * 16 cases and by Coulomb in 4, by the torsion term in NONE -- and the reference structure of an
     * independent search shares its torsion-state vector with two ensemble members while lying 2.7 A
     * away from them. A description that cannot distinguish those structures cannot be the basis for
     * proposing new ones. The NCI pattern adds exactly the missing half: which non-covalent contacts
     * a conformer forms. It is built with the same discrete logic as the torsion states (a binary
     * presence vector over a global list), so every tool downstream -- Hamming distance, matched
     * pairs, the cross-validated model -- applies to it unchanged.
     */
    struct NCIContact {
        enum Kind { HBond,   ///< D-H...A, directional
            XBond,           ///< C-X...B halogen bond, directional
            Ionic,           ///< non-bonded heavy pair with a large attractive charge product
            Contact };       ///< remaining close heavy-atom pair (dispersion/repulsion)
        Kind kind = Contact;
        int first = -1, second = -1; ///< heavy atoms (donor/acceptor, X/B, or the contact pair)
        int hydrogen = -1;           ///< bridging H for kind == HBond, else -1
        std::string label() const;   ///< e.g. "HB N38-H90...O12"
    };

    /// Detect the non-covalent interactions of ONE structure (geometry + partial charges).
    std::vector<NCIContact> detectNCI(const Molecule& mol, const std::vector<double>& charges) const;

    /// Union of all contacts observed anywhere in the ensemble -> per-frame presence vector.
    void buildNCISpace();

    /// Which energy term drives the spread of the ensemble: Cov(E_term, E_total) / Var(E_total).
    void reportTermVariance() const;

    /// Population, pattern count and Hamming statistics of the NCI space.
    void reportNCISpace() const;
    void writeNCITable(const std::string& path) const;

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
        /// False when every cross-validation fold had fewer training rows than the model has
        /// parameters. Such a fit was never tested; reporting its (zero) error as a score would
        /// invert the meaning of the number.
        bool evaluated = true;
    };

    /**
     * @brief Fit E(s) with indicator variables and score it OUT OF SAMPLE.
     *
     * level 0: constant only (the null model -- its error IS the energy spread)
     * level 1: + one coefficient per torsion state    E ~ c + sum_i h_i(s_i)
     * level 2: + one per state pair of two torsions   E ~ ... + sum_ij J_ij(s_i,s_j)
     * level 3: one coefficient per NCI contact        E ~ c + sum_k g_k n_k   (torsions NOT used)
     * level 4: torsion states AND NCI contacts        E ~ c + sum_i h_i + sum_k g_k
     *
     * Levels 3 and 4 answer the question the torsion-only analysis cannot: whether the energy of a
     * conformer is carried by its covalent skeleton or by its non-covalent pattern. Comparing 1 and 3
     * on the same folds is a like-for-like test of the two descriptions.
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
        /** Claude Generated (Aug 2026): smallest Hamming distance to ANY structure the run has
         *  already seen or already tried -- the coverage half of the selection. Ranking a proposal
         *  by its predicted energy alone uses the one quantity that three measurements found to be a
         *  poor predictor; ranking it by novelty alone throws away that the terms do relate to the
         *  final energy. Both enter the score, see PARAM proposal_ranking. */
        int novelty = 0;
        /**
         * Claude Generated (Aug 2026): an NCI move instead of a torsion move. Entries are
         * (index into m_nci_pairs, desired presence 0/1) relative to the template. Empty for the
         * torsion proposals; non-empty marks a proposal that is built by DISTANCE restraints on the
         * interacting atoms rather than by driving dihedrals. Motivated by the measurement that the
         * hydrogen-bond term carries +58 % of the ensemble's energy spread and the torsion term
         * -7 %: the move set has to act on the description that actually distinguishes conformers.
         */
        std::vector<std::pair<int, int>> nci_targets;
        bool nci_move() const { return !nci_targets.empty(); }
        std::string nci_label;     ///< human-readable move, e.g. "break HB 38-90...12, form HB 5-61...46"
        double predicted = 0.0;    ///< additive-model estimate, kJ/mol (ORDERING ONLY, see below)
        Molecule geometry;         ///< built structure (before optimisation)
        bool restrained_build = false; ///< rigid build clashed; geometry came from the restrained build
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
     * @brief Build a proposal by DRIVING its torsions instead of setting them rigidly (P0).
     *
     * Rigidly setting a torsion rotates a whole fragment on a frozen template and drops it wherever
     * the template happens to have atoms -- on a compact molecule that destroyed 72 % of all
     * proposals. Here the clash is never created: the optimisation starts from the clash-free
     * template and each target torsion carries a harmonic restraint towards its target state
     * (restraint_force, Eh/rad^2), so the torsions turn while the rest of the molecule relaxes out of
     * the way. The restraints are released afterwards -- optimiseProposals() runs the normal, free
     * optimisation on the result, so the reported energy is never a restrained one.
     *
     * @param p      proposal (target state vector + template)
     * @param driven output geometry, valid only when the function returns true
     * @return false when the restrained optimisation failed or left a torsion far from its target
     */
    /**
     * @param start  optional starting geometry. Default (nullptr) = the template, which makes the
     *               restrained optimisation perform the rotation itself. That is right for one or
     *               two torsions and WRONG for many: driving 29 dihedrals at once from the template
     *               stalls at 74 degrees worst deviation and never reaches the target (measured).
     *               Passing the RIGIDLY BUILT geometry instead turns the same machinery into a clash
     *               repair -- the fold is already correct, the restraints only hold it while the
     *               optimiser relieves the overlaps.
     */
    bool restrainedBuild(const Proposal& p, Molecule& driven, const Molecule* start = nullptr) const;

    /**
     * @brief Enumerate NCI moves: break a hydrogen bond the template has, form one it does not.
     *
     * The counterpart of generateProposals() on the second descriptor. Only bonds that OCCUR
     * SOMEWHERE in the ensemble are offered for forming -- like the torsion states, this recombines
     * what was observed and does not invent geometry. Candidates are ranked by the population of the
     * target bond (a bond realised in many structures is a plausible one to ask for), which needs no
     * energy model -- consistent with the measurement that the pattern separates but does not predict.
     */
    std::vector<Proposal> generateNCIProposals() const;

    /**
     * @brief Assemble a structure from the individually most favourable elements (de novo template).
     *
     * The mutation stages above are INCREMENTAL: they change one or two torsions of an existing
     * structure, and the measurements show why that is limiting -- the reference structure this work
     * chases shares its torsion vector with an ensemble member, so no small mutation points at it,
     * and it differs from the closest structure in at least seven hydrogen bonds at once.
     *
     * This stage does the opposite. For every torsion it takes the state that is on average the
     * BEST one -- the Boltzmann-weighted mean relative energy of all structures having that state,
     * the same statistic the state table reports -- and assembles all of them into ONE state vector.
     * That vector usually occurs nowhere in the ensemble and is far from every member, which is
     * exactly the point: it is reached in a single concerted build instead of a walk.
     *
     * The per-state energies are used ONLY to choose the geometry. They do not rank, predict or
     * filter anything -- the assembled structure is optimised and judged by the force field like
     * every other proposal. That distinction matters because those same numbers were measured to be
     * poor predictors (scatter 3.7 times their mean): as a template recipe they are still useful,
     * as an energy model they are not.
     *
     * Variants beyond the pure consensus are generated by flipping the torsion with the SMALLEST
     * margin between its best and second-best state -- those are the least certain choices.
     */
    std::vector<Proposal> generateConsensusProposals() const;

    /// Boltzmann-weighted mean relative energy per torsion state (kJ/mol); NaN for empty states.
    std::vector<std::vector<double>> stateEnergies() const;

    /**
     * @brief Build an NCI proposal with DISTANCE restraints instead of dihedral ones.
     *
     * Forming a bond restrains the H...acceptor distance to nci_form_distance, breaking one pushes it
     * to nci_break_distance -- i.e. outside the detection criterion. The rest of the molecule relaxes
     * around that, which is precisely the concerted motion a torsion move cannot express. Restraints
     * are released afterwards; optimiseProposals() reports a freely optimised energy as always.
     */
    bool restrainedBuildNCI(const Proposal& p, Molecule& driven) const;

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
    // Claude Generated (Aug 2026): bound on the candidate list, see generateProposals().
    int m_proposal_candidate_cap = 200000;
    int m_proposal_seed = 0;
    std::string m_analysis_file;
    std::string m_proposal_ranking = "mixed";
    int m_concerted_max = 5;
    mutable std::map<std::string, double> m_reference_terms; ///< terms of the deepest structure (see fitModel)
    mutable std::set<std::vector<int>> m_patterns_before; ///< contact patterns from the memory file
    double m_proposal_novelty_weight = 0.5;   ///< extra structures for the description (see PARAM analysis_file)
    double m_clash_factor = 1.2, m_new_rmsd = 1.0, m_topology_factor = 1.3;
    // Claude Generated (Jul 2026): restrained build (P0). See buildProposalGeometry().
    bool m_restrained_build = true;
    double m_restraint_force = 0.5;
    int m_restraint_max_iterations = 500;

    /**
     * ONE calculator for the whole run: analysis single points, proposal optimisations and the
     * reference re-optimisation. Beyond the physics argument (identical topology, parameters and
     * charge model, so every energy is comparable) this is also a robustness requirement -- creating
     * a second GFN-FF instance for the same molecule in one process reproducibly crashed inside
     * GFNFF::getGFNFFBondParameters (wild pointer in the parameter generation, Jul 2026).
     */
    mutable std::unique_ptr<EnergyCalculator> m_calculator;

    // Claude Generated (Jul 2026): NCI pattern as the second conformer descriptor.
    bool m_nci_analysis = true;
    double m_nci_hbond_distance = 2.60, m_nci_hbond_angle = 120.0;
    double m_nci_xbond_distance = 4.00, m_nci_xbond_angle = 140.0;
    double m_nci_contact_distance = 4.00, m_nci_charge_product = 0.05;
    int m_nci_min_population = 2;
    // NCI moves in the generation stage
    bool m_nci_generate = true;
    int m_nci_max_proposals = 10, m_nci_depth = 1;
    bool m_consensus_build = false;
    std::string m_proposal_memory_file;
    mutable std::set<std::vector<int>> m_proposed_before; ///< loaded from proposal_memory_file
    void loadProposalMemory();
    void appendProposalMemory(
        const std::vector<std::pair<std::vector<int>, std::vector<int>>>& entries) const;
    int m_consensus_max = 3;
    double m_nci_form_distance = 1.90, m_nci_break_distance = 3.50, m_nci_restraint_force = 1.0;

    std::vector<TorsionSpace::Torsion> m_torsions;
    std::vector<std::vector<double>> m_state_centres; ///< per torsion
    std::vector<Frame> m_frames;
    std::vector<NCIContact> m_nci_pairs;   ///< union of all contacts seen in the ensemble
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
    PARAM(concerted_max, Int, 5, "Number of CONCERTED proposals: one torsion change and one hydrogen-bond change in the SAME restrained optimisation, built on top of the two separate move sets. Reason from the measurements: the NCI move alone keeps its target pattern in 7 of 13 cases but lands in an already known minimum in 6 of those 7 -- changing one bond does not by itself change the basin, because the backbone has to follow. The torsion move alone ignores the term that carries +58 % of the energy spread. 0 disables them.", "NCI", {})
    PARAM(proposal_ranking, String, "mixed", "How the enumerated candidates are ordered before the budget cuts them off: energy = by the additive model alone (the behaviour up to Aug 2026), coverage = by the largest Hamming distance to everything already seen or tried, mixed = both, weighted by -proposal_novelty_weight. Measured background: the model predicts energies badly (cross-validated it explains ~30 % of the spread, and a delta model between the two surfaces even degrades the ranking), while the descriptions SEPARATE excellently -- 141 of 142 structures carry a distinct contact pattern. Ordering by energy alone therefore uses the weak property of the description and ignores the strong one.", "Proposals", {})
    PARAM(proposal_novelty_weight, Double, 0.5, "Weight of the novelty term in -proposal_ranking mixed: 0 = energy only, 1 = coverage only. Both contributions are standardised over the candidate set, so the weight is comparable across systems.", "Proposals", {})
    PARAM(analysis_file, String, "", "Additional structures used for the DESCRIPTION only -- torsion states, contact patterns, the additive model and the novelty check. Geometric templates still come from the file the run was called with. Motivated by the measured weakness of the statistics: a cycle of a search delivers a handful of structures (6 in one measured case), while 29 torsions with up to 11 states each and over 100 contacts need far more to be estimated at all. The cumulative pool of the whole run is one to two orders of magnitude larger and costs nothing to reuse.", "Proposals", {})
    PARAM(proposal_templates, Int, 5, "Number of lowest-energy ensemble members used as geometric templates.", "Generation", {})
    PARAM(proposal_depth, Int, 2, "Maximum number of torsions changed simultaneously relative to a template (Hamming distance). Depths beyond 3 are only usable together with -proposal_candidate_cap, which switches the enumeration to sampling.", "Generation", {})
    PARAM(proposal_candidate_cap, Int, 200000, "Upper bound on the number of candidate state vectors held in memory per call. The Hamming ball around a template grows combinatorially with -proposal_depth: measured on a 107-atom peptide with 29 torsions it holds about 4e3 combinations at depth 2, 1e5 at depth 3 and 5e7 at depth 5 -- enumerating the last one exhausted the memory (std::bad_alloc). Below the cap the ball is still enumerated EXACTLY, so nothing changes at the usual depths; above it the same ball is randomly SAMPLED down to the cap. Since the budget keeps only -max_proposals candidates anyway, the sample costs nothing but the guarantee of completeness -- and completeness at depth 5 is unattainable in any case.", "Generation", {})
    PARAM(proposal_seed, Int, 0, "Seed of the candidate sampling (only used when the ball exceeds -proposal_candidate_cap). 0 = derive a seed from the ensemble size and the number of remembered proposals, so a later repetition of the same run draws a DIFFERENT sample instead of the same one again. Any other value fixes the sample and makes the run reproducible.", "Generation", {})
    PARAM(clash_factor, Double, 1.2, "A built structure is rejected when a non-bonded atom pair comes closer than this factor times the sum of their covalent radii. The default is deliberately close to the BOND-DETECTION criterion (~1.3): a built structure that puts two atoms inside bonding distance makes the force field derive a new bond, and the optimisation then relaxes into a different molecule, not a conformer.", "Generation", {})
    PARAM(topology_factor, Double, 1.3, "Covalent-radius factor for the topology check of optimised proposals. A proposal whose bond list differs from the reference is a reaction product, not a conformer, and is rejected. Lower than Molecule's default 1.5, which counts compressed 1-3 contacts as bonds.", "Generation", {})
    PARAM(new_rmsd, Double, 1.0, "Best-fit RMSD in Angstrom above which an optimised proposal counts as a new conformer.", "Generation", {})
    PARAM(restrained_build, Bool, true, "When rigidly setting the torsions produces a clash or a changed bond topology, build the structure by a RESTRAINED optimisation instead: start from the clash-free template, hold the target torsions with a harmonic restraint and let the rest of the molecule relax out of the way, then release. Recovers proposals that the rigid build throws away (measured: 72 percent of them on a compact molecule). False = rigid build only.", "Generation", {})
    PARAM(restraint_force, Double, 0.5, "Force constant of the dihedral restraint in Eh/rad^2 during the restrained build. Larger holds the torsion closer to its target and pushes harder against the clash.", "Generation", {})
    PARAM(restraint_max_iterations, Int, 500, "Maximum optimisation steps of the restrained build stage.", "Generation", {})
    PARAM(cv_folds, Int, 5, "Number of cross-validation folds for the model comparison (additive vs. additive+couplings). Values below 2 disable the comparison.", "Analysis", {})
    PARAM(couplings, Bool, true, "Measure torsion-torsion couplings from double-mutant cycles (four ensemble members forming a rectangle in state space) and run the model comparison.", "Analysis", {})
    PARAM(report_threshold, Double, 1.0, "Only transitions whose mean total energy difference exceeds this many kJ/mol are printed in the summary. All of them are written to the CSV.", "Analysis", {})

    PARAM(nci_analysis, Bool, true, "Describe every conformer additionally by its pattern of non-covalent interactions (hydrogen bonds, halogen bonds, electrostatic contacts, close contacts) and compare that description with the torsion states. Needed because the per-term attribution shows the energy of a conformer is carried by the non-covalent terms, not by the torsion term.", "NCI", {})
    PARAM(nci_hbond_distance, Double, 2.60, "Maximum H...acceptor distance in Angstrom for a hydrogen bond of the NCI pattern.", "NCI", {})
    PARAM(nci_hbond_angle, Double, 120.0, "Minimum donor-H...acceptor angle in degrees for a hydrogen bond of the NCI pattern.", "NCI", {})
    PARAM(nci_xbond_distance, Double, 4.00, "Maximum halogen...acceptor distance in Angstrom for a halogen bond of the NCI pattern.", "NCI", {})
    PARAM(nci_xbond_angle, Double, 140.0, "Minimum C-halogen...acceptor angle in degrees for a halogen bond (sigma hole is directional).", "NCI", {})
    PARAM(nci_contact_distance, Double, 4.00, "Maximum distance in Angstrom between two heavy atoms that are at least four bonds apart for them to count as a non-covalent contact.", "NCI", {})
    PARAM(nci_charge_product, Double, 0.05, "A contact whose partial-charge product is more negative than minus this value is classified as an electrostatic (ionic) contact rather than a plain close contact.", "NCI", {})
    PARAM(nci_min_population, Int, 2, "A contact must occur in at least this many structures to enter the NCI pattern. Contacts seen once carry no contrast and only inflate the model.", "NCI", {})

    PARAM(nci_generate, Bool, true, "With -generate true, additionally propose NCI MOVES: break a hydrogen bond the template has and form one that occurs elsewhere in the ensemble, realised by distance restraints. This is the move set that acts on the description which actually distinguishes the conformers (the H-bond term carries the energy spread, the torsion term does not), and it reaches structures a torsion recombination cannot express.", "NCI", {})
    PARAM(nci_max_proposals, Int, 10, "Maximum number of NCI moves built and optimised, in addition to the torsion proposals.", "NCI", {})
    PARAM(nci_depth, Int, 1, "Number of hydrogen bonds changed simultaneously in one NCI move.", "NCI", {})
    PARAM(nci_form_distance, Double, 1.90, "Target H...acceptor distance in Angstrom when an NCI move FORMS a hydrogen bond.", "NCI", {})
    PARAM(nci_break_distance, Double, 3.50, "Target H...acceptor distance in Angstrom when an NCI move BREAKS a hydrogen bond. Must be clearly outside nci_hbond_distance, otherwise the bond re-forms during the free optimisation.", "NCI", {})
    PARAM(nci_restraint_force, Double, 1.0, "Force constant of the distance restraint in Eh/Angstrom^2 during an NCI move.", "NCI", {})

    PARAM(consensus_build, Bool, false, "With -generate true, additionally assemble structures DE NOVO from the individually most favourable torsion states instead of mutating an existing one, walking away from the ensemble one torsion at a time. Reaches state vectors that no sequence of one- or two-torsion mutations can reach: measured, every chemically valid assembly was a new conformer, but they sit high in energy (best +44 kJ/mol) and each costs a build plus two optimisations. Off by default for that cost.", "Generation", {})
    PARAM(proposal_memory_file, String, "", "Path to a file recording the state vectors that have already been proposed. When set, ConfGen skips combinations listed there and appends the ones it builds. ConfSearch passes one file per run, so a temperature stage does not rebuild what an earlier repetition already tried -- measured: repetitions 2 and 3 of a 600 K stage proposed the same structures again, and 32 of 113 proposals across a 7-cycle run were repeats.", "Generation", {})
    PARAM(consensus_max, Int, 3, "Number of de-novo assemblies built: the pure consensus plus variants that flip the torsions whose best and second-best state are closest in energy (the least certain choices).", "Generation", {})

    END_PARAMETER_DEFINITION
    // ^^^^^^^^^^^^ PARAMETER DEFINITION BLOCK ^^^^^^^^^^^^
};
