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
    /** Claude Generated (Sep 2026): number of hydrogen bonds (donor/acceptor pairs) present in exactly
     *  one of the two structures -- the symmetric difference of their H-bond sets. Geometric criterion
     *  only (detectNCI without charges). Used by the novelty rule of PARAM new_energy_gain. */
    int hbondHamming(const Molecule& a, const Molecule& b) const;

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
        /* Claude Generated (Sep 2026): the ensemble member the optimised proposal is closest to
         * (best-fit RMSD), and whether the proposal entered the ensemble through the energy/pattern
         * rule of PARAM new_energy_gain rather than the RMSD threshold. */
        int nearest_frame = -1;
        bool deeper = false;
        bool unconverged = false;  ///< kept although the optimiser hit its step cap (Sep 2026), like RELAX does
        /* Claude Generated (Sep 2026): ROUTE move -- the torsion states were chosen so that a rigid
         * rotation brings the intended donor/acceptor pair close (see generateRouteProposals()).
         * route_start is that rigidly rotated, clash-free geometry; the restrained build starts
         * from it instead of from the template. */
        Molecule route_start;
        double route_rigid_distance = 0.0; ///< H...A after the rigid rotation, Angstrom
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
        /**
         * Claude Generated (Aug 2026): explicit target ANGLES in degrees, as (torsion index,
         * angle). Needed because the isomerisation move drives a torsion to a value that is NOT
         * an observed rotamer state -- and a state index cannot express what the ensemble has
         * never shown. Measured on a 107-atom peptide: the two deepest known structures differ
         * from every one of 99 ensemble members in the guanidinium C-N torsion, which the whole
         * ensemble holds in a single state at +174 deg while they need -22 and +8 deg. A
         * recombination generator can only reassemble states it has seen, so no proposal depth
         * whatever can build them -- the BUILDING BLOCK is missing, not the reach.
         */
        std::vector<std::pair<int, double>> angle_targets;
        bool isomerisation() const { return !angle_targets.empty(); }
        std::string nci_label;     ///< human-readable move, e.g. "break HB 38-90...12, form HB 5-61...46"
        /** Claude Generated (Aug 2026): the geometry is already there and must not be rebuilt --
         *  a collective-mode displacement has no torsion target to drive towards. The clash and
         *  topology gates still apply. */
        bool prebuilt = false;
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
     * @brief Add the caller's polar X-H restraints to an optimiser config, next to its own ones.
     *
     * Claude Generated (Sep 2026): every optimisation this class runs can transfer a proton, and
     * until now none of them was guarded -- ConfSearch's -hold_polar_h reached RELAX and the
     * re-scoring pass, but not the proposal optimisations, because PerformConfGen never handed the
     * restraints down. Measured on a 107-atom peptide: 275 of 1049 optimised proposals of one
     * production run (26 %, 183 of them NCI moves) came out with a changed bond topology, H107
     * having moved from O57 to N25. They were correctly rejected afterwards -- .proposals.new.xyz
     * requires topology_ok -- so nothing entered the pool, but a quarter of the proposal budget
     * bought tautomers. A second run under the same code lost only 3 %, so the rate depends on how
     * strained the templates are, not on the move set alone.
     *
     * MERGE, never overwrite: an NCI move IS a set of distance restraints, and dropping them would
     * turn the move into a plain free optimisation. Own restraints win on a duplicated atom pair,
     * so a move that deliberately pulls on an X-H is not fought by this guard.
     */
    nlohmann::json withPolarHydrogenRestraints(const nlohmann::json& own) const;

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
     * @param calc   optional calculator to use instead of the shared m_calculator (Sep 2026: lets the
     *               parallel build loop give each worker its own instance -- see its call site).
     */
    bool restrainedBuild(const Proposal& p, Molecule& driven, const Molecule* start = nullptr,
        EnergyCalculator* calc = nullptr) const;

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

    /** Claude Generated (Aug 2026): crossover -- transfer a connected window of torsions from one
     *  conformer into another. The move for the measured 6-7-torsion gap that mutation cannot
     *  bridge; see the comment at the definition. */
    std::vector<Proposal> generateCrossoverProposals() const;

    /** Claude Generated (Aug 2026): collective modes -- displace a template along the principal
     *  components of the ensemble's own coordinate covariance. The complement to crossover: not
     *  restricted to observed combinations; see the comment at the definition. */
    std::vector<Proposal> generateModeProposals() const;

    /** Claude Generated (Aug 2026): path images -- the structures BETWEEN two known conformers that
     *  are far apart, rather than around one. See the comment at the definition. */
    std::vector<Proposal> generatePathProposals() const;

    /** Claude Generated (Aug 2026): the concerted move -- one torsion and one hydrogen bond in the
     *  same restrained optimisation, the torsion chosen by geometric coupling to the bridge. Until
     *  now -concerted_max was read and never used; see the definition. */
    std::vector<Proposal> generateConcertedProposals() const;

    /**
     * @brief Flip a torsion the ensemble only ever showed in ONE planar state (isomerise_max).
     *
     * Claude Generated (Aug 2026). Every other move set recombines OBSERVED rotamer states, so a
     * conformer that needs an unobserved state is unreachable at any depth. This is the only move
     * that can ADD a state to the space.
     *
     * OFF BY DEFAULT. It was built on a measurement that turned out to be an artefact of how the
     * measurement was made: the deepest known structure and an external reference both differ from
     * all 99 members of ONE CYCLE ensemble in the guanidinium C-N torsion. But ConfSearch hands
     * ConfGen the CUMULATIVE pool as its description basis (see analysis_file), and there that
     * torsion has ten states -- including -19 deg (n=109) and +31 deg (n=58), which is where those
     * two structures sit. Nothing was missing, and the earlier claim that this explains why a
     * larger proposal_depth never helped is withdrawn. Measured afterwards: 2 of 29 torsions are
     * single-state in the cumulative ensemble, none at all in a 223-structure force-field one.
     *
     * The mechanism itself is verified and cheap, so it stays available for the case it was meant
     * for -- a degree of freedom the WHOLE run never opened.
     *
     * Nothing here is keyed to a molecule. DETECTION is read off the rotamer analysis: a torsion
     * the ensemble showed in a SINGLE state is a degree of freedom the sampling never opened.
     * The TARGET is measured, not assumed -- the frozen torsion is scanned rigidly and every local
     * minimum of that profile at least isomerise_min_separation away from the known state becomes
     * a candidate. An earlier version proposed "the opposite planar value", which is chemical
     * knowledge brought in after seeing one system and blind to a non-planar second minimum.
     * Everything else -- restrained build, release, free optimisation, clash/topology/novelty
     * gates -- is the common path.
     */
    std::vector<Proposal> generateIsomerisationProposals() const;

    /** Claude Generated (Aug 2026): the calculator that judges a proposal (see -eval_method); the
     *  shared description calculator when the two methods are the same. */
    EnergyCalculator* evaluationCalculator() const;

    /** Claude Generated (Aug 2026): is a SEPARATE evaluation surface actually in use? The plain
     *  test `!m_eval_method.empty() && m_eval_method != m_method` was written at three call sites
     *  and is not enough: it stays true after the evaluation calculator has been rejected as
     *  unusable, so the method name still reached opt_config and the reference-energy calculator
     *  while the judge had already fallen back -- which segfaulted on an unknown method name.
     *  Consulting this instead keeps all sites on one answer, and it triggers the one-off probe. */
    bool useSplitEvaluation() const;

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
     * @param calc   optional calculator to use instead of the shared m_calculator (Sep 2026: lets the
     *               parallel build loop give each worker its own instance -- see its call site).
     */
    bool restrainedBuildNCI(const Proposal& p, Molecule& driven, const Molecule* start = nullptr,
        EnergyCalculator* calc = nullptr) const;
    /** Claude Generated (Sep 2026): the torsion route to a hydrogen bond -- for an absent D-H...A pair,
 *  the combination of observed rotamer states of the torsions on the bond path D->A that brings H and A
 *  closest by RIGID rotation, then the restrained build from that geometry. Unifies the two move sets:
 *  the goal is a contact, the means is a rotation around bonds (the physical path), not a pull through
 *  the molecule (measured: 0 of 45 contact pulls survive the clash gate). */
    std::vector<Proposal> generateRouteProposals() const;
    /** Claude Generated (Sep 2026): template order shared by every move set -- the lowest-energy input
 *  frame first, then (with -template_diversity) frames picked greedily for the largest H-bond-pattern
 *  distance to the ones already chosen, so the proposals do not all start in one basin. */
    std::vector<int> templateOrder() const;
    /** Claude Generated (Sep 2026): steered relaxation -- a short MD with the move's distance restraints
     *  (and the polar-H guard) still acting, so the rest of the molecule adapts to the new bond before
     *  the restraint is released. Returns the final MD geometry (the caller re-optimises restrained). */
    Molecule steerRelax(const Molecule& start, const nlohmann::json& distance_restraints, int seed) const;
    /** Claude Generated (Sep 2026): would forming pair k on top of pattern `have` over-saturate a partner?
 *  An acceptor already holding two bonds, or a donor hydrogen already donating, is skipped. */
    bool saturated(int k, const std::vector<int>& have) const;


    /** Gemessenes Ziel eines Kontaktzugs: Median des Abstands ueber die Traeger im Ensemble. */
    double contactTargetDistance(int pair_index) const;

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
    // Claude Generated (Aug 2026): isomerisation move, see generateIsomerisationProposals().
    int m_isomerise_max = 0, m_isomerise_scan_steps = 24;
    double m_isomerise_min_separation = 60.0, m_isomerise_max_rise = 100.0;
    // Claude Generated (Aug 2026): crossover move, see generateProposals().
    int m_crossover_max = 0, m_crossover_window = 6;
    // Claude Generated (Aug 2026): collective-mode move, see generateModeProposals().
    int m_mode_max = 0, m_mode_count = 3;
    double m_mode_amplitude = 2.0;
    // Claude Generated (Aug 2026): path images, see generatePathProposals().
    int m_path_max = 0, m_path_images = 3, m_path_min_distance = 6;
    // Claude Generated (Aug 2026): staged realisation of a multi-bridge NCI move.
    bool m_nci_staged = false;
    // Claude Generated (Aug 2026): the surface that JUDGES a proposal, see evaluationCalculator().
    std::string m_eval_method;
    mutable std::unique_ptr<EnergyCalculator> m_eval_calculator;
    /// Set once the evaluation method has been proven unusable; every split test then reads false.
    mutable bool m_eval_unusable = false;
    std::string m_analysis_file;
    std::string m_proposal_ranking = "mixed";
    std::string m_nci_ranking = "population";
    bool m_nci_contact_moves = false;         ///< Claude Generated (Sep 2026): Kompaktierungszug
    double m_nci_contact_break_distance = 8.0;
    mutable std::vector<bool> m_atom_is_hydrogen; ///< je Atom: ist es ein H (fuer die Stufenreihenfolge)
 ///< Claude Generated (Sep 2026): Richtung des Populationsterms im NCI-Zug
    int m_concerted_max = 5;
    mutable std::map<std::string, double> m_reference_terms; ///< terms of the deepest structure (see fitModel)
    mutable std::set<std::vector<int>> m_patterns_before; ///< contact patterns from the memory file
    double m_proposal_novelty_weight = 0.5;   ///< extra structures for the description (see PARAM analysis_file)
    double m_clash_factor = 1.2, m_new_rmsd = 1.0, m_topology_factor = 1.3;
    double m_new_energy_gain = 1.0; ///< see PARAM new_energy_gain
    int m_new_pattern_min = 1;      ///< see PARAM new_pattern_min
    // Claude Generated (Jul 2026): restrained build (P0). See buildProposalGeometry().
    bool m_restrained_build = true;
    double m_restraint_force = 0.5;
    int m_restraint_max_iterations = 500;
    /**
     * Claude Generated (Sep 2026): polar X-H distance restraints handed down by the caller, in the
     * json form GeometryRestraints::fromJson reads. Empty = none.
     *
     * Why the caller supplies them instead of this class deriving them: the reference bond lengths
     * belong to the SURFACE the proposals are optimised on, and only the caller knows which
     * structure defines it (ConfSearch keeps one reference per PES). See PolarHydrogenRestraints()
     * there and withPolarHydrogenRestraints() below.
     */
    nlohmann::json m_polar_h_restraints = nlohmann::json::array();

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
    // Claude Generated (Sep 2026): route move, saturation rule, clamp pairs, template diversity
    int m_route_max = 6, m_route_depth = 3;
    double m_route_reach = 3.0;
    bool m_nci_saturation = true, m_nci_clamp = true, m_template_diversity = true;
    double m_nci_steer_fs = 0.0, m_nci_steer_temperature = 300.0; ///< see PARAM nci_steer_fs

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
    PARAM(isomerise_max, Int, 0, "Number of ISOMERISATION proposals: flip a torsion the ensemble only ever showed in ONE state to a second minimum located by a scan. Every other move set recombines OBSERVED states, so a conformer needing an unobserved one is unreachable at any depth -- this is the only move that can ADD a state. OFF BY DEFAULT, because the case it repairs was not found in production. The move was built after measuring that the deepest known structure and an external reference both differ from all 99 members of one CYCLE ensemble in the guanidinium C-N torsion, which those 99 hold at +174 deg. That measurement was made on the wrong ensemble: ConfSearch hands ConfGen the CUMULATIVE pool as its description basis (analysis_file), and there the same torsion has ten states including -19 deg (n=109) and +31 deg (n=58) -- exactly where the two structures sit. Nothing was missing. In the cumulative ensemble only 2 of 29 torsions are single-state, and in a 223-structure force-field ensemble none at all. So the mechanism works (verified: the scan locates the second minimum, the build reaches it, the structure keeps it after free optimisation) but it addresses an impoverishment that only arises when ConfGen is run on a single cycle without analysis_file. Turn it on for that case, or when a rotation barrier is suspected to be uncrossed by the dynamics. 0 disables it.", "Proposals", {})
    PARAM(isomerise_scan_steps, Int, 24, "Points of the rigid torsion scan that LOCATES the second minimum of a single-state torsion (24 = every 15 degrees). The scan replaces an assumption with a measurement: an earlier version proposed 'the opposite planar value', which is chemical knowledge brought in after seeing one system, and it would miss a second minimum that is not planar. One energy per point on the description surface, per frozen torsion -- a handful of single points, not optimisations.", "Proposals", {})
    PARAM(isomerise_max_rise, Double, 100.0, "How far above the known state a scan minimum may sit, in kJ/mol, before it is dropped. The scan is RIGID -- it rotates the moving side and relaxes nothing -- so it overestimates badly and the threshold has to be generous. It still has to exist: measured on a 107-atom peptide, the same scan that correctly located the missing guanidinium state at +4.5 kJ/mol also reported minima at +2251 and +40475 kJ/mol, which are atoms driven into each other, not conformers. Those cost a restrained build each and cannot produce anything.", "Proposals", {})
    PARAM(isomerise_min_separation, Double, 60.0, "How far a scan minimum must lie from the state the ensemble already knows, in degrees, before it counts as a different state. Below this it is the same basin seen from its flank.", "Proposals", {})
    PARAM(concerted_max, Int, 5, "Number of CONCERTED proposals: one torsion change and one hydrogen-bond change in the SAME restrained optimisation. Implemented Aug 2026 -- before that the value was read and never used, and this help text described a move set that did not exist. The torsion is chosen by geometric coupling: its rotation must move exactly ONE of the two bridge partners, since only then does it change their relative position. Reason from the measurements: the NCI move alone keeps its target pattern in 7 of 13 cases but lands in an already known minimum in 6 of those 7 -- changing one bond does not by itself change the basin, the backbone has to follow. Measured with the move in place, same bridge budget: 4 new conformers instead of 2, bridge distance to the reference 11 instead of 13. 0 disables them.", "NCI", {})
    PARAM(proposal_ranking, String, "mixed", "How the enumerated candidates are ordered before the budget cuts them off: energy = by the additive model alone (the behaviour up to Aug 2026), coverage = by the largest Hamming distance to everything already seen or tried, mixed = both, weighted by -proposal_novelty_weight. Measured background: the model predicts energies badly (cross-validated it explains ~30 % of the spread, and a delta model between the two surfaces even degrades the ranking), while the descriptions SEPARATE excellently -- 141 of 142 structures carry a distinct contact pattern. Ordering by energy alone therefore uses the weak property of the description and ignores the strong one.", "Proposals", {})
    PARAM(proposal_novelty_weight, Double, 1.0, "Weight of the novelty (coverage) term in -proposal_ranking mixed; 1.0 = order the candidates by coverage alone, 0.0 = by the additive torsion-state model alone. Default 1.0 since Sep 2026: the model was measured at r = -0.02 between its ordering and the optimised energy of the same proposals (15 proposals, triose) and explains ~15 % of the energy variance in cross-validation on two peptide ensembles -- it has nothing to say about which candidate is worth an optimisation, while coverage (distance in state and contact space to everything tried) is the property the description is good at.", "Generation", {})
    PARAM(analysis_file, String, "", "Additional structures used for the DESCRIPTION only -- torsion states, contact patterns, the additive model and the novelty check. Geometric templates still come from the file the run was called with. Motivated by the measured weakness of the statistics: a cycle of a search delivers a handful of structures (6 in one measured case), while 29 torsions with up to 11 states each and over 100 contacts need far more to be estimated at all. The cumulative pool of the whole run is one to two orders of magnitude larger and costs nothing to reuse.", "Proposals", {})
    PARAM(proposal_templates, Int, 5, "Number of lowest-energy ensemble members used as geometric templates.", "Generation", {})
    PARAM(proposal_depth, Int, 2, "Maximum number of torsions changed simultaneously relative to a template (Hamming distance). Depths beyond 3 are only usable together with -proposal_candidate_cap, which switches the enumeration to sampling.", "Generation", {})
    PARAM(proposal_candidate_cap, Int, 200000, "Upper bound on the number of candidate state vectors held in memory per call. The Hamming ball around a template grows combinatorially with -proposal_depth: measured on a 107-atom peptide with 29 torsions it holds about 4e3 combinations at depth 2, 1e5 at depth 3 and 5e7 at depth 5 -- enumerating the last one exhausted the memory (std::bad_alloc). Below the cap the ball is still enumerated EXACTLY, so nothing changes at the usual depths; above it the same ball is randomly SAMPLED down to the cap. Since the budget keeps only -max_proposals candidates anyway, the sample costs nothing but the guarantee of completeness -- and completeness at depth 5 is unattainable in any case.", "Generation", {})
    PARAM(crossover_max, Int, 0, "Number of CROSSOVER proposals: transfer a connected window of torsions from one conformer of the ensemble into another, instead of changing torsions independently. 0 = off. Why it exists: the reference conformer of the measured 107-atom peptide differs from the nearest of ~4600 sampled structures in 6-7 of 29 torsions, further than any mutation reaches in one step -- and raising -proposal_depth that far is useless because NOT ONE deep-mutation structure kept its target state vector through the free optimisation (a random combination of seven torsions is almost never near a minimum). A transferred window carries a fold that already exists somewhere, so only its context changes. This is fragment insertion with the ensemble as the fragment library.", "Generation", {})
    PARAM(crossover_window, Int, 6, "Number of torsions transferred per crossover, grown breadth-first from a seed torsion over the CONNECTIVITY of the central bonds (a window in storage order would be an arbitrary set of unrelated dihedrals). The default matches the measured gap of 6-7 torsions to the reference conformer.", "Generation", {})
    PARAM(mode_max, Int, 0, "Number of COLLECTIVE-MODE proposals: displace a template along the principal components of the ensemble's own coordinate covariance (essential dynamics), then optimise freely. 0 = off. The complement to -crossover_max: a crossover can only transfer folds the ensemble already contains, a mode displacement is not restricted to observed combinations. Motivation from the same measurement -- the reference conformer sits 6-7 torsions away, and six or seven torsions do not move independently; the leading eigenvectors of the ensemble ARE the collective directions that rotate five to ten of them together.", "Generation", {})
    PARAM(mode_count, Int, 3, "How many of the leading collective modes are used. Each contributes two proposals per template (plus and minus displacement).", "Generation", {})
    PARAM(mode_amplitude, Double, 2.0, "Displacement along a mode in units of its own standard deviation sqrt(eigenvalue). 1 stays inside the sampled range, larger values extrapolate beyond it -- which is the point, but also where the clash and topology gates start rejecting.", "Generation", {})
    PARAM(eval_method, String, "", "Method used to OPTIMISE and JUDGE the proposals, when it should differ from the one that DESCRIBES the ensemble (-method). Empty = same for both. ConfGen used one method for everything, which mixes two surfaces that were measured not to agree: within one cycle of a 107-atom peptide GFN-FF and GFN2 rank the same structures at r = -0.32 ... -0.46, so every proposal was selected on a surface pointing elsewhere. Measured consequence: the collective-mode move produced the deepest GFN-FF structure this generator ever made -- 20.4 kJ/mol below the ensemble minimum -- and the same structures sit 87 to 134 kJ/mol ABOVE the reference on GFN2. With -eval_method the description keeps the force field it needs for the term decomposition while the proposals are optimised and compared on the accurate surface. The comparison base (the re-optimised templates) moves with it, so the reported gain stays on one scale.", "Generation", {})
    PARAM(nci_staged, Bool, false, "Realise an NCI move with more than one bridge in STAGES instead of pulling every distance restraint at once: close the bridge whose current distance is nearest its target, let the molecule relax around it, keep it restrained and add the next. Only relevant with -nci_depth > 1. Motivated by the measured failure of the simultaneous variant: taking the ensemble structure closest to a reference conformer (2.61 A) and closing all six of its target hydrogen bonds at once left the RMSD at 2.55-2.64 A and the optimisation stalled -- six distances are six conditions in 315 degrees of freedom, and the optimiser meets them with a local compromise instead of refolding. Costs one restrained optimisation per bridge.", "NCI", {})
    PARAM(path_max, Int, 0, "Number of PATH IMAGE proposals: structures that lie BETWEEN two known conformers instead of around one. 0 = off. Every other move set steps away from a single structure; this one exploits that two ensemble members which are themselves 6-7 torsions apart have that region between them -- and a point halfway along is only 3-4 torsions from either end, a distance the build and optimisation handle well. The path is the straight line in torsion space (circular interpolation per dihedral, snapped to populated states), not a minimum-energy path: it costs nothing and tests whether the missing structures lie between the known ones before the NEB machinery is brought in.", "Generation", {})
    PARAM(path_images, Int, 3, "Intermediate points generated per pair of endpoints.", "Generation", {})
    PARAM(path_min_distance, Int, 6, "Minimum torsion Hamming distance between two structures before the region between them is sampled. The default matches the measured gap to the reference conformer.", "Generation", {})
    PARAM(proposal_seed, Int, 0, "Seed of the candidate sampling (only used when the ball exceeds -proposal_candidate_cap). 0 = derive a seed from the ensemble size and the number of remembered proposals, so a later repetition of the same run draws a DIFFERENT sample instead of the same one again. Any other value fixes the sample and makes the run reproducible.", "Generation", {})
    PARAM(clash_factor, Double, 1.2, "A built structure is rejected when a non-bonded atom pair comes closer than this factor times the sum of their covalent radii. The default is deliberately close to the BOND-DETECTION criterion (~1.3): a built structure that puts two atoms inside bonding distance makes the force field derive a new bond, and the optimisation then relaxes into a different molecule, not a conformer.", "Generation", {})
    PARAM(topology_factor, Double, 1.3, "Covalent-radius factor for the topology check of optimised proposals. A proposal whose bond list differs from the reference is a reaction product, not a conformer, and is rejected. Lower than Molecule's default 1.5, which counts compressed 1-3 contacts as bonds.", "Generation", {})
    PARAM(new_rmsd, Double, 1.0, "Best-fit RMSD in Angstrom above which an optimised proposal counts as a new conformer.", "Generation", {})
    PARAM(new_energy_gain, Double, 1.0, "Second way for an optimised proposal to count as NEW, next to the RMSD threshold of -new_rmsd: it is a valid conformer (topology intact), it lies BELOW its closest ensemble member by more than this many kJ/mol on the SAME energy surface, and its hydrogen-bond pattern differs from that member by at least -new_pattern_min bonds. 0 = off (RMSD rule only, the behaviour up to Sep 2026). Why: a conformer is a minimum, and two topologically identical minima with different energies are different conformers no matter how close they are in Angstrom -- the RMSD rule alone throws the deeper one away when it sits inside the threshold. Measured on a 107-atom peptide over two production runs (35 + 20 repetitions): 61 and 38 valid proposals landed 0.3-1.25 A from a member AND more than 1 kJ/mol below it, 84-85 percent of them with a different H-bond pattern (control group under 0.3 A: 7-9 percent). One of them, built one repetition before the dynamics found it, was 5.2 kJ/mol below everything in the pool and 0.03 A heavy-atom RMSD from what later became that run's record -- discarded as not new. Requires that the input energies (XYZ comment lines) and the proposal optimisation are on the same surface with the same restraints; ConfSearch guarantees that and switches the rule off when its evaluation method differs from the RELAX method. A standalone -confgen run on a file with energies from another method must set 0. The deduplication at the end keeps the lower member of a basin, so an accepted deeper neighbour replaces rather than duplicates. Accepted structures carry the suffix _deeper in their name.", "Generation", {})
    PARAM(new_pattern_min, Int, 1, "Minimum number of hydrogen bonds (donor/acceptor pairs) in which a proposal accepted through -new_energy_gain must differ from its closest ensemble member. 0 = energy alone decides. 1 (default) demands a visibly different contact pattern, which is what separated the deeper neighbours from re-optimised copies in the measurement behind -new_energy_gain.", "Generation", {})
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
    PARAM(nci_contact_moves, Bool, false, "Admit CLOSE and IONIC contacts to the NCI move set, not only hydrogen bonds, restrained heavy-atom to heavy-atom with a MEASURED target (the median distance over exactly those ensemble members that carry the contact). Such a move runs staged and first, and contact moves are ordered by the FOLD they demand (current distance minus target) instead of by population. MEASURED RESULT: IT DOES NOT WORK, and the number says why. On WEKLQ, ranked by fold, 0 of 45 contact moves survive -- 44 are rejected by the clash/topology gate, 1 fails to reach its restraint. Merely admitting the contacts without the fold ranking builds 8 of 45 but comes out LESS compact than without them (R_gyr 4.275 against 4.202), because the population order picks some short contact among 467 movable ones instead of the one that folds the molecule. The reason is structural, not a tuning problem: a harmonic restraint between two distant heavy atoms satisfies its distance along the shortest path, which runs THROUGH the rest of the molecule. Compaction cannot be imposed by a static restrained optimisation; it has to come from dynamics under a global compaction bias. Kept as a documented negative and as the machinery for that experiment -- do not expect proposals from it.", "NCI", {})
    PARAM(nci_contact_break_distance, Double, 8.0, "Target distance in Angstrom when a contact move BREAKS a close contact. Must be clearly outside -nci_contact_distance (4.0), otherwise the contact re-forms during the free optimisation -- same reasoning as -nci_break_distance for hydrogen bonds, one length scale up because heavy-atom contacts are detected at a larger distance.", "NCI", {})
    PARAM(nci_ranking, String, "population", "Direction of the population term when NCI moves are ordered. population = prefer FORMING a hydrogen bond that many ensemble members realise (the behaviour up to Sep 2026). rarity = the opposite sign, i.e. prefer forming the bonds that are RARE in the ensemble. Why the switch exists: measured on a 166-structure WEKLQ ensemble constructed so that all nine hydrogen bonds of a known deep basin are present (7 to 80 carriers each), the templates carry up to 5 of those bonds and the generated proposals carry at most 2 at -nci_depth 1 and at most 3 at depth 3 with a four-fold budget. The move set therefore walks AWAY from that basin rather than merely approaching it slowly, and neither depth nor budget nor -nci_min_population 1 changes that. The mechanism is the sign in the term: breaking a RARE bond costs almost nothing while forming a COMMON one pays, so every structure drifts out of the rare-pattern corner into the population bulk. The counterweight already present (-proposal_ranking coverage / -proposal_novelty_weight) measures distance in contact space and does not distinguish rare from common. Not the default, because inverting it is untested on any system where the population term is the right prior.", "NCI", {})
    PARAM(nci_max_proposals, Int, 10, "Maximum number of NCI moves built and optimised, in addition to the torsion proposals.", "NCI", {})
    PARAM(nci_depth, Int, 1, "Number of hydrogen bonds changed simultaneously in one NCI move.", "NCI", {})
    PARAM(nci_form_distance, Double, 1.90, "Target H...acceptor distance in Angstrom when an NCI move FORMS a hydrogen bond.", "NCI", {})
    PARAM(nci_break_distance, Double, 3.50, "Target H...acceptor distance in Angstrom when an NCI move BREAKS a hydrogen bond. Must be clearly outside nci_hbond_distance, otherwise the bond re-forms during the free optimisation.", "NCI", {})
    PARAM(nci_restraint_force, Double, 1.0, "Force constant of the distance restraint in Eh/Angstrom^2 during an NCI move.", "NCI", {})
    PARAM(nci_steer_fs, Double, 0.0, "Steered relaxation of a built NCI or route proposal: after the restrained build, run this many femtoseconds of MD on the description method with the move's distance restraints (and the polar-H guard) still acting, then re-optimise restrained and hand the structure to the free optimisation. 0 = off. Why: measured on a 107-atom peptide, 25-55 percent of the built hydrogen-bond moves relax back into their template once the restraint is released -- the bond was closed, but the rest of the molecule never moved to accommodate it. A static restrained optimisation finds the nearest local compromise; a few hundred femtoseconds of thermal motion under the same restraint let the environment reorganise along a physical path. Costs one short force-field MD per proposal.", "Generation", {})
    PARAM(nci_steer_temperature, Double, 300.0, "Temperature in Kelvin of the steered relaxation (-nci_steer_fs).", "Generation", {})
    PARAM(route_max, Int, 6, "Number of ROUTE proposals per call (0 = off). A route move forms one hydrogen bond the template lacks by ROTATING: among the observed rotamer states of the torsions on the bond path between donor and acceptor (up to -route_depth of them at once) it takes the combination that brings H and acceptor closest to the form distance by rigid rotation, checks that geometry for clashes, and only then pulls the bond closed in the restrained build, starting from the rotated geometry. Why: a hydrogen bond between distant groups forms because a side chain or the backbone turns around bonds -- the physical path -- and not because two atoms are dragged through the molecule. Measured on a 125-structure peptide ensemble (gfnff build): 20 of 20 routes built with 0 clash rejections and 0 unreached restraints, against 7 of 15 for the pull-based NCI moves; 19-20 of 20 new by RMSD; 15-17 of 20 keep the bond after the free optimisation. Judged on gfn2 (12 routes, polar-H restraints): all 12 new, but 36 kJ/mol above the ensemble minimum, 6 of 12 snapped back, at most 2 of the 9 target bonds -- a reliable BUILDER, not yet a proven source of depth; n = 1 ensemble. Needs no target knowledge: pairs come from the ensemble's own contact vocabulary, states from its own rotamers.", "Generation", {})
    PARAM(route_depth, Int, 3, "Maximum number of path torsions changed together in one route move (-route_max). Measured on the peptide: depth 3 is where independent recombination stopped paying; the route changes only torsions between the two partners, so its combinations stay small.", "Generation", {})
    PARAM(route_reach, Double, 3.0, "A route counts only when the rigid rotation brings H...acceptor below this distance in Angstrom (and improves it by at least 1 A). From there the restraint closes the last stretch to -nci_form_distance without pulling through the molecule.", "Generation", {})
    PARAM(nci_saturation, Bool, true, "Do not propose to FORM a hydrogen bond whose acceptor already holds two bonds or whose donor hydrogen already donates one, in the template's own pattern. Plain chemistry -- a carbonyl oxygen accepts twice, an N-H donates once -- and it removes moves the free optimisation undoes anyway.", "Generation", {})
    PARAM(nci_clamp, Bool, false, "Admit CLAMP moves: two absent hydrogen bonds that share a heavy atom (the same side-chain nitrogen accepting one bond and donating another, or one acceptor taking a second donor) are formed together in one staged restrained build. Why it exists: the deepest basin of the measured peptide is defined by exactly such a bidirectional clamp of both basic side chains, and single-bond moves reached at most 2-3 of its 9 bonds. Why it is OFF: measured on a 125-structure ensemble carrying all nine bonds, 14 of 21 clamp builds died at the clash gate -- pulling two bonds through the molecule at once is the failure mode the route move avoids -- and the survivors displaced single moves inside nci_max_proposals without reaching more of the motif (max 2 of 9). The clamp that works has to be built by rotation (a two-target route), which is not implemented.", "Generation", {})
    PARAM(template_diversity, Bool, true, "Choose the proposal templates for pattern diversity, not only energy: the lowest-energy input structure first, then greedily the candidate (from the lowest 3 x proposal_templates by energy) whose hydrogen-bond pattern differs most from the templates already chosen. Otherwise every move set starts from the same basin, and the measured record chain shows that key steps came from seeds ranked 5-10.", "Generation", {})

    PARAM(consensus_build, Bool, false, "With -generate true, additionally assemble structures DE NOVO from the individually most favourable torsion states instead of mutating an existing one, walking away from the ensemble one torsion at a time. Reaches state vectors that no sequence of one- or two-torsion mutations can reach: measured, every chemically valid assembly was a new conformer, but they sit high in energy (best +44 kJ/mol) and each costs a build plus two optimisations. Off by default for that cost.", "Generation", {})
    PARAM(proposal_memory_file, String, "", "Path to a file recording the state vectors that have already been proposed. When set, ConfGen skips combinations listed there and appends the ones it builds. ConfSearch passes one file per run, so a temperature stage does not rebuild what an earlier repetition already tried -- measured: repetitions 2 and 3 of a 600 K stage proposed the same structures again, and 32 of 113 proposals across a 7-cycle run were repeats.", "Generation", {})
    PARAM(consensus_max, Int, 3, "Number of de-novo assemblies built: the pure consensus plus variants that flip the torsions whose best and second-best state are closest in energy (the least certain choices).", "Generation", {})

    END_PARAMETER_DEFINITION
    // ^^^^^^^^^^^^ PARAMETER DEFINITION BLOCK ^^^^^^^^^^^^
};
