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
    PARAM(cv_folds, Int, 5, "Number of cross-validation folds for the model comparison (additive vs. additive+couplings). Values below 2 disable the comparison.", "Analysis", {})
    PARAM(couplings, Bool, true, "Measure torsion-torsion couplings from double-mutant cycles (four ensemble members forming a rectangle in state space) and run the model comparison.", "Analysis", {})
    PARAM(report_threshold, Double, 1.0, "Only transitions whose mean total energy difference exceeds this many kJ/mol are printed in the summary. All of them are written to the CSV.", "Analysis", {})

    END_PARAMETER_DEFINITION
    // ^^^^^^^^^^^^ PARAMETER DEFINITION BLOCK ^^^^^^^^^^^^
};
