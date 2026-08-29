/*
 * <Conformational Search based on Molecular Dynamics>
 * Copyright (C) 2022 Conrad Hübler <Conrad.Huebler@gmx.net>
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
 */

#pragma once

#include <chrono>
#include <functional>
#include <limits>
#include <string>
#include <vector>

#include "src/tools/general.h"

#include "src/core/molecule.h"
#include "src/core/intra_parallel_context.h"

#include "json.hpp"

#include "external/CxxThreadPool/include/CxxThreadPool.hpp"

#include "src/capabilities/curcumamethod.h"
#include "src/capabilities/shared_bias_pool.h"
#include "src/capabilities/optimizer_factory.h"
#include "src/core/energycalculator.h"
#include "src/core/parameter_macros.h"
#include "src/core/parameter_registry.h"
#include "src/core/config_manager.h"
using namespace curcuma;

// Claude Generated (Apr 2026): Thread class for parallel geometry optimization
class OptThread : public CxxThread {
public:
    OptThread(const Molecule& molecule, const json& parameter)
        : m_molecule(molecule), m_parameter(parameter)
    {
        setAutoDelete(false);
    }

    /* Claude Generated (Jul 2026): completion hook for the live batch progress display.
     * Invoked from the worker thread when this structure is done, with its wall-time and its
     * result. Formatting lives in the caller (ConfSearch::PerformOptimisation) because the global
     * logger level is unreliable inside pool workers -- the same reason MDThread has this hook.
     * Without it the batch was completely silent between "Optimizing N structures" and the final
     * table, so a long batch gave no sign of which structure it was on. */
    void setOnComplete(std::function<void(double, const Optimization::OptimizationResult&)> cb)
    {
        m_on_complete = std::move(cb);
    }

    virtual int execute() override
    {
        const auto t0 = std::chrono::steady_clock::now();
        // This OptThread runs as one task among many under a molecule-level pool.
        // Suppress intra-molecule fan-out so methods that honor the flag (native
        // gfn1/gfn2) stay serial and the cores are not oversubscribed.
        curcuma::SuppressIntraParallel intra_guard;

        // Sync the global logger level with the json verbosity so GFN-FF init,
        // EnergyCalculator setup, and optimizer messages respect the requested level.
        // PerformOptimisation() restores the parent level after the pool finishes.
        CurcumaLogger::set_verbosity(m_parameter.value("verbosity", 0));

        std::string method = m_parameter.value("method", std::string("gfnff"));
        EnergyCalculator energy_calc(method, m_parameter);
        m_result = Optimization::OptimizationDispatcher::optimizeStructure(
            &m_molecule, Optimization::OptimizerType::LBFGSPP, &energy_calc, m_parameter);
        // Claude Generated (Jul 2026): did the METHOD itself fail on this structure (NaN energy or
        // gradient)? The OptimizationResult cannot say -- a NaN gradient reaches it as an ordinary
        // "did not converge". The calculator knows, and it dies with this thread, so ask it here.
        m_method_nan = energy_calc.HasNan() || energy_calc.Error();
        if (m_on_complete)
            m_on_complete(std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(), m_result);
        return m_result.success ? 0 : 1;
    }

    const Optimization::OptimizationResult& result() const { return m_result; }

    /** True when the energy method reported NaN/Inf for this structure (see execute()). */
    bool methodReportedNaN() const { return m_method_nan; }

private:
    Molecule m_molecule;
    json m_parameter;
    Optimization::OptimizationResult m_result;
    bool m_method_nan = false;
    std::function<void(double, const Optimization::OptimizationResult&)> m_on_complete;
};

// Claude Generated (May 2026, ICX-build): forward decl must be inside the namespace.
// ICX rejects `class curcuma::Molecule;` as a nested-name forward decl; GCC accepts it.
namespace curcuma { class Molecule; }

class ConfSearch : public CurcumaMethod {
public:
    ConfSearch(const json& controller, bool silent);
    ~ConfSearch();

    void setFile(const std::string& file);
    virtual bool Initialise() override;

    virtual void start() override;

    virtual void printHelp() const override
    {
        std::cout << "Usage: curcuma -confsearch input.xyz [parameters]\n\n";
        ParameterRegistry::getInstance().printHelp("confsearch");
        std::cout << "\nRMSD-MTD is enabled by default. Use -rmsd_mtd false to disable.\n";
    }

private:
    void PerformMolecularDynamics(const std::vector<Molecule*>& molecules, const nlohmann::json& parameter);

    /* Claude Generated (Aug 2026): `out_stem` gives the result its own name. Without it the output
     * is "<in_stem>.opt.xyz", which chained every step of a cycle onto the FIRST step's stem --
     * the deduplicated minima ended up in "...s1_explore.<method>.opt.accepted.xyz" although they
     * are the product of step 3. Each step now writes "<step>.<method>.xyz". */
    std::string PerformOptimisation(const std::string& filename, const nlohmann::json& parameter,
        const std::string& out_stem = std::string());

    /* Reads "<in_stem>.xyz"; ConfScan writes "<in_stem>.accepted.xyz", which is renamed to
     * "<out_stem>.xyz" when one is given (see PerformOptimisation). */
    std::string PerformFilter(const std::string& filename, const nlohmann::json& parameter,
        const std::string& out_stem = std::string());

    /* Claude Generated (Jul 2026): Phase 3c -- torsion recombination on this cycle's deduplicated
     * minima. Runs ConfGen on "<f>.xyz", which enumerates torsion-state vectors the ensemble does
     * NOT contain, builds them, optimises them at `method` and keeps the chemically valid, genuinely
     * new ones. Those are appended to "<f>.xyz", so everything downstream (Phase 3b re-optimisation,
     * the energy/topology filter, the bias feedback and the seed selection in Phase 4) treats them
     * exactly like a structure the metadynamics had found -- no special case anywhere.
     * @return number of new conformers appended */
    int PerformConfGen(const std::string& filename, const std::string& method);

    /* Claude Generated (Jul 2026): single source of truth for every child computation.
     * Every MD run, every optimisation batch and every nested ConfScan draws its config from
     * here, so the system identity (charge/spin), the runtime settings (threads/gpu/verbosity)
     * and the method sub-scopes (gfnff/xtb/tblite/...) can never silently differ between phases.
     * Replaces the ad-hoc {method, threads, gpu} JSONs that dropped charge, spin and solvation. */
    nlohmann::json ChildConfig(const std::string& method, int threads) const;

    /* Claude Generated (Aug 2026): cross-child topology lock. Writes the reference geometry,
     * runs one silent single point so gfnff derives and caches its topology, and returns the
     * .topo.json path (empty when locking does not apply to this method or the file failed).
     * ChildConfig then hands the path to every child via -gfnff.topology_file. */
    std::string PrepareTopologyLock(const std::string& method, const Molecule& ref, const std::string& stage);

    /* Claude Generated (Jul 2026): config for a nested ConfScan filter pass. Built from the
     * confscan registry defaults plus explicit, intentional overrides -- NOT from ConfSearch's
     * own defaults, which used to leak "method":"gfnff" into ConfScan's RMSD-alignment parameter
     * and pinned the dedup threshold to a ConfSearch default instead of the user's -rmsd. */
    nlohmann::json FilterConfig(const std::string& energy_method, int threads) const;

    /* Claude Generated (Jul 2026): Phase 3b -- re-optimisation of the cycle's deduplicated minima at
     * opt_method. Single-stage by default (one accurate optimisation per structure). With
     * phase3b_two_stage the work is split: a crude optimisation of EVERY structure first, then a
     * ConfScan dedup at that level (structures that collapse onto the same minimum during the crude
     * relaxation are removed before they are paid for), then the accurate optimisation of the
     * survivors only.
     * @return path of the file holding the opt_method ensemble (empty when nothing was produced) */
    std::string PerformHighLevelOptimisation(const std::string& basename);

    /* Claude Generated (Jul 2026): file stem of one stage, "<basename>.<purpose>.<method>".
     * Every intermediate file a run writes says WHAT it is and WHICH method produced it --
     * "input.initial.gfnff.opt.xyz" instead of "input.input.opt.xyz". The ".opt" / ".opt.accepted"
     * suffixes are still appended by PerformOptimisation / PerformFilter, so the chain stays
     * readable: <base>.explore.gfnff.opt.accepted.xyz is "the deduplicated gfnff minima of the
     * exploration phase". The final result keeps its established name
     * (<base>.cumulative.opt.accepted.xyz) -- that one is the deliverable, not an intermediate. */
    std::string stageBase(const std::string& purpose, const std::string& method) const
    {
        return Basename() + "." + purpose + "." + method;
    }

    /* Same, but for a file that belongs to ONE temperature cycle:
     * "<basename>.cycleNN_TxxxK.<purpose>.<method>". The cycle tag comes first so a directory
     * listing groups every file of a cycle together. Before this, all per-cycle working files were
     * overwritten by the next cycle -- what a given temperature actually did was gone. */
    std::string cycleStage(const std::string& purpose, const std::string& method) const
    {
        return Basename() + "." + m_cycle_tag + "." + purpose + "." + method;
    }

    /* Claude Generated (Jul 2026): drop MD snapshots whose bond topology differs from the reference,
     * BEFORE they are optimised.
     *
     * A conformer search explores one molecule; a snapshot that formed or broke a bond is not a
     * conformer of it and is rejected by the Phase 4 topology filter anyway -- after a full
     * optimisation has been paid for it. Worse, those geometries are the ones that make GFN-FF
     * return "energy is finite but the gradient contains NaN" (overlapping atoms in the 1/r hb and
     * three-body terms), which floods the log and wastes the optimiser's iterations.
     *
     * Rewrites the file in place with the surviving frames and reports what was dropped, split into
     * formed vs broken bonds (the two failure modes have different causes: a collision vs a
     * fragmentation).
     *
     * The criterion is the same covalent-radius rule as Phase 4 (factor 1.5), which is generous
     * enough that thermal bond stretching at the temperatures used here does not trip it -- a C-C
     * bond has to exceed 2.28 A to count as broken.
     *
     * @return number of frames kept */
    int FilterSnapshotsByTopology(const std::string& path) const;

    /** Claude Generated (Aug 2026): drop snapshots that are within -snapshot_dedup_rmsd of one
     *  already kept, BEFORE the optimisation. Cost saving, not deduplication -- see the definition. */
    int ThinSnapshots(const std::string& path) const;

    /** Claude Generated (Aug 2026): drop structures that changed SPECIES during the RELAX
     *  optimisation, right after it instead of at the end of the cycle. See the definition. */
    int FilterOptimisedByTopology(const std::string& path);

    /* Claude Generated (Jul 2026): try to REPAIR a snapshot instead of discarding it.
     *
     * A snapshot that broke a bond is not a conformer -- but its CONFORMATION may still be one, and
     * with GFN-FF a "broken" bond is often the force field losing a contact near its cutoff rather
     * than real chemistry. Instead of throwing the structure away, the offending atom pairs are
     * restrained back to a sane distance (a missing bond to the reference bond length, a spurious
     * contact to just outside bonding range), the structure is optimised under those restraints,
     * they are RELEASED, and it is optimised freely. It is kept only if the topology then matches --
     * the restraint is a way to reach the structure, never a licence to skip the check.
     *
     * Same mechanism as polymerbuild's interface-bond penalty (polymerbuild.cpp:181), as an energy.
     *
     * @param mol         snapshot, replaced by the repaired geometry on success
     * @param calculator  shared calculator (one instance for the whole batch -- creating a second
     *                    GFN-FF instance for the same molecule in one process has crashed before)
     * @return true when the repaired structure carries the reference topology */
    bool RepairSnapshot(Molecule& mol, EnergyCalculator& calculator) const;

    /* Claude Generated (Jul 2026): 0/1 bond matrix with an EXPLICIT covalent-radius factor.
     *
     * Replaces Molecule::DistanceMatrix().second, whose factor of 1.5 is far too generous to decide
     * whether a structure is still the same molecule: a compressed 1-3 contact counts as a bond.
     * Measured on the 107-atom peptide this was found on -- two carbons bonded to the same atom sit
     * 2.25 A apart, and the 1.5 criterion calls them BONDED at 2.28 A, i.e. as soon as their angle
     * closes from 107 deg to 109 deg. That is a tetrahedral angle: the "bond" appears and disappears
     * on a two-degree thermal fluctuation. Consequence, measured on a 10 ps 500 K trajectory:
     *   factor 1.5 -> 92 % of the frames "changed topology";  factor 1.3/1.4 -> 0 %.
     * The molecule never reacted. ConfGen already hit this and uses its own factor (1.3); ConfSearch
     * kept using 1.5 and therefore threw away ~90 % of perfectly good conformers in its Phase 4
     * topology filter. */
    Matrix TopologyMatrix(const Molecule& mol) const;

    /* Copy every frame of an XYZ file to another name. Used where a stage has to start from the
     * previous stage's OUTPUT: the copy is what lets the new stage own a file whose name states its
     * purpose and method instead of inheriting a chain of suffixes. Returns the frame count. */
    /* `append` keeps what is already in the destination -- used to accumulate the minima of every
     * repetition of a temperature stage into one pool (see the stage-repetition loop). */
    int CopyFrames(const std::string& source, const std::string& destination, bool append = false) const;

    /* Claude Generated (Jul 2026): energy summary of an XYZ ensemble on disk. Returns the energies
     * in ascending order (empty when the file is missing or holds no structures). */
    std::vector<double> EnsembleEnergies(const std::string& path) const;

    /* Claude Generated (Jul 2026): log "<label> [<method>]: N structures, min ... span ..." plus the
     * lowest `ensemble_report` conformers. Reads the file, so it always reports what is actually on
     * disk rather than what the phase believes it wrote. */
    void ReportEnsemble(const std::string& label, const std::string& method, const std::string& path) const;

    /* Claude Generated (Jul 2026): per-cycle ensemble output on one level of theory. Writes the
     * structures energy-sorted to "<base>.<cycle_tag>.<method>.xyz", appends the most stable one to
     * the "<base>.best_per_cycle.<method>.xyz" trajectory and reports both. Ordering is done by a
     * std::multimap keyed on the energy -- the container sorts, no explicit sort call. */
    void WriteCycleEnsemble(const std::string& cycle_tag, const std::string& method,
        const std::vector<Molecule*>& molecules) const;

    /* Claude Generated (Jun 2026): experimental adaptive bias calibration (Phase C). Clusters the
     * cycle's optimised structures onto the accepted distinct minima (best-fit RMSD + cached
     * permutations + energy tolerance), then derives the basin capture radius and per-atom RMSF.
     * In "cluster" mode it sets md["rmsd_mtd_alpha"] from the basin radius; in "weighted" mode it
     * sets flexibility weights (1/sigma^2) on the shared pool for the RMSF-weighted MTD bias. */
    void CalibrateBias(const std::string& basename, nlohmann::json& md);

    /* Permutation-aware best-fit (Kabsch) RMSD in Angstrom: min over identity + cached symmetry
     * permutations applied to the target. Both geometries must share the canonical atom order. */
    double PermRMSD(const Geometry& reference, const Geometry& target) const;

    /* Claude Generated (Jul 2026): pick the seeds for the next temperature cycle out of the
     * structures that already passed the seed energy window. "energy" keeps the seed_rank lowest
     * ones; "diverse" walks the same energy ranking but skips a candidate that sits closer than
     * seed_min_rmsd (permutation-aware best-fit RMSD) to a seed already chosen, so MD time is not
     * spent four times on the same basin. Sorts the vector in place, deletes the rejected
     * molecules and shrinks the vector to the kept seeds. Returns the number of rejects. */
    /**
     * @brief How many seeds this temperature gets (see seed_rank_cold_factor).
     *
     * Claude Generated (Aug 2026). Linear in the temperature between startT and endT, so the
     * factor is readable off the endpoints: 1.0 keeps the old constant count exactly.
     */
    int EffectiveSeedRank() const;

    /**
     * @brief Accumulated metadynamics hill height around each candidate (see seed_bias_penalty).
     *
     * Claude Generated (Aug 2026). Returns, per candidate, the sum of the counters of all pool
     * hills within seed_bias_radius. The counter IS the hill height (W = k*counter), so this is
     * "how much bias already sits here" -- the quantity that decides whether a trajectory started
     * here would be pushed straight out again. Screened by the radius of gyration first
     * (|Rg(a) - Rg(b)| <= RMSD(a,b) is a valid lower bound), which is what makes it affordable
     * against a pool of tens of thousands of hills.
     */
    std::vector<double> SeedBiasDensity(const std::vector<Molecule*>& candidates) const;

    int SelectSeeds(std::vector<Molecule*>& window_seeds, const std::string& method,
        double global_min) const;

    /* Per-seed RMSD spacing of the last diverse selection (diagnostic output only). */
    mutable std::vector<double> m_seed_spacing;

    /* Lets have this for all modules */
    virtual nlohmann::json WriteRestartInformation() override;

    /* Lets have this for all modules */
    virtual bool LoadRestartInformation() override;

    // Claude Generated (Jun 2026): ConfSearch restart/checkpoint.
    // Self-contained checkpoint of the whole search state (bias pool, cumulative pool, seeds,
    // energy progress, cycle/phase, symmetry permutations, topology reference). See
    // docs/CONFSEARCH_DUAL_METHOD.md / CONFSEARCH_RESTART.md.
    struct RestartState {
        bool valid = false;
        int entry_phase = 0; // where to re-enter the resumed cycle: 0=md,1=post_md,2=post_filter,3=post_refine
        double next_T = 0;   // temperature of the cycle to (re)start
        int temperature_cycle = 0; // number of cycles already completed
        int natoms = 0;
        std::string md_method, opt_method;
        double global_min = 0, best_energy = 0, initial_energy = 0;
        double initial_energy_opt = 0; // opt_method energy of the input structure (dual-method runs)
        std::vector<int> elements; // atomic numbers (shared by all frames)
        std::vector<BiasStructure> bias; // full bias pool (geometry + metadata)
        std::vector<Molecule> seeds; // m_in_stack for the resumed cycle
        std::vector<Molecule> cumulative; // accepted conformers from completed cycles
        std::vector<Molecule> accepted_md; // gfnff-filtered set (post_filter/post_refine)
        std::vector<Molecule> accepted_opt; // opt_method-reoptimised set (post_refine)
        Molecule topo_ref; // reference structure -> m_topo_matrix
        Molecule topo_ref_opt; // opt_method reference structure -> m_topo_matrix_opt (dual-method)
        std::vector<std::vector<int>> permutations; // m_permutation_cache
    };
    // Serialise one molecule (geometry + energy + name) given the shared element list is stored once.
    nlohmann::json molToJson(const Molecule& mol) const;
    nlohmann::json molVectorToJson(const std::vector<Molecule>& mols) const;
    nlohmann::json molPtrVectorToJson(const std::vector<Molecule*>& mols) const;
    nlohmann::json fileFramesToJson(const std::string& path) const; // read an xyz file into the geometry-only json form
    Molecule jsonToMol(const std::vector<int>& elements, const nlohmann::json& entry) const;
    std::vector<Molecule> jsonToMolVector(const std::vector<int>& elements, const nlohmann::json& arr) const;
    void writeMolVectorToFile(const std::vector<Molecule>& mols, const std::string& path) const; // first writes, rest append
    // Build the full checkpoint json (uses the m_ckpt_* staging members), write it to the BMT dir
    // and copy it back to the start directory (CWD) under the stable name.
    void writeCheckpoint(const std::string& phase, double next_T, int temperature_cycle);
    // Read + restore the CWD checkpoint into m_restart; returns true on a valid, matching checkpoint.
    bool loadCheckpoint();
    std::string restartFileName() const; // "<basename>.confsearch.restart.json"

    virtual StringList MethodName() const override
    {
        return { "ConfSearch" };
    }

    /* Lets have all methods read the input/control file */
    virtual void ReadControlFile() override;

    /* Read Controller has to be implemented for all */
    virtual void LoadControlJson() override;

    // Claude Generated (Jul 2026): registry-backed parameter access. ConfigManager resolves
    // intra-module aliases (config_manager.cpp:95), which plain Json2KeyWord(m_defaults, ...)
    // does not -- without it "-dt", "-velo", "-dump" etc. would silently stop working after
    // the canonical names were aligned with SimpleMD.
    ConfigManager m_config;

    StringList m_error_list;
    std::string m_method, m_md_method, m_opt_method, m_thermostat;
    bool m_silent = true;
    double m_dT = 4;
    std::vector<Molecule*> m_in_stack, m_final_stack;
    int m_spin = 0, m_charge = 0, m_repeat = 5, m_threads = 1, m_max_bias_export = 1000;
    double m_time = 1e4, m_startT = 500, m_endT = 300, m_deltaT = 50, m_currentT = 0, m_rmsd = 1.25, m_energy_window = 100;
    // Claude Generated (Jun 2026): efficiency/robustness controls (see the PARAM block below for meaning)
    double m_rattle_threshold_temp = 400, m_seed_energy_window = 50, m_seed_window_decay = 0.5, m_epot_abort_window = 250;
    int m_seed_rank = 1; // max lowest-energy seeds per cycle (0 = all in window; 1 = only most stable)
    int m_rattle_hot_mode = 2, m_topo_check_interval = 0, m_opt_feedback_height = 5;
    // Claude Generated (Aug 2026): surface-mismatch feedback + the manual bridge threshold.
    bool m_surface_feedback = false, m_hbond_excess_reject = false;
    int m_surface_feedback_min_structures = 20, m_hbond_excess_max = -1, m_reference_hbonds = -1;
    double m_surface_feedback_min_r = 0.4, m_surface_feedback_strength = 3.0;
    bool m_confgen_phase = false;
    int m_confgen_max_proposals = 20, m_confgen_templates = 3, m_confgen_depth = 3;
    // Claude Generated (Aug 2026): surface that judges ConfGen's proposals, see PerformConfGen().
    std::string m_confgen_eval_method = "none";
    double m_snapshot_dedup_rmsd = 0.0;   ///< pre-optimisation snapshot screen, see ThinSnapshots()
    std::string m_bias_reset = "never";
    double m_reduce_prefilter_window = 100.0;
    double m_reduce_prefilter_energy_tol = 1.0e-6;
    /** Cheap pre-filter between RELAX and REDUCE: energy window + identical-minimum collapse.
     *  @return number of structures written to `out` */
    int PrefilterForReduce(const std::string& in_file, const std::string& out_file) const;
    bool m_bias_rejected = false;
    /* Claude Generated (Aug 2026): optional per-repetition schedules for the two bias
     * parameters (see mtd_alpha_schedule). Empty = one value for the whole run. */
    std::vector<double> m_mtd_alpha_schedule, m_mtd_k_schedule;
    double m_verify_best = 1.0;   ///< tolerance in kJ/mol for the recomputed cycle minimum (0 = off)
    /** Recompute the energy of `mol` with a fresh calculator at `method`.
     *  @return the recomputed energy, or NaN if it could not be obtained. */
    double RecomputeEnergy(const Molecule& mol, const std::string& method) const;

    /* Claude Generated (Aug 2026): the two answers to a measured problem -- the exploration surface
     * and the ranking surface disagree systematically, and the exploration surface then drives the
     * search into a region the ranking surface rejects (measured: GFN-FF rewards every additional
     * hydrogen bond, so the cold cycles produced 10-11-bridge folds that GFN-FF put 230 kJ/mol below
     * everything else and GFN2 put 52-95 kJ/mol above its own best).
     *
     * MeasureSurfaceMismatch does NOT know that hydrogen bonds are the coordinate -- it measures
     * which of four generic shape descriptors correlates with the disagreement and returns a flag
     * per structure, so the same code finds a different coordinate on a different system. */
    std::vector<char> MeasureSurfaceMismatch(const std::vector<Molecule*>& structures) const;
    int HydrogenBondCount(const Molecule& mol) const;
    bool m_hold_polar_h = false;
    double m_hold_polar_h_force = 5.0;
    /** Harmonic distance restraints on the polar X-H bonds of a reference structure (json form).
     *  @param ref  the reference for THIS surface -- the md-level one for step 2 RELAX, the
     *              opt-level one for step 5 REFINE, because each side's topology filter
     *              compares against its own reference and the bond lengths differ. */
    nlohmann::json PolarHydrogenRestraints(const Molecule& ref) const;
    mutable bool m_polar_h_reported = false; ///< report the restraint list once, not per stage
    // Claude Generated (Aug 2026): the two additional move sets of Phase 3c.
    bool m_confgen_nci_moves = true, m_confgen_consensus = false;
    std::string m_confgen_method = "auto";
    /// Resolved method for the recombination step (see confgen_method).
    std::string ConfGenMethod() const;

    /**
     * @brief The method the per-repetition funnel computes with (RELAX, REDUCE, re-scoring).
     *
     * Claude Generated (Aug 2026). Everything between the MD and the stage pool -- optimising a
     * snapshot, ranking it inside the energy window, deciding whether it duplicates another
     * structure -- is one decision made on ONE surface, and this returns that surface. The MD
     * itself always runs at md_method; only what happens to its snapshots changes.
     * Identical to md_method in a single-method run, so nothing there is affected.
     */
    std::string RelaxMethod() const { return m_relax_pes == "opt" ? m_opt_method : m_md_method; }

    /// True when the funnel already produced opt_method minima, so REFINE has nothing left to do.
    bool FunnelIsRanking() const { return RelaxMethod() == m_opt_method; }

    /**
     * @brief Topology reference for structures the funnel produced.
     *
     * Claude Generated (Aug 2026). Two methods do not agree on bond lengths, so an atom pair near
     * the covalent-radius cutoff can be bonded for one and not for the other -- measured, checking
     * opt_method geometries against the md_method reference rejected an entire cycle. The funnel's
     * structures are optimised at RelaxMethod(), so they are compared against that method's own
     * optimised input. The MD snapshots keep m_topo_ref: they come from an md_method trajectory.
     */
    const Molecule& FunnelTopoRef() const
    {
        return (FunnelIsRanking() && m_topo_ref_opt.AtomCount() > 0) ? m_topo_ref_opt : m_topo_ref;
    }
    /// Bond matrix belonging to FunnelTopoRef().
    const Matrix& FunnelTopoMatrix() const
    {
        return (FunnelIsRanking() && m_topo_matrix_opt.rows() == m_topo_matrix.rows())
            ? m_topo_matrix_opt
            : m_topo_matrix;
    }

    /**
     * @brief Every structure the search has accepted, on each surface -- the pool the seeds come from.
     *
     * DESIGN, not a patch: the seeds of a cycle are the `seed_rank` best structures of ALL cycles,
     * every accepted structure is a bias hill, and a temperature stage is repeated via `repeat`.
     * Drawing the seeds from the current cycle alone breaks the first rule and was measured to matter:
     * a structure found once goes into the cumulative pool and into the bias pool as a PERSISTENT
     * hill, and persistent hills are deliberately never re-optimised, so it never re-enters a cycle
     * ensemble -- and therefore never again the seed list. On WEKLQ (7 cycles, seed_rank 10) cycle 3
     * found the global minimum -161.641649, cycle 4 still seeded from it, and from cycle 5 the search
     * seeded from structures 6-7 kJ/mol worse while the global minimum sat 3.1-3.5 A from anything
     * those cycles contained. The result was never lost -- only the seeding drifted downhill.
     *
     * The energy window (`seed_energy_window`) is applied at SELECTION time, not on insertion,
     * because the running minimum only ever drops: a structure inside the window today can be
     * outside it tomorrow, and the pool must not freeze that decision.
     */
    std::vector<Molecule> m_seed_pool_md, m_seed_pool_opt;

    /** One row per temperature stage for the closing summary -- what each stage actually
     *  contributed. Without it the only way to tell whether a cycle was worth its runtime was to
     *  read back through hundreds of log lines. */
    struct StageSummary {
        int cycle = 0;
        double temperature = 0.0;
        int repetitions = 0;
        int structures = 0;      ///< minima the stage handed to the accurate method
        int recombined = 0;      ///< of those, produced by the recombination step
        double best_md = 0.0;    ///< lowest md_method energy of the stage
        double best_opt = 0.0;   ///< lowest opt_method energy of the stage (NaN in single-method runs)
        double seconds = 0.0;
        bool new_global_best = false;
    };
    std::vector<StageSummary> m_stage_summary;
    void ReportStageSummary() const;

    /// Add this cycle's candidates to the cross-cycle pool (duplicates by energy dropped).
    void AccumulateSeedPool(std::vector<Molecule>& pool, const std::vector<Molecule*>& candidates) const;
    /// Offer the pool as seed candidates, filtered by the energy window (ownership passes on).
    int OfferSeedPool(std::vector<Molecule*>& seeds, const std::vector<Molecule>& pool,
        double global_min, double window) const;
    // Claude Generated (Jul 2026): two-stage Phase 3b (crude opt -> dedup -> accurate opt) and the
    // per-cycle ensemble/energy reporting.
    bool m_phase3b_two_stage = false, m_phase3b_filter = true, m_cycle_output = true;
    std::string m_phase3b_preopt_preset = "loose", m_phase3b_preset = "normal";
    int m_phase3b_preopt_max_iter = 0, m_ensemble_report = 3;
    bool m_snapshot_topology_gate = true; // Claude Generated (Jul 2026): pre-Phase-2 topology gate
    double m_snapshot_clash_ratio = 0.55;  // collapsed-contact screen, fraction of the covalent sum
    double m_topology_factor = 1.3;        // covalent-radius factor of every topology comparison
    // Claude Generated (Jul 2026): restrained repair of near-miss snapshots (see RepairSnapshot)
    bool m_repair_snapshots = false;
    bool m_topology_lock = true; ///< Claude Generated (Aug 2026): one gfnff parameterisation per search
    int m_stage_saturation_abort = 2; ///< Claude Generated (Aug 2026): dry repetitions before a stage is cut short
    int m_rmsd_mtd_cap_eff = 0; ///< Claude Generated (Aug 2026): the bias cap actually in force (resolved from -1 = adaptive)
    std::string m_topo_lock_md, m_topo_lock_opt; ///< reference .topo.json per PES (empty = unlocked)
    int m_repair_max = 20, m_repair_max_bonds = 2, m_repair_max_iterations = 300;
    double m_repair_force = 2.0;
    bool m_topo_check = false, m_epot_abort = false, m_opt_feedback_bias = true, m_opt_feedback_prune_snapshots = false, m_mtd_permutation = true;
    // Claude Generated (Jun 2026): temperature runaway abort + cross-run bias-height freeze.
    // ON by default for ConfSearch (bias-heating safety net + best conformer yield); see the PARAM block below.
    bool m_temp_abort = false, m_freeze_inherited = false;
    // Claude Generated (Jun 2026): initial energy at opt_method (dual-mode only)
    double m_initial_energy_opt = std::numeric_limits<double>::infinity();
    double m_temp_abort_factor = 1.5, m_temp_abort_delta = 300;
    int m_rmsd_mtd_max_height = 0;
    // Claude Generated (Jul 2026): bias-evaluation speedup controls (Gaussian-cutoff screen + pool cap)
    int m_rmsd_mtd_max_gaussians = -1;
    bool m_rmsd_mtd_screen = true;
    double m_rmsd_mtd_cutoff_tol = 1.0e-8, m_rmsd_mtd_screen_margin = 0.0;
    std::string m_seed_window_schedule = "static";
    // Claude Generated (Jul 2026): RMSD-aware seed selection (see SelectSeeds).
    std::string m_seed_selection = "diverse";
    /* Claude Generated (Jul 2026): which PES selects the next cycle's seeds ("md" | "opt"). */
    std::string m_seed_pes = "md";
    /* Claude Generated (Aug 2026): divergence guard for optimisation results, see
     * opt_divergence_factor. A diverged geometry is finite, so allFinite() does not catch it. */
    double m_opt_divergence_factor = 10.0;
    /* Claude Generated (Aug 2026): which PES the per-repetition funnel runs on ("md" | "opt").
     * See the relax_pes parameter for the measurement behind the "opt" default. */
    std::string m_relax_pes = "opt";
    double m_seed_min_rmsd = 0.0, m_seed_diversity_factor = 2.0;
    /* Claude Generated (Aug 2026): temperature-dependent seed count and bias-density-aware
     * ranking -- see seed_rank_cold_factor / seed_bias_penalty. Both default to "off". */
    double m_seed_rank_cold_factor = 1.0, m_seed_bias_penalty = 0.0, m_seed_bias_radius = 0.0;

    std::string m_bias_calibration = "off"; // adaptive MTD width mode: off | couple | cluster
    std::string m_bias_scale_mode = "global"; // global | weighted (RMSF-weighted RMSD)
    double m_bias_couple_factor = 1.0, m_bias_energy_tol = 4.0;
    double m_global_min = std::numeric_limits<double>::infinity(); // running lowest energy across all cycles
    double m_global_min_opt = std::numeric_limits<double>::infinity(); // the same on the opt_method PES
    std::vector<std::vector<int>> m_permutation_cache; // Claude Generated (Jun 2026): symmetry reorder rules from ConfScan, fed into MTD
    /* Claude Generated (Jul 2026): "cycleNN_TxxxK" of the cycle currently running. Set at the top of
     * every temperature cycle and used by cycleStage() for every file that cycle writes. */
    std::string m_cycle_tag = "cycle00_T0K";
    Matrix m_topo_matrix;
    /* Claude Generated (Jul 2026): SECOND topology reference, for the opt_method PES.
     * The refinement side used to check the opt_method geometries against the md_method reference.
     * Two methods do not agree on bond lengths, so a pair sitting near the covalent-radius cutoff
     * can be "bonded" for one and not for the other -- and then EVERY re-optimised structure of a
     * cycle is rejected as a reaction product (observed: 38 of 38 gfn2 structures lost, the whole
     * ranking side of the cycle silently discarded). Empty until the initial opt_method optimisation
     * has run; the refinement side falls back to m_topo_matrix while it is. */
    Matrix m_topo_matrix_opt;
    Molecule m_topo_ref_opt;
    SharedBiasPool* m_bias_pool = nullptr;  // Claude Generated (Apr 2026): shared bias pool for parallel ConfSearch
    // Claude Generated (Jun 2026): restart/checkpoint state.
    bool m_do_restart = false;                  // -restart: enable checkpoint write + resume
    RestartState m_restart;                     // populated by loadCheckpoint() on resume
    double m_best_energy = std::numeric_limits<double>::infinity();    // exploration best (md_method), persisted
    double m_initial_energy = std::numeric_limits<double>::infinity(); // cycle-1 exploration reference, persisted
    std::string m_cumulative_file;
    /** Claude Generated (Aug 2026): md-side collection of every topology-valid minimum of the run,
     *  appended per cycle and handed to the recombination as its description basis. Deliberately
     *  unfiltered -- the state and contact statistics need the higher-lying structures too, and the
     *  seed pool (which the energy window truncates at 50 kJ/mol) is the wrong source for that. */
    std::string m_cumulative_md_file;              // path to the cumulative output (set in start())
    std::vector<int> m_elements;                // atomic numbers of the system (set after pre-opt / on resume)
    Molecule m_topo_ref;                        // reference structure defining m_topo_matrix (persisted)
    std::string m_gpu = "none";                 // Claude Generated (Jul 2026): GPU backend, forwarded to every child

    // vvvvvvvvvvvv PARAMETER DEFINITION BLOCK vvvvvvvvvvvv
    // Claude Generated (Jul 2026): ConfSearch migrated off the static ConfSearchJson literal.
    // Owning these names is what stops main.cpp's auto-router from moving -opt_method into
    // polymerbuild, -thermostat/-coupling/-rattle/-wall_* into simplemd and -restart into
    // confscan. See docs/CONFSEARCH_DUAL_METHOD.md.
    //
    // DELIBERATELY NOT REGISTERED (they are registered by NO module today; claiming them here
    // would steal them from "-md", where SimpleMD does consume them):
    //   unique, rescue, respa, dipole, cleanenergy, MaxTopoDiff, impuls_scaling,
    //   rattle_tolerance, printOutput, wall_xl/yl/zl
    // They still work as unregistered flat flags and are forwarded to the md json as before.
    //
    // Every Double MUST carry a decimal-point literal: the generator emits the token verbatim,
    // an integer literal makes std::any hold int, getDefaultJson's any_cast<double> throws and
    // the key is SILENTLY dropped from the defaults -> Json2KeyWord then throws an uncaught -1.
    BEGIN_PARAMETER_DEFINITION(confsearch)

    // --- Methods ---
    PARAM(method, String, "gfnff", "Energy method used for both phases unless md_method or opt_method override it.", "Methods", {})
    PARAM(md_method, String, "", "Cheap method driving MD exploration and pre-optimisation. Empty falls back to method.", "Methods", {})
    PARAM(opt_method, String, "", "Accurate method for the per-cycle re-optimisation and the final ranking. Empty falls back to method.", "Methods", {})

    // --- System ---
    PARAM(charge, Int, 0, "Total charge of the system. Applied to every MD run, every optimisation and every ConfScan energy.", "System", {})
    PARAM(spin, Int, 0, "Total spin of the system.", "System", {"Spin"})
    PARAM(threads, Int, 1, "Number of ConfSearch cycles run in parallel. Each child MD or optimisation then runs single-threaded.", "System", {})

    // --- Search Schedule ---
    PARAM(startT, Double, 600.0, "Temperature of the first exploration cycle in Kelvin.", "Schedule", {})
    PARAM(endT, Double, 300.0, "Temperature of the last exploration cycle in Kelvin.", "Schedule", {})
    PARAM(deltaT, Double, 50.0, "Temperature decrement between cycles in Kelvin.", "Schedule", {})
    PARAM(repeat, Int, 4, "How often a temperature stage is repeated. Each repetition runs ONE MD per seed structure and then re-seeds from the best structures found so far, so repetition 2 starts from what repetition 1 discovered. The high-level re-optimisation runs once per temperature, at the end of the last repetition. Parallelism comes from the number of seeds (seed_rank), not from re-running one seed.", "Schedule", {})
    PARAM(time, Double, 2000.0, "Length of each MD run in femtoseconds.", "Schedule", {"max_time", "MaxTime"})

    // --- Filtering ---
    PARAM(rmsd, Double, 1.25, "RMSD threshold in Angstrom used to deduplicate conformers and to size the MTD hills.", "Filtering", {})
    PARAM(energy_window, Double, 100.0, "Energy window in kJ/mol above the running minimum for keeping conformers.", "Filtering", {})
    PARAM(seed_rank, Int, 1, "Maximum number of lowest-energy seeds carried into the next cycle. 0 keeps every structure inside seed_energy_window.", "Filtering", {})
    PARAM(seed_energy_window, Double, 50.0, "Energy window in kJ/mol versus the running global minimum for selecting next-cycle MD seeds.", "Filtering", {})
    PARAM(seed_pes, String, "md", "Which potential-energy surface picks the seeds of the next cycle in a dual-method run: md = the cheap exploration method (default; a basin it likes is never dropped because the accurate method ranks it higher), opt = the accurate ranking method, both = seed_rank seeds from EACH (an opt seed in the same basin as an md seed is dropped, since the opt structures are the re-optimised md ones). Motivated by measurement: on a 107-atom peptide the two surfaces correlate at only r = 0.40 and the gfn2 minimum sat at gfnff rank 59 of 141, so a gfnff-picked seed says little about where the accurate method has its minima. The MD still runs at md_method; only the starting geometries change. No effect in a single-method run.", "Filtering", {})
    PARAM(opt_divergence_factor, Double, 10.0, "An optimisation whose result spans more than this factor times its INPUT structure (never below 50 Angstrom, so small molecules keep room) is treated as diverged: the geometry is not written into the pool, the input geometry is used instead, and both are saved for inspection (<stage>.diverged.xyz for the blown-up result, <stage>.failed.xyz for the input that produced it). Why this exists next to the non-finite guard: a diverged optimisation is FINITE. Measured on a 107-atom peptide, gfn2 re-scoring of the recombination proposals: 72 of 2097 structures came back with coordinates of order 1e242 -- every isfinite() test passes, and the written energy (-161.62 Eh) sits in the middle of the conformer range, so such a structure was reported as the deepest of its cycle. Only the later species check removed it, after a full optimisation had been paid for. 0 disables the check.", "Optimisation", {})
    PARAM(relax_pes, String, "opt", "Which potential-energy surface the per-repetition selection funnel runs on in a dual-method run -- the RELAX optimisation of every snapshot, the energy window before the deduplication, the deduplication itself and the re-scoring of the recombination proposals. opt = the accurate ranking method (default); md = the cheap exploration method (the behaviour before Aug 2026). No effect in a single-method run, where both are the same surface. Why the default changed: selecting on the exploration surface is measurably wrong for this class of system -- the same reason -phase3b_two_stage carries, one level lower. -phase3b_two_stage only moved the ONCE-PER-STAGE re-optimisation onto the ranking surface; the funnel in front of it kept discarding on the exploration surface, once per repetition. Measured on a 107-atom peptide (WEKLQ): three gfnff/gfn2 runs stayed at +28.4/+28.6/+38.4 kJ/mol against a GOAT reference while gfn2/gfn2 reached -19.8, and the cause is visible in the budget -- the gfnff/gfn2 run spent 2345 QM optimisations, the gfn2/gfn2 run 7570. The target is NOT unreachable on the cheap surface (it sits at gfnff rank 5 of 649, a gfnff optimisation moves it 0.32 A, and the gfnff run came within 2.00 A of it); what the funnel discarded were the structures whose gfn2 basin was never determined. COST: with opt every snapshot is optimised on the accurate surface, so a dual-method run costs about what a single-method run of the accurate method costs, minus the cheaper dynamics. Set md to get the old, cheap behaviour back.", "Filtering", {"funnel_pes"})
    PARAM(seed_rank_cold_factor, Double, 1.0, "Scales seed_rank with the temperature: the effective number of seeds runs linearly from seed_rank at startT to seed_rank * this factor at endT. 1.0 (default) keeps the old constant behaviour exactly. Motivated by measurement on a 107-atom peptide: the colder the stage, the less the BEST seed contributes -- of the ten deepest structures per stage, the number NOT coming from the lowest-energy seed rises 1/10 (600 K), 4/10 (500 K), 5/10 (450 K), 7/10 (350 K), 10/10 (300 K). At 300 K the best seed produced none of them. Since a structure that once falls below rank seed_rank can never rise again (the pool only deepens), those basins are lost for the rest of the run. Raising the rank in the cold stages costs little, because the same stages produce far fewer snapshots anyway (measured 352 per repetition at 400 K against 113 at 300 K).", "Filtering", {})
    PARAM(seed_bias_penalty, Double, 0.0, "Energy penalty in kJ/mol applied to a seed CANDIDATE for the metadynamics hills already accumulated around it, used for the seed ranking only (never for reporting or for the ensemble). 0 (default) disables it and the selection is byte-identical to before. Why: the seed ranking is purely energetic, so the same deep basin is seeded again in every cycle -- and its surroundings fill with hills until the trajectory is pushed straight out of it. Measured on a 107-atom peptide: the lowest-energy seed dominates the early stages and contributes nothing to the ten deepest structures of the 300 K stage. The penalty is scaled by the accumulated hill height near the candidate (sum of counters within seed_bias_radius) relative to the most-covered candidate, so it only ever reorders candidates against each other.", "Filtering", {})
    PARAM(seed_bias_radius, Double, 0.0, "Radius in Angstrom within which bias hills count towards seed_bias_penalty. 0 derives it as the deduplication radius (rmsd), i.e. hills that the search itself would call the same structure. Only used when seed_bias_penalty > 0.", "Filtering", {})
    PARAM(seed_selection, String, "diverse", "How the next-cycle seeds are picked from the structures inside seed_energy_window: energy = strictly the seed_rank lowest-energy ones; diverse = lowest-energy first, then only structures at least seed_min_rmsd away (permutation-aware best-fit RMSD) from every seed already chosen.", "Filtering", {})
    PARAM(seed_min_rmsd, Double, 0.0, "Minimum RMSD in Angstrom between two seeds in the diverse selection. 0 derives it as seed_diversity_factor * rmsd.", "Filtering", {})
    PARAM(seed_diversity_factor, Double, 2.0, "Multiplier applied to rmsd when seed_min_rmsd is 0. Values around 2 keep the seeds one dedup radius apart from each other.", "Filtering", {})
    PARAM(seed_window_schedule, String, "static", "Seed window schedule: static or exp. exp shrinks the window each cycle by seed_window_decay.", "Filtering", {})
    PARAM(seed_window_decay, Double, 0.5, "Per-cycle multiplier applied to seed_energy_window in the exp schedule.", "Filtering", {})

    // --- Molecular Dynamics (canonical SimpleMD names, forwarded to the md json) ---
    PARAM(time_step, Double, 1.0, "MD integration time step in femtoseconds.", "MD", {"dt"})
    PARAM(thermostat, String, "csvr", "Thermostat: berendsen, andersen, nosehover, csvr or none.", "MD", {})
    PARAM(coupling, Double, 10.0, "Thermostat coupling time in femtoseconds.", "MD", {})
    PARAM(seed, Int, -1, "Random seed for the MD runs. -1 uses the clock.", "MD", {})
    PARAM(remove_com_motion, Double, 100.0, "Remove translation and rotation every N femtoseconds.", "MD", {"rm_COM"})
    PARAM(remove_com_mode, Int, 1, "Removal mode: 0 none, 1 translation, 2 rotation, 3 both.", "MD", {"rmrottrans"})
    PARAM(no_center, Bool, false, "Disable centering of the molecule at the origin.", "MD", {"nocenter"})
    PARAM(hydrogen_mass, Int, 1, "Hydrogen mass repartitioning factor.", "MD", {"hmass"})
    PARAM(initial_velocity_scale, Double, 1.0, "Scaling factor applied to the initial velocities.", "MD", {"velo"})

    // --- RATTLE ---
    PARAM(rattle, Int, 0, "RATTLE constraints in the baseline cycles: 0 off, 1 all bonds, 2 X-H only.", "RATTLE", {})
    PARAM(rattle_max_iterations, Int, 100, "Maximum RATTLE iterations per step.", "RATTLE", {"rattle_maxiter"})
    PARAM(rattle_threshold_temp, Double, 400.0, "Cycle temperature in Kelvin at or above which RATTLE is switched on automatically.", "RATTLE", {})
    PARAM(rattle_hot_mode, Int, 2, "RATTLE mode used above rattle_threshold_temp. 2 constrains X-H bonds only.", "RATTLE", {})

    // --- MD Output ---
    PARAM(dump_frequency, Int, 50, "Save MD coordinates every N steps.", "MD Output", {"dump"})
    PARAM(print_frequency, Int, 1000, "Print MD status every N steps.", "MD Output", {"print"})
    PARAM(write_xyz, Bool, true, "Write MD trajectories to XYZ files.", "MD Output", {"writeXYZ"})
    PARAM(write_initial_state, Bool, false, "Write the MD initial conditions to a .init.json file.", "MD Output", {"writeinit"})
    PARAM(restart_file, String, "none", "SimpleMD restart file used to seed the MD runs.", "MD Output", {"initfile"})
    PARAM(write_restart_frequency, Int, 1000, "Write the SimpleMD restart file every N steps.", "MD Output", {"writerestart"})
    PARAM(no_restart, Bool, false, "Disable automatic SimpleMD restart loading. ConfSearch manages its own checkpoint and forces this on regardless of the value here.", "MD Output", {"norestart"})

    // --- Robustness Gates ---
    PARAM(snapshot_topology_gate, Bool, true, "Drop MD snapshots whose bond topology differs from the reference structure BEFORE they are optimised. Such structures are not conformers of this molecule and are rejected by the Phase 4 filter anyway, after a full optimisation has been paid for them; they are also the geometries for which GFN-FF returns a finite energy together with a NaN gradient. Set false to optimise every snapshot as before.", "Robustness", {})
    PARAM(topology_factor, Double, 1.3, "Covalent-radius factor deciding whether two atoms count as bonded in EVERY topology comparison of the search (snapshot gate, Phase 4 filters, repair). Molecule's own default of 1.5 is far too generous here: a compressed 1-3 contact becomes a bond, so a two-degree bending fluctuation of a tetrahedral angle flips the topology. Measured on a 107-atom peptide, 10 ps at 500 K: 92 percent of the frames counted as topology changes at 1.5 and none at 1.3 or 1.4 -- with 1.5 the search discards the overwhelming majority of its own conformers.", "Robustness", {})
    PARAM(snapshot_clash_ratio, Double, 0.55, "A snapshot is rejected when any atom pair is closer than this fraction of the sum of their covalent radii. Independent of the topology test: a pair that is ALREADY BONDED forms no new bond when it collapses, so the topology check sees nothing, but the 1/r factors in the hydrogen-bond and three-body terms still blow up (observed as NaN in hb/atm/batm with a finite energy). 0 disables the screen.", "Robustness", {})
    PARAM(repair_snapshots, Bool, false, "Instead of discarding a snapshot whose topology changed, restrain the offending atom pairs back to a sane distance (missing bond -> its reference length, spurious contact -> outside bonding range), optimise, RELEASE the restraints and optimise freely; keep it only if the topology then matches. With GFN-FF a broken bond is often the force field losing a contact near its cutoff rather than real chemistry, and the conformation of such a snapshot can still be worth keeping. Costs two optimisations per attempt.", "Robustness", {})
    PARAM(repair_max, Int, 20, "Maximum repair attempts per cycle, lowest-energy candidates first.", "Robustness", {})
    PARAM(repair_max_bonds, Int, 2, "Only snapshots with at most this many changed bonds are repaired. More changed bonds means a different molecule, not a conformer with an artefact.", "Robustness", {})
    PARAM(repair_force, Double, 2.0, "Force constant of the distance restraints during a repair, in Eh/Angstrom^2.", "Robustness", {})
    PARAM(repair_max_iterations, Int, 300, "Maximum optimisation steps of the restrained repair stage.", "Robustness", {})
    PARAM(stage_saturation_abort, Int, 2, "Skip the remaining MIDDLE repetitions of a temperature stage after this many CONSECUTIVE dry repetitions -- dry meaning: no new best energy on either PES and not one structure added within seed_energy_window of the running global minimum. The LAST repetition of the stage always runs in full, because the once-per-stage finalisation (REFINE, seed selection) is gated on it. So with repeat 5 and saturation after repetition 2, repetitions 3-4 are skipped and repetition 5 finalises the stage. 0 disables. Measured motivation: the WEKLQ production run spent 3434 optimisations on three stages that never improved anything again, and a repeat-10 run stood still for 1500 optimisations; conversely cold stages DO deliver fine structure while their region still yields (four of the ten best conformers came from 300-400 K), which is why the criterion counts in-window additions and not only new minima.", "Efficiency", {})
    PARAM(topology_lock, Bool, true, "Lock the GFN-FF topology PERCEPTION of every child computation to the optimised input structure. After the initial optimisation, one reference .topo.json is written per PES whose method is gfnff, and every child (MD, snapshot optimisations, single points, ConfGen) adopts its hybridisation/itag and Phase-1 EEQ data instead of re-perceiving them from its own start geometry (-gfnff.topology_file; fingerprint-checked, so a structure whose BOND topology differs still derives fresh parameters and is caught by the species check as before). Why: without the lock every snapshot optimisation perceives its own topology -- a thermally distorted snapshot can cross a perception threshold, get parameters that legitimise its own distortion, and found a self-reinforcing artefact family (measured on WEKLQ: one =N-H at 179 deg -> sp via the GEODEP rule, 158 kJ/mol spurious depth, 75 percent of the pool within three temperature stages, all 10 deepest structures; the bonds-only species check passes the flip). It also removes the silently mixed energy scales of per-snapshot parameterisations (~1 kJ/mol routine spread between re-perceptions of the same conformer). Only affects force-field children (gfnff); QM methods have no perceived topology.", "Robustness", {})
    PARAM(topo_check, Bool, false, "Abort an MD run when the molecule fragments.", "Robustness", {})
    PARAM(topo_check_interval, Int, 0, "Steps between topology checks. 0 uses the MD dump frequency.", "Robustness", {})
    PARAM(epot_abort, Bool, false, "Abort an MD run when the running-mean potential climbs past epot_abort_window.", "Robustness", {})
    PARAM(epot_abort_window, Double, 250.0, "Energy window in kJ/mol above the run start energy for epot_abort.", "Robustness", {})
    PARAM(temp_abort, Bool, false, "Abort an MD run when the running-mean temperature runs away. Off by default since the strided RMSD-MTD scheme (soft residence counter) removes the unbounded bias-height growth that caused cross-cycle heating; set true to re-enable the safety net.", "Robustness", {})
    PARAM(temp_abort_factor, Double, 1.5, "Abort when the mean temperature exceeds this factor times the target. Values at or below 0 disable the check.", "Robustness", {})
    PARAM(temp_abort_delta, Double, 300.0, "Abort when the mean temperature exceeds the target plus this many Kelvin. Values at or below 0 disable the check.", "Robustness", {})

    // --- High-Level Re-Optimisation (Phase 3b) ---
    PARAM(phase3b_two_stage, Bool, true, "Split the high-level re-optimisation into a crude pre-optimisation of EVERY structure, a deduplication at that level and an accurate optimisation of the survivors. Default since Aug 2026, because deduplicating and selecting on the EXPLORATION surface was measured to be actively wrong for this class of system: on a 107-atom peptide the gfnff and gfn2 energies of one cycle correlate at r = -0.32 (rank -0.13), i.e. the structure gfnff likes best sits at gfn2 rank 10 of 14. A gfn2 single point on the gfnff geometry is better (r = +0.40, best pick at rank 2) but still weak, because the relaxation from the gfnff geometry to the gfn2 minimum spans 149-368 kJ/mol -- several times the conformer differences themselves. Only a crude optimisation on the ranking surface removes that term, which is what this mode does before anything is selected.", "Optimisation", {})
    PARAM(phase3b_preopt_preset, String, "loose", "Convergence preset of the crude stage in the two-stage mode: loose, normal, tight or verytight.", "Optimisation", {})
    PARAM(phase3b_preopt_max_iter, Int, 0, "Maximum optimisation steps in the crude stage. 0 keeps the preset value.", "Optimisation", {})
    PARAM(phase3b_preset, String, "normal", "Convergence preset of the accurate stage (and of the single-stage Phase 3b).", "Optimisation", {})
    PARAM(phase3b_filter, Bool, true, "Run the ConfScan dedup between the two stages. False chains the two optimisations directly.", "Optimisation", {})

    // --- Per-Cycle Output ---
    PARAM(cycle_output, Bool, true, "Write the ensemble and the most stable structure of every temperature cycle, on both levels of theory: <base>.cycleNN_TxxxK.<method>.xyz plus the <base>.best_per_cycle.<method>.xyz trajectory.", "Output", {})
    PARAM(ensemble_report, Int, 3, "Number of lowest conformers listed in the per-phase energy summaries. 0 prints the summary line only.", "Output", {})

    // --- Metadynamics Bias ---
    PARAM(confgen_phase, Bool, false, "Run the torsion-recombination phase (ConfGen) after the per-cycle dedup: build torsion-state combinations the cycle's minima do not contain, optimise them and add the genuinely new conformers to the cycle. Off by default -- it costs one geometry optimisation per proposal.", "Proposals", {})
    PARAM(confgen_max_proposals, Int, 20, "Maximum number of recombined structures built and optimised per cycle.", "Proposals", {})
    PARAM(confgen_templates, Int, 3, "Number of lowest-energy minima of the cycle used as geometric templates for the recombination.", "Proposals", {})
    PARAM(snapshot_dedup_rmsd, Double, 0.75, "Drop an MD snapshot before the RELAX optimisation when it is within this best-fit RMSD of a snapshot already kept. 0 = off. Measured on two repetitions of a gfn2 search: of 638 and 631 snapshots only 177 and 173 were pairwise further apart than 1.0 Angstrom, so 56-57 percent of the optimisations were spent on geometries that the later deduplication merges anyway. This is a COST SAVING, not a deduplication -- two snapshots half an Angstrom apart can still relax into different minima, and the screen cannot know that. Keep it well below -rmsd (the conformer threshold). Default 0.75 since Aug 2026: measured on a full production stage, it saved 49 percent of the optimisations without losing a single conformer within 20 kJ/mol of the minimum, and it had been carried manually in every production command line since.", "Filtering", {})
    PARAM(confgen_eval_method, String, "none", "Method with which the recombination step OPTIMISES and JUDGES its proposals, as opposed to the one it describes the ensemble with (-confgen_method, which must be a force field for the term decomposition). none = same for both, the behaviour up to Aug 2026. auto = the ranking method (-opt_method) whenever it differs. Any method name forces that one. Why it exists: measured on a 107-atom peptide, the exploration and ranking surfaces rank the same structures at r = -0.32 ... -0.46, so proposals chosen on the cheap surface are chosen against the search. The GFN2 cross-check of two new move sets is the clearest case -- their best structures sit +71.7 and +87 kJ/mol above the reference where the plain torsion move reaches +45.8. Costs one optimisation on the accurate surface per proposal; those structures then enter the refinement already relaxed.", "Proposals", {})
    PARAM(confgen_depth, Int, 3, "Maximum number of torsions changed simultaneously relative to a template. Raised from 2 to 3 in Aug 2026, when the depth stopped being a memory limit (ConfGen -proposal_candidate_cap): measured on a 398-structure ensemble with 30 proposals per depth, 3 is where the move set is best -- every built structure is a new conformer (29 of 29, against 28 of 32 at depth 2) and the best of them is 9.1 kJ/mol below the ensemble minimum, ten times the gain at depth 2. Deeper still works (4, 5, 7 all run, the ball is then sampled instead of enumerated) and gives the largest descriptor distance, but the structures get energetically worse with every step and none of them came closer to a reference conformer than the ensemble they were built from.", "Proposals", {})
    PARAM(confgen_method, String, "auto", "Energy method for the recombination step. auto = md_method when that is a force field, gfnff otherwise -- the analysis needs a per-term energy decomposition, which only a force field provides, so a search that EXPLORES with gfn2 would otherwise lose the step entirely. When the method differs from md_method the proposals are re-optimised at md_method before they enter the ensemble, unless the accurate re-optimisation runs anyway and does it.", "Proposals", {})
    PARAM(confgen_nci_moves, Bool, true, "Phase 3c also proposes hydrogen-bond moves (break one bridge, form another observed elsewhere) built with distance restraints, not only torsion recombinations.", "Proposals", {})
    PARAM(confgen_consensus, Bool, false, "Phase 3c additionally assembles structures DE NOVO from the individually most favourable torsion states, walking away from the ensemble one torsion at a time. Measured: every chemically valid assembly was a new conformer, but they sit high in energy and each costs a build plus two optimisations -- hence off by default.", "Proposals", {})
    PARAM(reduce_prefilter_window, Double, 100.0, "Energy window in kJ/mol above the cycle minimum that a RELAXED structure must fall into before it enters the RMSD deduplication. The deduplication is quadratic (133 structures take ~170 s, 10000 would take days), while a structure 400 kJ/mol above the minimum is not a relevant conformer. 0 disables the window.", "Filtering", {})
    PARAM(reduce_prefilter_energy_tol, Double, 1.0e-6, "Two relaxed structures whose energies agree to within this many Hartree ended in the SAME minimum and only one is passed to the deduplication. Snapshots from a biased trajectory collapse onto few minima, so this removes the bulk at O(N log N) instead of O(N^2). A tolerance this tight cannot merge two genuinely different conformers by accident; the RMSD pass still decides everything else. 0 disables it.", "Filtering", {})
    PARAM(verify_best, Double, 1.0, "Recompute the energy of the cycle's lowest structure with a FRESH calculator and correct it if it deviates by more than this many kJ/mol; 0 disables the check. The lowest structure sets the energy window and the running best, so a single wrong energy discards the whole cycle -- measured on a 107-atom peptide, a recombination proposal carried -18.968697 Eh while its geometry recomputes to -18.577015 Eh (1028 kJ/mol), passed the topology check, became the \"new best\" and left the cycle with 0 structures in the ensemble. The check costs one single point per cycle and is repeated up to five times, so several corrupted energies in a row are caught too.", "Filtering", {})
    PARAM(mtd_alpha_schedule, String, "", "Per-repetition schedule for the hill width, as a comma-separated list (e.g. \"0.325,0.5,0.7\"). Entry i is used in repetition i of every temperature stage; if the list is shorter than -repeat, the last entry is kept. Empty = the single value from -rmsd_mtd_alpha, unchanged behaviour. Motivated by a 43-run parameter sweep on a 107-atom peptide: the three repetitions have different jobs -- the first starts from ONE structure and should produce seeds cheaply (alpha 0.325 / k 0.1 gave 103 conformers per 1000 s), the last runs on a good seed pool and should buy depth (alpha 0.7 / k 0.05 reached -18.7869 Eh but cost 8072 s). A single value cannot do both.", "MD", {})
    PARAM(mtd_k_schedule, String, "", "Per-repetition schedule for the hill height, same syntax and rules as -mtd_alpha_schedule. Measured: k matters far less than alpha -- 0.05 and 0.1 are equivalent within noise, 0.01 loses depth, 0.2 doubles the runtime without gain.", "MD", {})
    PARAM(bias_rejected, Bool, false, "Deposit structures that FAILED the topology check as persistent hills in the bias pool instead of discarding them. They are chemically the wrong species (a tautomer, a reaction product), so they must not enter the ensemble -- but the metadynamics bias is purely geometric and does not care: a hill there tells the next trajectory that this region has been visited. Measured on a 107-atom peptide, one 600 K cycle: the 32 rejected zwitterions sit 2.24-3.54 A (median 2.82) from the nearest valid conformer, while valid conformers sit 2.20-4.80 A from each other -- so a hill on a reject is no closer to a good structure than the hills the search already places on its own minima. Note what this does NOT do: the proton transfer happens during the OPTIMISATION of a snapshot, so the trajectory never visits the rejected geometry itself; use -hold_polar_h to stop the transfer at its source.", "Filtering", {})
    PARAM(hold_polar_h, Bool, false, "Hold every polar X-H bond (X = N, O, F, S) of the reference structure at its reference length during EVERY optimisation of the search (step 2 RELAX and step 5 REFINE), with a harmonic distance restraint. Prevents a proton transfer from turning a conformer into a TAUTOMER, which the topology filter then discards after the optimisation has been paid for. Measured on a 107-atom peptide, one 600 K cycle: 32 of 96 structures were rejected for the identical transfer (H107 from O57 to N25), and those zwitterions were the DEEPEST structures of the cycle -- all ten lowest and 29 of the 32 lowest, the best of them 72.5 kJ/mol below the best neutral conformer and 17 kJ/mol below the (neutral) GOAT reference. The transfer happens during the optimisation, not in the MD snapshot, so the snapshot topology gate cannot see it. A bond sitting at its equilibrium feels no force from the restraint, so the reported energies are essentially unperturbed.", "Filtering", { "refine_hold_polar_h" })
    PARAM(hold_polar_h_force, Double, 5.0, "Force constant of that restraint in Eh/Angstrom^2.", "Filtering", { "refine_hold_polar_h_force" })
    PARAM(bias_reset, String, "never", "Opt-in: drop accumulated metadynamics hills at the start of every temperature stage. never = keep everything (default, unchanged behaviour); deposits = drop the MD deposits but KEEP the fed-back optimised minima, so the search still avoids basins it already knows while losing the accumulated hill mass; all = empty the pool completely. Motivation, measured on a 107-atom peptide: over five repetitions the pool grew to 15505 hills, and from repetition 2 on the accumulated bias tore one bond in more than 90 percent of all snapshots, so the topology gate rejected nearly everything the dynamics produced.", "Bias", {})
    PARAM(max_bias_export, Int, 0, "Maximum number of new bias structures handed to the optimisation per cycle. 0 = no limit (every new snapshot is offered). A positive value makes the export SUBSAMPLE with a stride, which silently discards snapshots: measured on a 107-atom peptide with a grown bias pool, a limit of 1000 meant a stride of 16, so 15332 of 16291 new deposits were never even checked. The topology gate and the optimisation cost money, but discarding before them costs structures.", "Bias", {})
    PARAM(rmsd_mtd_max_height, Int, -1, "Cap on the per-structure hill counter in the bias force (W = k*counter). -1 (default since Aug 2026) = ADAPTIVE: the cap is chosen per repetition so the tallest possible hill stays below 2x the retained conformer window, counter_cap = max(3, ceil(2*energy_window/k)) -- tall enough to flatten the whole retained energy landscape twice over (exploration), far below the measured catastrophe (hills of 2-98 Eh against a 0.03 Eh conformer spread tore every later-cycle MD apart within 10-20 fs; the fragment gate kept the pool clean but nothing bounded the height). 0 = explicitly unbounded (the pre-Aug-2026 behaviour). >0 = manual counter cap.", "Bias", {})
    PARAM(rmsd_mtd_freeze_inherited, Bool, false, "Freeze inherited bias heights each run so only new deposits grow. Off by default under the strided scheme (the soft counter bounds cross-run growth without freezing, which also keeps exploring inherited basins); set true for the legacy cross-run heating bound.", "Bias", {})
    PARAM(rmsd_mtd_max_gaussians, Int, -1, "Cap the shared bias pool to at most this many structures (dropping the lowest-counter non-persistent snapshots between cycles; optimised minima are always kept). -1 = unbounded. Bounds the per-step bias cost as the search accumulates structures.", "Bias", {"max_rmsd_N"})
    PARAM(rmsd_mtd_screen, Bool, true, "Skip bias hills whose Gaussian contribution is provably negligible (rotation-invariant RMSD lower bound + cutoff). Physics-preserving speedup for large pools; false = evaluate every hill.", "Bias", {})
    PARAM(rmsd_mtd_cutoff_tol, Double, 1.0e-8, "Gaussian tolerance for rmsd_mtd_screen: a hill is skipped when its lower-bound Gaussian falls below this. Smaller = more conservative.", "Bias", {})
    PARAM(rmsd_mtd_screen_margin, Double, 0.0, "Extra safety radius (RMSD length units) added to the rmsd_mtd_screen cutoff. 0 relies on the rigorous lower bound.", "Bias", {})
    PARAM(opt_feedback_bias, Bool, true, "Deposit optimised minima back into the shared bias pool.", "Bias", {})
    PARAM(opt_feedback_height, Int, 5, "Hill counter assigned to fed-back optimised minima.", "Bias", {})
    PARAM(surface_feedback, Bool, false, "Measure where the exploration surface and the ranking surface disagree, and bias the exploration away from that direction. Only active in a dual-method run. After each cycle every ranked structure also gets a single point on the EXPLORATION surface, so both energies exist for the same geometry; the difference of the two relative energies says which structures the cheap method over-rewards. That difference is then correlated with four generic shape descriptors (hydrogen bonds, close contacts, radius of gyration, buriedness of the polar atoms) and the strongest one -- if any passes -surface_feedback_min_r -- decides which structures get a taller bias hill. The coordinate is deliberately NOT hard-coded: 'too many hydrogen bonds' is true for the 107-atom peptide where this was measured (GFN-FF put 10-11-bridge folds 230 kJ/mol below everything, GFN2 put them 52-95 kJ/mol above its own best) and false for a polyol host. When no descriptor correlates, nothing happens.", "Bias", {})
    PARAM(surface_feedback_min_structures, Int, 20, "Minimum number of structures with both energies before the correlation is trusted. Below this the feedback stays inactive -- in practice it starts working from the second cycle.", "Bias", {})
    PARAM(surface_feedback_min_r, Double, 0.4, "Minimum |r| between the surface disagreement and a shape descriptor before that descriptor is used. Higher means fewer, better-founded interventions.", "Bias", {})
    PARAM(surface_feedback_strength, Double, 3.0, "Factor on the hill counter of a structure the exploration surface over-rewards (W = k*counter*strength). Bounded on purpose: an unbounded penalty repeats the very failure that made this necessary -- see the hill-height measurements in docs/CONFSEARCH_ROADMAP.md.", "Bias", {})
    PARAM(hbond_excess_max, Int, -1, "MANUAL alternative to -surface_feedback for a system whose behaviour is already known: a structure with more than (reference + this) hydrogen bonds is not used as a seed for the next cycle. -1 = off, and there is deliberately no sensible default -- the useful value depends on the molecule (the reference conformer of the measured peptide has 7 bridges, a polyol host has far more).", "Filtering", {})
    PARAM(hbond_excess_reject, Bool, false, "With -hbond_excess_max: also keep such structures OUT of the ensemble instead of only barring them from seeding. Off by default because they are real minima -- measured, they optimise to 52-95 kJ/mol above the best and the accurate method removes their excess bridges by itself.", "Filtering", {})
    PARAM(opt_feedback_prune_snapshots, Bool, false, "Remove raw MD snapshots after feeding back optimised minima.", "Bias", {})
    PARAM(mtd_permutation, Bool, true, "Feed the symmetry reorder rules found by ConfScan into the RMSD-MTD bias.", "Bias", {})
    PARAM(bias_calibration, String, "off", "Adaptive MTD hill width: off, couple or cluster. Experimental.", "Bias", {})
    PARAM(bias_couple_factor, Double, 1.0, "couple mode: place the hill half-max at this factor times rmsd.", "Bias", {})
    PARAM(bias_scale_mode, String, "global", "MTD RMSD scaling: global or weighted. Experimental.", "Bias", {})
    PARAM(bias_energy_tol, Double, 4.0, "Energy tolerance in kJ/mol when assigning optimised structures to the same minimum.", "Bias", {})

    // --- Wall Potentials (canonical SimpleMD semantics: wall_type is the GEOMETRY) ---
    PARAM(wall_type, String, "none", "Wall geometry: none, spheric or rect.", "Walls", {"wall"})
    PARAM(wall_potential, String, "harmonic", "Wall potential function: logfermi or harmonic.", "Walls", {})
    PARAM(wall_radius, Double, 0.0, "Radius of the spherical wall in Angstrom. Auto-sized when 0.", "Walls", {"wall_spheric_radius"})
    PARAM(wall_temp, Double, 298.15, "Wall strength expressed as a temperature in Kelvin.", "Walls", {})
    PARAM(wall_beta, Double, 6.0, "Steepness parameter of the wall potential.", "Walls", {})
    PARAM(wall_x_min, Double, 0.0, "Lower x boundary of the rectangular wall in Angstrom.", "Walls", {})
    PARAM(wall_x_max, Double, 0.0, "Upper x boundary of the rectangular wall in Angstrom.", "Walls", {})
    PARAM(wall_y_min, Double, 0.0, "Lower y boundary of the rectangular wall in Angstrom.", "Walls", {})
    PARAM(wall_y_max, Double, 0.0, "Upper y boundary of the rectangular wall in Angstrom.", "Walls", {})
    PARAM(wall_z_min, Double, 0.0, "Lower z boundary of the rectangular wall in Angstrom.", "Walls", {})
    PARAM(wall_z_max, Double, 0.0, "Upper z boundary of the rectangular wall in Angstrom.", "Walls", {})

    // --- Restart ---
    PARAM(restart, Bool, false, "Write a self-contained checkpoint after every sub-phase and resume from it on the next invocation.", "Restart", {})

    END_PARAMETER_DEFINITION
    // ^^^^^^^^^^^^ PARAMETER DEFINITION BLOCK ^^^^^^^^^^^^
};
