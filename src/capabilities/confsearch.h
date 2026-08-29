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

    virtual int execute() override
    {
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
        return m_result.success ? 0 : 1;
    }

    const Optimization::OptimizationResult& result() const { return m_result; }

private:
    Molecule m_molecule;
    json m_parameter;
    Optimization::OptimizationResult m_result;
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

    std::string PerformOptimisation(const std::string& filename, const nlohmann::json& parameter);

    std::string PerformFilter(const std::string& filename, const nlohmann::json& parameter);

    /* Claude Generated (Jul 2026): single source of truth for every child computation.
     * Every MD run, every optimisation batch and every nested ConfScan draws its config from
     * here, so the system identity (charge/spin), the runtime settings (threads/gpu/verbosity)
     * and the method sub-scopes (gfnff/xtb/tblite/...) can never silently differ between phases.
     * Replaces the ad-hoc {method, threads, gpu} JSONs that dropped charge, spin and solvation. */
    nlohmann::json ChildConfig(const std::string& method, int threads) const;

    /* Claude Generated (Jul 2026): config for a nested ConfScan filter pass. Built from the
     * confscan registry defaults plus explicit, intentional overrides -- NOT from ConfSearch's
     * own defaults, which used to leak "method":"gfnff" into ConfScan's RMSD-alignment parameter
     * and pinned the dedup threshold to a ConfSearch default instead of the user's -rmsd. */
    nlohmann::json FilterConfig(const std::string& energy_method, int threads) const;

    /* Claude Generated (Jun 2026): experimental adaptive bias calibration (Phase C). Clusters the
     * cycle's optimised structures onto the accepted distinct minima (best-fit RMSD + cached
     * permutations + energy tolerance), then derives the basin capture radius and per-atom RMSF.
     * In "cluster" mode it sets md["rmsd_mtd_alpha"] from the basin radius; in "weighted" mode it
     * sets flexibility weights (1/sigma^2) on the shared pool for the RMSF-weighted MTD bias. */
    void CalibrateBias(const std::string& basename, nlohmann::json& md);

    /* Permutation-aware best-fit (Kabsch) RMSD in Angstrom: min over identity + cached symmetry
     * permutations applied to the target. Both geometries must share the canonical atom order. */
    double PermRMSD(const Geometry& reference, const Geometry& target) const;

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
        std::vector<int> elements; // atomic numbers (shared by all frames)
        std::vector<BiasStructure> bias; // full bias pool (geometry + metadata)
        std::vector<Molecule> seeds; // m_in_stack for the resumed cycle
        std::vector<Molecule> cumulative; // accepted conformers from completed cycles
        std::vector<Molecule> accepted_md; // gfnff-filtered set (post_filter/post_refine)
        std::vector<Molecule> accepted_opt; // opt_method-reoptimised set (post_refine)
        Molecule topo_ref; // reference structure -> m_topo_matrix
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
    std::string m_bias_calibration = "off"; // adaptive MTD width mode: off | couple | cluster
    std::string m_bias_scale_mode = "global"; // global | weighted (RMSF-weighted RMSD)
    double m_bias_couple_factor = 1.0, m_bias_energy_tol = 4.0;
    double m_global_min = std::numeric_limits<double>::infinity(); // running lowest energy across all cycles
    std::vector<std::vector<int>> m_permutation_cache; // Claude Generated (Jun 2026): symmetry reorder rules from ConfScan, fed into MTD
    Matrix m_topo_matrix;
    SharedBiasPool* m_bias_pool = nullptr;  // Claude Generated (Apr 2026): shared bias pool for parallel ConfSearch
    // Claude Generated (Jun 2026): restart/checkpoint state.
    bool m_do_restart = false;                  // -restart: enable checkpoint write + resume
    RestartState m_restart;                     // populated by loadCheckpoint() on resume
    double m_best_energy = std::numeric_limits<double>::infinity();    // exploration best (md_method), persisted
    double m_initial_energy = std::numeric_limits<double>::infinity(); // cycle-1 exploration reference, persisted
    std::string m_cumulative_file;              // path to the cumulative output (set in start())
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
    PARAM(repeat, Int, 4, "Independent MD runs started from every seed structure per cycle.", "Schedule", {})
    PARAM(time, Double, 2000.0, "Length of each MD run in femtoseconds.", "Schedule", {"max_time", "MaxTime"})

    // --- Filtering ---
    PARAM(rmsd, Double, 1.25, "RMSD threshold in Angstrom used to deduplicate conformers and to size the MTD hills.", "Filtering", {})
    PARAM(energy_window, Double, 100.0, "Energy window in kJ/mol above the running minimum for keeping conformers.", "Filtering", {})
    PARAM(seed_rank, Int, 1, "Maximum number of lowest-energy seeds carried into the next cycle. 0 keeps every structure inside seed_energy_window.", "Filtering", {})
    PARAM(seed_energy_window, Double, 50.0, "Energy window in kJ/mol versus the running global minimum for selecting next-cycle MD seeds.", "Filtering", {})
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
    PARAM(topo_check, Bool, false, "Abort an MD run when the molecule fragments.", "Robustness", {})
    PARAM(topo_check_interval, Int, 0, "Steps between topology checks. 0 uses the MD dump frequency.", "Robustness", {})
    PARAM(epot_abort, Bool, false, "Abort an MD run when the running-mean potential climbs past epot_abort_window.", "Robustness", {})
    PARAM(epot_abort_window, Double, 250.0, "Energy window in kJ/mol above the run start energy for epot_abort.", "Robustness", {})
    PARAM(temp_abort, Bool, false, "Abort an MD run when the running-mean temperature runs away. Off by default since the strided RMSD-MTD scheme (soft residence counter) removes the unbounded bias-height growth that caused cross-cycle heating; set true to re-enable the safety net.", "Robustness", {})
    PARAM(temp_abort_factor, Double, 1.5, "Abort when the mean temperature exceeds this factor times the target. Values at or below 0 disable the check.", "Robustness", {})
    PARAM(temp_abort_delta, Double, 300.0, "Abort when the mean temperature exceeds the target plus this many Kelvin. Values at or below 0 disable the check.", "Robustness", {})

    // --- Metadynamics Bias ---
    PARAM(max_bias_export, Int, 1000, "Maximum number of bias structures written out per cycle.", "Bias", {})
    PARAM(rmsd_mtd_max_height, Int, 0, "Cap on the per-structure hill counter in the bias force. 0 is unbounded.", "Bias", {})
    PARAM(rmsd_mtd_freeze_inherited, Bool, false, "Freeze inherited bias heights each run so only new deposits grow. Off by default under the strided scheme (the soft counter bounds cross-run growth without freezing, which also keeps exploring inherited basins); set true for the legacy cross-run heating bound.", "Bias", {})
    PARAM(rmsd_mtd_max_gaussians, Int, -1, "Cap the shared bias pool to at most this many structures (dropping the lowest-counter non-persistent snapshots between cycles; optimised minima are always kept). -1 = unbounded. Bounds the per-step bias cost as the search accumulates structures.", "Bias", {"max_rmsd_N"})
    PARAM(rmsd_mtd_screen, Bool, true, "Skip bias hills whose Gaussian contribution is provably negligible (rotation-invariant RMSD lower bound + cutoff). Physics-preserving speedup for large pools; false = evaluate every hill.", "Bias", {})
    PARAM(rmsd_mtd_cutoff_tol, Double, 1.0e-8, "Gaussian tolerance for rmsd_mtd_screen: a hill is skipped when its lower-bound Gaussian falls below this. Smaller = more conservative.", "Bias", {})
    PARAM(rmsd_mtd_screen_margin, Double, 0.0, "Extra safety radius (RMSD length units) added to the rmsd_mtd_screen cutoff. 0 relies on the rigorous lower bound.", "Bias", {})
    PARAM(opt_feedback_bias, Bool, true, "Deposit optimised minima back into the shared bias pool.", "Bias", {})
    PARAM(opt_feedback_height, Int, 5, "Hill counter assigned to fed-back optimised minima.", "Bias", {})
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
