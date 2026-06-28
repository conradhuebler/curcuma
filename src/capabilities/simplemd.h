/*
 * <Simple MD Module for Curcuma. >
 * Copyright (C) 2023 - 2024 Conrad Hübler <Conrad.Huebler@gmx.net>
 *               2024 Gerd Gehrisch
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
#include <ctime>
#include <functional>
#include <random>
#include <ratio>
#include <unordered_map>

#ifdef USE_Plumed
#include "plumed2/src/wrapper/Plumed.h"
#endif

#include "src/capabilities/md_diagnostics.h"
#include "src/capabilities/rmsd.h"
#include "src/capabilities/rmsdtraj.h"

#include "src/core/energycalculator.h"
#include "src/core/intra_parallel_context.h"
#include "src/core/molecule.h"

#include "external/CxxThreadPool/include/CxxThreadPool.hpp"

#include "curcumamethod.h"
#include "src/core/parameter_macros.h"  // Claude Generated - For PARAM macro definitions
#include "src/core/parameter_registry.h"  // Claude Generated - For ParameterRegistry in getSimpleMDJson()
#include "src/core/config_manager.h"    // Claude Generated - Modern parameter access layer

struct BiasStructure {
    Geometry geometry;
    double time = 0;
    double rmsd_reference = 0;
    double energy = 0;
    double factor = 1;
    int index = 0;
    int counter = 0;
    double temperature = 0;  // Claude Generated (Apr 2026): deposition temperature for cross-T propagation
    bool persistent = false; // Claude Generated (Jun 2026): fed-back optimised minimum; exempt from counter pruning
};

class SharedBiasPool;  // Claude Generated (Apr 2026): forward declaration

class BiasThread : public CxxThread {
public:
    BiasThread(const Molecule& reference, const json& rmsdconfig, bool nocolvarfile, bool nohillsfile,
               const std::string& colvar_base = "COLVAR");
    ~BiasThread();

    virtual int execute() override;
    inline void addGeometry(const Geometry& geometry, double rmsd, double time, int index)
    {
        BiasStructure str;
        str.geometry = geometry;
        str.rmsd_reference = rmsd;
        str.time = time;
        str.counter = 1;
        str.index = index;
        m_biased_structures.push_back(str);
        if (m_nocolvarfile == false) {
            std::ofstream colvarfile;
            colvarfile.open(m_colvar_base + "_" + std::to_string(index));
            colvarfile << "#m_currentStep  rmsd  bias_energy   counter  factor" << std::endl;
            colvarfile.close();
        }
        /*
                std::ofstream hillsfile;
                hillsfile.open("HILLS_" + std::to_string(index));
                hillsfile << "#m_currentStep  rmsd  bias_energy   counter  factor" << std::endl;
                hillsfile.close();
         */
    }

    inline void addGeometry(const Geometry& geometry, const json& bias)
    {
        BiasStructure str;
        str.geometry = geometry;
        str.rmsd_reference = bias["rmsd_reference"];
        str.time = bias["time"];
        str.counter = bias["counter"];
        str.index = bias["index"];
        str.factor = bias["factor"];
        str.energy = bias["energy"];

        m_biased_structures.push_back(str);
    }

    inline void setCurrentGeometry(const Geometry& geometry, double currentStep)
    {
        m_reference.setGeometry(geometry);
        m_currentStep = currentStep;
    }

    inline Geometry Gradient() const { return m_gradient; }
    inline double RMSDReference() const { return m_rmsd_reference; }
    // Exploration bias V(x) = Sum_i k*counter_i*exp(-alpha*RMSD_i^2): drives force + deposition.
    inline double BiasEnergy() const { return m_current_bias; }
    // Optional well-tempered energy (opt-in, output only — never used for force/deposition).
    inline double BiasEnergyWT() const { return m_current_bias_wt; }
    inline void setk(double k) { m_k = k; }
    inline void setalpha(double alpha) { m_alpha = alpha; }
    inline void setDT(double DT) { m_DT = DT; } // well-tempered bias temperature Delta_T (K)

    inline void setEnergyConv(double rmsd_econv) { m_rmsd_econv = rmsd_econv; }
    inline void setWTMTD(bool wtmtd) { m_wtmtd = wtmtd; }
    inline int Counter() const { return m_counter; }
    std::vector<BiasStructure> getBiasStructure() const { return m_biased_structures; }
    std::vector<json> getBias() const;

private:
    std::vector<BiasStructure> m_biased_structures;
    RMSDDriver m_driver;
    json m_config, m_constrained;
    Molecule m_reference, m_target;
    Geometry m_gradient;
    double m_k, m_alpha, m_DT, m_currentStep, m_rmsd_reference, m_current_bias, m_rmsd_econv;
    double m_current_bias_wt = 0; // well-tempered bias energy (opt-in, output only)
    int m_counter = 0, m_atoms = 0;
    bool m_wtmtd = false, m_nocolvarfile = false, m_nohillsfile = false;
    std::string m_colvar_base = "COLVAR";
};

// Claude Generated 2025: CurcumaMDJson removed - replaced by ParameterRegistry + ConfigManager
// Legacy 80-line static JSON object removed - all defaults now managed through Parameter Registry System

// Claude Generated 2025: Type-safe thermostat selection
enum class ThermostatType {
    Berendsen,
    Andersen,
    NoseHover,
    CSVR,
    None
};

// Claude Generated 2025: Type-safe wall geometry
enum class WallGeometry {
    None,
    Spheric,
    Rect
};

// Claude Generated 2025: Type-safe wall potential
enum class WallPotentialType {
    LogFermi,
    Harmonic
};

class SimpleMD : public CurcumaMethod {
public:
    SimpleMD(const json& controller = json(), bool silent = true);  // Claude Generated 2025: Default to empty JSON, ParameterRegistry provides defaults
    ~SimpleMD();

    inline void setMolecule(const Molecule& molecule)
    {
        m_molecule = molecule;
    }
    /*
    void setBaseName(const std::string& name)
    {
        m_basename = name;
    }
    */
    bool Initialise() override;

    void start() override;

    std::vector<Molecule*> UniqueMolecules() const { return m_unique_structures; }
    const Molecule& CurrentMolecule() const { return m_molecule; }  ///< Claude Generated: access final MD geometry
    bool wasStable() const { return !m_unstable; }  ///< Claude Generated: check if MD completed without instability (NaN/Inf)
    void printHelp() const;

    /* ============================================================================
     * Stepwise MD API (Claude Generated 2026 - Qurcuma interactive simulation)
     *
     * Replaces the former callback-based integration. Callers drive the MD loop
     * themselves by invoking prepareRun() → step() repeatedly → finalizeRun().
     * Between two step() calls, applyExternalForces() can inject an additive
     * gradient contribution (e.g. from a user grabbing an atom in a GUI).
     * ============================================================================ */

    /** Setup phase before the MD loop (thermostat selection, initial energy,
     *  first trajectory frame, plumed init, console headers). Must be called
     *  after Initialise(). Safe to call from start(); automatically invoked there. */
    void prepareRun();

    /** Execute one MD step (position+velocity update, thermostat, rattle, output
     *  for the current step counter). Returns false when the simulation should
     *  terminate (m_unstable, maxtime reached, CheckStop()). */
    bool step();

    /** Final printout, final trajectory frame, plumed/metadynamics finalize,
     *  curcuma_final.json. Automatically invoked at end of start(). */
    void finalizeRun();

    /** Add an external per-atom force contribution that will be applied in the
     *  next step() call (added to m_eigen_gradient before the integrator). The
     *  contribution is cleared after application. Matrix shape: (natoms, 3) in
     *  Hartree/Bohr (same units as m_eigen_gradient). */
    void applyExternalForces(const Geometry& forces);

    // Getters for stepwise GUI feedback
    const Geometry& positions() const { return m_eigen_geometry; }
    const Geometry& velocities() const { return m_eigen_velocities; }
    const Geometry& gradient() const { return m_eigen_gradient; }
    double potentialEnergy() const { return m_Epot; }
    double kineticEnergy() const { return m_Ekin; }
    double currentTemperature() const { return m_T; }
    int stepCount() const { return m_step; }
    double currentTime() const { return m_currentStep; }
    const Molecule& currentMolecule() const { return m_molecule; }

    /** Claude Generated (2026): set the thermostat target temperature live (Kelvin).
     *  Safe to call between step() calls from the driving thread. Setting it marks the
     *  run as manually overridden, so any active temperature ramp stops touching m_T0
     *  for the rest of the run (manual control wins). */
    void setTargetTemperature(double T)
    {
        m_T0 = T;
        m_global_ramp.overridden = true;
    }
    /** Claude Generated (2026): current thermostat target temperature (Kelvin). Reflects
     *  the live setpoint, including the value driven by an active temperature ramp. */
    double targetTemperature() const { return m_T0; }

    /** Claude Generated (2026): live-set the wall potential energy/temperature scale (K).
     *  Applied before the next MD step without restarting. */
    void setWallTemp(double T) { m_wall_temp = T; }
    /** Claude Generated (2026): live-set the wall potential steepness parameter (β). */
    void setWallBeta(double beta) { m_wall_beta = beta; }

    // Claude Generated (Apr 2026): shared bias pool for parallel ConfSearch
    void setSharedBiasPool(SharedBiasPool* pool) { m_shared_pool = pool; }

private:
    std::function<void(void)> ThermostatFunction;
    void PrintStatus() const;

    // Claude Generated 2026: Build path inside snapshots subdirectory
    std::string snapshotPath(const std::string& filename) const {
        if (m_snapshots_dir.empty())
            return outputPath(filename);
        return m_snapshots_dir + "/" + filename;
    }

    /* Lets have this for all modules */
    virtual nlohmann::json WriteRestartInformation() override;

    /* Lets have this for all modules */
    virtual bool LoadRestartInformation() override;

    bool LoadRestartInformation(const json& state);

    virtual StringList MethodName() const override
    {
        return { "MD" };
    }

    /* Lets have all methods read the input/control file */
    virtual void ReadControlFile() override
    {
    }

    /* Read Controller has to be implemented for all */
    virtual void LoadControlJson() override;

    void InitVelocities(double scaling = 1.0);

    double FastEnergy();
    double CleanEnergy();

    void PrintMatrix(const double* matrix) const;

    bool WriteGeometry();
    void applyPeriodicBoundaryConditions();  // Claude Generated (Oct 2025): PBC wrapping
    void Verlet();
    void Rattle();
    void ApplyRMSDMTD();

    void Rattle_Verlet_First(double* coord, double* grad);
    void Rattle_Constrain_First(double* coord, double* grad);
    void Rattle_Verlet_Second(double* coord, double* grad);
    void Rattle_Constrain_Second(double* coord, double* grad);

    void AdjustRattleTolerance();

    void RemoveRotation();
    void RemoveRotations();

    void EKin();
    void AverageQuantities();

    void Berendson();
    void CSVR();
    void None();
    void Andersen();
    void NoseHover();

    // Claude Generated (2026): runtime temperature control / multi-stage ramp / regions.
    struct RampSegment {
        double target;                     ///< segment target temperature [K]
        enum Mode { Steps, Reach } mode;   ///< Steps: ramp over N steps; Reach: hold until <T> reached
        double value;                      ///< Steps: number of steps; Reach: tolerance [K]
    };
    /// One temperature schedule (global or per region) plus its running position.
    struct RampState {
        bool enabled = false;              ///< schedule active
        bool overridden = false;           ///< a live setTargetTemperature() cancelled it (global only)
        std::vector<RampSegment> schedule;
        int idx = 0;                       ///< active segment
        int seg_start_step = 0;            ///< m_step when the active segment began
        double seg_start_T = 298.15;       ///< setpoint at the start of the active segment
    };
    /// A subset of atoms thermostatted to their own (optionally ramped) target temperature.
    struct ThermalRegion {
        std::string atoms_spec;            ///< selection string (FragString2Indicies grammar)
        std::vector<int> atoms;            ///< resolved 0-based indices (set in prepareRun)
        double T0 = 298.15;                ///< current setpoint (driven by ramp)
        int dof = 0;                       ///< 3 * atoms.size()
        RampState ramp;
    };

    bool ParseSchedule(const std::string& spec, std::vector<RampSegment>& out, const std::string& ctx);  ///< parse 'T:mode:val;...'
    void StepRamp(RampState& rs, double& T0, double measuredT);  ///< advance one schedule by one step
    void UpdateTemperatureRamp();                                ///< drive global + region setpoints (called each step)
    void ParseThermalRegions();                                  ///< read temp_regions specs from the controller
    void ResolveThermalRegions();                                ///< resolve atom indices + default complement (needs molecule)
    double RegionTemperature(const std::vector<int>& atoms, int dof) const;  ///< instantaneous T of an atom subset
    void ApplyThermostat();                                      ///< per-region dispatch (or legacy global path)
    void ApplyThermostatRegion(const std::vector<int>& atoms, double T0, int dof, ThermostatType type);

    void InitialiseWalls();

    double ApplySphericLogFermiWalls();
    double ApplyRectLogFermiWalls();

    double ApplySphericHarmonicWalls();
    double ApplyRectHarmonicWalls();

    void InitConstrainedBonds();

    std::function<void()> Integrator;
    std::function<double()> Energy;
    std::function<double()> WallPotential;

    std::vector<std::pair<std::pair<int, int>, double>> m_bond_constrained, m_bond_13_constrained;
#ifdef USE_Plumed
    plumed m_plumedmain;
    std::vector<double> m_plumed_positions;  // Positions in Bohr for PLUMED
    double m_plumed_box[9] = {};             // Box vectors in Bohr for PLUMED (3x3 row-major)
                                             // NOTE: Currently never populated — PBC/PLUMED periodic box is NOT supported.
                                             // To enable: populate from lattice vectors when PBC is active, then call
                                             // plumed_cmd("setBox", m_plumed_box) during PLUMED initialization.
#endif

    std::string m_snapshots_dir;  // Claude Generated 2026: Snapshots subdirectory inside BMT

    int m_natoms = 0;
    int m_dump = 1;
    double m_T = 0, m_Epot = 0, m_aver_Epot = 0, m_Ekin = 0, m_aver_Ekin = 0, m_Etot = 0, m_aver_Etot = 0, m_aver_dipol_linear = 0;
    double m_rm_COM = 100;
    int m_rm_COM_step = -1;
    int m_hmass = 1;
    double m_single_step = 1;
    double m_dT = 0.5, m_currentStep = 0, m_maxtime = 1000;
    int m_spin = 0, m_charge = 0, m_print = 100;
    double m_T0 = 298.15, m_aver_Temp = 0, m_aver_rattle_Temp = 0, m_rmsd = 1.5;
    // Claude Generated (Jun 2026): initial velocity sampling temperature. Set to
    // -1 in the constructor and resolved against m_T0 in performMolecularDynamics
    // so the thermostat target (m_T0) and the MB-sampling temperature can be
    // controlled independently. -1 means "same as target temperature".
    double m_T_init = -1.0;
    double m_x0 = 0, m_y0 = 0, m_z0 = 0;
    double m_Ekin_exchange = 0.0;
//    std::vector<double> m_current_geometry, m_mass, m_velocities, m_gradient, m_rmass, m_virial, m_gradient_bias, m_scaling_vector_linear, m_scaling_vector_nonlinear, m_rt_geom_1, m_rt_geom_2, m_rt_velo;
    std::vector<double>  m_virial, m_gradient_bias, m_scaling_vector_linear, m_scaling_vector_nonlinear, m_rt_geom_1, m_rt_geom_2, m_rt_velo;

    std::vector<Position> m_curr_dipoles;
    std::vector<int> m_atomtype;
    Molecule m_molecule, m_reference, m_target, m_rmsd_mtd_molecule;
    ConfigManager m_config;  // Claude Generated - Modern type-safe parameter access
    bool m_initialised = false, m_restart = false, m_writeUnique = true, m_opt = false, m_rescue = false, m_writeXYZ = true, m_writeinit = false, m_norestart = false;
    int m_rmrottrans = 0, m_rattle_maxiter = 100;
    bool m_nocenter = false;
    bool m_COM = false;
    bool m_wall_render = false;
    EnergyCalculator* m_interface;
    RMSDTraj* m_unqiue;
    const std::vector<double> m_used_mass;
    std::vector<int> m_rmsd_indicies;
    std::vector<std::vector<int> > m_rmsd_fragments, m_start_fragments;

    Geometry m_eigen_geometry, m_eigen_geometry_old, m_eigen_gradient, m_eigen_gradient_old, m_eigen_velocities;
    Vector m_eigen_masses, m_eigen_inv_masses;

    // Stepwise-API state (Claude Generated 2026)
    Geometry m_external_forces;          // additive per-atom force contribution, cleared after use
    bool m_external_forces_pending = false;
    bool m_run_prepared = false;         // true after prepareRun(), false after finalizeRun()
    bool m_run_aborted = false;          // mirrors former local `aborted` flag in start()
    std::vector<json> m_run_states;      // rescue states carried across step() calls

    // Claude Generated (Jun 2026): ConfSearch robustness gates (opt-in, default off)
    bool m_topo_check = false;           // abort on fragmentation; also gates bias deposition
    int m_topo_check_interval = 0;       // user value (0 -> resolved to dump in prepareRun)
    int m_topo_check_every = 0;          // resolved interval actually used in step()
    int m_start_fragment_count = 1;      // connected-component count at run start (reference)
    bool m_epot_abort = false;           // abort when running-mean potential climbs too high
    double m_epot_abort_window = 250.0;  // kJ/mol above the starting energy
    double m_epot_ref = 0.0;             // bare potential energy at run start (reference)
    bool m_temp_abort = false;           // abort when running-mean temperature runs away from target
    double m_temp_abort_factor = 1.5;    // abort if <T> > factor * T0 (<= 0 disables)
    double m_temp_abort_delta = 300.0;   // abort if <T> > T0 + delta [K] (<= 0 disables)

    // Claude Generated (2026): runtime temperature control + multi-stage ramp + regions.
    // m_T0 is the global thermostat setpoint, read every step by the thermostats. The global
    // ramp (m_global_ramp) drives m_T0 across its segments; a live setTargetTemperature() sets
    // m_global_ramp.overridden, after which it leaves m_T0 alone for the rest of the run.
    // Each ThermalRegion thermostats its own atom subset to its own (optionally ramped) target;
    // atoms in no region (m_default_region_atoms) follow the global setpoint.
    bool m_temp_ramp = false;                       // global ramp enable (PARAM temp_ramp)
    RampState m_global_ramp;                         // global schedule state
    std::vector<ThermalRegion> m_thermal_regions;    // optional per-atom-subset thermostats
    std::vector<int> m_default_region_atoms;         // atoms in no region (complement; resolved in prepareRun)
    int m_default_region_dof = 0;                    // 3 * m_default_region_atoms.size()

    std::vector<Geometry> m_bias_structures;
    std::vector<BiasStructure> m_biased_structures;
    std::vector<BiasThread*> m_bias_threads;
    json m_bias_json;
    CxxThreadPool* m_bias_pool;
    SharedBiasPool* m_shared_pool = nullptr;  // Claude Generated (Apr 2026): shared bias pool for parallel ConfSearch
    RMSDDriver m_shared_pool_driver;  // Claude Generated (Apr 2026): local RMSDDriver for shared pool path
    Molecule m_shared_pool_target;  // Claude Generated (Apr 2026): target molecule for shared pool RMSD
    int m_unix_started = 0, m_prev_index = 0, m_max_rescue = 10, m_current_rescue = 0, m_currentTime = 0, m_max_top_diff = 15, m_step = 0;
    int m_writerestart = -1;
    int m_respa = 1;
    int m_rattle_dynamic_tol_iter = 100;
    double m_pos_conv = 0, m_scale_velo = 1.0, m_coupling = 10;
    double m_impuls = 0, m_impuls_scaling = 0.75, m_dt2 = 0;
    double m_rattle_tol_12 = 1e-4, m_rattle_tol_13 = 1e-3;
    double m_wall_spheric_radius = 6, m_wall_temp = 298.15, m_wall_beta = 6;
    double m_wall_x_min = 0, m_wall_x_max = 0, m_wall_y_min = 0, m_wall_y_max = 0, m_wall_z_min = 0, m_wall_z_max = 0;
    double m_wall_potential = 0, m_average_wall_potential = 0;
    double m_virial_correction = 0, m_average_virial_correction = 0;
    double m_deltaT = 0;
    double m_k_rmsd = 0.01;
    double m_alpha_rmsd = 10;
    double m_bias_energy = 1e8;
    double m_rmsd_rmsd = 1;
    double m_rmsd_econv = 1e8;
    double m_rmsd_DT = 2000; // well-tempered bias temperature Delta_T (K), used only when wtmtd=true
    double m_rattle_max = 10;
    double m_rattle_min = 1e-4;
    int m_max_rmsd_N = -1;
    int m_rmsd_mtd_max_height = 0;       // Claude Generated (Jun 2026): cap on counter used in W_i (0 = unbounded)
    bool m_freeze_inherited = false;     // Claude Generated (Jun 2026): freeze heights of structures inherited at run start
    std::unordered_map<int, int> m_frozen_height; // index -> frozen counter for inherited bias structures
    int m_mtd_steps = 10;
    int m_rattle = 0;
    int m_colvar_incr = 0;
    int m_threads = 0;
    int m_bias_structure_count = 0;
    int m_rmsd_fragment_count = 0;
    int m_wall_type = 0;
    int m_rattle_counter = 0;
    // RATTLE diagnostics: accumulated stats per print interval
    int m_rattle_iters_step1 = 0;    // iterations used in 1st step (position constraints)
    int m_rattle_iters_step2 = 0;    // iterations used in 2nd step (velocity constraints)
    int m_rattle_max_err_count = 0;  // how many steps hit max iterations
    double m_rattle_max_err_12 = 0;  // worst 1-2 constraint error this interval
    double m_rattle_max_err_13 = 0;  // worst 1-3 constraint error this interval
    int m_rattle_constrained_atoms = 0; // atoms involved in constraints
    std::vector<double> m_collected_dipole;
    Matrix m_topo_initial;
    std::vector<Molecule*> m_unique_structures;
    std::string m_method = "UFF", m_initfile = "none", m_thermostat = "csvr", m_plumed, m_rmsd_ref_file, m_rmsd_atoms = "-1", m_scaling_json = "none";
    bool m_unstable = false;
    bool m_dipole = false;
    bool m_clean_energy = false;
    bool m_mtd = false;
    bool m_eval_mtd = true;
    bool m_rmsd_mtd = false;
    bool m_wtmtd = false;
    bool m_rmsd_fix_structure = false;
    bool m_rattle_dynamic_tol = false;
    bool m_nocolvarfile = false;
    bool m_nohillsfile = false;
    bool m_rattle_12 = false;
    bool m_rattle_13 = false;

    // WP-S2 (May 2026): per-step diagnostics JSONL dump
    bool m_md_diagnostics = false;
    std::unique_ptr<MDDiagnosticsWriter> m_diag_writer;

    // WP-P1 (May 2026): per-phase wall-clock breakdown for diagnostics
    bool m_md_diagnostics_timing = false;
    double m_last_ff_ms = 0.0;     ///< wall-clock of last m_interface->CalculateEnergy()
    double m_last_hbxb_ms = 0.0;   ///< placeholder; HBXB-update lives inside Calculation() and is hard to isolate
    double m_last_integrator_ms = 0.0;  ///< wall-clock of last Integrator() call in step()

    int m_mtd_dT = -1;
    int m_seed = -1;
    int m_time_step = 0;
    int m_dof = 0;
    int m_mtd_time = 0, m_loop_time = 0;

    std::vector<std::vector<double>> m_atom_temp;
    std::vector<double> m_zeta; // Thermostatische Variablen
    std::vector<double> m_xi; // Zeitderivate von zeta
    std::vector<double> m_Q; // Trägheiten der Thermostatkette
    double m_eta; // Variable zur Speicherung der Thermostatenergie

    int m_chain_length = 3; // Länge der Thermostatkette

    double m_andersen = 0.01;

    std::vector<std::pair<double, double>> m_rattle_tol_temp;

    // Claude Generated: Wall configuration and statistics tracking
    bool m_wall_auto_configured = false;
    std::string m_wall_geometry = "none";
    std::string m_wall_potential_type = "harmonic";
    int m_wall_violation_count = 0;
    int m_wall_violation_last_reported = 0;
    double m_molecular_density = 0.0; // molecules/Å³

    // Claude Generated (Oct 2025): Coarse Graining system detection and optimization
    bool m_is_cg_system = false;           // True if molecule contains CG_ELEMENT (226)
    bool m_has_pbc = false;                // True if periodic boundary conditions active
    int m_cg_atom_count = 0;               // Number of CG atoms in system
    double m_cg_timestep_factor = 1.0;     // Timestep scaling for CG systems (10x for pure CG)

    // Claude Generated (Nov 2025): Orientational dynamics infrastructure (Phase 5 - prepared for Phase 6 ellipsoids)
    bool m_cg_enable_rotation = false;           // Future: enable orientational dynamics for ellipsoids
    std::vector<Eigen::Vector3d> m_cg_orientations;  // Future: orientation vectors (quaternions in Phase 6)
    std::vector<Eigen::Vector3d> m_cg_angular_velocities;  // Future: angular velocities
    bool m_cg_write_vtf = true;                   // Write VTF trajectory for CG systems

    // vvvvvvvvvvvv PARAMETER DEFINITION BLOCK vvvvvvvvvvvv
    // Claude Generated - Parameter Registry Integration (October 2025)
    BEGIN_PARAMETER_DEFINITION(simplemd)

    // --- Basic Simulation Parameters ---
    PARAM(method, String, "uff", "Energy calculation method (e.g., uff, gfn2).", "Basic", {})
    PARAM(temperature, Double, 298.15, "Target temperature in Kelvin.", "Basic", {"T"})
    PARAM(initial_temperature, Double, -1.0, "Initial temperature for velocity sampling (K). -1: same as 'temperature'. Use this to anneal into the target or to start cold/warm; the thermostat still drives toward 'temperature'. Ignored on restart (velocities come from the restart file).", "Basic", {"T_init", "T0", "initT"})
    PARAM(time_step, Double, 1.0, "Integration time step in femtoseconds.", "Basic", {"dt"})
    PARAM(max_time, Double, 1000.0, "Maximum simulation time in femtoseconds.", "Basic", {"MaxTime"})
    PARAM(charge, Int, 0, "Total charge of the system.", "Basic", {})
    PARAM(spin, Int, 0, "Total spin multiplicity of the system (0=singlet, 1=doublet).", "Basic", {"Spin"})
    PARAM(seed, Int, -1, "Random seed (-1: time, 0: auto).", "Basic", {})
    PARAM(threads, Int, 1, "Number of parallel threads.", "Basic", {})

    // --- Thermostat ---
    PARAM(thermostat, String, "csvr", "Thermostat type: berendsen|andersen|nosehover|csvr|none.", "Thermostat", {})
    PARAM(coupling, Double, 10.0, "Thermostat coupling time in fs.", "Thermostat", {})
    PARAM(andersen_probability, Double, 0.001, "Andersen thermostat collision probability.", "Thermostat", {})
    PARAM(chain_length, Int, 3, "Chain length for Nosé-Hoover thermostat.", "Thermostat", {"chainlength"})

    // --- System Control ---
    PARAM(remove_com_motion, Double, 100.0, "Remove translation/rotation every N fs.", "System", {"rm_COM"})
    PARAM(remove_com_mode, Int, 1, "Removal mode (0:none, 1:trans only, 2:rot only, 3:both). Rotation removal is opt-in (use 2 or 3).", "System", {"rmrottrans"})
    PARAM(no_center, Bool, false, "Disable centering of the molecule at the origin.", "System", {"nocenter"})
    PARAM(use_com, Bool, false, "Use center of mass (otherwise geometric center).", "System", {"COM"})
    PARAM(hydrogen_mass, Int, 1, "Hydrogen mass scaling factor for HMR.", "System", {"hmass"})
    PARAM(initial_velocity_scale, Double, 1.0, "Initial velocity scaling factor.", "System", {"velo"})

    // --- Output & Restart ---
    PARAM(dump_frequency, Int, 50, "Save coordinates every N steps.", "Output", {"dump"})
    PARAM(print_frequency, Int, 1000, "Print status every N steps.", "Output", {"print"})
    PARAM(write_xyz, Bool, true, "Write trajectory to XYZ file.", "Output", {"writeXYZ"})
    PARAM(write_initial_state, Bool, false, "Write initial conditions to a .init.json file.", "Output", {"writeinit"})
    PARAM(restart_file, String, "none", "Restart file to load initial state from.", "Restart", {"initfile"})
    PARAM(write_restart_frequency, Int, 1000, "Write restart file every N steps.", "Restart", {"writerestart"})
    PARAM(no_restart, Bool, false, "Disable automatic loading from restart files.", "Restart", {"norestart"})

    // --- RATTLE Constraints ---
    PARAM(rattle, Int, 0, "RATTLE constraint algorithm (0:off, 1:on, 2:H-only).", "RATTLE", {})
    PARAM(rattle_12, Bool, true, "Constrain 1-2 bond distances.", "RATTLE", {})
    PARAM(rattle_13, Bool, false, "Constrain 1-3 distances (angles).", "RATTLE", {})
    PARAM(rattle_tol_12, Double, 1e-4, "Tolerance for 1-2 constraints (squared Bohr).", "RATTLE", {})
    PARAM(rattle_tol_13, Double, 1e-3, "Tolerance for 1-3 constraints (squared Bohr).", "RATTLE", {})
    PARAM(rattle_max_iterations, Int, 100, "Maximum RATTLE iterations.", "RATTLE", {"rattle_maxiter"})

    // --- Wall Potentials ---
    PARAM(wall_type, String, "none", "Wall type: none|spheric|rect.", "Walls", {"wall"})
    PARAM(wall_potential, String, "harmonic", "Wall potential function: logfermi|harmonic.", "Walls", {})
    PARAM(wall_radius, Double, 0.0, "Radius for spherical wall (Å). Auto-sized if 0.", "Walls", {"wall_spheric_radius"})
    PARAM(wall_temp, Double, 298.15, "Wall temperature/strength in K.", "Walls", {})
    PARAM(wall_beta, Double, 6.0, "Steepness parameter for wall potential.", "Walls", {})
    PARAM(wall_x_min, Double, 0.0, "Min x-boundary for rectangular wall (Å).", "Walls", {})
    PARAM(wall_x_max, Double, 0.0, "Max x-boundary for rectangular wall (Å).", "Walls", {})
    PARAM(wall_y_min, Double, 0.0, "Min y-boundary for rectangular wall (Å).", "Walls", {})
    PARAM(wall_y_max, Double, 0.0, "Max y-boundary for rectangular wall (Å).", "Walls", {})
    PARAM(wall_z_min, Double, 0.0, "Min z-boundary for rectangular wall (Å).", "Walls", {})
    PARAM(wall_z_max, Double, 0.0, "Max z-boundary for rectangular wall (Å).", "Walls", {})

    // --- Metadynamics (PLUMED) ---
    PARAM(mtd, Bool, false, "Enable PLUMED metadynamics.", "Metadynamics", {})
    PARAM(plumed_file, String, "plumed.dat", "PLUMED input file.", "Metadynamics", {"plumed"})

    // --- RMSD-based Metadynamics (Internal) ---
    PARAM(rmsd_mtd, Bool, false, "Enable internal RMSD-based metadynamics.", "RMSD-MTD", {})
    PARAM(rmsd_mtd_k, Double, 0.01, "Hill-height constant: bias height W_i = k * counter_i (Eh). The force is the exact gradient of the bias, so k is ~100x smaller than the pre-2026 value.", "RMSD-MTD", {"k_rmsd"})
    PARAM(rmsd_mtd_alpha, Double, 10.0, "Width parameter for RMSD Gaussians.", "RMSD-MTD", {"alpha_rmsd"})
    PARAM(rmsd_mtd_pace, Int, 1, "Unused in the counter-based scheme (kept for compatibility); deposition is gated by the bias level, not a fixed pace.", "RMSD-MTD", {"mtd_steps"})
    PARAM(rmsd_mtd_max_gaussians, Int, -1, "Maximum number of stored bias structures.", "RMSD-MTD", {"max_rmsd_N"})
    PARAM(rmsd_mtd_ref_file, String, "none", "File with reference structures for RMSD-MTD.", "RMSD-MTD", {"rmsd_ref_file"})
    PARAM(rmsd_mtd_atoms, String, "-1", "Atom indices to use for RMSD calculation.", "RMSD-MTD", {"rmsd_atoms"})
    PARAM(rmsd_mtd_dt, Double, 2000.0, "Well-tempered bias temperature Delta_T (K). Only used when wtmtd=true, and only for the reported well-tempered energy -- it never affects the force or the exploration.", "RMSD-MTD", {"rmsd_DT"})
    PARAM(rmsd_mtd_max_height, Int, 0, "Cap the per-structure hill counter used in the bias force: W_i = k * min(counter_i, cap). 0 = unbounded (legacy). Stops the shared bias pool from heating the dynamics over many runs (counter_i grows on every visit).", "RMSD-MTD", {})
    PARAM(rmsd_mtd_freeze_inherited, Bool, false, "Freeze the hill heights of bias structures already present at this MD run's start; only structures deposited during this run gain height. Bounds the cumulative bias force across successive shared-pool runs (geometry sharing is preserved).", "RMSD-MTD", {})

    // --- Coarse Graining (CG) Parameters --- Claude Generated (Nov 2025)
    PARAM(cg_write_vtf, Bool, true, "Write VTF trajectory for CG systems.", "CG", {"write_vtf"})
    PARAM(cg_timestep_scaling, Bool, true, "Enable automatic timestep scaling for pure CG systems.", "CG", {"cg_dt_scaling"})
    PARAM(cg_timestep_factor, Double, 10.0, "Timestep multiplication factor for pure CG systems.", "CG", {})

    // --- WP-S2 Diagnostics (May 2026) ---
    PARAM(md_diagnostics, Bool, false, "Write per-step diagnostics to <basename>.diag.jsonl (energy decomposition, charges, CN, gradient norms, HB/XB counts). Frequency follows dump_frequency. One JSON object per line.", "Output", {})
    // --- WP-P1 Timing Instrumentation (May 2026) ---
    PARAM(md_diagnostics_timing, Bool, false, "Add a timing_ms block to each <basename>.diag.jsonl record (per-phase wall-clock: CN/EEQ/dcn/D4-weights/FF/integrator/HBXB/I-O). GPU runs add a gpu sub-block with per-kernel-category times. Requires md_diagnostics=true. ~1-2 us per hook.", "Output", {})

    // --- ConfSearch robustness gates (Jun 2026, Claude Generated) ---
    PARAM(topo_check, Bool, false, "Abort the MD run when the molecule fragments (number of connected components grows above the start value). Off by default so reactive paths are still sampled.", "ConfSearch", {})
    PARAM(topo_check_interval, Int, 0, "Steps between topology checks (0 -> use dump_frequency).", "ConfSearch", {})
    PARAM(epot_abort, Bool, false, "Abort the MD run when the running-mean potential energy climbs more than epot_abort_window above the run's starting energy.", "ConfSearch", {})
    PARAM(epot_abort_window, Double, 250.0, "Energy window (kJ/mol) above the starting energy for epot_abort. Must exceed the thermal baseline (~N_dof*kT/2) plus typical barriers.", "ConfSearch", {})
    PARAM(temp_abort, Bool, false, "Abort the MD run when the running-mean temperature runs away from the target (catches bias-driven heating). Uses temp_abort_factor and/or temp_abort_delta.", "ConfSearch", {})
    PARAM(temp_abort_factor, Double, 1.5, "Abort when <T> exceeds temp_abort_factor * target T. <= 0 disables this threshold. Only active when temp_abort=true.", "ConfSearch", {})
    PARAM(temp_abort_delta, Double, 300.0, "Abort when <T> exceeds (target T + temp_abort_delta) Kelvin. <= 0 disables this threshold. Only active when temp_abort=true.", "ConfSearch", {})

    // --- Temperature Ramp (Jun 2026, Claude Generated) ---
    PARAM(temp_ramp, Bool, false, "Enable a multi-stage temperature ramp schedule (see temp_schedule). A live GUI slider / setTargetTemperature() overrides it for the rest of the run.", "Temperature Ramp", {})
    PARAM(temp_schedule, String, "", "Ramp schedule 'target:mode:value;...'. mode=steps ramps the setpoint linearly from the previous target to <target> over <value> integration steps; mode=reach jumps the setpoint to <target> and advances once |<T>-target| < value Kelvin. Example: '500:steps:5000;500:steps:2000;300:reach:10'.", "Temperature Ramp", {})

    END_PARAMETER_DEFINITION
    // ^^^^^^^^^^^^ PARAMETER DEFINITION BLOCK ^^^^^^^^^^^^
};

// Claude Generated (October 2025): Replace static JSON with ParameterRegistry
// This ensures consistency with the parameter definitions above and enables
// automatic help generation, validation, and JSON export/import capabilities
inline json getSimpleMDJson() {
    return ParameterRegistry::getInstance().getDefaultJson("simplemd");
}

// Legacy static reference for backward compatibility - Claude Generated
// NOTE: This creates the JSON on first access, so it's safe for initialization order
#define SimpleMDJson getSimpleMDJson()

class MDThread : public CxxThread {
public:
    MDThread(const json& controller)
        : m_controller(controller)
    {
        setAutoDelete(true);
    }
    ~MDThread() { delete m_mddriver; }
    void setBasename(const std::string& basename) { m_basename = basename; }
    inline void setMolecule(const Molecule& molecule) { m_molecule = molecule; }
    SimpleMD* MDDriver() const { return m_mddriver; }

    // Claude Generated (Apr 2026): shared bias pool for parallel ConfSearch
    void setSharedBiasPool(SharedBiasPool* pool) { m_shared_pool = pool; }

    virtual int execute() override
    {
        // One MD run among many under a molecule-level pool: keep intra-molecule
        // fan-out suppressed so methods that honor the flag stay serial.
        curcuma::SuppressIntraParallel intra_guard;

        m_mddriver = new SimpleMD(m_controller, false);
        m_mddriver->setMolecule(m_molecule);
        m_mddriver->overrideBasename(m_basename + ".t" + std::to_string(getThreadId()));
        m_mddriver->setSharedBiasPool(m_shared_pool);
        m_mddriver->Initialise();
        m_mddriver->start();
        return 0;
    }

protected:
    Molecule m_molecule;
    std::string m_basename;
    std::string m_result;
    json m_controller;
    SimpleMD* m_mddriver;
    SharedBiasPool* m_shared_pool = nullptr;  // Claude Generated (Apr 2026)
};
