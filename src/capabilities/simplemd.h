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

#ifdef USE_Plumed
#include "plumed2/src/wrapper/Plumed.h"
#endif

#include "src/capabilities/rmsd.h"
#include "src/capabilities/rmsdtraj.h"

#include "src/core/energycalculator.h"
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
};

class BiasThread : public CxxThread {
public:
    BiasThread(const Molecule& reference, const json& rmsdconfig, bool nocolvarfile, bool nohillsfile);
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
            colvarfile.open("COLVAR_" + std::to_string(index));
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
    inline double BiasEnergy() const { return m_current_bias; }
    inline void setk(double k) { m_k = k; }
    inline void setalpha(double alpha) { m_alpha = alpha; }
    inline void setDT(double DT) { m_DT = DT; }
    inline void setdT(double dT) { m_dT = dT; }

    inline void setEnergyConv(double rmsd_econv) { m_rmsd_econv = rmsd_econv; }
    inline void setWTMTD(bool wtmtd) { m_wtmtd = wtmtd; }
    inline void setRamping(bool use_ramping) { m_use_ramping = use_ramping; }  // Claude Generated - Enable/disable ramping
    inline void setRampFactor(double ramp_factor) { m_ramp_factor = ramp_factor; }  // Claude Generated - Control ramping speed
    inline void incrementMetatime(double dt = 1.0) { m_metatime += dt; }  // Claude Generated - Increment ramping timer
    inline void resetMetatime() { m_metatime = 0.0; }  // Claude Generated - Reset timer when new structure added
    inline void removeOldestStructure() {  // Claude Generated - Rolling buffer support
        if (!m_biased_structures.empty()) {
            m_biased_structures.erase(m_biased_structures.begin());
        }
    }
    inline int Counter() const { return m_counter; }
    std::vector<BiasStructure> getBiasStructure() const { return m_biased_structures; }
    std::vector<json> getBias() const;

private:
    std::vector<BiasStructure> m_biased_structures;
    RMSDDriver m_driver;
    json m_config, m_constrained;
    Molecule m_reference, m_target;
    Geometry m_gradient;
    double m_k, m_alpha, m_DT, m_currentStep, m_rmsd_reference, m_current_bias, m_rmsd_econv, m_dT = 1;
    double m_metatime = 0.0, m_ramp_factor = 0.1;  // Claude Generated - Ramping parameters
    int m_counter = 0, m_atoms = 0;
    bool m_wtmtd = false, m_nocolvarfile = false, m_nohillsfile = false;
    bool m_use_ramping = true;  // Claude Generated - Enable ramping by default
};

// Claude Generated 2025: CurcumaMDJson removed - replaced by ParameterRegistry + ConfigManager
// Legacy 80-line static JSON object removed - all defaults now managed through Parameter Registry System

// Claude Generated 2025: Type-safe thermostat selection
enum class ThermostatType {
    Berendsen,
    Anderson,
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
    void printHelp() const;

private:
    std::function<void(void)> ThermostatFunction;
    void PrintStatus() const;

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
    void Anderson();
    void NoseHover();

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
#endif

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

    std::vector<Geometry> m_bias_structures;
    std::vector<BiasStructure> m_biased_structures;
    std::vector<BiasThread*> m_bias_threads;
    json m_bias_json;
    CxxThreadPool* m_bias_pool;
    int m_unix_started = 0, m_prev_index = 0, m_max_rescue = 10, m_current_rescue = 0, m_currentTime = 0, m_max_top_diff = 15, m_step = 0;
    int m_writerestart = -1;
    int m_respa = 1;
    int m_rattle_dynamic_tol_iter = 100;
    double m_pos_conv = 0, m_scale_velo = 1.0, m_coupling = 10;
    double m_impuls = 0, m_impuls_scaling = 0.75, m_dt2 = 0;
    double m_rattle_tol_12 = 0.1, m_rattle_tol_13 = 0.5;
    double m_wall_spheric_radius = 6, m_wall_temp = 298.15, m_wall_beta = 6;
    double m_wall_x_min = 0, m_wall_x_max = 0, m_wall_y_min = 0, m_wall_y_max = 0, m_wall_z_min = 0, m_wall_z_max = 0;
    double m_wall_potential = 0, m_average_wall_potential = 0;
    double m_virial_correction = 0, m_average_virial_correction = 0;
    double m_deltaT = 0;
    double m_k_rmsd = 0.1;
    double m_alpha_rmsd = 10;
    double m_bias_energy = 1e8;
    double m_rmsd_rmsd = 1;
    double m_rmsd_econv = 1e8;
    double m_rmsd_DT = 1000000;
    double m_rmsd_mtd_ramp_factor = 0.1;  // Claude Generated - Ramping speed parameter
    bool m_rmsd_mtd_ramping = true;  // Claude Generated - Enable ramping
    double m_rattle_max = 10;
    double m_rattle_min = 1e-4;
    int m_max_rmsd_N = -1;
    int m_mtd_steps = 10;
    int m_rattle = 0;
    int m_colvar_incr = 0;
    int m_threads = 0;
    int m_bias_structure_count = 0;
    int m_rmsd_fragment_count = 0;
    int m_wall_type = 0;
    int m_rattle_counter = 0;
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

    double m_anderson = 0.01;

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
    PARAM(time_step, Double, 1.0, "Integration time step in femtoseconds.", "Basic", {"dt"})
    PARAM(max_time, Double, 1000.0, "Maximum simulation time in femtoseconds.", "Basic", {"MaxTime"})
    PARAM(charge, Int, 0, "Total charge of the system.", "Basic", {})
    PARAM(spin, Int, 0, "Total spin multiplicity of the system (0=singlet, 1=doublet).", "Basic", {"Spin"})
    PARAM(seed, Int, -1, "Random seed (-1: time, 0: auto).", "Basic", {})
    PARAM(threads, Int, 1, "Number of parallel threads.", "Basic", {})

    // --- Thermostat ---
    PARAM(thermostat, String, "csvr", "Thermostat type: berendsen|anderson|nosehover|csvr|none.", "Thermostat", {})
    PARAM(coupling, Double, 10.0, "Thermostat coupling time in fs.", "Thermostat", {})
    PARAM(anderson_probability, Double, 0.001, "Anderson thermostat collision probability.", "Thermostat", {"anderson"})
    PARAM(chain_length, Int, 3, "Chain length for Nosé-Hoover thermostat.", "Thermostat", {"chainlength"})

    // --- System Control ---
    PARAM(remove_com_motion, Double, 100.0, "Remove translation/rotation every N fs.", "System", {"rm_COM"})
    PARAM(remove_com_mode, Int, 3, "Removal mode (0:none, 1:trans, 2:rot, 3:both).", "System", {"rmrottrans"})
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
    PARAM(rattle_tol_12, Double, 1e-1, "Tolerance for 1-2 constraints.", "RATTLE", {})
    PARAM(rattle_tol_13, Double, 2.0, "Tolerance for 1-3 constraints.", "RATTLE", {})
    PARAM(rattle_max_iterations, Int, 50, "Maximum RATTLE iterations.", "RATTLE", {"rattle_maxiter"})

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
    PARAM(rmsd_mtd_k, Double, 0.1, "Force constant for RMSD bias.", "RMSD-MTD", {"k_rmsd"})
    PARAM(rmsd_mtd_alpha, Double, 10.0, "Width parameter for RMSD Gaussians.", "RMSD-MTD", {"alpha_rmsd"})
    PARAM(rmsd_mtd_pace, Int, 1, "Add a new bias potential every N steps.", "RMSD-MTD", {"mtd_steps"})
    PARAM(rmsd_mtd_max_gaussians, Int, -1, "Maximum number of stored bias structures.", "RMSD-MTD", {"max_rmsd_N"})
    PARAM(rmsd_mtd_ref_file, String, "none", "File with reference structures for RMSD-MTD.", "RMSD-MTD", {"rmsd_ref_file"})
    PARAM(rmsd_mtd_atoms, String, "-1", "Atom indices to use for RMSD calculation.", "RMSD-MTD", {"rmsd_atoms"})
    PARAM(rmsd_mtd_dt, Double, 1000000.0, "RMSD-MTD bias deposition time.", "RMSD-MTD", {"rmsd_DT"})
    PARAM(rmsd_mtd_ramping, Bool, true, "Enable ramping function for smooth bias buildup (XTB-inspired).", "RMSD-MTD", {})  // Claude Generated
    PARAM(rmsd_mtd_ramp_factor, Double, 0.1, "Ramping speed parameter (sigmoid steepness, higher = faster).", "RMSD-MTD", {})  // Claude Generated

    // --- Coarse Graining (CG) Parameters --- Claude Generated (Nov 2025)
    PARAM(cg_write_vtf, Bool, true, "Write VTF trajectory for CG systems.", "CG", {"write_vtf"})
    PARAM(cg_timestep_scaling, Bool, true, "Enable automatic timestep scaling for pure CG systems.", "CG", {"cg_dt_scaling"})
    PARAM(cg_timestep_factor, Double, 10.0, "Timestep multiplication factor for pure CG systems.", "CG", {})

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
    ~MDThread() = default;
    void setBasename(const std::string& basename) { m_basename = basename; }
    inline void setMolecule(const Molecule& molecule) { m_molecule = molecule; }
    SimpleMD* MDDriver() const { return m_mddriver; }

    virtual int execute() override
    {
        json controller;
        controller["md"] = m_controller;
        m_mddriver = new SimpleMD(controller, false);
        m_mddriver->setMolecule(m_molecule);
        m_mddriver->overrideBasename(m_basename + ".t" + std::to_string(getThreadId()));
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
};
