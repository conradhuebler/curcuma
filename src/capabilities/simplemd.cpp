/*
 * <Simple MD Module for Cucuma. >
 * Copyright (C) 2020 - 2024 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include <cmath>
#include <cstring>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <set>
#include <string>
#include <vector>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "src/capabilities/optimizer_factory.h"
#include "src/capabilities/rmsd.h"
#include "src/capabilities/rmsdtraj.h"
#include "src/capabilities/shared_bias_pool.h"  // Claude Generated (Apr 2026)

#include "src/core/elements.h"
#include "src/core/energycalculator.h"
#include "src/core/fileiterator.h"
#include "src/core/global.h"
#include "src/core/molecule.h"
#include "src/core/parameter_registry.h"  // Claude Generated 2025: For ParameterRegistry::getInstance()

#include "src/tools/geometry.h"

#include "external/CxxThreadPool/include/CxxThreadPool.hpp"

#ifdef USE_Plumed
#include "plumed2/src/wrapper/Plumed.h"
#endif
#include "simplemd.h"
#include "src/core/citation_registry.h"

// Claude Generated: Unit conversion constants for wall statistics
const double au2eV = 1.0 / eV2Eh; // Convert Hartree to eV
const double au2N = 8.2387225e-8; // Convert atomic force units (Eh/bohr) to Newton

BiasThread::BiasThread(const Molecule& reference, const json& rmsdconfig, bool nocolvarfile, bool nohillsfile,
                       const std::string& colvar_base)
    : m_reference(reference)
    , m_target(reference)
    , m_nocolvarfile(nocolvarfile)
    , m_nohillsfile(nohillsfile)
    , m_colvar_base(colvar_base)
    , m_driver(rmsdconfig, true)
{
    m_config = rmsdconfig;
    setAutoDelete(true);
    m_current_bias = 0;
    m_counter = 0;
    m_atoms = m_reference.AtomCount();
    m_gradient = Eigen::MatrixXd::Zero(m_reference.AtomCount(), 3);
}

BiasThread::~BiasThread()
{
}

int BiasThread::execute()
{
    if (m_biased_structures.empty())
        return 0;
    m_current_bias = 0;    // exploration bias V(x): drives the force and deposition
    m_current_bias_wt = 0; // optional well-tempered energy (opt-in, output only)
    m_counter = 0;
    m_driver.setReference(m_reference); // reference = current walker geometry
    m_gradient = Eigen::MatrixXd::Zero(m_reference.AtomCount(), 3);

    // Non-textbook design (intentional): instead of depositing a fresh fixed-height
    // Gaussian every 'pace' steps, we keep ONE reference per region and raise its
    // height on every visit. Hill height W_i = m_k * counter_i, so the bias potential is
    //   V(x) = Sum_i W_i * exp(-alpha * RMSD(x, x_i)^2)
    // and the force below is its EXACT negative gradient. Claude Generated (Jun 2026).
    std::vector<int> visited;

    for (int i = 0; i < m_biased_structures.size(); ++i) {
        m_target.setGeometry(m_biased_structures[i].geometry);
        m_driver.setTarget(m_target);
        double rmsd = m_driver.BestFitRMSD();
        double expr = exp(-rmsd * rmsd * m_alpha);
        double height = m_k * m_biased_structures[i].counter; // W_i = k * counter_i
        double bias_energy = height * expr;

        m_current_bias += bias_energy;
        if (m_wtmtd) // separate well-tempered weight 'factor', output only
            m_current_bias_wt += m_k * m_biased_structures[i].factor * expr;

        if (i == 0)
            m_rmsd_reference = rmsd;

        // Force = -dV/dx for this hill. RMSDDriver::Gradient() returns dRMSD/dx of the
        // reference (= the walker), so dV/dx = W_i * exp * (-2*alpha*rmsd) * Gradient().
        double dEdR = -2.0 * m_alpha * height * rmsd * expr;
        m_gradient += m_driver.Gradient() * dEdR;

        if (m_nocolvarfile == false) {
            std::ofstream colvarfile;
            colvarfile.open(m_colvar_base + "_" + std::to_string(m_biased_structures[i].index), std::iostream::app);
            colvarfile << m_currentStep << " " << rmsd << " " << bias_energy << " "
                       << m_biased_structures[i].counter << " " << m_biased_structures[i].factor << std::endl;
            colvarfile.close();
        }

        // Visited if the walker sits inside this Gaussian (heuristic deposition gate).
        if (expr * m_rmsd_econv > static_cast<double>(m_biased_structures.size()))
            visited.push_back(i);

        m_counter += m_biased_structures[i].counter;
    }

    // Phase 2: update visited references. Exploration (counter) is ALWAYS undamped so
    // well-tempering never slows the search. WT (opt-in) only feeds the separate output
    // weight 'factor' with the standard W*exp(-V/(kB*Delta_T)) damping.
    for (int i : visited) {
        m_biased_structures[i].counter++;
        if (m_wtmtd)
            m_biased_structures[i].factor += exp(-m_current_bias / (kb_Eh * m_DT));
    }
    return 1;
}

std::vector<json> BiasThread::getBias() const
{
    std::vector<json> bias(m_biased_structures.size());
    for (int i = 0; i < m_biased_structures.size(); ++i) {
        json current;
        // current["geometry"] = Tools::Matrix2String(m_biased_structures[i].geometry);
        current["time"] = m_biased_structures[i].time;
        current["rmsd_reference"] = m_biased_structures[i].rmsd_reference;
        current["energy"] = m_biased_structures[i].energy;
        current["factor"] = m_biased_structures[i].factor;
        current["index"] = m_biased_structures[i].index;
        current["counter"] = m_biased_structures[i].counter;
        bias[i] = current;
    }
    return bias;
}

SimpleMD::SimpleMD(const json& controller, const bool silent)
    : CurcumaMethod(ParameterRegistry::getInstance().getDefaultJson("simplemd"), controller, silent)
    , m_config("simplemd", controller)  // Claude Generated - ConfigManager initialization
{
    UpdateController(controller);  // This calls LoadControlJson() internally
}

SimpleMD::~SimpleMD()
{
    delete m_interface;
    for (const auto & m_unique_structure : m_unique_structures)
        delete m_unique_structure;
    // delete m_bias_pool;
}

// Claude Generated 2025: String-to-Enum mappings for type-safe selection
static const std::map<std::string, ThermostatType> thermostat_map = {
    {"berendsen", ThermostatType::Berendsen},
    {"berendson", ThermostatType::Berendsen},  // Legacy typo support
    {"anderson", ThermostatType::Anderson},
    {"nosehover", ThermostatType::NoseHover},
    {"csvr", ThermostatType::CSVR},
    {"none", ThermostatType::None}
};

static const std::map<std::string, WallGeometry> wall_geometry_map = {
    {"none", WallGeometry::None},
    {"spheric", WallGeometry::Spheric},
    {"rect", WallGeometry::Rect}
};

static const std::map<std::string, WallPotentialType> wall_potential_map = {
    {"logfermi", WallPotentialType::LogFermi},
    {"harmonic", WallPotentialType::Harmonic}
};

void SimpleMD::LoadControlJson()
{
    try {
    // Claude Generated 2025: Basic Parameters - ConfigManager migration
    m_method = m_config.get<std::string>("method");
    m_thermostat = m_config.get<std::string>("thermostat");
    m_plumed = m_config.get<std::string>("plumed_file");

    m_spin = m_config.get<int>("spin");
    m_charge = m_config.get<int>("charge");
    m_dT = m_config.get<double>("time_step");
    m_maxtime = m_config.get<double>("max_time");
    m_T0 = m_config.get<double>("temperature");
    m_T_init = m_config.get<double>("initial_temperature");
    if (m_T_init < 0) m_T_init = m_T0;
    if (m_T_init <= 0) m_T_init = m_T0; // defensive: reject non-positive explicit values
    m_rmrottrans = m_config.get<int>("remove_com_mode");
    m_nocenter = m_config.get<bool>("no_center");
    m_COM = m_config.get<bool>("use_com");
    m_dump = m_config.get<int>("dump_frequency");
    m_print = m_config.get<int>("print_frequency");
    m_max_top_diff = m_config.get<int>("MaxTopoDiff", 15);  // Not in PARAM block, using default
    m_seed = m_config.get<int>("seed");
    m_threads = m_config.get<int>("threads");

    // Claude Generated 2025: System Control & Output Parameters
    m_rmsd = m_config.get<double>("rmsd", 1.5);  // Not in PARAM block - legacy parameter
    m_hmass = m_config.get<int>("hydrogen_mass");

    m_impuls = m_config.get<double>("impuls", 0.0);  // Not in PARAM block - legacy
    m_impuls_scaling = m_config.get<double>("impuls_scaling", 0.75);  // Not in PARAM block - legacy
    m_writeUnique = m_config.get<bool>("unique", false);  // Not in PARAM block - legacy
    m_opt = m_config.get<bool>("opt", false);  // Not in PARAM block - legacy
    m_scale_velo = m_config.get<double>("initial_velocity_scale");
    m_rescue = m_config.get<bool>("rescue", false);  // Not in PARAM block - legacy
    m_wall_render = m_config.get<bool>("wall_render", false);  // Not in PARAM block - legacy
    m_coupling = m_config.get<double>("coupling");
    m_anderson = m_config.get<double>("anderson_probability");

    if (m_coupling < m_dT)
        m_coupling = m_dT;

    // Claude Generated 2025: RMSD Metadynamics block
    // this one is used to recover https://doi.org/10.1021/acs.jctc.9b00143
    m_rmsd_mtd = m_config.get<bool>("rmsd_mtd");
    m_k_rmsd = m_config.get<double>("rmsd_mtd_k");
    m_alpha_rmsd = m_config.get<double>("rmsd_mtd_alpha");
    m_mtd_steps = m_config.get<int>("rmsd_mtd_pace");
    m_chain_length = m_config.get<int>("chain_length");
    m_rmsd_rmsd = m_config.get<double>("rmsd_rmsd", 1.0);  // Not in PARAM block - legacy
    m_max_rmsd_N = m_config.get<int>("rmsd_mtd_max_gaussians");
    m_rmsd_econv = m_config.get<double>("rmsd_econv", 1e8);  // Not in PARAM block - legacy
    m_rmsd_DT = m_config.get<double>("rmsd_mtd_dt");
    m_rmsd_mtd_max_height = m_config.get<int>("rmsd_mtd_max_height", 0);  // Claude Generated (Jun 2026): cap on hill height
    m_freeze_inherited = m_config.get<bool>("rmsd_mtd_freeze_inherited", false);  // Claude Generated (Jun 2026)
    m_wtmtd = m_config.get<bool>("wtmtd", false);  // Not in PARAM block - legacy
    m_rmsd_ref_file = m_config.get<std::string>("rmsd_mtd_ref_file");
    m_rmsd_fix_structure = m_config.get<bool>("rmsd_fix_structure", false);  // Not in PARAM block - legacy
    m_nocolvarfile = m_config.get<bool>("noCOLVARfile", false);  // Not in PARAM block - legacy
    m_nohillsfile = m_config.get<bool>("noHILSfile", false);  // Not in PARAM block - legacy

    m_rmsd_atoms = m_config.get<std::string>("rmsd_mtd_atoms");

    // Claude Generated (Jun 2026): ConfSearch robustness gates
    m_topo_check = m_config.get<bool>("topo_check", false);
    m_topo_check_interval = m_config.get<int>("topo_check_interval", 0);
    m_epot_abort = m_config.get<bool>("epot_abort", false);
    m_epot_abort_window = m_config.get<double>("epot_abort_window", 250.0);
    m_temp_abort = m_config.get<bool>("temp_abort", false);
    m_temp_abort_factor = m_config.get<double>("temp_abort_factor", 1.5);
    m_temp_abort_delta = m_config.get<double>("temp_abort_delta", 300.0);

    // Claude Generated 2025: Output & Restart Parameters
    m_writerestart = m_config.get<int>("write_restart_frequency");
    m_respa = m_config.get<int>("respa", 1);  // Not in PARAM block - legacy
    m_dipole = m_config.get<bool>("dipole", false);  // Not in PARAM block - legacy
    m_scaling_json = m_config.get<std::string>("scaling_json", "none");  // Not in PARAM block - legacy

    m_writeXYZ = m_config.get<bool>("write_xyz");
    m_writeinit = m_config.get<bool>("write_initial_state");
    m_mtd = m_config.get<bool>("mtd");
    m_mtd_dT = m_config.get<int>("mtd_dT", -1);  // Not in PARAM block - legacy
    if (m_mtd_dT < 0) {
        m_eval_mtd = true;
    } else {
        m_eval_mtd = false;
    }
    m_initfile = m_config.get<std::string>("restart_file");
    m_norestart = m_config.get<bool>("no_restart");
    m_dt2 = m_dT * m_dT;

    // Claude Generated (Nov 2025): CG-specific parameters
    m_cg_write_vtf = m_config.get<bool>("cg_write_vtf");
    bool cg_timestep_scaling_enabled = m_config.get<bool>("cg_timestep_scaling");
    double cg_timestep_factor_config = m_config.get<double>("cg_timestep_factor");

    // Override timestep factor if user disabled scaling
    if (!cg_timestep_scaling_enabled) {
        m_cg_timestep_factor = 1.0;
    } else {
        m_cg_timestep_factor = cg_timestep_factor_config;
    }

    // WP-S2 (May 2026): per-step diagnostics JSONL dump
    m_md_diagnostics = m_config.get<bool>("md_diagnostics");
    // WP-P1 (May 2026): per-phase wall-clock timing in the JSONL records
    m_md_diagnostics_timing = m_config.get<bool>("md_diagnostics_timing");

    // Claude Generated 2025: RATTLE Parameters
    m_rm_COM = m_config.get<double>("remove_com_motion");
    int rattle = m_config.get<int>("rattle");

    m_rattle_maxiter = m_config.get<int>("rattle_max_iterations");
    m_rattle_dynamic_tol_iter = m_config.get<int>("rattle_dynamic_tol_iter", 100);  // Not in PARAM block - legacy
    m_rattle_max = m_config.get<double>("rattle_max", 10.0);  // Not in PARAM block - legacy
    m_rattle_min = m_config.get<double>("rattle_min", 1e-4);  // Not in PARAM block - legacy

    m_rattle_dynamic_tol = m_config.get<bool>("rattle_dynamic_tol", false);  // Not in PARAM block - legacy

    // Claude Generated (Jun 2026): mode 2 (constrain X-H only) was documented (PARAM: 0:off,1:on,
    // 2:H-only) and handled by InitConstrainedBonds (m_rattle==2 path), but this activation gate
    // only fired for ==1, so -rattle 2 silently ran plain Verlet with 0 constraints. Activate the
    // RATTLE integrator for both modes.
    if (rattle == 1 || rattle == 2) {
        Integrator = [=]() {
            this->Rattle();
        };
        m_rattle_tol_12 = m_config.get<double>("rattle_tol_12");
        m_rattle_tol_13 = m_config.get<double>("rattle_tol_13");

        m_rattle_12 = m_config.get<bool>("rattle_12");
        m_rattle_13 = m_config.get<bool>("rattle_13");

        // m_coupling = m_dT;
        m_rattle = m_config.get<int>("rattle");
        std::cout << "Using rattle to constrain bonds!" << std::endl;
        if (m_rattle_12)
            std::cout << "Using rattle to constrain 1,2 distances!" << std::endl;
        if (m_rattle_13)
            std::cout << "Using rattle to constrain 1,3 distances between two bonds!" << std::endl;

    } else {
        Integrator = [=]() {
            this->Verlet();
        };
    }

    // Claude Generated 2025: Energy Calculator Selection
    if (m_config.get<bool>("cleanenergy", false)) {  // Not in PARAM block - legacy
        Energy = [=]() -> double {
            return this->CleanEnergy();
        };
        std::cout << "Energy Calculator will be set up for each step! Single steps are slower, but more reliable. Recommended for the combination of GFN2 and solvation." << std::endl;
    } else {
        Energy = [=]() -> double {
            return this->FastEnergy();
        };
        std::cout << "Energy Calculator will NOT be set up for each step! Fast energy calculation! This is the default way and should not be changed unless the energy and gradient calculation are unstable (happens with GFN2 and solvation)." << std::endl;
    }

    // Claude Generated 2025: Wall Potential Parameters - Enum-based selection
    std::string wall_geom_str = m_config.get<std::string>("wall_type");
    std::string wall_pot_str = m_config.get<std::string>("wall_potential");

    auto geom_it = wall_geometry_map.find(wall_geom_str);
    auto pot_it = wall_potential_map.find(wall_pot_str);

    WallGeometry wall_geom = (geom_it != wall_geometry_map.end()) ? geom_it->second : WallGeometry::None;
    WallPotentialType wall_pot = (pot_it != wall_potential_map.end()) ? pot_it->second : WallPotentialType::Harmonic;

    if (wall_geom == WallGeometry::Spheric) {
        switch (wall_pot) {
            case WallPotentialType::LogFermi:
                m_wall_type = 1;
                WallPotential = [=]() -> double {
                    this->m_wall_potential = this->ApplySphericLogFermiWalls();
                    return m_wall_potential;
                };
                break;
            case WallPotentialType::Harmonic:
                m_wall_type = 1;
                WallPotential = [=]() -> double {
                    this->m_wall_potential = this->ApplySphericHarmonicWalls();
                    return m_wall_potential;
                };
                break;
        }
        std::cout << "Setting up spherical potential" << std::endl;

    } else if (wall_geom == WallGeometry::Rect) {
        switch (wall_pot) {
            case WallPotentialType::LogFermi:
                m_wall_type = 2;
                WallPotential = [=]() -> double {
                    this->m_wall_potential = this->ApplyRectLogFermiWalls();
                    return m_wall_potential;
                };
                break;
            case WallPotentialType::Harmonic:
                m_wall_type = 2;
                WallPotential = [=]() -> double {
                    this->m_wall_potential = this->ApplyRectHarmonicWalls();
                    return m_wall_potential;
                };
                break;
        }
        std::cout << "Setting up rectangular potential" << std::endl;
    } else {
        WallPotential = [=]() -> double {
            return 0;
        };
    }
    m_rm_COM_step = static_cast<int>(m_rm_COM / m_dT);
    } catch (const std::exception& e) {
        std::cerr << "[SimpleMD] FATAL: Parameter conversion error: " << e.what() << std::endl;
        std::cerr << "[SimpleMD] Please check your configuration parameters for type mismatches." << std::endl;
        std::cerr << "[SimpleMD] Config dump: " << m_config.exportConfig().dump(2) << std::endl;
        throw;
    }
}

// Claude Generated 2025: Use ParameterRegistry for automatic help generation

void SimpleMD::printHelp() const
{
    ParameterRegistry::getInstance().printHelp("simplemd");

    std::cout << "\nExample usage:\n"
              << "  curcuma -md input.xyz -max_time 10000 -temperature 300 -thermostat csvr\n"
              << "  curcuma -md input.xyz -method gfn2 -wall_type spheric -wall_radius 10.0\n\n"
              << "Usage Tips:\n"
              << "- For stable dynamics, use time_step ≤ 1.0 fs\n"
              << "- The Berendsen thermostat is efficient but doesn't sample canonical ensemble\n"
              << "- For proper NVT sampling, use CSVR or Nosé-Hoover thermostats\n"
              << "- RATTLE constraints allow larger timesteps for bonds involving H atoms\n"
              << "- Wall potentials prevent molecules from drifting too far\n"
              << "- Metadynamics helps explore conformational space efficiently\n"
              << std::endl;
}

bool SimpleMD::Initialise()
{
    checkHelp();
    // Claude Generated 2026: Create .snapshots subdirectory inside BMT for JSON restarts,
    // but keep m_output_dir pointing at BMT root so trajectory files land there.
    // NOTE: With BMT default-on, OutputDir() is always set. The legacy CWD path (empty OutputDir)
    // is preserved for -no_bmt / -bmt false usage but is effectively unreachable in default config.
    if (OutputDir().empty()) {
        // Legacy path: no BMT directory, snapshots go to basename.snapshots/ in CWD
        setOutputDir(Basename() + ".snapshots");
    } else {
        // BMT is active: m_output_dir stays as BMT root, create snapshots subdirectory separately
        m_snapshots_dir = outputPath(Basename() + ".snapshots");
        ensureOutputDir();  // already ensured by setOutputDir, but ensure snapshots subdir too
#ifdef C17
#ifndef _WIN32
        std::filesystem::create_directories(m_snapshots_dir);
#endif
#endif
    }
    m_natoms = m_molecule.AtomCount();
    if(m_natoms == 0)
        return false;

    m_atomtype = std::vector<int>(m_natoms, 0);

    m_eigen_geometry = Eigen::MatrixXd::Zero(m_natoms, 3);
    m_eigen_geometry_old = Eigen::MatrixXd::Zero(m_natoms, 3);
    m_eigen_gradient = Eigen::MatrixXd::Zero(m_natoms, 3);
    m_eigen_gradient_old = Eigen::MatrixXd::Zero(m_natoms, 3);
    m_eigen_velocities = Eigen::MatrixXd::Zero(m_natoms, 3);

    m_eigen_masses = Eigen::VectorXd::Zero(3*m_natoms);
    m_eigen_inv_masses = Eigen::VectorXd::Zero(3*m_natoms);

    static std::random_device rd{};
    static std::mt19937 gen{ rd() };
    if (m_seed == -1) {
        const auto start = std::chrono::high_resolution_clock::now();
        m_seed = std::chrono::duration_cast<std::chrono::seconds>(start.time_since_epoch()).count();
    } else if (m_seed == 0)
        m_seed = m_natoms * m_T0;
    std::cout << "Random seed is " << m_seed << std::endl;
    gen.seed(m_seed);

    if (m_initfile != "none") {
        json md;
        std::ifstream restart_file(m_initfile);
        try {
            restart_file >> md;
        } catch ([[maybe_unused]] nlohmann::json::type_error& e) {
            throw 404;
        } catch ([[maybe_unused]] nlohmann::json::parse_error& e) {
            throw 404;
        }
        LoadRestartInformation(md);

    } else if (!m_restart && !m_norestart)
        LoadRestartInformation();

    if (m_molecule.AtomCount() == 0)
        return false;

    // Claude Generated (Oct 2025): Detect CG system and PBC for optimization
    m_is_cg_system = m_molecule.isCGSystem();
    m_has_pbc = m_molecule.hasPBC();
    std::vector<int> cg_atoms = m_molecule.getCGAtoms();
    m_cg_atom_count = cg_atoms.size();

    // CG systems can use larger timesteps (10x typical atomic MD)
    if (m_is_cg_system && m_cg_atom_count == m_natoms) {
        m_cg_timestep_factor = 10.0;
        int verbosity = m_config.get<int>("verbosity", 0);
        if (verbosity >= 1) {
            CurcumaLogger::success("Pure CG system detected with "
                                  + std::to_string(m_cg_atom_count) + " CG particles");
            CurcumaLogger::info("Timestep factor: " + std::to_string(m_cg_timestep_factor) + "x");
        }
        // Apply timestep scaling after detection
        if (m_cg_timestep_factor > 1.0) {
            double orig_dt = m_dT;
            m_dT *= m_cg_timestep_factor;
            m_dt2 = m_dT * m_dT;
            int verbosity_ts = m_config.get<int>("verbosity", 0);
            if (verbosity_ts >= 1) {
                CurcumaLogger::success("CG timestep scaling applied: "
                                      + std::to_string(orig_dt) + " fs → "
                                      + std::to_string(m_dT) + " fs");
            }
        }
    }

    if (m_has_pbc) {
        int verbosity = m_config.get<int>("verbosity", 0);
        if (verbosity >= 2) {
            Eigen::Matrix3d cell = m_molecule.getUnitCell();
            CurcumaLogger::param("Periodic Boundary Conditions",
                               "enabled (box: " + std::to_string(cell(0,0)) + " x " +
                               std::to_string(cell(1,1)) + " x " +
                               std::to_string(cell(2,2)) + " Å³)");
        }
    }

    // Claude Generated (Nov 2025): Initialize orientational dynamics infrastructure
    if (m_is_cg_system) {
        m_cg_orientations.resize(m_natoms, Eigen::Vector3d(1.0, 0.0, 0.0)); // Default: x-axis
        m_cg_angular_velocities.resize(m_natoms, Eigen::Vector3d::Zero());

        // Future Phase 6: Load orientations from VTF/JSON if ellipsoids detected
        // For now, spheres don't need orientation tracking (m_cg_enable_rotation = false)
        m_cg_enable_rotation = false;

        int verbosity = m_config.get<int>("verbosity", 0);
        if (verbosity >= 2) {
            CurcumaLogger::info("Orientational dynamics infrastructure initialized (spheres mode)");
        }
    }

    if (!m_restart) {
        std::ofstream result_file;
        result_file.open(outputPath(Basename() + ".trj.xyz"));
        result_file.close();
    }

    if (m_seed == -1) {
        const auto start = std::chrono::high_resolution_clock::now();
        m_seed = std::chrono::duration_cast<std::chrono::seconds>(start.time_since_epoch()).count();
    } else if (m_seed == 0)
        m_seed = m_T0 * m_natoms;
    std::cout << "Random seed is " << m_seed << std::endl;
    gen.seed(m_seed);


    m_start_fragments = m_molecule.GetFragments();
    m_start_fragment_count = static_cast<int>(m_start_fragments.size());  // Claude Generated (Jun 2026): topo-check reference
    m_scaling_vector_linear = std::vector<double>(m_natoms, 1);
    m_scaling_vector_nonlinear = std::vector<double>(m_natoms, 1);
    if (m_scaling_json != "none") {
        json scaling;
        std::ifstream file(m_scaling_json);
        try {
            file >> scaling;
        } catch ([[maybe_unused]] nlohmann::json::type_error& e) {
            throw 404;
        } catch ([[maybe_unused]] nlohmann::json::parse_error& e) {
            throw 404;
        }
        std::string scaling_vector_linear, scaling_vector_nonlinear;
        try {
            scaling_vector_linear = scaling["scaling_vector_linear"];
            scaling_vector_nonlinear = scaling["scaling_vector_nonlinear"];
        } catch ([[maybe_unused]] json::type_error& e) {
        }
        if (!scaling_vector_linear.empty()) {
            m_scaling_vector_linear = Tools::String2DoubleVec(scaling_vector_linear, "|");
        }
        if (!scaling_vector_nonlinear.empty()) {
            m_scaling_vector_nonlinear = Tools::String2DoubleVec(scaling_vector_nonlinear, "|");
        }
    }

    m_molecule.setCharge(0);
    if (!m_nocenter) {
        std::cout << "Move stucture to the origin ... " << std::endl;
        m_molecule.Center(m_COM);
    } else
        std::cout << "Move stucture NOT to the origin ... " << std::endl;



    if (!m_restart) {
        m_eigen_geometry = Eigen::MatrixXd::Zero(m_natoms, 3);
        m_eigen_velocities = Eigen::MatrixXd::Zero(m_natoms, 3);
        m_currentStep = 0;
    }
    /* */
    m_rt_geom_1 = std::vector<double>(3 * m_natoms, 0);
    m_rt_geom_2 = std::vector<double>(3 * m_natoms, 0);
    m_rt_velo = std::vector<double>(3 * m_natoms, 0);
    /* */

    //m_gradient = std::vector<double>(3 * m_natoms, 0);
    m_virial = std::vector<double>(3 * m_natoms, 0);
    m_atom_temp = std::vector<std::vector<double>>(m_natoms);
    if(m_opt)
    {
        // Claude Generated (Apr 2026): Use unified optimizer instead of legacy CurcumaOpt
        json opt_config;
        opt_config["method"] = m_method;
        opt_config["write_trajectory"] = false;

        EnergyCalculator energy_calc(m_method, m_defaults);
        auto result = Optimization::OptimizationDispatcher::optimizeStructure(
            &m_molecule, Optimization::OptimizerType::LBFGSPP, &energy_calc, opt_config);

        if (result.success) {
            m_molecule.setGeometry(result.final_molecule.getGeometry());
        }
        m_molecule.appendXYZFile(outputPath(Basename() + ".opt.xyz"));
    }
    double mass = 0;
    for (int i = 0; i < m_natoms; ++i) {
        m_atomtype[i] = m_molecule.Atom(i).first;
        if (!m_restart) {
            Position pos = m_molecule.Atom(i).second;
            // std::cout << pos << std::endl;
            m_eigen_geometry.data()[3 * i + 0] = pos(0);
            m_eigen_geometry.data()[3 * i + 1] = pos(1);
            m_eigen_geometry.data()[3 * i + 2] = pos(2);
            /*
            m_eigen_geometry(i, 0) = pos(0) / 1;
            m_eigen_geometry(i, 1) = pos(1) / 1;
            m_eigen_geometry(i, 2) = pos(2) / 1;
            */
        }
        if (m_atomtype[i] == 1) {
            m_eigen_masses.data()[3 * i + 0] = Elements::AtomicMass[m_atomtype[i]] * m_hmass;
            m_eigen_masses.data()[3 * i + 1] = Elements::AtomicMass[m_atomtype[i]] * m_hmass;
            m_eigen_masses.data()[3 * i + 2] = Elements::AtomicMass[m_atomtype[i]] * m_hmass;

            m_eigen_masses(3*i) = Elements::AtomicMass[m_atomtype[i]] * m_hmass;
            m_eigen_masses(3*i+1) = Elements::AtomicMass[m_atomtype[i]] * m_hmass;
            m_eigen_masses(3*i+2) = Elements::AtomicMass[m_atomtype[i]] * m_hmass;

            m_eigen_inv_masses(3*i) = 1 / (Elements::AtomicMass[m_atomtype[i]] * m_hmass);
            m_eigen_inv_masses(3*i + 1) = 1 / (Elements::AtomicMass[m_atomtype[i]] * m_hmass);
            m_eigen_inv_masses(3*i + 2) = 1 / (Elements::AtomicMass[m_atomtype[i]] * m_hmass);  

            mass += Elements::AtomicMass[m_atomtype[i]] * m_hmass;

            m_eigen_inv_masses.data()[3 * i + 0] = 1 / m_eigen_masses.data()[3 * i + 0];
            m_eigen_inv_masses.data()[3 * i + 1] = 1 / m_eigen_masses.data()[3 * i + 1];
            m_eigen_inv_masses.data()[3 * i + 2] = 1 / m_eigen_masses.data()[3 * i + 2];
        } else {
            m_eigen_masses.data()[3 * i + 0] = Elements::AtomicMass[m_atomtype[i]];
            m_eigen_masses.data()[3 * i + 1] = Elements::AtomicMass[m_atomtype[i]];
            m_eigen_masses.data()[3 * i + 2] = Elements::AtomicMass[m_atomtype[i]];
            mass += Elements::AtomicMass[m_atomtype[i]];

            m_eigen_inv_masses.data()[3 * i + 0] = 1 / m_eigen_masses.data()[3 * i + 0];
            m_eigen_inv_masses.data()[3 * i + 1] = 1 / m_eigen_masses.data()[3 * i + 1];
            m_eigen_inv_masses.data()[3 * i + 2] = 1 / m_eigen_masses.data()[3 * i + 2];

            m_eigen_masses(3*i) = Elements::AtomicMass[m_atomtype[i]];
            m_eigen_masses(3*i+1) = Elements::AtomicMass[m_atomtype[i]];
            m_eigen_masses(3*i+2) = Elements::AtomicMass[m_atomtype[i]];

            m_eigen_inv_masses(3*i) = 1 / (Elements::AtomicMass[m_atomtype[i]]);
            m_eigen_inv_masses(3*i + 1) = 1 / (Elements::AtomicMass[m_atomtype[i]]);
            m_eigen_inv_masses(3*i + 2) = 1 / (Elements::AtomicMass[m_atomtype[i]]);    
        }
    }
    // std::cout << m_eigen_geometry << std::endl;
    // std::cout << m_eigen_masses << std::endl;
    m_molecule.setCharge(m_charge);
    m_molecule.setSpin(m_spin);
    // Claude Generated (October 2025 - FIXED): Use controller["simplemd"] like CurcumaOpt uses controller["opt"]
    // Also merge in global parameters and defaults from ParameterRegistry
    json ec_config = m_controller.contains("simplemd") && m_controller["simplemd"].is_object()
        ? m_controller["simplemd"]
        : json::object();

    // Forward flat user parameters from controller to energy calculator config.
    // ConfSearch and other callers may pass parameters at the top level rather
    // than nested inside controller["simplemd"]. Without this forwarding, keys
    // like "gpu" or "charge" are silently dropped.
    for (auto& [key, value] : m_controller.items()) {
        if (key == "simplemd" || key == "global" || value.is_object())
            continue;
        if (!ec_config.contains(key) || ec_config[key].is_null()) {
            ec_config[key] = value;
        }
    }

    // Merge global parameters as fallback
    if (m_controller.contains("global") && m_controller["global"].is_object()) {
        for (auto& [key, value] : m_controller["global"].items()) {
            if (!ec_config.contains(key)) {
                ec_config[key] = value;
            }
        }
    }

    // WP-S/CLI-routing (May 2026): forward method-specific sub-scopes (gfnff, eeq_solver,
    // tblite, xtb, ...) routed by CLI2Json::findOwnerModules. Without this the GFNFF
    // constructor never sees flat flags like -static_all true or -eeq_distance_cutoff_auto
    // true that the registry routed into controller["gfnff"]. Mirrors the
    // EnergyCalculator::reattachMethodScopes fix used by the opt/sp path (WP6).
    static const std::vector<std::string> kMethodScopes = {
        "gfnff", "eeq_solver", "tblite", "xtb", "ulysses", "eht", "dftd3", "dftd4", "orca"
    };
    for (const auto& scope : kMethodScopes) {
        if (m_controller.contains(scope) && m_controller[scope].is_object()
            && !ec_config.contains(scope)) {
            ec_config[scope] = m_controller[scope];
        }
    }

    // Ensure critical parameters exist with defaults from ParameterRegistry
    json defaults = ParameterRegistry::getInstance().getDefaultJson("energycalculator");
    for (auto& [key, value] : defaults.items()) {
        if (!ec_config.contains(key) || ec_config[key].is_null()) {
            ec_config[key] = value;
        }
    }

    m_interface = new EnergyCalculator(m_method, ec_config, Basename());

    m_interface->setMolecule(m_molecule.getMolInfo());
    // Energy-method-setup boundary (Claude Generated, Jun 2026): EnergyCalculator construction +
    // setMolecule (e.g. GFN-FF parameter generation) leaves the global CurcumaLogger verbosity
    // clamped to 0 (it captures/restores around an already-clamped level), so the remaining setup
    // output below — notably the RATTLE report in InitConstrainedBonds — was silently dropped. The
    // CurcumaMethod base RAII can't help (EnergyCalculator is not a CurcumaMethod), so re-assert our
    // level here. (A deeper fix would stop the EnergyCalculator setup path leaking 0 in the first
    // place.)
    CurcumaLogger::set_verbosity(m_verbosity);
    // Iterative mode: raise SCF display threshold by one level so system verbosity
    // controls output (silent at default=1, visible at -v 2). Claude Generated.
    m_interface->setIterativeMode(true);
    // Warm-start: reuse converged shell charges from the previous geometry step.
    // Reduces SCF iterations for native GFN1/GFN2 in MD; no-op for other methods.
    if (m_method == "gfn1" || m_method == "gfn2")
        m_interface->setWarmStart(true);
    // m_interface->setGeometryFile(Basename() + ".xyz"); TODO this does not really work
    // m_interface->setBasename(Basename()); TODO this does not really work
    if (m_writeUnique) {
        // Claude Generated: Replace static RMSDTrajJson with ParameterRegistry
        json rmsdtraj = ParameterRegistry::getInstance().getDefaultJson("rmsdtraj");
        rmsdtraj["writeUnique"] = true;
        rmsdtraj["rmsd"] = m_rmsd;
        rmsdtraj["writeRMSD"] = false;
        m_unqiue = new RMSDTraj(rmsdtraj, true);
        m_unqiue->setBaseName(outputPath(Basename() + ".xyz"));
        m_unqiue->Initialise();
    }
    m_dof = 3 * m_natoms;
    InitialiseWalls();
    if (!m_restart) {

        InitConstrainedBonds();
        InitVelocities(m_scale_velo);
        m_xi.resize(m_chain_length, 0.0);
        m_Q.resize(m_chain_length, 100); // Setze eine geeignete Masse für jede Kette
        for (int i = 0; i < m_chain_length; ++i) {
            m_xi[i] = pow(10.0, static_cast<double>(i)) - 1;
            m_Q[i] = pow(10, i) * kb_Eh * m_T0 * m_dof * 100;
            // std::cout << m_xi[i] << "  " << m_Q[i] << std::endl;
        }
        m_eta = 0.0;
    }
    if (m_writeinit) {
        json init = WriteRestartInformation();
        std::ofstream result_file;
        result_file.open(snapshotPath(Basename() + ".init.json"));
        result_file << init;
        result_file.close();
    }
    /* Initialising MTD RMSD Threads */
    if (m_rmsd_mtd) {
        m_bias_pool = new CxxThreadPool;
        m_bias_pool->setProgressBar(CxxThreadPool::ProgressBarType::None);
        m_bias_pool->setActiveThreadCount(m_threads);
        m_molecule.GetFragments();
        m_rmsd_indicies = m_molecule.FragString2Indicies(m_rmsd_atoms);

        for(auto i : m_rmsd_indicies)
        {
            // std::cout << i << " ";
            m_rmsd_mtd_molecule.addPair(m_molecule.Atom(i));
        }
        m_rmsd_fragment_count = m_rmsd_mtd_molecule.GetFragments().size();

        json config = ParameterRegistry::getInstance().getDefaultJson("rmsd");  // Claude Generated 2025: Use ParameterRegistry instead of RMSDJson
        config["silent"] = true;
        config["reorder"] = false;

        // Claude Generated (Apr 2026): Initialize shared pool RMSDDriver for parallel ConfSearch
        m_shared_pool_driver = RMSDDriver(config, true);
        m_shared_pool_driver.setReference(m_rmsd_mtd_molecule);
        m_shared_pool_target = m_rmsd_mtd_molecule;

        for (int i = 0; i < m_threads; ++i) {
            auto* thread = new BiasThread(m_rmsd_mtd_molecule, config, m_nocolvarfile, m_nohillsfile, outputPath("COLVAR"));
            thread->setDT(m_rmsd_DT);
            thread->setk(m_k_rmsd);
            thread->setalpha(m_alpha_rmsd);
            thread->setEnergyConv(m_rmsd_econv);
            thread->setWTMTD(m_wtmtd);
            m_bias_threads.push_back(thread);
            m_bias_pool->addThread(thread);
        }
        if (m_restart) {
            std::cout << "Reading structure files from " << m_rmsd_ref_file << std::endl;
            for (const auto& i : m_bias_json)
                std::cout << i << std::endl;
            FileIterator file(m_rmsd_ref_file);
            int index = 0;
            while (!file.AtEnd()) {
                Molecule mol = file.Next();
                std::cout << m_bias_json[index] << std::endl;
                int thread_index = index % m_bias_threads.size();
                m_bias_threads[thread_index]->addGeometry(mol.getGeometry(), m_bias_json[index]);
                ++index;
            }
            m_bias_structure_count = index;
        } else {
            if (m_rmsd_ref_file != "none") {
                std::cout << "Reading structure files from " << m_rmsd_ref_file << std::endl;
                int index = 0;

                FileIterator file(m_rmsd_ref_file);
                while (!file.AtEnd()) {
                    Molecule mol = file.Next();
                    int thread_index = index % m_bias_threads.size();
                    m_bias_threads[thread_index]->addGeometry(mol.getGeometry(), 0, 0, index);
                    ++index;
                }
                m_bias_structure_count = index;
            }
        }
    }

    m_initialised = true;
    return true;
}

void SimpleMD::InitConstrainedBonds()
{
    int total_bonds = 0; // Claude Generated (Jun 2026): all 1-2 bonds in the topology (for the summary)
    if (m_rattle) {
        auto m = m_molecule.DistanceMatrix();
        m_topo_initial = m.second;
        for (int i = 0; i < m_molecule.AtomCount(); ++i) {
            for (int j = 0; j < i; ++j) {
                if (m.second(i, j)) {
                    ++total_bonds;
                    if (m_rattle == 2) {
                        if (m_molecule.Atom(i).first != 1 && m_molecule.Atom(j).first != 1)
                            continue;
                    }
                    std::pair<int, int> indicies(i, j);
                    std::pair<double, double> minmax(m_molecule.CalculateDistance(i, j) * m_molecule.CalculateDistance(i, j), m_molecule.CalculateDistance(i, j) * m_molecule.CalculateDistance(i, j));
                    std::pair<std::pair<int, int>, double> bond(indicies, m_molecule.CalculateDistance(i, j) * m_molecule.CalculateDistance(i, j));
                    if (m_rattle_12) {
                        m_bond_constrained.emplace_back(bond);
                    }

                    for (int k = 0; k < j; ++k) {
                        if (m.second(k, j)) {
                            std::pair<int, int> indicies(i, k);
                            std::pair<double, double> minmax(m_molecule.CalculateDistance(i, k) * m_molecule.CalculateDistance(i, k), m_molecule.CalculateDistance(i, k) * m_molecule.CalculateDistance(i, k));
                            std::pair<std::pair<int, int>, double> bond(indicies, m_molecule.CalculateDistance(i, k) * m_molecule.CalculateDistance(i, k));
                            if (m_rattle_13) {
                                m_bond_13_constrained.push_back(std::pair<std::pair<int, int>, double>(bond));
                            }
                        }
                    }
                }
                }
            }
    }

    // Step 1: subtract RATTLE constraints (each bond/angle constraint = 1 DOF)
    int n_constraints = static_cast<int>(m_bond_constrained.size() + m_bond_13_constrained.size());
    const int total_dof = m_dof; // 3N before any correction
    m_dof -= n_constraints;
    if (m_dof < 1) m_dof = 1;
    const int dof_after_rattle = m_dof;

    // Step 2: subtract frozen COM/rotation modes. Both RemoveRotation() and
    // RemoveRotations() zero translation AND rotation simultaneously (the
    // mode labels in remove_com_mode are misleading — all non-zero modes
    // remove both). Frozen DOF must be subtracted so that
    //   T = 2*Ekin / (kB * m_dof)
    // reflects only the active internal modes, otherwise the thermostat
    // overdrives kinetic energy and the reported temperature is wrong.
    // Non-linear assumption for 3+ atoms (linear-check is too expensive here).
    int dof_com_removed = 0;
    if (m_rmrottrans > 0) {
        auto dof_for_fragment = [](size_t n) -> int {
            if (n == 1) return 3;  // translation only (no rotational DOF)
            if (n == 2) return 5;  // 3 trans + 2 rot (linear diatomic)
            return 6;              // 3 trans + 3 rot (non-linear)
        };
        if (m_rmrottrans == 1) {
            dof_com_removed = dof_for_fragment(static_cast<size_t>(m_natoms));
        } else {
            for (const auto& frag : m_molecule.GetFragments())
                dof_com_removed += dof_for_fragment(frag.size());
        }
        m_dof -= dof_com_removed;
        if (m_dof < 1) m_dof = 1;
    }

    // Report
    if (m_rattle) {
        CurcumaLogger::result_fmt("RATTLE: {} constraints | {} of {} 1-2 bonds{} + {} 1-3 angles | DOF {} -> {} ({:+d})",
            n_constraints, m_bond_constrained.size(), total_bonds, m_rattle == 2 ? " (X-H only)" : "",
            m_bond_13_constrained.size(), total_dof, dof_after_rattle, dof_after_rattle - total_dof);
        if (CurcumaLogger::get_verbosity() >= 3) {
            if (!m_bond_constrained.empty()) {
                std::string line = "  1-2:";
                int col = 0;
                for (const auto& b : m_bond_constrained) {
                    int i = b.first.first, j = b.first.second;
                    line += fmt::format(" {}{}-{}{}({:.3f})",
                        Elements::ElementAbbr[m_molecule.Atom(i).first], i,
                        Elements::ElementAbbr[m_molecule.Atom(j).first], j, std::sqrt(b.second));
                    if (++col % 5 == 0 && col < static_cast<int>(m_bond_constrained.size()))
                        line += "\n     ";
                }
                CurcumaLogger::result(line);
            }
            if (!m_bond_13_constrained.empty()) {
                std::string line = "  1-3:";
                int col = 0;
                for (const auto& b : m_bond_13_constrained) {
                    line += fmt::format(" {}-{}({:.3f})", b.first.first, b.first.second, std::sqrt(b.second));
                    if (++col % 6 == 0 && col < static_cast<int>(m_bond_13_constrained.size()))
                        line += "\n     ";
                }
                CurcumaLogger::result(line);
            }
        }
    }
    if (dof_com_removed > 0)
        CurcumaLogger::result_fmt("COM/rot removal (mode {}): -{} DOF | effective DOF = {}",
            m_rmrottrans, dof_com_removed, m_dof);
    else if (!m_rattle)
        CurcumaLogger::result_fmt("{} degrees of freedom (no constraints)", m_dof);
}

void SimpleMD::InitVelocities(double scaling)
{
    static std::default_random_engine generator;
    for (size_t i = 0; i < m_natoms; ++i) {
        // Claude Generated (Jun 2026): sample from m_T_init (initial
        // temperature) rather than m_T0 (thermostat target) so callers can
        // anneal into the target temperature or start cold/warm without
        // touching the thermostat target. m_T_init defaults to m_T0 when
        // -initial_temperature is not set (backward compatible).
        std::normal_distribution<double> distribution(0.0, std::sqrt(kb_Eh * m_T_init * m_eigen_inv_masses.data()[3 * i]));
        m_eigen_velocities.data()[3 * i + 0] = distribution(generator);
        m_eigen_velocities.data()[3 * i + 1] = distribution(generator);
        m_eigen_velocities.data()[3 * i + 2] = distribution(generator);
    }

    // Match per-step removal logic exactly so initial velocities are
    // consistent with what the integrator loop enforces each step.
    if (m_rmrottrans == 1)
        RemoveRotation();
    else if (m_rmrottrans == 2)
        RemoveRotations();
    else if (m_rmrottrans == 3) {
        RemoveRotations();
        RemoveRotation();
    }
    EKin();
    // Normalize initial velocities to the requested sampling temperature.
    // When T_init == T0 (default): use two tight Berendson steps to remove
    // statistical fluctuations from the MB draw (original behavior).
    // When T_init != T0: simple velocity rescaling to exactly T_init — the
    // thermostat will then drive toward T0 during the run. Calling Berendson
    // here (which targets m_T0) would immediately destroy the T_init setting.
    if (m_T_init == m_T0) {
        double coupling = m_coupling;
        m_coupling = m_dT;
        Berendson();
        Berendson();
        EKin();
        m_coupling = coupling;
    } else if (m_T > 0.0) {
        double scale = std::sqrt(m_T_init / m_T);
        for (int i = 0; i < m_natoms; ++i) {
            m_eigen_velocities.data()[3 * i + 0] *= scale;
            m_eigen_velocities.data()[3 * i + 1] *= scale;
            m_eigen_velocities.data()[3 * i + 2] *= scale;
        }
        EKin();
    }

    // If RATTLE is active, project velocities onto constraint manifold
    // and rescale to target temperature using reduced DOF.
    // This prevents the initial temperature drop caused by SHAKE
    // removing kinetic energy from constrained DOF.
    if (m_rattle) {
        // Apply RATTLE velocity constraint: project out r_ij · v_ij for each constraint
        double max_mu = 10;
        int iter = 0;
        while (iter < m_rattle_maxiter) {
            iter++;
            int active = 0;
            for (auto bond : m_bond_constrained) {
                int i = bond.first.first, j = bond.first.second;
                double dx = m_eigen_geometry.data()[3 * i + 0] - m_eigen_geometry.data()[3 * j + 0];
                double dy = m_eigen_geometry.data()[3 * i + 1] - m_eigen_geometry.data()[3 * j + 1];
                double dz = m_eigen_geometry.data()[3 * i + 2] - m_eigen_geometry.data()[3 * j + 2];
                double dvx = m_eigen_velocities.data()[3 * i + 0] - m_eigen_velocities.data()[3 * j + 0];
                double dvy = m_eigen_velocities.data()[3 * i + 1] - m_eigen_velocities.data()[3 * j + 1];
                double dvz = m_eigen_velocities.data()[3 * i + 2] - m_eigen_velocities.data()[3 * j + 2];
                double distance_sq = dx * dx + dy * dy + dz * dz;
                double r = dx * dvx + dy * dvy + dz * dvz;
                double mu = -r / ((m_eigen_inv_masses.data()[3 * i] + m_eigen_inv_masses.data()[3 * j]) * distance_sq);
                while (std::abs(mu) > max_mu) mu /= 2;
                if (std::abs(mu) > m_rattle_tol_12) {
                    active = 1;
                    m_eigen_velocities.data()[3 * i + 0] += dx * mu * m_eigen_inv_masses.data()[3 * i];
                    m_eigen_velocities.data()[3 * i + 1] += dy * mu * m_eigen_inv_masses.data()[3 * i];
                    m_eigen_velocities.data()[3 * i + 2] += dz * mu * m_eigen_inv_masses.data()[3 * i];
                    m_eigen_velocities.data()[3 * j + 0] -= dx * mu * m_eigen_inv_masses.data()[3 * j];
                    m_eigen_velocities.data()[3 * j + 1] -= dy * mu * m_eigen_inv_masses.data()[3 * j];
                    m_eigen_velocities.data()[3 * j + 2] -= dz * mu * m_eigen_inv_masses.data()[3 * j];
                }
            }
            for (auto bond : m_bond_13_constrained) {
                int i = bond.first.first, j = bond.first.second;
                double dx = m_eigen_geometry.data()[3 * i + 0] - m_eigen_geometry.data()[3 * j + 0];
                double dy = m_eigen_geometry.data()[3 * i + 1] - m_eigen_geometry.data()[3 * j + 1];
                double dz = m_eigen_geometry.data()[3 * i + 2] - m_eigen_geometry.data()[3 * j + 2];
                double dvx = m_eigen_velocities.data()[3 * i + 0] - m_eigen_velocities.data()[3 * j + 0];
                double dvy = m_eigen_velocities.data()[3 * i + 1] - m_eigen_velocities.data()[3 * j + 1];
                double dvz = m_eigen_velocities.data()[3 * i + 2] - m_eigen_velocities.data()[3 * j + 2];
                double distance_sq = dx * dx + dy * dy + dz * dz;
                double r = dx * dvx + dy * dvy + dz * dvz;
                double mu = -r / ((m_eigen_inv_masses.data()[3 * i] + m_eigen_inv_masses.data()[3 * j]) * distance_sq);
                while (std::abs(mu) > max_mu) mu /= 2;
                if (std::abs(mu) > m_rattle_tol_13) {
                    active = 1;
                    m_eigen_velocities.data()[3 * i + 0] += dx * mu * m_eigen_inv_masses.data()[3 * i];
                    m_eigen_velocities.data()[3 * i + 1] += dy * mu * m_eigen_inv_masses.data()[3 * i];
                    m_eigen_velocities.data()[3 * i + 2] += dz * mu * m_eigen_inv_masses.data()[3 * i];
                    m_eigen_velocities.data()[3 * j + 0] -= dx * mu * m_eigen_inv_masses.data()[3 * j];
                    m_eigen_velocities.data()[3 * j + 1] -= dy * mu * m_eigen_inv_masses.data()[3 * j];
                    m_eigen_velocities.data()[3 * j + 2] -= dz * mu * m_eigen_inv_masses.data()[3 * j];
                }
            }
            if (active == 0) break;
        }
        // Rescale velocities to target temperature using reduced DOF
        EKin();
        double scale = std::sqrt(m_T0 / m_T);
        for (int i = 0; i < 3 * m_natoms; ++i) {
            m_eigen_velocities.data()[i] *= scale;
        }
        EKin();
    }
}

void SimpleMD::InitialiseWalls()
{
    /*
    { "wall_xl", 0},
    { "wall_yl", 0},
    { "wall_zl", 0},
    { "wall_x_min", 0},
    { "wall_x_max", 0},
    { "wall_y_min", 0},
    { "wall_y_max", 0},
    { "wall_z_min", 0},
    { "wall_z_max", 0},*/
    // Claude Generated 2025: Wall Boundary Parameters
    m_wall_spheric_radius = m_config.get<double>("wall_radius");
    m_wall_temp = m_config.get<double>("wall_temp");
    m_wall_beta = m_config.get<double>("wall_beta");

    m_wall_x_min = m_config.get<double>("wall_x_min");
    m_wall_x_max = m_config.get<double>("wall_x_max");
    m_wall_y_min = m_config.get<double>("wall_y_min");
    m_wall_y_max = m_config.get<double>("wall_y_max");
    m_wall_z_min = m_config.get<double>("wall_z_min");
    m_wall_z_max = m_config.get<double>("wall_z_max");
    // Claude Generated: Intelligent auto-sizing - only when boundaries are undefined or invalid
    double radius = 0;
    bool auto_configured = false;

    // Only auto-configure if bounds are completely undefined (0,0) or clearly invalid (max <= min)
    bool x_needs_config = (m_wall_x_min == 0 && m_wall_x_max == 0) || (m_wall_x_max <= m_wall_x_min);
    bool y_needs_config = (m_wall_y_min == 0 && m_wall_y_max == 0) || (m_wall_y_max <= m_wall_y_min);
    bool z_needs_config = (m_wall_z_min == 0 && m_wall_z_max == 0) || (m_wall_z_max <= m_wall_z_min);
    bool sphere_needs_config = (m_wall_spheric_radius == 0);

    if (x_needs_config || y_needs_config || z_needs_config || sphere_needs_config) {
        auto_configured = true;

        // Find actual molecular extent by analyzing all atom positions
        double min_x = 1e10, max_x = -1e10;
        double min_y = 1e10, max_y = -1e10;
        double min_z = 1e10, max_z = -1e10;
        double max_distance_from_origin = 0;

        for (int i = 0; i < m_natoms; ++i) {
            double x = m_eigen_geometry.data()[3 * i + 0];
            double y = m_eigen_geometry.data()[3 * i + 1];
            double z = m_eigen_geometry.data()[3 * i + 2];

            // Track coordinate ranges
            min_x = std::min(min_x, x);
            max_x = std::max(max_x, x);
            min_y = std::min(min_y, y);
            max_y = std::max(max_y, y);
            min_z = std::min(min_z, z);
            max_z = std::max(max_z, z);

            // Track maximum distance from origin for spherical walls
            double distance = std::sqrt(x * x + y * y + z * z);
            max_distance_from_origin = std::max(max_distance_from_origin, distance);
        }

        // Add safety margins (20% + 5 Å minimum)
        double margin_x = std::max(0.2 * (max_x - min_x), 5.0);
        double margin_y = std::max(0.2 * (max_y - min_y), 5.0);
        double margin_z = std::max(0.2 * (max_z - min_z), 5.0);
        double margin_sphere = std::max(0.2 * max_distance_from_origin, 5.0);

        // Only set boundaries that actually need configuration
        if (x_needs_config) {
            m_wall_x_min = min_x - margin_x;
            m_wall_x_max = max_x + margin_x;
        }
        if (y_needs_config) {
            m_wall_y_min = min_y - margin_y;
            m_wall_y_max = max_y + margin_y;
        }
        if (z_needs_config) {
            m_wall_z_min = min_z - margin_z;
            m_wall_z_max = max_z + margin_z;
        }

        // Set spherical radius with margin
        if (sphere_needs_config) {
            radius = max_distance_from_origin + margin_sphere;
        }
    }

    // Fallback to old box-based method if molecule geometry isn't available
    if (m_natoms == 0 && auto_configured) {
        std::vector<double> box = m_molecule.GetBox();
        if (x_needs_config) {
            m_wall_x_min = -box[0] * 0.75;
            m_wall_x_max = -1 * m_wall_x_min;
            radius = std::max(radius, box[0]);
        }
        if (y_needs_config) {
            m_wall_y_min = -box[1] * 0.75;
            m_wall_y_max = -1 * m_wall_y_min;
            radius = std::max(radius, box[1]);
        }
        if (z_needs_config) {
            m_wall_z_min = -box[2] * 0.75;
            m_wall_z_max = -1 * m_wall_z_min;
            radius = std::max(radius, box[2]);
        }
        radius += 5;
    }

    if (sphere_needs_config && m_wall_spheric_radius < radius) {
        m_wall_spheric_radius = radius;
    }
    if (m_wall_render) {
        std::cout << "render walls" << std::endl;
        if (m_wall_type == 1) {
            Position x0 = Position{ m_wall_spheric_radius, 0, 0 };
            Position x1 = Position{ -m_wall_spheric_radius, 0, 0 };
            Position y0 = Position{ 0, m_wall_spheric_radius, 0 };
            Position y1 = Position{ 0, -m_wall_spheric_radius, 0 };
            Position z0 = Position{ 0, 0, m_wall_spheric_radius };
            Position z1 = Position{ 0, 0, -m_wall_spheric_radius };
            m_molecule.addBorderPoint(x0);
            m_molecule.addBorderPoint(x1);
            m_molecule.addBorderPoint(y0);
            m_molecule.addBorderPoint(y1);
            m_molecule.addBorderPoint(z0);
            m_molecule.addBorderPoint(z1);

            double intermedia = 1 / sqrt(2.0) * m_wall_spheric_radius;
            x0 = Position{ intermedia, intermedia, 0 };
            y0 = Position{ 0, intermedia, intermedia };
            z0 = Position{ intermedia, 0, intermedia };

            m_molecule.addBorderPoint(x0);
            m_molecule.addBorderPoint(y0);
            m_molecule.addBorderPoint(z0);
            x0 = Position{ -intermedia, -intermedia, 0 };
            y0 = Position{ 0, -intermedia, -intermedia };
            z0 = Position{ -intermedia, 0, -intermedia };
            m_molecule.addBorderPoint(x0);
            m_molecule.addBorderPoint(y0);
            m_molecule.addBorderPoint(z0);
            x0 = Position{ -intermedia, intermedia, 0 };
            y0 = Position{ 0, -intermedia, intermedia };
            z0 = Position{ -intermedia, 0, intermedia };
            m_molecule.addBorderPoint(x0);
            m_molecule.addBorderPoint(y0);
            m_molecule.addBorderPoint(z0);
            x0 = Position{ intermedia, -intermedia, 0 };
            y0 = Position{ 0, intermedia, -intermedia };
            z0 = Position{ intermedia, 0, -intermedia };
            m_molecule.addBorderPoint(x0);
            m_molecule.addBorderPoint(y0);
            m_molecule.addBorderPoint(z0);
            intermedia = 1 / sqrt(3.0) * m_wall_spheric_radius;

            x0 = Position{ intermedia, intermedia, intermedia };
            m_molecule.addBorderPoint(x0);
            x0 = Position{ -intermedia, intermedia, intermedia };
            m_molecule.addBorderPoint(x0);
            x0 = Position{ intermedia, -intermedia, intermedia };
            m_molecule.addBorderPoint(x0);
            x0 = Position{ intermedia, intermedia, -intermedia };
            m_molecule.addBorderPoint(x0);
            x0 = Position{ -intermedia, intermedia, -intermedia };
            m_molecule.addBorderPoint(x0);
            x0 = Position{ intermedia, -intermedia, -intermedia };
            m_molecule.addBorderPoint(x0);
            x0 = Position{ -intermedia, -intermedia, intermedia };
            m_molecule.addBorderPoint(x0);
            x0 = Position{ -intermedia, -intermedia, -intermedia };
            m_molecule.addBorderPoint(x0);
        } else if (m_wall_type == 2) {
            Position x0 = Position{ m_wall_x_min, 0, 0 };
            Position x1 = Position{ m_wall_x_max, 0, 0 };
            Position y0 = Position{ 0, m_wall_y_min, 0 };
            Position y1 = Position{ 0, m_wall_y_max, 0 };
            Position z0 = Position{ 0, 0, m_wall_z_min };
            Position z1 = Position{ 0, 0, m_wall_z_max };

            m_molecule.addBorderPoint(x0);
            m_molecule.addBorderPoint(x1);
            m_molecule.addBorderPoint(y0);
            m_molecule.addBorderPoint(y1);
            m_molecule.addBorderPoint(z0);
            m_molecule.addBorderPoint(z1);

            x0 = Position{ m_wall_x_min, m_wall_y_min, 0 };
            x1 = Position{ m_wall_x_max, m_wall_y_max, 0 };
            y0 = Position{ m_wall_x_min, 0, m_wall_z_min };
            y1 = Position{ m_wall_x_max, 0, m_wall_z_min };
            z0 = Position{ 0, m_wall_y_min, m_wall_z_min };
            z1 = Position{ 0, m_wall_y_max, m_wall_z_max };

            m_molecule.addBorderPoint(x0);
            m_molecule.addBorderPoint(x1);
            m_molecule.addBorderPoint(y0);
            m_molecule.addBorderPoint(y1);
            m_molecule.addBorderPoint(z0);
            m_molecule.addBorderPoint(z1);

            x0 = Position{ m_wall_x_min, m_wall_y_min, m_wall_z_min };
            m_molecule.addBorderPoint(x0);
        }
    }
    // Claude Generated: Store wall configuration info for PrintStatus() display
    // Claude Generated 2025: Wall Configuration
    m_wall_geometry = m_config.get<std::string>("wall_type");
    m_wall_potential_type = m_config.get<std::string>("wall_potential");
    m_wall_auto_configured = auto_configured;

    // Calculate molecular density within wall boundaries
    if (m_wall_geometry == "rect" && (m_wall_x_max > m_wall_x_min) && (m_wall_y_max > m_wall_y_min) && (m_wall_z_max > m_wall_z_min)) {
        double volume = (m_wall_x_max - m_wall_x_min) * (m_wall_y_max - m_wall_y_min) * (m_wall_z_max - m_wall_z_min);
        m_molecular_density = 1.0 / volume; // molecules per Å³
    } else if (m_wall_geometry == "spheric" && m_wall_spheric_radius > 0) {
        double volume = (4.0 / 3.0) * pi * m_wall_spheric_radius * m_wall_spheric_radius * m_wall_spheric_radius;
        m_molecular_density = 1.0 / volume; // molecules per Å³
    }
    // Claude Generated: Wall configuration summary in PrintStatus() - show once every 10000 steps
    if (m_wall_geometry != "none" && m_wall_geometry != "") {
        std::cout << "\n--- Wall Setup ---\n";
        std::cout << "Geometry: " << m_wall_geometry << " | Potential: " << m_wall_potential_type;
        if (m_wall_auto_configured)
            std::cout << " (auto-sized)";
        std::cout << "\n";

        if (m_wall_geometry == "spheric") {
            std::cout << "Radius: " << m_wall_spheric_radius << " Å";
            if (m_molecular_density > 0) {
                std::cout << " | Density: " << m_molecular_density * 1e3 << " molecules/nm³";
            }
        } else if (m_wall_geometry == "rect") {
            double volume = (m_wall_x_max - m_wall_x_min) * (m_wall_y_max - m_wall_y_min) * (m_wall_z_max - m_wall_z_min);
            std::cout << "Bounds: [" << m_wall_x_min << "," << m_wall_x_max << "] x ["
                      << m_wall_y_min << "," << m_wall_y_max << "] x ["
                      << m_wall_z_min << "," << m_wall_z_max << "] Å";
            std::cout << " | Vol: " << volume << " Å³";
            if (m_molecular_density > 0) {
                std::cout << " | Density: " << m_molecular_density * 1e3 << " molecules/nm³";
            }
        }

        if (m_wall_violation_count > 0) {
            std::cout << " | Violations: " << m_wall_violation_count << "/" << m_natoms << " atoms";
        }
        std::cout << "\n---------------------------------\n";
    }
}

nlohmann::json SimpleMD::WriteRestartInformation()
{
    nlohmann::json restart;
    restart["method"] = m_method;
    restart["thermostat"] = m_thermostat;
    restart["dT"] = m_dT;
    restart["MaxTime"] = m_maxtime;
    restart["T"] = m_T0;
    restart["currentStep"] = m_currentStep;
    restart["seed"] = m_seed;
    restart["velocities"] = Tools::Geometry2String(m_eigen_velocities);
    restart["geometry"] = Tools::Geometry2String(m_eigen_geometry);
    restart["gradient"] = Tools::Geometry2String(m_eigen_gradient);
    restart["rmrottrans"] = m_rmrottrans;
    restart["nocenter"] = m_nocenter;
    restart["COM"] = m_COM;
    restart["average_T"] = m_aver_Temp;
    restart["average_Epot"] = m_aver_Epot;
    restart["average_Ekin"] = m_aver_Ekin;
    restart["average_Etot"] = m_aver_Etot;
    restart["average_Virial"] = m_average_virial_correction;
    restart["average_Wall"] = m_average_wall_potential;

    restart["rattle"] = m_rattle;
    restart["rattle_maxiter"] = m_rattle_maxiter;
    // restart["rattle_dynamic_tol"] = m_rattle_toler;
    restart["rattle_dynamic_tol_iter"] = m_rattle_dynamic_tol_iter;

    restart["coupling"] = m_coupling;
    restart["MaxTopoDiff"] = m_max_top_diff;
    restart["impuls"] = m_impuls;
    restart["impuls_scaling"] = m_impuls_scaling;
    restart["respa"] = m_respa;
    restart["rm_COM"] = m_rm_COM;
    restart["mtd"] = m_mtd;
    restart["rmsd_mtd"] = m_rmsd_mtd;
    restart["chainlength"] = m_chain_length;
    restart["eta"] = m_eta;
    restart["xi"] = Tools::DoubleVector2String(m_xi);
    restart["Q"] = Tools::DoubleVector2String(m_Q);

    if (m_rmsd_mtd) {
        restart["k_rmsd"] = m_k_rmsd;
        restart["alpha_rmsd"] = m_alpha_rmsd;
        restart["mtd_steps"] = m_mtd_steps;
        restart["rmsd_econv"] = m_rmsd_econv;
        restart["wtmtd"] = m_wtmtd;
        restart["rmsd_DT"] = m_rmsd_DT;
        restart["rmsd_ref_file"] = Basename() + ".mtd.xyz";
        restart["counter"] = m_bias_structure_count;
        restart["rmsd_atoms"] = m_rmsd_atoms;
        std::vector<json> bias(m_bias_structure_count);
        for (const auto & m_bias_thread : m_bias_threads) {
            for (const auto& stored_bias : m_bias_thread->getBias()) {
                bias[stored_bias["index"]] = stored_bias;
            }
        }
        json bias_restart;
        for (int i = 0; i < bias.size(); ++i) {
            bias_restart[i] = bias[i];
        }
        restart["bias"] = bias_restart;
    }
    if (m_rattle) {
        json constrains;
        if (m_rattle_12) {
            json constrains_12;
            for (int i = 0; i < m_bond_constrained.size(); ++i) {
                json element;
                element["i"] = m_bond_constrained[i].first.first;
                element["j"] = m_bond_constrained[i].first.second;
                element["d"] = m_bond_constrained[i].second;

                constrains_12[i] = element;
            }
            constrains["constrain_12"] = true;
            constrains["num_constrain_12"] = m_bond_constrained.size();
            constrains["constrains_12"] = constrains_12;
        }

        if (m_rattle_13) {
            json constrains_13;
            for (int i = 0; i < m_bond_13_constrained.size(); ++i) {
                json element;
                element["i"] = m_bond_13_constrained[i].first.first;
                element["j"] = m_bond_13_constrained[i].first.second;
                element["d"] = m_bond_13_constrained[i].second;

                constrains_13[i] = element;
            }
            constrains["constrain_13"] = true;
            constrains["num_constrain_13"] = m_bond_13_constrained.size();
            constrains["constrains_13"] = constrains_13;
        }
        restart["constrains"] = constrains;
    }
    return restart;
};

bool SimpleMD::LoadRestartInformation()
{
    if (!Restart())
        return false;
    StringList files = RestartFiles();
    int error = 0;
    for (const auto& f : files) {
        std::ifstream file(f);
        json restart;
        try {
            file >> restart;
        } catch ([[maybe_unused]] json::type_error& e) {
            error++;
            continue;
        } catch ([[maybe_unused]] json::parse_error& e) {
            error++;
            continue;
        }

        json md;
        try {
            md = restart[MethodName()[0]];
        } catch ([[maybe_unused]] json::type_error& e) {
            error++;
            continue;
        }
        return LoadRestartInformation(md);
    }
    return true;
};

bool SimpleMD::LoadRestartInformation(const json& state)
{
    // Claude Generated 2025: Validate restart data before parsing
    auto validation = validateRestartData(state,
        {"method", "geometry", "velocities"},  // required fields
        {"geometry", "velocities", "xi", "Q"}); // fields to validate for doubles

    if (!validation.valid) {
        std::cerr << "\033[1;31m[ERROR]\033[0m Restart file validation failed: "
                  << validation.error_message << std::endl;
        std::cerr << "\033[1;33m[WARNING]\033[0m Starting fresh simulation instead." << std::endl;
        return false;
    }

    std::string geometry, velocities, constrains, xi, Q;

    try {
        m_method = state["method"];
    } catch ([[maybe_unused]] json::type_error& e) {
    }

    try {
        m_dT = state["dT"];
    } catch ([[maybe_unused]] json::type_error& e) {
    }
    try {
        m_maxtime = state["MaxTime"];
    } catch ([[maybe_unused]] json::type_error& e) {
    }
    try {
        m_rmrottrans = state["rmrottrans"];
    } catch ([[maybe_unused]] json::type_error& e) {
    }
    try {
        m_nocenter = state["nocenter"];
    } catch ([[maybe_unused]] json::type_error& e) {
    }
    try {
        m_COM = state["COM"];
    } catch ([[maybe_unused]] json::type_error& e) {
    }
    try {
        m_T0 = state["T"];
    } catch ([[maybe_unused]] json::type_error& e) {
    }
    try {
        m_currentStep = state["currentStep"];
    } catch ([[maybe_unused]] json::type_error& e) {
    }

    try {
        m_aver_Epot = state["average_Epot"];
    } catch ([[maybe_unused]] json::type_error& e) {
    }

    try {
        m_aver_Ekin = state["average_Ekin"];
    } catch ([[maybe_unused]] json::type_error& e) {
    }

    try {
        m_aver_Etot = state["average_Etot"];
    } catch ([[maybe_unused]] json::type_error& e) {
    }
    try {
        m_aver_Temp = state["average_T"];
    } catch ([[maybe_unused]] json::type_error& e) {
    }

    try {
        m_average_virial_correction = state["average_Virial"];
    } catch ([[maybe_unused]] json::type_error& e) {
    }

    try {
        m_average_wall_potential = state["average_Wall"];
    } catch ([[maybe_unused]] json::type_error& e) {
    }

    try {
        m_coupling = state["coupling"];
    } catch ([[maybe_unused]] json::type_error& e) {
    }

    try {
        m_respa = state["respa"];
    } catch ([[maybe_unused]] json::type_error& e) {
    }

    try {
        m_eta = state["eta"];
    } catch ([[maybe_unused]] json::type_error& e) {
    }

    try {
        m_thermostat = state["thermostat"];
    } catch ([[maybe_unused]] json::type_error& e) {
    }

    try {
        geometry = state["geometry"];
    } catch ([[maybe_unused]] json::type_error& e) {
    }

    try {
        velocities = state["velocities"];
    } catch ([[maybe_unused]] json::type_error& e) {
    }

    try {
        xi = state["xi"];
    } catch ([[maybe_unused]] json::type_error& e) {
    }
    try {
        Q = state["Q"];
    } catch ([[maybe_unused]] json::type_error& e) {
    }

    try {
        m_mtd = state["mtd"];
    } catch ([[maybe_unused]] json::type_error& e) {
    }

    try {
        m_rattle = state["rattle"];
    } catch ([[maybe_unused]] json::type_error& e) {
    }
    if(m_rattle)
    {
        try {
            m_rattle_tol_12 = state["rattle_tol_12"];
        } catch ([[maybe_unused]] json::type_error& e) {
        }

        try {
            m_rattle_tol_13 = state["rattle_tol_13"];
        } catch (json::type_error& e) {
        }
        try {
            m_rattle_maxiter = state["rattle_maxiter"];
        } catch ([[maybe_unused]] json::type_error& e) {
        }

        try {
            m_rattle_dynamic_tol = state["rattle_dynamic_tol"];
        } catch ([[maybe_unused]] json::type_error& e) {
        }

        try {
            m_rattle_dynamic_tol_iter = state["rattle_dynamic_tol_iter"];
        } catch ([[maybe_unused]] json::type_error& e) {
        }
    }
    try {
        m_seed = state["seed"];
    } catch (json::type_error& e) {
    }
    try {
        m_rmsd_mtd = state["rmsd_mtd"];
        if (m_rmsd_mtd) {
            m_k_rmsd = state["k_rmsd"];
            m_alpha_rmsd = state["alpha_rmsd"];
            m_mtd_steps = state["mtd_steps"];
            m_rmsd_econv = state["rmsd_econv"];
            m_wtmtd = state["wtmtd"];
            m_rmsd_DT = state["rmsd_DT"];
            m_rmsd_ref_file = state["rmsd_ref_file"];
            m_bias_json = state["bias"];
        }
    } catch ([[maybe_unused]] json::type_error& e) {
    }

    // Claude Generated 2025: Robust parsing of restart data with error handling
    if (!geometry.empty()) {
        try {
            Tools::String2Geometry(m_eigen_geometry, geometry);
        } catch (const std::invalid_argument& e) {
            std::cerr << "\033[1;33m[WARNING]\033[0m Malformed geometry data in restart file (e.g., truncated NaN values). Ignoring restart geometry." << std::endl;
            return false;
        } catch (const std::out_of_range& e) {
            std::cerr << "\033[1;33m[WARNING]\033[0m Out-of-range geometry values in restart file. Ignoring restart geometry." << std::endl;
            return false;
        }
    }

    if (!velocities.empty()) {
        try {
            Tools::String2Geometry(m_eigen_velocities, velocities);
        } catch (const std::invalid_argument& e) {
            std::cerr << "\033[1;33m[WARNING]\033[0m Malformed velocity data in restart file (e.g., truncated NaN values). Ignoring restart velocities." << std::endl;
            return false;
        } catch (const std::out_of_range& e) {
            std::cerr << "\033[1;33m[WARNING]\033[0m Out-of-range velocity values in restart file. Ignoring restart velocities." << std::endl;
            return false;
        }
    }

    if (!xi.empty()) {
        try {
            m_xi = Tools::String2DoubleVec(xi, "|");
        } catch (const std::invalid_argument& e) {
            std::cerr << "\033[1;33m[WARNING]\033[0m Malformed thermostat xi data in restart file. Ignoring restart xi." << std::endl;
            return false;
        }
    }

    if (!Q.empty()) {
        try {
            m_Q = Tools::String2DoubleVec(Q, "|");
        } catch (const std::invalid_argument& e) {
            std::cerr << "\033[1;33m[WARNING]\033[0m Malformed thermostat Q data in restart file. Ignoring restart Q." << std::endl;
            return false;
        }
    }

    try {
        if (state.contains("constrains")) {
            const json& constrains = state["constrains"];

            if (constrains.contains("constrain_12") && constrains["constrain_12"].get<bool>()) {
                m_bond_constrained.clear();
                int num_constrain_12 = constrains["num_constrain_12"];
                const json& constrains_12 = constrains["constrains_12"];

                for (int i = 0; i < num_constrain_12; ++i) {
                    int i_index = constrains_12[i]["i"].get<int>();
                    int j_index = constrains_12[i]["j"].get<int>();
                    double distance = constrains_12[i]["d"].get<double>();
                    m_bond_constrained.emplace_back(std::make_pair(i_index, j_index), distance);
                    std::cout << "1,2: " << i_index << " " << j_index << " " << distance << " ";
                }
            }

            if (constrains.contains("constrain_13") && constrains["constrain_13"].get<bool>()) {
                m_bond_13_constrained.clear();
                int num_constrain_13 = constrains["num_constrain_13"];
                const json& constrains_13 = constrains["constrains_13"];

                for (int i = 0; i < num_constrain_13; ++i) {
                    int i_index = constrains_13[i]["i"].get<int>();
                    int j_index = constrains_13[i]["j"].get<int>();
                    double distance = constrains_13[i]["d"].get<double>();
                    m_bond_13_constrained.emplace_back(std::make_pair(i_index, j_index), distance);
                    std::cout << "1,3: " << i_index << " " << j_index << " " << distance << " ";
                }
            }
        }
    } catch ([[maybe_unused]] json::type_error& e) {
        std::cerr << e.what() << '\n';
    }

    m_restart = !geometry.empty() && !velocities.empty();

    return true;
}

// Claude Generated 2026 - Stepwise MD refactor (Qurcuma interactive simulation)
// start() previously contained setup + loop + cleanup in one monolithic block.
// It is now split into prepareRun() / step() / finalizeRun() so external
// callers (e.g. the qurcuma GUI worker) can drive the MD integration one step
// at a time and inject per-step perturbations between steps.
void SimpleMD::start()
{
    prepareRun();
    if (!m_run_prepared)
        return;
    while (step()) {
        // step() returns false when the simulation should terminate
    }
    finalizeRun();
}

void SimpleMD::prepareRun()
{
    if (m_initialised == false) {
        m_run_prepared = false;
        return;
    }
    m_run_aborted = false;
    m_run_states.clear();

    auto unix_timestamp = std::chrono::seconds(std::time(nullptr));
    m_unix_started = std::chrono::milliseconds(unix_timestamp).count();

    // Claude Generated 2025: Thermostat selection - Enum-based switch
    auto thermo_it = thermostat_map.find(m_thermostat);
    ThermostatType thermo = (thermo_it != thermostat_map.end()) ? thermo_it->second : ThermostatType::CSVR;

    switch (thermo) {
        case ThermostatType::CSVR:
            CitationRegistry::cite("csvr");
            ThermostatFunction = [this] { CSVR(); };
            break;
        case ThermostatType::Berendsen:
            CitationRegistry::cite("berendsen");
            ThermostatFunction = [this] { Berendson(); };
            break;
        case ThermostatType::Anderson:
            fmt::print(fg(fmt::color::green) | fmt::emphasis::bold, "\nUsing Anderson Thermostat\n ... \n\n");
            ThermostatFunction = [this] { Anderson(); };
            break;
        case ThermostatType::NoseHover:
            fmt::print(fg(fmt::color::green) | fmt::emphasis::bold, "\nUsing Nosé-Hoover-Chain Thermostat\n ... \n\n");
            ThermostatFunction = [this] { NoseHover(); };
            break;
        case ThermostatType::None:
            // No thermostat
            break;
    }

    if (thermo == ThermostatType::None) {
        ThermostatFunction = [this] { None(); };
        std::cout << "No Thermostat applied\n"
                  << std::endl;
    }

    m_Epot = Energy();
    EKin();
    m_Etot = m_Epot + m_Ekin;
    AverageQuantities();
    m_step = 0;

    // Claude Generated (Jun 2026): reference state for the opt-in robustness gates.
    // epot_ref is the bare starting potential; the topology check interval defaults to dump.
    m_epot_ref = m_Epot;
    m_topo_check_every = (m_topo_check_interval > 0) ? m_topo_check_interval : m_dump;

    // Claude Generated (Jun 2026): when freezing inherited hill heights, record the counters of
    // every bias structure already in the shared pool at this run's start. Those structures then
    // contribute at a fixed height for the whole run (and are not bumped by it), so only this run's
    // own deposits grow -> the cumulative bias force no longer escalates run after run.
    m_frozen_height.clear();
    if (m_shared_pool && m_freeze_inherited) {
        for (const auto& bs : m_shared_pool->snapshot())
            m_frozen_height[bs.index] = bs.counter;
    }

    // WP-S2 (May 2026): open per-step diagnostics JSONL file
    if (m_md_diagnostics) {
        m_diag_writer = std::make_unique<MDDiagnosticsWriter>(Basename() + ".diag.jsonl");
        if (!m_diag_writer->isOpen()) {
            CurcumaLogger::warn("MD diagnostics requested but JSONL file could not be opened — disabled");
            m_md_diagnostics = false;
            m_diag_writer.reset();
        }
    }

    // WP-P1 (May 2026): force per-phase timing collection when diagnostics-timing is on
    if (m_md_diagnostics && m_md_diagnostics_timing && m_interface) {
        m_interface->setForcePhaseTiming(true);
    }

    WriteGeometry();
#ifdef USE_Plumed
    if (m_mtd) {
        m_plumedmain = plumed_create();
        int real_precision = 8;
        double energyUnits = 2625.5;
        double lengthUnits = 10;
        double timeUnits = 1e-3;
        double massUnits = 1;
        double chargeUnit = 1;
        int restart = m_restart;
        plumed_cmd(m_plumedmain, "setRealPrecision", &real_precision); // Pass a pointer to an integer containing the size of a real number (4 or 8)
        plumed_cmd(m_plumedmain, "setMDEnergyUnits", &energyUnits); // Pass a pointer to the conversion factor between the energy unit used in your code and kJ mol-1
        plumed_cmd(m_plumedmain, "setMDLengthUnits", &lengthUnits); // Pass a pointer to the conversion factor between the length unit used in your code and nm
        plumed_cmd(m_plumedmain, "setMDTimeUnits", &timeUnits); // Pass a pointer to the conversion factor between the time unit used in your code and ps
        plumed_cmd(m_plumedmain, "setNatoms", &m_natoms); // Pass a pointer to the number of atoms in the system to plumed
        plumed_cmd(m_plumedmain, "setMDEngine", "curcuma");
        plumed_cmd(m_plumedmain, "setMDMassUnits", &massUnits); // Pass a pointer to the conversion factor between the mass unit used in your code and amu
        plumed_cmd(m_plumedmain, "setMDChargeUnits", &chargeUnit);
        plumed_cmd(m_plumedmain, "setTimestep", &m_dT); // Pass a pointer to the molecular dynamics timestep to plumed                       // Pass the name of your md engine to plumed (now it is just a label)
        plumed_cmd(m_plumedmain, "setKbT", &kb_Eh);
        plumed_cmd(m_plumedmain, "setLogFile", "plumed_log.out"); // Pass the file  on which to write out the plumed log (to be created)
        plumed_cmd(m_plumedmain, "setRestart", &restart); // Pointer to an integer saying if we are restarting (zero means no, one means yes)
        plumed_cmd(m_plumedmain, "init", NULL);
        plumed_cmd(m_plumedmain, "read", m_plumed.c_str());
        plumed_cmd(m_plumedmain, "setStep", &m_step);
        plumed_cmd(m_plumedmain, "setPositions", &m_eigen_geometry.data()[0]);
        plumed_cmd(m_plumedmain, "setEnergy", &m_Epot);
        plumed_cmd(m_plumedmain, "setForces", &m_eigen_gradient.data()[0]);
        plumed_cmd(m_plumedmain, "setVirial", &m_virial[0]);
        plumed_cmd(m_plumedmain, "setMasses", &m_eigen_masses.data()[0]);
        plumed_cmd(m_plumedmain, "prepareCalc", NULL);
        plumed_cmd(m_plumedmain, "performCalc", NULL);
    }
#endif
    std::vector<double> charge(0, m_natoms);

    // Claude Generated (May 2026): unified MD table header — matches the data rows
    // printed below. Drops the #ifdef GCC branches (the GCC-only fmt::format path was
    // dead because no compiler defines `GCC`; only `__GNUC__`) and the 5-column
    // tab-separated fallback that didn't line up with the 15+ column data.
    {
        std::string header = fmt::format(
            "{1: ^{0}} {2: ^{0}} {3: ^{0}} {4: ^{0}} {5: ^{0}} {6: ^{0}} {7: ^{0}} "
            "{8: ^{0}} {9: ^{0}} {10: ^{0}} {11: ^{0}} {12: ^{0}} {13: ^{0}} {14: ^{0}} "
            "{15: ^{0}}",
            15,
            "Time", "Epot", "<Epot>", "Ekin", "<Ekin>", "Etot", "<Etot>",
            "T", "<T>", "Wall", "<Wall>", "Virial", "<Virial>", "Remaining", "dt");
        std::string units = fmt::format(
            "{1: ^{0}} {2: ^{0}} {3: ^{0}} {4: ^{0}} {5: ^{0}} {6: ^{0}} {7: ^{0}} "
            "{8: ^{0}} {9: ^{0}} {10: ^{0}} {11: ^{0}} {12: ^{0}} {13: ^{0}} {14: ^{0}} "
            "{15: ^{0}}",
            15,
            "ps", "Eh", "Eh", "Eh", "Eh", "Eh", "Eh",
            "K", "K", "Eh", "Eh", "Eh", "Eh", "s", "ps");
        if (m_dipole) {
            header += fmt::format(" {: ^15}", "Dipole");
            units  += fmt::format(" {: ^15}", "Debye");
        }
        if (m_rmsd_mtd) {
            header += fmt::format(" {: ^15}", "nBias");
            units  += fmt::format(" {: ^15}", "#");
        } else if (m_writeUnique) {
            header += fmt::format(" {: ^15}", "nUnique");
            units  += fmt::format(" {: ^15}", "#");
        }
        std::cout << header << "\n" << units << "\n";
    }
    if (m_rmsd_mtd) {
        CurcumaLogger::result_fmt("RMSD-MTD: k={} Eh, alpha={} Bohr^-2, pace={} steps",
            m_k_rmsd, m_alpha_rmsd, m_mtd_steps);
        CurcumaLogger::result_fmt("RMSD-MTD: Econv={}, max_gaussians={}",
            m_rmsd_econv, m_max_rmsd_N);
        if (m_wtmtd)
            CurcumaLogger::result_fmt("RMSD-MTD: Well-tempered (dT={})", m_rmsd_DT);
        else
            CurcumaLogger::result("RMSD-MTD: Well-tempered Off");
        if (m_shared_pool)
            CurcumaLogger::result_fmt("RMSD-MTD: Shared bias pool active ({} structures)",
                m_shared_pool->biasStructureCount());
    }
    PrintStatus();
    m_run_prepared = true;
}

/* Claude Generated 2026 - One iteration of the MD loop.
 * Returns true while the simulation should continue, false when it should end.
 * Termination reasons: max_time reached, CheckStop() (stop file), unstable
 * dynamics, or exhausted rescue budget. The step() body mirrors the original
 * while-loop body in start() verbatim. External-force injection happens at the
 * very top so the contribution is folded into the upcoming integrator step. */
bool SimpleMD::step()
{
    if (!m_run_prepared)
        return false;

    if (!(m_maxtime <= 0 || m_currentStep < m_maxtime))
        return false;

    auto step0 = std::chrono::system_clock::now();

    if (CheckStop() == true) {
        TriggerWriteRestart();
        m_run_aborted = true;
#ifdef USE_Plumed
        if (m_mtd) {
            plumed_finalize(m_plumedmain);
        }
#endif
        return false;
    }

    // Claude Generated 2026 - Inject queued external forces before integration.
    // The contribution is consumed (cleared) in a single step; callers must
    // re-apply each step while the user is actively dragging an atom.
    if (m_external_forces_pending) {
        m_eigen_gradient += m_external_forces;
        m_external_forces.setZero();
        m_external_forces_pending = false;
    }

    if (m_rm_COM_step > 0 && m_step % m_rm_COM_step == 0) {
        if (m_rmrottrans == 1)
            RemoveRotation();
        else if (m_rmrottrans == 2)
            RemoveRotations();
        else if (m_rmrottrans == 3) {
            RemoveRotations();
            RemoveRotation();
        }
    }

    // WP-P1 (May 2026): wall-clock the integrator step for MD diagnostics
    double integrator_ms = 0.0;
    if (m_md_diagnostics_timing) {
        auto t0 = std::chrono::high_resolution_clock::now();
        Integrator();
        auto t1 = std::chrono::high_resolution_clock::now();
        integrator_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    } else {
        Integrator();
    }
    m_last_integrator_ms = integrator_ms;
    AverageQuantities();

    // Claude Generated (Jun 2026): ConfSearch robustness gates (opt-in, default off).
    // Both set m_run_aborted (like the stop-file path) so finalizeRun() still runs and
    // the bias pool / final frame are finalised cleanly.
    //
    // Early Epot abort: a sustained climb of the running-mean BARE potential (bias is
    // separate from m_Epot) means the walker left the relevant low-energy region. The
    // window must exceed the thermal baseline (~N_dof*kT/2) plus typical barriers, so a
    // single barrier crossing (transient spike) does not trip it. Warm up first so the
    // cumulative mean is not dominated by the initial transient.
    if (m_epot_abort && m_step > 100
        && (m_aver_Epot - m_epot_ref) * 2625.5 > m_epot_abort_window) {
        // fmt::print (not CurcumaLogger): abort diagnostics must stay visible even when the
        // global logger verbosity is clamped to 0 by the energy-setup path (see roadmap issue #3),
        // mirroring the "Simulation got unstable" message below.
        fmt::print(fg(fmt::color::orange) | fmt::emphasis::bold,
            "MD aborted: <Epot> climbed {:.1f} kJ/mol above start (window {:.1f})\n",
            (m_aver_Epot - m_epot_ref) * 2625.5, m_epot_abort_window);
        m_run_aborted = true;
        return false;
    }
    // Temperature runaway abort: the shared bias pool's hills (W_i = k*counter_i) grow over
    // successive runs and pump energy in faster than the thermostat removes it -> the running-mean
    // temperature climbs above the target. Either an over-factor (relative) or over-delta (absolute)
    // threshold trips it; a threshold with value <= 0 is disabled. Warm up first so the cumulative
    // mean is not dominated by the initial transient. Claude Generated (Jun 2026).
    if (m_temp_abort && m_step > 100) {
        const bool over_factor = (m_temp_abort_factor > 0 && m_aver_Temp > m_temp_abort_factor * m_T0);
        const bool over_delta = (m_temp_abort_delta > 0 && m_aver_Temp > m_T0 + m_temp_abort_delta);
        if (over_factor || over_delta) {
            // fmt::print: stay visible despite the verbosity clamp (see epot_abort note above).
            fmt::print(fg(fmt::color::orange) | fmt::emphasis::bold,
                "MD aborted: <T>={:.0f} K ran away from target {:.0f} K (factor limit {}x, delta limit {} K)\n",
                m_aver_Temp, m_T0, m_temp_abort_factor, m_temp_abort_delta);
            m_run_aborted = true;
            return false;
        }
    }
    // Topology abort: a growth in the connected-component count means the molecule
    // fragmented. Only the fragment count is used (robust against transient bond-length
    // fluctuations at high T). m_eigen_geometry holds the current coordinates.
    if (m_topo_check && m_topo_check_every > 0 && m_step > 0 && m_step % m_topo_check_every == 0) {
        m_molecule.setGeometry(m_eigen_geometry);
        int nfrag = static_cast<int>(m_molecule.GetFragments().size());
        if (nfrag > m_start_fragment_count) {
            // fmt::print: stay visible despite the verbosity clamp (see epot_abort note above).
            fmt::print(fg(fmt::color::orange) | fmt::emphasis::bold,
                "MD aborted: topology broke (fragments {} -> {})\n",
                m_start_fragment_count, nfrag);
            m_run_aborted = true;
            return false;
        }
    }

    if (m_mtd) {
        if (!m_eval_mtd) {
            if (std::abs(m_T0 - m_aver_Temp) < m_mtd_dT && m_step > 10) {
                m_eval_mtd = true;
                std::cout << "Starting with MetaDynamics ..." << std::endl;
            }
        }
    }

    /////////// Dipole
    if (m_dipole && m_method == "gfn2") {
        //linear Dipoles
        auto curr_dipoles_lin = m_molecule.CalculateDipoleMoments(m_scaling_vector_linear, m_start_fragments);
        std::ofstream file;
        file.open(outputPath(Basename() + "_dipole_linear.out"), std::ios_base::app);
        Position d = {0,0,0};
        for (const auto& dipole_lin : curr_dipoles_lin) {
            d += dipole_lin;
            file << dipole_lin[0] << " " << dipole_lin[1] << " " << dipole_lin[2] << " " << dipole_lin.norm() << ", ";
        }
        file << d[0] << " " << d[1] << " " << d[2] << ", " << m_molecule.getDipole()[0] << " " << m_molecule.getDipole()[1] << " " << m_molecule.getDipole()[2] << std::endl;
        file.close();
        //nonlinear Dipoles
        auto curr_dipoles_nlin = m_molecule.CalculateDipoleMoments(m_scaling_vector_nonlinear, m_start_fragments);
        std::ofstream file2;
        file2.open(outputPath(Basename() + "_dipole_nonlinear.out"), std::ios_base::app);
        Position sum = {0,0,0};
        for (const auto& dipole_nlin : curr_dipoles_nlin) {
            sum += dipole_nlin;
            file2 << dipole_nlin[0] << " " << dipole_nlin[1] << " " << dipole_nlin[2] << " " << dipole_nlin.norm() <<", ";
        }
        file2 << sum[0] << " " << sum[1] << " " << sum[2] << ", " << m_molecule.getDipole()[0] << " " << m_molecule.getDipole()[1] << " " << m_molecule.getDipole()[2] << std::endl;
        file2.close();
    }
    //////////// Dipole

    if (m_step % m_dump == 0) {
        if (bool write = WriteGeometry()) {
            m_run_states.push_back(WriteRestartInformation());
            m_current_rescue = 0;
            // WP-S2 (May 2026): append diagnostics record parallel to XYZ trajectory
            if (m_md_diagnostics && m_diag_writer) {
                // WP-P1 (May 2026): optional per-phase wall-clock breakdown
                json timing;
                if (m_md_diagnostics_timing) {
                    timing = m_interface->LastPrepTiming();
                    timing["ff_total"]    = m_last_ff_ms;
                    timing["integrator"]  = m_last_integrator_ms;
                    timing["hbxb_update"] = m_last_hbxb_ms;
                    auto t_now = std::chrono::system_clock::now();
                    timing["step_total"] = std::chrono::duration<double, std::milli>(t_now - step0).count();
                    json gpu_t = m_interface->StreamTimings();
                    if (!gpu_t.empty()) {
                        timing["gpu"] = gpu_t;
                    }
                }
                m_diag_writer->writeSnapshot(
                    m_step, m_currentStep,
                    m_interface->getEnergyDecomposition(),
                    m_interface->Charges(),
                    m_interface->CN(),
                    m_interface->Gradient(),
                    m_interface->HBCount(),
                    m_interface->XBCount(),
                    timing);
            }
        } else if (!write && m_rescue && m_run_states.size() > (1 - m_current_rescue)) {
            std::cout << "Molecule exploded, resetting to previous state ..." << std::endl;
            LoadRestartInformation(m_run_states[m_run_states.size() - 1 - m_current_rescue]);
            Geometry geometry = m_molecule.getGeometry();
            for (int i = 0; i < m_natoms; ++i) {
                geometry(i, 0) = m_eigen_geometry.data()[3 * i + 0];
                geometry(i, 1) = m_eigen_geometry.data()[3 * i + 1];
                geometry(i, 2) = m_eigen_geometry.data()[3 * i + 2];
            }
            m_molecule.setGeometry(geometry);
            m_molecule.GetFragments();
            InitVelocities(-1);
            Energy();
            EKin();
            m_Etot = m_Epot + m_Ekin;
            m_current_rescue++;
            PrintStatus();
            m_time_step = 0;
        }
    }

    if (m_unstable || m_interface->Error() || m_interface->HasNan()) {
        PrintStatus();
        fmt::print(fg(fmt::color::salmon) | fmt::emphasis::bold, "Simulation got unstable, exiting!\n");

        // Per-instance filename (Basename() carries the ConfSearch ".t<id>" suffix) so concurrent
        // MD workers do not clobber each other's crash dump during simultaneous instability cleanup.
        std::ofstream restart_file(snapshotPath(Basename() + ".unstable.json"));
        nlohmann::json restart;
        restart[MethodName()[0]] = WriteRestartInformation();
        restart_file << restart << std::endl;

        m_time_step = 0;
        m_run_aborted = true;

#ifdef USE_Plumed
        if (m_mtd) {
            plumed_finalize(m_plumedmain);
        }
#endif
        m_run_prepared = false;  // signal to start(): skip finalizeRun() cleanup
        return false;
    }

    if (m_writerestart > -1 && m_step % m_writerestart == 0) {
        std::ofstream restart_file(snapshotPath(Basename() + "_step_" + std::to_string(static_cast<int>(m_step * m_dT)) + ".json"));
        json restart;
        restart[MethodName()[0]] = WriteRestartInformation();
        restart_file << restart << std::endl;
    }
    if ((m_step && static_cast<int>(m_step * m_dT) % m_print == 0)) {
        m_Etot = m_Epot + m_Ekin;
        PrintStatus();
        // Reset RATTLE diagnostics after printing
        m_rattle_max_err_12 = 0;
        m_rattle_max_err_13 = 0;
        m_rattle_max_err_count = 0;
        m_time_step = 0;
    }
    if (m_rattle && m_rattle_dynamic_tol) {
        m_aver_rattle_Temp += m_T;
        m_rattle_counter++;
        if (m_rattle_counter == m_rattle_dynamic_tol_iter)
            AdjustRattleTolerance();
    }
    // Temporarily disabled: impuls re-initialization overrides thermostat-controlled
    // temperature ramps in ConfSearch. Re-enable after testing temperature stability.
    /*
    if (m_impuls > m_T) {
        InitVelocities(m_scale_velo * m_impuls_scaling);
        EKin();
        m_time_step = 0;
    }
    */

    if (m_current_rescue >= m_max_rescue) {
        fmt::print(fg(fmt::color::salmon) | fmt::emphasis::bold, "Nothing really helps");
        return false;
    }
    m_step++;
    m_currentStep += m_dT;
    m_time_step += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - step0).count();
    return true;
}

/* Claude Generated 2026 - Post-loop cleanup. Writes final trajectory frame,
 * prints summary, finalises plumed/RMSD-metadynamics, writes curcuma_final.json.
 * Skipped by start() when step() triggered an early abort (unstable dynamics).
 * Always resets m_run_prepared so the object can be re-run. */
void SimpleMD::finalizeRun()
{
    if (!m_run_prepared) {
        // early-exit path (unstable) already did its own cleanup inside step()
        return;
    }

    // Claude Generated (October 2025): Write final frame to ensure at least 2 frames (t=0, t=max_time)
    // Fixes trajectory file generation for short simulations where dump_frequency > total_steps
    WriteGeometry();

    PrintStatus();
    if (m_thermostat == "csvr")
        std::cout << "Exchange with heat bath " << m_Ekin_exchange << "Eh" << std::endl;
    if (m_dipole) {
        std::cout << "Calculated averaged dipole moment " << m_aver_dipol_linear * 2.5418 << " Debye and " << m_aver_dipol_linear * 2.5418 * 3.3356 << " Cm [e-30]" << std::endl;
    }

#ifdef USE_Plumed
    if (m_mtd) {
        plumed_finalize(m_plumedmain);
    }
#endif
    if (m_rmsd_mtd) {
        std::cout << "Sum of Energy of COLVARs:" << std::endl;
        for (int i = 0; i < m_bias_threads.size(); ++i) {
            auto structures = m_bias_threads[i]->getBiasStructure();
            for (int j = 0; j < structures.size(); ++j) {
                std::cout << structures[j].rmsd_reference << "\t" << structures[j].energy << "\t" << structures[j].counter / static_cast<double>(m_colvar_incr) * 100 << std::endl;

                m_rmsd_mtd_molecule.setGeometry(structures[j].geometry);
                m_rmsd_mtd_molecule.setEnergy(structures[j].energy);
                m_rmsd_mtd_molecule.setName(std::to_string(structures[j].index) + " " + std::to_string(structures[j].rmsd_reference));
                if (i == j && i == 0)
                    m_rmsd_mtd_molecule.writeXYZFile(outputPath(Basename() + ".mtd.xyz"));
                else
                    m_rmsd_mtd_molecule.appendXYZFile(outputPath(Basename() + ".mtd.xyz"));
            }
        }
    }
    // Per-instance filename so concurrent MD workers don't overwrite each other's final dump.
    std::ofstream restart_file(snapshotPath(Basename() + ".final.json"));
    nlohmann::json restart;
    restart[MethodName()[0]] = WriteRestartInformation();
    restart_file << restart << std::endl;
    if (m_run_aborted == false)
        std::remove(snapshotPath("curcuma_restart.json").c_str());

    m_run_prepared = false;
}

/* Claude Generated 2026 - Queue external per-atom force contribution for the
 * next step() call. Forces are in Hartree/Bohr (same units as m_eigen_gradient).
 * The contribution is added additively and then cleared, so callers must re-queue
 * each step while a user drag is active. Shape check: forces must match (natoms, 3). */
void SimpleMD::applyExternalForces(const Geometry& forces)
{
    if (forces.rows() != m_natoms || forces.cols() != 3)
        return;

    if (m_external_forces.rows() != m_natoms || m_external_forces.cols() != 3)
        m_external_forces = Geometry::Zero(m_natoms, 3);

    m_external_forces = forces;
    m_external_forces_pending = true;
}

void SimpleMD::AdjustRattleTolerance()
{
    m_aver_rattle_Temp /= static_cast<double>(m_rattle_counter);

    // std::pair<double, double> pair(m_rattle_tolerance, m_aver_Temp);

    if (m_aver_rattle_Temp > m_T0)
        m_rattle_tol_12 -= 0.01;
    else if (m_aver_rattle_Temp < m_T0)
        m_rattle_tol_12 += 0.01;
    std::cout << m_rattle_counter << " " << m_aver_rattle_Temp << " " << m_rattle_tol_12 << std::endl;
    m_rattle_tol_12 = std::abs(m_rattle_tol_12);
    m_rattle_counter = 0;
    m_aver_rattle_Temp = 0;
}

// Claude Generated (Oct 2025): Apply Periodic Boundary Conditions
// Wraps all atoms into the central unit cell using fractional coordinates
void SimpleMD::applyPeriodicBoundaryConditions()
{
    if (!m_has_pbc) return;

    // Get unit cell matrix and its inverse
    Eigen::Matrix3d cell = m_molecule.getUnitCell();
    Eigen::Matrix3d cell_inv = m_molecule.getUnitCellInverse();

    // Wrap each atom into central cell
    for (int i = 0; i < m_natoms; ++i) {
        Eigen::Vector3d pos = m_eigen_geometry.row(i);

        // Convert to fractional coordinates: r_frac = cell_inv * r_cart
        Eigen::Vector3d frac = cell_inv * pos;

        // Wrap to [0, 1) range
        frac[0] -= std::floor(frac[0]);
        frac[1] -= std::floor(frac[1]);
        frac[2] -= std::floor(frac[2]);

        // Convert back to Cartesian: r_cart = cell * r_frac
        pos = cell * frac;
        m_eigen_geometry.row(i) = pos;
    }

    // Update molecule geometry after PBC wrapping
    m_molecule.setGeometry(m_eigen_geometry);
}

void SimpleMD::Verlet()
{
    // CRITICAL FIX (Feb 2026): Check gradient for NaN/Inf BEFORE integration
    // Aromatic systems (benzene) can produce gradient singularities at φ=180° or ω=0
    // Early detection prevents NaN propagation through velocity/position updates
    for (int i = 0; i < 3 * m_natoms; ++i) {
        if (!std::isfinite(m_eigen_gradient.data()[i])) {
            CurcumaLogger::error("NaN/Inf gradient at coordinate " + std::to_string(i) +
                               " (atom " + std::to_string(i/3) + ") before integration");
            m_unstable = true;
            return;
        }
    }

    double ekin = 0;
    //std::cout << m_eigen_inv_masses << std::endl;
    for (int i = 0; i < m_natoms; ++i) {
        m_eigen_geometry.data()[3 * i + 0] = m_eigen_geometry.data()[3 * i + 0] + m_dT * m_eigen_velocities.data()[3 * i + 0] - 0.5 * m_eigen_gradient.data()[3 * i + 0] * m_eigen_inv_masses.data()[3 * i + 0] * m_dt2;
        m_eigen_geometry.data()[3 * i + 1] = m_eigen_geometry.data()[3 * i + 1] + m_dT * m_eigen_velocities.data()[3 * i + 1] - 0.5 * m_eigen_gradient.data()[3 * i + 1] * m_eigen_inv_masses.data()[3 * i + 1] * m_dt2;
        m_eigen_geometry.data()[3 * i + 2] = m_eigen_geometry.data()[3 * i + 2] + m_dT * m_eigen_velocities.data()[3 * i + 2] - 0.5 * m_eigen_gradient.data()[3 * i + 2] * m_eigen_inv_masses.data()[3 * i + 2] * m_dt2;

        m_eigen_velocities.data()[3 * i + 0] = m_eigen_velocities.data()[3 * i + 0] - 0.5 * m_dT * m_eigen_gradient.data()[3 * i + 0] * m_eigen_inv_masses.data()[3 * i + 0];
        m_eigen_velocities.data()[3 * i + 1] = m_eigen_velocities.data()[3 * i + 1] - 0.5 * m_dT * m_eigen_gradient.data()[3 * i + 1] * m_eigen_inv_masses.data()[3 * i + 1];
        m_eigen_velocities.data()[3 * i + 2] = m_eigen_velocities.data()[3 * i + 2] - 0.5 * m_dT * m_eigen_gradient.data()[3 * i + 2] * m_eigen_inv_masses.data()[3 * i + 2];
        ekin += m_eigen_masses.data()[3 * i] * (m_eigen_velocities.data()[3 * i] * m_eigen_velocities.data()[3 * i] + m_eigen_velocities.data()[3 * i + 1] * m_eigen_velocities.data()[3 * i + 1] + m_eigen_velocities.data()[3 * i + 2] * m_eigen_velocities.data()[3 * i + 2]);
    }

    // CRITICAL FIX (Feb 2026): Check velocities AFTER integration
    // Catch NaN/Inf that may result from extreme gradient values
    for (int i = 0; i < 3 * m_natoms; ++i) {
        if (!std::isfinite(m_eigen_velocities.data()[i])) {
            CurcumaLogger::error("NaN/Inf velocity at coordinate " + std::to_string(i) +
                               " (atom " + std::to_string(i/3) + ") after integration");
            m_unstable = true;
            return;
        }
    }
    ekin *= 0.5;
    m_T = 2.0 * ekin / (kb_Eh * m_dof);
    m_Ekin = ekin;
    // Thermostat NOT applied here — half-step velocities not representative.
    // Apply thermostat only after full velocity update (second half-step).
    m_Epot = Energy();
    if (m_rmsd_mtd) {
        if (m_step % m_mtd_steps == 0) {
            ApplyRMSDMTD();
        }
    }
#ifdef USE_Plumed
    if (m_mtd) {
        plumed_cmd(m_plumedmain, "setStep", &m_step);

        plumed_cmd(m_plumedmain, "setPositions", &m_eigen_geometry.data()[0]);

        plumed_cmd(m_plumedmain, "setEnergy", &m_Epot);
        plumed_cmd(m_plumedmain, "setForces", &m_eigen_gradient.data()[0]);
        plumed_cmd(m_plumedmain, "setVirial", &m_virial[0]);

        plumed_cmd(m_plumedmain, "setMasses", &m_eigen_masses.data()[0]);
        if (m_eval_mtd) {
            plumed_cmd(m_plumedmain, "prepareCalc", NULL);
            plumed_cmd(m_plumedmain, "performCalc", NULL);
        } else {
            if (std::abs(m_T0 - m_aver_Temp) < m_mtd_dT && m_step > 10) {
                m_eval_mtd = true;
                std::cout << "Starting with MetaDynamics ..." << std::endl;
            }
        }
    }
#endif
    WallPotential();
    ekin = 0.0;

    for (int i = 0; i < m_natoms; ++i) {
        m_eigen_velocities.data()[3 * i + 0] -= 0.5 * m_dT * m_eigen_gradient.data()[3 * i + 0] * m_eigen_inv_masses.data()[3 * i + 0];
        m_eigen_velocities.data()[3 * i + 1] -= 0.5 * m_dT * m_eigen_gradient.data()[3 * i + 1] * m_eigen_inv_masses.data()[3 * i + 1];
        m_eigen_velocities.data()[3 * i + 2] -= 0.5 * m_dT * m_eigen_gradient.data()[3 * i + 2] * m_eigen_inv_masses.data()[3 * i + 2];

        ekin += m_eigen_masses.data()[3 * i] * (m_eigen_velocities.data()[3 * i] * m_eigen_velocities.data()[3 * i] + m_eigen_velocities.data()[3 * i + 1] * m_eigen_velocities.data()[3 * i + 1] + m_eigen_velocities.data()[3 * i + 2] * m_eigen_velocities.data()[3 * i + 2]);
        //m_gradient[3 * i + 0] = m_eigen_gradient.data()[3 * i + 0];
        //m_gradient[3 * i + 1] = m_eigen_gradient.data()[3 * i + 1];
        //m_gradient[3 * i + 2] = m_eigen_gradient.data()[3 * i + 2];
    }
    ekin *= 0.5;
    double T = 2.0 * ekin / (kb_Eh * m_dof);
    m_unstable = T > 10000 * m_T || std::isnan(T);
    m_T = T;
    m_Ekin = ekin;
    ThermostatFunction();
    EKin();

    // Claude Generated (Oct 2025): Apply PBC wrapping after integration step
    applyPeriodicBoundaryConditions();
}

void SimpleMD::Rattle()
{
    /* this part was adopted from
     * Numerische Simulation in der Moleküldynamik
     * by
     * Griebel, Knapek, Zumbusch, Caglar
     * 2003, Springer-Verlag
     * and from
     * Molecular Simulation of Fluids
     * by Richard J. Sadus
     * some suff was just ignored or corrected
     * like dT^3 -> dT^2 and
     * updated velocities of the second atom (minus instead of plus)
     * and adjusted to some needs
     */
    TriggerWriteRestart();

    auto* coord = new double[3 * m_natoms];
    double m_dT_inverse = 1 / m_dT;
    std::set<int> constrained_atoms;
    bool move = false;
    double max_mu = 10;
    double max_err_12 = 0, max_err_13 = 0;
    for (int i = 0; i < m_natoms; ++i) {
        coord[3 * i + 0] = m_eigen_geometry.data()[3 * i + 0] + m_dT * m_eigen_velocities.data()[3 * i + 0] - 0.5 * m_eigen_gradient.data()[3 * i + 0] * m_eigen_inv_masses.data()[3 * i + 0] * m_dt2;
        coord[3 * i + 1] = m_eigen_geometry.data()[3 * i + 1] + m_dT * m_eigen_velocities.data()[3 * i + 1] - 0.5 * m_eigen_gradient.data()[3 * i + 1] * m_eigen_inv_masses.data()[3 * i + 1] * m_dt2;
        coord[3 * i + 2] = m_eigen_geometry.data()[3 * i + 2] + m_dT * m_eigen_velocities.data()[3 * i + 2] - 0.5 * m_eigen_gradient.data()[3 * i + 2] * m_eigen_inv_masses.data()[3 * i + 2] * m_dt2;

        m_eigen_velocities.data()[3 * i + 0] -= 0.5 * m_dT * m_eigen_gradient.data()[3 * i + 0] * m_eigen_inv_masses.data()[3 * i + 0];
        m_eigen_velocities.data()[3 * i + 1] -= 0.5 * m_dT * m_eigen_gradient.data()[3 * i + 1] * m_eigen_inv_masses.data()[3 * i + 1];
        m_eigen_velocities.data()[3 * i + 2] -= 0.5 * m_dT * m_eigen_gradient.data()[3 * i + 2] * m_eigen_inv_masses.data()[3 * i + 2];
    }

    double iter = 0;
    double max_violation_12 = 0, max_violation_13 = 0;
    while (iter < m_rattle_maxiter) {
        iter++;
        int active = 0;
        max_violation_12 = 0;
        max_violation_13 = 0;

        for (auto bond : m_bond_constrained) {
            int i = bond.first.first, j = bond.first.second;
            double distance = bond.second;
            double distance_current = ((coord[3 * i + 0] - coord[3 * j + 0]) * (coord[3 * i + 0] - coord[3 * j + 0])
                + (coord[3 * i + 1] - coord[3 * j + 1]) * (coord[3 * i + 1] - coord[3 * j + 1])
                + (coord[3 * i + 2] - coord[3 * j + 2]) * (coord[3 * i + 2] - coord[3 * j + 2]));
            double violation_12 = std::abs(distance - distance_current);
            if (violation_12 > max_violation_12) max_violation_12 = violation_12;
            if (violation_12 > m_rattle_tol_12) {
                if (violation_12 > max_err_12) max_err_12 = violation_12;
                move = true;
                double r = distance - distance_current;
                double dx = m_eigen_geometry.data()[3 * i + 0] - m_eigen_geometry.data()[3 * j + 0];
                double dy = m_eigen_geometry.data()[3 * i + 1] - m_eigen_geometry.data()[3 * j + 1];
                double dz = m_eigen_geometry.data()[3 * i + 2] - m_eigen_geometry.data()[3 * j + 2];

                double scalarproduct = (dx) * (coord[3 * i + 0] - coord[3 * j + 0])
                    + (dy) * (coord[3 * i + 1] - coord[3 * j + 1])
                    + (dz) * (coord[3 * i + 2] - coord[3 * j + 2]);

                constrained_atoms.insert(i);
                constrained_atoms.insert(j);
                active++;

                if (std::abs(scalarproduct) < m_rattle_min) {
                    if (scalarproduct < 0)
                        scalarproduct = -1 * m_rattle_min;
                    else
                        scalarproduct = m_rattle_min;
                }

                double lambda = r / ((m_eigen_inv_masses.data()[3 * i] + m_eigen_inv_masses.data()[3 * j]) * scalarproduct);
                if (std::isinf(lambda) || std::isnan(lambda)) {
                    std::cout << "RATTLE 1-2: " << i << " " << j << " lambda=" << lambda
                              << " r=" << r << " sp=" << scalarproduct << " dc=" << distance_current << std::endl;
                }

                while (std::abs(lambda) > max_mu)
                    lambda /= 2;

                coord[3 * i + 0] += dx * lambda * 0.5 * m_eigen_inv_masses.data()[3 * i];
                coord[3 * i + 1] += dy * lambda * 0.5 * m_eigen_inv_masses.data()[3 * i];
                coord[3 * i + 2] += dz * lambda * 0.5 * m_eigen_inv_masses.data()[3 * i];

                coord[3 * j + 0] -= dx * lambda * 0.5 * m_eigen_inv_masses.data()[3 * j];
                coord[3 * j + 1] -= dy * lambda * 0.5 * m_eigen_inv_masses.data()[3 * j];
                coord[3 * j + 2] -= dz * lambda * 0.5 * m_eigen_inv_masses.data()[3 * j];

                m_eigen_velocities.data()[3 * i + 0] += dx * lambda * 0.5 * m_eigen_inv_masses.data()[3 * i] * m_dT_inverse;
                m_eigen_velocities.data()[3 * i + 1] += dy * lambda * 0.5 * m_eigen_inv_masses.data()[3 * i] * m_dT_inverse;
                m_eigen_velocities.data()[3 * i + 2] += dz * lambda * 0.5 * m_eigen_inv_masses.data()[3 * i] * m_dT_inverse;

                m_eigen_velocities.data()[3 * j + 0] -= dx * lambda * 0.5 * m_eigen_inv_masses.data()[3 * j] * m_dT_inverse;
                m_eigen_velocities.data()[3 * j + 1] -= dy * lambda * 0.5 * m_eigen_inv_masses.data()[3 * j] * m_dT_inverse;
                m_eigen_velocities.data()[3 * j + 2] -= dz * lambda * 0.5 * m_eigen_inv_masses.data()[3 * j] * m_dT_inverse;
            }
        }

        for (auto bond : m_bond_13_constrained) {
            int i = bond.first.first, j = bond.first.second;
            double distance = bond.second;
            double distance_current = ((coord[3 * i + 0] - coord[3 * j + 0]) * (coord[3 * i + 0] - coord[3 * j + 0])
                + (coord[3 * i + 1] - coord[3 * j + 1]) * (coord[3 * i + 1] - coord[3 * j + 1])
                + (coord[3 * i + 2] - coord[3 * j + 2]) * (coord[3 * i + 2] - coord[3 * j + 2]));
            double violation_13 = std::abs(distance - distance_current);
            if (violation_13 > max_violation_13) max_violation_13 = violation_13;
            if (violation_13 > m_rattle_tol_13) {
                if (violation_13 > max_err_13) max_err_13 = violation_13;
                move = true;
                double r = distance - distance_current;
                double dx = m_eigen_geometry.data()[3 * i + 0] - m_eigen_geometry.data()[3 * j + 0];
                double dy = m_eigen_geometry.data()[3 * i + 1] - m_eigen_geometry.data()[3 * j + 1];
                double dz = m_eigen_geometry.data()[3 * i + 2] - m_eigen_geometry.data()[3 * j + 2];

                double scalarproduct = (dx) * (coord[3 * i + 0] - coord[3 * j + 0])
                    + (dy) * (coord[3 * i + 1] - coord[3 * j + 1])
                    + (dz) * (coord[3 * i + 2] - coord[3 * j + 2]);

                constrained_atoms.insert(i);
                constrained_atoms.insert(j);
                active++;

                if (std::abs(scalarproduct) < m_rattle_min) {
                    if (scalarproduct < 0)
                        scalarproduct = -1 * m_rattle_min;
                    else
                        scalarproduct = m_rattle_min;
                }

                double lambda = r / ((m_eigen_inv_masses.data()[3 * i] + m_eigen_inv_masses.data()[3 * j]) * scalarproduct);
                if (std::isinf(lambda) || std::isnan(lambda)) {
                    std::cout << "RATTLE 1-3: " << i << " " << j << " lambda=" << lambda
                              << " r=" << r << " sp=" << scalarproduct << " dc=" << distance_current << std::endl;
                }

                while (std::abs(lambda) > max_mu)
                    lambda /= 2;

                coord[3 * i + 0] += dx * lambda * 0.5 * m_eigen_inv_masses.data()[3 * i];
                coord[3 * i + 1] += dy * lambda * 0.5 * m_eigen_inv_masses.data()[3 * i];
                coord[3 * i + 2] += dz * lambda * 0.5 * m_eigen_inv_masses.data()[3 * i];

                coord[3 * j + 0] -= dx * lambda * 0.5 * m_eigen_inv_masses.data()[3 * j];
                coord[3 * j + 1] -= dy * lambda * 0.5 * m_eigen_inv_masses.data()[3 * j];
                coord[3 * j + 2] -= dz * lambda * 0.5 * m_eigen_inv_masses.data()[3 * j];

                m_eigen_velocities.data()[3 * i + 0] += dx * lambda * 0.5 * m_eigen_inv_masses.data()[3 * i] * m_dT_inverse;
                m_eigen_velocities.data()[3 * i + 1] += dy * lambda * 0.5 * m_eigen_inv_masses.data()[3 * i] * m_dT_inverse;
                m_eigen_velocities.data()[3 * i + 2] += dz * lambda * 0.5 * m_eigen_inv_masses.data()[3 * i] * m_dT_inverse;

                m_eigen_velocities.data()[3 * j + 0] -= dx * lambda * 0.5 * m_eigen_inv_masses.data()[3 * j] * m_dT_inverse;
                m_eigen_velocities.data()[3 * j + 1] -= dy * lambda * 0.5 * m_eigen_inv_masses.data()[3 * j] * m_dT_inverse;
                m_eigen_velocities.data()[3 * j + 2] -= dz * lambda * 0.5 * m_eigen_inv_masses.data()[3 * j] * m_dT_inverse;
            }
        }
        /*
                Geometry geometry = m_molecule.getGeometry();
                for (int i = 0; i < m_natoms; ++i) {
                    geometry(i, 0) = m_rt_geom_1[3 * i + 0];
                    geometry(i, 1) = m_rt_geom_1[3 * i + 1];
                    geometry(i, 2) = m_rt_geom_1[3 * i + 2];
                }
                m_molecule.setGeometry(geometry);

                m_molecule.appendXYZFile(outputPath(Basename() + ".rattle.trj.xyz"));
        */
        if (active == 0)
            break;
    }

    if (iter >= m_rattle_maxiter) {
        CurcumaLogger::info("RATTLE 1st step: max iterations reached");
        m_rattle_max_err_count++;
    }
    m_rattle_iters_step1 = static_cast<int>(iter);
    m_rattle_constrained_atoms = static_cast<int>(constrained_atoms.size());

    double ekin = 0;

    for (int i = 0; i < m_natoms; ++i) {
        m_eigen_geometry.data()[3 * i + 0] = coord[3 * i + 0];
        m_eigen_geometry.data()[3 * i + 1] = coord[3 * i + 1];
        m_eigen_geometry.data()[3 * i + 2] = coord[3 * i + 2];
        ekin += m_eigen_masses.data()[3 * i] * (m_eigen_velocities.data()[3 * i] * m_eigen_velocities.data()[3 * i] + m_eigen_velocities.data()[3 * i + 1] * m_eigen_velocities.data()[3 * i + 1] + m_eigen_velocities.data()[3 * i + 2] * m_eigen_velocities.data()[3 * i + 2]);
    }
    ekin *= 0.5;
    m_T = 2.0 * ekin / (kb_Eh * m_dof);
    m_Ekin = ekin;
    // Thermostat NOT applied here — half-step velocities include SHAKE corrections.
    // Apply thermostat only after full RATTLE step (velocity constraints + 2nd half-step).
    m_Epot = Energy();

    if (m_rmsd_mtd) {
        if (m_step % m_mtd_steps == 0) {
            ApplyRMSDMTD();
        }
    }
#ifdef USE_Plumed
    if (m_mtd) {
        plumed_cmd(m_plumedmain, "setStep", &m_step);

        plumed_cmd(m_plumedmain, "setPositions", &m_eigen_geometry.data()[0]);

        plumed_cmd(m_plumedmain, "setEnergy", &m_Epot);
        plumed_cmd(m_plumedmain, "setForces", &m_eigen_gradient.data()[0]);
        plumed_cmd(m_plumedmain, "setVirial", &m_virial[0]);

        plumed_cmd(m_plumedmain, "setMasses", &m_eigen_masses.data()[0]);
        if (m_eval_mtd) {
            plumed_cmd(m_plumedmain, "prepareCalc", NULL);
            plumed_cmd(m_plumedmain, "performCalc", NULL);
        } else {
            if (std::abs(m_T0 - m_aver_Temp) < m_mtd_dT && m_step > 10) {
                m_eval_mtd = true;
                std::cout << "Starting with MetaDynamics ..." << std::endl;
            }
        }
    }
#endif
    WallPotential();

    for (int i = 0; i < m_natoms; ++i) {
        m_eigen_velocities.data()[3 * i + 0] -= 0.5 * m_dT * m_eigen_gradient.data()[3 * i + 0] * m_eigen_inv_masses.data()[3 * i + 0];
        m_eigen_velocities.data()[3 * i + 1] -= 0.5 * m_dT * m_eigen_gradient.data()[3 * i + 1] * m_eigen_inv_masses.data()[3 * i + 1];
        m_eigen_velocities.data()[3 * i + 2] -= 0.5 * m_dT * m_eigen_gradient.data()[3 * i + 2] * m_eigen_inv_masses.data()[3 * i + 2];

        //m_gradient[3 * i + 0] = m_eigen_gradient.data()[3 * i + 0];
        //m_gradient[3 * i + 1] = m_eigen_gradient.data()[3 * i + 1];
        //m_gradient[3 * i + 2] = m_eigen_gradient.data()[3 * i + 2];
    }
    m_virial_correction = 0;
    iter = 0;
    ekin = 0.0;

    while (iter < m_rattle_maxiter) {
        iter++;
        int active = 0;
        for (auto bond : m_bond_constrained) {
            int i = bond.first.first, j = bond.first.second;
            if (constrained_atoms.count(i) && constrained_atoms.count(j)) {
                double distance_current = ((coord[3 * i + 0] - coord[3 * j + 0]) * (coord[3 * i + 0] - coord[3 * j + 0])
                    + (coord[3 * i + 1] - coord[3 * j + 1]) * (coord[3 * i + 1] - coord[3 * j + 1])
                    + (coord[3 * i + 2] - coord[3 * j + 2]) * (coord[3 * i + 2] - coord[3 * j + 2]));

                double distance = bond.second;

                double dx = coord[3 * i + 0] - coord[3 * j + 0];
                double dy = coord[3 * i + 1] - coord[3 * j + 1];
                double dz = coord[3 * i + 2] - coord[3 * j + 2];
                double dvx = m_eigen_velocities.data()[3 * i + 0] - m_eigen_velocities.data()[3 * j + 0];
                double dvy = m_eigen_velocities.data()[3 * i + 1] - m_eigen_velocities.data()[3 * j + 1];
                double dvz = m_eigen_velocities.data()[3 * i + 2] - m_eigen_velocities.data()[3 * j + 2];

                double r = (dx) * (dvx) + (dy) * (dvy) + (dz) * (dvz);

                double mu = -1 * r / ((m_eigen_inv_masses.data()[3 * i] + m_eigen_inv_masses.data()[3 * j]) * distance_current);
                while (std::abs(mu) > max_mu)
                    mu /= 2;
                if (std::abs(mu) > m_rattle_tol_12) {
                    active = 1;
                    m_virial_correction += mu * distance_current;
                    m_eigen_velocities.data()[3 * i + 0] += dx * mu * m_eigen_inv_masses.data()[3 * i];
                    m_eigen_velocities.data()[3 * i + 1] += dy * mu * m_eigen_inv_masses.data()[3 * i];
                    m_eigen_velocities.data()[3 * i + 2] += dz * mu * m_eigen_inv_masses.data()[3 * i];

                    m_eigen_velocities.data()[3 * j + 0] -= dx * mu * m_eigen_inv_masses.data()[3 * j];
                    m_eigen_velocities.data()[3 * j + 1] -= dy * mu * m_eigen_inv_masses.data()[3 * j];
                    m_eigen_velocities.data()[3 * j + 2] -= dz * mu * m_eigen_inv_masses.data()[3 * j];
                }
            }
        }

        for (auto bond : m_bond_13_constrained) {
            int i = bond.first.first, j = bond.first.second;
            if (constrained_atoms.count(i) && constrained_atoms.count(j)) {
                double distance_current = ((coord[3 * i + 0] - coord[3 * j + 0]) * (coord[3 * i + 0] - coord[3 * j + 0])
                    + (coord[3 * i + 1] - coord[3 * j + 1]) * (coord[3 * i + 1] - coord[3 * j + 1])
                    + (coord[3 * i + 2] - coord[3 * j + 2]) * (coord[3 * i + 2] - coord[3 * j + 2]));

                double distance = bond.second;

                double dx = coord[3 * i + 0] - coord[3 * j + 0];
                double dy = coord[3 * i + 1] - coord[3 * j + 1];
                double dz = coord[3 * i + 2] - coord[3 * j + 2];
                double dvx = m_eigen_velocities.data()[3 * i + 0] - m_eigen_velocities.data()[3 * j + 0];
                double dvy = m_eigen_velocities.data()[3 * i + 1] - m_eigen_velocities.data()[3 * j + 1];
                double dvz = m_eigen_velocities.data()[3 * i + 2] - m_eigen_velocities.data()[3 * j + 2];

                double r = (dx) * (dvx) + (dy) * (dvy) + (dz) * (dvz);

                double mu = -1 * r / ((m_eigen_inv_masses.data()[3 * i] + m_eigen_inv_masses.data()[3 * j]) * distance_current);
                while (std::abs(mu) > max_mu)
                    mu /= 2;
                if (std::abs(mu) > m_rattle_tol_13) {
                    active = 1;
                    m_virial_correction += mu * distance_current;
                    m_eigen_velocities.data()[3 * i + 0] += dx * mu * m_eigen_inv_masses.data()[3 * i];
                    m_eigen_velocities.data()[3 * i + 1] += dy * mu * m_eigen_inv_masses.data()[3 * i];
                    m_eigen_velocities.data()[3 * i + 2] += dz * mu * m_eigen_inv_masses.data()[3 * i];

                    m_eigen_velocities.data()[3 * j + 0] -= dx * mu * m_eigen_inv_masses.data()[3 * j];
                    m_eigen_velocities.data()[3 * j + 1] -= dy * mu * m_eigen_inv_masses.data()[3 * j];
                    m_eigen_velocities.data()[3 * j + 2] -= dz * mu * m_eigen_inv_masses.data()[3 * j];
                }
            }
        }
        if (active == 0)
            break;
    }

    if (iter >= m_rattle_maxiter) {
        CurcumaLogger::info("RATTLE 2nd step: max iterations reached (" + std::to_string(static_cast<int>(iter)) + ")");
        m_rattle_max_err_count++;
    }
    m_rattle_iters_step2 = static_cast<int>(iter);

    // Accumulate worst errors for periodic summary
    if (max_err_12 > m_rattle_max_err_12) m_rattle_max_err_12 = max_err_12;
    if (max_err_13 > m_rattle_max_err_13) m_rattle_max_err_13 = max_err_13;

    if (move)
        RemoveRotations();

    delete[] coord;
    for (int i = 0; i < m_natoms; ++i) {
        ekin += m_eigen_masses.data()[3 * i] * (m_eigen_velocities.data()[3 * i] * m_eigen_velocities.data()[3 * i] + m_eigen_velocities.data()[3 * i + 1] * m_eigen_velocities.data()[3 * i + 1] + m_eigen_velocities.data()[3 * i + 2] * m_eigen_velocities.data()[3 * i + 2]);
    }
    ekin *= 0.5;
    double T = 2.0 * ekin / (kb_Eh * m_dof);
    m_unstable = T > 10000 * m_T || std::isnan(T);
    m_T = T;
    m_Ekin = ekin;
    ThermostatFunction();
    EKin();

    // Claude Generated (Oct 2025): Apply PBC wrapping after integration step
    applyPeriodicBoundaryConditions();
}

void SimpleMD::ApplyRMSDMTD()
{
    std::chrono::time_point<std::chrono::system_clock> m_start, m_end;
    m_start = std::chrono::system_clock::now();
    m_colvar_incr = 0;

    // RMSD-subset geometry for bias evaluation
    Geometry current_geometry = m_rmsd_mtd_molecule.getGeometry();
    for (int i = 0; i < m_rmsd_indicies.size(); ++i) {
        current_geometry(i, 0) = m_eigen_geometry.data()[3 * m_rmsd_indicies[i] + 0];
        current_geometry(i, 1) = m_eigen_geometry.data()[3 * m_rmsd_indicies[i] + 1];
        current_geometry(i, 2) = m_eigen_geometry.data()[3 * m_rmsd_indicies[i] + 2];
    }

    // Full molecule geometry for storage and XYZ output (all atoms)
    Geometry full_geometry(m_natoms, 3);
    for (int i = 0; i < m_natoms; ++i) {
        full_geometry(i, 0) = m_eigen_geometry.data()[3 * i + 0];
        full_geometry(i, 1) = m_eigen_geometry.data()[3 * i + 1];
        full_geometry(i, 2) = m_eigen_geometry.data()[3 * i + 2];
    }

    double current_bias = 0;    // exploration bias: drives force + deposition
    double current_bias_wt = 0; // optional well-tempered energy (opt-in, COLVAR output only)
    double rmsd_reference = 0;

    // Claude Generated (Apr 2026): Shared bias pool path for parallel ConfSearch
    // When a shared pool is set, read bias structures from the pool and evaluate locally.
    // Deposit new structures back to the shared pool when the deposition criterion is met.
    if (m_shared_pool) {
        int global_count = m_shared_pool->biasStructureCount();

        if (global_count == 0) {
            // First structure: deposit initial reference with full geometry
            BiasStructure initial;
            initial.geometry = full_geometry;
            initial.time = m_currentStep;
            initial.rmsd_reference = 0;
            initial.counter = 1;
            initial.index = 0;
            initial.temperature = m_T0;
            int deposited = m_shared_pool->depositBiasStructure(initial);
            m_bias_structure_count++;
            CurcumaLogger::result_fmt("RMSD-MTD: Initial bias structure {} deposited (pool total: {})",
                deposited, m_shared_pool->biasStructureCount());
            // Write full molecule to per-thread .mtd.xyz for reference
            Molecule out_mol(m_molecule);
            out_mol.setGeometry(full_geometry);
            out_mol.setName(std::to_string(m_currentStep));
            out_mol.writeXYZFile(Basename() + ".mtd.xyz");
            if (m_nocolvarfile == false) {
                std::ofstream colvarfile;
                colvarfile.open(outputPath("COLVAR"));
                colvarfile.close();
            }
            m_end = std::chrono::system_clock::now();
            m_mtd_time += std::chrono::duration_cast<std::chrono::milliseconds>(m_end - m_start).count();
            return;
        }

        // Snapshot all current bias structures from the shared pool
        auto bias_snapshot = m_shared_pool->snapshot();

        // Reference for the RMSD driver must be the CURRENT walker geometry, set once
        // per step. Without this the driver keeps the initial geometry from Initialise()
        // (setReference copies), so RMSD would be measured against the start structure
        // (structure 0 -> rmsd~0 -> current_bias huge -> deposition never triggers and
        // the bias force is bogus). This mirrors the local BiasThread path (setReference
        // of the moving geometry each step). Claude Generated (Jun 2026).
        m_rmsd_mtd_molecule.setGeometry(current_geometry);
        m_shared_pool_driver.setReference(m_rmsd_mtd_molecule);

        // Visited references to bump after the loop (the WT weight needs the full V).
        std::vector<int> visited;

        // Symmetry/atom-permutation set discovered by ConfScan (full-atom reorder rules), set
        // on the shared pool between cycles. Empty -> identity only -> bit-identical to before.
        // We SUM a Gaussian over every symmetry image (NOT a hard min over images): the bias and
        // its force stay C-infinity smooth as the walker crosses a symmetry seam (a hard min would
        // make the force discontinuous there and break the MD integration). For well-separated
        // images only the nearest contributes appreciably, so it behaves like the min but smoothly.
        // Only applied when the RMSD subset is the full molecule (the reorder rules are full-atom).
        // Claude Generated (Jun 2026).
        std::vector<std::vector<int>> perms;
        if (static_cast<int>(m_rmsd_indicies.size()) == m_natoms)
            perms = m_shared_pool->permutations();

        // Claude Generated (Jun 2026): flexibility/RMSF weights (Phase C "weighted"). Empty ->
        // uniform -> standard best-fit RMSD (bit-identical). Only when the RMSD subset is the full
        // molecule (the weights are full-atom). Set once here; persists for every image below.
        if (static_cast<int>(m_rmsd_indicies.size()) == m_natoms) {
            std::vector<double> w = m_shared_pool->weights();
            if (static_cast<int>(w.size()) == m_natoms)
                m_shared_pool_driver.setRMSDWeights(w);
            else
                m_shared_pool_driver.clearRMSDWeights();
        }

        // Evaluate the bias from the snapshot — hill height W_i = k * counter_i,
        // V(x) = Sum_i Sum_p W_i * exp(-alpha*RMSD_{i,p}^2) (p over identity + symmetry images),
        // force = exact negative gradient. Well-tempered (opt-in) only feeds current_bias_wt.
        for (const auto& bs : bias_snapshot) {
            // Effective hill counter: frozen if inherited at run start (only this run's deposits
            // grow), then capped by rmsd_mtd_max_height. Both default-off -> eff = bs.counter
            // (legacy W_i = k * counter_i). Claude Generated (Jun 2026).
            int eff_counter = bs.counter;
            if (m_freeze_inherited) {
                auto it = m_frozen_height.find(bs.index);
                if (it != m_frozen_height.end())
                    eff_counter = it->second;
            }
            if (m_rmsd_mtd_max_height > 0)
                eff_counter = std::min(eff_counter, m_rmsd_mtd_max_height);
            const double height = m_k_rmsd * eff_counter; // W_i = k * counter_i
            double expr_sum = 0.0;      // Sum over images: drives deposition/visited bookkeeping
            double rmsd_identity = 0.0; // identity-image RMSD for COLVAR / rmsd_reference

            // Evaluate one image (a reordered copy of the bias structure subset) and accumulate
            // its Gaussian into the bias + its analytic force into the walker gradient.
            auto eval_image = [&](const Geometry& subset, bool is_identity) {
                m_shared_pool_target.setGeometry(subset);
                m_shared_pool_driver.setTarget(m_shared_pool_target);
                double rmsd = m_shared_pool_driver.BestFitRMSD();
                double expr = exp(-rmsd * rmsd * m_alpha_rmsd);
                current_bias += height * expr;
                if (m_wtmtd)
                    current_bias_wt += m_k_rmsd * bs.factor * expr;
                double dEdR = -2.0 * m_alpha_rmsd * height * rmsd * expr;
                Geometry grad = m_shared_pool_driver.Gradient();
                for (int j = 0; j < m_rmsd_indicies.size(); ++j) {
                    m_eigen_gradient.data()[3 * m_rmsd_indicies[j] + 0] += dEdR * grad(j, 0);
                    m_eigen_gradient.data()[3 * m_rmsd_indicies[j] + 1] += dEdR * grad(j, 1);
                    m_eigen_gradient.data()[3 * m_rmsd_indicies[j] + 2] += dEdR * grad(j, 2);
                }
                expr_sum += expr;
                if (is_identity)
                    rmsd_identity = rmsd;
            };

            // Identity image (bs.geometry is full-atom; project onto the RMSD subset).
            Geometry bs_rmsd_subset = m_shared_pool_target.getGeometry();
            for (int i = 0; i < m_rmsd_indicies.size(); ++i) {
                bs_rmsd_subset(i, 0) = bs.geometry(m_rmsd_indicies[i], 0);
                bs_rmsd_subset(i, 1) = bs.geometry(m_rmsd_indicies[i], 1);
                bs_rmsd_subset(i, 2) = bs.geometry(m_rmsd_indicies[i], 2);
            }
            eval_image(bs_rmsd_subset, true);

            // Symmetry images: position j holds atom rule[j] (same convention as ConfScan's
            // Rules2RMSD). Reference (walker) stays in canonical order, so this measures the
            // RMSD to the relabelled bias structure.
            for (const auto& rule : perms) {
                if (static_cast<int>(rule.size()) != m_natoms)
                    continue;
                Geometry bs_perm(m_natoms, 3);
                bool ok = true;
                for (int j = 0; j < m_natoms; ++j) {
                    int src = rule[j];
                    if (src < 0 || src >= m_natoms) { ok = false; break; }
                    bs_perm(j, 0) = bs.geometry(src, 0);
                    bs_perm(j, 1) = bs.geometry(src, 1);
                    bs_perm(j, 2) = bs.geometry(src, 2);
                }
                if (ok)
                    eval_image(bs_perm, false);
            }

            if (expr_sum * m_rmsd_econv > static_cast<double>(global_count))
                visited.push_back(bs.index);

            if (bs.index == 0)
                rmsd_reference = rmsd_identity;
        }

        // Phase 2: bump the visited references in the shared pool. counter++ always (drives
        // exploration); the well-tempered weight 'factor' grows only when opt-in (output-only).
        {
            std::vector<std::pair<int, double>> visit_updates;
            visit_updates.reserve(visited.size());
            double wt_inc = m_wtmtd ? exp(-current_bias / (kb_Eh * m_rmsd_DT)) : 0.0;
            for (int idx : visited) {
                // When freezing inherited heights, do not bump structures this run inherited:
                // their height stays fixed so the run's cumulative bias does not escalate.
                if (m_freeze_inherited && m_frozen_height.count(idx))
                    continue;
                visit_updates.emplace_back(idx, wt_inc);
            }
            m_shared_pool->registerVisits(visit_updates);
        }

        m_rmsd_mtd_molecule.setGeometry(current_geometry);

        // COLVAR output
        if (m_nocolvarfile == false) {
            std::ofstream colvarfile;
            colvarfile.open(outputPath("COLVAR"), std::iostream::app);
            colvarfile << m_currentStep << " ";
            if (m_rmsd_fragment_count < 2)
                colvarfile << rmsd_reference << " ";
            for (int i = 0; i < m_rmsd_fragment_count; ++i)
                for (int j = 0; j < i; ++j) {
                    colvarfile << (m_rmsd_mtd_molecule.Centroid(true, i) - m_rmsd_mtd_molecule.Centroid(true, j)).norm() << " ";
                }
            colvarfile << (m_wtmtd ? current_bias_wt : current_bias) << " " << std::endl;
            colvarfile.close();
        }
        m_bias_energy += current_bias;

        // Deposition criterion: same as local path — only deposit when in a genuinely new region.
        // current_bias * m_rmsd_econv < pool_count means: the region is under-biased relative
        // to the number of existing reference structures → it is a new conformation.
        // First structure (pool empty) is always accepted.
        int pool_count = m_shared_pool->biasStructureCount();
        bool deposit = !m_rmsd_fix_structure
            && (pool_count == 0 || current_bias * m_rmsd_econv < static_cast<double>(pool_count));
        // Claude Generated (Jun 2026): never deposit a fragmented structure into the shared
        // pool (only relevant when topo_check is on; the run aborts shortly after anyway).
        if (deposit && m_topo_check) {
            m_molecule.setGeometry(full_geometry);
            if (static_cast<int>(m_molecule.GetFragments().size()) > m_start_fragment_count)
                deposit = false;
        }
        if (deposit) {
            BiasStructure new_bs;
            new_bs.geometry = full_geometry;
            new_bs.rmsd_reference = rmsd_reference;
            new_bs.time = m_currentStep;
            new_bs.counter = 1;
            new_bs.temperature = m_T0;
            int new_count = m_shared_pool->depositBiasStructure(new_bs);
            m_bias_structure_count++;
            // Write full molecule to per-thread .mtd.xyz
            Molecule out_mol(m_molecule);
            out_mol.setGeometry(full_geometry);
            out_mol.setName(std::to_string(m_currentStep));
            out_mol.appendXYZFile(Basename() + ".mtd.xyz");
            if (CurcumaLogger::get_verbosity() >= 2)
                CurcumaLogger::result_fmt("RMSD-MTD: Deposited bias structure {} (pool total: {})",
                    new_count, m_shared_pool->biasStructureCount());
        }

        m_end = std::chrono::system_clock::now();
        m_mtd_time += std::chrono::duration_cast<std::chrono::milliseconds>(m_end - m_start).count();
        return;
    }

    // Original local-only bias path (unchanged)
    if (m_bias_structure_count == 0) {
        m_bias_threads[0]->addGeometry(current_geometry, 0, m_currentStep, 0);
        m_bias_structure_count++;
        m_rmsd_mtd_molecule.writeXYZFile(outputPath(Basename() + ".mtd.xyz"));
        if (m_nocolvarfile == false) {
            std::ofstream colvarfile;
            colvarfile.open(outputPath("COLVAR"));
            colvarfile.close();
        }
    }
    if (m_threads == 1 || m_bias_structure_count == 1) {
        for (auto & m_bias_thread : m_bias_threads) {
            m_bias_thread->setCurrentGeometry(current_geometry, m_currentStep);
            m_bias_thread->start();
            current_bias += m_bias_thread->BiasEnergy();
            current_bias_wt += m_bias_thread->BiasEnergyWT();
            for (int j = 0; j < m_rmsd_indicies.size(); ++j) {
                m_eigen_gradient.data()[3 * m_rmsd_indicies[j] + 0] += m_bias_thread->Gradient()(j, 0);
                m_eigen_gradient.data()[3 * m_rmsd_indicies[j] + 1] += m_bias_thread->Gradient()(j, 1);
                m_eigen_gradient.data()[3 * m_rmsd_indicies[j] + 2] += m_bias_thread->Gradient()(j, 2);
            }
            m_colvar_incr += m_bias_thread->Counter();
            m_loop_time += m_bias_thread->getExecutionTime();
        }
    } else {
        if (m_bias_structure_count < m_threads) {
            for (int i = 0; i < m_bias_structure_count; ++i) {
                m_bias_threads[i]->setCurrentGeometry(current_geometry, m_currentStep);
            }
        } else {
            for (auto & m_bias_thread : m_bias_threads) {
                m_bias_thread->setCurrentGeometry(current_geometry, m_currentStep);
            }
        }

        m_bias_pool->setActiveThreadCount(m_threads);
        m_bias_pool->StaticPool();
        m_bias_pool->StartAndWait();
        // m_bias_pool->setWakeUp(m_bias_pool->WakeUp() / 2);

        for (auto & m_bias_thread : m_bias_threads) {
            if (m_bias_thread->getReturnValue() == 1) {

                current_bias += m_bias_thread->BiasEnergy();
                current_bias_wt += m_bias_thread->BiasEnergyWT();
                for (int j = 0; j < m_rmsd_indicies.size(); ++j) {
                    m_eigen_gradient.data()[3 * m_rmsd_indicies[j] + 0] += m_bias_thread->Gradient()(j, 0);
                    m_eigen_gradient.data()[3 * m_rmsd_indicies[j] + 1] += m_bias_thread->Gradient()(j, 1);
                    m_eigen_gradient.data()[3 * m_rmsd_indicies[j] + 2] += m_bias_thread->Gradient()(j, 2);
                }
                m_colvar_incr += m_bias_thread->Counter();
            }
            m_loop_time += m_bias_thread->getExecutionTime();
        }
        m_bias_pool->Reset();
    }
    rmsd_reference = m_bias_threads[0]->RMSDReference();
    m_rmsd_mtd_molecule.setGeometry(current_geometry);

    if (m_nocolvarfile == false) {
        std::ofstream colvarfile;
        colvarfile.open(outputPath("COLVAR"), std::iostream::app);
        colvarfile << m_currentStep << " ";
        if (m_rmsd_fragment_count < 2)
            colvarfile << rmsd_reference << " ";

        for (int i = 0; i < m_rmsd_fragment_count; ++i)
            for (int j = 0; j < i; ++j) {
                colvarfile << (m_rmsd_mtd_molecule.Centroid(true, i) - m_rmsd_mtd_molecule.Centroid(true, j)).norm() << " ";
            }
        // Report the well-tempered energy when opt-in; the exploration bias otherwise.
        colvarfile << (m_wtmtd ? current_bias_wt : current_bias) << " " << std::endl;
        colvarfile.close();
    }
    m_bias_energy += current_bias;

    // Deposition uses the exploration bias only (well-tempering never gates the search).
    if (current_bias * m_rmsd_econv < m_bias_structure_count && m_rmsd_fix_structure == false) {
        int thread_index = m_bias_structure_count % m_bias_threads.size();
        m_bias_threads[thread_index]->addGeometry(current_geometry, rmsd_reference, m_currentStep, m_bias_structure_count);
        m_bias_structure_count++;
        m_rmsd_mtd_molecule.appendXYZFile(outputPath(Basename() + ".mtd.xyz"));
        std::cout << m_bias_structure_count << " stored structures currently" << std::endl;
    }
    m_end = std::chrono::system_clock::now();
    int m_time = std::chrono::duration_cast<std::chrono::milliseconds>(m_end - m_start).count();
    m_mtd_time += m_time;
}

void SimpleMD::Rattle_Verlet_First(double* coord, double* grad)
{
}

void SimpleMD::Rattle_Constrain_First(double* coord, double* grad)
{
}

void SimpleMD::Rattle_Verlet_Second(double* coord, double* grad)
{
}

double SimpleMD::ApplySphericLogFermiWalls()
{
    double potential = 0;
    double kbT = m_wall_temp * kb_Eh;
    int counter = 0;
    double sum_grad = 0; // Claude Generated: Track total wall force
    for (int i = 0; i < m_natoms; ++i) {
        double distance = sqrt(m_eigen_geometry.data()[3 * i + 0] * m_eigen_geometry.data()[3 * i + 0] + m_eigen_geometry.data()[3 * i + 1] * m_eigen_geometry.data()[3 * i + 1] + m_eigen_geometry.data()[3 * i + 2] * m_eigen_geometry.data()[3 * i + 2]);

        // Claude Generated: Add numerical stability - prevent exponential overflow and division by zero
        double beta_arg = m_wall_beta * (distance - m_wall_spheric_radius);
        double exp_expr;
        if (beta_arg > 700.0) { // Prevent overflow
            exp_expr = std::numeric_limits<double>::max() / 2.0;
        } else if (beta_arg < -700.0) { // Prevent underflow
            exp_expr = 0.0;
        } else {
            exp_expr = exp(beta_arg);
        }
        double curr_pot = kbT * log(1 + exp_expr);
        // counter += distance > m_wall_radius;
        // std::cout << m_wall_beta*m_eigen_geometry.data()[3 * i + 0]*exp_expr/(distance*(1-exp_expr)) << " ";
        // Claude Generated 2026: fx/fy/fz are dV/dr (the gradient of V = kbT·log(1+exp(β(s-R)))
        // projected radially: dV/dx = kbT·β·exp/(1+exp)·(x/s)). m_eigen_gradient holds dE/dr
        // (force = -gradient), so the wall gradient must be ADDED. The previous `gradient -= fx`
        // subtracted it — flipping the wall force outward (atoms outside the sphere were expelled
        // instead of confined). Mirrors ApplyRectLogFermiWalls, which already adds dV/dr.
        // Add numerical stability check for distance = 0
        if (distance > 1e-10) {
            double fx = kbT * m_wall_beta * m_eigen_geometry.data()[3 * i + 0] * exp_expr / (distance * (1 + exp_expr));
            double fy = kbT * m_wall_beta * m_eigen_geometry.data()[3 * i + 1] * exp_expr / (distance * (1 + exp_expr));
            double fz = kbT * m_wall_beta * m_eigen_geometry.data()[3 * i + 2] * exp_expr / (distance * (1 + exp_expr));

            m_eigen_gradient.data()[3 * i + 0] += fx;
            m_eigen_gradient.data()[3 * i + 1] += fy;
            m_eigen_gradient.data()[3 * i + 2] += fz;

            // Track wall force magnitude
            sum_grad += std::sqrt(fx * fx + fy * fy + fz * fz);
        }

        // Count atoms outside sphere
        if (distance > m_wall_spheric_radius)
            counter++;

        // std::cout << distance << " ";
        potential += curr_pot;
    }

    // Claude Generated: Smart wall violation reporting - prevent console spam
    m_wall_violation_count = counter;

    // Only report if violations exceed 5% of atoms OR it's been 1000 steps since last report
    bool should_report = (counter > m_natoms * 0.05) || (counter > 0 && (m_currentStep - m_wall_violation_last_reported) > 1000) || (sum_grad > 0.01); // Or if wall forces are very high

    if (should_report) {
        std::cout << "Wall stats - Atoms outside sphere: " << counter << "/" << m_natoms
                  << ", Total wall force: " << sum_grad * au2N << " N"
                  << ", Wall potential: " << potential * au2eV << " eV" << std::endl;
        m_wall_violation_last_reported = m_currentStep;
    }

    return potential;
    // std::cout << potential*kbT << std::endl;
}

double SimpleMD::ApplyRectLogFermiWalls()
{
    double potential = 0;
    double kbT = m_wall_temp * kb_Eh;
    int counter = 0;
    double b = m_wall_beta;
    double sum_grad = 0;
    for (int i = 0; i < m_natoms; ++i) {
        double exp_expr_xl = exp(b * (m_wall_x_min - m_eigen_geometry.data()[3 * i + 0]));
        double exp_expr_xu = exp(b * (m_eigen_geometry.data()[3 * i + 0] - m_wall_x_max));

        double exp_expr_yl = exp(b * (m_wall_y_min - m_eigen_geometry.data()[3 * i + 1]));
        double exp_expr_yu = exp(b * (m_eigen_geometry.data()[3 * i + 1] - m_wall_y_max));

        double exp_expr_zl = exp(b * (m_wall_z_min - m_eigen_geometry.data()[3 * i + 2]));
        double exp_expr_zu = exp(b * (m_eigen_geometry.data()[3 * i + 2] - m_wall_z_max));

        double curr_pot = kbT * (log(1 + exp_expr_xl) + log(1 + exp_expr_xu) + log(1 + exp_expr_yl) + log(1 + exp_expr_yu) + log(1 + exp_expr_zl) + log(1 + exp_expr_zu));
        counter += (m_eigen_geometry.data()[3 * i + 0] - m_wall_x_min) < 0 || (m_wall_x_max - m_eigen_geometry.data()[3 * i + 0]) < 0 || (m_eigen_geometry.data()[3 * i + 1] - m_wall_y_min) < 0 || (m_wall_y_max - m_eigen_geometry.data()[3 * i + 1]) < 0 || (m_eigen_geometry.data()[3 * i + 2] - m_wall_z_min) < 0 || (m_wall_z_max - m_eigen_geometry.data()[3 * i + 2]) < 0;
        // std::cout << i << " " << counter << std::endl;

        // std::cout << m_wall_beta*m_eigen_geometry.data()[3 * i + 0]*exp_expr/(distance*(1-exp_expr)) << " ";
        if (i == 81) {
            //    std::cout << std::endl;
            //    std::cout << m_eigen_geometry.data()[3 * i + 0] << " " << m_eigen_geometry.data()[3 * i + 1] << " " << m_eigen_geometry.data()[3 * i + 2] << std::endl;
            //    std::cout << m_eigen_gradient.data()[3 * i + 0] << " " << m_eigen_gradient.data()[3 * i + 1] << " " <<m_eigen_gradient.data()[3 * i + 2] << std::endl;
        }
        // Claude Generated: Fix rectangular log-Fermi forces - correct denominator (1 + exp) for derivative
        double fx = kbT * b * (exp_expr_xu / (1 + exp_expr_xu) - exp_expr_xl / (1 + exp_expr_xl));
        double fy = kbT * b * (exp_expr_yu / (1 + exp_expr_yu) - exp_expr_yl / (1 + exp_expr_yl));
        double fz = kbT * b * (exp_expr_zu / (1 + exp_expr_zu) - exp_expr_zl / (1 + exp_expr_zl));

        m_eigen_gradient.data()[3 * i + 0] += fx;
        m_eigen_gradient.data()[3 * i + 1] += fy;
        m_eigen_gradient.data()[3 * i + 2] += fz;

        // Track wall force magnitude
        sum_grad += std::abs(fx) + std::abs(fy) + std::abs(fz);
        // if( i == 81)
        {
            // std::cout << i << " " <<m_eigen_gradient.data()[3 * i + 0] << " " << m_eigen_gradient.data()[3 * i + 1] << " " <<m_eigen_gradient.data()[3 * i + 2] << std::endl;
        }
        // std::cout << distance << " ";
        potential += curr_pot;
    }

    // Claude Generated: Smart wall violation reporting - prevent console spam
    m_wall_violation_count = counter;

    // Only report if violations exceed 5% of atoms OR it's been 1000 steps since last report
    bool should_report = (counter > m_natoms * 0.05) || (counter > 0 && (m_currentStep - m_wall_violation_last_reported) > 1000) || (sum_grad > 0.01); // Or if wall forces are very high

    if (should_report) {
        std::cout << "Wall stats - Atoms outside rectangular: " << counter << "/" << m_natoms
                  << ", Total wall force: " << sum_grad * au2N << " N"
                  << ", Wall potential: " << potential * au2eV << " eV" << std::endl;
        m_wall_violation_last_reported = m_currentStep;
    }

    return potential;
    // std::cout << potential*kbT << std::endl;
}

double SimpleMD::ApplySphericHarmonicWalls()
{
    double potential = 0;
    double k = m_wall_temp * kb_Eh;
    int counter = 0;
    double sum_grad = 0; // Claude Generated: Track total wall force
    for (int i = 0; i < m_natoms; ++i) {
        double distance = sqrt(m_eigen_geometry.data()[3 * i + 0] * m_eigen_geometry.data()[3 * i + 0] + m_eigen_geometry.data()[3 * i + 1] * m_eigen_geometry.data()[3 * i + 1] + m_eigen_geometry.data()[3 * i + 2] * m_eigen_geometry.data()[3 * i + 2]);
        double curr_pot = 0.5 * k * (m_wall_spheric_radius - distance) * (m_wall_spheric_radius - distance) * (distance > m_wall_spheric_radius);
        double out = distance > m_wall_spheric_radius;
        counter += out;

        double diff = k * (m_wall_spheric_radius - distance) * (distance > m_wall_spheric_radius);

        double dx = diff * m_eigen_geometry.data()[3 * i + 0] / distance;
        double dy = diff * m_eigen_geometry.data()[3 * i + 1] / distance;
        double dz = diff * m_eigen_geometry.data()[3 * i + 2] / distance;

        m_eigen_gradient.data()[3 * i + 0] -= dx;
        m_eigen_gradient.data()[3 * i + 1] -= dy;
        m_eigen_gradient.data()[3 * i + 2] -= dz;

        // Claude Generated: Track wall force magnitude
        sum_grad += std::sqrt(dx * dx + dy * dy + dz * dz);

        /*
        if(out)
        {
            std::cout << m_eigen_geometry.data()[3 * i + 0]  << " " << m_eigen_geometry.data()[3 * i + 1]  << " " << m_eigen_geometry.data()[3 * i + 2] << std::endl;
            std::cout << dx << " " << dy << " " << dz << std::endl;
        }*/
        // std::cout << distance << " ";
        potential += curr_pot;
    }

    // Claude Generated: Smart wall violation reporting - prevent console spam
    m_wall_violation_count = counter;

    // Only report if violations exceed 5% of atoms OR it's been 1000 steps since last report
    bool should_report = (counter > m_natoms * 0.05) || (counter > 0 && (m_currentStep - m_wall_violation_last_reported) > 1000) || (sum_grad > 0.01); // Or if wall forces are very high

    if (should_report) {
        std::cout << "Wall stats - Atoms outside sphere: " << counter << "/" << m_natoms
                  << ", Total wall force: " << sum_grad * au2N << " N"
                  << ", Wall potential: " << potential * au2eV << " eV" << std::endl;
        m_wall_violation_last_reported = m_currentStep;
    }

    return potential;
    // std::cout << potential*kbT << std::endl;
}

double SimpleMD::ApplyRectHarmonicWalls()
{
    double potential = 0;
    double k = m_wall_temp * kb_Eh;
    int counter = 0;
    double b = m_wall_beta;
    double sum_grad = 0;
    for (int i = 0; i < m_natoms; ++i) {
        double Vx = (m_eigen_geometry.data()[3 * i + 0] - m_wall_x_min) * (m_eigen_geometry.data()[3 * i + 0] - m_wall_x_min) * (m_eigen_geometry.data()[3 * i + 0] < m_wall_x_min)
            + (m_eigen_geometry.data()[3 * i + 0] - m_wall_x_max) * (m_eigen_geometry.data()[3 * i + 0] - m_wall_x_max) * (m_eigen_geometry.data()[3 * i + 0] > m_wall_x_max);

        double Vy = (m_eigen_geometry.data()[3 * i + 1] - m_wall_y_min) * (m_eigen_geometry.data()[3 * i + 1] - m_wall_y_min) * (m_eigen_geometry.data()[3 * i + 1] < m_wall_y_min)
            + (m_eigen_geometry.data()[3 * i + 1] - m_wall_y_max) * (m_eigen_geometry.data()[3 * i + 1] - m_wall_y_max) * (m_eigen_geometry.data()[3 * i + 1] > m_wall_y_max);

        double Vz = (m_eigen_geometry.data()[3 * i + 2] - m_wall_z_min) * (m_eigen_geometry.data()[3 * i + 2] - m_wall_z_min) * (m_eigen_geometry.data()[3 * i + 2] < m_wall_z_min)
            + (m_eigen_geometry.data()[3 * i + 2] - m_wall_z_max) * (m_eigen_geometry.data()[3 * i + 2] - m_wall_z_max) * (m_eigen_geometry.data()[3 * i + 2] > m_wall_z_max);

        double curr_pot = 0.5 * k * (Vx + Vy + Vz);
        int out = (m_eigen_geometry.data()[3 * i + 0] - m_wall_x_min) < 0 || (m_wall_x_max - m_eigen_geometry.data()[3 * i + 0]) < 0 || (m_eigen_geometry.data()[3 * i + 1] - m_wall_y_min) < 0 || (m_wall_y_max - m_eigen_geometry.data()[3 * i + 1]) < 0 || (m_eigen_geometry.data()[3 * i + 2] - m_wall_z_min) < 0 || (m_wall_z_max - m_eigen_geometry.data()[3 * i + 2]) < 0;
        counter += out;

        // std::cout << i << " " << counter << std::endl;

        // Claude Generated 2026: Correct harmonic wall gradient.
        // m_eigen_gradient holds dE/dr (the integrator does v -= ½·dT·grad/m,
        // i.e. force = -gradient). The wall adds V = ½k·d² to the energy, so its
        // gradient contribution dV/dr = k·((r-r_min)·(r<r_min) + (r-r_max)·(r>r_max))
        // must be ADDED. The previous form used `gradient -= dx` with a minus between
        // the min/max terms: that added the max-wall gradient (correct) but SUBTRACTED
        // the min-wall gradient (sign error — atoms below r_min were pushed further out
        // instead of back in). Symmetric `+` with `gradient +=` fixes both walls.
        double gx = m_eigen_geometry.data()[3 * i + 0];
        double gy = m_eigen_geometry.data()[3 * i + 1];
        double gz = m_eigen_geometry.data()[3 * i + 2];
        double dx = k * ((gx - m_wall_x_min) * (gx < m_wall_x_min) + (gx - m_wall_x_max) * (gx > m_wall_x_max));
        double dy = k * ((gy - m_wall_y_min) * (gy < m_wall_y_min) + (gy - m_wall_y_max) * (gy > m_wall_y_max));
        double dz = k * ((gz - m_wall_z_min) * (gz < m_wall_z_min) + (gz - m_wall_z_max) * (gz > m_wall_z_max));
        m_eigen_gradient.data()[3 * i + 0] += dx;
        m_eigen_gradient.data()[3 * i + 1] += dy;
        m_eigen_gradient.data()[3 * i + 2] += dz;
        /* if(out)
         {
             std::cout << m_eigen_geometry.data()[3 * i + 0]  << " " << m_eigen_geometry.data()[3 * i + 1]  << " " << m_eigen_geometry.data()[3 * i + 2] << std::endl;
             std::cout << dx << " " << dy << " " << dz << std::endl;
         }*/
        sum_grad += std::abs(dx) + std::abs(dy) + std::abs(dz);

        potential += curr_pot;
    }

    // Claude Generated: Smart wall violation reporting - prevent console spam
    m_wall_violation_count = counter;

    // Only report if violations exceed 5% of atoms OR it's been 1000 steps since last report
    bool should_report = (counter > m_natoms * 0.05) || (counter > 0 && (m_currentStep - m_wall_violation_last_reported) > 1000) || (sum_grad > 0.01); // Or if wall forces are very high

    if (should_report) {
        std::cout << "Wall stats - Atoms outside rectangular: " << counter << "/" << m_natoms
                  << ", Total wall force: " << sum_grad * au2N << " N"
                  << ", Wall potential: " << potential * au2eV << " eV" << std::endl;
        m_wall_violation_last_reported = m_currentStep;
    }

    return potential;
    // std::cout << potential*kbT << std::endl;
}

void SimpleMD::RemoveRotations()
{
    /*
     * This code was taken and adopted from the xtb sources
     * https://github.com/grimme-lab/xtb/blob/main/src/rmrottr.f90
     * Special thanks to the developers
     */
    Geometry geom(m_natoms, 3);

    std::vector<std::vector<int>> fragments = m_molecule.GetFragments();
    for (auto & fragment : fragments) {
        // Reset per-fragment accumulators so each fragment's COM and angular
        // momentum are computed independently (pre-existing bug: these were
        // declared outside the loop and accumulated across fragments).
        double mass = 0;
        Position pos = { 0, 0, 0 }, angom{ 0, 0, 0 };
        for (const int i : fragment) {
            const double m = m_eigen_masses.data()[3 * i];
            mass += m;
            pos(0) += m * m_eigen_geometry.data()[3 * i + 0];
            pos(1) += m * m_eigen_geometry.data()[3 * i + 1];
            pos(2) += m * m_eigen_geometry.data()[3 * i + 2];

            geom(i, 0) = m_eigen_geometry.data()[3 * i + 0];
            geom(i, 1) = m_eigen_geometry.data()[3 * i + 1];
            geom(i, 2) = m_eigen_geometry.data()[3 * i + 2];
        }
        pos(0) /= mass;
        pos(1) /= mass;
        pos(2) /= mass;

        Geometry matrix = Geometry::Zero(3, 3);
        for (const int i : fragment) {
            const double m = m_eigen_masses.data()[3 * i];
            geom(i, 0) -= pos(0);
            geom(i, 1) -= pos(1);
            geom(i, 2) -= pos(2);

            const double x = geom(i, 0);
            const double y = geom(i, 1);
            const double z = geom(i, 2);
            angom(0) += m_eigen_masses.data()[3 * i] * (geom(i, 1) *  m_eigen_velocities.data()[3 * i + 2] - geom(i, 2) *  m_eigen_velocities.data()[3 * i + 1]);
            angom(1) += m_eigen_masses.data()[3 * i] * (geom(i, 2) *  m_eigen_velocities.data()[3 * i + 0] - geom(i, 0) *  m_eigen_velocities.data()[3 * i + 2]);
            angom(2) += m_eigen_masses.data()[3 * i] * (geom(i, 0) *  m_eigen_velocities.data()[3 * i + 1] - geom(i, 1) *  m_eigen_velocities.data()[3 * i + 0]);
            const double x2 = x * x;
            const double y2 = y * y;
            const double z2 = z * z;
            matrix(0, 0) += m * (y2 + z2);
            matrix(1, 1) += m * (x2 + z2);
            matrix(2, 2) += m * (x2 + y2);
            matrix(0, 1) -= m * x * y;
            matrix(0, 2) -= m * x * z;
            matrix(1, 2) -= m * y * z;
        }
        matrix(1, 0) = matrix(0, 1);
        matrix(2, 0) = matrix(0, 2);
        matrix(2, 1) = matrix(1, 2);

        // Robust solve for singular/near-singular inertia tensors (linear
        // molecules, single atoms). A naive matrix.inverse() returns NaN/Inf
        // when one or more principal moments of inertia vanish, which is
        // exactly the case for 2-atom systems (rotation about the bond axis
        // has zero moment). Use a JacobiSVD and clamp small singular values
        // so that the corresponding rotational DOF is correctly frozen
        // instead of producing NaN velocities.
        // Guard checks fragment size, not total atom count: a single-atom
        // fragment inside a larger system must also be skipped.
        Position omega = { 0, 0, 0 };
        if (fragment.size() > 1) {
            Eigen::JacobiSVD<Geometry> svd(matrix, Eigen::ComputeFullU | Eigen::ComputeFullV);
            const double sv0 = svd.singularValues()(0);
            const double sv_tol = 1e-12 * std::max(sv0, 1e-300);
            Eigen::Vector3d inv_sv = Eigen::Vector3d::Zero();
            for (int k = 0; k < 3; ++k) {
                const double sv = svd.singularValues()(k);
                inv_sv(k) = (sv > sv_tol) ? 1.0 / sv : 0.0;
            }
            omega = svd.matrixU() * inv_sv.asDiagonal() * svd.matrixV().transpose() * angom;
        }

        Position rlm = { 0, 0, 0 }, ram = { 0, 0, 0 };
        for (const int i : fragment) {
            rlm(0) = rlm(0) + m_eigen_masses.data()[3 * i] *  m_eigen_velocities.data()[3 * i + 0];
            rlm(1) = rlm(1) + m_eigen_masses.data()[3 * i] *  m_eigen_velocities.data()[3 * i + 1];
            rlm(2) = rlm(2) + m_eigen_masses.data()[3 * i] *  m_eigen_velocities.data()[3 * i + 2];
        }

        for (const int i : fragment) {
            ram(0) = (omega(1) * geom(i, 2) - omega(2) * geom(i, 1));
            ram(1) = (omega(2) * geom(i, 0) - omega(0) * geom(i, 2));
            ram(2) = (omega(0) * geom(i, 1) - omega(1) * geom(i, 0));

             m_eigen_velocities.data()[3 * i + 0] = m_eigen_velocities.data()[3 * i + 0] - rlm(0) / mass - ram(0);
             m_eigen_velocities.data()[3 * i + 1] = m_eigen_velocities.data()[3 * i + 1] - rlm(1) / mass - ram(1);
             m_eigen_velocities.data()[3 * i + 2] = m_eigen_velocities.data()[3 * i + 2] - rlm(2) / mass - ram(2);
        }
    }
}

void SimpleMD::RemoveRotation()
{
    /*
     * This code was taken and adopted from the xtb sources
     * https://github.com/grimme-lab/xtb/blob/main/src/rmrottr.f90
     * Special thanks to the developers
     */
    double mass = 0;
    Position pos = { 0, 0, 0 }, angom{ 0, 0, 0 };
    Geometry geom(m_natoms, 3);

    for (int i = 0; i < m_natoms; ++i) {
        double m = m_eigen_masses.data()[3 * i];
        mass += m;
        pos(0) += m * m_eigen_geometry.data()[3 * i + 0];
        pos(1) += m * m_eigen_geometry.data()[3 * i + 1];
        pos(2) += m * m_eigen_geometry.data()[3 * i + 2];

        geom(i, 0) = m_eigen_geometry.data()[3 * i + 0];
        geom(i, 1) = m_eigen_geometry.data()[3 * i + 1];
        geom(i, 2) = m_eigen_geometry.data()[3 * i + 2];
    }
    pos(0) /= mass;
    pos(1) /= mass;
    pos(2) /= mass;

    Geometry matrix = Geometry::Zero(3, 3);
    for (int i = 0; i < m_natoms; ++i) {
        double m = m_eigen_masses.data()[3 * i];
        geom(i, 0) -= pos(0);
        geom(i, 1) -= pos(1);
        geom(i, 2) -= pos(2);

        double x = geom(i, 0);
        double y = geom(i, 1);
        double z = geom(i, 2);
        angom(0) += m_eigen_masses.data()[3 * i] * (geom(i, 1) * m_eigen_velocities.data()[3 * i + 2] - geom(i, 2) *  m_eigen_velocities.data()[3 * i + 1]);
        angom(1) += m_eigen_masses.data()[3 * i] * (geom(i, 2) * m_eigen_velocities.data()[3 * i + 0] - geom(i, 0) *  m_eigen_velocities.data()[3 * i + 2]);
        angom(2) += m_eigen_masses.data()[3 * i] * (geom(i, 0) * m_eigen_velocities.data()[3 * i + 1] - geom(i, 1) *  m_eigen_velocities.data()[3 * i + 0]);
        double x2 = x * x;
        double y2 = y * y;
        double z2 = z * z;
        matrix(0, 0) += m * (y2 + z2);
        matrix(1, 1) += m * (x2 + z2);
        matrix(2, 2) += m * (x2 + y2);
        matrix(0, 1) -= m * x * y;
        matrix(0, 2) -= m * x * z;
        matrix(1, 2) -= m * y * z;
    }
    matrix(1, 0) = matrix(0, 1);
    matrix(2, 0) = matrix(0, 2);
    matrix(2, 1) = matrix(1, 2);

    // Robust solve for singular/near-singular inertia tensors (linear
    // molecules, single atoms). See SimpleMD::RemoveRotations above for the
    // rationale. Required so that 2-atom systems do not get NaN velocities
    // at step 0 and crash the integrator with "NaN/Inf velocity".
    Position omega = { 0, 0, 0 };
    if (m_natoms > 1) {
        Eigen::JacobiSVD<Geometry> svd(matrix, Eigen::ComputeFullU | Eigen::ComputeFullV);
        const double sv0 = svd.singularValues()(0);
        const double sv_tol = 1e-12 * std::max(sv0, 1e-300);
        Eigen::Vector3d inv_sv = Eigen::Vector3d::Zero();
        for (int k = 0; k < 3; ++k) {
            const double sv = svd.singularValues()(k);
            inv_sv(k) = (sv > sv_tol) ? 1.0 / sv : 0.0;
        }
        omega = svd.matrixU() * inv_sv.asDiagonal() * svd.matrixV().transpose() * angom;
    }

    Position rlm = { 0, 0, 0 }, ram = { 0, 0, 0 };
    for (int i = 0; i < m_natoms; ++i) {
        rlm(0) = rlm(0) + m_eigen_masses.data()[3 * i] *  m_eigen_velocities.data()[3 * i + 0];
        rlm(1) = rlm(1) + m_eigen_masses.data()[3 * i] *  m_eigen_velocities.data()[3 * i + 1];
        rlm(2) = rlm(2) + m_eigen_masses.data()[3 * i] *  m_eigen_velocities.data()[3 * i + 2];
    }

    for (int i = 0; i < m_natoms; ++i) {
        ram(0) = (omega(1) * geom(i, 2) - omega(2) * geom(i, 1));
        ram(1) = (omega(2) * geom(i, 0) - omega(0) * geom(i, 2));
        ram(2) = (omega(0) * geom(i, 1) - omega(1) * geom(i, 0));

         m_eigen_velocities.data()[3 * i + 0] =  m_eigen_velocities.data()[3 * i + 0] - rlm(0) / mass - ram(0);
         m_eigen_velocities.data()[3 * i + 1] =  m_eigen_velocities.data()[3 * i + 1] - rlm(1) / mass - ram(1);
         m_eigen_velocities.data()[3 * i + 2] =  m_eigen_velocities.data()[3 * i + 2] - rlm(2) / mass - ram(2);
    }
}

void SimpleMD::PrintStatus() const
{
    // Header is printed once before the MD loop starts (not here).

    const auto unix_timestamp = std::chrono::seconds(std::time(NULL));

    const int current = std::chrono::milliseconds(unix_timestamp).count();
    const double duration = (current - m_unix_started) / (1000.0 * static_cast<double>(m_currentStep));
    double remaining;
    if (const double tmp = (m_maxtime - m_currentStep) * duration / 60; tmp >= 1)
        remaining = tmp;
    else
        remaining = (m_maxtime - m_currentStep) * duration;
#pragma message("awfull, fix it ")
    // Claude Generated (May 2026): unified MD step printout — matches the header above.
    // Base 15 columns always present (Time, energies, T, wall, virial, remaining, dt).
    // Optional appendix columns: Dipole (when -dipole), nUnique (when -writeUnique).
    // Replaces three divergent layouts (writeUnique with nUnique-at-end, dipole with
    // dipole-mid-row, and a 5-col tab-separated fallback). The dropped #ifdef GCC
    // pre-processor guard was dead — no compiler defines `GCC`.
    {
        std::string line = fmt::format(
            "{1: ^{0}f} {2: ^{0}f} {3: ^{0}f} {4: ^{0}f} {5: ^{0}f} {6: ^{0}f} {7: ^{0}f} "
            "{8: ^{0}f} {9: ^{0}f} {10: ^{0}f} {11: ^{0}f} {12: ^{0}f} {13: ^{0}f} {14: ^{0}f} "
            "{15: ^{0}f}",
            15,
            m_currentStep / 1000.0,  // Time in ps (m_currentStep accumulates m_dT, both in fs)
            m_Epot, m_aver_Epot, m_Ekin, m_aver_Ekin, m_Etot, m_aver_Etot,
            m_T, m_aver_Temp, m_wall_potential, m_average_wall_potential,
            m_virial_correction, m_average_virial_correction,
            remaining, m_time_step / 1000.0);
        if (m_dipole)
            line += fmt::format(" {: ^15f}", m_aver_dipol_linear * 2.5418 * 3.3356);
        if (m_rmsd_mtd)
            line += fmt::format(" {: ^15}", m_bias_structure_count);
        else if (m_writeUnique)
            line += fmt::format(" {: ^15}", m_unqiue->StoredStructures());
        std::cout << line << "\n";
    }

    // RATTLE constraint summary (only when RATTLE is active)
    if (m_rattle) {
        CurcumaLogger::info(fmt::format("RATTLE: iters={}/{} max_err_12={:.2e} max_err_13={:.2e} failures={}",
            m_rattle_iters_step1, m_rattle_iters_step2,
            m_rattle_max_err_12, m_rattle_max_err_13, m_rattle_max_err_count));
    }

    //std::cout << m_mtd_time << " " << m_loop_time << std::endl;
}

void SimpleMD::PrintMatrix(const double* matrix) const
{
    std::cout << "Print Matrix" << std::endl;
    for (int i = 0; i < m_natoms; ++i) {
        std::cout << matrix[3 * i] << " " << matrix[3 * i + 1] << " " << matrix[3 * i + 2] << std::endl;
    }
    std::cout << std::endl;
}

double SimpleMD::CleanEnergy()
{
    // Claude Generated: Use new constructor with basename for parameter caching
    EnergyCalculator interface(m_method, m_defaults, Basename());
    interface.setMolecule(m_molecule.getMolInfo());
    interface.updateGeometry(m_eigen_geometry);
    interface.setIterativeMode(true);
    // WP-P1 (May 2026): wall-clock the FF call for MD diagnostics
    auto t_ff_start = std::chrono::high_resolution_clock::now();
    const double Energy = interface.CalculateEnergy(true);
    auto t_ff_end = std::chrono::high_resolution_clock::now();
    m_last_ff_ms = std::chrono::duration<double, std::milli>(t_ff_end - t_ff_start).count();
    m_eigen_gradient = interface.Gradient();
    if (m_dipole && m_method == "gfn2") {
        m_molecule.setDipole(interface.Dipole()*au);//in eA
        m_molecule.setPartialCharges(interface.Charges());
    }
    return Energy;
}

double SimpleMD::FastEnergy()
{
    m_interface->updateGeometry(m_eigen_geometry);

    // WP-P1 (May 2026): wall-clock the FF call for MD diagnostics
    auto t_ff_start = std::chrono::high_resolution_clock::now();
    const double Energy = m_interface->CalculateEnergy(true);
    auto t_ff_end = std::chrono::high_resolution_clock::now();
    m_last_ff_ms = std::chrono::duration<double, std::milli>(t_ff_end - t_ff_start).count();
    m_eigen_gradient = m_interface->Gradient();

    // Claude Generated (Feb 2026): Gradient sanity check for MD stability
    for (int i = 0; i < 3 * m_natoms; ++i) {
        if (!std::isfinite(m_eigen_gradient.data()[i])) {
            CurcumaLogger::error("NaN/Inf gradient detected in MD step - simulation unstable");
            m_unstable = true;
            return Energy;
        }
    }

    if (m_dipole && m_method == "gfn2") {
        m_molecule.setDipole(m_interface->Dipole()*au);// in eA
        m_molecule.setPartialCharges(m_interface->Charges());
    }
    return Energy;
}

void SimpleMD::EKin()
{
    double ekin = 0;
    for (int i = 0; i < m_natoms; ++i) {
        ekin += m_eigen_masses.data()[3 * i] * (m_eigen_velocities.data()[3 * i] * m_eigen_velocities.data()[3 * i] + m_eigen_velocities.data()[3 * i + 1] * m_eigen_velocities.data()[3 * i + 1] + m_eigen_velocities.data()[3 * i + 2] * m_eigen_velocities.data()[3 * i + 2]);
    }
    ekin *= 0.5;
    m_Ekin = ekin;
    m_T = 2.0 * ekin / (kb_Eh * m_dof);
}

void SimpleMD::AverageQuantities()
{
    m_aver_Temp = (m_T + (m_currentStep)*m_aver_Temp) / (m_currentStep + 1);
    m_aver_Epot = (m_Epot + (m_currentStep)*m_aver_Epot) / (m_currentStep + 1);
    m_aver_Ekin = (m_Ekin + (m_currentStep)*m_aver_Ekin) / (m_currentStep + 1);
    m_aver_Etot = (m_Etot + (m_currentStep)*m_aver_Etot) / (m_currentStep + 1);
    if (m_dipole) {
        //m_aver_dipol = (m_curr_dipoles + (m_currentStep)*m_aver_dipol) / (m_currentStep + 1);
    }
    m_average_wall_potential = (m_wall_potential + (m_currentStep)*m_average_wall_potential) / (m_currentStep + 1);
    m_average_virial_correction = (m_virial_correction + (m_currentStep)*m_average_virial_correction) / (m_currentStep + 1);
}

bool SimpleMD::WriteGeometry()
{
    bool result = true;
    Geometry geometry = m_molecule.getGeometry();
    for (int i = 0; i < m_natoms; ++i) {
        geometry(i, 0) = m_eigen_geometry.data()[3 * i + 0];
        geometry(i, 1) = m_eigen_geometry.data()[3 * i + 1];
        geometry(i, 2) = m_eigen_geometry.data()[3 * i + 2];
    }
    TriggerWriteRestart();
    m_molecule.setGeometry(geometry);

    // Claude Generated (Nov 2025): Write VTF trajectory for CG systems
    if (m_is_cg_system && m_cg_write_vtf) {
        m_molecule.setEnergy(m_Epot);
        m_molecule.setName(std::to_string(m_currentStep));
        m_molecule.appendVTFFile(outputPath(Basename() + ".trj.vtf"));
    }

    if (m_writeXYZ) {
        m_molecule.setEnergy(m_Epot);
        m_molecule.setName(std::to_string(m_currentStep));
        m_molecule.appendXYZFile(outputPath(Basename() + ".trj.xyz"));
    }
    if (m_writeUnique) {
        if (m_unqiue->CheckMolecule(new Molecule(m_molecule))) {
            std::cout << " ** new structure was added **" << std::endl;
            PrintStatus();
            m_time_step = 0;
            m_unique_structures.push_back(new Molecule(m_molecule));
        }
    }
    return result;
}

void SimpleMD::None()
{
}

void SimpleMD::Berendson()
{
    double lambda = sqrt(1 + (m_dT / 2.0 * (m_T0 - m_T)) / (m_T * m_coupling));
    for (int i = 0; i < m_natoms; ++i) {
        m_eigen_velocities.data()[3 * i + 0] *= lambda;
        m_eigen_velocities.data()[3 * i + 1] *= lambda;
        m_eigen_velocities.data()[3 * i + 2] *= lambda;
    }
}

void SimpleMD::CSVR()
{
    double Ekin_target = 0.5 * kb_Eh * (m_T0)*m_dof;
    double c = exp(-(m_dT / 2.0 * m_respa) / m_coupling);
    static std::default_random_engine rd{};
    static std::mt19937 gen{ rd() };
    static std::normal_distribution<> d{ 0, 1 };
    // Lazy-reinit when m_dof changes (e.g. after RATTLE constraint setup or
    // first call with a different molecule). A stale static distribution with
    // the wrong DOF produces wrong fluctuation widths in the CSVR rescaling.
    static int csvr_last_dof = -1;
    static std::chi_squared_distribution<float> dchi{ 1.0f };
    if (m_dof != csvr_last_dof) {
        dchi = std::chi_squared_distribution<float>(static_cast<float>(m_dof));
        csvr_last_dof = m_dof;
    }
    double R = d(gen);
    double SNf = dchi(gen);
    double alpha2 = c + (1 - c) * (SNf + R * R) * Ekin_target / (m_dof * m_Ekin) + 2 * R * sqrt(c * (1 - c) * Ekin_target / (m_dof * m_Ekin));
    m_Ekin_exchange += m_Ekin * (alpha2 - 1);
    double alpha = sqrt(alpha2);
    for (int i = 0; i < m_natoms; ++i) {
        m_eigen_velocities.data()[3 * i + 0] *= alpha;
        m_eigen_velocities.data()[3 * i + 1] *= alpha;
        m_eigen_velocities.data()[3 * i + 2] *= alpha;

        m_atom_temp[i].push_back(m_eigen_masses.data()[3 * i] * (m_eigen_velocities.data()[3 * i + 0] * m_eigen_velocities.data()[3 * i + 0] + m_eigen_velocities.data()[3 * i + 1] * m_eigen_velocities.data()[3 * i + 1] + m_eigen_velocities.data()[3 * i + 2] * m_eigen_velocities.data()[3 * i + 2]) / (kb_Eh * m_dof));
    }
    m_seed++;
}

void SimpleMD::Anderson()
{
    static std::default_random_engine generator;
    double probability = m_anderson * m_dT;
    std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);
    for (size_t i = 0; i < m_natoms; ++i) {
        if (uniform_dist(generator) < probability) {
            std::normal_distribution<double> distribution(0.0, std::sqrt(kb_Eh * m_T0 * m_eigen_inv_masses.data()[3 * i]));
            m_eigen_velocities.data()[3 * i + 0] = (m_eigen_velocities.data()[3 * i + 0] + distribution(generator)) / 2.0;
            m_eigen_velocities.data()[3 * i + 1] = (m_eigen_velocities.data()[3 * i + 1] + distribution(generator)) / 2.0;
            m_eigen_velocities.data()[3 * i + 2] = (m_eigen_velocities.data()[3 * i + 2] + distribution(generator)) / 2.0;
            m_seed += 3;
        }
    }
}
void SimpleMD::NoseHover()
{
    // Berechnung der kinetischen Energie
    double kinetic_energy = 0.0;
    for (int i = 0; i < m_natoms; ++i) {
        kinetic_energy += 0.5 * m_eigen_masses.data()[3 * i] * (m_eigen_velocities.data()[3 * i] * m_eigen_velocities.data()[3 * i] + m_eigen_velocities.data()[3 * i + 1] * m_eigen_velocities.data()[3 * i + 1] + m_eigen_velocities.data()[3 * i + 2] * m_eigen_velocities.data()[3 * i + 2]);
    }
    // Update der Thermostatkette
    m_xi[0] += 0.5 * m_dT * (2.0 * kinetic_energy - m_dof * m_T0 * kb_Eh) / m_Q[0];
    for (int j = 1; j < m_chain_length; ++j) {
        m_xi[j] += 0.5 * m_dT * (m_Q[j - 1] * m_xi[j - 1] * m_xi[j - 1] - m_T0 * kb_Eh) / m_Q[j];
    }

    // Update der Geschwindigkeiten
    double scale = exp(-m_xi[0] * m_dT);
    for (int i = 0; i < m_natoms; ++i) {
        m_eigen_velocities.data()[3 * i + 0] *= scale;
        m_eigen_velocities.data()[3 * i + 1] *= scale;
        m_eigen_velocities.data()[3 * i + 2] *= scale;
    }

    // Rückwärts-Update der Thermostatkette
    for (int j = m_chain_length - 1; j >= 1; --j) {
        m_xi[j] += 0.5 * m_dT * (m_Q[j - 1] * m_xi[j - 1] * m_xi[j - 1] - m_T0 * kb_Eh) / m_Q[j];
    }
    m_xi[0] += 0.5 * m_dT * (2.0 * kinetic_energy - m_dof * m_T0 * kb_Eh) / m_Q[0];
}
