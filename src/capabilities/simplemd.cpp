/*
 * <Simple MD Module for Cucuma. >
 * Copyright (C) 2020 - 2024 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include <cstring>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <vector>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "src/capabilities/curcumaopt.h"
#include "src/capabilities/rmsd.h"
#include "src/capabilities/rmsdtraj.h"

#include "src/core/elements.h"
#include "src/core/energycalculator.h"
#include "src/core/fileiterator.h"
#include "src/core/global.h"
#include "src/core/molecule.h"

#include "src/tools/geometry.h"

#include "external/CxxThreadPool/include/CxxThreadPool.h"

#ifdef USE_Plumed
#include "plumed2/src/wrapper/Plumed.h"
#endif
#include "simplemd.h"

BiasThread::BiasThread(const Molecule& reference, const json& rmsdconfig, bool nocolvarfile, bool nohillsfile)
    : m_reference(reference)
    , m_target(reference)
    , m_nocolvarfile(nocolvarfile)
    , m_nohillsfile(nohillsfile)
{
    m_driver = RMSDDriver(rmsdconfig, true);
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
    if (m_biased_structures.size() == 0)
        return 0;
    m_current_bias = 0;
    m_counter = 0;
    m_driver.setReference(m_reference);
    m_gradient = Eigen::MatrixXd::Zero(m_reference.AtomCount(), 3);

    for (int i = 0; i < m_biased_structures.size(); ++i) {
        double factor = 1;
        m_target.setGeometry(m_biased_structures[i].geometry);
        m_driver.setTarget(m_target);
        double rmsd = m_driver.BestFitRMSD();
        double expr = exp(-rmsd * rmsd * m_alpha);
        double bias_energy = expr * m_dT;
        factor = m_biased_structures[i].factor;

        if (!m_wtmtd)
            factor = m_biased_structures[i].counter;
        else
            factor += (exp(-(m_biased_structures[i].energy) / kb_Eh / m_DT));
        m_biased_structures[i].factor = factor;
        if (i == 0) {
            m_rmsd_reference = rmsd;
        }
        if (expr * m_rmsd_econv > 1 * m_biased_structures.size()) {
            m_biased_structures[i].counter++;
            m_biased_structures[i].energy += bias_energy;
        }
        bias_energy *= factor * m_k;

        m_current_bias += bias_energy;
        if (m_nocolvarfile == false) {
            std::ofstream colvarfile;
            colvarfile.open("COLVAR_" + std::to_string(m_biased_structures[i].index), std::iostream::app);
            colvarfile << m_currentStep << " " << rmsd << " " << bias_energy << " " << m_biased_structures[i].counter << " " << factor << std::endl;
            colvarfile.close();
        }
        /*
        if(nohillsfile == false)
        {
            std::ofstream hillsfile;
            if (i == 0) {
                hillsfile.open("HILLS", std::iostream::app);
            } else {
                hillsfile.open("HILLS_" + std::to_string(m_biased_structures[i].index), std::iostream::app);
            }
            hillsfile << m_currentStep << " " << rmsd << " " << m_alpha_rmsd << " " << m_k_rmsd << " " << "-1" << std::endl;
            hillsfile.close();
        }
        */

        double dEdR = -2 * m_alpha * m_k / m_atoms * exp(-rmsd * rmsd * m_alpha) * factor * m_dT;

        m_gradient += m_driver.Gradient() * dEdR;
        m_counter += m_biased_structures[i].counter;
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

SimpleMD::SimpleMD(const json& controller, bool silent)
    : CurcumaMethod(CurcumaMDJson, controller, silent)
{
    UpdateController(controller);
    m_interface = new EnergyCalculator(m_method, controller["md"]);
}

SimpleMD::~SimpleMD()
{
    for (int i = 0; i < m_unique_structures.size(); ++i)
        delete m_unique_structures[i];
    // delete m_bias_pool;
}

void SimpleMD::LoadControlJson()
{
    m_method = Json2KeyWord<std::string>(m_defaults, "method");
    m_thermostat = Json2KeyWord<std::string>(m_defaults, "thermostat");
    m_plumed = Json2KeyWord<std::string>(m_defaults, "plumed");

    m_spin = Json2KeyWord<int>(m_defaults, "spin");
    m_charge = Json2KeyWord<int>(m_defaults, "charge");
    m_dT = Json2KeyWord<double>(m_defaults, "dT");
    m_maxtime = Json2KeyWord<double>(m_defaults, "MaxTime");
    m_T0 = Json2KeyWord<double>(m_defaults, "T");
    m_rmrottrans = Json2KeyWord<int>(m_defaults, "rmrottrans");
    m_nocenter = Json2KeyWord<bool>(m_defaults, "nocenter");
    m_COM = Json2KeyWord<bool>(m_defaults, "COM");
    m_dump = Json2KeyWord<int>(m_defaults, "dump");
    m_print = Json2KeyWord<int>(m_defaults, "print");
    m_max_top_diff = Json2KeyWord<int>(m_defaults, "MaxTopoDiff");
    m_seed = Json2KeyWord<int>(m_defaults, "seed");
    m_threads = Json2KeyWord<double>(m_defaults, "threads");

    m_rmsd = Json2KeyWord<double>(m_defaults, "rmsd");
    m_hmass = Json2KeyWord<double>(m_defaults, "hmass");

    m_impuls = Json2KeyWord<double>(m_defaults, "impuls");
    m_impuls_scaling = Json2KeyWord<double>(m_defaults, "impuls_scaling");
    m_writeUnique = Json2KeyWord<bool>(m_defaults, "unique");
    m_opt = Json2KeyWord<bool>(m_defaults, "opt");
    m_scale_velo = Json2KeyWord<double>(m_defaults, "velo");
    m_rescue = Json2KeyWord<bool>(m_defaults, "rescue");
    m_wall_render = Json2KeyWord<bool>(m_defaults, "wall_render");
    m_coupling = Json2KeyWord<double>(m_defaults, "coupling");
    m_anderson = Json2KeyWord<double>(m_defaults, "anderson");

    m_threads = Json2KeyWord<double>(m_defaults, "threads");
    if (m_coupling < m_dT)
        m_coupling = m_dT;

    /* RMSD Metadynamik block */
    /* this one is used to recover https://doi.org/10.1021/acs.jctc.9b00143 */
    m_rmsd_mtd = Json2KeyWord<bool>(m_defaults, "rmsd_mtd");
    m_k_rmsd = Json2KeyWord<double>(m_defaults, "k_rmsd");
    m_alpha_rmsd = Json2KeyWord<double>(m_defaults, "alpha_rmsd");
    m_mtd_steps = Json2KeyWord<int>(m_defaults, "mtd_steps");
    m_chain_length = Json2KeyWord<int>(m_defaults, "chainlength");
    m_rmsd_rmsd = Json2KeyWord<double>(m_defaults, "rmsd_rmsd");
    m_max_rmsd_N = Json2KeyWord<int>(m_defaults, "max_rmsd_N");
    m_rmsd_econv = Json2KeyWord<double>(m_defaults, "rmsd_econv");
    m_rmsd_DT = Json2KeyWord<double>(m_defaults, "rmsd_DT");
    m_wtmtd = Json2KeyWord<bool>(m_defaults, "wtmtd");
    m_rmsd_ref_file = Json2KeyWord<std::string>(m_defaults, "rmsd_ref_file");
    m_rmsd_fix_structure = Json2KeyWord<bool>(m_defaults, "rmsd_fix_structure");
    m_nocolvarfile = Json2KeyWord<bool>(m_defaults, "noCOLVARfile");
    m_nohillsfile = Json2KeyWord<bool>(m_defaults, "noHILSfile");

    m_rmsd_atoms = Json2KeyWord<std::string>(m_defaults, "rmsd_atoms");

    m_writerestart = Json2KeyWord<int>(m_defaults, "writerestart");
    m_respa = Json2KeyWord<int>(m_defaults, "respa");
    m_dipole = Json2KeyWord<bool>(m_defaults, "dipole");

    m_writeXYZ = Json2KeyWord<bool>(m_defaults, "writeXYZ");
    m_writeinit = Json2KeyWord<bool>(m_defaults, "writeinit");
    m_mtd = Json2KeyWord<bool>(m_defaults, "mtd");
    m_mtd_dT = Json2KeyWord<int>(m_defaults, "mtd_dT");
    if (m_mtd_dT < 0) {
        m_eval_mtd = true;
    } else {
        m_eval_mtd = false;
    }
    m_initfile = Json2KeyWord<std::string>(m_defaults, "initfile");
    m_norestart = Json2KeyWord<bool>(m_defaults, "norestart");
    m_dt2 = m_dT * m_dT;
    m_rm_COM = Json2KeyWord<double>(m_defaults, "rm_COM");
    int rattle = Json2KeyWord<int>(m_defaults, "rattle");
    m_rattle_maxiter = Json2KeyWord<int>(m_defaults, "rattle_maxiter");
    m_rattle_dynamic_tol_iter = Json2KeyWord<int>(m_defaults, "rattle_dynamic_tol_iter");
    m_rattle_dynamic_tol = Json2KeyWord<bool>(m_defaults, "rattle_dynamic_tol");

    if (rattle == 1) {
        Integrator = [=](double* grad) {
            this->Rattle(grad);
        };
        m_rattle_tolerance = Json2KeyWord<double>(m_defaults, "rattle_tolerance");
        // m_coupling = m_dT;
        m_rattle = Json2KeyWord<int>(m_defaults, "rattle");
        std::cout << "Using rattle to constrain bonds!" << std::endl;
    } else {
        Integrator = [=](double* grad) {
            this->Verlet(grad);
        };
    }

    if (Json2KeyWord<bool>(m_defaults, "cleanenergy")) {
        Energy = [=](double* grad) -> double {
            return this->CleanEnergy(grad);
        };
        std::cout << "Energy Calculator will be set up for each step! Single steps are slower, but more reliable. Recommended for the combination of GFN2 and solvation." << std::endl;
    } else {
        Energy = [=](double* grad) -> double {
            return this->FastEnergy(grad);
        };
        std::cout << "Energy Calculator will NOT be set up for each step! Fast energy calculation! This is the default way and should not be changed unless the energy and gradient calculation are unstable (happens with GFN2 and solvation)." << std::endl;
    }

    if (Json2KeyWord<std::string>(m_defaults, "wall").compare("spheric") == 0) {
        if (Json2KeyWord<std::string>(m_defaults, "wall_type").compare("logfermi") == 0) {
            m_wall_type = 1;
            WallPotential = [=](double* grad) -> double {
                this->m_wall_potential = this->ApplySphericLogFermiWalls(grad);
                return m_wall_potential;
            };
        } else if (Json2KeyWord<std::string>(m_defaults, "wall_type").compare("harmonic") == 0) {
            m_wall_type = 1;
            WallPotential = [=](double* grad) -> double {
                this->m_wall_potential = this->ApplySphericHarmonicWalls(grad);
                return m_wall_potential;
            };
        } else {
            std::cout << "Did not understand wall potential input. Exit now!" << std::endl;
            exit(1);
        }
        std::cout << "Setting up spherical potential" << std::endl;

    } else if (Json2KeyWord<std::string>(m_defaults, "wall").compare("rect") == 0) {
        if (Json2KeyWord<std::string>(m_defaults, "wall_type").compare("logfermi") == 0) {
            m_wall_type = 2;
            WallPotential = [=](double* grad) -> double {
                this->m_wall_potential = this->ApplyRectLogFermiWalls(grad);
                return m_wall_potential;
            };
        } else if (Json2KeyWord<std::string>(m_defaults, "wall_type").compare("harmonic") == 0) {
            m_wall_type = 2;
            WallPotential = [=](double* grad) -> double {
                this->m_wall_potential = this->ApplyRectHarmonicWalls(grad);
                return m_wall_potential;
            };

        } else {
            std::cout << "Did not understand wall potential input. Exit now!" << std::endl;
            exit(1);
        }
        std::cout << "Setting up rectangular potential" << std::endl;
    } else
        WallPotential = [=](double* grad) -> double {
            return 0;
        };
    m_rm_COM_step = m_rm_COM / m_dT;
}

bool SimpleMD::Initialise()
{
    static std::random_device rd{};
    static std::mt19937 gen{ rd() };
    if (m_seed == -1) {
        const auto start = std::chrono::high_resolution_clock::now();
        m_seed = std::chrono::duration_cast<std::chrono::seconds>(start.time_since_epoch()).count();
    } else if (m_seed == 0)
        m_seed = m_T0 * m_mass.size();
    std::cout << "Random seed is " << m_seed << std::endl;
    gen.seed(m_seed);

    if (m_initfile.compare("none") != 0) {
        json md;
        std::ifstream restart_file(m_initfile);
        try {
            restart_file >> md;
        } catch (nlohmann::json::type_error& e) {
            throw 404;
        } catch (nlohmann::json::parse_error& e) {
            throw 404;
        }
        LoadRestartInformation(md);
        m_restart = true;
    } else if (!m_norestart)
        LoadRestartInformation();

    if (m_molecule.AtomCount() == 0)
        return false;

    if (!m_restart) {
        std::ofstream result_file;
        result_file.open(Basename() + ".trj.xyz");
        result_file.close();
    }
    m_natoms = m_molecule.AtomCount();
    m_molecule.setCharge(0);
    if (!m_nocenter) {
        std::cout << "Move stucture to the origin ... " << std::endl;
        m_molecule.Center(m_COM);
    } else
        std::cout << "Move stucture NOT to the origin ... " << std::endl;

    m_mass = std::vector<double>(3 * m_natoms, 0);
    m_rmass = std::vector<double>(3 * m_natoms, 0);
    m_atomtype = std::vector<int>(m_natoms, 0);

    if (!m_restart) {
        m_current_geometry = std::vector<double>(3 * m_natoms, 0);
        m_velocities = std::vector<double>(3 * m_natoms, 0);
        m_currentStep = 0;
    }

    m_gradient = std::vector<double>(3 * m_natoms, 0);
    m_virial = std::vector<double>(3 * m_natoms, 0);
    m_atom_temp = std::vector<std::vector<double>>(m_natoms);
    if(m_opt)
    {
        json js = CurcumaOptJson;
        js = MergeJson(js, m_defaults);
        js["writeXYZ"] = false;
        js["method"] = m_method;
        /*
        try {
            js["threads"] = m_defaults["threads"].get<int>();
        }
        catch (const nlohmann::detail::type_error& error) {

           }*/
        CurcumaOpt optimise(js, true);
        optimise.addMolecule(&m_molecule);
        optimise.start();
        auto mol = optimise.Molecules();

        auto molecule = ((*mol)[0]);
        m_molecule.setGeometry(molecule.getGeometry());
        m_molecule.appendXYZFile(Basename() + ".opt.xyz");
    }
    double mass = 0;
    for (int i = 0; i < m_natoms; ++i) {
        m_atomtype[i] = m_molecule.Atom(i).first;
        if (!m_restart) {
            Position pos = m_molecule.Atom(i).second;
            m_current_geometry[3 * i + 0] = pos(0) / 1;
            m_current_geometry[3 * i + 1] = pos(1) / 1;
            m_current_geometry[3 * i + 2] = pos(2) / 1;
        }
        if (m_atomtype[i] == 1) {
            m_mass[3 * i + 0] = Elements::AtomicMass[m_atomtype[i]] * m_hmass;
            m_mass[3 * i + 1] = Elements::AtomicMass[m_atomtype[i]] * m_hmass;
            m_mass[3 * i + 2] = Elements::AtomicMass[m_atomtype[i]] * m_hmass;
            mass += Elements::AtomicMass[m_atomtype[i]] * m_hmass;

            m_rmass[3 * i + 0] = 1 / m_mass[3 * i + 0];
            m_rmass[3 * i + 1] = 1 / m_mass[3 * i + 1];
            m_rmass[3 * i + 2] = 1 / m_mass[3 * i + 2];
        } else {
            m_mass[3 * i + 0] = Elements::AtomicMass[m_atomtype[i]];
            m_mass[3 * i + 1] = Elements::AtomicMass[m_atomtype[i]];
            m_mass[3 * i + 2] = Elements::AtomicMass[m_atomtype[i]];
            mass += Elements::AtomicMass[m_atomtype[i]];

            m_rmass[3 * i + 0] = 1 / m_mass[3 * i + 0];
            m_rmass[3 * i + 1] = 1 / m_mass[3 * i + 1];
            m_rmass[3 * i + 2] = 1 / m_mass[3 * i + 2];
        }
    }

    m_molecule.setCharge(m_charge);
    m_molecule.setSpin(m_spin);
    m_interface->setMolecule(m_molecule);

    if (m_writeUnique) {
        json rmsdtraj = RMSDTrajJson;
        rmsdtraj["writeUnique"] = true;
        rmsdtraj["rmsd"] = m_rmsd;
        rmsdtraj["writeRMSD"] = false;
        m_unqiue = new RMSDTraj(rmsdtraj, true);
        m_unqiue->setBaseName(Basename() + ".xyz");
        m_unqiue->Initialise();
    }
    m_dof = 3 * m_natoms;

    InitConstrainedBonds();
    InitialiseWalls();
    if (!m_restart) {
        InitVelocities(m_scale_velo);
        m_xi.resize(m_chain_length, 0.0);
        m_Q.resize(m_chain_length, 100); // Setze eine geeignete Masse für jede Kette
        for (int i = 0; i < m_chain_length; ++i) {
            m_xi[i] = pow(10.0, double(i)) - 1;
            m_Q[i] = pow(10, i) * kb_Eh * m_T0 * m_dof * 100;
            std::cout << m_xi[i] << "  " << m_Q[i] << std::endl;
        }
        m_eta = 0.0;
    }
    if (m_writeinit) {
        json init = WriteRestartInformation();
        std::ofstream result_file;
        result_file.open(Basename() + ".init.json");
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
            std::cout << i << " ";
            m_rmsd_mtd_molecule.addPair(m_molecule.Atom(i));
        }
        m_rmsd_fragment_count = m_rmsd_mtd_molecule.GetFragments().size();

        json config = RMSDJson;
        config["silent"] = true;
        config["reorder"] = false;
        for (int i = 0; i < m_threads; ++i) {
            BiasThread* thread = new BiasThread(m_rmsd_mtd_molecule, config, m_nocolvarfile, m_nohillsfile);
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
            if (m_rmsd_ref_file.compare("none") != 0) {
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

    if (m_rattle) {
        auto m = m_molecule.DistanceMatrix();
        m_topo_initial = m.second;
        for (int i = 0; i < m_molecule.AtomCount(); ++i) {
            for (int j = 0; j < i; ++j) {
                if (m.second(i, j)) {
                    if (m_rattle == 2) {
                        if (m_molecule.Atom(i).first != 1 && m_molecule.Atom(j).first != 1)
                            continue;
                    }
                    std::pair<int, int> indicies(i, j);
                    std::pair<double, double> minmax(m_molecule.CalculateDistance(i, j) * m_molecule.CalculateDistance(i, j), m_molecule.CalculateDistance(i, j) * m_molecule.CalculateDistance(i, j));
                    std::pair<std::pair<int, int>, double> bond(indicies, m_molecule.CalculateDistance(i, j) * m_molecule.CalculateDistance(i, j));
                    m_bond_constrained.push_back(std::pair<std::pair<int, int>, double>(bond));

                    // std::cout << i << " " << j << " " << bond.second << std::endl;
                }
            }
        }
    }

    std::cout << m_dof << " initial degrees of freedom " << std::endl;
    std::cout << m_bond_constrained.size() << " constrains active" << std::endl;
    m_dof -= m_bond_constrained.size();
    std::cout << m_dof << " degrees of freedom remaining ..." << std::endl;
}

void SimpleMD::InitVelocities(double scaling)
{
    static std::default_random_engine generator;
    for (size_t i = 0; i < m_natoms; ++i) {
        std::normal_distribution<double> distribution(0.0, std::sqrt(kb_Eh * m_T0 * m_rmass[i]));
        m_velocities[3 * i + 0] = distribution(generator);
        m_velocities[3 * i + 1] = distribution(generator);
        m_velocities[3 * i + 2] = distribution(generator);
    }
    RemoveRotation(m_velocities);
    EKin();
    double coupling = m_coupling;
    m_coupling = m_dT;
    Berendson();
    Berendson();
    EKin();
    m_coupling = coupling;
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
    m_wall_spheric_radius = Json2KeyWord<double>(m_defaults, "wall_spheric_radius");
    m_wall_temp = Json2KeyWord<double>(m_defaults, "wall_temp");
    m_wall_beta = Json2KeyWord<double>(m_defaults, "wall_beta");

    m_wall_x_min = Json2KeyWord<double>(m_defaults, "wall_x_min");
    m_wall_x_max = Json2KeyWord<double>(m_defaults, "wall_x_max");
    m_wall_y_min = Json2KeyWord<double>(m_defaults, "wall_y_min");
    m_wall_y_max = Json2KeyWord<double>(m_defaults, "wall_y_max");
    m_wall_z_min = Json2KeyWord<double>(m_defaults, "wall_z_min");
    m_wall_z_max = Json2KeyWord<double>(m_defaults, "wall_z_max");
    std::vector<double> box = m_molecule.GetBox();
    double radius = 0;
    if (m_wall_x_min - m_wall_x_max < 1) {
        m_wall_x_min = -box[0] * 0.75;
        m_wall_x_max = -1 * m_wall_x_min;
        radius = std::max(radius, box[0]);
    }

    if (m_wall_y_min - m_wall_y_max < 1) {
        m_wall_y_min = -box[1] * 0.75;
        m_wall_y_max = -1 * m_wall_y_min;
        radius = std::max(radius, box[1]);
    }

    if (m_wall_z_min - m_wall_z_max < 1) {
        m_wall_z_min = -box[2] * 0.75;
        m_wall_z_max = -1 * m_wall_z_min;
        radius = std::max(radius, box[2]);
    }

    if (m_wall_spheric_radius < radius) {
        m_wall_spheric_radius = radius + 5;
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
    std::cout << "Setting up potential walls " << std::endl;
    std::cout << "Radius " << m_wall_spheric_radius << std::endl;
    std::cout << "x-range " << m_wall_x_min << " ... " << m_wall_x_max << std::endl;
    std::cout << "y-range " << m_wall_y_min << " ... " << m_wall_y_max << std::endl;
    std::cout << "z-range " << m_wall_z_min << " ... " << m_wall_z_max << std::endl;
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
    restart["velocities"] = Tools::DoubleVector2String(m_velocities);
    restart["geometry"] = Tools::DoubleVector2String(m_current_geometry);
    restart["gradient"] = Tools::DoubleVector2String(m_gradient);
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
    restart["rattle_dynamic_tol"] = m_rattle_tolerance;
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
        for (int i = 0; i < m_bias_threads.size(); ++i) {
            for (const auto& stored_bias : m_bias_threads[i]->getBias()) {
                bias[stored_bias["index"]] = stored_bias;
            }
        }
        json bias_restart;
        for (int i = 0; i < bias.size(); ++i) {
            bias_restart[i] = bias[i];
        }
        restart["bias"] = bias_restart;
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
        } catch (json::type_error& e) {
            error++;
            continue;
        } catch (json::parse_error& e) {
            error++;
            continue;
        }

        json md;
        try {
            md = restart[MethodName()[0]];
        } catch (json::type_error& e) {
            error++;
            continue;
        }
        return LoadRestartInformation(md);
    }
    return true;
};

bool SimpleMD::LoadRestartInformation(const json& state)
{
    std::string geometry, velocities, constrains, xi, Q;

    try {
        m_method = state["method"];
    } catch (json::type_error& e) {
    }
    try {
        m_dT = state["dT"];
    } catch (json::type_error& e) {
    }
    try {
        m_maxtime = state["MaxTime"];
    } catch (json::type_error& e) {
    }
    try {
        m_rmrottrans = state["rmrottrans"];
    } catch (json::type_error& e) {
    }
    try {
        m_nocenter = state["nocenter"];
    } catch (json::type_error& e) {
    }
    try {
        m_COM = state["COM"];
    } catch (json::type_error& e) {
    }
    try {
        m_T0 = state["T"];
    } catch (json::type_error& e) {
    }
    try {
        m_currentStep = state["currentStep"];
    } catch (json::type_error& e) {
    }

    try {
        m_aver_Epot = state["average_Epot"];
    } catch (json::type_error& e) {
    }

    try {
        m_aver_Ekin = state["average_Ekin"];
    } catch (json::type_error& e) {
    }

    try {
        m_aver_Etot = state["average_Etot"];
    } catch (json::type_error& e) {
    }
    try {
        m_aver_Temp = state["average_T"];
    } catch (json::type_error& e) {
    }

    try {
        m_average_virial_correction = state["average_Virial"];
    } catch (json::type_error& e) {
    }

    try {
        m_average_wall_potential = state["average_Wall"];
    } catch (json::type_error& e) {
    }

    try {
        m_coupling = state["coupling"];
    } catch (json::type_error& e) {
    }

    try {
        m_respa = state["respa"];
    } catch (json::type_error& e) {
    }

    try {
        m_eta = state["eta"];
    } catch (json::type_error& e) {
    }

    try {
        m_thermostat = state["thermostat"];
    } catch (json::type_error& e) {
    }

    try {
        geometry = state["geometry"];
    } catch (json::type_error& e) {
    }

    try {
        velocities = state["velocities"];
    } catch (json::type_error& e) {
    }

    try {
        xi = state["xi"];
    } catch (json::type_error& e) {
    }
    try {
        Q = state["Q"];
    } catch (json::type_error& e) {
    }

    try {
        m_mtd = state["mtd"];
    } catch (json::type_error& e) {
    }

    try {
        m_rattle = state["rattle"];
    } catch (json::type_error& e) {
    }

    try {
        m_rattle_tolerance = state["rattle_tolerance"];
    } catch (json::type_error& e) {
    }

    try {
        m_rattle_maxiter = state["rattle_maxiter"];
    } catch (json::type_error& e) {
    }

    try {
        m_rattle_dynamic_tol = state["rattle_dynamic_tol"];
    } catch (json::type_error& e) {
    }

    try {
        m_rattle_dynamic_tol_iter = state["rattle_dynamic_tol_iter"];
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
    } catch (json::type_error& e) {
    }

    if (geometry.size()) {
        m_current_geometry = Tools::String2DoubleVec(geometry, "|");
    }
    if (velocities.size()) {
        m_velocities = Tools::String2DoubleVec(velocities, "|");
    }

    if (xi.size()) {
        m_xi = Tools::String2DoubleVec(xi, "|");
    }

    if (Q.size()) {
        m_Q = Tools::String2DoubleVec(Q, "|");
    }

    m_restart = geometry.size() && velocities.size();

    return true;
}

void SimpleMD::start()
{
    if (m_initialised == false)
        return;
    bool aborted = false;
    auto unix_timestamp = std::chrono::seconds(std::time(NULL));
    m_unix_started = std::chrono::milliseconds(unix_timestamp).count();
    double* gradient = new double[3 * m_natoms];
    std::vector<json> states;
    for (int i = 0; i < 3 * m_natoms; ++i) {
        gradient[i] = 0;
    }

    if (m_thermostat.compare("csvr") == 0) {
        fmt::print(fg(fmt::color::green) | fmt::emphasis::bold, "\nUsing Canonical sampling through velocity rescaling (CSVR) Thermostat\nJ. Chem. Phys. 126, 014101 (2007) - DOI: 10.1063/1.2408420\n\n");
        ThermostatFunction = std::bind(&SimpleMD::CSVR, this);
    } else if (m_thermostat.compare("berendson") == 0) {
        fmt::print(fg(fmt::color::green) | fmt::emphasis::bold, "\nUsing Berendson Thermostat\nJ. Chem. Phys. 81, 3684 (1984) - DOI: 10.1063/1.448118\n\n");
        ThermostatFunction = std::bind(&SimpleMD::Berendson, this);
    } else if (m_thermostat.compare("anderson") == 0) {
        fmt::print(fg(fmt::color::green) | fmt::emphasis::bold, "\nUsing Anderson Thermostat\n ... \n\n");
        ThermostatFunction = std::bind(&SimpleMD::Anderson, this);
    } else if (m_thermostat.compare("nosehover") == 0) {
        fmt::print(fg(fmt::color::green) | fmt::emphasis::bold, "\nUsing Nosé-Hoover-Chain Thermostat\n ... \n\n");
        ThermostatFunction = std::bind(&SimpleMD::NoseHover, this);
    } else {
        ThermostatFunction = std::bind(&SimpleMD::None, this);
        std::cout << "No Thermostat applied\n"
                  << std::endl;
    }

    m_Epot = Energy(gradient);
    EKin();
    m_Etot = m_Epot + m_Ekin;
    AverageQuantities();
    int m_step = 0;

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
        plumed_cmd(m_plumedmain, "setPositions", &m_current_geometry[0]);
        plumed_cmd(m_plumedmain, "setEnergy", &m_Epot);
        plumed_cmd(m_plumedmain, "setForces", &m_gradient[0]);
        plumed_cmd(m_plumedmain, "setVirial", &m_virial[0]);
        plumed_cmd(m_plumedmain, "setMasses", &m_mass[0]);
        plumed_cmd(m_plumedmain, "prepareCalc", NULL);
        plumed_cmd(m_plumedmain, "performCalc", NULL);
    }
#endif
    std::vector<double> charge(0, m_natoms);

#ifdef GCC
    //         std::cout << fmt::format("{0: ^{0}} {1: ^{1}} {2: ^{2}} {3: ^{3}} {4: ^{4}}\n", "Step", "Epot", "Ekin", "Etot", "T");
    // std::cout << fmt::format("{1: ^{0}} {1: ^{1}} {1: ^{2}} {1: ^{3}} {1: ^{4}}\n", "", "Eh", "Eh", "Eh", "K");
#else
    std::cout << "Step"
              << "\t"
              << "Epot"
              << "\t"
              << "Ekin"
              << "\t"
              << "Etot"
              << "\t"
              << "T" << std::endl;
    std::cout << "  "
              << "\t"
              << "Eh"
              << "\t"
              << "Eh"
              << "\t"
              << "Eh"
              << "\t"
              << "T" << std::endl;
#endif
    if (m_rmsd_mtd) {
        std::cout << "k\t" << m_k_rmsd << std::endl;
        std::cout << "alpha\t" << m_alpha_rmsd << std::endl;
        std::cout << "steps\t" << m_mtd_steps << std::endl;
        std::cout << "Ethresh\t" << m_rmsd_econv << std::endl;
        if (m_wtmtd)
            std::cout << "Well Tempered\tOn (" << m_rmsd_DT << ")" << std::endl;
        else
            std::cout << "Well Tempered\tOff" << std::endl;
    }
    PrintStatus();

    /* Start MD Lopp here */
    for (; m_currentStep < m_maxtime;) {
        auto step0 = std::chrono::system_clock::now();

        if (CheckStop() == true) {
            TriggerWriteRestart();
            aborted = true;
#ifdef USE_Plumed
            if (m_mtd) {
                plumed_finalize(m_plumedmain); // Call the plumed destructor
            }
#endif

            break;
        }

        if (m_rm_COM_step > 0 && m_step % m_rm_COM_step == 0) {
            // std::cout << "Removing COM motion." << std::endl;
            if (m_rmrottrans == 1)
                RemoveRotation(m_velocities);
            else if (m_rmrottrans == 2)
                RemoveRotations(m_velocities);
            else if (m_rmrottrans == 3) {
                RemoveRotations(m_velocities);
                RemoveRotation(m_velocities);
            }
        }

        Integrator(gradient);
        AverageQuantities();

        if (m_mtd) {
            if (!m_eval_mtd) {
                if (std::abs(m_T0 - m_aver_Temp) < m_mtd_dT && m_step > 10) {
                    m_eval_mtd = true;
                    std::cout << "Starting with MetaDynamics ..." << std::endl;
                }
            }
        }

/////////// Dipole calc




/////////// Dipole calc
        if (m_step % m_dump == 0) {
            bool write = WriteGeometry();
            if (write) {
                states.push_back(WriteRestartInformation());
                m_current_rescue = 0;
            } else if (!write && m_rescue && states.size() > (1 - m_current_rescue)) {
                std::cout << "Molecule exploded, resetting to previous state ..." << std::endl;
                LoadRestartInformation(states[states.size() - 1 - m_current_rescue]);
                Geometry geometry = m_molecule.getGeometry();
                for (int i = 0; i < m_natoms; ++i) {
                    geometry(i, 0) = m_current_geometry[3 * i + 0] * au;
                    geometry(i, 1) = m_current_geometry[3 * i + 1] * au;
                    geometry(i, 2) = m_current_geometry[3 * i + 2] * au;
                }
                m_molecule.setGeometry(geometry);
                m_molecule.GetFragments();
                InitVelocities(-1);
                Energy(gradient);
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

            std::ofstream restart_file("unstable_curcuma.json");
            nlohmann::json restart;
            restart[MethodName()[0]] = WriteRestartInformation();
            restart_file << restart << std::endl;

            m_time_step = 0;
            aborted = true;

#ifdef USE_Plumed
            if (m_mtd) {
                plumed_finalize(m_plumedmain); // Call the plumed destructor
            }
#endif
            return;
        }

        if (m_writerestart > -1 && m_step % m_writerestart == 0) {
            std::ofstream restart_file("curcuma_step_" + std::to_string(int(m_step * m_dT)) + ".json");
            nlohmann::json restart;
            restart[MethodName()[0]] = WriteRestartInformation();
            restart_file << restart << std::endl;
        }
        if ((m_step && int(m_step * m_dT) % m_print == 0)) {
            m_Etot = m_Epot + m_Ekin;
            PrintStatus();
            m_time_step = 0;
        }
        if (m_rattle && m_rattle_dynamic_tol) {
            m_aver_rattle_Temp += m_T;
            m_rattle_counter++;
            if (m_rattle_counter == m_rattle_dynamic_tol_iter)
                AdjustRattleTolerance();
        }
        if (m_impuls > m_T) {
            InitVelocities(m_scale_velo * m_impuls_scaling);
            EKin();
            // PrintStatus();
            m_time_step = 0;
        }

        if (m_current_rescue >= m_max_rescue) {
            fmt::print(fg(fmt::color::salmon) | fmt::emphasis::bold, "Nothing really helps");
            break;
        }
        m_step++;
        m_currentStep += m_dT;
        m_time_step += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - step0).count();
    }
    PrintStatus();
    if (m_thermostat.compare("csvr") == 0)
        std::cout << "Exchange with heat bath " << m_Ekin_exchange << "Eh" << std::endl;
    if (m_dipole) {
        /*
        double dipole = 0.0;
        for( auto d : m_collected_dipole)
            dipole += d;
        dipole /= m_collected_dipole.size();
        std::cout << dipole*2.5418 << " average dipole in Debye and " << dipole*2.5418*3.3356e-30 << " Cm" << std::endl;
        */
        std::cout << "Calculated averaged dipole moment " << m_aver_dipol * 2.5418 << " Debye and " << m_aver_dipol * 2.5418 * 3.3356 << " Cm [e-30]" << std::endl;
    }

#ifdef USE_Plumed
    if (m_mtd) {
        plumed_finalize(m_plumedmain); // Call the plumed destructor
    }
#endif
    if (m_rmsd_mtd) {
        std::cout << "Sum of Energy of COLVARs:" << std::endl;
        for (int i = 0; i < m_bias_threads.size(); ++i) {
            auto structures = m_bias_threads[i]->getBiasStructure();
            for (int j = 0; j < structures.size(); ++j) {
                std::cout << structures[j].rmsd_reference << "\t" << structures[j].energy << "\t" << structures[j].counter / double(m_colvar_incr) * 100 << std::endl;

                m_rmsd_mtd_molecule.setGeometry(structures[j].geometry);
                m_rmsd_mtd_molecule.setEnergy(structures[j].energy);
                m_rmsd_mtd_molecule.setName(std::to_string(structures[j].index) + " " + std::to_string(structures[j].rmsd_reference));
                if (i == j && i == 0)
                    m_rmsd_mtd_molecule.writeXYZFile(Basename() + ".mtd.xyz");
                else
                    m_rmsd_mtd_molecule.appendXYZFile(Basename() + ".mtd.xyz");
            }
        }
    }
    std::ofstream restart_file("curcuma_final.json");
    nlohmann::json restart;
    restart[MethodName()[0]] = WriteRestartInformation();
    restart_file << restart << std::endl;
    if (aborted == false)
        std::remove("curcuma_restart.json");
    delete[] gradient;
}

void SimpleMD::AdjustRattleTolerance()
{
    m_aver_rattle_Temp /= double(m_rattle_counter);

    // std::pair<double, double> pair(m_rattle_tolerance, m_aver_Temp);

    if (m_aver_rattle_Temp > m_T0)
        m_rattle_tolerance -= 0.01;
    else if (m_aver_rattle_Temp < m_T0)
        m_rattle_tolerance += 0.01;
    std::cout << m_rattle_counter << " " << m_aver_rattle_Temp << " " << m_rattle_tolerance << std::endl;
    m_rattle_tolerance = std::abs(m_rattle_tolerance);
    m_rattle_counter = 0;
    m_aver_rattle_Temp = 0;
}

void SimpleMD::Verlet(double* grad)
{
    double ekin = 0;

    for (int i = 0; i < m_natoms; ++i) {
        m_current_geometry[3 * i + 0] = m_current_geometry[3 * i + 0] + m_dT * m_velocities[3 * i + 0] - 0.5 * grad[3 * i + 0] * m_rmass[3 * i + 0] * m_dt2;
        m_current_geometry[3 * i + 1] = m_current_geometry[3 * i + 1] + m_dT * m_velocities[3 * i + 1] - 0.5 * grad[3 * i + 1] * m_rmass[3 * i + 1] * m_dt2;
        m_current_geometry[3 * i + 2] = m_current_geometry[3 * i + 2] + m_dT * m_velocities[3 * i + 2] - 0.5 * grad[3 * i + 2] * m_rmass[3 * i + 2] * m_dt2;

        m_velocities[3 * i + 0] = m_velocities[3 * i + 0] - 0.5 * m_dT * grad[3 * i + 0] * m_rmass[3 * i + 0];
        m_velocities[3 * i + 1] = m_velocities[3 * i + 1] - 0.5 * m_dT * grad[3 * i + 1] * m_rmass[3 * i + 1];
        m_velocities[3 * i + 2] = m_velocities[3 * i + 2] - 0.5 * m_dT * grad[3 * i + 2] * m_rmass[3 * i + 2];
        ekin += m_mass[i] * (m_velocities[3 * i] * m_velocities[3 * i] + m_velocities[3 * i + 1] * m_velocities[3 * i + 1] + m_velocities[3 * i + 2] * m_velocities[3 * i + 2]);
    }
    ekin *= 0.5;
    m_T = 2.0 * ekin / (kb_Eh * m_dof);
    m_Ekin = ekin;
    ThermostatFunction();
    m_Epot = Energy(grad);
    if (m_rmsd_mtd) {
        if (m_step % m_mtd_steps == 0) {
            ApplyRMSDMTD(grad);
        }
    }
#ifdef USE_Plumed
    if (m_mtd) {
        plumed_cmd(m_plumedmain, "setStep", &m_step);

        plumed_cmd(m_plumedmain, "setPositions", &m_current_geometry[0]);

        plumed_cmd(m_plumedmain, "setEnergy", &m_Epot);
        plumed_cmd(m_plumedmain, "setForces", &m_gradient[0]);
        plumed_cmd(m_plumedmain, "setVirial", &m_virial[0]);

        plumed_cmd(m_plumedmain, "setMasses", &m_mass[0]);
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
    WallPotential(grad);
    ekin = 0.0;

    for (int i = 0; i < m_natoms; ++i) {
        m_velocities[3 * i + 0] -= 0.5 * m_dT * grad[3 * i + 0] * m_rmass[3 * i + 0];
        m_velocities[3 * i + 1] -= 0.5 * m_dT * grad[3 * i + 1] * m_rmass[3 * i + 1];
        m_velocities[3 * i + 2] -= 0.5 * m_dT * grad[3 * i + 2] * m_rmass[3 * i + 2];

        ekin += m_mass[i] * (m_velocities[3 * i] * m_velocities[3 * i] + m_velocities[3 * i + 1] * m_velocities[3 * i + 1] + m_velocities[3 * i + 2] * m_velocities[3 * i + 2]);
        m_gradient[3 * i + 0] = grad[3 * i + 0];
        m_gradient[3 * i + 1] = grad[3 * i + 1];
        m_gradient[3 * i + 2] = grad[3 * i + 2];
    }
    ekin *= 0.5;
    double T = 2.0 * ekin / (kb_Eh * m_dof);
    m_unstable = T > 10000 * m_T || std::isnan(T);
    m_T = T;
    m_Ekin = ekin;
    ThermostatFunction();
    EKin();
}

void SimpleMD::Rattle(double* grad)
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
     */
    TriggerWriteRestart();
    double* coord = new double[3 * m_natoms];
    double m_dT_inverse = 1 / m_dT;
    std::vector<int> moved(m_natoms, 0);
    bool move = false;
    for (int i = 0; i < m_natoms; ++i) {
        coord[3 * i + 0] = m_current_geometry[3 * i + 0] + m_dT * m_velocities[3 * i + 0] - 0.5 * grad[3 * i + 0] * m_rmass[3 * i + 0] * m_dt2;
        coord[3 * i + 1] = m_current_geometry[3 * i + 1] + m_dT * m_velocities[3 * i + 1] - 0.5 * grad[3 * i + 1] * m_rmass[3 * i + 1] * m_dt2;
        coord[3 * i + 2] = m_current_geometry[3 * i + 2] + m_dT * m_velocities[3 * i + 2] - 0.5 * grad[3 * i + 2] * m_rmass[3 * i + 2] * m_dt2;

        m_velocities[3 * i + 0] -= 0.5 * m_dT * grad[3 * i + 0] * m_rmass[3 * i + 0];
        m_velocities[3 * i + 1] -= 0.5 * m_dT * grad[3 * i + 1] * m_rmass[3 * i + 1];
        m_velocities[3 * i + 2] -= 0.5 * m_dT * grad[3 * i + 2] * m_rmass[3 * i + 2];
    }
    double iter = 0;
    double max_mu = 10;

    while (iter < m_rattle_maxiter) {
        iter++;
        int active = 0;
        for (auto bond : m_bond_constrained) {
            int i = bond.first.first, j = bond.first.second;
            double distance = bond.second;
            double distance_current = ((coord[3 * i + 0] - coord[3 * j + 0]) * (coord[3 * i + 0] - coord[3 * j + 0])
                + (coord[3 * i + 1] - coord[3 * j + 1]) * (coord[3 * i + 1] - coord[3 * j + 1])
                + (coord[3 * i + 2] - coord[3 * j + 2]) * (coord[3 * i + 2] - coord[3 * j + 2]));

            if (std::abs(distance - distance_current) > m_rattle_tolerance) {
                move = true;
                double r = distance - distance_current;
                double dx = m_current_geometry[3 * i + 0] - m_current_geometry[3 * j + 0];
                double dy = m_current_geometry[3 * i + 1] - m_current_geometry[3 * j + 1];
                double dz = m_current_geometry[3 * i + 2] - m_current_geometry[3 * j + 2];

                double scalarproduct = (dx) * (coord[3 * i + 0] - coord[3 * j + 0])
                    + (dy) * (coord[3 * i + 1] - coord[3 * j + 1])
                    + (dz) * (coord[3 * i + 2] - coord[3 * j + 2]);
                if (scalarproduct >= m_rattle_tolerance * distance) {
                    moved[i] = 1;
                    moved[j] = 1;
                    active++;

                    double lambda = r / (1 * (m_rmass[i] + m_rmass[j]) * scalarproduct);
                    while (std::abs(lambda) > max_mu)
                        lambda /= 2;
                    coord[3 * i + 0] += dx * lambda * 0.5 * m_rmass[i];
                    coord[3 * i + 1] += dy * lambda * 0.5 * m_rmass[i];
                    coord[3 * i + 2] += dz * lambda * 0.5 * m_rmass[i];

                    coord[3 * j + 0] -= dx * lambda * 0.5 * m_rmass[j];
                    coord[3 * j + 1] -= dy * lambda * 0.5 * m_rmass[j];
                    coord[3 * j + 2] -= dz * lambda * 0.5 * m_rmass[j];

                    m_velocities[3 * i + 0] += dx * lambda * 0.5 * m_rmass[i] * m_dT_inverse;
                    m_velocities[3 * i + 1] += dy * lambda * 0.5 * m_rmass[i] * m_dT_inverse;
                    m_velocities[3 * i + 2] += dz * lambda * 0.5 * m_rmass[i] * m_dT_inverse;

                    m_velocities[3 * j + 0] -= dx * lambda * 0.5 * m_rmass[j] * m_dT_inverse;
                    m_velocities[3 * j + 1] -= dy * lambda * 0.5 * m_rmass[j] * m_dT_inverse;
                    m_velocities[3 * j + 2] -= dz * lambda * 0.5 * m_rmass[j] * m_dT_inverse;
                }
            }
        }
        if (active == 0)
            break;
    }
    if (iter >= m_rattle_maxiter) {
        std::cout << "numeric difficulties - 1st step in rattle velocity verlet" << std::endl;
        std::ofstream restart_file("unstable_curcuma_" + std::to_string(m_currentStep) + ".json");
        nlohmann::json restart;
        restart[MethodName()[0]] = WriteRestartInformation();
        restart_file << restart << std::endl;
    }
    double ekin = 0;

    for (int i = 0; i < m_natoms; ++i) {
        m_current_geometry[3 * i + 0] = coord[3 * i + 0];
        m_current_geometry[3 * i + 1] = coord[3 * i + 1];
        m_current_geometry[3 * i + 2] = coord[3 * i + 2];
        ekin += m_mass[i] * (m_velocities[3 * i] * m_velocities[3 * i] + m_velocities[3 * i + 1] * m_velocities[3 * i + 1] + m_velocities[3 * i + 2] * m_velocities[3 * i + 2]);
    }
    ekin *= 0.5;
    m_T = 2.0 * ekin / (kb_Eh * m_dof);
    m_Ekin = ekin;
    ThermostatFunction();
    m_Epot = Energy(grad);

    if (m_rmsd_mtd) {
        if (m_step % m_mtd_steps == 0) {
            ApplyRMSDMTD(grad);
        }
    }
#ifdef USE_Plumed
    if (m_mtd) {
        plumed_cmd(m_plumedmain, "setStep", &m_step);

        plumed_cmd(m_plumedmain, "setPositions", &m_current_geometry[0]);

        plumed_cmd(m_plumedmain, "setEnergy", &m_Epot);
        plumed_cmd(m_plumedmain, "setForces", &m_gradient[0]);
        plumed_cmd(m_plumedmain, "setVirial", &m_virial[0]);

        plumed_cmd(m_plumedmain, "setMasses", &m_mass[0]);
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
    WallPotential(grad);

    for (int i = 0; i < m_natoms; ++i) {
        m_velocities[3 * i + 0] -= 0.5 * m_dT * grad[3 * i + 0] * m_rmass[3 * i + 0];
        m_velocities[3 * i + 1] -= 0.5 * m_dT * grad[3 * i + 1] * m_rmass[3 * i + 1];
        m_velocities[3 * i + 2] -= 0.5 * m_dT * grad[3 * i + 2] * m_rmass[3 * i + 2];

        m_gradient[3 * i + 0] = grad[3 * i + 0];
        m_gradient[3 * i + 1] = grad[3 * i + 1];
        m_gradient[3 * i + 2] = grad[3 * i + 2];
    }
    m_virial_correction = 0;
    iter = 0;
    ekin = 0.0;
    while (iter < m_rattle_maxiter) {
        iter++;
        int active = 0;
        for (auto bond : m_bond_constrained) {
            int i = bond.first.first, j = bond.first.second;
            if (moved[i] == 1 && moved[j] == 1) {
                double distance = bond.second;

                double dx = coord[3 * i + 0] - coord[3 * j + 0];
                double dy = coord[3 * i + 1] - coord[3 * j + 1];
                double dz = coord[3 * i + 2] - coord[3 * j + 2];
                double dvx = m_velocities[3 * i + 0] - m_velocities[3 * j + 0];
                double dvy = m_velocities[3 * i + 1] - m_velocities[3 * j + 1];
                double dvz = m_velocities[3 * i + 2] - m_velocities[3 * j + 2];

                double r = (dx) * (dvx) + (dy) * (dvy) + (dz) * (dvz);

                double mu = -1 * r / ((m_rmass[i] + m_rmass[j]) * distance);
                while (std::abs(mu) > max_mu)
                    mu /= 2;
                if (std::abs(mu) > m_rattle_tolerance) {
                    active = 1;
                    m_virial_correction += mu * distance;
                    m_velocities[3 * i + 0] += dx * mu * m_rmass[i];
                    m_velocities[3 * i + 1] += dy * mu * m_rmass[i];
                    m_velocities[3 * i + 2] += dz * mu * m_rmass[i];

                    m_velocities[3 * j + 0] -= dx * mu * m_rmass[j];
                    m_velocities[3 * j + 1] -= dy * mu * m_rmass[j];
                    m_velocities[3 * j + 2] -= dz * mu * m_rmass[j];
                }
            }
        }
        if (active == 0)
            break;
    }

    if (iter >= m_rattle_maxiter) {
        std::cout << "numeric difficulties - 2nd in rattle velocity verlet" << iter << std::endl;
        std::ofstream restart_file("unstable_curcuma_" + std::to_string(m_currentStep) + ".json");
        nlohmann::json restart;
        restart[MethodName()[0]] = WriteRestartInformation();
        restart_file << restart << std::endl;
    }

    if (move)
        RemoveRotations(m_velocities);

    delete[] coord;
    for (int i = 0; i < m_natoms; ++i) {
        ekin += m_mass[i] * (m_velocities[3 * i] * m_velocities[3 * i] + m_velocities[3 * i + 1] * m_velocities[3 * i + 1] + m_velocities[3 * i + 2] * m_velocities[3 * i + 2]);
    }
    ekin *= 0.5;
    double T = 2.0 * ekin / (kb_Eh * m_dof);
    m_unstable = T > 10000 * m_T || std::isnan(T);
    m_T = T;
    ThermostatFunction();
    EKin();
}

void SimpleMD::ApplyRMSDMTD(double* grad)
{
    std::chrono::time_point<std::chrono::system_clock> m_start, m_end;
    m_start = std::chrono::system_clock::now();
    m_colvar_incr = 0;

    Geometry current_geometry = m_rmsd_mtd_molecule.getGeometry();
    for (int i = 0; i < m_rmsd_indicies.size(); ++i) {
        current_geometry(i, 0) = m_current_geometry[3 * m_rmsd_indicies[i] + 0];
        current_geometry(i, 1) = m_current_geometry[3 * m_rmsd_indicies[i] + 1];
        current_geometry(i, 2) = m_current_geometry[3 * m_rmsd_indicies[i] + 2];
    }

    double current_bias = 0;
    double rmsd_reference = 0;

    if (m_bias_structure_count == 0) {
        m_bias_threads[0]->addGeometry(current_geometry, 0, m_currentStep, 0);
        m_bias_structure_count++;
        m_rmsd_mtd_molecule.writeXYZFile(Basename() + ".mtd.xyz");
        std::ofstream colvarfile;
        colvarfile.open("COLVAR");
        colvarfile.close();
    }
    if (m_threads == 1 || m_bias_structure_count == 1) {
        for (int i = 0; i < m_bias_threads.size(); ++i) {
            m_bias_threads[i]->setCurrentGeometry(current_geometry, m_currentStep);
            m_bias_threads[i]->start();
            current_bias += m_bias_threads[i]->BiasEnergy();
            for (int j = 0; j < m_rmsd_indicies.size(); ++j) {
                grad[3 * m_rmsd_indicies[j] + 0] += m_bias_threads[i]->Gradient()(j, 0);
                grad[3 * m_rmsd_indicies[j] + 1] += m_bias_threads[i]->Gradient()(j, 1);
                grad[3 * m_rmsd_indicies[j] + 2] += m_bias_threads[i]->Gradient()(j, 2);
            }
            m_colvar_incr += m_bias_threads[i]->Counter();
            m_loop_time += m_bias_threads[i]->Time();
        }
    } else {
        if (m_bias_structure_count < m_threads) {
            for (int i = 0; i < m_bias_structure_count; ++i) {
                m_bias_threads[i]->setCurrentGeometry(current_geometry, m_currentStep);
            }
        } else {
            for (int i = 0; i < m_bias_threads.size(); ++i) {
                m_bias_threads[i]->setCurrentGeometry(current_geometry, m_currentStep);
            }
        }

        m_bias_pool->setActiveThreadCount(m_threads);
        m_bias_pool->StaticPool();
        m_bias_pool->StartAndWait();
        m_bias_pool->setWakeUp(m_bias_pool->WakeUp() / 2);

        for (int i = 0; i < m_bias_threads.size(); ++i) {
            if (m_bias_threads[i]->Return() == 1) {

                current_bias += m_bias_threads[i]->BiasEnergy();
                for (int j = 0; j < m_rmsd_indicies.size(); ++j) {
                    grad[3 * m_rmsd_indicies[j] + 0] += m_bias_threads[i]->Gradient()(j, 0);
                    grad[3 * m_rmsd_indicies[j] + 1] += m_bias_threads[i]->Gradient()(j, 1);
                    grad[3 * m_rmsd_indicies[j] + 2] += m_bias_threads[i]->Gradient()(j, 2);
                }
                m_colvar_incr += m_bias_threads[i]->Counter();
            }
            m_loop_time += m_bias_threads[i]->Time();
        }
        m_bias_pool->Reset();
    }
    rmsd_reference = m_bias_threads[0]->RMSDReference();
    if (m_nocolvarfile == false) {
        std::ofstream colvarfile;
        colvarfile.open("COLVAR", std::iostream::app);
        colvarfile << m_currentStep << " ";
        m_rmsd_mtd_molecule.setGeometry(current_geometry);
        if (m_rmsd_fragment_count < 2)
            colvarfile << rmsd_reference << " ";

        for (int i = 0; i < m_rmsd_fragment_count; ++i)
            for (int j = 0; j < i; ++j) {
                colvarfile << (m_rmsd_mtd_molecule.Centroid(true, i) - m_rmsd_mtd_molecule.Centroid(true, j)).norm() << " ";
            }
        colvarfile << current_bias << " " << std::endl;
        colvarfile.close();
    }
    m_bias_energy += current_bias;

    if (current_bias * m_rmsd_econv < m_bias_structure_count && m_rmsd_fix_structure == false) {
        int thread_index = m_bias_structure_count % m_bias_threads.size();
        m_bias_threads[thread_index]->addGeometry(current_geometry, rmsd_reference, m_currentStep, m_bias_structure_count);
        m_bias_structure_count++;
        m_rmsd_mtd_molecule.appendXYZFile(Basename() + ".mtd.xyz");
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

double SimpleMD::ApplySphericLogFermiWalls(double* grad)
{
    double potential = 0;
    double kbT = m_wall_temp * kb_Eh;
    // int counter = 0;
    for (int i = 0; i < m_natoms; ++i) {
        double distance = sqrt(m_current_geometry[3 * i + 0] * m_current_geometry[3 * i + 0] + m_current_geometry[3 * i + 1] * m_current_geometry[3 * i + 1] + m_current_geometry[3 * i + 2] * m_current_geometry[3 * i + 2]);
        double exp_expr = exp(m_wall_beta * (distance - m_wall_spheric_radius));
        double curr_pot = kbT * log(1 + exp_expr);
        // counter += distance > m_wall_radius;
        // std::cout << m_wall_beta*m_current_geometry[3 * i + 0]*exp_expr/(distance*(1-exp_expr)) << " ";
        grad[3 * i + 0] -= kbT * m_wall_beta * m_current_geometry[3 * i + 0] * exp_expr / (distance * (1 - exp_expr));
        grad[3 * i + 1] -= kbT * m_wall_beta * m_current_geometry[3 * i + 1] * exp_expr / (distance * (1 - exp_expr));
        grad[3 * i + 2] -= kbT * m_wall_beta * m_current_geometry[3 * i + 2] * exp_expr / (distance * (1 - exp_expr));

        // std::cout << distance << " ";
        potential += curr_pot;
    }
    //    std::cout << counter << " ";
    return potential;
    // std::cout << potential*kbT << std::endl;
}

double SimpleMD::ApplyRectLogFermiWalls(double* grad)
{
    double potential = 0;
    double kbT = m_wall_temp * kb_Eh;
    int counter = 0;
    double b = m_wall_beta;
    double sum_grad = 0;
    for (int i = 0; i < m_natoms; ++i) {
        double exp_expr_xl = exp(b * (m_wall_x_min - m_current_geometry[3 * i + 0]));
        double exp_expr_xu = exp(b * (m_current_geometry[3 * i + 0] - m_wall_x_max));

        double exp_expr_yl = exp(b * (m_wall_y_min - m_current_geometry[3 * i + 1]));
        double exp_expr_yu = exp(b * (m_current_geometry[3 * i + 1] - m_wall_y_max));

        double exp_expr_zl = exp(b * (m_wall_z_min - m_current_geometry[3 * i + 2]));
        double exp_expr_zu = exp(b * (m_current_geometry[3 * i + 2] - m_wall_z_max));

        double curr_pot = kbT * (log(1 + exp_expr_xl) + log(1 + exp_expr_xu) + log(1 + exp_expr_yl) + log(1 + exp_expr_yu) + log(1 + exp_expr_zl) + log(1 + exp_expr_zu));
        counter += (m_current_geometry[3 * i + 0] - m_wall_x_min) < 0 || (m_wall_x_max - m_current_geometry[3 * i + 0]) < 0 || (m_current_geometry[3 * i + 1] - m_wall_y_min) < 0 || (m_wall_y_max - m_current_geometry[3 * i + 1]) < 0 || (m_current_geometry[3 * i + 2] - m_wall_z_min) < 0 || (m_wall_z_max - m_current_geometry[3 * i + 2]) < 0;
        // std::cout << i << " " << counter << std::endl;

        // std::cout << m_wall_beta*m_current_geometry[3 * i + 0]*exp_expr/(distance*(1-exp_expr)) << " ";
        if (i == 81) {
            //    std::cout << std::endl;
            //    std::cout << m_current_geometry[3 * i + 0] << " " << m_current_geometry[3 * i + 1] << " " << m_current_geometry[3 * i + 2] << std::endl;
            //    std::cout << grad[3 * i + 0] << " " << grad[3 * i + 1] << " " <<grad[3 * i + 2] << std::endl;
        }
        grad[3 * i + 0] += kbT * b * (exp_expr_xu / (1 - exp_expr_xu) - exp_expr_xl / (1 - exp_expr_xl)); // m_current_geometry[3 * i + 0]*exp_expr/(distance*(1-exp_expr));
        grad[3 * i + 1] += kbT * b * (exp_expr_yu / (1 - exp_expr_yu) - exp_expr_yl / (1 - exp_expr_yl));
        grad[3 * i + 2] += kbT * b * (exp_expr_zu / (1 - exp_expr_zu) - exp_expr_zl / (1 - exp_expr_zl));
        sum_grad += grad[3 * i + 0] + grad[3 * i + 1] + grad[3 * i + 2];
        // if( i == 81)
        {
            // std::cout << i << " " <<grad[3 * i + 0] << " " << grad[3 * i + 1] << " " <<grad[3 * i + 2] << std::endl;
        }
        // std::cout << distance << " ";
        potential += curr_pot;
    }
    std::cout << counter << " " << sum_grad;
    return potential;
    // std::cout << potential*kbT << std::endl;
}

double SimpleMD::ApplySphericHarmonicWalls(double* grad)
{
    double potential = 0;
    double k = m_wall_temp * kb_Eh;
    int counter = 0;
    for (int i = 0; i < m_natoms; ++i) {
        double distance = sqrt(m_current_geometry[3 * i + 0] * m_current_geometry[3 * i + 0] + m_current_geometry[3 * i + 1] * m_current_geometry[3 * i + 1] + m_current_geometry[3 * i + 2] * m_current_geometry[3 * i + 2]);
        double curr_pot = 0.5 * k * (m_wall_spheric_radius - distance) * (m_wall_spheric_radius - distance) * (distance > m_wall_spheric_radius);
        double out = distance > m_wall_spheric_radius;
        counter += out;

        double diff = k * (m_wall_spheric_radius - distance) * (distance > m_wall_spheric_radius);

        double dx = diff * m_current_geometry[3 * i + 0] / distance;
        double dy = diff * m_current_geometry[3 * i + 1] / distance;
        double dz = diff * m_current_geometry[3 * i + 2] / distance;

        grad[3 * i + 0] -= dx;
        grad[3 * i + 1] -= dy;
        grad[3 * i + 2] -= dz;
        /*
        if(out)
        {
            std::cout << m_current_geometry[3 * i + 0]  << " " << m_current_geometry[3 * i + 1]  << " " << m_current_geometry[3 * i + 2] << std::endl;
            std::cout << dx << " " << dy << " " << dz << std::endl;
        }*/
        // std::cout << distance << " ";
        potential += curr_pot;
    }
    // std::cout << counter << " ";
    return potential;
    // std::cout << potential*kbT << std::endl;
}

double SimpleMD::ApplyRectHarmonicWalls(double* grad)
{
    double potential = 0;
    double k = m_wall_temp * kb_Eh;
    int counter = 0;
    double b = m_wall_beta;
    double sum_grad = 0;
    for (int i = 0; i < m_natoms; ++i) {
        double Vx = (m_current_geometry[3 * i + 0] - m_wall_x_min) * (m_current_geometry[3 * i + 0] - m_wall_x_min) * (m_current_geometry[3 * i + 0] < m_wall_x_min)
            + (m_current_geometry[3 * i + 0] - m_wall_x_max) * (m_current_geometry[3 * i + 0] - m_wall_x_max) * (m_current_geometry[3 * i + 0] > m_wall_x_max);

        double Vy = (m_current_geometry[3 * i + 1] - m_wall_y_min) * (m_current_geometry[3 * i + 1] - m_wall_y_min) * (m_current_geometry[3 * i + 1] < m_wall_y_min)
            + (m_current_geometry[3 * i + 1] - m_wall_y_max) * (m_current_geometry[3 * i + 1] - m_wall_y_max) * (m_current_geometry[3 * i + 1] > m_wall_y_max);

        double Vz = (m_current_geometry[3 * i + 2] - m_wall_z_min) * (m_current_geometry[3 * i + 2] - m_wall_z_min) * (m_current_geometry[3 * i + 2] < m_wall_z_min)
            + (m_current_geometry[3 * i + 2] - m_wall_z_max) * (m_current_geometry[3 * i + 2] - m_wall_z_max) * (m_current_geometry[3 * i + 2] > m_wall_z_max);

        double curr_pot = 0.5 * k * (Vx + Vy + Vz);
        int out = (m_current_geometry[3 * i + 0] - m_wall_x_min) < 0 || (m_wall_x_max - m_current_geometry[3 * i + 0]) < 0 || (m_current_geometry[3 * i + 1] - m_wall_y_min) < 0 || (m_wall_y_max - m_current_geometry[3 * i + 1]) < 0 || (m_current_geometry[3 * i + 2] - m_wall_z_min) < 0 || (m_wall_z_max - m_current_geometry[3 * i + 2]) < 0;
        counter += out;

        // std::cout << i << " " << counter << std::endl;

        double dx = k * (std::abs(m_current_geometry[3 * i + 0] - m_wall_x_min) * (m_current_geometry[3 * i + 0] < m_wall_x_min) - (m_current_geometry[3 * i + 0] - m_wall_x_max) * (m_current_geometry[3 * i + 0] > m_wall_x_max));

        double dy = k * (std::abs(m_current_geometry[3 * i + 1] - m_wall_y_min) * (m_current_geometry[3 * i + 1] < m_wall_y_min) - (m_current_geometry[3 * i + 1] - m_wall_y_max) * (m_current_geometry[3 * i + 1] > m_wall_y_max));

        double dz = k * (std::abs(m_current_geometry[3 * i + 2] - m_wall_z_min) * (m_current_geometry[3 * i + 2] < m_wall_z_min) - (m_current_geometry[3 * i + 2] - m_wall_z_max) * (m_current_geometry[3 * i + 2] > m_wall_z_max));
        grad[3 * i + 0] -= dx;
        grad[3 * i + 1] -= dy;
        grad[3 * i + 2] -= dz;
        /* if(out)
         {
             std::cout << m_current_geometry[3 * i + 0]  << " " << m_current_geometry[3 * i + 1]  << " " << m_current_geometry[3 * i + 2] << std::endl;
             std::cout << dx << " " << dy << " " << dz << std::endl;
         }*/
        sum_grad += dx + dy + dz;

        potential += curr_pot;
    }
    // std::cout << counter << " " << sum_grad;
    return potential;
    // std::cout << potential*kbT << std::endl;
}

void SimpleMD::RemoveRotations(std::vector<double>& velo)
{
    /*
     * This code was taken and adopted from the xtb sources
     * https://github.com/grimme-lab/xtb/blob/main/src/rmrottr.f90
     * Special thanks to the developers
     */
    double mass = 0;
    Position pos = { 0, 0, 0 }, angom{ 0, 0, 0 };
    Geometry geom(m_natoms, 3);

    std::vector<std::vector<int>> fragments = m_molecule.GetFragments();
    // std::cout << fragments.size() << std::endl;
    for (int f = 0; f < fragments.size(); ++f) {
        for (int i : fragments[f]) {
            double m = m_mass[i];
            mass += m;
            pos(0) += m * m_current_geometry[3 * i + 0];
            pos(1) += m * m_current_geometry[3 * i + 1];
            pos(2) += m * m_current_geometry[3 * i + 2];

            geom(i, 0) = m_current_geometry[3 * i + 0];
            geom(i, 1) = m_current_geometry[3 * i + 1];
            geom(i, 2) = m_current_geometry[3 * i + 2];
        }
        pos(0) /= mass;
        pos(1) /= mass;
        pos(2) /= mass;

        Geometry matrix = Geometry::Zero(3, 3);
        for (int i : fragments[f]) {
            double m = m_mass[i];
            geom(i, 0) -= pos(0);
            geom(i, 1) -= pos(1);
            geom(i, 2) -= pos(2);

            double x = geom(i, 0);
            double y = geom(i, 1);
            double z = geom(i, 2);
            angom(0) += m_mass[i] * (geom(i, 1) * velo[3 * i + 2] - geom(i, 2) * velo[3 * i + 1]);
            angom(1) += m_mass[i] * (geom(i, 2) * velo[3 * i + 0] - geom(i, 0) * velo[3 * i + 2]);
            angom(2) += m_mass[i] * (geom(i, 0) * velo[3 * i + 1] - geom(i, 1) * velo[3 * i + 0]);
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

        Position omega = matrix.inverse() * angom;

        Position rlm = { 0, 0, 0 }, ram = { 0, 0, 0 };
        for (int i : fragments[f]) {
            rlm(0) = rlm(0) + m_mass[i] * velo[3 * i + 0];
            rlm(1) = rlm(1) + m_mass[i] * velo[3 * i + 1];
            rlm(2) = rlm(2) + m_mass[i] * velo[3 * i + 2];
        }

        for (int i : fragments[f]) {
            ram(0) = (omega(1) * geom(i, 2) - omega(2) * geom(i, 1));
            ram(1) = (omega(2) * geom(i, 0) - omega(0) * geom(i, 2));
            ram(2) = (omega(0) * geom(i, 1) - omega(1) * geom(i, 0));

            velo[3 * i + 0] = velo[3 * i + 0] - rlm(0) / mass - ram(0);
            velo[3 * i + 1] = velo[3 * i + 1] - rlm(1) / mass - ram(1);
            velo[3 * i + 2] = velo[3 * i + 2] - rlm(2) / mass - ram(2);
        }
    }
}

void SimpleMD::RemoveRotation(std::vector<double>& velo)
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
        double m = m_mass[i];
        mass += m;
        pos(0) += m * m_current_geometry[3 * i + 0];
        pos(1) += m * m_current_geometry[3 * i + 1];
        pos(2) += m * m_current_geometry[3 * i + 2];

        geom(i, 0) = m_current_geometry[3 * i + 0];
        geom(i, 1) = m_current_geometry[3 * i + 1];
        geom(i, 2) = m_current_geometry[3 * i + 2];
    }
    pos(0) /= mass;
    pos(1) /= mass;
    pos(2) /= mass;

    Geometry matrix = Geometry::Zero(3, 3);
    for (int i = 0; i < m_natoms; ++i) {
        double m = m_mass[i];
        geom(i, 0) -= pos(0);
        geom(i, 1) -= pos(1);
        geom(i, 2) -= pos(2);

        double x = geom(i, 0);
        double y = geom(i, 1);
        double z = geom(i, 2);
        angom(0) += m_mass[i] * (geom(i, 1) * velo[3 * i + 2] - geom(i, 2) * velo[3 * i + 1]);
        angom(1) += m_mass[i] * (geom(i, 2) * velo[3 * i + 0] - geom(i, 0) * velo[3 * i + 2]);
        angom(2) += m_mass[i] * (geom(i, 0) * velo[3 * i + 1] - geom(i, 1) * velo[3 * i + 0]);
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

    Position omega = matrix.inverse() * angom;

    Position rlm = { 0, 0, 0 }, ram = { 0, 0, 0 };
    for (int i = 0; i < m_natoms; ++i) {
        rlm(0) = rlm(0) + m_mass[i] * velo[3 * i + 0];
        rlm(1) = rlm(1) + m_mass[i] * velo[3 * i + 1];
        rlm(2) = rlm(2) + m_mass[i] * velo[3 * i + 2];
    }

    for (int i = 0; i < m_natoms; ++i) {
        ram(0) = (omega(1) * geom(i, 2) - omega(2) * geom(i, 1));
        ram(1) = (omega(2) * geom(i, 0) - omega(0) * geom(i, 2));
        ram(2) = (omega(0) * geom(i, 1) - omega(1) * geom(i, 0));

        velo[3 * i + 0] = velo[3 * i + 0] - rlm(0) / mass - ram(0);
        velo[3 * i + 1] = velo[3 * i + 1] - rlm(1) / mass - ram(1);
        velo[3 * i + 2] = velo[3 * i + 2] - rlm(2) / mass - ram(2);
    }
}

void SimpleMD::PrintStatus() const
{
    auto unix_timestamp = std::chrono::seconds(std::time(NULL));

    int current = std::chrono::milliseconds(unix_timestamp).count();
    double duration = (current - m_unix_started) / (1000.0 * double(m_currentStep));
    double remaining;
    double tmp = (m_maxtime - m_currentStep) * duration / 60;
    if (tmp >= 1)
        remaining = tmp;
    else
        remaining = (m_maxtime - m_currentStep) * duration;
#pragma message("awfull, fix it ")
    if (m_writeUnique) {
#ifdef GCC
        std::cout << fmt::format("{1: ^{0}f} {2: ^{0}f} {3: ^{0}f} {4: ^{0}f} {5: ^{0}f} {6: ^{0}f} {7: ^{0}f} {8: ^{0}f} {9: ^{0}f} {10: ^{0}f} {11: ^{0}f} {12: ^{0}f} {13: ^{0}f} {14: ^{0}f} {15: ^{0}} {16: ^{0}}\n", 15,
            m_currentStep / 1000, m_Epot, m_aver_Epot, m_Ekin, m_aver_Ekin, m_Etot, m_aver_Etot, m_T, m_aver_Temp, m_wall_potential, m_average_wall_potential, m_virial_correction, m_average_virial_correction, remaining, m_time_step / 1000.0, m_unqiue->StoredStructures());
#else
        std::cout << m_currentStep * m_dT / fs2amu / 1000 << " " << m_Epot << " " << m_Ekin << " " << m_Epot + m_Ekin << m_T << std::endl;

#endif
    } else {
#ifdef GCC
        if (m_dipole)
            std::cout << fmt::format("{1: ^{0}f} {2: ^{0}f} {3: ^{0}f} {4: ^{0}f} {5: ^{0}f} {6: ^{0}f} {7: ^{0}f} {8: ^{0}f} {9: ^{0}f} {10: ^{0}f} {11: ^{0}f} {12: ^{0}f} {13: ^{0}f} {14: ^{0}f} {15: ^{0}f} {16: ^{0}f}\n", 15,
                m_currentStep / 1000, m_Epot, m_aver_Epot, m_Ekin, m_aver_Ekin, m_Etot, m_aver_Etot, m_T, m_aver_Temp, m_wall_potential, m_average_wall_potential, m_aver_dipol * 2.5418 * 3.3356, m_virial_correction, m_average_virial_correction, remaining, m_time_step / 1000.0);
        else
            std::cout << fmt::format("{1: ^{0}f} {2: ^{0}f} {3: ^{0}f} {4: ^{0}f} {5: ^{0}f} {6: ^{0}f} {7: ^{0}f} {8: ^{0}f} {9: ^{0}f} {10: ^{0}f} {11: ^{0}f} {12: ^{0}f} {13: ^{0}f} {14: ^{0}f} {15: ^{0}f}\n", 15,
                m_currentStep / 1000, m_Epot, m_aver_Epot, m_Ekin, m_aver_Ekin, m_Etot, m_aver_Etot, m_T, m_aver_Temp, m_wall_potential, m_average_wall_potential, m_virial_correction, m_average_virial_correction, remaining, m_time_step / 1000.0);
#else
        std::cout << m_currentStep * m_dT / fs2amu / 1000 << " " << m_Epot << " " << m_Ekin << " " << m_Epot + m_Ekin << m_T << std::endl;

#endif
    }
    //std::cout << m_mtd_time << " " << m_loop_time << std::endl;
}

void SimpleMD::PrintMatrix(const double* matrix)
{
    std::cout << "Print Matrix" << std::endl;
    for (int i = 0; i < m_natoms; ++i) {
        std::cout << matrix[3 * i] << " " << matrix[3 * i + 1] << " " << matrix[3 * i + 2] << std::endl;
    }
    std::cout << std::endl;
}

double SimpleMD::CleanEnergy(double* grad)
{
    EnergyCalculator interface(m_method, m_defaults);
    interface.setMolecule(m_molecule);
    interface.updateGeometry(m_current_geometry);

    double Energy = interface.CalculateEnergy(true);
    interface.getGradient(grad);
    if (m_dipole) {
        auto dipole = interface.Dipole();
        m_curr_dipole = sqrt(dipole[0] * dipole[0] + dipole[1] * dipole[1] + dipole[2] * dipole[2]);
        m_collected_dipole.push_back(sqrt(dipole[0] * dipole[0] + dipole[1] * dipole[1] + dipole[2] * dipole[2]));
    }
    return Energy;
}

double SimpleMD::FastEnergy(double* grad)
{
    m_interface->updateGeometry(m_current_geometry);

    double Energy = m_interface->CalculateEnergy(true);
    m_interface->getGradient(grad);
    if (m_dipole) {
        auto dipole = m_interface->Dipole();
        m_curr_dipole = sqrt(dipole[0] * dipole[0] + dipole[1] * dipole[1] + dipole[2] * dipole[2]);
        m_collected_dipole.push_back(sqrt(dipole[0] * dipole[0] + dipole[1] * dipole[1] + dipole[2] * dipole[2]));
    }
    return Energy;
}

void SimpleMD::EKin()
{
    double ekin = 0;
    for (int i = 0; i < m_natoms; ++i) {
        ekin += m_mass[i] * (m_velocities[3 * i] * m_velocities[3 * i] + m_velocities[3 * i + 1] * m_velocities[3 * i + 1] + m_velocities[3 * i + 2] * m_velocities[3 * i + 2]);
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
        m_aver_dipol = (m_curr_dipole + (m_currentStep)*m_aver_dipol) / (m_currentStep + 1);
    }
    m_average_wall_potential = (m_wall_potential + (m_currentStep)*m_average_wall_potential) / (m_currentStep + 1);
    m_average_virial_correction = (m_virial_correction + (m_currentStep)*m_average_virial_correction) / (m_currentStep + 1);
}

bool SimpleMD::WriteGeometry()
{
    bool result = true;
    Geometry geometry = m_molecule.getGeometry();
    for (int i = 0; i < m_natoms; ++i) {
        geometry(i, 0) = m_current_geometry[3 * i + 0];
        geometry(i, 1) = m_current_geometry[3 * i + 1];
        geometry(i, 2) = m_current_geometry[3 * i + 2];
    }
    TriggerWriteRestart();
    m_molecule.setGeometry(geometry);

    if (m_writeXYZ) {
        m_molecule.setEnergy(m_Epot);
        m_molecule.setName(std::to_string(m_currentStep));
        m_molecule.appendXYZFile(Basename() + ".trj.xyz");
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
        m_velocities[3 * i + 0] *= lambda;
        m_velocities[3 * i + 1] *= lambda;
        m_velocities[3 * i + 2] *= lambda;
    }
}

void SimpleMD::CSVR()
{
    double Ekin_target = 0.5 * kb_Eh * (m_T0)*m_dof;
    double c = exp(-(m_dT / 2.0 * m_respa) / m_coupling);
    static std::default_random_engine rd{};
    static std::mt19937 gen{ rd() };
    static std::normal_distribution<> d{ 0, 1 };
    static std::chi_squared_distribution<float> dchi{ m_dof };

    /*
        if(int(m_step * m_dT) % 1 == 0)
        {
        for (int i = 0; i < m_natoms; ++i) {
            m_atom_temp[i].push_back(m_mass[i]* (m_velocities[3 * i + 0] *m_velocities[3 * i + 0]
                + m_velocities[3 * i + 1] *m_velocities[3 * i + 1]
                                         + m_velocities[3 * i + 2]*m_velocities[3 * i + 2])/(kb_Eh * m_dof));
            double T = 0;
            double R = d(gen);
            double SNf = dchi(gen);
            for(int j = 0; j < m_atom_temp[i].size(); ++j)
                T += m_atom_temp[i][j];
            if(m_atom_temp[i].size() > 1000)
                m_atom_temp[i].erase(m_atom_temp[i].begin());
            double c = exp(-(T / 2.0 * m_respa) / m_coupling);
            double alpha2 = c + (1 - c) * (SNf + R * R) * Ekin_target / (m_dof * m_Ekin) + 2 * R * sqrt(c * (1 - c) * Ekin_target / (m_dof * m_Ekin));
            m_Ekin_exchange += m_Ekin * (alpha2 - 1);
            double alpha = sqrt(alpha2);
            m_velocities[3 * i + 0] *= alpha;
            m_velocities[3 * i + 1] *= alpha;
            m_velocities[3 * i + 2] *= alpha;
        }
        }else{ */
    double R = d(gen);
    double SNf = dchi(gen);
    double alpha2 = c + (1 - c) * (SNf + R * R) * Ekin_target / (m_dof * m_Ekin) + 2 * R * sqrt(c * (1 - c) * Ekin_target / (m_dof * m_Ekin));
    m_Ekin_exchange += m_Ekin * (alpha2 - 1);
    double alpha = sqrt(alpha2);
    for (int i = 0; i < m_natoms; ++i) {
        m_velocities[3 * i + 0] *= alpha;
        m_velocities[3 * i + 1] *= alpha;
        m_velocities[3 * i + 2] *= alpha;

        m_atom_temp[i].push_back(m_mass[i] * (m_velocities[3 * i + 0] * m_velocities[3 * i + 0] + m_velocities[3 * i + 1] * m_velocities[3 * i + 1] + m_velocities[3 * i + 2] * m_velocities[3 * i + 2]) / (kb_Eh * m_dof));
    }
    //    }
}

void SimpleMD::Anderson()
{
    static std::default_random_engine generator;
    double probability = m_anderson * m_dT;
    std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);
    for (size_t i = 0; i < m_natoms; ++i) {
        if (uniform_dist(generator) < probability) {
            std::normal_distribution<double> distribution(0.0, std::sqrt(kb_Eh * m_T0 * m_rmass[i]));
            m_velocities[3 * i + 0] = (m_velocities[3 * i + 0] + distribution(generator)) / 2.0;
            m_velocities[3 * i + 1] = (m_velocities[3 * i + 1] + distribution(generator)) / 2.0;
            m_velocities[3 * i + 2] = (m_velocities[3 * i + 2] + distribution(generator)) / 2.0;
        }
    }
}
void SimpleMD::NoseHover()
{
    // Berechnung der kinetischen Energie
    double kinetic_energy = 0.0;
    for (int i = 0; i < m_natoms; ++i) {
        kinetic_energy += 0.5 * m_mass[i] * (m_velocities[3 * i] * m_velocities[3 * i] + m_velocities[3 * i + 1] * m_velocities[3 * i + 1] + m_velocities[3 * i + 2] * m_velocities[3 * i + 2]);
    }
    // Update der Thermostatkette
    m_xi[0] += 0.5 * m_dT * (2.0 * kinetic_energy - m_dof * m_T0 * kb_Eh) / m_Q[0];
    for (int j = 1; j < m_chain_length; ++j) {
        m_xi[j] += 0.5 * m_dT * (m_Q[j - 1] * m_xi[j - 1] * m_xi[j - 1] - m_T0 * kb_Eh) / m_Q[j];
    }

    // Update der Geschwindigkeiten
    double scale = exp(-m_xi[0] * m_dT);
    for (int i = 0; i < m_natoms; ++i) {
        m_velocities[3 * i + 0] *= scale;
        m_velocities[3 * i + 1] *= scale;
        m_velocities[3 * i + 2] *= scale;
    }

    // Rückwärts-Update der Thermostatkette
    for (int j = m_chain_length - 1; j >= 1; --j) {
        m_xi[j] += 0.5 * m_dT * (m_Q[j - 1] * m_xi[j - 1] * m_xi[j - 1] - m_T0 * kb_Eh) / m_Q[j];
    }
    m_xi[0] += 0.5 * m_dT * (2.0 * kinetic_energy - m_dof * m_T0 * kb_Eh) / m_Q[0];
}
