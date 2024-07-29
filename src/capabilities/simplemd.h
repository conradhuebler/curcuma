/*
 * <Simple MD Module for Cucuma. >
 * Copyright (C) 2023 - 2024 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "src/capabilities/rmsdtraj.h"

#include "src/core/energycalculator.h"
#include "src/core/molecule.h"

#include "external/CxxThreadPool/include/CxxThreadPool.h"

#include "curcumamethod.h"

struct BiasStructure {
    Geometry geometry;
    double time = 0;
    double rmsd_reference = 0;
    double energy = 0;
    int counter = 0;
};

static json CurcumaMDJson{
    { "writeXYZ", true },
    { "printOutput", true },
    { "MaxTime", 5000 },
    { "T", 298.15 },
    { "dt", 1 }, // single step in fs
    { "rm_COM", 100 }, // remove translation and rotation every x fs
    { "charge", 0 },
    { "Spin", 0 },
    { "rmrottrans", 0 },
    { "nocenter", false },
    { "dump", 50 },
    { "print", 1000 },
    { "unique", false },
    { "rmsd", 1.5 },
    { "opt", false },
    { "hmass", 1 },
    { "velo", 1 },
    { "rescue", false },
    { "coupling", 10 },
    { "MaxTopoDiff", 15 },
    { "impuls", 0 },
    { "method", "uff" },
    { "impuls_scaling", 0.75 },
    { "writeinit", false },
    { "initfile", "none" },
    { "norestart", false },
    { "writerestart", 1000 },
    { "rattle", false },
    { "rattle_tolerance", 1e-6 },
    { "rattle_maxiter", 10 },
    { "thermostat", "csvr" },
    { "respa", 1 },
    { "dipole", false },
    { "seed", 1 },
    { "cleanenergy", false },
    { "wall", "none" }, // can be spheric or rect
    { "wall_type", "harmonic" }, // can be logfermi or harmonic
    { "wall_spheric_radius", 0 },
    { "wall_xl", 0 },
    { "wall_yl", 0 },
    { "wall_zl", 0 },
    { "wall_x_min", 0 },
    { "wall_x_max", 0 },
    { "wall_y_min", 0 },
    { "wall_y_max", 0 },
    { "wall_z_min", 0 },
    { "wall_z_max", 0 },
    { "wall_temp", 298.15 },
    { "wall_beta", 6 },
    { "mtd", false },
    { "plumed", "plumed.dat" },
    { "mtd_dT", -1 },
    { "rmsd_mtd", false },
    { "k_rmsd", 0.1 },
    { "alpha_rmsd", 10 },
    { "rmsd_rmsd", 1 },
    { "mtd_steps", 1 },
    { "max_rmsd_N", -1 },
    { "multi_rmsd", 1e8 }
};

class SimpleMD : public CurcumaMethod {
public:
    SimpleMD(const json& controller, bool silent);
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

    double FastEnergy(double* grad);
    double CleanEnergy(double* grad);

    void PrintMatrix(const double* matrix);

    bool WriteGeometry();
    void Verlet(double* grad);
    void Rattle(double* grad);
    void ApplyRMSDMTD(double* grad);

    void Rattle_Verlet_First(double* coord, double* grad);
    void Rattle_Constrain_First(double* coord, double* grad);
    void Rattle_Verlet_Second(double* coord, double* grad);
    void Rattle_Constrain_Second(double* coord, double* grad);

    void RemoveRotation(std::vector<double>& velo);
    void RemoveRotations(std::vector<double>& velo);

    double EKin();
    void Berendson();
    void CSVR();
    void None();

    void InitialiseWalls();

    double ApplySphericLogFermiWalls(double* grad);
    double ApplyRectLogFermiWalls(double* grad);

    double ApplySphericHarmonicWalls(double* grad);
    double ApplyRectHarmonicWalls(double* grad);

    void InitConstrainedBonds();

    std::function<void(double* grad)> Integrator;
    std::function<double(double* grad)> Energy;
    std::function<double(double* grad)> WallPotential;

    std::vector<std::pair<std::pair<int, int>, double>> m_bond_constrained;

    int m_natoms = 0;
    int m_dump = 1;
    double m_T = 0, m_Epot = 0, m_aver_Epot = 0, m_Ekin = 0, m_aver_Ekin = 0, m_Etot = 0, m_aver_Etot = 0, m_aver_dipol = 0, m_curr_dipole = 0;
    double m_rm_COM = 100;
    int m_rm_COM_step = -1;
    int m_hmass = 1;
    double m_single_step = 1;
    double m_dT = 0.5, m_currentStep = 0, m_maxtime = 1000;
    int m_spin = 0, m_charge = 0, m_print = 100;
    double m_T0 = 298.13, m_aver_Temp = 0, m_rmsd = 1.5;
    double m_x0 = 0, m_y0 = 0, m_z0 = 0;
    double m_Ekin_exchange = 0.0;
    std::vector<double> m_current_geometry, m_mass, m_velocities, m_gradient, m_rmass, m_virial, m_gradient_bias;
    std::vector<int> m_atomtype;
    Molecule m_molecule, m_reference, m_target;
    bool m_initialised = false, m_restart = false, m_writeUnique = true, m_opt = false, m_rescue = false, m_writeXYZ = true, m_writeinit = false, m_norestart = false;
    int m_rmrottrans = 0, m_rattle_maxiter = 100;
    bool m_nocenter = false;
    EnergyCalculator* m_interface;
    RMSDTraj* m_unqiue;
    const std::vector<double> m_used_mass;
    std::vector<Geometry> m_bias_structures;
    std::vector<BiasStructure> m_biased_structures;

    int m_unix_started = 0, m_prev_index = 0, m_max_rescue = 10, m_current_rescue = 0, m_currentTime = 0, m_max_top_diff = 15, m_step = 0;
    int m_writerestart = -1;
    int m_respa = 1;
    double m_pos_conv = 0, m_scale_velo = 1.0, m_coupling = 10;
    double m_impuls = 0, m_impuls_scaling = 0.75, m_dt2 = 0;
    double m_rattle_tolerance = 1;
    double m_wall_spheric_radius = 6, m_wall_temp = 298.15, m_wall_beta = 6;
    double m_wall_x_min = 0, m_wall_x_max = 0, m_wall_y_min = 0, m_wall_y_max = 0, m_wall_z_min = 0, m_wall_z_max = 0;
    double m_wall_potential = 0, m_average_wall_potential = 0;
    double m_virial_correction = 0, m_average_virial_correction = 0;
    double m_deltaT = 0;
    double m_k_rmsd = 0.04;
    double m_alpha_rmsd = 0.5;
    double m_bias_energy = 0;
    double m_rmsd_rmsd = 1;
    double m_mult_rmsd = 1e4;
    int m_max_rmsd_N = -1;
    int m_mtd_steps = 10;
    int m_rattle = 0;
    int m_colvar_incr = 0;
    std::vector<double> m_collected_dipole;
    Matrix m_topo_initial;
    std::vector<Molecule*> m_unique_structures;
    std::string m_method = "UFF", m_initfile = "none", m_thermostat = "csvr", m_plumed;
    bool m_unstable = false;
    bool m_dipole = false;
    bool m_clean_energy = false;
    bool m_mtd = false;
    bool m_eval_mtd = true;
    bool m_rmsd_mtd = false;
    int m_mtd_dT = -1;
    int m_seed = -1;
    int m_time_step = 0;
    int m_dof = 0;
};

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
        m_mddriver->overrideBasename(m_basename + ".t" + std::to_string(ThreadId()));
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
