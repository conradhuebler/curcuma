/*
 * <Simple MD Module for Cucuma. >
 * Copyright (C) 2023 - 2022 Conrad Hübler <Conrad.Huebler@gmx.net>
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
#include <ratio>

#include "src/capabilities/rmsdtraj.h"

#include "src/core/energycalculator.h"
#include "src/core/molecule.h"

#include "external/CxxThreadPool/include/CxxThreadPool.h"

#include "curcumamethod.h"

static json CurcumaMDJson{
    { "writeXYZ", true },
    { "printOutput", true },
    { "MaxTime", 5000 },
    { "T", 298.15 },
    { "dt", 1 },
    { "charge", 0 },
    { "Spin", 0 },
    { "centered", false },
    { "dump", 50 },
    { "print", 1000 },
    { "unique", false },
    { "rmsd", 1.5 },
    { "opt", false },
    { "hmass", 1 },
    { "velo", 1 },
    { "rescue", false },
    { "berendson", 10 },
    { "MaxTopoDiff", 15 },
    { "impuls", 0 },
    { "method", "uff" },
    { "impuls_scaling", 0.75 },
    { "writeinit", false },
    { "initfile", "none" },
    { "norestart", false },
    { "writerestart", 500 }
};

class SimpleMD : public CurcumaMethod {
public:
    SimpleMD(const json& controller, bool silent);
    ~SimpleMD();

    inline void setMolecule(const Molecule& molecule)
    {
        m_molecule = molecule;
    }

    void setBaseName(const std::string& name)
    {
        m_basename = name;
    }

    bool Initialise() override;

    void start() override;

    std::vector<Molecule*> UniqueMolecules() const { return m_unique_structures; }

private:
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

    double Gradient(const double* coord, double* grad);

    void PrintMatrix(const double* matrix);

    bool WriteGeometry();
    void Verlet(double* coord, double* grad);

    void RemoveRotation(std::vector<double>& velo);

    double EKin();
    void Berendson();

    void InitConstrainedBonds();

    std::vector<std::pair<int, int>> m_bond_constrained;
    std::string m_basename;
    int m_natoms = 0;
    int m_dumb = 1;
    double m_T = 0, m_Epot = 0, m_aver_Epot = 0, m_Ekin = 0, m_aver_Ekin = 0, m_Etot = 0, m_aver_Etot = 0;
    int m_hmass = 4;
    double m_single_step = 1;
    double m_timestep = 0.5, m_currentStep = 0, m_maxtime = 1000;
    int m_spin = 0, m_charge = 0, m_print = 100;
    double m_T0 = 298.13, m_aver_Temp = 0, m_rmsd = 1.5;
    std::vector<double> m_current_geometry, m_mass, m_velocities, m_gradient;
    std::vector<int> m_atomtype;
    Molecule m_molecule;
    bool m_initialised = false, m_restart = false, m_centered = false, m_writeUnique = true, m_opt = false, m_rescue = false, m_writeXYZ = true, m_writeinit = false, m_norestart = false;
    EnergyCalculator* m_interface;
    RMSDTraj* m_unqiue;
    const std::vector<double> m_used_mass;
    int m_unix_started = 0, m_prev_index = 0, m_max_rescue = 10, m_current_rescue = 0, m_currentTime = 0, m_max_top_diff = 15, m_step = 0;
    int m_writerestart = -1;
    double m_pos_conv = 0, m_scale_velo = 1.0, m_berendson = 1;
    double m_impuls = 0, m_impuls_scaling = 0.75, m_dt2 = 0;
    Matrix m_topo_initial;
    std::vector<Molecule*> m_unique_structures;
    std::string m_method = "UFF", m_initfile = "none";
    bool m_unstable = false;
};

class MDThread : public CxxThread {

public:
    MDThread(int thread, const json& controller)
        : m_thread(thread)
        , m_controller(controller)
    {
        setAutoDelete(true);
    }
    ~MDThread() = default;

    inline void setMolecule(const Molecule& molecule) { m_molecule = molecule; }
    SimpleMD* MDDriver() const { return m_mddriver; }

    virtual int execute() override
    {
        m_mddriver = new SimpleMD(m_controller, false);
        m_mddriver->setMolecule(m_molecule);
        m_mddriver->setBaseName("thread" + std::to_string(m_thread));
        m_mddriver->Initialise();
        m_mddriver->start();
        return 0;
    }

protected:
    int m_thread;
    Molecule m_molecule;
    std::string m_result;
    json m_controller;
    SimpleMD* m_mddriver;
};
