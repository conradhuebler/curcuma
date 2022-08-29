/*
 * <Simple MD Module for Cucuma. >
 * Copyright (C) 2020 - 2022 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "src/core/molecule.h"
#include "src/core/tbliteinterface.h"

#include "curcumamethod.h"

static json CurcumaMDJson{
    { "writeXYZ", true },
    { "printOutput", true },
    { "GFN", 2 },
    { "MaxSteps", 5000 },
    { "T", 298.15 },
    { "dt", 1 },
    { "charge", 0 },
    { "Spin", 0 },
    { "centered", false },
    { "dump", 50 },
    { "print", 1000 },
    { "unique", true },
    { "rmsd", 1.5 },
    { "opt", false },
    { "hmass", 4 },
    { "velo", 1 },
    { "rescue", false },
    { "thermostat_steps", 1 },
    { "berendson", 1 }
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

private:
    void PrintStatus() const;

    /* Lets have this for all modules */
    virtual nlohmann::json WriteRestartInformation();

    /* Lets have this for all modules */
    virtual bool LoadRestartInformation();

    bool LoadRestartInformation(const json& state);

    virtual StringList MethodName() const
    {
        return { "MD" };
    }

    /* Lets have all methods read the input/control file */
    virtual void ReadControlFile()
    {
    }

    /* Read Controller has to be implemented for all */
    virtual void LoadControlJson();

    void InitVelocities(double scaling = 1.0);

    double Gradient(const double* coord, double* grad);
    void UpdatePosition(const double* grad, double* coord);
    void UpdateVelocities(double* gradient_prev, const double* gradient_curr);

    void PrintMatrix(const double* matrix);

    bool WriteGeometry();
    void SimpleIntegrator(double* coord, double* grad_prev, double* grad_next);
    void VelocityVerlet(double* coord, double* grad_prev, double* grad_next);

    double EKin();
    void Berendson();

    std::string m_basename;
    int m_natoms = 0;
    int m_gfn = 2;
    int m_dumb = 1;
    double m_T = 0, m_Epot = 0, m_aver_Epot = 0, m_Ekin = 0, m_aver_Ekin = 0, m_Etot = 0, m_aver_Etot = 0;
    int m_single_step = 1, m_hmass = 4;
    double m_timestep = 0.5;
    int m_maxsteps = 10, m_currentStep = 0, m_spin = 0, m_charge = 0, m_print = 100;
    double m_temperatur = 298.13, m_aver_Temp = 0, m_rmsd = 1.5;
    std::vector<double> m_current_geometry, m_mass, m_velocities, m_gradient;
    std::vector<int> m_atomtype;
    Molecule m_molecule;
    bool m_initialised = false, m_restart = false, m_centered = false, m_writeUnique = true, m_opt = false, m_rescue = false, m_writeXYZ = true;
    TBLiteInterface* m_interface;
    RMSDTraj* m_unqiue;
    const std::vector<double> m_used_mass;
    int m_unix_started = 0, m_prev_index = 0, m_thermostat_steps = 10, m_max_rescue = 10, m_current_rescue = 0;
    double m_pos_conv = 0, m_scale_velo = 1.0, m_berendson = 1;
};
