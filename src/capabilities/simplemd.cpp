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
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "src/core/elements.h"
#include "src/core/global.h"
#include "src/core/molecule.h"
#include "src/core/tbliteinterface.h"

#include "src/tools/geometry.h"

#include "simplemd.h"

SimpleMD::SimpleMD(const json& controller, bool silent)
    : CurcumaMethod(CurcumaMDJson, controller, silent)
{
    UpdateController(controller);
    m_interface = new TBLiteInterface;
}

SimpleMD::~SimpleMD()
{
}

void SimpleMD::LoadControlJson()
{
    m_gfn = Json2KeyWord<int>(m_defaults, "GFN");
    m_spin = Json2KeyWord<int>(m_defaults, "spin");
    m_charge = Json2KeyWord<int>(m_defaults, "charge");
    m_timestep = Json2KeyWord<int>(m_defaults, "dT"); // * fs2amu;
    m_maxsteps = Json2KeyWord<int>(m_defaults, "MaxSteps"); // 1e4 * fs2amu;
    m_temperatur = Json2KeyWord<double>(m_defaults, "T");
    m_centered = Json2KeyWord<bool>(m_defaults, "centered");
}

bool SimpleMD::Initialise()
{
    LoadRestartInformation();

    if (m_molecule.AtomCount() == 0)
        return false;

    if (!m_restart) {
        std::ofstream result_file;
        result_file.open(m_basename + ".trj.xyz");
        result_file.close();
    }
    m_natoms = m_molecule.AtomCount();
    m_molecule.setCharge(0);

    m_mass = std::vector<double>(3 * m_natoms, 0);
    m_atomtype = std::vector<int>(m_natoms, 0);

    srand(time(NULL));
    if (!m_restart) {
        m_current_geometry = std::vector<double>(3 * m_natoms, 0);
        m_velocities = std::vector<double>(3 * m_natoms, 0);
        m_currentStep = 0;
    } else {
        //   std::cout << Tools::DoubleVector2String(m_current_geometry)  << std::endl;
        //   std::cout << Tools::DoubleVector2String(m_velocities)  << std::endl;
    }

    for (int i = 0; i < m_natoms; ++i) {
        m_atomtype[i] = m_molecule.Atom(i).first;

        if (!m_restart) {
            Position pos = m_molecule.Atom(i).second;
            m_current_geometry[3 * i + 0] = pos(0) / au;
            m_current_geometry[3 * i + 1] = pos(1) / au;
            m_current_geometry[3 * i + 2] = pos(2) / au;
            double factor = sqrt(kb * m_temperatur / (Elements::AtomicMass[m_atomtype[i]])) / m_natoms;
            m_velocities[3 * i + 0] = ((rand() / double(RAND_MAX) * 2 - 1) - 0.5) * factor;
            m_velocities[3 * i + 1] = ((rand() / double(RAND_MAX) * 2 - 1) - 0.5) * factor;
            m_velocities[3 * i + 2] = ((rand() / double(RAND_MAX) * 2 - 1) - 0.5) * factor;
        }

        m_mass[3 * i + 0] = Elements::AtomicMass[m_atomtype[i]] * amu2au;
        m_mass[3 * i + 1] = Elements::AtomicMass[m_atomtype[i]] * amu2au;
        m_mass[3 * i + 2] = Elements::AtomicMass[m_atomtype[i]] * amu2au;
    }
    m_molecule.setCharge(m_charge);
    m_molecule.setSpin(m_spin);
    m_interface->InitialiseMolecule(m_molecule);
    m_initialised = true;
    return true;
}

nlohmann::json SimpleMD::WriteRestartInformation()
{
    nlohmann::json restart;
    restart["GFN"] = m_gfn;
    restart["dT"] = m_timestep;
    restart["MaxSteps"] = m_maxsteps;
    restart["T"] = m_temperatur;
    restart["currentStep"] = m_currentStep;
    restart["velocities"] = Tools::DoubleVector2String(m_velocities);
    restart["geometry"] = Tools::DoubleVector2String(m_current_geometry);
    restart["centered"] = m_centered;
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
        std::string geometry, velocities;

        try {
            m_gfn = md["GFN"];
        } catch (json::type_error& e) {
        }
        try {
            m_timestep = md["dT"];
        } catch (json::type_error& e) {
        }
        try {
            m_maxsteps = md["MaxSteps"];
        } catch (json::type_error& e) {
        }
        try {
            m_centered = md["centered"];
        } catch (json::type_error& e) {
        }
        try {
            m_temperatur = md["T"];
        } catch (json::type_error& e) {
        }
        try {
            m_currentStep = md["currentStep"];
        } catch (json::type_error& e) {
        }
        try {
            geometry = md["geometry"];
        } catch (json::type_error& e) {
        }
        try {
            velocities = md["velocities"];
        } catch (json::type_error& e) {
        }
        if (geometry.size()) {
            m_current_geometry = Tools::String2DoubleVec(geometry);
        }
        if (velocities.size()) {
            m_velocities = Tools::String2DoubleVec(velocities);
        }
        m_restart = geometry.size() && velocities.size();
    }
    return true;
};

void SimpleMD::start()
{
    if (m_initialised == false)
        return;

    double coord[3 * m_natoms];
    double gradient_prev[3 * m_natoms], gradient_current[3 * m_natoms];

    for (int i = 0; i < 3 * m_natoms; ++i) {
        coord[i] = m_current_geometry[i];
        gradient_prev[i] = 0;
        gradient_current[i] = 0;
    }

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
    for (; m_currentStep < m_maxsteps; ++m_currentStep) {
        if (CheckStop() == true) {
            TriggerWriteRestart();
            return;
        }
        Gradient(coord, gradient_prev);

        WriteGeometry();

        UpdatePosition(gradient_prev, coord);
        m_Epot = Gradient(coord, gradient_current);
        m_molecule.setEnergy(m_Epot);

        UpdateVelocities(gradient_prev, gradient_current);
        m_Ekin = EKin();
        PrintStatus();
        Thermostat();
    }
}

void SimpleMD::PrintStatus() const
{
#ifdef GCC
    std::cout << fmt::format("{1: ^{0}} {2: ^{0}f} {3: ^{0}f} {4: ^{0}f} {5: ^{0}f}\n", 15, m_currentStep, m_Epot, m_Ekin, m_Epot + m_Ekin, m_T);
#else
    std::cout << m_currentStep << " " << m_Epot << " " << m_Ekin << " " << m_Epot + m_Ekin << m_T << std::endl;

#endif
}

void SimpleMD::PrintMatrix(const double* matrix)
{
    std::cout << "Print Matrix" << std::endl;
    for (int i = 0; i < m_natoms; ++i) {
        std::cout << matrix[3 * i] << " " << matrix[3 * i + 1] << " " << matrix[3 * i + 2] << std::endl;
    }
    std::cout << std::endl;
}

void SimpleMD::UpdatePosition(const double* gradient, double* coord)
{
    for (int i = 0; i < 3 * m_natoms; ++i) {
        m_current_geometry[i] -= m_timestep * m_velocities[i] + gradient[i] / (2 * m_mass[i]);
        coord[i] = m_current_geometry[i];
    }
}

void SimpleMD::UpdateVelocities(const double* gradient_prev, const double* gradient_curr)
{
    for (int i = 0; i < 3 * m_natoms; ++i) {
        m_velocities[i] += m_timestep * ((gradient_prev[i] + gradient_curr[i]) / (4 * m_mass[i]));
    }
}

double SimpleMD::Gradient(const double* coord, double* grad)
{
    m_interface->UpdateMolecule(coord);

    double Energy = m_interface->GFNCalculation(m_gfn, grad);
    return Energy;
}

double SimpleMD::EKin()
{
    double ekin = 0;
    for (int i = 0; i < 3 * m_natoms; ++i) {
        ekin += m_mass[i] * m_velocities[i] * m_velocities[i];
    }
    m_T = 2.0 * ekin / (kb * 3 * m_natoms);
    return ekin;
}

void SimpleMD::WriteGeometry()
{
    Geometry geometry = m_molecule.getGeometry();
    for (int i = 0; i < m_natoms; ++i) {
        Position pos = m_molecule.Atom(i).second;
        geometry(i, 0) = m_current_geometry[3 * i + 0] * au;
        geometry(i, 1) = m_current_geometry[3 * i + 1] * au;
        geometry(i, 2) = m_current_geometry[3 * i + 2] * au;
    }
    m_molecule.setGeometry(geometry);
    if (m_centered)
        m_molecule.setGeometry(GeometryTools::TranslateGeometry(m_molecule.getGeometry(), GeometryTools::Centroid(m_molecule.getGeometry()), Position{ 0, 0, 0 }));

    m_molecule.appendXYZFile(m_basename + ".trj.xyz");
}

void SimpleMD::Thermostat()
{
    double lambda = sqrt(m_temperatur / m_T);
    for (int i = 0; i < 3 * m_natoms; ++i) {
        m_velocities[i] *= lambda;
    }
}
