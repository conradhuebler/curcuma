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

#include "simplemd.h"

const int method = 66;

SimpleMD::SimpleMD()
{
    m_interface = new TBLiteInterface;
}

SimpleMD::~SimpleMD()
{
}

void SimpleMD::Initialise()
{
    if (m_molecule.AtomCount() == 0)
        return;

    std::ofstream result_file;
    result_file.open("trajectory.xyz");
    result_file.close();

    m_natoms = m_molecule.AtomCount();
    m_molecule.setCharge(0);

    m_timestep = 2; // * fs2amu;
    m_maxsteps = 20000; // 1e4 * fs2amu;
    m_current_geometry = std::vector<double>(3 * m_natoms, 0);
    m_velocities = std::vector<double>(3 * m_natoms, 0);
    m_mass = std::vector<double>(3 * m_natoms, 0);
    m_atomtype = std::vector<int>(m_natoms, 0);

    srand(time(NULL));
    for (int i = 0; i < m_natoms; ++i) {
        Position pos = m_molecule.Atom(i).second;
        m_atomtype[i] = m_molecule.Atom(i).first;
        m_current_geometry[3 * i + 0] = pos(0) / au;
        m_current_geometry[3 * i + 1] = pos(1) / au;
        m_current_geometry[3 * i + 2] = pos(2) / au;
        double factor = sqrt(kb * m_temperatur / (Elements::AtomicMass[m_atomtype[i]]) * amu2au) / m_natoms;
        m_velocities[3 * i + 0] = ((rand() / double(RAND_MAX) * 2 - 1) - 0.5) * factor;
        m_velocities[3 * i + 1] = ((rand() / double(RAND_MAX) * 2 - 1) - 0.5) * factor;
        m_velocities[3 * i + 2] = ((rand() / double(RAND_MAX) * 2 - 1) - 0.5) * factor;
        m_mass[3 * i + 0] = Elements::AtomicMass[m_atomtype[i]] * amu2au;
        m_mass[3 * i + 1] = Elements::AtomicMass[m_atomtype[i]] * amu2au;
        m_mass[3 * i + 2] = Elements::AtomicMass[m_atomtype[i]] * amu2au;
    }
    m_interface->InitialiseMolecule(m_molecule);
    m_initialised = true;
}

void SimpleMD::Dance()
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

    double energy = 0;
    double ekin = 0;
    for (int step = 0; step < m_maxsteps; ++step) {
        energy = Gradient(coord, gradient_prev);
        WriteGeometry();

        UpdatePosition(gradient_prev, coord);
        energy = Gradient(coord, gradient_current);
        m_molecule.setEnergy(energy);

        UpdateVelocities(gradient_prev, gradient_current);
        ekin = EKin();
        Thermostat();
        std::cout << energy << " " << ekin << " " << m_curr_temp << " " << step << std::endl;
    }
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

    double Energy = m_interface->GFNCalculation(method, grad);
    return Energy;
}

double SimpleMD::EKin()
{
    double ekin = 0;
    for (int i = 0; i < 3 * m_natoms; ++i) {
        ekin += m_mass[i] * m_velocities[i] * m_velocities[i];
    }
    m_curr_temp = 2.0 * ekin / (kb * 3 * m_natoms);
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
    m_molecule.appendXYZFile("trajectory.xyz");
}

void SimpleMD::Thermostat()
{
    double lambda = sqrt(m_temperatur / m_curr_temp);
    for (int i = 0; i < 3 * m_natoms; ++i) {
        m_velocities[i] *= lambda;
    }
}
