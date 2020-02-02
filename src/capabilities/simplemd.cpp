/*
 * <Simple MD Module for Cucuma. >
 * Copyright (C) 2020 Conrad Hübler <Conrad.Huebler@gmx.net>
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
#include "src/core/xtbinterface.h"

#include "simplemd.h"

SimpleMD::SimpleMD()
{
    m_interface = new XTBInterface;
}

SimpleMD::~SimpleMD()
{
}

void SimpleMD::Initialise()
{
    if (m_molecule.AtomCount() == 0)
        return;

    m_natoms = m_molecule.AtomCount();
    m_timestep = 0.01 * fs2amu;
    m_maxsteps = 1e4 * fs2amu;
    m_current_geometry = std::vector<double>(3 * m_natoms, 0);
    m_velocities = std::vector<double>(3 * m_natoms, 0);
    m_mass = std::vector<double>(3 * m_natoms, 0);
    m_atomtype = std::vector<int>(m_natoms, 0);
    // d
    // double dist_gradient[3 * m_natoms];
    srand(time(NULL));
    double temp = sqrt(m_temperatur);
    for (int i = 0; i < m_natoms; ++i) {
        Position pos = m_molecule.Atom(i).second;
        m_atomtype[i] = m_molecule.Atom(i).first;
        m_current_geometry[3 * i + 0] = pos(0) / au;
        m_current_geometry[3 * i + 1] = pos(1) / au;
        m_current_geometry[3 * i + 2] = pos(2) / au;
        m_velocities[3 * i + 0] = (rand() / double(RAND_MAX) * 2 - 1);
        m_velocities[3 * i + 1] = (rand() / double(RAND_MAX) * 2 - 1);
        m_velocities[3 * i + 2] = (rand() / double(RAND_MAX) * 2 - 1);
        m_mass[3 * i + 0] = Elements::AtomicMass[m_atomtype[i]] * amu2au;
        m_mass[3 * i + 1] = Elements::AtomicMass[m_atomtype[i]] * amu2au;
        m_mass[3 * i + 2] = Elements::AtomicMass[m_atomtype[i]] * amu2au;
    }

    m_initialised = true;
}

void SimpleMD::Dance()
{
    if (m_initialised == false)
        return;

    double coord[3 * m_natoms];
    for (int i = 0; i < 3 * m_natoms; ++i)
        coord[i] = m_current_geometry[i];

    int attyp[m_natoms];
    for (int i = 0; i < m_natoms; ++i)
        attyp[i] = m_atomtype[i];

    WriteGeometry();
    double gradient[3 * m_natoms];
    double energy = Gradient(attyp, coord, gradient);
    // EKin();
    // Thermostat();
    double ekin = 0;
    for (int step = 0; step < m_maxsteps; ++step) {
        std::cout << energy << " " << ekin << " " << m_curr_temp << " " << step / fs2amu << std::endl;
        for (int i = 0; i < 3 * m_natoms; ++i) {
            m_velocities[i] += +0.5 * m_timestep * gradient[i] / m_mass[i];
            m_current_geometry[i] += m_timestep * m_velocities[i];
            coord[i] = m_current_geometry[i];
        }
        energy = Gradient(attyp, coord, gradient);
        EKin();
        for (int i = 0; i < 3 * m_natoms; ++i) {
            m_velocities[i] += +0.5 * m_timestep * gradient[i] / m_mass[i];
        }
        WriteGeometry();
        Thermostat();
        ekin = EKin();
    }
}

double SimpleMD::Gradient(const int* attyp, const double* coord, double* grad)
{
    double Energy = m_interface->GFN2Energy(attyp, coord, m_natoms, m_charge, grad);
    return Energy;
}

void SimpleMD::Propagate()
{
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
        m_atomtype[i] = m_molecule.Atom(i).first;
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
        m_velocities[i] *= m_velocities[i] * lambda;
    }
}
