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

#pragma once

#include "src/core/molecule.h"
#include "src/core/xtbinterface.h"

class SimpleMD {
public:
    SimpleMD();
    ~SimpleMD();
    inline void setMolecule(const Molecule& molecule)
    {
        m_molecule = molecule;
    }
    void Initialise();

    void Dance();

private:
    double Gradient(const int* attyp, const double* coord, double* grad);
    void Propagate();
    void WriteGeometry();
    double EKin();
    void Thermostat();
    int m_natoms = 0;
    double m_curr_temp = 0;
    double m_timestep = 0.5;
    int m_maxsteps = 100;
    double m_charge = 0.0;
    double m_temperatur = 298.13;
    std::vector<double> m_current_geometry, m_mass, m_velocities, m_gradient;
    std::vector<int> m_atomtype;
    Molecule m_molecule;
    bool m_initialised = false;
    XTBInterface* m_interface;
};
