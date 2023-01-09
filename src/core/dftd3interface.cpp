/*
 * < C++ XTB and tblite Interface >
 * Copyright (C) 2020 - 2023 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "s-dftd3.h"

#include "src/core/global.h"
#include "src/tools/general.h"

#include "src/core/molecule.h"

#include <iostream>
#include <math.h>
#include <stdio.h>

#include "dftd3interface.h"

DFTD3Interface::DFTD3Interface(const json& controller)
{
    json parameter = MergeJson(DFTD3Settings, controller);
    m_d3_a1 = parameter["d_a1"];
    m_d3_a2 = parameter["d_a2"];
    m_d3_alp = parameter["d_alp"];

    m_d3_s6 = parameter["d_s6"];
    m_d3_s8 = parameter["d_s8"];
    m_d3_s9 = parameter["d_s9"];

    m_param = dftd3_new_rational_damping(m_error, m_d3_s6, m_d3_s8, m_d3_s9, m_d3_a1, m_d3_a2, m_d3_alp);

    m_error = dftd3_new_error();

    PrintParameter();
}

DFTD3Interface::~DFTD3Interface()
{
    dftd3_delete_param(&m_param);
    dftd3_delete_model(&m_disp);
    dftd3_delete_structure(&m_mol);
    dftd3_delete_error(&m_error);
}

void DFTD3Interface::PrintParameter() const
{
    // std::cout << m_d3_s6 << " " << m_d3_s8 << " " << m_d3_s9 << " " << m_d3_a1 << " " << m_d3_a2 << " " << m_d3_alp << std::endl;
}

void DFTD3Interface::UpdateParameters(const json& controller)
{
    json parameter = MergeJson(DFTD3Settings, controller);
    m_d3_a1 = parameter["d_a1"];
    m_d3_a2 = parameter["d_a2"];
    m_d3_alp = parameter["d_alp"];

    m_d3_s6 = parameter["d_s6"];
    m_d3_s8 = parameter["d_s8"];

    m_d3_s9 = parameter["d_s9"];

    dftd3_delete_param(&m_param);
    m_param = dftd3_new_rational_damping(m_error, m_d3_s6, m_d3_s8, m_d3_s9, m_d3_a1, m_d3_a2, m_d3_alp);

    PrintParameter();
}

bool DFTD3Interface::InitialiseMolecule(const std::vector<int>& atomtypes)
{
    int const natoms = atomtypes.size();

    int attyp[natoms];
    double coord[3 * natoms];
    m_geom = std::vector<double>(3 * natoms, 0);
    for (int i = 0; i < natoms; ++i) {
        coord[3 * i + 0] = 0 / au;
        coord[3 * i + 1] = 0 / au;
        coord[3 * i + 2] = 0 / au;
        attyp[i] = atomtypes[i];
    }

    m_mol = dftd3_new_structure(m_error, natoms, attyp, coord, NULL, NULL);
    m_disp = dftd3_new_d3_model(m_error, m_mol);

    return true;
}

double DFTD3Interface::DFTD3Calculation(double* grad)
{
    double energy = 0;
    double coord[m_geom.size()];
    for (int i = 0; i < m_geom.size(); ++i) {
        coord[i] = m_geom[i] / au;
    }
    double sigma[9];
    // std::cout << m_geom.size() << std::endl;
    // dftd3_delete_param(&m_param);
    // m_param = dftd3_load_zero_damping(m_error, "pbe", false);

    dftd3_update_structure(m_error, m_mol, coord, NULL);
    dftd3_get_dispersion(m_error, m_mol, m_disp, m_param, &energy, grad, sigma);
    return energy;
}

void DFTD3Interface::clear()
{
    dftd3_delete_param(&m_param);
    dftd3_delete_model(&m_disp);
    dftd3_delete_structure(&m_mol);
    dftd3_delete_error(&m_error);
}
