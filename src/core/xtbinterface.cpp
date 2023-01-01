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

#include "src/core/global.h"
#include "src/tools/general.h"

#include "external/xtb/include/xtb.h"

#include "src/core/molecule.h"

#include <iostream>
#include <math.h>
#include <stdio.h>

#include "xtbinterface.h"

XTBInterface::XTBInterface(const json& xtbsettings)
    : m_xtbsettings(xtbsettings)
{
    m_env = xtb_newEnvironment();
    m_xtb_calc = xtb_newCalculator();
    m_xtb_res = xtb_newResults();
}

XTBInterface::~XTBInterface()
{
    xtb_delResults(&m_xtb_res);
    xtb_delCalculator(&m_xtb_calc);
    xtb_delMolecule(&m_xtb_mol);
    xtb_delEnvironment(&m_env);
}

bool XTBInterface::InitialiseMolecule(const Molecule& molecule)
{
    if (m_initialised)
        UpdateMolecule(molecule);
    int const natoms = molecule.AtomCount();

    int attyp[natoms];
    std::vector<int> atoms = molecule.Atoms();
    double coord[3 * natoms];

    for (int i = 0; i < atoms.size(); ++i) {
        std::pair<int, Position> atom = molecule.Atom(i);
        coord[3 * i + 0] = atom.second(0) / au;
        coord[3 * i + 1] = atom.second(1) / au;
        coord[3 * i + 2] = atom.second(2) / au;
        attyp[i] = atoms[i];
    }
    return InitialiseMolecule(attyp, coord, natoms, molecule.Charge(), molecule.Spin());
}

bool XTBInterface::InitialiseMolecule(const Molecule* molecule)
{
    if (m_initialised)
        UpdateMolecule(molecule);
    int const natoms = molecule->AtomCount();

    int attyp[natoms];
    std::vector<int> atoms = molecule->Atoms();
    double coord[3 * natoms];

    for (int i = 0; i < atoms.size(); ++i) {
        std::pair<int, Position> atom = molecule->Atom(i);
        coord[3 * i + 0] = atom.second(0) / au;
        coord[3 * i + 1] = atom.second(1) / au;
        coord[3 * i + 2] = atom.second(2) / au;
        attyp[i] = atoms[i];
    }
    return InitialiseMolecule(attyp, coord, natoms, molecule->Charge(), molecule->Spin());
}

bool XTBInterface::InitialiseMolecule(const int* attyp, const double* coord, const int natoms, const double charge, const int spin)
{
    if (m_initialised)
        UpdateMolecule(coord);

    m_xtb_mol = xtb_newMolecule(m_env, &natoms, attyp, coord, &charge, &spin, NULL, NULL);

    m_initialised = true;
    return true;
}

bool XTBInterface::UpdateMolecule(const Molecule& molecule)
{
    int const natoms = molecule.AtomCount();
    double coord[3 * natoms];

    for (int i = 0; i < natoms; ++i) {
        std::pair<int, Position> atom = molecule.Atom(i);
        coord[3 * i + 0] = atom.second(0) / au;
        coord[3 * i + 1] = atom.second(1) / au;
        coord[3 * i + 2] = atom.second(2) / au;
    }

    return UpdateMolecule(coord);
}

bool XTBInterface::UpdateMolecule(const double* coord)
{
    xtb_updateMolecule(m_env, m_xtb_mol, coord, NULL);
    return true;
}

double XTBInterface::GFNCalculation(int parameter, double* grad)
{
    double energy = 0;
    xtb_setVerbosity(m_env, XTB_VERBOSITY_MUTED);
    if (parameter == 0) {
        xtb_loadGFN0xTB(m_env, m_xtb_mol, m_xtb_calc, NULL);
    } else if (parameter == 1) {
        xtb_loadGFN1xTB(m_env, m_xtb_mol, m_xtb_calc, NULL);
    } else if (parameter == 2) {
        xtb_loadGFN2xTB(m_env, m_xtb_mol, m_xtb_calc, NULL);
    } else if (parameter == 66) {
        xtb_loadGFNFF(m_env, m_xtb_mol, m_xtb_calc, NULL);
    }
    xtb_singlepoint(m_env, m_xtb_mol, m_xtb_calc, m_xtb_res);
    if (xtb_checkEnvironment(m_env)) {
        xtb_showEnvironment(m_env, NULL);
        return 4;
    }

    xtb_getEnergy(m_env, m_xtb_res, &energy);
    if (grad != NULL)
        xtb_getGradient(m_env, m_xtb_res, grad);
    return energy;
}

void XTBInterface::clear()
{
    xtb_delResults(&m_xtb_res);
    m_xtb_res = xtb_newResults();
}
