/*
 * < C++ XTB Interface >
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

#include "external/xtb/include/xtb.h"

#include "src/core/global.h"
#include "src/core/molecule.h"

#include <iostream>
#include <math.h>
#include <stdio.h>

#include "xtbinterface.h"

XTBInterface::XTBInterface()
{
#ifdef USE_XTB
    m_env = xtb_newEnvironment();
    m_calc = xtb_newCalculator();
    m_res = xtb_newResults();
#endif
}

XTBInterface::~XTBInterface()
{
#ifdef USE_XTB
    xtb_delResults(&m_res);
    xtb_delCalculator(&m_calc);
    xtb_delMolecule(&m_mol);
    xtb_delEnvironment(&m_env);
#endif
}

bool XTBInterface::InitialiseMolecule(const Molecule& molecule)
{
    if (m_initialised)
        UpdateMolecule(molecule);
    int const natoms = molecule.AtomCount();
    double const charge = 0.0;

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
    return InitialiseMolecule(attyp, coord, natoms, charge);
}

bool XTBInterface::InitialiseMolecule(const Molecule* molecule)
{
    if (m_initialised)
        UpdateMolecule(molecule);
    int const natoms = molecule->AtomCount();
    double const charge = 0.0;

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
    return InitialiseMolecule(attyp, coord, natoms, charge);
}

bool XTBInterface::InitialiseMolecule(const int* attyp, const double* coord, const int natoms, const double charge)
{
    if (m_initialised)
        UpdateMolecule(coord);

#ifdef USE_XTB
    m_mol = xtb_newMolecule(m_env, &natoms, attyp, coord, &charge, NULL, NULL, NULL);
    if (xtb_checkEnvironment(m_env)) {
        xtb_showEnvironment(m_env, NULL);
        return false;
    }
    m_initialised = true;
    return true;
#else
    return false;
#endif
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
#ifdef USE_XTB
    xtb_updateMolecule(m_env, m_mol, coord, NULL);
    if (xtb_checkEnvironment(m_env)) {
        xtb_showEnvironment(m_env, NULL);
        return false;
    }
    return true;
#else
    return false;
#endif
}

double XTBInterface::GFNCalculation(int parameter, double* grad)
{
    double energy = 0;
#ifdef USE_XTB
    xtb_setVerbosity(m_env, XTB_VERBOSITY_MUTED);
    if (xtb_checkEnvironment(m_env)) {
        xtb_showEnvironment(m_env, NULL);
        return 2;
    }

    if (parameter == 66) {
        xtb_loadGFNFF(m_env, m_mol, m_calc, NULL);
        if (xtb_checkEnvironment(m_env)) {
            xtb_showEnvironment(m_env, NULL);
            return 3;
        }
    } else if (parameter == 0) {
        xtb_loadGFN0xTB(m_env, m_mol, m_calc, NULL);
        if (xtb_checkEnvironment(m_env)) {
            xtb_showEnvironment(m_env, NULL);
            return 3;
        }
    } else if (parameter == 1) {
        xtb_loadGFN1xTB(m_env, m_mol, m_calc, NULL);
        if (xtb_checkEnvironment(m_env)) {
            xtb_showEnvironment(m_env, NULL);
            return 3;
        }
    }

    else if (parameter == 2) {
        xtb_loadGFN2xTB(m_env, m_mol, m_calc, NULL);
        if (xtb_checkEnvironment(m_env)) {
            xtb_showEnvironment(m_env, NULL);
            return 3;
        }
    }

    xtb_singlepoint(m_env, m_mol, m_calc, m_res);
    if (xtb_checkEnvironment(m_env)) {
        xtb_showEnvironment(m_env, NULL);
        return 4;
    }

    xtb_getEnergy(m_env, m_res, &energy);
    if (grad != NULL)
        xtb_getGradient(m_env, m_res, grad);
#else
    throw("XTB is not included, sorry for that");
#endif
    return energy;
}

void XTBInterface::clear()
{
    xtb_delResults(&m_res);
    m_res = xtb_newResults();
}
