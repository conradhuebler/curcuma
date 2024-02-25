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
    m_xtbsettings = MergeJson(XTBSettings, xtbsettings);
    // std::cout << m_xtbsettings << std::endl;

    m_accuracy = m_xtbsettings["xtb_ac"];
    m_maxiter = m_xtbsettings["xtb_maxiter"];
    m_temp = m_xtbsettings["xtb_temp"];
}

XTBInterface::~XTBInterface()
{
    xtb_delResults(&m_xtb_res);
    xtb_delCalculator(&m_xtb_calc);
    xtb_delMolecule(&m_xtb_mol);
    xtb_delEnvironment(&m_env);
    delete[] m_coord;
    delete[] m_attyp;
}

bool XTBInterface::InitialiseMolecule(const Molecule& molecule)
{
    if (m_initialised)
        UpdateMolecule(molecule);
    m_atomcount = molecule.AtomCount();

    m_attyp = new int[m_atomcount];
    std::vector<int> atoms = molecule.Atoms();
    m_coord = new double[3 * m_atomcount];

    for (int i = 0; i < m_atomcount; ++i) {
        std::pair<int, Position> atom = molecule.Atom(i);
        m_coord[3 * i + 0] = atom.second(0) / au;
        m_coord[3 * i + 1] = atom.second(1) / au;
        m_coord[3 * i + 2] = atom.second(2) / au;
        m_attyp[i] = atoms[i];
    }
    return InitialiseMolecule(m_attyp, m_coord, m_atomcount, molecule.Charge(), molecule.Spin());
}

bool XTBInterface::InitialiseMolecule(const Molecule* molecule)
{
    if (m_initialised)
        UpdateMolecule(molecule);
    m_atomcount = molecule->AtomCount();

    m_attyp = new int[m_atomcount];
    std::vector<int> atoms = molecule->Atoms();
    m_coord = new double[3 * m_atomcount];

    for (int i = 0; i < m_atomcount; ++i) {
        std::pair<int, Position> atom = molecule->Atom(i);
        m_coord[3 * i + 0] = atom.second(0) / au;
        m_coord[3 * i + 1] = atom.second(1) / au;
        m_coord[3 * i + 2] = atom.second(2) / au;
        m_attyp[i] = atoms[i];
    }
    return InitialiseMolecule(m_attyp, m_coord, m_atomcount, molecule->Charge(), molecule->Spin());
}

bool XTBInterface::InitialiseMolecule(const int* attyp, const double* coord, const int natoms, const double charge, const int spin)
{
    if (m_initialised)
        UpdateMolecule(coord);
    m_atomcount = natoms;

    m_xtb_mol = xtb_newMolecule(m_env, &natoms, attyp, coord, &charge, &spin, NULL, NULL);

    m_initialised = true;
    return true;
}

bool XTBInterface::UpdateMolecule(const Molecule& molecule)
{
    for (int i = 0; i < m_atomcount; ++i) {
        std::pair<int, Position> atom = molecule.Atom(i);
        m_coord[3 * i + 0] = atom.second(0) / au;
        m_coord[3 * i + 1] = atom.second(1) / au;
        m_coord[3 * i + 2] = atom.second(2) / au;
    }
    return UpdateMolecule(m_coord);
}

bool XTBInterface::UpdateMolecule(const double* coord)
{
    xtb_updateMolecule(m_env, m_xtb_mol, coord, NULL);
    return true;
}

double XTBInterface::GFNCalculation(int parameter, double* grad)
{
    double energy = 0;
    if (!m_setup) {
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
        xtb_setAccuracy(m_env, m_xtb_calc, m_accuracy);
        xtb_setMaxIter(m_env, m_xtb_calc, m_maxiter);
        xtb_setElectronicTemp(m_env, m_xtb_calc, m_temp);
        m_setup = true;
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

std::vector<double> XTBInterface::Charges() const
{
    std::vector<double> charges(m_atomcount);
    double* c = new double[m_atomcount];
    xtb_getCharges(m_env, m_xtb_res, c);
    for (int i = 0; i < m_atomcount; ++i)
        charges[i] = c[i];
    delete[] c;
    return charges;
}

std::vector<double> XTBInterface::Dipole() const
{
    std::vector<double> dipole(3);
    double* c = new double[3];
    xtb_getDipole(m_env, m_xtb_res, c);
    for (int i = 0; i < 3; ++i)
        dipole[i] = c[i];
    delete[] c;
    return dipole;
}

std::vector<std::vector<double>> XTBInterface::BondOrders() const
{
    std::vector<std::vector<double>> bond_orders(m_atomcount);
    double* bonds = new double[m_atomcount * m_atomcount];
    xtb_getBondOrders(m_env, m_xtb_res, bonds);
    for (int i = 0; i < m_atomcount; ++i) {
        std::vector<double> b(m_atomcount);
        for (int j = 0; j < m_atomcount; ++j)
            b[j] = bonds[i * j];
        bond_orders[i] = b;
    }
    delete[] bonds;
    return bond_orders;
}

void XTBInterface::clear()
{
    xtb_delResults(&m_xtb_res);
    m_xtb_res = xtb_newResults();
}
