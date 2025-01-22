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
#include "interface/abstract_interface.h"

#include "src/core/global.h"
#include "src/tools/general.h"

#include "external/xtb/include/xtb.h"

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
    m_SCFmaxiter = m_xtbsettings["SCFmaxiter"];
    m_Tele = m_xtbsettings["Tele"];
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

bool XTBInterface::InitialiseMolecule(const Mol& mol)
{
    if (m_initialised)
        UpdateMolecule(mol);
    m_atomcount = mol.m_number_atoms;
    m_attyp = new int[m_atomcount];
    m_coord = new double[3 * m_atomcount];

    std::vector<int> atoms = mol.m_atoms; // molecule.Atoms();

    for (int i = 0; i < m_atomcount; ++i) {
        m_coord[3 * i + 0] = mol.m_geometry(i, 0) / au;
        m_coord[3 * i + 1] = mol.m_geometry(i, 1) / au;
        m_coord[3 * i + 2] = mol.m_geometry(i, 2) / au;
        m_attyp[i] = atoms[i];
    }
    bool init = InitialiseMolecule(m_attyp, m_coord, m_atomcount, mol.m_charge, mol.m_spin);
    return init;
}

bool XTBInterface::InitialiseMolecule(const Mol* mol)
{
    if (m_initialised) {
        return UpdateMolecule(mol);
    }
    m_atomcount = mol->m_number_atoms;

    std::vector<int> atoms = mol->m_atoms;
    m_attyp = new int[m_atomcount];
    m_coord = new double[3 * m_atomcount];

    for (int i = 0; i < m_atomcount; ++i) {
        m_coord[3 * i + 0] = mol->m_geometry(i, 0) / au;
        m_coord[3 * i + 1] = mol->m_geometry(i, 1) / au;
        m_coord[3 * i + 2] = mol->m_geometry(i, 2) / au;
        m_attyp[i] = atoms[i];
    }
    return InitialiseMolecule();
}

bool XTBInterface::InitialiseMolecule(const int* attyp, const double* coord, const int natoms, const double charge, const int spin)
{
    if (m_initialised)
        UpdateMolecule(coord);

    m_atomcount = natoms;
    m_charge = charge;
    m_spin = spin;

    for (int i = 0; i < m_atomcount; ++i) {
        m_coord[3 * i + 0] = coord[3 * i + 0] / au;
        m_coord[3 * i + 1] = coord[3 * i + 1] / au;
        m_coord[3 * i + 2] = coord[3 * i + 2] / au;
        m_attyp[i] = attyp[i];
    }
    return InitialiseMolecule();
}

bool XTBInterface::InitialiseMolecule()
{
    if (m_initialised) {
        return UpdateMolecule(m_coord);
    }
    double charge = m_charge;
    int spin = m_spin;
    m_gradient = Matrix::Zero(m_atomcount, 3);

    m_xtb_mol = xtb_newMolecule(m_env, &m_atomcount, m_attyp, m_coord, &charge, &spin, NULL, NULL);
    /*
        char *error;
        int buffer;
        xtb_getError(m_env, &error, &buffer);
        if (error != NULL) {
            std::cerr << "[Error] " << error << std::endl;
            return false;
        }
    */
    m_initialised = true;
    return m_initialised;
}

bool XTBInterface::UpdateMolecule(const Mol& mol)
{
    for (int i = 0; i < m_atomcount; ++i) {
        m_coord[3 * i + 0] = mol.m_geometry(i, 0) / au;
        m_coord[3 * i + 1] = mol.m_geometry(i, 1) / au;
        m_coord[3 * i + 2] = mol.m_geometry(i, 2) / au;
    }
    return UpdateMolecule();
}

bool XTBInterface::UpdateMolecule(const Mol* mol)
{
    for (int i = 0; i < m_atomcount; ++i) {
        m_coord[3 * i + 0] = mol->m_geometry(i, 0) / au;
        m_coord[3 * i + 1] = mol->m_geometry(i, 1) / au;
        m_coord[3 * i + 2] = mol->m_geometry(i, 2) / au;
    }
    return UpdateMolecule();
}

bool XTBInterface::UpdateMolecule(const double* coord)
{
    for (int i = 0; i < m_atomcount; ++i) {
        m_coord[3 * i + 0] = coord[3 * i + 0] / au;
        m_coord[3 * i + 1] = coord[3 * i + 1] / au;
        m_coord[3 * i + 2] = coord[3 * i + 2] / au;
    }
    return UpdateMolecule();
}

bool XTBInterface::UpdateMolecule(const Matrix& geometry)
{
    for (int i = 0; i < m_atomcount; ++i) {
        m_coord[3 * i + 0] = geometry(i, 0) / au;
        m_coord[3 * i + 1] = geometry(i, 1) / au;
        m_coord[3 * i + 2] = geometry(i, 2) / au;
    }
    return UpdateMolecule();
}
bool XTBInterface::UpdateMolecule()
{
    xtb_updateMolecule(m_env, m_xtb_mol, m_coord, NULL);
    return true;
}

void XTBInterface::setMethod(const std::string& method)
{
    QMInterface::setMethod(method);
    if (method == "xtb-gfn0")
        m_method_switch = 0;
    else if (method == "xtb-gfn1")
        m_method_switch = 1;
    else if (method == "xtb-gfn2")
        m_method_switch = 2;
    else if (method == "gfnff")
        m_method_switch = 66;
}

double XTBInterface::Calculation(bool gradient, bool verbose)
{
    double energy = 0;
    if (!m_setup) {
        xtb_setVerbosity(m_env, XTB_VERBOSITY_MUTED);
        if (m_method_switch == 0) {
            xtb_loadGFN0xTB(m_env, m_xtb_mol, m_xtb_calc, NULL);
        } else if (m_method_switch == 1) {
            xtb_loadGFN1xTB(m_env, m_xtb_mol, m_xtb_calc, NULL);
        } else if (m_method_switch == 2) {
            xtb_loadGFN2xTB(m_env, m_xtb_mol, m_xtb_calc, NULL);
        } else if (m_method_switch == 66) {
            xtb_loadGFNFF(m_env, m_xtb_mol, m_xtb_calc, NULL);
        }
        xtb_setAccuracy(m_env, m_xtb_calc, m_accuracy);
        xtb_setMaxIter(m_env, m_xtb_calc, m_SCFmaxiter);
        xtb_setElectronicTemp(m_env, m_xtb_calc, m_Tele);
        m_setup = true;
    }
    xtb_singlepoint(m_env, m_xtb_mol, m_xtb_calc, m_xtb_res);
    if (xtb_checkEnvironment(m_env)) {
        xtb_showEnvironment(m_env, NULL);
        return 4;
    }

    xtb_getEnergy(m_env, m_xtb_res, &energy);
    if (gradient) {
        xtb_getGradient(m_env, m_xtb_res, m_gradient.data());
        m_gradient *= au;
    }
    return energy;
}

Vector XTBInterface::Charges() const
{
    Vector charges(m_atomcount);
    xtb_getCharges(m_env, m_xtb_res, charges.data());
    return charges;
}

Vector XTBInterface::Dipole() const
{
    Vector dipole(3);
    xtb_getDipole(m_env, m_xtb_res, dipole.data());

    return dipole;
}

Vector XTBInterface::BondOrders() const
{
    Vector bond_orders(m_atomcount);
    xtb_getBondOrders(m_env, m_xtb_res, bond_orders.data());
    return bond_orders;
}

Vector XTBInterface::OrbitalEnergies() const
{
    int num_orbitals;
    xtb_getNao(m_env, m_xtb_res, &num_orbitals);

    Vector orbital_energies(num_orbitals);
    xtb_getOrbitalEigenvalues(m_env, m_xtb_res, orbital_energies.data());

    return orbital_energies;
}

Vector XTBInterface::OrbitalOccupations() const
{
    int num_orbitals;
    xtb_getNao(m_env, m_xtb_res, &num_orbitals);

    Vector orbital_occupations(num_orbitals);
    xtb_getOrbitalOccupations(m_env, m_xtb_res, orbital_occupations.data());

    return orbital_occupations;
}

void XTBInterface::clear()
{
    xtb_delResults(&m_xtb_res);
    m_xtb_res = xtb_newResults();
}
