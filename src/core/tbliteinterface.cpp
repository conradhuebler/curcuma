/*
 * < C++ XTB and tblite Interface >
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

#ifdef USE_XTB
#include "external/tblite/include/tblite.h"
#include "external/xtb/include/xtb.h"

#endif

#include "src/core/global.h"
#include "src/core/molecule.h"

#include <iostream>
#include <math.h>
#include <stdio.h>

#include "tbliteinterface.h"

TBLiteInterface::TBLiteInterface()
{
#ifdef USE_XTB
    // tblite stuff
    m_error = tblite_new_error();
    m_ctx = tblite_new_context();
    m_tblite_res = tblite_new_result();

    // xtb stuff
    m_env = xtb_newEnvironment();
    m_xtb_calc = xtb_newCalculator();
    m_xtb_res = xtb_newResults();
#endif
}

TBLiteInterface::~TBLiteInterface()
{
#ifdef USE_XTB
    // tblite stuff
    delete m_error;
    delete m_ctx;
    delete m_tblite_res;
    delete m_tblite_mol;
    delete m_tblite_calc;

    // xtb stuff
    xtb_delResults(&m_xtb_res);
    xtb_delCalculator(&m_xtb_calc);
    xtb_delMolecule(&m_xtb_mol);
    xtb_delEnvironment(&m_env);
#endif
}

bool TBLiteInterface::InitialiseMolecule(const Molecule& molecule)
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

bool TBLiteInterface::InitialiseMolecule(const Molecule* molecule)
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

bool TBLiteInterface::InitialiseMolecule(const int* attyp, const double* coord, const int natoms, const double charge, const int spin)
{
    if (m_initialised)
        UpdateMolecule(coord);

#ifdef USE_XTB
    m_tblite_mol = tblite_new_structure(m_error, natoms, attyp, coord, &charge, &spin, NULL, NULL);
    m_xtb_mol = xtb_newMolecule(m_env, &natoms, attyp, coord, &charge, &spin, NULL, NULL);

    m_initialised = true;
    return true;
#else
    return false;
#endif
}

bool TBLiteInterface::UpdateMolecule(const Molecule& molecule)
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

bool TBLiteInterface::UpdateMolecule(const double* coord)
{
#ifdef USE_XTB
    tblite_update_structure_geometry(m_error, m_tblite_mol, coord, NULL);
    xtb_updateMolecule(m_env, m_xtb_mol, coord, NULL);
    return true;
#else
    return false;
#endif
}

double TBLiteInterface::GFNCalculation(int parameter, double* grad)
{
    double energy = 0;
#ifdef USE_XTB
    if (parameter == 66) {
        xtb_setVerbosity(m_env, XTB_VERBOSITY_MUTED);

        xtb_loadGFNFF(m_env, m_xtb_mol, m_xtb_calc, NULL);
        xtb_singlepoint(m_env, m_xtb_mol, m_xtb_calc, m_xtb_res);
        if (xtb_checkEnvironment(m_env)) {
            xtb_showEnvironment(m_env, NULL);
            return 4;
        }

        xtb_getEnergy(m_env, m_xtb_res, &energy);
        if (grad != NULL)
            xtb_getGradient(m_env, m_xtb_res, grad);
    } else {
        tblite_set_context_verbosity(m_ctx, 0);
        if (parameter == 0) {
            m_tblite_calc = tblite_new_ipea1_calculator(m_ctx, m_tblite_mol);
        } else if (parameter == 1) {
            m_tblite_calc = tblite_new_gfn1_calculator(m_ctx, m_tblite_mol);
        } else if (parameter == 2) {
            m_tblite_calc = tblite_new_gfn2_calculator(m_ctx, m_tblite_mol);
        }

        tblite_get_singlepoint(m_ctx, m_tblite_mol, m_tblite_calc, m_tblite_res);
        tblite_get_result_energy(m_error, m_tblite_res, &energy);

        if (grad != NULL)
            tblite_get_result_gradient(m_error, m_tblite_res, grad);
    }
#else
    throw("XTB is not included, sorry for that");
#endif
    return energy;
}

void TBLiteInterface::clear()
{
#ifdef USE_XTB
    xtb_delResults(&m_xtb_res);
    m_xtb_res = xtb_newResults();
#endif
}
