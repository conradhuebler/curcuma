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

#include "external/tblite/include/tblite.h"

#include "src/core/global.h"
#include "src/core/molecule.h"

#include <iostream>
#include <math.h>
#include <stdio.h>

#include "tbliteinterface.h"

void example_callback(char* msg, int len, void* udata)
{
    // printf("[callback] %.*s\n", len, msg);
}

TBLiteInterface::TBLiteInterface()
{
#ifdef USE_XTB
    m_error = tblite_new_error();
    //    mol = NULL;
    //    res = NULL;
    m_ctx = tblite_new_context();
    m_res = tblite_new_result();
//    cont = NULL;
//    xtb_setOutput(m_env, "/dev/null");
#endif
}

TBLiteInterface::~TBLiteInterface()
{
#ifdef USE_XTB
    delete m_error;
    delete m_ctx;
    delete m_res;
    delete m_mol;
    delete m_calc;
    /*
    tblite_delete(m_error);
    tblite_delete(m_ctx);
    //tblite_delete(m_cont);
    tblite_delete(m_mol);
    tblite_delete(m_res);
    tblite_delete(m_calc);*/
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
    m_mol = tblite_new_structure(m_error, natoms, attyp, coord, &charge, &spin, NULL, NULL);
    /* if (xtb_checkEnvironment(m_env)) {
         xtb_showEnvironment(m_env, NULL);
         return false;
     }*/
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
    tblite_update_structure_geometry(m_error, m_mol, coord, NULL);
    return true;
    /*  if (xtb_checkEnvironment(m_env)) {
          xtb_showEnvironment(m_env, NULL);
          return false;
      }
      return true;*/
#else
    return false;
#endif
}

double TBLiteInterface::GFNCalculation(int parameter, double* grad)
{
    double energy = 0;
#ifdef USE_XTB
    /*
    auto old_buffer = std::cout.rdbuf(nullptr);

    xtb_setVerbosity(m_env, XTB_VERBOSITY_MUTED);
    if (xtb_checkEnvironment(m_env)) {
        xtb_showEnvironment(m_env, NULL);
        return 2;
    }
*/
    /*
    if (parameter == 66) {

        //char* f = "filename";
        xtb_loadGFNFF(m_env, m_mol, m_calc, f);

        if (xtb_checkEnvironment(m_env)) {
            xtb_showEnvironment(m_env, f);
            return 3;
        }
    } else
*/
    tblite_logger_callback callback = example_callback;
    tblite_set_context_logger(m_ctx, callback, NULL);

    if (parameter == 0) {
        m_calc = tblite_new_ipea1_calculator(m_ctx, m_mol);
    } else if (parameter == 1) {
        m_calc = tblite_new_gfn1_calculator(m_ctx, m_mol);
    } else if (parameter == 2) {
        m_calc = tblite_new_gfn2_calculator(m_ctx, m_mol);
    }
    tblite_get_singlepoint(m_ctx, m_mol, m_calc, m_res);
    tblite_get_result_energy(m_error, m_res, &energy);

    if (grad != NULL)
        tblite_get_result_gradient(m_error, m_res, grad);

        // std::cout.rdbuf(old_buffer);

#else
    throw("XTB is not included, sorry for that");
#endif
    return energy;
}

void TBLiteInterface::clear()
{
#ifdef USE_XTB

#endif
}
