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

#ifndef tblite_delete
#include "tblite.h"
// #include "tblite/container.h"
#endif

#include "src/core/global.h"
#include "src/tools/general.h"

#include "src/core/molecule.h"

#include <iostream>
#include <math.h>
#include <stdio.h>

#include "tbliteinterface.h"

TBLiteInterface::TBLiteInterface(const json& tblitesettings)
{
    m_tblitesettings = MergeJson(TBLiteSettings, tblitesettings);

    m_acc = m_tblitesettings["tb_acc"];
    m_maxiter = m_tblitesettings["tb_max_iter"];
    m_damping = m_tblitesettings["tb_damping"];
    m_temp = m_tblitesettings["tb_temp"];
    m_verbose = m_tblitesettings["tb_verbose"];
    std::string guess = m_tblitesettings["tb_guess"];

    m_cpcm_eps = m_tblitesettings["cpcm_eps"];
    m_alpb_eps = m_tblitesettings["alpb_eps"];

    std::string tmp = m_tblitesettings["cpcm_solv"];
    m_cpcm = tmp.compare("none") != 0;
    m_cpcm_solv = new char[tmp.length() + 1];
    strcpy(m_cpcm_solv, tmp.c_str());

    tmp = m_tblitesettings["alpb_solv"];
    m_alpb = tmp.compare("none") != 0;

    m_alpb_solv = new char[tmp.length() + 1];
    strcpy(m_alpb_solv, tmp.c_str());

    if (guess.compare("SAD") == 0)
        m_guess = 0;
    else if (guess.compare("EEQ") == 0)
        m_guess = 1;
    else
        std::cout << "Dont know, ignoring ..." << std::endl;

    m_error = tblite_new_error();
    m_ctx = tblite_new_context();
    m_tblite_res = tblite_new_result();

    tblite_set_context_verbosity(m_ctx, m_verbose);
}

TBLiteInterface::~TBLiteInterface()
{
    delete m_error;
    delete m_ctx;
    delete m_tblite_res;
    delete m_tblite_mol;
    delete m_tblite_calc;
    delete m_tb_cont;

    delete[] m_coord;
    delete[] m_attyp;
}

bool TBLiteInterface::InitialiseMolecule(const Molecule& molecule)
{
    if (m_initialised)
        UpdateMolecule(molecule);
    m_atomcount = molecule.AtomCount();
    m_attyp = new int[m_atomcount];
    m_coord = new double[3 * m_atomcount];

    std::vector<int> atoms = molecule.Atoms();

    for (int i = 0; i < m_atomcount; ++i) {
        std::pair<int, Position> atom = molecule.Atom(i);
        m_coord[3 * i + 0] = atom.second(0) / au;
        m_coord[3 * i + 1] = atom.second(1) / au;
        m_coord[3 * i + 2] = atom.second(2) / au;
        m_attyp[i] = atoms[i];
    }
    bool init = InitialiseMolecule(m_attyp, m_coord, m_atomcount, molecule.Charge(), molecule.Spin());

    return init;
}

bool TBLiteInterface::InitialiseMolecule(const Molecule* molecule)
{
    if (m_initialised)
        UpdateMolecule(molecule);
    m_atomcount = molecule->AtomCount();

    std::vector<int> atoms = molecule->Atoms();
    m_attyp = new int[m_atomcount];
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

bool TBLiteInterface::InitialiseMolecule(const int* attyp, const double* coord, const int natoms, const double charge, const int spin)
{
    if (m_initialised)
        UpdateMolecule(coord);

    m_tblite_mol = tblite_new_structure(m_error, natoms, attyp, coord, &charge, &spin, NULL, NULL);

    m_initialised = true;
    return true;
}

bool TBLiteInterface::UpdateMolecule(const Molecule& molecule)
{
    m_atomcount = molecule.AtomCount();
    for (int i = 0; i < m_atomcount; ++i) {
        std::pair<int, Position> atom = molecule.Atom(i);
        m_coord[3 * i + 0] = atom.second(0) / au;
        m_coord[3 * i + 1] = atom.second(1) / au;
        m_coord[3 * i + 2] = atom.second(2) / au;
    }
    return UpdateMolecule(m_coord);
}

bool TBLiteInterface::UpdateMolecule(const double* coord)
{
    tblite_update_structure_geometry(m_error, m_tblite_mol, coord, NULL);
    return true;
}

double TBLiteInterface::GFNCalculation(int parameter, double* grad)
{
    double energy = 0;
    if (parameter == 0) {
        m_tblite_calc = tblite_new_ipea1_calculator(m_ctx, m_tblite_mol);
    } else if (parameter == 1) {
        m_tblite_calc = tblite_new_gfn1_calculator(m_ctx, m_tblite_mol);
    } else if (parameter == 2) {
        m_tblite_calc = tblite_new_gfn2_calculator(m_ctx, m_tblite_mol);
    }
    if (m_guess == 0)
        tblite_set_calculator_guess(m_ctx, m_tblite_calc, TBLITE_GUESS_SAD);
    else
        tblite_set_calculator_guess(m_ctx, m_tblite_calc, TBLITE_GUESS_EEQ);

    int count = m_cpcm + m_cpcm_eps != -1 + m_alpb + m_alpb_eps != -1;
    if (count == 1) {
        if (m_cpcm && m_cpcm_eps == -1) {
            m_tb_cont = tblite_new_cpcm_solvation_solvent(m_ctx, m_tblite_mol, m_tblite_calc, m_cpcm_solv);
            if (tblite_check_context(m_ctx)) {
                std::cout << "Error during CPCM calculation, ... Sorry for that" << std::endl;
                return 0;
            }
            tblite_calculator_push_back(m_ctx, m_tblite_calc, &m_tb_cont);
        } else if (!m_cpcm && m_cpcm_eps != -1) {
            m_tb_cont = tblite_new_cpcm_solvation_epsilon(m_ctx, m_tblite_mol, m_tblite_calc, m_cpcm_eps);
            if (tblite_check_context(m_ctx)) {
                std::cout << "Error during CPCM calculation, ... Sorry for that" << std::endl;
                return 0;
            }
            tblite_calculator_push_back(m_ctx, m_tblite_calc, &m_tb_cont);
        }

        if (m_alpb && m_alpb_eps == -1) {
            m_tb_cont = tblite_new_alpb_solvation_solvent(m_ctx, m_tblite_mol, m_tblite_calc, m_alpb_solv);
            if (tblite_check_context(m_ctx)) {
                std::cout << "Error during ALPB calculation, ... Sorry for that" << std::endl;
                return 0;
            }
            tblite_calculator_push_back(m_ctx, m_tblite_calc, &m_tb_cont);
        } else if (!m_alpb && m_alpb_eps != -1) {
            m_tb_cont = tblite_new_alpb_solvation_epsilon(m_ctx, m_tblite_mol, m_tblite_calc, m_alpb_eps);
            if (tblite_check_context(m_ctx)) {
                std::cout << "Error during ALPB calculation, ... Sorry for that" << std::endl;
                return 0;
            }
            tblite_calculator_push_back(m_ctx, m_tblite_calc, &m_tb_cont);
        }
    } else {
        std::cout << count << std::endl
                  << std::endl
                  << "If three witches had three watches, which witch would watch which watch?\n Ignoring epsilon and solvent or two solvation models given simultaneously" << std::endl
                  << std::endl
                  << std::endl;
    }
    tblite_set_calculator_accuracy(m_ctx, m_tblite_calc, m_acc);
    tblite_set_calculator_max_iter(m_ctx, m_tblite_calc, m_maxiter);
    tblite_set_calculator_mixer_damping(m_ctx, m_tblite_calc, m_damping);
    tblite_set_calculator_temperature(m_ctx, m_tblite_calc, m_temp);
    tblite_set_calculator_save_integrals(m_ctx, m_tblite_calc, 0);
    tblite_get_singlepoint(m_ctx, m_tblite_mol, m_tblite_calc, m_tblite_res);
    tblite_get_result_energy(m_error, m_tblite_res, &energy);

    if (grad != NULL)
        tblite_get_result_gradient(m_error, m_tblite_res, grad);

    return energy;
}

void TBLiteInterface::clear()
{
}

std::vector<double> TBLiteInterface::Charges() const
{
    std::vector<double> charges(m_atomcount);
    double* c = new double[m_atomcount];
    tblite_get_result_charges(m_error, m_tblite_res, c);
    for (int i = 0; i < m_atomcount; ++i)
        charges[i] = c[i];
    delete[] c;
    return charges;
}

std::vector<double> TBLiteInterface::Dipole() const
{
    std::vector<double> dipole(3);
    double* c = new double[3];
    tblite_get_result_dipole(m_error, m_tblite_res, c);
    for (int i = 0; i < 3; ++i)
        dipole[i] = c[i];
    delete[] c;
    return dipole;
}

std::vector<std::vector<double>> TBLiteInterface::BondOrders() const
{
    std::vector<std::vector<double>> bond_orders(m_atomcount);
    double* bonds = new double[m_atomcount * m_atomcount];
    tblite_get_result_bond_orders(m_error, m_tblite_res, bonds);
    for (int i = 0; i < m_atomcount; ++i) {
        std::vector<double> b(m_atomcount);
        for (int j = 0; j < m_atomcount; ++j)
            b[j] = bonds[i * j];
        bond_orders[i] = b;
    }
    delete[] bonds;
    return bond_orders;
}
