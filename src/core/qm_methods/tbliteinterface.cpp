/*
 * < C++ XTB and tblite Interface >
 * Copyright (C) 2020 - 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
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
// #include "tblite.h"
#include "tblite/container.h"
#include "tblite/context.h"
#include "tblite/error.h"
#include "tblite/result.h"
#include "tblite/solvation.h"
#endif

#include "src/core/global.h"

#include <iostream>
#include <math.h>
#include <stdio.h>

#include "tbliteinterface.h"

TBLiteInterface::TBLiteInterface(const json& tblitesettings)
{
    m_tblitesettings = MergeJson(TBLiteSettings, tblitesettings);
    // std::cout << "Initialising TBLite Interface" << std::endl;
    // std::cout << m_tblitesettings << std::endl;
    m_acc = m_tblitesettings["tb_acc"];
    m_SCFmaxiter = m_tblitesettings["SCFmaxiter"];
    m_damping = m_tblitesettings["tb_damping"];

    // convert kelvin to atomic units
    m_Tele = m_tblitesettings["Tele"];
    m_Tele /= 315775.326864009;
    m_verbose = m_tblitesettings["verbose"];
    m_spin = m_tblitesettings["spin"];
    std::string guess = m_tblitesettings["tb_guess"];
    m_solvent_eps = m_tblitesettings["solvent_eps"];
    m_solvent_model = m_tblitesettings["solvent_model"];

    std::string tmp = m_tblitesettings["solvent"];
    m_solvent = new char[tmp.length() + 1];
    strcpy(m_solvent, tmp.c_str());
    m_solvent_gb_version = m_tblitesettings["solvent_gb_version"];
    m_solvent_gb_kernel = m_tblitesettings["solvent_gb_kernel"];
    m_solvent_alpb_version = m_tblitesettings["solvent_alpb_version"];
    m_solvent_alpb_reference = m_tblitesettings["solvent_alpb_reference"];

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
    if (m_verbose > 0)
        std::cout << "Initialising tblite with accuracy " << m_acc << " and SCFmaxiter " << m_SCFmaxiter << std::endl;
}

TBLiteInterface::~TBLiteInterface()
{
    tblite_delete_error(&m_error);
    tblite_delete_context(&m_ctx);
    tblite_delete_result(&m_tblite_res);
    tblite_delete_structure(&m_tblite_mol);
    tblite_delete_calculator(&m_tblite_calc);
    tblite_delete_container(&m_tb_cont);

    delete[] m_coord;
    delete[] m_attyp;
}

bool TBLiteInterface::InitialiseMolecule(const Mol& mol)
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
    return InitialiseMolecule(m_attyp, m_coord, m_atomcount, mol.m_charge, mol.m_spin);
}

bool TBLiteInterface::InitialiseMolecule(const Mol* mol)
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

bool TBLiteInterface::InitialiseMolecule(const int* attyp, const double* coord, const int natoms, const double charge, const int spin)
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
    InitialiseMolecule();

    return m_initialised;
}

bool TBLiteInterface::InitialiseMolecule()
{
    if (m_initialised) {
        return UpdateMolecule(m_coord);
    }

    double charge = m_charge;
    int spin = m_spin;
    m_tblite_mol = tblite_new_structure(m_error, m_atomcount, m_attyp, m_coord, &charge, &spin, NULL, NULL);
    if (tblite_check_error(m_error)) {
        tbliteError();
        return false;
    }
    m_gradient = Matrix::Zero(m_atomcount, 3);
    m_initialised = true;
    return m_initialised;
}

bool TBLiteInterface::UpdateMolecule(const Mol& mol)
{
    for (int i = 0; i < m_atomcount; ++i) {
        m_coord[3 * i + 0] = mol.m_geometry(i, 0) / au;
        m_coord[3 * i + 1] = mol.m_geometry(i, 1) / au;
        m_coord[3 * i + 2] = mol.m_geometry(i, 2) / au;
    }
    return UpdateMolecule();
}

bool TBLiteInterface::UpdateMolecule(const Mol* mol)
{
    for (int i = 0; i < m_atomcount; ++i) {
        m_coord[3 * i + 0] = mol->m_geometry(i, 0) / au;
        m_coord[3 * i + 1] = mol->m_geometry(i, 1) / au;
        m_coord[3 * i + 2] = mol->m_geometry(i, 2) / au;
    }
    return UpdateMolecule();
}

bool TBLiteInterface::UpdateMolecule(const double* coord)
{
    for (int i = 0; i < m_atomcount; ++i) {
        m_coord[3 * i + 0] = coord[3 * i + 0] / au;
        m_coord[3 * i + 1] = coord[3 * i + 1] / au;
        m_coord[3 * i + 2] = coord[3 * i + 2] / au;
    }
    return UpdateMolecule();
}

bool TBLiteInterface::UpdateMolecule(const Matrix& geometry)
{
    for (int i = 0; i < m_atomcount; ++i) {
        m_coord[3 * i + 0] = geometry(i, 0) / au;
        m_coord[3 * i + 1] = geometry(i, 1) / au;
        m_coord[3 * i + 2] = geometry(i, 2) / au;
    }
    return UpdateMolecule();
}

bool TBLiteInterface::UpdateMolecule()
{
    tblite_update_structure_geometry(m_error, m_tblite_mol, m_coord, NULL);
    return true;
}

void TBLiteInterface::clear()
{
    tblite_delete_structure(&m_tblite_mol);
    m_initialised = false;
}

void TBLiteInterface::ApplySolvation()
{
    if (m_verbose > 0) {
        std::cout << "Applying Solvation Model" << std::endl;
        std::cout << "Solvent Model: " << m_solvent_model << std::endl;
    }
    if (m_solvent_model == 0)
        return;
    else if (m_solvent_model == 1 && m_solvent_eps > -1) {
        if (m_verbose > 0) {
            std::cout << "Using CPCM with epsilon " << m_solvent_eps << std::endl;
        }
        m_tb_cont = tblite_new_cpcm_solvation_epsilon(m_error, m_tblite_mol, m_solvent_eps);
        if (tblite_check_context(m_ctx)) {
            tbliteError();
            tbliteContextError();
            std::cout << "Error during CPCM calculation, ...  Retry with increased verbosity" << std::endl;
            tblite_set_context_verbosity(m_ctx, m_verbose + 1);
            tblite_delete_container(&m_tb_cont);
            m_tb_cont = tblite_new_cpcm_solvation_epsilon(m_error, m_tblite_mol, m_solvent_eps);
        }
        tblite_calculator_push_back(m_ctx, m_tblite_calc, &m_tb_cont);
        //  tbliteError();
        //   tbliteContextError();
    } else if (m_solvent_model == 2 && m_solvent_eps > -1) {
        if (m_verbose > 0) {
            std::cout << "Using GB with epsilon " << m_solvent_eps << " and born version " << m_solvent_gb_version << " and born kernel " << m_solvent_gb_kernel << std::endl;
            std::cout << "GB doesnt work for now" << std::endl;
        }
        exit(1);
        m_tb_cont = tblite_new_gb_solvation_epsilon(m_error, m_tblite_mol, m_solvent_eps, m_solvent_gb_version, m_solvent_gb_kernel);
        if (tblite_check_context(m_ctx)) {
            tbliteError();
            tbliteContextError();
            std::cout << "Error during ALPB calculation, ...  Retry with increased verbosity" << std::endl;
            tblite_set_context_verbosity(m_ctx, m_verbose + 1);
            tblite_delete_container(&m_tb_cont);
            m_tb_cont = tblite_new_gb_solvation_epsilon(m_error, m_tblite_mol, m_solvent_eps, m_solvent_gb_version, m_solvent_gb_kernel);
        }
        tblite_calculator_push_back(m_ctx, m_tblite_calc, &m_tb_cont);
        //  tbliteError();
        //  tbliteContextError();
    } else if (m_solvent_model == 3) {
        if (m_verbose > 0) {
            std::cout << "Using ALPB with solvent " << m_solvent << " and version " << m_solvent_alpb_version << " and reference " << m_solvent_alpb_reference << std::endl;
        }
        m_tb_cont = tblite_new_alpb_solvation_solvent(m_error, m_tblite_mol, m_solvent, m_solvent_alpb_version, m_solvent_alpb_reference);
        if (tblite_check_context(m_ctx)) {
            tbliteError();
            tbliteContextError();
            std::cout << "Error during ALPB calculation, ... Retry with increased verbosity" << std::endl;

            tblite_set_context_verbosity(m_ctx, m_verbose + 1);
            tblite_delete_container(&m_tb_cont);
            m_tb_cont = tblite_new_alpb_solvation_solvent(m_error, m_tblite_mol, m_solvent, m_solvent_alpb_version, m_solvent_alpb_reference);
        }
        tblite_calculator_push_back(m_ctx, m_tblite_calc, &m_tb_cont);
        //  tbliteError();
        //  tbliteContextError();
    }
}

double TBLiteInterface::Calculation(bool gradient, bool verbose)
{
    if (verbose)
        tblite_set_context_verbosity(m_ctx, 3);
    else
        tblite_set_context_verbosity(m_ctx, m_verbose);

    double energy = 0;
    int count = 0;
    if (!m_calculator) {
        if (m_method_switch == 0) {
            m_tblite_calc = tblite_new_ipea1_calculator(m_ctx, m_tblite_mol);
        } else if (m_method_switch == 1) {
            m_tblite_calc = tblite_new_gfn1_calculator(m_ctx, m_tblite_mol);
        } else if (m_method_switch == 2) {
            m_tblite_calc = tblite_new_gfn2_calculator(m_ctx, m_tblite_mol);
        }
        if (m_guess == 0)
            tblite_set_calculator_guess(m_ctx, m_tblite_calc, TBLITE_GUESS_SAD);
        else
            tblite_set_calculator_guess(m_ctx, m_tblite_calc, TBLITE_GUESS_EEQ);
        m_calculator = true;
    }

    ApplySolvation();

    tblite_set_calculator_accuracy(m_ctx, m_tblite_calc, m_acc);
    tblite_set_calculator_max_iter(m_ctx, m_tblite_calc, m_SCFmaxiter);
    tblite_set_calculator_mixer_damping(m_ctx, m_tblite_calc, m_damping);
    tblite_set_calculator_temperature(m_ctx, m_tblite_calc, m_Tele);
    tblite_set_calculator_save_integrals(m_ctx, m_tblite_calc, 0);
    tblite_get_singlepoint(m_ctx, m_tblite_mol, m_tblite_calc, m_tblite_res);
    tblite_get_result_energy(m_error, m_tblite_res, &energy);
    if (tblite_check_context(m_ctx)) {
        m_error_count++;
        std::cerr << "Error during SCF Calculation, reset wavefunction ..." << std::endl;
        tbliteContextError();
        tbliteError();

        tblite_set_context_verbosity(m_ctx, m_verbose + 1);
        if (count == 1) {
            tblite_delete_container(&m_tb_cont);
        }
        tblite_delete_error(&m_error);
        tblite_delete_context(&m_ctx);
        tblite_delete_result(&m_tblite_res);
        // tblite_delete_structure(  &m_tblite_mol );
        tblite_delete_calculator(&m_tblite_calc);
        tblite_delete_container(&m_tb_cont);

        m_error = tblite_new_error();
        m_ctx = tblite_new_context();
        m_tblite_res = tblite_new_result();

        tblite_set_context_verbosity(m_ctx, m_verbose);

        if (m_method_switch == 0) {
            m_tblite_calc = tblite_new_ipea1_calculator(m_ctx, m_tblite_mol);
        } else if (m_method_switch == 1) {
            m_tblite_calc = tblite_new_gfn1_calculator(m_ctx, m_tblite_mol);
        } else if (m_method_switch == 2) {
            m_tblite_calc = tblite_new_gfn2_calculator(m_ctx, m_tblite_mol);
        }
        if (m_guess == 0)
            tblite_set_calculator_guess(m_ctx, m_tblite_calc, TBLITE_GUESS_SAD);
        else
            tblite_set_calculator_guess(m_ctx, m_tblite_calc, TBLITE_GUESS_EEQ);
        ApplySolvation();
        tblite_get_singlepoint(m_ctx, m_tblite_mol, m_tblite_calc, m_tblite_res);
        tblite_set_context_verbosity(m_ctx, m_verbose);
    }
    // double* virial = 0;
    // tblite_get_result_virial(m_error, m_tblite_res, &virial);
    // std::cout << virial << std::endl;
    if (gradient) {
        tblite_get_result_gradient(m_error, m_tblite_res, m_gradient.data());
        m_gradient *= au;
        if (tblite_check_context(m_ctx)) {
            tbliteContextError();
            //     return 0;
        }
    }
    if (count == 1)
        tblite_delete_container(&m_tb_cont);
    return energy;
}

void TBLiteInterface::setMethod(const std::string& method)
{
    QMInterface::setMethod(method);

    if (m_method.compare("ipea1") == 0)
        m_method_switch = 0;
    else if (m_method.compare("gfn1") == 0)
        m_method_switch = 1;
    else if (m_method.compare("gfn2") == 0)
        m_method_switch = 2;
}

Vector TBLiteInterface::Charges() const
{
    Eigen::VectorXd charges(m_atomcount);
    tblite_get_result_charges(m_error, m_tblite_res, charges.data());
    return charges;
}

Vector TBLiteInterface::Dipole() const
{
    Vector dipole(3);
    tblite_get_result_dipole(m_error, m_tblite_res, dipole.data());
    return dipole;
}

Vector TBLiteInterface::BondOrders() const
{
    Eigen::VectorXd bonds(m_atomcount * m_atomcount);
    tblite_get_result_bond_orders(m_error, m_tblite_res, bonds.data());

    return bonds;
}

Vector TBLiteInterface::OrbitalEnergies() const
{
    int num_orbitals;
    tblite_get_result_number_of_orbitals(m_error, m_tblite_res, &num_orbitals);

    Vector orbital_energies(num_orbitals);
    tblite_get_result_orbital_energies(m_error, m_tblite_res, orbital_energies.data());

    return orbital_energies;
}

Vector TBLiteInterface::OrbitalOccupations() const
{
    int num_orbitals;
    tblite_get_result_number_of_orbitals(m_error, m_tblite_res, &num_orbitals);

    Vector orbital_occupations(num_orbitals);
    tblite_get_result_orbital_occupations(m_error, m_tblite_res, orbital_occupations.data());

    return orbital_occupations;
}

void TBLiteInterface::tbliteError()
{
    char message[512];
    tblite_get_error(m_error, message, NULL);
    std::cerr << "[Message] " << message << std::endl;

    // printf("[Message] %s\n", message);
    tblite_clear_error(m_error);
}

void TBLiteInterface::tbliteContextError()
{
    char message[512];
    tblite_get_context_error(m_ctx, message, NULL);
    std::cerr << "[Message] " << message << std::endl;
    // printf("[Message] %s\n", message);
    tblite_clear_error(m_error);
}
