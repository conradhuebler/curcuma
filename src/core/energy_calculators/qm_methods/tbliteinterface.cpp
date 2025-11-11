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

#ifdef USE_TBLITE
// #include "tblite.h"
#include "tblite/container.h"
#include "tblite/context.h"
#include "tblite/error.h"
#include "tblite/result.h"
#include "tblite/solvation.h"
#endif

#include "src/core/curcuma_logger.h"
#include "src/core/global.h"
#include <fmt/format.h>

#include <iostream>
#include <math.h>
#include <stdio.h>

#include "tbliteinterface.h"

TBLiteInterface::TBLiteInterface(const ConfigManager& config)
    : m_config(config)
{
    // Claude Generated: Type-safe parameter access via ConfigManager (Phase 3B)
    m_acc = m_config.get<int>("accuracy", 1);
    m_SCFmaxiter = m_config.get<int>("max_iterations", 100);
    m_damping = m_config.get<double>("damping", 0.4);

    // Convert Kelvin to atomic units
    m_Tele = m_config.get<double>("electronic_temperature", 300.0);
    m_Tele /= 315775.326864009;
    m_spin = m_config.get<double>("spin", 0.0);
    std::string guess = m_config.get<std::string>("initial_guess", "SAD");
    m_solvent_eps = m_config.get<double>("solvent_epsilon", -1.0);
    m_solvent_model = m_config.get<int>("solvent_model", 0);

    std::string tmp = m_config.get<std::string>("solvent", "none");
    m_solvent = new char[tmp.length() + 1];
    strcpy(m_solvent, tmp.c_str());
    m_solvent_gb_version = m_config.get<int>("solvent_gb_version", 0);
    m_solvent_gb_kernel = m_config.get<int>("solvent_gb_kernel", 1);
    m_solvent_alpb_version = m_config.get<int>("solvent_alpb_version", 12);
    m_solvent_alpb_reference = m_config.get<int>("solvent_alpb_reference", 1);

    if (guess.compare("SAD") == 0)
        m_guess = 0;
    else if (guess.compare("EEQ") == 0)
        m_guess = 1;
    else {
        CurcumaLogger::warn("Unknown guess method, using SAD default");
    }

    m_error = tblite_new_error();
    m_ctx = tblite_new_context();
    m_tblite_res = tblite_new_result();

    // Set TBLite verbosity based on our logging system
    if (CurcumaLogger::get_verbosity() >= 3) {
        tblite_set_context_verbosity(m_ctx, 3);
    } else if (CurcumaLogger::get_verbosity() >= 2) {
        tblite_set_context_verbosity(m_ctx, 1);
    } else {
        tblite_set_context_verbosity(m_ctx, 0);
    }

    // Verbosity Level 1+: TBLite initialization
    if (CurcumaLogger::get_verbosity() >= 1) {
        CurcumaLogger::info("Initializing TBLite quantum chemistry method");
        CurcumaLogger::param("accuracy", fmt::format("{}", m_acc));
        CurcumaLogger::param("SCF_maxiter", fmt::format("{}", m_SCFmaxiter));
        CurcumaLogger::param("temperature", fmt::format("{:.1f} K", m_Tele * 315775.326864009));
        CurcumaLogger::param("damping", fmt::format("{:.3f}", m_damping));
        if (m_solvent_model > 0) {
            CurcumaLogger::param("solvation", "enabled");
        }
    }
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

    // CRITICAL FIX: coord is already in atomic units (Bohr) - no conversion needed!
    // Previous bug: double conversion coord/au -> m_coord/au = coord/au^2
    for (int i = 0; i < m_atomcount; ++i) {
        m_coord[3 * i + 0] = coord[3 * i + 0]; // Already in Bohr from previous /au
        m_coord[3 * i + 1] = coord[3 * i + 1]; // Already in Bohr from previous /au
        m_coord[3 * i + 2] = coord[3 * i + 2]; // Already in Bohr from previous /au
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
    // Level 2+: Solvation setup
    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("Applying solvation model");
        CurcumaLogger::param("solvent_model", m_solvent_model);
    }

    if (m_solvent_model == 0)
        return;
    else if (m_solvent_model == 1 && m_solvent_eps > -1) {
        if (CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::param("solvation_method", "CPCM");
            CurcumaLogger::param("dielectric_constant", fmt::format("{:.3f}", m_solvent_eps));
        }
        m_tb_cont = tblite_new_cpcm_solvation_epsilon(m_error, m_tblite_mol, m_solvent_eps);
        if (tblite_check_context(m_ctx)) {
            tbliteError();
            tbliteContextError();
            CurcumaLogger::warn("Error during CPCM calculation - retrying with increased verbosity");
            tblite_set_context_verbosity(m_ctx, 3);
            tblite_delete_container(&m_tb_cont);
            m_tb_cont = tblite_new_cpcm_solvation_epsilon(m_error, m_tblite_mol, m_solvent_eps);
        }
        tblite_calculator_push_back(m_ctx, m_tblite_calc, &m_tb_cont);

    } else if (m_solvent_model == 2 && m_solvent_eps > -1) {
        if (CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::param("solvation_method", "GB");
            CurcumaLogger::param("dielectric_constant", fmt::format("{:.3f}", m_solvent_eps));
            CurcumaLogger::param("gb_version", m_solvent_gb_version);
            CurcumaLogger::param("gb_kernel", m_solvent_gb_kernel);
        }
        CurcumaLogger::error("GB solvation model not functional");
        exit(1);
        m_tb_cont = tblite_new_gb_solvation_epsilon(m_error, m_tblite_mol, m_solvent_eps, m_solvent_gb_version, m_solvent_gb_kernel);
        if (tblite_check_context(m_ctx)) {
            tbliteError();
            tbliteContextError();
            CurcumaLogger::warn("Error during GB calculation - retrying with increased verbosity");
            tblite_set_context_verbosity(m_ctx, 3);
            tblite_delete_container(&m_tb_cont);
            m_tb_cont = tblite_new_gb_solvation_epsilon(m_error, m_tblite_mol, m_solvent_eps, m_solvent_gb_version, m_solvent_gb_kernel);
        }
        tblite_calculator_push_back(m_ctx, m_tblite_calc, &m_tb_cont);

    } else if (m_solvent_model == 3) {
        if (CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::param("solvation_method", "ALPB");
            CurcumaLogger::param("solvent", std::string(m_solvent));
            CurcumaLogger::param("alpb_version", m_solvent_alpb_version);
            CurcumaLogger::param("alpb_reference", m_solvent_alpb_reference);
        }
        m_tb_cont = tblite_new_alpb_solvation_solvent(m_error, m_tblite_mol, m_solvent, m_solvent_alpb_version, m_solvent_alpb_reference);
        if (tblite_check_context(m_ctx)) {
            tbliteError();
            tbliteContextError();
            CurcumaLogger::warn("Error during ALPB calculation - retrying with increased verbosity");

            tblite_set_context_verbosity(m_ctx, 3);
            tblite_delete_container(&m_tb_cont);
            m_tb_cont = tblite_new_alpb_solvation_solvent(m_error, m_tblite_mol, m_solvent, m_solvent_alpb_version, m_solvent_alpb_reference);
        }
        tblite_calculator_push_back(m_ctx, m_tblite_calc, &m_tb_cont);
    }
}

double TBLiteInterface::Calculation(bool gradient)
{
    // Control TBLite verbosity based on our logging system
    if (CurcumaLogger::get_verbosity() >= 3) {
        tblite_set_context_verbosity(m_ctx, 3);
    } else if (CurcumaLogger::get_verbosity() >= 2) {
        tblite_set_context_verbosity(m_ctx, 1);
    } else {
        tblite_set_context_verbosity(m_ctx, 0);
    }

    double energy = 0;
    int count = 0;

    // Setup calculator on first call
    if (!m_calculator) {
        // Level 2+: Method setup info
        if (CurcumaLogger::get_verbosity() >= 2) {
            std::string method_name = "unknown";
            if (m_method_switch == 0)
                method_name = "iPEA1";
            else if (m_method_switch == 1)
                method_name = "GFN1-xTB";
            else if (m_method_switch == 2)
                method_name = "GFN2-xTB";

            CurcumaLogger::info("Starting TBLite calculation");
            CurcumaLogger::param("method", method_name);
            CurcumaLogger::param("guess", (m_guess == 0) ? "SAD" : "EEQ");
            if (gradient) {
                CurcumaLogger::param("gradient", "analytical");
            }
        }

        // Create appropriate calculator
        if (m_method_switch == 0) {
            m_tblite_calc = tblite_new_ipea1_calculator(m_ctx, m_tblite_mol);
        } else if (m_method_switch == 1) {
            m_tblite_calc = tblite_new_gfn1_calculator(m_ctx, m_tblite_mol);
        } else if (m_method_switch == 2) {
            m_tblite_calc = tblite_new_gfn2_calculator(m_ctx, m_tblite_mol);
        }

        // CRITICAL: Configure calculator parameters for convergence
        // Apply method-specific optimizations for GFN2 convergence issues
        double effective_damping = m_damping;
        int effective_maxiter = m_SCFmaxiter;
        double effective_accuracy = m_acc;

        if (m_method_switch == 2) { // GFN2-xTB needs stronger settings
            effective_damping = std::max(0.6, m_damping); // Stronger damping for GFN2
            effective_maxiter = std::max(150, m_SCFmaxiter); // More iterations for GFN2
            effective_accuracy = std::min(0.5, (double)m_acc); // Higher accuracy for GFN2

            if (CurcumaLogger::get_verbosity() >= 2) {
                CurcumaLogger::info("Applying GFN2-specific convergence optimizations");
                CurcumaLogger::param("enhanced_damping", std::to_string(effective_damping));
                CurcumaLogger::param("enhanced_maxiter", std::to_string(effective_maxiter));
                CurcumaLogger::param("enhanced_accuracy", std::to_string(effective_accuracy));
            }
        }

        tblite_set_calculator_accuracy(m_ctx, m_tblite_calc, effective_accuracy);
        tblite_set_calculator_max_iter(m_ctx, m_tblite_calc, effective_maxiter);
        tblite_set_calculator_mixer_damping(m_ctx, m_tblite_calc, effective_damping);
        tblite_set_calculator_temperature(m_ctx, m_tblite_calc, m_Tele);

        if (m_guess == 0)
            tblite_set_calculator_guess(m_ctx, m_tblite_calc, TBLITE_GUESS_SAD);
        else
            tblite_set_calculator_guess(m_ctx, m_tblite_calc, TBLITE_GUESS_EEQ);

        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info("TBLite Calculator configured:");
            CurcumaLogger::param("accuracy", std::to_string(m_acc));
            CurcumaLogger::param("max_iterations", std::to_string(m_SCFmaxiter));
            CurcumaLogger::param("damping", std::to_string(m_damping));
            CurcumaLogger::param("temperature", std::to_string(m_Tele));
            CurcumaLogger::param("guess", (m_guess == 0) ? "SAD" : "EEQ");
        }
        m_calculator = true;
    }

    ApplySolvation();

    // Parameters already set during calculator creation - no need to repeat
    tblite_set_calculator_save_integrals(m_ctx, m_tblite_calc, 0);
    // Perform calculation
    tblite_get_singlepoint(m_ctx, m_tblite_mol, m_tblite_calc, m_tblite_res);
    tblite_get_result_energy(m_error, m_tblite_res, &energy);

    // Handle calculation errors
    if (tblite_check_context(m_ctx)) {
        m_error_count++;
        CurcumaLogger::warn("Error during TBLite SCF calculation - resetting wavefunction");
        tbliteContextError();
        tbliteError();

        tblite_set_context_verbosity(m_ctx, 3);
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

        // CRITICAL: Re-apply TBLite verbosity after context recreation
        if (CurcumaLogger::get_verbosity() >= 3) {
            tblite_set_context_verbosity(m_ctx, 3);
        } else if (CurcumaLogger::get_verbosity() >= 2) {
            tblite_set_context_verbosity(m_ctx, 1);
        } else {
            tblite_set_context_verbosity(m_ctx, 0);
        }

        if (m_method_switch == 0) {
            m_tblite_calc = tblite_new_ipea1_calculator(m_ctx, m_tblite_mol);
        } else if (m_method_switch == 1) {
            m_tblite_calc = tblite_new_gfn1_calculator(m_ctx, m_tblite_mol);
        } else if (m_method_switch == 2) {
            m_tblite_calc = tblite_new_gfn2_calculator(m_ctx, m_tblite_mol);
        }

        // CRITICAL: Re-configure calculator parameters after recreation
        // Apply same method-specific optimizations as in main calculation
        double effective_damping = m_damping;
        int effective_maxiter = m_SCFmaxiter;
        double effective_accuracy = m_acc;

        if (m_method_switch == 2) { // GFN2-xTB needs stronger settings
            effective_damping = std::max(0.6, m_damping);
            effective_maxiter = std::max(150, m_SCFmaxiter);
            effective_accuracy = std::min(0.5, (double)m_acc);
        }

        tblite_set_calculator_accuracy(m_ctx, m_tblite_calc, effective_accuracy);
        tblite_set_calculator_max_iter(m_ctx, m_tblite_calc, effective_maxiter);
        tblite_set_calculator_mixer_damping(m_ctx, m_tblite_calc, effective_damping);
        tblite_set_calculator_temperature(m_ctx, m_tblite_calc, m_Tele);

        if (m_guess == 0)
            tblite_set_calculator_guess(m_ctx, m_tblite_calc, TBLITE_GUESS_SAD);
        else
            tblite_set_calculator_guess(m_ctx, m_tblite_calc, TBLITE_GUESS_EEQ);
        ApplySolvation();
        tblite_get_singlepoint(m_ctx, m_tblite_mol, m_tblite_calc, m_tblite_res);
        tblite_get_result_energy(m_error, m_tblite_res, &energy);

        // Reset to original verbosity
        if (CurcumaLogger::get_verbosity() >= 3) {
            tblite_set_context_verbosity(m_ctx, 3);
        } else if (CurcumaLogger::get_verbosity() >= 2) {
            tblite_set_context_verbosity(m_ctx, 1);
        } else {
            tblite_set_context_verbosity(m_ctx, 0);
        }
    }

    // Level 1+: Final energy result
    if (CurcumaLogger::get_verbosity() >= 1) {
        CurcumaLogger::energy_abs(energy, "TBLite Energy");
    }

    // Level 2+: Orbital analysis and molecular properties
    if (CurcumaLogger::get_verbosity() >= 2) {
        Vector orbital_energies = OrbitalEnergies();
        Vector orbital_occupations = OrbitalOccupations();

        if (orbital_energies.size() > 0 && orbital_occupations.size() > 0) {
            // Find HOMO/LUMO
            int homo_index = -1;
            for (int i = orbital_occupations.size() - 1; i >= 0; --i) {
                if (orbital_occupations(i) > 0.5) { // Occupied orbital
                    homo_index = i;
                    break;
                }
            }

            if (homo_index >= 0 && homo_index + 1 < orbital_energies.size()) {
                double homo = orbital_energies(homo_index);
                double lumo = orbital_energies(homo_index + 1);
                double gap = lumo - homo;

                CurcumaLogger::param("HOMO", fmt::format("{:.4f} Eh", homo));
                CurcumaLogger::param("LUMO", fmt::format("{:.4f} Eh", lumo));
                CurcumaLogger::param("HOMO-LUMO_gap", fmt::format("{:.4f} Eh ({:.2f} eV)", gap, gap * 27.211));
            }
        }

        // Molecular properties
        Vector charges = Charges();
        if (charges.size() > 0) {
            double total_charge = charges.sum();
            CurcumaLogger::param("total_charge", fmt::format("{:.6f}", total_charge));
        }

        Vector dipole = Dipole();
        if (dipole.size() >= 3) {
            double dipole_magnitude = dipole.norm();
            CurcumaLogger::param("dipole_moment", fmt::format("{:.4f} Debye", dipole_magnitude));
        }
    }

    // Level 3+: Complete orbital listing
    if (CurcumaLogger::get_verbosity() >= 3) {
        Vector orbital_energies = OrbitalEnergies();
        Vector orbital_occupations = OrbitalOccupations();

        if (orbital_energies.size() > 0) {
            CurcumaLogger::info("Complete TBLite orbital energy listing:");
            for (int i = 0; i < orbital_energies.size(); ++i) {
                double occ = (i < orbital_occupations.size()) ? orbital_occupations(i) : 0.0;
                CurcumaLogger::param(fmt::format("Orbital_{}", i + 1),
                    fmt::format("{:.6f} Eh (occ={:.3f})", orbital_energies(i), occ));
            }
        }
    }

    // Handle gradient calculation
    if (gradient) {
        tblite_get_result_gradient(m_error, m_tblite_res, m_gradient.data());
        m_gradient *= au;
        if (tblite_check_context(m_ctx)) {
            tbliteContextError();
        }

        if (CurcumaLogger::get_verbosity() >= 2) {
            double grad_norm = m_gradient.norm();
            CurcumaLogger::param("gradient_norm", fmt::format("{:.6f} Eh/Bohr", grad_norm));
        }
    }

    if (count == 1)
        tblite_delete_container(&m_tb_cont);
    return energy;
}

void TBLiteInterface::setMethod(const std::string& method)
{
    CurcumaLogger::info("TBLiteInterface::setMethod called");
    CurcumaLogger::param("requested_method", method);

    QMInterface::setMethod(method);
    CurcumaLogger::param("normalized_method", m_method);

    TBLiteMethod old_method = m_tblite_method;
    if (m_method.compare("ipea1") == 0) {
        m_tblite_method = TBLiteMethod::IPEA1;
        m_method_switch = static_cast<int>(TBLiteMethod::IPEA1);
        CurcumaLogger::info("Method resolved to: iPEA1 (TBLiteMethod::IPEA1)");
    } else if (m_method.compare("gfn1") == 0) {
        m_tblite_method = TBLiteMethod::GFN1;
        m_method_switch = static_cast<int>(TBLiteMethod::GFN1);
        CurcumaLogger::info("Method resolved to: GFN1-xTB (TBLiteMethod::GFN1)");
    } else if (m_method.compare("gfn2") == 0) {
        m_tblite_method = TBLiteMethod::GFN2;
        m_method_switch = static_cast<int>(TBLiteMethod::GFN2);
        CurcumaLogger::info("Method resolved to: GFN2-xTB (TBLiteMethod::GFN2)");
    } else {
        CurcumaLogger::warn("Unknown TBLite method '" + method + "', keeping current method");
    }

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::param("method_old", std::to_string(static_cast<int>(old_method)));
        CurcumaLogger::param("method_new", std::to_string(static_cast<int>(m_tblite_method)));
        CurcumaLogger::param("m_method_switch", std::to_string(m_method_switch));
    }
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
    CurcumaLogger::error(fmt::format("TBLite Error: {}", message));
    tblite_clear_error(m_error);
}

void TBLiteInterface::tbliteContextError()
{
    char message[512];
    tblite_get_context_error(m_ctx, message, NULL);
    CurcumaLogger::error(fmt::format("TBLite Context Error: {}", message));
    tblite_clear_error(m_error);
}
