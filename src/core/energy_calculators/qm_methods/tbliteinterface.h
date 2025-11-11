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

#pragma once

#include "interface/abstract_interface.h"
#include "src/core/parameter_macros.h"
#include "src/core/config_manager.h"

#ifdef USE_TBLITE
#include "tblite.h"
#endif

// Claude Generated 2025: TBLite Parameter Registry - replaces static TBLiteSettings JSON
BEGIN_PARAMETER_DEFINITION(tblite)
    // SCF Parameters
    PARAM(accuracy, Int, 1, "Accuracy level for TBLite calculations (0=crude, 1=normal, 2=tight).", "SCF", {"tb_acc"})
    PARAM(max_iterations, Int, 100, "Maximum number of SCF iterations.", "SCF", {"SCFmaxiter"})
    PARAM(damping, Double, 0.4, "SCF damping parameter for convergence.", "SCF", {"tb_damping"})
    PARAM(electronic_temperature, Double, 300.0, "Electronic temperature in Kelvin for Fermi smearing.", "SCF", {"Tele"})
    PARAM(initial_guess, String, "SAD", "Initial guess method (SAD=Superposition of Atomic Densities).", "SCF", {"tb_guess"})

    // Molecular Properties
    PARAM(spin, Double, 0.0, "Total spin of the system (0.0 = singlet).", "Molecular", {})

    // Solvent Model Parameters
    PARAM(solvent_model, Int, 0, "Solvent model type (0=none, 1=CPCM, 2=GB, 3=ALPB).", "Solvation", {})
    PARAM(solvent, String, "none", "Solvent name (e.g., 'water', 'acetone', 'none').", "Solvation", {})
    PARAM(solvent_epsilon, Double, -1.0, "Solvent dielectric constant (epsilon). -1 = auto-detect from name.", "Solvation", {"solvent_eps"})
    PARAM(solvent_gb_version, Int, 0, "Generalized Born model version (0=GBSA, 1=ALPB).", "Solvation", {})
    PARAM(solvent_gb_kernel, Int, 1, "GB kernel function type.", "Solvation", {})
    PARAM(solvent_alpb_version, Int, 12, "ALPB solvation model version.", "Solvation", {})
    PARAM(solvent_alpb_reference, Int, 1, "ALPB reference state.", "Solvation", {})
END_PARAMETER_DEFINITION

class TBLiteInterface : public QMInterface {
public:
    enum class TBLiteMethod {
        IPEA1 = 0,
        GFN1 = 1,
        GFN2 = 2
    };
    TBLiteInterface(const ConfigManager& config);
    ~TBLiteInterface();

    bool InitialiseMolecule(const Mol& mol) override;
    bool InitialiseMolecule(const Mol* mol) override;
    bool InitialiseMolecule(const int* attyp, const double* coord, const int natoms, const double charge, const int spin) override;
    bool InitialiseMolecule();

    bool UpdateMolecule(const Mol& mol) override;
    bool UpdateMolecule(const Mol* mol) override;
    bool UpdateMolecule(const Matrix& geometry) override;
    bool UpdateMolecule(const double* coord) override;
    bool UpdateMolecule();

    bool Error() override { return m_error_count >= 10; }
    double Calculation(bool gradient = 0);

    void clear() override;

    Vector Charges() const override;
    Vector Dipole() const override;

    Vector BondOrders() const;
    Vector OrbitalEnergies() const override;
    Vector OrbitalOccupations() const override;
    void setMethod(const std::string& method) override;
    virtual bool hasGradient() const { return true; }

private:
    void ApplySolvation();

    void tbliteError();
    void tbliteContextError();
    double* m_coord;
    int* m_attyp;

    double m_thr = 1.0e-10;
    int m_acc = 2;
    int m_SCFmaxiter = 100;
    int m_guess = 0;
    int m_error_count = 0;
    double m_damping = 0.5;
    double m_Tele = 298;
    double m_solvent_eps = -1;
    int m_solvent_model = 0;
    char* m_solvent = NULL;
    int m_solvent_gb_version = 0;
    int m_solvent_gb_kernel = 1;
    int m_solvent_alpb_version = 12;
    int m_solvent_alpb_reference = 1;

    tblite_error m_error = NULL;
    tblite_structure m_tblite_mol = NULL;
    tblite_result m_tblite_res = NULL;
    tblite_context m_ctx = NULL;
    tblite_calculator m_tblite_calc = NULL;
    tblite_container m_tb_cont = NULL;

    bool m_initialised = false, m_calculator = false;
    mutable ConfigManager m_config;

    // Method selection - Claude Generated improvements
    TBLiteMethod m_tblite_method = TBLiteMethod::IPEA1; // Type-safe enum
    int m_method_switch = 0; // Legacy compatibility
};
