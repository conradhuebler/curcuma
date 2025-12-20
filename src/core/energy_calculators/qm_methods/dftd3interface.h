/*
 * < C++ XTB and tblite Interface >
 * Copyright (C) 2020 - 2024 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#ifdef USE_D3

#include "src/tools/general.h"
#include "src/core/parameter_macros.h"
#include "src/core/config_manager.h"
#ifdef USE_D3
#include "s-dftd3.h"
#endif

#include "interface/abstract_interface.h"
#include "src/core/global.h"

// Claude Generated 2025: DFT-D3 Parameter Registry - replaces static DFTD3Settings JSON
BEGIN_PARAMETER_DEFINITION(dftd3)
    // Scaling Parameters (C6, C8, C9 terms)
    PARAM(s6, Double, 1.0, "Scaling factor for C6 dispersion term.", "Scaling", {"d_s6"})
    PARAM(s8, Double, 0.0, "Scaling factor for C8 dispersion term.", "Scaling", {"d_s8"})
    PARAM(s9, Double, 0.0, "Scaling factor for three-body ATM dispersion term.", "Scaling", {"d_s9"})

    // Damping Function Parameters
    PARAM(damping, String, "bj", "Damping function type (bj=Becke-Johnson, zero=zero-damping, op=optimized power).", "Damping", {"d_damping"})
    PARAM(a1, Double, 0.0, "Damping parameter a1 (distance scaling).", "Damping", {"d_a1"})
    PARAM(a2, Double, 0.0, "Damping parameter a2 (distance offset in Bohr).", "Damping", {"d_a2"})
    PARAM(alpha, Double, 14.0, "Alpha parameter for damping function.", "Damping", {"d_alp"})
    PARAM(beta, Double, 8.0, "Beta exponent for damping function.", "Damping", {"d_bet"})

    // Functional and Options
    PARAM(functional, String, "pbe0", "DFT functional name for dispersion parameters (e.g., pbe0, b3lyp, tpss).", "General", {"d_func"})
    PARAM(three_body, Bool, false, "Include Axilrod-Teller-Muto three-body dispersion term.", "General", {"d_atm"})
END_PARAMETER_DEFINITION

class DFTD3Interface : public QMInterface {
public:
    DFTD3Interface(const ConfigManager& config);
    DFTD3Interface();

    ~DFTD3Interface();

    bool InitialiseMolecule(const std::vector<int>& atomtypes);

    double Calculation(bool gradient = 0) override;

    // Claude Generated 2025: ConfigManager migration - Phase 2.1
    void UpdateParameters(const ConfigManager& config);
    void UpdateParametersD3(const ConfigManager& config);

    // Backward compatibility: JSON versions delegate to ConfigManager
    void UpdateParameters(const json& controller);
    void UpdateParametersD3(const json& controller);

    void clear();

    void UpdateGeometry(const double* coord);

    void PrintParameter() const;
    void UpdateAtom(int index, double x, double y, double z);

    inline double ParameterA1() const { return m_d3_a1; }
    inline double ParameterA2() const { return m_d3_a2; }
    inline double ParameterAlp() const { return m_d3_alp; }
    inline double ParameterS6() const { return m_d3_s6; }
    inline double ParameterS8() const { return m_d3_s8; }
    inline double ParameterS9() const { return m_d3_s9; }
    virtual bool hasGradient() const { return true; }

private:
    void CreateParameter();

    std::string m_damping;
    std::string m_functional;

    int m_charge = 0;
    int m_natoms = 0;  // Number of atoms (set in InitialiseMolecule)
    double* m_coord;
    int* m_attyp;

    double m_d3_a1 = 1;
    double m_d3_a2 = 1;
    double m_d3_alp = 1;

    double m_d3_s6 = 1;
    double m_d3_s8 = 1;
    double m_d3_s9 = 1;

    double m_bet = 8;
    bool m_atm = false;

#ifdef USE_D3
    dftd3_error m_error;
    dftd3_structure m_mol;
    dftd3_model m_disp;
    dftd3_param m_param;
#endif
    mutable ConfigManager m_config;
};

#endif // USE_D3
