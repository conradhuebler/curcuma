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

#include "src/tools/general.h"
#include "src/core/parameter_macros.h"
#include "src/core/config_manager.h"

#include "dftd_damping.h"
#include "dftd_dispersion.h"
#include "dftd_econv.h"
#include "dftd_geometry.h"

#include "interface/abstract_interface.h"

#include "src/core/global.h"

// Claude Generated 2025: DFT-D4 Parameter Registry - replaces static DFTD4Settings JSON
BEGIN_PARAMETER_DEFINITION(dftd4)
    // Scaling Parameters (C6, C8, C10 terms + 3-body)
    PARAM(s6, Double, 1.00, "Scaling factor for C6 dispersion term.", "Scaling", {"d4_s6"})
    PARAM(s8, Double, 1.20065498, "Scaling factor for C8 dispersion term.", "Scaling", {"d4_s8"})
    PARAM(s9, Double, 1.0, "Scaling factor for three-body ATM dispersion term.", "Scaling", {"d4_s9"})
    PARAM(s10, Double, 0.0, "Scaling factor for C10 dispersion term (rarely used).", "Scaling", {"d4_s10"})

    // Damping Function Parameters (optimized for D4)
    PARAM(a1, Double, 0.40085597, "Damping parameter a1 for Becke-Johnson damping.", "Damping", {"d4_a1"})
    PARAM(a2, Double, 5.02928789, "Damping parameter a2 in Bohr for Becke-Johnson damping.", "Damping", {"d4_a2"})
    PARAM(alpha, Double, 16.0, "Alpha parameter for damping function.", "Damping", {"d4_alp"})

    // Functional and Options
    PARAM(functional, String, "pbe0", "DFT functional name for D4 dispersion parameters (e.g., pbe0, b3lyp, wb97x).", "General", {"d4_func"})
    PARAM(three_body, Bool, true, "Include Axilrod-Teller-Muto three-body dispersion term (recommended for D4).", "General", {"d4_atm"})
END_PARAMETER_DEFINITION

class DFTD4Interface : public QMInterface {
public:
    DFTD4Interface(const ConfigManager& config);
    DFTD4Interface();

    ~DFTD4Interface();

    bool InitialiseMolecule(const Mol& mol, double factor = 1);
    bool InitialiseMolecule(const std::vector<int>& atomtype);

    double Calculation(bool gradient = 0) override;

    // Claude Generated 2025: ConfigManager migration - Phase 2.2
    void UpdateParameters(const ConfigManager& config);

    // Backward compatibility: JSON version delegates to ConfigManager
    void UpdateParameters(const json& controller);

    void clear();
    dftd4::dparam Parameter() const { return m_par; }

    inline void UpdateGeometry(const double* coord)
    {
        for (int i = 0; i < m_natoms; ++i) {
            m_mol.xyz(i, 0) = coord[3 * i + 0];
            m_mol.xyz(i, 1) = coord[3 * i + 1];
            m_mol.xyz(i, 2) = coord[3 * i + 2];
        }
    }

    inline void UpdateAtom(int index, double x, double y, double z, int element = -1)
    {
        m_mol.UpdateAtom(index, x, y, z, element);
    }

    void PrintParameter() const;

    virtual bool hasGradient() const { return true; }

private:
    dftd4::dparam m_par;
    dftd4::TMolecule m_mol;
    int m_charge = 0;
    int m_natoms;
    mutable ConfigManager m_config;
};
