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

#pragma once

#include "interface/abstract_interface.h"

#include "src/tools/general.h"
#include "src/core/parameter_macros.h"
#include "src/core/config_manager.h"

#include "external/xtb/include/xtb.h"

// Claude Generated 2025: XTB Parameter Registry - replaces static XTBSettings JSON
BEGIN_PARAMETER_DEFINITION(xtb)
    PARAM(accuracy, Int, 2, "Accuracy level for XTB calculations (0=crude, 1=sloppy, 2=normal, 3=tight, 4=vtight).", "SCF", {"xtb_ac"})
    PARAM(max_iterations, Int, 100, "Maximum number of SCF iterations.", "SCF", {"SCFmaxiter"})
    PARAM(electronic_temperature, Double, 300.0, "Electronic temperature in Kelvin for Fermi smearing.", "SCF", {"Tele"})
    PARAM(spin, Double, 0.0, "Total spin of the system (0.0 = singlet).", "Molecular", {})
END_PARAMETER_DEFINITION

class UFF;

class XTBInterface : public QMInterface {
public:
    XTBInterface(const ConfigManager& config);
    ~XTBInterface();

    bool InitialiseMolecule(const Mol& molecule);
    bool InitialiseMolecule(const Mol* molecule);
    bool InitialiseMolecule(const int* attyp, const double* coord, const int natoms, const double charge, const int spin);
    bool InitialiseMolecule();
    bool UpdateMolecule(const Mol& molecule);
    bool UpdateMolecule(const Mol* mol);
    bool UpdateMolecule(const double* coord);
    bool UpdateMolecule(const Matrix& geometry);
    bool UpdateMolecule();

    /* int parameter
     * 66 = xtb GFN FF
     * 0 = xtb GFN 0
     * 1 = xtb GFN 1
     * 2 = xtb GFN 2
     * */
    double Calculation(bool gradient = 0);

    void clear();

    Vector Charges() const;
    Vector Dipole() const;
    Vector BondOrders() const;
    Vector OrbitalEnergies() const;
    Vector OrbitalOccupations() const;

    void setMethod(const std::string& method) override;
    virtual bool hasGradient() const { return true; }

private:
    double m_thr = 1.0e-10;
    double* m_coord;
    int* m_attyp;
    int m_accuracy = 1, m_SCFmaxiter = 100;
    double m_Tele = 298;
    xtb_TEnvironment m_env = NULL;
    xtb_TMolecule m_xtb_mol = NULL;
    xtb_TCalculator m_xtb_calc = NULL;
    xtb_TResults m_xtb_res = NULL;
    bool m_initialised = false, m_setup = false;
    mutable ConfigManager m_config;
};
