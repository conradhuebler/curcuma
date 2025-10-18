/*
 * < C++ Ulysses Interface >
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include <string>
#include <vector>

// Claude Generated 2025: Ulysses Parameter Registry - replaces static UlyssesSettings JSON
BEGIN_PARAMETER_DEFINITION(ulysses)
    // Method Selection
    PARAM(method, String, "GFN2", "Semi-empirical method (GFN1, GFN2, iPEA1, PM3, AM1, MNDO, MNDO/D, PM6, OM2, OM3).", "Method", {})

    // SCF Parameters
    PARAM(max_iterations, Int, 100, "Maximum number of SCF iterations.", "SCF", {"SCFmaxiter"})
    PARAM(electronic_temperature, Double, 300.0, "Electronic temperature in Kelvin for Fermi smearing.", "SCF", {"Tele"})

    // Molecular Properties
    PARAM(multiplicity, Int, 1, "Spin multiplicity (1=singlet, 2=doublet, 3=triplet, etc.).", "Molecular", {"mult"})

    // Solvation
    PARAM(solvent, String, "none", "Solvent name for implicit solvation (water, acetone, dmso, none, etc.).", "Solvation", {})
END_PARAMETER_DEFINITION

#ifdef USE_ULYSSES
class UlyssesObject;
#endif

class UlyssesInterface : public QMInterface {
public:
    UlyssesInterface(const ConfigManager& config);
    ~UlyssesInterface();

    virtual bool InitialiseMolecule() override;
    virtual bool UpdateMolecule(const Geometry& geometry) override;

    double Calculation(bool gradient = false) override;
    virtual Vector Charges() const override;
    virtual Vector OrbitalEnergies() const override;
    virtual bool hasGradient() const { return true; }

private:
    Mol m_mol;
#ifdef USE_ULYSSES
    UlyssesObject* m_ulysses = nullptr;
#endif

    double m_Tele = 300;
    double m_SCFmaxiter = 100;
    int m_mult = 1;

    StringList m_solvents = { "acetone", "acetone", "acetonitrile", "aniline", "benzaldehyde", "benzene", "dichloromethane", "chloroform", "carbon disulfide", "dioxane", "dmf", "dmso", "ethanol", "diethyl ether", "ethyl acetate", "furane", "hexadecane", "hexane", "methanol", "nitromethane", "octanol", "phenol", "thf", "toluene", "water", "octanol wet" };
    std::string m_solvent = "none";
    std::string m_correction = "0"; // Claude Generated: Store correction parameter separately

    mutable ConfigManager m_config;
};
