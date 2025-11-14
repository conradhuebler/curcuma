/*
 * < C++ External GFN-FF Interface >
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
 * Claude Generated: External GFN-FF interface wrapper for C library
 */

#pragma once

#include "interface/abstract_interface.h"
#include "src/core/parameter_macros.h"
#include "src/core/config_manager.h"
#include "src/global_config.h"

#ifdef USE_GFNFF
#include "gfnff_interface_c.h"
#endif

// Claude Generated 2025: External GFN-FF Parameter Registry - replaces static GFNFFSettings JSON
BEGIN_PARAMETER_DEFINITION(gfnff_external)
    // Molecular Properties
    PARAM(charge, Int, 0, "Total molecular charge.", "Molecular", {})

    // Output Control
    PARAM(print_level, Int, 1, "Verbosity level for GFN-FF output (0=silent, 1=minimal, 2=verbose).", "Output", {"printlevel"})

    // Solvation
    PARAM(solvent, String, "none", "Solvent name for implicit solvation (none, water, etc.).", "Solvation", {})
END_PARAMETER_DEFINITION

class GFNFFInterface : public QMInterface {
public:
    GFNFFInterface(const ConfigManager& config);
    ~GFNFFInterface();

    bool InitialiseMolecule(const Mol& mol) override;
    bool InitialiseMolecule(const int* attyp, const double* coord, const int natoms, const double charge, const int spin) override;
    bool InitialiseMolecule() override;
    double Calculation(bool gradient = false) override;
    void clear() override;

    Vector BondOrders() const override;
    Vector Charges() const override;
    Geometry Gradient() const override;
    bool hasGradient() const override { return true; }

    void SetLogLevel(int level) { (void)level; } // Suppress unused parameter warning
    double Energy() const { return m_energy; }

private:
#ifdef USE_GFNFF
    c_gfnff_calculator m_calculator;
#endif
    bool m_initialized;
    double m_energy;
    Geometry m_gradient;
    Vector m_charges;
    int m_natoms;
    int m_charge;
    int m_printlevel;
    std::string m_solvent;
    std::vector<int> m_atom_types;
    std::vector<double> m_coordinates; // Current coordinates in Angstrom
    std::vector<double> m_coordinates_bohr; // Coordinates in Bohr for GFN-FF
    mutable ConfigManager m_config;

    void updateVerbosity();
};