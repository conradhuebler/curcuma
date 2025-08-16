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

#ifdef USE_GFNFF
#include "gfnff_interface_c.h"
#endif

static json GFNFFSettings{
    { "charge", 0 },
    { "printlevel", 1 },
    { "solvent", "none" },
    { "verbose", 0 }
};

class GFNFFInterface : public QMInterface {
public:
    GFNFFInterface(const json& gfnffsettings = GFNFFSettings);
    ~GFNFFInterface();

    bool InitialiseMolecule(const Mol& mol) override;
    bool InitialiseMolecule(const int* attyp, const double* coord, const int natoms, const double charge, const int spin) override;
    bool InitialiseMolecule() override;
    double Calculation(bool gradient = false, bool verbose = false) override;
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
    json m_settings;
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

    void updateVerbosity();
};