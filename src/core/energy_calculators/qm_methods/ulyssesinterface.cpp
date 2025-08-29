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

#include <string>
#include <vector>

#include "interface/abstract_interface.h"
#ifdef USE_ULYSSES
#include "interface/ulysses.h"
#endif

#include "src/core/curcuma_logger.h"
#include "ulyssesinterface.h"
#include <fmt/format.h>

UlyssesInterface::UlyssesInterface(const json& ulyssessettings)
{
    m_ulyssessettings = MergeJson(UlyssesSettings, ulyssessettings);
    m_Tele = m_ulyssessettings["Tele"];
    m_SCFmaxiter = m_ulyssessettings["SCFmaxiter"];
    m_solvent = m_ulyssessettings["solvent"];
    m_mult = m_ulyssessettings["mult"];

    // Claude Generated: Use parsed base_method and corecorrection from ulysses_method.cpp
    if (ulyssessettings.contains("base_method") && ulyssessettings.contains("corecorrection")) {
        m_method = ulyssessettings["base_method"];
        m_correction = ulyssessettings["corecorrection"];

        if (CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::info("UlyssesInterface using parsed parameters");
            CurcumaLogger::param("base_method", m_method);
            CurcumaLogger::param("correction", m_correction);
        }
    } else {
        // Fallback to original method for backward compatibility
        m_method = m_ulyssessettings["method"];
        m_correction = "0";

        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info("UlyssesInterface fallback mode");
            CurcumaLogger::param("original_method", m_method);
        }
    }

    // Validate solvent with CurcumaLogger
    if (std::find(m_solvents.begin(), m_solvents.end(), m_solvent) == m_solvents.end()) {
        CurcumaLogger::warn("Solvent '" + m_solvent + "' not supported by Ulysses - using 'none'");
        m_solvent = "none";
    }
#ifdef USE_ULYSSES
    m_ulysses = new UlyssesObject();
#endif
}

UlyssesInterface::~UlyssesInterface()
{
#ifdef USE_ULYSSES
    delete m_ulysses;
#endif
}

bool UlyssesInterface::InitialiseMolecule()
{
#ifdef USE_ULYSSES
    // Claude Generated: Pass base method only, no string parsing needed
    std::string base_method = m_method; // This is already the base method from parameter parsing
    m_ulysses->setMethod(base_method);

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("UlyssesObject method setup");
        CurcumaLogger::param("base_method_passed", base_method);
        CurcumaLogger::param("correction_to_apply", m_correction);
    }

    // Verbosity Level 1+: Method initialization info
    if (CurcumaLogger::get_verbosity() >= 1) {
        CurcumaLogger::info("Initializing Ulysses quantum chemistry method");
        CurcumaLogger::param("method", m_method);
        CurcumaLogger::param("SCF_maxiter", m_SCFmaxiter);
        CurcumaLogger::param("temperature", m_Tele);
        if (m_correction != "0") {
            CurcumaLogger::param("correction", m_correction);
        }
        if (m_solvent != "none") {
            CurcumaLogger::param("solvent", m_solvent);
        }
    }

    m_ulysses->setMolecule(m_geometry, m_atoms, m_charge, m_mult, "C1", m_correction);

    return true;
#else
    CurcumaLogger::error("Ulysses interface not compiled - USE_ULYSSES not defined");
    return false;
#endif
}

bool UlyssesInterface::UpdateMolecule(const Geometry& geometry)
{
#ifdef USE_ULYSSES
    m_geometry = geometry;
    m_ulysses->UpdateGeometry(geometry);
    return true;
#else
    return false;
#endif
}

double UlyssesInterface::Calculation(bool gradient)
{
#ifdef USE_ULYSSES
    m_ulysses->setTele(m_Tele);
    m_ulysses->setMaxIter(m_SCFmaxiter);
    m_ulysses->setSolvent(m_solvent);

    // Start calculation with verbosity control
    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("Starting Ulysses SCF calculation");
        if (gradient) {
            CurcumaLogger::param("gradient", "analytical");
        }
    }

    // Perform calculation - Ulysses handles its own output currently
    // TODO: Capture and filter Ulysses output based on verbosity
    m_ulysses->Calculate(gradient, CurcumaLogger::get_verbosity() >= 3);

    double energy = m_ulysses->Energy();

    // Final results - Level 1+
    if (CurcumaLogger::get_verbosity() >= 1) {
        CurcumaLogger::energy_abs(energy, "SCF Energy");
    }

    // Orbital properties - Level 2+
    if (CurcumaLogger::get_verbosity() >= 2) {
        Vector orbital_energies = OrbitalEnergies();
        if (orbital_energies.size() > 0) {
            // Find HOMO and LUMO - rough estimate: half-filled orbitals
            int homo_index = -1;
            int num_electrons = m_atomcount; // Rough estimate - neutral atoms
            int num_occupied = num_electrons / 2; // Closed shell approximation

            if (num_occupied > 0 && num_occupied <= orbital_energies.size()) {
                homo_index = num_occupied - 1; // HOMO is last occupied orbital
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
    }

    // Full orbital listing - Level 3+
    if (CurcumaLogger::get_verbosity() >= 3) {
        Vector orbital_energies = OrbitalEnergies();
        if (orbital_energies.size() > 0) {
            CurcumaLogger::info("Complete orbital energy listing:");
            for (int i = 0; i < orbital_energies.size(); ++i) {
                CurcumaLogger::param(fmt::format("Orbital_{}", i + 1),
                    fmt::format("{:.6f} Eh", orbital_energies(i)));
            }
        }
    }

    if (gradient) {
        m_gradient = m_ulysses->Gradient();
        if (CurcumaLogger::get_verbosity() >= 2) {
            double grad_norm = m_gradient.norm();
            CurcumaLogger::param("gradient_norm", fmt::format("{:.6f} Eh/Bohr", grad_norm));
        }
    }

    return energy;
#else
    CurcumaLogger::error("Ulysses calculation failed - not compiled");
    return 0;
#endif
}

Vector UlyssesInterface::Charges() const
{
#ifdef USE_ULYSSES
    return m_ulysses->Charges();
#else
    return Vector{};
#endif
}

Vector UlyssesInterface::OrbitalEnergies() const
{
#ifdef USE_ULYSSES
    return m_ulysses->OrbitalEnergies();
#else
    return Vector{};
#endif
}