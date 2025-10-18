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

#include "gfnffinterface.h"

#include "src/core/curcuma_logger.h"
#include "src/core/global.h"
#include <fmt/format.h>

#include <cstring>
#include <iostream>
#include <math.h>
#include <stdio.h>

GFNFFInterface::GFNFFInterface(const ConfigManager& config)
    : m_config(config)
    , m_initialized(false)
    , m_energy(0.0)
    , m_natoms(0)
    , m_charge(0)
    , m_printlevel(1)
    , m_solvent("none")
{
#ifdef USE_GFNFF
    // Initialize calculator struct
    m_calculator.ptr = nullptr;
#endif

    // Claude Generated 2025: ConfigManager migration - Phase 3B
    m_charge = m_config.get<int>("charge", 0);
    m_printlevel = m_config.get<int>("print_level", 1);
    m_solvent = m_config.get<std::string>("solvent", "none");

    updateVerbosity();

    // Level 1+: GFN-FF initialization
    if (CurcumaLogger::get_verbosity() >= 1) {
        CurcumaLogger::info("Initializing external GFN-FF quantum chemistry method");
        CurcumaLogger::param("charge", m_charge);
        CurcumaLogger::param("solvent", m_solvent);
        CurcumaLogger::param("printlevel", m_printlevel);
    }
}

GFNFFInterface::~GFNFFInterface()
{
#ifdef USE_GFNFF
    if (m_calculator.ptr != nullptr) {
        c_gfnff_calculator_deallocate(&m_calculator);
        m_calculator.ptr = nullptr;
    }
#endif
}

bool GFNFFInterface::InitialiseMolecule(const Mol& mol)
{
    // Following TBLiteInterface pattern - extract coordinates properly from Eigen matrix
    m_natoms = mol.m_number_atoms;
    std::vector<int> atoms = mol.m_atoms;
    std::vector<double> coords(3 * m_natoms);

    for (int i = 0; i < m_natoms; ++i) {
        coords[3 * i + 0] = mol.m_geometry(i, 0); // x coordinate
        coords[3 * i + 1] = mol.m_geometry(i, 1); // y coordinate
        coords[3 * i + 2] = mol.m_geometry(i, 2); // z coordinate
    }

    return InitialiseMolecule(atoms.data(), coords.data(), mol.m_number_atoms, mol.m_charge, mol.m_spin);
}

bool GFNFFInterface::InitialiseMolecule(const int* attyp, const double* coord, const int natoms, const double charge, const int spin)
{
    m_natoms = natoms;
    m_charge = charge;
    (void)spin; // Suppress unused parameter warning

    // Store coordinates (convert to Bohr for GFN-FF)
    m_gradient = Geometry::Zero(m_natoms, 3);
    m_charges = Vector::Zero(m_natoms);

    // Store atom types and coordinates (following gfnff_helper pattern)
    m_atom_types.assign(attyp, attyp + m_natoms);
    m_coordinates.assign(coord, coord + 3 * m_natoms); // Store in Angstrom

    // Convert and store coordinates in Bohr for GFN-FF (following gfnff_helper pattern)
    m_coordinates_bohr.resize(3 * m_natoms);
    for (int i = 0; i < m_natoms; ++i) {
        m_coordinates_bohr[3 * i + 0] = coord[3 * i + 0] / au;
        m_coordinates_bohr[3 * i + 1] = coord[3 * i + 1] / au;
        m_coordinates_bohr[3 * i + 2] = coord[3 * i + 2] / au;
    }

#ifdef USE_GFNFF
    // Clean up existing calculator if any
    if (m_calculator.ptr != nullptr) {
        c_gfnff_calculator_deallocate(&m_calculator);
        m_calculator.ptr = nullptr;
    }

    // Level 2+: Method setup info
    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("Setting up external GFN-FF calculator");
        CurcumaLogger::param("natoms", m_natoms);
        CurcumaLogger::param("charge", m_charge);
        CurcumaLogger::param("coordinates", "converted to Bohr");
    }

    // Initialize external GFN-FF calculator (using stored Bohr coordinates)
    // Note: API returns calculator object, not status code
    double (*coord_ptr)[3] = reinterpret_cast<double (*)[3]>(m_coordinates_bohr.data());
    m_calculator = c_gfnff_calculator_init(m_natoms, m_atom_types.data(), coord_ptr, m_charge, m_printlevel, m_solvent.c_str());

    // Check if calculator was initialized (assuming null ptr means failure)
    if (m_calculator.ptr == nullptr) {
        CurcumaLogger::error("Failed to initialize external GFN-FF calculator");
        return false;
    }

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::success("External GFN-FF calculator initialized successfully");
    }

    m_initialized = true;
    return true;
#else
    CurcumaLogger::error("External GFN-FF not available - USE_GFNFF not compiled");
    return false;
#endif
}

bool GFNFFInterface::InitialiseMolecule()
{
    CurcumaLogger::error("InitialiseMolecule() without parameters not supported for external GFN-FF");
    return false;
}

double GFNFFInterface::Calculation(bool gradient)
{
    // Level 3+: Entry point debug
    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("=== GFNFFInterface::Calculation() called ===");
        CurcumaLogger::param("gradient_requested", gradient);
        CurcumaLogger::param("curcuma_verbosity", CurcumaLogger::get_verbosity());
    }

#ifdef USE_GFNFF
    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("USE_GFNFF is defined");
    }

    if (!m_initialized || m_calculator.ptr == nullptr) {
        CurcumaLogger::error("External GFN-FF calculator not initialized");
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::param("m_initialized", m_initialized);
            CurcumaLogger::param("calculator_ptr_null", (m_calculator.ptr == nullptr));
        }
        return 0.0;
    }

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::success("External GFN-FF calculator appears to be initialized");
        CurcumaLogger::param("m_initialized", m_initialized);
        CurcumaLogger::param("calculator_ptr", fmt::format("{}", (void*)m_calculator.ptr));
        CurcumaLogger::param("m_natoms", m_natoms);
        CurcumaLogger::param("m_charge", m_charge);
        CurcumaLogger::param("atom_types_size", static_cast<int>(m_atom_types.size()));
        CurcumaLogger::param("coordinates_size", static_cast<int>(m_coordinates.size()));
    }

    updateVerbosity();

    // Level 2+: Calculation start
    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("Starting external GFN-FF calculation");
        if (gradient) {
            CurcumaLogger::param("gradient", "analytical");
        }
    }

    double energy = 0.0;
    // ALWAYS allocate gradient array - external GFN-FF requires it even if gradient not requested
    std::vector<double> grad_data(3 * m_natoms, 0.0);

    // Level 3+: Coordinate preparation debug
    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("Using stored Bohr coordinates for external GFN-FF calculation");
        CurcumaLogger::param("coordinate_source", "stored from initialization");
    }

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("Sample coordinates (first atom):");
        if (m_natoms > 0) {
            CurcumaLogger::param("atom_0_angstrom", fmt::format("[{:.6f}, {:.6f}, {:.6f}]", m_coordinates[0], m_coordinates[1], m_coordinates[2]));
            CurcumaLogger::param("atom_0_bohr", fmt::format("[{:.6f}, {:.6f}, {:.6f}]", m_coordinates_bohr[0], m_coordinates_bohr[1], m_coordinates_bohr[2]));
        }
        CurcumaLogger::param("atom_types_first_3", fmt::format("[{}, {}, {}]", m_atom_types.size() > 0 ? m_atom_types[0] : -1, m_atom_types.size() > 1 ? m_atom_types[1] : -1, m_atom_types.size() > 2 ? m_atom_types[2] : -1));
    }
    // Following gfnff_helper pattern exactly - use SAME coordinates as in initialization
    double (*coord_ptr)[3] = reinterpret_cast<double (*)[3]>(m_coordinates_bohr.data());
    // ALWAYS pass gradient array - external GFN-FF requires it even if not requested
    double (*grad_ptr)[3] = reinterpret_cast<double (*)[3]>(grad_data.data());

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("About to call c_gfnff_calculator_singlepoint");
        CurcumaLogger::param("calculator_ptr", fmt::format("{}", (void*)m_calculator.ptr));
        CurcumaLogger::param("natoms", m_natoms);
        CurcumaLogger::param("coord_ptr", fmt::format("{}", (void*)coord_ptr));
        CurcumaLogger::param("energy_ptr", fmt::format("{}", (void*)&energy));
        CurcumaLogger::param("gradient_requested", gradient);
        CurcumaLogger::param("grad_ptr", fmt::format("{}", (void*)grad_ptr));
    }

    int iostat = 0;
    CurcumaLogger::info("Calling external GFN-FF singlepoint calculation...");
    c_gfnff_calculator_singlepoint(&m_calculator, m_natoms, m_atom_types.data(), coord_ptr, &energy, grad_ptr, &iostat);

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("External GFN-FF singlepoint call completed");
        CurcumaLogger::param("iostat", iostat);
        CurcumaLogger::param("energy_result", fmt::format("{:.8f}", energy));
    }

    if (iostat != 0) {
        CurcumaLogger::error(fmt::format("External GFN-FF calculation failed (iostat: {})", iostat));
        return 0.0;
    }

    m_energy = energy;

    // Level 1+: Final energy result
    if (CurcumaLogger::get_verbosity() >= 1) {
        CurcumaLogger::energy_abs(energy, "External GFN-FF Energy");
    }

    // Handle gradient
    if (gradient && grad_ptr != nullptr) {
        for (int i = 0; i < m_natoms; ++i) {
            m_gradient(i, 0) = grad_data[3 * i + 0] * au; // Convert back to Angstrom
            m_gradient(i, 1) = grad_data[3 * i + 1] * au;
            m_gradient(i, 2) = grad_data[3 * i + 2] * au;
        }

        if (CurcumaLogger::get_verbosity() >= 2) {
            double grad_norm = m_gradient.norm();
            CurcumaLogger::param("gradient_norm", fmt::format("{:.6f} Eh/Bohr", grad_norm));
        }
    }

    // Level 2+: Molecular properties
    if (CurcumaLogger::get_verbosity() >= 2) {
        if (m_charges.size() > 0) {
            double total_charge = m_charges.sum();
            CurcumaLogger::param("total_charge", fmt::format("{:.6f}", total_charge));
        }
    }

    return energy;
#else
    CurcumaLogger::error("External GFN-FF not available - USE_GFNFF not compiled");
    return 0.0;
#endif
}

void GFNFFInterface::clear()
{
#ifdef USE_GFNFF
    if (m_calculator.ptr != nullptr) {
        c_gfnff_calculator_deallocate(&m_calculator);
        m_calculator.ptr = nullptr;
    }
#endif
    m_initialized = false;
}

Vector GFNFFInterface::BondOrders() const
{
    // External GFN-FF doesn't provide bond orders through C interface
    CurcumaLogger::warn("Bond orders not available from external GFN-FF interface");
    return Vector::Zero(m_natoms * m_natoms);
}

Vector GFNFFInterface::Charges() const
{
    return m_charges;
}

Geometry GFNFFInterface::Gradient() const
{
    return m_gradient;
}

void GFNFFInterface::updateVerbosity()
{
    // Update printlevel based on CurcumaLogger verbosity
    int verbosity = CurcumaLogger::get_verbosity();
    if (verbosity >= 3) {
        m_printlevel = 3;
    } else if (verbosity >= 2) {
        m_printlevel = 2;
    } else if (verbosity >= 1) {
        m_printlevel = 1;
    } else {
        m_printlevel = 0;
    }
}