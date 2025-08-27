/*
 * < External GFN-FF Method Wrapper >
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
 * Claude Generated: ComputationalMethod wrapper for external GFN-FF interface
 */

#include "external_gfnff_method.h"
#include "src/core/curcuma_logger.h"

#include <fmt/format.h>

ExternalGFNFFMethod::ExternalGFNFFMethod(const json& config)
    : m_config(config)
    , m_initialized(false)
{
    // Create the underlying GFN-FF interface
    m_interface = std::make_unique<GFNFFInterface>(config);

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("ExternalGFNFFMethod wrapper initialized");
        CurcumaLogger::param("method_wrapper", "External GFN-FF");
    }
}

bool ExternalGFNFFMethod::setMolecule(const Mol& mol)
{
    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("ExternalGFNFFMethod::setMolecule called");
        CurcumaLogger::param("natoms", mol.m_number_atoms);
        CurcumaLogger::param("charge", mol.m_charge);
    }

    m_current_molecule = mol;

    // Initialize the underlying interface
    bool success = m_interface->InitialiseMolecule(mol);
    m_initialized = success;

    if (!success) {
        CurcumaLogger::error("Failed to initialize external GFN-FF interface with molecule");
        return false;
    }

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::success("External GFN-FF method ready for calculations");
    }

    return true;
}

bool ExternalGFNFFMethod::updateGeometry(const Matrix& geometry)
{
    if (!m_initialized) {
        CurcumaLogger::error("External GFN-FF method not initialized - call setMolecule() first");
        return false;
    }

    // Update the molecule geometry and reinitialize interface
    m_current_molecule.m_geometry = geometry;
    return m_interface->InitialiseMolecule(m_current_molecule);
}

double ExternalGFNFFMethod::calculateEnergy(bool gradient)
{
    if (!m_initialized) {
        CurcumaLogger::error("External GFN-FF method not initialized - call setMolecule() first");
        return 0.0;
    }

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("ExternalGFNFFMethod::calculateEnergy called");
        CurcumaLogger::param("gradient_requested", gradient);
    }

    // Delegate to the underlying interface
    double energy = m_interface->Calculation(gradient);

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::param("final_energy", fmt::format("{:.8f} Eh", energy));
    }

    return energy;
}

Matrix ExternalGFNFFMethod::getGradient() const
{
    if (!m_initialized) {
        CurcumaLogger::warn("External GFN-FF method not initialized - returning empty gradient");
        return Matrix();
    }

    return m_interface->Gradient();
}

Vector ExternalGFNFFMethod::getCharges() const
{
    if (!m_initialized) {
        CurcumaLogger::warn("External GFN-FF method not initialized - returning empty charges");
        return Vector();
    }

    return m_interface->Charges();
}

Vector ExternalGFNFFMethod::getBondOrders() const
{
    if (!m_initialized) {
        CurcumaLogger::warn("External GFN-FF method not initialized - returning empty bond orders");
        return Vector();
    }

    return m_interface->BondOrders();
}

Position ExternalGFNFFMethod::getDipole() const
{
    if (!m_initialized) {
        CurcumaLogger::warn("External GFN-FF method not initialized - returning zero dipole");
        return Position::Zero();
    }

    // External GFN-FF doesn't provide dipole through C interface
    CurcumaLogger::warn("Dipole moment not available from external GFN-FF interface");
    return Position::Zero();
}

void ExternalGFNFFMethod::setParameters(const json& params)
{
    m_config = params;
    if (m_interface) {
        // Reinitialize interface with new parameters
        m_interface = std::make_unique<GFNFFInterface>(m_config);
        if (m_initialized) {
            m_interface->InitialiseMolecule(m_current_molecule);
        }
    }
}

json ExternalGFNFFMethod::getParameters() const
{
    json info = m_config;
    info["method_type"] = "External GFN-FF";
    info["initialized"] = m_initialized;
    info["supports_gradients"] = true;
    info["thread_safe"] = false;

    if (m_initialized) {
        info["natoms"] = m_current_molecule.m_number_atoms;
        info["charge"] = m_current_molecule.m_charge;
    }

    return info;
}