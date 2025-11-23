/*
 * < Native GFN-FF Method Wrapper Implementation >
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

#include "gfnff_method.h"
#include "src/tools/general.h"
#include "src/core/curcuma_logger.h"

#include <fmt/format.h>

// Minimal stub implementation for GFN-FF method
GFNFFMethod::GFNFFMethod(const json& config)
    : m_gfnff(std::make_unique<GFNFF>(config))
    , m_calculation_done(false)
    , m_last_energy(0.0)
{
    m_parameters = config;
}

bool GFNFFMethod::setMolecule(const Mol& mol) {
    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("=== GFNFFMethod::setMolecule() START ===");
        CurcumaLogger::param("atoms", std::to_string(mol.m_number_atoms));
        CurcumaLogger::param("charge", std::to_string(mol.m_charge));
    }

    m_molecule = mol;

    // CRITICAL FIX: Call base class method to set geometry/atoms, then call GFNFF::InitialiseMolecule()
    if (!m_gfnff->QMInterface::InitialiseMolecule(mol)) {
        CurcumaLogger::error("QMInterface::InitialiseMolecule failed");
        return false;
    }

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::success("QMInterface::InitialiseMolecule complete");
        CurcumaLogger::info("Calling GFNFF::InitialiseMolecule() for force field setup...");
    }

    // Now call the actual GFNFF initialization which sets up the force field
    if (!m_gfnff->InitialiseMolecule()) {
        CurcumaLogger::error("GFNFF::InitialiseMolecule failed");
        return false;
    }

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::success("GFNFF::InitialiseMolecule complete");
    }

    m_initialized = true;
    return true;
}

bool GFNFFMethod::updateGeometry(const Matrix& geometry) {
    m_molecule.m_geometry = geometry;
    return m_gfnff->QMInterface::UpdateMolecule(geometry);
}

double GFNFFMethod::calculateEnergy(bool gradient)
{
    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("=== GFNFFMethod::calculateEnergy() START ===");
        CurcumaLogger::param("gradient", gradient ? "true" : "false");
        CurcumaLogger::param("initialized", m_initialized ? "true" : "false");
    }

    if (!m_initialized) {
        CurcumaLogger::error("GFNFFMethod::calculateEnergy: Method not initialized!");
        return 0.0;
    }

    m_last_energy = m_gfnff->Calculation(gradient);
    m_calculation_done = true;

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::success("GFNFFMethod::calculateEnergy complete");
        CurcumaLogger::param("energy_hartree", fmt::format("{:.10f}", m_last_energy));
    }

    return m_last_energy;
}

Matrix GFNFFMethod::getGradient() const {
    return m_gfnff->Gradient();
}

Vector GFNFFMethod::getCharges() const {
    return m_gfnff->Charges();
}

Vector GFNFFMethod::getBondOrders() const {
    return m_gfnff->BondOrders();
}

Position GFNFFMethod::getDipole() const {
    return Position{0.0, 0.0, 0.0};
}

void GFNFFMethod::setThreadCount(int threads) {
    m_thread_count = threads;
}

void GFNFFMethod::setParameters(const json& params) {
    m_parameters = params;
    if (m_gfnff) {
        m_gfnff->setParameters(params);
    }
}

json GFNFFMethod::getParameters() const {
    return m_parameters;
}

bool GFNFFMethod::hasError() const {
    return m_has_error;
}

void GFNFFMethod::clearError() {
    m_has_error = false;
    m_error_message.clear();
}

std::string GFNFFMethod::getErrorMessage() const {
    return m_error_message;
}

// Claude Generated: Energy component getter methods (Nov 2025)
double GFNFFMethod::getBondEnergy() const {
    return m_gfnff ? m_gfnff->BondEnergy() : 0.0;
}

double GFNFFMethod::getAngleEnergy() const {
    return m_gfnff ? m_gfnff->AngleEnergy() : 0.0;
}

double GFNFFMethod::getDihedralEnergy() const {
    return m_gfnff ? m_gfnff->DihedralEnergy() : 0.0;
}

double GFNFFMethod::getInversionEnergy() const {
    return m_gfnff ? m_gfnff->InversionEnergy() : 0.0;
}

double GFNFFMethod::getVdWEnergy() const {
    return m_gfnff ? m_gfnff->VdWEnergy() : 0.0;
}

double GFNFFMethod::getRepulsionEnergy() const {
    return m_gfnff ? m_gfnff->RepulsionEnergy() : 0.0;
}

double GFNFFMethod::getDispersionEnergy() const {
    return m_gfnff ? m_gfnff->DispersionEnergy() : 0.0;
}

double GFNFFMethod::getCoulombEnergy() const {
    return m_gfnff ? m_gfnff->CoulombEnergy() : 0.0;
}