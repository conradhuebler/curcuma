/*
 * < GFN-FF Computational Method Wrapper >
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
#include "src/core/curcuma_logger.h"

GFNFFComputationalMethod::GFNFFComputationalMethod(const std::string& method_name, const json& config)
    : m_parameters(config)
{
    (void)method_name; // Unused parameter
    m_gfnff = std::make_unique<GFNFF>(config);
}

bool GFNFFComputationalMethod::setMolecule(const Mol& mol) {
    if (!m_gfnff) {
        CurcumaLogger::error("GFNFFComputationalMethod: m_gfnff is nullptr!");
        return false;
    }

    // Cache check is now done inside GFNFF::initializeForceField()
    // where it has access to the ForceField instance

    if (!m_gfnff->InitialiseMolecule(mol)) {
        CurcumaLogger::error("GFNFFComputationalMethod: InitialiseMolecule failed");
        return false;
    }

    return true;
}

bool GFNFFComputationalMethod::updateGeometry(const Matrix& geometry) {
    if (!m_gfnff) {
        CurcumaLogger::error("GFNFFComputationalMethod: m_gfnff is nullptr!");
        return false;
    }

    return m_gfnff->UpdateMolecule(geometry);
}

double GFNFFComputationalMethod::calculateEnergy(bool gradient) {
    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("=== GFNFFComputationalMethod::calculateEnergy() START ===");
    }

    m_last_energy = m_gfnff->Calculation(gradient);

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::success("GFNFFComputationalMethod::calculateEnergy complete");
        CurcumaLogger::energy_abs(m_last_energy, "GFN-FF Energy");
    }

    return m_last_energy;
}

Matrix GFNFFComputationalMethod::getGradient() const {
    return m_gfnff->Gradient();
}

Vector GFNFFComputationalMethod::getCharges() const {
    return m_gfnff->Charges();
}

Vector GFNFFComputationalMethod::getBondOrders() const {
    return m_gfnff->BondOrders();
}

Position GFNFFComputationalMethod::getDipole() const {
    return Position{0.0, 0.0, 0.0};
}

void GFNFFComputationalMethod::setThreadCount(int threads) {
    (void)threads; // Unused parameter - GFN-FF thread management handled internally
}

void GFNFFComputationalMethod::setParameters(const json& params) {
    m_parameters = params;
    if (m_gfnff) {
        m_gfnff->setParameters(params);
    }
}

json GFNFFComputationalMethod::getParameters() const {
    return m_parameters;
}

bool GFNFFComputationalMethod::hasError() const {
    return m_has_error;
}

void GFNFFComputationalMethod::clearError() {
    m_has_error = false;
    m_error_message.clear();
}

json GFNFFComputationalMethod::getEnergyDecomposition() const {
    json energy_json;

    if (!m_gfnff) {
        // Return zero JSON if GFNFF not initialized
        energy_json = {
            {"Bond", 0.0},
            {"Angle", 0.0},
            {"Torsion", 0.0},
            {"Inversion", 0.0},
            {"Dispersion", 0.0},
            {"Coulomb", 0.0},
            {"HBond", 0.0},
            {"XBond", 0.0},
            {"ATM", 0.0},
            {"BATM", 0.0}
        };
        return energy_json;
    }

    // Get all energy components from GFNFF
    energy_json["Bond"] = m_gfnff->BondEnergy();
    energy_json["Angle"] = m_gfnff->AngleEnergy();
    energy_json["Torsion"] = m_gfnff->DihedralEnergy();
    energy_json["Inversion"] = m_gfnff->InversionEnergy();
    energy_json["Dispersion"] = m_gfnff->DispersionEnergy();
    energy_json["Coulomb"] = m_gfnff->CoulombEnergy();
    energy_json["HBond"] = m_gfnff->HydrogenBondEnergy();
    energy_json["XBond"] = m_gfnff->HalogenBondEnergy();
    energy_json["ATM"] = m_gfnff->ATMEnergy();
    energy_json["BATM"] = m_gfnff->BatmEnergy();

    return energy_json;
}

std::string GFNFFComputationalMethod::getErrorMessage() const {
    return m_error_message;
}