/*
 * < EHT Method Wrapper Implementation >
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

#include "eht_method.h"
#include "src/tools/general.h"

#include <fmt/format.h>
#include <iostream>

// =================================================================================
// Constructor and Default Configuration
// =================================================================================

EHTMethod::EHTMethod(const json& config)
    : m_eht(nullptr)
    , m_calculation_done(false)
    , m_last_energy(0.0)
{
    m_parameters = MergeJson(getDefaultConfig(), config);
    updateEHTParameters();
    
    // Initialize EHT with Wolfsberg-Helmholz constant from parameters
    double K = m_parameters.value("wolfsberg_helmholz_constant", 1.75);
    m_eht = std::make_unique<EHT>(K);
}

json EHTMethod::getDefaultConfig() {
    return json{
        {"wolfsberg_helmholz_constant", 1.75},  // Standard WH constant
        {"verbose", false},                     // Enable detailed output
        {"threads", 1},                         // Single-threaded by default
        {"print_orbitals", false},             // Print orbital analysis
        {"num_orbitals_print", 5},             // Number of orbitals around gap
        {"save_orbitals", false},              // Save orbital data to file
        {"orbital_threshold", 1e-6}            // Orbital coefficient threshold
    };
}

// =================================================================================
// Core ComputationalMethod Interface Implementation
// =================================================================================

bool EHTMethod::setMolecule(const Mol& mol) {
    try {
        // Validate molecule for EHT compatibility
        if (!EHTMethodUtils::isMoleculeSupported(mol)) {
            m_has_error = true;
            m_error_message = "Molecule contains unsupported elements for EHT";
            return false;
        }
        
        m_molecule = mol;
        m_calculation_done = false;
        m_initialized = true;
        clearError();
        
        // Initialize EHT with molecule (use QMInterface base class method)
        if (!m_eht->QMInterface::InitialiseMolecule(mol)) {
            handleEHTError("molecule initialization");
            return false;
        }
        
        return true;
        
    } catch (const std::exception& e) {
        m_has_error = true;
        m_error_message = fmt::format("EHT setMolecule failed: {}", e.what());
        return false;
    }
}

bool EHTMethod::updateGeometry(const Matrix& geometry) {
    if (!m_initialized) {
        m_has_error = true;
        m_error_message = "Method not initialized - call setMolecule first";
        return false;
    }
    
    try {
        m_molecule.m_geometry = geometry;
        m_calculation_done = false;
        clearError();
        
        // Update geometry in EHT method
        if (!m_eht->UpdateMolecule(geometry)) {
            handleEHTError("geometry update");
            return false;
        }
        
        return true;
        
    } catch (const std::exception& e) {
        m_has_error = true;
        m_error_message = fmt::format("EHT updateGeometry failed: {}", e.what());
        return false;
    }
}

double EHTMethod::calculateEnergy(bool gradient, bool verbose) {
    if (!m_initialized) {
        m_has_error = true;
        m_error_message = "Method not initialized - call setMolecule first";
        return 0.0;
    }
    
    try {
        clearError();
        
        // Warn about gradient request (EHT doesn't have analytical gradients)
        if (gradient && verbose) {
            fmt::print("Warning: EHT does not provide analytical gradients\n");
        }
        
        // Perform EHT calculation
        m_last_energy = m_eht->Calculation(false, verbose);  // EHT doesn't support gradients
        m_calculation_done = true;
        
        // Print orbital analysis if requested
        if (m_parameters.value("print_orbitals", false)) {
            int num_orbitals = m_parameters.value("num_orbitals_print", 5);
            m_eht->printOrbitalAnalysis(num_orbitals);
        }
        
        return m_last_energy;
        
    } catch (const std::exception& e) {
        handleEHTError(fmt::format("energy calculation: {}", e.what()));
        return 0.0;
    }
}

// =================================================================================
// Property Access Methods
// =================================================================================

Matrix EHTMethod::getGradient() const {
    // EHT doesn't provide analytical gradients - return zero matrix
    if (!m_initialized) {
        return Matrix::Zero(1, 3);
    }
    return Matrix::Zero(m_molecule.m_number_atoms, 3);
}

Vector EHTMethod::getCharges() const {
    // EHT doesn't calculate charges - return zero vector
    if (!m_initialized) {
        return Vector::Zero(1);
    }
    return Vector::Zero(m_molecule.m_number_atoms);
}

Vector EHTMethod::getBondOrders() const {
    // EHT doesn't calculate bond orders - return empty vector
    return Vector{};
}

Position EHTMethod::getDipole() const {
    // EHT doesn't calculate dipole - return zero vector
    return Position{0.0, 0.0, 0.0};
}

// =================================================================================
// Method Information and Configuration
// =================================================================================

void EHTMethod::setThreadCount(int threads) {
    m_thread_count = std::max(1, threads);
    // EHT can potentially be parallelized at the matrix diagonalization level
    // For now, store thread count for future use
}

void EHTMethod::setParameters(const json& params) {
    m_parameters = MergeJson(m_parameters, params);
    updateEHTParameters();
}

json EHTMethod::getParameters() const {
    return m_parameters;
}

bool EHTMethod::hasError() const {
    return m_has_error;
}

void EHTMethod::clearError() {
    m_has_error = false;
    m_error_message.clear();
}

std::string EHTMethod::getErrorMessage() const {
    return m_error_message;
}

// =================================================================================
// EHT-specific Methods
// =================================================================================

Vector EHTMethod::getOrbitalEnergies() const {
    if (!m_calculation_done || !m_eht) {
        return Vector{};
    }
    
    try {
        return m_eht->Energies();
    } catch (const std::exception& e) {
        return Vector{};
    }
}

Matrix EHTMethod::getMolecularOrbitals() const {
    if (!m_calculation_done || !m_eht) {
        return Matrix{};
    }
    
    try {
        return m_eht->MolecularOrbitals();
    } catch (const std::exception& e) {
        return Matrix{};
    }
}

int EHTMethod::getNumElectrons() const {
    if (!m_calculation_done || !m_eht) {
        return 0;
    }
    
    try {
        return m_eht->NumElectrons();
    } catch (const std::exception& e) {
        return 0;
    }
}

double EHTMethod::getHOMOLUMOGap() const {
    if (!m_calculation_done || !m_eht) {
        return 0.0;
    }
    
    try {
        return m_eht->getHOMOLUMOGap();
    } catch (const std::exception& e) {
        return 0.0;
    }
}

double EHTMethod::getHOMOEnergy() const {
    if (!m_calculation_done || !m_eht) {
        return 0.0;
    }
    
    try {
        return m_eht->getHOMOEnergy();
    } catch (const std::exception& e) {
        return 0.0;
    }
}

double EHTMethod::getLUMOEnergy() const {
    if (!m_calculation_done || !m_eht) {
        return 0.0;
    }
    
    try {
        return m_eht->getLUMOEnergy();
    } catch (const std::exception& e) {
        return 0.0;
    }
}

void EHTMethod::setWolfsbergHelmholzConstant(double K) {
    if (m_eht) {
        m_eht->setWolfsbergHelmholzConstant(K);
        m_parameters["wolfsberg_helmholz_constant"] = K;
        m_calculation_done = false;  // Need to recalculate
    }
}

double EHTMethod::getWolfsbergHelmholzConstant() const {
    if (m_eht) {
        return m_eht->getWolfsbergHelmholzConstant();
    }
    return m_parameters.value("wolfsberg_helmholz_constant", 1.75);
}

void EHTMethod::printOrbitalAnalysis(int num_orbitals_around_gap) const {
    if (m_calculation_done && m_eht) {
        m_eht->printOrbitalAnalysis(num_orbitals_around_gap);
    } else {
        fmt::print("No EHT calculation results available for orbital analysis\n");
    }
}

bool EHTMethod::saveToFile(const std::string& filename) const {
    if (!m_calculation_done) {
        return false;
    }
    
    try {
        // Save EHT-specific data (orbital energies, coefficients, etc.)
        json output_data;
        output_data["method"] = "eht";
        output_data["energy"] = m_last_energy;
        output_data["orbital_energies"] = getOrbitalEnergies();
        output_data["num_electrons"] = getNumElectrons();
        output_data["homo_energy"] = getHOMOEnergy();
        output_data["lumo_energy"] = getLUMOEnergy();
        output_data["homo_lumo_gap"] = getHOMOLUMOGap();
        output_data["wolfsberg_helmholz_constant"] = getWolfsbergHelmholzConstant();
        output_data["parameters"] = m_parameters;
        
        std::ofstream file(filename);
        if (file.is_open()) {
            file << std::setw(2) << output_data << std::endl;
            file.close();
            return true;
        }
        
    } catch (const std::exception& e) {
        return false;
    }
    
    return false;
}

// =================================================================================
// Private Helper Methods
// =================================================================================

bool EHTMethod::initializeEHT() {
    if (!m_eht) {
        double K = m_parameters.value("wolfsberg_helmholz_constant", 1.75);
        m_eht = std::make_unique<EHT>(K);
    }
    return m_eht != nullptr;
}

void EHTMethod::updateEHTParameters() {
    if (m_eht && m_parameters.contains("wolfsberg_helmholz_constant")) {
        double K = m_parameters["wolfsberg_helmholz_constant"];
        m_eht->setWolfsbergHelmholzConstant(K);
        m_calculation_done = false;  // Parameters changed, need recalculation
    }
}

void EHTMethod::handleEHTError(const std::string& operation) {
    m_has_error = true;
    m_error_message = fmt::format("EHT error during {}", operation);
}

// =================================================================================
// Utility Functions
// =================================================================================

namespace EHTMethodUtils {
    
    bool validateEHTConfig(const json& config) {
        // Check for required parameters and validate ranges
        if (config.contains("wolfsberg_helmholz_constant")) {
            double K = config["wolfsberg_helmholz_constant"];
            if (K <= 0.0 || K > 3.0) {
                return false;  // Invalid WH constant
            }
        }
        
        if (config.contains("threads")) {
            int threads = config["threads"];
            if (threads < 1 || threads > 32) {
                return false;  // Invalid thread count
            }
        }
        
        return true;
    }
    
    std::vector<int> getSupportedElements() {
        // EHT supports common organic elements
        return {1, 6, 7, 8, 9, 16, 17};  // H, C, N, O, F, S, Cl
    }
    
    bool isMoleculeSupported(const Mol& mol) {
        auto supported = getSupportedElements();
        
        for (int atom : mol.m_atoms) {
            if (std::find(supported.begin(), supported.end(), atom) == supported.end()) {
                return false;  // Unsupported element found
            }
        }
        
        return true;
    }
}