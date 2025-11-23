/*
 * < ForceField Method Wrapper Implementation >
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

#include "forcefield_method.h"
#include "src/tools/general.h"
#include "src/core/curcuma_logger.h"
#include "src/core/config_manager.h"

#include <fmt/format.h>
#include <iostream>
#include <thread>

// =================================================================================
// Constructor and Configuration
// =================================================================================

ForceFieldMethod::ForceFieldMethod(const std::string& method_name, const json& config)
    : m_method_name(normalizeMethodName(method_name))
    , m_forcefield(nullptr)
    , m_calculation_done(false)
    , m_last_energy(0.0)
{
    CurcumaLogger::info("Creating ForceFieldMethod");
    CurcumaLogger::param("input_method", method_name);
    CurcumaLogger::param("normalized_method", m_method_name);
    
    if (!validateMethodName(m_method_name)) {
        CurcumaLogger::error("Unsupported ForceField method: " + method_name);
        m_has_error = true;
        m_error_message = fmt::format("Unsupported ForceField method: {}", method_name);
        return;
    }
    
    // Merge with default configuration for this method
    json default_config = getDefaultConfigForMethod(m_method_name);
    m_parameters = MergeJson(default_config, config);
    
    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::param_comparison_table(default_config, config, "ForceField Method Configuration");
    }
    
    // Initialize ForceField
    CurcumaLogger::info("Initializing underlying ForceField engine...");
    if (!initializeForceField()) {
        CurcumaLogger::error("Failed to initialize ForceField engine");
        return;
    }
    
    CurcumaLogger::success("ForceFieldMethod created successfully");
}

bool ForceFieldMethod::initializeForceField() {
    try {
        // Generate ForceField controller from parameters
        json controller = generateForceFieldController();
        
        // Create ForceField instance
        m_forcefield = std::make_unique<ForceField>(controller);
        
        if (!m_forcefield) {
            handleForceFieldError("ForceField creation");
            return false;
        }
        
        // Configure threading
        if (m_parameters.contains("threads")) {
            int threads = m_parameters["threads"];
            setThreadCount(threads);
        }
        
        // Configure parameter caching
        if (m_parameters.contains("parameter_caching")) {
            bool caching = m_parameters["parameter_caching"];
            m_forcefield->setParameterCaching(caching);
        }
        
        clearError();
        return true;
        
    } catch (const std::exception& e) {
        handleForceFieldError(fmt::format("initialization: {}", e.what()));
        return false;
    }
}

json ForceFieldMethod::generateForceFieldController() const {
    json controller;
    
    // Core ForceField parameters
    controller["threads"] = m_parameters.value("threads", 1);
    controller["gradient"] = m_parameters.value("gradient", 1);
    controller["method"] = m_method_name;
    
    // Parameter file handling
    if (m_parameters.contains("param_file")) {
        controller["param_file"] = m_parameters["param_file"];
    } else {
        controller["param_file"] = "none";  // Will trigger auto-generation
    }
    
    // Parameter caching
    if (m_parameters.contains("parameter_caching")) {
        controller["parameter_caching"] = m_parameters["parameter_caching"];
    }
    
    // Geometry file (for parameter auto-naming)
    if (m_parameters.contains("geometry_file")) {
        controller["geometry_file"] = m_parameters["geometry_file"];
    }
    
    // Method-specific parameters
    if (m_method_name == "uff-d3") {
        controller["d3_correction"] = true;
    } else if (m_method_name == "cgfnff") {
        controller["gfnff_mode"] = "native";
    }
    
    return controller;
}

// =================================================================================
// Core ComputationalMethod Interface Implementation
// =================================================================================

bool ForceFieldMethod::setMolecule(const Mol& mol) {
    if (!m_forcefield) {
        handleForceFieldError("ForceField not initialized");
        return false;
    }

    try {
        // Validate molecule for ForceField
        if (!ForceFieldMethodUtils::isMoleculeSuitableForFF(mol)) {
            m_has_error = true;
            m_error_message = "Molecule is not suitable for ForceField calculations";
            return false;
        }

        m_molecule = mol;
        m_calculation_done = false;
        clearError();

        // Set molecule in ForceField
        m_forcefield->setMolecule(mol);

        // Generate parameters if needed
        if (!generateParametersIfNeeded(mol)) {
            CurcumaLogger::error("Failed to generate ForceField parameters");
            handleForceFieldError("parameter generation");
            return false;
        }

        m_initialized = true;

        return true;

    } catch (const std::exception& e) {
        handleForceFieldError(fmt::format("setMolecule: {}", e.what()));
        return false;
    }
}

bool ForceFieldMethod::updateGeometry(const Matrix& geometry) {
    if (!m_initialized || !m_forcefield) {
        handleForceFieldError("ForceField not initialized - call setMolecule first");
        return false;
    }
    
    try {
        m_molecule.m_geometry = geometry;
        m_calculation_done = false;
        clearError();
        
        // Update geometry in ForceField
        m_forcefield->UpdateGeometry(geometry);
        
        return true;
        
    } catch (const std::exception& e) {
        handleForceFieldError(fmt::format("updateGeometry: {}", e.what()));
        return false;
    }
}

double ForceFieldMethod::calculateEnergy(bool gradient)
{
    if (!m_initialized || !m_forcefield) {
        handleForceFieldError("ForceField not initialized - call setMolecule first");
        return 0.0;
    }
    
    try {
        clearError();
        
        // Perform ForceField calculation
        m_last_energy = m_forcefield->Calculate(gradient);

        if (gradient) {
            m_last_gradient = m_forcefield->Gradient();
        }
        
        m_calculation_done = true;
        return m_last_energy;
        
    } catch (const std::exception& e) {
        handleForceFieldError(fmt::format("energy calculation: {}", e.what()));
        return 0.0;
    }
}

// =================================================================================
// Property Access Methods
// =================================================================================

Matrix ForceFieldMethod::getGradient() const {
    if (!m_calculation_done || !m_forcefield) {
        return Matrix::Zero(m_initialized ? m_molecule.m_number_atoms : 1, 3);
    }
    
    try {
        return m_forcefield->Gradient();
    } catch (const std::exception& e) {
        return m_last_gradient;  // Return cached gradient if direct access fails
    }
}

Vector ForceFieldMethod::getCharges() const {
    // ForceField methods typically don't provide charges
    // Return zero vector for consistency
    if (!m_initialized) {
        return Vector::Zero(1);
    }
    return Vector::Zero(m_molecule.m_number_atoms);
}

Vector ForceFieldMethod::getBondOrders() const {
    // ForceField methods don't provide bond orders
    return Vector{};
}

Position ForceFieldMethod::getDipole() const {
    // ForceField methods typically don't provide dipole moments
    return Position{0.0, 0.0, 0.0};
}

// =================================================================================
// Method Information and Configuration
// =================================================================================

void ForceFieldMethod::setThreadCount(int threads) {
    m_thread_count = std::max(1, threads);
    m_parameters["threads"] = m_thread_count;
    
    // ForceField doesn't have a direct setThreadCount method
    // Threading is configured via constructor - would need to recreate FF
    // For now, store for future calculations
    
    if (threads > 1 && m_parameters.value("parameter_caching", true)) {
        // Warn about parameter caching with threading
        fmt::print("Warning: Consider disabling parameter caching for multi-threaded calculations\n");
    }
}

void ForceFieldMethod::setParameters(const json& params) {
    m_parameters = MergeJson(m_parameters, params);
    updateForceFieldConfig();
}

json ForceFieldMethod::getParameters() const {
    // Include current ForceField parameters if available
    json current_params = m_parameters;
    
    if (m_forcefield && m_forcefield->hasParameters()) {
        current_params["ff_parameters"] = m_forcefield->exportCurrentParameters();
    }
    
    return current_params;
}

bool ForceFieldMethod::hasError() const {
    return m_has_error;
}

void ForceFieldMethod::clearError() {
    m_has_error = false;
    m_error_message.clear();
}

std::string ForceFieldMethod::getErrorMessage() const {
    return m_error_message;
}

// =================================================================================
// ForceField-specific Methods
// =================================================================================

void ForceFieldMethod::setParameterCaching(bool enable) {
    m_parameters["parameter_caching"] = enable;
    if (m_forcefield) {
        m_forcefield->setParameterCaching(enable);
    }
}

bool ForceFieldMethod::isParameterCachingEnabled() const {
    return m_parameters.value("parameter_caching", true);
}

bool ForceFieldMethod::loadParametersFromFile(const std::string& filename) {
    if (!m_forcefield) {
        handleForceFieldError("ForceField not initialized");
        return false;
    }
    
    try {
        bool success = m_forcefield->loadParametersFromFile(filename);
        if (success) {
            m_parameters["param_file"] = filename;
            m_calculation_done = false;  // Parameters changed
        }
        return success;
    } catch (const std::exception& e) {
        handleForceFieldError(fmt::format("parameter loading: {}", e.what()));
        return false;
    }
}

bool ForceFieldMethod::saveParametersToFile(const std::string& filename) const {
    if (!m_forcefield) {
        return false;
    }
    
    try {
        return m_forcefield->saveParametersToFile(filename);
    } catch (const std::exception& e) {
        return false;
    }
}

json ForceFieldMethod::exportCurrentParameters() const {
    if (!m_forcefield) {
        return json{};
    }
    
    try {
        return m_forcefield->exportCurrentParameters();
    } catch (const std::exception& e) {
        return json{};
    }
}

bool ForceFieldMethod::hasParameters() const {
    return m_forcefield && m_forcefield->hasParameters();
}

bool ForceFieldMethod::tryLoadAutoParameters(const std::string& method) {
    if (!m_forcefield) {
        return false;
    }
    
    std::string method_to_use = method.empty() ? m_method_name : method;
    return m_forcefield->tryLoadAutoParameters(method_to_use);
}

bool ForceFieldMethod::autoSaveParameters() const {
    if (!m_forcefield) {
        return false;
    }
    
    return m_forcefield->autoSaveParameters();
}

void ForceFieldMethod::printParameterSummary() const {
    if (m_forcefield) {
        m_forcefield->printParameterSummary();
    } else {
        fmt::print("ForceField not initialized - no parameters available\n");
    }
}

Matrix ForceFieldMethod::calculateNumericalGradient() {
    if (!m_forcefield) {
        return Matrix{};
    }
    
    try {
        return m_forcefield->NumGrad();
    } catch (const std::exception& e) {
        handleForceFieldError(fmt::format("numerical gradient: {}", e.what()));
        return Matrix{};
    }
}

void ForceFieldMethod::setParameterFile(const std::string& filename) {
    m_parameters["param_file"] = filename;
    if (m_forcefield) {
        m_forcefield->setParameterFile(filename);
        m_calculation_done = false;  // Parameters might change
    }
}

// =================================================================================
// Static Utility Methods
// =================================================================================

std::vector<std::string> ForceFieldMethod::getSupportedMethods() {
    return {"uff", "uff-d3", "qmdff", "cgfnff"};
}

bool ForceFieldMethod::isMethodSupported(const std::string& method_name) {
    auto supported = getSupportedMethods();
    std::string normalized = method_name;
    std::transform(normalized.begin(), normalized.end(), normalized.begin(), ::tolower);
    
    return std::find(supported.begin(), supported.end(), normalized) != supported.end();
}

// =================================================================================
// Private Helper Methods
// =================================================================================

void ForceFieldMethod::updateForceFieldConfig() {
    // Would need to recreate ForceField for major configuration changes
    // For now, update parameters that can be changed on the fly
    
    if (m_forcefield) {
        if (m_parameters.contains("parameter_caching")) {
            bool caching = m_parameters["parameter_caching"];
            m_forcefield->setParameterCaching(caching);
        }
        
        if (m_parameters.contains("param_file")) {
            std::string param_file = m_parameters["param_file"];
            m_forcefield->setParameterFile(param_file);
        }
    }
}

void ForceFieldMethod::handleForceFieldError(const std::string& operation) {
    m_has_error = true;
    m_error_message = fmt::format("ForceField error during {}", operation);
}

bool ForceFieldMethod::validateMethodName(const std::string& method_name) const {
    return isMethodSupported(method_name);
}

json ForceFieldMethod::getDefaultConfigForMethod(const std::string& method_name) {
    json config = {
        { "threads", ForceFieldMethodUtils::getOptimalThreadCount() },
        { "gradient", true },
        { "parameter_caching", true },
        { "param_file", "none" },
    };

    // Method-specific defaults
    if (method_name == "uff-d3") {
        config["d3_correction"] = true;
    } else if (method_name == "cgfnff") {
        config["gfnff_mode"] = "native";
        config["parameter_caching"] = false;  // cgfnff has parameter issues
    } else if (method_name == "qmdff") {
        config["qmdff_version"] = "latest";
    }
    
    return config;
}

std::string ForceFieldMethod::normalizeMethodName(const std::string& method_name) const {
    std::string normalized = method_name;
    std::transform(normalized.begin(), normalized.end(), normalized.begin(), ::tolower);
    return normalized;
}

bool ForceFieldMethod::generateParametersIfNeeded(const Mol& mol) {
    CurcumaLogger::info("Checking if ForceField parameters need to be generated");
    
    // Check if parameters are already set
    if (m_forcefield->hasParameters()) {
        CurcumaLogger::info("ForceField parameters already available - skipping generation");
        return true;
    }
    
    CurcumaLogger::info("Generating ForceField parameters using ForceFieldGenerator");
    CurcumaLogger::param("method", m_method_name);
    
    try {
        // Create generator configuration
        json generator_config = m_parameters;
        generator_config["method"] = m_method_name;

        // Claude Generated: Create ConfigManager for ForceFieldGenerator (Phase 3B)
        ConfigManager ff_config("forcefield", generator_config);

        // Create ForceFieldGenerator
        ForceFieldGenerator generator(ff_config);
        generator.setMolecule(mol);
        
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info("Calling ForceFieldGenerator::Generate()...");
        }
        
        // Generate parameters
        generator.Generate();
        
        // Get generated parameters
        json generated_params = generator.getParameter();
        
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info("ForceField parameter generation completed");
            CurcumaLogger::param("bonds", static_cast<int>(generated_params.value("bonds", json::array()).size()));
            CurcumaLogger::param("angles", static_cast<int>(generated_params.value("angles", json::array()).size()));
            CurcumaLogger::param("dihedrals", static_cast<int>(generated_params.value("dihedrals", json::array()).size()));
        }
        
        // Apply parameters to ForceField
        m_forcefield->setParameter(generated_params);
        
        CurcumaLogger::success("ForceField parameters generated and applied successfully");
        return true;
        
    } catch (const std::exception& e) {
        CurcumaLogger::error("ForceField parameter generation failed: " + std::string(e.what()));
        return false;
    }
}

// =================================================================================
// Utility Functions
// =================================================================================

namespace ForceFieldMethodUtils {
    
    bool validateForceFieldConfig(const json& config) {
        // Check thread count
        if (config.contains("threads")) {
            int threads = config["threads"];
            if (threads < 1 || threads > 64) {
                return false;
            }
        }
        
        // Check parameter file exists if specified
        if (config.contains("param_file")) {
            std::string param_file = config["param_file"];
            if (param_file != "none" && !std::filesystem::exists(param_file)) {
                return false;
            }
        }
        
        return true;
    }
    
    int getOptimalThreadCount() {
        unsigned int hw_threads = std::thread::hardware_concurrency();
        if (hw_threads == 0) {
            return 1;  // Fallback
        }
        
        // Use 75% of available threads for optimal performance
        return std::max(1u, static_cast<unsigned int>(hw_threads * 0.75));
    }
    
    std::string generateParameterFileName(const std::string& geometry_file) {
        return ForceField::generateParameterFileName(geometry_file);
    }
    
    bool isMoleculeSuitableForFF(const Mol& mol) {
        // Basic validation for ForceField suitability
        if (mol.m_number_atoms < 2) {
            return false;  // Need at least 2 atoms
        }

        if (mol.m_number_atoms > 10000) {
            return false;  // Very large systems might have performance issues
        }

        // Check for reasonable geometry (no NaN or infinite coordinates)
        for (int i = 0; i < mol.m_number_atoms; ++i) {
            for (int j = 0; j < 3; ++j) {
                double coord = mol.m_geometry(i, j);
                if (std::isnan(coord) || std::isinf(coord)) {
                    return false;
                }
            }
        }

        return true;
    }
    
    double estimateMemoryUsage(const Mol& mol, int threads) {
        // Rough estimate of memory usage in MB
        int natoms = mol.m_number_atoms;

        // Base memory for molecule data
        double base_memory = natoms * 0.001;  // ~1 KB per atom

        // Threading overhead
        double thread_memory = threads * 0.1;  // ~100 KB per thread

        // Force field matrices and parameters
        double ff_memory = natoms * natoms * 0.00001;  // Distance matrices etc.

        return base_memory + thread_memory + ff_memory;
    }
}

// Claude Generated: Energy component getter methods (Nov 2025)
double ForceFieldMethod::getBondEnergy() const {
    if (!m_forcefield) return 0.0;
    return m_forcefield->BondEnergy();
}

double ForceFieldMethod::getAngleEnergy() const {
    if (!m_forcefield) return 0.0;
    return m_forcefield->AngleEnergy();
}

double ForceFieldMethod::getDihedralEnergy() const {
    if (!m_forcefield) return 0.0;
    return m_forcefield->DihedralEnergy();
}

double ForceFieldMethod::getInversionEnergy() const {
    if (!m_forcefield) return 0.0;
    return m_forcefield->InversionEnergy();
}

double ForceFieldMethod::getVdWEnergy() const {
    if (!m_forcefield) return 0.0;
    return m_forcefield->VdWEnergy();
}

double ForceFieldMethod::getRepulsionEnergy() const {
    if (!m_forcefield) return 0.0;
    return m_forcefield->RepulsionEnergy();
}

double ForceFieldMethod::getDispersionEnergy() const {
    if (!m_forcefield) return 0.0;
    return m_forcefield->DispersionEnergy();
}

double ForceFieldMethod::getCoulombEnergy() const {
    if (!m_forcefield) return 0.0;
    return m_forcefield->CoulombEnergy();
}