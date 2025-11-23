/*
 * < Unified Energy and Gradient Calculator Implementation >
 * Copyright (C) 2022 - 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
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
 * Claude Generated: Big-Bang refactoring to unified ComputationalMethod interface
 */

#include "energycalculator.h"
#include "energy_calculators/ff_methods/forcefield_method.h"
#include "src/tools/general.h"
#include "src/core/curcuma_logger.h"
#include "config_manager.h"

#include <iostream>
#include <cmath>
#include <fmt/format.h>

// =================================================================================
// Constructors and Initialization
// =================================================================================

// JSON-based constructors (backward compatible) - delegate to ConfigManager versions
EnergyCalculator::EnergyCalculator(const std::string& method, const json& controller)
    : EnergyCalculator(method, ConfigManager("energycalculator", controller))
{
}

EnergyCalculator::EnergyCalculator(const std::string& method, const json& controller, const std::string& basename)
    : EnergyCalculator(method, ConfigManager("energycalculator", controller), basename)
{
}

// ConfigManager-based constructors (new, preferred) - Claude Generated: Phase 3C
EnergyCalculator::EnergyCalculator(const std::string& method, const ConfigManager& config)
    : m_method_name(method)
    , m_basename("")
    , m_energy(0.0)
{
    initializeCommonFromConfig(config);
}

EnergyCalculator::EnergyCalculator(const std::string& method, const ConfigManager& config, const std::string& basename)
    : m_method_name(method)
    , m_basename(basename)
    , m_energy(0.0)
{
    initializeCommonFromConfig(config);
}

EnergyCalculator::~EnergyCalculator() {
    // Smart pointers handle cleanup automatically
    // No need for manual method cleanup
}

// Claude Generated: Phase 3C - ConfigManager-based initialization
void EnergyCalculator::initializeCommonFromConfig(const ConfigManager& config) {
    // Only show initialization info if verbosity is high enough
    if (getEffectiveVerbosity() >= 2) {
        CurcumaLogger::info("Initializing EnergyCalculator (ConfigManager)");
        CurcumaLogger::param("method", m_method_name);
    }

    if (getEffectiveVerbosity() >= 3) {
        CurcumaLogger::info("=== EnergyCalculator Initialization Debug (ConfigManager) ===");
        CurcumaLogger::param("method_name", m_method_name);
        CurcumaLogger::param("basename", m_basename);
        CurcumaLogger::param("verbosity_override", m_verbosity_override);
    }

    // Store configuration
    m_controller = config.exportConfig();  // Export for compatibility with existing code

    if (getEffectiveVerbosity() >= 3) {
        CurcumaLogger::param_table(m_controller, "EnergyCalculator Configuration");
    }

    // Extract common parameters from ConfigManager
    if (m_controller.contains("param_file")) {
        m_parameter["param_file"] = m_controller["param_file"];
        if (getEffectiveVerbosity() >= 2) {
            CurcumaLogger::param("param_file", m_controller["param_file"].get<std::string>());
        }
    }

    if (m_controller.contains("geometry_file")) {
        m_geometry_file = m_controller["geometry_file"];
        if (getEffectiveVerbosity() >= 2) {
            CurcumaLogger::param("geometry_file", m_geometry_file);
        }
    }

    if (!m_basename.empty()) {
        m_geometry_file = m_basename + ".xyz";
        if (getEffectiveVerbosity() >= 2) {
            CurcumaLogger::param("auto_geometry_file", m_geometry_file);
        }
    }

    // Extract multiplicity and other settings
    m_mult = m_controller.value("multi", 1);
    if (getEffectiveVerbosity() >= 2) {
        CurcumaLogger::param("multiplicity", m_mult);
    }

    if (getEffectiveVerbosity() >= 3) {
        CurcumaLogger::info("Creating computational method using MethodFactory...");
    }

    // Create computational method using factory
    if (!createMethod(m_method_name, m_controller)) {
        m_error = true;
        m_error_message = fmt::format("Failed to create method: {}", m_method_name);
        CurcumaLogger::error("Failed to create computational method: " + m_method_name);
        return;
    }

    CurcumaLogger::success("EnergyCalculator initialized successfully");
    CurcumaLogger::param("actual_method", m_method_name);
}

// Backward compatible wrapper - delegates to ConfigManager version
void EnergyCalculator::initializeCommon(const json& controller) {
    initializeCommonFromConfig(ConfigManager("energycalculator", controller));
}

bool EnergyCalculator::createMethod(const std::string& method_name, const json& config) {
    CurcumaLogger::info("Creating computational method");
    CurcumaLogger::param("requested_method", method_name);
    
    try {
        // Prepare configuration for method creation
        json method_config = config;
        method_config["method"] = method_name;
        method_config["multiplicity"] = m_mult;
        
        // Add geometry file for parameter caching
        if (!m_geometry_file.empty()) {
            method_config["geometry_file"] = m_geometry_file;
        }
        if (!m_basename.empty()) {
            method_config["basename"] = m_basename;
        }
        
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info("Method configuration prepared:");
            CurcumaLogger::param_table(method_config, "Method Configuration");
        }
        
        // Use MethodFactory to create appropriate method
        CurcumaLogger::info("Calling MethodFactory::create...");
        m_method = MethodFactory::create(method_name, method_config);
        
        if (!m_method) {
            handleMethodError("method creation via factory");
            return false;
        }
        
        // Store actual method name (may differ from input due to resolution)
        m_method_name = m_method->getMethodName();
        
        // Configure threading if specified
        if (config.contains("threads")) {
            int threads = config["threads"];
            m_method->setThreadCount(threads);
        }
        
        ClearError();
        return true;
        
    } catch (const MethodCreationException& e) {
        m_error = true;
        m_error_message = fmt::format("Method creation failed: {}", e.what());
        return false;
    } catch (const std::exception& e) {
        handleMethodError(fmt::format("method creation: {}", e.what()));
        return false;
    }
}

// =================================================================================
// Molecular Setup and Geometry Management
// =================================================================================

void EnergyCalculator::setMolecule(const Mol& mol) {
    CurcumaLogger::info("Setting molecule in EnergyCalculator");
    CurcumaLogger::param("atoms", mol.m_number_atoms);
    CurcumaLogger::param("charge", mol.m_charge);
    CurcumaLogger::param("spin", mol.m_spin);
    
    if (!m_method) {
        CurcumaLogger::error("No computational method available for setMolecule");
        handleMethodError("no method available for setMolecule");
        return;
    }
    
    try {
        m_mol = mol;
        m_atoms = mol.m_number_atoms;
        
        // Initialize gradient matrix
        m_gradient = Matrix::Zero(m_atoms, 3);
        
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info("Calling computational method setMolecule...");
        }
        
        // Set molecule in computational method
        if (!m_method->setMolecule(mol)) {
            CurcumaLogger::error("Failed to set molecule in computational method");
            handleMethodError("setMolecule in computational method");
            return;
        }
        
        m_initialized = true;
        ClearError();

        if (CurcumaLogger::get_verbosity() >= 1) {
            CurcumaLogger::success(fmt::format("Molecule set: {} atoms, method: {}", m_atoms, m_method_name));
        }

    } catch (const std::exception& e) {
        handleMethodError(fmt::format("setMolecule: {}", e.what()));
    }
}

void EnergyCalculator::updateGeometry(const double* coord) {
    if (!m_initialized) {
        handleMethodError("updateGeometry called before setMolecule");
        return;
    }
    
    try {
        // Convert coordinates to Matrix format
        Matrix geometry = Matrix::Zero(m_atoms, 3);
        convertCoordinates(coord, geometry);
        
        // Update in computational method
        if (!m_method->updateGeometry(geometry)) {
            handleMethodError("updateGeometry in computational method");
            return;
        }
        
        // Store updated geometry
        m_mol.m_geometry = geometry;
        ClearError();
        
    } catch (const std::exception& e) {
        handleMethodError(fmt::format("updateGeometry: {}", e.what()));
    }
}

void EnergyCalculator::updateGeometry(const std::vector<double>& geometry) {
    if (!m_initialized) {
        handleMethodError("updateGeometry called before setMolecule");
        return;
    }
    
    try {
        // Convert coordinates to Matrix format
        Matrix coord_matrix = Matrix::Zero(m_atoms, 3);
        convertCoordinates(geometry, coord_matrix);
        
        // Update in computational method
        if (!m_method->updateGeometry(coord_matrix)) {
            handleMethodError("updateGeometry in computational method");
            return;
        }
        
        // Store updated geometry
        m_mol.m_geometry = coord_matrix;
        ClearError();
        
    } catch (const std::exception& e) {
        handleMethodError(fmt::format("updateGeometry: {}", e.what()));
    }
}

void EnergyCalculator::updateGeometry(const Matrix& geometry) {
    if (!m_initialized) {
        handleMethodError("updateGeometry called before setMolecule");
        return;
    }
    
    try {
        // Update in computational method
        if (!m_method->updateGeometry(geometry)) {
            handleMethodError("updateGeometry in computational method");
            return;
        }
        
        // Store updated geometry
        m_mol.m_geometry = geometry;
        ClearError();
        
    } catch (const std::exception& e) {
        handleMethodError(fmt::format("updateGeometry: {}", e.what()));
    }
}

void EnergyCalculator::updateGeometry(const Eigen::VectorXd& geometry) {
    if (!m_initialized) {
        handleMethodError("updateGeometry called before setMolecule");
        return;
    }
    
    try {
        // Convert coordinates to Matrix format
        Matrix coord_matrix = Matrix::Zero(m_atoms, 3);
        convertCoordinates(geometry, coord_matrix);
        
        // Update in computational method
        if (!m_method->updateGeometry(coord_matrix)) {
            handleMethodError("updateGeometry in computational method");
            return;
        }
        
        // Store updated geometry
        m_mol.m_geometry = coord_matrix;
        ClearError();
        
    } catch (const std::exception& e) {
        handleMethodError(fmt::format("updateGeometry: {}", e.what()));
    }
}

// =================================================================================
// Energy and Gradient Calculations
// =================================================================================

double EnergyCalculator::CalculateEnergy(bool gradient)
{
    if (getEffectiveVerbosity() >= 2) {
        CurcumaLogger::info("Starting energy calculation");
        CurcumaLogger::param("gradient", gradient);
        CurcumaLogger::param("method", m_method ? m_method->getMethodName() : "null");
    }

    if (!m_initialized || !m_method) {
        CurcumaLogger::error("CalculateEnergy called before proper initialization");
        CurcumaLogger::param("initialized", m_initialized);
        CurcumaLogger::param("method_available", m_method != nullptr);
        handleMethodError("CalculateEnergy called before proper initialization");
        return 0.0;
    }
    
    try {
        ClearError();

        if (getEffectiveVerbosity() >= 3) {
            CurcumaLogger::info("Calling computational method calculateEnergy...");
        }

        // Temporarily override global verbosity for submodule calculations
        int original_verbosity = CurcumaLogger::get_verbosity();
        if (m_verbosity_override >= 0) {
            CurcumaLogger::set_verbosity(m_verbosity_override);
        }

        // Perform calculation using computational method
        m_energy = m_method->calculateEnergy(gradient);

        // Restore original verbosity
        CurcumaLogger::set_verbosity(original_verbosity);

        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info("Energy calculation completed");
        }

        // Use effective verbosity to respect submodule silencing
        if (getEffectiveVerbosity() >= 1) {
            CurcumaLogger::energy_abs(m_energy, "Total Energy");
        }

        // Check for method-specific errors
        if (m_method->hasError()) {
            CurcumaLogger::error("Method reported error after calculation");
            CurcumaLogger::error("Method error: " + m_method->getErrorMessage());
            m_error = true;
            m_error_message = fmt::format("Method error: {}", m_method->getErrorMessage());
            return 0.0;
        }
        
        // Get gradient if requested
        if (gradient) {
            m_gradient = m_method->getGradient();
        }
        
        // Check for NaN values
        if (checkForNaN(m_energy, gradient ? m_gradient : Matrix{})) {
            m_containsNaN = true;
            handleMethodError("NaN values detected in calculation results");
            return 0.0;
        }

        if (getEffectiveVerbosity() >= 1) {
            CurcumaLogger::energy_abs(m_energy, fmt::format("{} Final Energy", m_method_name));
        }

        return m_energy;
        
    } catch (const std::exception& e) {
        handleMethodError(fmt::format("CalculateEnergy: {}", e.what()));
        return 0.0;
    }
}

Matrix EnergyCalculator::Gradient() const {
    if (!m_method) {
        return Matrix::Zero(m_atoms > 0 ? m_atoms : 1, 3);
    }
    
    try {
        return m_method->getGradient();
    } catch (const std::exception& e) {
        // Return cached gradient if direct access fails
        return m_gradient;
    }
}

Eigen::MatrixXd EnergyCalculator::NumGrad() {
    if (!m_initialized || !m_method) {
        return Eigen::MatrixXd::Zero(1, 3);
    }
    
    try {
        // Numerical gradient calculation using finite differences
        Eigen::MatrixXd gradient = Eigen::MatrixXd::Zero(m_atoms, 3);
        double dx = 1e-4;
        double E1, E2;
        
        Matrix current_geometry = m_mol.m_geometry;
        
        for (int i = 0; i < m_atoms; ++i) {
            for (int j = 0; j < 3; ++j) {
                // Forward step
                Matrix geom_forward = current_geometry;
                geom_forward(i, j) += dx;
                updateGeometry(geom_forward);
                E1 = CalculateEnergy(false);

                // Backward step
                Matrix geom_backward = current_geometry;
                geom_backward(i, j) -= dx;
                updateGeometry(geom_backward);
                E2 = CalculateEnergy(false);

                // Central difference
                gradient(i, j) = (E1 - E2) / (2 * dx);
            }
        }
        
        // Restore original geometry
        updateGeometry(current_geometry);
        
        return gradient;
        
    } catch (const std::exception& e) {
        handleMethodError(fmt::format("numerical gradient: {}", e.what()));
        return Eigen::MatrixXd::Zero(m_atoms, 3);
    }
}

// =================================================================================
// Property Access Methods
// =================================================================================

Vector EnergyCalculator::Charges() const {
    if (!m_method) {
        return Vector::Zero(m_atoms);
    }
    
    try {
        return m_method->getCharges();
    } catch (const std::exception& e) {
        return Vector::Zero(m_atoms);
    }
}

Position EnergyCalculator::Dipole() const {
    if (!m_method) {
        return Position{0.0, 0.0, 0.0};
    }
    
    try {
        return m_method->getDipole();
    } catch (const std::exception& e) {
        return Position{0.0, 0.0, 0.0};
    }
}

std::vector<std::vector<double>> EnergyCalculator::BondOrders() const {
    if (!m_method) {
        return std::vector<std::vector<double>>{{}};
    }
    
    try {
        Vector bond_orders = m_method->getBondOrders();
        // Convert Vector to nested vector format for compatibility
        std::vector<std::vector<double>> result;
        if (bond_orders.size() > 0) {
            result.push_back(std::vector<double>(bond_orders.data(), 
                                               bond_orders.data() + bond_orders.size()));
        }
        return result;
    } catch (const std::exception& e) {
        return std::vector<std::vector<double>>{{}};
    }
}

Vector EnergyCalculator::Energies() const {
    if (!m_method) {
        return Vector{};
    }
    
    try {
        return m_method->getOrbitalEnergies();
    } catch (const std::exception& e) {
        return Vector{};
    }
}

Vector EnergyCalculator::OrbitalOccuptations() const {
    if (!m_method) {
        return Vector{};
    }
    
    try {
        return m_method->getOrbitalOccupations();
    } catch (const std::exception& e) {
        return Vector{};
    }
}

int EnergyCalculator::NumElectrons() const {
    if (!m_method) {
        return 0;
    }
    
    try {
        return m_method->getNumElectrons();
    } catch (const std::exception& e) {
        return 0;
    }
}

// =================================================================================
// Error Handling and Status
// =================================================================================

std::string EnergyCalculator::ErrorMessage() const {
    if (m_method && m_method->hasError()) {
        return fmt::format("EnergyCalculator: {} | Method: {}", 
                          m_error_message, m_method->getErrorMessage());
    }
    return m_error_message;
}

void EnergyCalculator::ClearError() {
    m_error = false;
    m_error_message.clear();
    m_containsNaN = false;
    
    if (m_method) {
        m_method->clearError();
    }
}

// =================================================================================
// Configuration and Parameters
// =================================================================================

void EnergyCalculator::setGeometryFile(const std::string& filename) {
    m_geometry_file = filename;
    if (m_method) {
        json params = m_method->getParameters();
        params["geometry_file"] = filename;
        m_method->setParameters(params);
    }
}

void EnergyCalculator::setBasename(const std::string& basename) {
    m_basename = basename;
    m_geometry_file = basename + ".xyz";
}

void EnergyCalculator::setThreadCount(int threads) {
    if (m_method) {
        m_method->setThreadCount(threads);
    }
}

std::string EnergyCalculator::getMethodName() const {
    if (m_method) {
        return m_method->getMethodName();
    }
    return m_method_name;
}

bool EnergyCalculator::isThreadSafe() const {
    if (m_method) {
        return m_method->isThreadSafe();
    }
    return false;
}

// =================================================================================
// Static Methods
// =================================================================================

void EnergyCalculator::PrintAvailableMethods() {
    MethodFactory::printAvailableMethods();
}

// =================================================================================
// Advanced Features
// =================================================================================

bool EnergyCalculator::SaveToFile(const std::string& filename) const {
    if (!m_method) {
        return false;
    }
    
    try {
        return m_method->saveToFile(filename);
    } catch (const std::exception& e) {
        return false;
    }
}

json EnergyCalculator::GetMethodInfo() const {
    json info;
    info["method_name"] = m_method_name;
    info["initialized"] = m_initialized;
    info["has_error"] = m_error;
    info["error_message"] = m_error_message;
    info["atom_count"] = m_atoms;
    info["last_energy"] = m_energy;
    info["contains_nan"] = m_containsNaN;
    
    if (m_method) {
        info["thread_safe"] = m_method->isThreadSafe();
        info["has_gradient"] = m_method->hasGradient();
        info["method_parameters"] = m_method->getParameters();
    }
    
    return info;
}

// =================================================================================
// Submodule Verbosity Control (Claude Generated)
// =================================================================================

void EnergyCalculator::setVerbosity(int level)
{
    m_verbosity_override = level;
}

void EnergyCalculator::resetVerbosity()
{
    m_verbosity_override = -1;
}

int EnergyCalculator::getEffectiveVerbosity() const
{
    return (m_verbosity_override >= 0) ? m_verbosity_override : CurcumaLogger::get_verbosity();
}

// =================================================================================
// Internal Helper Methods
// =================================================================================

void EnergyCalculator::handleMethodError(const std::string& operation) {
    m_error = true;
    m_error_message = fmt::format("EnergyCalculator error during {}", operation);
}

bool EnergyCalculator::checkForNaN(double energy, const Matrix& gradient) {
    // Check energy for NaN
    if (std::isnan(energy) || std::isinf(energy)) {
        return true;
    }
    
    // Check gradient for NaN if provided
    if (gradient.size() > 0) {
        for (int i = 0; i < gradient.rows(); ++i) {
            for (int j = 0; j < gradient.cols(); ++j) {
                if (std::isnan(gradient(i, j)) || std::isinf(gradient(i, j))) {
                    return true;
                }
            }
        }
    }
    
    return false;
}

void EnergyCalculator::convertCoordinates(const double* coord, Matrix& geometry) {
    for (int i = 0; i < m_atoms; ++i) {
        geometry(i, 0) = coord[3 * i + 0];
        geometry(i, 1) = coord[3 * i + 1];
        geometry(i, 2) = coord[3 * i + 2];
    }
}

void EnergyCalculator::convertCoordinates(const std::vector<double>& coord, Matrix& geometry) {
    for (int i = 0; i < m_atoms; ++i) {
        geometry(i, 0) = coord[3 * i + 0];
        geometry(i, 1) = coord[3 * i + 1];
        geometry(i, 2) = coord[3 * i + 2];
    }
}

void EnergyCalculator::convertCoordinates(const Eigen::VectorXd& coord, Matrix& geometry) {
    for (int i = 0; i < m_atoms; ++i) {
        geometry(i, 0) = coord[3 * i + 0];
        geometry(i, 1) = coord[3 * i + 1];
        geometry(i, 2) = coord[3 * i + 2];
    }
}

// Claude Generated: Energy component getter methods for regression testing (Nov 2025)
// These delegate to the ComputationalMethod interface (which now has virtual methods)
double EnergyCalculator::getBondEnergy() const {
    return m_method ? m_method->getBondEnergy() : 0.0;
}

double EnergyCalculator::getAngleEnergy() const {
    return m_method ? m_method->getAngleEnergy() : 0.0;
}

double EnergyCalculator::getDihedralEnergy() const {
    return m_method ? m_method->getDihedralEnergy() : 0.0;
}

double EnergyCalculator::getInversionEnergy() const {
    return m_method ? m_method->getInversionEnergy() : 0.0;
}

double EnergyCalculator::getVdWEnergy() const {
    return m_method ? m_method->getVdWEnergy() : 0.0;
}

double EnergyCalculator::getRepulsionEnergy() const {
    return m_method ? m_method->getRepulsionEnergy() : 0.0;
}

double EnergyCalculator::getDispersionEnergy() const {
    return m_method ? m_method->getDispersionEnergy() : 0.0;
}

double EnergyCalculator::getCoulombEnergy() const {
    return m_method ? m_method->getCoulombEnergy() : 0.0;
}

double EnergyCalculator::getNonBondedEnergy() const {
    return getVdWEnergy() + getRepulsionEnergy();
}