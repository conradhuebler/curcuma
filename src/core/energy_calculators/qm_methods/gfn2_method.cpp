/*
 * < GFN2 Method Wrapper Implementation >
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Based on the GFN2-xTB method developed by:
 *   Stefan Grimme, Christoph Bannwarth, Sebastian Ehlert
 *   Mulliken Center for Theoretical Chemistry, University of Bonn
 *
 * This program is free software under GPL-3.0
 */

#include "gfn2_method.h"
#include "src/tools/general.h"
#include "src/core/curcuma_logger.h"

#include <fmt/format.h>
#include <fstream>
#include <iomanip>

// =================================================================================
// Constructor and Default Configuration
// =================================================================================

GFN2Method::GFN2Method(const json& config)
    : m_gfn2(nullptr)
    , m_calculation_done(false)
    , m_last_energy(0.0)
{
    m_parameters = MergeJson(getDefaultConfig(), config);
    updateGFN2Parameters();

    // Initialize GFN2 core implementation
    m_gfn2 = std::make_unique<GFN2>();
}

json GFN2Method::getDefaultConfig()
{
    return json{
        { "scf_max_iterations", 100 },       // Maximum SCF iterations
        { "scf_threshold", 1.0e-6 },         // Convergence threshold
        { "scf_damping", 0.4 },              // Density damping factor
        { "threads", 1 },                    // Single-threaded by default
        { "print_orbitals", false },         // Print orbital analysis
        { "save_orbitals", false },          // Save orbital data to file
        { "dispersion", "d4" }               // Dispersion correction (stub)
    };
}

// =================================================================================
// Core ComputationalMethod Interface Implementation
// =================================================================================

bool GFN2Method::setMolecule(const Mol& mol)
{
    try {
        // Validate molecule for GFN2 compatibility
        if (!GFN2MethodUtils::isMoleculeSupported(mol)) {
            m_has_error = true;
            m_error_message = "Molecule contains unsupported elements for GFN2";
            return false;
        }

        m_molecule = mol;
        m_calculation_done = false;
        m_initialized = true;
        clearError();

        // Initialize GFN2 with molecule
        if (!m_gfn2->QMInterface::InitialiseMolecule(mol)) {
            handleGFN2Error("molecule initialization");
            return false;
        }

        return true;

    } catch (const std::exception& e) {
        m_has_error = true;
        m_error_message = fmt::format("GFN2 setMolecule failed: {}", e.what());
        return false;
    }
}

bool GFN2Method::updateGeometry(const Matrix& geometry)
{
    if (!m_initialized) {
        m_has_error = true;
        m_error_message = "Method not initialized - call setMolecule first";
        return false;
    }

    try {
        m_molecule.m_geometry = geometry;
        m_calculation_done = false;
        clearError();

        // Update geometry in GFN2 method
        if (!m_gfn2->UpdateMolecule(geometry)) {
            handleGFN2Error("geometry update");
            return false;
        }

        return true;

    } catch (const std::exception& e) {
        m_has_error = true;
        m_error_message = fmt::format("GFN2 updateGeometry failed: {}", e.what());
        return false;
    }
}

double GFN2Method::calculateEnergy(bool gradient)
{
    if (!m_initialized) {
        m_has_error = true;
        m_error_message = "Method not initialized - call setMolecule first";
        return 0.0;
    }

    try {
        clearError();

        // Warn about gradient request if not implemented
        if (gradient) {
            if (CurcumaLogger::get_verbosity() >= 1) {
                CurcumaLogger::warn("GFN2 analytical gradients not yet implemented");
            }
        }

        // Perform GFN2 calculation
        m_last_energy = m_gfn2->Calculation(gradient);
        m_calculation_done = true;

        // Print orbital analysis if requested
        if (m_parameters.value("print_orbitals", false)) {
            if (CurcumaLogger::get_verbosity() >= 2) {
                json decomp = m_gfn2->getEnergyDecomposition();
                CurcumaLogger::info("GFN2 Energy Decomposition:");
                CurcumaLogger::param("electronic", fmt::format("{:.6f} Eh", decomp["electronic"].get<double>()));
                CurcumaLogger::param("repulsion", fmt::format("{:.6f} Eh", decomp["repulsion"].get<double>()));
                CurcumaLogger::param("coulomb", fmt::format("{:.6f} Eh", decomp["coulomb"].get<double>()));
                CurcumaLogger::param("dispersion", fmt::format("{:.6f} Eh", decomp["dispersion"].get<double>()));
            }
        }

        return m_last_energy;

    } catch (const std::exception& e) {
        handleGFN2Error(fmt::format("energy calculation: {}", e.what()));
        return 0.0;
    }
}

// =================================================================================
// Property Access Methods
// =================================================================================

Matrix GFN2Method::getGradient() const
{
    if (!m_calculation_done || !m_gfn2) {
        return Matrix::Zero(m_molecule.m_number_atoms, 3);
    }

    return m_gfn2->Gradient();  // Returns stub for now
}

Vector GFN2Method::getCharges() const
{
    if (!m_calculation_done || !m_gfn2) {
        return Vector::Zero(m_molecule.m_number_atoms);
    }

    return m_gfn2->getPartialCharges();
}

Vector GFN2Method::getBondOrders() const
{
    // GFN2 doesn't explicitly calculate bond orders
    // Could be computed from density matrix if needed
    return Vector{};
}

Position GFN2Method::getDipole() const
{
    // TODO: Calculate dipole moment from charges and positions
    return Position{0.0, 0.0, 0.0};
}

// =================================================================================
// Method Information and Configuration
// =================================================================================

void GFN2Method::setThreadCount(int threads)
{
    m_thread_count = std::max(1, threads);
    // GFN2 can be parallelized at the matrix operations level
}

void GFN2Method::setParameters(const json& params)
{
    m_parameters = MergeJson(m_parameters, params);
    updateGFN2Parameters();
}

json GFN2Method::getParameters() const
{
    return m_parameters;
}

bool GFN2Method::hasError() const
{
    return m_has_error;
}

void GFN2Method::clearError()
{
    m_has_error = false;
    m_error_message.clear();
}

std::string GFN2Method::getErrorMessage() const
{
    return m_error_message;
}

// =================================================================================
// GFN2-specific Methods
// =================================================================================

Vector GFN2Method::getOrbitalEnergies() const
{
    if (!m_calculation_done || !m_gfn2) {
        return Vector{};
    }

    try {
        return m_gfn2->Energies();
    } catch (const std::exception& e) {
        return Vector{};
    }
}

Matrix GFN2Method::getMolecularOrbitals() const
{
    if (!m_calculation_done || !m_gfn2) {
        return Matrix{};
    }

    try {
        return m_gfn2->MolecularOrbitals();
    } catch (const std::exception& e) {
        return Matrix{};
    }
}

int GFN2Method::getNumElectrons() const
{
    if (!m_calculation_done || !m_gfn2) {
        return 0;
    }

    try {
        return m_gfn2->NumElectrons();
    } catch (const std::exception& e) {
        return 0;
    }
}

double GFN2Method::getHOMOLUMOGap() const
{
    if (!m_calculation_done || !m_gfn2) {
        return 0.0;
    }

    try {
        return m_gfn2->getHOMOLUMOGap();
    } catch (const std::exception& e) {
        return 0.0;
    }
}

double GFN2Method::getHOMOEnergy() const
{
    if (!m_calculation_done || !m_gfn2) {
        return 0.0;
    }

    try {
        return m_gfn2->getHOMOEnergy();
    } catch (const std::exception& e) {
        return 0.0;
    }
}

double GFN2Method::getLUMOEnergy() const
{
    if (!m_calculation_done || !m_gfn2) {
        return 0.0;
    }

    try {
        return m_gfn2->getLUMOEnergy();
    } catch (const std::exception& e) {
        return 0.0;
    }
}

Vector GFN2Method::getCoordinationNumbers() const
{
    if (!m_calculation_done || !m_gfn2) {
        return Vector{};
    }

    try {
        return m_gfn2->getCoordinationNumbers();
    } catch (const std::exception& e) {
        return Vector{};
    }
}

json GFN2Method::getEnergyDecomposition() const
{
    if (!m_calculation_done || !m_gfn2) {
        return json{};
    }

    try {
        return m_gfn2->getEnergyDecomposition();
    } catch (const std::exception& e) {
        return json{};
    }
}

bool GFN2Method::saveToFile(const std::string& filename) const
{
    if (!m_calculation_done) {
        return false;
    }

    try {
        // Save GFN2-specific data
        json output_data;
        output_data["method"] = "gfn2";
        output_data["energy"] = m_last_energy;
        output_data["energy_decomposition"] = getEnergyDecomposition();
        output_data["orbital_energies"] = getOrbitalEnergies();
        output_data["coordination_numbers"] = getCoordinationNumbers();
        output_data["num_electrons"] = getNumElectrons();
        output_data["homo_energy"] = getHOMOEnergy();
        output_data["lumo_energy"] = getLUMOEnergy();
        output_data["homo_lumo_gap"] = getHOMOLUMOGap();
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

bool GFN2Method::initializeGFN2()
{
    if (!m_gfn2) {
        m_gfn2 = std::make_unique<GFN2>();
    }
    return m_gfn2 != nullptr;
}

void GFN2Method::updateGFN2Parameters()
{
    if (!m_gfn2) return;

    // Update SCF parameters if available
    // TODO: Pass parameters to GFN2 core implementation
}

void GFN2Method::handleGFN2Error(const std::string& operation)
{
    m_has_error = true;
    m_error_message = fmt::format("GFN2 error during {}", operation);
}

// =================================================================================
// Utility Functions
// =================================================================================

namespace GFN2MethodUtils {

    bool validateGFN2Config(const json& config)
    {
        // Check for required parameters and validate ranges
        if (config.contains("scf_max_iterations")) {
            int max_iter = config["scf_max_iterations"];
            if (max_iter < 1 || max_iter > 1000) {
                return false;
            }
        }

        if (config.contains("scf_threshold")) {
            double threshold = config["scf_threshold"];
            if (threshold <= 0.0 || threshold > 1.0) {
                return false;
            }
        }

        if (config.contains("scf_damping")) {
            double damping = config["scf_damping"];
            if (damping < 0.0 || damping > 1.0) {
                return false;
            }
        }

        return true;
    }

    std::vector<int> getSupportedElements()
    {
        // GFN2 supports Z=1-86 (H to Rn)
        std::vector<int> elements;
        for (int Z = 1; Z <= 86; ++Z) {
            elements.push_back(Z);
        }
        return elements;
    }

    bool isMoleculeSupported(const Mol& mol)
    {
        auto supported = getSupportedElements();

        for (int atom : mol.m_atoms) {
            if (std::find(supported.begin(), supported.end(), atom) == supported.end()) {
                return false;  // Unsupported element found
            }
        }

        return true;
    }
}
