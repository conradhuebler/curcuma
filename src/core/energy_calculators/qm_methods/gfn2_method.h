/*
 * < GFN2 Method Wrapper for ComputationalMethod Interface >
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Based on the GFN2-xTB method developed by:
 *   Stefan Grimme, Christoph Bannwarth, Sebastian Ehlert
 *   Mulliken Center for Theoretical Chemistry, University of Bonn
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 */

#pragma once

#include "../computational_method.h"
#include "gfn2.h"

#include <memory>

/**
 * @brief GFN2-xTB wrapper for ComputationalMethod interface
 *
 * This wrapper adapts the native GFN2 implementation to the unified
 * ComputationalMethod interface, providing consistent access to GFN2
 * calculations within the EnergyCalculator architecture.
 *
 * Features:
 * - Native C++ implementation of GFN2-xTB tight-binding method
 * - Excellent accuracy for geometries and noncovalent interactions
 * - ~1000x faster than DFT for large systems
 * - Educational transparency: full algorithm visibility
 * - Thread-safe by design
 *
 * Method Characteristics:
 * - Minimal valence basis (STOs)
 * - Self-consistent charge treatment
 * - Third-order electrostatics
 * - D4 dispersion correction (when integrated)
 * - Analytical gradients (when implemented)
 *
 * Reference:
 *   C. Bannwarth, S. Ehlert, S. Grimme
 *   J. Chem. Theory Comput. 2019, 15, 1652-1671
 *
 * Claude Generated: Native GFN2 method wrapper for Curcuma
 */
class GFN2Method : public ComputationalMethod {
public:
    /**
     * @brief Constructor with configuration
     * @param config JSON configuration (GFN2-specific parameters)
     */
    explicit GFN2Method(const json& config = json{});

    /**
     * @brief Destructor
     */
    virtual ~GFN2Method() = default;

    // =================================================================================
    // Core ComputationalMethod Interface Implementation
    // =================================================================================

    bool setMolecule(const Mol& mol) override;
    bool updateGeometry(const Matrix& geometry) override;
    double calculateEnergy(bool gradient = false) override;

    // Property access
    Matrix getGradient() const override;
    Vector getCharges() const override;
    Vector getBondOrders() const override;
    Position getDipole() const override;
    bool hasGradient() const override { return false; }  // TODO: implement gradients

    // Method information
    std::string getMethodName() const override { return "gfn2"; }
    bool isThreadSafe() const override { return true; }
    void setThreadCount(int threads) override;

    // Configuration
    void setParameters(const json& params) override;
    json getParameters() const override;
    bool hasError() const override;
    void clearError() override;
    std::string getErrorMessage() const override;

    // =================================================================================
    // GFN2-specific Methods (additional functionality)
    // =================================================================================

    /**
     * @brief Get orbital energies
     * @return Vector of orbital energies in eV
     */
    Vector getOrbitalEnergies() const override;

    /**
     * @brief Get molecular orbital coefficients
     * @return MO coefficient matrix
     */
    Matrix getMolecularOrbitals() const;

    /**
     * @brief Get number of electrons
     * @return Total number of electrons
     */
    int getNumElectrons() const override;

    /**
     * @brief Get HOMO-LUMO gap
     * @return Gap in eV
     */
    double getHOMOLUMOGap() const;

    /**
     * @brief Get HOMO energy
     * @return HOMO energy in eV
     */
    double getHOMOEnergy() const;

    /**
     * @brief Get LUMO energy
     * @return LUMO energy in eV
     */
    double getLUMOEnergy() const;

    /**
     * @brief Get coordination numbers
     * @return Vector of CN values
     */
    Vector getCoordinationNumbers() const;

    /**
     * @brief Get energy decomposition
     * @return JSON with energy components
     */
    json getEnergyDecomposition() const;

    /**
     * @brief Save calculation results to file
     * @param filename Output file name
     * @return true if save successful
     */
    bool saveToFile(const std::string& filename) const override;

private:
    std::unique_ptr<GFN2> m_gfn2;              ///< Wrapped GFN2 implementation
    Mol m_molecule;                             ///< Current molecule
    bool m_calculation_done;                    ///< Flag if calculation was performed
    double m_last_energy;                       ///< Last calculated energy

    /**
     * @brief Initialize GFN2 method with current parameters
     * @return true if initialization successful
     */
    bool initializeGFN2();

    /**
     * @brief Update GFN2 parameters from JSON configuration
     */
    void updateGFN2Parameters();

    /**
     * @brief Convert GFN2-specific errors to ComputationalMethod format
     */
    void handleGFN2Error(const std::string& operation);

    /**
     * @brief Default GFN2 configuration
     */
    static json getDefaultConfig();
};

/**
 * @brief Utility functions for GFN2 method
 */
namespace GFN2MethodUtils {
    /**
     * @brief Validate GFN2 parameters
     * @param config JSON configuration to validate
     * @return true if parameters are valid
     */
    bool validateGFN2Config(const json& config);

    /**
     * @brief Get supported elements for GFN2
     * @return Vector of atomic numbers supported by GFN2
     */
    std::vector<int> getSupportedElements();

    /**
     * @brief Check if molecule contains only supported elements
     * @param mol Molecule to check
     * @return true if all elements are supported
     */
    bool isMoleculeSupported(const Mol& mol);
}
