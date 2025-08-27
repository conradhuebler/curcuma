/*
 * < Dispersion Correction Method Wrapper for ComputationalMethod Interface >
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

#pragma once

#include "../computational_method.h"

#ifdef USE_D3
#include "dftd3interface.h"
#endif

#ifdef USE_D4
#include "dftd4interface.h"
#endif

#include <memory>

/**
 * @brief Dispersion correction method wrapper for ComputationalMethod interface
 * 
 * This wrapper adapts the existing DFT-D3 and DFT-D4 interfaces to the unified
 * ComputationalMethod interface. These methods provide dispersion corrections
 * for quantum chemical calculations.
 * 
 * Supported methods:
 * - d3: DFT-D3 dispersion correction
 * - d4: DFT-D4 dispersion correction
 * 
 * Claude Generated: Big-Bang EnergyCalculator refactoring wrapper
 */
class DispersionMethod : public ComputationalMethod {
public:
    DispersionMethod(const std::string& method_name, const json& config = json{});
    virtual ~DispersionMethod() = default;
    
    // ComputationalMethod interface
    bool setMolecule(const Mol& mol) override;
    bool updateGeometry(const Matrix& geometry) override;
    double calculateEnergy(bool gradient = false) override;

    Matrix getGradient() const override;
    Vector getCharges() const override;
    Vector getBondOrders() const override;
    Position getDipole() const override;
    bool hasGradient() const override { return true; }
    
    std::string getMethodName() const override { return m_method_name; }
    bool isThreadSafe() const override { return true; }
    void setThreadCount(int threads) override;
    
    void setParameters(const json& params) override;
    json getParameters() const override;
    bool hasError() const override;
    void clearError() override;
    std::string getErrorMessage() const override;

    static bool isD3Available();
    static bool isD4Available();
    static std::vector<std::string> getSupportedMethods();

private:
    std::unique_ptr<QMInterface> m_dispersion;  // Either DFTD3Interface or DFTD4Interface
    std::string m_method_name;
    Mol m_molecule;
    bool m_calculation_done;
    double m_last_energy;
};