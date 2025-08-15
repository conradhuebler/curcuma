/*
 * < Native GFN-FF Method Wrapper for ComputationalMethod Interface >
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
#include "src/core/qm_methods/gfnff.h"

#include <memory>

/**
 * @brief Native GFN-FF method wrapper for ComputationalMethod interface
 * 
 * This wrapper adapts Curcuma's native GFN-FF implementation (cgfnff)
 * to the unified ComputationalMethod interface. This is the "cgfnff" 
 * method - Curcuma's own implementation of the GFN-FF force field.
 * 
 * Status: WORK IN PROGRESS - Architecture complete, parameter generation debugging needed
 * 
 * Claude Generated: Big-Bang EnergyCalculator refactoring wrapper
 */
class GFNFFMethod : public ComputationalMethod {
public:
    explicit GFNFFMethod(const json& config = json{});
    virtual ~GFNFFMethod() = default;
    
    // ComputationalMethod interface
    bool setMolecule(const Mol& mol) override;
    bool updateGeometry(const Matrix& geometry) override;
    double calculateEnergy(bool gradient = false, bool verbose = false) override;
    
    Matrix getGradient() const override;
    Vector getCharges() const override;
    Vector getBondOrders() const override;
    Position getDipole() const override;
    bool hasGradient() const override { return true; }
    
    std::string getMethodName() const override { return "cgfnff"; }
    bool isThreadSafe() const override { return true; }
    void setThreadCount(int threads) override;
    
    void setParameters(const json& params) override;
    json getParameters() const override;
    bool hasError() const override;
    void clearError() override;
    std::string getErrorMessage() const override;

private:
    std::unique_ptr<GFNFF> m_gfnff;
    Mol m_molecule;
    bool m_calculation_done;
    double m_last_energy;
};