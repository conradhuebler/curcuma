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

#pragma once

#include "../computational_method.h"
#include "../ff_methods/gfnff.h"

#include <memory>

/**
 * @brief ComputationalMethod adapter for GFN-FF
 *
 * This wrapper provides the ComputationalMethod interface for the native
 * GFN-FF implementation, allowing it to be used through the EnergyCalculator
 * polymorphic system.
 */

  // TODO: Move GFNFF class to the force field methods directory
  // TODO: Is this obsolete now that GFNFF is located in ff_methods/gfnff.h or somewhere
class GFNFFComputationalMethod : public ComputationalMethod {
public:
    GFNFFComputationalMethod(const std::string& method_name, const json& config);
    ~GFNFFComputationalMethod() = default;

    // ComputationalMethod interface
    bool setMolecule(const Mol& mol) override;
    double calculateEnergy(bool gradient = false) override;
    Matrix getGradient() const override;
    Vector getCharges() const override;
    Vector getBondOrders() const override;
    Position getDipole() const override;
    bool updateGeometry(const Matrix& geometry) override;
    bool hasGradient() const override { return true; }
    bool isThreadSafe() const override { return false; }
    std::string getMethodName() const override { return "cgfnff"; }

    // Configuration
    void setThreadCount(int threads) override;
    void setParameters(const json& params) override;
    json getParameters() const override;

    // Error handling
    bool hasError() const override;
    void clearError() override;
    std::string getErrorMessage() const override;

private:
    std::unique_ptr<GFNFF> m_gfnff;
    json m_parameters;
    bool m_has_error = false;
    std::string m_error_message;
    double m_last_energy = 0.0;
};