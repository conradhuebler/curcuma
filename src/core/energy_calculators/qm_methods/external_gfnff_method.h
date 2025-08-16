/*
 * < External GFN-FF Method Wrapper >
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
 * Claude Generated: ComputationalMethod wrapper for external GFN-FF interface
 */

#pragma once

#include "../computational_method.h"
#include "gfnffinterface.h"

#include <memory>

class ExternalGFNFFMethod : public ComputationalMethod {
public:
    ExternalGFNFFMethod(const json& config = json{});
    virtual ~ExternalGFNFFMethod() = default;

    // ComputationalMethod interface
    bool setMolecule(const Mol& mol) override;
    bool updateGeometry(const Matrix& geometry) override;
    double calculateEnergy(bool gradient = false, bool verbose = false) override;
    Matrix getGradient() const override;
    Vector getCharges() const override;
    Vector getBondOrders() const override;
    Position getDipole() const override;
    bool hasGradient() const override { return true; }
    std::string getMethodName() const override { return "External GFN-FF"; }
    bool isThreadSafe() const override { return false; }
    void setThreadCount(int threads) override { (void)threads; }
    void setParameters(const json& params) override;
    json getParameters() const override;
    bool hasError() const override { return false; }

private:
    std::unique_ptr<GFNFFInterface> m_interface;
    json m_config;
    bool m_initialized;
    Mol m_current_molecule;
};