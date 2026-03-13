/*
 * < GFN-FF Method Wrapper Implementation >
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

#include "gfnff.h"
#include "src/core/parameter_macros.h"

#include <memory>

/**
 * @brief Native GFN-FF method wrapper for force field integration
 *
 * This class provides a wrapper around the native GFN-FF implementation
 * that integrates with the force field system in Curcuma.
 */
class GFNFFMethod {
public:
    GFNFFMethod(const json& config);
    ~GFNFFMethod() = default;

    // Core interface methods
    bool setMolecule(const Mol& mol);
    bool updateGeometry(const Matrix& geometry);
    double calculateEnergy(bool gradient = false);
    Matrix getGradient() const;
    Vector getCharges() const;
    Vector getBondOrders() const;
    Position getDipole() const;

    // Configuration methods
    void setThreadCount(int threads);
    void setParameters(const json& params);
    json getParameters() const;

    // Error handling
    bool hasError() const;
    void clearError();
    std::string getErrorMessage() const;

    // Energy component getters (for analysis/debugging)
    double getBondEnergy() const;
    double getAngleEnergy() const;
    double getDihedralEnergy() const;
    double getInversionEnergy() const;
    double getVdWEnergy() const;
    double getRepulsionEnergy() const;
    double getDispersionEnergy() const;
    double getCoulombEnergy() const;

private:
    std::unique_ptr<GFNFF> m_gfnff;
    Mol m_molecule;
    json m_parameters;
    bool m_initialized = false;
    bool m_has_error = false;
    std::string m_error_message;
    int m_thread_count = 1;
    double m_last_energy = 0.0;
    bool m_calculation_done = false;
};