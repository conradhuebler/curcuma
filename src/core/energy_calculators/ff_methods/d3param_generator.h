/*
 * DFT-D3 Parameter Generator for Curcuma
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
 * Claude Generated 2025 - Native D3 parameter generation from GFN-FF Fortran reference
 */

#pragma once

#include "src/core/global.h"
#include "src/core/parameter_macros.h"
#include "src/core/config_manager.h"

#include <Eigen/Dense>

#include "json.hpp"
using json = nlohmann::json;

BEGIN_PARAMETER_DEFINITION(d3param)
    // D3 Reference selection and scaling
    PARAM(d3_ref, Int, 0, "D3 reference coordination number set (0=standard).", "Reference", {})
    PARAM(d3_s6, Double, 1.0, "D3 global scaling factor for C6 term.", "Scaling", {})
    PARAM(d3_s8, Double, 1.0, "D3 global scaling factor for C8 term.", "Scaling", {})
    PARAM(d3_rs6, Double, 1.0, "D3 scaling for C6 damping function.", "Scaling", {})
    PARAM(d3_rs8, Double, 1.0, "D3 scaling for C8 damping function.", "Scaling", {})
    PARAM(d3_s9, Double, 1.0, "D3 scaling for three-body ATM term.", "Scaling", {})

    // D3 damping parameters (Becke-Johnson)
    PARAM(d3_a1, Double, 0.4, "D3 BJ damping parameter a1.", "Damping", {})
    PARAM(d3_a2, Double, 4.0, "D3 BJ damping parameter a2 (Bohr).", "Damping", {})
    PARAM(d3_alp, Double, 14.0, "D3 BJ alpha parameter.", "Damping", {})
END_PARAMETER_DEFINITION

class D3ParameterGenerator {
public:
    D3ParameterGenerator(const ConfigManager& config);
    ~D3ParameterGenerator() = default;

    // Main parameter generation interface
    void GenerateParameters(const std::vector<int>& atoms);
    json getParameters() const { return m_parameters; }

    // Individual parameter accessors
    double getC6(int atom_i, int atom_j, int ref_i = 0, int ref_j = 0) const;
    double getR6(int atom_i, int atom_j) const;
    double getReferenceCN(int atom, int ref_index) const;
    int getNumberofReferences(int atom) const;

private:
    void initializeReferenceData();
    void copyReferenceData();
    double calculateR6(int atom_i, int atom_j) const;

    // Reference data from GFN-FF Fortran implementation
    static const int MAX_ELEM = 94;
    static const int MAX_REF = 5;

    std::vector<std::vector<std::vector<double>>> m_reference_c6;
    std::vector<int> m_number_of_references;
    std::vector<std::vector<double>> m_reference_cn;
    std::vector<double> m_vdw_rad;

    ConfigManager m_config;
    json m_parameters;
    std::vector<int> m_atoms;
    bool m_data_initialized = false;
};