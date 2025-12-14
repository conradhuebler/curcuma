/*
 * DFT-D4 Parameter Generator for Curcuma
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
 * Claude Generated 2025 - Native D4 parameter generation from GFN-FF Fortran reference
 */

#pragma once

#include "src/core/global.h"
#include "src/core/parameter_macros.h"
#include "src/core/config_manager.h"

#include <Eigen/Dense>

#include "json.hpp"
using json = nlohmann::json;

BEGIN_PARAMETER_DEFINITION(d4param)
    // D4 reference selection and scaling
    PARAM(d4_refq, Int, 2, "D4 reference charges (0=gfn2xtb, 1=gasteiger, 2=hirshfeld).", "Reference", {})
    PARAM(d4_s6, Double, 1.0, "D4 global scaling factor for C6 term.", "Scaling", {})
    PARAM(d4_s8, Double, 1.0, "D4 global scaling factor for C8 term.", "Scaling", {})
    PARAM(d4_s10, Double, 1.0, "D4 global scaling factor for C10 term.", "Scaling", {})
    PARAM(d4_s12, Double, 1.0, "D4 global scaling factor for C12 term.", "Scaling", {})
    PARAM(d4_a1, Double, 0.43, "D4 damping parameter a1.", "Damping", {})
    PARAM(d4_a2, Double, 4.0, "D4 damping parameter a2 (Bohr).", "Damping", {})
    PARAM(d4_alp, Double, 14.0, "D4 alpha damping parameter.", "Damping", {})

    // D4 atomic polarizability data
    PARAM(d4_r4r2_model, Int, 1, "D4 r4/r2 model (0=static, 1=dynamic).", "Model", {})
END_PARAMETER_DEFINITION

class D4ParameterGenerator {
public:
    D4ParameterGenerator(const ConfigManager& config);
    ~D4ParameterGenerator() = default;

    // Main parameter generation interface
    void GenerateParameters(const std::vector<int>& atoms);
    json getParameters() const { return m_parameters; }

    // Individual parameter accessors
    double getC6(int atom_i, int atom_j, int ref_i = 0, int ref_j = 0) const;
    double getR4OverR2(int atom) const;
    double getSqrtZr4r2(int atom) const;
    double getAtomicPolarizability(int atom, int frequency_index = 1) const;

private:
    void initializeReferenceData();
    void calculateFrequencyDependentPolarizabilities();
    double getEffectiveC6(int atom_i, int atom_j) const;

    // Reference data from GFN-FF Fortran implementation
    static const int MAX_ELEM = 118;
    static const int MAX_REF = 7;
    static const int N_FREQ = 23;  // Frequency grid points
    static const int N_REFQ = 17;  // Reference charge states

    // Atomic properties from GFN-FF
    std::vector<double> m_r4_over_r2;
    std::vector<double> m_sqrt_z_r4_r2;

    // Frequency-dependent polarizability data
    std::vector<std::vector<std::vector<double>>> m_alpha_iw;
    std::vector<std::vector<double>> m_integrated_alpha;
    std::vector<int> m_refn;    // Number of reference systems
    std::vector<std::vector<double>> m_refq, m_refh, m_refcn;  // Reference charges/hardness/CN
    std::vector<std::vector<double>> m_ascale;  // Atomic scaling factors

    ConfigManager m_config;
    json m_parameters;
    std::vector<int> m_atoms;
    bool m_data_initialized = false;

    // Mathematical constants
    static constexpr double PI = 3.14159265358979323846;
    static constexpr double THREE_OVER_PI = 3.0 / PI;
    static constexpr double HALF_OVER_PI = 0.5 / PI;
};