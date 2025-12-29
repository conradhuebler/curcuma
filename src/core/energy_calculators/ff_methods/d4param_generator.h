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
#include "eeq_solver.h"  // EEQ charge calculation for D4
#include "cn_calculator.h"  // CN calculation for D4 (Claude Generated - December 2025)

#include <Eigen/Dense>
#include <memory>
#include <unordered_map>

#include "json.hpp"
using json = nlohmann::json;

BEGIN_PARAMETER_DEFINITION(d4param)
    // D4 reference selection and scaling
    PARAM(d4_refq, Int, 2, "D4 reference charges (0=gfn2xtb, 1=gasteiger, 2=hirshfeld).", "Reference", {})
    PARAM(d4_s6, Double, 1.0, "D4 global scaling factor for C6 term.", "Scaling", {})
    PARAM(d4_s8, Double, 1.0, "D4 global scaling factor for C8 term.", "Scaling", {})
    PARAM(d4_s10, Double, 1.0, "D4 global scaling factor for C10 term.", "Scaling", {})
    PARAM(d4_s12, Double, 1.0, "D4 global scaling factor for C12 term.", "Scaling", {})
    PARAM(d4_s9, Double, 1.0, "D4 scaling for three-body ATM term (default: enabled).", "Scaling", {})
    PARAM(d4_a1, Double, 0.44, "D4 damping parameter a1 (GFN-FF: 0.44, GFN2-xTB: 0.63).", "Damping", {})
    PARAM(d4_a2, Double, 4.60, "D4 damping parameter a2 (Bohr) - GFN-FF: 4.60, GFN2-xTB: 5.0.", "Damping", {})
    PARAM(d4_alp, Double, 14.0, "D4 alpha damping parameter.", "Damping", {})

    // D4 atomic polarizability data
    PARAM(d4_r4r2_model, Int, 1, "D4 r4/r2 model (0=static, 1=dynamic).", "Model", {})
END_PARAMETER_DEFINITION

class D4ParameterGenerator {
public:
    D4ParameterGenerator(const ConfigManager& config);
    ~D4ParameterGenerator() = default;

    // Main parameter generation interface
    void GenerateParameters(const std::vector<int>& atoms, const Matrix& geometry_bohr);
    json getParameters() const { return m_parameters; }

    // Individual parameter accessors
    double getC6(int atom_i, int atom_j, int ref_i = 0, int ref_j = 0) const;
    double getR4OverR2(int atom) const;
    double getSqrtZr4r2(int atom) const;
    double getAtomicPolarizability(int atom, int frequency_index = 1) const;

private:
    void initializeReferenceData();
    void calculateFrequencyDependentPolarizabilities();

    // OLD: Crude neutral-atom approximation (deprecated)
    double getEffectiveC6(int atom_i, int atom_j) const;

    // Claude Generated (Dec 27, 2025): Weight caching optimization
    void precomputeGaussianWeights();

    // Claude Generated (Dec 27, 2025): C6 reference matrix pre-computation
    void precomputeC6ReferenceMatrix();
    double computeC6Reference(int elem_i, int elem_j, int ref_i, int ref_j) const;

    // NEW: Charge-weighted C6 using EEQ charges and Gaussian weighting (Dec 2025)
    // Phase 2.2 (December 2025): CN+charge combined weighting
    // Claude Generated (Dec 27, 2025): Now uses cached weights and C6 reference for performance
    double getChargeWeightedC6(int Zi, int Zj, size_t atom_i, size_t atom_j) const;

    // ATM three-body helper (Claude Generated 2025)
    double calculateTripleScale(int i, int j, int k) const;

    // Reference data from GFN-FF Fortran implementation
    static const int MAX_ELEM = 118;
    static const int MAX_REF = 7;
    static const int N_FREQ = 23;  // Frequency grid points
    static const int N_REFQ = 17;  // Reference charge states

    // D4 reference data from external/gfnff/src/dftd4param.f90
    std::vector<int> m_refn;    // Number of reference systems per element
    std::vector<std::vector<double>> m_refq, m_refh;  // Reference charges and hydrogen counts
    std::vector<std::vector<double>> m_refcn;   // Reference coordination numbers (cpp-d4)

    // Legacy placeholders (deprecated - replaced with real data)
    std::vector<double> m_r4_over_r2;
    std::vector<double> m_sqrt_z_r4_r2;
    std::vector<std::vector<std::vector<double>>> m_alpha_iw;
    std::vector<std::vector<double>> m_integrated_alpha;

    ConfigManager m_config;
    json m_parameters;
    std::vector<int> m_atoms;
    bool m_data_initialized = false;

    // EEQ charge calculation (Dec 2025 - Phase 2 D4 integration)
    std::unique_ptr<EEQSolver> m_eeq_solver;
    Vector m_eeq_charges;  // Cached EEQ charges for current geometry

    // Molecular CN calculation (December 2025 Phase 2 - CN+Charge weighting)
    std::vector<double> m_cn_values;  // Cached GFNFFCN values for current geometry

    // Claude Generated (Dec 27, 2025): Cached Gaussian weights for C6 interpolation
    // Performance optimization: Pre-compute CN+charge weights once per atom instead of per pair
    // Expected speedup: 5-10x for molecules with 100+ atoms (same as D3)
    std::vector<std::vector<double>> m_gaussian_weights;  // [atom_idx][ref_idx]
    bool m_weights_cached = false;

    // Claude Generated (Dec 27, 2025): Cached C6 reference matrix
    // Performance optimization: Pre-compute Casimir-Polder integration ONCE per element-pair combination
    // Key: (elem_i << 24) | (elem_j << 16) | (ref_i << 8) | ref_j
    // Expected speedup: 50-100x for large molecules (eliminates ~340,000 integration operations)
    std::unordered_map<uint32_t, double> m_c6_reference_cache;
    bool m_c6_reference_cached = false;

    // Helper to generate cache key
    inline uint32_t c6CacheKey(int elem_i, int elem_j, int ref_i, int ref_j) const {
        return (static_cast<uint32_t>(elem_i) << 24) |
               (static_cast<uint32_t>(elem_j) << 16) |
               (static_cast<uint32_t>(ref_i) << 8) |
               static_cast<uint32_t>(ref_j);
    }

    // Mathematical constants
    static constexpr double PI = 3.14159265358979323846;
    static constexpr double THREE_OVER_PI = 3.0 / PI;
    static constexpr double HALF_OVER_PI = 0.5 / PI;
};