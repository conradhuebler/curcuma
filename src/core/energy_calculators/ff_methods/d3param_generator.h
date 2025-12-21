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

    // Factory methods for standard D3 parameter sets (Claude Generated Phase 3.1)
    /**
     * Create D3ParameterGenerator with GFN-FF parameters
     * Parameters: s6=1.0, s8=2.85, a1=0.80, a2=4.60 (Bohr)
     * Reference: Spicher & Grimme, J. Chem. Theory Comput. 2020
     */
    static D3ParameterGenerator createForGFNFF();

    /**
     * Create D3ParameterGenerator with UFF-D3 (PBE0/BJ) parameters
     * Parameters: s6=1.0, s8=1.2177, a1=0.4145, a2=4.8593 (Bohr)
     * Ideal for hybrid UFF bonded + D3 dispersion
     */
    static D3ParameterGenerator createForUFFD3();

    /**
     * Create D3ParameterGenerator with PBE0-D3-BJ parameters
     * Parameters: s6=1.0, s8=1.2177, a1=0.4145, a2=4.8593 (Bohr)
     * Standard PBE0 functional with D3-BJ damping
     */
    static D3ParameterGenerator createForPBE0();

    /**
     * Create D3ParameterGenerator with BLYP-D3-BJ parameters
     * Parameters: s6=1.0, s8=2.6996, a1=0.4298, a2=4.2359 (Bohr)
     * Reference: Grimme et al., J. Chem. Phys. 132, 154104 (2010), Table II
     */
    static D3ParameterGenerator createForBLYP();

    /**
     * Create D3ParameterGenerator with B3LYP-D3-BJ parameters
     * Parameters: s6=1.0, s8=1.9889, a1=0.3981, a2=4.4211 (Bohr)
     * Reference: Grimme et al., J. Chem. Phys. 132, 154104 (2010), Table II
     */
    static D3ParameterGenerator createForB3LYP();

    /**
     * Create D3ParameterGenerator with TPSS-D3-BJ parameters
     * Parameters: s6=1.0, s8=1.9435, a1=0.4535, a2=4.4752 (Bohr)
     * Reference: Grimme et al., J. Chem. Phys. 132, 154104 (2010), Table II
     */
    static D3ParameterGenerator createForTPSS();

    /**
     * Create D3ParameterGenerator with PBE-D3-BJ parameters
     * Parameters: s6=1.0, s8=0.7875, a1=0.4289, a2=4.4407 (Bohr)
     * Reference: Grimme et al., J. Chem. Phys. 132, 154104 (2010), Table II
     */
    static D3ParameterGenerator createForPBE();

    /**
     * Create D3ParameterGenerator with BP86-D3-BJ parameters
     * Parameters: s6=1.0, s8=3.2822, a1=0.3946, a2=4.8516 (Bohr)
     * Reference: Grimme et al., J. Chem. Phys. 132, 154104 (2010), Table II
     */
    static D3ParameterGenerator createForBP86();

    /**
     * Create D3ParameterGenerator with method-specific parameters
     * Supported methods: "pbe0", "blyp", "b3lyp", "tpss", "pbe", "bp86", "gfnff", "uff-d3"
     * Throws std::invalid_argument if method is unknown
     */
    static D3ParameterGenerator createForMethod(const std::string& method);

    // Main parameter generation interface
    void GenerateParameters(const std::vector<int>& atoms, const Eigen::MatrixXd& geometry);
    json getParameters() const { return m_parameters; }

    // Energy calculation
    double getTotalEnergy() const;  // Returns total D3 dispersion energy (Eh)

    // Individual parameter accessors
    // NOTE: Default ref indices use first valid reference (not 0,0 which is empty)
    double getC6(int atom_i, int atom_j, int ref_i = 0, int ref_j = 1) const;
    double getR6(int atom_i, int atom_j) const;
    double getReferenceCN(int atom, int ref_index) const;
    int getNumberofReferences(int atom) const;

private:
    void initializeReferenceData();
    void copyReferenceData();
    double calculateR6(int atom_i, int atom_j) const;

    // Geometry-dependent methods (NOT YET IMPLEMENTED)
    // TODO (Consolidation): Consider sharing CN calculation with GFN-FF/EEQ
    // - GFN-FF likely has CN calculation for dispersion terms
    // - EEQSolver may have geometry-dependent CN calculation
    // - Investigate creating shared CNCalculator utility class
    std::vector<double> calculateCoordinationNumbers(const std::vector<int>& atoms, const Eigen::MatrixXd& geometry) const;
    double interpolateC6(int elem_i, int elem_j, double cn_i, double cn_j) const;
    double calculateDistance(int i, int j, const Eigen::MatrixXd& geometry) const;

    // Reference data lookup methods
    double getC6Coefficient(int elem1, int elem2, double cn1, double cn2) const;
    int getElementPairIndex(int elem1, int elem2) const;

    // Reference data from s-dftd3 (COMPLETE MAX_REF=7)
    static const int MAX_ELEM = 103;  // Updated to match s-dftd3 (103 vs 94)
    static const int MAX_REF = 7;     // CRITICAL FIX: s-dftd3 uses 7 references

    std::vector<std::vector<std::vector<double>>> m_reference_c6;  // Legacy structure
    std::vector<double> m_reference_c6_flat;  // Flat storage for 2375 extracted values
    std::vector<int> m_number_of_references;
    std::vector<std::vector<double>> m_reference_cn;
    std::vector<double> m_vdw_rad;
    std::vector<double> m_sqrt_z_r4_over_r2;  // Transformed R8/R6 data from simple-dftd3

    ConfigManager m_config;
    json m_parameters;
    std::vector<int> m_atoms;
    Eigen::MatrixXd m_geometry;  // Store geometry for distance calculations
    bool m_data_initialized = false;
};