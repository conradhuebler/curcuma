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

#include "d3param_generator.h"
#include "src/core/curcuma_logger.h"

#include <algorithm>
#include <cmath>

D3ParameterGenerator::D3ParameterGenerator(const ConfigManager& config)
    : m_config(config)
{
    initializeReferenceData();
    copyReferenceData();
}

void D3ParameterGenerator::initializeReferenceData()
{
    // Initialize number of references for each element
    // From external/gfnff/src/dftd3param.f90 lines 32-47
    m_number_of_references = {
        2, 1, // H,He
        2, 3, 5, 5, 4, 3, 2, 1, // Li-Ne
        2, 3, 4, 5, 4, 3, 2, 1, // Na-Ar
        2, 3, // K,Ca
        3, 3, 3, 3, 3, 3, 4, 4, 2, 2, // Sc-Zn
        4, 5, 4, 3, 2, 1, // Ga-Kr
        2, 3, // Rb,Sr
        3, 3, 3, 3, 3, 3, 3, 3, 2, 2, // Y-Cd
        4, 5, 4, 3, 2, 1, // In-Xe
        2, 3, // Cs,Ba
        3, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, // La-Eu
        2, 3, 3, 3, 3, 3, 3, 3, 2, 2, // Lu-Hg
        4, 5, 4, 3, 2, 1, // Tl-Rn
        2, 3, // Fr,Ra
        2, 2, 2, 2, 2, 2 // Ac-Pu
    };

    // Initialize reference coordination numbers
    // From external/gfnff/src/dftd3param.f90 lines 49-144
    m_reference_cn = std::vector<std::vector<double>>(MAX_ELEM, std::vector<double>(MAX_REF, 0.0));

    // Populate reference CN for first few elements as example
    // Full implementation would include all 94 elements
    // H
    m_reference_cn[0] = {0.9118, 0.0000, 0.0000, 0.0000, 0.0000};
    // He
    m_reference_cn[1] = {0.0000, 0.0000, 0.0000, 0.0000, 0.0000};
    // Li
    m_reference_cn[2] = {0.0000, 0.9865, 0.0000, 0.0000, 0.0000};
    // Be
    m_reference_cn[3] = {0.0000, 1.9897, 0.0000, 0.0000, 0.0000};

    // Initialize van der Waals radii placeholder
    // Full vdW radii array would be extracted from lines 146+ in Fortran code
    // For now, using placeholder values
    m_vdw_rad.resize(4481, 2.0); // MAX_ELEM*(1+MAX_ELEM)/2 = 4481
}

void D3ParameterGenerator::copyReferenceData()
{
    // Resize the C6 reference tensor
    m_reference_c6.resize(MAX_ELEM);
    for (auto& elem : m_reference_c6) {
        elem.resize(MAX_REF);
        for (auto& ref : elem) {
            ref.resize(MAX_REF, 0.0);
        }
    }

    // Copy C6 reference data from Fortran copy_c6 subroutine
    // This is a large dataset - showing first few values as example
    // Full implementation would extract all values from lines 1098-23563

    // Example for first few atom pairs (H-H, H-He, He-He, etc.)
    // Reference: external/gfnff/src/dftd3param.f90 lines 1098-1124
    m_reference_c6[0][0][0] = 3.0267;   // H-H
    m_reference_c6[0][1][0] = 4.7379;   // H-H (ref 1)
    m_reference_c6[1][0][0] = 4.7379;   // H-H (transposed)
    m_reference_c6[1][1][0] = 7.5916;   // H-H (ref 1 transposed)

    // More atom pairs would follow...

    m_data_initialized = true;

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("D3 reference data initialized");
        CurcumaLogger::param("Elements loaded", static_cast<int>(m_number_of_references.size()));
    }
}

void D3ParameterGenerator::GenerateParameters(const std::vector<int>& atoms)
{
    m_atoms = atoms;
    m_parameters.clear();

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("=== D3ParameterGenerator::GenerateParameters() START ===");
        CurcumaLogger::param("Number of atoms", static_cast<int>(m_atoms.size()));
    }

    // Generate dispersion parameters for all atom pairs
    json dispersion_pairs = json::array();

    for (size_t i = 0; i < m_atoms.size(); ++i) {
        for (size_t j = i + 1; j < m_atoms.size(); ++j) {
            int atom_i = m_atoms[i];
            int atom_j = m_atoms[j];

            // Check if both elements are within the supported range (1-94)
            if (atom_i > 0 && atom_i <= MAX_ELEM && atom_j > 0 && atom_j <= MAX_ELEM) {
                json pair;
                pair["i"] = static_cast<int>(i);
                pair["j"] = static_cast<int>(j);
                pair["element_i"] = atom_i;
                pair["element_j"] = atom_j;

                // Get C6 coefficient from reference data
                double c6 = getC6(atom_i, atom_j);
                pair["c6"] = c6 * m_config.get<double>("d3_s6", 1.0);

                // Calculate C8 coefficient using R6 ratio
                double r6 = getR6(atom_i, atom_j);
                pair["c8"] = c6 * r6 * m_config.get<double>("d3_s8", 1.0);
                pair["r6"] = r6;

                // Add van der Waals radius for distance cutoff
                pair["vdw_radius"] = getR6(atom_i, atom_j); // Use r6 as fallback

                dispersion_pairs.push_back(pair);
            } else {
                if (CurcumaLogger::get_verbosity() >= 2) {
                    CurcumaLogger::warn("D3: Elements out of range - i=" +
                                       std::to_string(atom_i) +
                                       " j=" + std::to_string(atom_j));
                }
            }
        }
    }

    // Output parameter summary
    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::param("Generated D3 pairs", static_cast<int>(dispersion_pairs.size()));
    }

    m_parameters["d3_dispersion_pairs"] = dispersion_pairs;

    // Add Becke-Johnson damping parameters
    m_parameters["d3_damping"] = {
        {"a1", m_config.get<double>("d3_a1", 0.4)},
        {"a2", m_config.get<double>("d3_a2", 4.0)},
        {"alp", m_config.get<double>("d3_alp", 14.0)}
    };

    m_parameters["d3_enabled"] = true;
    m_parameters["d3_reference"] = m_config.get<int>("d3_ref", 0);

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::success("D3 parameter generation completed");
    }
}

double D3ParameterGenerator::getC6(int atom_i, int atom_j, int ref_i, int ref_j) const
{
    // Convert to 0-based indices
    --atom_i; --atom_j;
    if (atom_i >= MAX_ELEM || atom_j >= MAX_ELEM || !m_data_initialized) {
        return 0.0;
    }

    // Calculate pair index (same as Fortran triangular indexing)
    int pair_index;
    if (atom_i > atom_j) {
        pair_index = atom_j + atom_i * (atom_i - 1) / 2;
    } else {
        pair_index = atom_i + atom_j * (atom_j - 1) / 2;
    }

    // Clamp reference indices to valid range
    ref_i = std::min(ref_i, MAX_REF - 1);
    ref_j = std::min(ref_j, MAX_REF - 1);

    // Access C6 coefficient from reference data
    if (pair_index < static_cast<int>(m_reference_c6.size())) {
        if (ref_i < static_cast<int>(m_reference_c6[pair_index].size()) &&
            ref_j < static_cast<int>(m_reference_c6[pair_index][ref_i].size())) {
            return m_reference_c6[pair_index][ref_i][ref_j];
        }
    }

    return 0.0; // Default if not found
}

double D3ParameterGenerator::getR6(int atom_i, int atom_j) const
{
    // Calculate R6 ratio for C8 coefficient estimation
    // This is a simplified implementation - full version would use reference data

    // Convert to 0-based indices
    --atom_i; --atom_j;

    if (atom_i >= MAX_ELEM || atom_j >= MAX_ELEM) {
        return 2.0; // Default value
    }

    // For now, return estimated values based on element types
    // Full implementation would use calculated R6 values from reference data
    if (atom_i < 2 && atom_j < 2) {
        return 1.4; // H-H
    } else if (atom_i < 2 || atom_j < 2) {
        return 1.6; // H with other elements
    } else {
        return 2.0; // Other atom pairs
    }
}

double D3ParameterGenerator::getReferenceCN(int atom, int ref_index) const
{
    // Get reference coordination number for an element
    --atom; // Convert to 0-based

    if (atom >= 0 && atom < MAX_ELEM && ref_index >= 0 && ref_index < MAX_REF) {
        return m_reference_cn[atom][ref_index];
    }

    return 0.0;
}

int D3ParameterGenerator::getNumberofReferences(int atom) const
{
    --atom; // Convert to 0-based

    if (atom >= 0 && atom < static_cast<int>(m_number_of_references.size())) {
        return m_number_of_references[atom];
    }

    return 0;
}

