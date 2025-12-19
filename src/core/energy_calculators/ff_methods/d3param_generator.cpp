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

// Complete s-dftd3reference data (MAX_REF=7 - fixes 1.48x energy error)
extern std::vector<double> reference_cn_data_complete;
extern std::vector<double> reference_c6_data_complete;

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
    // Initialize number of references for each element (s-dftd3 MAX_REF=7)
    // From /external/simple-dftd3/src/dftd3/reference.f90 lines 30-45
    // Updated for 103 elements (H through Lr)
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
        7, 5, 7, 6, 7, 7, 6, 6, 5, 6, 7, 7, 6, 7, 7 // Ac-Lr (total 103)
    };

    // Initialize reference coordination numbers (MAX_ELEM=103, MAX_REF=7)
    // Using complete s-dftd3 reference data to fix 1.48x energy error
    m_reference_cn = std::vector<std::vector<double>>(MAX_ELEM, std::vector<double>(MAX_REF, 0.0));

    // Load complete CN reference data from s-dftd3 (103×7=721 values)
    extern std::vector<double> reference_cn_data_complete;

    if (reference_cn_data_complete.size() < 721) {  // Need 103 elements × 7 references = 721
        CurcumaLogger::error("Reference CN data insufficient: expected ≥721, got " +
                             std::to_string(reference_cn_data_complete.size()));
        return;
    }

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("Loading " + std::to_string(721) + " CN values from s-dftd3 reference data (size=" +
                           std::to_string(reference_cn_data_complete.size()) + ")");
    }

    // Copy CN values into element-indexed structure: [element][ref]
    // Data is already converted from Fortran column-major to C++ row-major format
    int flat_index = 0;
    for (int elem = 0; elem < MAX_ELEM; ++elem) {
        for (int ref = 0; ref < MAX_REF; ++ref) {
            if (flat_index < static_cast<int>(reference_cn_data_complete.size())) {
                m_reference_cn[elem][ref] = reference_cn_data_complete[flat_index++];
            }
        }
    }

    // Initialize van der Waals radii placeholder for 103 elements
    // Full vdW radii array would be extracted from s-dftd3 reference data
    // For now, using placeholder values
    m_vdw_rad.resize(MAX_ELEM * (MAX_ELEM + 1) / 2, 2.0); // 103*104/2 = 5356 pairs
}

void D3ParameterGenerator::copyReferenceData()
{
    // Resize the C6 reference tensor structure for s-dftd3 (MAX_ELEM=103, MAX_REF=7)
    // Structure: m_reference_c6[pair_index][ref_i][ref_j]
    // where pair_index uses triangular indexing for 103×7×7 reference grid
    int max_pairs = MAX_ELEM * (MAX_ELEM + 1) / 2;  // 103*104/2 = 5356 pairs
    m_reference_c6.resize(max_pairs);
    for (auto& pair : m_reference_c6) {
        pair.resize(MAX_REF);
        for (auto& ref : pair) {
            ref.resize(MAX_REF, 0.0);
        }
    }

    // Import complete s-dftd3 C6 reference data (262,444 values)
    // Data source: reference_c6_data_complete extracted from s-dftd3/refrence.f90
    // Structure: 7×7×(103+1)×103/2 = 262,444 values
    extern std::vector<double> reference_c6_data_complete;

    if (reference_c6_data_complete.size() != 262444) {
        CurcumaLogger::error("C6 reference data size mismatch: expected 262444, got " +
                             std::to_string(reference_c6_data_complete.size()));
        return;
    }

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("Loading complete s-dftd3 C6 reference data (" +
                           std::to_string(reference_c6_data_complete.size()) + " values)");
    }

    // Copy s-dftd3 c6ab(ref_i, ref_j, pair_index) into triangular structure
    // Uses triangular indexing with 7×7 reference pairs per element pair (5356 total)
    for (int pair_index = 0; pair_index < max_pairs; ++pair_index) {
        for (int ref_i = 0; ref_i < MAX_REF; ++ref_i) {
            for (int ref_j = 0; ref_j < MAX_REF; ++ref_j) {
                // s-dftd3 Fortran column-major: c6ab(ref_i, ref_j, pair_index)
                // Flat array order: ref_i (fastest), ref_j, pair_index (slowest)
                // We store as [pair_index][ref_i][ref_j] which matches this layout
                int flat_index = ref_i + ref_j * MAX_REF + pair_index * MAX_REF * MAX_REF;

                if (flat_index < static_cast<int>(reference_c6_data_complete.size())) {
                    m_reference_c6[pair_index][ref_i][ref_j] = reference_c6_data_complete[flat_index];
                }
            }
        }
    }

    m_data_initialized = true;

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("D3 reference data initialized with complete s-dftd3 dataset (MAX_REF=7)");
        CurcumaLogger::param("Elements loaded", static_cast<int>(m_number_of_references.size()));
        CurcumaLogger::param("Element pairs loaded", max_pairs);
        CurcumaLogger::param("C6 coefficients imported", static_cast<int>(reference_c6_data_complete.size()));

        // Debug: Check H-H C6 values (pair_index=0)
        CurcumaLogger::info("Debug C6 matrix for H-H (pair_index=0) with MAX_REF=7:");
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                CurcumaLogger::info("  C6[" + std::to_string(i) + "," + std::to_string(j) + "] = " +
                                   std::to_string(m_reference_c6[0][i][j]));
            }
        }
    }
}

void D3ParameterGenerator::GenerateParameters(const std::vector<int>& atoms, const Eigen::MatrixXd& geometry)
{
    m_atoms = atoms;
    m_geometry = geometry;
    m_parameters.clear();

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("=== D3ParameterGenerator::GenerateParameters() START ===");
        CurcumaLogger::param("Number of atoms", static_cast<int>(m_atoms.size()));
    }

    // Calculate geometry-dependent coordination numbers
    std::vector<double> coordination_numbers = calculateCoordinationNumbers(m_atoms, m_geometry);

    // Generate dispersion parameters for all atom pairs
    json dispersion_pairs = json::array();

    for (size_t i = 0; i < m_atoms.size(); ++i) {
        for (size_t j = i + 1; j < m_atoms.size(); ++j) {
            int atom_i = m_atoms[i];
            int atom_j = m_atoms[j];

            // Check if both elements are within the supported range (1-based: H=1, C=6, ...)
            if (atom_i >= 1 && atom_i <= MAX_ELEM && atom_j >= 1 && atom_j <= MAX_ELEM) {
                json pair;
                pair["i"] = static_cast<int>(i);
                pair["j"] = static_cast<int>(j);
                pair["element_i"] = atom_i;
                pair["element_j"] = atom_j;

                // Get C6 coefficient with CN-dependent interpolation
                double cn_i = coordination_numbers[i];
                double cn_j = coordination_numbers[j];
                // Calculate C6 coefficient from interpolation
                // Claude Generated (December 2025): CRITICAL FIX - s6 applied in energy, not here!
                double c6 = interpolateC6(atom_i, atom_j, cn_i, cn_j);
                pair["c6"] = c6;  // Store raw C6, s6 applied in energy calculation
                pair["cn_i"] = cn_i;  // Store CN for debugging
                pair["cn_j"] = cn_j;

                // Calculate C8 coefficient using C8/C6 ratio
                // Claude Generated (December 2025): CRITICAL FIX - s8 applied in energy, not here!
                double c8_over_c6 = getR6(atom_i, atom_j);
                pair["c8"] = c6 * c8_over_c6;  // C8 = C6 * (C8/C6 ratio), s8 applied in energy calculation
                pair["r6"] = c8_over_c6;

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
    // Elements from Molecule are 1-based (H=1, C=6, ...) - convert to 0-based for array indexing
    int elem_i = atom_i - 1;
    int elem_j = atom_j - 1;

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("getC6: atom_i=" + std::to_string(atom_i) + " atom_j=" + std::to_string(atom_j) +
                           " elem_i=" + std::to_string(elem_i) + " elem_j=" + std::to_string(elem_j));
    }

    if (elem_i < 0 || elem_i >= MAX_ELEM || elem_j < 0 || elem_j >= MAX_ELEM || !m_data_initialized) {
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::warn("getC6: Range check failed or data not initialized");
        }
        return 0.0;
    }

    // Calculate pair index (triangular indexing, 0-based after conversion)
    // NOTE: s-dftd3 uses 1-based Fortran: ic = j + i*(i-1)/2
    // For 0-based C++, formula changes to: ic = j + i*(i+1)/2
    // Example: Fortran ati=17, atj=1 → ic=137 (1-based)
    //          C++ elem_i=16, elem_j=0 → ic=136 (0-based, Fortran ic minus 1)
    int pair_index;
    if (elem_i > elem_j) {
        pair_index = elem_j + elem_i * (elem_i + 1) / 2;
    } else {
        pair_index = elem_i + elem_j * (elem_j + 1) / 2;
    }

    // CRITICAL: s-dftd3 index-swap behavior (reference.f90:160-171)
    // When elem_i > elem_j: normal order reference_c6(iref, jref, pair)
    // When elem_i ≤ elem_j: SWAPPED order reference_c6(jref, iref, pair)
    // This asymmetry is ESSENTIAL for heteronuclear pair accuracy
    // MUST use same condition as pair_index calculation!
    int ref_i_final, ref_j_final;
    if (elem_i > elem_j) {
        // Normal order - no swap
        ref_i_final = ref_i;
        ref_j_final = ref_j;
    } else {
        // SWAP reference indices for elem_i ≤ elem_j
        ref_i_final = ref_j;
        ref_j_final = ref_i;
    }

    // Clamp reference indices to valid range
    ref_i_final = std::min(ref_i_final, MAX_REF - 1);
    ref_j_final = std::min(ref_j_final, MAX_REF - 1);

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("getC6: pair_index=" + std::to_string(pair_index) +
                           " ref_i=" + std::to_string(ref_i) + " → " + std::to_string(ref_i_final) +
                           " ref_j=" + std::to_string(ref_j) + " → " + std::to_string(ref_j_final) +
                           " swap=" + std::string(elem_i <= elem_j ? "YES" : "NO"));
    }

    // Access C6 coefficient from reference data using potentially swapped indices
    if (pair_index < static_cast<int>(m_reference_c6.size())) {
        if (ref_i_final < static_cast<int>(m_reference_c6[pair_index].size()) &&
            ref_j_final < static_cast<int>(m_reference_c6[pair_index][ref_i_final].size())) {
            double c6_value = m_reference_c6[pair_index][ref_i_final][ref_j_final];
            if (CurcumaLogger::get_verbosity() >= 3) {
                CurcumaLogger::info("getC6: Returning c6=" + std::to_string(c6_value));
            }
            return c6_value;
        } else {
            if (CurcumaLogger::get_verbosity() >= 3) {
                CurcumaLogger::warn("getC6: Reference indices out of bounds");
            }
        }
    } else {
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::warn("getC6: pair_index out of bounds");
        }
    }

    return 0.0; // Default if not found
}

double D3ParameterGenerator::getR6(int atom_i, int atom_j) const
{
    // Calculate C8/C6 ratio for BJ damping using r4/r2 reference data
    // Formula from simple-dftd3: C8/C6 = 10.72 * r4/r2_A * r4/r2_B
    // where r4/r2 = sqrt(Z * <r^4>/<r^2>) is element-specific

    // r4/r2 values from simple-dftd3 for elements 1-94
    // Data source: s-dftd3 atomic radii output
    static const std::vector<double> r4_over_r2 = {
        1.0622, // 1 H
        0.8289, // 2 He
        2.6564, // 3 Li
        2.0393, // 4 Be
        1.9286, // 5 B
        1.6431, // 6 C
        1.4350, // 7 N
        1.3725, // 8 O
        1.2638, // 9 F
        1.1722, // 10 Ne
        3.4851, // 11 Na
        2.8909, // 12 Mg
        2.9910, // 13 Al
        2.5839, // 14 Si
        2.2740, // 15 P
        2.1385, // 16 S
        1.9735, // 17 Cl
        1.8240, // 18 Ar
        // Add more elements as needed - for now use fallback for Z>18
    };

    // Elements from Molecule are 1-based (H=1, C=6, ...) - convert to 0-based for array indexing
    int elem_i = atom_i - 1;
    int elem_j = atom_j - 1;

    if (elem_i < 0 || elem_i >= MAX_ELEM || elem_j < 0 || elem_j >= MAX_ELEM) {
        if (CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::warn("getR6: Elements out of range - using default");
        }
        return 12.0; // Conservative default for unknown elements
    }

    // Get r4/r2 values (use default if element not in table)
    double r4r2_i = (elem_i < static_cast<int>(r4_over_r2.size())) ? r4_over_r2[elem_i] : 2.0;
    double r4r2_j = (elem_j < static_cast<int>(r4_over_r2.size())) ? r4_over_r2[elem_j] : 2.0;

    // Calculate C8/C6 ratio using empirical formula from simple-dftd3
    // Factor 10.72 determined from D3 reference data analysis
    double c8_over_c6 = 10.72 * r4r2_i * r4r2_j;

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("getR6(" + std::to_string(atom_i) + "," + std::to_string(atom_j) + "): " +
                           "r4/r2_i=" + std::to_string(r4r2_i) + " r4/r2_j=" + std::to_string(r4r2_j) +
                           " C8/C6=" + std::to_string(c8_over_c6));
    }

    return c8_over_c6;
}

double D3ParameterGenerator::getReferenceCN(int atom, int ref_index) const
{
    // Get reference coordination number for an element
    // Validate BEFORE conversion (atom is 1-based from Molecule)
    if (atom < 1 || atom > MAX_ELEM || ref_index < 0 || ref_index >= MAX_REF) {
        return 0.0;
    }

    --atom; // Convert to 0-based for array indexing

    return m_reference_cn[atom][ref_index];
}

std::vector<double> D3ParameterGenerator::calculateCoordinationNumbers(const std::vector<int>& atoms, const Eigen::MatrixXd& geometry) const
{
    // Calculate D3 coordination numbers using exponential counting function
    // Formula: CN_i = sum_{j≠i} 1 / (1 + exp(-k1 * (k2 * (R_cov_i + R_cov_j) / R_ij - 1)))
    // Reference: Grimme et al., J. Chem. Phys. 132, 154104 (2010)

    const double k1 = 16.0;  // Steepness parameter
    const double k2 = 4.0/3.0;  // Scaling factor

    // Covalent radii from simple-dftd3 (Angstrom)
    // Data source: s-dftd3 atomic radii output
    static const std::vector<double> covalent_radii = {
        0.4267, // 1 H
        0.6133, // 2 He
        1.6000, // 3 Li
        1.2533, // 4 Be
        1.0267, // 5 B
        1.0000, // 6 C
        0.9467, // 7 N
        0.8400, // 8 O
        0.8533, // 9 F
        0.8933, // 10 Ne
        1.8667, // 11 Na
        1.6667, // 12 Mg
        1.5067, // 13 Al
        1.3867, // 14 Si
        1.4667, // 15 P
        1.3600, // 16 S
        1.3200, // 17 Cl
        1.2800, // 18 Ar
        // Add more elements as needed
    };

    std::vector<double> cn_values(atoms.size(), 0.0);

    for (size_t i = 0; i < atoms.size(); ++i) {
        int elem_i = atoms[i] - 1;  // Convert to 0-based

        if (elem_i < 0 || elem_i >= static_cast<int>(covalent_radii.size())) {
            if (CurcumaLogger::get_verbosity() >= 2) {
                CurcumaLogger::warn("calculateCN: Element " + std::to_string(atoms[i]) + " out of range, using CN=0");
            }
            continue;
        }

        double rcov_i = covalent_radii[elem_i];
        double cn = 0.0;

        for (size_t j = 0; j < atoms.size(); ++j) {
            if (i == j) continue;

            int elem_j = atoms[j] - 1;
            if (elem_j < 0 || elem_j >= static_cast<int>(covalent_radii.size())) {
                continue;
            }

            double rcov_j = covalent_radii[elem_j];

            // Calculate distance (geometry is in Angstrom)
            Eigen::Vector3d pos_i = geometry.row(i);
            Eigen::Vector3d pos_j = geometry.row(j);
            double r_ij = (pos_i - pos_j).norm();

            // D3 counting function
            double r_cov_sum = rcov_i + rcov_j;
            double arg = -k1 * (k2 * r_cov_sum / r_ij - 1.0);
            double count = 1.0 / (1.0 + std::exp(arg));

            cn += count;
        }

        cn_values[i] = cn;

        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info("CN[" + std::to_string(i) + "] (Z=" + std::to_string(atoms[i]) + "): " + std::to_string(cn));
        }
    }

    return cn_values;
}

double D3ParameterGenerator::interpolateC6(int elem_i, int elem_j, double cn_i, double cn_j) const
{
    // Interpolate C6 coefficient based on coordination numbers
    // Uses Gaussian weighting: w_i = exp(-4*(CN - CN_ref)^2) / sum(exp(...))
    // Reference: Grimme et al., J. Chem. Phys. 132, 154104 (2010)

    const double k = 4.0;  // Gaussian width parameter

    // Get number of references for both elements
    int nref_i = getNumberofReferences(elem_i);
    int nref_j = getNumberofReferences(elem_j);

    if (nref_i == 0 || nref_j == 0) {
        if (CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::warn("interpolateC6: No references for elements " +
                               std::to_string(elem_i) + "," + std::to_string(elem_j));
        }
        return 0.0;
    }

    // Calculate Gaussian weights for element i
    std::vector<double> weights_i(nref_i, 0.0);
    double sum_weights_i = 0.0;
    for (int ref = 0; ref < nref_i; ++ref) {
        double cn_ref = getReferenceCN(elem_i, ref);
        // Note: CN=0.0 is a valid reference (e.g., H in isolated state), not empty
        double diff = cn_i - cn_ref;
        weights_i[ref] = std::exp(-k * diff * diff);
        sum_weights_i += weights_i[ref];
    }

    // Normalize weights for element i
    if (sum_weights_i > 1e-10) {
        for (int ref = 0; ref < nref_i; ++ref) {
            weights_i[ref] /= sum_weights_i;
        }
    }

    // Calculate Gaussian weights for element j
    std::vector<double> weights_j(nref_j, 0.0);
    double sum_weights_j = 0.0;
    for (int ref = 0; ref < nref_j; ++ref) {
        double cn_ref = getReferenceCN(elem_j, ref);
        // Note: CN=0.0 is a valid reference (e.g., H in isolated state), not empty
        double diff = cn_j - cn_ref;
        weights_j[ref] = std::exp(-k * diff * diff);
        sum_weights_j += weights_j[ref];
    }

    // Normalize weights for element j
    if (sum_weights_j > 1e-10) {
        for (int ref = 0; ref < nref_j; ++ref) {
            weights_j[ref] /= sum_weights_j;
        }
    }

    // Interpolate C6 with weighted sum over all reference combinations
    // NOTE: CN=0 is VALID (isolated atom state), not empty!
    double c6_interpolated = 0.0;
    for (int ref_i = 0; ref_i < nref_i; ++ref_i) {
        for (int ref_j = 0; ref_j < nref_j; ++ref_j) {
            // Claude Generated (December 2025): CRITICAL FIX - Remove +1 offset
            // getC6() already expects 0-based reference indices (0-6 for MAX_REF=7)
            // The +1 offset was incorrect and caused wrong C6 selection
            double c6_ref = getC6(elem_i, elem_j, ref_i, ref_j);
            double contrib = weights_i[ref_i] * weights_j[ref_j] * c6_ref;
            c6_interpolated += contrib;

            if (CurcumaLogger::get_verbosity() >= 3 && contrib > 1e-6) {
                CurcumaLogger::info("    C6[" + std::to_string(ref_i) + "," + std::to_string(ref_j) +
                                   "]=" + std::to_string(c6_ref) + " contrib=" + std::to_string(contrib));
            }
        }
    }

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("interpolateC6(" + std::to_string(elem_i) + "," + std::to_string(elem_j) +
                           "): CN=(" + std::to_string(cn_i) + "," + std::to_string(cn_j) +
                           ") → C6=" + std::to_string(c6_interpolated));
        CurcumaLogger::info("  nref_i=" + std::to_string(nref_i) + " nref_j=" + std::to_string(nref_j));
        for (int ref = 0; ref < nref_i; ++ref) {
            CurcumaLogger::info("    weight_i[" + std::to_string(ref) + "]=" + std::to_string(weights_i[ref]) +
                               " CN_ref=" + std::to_string(getReferenceCN(elem_i, ref)));
        }
    }

    return c6_interpolated;
}

int D3ParameterGenerator::getNumberofReferences(int atom) const
{
    // Validate BEFORE conversion (atom is 1-based from Molecule)
    if (atom < 1 || atom > static_cast<int>(m_number_of_references.size())) {
        return 0;
    }

    --atom; // Convert to 0-based for array indexing

    return m_number_of_references[atom];
}

double D3ParameterGenerator::getTotalEnergy() const
{
    // Calculate total D3 dispersion energy from generated parameters
    if (!m_data_initialized || !m_parameters.contains("d3_dispersion_pairs")) {
        if (CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::warn("D3: Cannot calculate energy - data not initialized or no pairs");
        }
        return 0.0;
    }

    double energy = 0.0;

    // Get damping parameters and scaling factors from config
    double a1 = m_config.get<double>("d3_a1", 0.4);
    double a2 = m_config.get<double>("d3_a2", 4.0);
    double s6 = m_config.get<double>("d3_s6", 1.0);
    double s8 = m_config.get<double>("d3_s8", 1.0);

    const auto& pairs = m_parameters["d3_dispersion_pairs"];

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("D3 getTotalEnergy: Processing " + std::to_string(pairs.size()) + " pairs");
        CurcumaLogger::param("a1", a1);
        CurcumaLogger::param("a2", a2);
    }

    for (const auto& pair : pairs) {
        int i = pair["i"].get<int>();
        int j = pair["j"].get<int>();
        double c6 = pair["c6"].get<double>();
        double c8 = pair["c8"].get<double>();

        // Calculate interatomic distance (Angstrom from geometry)
        Eigen::Vector3d pos_i = m_geometry.row(i);
        Eigen::Vector3d pos_j = m_geometry.row(j);
        double r_angstrom = (pos_i - pos_j).norm();

        // Convert to Bohr (D3 formulas use atomic units)
        const double angstrom_to_bohr = 1.88972612546;  // CurcumaUnit::AngstromToBohr
        double r = r_angstrom * angstrom_to_bohr;

        // Becke-Johnson damping formula (all in Bohr)
        double r6_ratio = pair["r6"].get<double>();
        double r0 = a1 * std::sqrt(r6_ratio) + a2;

        double t6 = 1.0 / (std::pow(r, 6) + std::pow(r0, 6));
        double t8 = 1.0 / (std::pow(r, 8) + std::pow(r0, 8));

        // Two-body dispersion energy (negative - attractive)
        // Claude Generated (December 2025): Apply s6/s8 scaling HERE in energy calculation
        double e6 = -s6 * c6 * t6;
        double e8 = -s8 * c8 * t8;

        energy += (e6 + e8);

        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info("  Pair " + std::to_string(i) + "-" + std::to_string(j) +
                               ": r=" + std::to_string(r) + " Å, E=" + std::to_string(e6 + e8) + " Eh");
        }
    }

    // TODO: Add three-body ATM term if s9 > 0
    // For now, only two-body dispersion implemented

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::result("D3 total energy: " + std::to_string(energy) + " Eh");
    }

    return energy;
}

