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
#include "cn_calculator.h"
#include "../../curcuma_logger.h"

#include <stdexcept>
#include <cctype>
#include <map>

// Complete s-dftd3 reference data (MAX_REF=7 - fixes 1.48x energy error)
// Data split across:
// - test_cases/reference_data/d3_reference_cn.cpp (721 CN values)
// - test_cases/reference_data/d3_reference_c6.cpp (262,444 C6 values)
// Include declarations for automatic linking
#include "d3_reference_data.h"

#include <algorithm>
#include <cmath>

D3ParameterGenerator::D3ParameterGenerator(const ConfigManager& config)
    : m_config(config)
{
    // Phase 3.2: Parameter validation (Claude Generated)
    double s6 = m_config.get<double>("d3_s6", 1.0);
    double s8 = m_config.get<double>("d3_s8", 1.0);
    double a1 = m_config.get<double>("d3_a1", 0.4);
    double a2 = m_config.get<double>("d3_a2", 4.0);

    // Validate scaling factors (must be non-negative)
    if (s6 < 0.0 || s8 < 0.0) {
        throw std::invalid_argument(
            "D3 scaling factors must be non-negative: s6=" + std::to_string(s6) +
            ", s8=" + std::to_string(s8)
        );
    }

    // Warn on unusual parameters
    if (s6 > 2.0 || s8 > 5.0) {
        if (CurcumaLogger::get_verbosity() >= 1) {
            CurcumaLogger::warn(
                "Unusual D3 scaling factors: s6=" + std::to_string(s6) +
                " (typical 0.75-1.5), s8=" + std::to_string(s8) + " (typical 0.5-3.0)"
            );
        }
    }

    // Validate damping parameters
    if (a1 < 0.0 || a2 < 0.0) {
        throw std::invalid_argument(
            "D3 damping parameters must be non-negative: a1=" + std::to_string(a1) +
            ", a2=" + std::to_string(a2)
        );
    }

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
    auto t_start_total = std::chrono::high_resolution_clock::now();

    m_atoms = atoms;
    m_geometry = geometry;
    m_parameters.clear();

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("=== D3ParameterGenerator::GenerateParameters() START ===");
        CurcumaLogger::param("Number of atoms", static_cast<int>(m_atoms.size()));
    }

    // Calculate geometry-dependent coordination numbers
    auto t_cn_start = std::chrono::high_resolution_clock::now();
    std::vector<double> coordination_numbers = calculateCoordinationNumbers(m_atoms, m_geometry);
    auto t_cn_end = std::chrono::high_resolution_clock::now();
    double t_cn_ms = std::chrono::duration<double, std::milli>(t_cn_end - t_cn_start).count();

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info(fmt::format("D3: CN calculation took {:.2f} ms", t_cn_ms));
    }

    // Claude Generated (Dec 2025): Pre-compute Gaussian weights ONCE for all atoms
    // This eliminates redundant exp() calls in interpolateC6()
    // Performance: O(N×M) instead of O(N²×M) exp() calculations
    auto t_weights_start = std::chrono::high_resolution_clock::now();
    precomputeGaussianWeights(m_atoms, coordination_numbers);
    auto t_weights_end = std::chrono::high_resolution_clock::now();
    double t_weights_ms = std::chrono::duration<double, std::milli>(t_weights_end - t_weights_start).count();

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info(fmt::format("D3: Weight pre-computation took {:.2f} ms ({} atoms × avg {} refs)",
            t_weights_ms, m_atoms.size(), m_gaussian_weights.empty() ? 0 : m_gaussian_weights[0].size()));
    }

    // Generate dispersion parameters for all atom pairs
    auto t_pairs_start = std::chrono::high_resolution_clock::now();
    json dispersion_pairs = json::array();

    // Get global scaling factors from config
    double s6 = m_config.get<double>("d3_s6", 1.0);
    double s8 = m_config.get<double>("d3_s8", 1.0);
    double a1 = m_config.get<double>("d3_a1", 0.4);
    double a2 = m_config.get<double>("d3_a2", 4.0);

    int num_pairs = 0;
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
                // Claude Generated (Dec 2025): Use cached weights - pass atom indices, not CN values
                double c6 = interpolateC6(atom_i, atom_j, i, j);  // Pass atom indices for weight lookup
                pair["c6"] = c6;  // Store raw C6, s6 applied in energy calculation
                pair["cn_i"] = coordination_numbers[i];  // Store CN for debugging
                pair["cn_j"] = coordination_numbers[j];

                // Calculate C8 coefficient using C8/C6 ratio
                // Claude Generated (December 2025): CRITICAL FIX - s8 applied in energy, not here!
                double c8_over_c6 = getR6(atom_i, atom_j);
                pair["c8"] = c6 * c8_over_c6;  // C8 = C6 * (C8/C6 ratio), s8 applied in energy calculation
                pair["r6"] = c8_over_c6;

                // Add van der Waals radius for distance cutoff
                pair["vdw_radius"] = getR6(atom_i, atom_j); // Use r6 as fallback

                // Add missing fields expected by setGFNFFDispersions
                pair["C6"] = c6;  // Duplicate for compatibility
                pair["C8"] = c6 * c8_over_c6;  // Duplicate for compatibility
                pair["s6"] = s6;
                pair["s8"] = s8;
                pair["a1"] = a1;
                pair["a2"] = a2;
                // Use a reasonable default for r_cut (100 Bohr = ~52.9 Å)
                pair["r_cut"] = 100.0;

                dispersion_pairs.push_back(pair);
                num_pairs++;
            } else {
                if (CurcumaLogger::get_verbosity() >= 2) {
                    CurcumaLogger::warn("D3: Elements out of range - i=" +
                                       std::to_string(atom_i) +
                                       " j=" + std::to_string(atom_j));
                }
            }
        }
    }

    auto t_pairs_end = std::chrono::high_resolution_clock::now();
    double t_pairs_ms = std::chrono::duration<double, std::milli>(t_pairs_end - t_pairs_start).count();

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info(fmt::format("D3: C6 interpolation for {} pairs took {:.2f} ms ({:.4f} ms/pair)",
            num_pairs, t_pairs_ms, num_pairs > 0 ? t_pairs_ms / num_pairs : 0.0));
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

    // Phase 2.3 (December 2025): ATM three-body dispersion
    // Reference: external/cpp-d4/src/damping/atm.cpp:70-138
    json atm_triples = json::array();
    double t_atm_triples_ms = 0.0;  // Initialize for scope (used in final breakdown)

    if (m_config.get<double>("d3_s9", 0.0) > 1e-10) {
        double s9 = m_config.get<double>("d3_s9", 1.0);
        double a1 = m_config.get<double>("d3_a1", 0.40);
        double a2 = m_config.get<double>("d3_a2", 4.20);
        double alp = m_config.get<double>("d3_alp", 14.0);

        int n_atoms = static_cast<int>(m_atoms.size());

        // Claude Generated (January 2026): Build O(log N) C6 lookup map for ATM triples
        // Replaces O(N²) linear search per triple, improving performance from O(N⁶) to O(N³ log N)
        std::map<std::pair<int,int>, double> c6_lookup;

        for (const auto& pair : dispersion_pairs) {
            int i = pair["i"];
            int j = pair["j"];
            double c6 = pair["c6"];

            // Store both (i,j) and (j,i) for O(1) symmetric lookup
            c6_lookup[{i, j}] = c6;
            c6_lookup[{j, i}] = c6;
        }

        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::param("C6 lookup map size", static_cast<int>(c6_lookup.size()));
        }

        auto t_atm_triples_start = std::chrono::high_resolution_clock::now();

        for (int i = 0; i < n_atoms; ++i) {
            for (int j = 0; j < i; ++j) {
                for (int k = 0; k < j; ++k) {
                    json triple;
                    triple["i"] = i;
                    triple["j"] = j;
                    triple["k"] = k;

                    // Claude Generated (January 2026): O(log N) C6 lookup from pre-built map
                    // Extract C6 values from lookup map (fast vs. O(N²) linear search)
                    double c6_ij = c6_lookup[{i, j}];
                    double c6_ik = c6_lookup[{i, k}];
                    double c6_jk = c6_lookup[{j, k}];

                    triple["C6_ij"] = c6_ij;
                    triple["C6_ik"] = c6_ik;
                    triple["C6_jk"] = c6_jk;
                    triple["s9"] = s9;
                    triple["a1"] = a1;
                    triple["a2"] = a2;
                    triple["alp"] = alp;
                    triple["atm_method"] = "d3";

                    // Symmetry factor
                    triple["triple_scale"] = calculateTripleScale(i, j, k);

                    atm_triples.push_back(triple);
                }
            }
        }

        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::param("Generated D3 ATM triples", static_cast<int>(atm_triples.size()));
        }

        auto t_atm_triples_end = std::chrono::high_resolution_clock::now();
        double t_atm_triples_ms = std::chrono::duration<double, std::milli>(t_atm_triples_end - t_atm_triples_start).count();

        if (CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::info(fmt::format("D3: ATM triple generation took {:.2f} ms ({} triples)",
                t_atm_triples_ms, static_cast<int>(atm_triples.size())));
        }
    }

    m_parameters["atm_triples"] = atm_triples;

    auto t_end_total = std::chrono::high_resolution_clock::now();
    double t_total_ms = std::chrono::duration<double, std::milli>(t_end_total - t_start_total).count();

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::success(fmt::format("D3 parameter generation completed in {:.2f} ms", t_total_ms));
        CurcumaLogger::info(fmt::format("  └─ Breakdown: CN={:.1f}ms, Weights={:.1f}ms, C6Interp={:.1f}ms, ATM={:.1f}ms",
            t_cn_ms, t_weights_ms, t_pairs_ms, t_atm_triples_ms));
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
    // C8/C6 ratio for BJ damping — exact s-dftd3 / dftd4 form:
    //   rrij = C8/C6 = 3 * r4r2(i) * r4r2(j)            (s-dftd3 damping/rational.f90:173)
    //   r4r2(z) = sqrt(0.5 * <r^4>/<r^2>(z) * sqrt(z))  (s-dftd3 data/r4r2.f90:69-70)
    // The raw <r^4>/<r^2> table below is verbatim from s-dftd3 data/r4r2.f90
    // (PBE0/def2-QZVP, Grimme 2010; identical to D4ParameterGenerator::m_r4_over_r2).
    // Earlier code used an empirical "10.72 * r4r2" fit with a non-standard table
    // that was ~0.06% high on every pair — a size-extensive dispersion bias vs
    // tblite/s-dftd3. This reproduces s-dftd3's C8/C6 (and BJ R0) exactly.
    static const double raw_r4_over_r2[] = {
        8.0589,  3.4698,                                                                  // H,He
        29.0974, 14.8517, 11.8799, 7.8715, 5.5588, 4.7566, 3.8025, 3.1036,               // Li-Ne
        26.1552, 17.2304, 17.7210, 12.7442, 9.5361, 8.1652, 6.7463, 5.6004,              // Na-Ar
        29.2012, 22.3934,                                                                // K,Ca
        19.0598, 16.8590, 15.4023, 12.5589, 13.4788,                                     // Sc-Mn
        12.2309, 11.2809, 10.5569, 10.1428, 9.4907,                                      // Fe-Zn
        13.4606, 10.8544, 8.9386, 8.1350, 7.1251, 6.1971,                                // Ga-Kr
        30.0162, 24.4103,                                                                // Rb,Sr
        20.3537, 17.4780, 13.5528, 11.8451, 11.0355, 10.1997, 9.5414, 9.0061, 8.6417, 8.9975, // Y-Cd
        14.0834, 11.8333, 10.0179, 9.3844, 8.4110, 7.5152,                               // In-Xe
        32.7622, 27.5708,                                                                // Cs,Ba
        23.1671, 21.6003, 20.9615, 20.4562, 20.1010, 19.7475, 19.4828,                   // La-Eu
        15.6013, 19.2362, 17.4717, 17.8321, 17.4237, 17.1954, 17.1631,                   // Gd-Yb
        14.5716, 15.8758, 13.8989, 12.4834, 11.4421, 10.2671, 8.3549, 7.8496, 7.3278, 7.4820, // Lu-Hg
        13.5124, 11.6554, 10.0959, 9.7340, 8.8584, 8.0125,                               // Tl-Rn
        29.8135, 26.3157,                                                                // Fr,Ra
        19.1885, 15.8542, 16.1305, 15.6161, 15.1226, 16.1576, 14.6510,                   // Ac-Am
        14.7178, 13.9108, 13.5623, 13.2326, 12.9189, 12.6133, 12.3142,                   // Cm-No
        14.8326, 12.3771, 10.6378, 9.3638, 8.2297, 7.5667, 6.9456, 6.3946, 5.9159, 5.4929, // Lr-Cn
        6.7286, 6.5144, 10.9169, 10.3600, 9.4723, 8.6641                                 // Nh-Og
    };
    static const int n_r4r2 = static_cast<int>(sizeof(raw_r4_over_r2) / sizeof(double));

    // Elements from Molecule are 1-based (H=1, C=6, ...) - convert to 0-based for array indexing
    const int elem_i = atom_i - 1;
    const int elem_j = atom_j - 1;
    if (elem_i < 0 || elem_j < 0 || elem_i >= n_r4r2 || elem_j >= n_r4r2) {
        if (CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::warn("getR6: Elements out of range - using default");
        }
        return 12.0; // Conservative default for unknown elements
    }

    // Pre-scale: sqrtZr4r2(z) = sqrt(0.5 * <r^4>/<r^2> * sqrt(Z)).
    const double szr_i = std::sqrt(0.5 * raw_r4_over_r2[elem_i] * std::sqrt(static_cast<double>(atom_i)));
    const double szr_j = std::sqrt(0.5 * raw_r4_over_r2[elem_j] * std::sqrt(static_cast<double>(atom_j)));
    const double c8_over_c6 = 3.0 * szr_i * szr_j;

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("getR6(" + std::to_string(atom_i) + "," + std::to_string(atom_j) + "): " +
                           "sqrtZr4r2_i=" + std::to_string(szr_i) + " sqrtZr4r2_j=" + std::to_string(szr_j) +
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
    // Delegate to shared CNCalculator utility
    // Claude Generated (December 20, 2025): Consolidated CN calculation
    return CNCalculator::calculateD3CN(atoms, geometry);
}

void D3ParameterGenerator::precomputeGaussianWeights(
    const std::vector<int>& atoms,
    const std::vector<double>& coordination_numbers)
{
    // Claude Generated (Dec 2025): Pre-compute Gaussian weights once per atom
    // Reference: simple-dftd3/src/dftd3/model.f90:147-235 (weight_references)
    //
    // Performance optimization: Eliminates redundant exp() calls in interpolateC6()
    // - Old approach: O(N²×M) exp() calls (compute weights for every atom pair)
    // - New approach: O(N×M) exp() calls (compute weights once per atom)
    // Expected speedup: 5-10x for molecules with 100+ atoms

    const double k = 4.0;  // Gaussian width parameter (matches simple-dftd3)
    m_gaussian_weights.resize(atoms.size());
    m_cached_cn = coordination_numbers;

    for (size_t i = 0; i < atoms.size(); ++i) {
        int elem = atoms[i];
        int nref = getNumberofReferences(elem);

        if (nref == 0) {
            m_gaussian_weights[i].clear();
            continue;
        }

        // Compute Gaussian weights for all reference states of this atom
        std::vector<double> weights(nref, 0.0);
        double sum_weights = 0.0;

        for (int ref = 0; ref < nref; ++ref) {
            double cn_ref = getReferenceCN(elem, ref);
            double diff = coordination_numbers[i] - cn_ref;
            weights[ref] = std::exp(-k * diff * diff);  // COMPUTED ONCE PER ATOM
            sum_weights += weights[ref];
        }

        // Normalize weights
        if (sum_weights > 1e-10) {
            for (int ref = 0; ref < nref; ++ref) {
                weights[ref] /= sum_weights;
            }
        } else {
            // Exceptional case: all weights negligible
            // Set highest CN reference to 1.0 (matches Fortran behavior)
            double max_cn_ref = -1.0;
            int max_ref_idx = 0;
            for (int ref = 0; ref < nref; ++ref) {
                double cn_ref = getReferenceCN(elem, ref);
                if (cn_ref > max_cn_ref) {
                    max_cn_ref = cn_ref;
                    max_ref_idx = ref;
                }
            }
            std::fill(weights.begin(), weights.end(), 0.0);
            weights[max_ref_idx] = 1.0;
        }

        m_gaussian_weights[i] = std::move(weights);
    }

    m_weights_cached = true;

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("=== D3 Gaussian Weights Pre-computed ===");
        CurcumaLogger::param("Atoms processed", static_cast<int>(atoms.size()));
        for (size_t i = 0; i < atoms.size(); ++i) {
            std::string weights_str = "";
            for (size_t ref = 0; ref < m_gaussian_weights[i].size(); ++ref) {
                weights_str += std::to_string(m_gaussian_weights[i][ref]);
                if (ref < m_gaussian_weights[i].size() - 1) weights_str += ", ";
            }
            CurcumaLogger::info("  Atom " + std::to_string(i) + " (elem=" + std::to_string(atoms[i]) +
                               ", CN=" + std::to_string(coordination_numbers[i]) +
                               "): weights=[" + weights_str + "]");
        }
    }
}

double D3ParameterGenerator::interpolateC6(int elem_i, int elem_j, size_t atom_idx_i, size_t atom_idx_j) const
{
    // Claude Generated (Dec 2025): Optimized C6 interpolation using pre-computed weights
    // Reference: simple-dftd3/src/dftd3/model.f90:248-324 (get_atomic_c6)
    //
    // Performance: NO EXP() CALLS - uses cached Gaussian weights from precomputeGaussianWeights()
    // This function is called O(N²) times, but exp() is only computed O(N) times total

    if (!m_weights_cached) {
        CurcumaLogger::error("interpolateC6: Weights not cached! Call precomputeGaussianWeights() first.");
        return 0.0;
    }

    // Lookup pre-computed weights (NO COMPUTATION, only array access)
    const auto& weights_i = m_gaussian_weights[atom_idx_i];
    const auto& weights_j = m_gaussian_weights[atom_idx_j];

    if (weights_i.empty() || weights_j.empty()) {
        if (CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::warn("interpolateC6: Empty weights for atoms " +
                               std::to_string(atom_idx_i) + "," + std::to_string(atom_idx_j));
        }
        return 0.0;
    }

    // Weighted sum over reference combinations (NO EXP CALLS, only multiplications)
    double c6_interpolated = 0.0;
    for (size_t ref_i = 0; ref_i < weights_i.size(); ++ref_i) {
        for (size_t ref_j = 0; ref_j < weights_j.size(); ++ref_j) {
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
                           "): atoms=(" + std::to_string(atom_idx_i) + "," + std::to_string(atom_idx_j) +
                           "), CN=(" + std::to_string(m_cached_cn[atom_idx_i]) + "," +
                           std::to_string(m_cached_cn[atom_idx_j]) +
                           ") → C6=" + std::to_string(c6_interpolated));
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

// Claude Generated Phase 3.1: Factory Methods
// Provide convenient creation with standard D3 parameter sets

D3ParameterGenerator D3ParameterGenerator::createForGFN1()
{
    // Claude Generated: GFN1-xTB D3(BJ) parameters
    // Reference: tblite src/tblite/xtb/gfn1.f90:53 (Grimme et al. JCTC 2017, 13, 1989)
    //   s6 = 1.0, s8 = 2.4, a1 = 0.63, a2 = 5.0, s9 = 0.0
    // s9 = 0.0 pins GFN1-xTB's D3 to TWO-BODY ONLY (no Axilrod-Teller-Muto term),
    // matching tblite. This call already behaved as s9=0 (the PARAM default 1.0 is
    // not merged into this locally-built config), so this is an explicit guard, not
    // a behaviour change -- it prevents a future registry-default leak from silently
    // switching ATM on for GFN1. (The GFN1 dispersion residual vs tblite is in the
    // TWO-body term, not ATM; see SQM_WP2.)
    json config_json;
    config_json["d3_s6"] = 1.0;
    config_json["d3_s8"] = 2.4;
    config_json["d3_a1"] = 0.63;
    config_json["d3_a2"] = 5.0;
    config_json["d3_s9"] = 0.0;   // no three-body ATM (tblite gfn1.f90:53)
    config_json["d3_alp"] = 14.0;

    ConfigManager config("d3param", config_json);
    return D3ParameterGenerator(config);
}

D3ParameterGenerator D3ParameterGenerator::createForGFNFF()
{
    // GFN-FF D3-BJ parameters (Spicher & Grimme, J. Chem. Theory Comput. 2020)
    json config_json;
    config_json["d3_s6"] = 1.0;
    config_json["d3_s8"] = 2.85;
    config_json["d3_a1"] = 0.80;
    config_json["d3_a2"] = 4.60;
    config_json["d3_alp"] = 14.0;

    ConfigManager config("d3param", config_json);
    return D3ParameterGenerator(config);
}

D3ParameterGenerator D3ParameterGenerator::createForUFFD3()
{
    // PBE0-D3-BJ parameters (UFF bonded + D3 dispersion)
    json config_json;
    config_json["d3_s6"] = 1.0;
    config_json["d3_s8"] = 1.2177;
    config_json["d3_a1"] = 0.4145;
    config_json["d3_a2"] = 4.8593;
    config_json["d3_alp"] = 14.0;

    ConfigManager config("d3param", config_json);
    return D3ParameterGenerator(config);
}

D3ParameterGenerator D3ParameterGenerator::createForPBE0()
{
    // Standard PBE0-D3-BJ parameters
    // Same as UFF-D3 (PBE0 functional with D3-BJ damping)
    return createForUFFD3();
}

D3ParameterGenerator D3ParameterGenerator::createForBLYP()
{
    // BLYP-D3-BJ parameters (Grimme et al., J. Chem. Phys. 132, 154104 2010, Table II)
    json config_json;
    config_json["d3_s6"] = 1.0;
    config_json["d3_s8"] = 2.6996;
    config_json["d3_a1"] = 0.4298;
    config_json["d3_a2"] = 4.2359;
    config_json["d3_alp"] = 14.0;

    ConfigManager config("d3param", config_json);
    return D3ParameterGenerator(config);
}

D3ParameterGenerator D3ParameterGenerator::createForB3LYP()
{
    // B3LYP-D3-BJ parameters (Grimme et al., J. Chem. Phys. 132, 154104 2010, Table II)
    json config_json;
    config_json["d3_s6"] = 1.0;
    config_json["d3_s8"] = 1.9889;
    config_json["d3_a1"] = 0.3981;
    config_json["d3_a2"] = 4.4211;
    config_json["d3_alp"] = 14.0;

    ConfigManager config("d3param", config_json);
    return D3ParameterGenerator(config);
}

D3ParameterGenerator D3ParameterGenerator::createForTPSS()
{
    // TPSS-D3-BJ parameters (Grimme et al., J. Chem. Phys. 132, 154104 2010, Table II)
    json config_json;
    config_json["d3_s6"] = 1.0;
    config_json["d3_s8"] = 1.9435;
    config_json["d3_a1"] = 0.4535;
    config_json["d3_a2"] = 4.4752;
    config_json["d3_alp"] = 14.0;

    ConfigManager config("d3param", config_json);
    return D3ParameterGenerator(config);
}

D3ParameterGenerator D3ParameterGenerator::createForPBE()
{
    // PBE-D3-BJ parameters (Grimme et al., J. Chem. Phys. 132, 154104 2010, Table II)
    json config_json;
    config_json["d3_s6"] = 1.0;
    config_json["d3_s8"] = 0.7875;
    config_json["d3_a1"] = 0.4289;
    config_json["d3_a2"] = 4.4407;
    config_json["d3_alp"] = 14.0;

    ConfigManager config("d3param", config_json);
    return D3ParameterGenerator(config);
}

D3ParameterGenerator D3ParameterGenerator::createForBP86()
{
    // BP86-D3-BJ parameters (Grimme et al., J. Chem. Phys. 132, 154104 2010, Table II)
    json config_json;
    config_json["d3_s6"] = 1.0;
    config_json["d3_s8"] = 3.2822;
    config_json["d3_a1"] = 0.3946;
    config_json["d3_a2"] = 4.8516;
    config_json["d3_alp"] = 14.0;

    ConfigManager config("d3param", config_json);
    return D3ParameterGenerator(config);
}

D3ParameterGenerator D3ParameterGenerator::createForMethod(const std::string& method)
{
    std::string lower_method = method;
    std::transform(lower_method.begin(), lower_method.end(), lower_method.begin(), ::tolower);

    if (lower_method == "gfn1" || lower_method == "gfn1-xtb") {
        return createForGFN1();
    } else if (lower_method == "gfnff") {
        return createForGFNFF();
    } else if (lower_method == "uff-d3" || lower_method == "uffd3") {
        return createForUFFD3();
    } else if (lower_method == "pbe0") {
        return createForPBE0();
    } else if (lower_method == "blyp") {
        return createForBLYP();
    } else if (lower_method == "b3lyp") {
        return createForB3LYP();
    } else if (lower_method == "tpss") {
        return createForTPSS();
    } else if (lower_method == "pbe") {
        return createForPBE();
    } else if (lower_method == "bp86") {
        return createForBP86();
    } else {
        throw std::invalid_argument(
            "Unknown D3 preset: " + method +
            ". Supported: pbe0, blyp, b3lyp, tpss, pbe, bp86, gfnff, uff-d3"
        );
    }
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

// Claude Generated (May 2026): Analytical gradient for D3 dispersion
// -------------------------------------------------------------------------

void D3ParameterGenerator::computeDC6DCN() const
{
    // Reference: s-dftd3/src/dftd3/model.f90:248-324 (get_atomic_c6 derivatives)
    // Build the matrix m_dc6dcn where:
    //   m_dc6dcn(i,j) = dC6(i,j) / dCN(i)
    //   m_dc6dcn(j,i) = dC6(i,j) / dCN(j)
    // This is needed for the CN chain rule in the analytical gradient.

    if (m_dc6dcn_computed || !m_weights_cached)
        return;

    const size_t n = m_atoms.size();
    m_dc6dcn = Eigen::MatrixXd::Zero(n, n);

    const double wf = 4.0;  // Gaussian width parameter (same as precomputeGaussianWeights)

    // Pre-compute dgw/dcn for each atom's Gaussian weights
    // gw(ref) = exp(-wf*(cn - cn_ref)^2) / norm
    // d_gw/dcn = gw(ref) * 2*wf*(cn_ref - cn) / norm  [unnormalised raw derivative]
    // Actually: d/dcn[exp(-wf*(cn-cn_ref)^2)] = -2*wf*(cn-cn_ref)*exp(...)
    // Normalised derivative: d(gw/norm)/dcn = (d_gw*norm - gw*sum(d_gw)) / norm^2
    std::vector<std::vector<double>> dgw(n);

    for (size_t idx = 0; idx < n; ++idx) {
        int elem = m_atoms[idx];
        int nref = getNumberofReferences(elem);
        if (nref == 0) continue;

        double cn = m_cached_cn[idx];
        const auto& gw = m_gaussian_weights[idx];

        // Un-normalised weights and their derivatives
        std::vector<double> raw_gw(nref);
        std::vector<double> raw_dgw(nref);
        double norm = 0.0;
        double dnorm = 0.0;

        for (int ref = 0; ref < nref; ++ref) {
            double cn_ref = getReferenceCN(elem, ref);
            double diff = cn - cn_ref;
            raw_gw[ref] = std::exp(-wf * diff * diff);
            raw_dgw[ref] = -2.0 * wf * diff * raw_gw[ref];
            norm += raw_gw[ref];
            dnorm += raw_dgw[ref];
        }

        if (norm < 1e-10) {
            // Exceptional case: all weights negligible -> delta on max CN ref
            dgw[idx].resize(nref, 0.0);
            continue;
        }

        dgw[idx].resize(nref);
        for (int ref = 0; ref < nref; ++ref) {
            dgw[idx][ref] = (raw_dgw[ref] * norm - raw_gw[ref] * dnorm) / (norm * norm);
        }
    }

    // Now build dc6dcn using the weight derivatives
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            int elem_i = m_atoms[i];
            int elem_j = m_atoms[j];

            if (elem_i < 1 || elem_i > MAX_ELEM || elem_j < 1 || elem_j > MAX_ELEM)
                continue;

            const auto& gw_i = m_gaussian_weights[i];
            const auto& gw_j = m_gaussian_weights[j];
            const auto& dgw_i = dgw[i];
            const auto& dgw_j = dgw[j];

            if (gw_i.empty() || gw_j.empty()) continue;

            double dc6_di = 0.0;
            double dc6_dj = 0.0;

            for (size_t ref_i = 0; ref_i < gw_i.size(); ++ref_i) {
                for (size_t ref_j = 0; ref_j < gw_j.size(); ++ref_j) {
                    double c6_ref = getC6(elem_i, elem_j, ref_i, ref_j);
                    dc6_di += dgw_i[ref_i] * gw_j[ref_j] * c6_ref;
                    dc6_dj += gw_i[ref_i] * dgw_j[ref_j] * c6_ref;
                }
            }

            m_dc6dcn(i, j) = dc6_di;
            m_dc6dcn(j, i) = dc6_dj;
        }
    }

    m_dc6dcn_computed = true;
}

const Eigen::MatrixXd& D3ParameterGenerator::getDC6DCN() const
{
    if (!m_dc6dcn_computed)
        computeDC6DCN();
    return m_dc6dcn;
}

double D3ParameterGenerator::getEnergyAndGradient(bool need_gradient,
                                                  Matrix& gradient_out,
                                                  Vector& dEdcn_out) const
{
    // Reference: s-dftd3/src/dftd3/damping/rational.f90 (get_dispersion_derivs)
    // Computes D3(BJ) two-body dispersion energy and analytical gradient.
    //
    // The direct geometry gradient is folded into gradient_out immediately.
    // The CN-chain-rule term dE/dCN is returned in dEdcn_out; the caller must
    // multiply by dCN/dR (exponential counting function) and add to gradient.

    if (!m_data_initialized || !m_weights_cached || m_atoms.empty()) {
        return 0.0;
    }

    const double au = 1.88972612546;  // Angstrom -> Bohr
    const size_t n = m_atoms.size();

    double a1 = m_config.get<double>("d3_a1", 0.4);
    double a2 = m_config.get<double>("d3_a2", 4.0);
    double s6 = m_config.get<double>("d3_s6", 1.0);
    double s8 = m_config.get<double>("d3_s8", 1.0);

    double energy = 0.0;

    if (need_gradient) {
        computeDC6DCN();
        // Ensure outputs are sized (caller may or may not have pre-allocated)
        if (gradient_out.rows() != static_cast<int>(n) || gradient_out.cols() != 3)
            gradient_out = Eigen::MatrixXd::Zero(n, 3);
        if (dEdcn_out.size() != static_cast<int>(n))
            dEdcn_out = Eigen::VectorXd::Zero(n);
    }

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            int elem_i = m_atoms[i];
            int elem_j = m_atoms[j];

            if (elem_i < 1 || elem_i > MAX_ELEM || elem_j < 1 || elem_j > MAX_ELEM)
                continue;

            double c6 = interpolateC6(elem_i, elem_j, i, j);
            if (c6 <= 0.0) continue;

            double r6_ratio = getR6(elem_i, elem_j);
            double c8 = c6 * r6_ratio;

            Eigen::Vector3d pos_i = m_geometry.row(i);
            Eigen::Vector3d pos_j = m_geometry.row(j);
            Eigen::Vector3d rij_vec = pos_i - pos_j;          // Angstrom
            double r_angstrom = rij_vec.norm();
            if (r_angstrom < 1e-6) continue;

            double r = r_angstrom * au;                       // Bohr
            double r2 = r * r;

            // Becke-Johnson damping
            double r0 = a1 * std::sqrt(r6_ratio) + a2;        // Bohr

            double r6_0 = std::pow(r0, 6);
            double r8_0 = std::pow(r0, 8);
            double t6 = 1.0 / (std::pow(r, 6) + r6_0);
            double t8 = 1.0 / (std::pow(r, 8) + r8_0);

            // Two-body energy
            double e6 = -s6 * c6 * t6;
            double e8 = -s8 * c8 * t8;
            energy += (e6 + e8);

            if (!need_gradient) continue;

            // --- Direct geometry gradient (dE/dR) ---
            // d/dr [t6] = -6*r^5 * t6^2  ->  d6_contrib = -c6 * s6 * (-6*r^5*t6^2)
            // In s-dftd3 code: d6 = -6*r2^2*t6^2  (r2 = r^2)
            //   -> gdisp = s6*d6 + s8*rrij*d8
            //   -> dG = -c6ij * gdisp * vec   (vec = R_i - R_j in Bohr)
            //
            // But careful: vec in s-dftd3 is in Bohr, and dG is Eh/Bohr.
            // Our rij_vec is in Angstrom. We must convert to Bohr for consistency.

            double d6 = -6.0 * r2 * r2 * t6 * t6;   // = -6*r^4*t6^2
            double d8 = -8.0 * r2 * r2 * r2 * t8 * t8; // = -8*r^6*t8^2
            double gdisp = s6 * d6 + s8 * r6_ratio * d8;

            // rij in Bohr for gradient accumulation
            Eigen::Vector3d vec_bohr = rij_vec * au;
            Eigen::Vector3d dG = -c6 * gdisp * vec_bohr;   // Eh/Bohr

            gradient_out.row(i) += dG.transpose();
            gradient_out.row(j) -= dG.transpose();

            // --- CN chain-rule term (dE/dCN) ---
            // edisp = s6*t6 + s8*rrij*t8   (factor common to both dE/dCN terms)
            double edisp = s6 * t6 + s8 * r6_ratio * t8;
            dEdcn_out(i) -= m_dc6dcn(i, j) * edisp;
            dEdcn_out(j) -= m_dc6dcn(j, i) * edisp;
        }
    }

    return energy;
}

// Claude Generated (2025): ATM three-body symmetry factor calculation
double D3ParameterGenerator::calculateTripleScale(int i, int j, int k) const
{
    // Reference: external/cpp-d4/src/damping/atm.cpp:291-313
    if (i == j) {
        return (i == k) ? 1.0/6.0 : 0.5;  // iii: 1/6, iij: 1/2
    } else {
        return (i != k && j != k) ? 1.0 : 0.5;  // ijk: 1, ijj/iji: 1/2
    }
}
