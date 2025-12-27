/*
 * DFT-D4 Parameter Generator for Curcuma - Complete Data Integration
 * Integrates 118 elements of D4 reference data from GFN-FF Fortran implementation
 *
 * Claude Generated (December 2025): Phase 2.1 - D4 Reference Data Integration
 * December 25, 2025: Complete alphaiw (polarizability) data integration
 */

#include "d4param_generator.h"
#include "../../../../test_cases/reference_data/d4_reference_data_fixed.cpp"  // D4 reference data (365 lines)
#include "../../../../test_cases/reference_data/d4_reference_cn.cpp"          // D4 reference CN data (cpp-d4 - December 2025 Phase 1)
#include "../../../../test_cases/reference_data/d4_alphaiw_data.cpp"         // D4 alphaiw data (269 reference states)
#include "../../../../test_cases/reference_data/d4_corrections_data.cpp"     // D4 correction factors
#include "src/core/curcuma_logger.h"

#include <algorithm>
#include <cmath>
#include <limits>

// External declarations from d4_reference_cn.cpp (cpp-d4 CN data - December 2025 Phase 1)
extern const std::vector<std::vector<double>> D4ReferenceData::refcn_cppd4;

// External declarations from d4_alphaiw_data.cpp
extern std::vector<std::vector<std::vector<double>>> d4_alphaiw_data;
extern void initialize_d4_alphaiw();

// External declarations from d4_corrections_data.cpp
extern std::vector<std::vector<double>> d4_ascale_data;
extern std::map<int, double> d4_sscale_data;
extern std::map<int, std::vector<double>> d4_secaiw_data;
extern std::vector<std::vector<int>> d4_refsys_data;
extern void initialize_d4_corrections();

D4ParameterGenerator::D4ParameterGenerator(const ConfigManager& config)
    : m_config(config)
{
    // Initialize EEQ solver for charge calculation (Dec 2025 - Phase 2)
    ConfigManager eeq_config("eeq_solver", config.exportConfig());
    m_eeq_solver = std::make_unique<EEQSolver>(eeq_config);

    initializeReferenceData();

    // Initialize complete alphaiw data (Dec 25, 2025 - Phase 2.3)
    // 269 reference states with frequency-dependent polarizabilities
    initialize_d4_alphaiw();

    // Initialize correction factors (Dec 25, 2025 - Phase 2.4)
    // ascale, sscale, secaiw, refsys for accurate polarizability corrections
    initialize_d4_corrections();

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::success("D4: Complete data loaded (alphaiw + corrections)");
    }
}

void D4ParameterGenerator::initializeReferenceData()
{
    // Phase 2.1 (December 2025): D4 reference data integration
    // Import complete D4 reference data from d4_reference_data_fixed.cpp (365 lines)
    // - 118 elements (H through Og)
    // - Extracted from external/gfnff/src/dftd4param.f90 via extract_d4_data.py

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("=== D4ParameterGenerator::initializeReferenceData() ===");
    }

    // Load reference data from d4_reference_data_fixed.cpp
    m_refn = ::d4_refn;  // Number of reference states per element (118 elements)
    m_refq = ::d4_refq;  // Reference charges (118 × 7)
    m_refh = ::d4_refh;  // Reference hydrogen counts (118 × 7)

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::success(fmt::format("D4: Loaded {} elements from reference data", m_refn.size()));

        // Print first 5 elements as validation
        CurcumaLogger::info("First 5 elements validation:");
        for (int i = 0; i < 5 && i < static_cast<int>(m_refn.size()); ++i) {
            CurcumaLogger::result(fmt::format("  Element {} (Z={}): {} reference states, q[0]={:.6f}",
                                              i+1, i+1, m_refn[i], m_refq[i][0]));
        }
    }

    // Load reference coordination numbers (December 2025 Phase 2 - CN integration)
    // Use cpp-d4 data for consistency with CNCalculator (GFNFFCN)
    if (D4ReferenceData::refcn_cppd4.size() >= 7) {
        // cpp-d4 format: refcn[charge_state][element_number]
        // Transpose to our format: m_refcn[element_number][charge_state]
        m_refcn.resize(MAX_ELEM, std::vector<double>(MAX_REF, 0.0));
        for (int elem = 0; elem < MAX_ELEM && elem < 118; ++elem) {
            for (int charge_state = 0; charge_state < MAX_REF && charge_state < 7; ++charge_state) {
                // cpp-d4 has charge_state dimension first, then element
                if (charge_state < static_cast<int>(D4ReferenceData::refcn_cppd4.size()) &&
                    elem < static_cast<int>(D4ReferenceData::refcn_cppd4[charge_state].size())) {
                    m_refcn[elem][charge_state] = D4ReferenceData::refcn_cppd4[charge_state][elem];
                }
            }
        }
        if (CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::success("D4: Loaded reference CN data from cpp-d4 (86/118 elements)");
        }
    } else {
        if (CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::warn("D4: cpp-d4 CN data unavailable, using placeholders");
        }
        m_refcn.resize(MAX_ELEM, std::vector<double>(MAX_REF, 0.0));
    }

    // Initialize atomic scaling factors (WIP - Phase 2.2)
    // TODO: Extract ascale data from Fortran dftd4param.f90
    m_ascale.resize(MAX_ELEM, std::vector<double>(MAX_REF, 1.0));

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::warn("D4: Using placeholder refcn and ascale (Phase 2.2 pending)");
    }

    // Legacy data (retained for compatibility during transition)
    m_r4_over_r2.resize(MAX_ELEM, 0.0);
    m_sqrt_z_r4_r2.resize(MAX_ELEM, 0.0);

    // Basic atomic r4/r2 ratios from D4 reference (first few elements for now)
    if (MAX_ELEM >= 20) {
        m_r4_over_r2[0] = 5.1917362;  // H
        m_r4_over_r2[1] = 1.6325801;  // He
        m_r4_over_r2[5] = 32.6290215; // C
        m_r4_over_r2[6] = 24.9150620; // N
        m_r4_over_r2[7] = 20.0948582; // O
        m_r4_over_r2[10] = 90.3494164; // Na
        m_r4_over_r2[11] = 85.4616254; // Mg
        m_r4_over_r2[19] = 33.4697314; // Ca
    }

    for (int i = 0; i < MAX_ELEM; ++i) {
        m_sqrt_z_r4_r2[i] = std::sqrt((i + 1) * m_r4_over_r2[i]);
    }

    m_data_initialized = true;

    // Initialize frequency-dependent polarizabilities (simplified - real data extraction pending)
    m_alpha_iw.resize(MAX_ELEM, std::vector<std::vector<double>>(MAX_REF, std::vector<double>(N_FREQ, 0.0)));
    m_integrated_alpha.resize(MAX_ELEM, std::vector<double>(MAX_REF, 0.0));
}

void D4ParameterGenerator::calculateFrequencyDependentPolarizabilities()
{
    // Phase 2.1 (December 2025): Frequency-dependent polarizability calculation
    // Implements Casimir-Polder integration over 23-point frequency grid
    // Reference: Caldeweyher et al., J. Chem. Phys. 150, 154122 (2019)

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("=== D4: Calculating frequency-dependent polarizabilities ===");
    }

    // 23-point imaginary frequency grid (from GFN-FF Fortran dftd4param.f90)
    // Used for Casimir-Polder integration: α_static = (3/2π) ∫ α(iω) dω
    const std::vector<double> frequency_grid = {
        0.000001, 0.050000, 0.100000, 0.200000, 0.300000, 0.400000,
        0.500000, 0.600000, 0.700000, 0.800000, 0.900000, 1.000000,
        1.200000, 1.400000, 1.600000, 1.800000, 2.000000, 2.500000,
        3.000000, 4.000000, 5.000000, 7.500000, 10.000000
    };

    // Initialize 3D array: m_alpha_iw[element][reference][frequency]
    m_alpha_iw.resize(MAX_ELEM, std::vector<std::vector<double>>(MAX_REF, std::vector<double>(N_FREQ, 0.0)));
    m_integrated_alpha.resize(MAX_ELEM, std::vector<double>(MAX_REF, 0.0));

    // Simplified approach for Phase 2.1: Use approximate alpha_iw for key elements
    // Full extraction from Fortran dftd4param.f90 deferred to Phase 2.2

    // Hydrogen (Z=1) - 2 reference states (from Fortran secaiw data for H2)
    if (m_refn.size() > 0 && m_refn[0] >= 2) {
        std::vector<double> h_alpha = {
            5.4415160, 5.3912720, 5.2466780, 4.7462570, 4.1122050, 3.4827990,
            2.9256260, 2.4586020, 2.0763900, 1.7660350, 1.5138980, 1.3080740,
            0.9987770, 0.7833600, 0.6286810, 0.5145050, 0.4281480, 0.2867670,
            0.2047270, 0.1187560, 0.0772270, 0.0349350, 0.0197880
        };
        for (int ref = 0; ref < m_refn[0] && ref < MAX_REF; ++ref) {
            for (int iw = 0; iw < N_FREQ; ++iw) {
                m_alpha_iw[0][ref][iw] = h_alpha[iw];
            }
        }
    }

    // Carbon (Z=6) - 7 reference states (from Fortran secaiw data for C6H6)
    if (m_refn.size() > 5 && m_refn[5] >= 6) {
        std::vector<double> c_alpha = {
            68.5832590, 67.5115260, 64.6123080, 56.1286650, 47.4318310, 39.9459190,
            33.7814890, 28.7553020, 24.6561470, 21.2992860, 18.5340330, 16.2406480,
            12.7133690, 10.1832050, 8.3194640, 6.9133790, 5.8298100, 4.0106600,
            2.9230920, 1.7494800, 1.1654830, 0.5523060, 0.3242020
        };
        for (int ref = 0; ref < m_refn[5] && ref < MAX_REF; ++ref) {
            for (int iw = 0; iw < N_FREQ; ++iw) {
                m_alpha_iw[5][ref][iw] = c_alpha[iw];
            }
        }
    }

    // Integrate α(iω) over frequency grid using trapezoidal rule
    // Static polarizability: α_0 = (3/2π) ∫ α(iω) dω
    const double prefactor = THREE_OVER_PI;  // 3/(2π) with extra 1/2 absorbed

    for (int elem = 0; elem < MAX_ELEM; ++elem) {
        for (int ref = 0; ref < MAX_REF; ++ref) {
            double integral = 0.0;

            // Trapezoidal integration
            for (int iw = 0; iw < N_FREQ - 1; ++iw) {
                double dw = frequency_grid[iw + 1] - frequency_grid[iw];
                double avg_alpha = 0.5 * (m_alpha_iw[elem][ref][iw] + m_alpha_iw[elem][ref][iw + 1]);
                integral += avg_alpha * dw;
            }

            m_integrated_alpha[elem][ref] = prefactor * integral;
        }
    }

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::success("D4: Frequency-dependent polarizabilities calculated");

        // Print validation for H and C
        if (m_refn.size() > 0) {
            CurcumaLogger::result(fmt::format("  H (Z=1) α_0[0] = {:.4f}", m_integrated_alpha[0][0]));
        }
        if (m_refn.size() > 5) {
            CurcumaLogger::result(fmt::format("  C (Z=6) α_0[0] = {:.4f}", m_integrated_alpha[5][0]));
        }

        CurcumaLogger::warn("D4: Using simplified alpha_iw for H and C only (Phase 2.2: full extraction)");
    }
}

void D4ParameterGenerator::GenerateParameters(const std::vector<int>& atoms, const Matrix& geometry_bohr)
{
    m_atoms = atoms;

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("=== D4 Parameter Generation with EEQ Charges ===");
        CurcumaLogger::param("Number of atoms", static_cast<int>(m_atoms.size()));
    }

    // Claude Generated (2025): Calculate CN FIRST, then pass to EEQ to avoid duplicate calculation
    // This eliminates redundant O(n²) CN computation inside EEQ solver
    m_cn_values = CNCalculator::calculateGFNFFCN(m_atoms, geometry_bohr);

    if (m_cn_values.size() != m_atoms.size()) {
        CurcumaLogger::error("D4: CN calculation failed");
        return;
    }

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::success("D4: Molecular CN calculated (GFNFFCN)");
        if (CurcumaLogger::get_verbosity() >= 3) {
            for (size_t i = 0; i < std::min(size_t(5), m_atoms.size()); ++i) {
                CurcumaLogger::result(fmt::format("  Atom {} (Z={}) CN = {:.4f}",
                                                  i, m_atoms[i], m_cn_values[i]));
            }
        }
    }

    // STEP 1: Calculate EEQ charges from geometry (Dec 2025 - Phase 2)
    // Pass pre-calculated CN to avoid duplicate calculation
    // Convert std::vector<double> to Eigen::Vector for EEQ solver
    Vector cn_eigen = Eigen::Map<const Vector>(m_cn_values.data(), m_cn_values.size());
    m_eeq_charges = m_eeq_solver->calculateCharges(m_atoms, geometry_bohr, 0, &cn_eigen);

    if (m_eeq_charges.size() != m_atoms.size()) {
        CurcumaLogger::error("D4: EEQ charge calculation failed");
        return;
    }

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::success("D4: EEQ charges calculated");
        if (CurcumaLogger::get_verbosity() >= 3) {
            for (size_t i = 0; i < std::min(size_t(5), m_atoms.size()); ++i) {
                CurcumaLogger::result(fmt::format("  Atom {} (Z={}) q = {:.6f}",
                                                  i, m_atoms[i], m_eeq_charges(i)));
            }
        }
    }

    // STEP 2: Generate C6 pairs with CN+charge weighted C6 coefficients
    json dispersion_pairs = json::array();

    for (size_t i = 0; i < m_atoms.size(); ++i) {
        for (size_t j = i + 1; j < m_atoms.size(); ++j) {
            int atom_i = m_atoms[i];
            int atom_j = m_atoms[j];

            if (atom_i > 0 && atom_i <= MAX_ELEM && atom_j > 0 && atom_j <= MAX_ELEM) {
                json pair;
                pair["i"] = static_cast<int>(i);
                pair["j"] = static_cast<int>(j);
                pair["element_i"] = atom_i;
                pair["element_j"] = atom_j;

                // NEW: CN+charge weighted C6 using EEQ charges and GFNFFCN (Dec 2025 Phase 2.2)
                double c6 = getChargeWeightedC6(atom_i, atom_j, static_cast<int>(i), static_cast<int>(j));

                // D4 C6 coefficient is already charge-weighted via calculateChargeWeightedC6()
                // with Casimir-Polder integration over frequency-dependent polarizabilities

                // FIX (Dec 25, 2025): C8 IS used in pairwise D4!
                // BJ damping formula: E = -s6·C6/(r⁶+R0⁶) - s8·C8/(r⁸+R0⁸)
                // C8 = 3 * C6 * sqrt(Q_i * Q_j) where Q = <r⁴>/<r²> (atomic moment ratio)
                double sqrt_q_i = getSqrtZr4r2(atom_i);
                double sqrt_q_j = getSqrtZr4r2(atom_j);
                double c8 = 3.0 * c6 * sqrt_q_i * sqrt_q_j;

                // C10: Only needed for three-body corrections (not implemented yet)
                double c10 = 0.0;
                double c12 = 0.0;

                // Get D4 damping parameters from config (GFN-FF defaults)
                // Reference: Spicher/Grimme, Angew. Chem. Int. Ed. 2020, DOI: 10.1002/anie.202004239
                double s6 = m_config.get<double>("d4_s6", 1.0);
                double s8 = m_config.get<double>("d4_s8", 1.0);
                double a1 = m_config.get<double>("d4_a1", 0.44);  // GFN-FF D4 value
                double a2 = m_config.get<double>("d4_a2", 4.60);  // GFN-FF D4 value (Bohr)

                // CRITICAL FIX (Dec 2025): Match D3 JSON format for ForceField compatibility
                // Use uppercase C6/C8 and include s6/s8/a1/a2/r_cut fields
                // Claude Generated - Dec 25, 2025: Add dispersion_method tag to route to CalculateD4DispersionContribution()
                pair["dispersion_method"] = "d4";  // Route to native D4 in ForceFieldThread
                pair["C6"] = c6;  // Raw C6 (s6 applied in energy calculation)
                pair["C8"] = c8;  // Raw C8 (s8 applied in energy calculation)
                pair["s6"] = s6;
                pair["s8"] = s8;
                pair["a1"] = a1;
                pair["a2"] = a2;
                pair["r_cut"] = 100.0;  // Cutoff radius (Bohr)

                dispersion_pairs.push_back(pair);
            } else {
                if (CurcumaLogger::get_verbosity() >= 2) {
                    CurcumaLogger::warn("D4: Elements out of range - i=" +
                                       std::to_string(atom_i) +
                                       " j=" + std::to_string(atom_j));
                }
            }
        }
    }

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::param("Generated D4 pairs", static_cast<int>(dispersion_pairs.size()));
    }

    // Summary of C6 distribution (verbosity 2) - Fix 2
    if (CurcumaLogger::get_verbosity() >= 2) {
        double c6_min = std::numeric_limits<double>::max();
        double c6_max = 0.0;
        double c6_avg = 0.0;
        int zero_count = 0;

        for (const auto& pair : dispersion_pairs) {
            double c6 = pair["C6"];
            if (c6 < 1e-10) {
                zero_count++;
                if (CurcumaLogger::get_verbosity() >= 3) {
                    CurcumaLogger::warn(fmt::format("Zero C6 for pair [{},{}] (Zi={} Zj={})",
                        static_cast<int>(pair["i"]), static_cast<int>(pair["j"]),
                        static_cast<int>(pair["element_i"]), static_cast<int>(pair["element_j"])));
                }
            }
            c6_min = std::min(c6_min, c6);
            c6_max = std::max(c6_max, c6);
            c6_avg += c6;
        }

        if (!dispersion_pairs.empty()) {
            c6_avg /= dispersion_pairs.size();
            CurcumaLogger::success(fmt::format("D4 C6: {:.0f} pairs, avg={:.4f}, min={:.4f}, max={:.4f}, zeros={}",
                dispersion_pairs.size(), c6_avg, c6_min, c6_max, zero_count));
        }
    }

    m_parameters["d4_dispersion_pairs"] = dispersion_pairs;
    m_parameters["d4_damping"] = {
        {"a1", m_config.get<double>("d4_a1", 0.44)},  // GFN-FF D4 (Spicher/Grimme 2020)
        {"a2", m_config.get<double>("d4_a2", 4.60)},  // GFN-FF D4 (Bohr)
        {"alp", m_config.get<double>("d4_alp", 14.0)}
    };
    m_parameters["d4_enabled"] = true;
    m_parameters["d4_refq"] = m_config.get<int>("d4_refq", 2); // Hirshfeld charges default
    m_parameters["d4_r4r2_model"] = m_config.get<int>("d4_r4r2_model", 1);

    // Phase 2.4 (December 2025): D4 ATM three-body dispersion
    // Same formula as D3, but uses charge-weighted C6 coefficients
    json atm_triples = json::array();

    if (m_config.get<double>("d4_s9", 0.0) > 1e-10) {
        double s9 = m_config.get<double>("d4_s9", 1.0);
        double a1 = m_config.get<double>("d4_a1", 0.44);
        double a2 = m_config.get<double>("d4_a2", 4.60);
        double alp = m_config.get<double>("d4_alp", 14.0);

        int n_atoms = static_cast<int>(m_atoms.size());

        for (int i = 0; i < n_atoms; ++i) {
            for (int j = 0; j < i; ++j) {
                for (int k = 0; k < j; ++k) {
                    json triple;
                    triple["i"] = i;
                    triple["j"] = j;
                    triple["k"] = k;

                    // C6 from charge-weighted D4 coefficients
                    int Zi = m_atoms[i];
                    int Zj = m_atoms[j];
                    int Zk = m_atoms[k];

                    double c6_ij = getChargeWeightedC6(Zi, Zj, i, j);
                    double c6_ik = getChargeWeightedC6(Zi, Zk, i, k);
                    double c6_jk = getChargeWeightedC6(Zj, Zk, j, k);

                    triple["C6_ij"] = c6_ij;
                    triple["C6_ik"] = c6_ik;
                    triple["C6_jk"] = c6_jk;
                    triple["s9"] = s9;
                    triple["a1"] = a1;
                    triple["a2"] = a2;
                    triple["alp"] = alp;
                    triple["atm_method"] = "d4";

                    // Symmetry factor (same as D3)
                    triple["triple_scale"] = calculateTripleScale(i, j, k);

                    atm_triples.push_back(triple);
                }
            }
        }

        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::param("Generated D4 ATM triples", static_cast<int>(atm_triples.size()));
        }
    }

    m_parameters["atm_triples"] = atm_triples;

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::success("D4 parameter generation completed");
    }
}

double D4ParameterGenerator::getEffectiveC6(int atom_i, int atom_j) const
{
    // Convert to 0-based for internal arrays
    int elem_i = atom_i - 1;
    int elem_j = atom_j - 1;

    if (elem_i < 0 || elem_i >= MAX_ELEM || elem_j < 0 || elem_j >= MAX_ELEM || !m_data_initialized) {
        return 0.0;
    }

    // Use integrated polarizabilities for more accurate C6 estimation
    // This mimics the D4 approach with charge-state dependence
    double alpha_i = (elem_i < static_cast<int>(m_integrated_alpha.size())) ?
                     m_integrated_alpha[elem_i][0] : 1.0;
    double alpha_j = (elem_j < static_cast<int>(m_integrated_alpha.size())) ?
                     m_integrated_alpha[elem_j][0] : 1.0;

    // Weighted geometric mean accounting for atomic number differences
    double z_i = atom_i;
    double z_j = atom_j;
    double weight_factor = 2.0 * std::sqrt(z_i * z_j) / (z_i + z_j);

    // Simplified London dispersion formula
    double c6 = weight_factor * 1.5 * alpha_i * alpha_j;

    return c6;
}

double D4ParameterGenerator::getR4OverR2(int atom) const
{
    --atom; // Convert to 0-based
    if (atom >= 0 && atom < MAX_ELEM && atom < static_cast<int>(m_r4_over_r2.size())) {
        return m_r4_over_r2[atom];
    }
    return 10.0; // Default fallback
}

double D4ParameterGenerator::getSqrtZr4r2(int atom) const
{
    --atom; // Convert to 0-based
    if (atom >= 0 && atom < MAX_ELEM && atom < static_cast<int>(m_sqrt_z_r4_r2.size())) {
        return m_sqrt_z_r4_r2[atom];
    }
    return 1.0; // Default fallback
}

double D4ParameterGenerator::getAtomicPolarizability(int atom, int frequency_index) const
{
    --atom; // Convert to 0-based
    if (atom >= 0 && atom < MAX_ELEM && atom < static_cast<int>(m_alpha_iw.size()) &&
        frequency_index >= 0 && frequency_index < N_FREQ &&
        frequency_index < static_cast<int>(m_alpha_iw[atom][0].size())) {
        return m_alpha_iw[atom][0][frequency_index];
    }
    return 1.0; // Default value
}

/**
 * @brief Calculate charge-weighted C6 coefficient using Gaussian charge-state weighting
 *
 * Implements D4 core algorithm: Weight reference C6 values by Gaussian functions
 * based on atomic charges from EEQ calculation.
 *
 * Reference: E. Caldeweyher et al., J. Chem. Phys. 2019, 150, 154122 (D4 method)
 *
 * Formula: w_k = exp(-α * (q - q_ref_k)²) / Σ exp(...)
 *          C6 = Σ_k w_k * C6_ref_k
 *
 * Claude Generated - December 2025 (Phase 2: D4-EEQ Integration)
 *
 * @param Zi Atomic number of atom i
 * @param Zj Atomic number of atom j
 * @param qi EEQ charge of atom i
 * @param qj EEQ charge of atom j
 * @return Charge-weighted C6 coefficient (Hartree * Bohr^6)
 */
double D4ParameterGenerator::getChargeWeightedC6(int Zi, int Zj, int atom_i, int atom_j) const
{
    // Phase 2.2 (December 2025): Gaussian CN+charge-state weighting for C6 coefficients
    // Reference: E. Caldeweyher et al., J. Chem. Phys. 2019, 150, 154122 (D4 method)
    // cpp-d4 reference: dftd_model.h weight_cn() function
    //
    // Formula: C6(qi,qj,cni,cnj) = Σ_refi Σ_refj w(qi,cni,refi) * w(qj,cnj,refj) * C6_ref(refi,refj)
    // where: w(q,cn,ref) = exp(-wf * ((q - q_ref)² + (cn - cn_ref)²)) / norm
    //
    // Key change from charge-only (Phase 2.1): Added CN-dependent Gaussian term

    // Validate indices
    if (atom_i < 0 || atom_i >= static_cast<int>(m_eeq_charges.size()) ||
        atom_j < 0 || atom_j >= static_cast<int>(m_eeq_charges.size()) ||
        atom_i < 0 || atom_i >= static_cast<int>(m_cn_values.size()) ||
        atom_j < 0 || atom_j >= static_cast<int>(m_cn_values.size())) {
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::warn(fmt::format("D4: Invalid atom indices i={} j={}", atom_i, atom_j));
        }
        return 0.0;
    }

    double qi = m_eeq_charges(atom_i);
    double qj = m_eeq_charges(atom_j);
    double cni = m_cn_values[atom_i];
    double cnj = m_cn_values[atom_j];

    constexpr double wf = 4.0;  // Gaussian width parameter (from D4 cpp-d4)

    // Convert to 0-based indexing
    int elem_i = Zi - 1;
    int elem_j = Zj - 1;

    // Validate element range
    if (elem_i < 0 || elem_i >= MAX_ELEM || elem_j < 0 || elem_j >= MAX_ELEM) {
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::warn(fmt::format("D4: Invalid elements Zi={} Zj={}", Zi, Zj));
        }
        return 0.0;
    }

    // Get number of reference states for each element
    int nref_i = (elem_i < static_cast<int>(m_refn.size())) ? m_refn[elem_i] : 1;
    int nref_j = (elem_j < static_cast<int>(m_refn.size())) ? m_refn[elem_j] : 1;

    // Gaussian CN+charge-state weighting
    double c6_weighted = 0.0;
    double norm_i = 0.0;
    double norm_j = 0.0;

    // Calculate normalization factors with CN+charge combined weighting
    for (int refi = 0; refi < nref_i && refi < MAX_REF; ++refi) {
        double qi_ref = m_refq[elem_i][refi];
        double cni_ref = (elem_i < static_cast<int>(m_refcn.size()) && refi < static_cast<int>(m_refcn[elem_i].size()))
                        ? m_refcn[elem_i][refi] : 0.0;

        // KEY CHANGE: Combined CN+charge Gaussian weighting (December 2025 Phase 2.2)
        double wi = std::exp(-wf * ((qi - qi_ref) * (qi - qi_ref) + (cni - cni_ref) * (cni - cni_ref)));
        norm_i += wi;
    }

    for (int refj = 0; refj < nref_j && refj < MAX_REF; ++refj) {
        double qj_ref = m_refq[elem_j][refj];
        double cnj_ref = (elem_j < static_cast<int>(m_refcn.size()) && refj < static_cast<int>(m_refcn[elem_j].size()))
                        ? m_refcn[elem_j][refj] : 0.0;

        // KEY CHANGE: Combined CN+charge Gaussian weighting (December 2025 Phase 2.2)
        double wj = std::exp(-wf * ((qj - qj_ref) * (qj - qj_ref) + (cnj - cnj_ref) * (cnj - cnj_ref)));
        norm_j += wj;
    }

    // Avoid division by zero
    if (norm_i < 1e-10 || norm_j < 1e-10) {
        norm_i = std::max(norm_i, 1e-10);
        norm_j = std::max(norm_j, 1e-10);
    }

    // Frequency grid for Casimir-Polder integration (same as calculateFrequencyDependentPolarizabilities)
    const std::vector<double> frequency_grid = {
        0.000001, 0.050000, 0.100000, 0.200000, 0.300000, 0.400000,
        0.500000, 0.600000, 0.700000, 0.800000, 0.900000, 1.000000,
        1.200000, 1.400000, 1.600000, 1.800000, 2.000000, 2.500000,
        3.000000, 4.000000, 5.000000, 7.500000, 10.000000
    };

    // Weighted sum over reference states
    for (int refi = 0; refi < nref_i && refi < MAX_REF; ++refi) {
        double qi_ref = m_refq[elem_i][refi];
        double cni_ref = (elem_i < static_cast<int>(m_refcn.size()) && refi < static_cast<int>(m_refcn[elem_i].size()))
                        ? m_refcn[elem_i][refi] : 0.0;

        // KEY CHANGE: Combined CN+charge Gaussian weighting (December 2025 Phase 2.2)
        double wi = std::exp(-wf * ((qi - qi_ref) * (qi - qi_ref) + (cni - cni_ref) * (cni - cni_ref))) / norm_i;

        for (int refj = 0; refj < nref_j && refj < MAX_REF; ++refj) {
            double qj_ref = m_refq[elem_j][refj];
            double cnj_ref = (elem_j < static_cast<int>(m_refcn.size()) && refj < static_cast<int>(m_refcn[elem_j].size()))
                            ? m_refcn[elem_j][refj] : 0.0;

            // KEY CHANGE: Combined CN+charge Gaussian weighting (December 2025 Phase 2.2)
            double wj = std::exp(-wf * ((qj - qj_ref) * (qj - qj_ref) + (cnj - cnj_ref) * (cnj - cnj_ref))) / norm_j;

            // CRITICAL FIX (Dec 25, 2025): Correct C6 calculation from Casimir-Polder integration
            // Reference: XTB dftd4.F90 lines ~500: c6 = thopi * trapzd(alpha_i * alpha_j)
            //
            // Formula: C6_ij = (3/π) ∫ α_i(iω) * α_j(iω) dω
            // NOT: C6 = α_i * α_j (this was the 100x error bug!)

            double c6_ref = 0.0;

            // Debug: Track negative corrections for this reference state pair (Fix 3)
            bool negative_detected = false;

            // Check if we have frequency-dependent polarizabilities for both elements
            // Use d4_alphaiw_data (complete data from Fortran extraction - Dec 25, 2025)
            if (elem_i < static_cast<int>(d4_alphaiw_data.size()) &&
                elem_j < static_cast<int>(d4_alphaiw_data.size()) &&
                refi < static_cast<int>(d4_alphaiw_data[elem_i].size()) &&
                refj < static_cast<int>(d4_alphaiw_data[elem_j].size())) {

                // Get correction factors for reference states (Dec 25, 2025 - Phase 2.4)
                double ascale_i = (elem_i < static_cast<int>(d4_ascale_data.size()) && refi < static_cast<int>(d4_ascale_data[elem_i].size()))
                                ? d4_ascale_data[elem_i][refi] : 1.0;
                double ascale_j = (elem_j < static_cast<int>(d4_ascale_data.size()) && refj < static_cast<int>(d4_ascale_data[elem_j].size()))
                                ? d4_ascale_data[elem_j][refj] : 1.0;

                double hcount_i = (elem_i < static_cast<int>(m_refh.size()) && refi < static_cast<int>(m_refh[elem_i].size()))
                                ? m_refh[elem_i][refi] : 0.0;
                double hcount_j = (elem_j < static_cast<int>(m_refh.size()) && refj < static_cast<int>(m_refh[elem_j].size()))
                                ? m_refh[elem_j][refj] : 0.0;

                int refsys_i = (elem_i < static_cast<int>(d4_refsys_data.size()) && refi < static_cast<int>(d4_refsys_data[elem_i].size()))
                             ? d4_refsys_data[elem_i][refi] : 0;
                int refsys_j = (elem_j < static_cast<int>(d4_refsys_data.size()) && refj < static_cast<int>(d4_refsys_data[elem_j].size()))
                             ? d4_refsys_data[elem_j][refj] : 0;

                double sscale_i = (d4_sscale_data.find(refsys_i) != d4_sscale_data.end()) ? d4_sscale_data.at(refsys_i) : 0.0;
                double sscale_j = (d4_sscale_data.find(refsys_j) != d4_sscale_data.end()) ? d4_sscale_data.at(refsys_j) : 0.0;

                // Integrate product of CORRECTED polarizabilities using trapezoidal rule
                // Correction formula: α_corrected = ascale * (αᵢⱼw - hcount * sscale * secaiw)
                for (int iw = 0; iw < N_FREQ - 1; ++iw) {
                    // Get raw alphaiw values
                    double alphaiw_i_iw = d4_alphaiw_data[elem_i][refi][iw];
                    double alphaiw_i_next = d4_alphaiw_data[elem_i][refi][iw + 1];
                    double alphaiw_j_iw = d4_alphaiw_data[elem_j][refj][iw];
                    double alphaiw_j_next = d4_alphaiw_data[elem_j][refj][iw + 1];

                    // Get secaiw reference polarizabilities (if available)
                    double secaiw_i_iw = 0.0, secaiw_i_next = 0.0;
                    double secaiw_j_iw = 0.0, secaiw_j_next = 0.0;

                    if (d4_secaiw_data.find(refsys_i) != d4_secaiw_data.end() && iw < static_cast<int>(d4_secaiw_data.at(refsys_i).size())) {
                        secaiw_i_iw = d4_secaiw_data.at(refsys_i)[iw];
                        secaiw_i_next = d4_secaiw_data.at(refsys_i)[iw + 1];
                    }
                    if (d4_secaiw_data.find(refsys_j) != d4_secaiw_data.end() && iw < static_cast<int>(d4_secaiw_data.at(refsys_j).size())) {
                        secaiw_j_iw = d4_secaiw_data.at(refsys_j)[iw];
                        secaiw_j_next = d4_secaiw_data.at(refsys_j)[iw + 1];
                    }

                    // Apply correction formula
                    double alpha_i_iw = ascale_i * (alphaiw_i_iw - hcount_i * sscale_i * secaiw_i_iw);
                    double alpha_i_next = ascale_i * (alphaiw_i_next - hcount_i * sscale_i * secaiw_i_next);
                    double alpha_j_iw = ascale_j * (alphaiw_j_iw - hcount_j * sscale_j * secaiw_j_iw);
                    double alpha_j_next = ascale_j * (alphaiw_j_next - hcount_j * sscale_j * secaiw_j_next);

                    // DEBUG: Detect negative corrections before clamping (Fix 3)
                    if (alpha_i_iw < 0 || alpha_i_next < 0 || alpha_j_iw < 0 || alpha_j_next < 0) {
                        negative_detected = true;
                    }

                    // Ensure non-negative (as per Fortran: max(correction, 0.0))
                    alpha_i_iw = std::max(alpha_i_iw, 0.0);
                    alpha_i_next = std::max(alpha_i_next, 0.0);
                    alpha_j_iw = std::max(alpha_j_iw, 0.0);
                    alpha_j_next = std::max(alpha_j_next, 0.0);

                    // Product at current and next frequency points
                    double product_iw = alpha_i_iw * alpha_j_iw;
                    double product_next = alpha_i_next * alpha_j_next;

                    // Trapezoidal rule: ∫ f(x) dx ≈ Σ (f_i + f_{i+1})/2 * Δx
                    double dw = frequency_grid[iw + 1] - frequency_grid[iw];
                    c6_ref += 0.5 * (product_iw + product_next) * dw;
                }

                // Apply prefactor: 3/π (matches XTB's thopi constant)
                c6_ref *= THREE_OVER_PI;

            } else {
                // Fallback for elements without polarizability data
                // Use simple geometric estimate (will be inaccurate)
                if (CurcumaLogger::get_verbosity() >= 3) {
                    CurcumaLogger::warn(fmt::format("D4: No alpha_iw data for Zi={} Zj={}, using fallback", Zi, Zj));
                }
                c6_ref = 1.0;  // Minimal fallback value
            }

            // Add weighted contribution
            c6_weighted += wi * wj * c6_ref;

            // Warning for negative correction clamping (Fix 3)
            if (CurcumaLogger::get_verbosity() >= 3 && negative_detected) {
                CurcumaLogger::warn(fmt::format("Negative correction clamped: Zi={} Zj={}, refi={}, refj={}, c6={:.4f}",
                                                 Zi, Zj, refi, refj, c6_ref));
            }
        }
    }

    // Fix 1: Simplified debug output at level 3 (using info which has threshold >= 2)
    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info(fmt::format("D4 C6[{}][{}]: Zi={} Zj={} → C6={:.4f}",
                                         atom_i, atom_j, Zi, Zj, c6_weighted));
    }

    return c6_weighted;
}

// Claude Generated (2025): ATM three-body symmetry factor calculation
double D4ParameterGenerator::calculateTripleScale(int i, int j, int k) const
{
    // Reference: external/cpp-d4/src/damping/atm.cpp:291-313
    if (i == j) {
        return (i == k) ? 1.0/6.0 : 0.5;  // iii: 1/6, iij: 1/2
    } else {
        return (i != k && j != k) ? 1.0 : 0.5;  // ijk: 1, ijj/iji: 1/2
    }
}
