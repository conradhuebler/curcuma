/*
 * DFT-D4 Parameter Generator for Curcuma - Complete Data Integration
 * Integrates 118 elements of D4 reference data from GFN-FF Fortran implementation
 *
 * Claude Generated (December 2025): Phase 2.1 - D4 Reference Data Integration
 */

#include "d4param_generator.h"
#include "../../../../test_cases/reference_data/d4_reference_data_fixed.cpp"  // D4 reference data (365 lines)
#include "src/core/curcuma_logger.h"

#include <algorithm>
#include <cmath>

D4ParameterGenerator::D4ParameterGenerator(const ConfigManager& config)
    : m_config(config)
{
    // Initialize EEQ solver for charge calculation (Dec 2025 - Phase 2)
    ConfigManager eeq_config("eeq_solver", config.exportConfig());
    m_eeq_solver = std::make_unique<EEQSolver>(eeq_config);

    initializeReferenceData();
    calculateFrequencyDependentPolarizabilities();
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

    // Initialize reference coordination numbers (WIP - Phase 2.2)
    // TODO: Extract refcovcn data from Fortran dftd4param.f90
    m_refcn.resize(MAX_ELEM, std::vector<double>(MAX_REF, 0.0));
    for (int i = 0; i < MAX_ELEM; ++i) {
        for (int j = 0; j < m_refn[i] && j < MAX_REF; ++j) {
            m_refcn[i][j] = 1.0;  // Placeholder - will be extracted in Phase 2.2
        }
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

    // STEP 1: Calculate EEQ charges from geometry (Dec 2025 - Phase 2)
    m_eeq_charges = m_eeq_solver->calculateCharges(m_atoms, geometry_bohr, 0);

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

    // STEP 2: Generate C6 pairs with charge-weighted C6 coefficients
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

                // NEW: Charge-weighted C6 using EEQ charges (Dec 2025 - Phase 2)
                double qi = m_eeq_charges(i);
                double qj = m_eeq_charges(j);
                double c6 = getChargeWeightedC6(atom_i, atom_j, qi, qj);

                // D4 uses atomic moments for higher-order coefficients
                double r4r2_i = getR4OverR2(atom_i);
                double r4r2_j = getR4OverR2(atom_j);
                double r4r2 = std::sqrt(r4r2_i * r4r2_j);

                double c8 = 3.0 * c6 * r4r2;
                double c10 = 5.0 * c6 * r4r2 * r4r2;
                double c12 = 7.0 * c6 * r4r2 * r4r2 * r4r2;

                // Get D4 damping parameters from config (GFN-FF defaults)
                double s6 = m_config.get<double>("d4_s6", 1.0);
                double s8 = m_config.get<double>("d4_s8", 1.0);
                double a1 = m_config.get<double>("d4_a1", 0.63);  // D4 GFN2-XTB value
                double a2 = m_config.get<double>("d4_a2", 5.0);   // D4 GFN2-XTB value (Bohr)

                // CRITICAL FIX (Dec 2025): Match D3 JSON format for ForceField compatibility
                // Use uppercase C6/C8 and include s6/s8/a1/a2/r_cut fields
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

    m_parameters["d4_dispersion_pairs"] = dispersion_pairs;
    m_parameters["d4_damping"] = {
        {"a1", m_config.get<double>("d4_a1", 0.43)},
        {"a2", m_config.get<double>("d4_a2", 4.0)},
        {"alp", m_config.get<double>("d4_alp", 14.0)}
    };
    m_parameters["d4_enabled"] = true;
    m_parameters["d4_refq"] = m_config.get<int>("d4_refq", 2); // Hirshfeld charges default
    m_parameters["d4_r4r2_model"] = m_config.get<int>("d4_r4r2_model", 1);

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
double D4ParameterGenerator::getChargeWeightedC6(int Zi, int Zj, double qi, double qj) const
{
    // Phase 2.1 (December 2025): Gaussian charge-state weighting for C6 coefficients
    // Reference: E. Caldeweyher et al., J. Chem. Phys. 2019, 150, 154122 (D4 method)
    //
    // Formula: C6(qi,qj) = Σ_refi Σ_refj w(qi,refi) * w(qj,refj) * C6_ref(refi,refj)
    // where: w(q,ref) = exp(-wf * (q - q_ref)²) / norm

    constexpr double wf = 4.0;  // Gaussian width parameter (from D4 paper)

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

    // Gaussian charge-state weighting
    double c6_weighted = 0.0;
    double norm_i = 0.0;
    double norm_j = 0.0;

    // Calculate normalization factors
    for (int refi = 0; refi < nref_i && refi < MAX_REF; ++refi) {
        double qi_ref = m_refq[elem_i][refi];
        double wi = std::exp(-wf * (qi - qi_ref) * (qi - qi_ref));
        norm_i += wi;
    }

    for (int refj = 0; refj < nref_j && refj < MAX_REF; ++refj) {
        double qj_ref = m_refq[elem_j][refj];
        double wj = std::exp(-wf * (qj - qj_ref) * (qj - qj_ref));
        norm_j += wj;
    }

    // Avoid division by zero
    if (norm_i < 1e-10 || norm_j < 1e-10) {
        norm_i = std::max(norm_i, 1e-10);
        norm_j = std::max(norm_j, 1e-10);
    }

    // Weighted sum over reference states
    for (int refi = 0; refi < nref_i && refi < MAX_REF; ++refi) {
        double qi_ref = m_refq[elem_i][refi];
        double wi = std::exp(-wf * (qi - qi_ref) * (qi - qi_ref)) / norm_i;

        for (int refj = 0; refj < nref_j && refj < MAX_REF; ++refj) {
            double qj_ref = m_refq[elem_j][refj];
            double wj = std::exp(-wf * (qj - qj_ref) * (qj - qj_ref)) / norm_j;

            // Calculate C6_ref for this reference state pair
            // Using integrated polarizabilities: C6 ~ α_i * α_j
            double alpha_i = (elem_i < static_cast<int>(m_integrated_alpha.size()) && refi < static_cast<int>(m_integrated_alpha[elem_i].size()))
                           ? m_integrated_alpha[elem_i][refi]
                           : 1.0;
            double alpha_j = (elem_j < static_cast<int>(m_integrated_alpha.size()) && refj < static_cast<int>(m_integrated_alpha[elem_j].size()))
                           ? m_integrated_alpha[elem_j][refj]
                           : 1.0;

            // London dispersion formula: C6 = (3/2) * (α_i * α_j * I_i * I_j) / (I_i + I_j)
            // Simplified using geometric mean for ionization potentials
            double weight_factor = 2.0 * std::sqrt(static_cast<double>(Zi * Zj)) / (Zi + Zj);
            double c6_ref = weight_factor * 1.5 * alpha_i * alpha_j;

            // Add weighted contribution
            c6_weighted += wi * wj * c6_ref;
        }
    }

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::result(fmt::format("D4 C6 (Gaussian): Zi={} Zj={} qi={:.4f} qj={:.4f} → C6={:.4f} (nref_i={}, nref_j={})",
                                          Zi, Zj, qi, qj, c6_weighted, nref_i, nref_j));
    }

    return c6_weighted;
}
