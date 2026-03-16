/*
 * DFT-D4 Parameter Generator for Curcuma - Complete Data Integration
 * Integrates 118 elements of D4 reference data from GFN-FF Fortran implementation
 *
 * Claude Generated (December 2025): Phase 2.1 - D4 Reference Data Integration
 * December 25, 2025: Complete alphaiw (polarizability) data integration
 */

#include "d4param_generator.h"
#include "gfnff_par.h"
#include "d4_reference_data_fixed.cpp"  // D4 reference data (365 lines)
#include "d4_reference_cn_fortran.cpp"  // D4 reference CN data (Fortran dftd3param.f90 - January 2026)
#include "d4_alphaiw_data.cpp"  // D4 alphaiw polarizability data
#include "d4_corrections_data.cpp"  // D4 correction factors
#include "../../curcuma_logger.h"

#include <algorithm>
#include <cmath>
#include <chrono>
#include <fstream>  // Claude Generated (Mar 5, 2026): For diagnostic JSON output
#include <limits>
#include <thread>
#include <omp.h>  // Claude Generated (February 2026): Phase 3 - OpenMP parallelization

// External declarations from d4_reference_cn_fortran.cpp (Fortran dftd3param.f90 - January 2026)
// CRITICAL FIX: Use Fortran reference CN data (matches dftd3param.f90) instead of cpp-d4 data
// This fixes ~20 mEh dispersion error caused by wrong CN weighting
extern const std::vector<std::vector<double>> D4ReferenceData::refcn_fortran;

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

    // Claude Generated (Dec 27, 2025): Pre-compute C6 reference matrix
    // This is done ONCE at initialization, independent of molecular geometry
    // Expected to eliminate ~98% of D4 parameter generation time for large molecules
    precomputeC6ReferenceMatrix();
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

    // Load reference coordination numbers (January 2026 - Fortran dftd3param.f90 data)
    // CRITICAL FIX: Use Fortran reference CN data (matches dftd3param.f90) instead of cpp-d4 data
    // This fixes ~20 mEh dispersion error caused by wrong CN weighting
    if (D4ReferenceData::refcn_fortran.size() >= 7) {
        // Fortran format: refcn[charge_state][element_number]
        // Transpose to our format: m_refcn[element_number][charge_state]
        m_refcn.resize(MAX_ELEM, std::vector<double>(MAX_REF, 0.0));
        for (int elem = 0; elem < MAX_ELEM && elem < 118; ++elem) {
            for (int charge_state = 0; charge_state < MAX_REF && charge_state < 7; ++charge_state) {
                // Fortran data has charge_state dimension first, then element
                if (charge_state < static_cast<int>(D4ReferenceData::refcn_fortran.size()) &&
                    elem < static_cast<int>(D4ReferenceData::refcn_fortran[charge_state].size())) {
                    m_refcn[elem][charge_state] = D4ReferenceData::refcn_fortran[charge_state][elem];
                }
            }
        }
        if (CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::success("D4: Loaded reference CN data from Fortran dftd3param.f90");
        }
    } else {
        if (CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::warn("D4: Fortran CN data unavailable, using placeholders");
        }
        m_refcn.resize(MAX_ELEM, std::vector<double>(MAX_REF, 0.0));
    }

    // Legacy data (retained for compatibility during transition)
    m_r4_over_r2.resize(MAX_ELEM, 0.0);
    m_sqrt_z_r4_r2.resize(MAX_ELEM, 0.0);

    // GFN-FF Moment Ratios (r4Overr2) from dftd4param.f90:134-157
    // These are essential for correct damping radii R0 in GFN-FF
    m_r4_over_r2 = {
        8.0589, 3.4698,  // H-He
        29.0974, 14.8517, 11.8799, 7.8715, 5.5588, 4.7566, 3.8025, 3.1036,  // Li-Ne
        26.1552, 17.2304, 17.7210, 12.7442, 9.5361, 8.1652, 6.7463, 5.6004,  // Na-Ar
        29.2012, 22.3934, // K-Ca
        19.0598, 16.8590, 15.4023, 12.5589, 13.4788, // Sc-Mn
        12.2309, 11.2809, 10.5569, 10.1428, 9.4907,  // Fe-Zn
        13.4606, 10.8544, 8.9386, 8.1350, 7.1251, 6.1971 // Ga-Kr
    };
    if (m_r4_over_r2.size() < MAX_ELEM) m_r4_over_r2.resize(MAX_ELEM, 10.0);

    for (int i = 0; i < MAX_ELEM; ++i) {
        // Reference formula for sqrtZr4r2 (gfnff_param.f90:376 and dftd4param.f90:160):
        // sqrtZr4r2 = sqrt(0.5 * r4Overr2 * sqrt(Z))
        double Z = static_cast<double>(i + 1);
        m_sqrt_z_r4_r2[i] = std::sqrt(0.5 * m_r4_over_r2[i] * std::sqrt(Z));
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
    auto t_start_total = std::chrono::high_resolution_clock::now();

    m_atoms = atoms;

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("=== D4 Parameter Generation with EEQ Charges ===");
        CurcumaLogger::param("Number of atoms", static_cast<int>(m_atoms.size()));
    }

    // Claude Generated (2025): Calculate CN FIRST, then pass to EEQ to avoid duplicate calculation
    // This eliminates redundant O(n²) CN computation inside EEQ solver
    auto t_cn_start = std::chrono::high_resolution_clock::now();
    m_cn_values = CNCalculator::calculateGFNFFCN(m_atoms, geometry_bohr);
    auto t_cn_end = std::chrono::high_resolution_clock::now();
    double t_cn_ms = std::chrono::duration<double, std::milli>(t_cn_end - t_cn_start).count();

    if (m_cn_values.size() != m_atoms.size()) {
        CurcumaLogger::error("D4: CN calculation failed");
        return;
    }

    // Claude Generated (Feb 8, 2026): Diagnostic for dispersion CN values (verbosity >= 3 only)
    // Claude Generated (Mar 5, 2026): Removed atom count limit for large-system diagnostics
    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("D4_CN_DIAG: Dispersion CN values (from GFNFFCN):");
        for (size_t i = 0; i < m_atoms.size(); ++i) {
            CurcumaLogger::param(fmt::format("Atom {} (Z={})", i, m_atoms[i]),
                                fmt::format("{:.6f}", m_cn_values[i]));
        }
    }

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::success(fmt::format("D4: Molecular CN calculated in {:.2f} ms", t_cn_ms));
    }

    // Phase B debug output: Complete CN values for all atoms (Dec 2025)
    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("=== D4 Coordination Numbers ===");
        for (size_t i = 0; i < m_atoms.size(); ++i) {
            CurcumaLogger::param(fmt::format("Atom {} (Z={})", i, m_atoms[i]),
                                 fmt::format("CN={:.4f}", m_cn_values[i]));
        }
    }

    // STEP 1: Get EEQ charges - reuse topology charges if available, else compute fresh
    // Claude Generated (Feb 7, 2026): Phase 5 optimization - skip redundant EEQ solve
    // GFN-FF already computed charges during topology initialization; reuse them for D4
    auto t_eeq_start = std::chrono::high_resolution_clock::now();
    double t_eeq_ms = 0.0;

    if (m_topology_charges.size() == static_cast<Eigen::Index>(m_atoms.size())) {
        // Reuse topology charges (avoids O(N³) matrix decomposition)
        m_eeq_charges = m_topology_charges;
        auto t_eeq_end = std::chrono::high_resolution_clock::now();
        t_eeq_ms = std::chrono::duration<double, std::milli>(t_eeq_end - t_eeq_start).count();

        if (CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::success(fmt::format("D4: Reused topology charges in {:.2f} ms (skipped EEQ solve)", t_eeq_ms));
        }
    } else {
        // Fallback: compute fresh EEQ charges
        Vector cn_eigen = Eigen::Map<const Vector>(m_cn_values.data(), m_cn_values.size());
        m_eeq_charges = m_eeq_solver->calculateCharges(m_atoms, geometry_bohr, 0, &cn_eigen);
        auto t_eeq_end = std::chrono::high_resolution_clock::now();
        t_eeq_ms = std::chrono::duration<double, std::milli>(t_eeq_end - t_eeq_start).count();

        if (CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::success(fmt::format("D4: EEQ charges calculated in {:.2f} ms", t_eeq_ms));
        }
    }

    if (m_eeq_charges.size() != static_cast<Eigen::Index>(m_atoms.size())) {
        CurcumaLogger::error("D4: EEQ charge calculation failed");
        return;
    }

    // Phase B debug output: Complete EEQ charges for all atoms (Dec 2025)
    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("=== D4 EEQ Charges ===");
        for (size_t i = 0; i < m_atoms.size(); ++i) {
            CurcumaLogger::param(fmt::format("Atom {} (Z={})", i, m_atoms[i]),
                                 fmt::format("q={:.6f}", m_eeq_charges(i)));
        }
    }

    // Claude Generated (Dec 27, 2025): Pre-compute Gaussian weights ONCE for all atoms
    // This eliminates redundant exp() calls in getChargeWeightedC6()
    // Performance: O(N×M) instead of O(N²×M) exp() calculations
    auto t_weights_start = std::chrono::high_resolution_clock::now();
    precomputeGaussianWeights();
    auto t_weights_end = std::chrono::high_resolution_clock::now();
    double t_weights_ms = std::chrono::duration<double, std::milli>(t_weights_end - t_weights_start).count();

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info(fmt::format("D4: Weight pre-computation took {:.2f} ms ({} atoms × avg {} refs)",
            t_weights_ms, m_atoms.size(), m_gaussian_weights.empty() ? 0 : m_gaussian_weights[0].size()));
    }

    // STEP 2: Generate C6 pairs with CN+charge weighted C6 coefficients
    auto t_pairs_start = std::chrono::high_resolution_clock::now();
    int num_pairs = 0;

    // Claude Generated (Jan 3, 2026): Debug pair generation
    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info(fmt::format("D4: Generating pairs for {} atoms (expect {} pairs)",
            m_atoms.size(), (m_atoms.size() * (m_atoms.size() - 1)) / 2));
    }

    // Claude Generated (February 2026): Phase 3 - OpenMP D4 Pair Loop Parallelization
    //
    // PROBLEM: D4 C6 interpolation for 1410 atoms = 993,345 pair calculations (serial)
    //          Current timing: ~2.7 seconds with zero parallelization
    //
    // SOLUTION: Parallelize nested i-j loop using OpenMP with collapse(2)
    //           Expected speedup: 3-4× on 4 cores (2.7s → ~0.7s)
    //
    // ARCHITECTURE:
    //   - collapse(2): Parallelize BOTH i and j loops for better load balancing
    //   - Dynamic scheduling: Handles triangular iteration (j > i) efficiently
    //   - Thread-local storage: Each thread builds its own pair list
    //   - Critical section: Minimal synchronization for merging
    //
    // THREAD SAFETY:
    //   - m_atoms, m_cn_values, m_eeq_charges, m_topology_charges: Read-only ✅
    //   - getChargeWeightedC6(), getSqrtZr4r2(): Pure functions with cached data ✅
    //   - GFNFFParameters::zetaChargeScale(): Static data ✅
    //   - local_pairs: Thread-local ✅
    //   - pairs_vec: Protected by critical section ✅

    std::vector<json> pairs_vec;

    #pragma omp parallel
    {
        std::vector<json> local_pairs;

        // NOTE: Cannot use collapse(2) with triangular loops (j = i + 1)
        // Parallelize outer loop only, still provides good speedup
        #pragma omp for schedule(dynamic, 10)
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

                // CLAUDE GENERATED (January 25, 2026): GFN-FF modified dispersion formula
                // Reference: gfnff_gdisp0.f90:365-377, gfnff_param.f90:531-532
                //
                // GFN-FF uses a MODIFIED BJ damping formula (NOT standard D3/D4):
                //   E = -0.5 * C6 * (t6 + 2*r4r2ij*t8)
                // where:
                //   t6 = 1/(r^6 + R0^6)
                //   t8 = 1/(r^8 + R0^8)
                //   r4r2ij = 3 * sqrtZr4r2_i * sqrtZr4r2_j (implicit C8/C6 ratio)
                //   R0^2 = (a1*sqrt(r4r2ij) + a2)^2 with a1=0.58, a2=4.80
                //   sqrtZr4r2 = sqrt(0.5 * r4Overr2 * sqrt(Z)) [pre-computed in m_sqrt_z_r4_r2]

                // CN-weighted C6 using Casimir-Polder integration (Dec 2025 Phase 2.2)
                double c6 = getChargeWeightedC6(atom_i, atom_j, i, j);

                // Claude Generated (Feb 8, 2026): C6 reference diagnostic for first few pairs (verbosity >= 3 only)
                if (CurcumaLogger::get_verbosity() >= 3 && i == 0 && j == 1) {
                    int nref_i = (atom_i > 0 && atom_i <= MAX_ELEM) ? m_refn[atom_i - 1] : 0;
                    int nref_j = (atom_j > 0 && atom_j <= MAX_ELEM) ? m_refn[atom_j - 1] : 0;
                    CurcumaLogger::info(fmt::format("D4_C6REF_DIAG: pair ({},{}) Z=({},{}) c6_weighted={:.6f} nref=({},{})",
                                                   i, j, atom_i, atom_j, c6, nref_i, nref_j));
                    for (int ri = 0; ri < std::min(nref_i, MAX_REF); ++ri) {
                        for (int rj = 0; rj < std::min(nref_j, MAX_REF); ++rj) {
                            double c6ref = m_c6_flat_cache[c6FlatIndex(atom_i-1, atom_j-1, ri, rj)];
                            double wi = m_gaussian_weights[i][ri];
                            double wj = m_gaussian_weights[j][rj];
                            if (wi * wj > 1e-6) {
                                CurcumaLogger::info(fmt::format("D4_C6REF_DIAG: ref({},{}) c6ref={:.4f} wi={:.4f} wj={:.4f} contrib={:.6f}",
                                                               ri, rj, c6ref, wi, wj, wi*wj*c6ref));
                            }
                        }
                    }
                }

                // GFN-FF specific parameters (NOT standard D3/D4!)
                // sqrtZr4r2 values from pre-computed m_sqrt_z_r4_r2 array
                double sqrt_zr4r2_i = getSqrtZr4r2(atom_i);
                double sqrt_zr4r2_j = getSqrtZr4r2(atom_j);

                // r4r2ij = 3 * sqrtZr4r2_i * sqrtZr4r2_j (Fortran: gfnff_gdisp0.f90:365)
                double r4r2ij = 3.0 * sqrt_zr4r2_i * sqrt_zr4r2_j;

                // R0^2 = (a1*sqrt(r4r2ij) + a2)^2 (Fortran: gfnff_param.f90:532)
                // GFN-FF constants: a1=0.58, a2=4.80 (from gfnff_param.f90:841-842)
                double a1 = m_config.get<double>("d4_a1", 0.58);
                double a2 = m_config.get<double>("d4_a2", 4.80);
                double r0_squared = std::pow(a1 * std::sqrt(r4r2ij) + a2, 2);

                // Legacy C8 calculation (for backward compatibility, NOT used in GFN-FF formula)
                double c8 = 3.0 * c6 * sqrt_zr4r2_i * sqrt_zr4r2_j;

                // Claude Generated (Jan 31, 2026): GFN-FF zeta charge scaling
                // Reference: gfnff_ini.f90:789-806, gfnff_gdisp0.f90:374
                // The zeta function provides charge-dependent C6 scaling:
                //   zetac6_ij = zeta(Z_i, q_i) * zeta(Z_j, q_j)
                //
                // CRITICAL FIX (Jan 31, 2026): Use TOPOLOGY charges (topo%qa) for zeta scaling
                // Reference: Fortran gfnff_ini.f90:789 - f1 = zeta(ati, topo%qa(i))
                // GFN-FF computes zetac6 ONCE during initialization using topology charges,
                // which are calculated with INTEGER neighbor counts (neighbor_count),
                // NOT the fractional CN from geometry-dependent EEQ.
                double q_i, q_j;
                if (m_topology_charges.size() > 0) {
                    // Use topology charges (topo%qa equivalent) for zeta scaling
                    q_i = (i < static_cast<size_t>(m_topology_charges.size())) ? m_topology_charges(i) : 0.0;
                    q_j = (j < static_cast<size_t>(m_topology_charges.size())) ? m_topology_charges(j) : 0.0;
                    if (CurcumaLogger::get_verbosity() >= 3 && i == 0 && j == 1) {
                        CurcumaLogger::info("Zeta scaling: Using topology charges (topo%qa)");
                    }
                } else {
                    // Fallback to EEQ charges (geometry-dependent, less accurate for GFN-FF)
                    q_i = (i < static_cast<size_t>(m_eeq_charges.size())) ? m_eeq_charges(i) : 0.0;
                    q_j = (j < static_cast<size_t>(m_eeq_charges.size())) ? m_eeq_charges(j) : 0.0;
                    if (CurcumaLogger::get_verbosity() >= 3 && i == 0 && j == 1) {
                        CurcumaLogger::warn("Zeta scaling: Falling back to EEQ charges (no topology charges)");
                    }
                }
                double zeta_i = GFNFFParameters::zetaChargeScale(atom_i, q_i);
                double zeta_j = GFNFFParameters::zetaChargeScale(atom_j, q_j);
                double zetac6 = zeta_i * zeta_j;

                // GFN-FF dispersion pair parameters
                pair["dispersion_method"] = "d4";  // Route to native GFN-FF in ForceFieldThread
                pair["C6"] = c6;              // CN-weighted C6 from Casimir-Polder
                pair["r4r2ij"] = r4r2ij;      // GFN-FF: implicit C8/C6 factor
                pair["r0_squared"] = r0_squared;  // GFN-FF: pre-computed (a1*sqrt(r4r2ij)+a2)^2
                pair["zetac6"] = zetac6;      // GFN-FF: charge-dependent zeta scaling (Jan 31, 2026)

                // Legacy fields (for backward compatibility with standard D3/D4)
                pair["C8"] = c8;
                pair["s6"] = m_config.get<double>("d4_s6", 1.0);
                pair["s8"] = m_config.get<double>("d4_s8", 2.0);
                pair["a1"] = a1;
                pair["a2"] = a2;
                pair["r_cut"] = 100.0;  // Cutoff radius (Bohr)

                local_pairs.push_back(pair);
            } else {
                if (CurcumaLogger::get_verbosity() >= 2) {
                    CurcumaLogger::warn("D4: Elements out of range - i=" +
                                       std::to_string(atom_i) +
                                       " j=" + std::to_string(atom_j));
                }
            }
        }  // end j loop
    }  // end i loop

        // Merge thread-local results into global container
        #pragma omp critical
        {
            pairs_vec.insert(pairs_vec.end(), local_pairs.begin(), local_pairs.end());
        }
    }  // end omp parallel

    // Convert vector to JSON array and count pairs
    json dispersion_pairs = json::array();
    for (const auto& p : pairs_vec) {
        dispersion_pairs.push_back(p);
        num_pairs++;
    }

    auto t_pairs_end = std::chrono::high_resolution_clock::now();
    double t_pairs_ms = std::chrono::duration<double, std::milli>(t_pairs_end - t_pairs_start).count();

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info(fmt::format("D4: C6 interpolation for {} pairs took {:.2f} ms ({:.4f} ms/pair)",
            num_pairs, t_pairs_ms, num_pairs > 0 ? t_pairs_ms / num_pairs : 0.0));
    }

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::param("Generated D4 pairs", static_cast<int>(dispersion_pairs.size()));

        // Claude Generated (Mar 5, 2026): D4 diagnostics to JSON file
        {
            json diag;
            diag["type"] = "d4_dispersion";
            diag["n_pairs"] = static_cast<int>(dispersion_pairs.size());
            double sum_c6 = 0.0, sum_zetac6 = 0.0, sum_c6_zetac6 = 0.0;
            for (const auto& pair : dispersion_pairs) {
                double c6 = pair["C6"];
                double zetac6 = pair.value("zetac6", 1.0);
                sum_c6 += c6;
                sum_zetac6 += zetac6;
                sum_c6_zetac6 += c6 * zetac6;
            }
            diag["sum_C6"] = sum_c6;
            diag["sum_zetac6"] = sum_zetac6;
            diag["sum_C6_x_zetac6"] = sum_c6_zetac6;

            // Per-atom CN and zeta values
            json atoms = json::array();
            for (size_t i = 0; i < m_atoms.size(); ++i) {
                double q_i = (i < static_cast<size_t>(m_topology_charges.size())) ? m_topology_charges(i) : 0.0;
                double zeta_i = GFNFFParameters::zetaChargeScale(m_atoms[i], q_i);
                double cn_i = (i < m_cn_values.size()) ? m_cn_values[i] : 0.0;
                atoms.push_back({{"idx", static_cast<int>(i)}, {"Z", m_atoms[i]},
                                {"qa", q_i}, {"zeta", zeta_i}, {"CN", cn_i}});
            }
            diag["atoms"] = atoms;

            std::ofstream diag_file("gfnff_diag_d4.json");
            if (diag_file.is_open()) {
                diag_file << diag.dump(2) << std::endl;
                diag_file.close();
                CurcumaLogger::info("Wrote D4 dispersion diagnostics to gfnff_diag_d4.json");
            }
        }
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
            CurcumaLogger::success(fmt::format("D4 C6: {} pairs, avg={:.4f}, min={:.4f}, max={:.4f}, zeros={}",
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

    if (m_config.get<double>("d4_s9", 1.0) > 1e-10) {  // Claude Generated: Use PARAM default (1.0) not 0.0
        double s9 = m_config.get<double>("d4_s9", 1.0);
        double a1 = m_config.get<double>("d4_a1", 0.44);
        double a2 = m_config.get<double>("d4_a2", 4.60);
        double alp = m_config.get<double>("d4_alp", 14.0);

        int n_atoms = static_cast<int>(m_atoms.size());

        // Claude Generated (January 2025): Bonded triplet filtering for D4 (like XTB)
        // This optimization was reverted and is now being restored.
        // It reduces the number of triplets from O(N^3) to roughly O(N*<neighbors>^2),
        // which is a massive performance improvement for larger systems.
        std::vector<std::pair<int, int>> bonded_pairs;
        std::vector<std::vector<int>> adjacency(n_atoms);

        // Covalent radii in Angstrom (from GFN-FF method) for bond detection
        const double ANGSTROM_TO_BOHR = 1.8897261246257702;
        static const std::vector<double> rcov_angstrom = {
            0.32, 0.37, 1.30, 0.99, 0.84, 0.75, 0.71, 0.64, 0.60, 0.62, // H-Ne
            1.60, 1.40, 1.24, 1.14, 1.09, 1.04, 1.00, 1.01, // Na-Ar
            2.00, 1.74, 1.59, 1.48, 1.44, 1.30, 1.29, 1.24, 1.18, 1.17, 1.22, 1.20, // K-Zn
            1.23, 1.20, 1.20, 1.18, 1.17, 1.16, 2.15, 1.90, 1.76, 1.64, 1.56, 1.46, // Rb-Tc
            1.38, 1.36, 1.34, 1.30, 1.36, 1.40, 1.42, 1.40, 1.40, 1.37, 1.36, 1.36, // Ru-Xe
            2.38, 2.06, 1.94, 1.84, 1.90, 1.88, 1.86, 1.85, 1.83, 1.82, 1.81, 1.80, // Cs-Er
            1.79, 1.77, 1.77, 1.78, 1.74, 1.64, 1.58, 1.50, 1.41, 1.36, 1.32, 1.30, // Tm-Au
            1.30, 1.32, 1.44, 1.45, 1.50, 1.42, 1.48, 1.46 // Hg-Rn
        };

        for (int i = 0; i < n_atoms; ++i) {
            int Zi = m_atoms[i];
            double rcov_i = (Zi > 0 && Zi < rcov_angstrom.size()) ? rcov_angstrom[Zi-1] : 1.5;

            for (int j = i + 1; j < n_atoms; ++j) {
                int Zj = m_atoms[j];
                double rcov_j = (Zj > 0 && Zj < rcov_angstrom.size()) ? rcov_angstrom[Zj-1] : 1.5;

                double distance_bohr = (geometry_bohr.row(i) - geometry_bohr.row(j)).norm();
                // Bond threshold is defined in Angstrom, so compare in Angstrom
                double distance_angstrom = distance_bohr / ANGSTROM_TO_BOHR;
                double bond_threshold = 1.3 * (rcov_i + rcov_j);

                if (distance_angstrom < bond_threshold) {
                    bonded_pairs.push_back({i, j});
                    adjacency[i].push_back(j);
                    adjacency[j].push_back(i);
                }
            }
        }
        
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::param("D4 ATM bonded pairs identified", static_cast<int>(bonded_pairs.size()));
        }

        std::set<std::tuple<int, int, int>> unique_triplets;

        // Second pass: generate only bonded triplets
        for (const auto& bond : bonded_pairs) {
            int i = bond.first;
            int j = bond.second;

            // Find all neighbors k of i or j, forming a connected triplet
            std::set<int> neighbors_of_i_and_j;
            for (int k : adjacency[i]) {
                if (k != j) neighbors_of_i_and_j.insert(k);
            }
            for (int k : adjacency[j]) {
                if (k != i) neighbors_of_i_and_j.insert(k);
            }

            for (int k : neighbors_of_i_and_j) {
                // Ensure unique triplet ordering: i < j < k
                std::array<int, 3> t_sorted = {i, j, k};
                std::sort(t_sorted.begin(), t_sorted.end());
                unique_triplets.insert({t_sorted[0], t_sorted[1], t_sorted[2]});
            }
        }
        
        // Claude Generated (Feb 7, 2026): Phase 6 - Convert set to vector for OpenMP indexed iteration
        std::vector<std::tuple<int, int, int>> triplet_vec(unique_triplets.begin(), unique_triplets.end());

        // Claude Generated (Feb 7, 2026): Phase 6 - OpenMP parallelize ATM triple loop
        // Same proven pattern as C6 pair loop above: thread-local vectors + critical merge
        // Thread safety: m_atoms, m_gaussian_weights, m_c6_reference_cache are read-only
        //                getChargeWeightedC6() and calculateTripleScale() are pure functions
        #pragma omp parallel
        {
            std::vector<json> local_triples;

            #pragma omp for schedule(dynamic, 10)
            for (size_t idx = 0; idx < triplet_vec.size(); ++idx) {
                int i = std::get<0>(triplet_vec[idx]);
                int j = std::get<1>(triplet_vec[idx]);
                int k = std::get<2>(triplet_vec[idx]);

                json triple;
                triple["i"] = i;
                triple["j"] = j;
                triple["k"] = k;

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
                triple["triple_scale"] = calculateTripleScale(i, j, k);

                local_triples.push_back(std::move(triple));
            }

            #pragma omp critical
            {
                for (auto& t : local_triples) {
                    atm_triples.push_back(std::move(t));
                }
            }
        }

        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::param("Generated D4 ATM triples (bonded)", static_cast<int>(atm_triples.size()));
        }
    }

    m_parameters["atm_triples"] = atm_triples;

    auto t_end_total = std::chrono::high_resolution_clock::now();
    double t_total_ms = std::chrono::duration<double, std::milli>(t_end_total - t_start_total).count();

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::success(fmt::format("D4 parameter generation completed in {:.2f} ms", t_total_ms));
        CurcumaLogger::info(fmt::format("  └─ Breakdown: CN={:.1f}ms, EEQ={:.1f}ms, Weights={:.1f}ms, C6Interp={:.1f}ms",
            t_cn_ms, t_eeq_ms, t_weights_ms, t_pairs_ms));
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

void D4ParameterGenerator::precomputeC6ReferenceMatrix()
{
    // Claude Generated (Dec 27, 2025): Pre-compute C6 reference values via Casimir-Polder integration
    // Reference: XTB dftd4.F90 lines ~500: c6 = thopi * trapzd(alpha_i * alpha_j)
    //
    // Performance optimization: Eliminates redundant Casimir-Polder integrations
    // - Old approach: O(N²×M²×F) integrations (for every atom pair × ref states × frequencies)
    // - New approach: O(E²×M²×F) integrations (once per element-pair combination)
    // For triose (66 atoms): 2145 pairs × 49 ref combinations = 105,000 → ~1000 unique combinations
    // Expected speedup: 50-100x for large molecules

    auto t_start = std::chrono::high_resolution_clock::now();

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("=== D4 C6 Reference Matrix Pre-computation (PHASE 1 OPTIMIZED) ===");
    }

    // Frequency grid for Casimir-Polder integration
    const std::vector<double> frequency_grid = {
        0.000001, 0.050000, 0.100000, 0.200000, 0.300000, 0.400000,
        0.500000, 0.600000, 0.700000, 0.800000, 0.900000, 1.000000,
        1.200000, 1.400000, 1.600000, 1.800000, 2.000000, 2.500000,
        3.000000, 4.000000, 5.000000, 7.500000, 10.000000
    };

    int computed_count = 0;

    // Claude Generated (Feb 7, 2026): Phase 7a - Initialize flat dense array
    // Size: MAX_ELEM * MAX_ELEM * MAX_REF * MAX_REF = 118*118*7*7 = 406,952 entries (~3.1 MB)
    const size_t flat_size = static_cast<size_t>(MAX_ELEM) * MAX_ELEM * MAX_REF * MAX_REF;
    m_c6_flat_cache.assign(flat_size, 0.0);

    // Pre-compute C6 for all element-pair combinations that have alphaiw data
    for (int elem_i = 0; elem_i < MAX_ELEM && elem_i < static_cast<int>(d4_alphaiw_data.size()); ++elem_i) {
        int nref_i = (elem_i < static_cast<int>(m_refn.size())) ? m_refn[elem_i] : 0;
        if (nref_i == 0 || elem_i >= static_cast<int>(d4_alphaiw_data.size())) continue;

        for (int elem_j = 0; elem_j <= elem_i; ++elem_j) {  // Symmetric, only compute lower triangle
            int nref_j = (elem_j < static_cast<int>(m_refn.size())) ? m_refn[elem_j] : 0;
            if (nref_j == 0 || elem_j >= static_cast<int>(d4_alphaiw_data.size())) continue;

            // Compute C6 for all reference state combinations
            for (int ref_i = 0; ref_i < nref_i && ref_i < MAX_REF; ++ref_i) {
                if (ref_i >= static_cast<int>(d4_alphaiw_data[elem_i].size())) continue;

                for (int ref_j = 0; ref_j < nref_j && ref_j < MAX_REF; ++ref_j) {
                    if (ref_j >= static_cast<int>(d4_alphaiw_data[elem_j].size())) continue;

                    double c6 = computeC6Reference(elem_i, elem_j, ref_i, ref_j);

                    // Store in flat array (symmetric storage)
                    m_c6_flat_cache[c6FlatIndex(elem_i, elem_j, ref_i, ref_j)] = c6;
                    if (elem_i != elem_j || ref_i != ref_j) {
                        m_c6_flat_cache[c6FlatIndex(elem_j, elem_i, ref_j, ref_i)] = c6;
                    }

                    computed_count++;
                }
            }
        }
    }

    m_c6_reference_cached = true;

    auto t_end = std::chrono::high_resolution_clock::now();
    double t_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::success(fmt::format("D4: Pre-computed {} C6 reference values in {:.2f} ms",
            computed_count, t_ms));
        CurcumaLogger::info(fmt::format("  Cache size: {} entries (flat array, {:.1f} MB)",
            flat_size, flat_size * sizeof(double) / (1024.0 * 1024.0)));
    }
}

double D4ParameterGenerator::computeC6Reference(int elem_i, int elem_j, int ref_i, int ref_j) const
{
    // Claude Generated (Dec 27, 2025): Casimir-Polder integration for a single element-pair × ref-state combination
    // Reference: XTB dftd4.F90 lines ~500: c6 = thopi * trapzd(alpha_i * alpha_j)
    //
    // Formula: C6_ij = (3/π) ∫ α_i(iω) * α_j(iω) dω

    // Frequency grid for Casimir-Polder integration
    const std::vector<double> frequency_grid = {
        0.000001, 0.050000, 0.100000, 0.200000, 0.300000, 0.400000,
        0.500000, 0.600000, 0.700000, 0.800000, 0.900000, 1.000000,
        1.200000, 1.400000, 1.600000, 1.800000, 2.000000, 2.500000,
        3.000000, 4.000000, 5.000000, 7.500000, 10.000000
    };

    double c6_ref = 0.0;

    // Check if we have frequency-dependent polarizabilities for both elements
    if (elem_i >= static_cast<int>(d4_alphaiw_data.size()) ||
        elem_j >= static_cast<int>(d4_alphaiw_data.size()) ||
        ref_i >= static_cast<int>(d4_alphaiw_data[elem_i].size()) ||
        ref_j >= static_cast<int>(d4_alphaiw_data[elem_j].size())) {
        return 1.0;  // Fallback for elements without polarizability data
    }

    // Get correction factors for reference states
    double ascale_i = (elem_i < static_cast<int>(d4_ascale_data.size()) && ref_i < static_cast<int>(d4_ascale_data[elem_i].size()))
                    ? d4_ascale_data[elem_i][ref_i] : 1.0;
    double ascale_j = (elem_j < static_cast<int>(d4_ascale_data.size()) && ref_j < static_cast<int>(d4_ascale_data[elem_j].size()))
                    ? d4_ascale_data[elem_j][ref_j] : 1.0;

    double hcount_i = (elem_i < static_cast<int>(m_refh.size()) && ref_i < static_cast<int>(m_refh[elem_i].size()))
                    ? m_refh[elem_i][ref_i] : 0.0;
    double hcount_j = (elem_j < static_cast<int>(m_refh.size()) && ref_j < static_cast<int>(m_refh[elem_j].size()))
                    ? m_refh[elem_j][ref_j] : 0.0;

    int refsys_i = (elem_i < static_cast<int>(d4_refsys_data.size()) && ref_i < static_cast<int>(d4_refsys_data[elem_i].size()))
                 ? d4_refsys_data[elem_i][ref_i] : 0;
    int refsys_j = (elem_j < static_cast<int>(d4_refsys_data.size()) && ref_j < static_cast<int>(d4_refsys_data[elem_j].size()))
                 ? d4_refsys_data[elem_j][ref_j] : 0;

    double sscale_i = (d4_sscale_data.find(refsys_i) != d4_sscale_data.end()) ? d4_sscale_data.at(refsys_i) : 0.0;
    double sscale_j = (d4_sscale_data.find(refsys_j) != d4_sscale_data.end()) ? d4_sscale_data.at(refsys_j) : 0.0;

    // Debug output for ascale and hcount verification (Phase A - Dec 2025)
    static int debug_hcount = 0;
    if (CurcumaLogger::get_verbosity() >= 3 && debug_hcount < 5) {
        CurcumaLogger::info(fmt::format("D4 params: elem_i={} ref_i={} ascale={:.4f} hcount={:.1f}, elem_j={} ref_j={} ascale={:.4f} hcount={:.1f}",
                                         elem_i, ref_i, ascale_i, hcount_i, elem_j, ref_j, ascale_j, hcount_j));
        debug_hcount++;
    }

    // Integrate product of CORRECTED polarizabilities using trapezoidal rule
    // Correction formula: α_corrected = ascale * (αᵢⱼw - hcount * sscale * secaiw)
    for (int iw = 0; iw < N_FREQ - 1; ++iw) {
        // Get raw alphaiw values
        double alphaiw_i_iw = d4_alphaiw_data[elem_i][ref_i][iw];
        double alphaiw_i_next = d4_alphaiw_data[elem_i][ref_i][iw + 1];
        double alphaiw_j_iw = d4_alphaiw_data[elem_j][ref_j][iw];
        double alphaiw_j_next = d4_alphaiw_data[elem_j][ref_j][iw + 1];

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

    return c6_ref;
}

void D4ParameterGenerator::precomputeGaussianWeights(int num_threads)
{
    // Claude Generated (Dec 27, 2025): Pre-compute CN-only Gaussian weights once per atom
    // Claude Generated (Mar 2026): Internal std::thread parallelisation
    // Reference: external/gfnff/src/gfnff_gdisp0.f90:405 weight_cn() function

    constexpr double wf = 4.0;  // Gaussian width parameter (from D4 cpp-d4)

    int natoms_gw = static_cast<int>(m_atoms.size());
    m_gaussian_weights.resize(natoms_gw);

    // Worker lambda: each thread processes interleaved atom indices
    auto gw_worker = [&](int t_id, int T) {
    for (int i = t_id; i < natoms_gw; i += T) {
        int Zi = m_atoms[i];
        int elem_i = Zi - 1;  // Convert to 0-based

        // Validate element range
        if (elem_i < 0 || elem_i >= MAX_ELEM) {
            m_gaussian_weights[i].clear();
            continue;
        }

        // Get number of reference states for this element
        int nref = (elem_i < static_cast<int>(m_refn.size())) ? m_refn[elem_i] : 1;

        if (nref == 0) {
            m_gaussian_weights[i].clear();
            continue;
        }

        // Get atom properties
        // EEQ charges calculated but NOT used in GFN-FF weighting (CN-only model)
        // Charges remain available for potential future full D4 implementation
        double cni = m_cn_values[i];

        // Compute Gaussian weights for all reference states of this atom
        std::vector<double> weights(nref, 0.0);
        double sum_weights = 0.0;

        for (int ref = 0; ref < nref && ref < MAX_REF; ++ref) {
            double cni_ref = (elem_i < static_cast<int>(m_refcn.size()) &&
                             ref < static_cast<int>(m_refcn[elem_i].size()))
                            ? m_refcn[elem_i][ref] : 0.0;

            // KEY: CN-only Gaussian weighting (matches GFN-FF hybrid model)
            // Reference: external/gfnff/src/gfnff_gdisp0.f90:405 - cngw = exp(-wf * (cn - cnref)**2)
            double diff_cn = cni - cni_ref;
            weights[ref] = std::exp(-wf * diff_cn * diff_cn);
            sum_weights += weights[ref];
        }

        // Normalize weights
        if (sum_weights > 1e-10) {
            for (int ref = 0; ref < nref; ++ref) {
                weights[ref] /= sum_weights;
            }
        } else {
            // Exceptional case: set first reference to 1.0 (neutral state fallback)
            std::fill(weights.begin(), weights.end(), 0.0);
            weights[0] = 1.0;
        }

        m_gaussian_weights[i] = std::move(weights);
    }
    };  // end gw_worker lambda

    // Dispatch: parallel if beneficial, otherwise single-threaded
    if (num_threads > 1 && natoms_gw > 64) {
        int T = std::min(num_threads, natoms_gw);
        std::vector<std::thread> threads(T - 1);
        for (int t = 1; t < T; ++t)
            threads[t - 1] = std::thread(gw_worker, t, T);
        gw_worker(0, T);
        for (auto& th : threads) th.join();
    } else {
        gw_worker(0, 1);
    }

    // Claude Generated (Feb 8, 2026): Diagnostic for Gaussian weights (verbosity >= 3 only)
    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("D4_GW_DIAG: Gaussian weights per atom:");
        for (size_t i = 0; i < m_atoms.size(); ++i) {
            int elem_i = m_atoms[i];
            int nref_i = (elem_i > 0 && elem_i <= MAX_ELEM && (elem_i - 1) < static_cast<int>(m_refn.size()))
                         ? m_refn[elem_i - 1] : 0;
            std::string gw_str = fmt::format("atom {} (Z={}) CN={:.4f} nref={}", i, elem_i, m_cn_values[i], nref_i);
            for (size_t ref = 0; ref < m_gaussian_weights[i].size() && ref < 7; ++ref) {
                double cnref = ((elem_i - 1) < static_cast<int>(m_refcn.size()) &&
                               ref < m_refcn[elem_i - 1].size())
                              ? m_refcn[elem_i - 1][ref] : -1.0;
                gw_str += fmt::format(" w[{}]={:.4f}(cnref={:.2f})", ref, m_gaussian_weights[i][ref], cnref);
            }
            CurcumaLogger::param("D4_GW_DIAG", gw_str);
        }
    }

    // Claude Generated (Feb 7, 2026): Phase 7b - Pre-compute dominant reference indices
    // All refs with nonzero weight contribute to accuracy (Fortran uses all refs)
    constexpr double WEIGHT_THRESHOLD = 0.0;
    m_dominant_refs.resize(m_atoms.size());
    for (size_t i = 0; i < m_atoms.size(); ++i) {
        m_dominant_refs[i].clear();
        for (size_t ref = 0; ref < m_gaussian_weights[i].size(); ++ref) {
            if (m_gaussian_weights[i][ref] > WEIGHT_THRESHOLD) {
                m_dominant_refs[i].push_back(static_cast<int>(ref));
            }
        }
    }

    m_weights_cached = true;

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("=== D4 Gaussian Weights Pre-computed ===");
        CurcumaLogger::param("Atoms processed", static_cast<int>(m_atoms.size()));
        for (size_t i = 0; i < std::min<size_t>(m_atoms.size(), 231); ++i) {
            std::string weights_str = "";
            for (size_t ref = 0; ref < m_gaussian_weights[i].size(); ++ref) {
                weights_str += fmt::format("{:.6f}", m_gaussian_weights[i][ref]);
                if (ref < m_gaussian_weights[i].size() - 1) weights_str += ", ";
            }
            CurcumaLogger::info(fmt::format("  Atom {} (Z={}, q={:.4f}, CN={:.4f}): weights=[{}]",
                i, m_atoms[i], m_eeq_charges(i), m_cn_values[i], weights_str));
        }
    }
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
 * Claude Generated - December 27, 2025: Optimized to use cached weights
 *
 * @param Zi Atomic number of atom i
 * @param Zj Atomic number of atom j
 * @param atom_i Atom index i (for weight lookup)
 * @param atom_j Atom index j (for weight lookup)
 * @return Charge-weighted C6 coefficient (Hartree * Bohr^6)
 */
double D4ParameterGenerator::getChargeWeightedC6(int Zi, int Zj, size_t atom_i, size_t atom_j) const
{
    // Claude Generated (Feb 7, 2026): Phase 7a/7b/7c - Fully optimized C6 interpolation
    // - 7a: Dense flat array (no hash lookups)
    // - 7b: Dominant refs only (skip near-zero weights)
    // - 7c: No debug logging in hot path (moved to caller if needed)

    // Convert to 0-based indexing
    int elem_i = Zi - 1;
    int elem_j = Zj - 1;

    // Validate (rare path)
    if (elem_i < 0 || elem_i >= MAX_ELEM || elem_j < 0 || elem_j >= MAX_ELEM) {
        return 0.0;
    }

    const auto& weights_i = m_gaussian_weights[atom_i];
    const auto& weights_j = m_gaussian_weights[atom_j];

    if (weights_i.empty() || weights_j.empty()) {
        return 0.0;
    }

    // Phase 7b: Use dominant refs if available, else iterate all refs
    const auto& refs_i = (atom_i < m_dominant_refs.size() && !m_dominant_refs[atom_i].empty())
                         ? m_dominant_refs[atom_i]
                         : std::vector<int>();  // empty = iterate all
    const auto& refs_j = (atom_j < m_dominant_refs.size() && !m_dominant_refs[atom_j].empty())
                         ? m_dominant_refs[atom_j]
                         : std::vector<int>();

    double c6_weighted = 0.0;

    // Claude Generated (Feb 8, 2026): DEBUG: Log C6 computation for small molecules (verbosity >= 3 only)
    bool log_c6 = (CurcumaLogger::get_verbosity() >= 3 && atom_i == 0 && atom_j <= 1);

    // Base offset for elem_i and elem_j in flat cache
    const size_t base_ij = static_cast<size_t>(elem_i) * MAX_ELEM * MAX_REF * MAX_REF
                         + static_cast<size_t>(elem_j) * MAX_REF * MAX_REF;

    if (!refs_i.empty() && !refs_j.empty()) {
        // Fast path: only iterate dominant references (Phase 7b)
        for (int ri : refs_i) {
            double wi = weights_i[ri];
            size_t base_ri = base_ij + static_cast<size_t>(ri) * MAX_REF;
            for (int rj : refs_j) {
                double c6_ref = m_c6_flat_cache[base_ri + rj];
                double contrib = wi * weights_j[rj] * c6_ref;
                c6_weighted += contrib;
                if (log_c6) {
                    CurcumaLogger::info(fmt::format("C6_DEBUG: atom_i={} elem_i={} atom_j={} elem_j={} CN_i={:.4f} CN_j={:.4f}"
                                                     " ref_i={} ref_j={}: w_i={:.6f} w_j={:.6f} C6_ref={:.6f} contrib={:.6f}",
                                                     atom_i, Zi, atom_j, Zj, m_cn_values[atom_i], m_cn_values[atom_j],
                                                     ri, rj, wi, weights_j[rj], c6_ref, contrib));
                }
            }
        }
    } else {
        // Fallback: iterate all references
        int nref_i = static_cast<int>(weights_i.size());
        int nref_j = static_cast<int>(weights_j.size());
        for (int ri = 0; ri < nref_i; ++ri) {
            double wi = weights_i[ri];
            if (wi < 1e-12) continue;  // Micro-optimization: skip zero weights
            size_t base_ri = base_ij + static_cast<size_t>(ri) * MAX_REF;
            for (int rj = 0; rj < nref_j; ++rj) {
                double c6_ref = m_c6_flat_cache[base_ri + rj];
                double contrib = wi * weights_j[rj] * c6_ref;
                c6_weighted += contrib;
                if (log_c6) {
                    CurcumaLogger::info(fmt::format("C6_DEBUG: atom_i={} elem_i={} atom_j={} elem_j={} CN_i={:.4f} CN_j={:.4f}"
                                                     " ref_i={} ref_j={}: w_i={:.6f} w_j={:.6f} C6_ref={:.6f} contrib={:.6f}",
                                                     atom_i, Zi, atom_j, Zj, m_cn_values[atom_i], m_cn_values[atom_j],
                                                     ri, rj, wi, weights_j[rj], c6_ref, contrib));
                }
            }
        }
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

// Claude Generated (Feb 15, 2026): Compute derivatives of normalized Gaussian weights w.r.t. CN
// Reference: Fortran gfnff_gdisp0.f90:174-210 (weight_references_d4 subroutine)
//
// For normalized weight gw(ref) = expw(ref) / norm, the derivative is:
//   dgw/dCN = (d(expw)/dCN * norm - expw * d(norm)/dCN) / norm^2   (quotient rule)
// where:
//   expw(ref) = exp(-wf * (CN - CN_ref)^2)
//   d(expw)/dCN = -2*wf*(CN - CN_ref) * expw  = 2*wf*(CN_ref - CN) * expw
//   norm = sum_ref expw(ref)
//   d(norm)/dCN = sum_ref d(expw(ref))/dCN
void D4ParameterGenerator::computeGaussianWeightDerivatives(int num_threads)
{
    // Claude Generated (Mar 2026): Internal std::thread parallelisation
    constexpr double wf = 4.0;  // Gaussian width (matches precomputeGaussianWeights)

    int natoms_gwd = static_cast<int>(m_atoms.size());
    m_gaussian_weight_derivatives.resize(natoms_gwd);

    auto gwd_worker = [&](int t_id, int T) {
    for (int i = t_id; i < natoms_gwd; i += T) {
        int Zi = m_atoms[i];
        int elem_i = Zi - 1;

        if (elem_i < 0 || elem_i >= MAX_ELEM) {
            m_gaussian_weight_derivatives[i].clear();
            continue;
        }

        int nref = (elem_i < static_cast<int>(m_refn.size())) ? m_refn[elem_i] : 1;
        if (nref == 0) {
            m_gaussian_weight_derivatives[i].clear();
            continue;
        }

        double cni = m_cn_values[i];

        std::vector<double> expw(nref, 0.0);
        std::vector<double> dexpw(nref, 0.0);
        double norm = 0.0;
        double dnorm = 0.0;

        for (int ref = 0; ref < nref && ref < MAX_REF; ++ref) {
            double cni_ref = (elem_i < static_cast<int>(m_refcn.size()) &&
                             ref < static_cast<int>(m_refcn[elem_i].size()))
                            ? m_refcn[elem_i][ref] : 0.0;

            double diff_cn = cni - cni_ref;
            expw[ref] = std::exp(-wf * diff_cn * diff_cn);
            dexpw[ref] = 2.0 * wf * (cni_ref - cni) * expw[ref];
            norm += expw[ref];
            dnorm += dexpw[ref];
        }

        std::vector<double> dgwdcn(nref, 0.0);
        if (norm > 1e-10) {
            double norm2 = norm * norm;
            for (int ref = 0; ref < nref; ++ref) {
                dgwdcn[ref] = (dexpw[ref] * norm - expw[ref] * dnorm) / norm2;
            }
        }

        m_gaussian_weight_derivatives[i] = std::move(dgwdcn);
    }
    };  // end gwd_worker lambda

    if (num_threads > 1 && natoms_gwd > 64) {
        int T = std::min(num_threads, natoms_gwd);
        std::vector<std::thread> threads(T - 1);
        for (int t = 1; t < T; ++t)
            threads[t - 1] = std::thread(gwd_worker, t, T);
        gwd_worker(0, T);
        for (auto& th : threads) th.join();
    } else {
        gwd_worker(0, 1);
    }
}

// Claude Generated (Feb 15, 2026): Compute dc6dcn matrix from weight derivatives and C6 references
// Reference: Fortran gfnff_gdisp0.f90:262-305 (get_atomic_c6_d4 dc6dcn part)
//
// dc6dcn(i,j) = dC6(i,j)/dCN(i) = sum_{ri,rj} dgwdcn(ri,i) * gw(rj,j) * c6ref(ri,rj,Zi,Zj)
// dc6dcn(j,i) = dC6(i,j)/dCN(j) = sum_{ri,rj} gw(ri,i) * dgwdcn(rj,j) * c6ref(ri,rj,Zi,Zj)
void D4ParameterGenerator::computeDC6DCN(int num_threads)
{
    // Claude Generated (Mar 2026): Internal std::thread parallelisation for O(N²×M²)
    int natoms = static_cast<int>(m_atoms.size());
    m_dc6dcn = Matrix::Zero(natoms, natoms);

    // Thread safety proof: Each (i,j) pair with j >= i is processed by exactly the thread
    // owning row i. Writes: m_dc6dcn(i,j) = dc6_di (row i, thread's own row) and
    // m_dc6dcn(j,i) = dc6_dj (row j, column i). No other thread writes to (j,i) because
    // the pair (i,j) is only processed once by the thread owning row i.
    // Therefore direct writes to m_dc6dcn are safe without buffering.
    auto dc6_worker = [&](int t_id, int T) {
    for (int i = t_id; i < natoms; i += T) {
        int Zi = m_atoms[i];
        int elem_i = Zi - 1;
        if (elem_i < 0 || elem_i >= MAX_ELEM) continue;

        const auto& gw_i = m_gaussian_weights[i];
        const auto& dgw_i = m_gaussian_weight_derivatives[i];
        if (gw_i.empty() || dgw_i.empty()) continue;

        int nref_i = static_cast<int>(gw_i.size());

        for (int j = i; j < natoms; ++j) {
            int Zj = m_atoms[j];
            int elem_j = Zj - 1;
            if (elem_j < 0 || elem_j >= MAX_ELEM) continue;

            const auto& gw_j = m_gaussian_weights[j];
            const auto& dgw_j = m_gaussian_weight_derivatives[j];
            if (gw_j.empty() || dgw_j.empty()) continue;

            int nref_j = static_cast<int>(gw_j.size());

            const size_t base_ij = static_cast<size_t>(elem_i) * MAX_ELEM * MAX_REF * MAX_REF
                                 + static_cast<size_t>(elem_j) * MAX_REF * MAX_REF;

            double dc6_di = 0.0;
            double dc6_dj = 0.0;

            for (int ri = 0; ri < nref_i && ri < MAX_REF; ++ri) {
                size_t base_ri = base_ij + static_cast<size_t>(ri) * MAX_REF;
                for (int rj = 0; rj < nref_j && rj < MAX_REF; ++rj) {
                    double c6ref = m_c6_flat_cache[base_ri + rj];
                    if (std::abs(c6ref) < 1e-20) continue;

                    dc6_di += dgw_i[ri] * gw_j[rj] * c6ref;
                    dc6_dj += gw_i[ri] * dgw_j[rj] * c6ref;
                }
            }

            if (i == j) {
                m_dc6dcn(i, i) = dc6_di + dc6_dj;
            } else {
                m_dc6dcn(i, j) = dc6_di;  // dC6(i,j)/dCN(i)
                m_dc6dcn(j, i) = dc6_dj;  // dC6(i,j)/dCN(j) — safe, unique writer
            }
        }
    }
    };  // end dc6_worker lambda

    if (num_threads > 1 && natoms > 64) {
        int T = std::min(num_threads, natoms);
        std::vector<std::thread> threads(T - 1);
        for (int t = 1; t < T; ++t)
            threads[t - 1] = std::thread(dc6_worker, t, T);
        dc6_worker(0, T);
        for (auto& th : threads) th.join();
    } else {
        dc6_worker(0, 1);
    }

    m_dc6dcn_computed = true;
}

// Claude Generated (Feb 15, 2026): Update CN values and recompute weight derivatives + dc6dcn
// Called from GFNFF::Calculation() when gradient is requested
// Reference: Fortran gfnff_gdisp0.f90:382-395 - dc6dcn used for dispersion gradient
void D4ParameterGenerator::updateCNValuesForGradient(const std::vector<double>& cn, int num_threads)
{
    m_cn_values = cn;

    // Recompute gaussian weights with new CN values
    precomputeGaussianWeights(num_threads);

    // Compute weight derivatives and dc6dcn matrix
    computeGaussianWeightDerivatives(num_threads);
    computeDC6DCN(num_threads);
}

// Claude Generated (March 2026): Native dispersion pair generation — bypasses JSON entirely
// For 1410 atoms (993345 pairs), JSON overhead was ~10 seconds. Native: ~1 second.
std::vector<GFNFFDispersion> D4ParameterGenerator::GenerateDispersionPairsNative(
    const std::vector<int>& atoms, const Matrix& geometry_bohr)
{
    m_atoms = atoms;

    // Step 1: CN calculation
    m_cn_values = CNCalculator::calculateGFNFFCN(m_atoms, geometry_bohr);

    // Step 2: Reuse topology charges or compute fresh
    if (m_topology_charges.size() != static_cast<Eigen::Index>(m_atoms.size())) {
        Vector cn_eigen = Eigen::Map<const Vector>(m_cn_values.data(), m_cn_values.size());
        m_eeq_charges = m_eeq_solver->calculateCharges(m_atoms, geometry_bohr, 0, &cn_eigen);
    } else {
        m_eeq_charges = m_topology_charges;
    }

    // Step 3: Pre-compute Gaussian weights and C6 reference matrix
    if (!m_data_initialized) initializeReferenceData();
    if (!m_c6_reference_cached) precomputeC6ReferenceMatrix();
    precomputeGaussianWeights();

    // Step 4: Generate pairs as native structs (no JSON!)
    double a1 = m_config.get<double>("d4_a1", 0.58);
    double a2 = m_config.get<double>("d4_a2", 4.80);
    double s6 = m_config.get<double>("d4_s6", 1.0);
    double s8 = m_config.get<double>("d4_s8", 2.0);

    std::vector<GFNFFDispersion> all_pairs;

    #pragma omp parallel
    {
        std::vector<GFNFFDispersion> local_pairs;

        #pragma omp for schedule(dynamic, 10)
        for (size_t i = 0; i < m_atoms.size(); ++i) {
            for (size_t j = i + 1; j < m_atoms.size(); ++j) {
                int atom_i = m_atoms[i];
                int atom_j = m_atoms[j];

                if (atom_i < 1 || atom_i > MAX_ELEM || atom_j < 1 || atom_j > MAX_ELEM)
                    continue;

                double c6 = getChargeWeightedC6(atom_i, atom_j, i, j);
                double sqrt_zr4r2_i = getSqrtZr4r2(atom_i);
                double sqrt_zr4r2_j = getSqrtZr4r2(atom_j);
                double r4r2ij = 3.0 * sqrt_zr4r2_i * sqrt_zr4r2_j;
                double r0_sq = std::pow(a1 * std::sqrt(r4r2ij) + a2, 2);

                double q_i, q_j;
                if (m_topology_charges.size() > 0) {
                    q_i = (i < static_cast<size_t>(m_topology_charges.size())) ? m_topology_charges(i) : 0.0;
                    q_j = (j < static_cast<size_t>(m_topology_charges.size())) ? m_topology_charges(j) : 0.0;
                } else {
                    q_i = (i < static_cast<size_t>(m_eeq_charges.size())) ? m_eeq_charges(i) : 0.0;
                    q_j = (j < static_cast<size_t>(m_eeq_charges.size())) ? m_eeq_charges(j) : 0.0;
                }
                double zetac6 = GFNFFParameters::zetaChargeScale(atom_i, q_i)
                              * GFNFFParameters::zetaChargeScale(atom_j, q_j);

                GFNFFDispersion d;
                d.i = static_cast<int>(i);
                d.j = static_cast<int>(j);
                d.C6 = c6;
                d.r4r2ij = r4r2ij;
                d.r0_squared = r0_sq;
                d.r_cut = 50.0;
                d.zetac6 = zetac6;
                d.C8 = 3.0 * c6 * sqrt_zr4r2_i * sqrt_zr4r2_j;
                d.s6 = s6;
                d.s8 = s8;
                d.a1 = a1;
                d.a2 = a2;

                local_pairs.push_back(d);
            }
        }

        #pragma omp critical
        {
            all_pairs.insert(all_pairs.end(), local_pairs.begin(), local_pairs.end());
        }
    }

    // Also compute dc6dcn for gradient computation
    computeGaussianWeightDerivatives();
    computeDC6DCN();

    if (CurcumaLogger::get_verbosity() >= 1) {
        CurcumaLogger::result_fmt("D4 native dispersion: {} pairs generated", all_pairs.size());
    }

    return all_pairs;
}
