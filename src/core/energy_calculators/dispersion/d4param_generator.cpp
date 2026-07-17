/*
 * DFT-D4 Parameter Generator for Curcuma - Complete Data Integration
 * Integrates 118 elements of D4 reference data from GFN-FF Fortran implementation
 *
 * Claude Generated (December 2025): Phase 2.1 - D4 Reference Data Integration
 * December 25, 2025: Complete alphaiw (polarizability) data integration
 */

#include "d4param_generator.h"
#include "d4_charge_scaling.h"  // shared exact dftd4 zeta (AP6b)
#include "d4_ncoord.h"          // dftd4 EN-weighted covalent CN (GFN2 D4 path)
#include "src/core/energy_calculators/ff_methods/gfnff_par.h"
#include "d4_reference_data_fixed.cpp"  // D4 reference data (365 lines)
#include "d4_reference_cn_fortran.cpp"  // D4 reference CN data (Fortran dftd3param.f90 - January 2026)
#include "d4_alphaiw_data.cpp"  // D4 alphaiw polarizability data
#include "d4_corrections_data.cpp"  // D4 correction factors
#include "src/core/curcuma_logger.h"

#include <algorithm>
#include <cmath>
#include <chrono>
#include <fstream>  // Claude Generated (Mar 5, 2026): For diagnostic JSON output
#include <limits>
#include <thread>
#include <future>
#include <mutex>  // Claude Generated (Jun 2026): std::call_once guard for shared global D4 tables
#include <omp.h>  // Claude Generated (February 2026): Phase 3 - OpenMP parallelization
#include "external/CxxThreadPool/include/CxxThreadPool.hpp"

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

    // Lever 3 Opt B (Jun 2026): per-atom half-contraction fast path for C6/dc6dcn.
    // Plumbed down from the gfnff PARAM disp_half_contraction (default true).
    m_use_half_contraction = m_config.get<bool>("d4_disp_half_contraction", true);
    m_elem_slot.fill(-1);

    initializeReferenceData();

    // Initialize the shared global polarizability/correction tables exactly once for the
    // whole process. initialize_d4_alphaiw()/initialize_d4_corrections() populate MUTABLE
    // process-global containers (d4_alphaiw_data, d4_ascale_data, ...) by .resize()+assign.
    // Constructing several D4ParameterGenerators concurrently (e.g. parallel ConfSearch
    // OptThreads, each spinning up GFN-FF) otherwise races a resize/realloc against a read
    // in computeC6Reference -> SIGSEGV. The data is constant, so once-init is equivalent to
    // the previous per-construction re-init (and cheaper). call_once blocks late constructors
    // until the first finishes populating. (Claude Generated, Jun 2026)
    static std::once_flag d4_global_tables_once;
    std::call_once(d4_global_tables_once, []() {
        initialize_d4_alphaiw();      // 269 reference states, frequency-dependent polarizabilities
        initialize_d4_corrections();  // ascale, sscale, secaiw, refsys correction factors
    });

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
    m_refh = ::d4_refh;  // Reference hydrogen counts (118 × 7) — dftd4 calls this hcount
    m_refh_charges = ::d4_refh_charges;  // Reference H-charges (dftd4 refh) — used by GFN2 α-correction zeta scaling
    m_refcovcn = ::d4_refcovcn;  // dftd4 refcovcn — covalent reference CN for the CN-Gaussian weighting (set_refcn)

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

    // D4 Moment Ratios (r4Overr2), full 118-element table from
    // external/gfnff/src/dftd4param.f90:134-157 (identical to dftd4 / tblite).
    // These set the C8/C6 ratio and the BJ damping radius R0 for every element.
    // NOTE (2026-07): previously only H-Kr (Z<=36) were provided and everything
    // heavier was silently filled with a placeholder 10.0, which made native D4
    // dispersion wrong for I and all 4d/5d elements (too-large R0 -> underbinding).
    // The full array below fixes heavy elements without touching Z<=36.
    m_r4_over_r2 = {
        8.0589, 3.4698,  // H-He
        29.0974, 14.8517, 11.8799, 7.8715, 5.5588, 4.7566, 3.8025, 3.1036,  // Li-Ne
        26.1552, 17.2304, 17.7210, 12.7442, 9.5361, 8.1652, 6.7463, 5.6004,  // Na-Ar
        29.2012, 22.3934, // K-Ca
        19.0598, 16.8590, 15.4023, 12.5589, 13.4788, // Sc-Mn
        12.2309, 11.2809, 10.5569, 10.1428, 9.4907,  // Fe-Zn
        13.4606, 10.8544, 8.9386, 8.1350, 7.1251, 6.1971, // Ga-Kr
        30.0162, 24.4103, // Rb-Sr
        20.3537, 17.4780, 13.5528, 11.8451, 11.0355, // Y-Tc
        10.1997, 9.5414, 9.0061, 8.6417, 8.9975,     // Ru-Cd
        14.0834, 11.8333, 10.0179, 9.3844, 8.4110, 7.5152, // In-Xe
        32.7622, 27.5708, // Cs-Ba
        23.1671, 21.6003, 20.9615, 20.4562, 20.1010, 19.7475, 19.4828, // La-Eu
        15.6013, 19.2362, 17.4717, 17.8321, 17.4237, 17.1954, 17.1631, // Gd-Yb
        14.5716, 15.8758, 13.8989, 12.4834, 11.4421, // Lu-Re
        10.2671, 8.3549, 7.8496, 7.3278, 7.4820,     // Os-Hg
        13.5124, 11.6554, 10.0959, 9.7340, 8.8584, 8.0125, // Tl-Rn
        29.8135, 26.3157, // Fr-Ra
        19.1885, 15.8542, 16.1305, 15.6161, 15.1226, 16.1576, 0.0000, // Ac-Am
        0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, // Cm-No
        0.0000, 0.0000, 0.0000, 0.0000, 0.0000, // Lr-Bh
        0.0000, 0.0000, 0.0000, 0.0000, 5.4929, // Hs-Cn
        6.7286, 6.5144, 10.9169, 10.3600, 9.4723, 8.6641 // Nh-Og
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
    if (m_use_d4_covalent_cn) {
        // GFN2: dftd4 EN-weighted covalent CN (no log-cap) for the C6 interpolation.
        m_cn_values = curcuma::dispersion::computeD4CovalentCN(m_atoms, geometry_bohr);
    } else {
        m_cn_values = CNCalculator::calculateGFNFFCN(m_atoms, geometry_bohr);
    }
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
    } else if (m_use_d4_single_shot_eeq) {
        // Claude Generated (2026): canonical single-shot dftd4 EEQ (q-response path).
        // One smooth linear system → analytical dq/dx (D4ChargeModel). Used by
        // native GFN2 when d4_charge_source="eeq". Charges differ from the
        // GFN-FF two-phase solver, so D4 energies were re-validated for this mode.
        m_eeq_charges = m_d4_charge_model.computeCharges(m_atoms, geometry_bohr, 0.0);
        auto t_eeq_end = std::chrono::high_resolution_clock::now();
        t_eeq_ms = std::chrono::duration<double, std::milli>(t_eeq_end - t_eeq_start).count();

        if (CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::success(fmt::format("D4: single-shot EEQ charges in {:.2f} ms", t_eeq_ms));
        }
    } else {
        // Fallback: compute fresh EEQ charges via the GFN-FF two-phase solver.
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

    // Claude Generated (2026-05-28, GFN2-D4 C-path fix): rebuild the C6 reference
    // matrix if it was invalidated by a setD4CovalentCN()/setUseD4SingleShotEEQ()
    // toggle after construction. The constructor builds the cache with
    // m_use_d4_covalent_cn=false (unscaled); native GFN2 sets the flag true
    // afterwards to enable the set_refalpha_gfn2 alpha-zeta correction in
    // computeC6Reference. Without this rebuild, GenerateParameters/weightedC6Gfn2
    // read the stale UNSCALED reference C6 — bit-identical on hcount=0 references
    // but up to ~57% off on the zeta-corrected intermediate references, which is
    // the source of the CH4/triose C-path residual (see docs/GFN2_D4_STATUS.md,
    // diag_curcuma_d4_c6). GFN-FF (flag off) is unaffected: computeC6Reference's
    // zeta correction is gated by m_use_d4_covalent_cn.
    if (!m_c6_reference_cached) precomputeC6ReferenceMatrix();

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

    json dispersion_pairs = json::array();

    // Native GFN2 disables this whole block (setBuildPairLists(false)): it reads
    // the C6 reference cache through D4Evaluator and never consumes the JSON pair
    // list, which is otherwise rebuilt — a nlohmann::json object per atom pair
    // (12 string-keyed inserts each) plus an OpenMP spawn — on every geometry.
    // GFN-FF keeps it on; its ForceFieldThread consumes d4_dispersion_pairs.
    if (m_build_pair_lists) {
    std::vector<json> pairs_vec;

    #pragma omp parallel
    {
        std::vector<json> local_pairs;

        // NOTE: Cannot use collapse(2) with triangular loops (j = i + 1)
        // Parallelize outer loop only, still provides good speedup
        // Claude Generated 2026 - MSVC's OpenMP (2.0 semantics) requires a signed
        // loop index; GCC/Clang accept size_t fine, so this was never caught before.
        #pragma omp for schedule(dynamic, 10)
        for (int i = 0; i < static_cast<int>(m_atoms.size()); ++i) {
            for (int j = i + 1; j < static_cast<int>(m_atoms.size()); ++j) {
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
    for (const auto& p : pairs_vec) {
        dispersion_pairs.push_back(p);
        num_pairs++;
    }
    }  // end if (m_build_pair_lists)

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

            // Claude Generated 2026 - MSVC's OpenMP (2.0 semantics) requires a signed
            // loop index; GCC/Clang accept size_t fine, so this was never caught before.
            #pragma omp for schedule(dynamic, 10)
            for (int idx = 0; idx < static_cast<int>(triplet_vec.size()); ++idx) {
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

// Claude Generated (2026): D4 q-response chain rule (AP ∂q/∂x)
// Thin forwards to the closed-form zeta scaling and its derivative so the
// D4Evaluator can assemble dE_D4/dq without depending on gfnff_par.h directly.
double D4ParameterGenerator::getZeta(int Z, double q) const
{
    return GFNFFParameters::zetaChargeScale(Z, q);
}

double D4ParameterGenerator::getZetaDerivative(int Z, double q) const
{
    return GFNFFParameters::zetaChargeScaleDerivative(Z, q);
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

    // dftd4 GFN2 α-correction zeta scaling (set_refalpha_gfn2 in reference.f90:343-377):
    //   aiw[freq] = sscale·secaiw[freq] · zeta(ga, hardness(refsys)·gc, zeff(refsys), refh+zeff(refsys))
    // The zeta factor is frequency-independent → precompute once per (elem, ref).
    // Gated by m_use_d4_covalent_cn: GFN2 turns it on (dftd4 set_refalpha_gfn2 path);
    // GFN-FF keeps it off (its existing C6 reference is tuned against the Fortran
    // GFN-FF reference, which uses a different/no scaling — touching this regresses
    // the GFN-FF validation suite).
    auto zetaCorrection = [&](int elem, int ref, int refsys) -> double {
        if (!m_use_d4_covalent_cn) return 1.0;
        if (refsys <= 0 || refsys > MAX_ELEM) return 1.0;
        const int idx = refsys - 1;
        const double hardness = GFNFFParameters::zeta_c[idx];        // dftd4 chemical_hardness(refsys)
        const double zeff_rs  = GFNFFParameters::zeta_zeff[idx];     // dftd4 effective_nuclear_charge(refsys)
        const double refh_q   = (elem < static_cast<int>(m_refh_charges.size())
                                 && ref < static_cast<int>(m_refh_charges[elem].size()))
                                ? m_refh_charges[elem][ref] : 0.0;
        return curcuma::dispersion::d4_zeta(3.0, hardness * 2.0, zeff_rs, refh_q + zeff_rs);
    };
    const double zeta_corr_i = zetaCorrection(elem_i, ref_i, refsys_i);
    const double zeta_corr_j = zetaCorrection(elem_j, ref_j, refsys_j);

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

        // Apply correction formula (dftd4 set_refalpha_gfn2): the inner-shell
        // aiw is scaled by the precomputed zeta factor (zeta_corr_*).
        double alpha_i_iw = ascale_i * (alphaiw_i_iw - hcount_i * sscale_i * secaiw_i_iw * zeta_corr_i);
        double alpha_i_next = ascale_i * (alphaiw_i_next - hcount_i * sscale_i * secaiw_i_next * zeta_corr_i);
        double alpha_j_iw = ascale_j * (alphaiw_j_iw - hcount_j * sscale_j * secaiw_j_iw * zeta_corr_j);
        double alpha_j_next = ascale_j * (alphaiw_j_next - hcount_j * sscale_j * secaiw_j_next * zeta_corr_j);

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

void D4ParameterGenerator::precomputeGaussianWeights(CxxThreadPool* pool, int num_threads)
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

        // Claude Generated (WP5, May 2026): Stack-Eigen::Array pipeline for the inner
        // MAX_REF=7 loop. Replaces per-atom std::vector<double>(nref) heap allocation +
        // scalar std::exp loop. w_arr.exp() lets Eigen vectorise where possible while
        // remaining bit-identical to element-wise std::exp (Eigen's default exp() on a
        // scalar Array element delegates to libm, no approximation).
        const int n = std::min(nref, MAX_REF);
        const auto& refcn_row = m_refcn[elem_i];
        Eigen::Array<double, MAX_REF, 1> w_arr = Eigen::Array<double, MAX_REF, 1>::Zero();
        for (int ref = 0; ref < n; ++ref) {
            const double diff = cni - refcn_row[ref];
            w_arr(ref) = -wf * diff * diff;
        }
        w_arr.head(n) = w_arr.head(n).exp();
        const double sum_weights = w_arr.head(n).sum();

        // Normalize by the finite sum — faithful port of s-dftd4 weight_references:
        // 1/sum is finite for ANY nonzero sum (e.g. 1/2.3e-18 for an atom whose
        // reference CNs are far from the actual CN), so always divide, then guard
        // EACH normalized weight against NaN/Inf (true underflow: sum rounded to 0),
        // in which case the reference with the maximal CN gets weight 1. The old
        // `sum_weights > 1e-10` guard wrongly took the fallback for atoms with sparse
        // reference CNs, collapsing C6 onto one reference and biasing dispersion
        // (e.g. GFN2 Ti in MOR41 PR40 over-bound by ~2e-3 Eh).
        double max_cn_ref = refcn_row[0];
        for (int ref = 1; ref < n; ++ref) max_cn_ref = std::max(max_cn_ref, refcn_row[ref]);
        std::vector<double> weights(n, 0.0);
        const double inv = 1.0 / sum_weights;
        for (int ref = 0; ref < n; ++ref) {
            double gwk = w_arr(ref) * inv;
            if (!std::isfinite(gwk)) gwk = (refcn_row[ref] == max_cn_ref) ? 1.0 : 0.0;
            weights[ref] = gwk;
        }

        m_gaussian_weights[i] = std::move(weights);
    }
    };  // end gw_worker lambda

    // Dispatch: parallel if beneficial, otherwise single-threaded
    if (num_threads > 1 && natoms_gw > 64) {
        int T = std::min(num_threads, natoms_gw);
        if (pool) {
            std::vector<std::future<void>> futures;
            futures.reserve(T - 1);
            for (int t = 1; t < T; ++t)
                futures.push_back(pool->enqueue(gw_worker, t, T));
            gw_worker(0, T);
            for (auto& f : futures) f.get();
        } else {
            std::vector<std::thread> threads(T - 1);
            for (int t = 1; t < T; ++t)
                threads[t - 1] = std::thread(gw_worker, t, T);
            gw_worker(0, T);
            for (auto& th : threads) th.join();
        }
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

// Claude Generated (Lever 3 Opt B, Jun 2026): collect the distinct elements present
// in m_atoms and map Z -> compact slot. Used to size the per-atom × partner-element
// half-contraction tables. K = number of distinct elements (small in practice).
void D4ParameterGenerator::buildPresentElementMap()
{
    m_elem_slot.fill(-1);
    m_present_elems.clear();
    for (int Z : m_atoms) {
        if (Z < 1 || Z > MAX_ELEM) continue;
        if (m_elem_slot[Z] < 0) {
            m_elem_slot[Z] = static_cast<int>(m_present_elems.size());
            m_present_elems.push_back(Z);
        }
    }
    m_c6_half_K = static_cast<int>(m_present_elems.size());
}

// Claude Generated (Lever 3 Opt B, Jun 2026): per-atom half-contraction of the C6
// reference block over the FIRST ref index, for every partner element present.
//   g_i[e][b]  = Σ_a w_i[a]·c6ref[Zi][e][a][b]      (from m_gaussian_weights)
//   dg_i[e][b] = Σ_a dw_i[a]·c6ref[Zi][e][a][b]     (from m_gaussian_weight_derivatives)
// Requires both weight arrays. O(N·K·49) — ~1% of the O(N²·49) it removes from the
// energy loop and computeDC6DCN. Disabled (m_c6_half_valid=false) when the flag is off.
void D4ParameterGenerator::precomputeC6HalfContraction(CxxThreadPool* pool, int num_threads)
{
    m_c6_half_valid = false;
    if (!m_use_half_contraction) return;

    buildPresentElementMap();
    const int N = static_cast<int>(m_atoms.size());
    const int K = m_c6_half_K;
    if (N == 0 || K == 0) return;

    m_c6_half_g.assign(static_cast<size_t>(N) * K * MAX_REF, 0.0);
    m_c6_half_dg.assign(static_cast<size_t>(N) * K * MAX_REF, 0.0);
    m_c6_half_natoms = N;

    auto worker = [&](int t_id, int T) {
        for (int i = t_id; i < N; i += T) {
            int Zi = m_atoms[i];
            int elem_i = Zi - 1;
            if (elem_i < 0 || elem_i >= MAX_ELEM) continue;
            if (i >= static_cast<int>(m_gaussian_weights.size())) continue;
            const auto& wi  = m_gaussian_weights[i];
            const auto& dwi = m_gaussian_weight_derivatives[i];
            if (wi.empty() || dwi.empty()) continue;
            const int nref_i = std::min<int>(wi.size(), MAX_REF);

            for (int s = 0; s < K; ++s) {
                int elem_e = m_present_elems[s] - 1;
                const size_t base = static_cast<size_t>(elem_i) * MAX_ELEM * MAX_REF * MAX_REF
                                  + static_cast<size_t>(elem_e) * MAX_REF * MAX_REF;
                const int nref_e = std::min<int>(
                    (elem_e < static_cast<int>(m_refn.size())) ? m_refn[elem_e] : 0, MAX_REF);
                double* gout  = &m_c6_half_g [c6HalfIndex(i, s)];
                double* dgout = &m_c6_half_dg[c6HalfIndex(i, s)];
                for (int b = 0; b < nref_e; ++b) {
                    double gv = 0.0, dgv = 0.0;
                    for (int a = 0; a < nref_i; ++a) {
                        double c6ref = m_c6_flat_cache[base + static_cast<size_t>(a) * MAX_REF + b];
                        gv  += wi[a]  * c6ref;
                        dgv += dwi[a] * c6ref;
                    }
                    gout[b] = gv;
                    dgout[b] = dgv;
                }
            }
        }
    };

    if (num_threads > 1 && N > 64) {
        int T = std::min(num_threads, N);
        if (pool) {
            std::vector<std::future<void>> futures;
            futures.reserve(T - 1);
            for (int t = 1; t < T; ++t)
                futures.push_back(pool->enqueue(worker, t, T));
            worker(0, T);
            for (auto& f : futures) f.get();
        } else {
            std::vector<std::thread> threads(T - 1);
            for (int t = 1; t < T; ++t)
                threads[t - 1] = std::thread(worker, t, T);
            worker(0, T);
            for (auto& th : threads) th.join();
        }
    } else {
        worker(0, 1);
    }

    m_c6_half_valid = true;
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

    // Lever 3 Opt B fast path: c6(i,j) = Σ_b g_i[Zj][b]·w_j[b] using the per-atom
    // half-contraction (7 FMA instead of up to 49). Reassociates the FP sum vs the
    // flat loop below (~1e-16). Skipped at verbosity>=3 so the per-term debug print
    // in the flat path still works. Falls through to the exact path otherwise.
    if (m_c6_half_valid && CurcumaLogger::get_verbosity() < 3
        && static_cast<int>(atom_i) < m_c6_half_natoms) {
        const int slot_j = m_elem_slot[Zj];
        if (slot_j >= 0 && atom_j < m_gaussian_weights.size()) {
            const auto& wj = m_gaussian_weights[atom_j];
            const double* g = &m_c6_half_g[c6HalfIndex(static_cast<int>(atom_i), slot_j)];
            const int nb = std::min<int>(wj.size(), MAX_REF);
            double c6 = 0.0;
            for (int b = 0; b < nb; ++b) c6 += g[b] * wj[b];
            return c6;
        }
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

// Claude Generated (AP2 perf, 2026-06): per-atom reference weights, lifted out of
// weightedC6Gfn2 so the O(N) build runs once per atom instead of ~N times inside the
// O(N²) pair loop. Depends only on the atom's element/charge/CN. Mirrors dftd4
// weight_references exactly (numbers bit-identical to the old inline lambda).
D4ParameterGenerator::RefW
D4ParameterGenerator::buildAtomRefW(int Z, size_t atom_idx, double q,
                                    bool want_grad, bool want_hess) const
{
    using curcuma::dispersion::d4_zeta;
    using curcuma::dispersion::d4_dzeta;
    using curcuma::dispersion::d4_d2zeta;

    RefW r;
    const int elem = Z - 1;
    if (elem < 0 || elem >= MAX_ELEM) return r;
    if (atom_idx >= m_cn_values.size()) return r;
    const double cn = m_cn_values[atom_idx];

    constexpr double ga = 3.0;     // dftd4 ga_default (charge-scaling height)
    constexpr double gc = 2.0;     // dftd4 gc_default (charge-scaling steepness)
    constexpr double wf = 6.0;     // dftd4 wf_default (CN gaussian width)
    constexpr int    MAXCN = 19;   // dftd4 set_refgw max_cn

    const int nref = (elem < (int)m_refn.size()) ? m_refn[elem] : 0;
    r.nref = nref;
    if (nref <= 0) return r;
    const double eta  = GFNFFParameters::zeta_c[elem];      // = dftd4 chemical_hardness
    const double zeff = GFNFFParameters::zeta_zeff[elem];   // = dftd4 effective_nuclear_charge
    const double gi   = eta * gc;
    const double* refcn = m_refcn[elem].data();
    // dftd4 uses TWO reference-CN tables: refcn for the ngw bucketing
    // (set_refgw, below) and refcovcn (covalent CN) for the actual
    // CN-Gaussian (set_refcn -> model%cn -> weight_references). Using refcn
    // for the Gaussian over-broadens the weights (up to 14% wrong weighted
    // C6 on carbon -> the triose C-path residual). Fall back to refcn only
    // if the refcovcn table is unavailable. (Claude Generated 2026-05-28.)
    const double* refcovcn = (elem < (int)m_refcovcn.size() && !m_refcovcn[elem].empty())
                             ? m_refcovcn[elem].data() : refcn;

    // ngw[ref] from refcn (dftd4 set_refgw): count references sharing the
    // same rounded CN, then ngw = k(k+1)/2.
    int cnc[MAXCN + 1];
    for (int k = 0; k <= MAXCN; ++k) cnc[k] = 0;
    cnc[0] = 1;
    for (int ir = 0; ir < nref; ++ir) {
        int icn = (int)std::lround(refcn[ir]); if (icn < 0) icn = 0; if (icn > MAXCN) icn = MAXCN;
        cnc[icn] += 1;
    }
    int ngw[MAX_REF];
    for (int ir = 0; ir < nref; ++ir) {
        int icn = (int)std::lround(refcn[ir]); if (icn < 0) icn = 0; if (icn > MAXCN) icn = MAXCN;
        int k = cnc[icn];
        ngw[ir] = k * (k + 1) / 2;
    }

    // CN gaussian weights gwk[ref] (+ dgwk/dCN). dftd4 weight_cn = exp(-wf·(cn-cnref)²),
    // summed over igw=1..ngw with wf_eff = igw·wf.
    double expw[MAX_REF] = {0}, expd[MAX_REF] = {0};
    double norm = 0.0, dnorm = 0.0;
    for (int ir = 0; ir < nref; ++ir) {
        double ew = 0.0, ed = 0.0;
        for (int igw = 1; igw <= ngw[ir]; ++igw) {
            const double wfe = igw * wf;
            const double d   = cn - refcovcn[ir];        // Gaussian uses refcovcn (dftd4 model%cn)
            const double gw  = std::exp(-wfe * d * d);
            ew += gw;
            ed += 2.0 * wfe * (refcovcn[ir] - cn) * gw;  // d(gw)/dCN
        }
        expw[ir] = ew; expd[ir] = ed;
        norm += ew; dnorm += ed;
    }
    const double ninv = (norm > 0.0) ? 1.0 / norm : 0.0;

    for (int ir = 0; ir < nref; ++ir) {
        double gwk  = expw[ir] * ninv;
        double dgwk = ninv * (expd[ir] - expw[ir] * dnorm * ninv);
        if (!std::isfinite(gwk))  gwk  = 0.0;
        if (!std::isfinite(dgwk)) dgwk = 0.0;
        const double qref = m_refq[elem][ir] + zeff;
        const double qmod = q + zeff;
        const double z  = d4_zeta(ga, gi, qref, qmod);
        r.W[ir]   = gwk * z;
        if (want_grad) {
            r.dWq[ir] = gwk  * d4_dzeta(ga, gi, qref, qmod);
            r.dWc[ir] = dgwk * z;
        }
        if (want_hess) r.d2Wq[ir] = gwk * d4_d2zeta(ga, gi, qref, qmod);
    }
    return r;
}

// Claude Generated (AP2 perf, 2026-06): the bit-identical inner 7×7 contraction of
// weightedC6Gfn2 — ΣΣ ri.W[a]·rj.W[b]·c6ref over the cached element-pair block.
D4ParameterGenerator::C6Gfn2
D4ParameterGenerator::contractC6Gfn2(const RefW& ri, const RefW& rj, int Zi, int Zj,
                                     bool want_grad, bool want_hess) const
{
    C6Gfn2 out;
    const int ei = Zi - 1, ej = Zj - 1;
    if (ei < 0 || ei >= MAX_ELEM || ej < 0 || ej >= MAX_ELEM) return out;

    const size_t base = static_cast<size_t>(ei) * MAX_ELEM * MAX_REF * MAX_REF
                      + static_cast<size_t>(ej) * MAX_REF * MAX_REF;
    for (int a = 0; a < ri.nref; ++a) {
        const size_t base_a = base + static_cast<size_t>(a) * MAX_REF;
        for (int b = 0; b < rj.nref; ++b) {
            const double c6ref = m_c6_flat_cache[base_a + b];
            out.c6 += ri.W[a] * rj.W[b] * c6ref;
            if (want_grad) {
                out.dc6dqi  += ri.dWq[a] * rj.W[b]   * c6ref;
                out.dc6dqj  += ri.W[a]   * rj.dWq[b] * c6ref;
                out.dc6dcni += ri.dWc[a] * rj.W[b]   * c6ref;
                out.dc6dcnj += ri.W[a]   * rj.dWc[b] * c6ref;
            }
            if (want_hess) {
                out.d2c6dqi2 += ri.d2Wq[a] * rj.W[b]   * c6ref;
                out.d2c6dqj2 += ri.W[a]    * rj.d2Wq[b]* c6ref;
            }
        }
    }
    return out;
}

// Claude Generated (Stage 5 Part B2, 2026-06): flatten the per-atom reference
// weights (W = gwk·ζ, dWq = ∂W/∂q) the GPU D4-potential pair loop consumes. The
// CN-Gaussian + zeta math stays here on the validated CPU (buildAtomRefW); the
// device only does the O(N²) 7×7 contraction × BJ disp_sum.
void D4ParameterGenerator::buildRefWFlat(const std::vector<int>& atoms, const Vector& q,
                                         std::vector<double>& W_out,
                                         std::vector<double>& dWq_out) const
{
    const int nat = static_cast<int>(atoms.size());
    W_out.assign(static_cast<size_t>(nat) * MAX_REF, 0.0);
    dWq_out.assign(static_cast<size_t>(nat) * MAX_REF, 0.0);
    for (int a = 0; a < nat; ++a) {
        const double qa = (a < q.size()) ? q(a) : 0.0;
        const RefW r = buildAtomRefW(atoms[a], static_cast<size_t>(a), qa, true, false);
        const int nr = (r.nref < MAX_REF) ? r.nref : MAX_REF;
        for (int ref = 0; ref < nr; ++ref) {
            W_out[static_cast<size_t>(a) * MAX_REF + ref]   = r.W[ref];
            dWq_out[static_cast<size_t>(a) * MAX_REF + ref] = r.dWq[ref];
        }
    }
}

// Overload that also emits dWc = ∂W/∂CN (for the device 2-body D4 gradient kernel).
void D4ParameterGenerator::buildRefWFlat(const std::vector<int>& atoms, const Vector& q,
                                         std::vector<double>& W_out,
                                         std::vector<double>& dWq_out,
                                         std::vector<double>& dWc_out) const
{
    const int nat = static_cast<int>(atoms.size());
    W_out.assign(static_cast<size_t>(nat) * MAX_REF, 0.0);
    dWq_out.assign(static_cast<size_t>(nat) * MAX_REF, 0.0);
    dWc_out.assign(static_cast<size_t>(nat) * MAX_REF, 0.0);
    for (int a = 0; a < nat; ++a) {
        const double qa = (a < q.size()) ? q(a) : 0.0;
        const RefW r = buildAtomRefW(atoms[a], static_cast<size_t>(a), qa, true, false);
        const int nr = (r.nref < MAX_REF) ? r.nref : MAX_REF;
        for (int ref = 0; ref < nr; ++ref) {
            const size_t idx = static_cast<size_t>(a) * MAX_REF + ref;
            W_out[idx]   = r.W[ref];
            dWq_out[idx] = r.dWq[ref];
            dWc_out[idx] = r.dWc[ref];
        }
    }
}

// q=0 reference C6 (+ ∂C6/∂CN) matrices for the device ATM kernel. Mirrors the c6/dc6dcn
// build at the top of D4Evaluator::computeATM (buildAtomRefW at q=0 + contractC6Gfn2).
void D4ParameterGenerator::buildAtmC6Flat(const std::vector<int>& atoms,
                                          std::vector<double>& c6_out,
                                          std::vector<double>& dc6dcn_out) const
{
    const int nat = static_cast<int>(atoms.size());
    c6_out.assign(static_cast<size_t>(nat) * nat, 0.0);
    dc6dcn_out.assign(static_cast<size_t>(nat) * nat, 0.0);
    std::vector<RefW> refw0(nat);
    for (int a = 0; a < nat; ++a)
        refw0[a] = buildAtomRefW(atoms[a], static_cast<size_t>(a), /*q=*/0.0, true, false);
    for (int a = 0; a < nat; ++a) {
        for (int b = 0; b <= a; ++b) {
            const C6Gfn2 g = contractC6Gfn2(refw0[a], refw0[b], atoms[a], atoms[b], true, false);
            c6_out[static_cast<size_t>(a) * nat + b] = g.c6;
            c6_out[static_cast<size_t>(b) * nat + a] = g.c6;
            dc6dcn_out[static_cast<size_t>(a) * nat + b] = g.dc6dcni;   // ∂C6(a,b)/∂CN_a
            dc6dcn_out[static_cast<size_t>(b) * nat + a] = g.dc6dcnj;   // ∂C6(a,b)/∂CN_b
        }
    }
}

// Stage 6 (S6.2b, Claude Generated 2026-06): export the q-independent per-atom
// reference data for the device W/dWq rebuild (k_d4_build_refw). Mirrors the
// per-atom setup at the top of buildAtomRefW (gc=2 → gi=eta·gc; the refcovcn
// fallback to refcn). The device kernel then redoes the ngw bucketing + CN-
// Gaussian + zeta from these tables, so the resident SCF charges drive W/dWq
// without a host round-trip.
void D4ParameterGenerator::exportRefWDeviceData(const std::vector<int>& atoms,
                                                std::vector<double>& cn, std::vector<double>& gi,
                                                std::vector<double>& zeff, std::vector<double>& refcn,
                                                std::vector<double>& refcovcn, std::vector<double>& refq,
                                                std::vector<int>& nref) const
{
    constexpr double gc = 2.0;   // dftd4 gc_default (matches buildAtomRefW)
    const int nat = static_cast<int>(atoms.size());
    cn.assign(nat, 0.0);
    gi.assign(nat, 0.0);
    zeff.assign(nat, 0.0);
    nref.assign(nat, 0);
    refcn.assign(static_cast<size_t>(nat) * MAX_REF, 0.0);
    refcovcn.assign(static_cast<size_t>(nat) * MAX_REF, 0.0);
    refq.assign(static_cast<size_t>(nat) * MAX_REF, 0.0);
    for (int a = 0; a < nat; ++a) {
        const int elem = atoms[a] - 1;
        if (elem < 0 || elem >= MAX_ELEM) continue;
        cn[a]   = (static_cast<size_t>(a) < m_cn_values.size()) ? m_cn_values[a] : 0.0;
        gi[a]   = GFNFFParameters::zeta_c[elem] * gc;
        zeff[a] = GFNFFParameters::zeta_zeff[elem];
        const int nr = (elem < static_cast<int>(m_refn.size())) ? m_refn[elem] : 0;
        nref[a] = nr;
        const double* rcn = (elem < static_cast<int>(m_refcn.size())) ? m_refcn[elem].data() : nullptr;
        const double* rcov = (elem < static_cast<int>(m_refcovcn.size()) && !m_refcovcn[elem].empty())
                             ? m_refcovcn[elem].data() : rcn;
        const double* rq = (elem < static_cast<int>(m_refq.size())) ? m_refq[elem].data() : nullptr;
        for (int ir = 0; ir < nr && ir < MAX_REF; ++ir) {
            const size_t idx = static_cast<size_t>(a) * MAX_REF + ir;
            if (rcn)  refcn[idx]    = rcn[ir];
            if (rcov) refcovcn[idx] = rcov[ir];
            if (rq)   refq[idx]     = rq[ir];
        }
    }
}

// Claude Generated (AP6b exact D4 port, 2026): tblite/dftd4-exact per-reference
// charge-weighted C6 for native GFN2. Per-pair entry point — now a thin wrapper over
// buildAtomRefW + contractC6Gfn2 (numerics unchanged). Mirrors dftd4 model.f90
// weight_references (gwk with wf=6 + ngw multi-gaussian) and get_atomic_c6.
D4ParameterGenerator::C6Gfn2
D4ParameterGenerator::weightedC6Gfn2(int Zi, int Zj, size_t atom_i, size_t atom_j,
                                     double qi, double qj,
                                     bool want_grad, bool want_hess) const
{
    const RefW ri = buildAtomRefW(Zi, atom_i, qi, want_grad, want_hess);
    const RefW rj = buildAtomRefW(Zj, atom_j, qj, want_grad, want_hess);
    return contractC6Gfn2(ri, rj, Zi, Zj, want_grad, want_hess);
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
void D4ParameterGenerator::computeGaussianWeightDerivatives(CxxThreadPool* pool, int num_threads)
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

        // Claude Generated (WP5, May 2026): Stack-Eigen::Array pipeline. Same as
        // precomputeGaussianWeights but with the dgw/dCN extension (dexpw, dnorm).
        const int n = std::min(nref, MAX_REF);
        const auto& refcn_row = m_refcn[elem_i];
        Eigen::Array<double, MAX_REF, 1> diff_arr = Eigen::Array<double, MAX_REF, 1>::Zero();
        for (int ref = 0; ref < n; ++ref) diff_arr(ref) = cni - refcn_row[ref];
        Eigen::Array<double, MAX_REF, 1> expw = Eigen::Array<double, MAX_REF, 1>::Zero();
        expw.head(n) = (-wf * diff_arr.head(n).square()).exp();
        // dexpw = 2*wf * (cni_ref - cni) * expw = -2*wf * diff * expw
        Eigen::Array<double, MAX_REF, 1> dexpw = Eigen::Array<double, MAX_REF, 1>::Zero();
        dexpw.head(n) = (-2.0 * wf) * diff_arr.head(n) * expw.head(n);
        const double norm = expw.head(n).sum();
        const double dnorm = dexpw.head(n).sum();

        // Divide by the finite norm and guard each derivative against NaN/Inf (true
        // underflow), mirroring the energy-weight fix above so the CN chain-rule is
        // consistent for atoms with sparse reference CNs (the old `norm > 1e-10`
        // zeroed legitimate finite derivatives).
        std::vector<double> dgwdcn(n, 0.0);
        const double inv_norm2 = 1.0 / (norm * norm);
        for (int ref = 0; ref < n; ++ref) {
            double dgwk = (dexpw(ref) * norm - expw(ref) * dnorm) * inv_norm2;
            if (!std::isfinite(dgwk)) dgwk = 0.0;
            dgwdcn[ref] = dgwk;
        }

        m_gaussian_weight_derivatives[i] = std::move(dgwdcn);
    }
    };  // end gwd_worker lambda

    if (num_threads > 1 && natoms_gwd > 64) {
        int T = std::min(num_threads, natoms_gwd);
        if (pool) {
            std::vector<std::future<void>> futures;
            futures.reserve(T - 1);
            for (int t = 1; t < T; ++t)
                futures.push_back(pool->enqueue(gwd_worker, t, T));
            gwd_worker(0, T);
            for (auto& f : futures) f.get();
        } else {
            std::vector<std::thread> threads(T - 1);
            for (int t = 1; t < T; ++t)
                threads[t - 1] = std::thread(gwd_worker, t, T);
            gwd_worker(0, T);
            for (auto& th : threads) th.join();
        }
    } else {
        gwd_worker(0, 1);
    }
}

// Claude Generated (Feb 15, 2026): Compute dc6dcn matrix from weight derivatives and C6 references
// Reference: Fortran gfnff_gdisp0.f90:262-305 (get_atomic_c6_d4 dc6dcn part)
//
// dc6dcn(i,j) = dC6(i,j)/dCN(i) = sum_{ri,rj} dgwdcn(ri,i) * gw(rj,j) * c6ref(ri,rj,Zi,Zj)
// dc6dcn(j,i) = dC6(i,j)/dCN(j) = sum_{ri,rj} gw(ri,i) * dgwdcn(rj,j) * c6ref(ri,rj,Zi,Zj)
void D4ParameterGenerator::computeDC6DCN(CxxThreadPool* pool, int num_threads)
{
    // Claude Generated (Mar 2026): Internal std::thread parallelisation for O(N²×M²)
    int natoms = static_cast<int>(m_atoms.size());
    m_dc6dcn = Matrix::Zero(natoms, natoms);

    // Per-pair dc6/dCN kernel. For i<j writes (i,j)=dC6(i,j)/dCN(i) and
    // (j,i)=dC6(i,j)/dCN(j); for i==j (full-N² fallback only) writes the diagonal
    // sum. Each unordered pair is the sole writer of its two cells, so any
    // partition over pairs/rows is race-free without buffering.
    auto pair_kernel = [&](int i, int j) {
        int Zi = m_atoms[i];
        int elem_i = Zi - 1;
        if (elem_i < 0 || elem_i >= MAX_ELEM) return;
        const auto& gw_i = m_gaussian_weights[i];
        const auto& dgw_i = m_gaussian_weight_derivatives[i];
        if (gw_i.empty() || dgw_i.empty()) return;
        int nref_i = static_cast<int>(gw_i.size());

        int Zj = m_atoms[j];
        int elem_j = Zj - 1;
        if (elem_j < 0 || elem_j >= MAX_ELEM) return;
        const auto& gw_j = m_gaussian_weights[j];
        const auto& dgw_j = m_gaussian_weight_derivatives[j];
        if (gw_j.empty() || dgw_j.empty()) return;
        int nref_j = static_cast<int>(gw_j.size());

        double dc6_di = 0.0;
        double dc6_dj = 0.0;

        const int slot_j = m_c6_half_valid ? m_elem_slot[Zj] : -1;
        if (slot_j >= 0 && i < m_c6_half_natoms) {
            // Lever 3 Opt B half-contraction: 7 FMA instead of up to 49.
            //   dc6_di = Σ_b dg_i[Zj][b]·w_j[b];  dc6_dj = Σ_b g_i[Zj][b]·dw_j[b]
            const double* g  = &m_c6_half_g [c6HalfIndex(i, slot_j)];
            const double* dg = &m_c6_half_dg[c6HalfIndex(i, slot_j)];
            const int nb = std::min<int>(nref_j, MAX_REF);
            for (int b = 0; b < nb; ++b) {
                dc6_di += dg[b] * gw_j[b];
                dc6_dj += g[b]  * dgw_j[b];
            }
        } else {
            // Exact flat path (also the disp_half_contraction=false reproduction).
            const size_t base_ij = static_cast<size_t>(elem_i) * MAX_ELEM * MAX_REF * MAX_REF
                                 + static_cast<size_t>(elem_j) * MAX_REF * MAX_REF;
            for (int ri = 0; ri < nref_i && ri < MAX_REF; ++ri) {
                size_t base_ri = base_ij + static_cast<size_t>(ri) * MAX_REF;
                for (int rj = 0; rj < nref_j && rj < MAX_REF; ++rj) {
                    double c6ref = m_c6_flat_cache[base_ri + rj];
                    if (std::abs(c6ref) < 1e-20) continue;

                    dc6_di += dgw_i[ri] * gw_j[rj] * c6ref;
                    dc6_dj += gw_i[ri] * dgw_j[rj] * c6ref;
                }
            }
        }

        if (i == j) {
            m_dc6dcn(i, i) = dc6_di + dc6_dj;
        } else {
            m_dc6dcn(i, j) = dc6_di;  // dC6(i,j)/dCN(i)
            m_dc6dcn(j, i) = dc6_dj;  // dC6(i,j)/dCN(j) — unique writer
        }
    };

    // Lever 3 Opt A: iterate only the in-cutoff pair list — exactly the dc6dcn
    // entries the gradient consumer reads. Beyond-cutoff cells and the diagonal
    // stay 0 (never read). Bit-identical for every consumed entry. Falls back to
    // the full-N² double loop when the list is unavailable (GPU skip / first call).
    if (m_dc6dcn_pairs_valid) {
        const int P = static_cast<int>(m_dc6dcn_pairs.size());
        auto pl_worker = [&](int t_id, int T) {
            for (int p = t_id; p < P; p += T)
                pair_kernel(m_dc6dcn_pairs[p].first, m_dc6dcn_pairs[p].second);
        };
        if (num_threads > 1 && P > 64) {
            int T = std::min(num_threads, P);
            if (pool) {
                std::vector<std::future<void>> futures;
                futures.reserve(T - 1);
                for (int t = 1; t < T; ++t)
                    futures.push_back(pool->enqueue(pl_worker, t, T));
                pl_worker(0, T);
                for (auto& f : futures) f.get();
            } else {
                std::vector<std::thread> threads(T - 1);
                for (int t = 1; t < T; ++t)
                    threads[t - 1] = std::thread(pl_worker, t, T);
                pl_worker(0, T);
                for (auto& th : threads) th.join();
            }
        } else {
            pl_worker(0, 1);
        }
        m_dc6dcn_computed = true;
        return;
    }

    // Fallback: original full-N² row-partitioned loop. Each thread owns row i and
    // iterates j >= i, so it is the sole writer of cells (i,j) and (j,i).
    auto dc6_worker = [&](int t_id, int T) {
        for (int i = t_id; i < natoms; i += T) {
            int elem_i = m_atoms[i] - 1;
            if (elem_i < 0 || elem_i >= MAX_ELEM) continue;
            if (m_gaussian_weights[i].empty() || m_gaussian_weight_derivatives[i].empty()) continue;
            for (int j = i; j < natoms; ++j) pair_kernel(i, j);
        }
    };

    if (num_threads > 1 && natoms > 64) {
        int T = std::min(num_threads, natoms);
        if (pool) {
            std::vector<std::future<void>> futures;
            futures.reserve(T - 1);
            for (int t = 1; t < T; ++t)
                futures.push_back(pool->enqueue(dc6_worker, t, T));
            dc6_worker(0, T);
            for (auto& f : futures) f.get();
        } else {
            std::vector<std::thread> threads(T - 1);
            for (int t = 1; t < T; ++t)
                threads[t - 1] = std::thread(dc6_worker, t, T);
            dc6_worker(0, T);
            for (auto& th : threads) th.join();
        }
    } else {
        dc6_worker(0, 1);
    }

    m_dc6dcn_computed = true;
}

// Claude Generated (Feb 15, 2026): Update CN values and recompute weight derivatives + dc6dcn
// Called from GFNFF::Calculation() when gradient is requested
// Reference: Fortran gfnff_gdisp0.f90:382-395 - dc6dcn used for dispersion gradient
// Claude Generated (Apr 2026): P1a — Skip recomputation when CN change < threshold (MD optimization)
void D4ParameterGenerator::updateCNValuesForGradient(const std::vector<double>& cn, CxxThreadPool* pool, int num_threads, bool skip_dc6dcn)
{
    m_cn_values = cn;

    // P1a: Skip recomputation if CN change is below threshold
    if (canSkipGaussianWeightsUpdate(cn)) {
        return;  // gw, dgw, dc6dcn all unchanged — valid for this step
    }

    // Recompute gaussian weights with new CN values
    precomputeGaussianWeights(pool, num_threads);

    // Compute weight derivatives (needed for dc6dcn — CPU or GPU)
    computeGaussianWeightDerivatives(pool, num_threads);

    // Claude Generated (March 2026): Phase 2 GPU dc6dcn — skip O(N²) matrix when
    // GPU will compute dc6dcn per dispersion pair directly.
    if (!skip_dc6dcn) {
        // Lever 3 Opt B: rebuild the half-contraction (CN changed → weights changed)
        // before computeDC6DCN. Reuses m_dc6dcn_pairs from the last pair-list build.
        precomputeC6HalfContraction(pool, num_threads);
        computeDC6DCN(pool, num_threads);
    }

    // Store CN for next step's threshold check
    m_prev_cn_values = cn;
    m_cn_cached = true;
}

// Claude Generated (Apr 2026): P1a — CN-change threshold check for Gaussian weight caching
bool D4ParameterGenerator::canSkipGaussianWeightsUpdate(const std::vector<double>& new_cn) const
{
    double threshold = m_config.get<double>("d4_cn_cache_threshold", 0.01);
    if (threshold <= 0.0) return false;  // Caching disabled
    if (!m_cn_cached || m_prev_cn_values.size() != new_cn.size())
        return false;
    double max_change = 0.0;
    for (size_t i = 0; i < new_cn.size(); ++i) {
        double delta = std::abs(new_cn[i] - m_prev_cn_values[i]);
        if (delta > max_change) max_change = delta;
    }
    return max_change < threshold;
}

// Claude Generated (Apr 2026): P1a — Record CN values for next step's threshold check
void D4ParameterGenerator::recordCNValues(const std::vector<double>& cn)
{
    m_prev_cn_values = cn;
    m_cn_cached = true;
}

// Claude Generated (March 2026): Native dispersion pair generation — bypasses JSON entirely
// For 1410 atoms (993345 pairs), JSON overhead was ~10 seconds. Native: ~1 second.
std::vector<GFNFFDispersion> D4ParameterGenerator::GenerateDispersionPairsNative(
    const std::vector<int>& atoms, const Matrix& geometry_bohr)
{
    m_atoms = atoms;

    // Step 1: CN calculation
    if (m_use_d4_covalent_cn) {
        // GFN2: dftd4 EN-weighted covalent CN (no log-cap). This is the CN that
        // weightedC6Gfn2 reads, so it is the authoritative one for the GFN2 C6.
        m_cn_values = curcuma::dispersion::computeD4CovalentCN(m_atoms, geometry_bohr);
    } else {
        m_cn_values = CNCalculator::calculateGFNFFCN(m_atoms, geometry_bohr);
    }

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

    // WP-A (Jun 2026): with gpu_disp_pairs_on_device the GPU builds the pair list +
    // dc6dcn on the device, reusing only the CN + Gaussian weights computed above.
    // Skip the host O(N²) pair loop and the host dc6dcn entirely; return empty.
    if (m_skip_pair_loop) {
        m_dc6dcn_pairs_valid = false;  // GPU owns dc6dcn; keep host fallback exact (full N²)
        m_c6_half_valid = false;       // GPU path: no host half-contraction
        return {};
    }

    // Lever 3 Opt B: derivatives + half-contraction BEFORE the energy loop, because
    // both the energy loop's getChargeWeightedC6 fast path and computeDC6DCN consume
    // the per-atom half-contraction. computeGaussianWeightDerivatives moved up from
    // after the loop (it depends only on CN, always ran unconditionally).
    computeGaussianWeightDerivatives();
    precomputeC6HalfContraction();

    // Step 4: Generate pairs as native structs (no JSON!)
    double a1 = m_config.get<double>("d4_a1", 0.58);
    double a2 = m_config.get<double>("d4_a2", 4.80);
    double s6 = m_config.get<double>("d4_s6", 1.0);
    double s8 = m_config.get<double>("d4_s8", 2.0);

    std::vector<GFNFFDispersion> all_pairs;

    // Claude Generated (Apr 2026): Hard distance cutoff of 60 Bohr for dispersion pairs.
    // BJ-damped dispersion is negligible beyond this range; skipping saves N²/2 pairs.
    //
    // WP-C phase D (May 2026, REVERTED): tried aligning to D3's sqrt(dispthr=1500) ≈ 38.73 Bohr.
    // XTB-GFN-FF reference comparison on polymer.xyz showed this *worsens* the match
    // (Curcuma-vs-XTB diff: +0.92969 Eh @ cutoff=60.0  →  +0.92993 Eh @ cutoff=38.73).
    // Fortran's dispthr appears to apply only to D3, not D4 — D4's effective cutoff is larger.
    // 60 Bohr stays as the empirical Curcuma D4 default.
    constexpr double disp_cutoff_bohr = 60.0;
    constexpr double disp_cutoff_sq = disp_cutoff_bohr * disp_cutoff_bohr;

    #pragma omp parallel
    {
        std::vector<GFNFFDispersion> local_pairs;

        // Claude Generated 2026 - MSVC's OpenMP (2.0 semantics) requires a signed
        // loop index; GCC/Clang accept size_t fine, so this was never caught before.
        #pragma omp for schedule(dynamic, 10)
        for (int i = 0; i < static_cast<int>(m_atoms.size()); ++i) {
            for (int j = i + 1; j < static_cast<int>(m_atoms.size()); ++j) {
                double r2 = (geometry_bohr.row(i) - geometry_bohr.row(j)).squaredNorm();
                if (r2 > disp_cutoff_sq) continue;

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
                d.r_cut = 50.0;  // D4 pair-energy cutoff (struct default; verified by XTB comparison May 2026)
                d.zetac6 = zetac6;
                // P1c (Apr 2026): Legacy D3 fields removed from GFNFFDispersion

                local_pairs.push_back(d);
            }
        }

        #pragma omp critical
        {
            all_pairs.insert(all_pairs.end(), local_pairs.begin(), local_pairs.end());
        }
    }

    // Lever 3 Opt A (Jun 2026): capture the in-cutoff pair list. computeDC6DCN only
    // needs entries the gradient consumer reads — exactly the dispersion pairs (i<j,
    // within the 60 Bohr cutoff) collected above. updateCNValuesForGradient reuses it.
    m_dc6dcn_pairs.clear();
    m_dc6dcn_pairs.reserve(all_pairs.size());
    for (const auto& d : all_pairs) m_dc6dcn_pairs.emplace_back(d.i, d.j);
    m_dc6dcn_pairs_valid = true;

    // dc6dcn for the gradient (derivatives + half-contraction already built above).
    computeDC6DCN();

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::result_fmt("D4 native dispersion: {} pairs generated", all_pairs.size());
    }

    return all_pairs;
}
