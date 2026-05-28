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
#include "src/core/energy_calculators/ff_methods/eeq_solver.h"  // EEQ charge calculation for D4
#include "src/core/energy_calculators/ff_methods/gfnff_parameters.h"  // GFNFFDispersion struct
#include "src/core/energy_calculators/ff_methods/cn_calculator.h"  // CN calculation for D4
#include "d4_charge_model.h"  // single-shot EEQ + analytical dq/dx (q-response)

#include <Eigen/Dense>
#include <memory>

#include "json.hpp"
using json = nlohmann::json;

class CxxThreadPool;  // Forward declaration for pool-based parallelisation

BEGIN_PARAMETER_DEFINITION(d4param)
    // D4 reference selection and scaling
    PARAM(d4_refq, Int, 2, "D4 reference charges (0=gfn2xtb, 1=gasteiger, 2=hirshfeld).", "Reference", {})
    PARAM(d4_s6, Double, 1.0, "D4 global scaling factor for C6 term.", "Scaling", {})
    PARAM(d4_s8, Double, 2.0, "D4 global scaling factor for C8 term.", "Scaling", {})
    PARAM(d4_s10, Double, 1.0, "D4 global scaling factor for C10 term.", "Scaling", {})
    PARAM(d4_s12, Double, 1.0, "D4 global scaling factor for C12 term.", "Scaling", {})
    PARAM(d4_s9, Double, 1.0, "D4 scaling for three-body ATM term (default: enabled).", "Scaling", {})
    PARAM(d4_a1, Double, 0.58, "D4 damping parameter a1 (GFN-FF: 0.58, GFN2-xTB: 0.63).", "Damping", {})
    PARAM(d4_a2, Double, 4.80, "D4 damping parameter a2 (Bohr) - GFN-FF: 4.80, GFN2-xTB: 5.0.", "Damping", {})
    PARAM(d4_alp, Double, 14.0, "D4 alpha damping parameter.", "Damping", {})

    // D4 atomic polarizability data
    PARAM(d4_r4r2_model, Int, 1, "D4 r4/r2 model (0=static, 1=dynamic).", "Model", {})

    // Claude Generated (Apr 2026): P1a — CN-change threshold for Gaussian weight caching
    PARAM(d4_cn_cache_threshold, Double, 0.01, "CN change threshold for skipping Gaussian weight recomputation (MD optimization). Set to 0.0 to disable caching.", "Advanced", {})
END_PARAMETER_DEFINITION

class D4ParameterGenerator {
public:
    D4ParameterGenerator(const ConfigManager& config);
    ~D4ParameterGenerator() = default;

    // Main parameter generation interface
    void GenerateParameters(const std::vector<int>& atoms, const Matrix& geometry_bohr);
    json getParameters() const { return m_parameters; }

    /**
     * @brief Generate dispersion pairs as native structs (bypasses JSON entirely)
     *
     * Claude Generated (March 2026): Performance optimization for large systems.
     * For 1410 atoms (993345 pairs), JSON overhead was ~10 seconds.
     * Native struct generation reduces this to ~1 second.
     *
     * @param atoms Atomic numbers
     * @param geometry_bohr Coordinates in Bohr
     * @return Vector of GFNFFDispersion structs
     */
    std::vector<GFNFFDispersion> GenerateDispersionPairsNative(
        const std::vector<int>& atoms, const Matrix& geometry_bohr);

    // Individual parameter accessors
    double getC6(int atom_i, int atom_j, int ref_i = 0, int ref_j = 0) const;
    double getR4OverR2(int atom) const;
    double getSqrtZr4r2(int atom) const;
    double getAtomicPolarizability(int atom, int frequency_index = 1) const;

    // Claude Generated (Jan 31, 2026): Set topology charges for zeta scaling
    // Reference: Fortran gfnff_ini.f90:789 - f1 = zeta(ati, topo%qa(i))
    // GFN-FF computes zetac6 ONCE during initialization using topology charges (topo%qa),
    // which are calculated with INTEGER neighbor counts, NOT fractional CN from geometry.
    void setTopologyCharges(const Vector& charges) { m_topology_charges = charges; }

    // Claude Generated (Feb 15, 2026): dc6dcn computation for dispersion CN gradient
    // Reference: Fortran gfnff_gdisp0.f90:174-210, 262-305
    // @param skip_dc6dcn  If true, skip O(N²) dc6dcn matrix (GPU computes per-pair)
    void updateCNValuesForGradient(const std::vector<double>& cn, CxxThreadPool* pool = nullptr,
                                    int num_threads = 1, bool skip_dc6dcn = false);
    const Matrix& getDC6DCN() const { return m_dc6dcn; }

    // Claude Generated (March 2026): GPU dc6dcn Phase 2 — expose weight arrays
    const std::vector<std::vector<double>>& getGaussianWeights() const { return m_gaussian_weights; }
    const std::vector<std::vector<double>>& getGaussianWeightDerivatives() const { return m_gaussian_weight_derivatives; }
    const std::vector<double>& getC6FlatCache() const { return m_c6_flat_cache; }
    const std::vector<int>& getRefN() const { return m_refn; }
    const std::vector<std::vector<double>>& getRefCN() const { return m_refcn; }

    // Claude Generated: Made public for native D4 fallback in GFN2 (Dec 2025) and ATM triples (Mar 2026)
    double getChargeWeightedC6(int Zi, int Zj, size_t atom_i, size_t atom_j) const;
    double calculateTripleScale(int i, int j, int k) const;

    // Claude Generated (2026): D4 q-response chain rule (AP ∂q/∂x)
    // Per-atom zeta scaling and its closed-form charge derivative, used by
    // D4Evaluator to assemble dE_D4/dq. Both forward to GFNFFParameters.
    double getZeta(int Z, double q) const;
    double getZetaDerivative(int Z, double q) const;

    // Claude Generated (AP6b exact D4 port, 2026): tblite/dftd4-exact charge-weighted
    // C6 for native GFN2. Unlike getChargeWeightedC6 (CN-only, single per-atom zeta
    // prefactor — the GFN-FF approximation), this applies the charge-dependent zeta
    // PER REFERENCE STATE with the dftd4 weighting (wf=6, ngw multi-gaussian, gi=eta·gc,
    // gc=2, qref=refq[ref]+zeff). Mirrors dftd4 model.f90 weight_references+get_atomic_c6.
    //   C6  = ΣΣ (gw_i^ri·ζ_i^ri)(gw_j^rj·ζ_j^rj)·c6ref
    // Optional outputs (pass nullptr to skip): dC6/dq, dC6/dCN, and the second
    // charge derivative ∂²C6/∂q² (diagonal) for the self-consistent CPSCF kernel.
    // qi/qj are the SCF Mulliken atomic charges; cni/cnj the D4 coordination numbers.
    struct C6Gfn2 {
        double c6 = 0.0;
        double dc6dqi = 0.0, dc6dqj = 0.0;     // ∂C6/∂q
        double dc6dcni = 0.0, dc6dcnj = 0.0;   // ∂C6/∂CN
        double d2c6dqi2 = 0.0, d2c6dqj2 = 0.0; // ∂²C6/∂q² (diagonal, for CPSCF kernel)
    };
    // Per-atom CN is read internally from m_cn_values (same CN that drives the
    // existing CN-gaussian weights), so call after GenerateParameters().
    C6Gfn2 weightedC6Gfn2(int Zi, int Zj, size_t atom_i, size_t atom_j,
                          double qi, double qj,
                          bool want_grad = false, bool want_hess = false) const;

    // The charge vector that actually drives zetac6: topology charges if set
    // (GFN-FF path), otherwise the geometry-dependent EEQ charges (GFN2 path).
    // The dE_D4/dq term must use exactly these charges for consistency.
    const Vector& getZetaCharges() const {
        return (m_topology_charges.size() > 0) ? m_topology_charges : m_eeq_charges;
    }

    // Non-owning access to the EEQ solver that produced m_eeq_charges, so the
    // q-response path can reuse its cached factorisation for dq/dx (Phase 2).
    EEQSolver* getEEQSolver() const { return m_eeq_solver.get(); }
    const Vector& getEEQCharges() const { return m_eeq_charges; }

    // Claude Generated (2026): D4 q-response chain rule (AP ∂q/∂x, Phase 2)
    // Enable the canonical single-shot dftd4 EEQ charge model for zetac6. When
    // on (and no GFN-FF topology charges are set), GenerateParameters() fills
    // m_eeq_charges from D4ChargeModel — a single smooth linear system whose
    // dq/dx has a closed form. GFN-FF leaves this off (uses topology charges).
    void setUseD4SingleShotEEQ(bool on) { m_use_d4_single_shot_eeq = on; }
    bool usesD4SingleShotEEQ() const { return m_use_d4_single_shot_eeq; }

    // Claude Generated (2026): use the dftd4 EN-weighted covalent CN (no log-cap)
    // for the D4 C6 interpolation instead of the GFN-FF erf-CN. GFN2 turns this
    // on to match tblite's get_coordination_number(rcov, en); GFN-FF leaves it
    // off (its Fortran reference uses the log-capped, EN-free CN). See d4_ncoord.h.
    void setD4CovalentCN(bool on) {
        if (m_use_d4_covalent_cn != on) {
            m_use_d4_covalent_cn = on;
            // The α-zeta correction in computeC6Reference is gated by this flag,
            // so the cached C6 reference matrix must be rebuilt on toggle.
            m_c6_reference_cached = false;
        }
    }
    bool usesD4CovalentCN() const { return m_use_d4_covalent_cn; }

    // Accumulate the D4 charge-response gradient Σ_A dEdq(A)·∂q_A/∂R into
    // grad_out (Hartree/Bohr). Valid only when the single-shot EEQ model was
    // used for the current geometry; otherwise a no-op.
    void addChargeResponseGradient(const Vector& dEdq, Matrix& grad_out) const {
        if (m_use_d4_single_shot_eeq && m_d4_charge_model.valid())
            m_d4_charge_model.addChargeResponseGradient(dEdq, grad_out);
    }

    // Claude Generated (Apr 2026): P1a — CN-change threshold cache for Gaussian weights
    // Skip recomputation of gw/dgw/dc6dcn when CN changes < d4_cn_cache_threshold (MD optimization)
    bool canSkipGaussianWeightsUpdate(const std::vector<double>& new_cn) const;
    void recordCNValues(const std::vector<double>& cn);

private:
    void initializeReferenceData();
    void calculateFrequencyDependentPolarizabilities();

    // OLD: Crude neutral-atom approximation (deprecated)
    double getEffectiveC6(int atom_i, int atom_j) const;

    // Claude Generated (Dec 27, 2025): Weight caching optimization
    void precomputeGaussianWeights(CxxThreadPool* pool = nullptr, int num_threads = 1);

    // Claude Generated (Dec 27, 2025): C6 reference matrix pre-computation
    void precomputeC6ReferenceMatrix();
    double computeC6Reference(int elem_i, int elem_j, int ref_i, int ref_j) const;

    // getChargeWeightedC6 and calculateTripleScale moved to public section (March 2026)

    // Reference data from GFN-FF Fortran implementation
    // constexpr ensures inline definition (ODR-safe for C++14 and C++17)
    static constexpr int MAX_ELEM = 118;
    static constexpr int MAX_REF = 7;
    static constexpr int N_FREQ = 23;  // Frequency grid points
    static constexpr int N_REFQ = 17;  // Reference charge states

    // D4 reference data from external/gfnff/src/dftd4param.f90
    std::vector<int> m_refn;    // Number of reference systems per element
    std::vector<std::vector<double>> m_refq, m_refh;  // Reference charges and hydrogen counts (m_refh = dftd4 hcount, integer-valued)
    std::vector<std::vector<double>> m_refh_charges;  // dftd4 refh table: H-atom partial charges in compound reference states (used by GFN2 α-correction)
    std::vector<std::vector<double>> m_refcn;   // Reference CN — dftd4 'refcn', used ONLY for ngw bucketing (set_refgw)
    std::vector<std::vector<double>> m_refcovcn; // dftd4 'refcovcn' — covalent reference CN for the CN-Gaussian weighting (set_refcn -> model%cn -> weight_references). Distinct from m_refcn.

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

    // Claude Generated (2026): single-shot EEQ for D4 q-response (Phase 2)
    bool m_use_d4_single_shot_eeq = false;
    mutable curcuma::dispersion::D4ChargeModel m_d4_charge_model;

    // Claude Generated (2026): GFN2 uses the dftd4 EN-weighted covalent CN for
    // the C6 interpolation (see setD4CovalentCN / d4_ncoord.h).
    bool m_use_d4_covalent_cn = false;

    // Claude Generated (Jan 31, 2026): Topology charges for zeta scaling
    // Reference: Fortran gfnff_ini.f90:789 - f1 = zeta(ati, topo%qa(i))
    // GFN-FF uses topology-based charges (topo%qa) for zetac6 scaling, NOT geometry-dependent EEQ
    Vector m_topology_charges;  // Phase 1 topology charges from GFNFF::TopologyInfo

    // Molecular CN calculation (December 2025 Phase 2 - CN+Charge weighting)
    std::vector<double> m_cn_values;  // Cached GFNFFCN values for current geometry

    // Claude Generated (Dec 27, 2025): Cached Gaussian weights for C6 interpolation
    // Performance optimization: Pre-compute CN+charge weights once per atom instead of per pair
    // Expected speedup: 5-10x for molecules with 100+ atoms (same as D3)
    std::vector<std::vector<double>> m_gaussian_weights;  // [atom_idx][ref_idx]
    bool m_weights_cached = false;

    // Claude Generated (Feb 7, 2026): Phase 7a - Dense flat array for C6 reference values
    // Replaces unordered_map to eliminate hash overhead (40-100 cycles per find())
    // Layout: [elem_i * MAX_ELEM * MAX_REF * MAX_REF + elem_j * MAX_REF * MAX_REF + ref_i * MAX_REF + ref_j]
    // Size: 118 * 118 * 7 * 7 = 406,952 doubles = ~3.1 MB
    std::vector<double> m_c6_flat_cache;
    bool m_c6_reference_cached = false;

    // Dense array index for C6 reference lookup (no hash, just arithmetic)
    inline size_t c6FlatIndex(int elem_i, int elem_j, int ref_i, int ref_j) const {
        return static_cast<size_t>(elem_i) * MAX_ELEM * MAX_REF * MAX_REF
             + static_cast<size_t>(elem_j) * MAX_REF * MAX_REF
             + static_cast<size_t>(ref_i) * MAX_REF
             + static_cast<size_t>(ref_j);
    }

    // Claude Generated (Feb 7, 2026): Phase 7b - Dominant reference indices per atom
    // Pre-computed during weight calculation: only refs with weight > threshold
    // Eliminates near-zero weight iterations in inner loop of getChargeWeightedC6()
    std::vector<std::vector<int>> m_dominant_refs;  // [atom_idx] → list of significant ref indices

    // Claude Generated (Feb 15, 2026): dc6dcn computation helpers
    void computeGaussianWeightDerivatives(CxxThreadPool* pool = nullptr, int num_threads = 1);
    void computeDC6DCN(CxxThreadPool* pool = nullptr, int num_threads = 1);

    // Claude Generated (Feb 15, 2026): Gaussian weight derivatives and dc6dcn matrix
    // dgwdcn[atom_idx][ref_idx] = d(normalized_weight(ref))/d(CN(atom))
    std::vector<std::vector<double>> m_gaussian_weight_derivatives;
    // dc6dcn(i,j) = dC6(i,j)/dCN(i) - asymmetric matrix
    Matrix m_dc6dcn;
    bool m_dc6dcn_computed = false;

    // Claude Generated (Apr 2026): P1a — CN-change threshold cache for Gaussian weights
    // Skip recomputation of gw/dgw/dc6dcn when CN changes < threshold (MD optimization)
    std::vector<double> m_prev_cn_values;  // Previous step's CN values
    bool m_cn_cached = false;              // Whether m_prev_cn_values is valid

    // Mathematical constants
    static constexpr double PI = 3.14159265358979323846;
    static constexpr double THREE_OVER_PI = 3.0 / PI;
    static constexpr double HALF_OVER_PI = 0.5 / PI;
};