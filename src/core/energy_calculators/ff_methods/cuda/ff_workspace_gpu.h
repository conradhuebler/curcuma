/*
 * <FFWorkspaceGPU - CUDA-Accelerated GFN-FF Workspace>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (March 2026): GPU counterpart to FFWorkspace.
 * Accepts a GFNFFParameterSet, uploads interaction lists as SoA to GPU
 * once at construction, then for each step: uploads geometry/charges/CN,
 * launches all CUDA kernels, downloads results.
 *
 * All gradient contributions (including CN chain-rule) are computed on GPU
 * for full consistency.  No CPU postprocessing needed for gradients.
 *
 * Interface is deliberately CUDA-agnostic so it can be forward-declared
 * from non-CUDA translation units.  All CUDA types are hidden in the Pimpl
 * (FFWorkspaceGPUImpl, defined in ff_workspace_gpu.cu).
 *
 * Minimum compute capability: 6.0 (Pascal) — native double atomicAdd.
 *
 * Reference: Spicher/Grimme J. Chem. Theory Comput. 2020 (GFN-FF)
 */

#pragma once

#ifdef USE_CUDA

#include "src/core/global.h"
#include "../ff_workspace.h"           // FFEnergyComponents, Matrix, Vector, SpMatrix
#include "../gfnff_parameters.h"       // GFNFFParameterSet

#include <Eigen/Sparse>
#include <memory>
#include <vector>

// ---------------------------------------------------------------------------
// Forward declaration of CUDA-internal implementation (defined in .cu)
// ---------------------------------------------------------------------------
struct FFWorkspaceGPUImpl;

// ---------------------------------------------------------------------------
// FFWorkspaceGPU — GPU-accelerated GFN-FF energy and gradient
// ---------------------------------------------------------------------------
/**
 * @brief GPU-accelerated drop-in replacement for FFWorkspace (GFN-FF only).
 *
 * Claude Generated (March 2026): GPU version of FFWorkspace for the `ggfnff`
 * method.  The CPU FFWorkspace remains unchanged; this class provides the
 * same setter/calculate/getter API so GGFNFFComputationalMethod can swap
 * between GPU and CPU at construction time.
 *
 * Usage:
 *   FFWorkspaceGPU gpu_ws(params, natoms, atom_types);
 *   // Per optimisation step:
 *   gpu_ws.setGeometry(geom);           // Bohr row-major N×3
 *   gpu_ws.setEEQCharges(q);            // N EEQ charges
 *   gpu_ws.setD3CN(cn);                 // N D3 coordination numbers
 *   gpu_ws.setCNDerivatives(cn, cnf, dcn);
 *   double E = gpu_ws.calculate(true);
 *   Matrix G = gpu_ws.gradient();       // N×3 Hartree/Bohr
 */
class FFWorkspaceGPU {
public:
    /**
     * @brief Upload GFN-FF parameters to GPU.
     * @param params    Complete GFNFFParameterSet from GFNFF::generateGFNFFParameterSet()
     * @param natoms    Number of atoms
     * @param atom_types Element numbers (1-based, size natoms) for angle/dihedral damping
     * @throws std::runtime_error on CUDA allocation failure or if no CUDA device is available
     */
    explicit FFWorkspaceGPU(const GFNFFParameterSet& params,
                             int natoms,
                             const std::vector<int>& atom_types);

    ~FFWorkspaceGPU();

    // Non-copyable (owns GPU resources)
    FFWorkspaceGPU(const FFWorkspaceGPU&)             = delete;
    FFWorkspaceGPU& operator=(const FFWorkspaceGPU&)  = delete;

    // =========================================================================
    // Per-step state updates (mirror of FFWorkspace interface)
    // =========================================================================

    /// Set current geometry (Bohr, N×3 row-major Eigen matrix)
    void setGeometry(const Matrix& geom);

    /// Set dynamic EEQ charges (geometry-dependent, size N)
    void setEEQCharges(const Vector& q);

    /// Set topology charges (size N — used for BATM; forwarded to CPU path only)
    void setTopologyCharges(const Vector& q);

    /// Set D3 coordination numbers (size N, for CN-dependent bond r0)
    void setD3CN(const Vector& cn);

    /**
     * @brief Set coordination numbers + sparse dcn derivatives for gradient.
     * @param cn   D3 CN (size N)
     * @param cnf  CNF factors for Coulomb TERM 1b (size N)
     * @param dcn  3 sparse N×N matrices: ∂CN_i/∂x_j, ∂CN_i/∂y_j, ∂CN_i/∂z_j
     * @deprecated Use setCNPairList() + setDlogDCN() instead for GPU-only gradient
     */
    void setCNDerivatives(const Vector& cn, const Vector& cnf,
                          const std::vector<SpMatrix>& dcn);

    /// Set dc6/dcn pointer for D4 dispersion CN gradient (CPU chain-rule only)
    /// @deprecated Dispersion dEdcn will be computed on GPU
    void setDC6DCNPtr(const Matrix* ptr);

    // =========================================================================
    // GPU-only CN chain-rule data (replaces sparse dcn matrices)
    // Claude Generated (March 2026): Full GPU gradient consistency
    // =========================================================================

    /**
     * @brief Set CN pair list for GPU CN chain-rule gradient kernel.
     *
     * All atom pairs (i,j) within CN cutoff, with pre-computed scaled
     * covalent radius sums. Uploaded to GPU once; geometry checked at runtime.
     *
     * @param idx_i      Atom i indices (size n_pairs)
     * @param idx_j      Atom j indices (size n_pairs)
     * @param rcov_sum   Scaled covalent radius sum per pair: 4/3*(rcov_i+rcov_j) in Bohr
     */
    void setCNPairList(const std::vector<int>& idx_i,
                       const std::vector<int>& idx_j,
                       const std::vector<double>& rcov_sum);

    /**
     * @brief Set per-atom dlogdcn (logistic squashing factor for CN chain-rule).
     * dlogdcn[i] = exp(cnmax) / (exp(cnmax) + exp(cn_raw[i]))
     * Must be called every step after CN computation.
     * @param dlogdcn  Per-atom squashing factor (size N)
     */
    void setDlogDCN(const Vector& dlogdcn);

    /// Set baseline energy (e0 from parameter set, added to total)
    void setE0(double e0);

    /// Update per-bond HB coordination numbers (for egbond_hb alpha modulation)
    void updateBondHBCN(const std::vector<double>& hb_cn_values);

    /// Re-upload HBond SoA after dynamic re-detection (called after updateHBXBIfNeeded)
    void updateHBonds(const std::vector<GFNFFHydrogenBond>& hbonds,
                      const std::vector<int>& atom_types);

    /// Re-upload XBond SoA after dynamic re-detection
    void updateXBonds(const std::vector<GFNFFHalogenBond>& xbonds,
                      const std::vector<int>& atom_types);

    /**
     * @brief Provide Coulomb self-energy parameters (size N each).
     *
     * Called automatically from constructor if GFNFFParameterSet contains
     * Coulomb data.  Can be overridden by GFNFF::Calculation().
     */
    void setCoulombSelfEnergyParams(const Vector& chi_base, const Vector& gam,
                                     const Vector& alp,     const Vector& cnf,
                                     const Vector& chi_static);

    // =========================================================================
    // Term enable flags (match FFWorkspace API)
    // =========================================================================

    void setDispersionEnabled(bool v);
    void setHBondEnabled(bool v);
    void setRepulsionEnabled(bool v);
    void setCoulombEnabled(bool v);

    // =========================================================================
    // Main calculation
    // =========================================================================

    /**
     * @brief Run all GPU kernels and return total energy (Hartree).
     *
     * Execution order:
     *   1. Zero d_energy / d_grad / d_dEdcn on GPU
     *   2. Async upload: geometry, charges, CN
     *   3. k_dispersion, k_repulsion (×2), k_coulomb
     *   4. k_bonds (accumulates dEdcn), k_angles, k_dihedrals, k_inversions
     *   5. cudaStreamSynchronize
     *   6. Download energy scalar + gradient matrix + dEdcn vector
     *   7. postProcessCPU(): Coulomb TERM 2+3 self-energy + CN chain-rule gradient
     *   8. Return E_total + e0
     *
     * H-bonds, X-bonds, ATM, BATM are forwarded to the CPU residual workspace
     * (if set via setCPUResidualWorkspace) — Phase 1 leaves them on CPU.
     *
     * @param gradient  If true compute gradient (N×3 Hartree/Bohr)
     * @return Total GFN-FF energy in Hartree
     */
    double calculate(bool gradient);

    // =========================================================================
    // Results
    // =========================================================================

    const Matrix&             gradient()         const;
    const FFEnergyComponents& energyComponents() const;

    // dEdcn totals (for external CN chain-rule, if needed)
    const Vector& dEdcnTotal() const;

    // =========================================================================
    // Optional: CPU residual workspace for H/X-bonds and three-body terms
    // =========================================================================

    /**
     * @brief Attach a CPU FFWorkspace for terms not yet on GPU (Phase 1).
     *
     * When set, calculate() calls cpu_ws->calculate(gradient) after the GPU
     * kernels and adds the CPU energy/gradient to the GPU results.  The CPU
     * workspace must already have its per-step state set by the caller.
     */
    void setCPUResidualWorkspace(FFWorkspace* cpu_ws);

    // =========================================================================
    // Diagnostics
    // =========================================================================

    /// Number of dispersion pairs uploaded to GPU
    int dispersionCount() const;

    /// Number of bond pairs uploaded to GPU
    int bondCount() const;

private:
    std::unique_ptr<FFWorkspaceGPUImpl> m_impl;

    // CPU-side per-step state (needed for postProcessCPU)
    int     m_natoms = 0;
    double  m_e0     = 0.0;

    // Coulomb self-energy parameters (O(N), extracted at init)
    Vector  m_coul_chi_base, m_coul_gam, m_coul_alp, m_coul_cnf, m_coul_chi_static;

    // CN chain-rule state (legacy sparse matrix path, deprecated)
    Vector              m_cn, m_cnf;
    std::vector<SpMatrix> m_dcn;

    // GPU CN chain-rule state (replaces sparse dcn matrices)
    Vector m_dlogdcn;                       ///< [N] logistic squashing factor
    std::vector<int>    m_cn_pair_i;        ///< CN pair atom i indices
    std::vector<int>    m_cn_pair_j;        ///< CN pair atom j indices
    std::vector<double> m_cn_pair_rcov;     ///< CN pair scaled cov. radius sum
    int                 m_cn_n_pairs = 0;   ///< Number of CN pairs
    bool                m_cn_pairs_on_gpu = false; ///< Whether CN pairs uploaded to GPU
    double              m_kn = -7.5;        ///< CN exponential decay constant

    // Dynamic per-step state
    Vector  m_eeq_charges;
    Vector  m_topology_charges;

    // DC6/DCN pointer for D4 CN gradient (legacy, will be moved to GPU)
    const Matrix* m_dc6dcn_ptr = nullptr;

    // CPU-side dispersion pairs (for dEdcn chain-rule — temporary until GPU migration)
    std::vector<GFNFFDispersion> m_dispersions_cpu;
    Matrix m_geometry_cpu;  ///< Last geometry (Bohr) for CPU-side dispersion dEdcn

    // Term enable flags
    bool m_dispersion_enabled = true;
    bool m_hbond_enabled      = true;
    bool m_repulsion_enabled  = true;
    bool m_coulomb_enabled    = true;

    // Optional CPU residual for H/X-bonds and three-body terms
    FFWorkspace* m_cpu_residual = nullptr;

    // Pre-allocated staging buffers (avoid heap allocs on corrupted heap)
    // Claude Generated (March 2026): CUDA operations corrupt C++ heap metadata.
    // Repeated std::vector allocations in setGeometry/calculate cause segfaults
    // during iterative MD/Opt. Pre-allocating these buffers once avoids the issue.
    std::vector<double> m_h_coords;  ///< [3*N] staging buffer for geometry upload
    std::vector<double> m_h_grad;    ///< [3*N] staging buffer for gradient download

    // Results
    Matrix             m_result_gradient;
    FFEnergyComponents m_result_energy;
    Vector             m_dEdcn_total;

    /**
     * @brief CPU counterpart of FFWorkspace::postProcess().
     *
     * Computes:
     *   - Coulomb TERM 2+3: self-energy + electronegativity contribution
     *   - dEdcn chain-rule gradient: gradient += dcn * (dEdcn - qtmp)
     *
     * @param gradient  If true, apply CN chain-rule gradient
     */
    void postProcessCPU(bool gradient);
};

#endif // USE_CUDA
