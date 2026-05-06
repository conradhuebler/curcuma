/*
 * <FFWorkspaceGPU - CUDA-Accelerated GFN-FF Workspace>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (March 2026): GPU counterpart to FFWorkspace.
 * Accepts a GFNFFParameterSet, uploads interaction lists as SoA to GPU
 * once at construction, then for each step: uploads geometry/charges/CN,
 * launches all CUDA kernels, downloads results.
 *
 * All gradient contributions (including CN chain-rule and dispersion dEdcn)
 * are computed entirely on GPU.  No CPU postprocessing needed.
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
#include "../gfnff.h"                  // GFNFF::EEQGPUParams (WP2)

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
 * Claude Generated (March 2026): GPU version of FFWorkspace for GFN-FF GPU path.
 * All energy terms and gradient contributions run on GPU.
 *
 * Usage:
 *   FFWorkspaceGPU gpu_ws(params, natoms, atom_types);
 *   // Per optimisation step:
 *   gpu_ws.setGeometry(geom);           // Bohr row-major N×3
 *   gpu_ws.setEEQCharges(q);            // N EEQ charges
 *   gpu_ws.setD3CN(cn);                 // N D3 coordination numbers
 *   gpu_ws.updateDispersionDC6DCN(dc6dcn); // N×N dc6/dcn matrix
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
    // Per-step state updates
    // =========================================================================

    /// Set current geometry (Bohr, N×3 row-major Eigen matrix)
    void setGeometry(const Matrix& geom);

    /// Set dynamic EEQ charges (geometry-dependent, size N)
    void setEEQCharges(const Vector& q);

    /// WP5-A: D2D path — charges already on GPU after GPU Schur complement.
    /// Enqueues a D2D copy on the main workspace stream; no H2D upload happens
    /// in launchChargeDependentAndFinish() for this step.
    /// Caller must ensure d_src is valid (EEQ stream fully synced) before calling.
    void setEEQDeviceCharges(const double* d_src);

    /// Set topology charges (size N — used for BATM kernel)
    void setTopologyCharges(const Vector& q);

    /// Set D3 coordination numbers (size N, for CN-dependent bond r0)
    void setD3CN(const Vector& cn);

    /**
     * @brief Set CN factors needed for Coulomb TERM 1b gradient.
     * @param cn   D3 CN (size N) — stored for k_subtract_qtmp
     * @param cnf  CNF factors for Coulomb TERM 1b (size N)
     */
    void setCNDerivatives(const Vector& cn, const Vector& cnf,
                          const std::vector<SpMatrix>& dcn);

    // =========================================================================
    // GPU CN Computation (Phase 1: GPU migration)
    // Claude Generated (March 2026): Compute CN on GPU
    // =========================================================================

    /**
     * @brief Compute GFN-FF coordination numbers on GPU.
     *
     * Computes CN_raw and CN_final from current geometry and atom types.
     * Replaces CPU CNCalculator::calculateGFNFFCN().
     * Results are downloaded to pinned buffer — use getCNPinnedBuffer()
     * to access the data.  No Eigen allocations (heap-corruption safe).
     *
     * @param atom_types  Atomic numbers (1-based, size N)
     */
    void computeCN(const std::vector<int>& atom_types);

    /**
     * @brief Finalize deferred CN download after Phase 4 launch (G-P1, Apr 2026).
     * Syncs the main CUDA stream and copies CN_final from pinned buffer to out_cn.
     * Call after prepareAndLaunchChargeIndependent() so the GPU finishes CN during
     * the Phase 4 launch overhead (~0.3ms), eliminating the 4ms sleep penalty.
     * @param out_cn Output vector resized to N doubles (log-transformed CN)
     */
    void finalizeCNForCPU(Vector& out_cn);

    /**
     * @brief Get pointer to pinned CN_final buffer (valid after computeCN()).
     * Caller copies into their pre-allocated Vector via memcpy.
     * @return Pointer to N doubles (log-transformed CN values)
     */
    const double* getCNPinnedBuffer() const { return m_h_cn_final; }

    /**
     * @brief Get pointer to pinned CN_raw buffer (valid after computeCN()).
     * Raw (un-squashed) CN values needed for dlogdcn computation.
     * @return Pointer to N doubles (raw erf-sum CN values)
     */
    const double* getCNRawPinnedBuffer() const { return m_h_cn_raw; }

    /**
     * @brief Get pointer to pinned dlogdcn buffer (valid after computeCN() with gradient=true).
     * dlogdcn is computed on GPU by k_dlogdcn kernel (Apr 2026).
     * @return Pointer to N doubles (logistic derivative)
     */
    const double* getDlogdcnPinnedBuffer() const { return m_h_dlogdcn; }

    /**
     * @brief Check if GPU CN has been computed this step.
     */
    bool hasComputedCN() const { return m_cn_computed; }

    /**
     * @brief Reset CN computed flag (call before new geometry step).
     */
    void resetCNComputed() { m_cn_computed = false; }

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

    /**
     * @brief Generate CN pair list entirely on GPU (Apr 2026).
     *
     * Two-pass kernel: count valid pairs, allocate buffers, write (i,j,rcov_sum).
     * Replaces CPU generateCNPairList() O(N^2) loop.
     * Coordinates and atom types must already be on GPU.
     */
    void generateCNPairListOnGPU();

    /// Set baseline energy (e0 from parameter set, added to total)
    void setE0(double e0);

    /// Update per-bond HB coordination numbers (for egbond_hb alpha modulation)
    void updateBondHBCN(const std::vector<double>& hb_cn_values);

    /**
     * @brief Re-upload HB alpha (H,B) pair list after dynamic HB re-detection.
     *
     * Claude Generated (Apr 2026): Replaces the HBAlphaSoA built at construction time
     * with updated (H,B) pairs from the new HB list.  Must be called after
     * updateHBonds() when consumeHBXBUpdate() returns true.
     *
     * @param bond_hb_data  New flat BondHBEntry list (from GFNFF::rebuildBondHBData)
     * @param atom_types    Atomic numbers (1-based, size N) for covalent radius lookup
     */
    void updateHBAlphaPairs(const std::vector<BondHBEntry>& bond_hb_data,
                             const std::vector<int>& atom_types);

    /**
     * @brief Re-upload per-bond nr_hb and hb_H_atom arrays after dynamic HB re-detection.
     *
     * Claude Generated (Apr 2026): Overwrites the BondSoA nr_hb and hb_H_atom buffers
     * on the GPU so the k_bonds kernel uses updated HB participation counts.
     * Must be called after updateHBAlphaPairs() when consumeHBXBUpdate() returns true.
     *
     * @param nr_hb       Per-bond HB count (size = bond count)
     * @param hb_H_atom   Per-bond H atom index, -1 if not an HB bond (size = bond count)
     */
    void updateBondHBMetadata(const std::vector<int>& nr_hb,
                               const std::vector<int>& hb_H_atom);

    /// Compute HB CN per-bond entirely on GPU (Apr 2026).
    /// Uses HBAlphaSoA + BondSoA; no CPU loop or H2D upload.
    void computeBondHBCN();

    /// Re-upload HBond SoA after dynamic re-detection (called after updateHBXBIfNeeded)
    void updateHBonds(const std::vector<GFNFFHydrogenBond>& hbonds,
                      const std::vector<int>& atom_types);

    /// Re-upload XBond SoA after dynamic re-detection
    void updateXBonds(const std::vector<GFNFFHalogenBond>& xbonds,
                      const std::vector<int>& atom_types);

    /// Get last uploaded HBond list (for CPU vs GPU comparison debugging)
    const std::vector<GFNFFHydrogenBond>& getLastHBonds() const { return m_last_hbonds; }

    /// Get last uploaded XBond list (for CPU vs GPU comparison debugging)
    const std::vector<GFNFFHalogenBond>& getLastXBonds() const { return m_last_xbonds; }

    /**
     * @brief Provide Coulomb self-energy parameters (size N each).
     *
     * Called automatically from constructor if GFNFFParameterSet contains
     * Coulomb data.  Can be overridden by GFNFF::Calculation().
     */
    void setCoulombSelfEnergyParams(const Vector& chi_base, const Vector& gam,
                                     const Vector& alp,
                                     const Vector& chi_static);

    /**
     * @brief Upload per-pair dc6/dcn values for dispersion dEdcn on GPU.
     *
     * Extracts dc6dcn(i,j) and dc6dcn(j,i) for each dispersion pair from
     * the full N×N matrix and uploads to GPU SoA arrays.
     * Must be called every gradient step (dc6dcn is CN-dependent).
     *
     * @param dc6dcn  N×N matrix: dc6dcn(i,j) = dC6(i,j)/dCN(i)
     */
    void updateDispersionDC6DCN(const Matrix& dc6dcn);

    // =========================================================================
    // Phase 2: GPU dc6dcn per-pair computation (Claude Generated March 2026)
    // =========================================================================

    /**
     * @brief Upload C6 reference table and element refn counts (one-time at init).
     * @param c6_flat  Flat C6 reference array [MAX_ELEM² × MAX_REF²]
     * @param refn     Number of reference states per element [MAX_ELEM]
     */
    void uploadC6ReferenceTable(const std::vector<double>& c6_flat,
                                 const std::vector<int>& refn);

    /**
     * @brief Upload Gaussian weights + derivatives and compute dc6dcn per pair on GPU.
     * Replaces: CPU computeDC6DCN() O(N²) + updateDispersionDC6DCN() extraction.
     * @param gw   Nested [N][nref] Gaussian weights (will be flattened to [N×MAX_REF])
     * @param dgw  Nested [N][nref] weight derivatives
     */
    void computeDC6DCNOnGPU(const std::vector<std::vector<double>>& gw,
                             const std::vector<std::vector<double>>& dgw);

    /**
     * @brief Upload reference CN values (one-time at init).
     * @param refcn  Nested [MAX_ELEM][nref] reference CN values
     */
    void uploadRefCN(const std::vector<std::vector<double>>& refcn);

    /**
     * @brief Compute Gaussian weights + dc6dcn entirely on GPU (Phase 6).
     * Replaces: CPU precomputeGaussianWeights() + computeGaussianWeightDerivatives()
     * + computeDC6DCNOnGPU(). CN values must be on GPU (from computeCN() or setD3CN()).
     * @param use_cn_final If true, read CN from d_cn_final (current step, written by computeCN).
     *                     If false (default), read from d_cn (updated in prepareAndLaunch).
     *                     Set use_cn_final=true to call BEFORE prepareAndLaunchChargeIndependent()
     *                     so dc6dcn is ready for Phase 1 dispersion overlap with EEQ.
     */
    void computeGaussianWeightsOnGPU(bool use_cn_final = false);

    // =========================================================================
    // Term enable flags (match FFWorkspace API)
    // =========================================================================

    void setDispersionEnabled(bool v);
    void setHBondEnabled(bool v);
    void setRepulsionEnabled(bool v);
    void setCoulombEnabled(bool v);

    /// Claude Generated (April 2026): Set unit cell for PBC minimum image convention.
    /// cell_bohr: 3×3 matrix (column-major, in Bohr). Uploads to GPU constant memory.
    void setUnitCell(const double* cell_bohr_9, const double* cell_bohr_inv_9, bool has_pbc);

    /// Set verbosity for diagnostic snapshot downloads (only >= 3 triggers snapshot D2H)
    void setVerbosity(int v) { m_verbosity = v; }

    /// Enable/disable FP32 mixed precision for repulsion, BATM, and XBonds kernels.
    /// Default: enabled. FP64 fallback available for validation.
    void setMixedPrecision(bool enable);

    /// G3a (Apr 2026): Set fixed GPU kernel block size (0 = adaptive default).
    /// Valid: 0, 32, 64, 128, 256, 512. 0 uses adaptive logic (32/128/256/512).
    /// 512 maximizes occupancy but may be slower for small systems.
    void setBlockSize(int block_size) { m_block_size = block_size; }
    int getBlockSize() const { return m_block_size; }

    // =========================================================================
    // GPU Topology Displacement Check (Claude Generated March 2026)
    // =========================================================================

    /**
     * @brief Check if any atom moved more than threshold (Bohr) from reference geometry.
     *
     * Uses flag-based GPU kernel: each thread checks one atom, atomicOr sets flag
     * if any atom exceeds squared threshold. Single int D2H transfer.
     *
     * @param threshold  Maximum allowed displacement in Bohr (0.5 for topology)
     * @return true if any atom exceeded the threshold
     */
    bool checkDisplacement(double threshold);

    /// Copy current coords SoA to ref_coords SoA (device-to-device, call after topology update)
    void updateReferenceGeometry();

    // =========================================================================
    // Main calculation
    // =========================================================================

    /**
     * @brief Run all GPU kernels and return total energy (Hartree).
     *
     * Execution order:
     *   1. Zero d_energy / d_grad / d_dEdcn on GPU
     *   2. Upload: charges, CN, topology charges
     *   3. k_dispersion (energy + gradient + dEdcn), k_repulsion (×2), k_coulomb
     *   4. k_bonds (energy + gradient + dEdcn), k_angles, k_dihedrals, k_inversions
     *   5. k_storsions, k_batm, k_atm, k_xbonds, k_hbonds
     *   6. k_coulomb_self (energy only)
     *   7. k_subtract_qtmp + k_cn_chainrule (gradient postprocess)
     *   8. Download energy scalars + gradient matrix
     *   9. Return E_total + e0
     *
     * @param gradient  If true compute gradient (N×3 Hartree/Bohr)
     * @return Total GFN-FF energy in Hartree
     */
    double calculate(bool gradient);

    // =========================================================================
    // Phase 3: Split calculation for CPU/GPU overlap (Claude Generated March 2026)
    //
    // Architecture: Allows the O(N³) EEQ solver on CPU to run in parallel with
    // charge-independent GPU kernels. Only k_coulomb, k_coulomb_self, and
    // k_subtract_qtmp require EEQ charges — everything else can fire immediately.
    //
    // Usage (in gfnff_gpu_method.cpp):
    //   gpu->setGeometry(geom);
    //   gpu->computeCN(atom_types);
    //   gpu->setD3CN(cn);
    //   gpu->prepareAndLaunchChargeIndependent(gradient);  // fire-and-forget
    //   // ... CPU EEQ runs here in parallel ...
    //   gpu->setEEQCharges(charges);
    //   energy = gpu->launchChargeDependentAndFinish(gradient);
    // =========================================================================

    /**
     * @brief Zero accumulators, upload CN/geometry, launch all charge-independent kernels.
     *
     * Launches on 3 streams: dispersion+repulsion (sA), bonds+angles+dihedrals+
     * inversions+storsions+hb_alpha (sB), batm+atm+xbonds+hbonds (sC).
     * Does NOT upload or use EEQ charges. Returns immediately (non-blocking).
     *
     * @param gradient  If true, kernels compute gradient contributions
     */
    void prepareAndLaunchChargeIndependent(bool gradient);

    /**
     * @brief Upload EEQ charges, launch charge-dependent kernels, download results.
     *
     * Waits for charge-independent kernels to complete (via stream events),
     * then launches k_coulomb, k_coulomb_self, k_subtract_qtmp, k_cn_chainrule.
     * Downloads energy and gradient. Blocking call.
     *
     * @param gradient  If true, compute and download gradient
     * @return Total GFN-FF energy in Hartree
     */
    double launchChargeDependentAndFinish(bool gradient);

    // =========================================================================
    // Results
    // =========================================================================

    const Matrix&             gradient()         const;
    const FFEnergyComponents& energyComponents() const;

    // dEdcn totals (for external CN chain-rule, if needed)
    const Vector& dEdcnTotal() const;

    /// Gradient before CN chain-rule (for diagnostic comparison)
    const Matrix& gradientBeforeCN() const { return m_grad_before_cn; }

    // =========================================================================
    // Diagnostics
    // =========================================================================

    /// Number of dispersion pairs uploaded to GPU
    int dispersionCount() const;

    /// Number of bond pairs uploaded to GPU
    int bondCount() const;

    /// DEPRECATED: use getDeviceXPtr/YPtr/ZPtr() instead.
    /// Returns nullptr — coords are now in SoA layout (Phase 10, March 2026).
    const double* getDeviceCoordsPtr() const;

    /// Device pointer to x-coordinate array [N] (Bohr, SoA).
    /// Claude Generated (Phase 10, March 2026): SoA for coalesced warp reads.
    const double* getDeviceXPtr() const;
    const double* getDeviceYPtr() const;
    const double* getDeviceZPtr() const;

    /// Synchronize main stream to ensure all prior uploads (coords, CN) are visible.
    /// Claude Generated (March 2026): Required before cross-stream reads (e.g. EEQSolverGPU).
    void synchronizeMainStream();

    // =========================================================================
    // WP2: GPU-side EEQ topology parameters (Claude Generated May 2026)
    // =========================================================================

    /// Upload topology-constant EEQ parameters to GPU (once per topology build).
    /// Enables k_build_eeq_rhs to construct rhs_atoms[i] entirely on GPU each step.
    /// Call after InitialiseMolecule() and after each full topology recalculation.
    void uploadEEQTopologyParams(const GFNFF::EEQGPUParams& params);

    /// Returns true if uploadEEQTopologyParams() has been called at least once.
    bool isEEQTopoValid() const;

    /// Device pointer to EEQ RHS vector [N] (valid after k_build_eeq_rhs ran).
    const double* getDeviceRHSPtr() const;
    /// Device pointer to EEQ fragment target charges [nfrag] (topology-constant).
    const double* getDeviceRHSConstraintsPtr() const;
    /// Device pointer to alpha_corrected [N] (topology-constant).
    const double* getDeviceAlphaPtr() const;
    /// Device pointer to gam_corrected [N] (topology-constant).
    const double* getDeviceGamPtr() const;

    /// Invalidate the CUDA Graph — call whenever any SoA topology changes
    /// (full topology recalculation, HB/XB re-detection, or any n-value change).
    /// The next prepareAndLaunchChargeIndependent() call will re-capture the graph.
    /// Claude Generated (Phase 8, March 2026).
    void invalidateGraph();

    /// WP5-C (May 2026): Return the GPU-side dc6dcn skip flag from the previous step.
    /// True if max|ΔCN| < 0.01 between the last two steps, meaning dc6dcn recomputation
    /// can be skipped. Always false on the first step or after topology invalidation.
    bool shouldSkipDc6dcn() const { return m_dc6dcn_skip_pending; }

private:
    std::unique_ptr<FFWorkspaceGPUImpl> m_impl;

    // CPU-side per-step state
    int     m_natoms = 0;
    double  m_e0     = 0.0;

    // Coulomb self-energy parameters (O(N), extracted at init)
    Vector  m_coul_chi_base, m_coul_gam, m_coul_alp, m_coul_chi_static;

    // CN state for GPU upload and k_subtract_qtmp (Coulomb TERM 1b)
    // WP5-B (May 2026): m_cnf removed — CNF lives on GPU as impl.d_coul_cnf
    // (topology-constant). m_cn kept for legacy setD3CN() compatibility.
    Vector  m_cn;

    // GPU CN computation state
    // NOTE: No Eigen Vectors here — heap-corruption-safe.  CN data lives in pinned buffer only.
    double* m_h_cn_final = nullptr; ///< [N] pinned staging buffer for CN_final download
    double* m_h_cn_raw   = nullptr; ///< [N] pinned staging buffer for CN_raw download
    double* m_h_dlogdcn  = nullptr; ///< [N] pinned staging buffer for dlogdcn download (Apr 2026)
    bool    m_cn_computed = false; ///< GPU CN computed this step?
    std::vector<int> m_atom_types_cached; ///< Cached atom types for GPU CN

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
    bool    m_device_charges_ready = false;  ///< WP5-A: set by setEEQDeviceCharges(), skips H2D upload

    // WP5-C (May 2026): GPU-side D4 dc6dcn skip-check state
    bool    m_dc6dcn_skip_pending = false;   ///< skip flag from previous step (read in current step)

    // Last uploaded HB/XB bond lists (for CPU vs GPU comparison debugging)
    std::vector<GFNFFHydrogenBond> m_last_hbonds;
    std::vector<GFNFFHalogenBond> m_last_xbonds;

    // Host-side dispersion pair indices (for per-step dc6dcn extraction)
    std::vector<int> m_disp_idx_i_host;
    std::vector<int> m_disp_idx_j_host;

    // Per-step dc6dcn staging buffer (pre-allocated, reused)
    std::vector<double> m_h_dc6dcn_ij;
    std::vector<double> m_h_dc6dcn_ji;

    // Phase 2: GPU dc6dcn staging buffers (pre-allocated, reused per step)
    std::vector<double> m_h_gw_flat;   ///< [N * MAX_REF] flattened Gaussian weights
    std::vector<double> m_h_dgw_flat;  ///< [N * MAX_REF] flattened weight derivatives
    static constexpr int D4_MAX_REF = 7;

    // Term enable flags
    bool m_dispersion_enabled = true;
    bool m_hbond_enabled      = true;
    bool m_repulsion_enabled  = true;
    bool m_coulomb_enabled    = true;
    int  m_verbosity          = 0;

    // G3a (Apr 2026): GPU kernel block size override (0 = adaptive default)
    int  m_block_size         = 0;

    // Pre-allocated pinned staging buffers for async DMA transfers.
    // Claude Generated (March 2026): Pinned memory enables true async H2D/D2H via
    // cudaMemcpyAsync without CPU page-fault overhead. This is the primary per-step
    // performance optimization for iterative MD/Opt workloads.
    // Allocated via cudaMallocHost in constructor, freed via cudaFreeHost in destructor.
    double* m_h_x = nullptr;   ///< [N] pinned staging x-coords for geometry upload
    double* m_h_y = nullptr;   ///< [N] pinned staging y-coords for geometry upload
    double* m_h_z = nullptr;   ///< [N] pinned staging z-coords for geometry upload
    double* m_h_grad   = nullptr;   ///< [3*N] pinned staging for gradient download
    double* m_h_dEdcn_snap = nullptr; ///< [N] pinned staging for dEdcn snapshot download
    double* m_h_grad_snap  = nullptr; ///< [3*N] pinned staging for pre-CN grad snapshot
    // Per-kernel diagnostic pinned buffers (Claude Generated April 2026)
    double* m_h_grad_after_sA     = nullptr; ///< [3*N] grad after stream A (repulsion)
    double* m_h_grad_after_bonds  = nullptr; ///< [3*N] grad after k_bonds
    double* m_h_grad_after_angles = nullptr; ///< [3*N] grad after k_angles
    double* m_h_grad_after_sB     = nullptr; ///< [3*N] grad after all stream B
    double* m_h_grad_after_sC     = nullptr; ///< [3*N] grad after stream C (3-body)
    double* m_h_energies   = nullptr; ///< [14] pinned staging for per-term energy download

    // Results
    Matrix             m_result_gradient;
    Matrix             m_grad_before_cn;  ///< Gradient before CN chain-rule (diagnostic)
    FFEnergyComponents m_result_energy;
    Vector             m_dEdcn_total;
};

#endif // USE_CUDA
