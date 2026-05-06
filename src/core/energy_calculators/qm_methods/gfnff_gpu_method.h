/*
 * <GFNFFGPUComputationalMethod — GPU-accelerated GFN-FF wrapper>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (March 2026): ComputationalMethod adapter for gfnff GPU path.
 * Available only when compiled with USE_CUDA=ON.
 *
 * Usage:
 *   ./curcuma -sp mol.xyz -method gfnff -gpu cuda    # Explicit GPU
 *   ./curcuma -sp mol.xyz -method gfnff -gpu auto   # GPU if available
 *   ./curcuma -sp mol.xyz -method gfnff             # CPU (default)
 */

#pragma once

#ifdef USE_CUDA

#include "../computational_method.h"
#include "../ff_methods/gfnff.h"
#include "../ff_methods/eeq_solver.h"  // EEQSolveMethod enum
#include "../ff_methods/cuda/ff_workspace_gpu.h"
#include "../ff_methods/cuda/eeq_solver_gpu.h"

#include <memory>

/**
 * @brief GPU-accelerated GFN-FF via CUDA (method name: "gfnff" with -gpu cuda)
 *
 * Claude Generated (March 2026): Clean GPU/CPU separation architecture.
 * GFNFF is a pure CPU class (no GPU knowledge). This wrapper orchestrates:
 *
 * Architecture:
 *   - m_gfnff          : GFNFF instance (CPU topology + EEQ charges + CN)
 *   - m_gpu_workspace  : FFWorkspaceGPU (bonds, angles, dihedrals, inversions,
 *                        dispersion, repulsion, Coulomb on GPU)
 *   - m_cpu_residual   : FFWorkspace (HB, XB, ATM, BATM, sTors on CPU)
 *   - FFWorkspaceGPU::calculate() automatically adds m_cpu_residual results
 *
 * Initialization flow:
 *   1. setMolecule() → m_gfnff->InitialiseMolecule() (topology, params)
 *   2. consumeCachedParameterSet() → split into GPU params + CPU residual params
 *   3. FFWorkspaceGPU(full_params) + FFWorkspace(residual_params)
 *   4. m_gpu_workspace->setCPUResidualWorkspace(m_cpu_residual)
 *
 * Per-step calculation (orchestrated here, NOT delegated to GFNFF::Calculation):
 *   1. GPU: computeCN() — CN on GPU (k_cn_compute kernel)
 *   2. CPU: prepareCNAndEEQ(gradient, gpu_only, &gpu_cn) — EEQ only (CN from GPU)
 *   3. Distribute state (charges, CN, geometry) to GPU workspace
 *   4. m_gfnff->updateHBXBIfNeeded() — dynamic HB/XB re-detection
 *   5. m_gpu_workspace->calculate() — all energy terms on GPU
 */
class GFNFFGPUComputationalMethod : public ComputationalMethod {
public:
    explicit GFNFFGPUComputationalMethod(const std::string& method_name, const json& config);
    ~GFNFFGPUComputationalMethod();

    // === ComputationalMethod interface ===

    bool setMolecule(const Mol& mol) override;
    double calculateEnergy(bool gradient = false) override;
    Matrix getGradient() const override;
    void copyGradientTo(Matrix& target) const override;
    Vector getCharges() const override;
    Vector getBondOrders() const override;
    Position getDipole() const override;
    bool updateGeometry(const Matrix& geometry) override;
    bool hasGradient() const override { return true; }
    bool isThreadSafe() const override { return false; }
    std::string getMethodName() const override { return m_method_name; }

    // === Configuration ===

    void setThreadCount(int threads) override;
    void setParameters(const json& params) override;
    json getParameters() const override;

    // === Error handling ===

    bool hasError() const override;
    void clearError() override;
    std::string getErrorMessage() const override;

    // === Energy decomposition ===

    json getEnergyDecomposition() const override;

    /// Expose GPU dEdcn for diagnostics (valid after calculateEnergy with gradient)
    const Vector& getGPUdEdcn() const;

    /// Expose GPU workspace for gradient diagnostics (e.g. gradientBeforeCN)
    FFWorkspaceGPU* getGPUWorkspace() const { return m_gpu_workspace.get(); }

    /// Get GPU CN result (valid after calculateEnergy)
    const Vector& getGPUCN() const { return m_gpu_cn_final; }

private:
    /**
     * @brief Initialize GPU workspace from GFNFF parameter set.
     * Called after m_gfnff->InitialiseMolecule() succeeds.
     * @return true on success; sets m_has_error on CUDA failure.
     */
    bool initGPUWorkspace();

    std::unique_ptr<GFNFF>          m_gfnff;
    // NOTE: GPU params intentionally leaked (raw ptr, never freed).
    // FFWorkspaceGPU's CUDA allocations corrupt adjacent heap metadata, making
    // the GFNFFParameterSet unfreeable.  Cost: ~100 KB one-time leak.
    // TODO: Investigate CUDA driver heap corruption root cause.
    GFNFFParameterSet*              m_gpu_params_leaked = nullptr;
    std::unique_ptr<FFWorkspaceGPU> m_gpu_workspace;

    // Claude Generated (March 2026): GPU EEQ solver (cuSOLVER Cholesky)
    std::unique_ptr<EEQSolverGPU>   m_eeq_gpu;

    // Pre-allocated buffers for GPU EEQ Schur complement (avoid per-step heap allocs)
    std::vector<double> m_eeq_z1;           ///< [N] A⁻¹ · b_atoms
    std::vector<double> m_eeq_Z2;           ///< [N*nfrag] A⁻¹ · C^T (column-major)
    std::vector<double> m_eeq_charges_gpu;  ///< [N] final charges from GPU path
    std::vector<double> m_schur_workspace;  ///< [nfrag*(nfrag+2)] CPU Schur: S matrix + rhs + lambda (no heap after CUDA init)

    // WP2: cached topology-constant EEQ data (avoid prepareEEQParametersForGPU per step)
    std::vector<int>    m_eeq_fraglist;          ///< [N] fragment IDs (topology-constant)
    std::vector<double> m_eeq_rhs_constraints;   ///< [nfrag] target charges (topology-constant)
    int                 m_eeq_nfrag = 1;         ///< number of fragments (topology-constant)

    json             m_parameters;
    std::string      m_method_name;
    std::vector<int> m_atom_types;   ///< Element numbers (Z) — stored from setMolecule()
    bool             m_allow_unconverged_charges = false; ///< Respect allow_unconverged_charges for GPU→CPU fallback (Claude Generated Apr 2026)
    bool             m_skip_phase2 = false;               ///< Skip Phase 2 EEQ refinement, use Phase 1 topology charges (Claude Generated Apr 2026)
    bool             m_initialized   = false;
    bool             m_has_error     = false;
    std::string      m_error_message;
    double           m_last_energy   = 0.0;
    Matrix           m_cached_gradient; ///< Cached gradient (copied from GPU workspace after calculate)

    // CN chain-rule pair list (generated once at init, used every gradient step)
    // Claude Generated (March 2026): Full GPU gradient consistency
    std::vector<int>    m_cn_pair_i;        ///< atom i indices
    std::vector<int>    m_cn_pair_j;        ///< atom j indices
    std::vector<double> m_cn_pair_rcov;     ///< scaled cov. radius sum (Bohr)
    bool                m_cn_pairs_generated = false;

    // GPU CN result (always used — GPU CN is the default path)
    Vector              m_gpu_cn_final;        ///< Cached GPU CN result

    // RMSD-based EEQ lazy refactorization (Apr 2026)
    // Caller tracks geometry RMSD from last full Cholesky refactorization.
    // m_eeq_rmsd_threshold == 0.0 → always refactorize (default, matches reference).
    // m_eeq_rmsd_threshold >  0.0 → lazy: reuse L when per-atom RMSD < threshold (Bohr).
    double m_eeq_rmsd_threshold = 0.0;  ///< Bohr; from eeq_rmsd_threshold param
    Matrix m_eeq_ref_geom;              ///< geometry (Bohr) at last EEQ refactorization
    bool   m_eeq_has_ref_geom = false;

    // EEQ Coulomb-matrix distance cutoff (Bohr). Default 0 = no cutoff, matches Fortran
    // goed_gfnff. Set non-zero only for performance experiments — produces HF-inconsistent
    // gradients vs. the un-truncated Coulomb energy and degrades MD energy conservation.
    double m_eeq_distance_cutoff = 0.0;  ///< Bohr; from eeq_distance_cutoff param

    // WP7-B (May 2026): GPU EEQ solver strategy for nfrag>1 systems.
    // SchurCholesky → WP5-A (nfrag=1) / WP7-A (nfrag>1) — exact, full N×N Cholesky.
    // Batched      → WP7-B (nfrag>1)                   — per-fragment Cholesky, drops cross-fragment Coulomb.
    // PCG          → WP7-C (nfrag>1)                   — iterative O(k·N²), warm-started.
    // LU / Auto    → CPU concepts; Auto resolves to PCG above pcg_large_threshold else SchurCholesky.
    EEQSolveMethod m_eeq_strategy = EEQSolveMethod::SchurCholesky;
    double m_eeq_batched_min_distance_bohr = 15.0;  ///< warn if min inter-frag distance < this

    // WP7-C (May 2026): GPU PCG parameters (mirror CPU EEQSolver defaults).
    int    m_eeq_pcg_max_iter   = 200;     ///< per-PCG-call iteration cap
    double m_eeq_pcg_tolerance  = 1e-10;   ///< convergence tolerance on |r|
    int    m_eeq_pcg_threshold  = 500;     ///< Auto strategy: PCG for N>=this (else cholesky)

    int  m_calc_count = 0;  ///< counts calculateEnergy() calls; first 5 always print timing

    /**
     * @brief Generate CN pair list from geometry and covalent radii.
     * Called once after initGPUWorkspace(). Pairs with rcov_sum < 2*max contribution.
     */
    void generateCNPairList(const Matrix& geom_bohr);
};

#endif // USE_CUDA
