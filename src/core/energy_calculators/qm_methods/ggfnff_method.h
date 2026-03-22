/*
 * <GGFNFFComputationalMethod — GPU-accelerated GFN-FF wrapper>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (March 2026): ComputationalMethod adapter for ggfnff.
 * Available only when compiled with USE_CUDA=ON.
 *
 * Usage:
 *   ./curcuma -sp mol.xyz -method ggfnff
 */

#pragma once

#ifdef USE_CUDA

#include "../computational_method.h"
#include "../ff_methods/gfnff.h"
#include "../ff_methods/ff_workspace.h"
#include "../ff_methods/cuda/ff_workspace_gpu.h"

#include <memory>
#include <unordered_map>

/**
 * @brief GPU-accelerated GFN-FF via CUDA (method name: "ggfnff")
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
 *   1. m_gfnff->prepareCNAndEEQ(gradient)  — CN + EEQ on CPU
 *   2. Distribute state to GPU + CPU-residual workspaces
 *   3. m_gfnff->updateHBXBIfNeeded(m_cpu_residual)  — dynamic HB/XB
 *   4. m_gpu_workspace->calculate()  — GPU terms + CPU residual automatically
 */
class GGFNFFComputationalMethod : public ComputationalMethod {
public:
    explicit GGFNFFComputationalMethod(const std::string& method_name, const json& config);
    ~GGFNFFComputationalMethod();

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

private:
    /**
     * @brief Initialize GPU + CPU-residual workspaces from GFNFF parameter set.
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
    std::unique_ptr<FFWorkspace>    m_cpu_residual;  ///< CPU workspace for HB/XB/ATM/BATM/sTors
    std::unique_ptr<FFWorkspaceGPU> m_gpu_workspace; ///< GPU workspace (holds raw ptr to m_cpu_residual)

    json             m_parameters;
    std::string      m_method_name;
    std::vector<int> m_atom_types;   ///< Element numbers (Z) — stored from setMolecule()
    bool             m_initialized   = false;
    bool             m_has_error     = false;
    std::string      m_error_message;
    double           m_last_energy   = 0.0;
    Matrix           m_cached_gradient; ///< Cached gradient (copied from GPU workspace after calculate)

    // Pre-allocated buffers for per-step HB coordination number computation.
    // Avoids heap allocations on CUDA-corrupted heap during MD/Opt iterations.
    std::vector<double> m_hb_cn_values;                ///< [n_bonds] HB CN per bond
    std::unordered_map<int, double> m_hb_cn_map;       ///< H-atom → CN (cleared+reused each step)
};

#endif // USE_CUDA
