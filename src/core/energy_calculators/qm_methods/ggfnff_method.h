/*
 * <GGFNFFComputationalMethod — GPU-accelerated GFN-FF wrapper>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (March 2026): ComputationalMethod adapter for ggfnff.
 * Thin wrapper around GFNFF + FFWorkspaceGPU.
 * Available only when compiled with USE_CUDA=ON.
 *
 * Usage:
 *   ./curcuma -sp mol.xyz -method ggfnff
 */

#pragma once

#ifdef USE_CUDA

#include "../computational_method.h"
#include "../ff_methods/gfnff.h"
#include "../ff_methods/cuda/ff_workspace_gpu.h"

#include <memory>

/**
 * @brief GPU-accelerated GFN-FF via CUDA (method name: "ggfnff")
 *
 * Claude Generated (March 2026): Replaces FFWorkspace CPU calculation with
 * FFWorkspaceGPU while keeping all CPU-side logic (EEQ charges, CN, topology)
 * in the existing GFNFF class unchanged.
 *
 * Architecture:
 *   - m_gfnff          : GFNFF instance (CPU topology + EEQ + H/X-bonds)
 *   - m_gpu_workspace  : FFWorkspaceGPU (dispersion, repulsion, Coulomb TERM 1,
 *                        bonds, angles, dihedrals, inversions on GPU)
 *   - H-bonds, X-bonds, ATM, BATM: remain in CPU GFNFF (Phase 1)
 *
 * Initialization flow:
 *   1. setMolecule() → m_gfnff->InitialiseMolecule() (topology, params)
 *   2. m_gfnff->generateGFNFFParameterSet() → FFWorkspaceGPU ctor (upload to GPU)
 *   3. m_gfnff->setGPUWorkspace(m_gpu_workspace.get()) (inject non-owning ptr)
 *
 * Per-step calculation:
 *   calculateEnergy() → m_gfnff->Calculation() which:
 *     a. Computes CN (CPU)
 *     b. Computes EEQ charges (CPU)
 *     c. Calls m_gpu_workspace->setGeometry/setD3CN/setEEQCharges/setCNDerivatives
 *     d. Dispatches to m_gpu_workspace->calculate() instead of CPU workspace
 */
class GGFNFFComputationalMethod : public ComputationalMethod {
public:
    explicit GGFNFFComputationalMethod(const std::string& method_name, const json& config);
    ~GGFNFFComputationalMethod() = default;

    // === ComputationalMethod interface ===

    bool setMolecule(const Mol& mol) override;
    double calculateEnergy(bool gradient = false) override;
    Matrix getGradient() const override;
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
    std::unique_ptr<FFWorkspace>    m_cpu_residual;  ///< CPU workspace for HB/XB/ATM/BATM (must outlive m_gpu_workspace)
    std::unique_ptr<FFWorkspaceGPU> m_gpu_workspace; ///< Destroyed first (holds raw ptr to m_cpu_residual)

    json             m_parameters;
    std::string      m_method_name;
    std::vector<int> m_atom_types;   ///< Element numbers (Z) — stored from setMolecule()
    bool             m_initialized   = false;
    bool             m_has_error     = false;
    std::string      m_error_message;
    double           m_last_energy   = 0.0;
};

#endif // USE_CUDA
