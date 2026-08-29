/*
 * <Native xTB GPU Method Wrapper — GFN1 / GFN2 on CUDA>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0.
 *
 * Claude Generated (2026-06): ComputationalMethod adapter for the GPU path of
 * native GFN1/GFN2. Selected by the factory when `-gpu cuda` (or `-gpu auto`)
 * is given and the build has USE_CUDA.
 *
 * Design (minimal-disruption sibling of NativeXtbMethod):
 *   - OWNS a NativeXtbMethod (the validated CPU pipeline: config, large-system
 *     modes, error handling, property accessors). Not a subclass of XTB.
 *   - OWNS an XtbGpuContext (cuSOLVER/cuBLAS handles + stream).
 *
 * Stage 0 (current): every call forwards to the CPU NativeXtbMethod; the GPU
 * context is only initialized to prove the build/link and device handshake.
 * Later stages swap the hot numerics (eigensolve, SCF, integrals, gradient)
 * onto the device while keeping host-side setup in the owned NativeXtbMethod.
 *
 * Usage:
 *   ./curcuma -sp mol.xyz -method gfn2 -gpu cuda    # explicit GPU
 *   ./curcuma -sp mol.xyz -method gfn1 -gpu auto    # GPU if available
 *   ./curcuma -sp mol.xyz -method gfn2              # CPU (default)
 */

#pragma once

#ifdef USE_CUDA

#include "native_xtb_method.h"   // NativeXtbMethod + curcuma::xtb::MethodType

#include <memory>
#include <string>

namespace curcuma {
namespace xtb {
namespace gpu {
class XtbGpuContext;
}
}
}

/**
 * @brief GPU-accelerated native GFN1/GFN2 (method "gfn1"/"gfn2" with -gpu cuda).
 * Claude Generated.
 */
class XtbGpuComputationalMethod : public ComputationalMethod {
public:
    explicit XtbGpuComputationalMethod(curcuma::xtb::MethodType method,
                                       const json& config = json{});
    ~XtbGpuComputationalMethod() override;

    // ---- ComputationalMethod interface (Stage 0: forward to CPU) ----------
    bool setMolecule(const Mol& mol) override;
    bool updateGeometry(const Matrix& geometry) override;
    double calculateEnergy(bool gradient = false) override;

    Matrix getGradient() const override;
    Vector getCharges() const override;
    Vector getBondOrders() const override;
    Position getDipole() const override;
    bool hasGradient() const override;

    std::string getMethodName() const override;
    bool isThreadSafe() const override;
    void setThreadCount(int threads) override;

    void setParameters(const json& params) override;
    json getParameters() const override;

    bool hasError() const override;
    void clearError() override;
    std::string getErrorMessage() const override;

    Vector getOrbitalEnergies() const override;
    int getNumElectrons() const override;
    json getEnergyDecomposition() const override;
    bool saveToFile(const std::string& filename) const override;

    void setWarmStart(bool on) override;
    void setIterativeMode(bool on) override;

    /// True when the GPU context is live (else this object runs the CPU path).
    bool gpuActive() const;

private:
    curcuma::xtb::MethodType                          m_method;
    // Declared before m_cpu so it is destroyed AFTER m_cpu: m_cpu owns the XTB
    // that holds the eigensolver hook + resident-SCF backend pointer capturing
    // this context / the backend.
    std::unique_ptr<curcuma::xtb::gpu::XtbGpuContext> m_gpu;         ///< device handles
    std::unique_ptr<curcuma::xtb::GpuScfBackend>      m_scf_backend; ///< Stage-2 resident SCF (GFN1)
    std::unique_ptr<NativeXtbMethod>                  m_cpu; ///< validated CPU pipeline + hook holder
};

/// Build a device-resident SCF+gradient backend over an existing context
/// (Stage 4 validation, Claude Generated). Lets a bare XTB be driven on the GPU
/// path via XTB::setGpuScfBackend without the full XtbGpuComputationalMethod.
std::unique_ptr<curcuma::xtb::GpuScfBackend>
createXtbGpuScfBackend(curcuma::xtb::gpu::XtbGpuContext* ctx);

#endif // USE_CUDA
