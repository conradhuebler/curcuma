/*
 * <Native xTB Vulkan Method Wrapper — GFN1 / GFN2 on Vulkan>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0.
 *
 * Claude Generated (2026-06): ComputationalMethod adapter for the Vulkan path of
 * native GFN1/GFN2 — the Vulkan sibling of XtbGpuComputationalMethod (CUDA) /
 * XtbHipComputationalMethod (ROCm). Selected by the factory when `-gpu vulkan` (or
 * `-gpu auto` on a Vulkan-only build) is given and the build has USE_VULKAN_XTB.
 *
 * Design (identical to the CUDA/ROCm wrappers): OWNS a NativeXtbMethod (validated CPU
 * pipeline) + an XtbVulkanContext (device handles). Stage 0 forwards every call to the
 * CPU; later stages install the GPU eigensolver hook + device-resident SCF backend.
 *
 * Usage:
 *   ./curcuma -sp mol.xyz -method gfn2 -gpu vulkan   # explicit Vulkan
 *   ./curcuma -sp mol.xyz -method gfn1 -gpu auto     # GPU if available
 *   ./curcuma -sp mol.xyz -method gfn2               # CPU (default)
 */

#pragma once

#ifdef USE_VULKAN_XTB

#include "native_xtb_method.h"   // NativeXtbMethod + curcuma::xtb::MethodType

#include <memory>
#include <string>

namespace curcuma {
namespace xtb {
struct GpuScfBackend;
namespace gpu {
class XtbVulkanContext;
}
}
}

/**
 * @brief GPU-accelerated native GFN1/GFN2 (method "gfn1"/"gfn2" with -gpu vulkan).
 * Claude Generated.
 */
class XtbVulkanComputationalMethod : public ComputationalMethod {
public:
    explicit XtbVulkanComputationalMethod(curcuma::xtb::MethodType method,
                                          const json& config = json{});
    ~XtbVulkanComputationalMethod() override;

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

    /// True when the Vulkan context is live (else this object runs the CPU path).
    bool gpuActive() const;

private:
    curcuma::xtb::MethodType                             m_method;
    // m_gpu/m_scf_backend before m_cpu: the owned XTB holds the eigensolver hook +
    // resident-SCF backend pointer, so it must be destroyed (in m_cpu) first.
    std::unique_ptr<curcuma::xtb::gpu::XtbVulkanContext> m_gpu;          ///< device handles
    std::unique_ptr<curcuma::xtb::GpuScfBackend>         m_scf_backend;  ///< Stage-2 resident SCF (GFN1)
    std::unique_ptr<NativeXtbMethod>                     m_cpu;          ///< validated CPU pipeline
};

#endif // USE_VULKAN_XTB
