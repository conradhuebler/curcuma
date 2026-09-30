/*
 * <Native xTB Vulkan Method Wrapper — implementation>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0.
 *
 * Claude Generated (2026-06, restructured Sep 2026): everything the three GPU wrappers
 * share — the ComputationalMethod forwarding, the device handshake/logging/CPU fallback,
 * and the GpuScfBackend marshalling — now lives in xtb_gpu_adapter.h. This file keeps
 * ONLY what is specific to Vulkan:
 *   - the eigensolver hook (host Cholesky reduction + device Jacobi + back-transform;
 *     ROCm hands the generalized problem straight to rocSOLVER, CUDA to cuSOLVER),
 *   - three GpuScfBackend overrides (host Cholesky, the GFN2 gradient escape hatch),
 *   - the mixed-precision policy (OFF for Vulkan — measured net-negative),
 *   - the stage banner.
 * Host-compiled (no Vulkan headers — the context is pimpl).
 */

#ifdef USE_VULKAN

#include "xtb_vulkan_method.h"

#include <cstdlib>

namespace {
using curcuma::xtb::MethodType;
using curcuma::xtb::gpu::XtbVulkanContext;

/**
 * @brief Device-resident GFN1/GFN2 SCF backend over an XtbVulkanContext.
 *
 * All ~40 GpuScfBackend virtuals are forwarded by the shared XtbGpuResidentBackend
 * template (xtb_gpu_adapter.h). Only the two genuinely Vulkan-specific pieces are here.
 * Claude Generated (Stage 2a/V-AP2/V-AP3; restructured Sep 2026).
 */
class VulkanScfBackend
    : public curcuma::xtb::gpu::XtbGpuResidentBackend<XtbVulkanContext,
                                                     curcuma::xtb::gpu::XtbVulkanBasisData> {
public:
    explicit VulkanScfBackend(XtbVulkanContext* ctx)
        : XtbGpuResidentBackend(ctx) {}

    // ---- BACKEND-SPECIFIC 1: the Cholesky factor -------------------------
    // The device builds S but has no Cholesky shader; L = chol(S) is cheap host-side
    // (Eigen LLT) and matches the host m_X orthonormalizer. ROCm downloads a
    // device-built L instead, so this cannot live in the shared template.
    bool downloadCholesky(Eigen::MatrixXd& L_out) override
    {
        if (!m_ctx || m_n <= 0) return false;
        Eigen::MatrixXd S(m_n, m_n);
        if (!m_ctx->downloadOverlap(S.data())) return false;
        Eigen::LLT<Eigen::MatrixXd> llt(S);
        if (llt.info() != Eigen::Success) return false;
        L_out = llt.matrixL();
        return true;
    }

protected:
    // ---- BACKEND-SPECIFIC 2: the GFN2 gradient escape hatch ---------------
    // Both GFN1 and GFN2 gradients run on the device by default — the workgroup-per-atom
    // grad_pulay (V-PERF-2) made it fast (231-atom complex: GFN2 28 s -> 81 ms, GFN1
    // 18 s -> 98 ms), so it beats the host gradient. CURCUMA_VK_GFN2_CPUGRAD=1 forces the
    // validated host gradient for GFN2 (a Vulkan-only debug hatch). Claude Generated.
    bool deviceGradientEnabled(bool is_gfn2) const override
    {
        static const bool gfn2_cpu_grad = std::getenv("CURCUMA_VK_GFN2_CPUGRAD") != nullptr;
        return !(is_gfn2 && gfn2_cpu_grad);
    }
};
}  // namespace

XtbVulkanComputationalMethod::XtbVulkanComputationalMethod(MethodType method, const json& config)
    : XtbGpuAdapter(method, config, "Vulkan", "Vulkan compute device (FP64 required)")
{
    if (!gpuActive()) return;   // the adapter already warned; the CPU path stands

    curcuma::xtb::XTB* xtb = cpuSolver();
    if (!xtb) return;
    XtbVulkanContext* ctx = context();

    // Stage 1: install the GPU eigensolver. solveEigen() delegates the per-iteration
    // generalized eigenproblem (F, S=L·Lᵀ)→(C, eps) to the device; everything else
    // (integrals, Fock, density, gradient) stays on the CPU pipeline. VULKAN-SPECIFIC:
    // there is no LAPACK on Vulkan, so the generalized problem is reduced to a STANDARD
    // symmetric one on the host (two cheap triangular solves with the Cholesky factor L)
    // and only the dense symmetric eigensolve — the SCF hot path — runs on the GPU
    // (FP64 two-sided Jacobi, validated vs Eigen). ROCm/CUDA instead hand the generalized
    // pair straight to rocSOLVER/cuSOLVER.
    xtb->setExternalEigensolver(
        [ctx](const Matrix& F, const Eigen::MatrixXd& L,
              Matrix& C, Vector& eps) -> bool {
            const int n = static_cast<int>(F.rows());
            if (n <= 0 || L.rows() != n || L.cols() != n) return false;
            // Reduce F C = S C ε (S = L·Lᵀ) to the standard problem
            // Ã = L⁻¹ F L⁻ᵀ on the host (cheap triangular solves).
            Eigen::MatrixXd Fcm = F;  // → column-major (F symmetric)
            const Eigen::MatrixXd Y  = L.triangularView<Eigen::Lower>().solve(Fcm);          // L Y = F
            const Eigen::MatrixXd W  = L.triangularView<Eigen::Lower>().solve(Y.transpose()); // W = Ã (symmetric)
            Eigen::MatrixXd Atil = 0.5 * (W + W.transpose());                                 // symmetrize
            // GPU: standard symmetric eigensolve of Ã.
            eps.resize(n);
            Eigen::MatrixXd Ctil(n, n);
            if (!ctx->solveSymmetric(Atil.data(), n, eps.data(), Ctil.data())) return false;
            // Back-transform the generalized eigenvectors: C = L⁻ᵀ C̃ (Lᵀ C = C̃).
            C = L.transpose().triangularView<Eigen::Upper>().solve(Ctil);
            return true;
        });

    // Stage 2: install the device-resident SCF backend. Under the default Broyden mixing,
    // XTB::Calculation keeps H0/S (and the per-iteration density and MO coefficients)
    // resident on the device for the whole SCF — only length-nao vectors cross the bus per
    // iteration. GFN1 runs the isotropic loop, GFN2 the multipole loop (V-AP3).
    m_scf_backend = std::make_unique<VulkanScfBackend>(ctx);
    xtb->setGpuScfBackend(m_scf_backend.get());

    // VULKAN-SPECIFIC: mixed precision is deliberately NOT defaulted on (unlike
    // CUDA/ROCm). The FP32 Jacobi path exists (opt-in -scf_mixed_precision) and is
    // correct, but measured net-negative on the iGPU: the hand-written two-sided cyclic
    // Jacobi is dispatch/barrier-bound, not FP64-arithmetic-bound, so FP32 is ~equal per
    // iteration (complex/231: FP32 828 vs FP64 844 ms/iter) and the perturbed early
    // iterations occasionally cost an extra SCF cycle (GFN1 14->15), making it slower
    // overall. The real Vulkan lever is the eigensolve algorithm, not its precision.
    // See docs/SQM_GPU_ROADMAP.md X-AP3.
    //
    // This used to rely on the XTB member defaulting to false. That default is now TRUE
    // (it pays on the CPU/MKL eigensolve), so Vulkan must opt OUT explicitly.
    // applyXtbScfConfig already ran (inside the NativeXtbMethod constructor), so this call
    // would clobber a user-supplied flag — hence it only fires when the user did NOT ask
    // for mixed precision. Claude Generated.
    const bool mp_user_set
        = config.contains("scf_mixed_precision")
        || (config.contains("xtb") && config["xtb"].is_object()
            && config["xtb"].contains("scf_mixed_precision"));
    if (!mp_user_set)
        xtb->setMixedPrecision(false);

    if (CurcumaLogger::get_verbosity() >= 2) {
        const bool gfn1 = getMethodName() == "gfn1";
        CurcumaLogger::info(fmt::format(
            "{}: Vulkan on-device integral build (CN/S/H0/L/gamma) + GPU FP64 "
            "Jacobi/Lowdin eigensolve; {}",
            getMethodName(),
            gfn1 ? "device-resident isotropic SCF loop + nuclear gradient (Stage 4; "
                   "dispersion + CN chain-rule on CPU)"
                 : "device-resident multipole SCF loop (Stage 2b; Fock+moments on GPU, "
                   "gradient on CPU)"));
    }
}

XtbVulkanComputationalMethod::~XtbVulkanComputationalMethod() = default;

#endif // USE_VULKAN
