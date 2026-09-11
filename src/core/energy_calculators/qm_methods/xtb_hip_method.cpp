/*
 * <Native xTB ROCm/HIP Method Wrapper — implementation>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0.
 *
 * Claude Generated (2026-06, restructured Sep 2026): everything the three GPU wrappers
 * share — the ComputationalMethod forwarding, the device handshake/logging/CPU fallback,
 * and the GpuScfBackend marshalling — now lives in xtb_gpu_adapter.h. This file keeps
 * ONLY what is specific to ROCm:
 *   - the eigensolver hook (rocSOLVER dsygvd solves the GENERALIZED problem directly, so
 *     the host rebuilds S = L·Lᵀ; Vulkan instead reduces to a standard problem on the host),
 *   - two GpuScfBackend overrides (the device-built Cholesky download, the d-shell flag),
 *   - the mixed-precision policy (ON for ROCm),
 *   - the stage banner + the no-rocSOLVER Stage-0 note.
 * Host-compiled (no device code) — the HIP kernels live in the .hip TU.
 */

#ifdef USE_ROCM

#include "xtb_hip_method.h"

namespace {
using curcuma::xtb::MethodType;
using curcuma::xtb::gpu::XtbHipContext;

#ifdef HAVE_ROCSOLVER
/**
 * @brief Device-resident GFN1/GFN2 SCF backend over an XtbHipContext.
 *
 * All ~40 GpuScfBackend virtuals are forwarded by the shared XtbGpuResidentBackend
 * template (xtb_gpu_adapter.h). Only the two genuinely ROCm-specific pieces are here.
 * Claude Generated (Stage 2/R-AP1/R-AP2/R-AP3; restructured Sep 2026).
 */
class HipScfBackend
    : public curcuma::xtb::gpu::XtbGpuResidentBackend<XtbHipContext,
                                                     curcuma::xtb::gpu::XtbHipBasisData> {
public:
    explicit HipScfBackend(XtbHipContext* ctx)
        : XtbGpuResidentBackend(ctx) {}

    // ---- BACKEND-SPECIFIC 1: the Cholesky factor -------------------------
    // ROCm builds L = chol(S) on the device (rocSOLVER) and downloads it. Vulkan has no
    // Cholesky shader and computes it host-side from the device S, so this cannot live
    // in the shared template.
    bool downloadCholesky(Eigen::MatrixXd& L_out) override
    {
        if (!m_ctx || m_n <= 0) return false;
        L_out.resize(m_n, m_n);
        return m_ctx->downloadCholesky(L_out.data());
    }

    // ---- BACKEND-SPECIFIC 2: d shells ------------------------------------
    // X-I1 B6: the HIP integral/SCF/gradient kernels handle d shells (cartesian->spherical
    // dtrafo in rocm/xtb_hip_integrals.hiph + the dpair branch in k_overlap_h0 /
    // k_multipole_ints / k_grad_h0_pulay). The Vulkan GLSL shaders do not yet, so they
    // keep the base false (d systems route to the CPU path). Claude Generated.
    bool supportsDshell() const override { return true; }
};
#endif // HAVE_ROCSOLVER
}  // namespace

XtbHipComputationalMethod::XtbHipComputationalMethod(MethodType method, const json& config)
    : XtbGpuAdapter(method, config, "ROCm", "ROCm/HIP device")
{
    if (!gpuActive()) return;   // the adapter already warned; the CPU path stands

#ifdef HAVE_ROCSOLVER
    curcuma::xtb::XTB* xtb = cpuSolver();
    if (!xtb) return;
    XtbHipContext* ctx = context();

    // Stage 1: install the GPU eigensolver. solveEigen() delegates the per-iteration
    // generalized eigenproblem (F, S = L·Lᵀ) → (C, eps) to rocSOLVER (dsygvd) on the
    // device; everything else (integrals, Fock, density, gradient) stays on the CPU.
    // ROCm-SPECIFIC: rocSOLVER solves the GENERALIZED problem directly, so the host
    // reconstructs S from the Cholesky factor and hands both matrices over — no host
    // reduction/back-transform (that is the Vulkan route). Works for GFN1 and GFN2 alike
    // (the host hands it the complete Fock).
    xtb->setExternalEigensolver(
        [ctx](const Matrix& F, const Eigen::MatrixXd& L,
              Matrix& C, Vector& eps) -> bool {
            const int n = static_cast<int>(F.rows());
            if (n <= 0 || L.rows() != n || L.cols() != n) return false;
            Eigen::MatrixXd Fcm = F;                                   // column-major
            Eigen::MatrixXd Ll  = L.triangularView<Eigen::Lower>();    // lower Cholesky
            Eigen::MatrixXd Scm = Ll * Ll.transpose();                 // S = L·Lᵀ
            eps.resize(n);
            Eigen::MatrixXd Ccm(n, n);
            if (!ctx->solveGeneralized(Fcm.data(), Scm.data(), n, eps.data(), Ccm.data()))
                return false;
            C = Ccm;
            return true;
        });

    // Stage 2: install the device-resident SCF backend. Under the default Broyden mixing,
    // XTB::Calculation keeps H0/S (and the per-iteration density and MO coefficients)
    // resident on the device — the Fock build + populations are HIP kernels, the density a
    // rocBLAS GEMM, the eigensolve rocSOLVER. GFN1 runs the isotropic loop, GFN2 the
    // multipole loop (R-AP2).
    m_scf_backend = std::make_unique<HipScfBackend>(ctx);
    xtb->setGpuScfBackend(m_scf_backend.get());

    // ROCm-SPECIFIC: FP64 is ~1/16 of FP32 on this iGPU and the per-iteration eigensolve
    // dominates the resident SCF, so default mixed precision ON for the GPU path:
    // far-from-convergence iterations solve in FP32 (rocsolver_ssygvd), reverting to FP64
    // once max|dq| < scf_fp32_threshold so the converged fixed point and energy stay FP64.
    // Matches the CUDA backend (X-AP3); Vulkan deliberately opts out.
    xtb->setMixedPrecision(true);

    if (CurcumaLogger::get_verbosity() >= 2) {
        // The format string must be a compile-time constant (fmt >= 9 checks it in a
        // consteval context), so the method-dependent part is selected as plain text.
        const char* detail = getMethodName() == "gfn2"
            ? "integral build (CN/S/H0/L/gamma + dp/qp multipole) + multipole SCF + "
              "nuclear gradient incl. multipole-integral Pulay on the GPU (Stage 4 / "
              "R-AP3); only the multipole SD/DD/SQ interaction gradient + dispersion + "
              "CN chain-rule on CPU"
            : "integral build (CN/S/H0/L/gamma) + SCF + nuclear gradient on the GPU; "
              "only the dispersion gradient + CN chain-rule on CPU (Stage 4)";
        CurcumaLogger::info(fmt::format(
            "{}: ROCm fully device-resident (rocSOLVER + HIP kernels): {}",
            getMethodName(), detail));
    }
#else
    if (CurcumaLogger::get_verbosity() >= 2)
        CurcumaLogger::info(fmt::format(
            "{}: ROCm Stage 0 (device handshake); SCF/integrals/gradient on CPU "
            "(build without rocSOLVER)", getMethodName()));
#endif
}

XtbHipComputationalMethod::~XtbHipComputationalMethod() = default;

#endif // USE_ROCM
