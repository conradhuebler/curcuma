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

    // ---- BACKEND-SPECIFIC 3: multipole potential of the host-driven loop --------------
    // Claude Generated (Oct 2026), -xtb.rocm_device_multipole (default false). The GFN2
    // multipole potential and energy (XTB::addMultipolePotential / energyMultipole) run on the
    // device from the host interaction matrices, uploaded once per setupMultipole(). Vulkan has
    // no such kernels, so this is not in the shared template.
    void setDeviceMultipole(bool on) { m_device_mp = on; }
    bool supportsDeviceMultipolePotential() const override { return m_device_mp; }

    /// Device limit for the 15 nat^2 doubles of the matrices (4 GB, nat up to about 5770).
    static constexpr double kMpMaxBytes = 4.0e9;

    bool beginMultipolePotential(const std::array<Eigen::MatrixXd, 3>& sd,
                                 const std::array<std::array<Eigen::MatrixXd, 3>, 3>& dd,
                                 const std::array<Eigen::MatrixXd, 6>& sq,
                                 const std::vector<double>& dkernel,
                                 const std::vector<double>& qkernel,
                                 unsigned long long version) override
    {
        if (!m_ctx || m_nat <= 0) return false;
        const Eigen::Index n = m_nat;
        auto okm = [n](const Eigen::MatrixXd& m) { return m.rows() == n && m.cols() == n; };
        const double* mats[15] = {
            sd[0].data(), sd[1].data(), sd[2].data(),
            dd[0][0].data(), dd[0][1].data(), dd[0][2].data(), dd[1][1].data(), dd[1][2].data(), dd[2][2].data(),
            sq[0].data(), sq[1].data(), sq[2].data(), sq[3].data(), sq[4].data(), sq[5].data() };
        bool shapes = static_cast<Eigen::Index>(dkernel.size()) == n
                   && static_cast<Eigen::Index>(qkernel.size()) == n;
        for (int k = 0; k < 3; ++k) shapes = shapes && okm(sd[k]);
        for (int a = 0; a < 3; ++a) for (int b = a; b < 3; ++b) shapes = shapes && okm(dd[a][b]);
        for (int k = 0; k < 6; ++k) shapes = shapes && okm(sq[k]);
        if (!shapes) { m_mp_error = "matrix shapes do not match the atom count"; return false; }
        if (!m_ctx->beginMultipolePotential(m_nat, mats, dkernel.data(), qkernel.data(), version,
                                            kMpMaxBytes)) {
            m_mp_error = m_ctx->multipolePotentialError();
            return false;
        }
        return true;
    }
    bool multipolePotential(const Vector& q_at, const Eigen::MatrixXd& dp_at,
                            const Eigen::MatrixXd& qp_at, Vector& v_at,
                            Eigen::MatrixXd& v_dp, Eigen::MatrixXd& v_qp) override
    {
        if (!m_ctx || m_nat <= 0 || q_at.size() != m_nat
            || dp_at.rows() != 3 || dp_at.cols() != m_nat
            || qp_at.rows() != 6 || qp_at.cols() != m_nat) return false;
        v_at.resize(m_nat);
        v_dp.resize(3, m_nat);
        v_qp.resize(6, m_nat);
        return m_ctx->multipolePotential(q_at.data(), dp_at.data(), qp_at.data(),
                                         v_at.data(), v_dp.data(), v_qp.data());
    }
    bool multipoleEnergy(const Vector& q_at, const Eigen::MatrixXd& dp_at,
                         const Eigen::MatrixXd& qp_at, double& e) override
    {
        if (!m_ctx || m_nat <= 0 || q_at.size() != m_nat
            || dp_at.rows() != 3 || dp_at.cols() != m_nat
            || qp_at.rows() != 6 || qp_at.cols() != m_nat) return false;
        return m_ctx->multipoleEnergy(q_at.data(), dp_at.data(), qp_at.data(), &e);
    }
    bool solveMultipoleResident(const Eigen::VectorXd& v_ao, Vector& eps, bool fp32) override
    {
        if (!m_ctx || static_cast<int>(v_ao.size()) != m_n) return false;
        eps.resize(m_n);
        return m_ctx->solveMultipoleResident(v_ao.data(), eps.data(), fp32);
    }
    void noteDeviceMultipoleFallback(const std::string& why) override
    {
        if (!m_mp_fallback.empty()) return;   // keep the first event
        m_mp_fallback = why + (m_mp_error.empty() ? std::string() : " (" + m_mp_error + ")");
    }
    /// First device-multipole fallback of this backend ("" = none). Claude Generated (Oct 2026).
    const std::string& deviceMultipoleFallback() const { return m_mp_fallback; }

private:
    bool        m_device_mp = false;
    std::string m_mp_error;      ///< reason of the last failed upload
    std::string m_mp_fallback;   ///< first fallback, read by reportSolverFallbacks
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
    bool device_multipole = false;   // -xtb.rocm_device_multipole (Claude Generated, Oct 2026)
    // Resident-solver switches -xtb.rocm_cached_cholesky (default false),
    // -xtb.rocm_density_syrk (default false) and -xtb.rocm_explicit_linv (default false); the fallbacks here are the defaults of the
    // PARAM block in native_xtb_method.h. Claude Generated (Oct 2026).
    {
        auto cfgBool = [&config](const char* key, bool def) -> bool {
            const json* x = nullptr;
            if (config.contains("xtb") && config["xtb"].is_object() && config["xtb"].contains(key))
                x = &config["xtb"][key];
            else if (config.contains(key))
                x = &config[key];
            if (!x) return def;
            if (x->is_boolean()) return x->get<bool>();
            if (x->is_number()) return x->get<double>() != 0.0;
            if (x->is_string()) {
                const std::string v = x->get<std::string>();
                if (v == "true" || v == "1" || v == "on") return true;
                if (v == "false" || v == "0" || v == "off") return false;
            }
            return def;
        };
        const bool cached_chol  = cfgBool("rocm_cached_cholesky", false);
        const bool density_syrk = cfgBool("rocm_density_syrk", false);
        const bool explicit_linv = cfgBool("rocm_explicit_linv", false);    // Claude Generated (Oct 2026), default false since Oct 17
        ctx->setSolverOptions(cached_chol, density_syrk, explicit_linv);
        // Claude Generated (Oct 2026): -xtb.rocm_device_multipole, default false (PARAM default
        // in native_xtb_method.h; this fallback is what takes effect without the key).
        device_multipole = cfgBool("rocm_device_multipole", false);
        if (CurcumaLogger::get_verbosity() >= 2)
            CurcumaLogger::info(fmt::format(
                "{}: ROCm solver: {} eigensolve, density via {}", getMethodName(),
                !cached_chol ? "sygvd"
                             : (explicit_linv ? "cached-Cholesky with explicit L^-1 (trmm/syevd/trmm)"
                                              : "cached-Cholesky (sygst/syevd/trsm)"),
                density_syrk ? "dsyrk" : "dgemm"));
        if (CurcumaLogger::get_verbosity() >= 2 && device_multipole && getMethodName() == "gfn2")
            CurcumaLogger::info(fmt::format(
                "{}: ROCm multipole potential and energy on the device (-xtb.rocm_device_multipole)",
                getMethodName()));
    }

    {
        auto backend = std::make_unique<HipScfBackend>(ctx);
        backend->setDeviceMultipole(device_multipole);
        m_scf_backend = std::move(backend);
    }
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
        // Claude Generated (Oct 2026): text corrected after a code reading (docs/ROCM_HOST_PARTS_2026-10-17.md):
        // the SCF loop runs on the host with a device eigensolver; potentials, occupations and the mixer are host work.
        const char* detail = getMethodName() == "gfn2"
            ? "integral build (CN/S/H0/L/gamma + dp/qp multipole), Fock, eigensolve, density and "
              "nuclear gradient incl. multipole-integral Pulay and D4 pair terms on the GPU; "
              "SCF potentials (Coulomb, third order, multipole, D4 weights), occupations, mixer, "
              "multipole interaction gradient and CN chain-rule on the CPU"
            : "integral build (CN/S/H0/L/gamma), Fock, eigensolve, density and nuclear gradient on the GPU; "
              "SCF potentials, occupations, mixer, D3 dispersion, halogen bond and CN chain-rule on the CPU";
        CurcumaLogger::info(fmt::format(
            "{}: ROCm device kernels (rocSOLVER + HIP) with a host SCF loop: {}",
            getMethodName(), detail));
    }
#else
    if (CurcumaLogger::get_verbosity() >= 2)
        CurcumaLogger::info(fmt::format(
            "{}: ROCm Stage 0 (device handshake); SCF/integrals/gradient on CPU "
            "(build without rocSOLVER)", getMethodName()));
#endif
}

// The destructor reports what an optimiser or MD driver kept silent: they run the energy
// calls at logger verbosity 0, so a fallback seen there is still pending here and is printed
// at level 1, as printGpuFallbackSummary() does for the GPU fallback summary.
// Claude Generated (Oct 2026).
XtbHipComputationalMethod::~XtbHipComputationalMethod() { reportSolverFallbacks(true); }

// Claude Generated (Oct 2026): one notice per kind of solver fallback and method object.
double XtbHipComputationalMethod::calculateEnergy(bool gradient)
{
    const double e = XtbGpuAdapter::calculateEnergy(gradient);
    reportSolverFallbacks();
    return e;
}

// Claude Generated (Oct 2026). The resident solver (xtb_hip_context.hip) keeps running on a
// slower device path when the explicit L^-1 or the cached Cholesky factor cannot be formed;
// the results stay the same up to rounding. The plugin context has no logger, so the events
// are read here. These are GPU-internal path changes, not CPU fallbacks, so they are not
// counted in the gpu_fallback registry (which -gpu_strict turns into an exit).
void XtbHipComputationalMethod::reportSolverFallbacks(bool final_report)
{
    XtbHipContext* ctx = context();
    if (!ctx) return;
#ifdef HAVE_ROCSOLVER
    // Claude Generated (Oct 2026): -xtb.rocm_device_multipole fell back to the host potential.
    // Same once-per-object rule as the solver notices below.
    if (!m_reported_mp_fallback && m_scf_backend) {
        const std::string& why = static_cast<HipScfBackend*>(m_scf_backend.get())->deviceMultipoleFallback();
        const int v0 = CurcumaLogger::get_verbosity();
        if (!why.empty() && (v0 >= 1 || final_report)) {
            if (v0 < 1) CurcumaLogger::set_verbosity(1);
            m_reported_mp_fallback = true;
            CurcumaLogger::warn(fmt::format(
                "{}: ROCm device multipole potential failed at the {}; the multipole potential and "
                "energy are computed on the host for the rest of the run (reported once)",
                getMethodName(), why));
            if (v0 < 1) CurcumaLogger::set_verbosity(v0);
        }
    }
#endif
    const unsigned ev = ctx->solverFallbackEvents();
    const unsigned fresh = ev & ~m_reported_fallbacks
        & (XtbHipContext::kFallbackLinvAlloc | XtbHipContext::kFallbackLinvTrtri
           | XtbHipContext::kFallbackCholFailed | XtbHipContext::kFallbackFactorRetried);
    if (!fresh) return;
    // At verbosity 0 the warning would be swallowed: keep the bits pending for a later call;
    // the final report (destructor) raises the level to 1 for these lines.
    const int saved_verbosity = CurcumaLogger::get_verbosity();
    if (saved_verbosity < 1 && !final_report) return;
    if (saved_verbosity < 1) CurcumaLogger::set_verbosity(1);
    m_reported_fallbacks |= fresh;
    const char* src = (ev & XtbHipContext::kFallbackInjected)
        ? " [forced by CURCUMA_ROCM_FAIL_LINV]" : "";
    if (fresh & XtbHipContext::kFallbackLinvAlloc)
        CurcumaLogger::warn(fmt::format(
            "{}: ROCm solver: the explicit L^-1 buffer could not be allocated; the SCF uses the "
            "cached-Cholesky path (sygst/syevd/trsm) for the rest of the run{} (reported once)",
            getMethodName(), src));
    if (fresh & XtbHipContext::kFallbackLinvTrtri)
        CurcumaLogger::warn(fmt::format(
            "{}: ROCm solver: trtri could not invert the Cholesky factor of S; this geometry uses "
            "the cached-Cholesky path (sygst/syevd/trsm), later geometries try L^-1 again{} "
            "(reported once)", getMethodName(), src));
    if (fresh & XtbHipContext::kFallbackCholFailed)
        CurcumaLogger::warn(fmt::format(
            "{}: ROCm solver: the cached Cholesky factor of S could not be formed (allocation "
            "failed or potrf reported S not positive definite); the affected geometries use "
            "rocsolver sygvd instead{} (reported once)", getMethodName(), src));
    if (fresh & XtbHipContext::kFallbackFactorRetried)
        CurcumaLogger::warn(fmt::format(
            "{}: ROCm solver: rocSOLVER returned a Cholesky factor or inverse that failed the "
            "L L^T v = S v probe and it was recomputed{} (reported once)", getMethodName(), src));
    if (saved_verbosity < 1) CurcumaLogger::set_verbosity(saved_verbosity);
}

#endif // USE_ROCM
