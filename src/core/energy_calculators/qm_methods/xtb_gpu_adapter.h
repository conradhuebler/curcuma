/*
 * <Native xTB GPU adapter templates — shared CUDA / ROCm / Vulkan wrapper logic>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0.
 *
 * Claude Generated (Sep 2026): the three GPU wrappers (xtb_gpu_method.cpp = CUDA,
 * xtb_hip_method.cpp = ROCm, xtb_vulkan_method.cpp = Vulkan) were three copies of the
 * same two classes. This header holds the ONE copy of each; the per-backend files keep
 * only what is genuinely backend-specific.
 *
 *   1. XtbGpuAdapter<Context>        — the ComputationalMethod half. Owns the device
 *      context + the validated CPU NativeXtbMethod, does the device handshake and its
 *      logging, and forwards every ComputationalMethod call to the CPU pipeline
 *      (the GPU is reached through the hooks the derived constructor installs on the
 *      owned XTB, never through this interface). Used by ALL THREE backends.
 *
 *   2. XtbGpuResidentBackend<Context, BasisData> — the GpuScfBackend half: the ~40
 *      virtuals of the xtb_native.h seam forwarded to a device context whose methods
 *      take raw column-major pointers. Used by ROCm + Vulkan, whose contexts expose the
 *      SAME method set (residentBegin/Solve/Density/Finalize, beginBasis/beginComputed,
 *      download*, gradient, multipole, dispersion, EEQ). CUDA is deliberately NOT on
 *      this template: its context carries extra device-resident stages (the fused
 *      resident SCF loop, the device potential/solvation build, atomic/shell charges,
 *      the SCC energy) and passes `n` explicitly to most calls, so its backend class
 *      stays hand-written in xtb_gpu_method.cpp — forcing it in would mean a dozen
 *      more trait hooks for one user each. CUDA does share template 1.
 *
 * Educational note: the two templates are plain class templates with virtual methods —
 * a derived backend overrides the handful of virtuals where the hardware really differs
 * (see the "backend-specific" comments in the per-backend .cpp files). No CRTP, no
 * policy soup: the only template parameters are the context type and its basis-data
 * struct, because those types differ per backend but the CALLS on them do not.
 */

#pragma once

#include "native_xtb_method.h"   // NativeXtbMethod + curcuma::xtb::MethodType + XTB seam

#include "src/core/curcuma_logger.h"

#include <fmt/format.h>

#include <array>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

namespace curcuma {
namespace xtb {
namespace gpu {

/* ===================================================================== *
 *  1. XtbGpuAdapter — the shared ComputationalMethod wrapper.
 * ===================================================================== */

/**
 * @brief GPU-backed native GFN1/GFN2 ComputationalMethod, shared by all backends.
 *
 * Construction: create the device context (non-throwing; ok() is false when no usable
 * device), then the full validated CPU pipeline (NativeXtbMethod: config, large-system
 * modes, errors, properties). On a live device it logs the handshake; otherwise it warns
 * and the object simply runs the CPU path — the fallback that makes every GPU stage safe.
 *
 * The derived per-backend constructor then installs its own device hooks on
 * `cpuSolver()` (the external eigensolver, the resident GpuScfBackend, the mixed-
 * precision policy) and logs its stage banner. Everything else — all 20 forwarding
 * methods and gpuActive() — lives here once.
 *
 * @tparam Context device context class; needs Context(int device), ok() / deviceId() /
 *                deviceName() / bindDevice().
 *
 * Claude Generated (Sep 2026).
 */
template <class Context>
class XtbGpuAdapter : public ComputationalMethod {
public:
    /**
     * @param method       GFN1 or GFN2
     * @param config       method config (forwarded verbatim to NativeXtbMethod)
     * @param ctx_label    backend name in the ready message ("GPU" / "ROCm" / "Vulkan")
     * @param device_label device wording in the fallback warning ("CUDA device", ...)
     */
    XtbGpuAdapter(curcuma::xtb::MethodType method, const json& config,
                  const char* ctx_label, const char* device_label)
        : m_method(method)
    {
        // Device handshake (non-throwing; ok() is false when no usable device).
        // Claude Generated (Sep 2026, multi-GPU): `gpu_device` (global CLI key, set per
        // worker by the batch capabilities) selects the device; -1 = the backend default.
        m_device = gpuDeviceFromConfig(config);
        m_gpu = std::make_unique<Context>(m_device);

        // Full validated CPU pipeline (config, large-system modes, errors, properties).
        m_cpu = std::make_unique<NativeXtbMethod>(method, config);

        // m_cpu->getMethodName() rather than the virtual getMethodName(): this runs in a
        // BASE constructor, where virtual dispatch would not reach a derived override.
        if (m_gpu->ok()) {
            if (CurcumaLogger::get_verbosity() >= 1)
                CurcumaLogger::success(fmt::format(
                    "{}: {} context ready on device {} ({})",
                    m_cpu->getMethodName(), ctx_label, m_gpu->deviceId(), m_gpu->deviceName()));
        } else if (m_device >= 0) {
            CurcumaLogger::warn(fmt::format(
                "{}: {} {} is not usable (index out of range or init failed); running CPU path",
                m_cpu->getMethodName(), device_label, m_device));
        } else {
            CurcumaLogger::warn(fmt::format(
                "{}: no usable {}; running CPU path", m_cpu->getMethodName(), device_label));
        }
    }

    // Bind before the members (context, resident backend, CPU solver holding the hooks)
    // are destroyed, so device memory is freed on the device that owns it.
    ~XtbGpuAdapter() override { bindDevice(); }

    // ---- ComputationalMethod interface: forward to the CPU pipeline -------
    // The three calls that reach the device re-bind first: the object may be driven from a
    // different host thread than the one that built it (the CUDA current device is per thread).
    bool setMolecule(const Mol& mol) override { bindDevice(); return m_cpu->setMolecule(mol); }
    bool updateGeometry(const Matrix& g) override { bindDevice(); return m_cpu->updateGeometry(g); }
    double calculateEnergy(bool gradient = false) override { bindDevice(); return m_cpu->calculateEnergy(gradient); }

    Matrix getGradient() const override { return m_cpu->getGradient(); }
    Vector getCharges() const override { return m_cpu->getCharges(); }
    Vector getBondOrders() const override { return m_cpu->getBondOrders(); }
    Position getDipole() const override { return m_cpu->getDipole(); }
    bool hasGradient() const override { return m_cpu->hasGradient(); }

    std::string getMethodName() const override { return m_cpu->getMethodName(); }
    bool isThreadSafe() const override { return m_cpu->isThreadSafe(); }
    void setThreadCount(int threads) override { m_cpu->setThreadCount(threads); }

    void setParameters(const json& params) override { m_cpu->setParameters(params); }
    json getParameters() const override { return m_cpu->getParameters(); }

    bool hasError() const override { return m_cpu->hasError(); }
    void clearError() override { m_cpu->clearError(); }
    std::string getErrorMessage() const override { return m_cpu->getErrorMessage(); }

    Vector getOrbitalEnergies() const override { return m_cpu->getOrbitalEnergies(); }
    int getNumElectrons() const override { return m_cpu->getNumElectrons(); }
    json getEnergyDecomposition() const override { return m_cpu->getEnergyDecomposition(); }
    bool saveToFile(const std::string& f) const override { return m_cpu->saveToFile(f); }

    void setWarmStart(bool on) override { m_cpu->setWarmStart(on); }
    void setIterativeMode(bool on) override { m_cpu->setIterativeMode(on); }

    /// True when the device context is live (else this object runs the CPU path).
    bool gpuActive() const { return m_gpu && m_gpu->ok(); }

    /// Device index this object was asked to use (-1 = backend default). Claude Generated.
    int gpuDevice() const { return m_device; }

    /**
     * @brief Read the `gpu_device` key (int, or a numeric string) from a method config.
     * @return the index, or -1 when absent/invalid (= backend default device).
     * Claude Generated (Sep 2026, multi-GPU).
     */
    static int gpuDeviceFromConfig(const json& config)
    {
        if (!config.contains("gpu_device")) return -1;
        const auto& v = config["gpu_device"];
        try {
            // CLI2Json stores numeric flags as double ("-gpu_device 3" -> 3.0).
            if (v.is_number()) return static_cast<int>(v.get<double>());
            if (v.is_string() && !v.get<std::string>().empty()) return std::stoi(v.get<std::string>());
        } catch (...) {
        }
        return -1;
    }

protected:
    void bindDevice() const { if (m_gpu && m_gpu->ok()) m_gpu->bindDevice(); }

    /// The owned XTB the derived constructor hangs its device hooks on (may be null).
    curcuma::xtb::XTB* cpuSolver() { return m_cpu->solver(); }
    Context* context() { return m_gpu.get(); }

    curcuma::xtb::MethodType m_method;
    int                      m_device = -1;   ///< requested device index (-1 = default)
    // m_gpu/m_scf_backend before m_cpu: the owned XTB holds the eigensolver hook + the
    // resident-SCF backend pointer, so it must be destroyed (in m_cpu) FIRST — members
    // are destroyed in reverse declaration order.
    std::unique_ptr<Context>                     m_gpu;          ///< device handles
    std::unique_ptr<curcuma::xtb::GpuScfBackend> m_scf_backend;  ///< device-resident SCF
    std::unique_ptr<NativeXtbMethod>             m_cpu;          ///< validated CPU pipeline
};

/* ===================================================================== *
 *  2. XtbGpuResidentBackend — the shared GpuScfBackend forwarding.
 * ===================================================================== */

/**
 * @brief Device-resident SCF/integral/gradient backend over a pointer-API context.
 *
 * Adapts the project-typed GpuScfBackend seam (xtb_native.h) that the core XTB calls to
 * a device context whose methods take raw column-major `double*`. All the marshalling —
 * the row-major project `Matrix` → column-major copies (H0/S/P are symmetric, so the copy
 * is value-preserving), the flattened `GpuBasisFlat`/`GpuH0Flat` → the backend's
 * `BasisData` struct, the 3·nao²/6·nao² multipole unpacking, the (nat,3) gradient
 * transposition — is written here once.
 *
 * Any method returning false makes the core fall back to the CPU for that whole
 * calculation, so an unavailable stage never corrupts results.
 *
 * @tparam Context   device context (ROCm XtbHipContext / Vulkan XtbVulkanContext)
 * @tparam BasisData the context's flat basis-upload struct
 *
 * Claude Generated (Sep 2026).
 */
template <class Context, class BasisData>
class XtbGpuResidentBackend : public curcuma::xtb::GpuScfBackend {
public:
    explicit XtbGpuResidentBackend(Context* ctx) : m_ctx(ctx) {}

    // ---- Stage 2: device-resident SCF -------------------------------------
    bool begin(const Matrix& H0, const Matrix& S, const Eigen::MatrixXd& L) override
    {
        const int n = static_cast<int>(H0.rows());
        if (!m_ctx || n <= 0 || S.rows() != n) return false;
        // L unused: ROCm hands S straight to rocSOLVER dsygvd, Vulkan builds the Löwdin
        // S^-1/2 on-device from S — neither needs the host Cholesky factor.
        (void)L;
        m_n = n;
        Eigen::MatrixXd Hcm = H0, Scm = S;  // column-major copies (H0/S symmetric)
        return m_ctx->residentBegin(Hcm.data(), Scm.data(), n);
    }
    bool solve(const Eigen::VectorXd& v_ao, Vector& eps, bool fp32 = false, int n_eig = 0) override
    {
        (void)n_eig;  // full spectrum; fp32 → the backend's FP32 solver (X-AP3)
        if (!m_ctx || static_cast<int>(v_ao.size()) != m_n) return false;
        eps.resize(m_n);
        return m_ctx->residentSolve(v_ao.data(), eps.data(), fp32);
    }
    bool density(const Eigen::VectorXd& occ, int ncol, Eigen::VectorXd& pop_ao, double& band) override
    {
        if (!m_ctx || ncol < 0 || ncol > m_n) return false;
        pop_ao.resize(m_n); band = 0.0;
        return m_ctx->residentDensity(occ.data(), ncol, pop_ao.data(), &band);
    }
    bool finalize(Matrix& P, Matrix& C) override
    {
        if (!m_ctx || m_n <= 0) return false;
        Eigen::MatrixXd Pcm(m_n, m_n), Ccm(m_n, m_n);
        if (!m_ctx->residentFinalize(Pcm.data(), Ccm.data())) return false;
        P = Pcm; C = Ccm;
        return true;
    }

    // ---- Stage 3: device integral build -----------------------------------
    bool beginBasis(const curcuma::xtb::GpuBasisFlat& bf,
                    const curcuma::xtb::GpuH0Flat& hf) override
    {
        if (!m_ctx || bf.nao <= 0) return false;
        m_n = bf.nao; m_nsh = bf.nsh; m_nat = bf.nat;
        BasisData bd;
        bd.nat = bf.nat; bd.nsh = bf.nsh; bd.nao = bf.nao; bd.is_gfn2 = bf.is_gfn2;
        bd.nprim_total = static_cast<int>(bf.prim_alpha.size());
        bd.z           = bf.z.data();
        bd.sh2at       = bf.sh2at.data();
        bd.ang_sh      = bf.ang_sh.data();
        bd.iao_sh      = bf.iao_sh.data();
        bd.nao_sh      = bf.nao_sh.data();
        bd.sh_nprim    = bf.sh_nprim.data();
        bd.sh_prim_off = bf.sh_prim_off.data();
        bd.prim_alpha  = bf.prim_alpha.data();
        bd.prim_coeff  = bf.prim_coeff.data();
        bd.sh_zeta     = bf.sh_zeta.data();
        bd.valence     = bf.valence.empty() ? nullptr : bf.valence.data();
        bd.shell_hardness = bf.shell_hardness.data();
        bd.selfenergy  = hf.selfenergy.data();
        bd.kcn         = hf.kcn.data();
        bd.shpoly      = hf.shpoly.data();
        bd.ao2at       = bf.ao2at.empty() ? nullptr : bf.ao2at.data();   // Stage 4 (gradient)
        bd.ao2sh       = bf.ao2sh.empty() ? nullptr : bf.ao2sh.data();
        bd.rep_alpha   = bf.rep_alpha.empty() ? nullptr : bf.rep_alpha.data();
        bd.rep_zeff    = bf.rep_zeff.empty() ? nullptr : bf.rep_zeff.data();
        return m_ctx->beginBasis(bd);
    }
    bool beginComputed(const std::vector<double>& xyz_bohr) override
    {
        return m_ctx && m_ctx->beginComputed(xyz_bohr.data());
    }
    bool downloadOverlap(Eigen::MatrixXd& S_out) override
    {
        if (!m_ctx || m_n <= 0) return false;
        S_out.resize(m_n, m_n);
        return m_ctx->downloadOverlap(S_out.data());
    }
    bool downloadH0(Eigen::MatrixXd& H0_out) override
    {
        if (!m_ctx || m_n <= 0) return false;
        H0_out.resize(m_n, m_n);
        return m_ctx->downloadH0(H0_out.data());
    }
    bool downloadGamma(Eigen::MatrixXd& gamma_out) override
    {
        if (!m_ctx || m_nsh <= 0) return false;
        gamma_out.resize(m_nsh, m_nsh);
        return m_ctx->downloadGamma(gamma_out.data());
    }

    // ---- Stage 4: device nuclear gradient ---------------------------------
    // GFN1: v_dp/v_qp empty (isotropic). GFN2: non-empty v_dp/v_qp drive the on-device
    // multipole-integral Pulay term, so the device gradient is complete for GFN2 too (the
    // host then adds the multipole SD/DD/SQ interaction, the CN chain-rule and D4).
    // Density/MO coefficients are resident (pc_resident), so P/C are unused.
    bool supportsGradient() const override { return true; }
    bool gradient(const Matrix& P, const Eigen::MatrixXd& C, const Vector& eps,
                  int nocc_orbs, const Vector& v_ao, const Vector& q_sh,
                  const Eigen::MatrixXd& v_dp, const Eigen::MatrixXd& v_qp,
                  Matrix& grad_out, Vector& dEdcn_out, bool pc_resident) override
    {
        (void)P; (void)C; (void)pc_resident;
        if (!m_ctx || m_nat <= 0) return false;
        if (!deviceGradientEnabled(v_dp.size() > 0)) return false;   // backend escape hatch
        const double* vdp = v_dp.size() > 0 ? v_dp.data() : nullptr;
        const double* vqp = v_qp.size() > 0 ? v_qp.data() : nullptr;
        std::vector<double> grad(3 * static_cast<size_t>(m_nat), 0.0), dEdcn(m_nat, 0.0);
        if (!m_ctx->gradient(eps.data(), nocc_orbs, v_ao.data(), q_sh.data(),
                             vdp, vqp, grad.data(), dEdcn.data()))
            return false;
        grad_out.resize(m_nat, 3);
        for (int i = 0; i < m_nat; ++i) {
            grad_out(i, 0) = grad[3*i+0];
            grad_out(i, 1) = grad[3*i+1];
            grad_out(i, 2) = grad[3*i+2];
        }
        dEdcn_out.resize(m_nat);
        for (int i = 0; i < m_nat; ++i) dEdcn_out(i) = dEdcn[i];
        return true;
    }

    // ---- Stage 3m: GFN2 multipole integrals + resident multipole SCF -------
    // supportsMultipole()→true makes XTB::Calculation enter the resident multipole loop
    // (solveMultipole + density + multipoleMoments) instead of the host SCF, eliminating
    // the per-iteration nao² eigensolver transfer. dp_int/qp_int are built on the device
    // (beginBasis/beginComputed), so beginMultipoleComputed only confirms readiness. The
    // CPU-integral UPLOAD path (beginMultipole) is unsupported on both backends — when the
    // device integral build is unavailable the loop falls back to the host SCF
    // (use_gpu_resident stays false), which is correct.
    bool supportsMultipole() const override { return true; }
    bool beginMultipole(const std::array<Eigen::MatrixXd, 3>& dp_int,
                        const std::array<Eigen::MatrixXd, 6>& qp_int,
                        const std::vector<int>& ao2at) override
    { (void)dp_int; (void)qp_int; (void)ao2at; return false; }
    bool beginMultipoleComputed() override
    {
        return m_ctx && m_ctx->beginMultipoleComputed();
    }
    bool solveMultipole(const Eigen::VectorXd& v_ao, const Eigen::MatrixXd& v_dp,
                        const Eigen::MatrixXd& v_qp, Vector& eps, bool fp32 = false,
                        int n_eig = 0) override
    {
        (void)n_eig;   // full spectrum; fp32 → the backend's FP32 solver (X-AP3)
        if (!m_ctx || static_cast<int>(v_ao.size()) != m_n
            || v_dp.rows() != 3 || v_dp.cols() != m_nat
            || v_qp.rows() != 6 || v_qp.cols() != m_nat) return false;
        eps.resize(m_n);
        return m_ctx->solveMultipole(v_ao.data(), v_dp.data(), v_qp.data(), eps.data(), fp32);
    }
    bool multipoleMoments(Eigen::MatrixXd& dp_at, Eigen::MatrixXd& qp_at) override
    {
        if (!m_ctx || m_nat <= 0) return false;
        dp_at.resize(3, m_nat);
        qp_at.resize(6, m_nat);
        return m_ctx->multipoleMoments(dp_at.data(), qp_at.data());
    }
    bool downloadMultipoleInts(std::array<Eigen::MatrixXd, 3>& dp_int,
                               std::array<Eigen::MatrixXd, 6>& qp_int) override
    {
        if (!m_ctx || m_n <= 0) return false;
        const size_t nn = static_cast<size_t>(m_n) * m_n;
        std::vector<double> dp(3 * nn), qp(6 * nn);   // contiguous 3·nn / 6·nn, col-major
        if (!m_ctx->downloadMultipoleInts(dp.data(), qp.data())) return false;
        for (int k = 0; k < 3; ++k) {                 // Eigen MatrixXd is column-major →
            dp_int[k].resize(m_n, m_n);               // (mu,nu) at mu+nu*nao matches device
            std::memcpy(dp_int[k].data(), dp.data() + static_cast<size_t>(k) * nn, sizeof(double) * nn);
        }
        for (int k = 0; k < 6; ++k) {
            qp_int[k].resize(m_n, m_n);
            std::memcpy(qp_int[k].data(), qp.data() + static_cast<size_t>(k) * nn, sizeof(double) * nn);
        }
        return true;
    }

    // ---- In-SCF GFN2 D4 atom-potential (Stage 5, Part B2) -----------------
    // supportsDeviceDispersion()→true routes XTB::addDispersionPotential's per-iteration
    // dE_D4/dq to the device (xtb_native.cpp use_device_disp path), so the host O(N²) D4
    // contraction drops out of the GFN2 SCF loop.
    bool supportsDeviceDispersion() const override { return true; }
    bool beginDispersion(int nat, const int* Z, const double* sqrtZr4r2,
                         const int* nref, const double* xyz_bohr,
                         const double* c6_flat, int c6_flat_len,
                         double s6, double s8, double a1, double a2, double cutoff) override
    {
        if (!m_ctx) return false;
        return m_ctx->beginDispersion(nat, Z, sqrtZr4r2, nref, xyz_bohr,
                                      c6_flat, c6_flat_len, s6, s8, a1, a2, cutoff);
    }
    bool dispersionDedq(int nat, const double* W, const double* dWq, double* dEdq_out) override
    {
        if (!m_ctx) return false;
        return m_ctx->dispersionDedq(nat, W, dWq, dEdq_out);
    }
    bool dispersionGradient(int nat, const double* W, const double* dWq, const double* dWc,
                            double* e_atom_out, double* grad_out,
                            double* dEdcn_out, double* dEdq_out) override
    {
        if (!m_ctx) return false;
        return m_ctx->dispersionGradient(nat, W, dWq, dWc, e_atom_out, grad_out, dEdcn_out, dEdq_out);
    }
    bool dispersionATM(int nat, const double* c6, const double* dc6dcn,
                       double s9, double a1, double a2, double alp, double cutoff,
                       double* e_atom_out, double* grad_out, double* dEdcn_out) override
    {
        if (!m_ctx) return false;
        return m_ctx->dispersionATM(nat, c6, dc6dcn, s9, a1, a2, alp, cutoff,
                                    e_atom_out, grad_out, dEdcn_out);
    }

    // ---- Single-shot D4 EEQ charge model (Stage 5, Part A) ----------------
    // supportsDeviceEeq()→true routes the EEQ initial guess (scf_guess=eeq) and the D4
    // q-response (∂q/∂x in calcDispersionEnergy) onto the GPU, completing the device D4
    // gradient — the last host-resident piece.
    bool supportsDeviceEeq() const override { return true; }
    bool eeqCharges(int N, const double* xyz_bohr,
                    const double* chi, const double* gam, const double* alpha_sq,
                    const double* cnf, const double* rcov_bohr,
                    double total_charge, double* q_out) override
    {
        if (!m_ctx) return false;
        return m_ctx->eeqCharges(N, xyz_bohr, chi, gam, alpha_sq, cnf, rcov_bohr,
                                 total_charge, q_out);
    }
    bool eeqChargeResponse(int N, const double* dEdq, double* grad_add) override
    {
        if (!m_ctx) return false;
        return m_ctx->eeqChargeResponseGradient(N, dEdq, grad_add);
    }

protected:
    /**
     * @brief Backend veto on the device gradient (default: always run it on the device).
     *
     * The one behavioural difference between the ROCm and Vulkan gradient paths: Vulkan
     * keeps a debug escape hatch that forces the host gradient for GFN2. Overriding this
     * one predicate is cheaper than duplicating the whole 25-line gradient() marshalling.
     */
    virtual bool deviceGradientEnabled(bool is_gfn2) const { (void)is_gfn2; return true; }

    Context* m_ctx = nullptr;
    int      m_n   = 0;   ///< nao
    int      m_nsh = 0;   ///< number of shells
    int      m_nat = 0;   ///< number of atoms
};

} // namespace gpu
} // namespace xtb
} // namespace curcuma
