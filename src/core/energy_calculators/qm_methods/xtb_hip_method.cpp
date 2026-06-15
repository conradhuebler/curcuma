/*
 * <Native xTB ROCm/HIP Method Wrapper — implementation>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0.
 *
 * Claude Generated (2026-06): Stage-0 pass-through. Forwards the whole
 * ComputationalMethod interface to an owned NativeXtbMethod and initializes the HIP
 * context. Host-compiled (no device code) — the HIP kernels live in the .hip TU.
 */

#ifdef USE_ROCM_XTB

#include "xtb_hip_method.h"
#include "rocm/xtb_hip_context.h"

#include "src/core/curcuma_logger.h"

#include <fmt/format.h>
#include <array>
#include <cstring>
#include <vector>

namespace {
using curcuma::xtb::MethodType;
using curcuma::xtb::gpu::XtbHipContext;

#ifdef HAVE_ROCSOLVER
/**
 * @brief Device-resident GFN1 SCF backend (Stage 2) over an XtbHipContext.
 *
 * Adapts the column-major resident kernels to the project-typed GpuScfBackend the core
 * XTB calls. supportsMultipole() stays false, so GFN2 falls back to the Stage-1
 * ExternalEigensolver path. H0/S are symmetric, so the Eigen column-major copy preserves
 * values. Claude Generated (Stage 2).
 */
class HipScfBackend : public curcuma::xtb::GpuScfBackend {
public:
    explicit HipScfBackend(XtbHipContext* ctx) : m_ctx(ctx) {}

    bool begin(const Matrix& H0, const Matrix& S, const Eigen::MatrixXd& L) override
    {
        const int n = static_cast<int>(H0.rows());
        if (!m_ctx || n <= 0 || S.rows() != n) return false;
        (void)L;  // rocSOLVER dsygvd takes S directly (no host Cholesky needed)
        m_n = n;
        Eigen::MatrixXd Hcm = H0, Scm = S;  // column-major copies (H0/S symmetric)
        return m_ctx->residentBegin(Hcm.data(), Scm.data(), n);
    }
    bool solve(const Eigen::VectorXd& v_ao, Vector& eps, bool fp32 = false, int n_eig = 0) override
    {
        (void)n_eig;  // full spectrum; fp32 → rocsolver_ssygvd (X-AP3 mixed precision)
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
        curcuma::xtb::gpu::XtbHipBasisData bd;
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
        bd.ao2at       = bf.ao2at.empty() ? nullptr : bf.ao2at.data();   // Stage 4
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
    bool downloadCholesky(Eigen::MatrixXd& L_out) override
    {
        if (!m_ctx || m_n <= 0) return false;
        L_out.resize(m_n, m_n);
        return m_ctx->downloadCholesky(L_out.data());
    }
    bool downloadGamma(Eigen::MatrixXd& gamma_out) override
    {
        if (!m_ctx || m_nsh <= 0) return false;
        gamma_out.resize(m_nsh, m_nsh);
        return m_ctx->downloadGamma(gamma_out.data());
    }

    // ---- Stage 4: device nuclear gradient (GFN1) --------------------------
    bool supportsGradient() const override { return true; }
    bool gradient(const Matrix& P, const Eigen::MatrixXd& C, const Vector& eps,
                  int nocc_orbs, const Vector& v_ao, const Vector& q_sh,
                  const Eigen::MatrixXd& v_dp, const Eigen::MatrixXd& v_qp,
                  Matrix& grad_out, Vector& dEdcn_out, bool pc_resident) override
    {
        // GFN1 isotropic: density/MO coefficients are resident (pc_resident), no multipole.
        (void)P; (void)C; (void)pc_resident;
        if (!m_ctx || m_nat <= 0) return false;
        // GFN1 only: the device kernels omit the GFN2 multipole-integral Pulay term (R-AP3
        // pending). For GFN2 (non-empty v_dp/v_qp) return false so the caller falls back to
        // the full host gradient (the resident GFN2 path made use_gpu_resident=true).
        if (v_dp.size() > 0 || v_qp.size() > 0) return false;
        std::vector<double> grad(3 * static_cast<size_t>(m_nat), 0.0), dEdcn(m_nat, 0.0);
        if (!m_ctx->gradient(eps.data(), nocc_orbs, v_ao.data(), q_sh.data(),
                             grad.data(), dEdcn.data()))
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

    // ---- Stage 3m (R-AP1): GFN2 multipole integrals on device --------------
    // supportsMultipole() stays false (ROCm GFN2 runs the host SCF with the rocSOLVER
    // per-iteration eigensolver), so these are reached only via the non-resident host-
    // download path in xtb_native.cpp: build dp_int/qp_int on the device, download so the
    // CPU setupMultipole integral loop is skipped.
    // R-AP2: device-resident GFN2 multipole SCF. supportsMultipole()→true makes
    // XTB::Calculation enter the resident multipole loop (solveMultipole + density +
    // multipoleMoments) instead of the Stage-1 host SCF. dp_int/qp_int are built on the
    // device (computeIntegrals); beginMultipoleComputed confirms readiness. The CPU-
    // integral UPLOAD path (beginMultipole) is unsupported — when the device integral
    // build is unavailable the loop falls back to the host SCF (use_gpu_resident=false).
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
        (void)n_eig;   // full spectrum; fp32 → rocsolver_ssygvd (X-AP3 mixed precision)
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
        for (int k = 0; k < 3; ++k) {                 // Eigen MatrixXd column-major →
            dp_int[k].resize(m_n, m_n);               // (mu,nu) at mu+nu*nao matches device
            std::memcpy(dp_int[k].data(), dp.data() + static_cast<size_t>(k) * nn, sizeof(double) * nn);
        }
        for (int k = 0; k < 6; ++k) {
            qp_int[k].resize(m_n, m_n);
            std::memcpy(qp_int[k].data(), qp.data() + static_cast<size_t>(k) * nn, sizeof(double) * nn);
        }
        return true;
    }
private:
    XtbHipContext* m_ctx = nullptr;
    int            m_n   = 0;
    int            m_nsh = 0;
    int            m_nat = 0;
};
#endif // HAVE_ROCSOLVER
}  // namespace

XtbHipComputationalMethod::XtbHipComputationalMethod(MethodType method, const json& config)
    : m_method(method)
{
    // Device handshake (non-throwing; ok() is false when no usable device).
    m_gpu = std::make_unique<XtbHipContext>();

    // Full validated CPU pipeline (config, large-system modes, errors, properties).
    m_cpu = std::make_unique<NativeXtbMethod>(method, config);

    if (m_gpu->ok()) {
        if (CurcumaLogger::get_verbosity() >= 1)
            CurcumaLogger::success(fmt::format(
                "{}: ROCm context ready on device {} ({})",
                getMethodName(), m_gpu->deviceId(), m_gpu->deviceName()));

#ifdef HAVE_ROCSOLVER
        // Stage 1: install the GPU eigensolver. solveEigen() delegates the per-iteration
        // generalized eigenproblem (F, S = L·Lᵀ) → (C, eps) to rocSOLVER (dsygvd) on the
        // device; everything else (integrals, Fock, density, gradient) stays on the CPU.
        // rocSOLVER solves the generalized problem directly, so this works for both GFN1
        // and GFN2 (the host hands it the complete Fock).
        if (curcuma::xtb::XTB* xtb = m_cpu->solver()) {
            curcuma::xtb::gpu::XtbHipContext* ctx = m_gpu.get();
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
            // Stage 2: install the device-resident GFN1 SCF backend. Under the default
            // Broyden mixing, XTB::Calculation keeps H0/S (and the per-iteration density
            // and MO coefficients) resident on the device — the Fock build + populations
            // are HIP kernels, the density a rocBLAS GEMM, the eigensolve rocSOLVER. GFN2
            // (supportsMultipole()==false here) uses the Stage-1 eigensolver hook above.
            m_scf_backend = std::make_unique<HipScfBackend>(ctx);
            xtb->setGpuScfBackend(m_scf_backend.get());
            // FP64 is ~1/16 of FP32 on this iGPU and the per-iteration eigensolve
            // dominates the resident SCF, so default mixed precision ON for the GPU
            // path: far-from-convergence iterations solve in FP32 (rocsolver_ssygvd),
            // reverting to FP64 once max|dq| < scf_fp32_threshold so the converged
            // fixed point and energy stay FP64. Matches the CUDA backend (X-AP3).
            xtb->setMixedPrecision(true);

            if (CurcumaLogger::get_verbosity() >= 2)
                CurcumaLogger::info(fmt::format(
                    "{}: ROCm fully device-resident GFN1 (rocSOLVER + HIP kernels): integral "
                    "build (CN/S/H0/L/gamma) + SCF + nuclear gradient on the GPU; only the "
                    "dispersion gradient + CN chain-rule on CPU (Stage 4)", getMethodName()));
        }
#else
        if (CurcumaLogger::get_verbosity() >= 2)
            CurcumaLogger::info(fmt::format(
                "{}: ROCm Stage 0 (device handshake); SCF/integrals/gradient on CPU "
                "(build without rocSOLVER)", getMethodName()));
#endif
    } else {
        CurcumaLogger::warn(fmt::format(
            "{}: no usable ROCm/HIP device; running CPU path", getMethodName()));
    }
}

XtbHipComputationalMethod::~XtbHipComputationalMethod() = default;

// ---- forwarding (Stage 0) -------------------------------------------------
bool XtbHipComputationalMethod::setMolecule(const Mol& mol) { return m_cpu->setMolecule(mol); }
bool XtbHipComputationalMethod::updateGeometry(const Matrix& g) { return m_cpu->updateGeometry(g); }
double XtbHipComputationalMethod::calculateEnergy(bool gradient) { return m_cpu->calculateEnergy(gradient); }

Matrix XtbHipComputationalMethod::getGradient() const { return m_cpu->getGradient(); }
Vector XtbHipComputationalMethod::getCharges() const { return m_cpu->getCharges(); }
Vector XtbHipComputationalMethod::getBondOrders() const { return m_cpu->getBondOrders(); }
Position XtbHipComputationalMethod::getDipole() const { return m_cpu->getDipole(); }
bool XtbHipComputationalMethod::hasGradient() const { return m_cpu->hasGradient(); }

std::string XtbHipComputationalMethod::getMethodName() const { return m_cpu->getMethodName(); }
bool XtbHipComputationalMethod::isThreadSafe() const { return m_cpu->isThreadSafe(); }
void XtbHipComputationalMethod::setThreadCount(int threads) { m_cpu->setThreadCount(threads); }

void XtbHipComputationalMethod::setParameters(const json& params) { m_cpu->setParameters(params); }
json XtbHipComputationalMethod::getParameters() const { return m_cpu->getParameters(); }

bool XtbHipComputationalMethod::hasError() const { return m_cpu->hasError(); }
void XtbHipComputationalMethod::clearError() { m_cpu->clearError(); }
std::string XtbHipComputationalMethod::getErrorMessage() const { return m_cpu->getErrorMessage(); }

Vector XtbHipComputationalMethod::getOrbitalEnergies() const { return m_cpu->getOrbitalEnergies(); }
int XtbHipComputationalMethod::getNumElectrons() const { return m_cpu->getNumElectrons(); }
json XtbHipComputationalMethod::getEnergyDecomposition() const { return m_cpu->getEnergyDecomposition(); }
bool XtbHipComputationalMethod::saveToFile(const std::string& f) const { return m_cpu->saveToFile(f); }

void XtbHipComputationalMethod::setWarmStart(bool on) { m_cpu->setWarmStart(on); }
void XtbHipComputationalMethod::setIterativeMode(bool on) { m_cpu->setIterativeMode(on); }

bool XtbHipComputationalMethod::gpuActive() const { return m_gpu && m_gpu->ok(); }

#endif // USE_ROCM_XTB
