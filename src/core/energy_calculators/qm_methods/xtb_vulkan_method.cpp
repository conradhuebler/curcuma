/*
 * <Native xTB Vulkan Method Wrapper — implementation>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0.
 *
 * Claude Generated (2026-06): Stage-0 pass-through. Forwards the whole
 * ComputationalMethod interface to an owned NativeXtbMethod and initializes the Vulkan
 * context. Host-compiled (no Vulkan headers — the context is pimpl).
 */

#ifdef USE_VULKAN_XTB

#include "xtb_vulkan_method.h"
#include "vulkan/xtb_vulkan_context.h"

#include "src/core/curcuma_logger.h"

#include <fmt/format.h>
#include <array>
#include <cstring>
#include <vector>

namespace {
using curcuma::xtb::MethodType;
using curcuma::xtb::gpu::XtbVulkanContext;

/**
 * @brief Device-resident GFN1 SCF backend (Stage 2) over an XtbVulkanContext.
 *
 * Adapts the column-major resident kernels to the project-typed GpuScfBackend the core
 * XTB calls. supportsMultipole() stays false, so GFN2 falls back to the Stage-1
 * ExternalEigensolver path instead of this isotropic loop. H0/S are symmetric, so the
 * Eigen column-major copy preserves values. Claude Generated (Stage 2a).
 */
class VulkanScfBackend : public curcuma::xtb::GpuScfBackend {
public:
    explicit VulkanScfBackend(XtbVulkanContext* ctx) : m_ctx(ctx) {}

    bool begin(const Matrix& H0, const Matrix& S, const Eigen::MatrixXd& L) override
    {
        const int n = static_cast<int>(H0.rows());
        if (!m_ctx || n <= 0 || S.rows() != n) return false;
        (void)L;  // Löwdin S^-1/2 built on-device from S instead of the host Cholesky
        m_n = n;
        Eigen::MatrixXd Hcm = H0, Scm = S;  // column-major copies (H0/S symmetric)
        return m_ctx->residentBegin(Hcm.data(), Scm.data(), n);
    }
    bool solve(const Eigen::VectorXd& v_ao, Vector& eps, bool fp32 = false, int n_eig = 0) override
    {
        (void)fp32; (void)n_eig;  // always full FP64 spectrum
        if (!m_ctx || static_cast<int>(v_ao.size()) != m_n) return false;
        eps.resize(m_n);
        return m_ctx->residentSolve(v_ao.data(), eps.data());
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
        curcuma::xtb::gpu::XtbVulkanBasisData bd;
        bd.nat = bf.nat; bd.nsh = bf.nsh; bd.nao = bf.nao; bd.is_gfn2 = bf.is_gfn2;
        bd.nprim_total = static_cast<int>(bf.prim_alpha.size());
        bd.z = bf.z.data(); bd.sh2at = bf.sh2at.data(); bd.ang_sh = bf.ang_sh.data();
        bd.iao_sh = bf.iao_sh.data(); bd.nao_sh = bf.nao_sh.data();
        bd.sh_nprim = bf.sh_nprim.data(); bd.sh_prim_off = bf.sh_prim_off.data();
        bd.prim_alpha = bf.prim_alpha.data(); bd.prim_coeff = bf.prim_coeff.data();
        bd.sh_zeta = bf.sh_zeta.data();
        bd.valence = bf.valence.empty() ? nullptr : bf.valence.data();
        bd.shell_hardness = bf.shell_hardness.data();
        bd.selfenergy = hf.selfenergy.data(); bd.kcn = hf.kcn.data(); bd.shpoly = hf.shpoly.data();
        bd.ao2at     = bf.ao2at.empty() ? nullptr : bf.ao2at.data();   // Stage 4 (gradient)
        bd.ao2sh     = bf.ao2sh.empty() ? nullptr : bf.ao2sh.data();
        bd.rep_alpha = bf.rep_alpha.empty() ? nullptr : bf.rep_alpha.data();
        bd.rep_zeff  = bf.rep_zeff.empty() ? nullptr : bf.rep_zeff.data();
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
        // The device builds S; the lower Cholesky L = chol(S) is cheap host-side (Eigen LLT)
        // — no device Cholesky shader needed. Matches the host m_X (orthonormalizer).
        if (!m_ctx || m_n <= 0) return false;
        Eigen::MatrixXd S(m_n, m_n);
        if (!m_ctx->downloadOverlap(S.data())) return false;
        Eigen::LLT<Eigen::MatrixXd> llt(S);
        if (llt.info() != Eigen::Success) return false;
        L_out = llt.matrixL();
        return true;
    }
    bool downloadGamma(Eigen::MatrixXd& gamma_out) override
    {
        if (!m_ctx || m_nsh <= 0) return false;
        gamma_out.resize(m_nsh, m_nsh);
        return m_ctx->downloadGamma(gamma_out.data());
    }

    // ---- Stage 4: device nuclear gradient (GFN1) --------------------------
    // Only reached on the device-resident GFN1 path (XTB gates calculateGradientGpu on
    // use_gpu_resident); GFN2 on Vulkan runs the host SCF (Stage-1 eigensolver) so the
    // CPU gradient is used. Density/MO coefficients are resident (pc_resident), so P/C
    // are empty and v_dp/v_qp unused (GFN1 isotropic).
    bool supportsGradient() const override { return true; }
    bool gradient(const Matrix& P, const Eigen::MatrixXd& C, const Vector& eps,
                  int nocc_orbs, const Vector& v_ao, const Vector& q_sh,
                  const Eigen::MatrixXd& v_dp, const Eigen::MatrixXd& v_qp,
                  Matrix& grad_out, Vector& dEdcn_out, bool pc_resident) override
    {
        (void)P; (void)C; (void)pc_resident;
        if (!m_ctx || m_nat <= 0) return false;
        // GFN1 isotropic only: the device kernels do NOT compute the GFN2 multipole-integral
        // Pulay term (V-AP4 pending). For GFN2 (non-empty v_dp/v_qp) return false so the
        // caller falls back to the full host gradient. Without this, the resident GFN2 path
        // (V-AP3, use_gpu_resident=true) would silently use the isotropic-only gradient.
        if (v_dp.size() > 0 || v_qp.size() > 0) return false;
        std::vector<double> grad(3 * static_cast<size_t>(m_nat), 0.0), dEdcn(m_nat, 0.0);
        if (!m_ctx->gradient(eps.data(), nocc_orbs, v_ao.data(), q_sh.data(),
                             grad.data(), dEdcn.data()))
            return false;
        grad_out.resize(m_nat, 3);
        for (int i = 0; i < m_nat; ++i) {
            grad_out(i, 0) = grad[3*i+0]; grad_out(i, 1) = grad[3*i+1]; grad_out(i, 2) = grad[3*i+2];
        }
        dEdcn_out.resize(m_nat);
        for (int i = 0; i < m_nat; ++i) dEdcn_out(i) = dEdcn[i];
        return true;
    }

    // ---- Stage 3m (V-AP2): GFN2 multipole integrals on device --------------
    // supportsMultipole() stays false (GFN2 runs the host SCF via the Stage-1
    // eigensolver), so these are reached only via the non-resident host-download path
    // in xtb_native.cpp: build dp_int/qp_int on the device, download to the host so the
    // CPU setupMultipole integral loop is skipped.
    // V-AP3: device-resident GFN2 multipole SCF. supportsMultipole()→true makes
    // XTB::Calculation enter the resident multipole loop (solveMultipole + density +
    // multipoleMoments) instead of the Stage-1 host SCF, eliminating the per-iteration
    // nao² eigensolver transfer. dp_int/qp_int are built on the device (computeIntegrals);
    // the resident buffers/sets are set up in beginBasis, so beginMultipoleComputed only
    // confirms readiness. beginMultipole (the CPU-integral UPLOAD path) is unsupported on
    // Vulkan — when the device integral build is unavailable the resident loop falls back
    // to the host SCF (use_gpu_resident stays false), which is correct.
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
        (void)fp32; (void)n_eig;   // always full FP64 spectrum (like the GFN1 resident solve)
        if (!m_ctx || static_cast<int>(v_ao.size()) != m_n
            || v_dp.rows() != 3 || v_dp.cols() != m_nat
            || v_qp.rows() != 6 || v_qp.cols() != m_nat) return false;
        eps.resize(m_n);
        return m_ctx->solveMultipole(v_ao.data(), v_dp.data(), v_qp.data(), eps.data());
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
        if (!m_ctx->downloadMultipole(dp.data(), qp.data())) return false;
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
private:
    XtbVulkanContext* m_ctx = nullptr;
    int               m_n   = 0;
    int               m_nsh = 0;
    int               m_nat = 0;
};
}  // namespace

XtbVulkanComputationalMethod::XtbVulkanComputationalMethod(MethodType method, const json& config)
    : m_method(method)
{
    // Device handshake (non-throwing; ok() is false when no usable device).
    m_gpu = std::make_unique<XtbVulkanContext>();

    // Full validated CPU pipeline (config, large-system modes, errors, properties).
    m_cpu = std::make_unique<NativeXtbMethod>(method, config);

    if (m_gpu->ok()) {
        if (CurcumaLogger::get_verbosity() >= 1)
            CurcumaLogger::success(fmt::format(
                "{}: Vulkan context ready on device {} ({})",
                getMethodName(), m_gpu->deviceId(), m_gpu->deviceName()));

        // Stage 1: install the GPU eigensolver. solveEigen() delegates the per-iteration
        // generalized eigenproblem (F, S=L·Lᵀ)→(C, eps) to the device; everything else
        // (integrals, Fock, density, gradient) stays on the CPU pipeline. The generalized
        // problem reduces to a STANDARD symmetric one via the host Cholesky factor L
        // (two triangular solves + back-transform); only the dense symmetric eigensolve
        // — the SCF hot path — runs on the GPU (FP64 two-sided Jacobi, validated vs Eigen).
        if (curcuma::xtb::XTB* xtb = m_cpu->solver()) {
            curcuma::xtb::gpu::XtbVulkanContext* ctx = m_gpu.get();
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
            // Stage 2: install the device-resident GFN1 SCF backend. Under the default
            // Broyden mixing, XTB::Calculation keeps H0/S (and the per-iteration density
            // and MO coefficients) resident on the device for the whole SCF — only
            // length-nao vectors cross the bus per iteration. GFN1 runs the isotropic
            // loop; GFN2 (supportsMultipole()==false here) falls back to the Stage-1
            // eigensolver hook above.
            m_scf_backend = std::make_unique<VulkanScfBackend>(ctx);
            xtb->setGpuScfBackend(m_scf_backend.get());

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
    } else {
        CurcumaLogger::warn(fmt::format(
            "{}: no usable Vulkan compute device (FP64 required); running CPU path",
            getMethodName()));
    }
}

XtbVulkanComputationalMethod::~XtbVulkanComputationalMethod() = default;

// ---- forwarding (Stage 0) -------------------------------------------------
bool XtbVulkanComputationalMethod::setMolecule(const Mol& mol) { return m_cpu->setMolecule(mol); }
bool XtbVulkanComputationalMethod::updateGeometry(const Matrix& g) { return m_cpu->updateGeometry(g); }
double XtbVulkanComputationalMethod::calculateEnergy(bool gradient) { return m_cpu->calculateEnergy(gradient); }

Matrix XtbVulkanComputationalMethod::getGradient() const { return m_cpu->getGradient(); }
Vector XtbVulkanComputationalMethod::getCharges() const { return m_cpu->getCharges(); }
Vector XtbVulkanComputationalMethod::getBondOrders() const { return m_cpu->getBondOrders(); }
Position XtbVulkanComputationalMethod::getDipole() const { return m_cpu->getDipole(); }
bool XtbVulkanComputationalMethod::hasGradient() const { return m_cpu->hasGradient(); }

std::string XtbVulkanComputationalMethod::getMethodName() const { return m_cpu->getMethodName(); }
bool XtbVulkanComputationalMethod::isThreadSafe() const { return m_cpu->isThreadSafe(); }
void XtbVulkanComputationalMethod::setThreadCount(int threads) { m_cpu->setThreadCount(threads); }

void XtbVulkanComputationalMethod::setParameters(const json& params) { m_cpu->setParameters(params); }
json XtbVulkanComputationalMethod::getParameters() const { return m_cpu->getParameters(); }

bool XtbVulkanComputationalMethod::hasError() const { return m_cpu->hasError(); }
void XtbVulkanComputationalMethod::clearError() { m_cpu->clearError(); }
std::string XtbVulkanComputationalMethod::getErrorMessage() const { return m_cpu->getErrorMessage(); }

Vector XtbVulkanComputationalMethod::getOrbitalEnergies() const { return m_cpu->getOrbitalEnergies(); }
int XtbVulkanComputationalMethod::getNumElectrons() const { return m_cpu->getNumElectrons(); }
json XtbVulkanComputationalMethod::getEnergyDecomposition() const { return m_cpu->getEnergyDecomposition(); }
bool XtbVulkanComputationalMethod::saveToFile(const std::string& f) const { return m_cpu->saveToFile(f); }

void XtbVulkanComputationalMethod::setWarmStart(bool on) { m_cpu->setWarmStart(on); }
void XtbVulkanComputationalMethod::setIterativeMode(bool on) { m_cpu->setIterativeMode(on); }

bool XtbVulkanComputationalMethod::gpuActive() const { return m_gpu && m_gpu->ok(); }

#endif // USE_VULKAN_XTB
