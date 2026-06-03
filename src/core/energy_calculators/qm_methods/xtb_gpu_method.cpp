/*
 * <Native xTB GPU Method Wrapper — implementation>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0.
 *
 * Claude Generated (2026-06): Stage-0 pass-through. Forwards the whole
 * ComputationalMethod interface to an owned NativeXtbMethod and initializes the
 * GPU context. Host-compiled (no device code) — the CUDA kernels live in the
 * .cu translation units.
 */

#ifdef USE_CUDA_XTB

#include "xtb_gpu_method.h"
#include "cuda/xtb_gpu_context.h"

#include "src/core/curcuma_logger.h"

#include <fmt/format.h>

#include <algorithm>
#include <cstring>
#include <vector>

using curcuma::xtb::MethodType;
using curcuma::xtb::gpu::XtbGpuContext;

namespace {

/**
 * @brief Concrete device-resident SCF backend (GPU port Stage 2, Claude Generated).
 *
 * Adapts the column-major XtbGpuContext resident kernels to the project-typed
 * GpuScfBackend interface the core XTB calls. Owns nothing but a context pointer
 * (the wrapper owns both). H0/S are symmetric, so the row-major→column-major copy
 * (Eigen::MatrixXd = H0) preserves values; L (m_X) is already column-major.
 */
class XtbGpuScfBackend : public curcuma::xtb::GpuScfBackend {
public:
    explicit XtbGpuScfBackend(XtbGpuContext* ctx) : m_ctx(ctx) {}

    bool begin(const Matrix& H0, const Matrix& S, const Eigen::MatrixXd& L) override
    {
        const int n = static_cast<int>(H0.rows());
        if (!m_ctx || n <= 0 || S.rows() != n || L.rows() != n || L.cols() != n)
            return false;
        m_n = n;
        // Column-major copies of the symmetric H0/S; Hcm.data()[i+j*n] == H0(i,j).
        Eigen::MatrixXd Hcm = H0;
        Eigen::MatrixXd Scm = S;
        return m_ctx->residentBegin(Hcm.data(), Scm.data(), L.data(), n);
    }

    bool solve(const Eigen::VectorXd& v_ao, Vector& eps) override
    {
        if (!m_ctx || static_cast<int>(v_ao.size()) != m_n) return false;
        eps.resize(m_n);
        return m_ctx->residentSolve(v_ao.data(), m_n, eps.data());
    }

    bool density(const Eigen::VectorXd& occ, int ncol,
                 Eigen::VectorXd& pop_ao, double& band) override
    {
        if (!m_ctx || ncol < 0 || ncol > m_n) return false;
        pop_ao.resize(m_n);
        band = 0.0;
        return m_ctx->residentDensity(occ.data(), ncol, m_n, pop_ao.data(), &band);
    }

    bool finalize(Matrix& P, Matrix& C) override
    {
        if (!m_ctx) return false;
        Eigen::MatrixXd Pcm(m_n, m_n), Ccm(m_n, m_n);  // column-major device layout
        if (!m_ctx->residentFinalize(Pcm.data(), Ccm.data(), m_n)) return false;
        P = Pcm;   // → row-major project Matrix, values preserved (P symmetric)
        C = Ccm;   // col-major eigenvectors → row-major, values preserved
        return true;
    }

    // ---- GFN2 multipole (Stage 2b) ----------------------------------------
    bool supportsMultipole() const override { return true; }

    bool beginMultipole(const std::array<Eigen::MatrixXd, 3>& dp_int,
                        const std::array<Eigen::MatrixXd, 6>& qp_int,
                        const std::vector<int>& ao2at) override
    {
        if (!m_ctx || m_n <= 0 || static_cast<int>(ao2at.size()) != m_n) return false;
        const size_t nn = static_cast<size_t>(m_n) * static_cast<size_t>(m_n);
        // Pack the 3 dipole / 6 quadrupole integral matrices contiguously. Each
        // Eigen::MatrixXd is column-major n×n, so block k is its raw buffer.
        std::vector<double> dp(3 * nn), qp(6 * nn);
        for (int k = 0; k < 3; ++k) {
            if (dp_int[k].rows() != m_n || dp_int[k].cols() != m_n) return false;
            std::memcpy(dp.data() + static_cast<size_t>(k) * nn, dp_int[k].data(), nn * sizeof(double));
        }
        for (int k = 0; k < 6; ++k) {
            if (qp_int[k].rows() != m_n || qp_int[k].cols() != m_n) return false;
            std::memcpy(qp.data() + static_cast<size_t>(k) * nn, qp_int[k].data(), nn * sizeof(double));
        }
        m_nat = *std::max_element(ao2at.begin(), ao2at.end()) + 1;
        if (m_nat <= 0) return false;
        return m_ctx->residentBeginMultipole(dp.data(), qp.data(), ao2at.data(), m_n, m_nat);
    }

    bool solveMultipole(const Eigen::VectorXd& v_ao, const Eigen::MatrixXd& v_dp,
                        const Eigen::MatrixXd& v_qp, Vector& eps) override
    {
        if (!m_ctx || static_cast<int>(v_ao.size()) != m_n
            || v_dp.rows() != 3 || v_dp.cols() != m_nat
            || v_qp.rows() != 6 || v_qp.cols() != m_nat) return false;
        eps.resize(m_n);
        return m_ctx->residentSolveMultipole(v_ao.data(), v_dp.data(), v_qp.data(),
                                             m_n, eps.data());
    }

    bool multipoleMoments(Eigen::MatrixXd& dp_at, Eigen::MatrixXd& qp_at) override
    {
        if (!m_ctx || m_nat <= 0) return false;
        dp_at.resize(3, m_nat);
        qp_at.resize(6, m_nat);
        return m_ctx->residentMultipoleMoments(dp_at.data(), qp_at.data(), m_n, m_nat);
    }

    // ---- Device-side integral build (Stage 3) -----------------------------
    bool beginBasis(const curcuma::xtb::GpuBasisFlat& bf,
                    const curcuma::xtb::GpuH0Flat& hf) override
    {
        if (!m_ctx || bf.nao <= 0) return false;
        // Own copies so the device-upload pointers stay valid past this call.
        m_bf = bf;
        m_hf = hf;
        XtbGpuContext::XtbGpuBasisData bd;
        bd.nat         = m_bf.nat;
        bd.nsh         = m_bf.nsh;
        bd.nao         = m_bf.nao;
        bd.is_gfn2     = m_bf.is_gfn2;
        bd.z           = m_bf.z.data();
        bd.sh2at       = m_bf.sh2at.data();
        bd.ang_sh      = m_bf.ang_sh.data();
        bd.iao_sh      = m_bf.iao_sh.data();
        bd.nao_sh      = m_bf.nao_sh.data();
        bd.sh_nprim    = m_bf.sh_nprim.data();
        bd.sh_prim_off = m_bf.sh_prim_off.data();
        bd.nprim_total = static_cast<int>(m_bf.prim_alpha.size());
        bd.prim_alpha  = m_bf.prim_alpha.data();
        bd.prim_coeff  = m_bf.prim_coeff.data();
        bd.sh_zeta     = m_bf.sh_zeta.data();
        bd.selfenergy  = m_hf.selfenergy.data();
        bd.kcn         = m_hf.kcn.data();
        bd.shpoly      = m_hf.shpoly.data();
        bd.valence     = m_bf.valence.empty() ? nullptr : m_bf.valence.data();
        bd.shell_hardness = m_bf.shell_hardness.empty() ? nullptr : m_bf.shell_hardness.data();
        bd.ao2at       = m_bf.ao2at.empty() ? nullptr : m_bf.ao2at.data();
        bd.ao2sh       = m_bf.ao2sh.empty() ? nullptr : m_bf.ao2sh.data();
        bd.rep_alpha   = m_bf.rep_alpha.empty() ? nullptr : m_bf.rep_alpha.data();
        bd.rep_zeff    = m_bf.rep_zeff.empty() ? nullptr : m_bf.rep_zeff.data();
        m_n = m_bf.nao;
        m_nat = m_bf.nat;
        return m_ctx->beginBasis(bd);
    }

    bool beginComputed(const std::vector<double>& xyz_bohr) override
    {
        if (!m_ctx || m_n <= 0) return false;
        return m_ctx->computeIntegrals(xyz_bohr.data())
            && m_ctx->residentBeginComputed();
    }

    bool downloadGamma(Eigen::MatrixXd& gamma_out) override
    {
        if (!m_ctx) return false;
        const int nsh = static_cast<int>(m_bf.nsh);
        if (nsh <= 0) return false;
        gamma_out.resize(nsh, nsh);  // Eigen column-major matches the device layout
        return m_ctx->downloadGamma(gamma_out.data());
    }

    bool beginMultipoleComputed() override
    {
        if (!m_ctx || m_nat <= 0) return false;
        m_nat = m_bf.nat;  // matches the resident_nat set by residentBeginMultipoleComputed
        return m_ctx->residentBeginMultipoleComputed();
    }

    // ---- Device nuclear gradient (Stage 4) --------------------------------
    bool supportsGradient() const override { return true; }

    bool gradient(const Matrix& P, const Eigen::MatrixXd& C, const Vector& eps,
                  int nocc_orbs, const Vector& v_ao, const Vector& q_sh,
                  const Eigen::MatrixXd& v_dp, const Eigen::MatrixXd& v_qp,
                  Matrix& grad_out, Vector& dEdcn_out) override
    {
        if (!m_ctx || m_n <= 0 || m_nat <= 0) return false;
        if (P.rows() != m_n || C.rows() != m_n || static_cast<int>(eps.size()) < nocc_orbs
            || static_cast<int>(v_ao.size()) != m_n) return false;
        Eigen::MatrixXd Pcm = P;  // column-major copy (P symmetric → values preserved)
        // v_dp/v_qp are Eigen::MatrixXd (column-major) → contiguous [k + iat*rows].
        const bool with_mp = (v_dp.rows() == 3 && v_dp.cols() == m_nat
                           && v_qp.rows() == 6 && v_qp.cols() == m_nat);
        std::vector<double> grad(3 * m_nat, 0.0), dEdcn(m_nat, 0.0);
        if (!m_ctx->computeGradient(Pcm.data(), C.data(), eps.data(), nocc_orbs,
                                    v_ao.data(), q_sh.data(),
                                    with_mp ? v_dp.data() : nullptr,
                                    with_mp ? v_qp.data() : nullptr,
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

private:
    XtbGpuContext* m_ctx = nullptr;
    int            m_n   = 0;
    int            m_nat = 0;
    curcuma::xtb::GpuBasisFlat m_bf;
    curcuma::xtb::GpuH0Flat    m_hf;
};

} // namespace

// Public factory (Stage 4 validation, Claude Generated): build a device-resident
// SCF+gradient backend over an existing context, so a bare XTB can be driven on
// the GPU path from a test without the full XtbGpuComputationalMethod wrapper.
std::unique_ptr<curcuma::xtb::GpuScfBackend>
createXtbGpuScfBackend(curcuma::xtb::gpu::XtbGpuContext* ctx)
{
    return std::make_unique<XtbGpuScfBackend>(ctx);
}

XtbGpuComputationalMethod::XtbGpuComputationalMethod(MethodType method, const json& config)
    : m_method(method)
{
    // Device handshake (non-throwing; ok() is false when no usable device).
    m_gpu = std::make_unique<XtbGpuContext>();

    // Full validated CPU pipeline (config, large-system modes, errors, properties).
    m_cpu = std::make_unique<NativeXtbMethod>(method, config);

    if (m_gpu->ok()) {
        if (CurcumaLogger::get_verbosity() >= 1)
            CurcumaLogger::success(fmt::format(
                "{}: GPU context ready on device {} ({})",
                getMethodName(), m_gpu->deviceId(), m_gpu->deviceName()));

        // Stage 1: install the GPU eigensolver. solveEigen() delegates the
        // per-iteration generalized eigenproblem (F, S=L·Lᵀ)→(C, eps) to the
        // device; everything else (integrals, Fock build, density, gradient)
        // still runs on the CPU pipeline. The dense XTB exists at construction;
        // when a large_system_mode driver is active solver() is the (unused)
        // dense instance, so the hook is harmless there (GPU + large_system_mode
        // is not yet wired and runs on the CPU fragment driver).
        if (curcuma::xtb::XTB* xtb = m_cpu->solver()) {
            XtbGpuContext* ctx = m_gpu.get();
            xtb->setExternalEigensolver(
                [ctx](const Matrix& F, const Eigen::MatrixXd& L,
                      Matrix& C, Vector& eps) -> bool {
                    const int n = static_cast<int>(F.rows());
                    if (n <= 0 || L.rows() != n || L.cols() != n)
                        return false;
                    // cuSOLVER/cuBLAS are column-major. F (row-major Matrix) is
                    // symmetric, so its buffer already reads as column-major F.
                    // The eigenvectors come back column-major in Ccm; assigning to
                    // the row-major m_wfn.C (C) lets Eigen transpose-convert, values
                    // preserved. eps is a plain vector (no storage-order ambiguity).
                    Eigen::MatrixXd Fcm = F;     // guaranteed contiguous column-major
                    Eigen::MatrixXd Ccm(n, n);   // column-major eigenvector output
                    eps.resize(n);
                    if (!ctx->solveGeneralizedEigenF64(Fcm.data(), L.data(), n,
                                                       Ccm.data(), eps.data()))
                        return false;
                    C = Ccm;                     // col-major → row-major, values kept
                    return true;
                });
            if (CurcumaLogger::get_verbosity() >= 2)
                CurcumaLogger::info(fmt::format(
                    "{}: GPU eigensolver active (cuSOLVER Dsyevd, FP64); "
                    "SCF/integrals/gradient on CPU (Stage 1)", getMethodName()));

            // Stage 2: install the device-resident SCF backend. With the default
            // Broyden charge mixing, XTB::Calculation keeps H0/S/L and the
            // density/MO matrices on the device for the whole SCF (begin/solve/
            // density/finalize), so only length-nao vectors cross the bus each
            // iteration instead of the per-iteration Fock upload of Stage 1. GFN1
            // runs the isotropic loop (Stage 2a); GFN2 additionally keeps the
            // multipole integrals resident and adds the anisotropic Fock + atomic
            // moments on the device (Stage 2b). Non-Broyden modes keep the
            // Stage-1 eigensolver hook above.
            m_scf_backend = std::make_unique<XtbGpuScfBackend>(ctx);
            xtb->setGpuScfBackend(m_scf_backend.get());
            if (CurcumaLogger::get_verbosity() >= 2)
                CurcumaLogger::info(fmt::format(
                    "{}: GPU device-resident SCF backend active (Broyden; "
                    "GFN1 Stage 2a, GFN2 Stage 2b)", getMethodName()));
        }
    } else {
        CurcumaLogger::warn(fmt::format(
            "{}: no usable CUDA device; running CPU path", getMethodName()));
    }
}

XtbGpuComputationalMethod::~XtbGpuComputationalMethod() = default;

// ---- forwarding (Stage 0) -------------------------------------------------
bool XtbGpuComputationalMethod::setMolecule(const Mol& mol) { return m_cpu->setMolecule(mol); }
bool XtbGpuComputationalMethod::updateGeometry(const Matrix& g) { return m_cpu->updateGeometry(g); }
double XtbGpuComputationalMethod::calculateEnergy(bool gradient) { return m_cpu->calculateEnergy(gradient); }

Matrix XtbGpuComputationalMethod::getGradient() const { return m_cpu->getGradient(); }
Vector XtbGpuComputationalMethod::getCharges() const { return m_cpu->getCharges(); }
Vector XtbGpuComputationalMethod::getBondOrders() const { return m_cpu->getBondOrders(); }
Position XtbGpuComputationalMethod::getDipole() const { return m_cpu->getDipole(); }
bool XtbGpuComputationalMethod::hasGradient() const { return m_cpu->hasGradient(); }

std::string XtbGpuComputationalMethod::getMethodName() const { return m_cpu->getMethodName(); }
bool XtbGpuComputationalMethod::isThreadSafe() const { return m_cpu->isThreadSafe(); }
void XtbGpuComputationalMethod::setThreadCount(int threads) { m_cpu->setThreadCount(threads); }

void XtbGpuComputationalMethod::setParameters(const json& params) { m_cpu->setParameters(params); }
json XtbGpuComputationalMethod::getParameters() const { return m_cpu->getParameters(); }

bool XtbGpuComputationalMethod::hasError() const { return m_cpu->hasError(); }
void XtbGpuComputationalMethod::clearError() { m_cpu->clearError(); }
std::string XtbGpuComputationalMethod::getErrorMessage() const { return m_cpu->getErrorMessage(); }

Vector XtbGpuComputationalMethod::getOrbitalEnergies() const { return m_cpu->getOrbitalEnergies(); }
int XtbGpuComputationalMethod::getNumElectrons() const { return m_cpu->getNumElectrons(); }
json XtbGpuComputationalMethod::getEnergyDecomposition() const { return m_cpu->getEnergyDecomposition(); }
bool XtbGpuComputationalMethod::saveToFile(const std::string& f) const { return m_cpu->saveToFile(f); }

void XtbGpuComputationalMethod::setWarmStart(bool on) { m_cpu->setWarmStart(on); }
void XtbGpuComputationalMethod::setIterativeMode(bool on) { m_cpu->setIterativeMode(on); }

bool XtbGpuComputationalMethod::gpuActive() const { return m_gpu && m_gpu->ok(); }

#endif // USE_CUDA_XTB
