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
private:
    XtbVulkanContext* m_ctx = nullptr;
    int               m_n   = 0;
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

            if (CurcumaLogger::get_verbosity() >= 2)
                CurcumaLogger::info(fmt::format(
                    "{}: Vulkan GPU eigensolver + device-resident GFN1 SCF active "
                    "(FP64 Jacobi / Lowdin); integrals/gradient on CPU (Stage 2)",
                    getMethodName()));
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
