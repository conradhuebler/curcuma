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

namespace {
using curcuma::xtb::MethodType;
using curcuma::xtb::gpu::XtbHipContext;
}

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
            if (CurcumaLogger::get_verbosity() >= 2)
                CurcumaLogger::info(fmt::format(
                    "{}: ROCm GPU eigensolver active (rocSOLVER dsygvd, FP64); "
                    "SCF/integrals/gradient on CPU (Stage 1)", getMethodName()));
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
