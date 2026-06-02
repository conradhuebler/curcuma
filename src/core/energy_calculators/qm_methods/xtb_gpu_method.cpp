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

using curcuma::xtb::MethodType;
using curcuma::xtb::gpu::XtbGpuContext;

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
