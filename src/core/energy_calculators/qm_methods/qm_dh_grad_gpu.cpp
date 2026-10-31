/*
 * <Double-hybrid gradient: the parts that run on the GPU (WP12 stage 3)>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Device dispatch of the double-hybrid gradient of qm_dh_grad.cpp (equations: docs/DFT_ROADMAP/WP12_doppelhybride.md,
 * stage 2; device notes: stage 3 there and docs/QM_GPU_ROADMAP.md section 4.3):
 *   - F_xc[D] of the Z-vector equations              QMGpuXC::fxc            (cuda/qm_gpu_xc.cu)
 *   - P2-directional XC gradient (basis, points,      QMGpuXC::responseGradient
 *     Becke weights)
 *   - PT2 1e derivative integrals                     Q-G8 (cuda/qm_gpu_grad1e.cu)
 *   - PT2 2e separable part, by polarisation          Q-G7 (cuda/qm_gpu_grad.cu), 1/2 [G(P + P2) - G(P - P2)]
 *   - PT2 2e non-separable part (canonical Gamma)     MP2 2e kernel (cuda/qm_gpu_mp2_grad.cu) with P = P2 = 0
 *   - RI-PT2 3-index part                             cuda/qm_gpu_ri_grad.cu (RIJK::threeIndexGradientDevice)
 * Every part falls back to its CPU function on its own and counts that once (qm_dh_xc, qm_grad1e, qm_grad,
 * qm_mp2_grad, qm_ri_grad). The PT2 intermediates stay on the host (counted qm_dh_grad_pt2 in qm_dh_grad.cpp).
 *
 * Why a file of its own: qm_dh_grad.cpp.o emits out-of-line (weak) copies of the f_xc jets (qm_xc_fxc_kernels.h,
 * through the inline functionalCovered) whose machine code depends on the size of that translation unit (GCC
 * inlining budgets); the device code here would otherwise grow it. This TU does not odr-use the jets. Even so the
 * CPU B95 double-hybrid gradient differs from d0f21497 in the last bits (<= 8.4e-16 Eh/A on H2O/NH3, 2.0e-14 on
 * benzene; energies bitwise); the cause was not found (WP12 doc, stage 3).
 *
 * Claude Generated (Oct 31 2026, WP12 stage 3). AI-generated, machine-tested (RTX 5090 only); human production test
 * pending.
 *
 * This program is free software under GPL-3.0
 */

#include "qm_engine.h"
#include "qm_mp2.h"
#include "qm_ri.h"
#include "qm_xc_fxc.h"

#include "src/core/curcuma_logger.h"
#include "src/core/energy_calculators/gpu_plugin.h"
#include "src/core/gpu_fallback.h"

#include <fmt/format.h>

#include <chrono>
#include <cstdlib>
#include <memory>
#include <string>

namespace {
double msSinceGpu(std::chrono::steady_clock::time_point t0)
{
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}
bool envOn(const char* name)
{
    const char* v = std::getenv(name);
    return v && v[0] == '1';
}
}  // namespace

// Claude Generated (Oct 31 2026, WP12 stage 3): a device F_xc / response gradient failure, counted once.
void QMEngine::gpuDHXCFailed(const std::string& why) const
{
    m_gpu_dh_xc_failed = true;
    curcuma::reportGpuFallback("qm_dh_xc", "device double-hybrid XC response: " + why + "; CPU F_xc and XC response gradient");
}

// Claude Generated (Oct 31 2026, WP12 stage 3): F_xc[D] on the device XC object of the KS SCF (the quantity of
// qmxc::buildFxc). CURCUMA_QM_GPU_XC_CHECK=1 also forms the CPU F_xc of the same D and prints one "QM GPU FXC check"
// line (test_cases/qm_gpu_dft/check_qm_gpu_dft.py --dh). False: the caller runs the CPU builder.
bool QMEngine::dhFxcDevice(const Matrix& D, Matrix& V) const
{
    if (m_gpu_dh_xc_failed || m_xc_frozen_grid || !ensureGpuXC()) return false;
    const auto t0 = std::chrono::steady_clock::now();
    const Matrix P = m_density, Dc = D;  // contiguous
    Matrix Vd(m_nbf, m_nbf);
    if (!m_gpu_xc_obj->fxc(P.data(), Dc.data(), Vd.data())) {
        gpuDHXCFailed("F_xc kernel failed");
        return false;
    }
    V = Vd;
    if (envOn("CURCUMA_QM_GPU_XC_CHECK")) {
        const double dev_ms = msSinceGpu(t0);
        const auto c0 = std::chrono::steady_clock::now();
        Matrix Vc;
        qmxc::buildFxc(m_xc_grid, m_xc_ao, m_density, D, m_xc, m_xc_rho_threshold, m_threads, Vc, &m_xc_screen);
        const double vmax = Vc.cwiseAbs().maxCoeff(), dv = (Vd - Vc).cwiseAbs().maxCoeff();
        CurcumaLogger::result(fmt::format("QM GPU FXC check: n={} max|dV|={:.3e} max|V|={:.3e} rel_dV={:.3e} dev_ms={:.1f} "
                                          "cpu_ms={:.1f}", m_nbf, dv, vmax, vmax > 0.0 ? dv / vmax : dv, dev_ms, msSinceGpu(c0)));
    }
    return true;
}

// Claude Generated (Oct 31 2026, WP12 stage 3): basis, point and Becke weight parts of the P2-directional XC gradient
// on the device (QMGpuXC::responseGradient). CURCUMA_QM_GPU_XC_CHECK=1 also runs the CPU parts and prints one
// "QM GPU DH XC grad check" line. False: the caller runs the CPU path.
bool QMEngine::dhXCResponseGradientDevice(const Matrix& P2, Matrix& g_basis, Matrix& g_response, std::string& why) const
{
    if (m_gpu_dh_xc_failed || m_xc_frozen_grid || !ensureGpuXC()) return false;
    const auto t0 = std::chrono::steady_clock::now();
    Matrix gb = Matrix::Zero(m_atomcount, 3), gp = Matrix::Zero(m_atomcount, 3), gw = Matrix::Zero(m_atomcount, 3);
    long long used = 0;
    const Matrix P = m_density, P2c = P2;  // contiguous
    if (!m_gpu_xc_obj->responseGradient(P.data(), P2c.data(), m_grid_response, gb.data(), gp.data(), gw.data(), used)) {
        gpuDHXCFailed("XC response gradient kernel failed");
        return false;
    }
    g_basis = gb;
    g_response = m_grid_response ? Matrix(gp + gw) : Matrix(Matrix::Zero(m_atomcount, 3));
    if (envOn("CURCUMA_QM_GPU_XC_CHECK")) {
        const double dev_ms = msSinceGpu(t0);
        const auto c0 = std::chrono::steady_clock::now();
        qmxc::XCGradient xc;
        qmxc::buildXCResponseGradient(m_xc_grid, m_xc_ao, m_density, P2, m_xc, m_xc_rho_threshold, m_threads,
                                      activeFunctionAtoms(), m_atomcount, xc, &m_xc_screen);
        Matrix cw = Matrix::Zero(m_atomcount, 3);
        if (m_grid_response)
            qmgrid::gridWeightGradient(m_atoms, atomPositionsBohr(), m_grid_options, m_xc_grid, xc.edens, cw, why, m_threads);
        const double db = (gb - xc.basis).cwiseAbs().maxCoeff();
        const double dp = m_grid_response ? (gp - xc.points).cwiseAbs().maxCoeff() : 0.0;
        const double dw = m_grid_response ? (gw - cw).cwiseAbs().maxCoeff() : 0.0;
        CurcumaLogger::result(fmt::format(
            "QM GPU DH XC grad check: natoms={} max|dbasis|={:.3e} max|dpoints|={:.3e} max|dweights|={:.3e} "
            "max|basis|={:.3e} points_used dev={} cpu={} dev_ms={:.1f} cpu_ms={:.1f}",
            m_atomcount, db, dp, dw, xc.basis.cwiseAbs().maxCoeff(), used, xc.points_used, dev_ms, msSinceGpu(c0)));
    }
    return true;
}

// Claude Generated (Oct 31 2026, WP12 stage 3): PT2 two-electron forces on the device.
// The separable PT2 weight of gradientTwoElectronDHPT2, (P2_mn P_ls + P_mn P2_ls) - a_x/4 (4 exchange products), is
// the bilinear part of the Q-G7 weight Q(X) = X_mn X_ls - a_x/4 (X_ml X_ns + X_ms X_nl): Q(P + P2) - Q(P - P2) = 2 x
// that weight, so g_sep = 1/2 [G(P + P2) - G(P - P2)] with G the Q-G7 gradient (QMGpuGrad2e::compute, kx = a_x). Two
// device passes, no new kernel; the quartet screen of each pass bounds its own weight (max |P +- P2|). Rounding: both
// G are O(|g_KS|), the difference O(|g_PT2|), so the cancellation costs a factor |g_KS|/|g_PT2| (1e1..1e2) in relative
// precision, below the ~1e-13 run-to-run noise of the atomics in Q-G7. The canonical non-separable part 4 Gamma^sym
// runs through the MP2 2e kernel (QMGpuMP2Grad2e, cuda/qm_gpu_mp2_grad.cu) with P = P2 = 0, where its weight is Gamma.
bool QMEngine::dhPT2TwoElectronDevice(const qmint::EriGpuTables* tables, const std::string& why_tables,
                                      const std::string& backend, const Matrix& Pc, const Matrix& P2c, const double* Gamma,
                                      Matrix& g2, std::string& info) const
{
    if (!tables) return false;
    const double kx = m_xc.exact_exchange;
    const Matrix Pp = Pc + P2c, Pm = Pc - P2c;
    Matrix gp, gm;
    std::string ip, im;
    if (!gpuGradientTwoElectron(tables, why_tables, backend, Pp, kx, gp, ip)) return false;
    if (!gpuGradientTwoElectron(tables, why_tables, backend, Pm, kx, gm, im)) return false;
    g2 = 0.5 * (gp - gm);
    info = "separable part 1/2 [G(P + P2) - G(P - P2)] on " + ip;
    if (!Gamma) return true;
    const std::string be = gpuGradBackend(m_gpu_grad, "gpu_grad", "qm_mp2_grad", m_gpu_mp2_grad_failed);
    if (be.empty()) return false;
    std::string why = why_tables;
    qmint::QMGpuGradExtra extra;
    extra.natoms = m_atomcount;
    extra.sh_atom = tables->sh_atom.data();
    extra.nprimpair = (int)tables->pp_ea.size();
    extra.pp_ea = tables->pp_ea.data();
    extra.pp_eb = tables->pp_eb.data();
    std::unique_ptr<qmint::QMGpuMP2Grad2e> dev = gpu_plugin::createQmMp2Grad2e(be, tables->view(), extra, m_gpu_device, why);
    std::unique_ptr<double[]> Gc;
    const double* Gcart = Gamma;
    if (dev && m_Q.size() != 0) {
        qmint::mp2ActiveToCartesian(Gamma, m_Q, Gc, m_threads);
        Gcart = Gc.get();
    }
    const Matrix Z = Matrix::Zero(Pc.rows(), Pc.cols());
    Matrix gG = Matrix::Zero(m_atomcount, 3);
    if (!dev || !dev->compute(Z.data(), Z.data(), Gcart, gG.data())) {
        m_gpu_mp2_grad_failed = true;
        curcuma::reportGpuFallback("qm_mp2_grad", "device double-hybrid Gamma 2e gradient not possible: " +
                                                      (dev ? std::string("kernel failed") : why));
        return false;
    }
    g2 += gG;
    info += "; non-separable Gamma part on " + dev->name();
    return true;
}

// Claude Generated (Oct 31 2026, WP12 stage 3): step 5 of computeDoubleHybridGradient with -gpu (PT2 1e, 2e and
// RI-3-index derivative integrals), each part on the device or, if that is not possible, its CPU function.
void QMEngine::dhPT2ForcesDevice(const Matrix& posBohr, const Matrix& Pc, const Matrix& P2c, const Matrix& Wc,
                                 const double* Gamma, qmri::RIJK* ri, const Matrix& G3, const Eigen::MatrixXd& M2,
                                 bool ks2e_deferred, int zvector_iterations, qmint::MP2Result& res, Matrix& g1,
                                 Matrix& g2) const
{
    const std::string be1 = gpuGradBackend(m_gpu_grad1e, "gpu_grad1e", "qm_grad1e", m_gpu_grad1e_failed);
    // a deferred KS 2e part means the device 2e gradient was not possible in step 1: the CPU forms both in one pass
    const std::string be2 = ks2e_deferred ? std::string() : gpuGradBackend(m_gpu_grad, "gpu_grad", "qm_grad", m_gpu_grad_failed);
    qmint::EriGpuTables tables;
    std::string tab_why;
    bool tables_ok = false;
    if (!be1.empty() || !be2.empty()) tables_ok = qmint::buildEriGradGpuTables(m_gto_basis, m_eri_screening, tables, tab_why);
    const auto t_1e0 = std::chrono::steady_clock::now();
    std::string dev1_info;
    const bool on_device1 = !be1.empty()
        && gpuGradientOneElectron(tables_ok ? &tables : nullptr, tab_why, be1, posBohr, P2c, Wc, g1, dev1_info);
    if (!on_device1) g1 = qmint::gradientOneElectron(m_gto_basis, m_atoms, posBohr, P2c, Wc, m_threads);
    res.t_grad_deriv1e_ms = msSinceGpu(t_1e0);
    const auto t_2e0 = std::chrono::steady_clock::now();
    const Matrix* Qp = m_Q.size() != 0 ? &m_Q : nullptr;
    std::string dev2_info;
    const bool on_device2 = !be2.empty()
        && dhPT2TwoElectronDevice(tables_ok ? &tables : nullptr, tab_why, be2, Pc, P2c, Gamma, g2, dev2_info);
    if (on_device2 && envOn("CURCUMA_QM_GPU_GRAD_CHECK")) {
        const Matrix gc = qmint::gradientTwoElectronDHPT2(m_gto_basis, Pc, P2c, Gamma, Qp, m_xc.exact_exchange, m_atomcount,
                                                          m_threads, false);
        CurcumaLogger::result(fmt::format("QM GPU DH PT2 2e check: natoms={} max|dev-cpu|={:.3e} max|cpu|={:.3e} gamma={}",
                                          m_atomcount, (g2 - gc).cwiseAbs().maxCoeff(), gc.cwiseAbs().maxCoeff(), Gamma ? 1 : 0));
    }
    if (!on_device2)
        g2 = qmint::gradientTwoElectronDHPT2(m_gto_basis, Pc, P2c, Gamma, Qp, m_xc.exact_exchange, m_atomcount, m_threads,
                                             ks2e_deferred);
    bool on_device3 = false;
    if (ri) {
        // RI-PT2 3-index part (the quantity of RIJK::threeIndexGradient); the 2-index metric part stays on the host
        // (N_aux^2 shell pairs, as in the RI-JK gradient)
        Matrix g3;
        const std::string be3 = gpuGradBackend(m_gpu_ri_grad, "gpu_ri_grad", "qm_ri_grad", m_gpu_ri_grad_failed);
        if (!be3.empty()) {
            std::string info3, why3;
            on_device3 = ri->threeIndexGradientDevice(m_gto_basis, Qp, Eigen::VectorXd(), m_density, G3, m_atomcount, be3,
                                                      m_gpu_device, g3, info3, why3);
            if (!on_device3) {
                m_gpu_ri_grad_failed = true;
                curcuma::reportGpuFallback("qm_ri_grad", "device RI-PT2 3-index gradient not possible: " + why3);
            } else if (envOn("CURCUMA_QM_GPU_RI_GRAD_CHECK")) {
                const Matrix gc = ri->threeIndexGradient(m_gto_basis, Qp, Eigen::VectorXd(), m_density, G3, m_atomcount);
                CurcumaLogger::result(fmt::format("QM GPU DH RI grad check: natoms={} max|dev-cpu|={:.3e} max|cpu|={:.3e}",
                                                  m_atomcount, (g3 - gc).cwiseAbs().maxCoeff(), gc.cwiseAbs().maxCoeff()));
            }
        }
        if (!on_device3) g3 = ri->threeIndexGradient(m_gto_basis, Qp, Eigen::VectorXd(), m_density, G3, m_atomcount);
        g2 += g3 - ri->metricGradient(M2, m_atomcount);
    }
    res.t_grad_deriv2e_ms = msSinceGpu(t_2e0);
    if (on_device1) res.info += "; PT2 1e derivative integrals on the device (" + dev1_info + ")";
    if (on_device2) res.info += "; PT2 2e derivative integrals on the device (" + dev2_info + ")";
    if (on_device3) res.info += "; RI-PT2 3-index derivative integrals on the device";
    if (envOn("CURCUMA_QM_SCF_TRACE"))
        fmt::print(stderr, "QM DH gradient trace: zvector_ms={:.1f} ({} it) deriv1e_ms={:.1f} ({}) deriv2e_ms={:.1f} ({}{})\n",
                   res.t_grad_zvector_ms, zvector_iterations, res.t_grad_deriv1e_ms, on_device1 ? "device" : "cpu",
                   res.t_grad_deriv2e_ms, on_device2 ? "device" : "cpu",
                   ri ? (on_device3 ? ", RI 3-index device" : ", RI 3-index cpu") : "");
    if (CurcumaLogger::get_verbosity() >= 3)
        CurcumaLogger::info(fmt::format("double-hybrid gradient (device parts): PT2 |1e| {:.3e} ({}), |2e{}| {:.3e} ({}) Eh/Bohr",
                                        g1.norm(), on_device1 ? "device" : "CPU", ri ? " + RI" : "", g2.norm(),
                                        on_device2 ? "device" : "CPU"));
}
