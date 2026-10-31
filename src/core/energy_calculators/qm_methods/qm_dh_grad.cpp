/*
 * <Analytic nuclear gradient of the closed-shell double hybrids (WP12 stage 2, CPU)>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * E = E_KS[hybrid] + c_os E_os + c_ss E_ss (+ E_disp in DoubleHybridMethod), PT2 with the KS orbitals.
 * Equations: docs/DFT_ROADMAP/WP12_doppelhybride.md, section "Stufe 2"; the structure is the MP2 gradient of
 * qm_mp2_grad.cpp (derivation block there) with the HF orbital Hessian replaced by the Kohn-Sham one,
 *   G4(X) = C^T [4 J(C X C^T) - 2 a_x K(C X C^T) + 4 F_xc[C X C^T]] C,
 * the separable two-electron part with a_x-scaled exchange, and the P2-directional derivative of the XC gradient
 * (qm_xc_fxc.h). Contents:
 *   - QMEngine::dhRIIntermediates          RI-PT2 amplitudes, unrelaxed density, Lagrangian, 3-/2-index densities
 *   - QMEngine::dhG4, dhFxc, dhXCResponseGradient
 *   - QMEngine::computeDoubleHybridGradient KS gradient + Z-vector + PT2 force
 * References: N. C. Handy, H. F. Schaefer, J. Chem. Phys. 81, 5031 (1984) (Z-vector); F. Weigend, M. Haeser, Theor.
 * Chem. Acc. 97, 331 (1997) (RI-MP2 gradient); F. Neese, T. Schwabe, S. Grimme, J. Chem. Phys. 126, 124115 (2007),
 * DOI 10.1063/1.2712433 (double-hybrid gradient: KS orbital Hessian with the XC kernel in the Z-vector equations;
 * DOI checked against Crossref, Oct 31 2026; the derivation here was written out and checked numerically, not copied).
 *
 * Claude Generated (Oct 2026, WP12 stage 2). AI-generated, machine-tested; human production test pending.
 *
 * This program is free software under GPL-3.0
 */

#include "qm_citations.h"
#include "qm_engine.h"
#include "qm_mp2.h"
#include "qm_ri.h"
#include "qm_xc_fxc.h"

#include "src/core/curcuma_logger.h"
#include "src/core/gpu_fallback.h"
#include "src/core/units.h"

#include <fmt/format.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <new>
#include <vector>

namespace {
double msSince(std::chrono::steady_clock::time_point t0)
{
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}

using RowMat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

/// C = A B^T in fixed row blocks under OpenMP (each block one sequential product): independent of the thread count.
void gemmABtBlocked(const RowMat& A, const RowMat& B, RowMat& C, int threads)
{
    C.resize(A.rows(), B.rows());
    const Eigen::Index blk = 64, nb = (A.rows() + blk - 1) / blk;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(std::max(1, threads))
#endif
    for (Eigen::Index b = 0; b < nb; ++b) {
        const Eigen::Index r0 = b * blk, len = std::min(blk, A.rows() - r0);
        C.middleRows(r0, len).noalias() = A.middleRows(r0, len) * B.transpose();
    }
}
/// C = A B in fixed row blocks (as gemmABtBlocked).
void gemmABBlocked(const RowMat& A, const RowMat& B, RowMat& C, int threads)
{
    C.resize(A.rows(), B.cols());
    const Eigen::Index blk = 64, nb = (A.rows() + blk - 1) / blk;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(std::max(1, threads))
#endif
    for (Eigen::Index b = 0; b < nb; ++b) {
        const Eigen::Index r0 = b * blk, len = std::min(blk, A.rows() - r0);
        C.middleRows(r0, len).noalias() = A.middleRows(r0, len) * B;
    }
}
}  // namespace

// Claude Generated (Oct 2026, WP12 stage 2).
std::vector<int> QMEngine::activeFunctionAtoms() const
{
    std::vector<int> bf_atom(m_nbf, 0);
    if (m_Q.size() != 0) {
        for (int a = 0; a < m_nbf; ++a)
            for (int c = 0; c < m_Q.rows(); ++c)
                if (m_Q(c, a) != 0.0) {
                    bf_atom[a] = m_gto_basis[c].atom;
                    break;
                }
    } else {
        for (int a = 0; a < m_nbf; ++a) bf_atom[a] = m_gto_basis[a].atom;
    }
    return bf_atom;
}

// Claude Generated (Oct 2026, WP12 stage 2).
bool QMEngine::dhFxc(const Matrix& D, Matrix& V) const
{
    if (!ensureXCGrid()) return false;
    if (dhFxcDevice(D, V)) return true;  // WP12 stage 3: device F_xc (qm_dh_grad_gpu.cpp)
    return qmxc::buildFxc(m_xc_grid, m_xc_ao, m_density, D, m_xc, m_xc_rho_threshold, m_threads, V, &m_xc_screen);
}

// Claude Generated (Oct 2026, WP12 stage 2): KS orbital Hessian product (WP12 doc, stage 2, item 3).
Matrix QMEngine::dhG4(const Matrix& X_mo, const qmint::ERITensor* eri, bool& ok) const
{
    const Matrix Xao = m_mo * X_mo * m_mo.transpose();
    Matrix J, K;
    if (eri) {
        J = qmint::buildCoulomb(*eri, Xao, m_threads);
        K = qmint::buildExchange(*eri, Xao, m_threads);
    } else {
        buildJKExact(Xao, J, K);
    }
    Matrix Vx;
    if (!dhFxc(Xao, Vx)) {
        ok = false;
        return Matrix::Zero(m_nbf, m_nbf);
    }
    return m_mo.transpose() * (4.0 * J - 2.0 * m_xc.exact_exchange * K + 4.0 * Vx) * m_mo;
}

// Claude Generated (Oct 2026, WP12 stage 2): as xcGradient (CPU path) for the P2-directional derivative.
bool QMEngine::dhXCResponseGradient(const Matrix& P2, Matrix& g_basis, Matrix& g_response, std::string& why) const
{
    g_basis = Matrix::Zero(m_atomcount, 3);
    g_response = Matrix::Zero(m_atomcount, 3);
    if (!ensureXCGrid()) {
        why = "no DFT grid";
        return false;
    }
    if (dhXCResponseGradientDevice(P2, g_basis, g_response, why)) return true;  // WP12 stage 3 (qm_dh_grad_gpu.cpp)
    qmxc::XCGradient xg;
    if (!qmxc::buildXCResponseGradient(m_xc_grid, m_xc_ao, m_density, P2, m_xc, m_xc_rho_threshold, m_threads,
                                       activeFunctionAtoms(), m_atomcount, xg, &m_xc_screen)) {
        why = "no second functional derivative for a kernel of " + m_xc.name;
        return false;
    }
    g_basis = xg.basis;
    if (m_grid_response) {
        if (m_xc_frozen_grid) {
            why = "grid response with a frozen grid (validation hook) is undefined; use -qm.grid_response false";
            return false;
        }
        Matrix gw = Matrix::Zero(m_atomcount, 3);
        if (!qmgrid::gridWeightGradient(m_atoms, atomPositionsBohr(), m_grid_options, m_xc_grid, xg.edens, gw, why,
                                        m_threads))
            return false;
        g_response = xg.points + gw;
    }
    return true;
}

// Claude Generated (Oct 2026, WP12 stage 2): RI-PT2 intermediates (WP12 doc, stage 2, items 2 and 4).
bool QMEngine::dhRIIntermediates(const qmint::MP2Options& opt, int nf, double c_os, double c_ss, qmint::MP2Result& res,
                                 Matrix& Poo, Matrix& Pvv, Matrix& L1, Matrix& L2, Matrix& G3, Eigen::MatrixXd& M2,
                                 std::shared_ptr<qmri::RIJK>& ri, std::string& why)
{
    const int n = m_nbf, n_occ = m_num_electrons / 2, o = n_occ - nf, v = n - n_occ;
    const int threads = std::max(1, m_threads);
    std::string file;
    try {
        file = qmri::resolveAuxBasisFile(opt.ri_aux, "jk");
    } catch (const std::exception& e) {
        why = e.what();
        return false;
    }
    std::string src;
    double setup_ms = 0.0;
    ri = correlationRI(file, src, setup_ms, why);
    if (!ri) return false;
    res.t_ri_setup_ms = setup_ms;
    const int nrows = (int)ri->tensor().rows();
    const size_t n2 = (size_t)n * n, ov = (size_t)o * v;
    const int N = ri->naux();
    // memory: X (rows x n^2), amplitudes t, T~, K (3 ov^2), Gamma/c/d (3 ov naux), G3 (naux x n(n+1)/2)
    const double need_mb = 8.0e-6 * ((double)nrows * n2 + 3.0 * (double)ov * ov + 4.0 * (double)ov * nrows
                                     + (double)N * ri->npair() + (double)N * N);
    if (need_mb + ri->tensorMB() > riMemoryLimitMB()) {
        why = fmt::format("RI-PT2 gradient needs {:.0f} MB besides the {:.0f} MB tensor, above the RI limit ({}); use "
                          "-qm.dh_pt2 canonical or raise -qm.ri_max_memory_mb", need_mb, ri->tensorMB(), riMemoryLimitInfo());
        return false;
    }
    const auto t0 = std::chrono::steady_clock::now();
    // ---- B^R_pq for all MO pairs: X = C^T T^Q C per row, then the metric factor (fitRows) ----
    Matrix X(nrows, (Eigen::Index)n2);
    {
        const Matrix CT = m_mo.transpose();
#ifdef _OPENMP
#pragma omp parallel num_threads(threads)
#endif
        {
            std::vector<double> TQ(n2), Y(n2);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
            for (int Q = 0; Q < nrows; ++Q) {
                ri->unpackRow(Q, TQ.data());
                const Eigen::Map<const Matrix> Tq(TQ.data(), n, n);
                Eigen::Map<Matrix> Yq(Y.data(), n, n);
                Yq.noalias() = CT * Tq;
                Eigen::Map<Matrix>(X.data() + (size_t)Q * n2, n, n).noalias() = Yq * m_mo;  // [p][q]
            }
        }
    }
    ri->fitRows(X);
    const int nfit = (int)X.rows();
    res.n_aux = nfit;
    // Bov(i a, R) = B^R_ia (active i)
    RowMat Bov(ov, nfit);
    for (int R = 0; R < nfit; ++R) {
        const double* x = X.data() + (size_t)R * n2;
        for (int i = 0; i < o; ++i)
            for (int a = 0; a < v; ++a) Bov((size_t)i * v + a, R) = x[(size_t)(nf + i) * n + n_occ + a];
    }
    // ---- (ia|jb), amplitudes, energies ----
    RowMat Km;
    gemmABtBlocked(Bov, Bov, Km, threads);  // ov x ov
    RowMat tm(ov, ov), Tt(ov, ov);
    const double c1 = c_os + c_ss;
    std::vector<double> os_i(o, 0.0), ss_i(o, 0.0);
    const Vector& e = m_energies;
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
#endif
    for (int i = 0; i < o; ++i) {
        double os = 0.0, ss = 0.0;
        for (int a = 0; a < v; ++a)
            for (int j = 0; j < o; ++j)
                for (int b = 0; b < v; ++b) {
                    const size_t r = (size_t)i * v + a, c = (size_t)j * v + b;
                    const double k = Km(r, c), kx = Km((size_t)i * v + b, (size_t)j * v + a);
                    const double inv_d = 1.0 / (e(nf + i) + e(nf + j) - e(n_occ + a) - e(n_occ + b));
                    tm(r, c) = k * inv_d;
                    os += k * k * inv_d;
                    ss += k * (k - kx) * inv_d;
                }
        os_i[i] = os;
        ss_i[i] = ss;
    }
    res.e_os = res.e_ss = 0.0;
    for (int i = 0; i < o; ++i) {
        res.e_os += os_i[i];
        res.e_ss += ss_i[i];
    }
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif
    for (int i = 0; i < o; ++i)
        for (int a = 0; a < v; ++a)
            for (int j = 0; j < o; ++j)
                for (int b = 0; b < v; ++b)
                    Tt((size_t)i * v + a, (size_t)j * v + b) =
                        c1 * tm((size_t)i * v + a, (size_t)j * v + b) - c_ss * tm((size_t)i * v + b, (size_t)j * v + a);
    RowMat().swap(Km);
    // ---- unrelaxed density ----
    {
        const Eigen::Map<const RowMat> Tr(Tt.data(), o, (Eigen::Index)(ov * v)), tr(tm.data(), o, (Eigen::Index)(ov * v));
        Poo = Matrix(-2.0 * (Tr * tr.transpose()));
        std::vector<Eigen::MatrixXd> pv(o);
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
#endif
        for (int i = 0; i < o; ++i) {
            const Eigen::Map<const RowMat> Ti(Tt.data() + (size_t)i * v * ov, v, (Eigen::Index)ov);
            const Eigen::Map<const RowMat> ti(tm.data() + (size_t)i * v * ov, v, (Eigen::Index)ov);
            pv[i] = Ti * ti.transpose();
        }
        Pvv = Matrix::Zero(v, v);
        for (int i = 0; i < o; ++i) Pvv += 2.0 * pv[i];
    }
    RowMat().swap(tm);
    res.t_grad_amplitudes_ms = msSince(t0);
    // ---- Gamma^R_ia = sum_jb T~_ijab B^R_jb ----
    const auto t1 = std::chrono::steady_clock::now();
    RowMat Gam;
    gemmABBlocked(Tt, Bov, Gam, threads);  // ov x nfit
    RowMat().swap(Tt);
    // ---- L1 (n x o) and L2 (n x v): fixed chunks over R ----
    {
        const int nch = std::max(1, std::min(32, nfit));
        std::vector<Matrix> l1(nch, Matrix::Zero(n, o)), l2(nch, Matrix::Zero(n, v));
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
#endif
        for (int ch = 0; ch < nch; ++ch) {
            Eigen::MatrixXd Bpa(n, v), Bip(o, n), Gr(o, v);
            for (int R = (int)((long long)nfit * ch / nch); R < (int)((long long)nfit * (ch + 1) / nch); ++R) {
                const double* x = X.data() + (size_t)R * n2;
                for (int p = 0; p < n; ++p)
                    for (int a = 0; a < v; ++a) Bpa(p, a) = x[(size_t)p * n + n_occ + a];
                for (int i = 0; i < o; ++i)
                    for (int p = 0; p < n; ++p) Bip(i, p) = x[(size_t)(nf + i) * n + p];
                for (int i = 0; i < o; ++i)
                    for (int a = 0; a < v; ++a) Gr(i, a) = Gam((size_t)i * v + a, R);
                l1[ch].noalias() += Bpa * Gr.transpose();
                l2[ch].noalias() += Bip.transpose() * Gr;
            }
        }
        L1 = Matrix::Zero(n, o);
        L2 = Matrix::Zero(n, v);
        for (int ch = 0; ch < nch; ++ch) {
            L1 += 4.0 * l1[ch];
            L2 += 4.0 * l2[ch];
        }
    }
    // ---- raw-basis coefficients c = D F^T Gamma, d = D F^T B_ia; M2 = 2 sym(d c^T); G3 = 2 (H + H^T) ----
    Matrix cR(nfit, (Eigen::Index)ov), dR(nfit, (Eigen::Index)ov);
    for (int R = 0; R < nfit; ++R)
        for (size_t x = 0; x < ov; ++x) {
            cR(R, x) = Gam(x, R);
            dR(R, x) = Bov(x, R);
        }
    RowMat().swap(Gam);
    RowMat().swap(Bov);
    Matrix().swap(X);
    ri->fitTransposeRows(cR);  // N x ov
    ri->fitTransposeRows(dR);
    const Eigen::VectorXd& dn = ri->auxNormalisation();
    for (int P = 0; P < N; ++P) {
        cR.row(P) *= dn(P);
        dR.row(P) *= dn(P);
    }
    {
        const Eigen::MatrixXd gam = Eigen::MatrixXd(dR) * Eigen::MatrixXd(cR).transpose();  // N x N
        M2 = gam + gam.transpose();
    }
    const Matrix Co = m_mo.block(0, nf, n, o), Cv = m_mo.block(0, n_occ, n, v);
    G3.resize(N, ri->npair());
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
#endif
    for (int P = 0; P < N; ++P) {
        const Eigen::Map<const Matrix> cP(cR.data() + (size_t)P * ov, o, v);
        const Eigen::MatrixXd H = Co * cP * Cv.transpose();
        double* y = G3.data() + (size_t)P * G3.cols();
        for (int m = 0; m < n; ++m)
            for (int k = 0; k <= m; ++k) y[(size_t)m * (m + 1) / 2 + k] = 2.0 * (H(m, k) + H(k, m));
    }
    res.t_grad_backtrans_ms = msSince(t1);
    res.info = src + fmt::format("; RI-PT2 gradient on the CPU, {} threads", threads);
    return true;
}

// Claude Generated (Oct 2026, WP12 stage 2). Derivation: WP12 doc, stage 2; MP2 analogue: qm_mp2_grad.cpp.
bool QMEngine::computeDoubleHybridGradient(const qmint::MP2Options& opt, double c_os, double c_ss, qmint::MP2Result& res)
{
    res = qmint::MP2Result();
    res.gradient = true;
    res.ri = opt.ri;
    m_gradient_failed = false;
    m_gradient_error.clear();
    auto fail = [&](const std::string& msg) {
        m_gradient = Matrix::Zero(m_atomcount, 3);
        m_gradient_failed = true;
        m_gradient_error = msg;
        CurcumaLogger::error(msg);
        return false;
    };
    if (m_functional != QMFunctional::DoubleHybrid) return fail("double-hybrid gradient: needs the double-hybrid KS reference");
    if (!m_scf_converged) return fail("double-hybrid gradient: the KS SCF is not converged");
    if (riActive())
        return fail("double-hybrid gradient: needs the exact-integral KS SCF (-qm.ri off); the PT2 response with an "
                    "RI-J/RI-JK SCF is not implemented");
    std::string missing;
    if (!qmxc::fxc::functionalCovered(m_xc, &missing))
        return fail("double-hybrid gradient: no second functional derivative (f_xc) for kernel " + missing);
    // WP12 stage 3 (Oct 31 2026): with -gpu the KS gradient, F_xc, the XC response gradient and the PT2 1e/2e/RI
    // derivatives run on the device (qm_dh_grad_gpu.cpp, each part counts its own fallback); the PT2 intermediates
    // (amplitudes, unrelaxed density, Lagrangian, Z-vector algebra) are host code.
    if (gpuRequested())
        curcuma::reportGpuFallback("qm_dh_grad_pt2", "double-hybrid PT2 intermediates (amplitudes, unrelaxed density, "
                                                     "Lagrangian, Z-vector algebra) run on the CPU");
    const int n = m_nbf, n_occ = m_num_electrons / 2;
    if (m_mo.rows() != n || m_mo.cols() != n || m_energies.size() != n)
        return fail("double-hybrid gradient: no MO coefficients / orbital energies from the SCF");
    int nf = 0;
    if (opt.frozen_core)
        for (int Z : m_atoms) nf += qmint::frozenCoreOrbitals(Z);
    if (nf >= n_occ && n_occ > 0)
        return fail(fmt::format("double-hybrid gradient: frozen core ({} orbitals) leaves no active occupied orbital", nf));
    const int o = n_occ - nf, v = n - n_occ;
    res.n_frozen = nf;
    res.n_occ_active = o;
    res.n_virt = v;
    res.n_basis = n;
    qmcite::mp2(opt.ri ? "ri-mp2" : "mp2", true);  // Z-vector method
    CitationRegistry::cite("dh_gradient");

    // ---- 1. KS gradient of the hybrid (WP8/WP11 path, parts {1e, 2e, E_nn, XC basis, grid response}) ----
    // Oct 31 2026: a CPU KS 2e part is deferred and formed together with the PT2 2e part in step 5 (one pass over the
    // derivative integrals, qmint::gradientTwoElectronDHPT2 include_ks); a device KS 2e part (Q-G7) stays separate.
    bool ks2e_deferred = false;
    if (!computeGradient(true, ks2e_deferred)) return false;
    Matrix gKS = m_gradient;
    std::vector<Matrix> parts = m_gradient_parts;
    while (parts.size() < 5) parts.push_back(Matrix::Zero(m_atomcount, 3));
    auto ksTwoElectronOnly = [&]() {  // the deferred KS 2e part alone (no PT2 term to merge with)
        if (!ks2e_deferred) return;
        Matrix Pc0 = m_density;
        if (m_Q.size() != 0) Pc0 = m_Q * m_density * m_Q.transpose();
        const Matrix g2ks = qmint::gradientTwoElectronScaled(m_gto_basis, Pc0, m_xc.exact_exchange, m_atomcount, m_threads);
        gKS += g2ks;
        parts[1] += g2ks;
        m_gradient = gKS;
        m_gradient_parts = parts;
    };
    if (o == 0 || v == 0) {  // no PT2 term
        ksTwoElectronOnly();
        return true;
    }

    // ---- 2. PT2 intermediates ----
    Matrix Poo, Pvv, L1, L2, G3;
    Eigen::MatrixXd M2;
    std::shared_ptr<qmri::RIJK> ri;
    std::unique_ptr<double[]> GammaBuf;
    qmint::ERITensor local;
    const qmint::ERITensor* zeri = nullptr;
    const Matrix& C = m_mo;
    const Vector& e = m_energies;
    if (opt.ri) {
        std::string why;
        if (!dhRIIntermediates(opt, nf, c_os, c_ss, res, Poo, Pvv, L1, L2, G3, M2, ri, why))
            return fail("double-hybrid gradient: RI-PT2 part not possible: " + why);
    } else {
        const double gamma_mb = std::pow((double)n, 4) * 8.0 / 1.0e6;
        if (gamma_mb > opt.max_memory_mb)
            return fail(fmt::format("double-hybrid gradient: the AO two-particle density ({:.1f} MB) exceeds "
                                    "-qm.mp2_max_memory_mb {:.1f}", gamma_mb, opt.max_memory_mb));
        std::string src;
        const double* eri = mp2HostTensor(opt, res, local, src);
        if (!eri) return fail("double-hybrid gradient: no AO tensor (see the error above)");
        try {
            GammaBuf.reset(new double[(size_t)n * n * n * n]);
        } catch (const std::bad_alloc&) {
            return fail(fmt::format("double-hybrid gradient: cannot allocate the AO two-particle density ({:.1f} MB)", gamma_mb));
        }
        qmint::MP2GradIntermediates im;
        qmint::mp2GradientHost(eri, n, C, e, nf, n_occ, c_os, c_ss, m_threads, opt.batch_mb, GammaBuf.get(), im, res);
        res.info = src + fmt::format("; canonical PT2 gradient on the CPU, {} threads", m_threads);
        Poo = im.Poo;
        Pvv = im.Pvv;
        L1 = 4.0 * (C.transpose() * im.Y1);
        L2 = 4.0 * (C.transpose() * im.Y2);
        zeri = local.empty() ? nullptr : &local;
    }
    if (!std::isfinite(res.e_os) || !std::isfinite(res.e_ss))
        return fail("double-hybrid gradient: non-finite PT2 energy (degenerate HOMO/LUMO?)");

    // ---- 3. Lagrangian and Z-vector with the KS orbital Hessian ----
    const auto t_z0 = std::chrono::steady_clock::now();
    bool ok = true;
    Matrix P2 = Matrix::Zero(n, n);
    P2.block(nf, nf, o, o) = Poo;
    P2.block(n_occ, n_occ, v, v) = Pvv;
    Matrix w = Matrix::Zero(std::max(nf, 1), o);
    for (int I = 0; I < nf; ++I)
        for (int jj = 0; jj < o; ++jj) {
            const int j = nf + jj;
            w(I, jj) = L1(I, jj) / (e(j) - e(I));
            P2(I, j) = P2(j, I) = 0.5 * w(I, jj);
        }
    const Matrix G = dhG4(P2, zeri, ok);
    if (!ok) return fail("double-hybrid gradient: XC kernel response failed (grid or kernel)");
    Matrix Lag(v, n_occ);
    for (int a = 0; a < v; ++a)
        for (int k = 0; k < n_occ; ++k)
            Lag(a, k) = (k >= nf ? L1(n_occ + a, k - nf) : 0.0) - L2(k, a) + G(n_occ + a, k);
    auto toSym = [&](const Matrix& z) {
        Matrix Zs = Matrix::Zero(n, n);
        for (int a = 0; a < v; ++a)
            for (int k = 0; k < n_occ; ++k) Zs(n_occ + a, k) = Zs(k, n_occ + a) = 0.5 * z(a, k);
        return Zs;
    };
    auto applyM = [&](const Matrix& z) {
        const Matrix Gz = dhG4(toSym(z), zeri, ok);
        Matrix r(v, n_occ);
        for (int a = 0; a < v; ++a)
            for (int k = 0; k < n_occ; ++k) r(a, k) = (e(n_occ + a) - e(k)) * z(a, k) + Gz(n_occ + a, k);
        return r;
    };
    Matrix diag(v, n_occ);
    for (int a = 0; a < v; ++a)
        for (int k = 0; k < n_occ; ++k) diag(a, k) = e(n_occ + a) - e(k);
    const Matrix rhs = -Lag;
    Matrix z = rhs.cwiseQuotient(diag);
    Matrix r = rhs - applyM(z);
    Matrix pz = r.cwiseQuotient(diag);
    double rz = r.cwiseProduct(pz).sum();
    const double ztol = 1.0e-10;
    const int zmax = 100;
    int it = 0;
    double rmax = r.cwiseAbs().maxCoeff();
    while (ok && rmax > ztol && it < zmax) {
        const Matrix Ap = applyM(pz);
        const double alpha = rz / pz.cwiseProduct(Ap).sum();
        z += alpha * pz;
        r -= alpha * Ap;
        ++it;
        rmax = r.cwiseAbs().maxCoeff();
        if (CurcumaLogger::get_verbosity() >= 3)
            CurcumaLogger::info(fmt::format("double-hybrid Z-vector iteration {:3d}  max|R| = {:.2e}", it, rmax));
        if (rmax <= ztol) break;
        const Matrix s = r.cwiseQuotient(diag);
        const double rz_new = r.cwiseProduct(s).sum();
        pz = s + (rz_new / rz) * pz;
        rz = rz_new;
    }
    if (!ok) return fail("double-hybrid gradient: XC kernel response failed in the Z-vector iterations");
    res.zvector_iterations = it;
    res.zvector_residual = rmax;
    if (rmax > 1.0e-6)
        return fail(fmt::format("double-hybrid gradient: Z-vector equations not converged ({} iterations, max residual {:.2e})",
                                it, rmax));
    if (rmax > ztol && CurcumaLogger::get_verbosity() >= 1)
        CurcumaLogger::warn(fmt::format("double-hybrid gradient: Z-vector stopped at max residual {:.2e} after {} iterations",
                                        rmax, it));
    res.t_grad_zvector_ms = msSince(t_z0);

    // ---- 4. relaxed density and energy-weighted density (as MP2, KS G4) ----
    const auto t_d0 = std::chrono::steady_clock::now();
    P2 += toSym(z);
    const Matrix Gt = dhG4(P2, zeri, ok);
    if (!ok) return fail("double-hybrid gradient: XC kernel response failed (relaxed density)");
    Matrix Ws = Matrix::Zero(n, n);
    for (int k = 0; k < n_occ; ++k)
        for (int l = 0; l < n_occ; ++l) Ws(k, l) = -0.5 * Gt(k, l);
    for (int ii = 0; ii < o; ++ii)
        for (int jj = 0; jj < o; ++jj) {
            const int i = nf + ii, j = nf + jj;
            const double yij = L1(i, jj) + 2.0 * P2(i, j) * e(i);
            const double yji = L1(j, ii) + 2.0 * P2(j, i) * e(j);
            Ws(i, j) += -0.25 * (yij + yji);
        }
    for (int aa = 0; aa < v; ++aa)
        for (int bb = 0; bb < v; ++bb) {
            const int a = n_occ + aa, b = n_occ + bb;
            const double yab = L2(a, bb) + 2.0 * P2(a, b) * e(a);
            const double yba = L2(b, aa) + 2.0 * P2(b, a) * e(b);
            Ws(a, b) += -0.25 * (yab + yba);
        }
    for (int aa = 0; aa < v; ++aa)
        for (int k = 0; k < n_occ; ++k) {
            const double x = 0.5 * (-L2(k, aa) - z(aa, k) * e(k));
            Ws(n_occ + aa, k) += x;
            Ws(k, n_occ + aa) += x;
        }
    for (int I = 0; I < nf; ++I)
        for (int jj = 0; jj < o; ++jj) {
            const double x = -0.5 * w(I, jj) * e(nf + jj);
            Ws(I, nf + jj) += x;
            Ws(nf + jj, I) += x;
        }
    const Matrix P2ao = C * P2 * C.transpose();
    const Matrix Wpt = -(C * Ws * C.transpose());  // gradientOneElectron subtracts W dS: -(-Ws) dS = + Ws dS
    Matrix Pc = m_density, P2c = P2ao, Wc = Wpt;
    if (m_Q.size() != 0) {
        Pc = m_Q * m_density * m_Q.transpose();
        P2c = m_Q * P2ao * m_Q.transpose();
        Wc = m_Q * Wpt * m_Q.transpose();
    }
    if (GammaBuf) {  // non-separable canonical part D_ns = 4 Gamma^sym (8-fold average / 2, in place; as MP2)
        double* Gm = GammaBuf.get();
        const size_t Nn = (size_t)n;
        auto at = [&](size_t a, size_t b, size_t c, size_t d) -> double& { return Gm[((a * Nn + b) * Nn + c) * Nn + d]; };
#ifdef _OPENMP
#pragma omp parallel for num_threads(std::max(1, m_threads)) schedule(dynamic, 1)
#endif
        for (int m = 0; m < n; ++m)
            for (int nn = 0; nn <= m; ++nn) {
                const size_t mn = (size_t)m * (m + 1) / 2 + nn;
                for (int l = 0; l < n; ++l)
                    for (int s = 0; s <= l; ++s) {
                        const size_t ls = (size_t)l * (l + 1) / 2 + s;
                        if (ls > mn) break;
                        double& v0 = at(m, nn, l, s); double& v1 = at(nn, m, l, s);
                        double& v2 = at(m, nn, s, l); double& v3 = at(nn, m, s, l);
                        double& v4 = at(l, s, m, nn); double& v5 = at(s, l, m, nn);
                        double& v6 = at(l, s, nn, m); double& v7 = at(s, l, nn, m);
                        const double sum = 0.5 * (v0 + v1 + v2 + v3 + v4 + v5 + v6 + v7);
                        v0 = v1 = v2 = v3 = v4 = v5 = v6 = v7 = sum;
                    }
            }
    }
    res.t_grad_density_ms = msSince(t_d0);

    // ---- 5. PT2 forces ----
    Matrix posBohr = atomPositionsBohr();
    if (gpuRequested()) {  // WP12 stage 3: device 1e/2e/RI-3-index parts with per-part CPU fallback (qm_dh_grad_gpu.cpp)
        Matrix g1, g2;
        dhPT2ForcesDevice(posBohr, Pc, P2c, Wc, GammaBuf ? GammaBuf.get() : nullptr, ri.get(), G3, M2, ks2e_deferred, it, res,
                          g1, g2);
        Matrix gxb, gxr;
        std::string why;
        if (!dhXCResponseGradient(P2ao, gxb, gxr, why)) return fail("double-hybrid gradient: XC part failed: " + why);
        m_gradient = gKS + g1 + g2 + gxb + gxr;
        m_gradient_parts = { parts[0] + g1, parts[1] + g2, parts[2], parts[3] + gxb, parts[4] + gxr };
        return true;
    }
    const auto t_1e0 = std::chrono::steady_clock::now();
    const Matrix g1 = qmint::gradientOneElectron(m_gto_basis, m_atoms, posBohr, P2c, Wc, m_threads);
    res.t_grad_deriv1e_ms = msSince(t_1e0);
    const auto t_2e0 = std::chrono::steady_clock::now();
    const Matrix* Qp = m_Q.size() != 0 ? &m_Q : nullptr;
    // One pass over the 4-index derivative integrals: PT2 part (separable P2 x P, canonical non-separable Gamma) plus,
    // when deferred in step 1, the KS part P x P (same quartets, same kernel). The RI-PT2 3-index and metric
    // derivatives below are other integrals ((mn|P)^X, (P|Q)^X) and cannot join this pass; the KS SCF of a double
    // hybrid is exact (riActive() is an error above), so there is no KS 3-index pass to merge them with.
    Matrix g2 = qmint::gradientTwoElectronDHPT2(m_gto_basis, Pc, P2c, GammaBuf ? GammaBuf.get() : nullptr, Qp,
                                                m_xc.exact_exchange, m_atomcount, m_threads, ks2e_deferred);
    if (ri) g2 += ri->threeIndexGradient(m_gto_basis, Qp, Eigen::VectorXd(), m_density, G3, m_atomcount)
                - ri->metricGradient(M2, m_atomcount);
    res.t_grad_deriv2e_ms = msSince(t_2e0);
    Matrix gxb, gxr;
    std::string why;
    if (!dhXCResponseGradient(P2ao, gxb, gxr, why)) return fail("double-hybrid gradient: XC part failed: " + why);
    m_gradient = gKS + g1 + g2 + gxb + gxr;
    m_gradient_parts = { parts[0] + g1, parts[1] + g2, parts[2], parts[3] + gxb, parts[4] + gxr };
    if (CurcumaLogger::get_verbosity() >= 3)
        CurcumaLogger::info(fmt::format("double-hybrid gradient: PT2 parts |1e| {:.3e}, |2e{}{}| {:.3e}, |XC basis| {:.3e}, "
                                        "|XC grid response| {:.3e} Eh/Bohr; Z-vector {} iterations",
                                        g1.norm(), ks2e_deferred ? " incl. KS 2e (one pass)" : "", ri ? " + RI" : "",
                                        g2.norm(), gxb.norm(), gxr.norm(), it));
    return true;
}
