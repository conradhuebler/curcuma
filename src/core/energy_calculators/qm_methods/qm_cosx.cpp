/*
 * <Native QM: seminumerical exchange (chain-of-spheres, COSX / RIJCOSX) -- implementation>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Sections:
 *   1. SeminumericalExchange::setup    grid, AO blocks, overlap fit S~^-1 S, potential integrals
 *   2. SeminumericalExchange::exchange K(P) = sym( sum_g phi(g) [A(g) F(g)]^T ), F = w phi P'
 *   3. QMEngine glue                   -qm.exchange_method / -qm.cosx_*, refusals, buildJKCosx
 * Equations and literature: qm_cosx.h; measurements: docs/QM_COSX.md.
 *
 * Claude Generated (Nov 2026). AI-generated, machine-tested only.
 *
 * This program is free software under GPL-3.0
 */

#include "qm_cosx.h"
#include "qm_engine.h"
#include "qm_citations.h"
#include "qm_ri.h"

#include "src/core/curcuma_logger.h"
#include "src/core/gpu_fallback.h"

#include <Eigen/Dense>
#include <fmt/format.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <limits>

namespace qmcosx {

namespace {
double msSince(std::chrono::steady_clock::time_point t0)
{
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}
/// Number of fixed summation chunks over the blocks: independent of the thread count (bitwise 1 vs N threads).
int chunkCount(int nblocks) { return std::max(1, std::min(nblocks, 16)); }
}  // namespace

SeminumericalExchange::SeminumericalExchange() = default;
SeminumericalExchange::~SeminumericalExchange() = default;

// ---------------------------------------------------------------------------------------------
// 1. Setup (once per geometry)
// ---------------------------------------------------------------------------------------------
bool SeminumericalExchange::setup(const std::vector<GTO::Orbital>& basis, const Matrix& Q, const Matrix& S,
                                  const std::vector<int>& Z, const Matrix& posBohr, const qmgrid::GridOptions& gopt,
                                  const CosxOptions& opt, int threads, std::string& why)
{
    const auto t0 = std::chrono::steady_clock::now();
    m_ready = false;
    m_opt = opt;
    qmgrid::GridOptions g = gopt;
    g.level = opt.grid;
    if (!qmgrid::buildMolecularGrid(Z, posBohr, g, m_grid, why, threads)) {
        why = "COSX grid '" + opt.grid + "': " + why;
        return false;
    }
    m_ao = qmgrid::AOEvaluator(basis, Q);
    m_n = m_ao.nbf();
    if (m_n != (int)S.rows()) {
        why = fmt::format("COSX: AO evaluator has {} functions, overlap {}", m_n, S.rows());
        return false;
    }
    qmgrid::sortGridBlocks(m_grid, std::max(16, opt.block), m_begin, m_count);
    m_screen = qmgrid::BlockScreen();
    if (opt.ao_threshold > 0.0) {
        std::string w2;
        if (!qmgrid::buildBlockScreen(m_grid, m_ao, opt.ao_threshold, m_begin, m_count, m_screen, w2)) {
            m_screen = qmgrid::BlockScreen();
            if (CurcumaLogger::get_verbosity() >= 1) CurcumaLogger::warn("QM COSX: " + w2 + "; no AO screening");
        }
    }
    // Overlap fit (Neese 2009 sec. 2.3): proj = S~^-1 S with the numerical overlap of this grid
    m_proj = Matrix::Identity(m_n, m_n);
    m_fit_min_ev = 0.0;
    if (opt.fit_overlap) {
        const int nb = blocks(), nc = chunkCount(nb);
        std::vector<Eigen::MatrixXd> Sc(nc);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads > 0 ? threads : 1)
#endif
        for (int c = 0; c < nc; ++c) {
            Sc[c] = Eigen::MatrixXd::Zero(m_n, m_n);
            qmgrid::AOBlock blk;
            std::vector<int> cols;
            for (int b = (int)((long long)c * nb / nc); b < (int)((long long)(c + 1) * nb / nc); ++b) {
                const int begin = m_begin[b], count = m_count[b];
                if (m_screen.active()) {
                    const int off = m_screen.shell_off[b], nsh = m_screen.shell_off[b + 1] - off;
                    if (nsh == 0) continue;
                    m_ao.evaluateShells(m_grid.points, begin, count, 0, m_screen.shells.data() + off, nsh, blk, cols);
                } else {
                    m_ao.evaluate(m_grid.points, begin, count, 0, blk);
                    cols.resize(m_n);
                    for (int i = 0; i < m_n; ++i) cols[i] = i;
                }
                const Eigen::MatrixXd phi(blk.phi);
                Eigen::MatrixXd wphi = phi;
                for (int p = 0; p < count; ++p) wphi.row(p) *= m_grid.weights(begin + p);
                const Eigen::MatrixXd s = phi.transpose() * wphi;
                for (int i = 0; i < (int)cols.size(); ++i)
                    for (int j = 0; j < (int)cols.size(); ++j) Sc[c](cols[i], cols[j]) += s(i, j);
            }
        }
        Eigen::MatrixXd St = Eigen::MatrixXd::Zero(m_n, m_n);
        for (const auto& s : Sc) St += s;
        const Eigen::MatrixXd Sa(S);
        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> ges(St, Sa, Eigen::EigenvaluesOnly);
        m_fit_min_ev = ges.info() == Eigen::Success ? ges.eigenvalues()(0) : 0.0;
        if (!(m_fit_min_ev > 1.0e-8)) {
            why = fmt::format("COSX: numerical overlap of grid '{}' is (near) singular in the basis (smallest eigenvalue of "
                              "S~ relative to S {:.3e}); use a larger -qm.cosx_grid or -qm.cosx_fit_overlap false",
                              opt.grid, m_fit_min_ev);
            return false;
        }
        m_proj = St.partialPivLu().solve(Sa);
    }
    m_pot = qmint::GridPotentialIntegrals(basis, Q, opt.threshold > 0.0 ? opt.threshold * 1.0e-6 : 0.0);
    m_setup_ms = msSince(t0);
    m_build_ms = 0.0;
    m_builds = 0;
    m_ready = true;
    return true;
}

// ---------------------------------------------------------------------------------------------
// 2. K(P)
//   P' = S~^-1 S P;  per block: F = (w phi) P' (points x n), G = A F (G_gm = sum_l A_ml(g) F_gl, analytic),
//   K(cols, :) += phi^T G;  K = (K + K^T) / 2.
// ---------------------------------------------------------------------------------------------
Matrix SeminumericalExchange::exchange(const Matrix& P, int threads) const
{
    const auto t0 = std::chrono::steady_clock::now();
    const int n = m_n, nb = blocks(), nc = chunkCount(nb);
    const Eigen::MatrixXd Pp = m_opt.fit_overlap ? Eigen::MatrixXd(Eigen::MatrixXd(m_proj) * Eigen::MatrixXd(P)) : Eigen::MatrixXd(P);
    std::vector<Eigen::MatrixXd> Kc(nc);
    std::vector<long long> evals(nc, 0);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads > 0 ? threads : 1)
#endif
    for (int c = 0; c < nc; ++c) {
        Kc[c] = Eigen::MatrixXd::Zero(n, n);
        qmgrid::AOBlock blk;
        std::vector<int> cols;
        Matrix F, G;  // row-major points x n (the layout GridPotentialIntegrals::contract reads)
        for (int b = (int)((long long)c * nb / nc); b < (int)((long long)(c + 1) * nb / nc); ++b) {
            const int begin = m_begin[b], count = m_count[b];
            if (m_screen.active()) {
                const int off = m_screen.shell_off[b], nsh = m_screen.shell_off[b + 1] - off;
                if (nsh == 0) continue;
                m_ao.evaluateShells(m_grid.points, begin, count, 0, m_screen.shells.data() + off, nsh, blk, cols);
            } else {
                m_ao.evaluate(m_grid.points, begin, count, 0, blk);
                cols.resize(n);
                for (int i = 0; i < n; ++i) cols[i] = i;
            }
            const int nk = (int)cols.size();
            Eigen::MatrixXd wphi(blk.phi);
            for (int p = 0; p < count; ++p) wphi.row(p) *= m_grid.weights(begin + p);
            Eigen::MatrixXd Prow(nk, n);
            for (int i = 0; i < nk; ++i) Prow.row(i) = Pp.row(cols[i]);
            F = wphi * Prow;
            G = Matrix::Zero(count, n);
            evals[c] += m_pot.contract(m_grid.points.data() + 3 * (size_t)begin, count, F.data(), n, G.data(), n,
                                       m_opt.threshold);
            const Eigen::MatrixXd phi(blk.phi);
            const Eigen::MatrixXd kb = phi.transpose() * Eigen::MatrixXd(G);
            for (int i = 0; i < nk; ++i) Kc[c].row(cols[i]) += kb.row(i);
        }
    }
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(n, n);
    long long ev = 0;
    for (int c = 0; c < nc; ++c) {
        K += Kc[c];
        ev += evals[c];
    }
    const Eigen::MatrixXd Ks = 0.5 * (K + K.transpose());
    m_last_fraction = (double)ev / std::max(1.0, (double)m_pot.npairs() * (double)m_grid.size());
    m_build_ms += msSince(t0);
    ++m_builds;
    return Matrix(Ks);
}

}  // namespace qmcosx

// =============================================================================================
// 3. QMEngine glue
// =============================================================================================

namespace {
std::string lowerCopy(std::string s)
{
    for (char& c : s) c = (char)std::tolower((unsigned char)c);
    return s;
}
}  // namespace

void QMEngine::readCosxConfig(const json& config)
{
    if (config.contains("exchange_method") && config["exchange_method"].is_string())
        m_exchange_method = lowerCopy(config["exchange_method"].get<std::string>());
    if (m_exchange_method == "seminumerical" || m_exchange_method == "sgx") m_exchange_method = "cosx";
    if (m_exchange_method != "exact" && m_exchange_method != "cosx") {
        // an unknown name must not run silently with the exact exchange under another label: refuse in Calculation
        m_cosx_error = "unknown -qm.exchange_method '" + m_exchange_method + "' (exact | cosx; RI-K is -qm.ri jk)";
    }
    if (config.contains("cosx_grid") && config["cosx_grid"].is_string()) m_cosx_opt.grid = config["cosx_grid"].get<std::string>();
    if (config.contains("cosx_fit_overlap") && config["cosx_fit_overlap"].is_boolean())
        m_cosx_opt.fit_overlap = config["cosx_fit_overlap"].get<bool>();
    if (config.contains("cosx_threshold") && config["cosx_threshold"].is_number())
        m_cosx_opt.threshold = std::max(0.0, config["cosx_threshold"].get<double>());
    m_cosx_opt.ao_threshold = m_xc_ao_threshold;
}

bool QMEngine::cosxRequested() const { return m_exchange_method == "cosx"; }

bool QMEngine::cosxActive() const
{
    if (!cosxRequested()) return false;
    return !isKohnSham() || m_xc.exact_exchange != 0.0;  // a pure functional builds no exchange
}

// Refusals before any SCF (Claude Generated, Nov 2026). Returns false with m_cosx_error set.
bool QMEngine::cosxValidate()
{
    if (!m_cosx_error.empty() && !cosxRequested()) return false;  // unknown name
    if (!cosxRequested()) return true;
    m_cosx_error.clear();
    if (!cosxActive()) {
        static bool warned = false;
        if (!warned && CurcumaLogger::get_verbosity() >= 1) {
            warned = true;
            CurcumaLogger::warn("QM: -qm.exchange_method cosx ignored -- " + m_xc.name + " has no exact exchange");
        }
        return true;
    }
    auto refuse = [&](const std::string& why) {
        m_cosx_error = "-qm.exchange_method cosx: " + why;
        return false;
    };
    if (m_open_shell) return refuse("closed shell only (open-shell UHF/UKS is not implemented)");
    if (riActive() && m_ri_mode == "jk") return refuse("-qm.ri jk fits the exchange already; use -qm.ri j (RIJCOSX) or -qm.ri off");
    if (riActive() && m_ri_mixed_precision == "true") return refuse("-qm.ri_mixed_precision is not supported");
    if (hasECP()) return refuse("effective core potentials are not validated with COSX (use an all-electron basis)");
    if (m_scf_stability != "off") return refuse("-qm.scf_stability is not implemented for the COSX operator");
    if (gpuRequested()) curcuma::reportGpuFallback("qm_cosx", "COSX exchange and its J build run on the CPU");
    return true;
}

void QMEngine::invalidateCosx() const
{
    m_cosx.reset();
    m_cosx_failed = false;
}

bool QMEngine::ensureCosx() const
{
    if (m_cosx && m_cosx->ready()) return true;
    if (m_cosx_failed) return false;
    auto c = std::make_shared<qmcosx::SeminumericalExchange>();
    std::string why;
    qmcosx::CosxOptions opt = m_cosx_opt;
    opt.ao_threshold = m_xc_ao_threshold;
    if (!c->setup(m_gto_basis, m_Q, m_S, m_atoms, atomPositionsBohr(), m_grid_options, opt, m_threads, why)) {
        m_cosx_failed = true;
        CurcumaLogger::error("QM: " + why);
        return false;
    }
    qmcite::cosx();
    if (CurcumaLogger::get_verbosity() >= 2)
        CurcumaLogger::info(fmt::format("QM: COSX exchange grid {} ({} points, {} blocks), overlap fit {} (min eigenvalue "
                                        "S~/S {:.6f}), {} shell pairs, setup {:.1f} ms",
                                        opt.grid, c->gridPoints(), c->blocks(), opt.fit_overlap ? "on" : "off",
                                        c->fitMinEigenvalue(), c->pairs(), c->setupMs()));
    if (const char* dump = std::getenv("CURCUMA_QM_COSX_DUMP"); dump && dump[0] != '\0') {
        // grid points (Bohr) and weights, one line per point: diagnostic for the external reference
        std::ofstream f(dump);
        f.precision(17);
        const auto& g = c->grid();
        for (int p = 0; p < g.size(); ++p) f << g.points(p, 0) << " " << g.points(p, 1) << " " << g.points(p, 2) << " " << g.weights(p) << "\n";
    }
    m_cosx = c;
    return true;
}

// J from RI-J (-qm.ri j) or the exact CPU path (stored tensor / integral-direct), K from COSX. CPU only.
void QMEngine::buildJKCosx(const Matrix& P, Matrix& J, Matrix& K, bool need_k) const
{
    const double nan = std::numeric_limits<double>::quiet_NaN();
    if (riActive()) {
        if (!ensureRI()) {
            J = Matrix::Constant(m_nbf, m_nbf, nan);
            K = J;
            return;
        }
        const auto t0 = std::chrono::steady_clock::now();
        J = m_ri->coulomb(P);
        m_ri_build_ms += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
        ++m_ri_builds;
    } else if (!useDirectSCF()) {
        J = qmint::buildCoulomb(eriActive(), P, m_threads);
    } else {
        // integral-direct: the CPU digest forms J and K of the same quartets; K is discarded (no J-only kernel)
        if (!m_direct.ready()) m_direct = qmint::DirectJK(m_gto_basis, m_eri_screening, m_Q.size() != 0 ? &m_Q : nullptr);
        Matrix Kx;
        m_direct_quartets += m_direct.build(P, J, Kx, m_threads, m_direct_thr_now, m_eri_prim_density);
    }
    if (!need_k) {
        K = Matrix::Zero(m_nbf, m_nbf);
        return;
    }
    if (!ensureCosx()) {
        K = Matrix::Constant(m_nbf, m_nbf, nan);
        return;
    }
    K = m_cosx->exchange(P, m_threads);
}

void QMEngine::cosxReport() const
{
    if (!m_cosx || CurcumaLogger::get_verbosity() < 2) return;
    CurcumaLogger::info(fmt::format("QM: COSX {} K builds, {:.1f} ms ({:.1f} ms per build), last build evaluated {:.1f} % of "
                                    "the (shell pair, point) combinations",
                                    m_cosx->builds(), m_cosx->buildMs(), m_cosx->buildMs() / std::max<long>(1, m_cosx->builds()),
                                    100.0 * m_cosx->lastPairFraction()));
}
