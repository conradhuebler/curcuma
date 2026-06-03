/*
 * <Native xTB sparse + non-orthogonal purification SCF — implementation>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0.
 *
 * Claude Generated (June 2026). Implements XTB::calculateSparsePurification
 * (declared in xtb_native.h). It replaces the O(N^3) generalized eigensolve with
 * non-orthogonal canonical density-matrix purification:
 *
 *   Palser & Manolopoulos, PRB 58 (1998) 12704 — canonical (trace-conserving)
 *   purification. We purify the projector  M = P*S  (a plain nao x nao matrix,
 *   similar to the orthonormal idempotent S^{1/2} P S^{1/2}, so its eigenvalues
 *   are the occupations and Tr(M) = Nocc). The canonical PM map on M uses only
 *   plain matrix powers (M^2 = M*M) — the S metric enters only when forming the
 *   initial M0 = lambda*(mu*I - S^-1 F) + (Nocc/nao) I and when recovering the
 *   AO density D = 2 * M * S^-1. M is thresholded each step (drop |M_ij| <
 *   threshold) to expose the sparsity a gapped/insulating system supports
 *   (density-matrix nearsightedness, Kohn 1996).
 *
 * 0 K (integer occupation, inherent to purification), GAPPED systems only. If the
 * purification diverges (gapless / no HOMO-LUMO gap) it warns once and falls back
 * to the validated eigensolver for the rest of the SCF (energy stays correct;
 * the sparsity benefit is simply unavailable).
 *
 * First cut: dense Eigen storage + thresholding. It MEASURES the energy error
 * and the achievable nnz fraction vs the threshold (the validation contract);
 * true sparse storage / S^-1-free O(N) work is deferred. Energy-only.
 */

#include "xtb_native.h"
#include "src/core/curcuma_logger.h"

#include <Eigen/Dense>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fmt/format.h>

namespace curcuma::xtb {

namespace {
// Count of |A_ij| > tol, as a fraction of nao^2.
double nnzFraction(const Matrix& A, double tol)
{
    const long n = static_cast<long>(A.size());
    if (n == 0) return 0.0;
    long nz = 0;
    const double* p = A.data();
    for (long i = 0; i < n; ++i)
        if (std::abs(p[i]) > tol) ++nz;
    return static_cast<double>(nz) / static_cast<double>(n);
}

// Tight spectral bounds [emin,emax] of a SYMMETRIC matrix A via power iteration.
// Gershgorin is a safe outer bracket but grows far too loose with system size
// (wide radius -> the PM slope lambda collapses -> M0 degenerates -> purification
// stalls and is misread as "gapless"). Power iteration recovers the true extremal
// eigenvalues in O(n^2) per step: first the largest-magnitude eigenvalue, then
// (shifted by it) the one farthest away = the opposite extreme. The estimates are
// padded outward and clamped to the Gershgorin bracket so M0 always covers the
// true spectrum. Claude Generated.
void spectralBounds(const Eigen::MatrixXd& A, double& emin, double& emax)
{
    const int n = static_cast<int>(A.rows());
    double gmin = 1.0e30, gmax = -1.0e30;
    for (int i = 0; i < n; ++i) {
        double r = 0.0;
        for (int j = 0; j < n; ++j) if (j != i) r += std::abs(A(i, j));
        gmin = std::min(gmin, A(i, i) - r);
        gmax = std::max(gmax, A(i, i) + r);
    }
    auto rayleigh_extreme = [&](const Eigen::MatrixXd& B) {
        Eigen::VectorXd v = Eigen::VectorXd::Random(n).normalized();
        double lam = 0.0;
        for (int it = 0; it < 60; ++it) {
            Eigen::VectorXd w = B * v;
            const double nw = w.norm();
            if (nw < 1.0e-300) break;
            lam = v.dot(w);          // Rayleigh quotient
            v = w / nw;
        }
        return lam;                  // ~ eigenvalue of largest magnitude of B
    };
    const double lam1 = rayleigh_extreme(A);                 // largest-|lambda|
    // Shifted matrix A - lam1*I: its largest-|.| eigenvalue is the spectrum point
    // farthest from lam1, i.e. the opposite extreme of A.
    const Eigen::MatrixXd Ashift = A - lam1 * Eigen::MatrixXd::Identity(n, n);
    const double lam2 = lam1 + rayleigh_extreme(Ashift);
    emax = std::max(lam1, lam2);
    emin = std::min(lam1, lam2);
    const double pad = 0.05 * (emax - emin) + 1.0e-3;
    emin = std::max(gmin, emin - pad);
    emax = std::min(gmax, emax + pad);
}
} // namespace

/* -------------------------------------------------------------------------- *
 *  Non-orthogonal canonical PM purification of the projector M = P*S.
 *
 *  M0 is supplied with eigenvalues in [0,1] and Tr(M0)=Nocc. On success returns
 *  true and the idempotent M (Tr(M)=Nocc, M^2=M). Returns false if it fails to
 *  converge in max_steps (gapless symptom). NB: we purify WITHOUT in-loop
 *  thresholding — pruning a not-yet-idempotent M corrupts convergence and trips
 *  the watchdog. Sparsity is applied as a one-shot prune of the converged M by
 *  the caller, which both measures the achievable nnz and lets the pruned
 *  density feed the SCF (so the charges feel the sparsity).
 * -------------------------------------------------------------------------- */
static bool purifyProjector(Eigen::MatrixXd& M, int Nocc, int max_steps = 200)
{
    double best_err = 1.0e30;
    int no_improve = 0;
    for (int step = 0; step < max_steps; ++step) {
        const Eigen::MatrixXd M2 = M * M;
        const double trM  = M.trace();
        const double trM2 = M2.trace();
        const double idem_err = std::abs(trM - trM2);   // Tr(M-M^2) -> 0 at idempotency
        if (idem_err < 1.0e-11) return true;

        // Divergence / stall watchdog. Canonical PM is monotone in exact
        // arithmetic but thresholding adds noise, so only give up if the best
        // idempotency error has not improved for many consecutive steps
        // (gapless systems plateau well above zero).
        if (idem_err < best_err - 1.0e-13) { best_err = idem_err; no_improve = 0; }
        else if (++no_improve > 30) return false;

        const Eigen::MatrixXd M3 = M2 * M;
        const double trM3 = M3.trace();
        const double denom = trM - trM2;
        if (std::abs(denom) < 1.0e-14) return true;       // already idempotent
        const double cn = (trM2 - trM3) / denom;
        if (!std::isfinite(cn)) return false;

        if (cn <= 0.5)
            M = ((1.0 - 2.0 * cn) * M + (1.0 + cn) * M2 - M3) / (1.0 - cn);
        else
            M = ((1.0 + cn) * M2 - M3) / cn;
    }
    (void)Nocc;
    return false;
}

// Prune |M_ij| < threshold in place; return the surviving fraction (nnz/n^2).
static double pruneFraction(Eigen::MatrixXd& M, double threshold)
{
    const long n = static_cast<long>(M.size());
    if (n == 0) return 0.0;
    double* p = M.data();
    long nz = 0;
    for (long i = 0; i < n; ++i) {
        if (std::abs(p[i]) < threshold) p[i] = 0.0;
        else ++nz;
    }
    return static_cast<double>(nz) / static_cast<double>(n);
}

double XTB::calculateSparsePurification(double threshold,
                                        bool& converged_out, int& iters_out,
                                        double& nnz_frac_out)
{
    using clock = std::chrono::steady_clock;
    const int verb = CurcumaLogger::get_verbosity();
    const int nat  = m_atomcount;
    const int nsh  = m_basis.nsh;
    const int nao  = m_basis.nao;
    converged_out = false;
    iters_out = 0;
    nnz_frac_out = 1.0;

    const int Nocc = static_cast<int>(std::lround(m_wfn.nocc / 2.0));
    if (std::abs(m_wfn.nocc - 2.0 * Nocc) > 1.0e-6) {
        CurcumaLogger::warn("large_system_mode=sparse: purification needs a closed-shell (even-electron) "
                            "system; falling back to dense.");
        return 0.0;   // driver falls back to dense
    }

    // ---- build-once global setup (mirror Calculation steps 1-5) ----
    m_d4_prepared = false;
    m_d4_genparams_calls = 0;
    const auto t0 = clock::now();
    Vector cn = computeCoordinationNumbers();
    Vector se;
    getSelfEnergies(cn, se);
    Matrix S, H0;
    getHamiltonianH0(se, S, H0);
    m_S  = S;
    m_H0 = H0;
    buildOrthonormalizer();              // m_X = chol(S); for the gapless fallback solveEigen
    buildGammaMatrix();
    if (m_method == MethodType::GFN2)
        setupMultipole();
    m_coordination_numbers = cn;

    // S^-1 once (dense; first-cut — the S^-1-free O(N) route is deferred).
    const Eigen::MatrixXd Sinv = m_S.inverse();

    // ---- initial charge guess (EEQ) ----
    auto rebuildQatFromQsh = [&]() {
        m_wfn.q_at.setZero(nat);
        for (int s = 0; s < nsh; ++s) m_wfn.q_at(m_basis.sh2at[s]) += m_wfn.q_sh(s);
    };
    Vector q_sh_guess;
    m_wfn.q_sh = seedEEQGuess(q_sh_guess) ? q_sh_guess : Vector::Zero(nsh);
    rebuildQatFromQsh();
    m_wfn.dp_at = Eigen::MatrixXd::Zero(3, nat);
    m_wfn.qp_at = Eigen::MatrixXd::Zero(6, nat);

    auto packSCC = [&]() -> Vector {
        if (m_method == MethodType::GFN2) {
            Vector v = Vector::Zero(nsh + 9 * nat);
            v.head(nsh) = m_wfn.q_sh;
            for (int k = 0; k < 3; ++k) v.segment(nsh + k * nat, nat) = m_wfn.dp_at.row(k).transpose();
            for (int k = 0; k < 6; ++k) v.segment(nsh + 3 * nat + k * nat, nat) = m_wfn.qp_at.row(k).transpose();
            return v;
        }
        return m_wfn.q_sh;
    };
    auto unpackSCC = [&](const Vector& v) {
        m_wfn.q_sh = v.head(nsh);
        rebuildQatFromQsh();
        if (m_method == MethodType::GFN2) {
            m_wfn.dp_at = Eigen::MatrixXd::Zero(3, nat);
            m_wfn.qp_at = Eigen::MatrixXd::Zero(6, nat);
            for (int k = 0; k < 3; ++k) m_wfn.dp_at.row(k) = v.segment(nsh + k * nat, nat).transpose();
            for (int k = 0; k < 6; ++k) m_wfn.qp_at.row(k) = v.segment(nsh + 3 * nat + k * nat, nat).transpose();
        }
    };

    const double damp   = (m_scf_damping > 0.0 && m_scf_damping <= 1.0) ? m_scf_damping : 0.4;
    const double thresh = m_scf_threshold;
    const int    maxit  = m_scf_max_iter;

    if (verb >= 1) {
        CurcumaLogger::result(fmt::format(
            "large_system_mode=sparse: purification SCF (nao={}, Nocc={}, P-threshold={:.1e})", nao, Nocc, threshold));
        CurcumaLogger::result("  iter        max|dq|     nnz(M)     path     t/ms");
    }

    BroydenMixer mixer(damp, 20);
    Vector x = packSCC();
    bool   use_purify = true;     // drops to false (dense) once a gapless step is hit
    bool   warned_gapless = false;
    double final_dq = 1.0e30;
    double last_nnz = 1.0;
    int iter = 0;
    for (iter = 0; iter < maxit; ++iter) {
        const auto t_it0 = clock::now();
        unpackSCC(x);

        // Global potential + Fock.
        m_pot.reset();
        addCoulombShellPotential(m_pot);
        addThirdOrderPotential(m_pot);
        if (m_method == MethodType::GFN2) {
            addMultipolePotential(m_pot);
            addDispersionPotential(m_pot);
        }
        const Matrix F = buildFock(m_H0, m_S, m_pot);

        const char* path = "purify";
        bool got_density = false;
        if (use_purify) {
            // M0 = lambda*(mu I - Sinv F) + (Nocc/nao) I, eigenvalues in [0,1], Tr=Nocc.
            const Eigen::MatrixXd Heff = Sinv * F;     // eigenvalues = orbital energies
            const double mubar = Heff.trace() / nao;
            // Tight orbital-energy bounds from the SYMMETRIC reduced matrix
            // Atil = L^-1 F L^-T (same eigenvalues as Sinv F = the orbital
            // energies). Power-iteration bounds — plain Gershgorin grows far too
            // loose with system size, collapsing lambda so M0 degenerates and
            // purification is misread as gapless (e.g. the 1410-atom polymer).
            const Eigen::MatrixXd Yr   = m_X.triangularView<Eigen::Lower>().solve(F);
            const Eigen::MatrixXd Atil = m_X.triangularView<Eigen::Lower>().solve(Yr.transpose());
            double emin, emax;
            spectralBounds(Atil, emin, emax);
            // Palser-Manolopoulos initial slope (note the 1/nao): the line
            //   d(e) = Nocc/nao + lambda*(mubar - e)
            // pivots at the mean energy with value Nocc/nao, so Tr(M0)=Nocc, and
            // lambda is the largest slope keeping d(e) in [0,1] across [emin,emax].
            const double dpos = emax - mubar, dneg = mubar - emin;
            double lambda = std::min(dpos > 1e-12 ? Nocc / (nao * dpos) : 1e30,
                                     dneg > 1e-12 ? (nao - Nocc) / (nao * dneg) : 1e30);
            if (!std::isfinite(lambda) || lambda <= 0.0) lambda = 1.0 / std::max(1.0, emax - emin);

            Eigen::MatrixXd M = lambda * (mubar * Eigen::MatrixXd::Identity(nao, nao) - Heff)
                              + (static_cast<double>(Nocc) / nao) * Eigen::MatrixXd::Identity(nao, nao);

            if (purifyProjector(M, Nocc)) {
                // One-shot prune of the converged projector to expose sparsity;
                // the pruned M feeds the density (so the SCF charges feel it).
                last_nnz = (threshold > 0.0) ? pruneFraction(M, threshold)
                                             : nnzFraction(M, 1e-12);
                // AO density D = 2 * P = 2 * M * S^-1 (symmetrise to clean rounding).
                Eigen::MatrixXd D = 2.0 * (M * Sinv);
                D = 0.5 * (D + D.transpose());
                m_wfn.P = D;
                got_density = true;
            } else {
                if (!warned_gapless && verb >= 1)
                    CurcumaLogger::warn("large_system_mode=sparse: density purification did not converge (gapless / no "
                                        "HOMO-LUMO gap?). Falling back to the eigensolver.");
                warned_gapless = true;
                use_purify = false;
            }
        }
        if (!got_density) {
            // Dense fallback (gapless): the validated generalized eigensolve.
            path = "eigen";
            if (!solveEigen(F, m_S)) {
                CurcumaLogger::warn("large_system_mode=sparse: eigensolver fallback failed");
                return 0.0;
            }
            last_nnz = 1.0;
        }

        updatePopulations(m_S);
        const Vector g = packSCC();
        final_dq = (g.head(nsh) - x.head(nsh)).cwiseAbs().maxCoeff();

        if (verb >= 1)
            CurcumaLogger::result(fmt::format("  {:>4}   {:>12.3e}   {:>8.4f}   {:>6}  {:>7.1f}",
                iter, final_dq, last_nnz, path,
                std::chrono::duration<double, std::milli>(clock::now() - t_it0).count()));

        if (final_dq < thresh) { converged_out = true; ++iter; break; }
        x = mixer.update(x, g);
    }
    iters_out = iter;
    nnz_frac_out = last_nnz;

    if (!converged_out)
        CurcumaLogger::warn(fmt::format("large_system_mode=sparse: SCF did not converge in {} iterations (max|dq|={:.2e})",
                                        maxit, final_dq));

    evaluateComponentsAtFixedDensity(m_wfn.P, m_wfn.q_at, m_wfn.q_sh, m_wfn.dp_at, m_wfn.qp_at);
    m_scf_converged  = converged_out;
    m_scf_iterations = iters_out;

    if (verb >= 1)
        CurcumaLogger::result(fmt::format(
            "large_system_mode=sparse: total energy = {:.8f} Eh  ({} it, nnz(M)={:.3f}, {:.0f} ms{})",
            m_E_total, iters_out, nnz_frac_out,
            std::chrono::duration<double, std::milli>(clock::now() - t0).count(),
            use_purify ? "" : ", DENSE fallback (gapless)"));
    return m_E_total;
}

} // namespace curcuma::xtb
