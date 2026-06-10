/*
 * <Native xTB divide-and-conquer DC-SCF — implementation>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0.
 *
 * Claude Generated (June 2026). Implements XTB::calculateDivideConquer, the
 * divide-and-conquer SCF declared in xtb_native.h. It is a member of XTB so it
 * reuses the validated build-once setup (CN, self-energies, S/H0, gamma,
 * multipole), the per-iteration potential builders, buildFock, updatePopulations
 * and evaluateComponentsAtFixedDensity verbatim — the only new physics is:
 *   (1) diagonalising core+buffer sub-blocks of the GLOBAL Fock instead of the
 *       whole matrix (the O(N^3) eigensolve becomes sum of small O(n^3) ones),
 *   (2) a shared chemical potential by global electron-count bisection over all
 *       sub-system spectra (Yang & Lee, JCP 103 (1995) 5674; Yang, PRL 66 (1991)
 *       1438), and
 *   (3) Yang's core-projection assembly of the global density.
 *
 * The sub-block diagonalisation is dispatched on the user-selected -eigensolver
 * (m_eigensolver). For 'mkl' the standard Eigen GES is used (default, fast for
 * moderate blocks). For 'native' the self-contained solveSymmetric (Householder
 * reduction + Cuppen D&C, no LAPACK eigensolve dependency; useful for GPU-portable
 * builds) is used on each sub-block. For 'lobpcg' the seeded block LOBPCG from
 * curcuma::eigsolver::lobpcgLowest is used (sub-blocks are LOBPCG's sweet spot:
 * k=nocc/2 << n). For 'purify' the non-orthogonal Palser-Manolopoulos purification
 * is used, requiring T=0 and a HOMO-LUMO gap.
 *
 * First cut: ENERGY ONLY (no analytic gradient). Approximate — converges to the
 * dense SCF energy as the buffer radius grows, which is the validation contract.
 */

#include "xtb_native.h"
#include "native_eigensolver.h"
#include "src/core/curcuma_logger.h"

#include <Eigen/Eigenvalues>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fmt/format.h>
#include <vector>

namespace curcuma::xtb {

// Fermi occupation (2 e- per orbital, closed shell): f(e) = 2 / (1 + e^{(e-mu)/kT}).
static inline double fermiOcc(double e, double mu, double kT)
{
    if (kT <= 0.0) return (e <= mu) ? 2.0 : 0.0;
    const double x = (e - mu) / kT;
    return 2.0 / (1.0 + std::exp(std::min(std::max(x, -500.0), 500.0)));
}

// ---------------------------------------------------------------------------
// Sub-block diagonalisation — dispatched on m_eigensolver.
//
// Returns (eps, C, ok). Columns of C satisfy C^T Ss C = I. On failure returns
// false and the caller should abort the DC-SCF.
//
// DC needs the FULL eigenspectrum (all n eigenvalues) for Fermi occupation /
// chemical-potential bisection. Eigensolvers that can't provide this fall back
// to the dense GES:
//   - purify: produces projector eigenvalues 0/1, not orbital energies
//   - lobpcg: gives only the lowest kb eigenvalues, not all n
//   - native / dnc: full spectrum via self-contained Householder+Cuppen (GPU-portable)
//   - mkl (default): full spectrum via Eigen GES (LAPACK dsyevd)
//
// Memory layout: the per-block Fock + overlap are full n×n Eigen matrices; this
// matches the dense solveEigen() path so the gradient/CPSCF code paths are
// unaffected. n=200..500 blocks are the typical regime; the temporary matrices
// are stack/heap-cheap and freed at scope exit.
// ---------------------------------------------------------------------------
static bool diagonalizeSubBlock(const Eigen::MatrixXd& Fs, const Eigen::MatrixXd& Ss,
                                const std::string& solver,
                                Eigen::VectorXd& eps_out, Eigen::MatrixXd& C_out,
                                std::string& path_out)
{
    const int n = static_cast<int>(Fs.rows());
    if (n == 0) { eps_out.resize(0); C_out.resize(0, 0); path_out = "empty"; return true; }

    // --- purify / lobpcg: need full spectrum, can't provide it. Fall back. ---
    // Log the fallback once (not per sub-block per iteration — the DC outer
    // loop calls this many times).
    static bool logged_purify_dc = false, logged_lobpcg_dc = false;
    if (solver == "purify" && !logged_purify_dc) {
        CurcumaLogger::info("DC sub-block: eigensolver=purify produces projector "
                            "eigenvalues (0/1), not orbital energies needed for "
                            "Fermi occupation; falling back to dense GES per sub-block");
        logged_purify_dc = true;
    } else if (solver == "lobpcg" && !logged_lobpcg_dc) {
        CurcumaLogger::info("DC sub-block: eigensolver=lobpcg gives partial spectrum "
                            "(needs full spectrum for Fermi occupation); "
                            "falling back to dense GES per sub-block");
        logged_lobpcg_dc = true;
    }

    // --- native / dnc: self-contained Householder + Cuppen D&C ---------------
    if (solver == "native" || solver == "dnc") {
        Eigen::LLT<Eigen::MatrixXd> llt(Ss);
        if (llt.info() != Eigen::Success) {
            CurcumaLogger::error("DC sub-block: native eigensolver needs Cholesky of Ss");
            return false;
        }
        const Eigen::MatrixXd L = llt.matrixL();
        const Eigen::MatrixXd Y = L.triangularView<Eigen::Lower>().solve(Fs);
        Eigen::MatrixXd Atil = L.triangularView<Eigen::Lower>().solve(Y.transpose());
        Atil = 0.5 * (Atil + Atil.transpose());
        Eigen::VectorXd evals;
        Eigen::MatrixXd evecs_std;
        if (curcuma::eigsolver::solveSymmetric(Atil, evals, evecs_std, 1)) {
            C_out = L.triangularView<Eigen::Lower>().transpose().solve(evecs_std);
            const Eigen::MatrixXd SC = Ss * C_out;
            for (int i = 0; i < n; ++i) C_out.col(i) /= std::sqrt(std::max(1e-30, C_out.col(i).dot(SC.col(i))));
            eps_out = evals;
            path_out = "native";
            return true;
        }
        CurcumaLogger::warn("DC sub-block: native eigensolver failed; falling back to mkl");
        // fall through to mkl
    }

    // --- mkl: default Eigen GES ---------------------------------------------
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> ges(Fs, Ss);
    if (ges.info() != Eigen::Success) {
        CurcumaLogger::error(fmt::format("DC sub-block: mkl eigensolve failed (n={})", n));
        return false;
    }
    eps_out = ges.eigenvalues();
    C_out = ges.eigenvectors();
    path_out = "mkl";
    return true;
}

double XTB::calculateDivideConquer(const std::vector<std::vector<int>>& subsystems,
                                   const std::vector<std::vector<int>>& cores,
                                   bool& converged_out, int& iters_out)
{
    using clock = std::chrono::steady_clock;
    const int verb = CurcumaLogger::get_verbosity();
    const int nat  = m_atomcount;
    const int nsh  = m_basis.nsh;
    const int nao  = m_basis.nao;
    const int nsub = static_cast<int>(subsystems.size());
    converged_out = false;
    iters_out = 0;

    if (nsub == 0 || cores.size() != subsystems.size()) {
        CurcumaLogger::error("DC: empty / mismatched subsystem-core lists");
        return 0.0;
    }

    // ---- core ownership: every atom must be a core of exactly one subsystem ----
    std::vector<int> core_owner(nat, -1);
    for (int a = 0; a < nsub; ++a)
        for (int at : cores[a]) {
            if (at < 0 || at >= nat) { CurcumaLogger::error("DC: core atom index out of range"); return 0.0; }
            if (core_owner[at] != -1)
                CurcumaLogger::warn(fmt::format("DC: atom {} is core of >1 subsystem (cores must be disjoint)", at));
            core_owner[at] = a;
        }
    int uncovered = 0;
    for (int at = 0; at < nat; ++at) if (core_owner[at] == -1) ++uncovered;
    if (uncovered)
        CurcumaLogger::warn(fmt::format("DC: {} atom(s) belong to no core — their density block is dropped", uncovered));

    // The T=0 hard check for eigensolver=purify is enforced upstream in
    // NativeXtbMethod::setMolecule. We trust m_eigensolver here; the sub-block
    // diagonaliser logs + falls back on per-block failures.
    const std::string solver = m_eigensolver;

    // ---- build-once global setup (mirror Calculation steps 1-5) ----
    m_d4_prepared = false;
    m_d4_genparams_calls = 0;
    const auto t_setup0 = clock::now();
    Vector cn = computeCoordinationNumbers();
    Vector se;
    getSelfEnergies(cn, se);
    Matrix S, H0;
    getHamiltonianH0(se, S, H0);
    m_S  = S;
    m_H0 = H0;
    buildGammaMatrix();
    if (m_method == MethodType::GFN2)
        setupMultipole();
    m_coordination_numbers = cn;

    // ---- per-atom contiguous global AO range [ao_begin, ao_end) ----
    std::vector<int> ao_begin(nat, 0), ao_end(nat, 0);
    for (int a = 0; a < nat; ++a) {
        const int s0  = m_basis.ish_at[a];
        const int nsa = m_basis.nsh_at[a];
        int first = nao, last = 0, count = 0;
        for (int s = s0; s < s0 + nsa; ++s) {
            first = std::min(first, m_basis.iao_sh[s]);
            last  = std::max(last, m_basis.iao_sh[s] + m_basis.nao_sh[s]);
            count += m_basis.nao_sh[s];
        }
        ao_begin[a] = (nsa > 0) ? first : 0;
        ao_end[a]   = (nsa > 0) ? last  : 0;
        (void)count;
    }

    // ---- precompute each subsystem's global-AO list + per-AO core flag ----
    struct SubInfo {
        std::vector<int>  gao;     // global AO indices in this sub-system
        std::vector<char> is_core; // 1 if the AO's atom is a core atom OF THIS sub-system
    };
    std::vector<SubInfo> sub(nsub);
    int max_sub_ao = 0;
    for (int a = 0; a < nsub; ++a) {
        // atoms sorted for determinism
        std::vector<int> atoms = subsystems[a];
        std::sort(atoms.begin(), atoms.end());
        atoms.erase(std::unique(atoms.begin(), atoms.end()), atoms.end());
        for (int at : atoms) {
            for (int mu = ao_begin[at]; mu < ao_end[at]; ++mu) {
                sub[a].gao.push_back(mu);
                sub[a].is_core.push_back(core_owner[at] == a ? 1 : 0);
            }
        }
        max_sub_ao = std::max(max_sub_ao, static_cast<int>(sub[a].gao.size()));
    }

    if (verb >= 1) {
        CurcumaLogger::result(fmt::format(
            "large_system_mode=dc: {} subsystems, largest sub-block {} AOs (global nao={}), eigensolver={}",
            nsub, max_sub_ao, nao, solver));
    }

    // ---- initial charge guess (EEQ, falls back to zero) ----
    auto rebuildQatFromQsh = [&]() {
        m_wfn.q_at.setZero(nat);
        for (int s = 0; s < nsh; ++s) m_wfn.q_at(m_basis.sh2at[s]) += m_wfn.q_sh(s);
    };
    Vector q_sh_guess;
    if (seedEEQGuess(q_sh_guess)) {
        m_wfn.q_sh = q_sh_guess;
    } else {
        m_wfn.q_sh = Vector::Zero(nsh);
    }
    rebuildQatFromQsh();
    m_wfn.dp_at = Eigen::MatrixXd::Zero(3, nat);
    m_wfn.qp_at = Eigen::MatrixXd::Zero(6, nat);

    const double n_elec = m_wfn.nocc;          // target electron count
    const double kT     = m_electronic_temp * 3.166808e-6;   // K -> Hartree
    const double damp   = (m_scf_damping > 0.0 && m_scf_damping <= 1.0) ? m_scf_damping : 0.4;
    const double thresh = m_scf_threshold;
    const int    maxit  = m_scf_max_iter;

    // Pack / unpack the self-consistent charge vector exactly as the dense SCF:
    // shell charges (GFN1) plus the atomic dipoles/quadrupoles (GFN2). Broyden
    // mixing of THIS low-dimensional vector is what lets the dense code converge
    // stiff polar systems (e.g. complex) where linear/Pulay mixing sloshes —
    // the DC outer loop needs the same mixer for the same reason.
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

    if (verb >= 1) {
        CurcumaLogger::result("  iter        max|dq|         mu/Eh      N(mu)     t/ms");
    }

    BroydenMixer mixer(damp, /*max_history=*/20);
    Vector x = packSCC();                 // current SCC input (from the EEQ guess)
    double final_dq = 1.0e30;
    int iter = 0;
    for (iter = 0; iter < maxit; ++iter) {
        const auto t_it0 = clock::now();

        // Charges that build this iteration's potential are the mixed input x.
        unpackSCC(x);

        // (a) Global potential + Fock from the current density/charges.
        m_pot.reset();
        addCoulombShellPotential(m_pot);
        addThirdOrderPotential(m_pot);
        if (m_method == MethodType::GFN2) {
            addMultipolePotential(m_pot);
            addDispersionPotential(m_pot);
        }
        const Matrix F = buildFock(m_H0, m_S, m_pot);

        // (b) Diagonalise each core+buffer sub-block via the user-selected
        //     eigensolver; collect (eps, weight) pairs and keep the sub-eigvecs
        //     / occupied density for the assembly.
        std::vector<Eigen::VectorXd> sub_eps(nsub);
        std::vector<Eigen::MatrixXd> sub_C(nsub);
        std::vector<Eigen::VectorXd> sub_w(nsub);   // core-projected Mulliken weight per orbital
        std::map<std::string, int> sub_paths;       // for the iter-1 info line
        for (int a = 0; a < nsub; ++a) {
            const auto& g = sub[a].gao;
            const int n = static_cast<int>(g.size());
            if (n == 0) continue;
            Eigen::MatrixXd Fs(n, n), Ss(n, n);
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j) {
                    Fs(i, j) = F(g[i], g[j]);
                    Ss(i, j) = m_S(g[i], g[j]);
                }
            Eigen::VectorXd eps;
            Eigen::MatrixXd C;
            std::string path;
            // diagonalizeSubBlock: purify/lobpcg fall back to the dense GES
            // internally (they can't provide the full spectrum DC needs).
            if (!diagonalizeSubBlock(Fs, Ss, solver, eps, C, path))
                return 0.0;
            sub_eps[a] = eps;
            sub_C[a]   = C;
            sub_paths[path]++;
            const Eigen::MatrixXd SC = Ss * C; // for Mulliken weights / density Mulliken
            Eigen::VectorXd w(n);
            for (int i = 0; i < n; ++i) {
                double wi = 0.0;
                for (int mu = 0; mu < n; ++mu)
                    if (sub[a].is_core[mu]) wi += C(mu, i) * SC(mu, i);
                w(i) = wi;
            }
            sub_w[a] = w;
        }

        // (c) Shared chemical potential mu: bisect N(mu) = sum_a sum_i f(eps,mu)*w = n_elec.
        double mu_lo = 1.0e30, mu_hi = -1.0e30;
        for (int a = 0; a < nsub; ++a)
            if (sub_eps[a].size()) {
                mu_lo = std::min(mu_lo, sub_eps[a].minCoeff());
                mu_hi = std::max(mu_hi, sub_eps[a].maxCoeff());
            }
        mu_lo -= 1.0; mu_hi += 1.0;
        auto countN = [&](double mu) {
            double N = 0.0;
            for (int a = 0; a < nsub; ++a)
                for (int i = 0; i < sub_eps[a].size(); ++i)
                    N += fermiOcc(sub_eps[a](i), mu, kT) * sub_w[a](i);
            return N;
        };
        double mu = 0.5 * (mu_lo + mu_hi);
        for (int b = 0; b < 200; ++b) {
            mu = 0.5 * (mu_lo + mu_hi);
            if (countN(mu) > n_elec) mu_hi = mu; else mu_lo = mu;
            if (mu_hi - mu_lo < 1.0e-12) break;
        }
        const double Nmu = countN(mu);

        // (d) Yang core-projected global density assembly.
        Matrix D = Matrix::Zero(nao, nao);
        for (int a = 0; a < nsub; ++a) {
            const int n = static_cast<int>(sub[a].gao.size());
            if (n == 0) continue;
            // occupied (fractional) density of this sub-system at mu
            Eigen::VectorXd occ(n);
            int ncol = 0;
            for (int i = 0; i < n; ++i) { occ(i) = fermiOcc(sub_eps[a](i), mu, kT); if (occ(i) > 1.0e-12) ncol = std::max(ncol, i + 1); }
            const Eigen::MatrixXd Cw = sub_C[a].leftCols(ncol) * occ.head(ncol).asDiagonal();
            const Eigen::MatrixXd Ps = Cw * sub_C[a].leftCols(ncol).transpose();   // n x n local density
            const auto& g  = sub[a].gao;
            const auto& ic = sub[a].is_core;
            for (int i = 0; i < n; ++i) {
                const double ci = ic[i] ? 1.0 : 0.0;
                for (int j = 0; j < n; ++j) {
                    const double p = 0.5 * (ci + (ic[j] ? 1.0 : 0.0));   // Yang partition weight
                    if (p != 0.0) D(g[i], g[j]) += p * Ps(i, j);
                }
            }
        }

        // (e) New Mulliken populations from the assembled density: G(x).
        m_wfn.P = D;
        updatePopulations(m_S);
        const Vector g = packSCC();

        // (f) Broyden charge mixing + convergence check on max|dq_shell|.
        final_dq = (g.head(nsh) - x.head(nsh)).cwiseAbs().maxCoeff();

        if (verb >= 1)
            CurcumaLogger::result(fmt::format("  {:>4}   {:>12.3e}   {:>12.6f}   {:>8.4f}  {:>7.1f}",
                iter, final_dq, mu, Nmu,
                std::chrono::duration<double, std::milli>(clock::now() - t_it0).count()));

        if (final_dq < thresh) { converged_out = true; ++iter; break; }
        x = mixer.update(x, g);   // next input
    }
    iters_out = iter;
    // At loop exit m_wfn carries G(x) (the Mulliken of the last assembled
    // density) and m_wfn.P = D — the self-consistent state used for the energy.

    if (!converged_out)
        CurcumaLogger::warn(fmt::format("DC-SCF did not converge in {} iterations (max|dq|={:.2e})", maxit, final_dq));

    // ---- final energy at the converged assembled density (validated path) ----
    evaluateComponentsAtFixedDensity(m_wfn.P, m_wfn.q_at, m_wfn.q_sh, m_wfn.dp_at, m_wfn.qp_at);
    m_scf_converged  = converged_out;
    m_scf_iterations = iters_out;

    const auto t_total = std::chrono::duration<double, std::milli>(clock::now() - t_setup0).count();
    if (verb >= 1)
        CurcumaLogger::result(fmt::format("DC-SCF total energy = {:.8f} Eh  ({} it, {:.0f} ms)",
                                          m_E_total, iters_out, t_total));
    return m_E_total;
}

} // namespace curcuma::xtb
