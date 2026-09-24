/*
 * <Native KS-DFT -- WP3 closed-shell HF-SCF>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Closed-shell Roothaan-Hall SCF for the native DFT engine. This is the
 * "hard ERI gate" of the roadmap: it exercises the WP1 1e matrices and the
 * WP2 4-centre ERI end-to-end and, for the first time, makes the MO spectrum
 * + total energy vs ORCA a valid comparison (ORCA orbital energies are the SCF
 * Fock spectrum H + 2J - K, not the 1e Hcore that WP1 had alone).
 *
 * Convention (Szabo-Ostlund / PySCF, spin-summed density):
 *   P_munu  = 2 * sum_{i=1}^{n_occ} C_mui C_nui           (factor 2 = spin)
 *   J_munu = sum_{lam,sig} P_lamsig (mu nu | lam sig)     (Coulomb)
 *   K_munu = sum_{lam,sig} P_lamsig (mu lam | nu sig)     (exchange)
 *   F      = H + J - 1/2 K                                 (RHF Fock)
 *   E_elec = 1/2 * Tr(P (H + F))                           (total electronic)
 *   components: ET = Tr(P T), EV = Tr(P V), EJ = 1/2 Tr(P J),
 *               Ex = -1/4 Tr(P K);  E_elec = (ET + EV) + EJ + Ex.
 * (The roadmap note used the equivalent spatial-density form F = H + 2J - K;
 * both give the same Fock eigenvalues and the same total energy. The form here
 * is the standard spin-summed one and matches buildCoulomb/buildExchange with a
 * factor-2 density directly.)
 *
 * SCF driver mirrors the NDDO Fock-DIIS pattern (nddo.cpp runSCF): build Fock,
 * DIIS-extrapolate once the history is long enough, diagonalize, rebuild the
 * density, damp, check ||dP||. The generalized eigenproblem F C = S C eps is
 * reduced to the standard one via the one-time Lowdin orthonormalizer
 * X = S^{-1/2} (no generalized eigensolver exists in curcuma::eigsolver).
 *
 * Literature:
 *   - P. Pulay, Chem. Phys. Lett. 73, 393 (1980); J. Comput. Chem. 3, 556 (1982).
 *   - T. Helgaker, P. Jorgensen, J. Olsen, Molecular Electronic-Structure
 *     Theory (Wiley, 2000), ch. 10 (SCF) + ch. 9.9 (ERI).
 *   - A. Szabo, N. S. Ostlund, Modern Quantum Chemistry (Dover, 1996), ch. 3.
 *   - xcDFT RKS.f90 (TCCM winter school 2019: DFT) -- structural SCF template.
 *
 * Claude Generated: WP3 HF-SCF
 *
 * This program is free software under GPL-3.0
 */

#include "qm_engine.h"
#include "src/core/curcuma_logger.h"

#include "native_eigensolver.h"

#include <Eigen/Dense>
#include <fmt/format.h>
#include <chrono>
#include <cmath>
#include <string>

namespace {
// Symmetric Loewdin orthonormalizer X = S^{-1/2} via S = U diag(w) U^T,
// X = U diag(1/sqrt(w)) U^T. Returns true on success; leaves X empty on a
// singular / indefinite S (caller reports the error).
// Claude Generated: Lowdin S^{-1/2} for the generalized SCF eigenproblem.
bool lowdinOrthonormalizer(const Eigen::MatrixXd& S, Eigen::MatrixXd& X)
{
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(S);
    if (es.info() != Eigen::Success) return false;
    const Eigen::VectorXd w = es.eigenvalues();
    const Eigen::MatrixXd U = es.eigenvectors();
    // Reject a singular / indefinite overlap (linearly dependent basis).
    for (int i = 0; i < w.size(); ++i)
        if (w(i) < 1.0e-10) return false;
    Eigen::VectorXd inv_sqrt(w.size());
    for (int i = 0; i < w.size(); ++i) inv_sqrt(i) = 1.0 / std::sqrt(w(i));
    X = U * inv_sqrt.asDiagonal() * U.transpose();
    return true;
}
}  // namespace


// =================================================================================
// SCF convergence accelerator for the QM engine: Pulay DIIS and ADIIS
// (Claude Generated, Sep 2026). Self-contained here on purpose -- the shared
// DIISAccelerator (diis_accelerator.h) also drives the native xTB SCF, whose
// energies are gated at 1e-8 against tblite, and stays untouched.
//
//  DIIS  (P. Pulay, Chem. Phys. Lett. 73, 393 (1980); J. Comput. Chem. 3, 556 (1982)):
//        F = sum c_i F_i with sum c_i = 1 minimising |sum c_i e_i|^2, where
//        e_i = X^T (F_i P_i S - S P_i F_i) X is the commutator in the ORTHONORMAL
//        basis (X = S^-1/2) -- the AO-basis commutator is not a proper residual
//        when S is far from the identity, which is the Pulay-1982 recommendation.
//  ADIIS (X. Hu, W. Yang, J. Chem. Phys. 132, 054109 (2010)): far from convergence
//        the commutator is a poor guide; instead minimise the quadratic (ARH) model
//        of the energy over convex combinations of the stored densities,
//          f(c) = sum c_i Tr[(P_i - P_n) F_n] + 1/2 sum c_i c_j Tr[(P_i - P_n)(F_j - F_n)],
//        c_i >= 0, sum c_i = 1  (n = newest entry; spin-summed P, so dE = Tr(F dP)).
//  Blend (R. Garza, G. Scuseria, J. Chem. Phys. 137, 054110 (2012)): ADIIS while
//        max|e| > 1e-1, DIIS below 1e-4, linear mix of the two Fock matrices between.
// =================================================================================
namespace {

class QMScfAccelerator {
public:
    QMScfAccelerator(int max_size, bool adiis) : m_max(std::max(2, max_size)), m_adiis(adiis) {}

    /// Store F (built from P) and return max|e| of its orthonormal-basis error.
    double push(const Matrix& F, const Matrix& P, const Matrix& S, const Matrix& X)
    {
        const Matrix FPS = F * P * S;
        const Matrix e = X.transpose() * (FPS - FPS.transpose()) * X;
        m_F.push_back(F);
        m_P.push_back(P);
        m_e.push_back(e);
        if ((int)m_F.size() > m_max) {
            m_F.erase(m_F.begin());
            m_P.erase(m_P.begin());
            m_e.erase(m_e.begin());
        }
        m_last_err = e.cwiseAbs().maxCoeff();
        return m_last_err;
    }

    /// Fock matrix to diagonalise next.
    Matrix extrapolate()
    {
        if (m_F.size() < 2) return m_F.back();
        if (!m_adiis || m_last_err < 1e-4) return diis();
        if (m_last_err > 1e-1) return adiis();
        const double w = (std::log10(m_last_err) + 4.0) / 3.0;  // 0 at 1e-4 .. 1 at 1e-1
        return w * adiis() + (1.0 - w) * diis();
    }

private:
    int m_max;
    bool m_adiis;
    double m_last_err = 0.0;
    std::vector<Matrix> m_F, m_P, m_e;

    Matrix combine(const Eigen::VectorXd& c) const
    {
        Matrix F = Matrix::Zero(m_F[0].rows(), m_F[0].cols());
        for (int i = 0; i < c.size(); ++i) F += c(i) * m_F[m_F.size() - c.size() + i];
        return F;
    }

    Matrix diis()
    {
        // Solve the bordered Pulay system; if it is singular (nearly dependent error
        // vectors), drop the oldest entries until it is not.
        for (int first = 0; first + 2 <= (int)m_F.size(); ++first) {
            const int n = (int)m_F.size() - first;
            Eigen::MatrixXd B = Eigen::MatrixXd::Zero(n + 1, n + 1);
            for (int i = 0; i < n; ++i)
                for (int j = 0; j <= i; ++j)
                    B(i, j) = B(j, i) = m_e[first + i].cwiseProduct(m_e[first + j]).sum();
            const double scale = B.topLeftCorner(n, n).diagonal().maxCoeff();
            if (scale <= 0.0) return m_F.back();
            B.topLeftCorner(n, n) /= scale;  // conditioning only; the solution is scale-free
            for (int i = 0; i < n; ++i) B(i, n) = B(n, i) = -1.0;
            Eigen::VectorXd rhs = Eigen::VectorXd::Zero(n + 1);
            rhs(n) = -1.0;
            Eigen::FullPivLU<Eigen::MatrixXd> lu(B);
            if (lu.rank() < n + 1 || lu.rcond() < 1e-14) continue;
            const Eigen::VectorXd c = lu.solve(rhs).head(n);
            if (!c.allFinite()) continue;
            return combine(c);
        }
        return m_F.back();
    }

    Matrix adiis() const
    {
        const int n = (int)m_F.size();
        const Matrix& Pn = m_P.back();
        const Matrix& Fn = m_F.back();
        Eigen::VectorXd a(n);
        Eigen::MatrixXd M(n, n);
        for (int i = 0; i < n; ++i) {
            const Matrix dP = m_P[i] - Pn;
            a(i) = dP.cwiseProduct(Fn).sum();
            for (int j = 0; j < n; ++j)
                M(i, j) = dP.cwiseProduct(m_F[j] - Fn).sum();
        }
        M = 0.5 * (M + M.transpose()).eval();
        // c_i = t_i^2 / sum t^2 keeps c on the simplex; gradient descent with
        // backtracking on t (n <= subspace size, so this is cheap).
        auto fval = [&](const Eigen::VectorXd& c) { return a.dot(c) + 0.5 * c.dot(M * c); };
        auto toC = [](const Eigen::VectorXd& t) { return Eigen::VectorXd(t.cwiseAbs2() / t.squaredNorm()); };
        Eigen::VectorXd t = Eigen::VectorXd::Constant(n, 1.0);
        t(n - 1) = 2.0;  // lean towards the newest density
        Eigen::VectorXd c = toC(t);
        double f = fval(c), step = 1.0;
        for (int it = 0; it < 500; ++it) {
            const Eigen::VectorXd g = a + M * c;
            const double cg = c.dot(g);
            Eigen::VectorXd gt(n);
            const double S = t.squaredNorm();
            for (int k = 0; k < n; ++k) gt(k) = 2.0 * t(k) / S * (g(k) - cg);
            if (gt.norm() < 1e-12) break;
            bool accepted = false;
            for (int ls = 0; ls < 40; ++ls) {
                const Eigen::VectorXd tn = t - step * gt;
                const Eigen::VectorXd cn = toC(tn);
                const double fn = fval(cn);
                if (fn < f - 1e-4 * step * gt.squaredNorm()) {
                    t = tn; c = cn; f = fn; step *= 2.0; accepted = true;
                    break;
                }
                step *= 0.5;
            }
            if (!accepted) break;
        }
        return combine(c);
    }
};

}  // namespace

// =================================================================================
// SCF driver
// =================================================================================

bool QMEngine::runSCF()
{
    // Claude Generated: WP3 closed-shell RHF SCF (DIIS / plain damping).
    const int n = m_nbf;
    if (n <= 0) {
        CurcumaLogger::error("QM SCF: no basis functions (InitialiseMolecule failed?)");
        return false;
    }
    // Closed-shell only: even electron count.
    if (m_num_electrons < 0 || (m_num_electrons % 2) != 0) {
        CurcumaLogger::error(fmt::format(
            "QM SCF: closed-shell HF needs an even electron count, got {}",
            m_num_electrons));
        return false;
    }
    const int n_occ = m_num_electrons / 2;
    if (n_occ > n) {
        CurcumaLogger::error(fmt::format(
            "QM SCF: {} electrons need {} occupied orbitals but basis has only {}",
            m_num_electrons, n_occ, n));
        return false;
    }

    // One-time setup: Lowdin orthonormalizer + active-basis ERI.
    buildOrthonormalizer();
    if (m_X.size() == 0) {
        CurcumaLogger::error("QM SCF: S^{-1/2} construction failed (linear basis?)");
        return false;
    }
    if (CurcumaLogger::get_verbosity() >= 3) {
        Eigen::MatrixXd Scol(m_S), Hcol(m_H);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esS(Scol);
        const Eigen::VectorXd sw = esS.eigenvalues();
        std::string svals = "S eigvals:";
        for (int i = 0; i < sw.size(); ++i) svals += fmt::format(" {:.4e}", sw(i));
        CurcumaLogger::info(svals);
        // Cross-check: direct generalized eigensolver for the bare H (iter 0 Fock).
        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> ges(Hcol, Scol);
        const Eigen::VectorXd hv = ges.eigenvalues();
        std::string hvals = "H gen-eigvals:";
        for (int i = 0; i < std::min<int>(5, hv.size()); ++i)
            hvals += fmt::format(" {:.4f}", hv(i));
        CurcumaLogger::info(hvals);
    }
    // Two-electron part: stored tensor (built once here) or integral-direct.
    const bool direct = useDirectSCF();
    if (!direct) {
        const qmint::ERITensor& eri = eriActive();  // builds (and caches) on first call
        if (eri.n() != n) {
            CurcumaLogger::error(fmt::format(
                "QM SCF: active ERI dimension {} != basis dimension {}", eri.n(), n));
            return false;
        }
    }
    m_direct_quartets = 0;

    const bool use_diis = (m_scf_mode != "plain");   // diis | adiis
    if (m_scf_mode != "plain" && m_scf_mode != "diis" && m_scf_mode != "adiis" && CurcumaLogger::get_verbosity() >= 1)
        CurcumaLogger::warn(fmt::format("QM: unknown -qm.scf_mode '{}' (use diis|adiis|plain); using diis", m_scf_mode));
    const double damping = use_diis ? 0.0 : 0.5;  // plain mode needs density damping

    // Initial guess (m_scf_guess): SAD (default) or the bare core Hamiltonian
    // (zero density -> first Fock = H).
    m_density = buildInitialGuess();
    if (CurcumaLogger::get_verbosity() >= 3)
        CurcumaLogger::param("SCF guess",
                             fmt::format("{} (Tr(PS) = {:.6f})", m_last_warm_start ? "warm start" : m_scf_guess,
                                         m_density.cwiseProduct(m_S).sum()));

    QMScfAccelerator diis(m_diis_subspace, m_scf_mode == "adiis");
    const auto t_scf0 = std::chrono::steady_clock::now();
    double t_fock_ms = 0.0;

    m_scf_converged = false;
    m_scf_iterations = 0;

    // Incremental Fock build (integral-direct only): J and K are linear in P, so
    // J(P_n) = J(P_ref) + J(P_n - P_ref). The density change shrinks as the SCF
    // converges and the density-weighted screening then skips most quartets.
    // A full build every kDirectRebuild iterations bounds the accumulated
    // screening error; the final energy always comes from a full build.
    constexpr int kDirectRebuild = 8;
    Matrix J_acc, K_acc, P_ref;
    int since_full = 0;
    auto fockFromDensity = [&](const Matrix& P) -> Matrix {
        if (!direct || P_ref.size() == 0 || since_full >= kDirectRebuild) {
            buildJK(P, J_acc, K_acc);
            since_full = 0;
        } else {
            Matrix dJ, dK;
            buildJK(P - P_ref, dJ, dK);
            J_acc += dJ;
            K_acc += dK;
            ++since_full;
        }
        P_ref = P;
        return m_H + J_acc - 0.5 * K_acc;
    };

    for (int iter = 0; iter < m_scf_max_iter; ++iter) {
        m_scf_iterations = iter + 1;

        // Fock from the current density (F = H + J - 1/2 K, spin-summed convention).
        const auto t_f0 = std::chrono::steady_clock::now();
        Matrix fock = fockFromDensity(m_density);
        t_fock_ms += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t_f0).count();

        if (CurcumaLogger::get_verbosity() >= 3 && iter < 2) {
            // Inspect J/K in the MO basis of the *previous* orbitals (m_prev_C).
            Matrix Jm, Km;
            buildJK(m_density, Jm, Km);
            CurcumaLogger::info(fmt::format(
                "DBG Fock diag[0..3]={:.4} {:.4} {:.4} {:.4}  J_00(AO)={:.4} K_00(AO)={:.4} Tr(PJ)={:.6} Tr(PK)={:.6}",
                fock(0, 0), n > 1 ? fock(1, 1) : 0, n > 2 ? fock(2, 2) : 0, n > 3 ? fock(3, 3) : 0,
                Jm(0, 0), Km(0, 0),
                m_density.cwiseProduct(Jm).sum(), m_density.cwiseProduct(Km).sum()));
        }

        // DIIS: push the commutator error e = F P S - S P F every iter once the
        // history can support it; use the extrapolated Fock for the diagonalization.
        Matrix fock_use = fock;
        if (use_diis && iter >= m_diis_start) {
            diis.push(fock, m_density, m_S, m_X);
            fock_use = diis.extrapolate();
        }

        // Solve F C = S C eps (Lowdin-reduced standard problem).
        Matrix C;
        Vector eps;
        if (!solveFock(fock_use, C, eps)) {
            CurcumaLogger::error(fmt::format(
                "QM SCF: diagonalization failed at iteration {}", iter + 1));
            return false;
        }

        // New density from the occupied MOs.
        Matrix density_new = buildDensity(C, n_occ);

        // Convergence: Frobenius norm of the density change.
        const double delta_P = (density_new - m_density).norm();

        if (CurcumaLogger::get_verbosity() >= 3 && iter < 3) {
            const double trPS = (m_density.cwiseProduct(m_S).sum());
            const double trPSnew = (density_new.cwiseProduct(m_S).sum());
            const double cSc = (C.col(0).transpose() * m_S * C.col(0)).value();
            CurcumaLogger::info(fmt::format(
                "DBG iter {} Tr(PS)={} Tr(PS_new)={} c0Sc={} eps0={:.6} eps1={:.6} ||P||={:.3}",
                iter + 1, trPS, trPSnew, cSc, eps.size() > 0 ? eps(0) : 0.0,
                eps.size() > 1 ? eps(1) : 0.0, density_new.norm()));
        }

        if (CurcumaLogger::get_verbosity() >= 3) {
            double e_tmp = 0.5 * (density_new.cwiseProduct(m_H + fock).sum());
            CurcumaLogger::param(fmt::format("SCF iter {}", iter + 1),
                                 fmt::format("dP = {:.3e}, E_elec = {:.10f} Eh", delta_P, e_tmp));
        }

        if (delta_P < m_scf_threshold) {
            m_density = density_new;
            // Self-consistent final Fock from ONE full J/K build (never incremental).
            Matrix J, K;
            buildJK(m_density, J, K);
            m_fock = m_H + J - 0.5 * K;
            m_scf_converged = true;

            // Final energy + components from the converged density and Fock.
            m_et = m_density.cwiseProduct(m_T).sum();
            m_ev = m_density.cwiseProduct(m_V).sum();
            m_ej = 0.5 * m_density.cwiseProduct(J).sum();
            m_ex = -0.25 * m_density.cwiseProduct(K).sum();
            m_e_elec = 0.5 * m_density.cwiseProduct(m_H + m_fock).sum();

            // Store the MO spectrum for wrapper compatibility (QMDriver slots).
            m_mo = C;
            // Keep the occupied orbitals for the next geometry's warm start.
            m_warm_C = C.leftCols(n_occ);
            m_warm_atoms = m_atoms;
            m_energies = eps;

            if (CurcumaLogger::get_verbosity() >= 2) {
                CurcumaLogger::success(fmt::format(
                    "QM HF SCF converged in {} iterations (dP = {:.3e}, start: {})",
                    iter + 1, delta_P, m_last_warm_start ? "previous geometry" : m_scf_guess));
                CurcumaLogger::info(fmt::format("QM: SCF loop {:.1f} ms, of which Fock (J/K) builds {:.1f} ms{}",
                    std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t_scf0).count(),
                    t_fock_ms,
                    direct ? fmt::format(" (integral-direct, {} shell quartets computed)", m_direct_quartets) : std::string()));
            }
            return true;
        }

        // Density update: damped (plain) or direct (diis -- the extrapolation
        // already smoothed the Fock). Reuses NDDO's update form.
        m_density = damping * m_density + (1.0 - damping) * density_new;
    }

    // Ran out of iterations: keep the last density/Fock for diagnostics.
    m_fock = buildFock(m_density);
    m_e_elec = 0.5 * m_density.cwiseProduct(m_H + m_fock).sum();
    return false;
}

// =================================================================================
// SCF initial guess
// =================================================================================

// Atom index of each active basis function. The cartesian -> spherical transform is
// block-diagonal per shell (Q's column j has nonzeros only in its own shell's
// cartesian rows), so every active AO belongs to exactly one atom.
std::vector<int> QMEngine::activeAtomIndex() const
{
    const int n = m_nbf;
    std::vector<int> ao_atom(n, -1);
    const int ncart = static_cast<int>(m_gto_basis.size());
    if (m_Q.size() == 0) {
        for (int i = 0; i < n && i < ncart; ++i) ao_atom[i] = m_gto_basis[i].atom;
        return ao_atom;
    }
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < ncart; ++i) {
            if (std::abs(m_Q(i, j)) > 1.0e-10) { ao_atom[j] = m_gto_basis[i].atom; break; }
        }
    }
    return ao_atom;
}

// Claude Generated: SAD (superposition of atomic densities) initial guess.
// For every atom, diagonalize that atom's block of the core Hamiltonian in its own
// atomic basis and fill the atom's electrons into the resulting atomic orbitals
// (aufbau; the last orbital is fractionally occupied when Z is odd). Summing the
// atomic densities gives Tr(P S) = N and a start far closer to the physical
// solution than the bare core guess -- which is what keeps the HF SCF from settling
// on a secondary solution (BH is the case in the validation set).
Matrix QMEngine::buildAtomicGuess() const
{
    const int n = m_nbf;
    const std::vector<int> ao_atom = activeAtomIndex();
    Matrix P = Matrix::Zero(n, n);
    double occupied = 0.0;

    const Eigen::MatrixXd Hcol(m_H);   // RowMajor -> ColMajor for Eigen
    const Eigen::MatrixXd Scol(m_S);

    for (int a = 0; a < m_atomcount; ++a) {
        std::vector<int> idx;
        for (int i = 0; i < n; ++i)
            if (ao_atom[i] == a) idx.push_back(i);
        const int m = static_cast<int>(idx.size());
        if (m == 0) continue;

        Eigen::MatrixXd Hs(m, m), Ss(m, m);
        for (int p = 0; p < m; ++p)
            for (int q = 0; q < m; ++q) {
                Hs(p, q) = Hcol(idx[p], idx[q]);
                Ss(p, q) = Scol(idx[p], idx[q]);
            }
        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> ges(Hs, Ss);
        if (ges.info() != Eigen::Success) continue;   // skip: keep the atom's block empty
        const Eigen::MatrixXd C = ges.eigenvectors();

        // Aufbau over this atom's orbitals with fractional occupation of the last.
        double z = static_cast<double>(m_atoms[a]);
        for (int k = 0; k < m && z > 0.0; ++k) {
            const double f = std::min(2.0, z);
            z -= f;
            occupied += f;
            for (int p = 0; p < m; ++p)
                for (int q = 0; q < m; ++q)
                    P(idx[p], idx[q]) += f * C(p, k) * C(q, k);
        }
    }

    // Charged systems: scale so the guess carries exactly the right electron count.
    const double nelec = static_cast<double>(m_num_electrons);
    if (occupied > 1.0e-12 && std::abs(occupied - nelec) > 1.0e-12)
        P *= nelec / occupied;
    return P;
}

// Claude Generated (Sep 2026): warm start for the next geometry of the same molecule.
// The basis functions moved with the atoms, so the old occupied orbitals C are no
// longer orthonormal in the new metric S. Symmetric (Loewdin) re-orthonormalisation
//   C' = C (C^T S C)^{-1/2},   P = 2 C' C'^T
// restores Tr(P S) = N and idempotency (P S P = 2 P) exactly, and is the closest
// orthonormal set to C -- a much better start than superposed atomic densities once
// the steps are small (optimisation, MD).
Matrix QMEngine::buildWarmStartGuess() const
{
    const int n_occ = m_num_electrons / 2;
    if (!m_scf_warm_start || m_warm_C.size() == 0 || m_warm_atoms != m_atoms
        || m_warm_C.rows() != m_nbf || m_warm_C.cols() != n_occ || n_occ == 0)
        return Matrix();
    const Eigen::MatrixXd C(m_warm_C);
    const Eigen::MatrixXd M = C.transpose() * Eigen::MatrixXd(m_S) * C;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(M);
    if (es.info() != Eigen::Success || es.eigenvalues().minCoeff() < 1.0e-8)
        return Matrix();  // near-singular: the geometry step was too large, fall back
    const Eigen::MatrixXd Minvsqrt = es.eigenvectors()
        * es.eigenvalues().cwiseSqrt().cwiseInverse().asDiagonal() * es.eigenvectors().transpose();
    const Eigen::MatrixXd Cp = C * Minvsqrt;
    return Matrix(2.0 * Cp * Cp.transpose());
}

Matrix QMEngine::buildInitialGuess() const
{
    m_last_warm_start = false;
    const Matrix warm = buildWarmStartGuess();
    if (warm.size() != 0) {
        m_last_warm_start = true;
        return warm;
    }
    if (m_scf_guess == "h0")
        return Matrix::Zero(m_nbf, m_nbf);   // bare core: first Fock = H
    if (m_scf_guess != "sad" && CurcumaLogger::get_verbosity() >= 1)
        CurcumaLogger::warn(fmt::format(
            "QM: unknown -qm.scf_guess '{}' (use sad|h0); using sad", m_scf_guess));
    return buildAtomicGuess();
}

// =================================================================================
// Fock build: F = H + J - 1/2 K  (closed-shell RHF, spin-summed density)
// =================================================================================

Matrix QMEngine::buildFock(const Matrix& P) const
{
    Matrix J, K;
    buildJK(P, J, K);
    return m_H + J - 0.5 * K;
}

// Claude Generated (Sep 2026): -qm.scf_direct on|off|auto. auto keeps the stored
// tensor (fastest J/K) as long as it fits in -qm.eri_max_memory_mb.
bool QMEngine::useDirectSCF() const
{
    if (m_scf_direct == "on" || m_scf_direct == "true") return true;
    if (m_scf_direct == "off" || m_scf_direct == "false") return false;
    if (m_scf_direct != "auto" && CurcumaLogger::get_verbosity() >= 1)
        CurcumaLogger::warn(fmt::format("QM: unknown -qm.scf_direct '{}' (use auto|on|off); using auto", m_scf_direct));
    const double n = static_cast<double>(m_nbf);
    return n * n * n * n * 8.0 / 1.0e6 > m_eri_max_memory_mb;
}

// J(P), K(P) in the active basis: from the stored tensor, or integral-direct
// (shell-pair tables built on the first call per geometry).
void QMEngine::buildJK(const Matrix& P, Matrix& J, Matrix& K) const
{
    if (!useDirectSCF()) {
        const qmint::ERITensor& eri = eriActive();
        J = qmint::buildCoulomb(eri, P, m_threads);
        K = qmint::buildExchange(eri, P, m_threads);
        return;
    }
    if (!m_direct.ready()) {
        m_direct = qmint::DirectJK(m_gto_basis, m_eri_screening, m_Q.size() != 0 ? &m_Q : nullptr);
        if (CurcumaLogger::get_verbosity() >= 2)
            CurcumaLogger::info(fmt::format(
                "QM: integral-direct J/K ({} basis functions; stored tensor would need {:.1f} MB)",
                m_direct.n(), std::pow((double)m_direct.n(), 4) * 8.0 / 1.0e6));
    }
    m_direct_quartets += m_direct.build(P, J, K, m_threads);
}

// =================================================================================
// Generalized eigenproblem via the Lowdin reduce: F C = S C eps
//   F' = X^T F X  (symmetric),  F' C' = C' eps,  C = X C'.
// =================================================================================

bool QMEngine::solveFock(const Matrix& F, Matrix& C, Vector& eps) const
{
    const Eigen::MatrixXd Fcol(F);          // RowMajor -> ColMajor copy
    const Eigen::MatrixXd Fp = m_X.transpose() * Fcol * m_X;
    Eigen::VectorXd ev;
    Eigen::MatrixXd evecs;
    if (!curcuma::eigsolver::solveSymmetric(Fp, ev, evecs, m_threads))
        return false;
    eps = ev;                                // ascending orbital energies
    C = m_X * evecs;                          // back to the AO basis
    return true;
}

// =================================================================================
// Closed-shell density: P = 2 * sum_{i<n_occ} C_i C_i^T
// =================================================================================

Matrix QMEngine::buildDensity(const Matrix& C, int n_occ) const
{
    const int n = static_cast<int>(C.rows());
    Matrix P = Matrix::Zero(n, n);
    for (int i = 0; i < n_occ; ++i)
        P += 2.0 * (C.col(i) * C.col(i).transpose());
    return P;
}

// =================================================================================
// Lowdin orthonormalizer X = S^{-1/2}
// =================================================================================

void QMEngine::buildOrthonormalizer()
{
    const Eigen::MatrixXd Scol(m_S);         // RowMajor -> ColMajor copy
    Eigen::MatrixXd X;
    if (!lowdinOrthonormalizer(Scol, X)) {
        CurcumaLogger::error("QM SCF: overlap matrix is singular (basis linearly dependent)");
        m_X = Matrix();
        return;
    }
    m_X = X;                                  // store as RowMajor
}

// =================================================================================
// Active-basis ERI (cartesian, or the spherical-5d transform of it)
// =================================================================================

const qmint::ERITensor& QMEngine::eriActive() const
{
    // No d shells / cartesian_d: the cartesian tensor IS the active one -- return
    // it directly instead of keeping a second n^4 copy.
    if (m_Q.size() == 0)
        return cartesianERI();
    if (!m_eri_active_ready) {
        // Spherical: transform every shell quartet as it is computed, so the
        // cartesian tensor (1.7 GB for benzene/def2-SVP) is never formed.
        const auto t0 = std::chrono::steady_clock::now();
        m_eri_active = qmint::buildERI(m_gto_basis, m_threads, m_eri_screening, &m_Q);
        m_eri_active_ready = true;
        if (CurcumaLogger::get_verbosity() >= 2) {
            const int n = m_eri_active.n();
            CurcumaLogger::info(fmt::format(
                "QM: built 4-centre ERI (spherical, {} basis functions, {:.1f} MB) in {:.1f} ms ({} threads)",
                n, (double)n * n * n * n * 8.0 / 1.0e6,
                std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count(),
                m_threads));
        }
    }
    return m_eri_active;
}

// =================================================================================
// Energy components {ET, EV, EJ, Ex, E_elec} for the HF run (Hartree)
// =================================================================================

std::vector<double> QMEngine::energyComponents() const
{
    return { m_et, m_ev, m_ej, m_ex, m_e_elec };
}