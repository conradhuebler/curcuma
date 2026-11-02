/*
 * <Native QM: unrestricted Hartree-Fock and Kohn-Sham (UHF/UKS), stage 5b>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Pople-Nesbet equations (J. A. Pople, R. K. Nesbet, J. Chem. Phys. 22, 571 (1954)): one set of orbitals per spin,
 *   P_s = sum_{i <= n_s} C_s,i C_s,i^T                      (s = alpha, beta; no factor 2)
 *   F_s = h + J[P_a + P_b] - a_x K[P_s] + V_xc,s           (UHF: a_x = 1, no V_xc)
 *   E   = Tr(P_t h) + 1/2 Tr(P_t J[P_t]) - (a_x/2) sum_s Tr(P_s K[P_s]) + E_xc[rho_a, rho_b] + E_nn,  P_t = P_a + P_b
 *       = 1/2 sum_s Tr(P_s (h + F_s))   for UHF.
 * Each spin solves F_s C_s = S C_s eps_s (Lowdin reduction, the same solveFock as RHF). For P_a = P_b = P/2 every
 * expression is the RHF one of qm_scf.cpp.
 *
 * SCF driver (runUSCF): the RHF loop's building blocks (qm_scf_internal.h): Pulay DIIS / ADIIS with ONE coefficient
 * vector for both spins, error vector (e_a, e_b), e_s = X^T (F_s P_s S - S P_s F_s) X; the convergence test of
 * -qm.scf_convergence auto with |dE|, ||(e_a, e_b)||_F and ||(dP_a, dP_b)||_F. Every Fock build is a full J/K build
 * (no incremental integral-direct builds, no verified stop, no device loop: the open-shell path is CPU only).
 *
 * Spin contamination (Szabo-Ostlund eq. 3.226 in this notation):
 *   <S^2> = S_z (S_z + 1) + n_beta - Tr(P_a S P_b S),   S_z = (n_alpha - n_beta) / 2.
 *
 * Analytic gradient (computeGradientUHF): the RHF terms with P = P_a + P_b and W = sum_s sum_i eps_s,i C C^T in the
 * one-electron part (incl. ECP), the UHF two-electron weight (qmint::gradientTwoElectronUHF), and for UKS the spin XC
 * gradient (qmxc::buildXCGradientSpin) plus the grid response (Becke weight derivative with the total energy density).
 *
 * Claude Generated (Nov 2026, stage 5b). AI-generated, machine-tested only (docs/NATIVE_QM_IMPLEMENTATION.md, stage 5b).
 *
 * This program is free software under GPL-3.0
 */

#include "qm_citations.h"
#include "qm_engine.h"
#include "qm_ri.h"
#include "qm_scf_internal.h"

#include "src/core/curcuma_logger.h"
#include "src/core/gpu_fallback.h"
#include "src/core/units.h"

#include "native_eigensolver.h"

#include <Eigen/Dense>
#include <fmt/format.h>

#include <chrono>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace {

/// Orbital gradient e = X^T (F P S - S P F) X of one spin (zero at a stationary density).
Matrix spinOrbitalGradient(const Matrix& F, const Matrix& P, const Matrix& S, const Matrix& X, int threads)
{
    const Eigen::MatrixXd Fc(F), Pc(P), Sc(S), Xc(X);
    const Eigen::MatrixXd FPS = gemmThreaded(gemmThreaded(Fc, false, Pc, threads), false, Sc, threads);
    const Eigen::MatrixXd comm = FPS - FPS.transpose();
    return Matrix(gemmThreaded(gemmThreaded(Xc, true, comm, threads), false, Xc, threads));
}

/// Spin density P = C_occ C_occ^T (no factor 2).
Matrix spinDensity(const Matrix& C, int nocc, int threads)
{
    if (nocc <= 0) return Matrix::Zero(C.rows(), C.rows());
    const Eigen::MatrixXd CoT = Eigen::MatrixXd(C.leftCols(nocc)).transpose();
    const Eigen::MatrixXd P = gemmThreaded(CoT, true, CoT, threads);
    return Matrix(0.5 * (P + P.transpose()));
}

/// Lowdin re-orthonormalisation of the previous geometry's occupied orbitals, P = C' C'^T (as buildWarmStartGuess).
Matrix reorthonormalisedDensity(const Matrix& Cw, const Matrix& S)
{
    if (Cw.cols() == 0) return Matrix::Zero(S.rows(), S.rows());
    const Eigen::MatrixXd C(Cw);
    const Eigen::MatrixXd M = C.transpose() * Eigen::MatrixXd(S) * C;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(M);
    if (es.info() != Eigen::Success || es.eigenvalues().minCoeff() < 1.0e-8) return Matrix();
    const Eigen::MatrixXd Mis = es.eigenvectors() * es.eigenvalues().cwiseSqrt().cwiseInverse().asDiagonal()
        * es.eigenvectors().transpose();
    const Eigen::MatrixXd Cp = C * Mis;
    return Matrix(Cp * Cp.transpose());
}

/// DIIS history of the UHF loop: both spins per slot, one coefficient vector.
struct UhfDiisStore {
    std::vector<Matrix> Fa, Fb, Pa, Pb, ea, eb;
    explicit UhfDiisStore(int slots) : Fa(slots), Fb(slots), Pa(slots), Pb(slots), ea(slots), eb(slots) {}
    ~UhfDiisStore();  // out of line: -Winline (as HostDiisStore)
    Eigen::MatrixXd gram(const std::vector<int>& sl) const
    {
        const int m = static_cast<int>(sl.size());
        Eigen::MatrixXd G(m, m);
        for (int i = 0; i < m; ++i)
            for (int j = 0; j <= i; ++j)
                G(i, j) = G(j, i) = ea[sl[i]].cwiseProduct(ea[sl[j]]).sum() + eb[sl[i]].cwiseProduct(eb[sl[j]]).sum();
        return G;
    }
    /// ADIIS model summed over both spins: dE = sum_s Tr(F_s dP_s).
    void adiisData(const std::vector<int>& sl, Eigen::VectorXd& a, Eigen::MatrixXd& M) const
    {
        const int m = static_cast<int>(sl.size());
        const int nw = sl.back();
        a.resize(m);
        M.resize(m, m);
        for (int i = 0; i < m; ++i) {
            const Matrix dPa = Pa[sl[i]] - Pa[nw], dPb = Pb[sl[i]] - Pb[nw];
            a(i) = dPa.cwiseProduct(Fa[nw]).sum() + dPb.cwiseProduct(Fb[nw]).sum();
            for (int j = 0; j < m; ++j)
                M(i, j) = dPa.cwiseProduct(Fa[sl[j]] - Fa[nw]).sum() + dPb.cwiseProduct(Fb[sl[j]] - Fb[nw]).sum();
        }
    }
    void combine(const std::vector<int>& sl, const Eigen::VectorXd& c, Matrix& FaOut, Matrix& FbOut) const
    {
        FaOut = Matrix::Zero(Fa[sl[0]].rows(), Fa[sl[0]].cols());
        FbOut = FaOut;
        for (int i = 0; i < c.size(); ++i) {
            FaOut += c(i) * Fa[sl[i]];
            FbOut += c(i) * Fb[sl[i]];
        }
    }
};

UhfDiisStore::~UhfDiisStore() = default;

/// PySCF init_guess_mixed: rotate the HOMO/LUMO pair of one spin by `angle` (radians).
void mixHomoLumo(Matrix& C, int nocc, double angle)
{
    if (nocc <= 0 || nocc >= C.cols()) return;
    const Eigen::VectorXd h = C.col(nocc - 1), l = C.col(nocc);
    const double c = std::cos(angle), s = std::sin(angle);
    C.col(nocc - 1) = c * h + s * l;
    C.col(nocc) = -s * h + c * l;
}

}  // namespace

// =================================================================================
// Spin state: electrons, -spin, -qm.reference
// =================================================================================

bool QMEngine::resolveSpinState()
{
    m_spin_error.clear();
    const int N = m_num_electrons;
    const long s = std::lround(m_spin);
    auto fail = [&](const std::string& why) {
        m_spin_error = why;
        m_open_shell = false;
        CurcumaLogger::error("QM: " + why);
        return false;
    };
    if (m_reference != "auto" && m_reference != "rhf" && m_reference != "uhf")
        return fail("unknown -qm.reference '" + m_reference + "' (auto, rhf or uhf)");
    if (s < 0) return fail(fmt::format("-spin {} is negative (number of unpaired electrons, multiplicity - 1)", s));
    if (s > N || ((N - s) % 2) != 0)
        return fail(fmt::format("{} electrons and -spin {} (unpaired electrons = multiplicity - 1) do not fit: n_alpha - "
                                "n_beta must have the parity of the electron count{}",
                                N, s, (s == 0 && N % 2 != 0) ? "; for a radical pass -spin 1 (doublet)" : ""));
    if (m_reference == "rhf" && s > 0)
        return fail(fmt::format("-qm.reference rhf with -spin {}: ROHF is not implemented; use -qm.reference auto|uhf "
                                "(unrestricted)", s));
    const bool open = m_reference == "uhf" || (m_reference == "auto" && s > 0);
    if (open && !m_closed_shell_only.empty())
        return fail(m_closed_shell_only + " is closed shell only (RHF/RKS reference); open-shell input (-spin "
                    + std::to_string(s) + (m_reference == "uhf" ? ", -qm.reference uhf" : "")
                    + ") is rejected: no unrestricted " + m_closed_shell_only + " is implemented");
    if (open && m_functional == QMFunctional::DoubleHybrid)
        return fail("double hybrids are closed shell only; open-shell input is rejected");
    if (open && isKohnSham()) {
        std::string why;
        if (!qmxc::supportsSpinPolarised(m_xc, why))
            return fail("open-shell Kohn-Sham (UKS) is not implemented for " + m_xc.name + ": " + why
                        + "; closed-shell (-spin 0) runs are not affected");
    }
    m_open_shell = open;
    m_nalpha = static_cast<int>((N + s) / 2);
    m_nbeta = static_cast<int>((N - s) / 2);
    if (open && !(m_gpu == "none" || m_gpu == "cpu" || m_gpu.empty()) && !m_uhf_gpu_reported) {
        // GPU stage 5b: no open-shell device path; every part runs on the CPU (gpuRequested() is false), counted once
        m_uhf_gpu_reported = true;
        curcuma::reportGpuFallback("qm_uhf", "open-shell UHF/UKS has no device path; SCF, J/K, XC and gradient on the CPU");
    }
    return true;
}

// =================================================================================
// J[Pa + Pb], K[Pa], K[Pb]
// =================================================================================

void QMEngine::buildJKSpin(const Matrix& Pa, const Matrix& Pb, Matrix& J, Matrix& Ka, Matrix& Kb) const
{
    if (isKohnSham() && m_xc.exact_exchange == 0.0) {
        buildCoulombOnly(Pa + Pb, J);
        Ka = Matrix::Zero(m_nbf, m_nbf);
        Kb = Ka;
        return;
    }
    if (!riActive() && !useDirectSCF()) {
        // stored host tensor: J once from the total density, K per spin
        const qmint::ERITensor& eri = eriActive();
        J = qmint::buildCoulomb(eri, Pa + Pb, m_threads);
        Ka = qmint::buildExchange(eri, Pa, m_threads);
        Kb = qmint::buildExchange(eri, Pb, m_threads);
        return;
    }
    // integral-direct and RI: two builds, J is linear in P (J = J[Pa] + J[Pb])
    Matrix Ja, Jb;
    buildJK(Pa, Ja, Ka);
    buildJK(Pb, Jb, Kb);
    J = Ja + Jb;
}

const qmxc::XCResultSpin& QMEngine::xcForSpin(const Matrix& Pa, const Matrix& Pb) const
{
    if (m_xcs_last_valid && m_xcs_last_Pa.rows() == Pa.rows() && (m_xcs_last_Pa.array() == Pa.array()).all()
        && (m_xcs_last_Pb.array() == Pb.array()).all())
        return m_xcs_last;
    if (!ensureXCGrid()) {
        m_xcs_last = qmxc::XCResultSpin();
        m_xcs_last.Va = Matrix::Constant(m_nbf, m_nbf, std::numeric_limits<double>::quiet_NaN());
        m_xcs_last.Vb = m_xcs_last.Va;
        m_xcs_last.exc = std::numeric_limits<double>::quiet_NaN();
        m_xcs_last_valid = false;
        return m_xcs_last;
    }
    qmxc::buildVxcSpin(m_xc_grid, m_xc_ao, Pa, Pb, m_xc, m_xc_rho_threshold, m_threads, m_xcs_last);
    m_xcs_last_Pa = Pa;
    m_xcs_last_Pb = Pb;
    m_xcs_last_valid = true;
    ++m_xc_builds;
    m_xc_ms += m_xcs_last.ms;
    return m_xcs_last;
}

// =================================================================================
// UHF / UKS SCF
// =================================================================================

bool QMEngine::runUSCF()
{
    const int n = m_nbf;
    const int na = m_nalpha, nb = m_nbeta;
    if (n <= 0) {
        CurcumaLogger::error("QM UHF: no basis functions (InitialiseMolecule failed?)");
        return false;
    }
    if (na > n || nb > n) {
        CurcumaLogger::error(fmt::format("QM UHF: {} alpha / {} beta electrons but only {} basis functions", na, nb, n));
        return false;
    }
    buildOrthonormalizer();
    if (m_X.size() == 0) {
        CurcumaLogger::error("QM UHF: S^{-1/2} construction failed (linear basis?)");
        return false;
    }
    if (riActive() && !ensureRI()) return false;
    const bool ks = isKohnSham();
    if (ks) {
        if (!ensureXCGrid()) return false;
        qmcite::kohnSham(m_xc, m_grid_options.radial == "treutler");
    }
    qmcite::uhf(ks);
    const bool ri_jk = riActive() && m_ri_mode == "jk";
    const bool direct = !ri_jk && useDirectSCF();
    if (direct) qmcite::directScf();
    if (m_scf_mode == "adiis") CitationRegistry::cite("adiis");
    m_ri_occ_factor = Eigen::MatrixXd();  // the RHF orbital factor does not belong to a spin density
    m_ri_occ_density = Matrix();
    m_direct_thr_now = -1.0;
    m_direct_fp32_now = false;
    m_direct_quartets = 0;
    const double ax = ks ? m_xc.exact_exchange : 1.0;

    // Start: tightening a converged SCF of this geometry, the previous geometry's spin orbitals, else the RHF guess
    // density split evenly over the spins (PySCF's UHF start for every init_guess).
    Matrix Pa, Pb;
    m_last_warm_start = false;
    if (m_scf_continue && m_density_a.rows() == n && m_density_b.rows() == n) {
        Pa = m_density_a;
        Pb = m_density_b;
        m_last_warm_start = true;
    } else if (m_scf_warm_start && m_warm_atoms == m_atoms && m_warm_Ca.rows() == n && m_warm_Ca.cols() == na
               && m_warm_Cb.rows() == n && m_warm_Cb.cols() == nb) {
        Pa = reorthonormalisedDensity(m_warm_Ca, m_S);
        Pb = reorthonormalisedDensity(m_warm_Cb, m_S);
        m_last_warm_start = Pa.size() != 0 && Pb.size() != 0;
    }
    if (!m_last_warm_start) {
        const Matrix P0 = buildInitialGuess();  // RHF slots (warm C, density history) are never filled by UHF
        m_last_warm_start = false;
        Pa = 0.5 * P0;
        Pb = 0.5 * P0;
    }

    const bool use_diis = (m_scf_mode != "plain");
    const double damping = use_diis ? 0.0 : 1.0 - m_scf_damping;
    QMScfAccelerator diis(m_diis_subspace, m_scf_mode == "adiis");
    UhfDiisStore store(diis.capacity());
    ScfConvergence conv = makeConvergence(m_scf_convergence == "density", m_scf_threshold, m_scf_energy_threshold,
                                          m_scf_tight_now ? m_scf_gradient_threshold_tight : m_scf_gradient_threshold);
    const auto t_scf0 = std::chrono::steady_clock::now();
    double t_fock_ms = 0.0, t_eig_ms = 0.0;
    m_scf_converged = false;
    m_scf_iterations = 0;

    auto fockPair = [&](const Matrix& A, const Matrix& B, Matrix& Fa, Matrix& Fb, double& E) {
        Matrix J, Ka, Kb;
        buildJKSpin(A, B, J, Ka, Kb);
        Fa = m_H + J;
        Fb = m_H + J;
        if (ax != 0.0) {
            Fa -= ax * Ka;
            Fb -= ax * Kb;
        }
        E = 0.5 * ((A + B).cwiseProduct(m_H).sum() + A.cwiseProduct(Fa).sum() + B.cwiseProduct(Fb).sum());
        if (ks) {
            const qmxc::XCResultSpin& xc = xcForSpin(A, B);
            E += xc.exc;  // V_xc enters F only (1/2 Tr(P V_xc) is not part of E)
            Fa += xc.Va;
            Fb += xc.Vb;
        }
        return std::make_pair(J, std::make_pair(Ka, Kb));
    };

    Matrix Ca, Cb;
    Vector eps_a, eps_b;
    for (int iter = 0; iter < m_scf_max_iter; ++iter) {
        m_scf_iterations = iter + 1;
        const auto t_f0 = std::chrono::steady_clock::now();
        Matrix Fa, Fb;
        double E = 0.0;
        fockPair(Pa, Pb, Fa, Fb, E);
        t_fock_ms += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t_f0).count();

        Matrix Fa_use = Fa, Fb_use = Fb;
        double grad_norm = 0.0;
        if (use_diis && iter >= m_diis_start) {
            const int slot = diis.push(true);
            store.ea[slot] = spinOrbitalGradient(Fa, Pa, m_S, m_X, m_threads);
            store.eb[slot] = spinOrbitalGradient(Fb, Pb, m_S, m_X, m_threads);
            store.Fa[slot] = Fa;
            store.Fb[slot] = Fb;
            store.Pa[slot] = Pa;
            store.Pb[slot] = Pb;
            diis.setLastError(std::max(store.ea[slot].cwiseAbs().maxCoeff(), store.eb[slot].cwiseAbs().maxCoeff()));
            grad_norm = std::sqrt(store.ea[slot].squaredNorm() + store.eb[slot].squaredNorm());
            Eigen::MatrixXd G, M;
            Eigen::VectorXd a;
            if (diis.needsDiis()) G = store.gram(diis.slots());
            if (diis.needsAdiis()) store.adiisData(diis.slots(), a, M);
            store.combine(diis.slots(), diis.coefficients(G, a, M), Fa_use, Fb_use);
        } else if (!conv.density_only) {
            grad_norm = std::sqrt(spinOrbitalGradient(Fa, Pa, m_S, m_X, m_threads).squaredNorm()
                                  + spinOrbitalGradient(Fb, Pb, m_S, m_X, m_threads).squaredNorm());
        }
        conv.observe(E, grad_norm);

        const auto t_e0 = std::chrono::steady_clock::now();
        if (!solveFock(Fa_use, Ca, eps_a) || !solveFock(Fb_use, Cb, eps_b)) {
            CurcumaLogger::error(fmt::format("QM UHF: diagonalization failed at iteration {}", iter + 1));
            return false;
        }
        t_eig_ms += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t_e0).count();
        if (m_uhf_guess_mix && iter == 0 && !m_last_warm_start) {
            mixHomoLumo(Ca, na, 0.7853981633974483);  // +45 degrees
            mixHomoLumo(Cb, nb, -0.7853981633974483);
        }
        const Matrix Pa_new = spinDensity(Ca, na, m_threads), Pb_new = spinDensity(Cb, nb, m_threads);
        const double delta_P = std::sqrt((Pa_new - Pa).squaredNorm() + (Pb_new - Pb).squaredNorm());

        if (CurcumaLogger::get_verbosity() >= 1)
            CurcumaLogger::result(fmt::format("UHF iter {:3d}: E_elec = {:.10f} Eh, dE = {:.2e}, |g| = {:.2e}, dP = {:.3e}",
                                              iter + 1, conv.e_prev, conv.have_prev && iter > 0 ? conv.dE : 0.0,
                                              conv.grad, delta_P));
        if (conv.converged(delta_P)) {
            m_scf_last_dE = conv.dE;
            m_scf_last_grad = conv.grad;
            // final full build of the converged densities: Fock matrices and energy parts
            Matrix Fa_f, Fb_f;
            double E_f = 0.0;
            auto jk = fockPair(Pa_new, Pb_new, Fa_f, Fb_f, E_f);
            const Matrix& J = jk.first;
            const Matrix& Ka = jk.second.first;
            const Matrix& Kb = jk.second.second;
            m_density_a = Pa_new;
            m_density_b = Pb_new;
            m_density = Pa_new + Pb_new;
            m_fock_a = Fa_f;
            m_fock_b = Fb_f;
            m_fock = Fa_f;
            m_scf_converged = true;
            m_scf_tight_done = m_scf_tight_now;
            m_et = m_density.cwiseProduct(m_T).sum();
            m_ev = m_density.cwiseProduct(m_V).sum();
            m_ej = 0.5 * m_density.cwiseProduct(J).sum();
            m_ex = -0.5 * ax * (Pa_new.cwiseProduct(Ka).sum() + Pb_new.cwiseProduct(Kb).sum());
            m_e_elec = E_f;
            if (ks) {
                const qmxc::XCResultSpin& xc = xcForSpin(Pa_new, Pb_new);
                m_exc = xc.exc;
                m_exc_x = xc.ex;
                m_exc_c = xc.ec;
                m_exc_nelec = xc.nelec;
            }
            m_mo = Ca;
            m_energies = eps_a;
            m_mo_b = Cb;
            m_energies_b = eps_b;
            m_warm_Ca = Ca.leftCols(na);
            m_warm_Cb = Cb.leftCols(nb);
            m_warm_atoms = m_atoms;
            const double sz = 0.5 * (na - nb);
            m_s2_exact = sz * (sz + 1.0);
            const Matrix PaS = Pa_new * m_S, PbS = Pb_new * m_S;
            m_s2 = m_s2_exact + nb - PaS.cwiseProduct(PbS.transpose()).sum();
            if (const char* tr = std::getenv("CURCUMA_QM_SCF_TRACE"); tr && tr[0] == '1')
                fmt::print(stderr, "QM SCF trace: loop=host-uhf n={} iterations={} start={} E_elec={:.10f}\n", m_nbf,
                           iter + 1, m_last_warm_start ? "warm" : m_scf_guess, m_e_elec);
            if (CurcumaLogger::get_verbosity() >= 2) {
                CurcumaLogger::success(fmt::format(
                    "QM {} SCF converged in {} iterations (dP = {:.3e}, dE = {:.1e}, |g| = {:.1e}; n_alpha {}, n_beta {}, "
                    "start: {})",
                    ks ? "UKS-" + m_xc.name : std::string("UHF"), iter + 1, delta_P, m_scf_last_dE, m_scf_last_grad, na, nb,
                    m_scf_continue ? "tightened from the converged densities" : (m_last_warm_start ? "previous geometry" : m_scf_guess)));
                CurcumaLogger::param("<S^2>", fmt::format("{:.8f} (exact S(S+1) = {:.4f}, spin contamination {:.3e})", m_s2,
                                                          m_s2_exact, m_s2 - m_s2_exact));
                CurcumaLogger::info(fmt::format("QM: UHF loop {:.1f} ms, of which Fock (J/K{}) builds {:.1f} ms, eigensolve {:.1f} ms{}",
                                                std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t_scf0).count(),
                                                ks ? ", V_xc" : "", t_fock_ms, t_eig_ms,
                                                direct ? fmt::format(" (integral-direct, {} shell quartets)", m_direct_quartets)
                                                       : std::string()));
            }
            return true;
        }
        Pa = damping * Pa + (1.0 - damping) * Pa_new;
        Pb = damping * Pb + (1.0 - damping) * Pb_new;
    }
    // not converged: energy of the last densities for diagnostics
    Matrix Fa, Fb;
    double E = 0.0;
    fockPair(Pa, Pb, Fa, Fb, E);
    m_density_a = Pa;
    m_density_b = Pb;
    m_density = Pa + Pb;
    m_fock_a = Fa;
    m_fock_b = Fb;
    m_fock = Fa;
    m_e_elec = E;
    return false;
}

// =================================================================================
// UHF / UKS analytic gradient
// =================================================================================

bool QMEngine::computeGradientUHF()
{
    const int na = m_nalpha, nb = m_nbeta;
    Matrix W = Matrix::Zero(m_nbf, m_nbf);
    for (int i = 0; i < na; ++i) W += m_energies(i) * (m_mo.col(i) * m_mo.col(i).transpose());
    for (int i = 0; i < nb; ++i) W += m_energies_b(i) * (m_mo_b.col(i) * m_mo_b.col(i).transpose());
    Matrix Pc = m_density, Wc = W, Pac = m_density_a, Pbc = m_density_b;
    if (m_Q.size() != 0) {
        Pc = m_Q * m_density * m_Q.transpose();
        Wc = m_Q * W * m_Q.transpose();
        Pac = m_Q * m_density_a * m_Q.transpose();
        Pbc = m_Q * m_density_b * m_Q.transpose();
    }
    const Matrix posBohr = atomPositionsBohr();
    const bool ks = isKohnSham();
    const double kx = ks ? m_xc.exact_exchange : 1.0;
    const auto t0 = std::chrono::steady_clock::now();
    Matrix g1 = qmint::gradientOneElectron(m_gto_basis, m_core_charges, posBohr, Pc, Wc, m_threads);
    if (hasECP()) g1 += ecpGradient(Pc);
    const auto t1 = std::chrono::steady_clock::now();
    Matrix g2;
    std::string info = "CPU";
    if (riActive()) {
        if (!riGradientTwoElectronUHF(Pac, Pbc, g2, info)) {
            m_gradient = Matrix::Zero(m_atomcount, 3);
            m_gradient_failed = true;
            m_gradient_error = info.empty() ? std::string("RI gradient failed") : info;
            CurcumaLogger::error("QM: " + m_gradient_error);
            return false;
        }
    } else {
        g2 = qmint::gradientTwoElectronUHF(m_gto_basis, Pac, Pbc, kx, m_atomcount, m_threads);
    }
    const auto t2 = std::chrono::steady_clock::now();
    const Matrix gn = calculateCoreRepulsionGradient();
    if (!ks) {
        m_gradient = g1 + g2 + gn;
        m_gradient_parts = { g1, g2, gn };
    } else {
        qmcite::ksGradient(m_grid_response);
        if (!ensureXCGrid()) {
            m_gradient = Matrix::Zero(m_atomcount, 3);
            m_gradient_failed = true;
            m_gradient_error = "XC gradient failed: no DFT grid";
            return false;
        }
        const std::vector<int> bf_atom = activeFunctionAtoms();
        qmxc::XCGradient xg;
        qmxc::buildXCGradientSpin(m_xc_grid, m_xc_ao, m_density_a, m_density_b, m_xc, m_xc_rho_threshold, m_threads,
                                  bf_atom, m_atomcount, xg);
        Matrix gresp = Matrix::Zero(m_atomcount, 3);
        if (m_grid_response) {
            std::string why;
            if (m_xc_frozen_grid) why = "grid response with a frozen grid (validation hook) is undefined; use -qm.grid_response false";
            Matrix gw = Matrix::Zero(m_atomcount, 3);
            if (!why.empty() || !qmgrid::gridWeightGradient(m_atoms, posBohr, m_grid_options, m_xc_grid, xg.edens, gw, why,
                                                            m_threads)) {
                m_gradient = Matrix::Zero(m_atomcount, 3);
                m_gradient_failed = true;
                m_gradient_error = "XC gradient failed: " + why;
                CurcumaLogger::error("QM: " + m_gradient_error);
                return false;
            }
            gresp = xg.points + gw;
        }
        m_gradient = g1 + g2 + gn + xg.basis + gresp;
        m_gradient_parts = { g1, g2, gn, xg.basis, gresp };
    }
    if (CurcumaLogger::get_verbosity() >= 2) {
        auto ms = [](auto a, auto b) { return std::chrono::duration<double, std::milli>(b - a).count(); };
        CurcumaLogger::info(fmt::format("QM: analytic {} gradient |g| = {:.6e} Eh/Bohr (1e {:.1f} ms; 2e {:.1f} ms {}; {}{} threads, CPU)",
                                        ks ? "UKS" : "UHF", m_gradient.norm(), ms(t0, t1), ms(t1, t2), info,
                                        ks ? fmt::format("XC {:.1f} ms; ", ms(t2, std::chrono::steady_clock::now())) : std::string(),
                                        m_threads));
    }
    return true;
}
