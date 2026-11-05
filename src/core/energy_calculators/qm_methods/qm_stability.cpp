/*
 * <Native QM: wavefunction stability analysis of a converged HF/KS SCF>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * A converged SCF is a stationary point of E(kappa), kappa = real occupied-virtual orbital rotations. It is a minimum
 * only if the orbital Hessian is positive definite. For real rotations the Hessian is 4 (A + B) (RHF, spatial orbitals)
 * resp. (A + B) (UHF, spin orbitals) with (Seeger, Pople, J. Chem. Phys. 66, 3045 (1977); Bauernschmitt, Ahlrichs,
 * J. Chem. Phys. 104, 9047 (1996)):
 *
 *   RHF -> RHF (singlet):   (A + B)_ia,jb = d_ij F_ab - d_ab F_ij + 4 (ia|jb) - a_x [(ib|ja) + (ij|ab)] + 4 f_xc
 *   RHF -> UHF (triplet):   (A + B)_ia,jb = d_ij F_ab - d_ab F_ij        - a_x [(ib|ja) + (ij|ab)] + f_xc(triplet)
 *   UHF -> UHF (internal):  (A + B)_ia s,jb t = d_st (d_ij F_ab - d_ab F_ij) + 2 (ia|jb) - d_st a_x [(ib|ja) + (ij|ab)]
 *                                               + 2 f_xc^st
 *
 * with the Fock matrix in the MO basis (canonical orbitals: F_ab - F_ij = eps_a - eps_i). The matrix is never formed:
 * the product with a trial vector x (virtual x occupied) is one J/K build of the symmetric AO matrix
 *   D = C_v x C_o^T + (C_v x C_o^T)^T,
 * (A + B) x = F_vv x - x F_oo + C_v^T G[D] C_o,  G = the Fock response to the density change the rotation makes:
 *   singlet: G = 2 J[D] - a_x K[D] + F_xc[2 D];  triplet: G = -a_x K[D] + dV_xc,a for (P/2 + t D, P/2 - t D);
 *   UHF:     G_s = J[D_a + D_b] - a_x K[D_s] + dV_xc,s[D_a, D_b].
 * F_xc[D] of a closed-shell density is qmxc::buildFxc (exact second derivative); the spin-resolved kernels (triplet,
 * UKS) have no code of their own here and are formed as a central difference of the spin V_xc (step t with
 * t max|D| = 1e-4), so the KS triplet/UKS eigenvalues carry a difference error (about 1e-8 Eh, not measured on every
 * functional). The lowest eigenvalue comes from a Davidson solve with the diagonal F_aa - F_ii as preconditioner.
 *
 * Follow: the orbitals are rotated by U = exp(K), K_ai = -K_ia = v_ai (v = unstable eigenvector, unit norm as PySCF's
 * stability(), closed form through the SVD of v), and the SCF is re-run from that density. RHF -> UHF switches the
 * geometry to UHF/UKS (alpha +v, beta -v).
 *
 * Claude Generated (Nov 2026). AI-generated, machine-tested only (docs/QM_SCF_LEVEL_SHIFT_STABILITY_2026-11.md).
 *
 * This program is free software under GPL-3.0
 */

#include "qm_engine.h"
#include "qm_scf_internal.h"
#include "qm_xc_fxc.h"
#include "qm_xc_spin.h"

#include "src/core/curcuma_logger.h"

#include <Eigen/Dense>
#include <fmt/format.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <functional>
#include <limits>
#include <string>
#include <vector>

namespace {

/// An eigenvalue below this (Eh) counts as an instability (PySCF uses -1e-5 as well).
constexpr double kStabilityThreshold = -1.0e-5;

struct DavidsonResult {
    double value = std::numeric_limits<double>::quiet_NaN();
    Eigen::VectorXd vector;
    int iterations = 0;
    int products = 0;
    double residual = std::numeric_limits<double>::infinity();
    bool converged = false;
};

/**
 * @brief Lowest eigenpair of a symmetric operator by Davidson's method (E. R. Davidson, J. Comput. Phys. 17, 87 (1975)).
 *
 * Start block: unit vectors on the smallest diagonal entries plus one fixed pseudo-random vector (it has components in
 * every symmetry block, so a lowest mode of another irrep than the HOMO -> LUMO excitation is not missed). Correction
 * t = r / (theta - diag), orthogonalised twice against the subspace; collapse to the current vector at 40 vectors.
 * Claude Generated (Nov 2026).
 */
DavidsonResult davidsonLowest(const std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& op, const Eigen::VectorXd& diag,
                              double tol, int max_iter)
{
    DavidsonResult res;
    const int n = static_cast<int>(diag.size());
    if (n == 0) return res;
    const int max_sub = std::min(n, 40);
    std::vector<Eigen::VectorXd> V, AV;
    auto addVector = [&](Eigen::VectorXd t) {
        for (int pass = 0; pass < 2; ++pass)
            for (const auto& v : V) t -= v.dot(t) * v;
        const double nt = t.norm();
        if (nt < 1.0e-8 || static_cast<int>(V.size()) >= n) return false;
        V.push_back(t / nt);
        AV.push_back(op(V.back()));
        ++res.products;
        return true;
    };
    std::vector<int> order(n);
    for (int i = 0; i < n; ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&](int a, int b) { return diag(a) < diag(b); });
    for (int k = 0; k < std::min(n, 3); ++k) addVector(Eigen::VectorXd::Unit(n, order[k]));
    Eigen::VectorXd mix(n);
    for (int i = 0; i < n; ++i) mix(i) = std::sin(12.9898 * (i + 1)) * 0.5 + 0.1;  // deterministic, all components
    addVector(mix);

    for (int it = 0; it < max_iter; ++it) {
        res.iterations = it + 1;
        const int m = static_cast<int>(V.size());
        Eigen::MatrixXd H(m, m);
        for (int i = 0; i < m; ++i)
            for (int j = 0; j <= i; ++j) H(i, j) = H(j, i) = 0.5 * (V[i].dot(AV[j]) + V[j].dot(AV[i]));
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(H);
        const double theta = es.eigenvalues()(0);
        const Eigen::VectorXd y = es.eigenvectors().col(0);
        Eigen::VectorXd x = Eigen::VectorXd::Zero(n), Ax = Eigen::VectorXd::Zero(n);
        for (int i = 0; i < m; ++i) {
            x += y(i) * V[i];
            Ax += y(i) * AV[i];
        }
        const Eigen::VectorXd r = Ax - theta * x;
        res.value = theta;
        res.vector = x / x.norm();
        res.residual = r.norm();
        if (res.residual < tol) {
            res.converged = true;
            break;
        }
        Eigen::VectorXd t(n);
        for (int i = 0; i < n; ++i) {
            const double d = theta - diag(i);
            t(i) = r(i) / (std::abs(d) > 1.0e-4 ? d : (d < 0.0 ? -1.0e-4 : 1.0e-4));
        }
        if (m >= max_sub) {  // collapse: keep the current approximation and its product
            V.assign(1, x / x.norm());
            AV.assign(1, Ax / x.norm());
        }
        if (!addVector(t) && !addVector(r)) {  // no new direction: the subspace is exhausted
            res.converged = res.residual < 1.0e3 * tol;
            break;
        }
    }
    return res;
}

/// Orbitals C (columns: nocc occupied, then the virtuals) rotated by exp(K), K_ai = -K_ia = x_ai (x: nvir x nocc).
/// Closed form via the thin SVD x = U s W^T: the occupied block of exp(K) is W cos(s) W^T + (1 - W W^T), the
/// virtual-occupied block U sin(s) W^T, occupied-virtual -W sin(s) U^T, virtual-virtual U cos(s) U^T + (1 - U U^T).
/// Claude Generated (Nov 2026).
Matrix rotateOrbitals(const Matrix& C, int nocc, const Eigen::MatrixXd& x)
{
    const int n = static_cast<int>(C.cols()), nv = n - nocc;
    if (nocc <= 0 || nv <= 0) return C;
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(x, Eigen::ComputeThinU | Eigen::ComputeThinV);
    const Eigen::MatrixXd& U = svd.matrixU();
    const Eigen::MatrixXd& W = svd.matrixV();
    const Eigen::VectorXd s = svd.singularValues();
    const Eigen::VectorXd cs = s.array().cos().matrix(), sn = s.array().sin().matrix();
    Eigen::MatrixXd R = Eigen::MatrixXd::Identity(n, n);
    R.topLeftCorner(nocc, nocc) += W * (cs.array() - 1.0).matrix().asDiagonal() * W.transpose();
    R.bottomRightCorner(nv, nv) += U * (cs.array() - 1.0).matrix().asDiagonal() * U.transpose();
    R.bottomLeftCorner(nv, nocc) = U * sn.asDiagonal() * W.transpose();
    R.topRightCorner(nocc, nv) = -W * sn.asDiagonal() * U.transpose();
    return Matrix(Eigen::MatrixXd(C) * R);
}

/// Orbital blocks of one spin for the stability operator: C_occ, C_vir and the Fock matrix in the MO basis.
struct Spin {
    Eigen::MatrixXd Co, Cv, Foo, Fvv;
    int no = 0, nv = 0;
    Spin();
    Spin(const Spin&);
    ~Spin();  // out of line: -Winline (as HostDiisStore)
};
Spin::Spin() = default;
Spin::Spin(const Spin&) = default;
Spin::~Spin() = default;

/// Spin density C_occ C_occ^T (no factor 2).
Matrix occupiedProjector(const Matrix& C, int nocc)
{
    if (nocc <= 0) return Matrix::Zero(C.rows(), C.rows());
    const Eigen::MatrixXd Co(C.leftCols(nocc));
    const Eigen::MatrixXd P = Co * Co.transpose();
    return Matrix(0.5 * (P + P.transpose()));
}

}  // namespace

// ---------------------------------------------------------------------------------------------------------------
// XC response
// ---------------------------------------------------------------------------------------------------------------

// Claude Generated (Nov 2026): dV_xc,s = d/dt V_xc,s[Pa + t Da, Pb + t Db] at t = 0 as a central difference of the
// spin-polarised V_xc (CPU, uncached; t max|D| = 1e-4). Used for the RKS triplet and the UKS stability operators.
void QMEngine::xcResponseSpinFD(const Matrix& Pa, const Matrix& Pb, const Matrix& Da, const Matrix& Db, Matrix& dVa,
                                Matrix& dVb) const
{
    const double dmax = std::max(Da.cwiseAbs().maxCoeff(), Db.cwiseAbs().maxCoeff());
    dVa = Matrix::Zero(m_nbf, m_nbf);
    dVb = dVa;
    if (!(dmax > 0.0)) return;
    const double t = 1.0e-4 / dmax;
    qmxc::XCResultSpin p, m;
    qmxc::buildVxcSpin(m_xc_grid, m_xc_ao, Pa + t * Da, Pb + t * Db, m_xc, m_xc_rho_threshold, m_threads, p);
    qmxc::buildVxcSpin(m_xc_grid, m_xc_ao, Pa - t * Da, Pb - t * Db, m_xc, m_xc_rho_threshold, m_threads, m);
    dVa = (p.Va - m.Va) / (2.0 * t);
    dVb = (p.Vb - m.Vb) / (2.0 * t);
}

// ---------------------------------------------------------------------------------------------------------------
// Lowest eigenvalue of one stability operator
// ---------------------------------------------------------------------------------------------------------------

// Claude Generated (Nov 2026): kind 0 = RHF -> RHF (real singlet), 1 = RHF -> UHF (real triplet), 2 = UHF -> UHF
// (real internal). Returns the eigenvalue in Eh; `vec` is the eigenvector (layout: virtual x occupied, column-major,
// UHF: alpha block then beta block), `info` a short ASCII description of the solve.
double QMEngine::stabilityLowest(int kind, Eigen::VectorXd& vec, std::string& info) const
{
    const bool ks = isKohnSham();
    const double ax = ks ? m_xc.exact_exchange : 1.0;
    const bool need_k = ax != 0.0;
    auto makeSpin = [&](const Matrix& C, const Matrix& F, int no) {
        Spin s;
        s.no = no;
        s.nv = static_cast<int>(C.cols()) - no;
        const Eigen::MatrixXd Cc(C), Fc(F);
        s.Co = Cc.leftCols(no);
        s.Cv = Cc.rightCols(s.nv);
        s.Foo = s.Co.transpose() * Fc * s.Co;
        s.Fvv = s.Cv.transpose() * Fc * s.Cv;
        return s;
    };
    std::vector<Spin> sp;
    if (kind == 2) {
        sp.push_back(makeSpin(m_mo, m_fock_a, m_nalpha));
        sp.push_back(makeSpin(m_mo_b, m_fock_b, m_nbeta));
    } else {
        sp.push_back(makeSpin(m_mo, m_fock, m_num_electrons / 2));
    }
    int dim = 0;
    for (const auto& s : sp) dim += s.nv * s.no;
    Eigen::VectorXd diag(dim);
    {
        int off = 0;
        for (const auto& s : sp) {
            for (int i = 0; i < s.no; ++i)
                for (int a = 0; a < s.nv; ++a) diag(off + i * s.nv + a) = s.Fvv(a, a) - s.Foo(i, i);
            off += s.nv * s.no;
        }
    }
    auto aoMatrix = [&](const Spin& s, const double* x) {
        const Eigen::Map<const Eigen::MatrixXd> X(x, s.nv, s.no);
        const Eigen::MatrixXd D = s.Cv * X * s.Co.transpose();
        return Matrix(D + D.transpose());
    };
    auto project = [&](const Spin& s, const double* x, const Matrix& G, double* out) {
        const Eigen::Map<const Eigen::MatrixXd> X(x, s.nv, s.no);
        Eigen::Map<Eigen::MatrixXd> Y(out, s.nv, s.no);
        Y = s.Fvv * X - X * s.Foo + s.Cv.transpose() * Eigen::MatrixXd(G) * s.Co;
    };
    std::function<Eigen::VectorXd(const Eigen::VectorXd&)> op;
    if (kind == 2) {
        op = [&](const Eigen::VectorXd& x) {
            const Spin &a = sp[0], &b = sp[1];
            const Matrix Da = aoMatrix(a, x.data()), Db = aoMatrix(b, x.data() + a.nv * a.no);
            Matrix J, Ka, Kb;
            buildJKSpin(Da, Db, J, Ka, Kb);
            Matrix Ga = J, Gb = J;
            if (need_k) {
                Ga -= ax * Ka;
                Gb -= ax * Kb;
            }
            if (ks) {
                Matrix dVa, dVb;
                xcResponseSpinFD(m_density_a, m_density_b, Da, Db, dVa, dVb);
                Ga += dVa;
                Gb += dVb;
            }
            Eigen::VectorXd y(x.size());
            project(a, x.data(), Ga, y.data());
            project(b, x.data() + a.nv * a.no, Gb, y.data() + a.nv * a.no);
            return y;
        };
    } else {
        op = [&, kind](const Eigen::VectorXd& x) {
            const Spin& s = sp[0];
            const Matrix D = aoMatrix(s, x.data());
            Matrix J, K;
            if (need_k)
                buildJK(D, J, K);
            else
                buildCoulombOnly(D, J);
            Matrix G = Matrix::Zero(m_nbf, m_nbf);
            if (kind == 0) {
                G = 2.0 * J;
                if (need_k) G -= ax * K;
                if (ks) {
                    Matrix Vx;
                    if (!qmxc::buildFxc(m_xc_grid, m_xc_ao, m_density, Matrix(2.0 * D), m_xc, m_xc_rho_threshold, m_threads,
                                        Vx, &m_xc_screen)) {
                        Matrix dVa, dVb;  // no jet for this kernel: difference of the spin V_xc (singlet: Da = Db = D)
                        xcResponseSpinFD(Matrix(0.5 * m_density), Matrix(0.5 * m_density), D, D, dVa, dVb);
                        Vx = dVa;
                    }
                    G += Vx;
                }
            } else {
                if (need_k) G -= ax * K;
                if (ks) {
                    Matrix dVa, dVb;
                    xcResponseSpinFD(Matrix(0.5 * m_density), Matrix(0.5 * m_density), D, Matrix(-D), dVa, dVb);
                    G += dVa;
                }
            }
            Eigen::VectorXd y(x.size());
            project(s, x.data(), G, y.data());
            return y;
        };
    }
    const auto t0 = std::chrono::steady_clock::now();
    const DavidsonResult r = davidsonLowest(op, diag, 1.0e-6, 60);
    vec = r.vector;
    info = fmt::format("{} Davidson iterations, {} operator products, residual {:.1e}{}, {:.1f} ms", r.iterations,
                       r.products, r.residual, r.converged ? "" : " (NOT converged)",
                       std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count());
    return r.value;
}

// ---------------------------------------------------------------------------------------------------------------
// Check and follow
// ---------------------------------------------------------------------------------------------------------------

// Claude Generated (Nov 2026): -qm.scf_stability check|follow after a converged SCF (QMEngine::Calculation).
bool QMEngine::runStabilityAnalysis()
{
    const bool follow = m_scf_stability == "follow";
    const bool ks = isKohnSham();
    m_stability_report = json::object();
    m_stability_report["followed"] = 0;
    std::string why_triplet;
    bool triplet_ok = true;
    if (ks && !qmxc::supportsSpinPolarised(m_xc, why_triplet)) triplet_ok = false;
    // An RHF -> UHF instability can only be followed into an unrestricted SCF the method supports
    const bool may_break_spin = m_reference == "auto" && m_closed_shell_only.empty() && m_functional != QMFunctional::DoubleHybrid
        && triplet_ok;
    if (ks && !ensureXCGrid()) return false;
    const char* names[3] = { "RHF -> RHF (real, singlet)", "RHF -> UHF (real, triplet)", "UHF -> UHF (real, internal)" };
    const char* names_ks[3] = { "RKS -> RKS (real, singlet)", "RKS -> UKS (real, triplet)", "UKS -> UKS (real, internal)" };
    const char* keys[3] ={ "rhf_rhf", "rhf_uhf", "uhf_uhf" };
    for (int cycle = 0;; ++cycle) {
        std::vector<int> kinds;
        if (m_open_shell)
            kinds = { 2 };
        else if (triplet_ok)
            kinds = { 0, 1 };
        else
            kinds = { 0 };
        int unstable_kind = -1;
        Eigen::VectorXd unstable_vec;
        for (int kind : kinds) {
            Eigen::VectorXd v;
            std::string info;
            const double lam = stabilityLowest(kind, v, info);
            m_stability_report[keys[kind]] = lam;
            const bool unstable = std::isfinite(lam) && lam < kStabilityThreshold;
            if (CurcumaLogger::get_verbosity() >= 1) {
                const std::string txt = fmt::format("QM stability {}: lowest eigenvalue of (A + B) {:.8f} Eh ({}){}",
                                                    ks ? names_ks[kind] : names[kind], lam, unstable ? "UNSTABLE" : "stable",
                                                    CurcumaLogger::get_verbosity() >= 2 ? "; " + info : "");
                if (unstable)
                    CurcumaLogger::warn(txt);
                else
                    CurcumaLogger::result(txt);
            }
            if (unstable && unstable_kind < 0) {  // RHF: an internal (singlet) instability is followed first
                unstable_kind = kind;
                unstable_vec = v;
            }
        }
        if (!triplet_ok && !m_open_shell && CurcumaLogger::get_verbosity() >= 2)
            CurcumaLogger::info("QM stability: RKS -> UKS (triplet) not checked: " + why_triplet);
        m_stability_report["stable"] = unstable_kind < 0;
        if (unstable_kind < 0 || !follow) break;
        if (cycle >= m_scf_stability_max_follow) {
            if (CurcumaLogger::get_verbosity() >= 1)
                CurcumaLogger::warn(fmt::format("QM stability: still unstable after {} followed instabilities "
                                                "(-qm.scf_stability_max_follow)", cycle));
            break;
        }
        if (unstable_kind == 1 && !may_break_spin) {
            if (CurcumaLogger::get_verbosity() >= 1)
                CurcumaLogger::warn("QM stability: the RHF -> UHF instability is not followed (needs -qm.reference auto and "
                                    "an open-shell capable method); run with -qm.reference uhf -qm.uhf_guess_mix true instead");
            break;
        }
        // --- follow: rotate along the unstable mode, re-run the SCF from that density ---
        const bool was_open = m_open_shell;
        const Matrix P_old = m_density, Pa_old = m_density_a, Pb_old = m_density_b;
        const double e_old = m_e_elec;
        if (unstable_kind == 2) {
            const int nva = m_nbf - m_nalpha, nvb = m_nbf - m_nbeta;
            const Eigen::MatrixXd xa = Eigen::Map<const Eigen::MatrixXd>(unstable_vec.data(), nva, m_nalpha);
            const Eigen::MatrixXd xb = Eigen::Map<const Eigen::MatrixXd>(unstable_vec.data() + nva * m_nalpha, nvb, m_nbeta);
            m_density_a = occupiedProjector(rotateOrbitals(m_mo, m_nalpha, xa), m_nalpha);
            m_density_b = occupiedProjector(rotateOrbitals(m_mo_b, m_nbeta, xb), m_nbeta);
        } else {
            const int no = m_num_electrons / 2, nv = m_nbf - no;
            const Eigen::MatrixXd x = Eigen::Map<const Eigen::MatrixXd>(unstable_vec.data(), nv, no);
            if (unstable_kind == 0) {
                m_density = 2.0 * occupiedProjector(rotateOrbitals(m_mo, no, x), no);
            } else {  // RHF -> UHF: alpha along +v, beta along -v
                m_density_a = occupiedProjector(rotateOrbitals(m_mo, no, x), no);
                m_density_b = occupiedProjector(rotateOrbitals(m_mo, no, Eigen::MatrixXd(-x)), no);
                m_open_shell = true;
                m_nalpha = m_nbeta = no;
            }
        }
        m_scf_continue = true;
        bool ok = m_open_shell ? runUSCF() : runSCF();
        m_scf_continue = false;
        if (!ok || !(m_e_elec < e_old - 1.0e-10)) {
            if (CurcumaLogger::get_verbosity() >= 1)
                CurcumaLogger::warn(fmt::format("QM stability: following the instability {} (E_elec {:.10f} vs {:.10f} Eh); "
                                                "the previous solution is kept", ok ? "did not lower the energy" : "gave no converged SCF",
                                                m_e_elec, e_old));
            m_open_shell = was_open;
            m_density = P_old;
            m_density_a = Pa_old;
            m_density_b = Pb_old;
            m_scf_continue = true;
            ok = m_open_shell ? runUSCF() : runSCF();
            m_scf_continue = false;
            return ok;
        }
        if (unstable_kind == 1) m_broken_symmetry_serial = m_geom_serial;  // this geometry stays UHF/UKS
        m_stability_report["followed"] = m_stability_report["followed"].get<int>() + 1;
        if (CurcumaLogger::get_verbosity() >= 1)
            CurcumaLogger::result(fmt::format("QM stability: followed {} instability, E_elec {:.10f} -> {:.10f} Eh{}",
                                              keys[unstable_kind], e_old, m_e_elec,
                                              unstable_kind == 1 ? " (now unrestricted)" : ""));
    }
    return true;
}
