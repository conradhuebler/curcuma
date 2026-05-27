/*
 * <xTB SCF driver — Fock builder, eigensolver, population update>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Ported from test_xtb_scf_snapshot.cpp (tblite-validated kernels).
 *
 * Fock matrix assembly:
 *   F = H0 - 0.5 * S * diag(v_ao)  (diagonal expansion of shell potential)
 *
 * For GFN2, multipole contributions are added:
 *   vdp, vqp → F += -0.5 * (vdp contributions via dp_int + vqp contributions via qp_int)
 *
 * Generalized eigenproblem: F C = ε S C  (closed-shell)
 * Mulliken populations: P = 2 * C_occ * C_occ^T
 *
 * Claude Generated. GPL-3.0.
 */

#include "xtb_native.h"
#include "src/core/curcuma_logger.h"

#include <Eigen/Dense>

#include <cmath>

namespace curcuma::xtb {

/* ------------------------------------------------------------------ *
 *  Expand shell + atom potential to AO-resolved v_ao:                *
 *    v_ao(mu) = v_sh[shell(mu)] + v_at[atom(mu)]                     *
 * ------------------------------------------------------------------ */
static void expand_potential(const BasisMap& basis,
                             const Potential& pot,
                             Eigen::VectorXd& v_ao)
{
    v_ao.resize(basis.nao);
    for (int ao = 0; ao < basis.nao; ++ao) {
        const int sh = basis.ao2sh[ao];
        const int at = basis.ao2at[ao];
        v_ao(ao) = pot.v_sh(sh) + pot.v_at(at);
    }
}

/* ------------------------------------------------------------------ *
 *  Build Fock matrix from H0, S, and the current potential.          *
 *                                                                    *
 *  Isotropic part (tblite scf/potential.f90):                        *
 *    F_μν = H0_μν - 0.5 * S_μν * (v_ao(μ) + v_ao(ν))               *
 *                                                                    *
 *  GFN2 multipole (tblite add_vmp_to_h1):                            *
 *    F_μν -= 0.5 * [dp_int[:,μ,ν]·vdp(atom_ν) + dp_int[:,ν,μ]·vdp(atom_μ)
 *                   + qp_int[:,μ,ν]·vqp(atom_ν) + qp_int[:,ν,μ]·vqp(atom_μ)]
 * ------------------------------------------------------------------ */
Matrix XTB::buildFock(const Matrix& H0,
                       const Matrix& S,
                       const Potential& pot) const
{
    const int nao = m_basis.nao;
    Eigen::VectorXd v_ao;
    expand_potential(m_basis, pot, v_ao);

    Matrix F = H0;
    for (int mu = 0; mu < nao; ++mu)
        for (int nu = 0; nu < nao; ++nu)
            F(mu, nu) -= 0.5 * S(mu, nu) * (v_ao(mu) + v_ao(nu));

    // GFN2 multipole Fock contribution (tblite add_vmp_to_h1)
    if (m_method == MethodType::GFN2 && m_mp_initialized) {
        for (int mu = 0; mu < nao; ++mu) {
            const int iat = m_basis.ao2at[mu];
            for (int nu = 0; nu < nao; ++nu) {
                const int jat = m_basis.ao2at[nu];
                double dd = 0.0;
                for (int k = 0; k < 3; ++k)
                    dd += m_dp_int[k](mu, nu) * pot.v_dp(k, jat)
                        + m_dp_int[k](nu, mu) * pot.v_dp(k, iat);
                double qq = 0.0;
                for (int k = 0; k < 6; ++k)
                    qq += m_qp_int[k](mu, nu) * pot.v_qp(k, jat)
                        + m_qp_int[k](nu, mu) * pot.v_qp(k, iat);
                F(mu, nu) -= 0.5 * (dd + qq);
            }
        }
    }
    return F;
}

/* ------------------------------------------------------------------ *
 *  Solve generalized eigenvalue problem: F C = ε S C                 *
 *  using Eigen's GeneralizedSelfAdjointEigenSolver.                  *
 *                                                                    *
 *  Two occupation schemes controlled by m_electronic_temp:           *
 *    temp == 0  → integer closed-shell (2/0 per orbital)             *
 *    temp > 0   → Fermi-Dirac smearing (fractional occupation)       *
 *                                                                    *
 *  Default: 300 K (matches TBLite default).                          *
 * ------------------------------------------------------------------ */
bool XTB::solveEigen(const Matrix& F, const Matrix& S)
{
    const int nao = m_basis.nao;

    Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> solver(F, S);
    if (solver.info() != Eigen::Success) {
        return false;
    }

    m_wfn.eps = solver.eigenvalues();
    m_wfn.C   = solver.eigenvectors();

    if (m_electronic_temp > 0.0) {
        // Fermi-Dirac smearing: bisect for Fermi level, build fractional-occupation density
        // Pattern mirrors GFN2::buildDensityMatrix (gfn2.cpp:607-634).
        const double kT = m_electronic_temp * 3.166808e-6;  // K → Hartree
        const double n_elec = m_wfn.nocc;

        double mu_lo = m_wfn.eps.minCoeff() - 1.0;
        double mu_hi = m_wfn.eps.maxCoeff() + 1.0;
        for (int bisect = 0; bisect < 100; ++bisect) {
            const double mu = 0.5 * (mu_lo + mu_hi);
            double n_sum = 0.0;
            for (int i = 0; i < nao; ++i) {
                const double x = (m_wfn.eps(i) - mu) / kT;
                n_sum += 2.0 / (1.0 + std::exp(std::min(x, 500.0)));
            }
            if (n_sum > n_elec) mu_hi = mu;
            else                mu_lo = mu;
            if (mu_hi - mu_lo < 1e-14) break;
        }
        const double mu_f = 0.5 * (mu_lo + mu_hi);

        m_wfn.P = Matrix::Zero(nao, nao);
        for (int i = 0; i < nao; ++i) {
            const double x   = (m_wfn.eps(i) - mu_f) / kT;
            const double f_i = 2.0 / (1.0 + std::exp(std::min(x, 500.0)));
            if (f_i < 1e-12) continue;
            m_wfn.P += f_i * m_wfn.C.col(i) * m_wfn.C.col(i).transpose();
        }
    } else {
        // Integer occupation: closed-shell, 2 electrons per occupied orbital
        const int nocc = static_cast<int>(std::floor(m_wfn.nocc / 2.0));
        if (nocc < 0 || nocc > nao) return false;

        const Matrix C_occ = m_wfn.C.leftCols(nocc);
        m_wfn.P = 2.0 * C_occ * C_occ.transpose();
    }

    return true;
}

/* ------------------------------------------------------------------ *
 *  Update Mulliken populations from density and overlap.             *
 *    n_sh(s) = sum_{μ∈s, ν} P_μν * S_νμ                             *
 *    q_sh(s) = n0_sh(s) - n_sh(s)                                    *
 *    n_at(i) = sum_{s∈i} n_sh(s)                                     *
 * ------------------------------------------------------------------ */
void XTB::updatePopulations(const Matrix& S)
{
    const int nao = m_basis.nao;
    const int nsh = m_basis.nsh;

    // Shell populations via Mulliken: n_sh = sum_{μ∈sh, ν} P_μν * S_νμ
    Vector n_sh = Vector::Zero(nsh);
    for (int s = 0; s < nsh; ++s) {
        const int ao_start = m_basis.iao_sh[s];
        const int ao_nao   = m_basis.nao_sh[s];
        for (int mu = ao_start; mu < ao_start + ao_nao; ++mu) {
            for (int nu = 0; nu < nao; ++nu) {
                n_sh(s) += m_wfn.P(mu, nu) * S(nu, mu);
            }
        }
    }

    // Compute shell charges relative to reference
    m_wfn.q_sh = m_wfn.n0_sh - n_sh;

    // Accumulate atomic populations and charges
    m_wfn.n_at.setZero(m_atomcount);
    m_wfn.q_at.setZero(m_atomcount);
    for (int s = 0; s < nsh; ++s) {
        const int iat = m_basis.sh2at[s];
        m_wfn.n_at(iat) += n_sh(s);
    }
    for (int i = 0; i < m_atomcount; ++i) {
        m_wfn.q_at(i) = static_cast<double>(m_atoms[i]) - m_wfn.n_at(i);
        // valence-only charge (GFN convention): Z_valence - n_at
        // Actually, GFN uses q = n0_at - n_at (where n0 includes only valence)
        m_wfn.q_at(i) = m_wfn.n0_at(i) - m_wfn.n_at(i);
    }

    // GFN2 atomic multipoles via Mulliken on tblite-convention integrals
    if (m_method == MethodType::GFN2 && m_mp_initialized) {
        m_wfn.dp_at.setZero(3, m_atomcount);
        m_wfn.qp_at.setZero(6, m_atomcount);
        for (int mu = 0; mu < nao; ++mu) {
            const int iat = m_basis.ao2at[mu];
            double d0 = 0, d1 = 0, d2 = 0;
            double q0 = 0, q1 = 0, q2 = 0, q3 = 0, q4 = 0, q5 = 0;
            for (int nu = 0; nu < nao; ++nu) {
                const double Pji = m_wfn.P(nu, mu);
                d0 += Pji * m_dp_int[0](nu, mu);
                d1 += Pji * m_dp_int[1](nu, mu);
                d2 += Pji * m_dp_int[2](nu, mu);
                q0 += Pji * m_qp_int[0](nu, mu);
                q1 += Pji * m_qp_int[1](nu, mu);
                q2 += Pji * m_qp_int[2](nu, mu);
                q3 += Pji * m_qp_int[3](nu, mu);
                q4 += Pji * m_qp_int[4](nu, mu);
                q5 += Pji * m_qp_int[5](nu, mu);
            }
            // Sign: multipoles are valence deviations (like q_sh)
            m_wfn.dp_at(0, iat) -= d0;
            m_wfn.dp_at(1, iat) -= d1;
            m_wfn.dp_at(2, iat) -= d2;
            m_wfn.qp_at(0, iat) -= q0;
            m_wfn.qp_at(1, iat) -= q1;
            m_wfn.qp_at(2, iat) -= q2;
            m_wfn.qp_at(3, iat) -= q3;
            m_wfn.qp_at(4, iat) -= q4;
            m_wfn.qp_at(5, iat) -= q5;
        }
    }
}

/* ------------------------------------------------------------------ *
 *  SCF convergence check.                                            *
 *  Returns true when max |Δq| < threshold and |ΔE| < threshold.      *
 * ------------------------------------------------------------------ */
static bool checkConvergence(const Vector& q_sh_old,
                             const Vector& q_sh_new,
                             double e_old, double e_new,
                             double threshold)
{
    const double dq = (q_sh_new - q_sh_old).cwiseAbs().maxCoeff();
    const double de = std::fabs(e_new - e_old);
    return (dq < threshold && de < threshold);
}

/* ------------------------------------------------------------------ *
 *  buildFockFromPotential()                                          *
 *                                                                    *
 *  H0-free delta-Fock from a potential perturbation δpot:           *
 *    δF_μν = -0.5·S_μν·(δv_ao(μ) + δv_ao(ν))                        *
 *  + GFN2 multipole terms if δpot.v_dp/v_qp are populated.          *
 *                                                                    *
 *  Used in the CPSCF orbital-Hessian to compute (K·z)_ai:           *
 *    δF = buildFockFromPotential(δpot)  → (Cᵀ δF C)_ai             *
 * ------------------------------------------------------------------ */
Matrix XTB::buildFockFromPotential(const Potential& dpot) const
{
    const int nao = m_basis.nao;
    const int nat = m_atomcount;
    Eigen::VectorXd v_ao;
    expand_potential(m_basis, dpot, v_ao);

    Matrix dF = Matrix::Zero(nao, nao);
    for (int mu = 0; mu < nao; ++mu)
        for (int nu = 0; nu < nao; ++nu)
            dF(mu, nu) -= 0.5 * m_S(mu, nu) * (v_ao(mu) + v_ao(nu));

    if (m_method == MethodType::GFN2 && m_mp_initialized
        && dpot.v_dp.cols() == nat && dpot.v_qp.cols() == nat) {
        for (int mu = 0; mu < nao; ++mu) {
            const int iat = m_basis.ao2at[mu];
            for (int nu = 0; nu < nao; ++nu) {
                const int jat = m_basis.ao2at[nu];
                double dd = 0.0;
                for (int k = 0; k < 3; ++k)
                    dd += m_dp_int[k](mu, nu) * dpot.v_dp(k, jat)
                        + m_dp_int[k](nu, mu) * dpot.v_dp(k, iat);
                double qq = 0.0;
                for (int k = 0; k < 6; ++k)
                    qq += m_qp_int[k](mu, nu) * dpot.v_qp(k, jat)
                        + m_qp_int[k](nu, mu) * dpot.v_qp(k, iat);
                dF(mu, nu) -= 0.5 * (dd + qq);
            }
        }
    }
    return dF;
}

/* ------------------------------------------------------------------ *
 *  mullikenFromDensity()                                             *
 *                                                                    *
 *  Compute Mulliken charges and multipoles from a density            *
 *  perturbation δP without modifying any member state.              *
 *                                                                    *
 *    δn_sh(s) = Σ_{μ∈s,ν} δP_μν · S_νμ                              *
 *    δq_sh    = -δn_sh     (n0 constant)                             *
 *    δq_at    = -Σ_{s∈i} δn_sh(s)                                   *
 *    δdp_at, δqp_at  via tblite-convention integrals (GFN2 only)    *
 * ------------------------------------------------------------------ */
void XTB::mullikenFromDensity(const Matrix& dP,
                               Vector& dq_sh, Vector& dq_at,
                               Eigen::MatrixXd& ddp_at,
                               Eigen::MatrixXd& dqp_at) const
{
    const int nao = m_basis.nao;
    const int nsh = m_basis.nsh;
    const int nat = m_atomcount;

    Vector dn_sh = Vector::Zero(nsh);
    for (int s = 0; s < nsh; ++s) {
        const int ao_start = m_basis.iao_sh[s];
        const int ao_nao   = m_basis.nao_sh[s];
        for (int mu = ao_start; mu < ao_start + ao_nao; ++mu)
            for (int nu = 0; nu < nao; ++nu)
                dn_sh(s) += dP(mu, nu) * m_S(nu, mu);
    }
    dq_sh = -dn_sh;

    dq_at = Vector::Zero(nat);
    for (int s = 0; s < nsh; ++s)
        dq_at(m_basis.sh2at[s]) -= dn_sh(s);

    ddp_at = Eigen::MatrixXd::Zero(3, nat);
    dqp_at = Eigen::MatrixXd::Zero(6, nat);
    if (m_method == MethodType::GFN2 && m_mp_initialized) {
        for (int mu = 0; mu < nao; ++mu) {
            const int iat = m_basis.ao2at[mu];
            for (int nu = 0; nu < nao; ++nu) {
                const double dPji = dP(nu, mu);
                for (int k = 0; k < 3; ++k)
                    ddp_at(k, iat) -= dPji * m_dp_int[k](nu, mu);
                for (int k = 0; k < 6; ++k)
                    dqp_at(k, iat) -= dPji * m_qp_int[k](nu, mu);
            }
        }
    }
}

} // namespace curcuma::xtb
