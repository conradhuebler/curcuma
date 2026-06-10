/*
 * D4 single-shot EEQ charge model — implementation
 *
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 * Claude Generated 2026 — AP ∂q/∂x (Phase 2).
 */

#include "d4_charge_model.h"

#include "src/core/energy_calculators/ff_methods/cn_calculator.h"  // covalent radii + CN form
#include "src/core/energy_calculators/ff_methods/gfnff_par.h"      // chi_eeq, gam_eeq, alpha_eeq, cnf_eeq

#include <cmath>

namespace curcuma::dispersion {

namespace {
    constexpr double TSQRT2PI       = 0.797884560802866;   // sqrt(2/π)
    constexpr double TWO_OVER_SQRTPI = 1.1283791670955126; // 2/sqrt(π)
    constexpr double ANG2BOHR        = 1.8897259886;
    constexpr double K_SCALED        = 4.0 / 3.0;          // GFN-FF CN radius scaling
    constexpr double KN              = -7.5;               // GFN-FF erf-CN steepness
    constexpr double CNMAX           = 4.4;                // CN log-compression cap
    constexpr double CN_EPS          = 1.0e-10;            // guard for 1/sqrt(CN)
}

void D4ChargeModel::resolveParams(const std::vector<int>& atoms,
                                  std::vector<double>& chi,
                                  std::vector<double>& gam,
                                  std::vector<double>& alpha_sq,
                                  std::vector<double>& cnf,
                                  std::vector<double>& rcov_bohr)
{
    const int N = static_cast<int>(atoms.size());
    chi.assign(N, 1.0);
    gam.assign(N, 0.0);
    alpha_sq.assign(N, 1.0);
    cnf.assign(N, 0.0);
    rcov_bohr.assign(N, 0.0);
    for (int i = 0; i < N; ++i) {
        const int Z = atoms[i];
        if (Z >= 1 && Z <= 86) {
            chi[i]       = GFNFFParameters::chi_eeq[Z - 1];
            gam[i]       = GFNFFParameters::gam_eeq[Z - 1];
            alpha_sq[i]  = GFNFFParameters::alpha_eeq[Z - 1] * GFNFFParameters::alpha_eeq[Z - 1];
            cnf[i]       = GFNFFParameters::cnf_eeq[Z - 1];
            rcov_bohr[i] = K_SCALED * CNCalculator::getCovalentRadius(Z) * ANG2BOHR;
        }
    }
}

Vector D4ChargeModel::computeCharges(const std::vector<int>& atoms,
                                     const Matrix& geom_bohr,
                                     double total_charge)
{
    const int N = static_cast<int>(atoms.size());
    m_n = N;
    m_atoms = atoms;
    m_geom = geom_bohr;

    if (N == 0) { m_q = Vector(); return m_q; }

    // Resolve per-atom EEQ parameters (angewChem2020 set, α stored squared).
    // resolveParams is the single source of truth shared with the GPU port.
    std::vector<double> chi, gam;
    resolveParams(atoms, chi, gam, m_alp, m_cnf, m_rcov_bohr);

    // Coordination numbers (GFN-FF log-compressed erf form) + raw sum for the
    // log-compression chain factor used in the gradient.
    m_cn.assign(N, 0.0);
    m_cn_raw.assign(N, 0.0);
    const double log1p_ecnmax = std::log(1.0 + std::exp(CNMAX));
    for (int i = 0; i < N; ++i) {
        if (m_rcov_bohr[i] == 0.0) continue;
        double cn_raw = 0.0;
        for (int j = 0; j < N; ++j) {
            if (i == j || m_rcov_bohr[j] == 0.0) continue;
            const double r = (m_geom.row(i) - m_geom.row(j)).norm();
            const double rcov_ij = m_rcov_bohr[i] + m_rcov_bohr[j];
            const double dr = (r - rcov_ij) / rcov_ij;
            cn_raw += 0.5 * (1.0 + std::erf(KN * dr));
        }
        m_cn_raw[i] = cn_raw;
        m_cn[i] = log1p_ecnmax - std::log(1.0 + std::exp(CNMAX - cn_raw));
    }

    // Build augmented EEQ matrix M = [[A, 1], [1^T, 0]] and RHS c = [b; Q].
    const int m = N + 1;
    Matrix M = Matrix::Zero(m, m);
    Vector c = Vector::Zero(m);
    for (int i = 0; i < N; ++i) {
        M(i, i) = gam[i] + TSQRT2PI / std::sqrt(m_alp[i]);
        for (int j = 0; j < i; ++j) {
            const double r = (m_geom.row(i) - m_geom.row(j)).norm();
            const double gammij = 1.0 / std::sqrt(m_alp[i] + m_alp[j]);
            const double aij = std::erf(gammij * r) / r;
            M(i, j) = aij;
            M(j, i) = aij;
        }
        M(i, N) = 1.0;
        M(N, i) = 1.0;
        // b_i = -χ_i + κ_i·sqrt(CN_i)
        c(i) = -chi[i] + m_cnf[i] * std::sqrt(std::max(m_cn[i], 0.0));
    }
    c(N) = total_charge;

    m_lu.compute(M);
    Vector y = m_lu.solve(c);
    m_q = y.head(N);
    return m_q;
}

void D4ChargeModel::addChargeResponseGradient(const Vector& dEdq, Matrix& grad_out) const
{
    if (m_n <= 1) return;  // no pairwise geometry dependence for a single atom
    const int N = m_n;
    if (dEdq.size() != N) return;
    if (grad_out.rows() != N || grad_out.cols() != 3) return;

    // Adjoint (Z-vector): M·z = [dEdq; 0]. M is symmetric, reuse the cached LU.
    Vector rhs = Vector::Zero(N + 1);
    rhs.head(N) = dEdq;
    Vector z = m_lu.solve(rhs);
    const Vector zq = z.head(N);

    // Per-atom b-term weight: u_i = z_q(i)·κ_i/(2√CN_i)·g_i, where g_i is the
    // log-compression factor ∂CN_i/∂cn_raw_i = 1/(1 + e^(cn_raw_i − cnmax)).
    std::vector<double> u(N, 0.0);
    for (int i = 0; i < N; ++i) {
        const double cn = m_cn[i];
        if (cn <= CN_EPS || m_cnf[i] == 0.0) continue;
        const double g = 1.0 / (1.0 + std::exp(m_cn_raw[i] - CNMAX));
        u[i] = zq(i) * m_cnf[i] / (2.0 * std::sqrt(cn)) * g;
    }

    for (int a = 0; a < N; ++a) {
        if (m_rcov_bohr[a] == 0.0 && m_alp[a] <= 0.0) continue;
        for (int b = 0; b < a; ++b) {
            const Eigen::Vector3d rij = m_geom.row(a) - m_geom.row(b);
            const double r = rij.norm();
            if (r < 1.0e-10) continue;
            const Eigen::Vector3d uhat = rij / r;

            // --- A-term: −c_ab · A'(r) · û,  c_ab = z_q(a)q_b + z_q(b)q_a ---
            const double gammij = 1.0 / std::sqrt(m_alp[a] + m_alp[b]);
            const double gr = gammij * r;
            const double Aprime = gammij * TWO_OVER_SQRTPI * std::exp(-gr * gr) / r
                                  - std::erf(gr) / (r * r);
            const double c_ab = zq(a) * m_q(b) + zq(b) * m_q(a);

            // --- b-term (CN): (u_a + u_b) · d(cn_raw)/dr · û ---
            double Draw = 0.0;
            if (m_rcov_bohr[a] > 0.0 && m_rcov_bohr[b] > 0.0) {
                const double rcov_ij = m_rcov_bohr[a] + m_rcov_bohr[b];
                const double dr = (r - rcov_ij) / rcov_ij;
                const double earg = KN * dr;
                Draw = 0.5 * TWO_OVER_SQRTPI * std::exp(-earg * earg) * KN / rcov_ij;
            }

            const double coeff = (u[a] + u[b]) * Draw - c_ab * Aprime;
            const Eigen::Vector3d g_pair = coeff * uhat;
            grad_out.row(a) += g_pair.transpose();
            grad_out.row(b) -= g_pair.transpose();
        }
    }
}

}  // namespace curcuma::dispersion
