/*
 * D4 single-shot EEQ charge model — implementation
 *
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 * Claude Generated 2026 — AP ∂q/∂x (Phase 2).
 */

#include "src/core/units.h"
#include "d4_charge_model.h"

#include "src/core/energy_calculators/ff_methods/cn_calculator.h"  // covalent radii + CN form
#include "src/core/energy_calculators/ff_methods/gfnff_par.h"      // chi_eeq, gam_eeq, alpha_eeq, cnf_eeq
#include "src/core/math_compat.h"
#include "d4_eeq_tables.h"  // EEQ(2019) parameters (WP12 stage 4)
#include "d4_ncoord.h"      // D3 covalent radii x 4/3 (Bohr), shared with the D4 CN

#include <cmath>

namespace curcuma::dispersion {

namespace {
    constexpr double TSQRT2PI       = 0.797884560802866;   // sqrt(2/π)
    constexpr double TWO_OVER_SQRTPI = 1.1283791670955126; // 2/sqrt(π)
    constexpr double ANG2BOHR        = CurcumaUnit::Length::angstrom_to_bohr_or_legacy(1.8897259886);
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
    if (m_model == Model::Eeq2019) return computeChargesEeq2019(atoms, geom_bohr, total_charge);
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
    const double log1p_ecnmax = curcuma_log(1.0 + curcuma_exp(CNMAX));
    for (int i = 0; i < N; ++i) {
        if (m_rcov_bohr[i] == 0.0) continue;
        double cn_raw = 0.0;
        for (int j = 0; j < N; ++j) {
            if (i == j || m_rcov_bohr[j] == 0.0) continue;
            const double r = (m_geom.row(i) - m_geom.row(j)).norm();
            const double rcov_ij = m_rcov_bohr[i] + m_rcov_bohr[j];
            const double dr = (r - rcov_ij) / rcov_ij;
            cn_raw += 0.5 * (1.0 + curcuma_erf(KN * dr));
        }
        m_cn_raw[i] = cn_raw;
        m_cn[i] = log1p_ecnmax - curcuma_log(1.0 + curcuma_exp(CNMAX - cn_raw));
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
            const double aij = curcuma_erf(gammij * r) / r;
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
    if (m_model == Model::Eeq2019) { addChargeResponseGradientEeq2019(dEdq, grad_out); return; }
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
                                  - curcuma_erf(gr) / (r * r);
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

// ---------------------------------------------------------------------------------------------------------------
// EEQ(2019) of dftd4 (WP12 stage 4). Claude Generated (Oct 2026). Port of multicharge 0.6.0 model/eeq.f90 (0d:
// get_xvec, get_amat_0d, get_damat_0d, get_xvec_derivs) with param.f90 new_eeq2019_model (cutoff 25 Bohr,
// cn_exp 7.5, cn_max 8) and mctc-lib's erf counting function; the charges are what dftd4 4.3.1 feeds into
// weight_references for DFT-D4 (disp.f90 get_dispersion: get_charges(disp%mchrg, ...)).
// ---------------------------------------------------------------------------------------------------------------
namespace {
constexpr double kEeqCnCutoff = 25.0;  // Bohr, multicharge new_eeq2019_model
constexpr double kEeqKcn = 7.5;        // steepness of the erf count
constexpr double kEeqCnMax = 8.0;      // log cut of the CN
constexpr double kEeqReg = 1.0e-14;    // multicharge reg in sqrt(cn + reg)
}

int D4ChargeModel::eeq2019MaxZ() { return d4eeq::kEeqElements; }

Vector D4ChargeModel::computeChargesEeq2019(const std::vector<int>& atoms, const Matrix& geom_bohr, double total_charge)
{
    const int N = static_cast<int>(atoms.size());
    m_n = N;
    m_atoms = atoms;
    m_geom = geom_bohr;
    if (N == 0) { m_q = Vector(); return m_q; }

    const D4CovalentCNTables tab = d4CovalentCNTables();
    std::vector<double> chi(N), eta(N), rad(N);
    m_alp.assign(N, 0.0);
    m_cnf.assign(N, 0.0);
    m_rcov_bohr.assign(N, 0.0);
    for (int i = 0; i < N; ++i) {
        const int Z = atoms[i];
        if (Z < 1 || Z > d4eeq::kEeqElements) {
            m_n = 0;
            m_q = Vector();
            return m_q;  // caller checks the size (no EEQ(2019) parameters beyond Lr)
        }
        chi[i] = d4eeq::kEeq2019Chi[Z - 1];
        eta[i] = d4eeq::kEeq2019Eta[Z - 1];
        rad[i] = d4eeq::kEeq2019Rad[Z - 1];
        m_alp[i] = rad[i] * rad[i];
        m_cnf[i] = d4eeq::kEeq2019Kcn[Z - 1];
        m_rcov_bohr[i] = tab.rad[Z];
    }

    // CN: erf count 0.5 (1 + erf(-kcn (r/rc - 1))), rc = rcov_i + rcov_j, then the log cut at cn_max.
    m_cn_raw.assign(N, 0.0);
    m_cn.assign(N, 0.0);
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < i; ++j) {
            const double r = (m_geom.row(i) - m_geom.row(j)).norm();
            if (r > kEeqCnCutoff) continue;
            const double rr = r / (m_rcov_bohr[i] + m_rcov_bohr[j]);
            const double c = 0.5 * (1.0 + std::erf(-kEeqKcn * (rr - 1.0)));
            m_cn_raw[i] += c;
            m_cn_raw[j] += c;
        }
    for (int i = 0; i < N; ++i)
        m_cn[i] = std::log(1.0 + std::exp(kEeqCnMax)) - std::log(1.0 + std::exp(kEeqCnMax - m_cn_raw[i]));

    // [A 1; 1 0] [q; lambda] = [x; Q]
    const int m = N + 1;
    Matrix M = Matrix::Zero(m, m);
    Vector c = Vector::Zero(m);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < i; ++j) {
            const double r = (m_geom.row(i) - m_geom.row(j)).norm();
            const double g = 1.0 / std::sqrt(m_alp[i] + m_alp[j]);
            const double aij = std::erf(g * r) / r;
            M(i, j) = aij;
            M(j, i) = aij;
        }
        M(i, i) = eta[i] + TSQRT2PI / rad[i];
        M(i, N) = 1.0;
        M(N, i) = 1.0;
        c(i) = -chi[i] + m_cnf[i] / std::sqrt(m_cn[i] + kEeqReg) * m_cn[i];
    }
    c(N) = total_charge;
    m_lu.compute(M);
    const Vector y = m_lu.solve(c);
    m_q = y.head(N);
    return m_q;
}

void D4ChargeModel::addChargeResponseGradientEeq2019(const Vector& dEdq, Matrix& grad_out) const
{
    const int N = m_n;
    if (N <= 1 || dEdq.size() != N || grad_out.rows() != N || grad_out.cols() != 3) return;

    // Adjoint: M z = [dE/dq; 0]; dE/dR = z^T (dx/dR - (dA/dR) q).
    Vector rhs = Vector::Zero(N + 1);
    rhs.head(N) = dEdq;
    const Vector z = m_lu.solve(rhs);

    // u_i = z_i dx_i/dCN'_i dCN'_i/dCN_i (multicharge get_xvec_derivs, mctc cut derivative)
    std::vector<double> u(N, 0.0);
    for (int i = 0; i < N; ++i) {
        const double cn = m_cn[i];
        const double dxdcn = m_cnf[i] * (0.5 * cn + kEeqReg) / std::pow(cn + kEeqReg, 1.5);
        const double dcut = std::exp(kEeqCnMax) / (std::exp(kEeqCnMax) + std::exp(m_cn_raw[i]));
        u[i] = z(i) * dxdcn * dcut;
    }
    for (int a = 0; a < N; ++a)
        for (int b = 0; b < a; ++b) {
            const Eigen::Vector3d rab = (m_geom.row(a) - m_geom.row(b)).transpose();
            const double r = rab.norm();
            if (r < 1.0e-10) continue;
            const Eigen::Vector3d uhat = rab / r;
            // Coulomb element: d/dr [erf(g r)/r]
            const double g = 1.0 / std::sqrt(m_alp[a] + m_alp[b]);
            const double dA = TWO_OVER_SQRTPI * g * std::exp(-g * g * r * r) / r - std::erf(g * r) / (r * r);
            // CN count: d/dr 0.5 (1 + erf(-kcn (r/rc - 1))) = -kcn/sqrt(pi) exp(-(kcn (r/rc - 1))^2) / rc
            double dcount = 0.0;
            if (r <= kEeqCnCutoff) {
                const double rc = m_rcov_bohr[a] + m_rcov_bohr[b];
                const double t = kEeqKcn * (r / rc - 1.0);
                dcount = -kEeqKcn * 0.5 * TWO_OVER_SQRTPI * std::exp(-t * t) / rc;
            }
            const double coeff = (u[a] + u[b]) * dcount - (z(a) * m_q(b) + z(b) * m_q(a)) * dA;
            grad_out.row(a) += (coeff * uhat).transpose();
            grad_out.row(b) -= (coeff * uhat).transpose();
        }
}

}  // namespace curcuma::dispersion
