/*
 * <xTB CPSCF / Z-vector response — orbital Hessian and PCG solver>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Implements the orbital-Hessian matrix-vector product and a
 * preconditioned conjugate-gradient Z-vector solver needed for the
 * Mulliken charge-response gradient in GFN2-D4 (Phase 3b).
 *
 * References:
 *   Handy & Schaefer, J. Chem. Phys. 81 (1984) 5031 — Z-vector method.
 *   Pulay, Mol. Phys. 17 (1969) 197 — response-theory gradient.
 *
 * Claude Generated (Phase 3b-2, May 2026). GPL-3.0.
 */

#include "xtb_native.h"
#include "STO_CGTO.hpp"
#include "xtb_coulomb.hpp"
#include "xtb_multipole_ints.hpp"

#include "parameters/gfn1_params.hpp"
#include "parameters/gfn2_params.hpp"
#include "parameters/xtb_params_extra.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>

namespace curcuma::xtb {

// Local CGTO-shell converter (same as xtb_gradient.cpp / xtb_h0.cpp)
static CGTO::Shell as_cgto_shell_r(const CGTOShell& cg)
{
    CGTO::Shell s;
    s.ang   = cg.ang;
    s.nprim = static_cast<int>(cg.alpha.size());
    s.alpha = cg.alpha;
    s.coeff = cg.coeff;
    return s;
}

static inline int ao_to_type_r(int ang, int local_ao)
{
    if (ang == 0) return 0;
    if (ang == 1) { static const int p_map[3] = {2, 3, 1}; return p_map[local_ao]; }
    return -1;
}

/* ------------------------------------------------------------------ *
 *  mpPairGrad()  —  bilinear multipole-interaction pair gradient     *
 *                                                                    *
 *  Evaluates G(M_i, M_j), the per-pair Cartesian gradient of the     *
 *  GFN2 SD/DD/SQ multipole interaction energy, for atom iat (force   *
 *  on iat; jat receives −g). Faithful copy of the energy gradient    *
 *  (xtb_gradient.cpp Section 5) factored so that EVERY bilinear term *
 *  draws exactly one moment factor from M_i and one from M_j, and    *
 *  the geometry/damping prefactors (g3/g5/g7, fdmp*, ddmp*, rr_dmp)  *
 *  are moment-independent.                                           *
 *                                                                    *
 *  Because G is bilinear in (M_i, M_j), its linearization in the     *
 *  moments is  G(δM_i, M_j) + G(M_i, δM_j)  — i.e. call this twice   *
 *  with one atom's moment triple swapped to its δ-counterpart. This  *
 *  is the Mulliken charge-response of the multipole interaction      *
 *  (Phase 3b-4b); the δ-moments come from mullikenFromDensity(D_z).  *
 *                                                                    *
 *  vec = R_jat − R_iat (TBLite convention). qp packing: xx,xy,yy,    *
 *  xz,yz,zz. fddr is the ∂E/∂mrad scalar (added to both atoms).      *
 *  mask selects term families (bit0 SD, bit1 DD, bit2 SQ) for        *
 *  incremental FD validation.                                        *
 * ------------------------------------------------------------------ */
static void mpPairGrad(
    const double vec[3], double r2,
    double g1, double g3, double g5, double g7,
    double fdmp3, double fdmp5, double ddmp3, double ddmp5,
    double rr_dmp, double mp_dmp3v, double mp_dmp5v,
    double qi, const double dpi[3], const double qpi[6],
    double qj, const double dpj[3], const double qpj[6],
    double g_out[3], double& fddr_out, unsigned mask)
{
    const double vx = vec[0], vy = vec[1], vz = vec[2];
    g_out[0] = g_out[1] = g_out[2] = 0.0;
    fddr_out = 0.0;

    // ── SD (charge-dipole, g3) ────────────────────────────────────────
    if (mask & 1u) {
        const double dpiqj = (vx*dpi[0] + vy*dpi[1] + vz*dpi[2]) * qj;
        const double qidpj = (vx*dpj[0] + vy*dpj[1] + vz*dpj[2]) * qi;
        const double diff_sd = dpiqj - qidpj;
        g_out[0] += -ddmp3*vx*diff_sd + fdmp3*g3*(qi*dpj[0] - qj*dpi[0]);
        g_out[1] += -ddmp3*vy*diff_sd + fdmp3*g3*(qi*dpj[1] - qj*dpi[1]);
        g_out[2] += -ddmp3*vz*diff_sd + fdmp3*g3*(qi*dpj[2] - qj*dpi[2]);
        fddr_out += 3.0*diff_sd*mp_dmp3v*fdmp3*g3*(fdmp3/rr_dmp)*std::pow(rr_dmp*g1, mp_dmp3v);
    }

    // ── DD (dipole-dipole, g5) ────────────────────────────────────────
    if (mask & 2u) {
        const double dpidpj = dpi[0]*dpj[0] + dpi[1]*dpj[1] + dpi[2]*dpj[2];
        const double dpiv = vx*dpi[0] + vy*dpi[1] + vz*dpi[2];
        const double dpjv = vx*dpj[0] + vy*dpj[1] + vz*dpj[2];
        const double edd = dpidpj*r2 - 3.0*dpjv*dpiv;
        g_out[0] += -2.0*fdmp5*g5*dpidpj*vx + 3.0*fdmp5*g5*(dpiv*dpj[0] + dpjv*dpi[0]) - edd*ddmp5*g7*vx;
        g_out[1] += -2.0*fdmp5*g5*dpidpj*vy + 3.0*fdmp5*g5*(dpiv*dpj[1] + dpjv*dpi[1]) - edd*ddmp5*g7*vy;
        g_out[2] += -2.0*fdmp5*g5*dpidpj*vz + 3.0*fdmp5*g5*(dpiv*dpj[2] + dpjv*dpi[2]) - edd*ddmp5*g7*vz;
        fddr_out += 3.0*edd*mp_dmp5v*fdmp5*g5*(fdmp5/rr_dmp)*std::pow(rr_dmp*g1, mp_dmp5v);
    }

    // ── SQ (charge-quadrupole, g5) ────────────────────────────────────
    if (mask & 4u) {
        const double eq =
              (qj*qpi[0] + qpj[0]*qi) * vx*vx
            + 2.0*(qj*qpi[1] + qpj[1]*qi) * vx*vy
            + (qj*qpi[2] + qpj[2]*qi) * vy*vy
            + 2.0*(qj*qpi[3] + qpj[3]*qi) * vx*vz
            + 2.0*(qj*qpi[4] + qpj[4]*qi) * vy*vz
            + (qj*qpi[5] + qpj[5]*qi) * vz*vz;
        const double tivx = qpi[0]*vx + qpi[1]*vy + qpi[3]*vz;
        const double tivy = qpi[1]*vx + qpi[2]*vy + qpi[4]*vz;
        const double tivz = qpi[3]*vx + qpi[4]*vy + qpi[5]*vz;
        const double tjvx = qpj[0]*vx + qpj[1]*vy + qpj[3]*vz;
        const double tjvy = qpj[1]*vx + qpj[2]*vy + qpj[4]*vz;
        const double tjvz = qpj[3]*vx + qpj[4]*vy + qpj[5]*vz;
        g_out[0] += -eq*ddmp5*g7*vx - 2.0*fdmp5*g5*(qi*tjvx + qj*tivx);
        g_out[1] += -eq*ddmp5*g7*vy - 2.0*fdmp5*g5*(qi*tjvy + qj*tivy);
        g_out[2] += -eq*ddmp5*g7*vz - 2.0*fdmp5*g5*(qi*tjvz + qj*tivz);
        fddr_out += eq*3.0*mp_dmp5v*fdmp5*g5*(fdmp5/rr_dmp)*std::pow(rr_dmp*g1, mp_dmp5v);
    }
}

/* ------------------------------------------------------------------ *
 *  applyOrbitalHessian()                                             *
 *                                                                    *
 *  Computes A·z = (ε_a − ε_i)·z + K·z  where K is the SCC coupling  *
 *  kernel ∂²E_SCC/∂P².                                              *
 *                                                                    *
 *  Input:  z_ov — nocc × nvirt occ-virt amplitude matrix.            *
 *  Output: A·z  — same shape.                                        *
 *                                                                    *
 *  Kernel K·z steps:                                                 *
 *    1. δP = 2·(C_occ·z·C_virt^T + C_virt·z^T·C_occ^T)              *
 *    2. δq_sh, δq_at, δdp_at, δqp_at from mullikenFromDensity(δP)   *
 *    3. δpot: Coulomb γ·δq_sh + third-order Jacobian (2·Γ·q_SCF)·δq *
 *             + multipole maps (GFN2) from δq_at/δdp_at/δqp_at      *
 *    4. δF = buildFockFromPotential(δpot)                            *
 *    5. (K·z)_ia = (C_occ^T · δF · C_virt)_ia                       *
 * ------------------------------------------------------------------ */
Matrix XTB::applyOrbitalHessian(const Matrix& z_ov) const
{
    const int nao  = m_basis.nao;
    const int nsh  = m_basis.nsh;
    const int nat  = m_atomcount;
    const int nocc = static_cast<int>(std::floor(m_wfn.nocc / 2.0));
    const int nvirt = nao - nocc;

    const Matrix C_occ  = m_wfn.C.leftCols(nocc);
    const Matrix C_virt = m_wfn.C.rightCols(nvirt);

    // --- 1. δP from occ-virt rotation ---
    Matrix dP = 2.0 * (C_occ * z_ov * C_virt.transpose()
              + C_virt * z_ov.transpose() * C_occ.transpose());

    // --- 2. Mulliken populations from δP ---
    Vector dq_sh, dq_at;
    Eigen::MatrixXd ddp_at, dqp_at;
    mullikenFromDensity(dP, dq_sh, dq_at, ddp_at, dqp_at);

    // --- 3. Build response potential δpot ---
    Potential dpot;
    dpot.v_sh = Vector::Zero(nsh);
    dpot.v_at = Vector::Zero(nat);
    // v_dp / v_qp left unset here; addMultipolePotential fills them when called.

    // Coulomb kernel: δv_sh += γ · δq_sh  (linear, exact)
    addCoulombShellPotential(dpot, dq_sh);

    // Third-order Jacobian: linearization of v = q² · Γ around q_SCF
    if (m_method == MethodType::GFN2) {
        // Shell-resolved: δv_sh(s) += 2·Γ_s·q_sh_SCF(s) · δq_sh(s)
        Vector kernel = thirdOrderKernelDiag(m_wfn.q_sh, m_wfn.q_at);
        dpot.v_sh += kernel.cwiseProduct(dq_sh);
    } else {
        // GFN1: atom-resolved δv_at only. buildFockFromPotential() sums
        // v_sh + v_at, so broadcasting into v_sh as well would double-count the
        // third-order response (mirrors the addThirdOrderPotential fix).
        for (int i = 0; i < nat; ++i) {
            const double dv = 2.0 * coulomb::atom_hubbard_deriv_gfn1(m_basis.z[i])
                            * m_wfn.q_at(i) * dq_at(i);
            dpot.v_at(i) += dv;
        }
    }

    // Multipole interaction kernel (GFN2 only): δv_at, δv_dp, δv_qp
    // The multipole maps are linear — addMultipolePotential applies them to δ-quantities.
    if (m_method == MethodType::GFN2 && m_mp_initialized) {
        addMultipolePotential(dpot, dq_at, ddp_at, dqp_at);
    }

    // --- 4. δF from δpot (no H0 term) ---
    Matrix dF = buildFockFromPotential(dpot);

    // --- 5. (K·z)_ia = (C_occ^T · δF · C_virt)_ia  [nocc × nvirt] ---
    Matrix Kz = C_occ.transpose() * dF * C_virt;

    // --- Add diagonal orbital-energy differences: (A·z)_ia = (ε_a − ε_i)·z_ia + (K·z)_ia ---
    Matrix Az = Kz;
    for (int i = 0; i < nocc; ++i)
        for (int a = 0; a < nvirt; ++a)
            Az(i, a) += (m_wfn.eps(nocc + a) - m_wfn.eps(i)) * z_ov(i, a);

    return Az;
}

/* ------------------------------------------------------------------ *
 *  solveZVector()                                                    *
 *                                                                    *
 *  Solves A·z = rhs_ov using preconditioned conjugate gradient (PCG) *
 *  in the nocc × nvirt occ-virt space.                               *
 *                                                                    *
 *  Preconditioner: P_ia = 1/(ε_a − ε_i)   (always positive for a    *
 *  stable SCF solution with ε_a > ε_i).                              *
 *  Tolerance: 1e-10 on residual norm; max 200 iterations.            *
 * ------------------------------------------------------------------ */
Matrix XTB::solveZVector(const Matrix& rhs_ov) const
{
    const int nao  = m_basis.nao;
    const int nocc = static_cast<int>(std::floor(m_wfn.nocc / 2.0));
    const int nvirt = nao - nocc;

    // Preconditioner diagonal 1/(ε_a − ε_i)
    Matrix prec(nocc, nvirt);
    for (int i = 0; i < nocc; ++i)
        for (int a = 0; a < nvirt; ++a)
            prec(i, a) = 1.0 / std::max(m_wfn.eps(nocc + a) - m_wfn.eps(i), 1.0e-6);

    // Initial guess: z = P * rhs (first-order PT, ignores K)
    Matrix z = prec.cwiseProduct(rhs_ov);
    Matrix r = rhs_ov - applyOrbitalHessian(z);
    Matrix y = prec.cwiseProduct(r);
    Matrix p = y;
    double ry = r.cwiseProduct(y).sum();

    const int    max_iter  = 200;
    const double tol       = 1.0e-10;
    const double rhs_norm  = rhs_ov.norm();

    for (int iter = 0; iter < max_iter; ++iter) {
        if (r.norm() <= tol * (rhs_norm + 1.0e-30)) break;

        Matrix Ap = applyOrbitalHessian(p);
        const double pAp = p.cwiseProduct(Ap).sum();
        if (std::fabs(pAp) < 1.0e-30) break;

        const double alpha = ry / pAp;
        z  += alpha * p;
        r  -= alpha * Ap;
        y   = prec.cwiseProduct(r);

        const double ry_new = r.cwiseProduct(y).sum();
        const double beta   = ry_new / ry;
        ry  = ry_new;
        p   = y + beta * p;
    }
    return z;
}

/* ------------------------------------------------------------------ *
 *  computeMullikenChargeResponse()                                   *
 *                                                                    *
 *  Folds the D4 Mulliken charge-response  Σ_A w_A·∂q_A/∂x  into      *
 *  grad_out (Eh/Bohr), with w_A = dEdq(A) = ∂E_D4/∂q_A.             *
 *                                                                    *
 *  L = Σ_A w_A q_A = const − Tr(Λw·P·S),  Λw = diag(w_{atom(μ)}).   *
 *  dL/dx = Tr(D_z·Fˣ_skel) − Tr(W_z·Sˣ)    [Z-vector response]       *
 *          + linearized SCC charge products (Coulomb δq·q + q·δq).   *
 *                                                                    *
 *  The response is the LINEARIZATION of the energy gradient:         *
 *  P→D_z (Pulay sval), W→W_z, and the Coulomb charge product is      *
 *  linearized (q_is·q_js → δq_is·q_js + q_is·δq_js with δq from D_z).*
 *                                                                    *
 *  RHS_SIGN=+1 is FD-validated against the raw ∂q_A/∂x (test Test D).*
 *  Multipole charge-response terms exist but are gated OFF            *
 *  (MP_RESPONSE_ENABLED): correct for low-multipole systems, but     *
 *  mis-contracted for strongly polar molecules — see                 *
 *  docs/PHASE3B5_MULTIPOLE_RESPONSE_WP.md. The D4 gradient target     *
 *  (<5e-5) is met without them.                                      *
 * ------------------------------------------------------------------ */
void XTB::computeMullikenChargeResponse(const Vector& dEdq, Matrix& grad_out) const
{
    using namespace curcuma::xtb::gfn1_params;
    using namespace curcuma::xtb::gfn2_params;

    const int nat   = m_atomcount;
    const int nsh   = m_basis.nsh;
    const int nao   = m_basis.nao;
    const int nocc  = static_cast<int>(std::round(m_wfn.nocc / 2.0));
    const int nvirt = nao - nocc;
    if (nocc <= 0 || nvirt <= 0) return;

    // ---- response sign and multipole gating ----
    // RHS_SIGN=+1: sign of the property-gradient RHS. Validated against FD on the
    // RAW charge response ∂q_A/∂x (decoupled from D4) for H2O/HCN/CH4 — see
    // test_xtb_cpscf Test D, which calls this routine with dEdq=e_A. The previous
    // −1 was an artifact of a sign-blind max-abs validation metric and inverted
    // the entire response (anal ≈ −num); fixing it cut the error 5–8×.
    constexpr double RHS_SIGN = 1.0;

    // Multipole charge-response terms (3b-4a integral-Pulay + 7b interaction).
    // Correct and necessary for low-multipole systems (CH4 → ~1% vs FD) but
    // over-/mis-contracted for strongly multipolar molecules (H2O/HCN worsen).
    // Disabled pending a corrected property-gradient derivation; the D4 gradient
    // target (<5e-5) is met without them. See docs/PHASE3B5_MULTIPOLE_RESPONSE_WP.md.
    constexpr bool MP_RESPONSE_ENABLED = false;
    constexpr unsigned MP_RESP_MASK = 0x7u;  // bit0 SD, bit1 DD, bit2 SQ (when enabled)

    // ---- geometry in Bohr ----
    std::vector<double> xyz(3 * nat);
    for (int i = 0; i < nat; ++i) {
        xyz[3*i+0] = m_geometry(i, 0) * AA_TO_AU;
        xyz[3*i+1] = m_geometry(i, 1) * AA_TO_AU;
        xyz[3*i+2] = m_geometry(i, 2) * AA_TO_AU;
    }

    const Matrix C_occ  = m_wfn.C.leftCols(nocc);
    const Matrix C_virt = m_wfn.C.rightCols(nvirt);

    // ---- 1. property gradient G_P = ½ S_μν (Λw_μ + Λw_ν) ----
    Vector lam_ao(nao);
    for (int mu = 0; mu < nao; ++mu) lam_ao(mu) = dEdq(m_basis.ao2at[mu]);
    Matrix G_P(nao, nao);
    for (int mu = 0; mu < nao; ++mu)
        for (int nu = 0; nu < nao; ++nu)
            G_P(mu, nu) = 0.5 * m_S(mu, nu) * (lam_ao(mu) + lam_ao(nu));

    // ---- 2. Z-vector solve  A z = RHS_SIGN·(C_occ^T G_P C_virt) ----
    const Matrix rhs = RHS_SIGN * (C_occ.transpose() * G_P * C_virt);
    const Matrix z   = solveZVector(rhs);

    // ---- 3. relaxed density D_z and energy-weighted W_z (closed-shell ×2) ----
    const Matrix D_z = 2.0 * (C_occ * z * C_virt.transpose()
                     + C_virt * z.transpose() * C_occ.transpose());
    const Matrix Ez  = m_wfn.eps.head(nocc).asDiagonal() * z;   // ε_i · z_ia
    const Matrix W_z = 2.0 * (C_occ * Ez * C_virt.transpose()
                     + C_virt * Ez.transpose() * C_occ.transpose());

    // ---- 4. relaxed shell/atom charges from D_z (Mulliken) ----
    Vector dq_sh, dq_at;
    Eigen::MatrixXd ddp, dqp;
    mullikenFromDensity(D_z, dq_sh, dq_at, ddp, dqp);

    // ---- self-energies and valence flags (for h_av) ----
    // CN computed locally: m_coordination_numbers is only assigned AFTER
    // calcDispersionEnergy() in Calculation(), so it is not yet valid here.
    const Vector cn_resp = computeCoordinationNumbers();
    Vector se;
    getSelfEnergies(cn_resp, se);
    std::vector<bool> valence(nsh, false);
    if (m_method == MethodType::GFN1) {
        for (int iat = 0; iat < nat; ++iat) {
            bool ang_seen[3] = {false, false, false};
            for (int ish = 0; ish < m_basis.nsh_at[iat]; ++ish) {
                const int sh = m_basis.ish_at[iat] + ish;
                const int l  = m_basis.ang_sh[sh];
                if (!ang_seen[l]) { valence[sh] = true; ang_seen[l] = true; }
            }
        }
    }

    Vector dEdcn = Vector::Zero(nat);

    // ---- 5. H0-diagonal CN coupling from D_z (mirrors energy grad) ----
    for (int ish = 0; ish < nsh; ++ish) {
        const int iat = m_basis.sh2at[ish];
        const double kcn_s = m_h0.kcn[ish];
        for (int mu = m_basis.iao_sh[ish]; mu < m_basis.iao_sh[ish] + m_basis.nao_sh[ish]; ++mu)
            dEdcn(iat) += (-kcn_s) * D_z(mu, mu);
    }

    // ---- 6. Pulay/H0 + explicit-overlap gradient (off-site shell pairs) ----
    for (int iat = 0; iat < nat; ++iat) {
        for (int jat = iat + 1; jat < nat; ++jat) {
            const double dx_ij = xyz[3*iat+0] - xyz[3*jat+0];
            const double dy_ij = xyz[3*iat+1] - xyz[3*jat+1];
            const double dz_ij = xyz[3*iat+2] - xyz[3*jat+2];
            const double r2    = dx_ij*dx_ij + dy_ij*dy_ij + dz_ij*dz_ij;
            if (r2 < 1.0e-12) continue;
            const double r = std::sqrt(r2);

            const int zi = m_basis.z[iat];
            const int zj = m_basis.z[jat];
            const double rad_sum = atomic_rad_au(zi) + atomic_rad_au(zj);
            const double rr = std::sqrt(r / rad_sum);

            for (int ia = 0; ia < m_basis.nsh_at[iat]; ++ia) {
                const int ish_a    = m_basis.ish_at[iat] + ia;
                const int ia_start = m_basis.iao_sh[ish_a];
                const int ia_nao   = m_basis.nao_sh[ish_a];
                const double pi_a  = 1.0 + m_h0.shpoly[ish_a] * rr;
                const double zeta_a = m_basis.cgto[ish_a].slater_exp;
                const CGTO::Shell sh_a = as_cgto_shell_r(m_basis.cgto[ish_a]);

                for (int ib = 0; ib < m_basis.nsh_at[jat]; ++ib) {
                    const int ish_b    = m_basis.ish_at[jat] + ib;
                    const int jb_start = m_basis.iao_sh[ish_b];
                    const int jb_nao   = m_basis.nao_sh[ish_b];
                    const double pi_b  = 1.0 + m_h0.shpoly[ish_b] * rr;
                    const double zeta_b = m_basis.cgto[ish_b].slater_exp;
                    const CGTO::Shell sh_b = as_cgto_shell_r(m_basis.cgto[ish_b]);

                    double hs;
                    if (m_method == MethodType::GFN1) {
                        const bool vi = valence[ish_a], vj = valence[ish_b];
                        const double den = std::pow(pauling_en[zi-1] - pauling_en[zj-1], 2);
                        if (vi && vj)
                            hs = gfn1::kpair(zi, zj)
                               * gfn1::kshell(m_basis.ang_sh[ish_a], m_basis.ang_sh[ish_b])
                               * (1.0 + gfn1::enscale * den);
                        else if (vi && !vj)
                            hs = 0.5 * (gfn1::kshell(m_basis.ang_sh[ish_a], m_basis.ang_sh[ish_a]) + gfn1::kdiff);
                        else if (!vi && vj)
                            hs = 0.5 * (gfn1::kshell(m_basis.ang_sh[ish_b], m_basis.ang_sh[ish_b]) + gfn1::kdiff);
                        else
                            hs = gfn1::kdiff;
                    } else {
                        const double den = std::pow(pauling_en[zi-1] - pauling_en[zj-1], 2);
                        const double enp = 1.0 + gfn2::enscale * den;
                        const double km  = gfn2::kpair(zi, zj)
                                         * gfn2::kshell(m_basis.ang_sh[ish_a], m_basis.ang_sh[ish_b])
                                         * enp;
                        const double zij = std::pow(2.0 * std::sqrt(zeta_a * zeta_b)
                                                     / (zeta_a + zeta_b), gfn2::wexp);
                        hs = zij * km;
                    }
                    const double h_factor = hs * (pi_a * pi_b);
                    const double avg_eps  = 0.5 * (se(ish_a) + se(ish_b));
                    const double h_av     = avg_eps * h_factor;
                    const double dlog_pi_dr_r = (m_h0.shpoly[ish_a] / pi_a
                                                + m_h0.shpoly[ish_b] / pi_b) * rr / (2.0 * r2);
                    const double v_a = m_pot.v_sh(ish_a) + m_pot.v_at(iat);
                    const double v_b = m_pot.v_sh(ish_b) + m_pot.v_at(jat);

                    double G_sval[3]       = {0.0, 0.0, 0.0};
                    double G_shpoly_scalar = 0.0;
                    double cn_sum          = 0.0;

                    for (int mu = ia_start; mu < ia_start + ia_nao; ++mu) {
                        const int t_a = ao_to_type_r(m_basis.ang_sh[ish_a], mu - ia_start);
                        for (int nu = jb_start; nu < jb_start + jb_nao; ++nu) {
                            const int t_b = ao_to_type_r(m_basis.ang_sh[ish_b], nu - jb_start);
                            if (t_a < 0 || t_b < 0) continue;

                            const double Dmn = D_z(mu, nu);
                            const double Wmn = W_z(mu, nu);
                            const double Smn = m_S(mu, nu);
                            const double H0mn = m_H0(mu, nu);

                            double dS[3];
                            CGTO::cgto_overlap_grad(sh_a, sh_b,
                                xyz[3*iat+0], xyz[3*iat+1], xyz[3*iat+2],
                                xyz[3*jat+0], xyz[3*jat+1], xyz[3*jat+2],
                                t_a, t_b, dS);

                            // response sval = linearization of energy sval (P→D_z, W→W_z)
                            const double sval = 2.0*Dmn*h_av - 2.0*Wmn - Dmn*(v_a + v_b);
                            G_sval[0] += sval * dS[0];
                            G_sval[1] += sval * dS[1];
                            G_sval[2] += sval * dS[2];

                            // Multipole integral Pulay (AP5b) with D_z (converged v_dp/v_qp).
                            // Same contraction as the energy gradient, P→D_z.
                            if (m_method == MethodType::GFN2 && m_mp_initialized && MP_RESPONSE_ENABLED) {
                                double S_mp, D_mp[3], Q_mp[6];
                                double dD_dA[3][3], dD_dB[3][3];
                                double dQ_dA[3][6], dQ_dB[3][6];
                                using namespace curcuma::xtb::multipole_ints;
                                cgto_multipole_grad_transformed(
                                    sh_a, sh_b,
                                    xyz[3*iat+0], xyz[3*iat+1], xyz[3*iat+2],
                                    xyz[3*jat+0], xyz[3*jat+1], xyz[3*jat+2],
                                    t_a, t_b, S_mp, D_mp, Q_mp,
                                    dD_dA, dD_dB, dQ_dA, dQ_dB);
                                const double dR[3] = {
                                    xyz[3*jat+0]-xyz[3*iat+0],
                                    xyz[3*jat+1]-xyz[3*iat+1],
                                    xyz[3*jat+2]-xyz[3*iat+2]};
                                static const int qa6[6] = {0,0,1,0,1,2};
                                static const int qb6[6] = {0,1,1,2,2,2};
                                for (int l = 0; l < 3; ++l) {
                                    double term = 0.0;
                                    for (int k = 0; k < 3; ++k) {
                                        term += dD_dA[l][k] * m_pot.v_dp(k, jat);
                                        const double dDiat = dD_dA[l][k] + dR[k]*dS[l] - (k==l ? Smn : 0.0);
                                        term += dDiat * m_pot.v_dp(k, iat);
                                    }
                                    double Dqraw[6];
                                    for (int q = 0; q < 6; ++q) {
                                        const int a = qa6[q], b = qb6[q];
                                        Dqraw[q] = -(b==l ? D_mp[a] : 0.0) + dR[b]*dD_dA[l][a]
                                                 - (a==l ? D_mp[b] : 0.0) + dR[a]*dD_dA[l][b]
                                                 + (-(a==l ? dR[b] : 0.0) - (b==l ? dR[a] : 0.0))*Smn
                                                 + dR[a]*dR[b]*dS[l];
                                    }
                                    const double dtr_c = 0.5*(Dqraw[0] + Dqraw[2] + Dqraw[5]);
                                    for (int q = 0; q < 6; ++q) {
                                        const bool is_diag = (qa6[q] == qb6[q]);
                                        term += dQ_dA[l][q] * m_pot.v_qp(q, jat);
                                        const double dQiat = dQ_dA[l][q] + 1.5*Dqraw[q]
                                                           - (is_diag ? dtr_c : 0.0);
                                        term += dQiat * m_pot.v_qp(q, iat);
                                    }
                                    G_sval[l] -= Dmn * term;
                                }
                            }

                            G_shpoly_scalar += 2.0 * Dmn * H0mn * dlog_pi_dr_r;
                            cn_sum += Dmn * Smn;
                        }
                    }

                    const double Gx = G_sval[0] + G_shpoly_scalar * dx_ij;
                    const double Gy = G_sval[1] + G_shpoly_scalar * dy_ij;
                    const double Gz = G_sval[2] + G_shpoly_scalar * dz_ij;
                    grad_out(iat, 0) += Gx;  grad_out(jat, 0) -= Gx;
                    grad_out(iat, 1) += Gy;  grad_out(jat, 1) -= Gy;
                    grad_out(iat, 2) += Gz;  grad_out(jat, 2) -= Gz;

                    dEdcn(iat) += (-m_h0.kcn[ish_a]) * h_factor * cn_sum;
                    dEdcn(jat) += (-m_h0.kcn[ish_b]) * h_factor * cn_sum;
                }
            }
        }
    }

    // ---- 7. Coulomb (ES2) response: linearized charge product δq·q + q·δq ----
    {
        const double gexp = (m_method == MethodType::GFN1) ? gfn1_params::gexp : gfn2_params::gexp;
        const curcuma::xtb::coulomb::Method cm = (m_method == MethodType::GFN1)
            ? curcuma::xtb::coulomb::Method::GFN1 : curcuma::xtb::coulomb::Method::GFN2;
        std::vector<double> g_sh(nsh);
        for (int is = 0; is < nsh; ++is)
            g_sh[is] = curcuma::xtb::coulomb::shell_hardness(cm, m_basis.z[m_basis.sh2at[is]],
                                                             m_basis.ang_sh[is]);
        for (int is = 0; is < nsh; ++is) {
            const int iat = m_basis.sh2at[is];
            for (int js = 0; js < is; ++js) {
                const int jat = m_basis.sh2at[js];
                if (iat == jat) continue;
                const double dx = xyz[3*iat+0] - xyz[3*jat+0];
                const double dy = xyz[3*iat+1] - xyz[3*jat+1];
                const double dz = xyz[3*iat+2] - xyz[3*jat+2];
                const double r2 = dx*dx + dy*dy + dz*dz;
                if (r2 < 1.0e-12) continue;
                const double r1 = std::sqrt(r2);
                const double gam_bar = (m_method == MethodType::GFN1)
                    ? curcuma::xtb::coulomb::harmonic_average(g_sh[is], g_sh[js])
                    : curcuma::xtb::coulomb::arithmetic_average(g_sh[is], g_sh[js]);
                const double r1g     = std::pow(r1, gexp);
                const double gam_pow = std::pow(gam_bar, -gexp);
                const double gamma   = std::pow(r1g + gam_pow, -1.0 / gexp);
                const double dgamma_dr = -std::pow(r1, gexp - 2.0) * std::pow(gamma, gexp + 1.0);
                // linearized: δ(q_is·q_js) = δq_is·q_js + q_is·δq_js
                const double force = (dq_sh(is) * m_wfn.q_sh(js)
                                    + m_wfn.q_sh(is) * dq_sh(js)) * dgamma_dr;
                grad_out(iat, 0) += force * dx;  grad_out(jat, 0) -= force * dx;
                grad_out(iat, 1) += force * dy;  grad_out(jat, 1) -= force * dy;
                grad_out(iat, 2) += force * dz;  grad_out(jat, 2) -= force * dz;
            }
        }
    }

    // ---- 7b. Multipole-interaction response (Phase 3b-4b, GFN2 only) ----
    //  Linearization of the SD/DD/SQ interaction gradient (xtb_gradient.cpp
    //  Section 5) in the moments: G(δM_i,M_j)+G(M_i,δM_j), with δ-moments
    //  (dq_at, ddp, dqp) from D_z. The mrad/CN chain rule (fddr) is folded
    //  into the existing dEdcn accumulator (distributed in section 8).
    if (m_method == MethodType::GFN2 && m_mp_initialized && MP_RESPONSE_ENABLED) {
        Vector dEdr_mp = Vector::Zero(nat);
        for (int iat = 0; iat < nat; ++iat) {
            for (int jat = 0; jat < iat; ++jat) {
                const double vec[3] = {
                    xyz[3*jat+0] - xyz[3*iat+0],
                    xyz[3*jat+1] - xyz[3*iat+1],
                    xyz[3*jat+2] - xyz[3*iat+2]};
                const double r2 = vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2];
                if (r2 < 1.0e-12) continue;
                const double r1 = std::sqrt(r2);
                const double g1 = 1.0/r1, g3 = g1*g1*g1, g5 = g3*g1*g1, g7 = g5*g1*g1;
                const double rr_dmp = 0.5*(m_mp_mrad[iat] + m_mp_mrad[jat]);
                const double fdmp3 = 1.0/(1.0 + 6.0*std::pow(rr_dmp*g1, mp_dmp3));
                const double fdmp5 = 1.0/(1.0 + 6.0*std::pow(rr_dmp*g1, mp_dmp5));
                const double ddmp3 = g5*fdmp3*(-3.0 + mp_dmp3*(1.0 - fdmp3));
                const double ddmp5 = fdmp5*(-5.0 - mp_dmp5*(fdmp5 - 1.0));

                // SCF moments M
                const double qi = m_wfn.q_at(iat), qj = m_wfn.q_at(jat);
                const double dpi[3] = {m_wfn.dp_at(0,iat), m_wfn.dp_at(1,iat), m_wfn.dp_at(2,iat)};
                const double dpj[3] = {m_wfn.dp_at(0,jat), m_wfn.dp_at(1,jat), m_wfn.dp_at(2,jat)};
                double qpi[6], qpj[6];
                for (int k = 0; k < 6; ++k) { qpi[k] = m_wfn.qp_at(k,iat); qpj[k] = m_wfn.qp_at(k,jat); }

                // δ-moments from D_z (same convention as M)
                const double dqi = dq_at(iat), dqj = dq_at(jat);
                const double ddpi[3] = {ddp(0,iat), ddp(1,iat), ddp(2,iat)};
                const double ddpj[3] = {ddp(0,jat), ddp(1,jat), ddp(2,jat)};
                double dqpi[6], dqpj[6];
                for (int k = 0; k < 6; ++k) { dqpi[k] = dqp(k,iat); dqpj[k] = dqp(k,jat); }

                double gp1[3], gp2[3], f1 = 0.0, f2 = 0.0;
                // Pass 1: G(δM_i, M_j)
                mpPairGrad(vec, r2, g1, g3, g5, g7, fdmp3, fdmp5, ddmp3, ddmp5,
                           rr_dmp, mp_dmp3, mp_dmp5,
                           dqi, ddpi, dqpi, qj, dpj, qpj, gp1, f1, MP_RESP_MASK);
                // Pass 2: G(M_i, δM_j)
                mpPairGrad(vec, r2, g1, g3, g5, g7, fdmp3, fdmp5, ddmp3, ddmp5,
                           rr_dmp, mp_dmp3, mp_dmp5,
                           qi, dpi, qpi, dqj, ddpj, dqpj, gp2, f2, MP_RESP_MASK);

                const double gx = gp1[0] + gp2[0];
                const double gy = gp1[1] + gp2[1];
                const double gz = gp1[2] + gp2[2];
                grad_out(iat, 0) += gx;  grad_out(jat, 0) -= gx;
                grad_out(iat, 1) += gy;  grad_out(jat, 1) -= gy;
                grad_out(iat, 2) += gz;  grad_out(jat, 2) -= gz;

                const double fddr = f1 + f2;
                dEdr_mp(iat) += fddr;
                dEdr_mp(jat) += fddr;
            }
        }
        // mrad/CN chain rule (same dmrdcn as xtb_gradient.cpp Section 5)
        for (int i = 0; i < nat; ++i) {
            const int zi = m_atoms[i];
            const double rad = p_rad[zi - 1];
            const double vcn = p_vcn[zi - 1];
            const double arg = cn_resp(i) - vcn - mp_shift;
            const double t1  = std::exp(-mp_kexp * arg);
            const double t2  = (mp_rmax - rad) / (1.0 + t1);
            const double dmrdcn = -t2 * mp_kexp * t1 / (1.0 + t1);
            dEdcn(i) += dEdr_mp(i) * dmrdcn;
        }
    }

    // ---- 8. CN chain-rule distribution of dEdcn (mirrors energy grad section 4) ----
    {
        std::vector<double> rcov(nat);
        for (int i = 0; i < nat; ++i) rcov[i] = covalent_rad_d3_au(m_atoms[i]);
        for (int i = 0; i < nat; ++i) {
            for (int j = 0; j < i; ++j) {
                const double dx = xyz[3*i+0] - xyz[3*j+0];
                const double dy = xyz[3*i+1] - xyz[3*j+1];
                const double dz = xyz[3*i+2] - xyz[3*j+2];
                const double r2 = dx*dx + dy*dy + dz*dz;
                if (r2 < 1.0e-12) continue;
                const double r = std::sqrt(r2);
                const double rc = rcov[i] + rcov[j];
                double dcndr;
                if (m_method == MethodType::GFN1) {
                    constexpr double k = 16.0;
                    const double ef = std::exp(-k * (rc / r - 1.0));
                    const double c  = 1.0 / (1.0 + ef);
                    dcndr = -c * (1.0 - c) * k * rc / r2;
                } else {
                    constexpr double ka = 10.0, kb = 20.0, rs = 2.0;
                    const double c1 = 1.0 / (1.0 + std::exp(-ka * (rc / r - 1.0)));
                    const double c2 = 1.0 / (1.0 + std::exp(-kb * ((rc + rs) / r - 1.0)));
                    const double dc1dr = -c1 * (1.0 - c1) * ka * rc       / r2;
                    const double dc2dr = -c2 * (1.0 - c2) * kb * (rc + rs) / r2;
                    dcndr = dc1dr * c2 + c1 * dc2dr;
                }
                const double factor = (dEdcn(i) + dEdcn(j)) * dcndr / r;
                grad_out(i, 0) += factor * dx;  grad_out(j, 0) -= factor * dx;
                grad_out(i, 1) += factor * dy;  grad_out(j, 1) -= factor * dy;
                grad_out(i, 2) += factor * dz;  grad_out(j, 2) -= factor * dz;
            }
        }
    }
}

} // namespace curcuma::xtb
