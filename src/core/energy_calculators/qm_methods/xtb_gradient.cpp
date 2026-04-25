/*
 * <XTB analytical gradients — GFN1/GFN2>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Analytical gradients for GFN1/GFN2 tight-binding energy.
 * Fills m_gradient in Eh/Bohr; caller converts to Eh/Å via *= au.
 *
 * References:
 *   - external/tblite/src/tblite/xtb/h0.f90        (get_hamiltonian_gradient)
 *   - external/tblite/src/tblite/repulsion/effective.f90 (get_repulsion_derivs)
 *   - external/tblite/src/tblite/coulomb/charge/type.f90
 *
 * AP4: repulsion + H0/Pulay (s/p only) + Coulomb + CN chain-rule.
 * AP5 (TODO): GFN2 multipole gradient + overall energy accuracy.
 *
 * Claude Generated. GPL-3.0.
 */

#include "xtb_native.h"
#include "STO_CGTO.hpp"
#include "xtb_coulomb.hpp"

#include "parameters/gfn1_params.hpp"
#include "parameters/gfn2_params.hpp"
#include "parameters/xtb_params_extra.hpp"

#include "src/core/curcuma_logger.h"

#include <cmath>

namespace curcuma::xtb {

/* ---- local helpers shared with xtb_h0.cpp -------------------------------- */

// Convert internal CGTOShell → CGTO::Shell (duplicate of xtb_h0.cpp)
static CGTO::Shell as_cgto_shell_g(const CGTOShell& cg)
{
    CGTO::Shell s;
    s.ang   = cg.ang;
    s.nprim = static_cast<int>(cg.alpha.size());
    s.alpha = cg.alpha;
    s.coeff = cg.coeff;
    return s;
}

// AO-type within a shell (matching xtb_h0.cpp convention)
static inline int ao_to_type_g(int ang, int local_ao)
{
    if (ang == 0) return 0;
    if (ang == 1) {
        static const int p_map[3] = {2, 3, 1};  // py, pz, px (tblite ordering)
        return p_map[local_ao];
    }
    return -1;  // d not handled
}

/* ========================================================================== *
 *  XTB::calculateGradient()                                                   *
 *                                                                             *
 *  Fills m_gradient (nat×3) in Eh/Bohr.  Call after converged SCF.           *
 *                                                                             *
 *  Components:                                                                *
 *   1.  Repulsion gradient                                                    *
 *   2.  Hamiltonian/Pulay gradient (H0 + W + potential terms)                *
 *   3.  Coulomb (ES2) gradient                                                *
 *   4.  CN chain-rule gradient                                                *
 *  (GFN2 multipole gradient: deferred to Phase 4b)                           *
 * ========================================================================== */
void XTB::calculateGradient()
{
    const int nat = m_atomcount;
    const int nsh = m_basis.nsh;
    const int nao = m_basis.nao;

    m_gradient = Matrix::Zero(nat, 3);

    // ── Geometry in Bohr ────────────────────────────────────────────────────
    std::vector<double> xyz(3 * nat);
    for (int i = 0; i < nat; ++i) {
        xyz[3*i+0] = m_geometry(i, 0) * AA_TO_AU;
        xyz[3*i+1] = m_geometry(i, 1) * AA_TO_AU;
        xyz[3*i+2] = m_geometry(i, 2) * AA_TO_AU;
    }

    // ── Energy-weighted density matrix: W_μν = 2·Σ_{i<nocc} ε_i C_μi C_νi ─
    // Closed-shell: factor 2 for spin degeneracy; nocc counts electrons not pairs.
    const int nocc_orbs = static_cast<int>(std::round(m_wfn.nocc / 2.0));
    Matrix W = Matrix::Zero(nao, nao);
    for (int i = 0; i < nocc_orbs; ++i)
        W += 2.0 * m_wfn.eps(i) * m_wfn.C.col(i) * m_wfn.C.col(i).transpose();

    // ── v_ao: expand shell+atom potential to AO resolution ──────────────────
    // F_μν = H0_μν - 0.5·S_μν·(v_ao(μ)+v_ao(ν))  (xtb_scf.cpp:expand_potential)
    Vector v_ao = Vector::Zero(nao);
    for (int ao = 0; ao < nao; ++ao)
        v_ao(ao) = m_pot.v_sh(m_basis.ao2sh[ao]) + m_pot.v_at(m_basis.ao2at[ao]);

    // ── CN-shifted self-energies (needed for sval) ───────────────────────────
    Vector se;
    getSelfEnergies(m_coordination_numbers, se);

    // ── Valence flags (GFN1 only) ────────────────────────────────────────────
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

    // ── Element-Z per shell (shorthand) ─────────────────────────────────────
    std::vector<int> sh_z(nsh);
    for (int s = 0; s < nsh; ++s) sh_z[s] = m_basis.z[m_basis.sh2at[s]];

    // ==========================================================================
    //  1.  Repulsion gradient
    //      E_rep_ij = Z_i·Z_j/r · exp(-(α_i·α_j)^kexp · r^rexp)
    //      dE/dr    = -(rexp/r + a_pair·kexp·r^(kexp-1)) · E_pair   [rexp=1 here]
    //               = -(1/r + a_pair) · E_pair
    // ==========================================================================
    for (int i = 0; i < nat; ++i) {
        const int zi = m_atoms[i];
        double alfi, zeffi, kexp, rexp;
        if (m_method == MethodType::GFN1) {
            alfi  = gfn1_params::rep_alpha[zi - 1];
            zeffi = gfn1_params::rep_zeff [zi - 1];
            kexp  = gfn1_params::rep_kexp;
            rexp  = gfn1_params::rep_rexp;
        } else {
            alfi  = gfn2_params::rep_alpha[zi - 1];
            zeffi = gfn2_params::rep_zeff [zi - 1];
            kexp  = gfn2_params::rep_kexp;
            rexp  = gfn2_params::rep_rexp;
        }
        for (int j = 0; j < i; ++j) {
            const int zj = m_atoms[j];
            double alfj, zeffj;
            if (m_method == MethodType::GFN1) {
                alfj  = gfn1_params::rep_alpha[zj - 1];
                zeffj = gfn1_params::rep_zeff [zj - 1];
            } else {
                alfj  = gfn2_params::rep_alpha[zj - 1];
                zeffj = gfn2_params::rep_zeff [zj - 1];
            }
            const double dx = xyz[3*i+0] - xyz[3*j+0];
            const double dy = xyz[3*i+1] - xyz[3*j+1];
            const double dz = xyz[3*i+2] - xyz[3*j+2];
            const double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 < 1.0e-12) continue;
            const double r  = std::sqrt(r2);

            // a_pair = (α_i·α_j)^kexp  (constant for this atom pair)
            const double a_pair  = std::pow(alfi * alfj, kexp);
            const double E_pair  = zeffi * zeffj / r * std::exp(-a_pair * std::pow(r, rexp));
            // dE/dr = -(rexp/r + a_pair·rexp·r^(rexp-1)) · E_pair
            const double dEdr    = -(rexp / r + a_pair * rexp * std::pow(r, rexp - 1.0)) * E_pair;

            const double fx = dEdr * dx / r;
            const double fy = dEdr * dy / r;
            const double fz = dEdr * dz / r;
            m_gradient(i, 0) += fx;  m_gradient(j, 0) -= fx;
            m_gradient(i, 1) += fy;  m_gradient(j, 1) -= fy;
            m_gradient(i, 2) += fz;  m_gradient(j, 2) -= fz;
        }
    }

    // ==========================================================================
    //  2a.  On-site CN contribution: dEdcn[iat] += (-kcn[ish]) · P(μ,μ)
    //       (diagonal H0 elements: H0_μμ = se[ish]; dH0/dCN = -kcn[ish])
    // ==========================================================================
    Vector dEdcn = Vector::Zero(nat);

    for (int ish = 0; ish < nsh; ++ish) {
        const int iat      = m_basis.sh2at[ish];
        const int ia_start = m_basis.iao_sh[ish];
        const int ia_nao   = m_basis.nao_sh[ish];
        const double kcn_s = m_h0.kcn[ish];
        for (int ia = 0; ia < ia_nao; ++ia) {
            const int mu = ia_start + ia;
            dEdcn(iat) += (-kcn_s) * m_wfn.P(mu, mu);
        }
    }

    // ==========================================================================
    //  2b.  Hamiltonian/Pulay gradient (off-site shell pairs)
    //
    //  Loop over unique atom pairs (iat < jat).  For each shell pair
    //  (ish_a on iat, ish_b on jat) and each AO pair (μ∈ish_a, ν∈ish_b):
    //
    //    sval = 2·P_μν·h_av − 2·W_μν − P_μν·(v_ao(μ)+v_ao(ν))
    //    G_iat += sval · dS_μν/dR_iat  +  2·P_μν·H0_μν · dlog(pi_ij) · dir
    //    G_jat -= same  [Newton's 3rd law: dS/dR_jat = −dS/dR_iat]
    //
    //  CN accumulation (off-site):
    //    dEdcn[iat] += (-kcn[ish_a]) · h_factor · Σ_μν P_μν·S_μν
    //    dEdcn[jat] += (-kcn[ish_b]) · h_factor · Σ_μν P_μν·S_μν
    // ==========================================================================
    for (int iat = 0; iat < nat; ++iat) {
        for (int jat = iat + 1; jat < nat; ++jat) {

            const double dx_ij = xyz[3*iat+0] - xyz[3*jat+0];
            const double dy_ij = xyz[3*iat+1] - xyz[3*jat+1];
            const double dz_ij = xyz[3*iat+2] - xyz[3*jat+2];
            const double r2    = dx_ij*dx_ij + dy_ij*dy_ij + dz_ij*dz_ij;
            if (r2 < 1.0e-12) continue;
            const double r     = std::sqrt(r2);

            const int zi = m_basis.z[iat];
            const int zj = m_basis.z[jat];
            const double rad_sum = atomic_rad_au(zi) + atomic_rad_au(zj);
            // rr = sqrt( r / rad_sum ),  d(rr)/dr = rr/(2r)
            const double rr = std::sqrt(r / rad_sum);

            // Loop over all shell pairs (ish_a on iat, ish_b on jat)
            const int nsh_iat = m_basis.nsh_at[iat];
            const int nsh_jat = m_basis.nsh_at[jat];

            for (int ia = 0; ia < nsh_iat; ++ia) {
                const int ish_a    = m_basis.ish_at[iat] + ia;
                const int ia_start = m_basis.iao_sh[ish_a];
                const int ia_nao   = m_basis.nao_sh[ish_a];
                const double pi_a  = 1.0 + m_h0.shpoly[ish_a] * rr;
                const double zeta_a = m_basis.cgto[ish_a].slater_exp;
                const CGTO::Shell sh_a = as_cgto_shell_g(m_basis.cgto[ish_a]);

                for (int ib = 0; ib < nsh_jat; ++ib) {
                    const int ish_b    = m_basis.ish_at[jat] + ib;
                    const int jb_start = m_basis.iao_sh[ish_b];
                    const int jb_nao   = m_basis.nao_sh[ish_b];
                    const double pi_b  = 1.0 + m_h0.shpoly[ish_b] * rr;
                    const double zeta_b = m_basis.cgto[ish_b].slater_exp;
                    const CGTO::Shell sh_b = as_cgto_shell_g(m_basis.cgto[ish_b]);

                    // ── hscale hs ────────────────────────────────────────────
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
                    const double h_factor  = hs * (pi_a * pi_b);     // = H0_μν / (avg_eps · S_μν)
                    const double avg_eps   = 0.5 * (se(ish_a) + se(ish_b));
                    const double h_av      = avg_eps * h_factor;      // = H0_μν / S_μν

                    // ── d(log pi_ij)/dr · (1/r): scalar for shpoly term ─────
                    // pi_ij = pi_a · pi_b,  d(log pi_ij)/dr = (sp_a/pi_a + sp_b/pi_b) · rr/(2r)
                    // Component along ij: dlog_pi_dr_r = dlog_pi_dr / r (pre-divided by r for dir)
                    const double dlog_pi_dr_r = (m_h0.shpoly[ish_a] / pi_a
                                                + m_h0.shpoly[ish_b] / pi_b) * rr / (2.0 * r2);

                    // ── v_ao for this shell pair ─────────────────────────────
                    // All AOs in shell ish_a share the same v_ao value.
                    const double v_a = m_pot.v_sh(ish_a) + m_pot.v_at(iat);
                    const double v_b = m_pot.v_sh(ish_b) + m_pot.v_at(jat);

                    // ── Accumulate over AO pairs ──────────────────────────────
                    double G_sval[3]       = {0.0, 0.0, 0.0};
                    double G_shpoly_scalar = 0.0;
                    double cn_sum          = 0.0;  // Σ P_μν · S_μν for this shell pair

                    for (int iao = 0; iao < ia_nao; ++iao) {
                        const int mu   = ia_start + iao;
                        const int t_a  = ao_to_type_g(m_basis.ang_sh[ish_a], iao);
                        if (t_a < 0) continue;   // d-type: gradient = 0

                        for (int jao = 0; jao < jb_nao; ++jao) {
                            const int nu  = jb_start + jao;
                            const int t_b = ao_to_type_g(m_basis.ang_sh[ish_b], jao);
                            if (t_b < 0) continue;

                            const double Pmn  = m_wfn.P(mu, nu);
                            const double Smn  = m_S(mu, nu);
                            const double H0mn = m_H0(mu, nu);
                            const double Wmn  = W(mu, nu);

                            // Overlap gradient dS/dR_iat (Obara-Saika)
                            double dS[3];
                            CGTO::cgto_overlap_grad(sh_a, sh_b,
                                                    xyz[3*iat+0], xyz[3*iat+1], xyz[3*iat+2],
                                                    xyz[3*jat+0], xyz[3*jat+1], xyz[3*jat+2],
                                                    t_a, t_b, dS);

                            // sval: coefficient for the overlap-gradient term
                            // sval = 2·P·h_av - 2·W - P·(v_a+v_b)
                            const double sval = 2.0*Pmn*h_av - 2.0*Wmn - Pmn*(v_a + v_b);

                            G_sval[0] += sval * dS[0];
                            G_sval[1] += sval * dS[1];
                            G_sval[2] += sval * dS[2];

                            // shpoly contribution: 2·P·H0 · dlog(pi)/dR = 2·P·H0 · dlog_pi_dr_r · rij
                            G_shpoly_scalar += 2.0 * Pmn * H0mn * dlog_pi_dr_r;

                            // CN accumulation: Σ P·S (factor handled outside loop)
                            cn_sum += Pmn * Smn;
                        }
                    }

                    // ── Apply to atomic gradients (Newton's 3rd law) ─────────
                    const double Gx = G_sval[0] + G_shpoly_scalar * dx_ij;
                    const double Gy = G_sval[1] + G_shpoly_scalar * dy_ij;
                    const double Gz = G_sval[2] + G_shpoly_scalar * dz_ij;
                    m_gradient(iat, 0) += Gx;  m_gradient(jat, 0) -= Gx;
                    m_gradient(iat, 1) += Gy;  m_gradient(jat, 1) -= Gy;
                    m_gradient(iat, 2) += Gz;  m_gradient(jat, 2) -= Gz;

                    // ── CN accumulation (off-site) ───────────────────────────
                    // d(avg_eps·h_factor)/d(CN_iat) = 0.5·(-kcn_a)·h_factor
                    // Full-sum factor of 2 cancels the 0.5:
                    // dEdcn[iat] += (-kcn_a) · h_factor · Σ P·S
                    dEdcn(iat) += (-m_h0.kcn[ish_a]) * h_factor * cn_sum;
                    dEdcn(jat) += (-m_h0.kcn[ish_b]) * h_factor * cn_sum;
                }
            }
        }
    }

    // ==========================================================================
    //  3.  Coulomb (ES2) gradient
    //      E_iso = 0.5·q_sh^T·γ·q_sh
    //      dE/dR_iat = Σ_{is on iat, js on jat≠iat} q_is · dγ_is,js/dR_iat · q_js
    //      (factor 0.5·2 = 1 from combining both triangles with Newton's 3rd law)
    //
    //      Klopman-Ohno kernel: γ(r) = (r^g + gam^{-g})^{-1/g}
    //      dγ/dr = −r^{g−2} · γ^{g+1}
    // ==========================================================================
    {
        const double gexp = (m_method == MethodType::GFN1)
                          ? gfn1_params::gexp
                          : gfn2_params::gexp;

        using curcuma::xtb::coulomb::shell_hardness;
        using curcuma::xtb::coulomb::harmonic_average;
        using curcuma::xtb::coulomb::arithmetic_average;
        const curcuma::xtb::coulomb::Method coul_method =
            (m_method == MethodType::GFN1)
            ? curcuma::xtb::coulomb::Method::GFN1
            : curcuma::xtb::coulomb::Method::GFN2;

        // Precompute per-shell hardness for gam_bar
        std::vector<double> g_sh(nsh);
        for (int is = 0; is < nsh; ++is) {
            const int zs  = m_basis.z[m_basis.sh2at[is]];
            const int ang = m_basis.ang_sh[is];
            g_sh[is] = shell_hardness(coul_method, zs, ang);
        }

        // Loop over unique shell pairs on different atoms
        for (int is = 0; is < nsh; ++is) {
            const int iat = m_basis.sh2at[is];
            for (int js = 0; js < is; ++js) {
                const int jat = m_basis.sh2at[js];
                if (iat == jat) continue;   // on-site: no R-dependence

                const double dx = xyz[3*iat+0] - xyz[3*jat+0];
                const double dy = xyz[3*iat+1] - xyz[3*jat+1];
                const double dz = xyz[3*iat+2] - xyz[3*jat+2];
                const double r2 = dx*dx + dy*dy + dz*dz;
                if (r2 < 1.0e-12) continue;
                const double r1 = std::sqrt(r2);

                // gam_bar: harmonic (GFN1) or arithmetic (GFN2) average
                const double gam_bar = (m_method == MethodType::GFN1)
                    ? harmonic_average(g_sh[is], g_sh[js])
                    : arithmetic_average(g_sh[is], g_sh[js]);

                // γ(r) = (r^g + gam_bar^{-g})^{-1/g}
                const double r1g     = std::pow(r1, gexp);
                const double gam_pow = std::pow(gam_bar, -gexp);
                const double gamma   = std::pow(r1g + gam_pow, -1.0 / gexp);

                // dγ/dr = −r^{g−2}·γ^{g+1};  force scalar: q_is·q_js·dγ/dr
                const double dgamma_dr = -std::pow(r1, gexp - 2.0) * std::pow(gamma, gexp + 1.0);
                const double force     = m_wfn.q_sh(is) * m_wfn.q_sh(js) * dgamma_dr;

                // dγ/dR_iat[k] = dgamma_dr · (R_iat − R_jat)[k] / r1
                // G_iat += q_is·q_js·dgamma_dr · (R_iat−R_jat)
                // (here force already absorbed q_is·q_js·dgamma_dr; multiply by rij)
                m_gradient(iat, 0) += force * dx;  m_gradient(jat, 0) -= force * dx;
                m_gradient(iat, 1) += force * dy;  m_gradient(jat, 1) -= force * dy;
                m_gradient(iat, 2) += force * dz;  m_gradient(jat, 2) -= force * dz;
            }
        }
    }

    // ==========================================================================
    //  4.  CN chain-rule gradient
    //      G_iat += Σ_j (dEdcn[i] + dEdcn[j]) · dcount_ij/dr · (R_i−R_j)/r
    //
    //  GFN1 (exp counting):  count = 1/(1+exp(−k·(rc/r−1))), k=16
    //    dc/dr = −c·(1−c)·k·rc/r²
    //
    //  GFN2 (double-exp):    count = c1(ka,rc)·c2(kb,rc+rs)
    //    dc/dr = dc1/dr·c2 + c1·dc2/dr
    //
    //  Covalent radii: D3 parametrisation via covalent_rad_d3_au().
    // ==========================================================================
    {
        std::vector<double> rcov(nat);
        for (int i = 0; i < nat; ++i)
            rcov[i] = covalent_rad_d3_au(m_atoms[i]);

        for (int i = 0; i < nat; ++i) {
            for (int j = 0; j < i; ++j) {
                const double dx = xyz[3*i+0] - xyz[3*j+0];
                const double dy = xyz[3*i+1] - xyz[3*j+1];
                const double dz = xyz[3*i+2] - xyz[3*j+2];
                const double r2 = dx*dx + dy*dy + dz*dz;
                if (r2 < 1.0e-12) continue;
                const double r  = std::sqrt(r2);
                const double rc = rcov[i] + rcov[j];

                double dcndr;   // d(count_ij)/dr  (both CN_i and CN_j contain this pair)
                if (m_method == MethodType::GFN1) {
                    // count = 1/(1 + exp(−k·(rc/r−1))),  k=16
                    // dc/d(arg) = −c·(1−c);  d(arg)/dr = k·rc/r²
                    constexpr double k = 16.0;
                    const double arg = -k * (rc / r - 1.0);
                    const double ef  = std::exp(arg);
                    const double c   = 1.0 / (1.0 + ef);
                    dcndr = -c * (1.0 - c) * k * rc / r2;
                } else {
                    // count = c1(ka,rc)·c2(kb,rc+rs),  ka=10, kb=20, rs=2 Bohr
                    constexpr double ka = 10.0, kb = 20.0, rs = 2.0;
                    const double arg1 = -ka * (rc       / r - 1.0);
                    const double arg2 = -kb * ((rc + rs) / r - 1.0);
                    const double ef1  = std::exp(arg1);
                    const double ef2  = std::exp(arg2);
                    const double c1   = 1.0 / (1.0 + ef1);
                    const double c2   = 1.0 / (1.0 + ef2);
                    // dc1/dr = −c1·(1−c1)·ka·rc/r²
                    // dc2/dr = −c2·(1−c2)·kb·(rc+rs)/r²
                    const double dc1dr = -c1 * (1.0 - c1) * ka * rc       / r2;
                    const double dc2dr = -c2 * (1.0 - c2) * kb * (rc + rs) / r2;
                    dcndr = dc1dr * c2 + c1 * dc2dr;
                }

                // Both CN_i and CN_j receive a contribution from pair (i,j).
                // dCN_i/dR_i = dCN_j/dR_i = dcndr · rij/r  (r decreases as i→j)
                // dE/dR_i[k] = (dEdcn[i]+dEdcn[j]) · dcndr · rij[k]/r
                const double factor = (dEdcn(i) + dEdcn(j)) * dcndr / r;
                m_gradient(i, 0) += factor * dx;  m_gradient(j, 0) -= factor * dx;
                m_gradient(i, 1) += factor * dy;  m_gradient(j, 1) -= factor * dy;
                m_gradient(i, 2) += factor * dz;  m_gradient(j, 2) -= factor * dz;
            }
        }
    }

    // ==========================================================================
    //  5.  GFN2 multipole gradient
    //      TODO AP5: analytical gradient for sd/dd/sq interactions.
    //      Deferred together with GFN2 energy accuracy work (AP5).
    //      Impact: small relative to repulsion+H0+Coulomb; -opt converges.
    // ==========================================================================
}

} // namespace curcuma::xtb
