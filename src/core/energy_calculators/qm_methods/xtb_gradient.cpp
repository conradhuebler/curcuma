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
 *   - external/tblite/src/tblite/coulomb/multipole.f90  (get_multipole_gradient_0d)
 *
 * AP4: repulsion + H0/Pulay (s/p only) + Coulomb + CN chain-rule.
 * AP5: GFN2 direct multipole gradient (SD/DD/SQ) + mrad/CN chain-rule.
 * AP5 (TODO): integral Pulay term (multipole integral derivatives × v_dp/v_qp).
 *
 * Claude Generated. GPL-3.0.
 */

#include "xtb_native.h"
#include "STO_CGTO.hpp"
#include "xtb_ao_utils.hpp"
#include "xtb_coulomb.hpp"
#include "xtb_multipole_ints.hpp"

#include "parameters/gfn1_params.hpp"
#include "parameters/gfn2_params.hpp"
#include "parameters/xtb_params_extra.hpp"

#include "src/core/curcuma_logger.h"

#include <cmath>

namespace curcuma::xtb {

/* as_cgto_shell() and ao_to_type() now live in xtb_ao_utils.hpp (X-I5). */

/* ========================================================================== *
 *  XTB::calculateGradient()                                                   *
 *                                                                             *
 *  Fills m_gradient (nat×3) in Eh/Bohr.  Call after converged SCF.           *
 *                                                                             *
 *  Components:                                                                *
 *   1.  Repulsion gradient                                                    *
 *   2.  Hamiltonian/Pulay gradient (H0 + W + potential terms)                *
 *   3.  Coulomb (ES2) gradient                                                *
 *   5.  GFN2 multipole gradient (AP5) — must run before 4 to fill dEdcn      *
 *   4.  CN chain-rule gradient                                                *
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
    // W = C_occ · diag(2·ε_occ) · C_occᵀ as one BLAS-3 gemm (WP1, Claude Generated):
    // replaces the per-orbital rank-1 accumulation — faster serial and threads under MKL.
    // The purification density path (eigensolver="purify") has no eps/C and instead supplies W
    // directly (W = 2·L⁻ᵀ·P̃·Ã·P̃·L⁻¹); use it when present. Claude Generated.
    const int nocc_orbs = static_cast<int>(std::round(m_wfn.nocc / 2.0));
    Matrix W;
    if (m_wfn.W_valid && m_wfn.W.rows() == nao && m_wfn.W.cols() == nao) {
        W = m_wfn.W;
    } else if (m_wfn.focc.size() == nao && m_wfn.C.cols() == nao) {
        // X-G3 (Claude Generated): occupation-consistent W = Σ_i f_i·ε_i·C_i·C_iᵀ using the
        // same per-MO occupations that built the density. At electronic_temperature > 0 these
        // are fractional (Fermi) so the Pulay term matches the fractional density; at T = 0
        // f_i ∈ {2,0}, so W is bit-identical to the previous integer build. Restrict to the
        // leading columns with non-negligible occupation (eigenvalues ascending).
        int ncol = 0;
        for (int i = 0; i < nao; ++i)
            if (m_wfn.focc(i) > 1.0e-12) ncol = i + 1;
        const auto Cocc = m_wfn.C.leftCols(ncol);
        const Eigen::VectorXd we =
            (m_wfn.focc.head(ncol).array() * m_wfn.eps.head(ncol).array()).matrix();
        W = (Cocc * we.asDiagonal()) * Cocc.transpose();
    } else {
        // Fallback when focc was not populated (e.g. a device density path): integer occupation.
        const auto Cocc = m_wfn.C.leftCols(nocc_orbs);
        const Eigen::VectorXd w2 = 2.0 * m_wfn.eps.head(nocc_orbs);
        W = (Cocc * w2.asDiagonal()) * Cocc.transpose();
    }

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
    //      E_rep_ij = Z_i·Z_j / R^rexp * exp(-sqrt(α_i·α_j) * R^kexp_ij)
    //      dE/dr    = -(rexp/R + alpha_pair*kexp_ij*R^(kexp_ij-1)) * E_pair
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

            const double alpha_pair = std::sqrt(alfi * alfj);
            double kexp_pair = kexp;
            if (m_method == MethodType::GFN2) {
                if (zi <= 2 && zj <= 2)
                    kexp_pair = gfn2_params::rep_kexp_light;
            }
            const double r_kexp = std::pow(r, kexp_pair);
            const double E_pair = zeffi * zeffj / std::pow(r, rexp)
                                * std::exp(-alpha_pair * r_kexp);
            const double dEdr   = -(rexp / r + alpha_pair * kexp_pair
                                    * std::pow(r, kexp_pair - 1.0)) * E_pair;

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
    // Parallel over iat (Claude Generated). This is the dominant gradient cost
    // (per AO-pair overlap- and multipole-integral derivatives). Each thread writes
    // BOTH iat and jat (Newton's 3rd law) and accumulates dEdcn(iat)/dEdcn(jat),
    // which span thread boundaries — so threads write into private grad/dEdcn
    // partials that are summed afterwards. Bit-identical to the serial sum up to
    // floating-point reassociation (validated against the serial reference).
    const int grad_threads = effectiveIntraThreads(nat);
    std::vector<Matrix> grad_parts(grad_threads, Matrix::Zero(nat, 3));
    std::vector<Vector> dEdcn_parts(grad_threads, Vector::Zero(nat));
    parallelStripes(grad_threads, [&](int tid, int nth) {
    Matrix& g_loc    = grad_parts[tid];
    Vector& dEdcn_loc = dEdcn_parts[tid];
    for (int iat = tid; iat < nat; iat += nth) {
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
                const CGTO::Shell sh_a = as_cgto_shell(m_basis.cgto[ish_a]);

                for (int ib = 0; ib < nsh_jat; ++ib) {
                    const int ish_b    = m_basis.ish_at[jat] + ib;
                    const int jb_start = m_basis.iao_sh[ish_b];
                    const int jb_nao   = m_basis.nao_sh[ish_b];
                    const double pi_b  = 1.0 + m_h0.shpoly[ish_b] * rr;
                    const double zeta_b = m_basis.cgto[ish_b].slater_exp;
                    const CGTO::Shell sh_b = as_cgto_shell(m_basis.cgto[ish_b]);

                    // X-I1: d-touching shell pair uses the spherical-transform
                    // gradient blocks (computed once per shell pair); pure s/p
                    // pairs keep the scalar kernels below (byte-identical).
                    const int ang_a_g = m_basis.ang_sh[ish_a];
                    const int ang_b_g = m_basis.ang_sh[ish_b];
                    const bool dpair = (ang_a_g >= 2 || ang_b_g >= 2);
                    double dSblk[5*5*3];
                    double dDblk[5*5*9], dQblk[5*5*18];
                    if (dpair) {
                        sphericalOverlapGradBlock(sh_a, ang_a_g, sh_b, ang_b_g,
                            xyz[3*iat+0], xyz[3*iat+1], xyz[3*iat+2],
                            xyz[3*jat+0], xyz[3*jat+1], xyz[3*jat+2], dSblk, 5);
                        if (m_method == MethodType::GFN2 && m_mp_initialized)
                            sphericalMultipoleGradBlock(sh_a, ang_a_g, sh_b, ang_b_g,
                                xyz[3*iat+0], xyz[3*iat+1], xyz[3*iat+2],
                                xyz[3*jat+0], xyz[3*jat+1], xyz[3*jat+2], dDblk, dQblk, 5);
                    }

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
                        const int t_a  = dpair ? 0 : ao_to_type(m_basis.ang_sh[ish_a], iao);
                        if (!dpair && t_a < 0) continue;

                        for (int jao = 0; jao < jb_nao; ++jao) {
                            const int nu  = jb_start + jao;
                            const int t_b = dpair ? 0 : ao_to_type(m_basis.ang_sh[ish_b], jao);
                            if (!dpair && t_b < 0) continue;

                            const double Pmn  = m_wfn.P(mu, nu);
                            const double Smn  = m_S(mu, nu);
                            const double H0mn = m_H0(mu, nu);
                            const double Wmn  = W(mu, nu);

                            // Overlap gradient dS/dR_iat (Obara-Saika)
                            double dS[3];
                            if (dpair) {
                                const int o = (iao*5 + jao) * 3;
                                dS[0] = dSblk[o]; dS[1] = dSblk[o+1]; dS[2] = dSblk[o+2];
                            } else {
                                CGTO::cgto_overlap_grad(sh_a, sh_b,
                                                        xyz[3*iat+0], xyz[3*iat+1], xyz[3*iat+2],
                                                        xyz[3*jat+0], xyz[3*jat+1], xyz[3*jat+2],
                                                        t_a, t_b, dS);
                            }

                            // sval: coefficient for the overlap-gradient term
                            // sval = 2·P·h_av - 2·W - P·(v_a+v_b)
                            const double sval = 2.0*Pmn*h_av - 2.0*Wmn - Pmn*(v_a + v_b);

                            G_sval[0] += sval * dS[0];
                            G_sval[1] += sval * dS[1];
                            G_sval[2] += sval * dS[2];

                            // AP5b: GFN2 multipole integral Pulay term
                            // Fock: F_μν -= 0.5*(dp_jat·v_dp(jat) + dp_iat·v_dp(iat) + quad terms)
                            //   dp_jat = D with origin at jat (column atom)
                            //   dp_iat = dp_jat + ΔR·S  (ΔR = Rjat - Riat)
                            // Gradient (Pulay, iat-derivative, triangular loop):
                            //   term = d(dp_jat)/dRiat · v_dp(jat) + d(dp_iat)/dRiat · v_dp(iat)
                            //        + d(qp_jat)/dRiat · v_qp(jat) + d(qp_iat)/dRiat · v_qp(iat)
                            // Reference: tblite h0.f90:get_hamiltonian_gradient
                            if (m_method == MethodType::GFN2 && m_mp_initialized) {
                                double D_mp[3], Q_mp[6];
                                double dD_dA[3][3], dQ_dA[3][6];
                                if (dpair) {
                                    // Transformed integral values come from the stored
                                    // m_dp_int/m_qp_int (origin at column atom = jat);
                                    // their A-gradients from the d block (X-I1).
                                    for (int k = 0; k < 3; ++k) D_mp[k] = m_dp_int[k](mu, nu);
                                    for (int q = 0; q < 6; ++q) Q_mp[q] = m_qp_int[q](mu, nu);
                                    const int od = (iao*5 + jao) * 9;
                                    const int oq = (iao*5 + jao) * 18;
                                    for (int l = 0; l < 3; ++l) {
                                        for (int k = 0; k < 3; ++k) dD_dA[l][k] = dDblk[od + l*3 + k];
                                        for (int q = 0; q < 6; ++q) dQ_dA[l][q] = dQblk[oq + l*6 + q];
                                    }
                                } else {
                                    double S_mp, dD_dB[3][3], dQ_dB[3][6];
                                    using namespace curcuma::xtb::multipole_ints;
                                    cgto_multipole_grad_transformed(
                                        sh_a, sh_b,
                                        xyz[3*iat+0], xyz[3*iat+1], xyz[3*iat+2],
                                        xyz[3*jat+0], xyz[3*jat+1], xyz[3*jat+2],
                                        t_a, t_b,
                                        S_mp, D_mp, Q_mp,
                                        dD_dA, dD_dB, dQ_dA, dQ_dB);
                                }

                                // ΔR = Rjat - Riat (origin-shift vector)
                                const double dR[3] = {
                                    xyz[3*jat+0]-xyz[3*iat+0],
                                    xyz[3*jat+1]-xyz[3*iat+1],
                                    xyz[3*jat+2]-xyz[3*iat+2]
                                };
                                // Packed quadrupole index mapping (upper triangle)
                                static const int qa6[6] = {0,0,1,0,1,2};
                                static const int qb6[6] = {0,1,1,2,2,2};

                                for (int l = 0; l < 3; ++l) {
                                    double term = 0.0;

                                    // -- Dipole contributions --
                                    for (int k = 0; k < 3; ++k) {
                                        // d(dp_jat[k])/dRiat_l · v_dp(jat)
                                        term += dD_dA[l][k] * m_pot.v_dp(k, jat);
                                        // d(dp_iat[k])/dRiat_l · v_dp(iat)
                                        // = (dD_dA[l][k] + ΔR[k]·dS[l] - δ_{kl}·Smn) · v_dp(iat)
                                        const double dDiat = dD_dA[l][k] + dR[k]*dS[l] - (k==l ? Smn : 0.0);
                                        term += dDiat * m_pot.v_dp(k, iat);
                                    }

                                    // -- Quadrupole contributions --
                                    // d(qp_iat) = d(qp_jat) + traceless(Δqraw)
                                    // Δqraw[q] = origin-shift correction for iat-origin quad
                                    double Δqraw[6];
                                    for (int q = 0; q < 6; ++q) {
                                        const int a = qa6[q], b = qb6[q];
                                        // derivative of (ΔR[b]*D_mp[a] + ΔR[a]*D_mp[b] + ΔR[a]*ΔR[b]*S)
                                        // w.r.t. Riat_l  (dΔR[k]/dRiat_l = -δ_{kl})
                                        Δqraw[q] = -(b==l ? D_mp[a] : 0.0) + dR[b]*dD_dA[l][a]
                                                 - (a==l ? D_mp[b] : 0.0) + dR[a]*dD_dA[l][b]
                                                 + (-(a==l ? dR[b] : 0.0) - (b==l ? dR[a] : 0.0))*Smn
                                                 + dR[a]*dR[b]*dS[l];
                                    }
                                    const double dtr_c = 0.5*(Δqraw[0] + Δqraw[2] + Δqraw[5]);

                                    for (int q = 0; q < 6; ++q) {
                                        const bool is_diag = (qa6[q] == qb6[q]);
                                        // d(qp_jat[q])/dRiat_l · v_qp(jat)
                                        term += dQ_dA[l][q] * m_pot.v_qp(q, jat);
                                        // d(qp_iat[q])/dRiat_l · v_qp(iat)
                                        const double dQiat = dQ_dA[l][q]
                                                           + 1.5*Δqraw[q]
                                                           - (is_diag ? dtr_c : 0.0);
                                        term += dQiat * m_pot.v_qp(q, iat);
                                    }

                                    G_sval[l] -= Pmn * term;
                                }
                            }

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
                    g_loc(iat, 0) += Gx;  g_loc(jat, 0) -= Gx;
                    g_loc(iat, 1) += Gy;  g_loc(jat, 1) -= Gy;
                    g_loc(iat, 2) += Gz;  g_loc(jat, 2) -= Gz;

                    // ── CN accumulation (off-site) ───────────────────────────
                    // d(avg_eps·h_factor)/d(CN_iat) = 0.5·(-kcn_a)·h_factor
                    // Full-sum factor of 2 cancels the 0.5:
                    // dEdcn[iat] += (-kcn_a) · h_factor · Σ P·S
                    dEdcn_loc(iat) += (-m_h0.kcn[ish_a]) * h_factor * cn_sum;
                    dEdcn_loc(jat) += (-m_h0.kcn[ish_b]) * h_factor * cn_sum;
                }
            }
        }
    }
    });  // parallelStripes over iat
    // Reduce the per-thread partials into the shared accumulators. dEdcn must be
    // complete here: sections 5 (multipole mrad) and 4 (CN chain rule) read it.
    for (int t = 0; t < grad_threads; ++t) {
        m_gradient += grad_parts[t];
        dEdcn      += dEdcn_parts[t];
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
    //  5.  GFN2 multipole gradient (AP5)
    //      Port of tblite/coulomb/multipole.f90:get_multipole_gradient_0d
    //      Reference: Bannwarth et al., JCTC 2019; Grimme et al., JCTC 2017
    //
    //  Three pairwise terms (unique pairs iat > jat, vec = R_jat - R_iat):
    //
    //    SD (charge-dipole, g3):
    //      dpiqj = (vec·dp_iat) * q_jat
    //      qidpj = (vec·dp_jat) * q_iat
    //      ddmp3 = g5*fdmp3*(-3 + mp_dmp3*(1-fdmp3))          [includes g5]
    //      dg    = -ddmp3*vec*(dpiqj-qidpj) + fdmp3*g3*(q_iat*dp_jat - q_jat*dp_iat)
    //
    //    DD (dipole-dipole, g5):
    //      edd   = (dp_iat·dp_jat)*r² - 3*(dp_jat·vec)*(dp_iat·vec)
    //      ddmp5 = fdmp5*(-5 + mp_dmp5*(1-fdmp5))             [without g7]
    //      dg   += -2*fdmp5*g5*(dp_iat·dp_jat)*vec
    //              +3*fdmp5*g5*((dp_iat·vec)*dp_jat + (dp_jat·vec)*dp_iat)
    //              -edd*ddmp5*g7*vec
    //
    //    SQ (charge-quadrupole, g5):
    //      eq    = (q_jat*qp_iat + qp_jat*q_iat) : vec⊗vec   [packed sym tensor]
    //      θ·v   = symmetric matrix-vector product (packed xx,xy,yy,xz,yz,zz)
    //      dg   += -eq*ddmp5*g7*vec
    //              -2*fdmp5*g5*(q_iat*(θ_jat·vec) + q_jat*(θ_iat·vec))
    //
    //  mrad/CN chain rule:
    //    dEdr_mp[i] += fddr_{SD,DD,SQ}  (derivative w.r.t. mrad[i])
    //    dmrdcn[i]   = -t2*mp_kexp*t1/(1+t1)  (same formula as in setupMultipole)
    //    dEdcn[i]   += dEdr_mp[i] * dmrdcn[i]  → fed into existing CN pair loop (sect. 4)
    //
    //  Must run BEFORE section 4 so that dEdcn is complete when CN loop executes.
    // ==========================================================================
    if (m_method == MethodType::GFN2 && m_mp_initialized) {
        using namespace gfn2_params;

        Vector dEdr_mp = Vector::Zero(nat);

        for (int iat = 0; iat < nat; ++iat) {
            for (int jat = 0; jat < iat; ++jat) {
                // vec = R_jat - R_iat  (TBLite convention)
                const double vx = xyz[3*jat+0] - xyz[3*iat+0];
                const double vy = xyz[3*jat+1] - xyz[3*iat+1];
                const double vz = xyz[3*jat+2] - xyz[3*iat+2];
                const double r2 = vx*vx + vy*vy + vz*vz;
                if (r2 < 1.0e-12) continue;
                const double r1 = std::sqrt(r2);
                const double g1 = 1.0 / r1;
                const double g3 = g1 * g1 * g1;
                const double g5 = g3 * g1 * g1;
                const double g7 = g5 * g1 * g1;

                const double rr_dmp = 0.5 * (m_mp_mrad[iat] + m_mp_mrad[jat]);

                const double fdmp3 = 1.0 / (1.0 + 6.0 * std::pow(rr_dmp * g1, mp_dmp3));
                const double fdmp5 = 1.0 / (1.0 + 6.0 * std::pow(rr_dmp * g1, mp_dmp5));

                // ddmp3: d(g3*fdmp3)/dr divided by (vec/r), includes g5 factor
                const double ddmp3 = g5 * fdmp3 * (-3.0 + mp_dmp3 * (1.0 - fdmp3));
                // ddmp5: same for g5*fdmp5, WITHOUT g7 (g7 applied separately)
                const double ddmp5 = fdmp5 * (-5.0 - mp_dmp5 * (fdmp5 - 1.0));

                // ── SD (charge-dipole) ────────────────────────────────────────────
                const double dpiqj = (vx*m_wfn.dp_at(0,iat) + vy*m_wfn.dp_at(1,iat)
                                    + vz*m_wfn.dp_at(2,iat)) * m_wfn.q_at(jat);
                const double qidpj = (vx*m_wfn.dp_at(0,jat) + vy*m_wfn.dp_at(1,jat)
                                    + vz*m_wfn.dp_at(2,jat)) * m_wfn.q_at(iat);
                const double diff_sd = dpiqj - qidpj;

                double gx = -ddmp3 * vx * diff_sd
                          + fdmp3 * g3 * (m_wfn.q_at(iat)*m_wfn.dp_at(0,jat) - m_wfn.q_at(jat)*m_wfn.dp_at(0,iat));
                double gy = -ddmp3 * vy * diff_sd
                          + fdmp3 * g3 * (m_wfn.q_at(iat)*m_wfn.dp_at(1,jat) - m_wfn.q_at(jat)*m_wfn.dp_at(1,iat));
                double gz = -ddmp3 * vz * diff_sd
                          + fdmp3 * g3 * (m_wfn.q_at(iat)*m_wfn.dp_at(2,jat) - m_wfn.q_at(jat)*m_wfn.dp_at(2,iat));

                // fddr_SD: ∂E_SD/∂mrad (same value added to both iat and jat)
                const double fddr_sd = 3.0 * diff_sd * mp_dmp3 * fdmp3 * g3
                                     * (fdmp3 / rr_dmp) * std::pow(rr_dmp * g1, mp_dmp3);
                dEdr_mp(iat) += fddr_sd;
                dEdr_mp(jat) += fddr_sd;

                // ── DD (dipole-dipole) ────────────────────────────────────────────
                const double dpidpj = m_wfn.dp_at(0,iat)*m_wfn.dp_at(0,jat)
                                    + m_wfn.dp_at(1,iat)*m_wfn.dp_at(1,jat)
                                    + m_wfn.dp_at(2,iat)*m_wfn.dp_at(2,jat);
                const double dpiv   = vx*m_wfn.dp_at(0,iat) + vy*m_wfn.dp_at(1,iat) + vz*m_wfn.dp_at(2,iat);
                const double dpjv   = vx*m_wfn.dp_at(0,jat) + vy*m_wfn.dp_at(1,jat) + vz*m_wfn.dp_at(2,jat);
                const double edd    = dpidpj * r2 - 3.0 * dpjv * dpiv;

                gx += -2.0*fdmp5*g5*dpidpj*vx
                    + 3.0*fdmp5*g5*(dpiv*m_wfn.dp_at(0,jat) + dpjv*m_wfn.dp_at(0,iat))
                    - edd * ddmp5 * g7 * vx;
                gy += -2.0*fdmp5*g5*dpidpj*vy
                    + 3.0*fdmp5*g5*(dpiv*m_wfn.dp_at(1,jat) + dpjv*m_wfn.dp_at(1,iat))
                    - edd * ddmp5 * g7 * vy;
                gz += -2.0*fdmp5*g5*dpidpj*vz
                    + 3.0*fdmp5*g5*(dpiv*m_wfn.dp_at(2,jat) + dpjv*m_wfn.dp_at(2,iat))
                    - edd * ddmp5 * g7 * vz;

                const double fddr_dd = 3.0 * edd * mp_dmp5 * fdmp5 * g5
                                     * (fdmp5 / rr_dmp) * std::pow(rr_dmp * g1, mp_dmp5);
                dEdr_mp(iat) += fddr_dd;
                dEdr_mp(jat) += fddr_dd;

                // ── SQ (charge-quadrupole) ────────────────────────────────────────
                // qp_at packing: (xx=0, xy=1, yy=2, xz=3, yz=4, zz=5)
                // eq = Σ_{ab} (q_j*θ_i + θ_j*q_i)_{ab} * vec_a * vec_b  [packed sym]
                const double qi = m_wfn.q_at(iat);
                const double qj = m_wfn.q_at(jat);
                const double eq =
                      (qj*m_wfn.qp_at(0,iat) + m_wfn.qp_at(0,jat)*qi) * vx*vx
                    + 2.0*(qj*m_wfn.qp_at(1,iat) + m_wfn.qp_at(1,jat)*qi) * vx*vy
                    + (qj*m_wfn.qp_at(2,iat) + m_wfn.qp_at(2,jat)*qi) * vy*vy
                    + 2.0*(qj*m_wfn.qp_at(3,iat) + m_wfn.qp_at(3,jat)*qi) * vx*vz
                    + 2.0*(qj*m_wfn.qp_at(4,iat) + m_wfn.qp_at(4,jat)*qi) * vy*vz
                    + (qj*m_wfn.qp_at(5,iat) + m_wfn.qp_at(5,jat)*qi) * vz*vz;

                // θ_iat · vec  and  θ_jat · vec  (symmetric matrix-vector, packed)
                const double tivx = m_wfn.qp_at(0,iat)*vx + m_wfn.qp_at(1,iat)*vy + m_wfn.qp_at(3,iat)*vz;
                const double tivy = m_wfn.qp_at(1,iat)*vx + m_wfn.qp_at(2,iat)*vy + m_wfn.qp_at(4,iat)*vz;
                const double tivz = m_wfn.qp_at(3,iat)*vx + m_wfn.qp_at(4,iat)*vy + m_wfn.qp_at(5,iat)*vz;
                const double tjvx = m_wfn.qp_at(0,jat)*vx + m_wfn.qp_at(1,jat)*vy + m_wfn.qp_at(3,jat)*vz;
                const double tjvy = m_wfn.qp_at(1,jat)*vx + m_wfn.qp_at(2,jat)*vy + m_wfn.qp_at(4,jat)*vz;
                const double tjvz = m_wfn.qp_at(3,jat)*vx + m_wfn.qp_at(4,jat)*vy + m_wfn.qp_at(5,jat)*vz;

                gx += -eq * ddmp5 * g7 * vx
                    - 2.0*fdmp5*g5 * (qi*tjvx + qj*tivx);
                gy += -eq * ddmp5 * g7 * vy
                    - 2.0*fdmp5*g5 * (qi*tjvy + qj*tivy);
                gz += -eq * ddmp5 * g7 * vz
                    - 2.0*fdmp5*g5 * (qi*tjvz + qj*tivz);

                const double fddr_sq = eq * 3.0 * mp_dmp5 * fdmp5 * g5
                                     * (fdmp5 / rr_dmp) * std::pow(rr_dmp * g1, mp_dmp5);
                dEdr_mp(iat) += fddr_sq;
                dEdr_mp(jat) += fddr_sq;

                // ── Apply forces (Newton's 3rd law) ───────────────────────────────
                m_gradient(iat, 0) += gx;  m_gradient(jat, 0) -= gx;
                m_gradient(iat, 1) += gy;  m_gradient(jat, 1) -= gy;
                m_gradient(iat, 2) += gz;  m_gradient(jat, 2) -= gz;
            }
        }

        // mrad/CN chain rule:
        //   mrad[i] = p_rad[z-1] + (mp_rmax - p_rad[z-1]) / (1 + t1)
        //   dmrdcn[i] = -t2 * mp_kexp * t1 / (1 + t1)   (tblite get_mrad convention)
        //   dEdcn[i] += dEdr_mp[i] * dmrdcn[i]
        for (int i = 0; i < nat; ++i) {
            const int zi     = m_atoms[i];
            const double rad = p_rad[zi - 1];
            const double vcn = p_vcn[zi - 1];
            const double arg = m_coordination_numbers(i) - vcn - mp_shift;
            const double t1  = std::exp(-mp_kexp * arg);
            const double t2  = (mp_rmax - rad) / (1.0 + t1);
            const double dmrdcn = -t2 * mp_kexp * t1 / (1.0 + t1);
            dEdcn(i) += dEdr_mp(i) * dmrdcn;
        }
    }

    // ==========================================================================
    //  3b.  D4 dispersion gradient (GFN2 only) — Claude Generated 2026
    //
    //  The D4 energy was computed in xtb_native.cpp::calcDispersionEnergy().
    //  That function also produced two cached side-products with with_gradient=true:
    //    m_disp_gradient[i,:]  — geometry term  ∂E_D4/∂R_i  (Eh/Bohr)
    //    m_disp_dEdcn(i)       — chain-rule accumulator  dE_D4/dCN_i
    //
    //  The geometry term is added directly here; the CN-chain-rule contribution
    //  is folded into dEdcn so the existing CN-distribution loop in section 4
    //  carries it through ∂CN/∂x for free — no extra distance evaluations.
    //
    //  Caveat: the q-response chain rule (∂E_D4/∂q · ∂q/∂R via SCF response)
    //  is not yet implemented; zetac6 is treated as a static prefactor (sub-mEh
    //  residual on the AP test set). See dispersion/CLAUDE.md.
    // ==========================================================================
    if (m_disp_gradient_valid
        && m_disp_gradient.rows() == nat && m_disp_gradient.cols() == 3) {
        m_gradient += m_disp_gradient;
        if (m_disp_dEdcn.size() == nat) {
            dEdcn += m_disp_dEdcn;
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
    //  Note: dEdcn now includes H0/Pulay + GFN2 multipole mrad contributions.
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
}

/* ============================================================================ *
 *  GFN1 device nuclear gradient (Stage 4a, Claude Generated)
 *
 *  Sections 1 (repulsion), 2a/2b (on-site CN + H0/Pulay) and 3 (isotropic
 *  Coulomb) run on the GPU via the GpuScfBackend (device kernels mirror the
 *  corresponding host loops above). The dispersion gradient (3b, host-cached)
 *  and the CN chain-rule (4) stay on the host — identical formulas to
 *  calculateGradient. Fills m_gradient in Eh/Bohr; the caller applies /au.
 * ============================================================================ */
bool XTB::calculateGradientGpu()
{
    // X-I1: s/p-only device gradient kernels -> CPU gradient for d systems.
    if (!m_gpu_scf || !m_gpu_scf->supportsGradient() || m_has_dshell) return false;
    const int nat = m_atomcount;
    const int nao = m_basis.nao;

    // Converged AO potential (v_sh + v_at expanded to AO resolution).
    Vector v_ao = Vector::Zero(nao);
    for (int ao = 0; ao < nao; ++ao)
        v_ao(ao) = m_pot.v_sh(m_basis.ao2sh[ao]) + m_pot.v_at(m_basis.ao2at[ao]);

    const int nocc_orbs = static_cast<int>(std::round(m_wfn.nocc / 2.0));

    // GFN2 multipole potentials for the device multipole-integral Pulay term.
    Eigen::MatrixXd v_dp_empty, v_qp_empty;
    const Eigen::MatrixXd& v_dp = (m_method == MethodType::GFN2) ? m_pot.v_dp : v_dp_empty;
    const Eigen::MatrixXd& v_qp = (m_method == MethodType::GFN2) ? m_pot.v_qp : v_qp_empty;

    // AP8 (Claude Generated): calculateGradientGpu only runs on the device-resident
    // SCF path, so dP/dC are already the converged values on the GPU — reuse them
    // (pc_resident=true) and pass empty P/C, skipping the two nao²-sized uploads and
    // the host column-major copy. Bit-identical (the gradient only reads dP/dC).
    const Matrix P_empty;
    const Eigen::MatrixXd C_empty;
    Matrix grad_dev;   // nat×3, Eh/Bohr (sections 1/2/3 incl. GFN2 multipole Pulay)
    Vector dEdcn;      // nat
    if (!m_gpu_scf->gradient(P_empty, C_empty, m_wfn.eps, nocc_orbs, v_ao, m_wfn.q_sh,
                             v_dp, v_qp, grad_dev, dEdcn, /*pc_resident=*/true))
        return false;
    if (grad_dev.rows() != nat || grad_dev.cols() != 3 || dEdcn.size() != nat)
        return false;

    m_gradient = grad_dev;

    // Geometry in Bohr (for the host sections below).
    std::vector<double> xyz(3 * nat);
    for (int i = 0; i < nat; ++i) {
        xyz[3*i+0] = m_geometry(i, 0) * AA_TO_AU;
        xyz[3*i+1] = m_geometry(i, 1) * AA_TO_AU;
        xyz[3*i+2] = m_geometry(i, 2) * AA_TO_AU;
    }

    // 3b. Dispersion gradient (host-cached) + its CN chain-rule contribution.
    if (m_disp_gradient_valid
        && m_disp_gradient.rows() == nat && m_disp_gradient.cols() == 3) {
        m_gradient += m_disp_gradient;
        if (m_disp_dEdcn.size() == nat) dEdcn += m_disp_dEdcn;
    }

    // 5. GFN2 direct multipole interaction gradient (host; O(nat²) over the
    // atomic moments dp_at/qp_at — identical to section 5 in calculateGradient).
    // Must precede section 4 so dEdcn carries the mrad/CN chain contribution.
    if (m_method == MethodType::GFN2 && m_mp_initialized) {
        using namespace gfn2_params;
        Vector dEdr_mp = Vector::Zero(nat);
        for (int iat = 0; iat < nat; ++iat) {
            for (int jat = 0; jat < iat; ++jat) {
                const double vx = xyz[3*jat+0] - xyz[3*iat+0];
                const double vy = xyz[3*jat+1] - xyz[3*iat+1];
                const double vz = xyz[3*jat+2] - xyz[3*iat+2];
                const double r2 = vx*vx + vy*vy + vz*vz;
                if (r2 < 1.0e-12) continue;
                const double r1 = std::sqrt(r2);
                const double g1 = 1.0 / r1;
                const double g3 = g1*g1*g1, g5 = g3*g1*g1, g7 = g5*g1*g1;
                const double rr_dmp = 0.5 * (m_mp_mrad[iat] + m_mp_mrad[jat]);
                const double fdmp3 = 1.0 / (1.0 + 6.0 * std::pow(rr_dmp * g1, mp_dmp3));
                const double fdmp5 = 1.0 / (1.0 + 6.0 * std::pow(rr_dmp * g1, mp_dmp5));
                const double ddmp3 = g5 * fdmp3 * (-3.0 + mp_dmp3 * (1.0 - fdmp3));
                const double ddmp5 = fdmp5 * (-5.0 - mp_dmp5 * (fdmp5 - 1.0));

                const double dpiqj = (vx*m_wfn.dp_at(0,iat)+vy*m_wfn.dp_at(1,iat)+vz*m_wfn.dp_at(2,iat))*m_wfn.q_at(jat);
                const double qidpj = (vx*m_wfn.dp_at(0,jat)+vy*m_wfn.dp_at(1,jat)+vz*m_wfn.dp_at(2,jat))*m_wfn.q_at(iat);
                const double diff_sd = dpiqj - qidpj;
                double gx = -ddmp3*vx*diff_sd + fdmp3*g3*(m_wfn.q_at(iat)*m_wfn.dp_at(0,jat) - m_wfn.q_at(jat)*m_wfn.dp_at(0,iat));
                double gy = -ddmp3*vy*diff_sd + fdmp3*g3*(m_wfn.q_at(iat)*m_wfn.dp_at(1,jat) - m_wfn.q_at(jat)*m_wfn.dp_at(1,iat));
                double gz = -ddmp3*vz*diff_sd + fdmp3*g3*(m_wfn.q_at(iat)*m_wfn.dp_at(2,jat) - m_wfn.q_at(jat)*m_wfn.dp_at(2,iat));
                const double fddr_sd = 3.0*diff_sd*mp_dmp3*fdmp3*g3*(fdmp3/rr_dmp)*std::pow(rr_dmp*g1, mp_dmp3);
                dEdr_mp(iat) += fddr_sd; dEdr_mp(jat) += fddr_sd;

                const double dpidpj = m_wfn.dp_at(0,iat)*m_wfn.dp_at(0,jat)+m_wfn.dp_at(1,iat)*m_wfn.dp_at(1,jat)+m_wfn.dp_at(2,iat)*m_wfn.dp_at(2,jat);
                const double dpiv = vx*m_wfn.dp_at(0,iat)+vy*m_wfn.dp_at(1,iat)+vz*m_wfn.dp_at(2,iat);
                const double dpjv = vx*m_wfn.dp_at(0,jat)+vy*m_wfn.dp_at(1,jat)+vz*m_wfn.dp_at(2,jat);
                const double edd = dpidpj*r2 - 3.0*dpjv*dpiv;
                gx += -2.0*fdmp5*g5*dpidpj*vx + 3.0*fdmp5*g5*(dpiv*m_wfn.dp_at(0,jat)+dpjv*m_wfn.dp_at(0,iat)) - edd*ddmp5*g7*vx;
                gy += -2.0*fdmp5*g5*dpidpj*vy + 3.0*fdmp5*g5*(dpiv*m_wfn.dp_at(1,jat)+dpjv*m_wfn.dp_at(1,iat)) - edd*ddmp5*g7*vy;
                gz += -2.0*fdmp5*g5*dpidpj*vz + 3.0*fdmp5*g5*(dpiv*m_wfn.dp_at(2,jat)+dpjv*m_wfn.dp_at(2,iat)) - edd*ddmp5*g7*vz;
                const double fddr_dd = 3.0*edd*mp_dmp5*fdmp5*g5*(fdmp5/rr_dmp)*std::pow(rr_dmp*g1, mp_dmp5);
                dEdr_mp(iat) += fddr_dd; dEdr_mp(jat) += fddr_dd;

                const double qi = m_wfn.q_at(iat), qj = m_wfn.q_at(jat);
                const double eq =
                      (qj*m_wfn.qp_at(0,iat)+m_wfn.qp_at(0,jat)*qi)*vx*vx
                    + 2.0*(qj*m_wfn.qp_at(1,iat)+m_wfn.qp_at(1,jat)*qi)*vx*vy
                    + (qj*m_wfn.qp_at(2,iat)+m_wfn.qp_at(2,jat)*qi)*vy*vy
                    + 2.0*(qj*m_wfn.qp_at(3,iat)+m_wfn.qp_at(3,jat)*qi)*vx*vz
                    + 2.0*(qj*m_wfn.qp_at(4,iat)+m_wfn.qp_at(4,jat)*qi)*vy*vz
                    + (qj*m_wfn.qp_at(5,iat)+m_wfn.qp_at(5,jat)*qi)*vz*vz;
                const double tivx = m_wfn.qp_at(0,iat)*vx+m_wfn.qp_at(1,iat)*vy+m_wfn.qp_at(3,iat)*vz;
                const double tivy = m_wfn.qp_at(1,iat)*vx+m_wfn.qp_at(2,iat)*vy+m_wfn.qp_at(4,iat)*vz;
                const double tivz = m_wfn.qp_at(3,iat)*vx+m_wfn.qp_at(4,iat)*vy+m_wfn.qp_at(5,iat)*vz;
                const double tjvx = m_wfn.qp_at(0,jat)*vx+m_wfn.qp_at(1,jat)*vy+m_wfn.qp_at(3,jat)*vz;
                const double tjvy = m_wfn.qp_at(1,jat)*vx+m_wfn.qp_at(2,jat)*vy+m_wfn.qp_at(4,jat)*vz;
                const double tjvz = m_wfn.qp_at(3,jat)*vx+m_wfn.qp_at(4,jat)*vy+m_wfn.qp_at(5,jat)*vz;
                gx += -eq*ddmp5*g7*vx - 2.0*fdmp5*g5*(qi*tjvx+qj*tivx);
                gy += -eq*ddmp5*g7*vy - 2.0*fdmp5*g5*(qi*tjvy+qj*tivy);
                gz += -eq*ddmp5*g7*vz - 2.0*fdmp5*g5*(qi*tjvz+qj*tivz);
                const double fddr_sq = eq*3.0*mp_dmp5*fdmp5*g5*(fdmp5/rr_dmp)*std::pow(rr_dmp*g1, mp_dmp5);
                dEdr_mp(iat) += fddr_sq; dEdr_mp(jat) += fddr_sq;

                m_gradient(iat,0) += gx; m_gradient(jat,0) -= gx;
                m_gradient(iat,1) += gy; m_gradient(jat,1) -= gy;
                m_gradient(iat,2) += gz; m_gradient(jat,2) -= gz;
            }
        }
        for (int i = 0; i < nat; ++i) {
            const int zi = m_atoms[i];
            const double rad = p_rad[zi-1], vcn = p_vcn[zi-1];
            const double arg = m_coordination_numbers(i) - vcn - mp_shift;
            const double t1 = std::exp(-mp_kexp * arg);
            const double t2 = (mp_rmax - rad) / (1.0 + t1);
            dEdcn(i) += dEdr_mp(i) * (-t2 * mp_kexp * t1 / (1.0 + t1));
        }
    }

    // 4. CN chain-rule gradient (host; identical to section 4 in calculateGradient).
    std::vector<double> rcov(nat);
    for (int i = 0; i < nat; ++i) rcov[i] = covalent_rad_d3_au(m_atoms[i]);
    for (int i = 0; i < nat; ++i) {
        for (int j = 0; j < i; ++j) {
            const double dx = xyz[3*i+0] - xyz[3*j+0];
            const double dy = xyz[3*i+1] - xyz[3*j+1];
            const double dz = xyz[3*i+2] - xyz[3*j+2];
            const double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 < 1.0e-12) continue;
            const double r  = std::sqrt(r2);
            const double rc = rcov[i] + rcov[j];
            double dcndr;
            if (m_method == MethodType::GFN1) {
                constexpr double k = 16.0;
                const double c = 1.0 / (1.0 + std::exp(-k * (rc / r - 1.0)));
                dcndr = -c * (1.0 - c) * k * rc / r2;
            } else {
                constexpr double ka = 10.0, kb = 20.0, rs = 2.0;
                const double c1 = 1.0 / (1.0 + std::exp(-ka * (rc / r - 1.0)));
                const double c2 = 1.0 / (1.0 + std::exp(-kb * ((rc + rs) / r - 1.0)));
                const double dc1dr = -c1 * (1.0 - c1) * ka * rc / r2;
                const double dc2dr = -c2 * (1.0 - c2) * kb * (rc + rs) / r2;
                dcndr = dc1dr * c2 + c1 * dc2dr;
            }
            const double factor = (dEdcn(i) + dEdcn(j)) * dcndr / r;
            m_gradient(i, 0) += factor * dx;  m_gradient(j, 0) -= factor * dx;
            m_gradient(i, 1) += factor * dy;  m_gradient(j, 1) -= factor * dy;
            m_gradient(i, 2) += factor * dz;  m_gradient(j, 2) -= factor * dz;
        }
    }
    return true;
}

} // namespace curcuma::xtb
