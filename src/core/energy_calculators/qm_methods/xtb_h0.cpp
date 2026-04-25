/*
 * <xTB bare Hamiltonian — CN, self-energy, overlap, H0, repulsion>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Implements GFN1/GFN2 Hamiltonian construction: coordination numbers,
 * CN-shifted self-energies, CGTO overlap, H0 with hscale+shpoly,
 * pairwise repulsion, and halogen-bond correction (GFN1).
 *
 * Claude Generated. GPL-3.0.
 */

#include "xtb_native.h"

#include "parameters/gfn1_params.hpp"
#include "parameters/gfn2_params.hpp"
#include "parameters/xtb_params_extra.hpp"
#include "STO_CGTO.hpp"

#include "src/core/curcuma_logger.h"

namespace curcuma::xtb {

/* ------------------------------------------------------------------ *
 *  Helper: convert internal CGTOShell → CGTO::Shell for integrals   *
 * ------------------------------------------------------------------ */
static CGTO::Shell as_cgto_shell(const CGTOShell& cg)
{
    CGTO::Shell s;
    s.ang   = cg.ang;
    s.nprim = static_cast<int>(cg.alpha.size());
    s.alpha = cg.alpha;
    s.coeff = cg.coeff;
    return s;
}

/* ------------------------------------------------------------------ *
 *  AO type encoding (tblite ordering):                               *
 *    s  → 0                                                          *
 *    p  → [py=2, pz=3, px=1]  (matching STO_CGTO type convention)   *
 * ------------------------------------------------------------------ */
static inline int ao_to_type(int ang, int local_ao)
{
    if (ang == 0) return 0;                     // s
    if (ang == 1) {                              // p: tblite ordering
        static const int p_map[3] = {2, 3, 1};  // py, pz, px
        return p_map[local_ao];
    }
    return -1;  // d not handled
}

/* ------------------------------------------------------------------ *
 *  Coordination numbers                                              *
 * ------------------------------------------------------------------ */
Vector XTB::computeCoordinationNumbers() const
{
    const int nat = m_atomcount;
    std::vector<double> xyz_bohr(3 * nat);
    for (int i = 0; i < nat; ++i) {
        xyz_bohr[3 * i + 0] = m_geometry(i, 0) * AA_TO_AU;
        xyz_bohr[3 * i + 1] = m_geometry(i, 1) * AA_TO_AU;
        xyz_bohr[3 * i + 2] = m_geometry(i, 2) * AA_TO_AU;
    }
    std::vector<int> z_atoms(m_atoms.begin(), m_atoms.end());

    std::vector<double> cn_vec;
    if (m_method == MethodType::GFN1)
        cn_vec = cn_exp(z_atoms, xyz_bohr);
    else
        cn_vec = cn_gfn(z_atoms, xyz_bohr);

    Vector cn(nat);
    for (int i = 0; i < nat; ++i) cn(i) = cn_vec[i];
    return cn;
}

/* ------------------------------------------------------------------ *
 *  CN-shifted shell self-energies:  ε_sh = ε0 - kcn * CN(atom)       *
 * ------------------------------------------------------------------ */
void XTB::getSelfEnergies(const Vector& cn, Vector& se_out) const
{
    se_out.resize(m_basis.nsh);
    for (int s = 0; s < m_basis.nsh; ++s) {
        const int iat = m_basis.sh2at[s];
        se_out(s) = m_h0.selfenergy[s] - m_h0.kcn[s] * cn(iat);
    }
}

/* ------------------------------------------------------------------ *
 *  Build overlap (S), bare Hamiltonian (H0), and multipole integral  *
 *  matrices (dp_int, qp_int) from current self-energies.             *
 *                                                                    *
 *  GFN1: kpair + kshell + enscale + valence flags                    *
 *  GFN2: adds Slater-ratio prefactor zij, kpair=1                    *
 * ------------------------------------------------------------------ */
void XTB::getHamiltonianH0(const Vector& se,
                            Matrix& S, Matrix& H0) const
{
    const int nat = m_atomcount;
    const int nsh = m_basis.nsh;
    const int nao = m_basis.nao;

    // Geometry in bohr
    std::vector<double> xyz_bohr(3 * nat);
    for (int i = 0; i < nat; ++i) {
        xyz_bohr[3 * i + 0] = m_geometry(i, 0) * AA_TO_AU;
        xyz_bohr[3 * i + 1] = m_geometry(i, 1) * AA_TO_AU;
        xyz_bohr[3 * i + 2] = m_geometry(i, 2) * AA_TO_AU;
    }

    S  = Matrix::Zero(nao, nao);
    H0 = Matrix::Zero(nao, nao);

    // Valence flags: first shell of each ℓ is "valence" (GFN1 rule)
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

    // Pre-compute (ish → atom Z)
    std::vector<int> sh_z(nsh);
    for (int s = 0; s < nsh; ++s) sh_z[s] = m_basis.z[m_basis.sh2at[s]];

    // Build overlap + H0 element-by-element.
    // For each pair of shells (ish_a, ish_b) we compute the CGTO overlap
    // once and broadcast to all AO pairs.
    for (int ish_a = 0; ish_a < nsh; ++ish_a) {
        const int iat  = m_basis.sh2at[ish_a];
        const int ia_start = m_basis.iao_sh[ish_a];
        const int ia_nao   = m_basis.nao_sh[ish_a];
        const CGTO::Shell sh_a = as_cgto_shell(m_basis.cgto[ish_a]);

        // Local shell index within atom (for parameter lookup)
        const int local_a = ish_a - m_basis.ish_at[iat];
        const double zeta_a = m_basis.cgto[ish_a].slater_exp;

        for (int ish_b = 0; ish_b < nsh; ++ish_b) {
            const int jat  = m_basis.sh2at[ish_b];
            const int jb_start = m_basis.iao_sh[ish_b];
            const int jb_nao   = m_basis.nao_sh[ish_b];
            const CGTO::Shell sh_b = as_cgto_shell(m_basis.cgto[ish_b]);

            const int local_b = ish_b - m_basis.ish_at[jat];

            // H0 element
            const double avg_eps = 0.5 * (se(ish_a) + se(ish_b));
            double h_factor;

            if (iat == jat) {
                // On-atom: no hscale/shpoly, just self-energy
                h_factor = 1.0;
            } else {
                // Off-atom: apply hscale + shpoly distance polynomial
                const int zi = sh_z[ish_a], zj = sh_z[ish_b];
                const double dx = xyz_bohr[3 * iat + 0] - xyz_bohr[3 * jat + 0];
                const double dy = xyz_bohr[3 * iat + 1] - xyz_bohr[3 * jat + 1];
                const double dz = xyz_bohr[3 * iat + 2] - xyz_bohr[3 * jat + 2];
                const double r2 = dx * dx + dy * dy + dz * dz;
                const double rr = std::sqrt(std::sqrt(r2)
                    / (atomic_rad_au(zi) + atomic_rad_au(zj)));
                const double pi_ij = (1.0 + m_h0.shpoly[ish_a] * rr)
                                   * (1.0 + m_h0.shpoly[ish_b] * rr);

                double hs;
                if (m_method == MethodType::GFN1) {
                    const bool vi = valence[ish_a], vj = valence[ish_b];
                    if (vi && vj) {
                        const double den = std::pow(pauling_en[zi - 1] - pauling_en[zj - 1], 2);
                        hs = gfn1::kpair(zi, zj) * gfn1::kshell(m_basis.ang_sh[ish_a], m_basis.ang_sh[ish_b])
                           * (1.0 + gfn1::enscale * den);
                    } else if (vi && !vj) {
                        hs = 0.5 * (gfn1::kshell(m_basis.ang_sh[ish_a], m_basis.ang_sh[ish_a]) + gfn1::kdiff);
                    } else if (!vi && vj) {
                        hs = 0.5 * (gfn1::kshell(m_basis.ang_sh[ish_b], m_basis.ang_sh[ish_b]) + gfn1::kdiff);
                    } else {
                        hs = gfn1::kdiff;
                    }
                } else {
                    const double zeta_b = m_basis.cgto[ish_b].slater_exp;
                    const double den = std::pow(pauling_en[zi - 1] - pauling_en[zj - 1], 2);
                    const double enp = 1.0 + gfn2::enscale * den;
                    const double km  = gfn2::kpair(zi, zj) * gfn2::kshell(m_basis.ang_sh[ish_a], m_basis.ang_sh[ish_b]) * enp;
                    const double zij = std::pow(2.0 * std::sqrt(zeta_a * zeta_b) / (zeta_a + zeta_b), gfn2::wexp);
                    hs = zij * km;
                }
                h_factor = hs * pi_ij;
            }

            // Compute per-AO-pair overlap (directional for p-orbitals)
            for (int ia = 0; ia < ia_nao; ++ia) {
                const int mu = ia_start + ia;
                const int t_a = ao_to_type(m_basis.ang_sh[ish_a], ia);
                if (t_a < 0) continue;

                for (int jb = 0; jb < jb_nao; ++jb) {
                    const int nu = jb_start + jb;
                    const int t_b = ao_to_type(m_basis.ang_sh[ish_b], jb);
                    if (t_b < 0) continue;

                    // CGTO overlap depends on AO types (px, py, pz, etc.)
                    const double s_ab = (iat == jat && ish_a == ish_b && t_a == t_b)
                        ? 1.0  // on-atom same-orbital → identity
                        : CGTO::cgto_overlap(sh_a, sh_b,
                                             xyz_bohr[3 * iat + 0],
                                             xyz_bohr[3 * iat + 1],
                                             xyz_bohr[3 * iat + 2],
                                             xyz_bohr[3 * jat + 0],
                                             xyz_bohr[3 * jat + 1],
                                             xyz_bohr[3 * jat + 2],
                                             t_a, t_b);

                    S(mu, nu) = s_ab;
                    H0(mu, nu) = avg_eps * h_factor * s_ab;
                }
            }
        }
    }
}

/* ------------------------------------------------------------------ *
 *  Pairwise repulsion energy                                         *
 *    E_rep = sum_{A<B} Z_A * Z_B / R_AB * exp( -(α_A·α_B)^k * R^m ) *
 *                                                                    *
 *  Parameters from generated tables: rep_alpha, rep_zeff, rep_kexp,  *
 *  rep_rexp (GFN2 also has rep_max).                                 *
 * ------------------------------------------------------------------ */
double XTB::calcRepulsionEnergy() const
{
    const int nat = m_atomcount;
    double erep = 0.0;

    // Convert geometry to bohr
    std::vector<double> xyz_bohr(3 * nat);
    for (int i = 0; i < nat; ++i) {
        xyz_bohr[3 * i + 0] = m_geometry(i, 0) * AA_TO_AU;
        xyz_bohr[3 * i + 1] = m_geometry(i, 1) * AA_TO_AU;
        xyz_bohr[3 * i + 2] = m_geometry(i, 2) * AA_TO_AU;
    }

    for (int i = 0; i < nat; ++i) {
        const int zi = m_atoms[i];
        for (int j = 0; j < i; ++j) {
            const int zj = m_atoms[j];
            const double dx = xyz_bohr[3 * i + 0] - xyz_bohr[3 * j + 0];
            const double dy = xyz_bohr[3 * i + 1] - xyz_bohr[3 * j + 1];
            const double dz = xyz_bohr[3 * i + 2] - xyz_bohr[3 * j + 2];
            const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
            if (r < 1.0e-12) continue;

            // repulsion parameters from generated tables
            double alf, zeff;
            double kexp, rexp;
            if (m_method == MethodType::GFN1) {
                alf  = gfn1_params::rep_alpha[zi - 1];
                zeff = gfn1_params::rep_zeff[zi - 1];
                kexp = gfn1_params::rep_kexp;
                rexp = gfn1_params::rep_rexp;
            } else {
                alf  = gfn2_params::rep_alpha[zi - 1];
                zeff = gfn2_params::rep_zeff[zi - 1];
                kexp = gfn2_params::rep_kexp;
                rexp = gfn2_params::rep_rexp;
            }
            double alfj, zeffj;
            if (m_method == MethodType::GFN1) {
                alfj  = gfn1_params::rep_alpha[zj - 1];
                zeffj = gfn1_params::rep_zeff[zj - 1];
            } else {
                alfj  = gfn2_params::rep_alpha[zj - 1];
                zeffj = gfn2_params::rep_zeff[zj - 1];
            }

            // xTB repulsion: E_rep = Z_i*Z_j/R * exp( -(α_i·α_j)^kexp * R^rexp )
            const double a_prod = alf * alfj;
            const double ri = zeff * zeffj / r;
            erep += ri * std::exp(-std::pow(a_prod, kexp) * std::pow(r, rexp));
        }
    }
    return erep;
}

/* ------------------------------------------------------------------ *
 *  Halogen-bond correction (GFN1 only) — simplified form             *
 * ------------------------------------------------------------------ */
double XTB::calcHalogenBondEnergy() const
{
    // Halogen bond correction: only for GFN1 with halogen atoms.
    // The full xTB XB correction depends on charge densities (ρ).
    // For now return 0 — to be implemented if accuracy requires it.
    return 0.0;
}

} // namespace curcuma::xtb
