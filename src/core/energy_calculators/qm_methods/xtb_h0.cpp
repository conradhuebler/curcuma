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
#include "xtb_ao_utils.hpp"
#include "xtb_coulomb.hpp"

#include "src/core/curcuma_logger.h"

#include <limits>

namespace curcuma::xtb {

/* as_cgto_shell() and ao_to_type() now live in xtb_ao_utils.hpp (X-I5). */

/* ------------------------------------------------------------------ *
 *  Coordination numbers                                              *
 * ------------------------------------------------------------------ */
Vector XTB::computeCoordinationNumbers() const
{
    const int nat = m_atomcount;
    // X-I3 (Claude Generated): the H0/D4 coordination numbers depend only on the
    // current geometry but are requested up to 3× per Calculation (pre-SCF
    // self-energies, the GFN2 CPSCF gradient response, evaluateComponentsAtFixedDensity).
    // Memoise; the cache is invalidated on every geometry change (UpdateMolecule /
    // InitialiseMolecule), so the returned CN is bit-identical to recomputing.
    if (m_cn_cache_valid && m_cn_cache.size() == nat)
        return m_cn_cache;

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
    m_cn_cache = cn;            // X-I3: memoise for the rest of this geometry
    m_cn_cache_valid = true;
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
 *  Flatten the built basis + H0 parameters for the GPU integral build *
 *  (Stage 3, Claude Generated). Plain std::vector bundles only — no   *
 *  CUDA types — so the GPU wrapper can upload them and the core stays  *
 *  CUDA-free. The primitives are already post-orthogonalised on the   *
 *  host (CGTO::orthogonalize); we flatten them verbatim.              *
 * ------------------------------------------------------------------ */
void XTB::exportGpuBasis(GpuBasisFlat& bf, GpuH0Flat& hf) const
{
    const int nat = m_basis.nat;
    const int nsh = m_basis.nsh;
    const int nao = m_basis.nao;

    bf.nat = nat;
    bf.nsh = nsh;
    bf.nao = nao;

    bf.z = m_basis.z;
    bf.xyz_bohr.resize(3 * nat);
    for (int i = 0; i < nat; ++i) {
        bf.xyz_bohr[3 * i + 0] = m_geometry(i, 0) * AA_TO_AU;
        bf.xyz_bohr[3 * i + 1] = m_geometry(i, 1) * AA_TO_AU;
        bf.xyz_bohr[3 * i + 2] = m_geometry(i, 2) * AA_TO_AU;
    }

    bf.sh2at  = m_basis.sh2at;
    bf.ang_sh = m_basis.ang_sh;
    bf.iao_sh = m_basis.iao_sh;
    bf.nao_sh = m_basis.nao_sh;
    bf.ish_at = m_basis.ish_at;
    bf.nsh_at = m_basis.nsh_at;
    bf.ao2at  = m_basis.ao2at;
    bf.ao2sh  = m_basis.ao2sh;

    // Flatten the per-shell CGTO primitives with an exclusive-prefix-sum offset.
    bf.sh_nprim.resize(nsh);
    bf.sh_prim_off.resize(nsh);
    bf.sh_zeta.resize(nsh);
    int off = 0;
    for (int s = 0; s < nsh; ++s) {
        const int nprim = static_cast<int>(m_basis.cgto[s].alpha.size());
        bf.sh_nprim[s]    = nprim;
        bf.sh_prim_off[s] = off;
        bf.sh_zeta[s]     = m_basis.cgto[s].slater_exp;
        off += nprim;
    }
    bf.prim_alpha.resize(off);
    bf.prim_coeff.resize(off);
    for (int s = 0; s < nsh; ++s) {
        const int base = bf.sh_prim_off[s];
        const auto& cg = m_basis.cgto[s];
        for (int p = 0; p < bf.sh_nprim[s]; ++p) {
            bf.prim_alpha[base + p] = cg.alpha[p];
            bf.prim_coeff[base + p] = cg.coeff[p];
        }
    }

    // GFN1 valence flags: first shell of each ℓ per atom is "valence"
    // (mirrors getHamiltonianH0, xtb_h0.cpp:113-124). Zeroed for GFN2 (unused).
    bf.is_gfn2 = (m_method == MethodType::GFN2) ? 1 : 0;
    bf.valence.assign(nsh, 0);
    if (m_method == MethodType::GFN1) {
        for (int iat = 0; iat < nat; ++iat) {
            bool ang_seen[3] = {false, false, false};
            for (int ish = 0; ish < m_basis.nsh_at[iat]; ++ish) {
                const int sh = m_basis.ish_at[iat] + ish;
                const int l  = m_basis.ang_sh[sh];
                if (l >= 0 && l < 3 && !ang_seen[l]) { bf.valence[sh] = 1; ang_seen[l] = true; }
            }
        }
    }

    // Per-shell Coulomb hardness (molecule-constant): hubbard_parameter(Z) ×
    // shell_hubbard(ang, Z). Precomputed here so the device k_gamma needs no
    // element tables (mirrors coulomb::shell_hardness, xtb_coulomb.hpp:45).
    const coulomb::Method cm = (m_method == MethodType::GFN2)
        ? coulomb::Method::GFN2 : coulomb::Method::GFN1;
    bf.shell_hardness.resize(nsh);
    for (int s = 0; s < nsh; ++s)
        bf.shell_hardness[s] = coulomb::shell_hardness(cm, m_basis.z[m_basis.sh2at[s]],
                                                       m_basis.ang_sh[s]);

    // Per-atom repulsion parameters (molecule-constant) for the Stage-4 gradient.
    bf.rep_alpha.resize(nat);
    bf.rep_zeff.resize(nat);
    for (int i = 0; i < nat; ++i) {
        const int z = m_basis.z[i];
        if (m_method == MethodType::GFN1) {
            bf.rep_alpha[i] = gfn1_params::rep_alpha[z - 1];
            bf.rep_zeff[i]  = gfn1_params::rep_zeff[z - 1];
        } else {
            bf.rep_alpha[i] = gfn2_params::rep_alpha[z - 1];
            bf.rep_zeff[i]  = gfn2_params::rep_zeff[z - 1];
        }
    }

    hf.selfenergy = m_h0.selfenergy;
    hf.kcn        = m_h0.kcn;
    hf.shpoly     = m_h0.shpoly;
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
    //
    // Parallel over the outer shell ish_a (Claude Generated): each ish_a writes
    // only its own AO rows of S/H0, so the stripes touch disjoint matrix blocks —
    // no locking, bit-identical to the serial result. effectiveIntraThreads() keeps
    // this serial unless a single large molecule was granted a thread budget.
    // B1 (Jul 2026): as_cgto_shell() copies two std::vectors; called in the inner
    // ish_b loop it ran nsh^2 = 116k times on complex/231. Shells are
    // geometry-independent, so convert once and index (identical values).
    std::vector<CGTO::Shell> shells(nsh);
    for (int ish = 0; ish < nsh; ++ish)
        shells[ish] = as_cgto_shell(m_basis.cgto[ish]);

    const int n_threads = effectiveIntraThreads(nsh);
    parallelStripes(n_threads, [&](int tid, int nth) {
    for (int ish_a = tid; ish_a < nsh; ish_a += nth) {
        const int iat  = m_basis.sh2at[ish_a];
        const int ia_start = m_basis.iao_sh[ish_a];
        const int ia_nao   = m_basis.nao_sh[ish_a];
        const CGTO::Shell& sh_a = shells[ish_a];

        // Local shell index within atom (for parameter lookup)
        const int local_a = ish_a - m_basis.ish_at[iat];
        const double zeta_a = m_basis.cgto[ish_a].slater_exp;

        for (int ish_b = 0; ish_b < nsh; ++ish_b) {
            const int jat  = m_basis.sh2at[ish_b];
            const int jb_start = m_basis.iao_sh[ish_b];
            const int jb_nao   = m_basis.nao_sh[ish_b];
            const CGTO::Shell& sh_b = shells[ish_b];

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

            const int ang_a = m_basis.ang_sh[ish_a];
            const int ang_b = m_basis.ang_sh[ish_b];

            if (ang_a < 2 && ang_b < 2) {
                // ---- s/p path: B2 (Jul 2026) computes the whole component block
                // from one pass over primitive pairs instead of re-running the
                // primitive loop (and its pow/exp) once per AO component. ----
                double blk_s[9];
                CGTO::cgto_overlap_block(sh_a, ang_a, sh_b, ang_b,
                                         xyz_bohr[3 * iat + 0],
                                         xyz_bohr[3 * iat + 1],
                                         xyz_bohr[3 * iat + 2],
                                         xyz_bohr[3 * jat + 0],
                                         xyz_bohr[3 * jat + 1],
                                         xyz_bohr[3 * jat + 2],
                                         blk_s);
                for (int ia = 0; ia < ia_nao; ++ia) {
                    const int mu = ia_start + ia;
                    const int t_a = ao_to_type(ang_a, ia);
                    if (t_a < 0) continue;

                    for (int jb = 0; jb < jb_nao; ++jb) {
                        const int nu = jb_start + jb;
                        const int t_b = ao_to_type(ang_b, jb);
                        if (t_b < 0) continue;

                        const double s_ab = (iat == jat && ish_a == ish_b && t_a == t_b)
                            ? 1.0  // on-atom same-orbital → identity
                            : blk_s[ia * jb_nao + jb];

                        S(mu, nu) = s_ab;
                        H0(mu, nu) = avg_eps * h_factor * s_ab;
                    }
                }
            } else {
                // ---- X-I1: d-touching shell pair — cartesian(6)->spherical(5)
                // block via dtrafo. Local AO index maps to spherical m (tblite
                // order [-l..+l]), matching the s/p ordering above.
                double blk[6 * 6];
                sphericalOverlapBlock(sh_a, ang_a, sh_b, ang_b,
                                      xyz_bohr[3 * iat + 0],
                                      xyz_bohr[3 * iat + 1],
                                      xyz_bohr[3 * iat + 2],
                                      xyz_bohr[3 * jat + 0],
                                      xyz_bohr[3 * jat + 1],
                                      xyz_bohr[3 * jat + 2],
                                      blk, 6);
                for (int ia = 0; ia < ia_nao; ++ia) {
                    const int mu = ia_start + ia;
                    for (int jb = 0; jb < jb_nao; ++jb) {
                        const int nu = jb_start + jb;
                        double s_ab = blk[ia * 6 + jb];
                        if (iat == jat && ish_a == ish_b && ia == jb)
                            s_ab = 1.0;  // exact on-atom diagonal (orthonormal AOs)
                        S(mu, nu) = s_ab;
                        H0(mu, nu) = avg_eps * h_factor * s_ab;
                    }
                }
            }
        }
    }
    });  // parallelStripes over ish_a
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

            // xTB repulsion (TBLite effective.f90):
            //   alpha_ij = sqrt(alpha_i * alpha_j)
            //   kexp_ij  = kexp_light  if both Z<=2 (GFN2), else kexp
            //   E_rep    = Z_i*Z_j / R^rexp * exp(-alpha_ij * R^kexp_ij)
            const double alpha_pair = std::sqrt(alf * alfj);
            double kexp_pair = kexp;
            if (m_method == MethodType::GFN2) {
                if (zi <= 2 && zj <= 2)
                    kexp_pair = gfn2_params::rep_kexp_light;
            }
            const double r_kexp = std::pow(r, kexp_pair);
            erep += (zeff * zeffj / std::pow(r, rexp))
                  * std::exp(-alpha_pair * r_kexp);
        }
    }
    return erep;
}

/* ------------------------------------------------------------------ *
 *  Halogen-bond correction (GFN1 only)                               *
 *                                                                    *
 *  Classical three-body term for the B–X···A halogen bond: X is the  *
 *  halogen donor (Cl, Br, I, At), A a Lewis-base acceptor            *
 *  (N, O, P, S) and B the atom nearest to X (its covalent partner).  *
 *  It is a purely geometric correction — no density enters.          *
 *                                                                    *
 *    E_XB = sum_{X,A} f_ang * c_X * (t^12 - k_damp*t^6) / (1 + t^12) *
 *    t     = r0 / R_AX,   r0 = radScale * ( rad_X + rad_A )          *
 *    f_ang = ( 1/2 - 1/4 * cos(theta_BXA) )^6                        *
 *                                                                    *
 *  GFN1-xTB: Grimme, Bannwarth, Shushkov, JCTC 13 (2017) 1989.       *
 *  Ports external/xtb/src/xtb/halogen.f90 (xbpot) with the pair list *
 *  of scf_module.F90:370-404; external/tblite/src/tblite/classical/  *
 *  halogen.f90 is algebraically identical (checked line by line).    *
 *  c_X is zero for Cl in GFN1, so only Br, I and At contribute.      *
 *  Claude Generated.                                                 *
 * ------------------------------------------------------------------ */
namespace {

/// Halogen-bond donor (xtb scf_module.F90::xbond — F is deliberately excluded).
inline bool xb_is_halogen(int z) { return z == 17 || z == 35 || z == 53 || z == 85; }
/// Halogen-bond acceptor (Lewis base) — same source.
inline bool xb_is_acceptor(int z) { return z == 7 || z == 8 || z == 15 || z == 16; }

constexpr double xb_alp    = 6.0;   ///< exponent of the angular damping function
constexpr double xb_lj     = 12.0;  ///< LJ exponent (xtb: ljexp, GFN1 = 12)
constexpr double xb_lj2    = 0.5 * xb_lj;
constexpr double xb_cutoff = 20.0;  ///< Bohr; xtb tests sqrab < 400

/// One B–X···A triple: {halogen X, acceptor A, nearest neighbour B of X}.
struct XBTriple { int x, a, b; };

/// Build the halogen-bond list. B is the atom closest to X over the whole
/// molecule (xtb searches per pair, but the result depends on X alone).
std::vector<XBTriple> buildHalogenBondList(const std::vector<int>& z,
                                           const std::vector<double>& xyz)
{
    const int nat = static_cast<int>(z.size());
    std::vector<XBTriple> list;
    if (nat < 2) return list;

    // Nearest neighbour of every halogen (-1 = none / not a halogen).
    std::vector<int> nearest(nat, -1);
    for (int x = 0; x < nat; ++x) {
        if (!xb_is_halogen(z[x])) continue;
        double best = std::numeric_limits<double>::max();
        for (int m = 0; m < nat; ++m) {
            if (m == x) continue;
            const double dx = xyz[3*m+0] - xyz[3*x+0];
            const double dy = xyz[3*m+1] - xyz[3*x+1];
            const double dz = xyz[3*m+2] - xyz[3*x+2];
            const double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 < best) { best = r2; nearest[x] = m; }
        }
    }

    const double cut2 = xb_cutoff * xb_cutoff;
    for (int x = 0; x < nat; ++x) {
        if (nearest[x] < 0) continue;
        for (int a = 0; a < nat; ++a) {
            if (!xb_is_acceptor(z[a])) continue;
            const double dx = xyz[3*x+0] - xyz[3*a+0];
            const double dy = xyz[3*x+1] - xyz[3*a+1];
            const double dz = xyz[3*x+2] - xyz[3*a+2];
            if (dx*dx + dy*dy + dz*dz >= cut2) continue;
            list.push_back({x, a, nearest[x]});
        }
    }
    return list;
}

} // anonymous namespace

double XTB::calcHalogenBondEnergy() const
{
    if (m_method != MethodType::GFN1) return 0.0;   // GFN2 has no XB correction

    const int nat = m_atomcount;
    if (nat < 2) return 0.0;

    std::vector<int> z(m_atoms.begin(), m_atoms.end());
    std::vector<double> xyz(3 * nat);
    for (int i = 0; i < nat; ++i) {
        xyz[3*i+0] = m_geometry(i, 0) * AA_TO_AU;
        xyz[3*i+1] = m_geometry(i, 1) * AA_TO_AU;
        xyz[3*i+2] = m_geometry(i, 2) * AA_TO_AU;
    }

    const auto list = buildHalogenBondList(z, xyz);
    double exb = 0.0;

    for (const auto& t : list) {
        const double cc = gfn1_params::halogen_bond[z[t.x] - 1];
        if (cc == 0.0) continue;                    // Cl in GFN1
        const double r0ax = gfn1_params::halogen_radscale
                          * (atomic_rad_au(z[t.x]) + atomic_rad_au(z[t.a]));

        // dxa = A - X, dxb = B - X, dba = A - B
        double dxa[3], dxb[3], dba[3];
        for (int k = 0; k < 3; ++k) {
            dxa[k] = xyz[3*t.a+k] - xyz[3*t.x+k];
            dxb[k] = xyz[3*t.b+k] - xyz[3*t.x+k];
            dba[k] = xyz[3*t.a+k] - xyz[3*t.b+k];
        }
        const double d2ax = dxa[0]*dxa[0] + dxa[1]*dxa[1] + dxa[2]*dxa[2];
        const double d2bx = dxb[0]*dxb[0] + dxb[1]*dxb[1] + dxb[2]*dxb[2];
        const double d2ab = dba[0]*dba[0] + dba[1]*dba[1] + dba[2]*dba[2];
        const double rax  = std::sqrt(d2ax);

        // Angular part: term = cos(angle B-X-A) via the law of cosines.
        const double xy   = std::sqrt(d2bx * d2ax);
        const double term = (d2bx + d2ax - d2ab) / xy;
        const double aterm = std::pow(0.5 - 0.25 * term, xb_alp);

        const double t13 = r0ax / rax;
        const double t14 = std::pow(t13, xb_lj);
        exb += aterm * cc * (t14 - gfn1_params::halogen_damping * std::pow(t13, xb_lj2))
             / (1.0 + t14);
    }
    return exb;
}

/* ------------------------------------------------------------------ *
 *  Analytic gradient of the halogen-bond correction (GFN1 only).     *
 *  Verbatim port of the derivative block in xtb's xbpot; adds into   *
 *  `gradient` (nat×3, Eh/Bohr). Claude Generated.                    *
 * ------------------------------------------------------------------ */
void XTB::addHalogenBondGradient(Matrix& gradient) const
{
    if (m_method != MethodType::GFN1) return;

    const int nat = m_atomcount;
    if (nat < 2 || gradient.rows() != nat || gradient.cols() != 3) return;

    std::vector<int> z(m_atoms.begin(), m_atoms.end());
    std::vector<double> xyz(3 * nat);
    for (int i = 0; i < nat; ++i) {
        xyz[3*i+0] = m_geometry(i, 0) * AA_TO_AU;
        xyz[3*i+1] = m_geometry(i, 1) * AA_TO_AU;
        xyz[3*i+2] = m_geometry(i, 2) * AA_TO_AU;
    }

    const auto list = buildHalogenBondList(z, xyz);

    for (const auto& t : list) {
        const double cc = gfn1_params::halogen_bond[z[t.x] - 1];
        if (cc == 0.0) continue;
        const double damping = gfn1_params::halogen_damping;
        const double r0ax = gfn1_params::halogen_radscale
                          * (atomic_rad_au(z[t.x]) + atomic_rad_au(z[t.a]));

        double dxa[3], dxb[3], dba[3];
        for (int k = 0; k < 3; ++k) {
            dxa[k] = xyz[3*t.a+k] - xyz[3*t.x+k];
            dxb[k] = xyz[3*t.b+k] - xyz[3*t.x+k];
            dba[k] = xyz[3*t.a+k] - xyz[3*t.b+k];
        }
        const double d2ax = dxa[0]*dxa[0] + dxa[1]*dxa[1] + dxa[2]*dxa[2];
        const double d2bx = dxb[0]*dxb[0] + dxb[1]*dxb[1] + dxb[2]*dxb[2];
        const double d2ab = dba[0]*dba[0] + dba[1]*dba[1] + dba[2]*dba[2];
        const double rax  = std::sqrt(d2ax) + 1.0e-18;
        const double rbx  = std::sqrt(d2bx) + 1.0e-18;

        const double xy    = std::sqrt(d2bx * d2ax);
        const double term  = (d2bx + d2ax - d2ab) / xy;
        const double aterm = std::pow(0.5 - 0.25 * term, xb_alp);

        // Damped LJ in terms of t = (r0/r)^6, so the potential is (t²-k·t)/(1+t²).
        const double t14         = std::pow(r0ax / rax, xb_lj2);
        const double numerator   = t14 * t14 - damping * t14;
        const double denominator = 1.0 + t14 * t14;
        const double termlj      = numerator / denominator;

        // LJ radial derivative (denominator + numerator part), scaled by f_ang·c_X.
        double dtermlj = 2.0 * xb_lj2 * numerator * t14 * t14
                       / (rax * denominator * denominator);
        dtermlj += xb_lj2 * t14 * (damping - 2.0 * t14) / (rax * denominator);
        dtermlj *= aterm * cc / rax;
        for (int k = 0; k < 3; ++k) {
            gradient(t.a, k) += dtermlj * dxa[k];
            gradient(t.x, k) -= dtermlj * dxa[k];
        }

        // Derivative of the angular damping function.
        double prefactor = -0.25 * xb_alp * std::pow(0.5 - 0.25 * term, xb_alp - 1.0);
        prefactor *= cc * termlj;

        double dcosterm = (2.0 / rbx - term / rax) * prefactor / rax;   // A–X part
        for (int k = 0; k < 3; ++k) {
            gradient(t.a, k) += dcosterm * dxa[k];
            gradient(t.x, k) -= dcosterm * dxa[k];
        }
        dcosterm = (2.0 / rax - term / rbx) * prefactor / rbx;          // B–X part
        for (int k = 0; k < 3; ++k) {
            gradient(t.b, k) += dcosterm * dxb[k];
            gradient(t.x, k) -= dcosterm * dxb[k];
        }
        const double t13 = 2.0 * prefactor / xy;                        // A–B part
        for (int k = 0; k < 3; ++k) {
            gradient(t.a, k) -= t13 * dba[k];
            gradient(t.b, k) += t13 * dba[k];
        }
    }
}

} // namespace curcuma::xtb
