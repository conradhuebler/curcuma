/*
 * <xTB GFN2 multipole interactions — setup, potential, energy>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Ported from test_xtb_scf_snapshot.cpp (tblite-validated kernels).
 *
 * Provides the three public methods declared in xtb_native.h:
 *   setupMultipole()       — one-time geometry setup
 *   addMultipolePotential()— per-SCF-iteration potential assembly
 *   energyMultipole()      — multipole contribution to total energy
 *
 * Reference: coulomb/multipole.f90 + scf/potential.f90 in tblite.
 *
 * Claude Generated. GPL-3.0.
 */

#include "xtb_native.h"
#include "xtb_multipole_ints.hpp"
#include "parameters/gfn2_params.hpp"
#include "parameters/xtb_params_extra.hpp"

#include <Eigen/Dense>

#include <cmath>

namespace curcuma::xtb {
namespace MI = multipole_ints;

/* ------------------------------------------------------------------ *
 *  AO type encoding (same as xtb_h0.cpp):                            *
 *    s → 0,  p → [py=2, pz=3, px=1]                                 *
 * ------------------------------------------------------------------ */
static inline int ao_to_type(int ang, int local_ao)
{
    if (ang == 0) return 0;
    if (ang == 1) {
        static const int p_map[3] = {2, 3, 1};
        return p_map[local_ao];
    }
    return -1;
}

/* ------------------------------------------------------------------ *
 *  Convert internal CGTOShell → CGTO::Shell for integrals            *
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
 *  setupMultipole()                                                  *
 *                                                                    *
 *  Builds all multipole integrals and interaction matrices that      *
 *  depend only on geometry (not on density). Called once after       *
 *  buildBasis()/buildH0Data() for GFN2.                              *
 *                                                                    *
 *  1. AO-level dipole + quadrupole integrals with origin shift       *
 *     and traceless transform (→ m_dp_int, m_qp_int)                 *
 *  2. CN-dependent damping radii mrad                                *
 *  3. Interaction matrices amat_sd, amat_dd, amat_sq                 *
 * ------------------------------------------------------------------ */
void XTB::setupMultipole()
{
    const int nat = m_atomcount;
    const int nsh = m_basis.nsh;
    const int nao = m_basis.nao;

    if (nat == 0 || nao == 0) return;

    // Geometry in bohr
    std::vector<double> xyz_bohr(3 * nat);
    for (int i = 0; i < nat; ++i) {
        xyz_bohr[3 * i + 0] = m_geometry(i, 0) * AA_TO_AU;
        xyz_bohr[3 * i + 1] = m_geometry(i, 1) * AA_TO_AU;
        xyz_bohr[3 * i + 2] = m_geometry(i, 2) * AA_TO_AU;
    }

    // ---- 1. Global-origin raw dipole + raw-Cartesian quadrupole ----
    std::array<Eigen::MatrixXd, 3> dp_global;
    std::array<Eigen::MatrixXd, 6> qp_global_raw;
    for (int k = 0; k < 3; ++k) dp_global[k]     = Eigen::MatrixXd::Zero(nao, nao);
    for (int k = 0; k < 6; ++k) qp_global_raw[k] = Eigen::MatrixXd::Zero(nao, nao);

    for (int mu = 0; mu < nao; ++mu) {
        const int ish_a = m_basis.ao2sh[mu];
        const int iat   = m_basis.ao2at[mu];
        const int local_a = mu - m_basis.iao_sh[ish_a];
        const int t_a = ao_to_type(m_basis.ang_sh[ish_a], local_a);
        const CGTO::Shell sh_a = as_cgto_shell(m_basis.cgto[ish_a]);

        for (int nu = 0; nu < nao; ++nu) {
            const int ish_b = m_basis.ao2sh[nu];
            const int jat   = m_basis.ao2at[nu];
            const int local_b = nu - m_basis.iao_sh[ish_b];
            const int t_b = ao_to_type(m_basis.ang_sh[ish_b], local_b);
            if (t_a < 0 || t_b < 0) continue;
            const CGTO::Shell sh_b = as_cgto_shell(m_basis.cgto[ish_b]);

            double Sx, D[3], Q[6];
            MI::cgto_multipole(sh_a, sh_b,
                               xyz_bohr[3*iat+0], xyz_bohr[3*iat+1], xyz_bohr[3*iat+2],
                               xyz_bohr[3*jat+0], xyz_bohr[3*jat+1], xyz_bohr[3*jat+2],
                               t_a, t_b, Sx, D, Q);
            for (int k = 0; k < 3; ++k) dp_global[k](mu, nu)     = D[k];
            for (int k = 0; k < 6; ++k) qp_global_raw[k](mu, nu) = Q[k];
        }
    }

    // ---- 2. Shift to "origin at atom(column)" (tblite convention) ----
    //       and traceless quadrupole transform.
    for (int k = 0; k < 3; ++k) m_dp_int[k] = Eigen::MatrixXd::Zero(nao, nao);
    for (int k = 0; k < 6; ++k) m_qp_int[k] = Eigen::MatrixXd::Zero(nao, nao);

    for (int mu = 0; mu < nao; ++mu) {
        for (int nu = 0; nu < nao; ++nu) {
            const int iat = m_basis.ao2at[nu];   // column atom
            const double Rx = xyz_bohr[3*iat+0];
            const double Ry = xyz_bohr[3*iat+1];
            const double Rz = xyz_bohr[3*iat+2];
            const double Smn = m_S(mu, nu);
            const double dx = dp_global[0](mu, nu);
            const double dy = dp_global[1](mu, nu);
            const double dz = dp_global[2](mu, nu);

            m_dp_int[0](mu, nu) = dx - Rx * Smn;
            m_dp_int[1](mu, nu) = dy - Ry * Smn;
            m_dp_int[2](mu, nu) = dz - Rz * Smn;

            // Raw shifted quadrupole ⟨μ|(r-R)(r-R)|ν⟩
            const double qxx = qp_global_raw[0](mu, nu) - 2*Rx*dx + Rx*Rx*Smn;
            const double qxy = qp_global_raw[1](mu, nu) - Rx*dy - Ry*dx + Rx*Ry*Smn;
            const double qyy = qp_global_raw[2](mu, nu) - 2*Ry*dy + Ry*Ry*Smn;
            const double qxz = qp_global_raw[3](mu, nu) - Rx*dz - Rz*dx + Rx*Rz*Smn;
            const double qyz = qp_global_raw[4](mu, nu) - Ry*dz - Rz*dy + Ry*Rz*Smn;
            const double qzz = qp_global_raw[5](mu, nu) - 2*Rz*dz + Rz*Rz*Smn;

            // Traceless transform (tblite convention): Q_ab = 1.5·Q_raw_ab - 0.5·tr·δ_ab
            const double tr = 0.5 * (qxx + qyy + qzz);
            m_qp_int[0](mu, nu) = 1.5 * qxx - tr;
            m_qp_int[1](mu, nu) = 1.5 * qxy;
            m_qp_int[2](mu, nu) = 1.5 * qyy - tr;
            m_qp_int[3](mu, nu) = 1.5 * qxz;
            m_qp_int[4](mu, nu) = 1.5 * qyz;
            m_qp_int[5](mu, nu) = 1.5 * qzz - tr;
        }
    }

    // ---- 3. Coordination numbers (GFN2 double-exp form) ----
    std::vector<double> cn = cn_gfn(std::vector<int>(m_atoms.begin(), m_atoms.end()), xyz_bohr);

    // ---- 4. Damping radii and kernel parameters ----
    m_mp_mrad.resize(nat);
    m_mp_dkernel.resize(nat);
    m_mp_qkernel.resize(nat);
    using namespace curcuma::xtb::gfn2_params;
    for (int i = 0; i < nat; ++i) {
        const int z = m_atoms[i];
        const double vcn = p_vcn[z - 1];
        const double rad = p_rad[z - 1];
        const double arg = cn[i] - vcn - mp_shift;
        const double t1  = std::exp(-mp_kexp * arg);
        const double t2  = (mp_rmax - rad) / (1.0 + t1);
        m_mp_mrad[i]     = rad + t2;
        m_mp_dkernel[i]  = p_dkernel[z - 1];
        m_mp_qkernel[i]  = p_qkernel[z - 1];
    }

    // ---- 5. Interaction matrices ----
    for (int k = 0; k < 3; ++k) m_mp_amat_sd[k] = Eigen::MatrixXd::Zero(nat, nat);
    for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b) m_mp_amat_dd[a][b] = Eigen::MatrixXd::Zero(nat, nat);
    for (int k = 0; k < 6; ++k) m_mp_amat_sq[k] = Eigen::MatrixXd::Zero(nat, nat);

    for (int i = 0; i < nat; ++i) {
        for (int j = 0; j < nat; ++j) {
            if (i == j) continue;
            const double vx = xyz_bohr[3*i+0] - xyz_bohr[3*j+0];
            const double vy = xyz_bohr[3*i+1] - xyz_bohr[3*j+1];
            const double vz = xyz_bohr[3*i+2] - xyz_bohr[3*j+2];
            const double r1  = std::sqrt(vx*vx + vy*vy + vz*vz);
            const double g1  = 1.0 / r1;
            const double g3  = g1 * g1 * g1;
            const double g5  = g3 * g1 * g1;
            const double rr  = 0.5 * (m_mp_mrad[j] + m_mp_mrad[i]) * g1;
            const double fdmp3 = 1.0 / (1.0 + 6.0 * std::pow(rr, mp_dmp3));
            const double fdmp5 = 1.0 / (1.0 + 6.0 * std::pow(rr, mp_dmp5));

            // amat_sd (charge-dipole): row=target(j), col=source(i)
            m_mp_amat_sd[0](j, i) += vx * g3 * fdmp3;
            m_mp_amat_sd[1](j, i) += vy * g3 * fdmp3;
            m_mp_amat_sd[2](j, i) += vz * g3 * fdmp3;

            // amat_dd (dipole-dipole)
            const double dd_iso  = g3 * fdmp5;
            const double dd_anis = 3.0 * g5 * fdmp5;
            for (int a = 0; a < 3; ++a) {
                const double va = (a == 0) ? vx : ((a == 1) ? vy : vz);
                for (int b = 0; b < 3; ++b) {
                    const double vb = (b == 0) ? vx : ((b == 1) ? vy : vz);
                    m_mp_amat_dd[a][b](j, i) +=
                        ((a == b) ? dd_iso : 0.0) - va * vb * dd_anis;
                }
            }

            // amat_sq (charge-quadrupole), packed (xx, xy, yy, xz, yz, zz)
            m_mp_amat_sq[0](j, i) += vx * vx * g5 * fdmp5;
            m_mp_amat_sq[1](j, i) += 2.0 * vx * vy * g5 * fdmp5;
            m_mp_amat_sq[2](j, i) += vy * vy * g5 * fdmp5;
            m_mp_amat_sq[3](j, i) += 2.0 * vx * vz * g5 * fdmp5;
            m_mp_amat_sq[4](j, i) += 2.0 * vy * vz * g5 * fdmp5;
            m_mp_amat_sq[5](j, i) += vz * vz * g5 * fdmp5;
        }
    }

    m_mp_initialized = true;
}

/* ------------------------------------------------------------------ *
 *  addMultipolePotential()                                           *
 *                                                                    *
 *  Adds multipole contributions to the potential:                    *
 *    pot.v_at   += vat_extra  (scalar shift from dipole/quadrupole)  *
 *    pot.v_dp    = vdp        (dipole potential, 3 × nat)           *
 *    pot.v_qp    = vqp        (quadrupole potential, 6 × nat)        *
 *                                                                    *
 *  Must be called AFTER addCoulombShellPotential/addThirdOrder so    *
 *  those contributions are already in v_sh.                          *
 * ------------------------------------------------------------------ */
void XTB::addMultipolePotential(Potential& pot) const
{
    if (!m_mp_initialized) return;
    const int nat = m_atomcount;

    // --- vdp(k, iat) = Σ_jat amat_sd[k](iat,jat)·q_at(jat)
    //                 + Σ_a   amat_dd[k][a](iat,jat)·dpat(a,jat)
    //                 + 2·dkernel·dpat(k,iat)
    pot.v_dp.resize(3, nat);
    Eigen::MatrixXd& vdp = pot.v_dp;
    for (int i = 0; i < nat; ++i) {
        double vd0 = 0.0, vd1 = 0.0, vd2 = 0.0;
        for (int j = 0; j < nat; ++j) {
            vd0 += m_mp_amat_sd[0](i, j) * m_wfn.q_at(j);
            vd1 += m_mp_amat_sd[1](i, j) * m_wfn.q_at(j);
            vd2 += m_mp_amat_sd[2](i, j) * m_wfn.q_at(j);
            vd0 += m_mp_amat_dd[0][0](i, j) * m_wfn.dp_at(0, j)
                 + m_mp_amat_dd[0][1](i, j) * m_wfn.dp_at(1, j)
                 + m_mp_amat_dd[0][2](i, j) * m_wfn.dp_at(2, j);
            vd1 += m_mp_amat_dd[1][0](i, j) * m_wfn.dp_at(0, j)
                 + m_mp_amat_dd[1][1](i, j) * m_wfn.dp_at(1, j)
                 + m_mp_amat_dd[1][2](i, j) * m_wfn.dp_at(2, j);
            vd2 += m_mp_amat_dd[2][0](i, j) * m_wfn.dp_at(0, j)
                 + m_mp_amat_dd[2][1](i, j) * m_wfn.dp_at(1, j)
                 + m_mp_amat_dd[2][2](i, j) * m_wfn.dp_at(2, j);
        }
        vdp(0, i) = vd0 + 2.0 * m_mp_dkernel[i] * m_wfn.dp_at(0, i);
        vdp(1, i) = vd1 + 2.0 * m_mp_dkernel[i] * m_wfn.dp_at(1, i);
        vdp(2, i) = vd2 + 2.0 * m_mp_dkernel[i] * m_wfn.dp_at(2, i);
    }

    // --- vqp(k, iat) = Σ_jat amat_sq[k](iat,jat)·q_at(jat)
    //                 + 2·qkernel·qpat(k,iat)·mpscale_q[k]
    pot.v_qp.resize(6, nat);
    Eigen::MatrixXd& vqp = pot.v_qp;
    static const double mpscale_q[6] = {1.0, 2.0, 1.0, 2.0, 2.0, 1.0};
    for (int i = 0; i < nat; ++i) {
        for (int k = 0; k < 6; ++k) {
            double v = 0.0;
            for (int j = 0; j < nat; ++j)
                v += m_mp_amat_sq[k](i, j) * m_wfn.q_at(j);
            vqp(k, i) = v + 2.0 * m_mp_qkernel[i] * m_wfn.qp_at(k, i) * mpscale_q[k];
        }
    }

    // --- vat_extra(iat) = Σ_jat Σ_k amat_sd[k](jat,iat)·dpat(k,jat)
    //                    + Σ_jat Σ_k amat_sq[k](jat,iat)·qpat(k,jat)
    pot.v_at.resize(nat);
    pot.v_at.setZero();
    for (int i = 0; i < nat; ++i) {
        double acc = 0.0;
        for (int j = 0; j < nat; ++j) {
            acc += m_mp_amat_sd[0](j, i) * m_wfn.dp_at(0, j)
                 + m_mp_amat_sd[1](j, i) * m_wfn.dp_at(1, j)
                 + m_mp_amat_sd[2](j, i) * m_wfn.dp_at(2, j);
            for (int k = 0; k < 6; ++k)
                acc += m_mp_amat_sq[k](j, i) * m_wfn.qp_at(k, j);
        }
        pot.v_at(i) += acc;  // add to existing v_at (from third-order, etc.)
    }
}

/* ------------------------------------------------------------------ *
 *  energyMultipole()                                                *
 *                                                                    *
 *  Multipole contribution to total SCC energy:                       *
 *    E_mp = Σ_i Σ_j [                                                *
 *      Σ_k dpat(k,i) · amat_sd[k](i,j) · qat(j)                     *
 *    + Σ_a Σ_b 0.5 · dpat(a,i) · amat_dd[a][b](i,j) · dpat(b,j)    *
 *    + Σ_k qpat(k,i) · amat_sq[k](i,j) · qat(j) ]                   *
 *    + Σ_i [                                                         *
 *      Σ_k dkernel[i] · dpat(k,i)²                                   *
 *    + Σ_k qkernel[i] · qpat(k,i)² · mpscale_q[k] ]                 *
 * ------------------------------------------------------------------ */
double XTB::energyMultipole() const
{
    if (!m_mp_initialized) return 0.0;
    const int nat = m_atomcount;
    double e = 0.0;

    static const double mpscale_q[6] = {1.0, 2.0, 1.0, 2.0, 2.0, 1.0};

    for (int i = 0; i < nat; ++i) {
        for (int j = 0; j < nat; ++j) {
            // SD: dpat(k,i) · amat_sd[k](i,j) · qat(j)
            for (int k = 0; k < 3; ++k)
                e += m_wfn.dp_at(k, i) * m_mp_amat_sd[k](i, j) * m_wfn.q_at(j);

            // DD: 0.5 · dpat(a,i) · amat_dd[a][b](i,j) · dpat(b,j)
            for (int a = 0; a < 3; ++a)
                for (int b = 0; b < 3; ++b)
                    e += 0.5 * m_wfn.dp_at(a, i) * m_mp_amat_dd[a][b](i, j) * m_wfn.dp_at(b, j);

            // SQ: qpat(k,i) · amat_sq[k](i,j) · qat(j)
            for (int k = 0; k < 6; ++k)
                e += m_wfn.qp_at(k, i) * m_mp_amat_sq[k](i, j) * m_wfn.q_at(j);
        }
        // On-site dipole XC: dkernel · dpat(k,i)²
        for (int k = 0; k < 3; ++k)
            e += m_mp_dkernel[i] * m_wfn.dp_at(k, i) * m_wfn.dp_at(k, i);

        // On-site quadrupole XC: qkernel · qpat(k,i)² · mpscale_q[k]
        for (int k = 0; k < 6; ++k)
            e += m_mp_qkernel[i] * m_wfn.qp_at(k, i) * m_wfn.qp_at(k, i) * mpscale_q[k];
    }
    return e;
}

} // namespace curcuma::xtb
