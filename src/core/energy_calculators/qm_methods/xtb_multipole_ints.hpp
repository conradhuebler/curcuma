/*
 * STO/CGTO Multipole Integrals for Native xTB
 * Copyright (C) 2025 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Provides dipole and quadrupole moment integrals over the existing
 * STO->CGTO expansion (STO_CGTO.hpp), for the s/p subset used by the
 * GFN1/GFN2 basis on H, He, Li-Ne.
 *
 * Conventions chosen for simplicity:
 *   - Integrals are computed with the moment operator origin at 0
 *     (global Cartesian origin), so the resulting nao x nao matrices
 *     are symmetric (real-valued spherical basis).
 *   - Raw Cartesian quadrupole components are returned (xx, xy, yy,
 *     xz, yz, zz) with no traceless transform.  The shift to each
 *     atom's center and the tblite traceless convention
 *         Q_ab(traceless) = 1.5 * Q_ab(raw) - 0.5 * tr(Q_raw) * delta_ab
 *     are applied by the caller during Mulliken partitioning / Fock
 *     assembly.
 *   - AO ordering for p shells matches the caller's convention, which
 *     currently maps tblite's spherical (py, pz, px) onto the type
 *     codes 1=px, 2=py, 3=pz via an external map (see test_xtb_overlap.cpp).
 *     Only s and p shells are supported here.
 *
 * Reference: external/tblite/src/tblite/integral/multipole.f90 (routines
 * overlap_1d, multipole_3d, multipole_cgto).
 *
 * Claude Generated (Phase 3.5, Apr 2026). GPL-3.0.
 */

#pragma once

#include "STO_CGTO.hpp"

#include <cmath>

namespace curcuma::xtb::multipole_ints {

// -----------------------------------------------------------------------------
// 1D moment helper.  Returns the integral
//      I = int_{-inf}^{inf} (u+PA)^la * (u+PB)^lb * (u+P)^n * exp(-gamma u^2) du
// with la, lb in {0,1} and n in {0,1,2}, i.e. degree up to 4.
//
// Even-order moments:  M_0 = sqrt(pi/gamma),  M_2 = M_0/(2 gamma),
//                      M_4 = 3 M_0 / (4 gamma^2).
// Odd-order moments vanish.
// -----------------------------------------------------------------------------
inline double moment1d(int la, int lb, int n,
                       double PA, double PB, double P, double gamma)
{
    double a0 = (la == 0) ? 1.0 : PA;
    double a1 = (la == 0) ? 0.0 : 1.0;
    double b0 = (lb == 0) ? 1.0 : PB;
    double b1 = (lb == 0) ? 0.0 : 1.0;

    double p0, p1, p2;
    if (n == 0)      { p0 = 1.0;  p1 = 0.0;  p2 = 0.0; }
    else if (n == 1) { p0 = P;    p1 = 1.0;  p2 = 0.0; }
    else             { p0 = P*P;  p1 = 2.0*P; p2 = 1.0; }

    // ab(u) = a(u) * b(u), degree <= 2
    double ab0 = a0 * b0;
    double ab1 = a0 * b1 + a1 * b0;
    double ab2 = a1 * b1;

    // c(u) = ab(u) * p(u), degree <= 4 (only even powers contribute)
    double c0 = ab0 * p0;
    double c2 = ab0 * p2 + ab1 * p1 + ab2 * p0;
    double c4 = ab2 * p2;

    const double M0 = std::sqrt(M_PI / gamma);
    const double M2 = M0 / (2.0 * gamma);
    const double M4 = 3.0 * M0 / (4.0 * gamma * gamma);
    return c0 * M0 + c2 * M2 + c4 * M4;
}

// Cartesian (lx, ly, lz) indicator for a type code used by STO_CGTO.hpp.
//   0 = s (0,0,0), 1 = px (1,0,0), 2 = py (0,1,0), 3 = pz (0,0,1)
inline void type_to_cart(int type, int& lx, int& ly, int& lz)
{
    switch (type) {
    case 1: lx = 1; ly = 0; lz = 0; break;
    case 2: lx = 0; ly = 1; lz = 0; break;
    case 3: lx = 0; ly = 0; lz = 1; break;
    default: lx = 0; ly = 0; lz = 0; break;
    }
}

// -----------------------------------------------------------------------------
// Overlap + dipole + raw-Cartesian quadrupole between two primitive
// Gaussians at centres A, B with the indicator Cartesian powers encoded
// by type_a, type_b.  Origin for the moment operator is 0.
// Output ordering for the quadrupole: (xx, xy, yy, xz, yz, zz).
// -----------------------------------------------------------------------------
inline void primitive_multipole(double alpha_a, double alpha_b,
                                 double Ax, double Ay, double Az,
                                 double Bx, double By, double Bz,
                                 int type_a, int type_b,
                                 double& S, double D[3], double Q[6])
{
    int lxa, lya, lza, lxb, lyb, lzb;
    type_to_cart(type_a, lxa, lya, lza);
    type_to_cart(type_b, lxb, lyb, lzb);

    const double gamma = alpha_a + alpha_b;
    const double Px = (alpha_a * Ax + alpha_b * Bx) / gamma;
    const double Py = (alpha_a * Ay + alpha_b * By) / gamma;
    const double Pz = (alpha_a * Az + alpha_b * Bz) / gamma;
    const double dx = Ax - Bx, dy = Ay - By, dz = Az - Bz;
    const double R2 = dx*dx + dy*dy + dz*dz;
    const double K  = std::exp(-alpha_a * alpha_b / gamma * R2);

    const double PAx = Px - Ax, PAy = Py - Ay, PAz = Pz - Az;
    const double PBx = Px - Bx, PBy = Py - By, PBz = Pz - Bz;

    const double Sx0 = moment1d(lxa, lxb, 0, PAx, PBx, Px, gamma);
    const double Sy0 = moment1d(lya, lyb, 0, PAy, PBy, Py, gamma);
    const double Sz0 = moment1d(lza, lzb, 0, PAz, PBz, Pz, gamma);

    const double Sx1 = moment1d(lxa, lxb, 1, PAx, PBx, Px, gamma);
    const double Sy1 = moment1d(lya, lyb, 1, PAy, PBy, Py, gamma);
    const double Sz1 = moment1d(lza, lzb, 1, PAz, PBz, Pz, gamma);

    const double Sx2 = moment1d(lxa, lxb, 2, PAx, PBx, Px, gamma);
    const double Sy2 = moment1d(lya, lyb, 2, PAy, PBy, Py, gamma);
    const double Sz2 = moment1d(lza, lzb, 2, PAz, PBz, Pz, gamma);

    S    = K * Sx0 * Sy0 * Sz0;
    D[0] = K * Sx1 * Sy0 * Sz0;
    D[1] = K * Sx0 * Sy1 * Sz0;
    D[2] = K * Sx0 * Sy0 * Sz1;
    Q[0] = K * Sx2 * Sy0 * Sz0;   // xx
    Q[1] = K * Sx1 * Sy1 * Sz0;   // xy
    Q[2] = K * Sx0 * Sy2 * Sz0;   // yy
    Q[3] = K * Sx1 * Sy0 * Sz1;   // xz
    Q[4] = K * Sx0 * Sy1 * Sz1;   // yz
    Q[5] = K * Sx0 * Sy0 * Sz2;   // zz
}

// -----------------------------------------------------------------------------
// CGTO-level overlap + dipole + raw-Cartesian quadrupole for one AO pair.
// The CGTO coefficients are expected to be pre-normalised (see
// STO_CGTO::slater_to_gauss), matching the convention used by
// CGTO::cgto_overlap.
// -----------------------------------------------------------------------------
inline void cgto_multipole(const CGTO::Shell& shell_a, const CGTO::Shell& shell_b,
                           double xa, double ya, double za,
                           double xb, double yb, double zb,
                           int type_a, int type_b,
                           double& S, double D[3], double Q[6])
{
    S = 0.0;
    for (int k = 0; k < 3; ++k) D[k] = 0.0;
    for (int k = 0; k < 6; ++k) Q[k] = 0.0;

    for (int i = 0; i < shell_a.nprim; ++i) {
        const double ai = shell_a.alpha[i];
        const double ci = shell_a.coeff[i];
        for (int j = 0; j < shell_b.nprim; ++j) {
            const double aj = shell_b.alpha[j];
            const double cj = shell_b.coeff[j];
            double s, d[3], q[6];
            primitive_multipole(ai, aj, xa, ya, za, xb, yb, zb,
                                type_a, type_b, s, d, q);
            const double cc = ci * cj;
            S    += cc * s;
            D[0] += cc * d[0]; D[1] += cc * d[1]; D[2] += cc * d[2];
            for (int k = 0; k < 6; ++k) Q[k] += cc * q[k];
        }
    }
}

// =============================================================================
// AP5b: Gradient routines for multipole integrals (Pulay integral term)
//
// Provides gradients of the *transformed* multipole integrals used in the
// Fock matrix (origin at B = atom jat, traceless quadrupole), needed for the
// GFN2 Hamiltonian/Pulay gradient.
//
// Reference: tblite/integral/multipole.f90:multipole_grad_cgto (line 522)
//            tblite/xtb/h0.f90:get_hamiltonian_gradient (line 338)
//
// Math: For a primitive pair at A, B with γ = αa+αb:
//   G(nx,ny,nz) = K · Ix[nx] · Iy[ny] · Iz[nz]
//   where Iq[n] = moment1d(la_q, lb_q, n, ...)  and  K = exp(-αa·αb/γ · R²)
//
//   d/dA_l G = K · (-2αa·αb/γ · R_l) · Ix[nx]·Iy[ny]·Iz[nz]
//            + K · [la_l·(-αb/γ)·Ilo_l[n_l] + n_l·(αa/γ)·I_l[n_l-1]] · ∏_{q≠l} I_q[n_q]
//
//   d/dB_l G = K · (+2αa·αb/γ · R_l) · Ix[nx]·Iy[ny]·Iz[nz]
//            + K · [la_l·(αb/γ)·Ilo_l[n_l] + lb_l·(-αa/γ)·Jlo_l[n_l] + n_l·(αb/γ)·I_l[n_l-1]]
//              · ∏_{q≠l} I_q[n_q]
//
//   where Ilo_l = I with la_l lowered by 1; Jlo_l = I with lb_l lowered by 1.
//   Only lowered (not raised) terms appear — no d-type intermediates needed.
//
// The transformed integrals (origin at B, traceless Q) are:
//   dp_int[k]   = D_global[k] - B_k · S
//   qraw_ab     = Q_global[ab] - B_a·D_global[b] - B_b·D_global[a] + B_a·B_b·S
//   qp_int[k]   = 1.5·qraw_k - 0.5·tr(qraw)·δ_{k diag}
// =============================================================================

// -----------------------------------------------------------------------------
// Primitive-level gradient of global-origin integrals S, D[3], Q[6].
//
// Outputs:
//   S, D[3], Q[6]           — integrals (same as primitive_multipole)
//   dS_dA[3], dS_dB[3]     — dS/dA_l and dS/dB_l
//   dD_dA[3][3]             — d(D_global[k])/d(A_l):  [l][k]
//   dD_dB[3][3]             — d(D_global[k])/d(B_l):  [l][k]
//   dQ_dA[3][6]             — d(Q_global[q])/d(A_l):  [l][q]
//   dQ_dB[3][6]             — d(Q_global[q])/d(B_l):  [l][q]
// -----------------------------------------------------------------------------
inline void primitive_multipole_grad(
    double alpha_a, double alpha_b,
    double Ax, double Ay, double Az,
    double Bx, double By, double Bz,
    int type_a, int type_b,
    double& S, double D[3], double Q[6],
    double dS_dA[3], double dS_dB[3],
    double dD_dA[3][3], double dD_dB[3][3],
    double dQ_dA[3][6], double dQ_dB[3][6])
{
    int lxa, lya, lza, lxb, lyb, lzb;
    type_to_cart(type_a, lxa, lya, lza);
    type_to_cart(type_b, lxb, lyb, lzb);

    const double gamma = alpha_a + alpha_b;
    const double Px = (alpha_a*Ax + alpha_b*Bx) / gamma;
    const double Py = (alpha_a*Ay + alpha_b*By) / gamma;
    const double Pz = (alpha_a*Az + alpha_b*Bz) / gamma;
    const double PAx = Px - Ax, PAy = Py - Ay, PAz = Pz - Az;
    const double PBx = Px - Bx, PBy = Py - By, PBz = Pz - Bz;
    const double Rx = Ax - Bx, Ry = Ay - By, Rz = Az - Bz;
    const double R2 = Rx*Rx + Ry*Ry + Rz*Rz;
    const double K  = std::exp(-alpha_a * alpha_b / gamma * R2);

    // 1D moments: Iq[n] = moment1d(la_q, lb_q, n, ...)  for n=0,1,2
    double Ix[3], Iy[3], Iz[3];
    for (int n = 0; n <= 2; ++n) {
        Ix[n] = moment1d(lxa, lxb, n, PAx, PBx, Px, gamma);
        Iy[n] = moment1d(lya, lyb, n, PAy, PBy, Py, gamma);
        Iz[n] = moment1d(lza, lzb, n, PAz, PBz, Pz, gamma);
    }
    // Lowered A moments (la→la-1): only non-zero if la_q > 0
    double IxAlo[3] = {}, IyAlo[3] = {}, IzAlo[3] = {};
    if (lxa > 0) for (int n=0;n<=2;++n) IxAlo[n]=moment1d(lxa-1,lxb,n,PAx,PBx,Px,gamma);
    if (lya > 0) for (int n=0;n<=2;++n) IyAlo[n]=moment1d(lya-1,lyb,n,PAy,PBy,Py,gamma);
    if (lza > 0) for (int n=0;n<=2;++n) IzAlo[n]=moment1d(lza-1,lzb,n,PAz,PBz,Pz,gamma);
    // Lowered B moments (lb→lb-1): only non-zero if lb_q > 0
    double IxBlo[3] = {}, IyBlo[3] = {}, IzBlo[3] = {};
    if (lxb > 0) for (int n=0;n<=2;++n) IxBlo[n]=moment1d(lxa,lxb-1,n,PAx,PBx,Px,gamma);
    if (lyb > 0) for (int n=0;n<=2;++n) IyBlo[n]=moment1d(lya,lyb-1,n,PAy,PBy,Py,gamma);
    if (lzb > 0) for (int n=0;n<=2;++n) IzBlo[n]=moment1d(lza,lzb-1,n,PAz,PBz,Pz,gamma);

    // Assemble integrals
    S    = K * Ix[0] * Iy[0] * Iz[0];
    D[0] = K * Ix[1] * Iy[0] * Iz[0];
    D[1] = K * Ix[0] * Iy[1] * Iz[0];
    D[2] = K * Ix[0] * Iy[0] * Iz[1];
    Q[0] = K * Ix[2] * Iy[0] * Iz[0];  // xx
    Q[1] = K * Ix[1] * Iy[1] * Iz[0];  // xy
    Q[2] = K * Ix[0] * Iy[2] * Iz[0];  // yy
    Q[3] = K * Ix[1] * Iy[0] * Iz[1];  // xz
    Q[4] = K * Ix[0] * Iy[1] * Iz[1];  // yz
    Q[5] = K * Ix[0] * Iy[0] * Iz[2];  // zz

    // K gradient factors: dK/dA_l = K·(-2αa·αb/γ·R_l)
    const double kfA[3] = { -2.0*alpha_a*alpha_b/gamma*Rx,
                             -2.0*alpha_a*alpha_b/gamma*Ry,
                             -2.0*alpha_a*alpha_b/gamma*Rz };
    // dK/dB_l = +2αa·αb/γ·R_l  (= -dK/dA_l)

    // Helper: d/dA_l of K·Ix[nx]·Iy[ny]·Iz[nz]
    // = K·kfA[l]·Ix[nx]·Iy[ny]·Iz[nz]
    // + K · (l==0 ? [lxa·(-αb/γ)·IxAlo[nx] + nx·(αa/γ)·Ix[nx-1]] · Iy[ny]·Iz[nz]
    //      : l==1 ? Ix[nx] · [lya·(-αb/γ)·IyAlo[ny] + ny·(αa/γ)·Iy[ny-1]] · Iz[nz]
    //      :        Ix[nx]·Iy[ny] · [lza·(-αb/γ)·IzAlo[nz] + nz·(αa/γ)·Iz[nz-1]])
    // d/dA_l of K·Ix[nx]·Iy[ny]·Iz[nz]:
    // = K·kfA[l]·base
    // + K · d/dA_l[I_{l_q}[n_q]] · product of other I
    //
    // d/dAx[Ix[nx]] = dPAx/dAx · lxa · IxAlo  +  dPBx/dAx · lxb · IxBlo  +  dPx/dAx · nx · Ix[nx-1]
    //               = (-αb/γ) · lxa · IxAlo   +  (αa/γ) · lxb · IxBlo   +  (αa/γ) · nx · Ix[nx-1]
    auto dG_dA = [&](int nx, int ny, int nz, int l) -> double {
        const double base = Ix[nx] * Iy[ny] * Iz[nz];
        double r = K * kfA[l] * base;
        if (l == 0) {
            double t = (lxa > 0 ? lxa * (-alpha_b/gamma) * IxAlo[nx] : 0.0)
                     + (lxb > 0 ? lxb * ( alpha_a/gamma) * IxBlo[nx] : 0.0)
                     + (nx  > 0 ? nx  * ( alpha_a/gamma) *   Ix[nx-1] : 0.0);
            r += K * t * Iy[ny] * Iz[nz];
        } else if (l == 1) {
            double t = (lya > 0 ? lya * (-alpha_b/gamma) * IyAlo[ny] : 0.0)
                     + (lyb > 0 ? lyb * ( alpha_a/gamma) * IyBlo[ny] : 0.0)
                     + (ny  > 0 ? ny  * ( alpha_a/gamma) *   Iy[ny-1] : 0.0);
            r += K * Ix[nx] * t * Iz[nz];
        } else {
            double t = (lza > 0 ? lza * (-alpha_b/gamma) * IzAlo[nz] : 0.0)
                     + (lzb > 0 ? lzb * ( alpha_a/gamma) * IzBlo[nz] : 0.0)
                     + (nz  > 0 ? nz  * ( alpha_a/gamma) *   Iz[nz-1] : 0.0);
            r += K * Ix[nx] * Iy[ny] * t;
        }
        return r;
    };

    // d/dB_l of K·Ix[nx]·Iy[ny]·Iz[nz]:
    // K·(-kfA[l])·base + K·(la·(αb/γ)·Ilo + lb·(-αa/γ)·Jlo + n·(αb/γ)·I[n-1])·rest
    auto dG_dB = [&](int nx, int ny, int nz, int l) -> double {
        const double base = Ix[nx] * Iy[ny] * Iz[nz];
        double r = K * (-kfA[l]) * base;
        if (l == 0) {
            double t = (lxa > 0 ? lxa * ( alpha_b/gamma) * IxAlo[nx] : 0.0)
                     + (lxb > 0 ? lxb * (-alpha_a/gamma) * IxBlo[nx] : 0.0)
                     + (nx  > 0 ? nx  * ( alpha_b/gamma) *   Ix[nx-1] : 0.0);
            r += K * t * Iy[ny] * Iz[nz];
        } else if (l == 1) {
            double t = (lya > 0 ? lya * ( alpha_b/gamma) * IyAlo[ny] : 0.0)
                     + (lyb > 0 ? lyb * (-alpha_a/gamma) * IyBlo[ny] : 0.0)
                     + (ny  > 0 ? ny  * ( alpha_b/gamma) *   Iy[ny-1] : 0.0);
            r += K * Ix[nx] * t * Iz[nz];
        } else {
            double t = (lza > 0 ? lza * ( alpha_b/gamma) * IzAlo[nz] : 0.0)
                     + (lzb > 0 ? lzb * (-alpha_a/gamma) * IzBlo[nz] : 0.0)
                     + (nz  > 0 ? nz  * ( alpha_b/gamma) *   Iz[nz-1] : 0.0);
            r += K * Ix[nx] * Iy[ny] * t;
        }
        return r;
    };

    // n-tuples for S, D[3], Q[6]
    // S: (0,0,0); D: (1,0,0),(0,1,0),(0,0,1); Q: (2,0,0),(1,1,0),(0,2,0),(1,0,1),(0,1,1),(0,0,2)
    static const int qnx[6] = {2,1,0,1,0,0};
    static const int qny[6] = {0,1,2,0,1,0};
    static const int qnz[6] = {0,0,0,1,1,2};

    for (int l = 0; l < 3; ++l) {
        dS_dA[l] = dG_dA(0, 0, 0, l);
        dS_dB[l] = dG_dB(0, 0, 0, l);
        for (int k = 0; k < 3; ++k) {
            // D[k]: nx=delta_k0, ny=delta_k1, nz=delta_k2
            const int dnx = (k==0), dny = (k==1), dnz = (k==2);
            dD_dA[l][k] = dG_dA(dnx, dny, dnz, l);
            dD_dB[l][k] = dG_dB(dnx, dny, dnz, l);
        }
        for (int q = 0; q < 6; ++q) {
            dQ_dA[l][q] = dG_dA(qnx[q], qny[q], qnz[q], l);
            dQ_dB[l][q] = dG_dB(qnx[q], qny[q], qnz[q], l);
        }
    }
}

// -----------------------------------------------------------------------------
// CGTO-level gradient of transformed multipole integrals (origin at B, traceless Q).
//
// Sums primitive_multipole_grad over contracted pairs, then applies:
//   dp_int[k](A,B) = D_global[k] - B_k · S   (origin shift to B)
//   qp_int[q]      = 1.5·qraw[q] - 0.5·tr(qraw)·δ_{q diag}  (traceless, shifted)
//
// Outputs:
//   S_out, D_out[3], Q_out[6]  — transformed integral values
//   dD_dA[l][k]  = d(dp_int[k])/d(A_l)   (l=spatial dir, k=dipole comp)
//   dD_dB[l][k]  = d(dp_int[k])/d(B_l)
//   dQ_dA[l][q]  = d(qp_int[q])/d(A_l)   (q=0..5 packed)
//   dQ_dB[l][q]  = d(qp_int[q])/d(B_l)
// -----------------------------------------------------------------------------
inline void cgto_multipole_grad_transformed(
    const CGTO::Shell& shell_a, const CGTO::Shell& shell_b,
    double Ax, double Ay, double Az,  // center A = atom iat
    double Bx, double By, double Bz,  // center B = atom jat
    int type_a, int type_b,
    double& S_out, double D_out[3], double Q_out[6],
    double dD_dA[3][3], double dD_dB[3][3],
    double dQ_dA[3][6], double dQ_dB[3][6])
{
    // Accumulate over primitives
    double S_raw = 0.0, D_raw[3] = {}, Q_raw[6] = {};
    double dS_dA[3] = {}, dS_dB[3] = {};
    double dDg_dA[3][3] = {}, dDg_dB[3][3] = {};  // gradient of D_global
    double dQg_dA[3][6] = {}, dQg_dB[3][6] = {};  // gradient of Q_global

    for (int i = 0; i < shell_a.nprim; ++i) {
        const double ai = shell_a.alpha[i];
        const double ci = shell_a.coeff[i];
        for (int j = 0; j < shell_b.nprim; ++j) {
            const double aj = shell_b.alpha[j];
            const double cj = shell_b.coeff[j];
            const double cc = ci * cj;

            double S_p, D_p[3], Q_p[6];
            double dS_A[3], dS_B[3], dD_A[3][3], dD_B[3][3], dQ_A[3][6], dQ_B[3][6];
            primitive_multipole_grad(ai, aj, Ax, Ay, Az, Bx, By, Bz,
                                     type_a, type_b,
                                     S_p, D_p, Q_p,
                                     dS_A, dS_B, dD_A, dD_B, dQ_A, dQ_B);

            S_raw += cc * S_p;
            for (int k=0;k<3;++k) D_raw[k] += cc * D_p[k];
            for (int q=0;q<6;++q) Q_raw[q] += cc * Q_p[q];
            for (int l=0;l<3;++l) {
                dS_dA[l] += cc * dS_A[l];
                dS_dB[l] += cc * dS_B[l];
                for (int k=0;k<3;++k) {
                    dDg_dA[l][k] += cc * dD_A[l][k];
                    dDg_dB[l][k] += cc * dD_B[l][k];
                }
                for (int q=0;q<6;++q) {
                    dQg_dA[l][q] += cc * dQ_A[l][q];
                    dQg_dB[l][q] += cc * dQ_B[l][q];
                }
            }
        }
    }

    // Apply origin shift to B: dp_int[k] = D_global[k] - B_k · S
    const double B[3] = {Bx, By, Bz};
    for (int k=0;k<3;++k) D_out[k] = D_raw[k] - B[k] * S_raw;
    for (int l=0;l<3;++l) {
        for (int k=0;k<3;++k) {
            dD_dA[l][k] = dDg_dA[l][k] - B[k] * dS_dA[l];
            // d/dB_l: also d(-B_k)/dB_l = -δ_{kl}
            dD_dB[l][k] = dDg_dB[l][k] - (k==l ? S_raw : 0.0) - B[k] * dS_dB[l];
        }
    }

    // Apply shift and traceless transform to quadrupole
    // qraw_ab = Q_global[ab] - B_a·D_global[b] - B_b·D_global[a] + B_a·B_b·S
    // Packed indices: xx=0(a=0,b=0), xy=1(a=0,b=1), yy=2(a=1,b=1),
    //                xz=3(a=0,b=2), yz=4(a=1,b=2), zz=5(a=2,b=2)
    static const int qa[6] = {0,0,1,0,1,2};
    static const int qb[6] = {0,1,1,2,2,2};

    double qraw[6];
    for (int q=0;q<6;++q) {
        const int a = qa[q], b = qb[q];
        qraw[q] = Q_raw[q] - B[a]*D_raw[b] - B[b]*D_raw[a] + B[a]*B[b]*S_raw;
        // For diagonal (a==b): coefficient of Q_raw is 1; for off-diag: 2 in Q_raw[xy]?
        // Note: Q_raw[1]=xy = K·Ix1·Iy1·Iz0, so it's the full cross term, not 2x.
        // B[a]*D_raw[b] + B[b]*D_raw[a] handles both when a==b: 2*B[a]*D_raw[a].
    }
    const double tr_raw = 0.5 * (qraw[0] + qraw[2] + qraw[5]);
    Q_out[0] = 1.5*qraw[0] - tr_raw;
    Q_out[1] = 1.5*qraw[1];
    Q_out[2] = 1.5*qraw[2] - tr_raw;
    Q_out[3] = 1.5*qraw[3];
    Q_out[4] = 1.5*qraw[4];
    Q_out[5] = 1.5*qraw[5] - tr_raw;

    // Gradient of qraw[q] w.r.t. A_l:
    //   d(qraw_ab)/dA_l = dQ_global[ab]/dA_l - B_a·dD_global[b]/dA_l - B_b·dD_global[a]/dA_l + B_a·B_b·dS/dA_l
    for (int l=0;l<3;++l) {
        double dqraw_A[6], dqraw_B[6];
        for (int q=0;q<6;++q) {
            const int a = qa[q], b = qb[q];
            dqraw_A[q] = dQg_dA[l][q] - B[a]*dDg_dA[l][b] - B[b]*dDg_dA[l][a] + B[a]*B[b]*dS_dA[l];
            // d/dB_l: additional -δ_{al}·D_global[b] - δ_{bl}·D_global[a] + (δ_{al}·B_b+B_a·δ_{bl})·S
            dqraw_B[q] = dQg_dB[l][q]
                       - (a==l ? D_raw[b] : 0.0) - B[a]*dDg_dB[l][b]
                       - (b==l ? D_raw[a] : 0.0) - B[b]*dDg_dB[l][a]
                       + (a==l ? B[b] : 0.0)*S_raw + B[a]*(b==l ? S_raw : 0.0)
                       + B[a]*B[b]*dS_dB[l];
        }
        // Apply traceless transform to gradient: d(qp)/dR = 1.5·d(qraw)/dR - 0.5·d(tr)/dR·δ_diag
        const double dtr_A = 0.5*(dqraw_A[0] + dqraw_A[2] + dqraw_A[5]);
        const double dtr_B = 0.5*(dqraw_B[0] + dqraw_B[2] + dqraw_B[5]);
        dQ_dA[l][0] = 1.5*dqraw_A[0] - dtr_A;
        dQ_dA[l][1] = 1.5*dqraw_A[1];
        dQ_dA[l][2] = 1.5*dqraw_A[2] - dtr_A;
        dQ_dA[l][3] = 1.5*dqraw_A[3];
        dQ_dA[l][4] = 1.5*dqraw_A[4];
        dQ_dA[l][5] = 1.5*dqraw_A[5] - dtr_A;
        dQ_dB[l][0] = 1.5*dqraw_B[0] - dtr_B;
        dQ_dB[l][1] = 1.5*dqraw_B[1];
        dQ_dB[l][2] = 1.5*dqraw_B[2] - dtr_B;
        dQ_dB[l][3] = 1.5*dqraw_B[3];
        dQ_dB[l][4] = 1.5*dqraw_B[4];
        dQ_dB[l][5] = 1.5*dqraw_B[5] - dtr_B;
    }
}

} // namespace curcuma::xtb::multipole_ints
