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

} // namespace curcuma::xtb::multipole_ints
