/*
 * Shared AO/shell helpers + cartesian->spherical d transform for native xTB
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Single source of truth (resolves X-I5) for:
 *   - as_cgto_shell()  : internal CGTOShell -> CGTO::Shell converter
 *   - ao_to_type()     : s/p AO -> STO_CGTO scalar type code (-1 for d)
 *   - the d-shell machinery (X-I1): cartesian(6)->spherical(5) transform
 *     (tblite integral/trafo.f90 dtrafo), the general 1-D Gaussian moment,
 *     and a shell-pair spherical overlap block.
 *
 * Design rule for byte-identity: callers branch per shell-pair. Pairs with
 * no d shell keep using the existing scalar s/p kernels (cgto_overlap etc.);
 * only d-touching pairs use sphericalOverlapBlock(). Molecules without any d
 * shell therefore never touch this code path and stay bit-identical.
 *
 * Reference: external/tblite/src/tblite/integral/trafo.f90 (dtrafo),
 *            external/tblite/src/tblite/basis/slater.f90 (normalization).
 *
 * Claude Generated (X-I1, June 2026). GPL-3.0.
 */

#pragma once

#include "STO_CGTO.hpp"

#include <cmath>

namespace curcuma::xtb {

// ---------------------------------------------------------------------------
// Shared shell helpers (were duplicated in xtb_h0 / xtb_gradient / xtb_response
// / xtb_multipole as as_cgto_shell{,_g,_r} and ao_to_type{,_g,_r}).
// ---------------------------------------------------------------------------

// Convert an internal CGTOShell (xtb_native.h) to CGTO::Shell. Templated to
// avoid pulling the heavy xtb_native.h into this header; the only requirement
// is the .ang / .alpha / .coeff members.
template <class ShellT>
inline CGTO::Shell as_cgto_shell(const ShellT& cg)
{
    CGTO::Shell s;
    s.ang   = cg.ang;
    s.nprim = static_cast<int>(cg.alpha.size());
    s.alpha = cg.alpha;
    s.coeff = cg.coeff;
    return s;
}

// AO-within-shell -> STO_CGTO scalar type (tblite ordering):
//   s -> 0 ; p -> [py=2, pz=3, px=1] ; d -> -1 (use sphericalOverlapBlock).
inline int ao_to_type(int ang, int local_ao)
{
    if (ang == 0) return 0;
    if (ang == 1) {
        static const int p_map[3] = {2, 3, 1};  // py, pz, px
        return p_map[local_ao];
    }
    return -1;  // d not representable as a single cartesian type
}

// ===========================================================================
// X-I1: cartesian -> spherical d transform (tblite integral/trafo.f90).
//
// Cartesian d order : 0=xx 1=yy 2=zz 3=xy 4=xz 5=yz
// Spherical d order : 0=m(-2) 1=m(-1) 2=m(0) 3=m(+1) 4=m(+2)   [tblite -l..+l]
//   m=-2 : sqrt(3) * xy
//   m=-1 : sqrt(3) * yz
//   m= 0 : zz - 1/2 (xx + yy)
//   m=+1 : sqrt(3) * xz
//   m=+2 : sqrt(3)/2 (xx - yy)
// Operates on cartesian integrals computed with the single z^l-normalized
// CGTO coefficient from CGTO::slater_to_gauss (so the spherical functions are
// unit-normalized: e.g. <sqrt3 xy|sqrt3 xy> = 3 * 1/3 = 1).
// ===========================================================================
namespace dsh {
inline constexpr double s3   = 1.7320508075688772935; // sqrt(3)
inline constexpr double s3_4 = 0.8660254037844386468; // sqrt(3)/2
inline constexpr double dtrafo[5][6] = {
    //  xx     yy    zz   xy   xz   yz
    {  0.0,   0.0,  0.0,  s3, 0.0, 0.0 },  // m=-2
    {  0.0,   0.0,  0.0, 0.0, 0.0,  s3 },  // m=-1
    { -0.5,  -0.5,  1.0, 0.0, 0.0, 0.0 },  // m= 0
    {  0.0,   0.0,  0.0, 0.0,  s3, 0.0 },  // m=+1
    { s3_4, -s3_4,  0.0, 0.0, 0.0, 0.0 },  // m=+2
};
} // namespace dsh

// Number of cartesian / spherical functions for a shell.
inline int cartCount(int ang) { return (ang == 0) ? 1 : (ang == 1) ? 3 : 6; }
inline int sphCount(int ang)  { return 2 * ang + 1; }

// Cartesian powers (lx,ly,lz) for cartesian index ic of a shell.
//   p cartesian order: 0=x 1=y 2=z
//   d cartesian order: 0=xx 1=yy 2=zz 3=xy 4=xz 5=yz
inline void cartPower(int ang, int ic, int& lx, int& ly, int& lz)
{
    if (ang == 0) { lx = 0; ly = 0; lz = 0; return; }
    if (ang == 1) {
        lx = (ic == 0); ly = (ic == 1); lz = (ic == 2); return;
    }
    static const int L[6][3] = {{2,0,0},{0,2,0},{0,0,2},{1,1,0},{1,0,1},{0,1,1}};
    lx = L[ic][0]; ly = L[ic][1]; lz = L[ic][2];
}

// Transform coefficient T(isph, icart): spherical AO = sum_icart T * cartesian.
//   s : 1x1 identity
//   p : tblite spherical order (py, pz, px) -> cartesian (y, z, x)
//   d : dtrafo
inline double sphCartCoeff(int ang, int isph, int icart)
{
    if (ang == 0) return 1.0;
    if (ang == 1) {
        static const int p_cart[3] = {1, 2, 0};  // py->y, pz->z, px->x
        return (p_cart[isph] == icart) ? 1.0 : 0.0;
    }
    return dsh::dtrafo[isph][icart];
}

// General 1-D Gaussian moment over cartesian powers la,lb in {0,1,2,3}, moment
// order n in {0,1,2}:  I = int (u+PA)^la (u+PB)^lb (u+P)^n exp(-gamma u^2) du.
// la up to 3 is required for the gradient raise/lower (d -> d+1). Separate from
// multipole_ints::moment1d (which stays byte-identical for the s/p path); this
// one is only used on d-touching pairs.
inline double cartMoment1d(int la, int lb, int n,
                           double PA, double PB, double P, double gamma)
{
    // Polynomial coeffs of (u+X)^l (l<=3):  c[k] = C(l,k) * X^(l-k).
    auto poly = [](int l, double X, double* c) {
        static const int Bn[4][4] = {{1,0,0,0},{1,1,0,0},{1,2,1,0},{1,3,3,1}};
        const double xp[4] = {1.0, X, X * X, X * X * X};
        for (int k = 0; k <= l; ++k) c[k] = Bn[l][k] * xp[l - k];
    };
    double a[4] = {0,0,0,0}, b[4] = {0,0,0,0}, p[4] = {0,0,0,0};
    poly(la, PA, a); poly(lb, PB, b); poly(n, P, p);

    double ab[7] = {0,0,0,0,0,0,0};
    for (int i = 0; i <= la; ++i)
        for (int j = 0; j <= lb; ++j) ab[i + j] += a[i] * b[j];
    double c[10] = {0,0,0,0,0,0,0,0,0,0};
    for (int i = 0; i <= la + lb; ++i)
        for (int j = 0; j <= n; ++j) c[i + j] += ab[i] * p[j];

    // Even Gaussian moments: M_{2k} = (2k-1)!! / (2 gamma)^k * sqrt(pi/gamma).
    const double M0 = std::sqrt(M_PI / gamma);
    const double w  = 1.0 / (2.0 * gamma);
    static const double dfac[5] = {1.0, 1.0, 3.0, 15.0, 105.0};  // (2k-1)!!
    double sum = 0.0, wk = 1.0;
    for (int k = 0; k <= 4; ++k) { sum += c[2 * k] * (M0 * dfac[k] * wk); wk *= w; }
    return sum;
}

// Contracted cartesian overlap for explicit powers (la..) on shells A, B.
inline double cgtoCartOverlap(const CGTO::Shell& A, const CGTO::Shell& B,
                              double ax, double ay, double az,
                              double bx, double by, double bz,
                              int lxa, int lya, int lza,
                              int lxb, int lyb, int lzb)
{
    const double dx = bx - ax, dy = by - ay, dz = bz - az;
    const double R2 = dx * dx + dy * dy + dz * dz;
    double S = 0.0;
    for (int i = 0; i < A.nprim; ++i) {
        const double ai = A.alpha[i], ci = A.coeff[i];
        for (int j = 0; j < B.nprim; ++j) {
            const double aj = B.alpha[j], cj = B.coeff[j];
            const double g  = ai + aj;
            const double Px = (ai * ax + aj * bx) / g;
            const double Py = (ai * ay + aj * by) / g;
            const double Pz = (ai * az + aj * bz) / g;
            const double K  = std::exp(-ai * aj / g * R2);
            const double Mx = cartMoment1d(lxa, lxb, 0, Px - ax, Px - bx, Px, g);
            const double My = cartMoment1d(lya, lyb, 0, Py - ay, Py - by, Py, g);
            const double Mz = cartMoment1d(lza, lzb, 0, Pz - az, Pz - bz, Pz, g);
            S += ci * cj * K * Mx * My * Mz;
        }
    }
    return S;
}

// Spherical overlap block for a shell pair (ang_a x ang_b) for any combination
// involving d. Writes out[isph_a*ld + isph_b] = T_a . S_cart . T_b^T.
// out must hold at least sphCount(ang_a) rows of stride ld (>= sphCount(ang_b)).
inline void sphericalOverlapBlock(const CGTO::Shell& A, int ang_a,
                                  const CGTO::Shell& B, int ang_b,
                                  double ax, double ay, double az,
                                  double bx, double by, double bz,
                                  double* out, int ld)
{
    const int nca = cartCount(ang_a), ncb = cartCount(ang_b);
    const int nsa = sphCount(ang_a),  nsb = sphCount(ang_b);

    double cart[6 * 6];
    for (int ca = 0; ca < nca; ++ca) {
        int lxa, lya, lza; cartPower(ang_a, ca, lxa, lya, lza);
        for (int cb = 0; cb < ncb; ++cb) {
            int lxb, lyb, lzb; cartPower(ang_b, cb, lxb, lyb, lzb);
            cart[ca * 6 + cb] = cgtoCartOverlap(A, B, ax, ay, az, bx, by, bz,
                                                lxa, lya, lza, lxb, lyb, lzb);
        }
    }
    for (int sa = 0; sa < nsa; ++sa) {
        for (int sb = 0; sb < nsb; ++sb) {
            double v = 0.0;
            for (int ca = 0; ca < nca; ++ca) {
                const double ta = sphCartCoeff(ang_a, sa, ca);
                if (ta == 0.0) continue;
                for (int cb = 0; cb < ncb; ++cb) {
                    const double tb = sphCartCoeff(ang_b, sb, cb);
                    if (tb == 0.0) continue;
                    v += ta * tb * cart[ca * 6 + cb];
                }
            }
            out[sa * ld + sb] = v;
        }
    }
}

// Contracted cartesian dipole(3) + raw-cartesian quadrupole(6) for explicit
// powers, global moment origin (origin shift + traceless done by the caller).
// Q packing: 0=xx 1=xy 2=yy 3=xz 4=yz 5=zz  (matches xtb_multipole_ints qp order).
inline void cgtoCartMultipole(const CGTO::Shell& A, const CGTO::Shell& B,
                              double ax, double ay, double az,
                              double bx, double by, double bz,
                              int lxa, int lya, int lza,
                              int lxb, int lyb, int lzb,
                              double D[3], double Q[6])
{
    const double dx = bx - ax, dy = by - ay, dz = bz - az;
    const double R2 = dx * dx + dy * dy + dz * dz;
    for (int k = 0; k < 3; ++k) D[k] = 0.0;
    for (int k = 0; k < 6; ++k) Q[k] = 0.0;
    for (int i = 0; i < A.nprim; ++i) {
        const double ai = A.alpha[i], ci = A.coeff[i];
        for (int j = 0; j < B.nprim; ++j) {
            const double aj = B.alpha[j], cj = B.coeff[j];
            const double g  = ai + aj;
            const double Px = (ai * ax + aj * bx) / g;
            const double Py = (ai * ay + aj * by) / g;
            const double Pz = (ai * az + aj * bz) / g;
            const double K  = std::exp(-ai * aj / g * R2);
            const double kc = ci * cj * K;
            const double Mx0 = cartMoment1d(lxa, lxb, 0, Px - ax, Px - bx, Px, g);
            const double My0 = cartMoment1d(lya, lyb, 0, Py - ay, Py - by, Py, g);
            const double Mz0 = cartMoment1d(lza, lzb, 0, Pz - az, Pz - bz, Pz, g);
            const double Mx1 = cartMoment1d(lxa, lxb, 1, Px - ax, Px - bx, Px, g);
            const double My1 = cartMoment1d(lya, lyb, 1, Py - ay, Py - by, Py, g);
            const double Mz1 = cartMoment1d(lza, lzb, 1, Pz - az, Pz - bz, Pz, g);
            const double Mx2 = cartMoment1d(lxa, lxb, 2, Px - ax, Px - bx, Px, g);
            const double My2 = cartMoment1d(lya, lyb, 2, Py - ay, Py - by, Py, g);
            const double Mz2 = cartMoment1d(lza, lzb, 2, Pz - az, Pz - bz, Pz, g);
            D[0] += kc * Mx1 * My0 * Mz0;
            D[1] += kc * Mx0 * My1 * Mz0;
            D[2] += kc * Mx0 * My0 * Mz1;
            Q[0] += kc * Mx2 * My0 * Mz0;  // xx
            Q[1] += kc * Mx1 * My1 * Mz0;  // xy
            Q[2] += kc * Mx0 * My2 * Mz0;  // yy
            Q[3] += kc * Mx1 * My0 * Mz1;  // xz
            Q[4] += kc * Mx0 * My1 * Mz1;  // yz
            Q[5] += kc * Mx0 * My0 * Mz2;  // zz
        }
    }
}

// Spherical multipole block (global origin) for a d-touching shell pair.
// Writes per spherical AO pair the dipole Dout[(sa*ld+sb)*3 + k] and the raw
// cartesian quadrupole Qout[(sa*ld+sb)*6 + k]. The operator components k are
// cartesian and NOT transformed; only the bra/ket AO indices are (via dtrafo).
inline void sphericalMultipoleBlock(const CGTO::Shell& A, int ang_a,
                                    const CGTO::Shell& B, int ang_b,
                                    double ax, double ay, double az,
                                    double bx, double by, double bz,
                                    double* Dout, double* Qout, int ld)
{
    const int nca = cartCount(ang_a), ncb = cartCount(ang_b);
    const int nsa = sphCount(ang_a),  nsb = sphCount(ang_b);
    double cD[36][3], cQ[36][6];
    for (int ca = 0; ca < nca; ++ca) {
        int lxa, lya, lza; cartPower(ang_a, ca, lxa, lya, lza);
        for (int cb = 0; cb < ncb; ++cb) {
            int lxb, lyb, lzb; cartPower(ang_b, cb, lxb, lyb, lzb);
            cgtoCartMultipole(A, B, ax, ay, az, bx, by, bz,
                              lxa, lya, lza, lxb, lyb, lzb,
                              cD[ca * 6 + cb], cQ[ca * 6 + cb]);
        }
    }
    for (int sa = 0; sa < nsa; ++sa) {
        for (int sb = 0; sb < nsb; ++sb) {
            double d[3] = {0, 0, 0}, q[6] = {0, 0, 0, 0, 0, 0};
            for (int ca = 0; ca < nca; ++ca) {
                const double ta = sphCartCoeff(ang_a, sa, ca);
                if (ta == 0.0) continue;
                for (int cb = 0; cb < ncb; ++cb) {
                    const double tb = sphCartCoeff(ang_b, sb, cb);
                    if (tb == 0.0) continue;
                    const double w = ta * tb;
                    const int idx = ca * 6 + cb;
                    for (int k = 0; k < 3; ++k) d[k] += w * cD[idx][k];
                    for (int k = 0; k < 6; ++k) q[k] += w * cQ[idx][k];
                }
            }
            const int o = sa * ld + sb;
            for (int k = 0; k < 3; ++k) Dout[o * 3 + k] = d[k];
            for (int k = 0; k < 6; ++k) Qout[o * 6 + k] = q[k];
        }
    }
}

// ===========================================================================
// X-I1 gradients (B5). All blocks return the A-derivative (d/dR of the atom
// carrying shell A); the H0/Pulay loop uses Newton's 3rd law for the B atom.
// ===========================================================================

// Gradient of the contracted cartesian overlap wrt center A (Obara-Saika
// raise/lower):  dS/dA_k = sum_ij ci cj [2 ai S(a+1_k;b) - l_{a,k} S(a-1_k;b)].
inline void cgtoCartOverlapGrad(const CGTO::Shell& A, const CGTO::Shell& B,
                                double ax, double ay, double az,
                                double bx, double by, double bz,
                                int lxa, int lya, int lza,
                                int lxb, int lyb, int lzb,
                                double gA[3])
{
    gA[0] = gA[1] = gA[2] = 0.0;
    const double dx = bx - ax, dy = by - ay, dz = bz - az;
    const double R2 = dx * dx + dy * dy + dz * dz;
    for (int i = 0; i < A.nprim; ++i) {
        const double ai = A.alpha[i], ci = A.coeff[i];
        for (int j = 0; j < B.nprim; ++j) {
            const double aj = B.alpha[j], cj = B.coeff[j];
            const double g  = ai + aj, cc = ci * cj;
            const double Px = (ai*ax + aj*bx) / g, Py = (ai*ay + aj*by) / g, Pz = (ai*az + aj*bz) / g;
            const double K  = std::exp(-ai * aj / g * R2);
            const double PAx = Px-ax, PBx = Px-bx, PAy = Py-ay, PBy = Py-by, PAz = Pz-az, PBz = Pz-bz;
            const double Mx = cartMoment1d(lxa, lxb, 0, PAx, PBx, Px, g);
            const double My = cartMoment1d(lya, lyb, 0, PAy, PBy, Py, g);
            const double Mz = cartMoment1d(lza, lzb, 0, PAz, PBz, Pz, g);
            const double Mxp = cartMoment1d(lxa+1, lxb, 0, PAx, PBx, Px, g);
            const double Myp = cartMoment1d(lya+1, lyb, 0, PAy, PBy, Py, g);
            const double Mzp = cartMoment1d(lza+1, lzb, 0, PAz, PBz, Pz, g);
            const double Mxm = (lxa > 0) ? cartMoment1d(lxa-1, lxb, 0, PAx, PBx, Px, g) : 0.0;
            const double Mym = (lya > 0) ? cartMoment1d(lya-1, lyb, 0, PAy, PBy, Py, g) : 0.0;
            const double Mzm = (lza > 0) ? cartMoment1d(lza-1, lzb, 0, PAz, PBz, Pz, g) : 0.0;
            gA[0] += cc * K * My * Mz * (2.0*ai*Mxp - lxa*Mxm);
            gA[1] += cc * K * Mx * Mz * (2.0*ai*Myp - lya*Mym);
            gA[2] += cc * K * Mx * My * (2.0*ai*Mzp - lza*Mzm);
        }
    }
}

// dS_munu/dA[3] for each spherical AO pair of a d-touching shell pair.
// out[(sa*ld+sb)*3 + l] = dS/dA_l.  (dS/dB_l = -dS/dA_l.)
inline void sphericalOverlapGradBlock(const CGTO::Shell& A, int ang_a,
                                      const CGTO::Shell& B, int ang_b,
                                      double ax, double ay, double az,
                                      double bx, double by, double bz,
                                      double* out, int ld)
{
    const int nca = cartCount(ang_a), ncb = cartCount(ang_b);
    const int nsa = sphCount(ang_a),  nsb = sphCount(ang_b);
    double cg[36][3];
    for (int ca = 0; ca < nca; ++ca) {
        int lxa, lya, lza; cartPower(ang_a, ca, lxa, lya, lza);
        for (int cb = 0; cb < ncb; ++cb) {
            int lxb, lyb, lzb; cartPower(ang_b, cb, lxb, lyb, lzb);
            cgtoCartOverlapGrad(A, B, ax, ay, az, bx, by, bz,
                                lxa, lya, lza, lxb, lyb, lzb, cg[ca*6+cb]);
        }
    }
    for (int sa = 0; sa < nsa; ++sa)
        for (int sb = 0; sb < nsb; ++sb) {
            double g[3] = {0, 0, 0};
            for (int ca = 0; ca < nca; ++ca) {
                const double ta = sphCartCoeff(ang_a, sa, ca);
                if (ta == 0.0) continue;
                for (int cb = 0; cb < ncb; ++cb) {
                    const double tb = sphCartCoeff(ang_b, sb, cb);
                    if (tb == 0.0) continue;
                    const double w = ta * tb; const int idx = ca*6+cb;
                    for (int l = 0; l < 3; ++l) g[l] += w * cg[idx][l];
                }
            }
            const int o = (sa*ld + sb) * 3;
            out[o] = g[0]; out[o+1] = g[1]; out[o+2] = g[2];
        }
}

// Global-origin cartesian multipole + A-gradient for explicit powers (port of
// multipole_ints::primitive_multipole_grad, A-derivatives only, summed over
// primitives). Q packing 0=xx 1=xy 2=yy 3=xz 4=yz 5=zz.
inline void cgtoCartMultipoleGradGlobal(const CGTO::Shell& A, const CGTO::Shell& B,
                                        double ax, double ay, double az,
                                        double bx, double by, double bz,
                                        int lxa, int lya, int lza,
                                        int lxb, int lyb, int lzb,
                                        double& S, double D[3], double Q[6],
                                        double dS_dA[3], double dDg_dA[3][3],
                                        double dQg_dA[3][6])
{
    S = 0.0;
    for (int k = 0; k < 3; ++k) { D[k] = 0.0; dS_dA[k] = 0.0; }
    for (int k = 0; k < 6; ++k) Q[k] = 0.0;
    for (int l = 0; l < 3; ++l) { for (int k = 0; k < 3; ++k) dDg_dA[l][k] = 0.0;
                                  for (int q = 0; q < 6; ++q) dQg_dA[l][q] = 0.0; }
    static const int qnx[6] = {2,1,0,1,0,0};
    static const int qny[6] = {0,1,2,0,1,0};
    static const int qnz[6] = {0,0,0,1,1,2};
    const double Rx = ax - bx, Ry = ay - by, Rz = az - bz;
    const double R2 = Rx*Rx + Ry*Ry + Rz*Rz;
    for (int ip = 0; ip < A.nprim; ++ip) {
        const double ai = A.alpha[ip], ci = A.coeff[ip];
        for (int jp = 0; jp < B.nprim; ++jp) {
            const double aj = B.alpha[jp], cj = B.coeff[jp];
            const double g = ai + aj, cc = ci * cj;
            const double Px = (ai*ax + aj*bx) / g, Py = (ai*ay + aj*by) / g, Pz = (ai*az + aj*bz) / g;
            const double PAx = Px-ax, PBx = Px-bx, PAy = Py-ay, PBy = Py-by, PAz = Pz-az, PBz = Pz-bz;
            const double K = std::exp(-ai * aj / g * R2);
            double Ix[3], Iy[3], Iz[3], IxAlo[3] = {}, IyAlo[3] = {}, IzAlo[3] = {},
                   IxBlo[3] = {}, IyBlo[3] = {}, IzBlo[3] = {};
            for (int n = 0; n <= 2; ++n) {
                Ix[n] = cartMoment1d(lxa, lxb, n, PAx, PBx, Px, g);
                Iy[n] = cartMoment1d(lya, lyb, n, PAy, PBy, Py, g);
                Iz[n] = cartMoment1d(lza, lzb, n, PAz, PBz, Pz, g);
            }
            if (lxa > 0) for (int n=0;n<=2;++n) IxAlo[n] = cartMoment1d(lxa-1, lxb, n, PAx, PBx, Px, g);
            if (lya > 0) for (int n=0;n<=2;++n) IyAlo[n] = cartMoment1d(lya-1, lyb, n, PAy, PBy, Py, g);
            if (lza > 0) for (int n=0;n<=2;++n) IzAlo[n] = cartMoment1d(lza-1, lzb, n, PAz, PBz, Pz, g);
            if (lxb > 0) for (int n=0;n<=2;++n) IxBlo[n] = cartMoment1d(lxa, lxb-1, n, PAx, PBx, Px, g);
            if (lyb > 0) for (int n=0;n<=2;++n) IyBlo[n] = cartMoment1d(lya, lyb-1, n, PAy, PBy, Py, g);
            if (lzb > 0) for (int n=0;n<=2;++n) IzBlo[n] = cartMoment1d(lza, lzb-1, n, PAz, PBz, Pz, g);

            S    += cc * K * Ix[0]*Iy[0]*Iz[0];
            D[0] += cc * K * Ix[1]*Iy[0]*Iz[0];
            D[1] += cc * K * Ix[0]*Iy[1]*Iz[0];
            D[2] += cc * K * Ix[0]*Iy[0]*Iz[1];
            Q[0] += cc * K * Ix[2]*Iy[0]*Iz[0];
            Q[1] += cc * K * Ix[1]*Iy[1]*Iz[0];
            Q[2] += cc * K * Ix[0]*Iy[2]*Iz[0];
            Q[3] += cc * K * Ix[1]*Iy[0]*Iz[1];
            Q[4] += cc * K * Ix[0]*Iy[1]*Iz[1];
            Q[5] += cc * K * Ix[0]*Iy[0]*Iz[2];

            const double kfA[3] = { -2.0*ai*aj/g*Rx, -2.0*ai*aj/g*Ry, -2.0*ai*aj/g*Rz };
            auto dG_dA = [&](int nx, int ny, int nz, int l) -> double {
                const double base = Ix[nx]*Iy[ny]*Iz[nz];
                double r = K * kfA[l] * base;
                if (l == 0) {
                    double t = (lxa>0 ? lxa*(-aj/g)*IxAlo[nx] : 0.0)
                             + (lxb>0 ? lxb*( ai/g)*IxBlo[nx] : 0.0)
                             + (nx >0 ? nx *( ai/g)*Ix[nx-1] : 0.0);
                    r += K * t * Iy[ny] * Iz[nz];
                } else if (l == 1) {
                    double t = (lya>0 ? lya*(-aj/g)*IyAlo[ny] : 0.0)
                             + (lyb>0 ? lyb*( ai/g)*IyBlo[ny] : 0.0)
                             + (ny >0 ? ny *( ai/g)*Iy[ny-1] : 0.0);
                    r += K * Ix[nx] * t * Iz[nz];
                } else {
                    double t = (lza>0 ? lza*(-aj/g)*IzAlo[nz] : 0.0)
                             + (lzb>0 ? lzb*( ai/g)*IzBlo[nz] : 0.0)
                             + (nz >0 ? nz *( ai/g)*Iz[nz-1] : 0.0);
                    r += K * Ix[nx] * Iy[ny] * t;
                }
                return r;
            };
            for (int l = 0; l < 3; ++l) {
                dS_dA[l] += cc * dG_dA(0, 0, 0, l);
                for (int k = 0; k < 3; ++k) {
                    const int dnx = (k==0), dny = (k==1), dnz = (k==2);
                    dDg_dA[l][k] += cc * dG_dA(dnx, dny, dnz, l);
                }
                for (int q = 0; q < 6; ++q)
                    dQg_dA[l][q] += cc * dG_dA(qnx[q], qny[q], qnz[q], l);
            }
        }
    }
}

// A-gradient of the TRANSFORMED (origin at B, traceless) dipole/quadrupole
// integrals for a d-touching shell pair, per spherical AO pair.
//   dDout[(sa*ld+sb)*3 + k] = d(dp_int[k])/dA_l   (k=dipole comp, packed l outer)
//   dQout[(sa*ld+sb)*6 + q] = d(qp_int[q])/dA_l
// Matches cgto_multipole_grad_transformed's dD_dA / dQ_dA convention.
inline void sphericalMultipoleGradBlock(const CGTO::Shell& A, int ang_a,
                                        const CGTO::Shell& B, int ang_b,
                                        double ax, double ay, double az,
                                        double bx, double by, double bz,
                                        double* dDout, double* dQout, int ld)
{
    const int nca = cartCount(ang_a), ncb = cartCount(ang_b);
    const int nsa = sphCount(ang_a),  nsb = sphCount(ang_b);
    const double Bv[3] = {bx, by, bz};
    static const int qa6[6] = {0,0,1,0,1,2};
    static const int qb6[6] = {0,1,1,2,2,2};

    // Cartesian global integrals + A-gradients per cartesian pair.
    double cS[36], cDg[36][3][3], cQg[36][3][6], cdS[36][3];
    double cDv[36][3];   // global dipole value (for the origin-shift correction)
    for (int ca = 0; ca < nca; ++ca) {
        int lxa, lya, lza; cartPower(ang_a, ca, lxa, lya, lza);
        for (int cb = 0; cb < ncb; ++cb) {
            int lxb, lyb, lzb; cartPower(ang_b, cb, lxb, lyb, lzb);
            double S, D[3], Q[6], dS[3], dDg[3][3], dQg[3][6];
            cgtoCartMultipoleGradGlobal(A, B, ax, ay, az, bx, by, bz,
                                        lxa, lya, lza, lxb, lyb, lzb,
                                        S, D, Q, dS, dDg, dQg);
            const int idx = ca*6+cb;
            cS[idx] = S;
            for (int k=0;k<3;++k){ cDv[idx][k]=D[k]; cdS[idx][k]=dS[k];
                                   for(int l=0;l<3;++l) cDg[idx][l][k]=dDg[l][k]; }
            for (int l=0;l<3;++l) for (int q=0;q<6;++q) cQg[idx][l][q]=dQg[l][q];
        }
    }

    for (int sa = 0; sa < nsa; ++sa) {
        for (int sb = 0; sb < nsb; ++sb) {
            // Transform cartesian -> spherical (operator components untouched).
            double S_raw = 0.0, D_raw[3] = {}, dS_dA[3] = {}, dDg_dA[3][3] = {}, dQg_dA[3][6] = {};
            for (int ca = 0; ca < nca; ++ca) {
                const double ta = sphCartCoeff(ang_a, sa, ca);
                if (ta == 0.0) continue;
                for (int cb = 0; cb < ncb; ++cb) {
                    const double tb = sphCartCoeff(ang_b, sb, cb);
                    if (tb == 0.0) continue;
                    const double w = ta*tb; const int idx = ca*6+cb;
                    S_raw += w * cS[idx];
                    for (int k=0;k<3;++k){ D_raw[k] += w*cDv[idx][k]; dS_dA[k] += w*cdS[idx][k];
                        for(int l=0;l<3;++l) dDg_dA[l][k] += w*cDg[idx][l][k]; }
                    for (int l=0;l<3;++l) for (int q=0;q<6;++q) dQg_dA[l][q] += w*cQg[idx][l][q];
                }
            }
            // Origin shift to B + traceless (A-derivatives only).
            // Layout: dDout[(sa*ld+sb)*9 + l*3 + k], dQout[(sa*ld+sb)*18 + l*6 + q].
            (void)S_raw; (void)D_raw;
            const int oD = (sa*ld + sb) * 9;
            const int oQ = (sa*ld + sb) * 18;
            for (int l = 0; l < 3; ++l) {
                for (int k = 0; k < 3; ++k)
                    dDout[oD + l*3 + k] = dDg_dA[l][k] - Bv[k] * dS_dA[l];
                double dqraw[6];
                for (int q = 0; q < 6; ++q) {
                    const int a = qa6[q], b = qb6[q];
                    dqraw[q] = dQg_dA[l][q] - Bv[a]*dDg_dA[l][b] - Bv[b]*dDg_dA[l][a]
                             + Bv[a]*Bv[b]*dS_dA[l];
                }
                const double dtr = 0.5 * (dqraw[0] + dqraw[2] + dqraw[5]);
                double* dq = &dQout[oQ + l*6];
                dq[0] = 1.5*dqraw[0] - dtr;
                dq[1] = 1.5*dqraw[1];
                dq[2] = 1.5*dqraw[2] - dtr;
                dq[3] = 1.5*dqraw[3];
                dq[4] = 1.5*dqraw[4];
                dq[5] = 1.5*dqraw[5] - dtr;
            }
        }
    }
}

} // namespace curcuma::xtb
