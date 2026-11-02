/*
 * <Native QM: spin-polarised exchange-correlation point kernels (UKS), shared by the CPU (qm_xc_spin.cpp) and the
 *  GPU (cuda/qm_gpu_xc.cu)>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * ONE source of the spin-polarised formulas for both sides, as qm_xc_kernels.h is for the closed-shell ones: the CPU
 * wrapper qmxc::evaluatePointSpin (qm_xc_spin.cpp) and the device kernels of cuda/qm_gpu_xc.cu call
 * qmxc::kern::evaluatePointSpinTerms below. Formulas and literature: qm_xc_spin.h. Every correlation kernel was first
 * written as a short Python prototype and compared with libxc 7.0.0 (pyscf.dft.libxc.eval_xc, spin = 1); ctest
 * qm_uks_kernel repeats the comparison for eps, v_rho and v_sigma on a fixed point set, ctest qm_gpu_uhf compares the
 * device V_xc / E_xc / XC gradient with the CPU.
 *
 * Derivatives: forward-mode dual numbers in the five variables (rho_a, rho_b, sigma_aa, sigma_ab, sigma_bb), so the
 * potential is the exact derivative of the energy formula written here (no hand-derived v formulas to get wrong).
 *
 * Only plain C math from <cmath> in the global namespace (::pow, ::exp, ::log, ::sqrt, ::atan, ::cbrt), no std
 * containers: compiles as host code (g++) and as __host__ __device__ code (nvcc). Device math may differ from glibc
 * by an ulp, so the GPU agrees with the CPU to rounding, not bitwise.
 *
 * Claude Generated (Nov 2026, stage 6a: moved here unchanged from qm_xc_spin.cpp of stage 5b, std:: math replaced by
 * the same global-namespace functions). AI-generated, machine-tested only.
 *
 * This program is free software under GPL-3.0
 */

#pragma once

#include "qm_xc_kernels.h"

#include <cmath>

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Winline"
#endif

namespace qmxc {
namespace kern {
namespace spin {

// ---------------------------------------------------------------------------------------------------------------
// Dual numbers in five variables: value v and partial derivatives d[k] = dv / dx_k,
// x = (rho_a, rho_b, sigma_aa, sigma_ab, sigma_bb).
// ---------------------------------------------------------------------------------------------------------------
struct D5 {
    double v;
    double d[5];
};

QM_HD D5 cst(double c)
{
    D5 r;
    r.v = c;
    for (int k = 0; k < 5; ++k) r.d[k] = 0.0;
    return r;
}
QM_HD D5 var(double x, int k)
{
    D5 r = cst(x);
    r.d[k] = 1.0;
    return r;
}
/// f(a) with value fv and derivative fp = f'(a.v): chain rule.
QM_HD D5 chain(const D5& a, double fv, double fp)
{
    D5 r;
    r.v = fv;
    for (int k = 0; k < 5; ++k) r.d[k] = fp * a.d[k];
    return r;
}
QM_HD D5 operator+(const D5& a, const D5& b)
{
    D5 r;
    r.v = a.v + b.v;
    for (int k = 0; k < 5; ++k) r.d[k] = a.d[k] + b.d[k];
    return r;
}
QM_HD D5 operator-(const D5& a, const D5& b)
{
    D5 r;
    r.v = a.v - b.v;
    for (int k = 0; k < 5; ++k) r.d[k] = a.d[k] - b.d[k];
    return r;
}
QM_HD D5 operator-(const D5& a) { return chain(a, -a.v, -1.0); }
QM_HD D5 operator*(const D5& a, const D5& b)
{
    D5 r;
    r.v = a.v * b.v;
    for (int k = 0; k < 5; ++k) r.d[k] = a.d[k] * b.v + a.v * b.d[k];
    return r;
}
QM_HD D5 operator/(const D5& a, const D5& b)
{
    D5 r;
    r.v = a.v / b.v;
    const double ib2 = 1.0 / (b.v * b.v);
    for (int k = 0; k < 5; ++k) r.d[k] = (a.d[k] * b.v - a.v * b.d[k]) * ib2;
    return r;
}
QM_HD D5 operator+(const D5& a, double c) { return chain(a, a.v + c, 1.0); }
QM_HD D5 operator+(double c, const D5& a) { return chain(a, a.v + c, 1.0); }
QM_HD D5 operator-(const D5& a, double c) { return chain(a, a.v - c, 1.0); }
QM_HD D5 operator-(double c, const D5& a) { return chain(a, c - a.v, -1.0); }
QM_HD D5 operator*(double c, const D5& a) { return chain(a, a.v * c, c); }
QM_HD D5 operator/(const D5& a, double c) { return chain(a, a.v / c, 1.0 / c); }
QM_HD D5 operator/(double c, const D5& a) { return chain(a, c / a.v, -c / (a.v * a.v)); }
QM_HD D5 dexp(const D5& a)
{
    const double e = ::exp(a.v);
    return chain(a, e, e);
}
QM_HD D5 dlog(const D5& a) { return chain(a, ::log(a.v), 1.0 / a.v); }
QM_HD D5 dsqrt(const D5& a)
{
    const double s = ::sqrt(a.v);
    return chain(a, s, 0.5 / s);
}
QM_HD D5 datan(const D5& a) { return chain(a, ::atan(a.v), 1.0 / (1.0 + a.v * a.v)); }
/// a^p for a > 0 (p real); a == 0 with p > 1 gives 0 and a zero derivative.
QM_HD D5 dpow(const D5& a, double p)
{
    if (a.v == 0.0 && p > 1.0) return chain(a, 0.0, 0.0);
    const double vp = ::pow(a.v, p);
    return chain(a, vp, p * vp / a.v);
}

/// libxc's zeta_threshold (DBL_EPSILON): 1 + zeta and 1 - zeta are clamped to it (libxc my_piecewise in f_zeta, phi).
constexpr double kZetaThreshold = 2.220446049250313e-16;
/// (1 +- zeta)^p with libxc's clamp: below the threshold the constant threshold^p (no derivative).
QM_HD D5 opzPow(const D5& x, double p)
{
    if (x.v <= kZetaThreshold) return cst(::pow(kZetaThreshold, p));
    return dpow(x, p);
}

/// Spin interpolation f(z) = ((1+z)^(4/3) + (1-z)^(4/3) - 2) / (2^(4/3) - 2) (VWN 1980 eq. 2.4; PW92 eq. 9).
QM_HD D5 fZeta(const D5& z)
{
    return (opzPow(1.0 + z, 4.0 / 3.0) + opzPow(1.0 - z, 4.0 / 3.0) - 2.0) / (::pow(2.0, 4.0 / 3.0) - 2.0);
}

/// VWN fit function of x = sqrt(r_s) (VWN 1980 eq. 4.4) for one parameter set.
QM_HD D5 vwnAux(double A, double b, double c, double x0, const D5& x)
{
    const D5 X = x * x + b * x + c;
    const double X0 = x0 * x0 + b * x0 + c;
    const double Q = ::sqrt(4.0 * c - b * b);
    const D5 at = datan(Q / (2.0 * x + b));
    return A * (dlog(x * x / X) + (2.0 * b / Q) * at
                - (b * x0 / X0) * (dlog((x - x0) * (x - x0) / X) + (2.0 * (b + 2.0 * x0) / Q) * at));
}

/// VWN correlation energy per electron, spin-polarised. rpa = false: libxc LDA_C_VWN (VWN5 Monte-Carlo fits,
/// interpolation with the spin stiffness, VWN 1980 eq. 4.4/ "functional V"); rpa = true: libxc LDA_C_VWN_RPA (RPA fits,
/// plain f(z) interpolation between the para- and ferromagnetic fits; the closed-shell value is the vwn_rpa kernel).
QM_HD D5 vwnSpin(const D5& rs, const D5& z, bool rpa)
{
    const D5 x = dsqrt(rs);
    const D5 fz = fZeta(z);
    if (rpa) {
        const D5 eP = vwnAux(0.0310907, 13.0720, 42.7198, -0.409286, x);
        const D5 eF = vwnAux(0.01554535, 20.1231, 101.578, -0.743294, x);
        return eP * (1.0 - fz) + eF * fz;
    }
    const D5 eP = vwnAux(0.0310907, 3.72744, 12.9352, -0.10498, x);
    const D5 eF = vwnAux(0.01554535, 7.06042, 18.0578, -0.32500, x);
    const D5 aC = vwnAux(-1.0 / (6.0 * kPi * kPi), 1.13107, 13.0045, -0.0047584, x);
    const double fpp = 4.0 / (9.0 * (::cbrt(2.0) - 1.0));
    const D5 z4 = z * z * z * z;
    return eP + aC * fz * (1.0 - z4) / fpp + (eF - eP) * fz * z4;
}

/// PW92 G(r_s) of one fit (PW92 eq. 10).
QM_HD D5 pw92GSpin(const D5& rs, double A, double a1, double b1, double b2, double b3, double b4)
{
    const D5 srs = dsqrt(rs);
    const D5 Q1 = 2.0 * A * (b1 * srs + b2 * rs + b3 * rs * srs + b4 * rs * rs);
    return -2.0 * A * (1.0 + a1 * rs) * dlog(1.0 + 1.0 / Q1);
}

/// PW92 correlation per electron, spin-polarised (PW92 eq. 8), libxc LDA_C_PW_MOD constants and f''(0):
///   eps = G1 + z^4 f(z) (G2 - G1 + G3/f''(0)) - f(z) G3/f''(0)    (G3 = -alpha_c, libxc's form)
QM_HD D5 pw92Spin(const D5& rs, const D5& z)
{
    const D5 g1 = pw92GSpin(rs, 0.0310907, 0.21370, 7.5957, 3.5876, 1.6382, 0.49294);
    const D5 g2 = pw92GSpin(rs, 0.01554535, 0.20548, 14.1189, 6.1977, 3.3662, 0.62517);
    const D5 g3 = pw92GSpin(rs, 0.0168869, 0.11125, 10.357, 3.6231, 0.88026, 0.49671);
    const double fz20 = 1.709920934161365617563962776245;
    const D5 fz = fZeta(z);
    const D5 z4 = z * z * z * z;
    return g1 + z4 * fz * (g2 - g1 + g3 / fz20) - fz * g3 / fz20;
}

/// PBE correlation per electron, spin-polarised (PBE 1996 eqs. 3, 7, 8): eps_PW92(r_s, z) + H(r_s, z, t),
///   phi = ((1+z)^(2/3) + (1-z)^(2/3))/2, t^2 = sigma / (4 phi^2 k_s^2 rho^2), k_s = sqrt(4 k_F / pi).
QM_HD D5 pbeCSpin(const D5& rho, const D5& rs, const D5& z, const D5& sigma)
{
    const double beta = 0.06672455060314922, gam = (1.0 - ::log(2.0)) / (kPi * kPi);
    const D5 ec = pw92Spin(rs, z);
    const D5 phi = 0.5 * (opzPow(1.0 + z, 2.0 / 3.0) + opzPow(1.0 - z, 2.0 / 3.0));
    const D5 phi3 = phi * phi * phi;
    const D5 kf = dpow(3.0 * kPi * kPi * rho, 1.0 / 3.0);
    const D5 ks2 = 4.0 * kf / kPi;
    const D5 t2 = sigma / (4.0 * phi * phi * ks2 * rho * rho);
    const D5 A = (beta / gam) / (dexp(-ec / (gam * phi3)) - 1.0);
    const D5 At2 = A * t2;
    const D5 H = gam * phi3 * dlog(1.0 + (beta / gam) * t2 * (1.0 + At2) / (1.0 + At2 + At2 * At2));
    return ec + H;
}

/// LYP correlation ENERGY DENSITY (per volume), Miehlich et al. 1989 eq. 2 (no Laplacian; libxc GGA_C_LYP):
///   e = -4a/(1 + d r) ra rb / rho - a b w { ra rb [ 2^(11/3) C_F (ra^(8/3) + rb^(8/3)) + (47/18 - 7/18 delta) s
///       - (5/2 - delta/18)(saa + sbb) - (delta - 11)/9 (ra/rho saa + rb/rho sbb) ]
///       - 2/3 rho^2 s + (2/3 rho^2 - ra^2) sbb + (2/3 rho^2 - rb^2) saa },
///   r = rho^(-1/3), w = exp(-c r)/(1 + d r) rho^(-11/3), delta = c r + d r/(1 + d r), s = saa + 2 sab + sbb.
QM_HD D5 lypSpinDensity(const D5& ra, const D5& rb, const D5& saa, const D5& sab, const D5& sbb)
{
    const double a = 0.04918, b = 0.132, c = 0.2533, d = 0.349;
    const double CF = 0.3 * ::pow(3.0 * kPi * kPi, 2.0 / 3.0);
    const D5 rho = ra + rb;
    const D5 r = dpow(rho, -1.0 / 3.0);
    const D5 den = 1.0 + d * r;
    const D5 w = dexp(-c * r) / den * dpow(rho, -11.0 / 3.0);
    const D5 delta = c * r + d * r / den;
    const D5 s = saa + 2.0 * sab + sbb;
    const D5 rho2 = rho * rho;
    D5 t = ra * rb
        * (::pow(2.0, 11.0 / 3.0) * CF * (dpow(ra, 8.0 / 3.0) + dpow(rb, 8.0 / 3.0)) + (47.0 / 18.0 - 7.0 / 18.0 * delta) * s
           - (2.5 - delta / 18.0) * (saa + sbb) - (delta - 11.0) / 9.0 * (ra / rho * saa + rb / rho * sbb));
    t = t - (2.0 / 3.0) * rho2 * s + ((2.0 / 3.0) * rho2 - ra * ra) * sbb + ((2.0 / 3.0) * rho2 - rb * rb) * saa;
    return -4.0 * a / den * ra * rb / rho - a * b * w * t;
}

/// Correlation energy density of one kernel at (ra, rb, saa, sab, sbb) as a dual number; false: no spin form.
QM_HD bool correlationSpin(XCKernel k, const D5& ra, const D5& rb, const D5& saa, const D5& sab, const D5& sbb, D5& e)
{
    const D5 rho = ra + rb;
    switch (k) {
        case XCKernel::VWN5_C:
        case XCKernel::VWNRPA_C:
        case XCKernel::PW92_C:
        case XCKernel::PBE_C: {
            const D5 rs = dpow(3.0 / (4.0 * kPi * rho), 1.0 / 3.0);
            const D5 z = (ra - rb) / rho;
            D5 eps;
            if (k == XCKernel::PW92_C)
                eps = pw92Spin(rs, z);
            else if (k == XCKernel::PBE_C) {
                D5 sg = saa + 2.0 * sab + sbb;
                if (sg.v < 0.0) sg.v = 0.0;  // round-off of sigma_ab
                eps = pbeCSpin(rho, rs, z, sg);
            } else
                eps = vwnSpin(rs, z, k == XCKernel::VWNRPA_C);
            e = rho * eps;
            return true;
        }
        case XCKernel::LYP_C:
            e = lypSpinDensity(ra, rb, saa, sab, sbb);
            return true;
        default:
            return false;
    }
}

}  // namespace spin

/**
 * @brief Spin-polarised functional at one grid point (terms[0..nterms), weights applied); the body of
 * qmxc::evaluatePointSpin (qm_xc_spin.h, contract there) for host and device.
 *
 * ex, ec: exchange and correlation ENERGY DENSITIES (per volume); v[0..4] = d e / d (rho_a, rho_b, sigma_aa,
 * sigma_ab, sigma_bb). False (all outputs zero) below rho_min; a kernel without a spin form gives ex = ec = NaN and
 * false (never a silent zero). Exchange by spin scaling with the closed-shell kernels of qm_xc_kernels.h.
 */
QM_HD bool evaluatePointSpinTerms(const XCTerm* terms, int nterms, bool gga, double rho_min, double ra, double rb,
                                  double saa, double sab, double sbb, double& ex, double& ec, double v[5])
{
    using namespace spin;
    ex = ec = 0.0;
    for (int k = 0; k < 5; ++k) v[k] = 0.0;
    ra = ra > 0.0 ? ra : 0.0;  // negative round-off densities
    rb = rb > 0.0 ? rb : 0.0;
    const double rho = ra + rb;
    if (!(rho >= rho_min) || rho <= 0.0) return false;  // also NaN
    if (!gga) saa = sab = sbb = 0.0;
    saa = saa > 0.0 ? saa : 0.0;
    sbb = sbb > 0.0 ? sbb : 0.0;
    const D5 Ra = var(ra, 0), Rb = var(rb, 1), Saa = var(saa, 2), Sab = var(sab, 3), Sbb = var(sbb, 4);
    for (int it = 0; it < nterms; ++it) {
        const double w = terms[it].weight;
        if (w == 0.0) continue;
        if (kern::isExchange(terms[it].kernel)) {
            // spin scaling: e = sum_s rho_s eps(2 rho_s, 4 sigma_ss) with the closed-shell kernel
            const double rs[2] = { ra, rb }, ss[2] = { saa, sbb };
            for (int s = 0; s < 2; ++s) {
                const double r2 = 2.0 * rs[s];
                if (!(r2 >= rho_min) || r2 <= 0.0) continue;
                double eps = 0.0, vr = 0.0, vs = 0.0;
                kern::evaluateKernel(terms[it].kernel, r2, gga ? 4.0 * ss[s] : 0.0, eps, vr, vs);
                ex += w * rs[s] * eps;
                v[s] += w * vr;
                if (gga) v[s == 0 ? 2 : 4] += w * 2.0 * vs;
            }
            continue;
        }
        D5 e;
        if (!spin::correlationSpin(terms[it].kernel, Ra, Rb, Saa, Sab, Sbb, e)) {
            ex = ec = NAN;  // never a silent zero
            return false;
        }
        ec += w * e.v;
        v[0] += w * e.d[0];
        v[1] += w * e.d[1];
        if (gga) {
            v[2] += w * e.d[2];
            v[3] += w * e.d[3];
            v[4] += w * e.d[4];
        }
    }
    return true;
}

}  // namespace kern
}  // namespace qmxc

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif
