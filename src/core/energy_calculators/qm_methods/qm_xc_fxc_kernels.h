/*
 * <Native QM: second functional derivatives (XC kernel f_xc) of the double-hybrid exchange-correlation kernels,
 *  point functions shared by the CPU (qm_xc_response.cpp) and the GPU (cuda/qm_gpu_xc.cu)>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Why a second set of kernel formulas: the production kernels of qm_xc_kernels.h return e, v_rho, v_sigma, v_tau
 * with hand-written (GGA) or dual-number (meta-GGA) FIRST derivatives; they are shared with the CUDA build and their
 * bits are frozen by the regression checks. The double-hybrid gradient needs the SECOND derivatives
 * d^2 e / dx dy (x, y in rho, sigma, tau) for the coupled-perturbed Kohn-Sham (Z-vector) equations and for the
 * derivative of sum_mn P2_mn V_xc,mn. Hand-written second derivatives of LYP, P86 and B95 are long and error-prone,
 * so the energy densities e = rho eps are written here once more, as plain energy expressions, and differentiated by
 * forward-mode automatic differentiation to second order (Jet: value, 3 first and 6 second partials, exact to
 * rounding; no finite differences). The jet first derivatives are checked against the production v's and the second
 * derivatives against libxc (PySCF eval_xc deriv = 2) and against finite differences of the production potentials
 * (ctest qm_dh_grad_kernel), so a transcription error in this copy is caught.
 *
 * Kernels covered (the ones of the double-hybrid table): B88_X, LYP_C, PBE_X, P86_C (PZ81 local), P86VWN_C (VWN5
 * local), B95_C (meta-GGA). evaluatePointFxc returns false for any other kernel (the caller then refuses the
 * gradient instead of using a wrong kernel). Variables are the closed-shell total quantities of the production code
 * (rho, sigma = |grad rho|^2, tau = 1/2 sum_i |grad psi_i|^2), with the same clamps (sigma >= 0; B95: tau floor and
 * sigma <= 8 rho tau, derivatives at the clamped point as production). P86 uses sqrt(sigma): its second sigma
 * derivative diverges like sigma^(-1/2) at sigma -> 0 (a property of the functional, contracted with
 * sigma_1 = 2 grad rho . grad rho_1 it stays finite); sigma is raised to 1e-40 inside the P86 jet only.
 *
 * Like qm_xc_kernels.h: only plain C math in the global namespace, no Eigen, no std containers, so the same text
 * compiles as host code (g++) and as __host__ __device__ code (nvcc). The host arithmetic is that of the former copy
 * in qm_xc_fxc.h (std::exp etc. there are the same libm functions), so the CPU results did not change by the move
 * (checked with -dump_gradient md5 against d0f21497); device math may differ from glibc by an ulp.
 *
 * Claude Generated (Oct 2026, WP12 stage 2; moved here for the device in WP12 stage 3, Oct 31 2026,
 * docs/DFT_ROADMAP/WP12_doppelhybride.md). AI-generated, machine-tested; human production test pending.
 *
 * This program is free software under GPL-3.0
 */

#pragma once

#include "qm_xc_kernels.h"

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Winline"
#endif

namespace qmxc {
namespace fxc {

/// Value, gradient (3) and Hessian (6: 00 01 02 11 12 22) of a function of (rho, sigma, tau).
struct Jet {
    double v = 0.0;
    double g[3] = { 0.0, 0.0, 0.0 };
    double h[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
};

QM_HD Jet jconst(double c)
{
    Jet r;
    r.v = c;
    return r;
}
QM_HD Jet jvar(double x, int k)
{
    Jet r;
    r.v = x;
    r.g[k] = 1.0;
    return r;
}
/// f(a) with f = fv, f' = f1, f'' = f2: h_ij = f1 a.h_ij + f2 a.g_i a.g_j.
QM_HD Jet jchain(const Jet& a, double fv, double f1, double f2)
{
    Jet r;
    r.v = fv;
    for (int k = 0; k < 3; ++k) r.g[k] = f1 * a.g[k];
    const int I[6] = { 0, 0, 0, 1, 1, 2 }, J[6] = { 0, 1, 2, 1, 2, 2 };
    for (int k = 0; k < 6; ++k) r.h[k] = f1 * a.h[k] + f2 * a.g[I[k]] * a.g[J[k]];
    return r;
}
QM_HD Jet operator+(const Jet& a, const Jet& b)
{
    Jet r;
    r.v = a.v + b.v;
    for (int k = 0; k < 3; ++k) r.g[k] = a.g[k] + b.g[k];
    for (int k = 0; k < 6; ++k) r.h[k] = a.h[k] + b.h[k];
    return r;
}
QM_HD Jet operator-(const Jet& a, const Jet& b)
{
    Jet r;
    r.v = a.v - b.v;
    for (int k = 0; k < 3; ++k) r.g[k] = a.g[k] - b.g[k];
    for (int k = 0; k < 6; ++k) r.h[k] = a.h[k] - b.h[k];
    return r;
}
QM_HD Jet operator-(const Jet& a) { return jchain(a, -a.v, -1.0, 0.0); }
QM_HD Jet operator*(const Jet& a, const Jet& b)
{
    Jet r;
    r.v = a.v * b.v;
    for (int k = 0; k < 3; ++k) r.g[k] = a.g[k] * b.v + a.v * b.g[k];
    const int I[6] = { 0, 0, 0, 1, 1, 2 }, J[6] = { 0, 1, 2, 1, 2, 2 };
    for (int k = 0; k < 6; ++k)
        r.h[k] = a.h[k] * b.v + a.v * b.h[k] + a.g[I[k]] * b.g[J[k]] + a.g[J[k]] * b.g[I[k]];
    return r;
}
QM_HD Jet jinv(const Jet& a) { return jchain(a, 1.0 / a.v, -1.0 / (a.v * a.v), 2.0 / (a.v * a.v * a.v)); }
QM_HD Jet operator/(const Jet& a, const Jet& b) { return a * jinv(b); }
QM_HD Jet operator+(const Jet& a, double c) { return jchain(a, a.v + c, 1.0, 0.0); }
QM_HD Jet operator+(double c, const Jet& a) { return jchain(a, c + a.v, 1.0, 0.0); }
QM_HD Jet operator-(const Jet& a, double c) { return jchain(a, a.v - c, 1.0, 0.0); }
QM_HD Jet operator-(double c, const Jet& a) { return jchain(a, c - a.v, -1.0, 0.0); }
QM_HD Jet operator*(const Jet& a, double c) { return jchain(a, a.v * c, c, 0.0); }
QM_HD Jet operator*(double c, const Jet& a) { return jchain(a, c * a.v, c, 0.0); }
QM_HD Jet operator/(const Jet& a, double c) { return jchain(a, a.v / c, 1.0 / c, 0.0); }
QM_HD Jet operator/(double c, const Jet& a) { return jchain(a, c / a.v, -c / (a.v * a.v), 2.0 * c / (a.v * a.v * a.v)); }
QM_HD Jet jexp(const Jet& a)
{
    const double e = ::exp(a.v);
    return jchain(a, e, e, e);
}
QM_HD Jet jlog(const Jet& a) { return jchain(a, ::log(a.v), 1.0 / a.v, -1.0 / (a.v * a.v)); }
QM_HD Jet jsqrt(const Jet& a)
{
    const double s = ::sqrt(a.v);
    return jchain(a, s, 0.5 / s, -0.25 / (s * a.v));
}
QM_HD Jet jcbrt(const Jet& a)
{
    const double c = ::cbrt(a.v);
    return jchain(a, c, c / (3.0 * a.v), -2.0 * c / (9.0 * a.v * a.v));
}
QM_HD Jet jatan(const Jet& a)
{
    const double d = 1.0 + a.v * a.v;
    return jchain(a, ::atan(a.v), 1.0 / d, -2.0 * a.v / (d * d));
}
/// u(y) = sqrt(y) asinh(sqrt(y)) of B88 as one elementary function of y >= 0 (series below 1e-3: the closed-form
/// u'' cancels for small y).
QM_HD Jet jB88u(const Jet& y)
{
    const double Y = y.v;
    double u, u1, u2;
    if (Y < 1e-3) {
        // asinh x = x - x^3/6 + 3x^5/40 - 15x^7/336 + 105x^9/3456 - ...,  u = x asinh x
        u = Y * (1.0 + Y * (-1.0 / 6.0 + Y * (3.0 / 40.0 + Y * (-15.0 / 336.0 + Y * (105.0 / 3456.0)))));
        u1 = 1.0 + Y * (-2.0 / 6.0 + Y * (9.0 / 40.0 + Y * (-60.0 / 336.0 + Y * (525.0 / 3456.0))));
        u2 = -2.0 / 6.0 + Y * (18.0 / 40.0 + Y * (-180.0 / 336.0 + Y * (2100.0 / 3456.0)));
    } else {
        const double x = ::sqrt(Y), A = ::asinh(x), s = ::sqrt(1.0 + Y), B = x / s;
        u = x * A;
        u1 = (A + B) / (2.0 * x);
        // d/dy (A + B) = 1/(2 x s) + 1/(2 x s) - x/(2 s^3) ;  u'' = [ (A+B)' 2x - (A+B)/x ] / (4 x^2)
        const double dAB = 1.0 / (x * s) - x / (2.0 * s * s * s);
        u2 = (dAB * 2.0 * x - (A + B) / x) / (4.0 * Y);
    }
    return jchain(y, u, u1, u2);
}

// ------------------------------- energy densities e = rho eps (closed shell) -------------------------------------
// Same formulas, constants and equation references as qm_xc_kernels.h (comments there).

/// PW92 G(r_s) for one parameter set (PW92 eq. 10).
QM_HD Jet pw92G(const Jet& rs, double A, double a1, double b1, double b2, double b3, double b4)
{
    const Jet srs = jsqrt(rs);
    const Jet Q1 = 2.0 * A * (b1 * srs + b2 * rs + b3 * rs * srs + b4 * rs * rs);
    return -2.0 * A * (1.0 + a1 * rs) * jlog(1.0 + 1.0 / Q1);
}

QM_HD Jet b88(const Jet& R, const Jet& S)
{
    using namespace qmxc::kern;
    const double Cx = 1.5 * ::cbrt(3.0 / (4.0 * kPi));
    const Jet R13 = jcbrt(R);
    const Jet y = ::cbrt(4.0) * S / (R * R * R13 * R13);
    const Jet g = y / (1.0 + 6.0 * kB88_beta * jB88u(y));
    return R * (R13 / ::cbrt(2.0)) * (-Cx - kB88_beta * g);
}

QM_HD Jet lyp(const Jet& R, const Jet& S)
{
    using namespace qmxc::kern;
    const double CF = 0.3 * ::pow(3.0 * kPi * kPi, 2.0 / 3.0);
    const Jet r = 1.0 / jcbrt(R);
    const Jet den = 1.0 + kLYP_d * r;
    const Jet E = jexp(-kLYP_c * r);
    const Jet F = -kLYP_a / den;
    const Jet G = -kLYP_a * kLYP_b * E / den;
    const Jet delta = kLYP_c * r + kLYP_d * r / den;
    const Jet K = 1.0 / 24.0 + 7.0 * delta / 72.0;
    const Jet r2 = r * r;
    const Jet H = G * r2 * r2 * r * K;
    return R * (F + G * CF) - H * S;
}

QM_HD Jet pbeX(const Jet& R, const Jet& S)
{
    using namespace qmxc::kern;
    const Jet R13 = jcbrt(R);
    const Jet eps_lda = -0.75 * ::cbrt(3.0 / kPi) * R13;
    const Jet c = 1.0 / (4.0 * ::pow(3.0 * kPi * kPi, 2.0 / 3.0) * R * R * R13 * R13);
    const Jet s2 = c * S;
    const Jet F = 1.0 + kPBE_kappa - kPBE_kappa / (1.0 + kPBE_mu * s2 / kPBE_kappa);
    return R * eps_lda * F;
}

QM_HD Jet pz81(const Jet& R)
{
    using namespace qmxc::kern;
    const Jet rs = jcbrt(3.0 / (4.0 * kPi) / R);
    Jet eps;
    if (rs.v >= 1.0)
        eps = kPZ_gamma / (1.0 + kPZ_beta1 * jsqrt(rs) + kPZ_beta2 * rs);
    else {
        const Jet lr = jlog(rs);
        eps = kPZ_A * lr + kPZ_B + kPZ_C * rs * lr + kPZ_D * rs;
    }
    return R * eps;
}

QM_HD Jet vwn5(const Jet& R)
{
    const double A = 0.0310907, x0 = -0.10498, b = 3.72744, c = 12.9352;
    const Jet rs = jcbrt(3.0 / (4.0 * qmxc::kern::kPi) / R);
    const Jet x = jsqrt(rs);
    const Jet X = x * x + b * x + c;
    const double X0 = x0 * x0 + b * x0 + c;
    const double Q = ::sqrt(4.0 * c - b * b);
    const Jet at = jatan(Q / (2.0 * x + b));
    const double bx0 = b * x0 / X0;
    const Jet xm = x - x0;
    const Jet eps = A * (jlog(x * x / X) + 2.0 * b / Q * at - bx0 * (jlog(xm * xm / X) + 2.0 * (b + 2.0 * x0) / Q * at));
    return R * eps;
}

QM_HD Jet p86Gradient(const Jet& R, const Jet& Sin)
{
    using namespace qmxc::kern;
    const Jet S = Sin.v > 1e-40 ? Sin : jvar(1e-40, 1);  // sqrt(sigma) below (file comment)
    const Jet rs = jcbrt(3.0 / (4.0 * kPi) / R);
    const Jet N = kP86_b + kP86_alpha * rs + kP86_beta * rs * rs;
    const Jet Dn = 1.0 + kP86_gamma * rs + kP86_delta * rs * rs + 1.0e4 * kP86_beta * rs * rs * rs;
    const Jet C = kP86_a + N / Dn;
    const Jet R13 = jcbrt(R);
    const Jet rm43 = 1.0 / (R * R13);
    const Jet y = jsqrt(S) / (R * jsqrt(R13));
    const Jet Phi = kP86_ftilde * (kP86_a + kP86_b) / C * y;
    return S * rm43 * C * jexp(-Phi);
}

QM_HD Jet b95(const Jet& R, const Jet& S, const Jet& T)
{
    using namespace qmxc::kern;
    const Jet n = 0.5 * R;
    const Jet rs_n = jcbrt(3.0 / (4.0 * kPi) / n);
    const Jet rs_R = jcbrt(3.0 / (4.0 * kPi) / R);
    const Jet eF = pw92G(rs_n, kPWF_A, kPWF_alpha1, kPWF_beta1, kPWF_beta2, kPWF_beta3, kPWF_beta4);
    const Jet eP = pw92G(rs_R, kPW_A, kPW_alpha1, kPW_beta1, kPW_beta2, kPW_beta3, kPW_beta4);
    const Jet n13 = jcbrt(n);
    const Jet n53 = n * n13 * n13;
    const Jet xs2 = 0.25 * S / (n53 * n);
    const Jet ts = 0.5 * T / n53;
    const double KF = 0.3 * ::cbrt(36.0 * kPi * kPi * kPi * kPi);
    const Jet q = 1.0 + kB95_css * xs2;
    const Jet gpar = (ts - xs2 / 8.0) / (KF * q * q);
    const Jet gperp = 1.0 / (1.0 + 2.0 * kB95_copp * xs2);
    return R * (eF * gpar + (eP - eF) * gperp);
}

/// Energy density jet of one kernel at (rho, sigma, tau) after the production clamps; false if not covered.
QM_HD bool kernelJet(XCKernel k, double rho, double sigma, double tau, Jet& e)
{
    const double s = sigma > 0.0 ? sigma : 0.0;
    switch (k) {
        case XCKernel::B88_X: e = b88(jvar(rho, 0), jvar(s, 1)); return true;
        case XCKernel::LYP_C: e = lyp(jvar(rho, 0), jvar(s, 1)); return true;
        case XCKernel::PBE_X: e = pbeX(jvar(rho, 0), jvar(s, 1)); return true;
        case XCKernel::P86_C: e = pz81(jvar(rho, 0)) + p86Gradient(jvar(rho, 0), jvar(s, 1)); return true;
        case XCKernel::P86VWN_C: e = vwn5(jvar(rho, 0)) + p86Gradient(jvar(rho, 0), jvar(s, 1)); return true;
        case XCKernel::B95_C: {
            const double t = tau > kern::kMGGATauFloor ? tau : kern::kMGGATauFloor;
            const double sc = s > 8.0 * rho * t ? 8.0 * rho * t : s;
            e = b95(jvar(rho, 0), jvar(sc, 1), jvar(t, 2));
            return true;
        }
        default: return false;
    }
}

/**
 * @brief First and second derivatives of the weighted functional at one point.
 *
 * v[3] = (v_rho, v_sigma, v_tau), H[6] = second derivatives (rr, rs, rt, ss, st, tt), weights applied, terms in list
 * order. Points with rho < rho_min, rho <= 0 or NaN give zeros and false (they are left out, as in V_xc).
 */
QM_HD bool evaluatePointFxc(const XCTerm* terms, int nterms, double rho_min, double rho, double sigma, double tau,
                            double v[3], double H[6])
{
    for (int k = 0; k < 3; ++k) v[k] = 0.0;
    for (int k = 0; k < 6; ++k) H[k] = 0.0;
    if (!(rho >= rho_min) || rho <= 0.0) return false;
    for (int i = 0; i < nterms; ++i) {
        const double w = terms[i].weight;
        if (w == 0.0) continue;
        Jet e;
        if (!kernelJet(terms[i].kernel, rho, sigma, tau, e)) continue;
        for (int k = 0; k < 3; ++k) v[k] += w * e.g[k];
        for (int k = 0; k < 6; ++k) H[k] += w * e.h[k];
    }
    return true;
}

}  // namespace fxc
}  // namespace qmxc

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif
