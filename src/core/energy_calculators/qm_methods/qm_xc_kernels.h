/*
 * <Native QM: point kernels of the DFT exchange-correlation functionals, the Becke partition polynomials and
 *  one AO component on a grid point, shared by the CPU (qm_xc.cpp, qm_grid.cpp) and the GPU (cuda/qm_gpu_xc.cu)>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * ONE source of the formulas for both sides: the CPU wrappers in qm_xc.cpp / qm_grid.cpp and the device kernels of
 * cuda/qm_gpu_xc.cu call these inline functions, so CPU and GPU cannot drift apart. Formulas, derivations and
 * literature: qm_xc.h (XC kernels) and qm_grid.h (partition, AO values). A change here changes both sides; the
 * ctests qm_lda/qm_pbe/qm_b3lyp (CPU vs libxc) and qm_gpu_dft (GPU vs CPU) cover it.
 *
 * Only plain C math (cbrt, sqrt, log, exp, atan, asinh, pow from <cmath>, global namespace), no Eigen, no
 * std containers: compiles as host code (g++) and as __host__ __device__ code (nvcc; hipcc for a later ROCm port,
 * not built or tested). Device math (CUDA libdevice) may differ from glibc by an ulp; the GPU results agree with the
 * CPU to rounding, not bitwise.
 *
 * Claude Generated (Oct 2026, GPU KS-DFT G1: functions moved here from qm_xc.cpp / qm_grid.cpp with the same
 * arithmetic, so the CPU results are unchanged). AI-generated, machine-tested only.
 *
 * This program is free software under GPL-3.0
 */

#pragma once

#include <cmath>

#if defined(__CUDACC__) || defined(__HIPCC__)
#define QM_HD __host__ __device__ inline
#else
#define QM_HD inline
#endif

// The project builds with -Winline; GCC may decline to inline the larger kernels into evaluateKernel (only a
// diagnostic, the code is the same as before the move). Silenced for this header only.
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Winline"
#endif

namespace qmxc {

/**
 * @brief One closed-shell exchange or correlation kernel (a point function of rho and, for a GGA, sigma).
 *
 * A functional is a weighted sum of kernels (XCFunctional::terms). WP7 (B3LYP) combines B88 exchange and LYP
 * correlation with SlaterX and VWN5_C or VWNRPA_C.
 */
enum class XCKernel {
    SlaterX,   ///< LDA Dirac/Slater exchange (WP5)
    VWN5_C,    ///< LDA VWN5 correlation (WP5)
    VWNRPA_C,  ///< LDA VWN-RPA correlation (WP5)
    PW92_C,    ///< LDA Perdew-Wang 1992 correlation with libxc's LDA_C_PW_MOD constants (WP6, part of PBE)
    PBE_X,     ///< GGA PBE exchange (WP6), libxc GGA_X_PBE
    PBE_C,     ///< GGA PBE correlation = PW92 + H(rs, t) (WP6), libxc GGA_C_PBE
    B88_X,     ///< GGA Becke 1988 exchange (WP7), libxc GGA_X_B88
    LYP_C      ///< GGA Lee-Yang-Parr correlation, Miehlich form without Laplacian (WP7), libxc GGA_C_LYP
};

/// Weighted kernel: E_xc contribution = weight * integral rho eps_kernel.
struct XCTerm {
    XCKernel kernel = XCKernel::SlaterX;
    double weight = 1.0;
};

namespace kern {

constexpr double kPi = 3.14159265358979323846;

/// Paramagnetic VWN parameters {A, x0, b, c} (Hartree), libxc lda_c_vwn.c / VWN 1980 table.
struct VWNParams {
    double A, x0, b, c;
};
// vwn5: { 0.0310907, -0.10498, 3.72744, 12.9352 }, vwn_rpa: { 0.0310907, -0.409286, 13.0720, 42.7198 }; written out
// in vwnCorrelation (a namespace-scope constexpr struct is not readable in nvcc device code)

// PW92 paramagnetic correlation, libxc LDA_C_PW_MOD (lda_c_pw.c): A with the extra digits of PW92's own Fortran
// code (the paper prints 0.031091), p = 1 (the r_s^(p+1) term is beta4 r_s^2).
constexpr double kPW_A = 0.0310907, kPW_alpha1 = 0.21370;
constexpr double kPW_beta1 = 7.5957, kPW_beta2 = 3.5876, kPW_beta3 = 1.6382, kPW_beta4 = 0.49294;

// PBE (libxc gga_x_pbe.c / gga_c_pbe.c). mu = beta pi^2 / 3 with the PBE beta (PBE 1996 below eq. 12);
// gamma = (1 - ln 2) / pi^2 is formed where it is used (PBE 1996 below eq. 7).
constexpr double kPBE_kappa = 0.8040;
constexpr double kPBE_mu = 0.2195149727645171;
constexpr double kPBE_beta = 0.06672455060314922;

// B88 (Becke 1988, libxc gga_x_b88.c): beta, gamma = 6; C_x = (3/2) (3/(4 pi))^(1/3) is the per-spin Slater constant.
constexpr double kB88_beta = 0.0042;
// LYP (Lee-Yang-Parr 1988, libxc gga_c_lyp.c): a, b, c, d.
constexpr double kLYP_a = 0.04918, kLYP_b = 0.132, kLYP_c = 0.2533, kLYP_d = 0.349;

// Densities below this are never handed to a GGA kernel, whatever -qm.xc_density_threshold says (also 0): near
// rho ~ 1e-27 ln(1 + 1/Q1) rounds to 0, A = inf and the PBE formulas give NaN; rho^(8/3) underflows below ~1e-116.
// Far below any physical threshold (libxc: 1e-15 exchange, 1e-12 correlation). LDA kernels have no such floor.
constexpr double kGGADensityFloor = 1e-20;

/// True for an exchange kernel (its energy is reported as E_x), false for correlation.
QM_HD bool isExchange(XCKernel k) { return k == XCKernel::SlaterX || k == XCKernel::PBE_X || k == XCKernel::B88_X; }
/// True when the kernel depends on sigma = |grad rho|^2 (GGA).
QM_HD bool isGGA(XCKernel k)
{
    return k == XCKernel::PBE_X || k == XCKernel::PBE_C || k == XCKernel::B88_X || k == XCKernel::LYP_C;
}
/// Density below which a point is left out of E_xc and V_xc (gga: functional needs sigma).
QM_HD double effectiveThreshold(bool gga, double rho_threshold)
{
    return gga ? (rho_threshold > kGGADensityFloor ? rho_threshold : kGGADensityFloor) : rho_threshold;
}

/// PW92 G(r_s) of the paramagnetic gas and d G / d r_s (PW92 eq. 10 and its derivative, eq. A3 there).
///   G = -2A (1 + alpha1 r_s) ln(1 + 1/Q1),  Q1 = 2A (beta1 r_s^(1/2) + beta2 r_s + beta3 r_s^(3/2) + beta4 r_s^2)
///   dG/dr_s = -2A alpha1 ln(1 + 1/Q1) + 2A (1 + alpha1 r_s) Q1' / (Q1^2 + Q1)
/// ln(1 + x) rather than log1p(x), as libxc's generated code: at r_s ~ 6e3 (rho 1e-12) log1p would differ from libxc
/// by 1e-10 relative (libxc's rounding), everywhere else both agree to 1e-15.
QM_HD void pw92G(double rs, double& G, double& dG)
{
    const double srs = ::sqrt(rs);
    const double Q1 = 2.0 * kPW_A * (kPW_beta1 * srs + kPW_beta2 * rs + kPW_beta3 * rs * srs + kPW_beta4 * rs * rs);
    const double Q1p = kPW_A * (kPW_beta1 / srs + 2.0 * kPW_beta2 + 3.0 * kPW_beta3 * srs + 4.0 * kPW_beta4 * rs);
    const double L = ::log(1.0 + 1.0 / Q1);
    G = -2.0 * kPW_A * (1.0 + kPW_alpha1 * rs) * L;
    dG = -2.0 * kPW_A * kPW_alpha1 * L + 2.0 * kPW_A * (1.0 + kPW_alpha1 * rs) * Q1p / (Q1 * Q1 + Q1);
}

// Dirac (1930) / Slater (1951): e_x(rho) = rho eps_x = -(3/4) (3/pi)^(1/3) rho^(4/3), closed shell.
QM_HD void slaterExchange(double rho, double& eps, double& v)
{
    const double t = ::cbrt(3.0 * rho / kPi);  // (3 rho / pi)^(1/3)
    eps = -0.75 * t;
    v = -t;  // (4/3) eps_x
}

// Vosko-Wilk-Nusair (1980) eq. 4.4 with the paramagnetic parameters; derivative with respect to x = sqrt(r_s):
//   d eps/dx = A [ 2/x - (2x+b)/X - b/X - b x0/X0 ( 2/(x-x0) - (2x+b)/X - (b+2 x0)/X ) ]
// (the atan term gives -4b/((2x+b)^2 + Q^2) = -b/X since (2x+b)^2 + Q^2 = 4X), and
//   v = d(rho eps)/d rho = eps - (x/6) d eps/dx      (dx/d rho = -x/(6 rho)).
QM_HD void vwnCorrelation(double rho, bool rpa, double& eps, double& v)
{
    const VWNParams p = rpa ? VWNParams{ 0.0310907, -0.409286, 13.0720, 42.7198 }
                             : VWNParams{ 0.0310907, -0.10498, 3.72744, 12.9352 };
    const double rs = ::cbrt(3.0 / (4.0 * kPi * rho));
    const double x = ::sqrt(rs);
    const double X = x * x + p.b * x + p.c;
    const double X0 = p.x0 * p.x0 + p.b * p.x0 + p.c;
    const double Q = ::sqrt(4.0 * p.c - p.b * p.b);
    const double at = ::atan(Q / (2.0 * x + p.b));
    const double bx0 = p.b * p.x0 / X0;
    eps = p.A * (::log(x * x / X) + 2.0 * p.b / Q * at
                 - bx0 * (::log((x - p.x0) * (x - p.x0) / X) + 2.0 * (p.b + 2.0 * p.x0) / Q * at));
    const double dX = (2.0 * x + p.b) / X;
    const double deps_dx = p.A * (2.0 / x - dX - p.b / X - bx0 * (2.0 / (x - p.x0) - dX - (p.b + 2.0 * p.x0) / X));
    v = eps - x / 6.0 * deps_dx;
}

// PW92 (1992) paramagnetic correlation as an LDA kernel: v = eps - (r_s/3) d eps/d r_s  (d r_s/d rho = -r_s/(3 rho)).
QM_HD void pw92Correlation(double rho, double& eps, double& v)
{
    const double rs = ::cbrt(3.0 / (4.0 * kPi * rho));
    double dG = 0.0;
    pw92G(rs, eps, dG);
    v = eps - rs / 3.0 * dG;
}

// PBE exchange (PBE 1996 eqs. 10 and 14), closed shell, e_x = rho eps_x^LDA(rho) F_x(s^2):
//   s^2 = c sigma,  c = 1 / (4 (3 pi^2)^(2/3) rho^(8/3)),  F_x' = dF_x/d(s^2) = mu / (1 + mu s^2/kappa)^2
//   v_rho   = d e_x/d rho   = eps_x^LDA [ (4/3) F_x - (8/3) s^2 F_x' ]      (d s^2/d rho = -(8/3) s^2/rho)
//   v_sigma = d e_x/d sigma = rho eps_x^LDA F_x' c
QM_HD void pbeExchange(double rho, double sigma, double& eps, double& vrho, double& vsigma)
{
    const double eps_lda = -0.75 * ::cbrt(3.0 * rho / kPi);
    const double rho13 = ::cbrt(rho);
    const double c = 1.0 / (4.0 * ::pow(3.0 * kPi * kPi, 2.0 / 3.0) * rho * rho * rho13 * rho13);
    const double s2 = c * sigma;
    const double d = 1.0 + kPBE_mu * s2 / kPBE_kappa;
    const double F = 1.0 + kPBE_kappa - kPBE_kappa / d;   // eq. 14
    const double Fp = kPBE_mu / (d * d);
    eps = eps_lda * F;
    vrho = eps_lda * (4.0 / 3.0 * F - 8.0 / 3.0 * s2 * Fp);
    vsigma = rho * eps_lda * Fp * c;
}

// PBE correlation (PBE 1996 eqs. 3, 7, 8), zeta = 0 (phi = 1), e_c = rho (eps_c^PW92 + H(eps_c^PW92, u)), u = t^2:
//   u = sigma pi / (16 k_F rho^2)                        (t = |grad rho| / (2 k_s rho), k_s^2 = 4 k_F / pi)
//   A = (beta/gamma) / (E - 1),  E = exp(-eps_c^PW92 / gamma)                                         (eq. 8)
//   y = u (1 + z) / D,  z = A u,  D = 1 + z + z^2,  H = gamma ln(1 + (beta/gamma) y)                  (eq. 7)
// Partial derivatives (by hand, checked by finite differences in ctest qm_pbe_kernel):
//   dH/dy = beta / (1 + (beta/gamma) y),  dy/du = (1 + 2z) / D^2,  dy/dA = -u^2 z (2 + z) / D^2,
//   dA/d eps_c = A^2 E / beta
//   v_rho   = eps_c + H - (r_s/3) (1 + dH/d eps_c) d eps_c/d r_s - (7/3) u dH/du       (du/d rho = -(7/3) u / rho)
//   v_sigma = rho dH/du u / sigma = dH/du pi / (16 k_F rho)
// Large t (u -> inf): y -> 1/A and H -> -eps_c^PW92, so eps_c -> 0 by cancellation; the result there is only
// accurate to ~1e-16 * |eps_c^PW92| (absolute), as libxc's. The factors are grouped so that nothing overflows for
// z up to ~1e150.
QM_HD void pbeCorrelation(double rho, double sigma, double& eps, double& vrho, double& vsigma)
{
    const double gamma = (1.0 - ::log(2.0)) / (kPi * kPi);  // (1 - ln 2) / pi^2, PBE 1996 below eq. 7
    const double rs = ::cbrt(3.0 / (4.0 * kPi * rho));
    double ec = 0.0, dec = 0.0;
    pw92G(rs, ec, dec);
    const double kF = ::cbrt(3.0 * kPi * kPi * rho);
    const double u = sigma * kPi / (16.0 * kF * rho * rho);
    const double E = ::exp(-ec / gamma);
    const double A = (kPBE_beta / gamma) / (E - 1.0);
    const double z = A * u;
    const double D = 1.0 + z + z * z;
    const double y = u * ((1.0 + z) / D);
    const double arg = 1.0 + kPBE_beta / gamma * y;
    const double H = gamma * ::log(arg);
    const double H_y = kPBE_beta / arg;
    const double y_u = (1.0 + 2.0 * z) / D / D;
    const double y_A = -(u * z / D) * (u * (2.0 + z) / D);
    const double A_ec = A * A * E / kPBE_beta;
    const double H_ec = H_y * y_A * A_ec;
    const double H_u = H_y * y_u;
    eps = ec + H;
    vrho = ec + H - rs / 3.0 * (1.0 + H_ec) * dec - 7.0 / 3.0 * u * H_u;
    vsigma = H_u * kPi / (16.0 * kF * rho);
}

// Becke 1988 exchange (PRA 38, 3098, eq. 8), closed shell. With y = x_s^2 = 2^(2/3) sigma / rho^(8/3), x = sqrt(y),
//   g(y) = y / D,  D = 1 + 6 beta x asinh x,
//   e_x  = 2 (rho/2)^(4/3) [ -C_x - beta g ] = 2^(-1/3) rho^(4/3) [ -C_x - beta g ]
//   g'   = dg/dy = (D - 3 beta x (asinh x + x / sqrt(1 + y))) / D^2       (finite at y = 0: g' = 1)
//   v_rho   = 2^(-1/3) rho^(1/3) [ (4/3) (-C_x - beta g) + (8/3) beta y g' ]   (dy/d rho = -(8/3) y / rho)
//   v_sigma = -2^(1/3) beta g' rho^(-4/3)                                    (dy/d sigma = y / sigma)
// The derivative is taken in y, not in x = sqrt(y), so sigma = 0 needs no special case. Large x: x / sqrt(1 + y)
// is written as 1 / sqrt(1 + 1/y) (no overflow of y for any sigma, rho above the GGA floor).
QM_HD void b88Exchange(double rho, double sigma, double& eps, double& vrho, double& vsigma)
{
    const double Cx = 1.5 * ::cbrt(3.0 / (4.0 * kPi));
    const double rho13 = ::cbrt(rho);
    const double y = ::cbrt(4.0) * sigma / (rho * rho * rho13 * rho13);  // 2^(2/3) sigma / rho^(8/3)
    const double x = ::sqrt(y);
    const double ash = ::asinh(x);
    const double x_over = y > 1.0 ? 1.0 / ::sqrt(1.0 + 1.0 / y) : x / ::sqrt(1.0 + y);  // x / sqrt(1 + x^2)
    const double D = 1.0 + 6.0 * kB88_beta * x * ash;
    const double g = y / D;
    const double gp = (D - 3.0 * kB88_beta * x * (ash + x_over)) / (D * D);
    const double pre = rho13 / ::cbrt(2.0);  // 2^(-1/3) rho^(1/3)
    eps = pre * (-Cx - kB88_beta * g);
    vrho = pre * (4.0 / 3.0 * (-Cx - kB88_beta * g) + 8.0 / 3.0 * kB88_beta * y * gp);
    vsigma = -::cbrt(2.0) * kB88_beta * gp / (rho * rho13);
}

// LYP correlation, closed shell, in the Laplacian-free form of Miehlich, Savin, Stoll, Preuss (CPL 157, 200 (1989),
// eq. 2) with rho_a = rho_b = rho/2 and |grad rho_a|^2 = |grad rho_b|^2 = sigma/4 (all sums by hand):
//   the bracket of eq. 2 becomes 4 C_F rho^(8/3) + sigma (5/3 - 7 delta/18), times rho_a rho_b = rho^2/4, and the
//   last three terms give -(11/24) rho^2 sigma, so with r = rho^(-1/3), omega = exp(-c r)/(1 + d r) r^11:
//   e_c = F rho + G C_F rho - H sigma,
//   F = -a / (1 + d r),  G = -a b exp(-c r) / (1 + d r),  H = G r^5 K,  K = 1/24 + 7 delta / 72,
//   delta = c r + d r / (1 + d r).
// Derivatives in r (dr/d rho = -r / (3 rho)):
//   dF/dr = a d / (1 + d r)^2,  dG/dr = -G delta / r  (= -G (c + d / (1 + d r))),  d delta/dr = c + d / (1 + d r)^2,
//   dH/dr = r^5 K dG/dr + 5 r^4 K G + r^5 G (7/72) d delta/dr,
//   v_rho   = F + C_F G - (r/3) (dF/dr + C_F dG/dr) + sigma r / (3 rho) dH/dr,
//   v_sigma = -H.
// Small rho: exp(-c r) underflows to 0 below rho ~ 4e-11, then only the F term remains (no inf * 0: r^5 stays finite
// above the GGA floor 1e-20).
QM_HD void lypCorrelation(double rho, double sigma, double& eps, double& vrho, double& vsigma)
{
    const double CF = 0.3 * ::pow(3.0 * kPi * kPi, 2.0 / 3.0);
    const double r = 1.0 / ::cbrt(rho);
    const double den = 1.0 + kLYP_d * r;
    const double E = ::exp(-kLYP_c * r);
    const double F = -kLYP_a / den;
    const double G = -kLYP_a * kLYP_b * E / den;
    const double delta = kLYP_c * r + kLYP_d * r / den;
    const double K = 1.0 / 24.0 + 7.0 * delta / 72.0;
    const double r4 = r * r * r * r;
    const double H = G * r4 * r * K;
    const double dF = kLYP_a * kLYP_d / (den * den);
    const double dG = -G * (kLYP_c + kLYP_d / den);
    const double ddelta = kLYP_c + kLYP_d / (den * den);
    const double dH = dG * r4 * r * K + 5.0 * G * r4 * K + G * r4 * r * (7.0 / 72.0) * ddelta;
    eps = F + G * CF - H * sigma / rho;
    vrho = F + G * CF - r / 3.0 * (dF + CF * dG) + sigma * r / (3.0 * rho) * dH;
    vsigma = -H;
}

/// Any kernel by enum (sigma ignored by LDA kernels, vsigma = 0 for them). False for an unknown kernel.
QM_HD bool evaluateKernel(XCKernel k, double rho, double sigma, double& eps, double& vrho, double& vsigma)
{
    vsigma = 0.0;
    switch (k) {
        case XCKernel::SlaterX:  slaterExchange(rho, eps, vrho); return true;
        case XCKernel::VWN5_C:   vwnCorrelation(rho, false, eps, vrho); return true;
        case XCKernel::VWNRPA_C: vwnCorrelation(rho, true, eps, vrho); return true;
        case XCKernel::PW92_C:   pw92Correlation(rho, eps, vrho); return true;
        case XCKernel::PBE_X:    pbeExchange(rho, sigma, eps, vrho, vsigma); return true;
        case XCKernel::PBE_C:    pbeCorrelation(rho, sigma, eps, vrho, vsigma); return true;
        case XCKernel::B88_X:    b88Exchange(rho, sigma, eps, vrho, vsigma); return true;
        case XCKernel::LYP_C:    lypCorrelation(rho, sigma, eps, vrho, vsigma); return true;
    }
    return false;
}

/**
 * @brief The functional at one grid point: sum over the weighted terms (qmxc::evaluateXC per point).
 *
 * Points with rho < rho_min (= effectiveThreshold), rho <= 0 or NaN get ex = ec = vrho = vsigma = 0 and return
 * false. ex/ec: energy per electron of the exchange / correlation terms (weights applied); vrho, vsigma: weighted
 * sums (vsigma only accumulated for gga). Same summation order as the CPU always had (terms in list order).
 */
QM_HD bool evaluatePoint(const XCTerm* terms, int nterms, bool gga, double rho_min, double rho, double sigma,
                         double& ex, double& ec, double& vrho, double& vsigma)
{
    ex = ec = vrho = vsigma = 0.0;
    if (!(rho >= rho_min) || rho <= 0.0) return false;  // also drops NaN and negative round-off densities
    const double sg = gga ? (sigma > 0.0 ? sigma : 0.0) : 0.0;
    for (int t = 0; t < nterms; ++t) {
        const double w = terms[t].weight;
        if (w == 0.0) continue;
        double e = 0.0, v = 0.0, vs = 0.0;
        evaluateKernel(terms[t].kernel, rho, sg, e, v, vs);
        if (isExchange(terms[t].kernel))
            ex += w * e;
        else
            ec += w * e;
        vrho += w * v;
        if (gga) vsigma += w * vs;
    }
    return true;
}

}  // namespace kern
}  // namespace qmxc

namespace qmgrid {
namespace kern {

/// Becke's cell function argument after three iterations of f(x) = 3/2 x - 1/2 x^3 (Becke 1988, eq. 19-20).
QM_HD double beckeIterate(double g)
{
    g = (3.0 - g * g) * g * 0.5;
    g = (3.0 - g * g) * g * 0.5;
    g = (3.0 - g * g) * g * 0.5;
    return g;
}

/// Stratmann-Scuseria-Frisch polynomial (CPL 257, 213 (1996), eq. 14, a = 0.64), as pyscf gen_grid.stratmann.
QM_HD double stratmannIterate(double g)
{
    const double a = 0.64;
    if (g <= -a) return -1.0;
    if (g >= a) return 1.0;
    const double ma = g / a, ma2 = ma * ma;
    return (1.0 / 16.0) * (ma * (35.0 + ma2 * (-35.0 + ma2 * (21.0 - 5.0 * ma2))));
}

/// d h / d nu of the partition polynomial h (beckeIterate). Claude Generated (WP8).
QM_HD double beckeIterateDerivative(double g)
{
    // h = f(f(f(g))), f(x) = (3 - x^2) x / 2, f'(x) = 3/2 (1 - x^2)
    const double g1 = (3.0 - g * g) * g * 0.5;
    const double g2 = (3.0 - g1 * g1) * g1 * 0.5;
    return 1.5 * (1.0 - g2 * g2) * 1.5 * (1.0 - g1 * g1) * 1.5 * (1.0 - g * g);
}

/// d h / d nu of stratmannIterate. Claude Generated (WP8).
QM_HD double stratmannIterateDerivative(double g)
{
    const double a = 0.64;
    if (g <= -a || g >= a) return 0.0;
    const double z = g / a, z2 = z * z;
    return (1.0 / 16.0) * (35.0 + z2 * (-105.0 + z2 * (105.0 - 35.0 * z2))) / a;
}

/// x^k for small non-negative k (k = -1 returns 0; only used multiplied by a zero power).
QM_HD double ipow(double x, int k)
{
    if (k < 0) return 0.0;
    double v = 1.0;
    for (int i = 0; i < k; ++i) v *= x;
    return v;
}

/**
 * @brief One cartesian component x^l y^m z^n sum_k c_k exp(-a_k r^2) and its derivatives at one point.
 *
 * (x, y, z) relative to the shell centre; R0 = sum c e, R1 = sum c a e, R2 = sum c a^2 e (R2 only read for
 * deriv 2). out[0] value; deriv >= 1: out[1..3] d/dx, d/dy, d/dz; deriv >= 2: out[4..9] xx, xy, xz, yy, yz, zz.
 *   d/dx [x^l exp(-a r^2)]  = (l x^(l-1) - 2 a x^(l+1)) exp(-a r^2)
 *   d2/dx2                  = (l (l-1) x^(l-2) - 2 a (2l+1) x^l + 4 a^2 x^(l+2)) exp(-a r^2)
 *   d2/dxdy                 = (l x^(l-1) - 2 a x^(l+1)) (m y^(m-1) - 2 a y^(m+1)) exp(-a r^2), per power of a.
 */
QM_HD void aoComponent(int l, int m, int n, double x, double y, double z, double R0, double R1, double R2, int deriv,
                       double* out)
{
    const double xl = ipow(x, l), ym = ipow(y, m), zn = ipow(z, n);
    out[0] = xl * ym * zn * R0;
    if (deriv >= 1) {
        out[1] = (l * ipow(x, l - 1) * R0 - 2.0 * x * xl * R1) * ym * zn;
        out[2] = (m * ipow(y, m - 1) * R0 - 2.0 * y * ym * R1) * xl * zn;
        out[3] = (n * ipow(z, n - 1) * R0 - 2.0 * z * zn * R1) * xl * ym;
    }
    if (deriv >= 2) {
        const double sx = l * (l - 1) * ipow(x, l - 2) * R0 - 2.0 * (2 * l + 1) * xl * R1 + 4.0 * x * x * xl * R2;
        const double sy = m * (m - 1) * ipow(y, m - 2) * R0 - 2.0 * (2 * m + 1) * ym * R1 + 4.0 * y * y * ym * R2;
        const double sz = n * (n - 1) * ipow(z, n - 2) * R0 - 2.0 * (2 * n + 1) * zn * R1 + 4.0 * z * z * zn * R2;
        // coefficients of a^0 and a^1 of the first derivative factor per direction
        const double ax0 = l * ipow(x, l - 1), ax1 = -2.0 * x * xl;
        const double ay0 = m * ipow(y, m - 1), ay1 = -2.0 * y * ym;
        const double az0 = n * ipow(z, n - 1), az1 = -2.0 * z * zn;
        out[4] = sx * ym * zn;
        out[7] = sy * xl * zn;
        out[9] = sz * xl * ym;
        out[5] = (ax0 * ay0 * R0 + (ax0 * ay1 + ax1 * ay0) * R1 + ax1 * ay1 * R2) * zn;
        out[6] = (ax0 * az0 * R0 + (ax0 * az1 + ax1 * az0) * R1 + ax1 * az1 * R2) * ym;
        out[8] = (ay0 * az0 * R0 + (ay0 * az1 + ay1 * az0) * R1 + ay1 * az1 * R2) * xl;
    }
}

}  // namespace kern
}  // namespace qmgrid

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif
