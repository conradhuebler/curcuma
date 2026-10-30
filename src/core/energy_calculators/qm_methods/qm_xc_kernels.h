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
    LYP_C,     ///< GGA Lee-Yang-Parr correlation, Miehlich form without Laplacian (WP7), libxc GGA_C_LYP
    TPSS_X,    ///< meta-GGA TPSS exchange (WP11 stage 1, CPU only), libxc MGGA_X_TPSS
    TPSS_C,    ///< meta-GGA TPSS correlation (revPKZB + PBE, WP11 stage 1, CPU only), libxc MGGA_C_TPSS
    SCAN_X,    ///< meta-GGA SCAN exchange (WP11 stage 1, CPU only), libxc MGGA_X_SCAN
    SCAN_C     ///< meta-GGA SCAN correlation (WP11 stage 1, CPU only), libxc MGGA_C_SCAN
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
QM_HD bool isExchange(XCKernel k)
{
    return k == XCKernel::SlaterX || k == XCKernel::PBE_X || k == XCKernel::B88_X || k == XCKernel::TPSS_X
        || k == XCKernel::SCAN_X;
}
/// True when the kernel depends on tau = 1/2 sum_i |grad psi_i|^2 (meta-GGA, WP11); such a kernel needs sigma too.
QM_HD bool isMGGA(XCKernel k)
{
    return k == XCKernel::TPSS_X || k == XCKernel::TPSS_C || k == XCKernel::SCAN_X || k == XCKernel::SCAN_C;
}
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
        case XCKernel::TPSS_X:
        case XCKernel::TPSS_C:
        case XCKernel::SCAN_X:
        case XCKernel::SCAN_C:   eps = vrho = 0.0; return false;  // meta-GGA: evaluateKernelMGGA (needs tau)
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

// ===============================================================================================================
// Meta-GGA kernels (WP11 stage 1: TPSS and SCAN, closed shell, energies on the CPU; Claude Generated, Oct 2026)
// ===============================================================================================================
//
// Variables (closed shell, total quantities, as libxc's spin-unpolarised input):
//   rho, sigma = |grad rho|^2, tau = 1/2 sum_i^occ |grad psi_i|^2 summed over both spins
//       = 1/2 sum_mn P_mn grad phi_m . grad phi_n                         (TPSS 2003 eq. 2, summed over sigma)
//   tau_W = sigma / (8 rho),  tau_unif = (3/10) (3 pi^2)^(2/3) rho^(5/3),
//   p = s^2 = sigma / (4 (3 pi^2)^(2/3) rho^(8/3))                        (TPSS eq. 4)
//   z = tau_W / tau,  alpha = (tau - tau_W) / tau_unif                     (TPSS eq. 8; SCAN's alpha)
// libxc evaluates the same functions in per-spin variables (x_s, t_s = tau_s / rho_s^(5/3)); for rho_s = rho/2,
// sigma_ss = sigma/4, tau_s = tau/2 they reduce to the expressions above (checked in ctest qm_mgga_kernel).
//
// Derivatives: forward-mode dual numbers (Dual3 below: value and the three partials d/drho, d/dsigma, d/dtau,
// propagated through every operation by the chain rule). Reason: TPSS correlation nests three PBE evaluations, a max()
// switch and a min() cap, SCAN has piecewise interpolation functions; hand-written chain rules for that are ~200 lines
// in which a sign slip gives plausible but wrong potentials. Dual numbers are not an approximation: they give the
// derivative of the evaluated expression to rounding, from the SAME code that gives the energy (one source). They are
// checked against libxc's v_rho, v_sigma, v_tau and against finite differences (ctest qm_mgga_kernel). The existing
// analytic PBE correlation kernel (pbeCorrelation) enters TPSS correlation unchanged, lifted into a Dual3 by its
// partial derivatives.
//
// libxc thresholds (functionals.c, work_mgga_inc.c of libxc 7.0.0; dens_threshold 1e-15 for all four functionals):
//   tau  -> max(tau, 1e-20)                    (tau_threshold; "many functionals divide by tau")
//   sigma -> min(sigma, 8 rho tau)             (tau >= tau_W, i.e. z <= 1 and alpha >= 0)
//   sigma -> max(sigma, 1e-40)                 (sigma_threshold^2 = dens_threshold^(8/3)), rho < 1e-15 -> 0
// evaluateKernelMGGA applies the first two (they keep the formulas finite: z = tau_W/tau, alpha >= 0 for the SCAN
// interpolation); the derivatives are those at the clamped values, as libxc's (no chain rule through the clamp). The
// sigma floor is not applied (sigma = 0 is handled exactly, see tpssExchange and scanExchange), the density floor is
// kGGADensityFloor (1e-20) plus -qm.xc_density_threshold, as for the GGAs. For a positive semidefinite density matrix
// tau >= tau_W holds exactly (Cauchy-Schwarz), so the clamp only catches rounding.

/// Value and partial derivatives d/drho, d/dsigma, d/dtau (forward-mode automatic differentiation).
struct Dual3 {
    double v;
    double d[3];
};

QM_HD Dual3 dconst(double c)
{
    Dual3 r;
    r.v = c;
    r.d[0] = r.d[1] = r.d[2] = 0.0;
    return r;
}
/// Independent variable number k (0 rho, 1 sigma, 2 tau) with value x.
QM_HD Dual3 dvar(double x, int k)
{
    Dual3 r = dconst(x);
    r.d[k] = 1.0;
    return r;
}
/// f(a) with f(a.v) = fv and f'(a.v) = fp.
QM_HD Dual3 dchain(const Dual3& a, double fv, double fp)
{
    Dual3 r;
    r.v = fv;
    for (int k = 0; k < 3; ++k) r.d[k] = fp * a.d[k];
    return r;
}
QM_HD Dual3 operator+(const Dual3& a, const Dual3& b)
{
    Dual3 r;
    r.v = a.v + b.v;
    for (int k = 0; k < 3; ++k) r.d[k] = a.d[k] + b.d[k];
    return r;
}
QM_HD Dual3 operator-(const Dual3& a, const Dual3& b)
{
    Dual3 r;
    r.v = a.v - b.v;
    for (int k = 0; k < 3; ++k) r.d[k] = a.d[k] - b.d[k];
    return r;
}
QM_HD Dual3 operator-(const Dual3& a)
{
    Dual3 r;
    r.v = -a.v;
    for (int k = 0; k < 3; ++k) r.d[k] = -a.d[k];
    return r;
}
QM_HD Dual3 operator*(const Dual3& a, const Dual3& b)
{
    Dual3 r;
    r.v = a.v * b.v;
    for (int k = 0; k < 3; ++k) r.d[k] = a.d[k] * b.v + a.v * b.d[k];
    return r;
}
QM_HD Dual3 operator/(const Dual3& a, const Dual3& b)
{
    Dual3 r;
    r.v = a.v / b.v;
    for (int k = 0; k < 3; ++k) r.d[k] = (a.d[k] - r.v * b.d[k]) / b.v;
    return r;
}
QM_HD Dual3 operator+(const Dual3& a, double c) { return dchain(a, a.v + c, 1.0); }
QM_HD Dual3 operator+(double c, const Dual3& a) { return dchain(a, c + a.v, 1.0); }
QM_HD Dual3 operator-(const Dual3& a, double c) { return dchain(a, a.v - c, 1.0); }
QM_HD Dual3 operator-(double c, const Dual3& a) { return dchain(a, c - a.v, -1.0); }
QM_HD Dual3 operator*(const Dual3& a, double c) { return dchain(a, a.v * c, c); }
QM_HD Dual3 operator*(double c, const Dual3& a) { return dchain(a, c * a.v, c); }
QM_HD Dual3 operator/(const Dual3& a, double c) { return dchain(a, a.v / c, 1.0 / c); }
QM_HD Dual3 operator/(double c, const Dual3& a)
{
    const double q = c / a.v;
    return dchain(a, q, -q / a.v);
}
QM_HD Dual3 dexp(const Dual3& a)
{
    const double e = ::exp(a.v);
    return dchain(a, e, e);
}
QM_HD Dual3 dlog(const Dual3& a) { return dchain(a, ::log(a.v), 1.0 / a.v); }
QM_HD Dual3 dsqrt(const Dual3& a)
{
    const double s = ::sqrt(a.v);
    return dchain(a, s, 0.5 / s);
}
QM_HD Dual3 dcbrt(const Dual3& a)
{
    const double c = ::cbrt(a.v);
    return dchain(a, c, c / (3.0 * a.v));
}

// Constants of the meta-GGAs (libxc 7.0.0 parameter tables, which follow the papers).
// TPSS exchange (Tao, Perdew, Staroverov, Scuseria, PRL 91, 146401 (2003)): b (eq. 7), c, e (eq. 10), kappa (eq. 5),
// mu (large-p limit, text below eq. 9); libxc mgga_x_tpss.c tpss_values.
constexpr double kTPSS_b = 0.40, kTPSS_c = 1.59096, kTPSS_e = 1.537, kTPSS_kappa = 0.8040, kTPSS_mu = 0.21951;
// TPSS correlation: d (eq. 11, Eh^-1), C(0, 0) = 0.53 (eq. 13; closed shell zeta = 0 and xi = 0, eq. 14);
// beta of the PBE parts = kPBE_beta (libxc mgga_c_tpss.c, constant beta as PBE, not the r_s-dependent one of revTPSS).
constexpr double kTPSS_d = 2.8, kTPSS_C00 = 0.53;
// PW92 ferromagnetic set (zeta = 1) of libxc's LDA_C_PW_MOD (PW92 table I, A with the extra digits of PW92's code).
constexpr double kPWF_A = 0.01554535, kPWF_alpha1 = 0.20548;
constexpr double kPWF_beta1 = 14.1189, kPWF_beta2 = 6.1977, kPWF_beta3 = 3.3662, kPWF_beta4 = 0.62517;
// SCAN (Sun, Ruzsinszky, Perdew, PRL 115, 036402 (2015) and its Supplemental Material; libxc mgga_x_scan.c par_scan,
// mgga_x_scan.mpl, mgga_c_scan.mpl, gga_c_scan_e0.mpl). Exchange interpolation f_x(alpha): c1x, c2x, dx; h1x: k1;
// h0x, a1 (g_x). Correlation f_c(alpha): c1c, c2c, dc; eps_c^0: b1c, b2c, b3c, chi_inf; eps_c^1: beta(r_s) =
// beta_a (1 + 0.1 r_s) / (1 + 0.1778 r_s).
constexpr double kSCAN_c1x = 0.667, kSCAN_c2x = 0.8, kSCAN_dx = 1.24, kSCAN_k1 = 0.065;
constexpr double kSCAN_h0x = 1.174, kSCAN_a1 = 4.9479;
constexpr double kSCAN_c1c = 0.64, kSCAN_c2c = 1.5, kSCAN_dc = 0.7;
constexpr double kSCAN_b1c = 0.0285764, kSCAN_b2c = 0.0889, kSCAN_b3c = 0.125541, kSCAN_chi = 0.12802585262625815;
constexpr double kSCAN_beta_b = 0.1, kSCAN_beta_c = 0.1778;
constexpr double kMU_GE = 10.0 / 81.0;                 // gradient-expansion coefficient of exchange
constexpr double kDblEpsilon = 2.220446049250313e-16;  // DBL_EPSILON, used by libxc's SCAN interpolation cut-offs
constexpr double kMGGATauFloor = 1e-20;                // libxc tau_threshold

/// PW92 G(r_s) (PW92 eq. 10) for a parameter set, as a dual number.
QM_HD Dual3 pw92GDual(const Dual3& rs, double A, double a1, double b1, double b2, double b3, double b4)
{
    const Dual3 srs = dsqrt(rs);
    const Dual3 Q1 = 2.0 * A * (b1 * srs + b2 * rs + b3 * rs * srs + b4 * rs * rs);
    return -2.0 * A * (1.0 + a1 * rs) * dlog(1.0 + 1.0 / Q1);
}

/// Common meta-GGA variables: p = s^2, alpha = (tau - tau_W) / tau_unif, z = tau_W / tau (file comment above).
QM_HD void mggaVariables(const Dual3& R, const Dual3& S, const Dual3& T, Dual3& p, Dual3& alpha, Dual3& z)
{
    const double cp = 1.0 / (4.0 * ::cbrt(9.0 * kPi * kPi * kPi * kPi));  // 1 / (4 (3 pi^2)^(2/3))
    const double ctu = 0.3 * ::cbrt(9.0 * kPi * kPi * kPi * kPi);         // (3/10) (3 pi^2)^(2/3)
    const Dual3 R13 = dcbrt(R);
    const Dual3 R53 = R * R13 * R13;
    const Dual3 tauW = S / (8.0 * R);
    p = cp * S / (R53 * R);
    alpha = (T - tauW) / (ctu * R53);
    z = tauW / T;
}

/// eps (per electron) and v_rho = d(rho eps)/d rho, v_sigma = rho d eps/d sigma, v_tau = rho d eps/d tau.
QM_HD void mggaOutputs(double rho, const Dual3& eps, double& e, double& vrho, double& vsigma, double& vtau)
{
    e = eps.v;
    vrho = eps.v + rho * eps.d[0];
    vsigma = rho * eps.d[1];
    vtau = rho * eps.d[2];
}

/// TPSS exchange (PRL 91, 146401, eqs. 3-10): eps_x = eps_x^unif(rho) F_x(p, z), closed shell,
///   F_x = 1 + kappa - kappa / (1 + x/kappa)                                                          (eq. 5)
///   q_b = (9/20) (alpha - 1) / sqrt(1 + b alpha (alpha - 1)) + 2p/3                                    (eq. 7)
///   x = { [10/81 + c z^2 / (1 + z^2)^2] p + 146/2025 q_b^2 - 73/405 q_b sqrt(1/2 (9/25 z^2 + p^2))
///         + (1/kappa) (10/81)^2 p^2 + 2 sqrt(e) (10/81) (9/25) z^2 + e mu p^3 } / (1 + sqrt(e) p)^2       (eq. 10)
/// z and p are both proportional to sigma, so sqrt(1/2 (9/25 z^2 + p^2)) = sigma sqrt(1/2 (9/25 (z/sigma)^2 +
/// (p/sigma)^2)) is written with sigma outside the root: finite derivative at sigma = 0 (the root of zero would give
/// 0 * inf in the dual numbers). Claude Generated (WP11).
QM_HD Dual3 tpssExchangeDual(const Dual3& R, const Dual3& S, const Dual3& T)
{
    Dual3 p, alpha, z;
    mggaVariables(R, S, T, p, alpha, z);
    const double cp = 1.0 / (4.0 * ::cbrt(9.0 * kPi * kPi * kPi * kPi));
    const Dual3 R13 = dcbrt(R);
    const Dual3 zs = 1.0 / (8.0 * R * T);        // z / sigma
    const Dual3 ps = cp / (R * R * R13 * R13);   // p / sigma
    const Dual3 root = S * dsqrt(0.5 * (9.0 / 25.0 * zs * zs + ps * ps));
    const Dual3 qb = 9.0 / 20.0 * (alpha - 1.0) / dsqrt(1.0 + kTPSS_b * alpha * (alpha - 1.0)) + 2.0 / 3.0 * p;
    const Dual3 z2 = z * z;
    const Dual3 opz2 = 1.0 + z2;
    const double se = ::sqrt(kTPSS_e);
    const Dual3 num = (kMU_GE + kTPSS_c * z2 / (opz2 * opz2)) * p + 146.0 / 2025.0 * qb * qb
        - 73.0 / 405.0 * qb * root + kMU_GE * kMU_GE / kTPSS_kappa * p * p + 2.0 * se * kMU_GE * 9.0 / 25.0 * z2
        + kTPSS_e * kTPSS_mu * p * p * p;
    const Dual3 den1 = 1.0 + se * p;
    const Dual3 x = num / (den1 * den1);
    const Dual3 Fx = 1.0 + kTPSS_kappa * (1.0 - kTPSS_kappa / (kTPSS_kappa + x));
    const Dual3 eunif = -0.75 * dcbrt(3.0 / kPi * R);
    return eunif * Fx;
}

/// PBE correlation of a fully spin-polarised density n = rho/2 with |grad n|^2 = sigma/4 (zeta = 1, phi = 2^(-1/3),
/// phi^3 = 1/2): eps = eps_c^PW92(r_s(n), zeta = 1) + gamma phi^3 ln(1 + beta/gamma t^2 (1 + A t^2) / (1 + A t^2 +
/// A^2 t^4)), A = beta / (gamma (exp(-eps_c^PW92 / (gamma phi^3)) - 1)), t^2 = x_s^2 / (16 r_s(n)) with
/// x_s^2 = |grad n|^2 / n^(8/3) (PBE 1996 eqs. 3, 7, 8 with phi(zeta); libxc's tt() with 2^(2/3) phi^2 = 1).
/// The eps~_c^sigma = eps_c^PBE(n_sigma, 0, grad n_sigma, 0) of TPSS eq. 12. Claude Generated (WP11).
QM_HD Dual3 pbeCorrelationPolarisedDual(const Dual3& R, const Dual3& S)
{
    const double gamma = (1.0 - ::log(2.0)) / (kPi * kPi);
    const double phi3 = 0.5;
    const Dual3 n = 0.5 * R;
    const Dual3 rs = dcbrt(3.0 / (4.0 * kPi) / n);
    const Dual3 ec = pw92GDual(rs, kPWF_A, kPWF_alpha1, kPWF_beta1, kPWF_beta2, kPWF_beta3, kPWF_beta4);
    const Dual3 n13 = dcbrt(n);
    const Dual3 xs2 = 0.25 * S / (n * n * n13 * n13);
    const Dual3 t2 = xs2 / (16.0 * rs);
    const Dual3 A = kPBE_beta / (gamma * (dexp(-ec / (gamma * phi3)) - 1.0));
    const Dual3 f1 = t2 + A * t2 * t2;
    const Dual3 f2 = kPBE_beta * f1 / (gamma * (1.0 + A * f1));
    return ec + gamma * phi3 * dlog(1.0 + f2);
}

/// TPSS correlation (PRL 91, 146401, eqs. 11-14), closed shell (zeta = 0, xi = 0, so C = C(0, 0) = 0.53):
///   eps_c^revPKZB = eps_c^PBE(rho, sigma) [1 + C z^2] - (1 + C) z^2 eps~_c                              (eq. 12)
///   eps~_c = max(eps_c^PBE(rho/2, 0, sigma/4, 0), eps_c^PBE(rho, sigma))   (both spins equal, n_s/n = 1/2 each)
///   eps_c^TPSS = eps_c^revPKZB [1 + d eps_c^revPKZB z^3],  z = min(tau_W / tau, 1)                    (eq. 11)
/// eps_c^PBE(rho, sigma) is the WP6 kernel pbeCorrelation (analytic derivatives), lifted into a Dual3:
///   d eps/d rho = (v_rho - eps) / rho, d eps/d sigma = v_sigma / rho. max() and min() take the derivative of the
/// chosen branch (libxc's m_max/m_min). Claude Generated (WP11).
QM_HD Dual3 tpssCorrelationDual(const Dual3& R, const Dual3& S, const Dual3& T)
{
    double e = 0.0, vr = 0.0, vs = 0.0;
    pbeCorrelation(R.v, S.v, e, vr, vs);
    Dual3 epbe;
    epbe.v = e;
    for (int k = 0; k < 3; ++k) epbe.d[k] = (vr - e) / R.v * R.d[k] + vs / R.v * S.d[k];
    const Dual3 epol = pbeCorrelationPolarisedDual(R, S);
    const Dual3 etil = epol.v > epbe.v ? epol : epbe;
    const Dual3 zraw = S / (8.0 * R * T);
    const Dual3 z = zraw.v > 1.0 ? dconst(1.0) : zraw;
    const Dual3 z2 = z * z;
    const Dual3 f0 = epbe * (1.0 + kTPSS_C00 * z2) - (1.0 + kTPSS_C00) * z2 * etil;
    return f0 * (1.0 + kTPSS_d * f0 * z2 * z);
}

/// SCAN interpolation function f(alpha) (SCAN exchange f_x with c1x, c2x, dx; correlation f_c with c1c, c2c, dc):
///   alpha <= 1: exp(-c1 alpha / (1 - alpha)),   alpha > 1: -d exp(c2 / (1 - alpha)),
/// cut to 0 where the exponential falls below DBL_EPSILON, as libxc (left: alpha > L/(L + c1), L = -ln DBL_EPSILON;
/// right: alpha < (L_d + c2)/L_d, L_d = -ln(DBL_EPSILON/d)); f(1) = 0. Claude Generated (WP11).
QM_HD Dual3 scanInterpolation(const Dual3& a, double c1, double c2, double d)
{
    if (a.v <= 1.0) {
        const double L = -::log(kDblEpsilon);
        if (a.v > L / (L + c1)) return dconst(0.0);
        return dexp(-c1 * a / (1.0 - a));
    }
    const double Ld = -::log(kDblEpsilon / d);
    if (a.v < (Ld + c2) / Ld) return dconst(0.0);
    return -d * dexp(c2 / (1.0 - a));
}

/// SCAN exchange (Sun, Ruzsinszky, Perdew 2015; formulas of libxc mgga_x_scan.mpl), closed shell:
///   F_x = [h1x(y) (1 - f_x(alpha)) + h0x f_x(alpha)] g_x(p),  h1x = 1 + k1 - k1 / (1 + y/k1),
///   y = mu_GE p + b4 p^2 exp(-b4 p / mu_GE) + [b1 p + b2 (1 - alpha) exp(-b3 (1 - alpha)^2)]^2,
///   b2 = sqrt(5913/405000), b1 = (511/13500) / (2 b2), b3 = 1/2, b4 = mu_GE^2 / k1 - 1606/18225 - b1^2,
///   g_x = 1 - exp(-a1 / p^(1/4)),  eps_x = eps_x^unif F_x.
/// g_x = 1 to double precision for a1 / p^(1/4) > 700 (also p = 0): taken as the constant 1 there (the dual chain
/// rule would form 0 * inf). Claude Generated (WP11).
QM_HD Dual3 scanExchangeDual(const Dual3& R, const Dual3& S, const Dual3& T)
{
    Dual3 p, alpha, z;
    mggaVariables(R, S, T, p, alpha, z);
    const double b2 = ::sqrt(5913.0 / 405000.0);
    const double b1 = (511.0 / 13500.0) / (2.0 * b2);
    const double b3 = 0.5;
    const double b4 = kMU_GE * kMU_GE / kSCAN_k1 - 1606.0 / 18225.0 - b1 * b1;
    const Dual3 oma = 1.0 - alpha;
    const Dual3 w = b1 * p + b2 * oma * dexp(-b3 * oma * oma);
    const Dual3 y = kMU_GE * p + b4 * p * p * dexp(-b4 * p / kMU_GE) + w * w;
    const Dual3 h1 = 1.0 + kSCAN_k1 * (1.0 - kSCAN_k1 / (kSCAN_k1 + y));
    const Dual3 fa = scanInterpolation(alpha, kSCAN_c1x, kSCAN_c2x, kSCAN_dx);
    Dual3 gx = dconst(1.0);
    if (p.v > 0.0 && kSCAN_a1 / ::sqrt(::sqrt(p.v)) <= 700.0) gx = 1.0 - dexp(-kSCAN_a1 / dsqrt(dsqrt(p)));
    const Dual3 Fx = (h1 * (1.0 - fa) + kSCAN_h0x * fa) * gx;
    return -0.75 * dcbrt(3.0 / kPi * R) * Fx;
}

/// SCAN correlation (Sun, Ruzsinszky, Perdew 2015; libxc mgga_c_scan.mpl with gga_c_scan_e0.mpl), closed shell
/// (zeta = 0, phi = 1, G_c = 1):
///   eps_c = eps_c^1 + f_c(alpha) (eps_c^0 - eps_c^1)
///   eps_c^1 = eps_c^PW92 + gamma ln(1 + w1 (1 - g(A t^2))),  w1 = exp(-eps_c^PW92 / gamma) - 1,
///             A = beta(r_s) / (gamma w1),  g(y) = (1 + 4 y)^(-1/4),  beta(r_s) = beta (1 + 0.1 r_s)/(1 + 0.1778 r_s),
///             t^2 = sigma / (16 2^(2/3) r_s rho^(8/3))  (libxc tt(), phi = 1); written as libxc's
///             beta (1 - g) / (gamma A) for w1 (1 - g)
///   eps_c^0 = eps_c^LDA0 + b1c ln(1 + w0 (1 - g_inf(s))),  eps_c^LDA0 = -b1c / (1 + b2c sqrt(r_s) + b3c r_s),
///             w0 = exp(-eps_c^LDA0 / b1c) - 1,  g_inf(s) = (1 + 4 chi_inf s^2)^(-1/4),  s^2 = p.
/// Claude Generated (WP11).
QM_HD Dual3 scanCorrelationDual(const Dual3& R, const Dual3& S, const Dual3& T)
{
    Dual3 p, alpha, z;
    mggaVariables(R, S, T, p, alpha, z);
    const double gamma = (1.0 - ::log(2.0)) / (kPi * kPi);
    const Dual3 rs = dcbrt(3.0 / (4.0 * kPi) / R);
    const Dual3 ec = pw92GDual(rs, kPW_A, kPW_alpha1, kPW_beta1, kPW_beta2, kPW_beta3, kPW_beta4);
    const Dual3 beta = kPBE_beta * (1.0 + kSCAN_beta_b * rs) / (1.0 + kSCAN_beta_c * rs);
    const Dual3 w1 = dexp(-ec / gamma) - 1.0;
    const Dual3 A = beta / (gamma * w1);
    const Dual3 R13 = dcbrt(R);
    const Dual3 t2 = S / (R * R * R13 * R13) / (16.0 * ::cbrt(4.0) * rs);
    const Dual3 g = 1.0 / dsqrt(dsqrt(1.0 + 4.0 * A * t2));
    const Dual3 e1 = ec + gamma * dlog(1.0 + beta * (1.0 - g) / (gamma * A));
    const Dual3 el0 = -kSCAN_b1c / (1.0 + kSCAN_b2c * dsqrt(rs) + kSCAN_b3c * rs);
    const Dual3 w0 = dexp(-el0 / kSCAN_b1c) - 1.0;
    const Dual3 ginf = 1.0 / dsqrt(dsqrt(1.0 + 4.0 * kSCAN_chi * p));
    const Dual3 e0 = el0 + kSCAN_b1c * dlog(1.0 + w0 * (1.0 - ginf));
    const Dual3 fc = scanInterpolation(alpha, kSCAN_c1c, kSCAN_c2c, kSCAN_dc);
    return e1 + fc * (e0 - e1);
}

/**
 * @brief One meta-GGA kernel at one point: eps, v_rho = d(rho eps)/d rho, v_sigma = d(rho eps)/d sigma,
 * v_tau = d(rho eps)/d tau (closed shell, total rho, sigma, tau). Applies libxc's tau floor and sigma <= 8 rho tau
 * (file comment above). False (all zero) for a kernel that is not a meta-GGA. Claude Generated (WP11).
 */
QM_HD bool evaluateKernelMGGA(XCKernel k, double rho, double sigma, double tau, double& eps, double& vrho,
                              double& vsigma, double& vtau)
{
    eps = vrho = vsigma = vtau = 0.0;
    if (!isMGGA(k)) return false;
    const double t = tau > kMGGATauFloor ? tau : kMGGATauFloor;
    double s = sigma > 0.0 ? sigma : 0.0;
    if (s > 8.0 * rho * t) s = 8.0 * rho * t;
    const Dual3 R = dvar(rho, 0), S = dvar(s, 1), T = dvar(t, 2);
    Dual3 e = dconst(0.0);
    switch (k) {
        case XCKernel::TPSS_X: e = tpssExchangeDual(R, S, T); break;
        case XCKernel::TPSS_C: e = tpssCorrelationDual(R, S, T); break;
        case XCKernel::SCAN_X: e = scanExchangeDual(R, S, T); break;
        case XCKernel::SCAN_C: e = scanCorrelationDual(R, S, T); break;
        default: return false;
    }
    mggaOutputs(rho, e, eps, vrho, vsigma, vtau);
    return true;
}

/**
 * @brief A functional with meta-GGA terms at one grid point (the evaluatePoint of a meta-GGA functional).
 *
 * Meta-GGA terms through evaluateKernelMGGA, others (LDA/GGA, e.g. a later hybrid's parts) through evaluateKernel
 * (v_tau = 0). Points with rho < rho_min, rho <= 0 or NaN give zeros and false. Claude Generated (WP11).
 */
QM_HD bool evaluatePointMGGA(const XCTerm* terms, int nterms, double rho_min, double rho, double sigma, double tau,
                             double& ex, double& ec, double& vrho, double& vsigma, double& vtau)
{
    ex = ec = vrho = vsigma = vtau = 0.0;
    if (!(rho >= rho_min) || rho <= 0.0) return false;
    const double sg = sigma > 0.0 ? sigma : 0.0;
    for (int i = 0; i < nterms; ++i) {
        const double w = terms[i].weight;
        if (w == 0.0) continue;
        double e = 0.0, v = 0.0, vs = 0.0, vt = 0.0;
        if (isMGGA(terms[i].kernel))
            evaluateKernelMGGA(terms[i].kernel, rho, sg, tau, e, v, vs, vt);
        else
            evaluateKernel(terms[i].kernel, rho, sg, e, v, vs);
        if (isExchange(terms[i].kernel))
            ex += w * e;
        else
            ec += w * e;
        vrho += w * v;
        vsigma += w * vs;
        vtau += w * vt;
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
