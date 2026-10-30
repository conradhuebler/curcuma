/*
 * <Native QM: exchange-correlation functionals of closed-shell Kohn-Sham DFT (WP5: LDA, WP6: PBE, WP7: B3LYP)>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Formulas and literature: qm_xc.h. Checked point by point against libxc through PySCF
 * (pyscf.dft.libxc.eval_xc 'LDA_X', 'LDA_C_VWN', 'LDA_C_VWN_RPA'; ctest qm_lda_kernel; 'GGA_X_PBE', 'GGA_C_PBE',
 * 'LDA_C_PW_MOD'; ctest qm_pbe_kernel; 'GGA_X_B88', 'GGA_C_LYP'; ctest qm_b3lyp_kernel).
 *
 * xcDFT (TCCM winter school 2019) has a Slater exchange (lda_exchange_*.f90) and a "correlation" whose energy
 * formula Ec = -a sum rho w / (1 + b rho^(-1/3) w) puts the grid weight into the denominator; that is not a
 * functional of the density. It is not ported; VWN replaces it.
 *
 * Claude Generated (Oct 2026, WP5; WP6 PBE; WP7 B88, LYP, B3LYP). AI-generated, machine-tested only.
 *
 * This program is free software under GPL-3.0
 */

#include "qm_xc.h"

#include <Eigen/Dense>
#include <fmt/format.h>

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <limits>
#include <vector>

namespace qmxc {

namespace {
constexpr double kPi = 3.14159265358979323846;

/// Paramagnetic VWN parameters {A, x0, b, c} (Hartree), libxc lda_c_vwn.c / VWN 1980 table.
struct VWNParams {
    double A, x0, b, c;
};
constexpr VWNParams kVWN5{ 0.0310907, -0.10498, 3.72744, 12.9352 };
constexpr VWNParams kVWNRPA{ 0.0310907, -0.409286, 13.0720, 42.7198 };

// PW92 paramagnetic correlation, libxc LDA_C_PW_MOD (lda_c_pw.c): A with the extra digits of PW92's own Fortran
// code (the paper prints 0.031091), p = 1 (the r_s^(p+1) term is beta4 r_s^2).
constexpr double kPW_A = 0.0310907, kPW_alpha1 = 0.21370;
constexpr double kPW_beta1 = 7.5957, kPW_beta2 = 3.5876, kPW_beta3 = 1.6382, kPW_beta4 = 0.49294;

// PBE (libxc gga_x_pbe.c / gga_c_pbe.c). mu = beta pi^2 / 3 with the PBE beta (PBE 1996 below eq. 12).
constexpr double kPBE_kappa = 0.8040;
constexpr double kPBE_mu = 0.2195149727645171;
constexpr double kPBE_beta = 0.06672455060314922;
const double kPBE_gamma = (1.0 - std::log(2.0)) / (kPi * kPi);  // (1 - ln 2) / pi^2, PBE 1996 below eq. 7

// B88 (Becke 1988, libxc gga_x_b88.c): beta, gamma = 6; C_x = (3/2) (3/(4 pi))^(1/3) is the per-spin Slater constant.
constexpr double kB88_beta = 0.0042;
// LYP (Lee-Yang-Parr 1988, libxc gga_c_lyp.c): a, b, c, d.
constexpr double kLYP_a = 0.04918, kLYP_b = 0.132, kLYP_c = 0.2533, kLYP_d = 0.349;

// Densities below this are never handed to a GGA kernel, whatever -qm.xc_density_threshold says (also 0): near
// rho ~ 1e-27 ln(1 + 1/Q1) rounds to 0, A = inf and the PBE formulas give NaN; rho^(8/3) underflows below ~1e-116.
// Far below any physical threshold (libxc: 1e-15 exchange, 1e-12 correlation). LDA kernels have no such floor.
constexpr double kGGADensityFloor = 1e-20;

/// Density below which a point is left out of E_xc and V_xc for functional f.
double effectiveThreshold(const XCFunctional& f, double rho_threshold)
{
    return f.needsGradient() ? std::max(rho_threshold, kGGADensityFloor) : rho_threshold;
}

/// PW92 G(r_s) of the paramagnetic gas and d G / d r_s (PW92 eq. 10 and its derivative, eq. A3 there).
///   G = -2A (1 + alpha1 r_s) ln(1 + 1/Q1),  Q1 = 2A (beta1 r_s^(1/2) + beta2 r_s + beta3 r_s^(3/2) + beta4 r_s^2)
///   dG/dr_s = -2A alpha1 ln(1 + 1/Q1) + 2A (1 + alpha1 r_s) Q1' / (Q1^2 + Q1)
/// ln(1 + x) rather than log1p(x), as libxc's generated code: at r_s ~ 6e3 (rho 1e-12) log1p would differ from libxc
/// by 1e-10 relative (libxc's rounding), everywhere else both agree to 1e-15.
void pw92G(double rs, double& G, double& dG)
{
    const double srs = std::sqrt(rs);
    const double Q1 = 2.0 * kPW_A * (kPW_beta1 * srs + kPW_beta2 * rs + kPW_beta3 * rs * srs + kPW_beta4 * rs * rs);
    const double Q1p = kPW_A * (kPW_beta1 / srs + 2.0 * kPW_beta2 + 3.0 * kPW_beta3 * srs + 4.0 * kPW_beta4 * rs);
    const double L = std::log(1.0 + 1.0 / Q1);
    G = -2.0 * kPW_A * (1.0 + kPW_alpha1 * rs) * L;
    dG = -2.0 * kPW_A * kPW_alpha1 * L + 2.0 * kPW_A * (1.0 + kPW_alpha1 * rs) * Q1p / (Q1 * Q1 + Q1);
}
}  // namespace

bool parseVWNVariant(const std::string& name, VWNVariant& v)
{
    std::string s = name;
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::tolower(c); });
    if (s == "vwn5" || s == "vwn") {
        v = VWNVariant::VWN5;
        return true;
    }
    if (s == "vwn_rpa" || s == "vwn-rpa" || s == "rpa" || s == "vwn3") {
        v = VWNVariant::VWN_RPA;
        return true;
    }
    return false;
}

std::string vwnName(VWNVariant v) { return v == VWNVariant::VWN5 ? "VWN5" : "VWN-RPA"; }

bool isExchangeKernel(XCKernel k) { return k == XCKernel::SlaterX || k == XCKernel::PBE_X || k == XCKernel::B88_X; }

bool isGGAKernel(XCKernel k)
{
    return k == XCKernel::PBE_X || k == XCKernel::PBE_C || k == XCKernel::B88_X || k == XCKernel::LYP_C;
}

std::string kernelName(XCKernel k)
{
    switch (k) {
        case XCKernel::SlaterX:  return "Slater";
        case XCKernel::VWN5_C:   return "VWN5";
        case XCKernel::VWNRPA_C: return "VWN-RPA";
        case XCKernel::PW92_C:   return "PW92";
        case XCKernel::PBE_X:    return "PBE-X";
        case XCKernel::PBE_C:    return "PBE-C";
        case XCKernel::B88_X:    return "B88";
        case XCKernel::LYP_C:    return "LYP";
    }
    return "?";
}

XCFunctional::XCFunctional() = default;
XCFunctional::XCFunctional(const XCFunctional&) = default;
XCFunctional& XCFunctional::operator=(const XCFunctional&) = default;
XCFunctional::~XCFunctional() = default;

bool XCFunctional::needsGradient() const
{
    for (const XCTerm& t : terms)
        if (isGGAKernel(t.kernel) && t.weight != 0.0) return true;
    return false;
}

std::string XCFunctional::label() const
{
    if (!hasXC()) return "HF (no density functional)";
    if (!description.empty()) return description;
    std::string s;
    for (const XCTerm& t : terms) s += (s.empty() ? "" : " + ") + fmt::format("{:g} {}", t.weight, kernelName(t.kernel));
    return s;
}

XCFunctional ldaFunctional(VWNVariant v)
{
    XCFunctional f;
    f.name = "lda";
    f.exact_exchange = 0.0;
    f.terms = { { XCKernel::SlaterX, 1.0 }, { v == VWNVariant::VWN5 ? XCKernel::VWN5_C : XCKernel::VWNRPA_C, 1.0 } };
    f.description = std::string("SVWN") + (v == VWNVariant::VWN5 ? "5" : "-RPA") + " (Slater exchange + " + vwnName(v)
        + " correlation)";
    return f;
}

XCFunctional pbeFunctional()
{
    XCFunctional f;
    f.name = "pbe";
    f.exact_exchange = 0.0;
    f.terms = { { XCKernel::PBE_X, 1.0 }, { XCKernel::PBE_C, 1.0 } };
    f.description = "PBE-GGA (PBE exchange + PBE correlation, PW92 LDA part)";
    return f;
}

// WP7: B3LYP (Stephens et al. 1994) = 0.08 Slater + 0.72 B88 + 0.20 HF + 0.19 VWN + 0.81 LYP. The VWN variant is a
// choice of the caller (QMEngine: -qm.vwn_variant, default vwn5 as ORCA's B3LYP); it is always named in the label.
XCFunctional b3lypFunctional(VWNVariant v)
{
    XCFunctional f;
    f.name = "b3lyp";
    f.exact_exchange = 0.20;
    f.terms = { { XCKernel::SlaterX, 0.08 }, { XCKernel::B88_X, 0.72 },
                { v == VWNVariant::VWN5 ? XCKernel::VWN5_C : XCKernel::VWNRPA_C, 0.19 }, { XCKernel::LYP_C, 0.81 } };
    f.description = std::string("B3LYP-") + vwnName(v) + " (0.08 Slater + 0.72 B88 + 0.20 HF exchange, 0.19 "
        + vwnName(v) + " + 0.81 LYP correlation)";
    return f;
}

// Dirac (1930) / Slater (1951): e_x(rho) = rho eps_x = -(3/4) (3/pi)^(1/3) rho^(4/3), closed shell.
void slaterExchange(double rho, double& eps, double& v)
{
    const double t = std::cbrt(3.0 * rho / kPi);  // (3 rho / pi)^(1/3)
    eps = -0.75 * t;
    v = -t;  // (4/3) eps_x
}

// Vosko-Wilk-Nusair (1980) eq. 4.4 with the paramagnetic parameters; derivative with respect to x = sqrt(r_s):
//   d eps/dx = A [ 2/x - (2x+b)/X - b/X - b x0/X0 ( 2/(x-x0) - (2x+b)/X - (b+2 x0)/X ) ]
// (the atan term gives -4b/((2x+b)^2 + Q^2) = -b/X since (2x+b)^2 + Q^2 = 4X), and
//   v = d(rho eps)/d rho = eps - (x/6) d eps/dx      (dx/d rho = -x/(6 rho)).
void vwnCorrelation(double rho, VWNVariant variant, double& eps, double& v)
{
    const VWNParams& p = variant == VWNVariant::VWN5 ? kVWN5 : kVWNRPA;
    const double rs = std::cbrt(3.0 / (4.0 * kPi * rho));
    const double x = std::sqrt(rs);
    const double X = x * x + p.b * x + p.c;
    const double X0 = p.x0 * p.x0 + p.b * p.x0 + p.c;
    const double Q = std::sqrt(4.0 * p.c - p.b * p.b);
    const double at = std::atan(Q / (2.0 * x + p.b));
    const double bx0 = p.b * p.x0 / X0;
    eps = p.A * (std::log(x * x / X) + 2.0 * p.b / Q * at
                 - bx0 * (std::log((x - p.x0) * (x - p.x0) / X) + 2.0 * (p.b + 2.0 * p.x0) / Q * at));
    const double dX = (2.0 * x + p.b) / X;
    const double deps_dx = p.A * (2.0 / x - dX - p.b / X - bx0 * (2.0 / (x - p.x0) - dX - (p.b + 2.0 * p.x0) / X));
    v = eps - x / 6.0 * deps_dx;
}

// PW92 (1992) paramagnetic correlation as an LDA kernel: v = eps - (r_s/3) d eps/d r_s  (d r_s/d rho = -r_s/(3 rho)).
void pw92Correlation(double rho, double& eps, double& v)
{
    const double rs = std::cbrt(3.0 / (4.0 * kPi * rho));
    double dG = 0.0;
    pw92G(rs, eps, dG);
    v = eps - rs / 3.0 * dG;
}

// PBE exchange (PBE 1996 eqs. 10 and 14), closed shell, e_x = rho eps_x^LDA(rho) F_x(s^2):
//   s^2 = c sigma,  c = 1 / (4 (3 pi^2)^(2/3) rho^(8/3)),  F_x' = dF_x/d(s^2) = mu / (1 + mu s^2/kappa)^2
//   v_rho   = d e_x/d rho   = eps_x^LDA [ (4/3) F_x - (8/3) s^2 F_x' ]      (d s^2/d rho = -(8/3) s^2/rho)
//   v_sigma = d e_x/d sigma = rho eps_x^LDA F_x' c
void pbeExchange(double rho, double sigma, double& eps, double& vrho, double& vsigma)
{
    const double eps_lda = -0.75 * std::cbrt(3.0 * rho / kPi);
    const double rho13 = std::cbrt(rho);
    const double c = 1.0 / (4.0 * std::pow(3.0 * kPi * kPi, 2.0 / 3.0) * rho * rho * rho13 * rho13);
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
void pbeCorrelation(double rho, double sigma, double& eps, double& vrho, double& vsigma)
{
    const double rs = std::cbrt(3.0 / (4.0 * kPi * rho));
    double ec = 0.0, dec = 0.0;
    pw92G(rs, ec, dec);
    const double kF = std::cbrt(3.0 * kPi * kPi * rho);
    const double u = sigma * kPi / (16.0 * kF * rho * rho);
    const double E = std::exp(-ec / kPBE_gamma);
    const double A = (kPBE_beta / kPBE_gamma) / (E - 1.0);
    const double z = A * u;
    const double D = 1.0 + z + z * z;
    const double y = u * ((1.0 + z) / D);
    const double arg = 1.0 + kPBE_beta / kPBE_gamma * y;
    const double H = kPBE_gamma * std::log(arg);
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
void b88Exchange(double rho, double sigma, double& eps, double& vrho, double& vsigma)
{
    const double Cx = 1.5 * std::cbrt(3.0 / (4.0 * kPi));
    const double rho13 = std::cbrt(rho);
    const double y = std::cbrt(4.0) * sigma / (rho * rho * rho13 * rho13);  // 2^(2/3) sigma / rho^(8/3)
    const double x = std::sqrt(y);
    const double ash = std::asinh(x);
    const double x_over = y > 1.0 ? 1.0 / std::sqrt(1.0 + 1.0 / y) : x / std::sqrt(1.0 + y);  // x / sqrt(1 + x^2)
    const double D = 1.0 + 6.0 * kB88_beta * x * ash;
    const double g = y / D;
    const double gp = (D - 3.0 * kB88_beta * x * (ash + x_over)) / (D * D);
    const double pre = rho13 / std::cbrt(2.0);  // 2^(-1/3) rho^(1/3)
    eps = pre * (-Cx - kB88_beta * g);
    vrho = pre * (4.0 / 3.0 * (-Cx - kB88_beta * g) + 8.0 / 3.0 * kB88_beta * y * gp);
    vsigma = -std::cbrt(2.0) * kB88_beta * gp / (rho * rho13);
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
void lypCorrelation(double rho, double sigma, double& eps, double& vrho, double& vsigma)
{
    const double CF = 0.3 * std::pow(3.0 * kPi * kPi, 2.0 / 3.0);
    const double r = 1.0 / std::cbrt(rho);
    const double den = 1.0 + kLYP_d * r;
    const double E = std::exp(-kLYP_c * r);
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

void evaluateKernel(XCKernel k, double rho, double sigma, double& eps, double& vrho, double& vsigma)
{
    vsigma = 0.0;
    switch (k) {
        case XCKernel::SlaterX:  slaterExchange(rho, eps, vrho); return;
        case XCKernel::VWN5_C:   vwnCorrelation(rho, VWNVariant::VWN5, eps, vrho); return;
        case XCKernel::VWNRPA_C: vwnCorrelation(rho, VWNVariant::VWN_RPA, eps, vrho); return;
        case XCKernel::PW92_C:   pw92Correlation(rho, eps, vrho); return;
        case XCKernel::PBE_X:    pbeExchange(rho, sigma, eps, vrho, vsigma); return;
        case XCKernel::PBE_C:    pbeCorrelation(rho, sigma, eps, vrho, vsigma); return;
        case XCKernel::B88_X:    b88Exchange(rho, sigma, eps, vrho, vsigma); return;
        case XCKernel::LYP_C:    lypCorrelation(rho, sigma, eps, vrho, vsigma); return;
    }
    eps = vrho = std::numeric_limits<double>::quiet_NaN();
}

void evaluateXC(const XCFunctional& f, const Vector& rho, const Vector& sigma, double rho_threshold, Vector& ex,
                Vector& ec, Vector& vrho, Vector& vsigma)
{
    const int n = static_cast<int>(rho.size());
    const bool gga = f.needsGradient();
    const double rho_min = effectiveThreshold(f, rho_threshold);
    ex = Vector::Zero(n);
    ec = Vector::Zero(n);
    vrho = Vector::Zero(n);
    vsigma = Vector::Zero(n);
    for (int p = 0; p < n; ++p) {
        const double r = rho(p);
        if (!(r >= rho_min) || r <= 0.0) continue;  // also drops NaN and negative round-off densities
        const double sg = gga ? std::max(0.0, sigma(p)) : 0.0;
        for (const XCTerm& t : f.terms) {
            if (t.weight == 0.0) continue;
            double e = 0.0, v = 0.0, vs = 0.0;
            evaluateKernel(t.kernel, r, sg, e, v, vs);
            if (isExchangeKernel(t.kernel))
                ex(p) += t.weight * e;
            else
                ec(p) += t.weight * e;
            vrho(p) += t.weight * v;
            if (gga) vsigma(p) += t.weight * vs;
        }
    }
}

void buildVxc(const qmgrid::MolecularGrid& grid, const qmgrid::AOEvaluator& ao, const Matrix& P,
              const XCFunctional& f, double rho_threshold, int threads, XCResult& out, int block)
{
    const auto t0 = std::chrono::steady_clock::now();
    const int n = ao.nbf();
    const int npts = grid.size();
    const bool gga = f.needsGradient();
    const double rho_min = effectiveThreshold(f, rho_threshold);
    if (block < 1) block = 128;
    const int nblocks = (npts + block - 1) / block;
    // Fixed chunk count from the basis size only (memory budget ~128 MB for the per-chunk matrices), so the
    // summation order, and with it every bit of the result, is the same for any thread count.
    const double mat_mb = 8.0 * double(n) * n / 1.0e6;
    const int nchunks = std::max(1, std::min({ 32, nblocks, static_cast<int>(128.0 / std::max(mat_mb, 1e-3)) }));

    std::vector<Eigen::MatrixXd> Vc(nchunks);
    std::vector<double> exc_c(nchunks, 0.0), ec_c(nchunks, 0.0), ne_c(nchunks, 0.0);
    std::vector<long> used_c(nchunks, 0);

#pragma omp parallel for schedule(dynamic, 1) num_threads(threads > 0 ? threads : 1)
    for (int c = 0; c < nchunks; ++c) {
        const int b0 = static_cast<int>((long long)nblocks * c / nchunks);
        const int b1 = static_cast<int>((long long)nblocks * (c + 1) / nchunks);
        Eigen::MatrixXd V = Eigen::MatrixXd::Zero(n, n);
        double ex_sum = 0.0, ec_sum = 0.0, ne_sum = 0.0;
        long used = 0;
        qmgrid::AOBlock blk;
        Matrix grad;  // grad rho on the block (GGA only)
        for (int b = b0; b < b1; ++b) {
            const int begin = b * block;
            const int count = std::min(block, npts - begin);
            ao.evaluate(grid.points, begin, count, gga ? 1 : 0, blk);
            Vector rho, sigma;
            qmgrid::densityOnBlock(blk, P, rho, gga ? &grad : nullptr);
            if (gga) sigma = grad.rowwise().squaredNorm();
            Vector ex, ec, vrho, vsigma;
            evaluateXC(f, rho, sigma, rho_threshold, ex, ec, vrho, vsigma);
            // weighted potential per point; points below the density threshold carry vrho = 0
            Eigen::VectorXd wv(count);
            bool any = false;
            for (int q = 0; q < count; ++q) {
                const double w = grid.weights(begin + q);
                wv(q) = w * vrho(q);
                ne_sum += w * rho(q);
                if (rho(q) >= rho_min && rho(q) > 0.0) {
                    ex_sum += w * rho(q) * ex(q);
                    ec_sum += w * rho(q) * ec(q);
                    ++used;
                    any = true;
                }
            }
            if (!any) continue;
            const Eigen::MatrixXd phi(blk.phi);  // count x n, column-major
            if (!gga) {
                // LDA: V += phi^T diag(w v) phi  (one GEMM per block)
                const Eigen::MatrixXd wphi = wv.asDiagonal() * phi;
                V.noalias() += phi.transpose() * wphi;
                continue;
            }
            // GGA (WP6): Y = (w v_rho / 2) phi + 2 w v_sigma (grad rho . grad phi), V += phi^T Y + Y^T phi
            //            = sum_g w [v_rho phi_m phi_n + 2 v_sigma grad rho . (phi_m grad phi_n + grad phi_m phi_n)]
            Eigen::VectorXd gx(count), gy(count), gz(count);
            for (int q = 0; q < count; ++q) {
                const double w2s = 2.0 * grid.weights(begin + q) * vsigma(q);
                gx(q) = w2s * grad(q, 0);
                gy(q) = w2s * grad(q, 1);
                gz(q) = w2s * grad(q, 2);
            }
            Eigen::MatrixXd Y = (0.5 * wv).asDiagonal() * phi;
            Y.noalias() += gx.asDiagonal() * Eigen::MatrixXd(blk.dx);
            Y.noalias() += gy.asDiagonal() * Eigen::MatrixXd(blk.dy);
            Y.noalias() += gz.asDiagonal() * Eigen::MatrixXd(blk.dz);
            const Eigen::MatrixXd X = phi.transpose() * Y;
            V.noalias() += X;
            V.noalias() += X.transpose();
        }
        Vc[c] = std::move(V);
        exc_c[c] = ex_sum;
        ec_c[c] = ec_sum;
        ne_c[c] = ne_sum;
        used_c[c] = used;
    }

    Eigen::MatrixXd V = Eigen::MatrixXd::Zero(n, n);
    out.ex = out.ec = out.nelec = 0.0;
    out.points_used = 0;
    for (int c = 0; c < nchunks; ++c) {
        V += Vc[c];
        out.ex += exc_c[c];
        out.ec += ec_c[c];
        out.nelec += ne_c[c];
        out.points_used += used_c[c];
    }
    out.exc = out.ex + out.ec;
    out.V = Matrix(0.5 * (V + V.transpose()));
    out.ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}

// WP8 (Claude Generated, Oct 2026): nuclear gradient of E_xc, derivation in qm_xc.h (buildXCGradient).
void buildXCGradient(const qmgrid::MolecularGrid& grid, const qmgrid::AOEvaluator& ao, const Matrix& P,
                     const XCFunctional& f, double rho_threshold, int threads, const std::vector<int>& bf_atom,
                     int natoms, XCGradient& out, int block)
{
    const auto t0 = std::chrono::steady_clock::now();
    const int n = ao.nbf();
    const int npts = grid.size();
    const bool gga = f.needsGradient();
    const double rho_min = effectiveThreshold(f, rho_threshold);
    if (block < 1) block = 128;
    const int nblocks = (npts + block - 1) / block;
    const int nchunks = std::max(1, std::min(32, nblocks));
    const Eigen::MatrixXd Pd(P);
    std::vector<Matrix> gb(nchunks);
    std::vector<long> used_c(nchunks, 0);
    out.edens = Vector::Zero(npts);
    Matrix spoint = Matrix::Zero(npts, 3);  // per point: sum over all functions of the basis part

#pragma omp parallel for schedule(dynamic, 1) num_threads(threads > 0 ? threads : 1)
    for (int c = 0; c < nchunks; ++c) {
        const int b0 = static_cast<int>((long long)nblocks * c / nchunks);
        const int b1 = static_cast<int>((long long)nblocks * (c + 1) / nchunks);
        Matrix g = Matrix::Zero(natoms, 3);
        long used = 0;
        qmgrid::AOBlock blk;
        Matrix grad;
        for (int b = b0; b < b1; ++b) {
            const int begin = b * block;
            const int count = std::min(block, npts - begin);
            ao.evaluate(grid.points, begin, count, gga ? 2 : 1, blk);
            Vector rho, sigma;
            qmgrid::densityOnBlock(blk, P, rho, gga ? &grad : nullptr);
            if (gga) sigma = grad.rowwise().squaredNorm();
            Vector ex, ec, vrho, vsigma;
            evaluateXC(f, rho, sigma, rho_threshold, ex, ec, vrho, vsigma);
            Eigen::VectorXd wv(count), vx(count), vy(count), vz(count);
            bool any = false;
            for (int q = 0; q < count; ++q) {
                const double w = grid.weights(begin + q);
                wv(q) = w * vrho(q);
                vx(q) = vy(q) = vz(q) = 0.0;
                if (rho(q) >= rho_min && rho(q) > 0.0) {
                    out.edens(begin + q) = rho(q) * (ex(q) + ec(q));
                    ++used;
                    any = true;
                    if (gga) {
                        const double w2s = 2.0 * w * vsigma(q);
                        vx(q) = w2s * grad(q, 0);
                        vy(q) = w2s * grad(q, 1);
                        vz(q) = w2s * grad(q, 2);
                    }
                }
            }
            if (!any) continue;
            const Eigen::MatrixXd phi(blk.phi), dx(blk.dx), dy(blk.dy), dz(blk.dz);
            Eigen::MatrixXd Zm = wv.asDiagonal() * phi;
            if (gga) {
                Zm.noalias() += vx.asDiagonal() * dx;
                Zm.noalias() += vy.asDiagonal() * dy;
                Zm.noalias() += vz.asDiagonal() * dz;
            }
            const Eigen::MatrixXd ZP = Zm * Pd;  // (Z P)_pm, P symmetric
            Eigen::MatrixXd T[3] = { dx.cwiseProduct(ZP), dy.cwiseProduct(ZP), dz.cwiseProduct(ZP) };
            if (gga) {
                const Eigen::MatrixXd PhiP = phi * Pd;
                const Eigen::MatrixXd xx(blk.dxx), xy(blk.dxy), xz(blk.dxz), yy(blk.dyy), yz(blk.dyz), zz(blk.dzz);
                // sum_j v_j d_k d_j phi_m (P phi)_m
                T[0].noalias() += ((vx.asDiagonal() * xx) + (vy.asDiagonal() * xy) + (vz.asDiagonal() * xz)).cwiseProduct(PhiP);
                T[1].noalias() += ((vx.asDiagonal() * xy) + (vy.asDiagonal() * yy) + (vz.asDiagonal() * yz)).cwiseProduct(PhiP);
                T[2].noalias() += ((vx.asDiagonal() * xz) + (vy.asDiagonal() * yz) + (vz.asDiagonal() * zz)).cwiseProduct(PhiP);
            }
            for (int k = 0; k < 3; ++k) {
                const Eigen::VectorXd colsum = T[k].colwise().sum().transpose();  // per function
                for (int m = 0; m < n; ++m) g(bf_atom[m], k) += -2.0 * colsum(m);
                const Eigen::VectorXd rowsum = T[k].rowwise().sum();               // per point
                for (int q = 0; q < count; ++q) spoint(begin + q, k) = -2.0 * rowsum(q);
            }
        }
        gb[c] = std::move(g);
        used_c[c] = used;
    }
    out.basis = Matrix::Zero(natoms, 3);
    out.points_used = 0;
    for (int c = 0; c < nchunks; ++c) {
        out.basis += gb[c];
        out.points_used += used_c[c];
    }
    // points of atom B move with B: d e(r_p)/d R_B (point part) = d_r e = - (sum of all basis parts of the point)
    out.points = Matrix::Zero(natoms, 3);
    for (int p = 0; p < npts; ++p)
        for (int k = 0; k < 3; ++k) out.points(grid.atom[p], k) -= spoint(p, k);
    out.ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}

}  // namespace qmxc
