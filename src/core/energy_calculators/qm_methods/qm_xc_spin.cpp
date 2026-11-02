/*
 * <Native QM: spin-polarised exchange-correlation for open-shell Kohn-Sham (UKS), stage 5b>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Formulas and literature: qm_xc_spin.h. Every correlation kernel below was first written as a short Python
 * prototype and compared with libxc 7.0.0 (pyscf.dft.libxc.eval_xc, spin = 1) before this translation; the ctest
 * qm_uks_kernel repeats the comparison for eps, v_rho and v_sigma on a fixed point set.
 *
 * Derivatives: forward-mode dual numbers in the five variables (rho_a, rho_b, sigma_aa, sigma_ab, sigma_bb), so the
 * potential is the exact derivative of the energy formula written here (no hand-derived v formulas to get wrong).
 *
 * Claude Generated (Nov 2026, stage 5b). AI-generated, machine-tested only.
 *
 * This program is free software under GPL-3.0
 */

#include "qm_xc_spin.h"

#include <Eigen/Dense>
#include <fmt/format.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>

namespace qmxc {

namespace {

// ---------------------------------------------------------------------------------------------------------------
// Dual numbers in five variables: value v and partial derivatives d[k] = dv / dx_k,
// x = (rho_a, rho_b, sigma_aa, sigma_ab, sigma_bb).
// ---------------------------------------------------------------------------------------------------------------
struct D5 {
    double v = 0.0;
    double d[5] = { 0.0, 0.0, 0.0, 0.0, 0.0 };
};

D5 cst(double c)
{
    D5 r;
    r.v = c;
    return r;
}
D5 var(double x, int k)
{
    D5 r;
    r.v = x;
    r.d[k] = 1.0;
    return r;
}
/// f(a) with value fv and derivative fp = f'(a.v): chain rule.
D5 chain(const D5& a, double fv, double fp)
{
    D5 r;
    r.v = fv;
    for (int k = 0; k < 5; ++k) r.d[k] = fp * a.d[k];
    return r;
}
D5 operator+(const D5& a, const D5& b)
{
    D5 r;
    r.v = a.v + b.v;
    for (int k = 0; k < 5; ++k) r.d[k] = a.d[k] + b.d[k];
    return r;
}
D5 operator-(const D5& a, const D5& b)
{
    D5 r;
    r.v = a.v - b.v;
    for (int k = 0; k < 5; ++k) r.d[k] = a.d[k] - b.d[k];
    return r;
}
D5 operator-(const D5& a) { return chain(a, -a.v, -1.0); }
D5 operator*(const D5& a, const D5& b)
{
    D5 r;
    r.v = a.v * b.v;
    for (int k = 0; k < 5; ++k) r.d[k] = a.d[k] * b.v + a.v * b.d[k];
    return r;
}
D5 operator/(const D5& a, const D5& b)
{
    D5 r;
    r.v = a.v / b.v;
    const double ib2 = 1.0 / (b.v * b.v);
    for (int k = 0; k < 5; ++k) r.d[k] = (a.d[k] * b.v - a.v * b.d[k]) * ib2;
    return r;
}
D5 operator+(const D5& a, double c) { return chain(a, a.v + c, 1.0); }
D5 operator+(double c, const D5& a) { return chain(a, a.v + c, 1.0); }
D5 operator-(const D5& a, double c) { return chain(a, a.v - c, 1.0); }
D5 operator-(double c, const D5& a) { return chain(a, c - a.v, -1.0); }
D5 operator*(double c, const D5& a) { return chain(a, a.v * c, c); }
D5 operator/(const D5& a, double c) { return chain(a, a.v / c, 1.0 / c); }
D5 operator/(double c, const D5& a) { return chain(a, c / a.v, -c / (a.v * a.v)); }
D5 dexp(const D5& a)
{
    const double e = std::exp(a.v);
    return chain(a, e, e);
}
D5 dlog(const D5& a) { return chain(a, std::log(a.v), 1.0 / a.v); }
D5 dsqrt(const D5& a)
{
    const double s = std::sqrt(a.v);
    return chain(a, s, 0.5 / s);
}
D5 datan(const D5& a) { return chain(a, std::atan(a.v), 1.0 / (1.0 + a.v * a.v)); }
/// a^p for a > 0 (p real); a == 0 with p > 1 gives 0 and a zero derivative.
D5 dpow(const D5& a, double p)
{
    if (a.v == 0.0 && p > 1.0) return chain(a, 0.0, 0.0);
    const double vp = std::pow(a.v, p);
    return chain(a, vp, p * vp / a.v);
}

constexpr double kPi = 3.14159265358979323846;
/// libxc's zeta_threshold (DBL_EPSILON): 1 + zeta and 1 - zeta are clamped to it (libxc my_piecewise in f_zeta, phi).
constexpr double kZetaThreshold = 2.220446049250313e-16;
/// (1 +- zeta)^p with libxc's clamp: below the threshold the constant threshold^p (no derivative).
D5 opzPow(const D5& x, double p)
{
    if (x.v <= kZetaThreshold) return cst(std::pow(kZetaThreshold, p));
    return dpow(x, p);
}

/// Spin interpolation f(z) = ((1+z)^(4/3) + (1-z)^(4/3) - 2) / (2^(4/3) - 2) (VWN 1980 eq. 2.4; PW92 eq. 9).
D5 fZeta(const D5& z)
{
    return (opzPow(1.0 + z, 4.0 / 3.0) + opzPow(1.0 - z, 4.0 / 3.0) - 2.0) / (std::pow(2.0, 4.0 / 3.0) - 2.0);
}

/// VWN fit function of x = sqrt(r_s) (VWN 1980 eq. 4.4) for one parameter set.
D5 vwnAux(double A, double b, double c, double x0, const D5& x)
{
    const D5 X = x * x + b * x + c;
    const double X0 = x0 * x0 + b * x0 + c;
    const double Q = std::sqrt(4.0 * c - b * b);
    const D5 at = datan(Q / (2.0 * x + b));
    return A * (dlog(x * x / X) + (2.0 * b / Q) * at
                - (b * x0 / X0) * (dlog((x - x0) * (x - x0) / X) + (2.0 * (b + 2.0 * x0) / Q) * at));
}

/// VWN correlation energy per electron, spin-polarised. rpa = false: libxc LDA_C_VWN (VWN5 Monte-Carlo fits,
/// interpolation with the spin stiffness, VWN 1980 eq. 4.4/ "functional V"); rpa = true: libxc LDA_C_VWN_RPA (RPA fits,
/// plain f(z) interpolation between the para- and ferromagnetic fits; the closed-shell value is the vwn_rpa kernel).
D5 vwnSpin(const D5& rs, const D5& z, bool rpa)
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
    const double fpp = 4.0 / (9.0 * (std::cbrt(2.0) - 1.0));
    const D5 z4 = z * z * z * z;
    return eP + aC * fz * (1.0 - z4) / fpp + (eF - eP) * fz * z4;
}

/// PW92 G(r_s) of one fit (PW92 eq. 10).
D5 pw92GSpin(const D5& rs, double A, double a1, double b1, double b2, double b3, double b4)
{
    const D5 srs = dsqrt(rs);
    const D5 Q1 = 2.0 * A * (b1 * srs + b2 * rs + b3 * rs * srs + b4 * rs * rs);
    return -2.0 * A * (1.0 + a1 * rs) * dlog(1.0 + 1.0 / Q1);
}

/// PW92 correlation per electron, spin-polarised (PW92 eq. 8), libxc LDA_C_PW_MOD constants and f''(0):
///   eps = G1 + z^4 f(z) (G2 - G1 + G3/f''(0)) - f(z) G3/f''(0)    (G3 = -alpha_c, libxc's form)
D5 pw92Spin(const D5& rs, const D5& z)
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
D5 pbeCSpin(const D5& rho, const D5& rs, const D5& z, const D5& sigma)
{
    const double beta = 0.06672455060314922, gam = (1.0 - std::log(2.0)) / (kPi * kPi);
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
D5 lypSpinDensity(const D5& ra, const D5& rb, const D5& saa, const D5& sab, const D5& sbb)
{
    const double a = 0.04918, b = 0.132, c = 0.2533, d = 0.349;
    const double CF = 0.3 * std::pow(3.0 * kPi * kPi, 2.0 / 3.0);
    const D5 rho = ra + rb;
    const D5 r = dpow(rho, -1.0 / 3.0);
    const D5 den = 1.0 + d * r;
    const D5 w = dexp(-c * r) / den * dpow(rho, -11.0 / 3.0);
    const D5 delta = c * r + d * r / den;
    const D5 s = saa + 2.0 * sab + sbb;
    const D5 rho2 = rho * rho;
    D5 t = ra * rb
        * (std::pow(2.0, 11.0 / 3.0) * CF * (dpow(ra, 8.0 / 3.0) + dpow(rb, 8.0 / 3.0)) + (47.0 / 18.0 - 7.0 / 18.0 * delta) * s
           - (2.5 - delta / 18.0) * (saa + sbb) - (delta - 11.0) / 9.0 * (ra / rho * saa + rb / rho * sbb));
    t = t - (2.0 / 3.0) * rho2 * s + ((2.0 / 3.0) * rho2 - ra * ra) * sbb + ((2.0 / 3.0) * rho2 - rb * rb) * saa;
    return -4.0 * a / den * ra * rb / rho - a * b * w * t;
}

/// Correlation energy density of one kernel at (ra, rb, saa, sab, sbb) as a dual number; false: no spin form.
bool correlationSpin(XCKernel k, const D5& ra, const D5& rb, const D5& saa, const D5& sab, const D5& sbb, D5& e)
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

}  // namespace

XCResultSpin::XCResultSpin() = default;
XCResultSpin::XCResultSpin(const XCResultSpin&) = default;
XCResultSpin& XCResultSpin::operator=(const XCResultSpin&) = default;
XCResultSpin::~XCResultSpin() = default;

bool supportsSpinPolarised(const XCFunctional& f, std::string& why)
{
    for (const XCTerm& t : f.terms) {
        if (t.weight == 0.0) continue;
        switch (t.kernel) {
            case XCKernel::SlaterX:
            case XCKernel::PBE_X:
            case XCKernel::B88_X:
            case XCKernel::VWN5_C:
            case XCKernel::VWNRPA_C:
            case XCKernel::PW92_C:
            case XCKernel::PBE_C:
            case XCKernel::LYP_C:
                break;
            default:
                why = "no spin-polarised form of the kernel " + kernelName(t.kernel)
                    + " (open-shell Kohn-Sham: Slater, B88, PBE-X, VWN5, VWN-RPA, PW92, PBE-C, LYP, i.e. lda, pbe, b3lyp)";
                return false;
        }
    }
    return true;
}

bool evaluatePointSpin(const XCFunctional& f, bool gga, double rho_min, double ra, double rb, double saa, double sab,
                       double sbb, double& ex, double& ec, double v[5])
{
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
    for (const XCTerm& t : f.terms) {
        const double w = t.weight;
        if (w == 0.0) continue;
        if (kern::isExchange(t.kernel)) {
            // spin scaling: e = sum_s rho_s eps(2 rho_s, 4 sigma_ss) with the closed-shell kernel
            const double rs[2] = { ra, rb }, ss[2] = { saa, sbb };
            for (int s = 0; s < 2; ++s) {
                const double r2 = 2.0 * rs[s];
                if (!(r2 >= rho_min) || r2 <= 0.0) continue;
                double eps = 0.0, vr = 0.0, vs = 0.0;
                kern::evaluateKernel(t.kernel, r2, gga ? 4.0 * ss[s] : 0.0, eps, vr, vs);
                ex += w * rs[s] * eps;
                v[s] += w * vr;
                if (gga) v[s == 0 ? 2 : 4] += w * 2.0 * vs;
            }
            continue;
        }
        D5 e;
        if (!correlationSpin(t.kernel, Ra, Rb, Saa, Sab, Sbb, e)) {
            ex = ec = std::numeric_limits<double>::quiet_NaN();  // never a silent zero
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

void buildVxcSpin(const qmgrid::MolecularGrid& grid, const qmgrid::AOEvaluator& ao, const Matrix& Pa, const Matrix& Pb,
                  const XCFunctional& f, double rho_threshold, int threads, XCResultSpin& out, int block)
{
    const auto t0 = std::chrono::steady_clock::now();
    const int n = ao.nbf();
    const int npts = grid.size();
    const bool gga = f.needsGradient();
    const double rho_min = kern::effectiveThreshold(gga, rho_threshold);
    if (block < 1) block = 128;
    const int nblocks = (npts + block - 1) / block;
    // fixed chunk count (basis size only, two matrices per chunk): thread-count independent result
    const double mat_mb = 16.0 * double(n) * n / 1.0e6;
    const int nchunks = std::max(1, std::min({ 32, nblocks, static_cast<int>(128.0 / std::max(mat_mb, 1e-3)) }));
    std::vector<Eigen::MatrixXd> Vac(nchunks), Vbc(nchunks);
    std::vector<double> ex_c(nchunks, 0.0), ec_c(nchunks, 0.0), ne_c(nchunks, 0.0);
    std::vector<long> used_c(nchunks, 0);

#pragma omp parallel for schedule(dynamic, 1) num_threads(threads > 0 ? threads : 1)
    for (int c = 0; c < nchunks; ++c) {
        const int b0 = static_cast<int>((long long)nblocks * c / nchunks);
        const int b1 = static_cast<int>((long long)nblocks * (c + 1) / nchunks);
        Eigen::MatrixXd Va = Eigen::MatrixXd::Zero(n, n), Vb = Eigen::MatrixXd::Zero(n, n);
        double ex_sum = 0.0, ec_sum = 0.0, ne_sum = 0.0;
        long used = 0;
        qmgrid::AOBlock blk;
        Matrix ga, gb;
        for (int b = b0; b < b1; ++b) {
            const int begin = b * block;
            const int count = std::min(block, npts - begin);
            ao.evaluate(grid.points, begin, count, gga ? 1 : 0, blk);
            Vector ra, rb;
            qmgrid::densityOnBlock(blk, Pa, ra, gga ? &ga : nullptr);
            qmgrid::densityOnBlock(blk, Pb, rb, gga ? &gb : nullptr);
            Eigen::VectorXd wva = Eigen::VectorXd::Zero(count), wvb = Eigen::VectorXd::Zero(count);
            Eigen::MatrixXd ua = Eigen::MatrixXd::Zero(count, 3), ub = Eigen::MatrixXd::Zero(count, 3);
            bool any = false;
            for (int q = 0; q < count; ++q) {
                const double w = grid.weights(begin + q);
                ne_sum += w * (ra(q) + rb(q));
                double saa = 0.0, sab = 0.0, sbb = 0.0;
                if (gga) {
                    saa = ga.row(q).squaredNorm();
                    sbb = gb.row(q).squaredNorm();
                    sab = ga.row(q).dot(gb.row(q));
                }
                double ex = 0.0, ec = 0.0, v[5];
                if (!evaluatePointSpin(f, gga, rho_min, ra(q), rb(q), saa, sab, sbb, ex, ec, v)) continue;
                ex_sum += w * ex;
                ec_sum += w * ec;
                ++used;
                any = true;
                wva(q) = w * v[0];
                wvb(q) = w * v[1];
                if (gga)
                    for (int k = 0; k < 3; ++k) {
                        ua(q, k) = w * (2.0 * v[2] * ga(q, k) + v[3] * gb(q, k));
                        ub(q, k) = w * (2.0 * v[4] * gb(q, k) + v[3] * ga(q, k));
                    }
            }
            if (!any) continue;
            const Eigen::MatrixXd phi(blk.phi);
            if (!gga) {
                Va.noalias() += phi.transpose() * (wva.asDiagonal() * phi);
                Vb.noalias() += phi.transpose() * (wvb.asDiagonal() * phi);
                continue;
            }
            const Eigen::MatrixXd dx(blk.dx), dy(blk.dy), dz(blk.dz);
            // Y_s = (w v_rho_s / 2) phi + u_s . grad phi;  V_s += phi^T Y_s + Y_s^T phi
            Eigen::MatrixXd Ya = (0.5 * wva).asDiagonal() * phi;
            Ya.noalias() += ua.col(0).asDiagonal() * dx;
            Ya.noalias() += ua.col(1).asDiagonal() * dy;
            Ya.noalias() += ua.col(2).asDiagonal() * dz;
            Eigen::MatrixXd Yb = (0.5 * wvb).asDiagonal() * phi;
            Yb.noalias() += ub.col(0).asDiagonal() * dx;
            Yb.noalias() += ub.col(1).asDiagonal() * dy;
            Yb.noalias() += ub.col(2).asDiagonal() * dz;
            const Eigen::MatrixXd Xa = phi.transpose() * Ya, Xb = phi.transpose() * Yb;
            Va.noalias() += Xa;
            Va.noalias() += Xa.transpose();
            Vb.noalias() += Xb;
            Vb.noalias() += Xb.transpose();
        }
        Vac[c] = std::move(Va);
        Vbc[c] = std::move(Vb);
        ex_c[c] = ex_sum;
        ec_c[c] = ec_sum;
        ne_c[c] = ne_sum;
        used_c[c] = used;
    }
    Eigen::MatrixXd Va = Eigen::MatrixXd::Zero(n, n), Vb = Eigen::MatrixXd::Zero(n, n);
    out.ex = out.ec = out.nelec = 0.0;
    out.points_used = 0;
    for (int c = 0; c < nchunks; ++c) {  // chunk order: fixed summation order
        Va += Vac[c];
        Vb += Vbc[c];
        out.ex += ex_c[c];
        out.ec += ec_c[c];
        out.nelec += ne_c[c];
        out.points_used += used_c[c];
    }
    out.exc = out.ex + out.ec;
    out.Va = Matrix(0.5 * (Va + Va.transpose()));
    out.Vb = Matrix(0.5 * (Vb + Vb.transpose()));
    out.ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}

void buildXCGradientSpin(const qmgrid::MolecularGrid& grid, const qmgrid::AOEvaluator& ao, const Matrix& Pa,
                         const Matrix& Pb, const XCFunctional& f, double rho_threshold, int threads,
                         const std::vector<int>& bf_atom, int natoms, XCGradient& out, int block)
{
    const auto t0 = std::chrono::steady_clock::now();
    const int n = ao.nbf();
    const int npts = grid.size();
    const bool gga = f.needsGradient();
    const double rho_min = kern::effectiveThreshold(gga, rho_threshold);
    if (block < 1) block = 128;
    const int nblocks = (npts + block - 1) / block;
    const int nchunks = std::max(1, std::min(32, nblocks));
    const Eigen::MatrixXd Pad(Pa), Pbd(Pb);
    std::vector<Matrix> gbc(nchunks);
    std::vector<long> used_c(nchunks, 0);
    out.edens = Vector::Zero(npts);
    Matrix spoint = Matrix::Zero(npts, 3);

#pragma omp parallel for schedule(dynamic, 1) num_threads(threads > 0 ? threads : 1)
    for (int c = 0; c < nchunks; ++c) {
        const int b0 = static_cast<int>((long long)nblocks * c / nchunks);
        const int b1 = static_cast<int>((long long)nblocks * (c + 1) / nchunks);
        Matrix g = Matrix::Zero(natoms, 3);
        long used = 0;
        qmgrid::AOBlock blk;
        Matrix ga, gb;
        for (int b = b0; b < b1; ++b) {
            const int begin = b * block;
            const int count = std::min(block, npts - begin);
            ao.evaluate(grid.points, begin, count, gga ? 2 : 1, blk);
            Vector ra, rb;
            qmgrid::densityOnBlock(blk, Pa, ra, gga ? &ga : nullptr);
            qmgrid::densityOnBlock(blk, Pb, rb, gga ? &gb : nullptr);
            Eigen::VectorXd wv[2] = { Eigen::VectorXd::Zero(count), Eigen::VectorXd::Zero(count) };
            Eigen::MatrixXd u[2] = { Eigen::MatrixXd::Zero(count, 3), Eigen::MatrixXd::Zero(count, 3) };
            bool any = false;
            for (int q = 0; q < count; ++q) {
                const double w = grid.weights(begin + q);
                double saa = 0.0, sab = 0.0, sbb = 0.0;
                if (gga) {
                    saa = ga.row(q).squaredNorm();
                    sbb = gb.row(q).squaredNorm();
                    sab = ga.row(q).dot(gb.row(q));
                }
                double ex = 0.0, ec = 0.0, v[5];
                if (!evaluatePointSpin(f, gga, rho_min, ra(q), rb(q), saa, sab, sbb, ex, ec, v)) continue;
                out.edens(begin + q) = ex + ec;
                ++used;
                any = true;
                wv[0](q) = w * v[0];
                wv[1](q) = w * v[1];
                if (gga)
                    for (int k = 0; k < 3; ++k) {
                        u[0](q, k) = w * (2.0 * v[2] * ga(q, k) + v[3] * gb(q, k));
                        u[1](q, k) = w * (2.0 * v[4] * gb(q, k) + v[3] * ga(q, k));
                    }
            }
            if (!any) continue;
            const Eigen::MatrixXd phi(blk.phi), dx(blk.dx), dy(blk.dy), dz(blk.dz);
            Eigen::MatrixXd xx, xy, xz, yy, yz, zz;
            if (gga) {
                xx = blk.dxx;
                xy = blk.dxy;
                xz = blk.dxz;
                yy = blk.dyy;
                yz = blk.dyz;
                zz = blk.dzz;
            }
            Eigen::MatrixXd rowsum_all = Eigen::MatrixXd::Zero(count, 3);
            for (int s = 0; s < 2; ++s) {
                const Eigen::MatrixXd& Ps = s == 0 ? Pad : Pbd;
                Eigen::MatrixXd Zm = wv[s].asDiagonal() * phi;
                if (gga) {
                    Zm.noalias() += u[s].col(0).asDiagonal() * dx;
                    Zm.noalias() += u[s].col(1).asDiagonal() * dy;
                    Zm.noalias() += u[s].col(2).asDiagonal() * dz;
                }
                const Eigen::MatrixXd ZP = Zm * Ps;
                Eigen::MatrixXd T[3] = { dx.cwiseProduct(ZP), dy.cwiseProduct(ZP), dz.cwiseProduct(ZP) };
                if (gga) {
                    const Eigen::MatrixXd PhiP = phi * Ps;
                    const auto ux = u[s].col(0).asDiagonal(), uy = u[s].col(1).asDiagonal(), uz = u[s].col(2).asDiagonal();
                    T[0].noalias() += ((ux * xx) + (uy * xy) + (uz * xz)).cwiseProduct(PhiP);
                    T[1].noalias() += ((ux * xy) + (uy * yy) + (uz * yz)).cwiseProduct(PhiP);
                    T[2].noalias() += ((ux * xz) + (uy * yz) + (uz * zz)).cwiseProduct(PhiP);
                }
                for (int k = 0; k < 3; ++k) {
                    const Eigen::VectorXd colsum = T[k].colwise().sum().transpose();
                    for (int m = 0; m < n; ++m) g(bf_atom[m], k) += -2.0 * colsum(m);
                    rowsum_all.col(k) += T[k].rowwise().sum();
                }
            }
            for (int q = 0; q < count; ++q)
                for (int k = 0; k < 3; ++k) spoint(begin + q, k) = -2.0 * rowsum_all(q, k);
        }
        gbc[c] = std::move(g);
        used_c[c] = used;
    }
    out.basis = Matrix::Zero(natoms, 3);
    out.points_used = 0;
    for (int c = 0; c < nchunks; ++c) {
        out.basis += gbc[c];
        out.points_used += used_c[c];
    }
    out.points = Matrix::Zero(natoms, 3);
    for (int p = 0; p < npts; ++p)
        for (int k = 0; k < 3; ++k) out.points(grid.atom[p], k) -= spoint(p, k);
    out.ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}

}  // namespace qmxc
