/*
 * <Geometric counterpoise (gCP) + short-range basis (SRB) correction>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Ported from dftd3/simple-dftd3 src/dftd3/gcp.f90 (gcp_energy, gcp_deriv,
 * srb_energy, srb_deriv, dsovl, aaux_all, baux_all, bint_all) and
 * src/dftd3/gcp/param.f90 (emiss_hf_minix, nbas_minix, slater_s/p, the
 * p_minix_bas branch), LGPL-3.0-or-later. D3 pair radii from
 * src/dftd3/data/vdwrad.f90. See gcp.h for the model.
 *
 * Claude Generated: native gCP/SRB for HF-3c (Sep 2026)
 *
 * This program is free software under GPL-3.0
 */

#include "gcp.h"

#include "src/core/units.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace gcp {
namespace {

constexpr int kMaxZ = 10;  // H-Ne, the range of the shipped MINIX basis

// BSSE per atom for HF/MINIX: param.f90 emiss_hf_minix =
// [emiss_hf_minis(1:2), 0.177871, 0.171596, emiss_hf_minis(5:10), ...].
// NOTE: the digits matter -- a table rounded to 5 places (the earlier Python
// witness) is off by up to 2.6e-7 Eh on NH3.
constexpr double kEmissMinix[kMaxZ + 1] = {
    0.0,
    0.042400, 0.028324,
    0.177871, 0.171596, 0.224237, 0.279950, 0.357906, 0.479012, 0.638518, 0.832349
};

// Number of basis functions per atom in MINIX (param.f90 nbas_minix).
constexpr double kNbasMinix[kMaxZ + 1] = { 0, 1, 1, 5, 5, 5, 5, 5, 5, 5, 5 };

// Slater exponents: s for H/He, the s/p average for Li-Ne (param.f90 slater_exp).
// FLOAT on purpose: the reference declares `[real(wp) :: 1.2000, ...]` with no
// `_wp` kind suffix, so Fortran reads each literal as default (single-precision)
// REAL and only then widens it -- e.g. 2.5644 becomes 2.5643999576568604. Using
// the exact decimals instead shifts E_gCP by ~2e-8 relative (6e-10 Eh on F2).
// The s/p average is formed in double, like `(slater_s + slater_p)/2` in wp.
constexpr float kSlaterS[kMaxZ + 1] = {
    0.0f, 1.2000f, 1.6469f, 0.6534f, 1.0365f, 1.3990f, 1.7210f, 2.0348f, 2.2399f, 2.5644f, 2.8812f
};
constexpr float kSlaterP[kMaxZ + 1] = {
    0.0f, 0.0000f, 0.0000f, 0.5305f, 0.8994f, 1.2685f, 1.6105f, 1.9398f, 2.0477f, 2.4022f, 2.7421f
};

// D3 van-der-Waals pair radii R0_AB in Angstrom, packed lower triangle
// (vdwrad.f90: index num1 + num2*(num2-1)/2 with num1 <= num2, 1-based).
constexpr double kR0abAngstrom[kMaxZ * (kMaxZ + 1) / 2] = {
    2.1823, 1.8547, 1.7347, 2.9086, 2.5732, 3.4956, 2.3550, 2.5095, 2.9802, 3.0982,
    2.5141, 2.3917, 2.9977, 2.9484, 3.2160, 2.4492, 2.2527, 3.1933, 3.0214, 2.9531,
    2.9103, 2.3667, 2.1328, 2.8784, 2.7660, 2.7776, 2.7063, 2.6225, 2.1768, 2.0625,
    2.6395, 2.6648, 2.6482, 2.5697, 2.4846, 2.4817, 2.0646, 1.9891, 2.5086, 2.6908,
    2.6233, 2.4770, 2.3885, 2.3511, 2.2996, 1.9892, 1.9251, 2.4190, 2.5473, 2.4994,
    2.4091, 2.3176, 2.2571, 2.1946, 2.1374
};

double slaterExponent(int Z)
{
    return (Z <= 2) ? static_cast<double>(kSlaterS[Z])
                  : 0.5 * (static_cast<double>(kSlaterS[Z]) + static_cast<double>(kSlaterP[Z]));
}

double r0abBohr(int za, int zb)
{
    const int lo = std::min(za, zb), hi = std::max(za, zb);
    return CurcumaUnit::Length::angstrom_to_bohr(kR0abAngstrom[lo + hi * (hi - 1) / 2 - 1]);
}

// Principal quantum number of the valence s shell used in the overlap (1 for
// H/He, 2 for Li-Ne): the reference `shell` data.
int shellOf(int Z) { return (Z <= 2) ? 1 : 2; }

// A_k(x) = int_1^inf t^k e^{-xt} dt, upward recursion (aaux_all).
void aauxAll(double x, double* a, int kmax)
{
    const double ex = std::exp(-x), rx = 1.0 / x;
    a[0] = ex * rx;
    for (int k = 1; k <= kmax; ++k)
        a[k] = (k * a[k - 1] + ex) * rx;
}

// B_k(x) = int_-1^1 t^k e^{-xt} dt, closed form (baux_all).
void bauxAll(double x, double* b, int kmax)
{
    const double ep = std::exp(x), em = std::exp(-x), rx = 1.0 / x;
    for (int k = 0; k <= kmax; ++k) {
        double term = rx;
        double sgn = (k % 2 == 0) ? 1.0 : -1.0;
        double sp = sgn * term, sm = term;
        for (int j = 1; j <= k; ++j) {
            term *= (k - j + 1) * rx;
            sgn = -sgn;
            sp += sgn * term;
            sm += term;
        }
        b[k] = ep * sp - em * sm;
    }
}

// B_k(x) as the truncated 12-term Taylor series (bint_all). The reference uses it
// whenever the two Slater exponents differ by less than 0.1, where the closed
// form cancels badly -- this includes like-element pairs (x == 0) and, in H-Ne,
// the He-C pair. Ported with the same truncation so the energy matches exactly.
void bintAll(double x, double* b, int kmax)
{
    if (std::abs(x) < 1.0e-6) {
        for (int k = 0; k <= kmax; ++k)
            b[k] = (k % 2 == 0) ? 2.0 / (k + 1.0) : 0.0;
        return;
    }
    constexpr int nterm = 12;
    double pw[nterm + 1];
    pw[0] = 1.0;
    for (int i = 1; i <= nterm; ++i)
        pw[i] = pw[i - 1] * (-x) / i;
    for (int k = 0; k <= kmax; ++k) {
        double acc = 0.0;
        for (int i = k % 2; i <= nterm; i += 2)
            acc += pw[i] / (k + i + 1.0);
        b[k] = 2.0 * acc;
    }
}

// Overlap S(R) of two normalized ns Slater functions and dS/dR (dsovl), in the
// elliptic-coordinate form S = cnorm R^m sum_t w_t A_{p_t}(ax) B_{q_t}(bx).
// Za/Zb are the element numbers (they select 1s or 2s), za/zb the exponents.
void slaterOverlap(double r, int Za, int Zb, double za, double zb, double& s0, double& s1)
{
    const bool lsame = std::abs(za - zb) < 0.1;
    const int ii = shellOf(Za) * shellOf(Zb);

    int m = 0, nterm = 0;
    double wt[4] = {}, cnorm = 0.0;
    int pa[4] = {}, qb[4] = {};
    switch (ii) {
    case 1:  // <1s|1s>
        m = 3; nterm = 2;
        wt[0] = 1.0; wt[1] = -1.0;
        pa[0] = 2; pa[1] = 0;
        qb[0] = 0; qb[1] = 2;
        cnorm = 0.25 * std::sqrt(std::pow(za * zb, 3));
        break;
    case 2:  // <1s|2s>: the 1s function is always "a"
        if (shellOf(Za) >= shellOf(Zb)) std::swap(za, zb);
        m = 4; nterm = 4;
        wt[0] = 1.0; wt[1] = -1.0; wt[2] = 1.0; wt[3] = -1.0;
        pa[0] = 3; pa[1] = 0; pa[2] = 2; pa[3] = 1;
        qb[0] = 0; qb[1] = 3; qb[2] = 1; qb[3] = 2;
        cnorm = std::sqrt(1.0 / 3.0) * std::sqrt(std::pow(za, 3) * std::pow(zb, 5)) * 0.125;
        break;
    case 4:  // <2s|2s>
        m = 5; nterm = 3;
        wt[0] = 1.0; wt[1] = 1.0; wt[2] = -2.0;
        pa[0] = 4; pa[1] = 0; pa[2] = 2;
        qb[0] = 0; qb[1] = 4; qb[2] = 2;
        cnorm = std::sqrt(std::pow(za * zb, 5)) * 0.0625 / 3.0;
        break;
    default:
        throw std::runtime_error("gCP: shell combination outside the H-Ne scope");
    }

    const double ha = 0.5 * (za + zb), hb = 0.5 * (zb - za);
    double av[9], bv[9];
    aauxAll(ha * r, av, 8);
    if (lsame)
        bintAll(hb * r, bv, 8);
    else
        bauxAll(hb * r, bv, 8);

    double f0 = 0.0, f1 = 0.0;
    for (int t = 0; t < nterm; ++t) {
        const int p = pa[t], q = qb[t];
        f0 += wt[t] * av[p] * bv[q];
        f1 -= wt[t] * (ha * av[p + 1] * bv[q] + hb * av[p] * bv[q + 1]);
    }
    const double rm = std::pow(r, m);
    s0 = cnorm * rm * f0;
    s1 = cnorm * (m * std::pow(r, m - 1) * f0 + rm * f1);
}

// Number of virtual MINIX functions of atom Z, as 1/sqrt(nbas - nel/2) (0 if
// the atom has none, e.g. He with a single 1s function).
double inverseSqrtVirtuals(int Z)
{
    const double xv = kNbasMinix[Z] - 0.5 * Z;
    return (xv >= 0.5) ? 1.0 / std::sqrt(xv) : 0.0;
}

double srbPair(const std::vector<int>& atoms, const Matrix& xyz, const Parameters& p,
               int i, int j, Matrix* gradient)
{
    const double dx = xyz(i, 0) - xyz(j, 0), dy = xyz(i, 1) - xyz(j, 1), dz = xyz(i, 2) - xyz(j, 2);
    const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
    if (r > p.cutoff || r < 1.0e-14) return 0.0;
    const double r0 = p.rscal * std::pow(r0abBohr(atoms[i], atoms[j]), 0.75);
    const double ff = -std::pow(static_cast<double>(atoms[i]) * atoms[j], 1.5);
    const double expt = std::exp(-r0 * r);
    if (gradient) {
        const double g = -p.qscal * ff * r0 * expt / r;  // (dE/dR) / R
        (*gradient)(i, 0) += g * dx; (*gradient)(j, 0) -= g * dx;
        (*gradient)(i, 1) += g * dy; (*gradient)(j, 1) -= g * dy;
        (*gradient)(i, 2) += g * dz; (*gradient)(j, 2) -= g * dz;
    }
    return p.qscal * ff * expt;
}

} // namespace

Parameters hf3c()
{
    // gcp/param.f90, case(p_minix_bas) with method hf3c: HF/MINIX gCP parameters
    // and base = .true. with rscal = 0.7, qscal = 0.03 (Sure & Grimme 2013).
    Parameters p;
    p.sigma = 0.1290;
    p.eta = 1.1526;
    p.alpha = 1.1549;
    p.beta = 1.1763;
    p.base = true;
    p.rscal = 0.7;
    p.qscal = 0.03;
    return p;
}

bool supports(const std::vector<int>& atoms)
{
    for (int Z : atoms)
        if (Z < 1 || Z > kMaxZ) return false;
    return true;
}

double baseEnergy(const std::vector<int>& atoms, const Matrix& xyz, const Parameters& p)
{
    if (!p.base) return 0.0;
    double e = 0.0;
    const int n = static_cast<int>(atoms.size());
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < i; ++j)
            e += srbPair(atoms, xyz, p, i, j, nullptr);
    return e;
}

double energy(const std::vector<int>& atoms, const Matrix& xyz, const Parameters& p,
              Matrix* gradient)
{
    if (!supports(atoms))
        throw std::runtime_error("gCP: element outside the H-Ne scope of the MINIX parameters");

    const int n = static_cast<int>(atoms.size());
    if (gradient) *gradient = Matrix::Zero(n, 3);

    double e = 0.0;
    for (int i = 0; i < n; ++i) {
        const int Zi = atoms[i];
        const double xvi = inverseSqrtVirtuals(Zi);
        for (int j = 0; j < i; ++j) {
            const int Zj = atoms[j];
            const double xvj = inverseSqrtVirtuals(Zj);
            // BSSE that atom i suffers from j's virtuals, and vice versa.
            const double emij = p.sigma * (kEmissMinix[Zi] * xvj + kEmissMinix[Zj] * xvi);

            const double dx = xyz(i, 0) - xyz(j, 0), dy = xyz(i, 1) - xyz(j, 1), dz = xyz(i, 2) - xyz(j, 2);
            const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
            if (r > p.cutoff || r < 1.0e-14) continue;

            double s = 0.0, ds = 0.0;
            slaterOverlap(r, Zi, Zj, p.eta * slaterExponent(Zi), p.eta * slaterExponent(Zj), s, ds);
            const double expv = std::exp(-p.alpha * std::pow(r, p.beta));
            const double bsse = expv / std::sqrt(s);
            e += emij * bsse;

            if (gradient) {
                // d/dR [exp(-a R^b) S^{-1/2}] = bsse * (-a b R^{b-1} - S'/(2S)); times vec/R.
                const double dEdR = emij * bsse * (-p.alpha * p.beta * std::pow(r, p.beta - 1.0) - 0.5 * ds / s);
                const double g = dEdR / r;
                (*gradient)(i, 0) += g * dx; (*gradient)(j, 0) -= g * dx;
                (*gradient)(i, 1) += g * dy; (*gradient)(j, 1) -= g * dy;
                (*gradient)(i, 2) += g * dz; (*gradient)(j, 2) -= g * dz;
            }
        }
    }

    if (p.base)
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < i; ++j)
                e += srbPair(atoms, xyz, p, i, j, gradient);
    return e;
}

} // namespace gcp
