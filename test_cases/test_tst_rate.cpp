/*
 * <Unit test for the absolute transition-state-theory rate across s* (fe_profile.h tstRate).>
 * Copyright (C) 2019 - 2026 Conrad Huebler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated 2026 (AI-generated, machine-tested).
 *
 * (1) 1D classical TST: a particle of mass m (amu) in V(x) = Vb (x^2 - 1)^2 (Eh, x in Angstrom),
 *     s = x, grad s = 1, Z = 1/m, <Z^1/2> = m^-1/2. Exact one-dimensional TST rate
 *         k = sqrt(kT / (2 pi m)) exp(-V(0)/kT) / int_{x<=0} exp(-V/kT) dx ,
 *     evaluated here completely in SI (kT in J, m in kg, x in m -> 1/s -> 1/ps) with a
 *     composite Simpson rule on 2e6 intervals, i.e. independently of the MD time unit used by
 *     tstRate. tstRate on a 4001-point grid (G = V) must agree to rel. 1e-3, both directions
 *     (the well is symmetric, plus a tilted variant for k_fwd != k_rev).
 * (2) Eyring conversion: k = 1/ps at T = 300 K.
 *         kB T / h = 1.380649e-23 * 300 / 6.62607015e-34 = 6.250986e12 1/s = 6.250986 1/ps
 *         dG = kT ln(6.250986) = 0.5961613 kcal/mol * 1.8327392 = 1.0926081 kcal/mol
 *     (kT = R T = 8.314462618 / 4184 * 300 kcal/mol). Required: |dG - 1.0926081| < 1e-6.
 * (3) Scale invariance: the same physics on s' = a s (a = 0.37, 3.1): G'(s') = G(s'/a),
 *     <Z'^1/2> = a <Z^1/2>. The rate must be unchanged (rel. 1e-3), the basin barrier shifts by
 *     kT ln a (both printed).
 * (4) Z-mean identity: n = 10000 samples Z_i log-normal, w_i = Z_i^-1/2:
 *     N / sum w == sum w_i Z_i^1/2 / sum w_i to machine precision (rel. 1e-12).
 * (5) zHalfNear: window mean (sample-weighted) and the interpolating fallback.
 */

#include <cmath>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "src/capabilities/fe_profile.h"

using namespace curcuma::fe_profile;

static int failures = 0;
static void check(bool ok, const std::string& what)
{
    std::cout << (ok ? "  PASS  " : "  FAIL  ") << what << std::endl;
    if (!ok)
        ++failures;
}

// SI constants (CODATA-2018) used only for the independent reference
static const double kB_SI = 1.380649e-23, Eh_J = 4.3597447222071e-18, amu_kg = 1.66053906660e-27;

/*! Exact 1D TST rate in 1/ps, SI route. V in Eh as a function of x in Angstrom. */
static double exactRate1D(double m_amu, double T, double (*V)(double), double xlo, double xstar, bool forward,
                          double xhi)
{
    const double kT_J = kB_SI * T;
    const double kT_Eh = kT_J / Eh_J;
    const int N = 2000000; // even
    const double a = forward ? xlo : xstar, b = forward ? xstar : xhi;
    const double h = (b - a) / N;
    double sum = 0.0;
    for (int i = 0; i <= N; ++i) {
        const double x = a + i * h;
        const double w = (i == 0 || i == N) ? 1.0 : (i % 2 ? 4.0 : 2.0);
        sum += w * std::exp(-V(x) / kT_Eh);
    }
    const double I_ang = sum * h / 3.0;          // Angstrom
    const double I_m = I_ang * 1e-10;            // m
    const double v = std::sqrt(kT_J / (2.0 * M_PI * m_amu * amu_kg)); // m/s
    const double k_s = v * std::exp(-V(xstar) / kT_Eh) / I_m;         // 1/s
    return k_s * 1e-12;                                                // 1/ps
}

static double g_Vb = 0.0, g_tilt = 0.0;
static double Vdw(double x) { return g_Vb * (x * x - 1.0) * (x * x - 1.0) + g_tilt * x; }

int main()
{
    std::cout << "=== test_tst_rate (Claude Generated 2026) ===" << std::endl;
    const double T = 300.0;
    const double kT = CurcumaUnit::Constants::BOLTZMANN_HARTREE_K * T;
    const double kcal = 627.5094740631;

    // ---------- (1) 1D analytic ----------
    std::cout << "(1) 1D double well, s = x" << std::endl;
    struct Case { double m, Vb_kcal, tilt_kcal; };
    for (const Case c : { Case{ 15.035, 3.0, 0.0 }, Case{ 1.008, 5.0, 0.0 }, Case{ 30.07, 4.0, 0.8 } }) {
        g_Vb = c.Vb_kcal / kcal;
        g_tilt = c.tilt_kcal / kcal;
        const double xlo = -2.0, xhi = 2.0;
        const int n = 4001;
        std::vector<double> s(n), G(n);
        for (int i = 0; i < n; ++i) {
            s[i] = xlo + (xhi - xlo) * i / (n - 1.0);
            G[i] = Vdw(s[i]);
        }
        // barrier top: maximum of V between the wells (dV/dx = 0 near 0), found by a fine scan
        double xs = 0.0, vmax = -1e300;
        for (int i = 0; i <= 100000; ++i) { // x in [-0.5, 0.5]: between the wells only
            const double x = -0.5 + i * 1e-5;
            if (Vdw(x) > vmax) {
                vmax = Vdw(x);
                xs = x;
            }
        }
        const TstRate r = tstRate(s, G, kT, xs, 1.0 / std::sqrt(c.m));
        const double kf = exactRate1D(c.m, T, Vdw, xlo, xs, true, xhi);
        const double kr = exactRate1D(c.m, T, Vdw, xlo, xs, false, xhi);
        const double ef = std::abs(r.k_fwd / kf - 1.0), er = std::abs(r.k_rev / kr - 1.0);
        std::cout << "    m = " << c.m << " amu, Vb = " << c.Vb_kcal << " kcal/mol, tilt " << c.tilt_kcal
                  << ", x* = " << xs << ": k_fwd " << r.k_fwd << " vs exact " << kf << " (rel " << ef
                  << "), k_rev " << r.k_rev << " vs " << kr << " (rel " << er << ") 1/ps, dG_TST fwd "
                  << r.dG_fwd * kcal << " kcal/mol" << std::endl;
        check(r.valid && ef < 1e-3 && er < 1e-3, "1D TST rate fwd/rev within rel. 1e-3 of the exact SI value");
    }

    // ---------- (2) Eyring ----------
    std::cout << "(2) Eyring conversion" << std::endl;
    {
        const double dG = eyringFreeEnergy(1.0, kT) * kcal;
        const double pre = eyringPrefactorPerPs(kT);
        std::cout << "    kB T/h = " << pre << " 1/ps, dG(k = 1/ps, 300 K) = " << dG << " kcal/mol (expected 1.0926081)"
                  << std::endl;
        check(std::abs(pre - 6.250986) < 1e-5, "kB T / h at 300 K = 6.250986 1/ps");
        check(std::abs(dG - 1.0926081) < 1e-6, "Eyring dG(k = 1/ps, 300 K) = 1.0926081 kcal/mol");
    }

    // ---------- (3) scale invariance ----------
    std::cout << "(3) scale invariance s' = a s" << std::endl;
    {
        g_Vb = 3.0 / kcal;
        g_tilt = 0.5 / kcal;
        const int n = 4001;
        std::vector<double> s(n), G(n);
        for (int i = 0; i < n; ++i) {
            s[i] = -2.0 + 4.0 * i / (n - 1.0);
            G[i] = Vdw(s[i]);
        }
        const double m = 15.0, zh = 1.0 / std::sqrt(m);
        const BasinBarrier b0 = basinBarrier(s, G, kT);
        const TstRate r0 = tstRate(s, G, kT, b0.s_star, zh);
        for (double a : { 0.37, 3.1 }) {
            std::vector<double> sa(n);
            for (int i = 0; i < n; ++i)
                sa[i] = a * s[i];
            const BasinBarrier ba = basinBarrier(sa, G, kT);
            const TstRate ra = tstRate(sa, G, kT, ba.s_star, a * zh);
            const double e = std::abs(ra.k_fwd / r0.k_fwd - 1.0) + std::abs(ra.k_rev / r0.k_rev - 1.0);
            std::cout << "    a = " << a << ": k_fwd " << ra.k_fwd << " vs " << r0.k_fwd << " 1/ps, dG_TST fwd "
                      << ra.dG_fwd * kcal << " vs " << r0.dG_fwd * kcal << "; basin barrier fwd "
                      << ba.barrier_fwd * kcal << " vs " << b0.barrier_fwd * kcal << " kcal/mol (shift "
                      << (ba.barrier_fwd - b0.barrier_fwd) * kcal << ", kT ln a = " << kT * std::log(a) * kcal << ")"
                      << std::endl;
            check(ra.valid && e < 1e-3, "rate unchanged under s -> a s");
            check(std::abs((ba.barrier_fwd - b0.barrier_fwd) - kT * std::log(a)) < 1e-6 * kT,
                  "basin barrier shifts by kT ln a (scale dependent)");
        }
    }

    // ---------- (4) Z-mean identity ----------
    std::cout << "(4) <Z^1/2> = N / sum w" << std::endl;
    {
        std::mt19937 rng(20261002);
        std::lognormal_distribution<double> ln(-1.0, 0.8);
        const int n = 10000;
        double sw = 0.0, swz = 0.0;
        for (int i = 0; i < n; ++i) {
            const double Z = ln(rng);
            const double w = 1.0 / std::sqrt(Z);
            sw += w;
            swz += w * std::sqrt(Z);
        }
        const double a = n / sw, b = swz / sw;
        std::cout << "    n = " << n << ": N/sum w = " << a << ", sum w Z^1/2 / sum w = " << b << ", rel "
                  << std::abs(a / b - 1.0) << std::endl;
        check(std::abs(a / b - 1.0) < 1e-12, "N / sum w == sum w_i Z_i^1/2 / sum w_i");
    }

    // ---------- (5) zHalfNear ----------
    std::cout << "(5) zHalfNear" << std::endl;
    {
        // bins of width 0.01; window s* +- 0.015
        const std::vector<double> sb = { 0.40, 0.49, 0.50, 0.51, 0.53, 0.60 };
        const std::vector<long> nb = { 100, 40, 60, 100, 500, 100 };
        const std::vector<double> vb = { 1.0, 2.0, 3.0, 4.0, 9.0, 5.0 };
        const ZHalfEstimate w = zHalfNear(sb, nb, vb, 0.50, 0.015, 50);
        const double expect = (40 * 2.0 + 60 * 3.0 + 100 * 4.0) / 200.0;
        std::cout << "    window: " << w.value << " (expected " << expect << "), samples " << w.samples
                  << ", bins " << w.bins << std::endl;
        check(w.valid && !w.fallback && std::abs(w.value - expect) < 1e-14 && w.samples == 200 && w.bins == 3,
              "sample-weighted window mean");
        const ZHalfEstimate f = zHalfNear(sb, nb, vb, 0.45, 0.015, 50);
        // left side: bin 0.40 (100 >= 25 samples); right side: 0.49 (40 >= 25) -> interpolate
        const double ef = 1.0 + (0.45 - 0.40) / (0.49 - 0.40) * (2.0 - 1.0);
        std::cout << "    fallback: " << f.value << " (expected " << ef << "), warning: " << f.warning << std::endl;
        check(f.valid && f.fallback && std::abs(f.value - ef) < 1e-14 && !f.warning.empty(),
              "fallback interpolates between the nearest populated bins");
        // side accumulation: s* = 0.495, min 200 -> left needs 100: 0.49 (40) + 0.40 (100);
        // right needs 100: 0.50 (60) + 0.51 (100)
        const ZHalfEstimate g = zHalfNear(sb, nb, vb, 0.495, 0.001, 200);
        const double slm = (40 * 0.49 + 100 * 0.40) / 140.0, vlm = (40 * 2.0 + 100 * 1.0) / 140.0;
        const double srm = (60 * 0.50 + 100 * 0.51) / 160.0, vrm = (60 * 3.0 + 100 * 4.0) / 160.0;
        const double eg = vlm + (0.495 - slm) / (srm - slm) * (vrm - vlm);
        std::cout << "    fallback (side accumulation): " << g.value << " (expected " << eg << "), samples "
                  << g.samples << ", bins " << g.bins << std::endl;
        check(g.valid && g.fallback && std::abs(g.value - eg) < 1e-13 && g.samples == 300 && g.bins == 4,
              "fallback accumulates per side until min_samples/2");
    }

    std::cout << (failures == 0 ? "ALL PASSED" : "FAILURES: " + std::to_string(failures)) << std::endl;
    return failures == 0 ? 0 : 1;
}
