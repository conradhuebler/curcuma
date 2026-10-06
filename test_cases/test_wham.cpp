/*
 * <Validation of the WHAM solver against a known free-energy profile.>
 * Copyright (C) 2019 - 2026 Conrad Huebler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated 2026 (AI-generated, machine-tested).
 *
 * Histograms are drawn EXACTLY from p_i(b) ~ exp(-(G(b) + V_i(b))/kT) for a double-well
 * G(s) with windows of DIFFERENT stiffness. WHAM must recover G to statistical accuracy
 * when given the per-window k the windows really had. The negative control feeds it a
 * single wrong global k (the bug fixed in NebMD::writeWHAM) and must fail, which proves
 * the test is sensitive to it.
 */

#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>

#include "src/capabilities/wham.h"

static int failures = 0;
static void check(bool ok, const std::string& what)
{
    std::cout << (ok ? "  PASS  " : "  FAIL  ") << what << std::endl;
    if (!ok)
        ++failures;
}

int main()
{
    const double kT = 3.166811563e-6 * 300.0, kcal = 627.5094740631;
    const int B = 100, W = 20;
    const long N = 400000;
    // Double well, 10 kcal/mol barrier at s=0.5, asymmetric by a tilt (well B 1.5 kcal higher)
    auto G = [&](double s) {
        const double x = (s - 0.5) / 0.5;
        return (10.0 * (1.0 - x * x) + 1.5 * x) / kcal;
    };
    // (1-x^2) is a barrier hump at 0.5 of 10 kcal/mol; ends at 0 -> minima at the edges
    std::vector<double> s0(W), k(W), n(W, static_cast<double>(N));
    const double ds = 1.0 / (W - 1);
    for (int i = 0; i < W; ++i) {
        s0[i] = i * ds;
        // stiffness varies 1x..4x across windows (as after calibrateStiffness)
        k[i] = kT / std::pow(0.5 * ds, 2) * (1.0 + 3.0 * std::abs(std::sin(0.7 * i)));
    }
    std::mt19937 rng(12345);
    std::vector<std::vector<int>> hist(W, std::vector<int>(B, 0));
    for (int i = 0; i < W; ++i) {
        std::vector<double> w(B);
        for (int b = 0; b < B; ++b) {
            const double s = (b + 0.5) / B, d = s - s0[i];
            w[b] = std::exp(-(G(s) + 0.5 * k[i] * d * d) / kT);
        }
        std::discrete_distribution<int> dist(w.begin(), w.end());
        for (long t = 0; t < N; ++t)
            ++hist[i][dist(rng)];
    }

    auto maxerr = [&](const std::vector<double>& kk) {
        const std::vector<double> P = curcuma::whamSolve(hist, n, s0, kk, kT);
        std::vector<double> g(B);
        for (int b = 0; b < B; ++b)
            g[b] = -kT * std::log(P[b]) * kcal;
        // align at the profile minimum of the truth (bin 0 region is a minimum)
        const double off = g[0] - G(0.5 / B) * kcal;
        double e = 0.0;
        for (int b = 0; b < B; ++b)
            e = std::max(e, std::abs(g[b] - off - G((b + 0.5) / B) * kcal));
        return e;
    };

    const double e_ok = maxerr(k);
    std::vector<double> kwrong(W, k[0]);
    const double e_bad = maxerr(kwrong);
    std::cout << "  max |G_wham - G_true| per-window k : " << e_ok << " kcal/mol\n"
              << "  max |G_wham - G_true| global k     : " << e_bad << " kcal/mol\n";
    check(e_ok < 0.15, "per-window k recovers the profile within 0.15 kcal/mol");
    check(e_bad > 5.0 * e_ok, "a wrong global k is detected (error >5x larger)");
    return failures == 0 ? 0 : 1;
}
