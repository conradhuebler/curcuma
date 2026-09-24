/*
 * Accuracy of the Boys function used by every native-QM ERI.
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * qmint::boysFunction (grid + Taylor + recurrence, Sep 2026) and the long-double
 * qmint::boysFunctionReference are both checked against an independent series,
 *   F_n(T) = e^{-T} sum_k (2T)^k / ((2n+1)(2n+3)...(2n+2k+1)),
 * summed in long double until the terms vanish. Every n <= 16 (enough for d
 * functions and their gradient, Ltot <= 9) over 200000 random T in [0, 60] plus
 * the switch points. Pass: max relative error < 1e-13.
 *
 * Claude Generated (Sep 2026). GPL-3.0.
 */
#include "src/core/energy_calculators/qm_methods/qm_integrals.hpp"

#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

static long double seriesBoys(int n, long double T)
{
    long double term = 1.0L / (2.0L * n + 1.0L), sum = term;
    for (int k = 1; k < 100000; ++k) {
        term *= 2.0L * T / (2.0L * n + 2.0L * k + 1.0L);
        sum += term;
        if (term < sum * 1e-22L) break;
    }
    return expl(-T) * sum;
}

int main()
{
    const int maxN = 16;
    std::vector<double> Ts = { 0.0, 1e-15, 1e-10, 1e-6, 0.024999, 0.025, 0.5, 0.999999, 1.0, 1.000001,
                               17.3, 35.9999, 35.975, 36.0, 36.0001, 50.0, 60.0 };
    std::mt19937 rng(7);
    std::uniform_real_distribution<double> u(0.0, 60.0);
    for (int i = 0; i < 200000; ++i) Ts.push_back(u(rng));

    double err_fast = 0.0, err_ref = 0.0, T_fast = 0.0;
    int n_fast = 0;
    std::vector<double> F, Fr;
    for (double T : Ts) {
        qmint::boysFunction(maxN, T, F);
        qmint::boysFunctionReference(maxN, T, Fr);
        for (int n = 0; n <= maxN; ++n) {
            const long double exact = seriesBoys(n, T);
            const double ef = (double)(fabsl(F[n] - exact) / exact);
            const double er = (double)(fabsl(Fr[n] - exact) / exact);
            if (ef > err_fast) { err_fast = ef; T_fast = T; n_fast = n; }
            if (n <= 5) err_ref = std::max(err_ref, er);  // the reference is only claimed for n <= 5
        }
    }
    const bool ok = err_fast < 1e-13;
    std::printf("Boys F_n, n <= %d, %zu values of T in [0, 60]\n", maxN, Ts.size());
    std::printf("  fast (grid+Taylor):   max rel error %.2e (n=%d, T=%.6f)  %s\n", err_fast, n_fast, T_fast,
                ok ? "ok" : "FAIL");
    std::printf("  long-double reference (n <= 5): max rel error %.2e\n", err_ref);
    return ok ? 0 : 1;
}
