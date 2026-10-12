/*
 * test_frag_charge_selection.cpp - frag_charge_model ensemble: bounded placement selection
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (Oct 2026) - AI-generated, machine-tested; human production testing pending.
 *
 * GFNFF::fragSelectPlacements() (dynamic programming, O(G nq^2) states) replaced the full
 * enumeration of all C(nq+G-1, nq) charge distributions, which ran out of memory for many
 * fragments (497 groups, charge 3: 2e7 vectors, ~41 GB). The selection must return the SAME
 * placements in the SAME order as the enumeration (GFNFF::fragEnumeratePlacements(), the old
 * code) - including the tie-break by enumeration order, which decides who enters the
 * frag_charge_max_placements cap when free-charge scores are equal.
 *
 * Part 1: random equivalence, fixed seed, N_CASES cases (G <= 12, nq <= 4, both signs, cap
 *         1..8, electron counts that make radicals and bare nuclei, free charges drawn as
 *         continuous values, exact duplicates, values a few ulp apart, and zeros). 0 mismatches.
 * Part 2: scale - G = 497 / 2000 groups, nq up to 6: must finish quickly (the enumeration
 *         would need C(nq+G-1, nq) vectors); prints time and the size of the result.
 *
 * Usage: test_frag_charge_selection [n_cases]
 */
#include "src/core/energy_calculators/ff_methods/gfnff.h"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>

int main(int argc, char** argv)
{
    const long N_CASES = argc > 1 ? std::atol(argv[1]) : 20000;
    std::mt19937_64 rng(20261010ULL);
    auto uni = [&](int lo, int hi) { return std::uniform_int_distribution<int>(lo, hi)(rng); };
    auto real = [&](double lo, double hi) { return std::uniform_real_distribution<double>(lo, hi)(rng); };

    long mismatches = 0, n_capped = 0, n_ties = 0, n_bare = 0, n_inexact = 0;
    for (long t = 0; t < N_CASES; ++t) {
        const int G = uni(1, 12);
        const int nq = uni(1, 4);
        const double unit = uni(0, 1) ? 1.0 : -1.0;
        const int K = uni(1, 8);
        std::vector<long> nelec(G);
        for (auto& n : nelec) {
            const int mode = uni(0, 3);
            n = mode == 0 ? uni(0, 3) : (mode == 1 ? 2 * uni(1, 10) : uni(1, 40));   // small (bare nuclei), even, any
        }
        std::vector<double> q(G);
        const int qmode = uni(0, 4);
        std::vector<double> pool = { real(-1, 1), real(-1, 1), real(-1, 1) };
        for (int g = 0; g < G; ++g) {
            switch (qmode) {
            case 0: q[g] = real(-1.0, 1.0); break;                           // continuous
            case 1: q[g] = pool[uni(0, 2)]; break;                            // exact duplicates (symmetric fragments)
            case 2: {                                                          // a few ulp apart
                double x = pool[0];
                for (int s = uni(0, 4); s > 0; --s) x = std::nextafter(x, 2.0);
                q[g] = x;
                break;
            }
            case 3: q[g] = 0.0; break;                                         // all tied
            default: q[g] = std::ldexp(static_cast<double>(uni(-8, 8)), -3); break;   // exactly representable, many sums tie
            }
        }
        const auto ref = GFNFF::fragEnumeratePlacements(nq, unit, nelec, q, K);
        bool exact = true;
        const auto sel = GFNFF::fragSelectPlacements(nq, unit, nelec, q, K, &exact);
        if (!exact) ++n_inexact;
        if (qmode != 0) ++n_ties;
        for (const auto& v : ref)
            for (int g = 0; g < G; ++g)
                if (nelec[g] - static_cast<long>(unit * v[g]) <= 0) { ++n_bare; goto counted; }
    counted:
        // was the cap active? (more minimum-rank placements than K)
        if (static_cast<int>(GFNFF::fragEnumeratePlacements(nq, unit, nelec, q, 1000000).size()) > K) ++n_capped;
        if (ref != sel) {
            if (++mismatches <= 5) {
                std::printf("MISMATCH case %ld: G %d nq %d unit %+.0f K %d qmode %d (ref %zu, sel %zu)\n", t, G, nq, unit, K, qmode,
                    ref.size(), sel.size());
            }
        }
    }
    std::printf("equivalence: %ld cases, %ld mismatches (cap active in %ld, tie modes in %ld, bare-nucleus placements in %ld, "
                "inexact tie-break flag in %ld)\n",
        N_CASES, mismatches, n_capped, n_ties, n_bare, n_inexact);

    // Part 2: scale
    bool scale_ok = true;
    for (const auto& spec : std::vector<std::pair<int, int>>{ { 497, 1 }, { 497, 2 }, { 497, 3 }, { 497, 4 }, { 2000, 6 } }) {
        const int G = spec.first, nq = spec.second;
        std::vector<long> nelec(G);
        std::vector<double> q(G);
        for (int g = 0; g < G; ++g) {
            nelec[g] = (g % 11 == 0) ? 11 : 10;   // a few odd-electron groups (Na-like), the rest water-like
            q[g] = real(-0.5, 0.5);
        }
        const auto t0 = std::chrono::steady_clock::now();
        const auto sel = GFNFF::fragSelectPlacements(nq, 1.0, nelec, q, 6);
        const double dt = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        std::printf("scale: G %d nq %d -> %zu placements in %.3f s\n", G, nq, sel.size(), dt);
        if (sel.empty() || dt > 10.0) scale_ok = false;
    }
    const bool ok = mismatches == 0 && scale_ok && N_CASES > 0;
    std::printf("%s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
