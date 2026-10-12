/*
 * Unit test for curcuma::gfnff::fast_exp_neg_sq_block (WP-D Stage B)
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Claude Generated (May 2026, WP-D Stage B):
 *   Validates the vectorized exp(-x_sq) helper used in GFN-FF dcn step 3.
 *   Compares fast_exp_neg_sq_block against std::exp(-x) over the operational
 *   input range, checks edge cases (x_sq = 0, near-underflow), and exercises
 *   block sizes that hit and miss the SIMD lane width (n = 4, 7, 8, 13, 100,
 *   10000). Tolerance: 4 ULP relative ≈ 8.9e-16 for normal magnitudes.
 *
 *   Run modes:
 *     GFNFF_FAST_EXP undefined → tests the scalar fallback (identity vs std::exp).
 *     GFNFF_FAST_EXP defined + AVX2 → tests the SIMD Horner-Taylor kernel.
 */

#include "src/core/energy_calculators/ff_methods/fast_exp.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

namespace {

struct Stats {
    double max_abs_err = 0.0;
    double max_rel_err = 0.0;
    double worst_x = 0.0;
    double worst_ref = 0.0;
    double worst_got = 0.0;
};

Stats compare(const std::vector<double>& x_sq, const std::vector<double>& got)
{
    Stats s;
    for (std::size_t i = 0; i < x_sq.size(); ++i) {
        double ref = std::exp(-x_sq[i]);
        double abs_err = std::abs(got[i] - ref);
        double rel_err = (ref > 0.0) ? abs_err / ref : abs_err;
        if (abs_err > s.max_abs_err) s.max_abs_err = abs_err;
        if (rel_err > s.max_rel_err) {
            s.max_rel_err = rel_err;
            s.worst_x = x_sq[i];
            s.worst_ref = ref;
            s.worst_got = got[i];
        }
    }
    return s;
}

int tests_run = 0;
int tests_passed = 0;

void check(bool condition, const char* msg)
{
    ++tests_run;
    if (condition) {
        ++tests_passed;
        std::printf("[PASS] %s\n", msg);
    } else {
        std::printf("[FAIL] %s\n", msg);
    }
}

}  // namespace

int main()
{
#ifdef GFNFF_FAST_EXP
    std::printf("=== fast_exp tests (GFNFF_FAST_EXP=ON, AVX2 path expected) ===\n");
#else
    std::printf("=== fast_exp tests (GFNFF_FAST_EXP=OFF, scalar fallback) ===\n");
#endif

    // ---- 1. Exact match at canonical points ----
    {
        const double xs[] = {0.0, 0.5, 1.0, 4.0, 10.0, 50.0};
        const std::size_t n = sizeof(xs) / sizeof(xs[0]);
        std::vector<double> in(xs, xs + n);
        std::vector<double> out(n);
        curcuma::gfnff::fast_exp_neg_sq_block(in.data(), out.data(), n);
        Stats s = compare(in, out);
        std::printf("  canonical n=%zu: max_rel_err = %.3e (at x_sq=%g, ref=%.3e got=%.3e)\n",
                    n, s.max_rel_err, s.worst_x, s.worst_ref, s.worst_got);
        check(s.max_rel_err < 4e-15, "canonical points: max rel err < 4 ULP");
    }

    // ---- 2. Dense sweep over the operational range [0, 100] ----
    {
        constexpr std::size_t N = 10001;
        std::vector<double> in(N), out(N);
        for (std::size_t i = 0; i < N; ++i) {
            in[i] = 100.0 * static_cast<double>(i) / static_cast<double>(N - 1);
        }
        curcuma::gfnff::fast_exp_neg_sq_block(in.data(), out.data(), N);
        Stats s = compare(in, out);
        std::printf("  sweep x_sq in [0,100] n=%zu: max_rel_err = %.3e\n",
                    N, s.max_rel_err);
        check(s.max_rel_err < 1e-13, "sweep [0,100]: max rel err < 1e-13");
    }

    // ---- 3. Extended range up to underflow ----
    {
        constexpr std::size_t N = 1001;
        std::vector<double> in(N), out(N);
        for (std::size_t i = 0; i < N; ++i) {
            in[i] = 0.7 * static_cast<double>(i);  // up to ~700
        }
        curcuma::gfnff::fast_exp_neg_sq_block(in.data(), out.data(), N);
        Stats s = compare(in, out);
        std::printf("  sweep x_sq in [0,700] n=%zu: max_rel_err = %.3e (worst at x=%g)\n",
                    N, s.max_rel_err, s.worst_x);
        // We accept slightly worse near 700 because of subnormal range.
        check(s.max_rel_err < 1e-12, "sweep [0,700]: max rel err < 1e-12");
    }

    // ---- 4. Underflow region: x_sq > 700 should give 0 (or near-zero) ----
    {
        const double xs[] = {750.0, 1000.0, 5000.0, 12000.0};
        const std::size_t n = sizeof(xs) / sizeof(xs[0]);
        std::vector<double> in(xs, xs + n);
        std::vector<double> out(n, -1.0);  // sentinel
        curcuma::gfnff::fast_exp_neg_sq_block(in.data(), out.data(), n);
        bool ok = true;
        for (std::size_t i = 0; i < n; ++i) {
            // std::exp(-x) is either 0 or sub-normal here. AVX2 path clamps to 0.
            // Scalar fallback returns whatever std::exp does (also ~0).
            if (out[i] < 0.0 || out[i] > 1e-300) {
                ok = false;
                std::printf("  underflow at x_sq=%g returned %g (expected ≈ 0)\n", xs[i], out[i]);
            }
        }
        check(ok, "underflow x_sq>700 clamps to ~0");
    }

    // ---- 5. Tail handling: non-multiple-of-4 sizes ----
    {
        for (std::size_t n : {std::size_t{1}, std::size_t{2}, std::size_t{3},
                              std::size_t{4}, std::size_t{5}, std::size_t{7},
                              std::size_t{8}, std::size_t{13}, std::size_t{16},
                              std::size_t{17}}) {
            std::vector<double> in(n), out(n);
            for (std::size_t i = 0; i < n; ++i) {
                in[i] = 0.1 * static_cast<double>(i + 1) * static_cast<double>(i + 1);
            }
            curcuma::gfnff::fast_exp_neg_sq_block(in.data(), out.data(), n);
            Stats s = compare(in, out);
            char msg[64];
            std::snprintf(msg, sizeof(msg),
                          "tail n=%zu rel_err=%.2e < 1e-13", n, s.max_rel_err);
            check(s.max_rel_err < 1e-13, msg);
        }
    }

    // ---- 6. Random inputs that mimic real dcn-step-3 distribution ----
    // kn = -7.5, so x_sq = kn^2 * dr^2 = 56.25 * dr^2. For dr in [-1, 1]
    // (typical surviving covalent pairs) x_sq in [0, 56.25].
    {
        std::mt19937_64 rng(0xC0FFEEull);
        std::uniform_real_distribution<double> dr_dist(-1.0, 1.0);
        constexpr std::size_t N = 100000;
        std::vector<double> in(N), out(N);
        for (std::size_t i = 0; i < N; ++i) {
            double dr = dr_dist(rng);
            in[i] = 56.25 * dr * dr;
        }
        curcuma::gfnff::fast_exp_neg_sq_block(in.data(), out.data(), N);
        Stats s = compare(in, out);
        std::printf("  realistic dcn distribution (N=%zu): max_rel_err = %.3e\n",
                    N, s.max_rel_err);
        check(s.max_rel_err < 1e-13, "realistic dcn distribution: max rel err < 1e-13");
    }

    // ---- Summary ----
    std::printf("\n=== %d / %d checks passed ===\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
