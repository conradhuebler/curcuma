/*
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

// Claude Generated (Sep 2026): regression test for the vendored fdlibm ports in
// src/core/portable_{erf,acos,exp,log}.h (the CURCUMA_PORTABLE_MATH layer, see
// docs/PORTABLE_ERF.md). The ports exist so that GFN-FF's discrete classification
// thresholds (CN cutoffs, angle thresholds, CN log-compression) give the same
// branch on every CRT. This test pins them to the host libm (which is correctly
// rounded to within 1 ulp on glibc) so a transcription error or an accidental
// edit of a polynomial coefficient is caught, and checks the special values.
//
// Always built and run, independent of USE_PORTABLE_MATH, because the vendored
// functions must stay correct even when the default build does not use them.

#include "src/core/portable_acos.h"
#include "src/core/portable_erf.h"
#include "src/core/portable_exp.h"
#include "src/core/portable_log.h"

#include <cmath>
#include <cstdio>
#include <limits>
#include <random>

namespace {

int g_failures = 0;

double ulp_distance(double a, double b)
{
    if (a == b) return 0.0;
    const double scale = std::fmax(std::fabs(a), std::fabs(b));
    if (scale == 0.0) return 0.0;
    return std::fabs(a - b) / (scale * std::numeric_limits<double>::epsilon());
}

template <class F, class G>
void compare_range(const char* name, F portable, G reference, double lo, double hi,
                   int n, double max_ulps)
{
    std::mt19937_64 rng(20260906u);
    std::uniform_real_distribution<double> dist(lo, hi);
    double worst = 0.0, worst_x = 0.0;
    for (int i = 0; i < n; ++i) {
        const double x = dist(rng);
        const double u = ulp_distance(portable(x), reference(x));
        if (u > worst) { worst = u; worst_x = x; }
    }
    // Endpoints and the neighbourhood of 0 / 1 (where the fdlibm branches switch).
    for (double x : { lo, hi, 0.0, 1e-300, -1e-300, 0.5, 1.0, 0.84375, 6.0, -6.0 }) {
        if (x < lo || x > hi) continue;
        const double u = ulp_distance(portable(x), reference(x));
        if (u > worst) { worst = u; worst_x = x; }
    }
    const bool ok = worst <= max_ulps;
    std::printf("%-14s [%g, %g]  worst %.3f ulp at x=%.17g  %s\n", name, lo, hi, worst, worst_x,
                ok ? "OK" : "FAIL");
    if (!ok) ++g_failures;
}

void check(const char* what, bool ok)
{
    std::printf("%-46s %s\n", what, ok ? "OK" : "FAIL");
    if (!ok) ++g_failures;
}

} // namespace

int main()
{
    using namespace CurcumaMath;
    const double inf = std::numeric_limits<double>::infinity();
    const double nan = std::numeric_limits<double>::quiet_NaN();

    // Accuracy against the host libm: the fdlibm algorithms are within 1 ulp of the
    // correctly rounded result, and so is glibc, so 2 ulp is a safe combined bound.
    compare_range("erf",  [](double x) { return portable_erf(x); },  [](double x) { return std::erf(x); },  -6.0, 6.0, 200000, 2.0);
    compare_range("acos", [](double x) { return portable_acos(x); }, [](double x) { return std::acos(x); }, -1.0, 1.0, 200000, 2.0);
    compare_range("exp",  [](double x) { return portable_exp(x); },  [](double x) { return std::exp(x); },  -700.0, 700.0, 200000, 2.0);
    compare_range("exp (small)", [](double x) { return portable_exp(x); }, [](double x) { return std::exp(x); }, -2.0, 2.0, 200000, 2.0);
    compare_range("log",  [](double x) { return portable_log(x); },  [](double x) { return std::log(x); },  1e-300, 1e300, 200000, 2.0);
    compare_range("log (near 1)", [](double x) { return portable_log(x); }, [](double x) { return std::log(x); }, 0.5, 2.0, 200000, 2.0);

    // Special values and limits.
    check("erf(0) == 0",                 portable_erf(0.0) == 0.0);
    check("erf(+inf) == 1",              portable_erf(inf) == 1.0);
    check("erf(-inf) == -1",             portable_erf(-inf) == -1.0);
    check("erf(nan) is nan",             std::isnan(portable_erf(nan)));
    check("erf is odd",                  portable_erf(0.7) == -portable_erf(-0.7));
    check("acos(1) == 0",                portable_acos(1.0) == 0.0);
    check("acos(-1) == pi",              portable_acos(-1.0) == std::acos(-1.0));
    check("acos(0) == pi/2 (1 ulp)",     ulp_distance(portable_acos(0.0), std::acos(0.0)) <= 1.0);
    check("acos(1.5) is nan",            std::isnan(portable_acos(1.5)));
    check("exp(0) == 1",                 portable_exp(0.0) == 1.0);
    check("exp(1) == e (1 ulp)",         ulp_distance(portable_exp(1.0), std::exp(1.0)) <= 1.0);
    check("exp(+inf) == inf",            portable_exp(inf) == inf);
    check("exp(-inf) == 0",              portable_exp(-inf) == 0.0);
    check("exp(800) overflows to inf",   portable_exp(800.0) == inf);
    check("exp(-800) underflows to 0",   portable_exp(-800.0) == 0.0);
    check("exp(nan) is nan",             std::isnan(portable_exp(nan)));
    check("log(1) == 0",                 portable_log(1.0) == 0.0);
    check("log(0) == -inf",              portable_log(0.0) == -inf);
    check("log(-1) is nan",              std::isnan(portable_log(-1.0)));
    check("log(+inf) == inf",            portable_log(inf) == inf);
    check("log(e) == 1 (1 ulp)",         ulp_distance(portable_log(std::exp(1.0)), 1.0) <= 1.0);

    // The GFN-FF CN log-compression, end to end, at the values that decide a rounding
    // branch: log(1+exp(cnmax)) - log(1+exp(cnmax-cn)) must agree with libm to a few ulp.
    {
        const double cnmax = 4.4;
        double worst = 0.0;
        for (double cn = 0.0; cn <= 8.0; cn += 0.0625) {
            const double a = portable_log(1.0 + portable_exp(cnmax)) - portable_log(1.0 + portable_exp(cnmax - cn));
            const double b = std::log(1.0 + std::exp(cnmax)) - std::log(1.0 + std::exp(cnmax - cn));
            worst = std::fmax(worst, std::fabs(a - b));
        }
        std::printf("CN log-compression max |diff| = %.3e\n", worst);
        check("CN log-compression within 1e-14", worst < 1e-14);
    }

    if (g_failures == 0) {
        std::printf("test_portable_math: all checks passed\n");
        return 0;
    }
    std::printf("test_portable_math: %d check(s) FAILED\n", g_failures);
    return 1;
}
