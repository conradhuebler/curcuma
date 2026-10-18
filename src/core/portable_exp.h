/*
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This file also incorporates the public-domain-equivalent, freely
 * redistributable double-precision exp() algorithm originally developed
 * by Sun Microsystems, Inc. (fdlibm), reproduced here (via musl libc's
 * src/math/exp.c, itself "origin: FreeBSD /usr/src/lib/msun/src/e_exp.c")
 * under its original notice:
 *
 * ====================================================
 * Copyright (C) 2004 by Sun Microsystems, Inc. All rights reserved.
 *
 * Permission to use, copy, modify, and distribute this
 * software is freely granted, provided that this notice
 * is preserved.
 * ====================================================
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
 *
 */

#pragma once

// Claude Generated (Aug 2026): self-contained double-precision exp(),
// ported from fdlibm (via musl libc). See src/core/portable_erf.h for the
// full rationale (Wine vs. native-Windows CRT-DLL divergence). exp()
// specifically matters here because GFN-FF's coordination-number (CN)
// log-compression formula log(1+exp(cnmax))-log(1+exp(cnmax-cn_raw)) uses
// it, and the resulting CN is later rounded to an integer
// (std::round(cn)) to select a hybridisation branch -- so exp() must be
// deterministic for that rounding decision to be deterministic, even
// though round() itself is exact (see docs/PORTABLE_ERF.md).
// x86 / x86-64 little-endian only (curcuma's only supported targets).
//
// Two cosmetic deviations from the original fdlibm source, both provably
// inert (no effect on the returned VALUE, only on FPU exception-flag
// side effects we don't rely on): overflow returns
// std::numeric_limits<double>::infinity() directly instead of
// deliberately overflowing x*2^1023 to infinity, and the FORCE_EVAL(...)
// calls (which exist purely to force setting the FE_INEXACT/FE_UNDERFLOW
// flags) are omitted.

#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>

namespace CurcumaMath {
namespace exp_detail {

    inline std::uint32_t highWord(double x)
    {
        std::uint64_t bits;
        std::memcpy(&bits, &x, sizeof(bits));
        return static_cast<std::uint32_t>(bits >> 32);
    }

    inline constexpr double half[2] = { 0.5, -0.5 };
    inline constexpr double ln2hi = 6.93147180369123816490e-01;
    inline constexpr double ln2lo = 1.90821492927058770002e-10;
    inline constexpr double invln2 = 1.44269504088896338700e+00;
    inline constexpr double P1 = 1.66666666666666019037e-01;
    inline constexpr double P2 = -2.77777777770155933842e-03;
    inline constexpr double P3 = 6.61375632143793436117e-05;
    inline constexpr double P4 = -1.65339022054652515390e-06;
    inline constexpr double P5 = 4.13813679705723846039e-08;

} // namespace exp_detail

/// Claude Generated (Aug 2026): drop-in replacement for std::exp(), see
/// header comment above -- deterministic across OS/CRT (Wine vs native
/// Windows), ported from fdlibm/musl.
inline double portable_exp(double x)
{
    using namespace exp_detail;
    std::uint32_t hx = highWord(x);
    const int sign = hx >> 31;
    hx &= 0x7fffffffu; // high word of |x|

    if (hx >= 0x4086232bu) { // |x| >= 708.39...
        if (std::isnan(x))
            return x;
        if (x > 709.782712893383973096)
            return std::numeric_limits<double>::infinity();
        if (x < -708.39641853226410622) {
            if (x < -745.13321910194110842)
                return 0;
            // else: falls through to the normal reduction below, which
            // correctly produces a subnormal result via std::scalbn.
        }
    }

    double x_red = x;
    int k;
    double hi, lo;
    if (hx > 0x3fd62e42u) { // |x| > 0.5 ln2
        if (hx >= 0x3ff0a2b2u) // |x| >= 1.5 ln2
            k = static_cast<int>(invln2 * x_red + half[sign]);
        else
            k = 1 - sign - sign;
        hi = x_red - k * ln2hi; // k*ln2hi is exact here
        lo = k * ln2lo;
        x_red = hi - lo;
    } else if (hx > 0x3e300000u) { // |x| > 2^-28
        k = 0;
        hi = x_red;
        lo = 0;
    } else {
        // inexact if x!=0
        return 1 + x_red;
    }

    // x_red is now in the primary range
    const double xx = x_red * x_red;
    const double c = x_red - xx * (P1 + xx * (P2 + xx * (P3 + xx * (P4 + xx * P5))));
    const double y = 1 + (x_red * c / (2 - c) - lo + hi);
    if (k == 0)
        return y;
    return std::scalbn(y, k);
}

} // namespace CurcumaMath
