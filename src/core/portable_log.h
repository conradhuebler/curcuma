/*
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This file also incorporates the public-domain-equivalent, freely
 * redistributable double-precision log() algorithm originally developed
 * by Sun Microsystems, Inc. (fdlibm), reproduced here (via musl libc's
 * src/math/log.c, itself "origin: FreeBSD /usr/src/lib/msun/src/e_log.c")
 * under its original notice:
 *
 * ====================================================
 * Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
 *
 * Developed at SunSoft, a Sun Microsystems, Inc. business.
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

// Claude Generated (Aug 2026): self-contained double-precision log(),
// ported from fdlibm (via musl libc). See src/core/portable_erf.h for the
// full rationale (Wine vs. native-Windows CRT-DLL divergence), and
// src/core/portable_exp.h for why log() also needs vendoring: it is the
// other half of GFN-FF's CN log-compression formula that feeds a
// std::round()-based hybridisation branch.
// x86 / x86-64 little-endian only (curcuma's only supported targets).

#include <cstdint>
#include <cstring>

namespace CurcumaMath {
namespace log_detail {

    inline std::uint64_t bitsOf(double x)
    {
        std::uint64_t bits;
        std::memcpy(&bits, &x, sizeof(bits));
        return bits;
    }

    inline double fromBits(std::uint64_t bits)
    {
        double x;
        std::memcpy(&x, &bits, sizeof(x));
        return x;
    }

    inline constexpr double ln2_hi = 6.93147180369123816490e-01;
    inline constexpr double ln2_lo = 1.90821492927058770002e-10;
    inline constexpr double Lg1 = 6.666666666666735130e-01;
    inline constexpr double Lg2 = 3.999999999940941908e-01;
    inline constexpr double Lg3 = 2.857142874366239149e-01;
    inline constexpr double Lg4 = 2.222219843214978396e-01;
    inline constexpr double Lg5 = 1.818357216161805012e-01;
    inline constexpr double Lg6 = 1.531383769920937332e-01;
    inline constexpr double Lg7 = 1.479819860511658591e-01;

} // namespace log_detail

/// Claude Generated (Aug 2026): drop-in replacement for std::log(), see
/// header comment above -- deterministic across OS/CRT (Wine vs native
/// Windows), ported from fdlibm/musl.
inline double portable_log(double x)
{
    using namespace log_detail;
    std::uint64_t bits = bitsOf(x);
    std::uint32_t hx = static_cast<std::uint32_t>(bits >> 32);
    int k = 0;

    if (hx < 0x00100000u || (hx >> 31)) {
        if ((bits << 1) == 0)
            return -1 / (x * x); // log(+-0) = -inf
        if (hx >> 31)
            return (x - x) / 0.0; // log(-#) = NaN
        // subnormal number, scale x up
        k -= 54;
        x *= 18014398509481984.0; // 2^54
        bits = bitsOf(x);
        hx = static_cast<std::uint32_t>(bits >> 32);
    } else if (hx >= 0x7ff00000u) {
        return x;
    } else if (hx == 0x3ff00000u && static_cast<std::uint32_t>(bits) == 0) {
        return 0;
    }

    // reduce x into [sqrt(2)/2, sqrt(2)]
    hx += 0x3ff00000u - 0x3fe6a09eu;
    k += static_cast<int>(hx >> 20) - 0x3ff;
    hx = (hx & 0x000fffffu) + 0x3fe6a09eu;
    bits = (static_cast<std::uint64_t>(hx) << 32) | (bits & 0xffffffffULL);
    x = fromBits(bits);

    const double f = x - 1.0;
    const double hfsq = 0.5 * f * f;
    const double s = f / (2.0 + f);
    const double z = s * s;
    const double w = z * z;
    const double t1 = w * (Lg2 + w * (Lg4 + w * Lg6));
    const double t2 = z * (Lg1 + w * (Lg3 + w * (Lg5 + w * Lg7)));
    const double R = t2 + t1;
    const double dk = k;
    return s * (hfsq + R) + dk * ln2_lo - hfsq + f + dk * ln2_hi;
}

} // namespace CurcumaMath
