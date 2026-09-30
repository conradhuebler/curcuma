/*
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This file also incorporates the public-domain-equivalent, freely
 * redistributable double-precision acos() algorithm originally developed
 * by Sun Microsystems, Inc. (fdlibm), reproduced here (via musl libc's
 * src/math/acos.c, itself "origin: FreeBSD /usr/src/lib/msun/src/e_acos.c")
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

// Claude Generated (Aug 2026): self-contained double-precision acos(),
// ported from fdlibm (via musl libc). See src/core/portable_erf.h for the
// full rationale (Wine vs. native-Windows CRT-DLL divergence). acos()
// specifically matters here because GFN-FF derives bond angles from it and
// then compares against hard-coded degree thresholds (carbene 150°, linear
// N/dihedral 170°, planarity ~344-355°) to pick a discrete hybridisation /
// bond-type / HB-classification branch -- see docs/PORTABLE_ERF.md.
// x86 / x86-64 little-endian only (curcuma's only supported targets).

#include <cmath>
#include <cstdint>
#include <cstring>

namespace CurcumaMath {
namespace acos_detail {

    inline std::uint32_t highWord(double x)
    {
        std::uint64_t bits;
        std::memcpy(&bits, &x, sizeof(bits));
        return static_cast<std::uint32_t>(bits >> 32);
    }

    inline std::uint32_t lowWord(double x)
    {
        std::uint64_t bits;
        std::memcpy(&bits, &x, sizeof(bits));
        return static_cast<std::uint32_t>(bits);
    }

    // Zero the low 32 bits of x's mantissa (~29 significant bits left) --
    // used to split sqrt(z) into two exactly-representable parts so the
    // correction term c=(z-df*df)/(s+df) doesn't lose precision to
    // cancellation.
    inline double zeroLowWord(double x)
    {
        std::uint64_t bits;
        std::memcpy(&bits, &x, sizeof(bits));
        bits &= 0xFFFFFFFF00000000ULL;
        double z;
        std::memcpy(&z, &bits, sizeof(z));
        return z;
    }

    inline constexpr double pio2_hi = 1.57079632679489655800e+00;
    inline constexpr double pio2_lo = 6.12323399573676603587e-17;
    inline constexpr double pS0 = 1.66666666666666657415e-01;
    inline constexpr double pS1 = -3.25565818622400915405e-01;
    inline constexpr double pS2 = 2.01212532134862925881e-01;
    inline constexpr double pS3 = -4.00555345006794114027e-02;
    inline constexpr double pS4 = 7.91534994289814532176e-04;
    inline constexpr double pS5 = 3.47933107596021167570e-05;
    inline constexpr double qS1 = -2.40339491173441421878e+00;
    inline constexpr double qS2 = 2.02094576023350569471e+00;
    inline constexpr double qS3 = -6.88283971605453293030e-01;
    inline constexpr double qS4 = 7.70381505559019352791e-02;

    inline double R(double z)
    {
        double p = z * (pS0 + z * (pS1 + z * (pS2 + z * (pS3 + z * (pS4 + z * pS5)))));
        double q = 1.0 + z * (qS1 + z * (qS2 + z * (qS3 + z * qS4)));
        return p / q;
    }

} // namespace acos_detail

/// Claude Generated (Aug 2026): drop-in replacement for std::acos(), see
/// header comment above -- deterministic across OS/CRT (Wine vs native
/// Windows), ported from fdlibm/musl.
inline double portable_acos(double x)
{
    using namespace acos_detail;
    const std::uint32_t hx = highWord(x);
    const std::uint32_t ix = hx & 0x7fffffffu;

    if (ix >= 0x3ff00000u) { // |x| >= 1 or NaN
        const std::uint32_t lx = lowWord(x);
        if (((ix - 0x3ff00000u) | lx) == 0) {
            // acos(1) = 0, acos(-1) = pi
            if (hx >> 31)
                return 2 * pio2_hi;
            return 0;
        }
        return 0 / (x - x); // NaN
    }
    if (ix < 0x3fe00000u) { // |x| < 0.5
        if (ix <= 0x3c600000u) // |x| < 2^-57
            return pio2_hi;
        return pio2_hi - (x - (pio2_lo - x * R(x * x)));
    }
    if (hx >> 31) { // x < -0.5
        const double z = (1.0 + x) * 0.5;
        const double s = std::sqrt(z);
        const double w = R(z) * s - pio2_lo;
        return 2 * (pio2_hi - (s + w));
    }
    // x > 0.5
    const double z = (1.0 - x) * 0.5;
    const double s = std::sqrt(z);
    const double df = zeroLowWord(s);
    const double c = (z - df * df) / (s + df);
    const double w = R(z) * s + c;
    return 2 * (df + w);
}

} // namespace CurcumaMath
