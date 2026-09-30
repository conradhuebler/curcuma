/*
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This file also incorporates the public-domain-equivalent, freely
 * redistributable double-precision erf()/erfc() algorithm originally
 * developed by Sun Microsystems, Inc. (fdlibm), reproduced here (via
 * musl libc's src/math/erf.c, itself "origin: FreeBSD
 * /usr/src/lib/msun/src/s_erf.c") under its original notice:
 *
 * ====================================================
 * Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
 *
 * Developed at SunPro, a Sun Microsystems, Inc. business.
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

// Claude Generated (Aug 2026): self-contained double-precision erf(),
// ported from fdlibm (via musl libc). Statically compiled machine code
// -- unlike std::erf(), which on MinGW/Windows builds resolves to a
// dynamic import from api-ms-win-crt-math-l1-1-0.dll. Wine ships its
// own reimplementation of that DLL (it cannot redistribute Microsoft's
// binary), whose erf() rounds differently in the last bit than native
// Windows' ucrtbase erf(). GFN-FF builds several hard threshold
// decisions (hybridisation, bond type, pi-membership) directly on
// CN sums that use erf(), so a single ULP difference can flip a
// branch and change the resulting bond term (see docs/PORTABLE_ERF.md).
// Enable with -DUSE_PORTABLE_MATH=ON; see src/core/math_compat.h.
//
// x86 / x86-64 little-endian only (curcuma's only supported targets):
// the algorithm inspects a double's raw IEEE-754 bit pattern directly.

#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>

namespace CurcumaMath {
namespace erf_detail {

    inline std::uint32_t highWord(double x)
    {
        std::uint64_t bits;
        std::memcpy(&bits, &x, sizeof(bits));
        return static_cast<std::uint32_t>(bits >> 32);
    }

    // Returns x with the low 32 bits of its mantissa zeroed (~29 significant
    // bits left) -- used to split exp(-x*x) into two exactly-representable
    // factors without losing precision to cancellation (fdlibm Note1).
    inline double zeroLowWord(double x)
    {
        std::uint64_t bits;
        std::memcpy(&bits, &x, sizeof(bits));
        bits &= 0xFFFFFFFF00000000ULL;
        double z;
        std::memcpy(&z, &bits, sizeof(z));
        return z;
    }

    // Coefficients for approximation to erf on [0, 0.84375]
    inline constexpr double erx = 8.45062911510467529297e-01;
    inline constexpr double efx8 = 1.02703333676410069053e+00;
    inline constexpr double pp0 = 1.28379167095512558561e-01;
    inline constexpr double pp1 = -3.25042107247001499370e-01;
    inline constexpr double pp2 = -2.84817495755985104766e-02;
    inline constexpr double pp3 = -5.77027029648944159157e-03;
    inline constexpr double pp4 = -2.37630166566501626084e-05;
    inline constexpr double qq1 = 3.97917223959155352819e-01;
    inline constexpr double qq2 = 6.50222499887672944485e-02;
    inline constexpr double qq3 = 5.08130628187576562776e-03;
    inline constexpr double qq4 = 1.32494738004321644526e-04;
    inline constexpr double qq5 = -3.96022827877536812320e-06;

    // Coefficients for approximation to erf on [0.84375, 1.25]
    inline constexpr double pa0 = -2.36211856075265944077e-03;
    inline constexpr double pa1 = 4.14856118683748331666e-01;
    inline constexpr double pa2 = -3.72207876035701323847e-01;
    inline constexpr double pa3 = 3.18346619901161753674e-01;
    inline constexpr double pa4 = -1.10894694282396677476e-01;
    inline constexpr double pa5 = 3.54783043256182359371e-02;
    inline constexpr double pa6 = -2.16637559486879084300e-03;
    inline constexpr double qa1 = 1.06420880400844228286e-01;
    inline constexpr double qa2 = 5.40397917702171048937e-01;
    inline constexpr double qa3 = 7.18286544141962662868e-02;
    inline constexpr double qa4 = 1.26171219808761642112e-01;
    inline constexpr double qa5 = 1.36370839120290507362e-02;
    inline constexpr double qa6 = 1.19844998467991074170e-02;

    // Coefficients for approximation to erfc on [1.25, 1/0.35]
    inline constexpr double ra0 = -9.86494403484714822705e-03;
    inline constexpr double ra1 = -6.93858572707181764372e-01;
    inline constexpr double ra2 = -1.05586262253232909814e+01;
    inline constexpr double ra3 = -6.23753324503260060396e+01;
    inline constexpr double ra4 = -1.62396669462573470355e+02;
    inline constexpr double ra5 = -1.84605092906711035994e+02;
    inline constexpr double ra6 = -8.12874355063065934246e+01;
    inline constexpr double ra7 = -9.81432934416914548592e+00;
    inline constexpr double sa1 = 1.96512716674392571292e+01;
    inline constexpr double sa2 = 1.37657754143519042600e+02;
    inline constexpr double sa3 = 4.34565877475229228821e+02;
    inline constexpr double sa4 = 6.45387271733267880336e+02;
    inline constexpr double sa5 = 4.29008140027567833386e+02;
    inline constexpr double sa6 = 1.08635005541779435134e+02;
    inline constexpr double sa7 = 6.57024977031928170135e+00;
    inline constexpr double sa8 = -6.04244152148580987438e-02;

    // Coefficients for approximation to erfc on [1/0.35, 28]
    inline constexpr double rb0 = -9.86494292470009928597e-03;
    inline constexpr double rb1 = -7.99283237680523006574e-01;
    inline constexpr double rb2 = -1.77579549177547519889e+01;
    inline constexpr double rb3 = -1.60636384855821916062e+02;
    inline constexpr double rb4 = -6.37566443368389627722e+02;
    inline constexpr double rb5 = -1.02509513161107724954e+03;
    inline constexpr double rb6 = -4.83519191608651397019e+02;
    inline constexpr double sb1 = 3.03380607434824582924e+01;
    inline constexpr double sb2 = 3.25792512996573918826e+02;
    inline constexpr double sb3 = 1.53672958608443695994e+03;
    inline constexpr double sb4 = 3.19985821950859553908e+03;
    inline constexpr double sb5 = 2.55305040643316442583e+03;
    inline constexpr double sb6 = 4.74528541206955367215e+02;
    inline constexpr double sb7 = -2.24409524465858183362e+01;

    // erfc on [0.84375, 1.25] via Taylor expansion about x=1
    inline double erfc1(double x)
    {
        double s = std::fabs(x) - 1;
        double P = pa0 + s * (pa1 + s * (pa2 + s * (pa3 + s * (pa4 + s * (pa5 + s * pa6)))));
        double Q = 1 + s * (qa1 + s * (qa2 + s * (qa3 + s * (qa4 + s * (qa5 + s * qa6)))));
        return 1 - erx - P / Q;
    }

    // erfc on [1.25, 28] via the asymptotic series exp(-x*x)/(x*sqrt(pi))
    inline double erfc2(std::uint32_t ix, double x)
    {
        if (ix < 0x3ff40000u) // |x| < 1.25
            return erfc1(x);

        x = std::fabs(x);
        double s = 1 / (x * x);
        double R, S;
        if (ix < 0x4006db6du) { // |x| < 1/0.35 ~ 2.85714
            R = ra0 + s * (ra1 + s * (ra2 + s * (ra3 + s * (ra4 + s * (ra5 + s * (ra6 + s * ra7))))));
            S = 1.0 + s * (sa1 + s * (sa2 + s * (sa3 + s * (sa4 + s * (sa5 + s * (sa6 + s * (sa7 + s * sa8)))))));
        } else { // |x| > 1/0.35
            R = rb0 + s * (rb1 + s * (rb2 + s * (rb3 + s * (rb4 + s * (rb5 + s * rb6)))));
            S = 1.0 + s * (sb1 + s * (sb2 + s * (sb3 + s * (sb4 + s * (sb5 + s * (sb6 + s * sb7))))));
        }
        double z = zeroLowWord(x);
        return std::exp(-z * z - 0.5625) * std::exp((z - x) * (z + x) + R / S) / x;
    }

} // namespace erf_detail

/// Claude Generated (Aug 2026): drop-in replacement for std::erf(), see
/// header comment above -- deterministic across OS/CRT (Wine vs native
/// Windows), ported from fdlibm/musl, accurate to ~1e-16 vs std::erf().
inline double portable_erf(double x)
{
    using namespace erf_detail;
    std::uint32_t ix = highWord(x);
    const int sign = ix >> 31;
    ix &= 0x7fffffffu;
    if (ix >= 0x7ff00000u) {
        // erf(nan) = nan, erf(+-inf) = +-1
        return 1 - 2 * sign + 1 / x;
    }
    double y;
    if (ix < 0x3feb0000u) { // |x| < 0.84375
        if (ix < 0x3e300000u) // |x| < 2^-28: avoid underflow
            return 0.125 * (8 * x + efx8 * x);
        double z = x * x;
        double r = pp0 + z * (pp1 + z * (pp2 + z * (pp3 + z * pp4)));
        double s = 1.0 + z * (qq1 + z * (qq2 + z * (qq3 + z * (qq4 + z * qq5))));
        return x + x * (r / s);
    }
    if (ix < 0x40180000u) // 0.84375 <= |x| < 6
        y = 1 - erfc2(ix, x);
    else
        y = 1 - std::numeric_limits<double>::min();
    return sign ? -y : y;
}

} // namespace CurcumaMath
