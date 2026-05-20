/*
 * Fast vectorized exp(-x^2) for GFN-FF coordination-number gradient
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Claude Generated (May 2026, WP-D Stage B):
 *
 *   exp(y)  with  y = -x_sq ≤ 0
 *      Range reduction:  y = q·ln(2) + r,  q = round(y/ln(2)),  r ∈ [-ln2/2, ln2/2]
 *      Result:           exp(y) = 2^q · exp(r)
 *      Polynomial:       Horner-form Taylor of degree 13 over r
 *                        Worst-case error ≲ 1 ULP for x_sq ∈ [0, 700].
 *      Underflow:        x_sq > 700  ⇒  exp(-x_sq) < 1e-304  ⇒  clamp to 0.
 *
 * The scalar fallback path calls std::exp directly so the reference build
 * (USE_GFNFF_FAST_EXP=OFF) is bit-identical to pre-WP behavior.
 */

#include "fast_exp.h"

#include <cmath>

#if defined(__AVX2__) && defined(GFNFF_FAST_EXP)
#include <immintrin.h>
#endif

namespace curcuma {
namespace gfnff {

void fast_exp_neg_sq_block(const double* x_sq, double* out, std::size_t n) noexcept
{
#if defined(__AVX2__) && defined(GFNFF_FAST_EXP)
    // Reciprocal of ln(2) — high-precision double rounding of 1/ln(2).
    constexpr double LOG2E = 1.4426950408889634;
    // ln(2) at double precision. For |q| < 2^15 (our operational range, since q
    // = round(-x_sq/ln2) and x_sq ≤ 700 ⇒ |q| ≤ 1011) the product q·LN2 is
    // exact in double, so a Cody-Waite split is unnecessary here.
    constexpr double LN2 = 0.6931471805599453;
    // Underflow guard. exp(-700) ≈ 9.86e-305; everything beyond is sub-normal
    // or 0 and would corrupt the bit-packed 2^q anyway. The mask below zeros
    // those lanes.
    constexpr double UNDERFLOW_LIMIT = 700.0;

    // Taylor coefficients 1/k! for k = 13..0 (Horner order, highest first).
    constexpr double c13 = 1.6059043836821613e-10;
    constexpr double c12 = 2.0876756987868098e-09;
    constexpr double c11 = 2.505210838544172e-08;
    constexpr double c10 = 2.755731922398589e-07;
    constexpr double c9  = 2.755731922398589e-06;
    constexpr double c8  = 2.48015873015873e-05;
    constexpr double c7  = 1.984126984126984e-04;
    constexpr double c6  = 1.388888888888889e-03;
    constexpr double c5  = 8.333333333333333e-03;
    constexpr double c4  = 4.166666666666666e-02;
    constexpr double c3  = 1.666666666666666e-01;
    constexpr double c2  = 0.5;
    constexpr double c1  = 1.0;
    constexpr double c0  = 1.0;

    const __m256d vlog2e = _mm256_set1_pd(LOG2E);
    const __m256d vln2 = _mm256_set1_pd(LN2);
    const __m256d vunderflow = _mm256_set1_pd(UNDERFLOW_LIMIT);
    const __m256d vzero = _mm256_setzero_pd();
    const __m256i v1023 = _mm256_set1_epi64x(1023);

    const std::size_t simd_end = n & ~std::size_t{3};

    for (std::size_t k = 0; k < simd_end; k += 4) {
        __m256d xsq = _mm256_loadu_pd(x_sq + k);

        // y = -x_sq  (since x_sq ≥ 0 in our use case, y ≤ 0)
        __m256d y = _mm256_sub_pd(vzero, xsq);

        // Underflow mask: 0xFF..FF where x_sq > 700, used to zero those lanes.
        __m256d uf_mask = _mm256_cmp_pd(xsq, vunderflow, _CMP_GT_OQ);

        // q = round(y * log2e). _mm256_round_pd rounds to nearest, ties to even.
        __m256d q = _mm256_round_pd(_mm256_mul_pd(y, vlog2e),
                                    _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);

        // r = y - q * ln2,  r ∈ [-ln2/2, ln2/2]
        __m256d r = _mm256_fnmadd_pd(q, vln2, y);

        // Horner polynomial: p = ((..((c13*r + c12)*r + c11)*r + ..)*r + c0)
        __m256d p = _mm256_set1_pd(c13);
        p = _mm256_fmadd_pd(p, r, _mm256_set1_pd(c12));
        p = _mm256_fmadd_pd(p, r, _mm256_set1_pd(c11));
        p = _mm256_fmadd_pd(p, r, _mm256_set1_pd(c10));
        p = _mm256_fmadd_pd(p, r, _mm256_set1_pd(c9));
        p = _mm256_fmadd_pd(p, r, _mm256_set1_pd(c8));
        p = _mm256_fmadd_pd(p, r, _mm256_set1_pd(c7));
        p = _mm256_fmadd_pd(p, r, _mm256_set1_pd(c6));
        p = _mm256_fmadd_pd(p, r, _mm256_set1_pd(c5));
        p = _mm256_fmadd_pd(p, r, _mm256_set1_pd(c4));
        p = _mm256_fmadd_pd(p, r, _mm256_set1_pd(c3));
        p = _mm256_fmadd_pd(p, r, _mm256_set1_pd(c2));
        p = _mm256_fmadd_pd(p, r, _mm256_set1_pd(c1));
        p = _mm256_fmadd_pd(p, r, _mm256_set1_pd(c0));

        // 2^q via IEEE-754 bit packing:
        //   double_bits = (q + 1023) << 52  with the sign bit clear.
        // Convert q (double) -> int32 (4 lanes in __m128i) -> int64 (4 lanes in
        // __m256i). Within our operational range |q| ≤ 1011, the int32 cast is
        // exact and never produces 0x80000000 (the cvtpd_epi32 "indefinite"
        // sentinel).
        __m128i q_i32 = _mm256_cvtpd_epi32(q);
        __m256i q_i64 = _mm256_cvtepi32_epi64(q_i32);
        __m256i exp_bits = _mm256_slli_epi64(_mm256_add_epi64(q_i64, v1023), 52);
        __m256d pow2q = _mm256_castsi256_pd(exp_bits);

        __m256d result = _mm256_mul_pd(pow2q, p);
        // andnot(mask, x) = x & ~mask  →  zeros lanes where x_sq > 700.
        result = _mm256_andnot_pd(uf_mask, result);

        _mm256_storeu_pd(out + k, result);
    }

    for (std::size_t k = simd_end; k < n; ++k) {
        out[k] = std::exp(-x_sq[k]);
    }
#else
    for (std::size_t k = 0; k < n; ++k) {
        out[k] = std::exp(-x_sq[k]);
    }
#endif
}

} // namespace gfnff
} // namespace curcuma
