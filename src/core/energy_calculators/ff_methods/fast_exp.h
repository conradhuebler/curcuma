/*
 * Fast vectorized exp(-x^2) for GFN-FF coordination-number gradient
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
 * Claude Generated (May 2026, WP-D Stage B): Vectorized exp(-x^2) used in the
 * GFN-FF dcn step 3 inner loop. Range reduction via 2^q * exp(r), Remez/Taylor
 * polynomial of degree 13 over r in [-ln2/2, ln2/2]. AVX2 path is exercised
 * when both __AVX2__ and GFNFF_FAST_EXP are defined; otherwise the scalar
 * fallback (std::exp) is used so the reference path stays bit-identical.
 */

#pragma once

#include <cstddef>

namespace curcuma {
namespace gfnff {

// Vectorized batch: out[k] = exp(-x_sq[k]) for k in [0, n).
// x_sq must be non-negative (caller passes kn^2 * dr^2). Inputs above ~700
// underflow to 0 — the AVX2 path masks them out before forming the result.
// Buffers may alias (read x_sq[k] and write out[k] for the same k is safe).
void fast_exp_neg_sq_block(const double* x_sq, double* out, std::size_t n) noexcept;

} // namespace gfnff
} // namespace curcuma
