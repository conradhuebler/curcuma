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
 *
 */

#pragma once

// Claude Generated (Aug 2026): single call site to switch erf()/acos()/
// exp()/log() between the platform CRT/libm (default) and curcuma's
// self-contained portable implementations (-DUSE_PORTABLE_MATH=ON, see
// src/core/portable_{erf,acos,exp,log}.h).
//
// Background: on MinGW/Windows builds, std::erf()/acos()/exp()/log()
// dynamically import from api-ms-win-crt-math-l1-1-0.dll; Wine's
// reimplementation of that DLL (it cannot redistribute Microsoft's
// binary) rounds these functions differently in the last bit than native
// Windows' ucrtbase. GFN-FF/EEQ build many hard-coded classification
// thresholds directly on their output -- erf() feeds the CN sum compared
// against bond/pi-membership cutoffs; acos()-derived bond angles are
// compared against fixed degree thresholds (carbene 150 deg, linear N
// 170 deg, planarity ~344-355 deg) to pick sp/sp2/sp3 hybridisation,
// bond type, or HB donor/acceptor class; exp()/log() feed the CN
// log-compression formula whose result is later std::round()-ed to an
// integer that gates further hybridisation branches (round() itself is
// exact per IEEE 754 -- its *input* must be deterministic instead). A
// single ULP difference in any of these can flip a discrete branch and
// produce a genuinely different bond term for the same input geometry.
// See docs/PORTABLE_ERF.md.
//
// NOT vendored: std::round()/lround()/llround() -- unlike the
// transcendental functions above, round-to-nearest is a "basic
// operation" IEEE 754 requires to be exactly correctly rounded for every
// finite input, so no conformant libm (Wine's included) can legally
// disagree with another on it. sqrt() is likewise not vendored: on
// x86/x86-64 it compiles to the hardware SQRTSD instruction, which
// IEEE 754 also mandates to be exactly correctly rounded.

#ifdef CURCUMA_PORTABLE_MATH
#include "src/core/portable_erf.h"
#include "src/core/portable_acos.h"
#include "src/core/portable_exp.h"
#include "src/core/portable_log.h"
inline double curcuma_erf(double x) { return CurcumaMath::portable_erf(x); }
inline double curcuma_acos(double x) { return CurcumaMath::portable_acos(x); }
inline double curcuma_exp(double x) { return CurcumaMath::portable_exp(x); }
inline double curcuma_log(double x) { return CurcumaMath::portable_log(x); }
#else
#include <cmath>
inline double curcuma_erf(double x) { return std::erf(x); }
inline double curcuma_acos(double x) { return std::acos(x); }
inline double curcuma_exp(double x) { return std::exp(x); }
inline double curcuma_log(double x) { return std::log(x); }
#endif
