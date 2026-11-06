/*
 * Symbols for terminal and file output: ASCII by default, Unicode with -DCURCUMA_UNICODE_OUTPUT=ON.
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 */

#pragma once

// Each SYM_* is a string literal, so it concatenates with its neighbours at compile time:
//   fmt::format("{:.4f} " SYM_DEG, x)
// Letters with diacritics (u-umlaut, e-acute, the A with ring above, ...) are written as plain UTF-8 and need no symbol.
// Use it in output strings only. PARAM help texts stay ASCII (the parameter generator reads them as plain literals).
// The switch is a build option (CMake: CURCUMA_UNICODE_OUTPUT); there is no run-time switch.

#ifdef CURCUMA_UNICODE_OUTPUT
#define SYM_INV "⁻¹"
#define SYM_SQ "²"
#define SYM_CUBE "³"
#define SYM_SUP1 "¹"
#define SYM_TRANSPOSE "ᵀ"
#define SYM_SUB0 "₀"
#define SYM_SUB1 "₁"
#define SYM_SUB2 "₂"
#define SYM_SUB3 "₃"
#define SYM_DEG "°"
#define SYM_SIGMA "σ"
#define SYM_PM "±"
#define SYM_LE "≤"
#define SYM_GE "≥"
#define SYM_ARROW "→"
#define SYM_TIMES "×"
#define SYM_CDOT "·"
#define SYM_PI "π"
#define SYM_LAMBDA "λ"
#define SYM_DELTA "Δ"
#define SYM_ALPHA "α"
#define SYM_KAPPA "κ"
#define SYM_SUM "Σ"
#define SYM_CHI "χ"
#define SYM_GAMMA "γ"
#define SYM_EPS "ε"
#define SYM_ZETA "ζ"
#else
#define SYM_INV "^-1"
#define SYM_SQ "^2"
#define SYM_CUBE "^3"
#define SYM_SUP1 "^1"
#define SYM_TRANSPOSE "^T"
#define SYM_SUB0 "0"
#define SYM_SUB1 "1"
#define SYM_SUB2 "2"
#define SYM_SUB3 "3"
#define SYM_DEG " deg"
#define SYM_SIGMA "sigma"
#define SYM_PM "+/-"
#define SYM_LE "<="
#define SYM_GE ">="
#define SYM_ARROW "to"
#define SYM_TIMES "x"
#define SYM_CDOT "*"
#define SYM_PI "pi"
#define SYM_LAMBDA "lambda"
#define SYM_DELTA "Delta"
#define SYM_ALPHA "alpha"
#define SYM_KAPPA "kappa"
#define SYM_SUM "Sum"
#define SYM_CHI "chi"
#define SYM_GAMMA "gamma"
#define SYM_EPS "eps"
#define SYM_ZETA "zeta"
#endif
