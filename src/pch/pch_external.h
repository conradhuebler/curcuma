/*
 * External Libraries Precompiled Header
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This header precompiles external computational chemistry libraries
 * and threading utilities that are commonly used across the project.
 */

#pragma once

// Include base PCH first
#include "pch_base.h"

// Threading utilities
#include "external/CxxThreadPool/include/CxxThreadPool.hpp"

#ifdef USE_XTB
// XTB interface headers would go here
// Note: These may need to be included conditionally based on Fortran bindings
#endif

#ifdef USE_TBLITE
// TBLite interface headers would go here
// Note: These may need to be included conditionally based on Fortran bindings
#endif

// DFT-D3/D4 (s-dftd3 / cpp-d4) headers are NOT precompiled here.
// dftd_econv.h defines GLOBAL unit-conversion constants (kcaltoau, aatoau,
// autokcal, ...) that collide with names in ALPBParameters / other namespaces
// pulled into unrelated TUs via this PCH (e.g. alpb_solvation.cpp's
// `using namespace ALPBParameters;` → ambiguous reference, exposed by gcc-16).
// The D3/D4 interface TUs (dftd3interface.*, dftd4interface.*, dftd4_helper.cpp)
// include the dftd_* headers directly and are themselves USE_D3/USE_D4-gated,
// so dropping them from the PCH loses only precompilation, not functionality.
// See docs/TECHNICAL_DEBT.md Q-31. (Claude Generated fix, 2026-06)