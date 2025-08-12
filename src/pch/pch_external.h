/*
 * External Libraries Precompiled Header
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#ifdef USE_D3
// DFT-D3 dispersion correction headers
#include "dftd_damping.h"
#include "dftd_dispersion.h"
#include "dftd_econv.h"
#include "dftd_geometry.h"
#endif

#ifdef USE_D4
// DFT-D4 dispersion correction headers
#include "dftd_damping.h"
#include "dftd_dispersion.h"
#include "dftd_econv.h"
#include "dftd_geometry.h"
#endif