/*
 * Ulysses Precompiled Header - Header-Only Ulysses Library
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This header precompiles the heavy header-only Ulysses library
 * to significantly speed up compilation times for QM method builds.
 *
 * NOTE: Symbol conflicts with Curcuma exist but are acceptable
 * for the compilation speed benefit. Namespace wrapping doesn't work
 * due to global variable definitions in Ulysses headers.
 */

#pragma once

// Include base PCH first
#include "pch_base.h"

#ifdef USE_ULYSSES

// Core Ulysses headers (header-only library)
// These provide the biggest compilation speedup
#include "BSet.hpp"
#include "GFN.hpp"
#include "MNDOd.hpp"
#include "Molecule.hpp"
#include "QC.hpp"

#endif // USE_ULYSSES