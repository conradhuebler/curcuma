/*
 * Curcuma Core Precompiled Header
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This header includes safe Curcuma core headers that do not contain
 * inline functions that would cause redefinition errors.
 */

#pragma once

// Include base PCH first
#include "pch_base.h"

// Safe core headers (no inline functions, only declarations/definitions)
#include "src/core/elements.h" // Element constants and enums
#include "src/core/topology.h" // Topology structures

// Conditional parameter includes based on Ulysses usage
// Include parameter headers only if Ulysses is not used to avoid conflicts
#include "src/core/energy_calculators/ff_methods/qmdff_par.h" // QMDFF parameter definitions
#include "src/core/energy_calculators/ff_methods/uff_par.h" // UFF parameter definitions

// NOTE: When Ulysses is enabled, parameter headers must be included
// directly in source files due to symbol conflicts (kEN, etc.)

// NOTE: The following headers are NOT included due to inline functions:
// - src/tools/general.h    (has many inline utility functions)
// - src/tools/geometry.h   (has inline geometric calculations)
// - src/core/molecule.h    (may have inline methods)
// - src/core/energycalculator.h (complex header with dependencies)

// These headers should be included directly in source files where needed.