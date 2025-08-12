/*
 * Main Precompiled Header for Curcuma Project
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This is the main PCH entry point that includes modular PCH components
 * based on compilation options. This allows for flexible precompilation
 * of different library combinations while avoiding redefinition errors.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 */

#pragma once

// Always include base PCH (standard library + core external deps)
#include "pch/pch_base.h"

// Include safe Curcuma core headers
#include "pch/pch_core.h"

// Conditionally include external library PCH modules
#include "pch/pch_external.h"

// NOTE: Ulysses not included in PCH due to symbol conflicts
// Include Ulysses headers directly in source files where needed

// NOTE: Headers with inline functions (like general.h, geometry.h)
// are NOT precompiled to avoid redefinition errors. These should be
// included directly in source files where needed.