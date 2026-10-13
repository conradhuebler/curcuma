/*
 * < GPU fallback registry - make every GPU degradation visible >
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
 */

// Claude Generated (Sep 2026, docs/MULTI_GPU_GAPS.md G2-7/G2-11/F-19/G2-15).
//
// A GPU run can quietly degrade: a batch worker falls back to the CPU, the multi-GPU
// eigensolver picks a slower library, the device EEQ keeps the previous charges. The
// individual warnings are gated by verbosity, and batch workers run at verbosity 0, so
// such a run looks normal and is only slow (or, for the EEQ case, slightly different).
//
// Every such site calls reportGpuFallback() with a short, stable category string. The
// registry counts per category (thread-safe) and main() prints one summary at the end at
// warning level, independent of the verbosity the workers ran at. With `-gpu_strict true`
// the first report ends the process with exit code 3 instead.
//
// The GPU plugins resolve these symbols from the main executable (they already use
// CurcumaLogger and the device pool the same way).

#pragma once

#include <string>
#include <utility>
#include <vector>

namespace curcuma {

/// Exit code of a run ended by -gpu_strict.
constexpr int kGpuStrictExitCode = 3;

/**
 * @brief Record one GPU degradation.
 * @param category short stable text, e.g. "GFN-FF EEQ: no valid device solution"; the
 *                 summary counts per category
 * @param detail   optional context for the first occurrence (device, reason)
 *
 * Under -gpu_strict this prints an error and terminates the process (exit code 3).
 */
void reportGpuFallback(const std::string& category, const std::string& detail = std::string());

/// Enable/disable strict mode (`-gpu_strict`).
void setGpuStrict(bool strict);
bool gpuStrict();

/// Categories with their counts, in first-seen order.
std::vector<std::pair<std::string, int>> gpuFallbackCounts();

/// Print the summary as warnings (nothing when no fallback happened).
void printGpuFallbackSummary();

} // namespace curcuma
