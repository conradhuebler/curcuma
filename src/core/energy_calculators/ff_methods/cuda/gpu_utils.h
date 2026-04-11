/*
 * <GPU Memory Utilities>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#ifdef USE_CUDA

#include <cstddef>
#include <string>

namespace GPUUtils {

/**
 * @brief Get available GPU memory information.
 * @param free  Output: free memory in bytes
 * @param total Output: total memory in bytes
 * @return true if successful, false if CUDA error
 *
 * Claude Generated (March 2026): GPU memory management for Phase 2.
 */
bool getGPUMemoryInfo(size_t& free, size_t& total);

/**
 * @brief Estimate required GPU memory for a GFN-FF calculation.
 * @param natoms Number of atoms
 * @return Estimated memory in bytes
 *
 * Memory estimate includes:
 * - SoA coordinate buffers (N * 3 doubles for coords + gradient)
 * - Charge and CN buffers (N doubles each)
 * - EEQ matrix (N * N doubles during solve)
 * - Pair lists (~100 * N integers for dispersion, Coulomb, HB, XB)
 * - Temporary buffers (~50 * N doubles)
 */
size_t estimateGFNFFGPUMemory(int natoms);

/**
 * @brief Check if molecule fits in GPU memory.
 * @param natoms Number of atoms
 * @param threshold Fraction of free memory to use (default 0.8)
 * @return true if sufficient memory, false if insufficient or CUDA error
 *
 * Claude Generated (March 2026): Pre-allocation check to avoid OOM.
 */
bool checkGPUMemoryAvailable(int natoms, double threshold = 0.8);

/**
 * @brief Get human-readable memory size string.
 * @param bytes Memory in bytes
 * @return Formatted string (e.g., "4.2 GB", "256 MB")
 */
std::string formatMemorySize(size_t bytes);

/**
 * @brief Log GPU memory status via CurcumaLogger.
 * @param verbosity Minimum verbosity level to log at
 *
 * Logs free/total GPU memory at verbosity >= 2.
 */
void logGPUMemoryStatus(int verbosity = 2);

} // namespace GPUUtils

#endif // USE_CUDA