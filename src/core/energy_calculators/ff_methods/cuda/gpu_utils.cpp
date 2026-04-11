/*
 * <GPU Memory Utilities Implementation>
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

#ifdef USE_CUDA

#include "gpu_utils.h"
#include "src/core/curcuma_logger.h"

#include <cuda_runtime.h>
#include <cmath>
#include <sstream>
#include <iomanip>

namespace GPUUtils {

bool getGPUMemoryInfo(size_t& free, size_t& total) {
    cudaError_t err = cudaMemGetInfo(&free, &total);
    if (err != cudaSuccess) {
        CurcumaLogger::warn(std::string("cudaMemGetInfo failed: ") + cudaGetErrorString(err));
        return false;
    }
    return true;
}

size_t estimateGFNFFGPUMemory(int natoms) {
    // SoA buffers (see ff_workspace_gpu.cu):
    // - coordinates: N * 3 doubles
    // - gradient: N * 3 doubles
    // - charges: N doubles
    // - CN: N doubles
    // - dEdcn: N doubles
    // - EEQ matrix: N * N doubles (only during solve, but allocate once)
    // - Pair lists: ~100 * N integers (dispersion, Coulomb, HB, XB pairs)
    // - Temporary buffers: ~50 * N doubles

    // Base overhead for CUDA context, kernels, etc.
    size_t base_overhead = 2 * 1024 * 1024;  // 2 MB for kernel overhead

    // Per-atom memory (doubles)
    size_t per_atom_doubles = 3 + 3 + 1 + 1 + 1 + 50;  // coords, gradient, charges, CN, dEdcn, temps

    // Per-atom memory (integers for pair lists)
    size_t per_atom_ints = 100;  // conservative estimate for all pair types

    size_t memory = base_overhead;
    memory += natoms * per_atom_doubles * sizeof(double);
    memory += natoms * natoms * sizeof(double);  // EEQ matrix (worst case)
    memory += natoms * per_atom_ints * sizeof(int);

    // Add 20% safety margin
    memory = static_cast<size_t>(memory * 1.2);

    return memory;
}

bool checkGPUMemoryAvailable(int natoms, double threshold) {
    size_t free, total;
    if (!getGPUMemoryInfo(free, total)) {
        return false;  // CUDA not available or error
    }

    size_t required = estimateGFNFFGPUMemory(natoms);
    size_t available = static_cast<size_t>(free * threshold);

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("=== GPU Memory Check ===");
        CurcumaLogger::param("Required", formatMemorySize(required));
        CurcumaLogger::param("Available", formatMemorySize(available));
        CurcumaLogger::param("Free", formatMemorySize(free));
        CurcumaLogger::param("Total", formatMemorySize(total));
    }

    return required < available;
}

std::string formatMemorySize(size_t bytes) {
    const char* units[] = {"B", "KB", "MB", "GB", "TB"};
    int unit_idx = 0;
    double size = static_cast<double>(bytes);

    while (size >= 1024.0 && unit_idx < 4) {
        size /= 1024.0;
        unit_idx++;
    }

    std::ostringstream oss;
    oss << std::fixed << std::setprecision(1) << size << " " << units[unit_idx];
    return oss.str();
}

void logGPUMemoryStatus(int verbosity) {
    if (CurcumaLogger::get_verbosity() < verbosity) {
        return;
    }

    size_t free, total;
    if (!getGPUMemoryInfo(free, total)) {
        CurcumaLogger::warn("GPU not available for memory status query");
        return;
    }

    size_t used = total - free;

    CurcumaLogger::info("=== GPU Memory Status ===");
    CurcumaLogger::param("Total", formatMemorySize(total));
    CurcumaLogger::param("Used", formatMemorySize(used));
    CurcumaLogger::param("Free", formatMemorySize(free));

    std::ostringstream oss;
    oss << std::fixed << std::setprecision(1) << (100.0 * used / total) << "%";
    CurcumaLogger::param("Utilization", oss.str());
}

} // namespace GPUUtils

#endif // USE_CUDA