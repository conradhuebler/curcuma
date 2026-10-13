/*
 * < GPU device pool — distribute molecule-level batch workers over several GPUs >
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

// Claude Generated (Sep 2026, multi-GPU case A): "many structures, several GPUs".
//
// Curcuma already parallelises batch work at the MOLECULE level: the batch optimiser,
// CurcumaOpt::ProcessMolecules, ConfSearch and the numerical Hessian run one independent
// energy calculation per CxxThreadPool worker, and every worker builds its own
// EnergyCalculator (hence its own GPU context). Without this file all of those contexts
// land on device 0.
//
// The mechanism mirrors intra_parallel_context.h:
//   1. main() configures ONE process-wide pool from the CLI (`-gpu`, `-gpu_devices`,
//      `-gpu_workers_per_device`).
//   2. A batch worker constructs a `GpuDeviceLease` at the top of execute(). The lease
//      blocks until a device slot is free, takes the least-loaded device and publishes it
//      in a thread-local variable.
//   3. EnergyCalculator::createMethod() copies that thread-local device into the method
//      config as `gpu_device`, so the GPU plugin binds the context to it.
// Single calculations (no lease) are unaffected, and so is every run without `-gpu`.

#pragma once

#include <condition_variable>
#include <mutex>
#include <string>
#include <vector>

#include "json.hpp"

namespace curcuma {

/// Device assigned to the calling thread by a GpuDeviceLease, or -1 (no batch lease).
int& leasedGpuDevice();

/**
 * @brief Process-wide pool of GPU device slots for molecule-level batch workers.
 *
 * Capacity per device is `gpu_workers_per_device`; acquire() blocks when every slot is in
 * use, so a thread pool larger than the slot count simply queues on the GPUs instead of
 * overcommitting device memory.
 */
class GpuDevicePool {
public:
    static GpuDevicePool& instance();

    /**
     * @brief Configure from the top-level controller.
     *
     * Reads `gpu` (none|cpu|auto|cuda|rocm|vulkan), `gpu_devices` ("" or "all" = every
     * visible device, else a comma list like "0,2") and `gpu_workers_per_device` (>= 1).
     * Probing the device count loads the backend plugin, so this does nothing for CPU runs.
     */
    void configure(const nlohmann::json& controller);

    /// True when a GPU backend was requested and at least one device is usable.
    bool active() const;
    /// The device indices in use (runtime numbering, after CUDA_VISIBLE_DEVICES).
    std::vector<int> devices() const;
    /// Total slots = devices x workers per device (0 when inactive).
    int capacity() const;
    std::string backend() const;

    /// Block until a slot is free; return the least-loaded device index. -1 when inactive.
    int acquire();
    void release(int device);

private:
    GpuDevicePool() = default;

    mutable std::mutex m_mutex;
    std::condition_variable m_cv;
    std::string m_backend = "none";
    std::vector<int> m_devices;
    std::vector<int> m_load;   // active leases per entry of m_devices
    int m_per_device = 1;
};

/**
 * @brief RAII device lease for one batch task.
 *
 * Construct at the top of a molecule-level worker's execute(). No-op when the pool is
 * inactive. Restores the previous thread-local device on destruction (nesting-safe).
 */
class GpuDeviceLease {
public:
    GpuDeviceLease();
    ~GpuDeviceLease();
    GpuDeviceLease(const GpuDeviceLease&) = delete;
    GpuDeviceLease& operator=(const GpuDeviceLease&) = delete;

    /// Leased device, or -1 when no GPU pool is active.
    int device() const { return m_device >= 0 ? m_device : m_prev; }

private:
    int m_device = -1;
    int m_prev = -1;
    int m_prev_budget = 1;
};

} // namespace curcuma
