/*
 * < GPU device pool — implementation >
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

// Claude Generated (Sep 2026, multi-GPU case A). See gpu_device_pool.h.

#include "gpu_device_pool.h"

#include "src/core/intra_parallel_context.h"

#include "src/core/curcuma_logger.h"
#include "src/core/energy_calculators/gpu_plugin.h"

#include <fmt/format.h>

#include <algorithm>
#include <cctype>
#include <sstream>
#include <thread>

namespace curcuma {

int& leasedGpuDevice()
{
    static thread_local int device = -1;
    return device;
}

GpuDevicePool& GpuDevicePool::instance()
{
    static GpuDevicePool pool;
    return pool;
}

namespace {

// "0,2" / "2" / 2.0 / [0,2] -> indices; "" / "all" -> empty (= all visible).
std::vector<int> parseDeviceList(const nlohmann::json& v)
{
    std::vector<int> out;
    if (v.is_number()) {
        out.push_back(static_cast<int>(v.get<double>()));
    } else if (v.is_array()) {
        for (const auto& e : v)
            if (e.is_number()) out.push_back(static_cast<int>(e.get<double>()));
    } else if (v.is_string()) {
        std::string s = v.get<std::string>();
        if (s == "all") return out;
        std::stringstream ss(s);
        for (std::string tok; std::getline(ss, tok, ',');) {
            tok.erase(std::remove_if(tok.begin(), tok.end(), ::isspace), tok.end());
            if (!tok.empty()) out.push_back(std::stoi(tok));
        }
    }
    return out;
}

} // namespace

void GpuDevicePool::configure(const nlohmann::json& controller)
{
    std::lock_guard<std::mutex> lock(m_mutex);
    m_backend = "none";
    m_devices.clear();
    m_load.clear();

    std::string gpu = controller.value("gpu", std::string("none"));
    std::transform(gpu.begin(), gpu.end(), gpu.begin(), ::tolower);
    if (gpu.empty() || gpu == "none" || gpu == "cpu")
        return;
    if (gpu == "auto")
        gpu = gpu_plugin::firstAvailable();
    if (gpu == "none" || !gpu_plugin::available(gpu))
        return;   // the method factory warns about the missing plugin

    const int visible = gpu_plugin::deviceCount(gpu);
    if (visible <= 0)
        return;

    std::vector<int> wanted;
    try {
        if (controller.contains("gpu_devices"))
            wanted = parseDeviceList(controller["gpu_devices"]);
    } catch (const std::exception& e) {
        CurcumaLogger::warn(std::string("gpu_devices: cannot parse (") + e.what() + "); using all devices");
        wanted.clear();
    }
    if (wanted.empty()) {
        for (int i = 0; i < visible; ++i) wanted.push_back(i);
    }
    for (int d : wanted) {
        if (d < 0 || d >= visible) {
            // Claude Generated (Sep 2026): "-gpu_devices 2" reads as "two GPUs" but the option is
            // a list of device INDICES, which is the mistake this warning has to name.
            CurcumaLogger::warn(fmt::format(
                "gpu_devices: device {} ignored ({} visible, indices 0..{}). This option takes "
                "device INDICES (e.g. -gpu_devices 0,1) or 'all', not a device COUNT.",
                d, visible, visible - 1));
            continue;
        }
        if (std::find(m_devices.begin(), m_devices.end(), d) == m_devices.end())
            m_devices.push_back(d);
    }
    if (m_devices.empty())
        return;

    // A single visible device without explicit pool flags keeps the historical behaviour:
    // no lease, every worker shares device 0 concurrently. The pool only takes over when
    // there is something to distribute or the user asked for it.
    const bool explicit_pool = controller.contains("gpu_devices") || controller.contains("gpu_workers_per_device");
    if (m_devices.size() == 1 && !explicit_pool) {
        m_devices.clear();
        return;
    }

    const double per = controller.value("gpu_workers_per_device", 1.0);
    m_per_device = std::max(1, static_cast<int>(per));
    m_load.assign(m_devices.size(), 0);
    m_backend = gpu;

    if (m_devices.size() > 1 || m_per_device > 1) {
        std::string list;
        for (int d : m_devices) list += (list.empty() ? "" : ",") + std::to_string(d);
        CurcumaLogger::info(fmt::format("GPU batch pool: backend {}, devices [{}], {} worker(s) per device",
                                        m_backend, list, m_per_device));
    }
}

bool GpuDevicePool::active() const
{
    std::lock_guard<std::mutex> lock(m_mutex);
    return !m_devices.empty();
}

std::vector<int> GpuDevicePool::devices() const
{
    std::lock_guard<std::mutex> lock(m_mutex);
    return m_devices;
}

int GpuDevicePool::capacity() const
{
    std::lock_guard<std::mutex> lock(m_mutex);
    return static_cast<int>(m_devices.size()) * m_per_device;
}

std::string GpuDevicePool::backend() const
{
    std::lock_guard<std::mutex> lock(m_mutex);
    return m_backend;
}

int GpuDevicePool::acquire()
{
    std::unique_lock<std::mutex> lock(m_mutex);
    if (m_devices.empty())
        return -1;
    // Least-loaded device with a free slot; ties go to the lower list position, so a
    // two-worker run on "0,1,2,3" uses devices 0 and 1 (often a PCIe/NVLink-close pair).
    auto pick = [&]() -> int {
        int best = -1;
        for (size_t i = 0; i < m_devices.size(); ++i)
            if (m_load[i] < m_per_device && (best < 0 || m_load[i] < m_load[best]))
                best = static_cast<int>(i);
        return best;
    };
    int slot = pick();
    while (slot < 0) {
        m_cv.wait(lock);
        slot = pick();
    }
    ++m_load[slot];
    return m_devices[slot];
}

void GpuDevicePool::release(int device)
{
    {
        std::lock_guard<std::mutex> lock(m_mutex);
        for (size_t i = 0; i < m_devices.size(); ++i) {
            if (m_devices[i] == device && m_load[i] > 0) {
                --m_load[i];
                break;
            }
        }
    }
    m_cv.notify_one();
}

GpuDeviceLease::GpuDeviceLease()
    : m_prev(leasedGpuDevice())
{
    // Nested lease on a thread that already holds one (e.g. an optimisation launched from
    // inside a batch task): reuse it. Acquiring again would take a second slot for the same
    // work, or deadlock once every slot is held.
    if (m_prev >= 0)
        return;
    m_device = GpuDevicePool::instance().acquire();
    if (m_device >= 0) {
        leasedGpuDevice() = m_device;
        // Host-thread share for this GPU worker: the device does the heavy lifting, so split
        // the cores evenly over the GPU slots instead of running the host stages serially.
        m_prev_budget = intraThreadBudget();
        const unsigned hc = std::thread::hardware_concurrency();
        const int slots = std::max(1, GpuDevicePool::instance().capacity());
        intraThreadBudget() = std::max(1, static_cast<int>(hc > 0 ? hc : 1) / slots);
    }
}

GpuDeviceLease::~GpuDeviceLease()
{
    if (m_device >= 0) {
        GpuDevicePool::instance().release(m_device);
        leasedGpuDevice() = m_prev;
        intraThreadBudget() = m_prev_budget;
    }
}

} // namespace curcuma
