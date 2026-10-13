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

// Claude Generated (Sep 2026). See gpu_fallback.h.

#include "gpu_fallback.h"

#include "curcuma_logger.h"

#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <mutex>

#include <fmt/format.h>

namespace curcuma {

namespace {

struct Entry {
    std::string category;
    std::string first_detail;
    int count = 0;
};

std::mutex& registryMutex()
{
    static std::mutex m;
    return m;
}

std::vector<Entry>& registry()
{
    static std::vector<Entry> r;
    return r;
}

std::atomic<bool>& strictFlag()
{
    static std::atomic<bool> s{ false };
    return s;
}

} // namespace

void reportGpuFallback(const std::string& category, const std::string& detail)
{
    {
        std::lock_guard<std::mutex> lock(registryMutex());
        auto& r = registry();
        bool found = false;
        for (auto& e : r) {
            if (e.category == category) {
                ++e.count;
                found = true;
                break;
            }
        }
        if (!found)
            r.push_back({ category, detail, 1 });
    }

    if (strictFlag().load()) {
        CurcumaLogger::error(fmt::format("-gpu_strict: GPU fallback '{}'{} - stopping (exit code {})",
            category, detail.empty() ? std::string() : " (" + detail + ")", kGpuStrictExitCode));
        // _Exit, not exit: the report can come from a batch worker thread while other workers
        // hold GPU contexts; running static destructors under them can hang the process.
        std::fflush(nullptr);
        std::_Exit(kGpuStrictExitCode);
    }
}

void setGpuStrict(bool strict) { strictFlag().store(strict); }

bool gpuStrict() { return strictFlag().load(); }

std::vector<std::pair<std::string, int>> gpuFallbackCounts()
{
    std::lock_guard<std::mutex> lock(registryMutex());
    std::vector<std::pair<std::string, int>> out;
    for (const auto& e : registry())
        out.emplace_back(e.category, e.count);
    return out;
}

void printGpuFallbackSummary()
{
    std::vector<Entry> copy;
    {
        std::lock_guard<std::mutex> lock(registryMutex());
        copy = registry();
    }
    if (copy.empty())
        return;
    // warn() is suppressed only at verbosity 0 of the main thread; the summary is the one place
    // a silent run reports its degradations, so it is raised to level 1 for these lines.
    const int saved = CurcumaLogger::get_verbosity();
    if (saved < 1)
        CurcumaLogger::set_verbosity(1);
    CurcumaLogger::warn("GPU fallbacks during this run (use -gpu_strict true to stop at the first one):");
    for (const auto& e : copy)
        CurcumaLogger::warn(fmt::format("  {} x {}{}", e.count, e.category,
            e.first_detail.empty() ? std::string() : " - first: " + e.first_detail));
    if (saved < 1)
        CurcumaLogger::set_verbosity(saved);
}

} // namespace curcuma
