/*
 * < GPU backend plugin loader — implementation >
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

// Claude Generated (Jul 2026): see gpu_plugin.h for the rationale (fast CPU startup).

#include "gpu_plugin.h"

#include "src/core/curcuma_logger.h"

#include <dlfcn.h>
#include <limits.h>
#include <unistd.h>

#include <map>
#include <mutex>

namespace gpu_plugin {
namespace {

std::mutex g_mutex;
// dlopen handle cache per plugin .so. FAILED_SENTINEL marks a prior failed attempt so we
// do not re-probe (and re-log) every call.
std::map<std::string, void*> g_handles;
void* const FAILED_SENTINEL = reinterpret_cast<void*>(1);

// Load libcurcuma_<backend>.so, preferring the copy next to the curcuma executable (so no
// LD_LIBRARY_PATH is required), then the default loader search path. RTLD_GLOBAL so the
// plugin's undefined core symbols bind to the already-loaded executable (built -rdynamic).
void* pluginHandle(const std::string& backend)
{
    std::lock_guard<std::mutex> lock(g_mutex);
    const std::string soname = "libcurcuma_" + backend + ".so";
    auto it = g_handles.find(soname);
    if (it != g_handles.end())
        return it->second == FAILED_SENTINEL ? nullptr : it->second;

    void* handle = nullptr;
    char buf[PATH_MAX];
    const ssize_t n = ::readlink("/proc/self/exe", buf, sizeof(buf) - 1);
    if (n > 0) {
        buf[n] = '\0';
        std::string path(buf);
        const auto slash = path.find_last_of('/');
        if (slash != std::string::npos)
            handle = ::dlopen((path.substr(0, slash + 1) + soname).c_str(),
                              RTLD_NOW | RTLD_GLOBAL);
    }
    if (!handle)
        handle = ::dlopen(soname.c_str(), RTLD_NOW | RTLD_GLOBAL);

    if (!handle) {
        const char* err = ::dlerror();
        CurcumaLogger::warn("GPU plugin '" + soname + "' could not be loaded ("
                            + (err ? err : "unknown") + "); falling back to CPU.");
    }
    g_handles[soname] = handle ? handle : FAILED_SENTINEL;
    return handle;
}

template <typename FnPtr>
FnPtr resolveSymbol(void* handle, const std::string& name)
{
    ::dlerror();  // clear
    return reinterpret_cast<FnPtr>(::dlsym(handle, name.c_str()));
}

} // namespace

std::unique_ptr<ComputationalMethod> createNativeXtb(const std::string& backend,
                                                     int method_type, const json& config)
{
    void* handle = pluginHandle(backend);
    if (!handle)
        return nullptr;
    using create_fn = ComputationalMethod* (*)(int, const char*);
    const std::string sym = "curcuma_" + backend + "_create_native_xtb";
    create_fn fn = resolveSymbol<create_fn>(handle, sym);
    if (!fn) {
        CurcumaLogger::warn("GPU plugin missing entry point '" + sym + "'; falling back to CPU.");
        return nullptr;
    }
    const std::string cfg = config.dump();
    ComputationalMethod* method = fn(method_type, cfg.c_str());
    return std::unique_ptr<ComputationalMethod>(method);
}

std::unique_ptr<ComputationalMethod> createGfnff(const std::string& backend, const json& config)
{
    void* handle = pluginHandle(backend);
    if (!handle)
        return nullptr;
    using create_fn = ComputationalMethod* (*)(const char*);
    const std::string sym = "curcuma_" + backend + "_create_gfnff";
    create_fn fn = resolveSymbol<create_fn>(handle, sym);
    if (!fn) {
        CurcumaLogger::warn("GPU plugin missing entry point '" + sym + "'; falling back to CPU.");
        return nullptr;
    }
    const std::string cfg = config.dump();
    ComputationalMethod* method = fn(cfg.c_str());
    return std::unique_ptr<ComputationalMethod>(method);
}

} // namespace gpu_plugin
