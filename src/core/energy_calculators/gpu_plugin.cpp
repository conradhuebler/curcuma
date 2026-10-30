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

// Claude Generated (Aug 2026): USE_CUDA/ROCM/VULKAN are never enabled on the Windows CI
// runner (no GPU toolkit there), so no plugin is ever actually built or loaded on Windows —
// but this file is still part of curcuma_core unconditionally, so it must compile there.
// dlfcn.h/unistd.h are POSIX-only; use the LoadLibrary/GetProcAddress equivalents on Windows.
#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#else
#include <dlfcn.h>
#include <limits.h>
#include <unistd.h>
#endif

#include <map>
#include <mutex>

namespace gpu_plugin {
namespace {

std::mutex g_mutex;
// dlopen/LoadLibrary handle cache per plugin library. FAILED_SENTINEL marks a prior failed
// attempt so we do not re-probe (and re-log) every call.
std::map<std::string, void*> g_handles;
void* const FAILED_SENTINEL = reinterpret_cast<void*>(1);

#ifdef _WIN32
constexpr const char* kPluginExt = ".dll";
#else
constexpr const char* kPluginExt = ".so";
#endif

// Load libcurcuma_<backend>.{so,dll}, preferring the copy next to the curcuma executable (so
// no LD_LIBRARY_PATH/PATH entry is required), then the default loader search path. On POSIX,
// RTLD_GLOBAL so the plugin's undefined core symbols bind to the already-loaded executable
// (built -rdynamic).
void* pluginHandle(const std::string& backend, bool quiet = false)
{
    std::lock_guard<std::mutex> lock(g_mutex);
    const std::string soname = "libcurcuma_" + backend + kPluginExt;
    auto it = g_handles.find(soname);
    if (it != g_handles.end())
        return it->second == FAILED_SENTINEL ? nullptr : it->second;

    void* handle = nullptr;
#ifdef _WIN32
    char buf[MAX_PATH];
    const DWORD n = ::GetModuleFileNameA(nullptr, buf, sizeof(buf));
    if (n > 0 && n < sizeof(buf)) {
        std::string path(buf, n);
        const auto slash = path.find_last_of("\\/");
        if (slash != std::string::npos)
            handle = reinterpret_cast<void*>(::LoadLibraryA((path.substr(0, slash + 1) + soname).c_str()));
    }
    if (!handle)
        handle = reinterpret_cast<void*>(::LoadLibraryA(soname.c_str()));

    if (!handle) {
        if (!quiet) CurcumaLogger::warn("GPU plugin '" + soname + "' could not be loaded (error "
                            + std::to_string(::GetLastError()) + "); falling back to CPU.");
    }
#else
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
        if (!quiet) CurcumaLogger::warn("GPU plugin '" + soname + "' could not be loaded ("
                            + (err ? err : "unknown") + "); falling back to CPU.");
    }
#endif
    g_handles[soname] = handle ? handle : FAILED_SENTINEL;
    return handle;
}

template <typename FnPtr>
FnPtr resolveSymbol(void* handle, const std::string& name)
{
#ifdef _WIN32
    return reinterpret_cast<FnPtr>(::GetProcAddress(reinterpret_cast<HMODULE>(handle), name.c_str()));
#else
    ::dlerror();  // clear
    return reinterpret_cast<FnPtr>(::dlsym(handle, name.c_str()));
#endif
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

std::unique_ptr<qmint::QMGpuJK> createQmJK(const std::string& backend, const double* eri, int n, int device)
{
    void* handle = pluginHandle(backend, /*quiet=*/true);
    if (!handle)
        return nullptr;
    using create_fn = void* (*)(const double*, int, int);
    create_fn fn = resolveSymbol<create_fn>(handle, "curcuma_" + backend + "_create_qm_jk");
    if (!fn)
        return nullptr;
    return std::unique_ptr<qmint::QMGpuJK>(static_cast<qmint::QMGpuJK*>(fn(eri, n, device)));
}

// Claude Generated (Oct 2026, Q-G2): ERI tensor built on the device.
std::unique_ptr<qmint::QMGpuJK> createQmEriJK(const std::string& backend, const qmint::QMGpuEriTables& tables,
                                              int device, std::string& why)
{
    void* handle = pluginHandle(backend, /*quiet=*/true);
    if (!handle) {
        why = "plugin libcurcuma_" + backend + ".so not loadable";
        return nullptr;
    }
    using create_fn = void* (*)(const qmint::QMGpuEriTables*, int, char*, int);
    create_fn fn = resolveSymbol<create_fn>(handle, "curcuma_" + backend + "_create_qm_eri_jk");
    if (!fn) {
        why = "plugin has no entry point curcuma_" + backend + "_create_qm_eri_jk";
        return nullptr;
    }
    char buf[512] = { 0 };
    void* p = fn(&tables, device, buf, (int)sizeof buf);
    if (!p) why = buf;
    return std::unique_ptr<qmint::QMGpuJK>(static_cast<qmint::QMGpuJK*>(p));
}

// Claude Generated (Oct 2026, Q-G3): integral-direct J/K on the device.
std::unique_ptr<qmint::QMGpuJK> createQmDirectJK(const std::string& backend, const qmint::QMGpuEriTables& tables,
                                                 int device, std::string& why)
{
    void* handle = pluginHandle(backend, /*quiet=*/true);
    if (!handle) {
        why = "plugin libcurcuma_" + backend + ".so not loadable";
        return nullptr;
    }
    using create_fn = void* (*)(const qmint::QMGpuEriTables*, int, char*, int);
    create_fn fn = resolveSymbol<create_fn>(handle, "curcuma_" + backend + "_create_qm_direct_jk");
    if (!fn) {
        why = "plugin has no entry point curcuma_" + backend + "_create_qm_direct_jk";
        return nullptr;
    }
    char buf[512] = { 0 };
    void* p = fn(&tables, device, buf, (int)sizeof buf);
    if (!p) why = buf;
    return std::unique_ptr<qmint::QMGpuJK>(static_cast<qmint::QMGpuJK*>(p));
}

// Claude Generated (Oct 2026, Q-G4): device-resident SCF state.
std::unique_ptr<qmint::QMGpuSCF> createQmScf(const std::string& backend, int n, int device, std::string& why)
{
    void* handle = pluginHandle(backend, /*quiet=*/true);
    if (!handle) {
        why = "plugin libcurcuma_" + backend + ".so not loadable";
        return nullptr;
    }
    using create_fn = void* (*)(int, int, char*, int);
    // _ks (Oct 2026, GPU KS-DFT): QMGpuSCF gained setKohnSham(); the entry was renamed with the vtable, so an older
    // plugin is refused (counted host-loop fallback) instead of being called through a shorter vtable.
    create_fn fn = resolveSymbol<create_fn>(handle, "curcuma_" + backend + "_create_qm_scf_ks");
    if (!fn) {
        why = "plugin has no entry point curcuma_" + backend + "_create_qm_scf_ks";
        return nullptr;
    }
    char buf[512] = { 0 };
    void* p = fn(n, device, buf, (int)sizeof buf);
    if (!p) why = buf;
    return std::unique_ptr<qmint::QMGpuSCF>(static_cast<qmint::QMGpuSCF*>(p));
}

// Claude Generated (Oct 2026, Q-G7): two-electron gradient on the device.
std::unique_ptr<qmint::QMGpuGrad2e> createQmGrad2e(const std::string& backend, const qmint::QMGpuEriTables& tables,
                                                   const qmint::QMGpuGradExtra& extra, int device, std::string& why)
{
    void* handle = pluginHandle(backend, /*quiet=*/true);
    if (!handle) {
        why = "plugin libcurcuma_" + backend + ".so not loadable";
        return nullptr;
    }
    using create_fn = void* (*)(const qmint::QMGpuEriTables*, const qmint::QMGpuGradExtra*, int, char*, int);
    // _kx (Oct 2026, GPU KS-DFT): compute() gained the exchange weight kx; renamed with the vtable (see createQmScf).
    create_fn fn = resolveSymbol<create_fn>(handle, "curcuma_" + backend + "_create_qm_grad2e_kx");
    if (!fn) {
        why = "plugin has no entry point curcuma_" + backend + "_create_qm_grad2e_kx";
        return nullptr;
    }
    char buf[512] = { 0 };
    void* p = fn(&tables, &extra, device, buf, (int)sizeof buf);
    if (!p) why = buf;
    return std::unique_ptr<qmint::QMGpuGrad2e>(static_cast<qmint::QMGpuGrad2e*>(p));
}

// Claude Generated (Oct 2026, GPU KS-DFT): Kohn-Sham exchange-correlation on the device.
std::unique_ptr<qmint::QMGpuXC> createQmXC(const std::string& backend, const qmint::QMGpuXCSetup& setup, int device,
                                           std::string& why)
{
    void* handle = pluginHandle(backend, /*quiet=*/true);
    if (!handle) {
        why = "plugin libcurcuma_" + backend + ".so not loadable";
        return nullptr;
    }
    using create_fn = void* (*)(const qmint::QMGpuXCSetup*, int, char*, int);
    create_fn fn = resolveSymbol<create_fn>(handle, "curcuma_" + backend + "_create_qm_xc_screen");
    if (!fn) {
        why = "plugin has no entry point curcuma_" + backend + "_create_qm_xc_screen";
        return nullptr;
    }
    char buf[512] = { 0 };
    void* p = fn(&setup, device, buf, (int)sizeof buf);
    if (!p) why = buf;
    return std::unique_ptr<qmint::QMGpuXC>(static_cast<qmint::QMGpuXC*>(p));
}

// Claude Generated (Oct 2026, MP2 gradient): 2e derivative contraction of the MP2 gradient on the device.
std::unique_ptr<qmint::QMGpuMP2Grad2e> createQmMp2Grad2e(const std::string& backend, const qmint::QMGpuEriTables& tables,
                                                         const qmint::QMGpuGradExtra& extra, int device, std::string& why)
{
    void* handle = pluginHandle(backend, /*quiet=*/true);
    if (!handle) {
        why = "plugin libcurcuma_" + backend + ".so not loadable";
        return nullptr;
    }
    using create_fn = void* (*)(const qmint::QMGpuEriTables*, const qmint::QMGpuGradExtra*, int, char*, int);
    create_fn fn = resolveSymbol<create_fn>(handle, "curcuma_" + backend + "_create_qm_mp2_grad2e");
    if (!fn) {
        why = "plugin has no entry point curcuma_" + backend + "_create_qm_mp2_grad2e";
        return nullptr;
    }
    char buf[512] = { 0 };
    void* p = fn(&tables, &extra, device, buf, (int)sizeof buf);
    if (!p) why = buf;
    return std::unique_ptr<qmint::QMGpuMP2Grad2e>(static_cast<qmint::QMGpuMP2Grad2e*>(p));
}

// Claude Generated (Oct 2026, Q-G8): one-electron gradient on the device.
std::unique_ptr<qmint::QMGpuGrad1e> createQmGrad1e(const std::string& backend, const qmint::QMGpuEriTables& tables,
                                                   const qmint::QMGpuGrad1eExtra& extra, int device, std::string& why)
{
    void* handle = pluginHandle(backend, /*quiet=*/true);
    if (!handle) {
        why = "plugin libcurcuma_" + backend + ".so not loadable";
        return nullptr;
    }
    using create_fn = void* (*)(const qmint::QMGpuEriTables*, const qmint::QMGpuGrad1eExtra*, int, char*, int);
    create_fn fn = resolveSymbol<create_fn>(handle, "curcuma_" + backend + "_create_qm_grad1e");
    if (!fn) {
        why = "plugin has no entry point curcuma_" + backend + "_create_qm_grad1e";
        return nullptr;
    }
    char buf[512] = { 0 };
    void* p = fn(&tables, &extra, device, buf, (int)sizeof buf);
    if (!p) why = buf;
    return std::unique_ptr<qmint::QMGpuGrad1e>(static_cast<qmint::QMGpuGrad1e*>(p));
}

// Claude Generated (Oct 2026, WP10): density-fitting (RI) contractions on the device.
std::unique_ptr<qmint::QMGpuRI> createQmRI(const std::string& backend, const double* T, int nfit, int n, int device,
                                           std::string& why)
{
    void* handle = pluginHandle(backend, /*quiet=*/true);
    if (!handle) {
        why = "plugin libcurcuma_" + backend + ".so not loadable";
        return nullptr;
    }
    using create_fn = void* (*)(const double*, int, int, int, char*, int);
    create_fn fn = resolveSymbol<create_fn>(handle, "curcuma_" + backend + "_create_qm_ri_packed");
    if (!fn) {
        why = "plugin has no entry point curcuma_" + backend + "_create_qm_ri_packed (packed RI tensor)";
        return nullptr;
    }
    char buf[512] = { 0 };
    void* p = fn(T, nfit, n, device, buf, (int)sizeof buf);
    if (!p) why = buf;
    return std::unique_ptr<qmint::QMGpuRI>(static_cast<qmint::QMGpuRI*>(p));
}

// Claude Generated (Oct 2026, WP10 scaling): RI-JK metric solve T := L^-1 T on the device.
bool qmRiLowerSolve(const std::string& backend, double* T, int N, long ncol, const double* L, int device, int host_threads,
                    bool& touched, std::string& why)
{
    touched = false;
    void* handle = pluginHandle(backend, /*quiet=*/true);
    if (!handle) {
        why = "plugin libcurcuma_" + backend + ".so not loadable";
        return false;
    }
    using fn_t = int (*)(double*, int, long, const double*, int, int, int*, char*, int);
    fn_t fn = resolveSymbol<fn_t>(handle, "curcuma_" + backend + "_qm_ri_lower_solve");
    if (!fn) {
        why = "plugin has no entry point curcuma_" + backend + "_qm_ri_lower_solve";
        return false;
    }
    char buf[512] = { 0 };
    int t = 0;
    const int ok = fn(T, N, ncol, L, device, host_threads, &t, buf, (int)sizeof buf);
    touched = t != 0;
    if (!ok) why = buf;
    return ok != 0;
}

// Claude Generated (Oct 2026, WP10 low memory): resident or streamed RI tensor, FP64 or FP32.
std::unique_ptr<qmint::QMGpuRILowMem> createQmRILowMem(const std::string& backend, const void* T, bool fp32, int nfit, int n,
                                                       int device, int stream, double max_device_mb, double chunk_mb,
                                                       int host_threads, std::string& why)
{
    void* handle = pluginHandle(backend, /*quiet=*/true);
    if (!handle) {
        why = "plugin libcurcuma_" + backend + ".so not loadable";
        return nullptr;
    }
    using create_fn = void* (*)(const void*, int, int, int, int, int, double, double, int, char*, int);
    create_fn fn = resolveSymbol<create_fn>(handle, "curcuma_" + backend + "_create_qm_ri_lowmem");
    if (!fn) {
        why = "plugin has no entry point curcuma_" + backend + "_create_qm_ri_lowmem";
        return nullptr;
    }
    char buf[512] = { 0 };
    void* p = fn(T, fp32 ? 1 : 0, nfit, n, device, stream, max_device_mb, chunk_mb, host_threads, buf, (int)sizeof buf);
    if (!p) why = buf;
    return std::unique_ptr<qmint::QMGpuRILowMem>(static_cast<qmint::QMGpuRILowMem*>(p));
}

// Claude Generated (Oct 2026, WP10 low memory): batched RI-JK exchange X = L^-1 Y, K += sign X X^T on the device.
std::unique_ptr<qmint::QMGpuRIYSolver> createQmRIYSolver(const std::string& backend, const double* L, int naux, int n, int device,
                                                         double max_device_mb, int host_threads, std::string& why)
{
    void* handle = pluginHandle(backend, /*quiet=*/true);
    if (!handle) {
        why = "plugin libcurcuma_" + backend + ".so not loadable";
        return nullptr;
    }
    using create_fn = void* (*)(const double*, int, int, int, double, int, char*, int);
    create_fn fn = resolveSymbol<create_fn>(handle, "curcuma_" + backend + "_create_qm_ri_ysolver");
    if (!fn) {
        why = "plugin has no entry point curcuma_" + backend + "_create_qm_ri_ysolver";
        return nullptr;
    }
    char buf[512] = { 0 };
    void* p = fn(L, naux, n, device, max_device_mb, host_threads, buf, (int)sizeof buf);
    if (!p) why = buf;
    return std::unique_ptr<qmint::QMGpuRIYSolver>(static_cast<qmint::QMGpuRIYSolver*>(p));
}

// Claude Generated (Oct 2026, WP11): canonical MP2 energy on the device.
bool qmMp2Energy(const std::string& backend, const qmint::QMGpuMP2Input& in, qmint::QMGpuMP2Output& out,
                 std::string& why)
{
    void* handle = pluginHandle(backend, /*quiet=*/true);
    if (!handle) {
        why = "plugin libcurcuma_" + backend + ".so not loadable";
        return false;
    }
    using fn_t = int (*)(const qmint::QMGpuMP2Input*, qmint::QMGpuMP2Output*, char*, int);
    fn_t fn = resolveSymbol<fn_t>(handle, "curcuma_" + backend + "_qm_mp2_energy");
    if (!fn) {
        why = "plugin has no entry point curcuma_" + backend + "_qm_mp2_energy";
        return false;
    }
    char buf[512] = { 0 };
    const int ok = fn(&in, &out, buf, (int)sizeof buf);
    if (!ok) why = buf;
    return ok != 0;
}

// Claude Generated (Oct 2026): RI-MP2 pair energies on the device.
bool qmRiMp2Energy(const std::string& backend, const qmint::QMGpuRIMP2Input& in, qmint::QMGpuRIMP2Output& out,
                   std::string& why)
{
    void* handle = pluginHandle(backend, /*quiet=*/true);
    if (!handle) {
        why = "plugin libcurcuma_" + backend + ".so not loadable";
        return false;
    }
    using fn_t = int (*)(const qmint::QMGpuRIMP2Input*, qmint::QMGpuRIMP2Output*, char*, int);
    fn_t fn = resolveSymbol<fn_t>(handle, "curcuma_" + backend + "_qm_rimp2_energy");
    if (!fn) {
        why = "plugin has no entry point curcuma_" + backend + "_qm_rimp2_energy";
        return false;
    }
    char buf[512] = { 0 };
    const int ok = fn(&in, &out, buf, (int)sizeof buf);
    if (!ok) why = buf;
    return ok != 0;
}

// Claude Generated (Oct 2026): 3-index part of the RI gradient on the device.
bool qmRiGrad(const std::string& backend, const qmint::QMGpuRIGradInput& in, double* grad, qmint::QMGpuRIGradOutput& out,
              int device, std::string& why)
{
    void* handle = pluginHandle(backend, /*quiet=*/true);
    if (!handle) {
        why = "plugin libcurcuma_" + backend + ".so not loadable";
        return false;
    }
    using fn_t = int (*)(const qmint::QMGpuRIGradInput*, double*, qmint::QMGpuRIGradOutput*, int, char*, int);
    fn_t fn = resolveSymbol<fn_t>(handle, "curcuma_" + backend + "_qm_ri_grad");
    if (!fn) {
        why = "plugin has no entry point curcuma_" + backend + "_qm_ri_grad";
        return false;
    }
    char buf[512] = { 0 };
    const int ok = fn(&in, grad, &out, device, buf, (int)sizeof buf);
    if (!ok) why = buf;
    return ok != 0;
}

int deviceCount(const std::string& backend)
{
    void* handle = pluginHandle(backend, /*quiet=*/true);
    if (!handle)
        return 0;
    using count_fn = int (*)();
    count_fn fn = resolveSymbol<count_fn>(handle, "curcuma_" + backend + "_device_count");
    return fn ? fn() : 1;   // pre-multi-GPU plugin: one (default) device
}

// Claude Generated (Sep 2026): which distributed-eigensolver backends this plugin can reach
// ("mp", "mg", "mp,mg", "none ..."). Optional symbol - an older plugin returns "" and the
// caller stays silent. The point is that `curcuma -methods` answers "was this build linked
// against cuSOLVERMp, or is it the 15x-slower cusolverMg fallback?" without a calculation.
std::string mgpuBackends(const std::string& backend)
{
    void* handle = pluginHandle(backend, /*quiet=*/true);
    if (!handle)
        return {};
    using str_fn = const char* (*)();
    str_fn fn = resolveSymbol<str_fn>(handle, "curcuma_" + backend + "_mgpu_backends");
    if (!fn)
        return {};
    const char* s = fn();
    return s ? std::string(s) : std::string();
}

json deviceInfo(const std::string& backend, int index)
{
    void* handle = pluginHandle(backend, /*quiet=*/true);
    if (!handle)
        return json::object();
    using info_fn = int (*)(int, char*, int);
    info_fn fn = resolveSymbol<info_fn>(handle, "curcuma_" + backend + "_device_info");
    if (!fn)
        return json::object();
    std::string buf(1024, '\0');
    int need = fn(index, buf.data(), static_cast<int>(buf.size()));
    if (need < 0)
        return json::object();
    if (need >= static_cast<int>(buf.size())) {
        buf.assign(static_cast<size_t>(need) + 1, '\0');
        need = fn(index, buf.data(), static_cast<int>(buf.size()));
        if (need < 0)
            return json::object();
    }
    buf.resize(static_cast<size_t>(need));
    try {
        return json::parse(buf);
    } catch (...) {
        return json::object();
    }
}

bool available(const std::string& backend)
{
    return pluginHandle(backend, /*quiet=*/true) != nullptr;
}

const std::vector<std::string>& knownBackends()
{
    static const std::vector<std::string> backends = { "cuda", "rocm", "vulkan" };
    return backends;
}

} // namespace gpu_plugin
