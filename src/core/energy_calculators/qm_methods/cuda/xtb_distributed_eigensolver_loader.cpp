/*
 * < Loader for the optional multi-GPU eigensolver library libcurcuma_cuda_mgpu.so >
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

// Claude Generated (Sep 2026, multi-GPU step 3).
//
// Compiled into libcurcuma_cuda.so. The cuSOLVERMp/NCCL/cusolverMg code lives in
// libcurcuma_cuda_mgpu.so next to it and is dlopen'ed on the first request, so a host without
// those libraries still runs every single-GPU calculation; only -gpu_eigensolver_devices then
// reports why it is unavailable. The library is never closed (the solver objects' code lives there).

#ifdef USE_CUDA

#include "xtb_distributed_eigensolver.h"

#include <dlfcn.h>

#include <mutex>

namespace curcuma {
namespace xtb {
namespace gpu {

namespace {

using CreateFn = DistributedEigensolver* (*)(const std::string*, const std::vector<int>*, int, std::string*);
using BackendsFn = const char* (*)();

struct MgpuLibrary {
    void* handle = nullptr;
    CreateFn create = nullptr;
    BackendsFn backends = nullptr;
    std::string error;
};

const MgpuLibrary& library()
{
    static MgpuLibrary lib;
    static std::once_flag once;
    std::call_once(once, []() {
        // Look next to this plugin first (build tree / install dir), then the loader path.
        std::string path = "libcurcuma_cuda_mgpu.so";
        Dl_info info{};
        if (::dladdr(reinterpret_cast<void*>(&library), &info) && info.dli_fname) {
            const std::string self = info.dli_fname;
            const auto slash = self.rfind('/');
            if (slash != std::string::npos) path = self.substr(0, slash + 1) + path;
        }
        lib.handle = ::dlopen(path.c_str(), RTLD_NOW | RTLD_LOCAL);
        if (!lib.handle) lib.handle = ::dlopen("libcurcuma_cuda_mgpu.so", RTLD_NOW | RTLD_LOCAL);
        if (!lib.handle) {
            const char* e = ::dlerror();
            lib.error = std::string("libcurcuma_cuda_mgpu.so not loadable") + (e ? std::string(": ") + e : "");
            return;
        }
        lib.create = reinterpret_cast<CreateFn>(::dlsym(lib.handle, "curcuma_cuda_mgpu_create"));
        lib.backends = reinterpret_cast<BackendsFn>(::dlsym(lib.handle, "curcuma_cuda_mgpu_backends"));
        if (!lib.create || !lib.backends)
            lib.error = "libcurcuma_cuda_mgpu.so lacks its entry points";
    });
    return lib;
}

} // namespace

std::string DistributedEigensolver::availableBackends()
{
    const auto& lib = library();
    return lib.backends ? std::string(lib.backends()) : std::string("none (") + lib.error + ")";
}

std::unique_ptr<DistributedEigensolver> DistributedEigensolver::create(const std::string& backend,
                                                                       const std::vector<int>& devices,
                                                                       int block, std::string& why)
{
    if (devices.size() < 2) {
        why = "needs at least two devices";
        return nullptr;
    }
    const auto& lib = library();
    if (!lib.create) {
        why = lib.error;
        return nullptr;
    }
    return std::unique_ptr<DistributedEigensolver>(lib.create(&backend, &devices, block, &why));
}

} // namespace gpu
} // namespace xtb
} // namespace curcuma

#endif // USE_CUDA
