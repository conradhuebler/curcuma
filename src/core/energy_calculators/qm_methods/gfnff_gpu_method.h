/*
 * <GFNFFGPUComputationalMethod — CUDA backend of the shared GFN-FF GPU wrapper>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (March 2026): ComputationalMethod adapter for gfnff GPU path.
 * Claude Generated (Sep 2026): the ~1400 lines of wrapper logic moved to the shared
 * class template GFNFFGpuMethodImpl<Backend> (gfnff_gpu_method_impl.h), which the ROCm
 * wrapper uses as well.  This file is now only the CUDA backend traits plus the
 * concrete class that keeps the established public name.
 * Available only when compiled with USE_CUDA=ON.
 *
 * Usage:
 *   ./curcuma -sp mol.xyz -method gfnff -gpu cuda    # Explicit GPU
 *   ./curcuma -sp mol.xyz -method gfnff -gpu auto   # GPU if available
 *   ./curcuma -sp mol.xyz -method gfnff             # CPU (default)
 */

#pragma once

#ifdef USE_CUDA

#include "../ff_methods/cuda/ff_workspace_gpu.h"
#include "../ff_methods/cuda/eeq_solver_gpu.h"
#include "gfnff_gpu_method_impl.h"

/**
 * @brief CUDA backend traits for GFNFFGpuMethodImpl.
 *
 * Claude Generated (Sep 2026).  Everything the shared wrapper needs to know about the
 * backend it drives: the two device classes, the strings it puts in log lines, and the
 * one capability that actually differs between CUDA and ROCm (the device Schur EEQ).
 */
struct GFNFFCudaBackend {
    using Workspace = FFWorkspaceGPU;
    using EEQSolver = EEQSolverGPU;

    /// Backend name used in log prose.
    static constexpr const char* name = "CUDA";
    /// Device workspace class name, used in the init-failure message.
    static constexpr const char* workspace_name = "FFWorkspaceGPU";

    /// CUDA implements the device Schur EEQ solves (WP5-A / WP7-A / WP7-B / WP7-C),
    /// so a failure there is a real numerical failure and is warned about.
    static constexpr bool has_device_schur = true;

    /// 0 = never route the EEQ to the CPU PCG; the device Schur paths handle nfrag>1.
    static constexpr int default_eeq_cpu_fragment_threshold = 0;

    /// Blocking device->host copy of n doubles (GPU-Schur charge download).
    /// Defined in gfnff_gpu_method.cpp so <cuda_runtime.h> stays out of this header.
    static void downloadDoubles(double* host, const double* device, int n);

    /// Multi-GPU (Claude Generated, Sep 2026): number of visible devices, and make `device`
    /// current on the calling thread (the runtime's current device is per host thread).
    static int  deviceCount();
    static bool setDevice(int device);
};

// The wrapper is instantiated in exactly one translation unit (gfnff_gpu_method.cpp,
// compiled by nvcc) — same generated code as before the template refactor.
extern template class GFNFFGpuMethodImpl<GFNFFCudaBackend>;

/**
 * @brief GPU-accelerated GFN-FF via CUDA (method name: "gfnff" with -gpu cuda)
 *
 * Kept as a real class (not an alias): the plugin entry point creates it by name and
 * test_cases/test_gfnff_gpu.cpp dynamic_casts to it.
 */
// NOTE: deliberately NOT 'final' — marking it final lets the compiler
// devirtualise calls through a GFNFFGPUComputationalMethod* into direct references to the
// (extern-template) base members, which callers outside this plugin cannot resolve.
class GFNFFGPUComputationalMethod : public GFNFFGpuMethodImpl<GFNFFCudaBackend> {
public:
    using GFNFFGpuMethodImpl<GFNFFCudaBackend>::GFNFFGpuMethodImpl;
};

#endif // USE_CUDA
