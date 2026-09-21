/*
 * <GFNFFHipComputationalMethod — ROCm backend of the shared GFN-FF GPU wrapper>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (March 2026): ComputationalMethod adapter for gfnff GPU path.
 * Claude Generated (Sep 2026): the ~1400 lines of wrapper logic moved to the shared
 * class template GFNFFGpuMethodImpl<Backend> (gfnff_gpu_method_impl.h), which the CUDA
 * wrapper uses as well.  This file is now only the ROCm backend traits plus the
 * concrete class that keeps the established public name.
 * Available only when compiled with USE_ROCM=ON.
 *
 * Usage:
 *   ./curcuma -sp mol.xyz -method gfnff -gpu rocm    # Explicit GPU
 *   ./curcuma -sp mol.xyz -method gfnff -gpu auto   # GPU if available
 *   ./curcuma -sp mol.xyz -method gfnff             # CPU (default)
 */

#pragma once

#ifdef USE_ROCM

#include "../ff_methods/rocm/ff_workspace_hip.h"
#include "../ff_methods/rocm/eeq_solver_hip.h"
#include "gfnff_gpu_method_impl.h"

/**
 * @brief ROCm backend traits for GFNFFGpuMethodImpl.
 *
 * Claude Generated (Sep 2026).  Everything the shared wrapper needs to know about the
 * backend it drives: the two device classes, the strings it puts in log lines, and the
 * one capability that actually differs between CUDA and ROCm (the device Schur EEQ).
 */
struct GFNFFRocmBackend {
    using Workspace = FFWorkspaceHip;
    using EEQSolver = EEQSolverHip;

    /// Backend name used in log prose.
    static constexpr const char* name = "ROCm";
    /// Device workspace class name, used in the init-failure message.
    static constexpr const char* workspace_name = "FFWorkspaceHip";

    /// The device Schur EEQ solves (WP5-A / WP7-A / WP7-B / WP7-C) are NOT ported to
    /// HIP — they return false and the wrapper falls back to the WP2 GPU solve + CPU
    /// Schur complement.  That is the expected path here, not a numerical failure, so
    /// the wrapper reports it at verbosity 3 instead of warning.
    static constexpr bool has_device_schur = false;

    /// Deliverable 3 (Jun 2026): from 16 fragments on, the EEQ goes to the exact CPU
    /// PCG — the device path would do a dense N x N Cholesky (O(N^3)) for nfrag>1.
    static constexpr int default_eeq_cpu_fragment_threshold = 16;

    /// Blocking device->host copy of n doubles (GPU-Schur charge download).
    /// Defined in gfnff_hip_method.cpp so <hip/hip_runtime.h> stays out of this header.
    static void downloadDoubles(double* host, const double* device, int n);

    /// Multi-GPU (Claude Generated, Sep 2026): number of visible devices, and make `device`
    /// current on the calling thread (the runtime's current device is per host thread).
    static int  deviceCount();
    static bool setDevice(int device);
};

// The wrapper is instantiated in exactly one translation unit (gfnff_hip_method.cpp),
// same as before the template refactor.
extern template class GFNFFGpuMethodImpl<GFNFFRocmBackend>;

/**
 * @brief GPU-accelerated GFN-FF via HIP/ROCm (method name: "gfnff" with -gpu rocm)
 *
 * Kept as a real class (not an alias): the plugin entry point creates it by name.
 */
// NOTE: deliberately NOT 'final' — marking it final lets the compiler
// devirtualise calls through a GFNFFHipComputationalMethod* into direct references to the
// (extern-template) base members, which callers outside this plugin cannot resolve.
class GFNFFHipComputationalMethod : public GFNFFGpuMethodImpl<GFNFFRocmBackend> {
public:
    using GFNFFGpuMethodImpl<GFNFFRocmBackend>::GFNFFGpuMethodImpl;
};

#endif // USE_ROCM
