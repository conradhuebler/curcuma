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
 *   ./curcuma -sp mol.xyz -method gfnff -gpu rocm   # GPU
 *   ./curcuma -sp mol.xyz -method gfnff             # CPU (default)
 */

#pragma once

#ifdef USE_ROCM

#include "../ff_methods/rocm/ff_workspace_hip.h"
#include "../ff_methods/rocm/eeq_solver_hip.h"
#include "gfnff_gpu_method_impl.h"

/// Claude Generated (Oct 2026): defined in ff_methods/rocm/eeq_solver_hip.hiph (same plugin): copies the
/// charges of the last small-system device EEQ solve of this host thread when device and n match.
bool curcumaRocmEEQHostCharges(double* host, const double* device, int n);

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
    /// Oct 2026 (Claude Generated): the WP7-E projected PCG (solveWithDeviceRHSAndGPUProjectedPCG)
    /// is ported (diagonal Jacobi, or the local inverse below).
    static constexpr bool has_device_ppcg = true;
    /// Claude Generated (Oct 2026): the local-inverse preconditioner of the projected PCG, built on the
    /// host (EEQSolver::buildLocalInverse) and uploaded (setLocalInversePreconditioner).
    static constexpr bool has_localinv_upload = true;
    /// Claude Generated (Oct 2026): device build of the local inverse (requestLocalInverseDeviceBuild,
    /// localInverseDeviceStats); the default, CURCUMA_LOCALINV_HOST_BUILD=1 selects the host build.
    static constexpr bool has_localinv_device_build = true;
    /// No multi-GPU split on ROCm (CUDA only, Sep 2026).
    static constexpr bool has_multi_gpu_split = false;
    /// K1 (Oct 2026, Claude Generated): CN pair list as a Verlet list with skin
    /// (gpu_cn_pair_skin_bohr), rebuilt on its own displacement test (cnPairListStale()).
    static constexpr bool has_cn_pair_skin = true;
    /// K2/K3 (Oct 2026): dispersion and the device EEQ solve overlap dc6dcn on own streams
    /// (recordEEQInputsReady(), launchDispersionGradientEarly(), EEQ setInputReadyEvent()).
    /// Claude Generated (Oct 2026): on ROCm only K2, and only in the single-stream mode (below 500 atoms):
    /// the dispersion gather is enqueued right after the Gaussian weights, so the GPU works on it while the
    /// host solves the EEQ. K3 is a no-op here (no input event, setStreamPriorityHigh() returns false).
    static constexpr bool has_stream_overlap = true;
    /// Claude Generated (Oct 2026): -gfnff.gpu_staged_list_upload reaches the ROCm workspace; there it selects the
    /// pinned staged upload of the repulsion lists only (FFWorkspaceHip::updateRepulsion). The D4 list is built on
    /// the device, the HB list is packed by -threads host threads regardless of the switch.
    static constexpr bool has_staged_list_upload = true;
    /// Claude Generated (Oct 2026): H-bond triples enumerated on the device (FFWorkspaceHip::detectHBondsOnDevice(),
    /// rows by (case, A) instead of a key sort); -gfnff.gpu_hb_detect_on_device false restores the host detection.
    static constexpr bool has_device_hb_detect = true;
    /// Claude Generated (Oct 2026): -gfnff.gpu_disp_pairs_on_device defaults to true on ROCm as well (ordered
    /// device build with the gather CSR on the device, FFWorkspaceHip::generateDispersionPairListOnGPU);
    /// -gfnff.gpu_disp_pairs_on_device false restores the host build.
    static constexpr bool default_disp_pairs_on_device = true;
    /// WP9-A (Claude Generated, Oct 2026): -gfnff.gpu_rep_pairs_on_device reaches the ROCm workspace
    /// (FFWorkspaceHip::generateNonbondedRepulsionOnGPU, ordered row build with the gather CSR on the device);
    /// false keeps the staged host build (FFWorkspaceHip::updateRepulsion).
    static constexpr bool has_device_repulsion_build = true;

    /// Deliverable 3 (Jun 2026): from 16 fragments on, the dense device Cholesky (O(N^3) plus a
    /// host Schur step over nfrag+1 columns) is not used. Oct 2026: there the device projected PCG
    /// runs when the solver resolves to ppcg; otherwise, or when it fails, the exact CPU EEQ.
    static constexpr int default_eeq_cpu_fragment_threshold = 16;

    /// Blocking device->host copy of n doubles (GPU-Schur charge download).
    /// Defined in gfnff_hip_method.cpp so <hip/hip_runtime.h> stays out of this header.
    static void downloadDoubles(double* host, const double* device, int n);

    /// Claude Generated (Oct 2026): charge download after a device EEQ solve. The small-system solve
    /// (rocm/eeq_solver_hip.hiph) already read its charges back with its status; they are copied from
    /// that pinned host copy, anything else is downloaded. Detected by gfnffBackendEEQChargeDownload in
    /// gfnff_gpu_method_impl.h; a backend without the member (CUDA) calls downloadDoubles as before.
    static void downloadEEQCharges(double* host, const double* device, int n)
    {
        if (!curcumaRocmEEQHostCharges(host, device, n))
            downloadDoubles(host, device, n);
    }

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
    using Base = GFNFFGpuMethodImpl<GFNFFRocmBackend>;

public:
    /// Claude Generated (Oct 2026, WP9): reads the ROCm-only switches (RocmGfnffOptions) from the config.
    GFNFFHipComputationalMethod(const std::string& method_name, const json& config)
        : Base(method_name, config)
        , m_rocm_options(readRocmOptions(config))
    {
    }

    /// The device classes are created inside these two calls; the options are published around them.
    bool setMolecule(const Mol& mol) override
    {
        RocmOptionsScope scope(m_rocm_options);
        return Base::setMolecule(mol);
    }
    double calculateEnergy(bool gradient = false) override
    {
        RocmOptionsScope scope(m_rocm_options);
        return Base::calculateEnergy(gradient);
    }
    void setParameters(const json& params) override
    {
        Base::setParameters(params);
        m_rocm_options = readRocmOptions(params);
    }

    /// Claude Generated (Oct 2026, WP9): gfnff scope first, then the top level (as the wrapper's cfg_get_bool).
    static RocmGfnffOptions readRocmOptions(const json& config)
    {
        RocmGfnffOptions o;
        const json* g = (config.contains("gfnff") && config["gfnff"].is_object()) ? &config["gfnff"] : nullptr;
        auto get = [&](const char* k, bool d) {
            if (g && g->contains(k)) return g->value(k, d);
            return config.value(k, d);
        };
        o.dc6dcn_tables = get("rocm_dc6dcn_tables", o.dc6dcn_tables);
        o.device_schur  = get("rocm_device_schur", o.device_schur);
        o.fixed_point   = get("rocm_fixed_point", o.fixed_point);
        return o;
    }

private:
    /// Sets rocmGfnffOptionsSlot() for the lifetime of the scope and restores the previous value.
    struct RocmOptionsScope {
        const RocmGfnffOptions* prev;
        explicit RocmOptionsScope(const RocmGfnffOptions& o) : prev(rocmGfnffOptionsSlot()) { rocmGfnffOptionsSlot() = &o; }
        ~RocmOptionsScope() { rocmGfnffOptionsSlot() = prev; }
    };
    RocmGfnffOptions m_rocm_options;
};

#endif // USE_ROCM
