/*
 * <Native xTB ROCm/HIP Method Wrapper — GFN1 / GFN2 on ROCm>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0.
 *
 * Claude Generated (2026-06): ComputationalMethod adapter for the ROCm path of
 * native GFN1/GFN2 — the AMD/HIP sibling of XtbGpuComputationalMethod (CUDA).
 * Selected by the factory when `-gpu rocm` is
 * given and the build has USE_ROCM.
 *
 * Design: OWNS a NativeXtbMethod (the validated CPU pipeline: config, large-system
 * modes, error handling, property accessors) + an XtbHipContext (hipBLAS/rocSOLVER
 * handles + stream). Every ComputationalMethod call forwards to the CPU pipeline; the
 * GPU is reached through the hooks the constructor installs on the owned XTB. That
 * shared wrapper logic lives ONCE in XtbGpuAdapter (xtb_gpu_adapter.h, Claude Generated
 * Sep 2026); only the ROCm-specific eigensolver hook, the resident backend and the
 * mixed-precision policy are in the .cpp.
 *
 * Usage:
 *   ./curcuma -sp mol.xyz -method gfn2 -gpu rocm   # explicit ROCm
 *   ./curcuma -sp mol.xyz -method gfn1 -gpu rocm   # GPU
 *   ./curcuma -sp mol.xyz -method gfn2             # CPU (default)
 */

#pragma once

#ifdef USE_ROCM

#include "xtb_gpu_adapter.h"        // XtbGpuAdapter (shared wrapper logic)
#include "rocm/xtb_hip_context.h"   // XtbHipContext (pimpl — no HIP headers)

/**
 * @brief GPU-accelerated native GFN1/GFN2 (method "gfn1"/"gfn2" with -gpu rocm).
 *
 * A thin instantiation: the whole ComputationalMethod surface (forwarding, device
 * handshake + logging, CPU fallback, gpuActive()) comes from the shared adapter; the
 * constructor adds the rocSOLVER eigensolver hook and the device-resident SCF backend.
 * Claude Generated.
 */
class XtbHipComputationalMethod
    : public curcuma::xtb::gpu::XtbGpuAdapter<curcuma::xtb::gpu::XtbHipContext> {
public:
    explicit XtbHipComputationalMethod(curcuma::xtb::MethodType method,
                                       const json& config = json{});
    ~XtbHipComputationalMethod() override;

    /**
     * @brief Energy (and gradient) as in XtbGpuAdapter, then one warning per new solver
     *        fallback of the resident ROCm SCF (XtbHipContext::solverFallbackEvents()).
     *        Each kind of fallback is reported once per method object, not per geometry; while
     *        the logger is at verbosity 0 (optimiser/MD steps) it stays pending and is printed
     *        by a later call or, at level 1, by the destructor.
     *        Claude Generated (Oct 2026).
     */
    double calculateEnergy(bool gradient = false) override;

private:
    /// Prints a warning for every bit of the context's fallback mask not reported yet;
    /// below verbosity 1 only when final_report is true (then at level 1).
    void reportSolverFallbacks(bool final_report = false);
    unsigned m_reported_fallbacks = 0;   ///< SolverFallbackEvent bits already printed
    bool m_reported_mp_fallback = false; ///< -xtb.rocm_device_multipole fallback printed (Claude Generated, Oct 2026)
    bool m_reported_pot_fallback = false; ///< -xtb.rocm_device_potential fallback printed (Claude Generated, Oct 2026)
    bool m_reported_res_fallback = false; ///< -xtb.rocm_resident_loop fallback printed (Claude Generated, Oct 2026)
    bool m_reported_solv_fallback = false; ///< -xtb.rocm_device_solvation fallback printed (WP6, Claude Generated, Oct 2026)
    bool m_reported_pseudo_fallback = false; ///< -scf_pseudo_diag path turned off, printed (WP7, Claude Generated, Oct 2026)
    bool m_pseudo_diag = false;              ///< -scf_pseudo_diag active on ROCm (WP7): counts printed at verbosity 2
    bool m_reported_ps_fallback[3] = {false, false, false}; ///< WP8a post-SCF D4 / D3 / CN fallback printed (Claude Generated, Oct 2026)
    bool m_reported_gr_fallback[3] = {false, false, false}; ///< WP8b gradient potentials / tail / lazy copies fallback printed (Claude Generated, Oct 2026)
};

#endif // USE_ROCM
