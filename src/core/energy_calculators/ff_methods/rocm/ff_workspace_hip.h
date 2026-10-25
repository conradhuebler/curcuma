/*
 * <FFWorkspaceHip - ROCm include shim for the shared FFWorkspaceGPU declaration>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (Sep 2026): cuda/ff_workspace_gpu.h is backend-neutral (Pimpl: no
 * CUDA/HIP type appears in the interface), so the ROCm build uses that one declaration
 * instead of a hipified copy that has to be kept in sync by hand.
 *
 * Why the rename instead of one shared class name: the HIP port was generated with
 * hipify, so rocm/gfnff_rocm.hip defines its members as FFWorkspaceHip::... (7k lines,
 * only compilable with a ROCm SDK).  Keeping the name also keeps the two plugins'
 * exported symbols distinct - gpu_plugin.cpp dlopen's libcurcuma_cuda.so and
 * libcurcuma_rocm.so with RTLD_GLOBAL, so identical strong symbols in both would bind
 * across plugin boundaries.  The two object-like macros below give the shared
 * declaration its ROCm names; nothing else in this build spells "FFWorkspaceGPU".
 */

#pragma once

#define FFWorkspaceGPU     FFWorkspaceHip
#define FFWorkspaceGPUImpl FFWorkspaceHipImpl

#include "../cuda/ff_workspace_gpu.h"

/**
 * @brief ROCm-only GFN-FF switches of WP9 (Parameter Registry, gfnff scope).
 *
 * Claude Generated (Oct 2026). The shared wrapper template does not know them, so the ROCm
 * method class (GFNFFHipComputationalMethod) reads them from its configuration and publishes
 * them through rocmGfnffOptionsSlot() while setMolecule() and calculateEnergy() run; the
 * FFWorkspaceHip and EEQSolverHip constructors (called from there) copy them. The member
 * defaults are the fallback defaults and must equal the PARAM defaults in gfnff.h.
 */
struct RocmGfnffOptions {
    bool dc6dcn_tables = false;   ///< -gfnff.rocm_dc6dcn_tables: tabulated dc6dcn (k_d4_partner_tables)
    bool device_schur  = true;    ///< -gfnff.rocm_device_schur: Schur step of the dense EEQ on the device
    bool fixed_point   = true;    ///< -gfnff.rocm_fixed_point: order-independent fixed-point accumulation (WP10)
};

/// Claude Generated (Oct 2026): thread-local slot of the options above (nullptr: member defaults).
/// Defined in rocm/gfnff_rocm.hip.
const RocmGfnffOptions*& rocmGfnffOptionsSlot();
