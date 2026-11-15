/*
 * < GPU backend plugin ABI version >
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
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

// Claude Generated (Nov 2026, S0): one version number shared by the plugin loader (gpu_plugin.cpp) and every GPU
// plugin (libcurcuma_cuda.so, libcurcuma_rocm.so, libcurcuma_vulkan.so).
//
// The plugins return C++ objects (ComputationalMethod subclasses) to the executable and call back into it; both
// sides must agree on the vtable of ComputationalMethod and on the layout of every class or struct that crosses
// the boundary (GFNFF, Mol, json, Eigen types, ...). The compiler cannot check that across a dlopen. A plugin built
// from another source state can call the wrong virtual function without any message (example: A6 inserted
// ComputationalMethod::hasPeriodicEnergy() in the middle of the vtable).
//
// Each plugin exports `extern "C" int curcuma_<backend>_abi_version()` (CURCUMA_GPU_PLUGIN_DEFINE_ABI_VERSION below).
// The loader calls it right after dlopen, before any plugin object exists, and refuses the plugin (warning, counted
// GPU fallback, CPU; -gpu_strict true stops with exit code 3) when the symbol is missing or the value differs.
//
// RAISE CURCUMA_GPU_PLUGIN_ABI_VERSION BY ONE whenever a change touches what a plugin and the executable share:
//  - ComputationalMethod (any virtual function added, removed, reordered or with a changed signature; data members),
//  - the C entry points curcuma_<backend>_* (names, signatures, ownership rules),
//  - the layout of shared classes and structs used on both sides (GFNFF and its members, FFWorkspace, Mol,
//    EnergyCalculator-facing types, parameter structs passed by reference),
//  - compile definitions that change such layouts (both sides already share EIGEN_MAX_ALIGN_BYTES, checked separately).
// When in doubt, raise it: the cost is one rebuild of the plugin, the cost of not raising it is undefined behaviour.
//
// History: 0 = plugins built before this check (no symbol); 1 = Nov 2026, after A6 (hasPeriodicEnergy in the vtable);
// 2 = Nov 2026, S1 (GFNFF gained the data member m_gpu_periodic; GFNFF is shared across the boundary);
// 3 = Nov 2026, N1 (GFNFF gained m_want_virial; FFWorkspace and FFAccumulator gained the virial members);
// 4 = Nov 2026, N2 (GFNFF gained m_cutoffs_frozen, m_frozen_d4, m_lists_rebuild_pending, m_cell_update_count);
// 5 = Nov 2026, S-HB (FFWorkspace gained m_hb_scale, the diagnostic hb_scale factor);
// 6 = Nov 2026, V-LRC (FFVirialComponents gained tail_lrc; FFWorkspace gained m_disp_tail_e6/e8/lrc; GFNFF m_pressure_lrc);
// 7 = Nov 2026, T-SCALE (FFWorkspace gained m_repulsion_scale, m_dispersion_scale, m_coulomb_scale, the term-scale factors);
// 8 = Nov 2026, WP-PERIODIC-XTB P4/G1 (ComputationalMethod gained the virtual xtbInstance(); XTB gained m_pbc_device_scf,
//     m_want_virial, m_virial_valid, m_virial, m_disp_virial, m_pbc_mp_tblite_factors; XTB is shared with the xTB GPU
//     wrappers compiled into the plugins).

#pragma once

#define CURCUMA_GPU_PLUGIN_ABI_VERSION 8

/// Defines `extern "C" int curcuma_<backend>_abi_version()` in a plugin entry file (use at file scope).
#define CURCUMA_GPU_PLUGIN_DEFINE_ABI_VERSION(backend)                                                                \
    extern "C" int curcuma_##backend##_abi_version() { return CURCUMA_GPU_PLUGIN_ABI_VERSION; }
