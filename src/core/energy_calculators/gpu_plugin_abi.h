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
// 9 = Nov 2026, WP-PERIODIC-XTB G4a (GpuScfBackend gained the virtual dispersionATMPeriodic; the backend vtable is
//     called by XTB in the executable and implemented in the plugin).
// 10 = Nov 2026, WP-PERIODIC-XTB G3 (GpuScfBackend gained the virtual multipoleMatricesPeriodic).
// 11 = Nov 2026, WP-PERIODIC-XTB G5a (GpuScfBackend gained the virtual multipoleGradientPeriodic).
// 12 = Nov 2026, WP-PERIODIC-XTB G4b (GpuScfBackend gained beginDispersionImages/dispersionDedqPeriodic; XTB gained
//      m_pbc_d4_dev_ready).
// 13 = Nov 2026, WP-PERIODIC-XTB G2/G5b/G5c/G5d (GpuScfBackend gained the virtuals ewaldGammaPeriodic,
//      energyWeightedDensity, cnChainGradientPeriodic and repulsionPeriodic; XTB gained m_pbc_dev_wfn).
// 14 = Nov 2026, WP-PERIODIC-XTB G1b (GpuScfBackend::beginComputedPeriodic; XTB gained m_pbc_sx_deferred).
// 15 = Nov 2026, WP-PERIODIC-XTB G1c (GpuScfBackend::setHostLoopPseudo).
// 16 = Nov 2026, WP-PERIODIC-XTB E1 (GpuScfBackend::setHostLoopDensityComplement).
// 17 = Nov 2026, WP-PERIODIC-XTB C1 (GpuScfBackend::d4CnFoldPeriodic).
// 18 = Nov 2026, WP-PERIODIC-XTB H1 (GpuScfBackend::h0PulayGradientPeriodic).
// 19 = Nov 2026, WP-PERIODIC-XTB H2/H3 (GpuScfBackend::ewaldGradientPeriodic, koGradientPeriodic, cnPeriodic).
// 20 = Nov 2026, WP-PERIODIC-XTB H4 (GpuScfBackend::d3DispersionPeriodic).
// 21 = Nov 2026, WP-PERIODIC-XTB H5 (GpuScfBackend::d4PairsPeriodic).
// 22 = Nov 2026, WP-PERIODIC-XTB R2-R4 (GpuScfBackend::beginPotentialPeriodic, downloadConvergedPotentials).
// 23 = Nov 2026, WP-PERIODIC-XTB R5 (GpuScfBackend::periodicGammaResident).
// 24 = Nov 2026, WP-PERIODIC-XTB A1r/B1 (GpuScfBackend::setEigenRefinement, d4CovalentCnPeriodic).
// 25 = Nov 2026, merge of spin polarisation on CUDA (GpuScfBackend spin loop, finalizeSpin) and the EEQ sytrf solver
//      (Core symbol EEQSolver::hostSytrsUpper, EEQSolveMethod::Sytrf before Auto).

#pragma once

#define CURCUMA_GPU_PLUGIN_ABI_VERSION 25

/// Defines `extern "C" int curcuma_<backend>_abi_version()` in a plugin entry file (use at file scope).
#define CURCUMA_GPU_PLUGIN_DEFINE_ABI_VERSION(backend)                                                                \
    extern "C" int curcuma_##backend##_abi_version() { return CURCUMA_GPU_PLUGIN_ABI_VERSION; }
