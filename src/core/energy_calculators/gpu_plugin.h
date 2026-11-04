/*
 * < GPU backend plugin loader >
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

// Claude Generated (Jul 2026): runtime loader for the GPU backend plugins.
//
// The CUDA/ROCm/Vulkan backends are built into separate shared libraries
// (libcurcuma_cuda.so, ...) that are dlopen'd ON DEMAND — only when a GPU run is
// actually requested. This keeps the heavy GPU runtime libraries (cuBLAS/cuSOLVER,
// which run ~40 ms of DT_INIT table-building at load) OUT of the main binary's
// startup path, so CPU-only runs start in ~1 ms instead of ~50 ms.
//
// The plugin exposes C entry points (curcuma_<backend>_create_native_xtb /
// curcuma_<backend>_create_gfnff) that construct the concrete GPU ComputationalMethod
// and hand back a raw pointer the core wraps in a unique_ptr. The plugin .so is kept
// loaded for the process lifetime, so the returned object's vtable/destructor stay valid.

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "computational_method.h"   // ComputationalMethod, json (global namespace)
#include "qm_methods/qm_gpu_backend.h"  // qmint::QMGpuJK
#include "qm_methods/qm_gpu_scf_backend.h"  // qmint::QMGpuSCF (Q-G4)
#include "qm_methods/qm_gpu_xc_backend.h"  // qmint::QMGpuXC (GPU KS-DFT, Oct 2026)
#include "qm_methods/qm_gpu_grad_backend.h"  // qmint::QMGpuGrad2e (Q-G7)
#include "qm_methods/qm_gpu_mp2_grad_backend.h"  // qmint::QMGpuMP2Grad2e (MP2 gradient)
#include "qm_methods/qm_gpu_grad1e_backend.h"  // qmint::QMGpuGrad1e (Q-G8)
#include "qm_methods/qm_gpu_ecp_backend.h"     // qmecp::QMGpuECP (ECP on the device, Nov 2026)
#include "qm_methods/qm_gpu_ri_backend.h"  // qmint::QMGpuRI (WP10)
#include "qm_methods/qm_gpu_ri_lowmem_backend.h"  // qmint::QMGpuRILowMem, QMGpuRIYSolver (WP10 low memory)
#include "qm_methods/qm_gpu_mp2_backend.h"  // qmint::QMGpuMP2Input/Output (WP11)
#include "qm_methods/qm_gpu_rimp2_backend.h"  // qmint::QMGpuRIMP2Input/Output (RI-MP2)
#include "qm_methods/qm_gpu_ri_grad_backend.h"  // qmint::QMGpuRIGradInput/Output (RI gradient)

namespace gpu_plugin {

/**
 * Construct a native-xTB (GFN1/GFN2) GPU method from the given backend plugin.
 * @param backend    "cuda" | "rocm" | "vulkan"
 * @param method_type curcuma::xtb::MethodType value (GFN1=1, GFN2=2)
 * Returns nullptr if the plugin is unavailable or construction fails (caller must
 * fall back to the CPU method). Claude Generated.
 */
std::unique_ptr<ComputationalMethod> createNativeXtb(const std::string& backend,
                                                     int method_type, const json& config);

/**
 * Construct a GFN-FF GPU method from the given backend plugin.
 * Returns nullptr on failure (caller falls back to CPU). Claude Generated.
 */
std::unique_ptr<ComputationalMethod> createGfnff(const std::string& backend, const json& config);

/**
 * Native-QM J/K builder with the stored ERI tensor (n^4 doubles, row-major (mn|ls)) resident
 * on the device. nullptr when the plugin/entry point is missing or the upload fails (the
 * caller falls back to the CPU and reports it). Claude Generated (Oct 2026, Q-G1).
 */
std::unique_ptr<qmint::QMGpuJK> createQmJK(const std::string& backend, const double* eri, int n, int device);

/**
 * Native-QM J/K builder whose ERI tensor is computed ON the device from the flat tables
 * (qmint::buildEriGpuTables); nothing of size n^4 crosses the bus. nullptr with a reason in
 * `why` when the plugin/entry point is missing, the basis is unsupported (L > 2) or the
 * tensor does not fit into device memory. Claude Generated (Oct 2026, Q-G2).
 */
std::unique_ptr<qmint::QMGpuJK> createQmEriJK(const std::string& backend, const qmint::QMGpuEriTables& tables,
                                              int device, std::string& why);

/**
 * Native-QM integral-direct J/K builder (Q-G3): no stored tensor, every screened shell quartet
 * is recomputed per build on the device and contracted with P at once; O(n^2) device memory.
 * `tables` from qmint::buildEriGpuTables(..., list_quartets = false). nullptr with a reason in
 * `why` when the plugin/entry point is missing, the basis is unsupported or an allocation
 * fails. Claude Generated (Oct 2026).
 */
std::unique_ptr<qmint::QMGpuJK> createQmDirectJK(const std::string& backend, const qmint::QMGpuEriTables& tables,
                                                 int device, std::string& why);

/**
 * Native-QM device-resident SCF state (Q-G4) for n basis functions: device matrices plus
 * cuBLAS/cuSOLVER handles; the SCF driver stays in QMEngine::runSCF. nullptr with a reason in
 * `why` when the plugin/entry point is missing or an allocation fails. Claude Generated (Oct 2026).
 */
std::unique_ptr<qmint::QMGpuSCF> createQmScf(const std::string& backend, int n, int device, std::string& why);

/**
 * Native-QM Kohn-Sham exchange-correlation on the device (GPU KS-DFT, Oct 2026): grid resident, AO values streamed
 * per block of points; V_xc/E_xc and the XC gradient incl. grid response (qm_gpu_xc_backend.h). nullptr with a reason
 * in `why` when the plugin/entry point (curcuma_<backend>_create_qm_xc_fxc) is missing, the basis/functional is not
 * supported (shells above f, more than 32 primitives) or the device memory does not hold one block. Claude Generated.
 */
std::unique_ptr<qmint::QMGpuXC> createQmXC(const std::string& backend, const qmint::QMGpuXCSetup& setup, int device,
                                           std::string& why);

/**
 * Native-QM two-electron gradient builder (Q-G7): integral-direct derivative integrals contracted
 * with the density on the device. `tables`/`extra` from qmint::buildEriGradGpuTables (cartesian,
 * raised Hermite tables). nullptr with a reason in `why` when the plugin/entry point is missing,
 * the basis is unsupported or an allocation fails. Claude Generated (Oct 2026).
 */
std::unique_ptr<qmint::QMGpuGrad2e> createQmGrad2e(const std::string& backend, const qmint::QMGpuEriTables& tables,
                                                   const qmint::QMGpuGradExtra& extra, int device, std::string& why);

/**
 * MP2-gradient two-electron derivative contraction on the device (RHF + separable + non-separable 2PDM,
 * cartesian). Same tables as createQmGrad2e. nullptr with a reason in `why` when the plugin/entry point is
 * missing, the basis is unsupported or the n^4 2PDM does not fit. Claude Generated (Oct 2026).
 */
std::unique_ptr<qmint::QMGpuMP2Grad2e> createQmMp2Grad2e(const std::string& backend, const qmint::QMGpuEriTables& tables,
                                                         const qmint::QMGpuGradExtra& extra, int device, std::string& why);

/**
 * Native-QM one-electron gradient builder (Q-G8): S, T, V derivative integrals contracted with the
 * density and the energy-weighted density on the device. `tables` are the Q-G7 gradient tables
 * (qmint::buildEriGradGpuTables), `extra` adds shell centres, nuclei and charges. nullptr with a
 * reason in `why` when the plugin/entry point is missing, the basis is unsupported (L > 2) or an
 * allocation fails. Claude Generated (Oct 2026).
 */
std::unique_ptr<qmint::QMGpuGrad1e> createQmGrad1e(const std::string& backend, const qmint::QMGpuEriTables& tables,
                                                   const qmint::QMGpuGrad1eExtra& extra, int device, std::string& why);

/**
 * @brief ECP matrix / gradient device object (entry curcuma_<backend>_create_qm_ecp_v1) for the functions of
 * `setup` (qmecp::deviceFunctions); nullptr with `why` when the plugin, the entry or the layout check fails.
 * Claude Generated (Nov 2026).
 */
std::unique_ptr<qmecp::QMGpuECP> createQmECP(const std::string& backend, const qmecp::QMGpuEcpSetup& setup, int device,
                                             std::string& why);

/**
 * Native-QM density-fitting contractions (WP10): uploads the fitted 3-index tensor T (nfit x n(n+1)/2,
 * packed lower triangle per row, row-major, host pointer) once; J/K builds then run as cuBLAS GEMMs
 * (qm_gpu_ri_backend.h). Entry point curcuma_<backend>_create_qm_ri_packed (renamed with the packed
 * layout, Oct 2026, so an older plugin is refused instead of misreading the tensor). nullptr with a
 * reason in `why` when the plugin/entry point is missing or the tensor does not fit.
 * Claude Generated (Oct 2026).
 */
std::unique_ptr<qmint::QMGpuRI> createQmRI(const std::string& backend, const double* T, int nfit, int n, int device,
                                           std::string& why);

/**
 * RI-JK metric solve on the device (WP10 scaling): T := L^-1 T in place (T: N x ncol row-major host
 * array, L: N x N lower triangular, column-major), cublasDtrsm on streamed column chunks. False with a
 * reason; `touched` then says whether T was partly overwritten. Claude Generated (Oct 2026).
 */
bool qmRiLowerSolve(const std::string& backend, double* T, int N, long ncol, const double* L, int device, int host_threads,
                    bool& touched, std::string& why);

/**
 * Resource-saving RI contractions (WP10, Oct 2026): the fitted tensor T (nfit x n(n+1)/2 packed rows, host array of
 * double, or float when fp32) resident on the device or streamed from host memory per build. stream: 0 = auto
 * (resident when it fits, else streamed), 1 = always streamed, 2 = resident only (nullptr when it does not fit).
 * max_device_mb > 0 caps the device memory considered free (testing); chunk_mb > 0 sets the streamed chunk size.
 * Entry point curcuma_<backend>_create_qm_ri_lowmem. Claude Generated (Oct 2026).
 */
std::unique_ptr<qmint::QMGpuRILowMem> createQmRILowMem(const std::string& backend, const void* T, bool fp32, int nfit, int n,
                                                       int device, int stream, double max_device_mb, double chunk_mb,
                                                       int host_threads, std::string& why);

/**
 * Batched RI-JK exchange on the device (WP10, Oct 2026): keeps L (naux x naux, column-major lower Cholesky factor)
 * on the device, then per call X = L^-1 Y and K += sign X X^T on streamed column chunks of Y (qm_gpu_ri_lowmem_backend.h).
 * Entry point curcuma_<backend>_create_qm_ri_ysolver. Claude Generated (Oct 2026).
 */
std::unique_ptr<qmint::QMGpuRIYSolver> createQmRIYSolver(const std::string& backend, const double* L, int naux, int n, int device,
                                                         double max_device_mb, int host_threads, std::string& why);

/**
 * Native-QM canonical MP2 on the device (WP11): four-index transformation with DGEMM and the
 * pair-energy sum, FP64. False with a reason in `why` when the plugin/entry point is missing, the
 * device memory is too small or a device call fails. Claude Generated (Oct 2026).
 */
bool qmMp2Energy(const std::string& backend, const qmint::QMGpuMP2Input& in, qmint::QMGpuMP2Output& out,
                 std::string& why);

/**
 * Native-QM RI-MP2 pair energies on the device: (ia|jb) = sum_P B_ia^P B_jb^P as cuBLAS DGEMM per
 * occupied batch plus the pair-energy sum, FP64. False with a reason in `why` when the plugin/entry
 * point is missing, the device memory is too small or a device call fails. Claude Generated (Oct 2026).
 */
bool qmRiMp2Energy(const std::string& backend, const qmint::QMGpuRIMP2Input& in, qmint::QMGpuRIMP2Output& out,
                   std::string& why);

/**
 * Native-QM RI gradient, 3-index part (derivative integrals (mn|P)^x contracted with the host's
 * cartesian weight blocks), added into grad (natoms x 3, Eh/Bohr). False with a reason in `why` when
 * the plugin/entry point is missing, a limit is exceeded or a device call fails. Claude Generated (Oct 2026).
 */
bool qmRiGrad(const std::string& backend, const qmint::QMGpuRIGradInput& in, double* grad, qmint::QMGpuRIGradOutput& out,
              int device, std::string& why);

/**
 * Resolve a raw entry point from a backend plugin, or nullptr if the plugin or the symbol
 * is unavailable. Lets consumers other than the method factory (the delta-ML linear
 * algebra, for one) reuse the same dlopen cache instead of opening the library again.
 * Claude Generated (July 2026).
 */
void* resolveEntryPoint(const std::string& backend, const std::string& symbol);
/// True if the plugin library libcurcuma_<backend>.so can be loaded (cached; a failed probe
/// is not retried and, unlike the create* calls, is silent). Used for
/// the "was this backend built?" checks — the core no longer needs USE_CUDA/USE_ROCM/
/// USE_VULKAN for dispatch (Claude Generated, Sep 2026).
bool available(const std::string& backend);

/// The backend names the loader knows, in priority order.
const std::vector<std::string>& knownBackends();

/**
 * @brief Number of devices the backend's runtime can see (after CUDA_VISIBLE_DEVICES etc.).
 * Returns 0 when the plugin is not loadable. A plugin built before the multi-GPU ABI
 * (no `curcuma_<backend>_device_count` symbol) reports 1. Claude Generated (Sep 2026).
 */
int deviceCount(const std::string& backend);

/**
 * @brief Static properties of one device: {"index","name","memory_total_bytes",
 * "compute_capability", ...}. Empty object when unavailable. Creates no device context.
 * Claude Generated (Sep 2026).
 */
json deviceInfo(const std::string& backend, int index);

/**
 * @brief Distributed-eigensolver backends this plugin was built with ("mp", "mg", "mp,mg").
 * Empty when the plugin does not export the symbol. Claude Generated (Sep 2026).
 */
std::string mgpuBackends(const std::string& backend);

} // namespace gpu_plugin
