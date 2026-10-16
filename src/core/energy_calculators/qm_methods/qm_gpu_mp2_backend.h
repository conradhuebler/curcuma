/*
 * < Device seam of the native QM engine: canonical MP2 (WP11) >
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0
 */

#pragma once

// Claude Generated (Oct 2026, WP11 step 1 GPU variant). Separate small header, like
// qm_gpu_scf_backend.h: only cuda/qm_gpu_mp2.cu, gpu_plugin.cpp and qm_mp2.cpp see it, so a
// change here does not recompile the heavy ERI kernel TUs.

namespace qmint {

/// Input of one device MP2 evaluation. Host arrays are plain host memory; the plugin copies them.
struct QMGpuMP2Input {
    int n = 0;                         ///< active basis dimension
    int n_frozen = 0;                  ///< lowest MOs left out of the correlation treatment
    int n_occ = 0;                     ///< doubly occupied MOs (including the frozen ones)
    const double* C = nullptr;         ///< MO coefficients, n x n row-major (columns = MOs), host
    const double* eps = nullptr;       ///< orbital energies (n), host
    /// AO tensor (n^4, ERITensor layout): a DEVICE pointer on `device` (stored tensor of the SCF,
    /// QMGpuJK::deviceERI) or, if that is nullptr, a HOST pointer that the plugin uploads.
    const double* eri_device = nullptr;
    const double* eri_host = nullptr;
    int device = 0;
    double max_device_mb = 0.0;        ///< > 0: pretend the device has at most this much free memory
    double batch_mb = 0.0;             ///< > 0: upper bound for one batch, ib (n^3 + n^2 v) * 8 bytes (-qm.mp2_batch_mb)
};

/// Output (Hartree, milliseconds).
struct QMGpuMP2Output {
    double e_os = 0.0, e_ss = 0.0;
    double t_upload_ms = 0.0;          ///< tensor (if uploaded) + C, eps
    double t_transform_ms = 0.0;       ///< the four GEMM steps
    double t_energy_ms = 0.0;          ///< pair-energy reduction
    int batch = 0, nbatches = 0;       ///< occupied orbitals per batch
};

}  // namespace qmint
