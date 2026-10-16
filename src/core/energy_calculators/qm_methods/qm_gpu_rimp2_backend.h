/*
 * < Device seam of the native QM engine: RI-MP2 pair energies >
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0
 */

#pragma once

// Claude Generated (Oct 2026, RI-MP2). Separate small header like qm_gpu_mp2_backend.h: only
// cuda/qm_gpu_rimp2.cu, gpu_plugin.{h,cpp} and qm_rimp2.cpp see it.

namespace qmint {

/// Input of one device RI-MP2 energy evaluation. Host arrays; the plugin copies them.
struct QMGpuRIMP2Input {
    int o = 0;                    ///< active occupied orbitals
    int v = 0;                    ///< virtual orbitals
    int naux = 0;                 ///< fitted auxiliary functions
    const double* B = nullptr;    ///< fitted factors [i][a][P] (o*v rows of naux), host
    const double* eo = nullptr;   ///< active occupied orbital energies (o), host
    const double* ev = nullptr;   ///< virtual orbital energies (v), host
    int device = 0;
    double max_device_mb = 0.0;   ///< > 0: pretend the device has at most this much free memory
    double batch_mb = 0.0;        ///< > 0: upper bound for one (ia|jb) batch, ib * v * o * v * 8 bytes
};

/// Output (Hartree, milliseconds).
struct QMGpuRIMP2Output {
    double e_os = 0.0, e_ss = 0.0;
    double t_upload_ms = 0.0;     ///< B, energies
    double t_gemm_ms = 0.0;       ///< (ia|jb) = B_i B^T GEMMs
    double t_energy_ms = 0.0;     ///< pair-energy reduction
    int batch = 0, nbatches = 0;  ///< occupied orbitals per batch
};

}  // namespace qmint
