/*
 * < Device seam of the native QM engine: density-fitting contractions (WP10) >
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0
 */

#pragma once

// Claude Generated (Oct 2026, WP10). Own header (pattern of qm_gpu_scf_backend.h) so only the new
// cuda/qm_gpu_ri.cu and the core see it; the heavy ERI kernel TUs are not recompiled.

#include <string>

namespace qmint {

/// Device-resident fitted 3-index tensor T (nfit x n(n+1)/2: row Q = packed lower triangle of the
/// symmetric n x n block T_Q, column pairIndex(m, k) = m(m+1)/2 + k for m >= k, see qm_ri.h; packed since
/// Oct 2026) and the contractions of one RI J/K build (cuBLAS, FP64). The integrals, the metric and the
/// factorisation stay on the host (qmri::RIJK::setup); per build the host also forms the density
/// factor V (eigenvectors of P, n^3, or sqrt(2) C_occ) and, in mode j, c = J^-1 g (N^2). Every call
/// returns false on a device error; the caller then runs the CPU contraction.
class QMGpuRI {
public:
    virtual ~QMGpuRI() = default;
    virtual std::string name() const = 0;
    /// g_Q = sum_p T_Q,p w_p with the packed contraction weights w of P (qmri::packContractionWeights:
    /// P_mm, P_mk + P_km), i.e. sum_mk T_Q,mk P_mk  (w: n(n+1)/2; g: nfit)
    virtual bool contract(const double* w, double* g) = 0;
    /// J_p = sum_Q T_Q,p c_Q  (J: packed lower triangle, n(n+1)/2)
    virtual bool coulomb(const double* c, double* J) = 0;
    /// K = sum_Q (T_Q V+)(T_Q V+)^T - (T_Q V-)(T_Q V-)^T; V is n x r COLUMN-major, its first npos columns
    /// carry the positive eigenvalues of P (K: n x n full, row-major). T_Q is unpacked on the device per chunk.
    virtual bool exchange(const double* V, int r, int npos, double* K) = 0;
};

} // namespace qmint
