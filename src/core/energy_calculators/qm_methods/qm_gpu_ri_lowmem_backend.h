/*
 * < Device seam of the native QM engine: resource-saving density fitting (WP10) >
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0
 */

#pragma once

// Claude Generated (Oct 2026, WP10 "RI-JK ressourcensparend"). Own header so that only the new
// cuda/qm_gpu_ri_lowmem.cu and the core see it; qm_gpu_ri.cu and the ERI kernel units are not recompiled.

#include "qm_gpu_ri_backend.h"

#include <string>

namespace qmint {

/// Fitted RI tensor (nfit x n(n+1)/2, packed rows as in QMGpuRI) either RESIDENT on the device or STREAMED
/// from host memory in row chunks per call (-qm.gpu_ri_stream), in FP64 or FP32 (-qm.ri_mixed_precision; FP32
/// GEMM/SYRK per chunk, the chunk results accumulated in FP64). The host array must stay valid and unchanged for
/// the lifetime of a streamed object (it is page-locked with cudaHostRegister while the object lives).
class QMGpuRILowMem : public QMGpuRI {
public:
    /// One RI-JK build in ONE pass over the tensor (B already fitted, so c = g):
    ///   g_Q = T_Q . w,   J_p = sum_Q g_Q T_Q,p,   K = sum_Q (T_Q V+)(T_Q V+)^T - (T_Q V-)(T_Q V-)^T
    /// w: packed weights of P (n(n+1)/2); V: n x r column-major, first npos columns positive; outputs g (nfit),
    /// Jp (packed, n(n+1)/2), K (n x n full, row-major).
    virtual bool buildJK(const double* w, const double* V, int r, int npos, double* g, double* Jp, double* K) = 0;
    virtual bool streamed() const = 0;
    virtual bool fp32() const = 0;
    /// FP32 tensor only: true = rows widened to FP64 on the device and DGEMM/DSYRK (second mixed-precision phase),
    /// false = SGEMM/SSYRK (first phase, default). No effect on an FP64 tensor.
    virtual void setFP64Arithmetic(bool on) = 0;
    /// Accumulated copy / compute time and traffic since creation (diagnostic, one line).
    virtual std::string stats() const = 0;
};

/// Device part of the batched (no stored tensor) RI-JK exchange: X = L^-1 Y on the auxiliary index and
/// K += sign sum_P X^P X^P^T. L (naux x naux lower Cholesky factor of the Coulomb-normalised metric) stays on the
/// device; Y is streamed in column chunks (the solve acts on every column independently). FP64.
class QMGpuRIYSolver {
public:
    virtual ~QMGpuRIYSolver() = default;
    virtual std::string name() const = 0;
    /// Y: host row-major naux x (n rb); row P is the column-major n x rb block Y^P = T^P V_b. Adds
    /// sign * sum_P X^P (X^P)^T into K (n x n, full, row-major; only accumulated, never cleared).
    virtual bool exchange(const double* Y, int rb, double sign, double* K) = 0;
    virtual std::string stats() const = 0;
};

} // namespace qmint
