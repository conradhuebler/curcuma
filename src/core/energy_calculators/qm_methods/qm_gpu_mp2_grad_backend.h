/*
 * < Device seam of the native QM engine: two-electron derivative contraction of the MP2 gradient >
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0
 */

#pragma once

// Claude Generated (Oct 2026, MP2 gradient). Seen only by cuda/qm_gpu_mp2_grad.cu and the core (loader
// gpu_plugin.{h,cpp}); the tables and QMGpuGradExtra are those of the RHF gradient (Q-G7).

#include "qm_gpu_grad_backend.h"

#include <string>

namespace qmint {

/// Two-electron part of the MP2 nuclear gradient on a device:
///   g_X = 1/2 sum_{mnls} D_mnls d(mn|ls)/dX,
///   D = P P - 1/4 (P P + P P) + (P2 P + P P2) - 1/4 (P2 P + P P2 + P2 P + P P2) + G
/// (index pattern as qmint::gradientTwoElectronMP2), Eh/Bohr, canonical quartets, translational invariance,
/// every quartet (no screening). P, P2: cartesian n x n; G: cartesian n^4 (8-fold symmetric), uploaded per call.
class QMGpuMP2Grad2e {
public:
    virtual ~QMGpuMP2Grad2e() = default;
    virtual std::string name() const = 0;
    /// grad: natoms x 3 row-major, overwritten. False on a device error (caller: CPU path).
    virtual bool compute(const double* P, const double* P2, const double* G, double* grad) = 0;
    virtual long long lastQuartets() const { return 0; }
    virtual std::string stats() const { return std::string(); }
};

} // namespace qmint
