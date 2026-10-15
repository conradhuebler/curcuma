/*
 * < Device seam of the native QM engine (J/K builds) >
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0
 */

#pragma once

#include <string>

// Claude Generated (Oct 2026, Q-G1 of docs/QM_GPU_ROADMAP.md): backend-free interface the
// QMEngine drives a GPU through. The concrete class lives in a runtime plugin
// (libcurcuma_cuda.so, entry curcuma_cuda_create_qm_jk), so the core carries no CUDA symbol.
namespace qmint {

/// Coulomb/exchange builder holding the two-electron integrals on a device.
/// Same definitions as buildCoulomb/buildExchange in qm_integrals.hpp:
///   J_mn = sum_ls P_ls (mn|ls),  K_mn = sum_ls P_ls (ml|ns)
/// P, J, K are n x n row-major, contiguous.
class QMGpuJK {
public:
    virtual ~QMGpuJK() = default;
    /// Returns false on a device error; the caller then falls back to the CPU.
    virtual bool build(const double* P, double* J, double* K) = 0;
    virtual std::string name() const = 0;
    virtual int n() const = 0;
};

} // namespace qmint
