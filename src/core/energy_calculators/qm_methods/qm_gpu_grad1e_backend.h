/*
 * < Device seam of the native QM engine: one-electron part of the RHF gradient (Q-G8) >
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0
 */

#pragma once

// Claude Generated (Oct 2026, Q-G8 of docs/QM_GPU_ROADMAP.md). Own header, like
// qm_gpu_grad_backend.h, so that the ~9 min ERI kernel TUs (which include qm_gpu_backend.h) need
// not recompile when it changes; seen only by cuda/qm_gpu_grad1e.cu and the core.

#include "qm_gpu_grad_backend.h"

#include <string>

namespace qmint {

/// What the one-electron gradient kernel needs on top of the Q-G7 tables (QMGpuEriTables from
/// qmint::buildEriGradGpuTables, cartesian, Hermite tables raised by one on both shells, and
/// QMGpuGradExtra with the shell atoms and primitive exponents). All arrays are host memory.
struct QMGpuGrad1eExtra {
    QMGpuGradExtra base;              ///< natoms, sh_atom, pp_ea, pp_eb (as for the 2e gradient)
    const double* sh_xyz = nullptr;   ///< per shell: centre (Bohr), 3 per shell (= the basis-function centre)
    const double* atom_xyz = nullptr; ///< per atom: nuclear position (Bohr), 3 per atom
    const double* atom_Z = nullptr;   ///< per atom: nuclear charge of the attraction operator
};

/// One-electron part of the closed-shell RHF nuclear gradient on a device (Eh/Bohr):
///   g_X = sum_mn P_mn d(T_mn + V_mn)/dX - sum_mn W_mn dS_mn/dX,
/// with the basis-function (Pulay) parts on the function centres and the operator
/// (Hellmann-Feynman) part of V on the nucleus that moves; the same quantity, definitions and
/// units as qmint::gradientOneElectron (qm_integrals.cpp). The nuclear repulsion gradient is NOT
/// included (QMEngine::calculateCoreRepulsionGradient, as on the CPU).
class QMGpuGrad1e {
public:
    virtual ~QMGpuGrad1e() = default;
    virtual std::string name() const = 0;
    /// P, W: cartesian spin-summed density and energy-weighted density (ncart x ncart,
    /// row-major). grad: natoms x 3 row-major, overwritten. False on a device error (caller: CPU).
    virtual bool compute(const double* P, const double* W, double* grad) = 0;
    /// One ASCII line of device timings, empty if none.
    virtual std::string stats() const { return std::string(); }
};

} // namespace qmint
