/*
 * < Device seam of the native QM engine: two-electron part of the RHF gradient (Q-G7) >
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0
 */

#pragma once

// Claude Generated (Oct 2026, Q-G7 of docs/QM_GPU_ROADMAP.md). Separate from
// qm_gpu_backend.h on purpose: the heavy ERI kernel TUs (cuda/qm_gpu_eri.cu, cuda/qm_gpu_direct.cu)
// include that header and take ~9 min each to compile; this one is only seen by
// cuda/qm_gpu_grad.cu and the core.

#include "qm_gpu_backend.h"

#include <string>

namespace qmint {

/// What the gradient kernel needs on top of QMGpuEriTables. The tables themselves come from
/// qmint::buildEriGradGpuTables: CARTESIAN basis (no 6d->5d transform, n_out = number of
/// cartesian functions), no quartet list (pair classes only), and the Hermite tables E of every
/// primitive pair built ONE angular-momentum step higher on both shells (layout
/// (LA+LB+3)(LA+2)(LB+2) per direction), which the centre derivative
///   d/dA_x g_l = 2 alpha g_{l+1} - l g_{l-1}
/// needs. All arrays are host memory; the plugin copies them.
struct QMGpuGradExtra {
    int natoms = 0;
    const int* sh_atom = nullptr;    ///< per shell: atom index (0..natoms-1)
    int nprimpair = 0;
    const double* pp_ea = nullptr;   ///< per primitive pair: exponent of the bra-side (A) primitive
    const double* pp_eb = nullptr;   ///< per primitive pair: exponent of the ket-side (B) primitive
};

/// Two-electron part of the closed-shell RHF nuclear gradient on a device, integral-direct:
///   g_X = 1/2 sum_{mnls} D_mnls d(mn|ls)/dX,  D_mnls = P_mn P_ls - kx/4 (P_ml P_ns + P_ms P_nl),
/// the same definition, units (Eh/Bohr) and canonical-quartet / translational-invariance scheme
/// as qmint::gradientTwoElectron (qm_integrals.cpp). Device memory O(n^2) (P, Pmax) plus the
/// pair tables; no ERI or derivative-integral tensor is ever stored.
class QMGpuGrad2e {
public:
    virtual ~QMGpuGrad2e() = default;
    virtual std::string name() const = 0;
    /// P: cartesian spin-summed density (ncart x ncart, row-major). grad: natoms x 3 row-major,
    /// overwritten. `threshold`: density-weighted quartet screen (a bound on max |D| over the quartet),
    /// skip a canonical shell quartet if
    ///   Q_AB Q_CD (Pmax_AB Pmax_CD + 1/4 (Pmax_AC Pmax_BD + Pmax_AD Pmax_BC)) < threshold
    /// (0 = every quartet, as the CPU reference). False on a device error (caller: CPU fallback).
    /// kx (Claude Generated, Oct 2026, GPU KS-DFT): weight of the exchange part, D = P_mn P_ls - kx/4 (...), as
    /// qmint::gradientTwoElectronScaled; kx = 1 is the RHF expression (same code path), kx = 0 (LDA/PBE) drops the
    /// exchange products and their screening term. Entry point curcuma_<backend>_create_qm_grad2e_kx.
    virtual bool compute(const double* P, double threshold, double kx, double* grad) = 0;
    /// Shell quartets evaluated in the last compute().
    virtual long long lastQuartets() const { return 0; }
    /// One ASCII line of device timings, empty if none.
    virtual std::string stats() const { return std::string(); }
};

} // namespace qmint
