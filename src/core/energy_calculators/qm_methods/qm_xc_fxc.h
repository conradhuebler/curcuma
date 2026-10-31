/*
 * <Native QM: the XC parts of the double-hybrid gradient on the CPU (F_xc[D] and the P2-directional XC gradient,
 *  WP12 stage 2)>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * The second functional derivatives (jets of the kernels B88_X, LYP_C, PBE_X, P86_C, P86VWN_C, B95_C, the reason
 * for that second copy of the kernel formulas, and evaluatePointFxc) live in qm_xc_fxc_kernels.h, shared with the
 * device (cuda/qm_gpu_xc.cu, WP12 stage 3). This header adds the host-only parts: functionalCovered and the CPU
 * builders of qm_xc_response.cpp.
 *
 * Claude Generated (Oct 2026, WP12 stage 2, docs/DFT_ROADMAP/WP12_doppelhybride.md). AI-generated, machine-tested;
 * human production test pending.
 *
 * This program is free software under GPL-3.0
 */

#pragma once

#include "src/core/global.h"

#include "qm_grid.h"
#include "qm_xc.h"
#include "qm_xc_fxc_kernels.h"
#include "qm_xc_kernels.h"

#include <cmath>
#include <string>
#include <vector>

namespace qmxc {
namespace fxc {

/// True when every weighted term of f has a jet (the double-hybrid gradient refuses otherwise).
inline bool functionalCovered(const XCFunctional& f, std::string* missing = nullptr)
{
    for (const XCTerm& t : f.terms) {
        if (t.weight == 0.0) continue;
        Jet e;
        if (!kernelJet(t.kernel, 0.1, 0.01, 0.1, e)) {
            if (missing) *missing = kernelName(t.kernel);
            return false;
        }
    }
    return true;
}

}  // namespace fxc

/**
 * @brief Linear response of V_xc: F_xc[D]_mn = d/dt V_xc[P + t D]_mn at t = 0 (active AO basis, D symmetric).
 *
 *   F_xc[D]_mn = sum_g w_g [ dv_rho phi_m phi_n + 2 (dv_sigma grad rho + v_sigma grad rho_D) . grad(phi_m phi_n)
 *                + 1/2 dv_tau grad phi_m . grad phi_n ],
 *   dv_x = sum_y d^2e/dx dy y_D,  rho_D = sum D phi phi,  sigma_D = 2 grad rho . grad rho_D,  tau_D = 1/2 sum D grad phi grad phi.
 * Blocks as buildVxc (AO screen if active, else 128-point blocks), 32 fixed chunks (1 vs N threads bitwise).
 * Returns false (V untouched) if a kernel has no jet. Claude Generated (Oct 2026, WP12 stage 2).
 */
bool buildFxc(const qmgrid::MolecularGrid& grid, const qmgrid::AOEvaluator& ao, const Matrix& P, const Matrix& D,
              const XCFunctional& f, double rho_threshold, int threads, Matrix& V, const qmgrid::BlockScreen* screen);

/**
 * @brief Directional derivative, in direction P2, of the XC nuclear gradient of buildXCGradient at density P:
 *   out = d/dt [dE_xc[P + t P2]/dR]_{t=0} = the nuclear derivative of sum_mn P2_mn V_xc,mn[P] at fixed AO densities.
 *
 * Same parts as buildXCGradient: `basis`, `points` (points moving with their atom) and `edens` = the energy density of
 * the weight derivative, here the directional derivative v_rho rho2 + v_sigma sigma2 + v_tau tau2 (qmgrid::gridWeightGradient
 * gives the Becke weight part from it). Per point and function (T_k summed as in buildXCGradient, dE/dA_k = -2 sum T_k):
 *   T_k = d_k phi o (Z P2 + dZ P) + sum_j (v_j P2phi + dv_j Pphi) o d_k d_j phi + sum_j (wt dP2_j + dwt dP_j) o d_k d_j phi,
 *   Z = w v_rho phi + sum_j v_j d_j phi, v_j = 2 w v_sigma d_j rho, wt = w v_tau / 2, and the d-quantities their
 *   derivatives in direction P2 (dv_j = 2 w (dv_sigma d_j rho + v_sigma d_j rho2)).
 * Only second functional derivatives enter (WP12 doc, stage 2). Claude Generated (Oct 2026, WP12 stage 2).
 */
bool buildXCResponseGradient(const qmgrid::MolecularGrid& grid, const qmgrid::AOEvaluator& ao, const Matrix& P,
                             const Matrix& P2, const XCFunctional& f, double rho_threshold, int threads,
                             const std::vector<int>& bf_atom, int natoms, XCGradient& out,
                             const qmgrid::BlockScreen* screen);

}  // namespace qmxc
