/*
 * <RMSD-metadynamics core decisions for the strided scheme>
 * Copyright (C) 2019 - 2026 Conrad Huebler <Conrad.Huebler@gmx.net>
 *
 * Pure, side-effect-free decision helpers shared by the two RMSD-MTD code paths
 * (the local BiasThread path and the shared-pool path in SimpleMD::ApplyRMSDMTD).
 * Keeping the decisions here avoids duplicating them across the two mechanical
 * Kabsch loops. See docs/RMSD_MTD_TEXTBOOK.md sections 5.2-5.4 for the rationale.
 *
 * Claude Generated (Jul 2026)
 */

#pragma once

#include <algorithm>
#include <cmath>

namespace RMSDMTD {

// Auto hill spacing = FWHM of the hill exp(-alpha*RMSD^2): 2.35482 / sqrt(2*alpha).
// Neighbouring hills placed r_dep apart then overlap at ~half height. 0.5 A at alpha=10.
inline double autoRdep(double alpha)
{
    return 2.3548200450309493 / std::sqrt(2.0 * alpha);
}

// Deposition floor V_min: the bias one fresh unit-height hill produces at distance r_dep.
// Constant in pool size N (unlike the old current_bias*econv<N criterion).
inline double vMin(double k, double alpha, double r_dep)
{
    return k * std::exp(-alpha * r_dep * r_dep);
}

// Deposit a new hill when the pool is empty or the total bias at the walker is below V_min
// (i.e. no existing hill covers this region at spacing r_dep).
inline bool shouldDeposit(double bias, double v_min, int pool_count)
{
    return pool_count == 0 || bias < v_min;
}

// C1-continuous smoothstep ramp weight in [0,1] for the Milestone 2 held-force interpolation.
inline double smoothstep(double l)
{
    l = std::clamp(l, 0.0, 1.0);
    return l * l * (3.0 - 2.0 * l);
}

} // namespace RMSDMTD
