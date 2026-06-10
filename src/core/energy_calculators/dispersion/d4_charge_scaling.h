/*
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Shared D4 charge-scaling function (the "zeta" of Caldeweyher 2019).
 * Exact port of dftd4 `model.f90::zeta`/`dzeta` (release_tblite/_deps/dftd4-src).
 *
 * One source of truth for the charge-dependent C6 scaling, used by
 *   - native GFN2 D4 : per-reference weighting, gi = eta·gc (gc=2)
 *   - native GFN-FF D4: single per-atom prefactor, gi = eta (gc=1)
 * so both methods share identical math. Pure functions (no table dependency);
 * the caller supplies a = ga (height, 3.0), c = gi (eta·gc), qref, qmod.
 *
 * Claude Generated (AP6b exact D4 port, 2026).
 */
#pragma once

#include <cmath>

namespace curcuma {
namespace dispersion {

/// D4 charge-scaling factor. dftd4 model.f90:725.
///   zeta = exp(a·(1 - exp(c·(1 - qref/qmod))))   for qmod >= 0
///        = exp(a)                                 otherwise (clamped, anionic)
inline double d4_zeta(double a, double c, double qref, double qmod)
{
    if (qmod < 0.0)
        return std::exp(a);
    return std::exp(a * (1.0 - std::exp(c * (1.0 - qref / qmod))));
}

/// First derivative d(zeta)/d(qmod) = d(zeta)/d(q). dftd4 model.f90:743.
///   dzeta = -a·c·exp(c·(1 - qref/qmod))·zeta·qref/qmod²
inline double d4_dzeta(double a, double c, double qref, double qmod)
{
    if (qmod < 0.0)
        return 0.0;
    const double g = std::exp(c * (1.0 - qref / qmod));
    const double z = std::exp(a * (1.0 - g));
    return -a * c * g * z * qref / (qmod * qmod);
}

/// Second derivative d²(zeta)/d(qmod)². Needed for the self-consistent D4
/// CPSCF response kernel (∂²E_D4/∂q²). Derived from d4_dzeta by the product
/// rule; let u = c·qref/qmod², g = exp(c·(1-qref/qmod)), z = exp(a·(1-g)):
///   dzeta = -a·g·z·u
///   d(dzeta)/dq = -a·[ g'·z·u + g·z'·u + g·z·u' ]
/// with g' = g·(c·qref/qmod²) = g·u, z' = z·(-a·g·u) = -a·g·z·u,
///      u' = -2·c·qref/qmod³ = -2u/qmod.
inline double d4_d2zeta(double a, double c, double qref, double qmod)
{
    if (qmod < 0.0)
        return 0.0;
    const double g = std::exp(c * (1.0 - qref / qmod));
    const double z = std::exp(a * (1.0 - g));
    const double u = c * qref / (qmod * qmod);   // = -d/dq [c·(1-qref/qmod)] sign handled below
    // dzeta = -a·g·z·u
    const double gp = g * u;                     // dg/dq
    const double zp = -a * gp * z;               // dz/dq = z·(-a·dg/dq)
    const double up = -2.0 * u / qmod;           // du/dq
    return -a * (gp * z * u + g * zp * u + g * z * up);
}

} // namespace dispersion
} // namespace curcuma
