/*
 * <xTB third-order onsite correction — wrappers>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * GFN1: atom-resolved third-order.
 * GFN2: shell-resolved third-order.
 * Both port the logic in xtb_coulomb.hpp (coulomb/thirdorder.f90).
 *
 * Claude Generated. GPL-3.0.
 */

#include "xtb_native.h"
#include "xtb_coulomb.hpp"

namespace curcuma::xtb {

/* ------------------------------------------------------------------ *
 *  Add third-order potential.                                        *
 *    GFN1: v_at(i) += q_i^2 · Γ_i   (atom-resolved)                 *
 *    GFN2: v_sh(s) += qsh_s^2 · Γ_s (shell-resolved)                *
 * ------------------------------------------------------------------ */
void XTB::addThirdOrderPotential(Potential& pot) const
{
    const auto meth = (m_method == MethodType::GFN1)
        ? coulomb::Method::GFN1 : coulomb::Method::GFN2;
    const int nat = m_atomcount;

    if (m_method == MethodType::GFN1) {
        // GFN1: atom-resolved — broadcast to shells
        Vector v_at = coulomb::potential_third_order_atom(m_basis.z, m_wfn.q_at);
        pot.v_at += v_at;
        for (int s = 0; s < m_basis.nsh; ++s) {
            const int iat = m_basis.sh2at[s];
            pot.v_sh(s) += v_at(iat);
        }
    } else {
        // GFN2: shell-resolved
        pot.v_sh += coulomb::potential_third_order_shell(
            m_basis.z, m_basis.sh2at, m_basis.ang_sh, m_wfn.q_sh);
    }
}

/* ------------------------------------------------------------------ *
 *  Third-order energy.                                               *
 * ------------------------------------------------------------------ */
double XTB::energyThirdOrder() const
{
    return coulomb::energy_third_order(
        (m_method == MethodType::GFN1) ? coulomb::Method::GFN1 : coulomb::Method::GFN2,
        m_basis.z, m_basis.sh2at, m_basis.ang_sh, m_wfn.q_at, m_wfn.q_sh);
}

} // namespace curcuma::xtb
