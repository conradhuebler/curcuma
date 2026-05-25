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
 *  Parametrized overload takes q_sh / q_at explicitly.              *
 * ------------------------------------------------------------------ */
void XTB::addThirdOrderPotential(Potential& pot,
                                  const Vector& q_sh,
                                  const Vector& q_at) const
{
    if (m_method == MethodType::GFN1) {
        Vector v_at = coulomb::potential_third_order_atom(m_basis.z, q_at);
        pot.v_at += v_at;
        for (int s = 0; s < m_basis.nsh; ++s)
            pot.v_sh(s) += v_at(m_basis.sh2at[s]);
    } else {
        pot.v_sh += coulomb::potential_third_order_shell(
            m_basis.z, m_basis.sh2at, m_basis.ang_sh, q_sh);
    }
}

void XTB::addThirdOrderPotential(Potential& pot) const
{
    addThirdOrderPotential(pot, m_wfn.q_sh, m_wfn.q_at);
}

/* ------------------------------------------------------------------ *
 *  Linearized third-order Jacobian diagonal d(v_sh)/d(q_sh):        *
 *    GFN2: diag_s = 2·Γ_s·q_sh_s    (shell-resolved)                *
 *    GFN1: diag_s = 2·Γ_{at(s)}·q_at_{at(s)} (broadcast to shells)  *
 *  Used in CPSCF orbital-Hessian kernel K·z.                        *
 * ------------------------------------------------------------------ */
Vector XTB::thirdOrderKernelDiag(const Vector& q_sh, const Vector& q_at) const
{
    const int nsh = m_basis.nsh;
    Vector diag(nsh);
    if (m_method == MethodType::GFN1) {
        for (int s = 0; s < nsh; ++s) {
            const int iat = m_basis.sh2at[s];
            diag(s) = 2.0 * coulomb::atom_hubbard_deriv_gfn1(m_basis.z[iat]) * q_at(iat);
        }
    } else {
        for (int s = 0; s < nsh; ++s) {
            const int z = m_basis.z[m_basis.sh2at[s]];
            const int l = m_basis.ang_sh[s];
            diag(s) = 2.0 * coulomb::shell_hubbard_deriv_gfn2(z, l) * q_sh(s);
        }
    }
    return diag;
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
