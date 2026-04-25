/*
 * <xTB isotropic Coulomb kernel — gamma matrix wrappers>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Thin member-function wrappers around the header-only xtb_coulomb.hpp.
 *
 * Claude Generated. GPL-3.0.
 */

#include "xtb_native.h"
#include "xtb_coulomb.hpp"

#include "parameters/xtb_params_extra.hpp"

namespace curcuma::xtb {

/* ------------------------------------------------------------------ *
 *  Build the shell-resolved gamma matrix from current geometry.      *
 *  Called once per SCF cycle (or once per fixed geometry).           *
 * ------------------------------------------------------------------ */
void XTB::buildGammaMatrix()
{
    const double aa = AA_TO_AU;
    const int nat = m_atomcount;
    std::vector<double> xyz_bohr(3 * nat);
    for (int i = 0; i < nat; ++i) {
        xyz_bohr[3 * i + 0] = m_geometry(i, 0) * aa;
        xyz_bohr[3 * i + 1] = m_geometry(i, 1) * aa;
        xyz_bohr[3 * i + 2] = m_geometry(i, 2) * aa;
    }

    std::vector<int> z_atoms(m_atoms.begin(), m_atoms.end());

    const auto meth = (m_method == MethodType::GFN1)
        ? coulomb::Method::GFN1 : coulomb::Method::GFN2;

    m_gamma = coulomb::build_gamma_matrix(
        meth, z_atoms, xyz_bohr,
        m_basis.sh2at, m_basis.ang_sh);
}

/* ------------------------------------------------------------------ *
 *  Add isotropic Coulomb potential to v_sh: v_sh += γ * q_sh         *
 * ------------------------------------------------------------------ */
void XTB::addCoulombShellPotential(Potential& pot) const
{
    pot.v_sh += coulomb::potential_shell(m_gamma, m_wfn.q_sh);
}

/* ------------------------------------------------------------------ *
 *  Isotropic second-order Coulomb energy: E_iso = 0.5 * q_sh^T · γ · q_sh
 * ------------------------------------------------------------------ */
double XTB::energyCoulombShell() const
{
    return coulomb::energy_iso(m_gamma, m_wfn.q_sh);
}

} // namespace curcuma::xtb
