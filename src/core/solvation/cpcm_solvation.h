/*
 * <Native CPCM (ddCOSMO) implicit solvation model>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Claude Generated (June 2026): CPCM continuum-solvation model for the native
 * GFN1/GFN2/GFN-FF methods, a faithful port of tblite's CPCM
 * (external/tblite/src/tblite/solvation/cpcm.f90 + cpcm_dd.f90).
 *
 * tblite's CPCM is a domain-decomposition ddCOSMO model whose reaction field is
 * a LINEAR operator on the atomic charges. We exploit this: @ref update solves
 * the forward + adjoint ddCOSMO systems once per unit charge to assemble the
 * effective symmetric reaction matrix M (nat x nat) with v = M q and
 * E = 1/2 q^T M q (= tblite get_potential / get_energy). Downstream this behaves
 * exactly like the ALPB Born matrix, so the in-SCF potential, GFN-FF EEQ
 * coupling (A_eeq += M) and the GPU device path are reused unchanged. The
 * nuclear gradient is ported directly from tblite get_gradient (it needs the
 * geometry derivatives of the cavity, not expressible through M alone).
 *
 * CPCM is purely electrostatic (no CDS / state shift); only the dielectric
 * constant and the COSMO van-der-Waals radii are needed.
 */

#pragma once

#include "src/core/global.h"
#include "src/core/solvation/cpcm_dd.h"
#include "src/core/solvation/implicit_solvation.h"

#include <Eigen/Dense>
#include <string>
#include <vector>

namespace Curcuma {
namespace Solvation {

/// CPCM (ddCOSMO) implicit solvation model. See file header.
class CpcmSolvation : public ImplicitSolvationModel {
public:
    CpcmSolvation() = default;
    ~CpcmSolvation() override = default;

    /// Explicit dielectric constant; <= 0 (default) derives it from the solvent
    /// name (tblite epsv). Must be set before @ref init.
    void setEpsilon(double eps) { m_epsilon = eps; }

    bool init(const std::vector<int>& atomic_numbers, const std::string& solvent,
              const std::string& method = "gfn2") override;

    void update(const std::vector<int>& atomic_numbers, const Matrix& xyz_bohr) override;

    void addPotential(const Vector& q_at, Vector& v_at) const override;

    double energy(const Vector& q_at) const override;

    void addGradient(const std::vector<int>& atomic_numbers,
                     const Matrix& xyz_bohr,
                     const Vector& q_at,
                     Matrix& gradient) override;

    /// The effective symmetric reaction matrix M (== Born matrix analogue).
    const Eigen::MatrixXd& reactionMatrix() const override { return m_M; }

    /// CPCM is purely electrostatic: gborn = 1/2 q^T M q, all other parts zero.
    ALPBEnergyParts energyParts(const Vector& q_at) const override;

    /// GPU hook: M is contracted against the plain Mulliken charges (no CM5 for
    /// either GFN1 or GFN2 under CPCM), so the device path can use it directly.
    const double* deviceBornMatrix() const override {
        return m_initialized ? m_M.data() : nullptr;
    }

    bool isInitialized() const { return m_initialized; }

private:
    bool m_initialized = false;
    int m_nat = 0;
    double m_epsilon = -1.0;  ///< explicit dielectric; <=0 -> from solvent name
    double m_keps = 0.0;      ///< -1/2 (1/eps - 1) / (1 + alpha_alpb)
    std::string m_solvent = "none";
    std::string m_method = "gfn2";

    DomainDecomposition m_dd;
    Eigen::MatrixXd m_jmat;   ///< ncav x nat Coulomb matrix (cavity points vs atoms)
    Eigen::MatrixXd m_M;      ///< nat x nat effective reaction matrix
};

} // namespace Solvation
} // namespace Curcuma
