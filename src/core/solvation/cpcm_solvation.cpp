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
 * Claude Generated (June 2026): CPCM wrapper around the ddCOSMO core.
 */

#include "src/core/solvation/cpcm_solvation.h"

#include "src/core/curcuma_logger.h"
#include "src/core/solvation/cosmo_radii.h"
#include "src/core/solvation/lebedev_grid.h"
#include "src/core/solvation/tblite_solvation_params.h"

#include <cmath>

namespace Curcuma {
namespace Solvation {

namespace {

constexpr double kPi = 3.14159265358979323846;
const double kSqrt4pi = std::sqrt(4.0 * kPi);
constexpr double kAlphaAlpb = 0.571412; ///< CPCM/ALPB dielectric-scaling constant (tblite)

// Electric field of `nsrc` point charges `src` at `csrc` (3 x nsrc), evaluated
// at `ntrg` target points `ctrg` (3 x ntrg). Writes ef (3 x ntrg). tblite efld.
void efld(int nsrc, const double* src, const Eigen::MatrixXd& csrc,
          int ntrg, const Eigen::MatrixXd& ctrg, Eigen::MatrixXd& ef)
{
    ef.setZero();
    for (int j = 0; j < ntrg; ++j)
        for (int i = 0; i < nsrc; ++i) {
            const double vx = ctrg(0, j) - csrc(0, i);
            const double vy = ctrg(1, j) - csrc(1, i);
            const double vz = ctrg(2, j) - csrc(2, i);
            const double r2 = vx * vx + vy * vy + vz * vz;
            const double rr = std::sqrt(r2);
            const double r3 = r2 * rr;
            const double f = src[i] / r3;
            ef(0, j) += f * vx;
            ef(1, j) += f * vy;
            ef(2, j) += f * vz;
        }
}

} // namespace

bool CpcmSolvation::init(const std::vector<int>& atomic_numbers, const std::string& solvent,
                         const std::string& method)
{
    m_nat = static_cast<int>(atomic_numbers.size());
    m_solvent = solvent;
    m_method = method;

    // dielectric constant: explicit override, else tblite epsv for the solvent
    double eps = m_epsilon;
    if (eps <= 0.0) {
        const auto* p = curcuma::solvation::getTbliteSolvationParam(method, true, solvent);
        if (!p)
            p = curcuma::solvation::getTbliteSolvationParam(method, false, solvent);
        if (!p && method != "gfn2")
            p = curcuma::solvation::getTbliteSolvationParam("gfn2", true, solvent);
        if (!p) {
            CurcumaLogger::warn("CPCM solvation: no dielectric constant for solvent '" + solvent
                + "' (and no -xtb.solvent_epsilon given); solvation disabled.");
            return false;
        }
        eps = p->epsv;
    }
    if (eps <= 0.0)
        return false;

    // tblite keps = -1/2 (1/eps - 1) / (1 + alpha_alpb) (cpcm.f90:139)
    m_keps = -0.5 * (1.0 / eps - 1.0) / (1.0 + kAlphaAlpb);

    std::vector<double> rvdw(m_nat);
    for (int i = 0; i < m_nat; ++i)
        rvdw[i] = getVdwRadCosmo(atomic_numbers[i]); // Bohr

    // tblite cpcm_input defaults: nang = grid_size(6) = 74, lmax = 6, conv = 1e-8, eta = 2
    const int igrid = lebedevListBisection(74);
    bool ok = false;
    const std::vector<LebedevPoint> grid = lebedevAngularGrid(igrid, &ok);
    if (!ok || grid.empty()) {
        CurcumaLogger::warn("CPCM solvation: failed to build the Lebedev angular grid.");
        return false;
    }

    DomainDecompositionInput in;
    in.lmax = 6;
    // tblite's default cpcm conv is 1e-8, but it re-solves the surface every SCF
    // iteration with a warm restart, so its surface charges are refined far tighter
    // than a single cold solve. We build the reaction matrix M once from cold solves,
    // so we converge each one tighter (1e-10) to match tblite's effective accuracy.
    in.conv = 1e-10;
    in.eta = 2.0;
    m_dd.init(in, rvdw, grid);

    m_initialized = true;
    return true;
}

void CpcmSolvation::update(const std::vector<int>& /*atomic_numbers*/, const Matrix& xyz_bohr)
{
    if (!m_initialized)
        return;
    const int nat = m_nat;

    // copy the geometry into the dd column-major (3 x nat) layout
    DomainDecomposition::Mat ddxyz(3, nat);
    for (int i = 0; i < nat; ++i)
        for (int c = 0; c < 3; ++c)
            ddxyz(c, i) = xyz_bohr(i, c);
    m_dd.update(ddxyz);

    const int ncav = m_dd.ncav;
    const int nylm = m_dd.nylm;
    const auto& ccav = m_dd.cavityPoints(); // 3 x ncav

    // Coulomb matrix jmat(ic, j) = 1 / |ccav_ic - xyz_j| (tblite get_coulomb_matrix)
    m_jmat.resize(ncav, nat);
    for (int ic = 0; ic < ncav; ++ic)
        for (int j = 0; j < nat; ++j) {
            const double dx = ccav(0, ic) - ddxyz(0, j);
            const double dy = ccav(1, ic) - ddxyz(1, j);
            const double dz = ccav(2, ic) - ddxyz(2, j);
            m_jmat(ic, j) = 1.0 / std::sqrt(dx * dx + dy * dy + dz * dz);
        }

    // Assemble the effective reaction matrix M column-by-column from unit
    // charges: M(:, j) = tblite get_potential(e_j) (forward + adjoint ddCOSMO).
    m_M = Eigen::MatrixXd::Zero(nat, nat);
    DomainDecomposition::Mat sigma(nylm, nat), psi(nylm, nat), s(nylm, nat), dummy;
    std::vector<double> phi(ncav), zeta;
    for (int j = 0; j < nat; ++j) {
        for (int ic = 0; ic < ncav; ++ic)
            phi[ic] = m_jmat(ic, j);

        sigma.setZero();
        m_dd.solveCosmoDirect(true, phi, dummy, sigma, false);

        psi.setZero();
        psi(0, j) = kSqrt4pi;
        s.setZero();
        m_dd.solveCosmoAdjoint(psi, s, false);

        m_dd.getZeta(m_keps, s, zeta); // length ncav, includes keps

        for (int i = 0; i < nat; ++i) {
            double jtz = 0.0;
            for (int ic = 0; ic < ncav; ++ic)
                jtz += m_jmat(ic, i) * zeta[ic];
            m_M(i, j) = -jtz + m_keps * kSqrt4pi * sigma(0, i);
        }
    }
    // symmetrise (the ddCOSMO M is symmetric up to the solver tolerance)
    m_M = (0.5 * (m_M + m_M.transpose())).eval();
}

void CpcmSolvation::addPotential(const Vector& q_at, Vector& v_at) const
{
    if (!m_initialized)
        return;
    v_at.noalias() += m_M * q_at;
}

double CpcmSolvation::energy(const Vector& q_at) const
{
    if (!m_initialized)
        return 0.0;
    // Mirror tblite get_energy EXACTLY: forward-solve the surface charges for q and
    // E = keps * sum(sigma * psi) = keps * sqrt(4pi) * sum_i sigma(0,i) * q_i. This is
    // the forward-only contraction; using 1/2 q^T M q instead would expose the tiny
    // forward/adjoint inconsistency in the assembled M (grows to ~uEh on large polar
    // systems). The potential still uses M (= forward+adjoint), matching get_potential.
    const int nat = m_nat;
    const int ncav = m_dd.ncav;
    const int nylm = m_dd.nylm;

    std::vector<double> phi(ncav, 0.0);
    for (int ic = 0; ic < ncav; ++ic) {
        double sum = 0.0;
        for (int j = 0; j < nat; ++j)
            sum += m_jmat(ic, j) * q_at[j];
        phi[ic] = sum;
    }
    DomainDecomposition::Mat sigma(nylm, nat), dummy;
    sigma.setZero();
    m_dd.solveCosmoDirect(true, phi, dummy, sigma, false);

    double e = 0.0;
    for (int i = 0; i < nat; ++i)
        e += sigma(0, i) * q_at[i];
    return m_keps * kSqrt4pi * e;
}

ALPBEnergyParts CpcmSolvation::energyParts(const Vector& q_at) const
{
    ALPBEnergyParts parts;
    parts.gborn = energy(q_at);
    return parts;
}

void CpcmSolvation::addGradient(const std::vector<int>& /*atomic_numbers*/,
                                const Matrix& /*xyz_bohr*/,
                                const Vector& q_at,
                                Matrix& gradient)
{
    if (!m_initialized)
        return;
    const int nat = m_nat;
    const int ncav = m_dd.ncav;
    const int nylm = m_dd.nylm;
    const auto& ccav = m_dd.cavityPoints();  // 3 x ncav
    const auto& xyz = m_dd.atomCoords();      // 3 x nat
    const auto& cavAtom = m_dd.cavityAtom();

    std::vector<double> qv(nat);
    for (int i = 0; i < nat; ++i)
        qv[i] = q_at[i];

    // phi = jmat q  (potential at the cavity points from the atomic charges)
    std::vector<double> phi(ncav, 0.0);
    for (int ic = 0; ic < ncav; ++ic) {
        double sum = 0.0;
        for (int j = 0; j < nat; ++j)
            sum += m_jmat(ic, j) * qv[j];
        phi[ic] = sum;
    }

    // forward + adjoint solves at the actual charges (tblite get_gradient)
    DomainDecomposition::Mat sigma(nylm, nat), psi(nylm, nat), s(nylm, nat), dummy;
    sigma.setZero();
    m_dd.solveCosmoDirect(true, phi, dummy, sigma, false);
    psi.setZero();
    for (int iat = 0; iat < nat; ++iat)
        psi(0, iat) = kSqrt4pi * qv[iat];
    s.setZero();
    m_dd.solveCosmoAdjoint(psi, s, false);
    // refine the adjoint solution (tblite uses accuracy = conv * 1e-3 here); floor it
    // so the already-tight conv cannot drive the target below what the solver reaches.
    m_dd.solveCosmoAdjoint(psi, s, true, std::max(m_dd.conv * 1e-3, 1e-12));

    // ddCOSMO-specific force contribution
    DomainDecomposition::Mat gx(3, nat);
    m_dd.getDeriv(m_keps, phi, sigma, s, gx);

    std::vector<double> zeta;
    m_dd.getZeta(m_keps, s, zeta);

    // 1. solute electric field at the cavity points, contracted with zeta
    Eigen::MatrixXd ef1(3, ncav);
    efld(nat, qv.data(), xyz, ncav, ccav, ef1);
    for (int ii = 0; ii < ncav; ++ii) {
        const int iat = cavAtom[ii];
        gx(0, iat) += zeta[ii] * ef1(0, ii);
        gx(1, iat) += zeta[ii] * ef1(1, ii);
        gx(2, iat) += zeta[ii] * ef1(2, ii);
    }

    // 2. zeta's electric field at the nuclei, contracted with the charges
    Eigen::MatrixXd ef2(3, nat);
    efld(ncav, zeta.data(), ccav, nat, xyz, ef2);
    for (int iat = 0; iat < nat; ++iat) {
        gx(0, iat) += ef2(0, iat) * qv[iat];
        gx(1, iat) += ef2(1, iat) * qv[iat];
        gx(2, iat) += ef2(2, iat) * qv[iat];
    }

    // accumulate into the molecular gradient (Eh/Bohr, N x 3)
    for (int iat = 0; iat < nat; ++iat) {
        gradient(iat, 0) += gx(0, iat);
        gradient(iat, 1) += gx(1, iat);
        gradient(iat, 2) += gx(2, iat);
    }
}

} // namespace Solvation
} // namespace Curcuma
