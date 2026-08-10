/*
 * <QMDFF force-constant parametrisation from a QM Hessian>
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
 * Claude Generated (August 2026). See qmdff_parametrisation.h for the method,
 * the linearity argument and the literature references.
 */

#include "qmdff_parametrisation.h"

#include "gfnff_unit_terms.h"

#include "src/core/curcuma_logger.h"
#include "src/core/elements.h"

#include <fmt/core.h>
#include <fmt/format.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <set>

namespace curcuma::qmdff {

namespace {

/// Angstrom -> Bohr
constexpr double kAngstromToBohr = 1.889726125;

/**
 * @brief Central-difference Hessian block of a single term.
 *
 * @param n        number of atoms the term touches (2..4)
 * @param pos      the term's atom positions in Bohr (modified and restored)
 * @param h        displacement in Bohr
 * @param gradient callable (const std::vector<Vector3d>& pos, std::vector<Vector3d>& grad)
 *                 that fills `grad` with dE/dpos for the term with its constant set to 1
 * @param block    out: 12x12, only the leading 3n x 3n is written
 *
 * Differentiating the ALREADY-VALIDATED analytic term gradients (qmdff_terms.h, FD-checked
 * to <=7e-11 Eh/Angstrom by ctest qmdff_gradients) rather than deriving second derivatives
 * by hand: 6n kernel calls per term, and no new unverified math.
 */
template <typename GradFn>
void finiteDifferenceBlock(int n, std::vector<Eigen::Vector3d>& pos, double h,
    GradFn&& gradient, Eigen::Matrix<double, 12, 12>& block)
{
    block.setZero();
    std::vector<Eigen::Vector3d> gp(n), gm(n);

    for (int a = 0; a < n; ++a) {
        for (int al = 0; al < 3; ++al) {
            const double saved = pos[a](al);

            pos[a](al) = saved + h;
            gradient(pos, gp);
            pos[a](al) = saved - h;
            gradient(pos, gm);
            pos[a](al) = saved;

            for (int b = 0; b < n; ++b)
                for (int be = 0; be < 3; ++be)
                    block(3 * a + al, 3 * b + be) = (gp[b](be) - gm[b](be)) / (2.0 * h);
        }
    }

    // Symmetrise: the exact Hessian is symmetric, the finite difference is only nearly so.
    const int m = 3 * n;
    for (int p = 0; p < m; ++p)
        for (int q = p + 1; q < m; ++q) {
            const double s = 0.5 * (block(p, q) + block(q, p));
            block(p, q) = s;
            block(q, p) = s;
        }
}

/// Seminario projection of one 3x3 interaction block onto a direction.
/// Reference: Seminario, Int. J. Quantum Chem. 60 (1996) 1271, eq. 9.
double parameterSign(TermKind kind, TermSource source);

double seminarioProjection(const Eigen::Matrix3d& h_block, const Eigen::Vector3d& direction)
{
    // Modified Seminario (Allen/Payne/Cole, JCTC 14 (2018) 274): symmetrise first.
    const Eigen::Matrix3d k = -0.5 * (h_block + h_block.transpose());
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(k);
    if (solver.info() != Eigen::Success)
        return 0.0;

    double sum = 0.0;
    for (int a = 0; a < 3; ++a)
        sum += solver.eigenvalues()(a) * std::abs(direction.dot(solver.eigenvectors().col(a)));
    return sum;
}

/**
 * @brief Sign convention of a fitted constant.
 *
 * GFN-FF's bond term is E = fc * exp(-alpha dr^2), which is a well only for fc < 0 — the
 * constant IS the well depth. Every other constant here is a positive curvature or
 * barrier. The fit therefore works internally with s = sign * k >= 0 so that one
 * non-negativity constraint is correct for all of them, and undoes the flip on write-back.
 * Getting this wrong pins every GFN-FF bond at its floor (observed: 8 of 8 on CH3OCH3).
 */
double parameterSign(TermKind kind, TermSource source)
{
    return (source == TermSource::GFN_FF && kind == TermKind::Bond) ? -1.0 : 1.0;
}

} // anonymous namespace

const char* termKindName(TermKind kind)
{
    switch (kind) {
    case TermKind::Bond: return "bond";
    case TermKind::Angle: return "angle";
    case TermKind::Torsion: return "torsion";
    case TermKind::OutOfPlane: return "out-of-plane";
    case TermKind::ExtraTorsion: return "extra torsion";
    }
    return "unknown";
}

json FitReport::toJson() const
{
    json j;
    j["n_terms"] = n_terms;
    j["n_parameters"] = n_params;
    j["n_bonds"] = n_bonds;
    j["n_angles"] = n_angles;
    j["n_torsions"] = n_torsions;
    j["n_out_of_plane"] = n_oop;
    j["lambda_used"] = lambda_used;
    j["relative_residual_initial"] = rel_residual_initial;
    j["relative_residual_fitted"] = rel_residual_fitted;
    j["r_squared"] = r_squared;
    j["n_at_floor"] = n_at_floor;
    j["n_frozen"] = n_frozen;
    j["target_norm"] = target_norm;
    j["tr_content"] = tr_content;
    j["n_basins"] = n_basins;
    j["relative_residual_energy"] = rel_residual_energy;
    j["relative_residual_gradient"] = rel_residual_gradient;

    json ladder = json::array();
    for (const auto& [lam, res] : lambda_curve)
        ladder.push_back({ { "lambda", lam }, { "relative_residual", res } });
    j["lambda_ladder"] = ladder;

    json ki = json::array(), kf = json::array();
    for (int i = 0; i < k_initial.size(); ++i)
        ki.push_back(k_initial(i));
    for (int i = 0; i < k_fitted.size(); ++i)
        kf.push_back(k_fitted(i));
    j["k_initial"] = ki;
    j["k_fitted"] = kf;
    j["warnings"] = warnings;
    return j;
}

// =============================================================================
// Construction
// =============================================================================

QMDFFParametrisation::QMDFFParametrisation(const std::vector<int>& atoms,
    const Matrix& geometry_angstrom,
    const json& parameters,
    const FitOptions& options)
    : m_atoms(atoms)
    , m_parameters(parameters)
    , m_options(options)
    , m_natoms(static_cast<int>(atoms.size()))
{
    m_geometry_bohr = geometry_angstrom * kAngstromToBohr;

    m_weight = Vector::Ones(3 * m_natoms);
    if (m_options.mass_weighted) {
        for (int a = 0; a < m_natoms; ++a) {
            const double w = 1.0 / std::sqrt(Elements::AtomicMass[m_atoms[a]]);
            for (int c = 0; c < 3; ++c)
                m_weight(3 * a + c) = w;
        }
    }

    buildTermData(m_geometry_bohr, m_terms);
    assignParameterGroups();
    buildNonbonded(m_geometry_bohr, m_h_nonbonded, m_e_nonbonded, m_g_nonbonded);
}

Matrix QMDFFParametrisation::hessianAngstromToBohr(const Matrix& h)
{
    // d^2E/dx_Bohr^2 = d^2E/dx_Angstrom^2 / (Bohr per Angstrom)^2
    return h / (kAngstromToBohr * kAngstromToBohr);
}

void QMDFFParametrisation::setCNDerivatives(const Matrix& dcn)
{
    m_dcn = dcn;
    buildTermData(m_geometry_bohr, m_terms);
    assignParameterGroups();
}

void QMDFFParametrisation::setCoordinationNumbers(const Vector& cn)
{
    m_cn = cn;
    // The constructor built the terms before any coordination numbers existed, so the
    // GFN-FF bonds took the static r0_ij fallback instead of the dynamic
    // (r0_base + cnfak*cn + rabshift)*ff the force field actually evaluates. Rebuild.
    buildTermData(m_geometry_bohr, m_terms);
    assignParameterGroups();
}

void QMDFFParametrisation::buildTermData(const Matrix& geom_bohr,
    std::vector<TermHessian>& out_terms) const
{
    std::vector<TermHessian>& m_terms = out_terms;   // keep the body below unchanged
    m_terms.clear();
    const double h = m_options.fd_step_bohr;

    auto position = [&](int a) {
        return Eigen::Vector3d(geom_bohr(a, 0), geom_bohr(a, 1), geom_bohr(a, 2));
    };
    auto elementOf = [&](int a) {
        return (a >= 0 && a < m_natoms) ? m_atoms[a] : 1;
    };

    const bool gfnff = (m_options.source == TermSource::GFN_FF);

    // --- bonds ----------------------------------------------------------------
    // QMDFF : E = fc * f(r; r0, a)
    // GFN-FF: E = fc * exp(-alpha * dr^2), with the DYNAMIC r0 and the hydrogen-bridge
    //         alpha modulation reproduced exactly as FFWorkspace::calcBonds does them.
    const json bonds = m_parameters.value("bonds", json::array());
    for (int idx = 0; idx < static_cast<int>(bonds.size()); ++idx) {
        const auto& b = bonds[idx];
        TermHessian t;
        t.kind = TermKind::Bond;
        t.list_index = idx;
        t.natoms = 2;
        t.atom[0] = b["i"];
        t.atom[1] = b["j"];

        const double r0 = b.value("r0_ij", 0.0) * kAngstromToBohr;
        const double expo = b.value("exponent", 0.0);
        const auto potential = (b.value("qmdff_potential", 0) == 1)
            ? QMDFFTerms::BondPotential::Morse
            : QMDFFTerms::BondPotential::LennardJones;

        // GFN-FF bonds that take part in a hydrogen bridge have their exponent modulated by
        // alpha *= (1 - 0.1 hb_cn_H), and hb_cn_H is computed INSIDE FFWorkspace at
        // evaluation time (computeHBCoordinationNumbers) — it is 0 in the exported parameter
        // set and there is no accessor for the computed value. Fitting them would use the
        // wrong exponent, so they are left at their GFN-FF value instead; they remain fully
        // accounted for, because the non-fitted remainder is taken by difference.
        if (gfnff && b.value("nr_hb", 0) >= 1)
            continue;

        // GFN-FF: dynamic r0 from the coordination numbers
        double r0_gfnff = b.value("r0_ij", 0.0);
        double alpha_gfnff = expo;
        if (gfnff) {
            const int zi = b.value("z_i", 0);
            const int zj = b.value("z_j", 0);
            if (m_cn.size() > t.atom[0] && m_cn.size() > t.atom[1] && zi > 0 && zj > 0) {
                const double ra = b.value("r0_base_i", 0.0) + b.value("cnfak_i", 0.0) * m_cn(t.atom[0]);
                const double rb = b.value("r0_base_j", 0.0) + b.value("cnfak_j", 0.0) * m_cn(t.atom[1]);
                r0_gfnff = (ra + rb + b.value("rabshift", 0.0)) * b.value("ff", 1.0);
            }
            if (b.value("nr_hb", 0) >= 1) {
                constexpr double vbond_scale = 0.9;
                alpha_gfnff = (-(1.0 - vbond_scale) * b.value("hb_cn_H", 0.0) + 1.0) * expo;
            }
        }

        std::vector<Eigen::Vector3d> pos{ position(t.atom[0]), position(t.atom[1]) };
        finiteDifferenceBlock(2, pos, h,
            [&](const std::vector<Eigen::Vector3d>& p, std::vector<Eigen::Vector3d>& g) {
                g[0].setZero();
                g[1].setZero();
                if (gfnff)
                    GFNFFUnitTerms::bondEnergy(p[0], p[1], r0_gfnff, alpha_gfnff, g[0], g[1], true);
                else
                    QMDFFTerms::bondEnergy(p[0], p[1], r0, 1.0, expo, potential, g[0], g[1], true);
            },
            t.block);
        {   // energy and gradient at unit constant, same geometry as the block
            Eigen::Vector3d g0 = Eigen::Vector3d::Zero(), g1 = g0;
            t.energy = gfnff
                ? GFNFFUnitTerms::bondEnergy(pos[0], pos[1], r0_gfnff, alpha_gfnff, g0, g1, true)
                : QMDFFTerms::bondEnergy(pos[0], pos[1], r0, 1.0, expo, potential, g0, g1, true);
            t.grad.setZero();
            t.grad.segment<3>(0) = g0;
            t.grad.segment<3>(3) = g1;
        }

        // GFN-FF dynamic-r0 coordination-number response.
        //
        // FFWorkspace::calcBonds accumulates dE/dcn_i = -dE/dr * ff * cnfak_i and the
        // gradient assembly adds sum_a (dE/dcn_a) dcn_a/dx. With cn (and hence r0 and
        // dcn/dx) frozen at the reference geometry, the bond gradient is exactly
        //     g = dE/dr * (u - w),   u = dr/dx,   w = ff (cnfak_i dcn_i/dx + cnfak_j dcn_j/dx),
        // so the Hessian picks up  -d2E/dr2 * u w^T  on top of the local two-atom block.
        // The force field symmetrises its numerical Hessian, so the effective addition is
        // the symmetric part. Everything here is proportional to the force constant, which
        // is why the fit stays linear.
        if (gfnff && m_dcn.rows() == m_natoms && m_dcn.cols() == 3 * m_natoms) {
            const double ff = b.value("ff", 1.0);
            const double cnfak_i = b.value("cnfak_i", 0.0);
            const double cnfak_j = b.value("cnfak_j", 0.0);
            if (ff != 0.0 && (cnfak_i != 0.0 || cnfak_j != 0.0)) {
                const Eigen::Vector3d d = pos[0] - pos[1];
                const double r = d.norm();
                if (r > 1e-12) {
                    const double dr = r - r0_gfnff;
                    const double e = std::exp(-alpha_gfnff * dr * dr);
                    const double d2Edr2 = (4.0 * alpha_gfnff * alpha_gfnff * dr * dr
                                              - 2.0 * alpha_gfnff)
                        * e;
                    const double dEdr = -2.0 * alpha_gfnff * dr * e;

                    t.extra_u = Vector::Zero(3 * m_natoms);
                    t.extra_u.segment<3>(3 * t.atom[0]) = d / r;
                    t.extra_u.segment<3>(3 * t.atom[1]) = -d / r;

                    t.extra_w = ff
                        * (cnfak_i * m_dcn.row(t.atom[0]).transpose()
                            + cnfak_j * m_dcn.row(t.atom[1]).transpose());

                    t.extra_coeff = -0.5 * d2Edr2; // symmetrised -d2Edr2 * u w^T
                    t.extra_grad_coeff = -dEdr;    // g += (-dE/dr) * w
                }
            }
        }
        m_terms.push_back(t);
    }

    // --- angles: E = fc * (cos(theta) - cos(theta0))^2 * damp -----------------
    const json angles = m_parameters.value("angles", json::array());
    for (int idx = 0; idx < static_cast<int>(angles.size()); ++idx) {
        const auto& a = angles[idx];
        TermHessian t;
        t.kind = TermKind::Angle;
        t.list_index = idx;
        t.natoms = 3;
        t.atom[0] = a["i"];
        t.atom[1] = a["j"]; // centre, matching the Fortran angl(1,.) convention
        t.atom[2] = a["k"];

        const double theta0 = a.value("theta0_ijk", 0.0); // radians
        const int za = elementOf(t.atom[0]);
        const int zb = elementOf(t.atom[1]);
        const int zc = elementOf(t.atom[2]);

        std::vector<Eigen::Vector3d> pos{ position(t.atom[0]), position(t.atom[1]), position(t.atom[2]) };
        finiteDifferenceBlock(3, pos, h,
            [&](const std::vector<Eigen::Vector3d>& p, std::vector<Eigen::Vector3d>& g) {
                g[0].setZero();
                g[1].setZero();
                g[2].setZero();
                if (gfnff)
                    GFNFFUnitTerms::angleEnergy(p[0], p[1], p[2], za, zb, zc, theta0,
                        g[0], g[1], g[2], true);
                else
                    QMDFFTerms::angleEnergy(p[0], p[1], p[2], za, zb, zc, theta0, 1.0,
                        g[0], g[1], g[2], true);
            },
            t.block);
        {
            Eigen::Vector3d g0 = Eigen::Vector3d::Zero(), g1 = g0, g2 = g0;
            t.energy = gfnff
                ? GFNFFUnitTerms::angleEnergy(pos[0], pos[1], pos[2], za, zb, zc, theta0,
                    g0, g1, g2, true)
                : QMDFFTerms::angleEnergy(pos[0], pos[1], pos[2], za, zb, zc, theta0, 1.0,
                    g0, g1, g2, true);
            t.grad.setZero();
            t.grad.segment<3>(0) = g0;
            t.grad.segment<3>(3) = g1;
            t.grad.segment<3>(6) = g2;
        }
        m_terms.push_back(t);
    }

    // --- GFN-FF torsions and inversions ---------------------------------------
    // GFN-FF keeps them in separate lists (dihedrals / extra_dihedrals / inversions),
    // unlike QMDFF which puts out-of-plane terms into the torsion list.
    if (gfnff) {
        auto addTorsionList = [&](const char* key, TermKind kind) {
            int local_index = 0;
            for (const auto& d : m_parameters.value(key, json::array())) {
                if (!m_options.fit_torsions)
                    return;
                TermHessian t;
                t.kind = kind;
                t.list_index = local_index++;
                t.natoms = 4;
                t.atom[0] = d["i"]; t.atom[1] = d["j"]; t.atom[2] = d["k"]; t.atom[3] = d["l"];
                const double n = d.value("n", 0.0);
                const double phi0 = d.value("phi0", 0.0);
                const bool is_nci = d.value("is_nci", false);
                const int zi = elementOf(t.atom[0]), zj = elementOf(t.atom[1]);
                const int zk = elementOf(t.atom[2]), zl = elementOf(t.atom[3]);

                std::vector<Eigen::Vector3d> pos{ position(t.atom[0]), position(t.atom[1]),
                    position(t.atom[2]), position(t.atom[3]) };
                finiteDifferenceBlock(4, pos, h,
                    [&](const std::vector<Eigen::Vector3d>& p, std::vector<Eigen::Vector3d>& g) {
                        Eigen::Matrix<double, 4, 3> gm = Eigen::Matrix<double, 4, 3>::Zero();
                        GFNFFUnitTerms::torsionEnergy(p[0], p[1], p[2], p[3], zi, zj, zk, zl,
                            n, phi0, is_nci, gm, true);
                        for (int r = 0; r < 4; ++r)
                            g[r] = gm.row(r).transpose();
                    },
                    t.block);
                {
                    Eigen::Matrix<double, 4, 3> gm = Eigen::Matrix<double, 4, 3>::Zero();
                    t.energy = GFNFFUnitTerms::torsionEnergy(pos[0], pos[1], pos[2], pos[3],
                        zi, zj, zk, zl, n, phi0, is_nci, gm, true);
                    t.grad.setZero();
                    for (int r = 0; r < 4; ++r)
                        t.grad.segment<3>(3 * r) = gm.row(r).transpose();
                }
                m_terms.push_back(t);
            }
        };
        addTorsionList("dihedrals", TermKind::Torsion);
        addTorsionList("extra_dihedrals", TermKind::ExtraTorsion);

        if (m_options.fit_oop) {
            int inv_index = 0;
            for (const auto& inv : m_parameters.value("inversions", json::array())) {
                TermHessian t;
                t.kind = TermKind::OutOfPlane;
                t.list_index = inv_index++;
                t.natoms = 4;
                t.atom[0] = inv["i"]; t.atom[1] = inv["j"];
                t.atom[2] = inv["k"]; t.atom[3] = inv["l"];
                const int ptype = inv.value("potential_type", 0);
                const double omega0 = inv.value("omega0", 0.0);
                const int zc = elementOf(t.atom[0]), z1 = elementOf(t.atom[1]);
                const int z2 = elementOf(t.atom[2]), z3 = elementOf(t.atom[3]);

                std::vector<Eigen::Vector3d> pos{ position(t.atom[0]), position(t.atom[1]),
                    position(t.atom[2]), position(t.atom[3]) };
                finiteDifferenceBlock(4, pos, h,
                    [&](const std::vector<Eigen::Vector3d>& p, std::vector<Eigen::Vector3d>& g) {
                        Eigen::Matrix<double, 4, 3> gm = Eigen::Matrix<double, 4, 3>::Zero();
                        GFNFFUnitTerms::outOfPlaneEnergy(p[0], p[1], p[2], p[3],
                            zc, z1, z2, z3, ptype, omega0, gm, true);
                        for (int r = 0; r < 4; ++r)
                            g[r] = gm.row(r).transpose();
                    },
                    t.block);
                {
                    Eigen::Matrix<double, 4, 3> gm = Eigen::Matrix<double, 4, 3>::Zero();
                    t.energy = GFNFFUnitTerms::outOfPlaneEnergy(pos[0], pos[1], pos[2], pos[3],
                        zc, z1, z2, z3, ptype, omega0, gm, true);
                    t.grad.setZero();
                    for (int r = 0; r < 4; ++r)
                        t.grad.segment<3>(3 * r) = gm.row(r).transpose();
                }
                m_terms.push_back(t);
            }
        }
    }

    // --- QMDFF torsions and out-of-plane: E = scale * (...) * damp ------------
    const json torsions = gfnff ? json::array() : m_parameters.value("qmdff_torsions", json::array());
    for (int idx = 0; idx < static_cast<int>(torsions.size()); ++idx) {
        const auto& entry = torsions[idx];
        const bool oop = entry.value("out_of_plane", false);
        if (oop && !m_options.fit_oop)
            continue;
        if (!oop && !m_options.fit_torsions)
            continue;

        TermHessian t;
        t.kind = oop ? TermKind::OutOfPlane : TermKind::Torsion;
        t.list_index = idx;
        t.natoms = 4;
        t.atom[0] = entry["i"];
        t.atom[1] = entry["j"];
        t.atom[2] = entry["k"];
        t.atom[3] = entry["l"];

        // Unit term: same shape, scale set to 1
        QMDFFTerms::QMDFFTorsion tor;
        tor.i = t.atom[0];
        tor.j = t.atom[1];
        tor.k = t.atom[2];
        tor.l = t.atom[3];
        tor.out_of_plane = oop;
        tor.phi0 = entry.value("phi0", 0.0);
        tor.scale = 1.0;
        const json terms = entry.value("terms", json::array());
        tor.nterm = std::min<int>(static_cast<int>(terms.size()), QMDFFTerms::kMaxTorsionTerms);
        for (int n = 0; n < tor.nterm; ++n) {
            tor.rn[n] = terms[n].value("n", 0.0);
            tor.phase[n] = terms[n].value("phase", 0.0);
            tor.v[n] = terms[n].value("v", 0.0);
        }
        if (oop && tor.nterm == 0)
            tor.nterm = 1;

        const int zi = elementOf(t.atom[0]);
        const int zj = elementOf(t.atom[1]);
        const int zk = elementOf(t.atom[2]);
        const int zl = elementOf(t.atom[3]);

        std::vector<Eigen::Vector3d> pos{ position(t.atom[0]), position(t.atom[1]),
            position(t.atom[2]), position(t.atom[3]) };
        finiteDifferenceBlock(4, pos, h,
            [&](const std::vector<Eigen::Vector3d>& p, std::vector<Eigen::Vector3d>& g) {
                Eigen::Matrix<double, 4, 3> gm = Eigen::Matrix<double, 4, 3>::Zero();
                if (oop)
                    QMDFFTerms::outOfPlaneEnergy(p[0], p[1], p[2], p[3], zi, zj, zk, zl, tor, gm, true);
                else
                    QMDFFTerms::torsionEnergy(p[0], p[1], p[2], p[3], zi, zj, zk, zl, tor, gm, true);
                for (int r = 0; r < 4; ++r)
                    g[r] = gm.row(r).transpose();
            },
            t.block);
        {
            Eigen::Matrix<double, 4, 3> gm = Eigen::Matrix<double, 4, 3>::Zero();
            t.energy = oop
                ? QMDFFTerms::outOfPlaneEnergy(pos[0], pos[1], pos[2], pos[3], zi, zj, zk, zl, tor, gm, true)
                : QMDFFTerms::torsionEnergy(pos[0], pos[1], pos[2], pos[3], zi, zj, zk, zl, tor, gm, true);
            t.grad.setZero();
            for (int r = 0; r < 4; ++r)
                t.grad.segment<3>(3 * r) = gm.row(r).transpose();
        }
        m_terms.push_back(t);
    }

    // Apply the sign convention so the solver only ever sees non-negative constants.
    for (auto& t : m_terms) {
        const double sign = parameterSign(t.kind, m_options.source);
        if (sign < 0.0) {
            t.block = -t.block;
            t.energy = -t.energy;
            t.grad = -t.grad;
            t.extra_coeff = -t.extra_coeff;
            t.extra_grad_coeff = -t.extra_grad_coeff;
        }
    }

    // Cache <U,U> under the active inner product.
    const Vector mw = m_options.mass_weighted ? Vector(m_weight.cwiseProduct(m_weight)) : Vector();
    for (auto& t : m_terms)
        t.norm2 = t.dotTerm(t, mw);
}

void QMDFFParametrisation::assignParameterGroups()
{
    // Torsions about the same central bond have near-proportional unit Hessians, so their
    // individual scales are not identifiable — only their sum is. Tying them mirrors how
    // UFF and QMDFF both treat a rotatable bond and keeps the normal equations conditioned.
    std::map<std::pair<int, int>, int> central_bond_group;
    int next = 0;

    for (auto& t : m_terms) {
        // Tying is only valid when the tied terms share one constant. QMDFF's torsion
        // `scale` is 1 for all of them, so a group has a well-defined current value; GFN-FF
        // gives every dihedral its own barrier V, and collapsing those to a single number
        // corrupts both the starting guess and the remainder H_nb = H_FF - sum k U
        // (measured: a 1.9e-2 verification residual on a peptide). Fitting a shared
        // MULTIPLIER on the existing barriers would restore tying for GFN-FF; until then
        // each GFN-FF torsion is its own parameter, kept identifiable by lambda_torsion.
        if ((t.kind == TermKind::Torsion || t.kind == TermKind::ExtraTorsion)
            && m_options.tie_torsions_by_central_bond
            && m_options.source != TermSource::GFN_FF) {
            const std::pair<int, int> key{ std::min(t.atom[1], t.atom[2]),
                std::max(t.atom[1], t.atom[2]) };
            auto it = central_bond_group.find(key);
            if (it == central_bond_group.end()) {
                central_bond_group[key] = next;
                t.group = next;
                ++next;
            } else {
                t.group = it->second;
            }
        } else {
            t.group = next;
            ++next;
        }
    }
    m_nparams = next;

    // Lower bounds, per term kind
    m_floors = Vector::Zero(m_nparams);
    m_is_torsion_group.assign(m_nparams, 0);
    for (const auto& t : m_terms) {
        if (t.kind == TermKind::Torsion || t.kind == TermKind::ExtraTorsion
            || t.kind == TermKind::OutOfPlane)
            m_is_torsion_group[t.group] = 1;
        double floor_value = 0.0;
        switch (t.kind) {
        case TermKind::Bond: floor_value = m_options.bond_floor; break;
        case TermKind::Angle: floor_value = m_options.angle_floor; break;
        case TermKind::Torsion:
        case TermKind::ExtraTorsion: floor_value = m_options.torsion_floor; break;
        case TermKind::OutOfPlane: floor_value = m_options.oop_floor; break;
        }
        m_floors(t.group) = std::max(m_floors(t.group), floor_value);
    }
}

void QMDFFParametrisation::buildNonbonded(const Matrix& geom_bohr, Matrix& h_nb,
    double& e_nb, Matrix& g_nb) const
{
    const int n3 = 3 * m_natoms;
    h_nb = Matrix::Zero(n3, n3);
    g_nb = Matrix::Zero(m_natoms, 3);
    e_nb = 0.0;
    Matrix& m_h_nonbonded = h_nb;   // keep the scatter body below unchanged
    const double h = m_options.fd_step_bohr;

    auto position = [&](int a) {
        return Eigen::Vector3d(geom_bohr(a, 0), geom_bohr(a, 1), geom_bohr(a, 2));
    };
    auto scatter = [&](const Eigen::Matrix<double, 12, 12>& block,
                       const std::array<int, 4>& atom, int n) {
        for (int a = 0; a < n; ++a)
            for (int al = 0; al < 3; ++al)
                for (int b = 0; b < n; ++b)
                    for (int be = 0; be < 3; ++be)
                        m_h_nonbonded(3 * atom[a] + al, 3 * atom[b] + be)
                            += block(3 * a + al, 3 * b + be);
    };

    // --- non-covalent pairs (dispersion + electrostatics + repulsion) ---------
    for (const auto& nci : m_parameters.value("qmdff_ncis", json::array())) {
        const int i = nci["i"];
        const int j = nci["j"];
        const double c6 = nci.value("c6", 0.0);
        const double r0_bj = nci.value("r0_bj", 0.0);
        const double sr42 = nci.value("sr42", 0.0);
        const double zab = nci.value("zab", 0.0);
        const double alpha = nci.value("alpha", 0.0);
        const double qq = nci.value("qq", 0.0);
        const int nk = nci.value("nk", 5);

        std::vector<Eigen::Vector3d> pos{ position(i), position(j) };
        {   // energy and gradient of this pair, accumulated into the fixed remainder
            Eigen::Vector3d g0 = Eigen::Vector3d::Zero(), g1 = g0;
            double ed = 0.0, ee = 0.0, er = 0.0;
            QMDFFTerms::nonCovalentPair(pos[0], pos[1], c6, r0_bj, sr42, zab, alpha, qq, nk,
                ed, ee, er, g0, g1, true);
            e_nb += ed + ee + er;
            g_nb.row(i) += g0.transpose();
            g_nb.row(j) += g1.transpose();
        }
        Eigen::Matrix<double, 12, 12> block;
        finiteDifferenceBlock(2, pos, h,
            [&](const std::vector<Eigen::Vector3d>& p, std::vector<Eigen::Vector3d>& g) {
                g[0].setZero();
                g[1].setZero();
                double ed = 0.0, ee = 0.0, er = 0.0;
                QMDFFTerms::nonCovalentPair(p[0], p[1], c6, r0_bj, sr42, zab, alpha, qq, nk,
                    ed, ee, er, g[0], g[1], true);
            },
            block);
        scatter(block, { { i, j, -1, -1 } }, 2);
    }

    // --- hydrogen and halogen bonds ------------------------------------------
    for (const auto& hb : m_parameters.value("qmdff_hbonds", json::array())) {
        const int a = hb["a"];
        const int b = hb["b"];
        const int hh = hb["h"];
        const double c1 = hb.value("c1", 0.0);
        const double c2 = hb.value("c2", 0.0);
        const bool halogen = hb.value("halogen", false);

        // hydrogenBondEnergy bails out at E > -1e-8 (qmdff_terms.h, faithful to eabhag),
        // which makes the term discontinuous there. Differentiating across that step
        // produces garbage, so a negligible HB is simply left out of the target.
        if (!halogen) {
            Eigen::Vector3d g0 = Eigen::Vector3d::Zero(), g1 = g0, g2 = g0;
            const double e = QMDFFTerms::hydrogenBondEnergy(position(a), position(b), position(hh),
                c1, c2, g0, g1, g2, false);
            if (std::abs(e) < 1.0e-9)
                continue;
        }

        std::vector<Eigen::Vector3d> pos{ position(a), position(b), position(hh) };
        {
            Eigen::Vector3d g0 = Eigen::Vector3d::Zero(), g1 = g0, g2 = g0;
            e_nb += halogen
                ? QMDFFTerms::halogenBondEnergyGradient(pos[0], pos[1], pos[2], c1, g0, g1, g2, true)
                : QMDFFTerms::hydrogenBondEnergy(pos[0], pos[1], pos[2], c1, c2, g0, g1, g2, true);
            g_nb.row(a) += g0.transpose();
            g_nb.row(b) += g1.transpose();
            g_nb.row(hh) += g2.transpose();
        }
        Eigen::Matrix<double, 12, 12> block;
        finiteDifferenceBlock(3, pos, h,
            [&](const std::vector<Eigen::Vector3d>& p, std::vector<Eigen::Vector3d>& g) {
                g[0].setZero();
                g[1].setZero();
                g[2].setZero();
                if (halogen)
                    QMDFFTerms::halogenBondEnergyGradient(p[0], p[1], p[2], c1, g[0], g[1], g[2], true);
                else
                    QMDFFTerms::hydrogenBondEnergy(p[0], p[1], p[2], c1, c2, g[0], g[1], g[2], true);
            },
            block);
        scatter(block, { { a, b, hh, -1 } }, 3);
    }
}

// =============================================================================
// Stage 1: per-term starting guess
// =============================================================================


Vector QMDFFParametrisation::currentConstantsInternal() const
{
    Vector k = Vector::Zero(m_nparams);
    const json bonds = m_parameters.value("bonds", json::array());
    const json angles = m_parameters.value("angles", json::array());
    const json torsions = m_parameters.value("qmdff_torsions", json::array());
    const json dihedrals = m_parameters.value("dihedrals", json::array());
    const json extra = m_parameters.value("extra_dihedrals", json::array());
    const json inversions = m_parameters.value("inversions", json::array());
    const bool gfnff = (m_options.source == TermSource::GFN_FF);

    for (const auto& t : m_terms) {
        double value = 0.0;
        switch (t.kind) {
        case TermKind::Bond: value = bonds[t.list_index].value("fc", 0.0); break;
        case TermKind::Angle: value = angles[t.list_index].value("fc", 0.0); break;
        case TermKind::Torsion:
            value = gfnff ? dihedrals[t.list_index].value("V", 0.0)
                          : torsions[t.list_index].value("scale", 1.0);
            break;
        case TermKind::ExtraTorsion: value = extra[t.list_index].value("V", 0.0); break;
        case TermKind::OutOfPlane:
            value = gfnff ? inversions[t.list_index].value("fc", 0.0)
                          : torsions[t.list_index].value("scale", 1.0);
            break;
        }
        k(t.group) = parameterSign(t.kind, m_options.source) * value;
    }
    return k;
}

Vector QMDFFParametrisation::currentConstants() const
{
    Vector k = currentConstantsInternal();
    for (const auto& t : m_terms)
        k(t.group) *= parameterSign(t.kind, m_options.source);
    return k;
}

Matrix QMDFFParametrisation::currentBondedHessian() const
{
    const int n3 = 3 * m_natoms;
    Matrix h = Matrix::Zero(n3, n3);
    const Vector k = currentConstantsInternal();
    for (const auto& t : m_terms)
        t.addTo(h, k(t.group)); // internal convention, matching t.block
    return h;
}

Vector QMDFFParametrisation::seminarioGuess(const Matrix& h_target) const
{
    Vector k0 = Vector::Zero(m_nparams);
    Vector count = Vector::Zero(m_nparams);

    auto block3 = [&](int a, int b) {
        return Eigen::Matrix3d(h_target.block<3, 3>(3 * a, 3 * b));
    };
    auto position = [&](int a) {
        return Eigen::Vector3d(m_geometry_bohr(a, 0), m_geometry_bohr(a, 1), m_geometry_bohr(a, 2));
    };

    const json bonds = m_parameters.value("bonds", json::array());
    const json angles = m_parameters.value("angles", json::array());

    for (const auto& t : m_terms) {
        double guess = 0.0;
        bool have_seminario = false;

        if (t.kind == TermKind::Bond) {
            // k_AB from the interaction block, then map the harmonic curvature onto the
            // QMDFF shape: E''(r0) = fc * a^2 / (2 r0^2) for BOTH the LJ and the Morse
            // form, so the LJ->Morse promotion does not change the fitted constant.
            const double r0 = bonds[t.list_index].value("r0_ij", 0.0) * kAngstromToBohr;
            const double expo = bonds[t.list_index].value("exponent", 0.0);
            if (std::abs(expo) > 1e-8) {
                const Eigen::Vector3d u = (position(t.atom[1]) - position(t.atom[0])).normalized();
                const double k_ab = seminarioProjection(block3(t.atom[0], t.atom[1]), u);
                if (m_options.source == TermSource::GFN_FF) {
                    // E = fc exp(-alpha dr^2) -> E''(r0) = -2 alpha fc, and the solver works
                    // with s = -fc, so s = k_AB / (2 alpha).
                    guess = k_ab / (2.0 * expo);
                } else if (r0 > 1e-8) {
                    // QMDFF LJ/Morse: E''(r0) = fc a^2 / (2 r0^2)
                    guess = 2.0 * r0 * r0 * k_ab / (expo * expo);
                }
                have_seminario = std::isfinite(guess);
            }
        } else if (t.kind == TermKind::Angle) {
            const Eigen::Vector3d ra = position(t.atom[0]);
            const Eigen::Vector3d rb = position(t.atom[1]); // centre
            const Eigen::Vector3d rc = position(t.atom[2]);
            const double r_ab = (ra - rb).norm();
            const double r_cb = (rc - rb).norm();
            if (r_ab > 1e-8 && r_cb > 1e-8) {
                const Eigen::Vector3d u_ab = (ra - rb) / r_ab;
                const Eigen::Vector3d u_cb = (rc - rb) / r_cb;
                const Eigen::Vector3d nvec = u_cb.cross(u_ab);
                if (nvec.norm() > 1e-6) { // undefined for a linear arrangement
                    const Eigen::Vector3d nhat = nvec.normalized();
                    const Eigen::Vector3d u_pa = nhat.cross(u_ab);
                    const Eigen::Vector3d u_pc = u_cb.cross(nhat);

                    const double s_ab = r_ab * r_ab * seminarioProjection(block3(t.atom[0], t.atom[1]), u_pa);
                    const double s_cb = r_cb * r_cb * seminarioProjection(block3(t.atom[2], t.atom[1]), u_pc);
                    if (s_ab > 1e-12 && s_cb > 1e-12) {
                        const double k_theta = 1.0 / (1.0 / s_ab + 1.0 / s_cb);

                        // Map d2E/dtheta2 onto the QMDFF angle term. Use the SAME
                        // near-linear threshold the kernel uses, otherwise the guess and
                        // the evaluator disagree and k_theta/sin^2(theta0) explodes.
                        const double theta0 = angles[t.list_index].value("theta0_ijk", 0.0);
                        double damp_ab = 0.0, ddamp = 0.0, damp_cb = 0.0;
                        QMDFFTerms::abdamp(m_atoms[t.atom[0]], m_atoms[t.atom[1]], r_ab * r_ab, damp_ab, ddamp);
                        QMDFFTerms::abdamp(m_atoms[t.atom[2]], m_atoms[t.atom[1]], r_cb * r_cb, damp_cb, ddamp);
                        const double damp = damp_ab * damp_cb;

                        if (damp > 1e-12) {
                            const double s = std::sin(theta0);
                            const bool linear = (QMDFFTerms::kPi - theta0 < 1.0e-6);
                            guess = linear ? k_theta / (2.0 * damp)
                                           : k_theta / (2.0 * s * s * damp);
                            have_seminario = true;
                        }
                    }
                }
            }
        }

        if (have_seminario && std::isfinite(guess) && guess > 0.0) {
            k0(t.group) += guess;
            count(t.group) += 1.0;
        }
    }

    for (int p = 0; p < m_nparams; ++p)
        if (count(p) > 0.0)
            k0(p) /= count(p);

    return k0;
}

// =============================================================================
// Stage 2: global linear least squares
// =============================================================================

Matrix QMDFFParametrisation::projectorTR() const
{
    const int n3 = 3 * m_natoms;
    Matrix d = Matrix::Zero(n3, 6);

    Eigen::Vector3d com = Eigen::Vector3d::Zero();
    double mass_sum = 0.0;
    for (int a = 0; a < m_natoms; ++a) {
        const double m = Elements::AtomicMass[m_atoms[a]];
        com += m * m_geometry_bohr.row(a).transpose();
        mass_sum += m;
    }
    if (mass_sum > 0.0)
        com /= mass_sum;

    for (int a = 0; a < m_natoms; ++a) {
        for (int c = 0; c < 3; ++c)
            d(3 * a + c, c) = 1.0; // translations

        const Eigen::Vector3d r = m_geometry_bohr.row(a).transpose() - com;
        for (int c = 0; c < 3; ++c) {
            Eigen::Vector3d e = Eigen::Vector3d::Zero();
            e(c) = 1.0;
            const Eigen::Vector3d rot = r.cross(e);
            for (int k = 0; k < 3; ++k)
                d(3 * a + k, 3 + c) = rot(k);
        }
    }

    // Orthonormal basis of the T/R space; a linear molecule contributes only 5 vectors.
    Eigen::HouseholderQR<Matrix> qr(d);
    Matrix q = qr.householderQ() * Matrix::Identity(n3, 6);
    Matrix p = Matrix::Identity(n3, n3);
    for (int c = 0; c < 6; ++c) {
        const double nrm = d.col(c).norm();
        if (nrm < 1e-8)
            continue;
        const Vector v = q.col(c);
        p -= v * v.transpose();
    }
    return p;
}

Vector QMDFFParametrisation::boundedNNLS(const Matrix& A, const Vector& b,
    const Vector& floors, int max_iter)
{
    // Lawson-Hanson active set on the normal equations, shifted so the bounds are z >= 0.
    // Lower bounds rather than plain zero, so that a term is pinned at a small physical
    // floor instead of being silently deleted from the force field.
    const int n = static_cast<int>(A.rows());
    const Vector b_shift = b - A * floors;
    const double eps = 1e-12 * std::max(1.0, b_shift.cwiseAbs().maxCoeff());

    Vector z = Vector::Zero(n);
    std::vector<char> passive(n, 0);
    Vector w = b_shift;

    for (int outer = 0; outer < max_iter; ++outer) {
        int t = -1;
        double best = eps;
        for (int i = 0; i < n; ++i)
            if (!passive[i] && w(i) > best) {
                best = w(i);
                t = i;
            }
        if (t < 0)
            break; // KKT satisfied
        passive[t] = 1;

        for (int inner = 0; inner < 3 * n; ++inner) {
            std::vector<int> p_set;
            for (int i = 0; i < n; ++i)
                if (passive[i])
                    p_set.push_back(i);
            const int m = static_cast<int>(p_set.size());
            if (m == 0)
                break;

            Matrix ap(m, m);
            Vector bp(m);
            for (int a = 0; a < m; ++a) {
                bp(a) = b_shift(p_set[a]);
                for (int c = 0; c < m; ++c)
                    ap(a, c) = A(p_set[a], p_set[c]);
            }
            const Vector s = ap.ldlt().solve(bp);

            if (s.minCoeff() > 0.0) {
                z.setZero();
                for (int a = 0; a < m; ++a)
                    z(p_set[a]) = s(a);
                break;
            }

            double alpha = 1.0;
            for (int a = 0; a < m; ++a)
                if (s(a) <= 0.0) {
                    const double denom = z(p_set[a]) - s(a);
                    if (denom > 1e-300)
                        alpha = std::min(alpha, z(p_set[a]) / denom);
                }
            for (int a = 0; a < m; ++a)
                z(p_set[a]) += alpha * (s(a) - z(p_set[a]));
            for (int a = 0; a < m; ++a)
                if (z(p_set[a]) <= 1e-14) {
                    z(p_set[a]) = 0.0;
                    passive[p_set[a]] = 0;
                }
        }
        w = b_shift - A * z;
    }
    return floors + z;
}

Vector QMDFFParametrisation::solveGlobal(const Matrix& A, const Vector& b, const Vector& k0,
    const Vector& lambda, const Vector& floors) const
{
    const int n = static_cast<int>(A.rows());
    Vector diag = A.diagonal();
    for (int i = 0; i < n; ++i)
        if (diag(i) <= 0.0)
            diag(i) = 1.0; // an inactive parameter must still give a solvable row

    Matrix reg = A;
    Vector rhs = b;
    for (int i = 0; i < n; ++i) {
        reg(i, i) += lambda(i) * diag(i);
        rhs(i) += lambda(i) * diag(i) * k0(i);
    }

    if (m_options.nonnegative)
        return boundedNNLS(reg, rhs, floors, 4 * n + 16);

    // Rank-revealing rather than Cholesky. The diagonal of A spans many orders of
    // magnitude — a near-planar torsion has |U|^2 ~ 1e-5 while an angle has ~1e0 — and at
    // small lambda an LDLT on that loses several digits, which showed up as a 1e-3
    // residual on an exactly consistent synthetic system. The complete orthogonal
    // decomposition additionally returns the minimum-norm solution when A is rank
    // deficient (redundant torsions), instead of amplifying noise.
    return reg.completeOrthogonalDecomposition().solve(rhs);
}

// =============================================================================
// Driver
// =============================================================================

// =============================================================================
// Multi-basin fit: Hessian + relative-energy + gradient residuals
// =============================================================================
//
// Every residual type is linear in the SAME constants, so they stack into one system:
//   Hessian  at basin b :  sum_p k_p U_p(x_b)          ~ H_QM(x_b)  - H_nb(x_b)
//   energy   at basin b :  sum_p k_p [e_p(b)-e_p(0)]   ~ dE_QM(b,0) - dE_nb(b,0)
//   gradient at basin b :  sum_p k_p g_p(x_b)          ~ g_QM(x_b)  - g_nb(x_b)
//
// The energy block is the one that constrains conformer RANKING: a Hessian only sees
// curvature, which says nothing about the relative depth of two basins.
// Each block is normalised by its own target norm so the weights are dimensionless.

json QMDFFParametrisation::fitMultiBasin(const std::vector<BasinData>& basins, FitReport* report)
{
    FitReport local;
    FitReport& rep = report ? *report : local;
    rep.n_basins = static_cast<int>(basins.size());

    if (basins.empty()) {
        rep.warnings.push_back("no basins supplied — parameters returned unchanged");
        return m_parameters;
    }

    const int n3 = 3 * m_natoms;
    const int np = m_nparams;

    Matrix A_h = Matrix::Zero(np, np), A_e = Matrix::Zero(np, np), A_g = Matrix::Zero(np, np);
    Vector b_h = Vector::Zero(np), b_e = Vector::Zero(np), b_g = Vector::Zero(np);
    double t_h = 0.0, t_e = 0.0, t_g = 0.0;

    Matrix h_target_ref; // basin 0's Hessian target, used for the Seminario guess
    Vector e_ref = Vector::Zero(np);
    double e_qm_ref = 0.0;
    bool have_ref_energy = false;

    for (std::size_t ib = 0; ib < basins.size(); ++ib) {
        const BasinData& basin = basins[ib];
        const double w = basin.weight;
        const Matrix geom_bohr = basin.geometry_angstrom * kAngstromToBohr;

        std::vector<TermHessian> terms;
        buildTermData(geom_bohr, terms);
        Matrix h_nb;
        Matrix g_nb;
        double e_nb = 0.0;
        buildNonbonded(geom_bohr, h_nb, e_nb, g_nb);

        // Groups are assigned on the reference topology and are geometry independent.
        for (std::size_t i = 0; i < terms.size() && i < m_terms.size(); ++i)
            terms[i].group = m_terms[i].group;

        // Per-group energy and gradient at unit constant
        Vector e_unit = Vector::Zero(np);
        Matrix g_unit = Matrix::Zero(np, n3);
        for (const auto& t : terms) {
            e_unit(t.group) += t.energy;
            Vector row = g_unit.row(t.group).transpose();
            t.addGradTo(row, 1.0);
            g_unit.row(t.group) = row.transpose();
        }

        // ---- Hessian block ---------------------------------------------------
        if (basin.hessian_angstrom.rows() == n3) {
            Matrix h = hessianAngstromToBohr(basin.hessian_angstrom);
            h = 0.5 * (h + h.transpose()).eval();
            const Matrix h_target = h - h_nb;
            if (ib == 0)
                h_target_ref = h_target;

            const Vector no_weight;
            for (const auto& tp : terms)
                b_h(tp.group) += w * tp.dotHessian(h_target, no_weight);
            for (std::size_t ip = 0; ip < terms.size(); ++ip)
                for (std::size_t iq = 0; iq < terms.size(); ++iq)
                    A_h(terms[ip].group, terms[iq].group)
                        += w * terms[ip].dotTerm(terms[iq], no_weight);
            t_h += w * h_target.squaredNorm();
        }

        // ---- energy block (differences against basin 0) -----------------------
        if (basin.has_energy) {
            if (ib == 0) {
                e_ref = e_unit;
                e_qm_ref = basin.energy - e_nb;
                have_ref_energy = true;
            } else if (have_ref_energy) {
                const Vector d = e_unit - e_ref;
                const double target = (basin.energy - e_nb) - e_qm_ref;
                A_e += w * (d * d.transpose());
                b_e += w * target * d;
                t_e += w * target * target;
            }
        }

        // ---- gradient block ---------------------------------------------------
        if (basin.gradient_angstrom.rows() == m_natoms) {
            // dE/dx_Bohr = dE/dx_Angstrom / (Bohr per Angstrom)
            Vector gflat = Vector::Zero(n3);
            for (int a = 0; a < m_natoms; ++a)
                for (int c = 0; c < 3; ++c)
                    gflat(3 * a + c) = basin.gradient_angstrom(a, c) / kAngstromToBohr
                        - g_nb(a, c);
            A_g += w * (g_unit * g_unit.transpose());
            b_g += w * (g_unit * gflat);
            t_g += w * gflat.squaredNorm();
        }
    }

    // ---- combine, each block normalised by its own target norm ----------------
    auto scale = [](double weight, double target) {
        return (target > 1e-30) ? weight / target : 0.0;
    };
    const double s_h = scale(m_options.weight_hessian, t_h);
    const double s_e = scale(m_options.weight_energy, t_e);
    const double s_g = scale(m_options.weight_gradient, t_g);

    const Matrix A = s_h * A_h + s_e * A_e + s_g * A_g;
    const Vector b = s_h * b_h + s_e * b_e + s_g * b_g;
    const double target_norm2 = s_h * t_h + s_e * t_e + s_g * t_g;

    rep.target_norm = std::sqrt(std::max(0.0, target_norm2));
    if (t_e <= 1e-30)
        rep.warnings.push_back("no usable energy residuals — the fit constrains curvature "
                               "only, which cannot determine relative basin energies");

    const json out = solveAndWriteBack(A, b, target_norm2, h_target_ref, rep);

    // Per-block relative residuals, for diagnosis
    auto blockResidual = [&](const Matrix& Ax, const Vector& bx, double tx) {
        if (tx <= 1e-30) return 0.0;
        const Vector& k = rep.k_fitted;
        if (k.size() != np) return 0.0;
        return std::sqrt(std::max(0.0, tx - 2.0 * k.dot(bx) + k.dot(Ax * k)) / tx);
    };
    rep.rel_residual_energy = blockResidual(A_e, b_e, t_e);
    rep.rel_residual_gradient = blockResidual(A_g, b_g, t_g);

    if (m_options.verbosity >= 1) {
        CurcumaLogger::param("qmdff_fit_basins", fmt::format("{}", rep.n_basins));
        CurcumaLogger::param("qmdff_fit_block_residuals",
            fmt::format("energy {:.4f}, gradient {:.4f}",
                rep.rel_residual_energy, rep.rel_residual_gradient));
    }
    return out;
}

json QMDFFParametrisation::fit(const Matrix& hessian_qm_angstrom, FitReport* report)
{
    FitReport local;
    FitReport& rep = report ? *report : local;

    const int n3 = 3 * m_natoms;
    if (hessian_qm_angstrom.rows() != n3 || hessian_qm_angstrom.cols() != n3) {
        rep.warnings.push_back(fmt::format("Hessian has wrong dimension ({}x{}, expected {}x{}) "
                                           "— parameters returned unchanged",
            hessian_qm_angstrom.rows(), hessian_qm_angstrom.cols(), n3, n3));
        return m_parameters;
    }

    // --- target: QM Hessian minus everything we do not fit --------------------
    Matrix h = hessianAngstromToBohr(hessian_qm_angstrom);
    h = 0.5 * (h + h.transpose()).eval();

    // A plausibility guard on the magnitude. A missed Angstrom^2 -> Bohr^2 conversion is a
    // factor 3.5711 and a mass-weighted matrix lands inside the same band, so this only
    // catches gross errors — the synthetic round-trip ctest is the real check.
    const double hmax = h.diagonal().cwiseAbs().maxCoeff();
    if (hmax < 1e-4 || hmax > 10.0) {
        rep.warnings.push_back(fmt::format(
            "max |H_ii| = {:.4e} Eh/Bohr^2 is outside the plausible range [1e-4, 10] — "
            "wrong units, or a mass-weighted Hessian was passed instead of the raw one",
            hmax));
    }

    Matrix h_target = h - m_h_nonbonded;
    {
        // Diagnostic only — see FitOptions::project_tr for why this is NOT subtracted by
        // default. A nonzero value means the geometry is not a stationary point.
        const Matrix p = projectorTR();
        const Matrix projected = p * h_target * p;
        rep.tr_content = (h_target - projected).norm() / std::max(1e-30, h_target.norm());
        if (m_options.project_tr)
            h_target = projected;
    }
    if (rep.tr_content > 0.02) {
        rep.warnings.push_back(fmt::format(
            "{:.1f}% of the target Hessian lies in the translation/rotation subspace — the "
            "geometry is not a stationary point of the reference method, so part of the "
            "target cannot be represented by any internal-coordinate force field",
            100.0 * rep.tr_content));
    }

    // --- normal equations, assembled over shared-atom pairs only --------------
    // A_pq is non-zero only if terms p and q share an atom, so this is O(N) despite
    // looking like an all-pairs loop.
    std::vector<std::vector<int>> terms_of_atom(m_natoms);
    for (int p = 0; p < static_cast<int>(m_terms.size()); ++p)
        for (int a = 0; a < m_terms[p].natoms; ++a)
            terms_of_atom[m_terms[p].atom[a]].push_back(p);

    Matrix A = Matrix::Zero(m_nparams, m_nparams);
    Vector b = Vector::Zero(m_nparams);

    const Vector mw = m_options.mass_weighted ? Vector(m_weight.cwiseProduct(m_weight)) : Vector();
    auto overlap = [&](const TermHessian& tp, const TermHessian& tq) {
        return tp.dotTerm(tq, mw);
    };

    // Terms carrying the GFN-FF coordination-number correction are NOT local: their unit
    // Hessian reaches every atom the CN reaches, so they overlap each other even without a
    // shared atom. Collect them once and add them to every candidate set.
    std::vector<int> extra_terms;
    for (int q = 0; q < static_cast<int>(m_terms.size()); ++q)
        if (m_terms[q].hasExtra() && m_terms[q].extra_coeff != 0.0)
            extra_terms.push_back(q);

    for (int p = 0; p < static_cast<int>(m_terms.size()); ++p) {
        const TermHessian& tp = m_terms[p];

        // b_p = <U_p, H_target>
        b(tp.group) += tp.dotHessian(h_target, mw);

        // Candidate partners: every term sharing at least one atom
        std::set<int> candidates;
        for (int a = 0; a < tp.natoms; ++a)
            for (int q : terms_of_atom[tp.atom[a]])
                candidates.insert(q);
        candidates.insert(extra_terms.begin(), extra_terms.end());
        if (tp.hasExtra() && tp.extra_coeff != 0.0) {
            // p's rank-2 part has support wherever the coordination numbers reach, so it
            // overlaps terms it shares no atom with. Every term is a candidate.
            for (int q = 0; q < static_cast<int>(m_terms.size()); ++q)
                candidates.insert(q);
        }

        for (int q : candidates) {
            if (q < p)
                continue; // count each unordered pair once
            const double v = overlap(tp, m_terms[q]);
            A(tp.group, m_terms[q].group) += v;
            if (q != p)
                A(m_terms[q].group, tp.group) += v;
        }
    }

    // Weighted Frobenius norm of the target
    double target_norm2 = 0.0;
    for (int i = 0; i < n3; ++i)
        for (int j = 0; j < n3; ++j) {
            const double v = h_target(i, j) * m_weight(i) * m_weight(j);
            target_norm2 += v * v;
        }
    rep.target_norm = std::sqrt(target_norm2);
    return solveAndWriteBack(A, b, target_norm2, h_target, rep);
}

json QMDFFParametrisation::solveAndWriteBack(const Matrix& A, const Vector& b,
    double target_norm2, const Matrix& h_target, FitReport& rep) const
{

    // --- freeze numerically dead parameters -----------------------------------
    // A near-planar torsion has |sin(phi)| -> 0 and therefore a vanishing unit Hessian;
    // leaving it in the system lets the solver assign it arbitrary noise.
    const Vector k_current = currentConstantsInternal();

    const double max_diag = A.diagonal().maxCoeff();
    std::vector<char> frozen(m_nparams, 0);
    for (int p = 0; p < m_nparams; ++p)
        if (A(p, p) < 1e-14 * std::max(1.0, max_diag)) {
            frozen[p] = 1;
            ++rep.n_frozen;
        }

    // --- stage 1 --------------------------------------------------------------
    Vector k0 = seminarioGuess(h_target);
    for (int p = 0; p < m_nparams; ++p) {
        if (m_is_torsion_group[p]) {
            // Prior target for a torsion/out-of-plane scale: the rule-derived barrier it
            // came in with, not a Hessian projection (which carries no barrier information).
            k0(p) = k_current(p);
            continue;
        }
        if (frozen[p] || !(k0(p) > 0.0) || !std::isfinite(k0(p))) {
            // Rayleigh quotient — the exact single-parameter least-squares solution, well
            // defined for every term kind (Seminario has no torsion/out-of-plane analogue).
            k0(p) = (A(p, p) > 0.0) ? b(p) / A(p, p) : k_current(p);
        }
        if (!std::isfinite(k0(p)) || k0(p) < m_floors(p))
            k0(p) = std::max(m_floors(p), k_current(p));
    }

    auto residual = [&](const Vector& k) {
        const double r2 = std::max(0.0, target_norm2 - 2.0 * k.dot(b) + k.dot(A * k));
        return (target_norm2 > 0.0) ? std::sqrt(r2 / target_norm2) : 0.0;
    };
    rep.rel_residual_initial = residual(k0);

    // --- stage 2: pick lambda on the unconstrained problem, then solve ---------
    // Torsion and out-of-plane groups carry their own, much stronger prior: their constants
    // are barrier heights, which a Hessian at a single geometry does not determine. Without
    // it the solver drives them to zero and silently produces free internal rotation.
    auto lambdaVector = [&](double lam) {
        Vector v = Vector::Constant(m_nparams, lam);
        for (int i = 0; i < m_nparams; ++i)
            if (m_is_torsion_group[i])
                v(i) = m_options.lambda_torsion;
        return v;
    };

    double lambda = m_options.lambda;
    if (lambda < 0.0) {
        double best = std::numeric_limits<double>::max();
        for (double lam : m_options.lambda_ladder) {
            const Vector lam_vec = lambdaVector(lam);
            Matrix reg = A;
            Vector rhs = b;
            for (int i = 0; i < m_nparams; ++i) {
                const double d = (A(i, i) > 0.0) ? A(i, i) : 1.0;
                reg(i, i) += lam_vec(i) * d;
                rhs(i) += lam_vec(i) * d * k0(i);
            }
            const Vector k = reg.completeOrthogonalDecomposition().solve(rhs);
            const double r = residual(k);
            rep.lambda_curve.emplace_back(lam, r);
            best = std::min(best, r);
        }
        lambda = m_options.lambda_ladder.back();
        for (const auto& [lam, r] : rep.lambda_curve)
            if (r <= (1.0 + m_options.lambda_tol) * best) {
                lambda = lam;
                break;
            }
    }
    rep.lambda_used = lambda;

    Vector k = solveGlobal(A, b, k0, lambdaVector(lambda), m_floors);
    for (int p = 0; p < m_nparams; ++p) {
        if (frozen[p] || !std::isfinite(k(p)))
            k(p) = k_current(p);
        if (k(p) <= m_floors(p) + 1e-15)
            ++rep.n_at_floor;
    }

    rep.rel_residual_fitted = residual(k);
    rep.r_squared = 1.0 - rep.rel_residual_fitted * rep.rel_residual_fitted;
    rep.k_initial = k0;
    rep.k_fitted = k;
    rep.n_params = m_nparams;
    rep.n_terms = static_cast<int>(m_terms.size());
    rep.hessian_target = h_target;

    // --- write the constants back and rebuild the predicted Hessian ------------
    json out = m_parameters;
    Matrix h_pred = m_h_nonbonded;
    for (const auto& t : m_terms) {
        const double sign = parameterSign(t.kind, m_options.source);
        const double value = sign * k(t.group);   // back to the force field's convention
        switch (t.kind) {
        case TermKind::Bond:
            out["bonds"][t.list_index]["fc"] = value;
            ++rep.n_bonds;
            break;
        case TermKind::Angle:
            out["angles"][t.list_index]["fc"] = value;
            ++rep.n_angles;
            break;
        case TermKind::Torsion:
            if (m_options.source == TermSource::GFN_FF)
                out["dihedrals"][t.list_index]["V"] = value;
            else
                out["qmdff_torsions"][t.list_index]["scale"] = value;
            ++rep.n_torsions;
            break;
        case TermKind::ExtraTorsion:
            out["extra_dihedrals"][t.list_index]["V"] = value;
            ++rep.n_torsions;
            break;
        case TermKind::OutOfPlane:
            if (m_options.source == TermSource::GFN_FF)
                out["inversions"][t.list_index]["fc"] = value;
            else
                out["qmdff_torsions"][t.list_index]["scale"] = value;
            ++rep.n_oop;
            break;
        }
        // U_p already carries the sign flip, so the internal k enters here, not `value`.
        t.addTo(h_pred, k(t.group));
    }
    rep.hessian_predicted = h_pred;

    if (rep.n_at_floor > 0)
        rep.warnings.push_back(fmt::format("{} of {} force constants were pushed to their "
                                           "lower bound by the non-negativity constraint",
            rep.n_at_floor, m_nparams));
    if (rep.n_frozen > 0)
        rep.warnings.push_back(fmt::format("{} parameters have a numerically vanishing unit "
                                           "Hessian (near-planar torsions) and kept their "
                                           "initial value",
            rep.n_frozen));
    if (rep.r_squared < 0.9)
        rep.warnings.push_back(fmt::format("R^2 = {:.3f} — the QMDFF functional form cannot "
                                           "reproduce this Hessian well; the constants absorb "
                                           "physics QMDFF does not model and are not transferable",
            rep.r_squared));

    if (m_options.verbosity >= 1) {
        CurcumaLogger::param("qmdff_fit_parameters", fmt::format("{} ({} terms)", m_nparams, rep.n_terms));
        CurcumaLogger::param("qmdff_fit_lambda", fmt::format("{:.1e}", rep.lambda_used));
        CurcumaLogger::param("qmdff_fit_residual",
            fmt::format("{:.4f} -> {:.4f} (R^2 = {:.4f})",
                rep.rel_residual_initial, rep.rel_residual_fitted, rep.r_squared));
        for (const auto& w : rep.warnings)
            CurcumaLogger::warn(w);
    }

    return out;
}

} // namespace curcuma::qmdff
