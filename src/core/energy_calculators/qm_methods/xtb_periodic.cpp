/*
 * <Periodic (orthorhombic cell, Gamma point) kernels of the native xTB>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
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
 * Claude Generated (Nov 2026, WP-PERIODIC-XTB). Status: AI-generated, machine-tested, no human
 * production testing. Plan, decisions and tests: docs/WP-PERIODIC-XTB.md.
 *
 * Every function here is the periodic twin of a cell-free kernel in xtb_h0.cpp / xtb_native.cpp
 * and is reached only through `if (m_pbc)` at that kernel's entry, so the cell-free statements
 * stay verbatim (bitwise gate). Image sums use PBCUtils::forEachImage: lattice translations T
 * with |x_i - x_j - T| < cutoff, self images (i, i, T != 0) included, no L/2 restriction. The
 * cutoffs are those of tblite 0.4.0 (kPbc* in xtb_native.h), the reference of this work.
 *
 * P1 (this file): coordination number of the self-energies, repulsion, GFN1 halogen bond, GFN1 D3,
 * GFN2 D4 (CN, pairs, ATM, CN chain rule, dE/dq). The electronic part (S, H0, Coulomb, multipoles)
 * follows in P2.
 */

#include "xtb_native.h"

#include "parameters/gfn1_params.hpp"
#include "parameters/gfn2_params.hpp"
#include "parameters/xtb_params_extra.hpp"

#include "STO_CGTO.hpp"
#include "xtb_ao_utils.hpp"
#include "xtb_coulomb.hpp"

#include "src/core/curcuma_logger.h"
#include "src/core/energy_calculators/ff_methods/ewald.h"
#include "src/core/energy_calculators/dispersion/d4_evaluator.h"
#include "src/core/energy_calculators/dispersion/d4_ncoord.h"
#include "src/core/energy_calculators/dispersion/d4param_generator.h"
#include "src/core/energy_calculators/ff_methods/cn_calculator.h"
#include "src/core/energy_calculators/ff_methods/d3param_generator.h"

#include <cmath>
#include <limits>

namespace curcuma::xtb {

namespace {

/// Geometry (Angstrom rows) to a flat Bohr array, as the cell-free kernels do.
std::vector<double> geometryBohr(const Matrix& geometry, int nat)
{
    std::vector<double> xyz(3 * nat);
    for (int i = 0; i < nat; ++i) {
        xyz[3 * i + 0] = geometry(i, 0) * AA_TO_AU;
        xyz[3 * i + 1] = geometry(i, 1) * AA_TO_AU;
        xyz[3 * i + 2] = geometry(i, 2) * AA_TO_AU;
    }
    return xyz;
}

inline bool xb_is_halogen(int z) { return z == 17 || z == 35 || z == 53 || z == 85; }
inline bool xb_is_acceptor(int z) { return z == 7 || z == 8 || z == 15 || z == 16; }

} // namespace

/* ------------------------------------------------------------------ *
 *  Cell                                                               *
 * ------------------------------------------------------------------ */
bool XTB::setPeriodicCell(const Eigen::Matrix3d& cell_angstrom, std::string& error)
{
    error.clear();
    Eigen::Matrix3d cell_bohr = cell_angstrom * AA_TO_AU;
    PBCUtils::OrthoCell c = PBCUtils::makeOrthoCell(cell_bohr, error);
    if (!c.active) {
        if (error.empty()) error = "cell is not orthorhombic";
        error = "native " + methodName() + " with a unit cell: " + error
              + " (only orthorhombic cells are supported, docs/WP-PERIODIC-XTB.md)";
        return false;
    }
    if (m_method == MethodType::GFN2 && m_d4_charge_source != "mulliken") {
        error = "native gfn2 with a unit cell: d4_charge_source '" + m_d4_charge_source
              + "' is not implemented with periodic boundary conditions (the EEQ/CPSCF q-response "
                "is not periodic); use the default 'mulliken'";
        return false;
    }
    m_cell_bohr = c;
    m_pbc = true;
    m_cn_cache_valid = false;
    m_d4_prepared = false;
    return true;
}

/* ------------------------------------------------------------------ *
 *  Coordination number of the self-energies (tblite ncoord with trans) *
 *  GFN1: exp counting (kcn 16), GFN2: double exp (ka 10, kb 20, +2 Bohr) *
 * ------------------------------------------------------------------ */
Vector XTB::computeCoordinationNumbersPeriodic() const
{
    const int nat = m_atomcount;
    const std::vector<double> xyz = geometryBohr(m_geometry, nat);
    const PBCUtils::LatticeImageList images = PBCUtils::makeLatticeImageList(m_cell_bohr, kPbcCnCutoffBohr);
    std::vector<double> rcov(nat);
    for (int i = 0; i < nat; ++i) rcov[i] = covalent_rad_d3_au(m_atoms[i]);

    Vector cn = Vector::Zero(nat);
    const bool gfn1 = (m_method == MethodType::GFN1);
    for (int i = 0; i < nat; ++i) {
        for (int j = 0; j <= i; ++j) {
            const double rc = rcov[i] + rcov[j];
            const bool same = (i == j);
            PBCUtils::forEachImage(images, &xyz[3 * i], &xyz[3 * j], same,
                [&](double, double, double, double r2) {
                    const double r = std::sqrt(r2);
                    double c;
                    if (gfn1) {
                        c = 1.0 / (1.0 + std::exp(-16.0 * (rc / r - 1.0)));
                    } else {
                        constexpr double ka = 10.0, kb = 20.0, r_shift = 2.0;
                        c = (1.0 / (1.0 + std::exp(-ka * (rc / r - 1.0))))
                          * (1.0 / (1.0 + std::exp(-kb * ((rc + r_shift) / r - 1.0))));
                    }
                    cn(i) += c;
                    if (!same) cn(j) += c;
                });
        }
    }
    return cn;
}

/* ------------------------------------------------------------------ *
 *  Repulsion with images (tblite repulsion/effective.f90 with trans)  *
 * ------------------------------------------------------------------ */
double XTB::calcRepulsionEnergyPeriodic() const
{
    const int nat = m_atomcount;
    const std::vector<double> xyz = geometryBohr(m_geometry, nat);
    const PBCUtils::LatticeImageList images = PBCUtils::makeLatticeImageList(m_cell_bohr, kPbcRepulsionCutoffBohr);
    const bool g1 = (m_method == MethodType::GFN1);
    const double kexp = g1 ? gfn1_params::rep_kexp : gfn2_params::rep_kexp;
    const double rexp = g1 ? gfn1_params::rep_rexp : gfn2_params::rep_rexp;

    double erep = 0.0;
    for (int i = 0; i < nat; ++i) {
        const int zi = m_atoms[i];
        const double alfi = g1 ? gfn1_params::rep_alpha[zi - 1] : gfn2_params::rep_alpha[zi - 1];
        const double zeffi = g1 ? gfn1_params::rep_zeff[zi - 1] : gfn2_params::rep_zeff[zi - 1];
        for (int j = 0; j <= i; ++j) {
            const int zj = m_atoms[j];
            const double alfj = g1 ? gfn1_params::rep_alpha[zj - 1] : gfn2_params::rep_alpha[zj - 1];
            const double zeffj = g1 ? gfn1_params::rep_zeff[zj - 1] : gfn2_params::rep_zeff[zj - 1];
            const double alpha_pair = std::sqrt(alfi * alfj);
            double kexp_pair = kexp;
            if (!g1 && zi <= 2 && zj <= 2) kexp_pair = gfn2_params::rep_kexp_light;
            const bool same = (i == j);
            const double w = same ? 0.5 : 1.0;   // tblite: 0.5 dE to iat, (iat /= jat) 0.5 dE to jat
            PBCUtils::forEachImage(images, &xyz[3 * i], &xyz[3 * j], same,
                [&](double, double, double, double r2) {
                    const double r = std::sqrt(r2);
                    const double r_kexp = std::pow(r, kexp_pair);
                    erep += w * (zeffi * zeffj / std::pow(r, rexp)) * std::exp(-alpha_pair * r_kexp);
                });
        }
    }
    return erep;
}

/* ------------------------------------------------------------------ *
 *  GFN1 halogen bond with images (tblite classical/halogen.f90)       *
 *  X-A pairs over images within 20 Bohr; B = nearest atom or image of X *
 * ------------------------------------------------------------------ */
double XTB::calcHalogenBondEnergyPeriodic() const
{
    if (m_method != MethodType::GFN1) return 0.0;
    const int nat = m_atomcount;
    if (nat < 2) return 0.0;
    const std::vector<double> xyz = geometryBohr(m_geometry, nat);
    const PBCUtils::LatticeImageList images = PBCUtils::makeLatticeImageList(m_cell_bohr, kPbcHalogenCutoffBohr);

    double exb = 0.0;
    for (int x = 0; x < nat; ++x) {
        const int zx = m_atoms[x];
        if (!xb_is_halogen(zx)) continue;
        const double cc = gfn1_params::halogen_bond[zx - 1];
        if (cc == 0.0) continue;   // Cl in GFN1
        // nearest neighbour B of X over atoms and images (dxb = B - X)
        double best = std::numeric_limits<double>::max();
        double dxb[3] = { 0.0, 0.0, 0.0 };
        bool found = false;
        for (int m = 0; m < nat; ++m) {
            PBCUtils::forEachImage(images, &xyz[3 * x], &xyz[3 * m], m == x,
                [&](double dx, double dy, double dz, double r2) {
                    if (r2 < best) {
                        best = r2;
                        dxb[0] = -dx; dxb[1] = -dy; dxb[2] = -dz;   // d = X - (B + T)
                        found = true;
                    }
                });
        }
        if (!found) continue;
        const double d2bx = best;
        for (int a = 0; a < nat; ++a) {
            if (!xb_is_acceptor(m_atoms[a])) continue;
            const double r0ax = gfn1_params::halogen_radscale * (atomic_rad_au(zx) + atomic_rad_au(m_atoms[a]));
            PBCUtils::forEachImage(images, &xyz[3 * x], &xyz[3 * a], false,
                [&](double dx, double dy, double dz, double d2ax) {
                    const double dxa[3] = { -dx, -dy, -dz };   // A + T - X
                    const double dba[3] = { dxa[0] - dxb[0], dxa[1] - dxb[1], dxa[2] - dxb[2] };
                    const double d2ab = dba[0]*dba[0] + dba[1]*dba[1] + dba[2]*dba[2];
                    const double rax = std::sqrt(d2ax);
                    const double xy = std::sqrt(d2bx * d2ax);
                    const double term = (d2bx + d2ax - d2ab) / xy;
                    const double aterm = std::pow(0.5 - 0.25 * term, 6.0);
                    const double t13 = r0ax / rax;
                    const double t14 = std::pow(t13, 12.0);
                    exb += aterm * cc * (t14 - gfn1_params::halogen_damping * std::pow(t13, 6.0)) / (1.0 + t14);
                });
        }
    }
    return exb;
}

/* ------------------------------------------------------------------ *
 *  GFN1 D3 with images                                                *
 * ------------------------------------------------------------------ */
double XTB::calcD3DispersionPeriodic(bool need_gradient) const
{
    if (!m_d3_generator)
        m_d3_generator = std::make_unique<::D3ParameterGenerator>(::D3ParameterGenerator::createForGFN1());
    ::D3ParameterGenerator& d3 = *m_d3_generator;
    const int nat = m_atomcount;

    // D3 coordination number (exp counting, k1 16, k2 4/3 on the D3 covalent radii) over images within 25 Bohr
    const std::vector<double> xyz = geometryBohr(m_geometry, nat);
    const PBCUtils::LatticeImageList cn_images = PBCUtils::makeLatticeImageList(m_cell_bohr, kPbcD3CnCutoffBohr);
    const std::vector<double>& rcov_A = CNCalculator::covalentRadiiD3();
    std::vector<double> cn(nat, 0.0);
    for (int i = 0; i < nat; ++i) {
        const int ei = m_atoms[i] - 1;
        if (ei < 0 || ei >= static_cast<int>(rcov_A.size())) continue;
        for (int j = 0; j <= i; ++j) {
            const int ej = m_atoms[j] - 1;
            if (ej < 0 || ej >= static_cast<int>(rcov_A.size())) continue;
            const double rc_bohr = (4.0 / 3.0) * (rcov_A[ei] + rcov_A[ej]) * AA_TO_AU;
            const bool same = (i == j);
            PBCUtils::forEachImage(cn_images, &xyz[3 * i], &xyz[3 * j], same,
                [&](double, double, double, double r2) {
                    const double c = 1.0 / (1.0 + std::exp(-16.0 * (rc_bohr / std::sqrt(r2) - 1.0)));
                    cn[i] += c;
                    if (!same) cn[j] += c;
                });
        }
    }
    d3.prepareForEnergyGradientWithCN(m_atoms, m_geometry, cn);

    const PBCUtils::LatticeImageList pair_images = PBCUtils::makeLatticeImageList(m_cell_bohr, kPbcD3PairCutoffBohr);
    m_disp_gradient = Matrix::Zero(nat, 3);
    m_disp_dEdcn = need_gradient ? Vector::Zero(nat) : Vector();
    m_disp_dEdq = Vector();
    const double e_disp = d3.getEnergyAndGradientPeriodic(need_gradient, m_disp_gradient, m_disp_dEdcn, pair_images);
    if (!need_gradient) {
        m_disp_gradient_valid = false;
        return e_disp;
    }
    // CN chain rule with the same image sum (P3 gate): dE/dCN_k * dCN_k/dx over images (j, T) of each pair
    if (m_disp_dEdcn.size() == nat) {
        for (int i = 0; i < nat; ++i) {
            const int ei = m_atoms[i] - 1;
            if (ei < 0 || ei >= static_cast<int>(rcov_A.size())) continue;
            for (int j = 0; j < i; ++j) {   // self images: no Cartesian derivative
                const int ej = m_atoms[j] - 1;
                if (ej < 0 || ej >= static_cast<int>(rcov_A.size())) continue;
                const double rc_bohr = (4.0 / 3.0) * (rcov_A[ei] + rcov_A[ej]) * AA_TO_AU;
                PBCUtils::forEachImage(cn_images, &xyz[3 * i], &xyz[3 * j], false,
                    [&](double dx, double dy, double dz, double r2) {
                        const double r = std::sqrt(r2);
                        const double ef = std::exp(-16.0 * (rc_bohr / r - 1.0));
                        const double c = 1.0 / (1.0 + ef);
                        const double dcdr = -c * (1.0 - c) * 16.0 * rc_bohr / r2;
                        const double f = (m_disp_dEdcn(i) + m_disp_dEdcn(j)) * dcdr / r;
                        m_disp_gradient(i, 0) += f * dx; m_disp_gradient(j, 0) -= f * dx;
                        m_disp_gradient(i, 1) += f * dy; m_disp_gradient(j, 1) -= f * dy;
                        m_disp_gradient(i, 2) += f * dz; m_disp_gradient(j, 2) -= f * dz;
                    });
            }
        }
        m_disp_dEdcn = Vector();
    }
    m_disp_gradient_valid = true;
    return e_disp;
}

/* ------------------------------------------------------------------ *
 *  GFN2 D4 with images                                                *
 * ------------------------------------------------------------------ */
void XTB::primeD4CnPeriodic(const Matrix& geom_bohr) const
{
    if (!m_d4_generator) return;
    const PBCUtils::LatticeImageList images = PBCUtils::makeLatticeImageList(m_cell_bohr, kPbcD4CnCutoffBohr);
    m_d4_generator->setPrecomputedCovalentCN(
        curcuma::dispersion::computeD4CovalentCNPeriodic(m_atoms, geom_bohr, images));
}

double XTB::calcD4DispersionPeriodic(bool need_gradient) const
{
    (void)need_gradient;   // the gradient block is cheap and always filled, as on the cell-free path
    if (!m_d4_generator) {
        nlohmann::json d4_config;
        d4_config["d4_s6"] = 1.0;
        d4_config["d4_s8"] = 2.7;
        d4_config["d4_a1"] = 0.52;
        d4_config["d4_a2"] = 5.0;
        d4_config["d4_alp"] = 16.0;
        ConfigManager cfg("d4param", d4_config);
        m_d4_generator = std::make_unique<::D4ParameterGenerator>(cfg);
        m_d4_generator->setBuildPairLists(false);
    }
    if (!m_d4_evaluator) {
        curcuma::dispersion::D4Params p;
        p.s6 = 1.0; p.s8 = 2.7; p.a1 = 0.52; p.a2 = 5.0; p.s9 = 5.0; p.alpha = 16.0;
        p.damping = curcuma::dispersion::DampingFormula::StandardBJ_D4;
        p.per_reference_charge = true;
        m_d4_evaluator = std::make_unique<curcuma::dispersion::D4Evaluator>(m_d4_generator.get(), p);
    }
    const int nat = m_atomcount;
    const Matrix geom_bohr = m_geometry * AA_TO_AU;
    m_d4_generator->setTopologyCharges(m_wfn.q_at);
    if (!m_d4_prepared) {
        m_d4_generator->setUseD4SingleShotEEQ(false);
        m_d4_generator->setD4CovalentCN(true);
        primeD4CnPeriodic(geom_bohr);
        m_d4_generator->GenerateParameters(m_atoms, geom_bohr);
        ++m_d4_genparams_calls;
        m_d4_prepared = true;
        m_d4_evaluator->invalidatePairCache();
    }

    const PBCUtils::LatticeImageList pair_images = PBCUtils::makeLatticeImageList(m_cell_bohr, kPbcD4PairCutoffBohr);
    const PBCUtils::LatticeImageList atm_images = PBCUtils::makeLatticeImageList(m_cell_bohr, m_d4_atm_cutoff);
    const PBCUtils::LatticeImageList cn_images = PBCUtils::makeLatticeImageList(m_cell_bohr, kPbcD4CnCutoffBohr);

    double E = m_d4_evaluator->computeEnergyAndGradientPeriodic(m_atoms, geom_bohr, pair_images,
        /*with_gradient=*/true, m_disp_gradient, m_disp_dEdcn, m_disp_dEdq, /*with_dEdq=*/true);
    E += m_d4_evaluator->computeATMPeriodic(m_atoms, geom_bohr, atm_images, /*with_gradient=*/true,
        m_disp_gradient, m_disp_dEdcn);

    // CN chain rule with the D4 CN's own derivative over images (as the cell-free path, d4_ncoord)
    if (m_disp_dEdcn.size() == nat) {
        curcuma::dispersion::addD4CovalentCNGradientPeriodic(m_atoms, geom_bohr, m_disp_dEdcn, m_disp_gradient, cn_images);
        m_disp_dEdcn = Vector();
    }
    // q-response: only the variational route (dE_D4/dq in the Fock potential, charge Pulay) exists with a
    // cell; setPeriodicCell() refuses the other d4_charge_source values.
    m_disp_gradient_valid = true;
    if (CurcumaLogger::get_verbosity() >= 2)
        CurcumaLogger::param("dispersion_d4_native_periodic", fmt::format("{:.6f} Eh", E));
    return E;
}

/* ------------------------------------------------------------------ *
 *  P2a: overlap and bare Hamiltonian at the Gamma point               *
 *  S(Gamma) = sum_T S(T), H0 likewise (tblite xtb/h0.f90 with trans):   *
 *  every shell pair over the images within the integral cutoff, the   *
 *  self images (iat == jat, T != 0) with the off-site hscale/shpoly    *
 *  of Z_i = Z_j, the on-atom T = 0 block exactly as the vacuum kernel. *
 * ------------------------------------------------------------------ */
double XTB::pbcIntegralCutoffBohr() const
{
    if (m_pbc_integral_cutoff > 0.0) return m_pbc_integral_cutoff;
    // tblite basis/type.f90 get_cutoff: intcut = 25 at accuracy 1; cutoff = min(sqrt(2 intcut / min_alpha), 40)
    double min_alpha = std::numeric_limits<double>::max();
    for (const auto& sh : m_basis.cgto)
        for (double a : sh.alpha)
            if (a > 0.0 && a < min_alpha) min_alpha = a;
    if (!(min_alpha < std::numeric_limits<double>::max())) return 40.0;
    return std::min(std::sqrt(2.0 * 25.0 / min_alpha), 40.0);
}

void XTB::getHamiltonianH0Periodic(const Vector& se, Matrix& S, Matrix& H0) const
{
    const int nat = m_atomcount;
    const int nsh = m_basis.nsh;
    const int nao = m_basis.nao;
    const std::vector<double> xyz_bohr = geometryBohr(m_geometry, nat);
    const PBCUtils::LatticeImageList images = PBCUtils::makeLatticeImageList(m_cell_bohr, pbcIntegralCutoffBohr());

    S = Matrix::Zero(nao, nao);
    H0 = Matrix::Zero(nao, nao);

    std::vector<bool> valence(nsh, false);
    if (m_method == MethodType::GFN1) {
        for (int iat = 0; iat < nat; ++iat) {
            bool ang_seen[3] = { false, false, false };
            for (int ish = 0; ish < m_basis.nsh_at[iat]; ++ish) {
                const int sh = m_basis.ish_at[iat] + ish;
                const int l = m_basis.ang_sh[sh];
                if (!ang_seen[l]) { valence[sh] = true; ang_seen[l] = true; }
            }
        }
    }
    std::vector<int> sh_z(nsh);
    for (int s = 0; s < nsh; ++s) sh_z[s] = m_basis.z[m_basis.sh2at[s]];
    std::vector<CGTO::Shell> shells(nsh);
    for (int ish = 0; ish < nsh; ++ish) shells[ish] = as_cgto_shell(m_basis.cgto[ish]);

    const int n_threads = effectiveIntraThreads(nsh);
    parallelStripes(n_threads, [&](int tid, int nth) {
    for (int ish_a = tid; ish_a < nsh; ish_a += nth) {
        const int iat = m_basis.sh2at[ish_a];
        const int ia_start = m_basis.iao_sh[ish_a];
        const int ia_nao = m_basis.nao_sh[ish_a];
        const CGTO::Shell& sh_a = shells[ish_a];
        const double zeta_a = m_basis.cgto[ish_a].slater_exp;
        const int ang_a = m_basis.ang_sh[ish_a];
        const double* xi = &xyz_bohr[3 * iat];

        for (int ish_b = 0; ish_b < nsh; ++ish_b) {
            const int jat = m_basis.sh2at[ish_b];
            const int jb_start = m_basis.iao_sh[ish_b];
            const int jb_nao = m_basis.nao_sh[ish_b];
            const CGTO::Shell& sh_b = shells[ish_b];
            const int ang_b = m_basis.ang_sh[ish_b];
            const double avg_eps = 0.5 * (se(ish_a) + se(ish_b));
            const int zi = sh_z[ish_a], zj = sh_z[ish_b];

            // off-site scaling (distance independent part), also for the self images
            double hs;
            if (m_method == MethodType::GFN1) {
                const bool vi = valence[ish_a], vj = valence[ish_b];
                if (vi && vj) {
                    const double den = std::pow(pauling_en[zi - 1] - pauling_en[zj - 1], 2);
                    hs = gfn1::kpair(zi, zj) * gfn1::kshell(ang_a, ang_b) * (1.0 + gfn1::enscale * den);
                } else if (vi && !vj) {
                    hs = 0.5 * (gfn1::kshell(ang_a, ang_a) + gfn1::kdiff);
                } else if (!vi && vj) {
                    hs = 0.5 * (gfn1::kshell(ang_b, ang_b) + gfn1::kdiff);
                } else {
                    hs = gfn1::kdiff;
                }
            } else {
                const double zeta_b = m_basis.cgto[ish_b].slater_exp;
                const double den = std::pow(pauling_en[zi - 1] - pauling_en[zj - 1], 2);
                const double enp = 1.0 + gfn2::enscale * den;
                const double km = gfn2::kpair(zi, zj) * gfn2::kshell(ang_a, ang_b) * enp;
                const double zij = std::pow(2.0 * std::sqrt(zeta_a * zeta_b) / (zeta_a + zeta_b), gfn2::wexp);
                hs = zij * km;
            }
            const double rad_sum = atomic_rad_au(zi) + atomic_rad_au(zj);

            // one image block: positions of a at xi, of b at xj_img; h_factor given
            auto addBlock = [&](const double* xj_img, double h_factor, bool onsite) {
                if (ang_a < 2 && ang_b < 2) {
                    double blk_s[9];
                    CGTO::cgto_overlap_block(sh_a, ang_a, sh_b, ang_b, xi[0], xi[1], xi[2],
                                             xj_img[0], xj_img[1], xj_img[2], blk_s);
                    for (int ia = 0; ia < ia_nao; ++ia) {
                        const int mu = ia_start + ia;
                        const int t_a = ao_to_type(ang_a, ia);
                        if (t_a < 0) continue;
                        for (int jb = 0; jb < jb_nao; ++jb) {
                            const int nu = jb_start + jb;
                            const int t_b = ao_to_type(ang_b, jb);
                            if (t_b < 0) continue;
                            const double s_ab = (onsite && ish_a == ish_b && t_a == t_b)
                                ? 1.0 : blk_s[ia * jb_nao + jb];
                            S(mu, nu) += s_ab;
                            H0(mu, nu) += avg_eps * h_factor * s_ab;
                        }
                    }
                } else {
                    double blk[6 * 6];
                    sphericalOverlapBlock(sh_a, ang_a, sh_b, ang_b, xi[0], xi[1], xi[2],
                                          xj_img[0], xj_img[1], xj_img[2], blk, 6);
                    for (int ia = 0; ia < ia_nao; ++ia) {
                        const int mu = ia_start + ia;
                        for (int jb = 0; jb < jb_nao; ++jb) {
                            const int nu = jb_start + jb;
                            double s_ab = blk[ia * 6 + jb];
                            if (onsite && ish_a == ish_b && ia == jb) s_ab = 1.0;
                            S(mu, nu) += s_ab;
                            H0(mu, nu) += avg_eps * h_factor * s_ab;
                        }
                    }
                }
            };

            if (iat == jat)
                addBlock(&xyz_bohr[3 * jat], 1.0, true);   // T = 0 on-atom block (vacuum kernel)

            // images: d = x_i - (x_j + T), so the image of j sits at x_i - d
            PBCUtils::forEachImage(images, xi, &xyz_bohr[3 * jat], iat == jat,
                [&](double dx, double dy, double dz, double r2) {
                    const double xj_img[3] = { xi[0] - dx, xi[1] - dy, xi[2] - dz };
                    const double rr = std::sqrt(std::sqrt(r2) / rad_sum);
                    const double pi_ij = (1.0 + m_h0.shpoly[ish_a] * rr) * (1.0 + m_h0.shpoly[ish_b] * rr);
                    addBlock(xj_img, hs * pi_ij, false);
                });
        }
    }
    });
}

/* ------------------------------------------------------------------ *
 *  P2a: isotropic Coulomb gamma with Ewald summation                  *
 *  gamma(is, js) = A_pt(at(is), at(js))                              *
 *                + sum_T fcut(r) [KO_gam(r) - 1/r]   (r < switch)     *
 *                + on-atom hardness as in the vacuum kernel           *
 *  A_pt: point-charge Ewald matrix (real erfc(a r)/r under the         *
 *  minimum image, reciprocal, self -2a/sqrt(pi), background           *
 *  -pi/(a^2 V)); tblite effective.f90 get_amat_3d has the same split   *
 *  except for the background term (docs/WP-PERIODIC-XTB.md, P0.4).    *
 * ------------------------------------------------------------------ */
namespace {
// tblite effective.f90 fsmooth: 1 below rcut - 1, 0 above rcut, quintic in between
inline double koSwitch(double r, double rcut)
{
    constexpr double offset = 1.0;
    if (r < rcut - offset) return 1.0;
    if (r > rcut) return 0.0;
    const double x = (r - (rcut - offset)) / offset;
    return -6.0 * x * x * x * x * x + 15.0 * x * x * x * x - 10.0 * x * x * x + 1.0;
}
} // namespace

void XTB::buildGammaMatrixPeriodic()
{
    const int nat = m_atomcount;
    const int nsh = m_basis.nsh;
    const std::vector<double> xyz_bohr = geometryBohr(m_geometry, nat);
    const coulomb::Method meth = (m_method == MethodType::GFN1) ? coulomb::Method::GFN1 : coulomb::Method::GFN2;
    const double gexp = (m_method == MethodType::GFN1) ? gfn1_params::gexp : gfn2_params::gexp;
    const double inv_gexp = 1.0 / gexp;

    // 1. point-charge Ewald matrix of the atoms
    Eigen::MatrixXd pos(nat, 3);
    for (int i = 0; i < nat; ++i)
        for (int k = 0; k < 3; ++k) pos(i, k) = xyz_bohr[3 * i + k];
    const curcuma::ewald::EwaldParameters par = curcuma::ewald::chooseEwaldParameters(
        m_cell_bohr, m_pbc_ewald_tolerance, m_pbc_ewald_r_cut > 0.0 ? m_pbc_ewald_r_cut : -1.0);
    const curcuma::ewald::EwaldSum ewald(m_cell_bohr, par);
    const Eigen::MatrixXd A = ewald.matrix(pos);
    m_pbc_ewald_alpha_used = par.alpha;
    m_pbc_ewald_nk_used = ewald.numberOfKVectors();

    // 2. shell matrix: Ewald part per atom pair plus the short-range Klopman-Ohno remainder over images
    std::vector<double> g(nsh);
    for (int is = 0; is < nsh; ++is)
        g[is] = coulomb::shell_hardness(meth, m_basis.z[m_basis.sh2at[is]], m_basis.ang_sh[is]);
    const PBCUtils::LatticeImageList images = PBCUtils::makeLatticeImageList(m_cell_bohr, m_pbc_coulomb_switch);
    const double rcut = m_pbc_coulomb_switch;

    m_gamma = Eigen::MatrixXd::Zero(nsh, nsh);
    for (int is = 0; is < nsh; ++is) {
        const int iat = m_basis.sh2at[is];
        for (int js = 0; js <= is; ++js) {
            const int jat = m_basis.sh2at[js];
            const double gam = (is == js) ? g[is] : coulomb::average(meth, g[is], g[js]);
            const double gam_pow = std::pow(gam, -gexp);
            double ko = 0.0;
            PBCUtils::forEachImage(images, &xyz_bohr[3 * iat], &xyz_bohr[3 * jat], iat == jat,
                [&](double, double, double, double r2) {
                    const double r1 = std::sqrt(r2);
                    const double f = koSwitch(r1, rcut);
                    if (f == 0.0) return;
                    ko += f * (std::pow(std::pow(r1, gexp) + gam_pow, -inv_gexp) - 1.0 / r1);
                });
            double val = A(iat, jat) + ko;
            if (iat == jat) val += gam;   // on-atom: raw hardness (is == js) or the pair average
            m_gamma(is, js) = val;
            m_gamma(js, is) = val;
        }
    }
}

/* ------------------------------------------------------------------ *
 *  Test hook (P1): dispersion at given charges                        *
 * ------------------------------------------------------------------ */
double XTB::dispersionEnergyNow(const Vector& q_at, bool gradient)
{
    if (q_at.size() == m_atomcount) m_wfn.q_at = q_at;
    else m_wfn.q_at = Vector::Zero(m_atomcount);
    m_d4_prepared = false;
    m_d4_variational_qresp = true;   // no separate q-response in the probe
    return calcDispersionEnergy(gradient);
}

} // namespace curcuma::xtb
