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
#include "xtb_multipole_ints.hpp"

#include "src/core/curcuma_logger.h"
#include "src/core/energy_calculators/ff_methods/ewald.h"
#include "src/core/energy_calculators/dispersion/d4_evaluator.h"
#include "src/core/energy_calculators/dispersion/d4_ncoord.h"
#include "src/core/energy_calculators/dispersion/d4param_generator.h"
#include "src/core/energy_calculators/ff_methods/cn_calculator.h"
#include "src/core/energy_calculators/ff_methods/d3param_generator.h"

#include <chrono>
#include <cstdlib>
#include <string>
#include <cmath>
#include <limits>

namespace MI = curcuma::xtb::multipole_ints;

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

namespace {
/// G5d (Nov 2026): per-atom repulsion inputs and the image list as flat arrays for GpuScfBackend::repulsionPeriodic.
struct PbcRepulsionArrays {
    std::vector<double> alpha, zeff, trans;
    std::vector<int> light;
};
PbcRepulsionArrays pbcRepulsionArrays(const std::vector<int>& atoms, bool g1, const PBCUtils::LatticeImageList& images)
{
    PbcRepulsionArrays a;
    const int nat = static_cast<int>(atoms.size());
    a.alpha.resize(nat); a.zeff.resize(nat); a.light.resize(nat);
    for (int i = 0; i < nat; ++i) {
        const int z = atoms[i];
        a.alpha[i] = g1 ? gfn1_params::rep_alpha[z - 1] : gfn2_params::rep_alpha[z - 1];
        a.zeff[i] = g1 ? gfn1_params::rep_zeff[z - 1] : gfn2_params::rep_zeff[z - 1];
        a.light[i] = (!g1 && z <= 2) ? 1 : 0;
    }
    a.trans.resize(3 * images.trans.size());
    for (size_t t = 0; t < images.trans.size(); ++t)
        for (int k = 0; k < 3; ++k) a.trans[3 * t + k] = images.trans[t][k];
    return a;
}
} // namespace

/// G2-G6 (Nov 2026): CURCUMA_PBC_GPU_HOST=<list> keeps the named periodic device stages on the host (comparison
/// and timing in one binary; tokens: atm, mpmat, mpgrad, d4dedq, gamma, wmat, cnchain, rep, all). Diagnostic only, default empty = every ported stage on the device.
bool pbcGpuStageOnHost(const char* token)   // declared in xtb_native.h
{
    const char* e = std::getenv("CURCUMA_PBC_GPU_HOST");
    if (!e) return false;
    const std::string v(e);
    return v.find("all") != std::string::npos || v.find(token) != std::string::npos;
}


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
    if (m_has_dshell) {
        error = "native " + methodName() + " with a unit cell: d-shell elements are not implemented with periodic "
                "boundary conditions (the periodic gradient of the d multipole integrals, P3)";
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

    // G5d (Nov 2026): the pair loop on the CUDA device when the periodic device SCF is on; the host loop below runs on
    // any failure or with CURCUMA_PBC_GPU_HOST=rep.
    if (m_gpu_scf && m_pbc_device_scf && !pbcGpuStageOnHost("rep")) {
        const PbcRepulsionArrays ra = pbcRepulsionArrays(m_atoms, g1, images);
        double e_dev = 0.0;
        const bool ok_dev = m_gpu_scf->repulsionPeriodic(nat, xyz.data(), ra.alpha.data(), ra.zeff.data(),
            ra.light.data(), kexp, g1 ? kexp : gfn2_params::rep_kexp_light, rexp, m_cell_bohr.L,
            static_cast<int>(images.trans.size()), ra.trans.data(), images.cutoff2, false, false, &e_dev, nullptr, nullptr);
        if (CurcumaLogger::get_verbosity() >= 3)
            CurcumaLogger::info_fmt("  periodic repulsion energy on the {}", ok_dev ? "device" : "host (device call failed)");
        if (ok_dev) return e_dev;
    }

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
    m_disp_virial.setZero();
    Eigen::Matrix3d* vir = (need_gradient && m_want_virial) ? &m_disp_virial : nullptr;
    const double e_disp = d3.getEnergyAndGradientPeriodic(need_gradient, m_disp_gradient, m_disp_dEdcn, pair_images, vir);
    if (!need_gradient) {
        m_disp_gradient_valid = false;
        return e_disp;
    }
    // CN chain rule with the same image sum (P3 gate): dE/dCN_k * dCN_k/dx over images (j, T) of each pair
    if (m_disp_dEdcn.size() == nat) {
        for (int i = 0; i < nat; ++i) {
            const int ei = m_atoms[i] - 1;
            if (ei < 0 || ei >= static_cast<int>(rcov_A.size())) continue;
            for (int j = 0; j <= (vir ? i : i - 1); ++j) {   // self images: no Cartesian derivative, virial only (P4)
                const int ej = m_atoms[j] - 1;
                if (ej < 0 || ej >= static_cast<int>(rcov_A.size())) continue;
                const bool same = (i == j);
                const double rc_bohr = (4.0 / 3.0) * (rcov_A[ei] + rcov_A[ej]) * AA_TO_AU;
                PBCUtils::forEachImage(cn_images, &xyz[3 * i], &xyz[3 * j], same,
                    [&](double dx, double dy, double dz, double r2) {
                        const double r = std::sqrt(r2);
                        const double ef = std::exp(-16.0 * (rc_bohr / r - 1.0));
                        const double c = 1.0 / (1.0 + ef);
                        const double dcdr = -c * (1.0 - c) * 16.0 * rc_bohr / r2;
                        const double f = (same ? m_disp_dEdcn(i) : m_disp_dEdcn(i) + m_disp_dEdcn(j)) * dcdr / r;
                        if (!same) {
                            m_disp_gradient(i, 0) += f * dx; m_disp_gradient(j, 0) -= f * dx;
                            m_disp_gradient(i, 1) += f * dy; m_disp_gradient(j, 1) -= f * dy;
                            m_disp_gradient(i, 2) += f * dz; m_disp_gradient(j, 2) -= f * dz;
                        }
                        if (vir) {
                            const double d[3] = { dx, dy, dz };
                            for (int a = 0; a < 3; ++a)
                                for (int b = 0; b < 3; ++b) (*vir)(a, b) -= f * d[a] * d[b];
                        }
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

    // Stage T: the periodic pair and ATM loops run on fixed blocks over the gated intra-thread budget
    m_d4_evaluator->setThreads(effectiveIntraThreads(nat));
    using d4clk = std::chrono::steady_clock;
    const auto td0 = d4clk::now();
    m_disp_virial.setZero();
    Eigen::Matrix3d* vir = m_want_virial ? &m_disp_virial : nullptr;   // P4
    double E = m_d4_evaluator->computeEnergyAndGradientPeriodic(m_atoms, geom_bohr, pair_images,
        /*with_gradient=*/true, m_disp_gradient, m_disp_dEdcn, m_disp_dEdq, /*with_dEdq=*/true, vir);
    const auto td1 = d4clk::now();
    // G4a (Nov 2026): the ATM on the CUDA device when the periodic device SCF is on (G1 wrapper); every failure
    // (backend without the kernel, allocation) falls back to the host loop below.
    bool atm_device = false;
    double t_atm_c6 = 0.0;
    if (m_gpu_scf && m_pbc_device_scf && !pbcGpuStageOnHost("atm")) {
        std::vector<double> c6f, dc6f;
        m_d4_generator->buildAtmC6Flat(m_atoms, c6f, dc6f);
        t_atm_c6 = std::chrono::duration<double, std::milli>(d4clk::now() - td1).count();
        std::vector<double> xyz(3 * nat), r4r2(nat), trans(3 * atm_images.trans.size());
        for (int i = 0; i < nat; ++i) {
            for (int k = 0; k < 3; ++k) xyz[3 * i + k] = geom_bohr(i, k);
            r4r2[i] = m_d4_generator->getSqrtZr4r2(m_atoms[i]);
        }
        for (size_t t = 0; t < atm_images.trans.size(); ++t)
            for (int k = 0; k < 3; ++k) trans[3 * t + k] = atm_images.trans[t][k];
        std::vector<double> e_at(nat), g_at(3 * nat), dcn_at(nat), v_at(m_want_virial ? 6 * nat : 0);
        if (m_gpu_scf->dispersionATMPeriodic(nat, xyz.data(), r4r2.data(), c6f.data(), dc6f.data(), m_cell_bohr.L,
                static_cast<int>(atm_images.trans.size()), trans.data(), atm_images.cutoff2, /*s9=*/5.0, /*a1=*/0.52,
                /*a2=*/5.0, /*alp=*/16.0, m_want_virial, e_at.data(), g_at.data(), dcn_at.data(),
                m_want_virial ? v_at.data() : nullptr)) {
            double esum = 0.0;
            double vs[6] = { 0, 0, 0, 0, 0, 0 };
            for (int a = 0; a < nat; ++a) {
                esum += e_at[a];
                for (int k = 0; k < 3; ++k) m_disp_gradient(a, k) += g_at[3 * a + k];
                m_disp_dEdcn(a) += dcn_at[a];
                if (m_want_virial)
                    for (int k = 0; k < 6; ++k) vs[k] += v_at[6 * a + k];
            }
            E += esum / 3.0;   // every translation class of triples once per member
            if (vir) {
                const int ia[6] = { 0, 0, 0, 1, 1, 2 }, ib[6] = { 0, 1, 2, 1, 2, 2 };
                for (int k = 0; k < 6; ++k) {
                    (*vir)(ia[k], ib[k]) += vs[k] / 3.0;
                    if (ia[k] != ib[k]) (*vir)(ib[k], ia[k]) += vs[k] / 3.0;
                }
            }
            atm_device = true;
        }
    }
    if (!atm_device)
        E += m_d4_evaluator->computeATMPeriodic(m_atoms, geom_bohr, atm_images, /*with_gradient=*/true,
            m_disp_gradient, m_disp_dEdcn, vir);
    const auto td2 = d4clk::now();

    // CN chain rule with the D4 CN's own derivative over images (as the cell-free path, d4_ncoord)
    if (m_disp_dEdcn.size() == nat) {
        curcuma::dispersion::addD4CovalentCNGradientPeriodic(m_atoms, geom_bohr, m_disp_dEdcn, m_disp_gradient, cn_images, vir);
        m_disp_dEdcn = Vector();
    }
    if (CurcumaLogger::get_verbosity() >= 3) {
        auto ms = [](d4clk::time_point a, d4clk::time_point b) { return std::chrono::duration<double, std::milli>(b - a).count(); };
        CurcumaLogger::info_fmt("  periodic D4: pairs {:8.2f} ms, ATM {:8.2f} ms ({}, host C6 {:.1f} ms), CN fold {:8.2f} ms ({} threads)",
                                ms(td0, td1), ms(td1, td2), atm_device ? "device" : "host", t_atm_c6, ms(td2, d4clk::now()),
                                effectiveIntraThreads(nat));
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

    using gmclk = std::chrono::steady_clock;   // G2 profile (verbosity 3), no numerical effect
    const auto tgm0 = gmclk::now();
    // 1. point-charge Ewald matrix of the atoms
    Eigen::MatrixXd pos(nat, 3);
    for (int i = 0; i < nat; ++i)
        for (int k = 0; k < 3; ++k) pos(i, k) = xyz_bohr[3 * i + k];
    const curcuma::ewald::EwaldParameters par = curcuma::ewald::chooseEwaldParameters(
        m_cell_bohr, m_pbc_ewald_tolerance, m_pbc_ewald_r_cut > 0.0 ? m_pbc_ewald_r_cut : -1.0);
    const curcuma::ewald::EwaldSum ewald(m_cell_bohr, par);
    std::vector<double> g(nsh);
    for (int is = 0; is < nsh; ++is)
        g[is] = coulomb::shell_hardness(meth, m_basis.z[m_basis.sh2at[is]], m_basis.ang_sh[is]);
    const PBCUtils::LatticeImageList images = PBCUtils::makeLatticeImageList(m_cell_bohr, m_pbc_coulomb_switch);
    const double rcut = m_pbc_coulomb_switch;

    // G2 (Nov 2026): Ewald matrix, KO remainder and shell expansion on the CUDA device when the periodic device SCF
    // is on (G1 wrapper); the host statements below run on any failure or with CURCUMA_PBC_GPU_HOST=gamma.
    if (m_gpu_scf && m_pbc_device_scf && !pbcGpuStageOnHost("gamma")) {
        const std::vector<double> ktab = ewald.kVectorTable();
        std::vector<double> tr(3 * images.trans.size());
        for (size_t t = 0; t < images.trans.size(); ++t)
            for (int k = 0; k < 3; ++k) tr[3 * t + k] = images.trans[t][k];
        m_gamma.resize(nsh, nsh);
        const bool gm_device = m_gpu_scf->ewaldGammaPeriodic(nat, xyz_bohr.data(), m_cell_bohr.L, ewald.parameters().alpha,
            ewald.volume(), ewald.parameters().r_cut, ewald.numberOfKVectors(), ktab.data(), nsh, m_basis.sh2at.data(), g.data(),
            meth == coulomb::Method::GFN1, gexp, static_cast<int>(images.trans.size()), tr.data(), rcut, m_gamma.data());
        if (gm_device) {
            m_pbc_ewald_alpha_used = par.alpha;
            m_pbc_ewald_nk_used = ewald.numberOfKVectors();
            if (CurcumaLogger::get_verbosity() >= 3)
                CurcumaLogger::info_fmt("  periodic gamma on the device: {:.1f} ms ({} k vectors)",
                    std::chrono::duration<double, std::milli>(gmclk::now() - tgm0).count(), m_pbc_ewald_nk_used);
            return;
        }
        if (CurcumaLogger::get_verbosity() >= 3)
            CurcumaLogger::info("  periodic gamma: device call failed, host loop");
    }

    const Eigen::MatrixXd A = ewald.matrix(pos);
    m_pbc_ewald_alpha_used = par.alpha;
    m_pbc_ewald_nk_used = ewald.numberOfKVectors();
    const auto tgm1 = gmclk::now();

    // 2. shell matrix: Ewald part per atom pair plus the short-range Klopman-Ohno remainder over images
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
    if (CurcumaLogger::get_verbosity() >= 3) {
        auto ms = [](gmclk::time_point a, gmclk::time_point b) { return std::chrono::duration<double, std::milli>(b - a).count(); };
        CurcumaLogger::info_fmt("  periodic gamma: Ewald matrix {:.1f} ms ({} k vectors), KO remainder {:.1f} ms",
                                ms(tgm0, tgm1), m_pbc_ewald_nk_used, ms(tgm1, gmclk::now()));
    }
}

/* ------------------------------------------------------------------ *
 *  P2b: GFN2 anisotropic electrostatics (dipole, quadrupole) with a    *
 *  cell. The Ewald parameters and translation lists follow tblite     *
 *  exactly (coulomb/ewald.f90 search_alpha / search_cutoff with the    *
 *  multipole terms, conv = 100 eps, box-shaped lattice points of       *
 *  cutoff.f90) because the direct sum also truncates the damping       *
 *  remainder (fdmp - 1) g3 ~ r^-6 of the damped kernels: a different   *
 *  direct cutoff is a different (converging) model value, documented   *
 *  in docs/WP-PERIODIC-XTB.md P2b.                                      *
 * ------------------------------------------------------------------ */
namespace {

constexpr double kTbliteEps = 1.4901161193847656e-08;   // sqrt(epsilon(0.0_wp))
constexpr double kSqrtPi = 1.7724538509055159;

inline double tbDirMp(double dist, double alpha)
{
    const double arg = alpha * dist;
    return (2.0 / kSqrtPi * arg * std::exp(-arg * arg) + std::erfc(arg)) / (dist * dist * dist);
}
inline double tbRecMp(double dist, double alpha, double vol)
{
    return 4.0 * M_PI * std::exp(-0.25 * dist * dist / (alpha * alpha)) / vol;
}

// tblite search_cutoff: smallest x (by bisection) with term(x) <= conv
template <class F>
double tbSearchCutoff(F term, double conv)
{
    double x = kTbliteEps;
    double y = term(x);
    while (y > conv && x <= 1.0e300) { x *= 2.0; y = term(x); }
    double xl = 0.5 * x, yl = term(xl), xr = x, yr = y;
    for (int iter = 0; iter < 30; ++iter) {
        if (yl - yr <= conv) break;
        x = 0.5 * (xl + xr);
        y = term(x);
        if (y >= conv) { xl = x; yl = y; } else { xr = x; yr = y; }
    }
    return x;
}

// tblite search_alpha for the multipole terms (lattice orthorhombic: rlen = 2 pi / L_max, dlen = L_min)
double tbSearchAlphaMultipole(const PBCUtils::OrthoCell& cell)
{
    const double vol = cell.L[0] * cell.L[1] * cell.L[2];
    const double dlen = cell.minLength();
    const double rlen = 2.0 * M_PI / std::max(cell.L[0], std::max(cell.L[1], cell.L[2]));
    auto diff = [&](double a) {
        return (tbRecMp(4.0 * rlen, a, vol) - tbRecMp(5.0 * rlen, a, vol))
             - (tbDirMp(2.0 * dlen, a) - tbDirMp(3.0 * dlen, a));
    };
    const double tol = kTbliteEps;
    double alpha = kTbliteEps;
    double d = diff(alpha);
    while (d < -tol && alpha <= 1.0e300) { alpha *= 2.0; d = diff(alpha); }
    if (alpha > 1.0e300 || alpha == kTbliteEps) return 0.25;
    double alpl = 0.5 * alpha;
    while (d < tol && alpha <= 1.0e300) { alpha *= 2.0; d = diff(alpha); }
    if (alpha > 1.0e300) return 0.25;
    double alpr = alpha;
    alpha = 0.5 * (alpl + alpr);
    d = diff(alpha);
    int ibs = 0;
    while (std::abs(d) > tol && ibs <= 30) {
        if (d < 0) alpl = alpha; else alpr = alpha;
        alpha = 0.5 * (alpl + alpr);
        d = diff(alpha);
        ++ibs;
    }
    if (ibs > 30) return 0.25;
    return alpha;
}

// tblite get_lattice_points: every n with |n_a| <= ceil(rthr / L_a) (box, no spherical filter); origin optional
std::vector<std::array<double, 3>> tbBoxPoints(const double L[3], double rthr, bool with_origin)
{
    int rep[3];
    for (int a = 0; a < 3; ++a) rep[a] = static_cast<int>(std::ceil(rthr / L[a]));
    std::vector<std::array<double, 3>> pts;
    for (int ix = -rep[0]; ix <= rep[0]; ++ix)
        for (int iy = -rep[1]; iy <= rep[1]; ++iy)
            for (int iz = -rep[2]; iz <= rep[2]; ++iz) {
                if (!with_origin && ix == 0 && iy == 0 && iz == 0) continue;
                pts.push_back({ ix * L[0], iy * L[1], iz * L[2] });
            }
    return pts;
}

} // namespace

void XTB::buildMultipoleInteractionMatricesPeriodic() const
{
    using namespace curcuma::xtb::gfn2_params;
    const int nat = m_atomcount;
    const std::vector<double> xyz = geometryBohr(m_geometry, nat);
    for (int k = 0; k < 3; ++k) m_mp_amat_sd[k] = Eigen::MatrixXd::Zero(nat, nat);
    for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b) m_mp_amat_dd[a][b] = Eigen::MatrixXd::Zero(nat, nat);
    for (int k = 0; k < 6; ++k) m_mp_amat_sq[k] = Eigen::MatrixXd::Zero(nat, nat);

    const double vol = m_cell_bohr.L[0] * m_cell_bohr.L[1] * m_cell_bohr.L[2];
    double alpha = tbSearchAlphaMultipole(m_cell_bohr);
    double conv = 100.0 * kTbliteEps;
    if (const char* e = std::getenv("CURCUMA_PBC_MP_ALPHA")) alpha = std::atof(e);   // convergence experiments
    if (const char* e = std::getenv("CURCUMA_PBC_MP_CONV")) conv = std::atof(e);
    // charge-dipole reciprocal term: -grad of (4 pi/V) sum_k w_k cos(k.r) over the full k set is sum_k w_k k sin(k.r)
    // (factor 1); tblite's get_amat_sdq_rec_3d uses 2 k sin (factor 2) over the same full set. Experiment switch.
    double sd_rec_factor = m_pbc_mp_tblite_factors ? 2.0 : 1.0;   // -xtb.pbc_multipole_tblite_factors
    if (const char* e = std::getenv("CURCUMA_PBC_MP_SD_FACTOR")) sd_rec_factor = std::atof(e);
    double sq_rec_factor = m_pbc_mp_tblite_factors ? 1.0 : -1.0 / 3.0;   // tblite: +1 (see the comment at the reciprocal loop)
    if (const char* e = std::getenv("CURCUMA_PBC_MP_SQ_FACTOR")) sq_rec_factor = std::atof(e);
    auto build = [&](double alpha) {
    const double dir_cut = tbSearchCutoff([&](double x) { return tbDirMp(x, alpha); }, conv);
    const double rec_cut = tbSearchCutoff([&](double x) { return tbRecMp(x, alpha, vol); }, conv);
    const std::vector<std::array<double, 3>> dtrans = tbBoxPoints(m_cell_bohr.L, dir_cut, true);
    const double Lrec[3] = { 2.0 * M_PI / m_cell_bohr.L[0], 2.0 * M_PI / m_cell_bohr.L[1], 2.0 * M_PI / m_cell_bohr.L[2] };
    const std::vector<std::array<double, 3>> rtrans = tbBoxPoints(Lrec, rec_cut, false);
    if (CurcumaLogger::get_verbosity() >= 3)
        CurcumaLogger::info_fmt("  periodic multipoles: alpha {:.4f}, direct cutoff {:.2f} Bohr ({} translations), "
                                "reciprocal cutoff {:.3f} ({} k vectors)", alpha, dir_cut, dtrans.size(), rec_cut, rtrans.size());
    m_mp_pbc_alpha = alpha;      // the gradient (P3) uses the same split
    m_mp_pbc_dtrans = dtrans;
    m_mp_pbc_rtrans = rtrans;

    const double fac = 4.0 * M_PI / vol;
    const double alp2 = alpha * alpha;
    const double eps_r = kTbliteEps;

    // reciprocal weights per k (independent of the pair)
    std::vector<double> kexp(rtrans.size());
    for (size_t t = 0; t < rtrans.size(); ++t) {
        const double g2 = rtrans[t][0] * rtrans[t][0] + rtrans[t][1] * rtrans[t][1] + rtrans[t][2] * rtrans[t][2];
        kexp[t] = (g2 < eps_r) ? 0.0 : fac * std::exp(-0.25 * g2 / alp2) / g2;
    }

    // G3 (Nov 2026): the pair loop on the CUDA device when the periodic device SCF is on (G1 wrapper); it writes
    // every off-diagonal and diagonal element of the 18 blocks like the loop below, which runs on any failure.
    bool mp_device = false;
    if (m_gpu_scf && m_pbc_device_scf && !pbcGpuStageOnHost("mpmat")) {
        std::vector<double> dt(3 * dtrans.size()), rt(3 * rtrans.size());
        for (size_t t = 0; t < dtrans.size(); ++t)
            for (int k = 0; k < 3; ++k) dt[3 * t + k] = dtrans[t][k];
        for (size_t t = 0; t < rtrans.size(); ++t)
            for (int k = 0; k < 3; ++k) rt[3 * t + k] = rtrans[t][k];
        double* outs[18];
        for (int k = 0; k < 3; ++k) outs[k] = m_mp_amat_sd[k].data();
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b) outs[3 + a * 3 + b] = m_mp_amat_dd[a][b].data();
        for (int k = 0; k < 6; ++k) outs[12 + k] = m_mp_amat_sq[k].data();
        mp_device = m_gpu_scf->multipoleMatricesPeriodic(nat, xyz.data(), m_mp_mrad.data(), m_cell_bohr.L,
            static_cast<int>(dtrans.size()), dt.data(), static_cast<int>(rtrans.size()), rt.data(), kexp.data(),
            alpha, mp_dmp3, mp_dmp5, sd_rec_factor, sq_rec_factor, fac, outs);
        if (CurcumaLogger::get_verbosity() >= 3)
            CurcumaLogger::info_fmt("  periodic multipole matrices on the {}", mp_device ? "device" : "host (device call failed)");
    }
    const int mp_threads = mp_device ? 1 : effectiveIntraThreads(nat);
    if (!mp_device)
    parallelStripes(mp_threads, [&](int tid, int nth) {
    for (int iat = tid; iat < nat; iat += nth) {
        for (int jat = 0; jat <= iat; ++jat) {
            // minimum image of x_i - x_j (tblite: Wigner-Seitz image; identical for non-degenerate positions)
            double v0x = xyz[3 * iat] - xyz[3 * jat], v0y = xyz[3 * iat + 1] - xyz[3 * jat + 1], v0z = xyz[3 * iat + 2] - xyz[3 * jat + 2];
            PBCUtils::applyMinimumImageOrtho(m_cell_bohr, v0x, v0y, v0z);
            const double rr = 0.5 * (m_mp_mrad[jat] + m_mp_mrad[iat]);

            double sd[3] = { 0, 0, 0 }, dd[3][3] = { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } }, sq[6] = { 0, 0, 0, 0, 0, 0 };
            // reciprocal part (tblite get_amat_sdq_rec_3d), k.T = 2 pi n: the image choice does not matter
            dd[0][0] += fac / 6.0; dd[1][1] += fac / 6.0; dd[2][2] += fac / 6.0;
            sq[0] -= fac / 9.0; sq[2] -= fac / 9.0; sq[5] -= fac / 9.0;
            for (size_t t = 0; t < rtrans.size(); ++t) {
                if (kexp[t] == 0.0) continue;
                const double kx = rtrans[t][0], ky = rtrans[t][1], kz = rtrans[t][2];
                const double gv = v0x * kx + v0y * ky + v0z * kz;
                const double sink = std::sin(gv) * kexp[t];
                const double cosk = std::cos(gv) * kexp[t];
                sd[0] += sd_rec_factor * kx * sink; sd[1] += sd_rec_factor * ky * sink; sd[2] += sd_rec_factor * kz * sink;
                const double kv[3] = { kx, ky, kz };
                for (int a = 0; a < 3; ++a)
                    for (int b = 0; b < 3; ++b) dd[a][b] += kv[a] * kv[b] * cosk;
                // charge-quadrupole: the vacuum kernel v_a v_b / r^5 is (1/3) d_a d_b (1/r) for a traceless
                // quadrupole, so its reciprocal part is (1/3) d_a d_b sum_k w_k cos(k.r) = -(1/3) w_k k_a k_b cos(k.r)
                // (the direct part of tblite is exactly (1/3) d_a d_b erfc(a r)/r, see the derivation in
                // docs/WP-PERIODIC-XTB.md P2b). tblite's get_amat_sdq_rec_3d adds + k_a k_b cos (factor -3 off),
                // which made the block alpha dependent (probe CURCUMA_PBC_MP_DEBUG: 6.6e-3 of 3.3e-2 at L = 12).
                const double cq = sq_rec_factor * cosk;
                sq[0] += kx * kx * cq; sq[1] += 2.0 * kx * ky * cq; sq[2] += ky * ky * cq;
                sq[3] += 2.0 * kx * kz * cq; sq[4] += 2.0 * ky * kz * cq; sq[5] += kz * kz * cq;
            }
            // direct part (tblite get_amat_sdq_dir_3d): every translation of the box, damped kernel minus the
            // Ewald-screened bare kernel
            for (const auto& T : dtrans) {
                const double vx = v0x + T[0], vy = v0y + T[1], vz = v0z + T[2];
                const double r1 = std::sqrt(vx * vx + vy * vy + vz * vz);
                if (r1 < eps_r) continue;
                const double g1 = 1.0 / r1, g3 = g1 * g1 * g1, g5 = g3 * g1 * g1;
                const double fdmp3 = 1.0 / (1.0 + 6.0 * std::pow(rr / r1, mp_dmp3));
                const double fdmp5 = 1.0 / (1.0 + 6.0 * std::pow(rr / r1, mp_dmp5));
                const double arg = r1 * alpha, arg2 = arg * arg;
                const double expt = std::exp(-arg2) / kSqrtPi;
                const double erft = -std::erf(arg) * g1;
                const double e1 = g1 * g1 * (erft + 2.0 * expt * alpha);
                const double e2 = g1 * g1 * (e1 + 4.0 * expt * alp2 * alpha / 3.0);
                const double tmp3 = fdmp3 * g3 + e1;
                const double t5 = g5 * fdmp5 + e2;
                const double t3 = fdmp5 * g3 + e1;
                sd[0] += vx * tmp3; sd[1] += vy * tmp3; sd[2] += vz * tmp3;
                const double vv[3] = { vx, vy, vz };
                for (int a = 0; a < 3; ++a)
                    for (int b = 0; b < 3; ++b) dd[a][b] += ((a == b) ? t3 : 0.0) - vv[a] * vv[b] * 3.0 * t5;
                sq[0] += vx * vx * t5 - t3 / 3.0; sq[1] += 2.0 * vx * vy * t5; sq[2] += vy * vy * t5 - t3 / 3.0;
                sq[3] += 2.0 * vx * vz * t5; sq[4] += 2.0 * vy * vz * t5; sq[5] += vz * vz * t5 - t3 / 3.0;
            }
            // amat(:, jat, iat) with vec = x_i - x_j (target j, source i); the mirrored element by (anti)symmetry
            for (int k = 0; k < 3; ++k) {
                m_mp_amat_sd[k](jat, iat) = (iat == jat) ? 0.0 : sd[k];
                m_mp_amat_sd[k](iat, jat) = (iat == jat) ? 0.0 : -sd[k];
            }
            for (int a = 0; a < 3; ++a)
                for (int b = 0; b < 3; ++b) {
                    m_mp_amat_dd[a][b](jat, iat) = dd[a][b];
                    m_mp_amat_dd[a][b](iat, jat) = dd[a][b];
                }
            for (int k = 0; k < 6; ++k) {
                m_mp_amat_sq[k](jat, iat) = sq[k];
                m_mp_amat_sq[k](iat, jat) = sq[k];
            }
        }
    }
    });
    // self terms (tblite): dipole -2/3 alpha^3/sqrt(pi) (x2 on the diagonal), quadrupole trace 4/9 alpha^3/sqrt(pi)
    const double a3 = alpha * alpha * alpha / kSqrtPi;
    for (int iat = 0; iat < nat; ++iat) {
        for (int k = 0; k < 3; ++k) m_mp_amat_dd[k][k](iat, iat) += 2.0 * (-2.0 / 3.0 * a3);
        const double rq = 4.0 / 9.0 * a3;
        m_mp_amat_sq[0](iat, iat) += rq; m_mp_amat_sq[2](iat, iat) += rq; m_mp_amat_sq[5](iat, iat) += rq;
    }
    };   // build(alpha)
    build(alpha);
    if (std::getenv("CURCUMA_PBC_MP_DEBUG")) {
        // alpha-independence probe per block: rebuild with 1.5 alpha, print the largest element change
        const auto sd0 = m_mp_amat_sd;
        const auto dd0 = m_mp_amat_dd;
        const auto sq0 = m_mp_amat_sq;
        build(1.5 * alpha);
        double dsd = 0, ddd = 0, dsq = 0, msd = 0, mdd = 0, msq = 0;
        for (int k = 0; k < 3; ++k) { dsd = std::max(dsd, (m_mp_amat_sd[k] - sd0[k]).cwiseAbs().maxCoeff()); msd = std::max(msd, sd0[k].cwiseAbs().maxCoeff()); }
        for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) { ddd = std::max(ddd, (m_mp_amat_dd[a][b] - dd0[a][b]).cwiseAbs().maxCoeff()); mdd = std::max(mdd, dd0[a][b].cwiseAbs().maxCoeff()); }
        for (int k = 0; k < 6; ++k) { dsq = std::max(dsq, (m_mp_amat_sq[k] - sq0[k]).cwiseAbs().maxCoeff()); msq = std::max(msq, sq0[k].cwiseAbs().maxCoeff()); }
        CurcumaLogger::result(fmt::format("MP_DEBUG alpha {:.4f} vs {:.4f}: max|dSD| {:.3e} (max {:.3e}), max|dDD| {:.3e} (max {:.3e}), max|dSQ| {:.3e} (max {:.3e})",
                                          alpha, 1.5 * alpha, dsd, msd, ddd, mdd, dsq, msq));
        build(alpha);
    }
    m_mp_amat_deferred = false;
}

void XTB::setupMultipolePeriodic()
{
    using namespace curcuma::xtb::gfn2_params;
    const int nat = m_atomcount;
    const int nao = m_basis.nao;
    const int nsh = m_basis.nsh;
    const std::vector<double> xyz = geometryBohr(m_geometry, nat);
    const PBCUtils::LatticeImageList images = PBCUtils::makeLatticeImageList(m_cell_bohr, pbcIntegralCutoffBohr());
    std::vector<CGTO::Shell> shells(nsh);
    for (int ish = 0; ish < nsh; ++ish) shells[ish] = as_cgto_shell(m_basis.cgto[ish]);
    using mpclk = std::chrono::steady_clock;   // G2 profile (verbosity 3 setup table), no numerical effect
    auto mpms = [](mpclk::time_point a, mpclk::time_point b) { return std::chrono::duration<double, std::milli>(b - a).count(); };
    const auto tmp0 = mpclk::now();

    // 1 + 2: AO dipole and raw quadrupole integrals, per image shifted to the image position of the column atom
    // (origin at atom(nu) + T), summed over the images; the traceless transform is linear and applied at the end.
    for (int k = 0; k < 3; ++k) m_dp_int[k] = Eigen::MatrixXd::Zero(nao, nao);
    std::array<Eigen::MatrixXd, 6> qraw;
    for (int k = 0; k < 6; ++k) qraw[k] = Eigen::MatrixXd::Zero(nao, nao);

    const int mp_threads = effectiveIntraThreads(nsh);
    parallelStripes(mp_threads, [&](int tid, int nth) {
    for (int ish_a = tid; ish_a < nsh; ish_a += nth) {
        const int ang_a = m_basis.ang_sh[ish_a];
        const int iat = m_basis.sh2at[ish_a];
        const int ia0 = m_basis.iao_sh[ish_a];
        const int nca = m_basis.nao_sh[ish_a];
        const CGTO::Shell& sh_a = shells[ish_a];
        const double* xi = &xyz[3 * iat];
        for (int ish_b = 0; ish_b < nsh; ++ish_b) {
            const int ang_b = m_basis.ang_sh[ish_b];
            const int jat = m_basis.sh2at[ish_b];
            const int jb0 = m_basis.iao_sh[ish_b];
            const int ncb = m_basis.nao_sh[ish_b];
            const CGTO::Shell& sh_b = shells[ish_b];
            const bool dpair = (ang_a >= 2 || ang_b >= 2);

            auto addBlock = [&](const double* xj) {
                double bS[9], bD[27], bQ[54];
                double Db[5 * 5 * 3], Qb[5 * 5 * 6], Sb[6 * 6];
                if (!dpair) {
                    MI::cgto_multipole_block(sh_a, ang_a, sh_b, ang_b, xi[0], xi[1], xi[2], xj[0], xj[1], xj[2], bS, bD, bQ);
                } else {
                    sphericalMultipoleBlock(sh_a, ang_a, sh_b, ang_b, xi[0], xi[1], xi[2], xj[0], xj[1], xj[2], Db, Qb, 5);
                    sphericalOverlapBlock(sh_a, ang_a, sh_b, ang_b, xi[0], xi[1], xi[2], xj[0], xj[1], xj[2], Sb, 6);
                }
                const double Rx = xj[0], Ry = xj[1], Rz = xj[2];
                for (int ia = 0; ia < nca; ++ia) {
                    const int mu = ia0 + ia;
                    for (int jb = 0; jb < ncb; ++jb) {
                        const int nu = jb0 + jb;
                        double Smn, dx, dy, dz, q[6];
                        if (!dpair) {
                            const int c = ia * ncb + jb;
                            Smn = bS[c]; dx = bD[c * 3]; dy = bD[c * 3 + 1]; dz = bD[c * 3 + 2];
                            for (int k = 0; k < 6; ++k) q[k] = bQ[c * 6 + k];
                        } else {
                            const int o = ia * 5 + jb;
                            Smn = Sb[ia * 6 + jb]; dx = Db[o * 3]; dy = Db[o * 3 + 1]; dz = Db[o * 3 + 2];
                            for (int k = 0; k < 6; ++k) q[k] = Qb[o * 6 + k];
                        }
                        m_dp_int[0](mu, nu) += dx - Rx * Smn;
                        m_dp_int[1](mu, nu) += dy - Ry * Smn;
                        m_dp_int[2](mu, nu) += dz - Rz * Smn;
                        qraw[0](mu, nu) += q[0] - 2 * Rx * dx + Rx * Rx * Smn;
                        qraw[1](mu, nu) += q[1] - Rx * dy - Ry * dx + Rx * Ry * Smn;
                        qraw[2](mu, nu) += q[2] - 2 * Ry * dy + Ry * Ry * Smn;
                        qraw[3](mu, nu) += q[3] - Rx * dz - Rz * dx + Rx * Rz * Smn;
                        qraw[4](mu, nu) += q[4] - Ry * dz - Rz * dy + Ry * Rz * Smn;
                        qraw[5](mu, nu) += q[5] - 2 * Rz * dz + Rz * Rz * Smn;
                    }
                }
            };
            if (iat == jat) addBlock(&xyz[3 * jat]);   // T = 0 on-atom block
            PBCUtils::forEachImage(images, xi, &xyz[3 * jat], iat == jat,
                [&](double dx, double dy, double dz, double) {
                    const double xj_img[3] = { xi[0] - dx, xi[1] - dy, xi[2] - dz };
                    addBlock(xj_img);
                });
        }
    }
    });
    for (int k = 0; k < 6; ++k) m_qp_int[k] = Eigen::MatrixXd::Zero(nao, nao);
    for (int mu = 0; mu < nao; ++mu)
        for (int nu = 0; nu < nao; ++nu) {
            const double qxx = qraw[0](mu, nu), qxy = qraw[1](mu, nu), qyy = qraw[2](mu, nu);
            const double qxz = qraw[3](mu, nu), qyz = qraw[4](mu, nu), qzz = qraw[5](mu, nu);
            const double tr = 0.5 * (qxx + qyy + qzz);
            m_qp_int[0](mu, nu) = 1.5 * qxx - tr;
            m_qp_int[1](mu, nu) = 1.5 * qxy;
            m_qp_int[2](mu, nu) = 1.5 * qyy - tr;
            m_qp_int[3](mu, nu) = 1.5 * qxz;
            m_qp_int[4](mu, nu) = 1.5 * qyz;
            m_qp_int[5](mu, nu) = 1.5 * qzz - tr;
        }

    const auto tmp_ao = mpclk::now();
    // 3 + 4: periodic GFN2 CN (the same double-exponential CN as the self-energies) and the damping radii
    const Vector cn = computeCoordinationNumbers();
    m_mp_mrad.resize(nat);
    m_mp_dkernel.resize(nat);
    m_mp_qkernel.resize(nat);
    for (int i = 0; i < nat; ++i) {
        const int z = m_atoms[i];
        const double vcn = p_vcn[z - 1];
        const double rad = p_rad[z - 1];
        const double arg = cn(i) - vcn - mp_shift;
        const double t1 = std::exp(-mp_kexp * arg);
        const double t2 = (mp_rmax - rad) / (1.0 + t1);
        m_mp_mrad[i] = rad + t2;
        m_mp_dkernel[i] = p_dkernel[z - 1];
        m_mp_qkernel[i] = p_qkernel[z - 1];
    }
    const auto tmp_cn = mpclk::now();
    // 5: interaction matrices by Ewald summation
    buildMultipoleInteractionMatricesPeriodic();
    m_mp_t_d_block = m_mp_t_shift = 0.0;
    m_mp_t_ao_ints = mpms(tmp0, tmp_ao);
    m_mp_t_cn_mrad = mpms(tmp_ao, tmp_cn);
    m_mp_t_amat = mpms(tmp_cn, mpclk::now());
    m_mp_initialized = true;
    ++m_mp_version;
}

/* ------------------------------------------------------------------ *
 *  P3: nuclear gradient with a cell (Eh/Bohr into m_gradient)          *
 *  Sections as in calculateGradient() (xtb_gradient.cpp), each over    *
 *  the lattice images; a self image (i, i, T != 0) has no Cartesian    *
 *  derivative but contributes to dE/dCN (H0 of the self-image block).  *
 * ------------------------------------------------------------------ */
void XTB::addHalogenBondGradientPeriodic(Matrix& gradient, Eigen::Matrix3d* virial) const
{
    // P4: every pair of (gradient(p) += G, gradient(q) -= G) along d = x_p - x_q adds -G (x) d to the virial
    auto vir = [&](double gx, double gy, double gz, const double d[3]) {
        if (!virial) return;
        const double G[3] = { gx, gy, gz };
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b) (*virial)(a, b) -= G[a] * d[b];
    };
    if (m_method != MethodType::GFN1) return;
    const int nat = m_atomcount;
    if (nat < 2 || gradient.rows() != nat || gradient.cols() != 3) return;
    const std::vector<double> xyz = geometryBohr(m_geometry, nat);
    const PBCUtils::LatticeImageList images = PBCUtils::makeLatticeImageList(m_cell_bohr, kPbcHalogenCutoffBohr);
    constexpr double xb_alp = 6.0, xb_lj2 = 6.0;

    for (int x = 0; x < nat; ++x) {
        const int zx = m_atoms[x];
        if (!xb_is_halogen(zx)) continue;
        const double cc = gfn1_params::halogen_bond[zx - 1];
        if (cc == 0.0) continue;
        double best = std::numeric_limits<double>::max();
        double dxb[3] = { 0.0, 0.0, 0.0 };
        int b = -1;
        for (int m = 0; m < nat; ++m) {
            PBCUtils::forEachImage(images, &xyz[3 * x], &xyz[3 * m], m == x,
                [&](double dx, double dy, double dz, double r2) {
                    if (r2 < best) { best = r2; dxb[0] = -dx; dxb[1] = -dy; dxb[2] = -dz; b = m; }
                });
        }
        if (b < 0) continue;
        const double d2bx = best;
        const double rbx = std::sqrt(d2bx) + 1.0e-18;
        const double damping = gfn1_params::halogen_damping;
        for (int a = 0; a < nat; ++a) {
            if (!xb_is_acceptor(m_atoms[a])) continue;
            const double r0ax = gfn1_params::halogen_radscale * (atomic_rad_au(zx) + atomic_rad_au(m_atoms[a]));
            PBCUtils::forEachImage(images, &xyz[3 * x], &xyz[3 * a], false,
                [&](double dx, double dy, double dz, double d2ax) {
                    const double dxa[3] = { -dx, -dy, -dz };
                    const double dba[3] = { dxa[0] - dxb[0], dxa[1] - dxb[1], dxa[2] - dxb[2] };
                    const double d2ab = dba[0]*dba[0] + dba[1]*dba[1] + dba[2]*dba[2];
                    const double rax = std::sqrt(d2ax) + 1.0e-18;
                    const double xy = std::sqrt(d2bx * d2ax);
                    const double term = (d2bx + d2ax - d2ab) / xy;
                    const double aterm = std::pow(0.5 - 0.25 * term, xb_alp);
                    const double t14 = std::pow(r0ax / rax, xb_lj2);
                    const double numerator = t14 * t14 - damping * t14;
                    const double denominator = 1.0 + t14 * t14;
                    const double termlj = numerator / denominator;
                    double dtermlj = 2.0 * xb_lj2 * numerator * t14 * t14 / (rax * denominator * denominator);
                    dtermlj += xb_lj2 * t14 * (damping - 2.0 * t14) / (rax * denominator);
                    dtermlj *= aterm * cc / rax;
                    for (int k = 0; k < 3; ++k) { gradient(a, k) += dtermlj * dxa[k]; gradient(x, k) -= dtermlj * dxa[k]; }
                    vir(dtermlj * dxa[0], dtermlj * dxa[1], dtermlj * dxa[2], dxa);
                    double prefactor = -0.25 * xb_alp * std::pow(0.5 - 0.25 * term, xb_alp - 1.0);
                    prefactor *= cc * termlj;
                    double dcosterm = (2.0 / rbx - term / rax) * prefactor / rax;
                    for (int k = 0; k < 3; ++k) { gradient(a, k) += dcosterm * dxa[k]; gradient(x, k) -= dcosterm * dxa[k]; }
                    vir(dcosterm * dxa[0], dcosterm * dxa[1], dcosterm * dxa[2], dxa);
                    dcosterm = (2.0 / rax - term / rbx) * prefactor / rbx;
                    for (int k = 0; k < 3; ++k) { gradient(b, k) += dcosterm * dxb[k]; gradient(x, k) -= dcosterm * dxb[k]; }
                    vir(dcosterm * dxb[0], dcosterm * dxb[1], dcosterm * dxb[2], dxb);
                    const double t13 = 2.0 * prefactor / xy;
                    for (int k = 0; k < 3; ++k) { gradient(a, k) -= t13 * dba[k]; gradient(b, k) += t13 * dba[k]; }
                    vir(-t13 * dba[0], -t13 * dba[1], -t13 * dba[2], dba);
                });
        }
    }
}

// Port of tblite get_multipole_gradient_3d / get_damat_sdq_{dir,rec}_3d with the P2b reciprocal factors
// (charge-dipole 1, charge-quadrupole -1/3). vec = x_j - x_i + T (tblite gradient convention), gradient(i) += dG,
// gradient(j) -= dG; dE/d(mrad) into dEdcn through dmrad/dCN as in section 5 of calculateGradient().
void XTB::addMultipoleGradientPeriodic(Matrix& gradient, Vector& dEdcn, Eigen::Matrix3d* virial) const
{
    using namespace curcuma::xtb::gfn2_params;
    const int nat = m_atomcount;
    if (!m_mp_initialized || m_mp_pbc_alpha <= 0.0) return;
    const std::vector<double> xyz = geometryBohr(m_geometry, nat);
    const double vol = m_cell_bohr.L[0] * m_cell_bohr.L[1] * m_cell_bohr.L[2];
    const double alp = m_mp_pbc_alpha, alp2 = alp * alp;
    const double fac = 4.0 * M_PI / vol;
    const double eps_r = kTbliteEps;
    double sd_f = m_pbc_mp_tblite_factors ? 2.0 : 1.0, sq_f = m_pbc_mp_tblite_factors ? 1.0 : -1.0 / 3.0;
    if (const char* e = std::getenv("CURCUMA_PBC_MP_SD_FACTOR")) sd_f = std::atof(e);
    if (const char* e = std::getenv("CURCUMA_PBC_MP_SQ_FACTOR")) sq_f = std::atof(e);
    const auto& q = m_wfn.q_at;
    const auto& dp = m_wfn.dp_at;
    const auto& qp = m_wfn.qp_at;

    std::vector<double> kexp(m_mp_pbc_rtrans.size());
    for (size_t t = 0; t < m_mp_pbc_rtrans.size(); ++t) {
        const auto& k = m_mp_pbc_rtrans[t];
        const double g2 = k[0] * k[0] + k[1] * k[1] + k[2] * k[2];
        kexp[t] = (g2 < eps_r) ? 0.0 : fac * std::exp(-0.25 * g2 / alp2) / g2;
    }

    Vector dEdr_mp = Vector::Zero(nat);
    // pair (iat, jat): direct images and reciprocal part; iat == jat only for dE/d(mrad)
    // P4: dEde[a][b] = dE/d(eps_ab) of the pair at fixed fractional coordinates (direct images: dE/dvec (x) vec with
    // vec = x_j + T - x_i, i.e. -dG (x) vec; reciprocal: k_a -> k_a - eps_ba k_b, w_k -> w_k [1 - tr eps +
    // k_a k_b eps_ab (2/k^2 + 1/(2 alpha^2))], k.r invariant; the 1/V constants scale with -tr eps). The self pair
    // (iat == jat) carries the single-ordering weights (1/2 for dipole-dipole and the direct part).
    auto pairTerms = [&](int iat, int jat, double dG[3], double& dE, double dEde[3][3]) {
        dG[0] = dG[1] = dG[2] = 0.0;
        dE = 0.0;
        const bool want_v = (dEde != nullptr);
        if (want_v) for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) dEde[a][b] = 0.0;
        const double wself = (iat == jat) ? 0.5 : 1.0;
        double v0x = xyz[3 * jat] - xyz[3 * iat], v0y = xyz[3 * jat + 1] - xyz[3 * iat + 1], v0z = xyz[3 * jat + 2] - xyz[3 * iat + 2];
        PBCUtils::applyMinimumImageOrtho(m_cell_bohr, v0x, v0y, v0z);
        const double rr = 0.5 * (m_mp_mrad[jat] + m_mp_mrad[iat]);
        const double qi = q(iat), qj = q(jat);
        const double mi[3] = { dp(0, iat), dp(1, iat), dp(2, iat) };
        const double mj[3] = { dp(0, jat), dp(1, jat), dp(2, jat) };
        double ti[6], tj[6];
        for (int k = 0; k < 6; ++k) { ti[k] = qp(k, iat); tj[k] = qp(k, jat); }
        const double dpidpj = mi[0] * mj[0] + mi[1] * mj[1] + mi[2] * mj[2];

        // reciprocal (gradient for iat != jat; virial for every pair incl. iat == jat)
        if (iat != jat || want_v) {
            for (size_t t = 0; t < m_mp_pbc_rtrans.size(); ++t) {
                if (kexp[t] == 0.0) continue;
                const double* kv = m_mp_pbc_rtrans[t].data();
                const double gv = v0x * kv[0] + v0y * kv[1] + v0z * kv[2];
                const double cosk = std::cos(gv) * kexp[t];
                const double sink = std::sin(gv) * kexp[t];
                const double dpiv = kv[0] * mi[0] + kv[1] * mi[1] + kv[2] * mi[2];
                const double dpjv = kv[0] * mj[0] + kv[1] * mj[1] + kv[2] * mj[2];
                const double dpiqj = dpiv * qj, qidpj = dpjv * qi;
                const double qiqpj = qi * (tj[0]*kv[0]*kv[0] + tj[2]*kv[1]*kv[1] + tj[5]*kv[2]*kv[2]
                                          + 2*tj[1]*kv[0]*kv[1] + 2*tj[3]*kv[0]*kv[2] + 2*tj[4]*kv[1]*kv[2]);
                const double qpiqj = qj * (ti[0]*kv[0]*kv[0] + ti[2]*kv[1]*kv[1] + ti[5]*kv[2]*kv[2]
                                          + 2*ti[1]*kv[0]*kv[1] + 2*ti[3]*kv[0]*kv[2] + 2*ti[4]*kv[1]*kv[2]);
                if (iat != jat)
                    for (int c = 0; c < 3; ++c) {
                        dG[c] -= sd_f * kv[c] * cosk * (dpiqj - qidpj);        // tblite: 2 * ...
                        dG[c] += kv[c] * sink * dpiv * dpjv;
                        dG[c] += sq_f * kv[c] * sink * (qiqpj + qpiqj);        // tblite: +1 * ...
                    }
                if (want_v) {
                    const double g2 = kv[0] * kv[0] + kv[1] * kv[1] + kv[2] * kv[2];
                    const double X = 2.0 / g2 + 0.5 / alp2;
                    const double w = kexp[t];
                    const double cg = cosk / w, sg = sink / w;   // cos(gv), sin(gv)
                    // theta.k vectors
                    const double tjk[3] = { tj[0]*kv[0] + tj[1]*kv[1] + tj[3]*kv[2], tj[1]*kv[0] + tj[2]*kv[1] + tj[4]*kv[2], tj[3]*kv[0] + tj[4]*kv[1] + tj[5]*kv[2] };
                    const double tik[3] = { ti[0]*kv[0] + ti[1]*kv[1] + ti[3]*kv[2], ti[1]*kv[0] + ti[2]*kv[1] + ti[4]*kv[2], ti[3]*kv[0] + ti[4]*kv[1] + ti[5]*kv[2] };
                    const double A = dpiqj - qidpj;              // charge-dipole scalar (odd in k)
                    const double B = (iat != jat) ? (qiqpj + qpiqj) : qpiqj;   // charge-quadrupole scalar
                    const double DDw = (iat != jat) ? 1.0 : 0.5;  // dipole-dipole single-ordering weight for the self pair
                    for (int a = 0; a < 3; ++a)
                        for (int b = 0; b < 3; ++b) {
                            const double wprime = w * (-(a == b ? 1.0 : 0.0) + kv[a] * kv[b] * X);
                            double dE_ab = 0.0;
                            // charge-dipole: with v0 = x_j - x_i the pair energy is E = +sd_f sin(gv) w A (the matrix
                            // uses x_i - x_j and the energy -mu_i.SD(j,i) q_j; sin is odd), consistent with the
                            // gradient dG -= sd_f k cos A above
                            if (iat != jat)
                                dE_ab += sd_f * sg * (wprime * A + w * (-(mi[a] * qj - mj[a] * qi) * kv[b]));
                            // dipole-dipole: E = DDw cos(gv) w (k.mu_i)(k.mu_j)
                            dE_ab += DDw * cg * (wprime * dpiv * dpjv + w * (-mi[a] * kv[b] * dpjv - dpiv * mj[a] * kv[b]));
                            // charge-quadrupole: E = sq_f cos(gv) w B
                            const double Bprime = (iat != jat) ? -2.0 * (qi * tjk[a] + qj * tik[a]) * kv[b]
                                                                : -2.0 * qi * tik[a] * kv[b];
                            dE_ab += sq_f * cg * (wprime * B + w * Bprime);
                            dEde[a][b] += dE_ab;
                        }
                }
            }
            if (want_v) {   // dipole-dipole constant (fac/6) delta_ab, fac ~ 1/V: dE/d(eps_ab) = -(fac/6) mu_i.mu_j delta_ab
                const double cdd = (iat != jat) ? (fac / 6.0) * dpidpj : 0.5 * (fac / 6.0) * dpidpj;
                for (int a = 0; a < 3; ++a) dEde[a][a] -= cdd;
            }
        }
        // direct
        for (const auto& T : m_mp_pbc_dtrans) {
            const double vec[3] = { v0x + T[0], v0y + T[1], v0z + T[2] };
            const double r1 = std::sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
            if (r1 < eps_r) continue;
            const double r2 = r1 * r1, g1 = 1.0 / r1, g3 = g1 * g1 * g1, g5 = g3 * g1 * g1, g7 = g5 * g1 * g1;
            const double arg = r1 * alp, arg2 = arg * arg;
            const double erft = -std::erf(arg) * g1;
            const double expt = std::exp(-arg2) / kSqrtPi;
            const double e1 = g1 * g1 * (erft + expt * (2 * alp2) / alp);
            const double e2 = g1 * g1 * (e1 + expt * (2 * alp2) * (2 * alp2) / (3 * alp));
            const double e3 = g1 * g1 * (e2 + expt * (2 * alp2) * (2 * alp2) * (2 * alp2) / (15 * alp));
            const double fdmp3 = 1.0 / (1.0 + 6.0 * std::pow(rr * g1, mp_dmp3));
            const double ddmp3 = -3 * fdmp3 - mp_dmp3 * fdmp3 * (fdmp3 - 1.0);
            const double fdmp5 = 1.0 / (1.0 + 6.0 * std::pow(rr * g1, mp_dmp5));
            const double ddmp5 = -5 * fdmp5 - mp_dmp5 * (fdmp5 * fdmp5 - fdmp5);
            const double dpiv = vec[0] * mi[0] + vec[1] * mi[1] + vec[2] * mi[2];
            const double dpjv = vec[0] * mj[0] + vec[1] * mj[1] + vec[2] * mj[2];
            const double dpiqj = dpiv * qj, qidpj = dpjv * qi;
            const double edd = dpidpj * r2 - 3 * dpjv * dpiv;
            // charge - dipole
            double g_sd[3], tab[3][3];
            for (int c = 0; c < 3; ++c)
                g_sd[c] = -(ddmp3 * g5) * vec[c] * (dpiqj - qidpj) + fdmp3 * g3 * (qi * mj[c] - qj * mi[c]);
            for (int bb = 0; bb < 3; ++bb) {
                for (int a = 0; a < 3; ++a) tab[a][bb] = 3 * vec[a] * vec[bb] * e2;
                tab[bb][bb] -= e1;
            }
            for (int k = 0; k < 3; ++k)
                g_sd[k] += qj * (tab[0][k] * mi[0] + tab[1][k] * mi[1] + tab[2][k] * mi[2])
                         - qi * (tab[0][k] * mj[0] + tab[1][k] * mj[1] + tab[2][k] * mj[2]);
            dE += 3.0 * (dpiqj - qidpj) * mp_dmp3 * fdmp3 * g3 * (fdmp3 / rr) * std::pow(rr * g1, mp_dmp3);
            // dipole - dipole
            double g_dd[3], tabc[3][3][3];
            for (int c = 0; c < 3; ++c)
                g_dd[c] = -2.0 * fdmp5 * g5 * dpidpj * vec[c] + 3.0 * fdmp5 * g5 * (dpiv * mj[c] + dpjv * mi[c])
                        - edd * ddmp5 * g7 * vec[c];
            for (int c = 0; c < 3; ++c) {
                for (int bb = 0; bb < 3; ++bb)
                    for (int a = 0; a < 3; ++a) tabc[a][bb][c] = -15 * vec[a] * vec[bb] * vec[c] * e3;
                for (int a = 0; a < 3; ++a) {
                    tabc[a][a][c] += 3 * e2 * vec[c];
                    tabc[c][a][c] += 3 * e2 * vec[a];
                    tabc[a][c][c] += 3 * e2 * vec[a];
                }
            }
            for (int k = 0; k < 3; ++k)
                for (int a = 0; a < 3; ++a)
                    g_dd[k] += mi[a] * (tabc[0][a][k] * mj[0] + tabc[1][a][k] * mj[1] + tabc[2][a][k] * mj[2]);
            dE += 3.0 * edd * mp_dmp5 * fdmp5 * g5 * (fdmp5 / rr) * std::pow(rr * g1, mp_dmp5);
            // charge - quadrupole
            const double eq = 2 * (qj * ti[1] + tj[1] * qi) * vec[0] * vec[1] + 2 * (qj * ti[3] + tj[3] * qi) * vec[0] * vec[2]
                            + 2 * (qj * ti[4] + tj[4] * qi) * vec[1] * vec[2] + (qj * ti[0] + tj[0] * qi) * vec[0] * vec[0]
                            + (qj * ti[2] + tj[2] * qi) * vec[1] * vec[1] + (qj * ti[5] + tj[5] * qi) * vec[2] * vec[2];
            const double tjv[3] = { vec[0]*tj[0] + vec[1]*tj[1] + vec[2]*tj[3], vec[0]*tj[1] + vec[1]*tj[2] + vec[2]*tj[4], vec[0]*tj[3] + vec[1]*tj[4] + vec[2]*tj[5] };
            const double tiv[3] = { vec[0]*ti[0] + vec[1]*ti[1] + vec[2]*ti[3], vec[0]*ti[1] + vec[1]*ti[2] + vec[2]*ti[4], vec[0]*ti[3] + vec[1]*ti[4] + vec[2]*ti[5] };
            double g_sq[3];
            for (int c = 0; c < 3; ++c)
                g_sq[c] = -eq * ddmp5 * g7 * vec[c] - 2.0 * fdmp5 * g5 * qi * tjv[c] - 2.0 * fdmp5 * g5 * qj * tiv[c];
            for (int k = 0; k < 3; ++k)
                g_sq[k] += (-qi * (tabc[0][0][k] * tj[0] + 2 * tabc[1][0][k] * tj[1] + 2 * tabc[2][0][k] * tj[3]
                                   + tabc[1][1][k] * tj[2] + 2 * tabc[2][1][k] * tj[4] + tabc[2][2][k] * tj[5])
                            - qj * (tabc[0][0][k] * ti[0] + 2 * tabc[1][0][k] * ti[1] + 2 * tabc[2][0][k] * ti[3]
                                   + tabc[1][1][k] * ti[2] + 2 * tabc[2][1][k] * ti[4] + tabc[2][2][k] * ti[5])) / 3.0;
            dE += eq * 3.0 * mp_dmp5 * fdmp5 * g5 * fdmp5 / rr * std::pow(rr * g1, mp_dmp5);
            if (iat != jat)
                for (int c = 0; c < 3; ++c) dG[c] += g_sd[c] + g_dd[c] + g_sq[c];
            if (want_v)   // dE/d(eps_ab) = (dE/dvec_a) vec_b = -(g_a) vec_b (g = dE/dx_i of this image)
                for (int a = 0; a < 3; ++a)
                    for (int b = 0; b < 3; ++b) dEde[a][b] -= wself * (g_sd[a] + g_dd[a] + g_sq[a]) * vec[b];
        }
    };

    // G5a (Nov 2026): the pair loop on the CUDA device when the periodic device SCF is on (fixed-point sums,
    // order independent); on any failure the host blocks below run.
    bool mpg_device = false;
    if (m_gpu_scf && m_pbc_device_scf && !pbcGpuStageOnHost("mpgrad")) {
        std::vector<double> qv(nat), dpv(3 * nat), qpv(6 * nat);
        for (int i = 0; i < nat; ++i) {
            qv[i] = q(i);
            for (int k = 0; k < 3; ++k) dpv[3 * i + k] = dp(k, i);
            for (int k = 0; k < 6; ++k) qpv[6 * i + k] = qp(k, i);
        }
        std::vector<double> dt(3 * m_mp_pbc_dtrans.size()), rt(3 * m_mp_pbc_rtrans.size());
        for (size_t t = 0; t < m_mp_pbc_dtrans.size(); ++t)
            for (int k = 0; k < 3; ++k) dt[3 * t + k] = m_mp_pbc_dtrans[t][k];
        for (size_t t = 0; t < m_mp_pbc_rtrans.size(); ++t)
            for (int k = 0; k < 3; ++k) rt[3 * t + k] = m_mp_pbc_rtrans[t][k];
        std::vector<double> g3(3 * nat), der(nat), vir9(9);
        mpg_device = m_gpu_scf->multipoleGradientPeriodic(nat, xyz.data(), m_mp_mrad.data(), qv.data(), dpv.data(),
            qpv.data(), m_cell_bohr.L, static_cast<int>(m_mp_pbc_dtrans.size()), dt.data(),
            static_cast<int>(m_mp_pbc_rtrans.size()), rt.data(), kexp.data(), alp, fac, sd_f, sq_f, mp_dmp3, mp_dmp5,
            virial != nullptr, g3.data(), der.data(), vir9.data());
        if (mpg_device) {
            for (int i = 0; i < nat; ++i) {
                for (int k = 0; k < 3; ++k) gradient(i, k) += g3[3 * i + k];
                dEdr_mp(i) += der[i];
            }
            if (virial)
                for (int a = 0; a < 3; ++a)
                    for (int b = 0; b < 3; ++b) (*virial)(a, b) += vir9[3 * a + b];
        }
    }
    // Stage T: fixed row blocks (iat % NB == blk) with per-block accumulators reduced in block order, so the
    // result does not depend on the thread count; the same (iat, jat) order within a block as the plain loop.
    const int NB = mpg_device ? 0 : (nat < 32 ? nat : 32);
    std::vector<Matrix> g_b(NB, Matrix::Zero(nat, 3));
    std::vector<Vector> dEdr_b(NB, Vector::Zero(nat));
    std::vector<Eigen::Matrix3d> v_b(virial ? NB : 0, Eigen::Matrix3d::Zero());
    if (!mpg_device)
    parallelStripes(effectiveIntraThreads(nat), [&](int tid, int nth) {
        for (int blk = tid; blk < NB; blk += nth) {
            Matrix& gb = g_b[blk];
            Vector& eb = dEdr_b[blk];
            double dEde[3][3];
            for (int iat = blk; iat < nat; iat += NB) {
                for (int jat = 0; jat <= iat; ++jat) {
                    double dG[3], dE;
                    pairTerms(iat, jat, dG, dE, virial ? dEde : nullptr);
                    eb(iat) += dE;
                    if (jat != iat) {
                        eb(jat) += dE;
                        for (int c = 0; c < 3; ++c) { gb(iat, c) += dG[c]; gb(jat, c) -= dG[c]; }
                    }
                    if (virial)
                        for (int a = 0; a < 3; ++a)
                            for (int b = 0; b < 3; ++b) v_b[blk](a, b) -= dEde[a][b];   // Xi = -dE/d(eps)
                }
            }
        }
    });
    for (int blk = 0; blk < NB; ++blk) {
        gradient += g_b[blk];
        dEdr_mp += dEdr_b[blk];
        if (virial) *virial += v_b[blk];
    }
    const Vector cn = computeCoordinationNumbers();
    for (int i = 0; i < nat; ++i) {
        const int zi = m_atoms[i];
        const double rad = p_rad[zi - 1], vcn = p_vcn[zi - 1];
        const double arg = cn(i) - vcn - mp_shift;
        const double t1 = std::exp(-mp_kexp * arg);
        const double t2 = (mp_rmax - rad) / (1.0 + t1);
        const double dmrdcn = -t2 * mp_kexp * t1 / (1.0 + t1);
        dEdcn(i) += dEdr_mp(i) * dmrdcn;
    }
}

void XTB::calculateGradientPeriodic()
{
    ensureHostIntegrals(__func__);
    const int nat = m_atomcount;
    const int nsh = m_basis.nsh;
    const int nao = m_basis.nao;
    m_gradient = Matrix::Zero(nat, 3);
    const std::vector<double> xyz = geometryBohr(m_geometry, nat);
    // P4: atomic virial Xi_ab = -dE/d(eps_ab) at fixed fractional coordinates; every image term contributes
    // -(dE/dd_a) d_b with its pair vector d (self images included: |T| scales with the cell), Ewald terms their
    // own strain derivatives. Symmetrised at the end.
    const bool want_vir = m_want_virial;
    Eigen::Matrix3d Xi = Eigen::Matrix3d::Zero();
    m_virial_valid = false;
    auto virPair = [&](Eigen::Matrix3d& X, double gx, double gy, double gz, double dx, double dy, double dz) {
        const double G[3] = { gx, gy, gz }, d[3] = { dx, dy, dz };
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b) X(a, b) -= G[a] * d[b];
    };

    using gclk = std::chrono::steady_clock;   // G2 profile (verbosity 3), no numerical effect
    const auto tg0 = gclk::now();
    // energy-weighted density W (as calculateGradient)
    const int nocc_orbs = static_cast<int>(std::round(m_wfn.nocc / 2.0));
    Matrix W;
    if (m_wfn.W_valid && m_wfn.W.rows() == nao && m_wfn.W.cols() == nao) {
        W = m_wfn.W;
    } else if (m_wfn.focc.size() == nao && m_wfn.C.cols() == nao) {
        int ncol = 0;
        for (int i = 0; i < nao; ++i)
            if (m_wfn.focc(i) > 1.0e-12) ncol = i + 1;
        const Eigen::VectorXd we = (m_wfn.focc.head(ncol).array() * m_wfn.eps.head(ncol).array()).matrix();
        // G5b (Nov 2026): the GEMM on the CUDA device from its resident eigenvectors (same C as m_wfn.C) when the
        // periodic device SCF ran; host GEMM on any failure or with CURCUMA_PBC_GPU_HOST=wmat.
        bool w_device = false;
        if (m_gpu_scf && m_pbc_device_scf && m_pbc_dev_wfn && ncol > 0 && !pbcGpuStageOnHost("wmat")) {
            Eigen::MatrixXd Wcm(nao, nao);
            w_device = m_gpu_scf->energyWeightedDensity(we.data(), ncol, nao, Wcm.data());
            if (w_device) W = Wcm;   // element (i, j) kept (Matrix is row-major)
            if (CurcumaLogger::get_verbosity() >= 3)
                CurcumaLogger::info_fmt("  periodic gradient: W on the {}", w_device ? "device" : "host (device call failed)");
        }
        if (!w_device) {
            const auto Cocc = m_wfn.C.leftCols(ncol);
            W = (Cocc * we.asDiagonal()) * Cocc.transpose();
        }
    } else {
        const auto Cocc = m_wfn.C.leftCols(nocc_orbs);
        const Eigen::VectorXd w2 = 2.0 * m_wfn.eps.head(nocc_orbs);
        W = (Cocc * w2.asDiagonal()) * Cocc.transpose();
    }
    Vector se;
    getSelfEnergies(m_coordination_numbers, se);
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

    const auto tg_w = gclk::now();
    // 1. repulsion over images (j < i; self images have no derivative)
    {
        const PBCUtils::LatticeImageList images = PBCUtils::makeLatticeImageList(m_cell_bohr, kPbcRepulsionCutoffBohr);
        const bool g1m = (m_method == MethodType::GFN1);
        const double kexp = g1m ? gfn1_params::rep_kexp : gfn2_params::rep_kexp;
        const double rexp = g1m ? gfn1_params::rep_rexp : gfn2_params::rep_rexp;
        // G5d (Nov 2026): the pair loop on the CUDA device (host loop on failure or CURCUMA_PBC_GPU_HOST=rep)
        bool rep_device = false;
        if (m_gpu_scf && m_pbc_device_scf && !pbcGpuStageOnHost("rep")) {
            const PbcRepulsionArrays ra = pbcRepulsionArrays(m_atoms, g1m, images);
            std::vector<double> g(3 * nat, 0.0), v9(9, 0.0);
            rep_device = m_gpu_scf->repulsionPeriodic(nat, xyz.data(), ra.alpha.data(), ra.zeff.data(),
                ra.light.data(), kexp, g1m ? kexp : gfn2_params::rep_kexp_light, rexp, m_cell_bohr.L,
                static_cast<int>(images.trans.size()), ra.trans.data(), images.cutoff2, true, want_vir, nullptr,
                g.data(), v9.data());
            if (rep_device) {
                for (int i = 0; i < nat; ++i)
                    for (int k = 0; k < 3; ++k) m_gradient(i, k) += g[3 * i + k];
                if (want_vir)
                    for (int a = 0; a < 3; ++a)
                        for (int b = 0; b < 3; ++b) Xi(a, b) += v9[3 * a + b];
            }
        }
        if (!rep_device)
        for (int i = 0; i < nat; ++i) {
            const int zi = m_atoms[i];
            const double alfi = g1m ? gfn1_params::rep_alpha[zi - 1] : gfn2_params::rep_alpha[zi - 1];
            const double zeffi = g1m ? gfn1_params::rep_zeff[zi - 1] : gfn2_params::rep_zeff[zi - 1];
            for (int j = 0; j <= (want_vir ? i : i - 1); ++j) {   // self images only for the virial (weight 1/2)
                const int zj = m_atoms[j];
                const double alfj = g1m ? gfn1_params::rep_alpha[zj - 1] : gfn2_params::rep_alpha[zj - 1];
                const double zeffj = g1m ? gfn1_params::rep_zeff[zj - 1] : gfn2_params::rep_zeff[zj - 1];
                const double alpha_pair = std::sqrt(alfi * alfj);
                double kexp_pair = kexp;
                if (!g1m && zi <= 2 && zj <= 2) kexp_pair = gfn2_params::rep_kexp_light;
                const bool same = (i == j);
                PBCUtils::forEachImage(images, &xyz[3 * i], &xyz[3 * j], same,
                    [&](double dx, double dy, double dz, double r2) {
                        const double r = std::sqrt(r2);
                        const double r_kexp = std::pow(r, kexp_pair);
                        const double E_pair = zeffi * zeffj / std::pow(r, rexp) * std::exp(-alpha_pair * r_kexp);
                        const double dEdr = -(rexp / r + alpha_pair * kexp_pair * std::pow(r, kexp_pair - 1.0)) * E_pair;
                        const double fx = dEdr * dx / r, fy = dEdr * dy / r, fz = dEdr * dz / r;
                        if (!same) {
                            m_gradient(i, 0) += fx; m_gradient(j, 0) -= fx;
                            m_gradient(i, 1) += fy; m_gradient(j, 1) -= fy;
                            m_gradient(i, 2) += fz; m_gradient(j, 2) -= fz;
                        }
                        if (want_vir) {
                            const double w = same ? 0.5 : 1.0;
                            virPair(Xi, w * fx, w * fy, w * fz, dx, dy, dz);
                        }
                    });
            }
        }
    }
    // 1b. halogen bond
    addHalogenBondGradientPeriodic(m_gradient, want_vir ? &Xi : nullptr);

    const auto tg_rep = gclk::now();
    // 2a. on-site CN term
    Vector dEdcn = Vector::Zero(nat);
    for (int ish = 0; ish < nsh; ++ish) {
        const int iat = m_basis.sh2at[ish];
        for (int ia = 0; ia < m_basis.nao_sh[ish]; ++ia) {
            const int mu = m_basis.iao_sh[ish] + ia;
            dEdcn(iat) += (-m_h0.kcn[ish]) * m_wfn.P(mu, mu);
        }
    }

    // 2b. H0 / Pulay over images: pairs iat < jat (all images) and iat == jat (T != 0, dE/dCN only)
    {
        const PBCUtils::LatticeImageList images = PBCUtils::makeLatticeImageList(m_cell_bohr, pbcIntegralCutoffBohr());
        const bool do_mp = (m_method == MethodType::GFN2 && m_mp_initialized && !m_mp_grad_off);
        const int grad_threads = effectiveIntraThreads(nat);
        std::vector<Matrix> grad_parts(grad_threads, Matrix::Zero(nat, 3));
        std::vector<Vector> dEdcn_parts(grad_threads, Vector::Zero(nat));
        std::vector<Eigen::Matrix3d> vir_parts(grad_threads, Eigen::Matrix3d::Zero());
        parallelStripes(grad_threads, [&](int tid, int nth) {
        Matrix& g_loc = grad_parts[tid];
        Vector& dEdcn_loc = dEdcn_parts[tid];
        Eigen::Matrix3d& v_loc = vir_parts[tid];
        for (int iat = tid; iat < nat; iat += nth) {
            for (int jat = iat; jat < nat; ++jat) {
                const bool same = (iat == jat);
                const int zi = m_basis.z[iat], zj = m_basis.z[jat];
                const double rad_sum = atomic_rad_au(zi) + atomic_rad_au(zj);
                const double* xi = &xyz[3 * iat];
                PBCUtils::forEachImage(images, xi, &xyz[3 * jat], same,
                    [&](double dx_ij, double dy_ij, double dz_ij, double r2) {
                        // d = x_i - (x_j + T): image of j at x_i - d
                        const double xj[3] = { xi[0] - dx_ij, xi[1] - dy_ij, xi[2] - dz_ij };
                        const double r = std::sqrt(r2);
                        const double rr = std::sqrt(r / rad_sum);
                        for (int ia = 0; ia < m_basis.nsh_at[iat]; ++ia) {
                            const int ish_a = m_basis.ish_at[iat] + ia;
                            const int ia_start = m_basis.iao_sh[ish_a], ia_nao = m_basis.nao_sh[ish_a];
                            const double pi_a = 1.0 + m_h0.shpoly[ish_a] * rr;
                            const double zeta_a = m_basis.cgto[ish_a].slater_exp;
                            const CGTO::Shell& sh_a = shells[ish_a];
                            const int ang_a = m_basis.ang_sh[ish_a];
                            for (int ib = 0; ib < m_basis.nsh_at[jat]; ++ib) {
                                const int ish_b = m_basis.ish_at[jat] + ib;
                                const int jb_start = m_basis.iao_sh[ish_b], jb_nao = m_basis.nao_sh[ish_b];
                                const double pi_b = 1.0 + m_h0.shpoly[ish_b] * rr;
                                const double zeta_b = m_basis.cgto[ish_b].slater_exp;
                                const CGTO::Shell& sh_b = shells[ish_b];
                                const int ang_b = m_basis.ang_sh[ish_b];
                                double hs;
                                if (m_method == MethodType::GFN1) {
                                    const bool vi = valence[ish_a], vj = valence[ish_b];
                                    const double den = std::pow(pauling_en[zi - 1] - pauling_en[zj - 1], 2);
                                    if (vi && vj) hs = gfn1::kpair(zi, zj) * gfn1::kshell(ang_a, ang_b) * (1.0 + gfn1::enscale * den);
                                    else if (vi && !vj) hs = 0.5 * (gfn1::kshell(ang_a, ang_a) + gfn1::kdiff);
                                    else if (!vi && vj) hs = 0.5 * (gfn1::kshell(ang_b, ang_b) + gfn1::kdiff);
                                    else hs = gfn1::kdiff;
                                } else {
                                    const double den = std::pow(pauling_en[zi - 1] - pauling_en[zj - 1], 2);
                                    const double enp = 1.0 + gfn2::enscale * den;
                                    const double km = gfn2::kpair(zi, zj) * gfn2::kshell(ang_a, ang_b) * enp;
                                    const double zij = std::pow(2.0 * std::sqrt(zeta_a * zeta_b) / (zeta_a + zeta_b), gfn2::wexp);
                                    hs = zij * km;
                                }
                                const double h_factor = hs * (pi_a * pi_b);
                                const double avg_eps = 0.5 * (se(ish_a) + se(ish_b));
                                const double h_av = avg_eps * h_factor;
                                // image overlap block (the Gamma-point S is a sum; the shpoly and CN terms need S_T)
                                double blk_s[9];
                                CGTO::cgto_overlap_block(sh_a, ang_a, sh_b, ang_b, xi[0], xi[1], xi[2], xj[0], xj[1], xj[2], blk_s);
                                double cn_sum = 0.0;
                                // self image (iat == jat, T != 0): the loop visits both AO orderings and every T, so the
                                // pair formulas below (built for the half-visited iat < jat case) carry a weight 1/2 here;
                                // no Cartesian derivative (both centres move together), but the integrals depend on |T|
                                // and enter the virial (P4) and dE/dCN
                                const double wself = same ? 0.5 : 1.0;
                                if (same && !want_vir) {
                                    for (int iao = 0; iao < ia_nao; ++iao)
                                        for (int jao = 0; jao < jb_nao; ++jao)
                                            cn_sum += m_wfn.P(ia_start + iao, jb_start + jao) * blk_s[iao * jb_nao + jao];
                                    dEdcn_loc(iat) += wself * (-m_h0.kcn[ish_a]) * h_factor * cn_sum;
                                    dEdcn_loc(jat) += wself * (-m_h0.kcn[ish_b]) * h_factor * cn_sum;
                                    continue;
                                }
                                const double dlog_pi_dr_r = (m_h0.shpoly[ish_a] / pi_a + m_h0.shpoly[ish_b] / pi_b) * rr / (2.0 * r2);
                                const double v_a = m_pot.v_sh(ish_a) + m_pot.v_at(iat);
                                const double v_b = m_pot.v_sh(ish_b) + m_pot.v_at(jat);
                                CGTO::OverlapPrimPair pp[CGTO::kMaxPrimPairs];
                                const int npp = CGTO::buildOverlapPrimPairs(sh_a, sh_b, xi[0], xi[1], xi[2], xj[0], xj[1], xj[2], pp);
                                double G_sval[3] = { 0.0, 0.0, 0.0 };
                                double G_shpoly_scalar = 0.0;
                                for (int iao = 0; iao < ia_nao; ++iao) {
                                    const int mu = ia_start + iao;
                                    const int t_a = ao_to_type(ang_a, iao);
                                    if (t_a < 0) continue;
                                    for (int jao = 0; jao < jb_nao; ++jao) {
                                        const int nu = jb_start + jao;
                                        const int t_b = ao_to_type(ang_b, jao);
                                        if (t_b < 0) continue;
                                        const double Pmn = m_wfn.P(mu, nu);
                                        const double Wmn = W(mu, nu);
                                        const double S_T = blk_s[iao * jb_nao + jao];
                                        const double H0_T = h_av * S_T;
                                        double dS[3];
                                        if (npp >= 0) CGTO::cgto_overlap_grad_pre(pp, npp, t_a, t_b, dS);
                                        else CGTO::cgto_overlap_grad(sh_a, sh_b, xi[0], xi[1], xi[2], xj[0], xj[1], xj[2], t_a, t_b, dS);
                                        const double sval = 2.0 * Pmn * h_av - 2.0 * Wmn - Pmn * (v_a + v_b);
                                        G_sval[0] += sval * dS[0]; G_sval[1] += sval * dS[1]; G_sval[2] += sval * dS[2];
                                        if (do_mp) {
                                            double S_mp, D_mp[3], Q_mp[6], dD_dA[3][3], dD_dB[3][3], dQ_dA[3][6], dQ_dB[3][6];
                                            MI::cgto_multipole_grad_transformed(sh_a, sh_b, xi[0], xi[1], xi[2], xj[0], xj[1], xj[2],
                                                                                t_a, t_b, S_mp, D_mp, Q_mp, dD_dA, dD_dB, dQ_dA, dQ_dB);
                                            const double dR[3] = { xj[0] - xi[0], xj[1] - xi[1], xj[2] - xi[2] };
                                            static const int qa6[6] = { 0, 0, 1, 0, 1, 2 };
                                            static const int qb6[6] = { 0, 1, 1, 2, 2, 2 };
                                            for (int l = 0; l < 3; ++l) {
                                                double term = 0.0;
                                                for (int k = 0; k < 3; ++k) {
                                                    term += dD_dA[l][k] * m_pot.v_dp(k, jat);
                                                    const double dDiat = dD_dA[l][k] + dR[k] * dS[l] - (k == l ? S_T : 0.0);
                                                    term += dDiat * m_pot.v_dp(k, iat);
                                                }
                                                double dq[6];
                                                for (int qq = 0; qq < 6; ++qq) {
                                                    const int a = qa6[qq], b = qb6[qq];
                                                    dq[qq] = -(b == l ? D_mp[a] : 0.0) + dR[b] * dD_dA[l][a]
                                                           - (a == l ? D_mp[b] : 0.0) + dR[a] * dD_dA[l][b]
                                                           + (-(a == l ? dR[b] : 0.0) - (b == l ? dR[a] : 0.0)) * S_T
                                                           + dR[a] * dR[b] * dS[l];
                                                }
                                                const double dtr_c = 0.5 * (dq[0] + dq[2] + dq[5]);
                                                for (int qq = 0; qq < 6; ++qq) {
                                                    const bool is_diag = (qa6[qq] == qb6[qq]);
                                                    term += dQ_dA[l][qq] * m_pot.v_qp(qq, jat);
                                                    const double dQiat = dQ_dA[l][qq] + 1.5 * dq[qq] - (is_diag ? dtr_c : 0.0);
                                                    term += dQiat * m_pot.v_qp(qq, iat);
                                                }
                                                G_sval[l] -= Pmn * term;
                                            }
                                        }
                                        G_shpoly_scalar += 2.0 * Pmn * H0_T * dlog_pi_dr_r;
                                        cn_sum += Pmn * S_T;
                                    }
                                }
                                const double Gx = wself * (G_sval[0] + G_shpoly_scalar * dx_ij);
                                const double Gy = wself * (G_sval[1] + G_shpoly_scalar * dy_ij);
                                const double Gz = wself * (G_sval[2] + G_shpoly_scalar * dz_ij);
                                if (!same) {
                                    g_loc(iat, 0) += Gx; g_loc(jat, 0) -= Gx;
                                    g_loc(iat, 1) += Gy; g_loc(jat, 1) -= Gy;
                                    g_loc(iat, 2) += Gz; g_loc(jat, 2) -= Gz;
                                }
                                if (want_vir) virPair(v_loc, Gx, Gy, Gz, dx_ij, dy_ij, dz_ij);
                                dEdcn_loc(iat) += wself * (-m_h0.kcn[ish_a]) * h_factor * cn_sum;
                                dEdcn_loc(jat) += wself * (-m_h0.kcn[ish_b]) * h_factor * cn_sum;
                            }
                        }
                    });
            }
        }
        });
        for (int t = 0; t < grad_threads; ++t) {
            m_gradient += grad_parts[t];
            dEdcn += dEdcn_parts[t];
            Xi += vir_parts[t];
        }
    }

    const auto tg_h0 = gclk::now();
    auto tg_ew = tg_h0;
    // 3. isotropic Coulomb: point-charge Ewald gradient at q_at plus the Klopman-Ohno remainder over images
    if (!m_coulomb_grad_off) {
        Eigen::MatrixXd pos(nat, 3);
        for (int i = 0; i < nat; ++i)
            for (int k = 0; k < 3; ++k) pos(i, k) = xyz[3 * i + k];
        const curcuma::ewald::EwaldParameters par = curcuma::ewald::chooseEwaldParameters(
            m_cell_bohr, m_pbc_ewald_tolerance, m_pbc_ewald_r_cut > 0.0 ? m_pbc_ewald_r_cut : -1.0);
        const curcuma::ewald::EwaldSum ewald(m_cell_bohr, par);
        m_gradient += ewald.gradient(pos, m_wfn.q_at);
        if (want_vir) {   // P4: real, reciprocal and background strain derivatives of the point-charge Ewald sum (N1)
            curcuma::ewald::EwaldGeometryCache cache;
            ewald.buildGeometryCache(pos, nullptr, cache, static_cast<size_t>(1) << 30);
            ewald.realSpaceVirialRange(cache, m_wfn.q_at, 0, static_cast<int>(cache.pairs.size()), Xi);
            ewald.reciprocalVirialRange(cache.table, m_wfn.q_at, 0, ewald.numberOfKVectors(), Xi);
            curcuma::ewald::EwaldResult bg;
            ewald.selfAndBackground(m_wfn.q_at, nullptr, false, bg);
            Xi += curcuma::ewald::EwaldSum::backgroundVirial(bg.terms);
        }
        tg_ew = gclk::now();

        const double gexp = (m_method == MethodType::GFN1) ? gfn1_params::gexp : gfn2_params::gexp;
        const coulomb::Method cm = (m_method == MethodType::GFN1) ? coulomb::Method::GFN1 : coulomb::Method::GFN2;
        std::vector<double> g_sh(nsh);
        for (int is = 0; is < nsh; ++is) g_sh[is] = coulomb::shell_hardness(cm, m_basis.z[m_basis.sh2at[is]], m_basis.ang_sh[is]);
        const PBCUtils::LatticeImageList images = PBCUtils::makeLatticeImageList(m_cell_bohr, m_pbc_coulomb_switch);
        const double rcut = m_pbc_coulomb_switch;
        for (int is = 0; is < nsh; ++is) {
            const int iat = m_basis.sh2at[is];
            for (int js = 0; js <= is; ++js) {
                const int jat = m_basis.sh2at[js];
                const bool same_atom = (iat == jat);
                if (same_atom && !want_vir) continue;   // on-site and self images: no Cartesian derivative
                // self images of an on-site shell pair: energy 1/2 q q sum_T (is == js) or q q sum_T (is != js, both
                // orderings) -> weights for the virial
                const double wgt = same_atom ? ((is == js) ? 0.5 : 1.0) : 1.0;
                const double gam_bar = (is == js) ? g_sh[is] : coulomb::average(cm, g_sh[is], g_sh[js]);
                const double gam_pow = std::pow(gam_bar, -gexp);
                const double qq = m_wfn.q_sh(is) * m_wfn.q_sh(js);
                PBCUtils::forEachImage(images, &xyz[3 * iat], &xyz[3 * jat], same_atom,
                    [&](double dx, double dy, double dz, double r2) {
                        const double r1 = std::sqrt(r2);
                        const double f = koSwitch(r1, rcut);
                        if (f == 0.0) return;
                        const double gamma = std::pow(std::pow(r1, gexp) + gam_pow, -1.0 / gexp);
                        const double dgamma_dr = -std::pow(r1, gexp - 2.0) * std::pow(gamma, gexp + 1.0) * r1;
                        // d/dr [ f (gamma - 1/r) ] = f' (gamma - 1/r) + f (gamma' + 1/r^2)
                        double df = 0.0;
                        if (r1 > rcut - 1.0 && r1 < rcut) {
                            const double x = r1 - (rcut - 1.0);
                            df = -30.0 * x * x * x * x + 60.0 * x * x * x - 30.0 * x * x;
                        }
                        const double dko = df * (gamma - 1.0 / r1) + f * (dgamma_dr + 1.0 / r2);
                        const double force = qq * dko / r1;
                        if (!same_atom) {
                            m_gradient(iat, 0) += force * dx; m_gradient(jat, 0) -= force * dx;
                            m_gradient(iat, 1) += force * dy; m_gradient(jat, 1) -= force * dy;
                            m_gradient(iat, 2) += force * dz; m_gradient(jat, 2) -= force * dz;
                        }
                        if (want_vir) virPair(Xi, wgt * force * dx, wgt * force * dy, wgt * force * dz, dx, dy, dz);
                    });
            }
        }
    }

    const auto tg_coul = gclk::now();
    // 5. GFN2 multipoles
    if (m_method == MethodType::GFN2 && m_mp_initialized && !m_mp_grad_off)
        addMultipoleGradientPeriodic(m_gradient, dEdcn, want_vir ? &Xi : nullptr);
    const auto tg_mp = gclk::now();

    // 3b. dispersion (periodic D3/D4, CN chain rule already folded; its virial was accumulated with it)
    if (m_disp_gradient_valid && !m_d4_grad_off && m_disp_gradient.rows() == nat && m_disp_gradient.cols() == 3) {
        m_gradient += m_disp_gradient;
        if (m_disp_dEdcn.size() == nat) dEdcn += m_disp_dEdcn;
        if (want_vir) Xi += m_disp_virial;
    }

    // 4. CN chain rule of the H0 coordination number over images (j < i)
    if (!m_cnchain_off) {
        const PBCUtils::LatticeImageList images = PBCUtils::makeLatticeImageList(m_cell_bohr, kPbcCnCutoffBohr);
        std::vector<double> rcov(nat);
        for (int i = 0; i < nat; ++i) rcov[i] = covalent_rad_d3_au(m_atoms[i]);
        // G5c (Nov 2026): the pair loop on the CUDA device when the periodic device SCF is on; the host loop below
        // runs on any failure or with CURCUMA_PBC_GPU_HOST=cnchain.
        bool cn_device = false;
        if (m_gpu_scf && m_pbc_device_scf && !pbcGpuStageOnHost("cnchain")) {
            std::vector<double> tr(3 * images.trans.size()), g(3 * nat, 0.0), v9(9, 0.0);
            for (size_t t = 0; t < images.trans.size(); ++t)
                for (int k = 0; k < 3; ++k) tr[3 * t + k] = images.trans[t][k];
            cn_device = m_gpu_scf->cnChainGradientPeriodic(nat, xyz.data(), rcov.data(), dEdcn.data(), m_cell_bohr.L,
                static_cast<int>(images.trans.size()), tr.data(), images.cutoff2, m_method == MethodType::GFN1,
                want_vir, g.data(), v9.data());
            if (cn_device) {
                for (int i = 0; i < nat; ++i)
                    for (int k = 0; k < 3; ++k) m_gradient(i, k) += g[3 * i + k];
                if (want_vir)
                    for (int a = 0; a < 3; ++a)
                        for (int b = 0; b < 3; ++b) Xi(a, b) += v9[3 * a + b];
            }
            if (CurcumaLogger::get_verbosity() >= 3)
                CurcumaLogger::info_fmt("  periodic CN chain on the {}", cn_device ? "device" : "host (device call failed)");
        }
        if (!cn_device)
        for (int i = 0; i < nat; ++i)
            for (int j = 0; j <= (want_vir ? i : i - 1); ++j) {   // self images: virial only
                const double rc = rcov[i] + rcov[j];
                const bool same = (i == j);
                PBCUtils::forEachImage(images, &xyz[3 * i], &xyz[3 * j], same,
                    [&](double dx, double dy, double dz, double r2) {
                        const double r = std::sqrt(r2);
                        double dcndr;
                        if (m_method == MethodType::GFN1) {
                            constexpr double k = 16.0;
                            const double ef = std::exp(-k * (rc / r - 1.0));
                            const double c = 1.0 / (1.0 + ef);
                            dcndr = -c * (1.0 - c) * k * rc / r2;
                        } else {
                            constexpr double ka = 10.0, kb = 20.0, rs = 2.0;
                            const double ef1 = std::exp(-ka * (rc / r - 1.0));
                            const double ef2 = std::exp(-kb * ((rc + rs) / r - 1.0));
                            const double c1 = 1.0 / (1.0 + ef1), c2 = 1.0 / (1.0 + ef2);
                            const double dc1dr = -c1 * (1.0 - c1) * ka * rc / r2;
                            const double dc2dr = -c2 * (1.0 - c2) * kb * (rc + rs) / r2;
                            dcndr = dc1dr * c2 + c1 * dc2dr;
                        }
                        const double factor = (same ? dEdcn(i) : dEdcn(i) + dEdcn(j)) * dcndr / r;
                        if (!same) {
                            m_gradient(i, 0) += factor * dx; m_gradient(j, 0) -= factor * dx;
                            m_gradient(i, 1) += factor * dy; m_gradient(j, 1) -= factor * dy;
                            m_gradient(i, 2) += factor * dz; m_gradient(j, 2) -= factor * dz;
                        }
                        if (want_vir) virPair(Xi, factor * dx, factor * dy, factor * dz, dx, dy, dz);
                    });
            }
    }
    if (want_vir) {
        m_virial = 0.5 * (Xi + Xi.transpose());
        m_virial_valid = true;
    }
    if (CurcumaLogger::get_verbosity() >= 3) {
        auto ms = [](gclk::time_point a, gclk::time_point b) { return std::chrono::duration<double, std::milli>(b - a).count(); };
        CurcumaLogger::info_fmt("  periodic gradient: W {:.1f}, repulsion+XB {:.1f}, H0/Pulay {:.1f}, Ewald {:.1f}, "
                                "KO remainder {:.1f}, multipole {:.1f}, CN chain {:.1f} ms (virial {})",
                                ms(tg0, tg_w), ms(tg_w, tg_rep), ms(tg_rep, tg_h0), ms(tg_h0, tg_ew), ms(tg_ew, tg_coul),
                                ms(tg_coul, tg_mp), ms(tg_mp, gclk::now()), want_vir ? "yes" : "no");
    }
}

/* ------------------------------------------------------------------ *
 *  Test hooks (P4): multipole energy / virial at given moments        *
 * ------------------------------------------------------------------ */
double XTB::multipoleEnergyProbe(const Eigen::Matrix3d& cell_angstrom, const Matrix& geometry_angstrom,
                                 const Vector& q, const Matrix& dp, const Matrix& qp)
{
    std::string err;
    if (!setPeriodicCell(cell_angstrom, err)) return 0.0;
    m_geometry = geometry_angstrom;
    m_wfn.q_at = q;
    m_wfn.dp_at = dp;
    m_wfn.qp_at = qp;
    buildMultipoleInteractionMatricesPeriodic();   // current m_mp_mrad, alpha/boxes from the cell (or the env overrides)
    return energyMultipole();
}

Eigen::Matrix3d XTB::multipoleVirialProbe(const Vector& q, const Matrix& dp, const Matrix& qp)
{
    m_wfn.q_at = q;
    m_wfn.dp_at = dp;
    m_wfn.qp_at = qp;
    Matrix g = Matrix::Zero(m_atomcount, 3);
    Vector dEdcn = Vector::Zero(m_atomcount);
    Eigen::Matrix3d V = Eigen::Matrix3d::Zero();
    addMultipoleGradientPeriodic(g, dEdcn, &V);
    return 0.5 * (V + V.transpose());
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
    m_want_virial = gradient;        // P4 probe: dispersionVirialNow() holds the D3/D4 part of Xi
    return calcDispersionEnergy(gradient);
}

} // namespace curcuma::xtb
