/*
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 */
/**
 * Long-range pressure correction of the D4 tail (WP-PERIODIC-NONBONDED V-LRC) - Claude Generated (Nov 2026).
 *
 * The model virial (FFVirialComponents::total, FD-gated by pbc_virial) treats the homogeneous tail at a fixed pair set:
 * tr Xi_tail = 3 E_tail. The pressure of the untruncated pair potential needs tr = 6 E6 + 8 E8 (virial theorem with
 * g(r) = 1 beyond r_c); the difference tr Xi_lrc = 3 E6 + 5 E8 is the component tail_lrc (pressure_lrc, default true).
 * Criteria fixed before the first run:
 * (b) water boxes (20 and 60 H2O, L = 44 Bohr; dispersion_cutoff_bohr 20 and the default L/2):
 *     E6 + E8 = E_tail (rel 1e-12); tr tail_lrc = 3 E6 + 5 E8 (rel 1e-12); tail_lrc isotropic (off-diagonal 0);
 *     GFNFF::virial() = total() + tail_lrc (rel 1e-12); E6/E8 independently from E_tail at r_c = 15 and 20 Bohr
 *     (E_tail = a r^-3 + b r^-5) against the workspace split (rel 1e-9); pressure_lrc false: tail_lrc = 0,
 *     virial() == total() bitwise, every model component bitwise equal to the pressure_lrc true run.
 * (c) physical check, uncorrelated atoms: 800 Ne at random positions (minimum distance 5 Bohr) in L = 70 Bohr
 *     (rho = 2.33e-3 Bohr^-3), 5 configurations, r_c = 15, 20, 25, 30 Bohr, fixed Ewald parameters. The paired
 *     difference dP(r_c) = P(r_c) - P(30) per configuration, averaged; criterion
 *     |dP_lrc(20)| <= 0.2 |dP_model(20)| (variation between 20 and 30 Bohr with and without the correction).
 *     Charges: Ne is a closed-shell noble gas, every atom a neutral fragment; max |q| and E_coul are reported.
 * Exit 0 when every check passes.
 */
#include "src/core/curcuma_logger.h"
#include "src/core/energy_calculators/ff_methods/gfnff.h"
#include "src/core/molecule.h"
#include "src/core/units.h"
#include "test_molecule_registry.h"
#include "json.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

using json = nlohmann::json;

static int g_failures = 0;
static void check(bool ok, const std::string& what, double value, double limit)
{
    std::printf("  [%s] %-84s %.3e (limit %.1e)\n", ok ? "PASS" : "FAIL", what.c_str(), value, limit);
    if (!ok)
        ++g_failures;
}

static const double A2B = CurcumaUnit::Length::ANGSTROM_TO_BOHR;

static Mol withCell(Mol m, double l_ang)
{
    m.m_has_pbc = true;
    m.m_unit_cell = Eigen::Matrix3d::Identity() * l_ang;
    return m;
}

/// Rigid random packing of waters into a cube (Angstrom), every intermolecular pair >= min_any (as test_pbc_virial)
static Mol waterBox(int nwat, double L_bohr, unsigned seed, double min_any = 2.6)
{
    const Mol dimer = TestMolecules::TestMoleculeRegistry::createMolecule("H2O_dimer", false).getMolInfo();
    const double L = L_bohr / A2B;
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    Mol m = dimer;
    m.m_atoms.clear();
    m.m_geometry = Geometry(0, 3);
    m.m_charge = 0;
    const Geometry unit = dimer.m_geometry.topRows(3);
    const Eigen::Vector3d centre = unit.colwise().mean().transpose();
    auto mic = [&](Eigen::Vector3d v) {
        for (int k = 0; k < 3; ++k)
            v[k] -= L * std::rint(v[k] / L);
        return v;
    };
    for (int w = 0; w < nwat; ++w) {
        for (long tries = 0; tries < 200000; ++tries) {
            Eigen::Vector3d p(u(rng) * L, u(rng) * L, u(rng) * L);
            Eigen::Quaterniond qr(u(rng) - 0.5, u(rng) - 0.5, u(rng) - 0.5, u(rng) - 0.5);
            qr.normalize();
            Geometry x(3, 3);
            for (int i = 0; i < 3; ++i)
                x.row(i) = (p + qr.toRotationMatrix() * (unit.row(i).transpose() - centre)).transpose();
            bool ok = true;
            for (int a = 0; a < 3 && ok; ++a)
                for (int b = 0; b < m.m_geometry.rows(); ++b)
                    if (mic(x.row(a).transpose() - m.m_geometry.row(b).transpose()).norm() < min_any) {
                        ok = false;
                        break;
                    }
            if (!ok)
                continue;
            const int n0 = static_cast<int>(m.m_geometry.rows());
            m.m_geometry.conservativeResize(n0 + 3, 3);
            m.m_geometry.bottomRows(3) = x;
            for (int i = 0; i < 3; ++i)
                m.m_atoms.push_back(dimer.m_atoms[i]);
            break;
        }
    }
    m.m_number_atoms = static_cast<int>(m.m_geometry.rows());
    return withCell(m, L);
}

/// n noble-gas atoms (element z) at random positions in a cube of L_bohr, minimum-image distance >= dmin_bohr.
/// Generated, not a stored structure (the registry rule concerns molecular geometries).
static Mol randomAtoms(int z, int n, double L_bohr, double dmin_bohr, unsigned seed)
{
    const double L = L_bohr / A2B, dmin = dmin_bohr / A2B;
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> u(0.0, L);
    Mol m;
    m.m_geometry = Geometry(n, 3);
    int placed = 0;
    while (placed < n) {
        const Eigen::RowVector3d p(u(rng), u(rng), u(rng));
        bool ok = true;
        for (int b = 0; b < placed && ok; ++b) {
            Eigen::RowVector3d d = p - m.m_geometry.row(b);
            for (int k = 0; k < 3; ++k)
                d[k] -= L * std::rint(d[k] / L);
            ok = d.norm() >= dmin;
        }
        if (ok)
            m.m_geometry.row(placed++) = p;
    }
    m.m_atoms.assign(n, z);
    m.m_number_atoms = n;
    m.m_charge = 0;
    m.m_spin = 0;
    return withCell(m, L);
}

static json config(double rc, bool lrc, int threads)
{
    json c = json::object();
    c["verbosity"] = 0;
    c["threads"] = threads;
    c["cache_topology"] = false;
    c["eeq_ewald_tol"] = 1e-12;
    c["ewald_alpha"] = 0.32;    // fixed Ewald parameters: only the D4 radius changes between the runs
    c["ewald_r_cut"] = 19.0;
    c["frag_charge_model"] = "reference";
    if (rc > 0.0)
        c["dispersion_cutoff_bohr"] = rc;
    c["pressure_lrc"] = lrc;
    return c;
}

struct Run {
    bool ok = false;
    double e_tail = 0, e6 = 0, e8 = 0, e_coul = 0, qmax = 0, V = 0;
    FFVirialComponents comp;
    Eigen::Matrix3d virial = Eigen::Matrix3d::Zero();
};

static Run evaluate(const Mol& m, const json& cfg)
{
    Run r;
    try {
        GFNFF ff(cfg);
        ff.setWantVirial(true);
        if (!ff.InitialiseMolecule(m))
            return r;
        ff.Calculation(true);
        if (!ff.virialValid() || !ff.virialComponents() || !ff.getWorkspace())
            return r;
        r.e_tail = ff.DispersionTailEnergy();
        r.e6 = ff.getWorkspace()->dispersionTailE6();
        r.e8 = ff.getWorkspace()->dispersionTailE8();
        r.e_coul = ff.CoulombEnergy();
        const Vector q = ff.Charges();
        r.qmax = q.size() ? q.cwiseAbs().maxCoeff() : 0.0;
        r.comp = *ff.virialComponents();
        r.virial = ff.virial();
        r.V = ff.cellVolumeBohr3();
        r.ok = true;
    } catch (const std::exception& e) {
        std::printf("  evaluate threw: %s\n", e.what());
    }
    return r;
}

static double rel(double a, double b) { return std::abs(a - b) / std::max(std::abs(b), 1e-300); }

static bool bitwiseEqual(const FFVirialComponents& a, const FFVirialComponents& b)
{
    const Eigen::Matrix3d* pa[] = { &a.bond, &a.angle, &a.dihedral, &a.inversion, &a.stors, &a.dispersion, &a.disp_tail,
        &a.repulsion, &a.coulomb, &a.coulomb_bg, &a.hbond, &a.xbond, &a.batm, &a.atm, &a.cn_bond, &a.cn_dispersion, &a.cn_coulomb };
    const Eigen::Matrix3d* pb[] = { &b.bond, &b.angle, &b.dihedral, &b.inversion, &b.stors, &b.dispersion, &b.disp_tail,
        &b.repulsion, &b.coulomb, &b.coulomb_bg, &b.hbond, &b.xbond, &b.batm, &b.atm, &b.cn_bond, &b.cn_dispersion, &b.cn_coulomb };
    for (int k = 0; k < 17; ++k)
        if (!(*pa[k] == *pb[k]))
            return false;
    return true;
}

// ---------------------------------------------------------------------------------------------
// (b) the component on water boxes
static void testComponent(const char* label, const Mol& m, double rc)
{
    std::printf("(b) %s, r_c %s\n", label, rc > 0 ? std::to_string(rc).c_str() : "default (L/2)");
    const Run on = evaluate(m, config(rc, true, 1));
    const Run off = evaluate(m, config(rc, false, 1));
    check(on.ok && off.ok, "virial computed (pressure_lrc true and false)", on.ok && off.ok ? 1 : 0, 1);
    if (!on.ok || !off.ok)
        return;
    const double tr = on.comp.tail_lrc.trace();
    std::printf("  E_tail %.12e  E6 %.12e  E8 %.12e Eh (C8 share %.3f %%), tr Xi_lrc %.12e Eh\n", on.e_tail, on.e6, on.e8,
        100.0 * on.e8 / on.e_tail, tr);
    check(rel(on.e6 + on.e8, on.e_tail) <= 1e-12, "E6 + E8 = E_tail (workspace split)", rel(on.e6 + on.e8, on.e_tail), 1e-12);
    check(rel(tr, 3.0 * on.e6 + 5.0 * on.e8) <= 1e-12, "tr tail_lrc = 3 E6 + 5 E8", rel(tr, 3.0 * on.e6 + 5.0 * on.e8), 1e-12);
    Eigen::Matrix3d offd = on.comp.tail_lrc;
    offd.diagonal().setZero();
    check(offd.cwiseAbs().maxCoeff() == 0.0, "tail_lrc isotropic (off-diagonal exactly 0)", offd.cwiseAbs().maxCoeff(), 0.0);
    const double dv = (on.virial - on.comp.total() - on.comp.tail_lrc).cwiseAbs().maxCoeff() / on.virial.cwiseAbs().maxCoeff();
    check(dv <= 1e-12, "GFNFF::virial() = total() + tail_lrc", dv, 1e-12);
    check(off.comp.tail_lrc.cwiseAbs().maxCoeff() == 0.0, "pressure_lrc false: tail_lrc = 0", off.comp.tail_lrc.cwiseAbs().maxCoeff(), 0.0);
    check(off.virial == off.comp.total(), "pressure_lrc false: virial() == total() bitwise", (off.virial - off.comp.total()).cwiseAbs().maxCoeff(), 0.0);
    check(bitwiseEqual(on.comp, off.comp), "model components bitwise equal for pressure_lrc true / false", bitwiseEqual(on.comp, off.comp) ? 0 : 1, 0);
    std::printf("  pressure of the correction tr/(3V) = %.3f bar (model tail part tr Xi_tail/(3V) = %.3f bar)\n",
        CurcumaUnit::Pressure::au_to_bar(tr / (3.0 * on.V)), CurcumaUnit::Pressure::au_to_bar(on.comp.disp_tail.trace() / (3.0 * on.V)));
}

/// Independent split: E_tail(r) = a r^-3 + b r^-5 from two radii (same CN state, S6 and S8 do not depend on r_c)
static void testIndependentSplit(const Mol& m)
{
    std::printf("(b) independent C6/C8 split from E_tail at r_c = 15 and 20 Bohr (20 H2O)\n");
    const double r1 = 15.0, r2 = 20.0;
    const Run a = evaluate(m, config(r1, true, 1));
    const Run b = evaluate(m, config(r2, true, 1));
    if (!a.ok || !b.ok) {
        check(false, "instances", 0, 1);
        return;
    }
    // [r1^-3 r1^-5; r2^-3 r2^-5] [A; B] = [E1; E2]
    const double m11 = std::pow(r1, -3), m12 = std::pow(r1, -5), m21 = std::pow(r2, -3), m22 = std::pow(r2, -5);
    const double det = m11 * m22 - m12 * m21;
    const double A = (a.e_tail * m22 - m12 * b.e_tail) / det;
    const double B = (m11 * b.e_tail - m21 * a.e_tail) / det;
    std::printf("  E6(20): fit %.12e workspace %.12e; E8(20): fit %.12e workspace %.12e\n", A * m21, b.e6, B * m22, b.e8);
    check(rel(A * m21, b.e6) <= 1e-9, "E6 from the two-radius fit = workspace E6", rel(A * m21, b.e6), 1e-9);
    check(rel(B * m22, b.e8) <= 1e-9, "E8 from the two-radius fit = workspace E8", rel(B * m22, b.e8), 1e-9);
}

// ---------------------------------------------------------------------------------------------
// (c) r_c independence of the corrected pressure for uncorrelated atoms
static void testUncorrelatedGas()
{
    const int N = 800, ncfg = 5;
    const double L = 70.0;
    const double radii[4] = { 15.0, 20.0, 25.0, 30.0 };
    std::printf("(c) %d Ne at random positions (d_min 5 Bohr), L %.1f Bohr, rho %.4e Bohr^-3, %d configurations\n", N, L,
        N / (L * L * L), ncfg);
    double dP_model[4] = {}, dP_lrc[4] = {}, dP_model2[4] = {}, dP_lrc2[4] = {};
    double qmax = 0.0, ecoul = 0.0;
    for (int c = 0; c < ncfg; ++c) {
        const Mol m = randomAtoms(10, N, L, 5.0, 20261109u + 17u * c);
        double Pm[4], Pl[4];
        for (int k = 0; k < 4; ++k) {
            const Run r = evaluate(m, config(radii[k], true, 4));
            if (!r.ok) {
                check(false, "Ne instance", 0, 1);
                return;
            }
            Pm[k] = CurcumaUnit::Pressure::au_to_bar(r.comp.total().trace() / (3.0 * r.V));
            Pl[k] = CurcumaUnit::Pressure::au_to_bar(r.virial.trace() / (3.0 * r.V));
            qmax = std::max(qmax, r.qmax);
            ecoul = std::max(ecoul, std::abs(r.e_coul));
        }
        std::printf("  cfg %d: P_model(15/20/25/30) %.4f %.4f %.4f %.4f bar; P_lrc %.4f %.4f %.4f %.4f bar\n", c, Pm[0], Pm[1],
            Pm[2], Pm[3], Pl[0], Pl[1], Pl[2], Pl[3]);
        for (int k = 0; k < 4; ++k) {
            const double dm = Pm[k] - Pm[3], dl = Pl[k] - Pl[3];
            dP_model[k] += dm / ncfg;
            dP_lrc[k] += dl / ncfg;
            dP_model2[k] += dm * dm;
            dP_lrc2[k] += dl * dl;
        }
    }
    std::printf("  max |q| %.3e e, max |E_coul| %.3e Eh (noble gas: Coulomb and EEQ without effect)\n", qmax, ecoul);
    for (int k = 0; k < 3; ++k) {
        const double sm = std::sqrt(std::max(0.0, dP_model2[k] / ncfg - dP_model[k] * dP_model[k]) / (ncfg - 1));
        const double sl = std::sqrt(std::max(0.0, dP_lrc2[k] / ncfg - dP_lrc[k] * dP_lrc[k]) / (ncfg - 1));
        std::printf("  dP(r_c = %.0f) - dP(30): model %+.4f +- %.4f bar, with lrc %+.4f +- %.4f bar (mean +- SE)\n", radii[k],
            dP_model[k], sm, dP_lrc[k], sl);
    }
    check(qmax < 1e-8, "Ne charges vanish (max |q|)", qmax, 1e-8);
    check(std::abs(dP_model[1]) > 0.1, "model pressure depends on r_c (|dP_model(20)| in bar, sanity)", std::abs(dP_model[1]), 0.1);
    const double ratio = std::abs(dP_lrc[1]) / std::abs(dP_model[1]);
    check(ratio <= 0.2, "|dP_lrc(20 vs 30)| <= 0.2 |dP_model(20 vs 30)|", ratio, 0.2);
}

int main()
{
    CurcumaLogger::set_verbosity(0);
    const Mol box20 = waterBox(20, 44.0, 20261105);
    const Mol box60 = waterBox(60, 44.0, 20261108);
    testComponent("20 H2O, L 44 Bohr", box20, 20.0);
    testComponent("20 H2O, L 44 Bohr", box20, 0.0);
    testComponent("60 H2O, L 44 Bohr", box60, 20.0);
    testIndependentSplit(box20);
    testUncorrelatedGas();
    std::printf("%s: %d failure(s)\n", g_failures == 0 ? "PASS" : "FAIL", g_failures);
    return g_failures == 0 ? 0 : 1;
}
