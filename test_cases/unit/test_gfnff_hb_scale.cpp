/*
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 */
/**
 * Diagnostic switch gfnff.hb_scale (factor s on the GFN-FF hydrogen-bond term) - Claude Generated (Nov 2026).
 *
 * hb_scale multiplies E_HB (cases 1-4), its gradient and its virial component; detection, the HB-CN modulation of
 * the X-H bonds and XB stay unscaled (FFWorkspace::calcHydrogenBondsScaled). Checks:
 * (a) E_HB(s) = s E_HB(1) (relative 1e-12, per-case split too) for s = 0, 0.5, 1.5, 2 on the water dimer, a
 *     20-water cluster without cell and the 20-water box with cell; every other term and the EEQ charges are
 *     bitwise those of s = 1 (the charges do not depend on the HB term); E_tot(s) - E_tot(1) = (s-1) E_HB(1).
 * (b) central finite difference of the total gradient in Angstrom at s = 1.5, h = 1e-4 and 2e-4 A: quadratic in h
 *     (error ratio 2h/h within 3.5-4.5) and max |g_an - g_Richardson| <= 1e-6 * max(max|g|, 1e-3) (Eh/A). The raw
 *     h = 1e-4 difference carries the h^2 truncation error (2e-8 Eh/A on the clusters), so it is printed, not gated.
 * (c) virial: Xi_HB(s) = s Xi_HB(1) (1e-12) and the other components unchanged; FD gate of the total virial and
 *     of the HB term at s = 1.5 on the box (isotropic strain, eps = 1e-5, relative 1e-6 per term, as pbc_virial).
 * (d) hb_scale 0 against hbond false: equal per term except where the HB-CN modulation of the X-H bonds enters
 *     (bond energy; hbond false drops it with the detection, hb_scale 0 keeps it); the bond difference is printed.
 * Water monomer from the registry ("H2O_dimer", first molecule); the H-bonded dimer, the cluster and the box are built
 * from it by rigid placement (hbDimer, waterPack; not a minimum, not a liquid).
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
static const double kL = 44.0;  // Bohr

static Mol registryMol(const char* name)
{
    return TestMolecules::TestMoleculeRegistry::createMolecule(name, false).getMolInfo();
}

static Mol waterUnit()
{
    const Mol dimer = registryMol("H2O_dimer");
    Mol w = dimer;
    w.m_geometry = dimer.m_geometry.topRows(3);
    w.m_atoms.assign(dimer.m_atoms.begin(), dimer.m_atoms.begin() + 3);
    w.m_number_atoms = 3;
    return w;
}

/// Hydrogen-bonded water dimer built from two registry monomers (the registry "H2O_dimer" = water-dimer.roo2.8 faces
/// O to O with both H pairs pointing away: no H bond): the acceptor O at 2.9 A from the donor O along the donor O-H,
/// its H-H bisector pointing away from the donor (lone pairs towards the donor H). A start structure, not a minimum.
static Mol hbDimer()
{
    const Mol w = waterUnit();
    const Eigen::Vector3d o = w.m_geometry.row(0).transpose();
    const Eigen::Vector3d u = (w.m_geometry.row(1).transpose() - o).normalized();
    const Eigen::Vector3d b = 0.5 * (w.m_geometry.row(1) + w.m_geometry.row(2)).transpose() - o;
    const Eigen::Matrix3d R = Eigen::Quaterniond::FromTwoVectors(b, u).toRotationMatrix();
    Mol d = w;
    d.m_geometry.conservativeResize(6, 3);
    for (int i = 0; i < 3; ++i) {
        d.m_geometry.row(3 + i) = (o + 2.9 * u + R * (w.m_geometry.row(i).transpose() - o)).transpose();
        d.m_atoms.push_back(w.m_atoms[i]);
    }
    d.m_number_atoms = 6;
    d.m_charge = 0;
    return d;
}

/// Rigid random packing of n waters into a cube of L (Angstrom), intermolecular minimum-image distances >= min_any
static Mol waterPack(int n, double L_ang, unsigned seed, double min_any, bool cell)
{
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    const Mol unit = waterUnit();
    Mol m = unit;
    m.m_atoms.clear();
    m.m_geometry = Geometry(0, 3);
    auto mic = [&](Eigen::Vector3d v) {
        if (cell)
            for (int k = 0; k < 3; ++k)
                v[k] -= L_ang * std::rint(v[k] / L_ang);
        return v;
    };
    const Eigen::Vector3d centre = unit.m_geometry.colwise().mean().transpose();
    for (int w = 0; w < n; ++w) {
        for (long tries = 0; tries < 200000; ++tries) {
            Eigen::Vector3d p(u(rng) * L_ang, u(rng) * L_ang, u(rng) * L_ang);
            Eigen::Quaterniond qr(u(rng) - 0.5, u(rng) - 0.5, u(rng) - 0.5, u(rng) - 0.5);
            qr.normalize();
            const Eigen::Matrix3d R = qr.toRotationMatrix();
            Geometry x(3, 3);
            for (int i = 0; i < 3; ++i)
                x.row(i) = (p + R * (unit.m_geometry.row(i).transpose() - centre)).transpose();
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
                m.m_atoms.push_back(unit.m_atoms[i]);
            break;
        }
    }
    m.m_number_atoms = static_cast<int>(m.m_geometry.rows());
    m.m_charge = 0;
    m.m_has_pbc = cell;
    if (cell) {
        Eigen::Matrix3d c = Eigen::Matrix3d::Zero();
        c(0, 0) = c(1, 1) = c(2, 2) = L_ang;
        m.m_unit_cell = c;
    }
    return m;
}

static json config(double s, bool cell)
{
    json c = json::object();
    c["verbosity"] = 0;
    c["threads"] = 1;
    c["cache_topology"] = false;
    c["hb_update_force_every"] = 1;    // lists from the geometry at hand
    c["hb_min_pair_energy_eh"] = 0.0;  // no energy-based candidate filter (a discrete decision under FD)
    c["hb_scale"] = s;
    if (cell) {  // as pbc_virial: fixed Ewald parameters and an absolute D4 cutoff under the strain
        c["eeq_ewald_tol"] = 1e-14;
        c["ewald_tolerance"] = 1e-14;
        c["ewald_alpha"] = 0.32;
        c["ewald_r_cut"] = 19.0;
        c["dispersion_cutoff_bohr"] = 20.0;
    }
    return c;
}

enum Term { T_BOND, T_ANGLE, T_DIHEDRAL, T_INVERSION, T_STORS, T_DISP, T_REP, T_COULOMB, T_HBOND, T_XBOND, T_BATM, T_TOTAL, NT };
static const char* kTermName[NT] = { "bond", "angle", "dihedral", "inversion", "stors", "dispersion", "repulsion", "coulomb",
    "hbond", "xbond", "batm", "TOTAL" };

struct Eval {
    bool ok = false;
    double e[NT] = {};
    double hb_case[4] = {};
    Eigen::Matrix3d xi[NT];
    Matrix gradient;  // Eh/Bohr
    Vector charges;
};

static Eval evaluate(const Mol& m, const json& cfg, bool gradient, bool virial)
{
    Eval r;
    for (auto& x : r.xi)
        x.setZero();
    try {
        GFNFF ff(cfg);
        if (virial)
            ff.setWantVirial(true);
        if (!ff.InitialiseMolecule(m))
            return r;
        const double e_tot = ff.Calculation(gradient || virial);
        r.e[T_BOND] = ff.BondEnergy();
        r.e[T_ANGLE] = ff.AngleEnergy();
        r.e[T_DIHEDRAL] = ff.DihedralEnergy();
        r.e[T_INVERSION] = ff.InversionEnergy();
        r.e[T_STORS] = ff.STorsEnergy();
        r.e[T_DISP] = ff.DispersionEnergy() + ff.DispersionTailEnergy();
        r.e[T_REP] = ff.BondedRepulsionEnergy() + ff.NonbondedRepulsionEnergy();
        r.e[T_COULOMB] = ff.CoulombEnergy();
        r.e[T_HBOND] = ff.HydrogenBondEnergy();
        r.e[T_XBOND] = ff.HalogenBondEnergy();
        r.e[T_BATM] = ff.BatmEnergy();
        r.e[T_TOTAL] = e_tot;
        if (ff.getWorkspace()) {
            const FFEnergyComponents& ec = ff.getWorkspace()->energyComponents();
            r.hb_case[0] = ec.hbond_case1;
            r.hb_case[1] = ec.hbond_case2;
            r.hb_case[2] = ec.hbond_case3;
            r.hb_case[3] = ec.hbond_case4;
            if (gradient || virial)
                r.gradient = ff.getWorkspace()->gradient();
        }
        r.charges = ff.Charges();
        if (virial && ff.virialComponents()) {
            const FFVirialComponents& v = *ff.virialComponents();
            r.xi[T_BOND] = v.bondTotal();
            r.xi[T_ANGLE] = v.angle;
            r.xi[T_DIHEDRAL] = v.dihedral;
            r.xi[T_INVERSION] = v.inversion;
            r.xi[T_STORS] = v.stors;
            r.xi[T_DISP] = v.dispersionTotal();
            r.xi[T_REP] = v.repulsion;
            r.xi[T_COULOMB] = v.coulombTotal();
            r.xi[T_HBOND] = v.hbond;
            r.xi[T_XBOND] = v.xbond;
            r.xi[T_BATM] = v.batm;
            r.xi[T_TOTAL] = v.total();
        }
        r.ok = true;
    } catch (const std::exception& e) {
        std::printf("  evaluate threw: %s\n", e.what());
    }
    return r;
}

static double rel(double a, double b) { return std::abs(a - b) / std::max(std::abs(b), 1e-300); }

// (a) energy linearity, other terms and charges unchanged
static void testEnergy(const char* label, const Mol& m)
{
    const bool cell = m.m_has_pbc;
    std::printf("(a) energy: %s (%d atoms%s)\n", label, m.m_number_atoms, cell ? ", cell" : "");
    const Eval e1 = evaluate(m, config(1.0, cell), false, false);
    check(e1.ok && e1.e[T_HBOND] < 0.0, std::string(label) + ": E_HB(1) < 0", e1.e[T_HBOND], 0.0);
    std::printf("  E_HB(1) = %.15e Eh, cases %.6e %.6e %.6e %.6e\n", e1.e[T_HBOND], e1.hb_case[0], e1.hb_case[1], e1.hb_case[2], e1.hb_case[3]);
    for (double s : { 0.0, 0.5, 1.5, 2.0 }) {
        const Eval es = evaluate(m, config(s, cell), false, false);
        if (!es.ok) {
            check(false, "instance at s", s, 0);
            continue;
        }
        char buf[160];
        const double d_hb = std::abs(es.e[T_HBOND] - s * e1.e[T_HBOND]) / std::abs(e1.e[T_HBOND]);
        std::snprintf(buf, sizeof buf, "%s s=%.2f: |E_HB(s) - s E_HB(1)| / |E_HB(1)|", label, s);
        check(d_hb <= 1e-12, buf, d_hb, 1e-12);
        double d_case = 0.0;
        for (int c = 0; c < 4; ++c)
            d_case = std::max(d_case, std::abs(es.hb_case[c] - s * e1.hb_case[c]) / std::abs(e1.e[T_HBOND]));
        std::snprintf(buf, sizeof buf, "%s s=%.2f: per-case split scales (max over cases, rel. to |E_HB(1)|)", label, s);
        check(d_case <= 1e-12, buf, d_case, 1e-12);
        double d_other = 0.0;
        for (int t = 0; t < T_TOTAL; ++t)
            if (t != T_HBOND)
                d_other = std::max(d_other, std::abs(es.e[t] - e1.e[t]));
        std::snprintf(buf, sizeof buf, "%s s=%.2f: every other term bitwise unchanged (max abs diff, Eh)", label, s);
        check(d_other == 0.0, buf, d_other, 0.0);
        const double d_q = (es.charges.size() == e1.charges.size()) ? (es.charges - e1.charges).cwiseAbs().maxCoeff() : 1.0;
        std::snprintf(buf, sizeof buf, "%s s=%.2f: EEQ charges bitwise unchanged (max abs diff, e)", label, s);
        check(d_q == 0.0, buf, d_q, 0.0);
        const double d_tot = std::abs((es.e[T_TOTAL] - e1.e[T_TOTAL]) - (s - 1.0) * e1.e[T_HBOND]) / std::abs(e1.e[T_TOTAL]);
        std::snprintf(buf, sizeof buf, "%s s=%.2f: E_tot(s) - E_tot(1) = (s-1) E_HB(1) (rel. to |E_tot|)", label, s);
        check(d_tot <= 1e-13, buf, d_tot, 1e-13);
    }
}

// (b) central FD of the total gradient (Angstrom) at s = 1.5 on the first n_atoms atoms
static void testGradientFD(const char* label, const Mol& m, int n_atoms)
{
    const bool cell = m.m_has_pbc;
    const json cfg = config(1.5, cell);
    std::printf("(b) FD gradient at s = 1.5: %s (%d of %d atoms)\n", label, n_atoms, m.m_number_atoms);
    const Eval an = evaluate(m, cfg, true, false);
    if (!an.ok) {
        check(false, "analytic instance", 0, 0);
        return;
    }
    double gmax = 0.0, err_h = 0.0, err_2h = 0.0, err_rich = 0.0;
    const double h = 1e-4;  // Angstrom
    for (int i = 0; i < n_atoms; ++i)
        for (int k = 0; k < 3; ++k) {
            double e[4];
            const double steps[4] = { -2 * h, -h, h, 2 * h };
            for (int q = 0; q < 4; ++q) {
                Mol d = m;
                d.m_geometry(i, k) += steps[q];
                e[q] = evaluate(d, cfg, false, false).e[T_TOTAL];
            }
            const double fd_h = (e[2] - e[1]) / (2 * h);
            const double fd_2h = (e[3] - e[0]) / (4 * h);
            const double rich = (4 * fd_h - fd_2h) / 3;
            const double g = an.gradient(i, k) * A2B;  // Eh/Bohr -> Eh/A
            gmax = std::max(gmax, std::abs(fd_h));
            err_h = std::max(err_h, std::abs(g - fd_h));
            err_2h = std::max(err_2h, std::abs(g - fd_2h));
            err_rich = std::max(err_rich, std::abs(g - rich));
        }
    std::printf("  max|g| %.4e Eh/A; max|g_an - FD|: h %.3e, 2h %.3e, Richardson %.3e (ratio 2h/h %.2f)\n", gmax, err_h, err_2h,
        err_rich, err_2h / std::max(err_h, 1e-300));
    const double lim = 1e-6 * std::max(gmax, 1e-3);
    const double ratio = err_2h / std::max(err_h, 1e-300);
    check(ratio > 3.5 && ratio < 4.5, std::string(label) + ": FD error quadratic in h (ratio 2h/h, limit 3.5-4.5)", ratio, 4.0);
    check(err_rich <= lim, std::string(label) + ": max |g_an - g_Richardson|", err_rich, lim);
}

// (c) virial scaling and FD gate at s = 1.5
static Mol strained(const Mol& m, double eps)
{
    Mol s = m;
    s.m_geometry *= (1.0 + eps);
    s.m_unit_cell *= (1.0 + eps);
    return s;
}

static void testVirial(const char* label, const Mol& m)
{
    std::printf("(c) virial: %s\n", label);
    const Eval v1 = evaluate(m, config(1.0, true), true, true);
    const Eval v15 = evaluate(m, config(1.5, true), true, true);
    if (!v1.ok || !v15.ok) {
        check(false, "virial instances", 0, 0);
        return;
    }
    const double nrm = v1.xi[T_HBOND].cwiseAbs().maxCoeff();
    const double d_hb = (v15.xi[T_HBOND] - 1.5 * v1.xi[T_HBOND]).cwiseAbs().maxCoeff() / nrm;
    check(d_hb <= 1e-12, std::string(label) + ": Xi_HB(1.5) = 1.5 Xi_HB(1) (max abs / max |Xi_HB(1)|)", d_hb, 1e-12);
    double d_other = 0.0;
    for (int t = 0; t < T_TOTAL; ++t)
        if (t != T_HBOND)
            d_other = std::max(d_other, (v15.xi[t] - v1.xi[t]).cwiseAbs().maxCoeff());
    check(d_other == 0.0, std::string(label) + ": every other virial component bitwise unchanged (max abs, Eh)", d_other, 0.0);
    const double eps_list[4] = { -2e-5, -1e-5, 1e-5, 2e-5 };
    Eval ev[4];
    for (int q = 0; q < 4; ++q) {
        ev[q] = evaluate(strained(m, eps_list[q]), config(1.5, true), false, false);
        if (!ev[q].ok) {
            check(false, "strained instance", 0, 0);
            return;
        }
    }
    for (int t : { (int)T_HBOND, (int)T_TOTAL }) {
        const double fd = -(ev[2].e[t] - ev[1].e[t]) / (2e-5);
        const double fd2 = -(ev[3].e[t] - ev[0].e[t]) / (4e-5);
        const double rich = (4 * fd - fd2) / 3;
        const double an = v15.xi[t].trace();
        const double r1 = std::abs(an - fd) / std::max(std::abs(fd), 1e-6);
        std::printf("  %-8s Xi_an %.9e  FD %.9e  Richardson %.9e Eh\n", kTermName[t], an, fd, rich);
        check(r1 <= 1e-6, std::string(label) + ": FD gate s=1.5 " + kTermName[t] + ": |Xi_an - Xi_FD| / max(|Xi_FD|, 1e-6)", r1, 1e-6);
    }
}

// (d) hb_scale 0 against hbond false
static void testOff(const char* label, const Mol& m)
{
    const bool cell = m.m_has_pbc;
    std::printf("(d) hb_scale 0 vs hbond false: %s\n", label);
    const Eval z = evaluate(m, config(0.0, cell), false, false);
    json c = config(1.0, cell);
    c["hbond"] = false;
    const Eval off = evaluate(m, c, false, false);
    if (!z.ok || !off.ok) {
        check(false, "instances", 0, 0);
        return;
    }
    check(z.e[T_HBOND] == 0.0 && off.e[T_HBOND] == 0.0, std::string(label) + ": E_HB = 0 for both", std::abs(z.e[T_HBOND]) + std::abs(off.e[T_HBOND]), 0.0);
    std::printf("  per-term E(hb_scale 0) - E(hbond false), Eh:");
    double d_other = 0.0;
    for (int t = 0; t < T_TOTAL; ++t) {
        std::printf(" %s %.3e", kTermName[t], z.e[t] - off.e[t]);
        if (t != T_BOND && t != T_XBOND)
            d_other = std::max(d_other, std::abs(z.e[t] - off.e[t]));
    }
    std::printf("\n  total difference %.6e Eh (%.4f kcal/mol), bond part %.6e Eh: the HB-CN modulation of the X-H bonds, "
                "kept by hb_scale 0, dropped by hbond false\n",
        z.e[T_TOTAL] - off.e[T_TOTAL], (z.e[T_TOTAL] - off.e[T_TOTAL]) * CurcumaUnit::Energy::HARTREE_TO_KCALMOL, z.e[T_BOND] - off.e[T_BOND]);
    check(d_other <= 1e-12, std::string(label) + ": all terms except bond/xbond equal (max abs, Eh)", d_other, 1e-12);
    const double d_q = (z.charges - off.charges).cwiseAbs().maxCoeff();
    check(d_q == 0.0, std::string(label) + ": EEQ charges equal (max abs, e)", d_q, 0.0);
}

int main()
{
    CurcumaLogger::set_verbosity(0);
    const Mol dimer = hbDimer();
    const Mol cluster = waterPack(20, 9.0, 20261108, 2.4, false);
    const Mol box20 = waterPack(20, kL / A2B, 20261105, 2.6, true);  // as pbc_virial: sparse, few H bonds
    // HB-rich periodic case: the 20-water cluster in the L = 44 Bohr cell, shifted by -2 A so that molecules are torn
    // across the faces (exercises the A4b image frame of the HB kernel)
    Mol boxc = cluster;
    boxc.m_geometry.array() -= 2.0;
    boxc.m_has_pbc = true;
    boxc.m_unit_cell = Eigen::Matrix3d::Identity() * (kL / A2B);
    testEnergy("water dimer", dimer);
    testEnergy("20 H2O cluster", cluster);
    testEnergy("20 H2O box (sparse)", box20);
    testEnergy("20 H2O cluster in cell", boxc);
    testGradientFD("water dimer", dimer, dimer.m_number_atoms);
    testGradientFD("20 H2O cluster", cluster, 12);
    testGradientFD("20 H2O cluster in cell", boxc, 12);
    testVirial("20 H2O box (sparse)", box20);
    testVirial("20 H2O cluster in cell", boxc);
    testOff("water dimer", dimer);
    testOff("20 H2O cluster", cluster);
    testOff("20 H2O cluster in cell", boxc);
    std::printf("%s: %d failure(s)\n", g_failures == 0 ? "PASS" : "FAIL", g_failures);
    return g_failures == 0 ? 0 : 1;
}
