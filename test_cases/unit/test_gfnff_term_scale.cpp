/*
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 */
/**
 * Diagnostic term-scale switches gfnff.repulsion_scale, dispersion_scale, coulomb_scale - Claude Generated (Nov 2026).
 *
 * Each factor s multiplies one GFN-FF energy component with its gradient and virial (FFWorkspace::calcTermScaled,
 * setDispersionTail, postProcess); 1.0 is the reference. Pattern of test_gfnff_hb_scale. Checks, per switch:
 * (a) E_term(s) = s E_term(1) (relative 1e-12) for s = 0, 0.5, 1.5 on the water dimer, a 20-water cluster without
 *     cell and the same cluster in a cell (Ewald Coulomb, D4 tail); every other term and the EEQ charges bitwise those
 *     of s = 1; E_tot(s) - E_tot(1) = (s-1) E_term(1). Dispersion: pair term, tail and its E6/E8 split separately.
 * (b) central finite difference of the TOTAL gradient (Angstrom) at s = 1.5, every displaced geometry a fresh
 *     instance (charges solved anew): quadratic in h (ratio 2h/h within 3.5-4.5) and max |g_an - g_Richardson|
 *     <= 1e-8 Eh/A. For coulomb_scale this is the Hellmann-Feynman check: the analytic gradient has no dq/dx term.
 * (c) virial on the sparse 20-water box and the cluster in the cell: Xi_term(1.5) = 1.5 Xi_term(1) (1e-12), the
 *     other components bitwise; FD gate (isotropic strain, as pbc_virial) of the term and of the model total at
 *     s = 1.5 (relative 1e-6). dispersion_scale with pressure_lrc true and false: tail_lrc scales, is 0 when off,
 *     and the model components are bitwise equal between the two.
 * Structures are built as in test_gfnff_hb_scale (rigid placement of registry monomers; not minima, not a liquid).
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
static int g_checks = 0;
static void check(bool ok, const std::string& what, double value, double limit)
{
    ++g_checks;
    std::printf("  [%s] %-92s %.3e (limit %.1e)\n", ok ? "PASS" : "FAIL", what.c_str(), value, limit);
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

/// Hydrogen-bonded water dimer from two registry monomers (as test_gfnff_hb_scale); a start structure, not a minimum.
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

enum Switch { S_REP, S_DISP, S_COUL, NS };
static const char* kSwitchKey[NS] = { "repulsion_scale", "dispersion_scale", "coulomb_scale" };

static json config(int sw, double s, bool cell, bool lrc = true)
{
    json c = json::object();
    c["verbosity"] = 0;
    c["threads"] = 1;
    c["cache_topology"] = false;
    c["hb_update_force_every"] = 1;    // lists from the geometry at hand
    c["hb_min_pair_energy_eh"] = 0.0;  // no energy-based candidate filter (a discrete decision under FD)
    c[kSwitchKey[sw]] = s;
    c["pressure_lrc"] = lrc;
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
static const int kSwitchTerm[NS] = { T_REP, T_DISP, T_COULOMB };

struct Eval {
    bool ok = false;
    double e[NT] = {};
    double disp_pair = 0, disp_tail = 0, e6 = 0, e8 = 0, brep = 0, nbrep = 0;
    Eigen::Matrix3d xi[NT];
    Eigen::Matrix3d tail_lrc = Eigen::Matrix3d::Zero();
    FFVirialComponents comp;
    Matrix gradient;  // Eh/Bohr
    Vector charges;
};

static Eval evaluate(const Mol& m, const json& cfg, bool gradient, bool virial)
{
    Eval r;
    for (auto& x : r.xi)
        x.setZero();
    r.comp.reset();
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
        r.disp_pair = ff.DispersionEnergy();
        r.disp_tail = ff.DispersionTailEnergy();
        r.e[T_DISP] = r.disp_pair + r.disp_tail;
        r.brep = ff.BondedRepulsionEnergy();
        r.nbrep = ff.NonbondedRepulsionEnergy();
        r.e[T_REP] = r.brep + r.nbrep;
        r.e[T_COULOMB] = ff.CoulombEnergy();
        r.e[T_HBOND] = ff.HydrogenBondEnergy();
        r.e[T_XBOND] = ff.HalogenBondEnergy();
        r.e[T_BATM] = ff.BatmEnergy();
        r.e[T_TOTAL] = e_tot;
        if (ff.getWorkspace()) {
            r.e6 = ff.getWorkspace()->dispersionTailE6();
            r.e8 = ff.getWorkspace()->dispersionTailE8();
            if (gradient || virial)
                r.gradient = ff.getWorkspace()->gradient();
        }
        r.charges = ff.Charges();
        if (virial && ff.virialComponents()) {
            const FFVirialComponents& v = *ff.virialComponents();
            r.comp = v;
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
            r.tail_lrc = v.tail_lrc;
        }
        r.ok = true;
    } catch (const std::exception& e) {
        std::printf("  evaluate threw: %s\n", e.what());
    }
    return r;
}

static double relTo(double a, double b) { return std::abs(a) / std::max(std::abs(b), 1e-300); }

static bool bitwiseEqual(const FFVirialComponents& a, const FFVirialComponents& b)
{
    return a.bond == b.bond && a.angle == b.angle && a.dihedral == b.dihedral && a.inversion == b.inversion && a.stors == b.stors
        && a.dispersion == b.dispersion && a.disp_tail == b.disp_tail && a.repulsion == b.repulsion && a.coulomb == b.coulomb
        && a.coulomb_bg == b.coulomb_bg && a.hbond == b.hbond && a.xbond == b.xbond && a.batm == b.batm && a.atm == b.atm
        && a.cn_bond == b.cn_bond && a.cn_dispersion == b.cn_dispersion && a.cn_coulomb == b.cn_coulomb;
}

// (a) energy linearity, other terms and charges unchanged
static void testEnergy(int sw, const char* label, const Mol& m)
{
    const bool cell = m.m_has_pbc;
    const int T = kSwitchTerm[sw];
    std::printf("(a) %s energy: %s (%d atoms%s)\n", kSwitchKey[sw], label, m.m_number_atoms, cell ? ", cell" : "");
    const Eval e1 = evaluate(m, config(sw, 1.0, cell), false, false);
    check(e1.ok && e1.e[T] != 0.0, std::string(label) + ": E_term(1) != 0", e1.e[T], 0.0);
    std::printf("  E_%s(1) = %.15e Eh (disp pair %.6e tail %.6e E6 %.6e E8 %.6e; rep bonded %.6e nonbonded %.6e)\n",
        kTermName[T], e1.e[T], e1.disp_pair, e1.disp_tail, e1.e6, e1.e8, e1.brep, e1.nbrep);
    for (double s : { 0.0, 0.5, 1.5 }) {
        const Eval es = evaluate(m, config(sw, s, cell), false, false);
        if (!es.ok) {
            check(false, "instance at s", s, 0);
            continue;
        }
        char buf[200];
        const std::string pre = std::string(label) + " s=" + std::to_string(s).substr(0, 4) + ": ";
        const double d = relTo(es.e[T] - s * e1.e[T], e1.e[T]);
        check(d <= 1e-12, pre + "|E_term(s) - s E_term(1)| / |E_term(1)|", d, 1e-12);
        if (sw == S_DISP) {
            double d_parts = std::max(relTo(es.disp_pair - s * e1.disp_pair, e1.e[T]), relTo(es.disp_tail - s * e1.disp_tail, e1.e[T]));
            d_parts = std::max(d_parts, std::max(relTo(es.e6 - s * e1.e6, e1.e[T]), relTo(es.e8 - s * e1.e8, e1.e[T])));
            check(d_parts <= 1e-12, pre + "pair, tail, E6, E8 each scale (rel. to |E_disp(1)|)", d_parts, 1e-12);
        }
        if (sw == S_REP) {
            const double d_parts = std::max(relTo(es.brep - s * e1.brep, e1.e[T]), relTo(es.nbrep - s * e1.nbrep, e1.e[T]));
            check(d_parts <= 1e-12, pre + "bonded and non-bonded part each scale (rel. to |E_rep(1)|)", d_parts, 1e-12);
        }
        double d_other = 0.0;
        for (int t = 0; t < T_TOTAL; ++t)
            if (t != T)
                d_other = std::max(d_other, std::abs(es.e[t] - e1.e[t]));
        check(d_other == 0.0, pre + "every other term bitwise unchanged (max abs diff, Eh)", d_other, 0.0);
        const double d_q = (es.charges.size() == e1.charges.size()) ? (es.charges - e1.charges).cwiseAbs().maxCoeff() : 1.0;
        check(d_q == 0.0, pre + "EEQ charges bitwise unchanged (max abs diff, e)", d_q, 0.0);
        const double d_tot = relTo((es.e[T_TOTAL] - e1.e[T_TOTAL]) - (s - 1.0) * e1.e[T], e1.e[T_TOTAL]);
        std::snprintf(buf, sizeof buf, "E_tot(s) - E_tot(1) = (s-1) E_term(1) (rel. to |E_tot|)");
        check(d_tot <= 1e-13, pre + buf, d_tot, 1e-13);
    }
}

// (b) central FD of the total gradient (Angstrom) at s = 1.5 on the first n_atoms atoms
static void testGradientFD(int sw, const char* label, const Mol& m, int n_atoms)
{
    const bool cell = m.m_has_pbc;
    const json cfg = config(sw, 1.5, cell);
    std::printf("(b) %s FD gradient at s = 1.5: %s (%d of %d atoms)\n", kSwitchKey[sw], label, n_atoms, m.m_number_atoms);
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
    const double ratio = err_2h / std::max(err_h, 1e-300);
    const std::string pre = std::string(kSwitchKey[sw]) + " " + label + ": ";
    check(ratio > 3.5 && ratio < 4.5, pre + "FD error quadratic in h (ratio 2h/h, limit 3.5-4.5)", ratio, 4.0);
    check(err_rich <= 1e-8, pre + "max |g_an - g_Richardson| (Eh/A, charges solved anew)", err_rich, 1e-8);
}

// (c) virial scaling and FD gate at s = 1.5
static Mol strained(const Mol& m, double eps)
{
    Mol s = m;
    s.m_geometry *= (1.0 + eps);
    s.m_unit_cell *= (1.0 + eps);
    return s;
}

static void testVirial(int sw, const char* label, const Mol& m, bool lrc)
{
    const int T = kSwitchTerm[sw];
    std::printf("(c) %s virial: %s (pressure_lrc %s)\n", kSwitchKey[sw], label, lrc ? "true" : "false");
    const Eval v1 = evaluate(m, config(sw, 1.0, true, lrc), true, true);
    const Eval v15 = evaluate(m, config(sw, 1.5, true, lrc), true, true);
    if (!v1.ok || !v15.ok) {
        check(false, "virial instances", 0, 0);
        return;
    }
    const std::string pre = std::string(kSwitchKey[sw]) + " " + label + (lrc ? " lrc on" : " lrc off") + ": ";
    const double nrm = v1.xi[T].cwiseAbs().maxCoeff();
    const double d = (v15.xi[T] - 1.5 * v1.xi[T]).cwiseAbs().maxCoeff() / nrm;
    check(d <= 1e-12, pre + "Xi_term(1.5) = 1.5 Xi_term(1) (max abs / max |Xi_term(1)|)", d, 1e-12);
    // Other components: equal up to rounding, not bitwise. The kernels that run after the scaled one (HB, XB, sTors)
    // extract their virial from the difference of the accumulated gradient rows (VirialTermRows), whose base value
    // now differs by the scaled term, so the last bits of their virial move (energies and gradients stay bitwise).
    double d_other = 0.0;
    int t_worst = -1;
    for (int t = 0; t < T_TOTAL; ++t)
        if (t != T) {
            const double n = std::max(v1.xi[t].cwiseAbs().maxCoeff(), 1e-300);
            const double dt = (v15.xi[t] - v1.xi[t]).cwiseAbs().maxCoeff() / n;
            if (dt > d_other) {
                d_other = dt;
                t_worst = t;
            }
        }
    check(d_other <= 1e-12, pre + "other virial components unchanged (max rel., worst: " + (t_worst >= 0 ? kTermName[t_worst] : "none") + ")", d_other, 1e-12);
    if (sw == S_DISP) {
        if (lrc) {
            const double nl = v1.tail_lrc.cwiseAbs().maxCoeff();
            const double dl = (v15.tail_lrc - 1.5 * v1.tail_lrc).cwiseAbs().maxCoeff() / std::max(nl, 1e-300);
            check(nl > 0.0 && dl <= 1e-12, pre + "tail_lrc(1.5) = 1.5 tail_lrc(1)", dl, 1e-12);
            const Eval off = evaluate(m, config(sw, 1.5, true, false), true, true);
            check(off.ok && bitwiseEqual(off.comp, v15.comp), pre + "model components bitwise equal to pressure_lrc false",
                off.ok && bitwiseEqual(off.comp, v15.comp) ? 0.0 : 1.0, 0.0);
        } else {
            check(v15.tail_lrc.cwiseAbs().maxCoeff() == 0.0, pre + "tail_lrc = 0", v15.tail_lrc.cwiseAbs().maxCoeff(), 0.0);
        }
    } else {
        const double dl = (v15.tail_lrc - v1.tail_lrc).cwiseAbs().maxCoeff();
        check(dl == 0.0, pre + "tail_lrc bitwise unchanged", dl, 0.0);
    }
    const double eps_list[4] = { -2e-5, -1e-5, 1e-5, 2e-5 };
    Eval ev[4];
    for (int q = 0; q < 4; ++q) {
        ev[q] = evaluate(strained(m, eps_list[q]), config(sw, 1.5, true, lrc), false, false);
        if (!ev[q].ok) {
            check(false, "strained instance", 0, 0);
            return;
        }
    }
    for (int t : { T, (int)T_TOTAL }) {
        const double fd = -(ev[2].e[t] - ev[1].e[t]) / (2e-5);
        const double fd2 = -(ev[3].e[t] - ev[0].e[t]) / (4e-5);
        const double rich = (4 * fd - fd2) / 3;
        const double an = v15.xi[t].trace();  // model virial (total() without tail_lrc)
        const double r1 = std::abs(an - fd) / std::max(std::abs(fd), 1e-6);
        std::printf("  %-10s Xi_an %.9e  FD %.9e  Richardson %.9e Eh\n", kTermName[t], an, fd, rich);
        check(r1 <= 1e-6, pre + "FD gate s=1.5 " + kTermName[t] + ": |Xi_an - Xi_FD| / max(|Xi_FD|, 1e-6)", r1, 1e-6);
    }
}

int main()
{
    CurcumaLogger::set_verbosity(0);
    const Mol dimer = hbDimer();
    const Mol cluster = waterPack(20, 9.0, 20261108, 2.4, false);
    const Mol box20 = waterPack(20, kL / A2B, 20261105, 2.6, true);  // as pbc_virial: sparse
    Mol boxc = cluster;  // the 20-water cluster in the L = 44 Bohr cell, torn across the faces (as test_gfnff_hb_scale)
    boxc.m_geometry.array() -= 2.0;
    boxc.m_has_pbc = true;
    boxc.m_unit_cell = Eigen::Matrix3d::Identity() * (kL / A2B);
    for (int sw = 0; sw < NS; ++sw) {
        testEnergy(sw, "water dimer", dimer);
        testEnergy(sw, "20 H2O cluster", cluster);
        testEnergy(sw, "20 H2O cluster in cell", boxc);
        testGradientFD(sw, "water dimer", dimer, dimer.m_number_atoms);
        testGradientFD(sw, "20 H2O cluster", cluster, 12);
        testGradientFD(sw, "20 H2O cluster in cell", boxc, 12);
        testVirial(sw, "20 H2O box (sparse)", box20, true);
        testVirial(sw, "20 H2O cluster in cell", boxc, true);
        if (sw == S_DISP) {
            testVirial(sw, "20 H2O box (sparse)", box20, false);
            testVirial(sw, "20 H2O cluster in cell", boxc, false);
        }
    }
    std::printf("%s: %d of %d checks failed\n", g_failures == 0 ? "PASS" : "FAIL", g_failures, g_checks);
    return g_failures == 0 ? 0 : 1;
}
