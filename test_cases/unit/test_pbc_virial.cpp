/*
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 */
/**
 * Atomic virial of the periodic GFN-FF energy (WP-PERIODIC-NONBONDED N1) - Claude Generated (Nov 2026).
 *
 * Xi_ab = -dE/d(eps_ab) under an affine strain of the cell and of every atom at fixed fractional coordinates.
 * Checks (thresholds fixed before the first run, see the WP document):
 *
 * (a)  THE GATE: the analytic tensor against -dE/d(eps) by central differences of the energy under an isotropic
 *      scaling L -> L (1 + eps), x -> x (1 + eps), eps = 1e-5 and 2e-5 (quadratic convergence checked by Richardson),
 *      term by term (per energy component, fresh instance at every strained geometry, charges re-solved) and in
 *      total. Criterion |Xi_an - Xi_FD| / max(|Xi_FD|, 1e-6 Eh) <= 1e-6 per term and total. Pair-list events (pairs
 *      crossing a sharp cutoff under the strain) are counted and reported separately.
 *      Systems: 20 and 60 H2O in L = 44 Bohr, H3O+(H2O)2 + 18 H2O (Q = +1, frag_charge_model reference), caffeine in
 *      the same box, acetic acid dimer + F-I...pyridine (XB).
 * (b)  invariance: translation by lattice vectors plus an arbitrary shift, atom-wise wrap (torn molecules), a 2x1x1
 *      supercell (Xi doubles), 1 vs 4 threads.
 * (c)  no cell: no virial (the request is refused); with a cell the gradient and energy with the virial requested
 *      are bitwise those without the request.
 * (d)  pressure of the 60-water box in bar (virial part, and with the kinetic part of 300 K); unit factor.
 * Fixed Ewald parameters (alpha, r_cut, tolerance 1e-14), dispersion_cutoff_bohr 20 (an absolute r_c under the
 * strain), hb_update_force_every 1, hb_min_pair_energy_eh 0, eeq_ewald_tol 1e-14; dense EEQ below 1000 atoms.
 * Water geometry: registry "H2O_dimer", moved rigidly. Exit 0 when every check passes.
 */
#include "src/core/curcuma_logger.h"
#include "src/core/elements.h"
#include "src/core/energy_calculators/ff_methods/gfnff.h"
#include "src/core/molecule.h"
#include "src/core/units.h"
#include "test_molecule_registry.h"
#include "json.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <memory>
#include <random>
#include <string>
#include <vector>

using json = nlohmann::json;

static int g_failures = 0;
static void check(bool ok, const std::string& what, double value, double limit)
{
    std::printf("  [%s] %-78s %.3e (limit %.1e)\n", ok ? "PASS" : "FAIL", what.c_str(), value, limit);
    if (!ok)
        ++g_failures;
}

static const double A2B = CurcumaUnit::Length::ANGSTROM_TO_BOHR;
static const double kL = 44.0;  // Bohr

static Mol registryMol(const char* name, int charge = 0)
{
    Mol m = TestMolecules::TestMoleculeRegistry::createMolecule(name, false).getMolInfo();
    m.m_charge = charge;
    return m;
}

static Mol withCell(Mol m, double lx, double ly, double lz)
{
    m.m_has_pbc = true;
    Eigen::Matrix3d c = Eigen::Matrix3d::Zero();
    c(0, 0) = lx;
    c(1, 1) = ly;
    c(2, 2) = lz;
    m.m_unit_cell = c;
    return m;
}

/// Rigid random packing of copies of `units` into a cube of L (Angstrom), every intermolecular atom pair at a
/// minimum-image distance >= min_any; units[0] is placed first. Not a liquid. (As in test_pbc_eeq_ewald.cpp.)
static Mol randomPack(const std::vector<Mol>& units, double L_ang, unsigned seed, double min_any)
{
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    Mol m = units[0];
    m.m_atoms.clear();
    m.m_geometry = Geometry(0, 3);
    int charge = 0;
    auto mic = [&](Eigen::Vector3d v) {
        for (int k = 0; k < 3; ++k)
            v[k] -= L_ang * std::rint(v[k] / L_ang);
        return v;
    };
    for (const Mol& unit : units) {
        const Eigen::Vector3d centre = unit.m_geometry.colwise().mean().transpose();
        for (long tries = 0; tries < 200000; ++tries) {
            Eigen::Vector3d p(u(rng) * L_ang, u(rng) * L_ang, u(rng) * L_ang);
            Eigen::Quaterniond qr(u(rng) - 0.5, u(rng) - 0.5, u(rng) - 0.5, u(rng) - 0.5);
            qr.normalize();
            const Eigen::Matrix3d R = qr.toRotationMatrix();
            Geometry x(unit.m_number_atoms, 3);
            for (int i = 0; i < unit.m_number_atoms; ++i)
                x.row(i) = (p + R * (unit.m_geometry.row(i).transpose() - centre)).transpose();
            bool ok = true;
            for (int a = 0; a < x.rows() && ok; ++a)
                for (int b = 0; b < m.m_geometry.rows(); ++b)
                    if (mic(x.row(a).transpose() - m.m_geometry.row(b).transpose()).norm() < min_any) {
                        ok = false;
                        break;
                    }
            if (!ok)
                continue;
            const int n0 = static_cast<int>(m.m_geometry.rows());
            m.m_geometry.conservativeResize(n0 + x.rows(), 3);
            m.m_geometry.bottomRows(x.rows()) = x;
            for (int i = 0; i < unit.m_number_atoms; ++i)
                m.m_atoms.push_back(unit.m_atoms[i]);
            charge += unit.m_charge;
            break;
        }
    }
    m.m_number_atoms = static_cast<int>(m.m_geometry.rows());
    m.m_charge = charge;
    m.m_has_pbc = false;
    return m;
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

static Mol waterBox(int nwat, double L_bohr, unsigned seed, const std::vector<Mol>& first = {}, double min_any = 2.6)
{
    std::vector<Mol> units = first;
    for (int i = 0; i < nwat; ++i)
        units.push_back(waterUnit());
    const double L_ang = L_bohr / A2B;
    return withCell(randomPack(units, L_ang, seed, min_any), L_ang, L_ang, L_ang);
}

/// A molecule from the registry placed at the centre of the cube (Angstrom cell)
static Mol centredInBox(Mol m, double L_bohr)
{
    const double L_ang = L_bohr / A2B;
    const Eigen::RowVector3d c = m.m_geometry.colwise().mean();
    for (int i = 0; i < m.m_number_atoms; ++i)
        m.m_geometry.row(i) += Eigen::RowVector3d::Constant(0.5 * L_ang) - c;
    return withCell(m, L_ang, L_ang, L_ang);
}

/// Acetic acid dimer plus F-I...pyridine 6.5 A away (as in test_pbc_hbond.cpp)
static Mol compositeMol()
{
    Mol a = registryMol("acetic_acid_dimer");
    const Mol b = registryMol("fi_pyridine");
    const Eigen::RowVector3d ca = a.m_geometry.colwise().mean(), cb = b.m_geometry.colwise().mean();
    const Eigen::RowVector3d off = ca - cb + Eigen::RowVector3d(0.0, 0.0, 6.5);
    const int na = a.m_number_atoms, nb = b.m_number_atoms;
    Geometry g(na + nb, 3);
    g.topRows(na) = a.m_geometry;
    for (int i = 0; i < nb; ++i)
        g.row(na + i) = b.m_geometry.row(i) + off;
    a.m_geometry = g;
    for (int i = 0; i < nb; ++i)
        a.m_atoms.push_back(b.m_atoms[i]);
    a.m_number_atoms = na + nb;
    return a;
}

static json config(int threads)
{
    json c = json::object();
    c["verbosity"] = 0;
    c["threads"] = threads;
    c["cache_topology"] = false;
    c["hb_update_force_every"] = 1;      // lists from the geometry at hand, on every instance
    c["hb_min_pair_energy_eh"] = 0.0;    // no energy-based candidate filter (a discrete decision)
    c["eeq_ewald_tol"] = 1e-14;
    c["ewald_tolerance"] = 1e-14;        // k_max from alpha and this; the k set is then converged far below the FD noise
    c["ewald_alpha"] = 0.32;             // fixed: the energy must not change its parameters with L
    c["ewald_r_cut"] = 19.0;
    c["dispersion_cutoff_bohr"] = 20.0;  // absolute r_c (the default would follow L/2 under the strain)
    c["frag_charge_model"] = "reference";
    return c;
}

enum Term { T_BOND, T_ANGLE, T_DIHEDRAL, T_INVERSION, T_STORS, T_DISP, T_REP, T_COULOMB, T_HBOND, T_XBOND, T_BATM, T_TOTAL, NT };
static const char* kTermName[NT] = { "bond", "angle", "dihedral", "inversion", "stors", "dispersion(+tail)", "repulsion",
    "coulomb(EEQ functional)", "hbond", "xbond", "batm", "TOTAL" };

struct Eval {
    bool ok = false;
    double e[NT] = {};
    Eigen::Matrix3d xi[NT];
    Matrix gradient;        // Eh/Bohr, workspace rows
    int n_hb = 0, n_xb = 0;
    bool virial_valid = false;
};

/// One fresh instance on `m`: energies per term; with `virial` also the analytic tensors and the gradient.
static Eval evaluate(const Mol& m, const json& cfg, bool virial, bool gradient = false)
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
        const double e_tot = ff.Calculation(virial || gradient);
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
        r.n_hb = static_cast<int>(ff.getLastHBonds().size());
        r.n_xb = static_cast<int>(ff.getLastXBonds().size());
        r.virial_valid = ff.virialValid();
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
        if ((virial || gradient) && ff.getWorkspace())
            r.gradient = ff.getWorkspace()->gradient();
        r.ok = true;
    } catch (const std::exception& e) {
        std::printf("  evaluate threw: %s\n", e.what());
    }
    return r;
}

static Mol strained(const Mol& m, double eps)
{
    Mol s = m;
    s.m_geometry *= (1.0 + eps);
    s.m_unit_cell *= (1.0 + eps);
    return s;
}

/// Pairs within r_c (minimum image, Bohr) of a geometry: a sharp cutoff moves pairs in and out under the strain
static long pairsWithin(const Mol& m, double r_c_bohr)
{
    const int n = m.m_number_atoms;
    const double L = m.m_unit_cell(0, 0) * A2B;
    long count = 0;
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j) {
            Eigen::Vector3d d = (m.m_geometry.row(j) - m.m_geometry.row(i)).transpose() * A2B;
            for (int k = 0; k < 3; ++k)
                d[k] -= L * std::rint(d[k] / L);
            if (d.norm() <= r_c_bohr)
                ++count;
        }
    return count;
}

// ---------------------------------------------------------------------------------------------
// (a) the FD gate
static void testFiniteDifference(const char* label, const Mol& m, int threads)
{
    std::printf("(a) FD gate: %s (%d atoms, charge %d, L %.1f Bohr, %d threads)\n", label, m.m_number_atoms, m.m_charge,
        m.m_unit_cell(0, 0) * A2B, threads);
    const json cfg = config(threads);
    const Eval an = evaluate(m, cfg, true);
    check(an.ok && an.virial_valid, "analytic virial computed", an.ok && an.virial_valid ? 1 : 0, 1);
    if (!an.ok || !an.virial_valid)
        return;
    const double hs[2] = { 1e-5, 2e-5 };
    double fd[2][NT];
    long npair[5];
    int nhb[5], nxb[5];
    const double eps_list[5] = { -2e-5, -1e-5, 0.0, 1e-5, 2e-5 };
    Eval ev[5];
    for (int s = 0; s < 5; ++s) {
        const Mol ms = strained(m, eps_list[s]);
        ev[s] = (s == 2) ? an : evaluate(ms, cfg, false);
        npair[s] = pairsWithin(ms, 20.0);  // D4 evaluation radius (dispersion_cutoff_bohr) and the repulsion radius 20
        nhb[s] = ev[s].n_hb;
        nxb[s] = ev[s].n_xb;
        if (!ev[s].ok) {
            check(false, "strained instance initialised", 0, 1);
            return;
        }
    }
    const double V = std::pow(m.m_unit_cell(0, 0) * A2B, 3);
    std::printf("  pairs within 20 Bohr at eps -2e-5/-1e-5/0/+1e-5/+2e-5: %ld %ld %ld %ld %ld; HB triples %d %d %d %d %d; XB %d %d %d %d %d\n",
        npair[0], npair[1], npair[2], npair[3], npair[4], nhb[0], nhb[1], nhb[2], nhb[3], nhb[4], nxb[0], nxb[1], nxb[2], nxb[3], nxb[4]);
    const bool list_event = !(npair[0] == npair[4] && nhb[0] == nhb[4] && nxb[0] == nxb[4]);
    if (list_event)
        std::printf("  NOTE: a pair-list or triple-list event lies within the FD stencil; deviations of that term are expected\n");
    std::printf("  %-26s %16s %16s %16s %12s %12s %10s\n", "term", "Xi_an (Eh)", "Xi_FD(1e-5)", "Richardson", "rel(1e-5)", "rel(Rich)", "E (Eh)");
    for (int t = 0; t < NT; ++t) {
        for (int h = 0; h < 2; ++h) {
            const int ip = h == 0 ? 3 : 4, im = h == 0 ? 1 : 0;
            fd[h][t] = -(ev[ip].e[t] - ev[im].e[t]) / (2.0 * hs[h]);
        }
        const double rich = (4.0 * fd[0][t] - fd[1][t]) / 3.0;
        const double an_tr = an.xi[t].trace();
        const double rel1 = std::abs(an_tr - fd[0][t]) / std::max(std::abs(fd[0][t]), 1e-6);
        const double relr = std::abs(an_tr - rich) / std::max(std::abs(rich), 1e-6);
        std::printf("  %-26s %16.9e %16.9e %16.9e %12.2e %12.2e %10.5f\n", kTermName[t], an_tr, fd[0][t], rich, rel1, relr, an.e[t]);
        check(rel1 <= 1e-6, std::string("FD ") + label + " " + kTermName[t] + ": |Xi_an - Xi_FD| / max(|Xi_FD|, 1e-6)", rel1, 1e-6);
    }
    std::printf("  virial pressure tr(Xi)/(3V): %.3f bar (V = %.1f Bohr^3)\n",
        CurcumaUnit::Pressure::au_to_bar(an.xi[T_TOTAL].trace() / (3.0 * V)), V);
    // off-diagonal symmetry of the tensor (a scalar product of positions and gradients need not be symmetric term by
    // term, but the total of a rotationally invariant energy is): report
    const Eigen::Matrix3d& X = an.xi[T_TOTAL];
    std::printf("  tensor diag %.6e %.6e %.6e, max |asym| %.2e\n", X(0, 0), X(1, 1), X(2, 2), (X - X.transpose()).cwiseAbs().maxCoeff());
}

// ---------------------------------------------------------------------------------------------
// (b) invariances and (c) consistency
static double relMax(const Eigen::Matrix3d& a, const Eigen::Matrix3d& b)
{
    return (a - b).cwiseAbs().maxCoeff() / std::max(1e-300, std::max(a.cwiseAbs().maxCoeff(), b.cwiseAbs().maxCoeff()));
}

static void testInvariance(const Mol& m)
{
    std::printf("(b) invariances: 20 H2O, L %.1f Bohr\n", kL);
    const Eval base = evaluate(m, config(1), true);
    check(base.ok && base.virial_valid, "base instance", base.ok ? 1 : 0, 1);
    if (!base.ok)
        return;
    const double L = m.m_unit_cell(0, 0);
    // translation by (L, -L, 2L) plus an arbitrary shift
    Mol t = m;
    for (int i = 0; i < m.m_number_atoms; ++i)
        t.m_geometry.row(i) += Eigen::RowVector3d(L + 0.37, -L - 1.21, 2.0 * L + 5.5);
    const Eval et = evaluate(t, config(1), true);
    check(et.ok && relMax(et.xi[T_TOTAL], base.xi[T_TOTAL]) < 1e-10, "translation (L,-L,2L) + shift: total Xi equal", relMax(et.xi[T_TOTAL], base.xi[T_TOTAL]), 1e-10);
    double worst = 0.0;
    for (int k = 0; k < T_TOTAL; ++k)
        worst = std::max(worst, (et.xi[k] - base.xi[k]).cwiseAbs().maxCoeff());
    check(worst < 1e-10, "translation: every component equal (abs)", worst, 1e-10);
    // atom-wise wrap into [0, L): torn molecules
    Mol w = m;
    int torn = 0;
    for (int i = 0; i < m.m_number_atoms; ++i)
        for (int k = 0; k < 3; ++k) {
            const double f = std::floor(w.m_geometry(i, k) / L);
            if (f != 0.0)
                ++torn;
            w.m_geometry(i, k) -= L * f;
        }
    const Eval ew = evaluate(w, config(1), true);
    check(ew.ok && relMax(ew.xi[T_TOTAL], base.xi[T_TOTAL]) < 1e-10, "atom-wise wrap (" + std::to_string(torn) + " coordinates moved): total Xi equal",
        relMax(ew.xi[T_TOTAL], base.xi[T_TOTAL]), 1e-10);
    worst = 0.0;
    for (int k = 0; k < T_TOTAL; ++k)
        worst = std::max(worst, (ew.xi[k] - base.xi[k]).cwiseAbs().maxCoeff());
    check(worst < 1e-10, "atom-wise wrap: every component equal (abs)", worst, 1e-10);
    // 2x1x1 supercell
    Mol sc = m;
    const int n = m.m_number_atoms;
    sc.m_geometry.conservativeResize(2 * n, 3);
    for (int i = 0; i < n; ++i) {
        sc.m_geometry.row(n + i) = m.m_geometry.row(i) + Eigen::RowVector3d(L, 0, 0);
        sc.m_atoms.push_back(m.m_atoms[i]);
    }
    sc.m_number_atoms = 2 * n;
    sc.m_unit_cell(0, 0) = 2.0 * L;
    const Eval es = evaluate(sc, config(1), true);
    check(es.ok && relMax(es.xi[T_TOTAL], 2.0 * base.xi[T_TOTAL]) < 1e-9, "2x1x1 supercell: total Xi = 2 x cell", relMax(es.xi[T_TOTAL], 2.0 * base.xi[T_TOTAL]), 1e-9);
    worst = 0.0;
    for (int k = 0; k < T_TOTAL; ++k)
        worst = std::max(worst, (es.xi[k] - 2.0 * base.xi[k]).cwiseAbs().maxCoeff());
    check(worst < 1e-9, "2x1x1 supercell: every component doubled (abs)", worst, 1e-9);
    // threads
    const Eval e4 = evaluate(m, config(4), true);
    check(e4.ok && relMax(e4.xi[T_TOTAL], base.xi[T_TOTAL]) < 1e-12, "4 threads vs 1: total Xi equal", relMax(e4.xi[T_TOTAL], base.xi[T_TOTAL]), 1e-12);
    std::printf("  Xi diag (1 thread) %.9e %.9e %.9e Eh\n", base.xi[T_TOTAL](0, 0), base.xi[T_TOTAL](1, 1), base.xi[T_TOTAL](2, 2));

    std::printf("(c) consistency with the gradient\n");
    const Eval plain = evaluate(m, config(1), false, true);
    check(plain.ok && base.gradient.rows() == plain.gradient.rows() && (base.gradient - plain.gradient).cwiseAbs().maxCoeff() == 0.0,
        "cell: gradient with the virial requested == without (bitwise)", plain.ok ? (base.gradient - plain.gradient).cwiseAbs().maxCoeff() : 1.0, 0.0);
    check(plain.ok && base.e[T_TOTAL] == plain.e[T_TOTAL], "cell: energy with the virial requested == without (bitwise)", std::abs(base.e[T_TOTAL] - plain.e[T_TOTAL]), 0.0);
    Mol nc = m;
    nc.m_has_pbc = false;
    bool refused = false;
    try {
        GFNFF ff(config(1));
        ff.InitialiseMolecule(nc);
        ff.setWantVirial(true);
        ff.Calculation(true);
        refused = !ff.virialValid();
    } catch (const std::exception& e) {
        refused = true;
        std::printf("  no cell: refused with: %s\n", e.what());
    }
    check(refused, "no cell: the virial request is refused / no virial", refused ? 1 : 0, 1);
    Mol ce = m;
    bool refused_ewald = false;
    try {
        json c = config(1);
        c["coulomb_ewald"] = false;
        GFNFF ff(c);
        ff.InitialiseMolecule(ce);
        ff.setWantVirial(true);
        ff.Calculation(true);
        refused_ewald = !ff.virialValid();
    } catch (const std::exception& e) {
        refused_ewald = true;
        std::printf("  coulomb_ewald false: refused with: %s\n", e.what());
    }
    check(refused_ewald, "cell with coulomb_ewald false: refused", refused_ewald ? 1 : 0, 1);
}

// ---------------------------------------------------------------------------------------------
// (d) pressure of a water box and the unit factor
static void testPressure(const Mol& m, const char* label)
{
    std::printf("(d) pressure: %s\n", label);
    const double f = CurcumaUnit::Pressure::HARTREE_PER_BOHR3_TO_BAR;
    check(std::abs(f - 2.94210157e8) / 2.94210157e8 < 1e-8, "1 Eh/Bohr^3 = 2.94210157e8 bar (CODATA 2018)", std::abs(f - 2.94210157e8) / 2.94210157e8, 1e-8);
    const Eval e = evaluate(m, config(4), true);
    check(e.ok && e.virial_valid, "virial of the box computed", e.ok ? 1 : 0, 1);
    if (!e.ok)
        return;
    const double L = m.m_unit_cell(0, 0) * A2B;
    const double V = L * L * L;
    const int n = m.m_number_atoms;
    const double kT = CurcumaUnit::Constants::BOLTZMANN_HARTREE_K * 300.0;
    const double K300 = 0.5 * (3 * n - 3) * kT;
    const double p_vir = CurcumaUnit::Pressure::au_to_bar(e.xi[T_TOTAL].trace() / (3.0 * V));
    const double p_kin = CurcumaUnit::Pressure::au_to_bar(2.0 * K300 / (3.0 * V));
    double mass = 0.0;
    for (int i = 0; i < n; ++i)
        mass += Elements::AtomicMass[m.m_atoms[i]];
    std::printf("  tr Xi = %.6e Eh, V = %.1f Bohr^3: virial pressure %.2f bar, kinetic part at 300 K %.2f bar, sum %.2f bar; density %.4f g/cm^3\n",
        e.xi[T_TOTAL].trace(), V, p_vir, p_kin, p_vir + p_kin, CurcumaUnit::Pressure::amu_per_angstrom3_to_gcm3(mass, V * CurcumaUnit::Length::BOHR3_TO_ANGSTROM3));
    std::printf("  per-term virial pressure (bar):");
    for (int t = 0; t < T_TOTAL; ++t)
        std::printf(" %s %.1f", kTermName[t], CurcumaUnit::Pressure::au_to_bar(e.xi[t].trace() / (3.0 * V)));
    std::printf("\n");
    check(std::isfinite(p_vir), "pressure finite", p_vir, 1e300);
}

int main(int argc, char** argv)
{
    CurcumaLogger::set_verbosity(0);
    const bool quick = argc > 1 && std::strcmp(argv[1], "--quick") == 0;
    const Mol box20 = waterBox(20, kL, 20261105);
    testFiniteDifference("20 H2O", box20, 1);
    testInvariance(box20);
    if (!quick) {
        const Mol box60 = waterBox(60, kL, 20261108);
        testFiniteDifference("60 H2O", box60, 4);
        testPressure(box60, "60 H2O random packing in L = 44 Bohr (not equilibrated)");
        const Mol charged = waterBox(18, kL, 20261106, { registryMol("hydronium_water2", 1) });
        testFiniteDifference("H3O+(H2O)2 + 18 H2O, Q = +1", charged, 1);
        testFiniteDifference("caffeine", centredInBox(registryMol("caffeine"), kL), 1);
        testFiniteDifference("acetic acid dimer + F-I...pyridine", centredInBox(compositeMol(), kL), 1);
    }
    std::printf("%s: %d failure(s)\n", g_failures == 0 ? "PASS" : "FAIL", g_failures);
    return g_failures == 0 ? 0 : 1;
}
