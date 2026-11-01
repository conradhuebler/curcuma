/*
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 */
/**
 * Periodic GFN-FF repulsion and bonded terms (WP-PERIODIC-NONBONDED A4a) - Claude Generated (Nov 2026).
 *
 * Checks ONLY the terms that A4a makes periodic: bond, angle, torsion (incl. extra torsions and sTors),
 * inversion, bonded and non-bonded repulsion. They are read as single-term energies
 * (GFNFF::BondEnergy() ...) and workspace gradient components; HB/XB are switched off (not periodic,
 * A4b), Coulomb/EEQ and the topology are not periodic and enter no checked number.
 *
 * Topology is built from plain distances (A4b), so a molecule that is torn in the coordinates cannot
 * be initialised. Torn geometries are therefore handed to an instance initialised with the intact
 * molecule (UpdateMolecule, topology_mode constant), which keeps the bond list and term lists.
 *
 * (a) Water dimer and caffeine in a 44 Bohr cube: intact in the cell centre, whole molecule shifted
 *     by (L, -L, 2L), centred on the cell corner (intact, straddling three faces), torn by a per-atom
 *     wrap into [0, L), and one atom moved by L. Every term energy equals the centred one (< 1e-12 Eh)
 *     and the term gradients (bond, angle, torsion, bonded repulsion, total minus repulsion-off)
 *     agree per atom (< 1e-10 Eh/Bohr); without a cell the torn geometry differs (test not vacuous).
 * (b) 40 random waters in a 41 Bohr cube: E_nonbonded_rep equals a brute-force sum over the 27
 *     neighbouring images (r_c 20 < L/2) with the per-pair parameters of the workspace list, for the
 *     all-pairs list and the periodic cell list (nb_cell_list_min_atoms 0, build radius 22 > L/2),
 *     1 and 4 threads.
 * (c) 2x1x1 supercell of (b) and of the corner caffeine: repulsion and bonded terms exactly twice.
 * (d) Torn caffeine: central finite difference (Angstrom, h = 2e-4 and 1e-4 A) for the atoms nearest
 *     a face of E_angle, E_torsion + E_inversion + E_sTors and E_rep (bonded + non-bonded) against
 *     the analytic term gradients; error falls quadratically in h.
 * (e) Repulsion Verlet list (cell list): wrapping one molecule by +L does not rebuild and keeps the
 *     energy; a move beyond the half-skin rebuilds.
 *
 * Geometries: registry "H2O_dimer" (water-dimer.roo2.8) and "caffeine" (library caffeine), only
 * translated, rotated and imaged. Exit 0 when every check passes; each failure prints one line.
 */
#include "src/core/curcuma_logger.h"
#include "src/core/energy_calculators/ff_methods/gfnff.h"
#include "src/core/molecule.h"
#include "src/core/units.h"
#include "src/tools/pbc_utils.h"
#include "test_molecule_registry.h"
#include "json.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <memory>
#include <random>
#include <string>
#include <utility>

using json = nlohmann::json;

static int g_failures = 0;
static void check(bool ok, const std::string& what)
{
    if (!ok) {
        ++g_failures;
        std::printf("FAIL: %s\n", what.c_str());
    }
}
static std::string sci(double v)
{
    char b[32];
    std::snprintf(b, sizeof(b), "%.3e", v);
    return b;
}

static const double A2B = CurcumaUnit::Length::ANGSTROM_TO_BOHR;

static Mol registryMol(const char* name)
{
    return TestMolecules::TestMoleculeRegistry::createMolecule(name, false).getMolInfo();
}

/// Library molecule rigidly rotated off the axes and every atom displaced by up to 0.05 A (fixed
/// seed): the library caffeine is planar, where the inversion energy and its first derivative vanish.
static Mol distortedMol(const char* name)
{
    Mol m = registryMol(name);
    std::mt19937 rng(4242);
    std::uniform_real_distribution<double> u(-0.05, 0.05);
    const Eigen::Matrix3d R = Eigen::Quaterniond(0.9, 0.2, -0.3, 0.25).normalized().toRotationMatrix();
    for (int i = 0; i < m.m_number_atoms; ++i) {
        Eigen::Vector3d x = R * m.m_geometry.row(i).transpose();
        for (int k = 0; k < 3; ++k)
            m.m_geometry(i, k) = x[k] + u(rng);
    }
    return m;
}

// ---------------------------------------------------------------------------------------------
struct Terms {
    double bond = 0, angle = 0, torsion = 0, inversion = 0, stors = 0, brep = 0, nbrep = 0;
    double v(int k) const
    {
        const double a[7] = { bond, angle, torsion, inversion, stors, brep, nbrep };
        return a[k];
    }
};
static const char* kTermNames[7] = { "bond", "angle", "torsion", "inversion", "sTors", "bonded rep", "non-bonded rep" };

struct Grads {
    Matrix total, bond, angle, torsion, brep;
};

static json config(int threads, const json& extra)
{
    json c = json::object();
    c["verbosity"] = 0;
    c["threads"] = threads;
    c["cache_topology"] = false;
    c["hbond"] = false;  // HB/XB not periodic (A4b); also removes the HB modulation of the X-H bonds
    for (auto it = extra.begin(); it != extra.end(); ++it)
        c[it.key()] = it.value();
    return c;
}

static Mol withCell(Mol mol, const Eigen::Vector3d& L_ang)
{
    mol.m_has_pbc = L_ang.minCoeff() > 0.0;
    if (mol.m_has_pbc)
        mol.m_unit_cell = L_ang.asDiagonal();
    return mol;
}

/// One GFN-FF instance, initialised once; evaluate() hands it a geometry (Angstrom).
struct Instance {
    std::unique_ptr<GFNFF> ff;
    bool ok = false;
    Instance(const Mol& mol, const json& c)
    {
        try {
            ff = std::make_unique<GFNFF>(c);
            ok = ff->InitialiseMolecule(mol);
            if (ok && ff->getWorkspace())
                ff->getWorkspace()->setStoreGradientComponents(true);
        } catch (const std::exception& e) {
            std::printf("init threw: %s\n", e.what());
            ok = false;
        }
    }
    Terms evaluate(const Geometry& g_ang, bool gradient, Grads* grads = nullptr)
    {
        ff->UpdateMolecule(g_ang);
        ff->Calculation(gradient);
        Terms t;
        t.bond = ff->BondEnergy();
        t.angle = ff->AngleEnergy();
        t.torsion = ff->DihedralEnergy();
        t.inversion = ff->InversionEnergy();
        t.stors = ff->STorsEnergy();
        t.brep = ff->BondedRepulsionEnergy();
        t.nbrep = ff->NonbondedRepulsionEnergy();
        if (grads && gradient) {
            grads->total = ff->Gradient();
            const FFWorkspace* ws = ff->getWorkspace();
            grads->bond = ws->gradientBond();
            grads->angle = ws->gradientAngle();
            grads->torsion = ws->gradientTorsion();
            grads->brep = ws->gradientRepulsion();
        }
        return t;
    }
};

static double maxTermDiff(const Terms& a, const Terms& b, int* which = nullptr)
{
    double d = 0.0;
    for (int k = 0; k < 7; ++k) {
        const double x = std::abs(a.v(k) - b.v(k));
        if (x > d) {
            d = x;
            if (which)
                *which = k;
        }
    }
    return d;
}

static Geometry centred(const Mol& m, const Eigen::Vector3d& at)
{
    Geometry g = m.m_geometry;
    Eigen::RowVector3d c = g.colwise().mean();
    for (int i = 0; i < g.rows(); ++i)
        g.row(i) += at.transpose() - c;
    return g;
}

// ---------------------------------------------------------------------------------------------
// (a) image invariance of every term, torn input
static const double kCellBohr = 44.0;

static void testInvariance(const char* name)
{
    const Mol base = distortedMol(name);
    const int n = base.m_number_atoms;
    const double L = kCellBohr / A2B;
    const Eigen::Vector3d Lv(L, L, L);
    const Geometry g_mid = centred(base, 0.5 * Lv);
    const Geometry g_corner = centred(base, Eigen::Vector3d::Zero());

    Geometry g_shift = g_mid;
    for (int i = 0; i < n; ++i) {
        g_shift(i, 0) += L;
        g_shift(i, 1) -= L;
        g_shift(i, 2) += 2.0 * L;
    }
    Geometry g_torn = g_corner;  // every atom wrapped separately: bonds cross the faces in the coordinates
    for (int i = 0; i < n; ++i)
        for (int k = 0; k < 3; ++k)
            g_torn(i, k) -= L * std::floor(g_torn(i, k) / L);
    Geometry g_one = g_mid;  // one atom by one lattice vector
    g_one(n / 2, 1) += L;

    Mol m_mid = withCell(base, Lv);
    m_mid.m_geometry = g_mid;
    Instance ref(m_mid, config(1, json::object()));
    check(ref.ok, std::string(name) + ": init in the cell centre");
    if (!ref.ok)
        return;
    Grads gref;
    const Terms tref = ref.evaluate(g_mid, true, &gref);

    // no cell, same geometry: MIC of short vectors changes nothing (rounding only)
    {
        Mol m0 = base;
        m0.m_geometry = g_mid;
        Instance nc(m0, config(1, json::object()));
        const Terms t0 = nc.ok ? nc.evaluate(g_mid, false) : Terms();
        const double d = maxTermDiff(tref, t0);
        check(nc.ok && d < 1e-12, std::string(name) + ": centred with cell == without cell (" + sci(d) + ")");
    }

    // repulsion-off twin to isolate the non-bonded repulsion gradient
    json c_off = json::object();
    c_off["repulsion"] = false;
    Instance ref_off(m_mid, config(1, c_off));
    Grads gref_off;
    ref_off.evaluate(g_mid, true, &gref_off);
    const Matrix g_nb_ref = gref.total - gref_off.total - gref.brep;

    // fresh instances for intact geometries; the corner instance takes the torn ones
    Mol m_corner = withCell(base, Lv);
    m_corner.m_geometry = g_corner;
    Instance corner(m_corner, config(1, json::object()));
    Instance corner_off(m_corner, config(1, c_off));
    Mol m_shift = withCell(base, Lv);
    m_shift.m_geometry = g_shift;
    Instance shifted(m_shift, config(1, json::object()));
    check(corner.ok && corner_off.ok && shifted.ok, std::string(name) + ": init of the image variants");
    if (!corner.ok || !corner_off.ok || !shifted.ok)
        return;

    struct V {
        const char* what;
        Instance* inst;
        Instance* off;
        const Geometry* g;
    };
    const V variants[] = { { "shifted by (L,-L,2L)", &shifted, nullptr, &g_shift },
        { "centred on the corner", &corner, &corner_off, &g_corner },
        { "torn (per-atom wrap)", &corner, &corner_off, &g_torn },
        { "one atom moved by L", &corner, &corner_off, &g_one } };
    for (const auto& v : variants) {
        Grads gv, gv_off;
        const Terms t = v.inst->evaluate(*v.g, true, &gv);
        int k = 0;
        const double d = maxTermDiff(t, tref, &k);
        check(d < 1e-12, std::string(name) + " " + v.what + ": all terms == centred (max " + sci(d) + " in " + kTermNames[k] + ")");
        double dg = std::max({ (gv.bond - gref.bond).cwiseAbs().maxCoeff(), (gv.angle - gref.angle).cwiseAbs().maxCoeff(),
            (gv.torsion - gref.torsion).cwiseAbs().maxCoeff(), (gv.brep - gref.brep).cwiseAbs().maxCoeff() });
        if (v.off) {
            v.off->evaluate(*v.g, true, &gv_off);
            const Matrix g_nb = gv.total - gv_off.total - gv.brep;
            dg = std::max(dg, (g_nb - g_nb_ref).cwiseAbs().maxCoeff());
        }
        check(dg < 1e-10, std::string(name) + " " + v.what + ": term gradients per atom == centred (" + sci(dg) + ")");
        std::printf("%-9s %-22s max |dE| %.2e Eh (%s), max |dG| %.2e Eh/Bohr\n", name, v.what, d, kTermNames[k], dg);
    }
    std::printf("%-9s terms (centred): bond %.12f angle %.12f tors %.12f inv %.12f stors %.3e brep %.12f nbrep %.12f\n", name,
        tref.bond, tref.angle, tref.torsion, tref.inversion, tref.stors, tref.brep, tref.nbrep);

    // not vacuous: the same geometries without a cell (same topology, from the intact corner
    // geometry) are far off. The per-atom wrap tears a molecule only if it straddles a face, which
    // the corner placement does not guarantee for every molecule, so only "one atom by L" must differ.
    {
        Mol m0 = base;
        m0.m_geometry = g_corner;
        Instance nc(m0, config(1, json::object()));
        const double d_torn = nc.ok ? maxTermDiff(nc.evaluate(g_torn, false), tref) : 0.0;
        const double d_one = nc.ok ? maxTermDiff(nc.evaluate(g_one, false), tref) : 0.0;
        check(nc.ok && d_one > 1e-3, std::string(name) + ": one atom moved by L without a cell differs (" + sci(d_one) + ")");
        std::printf("%-9s without a cell: torn max |dE| %.3e Eh, one atom by L %.3e Eh (test not vacuous)\n", name, d_torn, d_one);
    }
}

// ---------------------------------------------------------------------------------------------
// random water box (intact molecules, minimum-image O-O >= 2.9 A)
static Mol waterBox(int nwat, double L_ang, unsigned seed)
{
    const Mol d = registryMol("H2O_dimer");
    Eigen::MatrixXd w(3, 3);
    for (int i = 0; i < 3; ++i)
        w.row(i) = d.m_geometry.row(i) - d.m_geometry.row(0);
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    std::vector<Eigen::Vector3d> O;
    Mol m = d;
    m.m_atoms.clear();
    m.m_geometry = Geometry(3 * nwat, 3);
    int placed = 0, tries = 0;
    while (placed < nwat && tries < 200000) {
        ++tries;
        Eigen::Vector3d p(u(rng) * L_ang, u(rng) * L_ang, u(rng) * L_ang);
        bool ok = true;
        for (const auto& q : O) {
            Eigen::Vector3d dv = p - q;
            for (int k = 0; k < 3; ++k)
                dv[k] -= L_ang * std::rint(dv[k] / L_ang);
            if (dv.norm() < 2.9) {
                ok = false;
                break;
            }
        }
        if (!ok)
            continue;
        Eigen::Quaterniond qr(u(rng) - 0.5, u(rng) - 0.5, u(rng) - 0.5, u(rng) - 0.5);
        qr.normalize();
        const Eigen::Matrix3d R = qr.toRotationMatrix();
        for (int i = 0; i < 3; ++i) {
            m.m_geometry.row(3 * placed + i) = (p + R * w.row(i).transpose()).transpose();
            m.m_atoms.push_back(d.m_atoms[i]);
        }
        O.push_back(p);
        ++placed;
    }
    m.m_number_atoms = 3 * placed;
    m.m_geometry.conservativeResize(3 * placed, 3);
    check(placed == nwat, "water box placed all molecules");
    return m;
}

/// Brute-force non-bonded repulsion over the 27 neighbouring images, radius rc (Bohr).
static double bruteForceRep(const Geometry& g_ang, const Eigen::Vector3d& L_bohr, double rc,
    const std::vector<GFNFFRepulsion>& pairs, double* central_only)
{
    const Eigen::MatrixXd g = g_ang * A2B;
    double e = 0.0, e0 = 0.0;
    for (const auto& p : pairs)
        for (int a = -1; a <= 1; ++a)
            for (int b = -1; b <= 1; ++b)
                for (int c = -1; c <= 1; ++c) {
                    const double dx = g(p.i, 0) - g(p.j, 0) - a * L_bohr[0];
                    const double dy = g(p.i, 1) - g(p.j, 1) - b * L_bohr[1];
                    const double dz = g(p.i, 2) - g(p.j, 2) - c * L_bohr[2];
                    const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
                    if (r > rc || r < 1e-8)
                        continue;
                    const double ep = p.repab * std::exp(-p.alpha * r * std::sqrt(r)) / r;
                    e += ep;
                    if (a == 0 && b == 0 && c == 0)
                        e0 += ep;
                }
    if (central_only)
        *central_only = e0;
    return e;
}

static const double kBoxBohr = 41.0;
static const int kWaters = 40;

// (b) brute force and (c) supercell
static void testBox()
{
    const double L = kBoxBohr / A2B;
    const Eigen::Vector3d Lv(L, L, L);
    const Mol box = withCell(waterBox(kWaters, L, 20261104), Lv);
    // all-pairs reference list (N < nb_cell_list_min_atoms): every non-bonded pair once
    Instance all1(box, config(1, json::object()));
    check(all1.ok, "box: init (all-pairs list)");
    if (!all1.ok)
        return;
    const Terms t_all = all1.evaluate(box.m_geometry, false);
    const std::vector<GFNFFRepulsion> pairs = all1.ff->getWorkspace()->nonbondedRepulsions();
    double e_central = 0.0;
    const double eb = bruteForceRep(box.m_geometry, Lv * A2B, 20.0, pairs, &e_central);
    const double d = std::abs(t_all.nbrep - eb);
    check(d < 1e-12, "box: E_nbrep (all pairs) == 27-image brute force (" + sci(d) + ")");
    check(std::abs(eb - e_central) > 1e-6, "box: images contribute (test not vacuous): " + sci(eb - e_central));

    json cl = json::object();
    cl["nb_cell_list_min_atoms"] = 0;  // periodic cell list, build radius 20 + 2 > L/2
    double e_cl[2] = { 0.0, 0.0 };
    size_t n_cl = 0;
    const int threads[2] = { 1, 4 };
    for (int t = 0; t < 2; ++t) {
        Instance inst(box, config(threads[t], cl));
        check(inst.ok, "box: init (periodic cell list)");
        if (!inst.ok)
            return;
        e_cl[t] = inst.evaluate(box.m_geometry, false).nbrep;
        n_cl = inst.ff->getWorkspace()->nonbondedRepulsions().size();
    }
    const double dcl = std::abs(e_cl[0] - eb), dt = std::abs(e_cl[0] - e_cl[1]);
    check(dcl < 1e-12, "box: E_nbrep (periodic cell list) == brute force (" + sci(dcl) + ")");
    check(dt < 1e-14, "box: cell list 1 and 4 threads agree (" + sci(dt) + ")");
    std::printf("box %d atoms, L %.1f Bohr: E_nbrep all-pairs %.15f, cell list %.15f, brute force %.15f, |d| %.2e / %.2e, "
                "image share %.3e Eh, pairs %zu (all) / %zu (cell list), t1-t4 %.1e\n",
        box.m_number_atoms, kBoxBohr, t_all.nbrep, e_cl[0], eb, d, dcl, eb - e_central, pairs.size(), n_cl, dt);

    // cell check: with a cell only the repulsion evaluation radius 20 is checked (the cell-list
    // build above used 22 > L/2 = 20.5); L = 39.8 Bohr (L/2 < 20) is refused
    {
        const Mol small = withCell(box, Eigen::Vector3d::Constant(39.8 / A2B));
        Instance refused(small, config(1, cl));  // prints the cell-check message ("init threw")
        check(!refused.ok, "cell with L/2 < 20 Bohr (repulsion radius) is refused");
    }

    // (c) 2x1x1 supercell of the box
    {
        Mol sc = box;
        const int n = box.m_number_atoms;
        sc.m_geometry = Geometry(2 * n, 3);
        sc.m_atoms.clear();
        for (int r = 0; r < 2; ++r)
            for (int i = 0; i < n; ++i) {
                sc.m_geometry.row(r * n + i) = box.m_geometry.row(i);
                sc.m_geometry(r * n + i, 0) += r * L;
                sc.m_atoms.push_back(box.m_atoms[i]);
            }
        sc.m_number_atoms = 2 * n;
        sc = withCell(sc, Eigen::Vector3d(2.0 * L, L, L));
        for (int variant = 0; variant < 2; ++variant) {
            Instance s(sc, config(1, variant == 0 ? json::object() : cl));
            check(s.ok, "supercell runs");
            if (!s.ok)
                continue;
            const Terms ts = s.evaluate(sc.m_geometry, false);
            const double dd = std::max({ std::abs(ts.nbrep - 2.0 * t_all.nbrep), std::abs(ts.brep - 2.0 * t_all.brep),
                std::abs(ts.bond - 2.0 * t_all.bond), std::abs(ts.angle - 2.0 * t_all.angle) });
            check(dd < 1e-12, std::string("water supercell 2x1x1 (") + (variant ? "cell list" : "all pairs") + "): terms == 2 x cell (" + sci(dd) + ")");
            std::printf("water supercell 2x1x1 (%s): E_nbrep %.15f, 2 x cell %.15f, max |d| over nbrep/brep/bond/angle %.2e\n",
                variant ? "cell list" : "all pairs", ts.nbrep, 2.0 * t_all.nbrep, dd);
        }
    }
}

// (c) caffeine 2x1x1 supercell (torsions, inversions), molecule on the corner of the small cell
static void testCaffeineSupercell()
{
    const Mol base = distortedMol("caffeine");
    const int n = base.m_number_atoms;
    const double L = kCellBohr / A2B;
    Mol cell = withCell(base, Eigen::Vector3d(L, L, L));
    cell.m_geometry = centred(base, Eigen::Vector3d::Zero());
    Instance c1(cell, config(1, json::object()));
    Mol sc = cell;
    sc.m_geometry = Geometry(2 * n, 3);
    sc.m_atoms.clear();
    for (int r = 0; r < 2; ++r)
        for (int i = 0; i < n; ++i) {
            sc.m_geometry.row(r * n + i) = cell.m_geometry.row(i);
            sc.m_geometry(r * n + i, 0) += r * L;
            sc.m_atoms.push_back(base.m_atoms[i]);
        }
    sc.m_number_atoms = 2 * n;
    sc = withCell(sc, Eigen::Vector3d(2.0 * L, L, L));
    Instance c2(sc, config(1, json::object()));
    check(c1.ok && c2.ok, "caffeine supercell runs");
    if (!c1.ok || !c2.ok)
        return;
    const Terms t1 = c1.evaluate(cell.m_geometry, false), t2 = c2.evaluate(sc.m_geometry, false);
    double dd = 0.0;
    int which = 0;
    for (int k = 0; k < 7; ++k) {
        const double x = std::abs(t2.v(k) - 2.0 * t1.v(k));
        if (x > dd) {
            dd = x;
            which = k;
        }
    }
    check(dd < 1e-12, "caffeine supercell 2x1x1: every term == 2 x cell (" + sci(dd) + ")");
    std::printf("caffeine supercell 2x1x1: max |E - 2 E_cell| %.2e (%s)\n", dd, kTermNames[which]);
}

// ---------------------------------------------------------------------------------------------
// (d) finite difference on the torn caffeine
static void testFiniteDifference()
{
    const Mol base = distortedMol("caffeine");
    const int n = base.m_number_atoms;
    const double L = kCellBohr / A2B;
    const Eigen::Vector3d Lv(L, L, L);
    Mol m = withCell(base, Lv);
    m.m_geometry = centred(base, Eigen::Vector3d::Zero());
    Geometry g_torn = m.m_geometry;
    for (int i = 0; i < n; ++i)
        for (int k = 0; k < 3; ++k)
            g_torn(i, k) -= L * std::floor(g_torn(i, k) / L);

    Instance on(m, config(1, json::object()));
    json c_off = json::object();
    c_off["repulsion"] = false;
    Instance off(m, config(1, c_off));
    check(on.ok && off.ok, "FD: init");
    if (!on.ok || !off.ok)
        return;
    Grads ga, gb;
    on.evaluate(g_torn, true, &ga);
    off.evaluate(g_torn, true, &gb);
    const Matrix g_rep = ga.total - gb.total;  // bonded + non-bonded repulsion, Eh/Bohr

    std::vector<std::pair<double, int>> dist;
    for (int i = 0; i < n; ++i) {
        double dmin = 1e300;
        for (int k = 0; k < 3; ++k)
            dmin = std::min({ dmin, g_torn(i, k), L - g_torn(i, k) });
        dist.push_back({ dmin, i });
    }
    std::sort(dist.begin(), dist.end());
    const int nfd = 8;

    // term groups: 0 angle, 1 torsion + inversion + sTors, 2 repulsion
    const char* group[3] = { "angle", "torsion+inversion", "repulsion" };
    auto groupEnergy = [](const Terms& t, int gidx) {
        if (gidx == 0)
            return t.angle;
        if (gidx == 1)
            return t.torsion + t.inversion + t.stors;
        return t.brep + t.nbrep;
    };
    const Matrix* an[3] = { &ga.angle, &ga.torsion, &g_rep };
    double err[2][3] = { { 0, 0, 0 }, { 0, 0, 0 } }, gmax[3] = { 0, 0, 0 }, inv_fd = 0.0;
    const double hs[2] = { 2e-4, 1e-4 };
    for (int hi = 0; hi < 2; ++hi) {
        const double h = hs[hi];
        for (int s = 0; s < nfd; ++s) {
            const int a = dist[s].second;
            for (int k = 0; k < 3; ++k) {
                Geometry gp = g_torn, gm = g_torn;
                gp(a, k) += h;
                gm(a, k) -= h;
                const Terms tp = on.evaluate(gp, false), tm = on.evaluate(gm, false);
                for (int gi = 0; gi < 3; ++gi) {
                    const double fd = (groupEnergy(tp, gi) - groupEnergy(tm, gi)) / (2.0 * h);  // Eh/A
                    const double ana = (*an[gi])(a, k) * A2B;
                    err[hi][gi] = std::max(err[hi][gi], std::abs(fd - ana));
                    gmax[gi] = std::max(gmax[gi], std::abs(ana));
                }
                inv_fd = std::max(inv_fd, std::abs((tp.inversion - tm.inversion) / (2.0 * h)));
            }
        }
    }
    for (int gi = 0; gi < 3; ++gi) {
        check(err[1][gi] < 1e-7, std::string("FD ") + group[gi] + ": max error " + sci(err[1][gi]) + " Eh/A");
        check(err[1][gi] < err[0][gi] / 3.0 || err[1][gi] < 1e-10,
            std::string("FD ") + group[gi] + ": error falls quadratically (" + sci(err[0][gi]) + " -> " + sci(err[1][gi]) + ")");
        check(gmax[gi] > 1e-5, std::string("FD ") + group[gi] + ": gradient non-zero (test not vacuous)");
        std::printf("FD %-18s (%d atoms %.2f..%.2f A from a face): max |g| %.3e Eh/A, |an - FD| %.2e (h 2e-4 A), %.2e (h 1e-4 A)\n",
            group[gi], nfd, dist[0].first, dist[nfd - 1].first, gmax[gi], err[0][gi], err[1][gi]);
    }
    check(inv_fd > 1e-6, "FD: inversion energy moves with the boundary atoms (" + sci(inv_fd) + " Eh/A)");
    std::printf("FD inversion part alone: max |dE_inv/dx| %.3e Eh/A\n", inv_fd);
}

// (e) repulsion Verlet list with the periodic cell list
static void testVerlet()
{
    const double L = kBoxBohr / A2B;
    const Eigen::Vector3d Lv(L, L, L);
    const Mol box = withCell(waterBox(kWaters, L, 20261104), Lv);
    json cl = json::object();
    cl["nb_cell_list_min_atoms"] = 0;
    Instance inst(box, config(1, cl));
    check(inst.ok, "Verlet: init");
    if (!inst.ok)
        return;
    Geometry g = box.m_geometry;
    const double e0 = inst.evaluate(g, true).nbrep;
    const long c0 = inst.ff->repulsionRebuildCount();
    for (int i = 0; i < 3; ++i)
        g(i, 0) += L;
    const double e1 = inst.evaluate(g, true).nbrep;
    const long c1 = inst.ff->repulsionRebuildCount();
    check(std::abs(e1 - e0) < 1e-12, "wrap by +L: E_nbrep unchanged (" + sci(std::abs(e1 - e0)) + ")");
    check(c1 == c0, "wrap by +L: no repulsion list rebuild (" + std::to_string(c0) + " -> " + std::to_string(c1) + ")");
    for (int i = 3; i < 6; ++i)
        g(i, 1) += 1.0;  // 1.9 Bohr > half-skin 1 Bohr
    const double e2 = inst.evaluate(g, true).nbrep;
    const long c2 = inst.ff->repulsionRebuildCount();
    check(c2 == c0 + 1, "move beyond the half-skin rebuilds once (" + std::to_string(c2 - c0) + ")");
    Mol moved = box;
    moved.m_geometry = g;
    Instance fresh(moved, config(1, cl));
    const double ef = fresh.ok ? fresh.evaluate(g, false).nbrep : 0.0;
    check(fresh.ok && std::abs(e2 - ef) < 1e-12, "after the rebuild E_nbrep == fresh instance (" + sci(std::abs(e2 - ef)) + ")");
    std::printf("Verlet: wrap |dE| %.2e, rebuilds %ld -> %ld; after 1 A move rebuilds %ld, |E - fresh| %.2e\n",
        std::abs(e1 - e0), c0, c1, c2, std::abs(e2 - ef));
}

int main()
{
    CurcumaLogger::set_verbosity(0);
    testInvariance("H2O_dimer");
    testInvariance("caffeine");
    testBox();
    testCaffeineSupercell();
    testFiniteDifference();
    testVerlet();
    if (g_failures == 0)
        std::printf("pbc_bonded_repulsion: all checks passed\n");
    return g_failures == 0 ? 0 : 1;
}
