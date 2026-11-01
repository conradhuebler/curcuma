/*
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 */
/**
 * Periodic GFN-FF topology generation and reactive topology (WP-PERIODIC-NONBONDED A4c) - Claude Generated (Nov 2026).
 *
 * A4c makes the topology build periodic: the bond criterion takes minimum-image distances, every connected fragment
 * is reassembled (TopologyGeometryScope, lowest-index atom anchored) before hybridisation, angle/torsion/inversion
 * generation and the Phase-1 EEQ run on it, the reactive bond scan takes minimum-image distances, the charge-ensemble
 * contact window likewise, and the topology cache fingerprint carries the cell. Coulomb and the per-step Phase-2 EEQ
 * stay non-periodic (plain distances on the input geometry), so the Coulomb term is reported but not compared.
 *
 * (a) caffeine, hydronium-water2 (+1, default charge ensemble), acetic acid dimer in a 44 Bohr cube: a NEW instance is
 *     initialised directly with the molecule centred on the cell corner and wrapped per atom (torn). Two variants:
 *     "anchored" (the lowest-index atom of every fragment keeps its position, so the reassembled geometry is the corner
 *     geometry itself: topology charges equal too) and "fully torn" (all atoms wrapped; periodic terms and the graph
 *     topology equal, the Phase-1 charges of a multi-fragment system may differ by the non-periodic inter-fragment
 *     Coulomb of the Phase-1 EEQ, reported). Compared with an instance initialised intact in the cell centre: bond list,
 *     angle/dihedral/inversion counts, hybridisation, fragment count, every periodic energy term (bond, angle, torsion,
 *     inversion, sTors, repulsion, D4 pair + tail, HB, XB) <= 1e-12 Eh; and the graph topology, bonded terms and HB
 *     against an instance without a cell.
 * (b) 2x1x1 supercell of the centred caffeine and hydronium: counts exactly doubled, periodic terms twice the cell.
 * (c) 60 random waters (L 43.5 Bohr, liquid-like density), wrapped per atom: bond list == the one expected physically
 *     (two O-H per water, nothing else) and == an independent 27-image enumeration of the charge-shrunk bond criterion
 *     (gfnff_ini2.f90:111-126, 361-419) with the instance's own Phase-1 charges; threads 1 and 4 identical; the same
 *     geometry without a cell loses the bonds of the torn waters (test not vacuous).
 * (d) reactive topology (topology_mode react, react_check_every 1): two H atoms approach from 4.0 to 1.2 Bohr, once in
 *     the cell centre and once across the x face (minimum-image distance identical); the formation event happens at the
 *     same step, the energy after the event, the dE_jump and the bond list are equal.
 * (e) topology cache: the fingerprint of a topology written with a cell differs from the one written without a cell
 *     for the same intact geometry and bond list, and is reproduced by the same cell; a cell-made cache file is thus
 *     not read without a cell and vice versa.
 * (f) smallest cell: water dimer, hbond false: L 40.5 Bohr accepted, 39.9 refused (repulsion 20); the topology
 *     thresholds add no blocker.
 *
 * Exit 0 when every check passes; each failure prints one line.
 */
#include "src/core/curcuma_logger.h"
#include "src/core/energy_calculators/ff_methods/gfnff.h"
#include "src/core/energy_calculators/ff_methods/gfnff_par.h"
#include "src/core/molecule.h"
#include "src/core/units.h"
#include "src/tools/pbc_utils.h"
#include "test_molecule_registry.h"
#include "json.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <memory>
#include <random>
#include <set>
#include <string>
#include <utility>
#include <vector>

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

static Mol registryMol(const char* name, int charge = 0)
{
    Mol m = TestMolecules::TestMoleculeRegistry::createMolecule(name, false).getMolInfo();
    m.m_charge = charge;
    return m;
}

static Mol withCell(Mol mol, const Eigen::Vector3d& L_ang)
{
    mol.m_has_pbc = L_ang.minCoeff() > 0.0;
    if (mol.m_has_pbc)
        mol.m_unit_cell = L_ang.asDiagonal();
    return mol;
}

static json config(int threads, const json& extra = json::object())
{
    json c = json::object();
    c["verbosity"] = 0;
    c["threads"] = threads;
    c["cache_topology"] = false;
    for (auto it = extra.begin(); it != extra.end(); ++it)
        c[it.key()] = it.value();
    return c;
}

using Pair = std::pair<int, int>;

struct Topo {
    std::vector<Pair> bonds;
    size_t nangles = 0, ndihedrals = 0, ninversions = 0;
    std::vector<int> hyb;
    int nfrag = 0;
    Vector qtopo;
};

struct Terms {
    double bond = 0, angle = 0, dihedral = 0, inversion = 0, stors = 0, rep = 0, disp = 0, tail = 0, hb = 0, xb = 0, coulomb = 0, total = 0;
    double periodicSum() const { return bond + angle + dihedral + inversion + stors + rep + disp + tail + hb + xb; }
};

struct Instance {
    std::unique_ptr<GFNFF> ff;
    bool ok = false;
    std::string error;
    Instance(const Mol& mol, const json& c)
    {
        try {
            ff = std::make_unique<GFNFF>(c);
            ok = ff->InitialiseMolecule(mol);
        } catch (const std::exception& e) {
            error = e.what();
            ok = false;
        }
    }
    Topo topology() const
    {
        Topo t;
        const FFWorkspace* ws = ff->getWorkspace();
        for (const auto& b : ws->bonds())
            t.bonds.emplace_back(std::min(b.i, b.j), std::max(b.i, b.j));
        std::sort(t.bonds.begin(), t.bonds.end());
        t.nangles = ws->angles().size();
        t.ndihedrals = ws->dihedrals().size();
        t.ninversions = ws->inversions().size();
        const json tj = ff->exportTopology();
        if (tj.contains("hybridization"))
            t.hyb = tj["hybridization"].get<std::vector<int>>();
        t.nfrag = tj.value("nfrag", 0);
        t.qtopo = ff->getTopologyCharges();
        return t;
    }
    Terms evaluate(const Geometry& g_ang, bool gradient = false)
    {
        ff->UpdateMolecule(g_ang);
        Terms r;
        r.total = ff->Calculation(gradient);
        r.bond = ff->BondEnergy();
        r.angle = ff->AngleEnergy();
        r.dihedral = ff->DihedralEnergy();
        r.inversion = ff->InversionEnergy();
        r.stors = ff->STorsEnergy();
        r.rep = ff->RepulsionEnergy();
        r.disp = ff->DispersionEnergy();
        r.tail = ff->DispersionTailEnergy();
        r.hb = ff->HydrogenBondEnergy();
        r.xb = ff->HalogenBondEnergy();
        r.coulomb = ff->CoulombEnergy();
        return r;
    }
};

static double maxTermDiff(const Terms& a, const Terms& b, double scale_b = 1.0)
{
    return std::max({ std::abs(a.bond - scale_b * b.bond), std::abs(a.angle - scale_b * b.angle), std::abs(a.dihedral - scale_b * b.dihedral),
        std::abs(a.inversion - scale_b * b.inversion), std::abs(a.stors - scale_b * b.stors), std::abs(a.rep - scale_b * b.rep),
        std::abs(a.disp - scale_b * b.disp), std::abs(a.tail - scale_b * b.tail), std::abs(a.hb - scale_b * b.hb), std::abs(a.xb - scale_b * b.xb) });
}

static bool sameGraph(const Topo& a, const Topo& b, std::string& why)
{
    if (a.bonds != b.bonds) { why = "bond list (" + std::to_string(a.bonds.size()) + " vs " + std::to_string(b.bonds.size()) + ")"; return false; }
    if (a.nangles != b.nangles) { why = "angle count (" + std::to_string(a.nangles) + " vs " + std::to_string(b.nangles) + ")"; return false; }
    if (a.ndihedrals != b.ndihedrals) { why = "dihedral count (" + std::to_string(a.ndihedrals) + " vs " + std::to_string(b.ndihedrals) + ")"; return false; }
    if (a.ninversions != b.ninversions) { why = "inversion count (" + std::to_string(a.ninversions) + " vs " + std::to_string(b.ninversions) + ")"; return false; }
    if (a.hyb != b.hyb) { why = "hybridisation"; return false; }
    if (a.nfrag != b.nfrag) { why = "fragment count (" + std::to_string(a.nfrag) + " vs " + std::to_string(b.nfrag) + ")"; return false; }
    return true;
}

static double maxAbsDiff(const Vector& a, const Vector& b)
{
    if (a.size() != b.size())
        return 1e300;
    return a.size() ? (a - b).cwiseAbs().maxCoeff() : 0.0;
}

static Geometry centred(const Geometry& g0, const Eigen::Vector3d& at)
{
    Geometry g = g0;
    Eigen::RowVector3d c = g.colwise().mean();
    for (int i = 0; i < g.rows(); ++i)
        g.row(i) += at.transpose() - c;
    return g;
}

static Geometry wrapped(const Geometry& g0, double L)
{
    Geometry g = g0;
    for (int i = 0; i < g.rows(); ++i)
        for (int k = 0; k < 3; ++k)
            g(i, k) -= L * std::floor(g(i, k) / L);
    return g;
}

static std::vector<std::vector<int>> fragmentsOf(const Mol& m)
{
    curcuma::Molecule mol(m);
    return mol.GetFragments();
}

static const double kCellBohr = 44.0;

// ---------------------------------------------------------------------------------------------
// (a) torn input initialised directly
static void testTornInit(const char* label, const Mol& base, const json& extra = json::object())
{
    const int n = base.m_number_atoms;
    const double L = kCellBohr / A2B;
    const Eigen::Vector3d Lv(L, L, L);
    const Geometry g_mid = centred(base.m_geometry, 0.5 * Lv);
    const Geometry g_corner = centred(base.m_geometry, Eigen::Vector3d::Zero());
    const Geometry g_torn = wrapped(g_corner, L);
    // anchored tear: every fragment moved by the lattice vector that brings its lowest-index atom back
    Geometry g_anchored = g_torn;
    const auto frags = fragmentsOf(base);
    for (const auto& f : frags) {
        const int a0 = *std::min_element(f.begin(), f.end());
        const Eigen::RowVector3d back = g_corner.row(a0) - g_torn.row(a0);
        for (int a : f)
            g_anchored.row(a) += back;
    }

    Mol m_mid = withCell(base, Lv);
    m_mid.m_geometry = g_mid;
    Instance ref(m_mid, config(1, extra));
    check(ref.ok, std::string(label) + ": init intact in the cell centre (" + ref.error + ")");
    if (!ref.ok)
        return;
    const Topo t_ref = ref.topology();
    const Terms e_ref = ref.evaluate(g_mid);

    {  // no cell, same geometry: graph topology, bonded terms and HB must agree (D4 radius differs with a cell)
        Mol m0 = base;
        m0.m_geometry = g_mid;
        Instance nc(m0, config(1, extra));
        check(nc.ok, std::string(label) + ": init without a cell");
        if (nc.ok) {
            const Topo t0 = nc.topology();
            const Terms e0 = nc.evaluate(g_mid);
            std::string why;
            check(sameGraph(t0, t_ref, why), std::string(label) + ": topology with cell == without cell (" + why + ")");
            const double dq = maxAbsDiff(t0.qtopo, t_ref.qtopo);
            const double de = std::max({ std::abs(e0.bond - e_ref.bond), std::abs(e0.angle - e_ref.angle), std::abs(e0.dihedral - e_ref.dihedral),
                std::abs(e0.inversion - e_ref.inversion), std::abs(e0.hb - e_ref.hb), std::abs(e0.rep - e_ref.rep) });
            check(dq < 1e-10 && de < 1e-12, std::string(label) + ": charges/bonded/HB/rep with cell == without cell (" + sci(dq) + ", " + sci(de) + ")");
        }
    }

    struct V {
        const char* what;
        const Geometry* g;
        bool charges_equal;
    };
    const V variants[] = { { "anchored tear", &g_anchored, true }, { "fully torn", &g_torn, frags.size() == 1 } };
    for (const auto& v : variants) {
        Mol m = withCell(base, Lv);
        m.m_geometry = *v.g;
        int nmoved = 0;
        for (int i = 0; i < n; ++i)
            if ((v.g->row(i) - g_corner.row(i)).norm() > 1e-9)
                ++nmoved;
        Instance inst(m, config(1, extra));
        check(inst.ok, std::string(label) + " " + v.what + ": direct initialisation (" + inst.error + ")");
        if (!inst.ok)
            continue;
        const Topo t = inst.topology();
        std::string why;
        check(sameGraph(t, t_ref, why), std::string(label) + " " + v.what + ": topology == intact (" + why + ")");
        const double dq = maxAbsDiff(t.qtopo, t_ref.qtopo);
        if (v.charges_equal)
            check(dq < 1e-10, std::string(label) + " " + v.what + ": topology charges == intact (" + sci(dq) + ")");
        const Terms e = inst.evaluate(*v.g);
        const double de = maxTermDiff(e, e_ref);
        check(de < 1e-12, std::string(label) + " " + v.what + ": periodic terms == intact (" + sci(de) + " Eh)");
        std::printf("%-18s %-14s %3d/%3d atoms moved: %zu bonds, %zu/%zu/%zu ang/dih/inv, %d frag; |dq_topo| %.1e, |dE_periodic| %.1e, "
                    "E_coul %.9f vs %.9f (not periodic), E_hb %.3e\n",
            label, v.what, nmoved, n, t.bonds.size(), t.nangles, t.ndihedrals, t.ninversions, t.nfrag, dq, de, e.coulomb, e_ref.coulomb, e.hb);
    }

    {  // not vacuous: the torn geometry without a cell has a different topology
        Mol m0 = base;
        m0.m_geometry = g_torn;
        Instance nc(m0, config(1, extra));
        if (nc.ok) {
            const Topo t0 = nc.topology();
            check(t0.bonds != t_ref.bonds, std::string(label) + ": torn geometry without a cell differs (" + std::to_string(t0.bonds.size())
                    + " vs " + std::to_string(t_ref.bonds.size()) + " bonds; test not vacuous)");
        }
    }
}

// ---------------------------------------------------------------------------------------------
// (b) 2x1x1 supercell
static void testSupercell(const char* label, const Mol& base, const json& extra = json::object())
{
    const double L = kCellBohr / A2B;
    const int n = base.m_number_atoms;
    Mol cell = withCell(base, Eigen::Vector3d(L, L, L));
    cell.m_geometry = centred(base.m_geometry, 0.5 * Eigen::Vector3d(L, L, L));
    Instance one(cell, config(1, extra));
    check(one.ok, std::string(label) + ": single cell init");
    if (!one.ok)
        return;
    const Topo t1 = one.topology();
    const Terms e1 = one.evaluate(cell.m_geometry);

    Mol sup = withCell(base, Eigen::Vector3d(2.0 * L, L, L));
    sup.m_number_atoms = 2 * n;
    sup.m_charge = 2 * base.m_charge;
    sup.m_geometry = Geometry(2 * n, 3);
    sup.m_atoms.clear();
    for (int c = 0; c < 2; ++c)
        for (int i = 0; i < n; ++i) {
            sup.m_geometry.row(c * n + i) = cell.m_geometry.row(i) + Eigen::RowVector3d(c * L, 0.0, 0.0);
            sup.m_atoms.push_back(base.m_atoms[i]);
        }
    Instance two(sup, config(1, extra));
    check(two.ok, std::string(label) + ": supercell init (" + two.error + ")");
    if (!two.ok)
        return;
    const Topo t2 = two.topology();
    const Terms e2 = two.evaluate(sup.m_geometry);
    check(t2.bonds.size() == 2 * t1.bonds.size() && t2.nangles == 2 * t1.nangles && t2.ndihedrals == 2 * t1.ndihedrals
            && t2.ninversions == 2 * t1.ninversions,
        std::string(label) + ": supercell counts exactly doubled");
    const double de = maxTermDiff(e2, e1, 2.0);
    check(de < 2e-12, std::string(label) + ": supercell periodic terms == 2 x cell (" + sci(de) + " Eh)");
    std::printf("%-18s 2x1x1: bonds %zu/%zu, ang %zu/%zu, dih %zu/%zu, inv %zu/%zu, max |E2 - 2 E1| %.2e Eh (E_coul %.6f vs 2 x %.6f, not periodic)\n",
        label, t2.bonds.size(), t1.bonds.size(), t2.nangles, t1.nangles, t2.ndihedrals, t1.ndihedrals, t2.ninversions, t1.ninversions, de,
        e2.coulomb, e1.coulomb);
}

// ---------------------------------------------------------------------------------------------
// (c) random water box: bond detection against an explicit image enumeration
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
    while (placed < nwat && tries < 400000) {
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

/// Pass-2 bond criterion of gfnff_ini2.f90:111-126 / 361-419 (icase 1, non-metals) over all 27 images.
static std::vector<Pair> imageBonds(const Mol& m, double L_bohr, const Vector& qa)
{
    using GFNFFParameters::normcn;
    constexpr double rthr = 1.25, rqshrink = 0.23;
    const int n = m.m_number_atoms;
    std::vector<Pair> out;
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j) {
            const int zi = m.m_atoms[i], zj = m.m_atoms[j];
            double rco = GFNFFParameters::computeRabEstimate(zi, zj, normcn[zi - 1], normcn[zj - 1]);
            rco -= qa(i) * rqshrink + qa(j) * rqshrink;
            rco *= GFNFF::bondDetectionElementFactor(zi) * GFNFF::bondDetectionElementFactor(zj);
            const double thr = rthr * rco;
            bool bonded = false;
            for (int a = -1; a <= 1 && !bonded; ++a)
                for (int b = -1; b <= 1 && !bonded; ++b)
                    for (int c = -1; c <= 1 && !bonded; ++c) {
                        Eigen::RowVector3d d = (m.m_geometry.row(i) - m.m_geometry.row(j)) * A2B + Eigen::RowVector3d(a, b, c) * L_bohr;
                        if (d.norm() < thr)
                            bonded = true;
                    }
            if (bonded)
                out.emplace_back(i, j);
        }
    return out;
}

static void testWaterBox()
{
    const int nwat = 60;
    const double L = 43.5 / A2B;
    Mol box = waterBox(nwat, L, 20261101u);
    const int n = box.m_number_atoms;
    if (n != 3 * nwat)
        return;
    box.m_geometry = wrapped(box.m_geometry, L);
    int ntorn = 0;
    for (int w = 0; w < nwat; ++w) {
        bool torn = false;
        for (int h = 1; h < 3; ++h)
            if ((box.m_geometry.row(3 * w + h) - box.m_geometry.row(3 * w)).norm() > 2.0)
                torn = true;
        ntorn += torn;
    }
    std::vector<Pair> expected;
    for (int w = 0; w < nwat; ++w) {
        expected.emplace_back(3 * w, 3 * w + 1);
        expected.emplace_back(3 * w, 3 * w + 2);
    }
    Mol cell = withCell(box, Eigen::Vector3d(L, L, L));
    Instance t1(cell, config(1));
    check(t1.ok, "water box: init with cell (" + t1.error + ")");
    if (!t1.ok)
        return;
    const Topo a = t1.topology();
    check(a.bonds == expected, "water box: bonds are exactly the O-H pairs of each water (" + std::to_string(a.bonds.size()) + " vs 120)");
    const std::vector<Pair> img = imageBonds(box, L * A2B, a.qtopo);
    check(a.bonds == img, "water box: bonds == 27-image enumeration of the bond criterion (" + std::to_string(img.size()) + ")");
    check(a.nfrag == nwat && a.nangles == static_cast<size_t>(nwat), "water box: " + std::to_string(a.nfrag) + " fragments, " + std::to_string(a.nangles) + " angles");
    const Terms e1 = t1.evaluate(box.m_geometry);
    Instance t4(cell, config(4));
    check(t4.ok, "water box: init with 4 threads");
    if (t4.ok) {
        const Topo b = t4.topology();
        const Terms e4 = t4.evaluate(box.m_geometry);
        std::string why;
        check(sameGraph(a, b, why) && maxAbsDiff(a.qtopo, b.qtopo) < 1e-10, "water box: topology 1 == 4 threads (" + why + ")");
        check(maxTermDiff(e1, e4) < 1e-12, "water box: periodic terms 1 == 4 threads (" + sci(maxTermDiff(e1, e4)) + ")");
    }
    {  // not vacuous
        Instance nc(box, config(1));
        if (nc.ok) {
            const Topo c = nc.topology();
            check(ntorn == 0 || c.bonds.size() < expected.size(), "water box: without a cell the torn waters lose bonds (" + std::to_string(c.bonds.size()) + ")");
            std::printf("water box: %d waters, %d torn by the wrap, %zu bonds with cell (expected 120, 27-image reference %zu), %zu without a cell; "
                        "E_bond %.9f E_hb %.6e E_rep %.6e, 1/4 threads\n",
                nwat, ntorn, a.bonds.size(), img.size(), c.bonds.size(), e1.bond, e1.hb, e1.rep);
        }
    }
}

// ---------------------------------------------------------------------------------------------
// (d) reactive topology across the boundary
struct ReactTrace {
    int event_step = -1;
    double e_after = 0.0, de_jump = 0.0, e_last = 0.0;
    std::vector<Pair> formed;
    std::vector<Pair> bonds_after;
    bool ok = false;
};

static ReactTrace runRecombination(double L_bohr, bool across_face)
{
    Mol h2 = registryMol("H2");
    const double L = L_bohr / A2B;
    h2 = withCell(h2, Eigen::Vector3d(L, L, L));
    const double d0 = 4.0 / A2B, d1 = 1.2 / A2B;  // Angstrom; H-H from 4.0 to 1.2 Bohr
    auto place = [&](double d) {
        Geometry g(2, 3);
        if (across_face) {
            g.row(0) = Eigen::RowVector3d(0.3 / A2B, 0.5 * L, 0.5 * L);
            g.row(1) = Eigen::RowVector3d(L - (d - 0.3 / A2B), 0.5 * L, 0.5 * L);
        } else {
            g.row(0) = Eigen::RowVector3d(0.5 * L - 0.5 * d, 0.5 * L, 0.5 * L);
            g.row(1) = Eigen::RowVector3d(0.5 * L + 0.5 * d, 0.5 * L, 0.5 * L);
        }
        return g;
    };
    h2.m_geometry = place(d0);
    json extra = { { "topology_mode", "react" }, { "react_check_every", 1 }, { "hbond", false } };
    ReactTrace tr;
    Instance inst(h2, config(1, extra));
    if (!inst.ok) {
        std::printf("react init failed: %s\n", inst.error.c_str());
        return tr;
    }
    const int nsteps = 28;
    for (int s = 0; s <= nsteps; ++s) {
        const double d = d0 + (d1 - d0) * s / nsteps;
        inst.ff->UpdateMolecule(place(d));
        tr.e_last = inst.ff->Calculation(false);
        auto events = inst.ff->consumeReactEvents();
        for (const auto& ev : events)
            if (!ev.formed.empty() && tr.event_step < 0) {
                tr.event_step = s;
                tr.formed = ev.formed;
                tr.de_jump = ev.de_jump_eh;
                tr.e_after = tr.e_last;
            }
    }
    tr.bonds_after = inst.ff->reactiveBonds();
    tr.ok = true;
    return tr;
}

static void testReactive()
{
    const ReactTrace c = runRecombination(kCellBohr, false);
    const ReactTrace f = runRecombination(kCellBohr, true);
    check(c.ok && f.ok, "react: both runs initialised");
    if (!(c.ok && f.ok))
        return;
    check(c.event_step >= 0, "react: formation event in the centre (step " + std::to_string(c.event_step) + ")");
    check(f.event_step == c.event_step, "react: formation event across the face at the same step (" + std::to_string(f.event_step) + " vs "
            + std::to_string(c.event_step) + ")");
    check(f.formed == c.formed && f.bonds_after == c.bonds_after && c.bonds_after.size() == 1, "react: same bond formed across the face");
    check(std::abs(f.e_after - c.e_after) < 1e-12, "react: energy after the event equal across the face (" + sci(std::abs(f.e_after - c.e_after)) + ")");
    check(std::abs(f.e_last - c.e_last) < 1e-12, "react: final energy equal across the face (" + sci(std::abs(f.e_last - c.e_last)) + ")");
    check(std::isfinite(c.de_jump) && std::abs(f.de_jump - c.de_jump) < 1e-12, "react: dE_jump equal across the face (" + sci(std::abs(f.de_jump - c.de_jump)) + ")");
    std::printf("react H + H: centre event at step %d, across the face at step %d; E_after %.12f / %.12f, dE_jump %.6e / %.6e, E_final %.12f / %.12f\n",
        c.event_step, f.event_step, c.e_after, f.e_after, c.de_jump, f.de_jump, c.e_last, f.e_last);
}

// ---------------------------------------------------------------------------------------------
// (e) topology cache fingerprint carries the cell
static std::string fingerprintOf(const std::string& topo_file)
{
    std::ifstream in(topo_file);
    if (!in.good())
        return "";
    json j;
    in >> j;
    return j.value("fingerprint", std::string());
}

static void testCacheKey()
{
    const std::string base = "pbc_topology_cache_probe";
    const std::string xyz = base + ".xyz", topo = base + ".topo.json";
    std::remove(topo.c_str());
    Mol caf = registryMol("caffeine");
    const double L = kCellBohr / A2B;
    caf.m_geometry = centred(caf.m_geometry, 0.5 * Eigen::Vector3d(L, L, L));
    json c = config(1);
    c["cache_topology"] = true;
    c["geometry_file"] = xyz;

    Instance with_cell(withCell(caf, Eigen::Vector3d(L, L, L)), c);
    const std::string fp_cell = fingerprintOf(topo);
    Instance no_cell(caf, c);
    const std::string fp_plain = fingerprintOf(topo);
    Instance with_cell2(withCell(caf, Eigen::Vector3d(L, L, L)), c);
    const std::string fp_cell2 = fingerprintOf(topo);
    Instance other_cell(withCell(caf, Eigen::Vector3d(L, L, 1.1 * L)), c);
    const std::string fp_other = fingerprintOf(topo);
    std::remove(topo.c_str());
    check(with_cell.ok && no_cell.ok && with_cell2.ok && other_cell.ok, "cache: all instances initialised");
    check(!fp_cell.empty() && !fp_plain.empty(), "cache: topology cache files were written");
    check(fp_cell != fp_plain, "cache: fingerprint with a cell differs from the one without (same geometry and bond list)");
    check(fp_cell == fp_cell2, "cache: fingerprint reproduced by the same cell");
    check(fp_other != fp_cell && fp_other != fp_plain, "cache: another cell gives another fingerprint");
    std::printf("cache: fingerprints cell %s, plain %s, cell again %s, other cell %s\n", fp_cell.c_str(), fp_plain.c_str(), fp_cell2.c_str(), fp_other.c_str());
}

// ---------------------------------------------------------------------------------------------
// (f) smallest accepted cell with hbond false
static void testSmallestCell()
{
    const Mol dimer = registryMol("H2O_dimer");
    for (double Lb : { 40.5, 39.9 }) {
        const double L = Lb / A2B;
        Mol m = withCell(dimer, Eigen::Vector3d(L, L, L));
        m.m_geometry = centred(dimer.m_geometry, 0.5 * Eigen::Vector3d(L, L, L));
        Instance inst(m, config(1, { { "hbond", false } }));
        const bool expect_ok = Lb > 40.0;
        check(inst.ok == expect_ok, "smallest cell: L " + std::to_string(Lb) + " Bohr with hbond false " + (expect_ok ? "accepted" : "refused")
                + (inst.ok ? "" : " (" + inst.error + ")"));
    }
}

int main()
{
    CurcumaLogger::set_verbosity(0);
    testTornInit("caffeine", registryMol("caffeine"));
    testTornInit("hydronium-water2", registryMol("hydronium_water2", 1));
    testTornInit("acetic acid dimer", registryMol("acetic_acid_dimer"));
    testSupercell("caffeine", registryMol("caffeine"));
    testSupercell("hydronium-water2", registryMol("hydronium_water2", 1));
    testWaterBox();
    testReactive();
    testCacheKey();
    testSmallestCell();
    if (g_failures == 0)
        std::printf("pbc_topology: all checks passed\n");
    else
        std::printf("pbc_topology: %d check(s) failed\n", g_failures);
    return g_failures == 0 ? 0 : 1;
}
