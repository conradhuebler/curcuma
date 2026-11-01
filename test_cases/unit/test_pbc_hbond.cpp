/*
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 */
/**
 * Periodic GFN-FF hydrogen- and halogen-bond terms (WP-PERIODIC-NONBONDED A4b) - Claude Generated (Nov 2026).
 *
 * Checks ONLY what A4b makes periodic: the HB and XB kernels, their (re-)detection, the HB coordination
 * number that modulates the X-H bond term (read through the bond energy), and the RMSD trigger of the
 * re-detection. Read as single-term energies (GFNFF::HydrogenBondEnergy(), HalogenBondEnergy(),
 * BondEnergy()) and workspace gradient components. The HB/XB terms use the Phase-1 topology charges,
 * so Coulomb/EEQ (not periodic) enter no checked number.
 *
 * The topology is built from plain distances (A4c), so every image variant is handed to ONE instance
 * initialised with the intact system in the cell centre (UpdateMolecule, topology_mode constant);
 * hb_update_force_every 1 re-detects the HB/XB lists on every call, i.e. on the variant geometry.
 *
 * (a) hydronium-water2 (+1, real O-H...O), cyclic acetic acid dimer (O-H...O=C, case 3) and a composite of
 *     the dimer and F-I...pyridine (case 4 to the pyridine N, the XB F-I...N) in a 44 Bohr cube: centred,
 *     shifted by (L,-L,2L), each fragment by a different lattice vector (acceptor across the boundary),
 *     centred on the corner and wrapped per atom (torn). E_HB, E_XB, E_bond equal (< 1e-12 Eh), the HB,
 *     XB and bond gradients per atom (< 1e-12 Eh/Bohr), the triple lists identical (indices and case).
 *     HB case 4 (sp2 N with two neighbours) occurs in none of them (pyridine N is bonded to I in the topology,
 *     caffeine gives no case-4 triple either): not tested. Hydronium runs with frag_charge_model reference: the
 *     default ensemble model blends charge variants whose fragmentation uses plain distances (not periodic).
 * (b) XB: fi-pyridine (in the composite of (a)); the only library structure with a halogen bond.
 * (c) 40 random waters (L 43 Bohr, hb_min_pair_energy_eh 0): HB triple set and E_HB against an independent
 *     27-image enumeration (image of A and B next to H, chosen by explicit search over the images),
 *     threads 1 and 4; L 42.3 Bohr (sqrt(hbthr2) 21.2 > L/2) refused.
 * (d) 2x1x1 supercell of (c): E_HB exactly twice.
 * (e) torn composite: central finite difference (Angstrom, h 2e-4 and 1e-4 A) of E_HB and E_XB for the
 *     atoms nearest a face against the analytic HB/XB gradients; quadratic in h.
 * (f) re-detection trigger (hb_update_force_every 0): a wrap by +L does not re-detect and keeps E_HB; a
 *     move with RMSD > 0.3 Bohr re-detects and then equals an instance that re-detects every call.
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
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <memory>
#include <random>
#include <set>
#include <string>
#include <tuple>
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

static json config(int threads, const json& extra)
{
    json c = json::object();
    c["verbosity"] = 0;
    c["threads"] = threads;
    c["cache_topology"] = false;
    c["hb_update_force_every"] = 1;  // re-detect on every call: the lists come from the variant geometry
    for (auto it = extra.begin(); it != extra.end(); ++it)
        c[it.key()] = it.value();
    return c;
}

using Triple = std::tuple<int, int, int, int>;  // A, H, B, case

struct Result {
    double hb = 0, xb = 0, bond = 0;
    Matrix g_hb, g_xb, g_bond;
    std::vector<Triple> hbs, xbs;
};

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
    Result evaluate(const Geometry& g_ang, bool gradient)
    {
        ff->UpdateMolecule(g_ang);
        ff->Calculation(gradient);
        Result r;
        r.hb = ff->HydrogenBondEnergy();
        r.xb = ff->HalogenBondEnergy();
        r.bond = ff->BondEnergy();
        if (gradient) {
            const FFWorkspace* ws = ff->getWorkspace();
            r.g_hb = ws->gradientHB();
            r.g_xb = ws->gradientXB();
            r.g_bond = ws->gradientBond();
        }
        for (const auto& h : ff->getLastHBonds())
            r.hbs.emplace_back(h.i, h.j, h.k, h.case_type);
        for (const auto& x : ff->getLastXBonds())
            r.xbs.emplace_back(x.i, x.j, x.k, 0);
        std::sort(r.hbs.begin(), r.hbs.end());
        std::sort(r.xbs.begin(), r.xbs.end());
        return r;
    }
};

static double maxAbs(const Matrix& a, const Matrix& b)
{
    if (a.rows() != b.rows() || a.cols() != b.cols())
        return 1e300;
    return a.rows() ? (a - b).cwiseAbs().maxCoeff() : 0.0;
}

static Geometry centred(const Geometry& g0, const Eigen::Vector3d& at)
{
    Geometry g = g0;
    Eigen::RowVector3d c = g.colwise().mean();
    for (int i = 0; i < g.rows(); ++i)
        g.row(i) += at.transpose() - c;
    return g;
}

/// Acetic acid dimer plus F-I...pyridine 6.5 A away (rigid translation of the library geometries).
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

static std::vector<std::vector<int>> fragmentsOf(const Mol& m)
{
    curcuma::Molecule mol(m);
    return mol.GetFragments();
}

static const double kCellBohr = 44.0;

// ---------------------------------------------------------------------------------------------
// (a) image invariance, torn input, triple identity
static void testInvariance(const char* label, const Mol& base, const json& extra = json::object())
{
    const int n = base.m_number_atoms;
    const double L = kCellBohr / A2B;
    const Eigen::Vector3d Lv(L, L, L);
    const Geometry g_mid = centred(base.m_geometry, 0.5 * Lv);
    const Geometry g_corner = centred(base.m_geometry, Eigen::Vector3d::Zero());
    const auto frags = fragmentsOf(base);

    Geometry g_shift = g_mid;
    for (int i = 0; i < n; ++i) {
        g_shift(i, 0) += L;
        g_shift(i, 1) -= L;
        g_shift(i, 2) += 2.0 * L;
    }
    Geometry g_frag = g_mid;  // fragment f by a different lattice vector; f = 1 (acceptor side) by +L in x
    const int shifts[4][3] = { { 0, 0, 0 }, { 1, 0, 0 }, { -1, 1, 0 }, { 0, -1, 2 } };
    for (size_t f = 0; f < frags.size(); ++f)
        for (int a : frags[f])
            for (int k = 0; k < 3; ++k)
                g_frag(a, k) += shifts[f % 4][k] * L;
    Geometry g_torn = g_corner;
    for (int i = 0; i < n; ++i)
        for (int k = 0; k < 3; ++k)
            g_torn(i, k) -= L * std::floor(g_torn(i, k) / L);

    Mol m_mid = withCell(base, Lv);
    m_mid.m_geometry = g_mid;
    Instance inst(m_mid, config(1, extra));
    check(inst.ok, std::string(label) + ": init in the cell centre");
    if (!inst.ok)
        return;
    const Result ref = inst.evaluate(g_mid, true);
    check(!ref.hbs.empty() && std::abs(ref.hb) > 1e-6, std::string(label) + ": has hydrogen bonds (" + sci(ref.hb) + " Eh)");
    int ncase[5] = { 0, 0, 0, 0, 0 };
    for (const auto& t : ref.hbs)
        ++ncase[std::clamp(std::get<3>(t), 0, 4)];

    {  // no cell, same geometry (the cell is 44 Bohr, so contacts beyond L/2 = 22 Bohr may differ)
        Mol m0 = base;
        m0.m_geometry = g_mid;
        Instance nc(m0, config(1, extra));
        const Result r0 = nc.ok ? nc.evaluate(g_mid, true) : Result();
        const double d = std::max({ std::abs(r0.hb - ref.hb), std::abs(r0.xb - ref.xb), std::abs(r0.bond - ref.bond) });
        const double dg = std::max({ maxAbs(r0.g_hb, ref.g_hb), maxAbs(r0.g_xb, ref.g_xb), maxAbs(r0.g_bond, ref.g_bond) });
        // The XB list may legitimately differ: a system wider than L/2 (the composite, 18 A) sees the
        // periodic image of its own atoms within the 21.2 Bohr XB search radius (F...O image at
        // r_XB > 20 Bohr, energy 0 in the kernel). Energies and gradients must agree.
        if (r0.xbs != ref.xbs)
            std::printf("%-18s centred: XB list with a cell %zu triples, without %zu (own periodic image within the search radius)\n",
                label, ref.xbs.size(), r0.xbs.size());
        check(nc.ok && d < 1e-12 && dg < 1e-12 && r0.hbs == ref.hbs,
            std::string(label) + ": centred with cell == without cell (" + sci(d) + ", " + sci(dg) + ", triples "
                + std::to_string(r0.hbs.size()) + "/" + std::to_string(ref.hbs.size()) + " HB, " + std::to_string(r0.xbs.size()) + "/"
                + std::to_string(ref.xbs.size()) + " XB)");
    }

    struct V {
        const char* what;
        const Geometry* g;
    };
    const V variants[] = { { "shifted by (L,-L,2L)", &g_shift }, { "fragments by different lattice vectors", &g_frag },
        { "centred on the corner", &g_corner }, { "torn (per-atom wrap)", &g_torn } };
    for (const auto& v : variants) {
        const Result r = inst.evaluate(*v.g, true);
        const double d = std::max({ std::abs(r.hb - ref.hb), std::abs(r.xb - ref.xb), std::abs(r.bond - ref.bond) });
        const double dg = std::max({ maxAbs(r.g_hb, ref.g_hb), maxAbs(r.g_xb, ref.g_xb), maxAbs(r.g_bond, ref.g_bond) });
        check(d < 1e-12, std::string(label) + " " + v.what + ": E_HB, E_XB, E_bond == centred (" + sci(std::abs(r.hb - ref.hb)) + ", "
                + sci(std::abs(r.xb - ref.xb)) + ", " + sci(std::abs(r.bond - ref.bond)) + ")");
        check(dg < 1e-12, std::string(label) + " " + v.what + ": HB/XB/bond gradients == centred (" + sci(dg) + " Eh/Bohr)");
        check(r.hbs == ref.hbs, std::string(label) + " " + v.what + ": same HB triples (" + std::to_string(r.hbs.size()) + " vs "
                + std::to_string(ref.hbs.size()) + ")");
        check(r.xbs == ref.xbs, std::string(label) + " " + v.what + ": same XB triples");
        std::printf("%-18s %-40s max |dE| %.2e Eh, max |dG| %.2e Eh/Bohr, triples %zu HB / %zu XB %s\n", label, v.what, d, dg,
            r.hbs.size(), r.xbs.size(), (r.hbs == ref.hbs && r.xbs == ref.xbs) ? "identical" : "DIFFERENT");
    }
    std::printf("%-18s centred: E_HB %.12e E_XB %.12e E_bond %.12f; %zu fragments; HB cases 1/2/3/4: %d/%d/%d/%d, XB %zu\n", label,
        ref.hb, ref.xb, ref.bond, frags.size(), ncase[1], ncase[2], ncase[3], ncase[4], ref.xbs.size());

    {  // not vacuous: without a cell the fragment-shifted geometry differs
        Mol m0 = base;
        m0.m_geometry = g_mid;
        Instance nc(m0, config(1, extra));
        const Result r = nc.ok ? nc.evaluate(g_frag, false) : Result();
        const double d = std::abs(r.hb - ref.hb) + std::abs(r.xb - ref.xb);
        check(nc.ok && (frags.size() < 2 || d > 1e-6), std::string(label) + ": fragment-shifted geometry without a cell differs (" + sci(d) + ")");
        std::printf("%-18s without a cell, fragments shifted: |dE_HB| + |dE_XB| %.3e Eh (test not vacuous)\n", label, d);
    }
}

// ---------------------------------------------------------------------------------------------
// random water box (intact molecules, minimum-image O-O >= 2.9 A), as in pbc_bonded_repulsion
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

/// Image of atom a next to ref, by explicit search over the 27 neighbouring images (independent of
/// PBCUtils). Bohr.
static Eigen::Vector3d nearestImage(const Eigen::MatrixXd& g, int a, const Eigen::Vector3d& ref, const Eigen::Vector3d& L)
{
    Eigen::Vector3d best = g.row(a).transpose();
    double best_d = 1e300;
    for (int i = -1; i <= 1; ++i)
        for (int j = -1; j <= 1; ++j)
            for (int k = -1; k <= 1; ++k) {
                const Eigen::Vector3d x = g.row(a).transpose() + Eigen::Vector3d(i * L[0], j * L[1], k * L[2]);
                const double d = (x - ref).squaredNorm();
                if (d < best_d) {
                    best_d = d;
                    best = x;
                }
            }
    return best;
}

/// HB energy of one water triple (case 1 or 2) in the frame of explicit images; formula of
/// gfnff_engrad.F90 abhgfnff_eg1 / eg2new (energy only), parameters from the curcuma list entry.
static double tripleEnergy(const GFNFFHydrogenBond& hb, const Eigen::Vector3d& pA, const Eigen::Vector3d& pH,
    const Eigen::Vector3d& pB, const std::vector<Eigen::Vector3d>& nbB, int zA, int zB)
{
    using namespace GFNFFParameters;
    const double r_AH = (pH - pA).norm(), r_HB = (pB - pH).norm(), r_AB = (pB - pA).norm();
    const double r_AH_4 = std::pow(r_AH, 4), r_HB_4 = std::pow(r_HB, 4);
    const double Q_H = ws_charge_scaling(hb.q_H, HB_ST, HB_SF);
    const double Q_A = ws_charge_scaling(-hb.q_A, HB_ST, HB_SF);
    const double Q_B = ws_charge_scaling(-hb.q_B, HB_ST, HB_SF);
    const double r_vdw = covalent_radii[zA - 1] + covalent_radii[zB - 1];
    const double damp_env = ws_damping_short_range(r_AB, r_vdw, HB_SCUT, HB_ALP) * ws_damping_long_range(r_AB, HB_LONGCUT, HB_ALP);
    const double damp_outl = ws_damping_out_of_line(r_AH, r_HB, r_AB, r_vdw, HB_BACUT);
    if (r_AB > hb.r_cut)
        return 0.0;
    if (hb.case_type == 1) {
        const double bas = (Q_A * hb.basicity_A * r_AH_4 + Q_B * hb.basicity_B * r_HB_4) / (r_AH_4 + r_HB_4);
        const double aci = (hb.acidity_B * r_AH_4 + hb.acidity_A * r_HB_4) / (r_AH_4 + r_HB_4);
        return -bas * aci * damp_env / std::pow(r_AB, 3) * Q_H * damp_outl;
    }
    const double nbcut = (zB == 7 && nbB.size() == 1) ? 2.0 : HB_NBCUT;
    double outl_nb = 1.0;
    for (const auto& x : nbB)
        outl_nb *= 2.0 / (1.0 + std::exp(-(nbcut / r_vdw) * (((x - pA).norm() + (x - pB).norm()) / r_AB - 1.0))) - 1.0;
    const double rdamp = damp_env * (1.8 / std::pow(r_HB, 3) - 0.8 / std::pow(r_AB, 3));
    return -rdamp * Q_H * damp_outl * outl_nb * hb.acidity_A * hb.basicity_B * Q_A * Q_B * XHACI_GLOBABH;
}

/// Case-1 triples are stored as (i, H, j) with the AB-pair order of the detection branch (legacy j < i,
/// cell list i < j); normalised to A > B so that two branches compare as sets.
static std::set<Triple> normalisedSet(const std::vector<Triple>& v)
{
    std::set<Triple> out;
    for (auto t : v) {
        if (std::get<3>(t) == 1 && std::get<0>(t) < std::get<2>(t))
            t = { std::get<2>(t), std::get<1>(t), std::get<0>(t), 1 };
        out.insert(t);
    }
    return out;
}

static const double kBoxBohr = 43.0;
static const int kWaters = 40;

// (c) brute force and (d) supercell
static void testBox()
{
    const double L = kBoxBohr / A2B;
    const Eigen::Vector3d Lv(L, L, L), Lb = Lv * A2B;
    const Mol box = withCell(waterBox(kWaters, L, 20261105), Lv);
    const int n = box.m_number_atoms;
    json nopr = json::object();
    nopr["hb_min_pair_energy_eh"] = 0.0;  // keep every case-1 triple, so the set is purely geometric
    Instance i1(box, config(1, nopr));
    check(i1.ok, "box: init");
    if (!i1.ok)
        return;
    const Result r1 = i1.evaluate(box.m_geometry, false);
    const auto& list = i1.ff->getLastHBonds();

    // independent reference: triple set from an explicit image search (the hydrogen fixes the frame)
    const Eigen::MatrixXd g = box.m_geometry * A2B;
    const double hbthr1 = 250.0, hbthr2 = 450.0;  // 200/400 - 50 log10(0.1), the defaults
    std::set<Triple> ref_set;
    for (int oi = 0; oi < n; oi += 3)
        for (int oj = 0; oj < n; oj += 3) {
            if (oi == oj)
                continue;
            // A-B distance: nearest image of oj next to oi
            const Eigen::Vector3d pi = g.row(oi).transpose();
            const double rAB2 = (nearestImage(g, oj, pi, Lb) - pi).squaredNorm();
            if (rAB2 > hbthr1)
                continue;
            for (int h = 0; h < n; ++h) {
                if (box.m_atoms[h] != 1)
                    continue;
                const Eigen::Vector3d pH = g.row(h).transpose();
                if (h == oi + 1 || h == oi + 2) {  // case 2, donor oi
                    ref_set.insert({ oi, h, oj, 2 });
                } else if (h != oj + 1 && h != oj + 2 && oi > oj) {  // case 1, each O-O pair once (A > B as curcuma)
                    const Eigen::Vector3d a = nearestImage(g, oi, pH, Lb), b = nearestImage(g, oj, pH, Lb);
                    if ((b - a).squaredNorm() + (a - pH).squaredNorm() + (b - pH).squaredNorm() < hbthr2)
                        ref_set.insert({ oi, h, oj, 1 });
                }
            }
        }
    const std::set<Triple> cur_norm = normalisedSet(r1.hbs);
    check(cur_norm == ref_set, "box: HB triple set == 27-image enumeration (" + std::to_string(cur_norm.size()) + " vs "
            + std::to_string(ref_set.size()) + ")");

    double e_ref = 0.0, e_central = 0.0;
    int n_imaged = 0;
    for (const auto& hb : list) {
        const Eigen::Vector3d pH = g.row(hb.j).transpose();
        const Eigen::Vector3d pA = nearestImage(g, hb.i, pH, Lb), pB = nearestImage(g, hb.k, pH, Lb);
        std::vector<Eigen::Vector3d> nb, nb0;
        for (int x : hb.neighbors_B) {
            nb.push_back(nearestImage(g, x, pB, Lb));
            nb0.push_back(g.row(x).transpose());
        }
        const double e = tripleEnergy(hb, pA, pH, pB, nb, box.m_atoms[hb.i], box.m_atoms[hb.k]);
        e_ref += e;
        const bool central = (pA - g.row(hb.i).transpose()).norm() < 1e-9 && (pB - g.row(hb.k).transpose()).norm() < 1e-9;
        if (central)
            e_central += e;
        else
            ++n_imaged;
    }
    const double d = std::abs(r1.hb - e_ref);
    check(d < 1e-12, "box: E_HB == 27-image reference (" + sci(d) + ")");
    check(n_imaged > 0 && std::abs(e_ref - e_central) > 1e-6, "box: image triples contribute (test not vacuous)");

    Instance i4(box, config(4, nopr));
    const Result r4 = i4.ok ? i4.evaluate(box.m_geometry, false) : Result();
    check(i4.ok && std::abs(r4.hb - r1.hb) < 1e-14 && r4.hbs == r1.hbs, "box: 1 and 4 threads agree (" + sci(std::abs(r4.hb - r1.hb)) + ")");

    json cl = nopr;
    cl["nb_cell_list_min_atoms"] = 0;  // periodic cell-list branch of the detection
    Instance ic(box, config(1, cl));
    const Result rc = ic.ok ? ic.evaluate(box.m_geometry, false) : Result();
    check(ic.ok && std::abs(rc.hb - r1.hb) < 1e-12 && normalisedSet(rc.hbs) == cur_norm, "box: cell-list detection == legacy detection (" + sci(std::abs(rc.hb - r1.hb)) + ")");
    std::printf("box %d atoms, L %.1f Bohr: E_HB %.15f, reference %.15f, |d| %.2e; triples %zu (ref %zu), %d through images "
                "(%.3e Eh); 1-4 threads %.1e, cell list %.1e\n",
        n, kBoxBohr, r1.hb, e_ref, d, cur_norm.size(), ref_set.size(), n_imaged, e_ref - e_central, std::abs(r4.hb - r1.hb),
        std::abs(rc.hb - r1.hb));

    {  // largest triple energy with r_HB near L/2: size of a possible image switch between re-detections
        double emax = 0.0, rmax = 0.0;
        const auto energies = i1.ff->hydrogenBondEnergies();
        for (const auto& hb : list) {
            const Eigen::Vector3d pH = g.row(hb.j).transpose();
            const double r = (nearestImage(g, hb.k, pH, Lb) - pH).norm();
            rmax = std::max(rmax, r);
            if (r > 18.0)
                emax = std::max(emax, std::abs(tripleEnergy(hb, nearestImage(g, hb.i, pH, Lb), pH, nearestImage(g, hb.k, pH, Lb),
                                                     std::vector<Eigen::Vector3d>(), box.m_atoms[hb.i], box.m_atoms[hb.k])));
        }
        (void)energies;
        std::printf("box: largest r_HB %.2f Bohr (L/2 %.2f); max |E_triple| with r_HB > 18 Bohr (neighbour factor 1) %.2e Eh\n",
            rmax, 0.5 * kBoxBohr, emax);
    }

    {  // smallest cube: sqrt(hbthr2) = 21.2 Bohr must be below L/2
        const Mol small = withCell(box, Eigen::Vector3d::Constant(42.3 / A2B));
        Instance refused(small, config(1, nopr));
        check(!refused.ok, "cell with L/2 < sqrt(hbthr2) is refused");
        json off = nopr;
        off["hbond"] = false;
        Instance accepted(small, config(1, off));
        check(accepted.ok, "the same cell with hbond false is accepted");
    }

    // (d) 2x1x1 supercell
    {
        Mol sc = box;
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
        Instance s(sc, config(1, nopr));
        check(s.ok, "supercell runs");
        if (s.ok) {
            const Result rs = s.evaluate(sc.m_geometry, false);
            const double dd = std::max(std::abs(rs.hb - 2.0 * r1.hb), std::abs(rs.bond - 2.0 * r1.bond));
            check(dd < 1e-12 && rs.hbs.size() == 2 * r1.hbs.size(), "water supercell 2x1x1: E_HB, E_bond == 2 x cell (" + sci(dd) + ")");
            std::printf("water supercell 2x1x1: E_HB %.15f, 2 x cell %.15f, |d| (HB, bond) %.2e, triples %zu = 2 x %zu\n", rs.hb,
                2.0 * r1.hb, dd, rs.hbs.size(), r1.hbs.size());
        }
    }
}

// ---------------------------------------------------------------------------------------------
// (e) finite difference on the torn composite
static void testFiniteDifference()
{
    const Mol base = compositeMol();
    const int n = base.m_number_atoms;
    const double L = kCellBohr / A2B;
    Mol m = withCell(base, Eigen::Vector3d(L, L, L));
    m.m_geometry = centred(base.m_geometry, Eigen::Vector3d::Zero());
    Geometry g_torn = m.m_geometry;
    for (int i = 0; i < n; ++i)
        for (int k = 0; k < 3; ++k)
            g_torn(i, k) -= L * std::floor(g_torn(i, k) / L);
    json fixed = json::object();
    fixed["hb_update_force_every"] = 0;  // fixed lists (no re-detection between the displaced points)
    Instance inst(m, config(1, fixed));
    check(inst.ok, "FD: init");
    if (!inst.ok)
        return;
    const Result r = inst.evaluate(g_torn, true);
    // the same molecule intact and without a cell: reference for the XB residual (see below)
    Mol m0 = base;
    m0.m_geometry = m.m_geometry;
    Instance nc(m0, config(1, fixed));
    check(nc.ok, "FD: init without a cell");
    if (!nc.ok)
        return;
    const Result r0 = nc.evaluate(m.m_geometry, true);
    double xb_res_diff = 0.0, xb_res_max = 0.0;
    std::vector<std::pair<double, int>> dist;
    for (int i = 0; i < n; ++i) {
        double dmin = 1e300;
        for (int k = 0; k < 3; ++k)
            dmin = std::min({ dmin, g_torn(i, k), L - g_torn(i, k) });
        dist.push_back({ dmin, i });
    }
    std::sort(dist.begin(), dist.end());
    const int nfd = 10;
    const char* names[2] = { "HB", "XB" };
    const Matrix* an[2] = { &r.g_hb, &r.g_xb };
    double err[2][2] = { { 0, 0 }, { 0, 0 } }, gmax[2] = { 0, 0 };
    const double hs[2] = { 2e-4, 1e-4 };
    for (int hi = 0; hi < 2; ++hi)
        for (int s = 0; s < nfd; ++s) {
            const int a = dist[s].second;
            for (int k = 0; k < 3; ++k) {
                Geometry gp = g_torn, gm = g_torn;
                gp(a, k) += hs[hi];
                gm(a, k) -= hs[hi];
                const Result rp = inst.evaluate(gp, false), rm = inst.evaluate(gm, false);
                const double fd[2] = { (rp.hb - rm.hb) / (2.0 * hs[hi]), (rp.xb - rm.xb) / (2.0 * hs[hi]) };
                for (int t = 0; t < 2; ++t) {
                    const double ana = (*an[t])(a, k) * A2B;
                    err[hi][t] = std::max(err[hi][t], std::abs(fd[t] - ana));
                    gmax[t] = std::max(gmax[t], std::abs(ana));
                }
                Geometry hp = m.m_geometry, hm = m.m_geometry;
                hp(a, k) += hs[hi];
                hm(a, k) -= hs[hi];
                const double fd0 = (nc.evaluate(hp, false).xb - nc.evaluate(hm, false).xb) / (2.0 * hs[hi]);
                const double res = fd[1] - (*an[1])(a, k) * A2B, res0 = fd0 - r0.g_xb(a, k) * A2B;
                xb_res_diff = std::max(xb_res_diff, std::abs(res - res0));
                xb_res_max = std::max(xb_res_max, std::abs(res0));
            }
        }
    // HB: analytic against FD. XB: the analytic XB gradient has a defect that does not depend on the
    // cell (FFWorkspace::calcHalogenBonds, damping derivatives: the short-range sign and r/c instead of
    // 1/r in the long-range one; residual 1.2e-4 of 5.4e-4 Eh/A without a cell as well, TODO.md). So for
    // XB the check is that the cell changes nothing: the FD residual of the torn system equals the
    // residual of the intact system without a cell, per atom and component.
    check(xb_res_diff < 1e-9, "FD XB: residual torn with cell == residual intact without cell (" + sci(xb_res_diff) + ")");
    std::printf("FD XB: max residual without a cell %.3e Eh/A (pre-existing XB gradient defect), torn - intact %.2e\n", xb_res_max,
        xb_res_diff);
    for (int t = 0; t < 2; ++t) {
        if (t == 1) {
            std::printf("FD %s (%d atoms %.2f..%.2f A from a face): max |g| %.3e Eh/A, |an - FD| %.2e (h 2e-4 A), %.2e (h 1e-4 A), not checked\n",
                names[t], nfd, dist[0].first, dist[nfd - 1].first, gmax[t], err[0][t], err[1][t]);
            continue;
        }
        check(err[1][t] < 1e-8, std::string("FD ") + names[t] + ": max error " + sci(err[1][t]) + " Eh/A");
        check(err[1][t] < err[0][t] / 3.0 || err[1][t] < 1e-11,
            std::string("FD ") + names[t] + ": error falls quadratically (" + sci(err[0][t]) + " -> " + sci(err[1][t]) + ")");
        check(gmax[t] > 1e-6, std::string("FD ") + names[t] + ": gradient non-zero on the boundary atoms (test not vacuous)");
        std::printf("FD %s (%d atoms %.2f..%.2f A from a face): max |g| %.3e Eh/A, |an - FD| %.2e (h 2e-4 A), %.2e (h 1e-4 A)\n",
            names[t], nfd, dist[0].first, dist[nfd - 1].first, gmax[t], err[0][t], err[1][t]);
    }
}

// (f) re-detection trigger
static void testRedetection()
{
    const Mol base = registryMol("hydronium_water2", 1);
    const double L = kCellBohr / A2B;
    Mol m = withCell(base, Eigen::Vector3d(L, L, L));
    m.m_geometry = centred(base.m_geometry, 0.5 * Eigen::Vector3d(L, L, L));
    json rmsd_only = json::object();
    rmsd_only["hb_update_force_every"] = 0;
    rmsd_only["frag_charge_model"] = "reference";
    Instance inst(m, config(1, rmsd_only));
    check(inst.ok, "re-detection: init");
    if (!inst.ok)
        return;
    Geometry g = m.m_geometry;
    const double e0 = inst.evaluate(g, true).hb;
    inst.ff->consumeHBXBUpdate();
    for (int i = 0; i < 3; ++i)
        g(i, 0) += L;  // first water by +L
    const double e1 = inst.evaluate(g, true).hb;
    const bool upd1 = inst.ff->consumeHBXBUpdate();
    check(!upd1, "wrap by +L: no HB/XB re-detection");
    check(std::abs(e1 - e0) < 1e-12, "wrap by +L: E_HB unchanged (" + sci(std::abs(e1 - e0)) + ")");
    // RMSD = sqrt(sum d^2) / N (reference formula, gfnff_ini2.f90:717): 3 atoms by 1.0 A = 1.89 Bohr,
    // sqrt(3) * 1.89 / 10 = 0.33 > 0.3
    for (int i = 0; i < 3; ++i)
        g(i, 1) += 1.0;
    const double e2 = inst.evaluate(g, true).hb;
    const bool upd2 = inst.ff->consumeHBXBUpdate();
    check(upd2, "move with RMSD > 0.3 Bohr re-detects");
    json every_cfg = json::object();
    every_cfg["frag_charge_model"] = "reference";
    Instance every(m, config(1, every_cfg));
    const double ef = every.ok ? every.evaluate(g, false).hb : 0.0;
    check(every.ok && std::abs(e2 - ef) < 1e-12, "after the re-detection E_HB == instance re-detecting every call (" + sci(std::abs(e2 - ef)) + ")");
    std::printf("re-detection: wrap |dE_HB| %.2e, re-detected %s; RMSD 0.33 Bohr move re-detected %s, |E - every-call| %.2e\n",
        std::abs(e1 - e0), upd1 ? "yes" : "no", upd2 ? "yes" : "no", std::abs(e2 - ef));
}

// Diagnostic (PBC_HB_DEBUG=1): the same FD without a cell on the intact composite (is an XB gradient error periodic?)
static void fdNoCell()
{
    Mol m = compositeMol();
    json fixed = json::object();
    fixed["hb_update_force_every"] = 0;
    Instance inst(m, config(1, fixed));
    if (!inst.ok)
        return;
    const Result r = inst.evaluate(m.m_geometry, true);
    double err[2] = { 0, 0 }, gmax[2] = { 0, 0 };
    int worst = -1;
    for (int a = 0; a < m.m_number_atoms; ++a)
        for (int k = 0; k < 3; ++k) {
            Geometry gp = m.m_geometry, gm = m.m_geometry;
            gp(a, k) += 1e-4;
            gm(a, k) -= 1e-4;
            const Result rp = inst.evaluate(gp, false), rm = inst.evaluate(gm, false);
            const double fd[2] = { (rp.hb - rm.hb) / 2e-4, (rp.xb - rm.xb) / 2e-4 };
            const double an[2] = { r.g_hb(a, k) * A2B, r.g_xb(a, k) * A2B };
            for (int t = 0; t < 2; ++t) {
                if (std::abs(fd[t] - an[t]) > err[t] && t == 1) worst = a;
                err[t] = std::max(err[t], std::abs(fd[t] - an[t]));
                gmax[t] = std::max(gmax[t], std::abs(an[t]));
            }
        }
    std::printf("FD without a cell (intact composite, all atoms, h 1e-4 A): HB |an - FD| %.2e (max |g| %.2e), XB %.2e (max |g| %.2e, worst atom %d)\n",
        err[0], gmax[0], err[1], gmax[1], worst);
}

int main()
{
    CurcumaLogger::set_verbosity(0);
    // net charge with several fragments: frag_charge_model ensemble blends charge variants whose
    // fragmentation and window use plain distances (not periodic, Block B); reference = master charges
    json refq = json::object();
    refq["frag_charge_model"] = "reference";
    testInvariance("hydronium-water2", registryMol("hydronium_water2", 1), refq);
    testInvariance("acetic-acid-dimer", registryMol("acetic_acid_dimer"));
    testInvariance("aa-dimer+fi-pyr", compositeMol());
    testBox();
    testFiniteDifference();
    if (std::getenv("PBC_HB_DEBUG"))
        fdNoCell();
    testRedetection();
    if (g_failures == 0)
        std::printf("pbc_hbond: all checks passed\n");
    return g_failures == 0 ? 0 : 1;
}
