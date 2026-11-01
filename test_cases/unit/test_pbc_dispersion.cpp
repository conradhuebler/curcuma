/*
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 */
/**
 * Periodic GFN-FF D4 dispersion pair term (WP-PERIODIC-NONBONDED A3b) - Claude Generated (Nov 2026).
 *
 * Checks ONLY the D4 dispersion term (GFNFF::DispersionEnergy(), the workspace dispersion energy, and
 * its gradient taken as G(dispersion on) - G(dispersion off) at the same geometry, which contains the
 * pair force and the dC6/dCN chain rule). Repulsion, Coulomb/EEQ, HB/XB and the bonded terms are not
 * periodic yet and are not tested. Topology (bonds) is not periodic either, so every molecule is
 * kept intact: molecules are moved by lattice vectors, never torn.
 *
 * (a) Water dimer in a 130 Bohr cell: E_disp with the second molecule shifted by L, and with the
 *     dimer contact across the y boundary, equals the non-periodic dimer (< 1e-12 Eh).
 * (b) 21 randomly placed and oriented waters (63 atoms) in a 24 Bohr cube, r_c = 11.5 Bohr
 *     (hbond and repulsion off so the A1 cell check passes): E_disp equals a brute-force sum over
 *     the 27 neighbouring images, using the per-pair parameters of the workspace list; every pair
 *     with an image within r_c must be in the list. 1 and 4 threads agree.
 * (c) 2x1x1 supercell of (b), same r_c: E_disp exactly twice the single cell (< 1e-12 Eh).
 * (d) Central finite difference of E_disp (Angstrom displacement, h = 2e-4 and 1e-4 A) for the atoms
 *     next to a cell face against the analytic dispersion gradient, error quadratic in h; the
 *     dE/dCN part (analytic minus the pair-only workspace gradient) is reported and must be non-zero.
 * (e) Cut-off reduction: default cut-off in the 24 Bohr cell is reduced to just below L/2 and the
 *     energy equals the brute-force sum at that radius; an explicit 12.5 or 50 Bohr is refused.
 * (f) Verlet list: wrapping one molecule by +L does not rebuild the D4 list and leaves E_disp
 *     unchanged; small random moves below the half-skin keep E_disp equal to a fresh instance; a
 *     move beyond the half-skin rebuilds.
 *
 * Geometry: water-dimer.roo2.8 from the structure library (registry "H2O_dimer"), only translated and
 * rotated rigidly. Exit 0 when every check passes; each failure prints one line.
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

static Mol dimerAngstrom()
{
    return TestMolecules::TestMoleculeRegistry::createMolecule("H2O_dimer", false).getMolInfo();
}

// ---------------------------------------------------------------------------------------------
struct Run {
    bool ok = false;
    bool threw = false;
    double e_disp = 0.0;
    Matrix grad;        // full gradient (Eh/Bohr, GFNFF::Gradient)
    Matrix grad_ws;     // workspace dispersion gradient (pair + dC6/dCN chain rule, Eh/Bohr)
    std::vector<GFNFFDispersion> pairs;
    GFNFF::D4RadiiBohr radii;
};

static json config(int threads, const json& extra)
{
    json c = json::object();
    c["verbosity"] = 0;
    c["threads"] = threads;
    c["cache_topology"] = false;
    for (auto it = extra.begin(); it != extra.end(); ++it)
        c[it.key()] = it.value();
    return c;
}

static Run runGFNFF(Mol mol, double L_ang, const json& extra, int threads = 1, bool dispersion = true)
{
    Run r;
    mol.m_has_pbc = (L_ang > 0.0);
    if (L_ang > 0.0)
        mol.m_unit_cell = Eigen::Matrix3d::Identity() * L_ang;
    json c = config(threads, extra);
    c["dispersion"] = dispersion;
    try {
        GFNFF ff(c);
        if (!ff.InitialiseMolecule(mol))
            return r;
        if (ff.getWorkspace())
            ff.getWorkspace()->setStoreGradientComponents(true);
        ff.Calculation(true);
        r.e_disp = ff.DispersionEnergy();
        r.grad = ff.Gradient();
        if (ff.getWorkspace()) {
            r.grad_ws = ff.getWorkspace()->gradientDispersion();
            r.pairs = ff.getWorkspace()->d4Dispersions();
        }
        r.radii = ff.dispersionRadiiBohr();
        r.ok = true;
    } catch (const std::exception&) {
        r.threw = true;
    }
    return r;
}

// ---------------------------------------------------------------------------------------------
// (a) water dimer, shift invariance
static void testDimer()
{
    const Mol base = dimerAngstrom();
    const int n = base.m_number_atoms;
    const double L = 130.0 / A2B;  // 130 Bohr: every default cutoff below L/2
    Mol centred = base;
    for (int i = 0; i < n; ++i)
        for (int k = 0; k < 3; ++k)
            centred.m_geometry(i, k) += 0.5 * L;
    const Run ref = runGFNFF(centred, 0.0, json::object());
    check(ref.ok && ref.e_disp < 0.0, "dimer without cell runs");

    // atoms 0-2 first water, 3-5 second water (water-dimer.roo2.8)
    Mol shifted = centred;
    for (int i = 3; i < 6; ++i)
        shifted.m_geometry(i, 2) += L;
    // contact across y = 0: the plane through the midpoint of the two oxygens, then the molecule
    // below it moved by +L (intact molecules, the O...O contact crosses the boundary)
    Mol across = centred;
    {
        const double ymid = 0.5 * (centred.m_geometry(0, 1) + centred.m_geometry(3, 1));
        for (int i = 0; i < n; ++i)
            across.m_geometry(i, 1) -= ymid;
        const int low = across.m_geometry(0, 1) < across.m_geometry(3, 1) ? 0 : 3;
        for (int i = low; i < low + 3; ++i)
            across.m_geometry(i, 1) += L;
    }
    const std::pair<const char*, Mol*> variants[] = { { "in the cell", &centred }, { "second molecule shifted by L", &shifted },
        { "contact across y = 0", &across } };
    for (const auto& v : variants) {
        const Run per = runGFNFF(*v.second, L, json::object());
        check(per.ok, std::string("dimer with cell runs: ") + v.first);
        const double d = std::abs(per.e_disp - ref.e_disp);
        check(d < 1e-12, std::string("dimer ") + v.first + ": E_disp == non-periodic (" + sci(d) + ")");
        std::printf("dimer %-30s E_disp %.15f, |dE| vs no cell %.2e\n", v.first, per.e_disp, d);
    }
    const Run flat = runGFNFF(across, 0.0, json::object());
    check(flat.ok && std::abs(flat.e_disp - ref.e_disp) > 1e-5, "contact across y=0 without a cell differs (test not vacuous)");
    std::printf("contact across y=0 without cell: |dE| %.3e Eh\n", std::abs(flat.e_disp - ref.e_disp));
}

// ---------------------------------------------------------------------------------------------
// random water box
static Mol waterBox(int nwat, double L_ang, unsigned seed)
{
    const Mol d = dimerAngstrom();
    Eigen::MatrixXd w(3, 3);
    for (int i = 0; i < 3; ++i)
        w.row(i) = d.m_geometry.row(i) - d.m_geometry.row(0);  // O at the origin
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    const double min_oo = 3.0;  // Angstrom, minimum image
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
            if (dv.norm() < min_oo) {
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
            const Eigen::Vector3d x = p + R * w.row(i).transpose();
            m.m_geometry.row(3 * placed + i) = x.transpose();
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

/// Brute-force D4 energy over the 27 neighbouring images with radius rc (Bohr), per-pair data from
/// the workspace list. Counts the images found within rc for pairs absent from the list.
static double bruteForce(const Mol& m, double L_bohr, double rc, const std::vector<GFNFFDispersion>& pairs,
    int& missing, int& self_images, double* central_only = nullptr)
{
    std::map<std::pair<int, int>, const GFNFFDispersion*> lut;
    for (const auto& p : pairs)
        lut[{ std::min(p.i, p.j), std::max(p.i, p.j) }] = &p;
    const Eigen::MatrixXd g = m.m_geometry * A2B;
    const int n = m.m_number_atoms;
    double e = 0.0, e0 = 0.0;
    missing = 0;
    self_images = 0;
    for (int i = 0; i < n; ++i)
        for (int j = i; j < n; ++j)
            for (int a = -1; a <= 1; ++a)
                for (int b = -1; b <= 1; ++b)
                    for (int c = -1; c <= 1; ++c) {
                        if (i == j && a == 0 && b == 0 && c == 0)
                            continue;
                        const double dx = g(i, 0) - g(j, 0) - a * L_bohr;
                        const double dy = g(i, 1) - g(j, 1) - b * L_bohr;
                        const double dz = g(i, 2) - g(j, 2) - c * L_bohr;
                        const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
                        if (r > rc)
                            continue;
                        if (i == j) {
                            ++self_images;
                            continue;
                        }
                        auto it = lut.find({ i, j });
                        if (it == lut.end()) {
                            ++missing;
                            continue;
                        }
                        const GFNFFDispersion& d = *it->second;
                        const double r2 = r * r, r6 = r2 * r2 * r2;
                        const double r06 = d.r0_squared * d.r0_squared * d.r0_squared;
                        const double t6 = 1.0 / (r6 + r06);
                        const double t8 = 1.0 / (r6 * r2 + r06 * d.r0_squared);
                        const double ep = -d.C6 * d.zetac6 * (t6 + 2.0 * d.r4r2ij * t8);
                        e += ep;
                        if (a == 0 && b == 0 && c == 0)
                            e0 += ep;
                    }
    if (central_only)
        *central_only = e0;
    return e;
}

/// Pair-only analytic gradient (C6 held fixed, minimum image, rc < L/2), Eh/Bohr: the dispersion
/// gradient minus this is the dC6/dCN chain-rule part.
static Matrix pairGradient(const Mol& m, double L_bohr, double rc, const std::vector<GFNFFDispersion>& pairs)
{
    const Eigen::MatrixXd g = m.m_geometry * A2B;
    Matrix out = Matrix::Zero(m.m_number_atoms, 3);
    for (const auto& d : pairs) {
        double v[3];
        for (int k = 0; k < 3; ++k) {
            v[k] = g(d.i, k) - g(d.j, k);
            v[k] -= L_bohr * std::rint(v[k] / L_bohr);
        }
        const double r2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
        const double r = std::sqrt(r2);
        if (r > rc)
            continue;
        const double r06 = d.r0_squared * d.r0_squared * d.r0_squared;
        const double t6 = 1.0 / (r2 * r2 * r2 + r06), t8 = 1.0 / (r2 * r2 * r2 * r2 + r06 * d.r0_squared);
        // dE/dr = -C6 zeta [ -6 r^5 t6^2 - 2 r4r2 8 r^7 t8^2 ]
        const double dEdr = -d.C6 * d.zetac6 * (-6.0 * r2 * r2 * r * t6 * t6 - 16.0 * d.r4r2ij * r2 * r2 * r2 * r * t8 * t8);
        for (int k = 0; k < 3; ++k) {
            out(d.i, k) += dEdr * v[k] / r;
            out(d.j, k) -= dEdr * v[k] / r;
        }
    }
    return out;
}

static json boxFlags(double rc)
{
    json e = json::object();
    e["hbond"] = false;      // HB/XB search radii 15.8/21.2 Bohr exceed L/2 = 12 Bohr
    e["repulsion"] = false;  // repulsion list 22 Bohr likewise
    e["dispersion_tail"] = false;  // Nov 2026 (A3c): pair term only here; the tail is tested in pbc_dispersion_tail
    if (rc > 0.0)
        e["dispersion_cutoff_bohr"] = rc;
    return e;
}

static const double kBoxBohr = 24.0;
static const int kWaters = 21;

// (b) brute force, (c) supercell, (e) cut-off reduction
static void testBox()
{
    const double L = kBoxBohr / A2B;
    const Mol box = waterBox(kWaters, L, 20261103);
    const double rc = 11.5;
    const Run r1 = runGFNFF(box, L, boxFlags(rc), 1);
    const Run r4 = runGFNFF(box, L, boxFlags(rc), 4);
    check(r1.ok && r4.ok, "water box with cell runs");
    if (!r1.ok)
        return;
    int missing = 0, self = 0;
    double e_central = 0.0;
    const double eb = bruteForce(box, kBoxBohr, rc, r1.pairs, missing, self, &e_central);
    const double d = std::abs(r1.e_disp - eb);
    check(d < 1e-12, "box: E_disp == 27-image brute force (" + sci(d) + ")");
    check(missing == 0, "box: every pair with an image within r_c is in the list (" + std::to_string(missing) + " missing)");
    check(self == 0, "box: no self-image within r_c");
    check(std::abs(eb - e_central) > 1e-4, "box: images contribute (test not vacuous)");
    check(std::abs(r1.e_disp - r4.e_disp) < 1e-14, "box: 1 and 4 threads agree");
    check(std::abs(r1.radii.eval - rc) == 0.0 && !r1.radii.capped, "box: explicit r_c is used unchanged");
    std::printf("box %d atoms, L %.1f Bohr, r_c %.1f: E_disp %.15f, brute force %.15f, |d| %.2e, image share %.3e Eh, %zu pairs, t1-t4 %.1e\n",
        box.m_number_atoms, kBoxBohr, rc, r1.e_disp, eb, d, eb - e_central, r1.pairs.size(), std::abs(r1.e_disp - r4.e_disp));

    // (c) 2x1x1 supercell (only x doubled; the cell is non-cubic, L_min/2 stays 12 Bohr)
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
        sc.m_has_pbc = true;
        Eigen::Matrix3d cell = Eigen::Matrix3d::Identity() * L;
        cell(0, 0) = 2.0 * L;
        sc.m_unit_cell = cell;
        Run rs;
        try {
            GFNFF ff(config(1, boxFlags(rc)));
            if (ff.InitialiseMolecule(sc)) {
                ff.Calculation(false);
                rs.e_disp = ff.DispersionEnergy();
                rs.ok = true;
            }
        } catch (const std::exception& e) {
            std::printf("supercell threw: %s\n", e.what());
        }
        check(rs.ok, "supercell runs");
        const double ds = std::abs(rs.e_disp - 2.0 * r1.e_disp);
        check(ds < 1e-12, "supercell 2x1x1: E_disp == 2 x cell (" + sci(ds) + ")");
        std::printf("supercell 2x1x1: E_disp %.15f, 2 x cell %.15f, |d| %.2e\n", rs.e_disp, 2.0 * r1.e_disp, ds);
    }

    // (e) default cut-off reduced to L/2; explicit too large refused
    {
        const Run rd = runGFNFF(box, L, boxFlags(0.0), 1);
        check(rd.ok, "default cut-off in a small cell runs (reduced, not refused)");
        if (rd.ok) {
            check(rd.radii.capped && rd.radii.eval < 0.5 * kBoxBohr && rd.radii.eval > 0.5 * kBoxBohr - 1e-12,
                "default cut-off reduced to just below L/2 (" + std::to_string(rd.radii.eval) + ")");
            check(std::abs(rd.radii.build - (rd.radii.eval + 10.0)) < 1e-12, "reduced list keeps the 10 Bohr skin");
            int mi = 0, si = 0;
            const double ebd = bruteForce(box, kBoxBohr, rd.radii.eval, rd.pairs, mi, si);
            const double dd = std::abs(rd.e_disp - ebd);
            check(dd < 1e-12 && mi == 0 && si == 0, "reduced cut-off: E_disp == brute force at r_c = L/2 (" + sci(dd) + ")");
            std::printf("default cut-off: reduced 50 -> %.12f Bohr, list %.4f, E_disp %.15f, brute force |d| %.2e\n",
                rd.radii.eval, rd.radii.build, rd.e_disp, dd);
        }
        for (double big : { 12.5, 50.0 }) {
            const Run rb = runGFNFF(box, L, boxFlags(big), 1);
            check(rb.threw && !rb.ok, "explicit dispersion_cutoff_bohr " + std::to_string(big) + " >= L/2 is refused");
        }
        // without a cell the radii are the former ones
        const Run rn = runGFNFF(box, 0.0, boxFlags(0.0), 1);
        check(rn.ok && rn.radii.eval == 50.0 && rn.radii.build == 60.0 && !rn.radii.capped, "no cell: radii 50 / 60");
    }
}

// (d) finite difference near the cell faces
static void testFiniteDifference()
{
    const double L = kBoxBohr / A2B;
    const Mol box = waterBox(kWaters, L, 20261103);
    const json flags = boxFlags(0.0);  // reduced default cut-off, positive skin
    const Run on = runGFNFF(box, L, flags, 1, true);
    const Run off = runGFNFF(box, L, flags, 1, false);
    check(on.ok && off.ok, "FD reference runs");
    if (!on.ok || !off.ok)
        return;
    const Matrix g_disp = on.grad - off.grad;  // Eh/Bohr, pair + dC6/dCN chain rule
    const Matrix g_pair = pairGradient(box, kBoxBohr, on.radii.eval, on.pairs);
    const double ws_diff = on.grad_ws.rows() == g_disp.rows() ? (on.grad_ws - g_disp).cwiseAbs().maxCoeff() : 1e300;
    check(ws_diff < 1e-14, "G(disp on) - G(disp off) == workspace dispersion gradient (" + sci(ws_diff) + ")");
    // atoms closest to a face (wrapped into [0, L))
    std::vector<std::pair<double, int>> dist;
    for (int i = 0; i < box.m_number_atoms; ++i) {
        double dmin = 1e300;
        for (int k = 0; k < 3; ++k) {
            double x = box.m_geometry(i, k) - L * std::floor(box.m_geometry(i, k) / L);
            dmin = std::min({ dmin, x, L - x });
        }
        dist.push_back({ dmin, i });
    }
    std::sort(dist.begin(), dist.end());
    const int natoms_fd = 6;
    auto eAt = [&](const Mol& m) {
        const Run r = runGFNFF(m, L, flags, 1, true);
        return r.e_disp;
    };
    double max_g = 0.0, max_cn = 0.0;
    auto fdError = [&](double h) {
        double err = 0.0;
        for (int s = 0; s < natoms_fd; ++s) {
            const int a = dist[s].second;
            for (int k = 0; k < 3; ++k) {
                Mol p = box, mm = box;
                p.m_geometry(a, k) += h;
                mm.m_geometry(a, k) -= h;
                const double fd = (eAt(p) - eAt(mm)) / (2.0 * h);  // Eh/A
                const double an = g_disp(a, k) * A2B;              // Eh/Bohr -> Eh/A
                err = std::max(err, std::abs(fd - an));
                max_g = std::max(max_g, std::abs(an));
                max_cn = std::max(max_cn, std::abs((g_disp(a, k) - g_pair(a, k)) * A2B));
            }
        }
        return err;
    };
    const double e2 = fdError(2e-4), e1 = fdError(1e-4);
    check(e1 < 1e-8, "dispersion gradient vs central FD near the faces, max error " + sci(e1) + " Eh/A");
    check(e1 < e2 / 3.0, "FD error falls quadratically with h (" + sci(e2) + " -> " + sci(e1) + ")");
    check(max_cn > 1e-7, "dE/dCN part of the gradient is non-zero (test not vacuous): " + sci(max_cn));
    std::printf("FD (6 atoms nearest a face, %.2f..%.2f A): max |g_disp| %.3e Eh/A, dE/dCN part max %.3e, |an - FD| %.2e (h 2e-4 A), %.2e (h 1e-4 A)\n",
        dist[0].first, dist[natoms_fd - 1].first, max_g, max_cn, e2, e1);
}

// (f) Verlet list under wrapping and small moves
static void testVerlet()
{
    const double L = kBoxBohr / A2B;
    Mol box = waterBox(kWaters, L, 20261103);
    box.m_has_pbc = true;
    box.m_unit_cell = Eigen::Matrix3d::Identity() * L;
    try {
        GFNFF ff(config(1, boxFlags(0.0)));
        check(ff.InitialiseMolecule(box), "Verlet: init");
        ff.Calculation(true);
        const double e0 = ff.DispersionEnergy();
        const long c0 = ff.dispersionRebuildCount();
        Matrix g = box.m_geometry;
        for (int i = 0; i < 3; ++i)
            g(i, 0) += L;  // first water by one lattice vector
        ff.UpdateMolecule(g);
        ff.Calculation(true);
        const double e1 = ff.DispersionEnergy();
        const long c1 = ff.dispersionRebuildCount();
        check(std::abs(e1 - e0) < 1e-12, "wrap by +L: E_disp unchanged (" + sci(std::abs(e1 - e0)) + ")");
        check(c1 == c0, "wrap by +L: no D4 list rebuild (" + std::to_string(c0) + " -> " + std::to_string(c1) + ")");
        // small random moves below the half-skin (5 Bohr): no rebuild, energy equals a fresh instance
        std::mt19937 rng(7);
        std::uniform_real_distribution<double> u(-0.02, 0.02);  // Angstrom per step (topology stays the same)
        for (int step = 0; step < 4; ++step) {
            for (int i = 0; i < g.rows(); ++i)
                for (int k = 0; k < 3; ++k)
                    g(i, k) += u(rng);
            ff.UpdateMolecule(g);
            ff.Calculation(true);
        }
        Mol moved = box;
        moved.m_geometry = g;
        const Run fresh = runGFNFF(moved, L, boxFlags(0.0), 1);
        const double df = std::abs(ff.DispersionEnergy() - fresh.e_disp);
        check(fresh.ok && df < 1e-12, "moves below the half-skin: E_disp == fresh instance (" + sci(df) + ")");
        check(ff.dispersionRebuildCount() == c0, "moves below the half-skin: no rebuild");
        // one water by 3 Angstrom (5.7 Bohr > half-skin 5 Bohr): rebuild
        for (int i = 3; i < 6; ++i)
            g(i, 1) += 3.0;
        ff.UpdateMolecule(g);
        ff.Calculation(true);
        check(ff.dispersionRebuildCount() == c0 + 1, "move beyond the half-skin rebuilds once");
        std::printf("Verlet: wrap |dE| %.2e, rebuilds %ld -> %ld; after 4 small moves |E - fresh| %.2e; after 5.7 Bohr move rebuilds %ld\n",
            std::abs(e1 - e0), c0, c1, df, ff.dispersionRebuildCount());
    } catch (const std::exception& e) {
        check(false, std::string("Verlet test threw: ") + e.what());
    }
}

int main()
{
    CurcumaLogger::set_verbosity(0);
    testDimer();
    testBox();
    testFiniteDifference();
    testVerlet();
    if (g_failures == 0)
        std::printf("pbc_dispersion: all checks passed\n");
    return g_failures == 0 ? 0 : 1;
}
