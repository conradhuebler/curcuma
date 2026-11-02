/*
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 */
/**
 * Periodic GFN-FF coordination number (WP-PERIODIC-NONBONDED A2) - Claude Generated (Nov 2026).
 *
 * Checks ONLY the coordination number (value and derivative) under the minimum-image convention.
 * Dispersion, repulsion, Coulomb/EEQ, HB/XB and the bonded terms have no minimum image yet, so the
 * total energy and gradient of a periodic system are NOT tested here (they are not periodic yet).
 *
 * (1) Periodic SpatialCellList: pairs within the cutoff equal a brute-force minimum-image scan for
 *     random points (grids of 1, 2 and >= 3 cells per axis).
 * (2) Shift invariance (GFN-FF path, prepareCNAndEEQ): CN of the water dimer without a cell equals
 *     the CN with a 130 Bohr cell (a) after shifting one molecule by L, (b) after tearing one molecule
 *     across the y boundary and (c) after both; difference < 1e-12. The torn variant differs from the
 *     non-periodic CN of the same coordinates by O(1), so the check is not vacuous.
 * (3) Supercell: water dimer in a small (8, 8, 9) Bohr cell with a 3.9 Bohr CN cutoff (images do
 *     interact); 2x1x1 and 1x1x2 replicas give the same CN per atom (< 1e-12), and both equal a
 *     brute-force sum over the 27 neighbouring images. A 4x4x4 replica (384 atoms) in a 128 Bohr cell
 *     through GFN-FF: CN equal to the single cell per atom and bitwise equal for 1 and 4 threads.
 * (4) Derivative: the CN-derivative store of GFN-FF (getLastCNDerivatives, Bohr) for the torn dimer
 *     against a central finite difference of the periodic CN (Angstrom displacement, converted to
 *     per Angstrom), and against the derivative of the intact dimer without a cell.
 *
 * Geometry from the structure library (water-dimer.roo2.8 via the registry), only translated.
 * Exit 0 when every check passes; each failure prints one line.
 */
#include "src/core/curcuma_logger.h"
#include "src/core/energy_calculators/ff_methods/cn_calculator.h"
#include "src/core/energy_calculators/ff_methods/gfnff.h"
#include "src/core/energy_calculators/ff_methods/gfnff_par.h"
#include "src/core/molecule.h"
#include "src/core/units.h"
#include "src/tools/pbc_utils.h"
#include "src/tools/spatial_cell_list.h"
#include "test_molecule_registry.h"
#include "json.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <set>
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

static const double A2B = CurcumaUnit::Length::ANGSTROM_TO_BOHR;

static PBCUtils::OrthoCell orthoCell(double lx, double ly, double lz)
{
    Eigen::Matrix3d m = Eigen::Matrix3d::Zero();
    m(0, 0) = lx;
    m(1, 1) = ly;
    m(2, 2) = lz;
    std::string err;
    PBCUtils::OrthoCell c = PBCUtils::makeOrthoCell(m, err);
    check(c.active, "makeOrthoCell");
    return c;
}

static Mol dimerAngstrom()
{
    return TestMolecules::TestMoleculeRegistry::createMolecule("H2O_dimer", false).getMolInfo();
}

static double maxAbsDiff(const Vector& a, const Vector& b)
{
    if (a.size() != b.size())
        return 1e300;
    return (a - b).cwiseAbs().maxCoeff();
}

// ---------------------------------------------------------------------------------------------
// (1) periodic cell list against brute force
static void testCellList()
{
    std::mt19937 rng(20261102);
    struct Case {
        double L[3];
        double cut;
    };
    // cells per axis: (1,1,1) / (2,2,2) / (4,3,5) / (2,6,3)
    const Case cases[] = { { { 9.0, 9.0, 9.0 }, 4.4 }, { { 10.0, 10.0, 10.0 }, 4.9 },
        { { 20.0, 15.0, 25.0 }, 4.9 }, { { 12.0, 36.0, 18.0 }, 5.9 } };
    for (const Case& cs : cases) {
        const PBCUtils::OrthoCell cell = orthoCell(cs.L[0], cs.L[1], cs.L[2]);
        const int n = 300;
        Eigen::MatrixXd x(n, 3);
        for (int k = 0; k < 3; ++k) {
            // points also outside [0, L): the list must wrap them itself
            std::uniform_real_distribution<double> u(-1.3 * cs.L[k], 2.2 * cs.L[k]);
            for (int i = 0; i < n; ++i)
                x(i, k) = u(rng);
        }
        std::set<std::pair<int, int>> ref, got, nb;
        const double c2 = cs.cut * cs.cut;
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j) {
                double dx = x(i, 0) - x(j, 0), dy = x(i, 1) - x(j, 1), dz = x(i, 2) - x(j, 2);
                PBCUtils::applyMinimumImageOrtho(cell, dx, dy, dz);
                if (dx * dx + dy * dy + dz * dz <= c2)
                    ref.insert({ i, j });
            }
        SpatialCellList list;
        list.build(x, cs.cut, cell);
        check(list.periodic(), "cell list built periodic");
        int dup = 0;
        list.forEachPair([&](int i, int j, double) {
            if (!got.insert({ i, j }).second)
                ++dup;
        });
        for (int i = 0; i < n; ++i)
            list.forEachNeighbor(i, c2, [&](int j, double) {
                if (j > i && !nb.insert({ i, j }).second)
                    ++dup;
            });
        check(dup == 0, "periodic cell list: no pair visited twice");
        check(got == ref, "periodic forEachPair == brute-force minimum image (L " + std::to_string(cs.L[0]) + ")");
        check(nb == ref, "periodic forEachNeighbor == brute-force minimum image (L " + std::to_string(cs.L[0]) + ")");
        std::printf("cell list L=(%.0f,%.0f,%.0f) cut %.1f: %zu pairs, brute force %zu\n", cs.L[0], cs.L[1],
            cs.L[2], cs.cut, got.size(), ref.size());
    }
    // inactive cell: identical to the non-periodic build
    {
        Eigen::MatrixXd x(50, 3);
        std::uniform_real_distribution<double> u(-10.0, 10.0);
        for (int i = 0; i < 50; ++i)
            for (int k = 0; k < 3; ++k)
                x(i, k) = u(rng);
        SpatialCellList a, b;
        a.build(x, 3.0);
        b.build(x, 3.0, PBCUtils::OrthoCell());
        std::vector<std::pair<int, int>> pa, pb;
        a.forEachPair([&](int i, int j, double) { pa.push_back({ i, j }); });
        b.forEachPair([&](int i, int j, double) { pb.push_back({ i, j }); });
        check(!b.periodic() && pa == pb, "inactive cell: same pairs in the same order as build(coords, cutoff)");
    }
}

// ---------------------------------------------------------------------------------------------
// GFN-FF set-up helpers
static json gfnffConfig(int threads)
{
    json c = json::object();
    c["verbosity"] = 0;
    c["threads"] = threads;
    c["cache_topology"] = false;
    c["coulomb_ewald"] = false;  // Nov 2026 (B3): this test isolates terms other than Coulomb/EEQ; keep the plain (non-periodic) Coulomb and EEQ of A1-A4
    return c;
}

struct CNRun {
    bool ok = false;
    Vector cn;
    Matrix dcn_per_bohr;  // (3N x N): column k = dCN_k / dx (Bohr), row 3*i+d
};

static Matrix denseDCN(const CNDerivStore& s, int n)
{
    Matrix out = Matrix::Zero(3 * n, n);
    for (int k = 0; k < n; ++k) {
        Vector e = Vector::Zero(n);
        e(k) = 1.0;
        GeoGradMatrix g = GeoGradMatrix::Zero(n, 3);
        s.applyAdd(e, g);
        for (int i = 0; i < n; ++i)
            for (int d = 0; d < 3; ++d)
                out(3 * i + d, k) = g(i, d);
    }
    return out;
}

static CNRun gfnffCN(Mol mol, const Eigen::Matrix3d* cell_ang, int threads = 1)
{
    CNRun r;
    mol.m_has_pbc = (cell_ang != nullptr);
    if (cell_ang)
        mol.m_unit_cell = *cell_ang;
    try {
        GFNFF ff(gfnffConfig(threads));
        if (!ff.InitialiseMolecule(mol))
            return r;
        ff.Calculation(true);
        r.cn = ff.getLastCN();
        r.dcn_per_bohr = denseDCN(ff.getLastCNDerivatives(), mol.m_number_atoms);
        r.ok = (r.cn.size() == mol.m_number_atoms);
    } catch (const std::exception& e) {
        std::printf("GFN-FF threw: %s\n", e.what());
    }
    return r;
}

// ---------------------------------------------------------------------------------------------
// (2) shift invariance through the GFN-FF path; (4) derivative
static void testShiftAndDerivative()
{
    const Mol base = dimerAngstrom();
    const int n = base.m_number_atoms;
    const double L_ang = 130.0 / A2B;  // 130 Bohr: above 2x every default cutoff (D4 build 60)
    const Eigen::Matrix3d cell_ang = Eigen::Matrix3d::Identity() * L_ang;
    const PBCUtils::OrthoCell cell_bohr = orthoCell(130.0, 130.0, 130.0);

    // reference: intact dimer, no cell, placed in the middle of the cell
    Mol centred = base;
    for (int i = 0; i < n; ++i)
        for (int k = 0; k < 3; ++k)
            centred.m_geometry(i, k) += 0.5 * L_ang;
    const CNRun ref = gfnffCN(centred, nullptr);
    check(ref.ok, "GFN-FF no cell runs");

    // atoms 0-2 first water, 3-5 second water (water-dimer.roo2.8)
    auto shiftedSecond = [&](Mol m) {  // (a) second molecule moved by +L along z
        for (int i = 3; i < 6; ++i)
            m.m_geometry(i, 2) += L_ang;
        return m;
    };
    auto tornFirst = [&](Mol m) {  // (b) cell boundary y = 0 through the first water (O-H across it)
        const double y0 = m.m_geometry(0, 1) + 0.3;  // O at y0 - 0.3 -> -0.3 A, H(+0.757) -> +0.457 A
        for (int i = 0; i < n; ++i)
            m.m_geometry(i, 1) -= y0;
        for (int i = 0; i < 3; ++i)
            if (m.m_geometry(i, 1) < 0.0)
                m.m_geometry(i, 1) += L_ang;  // per-atom wrap: the molecule is torn
        return m;
    };

    struct Variant {
        std::string name;
        Mol mol;
    };
    const Variant variants[] = { { "second molecule shifted by L", shiftedSecond(centred) },
        { "first molecule torn across y=0", tornFirst(centred) },
        { "torn and shifted", shiftedSecond(tornFirst(centred)) } };
    for (const Variant& v : variants) {
        const CNRun per = gfnffCN(v.mol, &cell_ang);
        check(per.ok, "GFN-FF with cell runs: " + v.name);
        if (!per.ok || !ref.ok)
            continue;
        const double d = maxAbsDiff(per.cn, ref.cn);
        check(d < 1e-12, v.name + ": periodic CN == non-periodic CN of the intact dimer (" + std::to_string(d) + ")");
        const double dd = (per.dcn_per_bohr - ref.dcn_per_bohr).cwiseAbs().maxCoeff();
        check(dd < 1e-12, v.name + ": periodic dCN == dCN of the intact dimer (" + std::to_string(dd) + ")");
        // the same coordinates without a cell
        const CNRun flat = gfnffCN(v.mol, nullptr);
        const double dflat = flat.ok ? maxAbsDiff(flat.cn, ref.cn) : -1.0;
        std::printf("%-32s |CN_pbc - CN_ref| = %.2e, |dCN| %.2e, same coords without cell: |CN - CN_ref| = %.2e\n",
            v.name.c_str(), d, dd, dflat);
    }
    // the torn variant must be sensitive to the minimum image at all
    {
        const CNRun flat = gfnffCN(variants[1].mol, nullptr);
        check(flat.ok && maxAbsDiff(flat.cn, ref.cn) > 0.1, "torn dimer without a cell has a different CN (test not vacuous)");
    }

    // (4) central finite difference of the periodic CN (CNCalculator, the per-step function) for
    // the torn dimer against the analytic derivative of GFN-FF. h in Angstrom, derivative per A.
    const Mol& torn = variants[2].mol;
    const CNRun per = gfnffCN(torn, &cell_ang);
    if (!per.ok)
        return;
    const std::vector<int> atoms(torn.m_atoms.begin(), torn.m_atoms.end());
    auto cnAt = [&](const Eigen::MatrixXd& geom_ang) {
        const Eigen::MatrixXd gb = geom_ang * A2B;
        const auto res = CNCalculator::calculateGFNFFCNWithNeighbors(atoms, gb, 10.0, -7.5, 4.4, nullptr, 1, &cell_bohr);
        return Vector(Vector::Map(res.cn_values.data(), n));
    };
    const Eigen::MatrixXd g0 = torn.m_geometry;
    check(maxAbsDiff(cnAt(g0), per.cn) == 0.0, "CNCalculator periodic CN == GFN-FF per-step CN bitwise");
    double max_val = 0.0, cross_val = 0.0;
    auto fdError = [&](double h) {
        double max_err = 0.0;
        for (int a = 0; a < n; ++a)
            for (int d = 0; d < 3; ++d) {
                Eigen::MatrixXd gp = g0, gm = g0;
                gp(a, d) += h;
                gm(a, d) -= h;
                const Vector fd = (cnAt(gp) - cnAt(gm)) / (2.0 * h);
                for (int k = 0; k < n; ++k) {
                    const double an = per.dcn_per_bohr(3 * a + d, k) * A2B;  // per Bohr -> per Angstrom
                    max_err = std::max(max_err, std::abs(an - fd(k)));
                    max_val = std::max(max_val, std::abs(an));
                    // O(0)-H(1) of the torn water: O was wrapped by +L, H(1) was not; the pair crosses y = 0
                    if ((a == 1 && k == 0) || (a == 0 && k == 1))
                        cross_val = std::max(cross_val, std::abs(an));
                }
            }
        return max_err;
    };
    const double err2 = fdError(2e-4), err1 = fdError(1e-4);
    check(err1 < 1e-7, "dCN analytic vs central FD (torn dimer, cell) max error " + std::to_string(err1));
    // central difference: error ~ h^2, so halving h must cut it by about 4 (not a constant offset)
    check(err1 < err2 / 3.0, "FD error falls quadratically with h (" + std::to_string(err2) + " -> " + std::to_string(err1) + ")");
    check(cross_val > 1e-2, "derivative of the pair across the boundary is large (test not vacuous)");
    std::printf("dCN/dx torn dimer: max |analytic| %.3e /A, cross-boundary O-H %.3e /A, max |analytic - FD| %.2e (h = 2e-4 A), %.2e (h = 1e-4 A)\n",
        max_val, cross_val, err2, err1);
}

// ---------------------------------------------------------------------------------------------
// (3) supercells
static Eigen::MatrixXd replicate(const Eigen::MatrixXd& g, const double L[3], int nx, int ny, int nz)
{
    const int n = static_cast<int>(g.rows());
    Eigen::MatrixXd out(n * nx * ny * nz, 3);
    int r = 0;
    for (int a = 0; a < nx; ++a)
        for (int b = 0; b < ny; ++b)
            for (int c = 0; c < nz; ++c)
                for (int i = 0; i < n; ++i, ++r) {
                    out(r, 0) = g(i, 0) + a * L[0];
                    out(r, 1) = g(i, 1) + b * L[1];
                    out(r, 2) = g(i, 2) + c * L[2];
                }
    return out;
}

static void testSupercell()
{
    const Mol base = dimerAngstrom();
    const int n = base.m_number_atoms;
    const std::vector<int> atoms(base.m_atoms.begin(), base.m_atoms.end());
    const Eigen::MatrixXd gb = base.m_geometry * A2B;

    // small cell, images interact: z extent of the dimer is ~7.5 Bohr in a 9 Bohr cell
    const double L[3] = { 8.0, 8.0, 9.0 };
    const double cut = 3.9;  // < L_min/2 = 4
    const PBCUtils::OrthoCell cell = orthoCell(L[0], L[1], L[2]);
    const auto one = CNCalculator::calculateGFNFFCNWithNeighbors(atoms, gb, cut, -7.5, 4.4, nullptr, 1, &cell);

    // brute force over the 27 neighbouring images (cut < L/2: at most one image within the cutoff)
    std::vector<double> brute(n, 0.0);
    for (int i = 0; i < n; ++i) {
        const double ri = GFNFFParameters::gfnff_cn_rcov_bohr(atoms[i]);
        double raw = 0.0;
        for (int j = 0; j < n; ++j)
            for (int a = -1; a <= 1; ++a)
                for (int b = -1; b <= 1; ++b)
                    for (int c = -1; c <= 1; ++c) {
                        if (i == j && a == 0 && b == 0 && c == 0)
                            continue;
                        const double dx = gb(i, 0) - gb(j, 0) - a * L[0];
                        const double dy = gb(i, 1) - gb(j, 1) - b * L[1];
                        const double dz = gb(i, 2) - gb(j, 2) - c * L[2];
                        const double r2 = dx * dx + dy * dy + dz * dz;
                        if (!(r2 < cut * cut))
                            continue;
                        const double rc = ri + GFNFFParameters::gfnff_cn_rcov_bohr(atoms[j]);
                        raw += 0.5 * (1.0 + std::erf(-7.5 * (std::sqrt(r2) - rc) / rc));
                    }
        brute[i] = std::log(1.0 + std::exp(4.4)) - std::log(1.0 + std::exp(4.4 - raw));
    }
    // reference without a cell, to show the images contribute
    const auto flat = CNCalculator::calculateGFNFFCNWithNeighbors(atoms, gb, cut, -7.5, 4.4);
    double d_brute = 0.0, d_flat = 0.0;
    for (int i = 0; i < n; ++i) {
        d_brute = std::max(d_brute, std::abs(one.cn_values[i] - brute[i]));
        d_flat = std::max(d_flat, std::abs(one.cn_values[i] - flat.cn_values[i]));
    }
    check(d_brute < 1e-12, "small cell: periodic CN == brute-force image sum (" + std::to_string(d_brute) + ")");
    check(d_flat > 1e-3, "small cell: images change the CN (test not vacuous)");
    std::printf("small cell (8,8,9) Bohr, cut 3.9: |CN - brute force| %.2e, image contribution max %.3e\n", d_brute, d_flat);

    // the threshold (N^2) and mode functions agree with the neighbour-list pass
    {
        const auto thr = CNCalculator::calculateGFNFFCN(atoms, gb, cut * cut, -7.5, 4.4, nullptr, &cell);
        const auto mode = CNCalculator::calculateGFNFFCN(atoms, gb, cut, 1.0, -7.5, 4.4, &cell);
        double d1 = 0.0, d2 = 0.0;
        for (int i = 0; i < n; ++i) {
            d1 = std::max(d1, std::abs(thr[i] - one.cn_values[i]));
            d2 = std::max(d2, std::abs(mode[i] - one.cn_values[i]));
        }
        check(d1 < 1e-12, "threshold CN with cell == neighbour-list CN (" + std::to_string(d1) + ")");
        check(d2 < 1e-12, "mode CN with cell == neighbour-list CN (" + std::to_string(d2) + ")");
    }

    const int reps[2][3] = { { 2, 1, 1 }, { 1, 1, 2 } };
    for (const auto& rp : reps) {
        const Eigen::MatrixXd gs = replicate(gb, L, rp[0], rp[1], rp[2]);
        const PBCUtils::OrthoCell sc = orthoCell(L[0] * rp[0], L[1] * rp[1], L[2] * rp[2]);
        std::vector<int> as;
        for (int r = 0; r < rp[0] * rp[1] * rp[2]; ++r)
            as.insert(as.end(), atoms.begin(), atoms.end());
        const auto big = CNCalculator::calculateGFNFFCNWithNeighbors(as, gs, cut, -7.5, 4.4, nullptr, 1, &sc);
        double d = 0.0;
        for (int i = 0; i < static_cast<int>(as.size()); ++i)
            d = std::max(d, std::abs(big.cn_values[i] - one.cn_values[i % n]));
        check(d < 1e-12, "supercell " + std::to_string(rp[0]) + "x" + std::to_string(rp[1]) + "x" + std::to_string(rp[2]) + ": CN per atom equal (" + std::to_string(d) + ")");
        std::printf("supercell %dx%dx%d: max |CN_super - CN_cell| %.2e\n", rp[0], rp[1], rp[2], d);
    }

    // 4x4x4 replica through GFN-FF (cell check passes: 128 Bohr), threads 1 and 4
    {
        const double L32[3] = { 32.0, 32.0, 32.0 };
        const PBCUtils::OrthoCell c32 = orthoCell(32.0, 32.0, 32.0);
        const auto cn32 = CNCalculator::calculateGFNFFCNWithNeighbors(atoms, gb, 10.0, -7.5, 4.4, nullptr, 1, &c32);
        Mol big = base;
        big.m_geometry = replicate(gb, L32, 4, 4, 4) / A2B;
        big.m_number_atoms = static_cast<int>(big.m_geometry.rows());
        big.m_atoms.clear();
        for (int r = 0; r < 64; ++r)
            big.m_atoms.insert(big.m_atoms.end(), base.m_atoms.begin(), base.m_atoms.end());
        const Eigen::Matrix3d cell128 = Eigen::Matrix3d::Identity() * (128.0 / A2B);
        const CNRun t1 = gfnffCN(big, &cell128, 1);
        const CNRun t4 = gfnffCN(big, &cell128, 4);
        check(t1.ok && t4.ok, "GFN-FF 4x4x4 replica runs");
        if (t1.ok && t4.ok) {
            check(t1.cn == t4.cn, "4x4x4 replica: CN bitwise equal for 1 and 4 threads");
            double d = 0.0;
            for (int i = 0; i < big.m_number_atoms; ++i)
                d = std::max(d, std::abs(t1.cn(i) - cn32.cn_values[i % n]));
            check(d < 1e-12, "4x4x4 replica: CN per atom == single 32 Bohr cell (" + std::to_string(d) + ")");
            std::printf("4x4x4 replica (384 atoms, 128 Bohr): max |CN - CN_cell| %.2e, threads 1 vs 4 bitwise %s\n",
                d, t1.cn == t4.cn ? "equal" : "DIFFERENT");
        }
    }
}

// Claude Generated (Nov 2026, A2): cn_cutoff_bohr = 0 means the reference threshold
// cnthr = 100 - 50 log10(acc) Bohr^2 (gfnff_param.f90:551) in every call; the mode overload
// and the threshold overload fed by gfnffAccuracyThresholdSqBohr() must agree bitwise.
static void testAccuracyThreshold()
{
    check(CNCalculator::gfnffAccuracyThresholdSqBohr(1.0) == 100.0, "cnthr(1) = 100 Bohr^2");
    check(std::abs(CNCalculator::gfnffAccuracyThresholdSqBohr(0.1) - 150.0) < 1e-12, "cnthr(0.1) = 150 Bohr^2");
    check(CNCalculator::gfnffAccuracyThresholdSqBohr(0.0) >= 1e30, "cnthr(0) = all pairs");
    const Mol m = dimerAngstrom();
    const std::vector<int> atoms(m.m_atoms.begin(), m.m_atoms.end());
    const Eigen::MatrixXd gb = m.m_geometry * A2B;
    for (double acc : { 1.0, 0.1, 0.0 }) {
        const auto mode = CNCalculator::calculateGFNFFCN(atoms, gb, 0.0, acc, -7.5, 4.4);
        const auto thr = CNCalculator::calculateGFNFFCN(atoms, gb,
            CNCalculator::gfnffAccuracyThresholdSqBohr(acc), -7.5, 4.4, nullptr);
        check(mode == thr, "cn_cutoff_bohr 0: mode CN == threshold CN bitwise, acc " + std::to_string(acc));
    }
}

int main()
{
    CurcumaLogger::set_verbosity(0);
    testCellList();
    testShiftAndDerivative();
    testSupercell();
    testAccuracyThreshold();
    if (g_failures == 0)
        std::printf("pbc_cn: all checks passed\n");
    return g_failures == 0 ? 0 : 1;
}
