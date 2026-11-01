/*
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 */
/**
 * Periodic GFN-FF foundation (WP-PERIODIC-NONBONDED A1) - Claude Generated (Nov 2026).
 *
 * (a) PBCUtils::minimumImageOrtho: pair across the boundary, on the boundary, exact magnitude,
 *     symmetry mic(-d) = -mic(d), agreement with the general applyMinimumImage().
 * (b) GFN-FF cell check: a cell too small for the cutoffs (and a triclinic cell) raises an error,
 *     no cell and a large cell do not; without minimum image in the energy terms (A1) the energy
 *     with a large cell equals the energy without a cell bit for bit.
 * (c) PBCUtils::wrapFragmentsIntoCell on a water dimer: every intramolecular distance kept to
 *     1e-12 Angstrom, each centre of mass in fractional [0,1).
 *
 * Exit 0 when every check passes; each failure prints one line.
 */
#include "src/core/energycalculator.h"
#include "src/core/molecule.h"
#include "src/core/curcuma_logger.h"
#include "src/core/elements.h"
#include "src/tools/pbc_utils.h"
#include "test_molecule_registry.h"
#include "json.hpp"

#include <cmath>
#include <cstdio>
#include <random>
#include <string>

using json = nlohmann::json;

static int g_failures = 0;
static void check(bool ok, const std::string& what)
{
    if (!ok) {
        ++g_failures;
        std::printf("FAIL: %s\n", what.c_str());
    }
}

static PBCUtils::OrthoCell cubicCell(double L)
{
    std::string err;
    PBCUtils::OrthoCell c = PBCUtils::makeOrthoCell(Eigen::Matrix3d::Identity() * L, err);
    check(c.active && err.empty(), "makeOrthoCell accepts a cube");
    return c;
}

// ---------------------------------------------------------------------------------------------
static void testMinimumImage()
{
    // L = 16 is a power of two: d/L and L*s are exact, so the expected values are exact too.
    const PBCUtils::OrthoCell c = cubicCell(16.0);
    const double L = c.L[0], iL = c.inv_L[0];
    check(iL == 0.0625, "inv_L of a 16 cube is exactly 1/16");

    // across the boundary: +10 -> -6, -10 -> +6, +23 -> +7, -41 -> +7 (-41 + 3*16)
    check(PBCUtils::minimumImageOrtho(10.0, L, iL) == -6.0, "mic(+10) == -6 in L=16");
    check(PBCUtils::minimumImageOrtho(-10.0, L, iL) == 6.0, "mic(-10) == +6 in L=16");
    check(PBCUtils::minimumImageOrtho(23.0, L, iL) == 7.0, "mic(+23) == +7 in L=16");
    check(PBCUtils::minimumImageOrtho(-41.0, L, iL) == 7.0, "mic(-41) == +7 in L=16");
    // inside: unchanged
    check(PBCUtils::minimumImageOrtho(3.25, L, iL) == 3.25, "mic(3.25) unchanged");
    check(PBCUtils::minimumImageOrtho(0.0, L, iL) == 0.0, "mic(0) == 0");
    // on the boundary: |s| = 0.5 keeps the sign (rint(+-0.5) = 0), magnitude exactly L/2
    check(PBCUtils::minimumImageOrtho(8.0, L, iL) == 8.0, "mic(+L/2) == +L/2");
    check(PBCUtils::minimumImageOrtho(-8.0, L, iL) == -8.0, "mic(-L/2) == -L/2");
    check(std::abs(PBCUtils::minimumImageOrtho(24.0, L, iL)) == 8.0, "|mic(1.5 L)| == L/2");
    check(std::abs(PBCUtils::minimumImageOrtho(-40.0, L, iL)) == 8.0, "|mic(-2.5 L)| == L/2");

    // inactive cell: no-op
    {
        PBCUtils::OrthoCell none;
        double dx = 123.4, dy = -77.0, dz = 0.5;
        PBCUtils::applyMinimumImageOrtho(none, dx, dy, dz);
        check(dx == 123.4 && dy == -77.0 && dz == 0.5, "inactive cell leaves the vector unchanged");
    }

    // random displacements in a non-power-of-two box (Bohr scale of a water box)
    std::string err;
    Eigen::Matrix3d m = Eigen::Matrix3d::Zero();
    m(0, 0) = 37.9 * 1.8897261246;
    m(1, 1) = 23.1;
    m(2, 2) = 51.7;
    const PBCUtils::OrthoCell box = PBCUtils::makeOrthoCell(m, err);
    check(box.active, "makeOrthoCell accepts a rectangular box");
    const Eigen::Matrix3d m_inv = m.inverse();
    std::mt19937 rng(20261101);
    std::uniform_real_distribution<double> u(-200.0, 200.0);
    int bad_sym = 0, bad_range = 0, bad_lattice = 0, bad_general = 0;
    for (int n = 0; n < 20000; ++n) {
        double d[3] = { u(rng), u(rng), u(rng) };
        double p[3] = { d[0], d[1], d[2] };
        double q[3] = { -d[0], -d[1], -d[2] };
        PBCUtils::applyMinimumImageOrtho(box, p[0], p[1], p[2]);
        PBCUtils::applyMinimumImageOrtho(box, q[0], q[1], q[2]);
        const Eigen::Vector3d g = PBCUtils::applyMinimumImage(Eigen::Vector3d(d[0], d[1], d[2]), m, m_inv);
        for (int k = 0; k < 3; ++k) {
            if (q[k] != -p[k]) ++bad_sym;
            if (std::abs(p[k]) > 0.5 * box.L[k] * (1.0 + 1e-15)) ++bad_range;
            const double nimg = (d[k] - p[k]) / box.L[k];
            if (std::abs(nimg - std::round(nimg)) > 1e-12) ++bad_lattice;
            if (std::abs(g(k) - p[k]) > 1e-12 * box.L[k]) ++bad_general;
        }
    }
    check(bad_sym == 0, "mic(-d) == -mic(d) bitwise (" + std::to_string(bad_sym) + " violations)");
    check(bad_range == 0, "|mic(d)| <= L/2 (" + std::to_string(bad_range) + " violations)");
    check(bad_lattice == 0, "d - mic(d) is a lattice vector (" + std::to_string(bad_lattice) + " violations)");
    check(bad_general == 0, "agrees with applyMinimumImage (" + std::to_string(bad_general) + " violations)");

    // triclinic refused, 90/90/90 from buildLatticeVectors accepted
    check(!PBCUtils::isOrthorhombic(PBCUtils::buildLatticeVectors(10, 10, 10, 90, 90, 60)), "triclinic cell refused");
    check(PBCUtils::isOrthorhombic(PBCUtils::buildLatticeVectors(10, 11, 12, 90, 90, 90)), "90/90/90 cell accepted");
    PBCUtils::OrthoCell tri = PBCUtils::makeOrthoCell(PBCUtils::buildLatticeVectors(10, 10, 10, 90, 90, 60), err);
    check(!tri.active && !err.empty(), "makeOrthoCell reports a triclinic cell");

    // the check function itself
    const PBCUtils::OrthoCell c20 = cubicCell(20.0);
    check(PBCUtils::checkCutoffsAgainstCell(c20, { { "a", 9.99 } }, "Bohr").empty(), "9.99 < 10 passes");
    check(!PBCUtils::checkCutoffsAgainstCell(c20, { { "a", 10.0 } }, "Bohr").empty(), "10 == L/2 fails");
    check(PBCUtils::checkCutoffsAgainstCell(PBCUtils::OrthoCell(), { { "a", 1e9 } }, "Bohr").empty(),
        "no cell: nothing is checked");
}

// ---------------------------------------------------------------------------------------------
static curcuma::Molecule waterDimer()
{
    // Structure library id water-dimer.roo2.8 (registry rule: no hardcoded geometry), Angstrom.
    return TestMolecules::TestMoleculeRegistry::createMolecule("H2O_dimer", false);
}

// Returns (error flag of the calculator, energy); the cell is in Angstrom.
static std::pair<bool, double> gfnffEnergy(const Eigen::Matrix3d* cell)
{
    json config = { { "verbosity", 0 }, { "threads", 1 }, { "gfnff", json::object() }, { "gpu", "none" } };
    curcuma::Molecule mol = waterDimer();
    if (cell)
        mol.setUnitCell(*cell, true);
    EnergyCalculator calc("gfnff", config);
    calc.setMolecule(mol.getMolInfo());
    if (calc.Error())
        return { true, 0.0 };
    const double e = calc.CalculateEnergy(true);
    return { calc.Error(), e };
}

static void testCutoffCheck()
{
    const auto none = gfnffEnergy(nullptr);
    check(!none.first, "no cell: GFN-FF runs without error");

    // 20 A cube: L/2 = 18.9 Bohr, below the CN-derivative range (40) and the D4 radii (50/60)
    const Eigen::Matrix3d small = Eigen::Matrix3d::Identity() * 20.0;
    check(gfnffEnergy(&small).first, "20 A cell: cutoff check raises an error");

    // triclinic cell: refused
    const Eigen::Matrix3d tri = PBCUtils::buildLatticeVectors(100, 100, 100, 90, 90, 60);
    check(gfnffEnergy(&tri).first, "triclinic cell: refused with an error");

    // 100 A cube: L/2 = 94.5 Bohr > every default cutoff (largest: D4 build 60)
    const Eigen::Matrix3d large = Eigen::Matrix3d::Identity() * 100.0;
    const auto big = gfnffEnergy(&large);
    check(!big.first, "100 A cell: no error");
    // A1 applies no minimum image yet: the energy must be the non-periodic one, bit for bit.
    check(big.second == none.second, "100 A cell: energy equals the no-cell energy bit for bit");
    std::printf("water dimer GFN-FF: no cell %.15f Eh, 100 A cell %.15f Eh\n", none.second, big.second);
}

// ---------------------------------------------------------------------------------------------
static void testFragmentWrap()
{
    curcuma::Molecule mol = waterDimer();
    const double L = 10.0;  // Angstrom
    const Eigen::Matrix3d cell = Eigen::Matrix3d::Identity() * L;
    mol.setUnitCell(cell, true);

    // Move the dimer so that both waters sit outside the cell, the second one straddling
    // the z face after the wrap would apply per atom: shift (+13.3, -4.2, +7.1).
    Geometry geom = mol.getGeometry();
    for (int i = 0; i < geom.rows(); ++i) {
        geom(i, 0) += 13.3;
        geom(i, 1) -= 4.2;
        geom(i, 2) += 7.1;
    }
    mol.setGeometry(geom);
    const std::vector<std::vector<int>> fragments = mol.GetFragments();
    check(fragments.size() == 2, "water dimer: two fragments (got " + std::to_string(fragments.size()) + ")");

    std::vector<double> masses(geom.rows());
    for (int i = 0; i < geom.rows(); ++i)
        masses[i] = Elements::AtomicMass[mol.Atom(i).first];

    const Geometry before = geom;
    const int moved = PBCUtils::wrapFragmentsIntoCell(geom, masses, fragments, cell, cell.inverse());
    check(moved == 2, "both fragments moved (got " + std::to_string(moved) + ")");

    double max_dev = 0.0;
    for (const auto& frag : fragments) {
        for (std::size_t a = 0; a < frag.size(); ++a)
            for (std::size_t b = a + 1; b < frag.size(); ++b) {
                const double r0 = (before.row(frag[a]) - before.row(frag[b])).norm();
                const double r1 = (geom.row(frag[a]) - geom.row(frag[b])).norm();
                max_dev = std::max(max_dev, std::abs(r1 - r0));
            }
        double msum = 0.0;
        Eigen::Vector3d com = Eigen::Vector3d::Zero();
        for (int a : frag) {
            msum += masses[a];
            com += masses[a] * Eigen::Vector3d(geom(a, 0), geom(a, 1), geom(a, 2));
        }
        com /= msum;
        for (int k = 0; k < 3; ++k)
            check(com(k) >= 0.0 && com(k) < L, "centre of mass inside the cell, axis " + std::to_string(k));
    }
    check(max_dev <= 1e-12, "intramolecular distances kept to 1e-12 A (max dev " + std::to_string(max_dev) + ")");
    std::printf("fragment wrap: %d fragments moved, max intramolecular deviation %.3e A\n", moved, max_dev);

    // A second wrap is a no-op.
    const Geometry once = geom;
    check(PBCUtils::wrapFragmentsIntoCell(geom, masses, fragments, cell, cell.inverse()) == 0 && geom == once,
        "second wrap moves nothing");
}

int main()
{
    CurcumaLogger::set_verbosity(0);
    testMinimumImage();
    testCutoffCheck();
    testFragmentWrap();
    if (g_failures == 0)
        std::printf("pbc_foundation: all checks passed\n");
    return g_failures == 0 ? 0 : 1;
}
