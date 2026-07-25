/*
 * test_torsion_space.cpp - Unit tests for the torsion-space reduction
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Tests (no energy evaluation, pure topology + geometry):
 *   1. detectTorsions on butane           -> exactly the central C-C bond
 *   2. detectTorsions on a ring molecule  -> ring bonds are rejected
 *   3. setDihedral                        -> hits the target angle, keeps bond lengths/angles,
 *                                            moves only the atoms in `moving`
 *   4. round trip                         -> read angle, set another, set it back == original
 *   5. clusterStates                      -> periodic clustering incl. the +-180 wrap-around
 *   6. stateVector / hammingDistance
 *
 * Claude Generated (Jul 2026)
 */

#include "src/capabilities/torsion_space.h"
#include "core/test_molecule_registry.h"
#include "src/core/curcuma_logger.h"
#include "src/core/molecule.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

static int tests_run = 0;
static int tests_passed = 0;
static int tests_failed = 0;

#define TEST_ASSERT(condition, msg)                        \
    do {                                                   \
        tests_run++;                                       \
        if (condition) {                                   \
            std::cout << "  [PASS] " << msg << std::endl;  \
            tests_passed++;                                \
        } else {                                           \
            std::cout << "  [FAIL] " << msg << std::endl;  \
            tests_failed++;                                \
        }                                                  \
    } while (0)

using namespace curcuma;

// Geometries come from the shared molecule registry (test_cases/CLAUDE.md: NEVER hardcode geometry
// in a test file). scale_coordinates=false keeps Angstrom, which is what Molecule expects.
static Molecule fromRegistry(const std::string& name)
{
    return TestMolecules::TestMoleculeRegistry::createMolecule(name, false);
}

// Heavy-neighbour count from the bond matrix -- lets the assertions below identify "the central
// bond" without hardcoding atom indices (the registry order is not the order of the xyz file).
static int heavyNeighbours(const Molecule& mol, int atom)
{
    const Matrix bonds = mol.DistanceMatrix().second;
    int count = 0;
    for (int b = 0; b < mol.AtomCount(); ++b)
        if (b != atom && bonds(atom, b) > 0.5 && mol.Atom(b).first != 1)
            count++;
    return count;
}

static void test_detect_butane()
{
    std::cout << "\n=== detectTorsions (butane) ===" << std::endl;
    Molecule mol = fromRegistry("butane");
    auto torsions = TorsionSpace::detectTorsions(mol);

    TEST_ASSERT(torsions.size() == 1, "butane has exactly 1 rotatable torsion (found "
            + std::to_string(torsions.size()) + ")");
    if (torsions.empty())
        return;

    const auto& t = torsions.front();
    const bool central = mol.Atom(t.j).first == 6 && mol.Atom(t.k).first == 6
        && heavyNeighbours(mol, t.j) == 2 && heavyNeighbours(mol, t.k) == 2;
    TEST_ASSERT(central, "the rotatable bond connects the two inner carbons (" + t.label(mol) + ")");
    TEST_ASSERT(t.i != t.k && t.l != t.j && t.i >= 0 && t.l >= 0, "reference atoms i/l are valid");
    // Rotating about C1-C2 moves everything on the C2 side except C2 itself (it sits on the axis):
    // the terminal CH3 (C + 3 H) plus the two hydrogens on C2 -> 6 atoms, 1 of them heavy.
    TEST_ASSERT(t.moving.size() == 6, "the whole ethyl side moves without the axis atom (6 atoms, found "
            + std::to_string(t.moving.size()) + ")");
    TEST_ASSERT(t.heavy_moving == 1, "one heavy atom moves");
    TEST_ASSERT(std::find(t.moving.begin(), t.moving.end(), t.k) == t.moving.end(),
        "the axis atom k is not in the moving set");

    // Note: the stored geometry is NOT anti (its backbone dihedral is ~116 deg -- the "Anti
    // conformation" claim in test_cases/validation/butane.xyz does not match its coordinates), so
    // the reference angle is constructed here instead of assumed.
    const Geometry anti = TorsionSpace::setDihedral(mol.getGeometry(), t, 180.0);
    TEST_ASSERT(std::abs(std::abs(TorsionSpace::dihedral(anti, t)) - 180.0) < 1e-6,
        "an explicitly built anti geometry reads back as 180 deg");
}

static void test_detect_ring()
{
    std::cout << "\n=== detectTorsions (rings) ===" << std::endl;

    // Benzene: every bond is either inside the ring or a C-H bond -> nothing rotatable. A ring bond
    // slipping through would show up here immediately.
    Molecule benzene = fromRegistry("C6H6");
    auto benzene_torsions = TorsionSpace::detectTorsions(benzene);
    TEST_ASSERT(benzene_torsions.empty(), "benzene has no rotatable torsion (found "
            + std::to_string(benzene_torsions.size()) + ")");

    // A real ring-containing molecule with side chains: torsions must be found, and every one of
    // them must move a proper, non-empty subset of the molecule (a ring bond would have to move
    // "everything", which the BFS cannot produce).
    Molecule sugar = fromRegistry("monosaccharide");
    auto sugar_torsions = TorsionSpace::detectTorsions(sugar);
    TEST_ASSERT(!sugar_torsions.empty(), "monosaccharide has rotatable torsions (found "
            + std::to_string(sugar_torsions.size()) + ")");

    bool subsets_ok = true, rigid_ok = true, targets_ok = true;
    const Geometry start = sugar.getGeometry();
    for (const auto& t : sugar_torsions) {
        if (t.moving.empty() || static_cast<int>(t.moving.size()) >= sugar.AtomCount() - 1)
            subsets_ok = false;
        Geometry moved = TorsionSpace::setDihedral(start, t, TorsionSpace::dihedral(start, t) + 35.0);
        if (std::abs(TorsionSpace::angleDifference(TorsionSpace::dihedral(moved, t),
                TorsionSpace::dihedral(start, t) + 35.0))
            > 1e-6)
            targets_ok = false;
        for (int a : t.moving)
            for (int b : t.moving) {
                if (a >= b)
                    continue;
                const double d0 = (Eigen::Vector3d(start.row(a)) - Eigen::Vector3d(start.row(b))).norm();
                const double d1 = (Eigen::Vector3d(moved.row(a)) - Eigen::Vector3d(moved.row(b))).norm();
                if (std::abs(d1 - d0) > 1e-10)
                    rigid_ok = false;
            }
    }
    TEST_ASSERT(subsets_ok, "every torsion moves a proper, non-empty part of the molecule");
    TEST_ASSERT(targets_ok, "setDihedral hits its target on all of them (real geometry with rings)");
    TEST_ASSERT(rigid_ok, "the moved fragments stay rigid");
}

static void test_set_dihedral()
{
    std::cout << "\n=== setDihedral ===" << std::endl;
    Molecule mol = fromRegistry("butane");
    auto torsions = TorsionSpace::detectTorsions(mol);
    if (torsions.empty()) {
        TEST_ASSERT(false, "no torsion to test setDihedral on");
        return;
    }
    const auto& t = torsions.front();
    const Geometry start = mol.getGeometry();

    for (double target : { 60.0, -60.0, 90.0, 179.0, -120.0 }) {
        Geometry moved = TorsionSpace::setDihedral(start, t, target);
        const double reached = TorsionSpace::dihedral(moved, t);
        TEST_ASSERT(std::abs(TorsionSpace::angleDifference(reached, target)) < 1e-6,
            "target " + std::to_string(static_cast<int>(target)) + " deg reached (got "
                + std::to_string(reached) + ")");
    }

    // Rigid-body property: bond lengths inside the moved fragment must not change, and atoms
    // outside `moving` must not move at all.
    Geometry moved = TorsionSpace::setDihedral(start, t, 60.0);
    double max_shift_static = 0.0;
    std::vector<char> is_moving(mol.AtomCount(), 0);
    for (int idx : t.moving)
        is_moving[idx] = 1;
    for (int a = 0; a < mol.AtomCount(); ++a)
        if (!is_moving[a])
            max_shift_static = std::max(max_shift_static,
                (Eigen::Vector3d(moved.row(a)) - Eigen::Vector3d(start.row(a))).norm());
    TEST_ASSERT(max_shift_static < 1e-12, "atoms outside `moving` stay exactly in place");

    double max_len_change = 0.0;
    for (int a = 0; a < mol.AtomCount(); ++a)
        for (int b = a + 1; b < mol.AtomCount(); ++b) {
            if (is_moving[a] != is_moving[b])
                continue; // distances across the axis are supposed to change
            const double d0 = (Eigen::Vector3d(start.row(a)) - Eigen::Vector3d(start.row(b))).norm();
            const double d1 = (Eigen::Vector3d(moved.row(a)) - Eigen::Vector3d(moved.row(b))).norm();
            max_len_change = std::max(max_len_change, std::abs(d1 - d0));
        }
    TEST_ASSERT(max_len_change < 1e-10,
        "the moved fragment stays rigid (max internal distance change "
            + std::to_string(max_len_change) + " A)");
}

static void test_round_trip()
{
    std::cout << "\n=== round trip (set away and back) ===" << std::endl;
    Molecule mol = fromRegistry("butane");
    auto torsions = TorsionSpace::detectTorsions(mol);
    if (torsions.empty()) {
        TEST_ASSERT(false, "no torsion for the round trip");
        return;
    }
    const auto& t = torsions.front();
    const Geometry start = mol.getGeometry();
    const double phi0 = TorsionSpace::dihedral(start, t);

    Geometry there = TorsionSpace::setDihedral(start, t, 65.0);
    Geometry back = TorsionSpace::setDihedral(there, t, phi0);

    double max_dev = 0.0;
    for (int a = 0; a < mol.AtomCount(); ++a)
        max_dev = std::max(max_dev,
            (Eigen::Vector3d(back.row(a)) - Eigen::Vector3d(start.row(a))).norm());
    TEST_ASSERT(max_dev < 1e-9, "geometry restored to the original (max deviation "
            + std::to_string(max_dev) + " A)");
}

static void test_cluster_states()
{
    std::cout << "\n=== clusterStates / assignState ===" << std::endl;

    // three well-separated rotamer states
    std::vector<double> angles = { 59.0, 61.0, 62.0, -58.0, -61.0, 178.0, 179.5 };
    auto centres = TorsionSpace::clusterStates(angles, 40.0);
    TEST_ASSERT(centres.size() == 3, "gauche+/gauche-/anti recognised as 3 states (found "
            + std::to_string(centres.size()) + ")");

    // wrap-around: anti scattered across +-180 is ONE state, not two
    std::vector<double> wrapped = { 176.0, 178.0, -179.0, -176.0 };
    auto wrapped_centres = TorsionSpace::clusterStates(wrapped, 40.0);
    TEST_ASSERT(wrapped_centres.size() == 1,
        "values around +-180 form a single state (found " + std::to_string(wrapped_centres.size()) + ")");
    if (wrapped_centres.size() == 1)
        TEST_ASSERT(std::abs(std::abs(wrapped_centres[0]) - 179.75) < 1.0,
            "its centre is near +-180 (got " + std::to_string(wrapped_centres[0]) + ")");

    TEST_ASSERT(TorsionSpace::assignState(58.0, centres) == TorsionSpace::assignState(63.0, centres),
        "two values in the same basin get the same state");
    TEST_ASSERT(TorsionSpace::assignState(58.0, centres) != TorsionSpace::assignState(-58.0, centres),
        "gauche+ and gauche- get different states");
    TEST_ASSERT(TorsionSpace::clusterStates({}, 40.0).empty(), "empty input -> no states");
}

static void test_state_vector()
{
    std::cout << "\n=== stateVector / hammingDistance ===" << std::endl;
    Molecule mol = fromRegistry("butane");
    auto torsions = TorsionSpace::detectTorsions(mol);
    if (torsions.empty()) {
        TEST_ASSERT(false, "no torsion for the state vector");
        return;
    }
    std::vector<std::vector<double>> centres = { { -60.0, 60.0, 180.0 } };

    const Geometry anti = TorsionSpace::setDihedral(mol.getGeometry(), torsions.front(), 180.0);
    const Geometry gauche = TorsionSpace::setDihedral(mol.getGeometry(), torsions.front(), 60.0);

    auto s_anti = TorsionSpace::stateVector(anti, torsions, centres);
    auto s_gauche = TorsionSpace::stateVector(gauche, torsions, centres);

    TEST_ASSERT(s_anti.size() == 1 && s_gauche.size() == 1, "one state per torsion");
    TEST_ASSERT(s_anti[0] == 2, "anti maps onto the 180 deg state (got " + std::to_string(s_anti[0]) + ")");
    TEST_ASSERT(s_gauche[0] == 1, "60 deg maps onto the +60 state (got " + std::to_string(s_gauche[0]) + ")");
    TEST_ASSERT(TorsionSpace::hammingDistance(s_anti, s_gauche) == 1, "Hamming distance 1");
    TEST_ASSERT(TorsionSpace::hammingDistance(s_anti, s_anti) == 0, "Hamming distance to itself is 0");
    TEST_ASSERT(TorsionSpace::hammingDistance(s_anti, { 1, 2 }) == -1, "size mismatch reported as -1");
}

int main()
{
    CurcumaLogger::set_verbosity(0);
    std::cout << "========================================" << std::endl;
    std::cout << "Torsion space unit tests" << std::endl;
    std::cout << "========================================" << std::endl;

    test_detect_butane();
    test_detect_ring();
    test_set_dihedral();
    test_round_trip();
    test_cluster_states();
    test_state_vector();

    std::cout << "\n=== Summary ===" << std::endl;
    std::cout << "Tests run:    " << tests_run << std::endl;
    std::cout << "Tests passed: " << tests_passed << std::endl;
    std::cout << "Tests failed: " << tests_failed << std::endl;
    return tests_failed > 0 ? 1 : 0;
}
