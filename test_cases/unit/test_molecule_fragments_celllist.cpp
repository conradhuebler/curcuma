/*
 * < Molecule::GetFragments: cell-grid search equals the all-pairs scan >
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * Claude Generated (Nov 2026, WP-PERIODIC-NONBONDED M1).
 *
 * GetFragments() runs a cell-grid breadth-first search from Molecule::fragmentCellListMinAtoms() atoms on
 * (SimpleMD calls it every step of a periodic run). Its result must be the all-pairs scan bit for bit:
 * the same fragments, the same order of fragments (multimap keyed by a mass sum that depends on the join
 * order) and the same FragmentMass() values. Systems (built from registry molecules, positions randomised
 * by a fixed-seed engine, no geometry literals):
 *   (a) water boxes, cubic cell, liquid density, molecules torn across faces, atoms outside the cell,
 *       whole molecules shifted by lattice vectors; 10 seeds;
 *   (b) the same coordinates without a cell (Cartesian grid);
 *   (c) a triclinic cell (rounding-based minimum image of Molecule::CalculateDistancePBC);
 *   (d) a dense random gas of caffeine atoms (chains merge, long breadth-first join orders), with and
 *       without a cell;
 *   (e) a cell too small for 3 bins per axis (falls back to the scan) and an exploded, very sparse cluster
 *       (bin count cap);
 *   (f) a timing line for 2910 atoms (informative, not a gate).
 */

#include "src/core/elements.h"
#include "src/core/molecule.h"
#include "test_molecule_registry.h"

#include <chrono>
#include <climits>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

#include "generated/parameter_registry.h"
static struct RegistryInitializer {
    RegistryInitializer() { initialize_generated_registry(); }
} registry_initializer;

static int g_failures = 0;
static void check(bool ok, const std::string& what)
{
    std::printf("%s %s\n", ok ? "[PASS]" : "[FAIL]", what.c_str());
    if (!ok)
        ++g_failures;
}

struct FragResult {
    std::vector<std::vector<int>> fragments;
    std::vector<double> masses;
};

static FragResult fragmentsWith(curcuma::Molecule& mol, int min_atoms)
{
    curcuma::Molecule::setFragmentCellListMinAtoms(min_atoms);
    mol.setGeometry(mol.getGeometry()); // invalidates the fragment cache
    FragResult r;
    r.fragments = mol.GetFragments();
    r.masses = mol.FragmentMass();
    return r;
}

// returns true when scan and grid agree bit for bit; prints the fragment count
static bool compare(curcuma::Molecule& mol, const std::string& label)
{
    const FragResult scan = fragmentsWith(mol, INT_MAX);
    const FragResult grid = fragmentsWith(mol, 0);
    bool same = scan.fragments == grid.fragments && scan.masses.size() == grid.masses.size();
    for (std::size_t i = 0; same && i < scan.masses.size(); ++i)
        same = scan.masses[i] == grid.masses[i];
    check(same, label + ": " + std::to_string(mol.AtomCount()) + " atoms, " + std::to_string(scan.fragments.size()) + " fragments, grid == scan");
    return same;
}

// n waters (registry H2O) at random positions and orientations in [0, L)^3, then perturbed: a share of
// the molecules shifted by a lattice vector, a share of single atoms shifted by a lattice vector (torn)
static curcuma::Molecule waterBox(int n, double L, unsigned seed, bool cell, bool tear)
{
    const curcuma::Molecule w = TestMolecules::TestMoleculeRegistry::createMolecule("H2O", false);
    const Geometry g = w.getGeometry();
    Eigen::Vector3d c = Eigen::Vector3d::Zero();
    for (int k = 0; k < g.rows(); ++k)
        c += g.row(k).transpose();
    c /= g.rows();
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> u01(0.0, 1.0);
    std::normal_distribution<double> nrm(0.0, 1.0);
    curcuma::Molecule mol;
    for (int m = 0; m < n; ++m) {
        const Eigen::Quaterniond q = Eigen::Quaterniond(nrm(rng), nrm(rng), nrm(rng), nrm(rng)).normalized();
        const Eigen::Vector3d centre(L * u01(rng), L * u01(rng), L * u01(rng));
        Eigen::Vector3d shift = Eigen::Vector3d::Zero();
        if (tear && u01(rng) < 0.1)
            shift = L * Eigen::Vector3d(std::floor(3 * u01(rng)) - 1, std::floor(3 * u01(rng)) - 1, std::floor(3 * u01(rng)) - 1);
        for (int k = 0; k < g.rows(); ++k) {
            Eigen::Vector3d p = q * (g.row(k).transpose() - c) + centre + shift;
            if (tear && u01(rng) < 0.05)
                p(int(3 * u01(rng)) % 3) += (u01(rng) < 0.5 ? -L : L);
            mol.addPair({ w.Atom(k).first, Position(p) });
        }
    }
    if (cell)
        mol.setUnitCell(Eigen::Matrix3d::Identity() * L, true);
    return mol;
}

// atoms of caffeine (registry) repeated, every atom at a random position in [0, L)^3
static curcuma::Molecule caffeineGas(int copies, double L, unsigned seed, bool cell)
{
    const curcuma::Molecule caf = TestMolecules::TestMoleculeRegistry::createMolecule("caffeine", false);
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> u01(0.0, 1.0);
    curcuma::Molecule mol;
    for (int m = 0; m < copies; ++m)
        for (int k = 0; k < caf.AtomCount(); ++k)
            mol.addPair({ caf.Atom(k).first, Position(L * u01(rng), L * u01(rng), L * u01(rng)) });
    if (cell)
        mol.setUnitCell(Eigen::Matrix3d::Identity() * L, true);
    return mol;
}

int main()
{
    const int default_min = curcuma::Molecule::fragmentCellListMinAtoms();
    check(default_min == 256, "default switch-over at 256 atoms (got " + std::to_string(default_min) + ")");

    // (a), (b): 300 waters at liquid density (L = 20.8 A: 0.1 atoms/A^3)
    for (unsigned seed = 1; seed <= 10; ++seed) {
        curcuma::Molecule box = waterBox(300, 20.8, seed, true, true);
        compare(box, "(a) water box, cell, torn, seed " + std::to_string(seed));
        curcuma::Molecule nocell = waterBox(300, 20.8, seed, false, true);
        compare(nocell, "(b) same coordinates, no cell, seed " + std::to_string(seed));
    }
    // (c) triclinic cell (rows = lattice vectors, as setUnitCell stores them)
    {
        curcuma::Molecule box = waterBox(300, 20.8, 11, false, false);
        Eigen::Matrix3d tri;
        tri << 20.8, 0.0, 0.0,
            4.0, 20.0, 0.0,
            -3.0, 2.5, 21.5;
        box.setUnitCell(tri, true);
        compare(box, "(c) triclinic cell");
        box.setUnitCell(tri.transpose(), true);
        compare(box, "(c) triclinic cell, transposed");
    }
    // (d) dense random gas: chains of bonded atoms merge into large fragments
    for (unsigned seed = 1; seed <= 3; ++seed) {
        curcuma::Molecule gas = caffeineGas(20, 14.0, seed, true);
        compare(gas, "(d) caffeine-atom gas, cell, seed " + std::to_string(seed));
        curcuma::Molecule gas_nc = caffeineGas(20, 14.0, seed, false);
        compare(gas_nc, "(d) caffeine-atom gas, no cell, seed " + std::to_string(seed));
    }
    // (e) small cell: fewer than 3 bins on an axis -> scan; sparse exploded cluster -> bin cap
    {
        curcuma::Molecule small = waterBox(20, 5.0, 3, true, true);
        compare(small, "(e) 5 A cell (scan fallback)");
        curcuma::Molecule sparse = waterBox(300, 5000.0, 4, false, false);
        compare(sparse, "(e) 300 waters spread over 5000 A, no cell (bin cap)");
        curcuma::Molecule mixed = waterBox(300, 20.8, 5, false, false);
        Geometry g = mixed.getGeometry();
        g.row(0) += Eigen::RowVector3d(1e5, -2e5, 3e5); // one atom very far away
        mixed.setGeometry(g);
        compare(mixed, "(e) one atom 3.7e5 A away, no cell");
    }
    // (f) timing at the MD box size (informative)
    {
        curcuma::Molecule box = waterBox(970, 30.692278, 7, true, true);
        const auto t0 = std::chrono::steady_clock::now();
        const FragResult scan = fragmentsWith(box, INT_MAX);
        const auto t1 = std::chrono::steady_clock::now();
        const FragResult grid = fragmentsWith(box, 0);
        const auto t2 = std::chrono::steady_clock::now();
        check(scan.fragments == grid.fragments && scan.masses == grid.masses, "(f) 970 waters, cell, grid == scan");
        std::printf("      timing 2910 atoms: scan %.2f ms, grid %.2f ms\n",
            std::chrono::duration<double, std::milli>(t1 - t0).count(),
            std::chrono::duration<double, std::milli>(t2 - t1).count());
    }
    curcuma::Molecule::setFragmentCellListMinAtoms(default_min);

    std::printf("%s: %d failure(s)\n", g_failures ? "FAILED" : "PASSED", g_failures);
    return g_failures ? 1 : 0;
}
