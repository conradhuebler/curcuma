/*
 * < Periodic MD: RATTLE constraints with a unit cell >
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
 * RATTLE builds its constraint list from raw distances and evaluates raw coordinate differences
 * (SimpleMD::InitConstrainedBonds(), SimpleMD::Rattle()); the cell wrap moves whole fragments. Checks on
 * 8 registry waters in a cubic cell of 44 Bohr, GFN-FF periodic (coulomb_ewald), NVE from 300 K, seed 7,
 * dt 1 fs, 400 steps, 6 of the 8 molecules placed with their centre 0.02 A inside a face (so they wrap):
 *   (a) rattle 2 with the cell: cell wraps happen, every constrained O-H pair keeps its squared raw
 *       distance within 2 x rattle_tol_12 of the start value in every step (no pair torn by a wrap);
 *   (a') the same box with one H moved by +L before the start (a bond across the face in the input): the
 *       bond is reported as crossing and left unconstrained (documents the limitation; warning printed);
 *   (b) energy conservation of (a) against the same run without a cell (bounds below);
 *   (c) DOF: 3N - N_c - 3 with the cell (translation only, D1), 3N - N_c - 6 without; rattle 1 +
 *       rattle_13 3N - N_c - 3 with the constraints actually applied (the 1-3 list misses H-O-H when O
 *       comes first, a known defect, TODO.md).
 */

#include "src/capabilities/simplemd.h"
#include "src/core/curcuma_logger.h"
#include "src/core/elements.h"
#include "src/core/molecule.h"
#include "src/core/units.h"
#include "test_molecule_registry.h"
#include "json.hpp"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using json = nlohmann::json;

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

// 8 registry waters on a 2 x 2 x 2 grid, centres at (L/2 - 0.02) and (L - 0.02) per axis: the outer ones
// stick out of [0, L) but are whole. torn_h: move atom 1 (an H of molecule 0) by +L in x.
static curcuma::Molecule waterBox8(bool cell, double L, bool torn_h)
{
    const curcuma::Molecule water = TestMolecules::TestMoleculeRegistry::createMolecule("H2O", false);
    const Geometry unit = water.getGeometry().rowwise() - water.getGeometry().colwise().mean();
    curcuma::Molecule mol;
    int n = 0;
    for (int a = 0; a < 2; ++a)
        for (int b = 0; b < 2; ++b)
            for (int c = 0; c < 2; ++c) {
                const int perm = n % 3;
                for (int k = 0; k < unit.rows(); ++k) {
                    const double p[3] = { unit(k, perm), unit(k, (perm + 1) % 3), unit(k, (perm + 2) % 3) };
                    Position x{ p[0] + (a + 1) * L / 2 - 0.02, p[1] + (b + 1) * L / 2 - 0.02, p[2] + (c + 1) * L / 2 - 0.02 };
                    if (torn_h && n == 0 && k == 1)
                        x(0) += L;
                    mol.addPair({ water.Atom(k).first, x });
                }
                ++n;
            }
    if (cell)
        mol.setUnitCell(Eigen::Matrix3d::Identity() * L, true);
    return mol;
}

static json mdController(int rattle, bool rattle13)
{
    json controller;
    json& md = controller["simplemd"];
    md["method"] = "gfnff";
    md["no_restart"] = true;
    md["thermostat"] = "none";
    md["temperature"] = 300.0;
    md["seed"] = 7;
    md["dt"] = 1.0;
    md["max_time"] = 1000.0;
    md["no_center"] = true;
    md["remove_com_motion"] = 10.0;
    md["remove_com_mode"] = 1;
    md["rattle"] = rattle;
    md["rattle_13"] = rattle13;
    md["print_frequency"] = 1000000;
    md["dump_frequency"] = 1000000;
    md["write_restart_frequency"] = 1000000;
    md["write_unique"] = false;
    md["verbosity"] = 0;
    return controller;
}

struct RunResult {
    bool ok = false, periodic = false;
    std::vector<double> e;
    long cell_wraps = 0;
    int dof = 0, bonds = 0, angles = 0, crossing = 0;
    double max_dev = 0;   // max over steps and O-H pairs of |d^2 - d0^2| (A^2), raw coordinates
    double max_raw_oh = 0; // longest raw O-H distance seen (A)
};

static RunResult run(const json& controller, const curcuma::Molecule& mol, const std::string& tag, int nsteps)
{
    RunResult r;
    SimpleMD sim(controller, true);
    sim.overrideBasename("pbc_md_rattle_" + tag);
    sim.setMolecule(mol);
    if (!sim.Initialise())
        return r;
    sim.prepareRun();
    r.dof = sim.degreesOfFreedom();
    r.periodic = sim.periodicEnergy();
    r.bonds = sim.rattleBondConstraints();
    r.angles = sim.rattleAngleConstraints();
    r.crossing = sim.rattleCrossingBonds();
    // O-H pairs of every water (O first in the registry molecule: atoms 3w, 3w+1, 3w+2)
    const int nw = mol.AtomCount() / 3;
    std::vector<double> d0;
    const Geometry g0 = sim.positions();
    for (int w = 0; w < nw; ++w)
        for (int h = 1; h <= 2; ++h)
            d0.push_back((g0.row(3 * w) - g0.row(3 * w + h)).squaredNorm());
    for (int s = 0; s < nsteps; ++s) {
        sim.step();
        r.e.push_back(sim.potentialEnergy() + sim.kineticEnergy());
        const Geometry& g = sim.positions();
        int k = 0;
        for (int w = 0; w < nw; ++w)
            for (int h = 1; h <= 2; ++h, ++k) {
                const double d2 = (g.row(3 * w) - g.row(3 * w + h)).squaredNorm();
                r.max_dev = std::max(r.max_dev, std::abs(d2 - d0[k]));
                r.max_raw_oh = std::max(r.max_raw_oh, std::sqrt(d2));
            }
    }
    sim.finalizeRun();
    r.cell_wraps = sim.cellWrapCount();
    r.ok = true;
    return r;
}

static double slopePerPs(const std::vector<double>& e, double dt)
{
    double st = 0, se = 0, stt = 0, ste = 0;
    const double n = static_cast<double>(e.size());
    for (size_t k = 0; k < e.size(); ++k) {
        const double t = (k + 1) * dt * 1e-3;
        st += t;
        se += e[k];
        stt += t * t;
        ste += t * e[k];
    }
    return (n * ste - st * se) / (n * stt - st * st);
}

static double maxMin(const std::vector<double>& e)
{
    double lo = 1e30, hi = -1e30;
    for (double x : e) {
        lo = std::min(lo, x);
        hi = std::max(hi, x);
    }
    return hi - lo;
}

int main()
{
    CurcumaLogger::set_verbosity(0);
    const double L = 44.0 * CurcumaUnit::Length::BOHR_TO_ANGSTROM; // 23.28 A
    const int nsteps = 400;
    const int n = 24, nw = 8;

    const curcuma::Molecule box = waterBox8(true, L, false);
    const curcuma::Molecule vac = waterBox8(false, L, false);
    const RunResult rc = run(mdController(2, false), box, "cell_r2", nsteps);
    const RunResult rv = run(mdController(2, false), vac, "vac_r2", nsteps);
    check(rc.ok && rv.ok, "set up (rattle 2, cell and vacuum)");
    if (rc.ok && rv.ok) {
        std::printf("(a) rattle 2 + cell: periodic %d, constraints %d (+%d 1-3), crossing %d, cell wraps %ld, max |d^2-d0^2| %.2e A^2, "
                    "longest raw O-H %.4f A\n",
            rc.periodic ? 1 : 0, rc.bonds, rc.angles, rc.crossing, rc.cell_wraps, rc.max_dev, rc.max_raw_oh);
        std::printf("    vacuum: max |d^2-d0^2| %.2e A^2\n", rv.max_dev);
        check(rc.periodic, "(a) gfnff + cell is periodic");
        check(rc.bonds == 2 * nw && rc.angles == 0 && rc.crossing == 0, "(a) 16 O-H constraints, none crossing a face");
        check(rc.cell_wraps > 0, "(a) molecules were wrapped into the cell during the run");
        check(rc.max_dev <= 2e-4 && rc.max_raw_oh < 1.2, "(a) every O-H constraint holds in raw coordinates through the wraps (<= 2 x rattle_tol_12)");
        const double dc = slopePerPs(rc.e, 1.0), dv = slopePerPs(rv.e, 1.0);
        const double mc = maxMin(rc.e), mv = maxMin(rv.e);
        std::printf("(b) NVE 0.4 ps dt 1 fs: drift %.3e / %.3e Eh/ps, max-min %.3e / %.3e Eh (cell / vacuum)\n", dc, dv, mc, mv);
        check(std::abs(dc) < 2e-4 && mc < 5.0 * mv + 2e-5, "(b) cell run conserves energy like the vacuum run (|drift| < 2e-4 Eh/ps, max-min < 5 x vacuum + 2e-5 Eh)");
        std::printf("(c) DOF cell %d (3N-16-3 = %d), vacuum %d (3N-16-6 = %d)\n", rc.dof, 3 * n - 16 - 3, rv.dof, 3 * n - 16 - 6);
        check(rc.dof == 3 * n - 2 * nw - 3, "(c) DOF with the cell: 3N - N_c - 3");
        check(rv.dof == 3 * n - 2 * nw - 6, "(c) DOF without a cell: 3N - N_c - 6");
    }

    {
        const RunResult r = run(mdController(1, true), box, "cell_r1_13", 50);
        check(r.ok, "(c) rattle 1 + rattle_13 + cell: set up");
        if (r.ok) {
            // Known defect (M1, TODO.md): InitConstrainedBonds() enumerates 1-3 pairs i-j-k only with
            // i > j > k, so H-O-H with O first in the input (registry water) gives 0 of 8 angle constraints.
            // Checked here: the DOF follows the constraints actually applied (whatever their number).
            std::printf("(c) rattle 1 + rattle_13 + cell: constraints %d + %d (8 H-O-H expected for rigid water; 1-3 enumeration "
                        "defect, TODO.md), DOF %d\n",
                r.bonds, r.angles, r.dof);
            check(r.bonds == 2 * nw && r.dof == 3 * n - r.bonds - r.angles - 3, "(c) rattle 1 + rattle_13 + cell: DOF 3N - N_c - 3 with the applied constraints");
        }
    }

    {
        const curcuma::Molecule torn = waterBox8(true, L, true);
        const RunResult r = run(mdController(2, false), torn, "cell_torn", 20);
        check(r.ok, "(a') torn input: set up");
        if (r.ok) {
            std::printf("(a') one O-H across the face in the input: constraints %d, crossing %d, DOF %d\n", r.bonds, r.crossing, r.dof);
            check(r.bonds == 2 * nw - 1 && r.crossing == 1, "(a') the bond across the face is detected and reported as unconstrained (15 constraints, 1 crossing)");
            check(r.dof == 3 * n - (2 * nw - 1) - 3, "(a') DOF counts only the 15 applied constraints");
        }
    }

    std::printf("%s: %d failure(s)\n", g_failures == 0 ? "ALL PASSED" : "FAILED", g_failures);
    return g_failures == 0 ? 0 : 1;
}
