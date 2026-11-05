/*
 * < Periodic MD: barostat cell update, cutoff limit and restart (WP-PERIODIC-NONBONDED N2) >
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
 * Claude Generated (Nov 2026, WP-PERIODIC-NONBONDED N2).
 *
 * (c) Cell update consistency: a GFN-FF instance whose cell and positions are scaled by mu = 1.001 through
 *     GFNFF::updateCell() + UpdateMolecule() must give the energy, gradient and virial of a FRESH instance set up
 *     on the scaled structure with the same Ewald parameter set (alpha, r_cut, k_max fixed by ewald_alpha,
 *     ewald_r_cut, ewald_tolerance; the same absolute D4 radius): |dE| <= 1e-12 Eh, max |dG| <= 1e-12 Eh/Bohr,
 *     |d tr Xi| <= 1e-12 Eh. Round trip back to mu = 1 reproduces the start energy. A cell below the cutoff limit
 *     is refused and leaves the instance unchanged.
 * (d) Limit abort: SimpleMD with a Berendsen barostat driven by an absurd target pressure shrinks the 44 Bohr box
 *     below the HB/XB search radius limit (L_min > 42.43 Bohr): the run stops with StopReason::CellLimit, the
 *     message names the limit, and <basename>.npt_limit.json holds the last valid state.
 * (e) Restart: a crescale run writes a restart file; a new SimpleMD started from it carries the saved cell and
 *     barostat state (events counter) and continues to change the cell.
 * Both MD parts use 20 registry waters packed randomly (rigid, >= 2.6 A between molecules) into L = 44 Bohr:
 * a dilute gas, not a liquid; only the mechanics are tested here.
 */

#include "src/capabilities/simplemd.h"
#include "src/core/curcuma_logger.h"
#include "src/core/energy_calculators/ff_methods/gfnff.h"
#include "src/core/molecule.h"
#include "src/core/units.h"
#include "test_molecule_registry.h"
#include "json.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <random>
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

/// nwat rigid copies of the registry water in a cube of L (Angstrom), every intermolecular atom pair at a
/// minimum-image distance >= min_any (random sequential placement; not a liquid).
static curcuma::Molecule waterBox(int nwat, double L, double min_any, unsigned seed)
{
    const curcuma::Molecule water = TestMolecules::TestMoleculeRegistry::createMolecule("H2O", false);
    const Geometry unit = water.getGeometry().rowwise() - water.getGeometry().colwise().mean();
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    std::vector<Eigen::Vector3d> placed;
    auto mic = [&](Eigen::Vector3d v) {
        for (int k = 0; k < 3; ++k)
            v[k] -= L * std::rint(v[k] / L);
        return v;
    };
    curcuma::Molecule mol;
    for (int n = 0; n < nwat; ++n) {
        for (long tries = 0; tries < 100000; ++tries) {
            Eigen::Vector3d p(u(rng) * L, u(rng) * L, u(rng) * L);
            Eigen::Quaterniond qr(u(rng) - 0.5, u(rng) - 0.5, u(rng) - 0.5, u(rng) - 0.5);
            qr.normalize();
            const Eigen::Matrix3d R = qr.toRotationMatrix();
            std::vector<Eigen::Vector3d> x(unit.rows());
            for (int i = 0; i < unit.rows(); ++i)
                x[i] = p + R * unit.row(i).transpose();
            bool ok = true;
            for (const auto& a : x)
                for (const auto& b : placed)
                    if (mic(a - b).norm() < min_any)
                        ok = false;
            if (!ok)
                continue;
            for (int i = 0; i < unit.rows(); ++i) {
                mol.addPair({ water.Atom(i).first, Position{ x[i](0), x[i](1), x[i](2) } });
                placed.push_back(x[i]);
            }
            break;
        }
    }
    mol.setUnitCell(Eigen::Matrix3d::Identity() * L, true);
    return mol;
}

static json gfnffConfig()
{
    json c = json::object();
    c["verbosity"] = 0;
    c["threads"] = 1;
    c["cache_topology"] = false;
    c["ewald_alpha"] = 0.32;      // fixed parameter set: alpha, r_cut and k_max do not follow L
    c["ewald_r_cut"] = 19.0;
    c["ewald_tolerance"] = 1e-10;
    c["frag_charge_model"] = "reference";
    return c;
}

struct Eval {
    double e = 0.0, tr = 0.0;
    Matrix g;
};

static Eval evalInstance(GFNFF& ff, bool want_virial)
{
    Eval r;
    r.e = ff.Calculation(true);
    r.g = ff.getWorkspace()->gradient();
    r.tr = want_virial ? ff.virial().trace() : 0.0;
    return r;
}

static double maxAbsDiff(const Matrix& a, const Matrix& b)
{
    if (a.rows() != b.rows() || a.cols() != b.cols())
        return 1e30;
    return (a - b).cwiseAbs().maxCoeff();
}

static Mol scaledMol(const Mol& m, double mu)
{
    Mol s = m;
    s.m_geometry *= mu;
    s.m_unit_cell *= mu;
    return s;
}

static void testCellUpdate(const curcuma::Molecule& box)
{
    const Mol m0 = box.getMolInfo();
    const double mu = 1.001;
    GFNFF a(gfnffConfig());
    a.setWantVirial(true);
    check(a.InitialiseMolecule(m0), "(c) instance A set up with the cell");
    a.freezeCutoffsForCellChanges(0.15);
    const GFNFF::D4RadiiBohr frozen = a.dispersionRadiiBohr();
    const double half0 = 0.5 * 44.0;
    std::printf("(c) frozen D4 radius %.6f Bohr (list %.6f), start L/2 %.4f, allowed by 15 %% density gain %.6f\n",
        frozen.eval, frozen.build, half0, half0 * std::cbrt(1.0 / 1.15));
    check(std::abs(frozen.eval - half0 * std::cbrt(1.0 / 1.15)) < 1e-9 && frozen.build - frozen.eval == 10.0,
        "(c) frozen D4 evaluation radius = L/2 (1.15)^(-1/3) (just below), list radius + 10");
    const Eval e0 = evalInstance(a, true);
    // fresh reference at the start geometry with the same absolute D4 radius: equal to A before any update
    json cfg_fixed = gfnffConfig();
    cfg_fixed["dispersion_cutoff_bohr"] = frozen.eval;
    {
        GFNFF r0(cfg_fixed);
        r0.setWantVirial(true);
        check(r0.InitialiseMolecule(m0), "(c) fresh instance at the start cell set up");
        const Eval f0 = evalInstance(r0, true);
        std::printf("(c) start: E_A %.12f, E_fresh %.12f, dE %.2e, max|dG| %.2e, d trXi %.2e\n", e0.e, f0.e, e0.e - f0.e,
            maxAbsDiff(e0.g, f0.g), e0.tr - f0.tr);
        check(std::abs(e0.e - f0.e) <= 1e-12 && maxAbsDiff(e0.g, f0.g) <= 1e-12 && std::abs(e0.tr - f0.tr) <= 1e-12,
            "(c) frozen radius on A == explicit dispersion_cutoff_bohr on a fresh instance (1e-12)");
    }
    // scale A: cell first (checks the cutoffs, rebuilds the Ewald lattice, marks the lists), then the positions
    const Mol m1 = scaledMol(m0, mu);
    const long disp_before = a.dispersionRebuildCount(), rep_before = a.repulsionRebuildCount();
    a.updateCell(m1.m_unit_cell);
    a.UpdateMolecule(m1.m_geometry);
    const Eval e1 = evalInstance(a, true);
    check(a.cellUpdateCount() == 1, "(c) one cell update counted");
    // the repulsion list of 60 atoms holds every pair (not distance-filtered below nb_cell_list_min_atoms), so only
    // the D4 rebuild is counted here
    std::printf("(c) rebuilds after the cell change: D4 %ld -> %ld, repulsion %ld -> %ld\n", disp_before, a.dispersionRebuildCount(),
        rep_before, a.repulsionRebuildCount());
    check(a.dispersionRebuildCount() > disp_before, "(c) D4 list rebuilt after the cell change");
    GFNFF b(cfg_fixed);
    b.setWantVirial(true);
    check(b.InitialiseMolecule(m1), "(c) fresh instance B at the scaled cell set up");
    const Eval e2 = evalInstance(b, true);
    const double v1 = a.cellVolumeBohr3(), v2 = b.cellVolumeBohr3();
    std::printf("(c) mu = %.3f: E_updated %.12f, E_fresh %.12f, dE %.2e Eh, max|dG| %.2e Eh/Bohr, trXi %.9f / %.9f (d %.2e), V %.3f / %.3f Bohr^3\n",
        mu, e1.e, e2.e, e1.e - e2.e, maxAbsDiff(e1.g, e2.g), e1.tr, e2.tr, e1.tr - e2.tr, v1, v2);
    check(std::abs(e1.e - e2.e) <= 1e-12, "(c) energy after updateCell == fresh instance (1e-12 Eh)");
    check(maxAbsDiff(e1.g, e2.g) <= 1e-12, "(c) gradient after updateCell == fresh instance (1e-12 Eh/Bohr)");
    check(std::abs(e1.tr - e2.tr) <= 1e-12, "(c) virial trace after updateCell == fresh instance (1e-12 Eh)");
    check(std::abs(v1 - v2) <= 1e-9 * v2, "(c) cell volume updated");
    check(std::abs(e1.e - e0.e) > 1e-6, "(c) the scaled structure has a different energy (the test is not trivial)");
    // the same with the active Ewald parameters printed: fixed by construction
    const auto pa = a.activeEwaldParameters();
    const auto pb = b.activeEwaldParameters();
    check(pa.alpha == pb.alpha && pa.r_cut == pb.r_cut && pa.k_max == pb.k_max, "(c) identical Ewald parameter set (alpha, r_cut, k_max)");
    // round trip
    a.updateCell(m0.m_unit_cell);
    a.UpdateMolecule(m0.m_geometry);
    const Eval e3 = evalInstance(a, true);
    std::printf("(c) round trip to mu = 1: dE %.2e, max|dG| %.2e, d trXi %.2e\n", e3.e - e0.e, maxAbsDiff(e3.g, e0.g), e3.tr - e0.tr);
    check(std::abs(e3.e - e0.e) <= 1e-12 && maxAbsDiff(e3.g, e0.g) <= 1e-12 && std::abs(e3.tr - e0.tr) <= 1e-12,
        "(c) round trip reproduces the start (1e-12)");
    // a cell below the limit is refused and nothing changes
    bool threw = false;
    std::string why;
    try {
        a.updateCell(m0.m_unit_cell * 0.9);  // L = 39.6 Bohr < 42.43
    } catch (const std::exception& e) {
        threw = true;
        why = e.what();
    }
    std::printf("(c) refused cell: %s\n", why.c_str());
    check(threw, "(c) a cell below the cutoff limit is refused");
    a.UpdateMolecule(m0.m_geometry);
    const Eval e4 = evalInstance(a, true);
    check(a.cellUpdateCount() == 2 && std::abs(e4.e - e0.e) <= 1e-12, "(c) the refused update left the instance unchanged");
    // 4 threads: the updated instance equals the single-thread one
    {
        json c4 = gfnffConfig();
        c4["threads"] = 4;
        GFNFF a4(c4);
        a4.setWantVirial(true);
        a4.InitialiseMolecule(m0);
        a4.freezeCutoffsForCellChanges(0.15);
        a4.Calculation(true);
        a4.updateCell(m1.m_unit_cell);
        a4.UpdateMolecule(m1.m_geometry);
        const Eval e5 = evalInstance(a4, true);
        std::printf("(c) 4 threads after the update: dE %.2e, max|dG| %.2e, d trXi %.2e vs 1 thread\n", e5.e - e1.e, maxAbsDiff(e5.g, e1.g), e5.tr - e1.tr);
        check(std::abs(e5.e - e1.e) <= 1e-11 && maxAbsDiff(e5.g, e1.g) <= 1e-11, "(c) 4 threads == 1 thread after the update (1e-11)");
    }
}

static json mdController(const std::string& barostat, double p0, double tau_p, int interval)
{
    json controller;
    json& md = controller["simplemd"];
    md["method"] = "gfnff";
    md["no_restart"] = true;
    md["thermostat"] = "csvr";
    md["coupling"] = 100.0;
    md["temperature"] = 300.0;
    md["seed"] = 7;
    md["dt"] = 0.5;
    md["max_time"] = 50.0;
    md["no_center"] = true;
    md["print_frequency"] = 1000000;
    md["dump_frequency"] = 10;
    md["write_restart_frequency"] = 1000000;
    md["write_unique"] = false;
    md["verbosity"] = 0;
    md["barostat"] = barostat;
    md["pressure"] = p0;
    md["tau_p"] = tau_p;
    md["compressibility"] = 4.5e-5;
    md["pcoupl_interval"] = interval;
    controller["gfnff"]["frag_charge_model"] = "reference";
    controller["gfnff"]["cache_topology"] = false;
    return controller;
}

static bool fileExists(const std::string& name)
{
    std::ifstream f(name);
    return f.good();
}

static void testLimitAbort(const curcuma::Molecule& box)
{
    // Berendsen: mu^3 = 1 - (dt/tau) beta (P0 - P) = 1 - 0.005 * 4.5e-5 * 2e5 = 0.955 per step (dt 0.5, tau 100,
    // every step): L shrinks 1.5 % per step, the limit L_min > 42.43 Bohr (HB/XB search radius 21.21) is reached
    // at the third coupling (44 * 0.985^3 = 42.05).
    json controller = mdController("berendsen", 2.0e5, 100.0, 1);
    SimpleMD sim(controller, true);
    sim.overrideBasename("pbc_npt_limit");
    sim.setMolecule(box);
    check(sim.Initialise(), "(d) set up (berendsen, 2e5 bar)");
    sim.prepareRun();
    const double l0 = sim.unitCell()(0, 0);
    int steps = 0;
    while (sim.step())
        ++steps;
    const double l1 = sim.unitCell()(0, 0);
    const std::string restart_name = "pbc_npt_limit.snapshots/pbc_npt_limit.npt_limit.json";  // snapshotPath() without BMT
    std::printf("(d) stopped after %d steps, reason %d (CellLimit = %d), events %ld, L %.4f -> %.4f A (%.4f Bohr), restart file %s\n",
        steps, static_cast<int>(sim.stopReason()), static_cast<int>(SimpleMD::StopReason::CellLimit), sim.barostatEvents(), l0, l1,
        l1 * CurcumaUnit::Length::ANGSTROM_TO_BOHR, fileExists(restart_name) ? "written" : "MISSING");
    check(sim.stopReason() == SimpleMD::StopReason::CellLimit, "(d) run stopped with StopReason::CellLimit");
    check(steps >= 2 && steps <= 4, "(d) the limit was hit at the expected coupling (steps 2..4)");
    check(l1 * CurcumaUnit::Length::ANGSTROM_TO_BOHR > 42.43, "(d) the kept cell is the last valid one (above the limit)");
    check(fileExists(restart_name), "(d) restart file of the last valid state written");
    if (fileExists(restart_name)) {
        std::ifstream f(restart_name);
        json r;
        f >> r;
        const bool has = r.contains("MD") && r["MD"].contains("cell") && r["MD"].contains("barostat");
        check(has, "(d) the restart file carries the cell and the barostat state");
        if (has)
            check(std::abs(r["MD"]["cell"][0].get<double>() - l1) < 1e-9, "(d) the saved cell is the kept one");
    }
    sim.finalizeRun();
}

static void testRestart(const curcuma::Molecule& box)
{
    const std::string file = "pbc_npt_restart_state.json";
    double l_saved = 0.0;
    long events_saved = 0;
    {
        json controller = mdController("crescale", 1.01325, 100.0, 2);
        SimpleMD sim(controller, true);
        sim.overrideBasename("pbc_npt_restart_a");
        sim.setMolecule(box);
        check(sim.Initialise(), "(e) set up (crescale)");
        sim.prepareRun();
        for (int s = 0; s < 40; ++s)
            sim.step();
        l_saved = sim.unitCell()(0, 0);
        events_saved = sim.barostatEvents();
        std::printf("(e) run A: 40 steps, %ld cell changes, L %.6f A, last mu %.6f, P %.2f bar, rho %.5f\n", events_saved, l_saved,
            sim.lastCellScale(), sim.pressureBar(), sim.densityGcm3());
        check(events_saved == 20 && std::abs(l_saved - 44.0 * CurcumaUnit::Length::BOHR_TO_ANGSTROM) > 1e-6, "(e) crescale changed the cell every 2nd step");
        json restart;
        restart["MD"] = sim.restartState();
        std::ofstream out(file);
        out << restart << std::endl;
        sim.finalizeRun();
    }
    {
        json controller = mdController("none", 1.01325, 100.0, 2);  // the file decides: barostat and cell come from it
        controller["simplemd"]["restart_file"] = file;
        controller["simplemd"]["no_restart"] = false;
        SimpleMD sim(controller, true);
        sim.overrideBasename("pbc_npt_restart_b");
        sim.setMolecule(box);  // start cell 44 Bohr; the restart overrides it
        check(sim.Initialise(), "(e) set up from the restart file");
        const double l_loaded = sim.unitCell()(0, 0);
        std::printf("(e) run B: loaded L %.6f A (saved %.6f), events %ld (saved %ld)\n", l_loaded, l_saved, sim.barostatEvents(), events_saved);
        check(std::abs(l_loaded - l_saved) < 1e-9, "(e) the cell of the restart file is used");
        check(sim.barostatEvents() == events_saved, "(e) the barostat state is restored");
        sim.prepareRun();
        for (int s = 0; s < 20; ++s)
            sim.step();
        std::printf("(e) run B: +20 steps, events %ld, L %.6f A, P %.2f bar\n", sim.barostatEvents(), sim.unitCell()(0, 0), sim.pressureBar());
        check(sim.barostatEvents() == events_saved + 10 && std::abs(sim.unitCell()(0, 0) - l_loaded) > 1e-9,
            "(e) the barostat continues from the restart (10 more cell changes, cell moved)");
        check(sim.stopReason() == SimpleMD::StopReason::Running, "(e) the continued run is healthy");
        sim.finalizeRun();
    }
}

int main()
{
    CurcumaLogger::set_verbosity(0);
    const double L = 44.0 * CurcumaUnit::Length::BOHR_TO_ANGSTROM; // 23.28 A
    const curcuma::Molecule box = waterBox(20, L, 2.6, 11);
    std::printf("box: %d atoms in L = %.4f A (%.1f Bohr)\n", box.AtomCount(), L, 44.0);
    check(box.AtomCount() == 60, "20 waters packed");
    testCellUpdate(box);
    testLimitAbort(box);
    testRestart(box);
    std::printf("%s: %d failure(s)\n", g_failures ? "FAILED" : "OK", g_failures);
    return g_failures ? 1 : 0;
}
