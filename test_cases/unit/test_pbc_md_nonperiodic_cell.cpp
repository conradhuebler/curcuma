/*
 * < Periodic MD: a unit cell with a method whose energy is not periodic >
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
 * Claude Generated (Nov 2026, WP-PERIODIC-NONBONDED A6b).
 *
 * With a unit cell and a method whose energy does not read the cell (hasPeriodicEnergy() false: uff,
 * eht, GFN-FF with coulomb_ewald false, ...) SimpleMD no longer wraps into the cell and removes
 * COM translation AND rotation as in vacuum (-6 DOF). Since WP-PERIODIC-XTB P5 (Nov 2026) the native
 * gfn1/gfn2 with a cell are periodic energies (check (d)); (a) therefore uses eht instead of gfn2. Checks:
 *   (a) uff and eht, registry water dimer around the origin, cell L = 44 Bohr ([0, L), so the dimer
 *       starts outside it), wall_potential pbc with the box [-L/2, L/2), fresh 300 K start, seed 7,
 *       NVE 200 steps of 0.5 fs, COM/rotation removal every 10 fs: 0 cell wraps, 0 container actions
 *       (the dimer stays in the box), DOF 3N-6, and the E_tot series bitwise equal to the same run
 *       without a cell (so the drift is the vacuum reference by construction; both are printed);
 *   (b) uff, box = cell, the dimer drifting through the +x face: the container wrap acts (it really
 *       leaves the box), the cell wrap does not;
 *   (c) the other corner of the matrix: gfnff + cell (periodic) has DOF 3N-3, gfnff + cell +
 *       coulomb_ewald false has DOF 3N-6;
 *   (d) gfn2 + cell (WP-PERIODIC-XTB): periodic energy, DOF 3N-3.
 */

#include "src/capabilities/simplemd.h"
#include "src/core/curcuma_logger.h"
#include "src/core/elements.h"
#include "src/core/molecule.h"
#include "src/core/units.h"
#include "src/tools/general.h"
#include "test_molecule_registry.h"
#include "json.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
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

// Registry water dimer (Angstrom), its centre of mass moved to `centre`; optional cubic cell of edge L.
static curcuma::Molecule dimerAt(const Position& centre, bool cell, double L)
{
    const curcuma::Molecule dimer = TestMolecules::TestMoleculeRegistry::createMolecule("H2O_dimer", false);
    const Geometry g = dimer.getGeometry();
    double mass = 0;
    Position com = Position::Zero();
    for (int k = 0; k < g.rows(); ++k) {
        const double m = Elements::AtomicMass[dimer.Atom(k).first];
        mass += m;
        com += m * g.row(k).transpose();
    }
    com /= mass;
    curcuma::Molecule mol;
    for (int k = 0; k < g.rows(); ++k)
        mol.addPair({ dimer.Atom(k).first, Position(g.row(k).transpose() - com + centre) });
    if (cell)
        mol.setUnitCell(Eigen::Matrix3d::Identity() * L, true);
    return mol;
}

static json mdController(const std::string& method, double dt, double max_time, double box_min, double box_max,
    const std::string& restart_file = "", bool ewald = true)
{
    json controller;
    json& md = controller["simplemd"];
    md["method"] = method;
    if (!restart_file.empty())
        md["restart_file"] = restart_file;
    md["no_restart"] = restart_file.empty();
    md["thermostat"] = "none";
    md["temperature"] = 300.0;
    md["seed"] = 7;
    md["dt"] = dt;
    md["max_time"] = max_time;
    md["no_center"] = true;
    md["remove_com_motion"] = 10.0;
    md["remove_com_mode"] = 1;
    md["rattle"] = 0;
    md["print_frequency"] = 1000000;
    md["dump_frequency"] = 1000000;
    md["write_restart_frequency"] = 1000000;
    md["write_unique"] = false;
    md["verbosity"] = 0;
    md["wall_type"] = "rect";
    md["wall_potential"] = "pbc";
    md["wall_x_min"] = box_min;
    md["wall_y_min"] = box_min;
    md["wall_z_min"] = box_min;
    md["wall_x_max"] = box_max;
    md["wall_y_max"] = box_max;
    md["wall_z_max"] = box_max;
    if (!ewald)
        controller["gfnff"]["coulomb_ewald"] = false;
    return controller;
}

struct RunResult {
    bool ok = false;
    std::vector<double> e;
    long cell_wraps = 0, container = 0;
    int dof = 0;
    bool periodic = false;
};

static RunResult run(const json& controller, const curcuma::Molecule& mol, const std::string& tag, int nsteps)
{
    RunResult r;
    SimpleMD sim(controller, true);
    sim.overrideBasename("pbc_md_nonperiodic_cell_" + tag);
    sim.setMolecule(mol);
    if (!sim.Initialise())
        return r;
    sim.prepareRun();
    r.dof = sim.degreesOfFreedom();
    r.periodic = sim.periodicEnergy();
    for (int s = 0; s < nsteps; ++s) {
        sim.step();
        r.e.push_back(sim.potentialEnergy() + sim.kineticEnergy());
    }
    sim.finalizeRun();
    r.cell_wraps = sim.cellWrapCount();
    r.container = sim.containerWrapCount() + sim.containerRefusedCount() + sim.containerReflectionCount();
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
    const double dt = 0.5;
    const int nsteps = 200;
    const double max_time = nsteps * dt + 10;

    // (a) cell [0, L), box [-L/2, L/2), dimer at the origin: before A6b cell and container wrap fought (+L/-L per step)
    for (const std::string method : { std::string("uff"), std::string("eht") }) {
        const curcuma::Molecule with_cell = dimerAt(Position::Zero(), true, L);
        const curcuma::Molecule no_cell = dimerAt(Position::Zero(), false, L);
        const RunResult rc = run(mdController(method, dt, max_time, -L / 2, L / 2), with_cell, method + "_cell", nsteps);
        const RunResult rv = run(mdController(method, dt, max_time, -L / 2, L / 2), no_cell, method + "_vac", nsteps);
        check(rc.ok && rv.ok, "(a) " + method + ": both runs set up");
        if (!rc.ok || !rv.ok)
            continue;
        bool identical = rc.e.size() == rv.e.size();
        for (size_t k = 0; identical && k < rc.e.size(); ++k)
            identical = (rc.e[k] == rv.e[k]);
        std::printf("(a) %s: periodic %d, DOF %d / %d, cell wraps %ld, container actions %ld / %ld; drift %.3e / %.3e Eh/ps, "
                    "max-min %.3e / %.3e Eh (cell / vacuum)\n",
            method.c_str(), rc.periodic ? 1 : 0, rc.dof, rv.dof, rc.cell_wraps, rc.container, rv.container,
            slopePerPs(rc.e, dt), slopePerPs(rv.e, dt), maxMin(rc.e), maxMin(rv.e));
        check(!rc.periodic, "(a) " + method + " + cell: energy not periodic");
        check(rc.cell_wraps == 0, "(a) " + method + " + cell: no cell wrap");
        check(rc.container == 0, "(a) " + method + " + cell: no container action (the dimer stays in the box)");
        check(rc.dof == 3 * with_cell.AtomCount() - 6 && rv.dof == rc.dof, "(a) " + method + " + cell: DOF 3N-6 as in vacuum");
        check(identical, "(a) " + method + ": E_tot series bitwise equal to the run without a cell");
    }

    // (b) box = cell, the dimer leaves through the +x face: the container wrap acts, the cell wrap does not
    {
        const curcuma::Molecule mol = dimerAt(Position(L - 0.5, L / 2, L / 2), true, L);
        Geometry v = Geometry::Zero(mol.AtomCount(), 3);
        for (int i = 0; i < mol.AtomCount(); ++i)
            v(i, 0) = 0.04; // ~0.02 A/fs: 2 A in 100 fs
        json s;
        s["method"] = "uff";
        s["thermostat"] = "none";
        s["dT"] = dt;
        s["MaxTime"] = max_time;
        s["T"] = 300.0;
        s["currentStep"] = 0.0;
        s["seed"] = 1;
        s["rmrottrans"] = 0;
        s["nocenter"] = true;
        s["COM"] = false;
        s["geometry"] = Tools::Geometry2String(mol.getGeometry());
        s["velocities"] = Tools::Geometry2String(v);
        {
            std::ofstream f("pbc_md_nonperiodic_cell_b.json");
            f << s.dump(2);
        }
        json c = mdController("uff", dt, max_time, 0.0, L, "pbc_md_nonperiodic_cell_b.json");
        c["simplemd"]["remove_com_motion"] = 0.0;
        const RunResult r = run(c, mol, "b", nsteps);
        check(r.ok, "(b) set up");
        if (r.ok) {
            std::printf("(b) uff, box = cell, leaving through +x: cell wraps %ld, container actions %ld\n", r.cell_wraps, r.container);
            check(r.cell_wraps == 0, "(b) no cell wrap");
            check(r.container > 0, "(b) the container wrap acts when the dimer leaves the box");
        }
    }

    // (c) DOF of the periodic corner and of gfnff without the periodic Coulomb
    for (const bool ewald : { true, false }) {
        const curcuma::Molecule mol = dimerAt(Position(L / 2, L / 2, L / 2), true, L);
        const RunResult r = run(mdController("gfnff", dt, max_time, 0.0, L, "", ewald), mol, ewald ? "c_per" : "c_nonper", 2);
        check(r.ok, std::string("(c) gfnff + cell") + (ewald ? "" : " + coulomb_ewald false") + ": set up");
        if (!r.ok)
            continue;
        std::printf("(c) gfnff + cell%s: periodic %d, DOF %d\n", ewald ? "" : " + coulomb_ewald false", r.periodic ? 1 : 0, r.dof);
        const int n = mol.AtomCount();
        if (ewald)
            check(r.periodic && r.dof == 3 * n - 3, "(c) gfnff + cell: periodic, translation only removed, DOF 3N-3");
        else
            check(!r.periodic && r.dof == 3 * n - 6, "(c) gfnff + cell + coulomb_ewald false: not periodic, DOF 3N-6");
    }

    // (d) native gfn2 with a cell is a periodic energy since WP-PERIODIC-XTB P5: wrap rules of the periodic corner
    {
        const curcuma::Molecule mol = dimerAt(Position(L / 2, L / 2, L / 2), true, L);
        const RunResult r = run(mdController("gfn2", dt, max_time, 0.0, L), mol, "d_gfn2", 2);
        check(r.ok, "(d) gfn2 + cell: set up");
        if (r.ok) {
            std::printf("(d) gfn2 + cell: periodic %d, DOF %d\n", r.periodic ? 1 : 0, r.dof);
            check(r.periodic && r.dof == 3 * mol.AtomCount() - 3, "(d) gfn2 + cell: periodic, translation only removed, DOF 3N-3");
        }
    }

    std::printf("%s: %d failure(s)\n", g_failures == 0 ? "ALL PASSED" : "FAILED", g_failures);
    return g_failures == 0 ? 0 : 1;
}
