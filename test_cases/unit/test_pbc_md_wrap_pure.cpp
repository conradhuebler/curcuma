/*
 * < Periodic MD: wall_potential pbc is a pure cell wrap when the energy is periodic >
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
 * Claude Generated (Nov 2026, WP-PERIODIC-NONBONDED A6).
 *
 * wall_potential pbc used to refuse a wrap into occupied space (kWrapClearance) and reflect the
 * fragment instead, because the energy had no minimum image. With a unit cell the native CPU
 * GFN-FF energy is periodic (A2-B3), so the wrap is pure bookkeeping. The test checks:
 *   (a) the source of truth EnergyCalculator::hasPeriodicEnergy(): GFN-FF + cell true,
 *       GFN-FF + cell + coulomb_ewald false, GFN-FF without cell and UFF false;
 *   (b) pbc_md_wrap_pure: the registry water dimer in L = 44 Bohr, split across the +x face (monomer A
 *       0.2 A inside, B through the face, MIC O-O 2.8 A) and drifting along +x, so both monomers
 *       cross the face: GFN-FF, rect container = cell, NVE 2 ps,
 *       dt 0.5 fs, COM removal off. The container wrap is the cell wrap (wallWrapByCell), no
 *       refusal, no reflection, cell wraps > 0, energy drift < 1e-4 Eh/ps (D1 level), total
 *       momentum conserved;
 *   (c) UFF and GFN-FF without a cell, same container, a water whose wrap destination is
 *       occupied: the clearance/reflection fallback still runs (refusals and reflections > 0);
 *   (d) container bounds that differ from the cell: with a periodic energy only the cell wraps.
 */

#include "src/capabilities/simplemd.h"
#include "src/core/curcuma_logger.h"
#include "src/core/elements.h"
#include "src/core/energycalculator.h"
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

// Two registry waters (Angstrom), centres of mass at x_a and x_b, both at y = z = L/2.
static curcuma::Molecule waterPair(double L, double x_a, double x_b, bool cell)
{
    const curcuma::Molecule water = TestMolecules::TestMoleculeRegistry::createMolecule("H2O", false);
    const Geometry g = water.getGeometry();
    double mass = 0;
    Position com = Position::Zero();
    for (int k = 0; k < g.rows(); ++k) {
        const double m = Elements::AtomicMass[water.Atom(k).first];
        mass += m;
        com += m * g.row(k).transpose();
    }
    com /= mass;
    curcuma::Molecule mol;
    for (double xc : { x_a, x_b })
        for (int k = 0; k < g.rows(); ++k) {
            const Position d = g.row(k).transpose() - com;
            mol.addPair({ water.Atom(k).first, Position{ d(0) + xc, d(1) + L / 2, d(2) + L / 2 } });
        }
    if (cell)
        mol.setUnitCell(Eigen::Matrix3d::Identity() * L, true);
    return mol;
}

// The registry water dimer (O-O 2.8 A along z) turned so that O-O lies along x and placed across the +x face:
// monomer A (atoms 0-2) with its centre of mass at x_a, monomer B (atoms 3-5) one O-O distance further on, i.e.
// through the face at x_a + 2.8 - L. Raw coordinates put them L - 2.8 A apart; only the minimum image joins them.
static curcuma::Molecule dimerAcrossFace(double L, double x_a)
{
    const curcuma::Molecule dimer = TestMolecules::TestMoleculeRegistry::createMolecule("H2O_dimer", false);
    const Geometry g = dimer.getGeometry();
    Position com[2] = { Position::Zero(), Position::Zero() };
    double mass[2] = { 0, 0 };
    for (int k = 0; k < g.rows(); ++k) {
        const int f = k < 3 ? 0 : 1;
        const double m = Elements::AtomicMass[dimer.Atom(k).first];
        mass[f] += m;
        com[f] += m * g.row(k).transpose();
    }
    com[0] /= mass[0];
    com[1] /= mass[1];
    const double d = (com[1] - com[0]).norm();
    curcuma::Molecule mol;
    for (int k = 0; k < g.rows(); ++k) {
        const int f = k < 3 ? 0 : 1;
        const Position r = g.row(k).transpose() - com[f];
        // rotation z -> x, x -> -z
        const double xc = (f == 0) ? x_a : x_a + d - L;
        mol.addPair({ dimer.Atom(k).first, Position{ r(2) + xc, r(1) + L / 2, -r(0) + L / 2 } });
    }
    mol.setUnitCell(Eigen::Matrix3d::Identity() * L, true);
    return mol;
}

static json restartState(const curcuma::Molecule& mol, const Geometry& v, double dt, double max_time, const std::string& method)
{
    json s;
    s["method"] = method;
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
    return s;
}

static json mdController(const std::string& method, const std::string& restart_file, double dt, double max_time, double L,
    bool ewald = true, double box_shift = 0.0)
{
    json controller;
    json& md = controller["simplemd"];
    md["method"] = method;
    md["restart_file"] = restart_file;
    md["thermostat"] = "none";
    md["dt"] = dt;
    md["max_time"] = max_time;
    md["no_center"] = true;
    md["remove_com_motion"] = 0.0;
    md["rattle"] = 0;
    md["print_frequency"] = 1000000;
    md["dump_frequency"] = 1000000;
    md["write_restart_frequency"] = 1000000;
    md["write_unique"] = false;
    md["verbosity"] = 0;
    md["wall_type"] = "rect";
    md["wall_potential"] = "pbc";
    md["wall_x_min"] = box_shift;
    md["wall_y_min"] = box_shift;
    md["wall_z_min"] = box_shift;
    md["wall_x_max"] = L + box_shift;
    md["wall_y_max"] = L + box_shift;
    md["wall_z_max"] = L + box_shift;
    if (!ewald)
        controller["gfnff"]["coulomb_ewald"] = false;
    return controller;
}

static void writeJson(const std::string& path, const json& j)
{
    std::ofstream f(path);
    f << j.dump(2);
}

// Velocities: water A drifts along +x (internal unit sqrt(Eh/amu), 0.01 = 0.005 A/fs), water B at rest.
static Geometry drift(int natoms, double vx)
{
    Geometry v = Geometry::Zero(natoms, 3);
    for (int i = 0; i < 3; ++i)
        v(i, 0) = vx;
    return v;
}

static Position momentum(const Geometry& v, const curcuma::Molecule& mol)
{
    Position p = Position::Zero();
    for (int i = 0; i < v.rows(); ++i)
        p += Elements::AtomicMass[mol.Atom(i).first] * v.row(i).transpose();
    return p;
}

static bool periodicEnergy(const std::string& method, bool cell, bool ewald)
{
    const double L = 44.0 * CurcumaUnit::Length::BOHR_TO_ANGSTROM;
    curcuma::Molecule mol = waterPair(L, L / 2 - 1.5, L / 2 + 1.5, cell);
    json config;
    config["method"] = method;
    config["verbosity"] = 0;
    if (!ewald)
        config["gfnff"]["coulomb_ewald"] = false;
    EnergyCalculator ec(method, config);
    ec.setMolecule(mol.getMolInfo());
    return ec.hasPeriodicEnergy();
}

int main()
{
    CurcumaLogger::set_verbosity(0);
    const double L = 44.0 * CurcumaUnit::Length::BOHR_TO_ANGSTROM; // 23.28 A
    const double dt = 0.5; // fs

    // (a) the source of truth
    check(periodicEnergy("gfnff", true, true), "(a) gfnff + cell: hasPeriodicEnergy() true");
    check(!periodicEnergy("gfnff", true, false), "(a) gfnff + cell + coulomb_ewald false: hasPeriodicEnergy() false");
    check(!periodicEnergy("gfnff", false, true), "(a) gfnff without cell: hasPeriodicEnergy() false");
    check(!periodicEnergy("uff", true, true), "(a) uff + cell: hasPeriodicEnergy() false");

    // (b) pure wrap, GFN-FF with the cell, NVE 2 ps
    {
        const double max_time = 2000.0;
        curcuma::Molecule mol = dimerAcrossFace(L, L - 0.2);
        Geometry v_in = Geometry::Zero(mol.AtomCount(), 3);
        for (int i = 0; i < mol.AtomCount(); ++i)
            v_in(i, 0) = 0.02; // the whole dimer drifts along +x (about 0.01 A/fs), so A and later B cross the face
        writeJson("pbc_md_wrap_pure_b.json", restartState(mol, v_in, dt, max_time + 10, "gfnff"));
        SimpleMD sim(mdController("gfnff", "pbc_md_wrap_pure_b.json", dt, max_time + 10, L), true);
        sim.overrideBasename("pbc_md_wrap_pure_b");
        sim.setMolecule(mol);
        const bool ok = sim.Initialise();
        check(ok, "(b) SimpleMD set-up");
        if (ok) {
            sim.prepareRun();
            check(sim.wallWrapByCell(), "(b) wall_potential pbc with a periodic energy: wrap by the cell");
            const int nsteps = static_cast<int>(max_time / dt);
            std::vector<double> t, e;
            Position p0 = Position::Zero();
            double max_dp = 0.0, min_contact = 1e9;
            long last_wraps = 0;
            bool finite = true;
            for (int s = 0; s < nsteps; ++s) {
                sim.step();
                const double et = sim.potentialEnergy() + sim.kineticEnergy();
                if (!std::isfinite(et)) {
                    finite = false;
                    break;
                }
                t.push_back((s + 1) * dt * 1e-3);
                e.push_back(et);
                const Position p = momentum(sim.velocities(), mol);
                if (s == 0)
                    p0 = p;
                max_dp = std::max(max_dp, (p - p0).norm());
                if (sim.cellWrapCount() > last_wraps) {
                    // raw (not minimum-image) contact right after the wrap: what the old clearance test saw
                    const Geometry& r = sim.positions();
                    for (int i = 0; i < 3; ++i)
                        for (int j = 3; j < 6; ++j)
                            min_contact = std::min(min_contact, (r.row(i) - r.row(j)).norm());
                    last_wraps = sim.cellWrapCount();
                }
            }
            sim.finalizeRun();
            check(finite, "(b) finite total energy over the run");
            // linear fit of E_tot(t) over the whole run
            double st = 0, se = 0, stt = 0, ste = 0;
            const double n = static_cast<double>(t.size());
            for (size_t k = 0; k < t.size(); ++k) {
                st += t[k];
                se += e[k];
                stt += t[k] * t[k];
                ste += t[k] * e[k];
            }
            const double slope = (n * ste - st * se) / (n * stt - st * st);
            double emin = 1e30, emax = -1e30;
            for (double x : e) {
                emin = std::min(emin, x);
                emax = std::max(emax, x);
            }
            std::printf("(b) steps %zu, cell wraps %ld, container wraps %ld, refused %ld, reflections %ld, "
                        "drift %.3e Eh/ps, max-min %.3e Eh, |P0| %.4e, max |P - P0| %.3e, "
                        "min raw A-B contact after a wrap %.3f A\n",
                t.size(), sim.cellWrapCount(), sim.containerWrapCount(), sim.containerRefusedCount(),
                sim.containerReflectionCount(), slope, emax - emin, p0.norm(), max_dp, min_contact);
            check(sim.cellWrapCount() >= 2, "(b) both monomers wrapped by the cell");
            check(sim.containerRefusedCount() == 0, "(b) no wrap refused");
            check(sim.containerReflectionCount() == 0, "(b) no reflection");
            check(sim.containerWrapCount() == 0, "(b) no container (fallback) wrap");
            check(std::abs(slope) < 1e-4, "(b) NVE drift below 1e-4 Eh/ps");
            check(max_dp < 1e-6 * std::max(p0.norm(), 1e-12) + 1e-9, "(b) total momentum conserved");
        }
    }

    // (c) no periodic energy: the fallback stays (a water lands next to another one)
    for (const std::string method : { std::string("uff"), std::string("gfnff") }) {
        const double max_time = 200.0;
        curcuma::Molecule mol = waterPair(L, L - 0.3, 1.0, false);
        const Geometry v_in = drift(mol.AtomCount(), 0.01);
        const std::string file = "pbc_md_wrap_pure_c_" + method + ".json";
        writeJson(file, restartState(mol, v_in, dt, max_time + 10, method));
        SimpleMD sim(mdController(method, file, dt, max_time + 10, L), true);
        sim.overrideBasename("pbc_md_wrap_pure_c_" + method);
        sim.setMolecule(mol);
        const bool ok = sim.Initialise();
        check(ok, "(c) " + method + ": SimpleMD set-up");
        if (!ok)
            continue;
        sim.prepareRun();
        check(!sim.wallWrapByCell(), "(c) " + method + " without cell: fallback path");
        for (int s = 0; s < static_cast<int>(max_time / dt); ++s)
            sim.step();
        sim.finalizeRun();
        std::printf("(c) %s: container wraps %ld, refused %ld, reflections %ld, cell wraps %ld\n", method.c_str(),
            sim.containerWrapCount(), sim.containerRefusedCount(), sim.containerReflectionCount(), sim.cellWrapCount());
        check(sim.containerRefusedCount() > 0, "(c) " + method + ": a wrap into occupied space is refused");
        check(sim.containerReflectionCount() > 0, "(c) " + method + ": the refused fragment is reflected");
        check(sim.cellWrapCount() == 0, "(c) " + method + ": no cell wrap without a cell");
    }

    // (d) container bounds that differ from the cell ([-L/2, L/2) against [0, L)): the cell wrap and the old container
    // wrap disagree about where a fragment belongs. With a periodic energy only the cell wraps; with the
    // non-periodic fallback (coulomb_ewald false) both run, as before A6, and move the same fragment back and forth.
    for (const bool ewald : { true, false }) {
        const double max_time = 100.0;
        curcuma::Molecule mol = dimerAcrossFace(L, L - 0.2);
        Geometry v_in = Geometry::Zero(mol.AtomCount(), 3);
        for (int i = 0; i < mol.AtomCount(); ++i)
            v_in(i, 0) = 0.02;
        const std::string tag = ewald ? "periodic" : "fallback";
        const std::string file = "pbc_md_wrap_pure_d_" + tag + ".json";
        writeJson(file, restartState(mol, v_in, dt, max_time + 10, "gfnff"));
        SimpleMD sim(mdController("gfnff", file, dt, max_time + 10, L, ewald, -L / 2), true);
        sim.overrideBasename("pbc_md_wrap_pure_d_" + tag);
        sim.setMolecule(mol);
        const bool ok = sim.Initialise();
        check(ok, "(d) " + tag + ": SimpleMD set-up");
        if (!ok)
            continue;
        sim.prepareRun();
        for (int s = 0; s < static_cast<int>(max_time / dt); ++s)
            sim.step();
        sim.finalizeRun();
        std::printf("(d) %s, box [-L/2, L/2): wrap by cell %d, cell wraps %ld, container wraps %ld, refused %ld, reflections %ld\n",
            tag.c_str(), sim.wallWrapByCell() ? 1 : 0, sim.cellWrapCount(), sim.containerWrapCount(), sim.containerRefusedCount(),
            sim.containerReflectionCount());
        if (ewald) {
            check(sim.wallWrapByCell(), "(d) periodic: the cell defines the wrap, the container bounds are ignored");
            check(sim.containerWrapCount() == 0 && sim.containerRefusedCount() == 0, "(d) periodic: no container wrap, no refusal");
            check(sim.cellWrapCount() <= 2, "(d) periodic: no back-and-forth wrapping");
        } else {
            check(!sim.wallWrapByCell(), "(d) fallback: coulomb_ewald false keeps the container path");
            check(sim.containerWrapCount() + sim.containerRefusedCount() > 0, "(d) fallback: the container path acts (pre-A6 behaviour)");
        }
    }

    std::printf("%s: %d failure(s)\n", g_failures == 0 ? "ALL PASSED" : "FAILED", g_failures);
    return g_failures == 0 ? 0 : 1;
}
