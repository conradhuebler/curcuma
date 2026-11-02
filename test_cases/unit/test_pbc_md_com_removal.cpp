/*
 * < Periodic MD: COM/rotation removal must not take energy out of a periodic system >
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
 * Claude Generated (Nov 2026, WP-PERIODIC-NONBONDED D1).
 *
 * SimpleMD::RemoveRotation() (remove_com_mode 1, the default, applied every remove_com_motion fs)
 * subtracts the centre-of-mass velocity AND the rigid rotation omega = I^-1 L of the whole system.
 * In a cluster in vacuum L is a constant of motion, so after the first removal nothing is taken
 * out any more. With a unit cell L about the centre of mass of the wrapped coordinates is not
 * conserved (no rotational symmetry, wraps move molecules by lattice vectors), and every removal
 * took 1/2 L^T I^-1 L of kinetic energy out of the run: 60 H2O in L = 44 Bohr lost 6e-3 Eh/ps in
 * NVE, 8 H2O 3e-3 Eh/ps, nothing with the removal switched off (docs/WP-PERIODIC-NONBONDED.md D1).
 *
 * The test drives one MD step of 1e-6 fs (the integrator changes nothing at the 1e-8 level) from
 * a restart state whose velocity field is thermal noise + a centre-of-mass drift + a rigid
 * rotation, with the removal due at step 0:
 *   (a) with a cell the velocities after the step equal v_in minus the COM velocity (rotation kept),
 *   (b) without a cell they equal v_in minus the COM velocity minus omega x (r - r_com) (unchanged
 *       behaviour of the vacuum path),
 *   (c) the kinetic energy after the step matches the respective expectation,
 *   (d) a restart file in the {"MD": {...}} layout that curcuma writes is accepted by -restart_file.
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
#include <random>
#include <string>

using json = nlohmann::json;

// SimpleMD and its RMSDDriver members read their defaults from the registry, which only the
// generated table fills (as main.cpp and confscan_test.cpp do).
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

// 8 copies of the registry water (library structure, Angstrom) on a 2 x 2 x 2 grid in a cubic
// cell of 44 Bohr (23.28 A); the molecules are far apart, so the GFN-FF forces are small and the
// one tiny step is inert.
static curcuma::Molecule waterBox8(bool cell, double L)
{
    const curcuma::Molecule water = TestMolecules::TestMoleculeRegistry::createMolecule("H2O", false);
    const Geometry unit = water.getGeometry().rowwise() - water.getGeometry().colwise().mean();
    curcuma::Molecule mol;
    int n = 0;
    for (int a = 0; a < 2; ++a)
        for (int b = 0; b < 2; ++b)
            for (int c = 0; c < 2; ++c) {
                // a different orientation per molecule (cyclic permutation of the axes)
                const int perm = n % 3;
                for (int k = 0; k < unit.rows(); ++k) {
                    const double p[3] = { unit(k, perm), unit(k, (perm + 1) % 3), unit(k, (perm + 2) % 3) };
                    mol.addPair({ water.Atom(k).first, Position{ p[0] + (a + 0.5) * L / 2, p[1] + (b + 0.5) * L / 2, p[2] + (c + 0.5) * L / 2 } });
                }
                ++n;
            }
    if (cell)
        mol.setUnitCell(Eigen::Matrix3d::Identity() * L, true);
    return mol;
}

struct Expectation {
    Geometry v_trans_only; // v_in - v_com
    Geometry v_trans_rot; // v_in - v_com - omega x (r - r_com)
    double ekin_trans_only = 0, ekin_trans_rot = 0;
};

static double kinetic(const Geometry& v, const std::vector<double>& m)
{
    double e = 0;
    for (int i = 0; i < v.rows(); ++i)
        e += 0.5 * m[i] * v.row(i).squaredNorm();
    return e;
}

// Velocity field and the two reference results, built independently of SimpleMD.
static Geometry velocityField(const curcuma::Molecule& mol, Expectation& exp)
{
    const int n = mol.AtomCount();
    const Geometry r = mol.getGeometry();
    std::vector<double> m(n);
    double mass = 0;
    Position r_com = Position::Zero();
    for (int i = 0; i < n; ++i) {
        m[i] = Elements::AtomicMass[mol.Atom(i).first];
        mass += m[i];
        r_com += m[i] * r.row(i).transpose();
    }
    r_com /= mass;

    std::mt19937 rng(1);
    std::normal_distribution<double> gauss(0.0, 1.0);
    const Position v_drift(0.003, -0.002, 0.001);
    const Position omega(0.002, -0.0015, 0.001);
    Geometry v(n, 3);
    for (int i = 0; i < n; ++i) {
        const double sigma = std::sqrt(0.00095 / m[i]); // kT/m at 300 K in Eh/amu
        const Position d = r.row(i).transpose() - r_com;
        v.row(i) = (Position(gauss(rng), gauss(rng), gauss(rng)) * sigma + v_drift + omega.cross(d)).transpose();
    }

    Position p = Position::Zero();
    for (int i = 0; i < n; ++i)
        p += m[i] * v.row(i).transpose();
    const Position v_com = p / mass;
    exp.v_trans_only = v;
    for (int i = 0; i < n; ++i)
        exp.v_trans_only.row(i) -= v_com.transpose();

    Eigen::Matrix3d inertia = Eigen::Matrix3d::Zero();
    Position L = Position::Zero();
    for (int i = 0; i < n; ++i) {
        const Position d = r.row(i).transpose() - r_com;
        L += m[i] * d.cross(Position(exp.v_trans_only.row(i)));
        inertia += m[i] * (d.squaredNorm() * Eigen::Matrix3d::Identity() - d * d.transpose());
    }
    const Position w = inertia.inverse() * L;
    exp.v_trans_rot = exp.v_trans_only;
    for (int i = 0; i < n; ++i) {
        const Position d = r.row(i).transpose() - r_com;
        exp.v_trans_rot.row(i) -= w.cross(d).transpose();
    }
    exp.ekin_trans_only = kinetic(exp.v_trans_only, m);
    exp.ekin_trans_rot = kinetic(exp.v_trans_rot, m);
    return v;
}

static json restartState(const curcuma::Molecule& mol, const Geometry& v, double dt)
{
    json s;
    s["method"] = "gfnff";
    s["thermostat"] = "none";
    s["dT"] = dt;
    s["MaxTime"] = 1.0;
    s["T"] = 300.0;
    s["currentStep"] = 0.0;
    s["seed"] = 1;
    s["rmrottrans"] = 1;
    s["nocenter"] = true;
    s["COM"] = false;
    s["geometry"] = Tools::Geometry2String(mol.getGeometry());
    s["velocities"] = Tools::Geometry2String(v);
    return s;
}

// One MD step from the restart file; returns false when the set-up failed.
static bool oneStep(const curcuma::Molecule& mol, const std::string& restart_file, const std::string& tag,
    double dt, Geometry& v_out, double& ekin_out)
{
    json controller;
    json& md = controller["simplemd"];
    md["method"] = "gfnff";
    md["restart_file"] = restart_file;
    md["thermostat"] = "none";
    md["dt"] = dt;
    md["max_time"] = 1.0;
    md["no_center"] = true;
    md["remove_com_motion"] = 100.0;
    md["remove_com_mode"] = 1;
    md["print_frequency"] = 1000000;
    md["dump_frequency"] = 1000000;
    md["write_restart_frequency"] = 1000000;
    md["write_unique"] = false;
    md["verbosity"] = 0;
    SimpleMD sim(controller, true);
    sim.overrideBasename("pbc_md_com_removal_" + tag);
    sim.setMolecule(mol);
    if (!sim.Initialise())
        return false;
    sim.prepareRun();
    sim.step();
    v_out = sim.velocities();
    ekin_out = sim.kineticEnergy();
    return true;
}

static void writeJson(const std::string& path, const json& j)
{
    std::ofstream f(path);
    f << j.dump(2);
}

int main()
{
    CurcumaLogger::set_verbosity(0);
    const double L = 44.0 * CurcumaUnit::Length::BOHR_TO_ANGSTROM;
    const double dt = 1e-6; // fs
    // Velocity units: A per MD time unit. The restart file stores 6 decimals (std::to_string), so the
    // round trip of a 1e-2 velocity is accurate to ~1e-5; the step itself moves ~1e-8. The removed
    // rotation, if any, is ~2e-2 (measured 2.4e-2 before the fix), far above the tolerance.
    const double tol = 3e-5;
    const double ekin_rtol = 1e-4;

    // (a) cell: translation removed, rotation kept
    {
        curcuma::Molecule mol = waterBox8(true, L);
        Expectation exp;
        const Geometry v_in = velocityField(mol, exp);
        writeJson("pbc_md_com_removal_cell.json", restartState(mol, v_in, dt));
        Geometry v_out;
        double ekin = 0;
        const bool ok = oneStep(mol, "pbc_md_com_removal_cell.json", "cell", dt, v_out, ekin);
        check(ok, "cell: SimpleMD set-up from the restart file");
        if (ok) {
            const double d_trans = (v_out - exp.v_trans_only).cwiseAbs().maxCoeff();
            const double d_rot = (v_out - exp.v_trans_rot).cwiseAbs().maxCoeff();
            std::printf("cell: |v - (v_in - v_com)| max %.2e, |v - (v_in - v_com - w x r)| max %.2e; Ekin %.8f, "
                        "expected (translation only) %.8f, (translation + rotation) %.8f Eh\n",
                d_trans, d_rot, ekin, exp.ekin_trans_only, exp.ekin_trans_rot);
            check(d_trans < tol, "cell: velocities after the step = v_in - v_com (rigid rotation NOT removed)");
            check(std::abs(ekin - exp.ekin_trans_only) < ekin_rtol * exp.ekin_trans_only,
                "cell: kinetic energy after the step = Ekin(v_in - v_com)");
        }
    }

    // (b) no cell: the vacuum path still removes translation and rotation
    {
        curcuma::Molecule mol = waterBox8(false, L);
        Expectation exp;
        const Geometry v_in = velocityField(mol, exp);
        writeJson("pbc_md_com_removal_nocell.json", restartState(mol, v_in, dt));
        Geometry v_out;
        double ekin = 0;
        const bool ok = oneStep(mol, "pbc_md_com_removal_nocell.json", "nocell", dt, v_out, ekin);
        check(ok, "no cell: SimpleMD set-up from the restart file");
        if (ok) {
            const double d_rot = (v_out - exp.v_trans_rot).cwiseAbs().maxCoeff();
            std::printf("no cell: |v - (v_in - v_com - w x r)| max %.2e; Ekin %.8f, expected %.8f Eh\n",
                d_rot, ekin, exp.ekin_trans_rot);
            check(d_rot < tol, "no cell: velocities after the step = v_in - v_com - w x (r - r_com)");
            check(std::abs(ekin - exp.ekin_trans_rot) < ekin_rtol * exp.ekin_trans_rot,
                "no cell: kinetic energy after the step = Ekin(v_in - v_com - w x r)");
        }
    }

    // (d) the wrapped layout {"MD": {...}} of curcuma's own restart files through -restart_file
    {
        curcuma::Molecule mol = waterBox8(true, L);
        Expectation exp;
        const Geometry v_in = velocityField(mol, exp);
        json wrapped;
        wrapped["MD"] = restartState(mol, v_in, dt);
        writeJson("pbc_md_com_removal_wrapped.json", wrapped);
        Geometry v_out;
        double ekin = 0;
        const bool ok = oneStep(mol, "pbc_md_com_removal_wrapped.json", "wrapped", dt, v_out, ekin);
        check(ok, "wrapped restart: SimpleMD set-up");
        if (ok) {
            // A fresh start would draw Maxwell-Boltzmann velocities at 300 K instead of v_in.
            const double d_trans = (v_out - exp.v_trans_only).cwiseAbs().maxCoeff();
            std::printf("wrapped restart: |v - (v_in - v_com)| max %.2e\n", d_trans);
            check(d_trans < tol, "wrapped restart {\"MD\": {...}}: the velocities of the file were used");
        }
    }

    std::printf("%s: %d failure(s)\n", g_failures == 0 ? "ALL PASSED" : "FAILED", g_failures);
    return g_failures == 0 ? 0 : 1;
}
