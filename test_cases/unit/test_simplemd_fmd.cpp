/*
 * < SimpleMD: print_frequency in fs, and the state after a restart (seed, DOF) >
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
 * Claude Generated (Nov 2026, F-MD, docs/WP-PERIODIC-NONBONDED.md section F-MD). Two modes:
 *
 *   print   (ctest simplemd_print_frequency): print_frequency is a time in fs. Registry water dimer, uff,
 *           NVE; the status lines on stdout are counted. Before F-MD print_frequency was read as an int:
 *           0.5 became 0 and the modulus by 0 killed the run with SIGFPE. Expected: dt 1 / freq 0.5 ->
 *           every step; dt 0.5 / freq 1 -> every step (integer frequencies unchanged: the truncated time repeats,
 *           not fixed, see F-MD item 2); dt 1 / freq 2 -> every 2nd step; freq 0 -> every step.
 *
 *   restart (ctest simplemd_restart_state): a fresh run of 5 steps, its restartState() written to a file,
 *           then restarted with -restart_file.
 *           (1) DOF after the restart == DOF of the fresh run for (a) vacuum uff, (b) cell + GFN-FF
 *               (periodic, coulomb_ewald), (c) RATTLE (uff, rattle 1), (d) cell + uff (non-periodic).
 *               Before F-MD a restart kept m_dof = 3N.
 *           (2) CSVR, 20 steps after the restart: -seed 202 and -seed 303 give different kinetic energy
 *               series, -seed 202 twice gives the same; without -seed the file's seed is used (two runs
 *               identical and equal to an explicit -seed of the file's value). Before F-MD the generator
 *               was not seeded after a restart: 202 and 303 were bit-identical.
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
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unistd.h>
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
    std::fflush(stdout);
    if (!ok)
        ++g_failures;
}

static curcuma::Molecule dimer(bool cell)
{
    const double L = 44.0 * CurcumaUnit::Length::BOHR_TO_ANGSTROM;
    const curcuma::Molecule d = TestMolecules::TestMoleculeRegistry::createMolecule("H2O_dimer", false);
    const Geometry g = d.getGeometry();
    curcuma::Molecule mol;
    for (int k = 0; k < g.rows(); ++k)
        mol.addPair({ d.Atom(k).first, Position(g.row(k).transpose() + Position(L / 2, L / 2, L / 2)) });
    if (cell)
        mol.setUnitCell(Eigen::Matrix3d::Identity() * L, true);
    return mol;
}

static json mdController(const std::string& method, const std::string& thermostat, double dt, double print,
    int verbosity, int rattle = 0, const std::string& restart_file = "")
{
    json controller;
    json& md = controller["simplemd"];
    md["method"] = method;
    if (!restart_file.empty())
        md["restart_file"] = restart_file;
    md["no_restart"] = restart_file.empty();
    md["thermostat"] = thermostat;
    md["temperature"] = 300.0;
    md["seed"] = 1;
    md["dt"] = dt;
    md["max_time"] = 1000.0;
    md["no_center"] = true;
    md["remove_com_motion"] = 10.0;
    md["remove_com_mode"] = 1;
    md["rattle"] = rattle;
    md["print_frequency"] = print;
    md["dump_frequency"] = 1000000;
    md["write_restart_frequency"] = 1000000;
    md["write_unique"] = false;
    md["write_xyz"] = false;
    md["verbosity"] = verbosity;
    return controller;
}

// ---------------------------------------------------------------- print mode
static bool isStatusLine(const std::string& line)
{
    std::istringstream in(line);
    std::string tok;
    int n = 0;
    while (in >> tok) {
        char* end = nullptr;
        std::strtod(tok.c_str(), &end);
        if (end == tok.c_str() || *end != '\0')
            return false;
        ++n;
    }
    return n >= 15;
}

// Runs nsteps step() calls with stdout redirected to a file; returns the number of status lines.
static int countStatusLines(double dt, double print, int nsteps, const std::string& tag)
{
    const std::string file = "simplemd_fmd_print_" + tag + ".out";
    std::fflush(stdout);
    std::cout.flush();
    const int saved = dup(fileno(stdout));
    if (!std::freopen(file.c_str(), "w", stdout))
        return -1;
    {
        CurcumaLogger::set_verbosity(1);
        SimpleMD sim(mdController("uff", "none", dt, print, 1), false); // not silent: the status lines are the probe
        sim.overrideBasename("simplemd_fmd_print_" + tag);
        sim.setMolecule(dimer(false));
        if (sim.Initialise()) {
            sim.prepareRun();
            for (int s = 0; s < nsteps; ++s)
                sim.step();
            sim.finalizeRun();
        }
    }
    std::cout.flush();
    std::fflush(stdout);
    dup2(saved, fileno(stdout));
    close(saved);
    std::ifstream in(file);
    std::string line;
    int lines = 0;
    while (std::getline(in, line))
        if (isStatusLine(line))
            ++lines;
    return lines;
}

static int printMode()
{
    struct Case {
        double dt, print;
        int nsteps, expected;
        const char* what;
    };
    // step() k = 0 .. nsteps-1 prints for k >= 1 when due, so nsteps steps give at most nsteps-1 lines; plus the
    // line at t = 0 (prepareRun) and the final line (finalizeRun), which do not depend on print_frequency: +2
    const Case cases[] = {
        { 1.0, 0.5, 7, 6 + 2, "dt 1 fs, print_frequency 0.5 fs: every step, no SIGFPE" },
        { 0.5, 1.0, 21, 20 + 2, "dt 0.5 fs, print_frequency 1 fs: every step, as before F-MD (repeated times kept)" },
        { 1.0, 2.0, 11, 5 + 2, "dt 1 fs, print_frequency 2 fs: every 2nd step (unchanged)" },
        { 1.0, 0.0, 5, 4 + 2, "print_frequency 0: every step, no SIGFPE" },
    };
    int k = 0;
    for (const Case& c : cases) {
        const int n = countStatusLines(c.dt, c.print, c.nsteps, std::to_string(k++));
        std::printf("print: %s -> %d status lines (expected %d)\n", c.what, n, c.expected);
        check(n == c.expected, c.what);
    }
    std::printf("%s: %d failure(s)\n", g_failures == 0 ? "ALL PASSED" : "FAILED", g_failures);
    return g_failures == 0 ? 0 : 1;
}

// ---------------------------------------------------------------- restart mode
struct Run {
    bool ok = false;
    int dof = 0;
    std::vector<double> ekin;
};

static Run runMD(const json& controller, const curcuma::Molecule& mol, const std::string& tag, int nsteps,
    json* state = nullptr)
{
    Run r;
    SimpleMD sim(controller, true);
    sim.overrideBasename("simplemd_fmd_" + tag);
    sim.setMolecule(mol);
    if (!sim.Initialise())
        return r;
    sim.prepareRun();
    r.dof = sim.degreesOfFreedom();
    for (int s = 0; s < nsteps; ++s) {
        sim.step();
        r.ekin.push_back(sim.kineticEnergy());
    }
    if (state)
        *state = sim.restartState();
    sim.finalizeRun();
    r.ok = true;
    return r;
}

static std::string writeState(const json& state, const std::string& tag)
{
    const std::string file = "simplemd_fmd_" + tag + ".restart.json";
    std::ofstream f(file);
    f << state.dump(1);
    return file;
}

static int restartMode()
{
    CurcumaLogger::set_verbosity(0);

    // (1) DOF after a restart, the four corners
    struct DofCase {
        const char* tag;
        std::string method;
        bool cell;
        int rattle;
        const char* what;
    };
    const DofCase dcases[] = {
        { "a_vac", "uff", false, 0, "(a) vacuum uff" },
        { "b_gfnff_cell", "gfnff", true, 0, "(b) cell + GFN-FF (periodic)" },
        { "c_rattle", "uff", false, 1, "(c) RATTLE (uff, rattle 1)" },
        { "d_uff_cell", "uff", true, 0, "(d) cell + uff (non-periodic)" },
    };
    for (const DofCase& c : dcases) {
        const curcuma::Molecule mol = dimer(c.cell);
        json state;
        const Run fresh = runMD(mdController(c.method, "none", 0.5, 1e6, 0, c.rattle), mol, std::string(c.tag) + "_fresh", 5, &state);
        check(fresh.ok, std::string(c.what) + ": fresh run set up");
        if (!fresh.ok)
            continue;
        const std::string file = writeState(state, c.tag);
        const Run rest = runMD(mdController(c.method, "none", 0.5, 1e6, 0, c.rattle, file), mol, std::string(c.tag) + "_restart", 2);
        std::printf("%s: DOF fresh %d, after restart %d (3N = %d)\n", c.what, fresh.dof, rest.dof, 3 * mol.AtomCount());
        check(rest.ok && rest.dof == fresh.dof, std::string(c.what) + ": DOF after restart == DOF of the fresh run");
    }

    // (2) CSVR noise after a restart follows -seed
    const curcuma::Molecule mol = dimer(false);
    json state;
    const Run fresh = runMD(mdController("uff", "csvr", 0.5, 1e6, 0), mol, "seed_fresh", 5, &state);
    check(fresh.ok, "(2) fresh CSVR run set up");
    const std::string file = writeState(state, "seed");
    auto restarted = [&](int seed, const std::string& tag) {
        json c = mdController("uff", "csvr", 0.5, 1e6, 0, 0, file);
        if (seed == INT32_MIN)
            c["simplemd"].erase("seed");
        else
            c["simplemd"]["seed"] = seed;
        return runMD(c, mol, tag, 20);
    };
    const Run s202a = restarted(202, "seed_202a"), s202b = restarted(202, "seed_202b"), s303 = restarted(303, "seed_303");
    const Run n1 = restarted(INT32_MIN, "seed_none1"), n2 = restarted(INT32_MIN, "seed_none2");
    const Run file_seed = restarted(state.value("seed", 1), "seed_file");
    check(s202a.ok && s202b.ok && s303.ok && n1.ok && n2.ok && file_seed.ok, "(2) restarts set up");
    auto maxdiff = [](const Run& a, const Run& b) {
        double d = a.ekin.size() == b.ekin.size() ? 0.0 : 1e30;
        for (size_t k = 0; k < a.ekin.size() && k < b.ekin.size(); ++k)
            d = std::max(d, std::abs(a.ekin[k] - b.ekin[k]));
        return d;
    };
    std::printf("(2) max |dEkin| over 20 CSVR steps: 202 vs 202 %.3e, 202 vs 303 %.3e, no seed vs no seed %.3e, "
                "no seed vs file seed %lld %.3e Eh\n",
        maxdiff(s202a, s202b), maxdiff(s202a, s303), maxdiff(n1, n2), static_cast<long long>(state.value("seed", 1)),
        maxdiff(n1, file_seed));
    check(maxdiff(s202a, s202b) == 0.0, "(2) -seed 202 twice after a restart: identical");
    check(maxdiff(s202a, s303) > 0.0, "(2) -seed 202 vs 303 after a restart: different thermostat noise");
    check(maxdiff(n1, n2) == 0.0, "(2) no -seed after a restart: two runs identical (seed from the file)");
    check(maxdiff(n1, file_seed) == 0.0, "(2) no -seed == explicit -seed of the file's value");

    std::printf("%s: %d failure(s)\n", g_failures == 0 ? "ALL PASSED" : "FAILED", g_failures);
    return g_failures == 0 ? 0 : 1;
}

int main(int argc, char** argv)
{
    const std::string mode = argc > 1 ? argv[1] : "";
    if (mode == "print")
        return printMode();
    if (mode == "restart")
        return restartMode();
    std::fprintf(stderr, "usage: test_simplemd_fmd print|restart\n");
    return 2;
}
