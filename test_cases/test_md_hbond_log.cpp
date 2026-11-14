/*
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated 2026 - the opt-in hydrogen-bond log of SimpleMD (hb_log_frequency).
 *
 * The properties checked hold by construction:
 *  - with the option off nothing is written and the run is the reference run;
 *  - with the option on the trajectory is identical to the reference (the log only reads);
 *  - one row per logged step, total = a_donates + b_donates;
 *  - the same trajectory counted with a looser limit never has fewer bonds than with a stricter one;
 *  - with a method that has no hydrogen-bond triples the log is refused and the run goes on.
 * Run on stacked water dimers (18 atoms) for 60 fs; fixed seed, one thread.
 */

#include "src/capabilities/simplemd.h"
#include "src/core/molecule.h"

#include "generated/parameter_registry.h"

#include <cstdint>
#include <cstdlib>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

struct RegistryInitializer {
    RegistryInitializer() { initialize_generated_registry(); }
} registry_initializer;

static int g_failed = 0;

static void check(bool ok, const std::string& what)
{
    std::printf("  %s  %s\n", ok ? "PASS" : "FAIL", what.c_str());
    if (!ok)
        ++g_failed;
}

static curcuma::Molecule waterDimers(int copies)
{
    curcuma::Molecule m;
    for (int c = 0; c < copies; ++c) {
        const double z = 3.2 * c;
        m.addPair({ 8, Position(0.0, 0.0, z) });
        m.addPair({ 1, Position(0.9572, 0.0, z) });
        m.addPair({ 1, Position(-0.2400, 0.9266, z) });
        m.addPair({ 8, Position(2.9, 0.0, z) });
        m.addPair({ 1, Position(3.14, 0.93, z) });
        m.addPair({ 1, Position(3.14, -0.93, z) });
    }
    return m;
}

/// Where a run of basename @p name leaves @p file (SimpleMD writes into <name>.snapshots without a BMT directory).
static std::string where(const std::string& name, const std::string& file)
{
    return name + ".snapshots/" + name + file;
}

/// Content with the run's own name masked, so two runs under different names can be compared.
static std::string slurpMasked(const std::string& name, const std::string& file)
{
    std::ifstream in(where(name, file));
    std::stringstream ss;
    ss << in.rdbuf();
    std::string text = ss.str();
    for (std::size_t pos = text.find(name); pos != std::string::npos; pos = text.find(name, pos))
        text.replace(pos, name.size(), "RUN");
    return text;
}

static std::string slurp(const std::string& path)
{
    std::ifstream in(path);
    std::stringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

/// Run 60 fs of MD under the name @p name; @p extra adds simplemd parameters. Returns false if setup failed.
static bool runHere(const std::string& name, const json& extra)
{
    json md = { { "method", "gfnff" }, { "max_time", 60.0 }, { "time_step", 1.0 }, { "temperature", 300.0 },
        { "thermostat", "none" }, { "seed", 1 }, { "threads", 1 }, { "dump_frequency", 5 }, { "print_frequency", 1000000 },
        { "no_restart", true }, { "verbosity", 0 } };
    for (auto it = extra.begin(); it != extra.end(); ++it)
        md[it.key()] = it.value();
    json controller = { { "simplemd", md }, { "verbosity", 0 } };
    SimpleMD simulation(controller, true);
    simulation.overrideBasename(name);
    simulation.setMolecule(waterDimers(3));
    if (!simulation.Initialise())
        return false;
    simulation.start();
    return true;
}

static std::string g_self;

/// Every run is a separate process of this program: two runs inside one process do not start from the same
/// random state (the second run's initial velocities differ), which would make any comparison meaningless.
static bool run(const std::string& name, const json& extra)
{
    {
        std::ofstream args(name + ".args.json");
        args << extra.dump();
    }
    const std::string command = "\"" + g_self + "\" --run " + name + " > " + name + ".out 2>&1";
    return std::system(command.c_str()) == 0;
}

static std::vector<std::vector<double>> readRows(const std::string& path)
{
    std::vector<std::vector<double>> rows;
    std::ifstream in(path);
    std::string line;
    std::getline(in, line);   // header
    while (std::getline(in, line)) {
        std::vector<double> row;
        std::stringstream ss(line);
        std::string cell;
        while (std::getline(ss, cell, ','))
            row.push_back(std::stod(cell));
        rows.push_back(row);
    }
    return rows;
}

int main(int argc, char** argv)
{
    if (argc >= 3 && std::string(argv[1]) == "--run") {
        std::ifstream args(std::string(argv[2]) + ".args.json");
        std::stringstream text;
        text << args.rdbuf();
        const json extra = json::parse(text.str(), nullptr, false);
        return runHere(argv[2], extra.is_object() ? extra : json::object()) ? 0 : 1;
    }
    g_self = std::filesystem::absolute(argv[0]).string();
    std::printf("md hydrogen-bond log\n");
    const std::filesystem::path dir = std::filesystem::temp_directory_path() / "curcuma_test_md_hbond_log";
    std::filesystem::remove_all(dir);
    std::filesystem::create_directories(dir);
    const std::filesystem::path previous = std::filesystem::current_path();
    std::filesystem::current_path(dir);

    check(run("ref", json::object()), "reference run (option off) works");
    check(!std::filesystem::exists(where("ref", ".hbonds.csv")), "option off: no log file");
    const std::string reference = slurpMasked("ref", ".trj.xyz");
    check(reference.size() > 1000, "reference trajectory written (" + std::to_string(reference.size()) + " bytes)");
    // The settings must have reached SimpleMD: 60 fs written every 5 steps is 13 frames (step 0 included) of 18 atoms.
    std::size_t frames = 0;
    for (std::size_t pos = reference.find("\n18\n"); pos != std::string::npos; pos = reference.find("\n18\n", pos + 1))
        ++frames;
    check(frames == 12 || frames == 13, "the run used max_time 60 and dump_frequency 5 (" + std::to_string(frames) + " frame headers after the first)");

    // Control: two runs without the option must agree, otherwise a difference below proves nothing.
    check(run("ref2", json::object()), "second reference run");
    const std::string reference2 = slurpMasked("ref2", ".trj.xyz");
    if (reference2 != reference) {
        std::stringstream a(reference), b(reference2);
        std::string la, lb;
        int line = 0;
        while (std::getline(a, la) && std::getline(b, lb)) {
            ++line;
            if (la != lb) {
                std::printf("    first difference at line %d:\n      %s\n      %s\n", line, la.c_str(), lb.c_str());
                break;
            }
        }
    }
    check(reference2 == reference, "two reference runs give the same trajectory (the control)");

    check(run("strict", { { "hb_log_frequency", 10 }, { "hb_log_energy_limit", -5.0 }, { "hb_log_selection", "1:6" },
              { "hb_log_selection_b", "7:18" } }), "run with the log on (limit -5, group 1:6 against 7:18)");
    check(slurpMasked("strict", ".trj.xyz") == reference, "log on: the trajectory is identical to the reference run");
    check(std::filesystem::exists(where("strict", ".hbonds.csv")), "log file written");
    const auto strict = readRows(where("strict", ".hbonds.csv"));
    check(strict.size() == 6, "60 steps, every 10th: six rows (got " + std::to_string(strict.size()) + ")");
    bool consistent = !strict.empty();
    for (const auto& r : strict)
        consistent = consistent && r.size() == 6 && r[2] == r[3] + r[4] && r[2] >= 0;
    check(consistent, "six columns, total = a_donates + b_donates, not negative");

    check(run("loose", { { "hb_log_frequency", 10 }, { "hb_log_energy_limit", 0.0 }, { "hb_log_selection", "1:6" },
              { "hb_log_selection_b", "7:18" } }), "run with limit 0");
    const auto loose = readRows(where("loose", ".hbonds.csv"));
    bool monotone = loose.size() == strict.size() && !loose.empty();
    bool sameSteps = monotone;
    for (std::size_t i = 0; monotone && i < loose.size(); ++i) {
        monotone = loose[i][2] >= strict[i][2];
        sameSteps = sameSteps && loose[i][0] == strict[i][0];
    }
    check(sameSteps && monotone, "same steps; limit 0 counts at least as many bonds as limit -5 in every row");
    check(slurpMasked("loose", ".trj.xyz") == reference, "limit 0: trajectory identical to the reference as well");

    check(run("notriples", { { "method", "uff" }, { "hb_log_frequency", 10 } }), "method uff with the log on still runs");
    check(!std::filesystem::exists(where("notriples", ".hbonds.csv")), "method without triples: the log is refused, no file");

    if (g_failed) {
        std::printf("  files in the run directory:\n");
        for (const auto& entry : std::filesystem::recursive_directory_iterator(dir))
            std::printf("    %s (%ju bytes)\n", std::filesystem::relative(entry.path(), dir).c_str(),
                static_cast<std::uintmax_t>(entry.is_regular_file() ? entry.file_size() : 0));
    }
    std::filesystem::current_path(previous);
    std::filesystem::remove_all(dir);
    std::printf("%s\n", g_failed ? "FAILED" : "all passed");
    return g_failed ? 1 : 0;
}
