/*
 * Golden test for CLI2Json (argv -> controller JSON).
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 * Claude Generated (Nov 2026): fixes today's routing so that later refactoring steps
 * (docs/debt_2026-11/PARAMS.md) show exactly which invocations change.
 *
 * Usage: cli_parser_test <golden.json> [--update]
 * Each case is a command line; the expected controller is stored in the golden file (today's
 * behaviour, including known defects marked in the case list; regenerate with --update).
 */
#include "generated/parameter_registry.h"
#include "src/core/cli_parser.h"
#include "src/core/parameter_registry.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

static const std::vector<std::vector<std::string>> kCases = {
    { "-sp", "w.xyz", "-method", "gfn2" },
    { "-sp", "w.xyz", "-max_iterations", "3", "-method", "gfn2" },
    { "-sp", "w.xyz", "-gfnff.cn_cutoff_bohr", "5.5", "-cn_cutoff_bohr", "6" },
    { "-sp", "w.xyz", "-charge", "-1", "-spin", "1" },
    { "-sp", "w.xyz", "-verbose" },
    { "-sp", "w.xyz", "-threads", "4", "5" },
    { "-opt", "w.xyz", "-write_trajectory", "0" },
    { "-md", "w.xyz", "-temperature", "400", "-max_time", "10" },
    { "-confscan", "c.xyz", "-method", "gfn2" }, // method stays global; H2 (copy into rmsd) happens later in ConfigManager
    { "-confsearch", "w.xyz", "-temperature", "400" },
    { "-casino", "w.xyz", "-temperature", "400" },
    { "-hessian", "w.xyz", "-method", "gfn2", "-electronic_temperature", "500" },
    { "-rmsdtraj", "t.xyz", "-rmsd", "0.5", "-heavy" },
    { "-traj", "t.xyz", "-stride", "5" },
    { "-dock", "-host", "A.xyz", "-guest", "B.xyz", "-step_x", "4" },
    { "-sp", "w.xyz", "-temprature", "400" }, // misspelt key is accepted silently (M-class)
    { "-sp", "w.xyz", "-export_run", "2026run.json" }, // numeric prefix stays a string
    { "-sp", "w.xyz", "-bak", "1abc.log", "-bak", "x.log" },
    { "-sp", "w.xyz", "-threshold", "-1e-6" },
};

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::cerr << "usage: cli_parser_test <golden.json> [--update]\n";
        return 2;
    }
    initialize_generated_registry();
    const std::string golden_path = argv[1];
    const bool update = argc > 2 && std::string(argv[2]) == "--update";

    json actual = json::array();
    for (const auto& args : kCases) {
        std::vector<std::string> store = { "curcuma" };
        store.insert(store.end(), args.begin(), args.end());
        std::vector<char*> v;
        for (auto& s : store)
            v.push_back(&s[0]);
        json entry;
        entry["argv"] = args;
        entry["controller"] = CLI2Json(static_cast<int>(v.size()), v.data());
        actual.push_back(entry);
    }

    if (update) {
        std::ofstream(golden_path) << actual.dump(1) << "\n";
        std::cout << "golden written: " << golden_path << "\n";
        return 0;
    }
    std::ifstream in(golden_path);
    if (!in) {
        std::cerr << "golden file missing: " << golden_path << "\n";
        return 2;
    }
    json expected = json::parse(in);
    int bad = 0;
    for (size_t i = 0; i < actual.size(); ++i) {
        if (i >= expected.size() || expected[i] != actual[i]) {
            ++bad;
            std::cerr << "MISMATCH case " << i << ": " << actual[i]["argv"].dump() << "\n  expected: "
                      << (i < expected.size() ? expected[i]["controller"].dump() : "(none)")
                      << "\n  actual:   " << actual[i]["controller"].dump() << "\n";
        }
    }
    if (expected.size() != actual.size())
        ++bad;
    std::cout << (kCases.size() - bad) << "/" << kCases.size() << " cases match\n";
    return bad == 0 ? 0 : 1;
}
