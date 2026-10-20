/*
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated 2026 - the streaming trajectory properties.
 *
 * Expected values are worked out from the geometry, not taken from a run: two
 * equal masses d apart have R_g = d/2; two unequal masses m1, m2 have
 * R_g^2 = m1 m2 d^2 / (m1 + m2)^2; the corners of a square of side a have
 * R_g = a / sqrt(2). Block statistics are checked on 1..10 where the block means
 * (1.5, 3.5, ..., 9.5) are known by hand.
 */

#include "src/capabilities/trajprops.h"
#include "src/core/elements.h"
#include "src/core/molecule.h"

#include "generated/parameter_registry.h"

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>

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

static bool near(double a, double b, double tol = 1e-9)
{
    return std::fabs(a - b) < tol;
}

int main()
{
    std::printf("trajprops\n");

    // --- gyration radius -----------------------------------------------------------
    curcuma::Molecule pair;
    pair.addPair({ 6, Position(0, 0, 0) });
    pair.addPair({ 6, Position(2, 0, 0) });
    check(near(curcuma::gyrationRadius(pair, {}, true), 1.0), "two equal masses 2 A apart: R_g = 1");
    check(near(curcuma::gyrationRadius(pair, {}, false), 1.0), "geometric, the same");

    curcuma::Molecule ch;
    ch.addPair({ 6, Position(0, 0, 0) });
    ch.addPair({ 1, Position(1, 0, 0) });
    const double mC = Elements::AtomicMass[6], mH = Elements::AtomicMass[1];
    check(near(curcuma::gyrationRadius(ch, {}, true), std::sqrt(mC * mH) / (mC + mH)),
        "C and H 1 A apart: R_g^2 = m1 m2 d^2 / (m1 + m2)^2");
    check(near(curcuma::gyrationRadius(ch, {}, false), 0.5), "geometric C-H: 0.5");

    curcuma::Molecule square;
    square.addPair({ 6, Position(0, 0, 0) });
    square.addPair({ 6, Position(2, 0, 0) });
    square.addPair({ 6, Position(2, 2, 0) });
    square.addPair({ 6, Position(0, 2, 0) });
    check(near(curcuma::gyrationRadius(square, {}, true), std::sqrt(2.0)), "square of side 2: R_g = sqrt(2)");
    check(near(curcuma::gyrationRadius(square, { 0, 1 }, true), 1.0), "subset of two corners: R_g = 1");

    // --- frame range -----------------------------------------------------------------
    check(curcuma::resolveFrames(10, 0, -1, 0.0, 1).size() == 10, "all frames");
    const auto lastThird = curcuma::resolveFrames(30, 0, -1, 1.0 / 3.0, 1);
    check(lastThird.size() == 10 && lastThird.front() == 20 && lastThird.back() == 29, "last third of 30 frames: 20..29");
    const auto strided = curcuma::resolveFrames(10, 2, 8, 0.0, 3);
    check(strided.size() == 3 && strided[0] == 2 && strided[1] == 5 && strided[2] == 8, "2..8 every third: 2, 5, 8");
    check(curcuma::resolveFrames(10, 12, -1, 0.0, 1).empty(), "first frame beyond the end: nothing");
    check(curcuma::resolveFrames(0, 0, -1, 0.0, 1).empty(), "empty trajectory: nothing");

    // --- statistics --------------------------------------------------------------------
    std::vector<double> ramp;
    for (int i = 1; i <= 10; ++i)
        ramp.push_back(i);
    const curcuma::SeriesStatistics s = curcuma::seriesStatistics(ramp, 5);
    check(s.n == 10 && near(s.mean, 5.5) && near(s.median, 5.5) && near(s.min, 1) && near(s.max, 10), "mean, median, min, max of 1..10");
    check(near(s.stddev, std::sqrt(82.5 / 9.0)), "sample standard deviation of 1..10");
    // block means 1.5, 3.5, 5.5, 7.5, 9.5: sample variance 10, error = sqrt(10 / 4 / 5)... = sqrt(10)/sqrt(5)
    check(s.blocks == 5 && near(s.blockError, std::sqrt(2.0)), "block error of 1..10 in 5 blocks: sqrt(2)");
    const curcuma::SeriesStatistics one = curcuma::seriesStatistics({ 3.0 }, 5);
    check(one.n == 1 && one.stddev == 0.0 && one.blockError == 0.0, "one value: no spread, no block error");

    // --- the capability on a file -------------------------------------------------------
    const std::filesystem::path dir = std::filesystem::temp_directory_path() / "curcuma_test_trajprops";
    std::filesystem::create_directories(dir);
    const std::filesystem::path previous = std::filesystem::current_path();
    std::filesystem::current_path(dir);
    {
        std::ofstream xyz("traj.xyz");
        for (int k = 0; k < 10; ++k) {
            const double d = 1.0 + 0.1 * k;
            xyz << "3\nframe " << k << "\n";
            xyz << "C 0.0 0.0 0.0\nC " << d << " 0.0 0.0\nH 50.0 50.0 50.0\n";
        }
    }
    json controller = json::object();
    controller["trajprops"] = { { "selection", "1:2" }, { "mass_weighted", false }, { "last_fraction", 0.5 }, { "blocks", 1 } };
    curcuma::TrajectoryProperties props(controller, true);
    props.setFile("traj.xyz");
    props.start();
    const json r = props.Results();
    check(props.error().empty(), "capability runs: " + props.error());
    if (props.error().empty()) {
        check(r["frames_in_file"] == 10 && r["frames_used"] == 5 && r["first_frame_used"] == 5 && r["last_frame_used"] == 9,
            "last half of 10 frames: 5..9");
        // d_k = 1 + 0.1 k, R_g = d_k / 2 = 0.5 + 0.05 k, k = 5..9 -> mean 0.5 + 0.05 * 7
        check(near(r["mean"].get<double>(), 0.85, 1e-9), "mean R_g of frames 5..9 is 0.85 (the far hydrogen is excluded by the selection)");
        check(r["atoms"] == 2, "selection \"1:2\" gives two atoms");
        check(std::filesystem::exists(r["series_file"].get<std::string>()), "per-frame CSV written");
    }
    json bad = json::object();
    bad["trajprops"] = { { "selection", "9999" } };
    curcuma::TrajectoryProperties badProps(bad, true);
    badProps.setFile("traj.xyz");
    badProps.start();
    check(!badProps.error().empty(), "a selection that matches nothing is an error");
    std::filesystem::current_path(previous);
    std::filesystem::remove_all(dir);

    std::printf("%s\n", g_failed ? "FAILED" : "all passed");
    return g_failed ? 1 : 0;
}
