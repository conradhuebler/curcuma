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


    // --- hydrogen bonds ---------------------------------------------------------------
    // Water dimer along x. O1 at the origin donates H1 to O2; the other three hydrogens point away.
    // Linear: D-H...A = 180 deg, O...O = 2.8 A. Bent: H1 at 60 deg from the O-O axis, 100.3 deg at H.
    const auto dimer = [](double ooDistance, bool bent) {
        curcuma::Molecule m;
        m.addPair({ 8, Position(0, 0, 0) });
        m.addPair({ 1, bent ? Position(0.48, 0.831384, 0) : Position(0.96, 0, 0) });
        m.addPair({ 1, Position(-0.24, 0.93, 0) });
        m.addPair({ 8, Position(ooDistance, 0, 0) });
        m.addPair({ 1, Position(ooDistance + 0.24, 0.93, 0) });
        m.addPair({ 1, Position(ooDistance + 0.24, -0.93, 0) });
        return m;
    };
    const std::vector<int> waterA { 0, 1, 2 }, waterB { 3, 4, 5 };
    curcuma::HBondCriteria standard;
    const curcuma::HBondCount linear = curcuma::countHydrogenBonds(dimer(2.8, false), waterA, waterB, standard);
    check(linear.aDonates == 1 && linear.bDonates == 0, "linear dimer: one bond, donated by the first water");
    const curcuma::HBondCount reversed = curcuma::countHydrogenBonds(dimer(2.8, false), waterB, waterA, standard);
    check(reversed.aDonates == 0 && reversed.bDonates == 1, "groups swapped: the donor is now in group B");
    check(curcuma::countHydrogenBonds(dimer(2.8, false), {}, {}, standard).total() == 1, "one group (all atoms): one bond inside it");
    check(curcuma::countHydrogenBonds(dimer(2.8, true), waterA, waterB, standard).total() == 0, "bent to 100 deg: no bond at the 150 deg limit");
    curcuma::HBondCriteria loose = standard;
    loose.angle = 90.0;
    check(curcuma::countHydrogenBonds(dimer(2.8, true), waterA, waterB, loose).total() == 1, "bent to 100 deg: a bond with the limit at 90 deg");
    check(curcuma::countHydrogenBonds(dimer(3.6, false), waterA, waterB, standard).total() == 0, "O...O 3.6 A: none at the 3.5 A limit");
    curcuma::HBondCriteria wide = standard;
    wide.distance = 4.0;
    check(curcuma::countHydrogenBonds(dimer(3.6, false), waterA, waterB, wide).total() == 1, "O...O 3.6 A: a bond with the limit at 4.0 A");
    curcuma::HBondCriteria nitrogen = standard;
    nitrogen.elements = { 7 };
    check(curcuma::countHydrogenBonds(dimer(2.8, false), waterA, waterB, nitrogen).total() == 0, "only N as donor and acceptor: none in a water dimer");
    check(curcuma::countHydrogenBonds(dimer(2.8, false), waterA, { 3 }, standard).total() == 1
            && curcuma::countHydrogenBonds(dimer(2.8, false), waterA, { 4, 5 }, standard).total() == 0,
        "group B restricted to the acceptor oxygen: found; to its hydrogens only: not");


    // --- hydrogen bonds by GFN-FF triple energy ---------------------------------------
    // Water dimer as in test_hbond_energies: O1 donates to O2. The strongest triple of the bonded dimer is
    // about -12 kJ/mol; "clearly attractive" is taken as below -5 kJ/mol, and the dimer stretched to 6 A as above it.
    const auto gfnffDimer = [](double oo) {
        curcuma::Molecule m;
        m.addPair({ 8, Position(0.0, 0.0, 0.0) });
        m.addPair({ 1, Position(0.9572, 0.0, 0.0) });
        m.addPair({ 1, Position(-0.2400, 0.9266, 0.0) });
        m.addPair({ 8, Position(oo, 0.0, 0.0) });
        m.addPair({ 1, Position(oo + 0.24, 0.93, 0.0) });
        m.addPair({ 1, Position(oo + 0.24, -0.93, 0.0) });
        return m;
    };
    curcuma::HBondEnergyOptions byEnergy;
    byEnergy.limitKJ = -5.0;
    const curcuma::HBondEnergyCount bonded = curcuma::countHydrogenBondsByEnergy(gfnffDimer(2.9), waterA, waterB, byEnergy);
    check(bonded.ok, "energy criterion runs: " + bonded.error);
    check(bonded.count.aDonates == 1 && bonded.count.bDonates == 0, "bonded dimer, limit -5 kJ/mol: one bond, donated by the first water");
    check(bonded.energySumKJ < -5.0 && !bonded.energiesKJ.empty(), "its triple energy is below -5 kJ/mol");
    curcuma::HBondEnergyOptions strict = byEnergy;
    strict.limitKJ = bonded.energySumKJ - 1.0;
    check(curcuma::countHydrogenBondsByEnergy(gfnffDimer(2.9), waterA, waterB, strict).count.total() == 0,
        "a limit 1 kJ/mol below the triple energy counts nothing");
    curcuma::HBondEnergyOptions all = byEnergy;
    all.limitKJ = 0.0;
    check(curcuma::countHydrogenBondsByEnergy(gfnffDimer(2.9), waterA, waterB, all).count.total() >= 1, "limit 0 counts every attractive triple");
    const curcuma::HBondEnergyCount swapped = curcuma::countHydrogenBondsByEnergy(gfnffDimer(2.9), waterB, waterA, byEnergy);
    check(swapped.count.aDonates == 0 && swapped.count.bDonates == 1, "groups swapped: the donor is in group B");
    const curcuma::HBondEnergyCount stretched = curcuma::countHydrogenBondsByEnergy(gfnffDimer(6.0), waterA, waterB, byEnergy);
    check(stretched.ok && stretched.count.total() == 0, "stretched to 6 A: no triple below -5 kJ/mol");

    // Selection with a complement: two waters 6 A apart are two fragments.
    curcuma::Molecule far = dimer(6.0, false);
    check(curcuma::resolveSelection(far, "F1").size() == 3, "F1 is the first water");
    const std::vector<int> rest = curcuma::resolveSelection(far, "!F1");
    check(rest.size() == 3 && rest[0] == 3 && rest[2] == 5, "!F1 is everything but the first water");

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

    // hbonds over a three-frame file: linear, bent, linear
    {
        std::ofstream xyz("dimer.xyz");
        const bool pattern[3] = { false, true, false };
        for (int k = 0; k < 3; ++k) {
            const curcuma::Molecule m = dimer(2.8, pattern[k]);
            xyz << "6\nframe " << k << "\n";
            for (int i = 0; i < 6; ++i) {
                const auto a = m.Atom(i);
                xyz << (a.first == 8 ? "O " : "H ") << a.second.x() << " " << a.second.y() << " " << a.second.z() << "\n";
            }
        }
    }
    json hbController = json::object();
    hbController["trajprops"] = { { "quantity", "hbonds" }, { "selection", "1:3" }, { "selection_b", "4:6" } };
    curcuma::TrajectoryProperties hbProps(hbController, true);
    hbProps.setFile("dimer.xyz");
    hbProps.start();
    const json hb = hbProps.Results();
    check(hbProps.error().empty(), "hbonds capability runs: " + hbProps.error());
    if (hbProps.error().empty()) {
        check(hb["frames_used"] == 3 && near(hb["mean"].get<double>(), 2.0 / 3.0) && near(hb["max"].get<double>(), 1.0)
                && near(hb["min"].get<double>(), 0.0),
            "linear, bent, linear: mean 2/3, min 0, max 1");
        check(near(hb["a_donates"]["mean"].get<double>(), 2.0 / 3.0) && near(hb["b_donates"]["mean"].get<double>(), 0.0),
            "all bonds are donated by group A");
    }
    json hbLoose = json::object();
    hbLoose["trajprops"] = { { "quantity", "hbonds" }, { "selection", "1:3" }, { "selection_b", "4:6" }, { "hb_angle", 90.0 } };
    curcuma::TrajectoryProperties looseProps(hbLoose, true);
    looseProps.setFile("dimer.xyz");
    looseProps.start();
    check(looseProps.error().empty() && near(looseProps.Results()["mean"].get<double>(), 1.0), "angle limit 90 deg: the bent frame counts, mean 1");

    // energy criterion over a three-frame file: bonded, stretched, bonded
    {
        std::ofstream xyz("gfnffdimer.xyz");
        const double oo[3] = { 2.9, 6.0, 2.9 };
        for (int k = 0; k < 3; ++k) {
            const curcuma::Molecule m = gfnffDimer(oo[k]);
            xyz << "6\nframe " << k << "\n";
            for (int i = 0; i < 6; ++i) {
                const auto a = m.Atom(i);
                xyz << (a.first == 8 ? "O " : "H ") << a.second.x() << " " << a.second.y() << " " << a.second.z() << "\n";
            }
        }
    }
    json enController = json::object();
    enController["trajprops"] = { { "quantity", "hbonds" }, { "hb_criterion", "energy" }, { "hb_energy_limit", -5.0 },
        { "selection", "1:3" }, { "selection_b", "4:6" }, { "hb_hist_bins", 10 }, { "hb_hist_max", 20.0 } };
    curcuma::TrajectoryProperties enProps(enController, true);
    enProps.setFile("gfnffdimer.xyz");
    enProps.start();
    const json en = enProps.Results();
    check(enProps.error().empty(), "energy-criterion capability runs: " + enProps.error());
    if (enProps.error().empty()) {
        check(en["frames_used"] == 3 && near(en["mean"].get<double>(), 2.0 / 3.0) && near(en["max"].get<double>(), 1.0)
                && near(en["min"].get<double>(), 0.0),
            "bonded, stretched, bonded: mean 2/3, min 0, max 1");
        check(en["energy_sum_kj_mol"]["mean"].get<double>() < 0.0, "the mean triple-energy sum is attractive");
        check(std::filesystem::exists(en["histogram_file"].get<std::string>()), "histogram file written");
        std::ifstream histogram(en["histogram_file"].get<std::string>());
        std::string line;
        std::getline(histogram, line);   // header
        double total = 0.0;
        int rows = 0;
        while (std::getline(histogram, line)) {
            ++rows;
            total += std::stod(line.substr(line.rfind(',') + 1));
        }
        check(rows == 10, "ten histogram rows");
        // Triples per frame: the two bonded frames contribute their nonzero triples inside [-20, 0], the stretched one is
        // below the printed precision of nothing; so the total is at least 2/3 of one triple per frame.
        check(total >= 2.0 / 3.0 - 1e-9, "the histogram holds at least the bonded triple (2 of 3 frames)");
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
