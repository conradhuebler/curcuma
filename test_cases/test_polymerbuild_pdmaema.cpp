/*
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated 2026 - PDMAEMA built by polymerbuild from src/molecules/pdmaema.xyz.
 *
 * What must hold for a 4-mer, a 10-mer and an alternating (AB)3 with protonated B built with the real optimisation (nothing here is a number from an earlier run):
 *  - composition: n repeat units C8H15NO2 plus two chain-end hydrogens, 26n+2 atoms (4-mer: C32H62N4O8, 106 atoms);
 *  - every atom is bonded (distance criterion): the chain-end caps are hydrogens ON a carbon, not hydrogens
 *    floating a few Angstrom away (the defect this test was written for: an open Xx became an H for the
 *    optimisation without a bond to its carbon and drifted away);
 *  - the bond graph is a tree: atoms-1 bonds (the C=O bonds are single pairs in a distance graph);
 *  - the backbone alternates quaternary carbon and CH2: n carbons carrying the ester, each with a methyl
 *    (n+1 CH3 on quaternary carbons in all: n side methyls and the chain-end CH3), joined by n-1 CH2 groups;
 *  - C-H 1.05 to 1.15 A, C-C 1.45 to 1.60 A, no atom pair closer than 1.0 A (with H) or 1.1 A (heavy);
 *  - bond lengths per type (C-H, C-C, C-N, C-O) in chemical windows; angles per centre (C with 4 or 3 neighbours, N, O)
 *    in windows, mean of the tetrahedral carbons 109.5 +- 3 deg; the backbone C-CH2-C angle is printed, not judged.
 */

#include "src/capabilities/polymerbuild.h"
#include "src/core/molecule.h"

#include "generated/parameter_registry.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
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

struct Atom {
    std::string element;
    double x, y, z;
};

static double distance(const Atom& a, const Atom& b)
{
    return std::sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z));
}

/// One chain of `units` repeat units: every expected count follows from the repeat unit C8H15NO2 plus two chain-end H.
static void checkChain(const std::string& sequence, int units, int protonated, int subchain)
{
    std::printf("polymerbuild PDMAEMA %s%s\n", sequence.c_str(), subchain > 0 ? (" subchain_size " + std::to_string(subchain)).c_str() : "");
    const std::filesystem::path dir = std::filesystem::temp_directory_path() / ("curcuma_test_polymerbuild_pdmaema_" + std::to_string(units) + "_" + std::to_string(protonated) + "_" + std::to_string(subchain));
    std::filesystem::remove_all(dir);
    std::filesystem::create_directories(dir);
    const std::filesystem::path previous = std::filesystem::current_path();
    std::filesystem::current_path(dir);

    json controller;
    controller["fragments"] = { { "A", std::string(CURCUMA_SOURCE_DIR) + "/src/molecules/pdmaema.xyz" },
        { "B", std::string(CURCUMA_SOURCE_DIR) + "/src/molecules/pdmaema_h.xyz" } };
    controller["sequence"] = sequence;
    if (subchain > 0)
        controller["subchain_size"] = subchain; // divide-and-conquer assembly: sub-chains built, cloned, joined
    controller["verbosity"] = 0;
    {
        PolymerBuild builder(controller, true);
        builder.start();
    }

    std::string found;
    for (const auto& entry : std::filesystem::recursive_directory_iterator(dir))
        if (entry.path().filename() == "polymer_A_final.xyz")
            found = entry.path().string();
    check(!found.empty(), "polymer_A_final.xyz written");

    std::vector<Atom> atoms;
    if (!found.empty()) {
        std::ifstream in(found);
        std::string line;
        std::getline(in, line);
        std::getline(in, line);
        while (std::getline(in, line)) {
            std::stringstream ss(line);
            Atom a;
            if (ss >> a.element >> a.x >> a.y >> a.z)
                atoms.push_back(a);
        }
    }
    const int n = static_cast<int>(atoms.size());
    std::map<std::string, int> composition;
    for (const Atom& a : atoms)
        ++composition[a.element];
    // Unit B is unit A plus one H on the nitrogen (C8H16NO2).
    const int expectedAtoms = 26 * units + protonated + 2;
    check(n == expectedAtoms && composition["C"] == 8 * units && composition["H"] == 15 * units + protonated + 2 && composition["N"] == units && composition["O"] == 2 * units,
        "composition of " + std::to_string(units) + " repeat units (" + std::to_string(protonated) + " protonated) plus two chain-end H, " + std::to_string(expectedAtoms) + " atoms (got " + std::to_string(n) + " atoms)");

    const std::map<std::string, double> radius { { "H", 0.31 }, { "C", 0.76 }, { "N", 0.71 }, { "O", 0.66 } };
    std::vector<std::vector<int>> neighbours(n);
    int bonds = 0;
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
            if (radius.count(atoms[i].element) && radius.count(atoms[j].element)
                && distance(atoms[i], atoms[j]) < 1.2 * (radius.at(atoms[i].element) + radius.at(atoms[j].element))) {
                neighbours[i].push_back(j);
                neighbours[j].push_back(i);
                ++bonds;
            }
    std::vector<std::string> unbound;
    for (int i = 0; i < n; ++i)
        if (neighbours[i].empty())
            unbound.push_back(atoms[i].element + std::to_string(i + 1));
    std::string unboundList;
    for (const auto& u : unbound)
        unboundList += u + " ";
    check(unbound.empty(), "no atom without a bond (unbound: " + unboundList + ")");
    check(bonds == n - 1, "the bond graph is a tree: " + std::to_string(bonds) + " bonds for " + std::to_string(n) + " atoms");

    double chMin = 9, chMax = 0, ccMin = 9, ccMax = 0, closest = 9;
    for (int i = 0; i < n; ++i) {
        for (int j : neighbours[i]) {
            if (j < i)
                continue;
            const double d = distance(atoms[i], atoms[j]);
            if (atoms[i].element == "C" && atoms[j].element == "H") { chMin = std::min(chMin, d); chMax = std::max(chMax, d); }
            if (atoms[j].element == "C" && atoms[i].element == "H") { chMin = std::min(chMin, d); chMax = std::max(chMax, d); }
            if (atoms[i].element == "C" && atoms[j].element == "C") { ccMin = std::min(ccMin, d); ccMax = std::max(ccMax, d); }
        }
        for (int j = i + 1; j < n; ++j) {
            const double d = distance(atoms[i], atoms[j]);
            const bool hydrogen = atoms[i].element == "H" || atoms[j].element == "H";
            if (d < (hydrogen ? 1.0 : 1.1))
                closest = std::min(closest, d);
        }
    }
    check(chMin > 1.05 && chMax < 1.15, "C-H bonds " + std::to_string(chMin) + " to " + std::to_string(chMax) + " A");
    check(ccMin > 1.45 && ccMax < 1.60, "C-C bonds " + std::to_string(ccMin) + " to " + std::to_string(ccMax) + " A");
    check(closest == 9, "no pair closer than 1.0 A (with H) or 1.1 A (heavy atoms)");

    // Backbone: carbons that carry the ester (bonded to a C that has two O) and their CH2 partners.
    auto isCarbonyl = [&](int i) {
        int oxygens = 0;
        for (int j : neighbours[i])
            oxygens += atoms[j].element == "O";
        return atoms[i].element == "C" && oxygens == 2;
    };
    auto hydrogens = [&](int i) {
        int h = 0;
        for (int j : neighbours[i])
            h += atoms[j].element == "H";
        return h;
    };
    std::vector<int> quaternary;
    for (int i = 0; i < n; ++i) {
        if (atoms[i].element != "C")
            continue;
        for (int j : neighbours[i])
            if (isCarbonyl(j))
                quaternary.push_back(i);
    }
    int methylOnQuaternary = 0, joiningCH2 = 0;
    for (int q : quaternary) {
        for (int j : neighbours[q])
            if (atoms[j].element == "C" && !isCarbonyl(j) && hydrogens(j) == 3)
                ++methylOnQuaternary;
    }
    for (int i = 0; i < n; ++i) {
        if (atoms[i].element != "C" || hydrogens(i) != 2)
            continue;
        int attached = 0;
        for (int j : neighbours[i])
            attached += std::find(quaternary.begin(), quaternary.end(), j) != quaternary.end();
        joiningCH2 += attached == 2;
    }
    check(static_cast<int>(quaternary.size()) == units, "one ester-bearing carbon per unit (got " + std::to_string(quaternary.size()) + ")");
    // n side methyls, plus the CH3 at the chain end (the last CH2 with its cap H), which also hangs on a quaternary carbon.
    check(methylOnQuaternary == units + 1, "one side methyl per unit and the chain-end CH3 hang on the quaternary carbons (got " + std::to_string(methylOnQuaternary) + ")");
    check(joiningCH2 == units - 1, "CH2 groups join the units (got " + std::to_string(joiningCH2) + ")");

    // Bond-length windows per bond type (generous chemical windows, not fitted to a run): a broken or
    // stretched junction falls outside, a normal GFN-FF geometry inside.
    struct Window { double lo, hi; };
    const std::map<std::string, Window> bondWindow {
        { "C-H", { 1.06, 1.13 } }, { "H-N", { 0.98, 1.08 } }, { "C-C", { 1.48, 1.60 } }, { "C-N", { 1.43, 1.52 } }, { "C-O", { 1.18, 1.45 } }
    };
    std::map<std::string, std::pair<double, double>> bondSeen;
    std::map<std::string, int> bondCount;
    for (int i = 0; i < n; ++i)
        for (int j : neighbours[i]) {
            if (j < i)
                continue;
            std::string a = atoms[i].element, b = atoms[j].element;
            if (a > b)
                std::swap(a, b);
            const std::string type = a + "-" + b;
            const double d = distance(atoms[i], atoms[j]);
            auto it = bondSeen.find(type);
            if (it == bondSeen.end())
                bondSeen[type] = { d, d };
            else {
                it->second.first = std::min(it->second.first, d);
                it->second.second = std::max(it->second.second, d);
            }
            ++bondCount[type];
        }
    for (const auto& [type, range] : bondSeen) {
        auto w = bondWindow.find(type);
        if (w == bondWindow.end())
            check(false, "unexpected bond type " + type + " (" + std::to_string(bondCount[type]) + " bonds)");
        else
            check(range.first >= w->second.lo && range.second <= w->second.hi,
                type + " (" + std::to_string(bondCount[type]) + " bonds) " + std::to_string(range.first) + " to "
                    + std::to_string(range.second) + " A, window " + std::to_string(w->second.lo) + " to "
                    + std::to_string(w->second.hi));
    }

    // Angle distribution by centre: element and number of neighbours (4 = tetrahedral, 3 = planar ester C or amine N).
    // The window rejects distorted junctions; the mean of the tetrahedral centres must sit near 109.5 degrees.
    struct AngleStat { double lo = 999, hi = 0, sum = 0; int count = 0; };
    std::map<std::string, AngleStat> angleStat;
    AngleStat backbone; // C-CH2-C between two quaternary carbons, reported on its own
    auto angleAt = [&](int c, int a, int b) {
        const double ux = atoms[a].x - atoms[c].x, uy = atoms[a].y - atoms[c].y, uz = atoms[a].z - atoms[c].z;
        const double vx = atoms[b].x - atoms[c].x, vy = atoms[b].y - atoms[c].y, vz = atoms[b].z - atoms[c].z;
        const double cosine = (ux * vx + uy * vy + uz * vz) / (std::sqrt(ux * ux + uy * uy + uz * uz) * std::sqrt(vx * vx + vy * vy + vz * vz));
        return std::acos(std::max(-1.0, std::min(1.0, cosine))) * 180.0 / 3.14159265358979323846;
    };
    for (int c = 0; c < n; ++c) {
        const auto& nb = neighbours[c];
        const std::string centre = atoms[c].element + std::to_string(nb.size());
        const bool isBackbone = atoms[c].element == "C" && nb.size() == 4 && hydrogens(c) == 2
            && std::count_if(nb.begin(), nb.end(), [&](int j) { return std::find(quaternary.begin(), quaternary.end(), j) != quaternary.end(); }) == 2;
        for (size_t x = 0; x < nb.size(); ++x)
            for (size_t y = x + 1; y < nb.size(); ++y) {
                const double t = angleAt(c, nb[x], nb[y]);
                AngleStat& st = angleStat[centre];
                st.lo = std::min(st.lo, t);
                st.hi = std::max(st.hi, t);
                st.sum += t;
                ++st.count;
                if (isBackbone && std::find(quaternary.begin(), quaternary.end(), nb[x]) != quaternary.end()
                    && std::find(quaternary.begin(), quaternary.end(), nb[y]) != quaternary.end()) {
                    backbone.lo = std::min(backbone.lo, t);
                    backbone.hi = std::max(backbone.hi, t);
                    backbone.sum += t;
                    ++backbone.count;
                }
            }
    }
    const std::map<std::string, Window> angleWindow {
        { "C4", { 100.0, 128.0 } }, { "C3", { 108.0, 130.0 } }, { "N3", { 105.0, 120.0 } }, { "N4", { 100.0, 120.0 } }, { "O2", { 105.0, 120.0 } }
    };
    for (const auto& [centre, st] : angleStat) {
        auto w = angleWindow.find(centre);
        if (w == angleWindow.end()) {
            check(false, "unexpected centre " + centre);
            continue;
        }
        char text[200];
        std::snprintf(text, sizeof text, "angles at %s: %d angles, %.1f to %.1f deg (mean %.1f), window %.0f to %.0f",
            centre.c_str(), st.count, st.lo, st.hi, st.sum / st.count, w->second.lo, w->second.hi);
        check(st.lo >= w->second.lo && st.hi <= w->second.hi, text);
    }
    if (angleStat.count("C4"))
        check(std::abs(angleStat["C4"].sum / angleStat["C4"].count - 109.5) < 3.0, "mean angle at tetrahedral carbon within 3 deg of 109.5");
    if (backbone.count > 0)
        std::printf("  INFO  backbone C-CH2-C angle (%d): %.1f to %.1f deg, mean %.1f (not tested against a reference)\n",
            backbone.count, backbone.lo, backbone.hi, backbone.sum / backbone.count);

    std::filesystem::current_path(previous);
    std::filesystem::remove_all(dir);
}

int main(int argc, char** argv)
{
    // Chains from the command line: "A30" = (A)30, "AB5" = (AB)5 (B is the protonated unit), ":4" appended = sub-chains of
    // 4 units. Default: (A)4, (A)10, (AB)3, (AB)5 in sub-chains of 4.
    // (AB)3 is a regression case: the placement of the fifth unit (A after ABAB) used to stay at 2.5 A and tear the chain.
    std::vector<std::string> specs { "A4", "A10", "AB3", "AB5:4" };
    if (argc > 1)
        specs.assign(argv + 1, argv + argc);
    for (std::string spec : specs) {
        int subchain = 0;
        const size_t colon = spec.find(':');
        if (colon != std::string::npos) {
            subchain = std::max(2, std::atoi(spec.c_str() + colon + 1));
            spec.erase(colon);
        }
        if (spec.rfind("AB", 0) == 0) {
            const int pairs = std::max(1, std::atoi(spec.c_str() + 2));
            checkChain("(AB)" + std::to_string(pairs), 2 * pairs, pairs, subchain);
        } else {
            const int units = std::max(2, std::atoi(spec.c_str() + 1));
            checkChain("(A)" + std::to_string(units), units, 0, subchain);
        }
    }
    std::printf("%s\n", g_failed ? "FAILED" : "all passed");
    return g_failed ? 1 : 0;
}
