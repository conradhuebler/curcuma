/*
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated 2026 - PDMAEMA built by polymerbuild from src/molecules/pdmaema.xyz.
 *
 * What must hold for a 4-mer built with the real optimisation (nothing here is a number from an earlier run):
 *  - composition: four repeat units C8H15NO2 plus two chain-end hydrogens = C32H62N4O8, 106 atoms;
 *  - every atom is bonded (distance criterion): the chain-end caps are hydrogens ON a carbon, not hydrogens
 *    floating a few Angstrom away (the defect this test was written for: an open Xx became an H for the
 *    optimisation without a bond to its carbon and drifted away);
 *  - the bond graph is a tree: 105 bonds for 106 atoms (the C=O bonds are single pairs in a distance graph);
 *  - the backbone alternates quaternary carbon and CH2: four carbons carrying the ester, each with a methyl
 *    (five CH3 on quaternary carbons in all: four side methyls and the chain-end CH3), joined by three CH2 groups;
 *  - C-H 1.05 to 1.15 A, C-C 1.45 to 1.60 A, no atom pair closer than 1.0 A (with H) or 1.1 A (heavy).
 */

#include "src/capabilities/polymerbuild.h"
#include "src/core/molecule.h"

#include "generated/parameter_registry.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
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

int main()
{
    std::printf("polymerbuild PDMAEMA\n");
    const std::filesystem::path dir = std::filesystem::temp_directory_path() / "curcuma_test_polymerbuild_pdmaema";
    std::filesystem::remove_all(dir);
    std::filesystem::create_directories(dir);
    const std::filesystem::path previous = std::filesystem::current_path();
    std::filesystem::current_path(dir);

    json controller;
    controller["fragments"] = { { "A", std::string(CURCUMA_SOURCE_DIR) + "/src/molecules/pdmaema.xyz" } };
    controller["sequence"] = "(A)4";
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
    check(n == 106 && composition["C"] == 32 && composition["H"] == 62 && composition["N"] == 4 && composition["O"] == 8,
        "composition C32H62N4O8, 106 atoms (got " + std::to_string(n) + " atoms)");

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
    check(quaternary.size() == 4, "four carbons carry the ester (got " + std::to_string(quaternary.size()) + ")");
    // Four side methyls, plus the CH3 at the chain end (the last CH2 with its cap H), which also hangs on a quaternary carbon.
    check(methylOnQuaternary == 5, "four side methyls and the chain-end CH3 hang on the quaternary carbons (got " + std::to_string(methylOnQuaternary) + ")");
    check(joiningCH2 == 3, "three CH2 groups join the four units (got " + std::to_string(joiningCH2) + ")");

    std::filesystem::current_path(previous);
    std::filesystem::remove_all(dir);
    std::printf("%s\n", g_failed ? "FAILED" : "all passed");
    return g_failed ? 1 : 0;
}
