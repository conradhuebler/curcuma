/*
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated 2026 — external potentials: parsing and selection resolution.
 *
 * A configuration error here is the dangerous kind: a bias on a selection that
 * matched nothing does not crash and does not warn, it simply does not act, and
 * the run looks like the potential was too weak. So the parser refuses rather
 * than shrugs, and that is what these check.
 */

#include "src/capabilities/external_potentials.h"
#include "src/core/molecule.h"

// The parameter registry is filled by a generated function, and ConfigManager
// throws from deep inside a driver without it -- the same trap every embedder of
// curcuma meets once.
#include "generated/parameter_registry.h"

#include <cmath>
#include <cstdio>
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

/// Two well-separated water molecules: two fragments, so "F1"/"F2" mean something
/// and a distance restraint between their centroids has room to act.
static curcuma::Molecule twoWaters(double separation)
{
    curcuma::Molecule mol;
    const double geom[3][3] = { { 0.0, 0.0, 0.0 }, { 0.9572, 0.0, 0.0 },
        { -0.2400, 0.9266, 0.0 } };
    const int z[3] = { 8, 1, 1 };
    for (int copy = 0; copy < 2; ++copy) {
        for (int a = 0; a < 3; ++a) {
            mol.addPair({ z[a], Position(geom[a][0] + copy * separation, geom[a][1],
                                   geom[a][2]) });
        }
    }
    return mol;
}

int main()
{
    std::printf("external potentials\n");

    // --- the parser ---------------------------------------------------------
    {
        curcuma::Molecule mol = twoWaters(6.0);
        std::string error;

        json list = json::array();
        json entry;
        entry["kind"] = "distance_harmonic";
        entry["atoms"] = "F1";
        entry["atoms_b"] = "F2";
        entry["k"] = 1.0;
        entry["r0"] = 4.0;
        list.push_back(entry);

        std::vector<curcuma::ExternalPotential> parsed
            = curcuma::parseExternalPotentials(list, mol, &error);
        check(parsed.size() == 1 && error.empty(),
            "a distance_harmonic parses" + (error.empty() ? std::string() : " [" + error + "]"));
        check(parsed.size() == 1 && parsed[0].atoms.size() == 3 && parsed[0].atoms_b.size() == 3,
            "and both selections resolve through the fragment cache");

        json bad = json::array();
        json missing;
        missing["kind"] = "distance_harmonic";
        missing["atoms"] = "F9";
        missing["atoms_b"] = "F2";
        missing["k"] = 1.0;
        missing["r0"] = 4.0;
        bad.push_back(missing);
        curcuma::parseExternalPotentials(bad, mol, &error);
        check(!error.empty() && error.find("matched no atoms") != std::string::npos,
            "a selection that matches nothing is refused, not quietly ignored");

        json noKind = json::array();
        json unknown;
        unknown["kind"] = "spring";
        unknown["atoms"] = "F1";
        noKind.push_back(unknown);
        curcuma::parseExternalPotentials(noKind, mol, &error);
        check(!error.empty() && error.find("kind must be") != std::string::npos,
            "an unknown kind is named in the error");
    }

    // The optimiser half of this test is written and was run: it drives a
    // distance_harmonic on two waters through every backend and checks that the
    // restrained centroid distance lands on r0, which is the measurement that
    // catches a wrong Angstrom/Bohr conversion. It is not here because the
    // wiring it needs is not in the tree -- see TOOL_API_WP.md, "the same
    // potentials in a geometry optimisation", for what was found.

    std::printf("%s (%d failed)\n", g_failed ? "FAIL" : "PASS", g_failed);
    return g_failed ? 1 : 0;
}
