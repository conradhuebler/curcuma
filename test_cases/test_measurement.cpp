/*
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated 2026 — the measurement capability.
 *
 * The numbers are checked against geometry worked out by hand, not against a
 * previous run of the same code: a capability that agrees with itself proves
 * nothing. Where a value is a known constant of the test geometry (a right angle,
 * a 90-degree dihedral, the gyration radius of a square) that constant is what it
 * is compared against.
 */

#include "src/capabilities/measurement.h"
#include "src/core/molecule.h"

#include "generated/parameter_registry.h"

#include <cmath>
#include <cstdio>
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

static bool near(double a, double b, double tol = 1e-6)
{
    return std::fabs(a - b) < tol;
}

/// A right angle with legs of 1 and 2 Angstrom, and a fourth atom lifted out of
/// the plane so the dihedral is exactly 90 degrees.
static curcuma::Molecule bentChain()
{
    curcuma::Molecule mol;
    mol.addPair({ 6, Position(1.0, 0.0, 0.0) });   // 0
    mol.addPair({ 6, Position(0.0, 0.0, 0.0) });   // 1, the vertex
    mol.addPair({ 6, Position(0.0, 2.0, 0.0) });   // 2
    mol.addPair({ 6, Position(0.0, 2.0, 3.0) });   // 3, straight up out of the plane
    return mol;
}

/// Four atoms on a unit square in the xy plane: every point is 1/sqrt(2) from the
/// centre, so the unweighted radius of gyration is exactly that.
static curcuma::Molecule square()
{
    curcuma::Molecule mol;
    mol.addPair({ 6, Position(0.5, 0.5, 0.0) });
    mol.addPair({ 6, Position(-0.5, 0.5, 0.0) });
    mol.addPair({ 6, Position(-0.5, -0.5, 0.0) });
    mol.addPair({ 6, Position(0.5, -0.5, 0.0) });
    return mol;
}

static json run(const curcuma::Molecule& mol, const json& controller)
{
    curcuma::Measurement measurement(controller, true);
    measurement.setMolecule(mol);
    measurement.start();
    return measurement.Results();
}

int main()
{
    std::printf("measurement capability\n");

    {
        json controller;
        controller["kind"] = "distance";
        controller["atoms"] = "1,2";     // one-based grammar: atoms 0 and 1
        const json r = run(bentChain(), controller);
        check(r.contains("values") && near(r["values"][0].get<double>(), 1.0),
            "a distance is the plain separation (1.0 A)");
        check(r.value("unit", std::string()) == "Angstrom", "and is reported in Angstrom");
    }
    {
        json controller;
        controller["kind"] = "angle";
        controller["atoms"] = "1,2,3";
        const json r = run(bentChain(), controller);
        check(r.contains("values") && near(r["values"][0].get<double>(), 90.0, 1e-4),
            "a right angle comes out at 90 degrees");
    }
    {
        json controller;
        controller["kind"] = "angle";
        controller["atoms"] = "1,2,3";
        controller["unit"] = "radians";
        const json r = run(bentChain(), controller);
        check(r.contains("values") && near(r["values"][0].get<double>(), M_PI / 2.0, 1e-4),
            "and the same angle in radians when asked");
        check(r.value("unit", std::string()) == "radians", "with the unit said out loud");
    }
    {
        json controller;
        controller["kind"] = "dihedral";
        controller["atoms"] = "1,2,3,4";
        const json r = run(bentChain(), controller);
        check(r.contains("values") && near(std::fabs(r["values"][0].get<double>()), 90.0, 1e-4),
            "a fourth atom lifted straight out of the plane gives a 90 degree dihedral");
    }
    {
        json controller;
        controller["kind"] = "gyration";
        const json r = run(square(), controller);
        check(r.contains("values")
                && near(r["values"][0].get<double>(), 1.0 / std::sqrt(2.0), 1e-6),
            "the unweighted gyration radius of a unit square is 1/sqrt(2)");
    }
    {
        json controller;
        controller["kind"] = "centroid";
        const json r = run(square(), controller);
        check(r.contains("positions") && !r.contains("values"),
            "a centroid is reported as a position, not squeezed into a scalar");
        check(r.contains("positions") && near(r["positions"][0][0].get<double>(), 0.0)
                && near(r["positions"][0][1].get<double>(), 0.0),
            "and the square's centroid is the origin");
    }
    {
        // The kind decides how many atoms it needs, and says so rather than
        // measuring whatever it was handed.
        json controller;
        controller["kind"] = "angle";
        controller["atoms"] = "1,2";
        curcuma::Measurement measurement(controller, true);
        measurement.setMolecule(bentChain());
        measurement.start();
        check(measurement.frameCount() == 0
                && measurement.error().find("needs exactly 3") != std::string::npos,
            "an angle over two atoms is refused, with the count in the message");
    }
    {
        json controller;
        controller["kind"] = "distance";
        controller["atoms"] = "1,2";
        curcuma::Measurement measurement(controller, true);
        measurement.start();
        check(measurement.frameCount() == 0 && !measurement.error().empty(),
            "and so is a measurement with nothing to measure");
    }
    {
        check(curcuma::requiredAtoms(curcuma::MeasurementKind::Dihedral) == 4
                && curcuma::requiredAtoms(curcuma::MeasurementKind::Gyration) == 0,
            "how many atoms a kind needs is a property of the kind");
        curcuma::MeasurementKind kind;
        check(curcuma::parseMeasurementKind("torsion", kind)
                && kind == curcuma::MeasurementKind::Dihedral,
            "the CLI's \"torsion\" is the same kind as \"dihedral\"");
    }

    std::printf("%s (%d failed)\n", g_failed ? "FAIL" : "PASS", g_failed);
    return g_failed ? 1 : 0;
}
