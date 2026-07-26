/*
 * test_geometry_restraints.cpp - Unit tests for the harmonic geometry restraints
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * The restraint is what lets a torsion be DRIVEN to a target while the rest of the molecule relaxes
 * (docs/CONFSEARCH_PROPOSALS.md, P0). Its gradient is added to a force-field gradient, so a wrong
 * sign or a wrong unit conversion would not crash -- it would quietly optimise the wrong surface.
 * Hence these tests check numbers, not behaviour:
 *
 *   1. energy is zero at the target and grows quadratically away from it
 *   2. the deviation wraps the short way round (+170 deg vs -175 deg is 15 deg, not 345)
 *   3. ANALYTIC GRADIENT vs FINITE DIFFERENCES, including the Angstrom -> Bohr conversion
 *   4. an empty restraint set is an exact no-op
 *   5. json round trip (degrees in the file, radians internally)
 *   6. DISTANCE restraint: energy, analytic gradient vs finite differences, action=reaction,
 *      json round trip (used to pull a snapshot with a broken/formed bond back to the reference
 *      topology -- same idea as polymerbuild's interface-bond penalty)
 *
 * Claude Generated (Jul 2026)
 */

#include "src/capabilities/optimisation/geometry_restraints.h"
#include "src/core/curcuma_logger.h"
#include "src/core/molecule.h"

#include <cmath>
#include <iostream>
#include <vector>

static int tests_run = 0;
static int tests_passed = 0;
static int tests_failed = 0;

#define TEST_ASSERT(condition, msg)                       \
    do {                                                  \
        tests_run++;                                      \
        if (condition) {                                  \
            std::cout << "  [PASS] " << msg << std::endl; \
            tests_passed++;                               \
        } else {                                          \
            std::cout << "  [FAIL] " << msg << std::endl; \
            tests_failed++;                               \
        }                                                 \
    } while (0)

using namespace Optimization;

// n-butane-like backbone: 4 atoms is all a dihedral needs. Positions in Angstrom, deliberately
// asymmetric so no gradient component is zero by symmetry.
static Geometry testGeometry()
{
    Geometry g(5, 3);
    g << 0.000, 0.000, 0.000,
        1.520, 0.000, 0.000,
        2.100, 1.420, 0.150,
        3.590, 1.510, -0.320,
        1.900, -0.550, 0.930; // a fifth atom, not part of the dihedral -> must stay gradient-free
    return g;
}

static double dihedralDegrees(const Geometry& g, int i, int j, int k, int l)
{
    Matrix dummy(4, 3);
    return GFNFF_Geometry::calculateDihedralAngle(Eigen::Vector3d(g.row(i)), Eigen::Vector3d(g.row(j)),
               Eigen::Vector3d(g.row(k)), Eigen::Vector3d(g.row(l)), dummy, false)
        * 180.0 / M_PI;
}

int main()
{
    CurcumaLogger::set_verbosity(1);
    std::cout << "=== Dihedral restraint tests ===" << std::endl;

    const Geometry geom = testGeometry();
    const double phi_deg = dihedralDegrees(geom, 0, 1, 2, 3);
    std::cout << "  test dihedral 0-1-2-3 = " << phi_deg << " deg" << std::endl;

    // --- 1. energy at / away from the target -----------------------------------------------------
    {
        GeometryRestraints at_target;
        at_target.add({ 0, 1, 2, 3, phi_deg * M_PI / 180.0, 0.5 });
        Vector grad = Vector::Zero(15);
        const double e0 = at_target.apply(geom, grad);
        TEST_ASSERT(std::abs(e0) < 1e-12, "energy is zero when the torsion sits at its target");
        TEST_ASSERT(grad.norm() < 1e-12, "gradient vanishes at the target");

        // 30 degrees off, k = 0.5 Eh/rad^2 -> E = 0.5 * 0.5 * (pi/6)^2
        GeometryRestraints off;
        off.add({ 0, 1, 2, 3, (phi_deg - 30.0) * M_PI / 180.0, 0.5 });
        Vector g2 = Vector::Zero(15);
        const double e1 = off.apply(geom, g2);
        const double expected = 0.5 * 0.5 * std::pow(30.0 * M_PI / 180.0, 2);
        TEST_ASSERT(std::abs(e1 - expected) < 1e-10,
            "energy is 1/2 k dphi^2 (" + std::to_string(e1) + " vs " + std::to_string(expected) + ")");
        TEST_ASSERT(g2.norm() > 1e-8, "gradient is non-zero away from the target");
    }

    // --- 2. shortest-path wrap-around ------------------------------------------------------------
    {
        // Restrain to a target 350 deg "away" the long way = 10 deg the short way.
        double target = phi_deg + 350.0;
        GeometryRestraints wrapped;
        wrapped.add({ 0, 1, 2, 3, target * M_PI / 180.0, 1.0 });
        Vector grad = Vector::Zero(15);
        const double e = wrapped.apply(geom, grad);
        const double expected = 0.5 * 1.0 * std::pow(10.0 * M_PI / 180.0, 2);
        TEST_ASSERT(std::abs(e - expected) < 1e-9,
            "deviation wraps the short way round (10 deg, not 350 deg)");
    }

    // --- 3. analytic gradient vs finite differences ----------------------------------------------
    {
        GeometryRestraints r;
        r.add({ 0, 1, 2, 3, (phi_deg - 47.0) * M_PI / 180.0, 0.37 }); // arbitrary, nothing special
        Vector analytic = Vector::Zero(15);
        r.apply(geom, analytic);

        // Finite differences in ANGSTROM; the analytic gradient is in Eh/Bohr, so the numerical
        // derivative dE/dAngstrom has to be converted the same way the implementation does.
        const double kAngstromPerBohr = 0.52917721067;
        const double h = 1e-5; // Angstrom
        double max_error = 0.0;
        int worst_component = -1;
        for (int atom = 0; atom < 5; ++atom) {
            for (int c = 0; c < 3; ++c) {
                Geometry plus = geom, minus = geom;
                plus(atom, c) += h;
                minus(atom, c) -= h;
                Vector dummy_p = Vector::Zero(15), dummy_m = Vector::Zero(15);
                const double e_plus = r.apply(plus, dummy_p);
                const double e_minus = r.apply(minus, dummy_m);
                const double numeric = (e_plus - e_minus) / (2.0 * h) * kAngstromPerBohr;
                const double error = std::abs(numeric - analytic[3 * atom + c]);
                if (error > max_error) {
                    max_error = error;
                    worst_component = 3 * atom + c;
                }
            }
        }
        std::cout << "  max |analytic - numerical| = " << max_error
                  << " Eh/Bohr (component " << worst_component << ")" << std::endl;
        TEST_ASSERT(max_error < 1e-7, "analytic gradient matches finite differences (< 1e-7 Eh/Bohr)");

        // The fifth atom is not part of the dihedral and must not feel the restraint at all.
        const double untouched = std::abs(analytic[12]) + std::abs(analytic[13]) + std::abs(analytic[14]);
        TEST_ASSERT(untouched < 1e-14, "atoms outside the dihedral get no restraint gradient");
    }

    // --- 4. empty set is a no-op -----------------------------------------------------------------
    {
        GeometryRestraints none;
        Vector grad = Vector::Zero(15);
        const double e = none.apply(geom, grad);
        TEST_ASSERT(none.empty() && e == 0.0 && grad.norm() == 0.0,
            "an empty restraint set changes neither energy nor gradient");
    }

    // --- 5. json round trip ----------------------------------------------------------------------
    {
        std::vector<DihedralRestraint> in = { { 0, 1, 2, 3, 45.0 * M_PI / 180.0, 0.25 } };
        nlohmann::json config;
        config["dihedral_restraints"] = GeometryRestraints::toJson(in);
        GeometryRestraints out = GeometryRestraints::fromJson(config);
        TEST_ASSERT(out.size() == 1, "json round trip keeps the restraint");
        if (out.size() == 1) {
            const DihedralRestraint& r = out.dihedrals().front();
            TEST_ASSERT(r.i == 0 && r.j == 1 && r.k == 2 && r.l == 3, "indices survive the round trip");
            TEST_ASSERT(std::abs(r.target - 45.0 * M_PI / 180.0) < 1e-12,
                "target survives the degree/radian conversion");
            TEST_ASSERT(std::abs(r.force - 0.25) < 1e-12, "force constant survives the round trip");
        }
        // A malformed entry must be dropped, not crash.
        nlohmann::json bad;
        bad["dihedral_restraints"] = nlohmann::json::array({ nlohmann::json::object() });
        TEST_ASSERT(GeometryRestraints::fromJson(bad).empty(), "entries without indices are dropped");
    }

    // --- 6. distance restraints ------------------------------------------------------------------
    {
        const double target = 1.10, k = 2.0; // the 0-1 pair sits at 1.52 A
        GeometryRestraints r;
        r.add(DistanceRestraint{ 0, 1, target, k });
        Vector analytic = Vector::Zero(15);
        const double e = r.apply(geom, analytic);
        const double dist = (Eigen::Vector3d(geom.row(0)) - Eigen::Vector3d(geom.row(1))).norm();
        const double expected = 0.5 * k * std::pow(dist - target, 2);
        TEST_ASSERT(std::abs(e - expected) < 1e-12,
            "distance restraint energy is 1/2 k (r - r0)^2 (" + std::to_string(e) + ")");

        const double kAngstromPerBohr = 0.52917721067;
        const double h = 1e-6;
        double max_error = 0.0;
        for (int atom = 0; atom < 5; ++atom) {
            for (int c = 0; c < 3; ++c) {
                Geometry plus = geom, minus = geom;
                plus(atom, c) += h;
                minus(atom, c) -= h;
                Vector dp = Vector::Zero(15), dm = Vector::Zero(15);
                const double numeric = (r.apply(plus, dp) - r.apply(minus, dm)) / (2.0 * h) * kAngstromPerBohr;
                max_error = std::max(max_error, std::abs(numeric - analytic[3 * atom + c]));
            }
        }
        std::cout << "  distance restraint: max |analytic - numerical| = " << max_error << " Eh/Bohr" << std::endl;
        TEST_ASSERT(max_error < 1e-7, "distance gradient matches finite differences (< 1e-7 Eh/Bohr)");

        double sum = 0.0;
        for (int c = 0; c < 3; ++c)
            sum += std::abs(analytic[c] + analytic[3 + c]);
        TEST_ASSERT(sum < 1e-14, "the pair's gradient contributions cancel (equal and opposite)");
        double others = 0.0;
        for (int idx = 6; idx < 15; ++idx)
            others += std::abs(analytic[idx]);
        TEST_ASSERT(others < 1e-14, "atoms outside the pair get no restraint gradient");

        std::vector<DistanceRestraint> in = { { 2, 3, 1.45, 3.5 } };
        nlohmann::json config;
        config["distance_restraints"] = GeometryRestraints::toJson(in);
        GeometryRestraints out = GeometryRestraints::fromJson(config);
        TEST_ASSERT(out.size() == 1 && !out.distances().empty(), "distance json round trip keeps the restraint");
        if (!out.distances().empty()) {
            const DistanceRestraint& d = out.distances().front();
            TEST_ASSERT(d.i == 2 && d.j == 3 && std::abs(d.target - 1.45) < 1e-12
                    && std::abs(d.force - 3.5) < 1e-12,
                "distance restraint survives the round trip unchanged");
        }
    }

    std::cout << "\n=== Summary: " << tests_passed << "/" << tests_run << " passed, "
              << tests_failed << " failed ===" << std::endl;
    return tests_failed == 0 ? 0 : 1;
}
