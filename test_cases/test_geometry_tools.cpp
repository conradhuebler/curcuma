/*
 * test_geometry_tools.cpp — GeometryTools angle, dihedral and radius of gyration
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (Sep 2026): reference values worked out by hand, so the sign
 * convention of the dihedral is pinned deliberately rather than by accident.
 */
#include <cmath>
#include <cstdio>
#include <string>

#include "src/capabilities/rmsd/rmsd_functions.h"
#include "src/tools/geometry.h"

static int g_failed = 0;

static void near(double got, double want, double tol, const std::string& what)
{
    const bool ok = std::isfinite(got) && std::abs(got - want) <= tol;
    std::printf("  %s  %s (got %.4f, want %.4f)\n", ok ? "PASS" : "FAIL", what.c_str(), got, want);
    if (!ok)
        ++g_failed;
}

int main()
{
    // --- angle --------------------------------------------------------------
    near(GeometryTools::Angle({ 1, 0, 0 }, { 0, 0, 0 }, { 0, 1, 0 }), 90.0, 1e-6, "angle 90 deg");
    near(GeometryTools::Angle({ 1, 0, 0 }, { 0, 0, 0 }, { -1, 0, 0 }), 180.0, 1e-6,
        "angle 180 deg (straight, clamped acos)");
    near(GeometryTools::Angle({ 1, 0, 0 }, { 0, 0, 0 }, { 1, 0, 0 }), 0.0, 1e-6,
        "angle 0 deg (folded, clamped acos)");
    // Water at the experimental geometry.
    near(GeometryTools::Angle({ 0.757, 0.586, 0.0 }, { 0, 0, 0 }, { -0.757, 0.586, 0.0 }),
        104.48, 0.05, "water H-O-H");
    near(GeometryTools::Angle({ 1, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }), 0.0, 1e-9,
        "degenerate angle (coincident atoms) returns 0, not NaN");

    // --- dihedral -----------------------------------------------------------
    // Backbone a-b-c-d with b-c along x.
    const Position b { 0, 0, 0 }, c { 1.5, 0, 0 };
    near(GeometryTools::Dihedral({ -0.5, 1.0, 0.0 }, b, c, { 2.0, 1.0, 0.0 }), 0.0, 1e-6,
        "dihedral 0 deg (syn)");
    near(GeometryTools::Dihedral({ -0.5, 1.0, 0.0 }, b, c, { 2.0, -1.0, 0.0 }), 180.0, 1e-6,
        "dihedral 180 deg (anti)");
    // Hand-derived: b1 = (0.5,-1,0), b2 = (1.5,0,0), b3 = (0.5,0,1)
    //   n1 = b1 x b2 = (0,0,1.5), n2 = b2 x b3 = (0,-1.5,0), m = n1 x b2_hat = (0,1.5,0)
    //   atan2(m.n2, n1.n2) = atan2(-2.25, 0) = -90 deg
    near(GeometryTools::Dihedral({ -0.5, 1.0, 0.0 }, b, c, { 2.0, 0.0, 1.0 }), -90.0, 1e-6,
        "dihedral -90 deg (d above the bc axis)");
    near(GeometryTools::Dihedral({ -0.5, 1.0, 0.0 }, b, c, { 2.0, 0.0, -1.0 }), 90.0, 1e-6,
        "dihedral +90 deg (mirrored d flips the sign)");

    // --- radius of gyration -------------------------------------------------
    // Four points on a square of half-diagonal r: every point sits at r from the
    // centroid, so Rg is exactly r.
    Geometry square(4, 3);
    square << 1, 0, 0,
        -1, 0, 0,
        0, 1, 0,
        0, -1, 0;
    near(GeometryTools::GyrationRadius(square), 1.0, 1e-9, "Rg of a unit square is 1");

    // Translating the whole set must not change Rg.
    Geometry moved = GeometryTools::TranslateGeometry(square, Position { 5, -3, 7 });
    near(GeometryTools::GyrationRadius(moved), 1.0, 1e-9, "Rg is translation invariant");

    // Scaling by 2 doubles it.
    near(GeometryTools::GyrationRadius(Geometry(square * 2.0)), 2.0, 1e-9, "Rg scales linearly");

    Geometry single(1, 3);
    single << 3, 4, 5;
    near(GeometryTools::GyrationRadius(single), 0.0, 1e-9, "Rg of one point is 0");
    near(GeometryTools::GyrationRadius(Geometry(0, 3)), 0.0, 1e-9, "Rg of an empty set is 0");

    // --- best-fit rotation on a PLANAR set ----------------------------------
    // A planar structure gives a rank-deficient covariance matrix (det = 0) while
    // its Kabsch rotation is perfectly well defined. The guard in BestFitRotation
    // must key on the matrix norm, not the determinant, or every flat molecule
    // silently stays unaligned.
    Geometry rotated(4, 3);   // the square above, turned 90 deg about z
    for (int i = 0; i < 4; ++i) {
        rotated(i, 0) = -square(i, 1);
        rotated(i, 1) = square(i, 0);
        rotated(i, 2) = square(i, 2);
    }
    const Eigen::Matrix3d R = RMSDFunctions::BestFitRotation(square, rotated, 1);
    near(RMSDFunctions::getRMSD(square, RMSDFunctions::applyRotation(rotated, R)), 0.0, 1e-9,
        "planar set: best-fit rotation recovers a rigid 90 deg turn (RMSD 0)");
    near(RMSDFunctions::getRMSD(square, rotated), std::sqrt(2.0), 1e-9,
        "and the unaligned RMSD of that pair really is sqrt(2)");

    std::printf("%s (%d failed)\n", g_failed == 0 ? "PASS" : "FAIL", g_failed);
    return g_failed == 0 ? 0 : 1;
}
