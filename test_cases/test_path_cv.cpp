/*
 * <Validation of the path collective variables (Branduardi s,z).>
 * Copyright (C) 2019 - 2026 Conrad Huebler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated 2026 (AI-generated, machine-tested).
 *
 * Checks that (a) s runs monotonically 0 -> 1 along the reference path, (b) s is an
 * explicit function of the geometry alone (the property the NEB-MD arclength lacks),
 * and (c) the analytic ds/dR matches finite differences. Without (c) the CV cannot be
 * used to bias a simulation.
 */

#include <cmath>
#include <iostream>
#include <vector>

#include "generated/parameter_registry.h"
#include "src/capabilities/path_cv.h"
#include "src/core/molecule.h"

using curcuma::Molecule;

// The registry is normally populated by main.cpp; a standalone test has to do it
// itself or ConfigManager throws on the first missing rmsd parameter.
static struct RegistryInitializer {
    RegistryInitializer() { initialize_generated_registry(); }
} registry_initializer;

static int failures = 0;

static void check(bool ok, const std::string& what)
{
    std::cout << (ok ? "  PASS  " : "  FAIL  ") << what << std::endl;
    if (!ok)
        ++failures;
}

/*! \brief Build a trivial two-atom path: the second atom slides along z. */
static std::vector<Geometry> makePath(int n, Molecule& templ)
{
    templ = Molecule();
    Position p;
    p << 0.0, 0.0, 0.0;
    templ.addPair({ 6, p });
    p << 0.0, 0.0, 1.5;
    templ.addPair({ 6, p });

    std::vector<Geometry> path;
    for (int i = 0; i < n; ++i) {
        Geometry g = Geometry::Zero(2, 3);
        g(1, 2) = 1.5 + 1.5 * static_cast<double>(i) / (n - 1); // 1.5 -> 3.0 A
        path.push_back(g);
    }
    return path;
}

int main()
{
    std::cout << "=== PathCV validation ===" << std::endl;

    Molecule templ;
    const int N = 6;
    std::vector<Geometry> path = makePath(N, templ);
    PathCV cv(path, templ);

    std::cout << "lambda = " << cv.lambda() << ", frames = " << cv.frames() << std::endl;

    // (a) s must increase monotonically along the reference frames.
    bool monotone = true;
    double prev = -1.0;
    for (int i = 0; i < N; ++i) {
        const double s = cv.s(path[i]);
        std::cout << "   frame " << i << ": s = " << s << std::endl;
        if (s < prev - 1e-9)
            monotone = false;
        prev = s;
    }
    check(monotone, "s increases monotonically along the path");

    // (b) endpoints must map near 0 and 1.
    check(std::abs(cv.s(path[0])) < 0.25, "s(first frame) is near 0");
    check(std::abs(cv.s(path[N - 1]) - 1.0) < 0.25, "s(last frame) is near 1");

    // (c) analytic gradient vs finite differences at an off-path point.
    Geometry R = path[2];
    R(1, 2) += 0.13;   // displace so we are not sitting exactly on a frame
    R(0, 0) += 0.05;

    double s0 = 0.0, z0 = 0.0;
    Geometry ds;
    cv.evaluate(R, &s0, &z0, &ds);
    std::cout << "test point: s = " << s0 << ", z = " << z0 << std::endl;

    const double h = 1e-5;
    double max_err = 0.0, max_ref = 0.0;
    for (int a = 0; a < R.rows(); ++a) {
        for (int k = 0; k < 3; ++k) {
            Geometry Rp = R, Rm = R;
            Rp(a, k) += h;
            Rm(a, k) -= h;
            const double num = (cv.s(Rp) - cv.s(Rm)) / (2.0 * h);
            const double err = std::abs(num - ds(a, k));
            max_err = std::max(max_err, err);
            max_ref = std::max(max_ref, std::abs(num));
        }
    }
    std::cout << "max |analytic - numeric| = " << max_err
              << "  (largest numeric component " << max_ref << ")" << std::endl;
    // The RMSD gradient ignores the derivative of the best-fit rotation, so an exact
    // match is not expected; require agreement well below the gradient magnitude.
    check(max_err < 0.05 * std::max(max_ref, 1e-6) + 1e-4,
          "analytic ds/dR agrees with finite differences");

    // (c2) analytic dz/dR vs finite differences - needed for the tube restraint that
    // keeps the sampling on the path.
    {
        double sz = 0.0, zz = 0.0;
        Geometry dsz, dzz;
        cv.evaluate(R, &sz, &zz, &dsz, &dzz);
        double zmax_err = 0.0, zmax_ref = 0.0;
        for (int a = 0; a < R.rows(); ++a) {
            for (int k = 0; k < 3; ++k) {
                Geometry Rp = R, Rm = R;
                Rp(a, k) += h;
                Rm(a, k) -= h;
                const double num = (cv.z(Rp) - cv.z(Rm)) / (2.0 * h);
                zmax_err = std::max(zmax_err, std::abs(num - dzz(a, k)));
                zmax_ref = std::max(zmax_ref, std::abs(num));
            }
        }
        std::cout << "max |analytic - numeric| dz/dR = " << zmax_err
                  << "  (largest numeric component " << zmax_ref << ")" << std::endl;
        check(zmax_err < 0.05 * std::max(zmax_ref, 1e-6) + 1e-4,
              "analytic dz/dR agrees with finite differences");
    }

    // (d) the CV must be a pure function of R - evaluating twice gives the same value.
    check(std::abs(cv.s(R) - s0) < 1e-12, "s is deterministic (pure function of R)");

    std::cout << (failures == 0 ? "ALL TESTS PASSED" : "SOME TESTS FAILED") << std::endl;
    return failures == 0 ? 0 : 1;
}
