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
#include <random>
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

    // (e) gradient in a ROTATED frame (Claude Generated 2026). The 2-atom collinear case above
    // cannot see a missing back-rotation of the best-fit gradient: the target is aligned
    // onto the reference before the RMSD is taken, so dRMSD/dR is only valid in the
    // laboratory frame if it is rotated back. Use a non-collinear 5-atom molecule, a
    // bent path and a target rotated by 40 degrees; compare analytic and FD gradients.
    {
        Molecule t5;
        Position q;
        const double X[5][3] = { { 0, 0, 0 }, { 1.4, 0.1, 0.0 }, { -0.5, 1.3, 0.2 }, { -0.4, -0.6, 1.2 }, { 0.3, 0.4, -1.3 } };
        for (const auto& x : X) {
            q << x[0], x[1], x[2];
            t5.addPair({ 6, q });
        }
        std::vector<Geometry> p5;
        for (int i = 0; i < 6; ++i) {
            Geometry g(5, 3);
            for (int a = 0; a < 5; ++a)
                for (int c = 0; c < 3; ++c)
                    g(a, c) = X[a][c];
            const double f = i / 5.0;
            g(1, 0) += 0.8 * f;
            g(2, 1) -= 0.6 * f;
            g(3, 2) += 0.5 * f * f;
            g(4, 0) += 0.4 * std::sin(3.0 * f);
            p5.push_back(g);
        }
        PathCV cv5(p5, t5);
        const double th = 40.0 * M_PI / 180.0;
        Eigen::Matrix3d Rz;
        Rz << std::cos(th), -std::sin(th), 0, std::sin(th), std::cos(th), 0, 0, 0, 1;
        Geometry R5 = 0.5 * (p5[2] + p5[3]);
        R5(1, 1) += 0.07;
        R5(4, 2) -= 0.05;
        R5 = (R5 * Rz.transpose()).eval();
        R5.rowwise() += Eigen::RowVector3d(0.3, -0.2, 0.1);
        double sv = 0, zv = 0;
        Geometry ds, dz;
        cv5.evaluate(R5, &sv, &zv, &ds, &dz);
        double err_s = 0, err_z = 0, nrm_s = 0, nrm_z = 0;
        const double h = 1e-5;
        for (int a = 0; a < 5; ++a)
            for (int c = 0; c < 3; ++c) {
                Geometry Rp = R5, Rm = R5;
                Rp(a, c) += h;
                Rm(a, c) -= h;
                double sp = 0, zp = 0, sm = 0, zm = 0;
                cv5.evaluate(Rp, &sp, &zp, nullptr, nullptr);
                cv5.evaluate(Rm, &sm, &zm, nullptr, nullptr);
                const double fs = (sp - sm) / (2 * h), fz = (zp - zm) / (2 * h);
                err_s += (fs - ds(a, c)) * (fs - ds(a, c));
                err_z += (fz - dz(a, c)) * (fz - dz(a, c));
                nrm_s += fs * fs;
                nrm_z += fz * fz;
            }
        const double rel_s = std::sqrt(err_s / std::max(nrm_s, 1e-30));
        const double rel_z = std::sqrt(err_z / std::max(nrm_z, 1e-30));
        std::cout << "   rotated-frame FD: rel err ds = " << rel_s << ", dz = " << rel_z << std::endl;
        check(rel_s < 1e-5, "ds/dR correct for a rotated (non-collinear) geometry");
        check(rel_z < 1e-5, "dz/dR correct for a rotated (non-collinear) geometry");

        // (z wall, Claude Generated 2026) flat-bottom wall W = k/2 (z - z_wall)^2 for z > z_wall:
        // energy and lab-frame gradient vs a central finite difference (Angstrom displacements) on
        // a randomly displaced structure, once with the wall active (z_wall below z) and once
        // inactive (z_wall above z: W and dW/dR must be exactly zero).
        std::mt19937 gen(4711);
        std::normal_distribution<double> nd(0.0, 0.1);
        Geometry Rr = R5;
        for (int a = 0; a < 5; ++a)
            for (int c = 0; c < 3; ++c)
                Rr(a, c) += nd(gen);
        double zr = 0;
        cv5.evaluate(Rr, nullptr, &zr, nullptr, nullptr);
        std::cout << "   z wall test: z of the random structure = " << zr << " A^2" << std::endl;
        const double kw = 0.8;
        for (double zw : { zr - 0.05, zr + 0.05 }) { // z can be < 0 (log-sum-exp), so shift, not scale
            Geometry gw;
            const double W = cv5.zWallEnergy(Rr, zw, kw, &gw);
            const double Wref = (zr > zw) ? 0.5 * kw * (zr - zw) * (zr - zw) : 0.0;
            double err = 0, nrm = 0;
            for (int a = 0; a < 5; ++a)
                for (int c = 0; c < 3; ++c) {
                    Geometry Rp = Rr, Rm = Rr;
                    Rp(a, c) += h;
                    Rm(a, c) -= h;
                    const double fd = (cv5.zWallEnergy(Rp, zw, kw, nullptr) - cv5.zWallEnergy(Rm, zw, kw, nullptr)) / (2 * h);
                    err += (fd - gw(a, c)) * (fd - gw(a, c));
                    nrm += fd * fd;
                }
            if (zr > zw) {
                const double rel = std::sqrt(err / std::max(nrm, 1e-30));
                std::cout << "   z_wall = " << zw << ": W = " << W << " Eh, |dW/dR| = " << gw.norm()
                          << ", rel err vs FD = " << rel << std::endl;
                check(std::abs(W - Wref) < 1e-14 && W > 0.0 && rel < 1e-6,
                      "z wall active: W = k/2 (z-z_wall)^2 and dW/dR vs FD (rel < 1e-6)");
            } else {
                std::cout << "   z_wall = " << zw << ": W = " << W << ", |dW/dR| = " << gw.norm()
                          << ", |FD| = " << std::sqrt(nrm) << std::endl;
                check(W == 0.0 && gw.norm() == 0.0 && nrm == 0.0, "z wall inactive inside the tube: W = 0, dW/dR = 0");
            }
        }
        // Default stiffness: W = 3 kT at z = 1.2 z_wall.
        const double kT = 9.500e-4, zw0 = 0.5;
        const double k0 = PathCV::zWallDefaultK(kT, zw0);
        check(std::abs(0.5 * k0 * (0.2 * zw0) * (0.2 * zw0) - 3.0 * kT) < 1e-15,
              "zWallDefaultK: W(1.2 z_wall) = 3 kT");
    }

    std::cout << (failures == 0 ? "ALL TESTS PASSED" : "SOME TESTS FAILED") << std::endl;
    return failures == 0 ? 0 : 1;
}
