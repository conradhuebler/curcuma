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

#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <string>
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

    // (f) atom subset (cv_atoms), Claude Generated 2026. A random 10-atom molecule
    // (C O C N H H H C H H), a random smooth 8-frame path, a test point displaced off the path
    // and rotated. Checks: (i) s, z do not change when unselected atoms move; (ii) analytic
    // lab-frame ds/dR, dz/dR vs central FD (Angstrom displacements), rel < 1e-6; (iii) the
    // gradient rows of unselected atoms are exactly 0; (iv) the full explicit list equals the
    // 'all' code path; (v) s, z equal an independent subset-only reference (a PathCV built on
    // the selected rows alone); parser checks.
    {
        const int Z[10] = { 6, 8, 6, 7, 1, 1, 1, 6, 1, 1 };
        std::mt19937 gen(1234);
        std::normal_distribution<double> nd(0.0, 1.0);
        Molecule m10;
        Geometry X0(10, 3);
        for (int a = 0; a < 10; ++a) {
            Position q;
            q << 1.2 * nd(gen), 1.2 * nd(gen), 1.2 * nd(gen);
            X0.row(a) = q.transpose();
            m10.addPair({ Z[a], q });
        }
        Geometry D1(10, 3), D2(10, 3);
        for (int a = 0; a < 10; ++a)
            for (int c = 0; c < 3; ++c) {
                D1(a, c) = 0.5 * nd(gen);
                D2(a, c) = 0.3 * nd(gen);
            }
        std::vector<Geometry> p10;
        for (int i = 0; i < 8; ++i) {
            const double f = i / 7.0;
            p10.push_back(X0 + f * D1 + std::sin(M_PI * f) * D2);
        }
        const double th = 35.0 * M_PI / 180.0;
        Eigen::Matrix3d Ry;
        Ry << std::cos(th), 0, std::sin(th), 0, 1, 0, -std::sin(th), 0, std::cos(th);
        Geometry Rt = 0.5 * (p10[3] + p10[4]);
        for (int a = 0; a < 10; ++a)
            for (int c = 0; c < 3; ++c)
                Rt(a, c) += 0.08 * nd(gen);
        Rt = (Rt * Ry.transpose()).eval();
        Rt.rowwise() += Eigen::RowVector3d(0.4, -0.3, 0.2);

        // parser
        std::vector<int> heavy, list, tmp;
        std::string err;
        check(PathCV::parseAtomSelection("heavy", m10, heavy, err)
                  && heavy == std::vector<int>({ 0, 1, 2, 3, 7 }),
              "cv_atoms 'heavy' = atoms with Z > 1 (0-based 0,1,2,3,7)");
        check(PathCV::parseAtomSelection("1-6,9", m10, list, err)
                  && list == std::vector<int>({ 0, 1, 2, 3, 4, 5, 8 }),
              "cv_atoms '1-6,9' = 0-based 0..5,8");
        check(PathCV::parseAtomSelection("all", m10, tmp, err) && tmp.empty(), "cv_atoms 'all' = empty selection");
        check(PathCV::parseAtomSelection(" 9, 2 ,4-5 ", m10, tmp, err) && tmp == std::vector<int>({ 1, 3, 4, 8 }),
              "cv_atoms list with blanks, unsorted, is sorted");
        for (const std::string bad : { "0", "11", "1-3,3", "1,2", "x", "3-1", "1-", "2-11", "1,,3" }) {
            const bool ok = PathCV::parseAtomSelection(bad, m10, tmp, err);
            std::cout << "   '" << bad << "' -> " << (ok ? "accepted" : err) << std::endl;
            check(!ok && !err.empty(), "cv_atoms '" + bad + "' is rejected with a message");
        }
        {
            Molecule eth; // C2H6: 'heavy' = 2 atoms -> rejected
            Position q = Position::Zero();
            for (int a = 0; a < 8; ++a) {
                q(0) = a;
                eth.addPair({ a < 2 ? 6 : 1, q });
            }
            check(!PathCV::parseAtomSelection("heavy", eth, tmp, err), "cv_atoms 'heavy' on ethane (2 atoms) is rejected");
        }

        const PathCV cvall(p10, m10);
        std::vector<int> full(10);
        for (int a = 0; a < 10; ++a)
            full[a] = a;
        const PathCV cvfull(p10, m10, -1.0, full);
        {
            double s1, z1, s2, z2;
            Geometry ds1, dz1, ds2, dz2;
            cvall.evaluate(Rt, &s1, &z1, &ds1, &dz1);
            cvfull.evaluate(Rt, &s2, &z2, &ds2, &dz2);
            const double dmax = std::max({ std::abs(s1 - s2), std::abs(z1 - z2), (ds1 - ds2).cwiseAbs().maxCoeff(),
                                           (dz1 - dz2).cwiseAbs().maxCoeff(), std::abs(cvall.lambda() - cvfull.lambda()) });
            std::cout << "   all vs full explicit list: max |diff| = " << dmax << std::endl;
            check(dmax < 1e-12, "(iv) full explicit list reproduces the 'all' code path");
        }

        for (const auto* sel : { &heavy, &list }) {
            const std::string tag = (sel == &heavy) ? "heavy" : "1-6,9";
            const PathCV cv(p10, m10, -1.0, *sel);
            double s0 = 0, z0 = 0;
            Geometry ds, dz;
            cv.evaluate(Rt, &s0, &z0, &ds, &dz);

            // (v) independent reference: PathCV built on the selected rows only
            std::vector<Geometry> psub;
            Molecule msub;
            for (const auto& g : p10) {
                Geometry h(static_cast<int>(sel->size()), 3);
                for (size_t k = 0; k < sel->size(); ++k)
                    h.row(static_cast<int>(k)) = g.row((*sel)[k]);
                psub.push_back(h);
            }
            Geometry Rsub(static_cast<int>(sel->size()), 3);
            for (size_t k = 0; k < sel->size(); ++k) {
                Rsub.row(static_cast<int>(k)) = Rt.row((*sel)[k]);
                const Position q = Rt.row((*sel)[k]).transpose();
                msub.addPair({ Z[(*sel)[k]], q });
            }
            const PathCV cvsub(psub, msub);
            double ss = 0, zs = 0;
            cvsub.evaluate(Rsub, &ss, &zs, nullptr, nullptr);
            check(std::abs(ss - s0) < 1e-12 && std::abs(zs - z0) < 1e-12 && std::abs(cv.lambda() - cvsub.lambda()) < 1e-9,
                  "(v) " + tag + ": s, z, lambda equal a CV built on the selected atoms alone");

            // (i) invariance under moving unselected atoms
            Geometry Rmov = Rt;
            for (int a = 0; a < 10; ++a)
                if (std::find(sel->begin(), sel->end(), a) == sel->end())
                    for (int c = 0; c < 3; ++c)
                        Rmov(a, c) += 0.7 * nd(gen);
            double sm = 0, zm = 0;
            cv.evaluate(Rmov, &sm, &zm, nullptr, nullptr);
            check(sm == s0 && zm == z0, "(i) " + tag + ": s, z unchanged when unselected atoms move");

            // (ii) FD, (iii) zero rows
            const double h = 1e-5;
            double es = 0, ez = 0, ns = 0, nz = 0;
            bool zero_rows = true;
            for (int a = 0; a < 10; ++a) {
                const bool selected = std::find(sel->begin(), sel->end(), a) != sel->end();
                if (!selected && (ds.row(a).norm() != 0.0 || dz.row(a).norm() != 0.0))
                    zero_rows = false;
                for (int c = 0; c < 3; ++c) {
                    Geometry Rp = Rt, Rm = Rt;
                    Rp(a, c) += h;
                    Rm(a, c) -= h;
                    double sp, zp, sq, zq;
                    cv.evaluate(Rp, &sp, &zp, nullptr, nullptr);
                    cv.evaluate(Rm, &sq, &zq, nullptr, nullptr);
                    const double fs = (sp - sq) / (2 * h), fz = (zp - zq) / (2 * h);
                    es += (fs - ds(a, c)) * (fs - ds(a, c));
                    ez += (fz - dz(a, c)) * (fz - dz(a, c));
                    ns += fs * fs;
                    nz += fz * fz;
                }
            }
            const double rs = std::sqrt(es / std::max(ns, 1e-30)), rz = std::sqrt(ez / std::max(nz, 1e-30));
            std::cout << "   " << tag << " (" << sel->size() << " atoms): s = " << s0 << ", z = " << z0
                      << ", lambda = " << cv.lambda() << ", FD rel err ds = " << rs << ", dz = " << rz << std::endl;
            check(rs < 1e-6 && rz < 1e-6 && ns > 0.0 && nz > 0.0, "(ii) " + tag + ": ds/dR, dz/dR vs FD rel < 1e-6");
            check(zero_rows, "(iii) " + tag + ": gradient rows of unselected atoms exactly 0");
        }
    }

    // Distance-space metric (cv_metric distances), Claude Generated 2026. A random 8-atom molecule
    // (C O C N H H C H), a smooth 10-frame path with a dominant linear part, CV pairs = heavy-atom
    // pairs whose distance spans > 0.2 A. Blocks:
    //   (dist-d) s strictly monotone over the path frames (so s(i) lies between s(i-1) and s(i+1))
    //   (dist-e) analytic ds/dR, dz/dR vs central FD (1e-4 A) at an off-path random point, rel < 1e-6,
    //            for the transforms d and inv
    //   (dist-f) s, z invariant under a random rotation + translation (< 1e-12); gradient covariant
    //   (dist-g) pair-list parser (valid/invalid) and the 'auto' rule on a path with known spans
    {
        const int Z[8] = { 6, 8, 6, 7, 1, 1, 6, 1 };
        std::mt19937 gen(2718);
        std::normal_distribution<double> nd(0.0, 1.0);
        Molecule m8;
        Geometry X0(8, 3), D1(8, 3), D2(8, 3);
        for (int a = 0; a < 8; ++a) {
            Position q;
            q << 1.3 * nd(gen), 1.3 * nd(gen), 1.3 * nd(gen);
            X0.row(a) = q.transpose();
            m8.addPair({ Z[a], q });
            for (int c = 0; c < 3; ++c) {
                D1(a, c) = 0.9 * nd(gen);
                D2(a, c) = 0.15 * nd(gen);
            }
        }
        std::vector<Geometry> p8;
        for (int i = 0; i < 10; ++i) {
            const double f = i / 9.0;
            p8.push_back(X0 + f * D1 + std::sin(M_PI * f) * D2);
        }
        const std::vector<int> heavy8 = { 0, 1, 2, 3, 6 };
        std::vector<double> spans;
        const auto pairs = PathCV::selectPairsBySpan(p8, heavy8, 0.2, &spans);
        std::cout << "   distance metric: " << pairs.size() << " heavy-atom pairs with span > 0.2 A" << std::endl;
        check(pairs.size() >= 2, "(dist) at least 2 auto pairs on the random path");

        for (auto tr : { PathCV::PairTransform::Distance, PathCV::PairTransform::Inverse }) {
            const std::string tag = (tr == PathCV::PairTransform::Inverse) ? "inv" : "d";
            const PathCV cvd(p8, m8, -1.0, {}, PathCV::Metric::Distances, pairs, tr);
            // lambda = 1/<d_k,k+1^2>
            double m2 = 0.0;
            for (int k = 1; k < 10; ++k)
                m2 += PathCV::pairSpacing2(p8[k], p8[k - 1], pairs, tr);
            m2 /= 9.0;
            check(std::abs(cvd.lambda() * m2 - 1.0) < 1e-12, "(dist) " + tag + ": lambda = 1/<d_k,k+1^2>");

            // (dist-d)
            bool mono = true;
            std::vector<double> sv(10);
            for (int i = 0; i < 10; ++i)
                sv[i] = cvd.s(p8[i]);
            for (int i = 1; i < 10; ++i)
                if (!(sv[i] > sv[i - 1]))
                    mono = false;
            bool between = true;
            for (int i = 1; i < 9; ++i)
                if (!(sv[i] > sv[i - 1] && sv[i] < sv[i + 1]))
                    between = false;
            std::cout << "   " << tag << ": s(frames) = " << sv.front() << " .. " << sv.back() << ", lambda = " << cvd.lambda()
                      << std::endl;
            check(mono && between, "(dist-d) " + tag + ": s strictly monotone, s(i) between s(i-1) and s(i+1)");
            // nearest frame and z consistency
            {
                double z1 = 0, z2 = 0;
                cvd.evaluate(p8[4], nullptr, &z1, nullptr, nullptr);
                const int nn = cvd.nearestFrame(p8[4], &z2);
                check(nn == 4 && z1 == z2, "(dist) " + tag + ": nearestFrame(frame 4) = 4, z identical to evaluate()");
            }

            // (dist-e) FD at an off-path random point
            Geometry Rt = 0.5 * (p8[4] + p8[5]);
            for (int a = 0; a < 8; ++a)
                for (int c = 0; c < 3; ++c)
                    Rt(a, c) += 0.1 * nd(gen);
            double s0 = 0, z0 = 0;
            Geometry ds, dz;
            cvd.evaluate(Rt, &s0, &z0, &ds, &dz);
            const double h = 1e-4;
            double es = 0, ez = 0, ns = 0, nz = 0;
            for (int a = 0; a < 8; ++a)
                for (int c = 0; c < 3; ++c) {
                    Geometry Rp = Rt, Rm = Rt;
                    Rp(a, c) += h;
                    Rm(a, c) -= h;
                    double sp, zp, sq, zq;
                    cvd.evaluate(Rp, &sp, &zp, nullptr, nullptr);
                    cvd.evaluate(Rm, &sq, &zq, nullptr, nullptr);
                    const double fs = (sp - sq) / (2 * h), fz = (zp - zq) / (2 * h);
                    es += (fs - ds(a, c)) * (fs - ds(a, c));
                    ez += (fz - dz(a, c)) * (fz - dz(a, c));
                    ns += fs * fs;
                    nz += fz * fz;
                }
            const double rs = std::sqrt(es / std::max(ns, 1e-30)), rz = std::sqrt(ez / std::max(nz, 1e-30));
            std::cout << "   " << tag << ": off-path s = " << s0 << ", z = " << z0 << " " << cvd.zUnit()
                      << ", FD rel err ds = " << rs << ", dz = " << rz << std::endl;
            check(rs < 1e-6 && rz < 1e-6 && ns > 0.0 && nz > 0.0, "(dist-e) " + tag + ": ds/dR, dz/dR vs FD rel < 1e-6");
            // rows of atoms in no pair are exactly zero (atom 7 is H, never in a heavy pair)
            check(ds.row(7).norm() == 0.0 && dz.row(7).norm() == 0.0, "(dist-e) " + tag + ": gradient of an atom in no pair is 0");

            // (dist-f) rigid motion
            Eigen::Matrix3d A;
            for (int r = 0; r < 3; ++r)
                for (int c = 0; c < 3; ++c)
                    A(r, c) = nd(gen);
            Eigen::HouseholderQR<Eigen::Matrix3d> qr(A);
            Eigen::Matrix3d Q = qr.householderQ();
            if (Q.determinant() < 0)
                Q.col(0) *= -1.0;
            Geometry Rr = (Rt * Q.transpose()).eval();
            Rr.rowwise() += Eigen::RowVector3d(1.7, -0.4, 2.2);
            double s1 = 0, z1 = 0;
            Geometry ds1, dz1;
            cvd.evaluate(Rr, &s1, &z1, &ds1, &dz1);
            const double gerr = std::max((ds1 - ds * Q.transpose()).cwiseAbs().maxCoeff(),
                                         (dz1 - dz * Q.transpose()).cwiseAbs().maxCoeff());
            std::cout << "   " << tag << ": rigid motion |ds| = " << std::abs(s1 - s0) << ", |dz| = " << std::abs(z1 - z0)
                      << ", gradient covariance err = " << gerr << std::endl;
            check(std::abs(s1 - s0) < 1e-12 && std::abs(z1 - z0) < 1e-12, "(dist-f) " + tag + ": s, z invariant under rotation + translation");
            check(gerr < 1e-10, "(dist-f) " + tag + ": gradient rotates covariantly (grad(QR) = grad(R) Q^T)");
        }

        // resampling at equal distance-space arclength: chord spacing nearly uniform
        {
            const auto rs = PathCV::resampleByPairs(p8, 15, pairs, PathCV::PairTransform::Distance);
            double dmin = 1e300, dmax = 0;
            for (int k = 1; k < 15; ++k) {
                const double d = std::sqrt(PathCV::pairSpacing2(rs[k], rs[k - 1], pairs, PathCV::PairTransform::Distance));
                dmin = std::min(dmin, d);
                dmax = std::max(dmax, d);
            }
            std::cout << "   resampleByPairs 10 -> 15 frames: spacing min " << dmin << " max " << dmax << std::endl;
            check(rs.size() == 15 && (rs.front() - p8.front()).norm() == 0.0 && (rs.back() - p8.back()).norm() == 0.0
                      && dmax / dmin < 1.02,
                  "(dist) resampleByPairs: endpoints kept, distance-space spacing uniform within 2 %");
        }

        // (dist-g) parser
        std::vector<PathCV::AtomPair> pl;
        std::string err;
        check(PathCV::parsePairList("1-4,1-9,4-9", 10, pl, err)
                  && pl == std::vector<PathCV::AtomPair>({ { 0, 3 }, { 0, 8 }, { 3, 8 } }),
              "(dist-g) cv_pairs '1-4,1-9,4-9' = pairs (0,3),(0,8),(3,8) (a pair, not a range)");
        check(PathCV::parsePairList(" 4-1 , 2 - 3 ", 10, pl, err) && pl == std::vector<PathCV::AtomPair>({ { 0, 3 }, { 1, 2 } }),
              "(dist-g) blanks and reversed order are normalised");
        for (const std::string bad : { "1-4,4-1", "1-4,1-4", "1-1,2-3", "1-11,2-3", "0-2,2-3", "1-4", "", "1-4,x",
                                       "1-4,1", "1-4,1-2-3", "1-4,,2-3", "a-b,1-2", "1-4,2-" }) {
            const bool ok = PathCV::parsePairList(bad, 10, pl, err);
            std::cout << "   '" << bad << "' -> " << (ok ? "accepted" : err) << std::endl;
            check(!ok && !err.empty() && pl.empty(), "(dist-g) cv_pairs '" + bad + "' is rejected with a message");
        }
        // 'auto' on a path with known spans: C0 fixed, C1 (1,0,0)->(2,0,0), C2 (0,3,0)->(0,3.8,0),
        // H3 (0,0,1)->(0,0,3). Spans (all monotone): 0-1 1.0, 0-2 0.8, 1-2 sqrt(4+14.44)-sqrt(10) = 1.1325,
        // 0-3 2.0 (H).
        {
            std::vector<Geometry> pk;
            for (int i = 0; i < 5; ++i) {
                const double t = i / 4.0;
                Geometry g = Geometry::Zero(4, 3);
                g(1, 0) = 1.0 + t;
                g(2, 1) = 3.0 + 0.8 * t;
                g(3, 2) = 1.0 + 2.0 * t;
                pk.push_back(g);
            }
            std::vector<double> sp;
            const auto h05 = PathCV::selectPairsBySpan(pk, { 0, 1, 2 }, 0.5, &sp);
            const double s12 = std::sqrt(4.0 + 3.8 * 3.8) - std::sqrt(10.0);
            check(h05 == std::vector<PathCV::AtomPair>({ { 0, 1 }, { 0, 2 }, { 1, 2 } }) && sp.size() == 3
                      && std::abs(sp[0] - 1.0) < 1e-12 && std::abs(sp[1] - 0.8) < 1e-12 && std::abs(sp[2] - s12) < 1e-12,
                  "(dist-g) auto, heavy atoms, threshold 0.5: pairs 1-2, 1-3, 2-3 with spans 1.0, 0.8, 1.1325");
            const auto h09 = PathCV::selectPairsBySpan(pk, { 0, 1, 2 }, 0.9, &sp);
            check(h09 == std::vector<PathCV::AtomPair>({ { 0, 1 }, { 1, 2 } }), "(dist-g) auto, threshold 0.9: pairs 1-2, 2-3");
            const auto all15 = PathCV::selectPairsBySpan(pk, { 0, 1, 2, 3 }, 1.5, &sp);
            check(all15.size() >= 1 && all15.front() == PathCV::AtomPair(0, 3),
                  "(dist-g) auto with H candidates, threshold 1.5: H pair 1-4 (span 2.0) selected");
            const auto none = PathCV::selectPairsBySpan(pk, { 0, 1, 2 }, 1.2, &sp);
            check(none.empty(), "(dist-g) auto, threshold 1.2: no heavy pair (caller reports < 2 pairs as an error)");
        }
    }

    std::cout << (failures == 0 ? "ALL TESTS PASSED" : "SOME TESTS FAILED") << std::endl;
    return failures == 0 ? 0 : 1;
}
