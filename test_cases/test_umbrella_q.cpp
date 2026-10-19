/*
 * <Umbrella restraint on a distance difference (SimpleMD umbrella_q): gradient and parser.>
 * Copyright (C) 2026 Conrad Huebler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated 2026 (AI-generated, machine-tested).
 *
 * W = kappa/2 (q - q0)^2 with q = d(i,j) - d(i,k): analytic gradient vs central finite
 * difference (Angstrom displacements), W = 0 and grad = 0 at q = q0, parser rejects bad input.
 * umbrella_cv (sum, dihedral, diff synonym): the same checks on random 6-atom structures, the IUPAC
 * dihedral sign (vs gfnff_geometry.h), periodic wrap of phi - phi0, and the collinear guard.
 */
#include <cmath>
#include <functional>
#include <iostream>
#include <random>

#include "src/capabilities/simplemd.h"
#include "src/core/energy_calculators/ff_methods/gfnff_geometry.h"
#include "src/core/units.h"

// Claude Generated 2026: umbrella_cv (sum, dihedral, diff synonym) checks, see docs/UMBRELLA_Q.md
static int testUmbrellaCV(const std::function<void(bool, const std::string&)>& check)
{
    const double pi = CurcumaUnit::Constants::PI, deg = pi / 180.0;
    std::string err;
    // (i) central finite difference (Angstrom displacements) on several random 6-atom structures
    std::mt19937 rng(11);
    std::uniform_real_distribution<double> U(-1.8, 1.8);
    for (const std::string spec : { "sum,1,2,3,4,2.5,100", "sum,2,5,2,6,3.0,40", "dihedral,1,2,3,4,60,20",
             "dihedral,6,3,2,5,-150,35", "diff,1,2,3,0.2,80" }) {
        SimpleMD::UmbrellaCV u;
        check(SimpleMD::parseUmbrellaCV(spec, 6, u, err), "parse '" + spec + "'" + (err.empty() ? "" : " " + err));
        double maxrel = 0.0;
        for (int trial = 0; trial < 6; ++trial) {
            Geometry R(6, 3);
            for (int a = 0; a < 6; ++a)
                for (int c = 0; c < 3; ++c)
                    R(a, c) = U(rng);
            Geometry g, gp, gm;
            SimpleMD::umbrellaCVEnergy(u, R, g);
            const double h = 1e-5;
            for (int a = 0; a < 6; ++a)
                for (int c = 0; c < 3; ++c) {
                    Geometry Rp = R, Rm = R;
                    Rp(a, c) += h;
                    Rm(a, c) -= h;
                    const double fd = (SimpleMD::umbrellaCVEnergy(u, Rp, gp) - SimpleMD::umbrellaCVEnergy(u, Rm, gm)) / (2 * h);
                    maxrel = std::max(maxrel, std::abs(fd - g(a, c)) / std::max(g.cwiseAbs().maxCoeff(), 1e-12));
                }
        }
        std::cout << "  " << spec << ": max rel. error analytic vs FD over 6 structures: " << maxrel << std::endl;
        check(maxrel < 1e-6, spec + " gradient == central FD (rel. < 1e-6)");

        // (ii) W = 0 and gradient exactly 0 at the minimum
        Geometry R(6, 3), g;
        for (int a = 0; a < 6; ++a)
            for (int c = 0; c < 3; ++c)
                R(a, c) = U(rng);
        double cv = 0.0;
        SimpleMD::umbrellaCVEnergy(u, R, g, &cv);
        SimpleMD::UmbrellaCV u0 = u;
        u0.q0 = cv;
        const double W0 = SimpleMD::umbrellaCVEnergy(u0, R, g);
        check(W0 == 0.0 && g.cwiseAbs().maxCoeff() == 0.0, spec + ": W = 0 and gradient = 0 at q = q0");
    }

    // diff synonym gives exactly the umbrella_q numbers
    {
        SimpleMD::UmbrellaCV c;
        SimpleMD::UmbrellaQ q;
        SimpleMD::parseUmbrellaCV("diff,2,4,5,0.3,100", 6, c, err);
        SimpleMD::parseUmbrellaQ("2,4,5,0.3,100", 6, q, err);
        Geometry R = Geometry::Random(6, 3), g1, g2;
        const double w1 = SimpleMD::umbrellaCVEnergy(c, R, g1), w2 = SimpleMD::umbrellaQEnergy(q, R, g2);
        check(w1 == w2 && (g1 - g2).cwiseAbs().maxCoeff() == 0.0, "diff,... is bit-identical to umbrella_q");
    }

    // (iv) convention on a hand-built geometry: j = origin, k on +z, i on +x, l rotated by a about +z
    auto build = [](double a) {
        Geometry R(4, 3);
        R << 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.5, std::cos(a), std::sin(a), 1.5;
        return R;
    };
    for (const double a : { 0.0, 60.0, -60.0, 180.0, 120.0 }) {
        const Geometry R = build(a * deg);
        const double phi = SimpleMD::dihedralAngle(R, 0, 1, 2, 3) / deg;
        Matrix dummy;
        const double ref = GFNFF_Geometry::calculateDihedralAngle(R.row(0).transpose(), R.row(1).transpose(), R.row(2).transpose(),
                               R.row(3).transpose(), dummy, false) / deg;
        const double expect = a == 180.0 ? 180.0 : a;
        std::cout << "  rotation " << a << " deg: phi = " << phi << ", gfnff_geometry = " << ref << std::endl;
        check(std::abs(std::abs(std::remainder(phi - expect, 360.0))) < 1e-9 && std::abs(std::remainder(phi - ref, 360.0)) < 1e-9,
            "convention: rotation " + std::to_string(a) + " deg gives phi = " + std::to_string(expect) + " (IUPAC, = gfnff_geometry)");
    }
    // gauche+ in IUPAC terms: looking along j->k the front bond turns clockwise onto the back bond.
    // Viewer looks along +z; a positive rotation about +z appears clockwise to that viewer -> phi = +60.

    // (iii) wrap-around: phi0 = 175, phi = -175 -> dphi = +10 deg; also phi0 = 180 vs phi = -180 + eps
    {
        SimpleMD::UmbrellaCV u;
        SimpleMD::parseUmbrellaCV("dihedral,1,2,3,4,175,20", 4, u, err);
        const Geometry R = build(-175.0 * deg);
        Geometry g;
        double phi = 0.0, dq = 0.0;
        const double W = SimpleMD::umbrellaCVEnergy(u, R, g, &phi, &dq);
        const double Wexp = 0.5 * u.kappa * std::pow(10.0 * deg, 2);
        std::cout << "  wrap: phi = " << phi / deg << " dphi = " << dq / deg << " W = " << W << " expected " << Wexp << std::endl;
        check(std::abs(dq / deg - 10.0) < 1e-9 && std::abs(W - Wexp) < 1e-15, "wrap: phi0 175, phi -175 -> dphi = +10 deg");
        // gradient sign: moving along -grad must lower W, i.e. move phi towards 175 (decrease |dphi|)
        double phi2 = 0.0, dq2 = 0.0;
        Geometry g2;
        SimpleMD::umbrellaCVEnergy(u, R - 1e-4 * g / g.norm(), g2, &phi2, &dq2);
        check(std::abs(dq2) < std::abs(dq), "wrap: step along -grad reduces |dphi|");
        SimpleMD::UmbrellaCV v;
        SimpleMD::parseUmbrellaCV("dihedral,1,2,3,4,180,20", 4, v, err);
        for (const double a : { 179.9, -179.9 }) {
            const double Wa = SimpleMD::umbrellaCVEnergy(v, build(a * deg), g, &phi, &dq);
            check(std::abs(std::abs(dq / deg) - 0.1) < 1e-9 && std::abs(Wa - 0.5 * v.kappa * std::pow(0.1 * deg, 2)) < 1e-15,
                "wrap: phi0 180, phi " + std::to_string(a) + " -> |dphi| = 0.1 deg");
        }
        SimpleMD::parseUmbrellaCV("dihedral,1,2,3,4,-180,20", 4, v, err);
        SimpleMD::umbrellaCVEnergy(v, build(179.9 * deg), g, &phi, &dq);
        check(std::abs(dq / deg + 0.1) < 1e-9, "wrap: phi0 -180, phi 179.9 -> dphi = -0.1 deg");
    }

    // (vi) collinear guard: i-j-k linear
    {
        SimpleMD::UmbrellaCV u;
        SimpleMD::parseUmbrellaCV("dihedral,1,2,3,4,60,20", 4, u, err);
        Geometry R = build(60.0 * deg);
        R.row(0) << 0.0, 0.0, -1.0;
        Geometry g;
        bool thrown = false;
        try {
            SimpleMD::umbrellaCVEnergy(u, R, g);
        } catch (const std::runtime_error& e) {
            thrown = true;
            std::cout << "  guard message: " << e.what() << std::endl;
        }
        check(thrown, "collinear i-j-k raises an error");
        R = build(60.0 * deg);
        R.row(3) << 0.0, 0.0, 2.7;
        thrown = false;
        try {
            SimpleMD::umbrellaCVEnergy(u, R, g);
        } catch (const std::runtime_error&) {
            thrown = true;
        }
        check(thrown, "collinear j-k-l raises an error");
    }

    // (v) parser rejections
    SimpleMD::UmbrellaCV d;
    for (const std::string bad : { "sum,1,2,3,4,2.5", "sum,1,2,3,4,2.5,10,1", "dihedral,1,2,3,60,20", "angle,1,2,3,4,60,20",
             "sum", "", "sum,0,2,3,4,2.5,10", "sum,1,2,3,7,2.5,10", "sum,1,1,3,4,2.5,10", "sum,1,2,4,4,2.5,10",
             "dihedral,1,2,3,1,60,20", "dihedral,1,2,2,4,60,20", "sum,1,2,3,4,2.5,0", "dihedral,1,2,3,4,60,-1",
             "sum,1,2,a,4,2.5,10", "dihedral,1,2,3,4,x,20", "sum,1,2,3,4,2.5,10kcal", "diff,1,2,2,0.1,10",
             "diff,1,2,3,4,0.1,10" }) {
        const bool rejected = !SimpleMD::parseUmbrellaCV(bad, 6, d, err);
        check(rejected, "reject '" + bad + "'" + (rejected ? " (" + err + ")" : ""));
    }
    check(SimpleMD::parseUmbrellaCV("sum,1,2,1,3,2.5,10", 6, d, err), "sum accepts i == k (shared atom)");
    check(std::abs(d.kappa - 10.0 / CurcumaUnit::Energy::HARTREE_TO_KCALMOL) < 1e-15, "kappa kcal/mol -> Eh");
    check(SimpleMD::parseUmbrellaCV("dihedral,1,2,3,4,90,10", 6, d, err) && std::abs(d.q0 - pi / 2) < 1e-15, "phi0 degrees -> rad");
    return 0;
}

int main()
{
    int failures = 0;
    auto check = [&](bool ok, const std::string& w) {
        std::cout << (ok ? "  PASS  " : "  FAIL  ") << w << std::endl;
        if (!ok) ++failures;
    };
    std::mt19937 rng(7);
    std::uniform_real_distribution<double> U(-2.0, 2.0);
    Geometry R(5, 3);
    for (int a = 0; a < 5; ++a)
        for (int c = 0; c < 3; ++c)
            R(a, c) = U(rng);

    SimpleMD::UmbrellaQ u;
    std::string err;
    check(SimpleMD::parseUmbrellaQ("2,4,5,0.3,100", 5, u, err), "parse '2,4,5,0.3,100'");
    check(u.i == 1 && u.j == 3 && u.k == 4 && u.q0 == 0.3, "indices 1-based -> 0-based, q0 kept");
    check(std::abs(u.kappa - 100.0 / CurcumaUnit::Energy::HARTREE_TO_KCALMOL) < 1e-15, "kappa kcal/mol/A^2 -> Eh/A^2");

    // (1) gradient vs central finite difference
    Geometry g, gp, gm;
    SimpleMD::umbrellaQEnergy(u, R, g);
    const double h = 1e-5;
    double maxrel = 0.0;
    for (int a = 0; a < 5; ++a)
        for (int c = 0; c < 3; ++c) {
            Geometry Rp = R, Rm = R;
            Rp(a, c) += h;
            Rm(a, c) -= h;
            const double fd = (SimpleMD::umbrellaQEnergy(u, Rp, gp) - SimpleMD::umbrellaQEnergy(u, Rm, gm)) / (2 * h);
            maxrel = std::max(maxrel, std::abs(fd - g(a, c)) / std::max(g.cwiseAbs().maxCoeff(), 1e-12));
        }
    std::cout << "  max rel. error analytic vs FD gradient: " << maxrel << std::endl;
    check(maxrel < 1e-6, "gradient == central FD (rel. < 1e-6)");
    check(g.row(0).norm() == 0.0 && g.row(2).norm() == 0.0, "atoms outside i,j,k carry no gradient");

    // (2) W = 0 and grad exactly 0 at q = q0
    double q = 0.0;
    SimpleMD::umbrellaQEnergy(u, R, g, &q);
    SimpleMD::UmbrellaQ u0 = u;
    u0.q0 = q;
    const double W0 = SimpleMD::umbrellaQEnergy(u0, R, g);
    check(W0 == 0.0 && g.cwiseAbs().maxCoeff() == 0.0, "W = 0 and gradient = 0 at q = q0");

    // (3) parser rejects bad strings
    SimpleMD::UmbrellaQ d;
    for (const std::string bad : { "1,2,3,0.1", "1,2,3,0.1,10,5", "0,2,3,0.1,10", "1,2,6,0.1,10", "1,2,2,0.1,10",
             "1,1,3,0.1,10", "1,2,3,0.1,0", "1,2,3,0.1,-5", "1,a,3,0.1,10", "1,2,3,x,10", "1,2,3,0.1,10kcal", "" }) {
        const bool rejected = !SimpleMD::parseUmbrellaQ(bad, 5, d, err);
        check(rejected, "reject '" + bad + "'" + (rejected ? " (" + err + ")" : ""));
    }
    testUmbrellaCV(check);
    return failures == 0 ? 0 : 1;
}
