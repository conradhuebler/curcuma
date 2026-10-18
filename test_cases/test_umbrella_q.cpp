/*
 * <Umbrella restraint on a distance difference (SimpleMD umbrella_q): gradient and parser.>
 * Copyright (C) 2026 Conrad Huebler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated 2026 (AI-generated, machine-tested).
 *
 * W = kappa/2 (q - q0)^2 with q = d(i,j) - d(i,k): analytic gradient vs central finite
 * difference (Angstrom displacements), W = 0 and grad = 0 at q = q0, parser rejects bad input.
 */
#include <cmath>
#include <iostream>
#include <random>

#include "src/capabilities/simplemd.h"
#include "src/core/units.h"

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
    return failures == 0 ? 0 : 1;
}
