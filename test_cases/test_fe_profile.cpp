/*
 * <Unit test for the grid-free (natural cubic spline) integration of a mean-force profile.>
 * Copyright (C) 2019 - 2026 Conrad Huebler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated 2026 (AI-generated, machine-tested).
 *
 * (1) Barrier G(s) = 10 kcal/mol (1 - x^2) + 1.5 kcal/mol x, x = (s - 0.5)/0.5, sampled by its
 *     exact dG/ds on 6, 10 and 22 equidistant support points in [0.05, 0.95]. The spline
 *     barrier must lie within 0.05 kcal/mol of the exact value (10, 22 points) and within 0.15
 *     (6 points); the support-point maximum is printed next to it. Both barrier definitions
 *     (max G - G(first) and max G - mean of the end values) are checked.
 *     Note: dG/ds is linear here, which a natural spline reproduces exactly.
 * (2) Non-polynomial profile G(s) = 5 kcal/mol (1 - cos(2 pi s)) + 1.5 kcal/mol s (dG/ds not
 *     a polynomial): spline barrier within 0.05 kcal/mol at 10 and 22 points.
 * (3) Symmetric profile (no tilt): reaction dG = 0 +- 1e-6 kcal/mol.
 * (4) Fewer than 3 support points: trapezoid fallback with a warning.
 * (5) Basin barrier (basinBarrier): double well G = 5 (u^2-1)^2 + c s kcal/mol, u = (s-0.5)/0.3,
 *     exact dG/ds on 6/10/22 points in [0,1], spline-integrated: c = 0 -> forward = reverse
 *     (+- 1e-3); c = 1.5 -> fwd - rev = reaction (+- 1e-6); n >= 10 -> within 0.05 kcal/mol of the
 *     same basin integrals on the exact G (2e6+1 points). basinCoverage: [0,1] covered, flanks
 *     [0.3,0.7] not covered; maximum at the edge -> invalid.
 * (6) splitMeanForce: gradE + metric + wall = total dG/ds to 1e-12 Eh (both estimator-I forms);
 *     passive wall -> wall part exactly 0.
 * Optional arguments: stage-2 *.neb.fe.csv files; the estimator-I rows are integrated and
 * both barrier definitions printed (offline cross-check, not part of the ctest run).
 */

#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "src/capabilities/fe_profile.h"

static int failures = 0;
static void check(bool ok, const std::string& what)
{
    std::cout << (ok ? "  PASS  " : "  FAIL  ") << what << std::endl;
    if (!ok)
        ++failures;
}

static const double kcal = 627.5094740631;
static const double PI = 3.14159265358979323846;

struct Case {
    double node_fwd, node_mean, spl_fwd, spl_mean, reaction;
};

static Case run(int n, const std::function<double(double)>& dG_kcal)
{
    std::vector<double> s(n), y(n);
    for (int k = 0; k < n; ++k) {
        s[k] = 0.05 + 0.9 * k / (n - 1.0);
        y[k] = dG_kcal(s[k]) / kcal;
    }
    const auto p = curcuma::fe_profile::integrateSpline(s, y);
    const auto Gn = curcuma::fe_profile::integrateTrapezoid(s, y);
    double gm = Gn[0];
    for (double g : Gn)
        gm = std::max(gm, g);
    return { (gm - Gn.front()) * kcal, (gm - 0.5 * (Gn.front() + Gn.back())) * kcal, p.barrier_forward * kcal,
        p.barrier_mean_ends * kcal, p.reaction * kcal };
}

static void csvCheck(const std::string& file)
{
    std::ifstream in(file);
    std::string line;
    std::vector<double> s, y, gk;
    while (std::getline(in, line)) {
        if (line.rfind("I,", 0) != 0)
            continue;
        std::stringstream ss(line.substr(2));
        std::string a, b, c;
        std::getline(ss, a, ',');
        std::getline(ss, b, ',');
        std::getline(ss, c, ',');
        s.push_back(std::stod(a));
        y.push_back(std::stod(b));
        gk.push_back(std::stod(c));
    }
    if (s.size() < 2) {
        std::cout << "  " << file << ": fewer than 2 estimator-I rows" << std::endl;
        return;
    }
    const auto p = curcuma::fe_profile::integrateSpline(s, y);
    double gm = gk[0];
    for (double g : gk)
        gm = std::max(gm, g);
    std::printf("  CSV %s\n      n=%zu  spline: fwd %.4f  mean-ends %.4f  rxn %.4f  s_max %.4f | "
                "nodes (G_kcal column): fwd %.4f  mean-ends %.4f kcal/mol\n",
        file.c_str(), s.size(), p.barrier_forward * kcal, p.barrier_mean_ends * kcal, p.reaction * kcal, p.s_max,
        gm - gk.front(), gm - 0.5 * (gk.front() + gk.back()));
}

int main(int argc, char** argv)
{
    // (1) tilted barrier, exact values
    auto dG1 = [](double s) { const double x = (s - 0.5) / 0.5; return (-20.0 * x + 1.5) * 2.0; };
    auto G1 = [](double s) { const double x = (s - 0.5) / 0.5; return 10.0 * (1.0 - x * x) + 1.5 * x; };
    const double xm = 0.075, gmax1 = 10.0 * (1.0 - xm * xm) + 1.5 * xm;
    const double ex_fwd = gmax1 - G1(0.05), ex_mean = gmax1 - 0.5 * (G1(0.05) + G1(0.95));
    std::printf("(1) G = 10(1-x^2) + 1.5x kcal/mol: exact barrier fwd %.5f, mean-ends %.5f\n", ex_fwd, ex_mean);
    for (int n : { 6, 10, 22 }) {
        const Case c = run(n, dG1);
        const double tol = (n == 6) ? 0.15 : 0.05;
        std::printf("    n=%2d  nodes: fwd %.5f (err %+.5f)  mean %.5f (err %+.5f) | spline: fwd %.5f (err %+.5f)  mean %.5f (err %+.5f)\n",
            n, c.node_fwd, c.node_fwd - ex_fwd, c.node_mean, c.node_mean - ex_mean, c.spl_fwd, c.spl_fwd - ex_fwd,
            c.spl_mean, c.spl_mean - ex_mean);
        check(std::abs(c.spl_fwd - ex_fwd) <= tol && std::abs(c.spl_mean - ex_mean) <= tol,
            "spline barrier within " + std::to_string(tol) + " kcal/mol at n=" + std::to_string(n));
        check(std::abs(c.reaction - (G1(0.95) - G1(0.05))) <= 1e-6, "reaction dG exact at n=" + std::to_string(n));
        if (n == 6)
            check(std::abs(c.node_fwd - ex_fwd) > 3.0 * std::abs(c.spl_fwd - ex_fwd) + 0.05,
                "n=6: support-point maximum clearly worse than the spline");
    }
    // (2) non-polynomial profile
    auto dG2 = [](double s) { return 5.0 * 2.0 * PI * std::sin(2.0 * PI * s) + 1.5; };
    auto G2 = [](double s) { return 5.0 * (1.0 - std::cos(2.0 * PI * s)) + 1.5 * s; };
    // maximum: 10 pi sin(2 pi s) = -1.5 -> s slightly above 0.5
    const double s2 = 0.5 + std::asin(1.5 / (10.0 * PI)) / (2.0 * PI);
    const double ex2 = G2(s2) - G2(0.05);
    std::printf("(2) G = 5(1-cos 2 pi s) + 1.5 s kcal/mol: exact barrier fwd %.5f\n", ex2);
    for (int n : { 6, 10, 22 }) {
        const Case c = run(n, dG2);
        std::printf("    n=%2d  nodes fwd %.5f (err %+.5f) | spline fwd %.5f (err %+.5f)\n", n, c.node_fwd,
            c.node_fwd - ex2, c.spl_fwd, c.spl_fwd - ex2);
        if (n >= 10)
            check(std::abs(c.spl_fwd - ex2) <= 0.05, "non-polynomial: spline barrier within 0.05 at n=" + std::to_string(n));
    }
    // (3) symmetric profile -> reaction dG = 0
    auto dG3 = [](double s) { const double x = (s - 0.5) / 0.5; return -20.0 * x * 2.0; };
    for (int n : { 6, 10, 22 }) {
        const Case c = run(n, dG3);
        check(std::abs(c.reaction) <= 1e-6, "symmetric profile: reaction dG " + std::to_string(c.reaction)
                + " = 0 +- 1e-6 at n=" + std::to_string(n));
    }
    // (4) fallback
    const auto f = curcuma::fe_profile::integrateSpline({ 0.2, 0.8 }, { 1.0, -1.0 });
    check(!f.spline && !f.warning.empty() && std::abs(f.reaction) < 1e-15, "2 support points -> trapezoid fallback with warning");

    // (5) basin free energies on a double well (Claude Generated 2026)
    //     G = A (u^2 - 1)^2 + c s, u = (s - 0.5)/0.3 -> minima near s = 0.2 / 0.8, maximum near 0.5.
    //     Support points: exact dG/ds on n = 6/10/22 equidistant points in [0, 1] (the walls at the
    //     ends lie ~3 A above the minima, so both basins are inside the domain). Reference: the
    //     same basin integrals on the EXACT G at 2,000,001 points.
    {
        const double A = 5.0, kT = 0.0019872043 * 300.0; // kcal/mol
        for (double c : { 0.0, 1.5 }) {
            auto Gx = [A, c](double s) { const double u = (s - 0.5) / 0.3; return A * (u * u - 1.0) * (u * u - 1.0) + c * s; };
            auto dGx = [A, c](double s) { const double u = (s - 0.5) / 0.3; return 4.0 * A * u * (u * u - 1.0) / 0.3 + c; };
            const int nf = 2000001;
            std::vector<double> sf(nf), gf(nf);
            for (int i = 0; i < nf; ++i) {
                sf[i] = static_cast<double>(i) / (nf - 1);
                gf[i] = Gx(sf[i]);
            }
            const auto ref = curcuma::fe_profile::basinBarrier(sf, gf, kT);
            std::printf("(5) double well A=%.1f, tilt c=%.1f kcal/mol: reference basin barrier fwd %.5f rev %.5f, "
                        "reaction %.5f (s* = %.5f, G(s*) - G(0.2) = %.5f)\n",
                A, c, ref.barrier_fwd, ref.barrier_rev, ref.reaction, ref.s_star, ref.g_star - Gx(0.2));
            for (int n : { 6, 10, 22 }) {
                std::vector<double> s(n), y(n);
                for (int k = 0; k < n; ++k) {
                    s[k] = static_cast<double>(k) / (n - 1.0);
                    y[k] = dGx(s[k]);
                }
                const auto p = curcuma::fe_profile::integrateSpline(s, y); // units: kcal/mol here
                const auto b = curcuma::fe_profile::basinBarrier(p.s, p.G, kT);
                const auto cov = curcuma::fe_profile::basinCoverage(s, y);
                std::printf("    n=%2d  basin: fwd %.5f (err %+.5f)  rev %.5f (err %+.5f)  rxn %.6f | fwd-rev-rxn %.2e | "
                            "covered L/R %d/%d\n",
                    n, b.barrier_fwd, b.barrier_fwd - ref.barrier_fwd, b.barrier_rev, b.barrier_rev - ref.barrier_rev,
                    b.reaction, b.barrier_fwd - b.barrier_rev - b.reaction, cov.left ? 1 : 0, cov.right ? 1 : 0);
                const std::string tag = " (c=" + std::to_string(c).substr(0, 3) + ", n=" + std::to_string(n) + ")";
                check(b.valid, "basin barrier valid" + tag);
                if (c == 0.0)
                    check(std::abs(b.barrier_fwd - b.barrier_rev) <= 1e-3, "symmetric: forward = reverse +- 1e-3" + tag);
                else
                    check(std::abs((b.barrier_fwd - b.barrier_rev) - b.reaction) <= 1e-6,
                        "tilted: fwd - rev = reaction dG (basin) +- 1e-6" + tag);
                if (n >= 10)
                    check(std::abs(b.barrier_fwd - ref.barrier_fwd) <= 0.05 && std::abs(b.barrier_rev - ref.barrier_rev) <= 0.05,
                        "basin barriers within 0.05 kcal/mol of the fine reference" + tag);
                check(cov.left && cov.right, "support points on [0,1] cover both basins" + tag);
            }
        }
        // coverage: support points only on the barrier flanks [0.3, 0.7] -> basins not reached
        std::vector<double> s(10), y(10);
        for (int k = 0; k < 10; ++k) {
            s[k] = 0.3 + 0.4 * k / 9.0;
            const double u = (s[k] - 0.5) / 0.3;
            y[k] = 4.0 * 5.0 * u * (u * u - 1.0) / 0.3;
        }
        const auto cov = curcuma::fe_profile::basinCoverage(s, y);
        std::printf("    flanks only [0.3,0.7]: slope L %.3f R %.3f max %.3f -> covered L/R %d/%d\n", cov.slope_left,
            cov.slope_right, cov.max_abs, cov.left ? 1 : 0, cov.right ? 1 : 0);
        check(!cov.left && !cov.right, "support points on the flanks only -> 'basin not covered' on both sides");
        // edge maximum -> invalid
        const auto e = curcuma::fe_profile::basinBarrier({ 0.0, 0.5, 1.0 }, { 0.0, 1.0, 2.0 }, 0.6);
        check(!e.valid && !e.warning.empty(), "maximum at the edge -> invalid with warning");
    }

    // (6) Claude Generated 2026: splitMeanForce (fe.csv columns dGds_gradE/metric/wall). The three
    //     parts must add up to the total exactly as NebMD forms it, in both estimator-I paths:
    //     constrained (wf - kT wd) / w and binned wf/w - kT wd/w, to 1e-12 Eh, for a passive wall
    //     (sum_wwall = 0 -> wall column 0, gradE = wf/w) and an active wall of either sign.
    {
        std::cout << "(6) mean-force split gradE + metric + wall = total" << std::endl;
        const double kT = 9.500366e-4; // 300 K in Eh
        const double cases[][4] = { // sum_wf, sum_wwall, sum_wd, sum_w (weighted sums of one support point)
            { 0.0, 0.0, 0.0, 1.0 }, { 12.345678901, 0.0, -3.21, 517.25 }, { -4.2e2, 3.7e1, 1.9e3, 2.4e4 },
            { 7.77e-3, -1.5e-3, 4.4e2, 3.3e-1 }, { 1.0e5, 9.99e4, -2.0e5, 7.0e3 }
        };
        double worst = 0.0;
        bool wall_zero_ok = true;
        for (const auto& c : cases) {
            const auto sp = curcuma::fe_profile::splitMeanForce(c[0], c[1], c[2], c[3], kT);
            const double sum = sp.grad_e + sp.metric + sp.wall;
            worst = std::max(worst, std::abs(sum - (c[0] - kT * c[2]) / c[3]));
            worst = std::max(worst, std::abs(sum - (c[0] / c[3] - kT * c[2] / c[3])));
            if (c[1] == 0.0 && (sp.wall != 0.0 || sp.grad_e != c[0] / c[3]))
                wall_zero_ok = false;
        }
        std::printf("    max |gradE + metric + wall - total| = %.3e Eh\n", worst);
        check(worst <= 1e-12, "gradE + metric + wall = dG/ds (both estimator-I forms) to 1e-12 Eh");
        check(wall_zero_ok, "passive wall: wall column exactly 0, gradE = <b.grad(E+W)>");
    }

    for (int i = 1; i < argc; ++i)
        csvCheck(argv[i]);
    std::cout << (failures ? "FAILED" : "ALL PASSED") << " (" << failures << " failures)" << std::endl;
    return failures ? 1 : 0;
}
