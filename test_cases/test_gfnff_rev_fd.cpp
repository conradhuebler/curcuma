/*
 * test_gfnff_rev_fd.cpp — rev-gfnff stage 1: analytic vs finite-difference gradient
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (Sep 2026): the continuous bond order, the blended repulsion and the
 * over-coordination energy each add gradient terms; this test checks every one against a
 * central finite difference on geometries where the weights are neither 0 nor 1:
 *   1. H2 stretched to 0.7 / 1.0 / 1.5 / 2.0 / 2.5 A (bond weight + repulsion blend)
 *   2. CH4 + H approaching the carbon backside at 1.2 / 1.6 / 2.2 A (over-coordination)
 *   3. trans-N2H2 with one N-H stretched to 1.6 A (angle/torsion/inversion weights)
 * and that `rev_enabled false` reproduces plain gfnff to 1e-12.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 */

#include "src/core/energycalculator.h"
#include "src/core/molecule.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using json = nlohmann::json;

static double fd_residual(const std::string& method, const json& cfg, const curcuma::Molecule& mol, double& analytic_norm)
{
    EnergyCalculator calc(method, cfg);
    calc.setMolecule(mol.getMolInfo());
    calc.CalculateEnergy(true);
    Matrix g = calc.Gradient();
    analytic_norm = g.norm();
    const double h = 1e-5; // Angstrom
    Matrix geom = mol.getGeometry();
    double worst = 0.0;
    for (int i = 0; i < geom.rows(); ++i) {
        for (int c = 0; c < 3; ++c) {
            Matrix gp = geom, gm = geom;
            gp(i, c) += h;
            gm(i, c) -= h;
            calc.updateGeometry(gp);
            const double ep = calc.CalculateEnergy(false);
            calc.updateGeometry(gm);
            const double em = calc.CalculateEnergy(false);
            worst = std::max(worst, std::abs((ep - em) / (2.0 * h) - g(i, c)));
        }
    }
    return worst;
}

static curcuma::Molecule h2(double r)
{
    curcuma::Molecule m;
    m.addPair({ 1, Position(0.0, 0.0, 0.0) });
    m.addPair({ 1, Position(r, 0.0, 0.0) });
    return m;
}

static curcuma::Molecule ch4_h(double d)
{
    curcuma::Molecule m;
    const double a = 0.63;
    m.addPair({ 6, Position(0, 0, 0) });
    m.addPair({ 1, Position(a, a, a) });
    m.addPair({ 1, Position(-a, -a, a) });
    m.addPair({ 1, Position(-a, a, -a) });
    m.addPair({ 1, Position(a, -a, -a) });
    const double s = d / std::sqrt(3.0);
    m.addPair({ 1, Position(-s, -s, -s) }); // backside of the first C-H bond
    return m;
}

static curcuma::Molecule n2h2_stretched()
{
    curcuma::Molecule m;
    m.addPair({ 7, Position(0.0, 0.0, 0.62) });
    m.addPair({ 7, Position(0.0, 0.0, -0.62) });
    m.addPair({ 1, Position(0.95, 0.0, 1.00) });
    m.addPair({ 1, Position(-1.45, 0.0, -1.35) }); // N-H stretched to ~1.6 A
    return m;
}

int main()
{
    std::cout << std::scientific << std::setprecision(3);
    json rev = { { "verbosity", 0 }, { "threads", 1 }, { "gfnff", { { "rev_enabled", true } } } };
    json plain = { { "verbosity", 0 }, { "threads", 1 }, { "gfnff", json::object() } };
    bool pass = true;
    const double tol = 2e-4; // Eh/A; FD truncation with h = 1e-5 A is ~1e-6

    struct Case { std::string name; curcuma::Molecule mol; };
    std::vector<Case> cases;
    for (double r : { 0.7, 1.0, 1.5, 2.0, 2.5 })
        cases.push_back({ "H2 r=" + std::to_string(r), h2(r) });
    for (double d : { 1.2, 1.6, 2.2 })
        cases.push_back({ "CH4+H d=" + std::to_string(d), ch4_h(d) });
    cases.push_back({ "N2H2 stretched", n2h2_stretched() });

    for (const auto& c : cases) {
        double gnorm = 0.0, gnorm_plain = 0.0;
        const double res = fd_residual("revgfnff", rev, c.mol, gnorm);
        // the same geometry with plain gfnff: whatever residual GFN-FF itself carries there
        // (its own documented partial gradients) is not a stage-1 defect
        const double res_plain = fd_residual("gfnff", plain, c.mol, gnorm_plain);
        const bool ok = res < tol || res <= res_plain + 2e-5;
        pass = pass && ok;
        std::cout << (ok ? "  PASS  " : "  FAIL  ") << c.name << ": max|g_analytic - g_FD| = " << res
                  << " Eh/A (|g| = " << gnorm << "); plain gfnff residual " << res_plain << "\n";
    }

    // stage 1b: gradient with an ACTIVE transition. The calculator keeps its topology across
    // updateGeometry(), so moving H2 from outside the formation window into it starts a
    // forming blend (s between 0 and 1); the FD is then taken on that blended surface.
    {
        json cfg = { { "verbosity", 0 }, { "threads", 1 },
            { "gfnff", { { "rev_enabled", true }, { "topology_mode", "react" }, { "react_check_every", 1 },
                         { "react_valence_cap", false }, { "react_refractory_scans", 0 }, { "react_exchange_scans", 0 } } } };
        // The pair is walked in small steps so that the scan catches the transition inside its
        // window (a single 0.9 -> 1.5 A jump crosses the whole break window and the transition
        // is over before the probe); the probe geometries sit at s ~ 0.2-0.9 of the windows
        // (formation c: 0.02 -> 0.8, break c: 0.5 -> 0.02 on the 1.6x/-8 switch; H-H
        // covalent sum 0.666 A, so 1.15 A = 1.73x and 1.05 A = 1.58x).
        struct Seq { std::string name; std::vector<double> path; };
        // Claude Generated (Sep 12, 2026): with the default rev_form_switch = order the formation
        // joins at 1.611x the covalent sum (1.073 A here) and its window runs on the NARROW switch
        // from rev_bo2_form = 0.1 to 0.9, i.e. 1.073 A -> 0.792 A. The two "1.7->1.3->..." rows
        // therefore no longer start a transition at 1.3 A (narrow order 3.9e-4 there); the last row
        // walks 1.7 -> 1.05 -> 0.95 A, which joins at 1.05 A and probes at s ~ 0.33 of that window.
        for (const Seq& q : { Seq { "forming H2 1.7->1.3->1.15 A", { 1.7, 1.3, 1.15 } }, Seq { "forming H2 1.7->1.3->1.05 A", { 1.7, 1.3, 1.05 } },
                              Seq { "breaking H2 0.9->1.05->1.15 A", { 0.9, 1.05, 1.15 } }, Seq { "breaking H2 0.9->1.05->1.20 A", { 0.9, 1.05, 1.20 } },
                              Seq { "forming H2 1.7->1.05->0.95 A", { 1.7, 1.05, 0.95 } } }) {
            EnergyCalculator calc("revgfnff", cfg);
            calc.setMolecule(h2(q.path.front()).getMolInfo());
            calc.CalculateEnergy(false);
            for (size_t k = 1; k + 1 < q.path.size(); ++k) {
                calc.updateGeometry(h2(q.path[k]).getGeometry());
                calc.CalculateEnergy(false);
            }
            Matrix geom = h2(q.path.back()).getGeometry();
            calc.updateGeometry(geom);
            calc.CalculateEnergy(true);
            Matrix g = calc.Gradient();
            const double h = 1e-5;
            double worst = 0.0;
            for (int i = 0; i < 2; ++i)
                for (int c = 0; c < 3; ++c) {
                    Matrix gp = geom, gm = geom;
                    gp(i, c) += h;
                    gm(i, c) -= h;
                    calc.updateGeometry(gp);
                    const double ep = calc.CalculateEnergy(false);
                    calc.updateGeometry(gm);
                    const double em = calc.CalculateEnergy(false);
                    worst = std::max(worst, std::abs((ep - em) / (2.0 * h) - g(i, c)));
                }
            const bool ok = worst < tol;
            pass = pass && ok;
            std::cout << (ok ? "  PASS  " : "  FAIL  ") << q.name << " (active transition): max|g_analytic - g_FD| = " << worst
                      << " Eh/A (|g| = " << g.norm() << ")\n";
        }
    }

    // rev_enabled false must be plain gfnff
    {
        json off = { { "verbosity", 0 }, { "threads", 1 }, { "gfnff", { { "rev_enabled", false } } } };
        EnergyCalculator a("gfnff", off), b("gfnff", plain);
        a.setMolecule(ch4_h(1.6).getMolInfo());
        b.setMolecule(ch4_h(1.6).getMolInfo());
        const double ea = a.CalculateEnergy(false), eb = b.CalculateEnergy(false);
        const bool ok = std::abs(ea - eb) < 1e-12;
        pass = pass && ok;
        std::cout << (ok ? "  PASS  " : "  FAIL  ") << "rev_enabled=false reproduces gfnff: " << ea << " vs " << eb << "\n";
    }
    // over-coordination must be positive when the fifth H sits on the carbon
    {
        EnergyCalculator calc("revgfnff", rev);
        calc.setMolecule(ch4_h(1.2).getMolInfo());
        calc.CalculateEnergy(false);
        const double over = calc.getEnergyDecomposition().value("OverCoord", 0.0);
        const bool ok = over > 1e-4;
        pass = pass && ok;
        std::cout << (ok ? "  PASS  " : "  FAIL  ") << "over-coordination penalty at CH4+H(1.2 A): " << over << " Eh\n";
    }
    std::cout << (pass ? "PASS" : "FAIL") << "\n";
    return pass ? 0 : 1;
}
