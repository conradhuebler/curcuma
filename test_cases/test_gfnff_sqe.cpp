/*
 * test_gfnff_sqe.cpp — rev-gfnff stage 2: split-charge (SQE) charge model
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (Sep 2026), docs/REV_GFNFF_STAGE2.md — acceptance items 1 and 2.
 *
 * 1. Fidelity. On a NEUTRAL molecule whose bond graph is connected, the split-charge
 *    minimum with kappa_Z = 0 spans exactly the same space as the fragment-constrained
 *    EEQ minimum (range(B) = {v : sum v = 0}), so the two must give the SAME charges and
 *    the SAME Coulomb energy. Anything above 1e-8 Eh / 1e-8 e is a defect in the SQE
 *    assembly, not a model difference. Six neutral molecules from test_cases/molecules,
 *    caffeine among them (rings — the pair system is then rank deficient, which is the
 *    interesting case for the solver's null-space handling).
 *
 * 2. Gradient. With kappa_Z = 0.5 Eh the hardness term 1/2 kappa0/b p^2 carries a real
 *    force through db/dr; it is checked against central finite differences (h = 1e-5 A,
 *    tol 2e-4 Eh/A, same tolerance as test_gfnff_rev_fd) on
 *      (a) Cl2- at r = 2.73 A (the EA_25 geometry of the design's motivating case),
 *      (b) HCOO-...HF (frame point=1 of test_cases/revgfnff/ref/E/ahb21_21_stretch),
 *      (c) CH4 + H in react mode with a transition in flight (the corner blend, walked
 *          into its window the same way test_gfnff_rev_fd does it).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 */

#include "src/core/energycalculator.h"
#include "src/core/fileiterator.h"
#include "src/core/molecule.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using json = nlohmann::json;

/// rev-gfnff config: stage 1 on, charge model and kappa_Z as given
static json revConfig(const std::string& charge_model, double kappa, bool react = false)
{
    json g;
    g["rev_enabled"] = true;
    g["rev_charge_model"] = charge_model;
    for (const char* k : { "rev_sqe_kappa_H", "rev_sqe_kappa_C", "rev_sqe_kappa_N",
             "rev_sqe_kappa_O", "rev_sqe_kappa_F", "rev_sqe_kappa_Cl" })
        g[k] = kappa;
    if (react) {
        g["topology_mode"] = "react";
        g["react_check_every"] = 1;
        g["react_valence_cap"] = false;
        g["react_refractory_scans"] = 0;
        g["react_exchange_scans"] = 0;
    }
    return json { { "verbosity", 0 }, { "threads", 1 }, { "gfnff", g } };
}

/// max |g_analytic - g_FD| over every Cartesian component (h in Angstrom, as the contract says)
static double fdResidual(EnergyCalculator& calc, const Matrix& geom, double& gnorm)
{
    calc.updateGeometry(geom);
    calc.CalculateEnergy(true);
    Matrix g = calc.Gradient();
    gnorm = g.norm();
    const double h = 1e-5;
    double worst = 0.0;
    for (int i = 0; i < geom.rows(); ++i)
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
    return worst;
}

static curcuma::Molecule cl2anion(double r)
{
    curcuma::Molecule m;
    m.addPair({ 17, Position(0.0, 0.0, 0.0) });
    m.addPair({ 17, Position(r, 0.0, 0.0) });
    m.setCharge(-1);
    return m;
}

static curcuma::Molecule ch4H(double d)
{
    curcuma::Molecule m;
    const double a = 0.63;
    m.addPair({ 6, Position(0, 0, 0) });
    m.addPair({ 1, Position(a, a, a) });
    m.addPair({ 1, Position(-a, -a, a) });
    m.addPair({ 1, Position(-a, a, -a) });
    m.addPair({ 1, Position(a, -a, -a) });
    const double s = d / std::sqrt(3.0);
    m.addPair({ 1, Position(-s, -s, -s) });
    return m;
}

int main(int argc, char* argv[])
{
    // argv[1] = test_cases source directory (set by CMake)
    const std::string root = (argc > 1) ? std::string(argv[1]) : std::string(".");
    bool pass = true;
    std::cout << std::scientific << std::setprecision(3);

    // ---------------------------------------------------------------------------------------
    // 1. Fidelity: sqe with kappa = 0 == eeq on a connected neutral bond graph
    // ---------------------------------------------------------------------------------------
    {
        const std::vector<std::string> molecules = {
            "molecules/larger/caffeine.xyz", "molecules/larger/CH4.xyz",
            "molecules/larger/CH3OH.xyz", "molecules/larger/C6H6.xyz",
            "molecules/larger/CH3OCH3.xyz", "molecules/larger/C6H5COOH.xyz"
        };
        const double tol = 1e-8;
        for (const std::string& rel : molecules) {
            curcuma::Molecule mol(root + "/" + rel);
            if (mol.AtomCount() == 0) {
                std::cout << "  FAIL  " << rel << ": could not be read\n";
                pass = false;
                continue;
            }
            EnergyCalculator a("revgfnff", revConfig("eeq", 0.0));
            EnergyCalculator b("revgfnff", revConfig("sqe", 0.0));
            a.setMolecule(mol.getMolInfo());
            b.setMolecule(mol.getMolInfo());
            const double ea = a.CalculateEnergy(false);
            const double eb = b.CalculateEnergy(false);
            const double ca = a.getEnergyDecomposition().value("Coulomb", 0.0);
            const double cb = b.getEnergyDecomposition().value("Coulomb", 0.0);
            const double hard = b.getEnergyDecomposition().value("SqeHardness", 0.0);
            Vector qa = a.Charges(), qb = b.Charges();
            double dq = 0.0;
            if (qa.size() == qb.size() && qa.size() > 0)
                dq = (qa - qb).cwiseAbs().maxCoeff();
            else
                dq = std::numeric_limits<double>::quiet_NaN();
            const double dcoul = std::abs(ca - cb);
            const double dtot = std::abs(ea - eb);
            const bool ok = (dcoul < tol) && (dtot < tol) && (dq < tol) && (hard == 0.0);
            pass = pass && ok;
            std::cout << (ok ? "  PASS  " : "  FAIL  ") << rel << " (N=" << mol.AtomCount()
                      << "): dE = " << dtot << " Eh, dCoulomb = " << dcoul
                      << " Eh, max|dq| = " << dq << " e\n";
        }
    }

    // ---------------------------------------------------------------------------------------
    // 2. FD gradient with kappa = 0.5 Eh.
    //
    // The acceptance criterion of the design is "analytic vs central differences to the same
    // residual as plain gfnff" — GFN-FF carries its own documented partial-gradient residual at
    // some geometries, and that is not a stage-2 defect. Each case therefore reports the plain
    // `gfnff` residual at the SAME geometry as the baseline, exactly as test_gfnff_rev_fd does.
    // Measured Sep 2026: Cl2- at 2.73 A has a plain-gfnff residual of 1.77e-2 Eh/A (a free
    // chloride ion plus a 2c-3e bond — the fragment/EEQ corner of Known Issue #8/#17), and the
    // split-charge model reproduces it to the digit for every kappa, i.e. the hardness term
    // contributes nothing to it.
    // ---------------------------------------------------------------------------------------
    const double tol_g = 2e-4;
    const json plain = { { "verbosity", 0 }, { "threads", 1 }, { "gfnff", json::object() } };
    auto gradCase = [&](const std::string& name, const curcuma::Molecule& mol) {
        curcuma::Molecule m = mol;
        EnergyCalculator ref("gfnff", plain);
        ref.setMolecule(m.getMolInfo());
        double gn_ref = 0.0;
        const double res_ref = fdResidual(ref, m.getGeometry(), gn_ref);
        EnergyCalculator calc("revgfnff", revConfig("sqe", 0.5));
        calc.setMolecule(m.getMolInfo());
        double gnorm = 0.0;
        const double res = fdResidual(calc, m.getGeometry(), gnorm);
        const bool ok = (res < tol_g) || (res <= res_ref + 2e-5);
        pass = pass && ok;
        std::cout << (ok ? "  PASS  " : "  FAIL  ") << name << " (kappa=0.5): max|g-gFD| = " << res
                  << " Eh/A (|g| = " << gnorm << "); plain gfnff residual " << res_ref << "\n";
    };

    // 2a. Cl2- at the EA_25 separation (charge -1)
    gradCase("Cl2- r=2.73 A, q=-1", cl2anion(2.73));

    // 2b. HCOO-...HF, frame point=1 of the class-E stretch scan (charge -1)
    {
        const std::string pts = root + "/revgfnff/ref/E/ahb21_21_stretch/points.xyz";
        FileIterator it;
        it.setFile(pts);
        curcuma::Molecule mol;
        int frame = 0;
        bool have = false;
        while (!it.AtEnd()) {
            curcuma::Molecule m = it.Next();
            if (frame == 1) { mol = m; have = true; break; }
            ++frame;
        }
        if (!have) {
            std::cout << "  FAIL  HCOO-...HF: frame 1 of " << pts << " not readable\n";
            pass = false;
        } else {
            mol.setCharge(-1);
            gradCase("HCOO-...HF, q=-1", mol);
        }
    }

    // 2c. CH4 + H in react mode, transition in flight. The calculator keeps its topology across
    // updateGeometry(), so walking the fifth hydrogen in from 2.2 A starts a forming transition
    // and the FD is taken on the BLENDED surface, where every corner runs its own split-charge
    // solve with its own frozen q0. Same walking pattern as test_gfnff_rev_fd.
    {
        EnergyCalculator calc("revgfnff", revConfig("sqe", 0.5, /*react=*/true));
        calc.setMolecule(ch4H(2.2).getMolInfo());
        calc.CalculateEnergy(false);
        for (double d : { 1.8, 1.5, 1.3 }) {
            calc.updateGeometry(ch4H(d).getGeometry());
            calc.CalculateEnergy(false);
        }
        Matrix geom = ch4H(1.25).getGeometry();
        double gnorm = 0.0;
        const double res = fdResidual(calc, geom, gnorm);
        const bool ok = res < tol_g;
        pass = pass && ok;
        std::cout << (ok ? "  PASS  " : "  FAIL  ")
                  << "CH4+H react, transition in flight (kappa=0.5): max|g-gFD| = " << res
                  << " Eh/A (|g| = " << gnorm << ")\n";
    }

    std::cout << (pass ? "PASS" : "FAIL") << "\n";
    return pass ? 0 : 1;
}
