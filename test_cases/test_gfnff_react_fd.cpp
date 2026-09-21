/*
 * test_gfnff_react_fd.cpp — react-topology rebuild correctness
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (Aug 2026): validates the react-mode bond-formation rebuild.
 *
 * 1. Equivalence: after a react-mode formation rebuild, energy and analytic
 *    gradient must equal those of a FRESH initialisation with the same bond
 *    topology (forced bonds) at the same geometry — the rebuild must reproduce
 *    exactly the state a from-scratch setup would produce.
 * 2. FD consistency: the analytic-vs-finite-difference residual on the
 *    post-rebuild surface must match the residual of the fresh-init reference
 *    surface. (The residual itself is a pre-existing property of the GFN-FF
 *    gradient for strongly CN-dependent bonds like H-H and is reported, not
 *    gated here.)
 *
 * Geometry: two H-H pairs at 1.00 A — outside plain 1.3x geometric detection
 * (0.87 A for H-H), inside the react formation radius (1.6 x threshold = 1.06 A),
 * so react mode starts with 0 bonds and forms both on the first energy call.
 */
#include <cmath>
#include <iomanip>
#include <iostream>

#include "src/core/energycalculator.h"
#include "src/core/molecule.h"

#include "json.hpp"
using json = nlohmann::json;

static double fd_check(EnergyCalculator& calc, const Matrix& geom, Matrix& g_analytic)
{
    const double h = 1e-5; // Angstrom
    Matrix g_fd = Matrix::Zero(geom.rows(), 3);
    for (int i = 0; i < geom.rows(); ++i) {
        for (int c = 0; c < 3; ++c) {
            Matrix gp = geom, gm = geom;
            gp(i, c) += h;
            gm(i, c) -= h;
            calc.updateGeometry(gp);
            const double ep = calc.CalculateEnergy(false);
            calc.updateGeometry(gm);
            const double em = calc.CalculateEnergy(false);
            g_fd(i, c) = (ep - em) / (2.0 * h); // Eh/Angstrom
        }
    }
    calc.updateGeometry(geom);
    calc.CalculateEnergy(true);
    g_analytic = calc.Gradient();
    return (g_analytic - g_fd).cwiseAbs().maxCoeff();
}

int main()
{
    curcuma::Molecule mol;
    mol.addPair({ 1, Position(0.00, 0.0, 0.0) });
    mol.addPair({ 1, Position(1.00, 0.0, 0.0) });
    mol.addPair({ 1, Position(0.00, 3.5, 0.0) });
    mol.addPair({ 1, Position(1.00, 3.5, 0.0) });
    const Matrix geom = mol.getGeometry();

    std::cout << std::scientific << std::setprecision(4);

    // --- React path: bonds form via the hysteresis scan on the first call ---
    json cfg_react = { { "verbosity", 0 }, { "threads", 1 },
        { "gfnff", { { "topology_mode", "react" }, { "react_check_every", 1 } } } };
    EnergyCalculator calc_react("gfnff", cfg_react);
    calc_react.setMolecule(mol.getMolInfo());
    calc_react.CalculateEnergy(false); // triggers formation rebuild
    const double e_react = calc_react.CalculateEnergy(true);
    Matrix g_react = calc_react.Gradient();

    // --- Reference path: fresh initialisation with the same topology (forced bonds) ---
    Mol forced = mol.getMolInfo();
    forced.m_bonds = { { 0, 1 }, { 2, 3 } };
    json cfg_auto = { { "verbosity", 0 }, { "threads", 1 }, { "gfnff", json::object() } };
    EnergyCalculator calc_ref("gfnff", cfg_auto);
    calc_ref.setMolecule(forced);
    const double e_ref = calc_ref.CalculateEnergy(true);
    Matrix g_ref = calc_ref.Gradient();

    const double de = std::abs(e_react - e_ref);
    const double dg = (g_react - g_ref).cwiseAbs().maxCoeff();
    std::cout << "rebuild vs fresh-init:  |dE| = " << de << " Eh, max|dg| = " << dg << " Eh/A\n";

    // --- FD residual on both surfaces (informational + consistency gate) ---
    Matrix ga, gb;
    const double fd_react = fd_check(calc_react, geom, ga);
    EnergyCalculator calc_ref2("gfnff", cfg_auto);
    calc_ref2.setMolecule(forced);
    calc_ref2.CalculateEnergy(false);
    const double fd_ref = fd_check(calc_ref2, geom, gb);
    std::cout << "FD residual: react surface " << fd_react
              << " Eh/A, fresh-init surface " << fd_ref << " Eh/A\n";

    // A bound pair at 1.00 A must be attractive along the bond axis.
    const bool bonded_surface = (g_react(0, 0) < -0.01);
    std::cout << "bond formed (attractive surface): " << (bonded_surface ? "yes" : "NO") << "\n";

    const bool pass = (de < 1e-9) && (dg < 1e-9)
        && (std::abs(fd_react - fd_ref) < 1e-6) && bonded_surface;
    std::cout << (pass ? "PASS" : "FAIL") << "\n";
    return pass ? 0 : 1;
}
