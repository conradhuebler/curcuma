/*
 * < GFN1/GFN2 Gradient Finite-Difference Validation >
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Validates the analytical gradients of the native GFN1-xTB and GFN2-xTB
 * implementations against central-finite-difference numerical gradients.
 * Tests AP4 (repulsion, H0/Pulay, Coulomb, CN chain-rule) and
 * AP5 (GFN2 multipole direct gradient, Schritt 1).
 *
 * Usage:
 *   test_xtb_gradient <molecule.xyz> [...]
 *
 *   Returns 0 if all molecules pass within tolerance, 1 otherwise.
 *
 * Tolerance:  5e-4 Eh/Å (roughly 1e-3 Eh/Bohr) per gradient component.
 *             Both methods should pass; GFN2 tests the AP5 multipole term.
 *
 * Claude Generated (AP5, Apr 2026). GPL-3.0.
 */

#include "src/core/energycalculator.h"
#include "src/core/molecule.h"
#include "src/core/curcuma_logger.h"
#include "src/core/global.h"

#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

// Central-finite-difference step in Angstroms
static constexpr double FD_H = 1.0e-3;

// Pass threshold: max |G_anal - G_num| per component in Eh/Å
static constexpr double GRAD_TOL = 5.0e-4;

struct GradResult {
    std::string mol_name;
    std::string method;
    int         nat;
    double      max_err;
    double      rms_err;
    bool        passed;
    double      wall_s;
};

static Mol buildMol(const std::string& filename)
{
    curcuma::Molecule m(filename);
    return m.getMolInfo();
}

/*
 * Compute the analytical gradient and then validate it against the
 * 3N central-finite-difference numerical gradient.
 */
static GradResult validateGradient(const std::string& xyz_path,
                                   const std::string& mol_name,
                                   const std::string& method)
{
    auto t0 = std::chrono::steady_clock::now();

    GradResult res;
    res.mol_name = mol_name;
    res.method   = method;
    res.passed   = false;
    res.max_err  = 0.0;
    res.rms_err  = 0.0;

    Mol mol = buildMol(xyz_path);
    res.nat = mol.m_number_atoms;

    const json cfg = {{"verbosity", 0}, {"threads", 1}};

    // ─── Analytical gradient ─────────────────────────────────────────────────
    EnergyCalculator calc(method, cfg);
    calc.setMolecule(mol);
    calc.CalculateEnergy(/*gradient=*/true);
    const Matrix G_anal = calc.Gradient();   // nat×3, units: same as FD below

    // ─── Numerical gradient (central differences in Å) ───────────────────────
    Geometry geo_ref = mol.m_geometry;       // nat×3, rows=atoms, cols=xyz in Å
    Matrix G_num = Matrix::Zero(res.nat, 3);

    for (int at = 0; at < res.nat; ++at) {
        for (int c = 0; c < 3; ++c) {
            mol.m_geometry = geo_ref;
            mol.m_geometry(at, c) += FD_H;
            EnergyCalculator cp(method, cfg);
            cp.setMolecule(mol);
            const double Ep = cp.CalculateEnergy(false);

            mol.m_geometry = geo_ref;
            mol.m_geometry(at, c) -= FD_H;
            EnergyCalculator cm(method, cfg);
            cm.setMolecule(mol);
            const double Em = cm.CalculateEnergy(false);

            G_num(at, c) = (Ep - Em) / (2.0 * FD_H);
        }
    }
    mol.m_geometry = geo_ref;

    // ─── Compare ──────────────────────────────────────────────────────────────
    double sum_sq = 0.0;
    for (int at = 0; at < res.nat; ++at) {
        for (int c = 0; c < 3; ++c) {
            const double d = std::abs(G_anal(at, c) - G_num(at, c));
            if (d > res.max_err) res.max_err = d;
            sum_sq += d * d;
        }
    }
    res.rms_err = std::sqrt(sum_sq / (res.nat * 3));
    res.passed  = res.max_err < GRAD_TOL;

    auto t1 = std::chrono::steady_clock::now();
    res.wall_s = std::chrono::duration<double>(t1 - t0).count();
    return res;
}

static void printHeader()
{
    std::cout << std::string(72, '=') << "\n"
              << "  GFN1/GFN2 Gradient Finite-Difference Validation (AP4+AP5)\n"
              << "  FD step: " << FD_H << " Å   tolerance: " << GRAD_TOL << " Eh/Å\n"
              << std::string(72, '=') << "\n";
    std::cout << std::left
              << std::setw(10) << "Molecule"
              << std::setw(8)  << "Method"
              << std::setw(8)  << "Atoms"
              << std::setw(14) << "MaxErr(Eh/Å)"
              << std::setw(14) << "RMSErr(Eh/Å)"
              << std::setw(8)  << "Time(s)"
              << "Result\n"
              << std::string(72, '-') << "\n";
}

static void printResult(const GradResult& r)
{
    std::cout << std::left << std::fixed << std::setprecision(3)
              << std::setw(10) << r.mol_name
              << std::setw(8)  << r.method
              << std::setw(8)  << r.nat
              << std::setw(14) << std::setprecision(2) << std::scientific << r.max_err
              << std::setw(14) << r.rms_err
              << std::setw(8)  << std::fixed << std::setprecision(1) << r.wall_s
              << (r.passed ? "[PASS]" : "[FAIL]") << "\n";
}

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::cerr << "Usage: test_xtb_gradient <mol1.xyz> [mol2.xyz ...]\n";
        return 1;
    }

    CurcumaLogger::set_verbosity(0);
    printHeader();

    bool all_pass = true;
    for (int i = 1; i < argc; ++i) {
        const std::string path = argv[i];
        // Extract bare filename as label
        const std::size_t sl = path.rfind('/');
        std::string name = (sl == std::string::npos) ? path : path.substr(sl + 1);
        // Strip extension
        const std::size_t dot = name.rfind('.');
        if (dot != std::string::npos) name = name.substr(0, dot);

        for (const std::string& meth : {"gfn1", "gfn2"}) {
            GradResult r = validateGradient(path, name, meth);
            printResult(r);
            all_pass = all_pass && r.passed;
        }
    }

    std::cout << std::string(72, '=') << "\n"
              << (all_pass ? "All tests PASSED.\n" : "Some tests FAILED.\n");
    return all_pass ? 0 : 1;
}
