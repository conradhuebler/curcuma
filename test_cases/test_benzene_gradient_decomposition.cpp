/**
 * Benzene Gradient Term-by-Term Decomposition
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Diagnostic tool to identify which GFN-FF gradient term causes sign flip in benzene
 *
 * Strategy:
 * 1. Test each energy term in isolation
 * 2. Identify which term(s) have wrong sign
 * 3. Compare with Fortran reference implementation
 *
 * Claude Generated - February 2026
 */

#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "src/core/molecule.h"
#include "src/core/global.h"
#include "src/core/curcuma_logger.h"
#include "src/tools/formats.h"
#include "src/core/units.h"

#include "src/core/energy_calculators/ff_methods/forcefield.h"
#include "src/core/energy_calculators/ff_methods/gfnff.h"

using namespace curcuma;

const double FINITE_DIFF_STEP = 1e-4;  // Bohr

// Benzene geometry (from validation directory)
std::vector<int> benzene_atoms = {6,6,6,6,6,6,1,1,1,1,1,1};  // 6 C + 6 H
Matrix benzene_geometry_angstrom(12, 3);

void initializeBenzeneGeometry() {
    benzene_geometry_angstrom <<
        1.398420,    0.000000,    0.000000,
        0.699210,    1.211530,    0.000000,
       -0.699210,    1.211530,    0.000000,
       -1.398420,    0.000000,    0.000000,
       -0.699210,   -1.211530,    0.000000,
        0.699210,   -1.211530,    0.000000,
        2.487520,    0.000000,    0.000000,
        1.243760,    2.154210,    0.000000,
       -1.243760,    2.154210,    0.000000,
       -2.487520,    0.000000,    0.000000,
       -1.243760,   -2.154210,    0.000000,
        1.243760,   -2.154210,    0.000000;
}

Matrix calculateNumericalGradient(GFNFF& gfnff, const Matrix& geometry_bohr) {
    const int natoms = geometry_bohr.rows();
    Matrix numerical_gradient = Matrix::Zero(natoms, 3);

    Matrix geom_plus = geometry_bohr;
    Matrix geom_minus = geometry_bohr;

    for (int i = 0; i < natoms; ++i) {
        for (int coord = 0; coord < 3; ++coord) {
            geom_plus(i, coord) = geometry_bohr(i, coord) + FINITE_DIFF_STEP;
            gfnff.UpdateMolecule(geom_plus * CurcumaUnit::Length::BOHR_TO_ANGSTROM);
            double E_plus = gfnff.Calculation(false);

            geom_minus(i, coord) = geometry_bohr(i, coord) - FINITE_DIFF_STEP;
            gfnff.UpdateMolecule(geom_minus * CurcumaUnit::Length::BOHR_TO_ANGSTROM);
            double E_minus = gfnff.Calculation(false);

            numerical_gradient(i, coord) = (E_plus - E_minus) / (2.0 * FINITE_DIFF_STEP);

            geom_plus(i, coord) = geometry_bohr(i, coord);
            geom_minus(i, coord) = geometry_bohr(i, coord);
        }
    }

    // Convert to Hartree/Angstrom to match analytical gradient
    return numerical_gradient / CurcumaUnit::Length::BOHR_TO_ANGSTROM;
}

struct TermTestResult {
    std::string term_name;
    double energy;
    Matrix analytical_grad;
    Matrix numerical_grad;
    double max_error;
    double sign_ratio;  // analytical/numerical for first component
};

TermTestResult testGradientTerm(const std::string& term_name, const json& gfnff_params) {
    TermTestResult result;
    result.term_name = term_name;

    // Create molecule
    Molecule mol(benzene_atoms.size());
    for (size_t i = 0; i < benzene_atoms.size(); ++i) {
        mol.setAtom(std::to_string(benzene_atoms[i]), (int)i);
    }
    mol.setGeometry(benzene_geometry_angstrom);

    // Create GFN-FF instance
    GFNFF gfnff(gfnff_params);
    gfnff.InitialiseMolecule(mol.getMolInfo());

    // Get geometry in Bohr for numerical gradient
    Matrix geometry_bohr = benzene_geometry_angstrom * CurcumaUnit::Length::ANGSTROM_TO_BOHR;

    // Calculate analytical gradient
    gfnff.UpdateMolecule(benzene_geometry_angstrom);
    result.energy = gfnff.Calculation(true);
    result.analytical_grad = gfnff.Gradient();

    // Calculate numerical gradient
    result.numerical_grad = calculateNumericalGradient(gfnff, geometry_bohr);

    // Calculate errors
    Matrix diff = result.analytical_grad - result.numerical_grad;
    result.max_error = diff.array().abs().maxCoeff();

    // Calculate sign ratio for first atom X component (to detect sign flips)
    if (std::abs(result.numerical_grad(0,0)) > 1e-10) {
        result.sign_ratio = result.analytical_grad(0,0) / result.numerical_grad(0,0);
    } else {
        result.sign_ratio = 0.0;
    }

    return result;
}

void printTermResult(const TermTestResult& result) {
    std::cout << "\n========================================" << std::endl;
    std::cout << "Term: " << result.term_name << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "  Energy: " << std::scientific << std::setprecision(8) << result.energy << " Eh" << std::endl;
    std::cout << "  Max gradient error: " << result.max_error << " Eh/Bohr" << std::endl;
    std::cout << "  Sign ratio (atom 0, x): " << std::fixed << std::setprecision(4) << result.sign_ratio << std::endl;

    std::cout << "\n  First atom gradient:" << std::endl;
    std::cout << "    Analytical: [" << std::scientific
              << result.analytical_grad(0,0) << ", "
              << result.analytical_grad(0,1) << ", "
              << result.analytical_grad(0,2) << "]" << std::endl;
    std::cout << "    Numerical:  ["
              << result.numerical_grad(0,0) << ", "
              << result.numerical_grad(0,1) << ", "
              << result.numerical_grad(0,2) << "]" << std::endl;

    // Interpret sign ratio
    if (std::abs(result.sign_ratio + 1.0) < 0.3) {  // Close to -1
        std::cout << "\n  ⚠️  **SIGN FLIP DETECTED** - Analytical gradient has opposite sign!" << std::endl;
    } else if (std::abs(result.sign_ratio - 1.0) < 0.3) {  // Close to +1
        std::cout << "\n  ✓  Sign appears correct" << std::endl;
    } else if (std::abs(result.sign_ratio) < 0.1) {  // Close to 0
        std::cout << "\n  ℹ️  Gradient too small to determine sign" << std::endl;
    } else {
        std::cout << "\n  ⚠️  Unexpected sign ratio - needs investigation" << std::endl;
    }
}

int main() {
    std::cout << "\n";
    std::cout << "╔════════════════════════════════════════════════════════════╗\n";
    std::cout << "║   Benzene Gradient Term-by-Term Decomposition              ║\n";
    std::cout << "║   Diagnostic tool to identify sign flip source             ║\n";
    std::cout << "╚════════════════════════════════════════════════════════════╝\n";
    std::cout << "\n";

    CurcumaLogger::set_verbosity(0);  // Suppress method output

    initializeBenzeneGeometry();

    // Test 1: Bonds only
    json params_bonds = {
        {"method", "cgfnff"},
        {"verbosity", 0},
        {"threads", 1},
        {"gfnff", {
            {"enable_angles", false},
            {"enable_torsions", false},
            {"enable_repulsion", false},
            {"enable_dispersion", false},
            {"enable_coulomb", false},
            {"enable_hbond", false}
        }}
    };
    auto result_bonds = testGradientTerm("Bonds only", params_bonds);
    printTermResult(result_bonds);

    // Test 2: Angles only
    json params_angles = {
        {"method", "cgfnff"},
        {"verbosity", 0},
        {"threads", 1},
        {"gfnff", {
            {"enable_bonds", false},
            {"enable_torsions", false},
            {"enable_repulsion", false},
            {"enable_dispersion", false},
            {"enable_coulomb", false},
            {"enable_hbond", false}
        }}
    };
    auto result_angles = testGradientTerm("Angles only", params_angles);
    printTermResult(result_angles);

    // Test 3: Coulomb only
    json params_coulomb = {
        {"method", "cgfnff"},
        {"verbosity", 0},
        {"threads", 1},
        {"gfnff", {
            {"enable_bonds", false},
            {"enable_angles", false},
            {"enable_torsions", false},
            {"enable_repulsion", false},
            {"enable_dispersion", false},
            {"enable_hbond", false}
        }}
    };
    auto result_coulomb = testGradientTerm("Coulomb only", params_coulomb);
    printTermResult(result_coulomb);

    // Test 4: Dispersion only
    json params_disp = {
        {"method", "cgfnff"},
        {"verbosity", 0},
        {"threads", 1},
        {"gfnff", {
            {"enable_bonds", false},
            {"enable_angles", false},
            {"enable_torsions", false},
            {"enable_repulsion", false},
            {"enable_coulomb", false},
            {"enable_hbond", false}
        }}
    };
    auto result_disp = testGradientTerm("Dispersion only", params_disp);
    printTermResult(result_disp);

    // Test 5: Repulsion only
    json params_rep = {
        {"method", "cgfnff"},
        {"verbosity", 0},
        {"threads", 1},
        {"gfnff", {
            {"enable_bonds", false},
            {"enable_angles", false},
            {"enable_torsions", false},
            {"enable_dispersion", false},
            {"enable_coulomb", false},
            {"enable_hbond", false}
        }}
    };
    auto result_rep = testGradientTerm("Repulsion only", params_rep);
    printTermResult(result_rep);

    // Test 6: Full gradient
    json params_full = {
        {"method", "cgfnff"},
        {"verbosity", 0},
        {"threads", 1}
    };
    auto result_full = testGradientTerm("Full GFN-FF", params_full);
    printTermResult(result_full);

    // Summary
    std::cout << "\n";
    std::cout << "═══════════════════════════════════════════════════════════════" << std::endl;
    std::cout << "SUMMARY: Sign Ratio Analysis" << std::endl;
    std::cout << "═══════════════════════════════════════════════════════════════" << std::endl;
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "  Bonds:      " << result_bonds.sign_ratio << std::endl;
    std::cout << "  Angles:     " << result_angles.sign_ratio << std::endl;
    std::cout << "  Coulomb:    " << result_coulomb.sign_ratio << std::endl;
    std::cout << "  Dispersion: " << result_disp.sign_ratio << std::endl;
    std::cout << "  Repulsion:  " << result_rep.sign_ratio << std::endl;
    std::cout << "  Full:       " << result_full.sign_ratio << std::endl;
    std::cout << "\nInterpretation:" << std::endl;
    std::cout << "  Ratio ≈ +1.0: Correct sign" << std::endl;
    std::cout << "  Ratio ≈ -1.0: SIGN FLIP (needs fixing)" << std::endl;
    std::cout << "  Ratio ≈  0.0: Gradient too small" << std::endl;

    return 0;
}
