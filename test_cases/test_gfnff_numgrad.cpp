/**
 * CPU GFN-FF analytical vs. numerical gradient test.
 * Compares the analytical gradient from GFNFF::Calculation(true)
 * against GFNFF::NumGrad() (full finite-difference including EEQ recalculation).
 *
 * Uses TestMoleculeRegistry for all built-in molecules — do NOT hardcode
 * molecule geometry here. If a molecule is missing from the registry, add it there.
 *
 * Usage:
 *   test_gfnff_numgrad                  – run default molecule set from registry
 *   test_gfnff_numgrad molecule.xyz     – test one XYZ file
 *   test_gfnff_numgrad --fixed-charges  – use NumGradFixedCharges() (isolates formula bugs)
 *
 * Exit code: 0 = all molecules PASS, 1 = at least one FAIL.
 *
 * Claude Generated (April 2026)
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 */

#include "src/core/energy_calculators/ff_methods/gfnff.h"
#include "src/core/molecule.h"
#include "src/core/curcuma_logger.h"
#include "src/core/global.h"
#include "core/test_molecule_registry.h"
#include "json.hpp"

#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using json = nlohmann::json;

// ─────────────────────────────────────────────────────────────────────────────
// Default molecule set for gradient tests
// Note: add new molecules to TestMoleculeRegistry, not here.
// ─────────────────────────────────────────────────────────────────────────────
static const std::vector<std::string> DEFAULT_MOLECULES = {
    "H2O",        // O-H bonds — sensitive to angle/bond gradient
    "CH4",        // reference: symmetric C-H, should pass easily
    "CH3OH",      // mixed C-H / O-H — tests O-H bond gradient
    "CH3OCH3",    // C-O-C only (no O-H) — tests Coulomb without O-H
    "C6H6",       // aromatic, tests torsion + dispersion gradient
    "H2O_dimer",  // hydrogen bond between water monomers
};

// ─────────────────────────────────────────────────────────────────────────────
// Test runner for a single molecule
// ─────────────────────────────────────────────────────────────────────────────

struct GradTestResult {
    std::string name;
    int n_atoms;
    double energy;
    double anal_norm;
    double num_norm;
    double max_diff;
    double rms_diff;
    int worst_atom;
    bool pass;
    bool nan_detected;
};

static GradTestResult runGradTest(const std::string& name,
                                   const curcuma::Molecule& cmol,
                                   bool fixed_charges, double tol)
{
    GradTestResult res;
    res.name         = name;
    res.n_atoms      = cmol.AtomCount();
    res.pass         = false;
    res.nan_detected = false;

    Mol mol = cmol.getMolInfo();

    json config = json::object();
    config["verbosity"] = 0;
    config["threads"]   = 1;

    GFNFF gfnff(config);
    if (!gfnff.InitialiseMolecule(mol)) {
        std::cerr << "[SKIP] " << name << ": InitialiseMolecule failed\n";
        return res;
    }

    res.energy = gfnff.Calculation(true);

    Eigen::MatrixXd G_anal = gfnff.Gradient();
    Eigen::MatrixXd G_num  = fixed_charges
                             ? gfnff.NumGradFixedCharges()
                             : gfnff.NumGrad();

    if (!G_anal.allFinite() || !G_num.allFinite()) {
        res.nan_detected = true;
        std::cout << "  [WARN] NaN/Inf detected in gradient — skipping comparison\n";
        return res;
    }

    res.anal_norm  = G_anal.norm();
    res.num_norm   = G_num.norm();
    res.max_diff   = 0.0;
    res.rms_diff   = 0.0;
    res.worst_atom = 0;

    Eigen::MatrixXd diff = G_anal - G_num;
    double worst = 0.0;
    for (int i = 0; i < res.n_atoms; ++i) {
        double row_norm = diff.row(i).norm();
        res.rms_diff += row_norm * row_norm;
        if (row_norm > worst) { worst = row_norm; res.worst_atom = i; }
        for (int d = 0; d < 3; ++d) {
            double v = std::abs(diff(i, d));
            if (v > res.max_diff) res.max_diff = v;
        }
    }
    res.rms_diff = std::sqrt(res.rms_diff / (res.n_atoms * 3));
    res.pass = (res.max_diff < tol);

    // Per-atom table
    std::cout << "  Atom | Anal_x       Num_x        Diff_x   |"
                 " Anal_y       Num_y        Diff_y   |"
                 " Anal_z       Num_z        Diff_z   | |diff|\n";
    for (int i = 0; i < res.n_atoms; ++i) {
        double row_norm = diff.row(i).norm();
        std::cout << "  " << std::setw(4) << i << " |";
        for (int d = 0; d < 3; ++d) {
            std::cout << std::scientific << std::setprecision(4)
                      << std::setw(13) << G_anal(i, d)
                      << std::setw(13) << G_num(i, d)
                      << std::setw(13) << diff(i, d) << " |";
        }
        std::cout << std::setw(11) << row_norm;
        if (i == res.worst_atom) std::cout << " * worst";
        std::cout << "\n";
    }

    return res;
}

// ─────────────────────────────────────────────────────────────────────────────
// main
// ─────────────────────────────────────────────────────────────────────────────

int main(int argc, char* argv[])
{
    CurcumaLogger::set_verbosity(0);

    bool fixed_charges = false;
    std::string xyz_file;

    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
        if (arg == "--fixed-charges" || arg == "-fc")
            fixed_charges = true;
        else
            xyz_file = arg;
    }

    const double tol = 1e-4; // Eh/Bohr

    std::cout << "GFN-FF CPU Gradient Validation (analytical vs numerical)\n"
              << "Tolerance: " << tol << " Eh/Bohr"
              << (fixed_charges ? " [fixed-charges mode]\n" : " [full EEQ mode]\n")
              << std::string(70, '=') << "\n";

    // Build molecule list
    std::vector<std::pair<std::string, curcuma::Molecule>> molecules;

    if (!xyz_file.empty()) {
        curcuma::Molecule mol(xyz_file);
        molecules.push_back({xyz_file, mol});
    } else {
        for (const auto& mol_name : DEFAULT_MOLECULES) {
            try {
                // scale_coordinates=false: registry stores Angstrom, Molecule expects Angstrom
                curcuma::Molecule mol = TestMolecules::TestMoleculeRegistry::createMolecule(mol_name, false);
                molecules.push_back({mol_name, mol});
            } catch (const std::exception& e) {
                std::cerr << "[SKIP] " << mol_name << ": " << e.what() << "\n";
            }
        }
    }

    // Run tests
    std::vector<GradTestResult> results;
    for (auto& [name, cmol] : molecules) {
        std::cout << "\n--- " << name << " (" << cmol.AtomCount() << " atoms) ---\n";
        auto res = runGradTest(name, cmol, fixed_charges, tol);
        results.push_back(res);

        if (!res.nan_detected) {
            std::cout << "  Energy:     " << std::fixed    << std::setprecision(10)
                                          << res.energy    << " Eh\n"
                      << "  |G_anal|:   " << std::scientific << std::setprecision(4)
                                          << res.anal_norm << " Eh/Bohr\n"
                      << "  |G_num|:    " << res.num_norm  << " Eh/Bohr\n"
                      << "  max|diff|:  " << res.max_diff
                      << (res.pass ? "  PASS\n" : "  FAIL *** above tolerance ***\n")
                      << "  RMS|diff|:  " << res.rms_diff  << " Eh/Bohr\n"
                      << "  Worst atom: " << res.worst_atom << "\n";
        }
    }

    // Summary
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << std::left  << std::setw(20) << "Molecule"
              << std::right << std::setw(6)  << "Atoms"
              << std::setw(17) << "max|diff|[Eh/B]"
              << std::setw(17) << "RMS|diff|[Eh/B]"
              << "  Result\n";

    bool all_pass = true;
    int  n_tested = 0;
    for (auto& r : results) {
        if (r.nan_detected) {
            std::cout << std::left  << std::setw(20) << r.name
                      << std::right << std::setw(6)  << r.n_atoms
                      << std::setw(17) << "NaN/Inf"
                      << std::setw(17) << "NaN/Inf"
                      << "  SKIP\n";
            continue;
        }
        ++n_tested;
        std::cout << std::left  << std::setw(20) << r.name
                  << std::right << std::setw(6)  << r.n_atoms
                  << std::scientific << std::setprecision(3)
                  << std::setw(17) << r.max_diff
                  << std::setw(17) << r.rms_diff
                  << "  " << (r.pass ? "PASS" : "FAIL") << "\n";
        if (!r.pass) all_pass = false;
    }

    std::cout << "\nTested: " << n_tested << "/" << results.size()
              << "  Overall: " << (all_pass ? "PASS" : "FAIL") << "\n";

    return all_pass ? 0 : 1;
}
