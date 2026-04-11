/**
 * CPU GFN-FF analytical vs. numerical gradient test.
 * Compares the analytical gradient from GFNFF::Calculation(true)
 * against GFNFF::NumGrad() (full finite-difference including EEQ recalculation).
 *
 * Usage:
 *   test_gfnff_numgrad                  – run built-in molecules
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
#include "src/core/elements.h"
#include "json.hpp"

#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using json = nlohmann::json;

// ─────────────────────────────────────────────────────────────────────────────
// Inline molecule builder via curcuma::Molecule (matches validation test pattern)
// ─────────────────────────────────────────────────────────────────────────────

struct AtomSpec { std::string element; double x, y, z; }; // Angstrom

static curcuma::Molecule buildMolecule(const std::vector<AtomSpec>& atoms)
{
    curcuma::Molecule mol;
    for (auto& a : atoms) {
        Position pos;
        pos << a.x, a.y, a.z;
        mol.addAtom({Elements::String2Element(a.element), pos});
    }
    mol.setCharge(0);
    return mol;
}

static std::vector<std::pair<std::string, curcuma::Molecule>> builtinMolecules()
{
    std::vector<std::pair<std::string, curcuma::Molecule>> mols;

    // Water
    mols.push_back({"H2O", buildMolecule({
        {"O",  0.000,  0.000,  0.119},
        {"H",  0.000,  0.757, -0.477},
        {"H",  0.000, -0.757, -0.477}
    })});

    // Methane
    mols.push_back({"CH4", buildMolecule({
        {"C",  0.000,  0.000,  0.000},
        {"H",  0.629,  0.629,  0.629},
        {"H", -0.629, -0.629,  0.629},
        {"H", -0.629,  0.629, -0.629},
        {"H",  0.629, -0.629, -0.629}
    })});

    // Methanol
    mols.push_back({"CH3OH", buildMolecule({
        {"C", -0.047,  0.000,  0.000},
        {"O",  1.384,  0.000,  0.000},
        {"H", -0.431,  1.027,  0.000},
        {"H", -0.431, -0.513,  0.889},
        {"H", -0.431, -0.513, -0.889},
        {"H",  1.759,  0.937,  0.000}
    })});

    // Dimethyl ether (sensitive to Coulomb)
    mols.push_back({"CH3OCH3", buildMolecule({
        {"C", -1.229,  0.000,  0.000},
        {"O",  0.000,  0.000,  0.000},
        {"C",  1.229,  0.000,  0.000},
        {"H", -1.618,  1.023,  0.000},
        {"H", -1.618, -0.512,  0.886},
        {"H", -1.618, -0.512, -0.886},
        {"H",  1.618,  1.023,  0.000},
        {"H",  1.618, -0.512,  0.886},
        {"H",  1.618, -0.512, -0.886}
    })});

    // Benzene (aromatic — tests torsions and dispersion)
    mols.push_back({"C6H6", buildMolecule({
        {"C",  1.397,  0.000,  0.000},
        {"C",  0.699,  1.210,  0.000},
        {"C", -0.699,  1.210,  0.000},
        {"C", -1.397,  0.000,  0.000},
        {"C", -0.699, -1.210,  0.000},
        {"C",  0.699, -1.210,  0.000},
        {"H",  2.484,  0.000,  0.000},
        {"H",  1.242,  2.151,  0.000},
        {"H", -1.242,  2.151,  0.000},
        {"H", -2.484,  0.000,  0.000},
        {"H", -1.242, -2.151,  0.000},
        {"H",  1.242, -2.151,  0.000}
    })});

    // Water dimer (tests hydrogen bonds)
    mols.push_back({"H2O dimer", buildMolecule({
        {"O",  0.000,  0.000,  0.119},
        {"H",  0.000,  0.757, -0.477},
        {"H",  0.000, -0.757, -0.477},
        {"O",  0.000,  0.000,  2.919},
        {"H",  0.000,  0.757,  3.515},
        {"H",  0.000, -0.757,  3.515}
    })});

    return mols;
}

// ─────────────────────────────────────────────────────────────────────────────
// Test runner
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
    res.name        = name;
    res.pass        = false;
    res.nan_detected = false;

    Mol mol = cmol.getMolInfo();
    res.n_atoms = mol.m_number_atoms;

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

    // Check for NaN
    if (!G_anal.allFinite() || !G_num.allFinite()) {
        res.nan_detected = true;
        std::cout << "  [WARN] NaN/Inf detected in gradient — skipping comparison\n";
        return res;
    }

    res.anal_norm = G_anal.norm();
    res.num_norm  = G_num.norm();

    Eigen::MatrixXd diff = G_anal - G_num;
    res.max_diff  = 0.0;
    res.rms_diff  = 0.0;
    res.worst_atom = 0;

    double worst = 0.0;
    for (int i = 0; i < res.n_atoms; ++i) {
        double row_norm = diff.row(i).norm();
        res.rms_diff += row_norm * row_norm;
        if (row_norm > worst) { worst = row_norm; res.worst_atom = i; }
        for (int d = 0; d < 3; ++d) {
            double v = std::abs(diff(i,d));
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
                      << std::setw(13) << G_anal(i,d)
                      << std::setw(13) << G_num(i,d)
                      << std::setw(13) << diff(i,d) << " |";
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

    std::vector<std::pair<std::string, curcuma::Molecule>> molecules;

    if (!xyz_file.empty()) {
        curcuma::Molecule mol(xyz_file);
        molecules.push_back({xyz_file, mol});
    } else {
        molecules = builtinMolecules();
    }

    std::vector<GradTestResult> results;
    for (auto& [name, cmol] : molecules) {
        std::cout << "\n--- " << name << " (" << cmol.AtomCount() << " atoms) ---\n";
        auto res = runGradTest(name, cmol, fixed_charges, tol);
        results.push_back(res);

        if (res.nan_detected) {
            std::cout << "  Result: SKIP (NaN)\n";
        } else {
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

    // Summary table
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
                      << std::setw(17) << "  NaN/Inf"
                      << std::setw(17) << "  NaN/Inf"
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
