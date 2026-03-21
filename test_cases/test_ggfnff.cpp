/**
 * GPU GFN-FF Validation: ggfnff vs gfnff (CPU)
 *
 * Claude Generated (March 2026): Validates that the GPU-accelerated GFN-FF
 * (ggfnff) produces energies and gradients numerically identical to the CPU
 * reference implementation (gfnff) within tight tolerances.
 *
 * Usage:
 *   test_ggfnff [reference_molecule.ref.json]
 *   test_ggfnff  (no arg: uses built-in molecules)
 *
 * Reference tolerances (from plan):
 *   |E_GPU - E_CPU| < 1e-6 Eh
 *   ||G_GPU - G_CPU||_max < 1e-5 Eh/Bohr
 *
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 */

#ifdef USE_CUDA

#include "src/core/molecule.h"
#include "src/core/curcuma_logger.h"
#include "src/core/global.h"
#include "src/core/energy_calculators/ff_methods/gfnff.h"
#include "src/core/energy_calculators/ff_methods/cuda/ff_workspace_gpu.h"
#include "src/tools/formats.h"
#include "json.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using json = nlohmann::json;
// Matrix and Vector are already defined in global.h as RowMajor types
// using Matrix = Eigen::MatrixXd;  // Would conflict (ColMajor vs RowMajor)
// using Vector = Eigen::VectorXd;

// ============================================================================
// Test utilities
// ============================================================================

static int g_passed = 0;
static int g_failed = 0;

static bool check(const std::string& name, double value, double ref,
                  double tol, const std::string& unit = "Eh")
{
    const double err = std::abs(value - ref);
    const bool ok = (err <= tol);
    if (ok) {
        std::cout << "  [PASS] " << name << ": " << std::setprecision(10) << value
                  << " (ref=" << ref << ", err=" << std::scientific << err << " " << unit << ")\n";
        ++g_passed;
    } else {
        std::cout << "  [FAIL] " << name << ": " << std::setprecision(10) << value
                  << " (ref=" << ref << ", err=" << std::scientific << err
                  << " " << unit << ", tol=" << tol << ")\n";
        ++g_failed;
    }
    return ok;
}

static void section(const std::string& title)
{
    std::cout << "\n=== " << title << " ===\n";
}

// ============================================================================
// Run one molecule: compare CPU gfnff vs GPU ggfnff
// ============================================================================

struct TestResult {
    std::string mol_name;
    bool        passed_energy   = false;
    bool        passed_gradient = false;
    double      energy_error    = 0.0;
    double      gradient_error  = 0.0;
};

TestResult runComparisonTest(const std::string& mol_name, const Mol& mol)
{
    TestResult result;
    result.mol_name = mol_name;

    constexpr double E_TOL  = 1e-6;   // Eh — total energy tolerance
    constexpr double G_TOL  = 1e-5;   // Eh/Bohr — per-component gradient tolerance

    section("CPU gfnff: " + mol_name);

    // =========================================================================
    // CPU path: run gfnff via GFNFF class
    // =========================================================================
    std::cout << "[DEBUG] Step 1: Creating CPU GFNFF instance\n" << std::flush;
    json cpu_config;
    GFNFF cpu_gfnff(cpu_config);

    std::cout << "[DEBUG] Step 2: CPU InitialiseMolecule\n" << std::flush;
    if (!cpu_gfnff.InitialiseMolecule(mol)) {
        std::cout << "  [SKIP] CPU InitialiseMolecule failed\n";
        ++g_failed;
        return result;
    }

    std::cout << "[DEBUG] Step 3: CPU Calculation\n" << std::flush;
    const double E_cpu = cpu_gfnff.Calculation(/*gradient=*/true);
    const Matrix G_cpu = cpu_gfnff.Gradient();

    std::cout << "  CPU energy: " << std::setprecision(12) << E_cpu << " Eh\n";
    std::cout << "  CPU gradient norm: " << G_cpu.norm() << " Eh/Bohr\n";

    section("GPU ggfnff: " + mol_name);

    // =========================================================================
    // GPU path: generate parameter set from CPU GFNFF, upload to GPU workspace
    // =========================================================================
    std::cout << "[DEBUG] Step 4: Creating GPU GFNFF instance\n" << std::flush;
    json gpu_config;
    GFNFF gpu_gfnff(gpu_config);

    std::cout << "[DEBUG] Step 5: GPU InitialiseMolecule\n" << std::flush;
    if (!gpu_gfnff.InitialiseMolecule(mol)) {
        std::cout << "  [SKIP] GPU InitialiseMolecule failed\n";
        ++g_failed;
        return result;
    }

    std::cout << "[DEBUG] Step 6: Consuming pending GPU params\n" << std::flush;
    std::unique_ptr<GFNFFParameterSet> pending = gpu_gfnff.consumePendingGPUParams();
    if (!pending) {
        std::cout << "  [FAIL] No pending GPU params after InitialiseMolecule\n";
        ++g_failed;
        return result;
    }
    const GFNFFParameterSet& params = *pending;
    const int natoms = mol.m_number_atoms;

    std::cout << "[DEBUG]   Parameter counts: bonds=" << params.bonds.size()
              << " angles=" << params.angles.size()
              << " dihedrals=" << params.dihedrals.size()
              << " inversions=" << params.inversions.size() << "\n" << std::flush;
    std::cout << "[DEBUG]   Non-bonded: disp=" << params.dispersions.size()
              << " coul=" << params.coulombs.size()
              << " brep=" << params.bonded_repulsions.size()
              << " nbrep=" << params.nonbonded_repulsions.size() << "\n" << std::flush;

    // Build atom type vector
    std::vector<int> atom_types = mol.m_atoms;

    std::cout << "[DEBUG] Step 7: Creating FFWorkspaceGPU (natoms=" << natoms << ")\n" << std::flush;
    std::unique_ptr<FFWorkspaceGPU> gpu_ws;
    try {
        gpu_ws = std::make_unique<FFWorkspaceGPU>(params, natoms, atom_types);
        std::cout << "[DEBUG]   GPU workspace created: disp=" << gpu_ws->dispersionCount()
                  << " bonds=" << gpu_ws->bondCount() << "\n" << std::flush;
    } catch (const std::exception& e) {
        std::cout << "  [FAIL] FFWorkspaceGPU init: " << e.what() << "\n";
        ++g_failed;
        return result;
    }

    // Inject GPU workspace and run
    std::cout << "[DEBUG] Step 8: Setting GPU workspace\n" << std::flush;
    gpu_gfnff.setGPUWorkspace(gpu_ws.get());

    std::cout << "[DEBUG] Step 9: GPU Calculation\n" << std::flush;
    const double E_gpu = gpu_gfnff.Calculation(/*gradient=*/true);
    std::cout << "[DEBUG] Step 10: Getting GPU gradient\n" << std::flush;
    const Matrix G_gpu = gpu_gfnff.Gradient();

    std::cout << "  GPU energy: " << std::setprecision(12) << E_gpu << " Eh\n";
    std::cout << "  GPU gradient norm: " << G_gpu.norm() << " Eh/Bohr\n";

    // =========================================================================
    // Compare
    // =========================================================================
    section("Comparison: " + mol_name);

    const double E_err = std::abs(E_gpu - E_cpu);
    result.energy_error = E_err;
    result.passed_energy = check("Total energy", E_gpu, E_cpu, E_TOL, "Eh");

    // Max absolute gradient component error
    double G_err = 0.0;
    if (G_gpu.rows() == G_cpu.rows() && G_gpu.cols() == G_cpu.cols()) {
        G_err = (G_gpu - G_cpu).cwiseAbs().maxCoeff();
    } else {
        std::cout << "  [FAIL] Gradient dimension mismatch\n";
        ++g_failed;
    }
    result.gradient_error = G_err;
    result.passed_gradient = check("Gradient (max comp.)", G_err, 0.0, G_TOL, "Eh/Bohr");

    return result;
}

// ============================================================================
// Inline test molecules
// ============================================================================

static Mol makeWater()
{
    Mol mol;
    mol.m_number_atoms = 3;
    mol.m_charge = 0;
    mol.m_atoms = {8, 1, 1};  // O, H, H
    mol.m_geometry.resize(3, 3);
    // Geometry in Angstrom
    mol.m_geometry << 0.000,  0.000,  0.000,
                      0.000,  0.757,  0.586,
                      0.000, -0.757,  0.586;
    mol.m_commentline = "H2O";
    return mol;
}

static Mol makeMethane()
{
    Mol mol;
    mol.m_number_atoms = 5;
    mol.m_charge = 0;
    mol.m_atoms = {6, 1, 1, 1, 1};  // C, H×4
    mol.m_geometry.resize(5, 3);
    // Geometry in Angstrom
    mol.m_geometry << 0.000,  0.000,  0.000,
                      0.629,  0.629,  0.629,
                     -0.629, -0.629,  0.629,
                     -0.629,  0.629, -0.629,
                      0.629, -0.629, -0.629;
    mol.m_commentline = "CH4";
    return mol;
}

static Mol makeEthanol()
{
    Mol mol;
    mol.m_number_atoms = 9;
    mol.m_charge = 0;
    mol.m_atoms = {6, 6, 8, 1, 1, 1, 1, 1, 1};  // C, C, O, H×6
    mol.m_geometry.resize(9, 3);
    // Geometry in Angstrom (approximate ethanol structure)
    mol.m_geometry <<  0.000,  0.000,  0.000,   // C1
                       1.530,  0.000,  0.000,   // C2
                       2.028,  1.312,  0.000,   // O
                      -0.383,  1.026,  0.000,   // H
                      -0.383, -0.513,  0.887,   // H
                      -0.383, -0.513, -0.887,   // H
                       1.913, -0.513,  0.887,   // H
                       1.913, -0.513, -0.887,   // H
                       2.956,  1.278,  0.000;   // OH
    mol.m_commentline = "ethanol";
    return mol;
}

// ============================================================================
// Load molecule from reference JSON (optional)
// ============================================================================

static Mol loadFromRefJSON(const std::string& path)
{
    std::ifstream f(path);
    if (!f.is_open())
        throw std::runtime_error("Cannot open: " + path);

    json j;
    f >> j;

    Mol mol;
    mol.m_number_atoms = j["molecule"]["natoms"].get<int>();
    mol.m_charge       = j["molecule"].value("charge", 0);
    mol.m_commentline  = j["molecule"].value("name", path);

    const auto& atoms  = j["molecule"]["atoms"];
    const auto& coords = j["molecule"]["coords"];  // flat [x0,y0,z0,x1,...]

    mol.m_atoms.resize(mol.m_number_atoms);
    mol.m_geometry.resize(mol.m_number_atoms, 3);

    for (int i = 0; i < mol.m_number_atoms; ++i) {
        mol.m_atoms[i] = atoms[i].get<int>();
        mol.m_geometry(i, 0) = coords[3*i+0].get<double>();
        mol.m_geometry(i, 1) = coords[3*i+1].get<double>();
        mol.m_geometry(i, 2) = coords[3*i+2].get<double>();
    }

    return mol;
}

// ============================================================================
// main
// ============================================================================

int main(int argc, char* argv[])
{
    CurcumaLogger::set_verbosity(0);  // Silent: suppress calculation output

    std::cout << "==================================================\n";
    std::cout << "  GFN-FF GPU Validation: ggfnff vs gfnff (CPU)\n";
    std::cout << "==================================================\n\n";

    std::vector<TestResult> results;

    if (argc > 1) {
        // Load molecules from reference JSON files passed as arguments
        for (int i = 1; i < argc; ++i) {
            try {
                Mol mol = loadFromRefJSON(argv[i]);
                results.push_back(runComparisonTest(mol.m_commentline, mol));
            } catch (const std::exception& e) {
                std::cout << "[SKIP] " << argv[i] << ": " << e.what() << "\n";
                ++g_failed;
            }
        }
    } else {
        // Use built-in test molecules
        results.push_back(runComparisonTest("H2O (water)", makeWater()));
        results.push_back(runComparisonTest("CH4 (methane)", makeMethane()));
        results.push_back(runComparisonTest("C2H5OH (ethanol)", makeEthanol()));
    }

    // =========================================================================
    // Final summary
    // =========================================================================
    std::cout << "\n==================================================\n";
    std::cout << "  SUMMARY\n";
    std::cout << "==================================================\n";
    for (const auto& r : results) {
        std::cout << "  " << r.mol_name << ": "
                  << (r.passed_energy   ? "E✓ " : "E✗ ")
                  << (r.passed_gradient ? "G✓ " : "G✗ ")
                  << "  [Eerr=" << std::scientific << std::setprecision(2) << r.energy_error
                  << " Eh, Gerr=" << r.gradient_error << " Eh/Bohr]\n";
    }
    std::cout << "\nPassed: " << g_passed << "/" << (g_passed + g_failed) << "\n";

    return (g_failed == 0) ? 0 : 1;
}

#else  // not USE_CUDA

#include <iostream>
int main() {
    std::cout << "test_ggfnff: built without USE_CUDA — skipping GPU tests\n";
    return 0;
}

#endif
