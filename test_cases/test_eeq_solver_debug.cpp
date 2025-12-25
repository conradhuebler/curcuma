/*
 * EEQSolver Debugging Test
 * Tests EEQ charge calculation in isolation to identify zero charge issue
 * Claude Generated - December 25, 2025
 */

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

#include "src/core/config_manager.h"
#include "src/core/energy_calculators/ff_methods/eeq_solver.h"
#include "src/core/energy_calculators/ff_methods/cn_calculator.h"
#include "src/core/molecule.h"
#include "src/core/curcuma_logger.h"
#include "json.hpp"

using json = nlohmann::json;
using Eigen::MatrixXd;
using Eigen::VectorXd;

int main(int argc, char *argv[]) {
    // Test 1: H₂ molecule (simplest case)
    std::cout << "\n=== TEST 1: H2 Molecule (H at 0,0,0 and 0,0,0.74) ===" << std::endl;

    std::vector<int> atoms_h2 = {1, 1};  // Two hydrogens

    // Geometry in Bohr (0.74 Å ≈ 1.40 Bohr)
    MatrixXd geom_h2(2, 3);
    geom_h2 << 0.0, 0.0, 0.0,
              0.0, 0.0, 1.40;

    std::cout << "Atoms: " << atoms_h2[0] << ", " << atoms_h2[1] << std::endl;
    std::cout << "Geometry (Bohr):\n" << geom_h2 << std::endl;

    // Test CN calculation first
    std::cout << "\n--- CN Calculation ---" << std::endl;
    std::vector<double> cn_values = CNCalculator::calculateGFNFFCN(atoms_h2, geom_h2);
    std::cout << "CN[0]: " << cn_values[0] << std::endl;
    std::cout << "CN[1]: " << cn_values[1] << std::endl;

    // Test EEQSolver
    std::cout << "\n--- EEQSolver Test ---" << std::endl;
    json eeq_params = {
        {"max_iterations", 50},
        {"convergence_threshold", 1e-6},
        {"verbosity", 3},
        {"calculate_cn", true}
    };

    ConfigManager eeq_config("eeq_solver", eeq_params);
    EEQSolver solver(eeq_config);

    // Convert CN to Eigen::VectorXd
    VectorXd cn_eigen(2);
    cn_eigen << cn_values[0], cn_values[1];

    std::cout << "\nCalling EEQSolver::calculateCharges for H2..." << std::endl;
    VectorXd charges = solver.calculateCharges(atoms_h2, geom_h2, 0, &cn_eigen);

    std::cout << "\nResults:" << std::endl;
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Charge[0]: " << charges(0) << std::endl;
    std::cout << "Charge[1]: " << charges(1) << std::endl;
    std::cout << "Total charge: " << charges.sum() << std::endl;

    // Check if all charges are zero
    double max_abs_charge = charges.cwiseAbs().maxCoeff();
    std::cout << "Max absolute charge: " << max_abs_charge << std::endl;

    if (max_abs_charge < 1e-10) {
        std::cerr << "\n❌ ERROR: All charges are ZERO! This is the bug." << std::endl;
        std::cerr << "Possible causes:" << std::endl;
        std::cerr << "  1. EEQ parameters (chi, gam, alpha) not loaded correctly" << std::endl;
        std::cerr << "  2. Linear system is singular or near-singular" << std::endl;
        std::cerr << "  3. CN calculation returns unexpected values" << std::endl;
        std::cerr << "  4. Matrix setup has a bug (all zeros in matrix)" << std::endl;
    } else {
        std::cout << "\n✅ SUCCESS: Charges are non-zero!" << std::endl;
    }

    // Test 2: CH₄ molecule
    std::cout << "\n\n=== TEST 2: CH4 Molecule ===" << std::endl;

    std::vector<int> atoms_ch4 = {6, 1, 1, 1, 1};  // C + 4 H
    MatrixXd geom_ch4(5, 3);
    geom_ch4 << 0.000000,  0.000000,  0.000000,   // C
               1.089000,  1.089000,  1.089000,   // H
              -1.089000, -1.089000,  1.089000,   // H
              -1.089000,  1.089000, -1.089000,   // H
               1.089000, -1.089000, -1.089000;   // H

    std::cout << "Testing CH4 (5 atoms)..." << std::endl;

    std::vector<double> cn_ch4 = CNCalculator::calculateGFNFFCN(atoms_ch4, geom_ch4);
    std::cout << "CN values: ";
    for (auto cn : cn_ch4) std::cout << cn << " ";
    std::cout << std::endl;

    VectorXd cn_ch4_eigen = VectorXd::Map(cn_ch4.data(), cn_ch4.size());
    VectorXd charges_ch4 = solver.calculateCharges(atoms_ch4, geom_ch4, 0, &cn_ch4_eigen);

    std::cout << "Charges (C + 4H): ";
    for (int i = 0; i < charges_ch4.size(); ++i) {
        std::cout << charges_ch4(i) << " ";
    }
    std::cout << std::endl;
    std::cout << "Total charge: " << charges_ch4.sum() << std::endl;

    double max_abs_ch4 = charges_ch4.cwiseAbs().maxCoeff();
    if (max_abs_ch4 < 1e-10) {
        std::cerr << "❌ CH4: All charges ZERO" << std::endl;
    } else {
        std::cout << "✅ CH4: Charges OK" << std::endl;
    }

    return 0;
}
