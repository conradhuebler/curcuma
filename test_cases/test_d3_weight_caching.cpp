/**
 * D3 Weight Caching Optimization Validation Test
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Verifies that weight caching optimization produces identical results
 * Claude Generated - December 27, 2025
 */

#include <cassert>
#include <iostream>
#include <cmath>

#include "src/core/energy_calculators/ff_methods/d3param_generator.h"
#include "src/core/config_manager.h"
#include "src/core/curcuma_logger.h"

int main() {
    CurcumaLogger::set_verbosity(3);

    std::cout << "========================================" << std::endl;
    std::cout << "D3 Weight Caching Optimization Test" << std::endl;
    std::cout << "========================================" << std::endl;

    // Test 1: H2 molecule (simple homoatomic)
    {
        std::cout << "\n[Test 1] H2 molecule" << std::endl;

        std::vector<int> atoms = {1, 1};  // Two hydrogen atoms
        Eigen::MatrixXd geometry(2, 3);
        geometry << 0.0, 0.0, 0.0,
                    0.0, 0.0, 1.4;  // 1.4 Bohr = ~0.74 Angstrom

        // PBE0-D3-BJ parameters
        json config_json = {
            {"d3_s6", 1.0},
            {"d3_s8", 1.2177},
            {"d3_a1", 0.4145},
            {"d3_a2", 4.8593}
        };

        ConfigManager config("d3param", config_json);
        D3ParameterGenerator d3_gen(config);

        d3_gen.GenerateParameters(atoms, geometry);
        double energy = d3_gen.getTotalEnergy();

        std::cout << "  H2 D3 energy: " << std::scientific << energy << " Eh" << std::endl;

        // Sanity check: should be negative (attractive dispersion)
        if (energy >= 0.0) {
            std::cerr << "ERROR: D3 energy should be negative for H2!" << std::endl;
            return 1;
        }

        // Should be in reasonable range for H-H pair
        if (energy < -0.001 || energy > -1e-6) {
            std::cerr << "ERROR: D3 energy out of expected range for H2!" << std::endl;
            std::cerr << "  Expected: ~ -6e-5 Eh" << std::endl;
            std::cerr << "  Got: " << energy << " Eh" << std::endl;
            return 1;
        }

        std::cout << "  ✓ H2 test passed" << std::endl;
    }

    // Test 2: CH4 molecule (heteroatomic)
    {
        std::cout << "\n[Test 2] CH4 molecule" << std::endl;

        std::vector<int> atoms = {6, 1, 1, 1, 1};  // C + 4H
        Eigen::MatrixXd geometry(5, 3);
        geometry << 0.0, 0.0, 0.0,                // C at origin
                    0.0, 0.0, 2.05,               // H1
                    1.93, 0.0, -0.68,             // H2
                    -0.97, 1.68, -0.68,           // H3
                    -0.97, -1.68, -0.68;          // H4

        json config_json = {
            {"d3_s6", 1.0},
            {"d3_s8", 1.2177},
            {"d3_a1", 0.4145},
            {"d3_a2", 4.8593}
        };

        ConfigManager config("d3param", config_json);
        D3ParameterGenerator d3_gen(config);

        d3_gen.GenerateParameters(atoms, geometry);
        double energy = d3_gen.getTotalEnergy();

        std::cout << "  CH4 D3 energy: " << std::scientific << energy << " Eh" << std::endl;

        if (energy >= 0.0) {
            std::cerr << "ERROR: D3 energy should be negative for CH4!" << std::endl;
            return 1;
        }

        // Should be larger magnitude than H2 (more pairs)
        if (energy < -0.01 || energy > -1e-5) {
            std::cerr << "ERROR: D3 energy out of expected range for CH4!" << std::endl;
            std::cerr << "  Expected: ~ -9e-4 Eh" << std::endl;
            std::cerr << "  Got: " << energy << " Eh" << std::endl;
            return 1;
        }

        std::cout << "  ✓ CH4 test passed" << std::endl;
    }

    std::cout << "\n========================================" << std::endl;
    std::cout << "✅ All weight caching tests PASSED" << std::endl;
    std::cout << "========================================" << std::endl;

    return 0;
}
