/*
 * <Test Optimizer Factory - Direct Test of New Optimizer Infrastructure>
 * Copyright (C) 2025 Claude AI - Generated Code
 *
 * Quick test to verify OptimizerFactory, ANCOPT, and LBFGSPP work correctly
 */

#include "src/capabilities/optimizer_factory.h"
#include "src/core/energycalculator.h"
#include "src/core/molecule.h"
#include "src/core/curcuma_logger.h"
#include <iostream>

using namespace Optimization;
using namespace curcuma;

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cout << "Usage: " << argv[0] << " <xyz_file> <method>" << std::endl;
        std::cout << "Example: " << argv[0] << " molecule.xyz uff" << std::endl;
        return 1;
    }

    std::string xyz_file = argv[1];
    std::string method = argv[2];

    try {
        // Set verbosity
        CurcumaLogger::set_verbosity(3);

        CurcumaLogger::header("Optimizer Factory Test");
        CurcumaLogger::info(fmt::format("Testing file: {}", xyz_file));
        CurcumaLogger::info(fmt::format("QM Method: {}", method));

        // Load molecule
        Molecule mol(xyz_file);
        CurcumaLogger::success(fmt::format("Loaded molecule: {} atoms", mol.AtomCount()));

        // Create energy calculator
        json calc_config;
        calc_config["method"] = method;
        EnergyCalculator energy_calc(method, calc_config);
        energy_calc.setMolecule(mol.getMolInfo());

        // Test initial energy
        double initial_energy = energy_calc.CalculateEnergy(true);
        CurcumaLogger::energy_abs(initial_energy, "Initial energy");

        // Test all available optimizers
        std::vector<OptimizerType> optimizers_to_test = {
            OptimizerType::LBFGSPP,
            OptimizerType::ANCOPT
        };

        for (auto opt_type : optimizers_to_test) {
            CurcumaLogger::header(fmt::format("Testing {}", optimizerTypeToString(opt_type)));

            // Create optimizer copy for this test
            Molecule test_mol = mol;

            // Create optimizer via factory
            auto optimizer = OptimizerFactory::createOptimizer(opt_type, &energy_calc);

            if (!optimizer) {
                CurcumaLogger::error(fmt::format("Failed to create {} optimizer",
                    optimizerTypeToString(opt_type)));
                continue;
            }

            // Configure optimizer
            json opt_config;
            opt_config["max_iterations"] = 50;
            opt_config["energy_threshold"] = 0.1; // kJ/mol
            opt_config["gradient_threshold"] = 5e-4; // Eh/Bohr
            opt_config["write_trajectory"] = false;
            opt_config["verbose"] = true;

            optimizer->LoadConfiguration(opt_config);

            // Initialize optimization
            if (!optimizer->InitializeOptimization(test_mol)) {
                CurcumaLogger::error("Optimizer initialization failed");
                continue;
            }

            // Perform optimization
            auto result = optimizer->Optimize(false, true);

            // Report results
            if (result.success) {
                CurcumaLogger::success(fmt::format("{} optimization completed!",
                    optimizerTypeToString(opt_type)));
                CurcumaLogger::energy_abs(result.final_energy, "Final energy");
                CurcumaLogger::param("Iterations", result.iterations_performed);
                CurcumaLogger::param("Time", fmt::format("{:.3f} seconds", result.optimization_time_seconds));
                CurcumaLogger::param("Final gradient norm", fmt::format("{:.6e} Eh/Bohr", result.final_gradient_norm));
            } else {
                CurcumaLogger::error(fmt::format("{} optimization failed: {}",
                    optimizerTypeToString(opt_type), result.error_message));
            }

            std::cout << "\n";
        }

        CurcumaLogger::success("All optimizer tests completed!");
        return 0;

    } catch (const std::exception& e) {
        CurcumaLogger::error(fmt::format("Test failed with exception: {}", e.what()));
        return 1;
    }
}
