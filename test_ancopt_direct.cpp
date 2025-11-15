#include "src/capabilities/optimizer_factory.h"
#include "src/core/energycalculator.h"
#include "src/core/molecule.h"
#include <iostream>

int main() {
    using namespace Optimization;
    using namespace curcuma;

    // Load water molecule
    Molecule mol("/home/user/curcuma/test_cases/cli/curcumaopt/01_default_uff_opt/input.xyz");
    std::cout << "Loaded " << mol.AtomCount() << " atoms\n";

    // Create energy calculator
    json calc_config = {{"method", "uff"}};
    EnergyCalculator calc("uff", calc_config);
    calc.setMolecule(mol.getMolInfo());
    
    double initial_energy = calc.CalculateEnergy(true);
    std::cout << "Initial energy: " << initial_energy << " Eh\n";

    // Test ANCOPT via factory
    std::cout << "\n=== Testing ANCOPT via Factory ===\n";
    auto ancopt = OptimizerFactory::createOptimizer(OptimizerType::ANCOPT, &calc);
    std::cout << "Created optimizer: " << ancopt->getName() << "\n";
    
    json config = {
        {"max_iterations", 20},
        {"energy_threshold", 0.1},
        {"gradient_threshold", 5e-4}
    };
    ancopt->LoadConfiguration(config);
    
    if (ancopt->InitializeOptimization(mol)) {
        std::cout << "Initialization successful!\n";
        auto result = ancopt->Optimize(false, false);
        
        if (result.success) {
            std::cout << "Optimization succeeded!\n";
            std::cout << "Final energy: " << result.final_energy << " Eh\n";
            std::cout << "Iterations: " << result.iterations_performed << "\n";
        } else {
            std::cout << "Optimization failed: " << result.error_message << "\n";
        }
    } else {
        std::cout << "Initialization failed!\n";
    }

    return 0;
}
