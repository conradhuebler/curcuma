#include "src/core/energy_calculators/qm_methods/gfnff.h"
#include "src/core/global.h"
#include <iostream>

int main() {
    std::cout << "Testing GFN-FF Implementation..." << std::endl;
    
    // Create simple H2O molecule for testing
    json gfnff_params = {
        {"method", "gfnff"},
        {"threads", 1},
        {"gradient", true},
        {"dispersion", false},  // Disable for simple test
        {"hbond", false}        // Disable for simple test
    };
    
    GFNFF gfnff(gfnff_params);
    
    // Set up H2O geometry (3 atoms: H, O, H)
    std::vector<int> atoms = {1, 8, 1};  // H, O, H
    Matrix geometry(3, 3);
    geometry << 0.0,  0.757,  0.587,   // H
                0.0,  0.0,    0.0,     // O  
                0.0, -0.757,  0.587;   // H
    
    // Initialize molecule
    if (!gfnff.InitialiseMolecule(atoms.data(), geometry.data(), 3, 0.0, 0)) {
        std::cerr << "Failed to initialize H2O molecule" << std::endl;
        return 1;
    }
    
    std::cout << "H2O molecule initialized successfully!" << std::endl;
    
    // Run calculation
    double energy = gfnff.Calculation(true, true);
    
    std::cout << "GFN-FF Energy: " << energy << " Hartree" << std::endl;
    
    // Get gradient
    Matrix gradient = gfnff.Gradient();
    std::cout << "Gradient shape: " << gradient.rows() << "x" << gradient.cols() << std::endl;
    
    std::cout << "GFN-FF Proof-of-Concept test completed!" << std::endl;
    
    return 0;
}