#include <iostream>
#include <iomanip>
#include "src/core/energycalculator.h"
#include "src/core/molecule.h"
#include "src/core/curcuma_logger.h"

int main() {
    // Set verbosity
    CurcumaLogger::set_verbosity(2);

    // Create water molecule
    std::string xyz_content = "3\nWater\nH    0.0000    0.0000    0.9584\nO    0.0000    0.0000    0.0000\nH    0.0000    0.9277   -0.2398\n";

    // Create molecule from string
    Molecule molecule;
    molecule.setXYZ(xyz_content);

    std::cout << "Testing GFN-FF energy components for water..." << std::endl;
    std::cout << "Atoms: " << molecule.AtomCount() << std::endl;

    try {
        // Create EnergyCalculator with GFN-FF
        json config = json::object();
        config["method"] = "cgfnff";
        config["verbosity"] = 3;

        EnergyCalculator ec("cgfnff", config);

        // Set molecule
        Mol mol = molecule.getMolInfo();
        ec.setMolecule(mol);

        // Calculate energy (no gradient)
        double energy = ec.CalculateEnergy(false);

        std::cout << "\n=== GFN-FF Energy Components ===" << std::endl;
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "Total Energy:    " << energy << " Eh" << std::endl;
        std::cout << "Bond Energy:     " << ec.getBondEnergy() << " Eh" << std::endl;
        std::cout << "Angle Energy:    " << ec.getAngleEnergy() << " Eh" << std::endl;
        std::cout << "Dihedral Energy: " << ec.getDihedralEnergy() << " Eh" << std::endl;
        std::cout << "Repulsion:       " << ec.getRepulsionEnergy() << " Eh" << std::endl;
        std::cout << "Coulomb Energy:  " << ec.getCoulombEnergy() << " Eh" << std::endl;
        std::cout << "Dispersion:       " << ec.getDispersionEnergy() << " Eh" << std::endl;

        // Expected values from XTB reference (approximate)
        std::cout << "\n=== Expected XTB Reference (Water) ===" << std::endl;
        std::cout << "Bond Energy:     -0.268291 Eh (expected)" << std::endl;
        std::cout << "Angle Energy:    +0.000833 Eh (expected)" << std::endl;
        std::cout << "Dihedral Energy: +0.000000 Eh (expected)" << std::endl;
        std::cout << "Repulsion:       +0.027517 Eh (expected)" << std::endl;
        std::cout << "Coulomb Energy:  -0.087576 Eh (expected)" << std::endl;
        std::cout << "Dispersion:       -0.000139 Eh (expected)" << std::endl;
        std::cout << "Total:           -0.327656 Eh (expected)" << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}