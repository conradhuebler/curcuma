#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "src/core/molecule.h"
#include "src/core/energy_calculators/qm_methods/gfnff.h"

int main() {
    // Create OH molecule
    Molecule mol;
    mol.addAtom(8, 0.0, 0.0, 0.0);  // Oxygen
    mol.addAtom(1, 0.0, 0.0, 0.9700);  // Hydrogen
    mol.setCharge(0.0);

    std::cout << "=== Testing OH EEQ Charges ===" << std::endl;

    // Create GFNFF instance
    json controller = json{};
    controller["verbosity"] = 3;
    GFNFF gfnff(controller, mol);

    // Generate topology parameters
    Vector cn = gfnff.calculateCoordinationNumbers();
    std::vector<int> hyb = gfnff.determineHybridization(cn);
    std::vector<int> rings = gfnff.findSmallestRings();

    std::cout << "Coordination Numbers:" << std::endl;
    std::cout << "  O: " << cn[0] << std::endl;
    std::cout << "  H: " << cn[1] << std::endl;

    // Calculate EEQ charges
    Vector charges = gfnff.calculateEEQCharges(cn, hyb, rings);

    std::cout << "\nEEQ Charges:" << std::endl;
    std::cout << "  O: q = " << charges[0] << " e⁻" << std::endl;
    std::cout << "  H: q = " << charges[1] << " e⁻" << std::endl;
    std::cout << "  Sum: " << charges.sum() << std::endl;

    // Expected values from reference
    std::cout << "\nExpected reference values:" << std::endl;
    std::cout << "  O: q = -0.385 e⁻" << std::endl;
    std::cout << "  H: q = +0.385 e⁻" << std::endl;

    // Check if fixed
    bool o_negative = charges[0] < 0;
    bool h_positive = charges[1] > 0;
    bool magnitude_ok = (std::abs(charges[0]) > 0.3 && std::abs(charges[1]) > 0.3);

    std::cout << "\n=== FIX STATUS ===" << std::endl;
    if (o_negative && h_positive && magnitude_ok) {
        std::cout << "✅ FIXED: Charges have correct sign AND magnitude!" << std::endl;
    } else {
        std::cout << "❌ STILL BROKEN:" << std::endl;
        if (!o_negative) std::cout << "   - Oxygen charge is not negative" << std::endl;
        if (!h_positive) std::cout << "   - Hydrogen charge is not positive" << std::endl;
        if (!magnitude_ok) std::cout << "   - Magnitude is too small (< 0.3)" << std::endl;
    }

    return 0;
}