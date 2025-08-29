/*
 * Simple Ulysses Method Test
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Simplified test focusing on the core issue: why PM3/AM1 don't work correctly
 * Tests only PM6, AM1, and PM3 to debug the specific problems identified
 */

#include "../external/ulysses-main/src/MNDOd.hpp"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <string>

void testMethod(const std::string& method_name, const std::string& correction = "0")
{
    std::cout << "\n=== Testing " << method_name;
    if (correction != "0")
        std::cout << " with " << correction;
    std::cout << " ===" << std::endl;

    try {
        auto start = std::chrono::high_resolution_clock::now();

        // Load molecule - A.xyz with charge 0, singlet
        Molecule Mol1("A.xyz", 0, 1);

        std::cout << "Molecule loaded: " << Mol1.Atoms().size() << " atoms" << std::endl;

        // Create basis set
        BSet basis(Mol1, method_name.c_str());
        std::cout << "Basis set created for method: " << method_name << std::endl;

        // Create method object
        QCbasis* electron = nullptr;

        if (method_name == "pm6") {
            electron = new PM6(basis, Mol1, "0", correction);
            std::cout << "PM6 object created with correction: " << correction << std::endl;
        } else if (method_name == "am1") {
            electron = new AM1(basis, Mol1, "0", correction);
            std::cout << "AM1 object created with correction: " << correction << std::endl;
        } else if (method_name == "pm3") {
            electron = new PM3(basis, Mol1, "0", correction);
            std::cout << "PM3 object created with correction: " << correction << std::endl;
        } else {
            std::cout << "ERROR: Unknown method " << method_name << std::endl;
            return;
        }

        // Set level shift
        electron->setEpsilonS(0.1);

        std::cout << "Starting SCF calculation..." << std::endl;

        // Calculate
        electron->Calculate(0);

        auto end = std::chrono::high_resolution_clock::now();
        double time_s = std::chrono::duration<double>(end - start).count();

        // Get results
        double energy = electron->getEnergy();

        std::cout << std::fixed << std::setprecision(8);
        std::cout << "✓ Energy: " << energy << " Eh" << std::endl;
        std::cout << "✓ Time: " << std::setprecision(3) << time_s << " seconds" << std::endl;

        // Get orbital info if possible
        try {
            std::vector<double> HOMOorb;
            double Ehomo = electron->getHOMO(HOMOorb);
            std::vector<double> LUMOorb;
            double Elumo = electron->getLUMO(LUMOorb);

            std::cout << "  HOMO: " << std::setprecision(6) << Ehomo << " Eh" << std::endl;
            std::cout << "  LUMO: " << std::setprecision(6) << Elumo << " Eh" << std::endl;
            std::cout << "  Gap:  " << std::setprecision(4) << (Elumo - Ehomo) * 27.211 << " eV" << std::endl;

        } catch (...) {
            std::cout << "  (Orbital info not available)" << std::endl;
        }

        delete electron;

    } catch (const std::exception& e) {
        std::cout << "✗ FAILED: " << e.what() << std::endl;
    } catch (...) {
        std::cout << "✗ FAILED: Unknown error" << std::endl;
    }
}

int main()
{
    std::cout << "=== Simple Ulysses Method Test ===" << std::endl;
    std::cout << "Debugging PM3/AM1 issues identified in Curcuma integration" << std::endl;

    // Test the three key methods that are causing problems
    testMethod("pm6", "0"); // Reference method (should work)
    testMethod("pm6", "D3H4X"); // PM6 with correction (should differ from above)
    testMethod("am1", "0"); // AM1 basic (falls back to UFF in Curcuma)
    testMethod("pm3", "0"); // PM3 basic (returns same energy as PM6 in Curcuma)
    testMethod("am1", "D3H4X"); // AM1 with correction
    testMethod("pm3", "D3H4X"); // PM3 with correction

    std::cout << "\n=== Analysis ===" << std::endl;
    std::cout << "Compare these results with Curcuma test_energy_methods output:" << std::endl;
    std::cout << "- PM6 should match: -328.41882169 Eh" << std::endl;
    std::cout << "- AM1 in Curcuma falls back to UFF: 1.25494377 Eh" << std::endl;
    std::cout << "- PM3 in Curcuma returns identical to PM6: -328.41882169 Eh" << std::endl;
    std::cout << "\nIf direct Ulysses gives different results, the issue is in Curcuma integration." << std::endl;

    return 0;
}