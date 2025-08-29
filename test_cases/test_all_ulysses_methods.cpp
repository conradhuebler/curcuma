/*
 * Complete Ulysses Method Test Suite
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Based on ulysses-main/examples/pm6-corrected.cpp
 * Tests ALL available Ulysses methods directly using the Ulysses API
 * This bypasses the Curcuma wrapper to test the methods at their source
 */

#include "../external/ulysses-main/src/GFN.hpp"
#include "../external/ulysses-main/src/MNDOd.hpp"
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

struct MethodResult {
    std::string method_name;
    std::string correction;
    double energy;
    double time_seconds;
    bool success;
    std::string error_message;
};

void testUlyssesMethod(const std::string& method_name, const std::string& correction_type,
    const std::string& xyz_file, int charge, double level_shift,
    std::vector<MethodResult>& results)
{

    MethodResult result;
    result.method_name = method_name;
    result.correction = correction_type;
    result.success = false;

    std::cout << "\n=== Testing " << method_name;
    if (!correction_type.empty() && correction_type != "0") {
        std::cout << " with " << correction_type << " correction";
    }
    std::cout << " ===" << std::endl;

    auto start_time = std::chrono::high_resolution_clock::now();

    try {
        // Create molecule object - Ulysses can read XYZ files directly
        Molecule Mol1(xyz_file.c_str(), charge, 1); // singlet

        std::vector<size_t> atoms = Mol1.Atoms();
        matrixE geometry = Mol1.Geometry();

        std::cout << "Molecule loaded: " << atoms.size() << " atoms" << std::endl;
        std::cout << "Method: " << method_name << ", Correction: " << correction_type << std::endl;

        // Create basis set
        BSet basis(Mol1, method_name.c_str());

        // Create method object based on method name
        QCbasis* electron = nullptr;

        if (method_name == "pm6") {
            electron = new PM6(basis, Mol1, "0", correction_type);
        } else if (method_name == "am1") {
            electron = new AM1(basis, Mol1, "0", correction_type);
        } else if (method_name == "pm3") {
            electron = new PM3(basis, Mol1, "0", correction_type);
        } else if (method_name == "mndo") {
            electron = new MNDO(basis, Mol1, "0", correction_type);
        } else if (method_name == "mndod") {
            electron = new MNDOd(basis, Mol1, "0", correction_type);
        } else if (method_name == "rm1") {
            electron = new RM1(basis, Mol1, "0", correction_type);
        } else if (method_name == "pm3pddg") {
            electron = new PM3PDDG(basis, Mol1, "0", correction_type);
        } else if (method_name == "mndopddg") {
            electron = new MNDOPDDG(basis, Mol1, "0", correction_type);
        } else if (method_name == "pm3bp") {
            electron = new PM3BP(basis, Mol1, "0", correction_type);
        } else if (method_name == "ugfn2" || method_name == "gfn2") {
            electron = new GFN2(basis, Mol1); // GFN2 doesn't use corrections
        } else {
            result.error_message = "Unknown method: " + method_name;
            results.push_back(result);
            return;
        }

        // Set level shift
        electron->setEpsilonS(level_shift);

        std::cout << "Starting SCF calculation..." << std::endl;

        // Perform calculation
        electron->Calculate(0);

        // Get energy
        result.energy = electron->getEnergy();

        auto end_time = std::chrono::high_resolution_clock::now();
        result.time_seconds = std::chrono::duration<double>(end_time - start_time).count();

        result.success = true;

        std::cout << std::fixed << std::setprecision(8);
        std::cout << "✓ SUCCESS: Energy = " << result.energy << " Eh" << std::endl;
        std::cout << "  Calculation time: " << std::setprecision(3) << result.time_seconds << "s" << std::endl;

        // Get additional properties for successful calculations
        try {
            std::vector<double> HOMOorb;
            double Ehomo = electron->getHOMO(HOMOorb);
            std::vector<double> LUMOorb;
            double Elumo = electron->getLUMO(LUMOorb);

            std::cout << "  HOMO: " << std::setprecision(6) << Ehomo << " Eh" << std::endl;
            std::cout << "  LUMO: " << std::setprecision(6) << Elumo << " Eh" << std::endl;
            std::cout << "  HOMO-LUMO gap: " << std::setprecision(4) << (Elumo - Ehomo) * 27.211 << " eV" << std::endl;

        } catch (...) {
            std::cout << "  (Unable to get orbital properties)" << std::endl;
        }

        // Cleanup
        delete electron;

    } catch (const std::exception& e) {
        auto end_time = std::chrono::high_resolution_clock::now();
        result.time_seconds = std::chrono::duration<double>(end_time - start_time).count();
        result.error_message = std::string("Exception: ") + e.what();
        std::cout << "✗ FAILED: " << result.error_message << std::endl;
    } catch (...) {
        auto end_time = std::chrono::high_resolution_clock::now();
        result.time_seconds = std::chrono::duration<double>(end_time - start_time).count();
        result.error_message = "Unknown exception occurred";
        std::cout << "✗ FAILED: " << result.error_message << std::endl;
    }

    results.push_back(result);
}

int main(int argc, char* argv[])
{
    std::cout << "=== Complete Ulysses Method Test Suite ===" << std::endl;
    std::cout << "Testing all available Ulysses semi-empirical methods" << std::endl;
    std::cout << "Based on ulysses-main/examples/pm6-corrected.cpp" << std::endl;

    // Default parameters (can be overridden by command line)
    std::string xyz_file = "A.xyz";
    int charge = 0;
    double level_shift = 0.1; // Default from example

    if (argc > 1)
        xyz_file = argv[1];
    if (argc > 2)
        charge = std::atoi(argv[2]);
    if (argc > 3)
        level_shift = std::atof(argv[3]);

    std::cout << "Test molecule: " << xyz_file << std::endl;
    std::cout << "Charge: " << charge << std::endl;
    std::cout << "Level shift: " << level_shift << std::endl;

    std::vector<MethodResult> results;

    // Define all base methods (from MNDO.hpp and MNDOd.hpp analysis)
    std::vector<std::string> base_methods = {
        "gfn2", // GFN2-xTB via Ulysses
        "pm6", // PM6 method
        "am1", // AM1 method
        "pm3", // PM3 method
        "mndo", // MNDO method
        "mndod", // MNDO-d method
        "rm1", // RM1 method
        "pm3pddg", // PM3-PDDG method
        "mndopddg", // MNDO-PDDG method
        "pm3bp" // PM3-BP method
    };

    // Define correction types
    std::vector<std::string> correction_types = {
        "0", // No correction
        "D3H4X", // D3H4X dispersion correction
        "D3H+" // D3H+ correction
    };

    std::cout << "\nTesting " << base_methods.size() << " base methods × "
              << correction_types.size() << " correction modes = "
              << (base_methods.size() * correction_types.size()) << " total combinations\n"
              << std::endl;

    // Test all combinations
    for (const auto& method : base_methods) {
        for (const auto& correction : correction_types) {
            // Skip corrections for GFN2 as it doesn't support them
            if (method == "ugfn2" && correction != "0") {
                continue;
            }

            testUlyssesMethod(method, correction, xyz_file, charge, level_shift, results);
        }
    }

    // Generate summary report
    std::cout << "\n"
              << std::string(80, '=') << std::endl;
    std::cout << "=== FINAL SUMMARY REPORT ===" << std::endl;
    std::cout << std::string(80, '=') << std::endl;

    int successful = 0;
    int total = results.size();
    double total_time = 0.0;

    std::cout << std::left << std::setw(15) << "Method"
              << std::setw(10) << "Correction"
              << std::setw(18) << "Energy (Eh)"
              << std::setw(10) << "Time (s)"
              << std::setw(10) << "Status" << std::endl;
    std::cout << std::string(80, '-') << std::endl;

    for (const auto& result : results) {
        std::cout << std::left << std::setw(15) << result.method_name
                  << std::setw(10) << (result.correction == "0" ? "none" : result.correction);

        if (result.success) {
            std::cout << std::right << std::setw(18) << std::fixed << std::setprecision(8) << result.energy;
            successful++;
        } else {
            std::cout << std::right << std::setw(18) << "FAILED";
        }

        std::cout << std::right << std::setw(10) << std::fixed << std::setprecision(3) << result.time_seconds
                  << std::setw(10) << (result.success ? "✓" : "✗") << std::endl;

        total_time += result.time_seconds;

        if (!result.success && !result.error_message.empty()) {
            std::cout << "    Error: " << result.error_message << std::endl;
        }
    }

    std::cout << std::string(80, '-') << std::endl;
    std::cout << "Total: " << successful << "/" << total << " methods successful ("
              << std::fixed << std::setprecision(1) << (100.0 * successful / total) << "%)" << std::endl;
    std::cout << "Total execution time: " << std::setprecision(2) << total_time << " seconds" << std::endl;

    // Identify working methods for Curcuma integration
    std::cout << "\n=== WORKING METHODS (for Curcuma integration) ===" << std::endl;
    for (const auto& result : results) {
        if (result.success) {
            std::string full_name = result.method_name;
            if (result.correction != "0") {
                full_name += "-" + result.correction;
                // Convert to lowercase for consistency with Curcuma naming
                std::transform(full_name.begin(), full_name.end(), full_name.begin(), ::tolower);
            }
            std::cout << "\"" << full_name << "\": " << std::setprecision(8) << result.energy << " Eh" << std::endl;
        }
    }

    return successful == total ? 0 : 1;
}