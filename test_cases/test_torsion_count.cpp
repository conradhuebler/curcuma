/*
 * Torsion Count Comparison Test
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 * Claude Generated (2026-01-12) for GFN-FF torsion debugging
 */

#include <iostream>
#include <fstream>
#include "json.hpp"
#include "src/core/molecule.h"
#include "src/core/energy_calculators/ff_methods/gfnff.h"
#include "src/core/curcuma_logger.h"

using json = nlohmann::json;

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <molecule.xyz>" << std::endl;
        return 1;
    }

    // Set verbosity to maximum for detailed output
    CurcumaLogger::set_verbosity(3);

    // Load molecule
    Molecule mol(argv[1], FileIterator::Iterate::First);
    std::cout << "\n=== Molecule: " << argv[1] << " ===" << std::endl;
    std::cout << "Atoms: " << mol.AtomCount() << std::endl;

    // Create GFN-FF instance with default parameters
    json config = json::object();
    config["method"] = "gfnff";
    config["charge"] = 0;
    config["spin"] = 0;
    config["threads"] = 1;

    GFNFF gfnff(mol.Atoms(), mol.Geometry(), config);

    // Generate torsions
    std::cout << "\n=== Generating Torsions ===" << std::endl;
    json torsions = gfnff.generateGFNFFTorsions();

    std::cout << "\n=== Torsion Count ===" << std::endl;
    std::cout << "Primary torsions generated: " << torsions.size() << std::endl;

    // Print torsion details
    if (torsions.size() > 0 && torsions.size() <= 20) {
        std::cout << "\n=== Torsion Details ===" << std::endl;
        std::cout << std::setw(4) << "Idx" << "  "
                  << std::setw(3) << "i" << "  "
                  << std::setw(3) << "j" << "  "
                  << std::setw(3) << "k" << "  "
                  << std::setw(3) << "l" << "  "
                  << std::setw(2) << "n" << "  "
                  << std::setw(12) << "V (Eh)" << std::endl;

        for (size_t idx = 0; idx < torsions.size(); ++idx) {
            const auto& tor = torsions[idx];
            std::cout << std::setw(4) << idx << "  "
                      << std::setw(3) << tor["i"] << "  "
                      << std::setw(3) << tor["j"] << "  "
                      << std::setw(3) << tor["k"] << "  "
                      << std::setw(3) << tor["l"] << "  "
                      << std::setw(2) << tor["n"] << "  "
                      << std::setw(12) << std::scientific << std::setprecision(6) << tor["V"].get<double>() << std::endl;
        }
    }

    // Generate extra SP3-SP3 torsions
    std::cout << "\n=== Generating Extra SP3-SP3 Torsions ===" << std::endl;
    json extra_torsions = gfnff.generateExtraTorsions();

    std::cout << "\n=== Extra Torsion Count ===" << std::endl;
    std::cout << "Extra torsions generated: " << extra_torsions.size() << std::endl;

    // Total count
    std::cout << "\n=== TOTAL ===" << std::endl;
    std::cout << "Primary + Extra: " << torsions.size() << " + " << extra_torsions.size()
              << " = " << (torsions.size() + extra_torsions.size()) << std::endl;

    // Compare with XTB reference
    std::cout << "\n=== Comparison with XTB Reference ===" << std::endl;
    std::cout << "XTB generates: 6 primary torsions for CH3OCH3" << std::endl;
    std::cout << "Curcuma generates: " << torsions.size() << " primary torsions" << std::endl;

    if (torsions.size() == 6) {
        std::cout << "✅ MATCH: Same number of torsions as XTB!" << std::endl;
    } else if (torsions.size() > 6) {
        std::cout << "⚠️ MISMATCH: Curcuma generates " << (torsions.size() - 6) << " MORE torsions than XTB" << std::endl;
    } else {
        std::cout << "⚠️ MISMATCH: Curcuma generates " << (6 - torsions.size()) << " FEWER torsions than XTB" << std::endl;
    }

    return 0;
}
