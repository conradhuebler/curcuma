/*
 * GFN-FF Term Validator - Individual Energy Term Validation Tool
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated - GFN-FF 100% Accuracy Implementation Plan
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <nlohmann/json.hpp>

#include "src/core/energy_calculators/qm_methods/gfnff.h"
#include "src/core/molecule.h"
#include "src/core/curcuma_logger.h"

using json = nlohmann::json;

struct EnergyResult {
    std::string system;
    std::string term;
    double calculated;
    double reference;
    double difference;
    double relative_error;
    bool passed;
};

class GFNFFTValidator {
private:
    std::vector<EnergyResult> m_results;
    const double EPSILON = 1e-6;  // Target tolerance: 1e-6 Eh

public:
    struct TestMolecule {
        std::string name;
        std::string xyz_file;
        std::string reference_file;
    };

    void validateSystem(const TestMolecule& test) {
        std::cout << "=== " << test.name << " Validation ===" << std::endl;

        // Load reference data
        std::ifstream ref_file(test.reference_file);
        json reference;
        ref_file >> reference;

        // Load molecule
        Mol mol("molecules/" + test.xyz_file);

        // Run GFN-FF calculation
        GFNFF gfnff;
        gfnff.setVerbosity(0);  // Silent mode for clean output
        if (!gfnff.InitialiseMolecule(mol)) {
            std::cout << "ERROR: Failed to initialize " << test.name << std::endl;
            return;
        }

        double total_energy = gfnff.CalculateEnergy(false);

        // Get energy breakdown (need to extract from logging or add getter methods)
        json breakdown = extractEnergyBreakdown(gfnff);

        // Validate each term
        const std::map<std::string, std::string> term_map = {
            {"bond", "Bond energy"},
            {"angle", "Angle energy"},
            {"torsion", "Torsion energy"},
            {"repulsion", "Repulsion energy"},
            {"electrostatic", "Electrostatic"},
            {"dispersion", "Dispersion"}
        };

        for (const auto& [term_key, term_name] : term_map) {
            validateTerm(test.name, term_key, breakdown[term_key],
                         reference["reference_data"]["energy_decomposition"][term_key]);
        }

        std::cout << std::endl;
    }

    void validateTerm(const std::string& system, const std::string& term,
                      double calculated, const json& reference) {
        EnergyResult result;
        result.system = system;
        result.term = term;
        result.calculated = calculated;
        result.reference = reference["value"];
        result.difference = std::abs(calculated - result.reference);
        result.relative_error = result.difference / std::abs(result.reference);
        result.passed = result.difference < EPSILON;

        m_results.push_back(result);
    }

    void printValidationReport() {
        std::cout << "\n=== GFN-FF VALIDATION REPORT ===" << std::endl;
        std::cout << "Target tolerance: 1e-6 Hartree per term" << std::endl;
        std::cout << "Tested systems: HH, CH4, H2O" << std::endl << std::endl;

        // Summary by system
        std::map<std::string, int> passed_by_system = {
            {"HH", 0}, {"CH4", 0}, {"H2O", 0}
        };
        std::map<std::string, int> total_by_system = {
            {"HH", 0}, {"CH4", 0}, {"H2O", 0}
        };

        for (const auto& result : m_results) {
            if (result.passed) {
                passed_by_system[result.system]++;
            }
            total_by_system[result.system]++;
        }

        bool all_passed = true;
        for (const auto& [system, passed] : passed_by_system) {
            int total = total_by_system[system];
            bool system_passed = (passed == total);
            if (system_passed) {
                std::cout << "✅ " << system << ": " << passed << "/" << total << " tests passed" << std::endl;
            } else {
                std::cout << "❌ " << system << ": " << passed << "/" << total << " tests passed" << std::endl;
                all_passed = false;
            }
        }

        std::cout << std::endl;
        if (all_passed) {
            std::cout << "🎯 SUCCESS: All energy terms within 1e-6 Hartree tolerance!" << std::endl;
            std::cout << "💰 Ready to claim the $200 bounty!" << std::endl;
        } else {
            std::cout << "🔧 ACTION NEEDED: Some terms exceed tolerance limits" << std::endl;
            printDetailedBreakdown();
        }
    }

    void printDetailedBreakdown() {
        std::cout << "\n=== DETAILED BREAKDOWN ===" << std::endl;
        std::cout << std::setw(8) << "System"
                  << std::setw(12) << "Term"
                  << std::setw(15) << "Calculated"
                  << std::setw(15) << "Reference"
                  << std::setw(15) << "Difference"
                  << std::setw(12) << "Rel Error"
                  << "  Status" << std::endl;
        std::cout << std::string(84, '-') << std::endl;

        for (const auto& result : m_results) {
            std::cout << std::setw(8) << result.system
                      << std::setw(12) << result.term
                      << std::setw(15) << std::fixed << std::setprecision(9) << result.calculated
                      << std::setw(15) << result.reference
                      << std::setw(15) << result.difference
                      << std::setw(12) << std::scientific << result.relative_error;

            if (result.passed) {
                std::cout << "  ✅ PASS" << std::endl;
            } else {
                std::cout << "  ❌ FAIL" << std::endl;
            }
        }
        std::cout << std::defaultfloat << std::endl;

        // Identify largest issues
        auto max_diff = std::max_element(m_results.begin(), m_results.end(),
            [](const EnergyResult& a, const EnergyResult& b) {
                return a.difference < b.difference;
            });

        std::cout << "🔥 LARGEST DISCREPANCY:" << std::endl;
        std::cout << "   System: " << max_diff->system << std::endl;
        std::cout << "   Term: " << max_diff->term << std::endl;
        std::cout << "   Difference: " << max_diff->difference << " Hartree" << std::endl;
        std::cout << "   Factor needed: " << (max_diff->calculated / max_diff->reference) << "x" << std::endl;
    }

private:
    json extractEnergyBreakdown(GFNFF& gfnff) {
        // TODO: Add getter methods to GFNFF class for individual energy terms
        // For now, we need to parse this from logging or modify GFNFF to return breakdown
        json breakdown;

        // Placeholder - will need to implement based on actual GFNFF interface
        breakdown["bond"] = 0.0;           // Get from gfnff.getBondEnergy()
        breakdown["angle"] = 0.0;          // Get from gfnff.getAngleEnergy()
        breakdown["torsion"] = 0.0;        // Get from gfnff.getTorsionEnergy()
        breakdown["repulsion"] = 0.0;      // Get from gfnff.getRepulsionEnergy()
        breakdown["electrostatic"] = 0.0;  // Get from gfnff.getElectrostaticEnergy()
        breakdown["dispersion"] = 0.0;     // Get from gfnff.getDispersionEnergy()

        return breakdown;
    }
};

int main() {
    std::cout << "🔬 GFN-FF 100% Accuracy Validator" << std::endl;
    std::cout << "Validating individual energy terms against XTB 6.6.1 reference" << std::endl;
    std::cout << "Target: 1e-6 Hartree tolerance per term" << std::endl << std::endl;

    GFNFFTValidator validator;

    std::vector<GFNFFTValidator::TestMolecule> tests = {
        {"HH", "dimers/HH.xyz", "golden_references/gfnff_hh.json"},
        {"CH4", "larger/CH4.xyz", "golden_references/gfnff_ch4.json"},
        {"H2O", "trimers/water.xyz", "golden_references/gfnff_h2o.json"}
    };

    for (const auto& test : tests) {
        validator.validateSystem(test);
    }

    validator.printValidationReport();

    return 0;
}