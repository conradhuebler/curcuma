/**
 * GFN-FF Unified Validation Runner
 *
 * Standardized test suite for validating GFN-FF implementation against
 * reference data generated from the Fortran gfnff analyzer.
 *
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 */

#include "src/core/energycalculator.h"
#include "src/core/molecule.h"
#include "src/core/curcuma_logger.h"
#include "src/core/energy_calculators/ff_methods/gfnff.h"
#include "src/core/elements.h"
#include "src/tools/formats.h"
#include "json.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>
#include <filesystem>
#include <memory>

using json = nlohmann::json;
using curcuma::Molecule;

struct ValidationResult {
    std::string category;
    std::string subcategory;
    bool passed;
    double value;
    double reference;
    double error;
    std::string message;
};

class GFNFFValidator {
public:
    GFNFFValidator(const std::string& ref_json_path) {
        loadReferenceData(ref_json_path);
        setupMolecule();
    }

    void runValidation() {
        std::cout << "=== GFN-FF Validation Runner ===" << std::endl;
        std::cout << "Molecule: " << m_ref_data["molecule"]["name"] << std::endl;
        std::cout << "Reference File: " << m_ref_path << std::endl;
        std::cout << std::endl;

        validateTopology();
        validateCharges();
        validateEnergyComponents();
        validateGradients();

        printSummary();
    }

    bool allPassed() const {
        for (const auto& res : m_results) {
            if (!res.passed) return false;
        }
        return true;
    }

private:
    void loadReferenceData(const std::string& path) {
        std::ifstream f(path);
        if (!f.is_open()) {
            throw std::runtime_error("Could not open reference JSON: " + path);
        }
        f >> m_ref_data;
        m_ref_path = path;
    }

    void setupMolecule() {
        m_molecule = Molecule();
        // Reference data from analyzer is in BOHR
        const double BOHR_TO_ANGSTROM = 0.529177210903;
        for (const auto& atom : m_ref_data["molecule"]["atoms"]) {
            Position pos;
            pos << (double)atom["x"] * BOHR_TO_ANGSTROM,
                   (double)atom["y"] * BOHR_TO_ANGSTROM,
                   (double)atom["z"] * BOHR_TO_ANGSTROM;
            std::string element_str = atom["element"];
            m_molecule.addAtom({Elements::String2Element(element_str), pos});
        }
        m_molecule.setCharge(0); // Default to neutral for now

        m_gfnff = std::make_unique<GFNFF>();
        m_gfnff->InitialiseMolecule(m_molecule.getMolInfo());
    }

    void validateTopology() {
        const auto& topo = m_gfnff->getTopologyInfo();
        const auto& ref_cn = m_ref_data["topology"]["coordination_numbers"];
        const auto& ref_hyb = m_ref_data["topology"]["hybridizations"];

        for (size_t i = 0; i < (size_t)topo.coordination_numbers.size(); ++i) {
            double cn = topo.coordination_numbers(i);
            double ref = ref_cn[i];
            addResult("Topology", "CN_" + std::to_string(i), std::abs(cn - ref) < 0.01, cn, ref, cn - ref);

            int hyb = topo.hybridization[i];
            int rhyb = ref_hyb[i];
            addResult("Topology", "Hyb_" + std::to_string(i), hyb == rhyb, hyb, rhyb, hyb - rhyb);
        }
    }

    void validateCharges() {
        m_gfnff->Calculation(false); // Triggers EEQ
        Vector charges = m_gfnff->Charges();
        const auto& ref_charges = m_ref_data["topology"]["charges"];

        for (int i = 0; i < (int)charges.size(); ++i) {
            double q = charges(i);
            double ref = ref_charges[i];
            // Higher tolerance for charges (0.01 e) as they are sensitive
            addResult("Charges", "q_" + std::to_string(i), std::abs(q - ref) < 0.01, q, ref, q - ref);
        }
    }

    void validateEnergyComponents() {
        double total = m_gfnff->Calculation(false);
        const auto& ref_e = m_ref_data["energy_components"];

        double tol_comp = 1e-6; // 1 µEh tolerance for individual terms
        double tol_total = 1e-4; // 100 µEh tolerance for total energy (allow some cumulative error)

        // Calculate total reference energy from components if "total" is missing
        double ref_total = 0.0;
        if (ref_e.contains("total")) {
            ref_total = ref_e["total"];
        } else {
            for (auto it = ref_e.begin(); it != ref_e.end(); ++it) {
                ref_total += (double)it.value();
            }
        }

        addResult("Energy", "Total", std::abs(total - ref_total) < tol_total, total, ref_total, total - ref_total);
        addResult("Energy", "Bond", std::abs(m_gfnff->BondEnergy() - (double)ref_e["bond"]) < tol_comp, m_gfnff->BondEnergy(), ref_e["bond"], m_gfnff->BondEnergy() - (double)ref_e["bond"]);
        addResult("Energy", "Angle", std::abs(m_gfnff->AngleEnergy() - (double)ref_e["angle"]) < tol_comp, m_gfnff->AngleEnergy(), ref_e["angle"], m_gfnff->AngleEnergy() - (double)ref_e["angle"]);
        addResult("Energy", "Torsion", std::abs(m_gfnff->DihedralEnergy() - (double)ref_e["torsion"]) < tol_comp, m_gfnff->DihedralEnergy(), ref_e["torsion"], m_gfnff->DihedralEnergy() - (double)ref_e["torsion"]);
        addResult("Energy", "Repulsion", std::abs(m_gfnff->RepulsionEnergy() - (double)ref_e["repulsion"]) < tol_comp, m_gfnff->RepulsionEnergy(), ref_e["repulsion"], m_gfnff->RepulsionEnergy() - (double)ref_e["repulsion"]);
        addResult("Energy", "Coulomb", std::abs(m_gfnff->CoulombEnergy() - (double)ref_e["electrostatic"]) < tol_comp, m_gfnff->CoulombEnergy(), ref_e["electrostatic"], m_gfnff->CoulombEnergy() - (double)ref_e["electrostatic"]);
        addResult("Energy", "Dispersion", std::abs(m_gfnff->DispersionEnergy() - (double)ref_e["dispersion"]) < tol_comp, m_gfnff->DispersionEnergy(), ref_e["dispersion"], m_gfnff->DispersionEnergy() - (double)ref_e["dispersion"]);
    }

    void validateGradients() {
        if (!m_ref_data.contains("gradients")) return;

        m_gfnff->Calculation(true);
        Geometry grad = m_gfnff->Gradient();

        // Compare norm
        double norm = 0.0;
        for (int i = 0; i < grad.rows(); ++i) {
            norm += grad.row(i).squaredNorm();
        }
        norm = std::sqrt(norm);

        double ref_norm = m_ref_data["gradients"]["norm"];
        addResult("Gradient", "Norm", std::abs(norm - ref_norm) < 1e-4, norm, ref_norm, norm - ref_norm);
    }

    void addResult(const std::string& cat, const std::string& sub, bool passed, double val, double ref, double err) {
        m_results.push_back({cat, sub, passed, val, ref, err, ""});
    }

    void printSummary() {
        std::cout << std::left << std::setw(15) << "Category" << std::setw(15) << "Subcategory"
                  << std::setw(10) << "Status" << std::setw(15) << "Value" << std::setw(15) << "Reference"
                  << std::setw(15) << "Error" << std::endl;
        std::cout << std::string(85, '-') << std::endl;

        int passed_count = 0;
        for (const auto& res : m_results) {
            std::cout << std::left << std::setw(15) << res.category << std::setw(15) << res.subcategory
                      << std::setw(10) << (res.passed ? "PASSED" : "FAILED")
                      << std::fixed << std::setprecision(8) << std::setw(15) << res.value
                      << std::setw(15) << res.reference
                      << std::setw(15) << res.error << std::endl;
            if (res.passed) passed_count++;
        }

        std::cout << std::endl;
        std::cout << "Summary: " << passed_count << "/" << m_results.size() << " checks passed." << std::endl;
    }

    json m_ref_data;
    std::string m_ref_path;
    Molecule m_molecule;
    std::unique_ptr<GFNFF> m_gfnff;
    std::vector<ValidationResult> m_results;
};

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <reference_json> [reference_json2 ...]" << std::endl;
        return 1;
    }

    bool all_files_passed = true;
    for (int i = 1; i < argc; ++i) {
        try {
            GFNFFValidator validator(argv[i]);
            validator.runValidation();
            if (!validator.allPassed()) all_files_passed = false;
        } catch (const std::exception& e) {
            std::cerr << "Error processing " << argv[i] << ": " << e.what() << std::endl;
            all_files_passed = false;
        }
        std::cout << std::endl;
    }

    return all_files_passed ? 0 : 1;
}
