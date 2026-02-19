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
#include <functional>
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

    void runValidation(bool run_fd = false) {
        std::cout << "=== GFN-FF Validation Runner ===" << std::endl;
        std::cout << "Molecule: " << m_ref_data["molecule"]["name"] << std::endl;
        std::cout << "Reference File: " << m_ref_path << std::endl;
        std::cout << std::endl;

        validateTopology();
        validateCharges();
        validateEnergyComponents();
        validateGradients();

        if (run_fd) {
            validateFiniteDifference();
        }

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
        const auto& ref_cn = m_ref_data["topology"]["coordination_numbers"];
        const auto& ref_hyb = m_ref_data["topology"]["hybridizations"];

        // Skip if reference topology is empty (parser didn't capture it)
        if (ref_cn.empty() || ref_hyb.empty()) {
            std::cout << "  [SKIP] Topology validation: reference data empty" << std::endl;
            return;
        }

        const auto& topo = m_gfnff->getTopologyInfo();
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
        // Prefer Phase 2 energy charges (from goed_gfnff) over Phase 1 topology charges
        // Phase 2 charges are what Curcuma returns from Charges() method
        bool has_energy_charges = m_ref_data.contains("energy_charges") && !m_ref_data["energy_charges"].empty();
        const auto& ref_charges = has_energy_charges
            ? m_ref_data["energy_charges"]
            : m_ref_data["topology"]["charges"];

        // Skip if reference charges are empty (parser didn't capture them)
        if (ref_charges.empty()) {
            std::cout << "  [SKIP] Charge validation: reference data empty" << std::endl;
            return;
        }

        if (has_energy_charges) {
            std::cout << "  [INFO] Using Phase 2 energy charges for comparison" << std::endl;
        } else {
            std::cout << "  [WARN] Using Phase 1 topology charges (Phase 2 not available)" << std::endl;
        }

        m_gfnff->Calculation(false); // Triggers EEQ
        Vector charges = m_gfnff->Charges();

        for (int i = 0; i < (int)charges.size(); ++i) {
            double q = charges(i);
            double ref = ref_charges[i];
            // Tight tolerance for energy charges (0.001 e), loose for topology charges (0.01 e)
            double tol = has_energy_charges ? 0.001 : 0.01;
            addResult("Charges", "q_" + std::to_string(i), std::abs(q - ref) < tol, q, ref, q - ref);
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
        // Fortran "etors" includes both torsions AND inversions in one loop
        double torsion_plus_inversion = m_gfnff->DihedralEnergy() + m_gfnff->InversionEnergy();
        addResult("Energy", "Torsion", std::abs(torsion_plus_inversion - (double)ref_e["torsion"]) < tol_comp, torsion_plus_inversion, ref_e["torsion"], torsion_plus_inversion - (double)ref_e["torsion"]);
        addResult("Energy", "Repulsion", std::abs(m_gfnff->RepulsionEnergy() - (double)ref_e["repulsion"]) < tol_comp, m_gfnff->RepulsionEnergy(), ref_e["repulsion"], m_gfnff->RepulsionEnergy() - (double)ref_e["repulsion"]);
        addResult("Energy", "Coulomb", std::abs(m_gfnff->CoulombEnergy() - (double)ref_e["electrostatic"]) < tol_comp, m_gfnff->CoulombEnergy(), ref_e["electrostatic"], m_gfnff->CoulombEnergy() - (double)ref_e["electrostatic"]);
        addResult("Energy", "Dispersion", std::abs(m_gfnff->DispersionEnergy() - (double)ref_e["dispersion"]) < tol_comp, m_gfnff->DispersionEnergy(), ref_e["dispersion"], m_gfnff->DispersionEnergy() - (double)ref_e["dispersion"]);

        // HB energy validation (reference key can be "hbond" or "hb")
        std::string hb_key = ref_e.contains("hbond") ? "hbond" : (ref_e.contains("hb") ? "hb" : "");
        if (!hb_key.empty()) {
            double ref_hb = ref_e[hb_key];
            double cur_hb = m_gfnff->HydrogenBondEnergy();
            addResult("Energy", "HBond", std::abs(cur_hb - ref_hb) < tol_comp, cur_hb, ref_hb, cur_hb - ref_hb);
        }
    }

    void validateGradients() {
        if (!m_ref_data.contains("gradients")) return;

        // Enable per-component gradient storage before calculation
        bool has_decomp = m_ref_data.contains("gradient_decomposition") &&
                          !m_ref_data["gradient_decomposition"].empty();
        if (has_decomp) {
            m_gfnff->setStoreGradientComponents(true);
        }

        m_gfnff->Calculation(true);
        Geometry grad = m_gfnff->Gradient();

        // Compare total norm
        double norm = 0.0;
        for (int i = 0; i < grad.rows(); ++i) {
            norm += grad.row(i).squaredNorm();
        }
        norm = std::sqrt(norm);

        double ref_norm = m_ref_data["gradients"]["norm"];
        addResult("Gradient", "TotalNorm", std::abs(norm - ref_norm) < 1e-4, norm, ref_norm, norm - ref_norm);

        // Per-atom and per-component gradient validation
        if (has_decomp) {
            validateGradientDecomposition();
            validateGradientComponents();
        }
    }

    /**
     * @brief Per-atom gradient validation against Fortran reference decomposition
     *
     * Compares per-atom TOTAL gradient (and optionally SUM) against reference data.
     * Also prints per-component diagnostic information for debugging.
     *
     * Claude Generated (February 2026): Per-atom gradient validation
     */
    void validateGradientDecomposition() {
        std::cout << std::endl << "=== Per-Atom Gradient Validation ===" << std::endl;
        const auto& ref_decomp = m_ref_data["gradient_decomposition"];

        Geometry grad = m_gfnff->Gradient();
        int natoms = grad.rows();

        // Use SUM (sum of all components) instead of TOTAL (Fortran analyzer bug: always [0,0,0])
        double tol_sum = 1e-4; // Eh/Bohr
        double max_atom_error = 0.0;
        int max_error_atom = -1;

        for (auto it = ref_decomp.begin(); it != ref_decomp.end(); ++it) {
            int atom_idx = std::stoi(it.key());
            if (atom_idx >= natoms) continue;

            const auto& comps = it.value();

            // Validate per-atom gradient against Fortran SUM (correct sum of all components)
            if (comps.contains("SUM")) {
                double rx = comps["SUM"]["x"];
                double ry = comps["SUM"]["y"];
                double rz = comps["SUM"]["z"];

                double cx = grad(atom_idx, 0);
                double cy = grad(atom_idx, 1);
                double cz = grad(atom_idx, 2);

                double err_x = cx - rx;
                double err_y = cy - ry;
                double err_z = cz - rz;
                double err_norm = std::sqrt(err_x*err_x + err_y*err_y + err_z*err_z);

                if (err_norm > max_atom_error) {
                    max_atom_error = err_norm;
                    max_error_atom = atom_idx;
                }

                bool passed = err_norm < tol_sum;
                addResult("Gradient", "Atom_" + std::to_string(atom_idx) + "_SUM",
                         passed, err_norm, 0.0, err_norm);

                if (!passed || err_norm > 1e-5) {
                    std::string elem = m_ref_data["molecule"]["atoms"][atom_idx]["element"].get<std::string>();
                    std::cout << "  Atom " << atom_idx << " (" << elem << "): "
                              << "err=" << std::scientific << std::setprecision(4) << err_norm
                              << " Curcuma=[" << std::fixed << std::setprecision(8) << cx << "," << cy << "," << cz << "]"
                              << " Ref_SUM=[" << rx << "," << ry << "," << rz << "]" << std::endl;
                }
            }

            // Print per-component diagnostic for atoms with significant gradients
            double ref_sum_norm = 0.0;
            if (comps.contains("SUM")) {
                double sx = (double)comps["SUM"]["x"], sy = (double)comps["SUM"]["y"], sz = (double)comps["SUM"]["z"];
                ref_sum_norm = std::sqrt(sx*sx + sy*sy + sz*sz);
            }
            if (ref_sum_norm > 1e-4) {
                std::string elem = m_ref_data["molecule"]["atoms"][atom_idx]["element"].get<std::string>();
                std::cout << "  Atom " << atom_idx << " (" << elem << ") component breakdown:" << std::endl;
                for (auto cit = comps.begin(); cit != comps.end(); ++cit) {
                    std::string name = cit.key();
                    double gx = cit.value()["x"], gy = cit.value()["y"], gz = cit.value()["z"];
                    double gnorm = std::sqrt(gx*gx + gy*gy + gz*gz);
                    if (gnorm > 1e-6) {
                        std::cout << "    " << std::left << std::setw(12) << name
                                  << " Norm: " << std::scientific << std::setprecision(4) << gnorm
                                  << "  [" << std::fixed << std::setprecision(8) << gx << ", " << gy << ", " << gz << "]"
                                  << std::endl;
                    }
                }
            }
        }

        std::cout << "  Max per-atom SUM error: " << std::scientific << std::setprecision(4)
                  << max_atom_error << " Eh/Bohr (atom " << max_error_atom << ")" << std::endl;
    }

    /**
     * @brief Finite-difference gradient validation (internal consistency)
     *
     * Compares analytical gradient against numerical central-difference gradient.
     * This validates that the analytical gradient implementation is correct,
     * independent of the Fortran reference.
     *
     * Claude Generated (February 2026)
     */
    void validateFiniteDifference() {
        std::cout << std::endl << "=== Finite-Difference Gradient Validation ===" << std::endl;

        // Step size for central differences
        const double h = 1e-5; // Bohr (appropriate for force fields)
        // Note: FD creates new GFNFF instances per perturbation (full topology rebuild),
        // so CN-dependent terms (bond r0, dispersion C6) are fully self-consistent.
        // Analytical gradients use CN chain-rule approximation from fixed topology.
        // Expected discrepancy: ~1e-3 for CN-dependent terms, ~1e-6 for others.
        const double tol = 1e-3; // Eh/Bohr (relaxed due to CN chain-rule approximation)

        // Get analytical gradient
        m_gfnff->Calculation(true);
        Geometry anal_grad = m_gfnff->Gradient();
        int natoms = anal_grad.rows();

        // Get reference geometry in Bohr (the GFNFF uses Bohr internally)
        const double BOHR_TO_ANGSTROM = 0.529177210903;

        double max_err = 0.0;
        int max_err_atom = -1;
        int max_err_dim = -1;

        for (int i = 0; i < natoms; ++i) {
            for (int d = 0; d < 3; ++d) {
                // Perturb atom i in dimension d by +h and -h
                // GFNFF stores geometry in Angstrom but parameters are in Bohr
                // We need to perturb in the same units as the gradient
                double h_ang = h * BOHR_TO_ANGSTROM; // Convert Bohr step to Angstrom

                // Save original coordinate
                Mol mol_plus = m_molecule.getMolInfo();
                Mol mol_minus = m_molecule.getMolInfo();

                mol_plus.m_geometry(i, d) += h_ang;
                mol_minus.m_geometry(i, d) -= h_ang;

                // Create fresh GFNFF instances for perturbed geometries
                // (reuse topology from the base calculation via geometry update)
                GFNFF gfnff_plus, gfnff_minus;
                gfnff_plus.InitialiseMolecule(mol_plus);
                gfnff_minus.InitialiseMolecule(mol_minus);

                double E_plus = gfnff_plus.Calculation(false);
                double E_minus = gfnff_minus.Calculation(false);

                double fd_grad = (E_plus - E_minus) / (2.0 * h);
                double anal = anal_grad(i, d);

                double err = std::abs(fd_grad - anal);
                if (err > max_err) {
                    max_err = err;
                    max_err_atom = i;
                    max_err_dim = d;
                }

                if (err > tol) {
                    std::string elem = m_ref_data["molecule"]["atoms"][i]["element"].get<std::string>();
                    const char* dim_name[] = {"x", "y", "z"};
                    std::cout << "  Atom " << i << " (" << elem << ") " << dim_name[d]
                              << ": FD=" << std::scientific << std::setprecision(6) << fd_grad
                              << " Anal=" << anal
                              << " err=" << err << std::endl;
                }
            }
        }

        bool passed = max_err < tol;
        addResult("FD_Gradient", "MaxError", passed, max_err, 0.0, max_err);

        const char* dim_name[] = {"x", "y", "z"};
        std::cout << "  Max FD error: " << std::scientific << std::setprecision(4) << max_err
                  << " Eh/Bohr (atom " << max_err_atom;
        if (max_err_atom >= 0)
            std::cout << " " << dim_name[max_err_dim];
        std::cout << ")" << std::endl;
        std::cout << "  " << natoms << " atoms × 3 dims = " << natoms * 3
                  << " FD evaluations (" << natoms * 6 << " energy calls)" << std::endl;
    }

    /**
     * @brief Per-component gradient validation against Fortran reference
     *
     * Maps Fortran component names to Curcuma gradient getters and compares
     * per-atom gradient vectors for each energy term.
     *
     * Claude Generated (February 2026)
     */
    void validateGradientComponents() {
        std::cout << std::endl << "=== Per-Component Gradient Validation ===" << std::endl;
        const auto& ref_decomp = m_ref_data["gradient_decomposition"];
        int natoms = m_gfnff->Gradient().rows();

        // Map Fortran component names to Curcuma getter functions
        // Fortran names: Bond, Angle, Torsion, Repulsion, ES, Disp
        struct ComponentMap {
            std::string fortran_name;
            std::string display_name;
            std::function<Matrix()> getter;
            double tolerance;
        };

        std::vector<ComponentMap> components = {
            {"Bond", "Bond", [this]() { return m_gfnff->GradientBond(); }, 1e-5},
            {"Angle", "Angle", [this]() { return m_gfnff->GradientAngle(); }, 1e-5},
            {"Torsion", "Torsion", [this]() { return m_gfnff->GradientTorsion(); }, 1e-5},
            {"Repulsion", "Repulsion", [this]() { return m_gfnff->GradientRepulsion(); }, 1e-4},
            {"ES", "Coulomb", [this]() { return m_gfnff->GradientCoulomb(); }, 1e-4},
            {"Disp", "Dispersion", [this]() { return m_gfnff->GradientDispersion(); }, 1e-4},
            {"HB", "HBond", [this]() { return m_gfnff->GradientHB(); }, 1e-4},
        };

        for (const auto& comp : components) {
            Matrix curcuma_grad = comp.getter();
            if (curcuma_grad.rows() == 0) continue;

            double max_err = 0.0;
            int max_err_atom = -1;
            int atoms_checked = 0;

            for (auto it = ref_decomp.begin(); it != ref_decomp.end(); ++it) {
                int atom_idx = std::stoi(it.key());
                if (atom_idx >= natoms) continue;
                const auto& comps = it.value();

                if (!comps.contains(comp.fortran_name)) continue;

                double rx = comps[comp.fortran_name]["x"];
                double ry = comps[comp.fortran_name]["y"];
                double rz = comps[comp.fortran_name]["z"];

                double cx = curcuma_grad(atom_idx, 0);
                double cy = curcuma_grad(atom_idx, 1);
                double cz = curcuma_grad(atom_idx, 2);

                double err = std::sqrt((cx-rx)*(cx-rx) + (cy-ry)*(cy-ry) + (cz-rz)*(cz-rz));
                atoms_checked++;

                if (err > max_err) {
                    max_err = err;
                    max_err_atom = atom_idx;
                }
            }

            bool passed = max_err < comp.tolerance;
            addResult("GradComp", comp.display_name,
                     passed, max_err, 0.0, max_err);

            std::cout << "  " << std::left << std::setw(12) << comp.display_name
                      << " max_err=" << std::scientific << std::setprecision(4) << max_err
                      << " (atom " << max_err_atom << ")"
                      << " tol=" << comp.tolerance
                      << (passed ? " PASS" : " FAIL") << std::endl;

            // Print details for failed components
            if (!passed && max_err_atom >= 0) {
                std::string elem = m_ref_data["molecule"]["atoms"][max_err_atom]["element"].get<std::string>();
                auto atom_key = std::to_string(max_err_atom);
                if (ref_decomp.contains(atom_key) && ref_decomp[atom_key].contains(comp.fortran_name)) {
                    double rx = ref_decomp[atom_key][comp.fortran_name]["x"];
                    double ry = ref_decomp[atom_key][comp.fortran_name]["y"];
                    double rz = ref_decomp[atom_key][comp.fortran_name]["z"];
                    double cx = curcuma_grad(max_err_atom, 0);
                    double cy = curcuma_grad(max_err_atom, 1);
                    double cz = curcuma_grad(max_err_atom, 2);
                    std::cout << "    Atom " << max_err_atom << " (" << elem << "): "
                              << "Curcuma=[" << std::fixed << std::setprecision(8) << cx << "," << cy << "," << cz << "] "
                              << "Ref=[" << rx << "," << ry << "," << rz << "]" << std::endl;
                }
            }
        }
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
        std::cerr << "Usage: " << argv[0] << " [--fd] <reference_json> [reference_json2 ...]" << std::endl;
        std::cerr << "  --fd  Enable finite-difference gradient validation (slow)" << std::endl;
        return 1;
    }

    bool run_fd = false;
    bool all_files_passed = true;
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--fd") {
            run_fd = true;
            continue;
        }
        try {
            GFNFFValidator validator(argv[i]);
            validator.runValidation(run_fd);
            if (!validator.allPassed()) all_files_passed = false;
        } catch (const std::exception& e) {
            std::cerr << "Error processing " << argv[i] << ": " << e.what() << std::endl;
            all_files_passed = false;
        }
        std::cout << std::endl;
    }

    return all_files_passed ? 0 : 1;
}
