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
    GFNFFValidator(const std::string& ref_json_path, bool rep_diag = false,
                   const std::string& solve_method = "")
        : m_rep_diag(rep_diag), m_solve_method(solve_method)
    {
        loadReferenceData(ref_json_path);
        setupMolecule();
    }

    void runValidation(bool run_fd = false, bool run_charge_inject = false,
                       bool run_torsion_detail = false, bool run_disp_diag = false) {
        std::cout << "=== GFN-FF Validation Runner ===" << std::endl;
        std::cout << "Molecule: " << m_ref_data["molecule"]["name"] << std::endl;
        std::cout << "Reference File: " << m_ref_path << std::endl;
        std::cout << std::endl;

        validateTopology();
        validateCharges();
        validateEnergyComponents();
        validateBondParameters();
        validateGradients();

        if (run_fd) {
            validateFiniteDifference();
        }

        if (run_charge_inject) {
            validateChargeInjection();
        }

        if (run_torsion_detail) {
            validateTorsionDetail();
        }

        if (run_disp_diag) {
            validateDispersionDiagnostic();
        }

        printSummary();
    }

    bool allPassed() const {
        for (const auto& res : m_results) {
            if (!res.passed) return false;
        }
        return true;
    }

    GFNFF* getGFNFF() { return m_gfnff.get(); }

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

        // Claude Generated (March 2026): Pass solve_method to GFNFF via JSON parameters
        if (!m_solve_method.empty()) {
            json params;
            params["eeq_solver"]["solve_method"] = m_solve_method;
            m_gfnff = std::make_unique<GFNFF>(params);
        } else {
            m_gfnff = std::make_unique<GFNFF>();
        }
        if (m_rep_diag) m_gfnff->setRepDiag(true);
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
            // Charge tolerance scales with system size (EEQ precision limit for large molecules)
            // Base: 0.001 e for small molecules, scales as max(0.001, natoms * 1.2e-5)
            int natoms = (int)charges.size();
            double tol = has_energy_charges ? std::max(0.001, natoms * 1.2e-5) : 0.01;
            addResult("Charges", "q_" + std::to_string(i), std::abs(q - ref) < tol, q, ref, q - ref);
        }
    }

    void validateEnergyComponents() {
        double total = m_gfnff->Calculation(false);
        const auto& ref_e = m_ref_data["energy_components"];

        // Energy tolerance scales with system size: larger molecules accumulate
        // more numerical precision loss in EEQ charges and CN derivatives.
        // Base: 1 µEh for small molecules, scales as max(1e-6, natoms * 1e-8)
        int natoms = m_gfnff->Charges().size();
        double tol_comp = std::max(1e-6, natoms * 1e-8);
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
        std::cout << "  Bonded Rep:    " << std::setw(14) << std::fixed << std::setprecision(8) << m_gfnff->BondedRepulsionEnergy() << " Eh\n";
        std::cout << "  Nonbonded Rep: " << std::setw(14) << std::fixed << std::setprecision(8) << m_gfnff->NonbondedRepulsionEnergy() << " Eh\n";
        std::cout << "  Rep sum check: " << std::setw(14) << std::fixed << std::setprecision(8) << (m_gfnff->BondedRepulsionEnergy() + m_gfnff->NonbondedRepulsionEnergy()) << " Eh (should == RepulsionEnergy)\n";
        // Coulomb uses same tolerance as other components (piadr fix Mar 2026 eliminated large-system errors)
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

    /**
     * @brief Per-bond parameter comparison for diagnosing bond energy errors
     *
     * Compares per-bond fc, alpha, r0 between Curcuma and Fortran reference,
     * then computes per-bond energy decomposition to identify error sources.
     *
     * Claude Generated (March 2026): Bond energy error diagnosis
     */
    void validateBondParameters() {
        if (!m_ref_data.contains("bonds") || m_ref_data["bonds"].empty()) {
            std::cout << "  [SKIP] Bond parameter validation: no reference bond data" << std::endl;
            return;
        }

        std::cout << std::endl << "=== Per-Bond Parameter Comparison ===" << std::endl;

        // Get Curcuma's bond parameters
        json curcuma_bonds = m_gfnff->getBondParameters();
        const auto& ref_bonds = m_ref_data["bonds"];

        // Build lookup: (min_atom, max_atom) -> curcuma bond index
        std::map<std::pair<int,int>, int> curcuma_bond_map;
        for (int b = 0; b < (int)curcuma_bonds.size(); ++b) {
            int ai = curcuma_bonds[b]["i"];
            int aj = curcuma_bonds[b]["j"];
            auto key = std::make_pair(std::min(ai, aj), std::max(ai, aj));
            curcuma_bond_map[key] = b;
        }

        // Get D3 CN for dynamic r0 computation
        // The GFNFF object stores these internally after Calculation()
        // We'll use the bond's dynamic r0 parameters to compute r0

        // Get current geometry for distance computation
        // NOTE: GFN-FF works in Bohr internally, so convert Angstrom geometry
        const auto& mol = m_molecule.getMolInfo();
        const double BOHR_TO_ANGSTROM = 0.529177210903;
        const double ANG_TO_BOHR = 1.0 / BOHR_TO_ANGSTROM;

        double total_energy_ref = 0.0;
        double total_energy_curcuma = 0.0;
        double total_fc_diff = 0.0;
        double total_r0_diff = 0.0;

        int bonds_matched = 0;
        int bonds_with_fc_diff = 0;
        int bonds_with_r0_diff = 0;

        // Header for CSV diagnostic
        std::cout << std::setw(4) << "Bond"
                  << std::setw(6) << "At_i"
                  << std::setw(6) << "At_j"
                  << std::setw(5) << "Type"
                  << std::setw(13) << "fc_cur"
                  << std::setw(13) << "fc_ref"
                  << std::setw(12) << "fc_diff"
                  << std::setw(12) << "alp_cur"
                  << std::setw(12) << "alp_ref"
                  << std::setw(12) << "alp_diff"
                  << std::setw(11) << "r0_cur"
                  << std::setw(11) << "r0_ref"
                  << std::setw(11) << "r0_diff"
                  << std::setw(13) << "E_diff"
                  << std::endl;
        std::cout << std::string(140, '-') << std::endl;

        for (int rb = 0; rb < (int)ref_bonds.size(); ++rb) {
            int ref_i = ref_bonds[rb]["atoms"][0];
            int ref_j = ref_bonds[rb]["atoms"][1];
            double ref_r0_ang = (double)ref_bonds[rb]["r0"];        // Å (from Fortran)
            double ref_r0 = ref_r0_ang * ANG_TO_BOHR;               // Convert to Bohr for comparison
            double ref_alpha = ref_bonds[rb]["kbond"];       // vbond(2) = alpha
            double ref_fc = ref_bonds[rb]["prefactor"];      // vbond(3) = fc

            auto key = std::make_pair(std::min(ref_i, ref_j), std::max(ref_i, ref_j));
            auto it = curcuma_bond_map.find(key);
            if (it == curcuma_bond_map.end()) {
                std::cout << "  [MISS] Reference bond " << ref_i << "-" << ref_j << " not found in Curcuma" << std::endl;
                continue;
            }

            int cb = it->second;
            double cur_fc = curcuma_bonds[cb]["fc"];          // force constant (prefactor)
            double cur_alpha = curcuma_bonds[cb]["exponent"];  // alpha
            double cur_r0 = curcuma_bonds[cb]["r0_ij"];        // static initialization r0 (Å)

            // Compute actual distance from geometry in Bohr (GFN-FF internal units)
            int ai = curcuma_bonds[cb]["i"];
            int aj = curcuma_bonds[cb]["j"];
            double dx = (mol.m_geometry(ai, 0) - mol.m_geometry(aj, 0)) * ANG_TO_BOHR;
            double dy = (mol.m_geometry(ai, 1) - mol.m_geometry(aj, 1)) * ANG_TO_BOHR;
            double dz = (mol.m_geometry(ai, 2) - mol.m_geometry(aj, 2)) * ANG_TO_BOHR;
            double rij = std::sqrt(dx*dx + dy*dy + dz*dz); // Bohr

            // Compute per-bond energy with each set of parameters
            double dr_cur = rij - cur_r0;
            double E_cur = cur_fc * std::exp(-cur_alpha * dr_cur * dr_cur);

            double dr_ref = rij - ref_r0;
            double E_ref = ref_fc * std::exp(-ref_alpha * dr_ref * dr_ref);

            double E_diff = E_cur - E_ref;
            total_energy_ref += E_ref;
            total_energy_curcuma += E_cur;

            double fc_diff = cur_fc - ref_fc;
            double alpha_diff = cur_alpha - ref_alpha;
            double r0_diff = cur_r0 - ref_r0;

            if (std::abs(fc_diff) > 1e-10) {
                bonds_with_fc_diff++;
                total_fc_diff += fc_diff;
            }
            if (std::abs(r0_diff) > 1e-6) {
                bonds_with_r0_diff++;
                total_r0_diff += r0_diff;
            }
            bonds_matched++;

            // Determine bond type for display
            int z_i = curcuma_bonds[cb]["z_i"];
            int z_j = curcuma_bonds[cb]["z_j"];
            std::string btype;
            if (z_i == 1 || z_j == 1) {
                int heavy_z = (z_i == 1) ? z_j : z_i;
                if (heavy_z == 6) btype = "H-C";
                else if (heavy_z == 8) btype = "H-O";
                else btype = "H-?";
            } else if (z_i == 6 && z_j == 6) btype = "C-C";
            else if ((z_i == 6 && z_j == 8) || (z_i == 8 && z_j == 6)) btype = "C-O";
            else btype = "?-?";

            // Print only bonds with significant differences or all if few bonds
            bool significant = std::abs(E_diff) > 1e-6 || std::abs(fc_diff) > 1e-8 || std::abs(r0_diff) > 1e-5;
            if (significant || ref_bonds.size() <= 80) {
                std::cout << std::setw(4) << rb
                          << std::setw(6) << ref_i
                          << std::setw(6) << ref_j
                          << std::setw(5) << btype
                          << std::fixed << std::setprecision(9)
                          << std::setw(13) << cur_fc
                          << std::setw(13) << ref_fc
                          << std::scientific << std::setprecision(3)
                          << std::setw(12) << fc_diff
                          << std::fixed << std::setprecision(6)
                          << std::setw(12) << cur_alpha
                          << std::setw(12) << ref_alpha
                          << std::scientific << std::setprecision(3)
                          << std::setw(12) << alpha_diff
                          << std::fixed << std::setprecision(6)
                          << std::setw(11) << cur_r0
                          << std::setw(11) << ref_r0
                          << std::scientific << std::setprecision(3)
                          << std::setw(11) << r0_diff
                          << std::setw(13) << E_diff
                          << std::endl;
            }
        }

        // Summary
        std::cout << std::string(140, '-') << std::endl;
        std::cout << std::fixed << std::setprecision(9)
                  << "Total bond energy (Curcuma params): " << total_energy_curcuma << " Eh" << std::endl
                  << "Total bond energy (Fortran params): " << total_energy_ref << " Eh" << std::endl
                  << "Difference (Curcuma - Fortran):     " << std::scientific << std::setprecision(6)
                  << (total_energy_curcuma - total_energy_ref) << " Eh ("
                  << std::fixed << std::setprecision(3) << (total_energy_curcuma - total_energy_ref) * 1000.0 << " mEh)" << std::endl;
        std::cout << "Bonds matched: " << bonds_matched << "/" << ref_bonds.size() << std::endl;
        std::cout << "Bonds with fc diff:  " << bonds_with_fc_diff << std::endl;
        std::cout << "Bonds with r0 diff:  " << bonds_with_r0_diff << std::endl;

        // Error attribution: compute energy with mixed parameters
        // E_fc_only: use Curcuma fc + Fortran alpha + Fortran r0
        // E_r0_only: use Fortran fc + Fortran alpha + Curcuma r0
        double E_fc_attribution = 0.0;
        double E_r0_attribution = 0.0;
        double E_alpha_attribution = 0.0;

        for (int rb = 0; rb < (int)ref_bonds.size(); ++rb) {
            int ref_i = ref_bonds[rb]["atoms"][0];
            int ref_j = ref_bonds[rb]["atoms"][1];
            double ref_r0 = (double)ref_bonds[rb]["r0"] * ANG_TO_BOHR; // Convert Å→Bohr
            double ref_alpha = ref_bonds[rb]["kbond"];
            double ref_fc = ref_bonds[rb]["prefactor"];

            auto key = std::make_pair(std::min(ref_i, ref_j), std::max(ref_i, ref_j));
            auto it = curcuma_bond_map.find(key);
            if (it == curcuma_bond_map.end()) continue;

            int cb = it->second;
            double cur_fc = curcuma_bonds[cb]["fc"];
            double cur_alpha = curcuma_bonds[cb]["exponent"];
            double cur_r0 = curcuma_bonds[cb]["r0_ij"];

            int ai = curcuma_bonds[cb]["i"];
            int aj = curcuma_bonds[cb]["j"];
            double dx = (mol.m_geometry(ai, 0) - mol.m_geometry(aj, 0)) * ANG_TO_BOHR;
            double dy = (mol.m_geometry(ai, 1) - mol.m_geometry(aj, 1)) * ANG_TO_BOHR;
            double dz = (mol.m_geometry(ai, 2) - mol.m_geometry(aj, 2)) * ANG_TO_BOHR;
            double rij = std::sqrt(dx*dx + dy*dy + dz*dz); // Bohr

            double dr_ref = rij - ref_r0;
            double dr_cur = rij - cur_r0;
            double E_ref = ref_fc * std::exp(-ref_alpha * dr_ref * dr_ref);

            // fc attribution: change only fc, keep Fortran alpha + r0
            double E_fc_only = cur_fc * std::exp(-ref_alpha * dr_ref * dr_ref);
            E_fc_attribution += (E_fc_only - E_ref);

            // r0 attribution: change only r0, keep Fortran fc + alpha
            double E_r0_only = ref_fc * std::exp(-ref_alpha * dr_cur * dr_cur);
            E_r0_attribution += (E_r0_only - E_ref);

            // alpha attribution: change only alpha, keep Fortran fc + r0
            double E_alpha_only = ref_fc * std::exp(-cur_alpha * dr_ref * dr_ref);
            E_alpha_attribution += (E_alpha_only - E_ref);
        }

        std::cout << std::endl << "Error Attribution (mEh):" << std::endl;
        std::cout << "  fc (prefactor) contribution:   " << std::fixed << std::setprecision(4) << E_fc_attribution * 1000.0 << " mEh" << std::endl;
        std::cout << "  alpha (exponent) contribution: " << E_alpha_attribution * 1000.0 << " mEh" << std::endl;
        std::cout << "  r0 (eq. dist.) contribution:   " << E_r0_attribution * 1000.0 << " mEh" << std::endl;
        std::cout << "  Total attributed:              " << (E_fc_attribution + E_alpha_attribution + E_r0_attribution) * 1000.0 << " mEh" << std::endl;
        std::cout << "  Actual total error:            " << (total_energy_curcuma - total_energy_ref) * 1000.0 << " mEh" << std::endl;
        std::cout << "  Cross-terms:                   "
                  << ((total_energy_curcuma - total_energy_ref) - E_fc_attribution - E_alpha_attribution - E_r0_attribution) * 1000.0
                  << " mEh" << std::endl;

        // === CORRECTED COMPARISON ===
        // The reference JSON r0 uses INIT formula: (ra+rb)*ff + rabshift
        // Fortran RUNTIME uses: (ra+rb+rabshift)*ff (same as Curcuma)
        // Correction: Fortran_runtime_r0 = ref_r0_Bohr + rabshift*(ff-1)
        std::cout << std::endl << "=== Corrected Runtime r0 Comparison ===" << std::endl;
        std::cout << "Init formula: r0 = (ra+rb)*ff + rabshift  [reference JSON]" << std::endl;
        std::cout << "Runtime formula: r0 = (ra+rb+rabshift)*ff  [both Curcuma & Fortran]" << std::endl;
        std::cout << "Correction: r0_runtime = r0_init + rabshift*(ff-1)" << std::endl;
        std::cout << std::endl;

        double total_E_corrected_ref = 0.0;
        double total_E_curcuma_dynamic = 0.0;
        int bonds_with_runtime_r0_diff = 0;

        std::cout << std::setw(4) << "Bond" << std::setw(5) << "Type"
                  << std::setw(12) << "r0_cur" << std::setw(12) << "r0_f_rt"
                  << std::setw(12) << "r0_diff" << std::setw(10) << "rij"
                  << std::setw(13) << "E_cur" << std::setw(13) << "E_f_rt"
                  << std::setw(13) << "dE" << std::endl;
        std::cout << std::string(105, '-') << std::endl;

        for (int rb = 0; rb < (int)ref_bonds.size(); ++rb) {
            int ref_i = ref_bonds[rb]["atoms"][0];
            int ref_j = ref_bonds[rb]["atoms"][1];
            double ref_r0_init_Bohr = (double)ref_bonds[rb]["r0"] * ANG_TO_BOHR;
            double ref_alpha = ref_bonds[rb]["kbond"];
            double ref_fc = ref_bonds[rb]["prefactor"];

            auto key = std::make_pair(std::min(ref_i, ref_j), std::max(ref_i, ref_j));
            auto it = curcuma_bond_map.find(key);
            if (it == curcuma_bond_map.end()) continue;

            int cb = it->second;
            double cur_r0 = curcuma_bonds[cb]["r0_ij"];
            double cur_fc = curcuma_bonds[cb]["fc"];
            double cur_alpha = curcuma_bonds[cb]["exponent"];
            double rabshift = curcuma_bonds[cb]["rabshift"];
            double ff = curcuma_bonds[cb]["ff"];

            // Compute corrected Fortran runtime r0
            double ref_r0_runtime = ref_r0_init_Bohr + rabshift * (ff - 1.0);

            int ai = curcuma_bonds[cb]["i"];
            int aj = curcuma_bonds[cb]["j"];
            double dx = (mol.m_geometry(ai, 0) - mol.m_geometry(aj, 0)) * ANG_TO_BOHR;
            double dy = (mol.m_geometry(ai, 1) - mol.m_geometry(aj, 1)) * ANG_TO_BOHR;
            double dz = (mol.m_geometry(ai, 2) - mol.m_geometry(aj, 2)) * ANG_TO_BOHR;
            double rij = std::sqrt(dx*dx + dy*dy + dz*dz);

            double dr_cur = rij - cur_r0;
            double E_cur = cur_fc * std::exp(-cur_alpha * dr_cur * dr_cur);

            double dr_ref_rt = rij - ref_r0_runtime;
            double E_ref_rt = ref_fc * std::exp(-ref_alpha * dr_ref_rt * dr_ref_rt);

            total_E_curcuma_dynamic += E_cur;
            total_E_corrected_ref += E_ref_rt;

            double r0_diff = cur_r0 - ref_r0_runtime;
            if (std::abs(r0_diff) > 1e-6) bonds_with_runtime_r0_diff++;

            // Determine bond type
            int z_i = curcuma_bonds[cb]["z_i"];
            int z_j = curcuma_bonds[cb]["z_j"];
            std::string btype;
            if (z_i == 1 || z_j == 1) {
                int heavy_z = (z_i == 1) ? z_j : z_i;
                if (heavy_z == 6) btype = "H-C";
                else if (heavy_z == 8) btype = "H-O";
                else btype = "H-?";
            } else if (z_i == 6 && z_j == 6) btype = "C-C";
            else if ((z_i == 6 && z_j == 8) || (z_i == 8 && z_j == 6)) btype = "C-O";
            else btype = "?-?";

            double dE = E_cur - E_ref_rt;
            bool show = std::abs(dE) > 1e-6 || std::abs(r0_diff) > 1e-5;
            if (show || ref_bonds.size() <= 80) {
                std::cout << std::setw(4) << rb << std::setw(5) << btype
                          << std::fixed << std::setprecision(6)
                          << std::setw(12) << cur_r0
                          << std::setw(12) << ref_r0_runtime
                          << std::scientific << std::setprecision(3)
                          << std::setw(12) << r0_diff
                          << std::fixed << std::setprecision(4)
                          << std::setw(10) << rij
                          << std::fixed << std::setprecision(9)
                          << std::setw(13) << E_cur
                          << std::setw(13) << E_ref_rt
                          << std::scientific << std::setprecision(4)
                          << std::setw(13) << dE
                          << std::endl;
            }
        }

        std::cout << std::string(105, '-') << std::endl;
        std::cout << std::fixed << std::setprecision(9)
                  << "Total E (Curcuma):          " << total_E_curcuma_dynamic << " Eh" << std::endl
                  << "Total E (Fortran runtime):  " << total_E_corrected_ref << " Eh" << std::endl
                  << "Difference:                 " << std::scientific << std::setprecision(6)
                  << (total_E_curcuma_dynamic - total_E_corrected_ref) << " Eh ("
                  << std::fixed << std::setprecision(3) << (total_E_curcuma_dynamic - total_E_corrected_ref) * 1000.0
                  << " mEh)" << std::endl;
        std::cout << "Actual bond energy error:   "
                  << std::fixed << std::setprecision(3)
                  << (m_gfnff->BondEnergy() - (double)m_ref_data["energy_components"]["bond"]) * 1000.0
                  << " mEh" << std::endl;
        std::cout << "Bonds with runtime r0 diff: " << bonds_with_runtime_r0_diff << std::endl;
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

        // Sum only the named components that match the Fortran SUM scope
        // (excludes nonbonded repulsion, sTorsions, XB, BATM, ATM, CN chain-rule)
        int natoms = m_gfnff->GradientBond().rows();
        Matrix grad = Matrix::Zero(natoms, 3);
        grad += m_gfnff->GradientBond();
        grad += m_gfnff->GradientAngle();
        grad += m_gfnff->GradientTorsion();
        grad += m_gfnff->GradientRepulsion();
        grad += m_gfnff->GradientCoulomb();
        grad += m_gfnff->GradientDispersion();
        grad += m_gfnff->GradientHB();

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

            // Accumulate Frobenius norms: ||g_curcuma|| and ||g_ref|| for this component
            double curcuma_norm_sq = curcuma_grad.squaredNorm();
            double ref_norm_sq = 0.0;

            for (auto it = ref_decomp.begin(); it != ref_decomp.end(); ++it) {
                int atom_idx = std::stoi(it.key());
                if (atom_idx >= natoms) continue;
                const auto& comps = it.value();

                if (!comps.contains(comp.fortran_name)) continue;

                double rx = comps[comp.fortran_name]["x"];
                double ry = comps[comp.fortran_name]["y"];
                double rz = comps[comp.fortran_name]["z"];
                ref_norm_sq += rx*rx + ry*ry + rz*rz;

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

            double curcuma_norm = std::sqrt(curcuma_norm_sq);
            double ref_norm = std::sqrt(ref_norm_sq);
            double norm_err = std::abs(curcuma_norm - ref_norm);

            bool passed = max_err < comp.tolerance;
            addResult("GradComp", comp.display_name,
                     passed, max_err, 0.0, max_err);

            // Claude Generated (Mar 2026): Also validate component gradient norm
            bool norm_passed = norm_err < comp.tolerance;
            addResult("GradNorm", comp.display_name,
                     norm_passed, curcuma_norm, ref_norm, norm_err);

            std::cout << "  " << std::left << std::setw(12) << comp.display_name
                      << " max_err=" << std::scientific << std::setprecision(4) << max_err
                      << " (atom " << max_err_atom << ")"
                      << " |g|=" << curcuma_norm
                      << " |g_ref|=" << ref_norm
                      << " d|g|=" << norm_err
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

    /**
     * @brief Charge injection diagnostic: compare native vs Fortran-reference charges
     *
     * Claude Generated (March 2026): Isolates charge-attributable errors per energy component.
     * 1. Run Calculation() with native EEQ charges → save component energies
     * 2. Inject Fortran reference charges, skip EEQ recalc → run Calculation() again
     * 3. Print comparison table: native error, injected error, charge attribution
     */
    void validateChargeInjection() {
        if (!m_ref_data.contains("energy_charges") || m_ref_data["energy_charges"].empty()) {
            std::cout << "  [SKIP] Charge injection: no energy_charges in reference JSON" << std::endl;
            return;
        }

        std::cout << std::endl << "=== Charge Injection Diagnostic ===" << std::endl;

        // Step 1: Normal calculation with native EEQ charges
        m_gfnff->Calculation(false);

        struct ComponentEnergy {
            std::string name;
            double native;
            double injected;
            double reference;
        };

        const auto& ref_e = m_ref_data["energy_components"];

        // Collect Fortran reference energies
        double ref_bond = ref_e["bond"];
        double ref_angle = ref_e["angle"];
        double ref_torsion = ref_e["torsion"];
        double ref_repulsion = ref_e["repulsion"];
        double ref_coulomb = ref_e["electrostatic"];
        double ref_dispersion = ref_e["dispersion"];
        std::string hb_key = ref_e.contains("hbond") ? "hbond" : (ref_e.contains("hb") ? "hb" : "");
        double ref_hb = hb_key.empty() ? 0.0 : (double)ref_e[hb_key];
        double ref_batm = ref_e.contains("batm") ? (double)ref_e["batm"] : 0.0;

        // Collect native energies
        double torsion_native = m_gfnff->DihedralEnergy() + m_gfnff->InversionEnergy();
        std::vector<ComponentEnergy> components = {
            {"Bond",       m_gfnff->BondEnergy(),       0.0, ref_bond},
            {"Angle",      m_gfnff->AngleEnergy(),      0.0, ref_angle},
            {"Torsion",    torsion_native,               0.0, ref_torsion},
            {"Repulsion",  m_gfnff->RepulsionEnergy(),  0.0, ref_repulsion},
            {"Coulomb",    m_gfnff->CoulombEnergy(),     0.0, ref_coulomb},
            {"Dispersion", m_gfnff->DispersionEnergy(), 0.0, ref_dispersion},
            {"HBond",      m_gfnff->HydrogenBondEnergy(), 0.0, ref_hb},
            {"BATM",       m_gfnff->BatmEnergy(),       0.0, ref_batm},
        };

        // Step 2: Inject Fortran reference charges
        const auto& ref_charges_json = m_ref_data["energy_charges"];
        Eigen::VectorXd ref_charges(ref_charges_json.size());
        for (int i = 0; i < (int)ref_charges_json.size(); ++i) {
            ref_charges(i) = ref_charges_json[i];
        }

        m_gfnff->setCharges(ref_charges);
        m_gfnff->setSkipEEQRecalc(true);

        // Step 3: Recalculate with injected charges
        m_gfnff->Calculation(false);

        double torsion_injected = m_gfnff->DihedralEnergy() + m_gfnff->InversionEnergy();
        components[0].injected = m_gfnff->BondEnergy();
        components[1].injected = m_gfnff->AngleEnergy();
        components[2].injected = torsion_injected;
        components[3].injected = m_gfnff->RepulsionEnergy();
        components[4].injected = m_gfnff->CoulombEnergy();
        components[5].injected = m_gfnff->DispersionEnergy();
        components[6].injected = m_gfnff->HydrogenBondEnergy();
        components[7].injected = m_gfnff->BatmEnergy();

        // Restore normal operation
        m_gfnff->setSkipEEQRecalc(false);

        // Step 4: Print comparison table
        std::cout << std::left
                  << std::setw(12) << "Component"
                  << std::right
                  << std::setw(14) << "Native(Eh)"
                  << std::setw(14) << "Injected(Eh)"
                  << std::setw(14) << "Fortran(Eh)"
                  << std::setw(16) << "Err_Nat(mEh)"
                  << std::setw(16) << "Err_Inj(mEh)"
                  << std::setw(16) << "ChgAttr(mEh)"
                  << std::endl;
        std::cout << std::string(102, '-') << std::endl;

        double total_native = 0, total_injected = 0, total_ref = 0;
        for (const auto& c : components) {
            double err_native = (c.native - c.reference) * 1000.0;
            double err_injected = (c.injected - c.reference) * 1000.0;
            double charge_attr = err_native - err_injected;

            std::cout << std::left << std::setw(12) << c.name
                      << std::right << std::fixed
                      << std::setprecision(6) << std::setw(14) << c.native
                      << std::setw(14) << c.injected
                      << std::setw(14) << c.reference
                      << std::setprecision(3) << std::setw(16) << err_native
                      << std::setw(16) << err_injected
                      << std::setw(16) << charge_attr
                      << std::endl;

            total_native += c.native;
            total_injected += c.injected;
            total_ref += c.reference;
        }

        std::cout << std::string(102, '-') << std::endl;
        double err_total_nat = (total_native - total_ref) * 1000.0;
        double err_total_inj = (total_injected - total_ref) * 1000.0;
        std::cout << std::left << std::setw(12) << "TOTAL"
                  << std::right << std::fixed
                  << std::setprecision(6) << std::setw(14) << total_native
                  << std::setw(14) << total_injected
                  << std::setw(14) << total_ref
                  << std::setprecision(3) << std::setw(16) << err_total_nat
                  << std::setw(16) << err_total_inj
                  << std::setw(16) << (err_total_nat - err_total_inj)
                  << std::endl;
        std::cout << std::endl;
        std::cout << "ChgAttr = Err_Native - Err_Injected (positive = error caused by charge differences)" << std::endl;
    }

    /**
     * @brief Per-torsion diagnostic: compare C++ vs Fortran torsion parameters
     *
     * Claude Generated (March 2026): Pinpoints which torsions have wrong V/n/phi0
     * or are missing/extra between C++ and Fortran implementations.
     */
    void validateTorsionDetail() {
        std::cout << std::endl << "=== Torsion Detail Diagnostic ===" << std::endl;

        // Run calculation to populate parameters
        m_gfnff->Calculation(false);

        // Get C++ torsion data
        json cpp_torsions = m_gfnff->getTorsionParameters();
        json cpp_primary = cpp_torsions.value("primary", json::array());
        json cpp_extra = cpp_torsions.value("extra", json::array());
        json cpp_inversions = m_gfnff->getInversionParameters();

        // Reference has all torsions (primary+extra) in "torsions" array
        json ref_all = m_ref_data.value("torsions", json::array());

        std::cout << "  C++ primary: " << cpp_primary.size()
                  << "  C++ extra: " << cpp_extra.size()
                  << "  C++ total: " << (cpp_primary.size() + cpp_extra.size())
                  << "  Ref total: " << ref_all.size()
                  << "  C++ inversions: " << cpp_inversions.size() << std::endl;

        // Canonical key: sort forward/reverse to match regardless of direction
        auto canonical_key = [](int a, int b, int c, int d) -> std::array<int,4> {
            std::array<int,4> fwd = {a, b, c, d};
            std::array<int,4> rev = {d, c, b, a};
            return (fwd < rev) ? fwd : rev;
        };

        auto key_str = [](const std::array<int,4>& k) -> std::string {
            return std::to_string(k[0]) + "-" + std::to_string(k[1]) + "-"
                 + std::to_string(k[2]) + "-" + std::to_string(k[3]);
        };

        struct TorsionInfo {
            double V = 0, phi0 = 0, energy = 0;
            int n = 0, idx = -1;
            std::array<int,4> atoms;
            bool is_extra = false;
        };

        // Build combined C++ map (primary + extra)
        std::map<std::array<int,4>, TorsionInfo> cpp_map;
        for (int idx = 0; idx < (int)cpp_primary.size(); ++idx) {
            const auto& t = cpp_primary[idx];
            int i = t["i"].get<int>(), j = t["j"].get<int>();
            int k = t["k"].get<int>(), l = t["l"].get<int>();
            auto key = canonical_key(i, j, k, l);
            TorsionInfo info;
            info.V = t.value("V", 0.0);
            info.n = t.value("n", 0);
            info.phi0 = t.value("phi0", 0.0);
            info.idx = idx;
            info.atoms = {i, j, k, l};
            info.is_extra = false;
            cpp_map[key] = info;
        }
        for (int idx = 0; idx < (int)cpp_extra.size(); ++idx) {
            const auto& t = cpp_extra[idx];
            int i = t["i"].get<int>(), j = t["j"].get<int>();
            int k = t["k"].get<int>(), l = t["l"].get<int>();
            auto key = canonical_key(i, j, k, l);
            TorsionInfo info;
            info.V = t.value("V", 0.0);
            info.n = t.value("n", 0);
            info.phi0 = t.value("phi0", 0.0);
            info.idx = idx;
            info.atoms = {i, j, k, l};
            info.is_extra = true;
            cpp_map[key] = info;
        }

        // Build reference map
        std::map<std::array<int,4>, TorsionInfo> ref_map;
        for (int idx = 0; idx < (int)ref_all.size(); ++idx) {
            const auto& t = ref_all[idx];
            auto atoms_arr = t["atoms"];
            int i = atoms_arr[0].get<int>(), j = atoms_arr[1].get<int>();
            int k = atoms_arr[2].get<int>(), l = atoms_arr[3].get<int>();
            auto key = canonical_key(i, j, k, l);
            TorsionInfo info;
            info.V = t.value("V", 0.0);
            info.n = t.value("n", 0);
            info.phi0 = t.value("phi0", 0.0);
            info.energy = t.value("energy", 0.0);
            info.idx = idx;
            info.atoms = {i, j, k, l};
            ref_map[key] = info;
        }

        // Compare
        std::cout << std::endl << std::left << std::setw(22) << "Atoms"
                  << std::setw(5) << "Type"
                  << std::right
                  << std::setw(4) << "n_C" << std::setw(4) << "n_R"
                  << std::setw(14) << "V_C++"
                  << std::setw(14) << "V_Ref"
                  << std::setw(12) << "dV(mEh)"
                  << std::setw(10) << "phi0_C"
                  << std::setw(10) << "phi0_R"
                  << std::endl;
        std::cout << std::string(95, '-') << std::endl;

        int matched = 0, n_mismatch = 0, phi0_mismatch = 0;
        double sum_dV = 0.0, max_dV = 0.0;
        std::string max_dV_key;

        for (const auto& [key, cpp_info] : cpp_map) {
            auto it = ref_map.find(key);
            if (it != ref_map.end()) {
                matched++;
                const auto& ref_info = it->second;
                double dV = (cpp_info.V - ref_info.V) * 1000.0;
                sum_dV += dV;
                if (std::abs(dV) > std::abs(max_dV)) {
                    max_dV = dV;
                    max_dV_key = key_str(key);
                }
                bool n_diff = (cpp_info.n != ref_info.n);
                bool phi0_diff = (std::abs(cpp_info.phi0 - ref_info.phi0) > 0.01);
                if (n_diff) n_mismatch++;
                if (phi0_diff) phi0_mismatch++;

                // Print entries with significant V difference, or n/phi0 mismatch
                if (std::abs(dV) > 0.001 || n_diff || phi0_diff) {
                    std::cout << std::left << std::setw(22) << key_str(key)
                              << std::setw(5) << (cpp_info.is_extra ? "extra" : "prim")
                              << std::right
                              << std::setw(4) << cpp_info.n
                              << std::setw(4) << ref_info.n
                              << std::fixed << std::setprecision(6)
                              << std::setw(14) << cpp_info.V
                              << std::setw(14) << ref_info.V
                              << std::setprecision(3) << std::setw(12) << dV
                              << std::setprecision(1) << std::setw(10) << (cpp_info.phi0 * 180.0 / M_PI)
                              << std::setw(10) << (ref_info.phi0 * 180.0 / M_PI)
                              << std::endl;
                }
            }
        }

        // Unmatched
        std::vector<std::string> unmatched_cpp, unmatched_ref;
        for (const auto& [key, info] : cpp_map) {
            if (ref_map.find(key) == ref_map.end()) {
                unmatched_cpp.push_back(key_str(key) + " n=" + std::to_string(info.n)
                    + " V=" + std::to_string(info.V) + (info.is_extra ? " [extra]" : " [prim]"));
            }
        }
        for (const auto& [key, info] : ref_map) {
            if (cpp_map.find(key) == cpp_map.end()) {
                unmatched_ref.push_back(key_str(key) + " n=" + std::to_string(info.n)
                    + " V=" + std::to_string(info.V));
            }
        }

        if (!unmatched_cpp.empty()) {
            std::cout << std::endl << "  UNMATCHED C++ (" << unmatched_cpp.size() << " not in reference):" << std::endl;
            for (const auto& s : unmatched_cpp) std::cout << "    " << s << std::endl;
        }
        if (!unmatched_ref.empty()) {
            std::cout << std::endl << "  UNMATCHED REF (" << unmatched_ref.size() << " not in C++):" << std::endl;
            for (const auto& s : unmatched_ref) std::cout << "    " << s << std::endl;
        }

        std::cout << std::endl << "  Matched=" << matched
                  << "  n_mismatch=" << n_mismatch
                  << "  phi0_mismatch=" << phi0_mismatch
                  << "  Unmatched_C++=" << unmatched_cpp.size()
                  << "  Unmatched_Ref=" << unmatched_ref.size() << std::endl;
        std::cout << "  Sum(dV)=" << std::fixed << std::setprecision(3) << sum_dV
                  << " mEh  Max(|dV|)=" << std::abs(max_dV) << " mEh at " << max_dV_key << std::endl;

        // Inversion summary
        std::cout << std::endl << "=== INVERSIONS ===" << std::endl;
        std::cout << "  C++ inversions: " << cpp_inversions.size() << std::endl;
        for (int idx = 0; idx < (int)cpp_inversions.size(); ++idx) {
            const auto& inv = cpp_inversions[idx];
            std::cout << "    " << inv["i"].get<int>() << "-" << inv["j"].get<int>()
                      << "-" << inv["k"].get<int>() << "-" << inv["l"].get<int>()
                      << "  fc=" << std::fixed << std::setprecision(6) << inv.value("fc", 0.0)
                      << "  type=" << inv.value("potential_type", 0)
                      << "  omega0=" << std::setprecision(4) << inv.value("omega0", 0.0)
                      << std::endl;
        }

        // Energy summary
        double dihedral_energy = m_gfnff->DihedralEnergy();
        double inversion_energy = m_gfnff->InversionEnergy();
        double combined = dihedral_energy + inversion_energy;
        double ref_torsion = m_ref_data.contains("energy_components")
                           && m_ref_data["energy_components"].contains("torsion")
                           ? (double)m_ref_data["energy_components"]["torsion"] : 0.0;
        double diff = (combined - ref_torsion) * 1000.0;

        std::cout << std::endl << "=== ENERGY SUMMARY ===" << std::endl;
        std::cout << "  DihedralEnergy  = " << std::fixed << std::setprecision(8) << dihedral_energy << " Eh" << std::endl;
        std::cout << "  InversionEnergy = " << inversion_energy << " Eh" << std::endl;
        std::cout << "  Combined        = " << combined << " Eh" << std::endl;
        std::cout << "  Ref torsion     = " << ref_torsion << " Eh" << std::endl;
        std::cout << "  Diff            = " << std::setprecision(3) << diff << " mEh" << std::endl;
    }

    /**
     * @brief Dispersion gradient diagnostic: split direct vs CN chain-rule
     *
     * Claude Generated (March 2026): Splits dispersion gradient into:
     *   1. Direct pairwise component (from thread m_gradient_dispersion)
     *   2. CN chain-rule correction (dEdcn_disp * dcn)
     * Reports per-atom norms and compares against Fortran reference.
     */
    void validateDispersionDiagnostic() {
        std::cout << std::endl << "=== Dispersion Gradient Diagnostic ===" << std::endl;

        // Get total dispersion gradient (direct + CN chain-rule)
        Matrix grad_disp_total = m_gfnff->GradientDispersion();
        Matrix grad_disp_cn = m_gfnff->getDispCNCorrection();

        if (grad_disp_total.rows() == 0) {
            std::cout << "  [SKIP] No dispersion gradient available" << std::endl;
            return;
        }

        int natoms = grad_disp_total.rows();

        // Direct = Total - CN correction
        Matrix grad_disp_direct = grad_disp_total;
        bool has_cn_correction = (grad_disp_cn.rows() == natoms);
        if (has_cn_correction) {
            grad_disp_direct = grad_disp_total - grad_disp_cn;
        }

        // Compare against Fortran reference if available
        const auto& ref_decomp = m_ref_data["gradient_decomposition"];

        double max_total_err = 0.0, max_direct_err = 0.0, max_cn_norm = 0.0;
        int max_total_atom = -1, max_direct_atom = -1, max_cn_atom = -1;

        std::cout << std::left << std::setw(8) << "Atom"
                  << std::right
                  << std::setw(14) << "|Direct|"
                  << std::setw(14) << "|CN_corr|"
                  << std::setw(14) << "|Total|"
                  << std::setw(14) << "|Ref_Disp|"
                  << std::setw(14) << "|Err_Total|"
                  << std::setw(14) << "|Err_Direct|"
                  << std::endl;
        std::cout << std::string(92, '-') << std::endl;

        for (int i = 0; i < natoms; ++i) {
            double direct_norm = grad_disp_direct.row(i).norm();
            double cn_norm = has_cn_correction ? grad_disp_cn.row(i).norm() : 0.0;
            double total_norm = grad_disp_total.row(i).norm();

            if (cn_norm > max_cn_norm) { max_cn_norm = cn_norm; max_cn_atom = i; }

            // Compare against Fortran Disp reference
            std::string atom_key = std::to_string(i);
            double ref_norm = 0.0;
            double total_err = 0.0, direct_err = 0.0;

            if (ref_decomp.contains(atom_key) && ref_decomp[atom_key].contains("Disp")) {
                double rx = ref_decomp[atom_key]["Disp"]["x"];
                double ry = ref_decomp[atom_key]["Disp"]["y"];
                double rz = ref_decomp[atom_key]["Disp"]["z"];
                ref_norm = std::sqrt(rx*rx + ry*ry + rz*rz);

                double tx = grad_disp_total(i, 0), ty = grad_disp_total(i, 1), tz = grad_disp_total(i, 2);
                total_err = std::sqrt((tx-rx)*(tx-rx) + (ty-ry)*(ty-ry) + (tz-rz)*(tz-rz));

                double dx = grad_disp_direct(i, 0), dy = grad_disp_direct(i, 1), dz = grad_disp_direct(i, 2);
                direct_err = std::sqrt((dx-rx)*(dx-rx) + (dy-ry)*(dy-ry) + (dz-rz)*(dz-rz));

                if (total_err > max_total_err) { max_total_err = total_err; max_total_atom = i; }
                if (direct_err > max_direct_err) { max_direct_err = direct_err; max_direct_atom = i; }
            }

            // Print atoms with significant gradients or errors
            if (total_norm > 1e-5 || total_err > 1e-5 || cn_norm > 1e-6) {
                std::string elem = "";
                if (m_ref_data.contains("molecule") && m_ref_data["molecule"].contains("atoms")
                    && i < (int)m_ref_data["molecule"]["atoms"].size()) {
                    elem = m_ref_data["molecule"]["atoms"][i]["element"].get<std::string>();
                }
                std::cout << std::left << std::setw(4) << i << std::setw(4) << elem
                          << std::right << std::scientific << std::setprecision(4)
                          << std::setw(14) << direct_norm
                          << std::setw(14) << cn_norm
                          << std::setw(14) << total_norm
                          << std::setw(14) << ref_norm
                          << std::setw(14) << total_err
                          << std::setw(14) << direct_err
                          << std::endl;
            }
        }

        std::cout << std::endl;
        std::cout << "  Max |CN correction|: " << std::scientific << std::setprecision(4)
                  << max_cn_norm << " (atom " << max_cn_atom << ")" << std::endl;
        std::cout << "  Max |Total err|:     " << max_total_err << " (atom " << max_total_atom << ")" << std::endl;
        std::cout << "  Max |Direct-only err|: " << max_direct_err << " (atom " << max_direct_atom << ")" << std::endl;
        std::cout << "  CN correction " << (has_cn_correction ? "PRESENT" : "MISSING (all zero)")
                  << std::endl;
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
    bool m_rep_diag = false;
    std::string m_solve_method;  ///< EEQ solve method override (empty = default)
};

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " [--fd] [--charge-inject] [--torsion-detail] [--disp-diag] [--solve-method lu|schur_cholesky|pcg] <reference_json> [reference_json2 ...]" << std::endl;
        std::cerr << "  --fd              Enable finite-difference gradient validation (slow)" << std::endl;
        std::cerr << "  --charge-inject   Inject Fortran reference charges to isolate charge-attributable errors" << std::endl;
        std::cerr << "  --torsion-detail  Per-torsion parameter comparison between C++ and Fortran" << std::endl;
        std::cerr << "  --disp-diag       Split dispersion gradient into direct vs CN chain-rule parts" << std::endl;
        std::cerr << "  --solve-method    EEQ linear solve: lu, schur_cholesky (default), pcg" << std::endl;
        return 1;
    }

    // Two-pass argument parsing: flags first, then files
    bool run_fd = false;
    bool run_charge_inject = false;
    bool run_torsion_detail = false;
    bool run_rep_diag = false;
    bool run_disp_diag = false;
    bool all_files_passed = true;
    std::string solve_method;  // Empty = use default (schur_cholesky)
    std::vector<std::string> files;

    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
        if (arg == "--fd") run_fd = true;
        else if (arg == "--charge-inject") run_charge_inject = true;
        else if (arg == "--torsion-detail") run_torsion_detail = true;
        else if (arg == "--rep-diag") run_rep_diag = true;
        else if (arg == "--disp-diag") run_disp_diag = true;
        else if (arg == "--solve-method" && i + 1 < argc) {
            solve_method = argv[++i];
            std::cout << "EEQ solve method: " << solve_method << std::endl;
        }
        else files.push_back(arg);
    }

    for (const auto& file : files) {
        try {
            GFNFFValidator validator(file.c_str(), run_rep_diag, solve_method);
            validator.runValidation(run_fd, run_charge_inject, run_torsion_detail, run_disp_diag);
            if (!validator.allPassed()) all_files_passed = false;
        } catch (const std::exception& e) {
            std::cerr << "Error processing " << file << ": " << e.what() << std::endl;
            all_files_passed = false;
        }
        std::cout << std::endl;
    }

    return all_files_passed ? 0 : 1;
}
