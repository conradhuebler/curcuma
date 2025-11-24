/*
 * <GFN-FF Implementation for Curcuma>
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "gfnff.h"
#include "src/core/curcuma_logger.h"
#include "src/core/units.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <queue>
#include <stack>
#include <string>

#include <fmt/format.h>

GFNFF::GFNFF()
    : m_forcefield(nullptr)
    , m_initialized(false)
    , m_energy_total(0.0)
{
    json default_parameters = {
        { "method", "gfnff" },
        { "threads", 1 },
        { "gradient", 1 },
        { "dispersion", true },
        { "hbond", true },
        { "repulsion_scaling", 1.0 }
    };
    m_parameters = default_parameters;
}

GFNFF::GFNFF(const json& parameters)
    : m_forcefield(nullptr)
    , m_initialized(false)
    , m_energy_total(0.0)
{
    json default_parameters = {
        { "method", "gfnff" },
        { "threads", 1 },
        { "gradient", 1 },
        { "dispersion", true },
        { "hbond", true },
        { "repulsion_scaling", 1.0 }
    };

    m_parameters = MergeJson(default_parameters, parameters);
}

GFNFF::~GFNFF()
{
    if (m_forcefield) {
        delete m_forcefield;
        m_forcefield = nullptr;
    }
}

bool GFNFF::InitialiseMolecule()
{
    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("=== GFNFF::InitialiseMolecule() START ===");
        CurcumaLogger::param("atom_count", std::to_string(m_atomcount));
        CurcumaLogger::param("geometry_rows", std::to_string(m_geometry_bohr.rows()));
        CurcumaLogger::param("geometry_cols", std::to_string(m_geometry.cols()));
    }

    if (m_atomcount == 0) {
        CurcumaLogger::error("GFN-FF initialization failed: No atoms in molecule");
        return false;
    }

    // Convert geometry from Angström to Bohr (GFN-FF parameters are in Bohr)
    // TODO: Alternative approach (Option B) would convert parameters to Angström instead.
    // See documentation block at parameter arrays for details on parameter conversion.
    m_geometry_bohr = m_geometry * CurcumaUnit::Length::ANGSTROM_TO_BOHR;

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("Validating molecule structure...");
    }

    if (!validateMolecule()) {
        CurcumaLogger::error("GFN-FF initialization failed: Molecule validation failed");
        return false;
    }

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::success("Molecule validation passed");
        CurcumaLogger::info("Initializing force field...");
    }

    if (!initializeForceField()) {
        CurcumaLogger::error("GFN-FF initialization failed: Force field initialization failed");
        return false;
    }

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::success("Force field initialization successful");
        CurcumaLogger::info("Allocating gradient and charge arrays...");
    }

    m_gradient = Matrix::Zero(m_atomcount, 3);
    m_charges = Vector::Zero(m_atomcount);
    m_bond_orders = Vector::Zero(m_atomcount * (m_atomcount - 1) / 2);

    m_initialized = true;

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::success("GFN-FF initialization complete");
        CurcumaLogger::param("initialized", "true");
    }

    return true;
}

bool GFNFF::UpdateMolecule()
{
    if (!m_initialized) {
        return InitialiseMolecule();
    }

    // Update Bohr geometry for GFN-FF
    m_geometry_bohr = m_geometry * CurcumaUnit::Length::ANGSTROM_TO_BOHR;

    if (m_forcefield) {
        m_forcefield->UpdateGeometry(m_geometry_bohr);  // Pass Bohr geometry
    }

    return true;
}

double GFNFF::Calculation(bool gradient)
{
    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("=== GFNFF::Calculation() START ===");
        CurcumaLogger::param("gradient_requested", gradient ? "true" : "false");
        CurcumaLogger::param("initialized", m_initialized ? "true" : "false");
        CurcumaLogger::param("forcefield_ptr", m_forcefield ? "valid" : "null");
    }

    if (!m_initialized) {
        CurcumaLogger::error("GFN-FF calculation failed: Not initialized");
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::param("atom_count", std::to_string(m_atomcount));
            CurcumaLogger::param("geometry_size", std::to_string(m_geometry_bohr.rows()) + "x" + std::to_string(m_geometry.cols()));
        }
        return 0.0;
    }

    if (!m_forcefield) {
        CurcumaLogger::error("GFN-FF calculation failed: Force field not available");
        return 0.0;
    }

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("Calling ForceField::Calculate()...");
    }

    // ForceField returns energy in Hartree (m_final_factor = 1)
    double energy_hartree = m_forcefield->Calculate(gradient);

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::param("energy_hartree_raw", fmt::format("{:.8f}", energy_hartree));
    }

    if (gradient) {
        // ForceField returns gradient in Hartree/Bohr (m_final_factor = 1)
        Matrix grad_hartree = m_forcefield->Gradient();
        m_gradient = grad_hartree;  // No conversion needed

        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::param("gradient_norm", fmt::format("{:.8f}", m_gradient.norm()));
        }
    }

    // No unit conversion needed - already in Hartree
    m_energy_total = energy_hartree;

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::energy_abs(m_energy_total, "GFN-FF Energy");
        CurcumaLogger::param("energy_hartree", fmt::format("{:.10f}", energy_hartree));
    }

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::success("GFN-FF calculation complete");
        CurcumaLogger::param("energy_hartree", fmt::format("{:.10f}", m_energy_total));
    }

    return m_energy_total;
}

Vector GFNFF::Charges() const
{
    return m_charges;
}

Vector GFNFF::BondOrders() const
{
    return m_bond_orders;
}

void GFNFF::setParameters(const json& parameters)
{
    m_parameters = MergeJson(m_parameters, parameters);

    if (m_forcefield && m_initialized) {
        json ff_params = generateGFNFFParameters();
        m_forcefield->setParameter(ff_params);
    }
}

bool GFNFF::initializeForceField()
{
    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("=== GFNFF::initializeForceField() START ===");
        CurcumaLogger::param("forcefield_exists", m_forcefield ? "yes (will delete)" : "no");
    }

    if (m_forcefield) {
        delete m_forcefield;
    }

    json ff_config = {
        { "threads", m_parameters["threads"] },
        { "gradient", m_parameters["gradient"] },
        { "method", "gfnff" }
    };

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("Creating ForceField instance...");
        CurcumaLogger::param("threads", std::to_string(m_parameters.value("threads", 1)));
        CurcumaLogger::param("gradient", std::to_string(m_parameters.value("gradient", 1)));
    }

    m_forcefield = new ForceField(ff_config);
    m_forcefield->setAtomTypes(m_atoms);

    // CRITICAL FIX: Set geometry in ForceField! setAtomTypes() only sets atom types, not geometry
    // Without this, m_geometry in ForceField is empty, causing out-of-bounds access in threads
    // NOTE: Use Bohr geometry (m_geometry_bohr) because GFN-FF parameters are in Bohr
    m_forcefield->UpdateGeometry(m_geometry_bohr);

    // TEMPORARY DEBUG: Disable caching to isolate the problem
    m_forcefield->setParameterCaching(false);

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::success("ForceField instance created");
        CurcumaLogger::param("geometry_set", std::to_string(m_geometry_bohr.rows()) + " atoms");
        CurcumaLogger::warn("TEMPORARY: Parameter caching disabled for debugging");
        CurcumaLogger::info("Calculating topology (bonds, angles, torsions, inversions)...");
    }

    if (!calculateTopology()) {
        CurcumaLogger::error("Topology calculation failed");
        return false;
    }

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::success("Topology calculation complete");
        CurcumaLogger::info("Generating GFN-FF parameters...");
    }

    json ff_params = generateGFNFFParameters();

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::success("GFN-FF parameters generated");
        CurcumaLogger::param("bonds_count", std::to_string(ff_params.value("bonds", json::array()).size()));
        CurcumaLogger::param("angles_count", std::to_string(ff_params.value("angles", json::array()).size()));
        CurcumaLogger::param("torsions_count", std::to_string(ff_params.value("dihedrals", json::array()).size()));
        CurcumaLogger::param("inversions_count", std::to_string(ff_params.value("inversions", json::array()).size()));
        CurcumaLogger::info("Setting parameters in ForceField...");
    }

    m_forcefield->setParameter(ff_params);

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::success("ForceField initialization complete");
    }

    return true;
}

json GFNFF::generateGFNFFParameters()
{
    json parameters;
    parameters["method"] = "gfnff";
    parameters["e0"] = 0.0;

    // Check if advanced parametrization is enabled
    bool use_advanced = m_parameters.value("use_advanced_parametrization", false);

    if (use_advanced) {
        std::cout << "Using advanced GFN-FF parametrization (experimental)" << std::endl;

        // Calculate topology information for advanced parametrization
        TopologyInfo topo_info = calculateTopologyInfo();

        // Generate advanced parameters
        json bonds = generateTopologyAwareBonds(topo_info.coordination_numbers,
            topo_info.hybridization,
            topo_info.eeq_charges,
            topo_info.ring_sizes);
        json angles = generateTopologyAwareAngles(topo_info.coordination_numbers,
            topo_info.hybridization,
            topo_info.eeq_charges,
            topo_info.ring_sizes);

        parameters["bonds"] = bonds;
        parameters["angles"] = angles;
        parameters["dihedrals"] = generateGFNFFTorsions(); // ✅ Phase 1.1 implemented
        parameters["inversions"] = generateGFNFFInversions(); // ✅ Phase 1.2 implemented

        // Phase 4.2: Generate pairwise non-bonded parameters
        parameters["gfnff_coulombs"] = generateGFNFFCoulombPairs();
        parameters["gfnff_repulsions"] = generateGFNFFRepulsionPairs();
        parameters["gfnff_dispersions"] = generateGFNFFDispersionPairs();

        parameters["vdws"] = json::array(); // Legacy vdW (will be replaced by pairwise)
        parameters["hbonds"] = detectHydrogenBonds(topo_info.eeq_charges);

        // Store topology information for debugging
        parameters["topology_info"] = {
            { "coordination_numbers", std::vector<double>(topo_info.coordination_numbers.data(), topo_info.coordination_numbers.data() + topo_info.coordination_numbers.size()) },
            { "hybridization", topo_info.hybridization },
            { "ring_sizes", topo_info.ring_sizes },
            { "eeq_charges", std::vector<double>(topo_info.eeq_charges.data(), topo_info.eeq_charges.data() + topo_info.eeq_charges.size()) }
        };

        // Use calculated charges instead of loading from file
        m_charges = topo_info.eeq_charges;

    } else {
        std::cout << "Using basic GFN-FF parametrization" << std::endl;

        // Calculate topology information for Phase 2 angle parametrization
        TopologyInfo topo_info = calculateTopologyInfo();

        // Generate GFN-FF bonds with real parameters
        json bonds = generateGFNFFBonds();
        json angles = generateGFNFFAngles(topo_info);
        json torsions = generateGFNFFTorsions(); // ✅ Phase 1.1 implemented
        json inversions = generateGFNFFInversions(); // ✅ Phase 1.2 implemented

        parameters["bonds"] = bonds;
        parameters["angles"] = angles;
        parameters["dihedrals"] = torsions;
        parameters["inversions"] = inversions;

        // Store topology charges for use in other functions
        m_charges = topo_info.eeq_charges;

        // Phase 4.2: Generate pairwise non-bonded parameters
        parameters["gfnff_coulombs"] = generateGFNFFCoulombPairs();
        parameters["gfnff_repulsions"] = generateGFNFFRepulsionPairs();
        parameters["gfnff_dispersions"] = generateGFNFFDispersionPairs();

        parameters["vdws"] = json::array(); // Legacy vdW (will be replaced by pairwise)
    }

    // Integration with existing corrections (H4, D3/D4)
    if (m_parameters.contains("dispersion") && m_parameters["dispersion"]) {
        parameters["use_dispersion"] = true;
        parameters["use_d4"] = true; // Use existing D4 implementation
    }

    if (m_parameters.contains("hbond") && m_parameters["hbond"]) {
        parameters["use_hbond"] = true;
        parameters["use_h4"] = true; // Use existing H4 implementation
    }

    // GFN-FF specific settings
    parameters["repulsion_scaling"] = m_parameters.value("repulsion_scaling", 1.0);

    return parameters;
}

json GFNFF::generateGFNFFBonds() const
{
    json bonds = json::array();

    // Phase 9: Calculate topology information for complete bond parameters
    TopologyInfo topo_info = calculateTopologyInfo();

    // GFN-FF bond detection with connectivity threshold
    double bond_threshold = 1.3; // Factor for covalent radii sum

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info(fmt::format("GFN-FF bond detection: {} atoms, threshold {:.2f}",
                                         m_atomcount, bond_threshold));
    } else {
        std::cout << "GFN-FF bond detection: " << m_atomcount << " atoms, threshold " << bond_threshold << std::endl;
    }

    for (int i = 0; i < m_atomcount; ++i) {
        for (int j = i + 1; j < m_atomcount; ++j) {
            Vector ri = m_geometry_bohr.row(i);
            Vector rj = m_geometry_bohr.row(j);
            double distance = (ri - rj).norm();

            // Get covalent radii for atoms i and j
            double rcov_i = getCovalentRadius(m_atoms[i]);
            double rcov_j = getCovalentRadius(m_atoms[j]);

            if (distance < bond_threshold * (rcov_i + rcov_j)) {
                // This is a bond - generate GFN-FF parameters
                json bond;
                bond["type"] = 3; // GFN-FF type
                bond["i"] = i;
                bond["j"] = j;
                bond["k"] = 0; // Not used in GFN-FF but required by ForceField
                bond["distance"] = distance; // Current bond distance

                // Phase 9: GFN-FF bond parameters with full topology awareness
                auto bond_params = getGFNFFBondParameters(i, j, m_atoms[i], m_atoms[j], distance, topo_info);
                bond["fc"] = bond_params.force_constant;
                bond["r0_ij"] = bond_params.equilibrium_distance;
                bond["r0_ik"] = 0.0; // Not used in GFN-FF but required by ForceField
                bond["exponent"] = bond_params.alpha;  // Phase 1.3: store α in exponent field

                bonds.push_back(bond);
            }
        }
    }

    if (bonds.empty()) {
        CurcumaLogger::warn("No bonds detected in GFN-FF");
    } else {
        if (CurcumaLogger::get_verbosity() >= 1) {
            CurcumaLogger::success(fmt::format("GFN-FF detected {} bonds", bonds.size()));
        } else {
            std::cout << "GFN-FF detected " << bonds.size() << " bonds" << std::endl;
        }
    }

    return bonds;
}

json GFNFF::generateGFNFFAngles(const TopologyInfo& topo_info) const
{
    json angles = json::array();

    // First, collect all bonds for angle detection
    std::vector<std::pair<int, int>> bond_list;
    json bonds = generateGFNFFBonds();

    for (const auto& bond : bonds) {
        bond_list.push_back({ bond["i"], bond["j"] });
    }

    // Generate angles from bonded topology
    for (int center = 0; center < m_atomcount; ++center) {
        std::vector<int> neighbors;

        // Find all atoms bonded to center
        for (const auto& bond : bond_list) {
            if (bond.first == center)
                neighbors.push_back(bond.second);
            if (bond.second == center)
                neighbors.push_back(bond.first);
        }

        // Generate all possible angles with center as middle atom
        for (int i = 0; i < neighbors.size(); ++i) {
            for (int j = i + 1; j < neighbors.size(); ++j) {
                json angle;
                angle["type"] = 3; // GFN-FF type
                angle["i"] = neighbors[i];
                angle["j"] = center;
                angle["k"] = neighbors[j];

                // Calculate current angle for reference
                Vector ri = m_geometry_bohr.row(neighbors[i]);
                Vector rj = m_geometry_bohr.row(center);
                Vector rk = m_geometry_bohr.row(neighbors[j]);

                Vector v1 = ri - rj;
                Vector v2 = rk - rj;

                // Safely calculate angle with bounds checking
                double v1_norm = v1.norm();
                double v2_norm = v2.norm();

                // Skip if vectors are too small (linear geometry or duplicate atoms)
                if (v1_norm < 1e-10 || v2_norm < 1e-10) {
                    continue;
                }

                double cos_angle = v1.dot(v2) / (v1_norm * v2_norm);
                // Clamp to valid acos range [-1, 1] to avoid NaN
                cos_angle = std::max(-1.0, std::min(1.0, cos_angle));
                double current_angle = acos(cos_angle);

                // GFN-FF angle parameters (Claude Generated Nov 2025: Phase 2 with topology info)
                auto angle_params = getGFNFFAngleParameters(neighbors[i],
                    center,
                    neighbors[j],
                    current_angle,
                    topo_info);

                angle["fc"] = angle_params.force_constant;
                angle["theta0_ijk"] = angle_params.equilibrium_angle;
                angle["r0_ij"] = (ri - rj).norm(); // Distance i-j
                angle["r0_ik"] = (rk - rj).norm(); // Distance k-j
                // Phase 1.3: No longer using Fourier coefficients (C0/C1/C2)
                // GFN-FF uses simple angle bending formula

                angles.push_back(angle);
            }
        }
    }

    std::cout << "GFN-FF detected " << angles.size() << " angles" << std::endl;

    return angles;
}

bool GFNFF::calculateTopology()
{
    // TODO: Implement GFN-FF specific topology calculation
    // This should include:
    // 1. Bond detection based on covalent radii and GFN-FF rules
    // 2. Ring detection for proper parameter assignment
    // 3. Hybridization state determination
    // 4. Coordination number calculation
    // 5. Formal charge assignment

    // For now, return true - actual implementation would go here
    return true;
}

bool GFNFF::validateMolecule() const
{
    // Check if all atoms are supported by GFN-FF (Z <= 86)
    for (int atom : m_atoms) {
        if (atom < 1 || atom > 86) {
            std::cerr << "Error: Atom type " << atom << " not supported by GFN-FF" << std::endl;
            return false;
        }
    }

    // Check for reasonable geometry
    if (m_geometry_bohr.rows() != m_atomcount || m_geometry.cols() != 3) {
        std::cerr << "Error: Invalid geometry dimensions for GFN-FF" << std::endl;
        return false;
    }

    return true;
}

double GFNFF::convertToHartree(double energy) const
{
    return energy * KCAL_TO_HARTREE;
}

Matrix GFNFF::convertGradientToHartree(const Matrix& gradient) const
{
    return gradient * KCAL_TO_HARTREE * ANGSTROM_TO_BOHR;
}

double GFNFF::getCovalentRadius(int atomic_number) const
{
    // Original GFN-FF covalent radii from gfnff_param.f90 (rad array) in Angström
    // NOTE: Converted to Bohr on return to match m_geometry_bohr units
    static const std::vector<double> covalent_radii = {
        0.32, 0.37, 1.30, 0.99, 0.84, 0.75, 0.71, 0.64, 0.60, 0.62, // H-Ne
        1.60, 1.40, 1.24, 1.14, 1.09, 1.04, 1.00, 1.01, // Na-Ar
        2.00, 1.74, 1.59, 1.48, 1.44, 1.30, 1.29, 1.24, 1.18, // K-Ni
        1.17, 1.22, 1.20, 1.23, 1.20, 1.20, 1.18, 1.17, 1.16, // Cu-Kr
        2.15, 1.90, 1.76, 1.64, 1.56, 1.46, 1.38, 1.36, 1.34, // Rb-Pd
        1.30, 1.36, 1.40, 1.42, 1.40, 1.40, 1.37, 1.36, 1.36, // Ag-Xe
        2.38, 2.06, 1.94, 1.84, 1.90, 1.88, 1.86, 1.85, 1.83, // Cs-Eu
        1.82, 1.81, 1.80, 1.79, 1.77, 1.77, 1.78, 1.74, 1.64, // Gd-Lu
        1.58, 1.50, 1.41, 1.36, 1.32, 1.30, 1.30, 1.32, 1.44, // Hf-Au
        1.45, 1.50, 1.42, 1.48, 1.46, // Hg-Po
        2.42, 2.11, 2.01, 1.90, 1.84, 1.83, 1.80, 1.80, 1.73, // Fr-Am
        1.68, 1.68, 1.68, 1.65, 1.67, 1.73, 1.76, 1.61 // Cm-Lr
    };

    if (atomic_number >= 1 && atomic_number <= static_cast<int>(covalent_radii.size())) {
        // Convert Angström → Bohr to match m_geometry_bohr units
        return covalent_radii[atomic_number - 1] * CurcumaUnit::Length::ANGSTROM_TO_BOHR;
    } else {
        // Fallback for unknown elements (convert to Bohr)
        std::cerr << "Warning: No covalent radius for element " << atomic_number
                  << ", using default 1.0 Å → Bohr" << std::endl;
        return 1.0 * CurcumaUnit::Length::ANGSTROM_TO_BOHR;
    }
}

// Phase 3: EEQ (Electronegativity Equalization) parameters
// Reference: gfnff_param.f90 chi/gam/alp/cnf_angewChem2020 arrays
struct EEQParameters {
    double chi;  // Electronegativity
    double gam;  // Chemical hardness
    double alp;  // Damping parameter
    double cnf;  // CN correction factor
};

GFNFF::EEQParameters GFNFF::getEEQParameters(int atomic_number) const
{
    GFNFF::EEQParameters params;

    // Phase 3: EEQ parameters from gfnff_param.f90 (angewChem2020 parameter set)
    // Reference: S. Spicher, S. Grimme, Angew. Chem. Int. Ed. 2020, 59, 15665-15673

    // chi_angewChem2020(86) - Electronegativity
    static const std::vector<double> chi_eeq = {
        1.227054, 1.451412, 0.813363, 1.062841, 1.186499, // H-B
        1.311555, 1.528485, 1.691201, 1.456784, 1.231037, // C-Ne
        0.772989, 1.199092, 1.221576, 1.245964, 1.248942, // Na-P
        1.301708, 1.312474, 1.247701, 0.781237, 0.940834, // S-Ca
        0.950000, 0.974455, 0.998911, 1.023366, 1.047822, // Sc-Mn
        1.072277, 1.096733, 1.121188, 1.145644, 1.170099, // Fe-Zn
        1.205357, 1.145447, 1.169499, 1.253293, 1.329909, // Ga-Br
        1.116527, 0.950975, 0.964592, 0.897786, 0.932824, // Kr-Zr
        0.967863, 1.002901, 1.037940, 1.072978, 1.108017, // Nb-Rh
        1.143055, 1.178094, 1.213132, 1.205076, 1.075529, // Pd-Sn
        1.206919, 1.303658, 1.332656, 1.179317, 0.789115, // Sb-Cs
        0.798704, 1.127797, 1.127863, 1.127928, 1.127994, // Ba-Nd
        1.128059, 1.128125, 1.128190, 1.128256, 1.128322, // Pm-Tb
        1.128387, 1.128453, 1.128518, 1.128584, 1.128649, // Dy-Yb
        1.128715, 1.128780, 1.129764, 1.130747, 1.131731, // Lu-Re
        1.132714, 1.133698, 1.134681, 1.135665, 1.136648, // Os-Hg
        1.061832, 1.053084, 1.207830, 1.236314, 1.310129, // Tl-At
        1.157380 // Rn
    };

    // gam_angewChem2020(86) - Chemical hardness (FIXED: use correct positive values!)
    // Claude Generated (2025): Replaced incorrect negative values with proper gam_angewChem2020
    static const std::vector<double> gam_eeq = {
        0.473762, 0.455827, // H-He
        0.460630, 0.442857, 0.418838, 0.392544, 0.360555, 0.323150, 0.277053, 0.242707, // Li-Ne
        0.517690, 0.479682, 0.397478, 0.410702, 0.365801, // Na-P
        0.339026, 0.297061, 0.280301, 0.521750, 0.445320, // S-Ca
        0.295042, 0.297024, 0.299006, 0.300987, 0.302969, // Sc-Mn
        0.304951, 0.306933, 0.308915, 0.310896, 0.312878, // Fe-Zn
        0.358152, 0.352705, 0.311348, 0.312305, 0.328874, // Ga-Br
        0.347124, 0.577848, 0.621881, 0.292455, 0.293407, // Kr-Zr
        0.294358, 0.295310, 0.296261, 0.297213, 0.298164, // Nb-Rh
        0.299116, 0.300068, 0.301019, 0.288934, 0.305239, // Pd-Sn
        0.304503, 0.279467, 0.310037, 0.353098, 0.614952, // Sb-Cs
        0.627408, 0.221291, 0.252111, 0.260020, 0.267928, // Ba-Nd
        0.275837, 0.283745, 0.291654, 0.299562, 0.307471, // Pm-Tb
        0.315379, 0.323288, 0.331197, 0.339105, 0.347014, // Dy-Yb
        0.354922, 0.362831, 0.298030, 0.293780, 0.289531, // Lu-Re
        0.285282, 0.281032, 0.276783, 0.272533, 0.268284, // Os-Hg
        0.293447, 0.279865, 0.324131, 0.349663, 0.627862, // Tl-At
        0.646316 // Rn
    };

    // alp_angewChem2020(86) - Damping parameter
    static const std::vector<double> alp_eeq = {
        0.585069, 0.432382, 0.628636, 0.743646, 1.167323, // H-B
        0.903430, 1.278388, 0.905347, 1.067014, 2.941513, // C-Ne
        0.687680, 0.792170, 1.337040, 1.251409, 1.068295, // Na-P
        1.186476, 1.593532, 2.056749, 0.674196, 0.868052, // S-Ca
        0.575052, 0.613424, 0.651796, 0.690169, 0.728541, // Sc-Mn
        0.766913, 0.805285, 0.843658, 0.882030, 0.920402, // Fe-Zn
        0.877178, 1.422350, 1.405901, 1.646860, 2.001970, // Ga-Br
        2.301695, 1.020617, 0.634141, 0.652752, 0.668845, // Kr-Zr
        0.684938, 0.701032, 0.717125, 0.733218, 0.749311, // Nb-Rh
        0.765405, 0.781498, 0.797591, 1.296844, 1.534068, // Pd-Sn
        1.727781, 1.926871, 2.175548, 2.177702, 0.977079, // Sb-Cs
        0.770260, 0.757372, 0.757352, 0.757332, 0.757313, // Ba-Nd
        0.757293, 0.757273, 0.757253, 0.757233, 0.757213, // Pm-Tb
        0.757194, 0.757174, 0.757154, 0.757134, 0.757114, // Dy-Yb
        0.757095, 0.757075, 0.756778, 0.756480, 0.756183, // Lu-Re
        0.755886, 0.755589, 0.755291, 0.754994, 0.754697, // Os-Hg
        0.868029, 1.684375, 2.001040, 2.067331, 2.228923, // Tl-At
        1.874218 // Rn
    };

    // cnf_angewChem2020(86) - CN correction factor
    static const std::vector<double> cnf_eeq = {
        0.008904, 0.004641, 0.048324, 0.080316, -0.051990, // H-B
        0.031779, 0.132184, 0.157353, 0.064120, 0.036540, // C-Ne
        -0.000627, 0.005412, 0.018809, 0.016329, 0.012149, // Na-P
        0.021484, 0.014212, 0.014939, 0.003597, 0.032921, // S-Ca
        -0.021804, -0.022797, -0.023789, -0.024782, -0.025775, // Sc-Mn
        -0.026767, -0.027760, -0.028753, -0.029745, -0.030738, // Fe-Zn
        -0.004189, -0.011113, -0.021305, -0.012311, 0.049781, // Ga-Br
        -0.040533, 0.012872, 0.021056, -0.003395, 0.000799, // Kr-Zr
        0.004992, 0.009186, 0.013379, 0.017573, 0.021766, // Nb-Rh
        0.025960, 0.030153, 0.034347, -0.000052, -0.039776, // Pd-Sn
        0.006661, 0.050424, 0.068985, 0.023470, -0.024950, // Sb-Cs
        -0.033006, 0.058973, 0.058595, 0.058217, 0.057838, // Ba-Nd
        0.057460, 0.057082, 0.056704, 0.056326, 0.055948, // Pm-Tb
        0.055569, 0.055191, 0.054813, 0.054435, 0.054057, // Dy-Yb
        0.053679, 0.053300, 0.047628, 0.041955, 0.036282, // Lu-Re
        0.030610, 0.024937, 0.019264, 0.013592, 0.007919, // Os-Hg
        0.006383, -0.089155, -0.001293, 0.019269, 0.074803, // Tl-At
        0.016657 // Rn
    };

    // Validate atomic number and return parameters
    if (atomic_number >= 1 && atomic_number <= 86) {
        int idx = atomic_number - 1;  // Convert to 0-based indexing
        params.chi = chi_eeq[idx];
        params.gam = gam_eeq[idx];
        params.alp = alp_eeq[idx];
        params.cnf = cnf_eeq[idx];
        params.xi_corr = 0.0;  // No environment correction in simple version
    } else {
        // Fallback for unknown elements (use default values)
        std::cerr << "Warning: No EEQ parameters for element " << atomic_number
                  << ", using default values" << std::endl;
        params.chi = 1.0;
        params.gam = 0.0;
        params.alp = 1.0;
        params.cnf = 0.0;
        params.xi_corr = 0.0;
    }

    return params;
}

GFNFF::GFNFFBondParams GFNFF::getGFNFFBondParameters(int atom1, int atom2, int z1, int z2,
                                                      double distance, const TopologyInfo& topo) const
{
    // Phase 9: Complete assembly with all topology corrections
    // Reference: Fortran gfnff_ini.f90:1087-1285 (complete bond parameter generation)
    //
    // Educational Documentation:
    // ==========================
    // GFN-FF bond parameters with full topology awareness:
    //
    //   r_eq = (r0 + cnfak*CN) * (1 - k1*|ΔEN| - k2*|ΔEN|²) + shift
    //   fc = bond_i * bond_j * ringf * bstrength * fqq * fheavy * fpi * fxh * fcn
    //   α = srb1 * (1 + fsrb2*ΔEN² + srb3*(bstrength-1))
    //
    // Where all corrections use actual topology data:
    //   CN from topo.coordination_numbers
    //   hyb from topo.hybridization
    //   qa from topo.eeq_charges
    //   rings from topo.ring_sizes
    //
    // Literature: Spicher, S.; Grimme, S. Angew. Chem. Int. Ed. 2020, 59, 15665–15673

    GFNFFBondParams params;

    // Original GFN-FF bond force constant parameters (bond_angewChem2020 array)
    // Used for force constant calculation: k = bond(i) * bond(j) * corrections
    // TODO: Bond energy is ~1479x too small - root cause unknown
    // Hypothesis: Missing initialization factor or unit conversion in parameter generation
    static const std::vector<double> bond_params = {
        0.417997, 0.258490, 0.113608, 0.195935, 0.231217, // H-B
        0.385248, 0.379257, 0.339249, 0.330706, 0.120319, // C-Ne
        0.127255, 0.173647, 0.183796, 0.273055, 0.249044, // Na-P
        0.290653, 0.218744, 0.034706, 0.136353, 0.192467, // S-Ca
        0.335860, 0.314452, 0.293044, 0.271636, 0.250228, // Sc-Mn
        0.228819, 0.207411, 0.186003, 0.164595, 0.143187, // Fe-Zn
        0.212434, 0.210451, 0.219870, 0.224618, 0.272206, // Ga-Br
        0.147864, 0.150000, 0.150000, 0.329501, 0.309632, // Kr-Zr
        0.289763, 0.269894, 0.250025, 0.230155, 0.210286, // Nb-Rh
        0.190417, 0.170548, 0.150679, 0.192977, 0.173411, // Pd-Sn
        0.186907, 0.192891, 0.223202, 0.172577, 0.150000, // Sb-Cs
        0.150000, 0.370682, 0.368511, 0.366339, 0.364168, // Ba-Nd
        0.361996, 0.359825, 0.357654, 0.355482, 0.353311, // Pm-Tb
        0.351139, 0.348968, 0.346797, 0.344625, 0.342454, // Dy-Yb
        0.340282, 0.338111, 0.305540, 0.272969, 0.240398, // Lu-Re
        0.207828, 0.175257, 0.142686, 0.110115, 0.077544, // Os-Hg
        0.108597, 0.148422, 0.183731, 0.192274, 0.127706, // Tl-At
        0.086756, 0.150000, 0.150000, 0.370682, 0.368511, // Rn-Ra
        0.366339, 0.364168, 0.361996, 0.359825, 0.357654, // Ac-Am
        0.355482, 0.353311, 0.351139, 0.348968, 0.346797, // Cm-Cf
        0.344625, 0.342454, 0.340282 // Es-Lr
    };

    // Phase 2 NEW: CN-independent base covalent radii (Bohr)
    // Fortran gfnff_rab.f90:82-102
    static const std::vector<double> r0_gfnff = {
        0.55682207, 0.80966997, 2.49092101, 1.91705642, 1.35974851, // H-B
        0.98310699, 0.98423007, 0.76716063, 1.06139799, 1.17736822, // C-Ne
        2.85570926, 2.56149012, 2.31673425, 2.03181740, 1.82568535, // Na-P
        1.73685958, 1.97498207, 2.00136196, 3.58772537, 2.68096221, // S-Ca
        2.23355957, 2.33135502, 2.15870365, 2.10522128, 2.16376162, // Sc-Mn
        2.10804037, 1.96460045, 2.00476257, 2.22628712, 2.43846700, // Fe-Zn
        2.39408483, 2.24245792, 2.05751204, 2.15427677, 2.27191920, // Ga-Br
        2.19722638, 3.80910350, 3.26020971, 2.99716916, 2.71707818, // Kr-Zr
        2.34950167, 2.11644818, 2.47180659, 2.32198800, 2.32809515, // Nb-Rh
        2.15244869, 2.55958313, 2.59141300, 2.62030465, 2.39935278, // Pd-Sn
        2.56912355, 2.54374096, 2.56914830, 2.53680807, 4.24537037, // Sb-Cs
        3.66542289, 3.19903011, 2.80000000, 2.80000000, 2.80000000, // Ba-Nd (58-60 placeholder)
        2.80000000, 2.80000000, 2.80000000, 2.80000000, 2.80000000, // Pm-Tb (61-65)
        2.80000000, 2.80000000, 2.80000000, 2.80000000, 2.80000000, // Dy-Yb (66-70)
        2.80000000, 2.34880037, 2.37597108, 2.49067697, 2.14100577, // Lu-Re (71-75)
        2.33473532, 2.19498900, 2.12678348, 2.34895048, 2.33422774, // Os-Hg (76-80)
        2.86560827, 2.62488837, 2.88376127, 2.75174124, 2.83054552, // Tl-At (81-85)
        2.63264944 // Rn (86)
    };

    // Phase 2 NEW: CN-dependent correction factors
    // Fortran gfnff_rab.f90:103-122
    static const std::vector<double> cnfak_gfnff = {
        0.17957827,  0.25584045, -0.02485871,  0.00374217,  0.05646607, // H-B
        0.10514203,  0.09753494,  0.30470380,  0.23261783,  0.36752208, // C-Ne
        0.00131819, -0.00368122, -0.01364510,  0.04265789,  0.07583916, // Na-P
        0.08973207, -0.00589677,  0.13689929, -0.01861307,  0.11061699, // S-Ca
        0.10201137,  0.05426229,  0.06014681,  0.05667719,  0.02992924, // Sc-Mn
        0.03764312,  0.06140790,  0.08563465,  0.03707679,  0.03053526, // Fe-Zn
       -0.00843454,  0.01887497,  0.06876354,  0.01370795, -0.01129196, // Ga-Br
        0.07226529,  0.01005367,  0.01541506,  0.05301365,  0.07066571, // Kr-Zr
        0.07637611,  0.07873977,  0.02997732,  0.04745400,  0.04582912, // Nb-Rh
        0.10557321,  0.02167468,  0.05463616,  0.05370913,  0.05985441, // Pd-Sn
        0.02793994,  0.02922983,  0.02220438,  0.03340460, -0.04110969, // Sb-Cs
       -0.01987240,  0.07260201,  0.07700000,  0.07700000,  0.07700000, // Ba-Nd (58-60)
        0.07700000,  0.07700000,  0.07700000,  0.07700000,  0.07700000, // Pm-Tb (61-65)
        0.07700000,  0.07700000,  0.07700000,  0.07700000,  0.07700000, // Dy-Yb (66-70)
        0.07700000,  0.08379100,  0.07314553,  0.05318438,  0.06799334, // Lu-Re (71-75)
        0.04671159,  0.06758819,  0.09488437,  0.07556405,  0.13384502, // Os-Hg (76-80)
        0.03203572,  0.04235009,  0.03153769, -0.00152488,  0.02714675, // Tl-At (81-85)
        0.04800662 // Rn (86)
    };

    // =============================================================================
    // GFN-FF UNIT SYSTEM NOTES
    // =============================================================================
    //
    // CURRENT IMPLEMENTATION (Option A - Geometry Conversion):
    // - All GFN-FF literature parameters are in BOHR (r0, cnfak, distances)
    // - Curcuma input geometry is in ANGSTRÖM (XYZ standard)
    // - Solution: Convert geometry to Bohr in InitialiseMolecule()
    // - Member variable m_geometry_bohr stores converted geometry
    // - Gradients are converted back to Angström for output
    //
    // ALTERNATIVE APPROACH (Option B - Parameter Conversion):
    // - Convert all parameters from Bohr → Angström once at startup
    // - All calculations run in Angström (no geometry conversion needed)
    // - Advantages: Consistent units, simpler debugging, no dual geometry storage
    // - Disadvantages: Parameter values differ from literature (needs documentation)
    // - To implement: Multiply r0_gfnff by 0.529177, divide cnfak_gfnff by 0.529177
    // - Arrays to convert: r0_gfnff[86], cnfak_gfnff[86]
    //
    // Current choice: Option A keeps parameters identical to Fortran reference
    // =============================================================================

    // Phase 2 NEW: GFN-FF electronegativities (different from Pauling!)
    // Fortran gfnff_rab.f90:62-81
    static const std::vector<double> en_gfnff = {
        2.30085633, 2.78445145, 1.52956084, 1.51714704, 2.20568300, // H-B
        2.49640820, 2.81007174, 4.51078438, 4.67476223, 3.29383610, // C-Ne
        2.84505365, 2.20047950, 2.31739628, 2.03636974, 1.97558064, // Na-P
        2.13446570, 2.91638164, 1.54098156, 2.91656301, 2.26312147, // S-Ca
        2.25621439, 1.32628677, 2.27050569, 1.86790977, 2.44759456, // Sc-Mn
        2.49480042, 2.91545568, 3.25897750, 2.68723778, 1.86132251, // Fe-Zn
        2.01200832, 1.97030722, 1.95495427, 2.68920990, 2.84503857, // Ga-Br
        2.61591858, 2.64188286, 2.28442252, 1.33011187, 1.19809388, // Kr-Zr
        1.89181390, 2.40186898, 1.89282464, 3.09963488, 2.50677823, // Nb-Rh
        2.61196704, 2.09943450, 2.66930105, 1.78349472, 2.09634533, // Pd-Sn
        2.00028974, 1.99869908, 2.59072029, 2.54497829, 2.52387890, // Sb-Cs
        2.30204667, 1.60119300, 2.00000000, 2.00000000, 2.00000000, // Ba-Nd
        2.00000000, 2.00000000, 2.00000000, 2.00000000, 2.00000000, // Pm-Tb
        2.00000000, 2.00000000, 2.00000000, 2.00000000, 2.00000000, // Dy-Yb
        2.00000000, 2.30089349, 1.75039077, 1.51785130, 2.62972945, // Lu-Re
        2.75372921, 2.62540906, 2.55860939, 3.32492356, 2.65140898, // Os-Hg
        1.52014458, 2.54984804, 1.72021963, 2.69303422, 1.81031095, // Tl-At
        2.34224386 // Rn
    };

    // Phase 2 NEW: Row-dependent EN polynomial coefficients (scaled by 10^-3 in Fortran)
    // Fortran gfnff_rab.f90:125-136
    // Index: [row-1][0 or 1]  (row = 1..6 for H-He, Li-Ne, Na-Ar, K-Kr, Rb-Xe, Cs-Rn)
    static const double p_enpoly[6][2] = {
        { 29.84522887, -8.87843763 },  // Row 1: H, He
        { -1.70549806,  2.10878369 },  // Row 2: Li-Ne
        {  6.54013762,  0.08009374 },  // Row 3: Na-Ar
        {  6.39169003, -0.85808076 },  // Row 4: K-Kr
        {  6.00000000, -1.15000000 },  // Row 5: Rb-Xe (extrapolated)
        {  5.60000000, -1.30000000 }   // Row 6: Cs-Rn (extrapolated)
    };

    // Helper lambda: Get periodic table row (1-6)
    // Fortran gfnff_rab.f90:165-185 (iTabRow6 function)
    auto getPeriodicTableRow = [](int z) -> int {
        if (z <= 0) return 0;
        if (z <= 2) return 1;        // H, He
        if (z <= 10) return 2;       // Li-Ne
        if (z <= 18) return 3;       // Na-Ar
        if (z <= 36) return 4;       // K-Kr
        if (z <= 54) return 5;       // Rb-Xe
        return 6;                    // Cs-Rn and beyond
    };

    // Step 1: Get bond force constant parameters (will be used in Phase 9)
    double bond_param_1 = (z1 >= 1 && z1 <= static_cast<int>(bond_params.size())) ? bond_params[z1 - 1] : 0.15;
    double bond_param_2 = (z2 >= 1 && z2 <= static_cast<int>(bond_params.size())) ? bond_params[z2 - 1] : 0.15;

    // Step 2: Phase 9 - Use actual CN-dependent radii from topology
    // Fortran gfnff_rab.f90:147-148
    double cn1 = topo.coordination_numbers[atom1];
    double cn2 = topo.coordination_numbers[atom2];

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info(fmt::format("=== GFN-FF Bond Parameters: Atom {} (Z={}) - Atom {} (Z={}) ===",
                                         atom1, z1, atom2, z2));
        CurcumaLogger::info(fmt::format("  Topology: CN1={:.2f}, CN2={:.2f}, hyb1={}, hyb2={}, ring1={}, ring2={}",
                                         cn1, cn2, topo.hybridization[atom1], topo.hybridization[atom2],
                                         topo.ring_sizes[atom1], topo.ring_sizes[atom2]));
    }

    double r0_1 = (z1 >= 1 && z1 <= static_cast<int>(r0_gfnff.size())) ? r0_gfnff[z1 - 1] : 2.0;
    double r0_2 = (z2 >= 1 && z2 <= static_cast<int>(r0_gfnff.size())) ? r0_gfnff[z2 - 1] : 2.0;
    double cnfak_1 = (z1 >= 1 && z1 <= static_cast<int>(cnfak_gfnff.size())) ? cnfak_gfnff[z1 - 1] : 0.0;
    double cnfak_2 = (z2 >= 1 && z2 <= static_cast<int>(cnfak_gfnff.size())) ? cnfak_gfnff[z2 - 1] : 0.0;

    double ra = r0_1 + cnfak_1 * cn1;  // CN-dependent radius A
    double rb = r0_2 + cnfak_2 * cn2;  // CN-dependent radius B

    std::cout << "Summ of atomic radius and correction for atom 1: " << r0_1 <<  " " << cnfak_1 << " " << cn1 << " " << ra << std::endl;
    std::cout << "Summ of atomic radius and correction for atom 2: " << r0_2 <<  " " << cnfak_2 << " " << cn2 << " " << rb << std::endl;

    // Step 3: Phase 2 - Row-dependent EN correction
    // Fortran gfnff_rab.f90:144-152
    double en1 = (z1 >= 1 && z1 <= static_cast<int>(en_gfnff.size())) ? en_gfnff[z1 - 1] : 2.0;
    double en2 = (z2 >= 1 && z2 <= static_cast<int>(en_gfnff.size())) ? en_gfnff[z2 - 1] : 2.0;

    int row1 = getPeriodicTableRow(z1);
    int row2 = getPeriodicTableRow(z2);

    // EN polynomial coefficients (Fortran multiplies by 0.005)
    double k1 = 0.005 * (p_enpoly[row1 - 1][0] + p_enpoly[row2 - 1][0]);
    double k2 = 0.005 * (p_enpoly[row1 - 1][1] + p_enpoly[row2 - 1][1]);

    double en_diff = std::abs(en1 - en2);
    double ff = 1.0 - k1 * en_diff - k2 * en_diff * en_diff;

    // Step 4: Phase 2 - Final equilibrium distance with shift correction
    // Fortran gfnff_ini.f90:1235 & 1248: rab(k) = (ra + rb + shift) * ff, then r0 = rab(k)*0.529167
    //
    // Shift calculation (Fortran gfnff_ini.f90:1063-1231):
    //   1. Base shift: gen%rabshift = -0.110 (general shift in Bohr)
    //   2. XH correction: if(ia.eq.1 or ja.eq.1) shift += gen%rabshifth = -0.050
    //   3. Hypervalent: if(bbtyp==4) shift = gen%hyper_shift
    //   4. Ring effects: if X-sp3 in ring, shift -= 0.022
    //   5. Heavy atom effects (Z>36): shift += gen%hshift5
    //
    // For now: implement basic rabshift + XH correction
    // Full ring/hypervalent/heavy atom corrections will be added in Phase 10
    double gen_rabshift = -0.110;    // Fortran: gen%rabshift
    double gen_rabshifth = -0.050;   // Fortran: gen%rabshifth (XH bonds)

    double shift = 0.0;
    // XH bond correction
    if (z1 == 1 || z2 == 1) {
        shift = gen_rabshifth;  // XH uses special shift
    }
    // Additional shifts (hypervalent, heavy atoms, rings) - TODO Phase 10

    double rabshift = gen_rabshift + shift;  // Total shift in Bohr
    // NOTE: Result is in Bohr (ra, rb are in Bohr) - NO conversion needed since we use m_geometry_bohr
    params.equilibrium_distance = (ra + rb + rabshift) * ff;

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info(fmt::format("  Equilibrium Distance: r_A={:.3f}, r_B={:.3f}, EN_corr={:.3f} -> r_eq={:.3f} Bohr",
                                         ra, rb, ff, params.equilibrium_distance));
    }

    // Step 5: Phase 3 - Hybridization-based bond strength
    // Reference: Fortran gfnff_param.f90:733-815 and gfnff_ini.f90:1127-1142
    //
    // Educational Documentation:
    // ==========================
    // GFN-FF bond strength depends on hybridization of both atoms:
    //   bstrength = bsmat[max(hyb_i, hyb_j), min(hyb_i, hyb_j)]
    //
    // Hybridization codes (Fortran):
    //   0 = unknown/default (treated as sp3)
    //   1 = sp (linear)
    //   2 = sp2 (trigonal planar)
    //   3 = sp3 (tetrahedral)
    //   5 = hypervalent (sp3d, sp3d2)
    //
    // Bond strength values (bstren):
    //   1.00 = single bond
    //   1.24 = double bond (sp2-sp2)
    //   1.98 = triple bond (sp-sp)
    //   1.22 = hypervalent
    //
    // Literature: Spicher, S.; Grimme, S. Angew. Chem. Int. Ed. 2020

    // Phase 3 NEW: Bond strength base values
    // Fortran gfnff_param.f90:733-741
    static const double bstren[9] = {
        0.0,   // Index 0 unused
        1.00,  // [1] single bond
        1.24,  // [2] double bond
        1.98,  // [3] triple bond
        1.22,  // [4] hypervalent bond
        1.00,  // [5] M-X (metal-ligand)
        0.78,  // [6] M eta
        3.40,  // [7] M-M
        3.40   // [8] M-M
    };

    // Phase 3 NEW: Hybridization-dependent bond strength matrix (4×4)
    // Fortran gfnff_param.f90:804-814
    // split0 = 0.67, split1 = 0.33 for mixed hybridizations
    static const double bsmat[4][4] = {
        // hyb=0 (unknown/sp3)
        { 1.000, 1.323, 1.079, 1.000 },  // vs. hyb=0,1,2,3
        // hyb=1 (sp)
        { 1.323, 1.980, 1.484, 1.323 },  // vs. hyb=0,1,2,3
        // hyb=2 (sp2)
        { 1.079, 1.484, 1.240, 1.079 },  // vs. hyb=0,1,2,3
        // hyb=3 (sp3)
        { 1.000, 1.323, 1.079, 1.000 }   // vs. hyb=0,1,2,3
    };
    // Computed as: bsmat[i][j] = split0*bstren[bond_i] + split1*bstren[bond_j]
    // Example: bsmat[1][0] = 0.67*1.00 + 0.33*1.98 = 1.323

    // Phase 9: Use actual hybridization from topology
    // Curcuma uses hyb=1,2,3 (sp, sp2, sp3), Fortran uses hyb=0,1,2,3
    int hyb1 = topo.hybridization[atom1];
    int hyb2 = topo.hybridization[atom2];

    // Map Curcuma hybridization (1,2,3) to Fortran (0,1,2,3)
    // Curcuma: 1=sp, 2=sp2, 3=sp3 → Fortran: 1=sp, 2=sp2, 3=sp3, 0=unknown
    // For safety, map Curcuma 3 → Fortran 3 (sp3)
    int hyb1_fortran = (hyb1 >= 1 && hyb1 <= 3) ? hyb1 : 3;
    int hyb2_fortran = (hyb2 >= 1 && hyb2 <= 3) ? hyb2 : 3;

    // Get hybridization indices for lookup (needed for both H and non-H bonds)
    int hybi = std::max(hyb1_fortran, hyb2_fortran);  // Max hyb index
    int hybj = std::min(hyb1_fortran, hyb2_fortran);  // Min hyb index

    // CRITICAL FIX (Phase 11): Special handling for hydrogen bonds!
    // Fortran gfnff_param.f90 has bsmat(1,1)=1.98 for sp-sp (triple bond)
    // But Fortran gfnff_ini.f90:1069 sets bbtyp=3 for ANY sp-X bond!
    // For H-H (both atoms Z=1), this incorrectly uses triple bond strength
    // Solution: H-H and H-X bonds should be SINGLE BONDS (bstrength=1.0)
    double bstrength;

    if (z1 == 1 || z2 == 1) {
        // Hydrogen bonds are always single bonds in GFN-FF
        bstrength = bstren[1];  // 1.00 (single bond)
    } else {
        // Get bond strength from hybridization matrix (Fortran gfnff_ini.f90:1127-1133)
        if (hybi == 5 || hybj == 5) {
            // Hypervalent atoms
            bstrength = bstren[4];  // 1.22
        } else {
            // Normal hybridization: lookup in bsmat
            // Since Curcuma uses 1-3, and bsmat is indexed 0-3, we need to map:
            // Curcuma hyb=1,2,3 → bsmat indices 1,2,3 (Fortran hyb=1,2,3)
            // But bsmat also has index 0 for unknown → use hyb=3 for unknown
            bstrength = bsmat[hybi][hybj];
        }
    }

    // Phase 3: Special cases (Fortran gfnff_ini.f90:1134-1142)
    // N-sp2 correction: if one atom is sp3 (hyb=3) and other is sp2 (hyb=2) and it's nitrogen
    if ((hybi == 3 && hybj == 2 && z1 == 7) || (hybi == 3 && hybj == 2 && z2 == 7)) {
        bstrength = bstren[2] * 1.04;  // N-sp2: 1.24 * 1.04 = 1.2896
    }

    // Bridging atoms (will be detected in Phase 6)
    bool is_bridge = false;  // Placeholder
    if (is_bridge) {
        // Group 7 (halogens): bstrength *= 0.50
        // H or F bridging: bstrength *= 0.30
        // (Will be implemented in Phase 6)
    }

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info(fmt::format("  Hybridization: hyb1={}, hyb2={} -> bstrength={:.3f}",
                                         hyb1, hyb2, bstrength));
    }

    // Step 6: Phase 5 - EEQ charge-dependent force constant correction
    // Reference: Fortran gfnff_ini.f90:1185-1186
    //
    // Educational Documentation:
    // ==========================
    // GFN-FF force constant depends on EEQ atomic partial charges:
    //   fqq = 1 + qfacbm0 * tanh(15 * qa_i * qa_j * 70)
    //
    // Where:
    //   qa_i, qa_j = EEQ atomic partial charges (electron units)
    //   qfacbm0 = 0.047 (empirical scaling factor)
    //   70.0 = charge product amplification factor
    //   15.0 = sigmoid steepness parameter
    //
    // Physical meaning: Opposite charges strengthen bonds, like charges weaken them
    // Literature: Spicher, S.; Grimme, S. Angew. Chem. Int. Ed. 2020

    double fqq = 1.0;  // Default: no charge correction

    // Phase 9: Use actual EEQ charges from topology
    double qa1 = topo.eeq_charges[atom1];
    double qa2 = topo.eeq_charges[atom2];

    // Fortran formula (sigmoid function for smooth charge-dependence)
    // fqq = 1.0 + qfacbm0 * exp(-15*qafac) / (1 + exp(-15*qafac))
    // This is equivalent to: fqq = 1.0 + qfacbm0 * tanh(15*qafac/2)
    double qfacbm0 = 0.047;  // Fortran gfnff_param.f90:772
    double qafac = qa1 * qa2 * 70.0;
    double exp_term = std::exp(-15.0 * qafac);
    fqq = 1.0 + qfacbm0 * exp_term / (1.0 + exp_term);

    // Step 7: Phase 6 - Ring strain corrections and XH special cases
    // Reference: Fortran gfnff_ini.f90:1150-1165, 1278-1279
    //
    // Educational Documentation:
    // ==========================
    // GFN-FF bond strength modified by ring strain:
    //   ringf = 1 + fringbo * (6 - ring_size)²
    //
    // Where:
    //   ring_size = smallest ring containing this bond (3, 4, 5, 6+)
    //   fringbo = 0.020 (ring strain scaling factor)
    //   Reference size = 6 (benzene-like, no strain)
    //
    // XH bond corrections (fxh):
    //   3-ring CH: +5%  (strained C-H bonds)
    //   Aldehyde CH: -5% (weak C-H in CHO group)
    //   BH: +10%  (strong B-H bonds)
    //   NH: +6%   (strong N-H bonds)
    //   OH: -7%   (weak O-H bonds)
    //
    // Literature: Spicher, S.; Grimme, S. Angew. Chem. Int. Ed. 2020

    // Phase 9: Ring strain correction with actual topology data
    double ringf = 1.0;  // Default: no ring strain

    // Use ring size from topology (smallest ring containing either atom)
    int ring_size = std::max(topo.ring_sizes[atom1], topo.ring_sizes[atom2]);

    if (ring_size > 0) {
        // Fortran gfnff_ini.f90:1279
        double fringbo = 0.020;  // Fortran gfnff_param.f90:800
        ringf = 1.0 + fringbo * std::pow(6.0 - ring_size, 2);
    }
    // Examples:
    //   3-ring: ringf = 1.0 + 0.020*(6-3)² = 1.180 (18% stronger)
    //   4-ring: ringf = 1.0 + 0.020*(6-4)² = 1.080 (8% stronger)
    //   5-ring: ringf = 1.0 + 0.020*(6-5)² = 1.020 (2% stronger)
    //   6-ring: ringf = 1.0 + 0.020*(6-6)² = 1.000 (no strain)

    // Phase 6 NEW: XH bond special cases
    // Fortran gfnff_ini.f90:1150-1165
    double fxh = 1.0;  // Default: no XH correction

    // Check if this is an X-H bond (one atom is hydrogen)
    if (z1 == 1 || z2 == 1) {
        int heavy_atom = (z1 == 1) ? z2 : z1;  // The non-hydrogen atom

        // TEMPORARY: Ring and functional group detection needed for full XH corrections
        // For now, implement element-specific corrections only
        bool is_3ring = false;  // Placeholder (will be from ring detection in Phase 9)
        bool is_aldehyde = false;  // Placeholder (requires functional group detection)

        if (heavy_atom == 6 && is_3ring) {
            // 3-ring CH: stronger due to ring strain
            fxh = 1.05;  // +5%
        } else if (heavy_atom == 6 && is_aldehyde) {
            // Aldehyde CH: weaker
            fxh = 0.95;  // -5%
        } else if (heavy_atom == 5) {
            // B-H bond: strong
            fxh = 1.10;  // +10%
        } else if (heavy_atom == 7) {
            // N-H bond: moderate strength
            fxh = 1.06;  // +6%
        } else if (heavy_atom == 8) {
            // O-H bond: weak (hydrogen bonding)
            fxh = 0.93;  // -7%
        }
    }

    // Step 8: Phase 8 - CN-dependent heavy atom corrections
    // Reference: Fortran gfnff_ini.f90:1181-1184
    //
    // Educational Documentation:
    // ==========================
    // GFN-FF weakens bonds for highly coordinated heavy atoms:
    //   fcn = 1 / (1 + 0.007 * nb20_A²) / (1 + 0.007 * nb20_B²)
    //
    // Where:
    //   nb20 = number of neighbors within 20 Bohr cutoff
    //   0.007 = empirical scaling factor
    //   Only applied for heavy atoms (Z > 10)
    //
    // Physical meaning: Highly coordinated atoms have weaker individual bonds
    // Example: 6-coordinate metal complex has weaker M-L bonds than 4-coordinate
    //
    // Literature: Spicher, S.; Grimme, S. Angew. Chem. Int. Ed. 2020

    double fcn = 1.0;  // Default: no CN correction

    // Phase 9: Use actual coordination numbers as approximation for nb20
    // TODO: Implement exact nb20 (neighbors within 20 Bohr cutoff) for precision
    // For now, use CN as reasonable approximation
    int nb20_1 = static_cast<int>(std::round(cn1));
    int nb20_2 = static_cast<int>(std::round(cn2));

    // Only apply to heavy atoms (Z > 10, i.e., beyond neon)
    // Fortran gfnff_ini.f90:1181-1184
    if (z1 > 10 && z2 > 10) {
        fcn /= (1.0 + 0.007 * nb20_1 * nb20_1);
        fcn /= (1.0 + 0.007 * nb20_2 * nb20_2);
    }
    // Examples:
    //   nb20=4: fcn = 1/(1+0.007*16) = 0.898 (10% weaker)
    //   nb20=6: fcn = 1/(1+0.007*36) = 0.799 (20% weaker)
    //   nb20=8: fcn = 1/(1+0.007*64) = 0.689 (31% weaker)

    // Step 9: Phase 7 - Metal-specific logic (fheavy, metal-specific shifts, alpha sign flip)
    // Reference: Fortran gfnff_ini.f90:1188-1265
    //
    // Educational Documentation:
    // ==========================
    // GFN-FF has extensive metal chemistry corrections:
    //
    // Metal classification (Fortran param_gfnff.f90:408-415):
    //   metal[Z] = 0: non-metal
    //   metal[Z] = 1: main group metal (Li, Na, Mg, Al, Sn, Pb, etc.)
    //   metal[Z] = 2: transition metal (TM: Sc-Zn, Y-Cd, Hf-Hg)
    //
    // Metal types (mtyp):
    //   0: non-metal
    //   1: Group 1 (alkali metals)
    //   2: Group 2 (alkaline earth)
    //   3: Main group metal (Al, Ga, In, Sn, Pb, Bi, Po)
    //   4: Transition metal
    //
    // fheavy corrections (TM-ligand bond strength):
    //   TM-heavy (Z>10): 0.65 (weaken metal-heavy atom bonds)
    //   TM-P: 1.60 (strengthen P ligands)
    //   TM-chalcogen: 0.85 (S, Se, Te ligands)
    //   TM-halogen: 1.30 (F, Cl, Br, I ligands)
    //
    // Literature: Spicher, S.; Grimme, S. Angew. Chem. Int. Ed. 2020

    // Phase 7 NEW: Metal arrays (Fortran param_gfnff.f90)
    static const int metal_type[86] = {
        // 0=non-metal, 1=main group metal, 2=transition metal
        0,                                                                0, // H-He
        1,1,                                               0, 0, 0, 0, 0, 0, // Li-Ne
        1,1,                                               1, 0, 0, 0, 0, 0, // Na-Ar
        1,1,2,                2, 2, 2, 2, 2, 2, 2, 2, 2,   1, 0, 0, 0, 0, 0, // K-Kr
        1,2,2,                2, 2, 2, 2, 2, 2, 2, 2, 2,   1, 1, 0, 0, 0, 0, // Rb-Xe
        1,2,2, 2,2,2,2,2,2,2,2,2,2,2,2,2,2, 2, 2, 2, 2, 2, 2, 2, 2, 2,   1, 1, 1, 1, 0, 0  // Cs-Rn
    };

    static const int periodic_group[86] = {
        // 1-8: main group, negative: transition metals (d-block)
        1,                                                                   8, // H-He
        1,2,                                                  3, 4, 5, 6, 7, 8, // Li-Ne
        1,2,                                                  3, 4, 5, 6, 7, 8, // Na-Ar
        1,2,-3,                 -4,-5,-6,-7,-8,-9,-10,-11,-12,3, 4, 5, 6, 7, 8, // K-Kr
        1,2,-3,                 -4,-5,-6,-7,-8,-9,-10,-11,-12,3, 4, 5, 6, 7, 8, // Rb-Xe
        1,2,-3,  -3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3, -4,-5,-6,-7,-8,-9,-10,-11,-12,3, 4, 5, 6, 7, 8  // Cs-Rn
    };

    double fheavy = 1.0;  // Default: no metal-ligand correction
    double fpi = 1.0;     // Default: no pi-bond order correction (Phase 4 will override)
    double metal_shift = 0.0;  // Additional equilibrium distance shift for metals

    // Get metal types for both atoms
    int imetal1 = (z1 >= 1 && z1 <= 86) ? metal_type[z1 - 1] : 0;
    int imetal2 = (z2 >= 1 && z2 <= 86) ? metal_type[z2 - 1] : 0;
    int group1 = (z1 >= 1 && z1 <= 86) ? periodic_group[z1 - 1] : 0;
    int group2 = (z2 >= 1 && z2 <= 86) ? periodic_group[z2 - 1] : 0;

    // Metal type classification (Fortran gfnff_ini.f90:1158-1167)
    // CRITICAL: H must not be treated as alkali metal even though it's in group 1!
    int mtyp1 = 0;  // 0=non-metal, 1=Group1, 2=Group2, 3=main group metal, 4=TM
    int mtyp2 = 0;

    // Hydrogen is explicitly non-metal (Z=1 → mtyp=0)
    if (z1 > 1 && group1 == 1) mtyp1 = 1;  // Li, Na, K, Rb, Cs (NOT hydrogen!)
    else if (z1 > 1 && group1 == 2) mtyp1 = 2;  // Be, Mg, Ca, Sr, Ba
    else if (group1 > 2 && imetal1 == 1) mtyp1 = 3;  // Al, Ga, In, Sn, Pb, Bi, Po
    else if (imetal1 == 2) mtyp1 = 4;  // Transition metals

    if (z2 > 1 && group2 == 1) mtyp2 = 1;  // Li, Na, K, Rb, Cs (NOT hydrogen!)
    else if (z2 > 1 && group2 == 2) mtyp2 = 2;
    else if (group2 > 2 && imetal2 == 1) mtyp2 = 3;
    else if (imetal2 == 2) mtyp2 = 4;

    // Phase 7: Metal-ligand corrections (Fortran gfnff_ini.f90:1212-1245)
    if (mtyp1 == 4 || mtyp2 == 4) {  // At least one atom is a transition metal
        // fheavy corrections for TM-ligand bonds
        if (imetal1 == 2 && z2 > 10) fheavy = 0.65;  // TM with heavy ligand
        if (imetal2 == 2 && z1 > 10) fheavy = 0.65;

        if (imetal1 == 2 && z2 == 15) fheavy = 1.60;  // TM-P (phosphine ligands)
        if (imetal2 == 2 && z1 == 15) fheavy = 1.60;

        if (imetal1 == 2 && group2 == 6) fheavy = 0.85;  // TM-chalcogen (S, Se, Te)
        if (imetal2 == 2 && group1 == 6) fheavy = 0.85;

        if (imetal1 == 2 && group2 == 7) fheavy = 1.30;  // TM-halogen (F, Cl, Br, I)
        if (imetal2 == 2 && group1 == 7) fheavy = 1.30;

        // M-CO and M-CN corrections (Fortran gfnff_ini.f90:1226-1245)
        // Requires sp hybridization check
        if (imetal2 == 2 && hyb1 == 1) {  // TM bonded to sp atom
            if (z1 == 6) {  // M-CO (carbon monoxide ligand)
                fpi = 1.5;
                metal_shift = -0.45;
            }
            if (z1 == 7 && nb20_1 != 1) {  // M-CN (cyanide, not terminal N)
                fpi = 0.4;
                metal_shift = 0.47;
            }
        }
        if (imetal1 == 2 && hyb2 == 1) {  // sp atom bonded to TM
            if (z2 == 6) {  // M-CO
                fpi = 1.5;
                metal_shift = -0.45;
            }
            if (z2 == 7 && nb20_2 != 1) {  // M-CN
                fpi = 0.4;
                metal_shift = 0.47;
            }
        }

        // fxh corrections for M-H bonds (Fortran gfnff_ini.f90:1220-1225)
        int row1 = getPeriodicTableRow(z1);
        int row2 = getPeriodicTableRow(z2);

        if (imetal1 == 2 && z2 == 1) {  // TM-H
            if (row1 <= 5) fxh = 0.80;  // 3d/4d metals (Sc-Cd, Y-Cd)
            else fxh = 1.00;            // 5d metals (Hf-Hg)
        }
        if (imetal2 == 2 && z1 == 1) {  // H-TM
            if (row2 <= 5) fxh = 0.80;
            else fxh = 1.00;
        }
    }

    // Main group metal-H corrections
    if ((imetal1 == 1 && z2 == 1) || (imetal2 == 1 && z1 == 1)) {
        fxh = 1.20;  // Main group metal-H (e.g., Al-H, Sn-H)
    }

    // Metal-specific equilibrium distance shifts (Fortran gfnff_ini.f90:1246-1253)
    if (imetal1 == 2) metal_shift += 0.15;  // TM shift (metal2_shift)
    if (imetal2 == 2) metal_shift += 0.15;

    if (imetal1 == 1 && group1 <= 2) metal_shift += 0.20;  // Group 1+2 shift (metal1_shift)
    if (imetal2 == 1 && group2 <= 2) metal_shift += 0.20;

    if (mtyp1 == 3) metal_shift += 0.05;  // Main group metal shift (metal3_shift)
    if (mtyp2 == 3) metal_shift += 0.05;

    // Apply metal shift to equilibrium distance (metal_shift is already in Bohr)
    params.equilibrium_distance += metal_shift;

    // Metal-specific fcn corrections (Fortran gfnff_ini.f90:1254-1259)
    // Different CN-dependence for metals vs. non-metals
    if (mtyp1 > 0 && mtyp1 < 3) {  // Group 1+2 metals
        fcn /= (1.0 + 0.100 * nb20_1 * nb20_1);  // Stronger CN dependence
    } else if (mtyp1 == 3) {  // Main group metal
        fcn /= (1.0 + 0.030 * nb20_1 * nb20_1);
    } else if (mtyp1 == 4) {  // Transition metal
        fcn /= (1.0 + 0.036 * nb20_1 * nb20_1);
    }

    if (mtyp2 > 0 && mtyp2 < 3) {  // Group 1+2 metals
        fcn /= (1.0 + 0.100 * nb20_2 * nb20_2);
    } else if (mtyp2 == 3) {  // Main group metal
        fcn /= (1.0 + 0.030 * nb20_2 * nb20_2);
    } else if (mtyp2 == 4) {  // Transition metal
        fcn /= (1.0 + 0.036 * nb20_2 * nb20_2);
    }

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info(fmt::format("  Correction Factors: fqq={:.3f}, ringf={:.3f}, fxh={:.3f}, fcn={:.3f}, fheavy={:.3f}, fpi={:.3f}",
                                         fqq, ringf, fxh, fcn, fheavy, fpi));
        if (mtyp1 > 0 || mtyp2 > 0) {
            CurcumaLogger::info(fmt::format("  Metal Types: mtyp1={} (imetal={}), mtyp2={} (imetal={}), metal_shift={:.3f}",
                                             mtyp1, imetal1, mtyp2, imetal2, metal_shift));
        }
    }

    // Step 10: Force constant (COMPLETE with all 9 factors!)
    // Fortran gfnff_ini.f90:1285: fc = -bond(i)*bond(j) * ringf * bstrength * fqq * fheavy * fpi * fxh * fcn
    // CRITICAL: Negative sign required! Bond energy formula: E = k_b * exp(-α*(r-r₀)²)
    // With k_b < 0, bond energy becomes negative (attractive) at equilibrium
    // Phase 7 completes this formula with fheavy and updated fxh/fcn
    params.force_constant = -(bond_param_1 * bond_param_2 * bstrength * fqq * ringf * fheavy * fpi * fxh * fcn);

    // Step 11: Alpha parameter with metal-specific sign flip (Fortran gfnff_param.f90:642-644, gfnff_ini.f90:1240)
    // CRITICAL PARAMETERS FROM FORTRAN (these were WRONG in earlier implementation!)
    double srb1 = 0.3731;   // Fortran: gen%srb1
    double srb2 = 0.3171;   // Fortran: gen%srb2
    double srb3 = 0.2538;   // Fortran: gen%srb3

    // EN-dependence scaling with metal-specific SIGN FLIP (Fortran gfnff_ini.f90:1227-1234)
    double fsrb2;
    if (mtyp1 == 4 || mtyp2 == 4) {
        // Transition metals: INVERSE EN dependence (negative sign!)
        fsrb2 = -srb2 * 0.22;  // Fortran: -gen.srb2*0.22
    } else if (mtyp1 > 0 || mtyp2 > 0) {
        // Other metals: normal but scaled EN dependence
        fsrb2 = srb2 * 0.28;   // Fortran: gen.srb2*0.28
    } else {
        // Non-metals: standard EN dependence (Fortran gfnff_ini.f90:1233)
        fsrb2 = srb2;  // Just srb2 directly
    }

    params.alpha = srb1 * (1.0 + fsrb2 * en_diff * en_diff + srb3 * (bstrength - 1.0));

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::success(fmt::format("  FINAL: fc={:.4f}, r_eq={:.4f} Bohr, alpha={:.4f} (fsrb2={:.4f})",
                                            params.force_constant, params.equilibrium_distance, params.alpha, fsrb2));
        // Detailed factor breakdown for debugging
        CurcumaLogger::info(fmt::format("  Force Constant Factors: bond_i={:.4f}, bond_j={:.4f}, bstrength={:.4f}",
                                        bond_param_1, bond_param_2, bstrength));
        CurcumaLogger::info(fmt::format("  Corrections: fqq={:.4f}, ringf={:.4f}, fheavy={:.4f}, fpi={:.4f}, fxh={:.4f}, fcn={:.4f}",
                                        fqq, ringf, fheavy, fpi, fxh, fcn));
    }

    return params;
}

GFNFF::GFNFFAngleParams GFNFF::getGFNFFAngleParameters(int atom_i, int atom_j, int atom_k,
                                                        double current_angle, const TopologyInfo& topo_info) const
{
    GFNFFAngleParams params;

    // Original GFN-FF angle parameters from gfnff_param.f90 (angl_angewChem2020 array)
    static const std::vector<double> angle_params = {
        1.661808, 0.300000, 0.018158, 0.029224, 0.572683, // H-B
        0.771055, 1.053577, 2.159889, 1.525582, 0.400000, // C-Ne
        0.041070, 0.028889, 0.086910, 0.494456, 0.409204, // Na-P
        0.864972, 1.986025, 0.491537, 0.050168, 0.072745, // S-Ca
        0.378334, 0.346400, 0.314466, 0.282532, 0.250598, // Sc-Mn
        0.218663, 0.186729, 0.154795, 0.122861, 0.090927, // Fe-Zn
        0.140458, 0.653971, 0.528465, 0.420379, 2.243492, // Ga-Br
        0.400000, 0.035341, 0.022704, 0.195060, 0.188476, // Kr-Zr
        0.181892, 0.175308, 0.168724, 0.162139, 0.155555, // Nb-Rh
        0.148971, 0.142387, 0.135803, 0.169779, 0.265730, // Pd-Sn
        0.505495, 0.398254, 2.640752, 0.568026, 0.032198, // Sb-Cs
        0.036663, 0.281449, 0.280526, 0.279603, 0.278680, // Ba-Nd
        0.277757, 0.276834, 0.275911, 0.274988, 0.274065, // Pm-Tb
        0.273142, 0.272219, 0.271296, 0.270373, 0.269450, // Dy-Yb
        0.268528, 0.267605, 0.253760, 0.239916, 0.226071, // Lu-Re
        0.212227, 0.198382, 0.184538, 0.170693, 0.156849, // Os-Hg
        0.104547, 0.313474, 0.220185, 0.415042, 1.259822, // Tl-At
        0.400000, 0.032198, 0.036663, 0.281449, 0.280526, // Rn-Ra
        0.279603, 0.278680, 0.277757, 0.276834, 0.275911, // Ac-Am
        0.274988, 0.274065, 0.273142, 0.272219, 0.271296, // Cm-Cf
        0.270373, 0.269450, 0.268528 // Es-Lr
    };

    // Claude Generated (Nov 2025): Phase 2 - Neighbor scaling factors for angle force constants
    // CORRECTED: Reference: gfnff_param.f90:247-265 (angl2_angewChem2020 array)
    // These are the CORRECT Fortran values, not the truncated version I had before!
    static const std::vector<double> angl2_neighbors = {
        0.624197, 0.600000, 0.050000, 0.101579, 0.180347, // H-B
        0.755851, 0.761551, 0.813653, 0.791274, 0.400000, // C-Ne
        0.000000, 0.022706, 0.100000, 0.338514, 0.453023, // Na-P
        0.603722, 1.051121, 0.547904, 0.000000, 0.059059, // S-Ca
        0.117040, 0.118438, 0.119836, 0.121234, 0.122632, // Sc-Mn (CORRECTED!)
        0.124031, 0.125429, 0.126827, 0.128225, 0.129623, // Fe-Zn (CORRECTED!)
        0.206779, 0.466678, 0.496442, 0.617321, 0.409933, // Ga-Br (CORRECTED!)
        0.400000, 0.000000, 0.000000, 0.119120, 0.118163, // Kr-Zr (CORRECTED!)
        0.117206, 0.116249, 0.115292, 0.114336, 0.113379, // Nb-Rh (CORRECTED!)
        0.112422, 0.111465, 0.110508, 0.149917, 0.308383, // Pd-Sn (CORRECTED!)
        0.527398, 0.577885, 0.320371, 0.568026, 0.000000, // Sb-Cs (CORRECTED!)
        0.000000, 0.078710, 0.079266, 0.079822, 0.080379, // Ba-Nd (CORRECTED!)
        0.080935, 0.081491, 0.082047, 0.082603, 0.083159, // Pm-Tb (CORRECTED!)
        0.083716, 0.084272, 0.084828, 0.085384, 0.085940, // Dy-Yb (CORRECTED!)
        0.086496, 0.087053, 0.095395, 0.103738, 0.112081, // Lu-Re (CORRECTED!)
        0.120423, 0.128766, 0.137109, 0.145451, 0.153794, // Os-Hg (CORRECTED!)
        0.323570, 0.233450, 0.268137, 0.307481, 0.316447, // Tl-At (CORRECTED!)
        0.400000, 0.000000, 0.000000, 0.119120, 0.118163, // Rn-Ra (CORRECTED!)
        0.117206, 0.116249, 0.115292, 0.114336, 0.113379, // Ac-Am (CORRECTED!)
        0.112422, 0.111465, 0.110508, 0.149917, 0.308383, // Cm-Cf (CORRECTED!)
        0.527398, 0.577885, 0.320371, 0.568026, 0.000000, // Es-Lr (CORRECTED!)
        0.000000, 0.078710, 0.079266, 0.079822, 0.080379  // Last elements
    };

    // Get center atom data
    int z_center = m_atoms[atom_j];

    // Claude Generated (Nov 2025): Determine hybridization of center atom for θ₀ lookup
    // Count neighbors to determine hybridization
    int neighbor_count = 0;
    for (int i = 0; i < m_atomcount; ++i) {
        if (i == atom_j) continue;
        double distance = (m_geometry_bohr.row(atom_j) - m_geometry_bohr.row(i)).norm();
        // Bond threshold: sum of covalent radii * 1.3
        if (distance < 2.0) neighbor_count++;  // Bohr units (~1.06 Å)
    }

    // Assign hybridization based on neighbor count and Z
    int hyb_center = 3;  // Default sp³
    if (neighbor_count <= 1) hyb_center = 1;    // sp (linear/terminal)
    else if (neighbor_count == 2) hyb_center = 2;  // sp² (for C, N, O) or linear for others
    else if (neighbor_count == 3) hyb_center = 3;  // sp³
    else if (neighbor_count >= 5) hyb_center = 5;  // hypervalent

    // Get angle parameters for center atom and neighbors
    int z_i = m_atoms[atom_i];
    int z_k = m_atoms[atom_k];

    double angle_param = (z_center >= 1 && z_center <= static_cast<int>(angle_params.size())) ? angle_params[z_center - 1] : 0.1;
    double angl2_i = (z_i >= 1 && z_i <= static_cast<int>(angl2_neighbors.size())) ? angl2_neighbors[z_i - 1] : 0.1;
    double angl2_k = (z_k >= 1 && z_k <= static_cast<int>(angl2_neighbors.size())) ? angl2_neighbors[z_k - 1] : 0.1;

    // Claude Generated (Nov 2025): Phase 2 implementation - Full angle force constant formula
    // Reference: gfnff_ini.f90:1359-1621
    // Formula: k_ijk = fijk * fqq * f2 * fn * fbsmall * feta

    // Claude Generated (Nov 2025): Phase 2 - Full angle force constant formula infrastructure
    // Reference: gfnff_ini.f90:1359-1621
    // Formula: k_ijk = fijk * fqq * f2 * fn * fbsmall * feta
    //
    // NOTE: The angl2 neighbor arrays have a different meaning in Fortran than simple scaling.
    // In Fortran, angl2 is used with specific logic for parameter assignment, not as direct
    // multipliers. For now, using Phase 1 baseline of 0.01 and deferring full Phase 2b
    // implementation to when element-specific f2 corrections are added.

    // Phase 2: Calculate fijk directly without extra 0.01 scaling
    // The 0.01 in Phase 1 was arbitrary and wrong!
    // Real formula: fijk = angl(center) * angl2(i) * angl2(k)
    double fijk_calc = angle_param * angl2_i * angl2_k;

    // Factor 1: fijk calculation (DEFERRED - needs deeper Fortran analysis)
    // TODO Phase 2b: Proper fijk with angl2 logic from gfnff_param.f90:1359

    // Factor 2: fqq = charge-dependent correction for angles
    // TODO Phase 2b: Re-enable and test fqq after getting fn + fijk working
    // For now: fqq = 1.0 to test other factors first
    double fqq = 1.0;
    // if (atom_i < static_cast<int>(topo_info.eeq_charges.size()) &&
    //     atom_j < static_cast<int>(topo_info.eeq_charges.size()) &&
    //     atom_k < static_cast<int>(topo_info.eeq_charges.size())) {
    //     double q_center = topo_info.eeq_charges[atom_j];
    //     double q_i = topo_info.eeq_charges[atom_i];
    //     double q_k = topo_info.eeq_charges[atom_k];
    //     double qfacBEN = 0.002;
    //     fqq = 1.0 - (q_center * q_i + q_center * q_k) * qfacBEN;
    // }

    // Factor 3: f2 = element-specific correction factor
    // PHASE 4: Implement element-specific f2 for water (gfnff_ini.f90:1486-1491)
    // Water with 2 hydrogens: f2 = 1.20
    // TODO: Complete all element-specific corrections from gfnff_ini.f90:1463-1599
    double f2 = 1.0;  // Default
    // Water case: O with H-O-H angle and both neighbors are H
    if (m_atoms[atom_j] == 8) {  // Central atom is oxygen
        int nh = 0;  // Count hydrogens
        if (m_atoms[atom_i] == 1) nh++;
        if (m_atoms[atom_k] == 1) nh++;
        if (nh == 2) f2 = 1.20;  // H2O specific: f2 = 1.20 (gfnff_ini.f90:1491)
    }

    // Factor 4: fn = coordination number dependence
    // Reference: gfnff_ini.f90:1612
    // Fortran: fn = 1.0d0 - 2.36d0 / dble(nn)**2  where nn = neighbor count
    int neighbor_count_angle = 0;
    for (int i = 0; i < m_atomcount; ++i) {
        if (i == atom_j) continue;
        double distance = (m_geometry_bohr.row(atom_j) - m_geometry_bohr.row(i)).norm();
        if (distance < 2.0) neighbor_count_angle++;  // Bond threshold
    }
    // Factor 4: fn - coordination number dependence
    // PHASE 4: Use REAL coordination number from D3-style calculation
    // Formula (Fortran gfnff_ini.f90:1612): fn = 1.0 - 2.36 / nn²
    // where nn = topo%nb(20,i) = Coordination number of central atom
    // Reference: gfnff_ini.f90:1377,1612
    const double threshold_cn_squared = 40.0 * 40.0;  // ~40 Bohr cutoff (standard GFN-FF)
    Vector coord_numbers = calculateCoordinationNumbers(threshold_cn_squared);
    double nn = static_cast<int>(std::round(coord_numbers[atom_j]));  // Central atom coordination number
    nn = std::max(1.0, nn);  // Ensure nn >= 1
    double fn = 1.0 - 2.36 / (nn * nn);
    fn = std::max(0.05, fn);  // Ensure fn stays positive

    // Factor 5: fbsmall = small-angle correction
    // PHASE 4: Implement small-angle correction from gfnff_ini.f90:1618
    // Formula: fbsmall = 1.0 - fbs1 * exp(-0.64*(theta - π)²)
    // where theta is the equilibrium angle in RADIANS, fbs1 is a parameter
    // Reference: gfnff_param.f90:754 gen%fbs1 = 0.50
    const double pi = 3.14159265358979323846;
    const double fbs1 = 0.50;  // GFN-FF parameter from gfnff_param.f90:754
    double theta_eq_rad = params.equilibrium_angle;  // Already in radians from topology
    double fbsmall = 1.0 - fbs1 * std::exp(-0.64 * (theta_eq_rad - pi) * (theta_eq_rad - pi));

    // Factor 6: feta = metal eta-coordination correction
    // TODO Phase 2.5: Implement feta correction for metals
    double feta = 1.0;

    // PHASE 2D: Test hypothesis - angl2 is only for FILTERING, not force constant!
    // Discovery: angl2 may only be used for threshold check (skip angles)
    // The force constant is purely from angle_param with 0.01 scaling + correction factors

    double fijk = fijk_calc;  // Only for threshold check: angle_param * angl2_i * angl2_k

    // Check threshold: if fijk is too small, skip this angle
    static const double THRESHOLD = 0.001;  // gen%fcthr in Fortran
    if (fijk < THRESHOLD) {
        // Skip angles with too small fijk - this matches Fortran filtering!
        params.force_constant = 0.0;
    } else {
        // PHASE 4: CORRECT Force Constant Formula from gfnff_ini.f90:1621
        // Previous WRONG: base_k = angle_param * 0.01; corrected_k = base_k * corrections
        // CORRECT: fc = fijk * fqq * f2 * fn * fbsmall * feta
        // where fijk = angle_param * angl2_i * angl2_k
        params.force_constant = fijk * fqq * f2 * fn * fbsmall * feta;
    }

    // ===========================================================================================
    // CRITICAL FIX: Equilibrium angle from TOPOLOGY (hybridization), NOT current geometry!
    // ===========================================================================================
    // Reference: gfnff_ini.f90:1441-1617 shows θ₀ is topology-dependent:
    // sp=180°, sp²=120°, sp³=109.5°, hypervalent=90°
    // This is NOT dependent on current geometry angle.
    // Using current_angle makes restoring force = 0, breaking optimization!

    double r0_deg = 100.0;  // Fallback default

    // Determine base equilibrium angle from hybridization
    // Reference: gfnff_ini.f90:1441-1617
    switch (hyb_center) {
        case 1:   r0_deg = 180.0;   // sp   - linear
            break;
        case 2:   r0_deg = 120.0;   // sp²  - trigonal planar
            break;
        case 3:   r0_deg = 109.5;   // sp³  - tetrahedral
            break;
        case 5:   r0_deg = 90.0;    // hypervalent - square planar
            break;
        default:  r0_deg = 100.0;   // Fallback
    }

    // TODO (Phase 2+): Element-specific corrections
    // Reference: gfnff_ini.f90:1486-1599
    // - Water: θ₀=100°, f2=1.20 for H-O-H
    // - Aromatic: θ₀=120° for aromatic C and Si
    // - Ring corrections: 3-ring=82°, 4-ring=96°, etc.
    // - Metal coordination: Special cases for transition metals

    // Convert to radians
    params.equilibrium_angle = r0_deg * M_PI / 180.0;

    return params;
}

bool GFNFF::loadGFNFFCharges()
{
    // Try to load charges from reference GFN-FF calculation
    std::string charges_file = "releaseX/gfnff_charges";
    std::ifstream file(charges_file);

    if (!file.is_open()) {
        std::cerr << "Warning: Could not open " << charges_file << " for reading charges" << std::endl;
        return false;
    }

    m_charges = Vector::Zero(m_atomcount);

    std::string line;
    int atom_idx = 0;

    while (std::getline(file, line) && atom_idx < m_atomcount) {
        // Remove leading/trailing whitespace
        line.erase(0, line.find_first_not_of(" \t"));
        line.erase(line.find_last_not_of(" \t") + 1);

        if (!line.empty()) {
            try {
                double charge = std::stod(line);
                m_charges[atom_idx] = charge;
                atom_idx++;
            } catch (const std::exception& e) {
                std::cerr << "Error parsing charge on line " << atom_idx + 1 << ": " << e.what() << std::endl;
                return false;
            }
        }
    }

    file.close();

    if (atom_idx != m_atomcount) {
        std::cerr << "Warning: Expected " << m_atomcount << " charges, got " << atom_idx << std::endl;
        return false;
    }

    std::cout << "Loaded " << atom_idx << " GFN-FF charges from " << charges_file << std::endl;
    return true;
}

// =================================================================================
// Advanced GFN-FF Parameter Generation (Placeholder implementations)
// =================================================================================

Vector GFNFF::calculateCoordinationNumbers(double threshold) const
{
    // GFN-FF Coordination Number Calculation using D3-style Error Function
    // Reference: external/gfnff/src/gfnff_cn.f90:66-126
    //
    // This implements the create_erfCN and create_logCN functions from Fortran:
    // 1. Calculate raw CN using error function: erfCN = 0.5 * (1 + erf(kn * dr))
    // 2. Apply logarithmic transformation: logCN = log(1+e^cnmax) - log(1+e^(cnmax-cn))
    //
    // NOTE: This replaces the old create_expCN(16.0) formula that was commented out
    // in the Fortran reference (gfnff_cn.f90:88)

    // D3 coordination number parameters (gfnff_cn.f90:66, gfnff_param.f90:462)
    const double kn = -7.5;      // Error function steepness parameter
    const double cnmax = 4.4;    // Maximum coordination number cutoff

    Vector cn = Vector::Zero(m_atomcount);

    // Step 1: Calculate raw coordination numbers using error function
    for (int i = 0; i < m_atomcount; ++i) {
        double cn_i = 0.0;
        for (int j = 0; j < m_atomcount; ++j) {
            if (i == j)
                continue;

            Vector ri = m_geometry_bohr.row(i);
            Vector rj = m_geometry_bohr.row(j);
            double distance = (ri - rj).norm();

            // Distance threshold check (Fortran: thr = sqrt(thr2), typically ~40 Bohr)
            if (distance > threshold)
                continue;

            double rcov_i = getCovalentRadius(m_atoms[i]);
            double rcov_j = getCovalentRadius(m_atoms[j]);
            double r_cov = rcov_i + rcov_j;

            // D3-style error function CN contribution
            // Formula: erfCN = 0.5 * (1 + erf(kn * dr))
            // where dr = (r - r0) / r0
            double dr = (distance - r_cov) / r_cov;
            double erfCN = 0.5 * (1.0 + std::erf(kn * dr));

            cn_i += erfCN;
        }

        // Step 2: Apply logarithmic transformation (Fortran create_logCN)
        // logCN = log(1 + e^cnmax) - log(1 + e^(cnmax - cn))
        // This smoothly caps CN at cnmax and provides better numerical behavior
        cn[i] = std::log(1.0 + std::exp(cnmax)) - std::log(1.0 + std::exp(cnmax - cn_i));
    }

    return cn;
}

std::vector<Matrix> GFNFF::calculateCoordinationNumberDerivatives(const Vector& cn, double threshold) const
{
    // GFN-FF Coordination Number Derivatives using D3-style Error Function
    // Reference: external/gfnff/src/gfnff_cn.f90:94-117
    //
    // Mathematical formulation (two-step chain rule):
    // 1. Raw CN: cn_raw_i = Σ_j erfCN(r_ij) where erfCN = 0.5 * (1 + erf(kn * dr))
    // 2. Log CN: logCN_i = log(1 + e^cnmax) - log(1 + e^(cnmax - cn_raw))
    //
    // Derivatives:
    // ∂logCN_i/∂r_ij = (∂logCN/∂cn_raw) * (∂erfCN/∂r_ij)
    //
    // where:
    // - ∂logCN/∂cn_raw = e^cnmax / (e^cnmax + e^cn_raw)  [create_dlogCN]
    // - ∂erfCN/∂r = (kn / sqrt(π)) * exp(-kn² * dr²) / r0  [create_derfCN]
    //
    // NOTE: The input 'cn' vector already contains logCN values from calculateCoordinationNumbers()
    // We need to invert the log transformation to get raw CN for the derivative calculation.

    // D3 coordination number parameters (matching calculateCoordinationNumbers)
    const double kn = -7.5;         // Error function steepness parameter
    const double cnmax = 4.4;       // Maximum coordination number cutoff
    const double sqrtpi = 1.77245385091;  // sqrt(π) for derivative formula

    // Initialize 3D tensor as vector of matrices
    // dcn[0] = ∂logCN/∂x, dcn[1] = ∂logCN/∂y, dcn[2] = ∂logCN/∂z
    std::vector<Matrix> dcn(3, Matrix::Zero(m_atomcount, m_atomcount));

    // Step 1: Compute raw CN from input logCN by inverting the transformation
    // logCN = log(1 + e^cnmax) - log(1 + e^(cnmax - cn_raw))
    // Solving for cn_raw is complex, so we'll recalculate raw CN
    Vector cn_raw = Vector::Zero(m_atomcount);
    for (int i = 0; i < m_atomcount; ++i) {
        double cn_i = 0.0;
        for (int j = 0; j < m_atomcount; ++j) {
            if (i == j) continue;

            Vector ri = m_geometry_bohr.row(i);
            Vector rj = m_geometry_bohr.row(j);
            double distance = (ri - rj).norm();

            if (distance > threshold) continue;

            double rcov_i = getCovalentRadius(m_atoms[i]);
            double rcov_j = getCovalentRadius(m_atoms[j]);
            double r_cov = rcov_i + rcov_j;

            double dr = (distance - r_cov) / r_cov;
            double erfCN = 0.5 * (1.0 + std::erf(kn * dr));
            cn_i += erfCN;
        }
        cn_raw[i] = cn_i;
    }

    // Step 2: Calculate dlogCN/dcn for each atom (Fortran create_dlogCN)
    Vector dlogdcn = Vector::Zero(m_atomcount);
    for (int i = 0; i < m_atomcount; ++i) {
        // dlogCN/dcn = e^cnmax / (e^cnmax + e^cn_raw)
        dlogdcn[i] = std::exp(cnmax) / (std::exp(cnmax) + std::exp(cn_raw[i]));
    }

    // Step 3: Calculate derivatives for all atom pairs
    for (int i = 0; i < m_atomcount; ++i) {
        Vector ri = m_geometry_bohr.row(i);
        double rcov_i = getCovalentRadius(m_atoms[i]);
        double dlogdcn_i = dlogdcn[i];

        for (int j = 0; j < i; ++j) {  // Only j < i to avoid double counting
            Vector rj = m_geometry_bohr.row(j);
            double rcov_j = getCovalentRadius(m_atoms[j]);

            // Distance and direction
            Vector r_ij_vec = rj - ri;
            double r_ij = r_ij_vec.norm();

            if (r_ij > threshold) continue;

            double rcov_sum = rcov_i + rcov_j;

            // Error function CN derivative (Fortran create_derfCN)
            // derfCN/dr = (kn / sqrt(π)) * exp(-kn² * dr²) / r0
            double dr = (r_ij - rcov_sum) / rcov_sum;
            double derfCN_dr = (kn / sqrtpi) * std::exp(-kn * kn * dr * dr) / rcov_sum;

            // Unit vector from i to j
            Vector grad_direction = r_ij_vec / r_ij;

            // Chain rule: d(logCN)/dr = (dlogCN/dcn) * (derfCN/dr)
            double dlogdcn_j = dlogdcn[j];

            // Apply derivatives (Fortran lines 110-115)
            // Note: Both CN_i and CN_j are affected by r_ij
            for (int dim = 0; dim < 3; ++dim) {
                double rij_component = derfCN_dr * grad_direction[dim];

                // ∂logCN_j/∂r_j and ∂logCN_i/∂r_i (diagonal terms)
                dcn[dim](j, j) += dlogdcn_j * rij_component;
                dcn[dim](i, i) -= dlogdcn_i * rij_component;

                // ∂logCN_j/∂r_i and ∂logCN_i/∂r_j (off-diagonal terms)
                dcn[dim](j, i) = -dlogdcn_j * rij_component;
                dcn[dim](i, j) = dlogdcn_i * rij_component;
            }
        }
    }

    return dcn;
}

std::vector<int> GFNFF::determineHybridization() const
{
    // Phase 2.3: Enhanced hybridization detection (geometry-based)
    // Reference: gfnff_ini.f90:400-550 (hybridization assignment)
    // Analyzes bond angles and geometry, not just neighbor count

    std::vector<int> hyb(m_atomcount, 3); // Default to sp3
    double bond_threshold = 1.3;

    for (int i = 0; i < m_atomcount; ++i) {
        int z = m_atoms[i];
        Vector ri = m_geometry_bohr.row(i);

        // Step 1: Find bonded neighbors
        std::vector<int> neighbors;
        std::vector<Vector> bond_vectors;

        for (int j = 0; j < m_atomcount; ++j) {
            if (i == j) continue;

            Vector rj = m_geometry_bohr.row(j);
            double distance = (ri - rj).norm();
            double rcov_sum = getCovalentRadius(z) + getCovalentRadius(m_atoms[j]);

            if (distance < bond_threshold * rcov_sum) {
                neighbors.push_back(j);
                bond_vectors.push_back((rj - ri).normalized());
            }
        }

        int neighbor_count = neighbors.size();

        // Step 2: Geometry-based hybridization assignment
        if (neighbor_count == 0 || neighbor_count == 1) {
            hyb[i] = 1; // sp (terminal or diatomic)

        } else if (neighbor_count == 2) {
            // Check if linear (sp) or bent (sp2)
            double dot_product = bond_vectors[0].dot(bond_vectors[1]);
            double angle = std::acos(std::max(-1.0, std::min(1.0, dot_product)));

            if (angle > 2.8) { // ~160° - linear geometry
                hyb[i] = 1; // sp
            } else {
                hyb[i] = 2; // sp2 (bent, like H2O oxygen)
            }

        } else if (neighbor_count == 3) {
            // Check if planar (sp2) or pyramidal (sp3)
            // Calculate sum of bond angles (should be ~360° for planar)
            double angle_sum = 0.0;
            for (int j = 0; j < 3; ++j) {
                int k = (j + 1) % 3;
                double dot = bond_vectors[j].dot(bond_vectors[k]);
                angle_sum += std::acos(std::max(-1.0, std::min(1.0, dot)));
            }

            // Planar: sum ≈ 2π, Pyramidal: sum < 2π
            if (angle_sum > 6.0) { // ~345° - nearly planar
                hyb[i] = 2; // sp2
            } else {
                hyb[i] = 3; // sp3
            }

        } else if (neighbor_count >= 4) {
            hyb[i] = 3; // sp3 (tetrahedral or higher coordination)

            // Special case: sp3d, sp3d2 for main group elements (P, S, etc.)
            if ((z == 15 || z == 16 || z == 17) && neighbor_count >= 5) {
                hyb[i] = 3; // Still mark as sp3 for GFN-FF purposes
            }
        }
    }

    return hyb;
}

std::vector<int> GFNFF::detectPiSystems(const std::vector<int>& hyb) const
{
    // Phase 2.2: Pi-system detection (conjugated fragments)
    // Reference: gfnff_ini.f90:1100-1300 (pi-system setup)
    // Identifies conjugated chains/rings and marks aromatic systems

    std::vector<int> pi_fragments(m_atomcount, 0); // 0 = not in pi-system
    double bond_threshold = 1.3;

    // Step 1: Find all sp2 and sp atoms (potential pi-system members)
    std::vector<bool> is_pi_atom(m_atomcount, false);
    for (int i = 0; i < m_atomcount; ++i) {
        if (hyb[i] == 1 || hyb[i] == 2) { // sp or sp2
            is_pi_atom[i] = true;
        }
    }

    // Step 2: Build adjacency for pi-atoms only
    std::vector<std::vector<int>> pi_neighbors(m_atomcount);
    for (int i = 0; i < m_atomcount; ++i) {
        if (!is_pi_atom[i]) continue;

        for (int j = i + 1; j < m_atomcount; ++j) {
            if (!is_pi_atom[j]) continue;

            double distance = (m_geometry_bohr.row(i) - m_geometry_bohr.row(j)).norm();
            double rcov_sum = getCovalentRadius(m_atoms[i]) + getCovalentRadius(m_atoms[j]);

            if (distance < bond_threshold * rcov_sum) {
                pi_neighbors[i].push_back(j);
                pi_neighbors[j].push_back(i);
            }
        }
    }

    // Step 3: Connected component analysis (union-find/DFS)
    // Assign same fragment ID to all connected pi-atoms
    std::vector<bool> visited(m_atomcount, false);
    int fragment_id = 1;

    for (int start = 0; start < m_atomcount; ++start) {
        if (!is_pi_atom[start] || visited[start]) continue;

        // DFS to mark all atoms in this conjugated system
        std::stack<int> stack;
        stack.push(start);
        visited[start] = true;

        while (!stack.empty()) {
            int current = stack.top();
            stack.pop();
            pi_fragments[current] = fragment_id;

            for (int neighbor : pi_neighbors[current]) {
                if (!visited[neighbor]) {
                    visited[neighbor] = true;
                    stack.push(neighbor);
                }
            }
        }

        fragment_id++;
    }

    return pi_fragments;
}

std::vector<int> GFNFF::findSmallestRings() const
{
    // Phase 2.1: Ring detection algorithm (DFS-based)
    // Finds the smallest ring each atom belongs to (3-8 membered rings)
    // Reference: gfnff_helpers.f90:99-316 (getring36 subroutine)
    // Educational simplification: DFS instead of exhaustive enumeration

    std::vector<int> ring_sizes(m_atomcount, 0); // 0 = not in ring

    // Step 1: Build adjacency list from bonds
    std::vector<std::vector<int>> neighbors(m_atomcount);
    double bond_threshold = 1.3; // Same as bond detection

    for (int i = 0; i < m_atomcount; ++i) {
        for (int j = i + 1; j < m_atomcount; ++j) {
            double distance = (m_geometry_bohr.row(i) - m_geometry_bohr.row(j)).norm();
            double rcov_sum = getCovalentRadius(m_atoms[i]) + getCovalentRadius(m_atoms[j]);

            if (distance < bond_threshold * rcov_sum) {
                neighbors[i].push_back(j);
                neighbors[j].push_back(i);
            }
        }
    }

    // Step 2: For each atom, find shortest cycle using BFS
    // (BFS guarantees shortest path = smallest ring)
    for (int start_atom = 0; start_atom < m_atomcount; ++start_atom) {
        std::vector<int> distance(m_atomcount, -1);  // -1 = not visited
        std::vector<int> parent(m_atomcount, -1);
        std::queue<int> queue;

        distance[start_atom] = 0;
        queue.push(start_atom);

        int smallest_ring = 0;

        while (!queue.empty() && smallest_ring == 0) {
            int current = queue.front();
            queue.pop();

            for (int neighbor : neighbors[current]) {
                if (distance[neighbor] == -1) {
                    // Unvisited neighbor: explore
                    distance[neighbor] = distance[current] + 1;
                    parent[neighbor] = current;
                    queue.push(neighbor);

                } else if (parent[current] != neighbor && distance[neighbor] <= distance[current]) {
                    // Found cycle! Ring size = distance[current] + distance[neighbor] + 1
                    int ring_size = distance[current] + distance[neighbor] + 1;

                    // Only consider rings up to size 8 (common in organic chemistry)
                    if (ring_size >= 3 && ring_size <= 8) {
                        smallest_ring = ring_size;
                        break; // BFS ensures this is the smallest
                    }
                }
            }
        }

        ring_sizes[start_atom] = smallest_ring;
    }

    return ring_sizes;
}

Vector GFNFF::calculateEEQCharges(const Vector& cn, const std::vector<int>& hyb, const std::vector<int>& rings) const
{
    // Phase 3: Full EEQ (Electronegativity Equalization) implementation
    // Reference: gfnff_engrad.F90:1274-1391 (goed_gfnff subroutine)
    //
    // Solves linear system: A·q = b
    // where A = hardness matrix + Coulomb interaction + constraint
    //       b = -electronegativity - CN corrections
    //       q = atomic charges (unknown)
    //
    // Mathematical formulation:
    // A(i,i) = γ_i + sqrt(2π)/sqrt(α_i)  (self-interaction)
    // A(i,j) = erf(γ_ij * r_ij) / r_ij    (Coulomb interaction with damping)
    // A(n+1, :) = 1 (charge constraint: Σq = total_charge)
    //
    // Literature: S. Spicher, S. Grimme, Angew. Chem. Int. Ed. 2020, 59, 15665-15673

    const int n = m_atomcount;
    const int m = n + 1;  // n atoms + 1 charge constraint

    // Constants
    const double sqrt_2pi = 0.79788456080287;  // sqrt(2/π)

    // Setup EEQ matrix A and RHS vector b
    Matrix A = Matrix::Zero(m, m);
    Vector b = Vector::Zero(m);

    // Phase 3.1: Build diagonal and off-diagonal elements
    for (int i = 0; i < n; ++i) {
        EEQParameters params_i = getEEQParameters(m_atoms[i]);

        // Diagonal: A(i,i) = gamma_i + sqrt(2π)/sqrt(alpha_i)
        // This is the self-interaction (chemical hardness + Coulomb self-energy)
        A(i, i) = params_i.gam + sqrt_2pi / std::sqrt(params_i.alp);

        // RHS: b(i) = -chi_i - cnf_i * sqrt(CN_i)
        // Electronegativity with CN correction
        b(i) = -params_i.chi - params_i.cnf * std::sqrt(cn[i]);

        // Off-diagonal: Coulomb interaction with error function damping
        for (int j = 0; j < i; ++j) {
            EEQParameters params_j = getEEQParameters(m_atoms[j]);

            // Distance between atoms
            Vector r_ij_vec = m_geometry_bohr.row(i) - m_geometry_bohr.row(j);
            double r_ij = r_ij_vec.norm();

            // Damping parameter: γ_ij = 1/sqrt(α_i + α_j)
            double gamma_ij = 1.0 / std::sqrt(params_i.alp + params_j.alp);

            // Coulomb interaction with error function damping
            // A(i,j) = erf(γ_ij * r_ij) / r_ij
            double erf_arg = gamma_ij * r_ij;
            double erf_val = std::erf(erf_arg);
            double coulomb = erf_val / r_ij;

            A(i, j) = coulomb;
            A(j, i) = coulomb;  // Symmetric matrix
        }
    }

    // Phase 3.2: Add charge constraint (Lagrange multiplier method)
    // Σq_i = total_charge
    // Last row and column enforce this constraint
    for (int i = 0; i < n; ++i) {
        A(n, i) = 1.0;  // ∂(Σq)/∂q_i = 1
        A(i, n) = 1.0;  // Symmetric
    }
    b(n) = m_charge;  // Total charge on system

    // Phase 3.3: Solve linear system A·q = b
    // Use Eigen's LDLT decomposition (symmetric indefinite solver)
    // This handles the constraint properly via the Lagrange multiplier
    Vector q_extended = A.ldlt().solve(b);

    // Extract charges (first n elements, last element is Lagrange multiplier)
    Vector charges = q_extended.head(n);

    // Phase 3.4: Safety check - verify charge conservation
    double total_charge_actual = charges.sum();
    double charge_error = std::abs(total_charge_actual - m_charge);
    if (charge_error > 1e-6) {
        std::cerr << "Warning: EEQ charge constraint not satisfied. "
                  << "Expected: " << m_charge << ", Got: " << total_charge_actual
                  << ", Error: " << charge_error << std::endl;
    }

    // Claude Generated (2025): Debug EEQ charges
    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("=== EEQ Charges Calculated ===");
        for (int i = 0; i < n; ++i) {
            EEQParameters params = getEEQParameters(m_atoms[i]);
            CurcumaLogger::param(fmt::format("atom_{}_Z{}", i, m_atoms[i]),
                fmt::format("q={:.6f}, χ={:.6f}, γ={:.6f}, α={:.6f}, CN={:.3f}",
                    charges[i], params.chi, params.gam, params.alp, cn[i]));
        }
        CurcumaLogger::param("total_charge", fmt::format("{:.6f}", total_charge_actual));
    }

    return charges;
}

double GFNFF::calculateEEQEnergy(const Vector& charges, const Vector& cn) const
{
    // Phase 3.2: EEQ electrostatic energy calculation
    // Reference: gfnff_engrad.F90:1378-1389 (EEQ energy in goed_gfnff)
    //
    // Energy formula:
    // E_EEQ = Σ_i<j q_i*q_j*erf(γ_ij*r_ij)/r_ij
    //       + Σ_i [-q_i*(χ_i + cnf_i*√CN_i) + 0.5*q_i²*(γ_i + √(2π)/√α_i)]
    //
    // This is the electrostatic energy from the EEQ method.
    // Note: In GFN-FF, this uses "frozen charge" approximation for gradients.

    const int n = m_atomcount;
    const double sqrt_2pi = 0.79788456080287;  // sqrt(2/π)
    double energy = 0.0;

    // Pairwise Coulomb interactions: Σ_i<j q_i*q_j*γ_ij(r_ij)
    for (int i = 0; i < n; ++i) {
        EEQParameters params_i = getEEQParameters(m_atoms[i]);

        for (int j = 0; j < i; ++j) {
            EEQParameters params_j = getEEQParameters(m_atoms[j]);

            // Distance between atoms
            Vector r_ij_vec = m_geometry_bohr.row(i) - m_geometry_bohr.row(j);
            double r_ij = r_ij_vec.norm();

            // Damping parameter: γ_ij = 1/sqrt(α_i + α_j)
            double gamma_ij = 1.0 / std::sqrt(params_i.alp + params_j.alp);

            // Coulomb interaction with error function damping
            double erf_arg = gamma_ij * r_ij;
            double erf_val = std::erf(erf_arg);
            double coulomb = erf_val / r_ij;

            // Pairwise energy contribution
            energy += charges[i] * charges[j] * coulomb;
        }

        // Self-energy and electronegativity terms for atom i
        // E_i = -q_i*(χ_i + cnf_i*√CN_i) + 0.5*q_i²*(γ_i + √(2π)/√α_i)
        double self_interaction = params_i.gam + sqrt_2pi / std::sqrt(params_i.alp);
        double en_term = params_i.chi + params_i.cnf * std::sqrt(cn[i]);

        energy += -charges[i] * en_term + 0.5 * charges[i] * charges[i] * self_interaction;
    }

    return energy;
}

json GFNFF::generateTopologyAwareBonds(const Vector& cn, const std::vector<int>& hyb,
    const Vector& charges, const std::vector<int>& rings) const
{
    json bonds = json::array();

    // Phase 9: Create TopologyInfo structure from separate parameters
    TopologyInfo topo_info;
    topo_info.coordination_numbers = cn;
    topo_info.hybridization = hyb;
    topo_info.eeq_charges = charges;
    topo_info.ring_sizes = rings;
    // pi_fragments, is_metal, is_aromatic not used in getGFNFFBondParameters

    // Phase 2: Topology-aware bond parameter generation
    // Start with basic GFN-FF bond detection
    double bond_threshold = 1.3; // Factor for covalent radii sum

    for (int i = 0; i < m_atomcount; ++i) {
        for (int j = i + 1; j < m_atomcount; ++j) {
            Vector ri = m_geometry_bohr.row(i);
            Vector rj = m_geometry_bohr.row(j);
            double distance = (ri - rj).norm();

            // Get covalent radii for atoms i and j
            double rcov_i = getCovalentRadius(m_atoms[i]);
            double rcov_j = getCovalentRadius(m_atoms[j]);

            if (distance < bond_threshold * (rcov_i + rcov_j)) {
                // This is a bond - generate topology-aware GFN-FF parameters
                json bond;
                bond["type"] = 3; // GFN-FF type
                bond["i"] = i;
                bond["j"] = j;
                bond["k"] = 0;
                bond["distance"] = distance;

                // Phase 9: Get bond parameters with full topology awareness
                auto bond_params = getGFNFFBondParameters(i, j, m_atoms[i], m_atoms[j], distance, topo_info);

                // Phase 2: Apply topology corrections to force constant
                double topology_factor = 1.0;

                // Ring strain correction (small rings are stiffer)
                int ring_i = rings[i];
                int ring_j = rings[j];
                if (ring_i > 0 && ring_j > 0) {
                    // Both atoms in rings - assume they're in the same ring if bonded
                    int ring_size = std::min(ring_i, ring_j);
                    if (ring_size == 3) {
                        topology_factor *= 1.25; // Cyclopropane +25% strain
                    } else if (ring_size == 4) {
                        topology_factor *= 1.15; // Cyclobutane +15% strain
                    } else if (ring_size == 5) {
                        topology_factor *= 1.05; // Cyclopentane slight strain
                    }
                    // 6+ membered rings have normal parameters
                }

                // Pi-system correction (conjugated bonds are stiffer)
                int pi_i = static_cast<int>(cn[i]); // Using CN as proxy for pi_fragments
                int pi_j = static_cast<int>(cn[j]);
                if (pi_i > 0 && pi_j > 0 && (hyb[i] == 2 || hyb[i] == 1) && (hyb[j] == 2 || hyb[j] == 1)) {
                    // Both atoms are sp2/sp and in conjugated system
                    topology_factor *= 1.15; // Conjugated bonds +15%
                }

                // Apply topology corrections
                bond_params.force_constant *= topology_factor;

                bond["fc"] = bond_params.force_constant;
                bond["r0_ij"] = bond_params.equilibrium_distance;
                bond["r0_ik"] = 0.0;
                bond["exponent"] = bond_params.alpha;

                bonds.push_back(bond);
            }
        }
    }

    return bonds;
}

json GFNFF::generateTopologyAwareAngles(const Vector& cn, const std::vector<int>& hyb,
    const Vector& charges, const std::vector<int>& rings) const
{
    json angles = json::array();

    // Phase 2: Topology-aware angle parameter generation
    // Build bond list first
    std::vector<std::pair<int, int>> bond_list;
    json bonds = generateTopologyAwareBonds(cn, hyb, charges, rings);

    for (const auto& bond : bonds) {
        bond_list.push_back({ bond["i"], bond["j"] });
    }

    // Generate angles from bonded topology
    for (int center = 0; center < m_atomcount; ++center) {
        std::vector<int> neighbors;

        // Find all atoms bonded to center
        for (const auto& bond : bond_list) {
            if (bond.first == center)
                neighbors.push_back(bond.second);
            if (bond.second == center)
                neighbors.push_back(bond.first);
        }

        // Generate all possible angles with center as middle atom
        for (size_t i = 0; i < neighbors.size(); ++i) {
            for (size_t j = i + 1; j < neighbors.size(); ++j) {
                json angle;
                angle["type"] = 3; // GFN-FF type
                angle["i"] = neighbors[i];
                angle["j"] = center;
                angle["k"] = neighbors[j];

                // Calculate current angle
                Vector ri = m_geometry_bohr.row(neighbors[i]);
                Vector rj = m_geometry_bohr.row(center);
                Vector rk = m_geometry_bohr.row(neighbors[j]);

                Vector v1 = ri - rj;
                Vector v2 = rk - rj;

                double v1_norm = v1.norm();
                double v2_norm = v2.norm();

                // Skip if vectors are too small
                if (v1_norm < 1e-10 || v2_norm < 1e-10) {
                    continue;
                }

                double cos_angle = v1.dot(v2) / (v1_norm * v2_norm);
                cos_angle = std::max(-1.0, std::min(1.0, cos_angle));
                double current_angle = acos(cos_angle);

                // Get basic angle parameters
                // Note: generateTopologyAwareAngles is called with raw charges vector
                // Create a minimal TopologyInfo for compatibility with Phase 2
                TopologyInfo topo_compat;
                topo_compat.eeq_charges = charges;  // Use provided charges
                topo_compat.coordination_numbers = cn;  // Use provided CN
                topo_compat.hybridization = hyb;  // Use provided hybridization

                auto angle_params = getGFNFFAngleParameters(m_atoms[neighbors[i]],
                    m_atoms[center],
                    m_atoms[neighbors[j]],
                    current_angle,
                    topo_compat);

                // Phase 2: Apply topology corrections to force constant
                double topology_factor = 1.0;

                // Ring strain correction (small ring angles are stiffer)
                int ring_center = rings[center];
                int ring_i = rings[neighbors[i]];
                int ring_k = rings[neighbors[j]];

                // If all three atoms are in rings, assume they form a ring angle
                if (ring_center > 0 && ring_i > 0 && ring_k > 0) {
                    int ring_size = std::min({ring_center, ring_i, ring_k});
                    if (ring_size == 3) {
                        topology_factor *= 1.30; // Cyclopropane +30% angle strain
                    } else if (ring_size == 4) {
                        topology_factor *= 1.20; // Cyclobutane +20% angle strain
                    } else if (ring_size == 5) {
                        topology_factor *= 1.08; // Cyclopentane +8% strain
                    }
                }

                // Hybridization correction (linear/planar geometries)
                int hyb_center = hyb[center];
                if (hyb_center == 1) {
                    // sp center - expects linear geometry (180°)
                    topology_factor *= 1.25; // Stiffer angle bending for sp
                } else if (hyb_center == 2) {
                    // sp2 center - expects planar geometry (120°)
                    topology_factor *= 1.10; // Moderate stiffness for sp2
                }
                // sp3 uses default parameters

                // Apply topology corrections
                angle_params.force_constant *= topology_factor;

                angle["fc"] = angle_params.force_constant;
                angle["theta0_ijk"] = angle_params.equilibrium_angle;
                angle["r0_ij"] = v1_norm; // Distance i-j
                angle["r0_ik"] = v2_norm; // Distance k-j

                angles.push_back(angle);
            }
        }
    }

    return angles;
}

json GFNFF::detectHydrogenBonds(const Vector& charges) const
{
    json hbonds = json::array();

    // TODO: Implement hydrogen bond detection
    // This should identify A-H...B interactions
    // Based on gfnff_hbset functions

    return hbonds;
}

GFNFF::TopologyInfo GFNFF::calculateTopologyInfo() const
{
    TopologyInfo topo_info;

    // Calculate all topology information (Phase 2 implementations)
    topo_info.coordination_numbers = calculateCoordinationNumbers();
    topo_info.hybridization = determineHybridization();         // Phase 2.3 ✅
    topo_info.pi_fragments = detectPiSystems(topo_info.hybridization);  // Phase 2.2 ✅
    topo_info.ring_sizes = findSmallestRings();                 // Phase 2.1 ✅
    topo_info.eeq_charges = calculateEEQCharges(topo_info.coordination_numbers,
        topo_info.hybridization,
        topo_info.ring_sizes);

    // Initialize metal and aromatic flags
    topo_info.is_metal.resize(m_atomcount, false);
    topo_info.is_aromatic.resize(m_atomcount, false);

    for (int i = 0; i < m_atomcount; ++i) {
        int z = m_atoms[i];

        // Mark metals (simplified: transition metals and lanthanides/actinides)
        if ((z >= 21 && z <= 30) || (z >= 39 && z <= 48) ||
            (z >= 57 && z <= 80) || (z >= 89 && z <= 103)) {
            topo_info.is_metal[i] = true;
        }

        // Phase 2.2: Aromaticity detection (Hückel 4n+2 rule)
        // Simplified: 6-membered rings in pi-systems are aromatic
        if (topo_info.ring_sizes[i] == 6 && topo_info.pi_fragments[i] > 0 && topo_info.hybridization[i] == 2) {
            topo_info.is_aromatic[i] = true;
        }
        // 5-membered rings with heteroatoms (pyrrole, furan) are also aromatic
        else if (topo_info.ring_sizes[i] == 5 && topo_info.pi_fragments[i] > 0) {
            if ((z == 7 || z == 8 || z == 16) && topo_info.hybridization[i] == 2) {
                topo_info.is_aromatic[i] = true;
            }
        }
    }

    return topo_info;
}

GFNFF::EEQParameters GFNFF::getEEQParameters(int atom_idx, const TopologyInfo& topo_info) const
{
    EEQParameters params;

    int z = m_atoms[atom_idx];

    // EEQ parameters from gfnff_param.f90 (chi_angewChem2020, gam_angewChem2020)
    static const std::vector<double> chi_angewChem2020 = {
        1.227054, 1.451412, 0.813363, 1.062841, 1.186499, // H-B
        1.311555, 1.528485, 1.691201, 1.456784, 1.231037, // C-Ne
        0.772989, 1.199092, 1.221576, 1.245964, 1.248942, // Na-P
        1.301708, 1.312474, 1.247701, 0.781237, 0.940834, // S-Ca
        0.950000, 0.974455, 0.998911, 1.023366, 1.047822, // Sc-Mn
        1.072277, 1.096733, 1.121188, 1.145644, 1.170099, // Fe-Zn
        1.205357, 1.145447, 1.169499, 1.253293, 1.329909, // Ga-Br
        1.116527, 0.950975, 0.964592, 0.897786, 0.932824, // Kr-Zr
        0.967863, 1.002901, 1.037940, 1.072978, 1.108017, // Nb-Rh
        1.143055, 1.178094, 1.213132, 1.205076, 1.075529, // Pd-Sn
        1.206919, 1.303658, 1.332656, 1.179317, 0.789115, // Sb-Cs
        0.798704, 0.993208, 0.907847, 0.836114, 0.778008, // Ba-Nd
        0.733529, 0.702678, 0.685455, 0.681858, 0.691889, // Pm-Tb
        0.715548, 0.752834, 0.803747, 0.868288, 0.946457, // Dy-Yb
        1.038252, 1.128780, 1.129764, 1.130747, 1.131731, // Lu-Re
        1.132714, 1.133698, 1.134681, 1.135665, 1.136648, // Os-Hg
        1.061832, 1.053084, 1.207830, 1.236314, 1.310129, // Tl-At
        1.157380, 0.789115, 0.798704, 1.053384, 1.056040, // Rn-Ra
        1.058772, 1.061580, 1.064463, 1.067422, 1.070456, // Ac-Am
        1.073566, 1.076751, 1.080012, 1.083349, 1.086761, // Cm-Cf
        1.090249, 1.093812, 1.097451 // Es-Lr
    };

    static const std::vector<double> gam_angewChem2020 = {
        -0.448428, 0.131022, 0.571431, 0.334622, -0.089208, // H-B
        -0.025895, -0.027280, -0.031236, -0.159892, 0.074198, // C-Ne
        0.316829, 0.326072, 0.069748, -0.120184, -0.193159, // Na-P
        -0.182428, -0.064093, 0.061914, 0.318112, 0.189248, // S-Ca
        -0.104172, -0.082038, -0.059903, -0.037769, -0.015635, // Sc-Mn
        0.006500, 0.028634, 0.050768, 0.072903, 0.095037, // Fe-Zn
        0.131140, 0.097006, -0.065744, -0.058394, 0.063307, // Ga-Br
        0.091652, 0.386337, 0.530677, -0.030705, -0.020787, // Kr-Zr
        -0.010869, -0.000951, 0.008967, 0.018884, 0.028802, // Nb-Rh
        0.038720, 0.048638, 0.058556, 0.036488, 0.077711, // Pd-Sn
        0.077025, 0.004547, 0.039909, 0.082630, 0.485375, // Sb-Cs
        0.498677, 0.192222, 0.221806, 0.229117, 0.236428, // Ba-Nd
        0.243740, 0.251051, 0.258362, 0.265673, 0.272984, // Pm-Tb
        0.280296, 0.287607, 0.294918, 0.302229, 0.309540, // Dy-Yb
        0.316851, 0.324163, 0.068830, 0.064240, 0.059650, // Lu-Re
        0.055060, 0.050471, 0.045881, 0.041291, 0.036701, // Os-Hg
        0.032111, 0.027521, 0.010600, 0.004800, 0.053600, // Tl-At
        0.072000, 0.485375, 0.498677, 0.192222, 0.221806, // Rn-Ra
        0.229117, 0.236428, 0.243740, 0.251051, 0.258362, // Ac-Am
        0.265673, 0.272984, 0.280296, 0.287607, 0.294918, // Cm-Cf
        0.302229, 0.309540, 0.316851 // Es-Lr
    };

    // Phase 4.3: Complete alp (polarizability) array from Fortran gfnff_param.f90
    static const std::vector<double> alp_angewChem2020 = {
        0.585069, 0.432382, 0.628636, 0.743646, 1.167323, // H-He, Li-B
        0.903430, 1.278388, 0.905347, 1.067014, 2.941513, // C-Ne
        0.687680, 0.792170, 1.337040, 1.251409, 1.068295, // Na-P
        1.186476, 1.593532, 2.056749, 0.674196, 0.868052, // S-Ca
        0.575052, 0.613424, 0.651796, 0.690169, 0.728541, // Sc-Mn
        0.766913, 0.805285, 0.843658, 0.882030, 0.920402, // Fe-Zn
        0.877178, 1.422350, 1.405901, 1.646860, 2.001970, // Ga-Br
        2.301695, 1.020617, 0.634141, 0.652752, 0.668845, // Kr-Zr
        0.684938, 0.701032, 0.717125, 0.733218, 0.749311, // Nb-Rh
        0.765405, 0.781498, 0.797591, 1.296844, 1.534068, // Pd-Sn
        1.727781, 1.926871, 2.175548, 2.177702, 0.977079, // Sb-Cs
        0.770260, 0.757372, 0.757352, 0.757332, 0.757313, // Ba-Nd
        0.757293, 0.757273, 0.757253, 0.757233, 0.757213, // Pm-Tb
        0.757194, 0.757174, 0.757154, 0.757134, 0.757114, // Dy-Yb
        0.757095, 0.757075, 0.756778, 0.756480, 0.756183, // Lu-Re
        0.755886, 0.755589, 0.755291, 0.754994, 0.754697, // Os-Hg
        0.868029, 1.684375, 2.001040, 2.067331, 2.228923, // Tl-At
        1.874218, 0.977079, 0.770260, 0.757372, 0.757352, // Rn-Th (approx)
        0.757332, 0.757313, 0.757293, 0.757273, 0.757253, // Pa-Am
        0.757233, 0.757213, 0.757194, 0.757174, 0.757154, // Cm-Cf
        0.757134, 0.757114, 0.757095 // Es-Lr
    };

    // Phase 4.3: Complete cnf_angewChem2020 array - CN correction factor (86 elements)
    static const std::vector<double> cnf_angewChem2020 = {
        0.008904, 0.004641, 0.048324, 0.080316, -0.051990, // H-B
        0.031779, 0.132184, 0.157353, 0.064120, 0.036540, // C-Ne
        -0.000627, 0.005412, 0.018809, 0.016329, 0.012149, // Na-P
        0.021484, 0.014212, 0.014939, 0.003597, 0.032921, // S-Ca
        -0.021804, -0.022797, -0.023789, -0.024782, -0.025775, // Sc-Mn
        -0.026767, -0.027760, -0.028753, -0.029745, -0.030738, // Fe-Zn
        -0.004189, -0.011113, -0.021305, -0.012311, 0.049781, // Ga-Br
        -0.040533, 0.012872, 0.021056, -0.003395, 0.000799, // Kr-Zr
        0.004992, 0.009186, 0.013379, 0.017573, 0.021766, // Nb-Rh
        0.025960, 0.030153, 0.034347, -0.000052, -0.039776, // Pd-Sn
        0.006661, 0.050424, 0.068985, 0.023470, -0.024950, // Sb-Cs
        -0.033006, 0.058973, 0.058595, 0.058217, 0.057838, // Ba-Nd
        0.057460, 0.057082, 0.056704, 0.056326, 0.055948, // Pm-Tb
        0.055569, 0.055191, 0.054813, 0.054435, 0.054057, // Dy-Yb
        0.053679, 0.053300, 0.047628, 0.041955, 0.036282, // Lu-Re
        0.030610, 0.024937, 0.019264, 0.013592, 0.007919, // Os-Hg
        0.006383, -0.089155, -0.001293, 0.019269, 0.074803, // Tl-At
        0.016657 // Rn
    };

    // Get base parameters
    if (z >= 1 && z <= static_cast<int>(chi_angewChem2020.size())) {
        params.chi = chi_angewChem2020[z - 1];
        params.gam = gam_angewChem2020[z - 1];
        params.alp = alp_angewChem2020[z - 1];  // Phase 4.3: Real polarizability (not fixed 5.0)
        params.cnf = cnf_angewChem2020[z - 1];  // Phase 4.3: CN correction factor
    } else {
        // Fallback values
        params.chi = 1.0;
        params.gam = 0.5;
        params.alp = 5.0;  // Fallback for undefined elements
        params.cnf = 0.0;
    }

    // TODO: Add environment-dependent corrections (dxi, dgam)
    // This should include:
    // - Coordination number corrections
    // - Hybridization corrections
    // - Ring corrections
    // - Charge corrections

    params.xi_corr = 0.0; // TODO: Calculate environment correction

    return params;
}

// ============================================================================
// Phase 4.2: GFN-FF Pairwise Non-Bonded Parameter Generation (Claude Generated 2025)
// ============================================================================

json GFNFF::generateGFNFFCoulombPairs() const
{
    /**
     * @brief Generate EEQ-based Coulomb electrostatics pairwise parameters
     *
     * Reference: Phase 3 EEQ charge calculation
     * Formula: E_coul = q_i * q_j * erf(γ_ij * r_ij) / r_ij
     *
     * Claude Generated (2025): Phase 4.2 parameter generation
     */

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("=== generateGFNFFCoulombPairs() START ===");
        CurcumaLogger::param("m_atomcount", std::to_string(m_atomcount));
    }

    json coulomb_pairs = json::array();

    // Calculate EEQ charges (Phase 3 implementation)
    Vector cn = calculateCoordinationNumbers();
    std::vector<int> hyb = determineHybridization();
    std::vector<int> rings = findSmallestRings();
    Vector charges = calculateEEQCharges(cn, hyb, rings);

    // Generate all pairwise Coulomb interactions (i<j to avoid double-counting)
    for (int i = 0; i < m_atomcount; ++i) {
        for (int j = i + 1; j < m_atomcount; ++j) {
            json coulomb;
            coulomb["i"] = i;
            coulomb["j"] = j;
            coulomb["q_i"] = charges[i];
            coulomb["q_j"] = charges[j];

            // Calculate damping parameter: γ_ij = 1 / sqrt(α_i + α_j)
            EEQParameters params_i = getEEQParameters(m_atoms[i]);
            EEQParameters params_j = getEEQParameters(m_atoms[j]);
            double gamma_ij = 1.0 / std::sqrt(params_i.alp + params_j.alp);
            coulomb["gamma_ij"] = gamma_ij;

            // Cutoff radius (50 Bohr ~ 26 Å, typical for electrostatics)
            coulomb["r_cut"] = 50.0;

            coulomb_pairs.push_back(coulomb);

            if (CurcumaLogger::get_verbosity() >= 3) {
                CurcumaLogger::param(fmt::format("coulomb_{}-{}", i, j),
                    fmt::format("q_i={:.6f}, q_j={:.6f}, γ={:.6f}", charges[i], charges[j], gamma_ij));
            }
        }
    }

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::success(fmt::format("Generated {} Coulomb pairs", coulomb_pairs.size()));
    }

    return coulomb_pairs;
}

json GFNFF::generateGFNFFRepulsionPairs() const
{
    /**
     * @brief Generate GFN-FF repulsion pairwise parameters
     *
     * Reference: Fortran gfnff_engrad.F90:407-439 (bonded repulsion)
     * Formula: E_rep = repab * exp(-α*r^1.5) / r
     * Parameters: alpha = sqrt(repa_i * repa_j), repab = repz_i * repz_j * scale
     *
     * Claude Generated (2025): Phase 4.2 parameter generation
     */

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("=== generateGFNFFRepulsionPairs() START ===");
        CurcumaLogger::param("m_atomcount", std::to_string(m_atomcount));
    }

    json repulsion_pairs = json::array();

    // Phase 4.3: Complete GFN-FF repulsion parameters from gfnff_param.f90
    // repa_angewChem2020 - Repulsion exponent parameter (86 elements)
    static const std::vector<double> repa_angewChem2020 = {
        2.639785, 3.575012, 0.732142, 1.159621, 1.561585, // H-B
        1.762895, 2.173015, 2.262269, 2.511112, 3.577220, // C-Ne
        0.338845, 0.693023, 0.678792, 0.804784, 1.012178, // Na-P
        1.103469, 1.209798, 1.167791, 0.326946, 0.595242, // S-Ca
        1.447860, 1.414501, 1.381142, 1.347783, 1.314424, // Sc-Mn
        1.281065, 1.247706, 1.214347, 1.180988, 1.147629, // Fe-Zn
        0.700620, 0.721266, 0.741789, 0.857434, 0.875583, // Ga-Br
        0.835876, 0.290625, 0.554446, 0.623980, 0.696005, // Kr-Zr
        0.768030, 0.840055, 0.912081, 0.984106, 1.056131, // Nb-Rh
        1.128156, 1.200181, 1.272206, 0.478807, 0.479759, // Pd-Sn
        0.579840, 0.595241, 0.644458, 0.655289, 0.574626, // Sb-Cs
        0.560506, 0.682723, 0.684824, 0.686925, 0.689026, // Ba-Nd
        0.691127, 0.693228, 0.695329, 0.697430, 0.699531, // Pm-Tb
        0.701631, 0.703732, 0.705833, 0.707934, 0.710035, // Dy-Yb
        0.712136, 0.714237, 0.745751, 0.777265, 0.808779, // Lu-Re
        0.840294, 0.871808, 0.903322, 0.934836, 0.966350, // Os-Hg
        0.467729, 0.486102, 0.559176, 0.557520, 0.563373, // Tl-At
        0.484713  // Rn (last element)
    };

    // repz - Effective nuclear charges for repulsion (86 elements)
    // NOTE: Fortran uses (/ /) notation, values are simply 1,2,3... for H,He,Li... with some variations
    static const std::vector<double> repz = {
        1., 2.,                                                   // H-He
        1., 2., 3., 4., 5., 6., 7., 8.,                          // Li-Ne
        1., 2., 3., 4., 5., 6., 7., 8.,                          // Na-Ar
        1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12.,       // K-Zn
        3., 4., 5., 6., 7., 8.,                                  // Ga-Kr
        1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12.,       // Rb-Cd
        3., 4., 5., 6., 7., 8.,                                  // In-Xe
        1., 2., 3., 3., 3., 3., 3., 3., 3., 3., 3., 3., 3., 3., 3., 3., 3.,  // Cs-Lu (lanthanides use Z=3)
        4., 5., 6., 7., 8., 9., 10., 11., 12.,                   // Hf-Hg
        3., 4., 5., 6., 7., 8.                                   // Tl-Rn
    };

    // Scaling factors from Fortran gfnff_param.f90:373-374
    // Claude Generated (Nov 2025): Apply bonded/non-bonded repulsion scaling
    static const double REPSCALB = 1.7583;  // Bonded repulsion scaling
    static const double REPSCALN = 0.4270;  // Non-bonded repulsion scaling

    // Reuse existing bond topology (same pattern as generateGFNFFAngles)
    // Claude Generated (Nov 2025): Fix duplicate bond detection - use existing topology
    std::set<std::pair<int, int>> bonded_pairs;
    json bonds = generateGFNFFBonds();  // Get already-detected bonds (consistent with rest of GFN-FF)

    for (const auto& bond : bonds) {
        int i = bond["i"];
        int j = bond["j"];
        bonded_pairs.insert({i, j});
    }

    // Generate all pairwise repulsion interactions
    for (int i = 0; i < m_atomcount; ++i) {
        for (int j = i + 1; j < m_atomcount; ++j) {
            json repulsion;
            repulsion["i"] = i;
            repulsion["j"] = j;

            // Get atomic parameters (with bounds checking)
            int zi = m_atoms[i] - 1; // 0-indexed
            int zj = m_atoms[j] - 1;

            // Phase 4.3: Use complete parameter arrays (86 elements)
            double repa_i = (zi >= 0 && zi < static_cast<int>(repa_angewChem2020.size())) ? repa_angewChem2020[zi] : 2.0;
            double repa_j = (zj >= 0 && zj < static_cast<int>(repa_angewChem2020.size())) ? repa_angewChem2020[zj] : 2.0;
            double repz_i = (zi >= 0 && zi < static_cast<int>(repz.size())) ? repz[zi] : 1.0;
            double repz_j = (zj >= 0 && zj < static_cast<int>(repz.size())) ? repz[zj] : 1.0;

            // Calculate pairwise parameters (match Fortran reference gfnff_engrad.F90:407-439)
            // Claude Generated (Nov 2025): Apply bonded/non-bonded repulsion scaling
            bool is_bonded = (bonded_pairs.count({i, j}) > 0);
            double repulsion_scale = is_bonded ? REPSCALB : REPSCALN;

            repulsion["alpha"] = std::sqrt(repa_i * repa_j);
            repulsion["repab"] = repz_i * repz_j * repulsion_scale;  // Apply scaling based on bonded status
            repulsion["r_cut"] = 50.0; // Cutoff radius (Bohr)

            if (CurcumaLogger::get_verbosity() >= 3) {
                CurcumaLogger::param(fmt::format("repulsion_{}-{}", i, j),
                    fmt::format("bonded=%s, repab=%.6f, alpha=%.6f",
                    is_bonded ? "true" : "false", repulsion["repab"].get<double>(), repulsion["alpha"].get<double>()));
            }

            repulsion_pairs.push_back(repulsion);
        }
    }

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::success(fmt::format("Generated {} repulsion pairs", repulsion_pairs.size()));
    }

    return repulsion_pairs;
}

json GFNFF::generateGFNFFDispersionPairs() const
{
    /**
     * @brief Generate D3/D4 dispersion pairwise parameters with BJ damping
     *
     * Reference: Grimme et al., J. Chem. Phys. 132, 154104 (2010) [D3-BJ]
     * Formula: E_disp = -Σ_ij f_damp(r) * (s6*C6/r^6 + s8*C8/r^8)
     *
     * NOTE: This is a simplified implementation using D3 parameters.
     * Full GFN-FF uses D4 with geometry-dependent C6 coefficients.
     *
     * Claude Generated (2025): Phase 4.2 parameter generation
     */

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("=== generateGFNFFDispersionPairs() START ===");
        CurcumaLogger::param("m_atomcount", std::to_string(m_atomcount));
    }

    json dispersion_pairs = json::array();

    // Phase 4.3: Complete free-atom C6 coefficients (Hartree*Bohr^6)
    // Approximate values based on D3/D4 references (Grimme et al.)
    // NOTE: GFN-FF uses geometry-dependent D4 C6 (future: implement full D4 calculation)
    // For now: simplified free-atom values for initial implementation
    static const std::vector<double> C6_atomic = {
        // H-He
        6.50, 1.42,
        // Li-Ne
        1387.0, 214.0, 99.5, 46.6, 24.2, 15.6, 9.52, 6.38,
        // Na-Ar
        1556.0, 627.0, 528.0, 305.0, 185.0, 134.0, 94.6, 64.3,
        // K-Ca
        3897.0, 2221.0,
        // Sc-Zn
        1383.0, 1044.0, 832.0, 602.0, 552.0, 482.0, 408.0, 373.0, 253.0, 284.0,
        // Ga-Kr
        498.0, 354.0, 246.0, 210.0, 162.0, 129.5,
        // Rb-Sr
        4691.0, 3170.0,
        // Y-Cd
        1968.0, 1677.0, 1263.0, 1028.0, 1390.0, 1029.0, 1118.0, 1251.0, 1225.0, 1225.0,
        // In-Xe
        2896.0, 2290.0, 1896.0, 1830.0, 1612.0, 1416.0,
        // Cs-Ba
        6582.0, 5727.0,
        // La-Lu (lanthanides - approximate values)
        3884.0, 3708.0, 3551.0, 3410.0, 3280.0, 3163.0, 3056.0, 2958.0, 2868.0, 2785.0,
        2708.0, 2638.0, 2573.0, 2513.0, 2458.0,
        // Hf-Hg
        2051.0, 1877.0, 1659.0, 1529.0, 1414.0, 1305.0, 1206.0, 1118.0, 1037.0, 1185.0,
        // Tl-Rn
        3292.0, 3135.0, 2762.0, 2600.0, 2452.0, 2318.0
    };

    // GFN-FF specific parameters (from gfnff_param.f90)
    const double s6 = 1.0;  // C6 scaling factor
    const double s8 = 2.4;  // C8 scaling factor (typical GFN-FF)
    const double a1 = 0.48; // BJ damping parameter 1
    const double a2 = 4.80; // BJ damping parameter 2

    // Generate all pairwise dispersion interactions
    for (int i = 0; i < m_atomcount; ++i) {
        for (int j = i + 1; j < m_atomcount; ++j) {
            json dispersion;
            dispersion["i"] = i;
            dispersion["j"] = j;

            // Get atomic C6 coefficients
            int zi = m_atoms[i] - 1; // 0-indexed
            int zj = m_atoms[j] - 1;

            double C6_i = (zi >= 0 && zi < C6_atomic.size()) ? C6_atomic[zi] : 50.0;
            double C6_j = (zj >= 0 && zj < C6_atomic.size()) ? C6_atomic[zj] : 50.0;

            // Combine rule for C6_ij: geometric mean
            double C6_ij = std::sqrt(C6_i * C6_j);

            // Estimate C8 from C6 (typical ratio C8/C6 ≈ 50 Bohr^2)
            double C8_ij = C6_ij * 50.0;

            dispersion["C6"] = C6_ij;
            dispersion["C8"] = C8_ij;
            dispersion["s6"] = s6;
            dispersion["s8"] = s8;
            dispersion["a1"] = a1;
            dispersion["a2"] = a2;
            dispersion["r_cut"] = 50.0; // Cutoff radius (Bohr)

            dispersion_pairs.push_back(dispersion);
        }
    }

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::success(fmt::format("Generated {} dispersion pairs", dispersion_pairs.size()));
    }

    return dispersion_pairs;
}

// Claude Generated: Energy component getter methods for regression testing (Nov 2025)
double GFNFF::BondEnergy() const {
    if (!m_forcefield) return 0.0;
    return m_forcefield->BondEnergy();
}

double GFNFF::AngleEnergy() const {
    if (!m_forcefield) return 0.0;
    return m_forcefield->AngleEnergy();
}

double GFNFF::DihedralEnergy() const {
    if (!m_forcefield) return 0.0;
    return m_forcefield->DihedralEnergy();
}

double GFNFF::InversionEnergy() const {
    if (!m_forcefield) return 0.0;
    return m_forcefield->InversionEnergy();
}

double GFNFF::VdWEnergy() const {
    if (!m_forcefield) return 0.0;
    return m_forcefield->VdWEnergy();
}

double GFNFF::RepulsionEnergy() const {
    if (!m_forcefield) return 0.0;
    return m_forcefield->RepulsionEnergy();
}

double GFNFF::DispersionEnergy() const {
    if (!m_forcefield) return 0.0;
    return m_forcefield->DispersionEnergy();
}

double GFNFF::CoulombEnergy() const {
    if (!m_forcefield) return 0.0;
    return m_forcefield->CoulombEnergy();
}