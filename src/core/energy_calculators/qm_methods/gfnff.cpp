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
        CurcumaLogger::param("geometry_rows", std::to_string(m_geometry.rows()));
        CurcumaLogger::param("geometry_cols", std::to_string(m_geometry.cols()));
    }

    if (m_atomcount == 0) {
        CurcumaLogger::error("GFN-FF initialization failed: No atoms in molecule");
        return false;
    }

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

    if (m_forcefield) {
        m_forcefield->UpdateGeometry(m_geometry);
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
            CurcumaLogger::param("geometry_size", std::to_string(m_geometry.rows()) + "x" + std::to_string(m_geometry.cols()));
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

    double energy_kcal = m_forcefield->Calculate(gradient);

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::param("energy_kcal/mol_raw", fmt::format("{:.8f}", energy_kcal));
    }

    if (gradient) {
        Matrix grad_kcal = m_forcefield->Gradient();
        m_gradient = convertGradientToHartree(grad_kcal);

        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::param("gradient_norm", fmt::format("{:.8f}", m_gradient.norm()));
        }
    }

    m_energy_total = convertToHartree(energy_kcal);

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::energy_abs(m_energy_total, "GFN-FF Energy");
        CurcumaLogger::param("energy_kcal/mol", fmt::format("{:.6f}", energy_kcal));
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

    // TEMPORARY DEBUG: Disable caching to isolate the problem
    m_forcefield->setParameterCaching(false);

    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::success("ForceField instance created");
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

        // Load atomic charges from reference file
        loadGFNFFCharges();

        // Generate GFN-FF bonds with real parameters
        json bonds = generateGFNFFBonds();
        json angles = generateGFNFFAngles();
        json torsions = generateGFNFFTorsions(); // ✅ Phase 1.1 implemented
        json inversions = generateGFNFFInversions(); // ✅ Phase 1.2 implemented

        parameters["bonds"] = bonds;
        parameters["angles"] = angles;
        parameters["dihedrals"] = torsions;
        parameters["inversions"] = inversions;

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

    // GFN-FF bond detection with connectivity threshold
    double bond_threshold = 1.3; // Factor for covalent radii sum

    std::cout << "GFN-FF bond detection: " << m_atomcount << " atoms, threshold " << bond_threshold << std::endl;

    for (int i = 0; i < m_atomcount; ++i) {
        for (int j = i + 1; j < m_atomcount; ++j) {
            Vector ri = m_geometry.row(i);
            Vector rj = m_geometry.row(j);
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

                // GFN-FF bond parameters based on element types
                auto bond_params = getGFNFFBondParameters(m_atoms[i], m_atoms[j], distance);
                bond["fc"] = bond_params.force_constant;
                bond["r0_ij"] = bond_params.equilibrium_distance;
                bond["r0_ik"] = 0.0; // Not used in GFN-FF but required by ForceField
                bond["exponent"] = bond_params.alpha;  // Phase 1.3: store α in exponent field

                bonds.push_back(bond);
            }
        }
    }

    if (bonds.empty()) {
        std::cerr << "Warning: No bonds detected in GFN-FF" << std::endl;
    } else {
        std::cout << "GFN-FF detected " << bonds.size() << " bonds" << std::endl;
    }

    return bonds;
}

json GFNFF::generateGFNFFAngles() const
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
                Vector ri = m_geometry.row(neighbors[i]);
                Vector rj = m_geometry.row(center);
                Vector rk = m_geometry.row(neighbors[j]);

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

                // GFN-FF angle parameters
                auto angle_params = getGFNFFAngleParameters(m_atoms[neighbors[i]],
                    m_atoms[center],
                    m_atoms[neighbors[j]],
                    current_angle);

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
    if (m_geometry.rows() != m_atomcount || m_geometry.cols() != 3) {
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
    // Original GFN-FF covalent radii from gfnff_param.f90 (rad array)
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
        return covalent_radii[atomic_number - 1]; // Convert to 0-based indexing
    } else {
        // Fallback for unknown elements
        std::cerr << "Warning: No covalent radius for element " << atomic_number
                  << ", using default 1.0 Å" << std::endl;
        return 1.0;
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

    // gam_angewChem2020(86) - Chemical hardness
    static const std::vector<double> gam_eeq = {
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
        0.416264, -0.011212, -0.011046, -0.010879, -0.010713, // Ba-Nd
        -0.010546, -0.010380, -0.010214, -0.010047, -0.009881, // Pm-Tb
        -0.009714, -0.009548, -0.009382, -0.009215, -0.009049, // Dy-Yb
        -0.008883, -0.008716, -0.006220, -0.003724, -0.001228, // Lu-Re
        0.001267, 0.003763, 0.006259, 0.008755, 0.011251, // Os-Hg
        0.020477, -0.056566, 0.051943, 0.076708, 0.000273, // Tl-At
        -0.068929 // Rn
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

GFNFF::GFNFFBondParams GFNFF::getGFNFFBondParameters(int z1, int z2, double distance) const
{
    GFNFFBondParams params;

    // Original GFN-FF bond parameters from gfnff_param.f90 (bond_angewChem2020 array)
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

    // Get bond parameter for both elements
    double bond_param_1 = (z1 >= 1 && z1 <= static_cast<int>(bond_params.size())) ? bond_params[z1 - 1] : 0.15;
    double bond_param_2 = (z2 >= 1 && z2 <= static_cast<int>(bond_params.size())) ? bond_params[z2 - 1] : 0.15;

    // Use geometric mean of bond parameters as force constant
    params.force_constant = std::sqrt(bond_param_1 * bond_param_2);

    // Pauling electronegativities from gfnff_param.f90 (en array)
    static const std::vector<double> electronegativities = {
        2.200, 3.000, 0.980, 1.570, 2.040, 2.550, 3.040, 3.440, 3.980, // H-F
        4.500, 0.930, 1.310, 1.610, 1.900, 2.190, 2.580, 3.160, 3.500, // Ne-Ar
        0.820, 1.000, 1.360, 1.540, 1.630, 1.660, 1.550, 1.830, 1.880, // K-Co
        1.910, 1.900, 1.650, 1.810, 2.010, 2.180, 2.550, 2.960, 3.000, // Ni-Kr
        0.820, 0.950, 1.220, 1.330, 1.600, 2.160, 1.900, 2.200, 2.280, // Rb-Rh
        2.200, 1.930, 1.690, 1.780, 1.960, 2.050, 2.100, 2.660, 2.600, // Pd-Xe
        0.79, 0.89, 1.10, 1.12, 1.13, 1.14, 1.15, 1.17, 1.18, 1.20, 1.21, 1.22, // Cs-Gd
        1.23, 1.24, 1.25, 1.26, 1.27, 1.3, 1.5, 1.7, 1.9, 2.1, 2.2, 2.2, 2.2, // Tb-Au
        2.00, 1.62, 2.33, 2.02, 2.0, 2.2, 2.2 // Hg-Lr
    };

    // Get electronegativities
    double en1 = (z1 >= 1 && z1 <= static_cast<int>(electronegativities.size())) ? electronegativities[z1 - 1] : 2.0;
    double en2 = (z2 >= 1 && z2 <= static_cast<int>(electronegativities.size())) ? electronegativities[z2 - 1] : 2.0;

    // Force constant: geometric mean of bond parameters
    // Full GFN-FF: fc = bond(i)*bond(j) * ringf * bstrength * fqq * fheavy * fpi * fxh * fcn
    params.force_constant = std::sqrt(bond_param_1 * bond_param_2);

    // Phase 2: Apply basic topology corrections
    // Note: Full topology corrections (ring, pi-system) applied in generateTopologyAwareBonds()

    // Heavy atom correction (elements > period 2)
    double heavy_factor = 1.0;
    if (z1 > 10 || z2 > 10) {
        heavy_factor = 0.95; // Slightly weaker for heavy atoms
    }
    params.force_constant *= heavy_factor;

    // Equilibrium distance: use current distance (topology-dependent)
    params.equilibrium_distance = distance;

    // Alpha parameter (exponential decay): Fortran vbond(2,i) = srb1*(1.0 + fsrb2*ΔEN² + srb3*bstrength)
    double srb1 = 16.0;  // Base exponential decay parameter
    double fsrb2 = 0.1;   // Electronegativity scaling factor
    double en_diff = en1 - en2;

    // Phase 2: Add bond order correction (simplified)
    // Double bonds: α *= 1.2, Triple bonds: α *= 1.4
    // Will be refined with actual bond order detection
    double srb3 = 0.5;  // Bond strength scaling
    double bstrength = 1.0; // Assume single bond (will be improved with topology)

    params.alpha = srb1 * (1.0 + fsrb2 * en_diff * en_diff + srb3 * (bstrength - 1.0));

    return params;
}

GFNFF::GFNFFAngleParams GFNFF::getGFNFFAngleParameters(int z1, int z2, int z3, double current_angle) const
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

    // Get angle parameter for center atom (z2)
    double angle_param = (z2 >= 1 && z2 <= static_cast<int>(angle_params.size())) ? angle_params[z2 - 1] : 0.1;

    // Phase 1.3: Simplified force constant (without full topology corrections)
    // Full GFN-FF: k_ijk = angl(center)*angl2(i)*angl2(k) * fqq * f2 * fn * fbsmall * feta
    // For now: use only center atom parameter
    params.force_constant = angle_param * 0.001;  // Scale to match kcal/mol units

    // Use current angle as equilibrium angle
    // NOTE: Full GFN-FF calculates θ₀ from ideal geometries (sp/sp²/sp³)
    params.equilibrium_angle = current_angle;

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
    Vector cn = Vector::Zero(m_atomcount);

    // TODO: Implement coordination number calculation based on gfnff_dlogcoord
    // For now, use simple distance-based approach
    for (int i = 0; i < m_atomcount; ++i) {
        double cn_i = 0.0;
        for (int j = 0; j < m_atomcount; ++j) {
            if (i == j)
                continue;

            Vector ri = m_geometry.row(i);
            Vector rj = m_geometry.row(j);
            double distance = (ri - rj).norm();

            double rcov_i = getCovalentRadius(m_atoms[i]);
            double rcov_j = getCovalentRadius(m_atoms[j]);

            // Coordination number contribution with exponential decay
            double r_cov = rcov_i + rcov_j;
            double exp_arg = -16.0 * (distance / r_cov - 1.0);
            if (exp_arg > -threshold) {
                cn_i += 1.0 / (1.0 + exp(exp_arg));
            }
        }
        cn[i] = cn_i;
    }

    return cn;
}

std::vector<Matrix> GFNFF::calculateCoordinationNumberDerivatives(const Vector& cn, double threshold) const
{
    // Phase 3.2: CN derivatives for gradient calculations
    // Reference: gfnff_engrad.F90:802-853 (dncoord_erf subroutine)
    //
    // Mathematical formulation:
    // CN_i = Σ_j f(r_ij) where f(r) = 1 / (1 + exp(-16*(r/rcov - 1)))
    //
    // Derivatives:
    // ∂CN_i/∂r_ij = df/dr = 16/rcov * exp(a) / (1 + exp(a))^2
    // where a = -16*(r/rcov - 1)
    //
    // Chain rule for Cartesian coordinates:
    // ∂CN_i/∂x_k = Σ_j ∂CN_i/∂r_ij * ∂r_ij/∂x_k
    //
    // Tensor structure: dcn[dim][i][j] = ∂CN_i/∂coord_dim(atom_j)

    // Initialize 3D tensor as vector of matrices
    // dcn[0] = ∂CN/∂x, dcn[1] = ∂CN/∂y, dcn[2] = ∂CN/∂z
    std::vector<Matrix> dcn(3, Matrix::Zero(m_atomcount, m_atomcount));

    // Calculate derivatives for all atom pairs
    for (int i = 0; i < m_atomcount; ++i) {
        Vector ri = m_geometry.row(i);
        double rcov_i = getCovalentRadius(m_atoms[i]);

        for (int j = 0; j < m_atomcount; ++j) {
            if (i == j) continue;

            Vector rj = m_geometry.row(j);
            double rcov_j = getCovalentRadius(m_atoms[j]);

            // Distance and direction
            Vector r_ij_vec = rj - ri;
            double r_ij = r_ij_vec.norm();
            double rcov_sum = rcov_i + rcov_j;

            // Exponential argument: a = -16*(r/rcov - 1)
            double exp_arg = -16.0 * (r_ij / rcov_sum - 1.0);

            // Only compute if within threshold (same as CN calculation)
            if (exp_arg > -threshold) {
                // Derivative: dCN/dr = 16/rcov * exp(a) / (1 + exp(a))^2
                double exp_val = std::exp(exp_arg);
                double denom = 1.0 + exp_val;
                double dCN_dr = (16.0 / rcov_sum) * exp_val / (denom * denom);

                // Chain rule: dCN/dx_k = dCN/dr * dr/dx_k
                // where dr/dx_k = (x_j - x_i) / r for k=j
                //       dr/dx_k = -(x_j - x_i) / r for k=i
                Vector grad_direction = r_ij_vec / r_ij;  // Unit vector i→j

                // ∂CN_i/∂x_j (moving atom j affects CN_i)
                for (int dim = 0; dim < 3; ++dim) {
                    dcn[dim](i, j) += dCN_dr * grad_direction[dim];
                }

                // ∂CN_i/∂x_i (moving atom i affects CN_i)
                for (int dim = 0; dim < 3; ++dim) {
                    dcn[dim](i, i) -= dCN_dr * grad_direction[dim];
                }
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
        Vector ri = m_geometry.row(i);

        // Step 1: Find bonded neighbors
        std::vector<int> neighbors;
        std::vector<Vector> bond_vectors;

        for (int j = 0; j < m_atomcount; ++j) {
            if (i == j) continue;

            Vector rj = m_geometry.row(j);
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

            double distance = (m_geometry.row(i) - m_geometry.row(j)).norm();
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
            double distance = (m_geometry.row(i) - m_geometry.row(j)).norm();
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
            Vector r_ij_vec = m_geometry.row(i) - m_geometry.row(j);
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
            Vector r_ij_vec = m_geometry.row(i) - m_geometry.row(j);
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

    // Phase 2: Topology-aware bond parameter generation
    // Start with basic GFN-FF bond detection
    double bond_threshold = 1.3; // Factor for covalent radii sum

    for (int i = 0; i < m_atomcount; ++i) {
        for (int j = i + 1; j < m_atomcount; ++j) {
            Vector ri = m_geometry.row(i);
            Vector rj = m_geometry.row(j);
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

                // Get basic bond parameters
                auto bond_params = getGFNFFBondParameters(m_atoms[i], m_atoms[j], distance);

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
                Vector ri = m_geometry.row(neighbors[i]);
                Vector rj = m_geometry.row(center);
                Vector rk = m_geometry.row(neighbors[j]);

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
                auto angle_params = getGFNFFAngleParameters(m_atoms[neighbors[i]],
                    m_atoms[center],
                    m_atoms[neighbors[j]],
                    current_angle);

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
        }
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

            // Calculate pairwise parameters
            repulsion["alpha"] = std::sqrt(repa_i * repa_j);
            repulsion["repab"] = repz_i * repz_j * 0.5; // Scale factor 0.5 (typical)
            repulsion["r_cut"] = 50.0; // Cutoff radius (Bohr)

            repulsion_pairs.push_back(repulsion);
        }
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

    return dispersion_pairs;
}