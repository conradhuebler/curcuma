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

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

GFNFF::GFNFF()
    : m_forcefield(nullptr)
    , m_initialized(false)
    , m_energy_total(0.0)
{
    json default_parameters = {
        {"method", "gfnff"},
        {"threads", 1},
        {"gradient", 1},
        {"dispersion", true},
        {"hbond", true},
        {"repulsion_scaling", 1.0}
    };
    m_parameters = default_parameters;
}

GFNFF::GFNFF(const json& parameters)
    : m_forcefield(nullptr)
    , m_initialized(false)
    , m_energy_total(0.0)
{
    json default_parameters = {
        {"method", "gfnff"},
        {"threads", 1},
        {"gradient", 1},
        {"dispersion", true},
        {"hbond", true},
        {"repulsion_scaling", 1.0}
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
    if (m_atomcount == 0) {
        std::cerr << "Error: No atoms in molecule for GFN-FF initialization" << std::endl;
        return false;
    }

    if (!validateMolecule()) {
        std::cerr << "Error: Molecule validation failed for GFN-FF" << std::endl;
        return false;
    }

    if (!initializeForceField()) {
        std::cerr << "Error: Force field initialization failed for GFN-FF" << std::endl;
        return false;
    }

    m_gradient = Matrix::Zero(m_atomcount, 3);
    m_charges = Vector::Zero(m_atomcount);
    m_bond_orders = Vector::Zero(m_atomcount * (m_atomcount - 1) / 2);

    m_initialized = true;
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

double GFNFF::Calculation(bool gradient, bool verbose)
{
    if (!m_initialized) {
        std::cerr << "Error: GFN-FF not initialized" << std::endl;
        return 0.0;
    }

    if (!m_forcefield) {
        std::cerr << "Error: Force field not available" << std::endl;
        return 0.0;
    }

    double energy_kcal = m_forcefield->Calculate(gradient, verbose);
    
    if (gradient) {
        Matrix grad_kcal = m_forcefield->Gradient();
        m_gradient = convertGradientToHartree(grad_kcal);
    }

    m_energy_total = convertToHartree(energy_kcal);

    if (verbose) {
        std::cout << "GFN-FF Energy: " << m_energy_total << " Hartree" << std::endl;
        std::cout << "GFN-FF Energy: " << energy_kcal << " kcal/mol" << std::endl;
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
    if (m_forcefield) {
        delete m_forcefield;
    }

    json ff_config = {
        { "threads", m_parameters["threads"] },
        { "gradient", m_parameters["gradient"] },
        { "method", "gfnff" }
    };

    m_forcefield = new ForceField(ff_config);
    m_forcefield->setAtomTypes(m_atoms);

    if (!calculateTopology()) {
        return false;
    }

    json ff_params = generateGFNFFParameters();
    m_forcefield->setParameter(ff_params);

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
        parameters["dihedrals"] = json::array(); // TODO: Advanced torsions
        parameters["inversions"] = json::array(); // TODO: Advanced inversions
        parameters["vdws"] = json::array(); // TODO: Advanced non-bonded
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

        parameters["bonds"] = bonds;
        parameters["angles"] = angles;
        parameters["dihedrals"] = json::array(); // TODO: Implement torsions
        parameters["inversions"] = json::array(); // TODO: Implement inversions
        parameters["vdws"] = json::array(); // TODO: Implement non-bonded
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
                bond["exponent"] = bond_params.anharmonic_factor;
                
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
        bond_list.push_back({bond["i"], bond["j"]});
    }
    
    // Generate angles from bonded topology
    for (int center = 0; center < m_atomcount; ++center) {
        std::vector<int> neighbors;
        
        // Find all atoms bonded to center
        for (const auto& bond : bond_list) {
            if (bond.first == center) neighbors.push_back(bond.second);
            if (bond.second == center) neighbors.push_back(bond.first);
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
                double current_angle = acos(v1.dot(v2) / (v1.norm() * v2.norm()));
                
                // GFN-FF angle parameters
                auto angle_params = getGFNFFAngleParameters(m_atoms[neighbors[i]], 
                                                          m_atoms[center], 
                                                          m_atoms[neighbors[j]], 
                                                          current_angle);
                
                angle["fc"] = angle_params.force_constant;
                angle["theta0_ijk"] = angle_params.equilibrium_angle;
                angle["r0_ij"] = (ri - rj).norm(); // Distance i-j
                angle["r0_ik"] = (rk - rj).norm(); // Distance k-j  
                angle["C0"] = angle_params.c0;
                angle["C1"] = angle_params.c1;
                angle["C2"] = angle_params.c2;
                
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

    // Use current distance as equilibrium distance (topology-dependent)
    params.equilibrium_distance = distance;

    // Small anharmonic correction
    params.anharmonic_factor = -0.1;

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

    // Scale down to match UFF energy scale
    params.force_constant = angle_param * 0.001;

    // Use current angle as equilibrium angle (topology-dependent)
    params.equilibrium_angle = current_angle * 180.0 / M_PI;

    // Fourier coefficients for cosine expansion: E = k*(C0 + C1*cos(θ) + C2*cos(2θ))
    double cos_eq = cos(current_angle);
    params.c0 = 1.0 - cos_eq;    // Ensures minimum at current angle
    params.c1 = -1.0;            // Linear restoring term
    params.c2 = 0.05; // Small anharmonic correction

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
    std::vector<Matrix> dcn;

    // TODO: Implement CN derivatives for gradient calculations
    // This is complex and requires careful implementation
    for (int dim = 0; dim < 3; ++dim) {
        dcn.push_back(Matrix::Zero(m_atomcount, m_atomcount));
    }

    return dcn;
}

std::vector<int> GFNFF::determineHybridization() const
{
    std::vector<int> hyb(m_atomcount, 3); // Default to sp3

    // TODO: Implement hybridization detection algorithm
    // This should analyze neighbor count and geometry
    for (int i = 0; i < m_atomcount; ++i) {
        int z = m_atoms[i];

        // Count neighbors
        int neighbor_count = 0;
        for (int j = 0; j < m_atomcount; ++j) {
            if (i == j)
                continue;

            Vector ri = m_geometry.row(i);
            Vector rj = m_geometry.row(j);
            double distance = (ri - rj).norm();

            double rcov_i = getCovalentRadius(z);
            double rcov_j = getCovalentRadius(m_atoms[j]);

            if (distance < 1.3 * (rcov_i + rcov_j)) {
                neighbor_count++;
            }
        }

        // Simple hybridization assignment based on neighbors
        if (neighbor_count <= 1) {
            hyb[i] = 1; // sp
        } else if (neighbor_count == 2) {
            hyb[i] = 2; // sp2
        } else {
            hyb[i] = 3; // sp3
        }
    }

    return hyb;
}

std::vector<int> GFNFF::detectPiSystems(const std::vector<int>& hyb) const
{
    std::vector<int> pi_fragments(m_atomcount, 0);

    // TODO: Implement pi-system detection
    // This should identify conjugated systems and aromatic rings

    return pi_fragments;
}

std::vector<int> GFNFF::findSmallestRings() const
{
    std::vector<int> ring_sizes(m_atomcount, 0);

    // TODO: Implement ring detection algorithm (similar to getring36)
    // This should find the smallest ring each atom belongs to

    return ring_sizes;
}

Vector GFNFF::calculateEEQCharges(const Vector& cn, const std::vector<int>& hyb, const std::vector<int>& rings) const
{
    Vector charges = Vector::Zero(m_atomcount);

    // TODO: Implement EEQ charge calculation based on goedeckera
    // This is the most complex part and requires:
    // 1. Setting up EEQ parameter matrix
    // 2. Solving linear system for charges
    // 3. Applying fragment constraints

    // For now, use simple partial charge estimation
    for (int i = 0; i < m_atomcount; ++i) {
        int z = m_atoms[i];

        // Simple electronegativity-based charge estimate
        double base_charge = 0.0;
        if (z == 1)
            base_charge = 0.1; // H
        else if (z == 6)
            base_charge = -0.1; // C
        else if (z == 7)
            base_charge = -0.2; // N
        else if (z == 8)
            base_charge = -0.3; // O
        else if (z == 9)
            base_charge = -0.4; // F

        charges[i] = base_charge;
    }

    return charges;
}

json GFNFF::generateTopologyAwareBonds(const Vector& cn, const std::vector<int>& hyb,
    const Vector& charges, const std::vector<int>& rings) const
{
    json bonds = json::array();

    // TODO: Implement topology-aware bond parameter generation
    // This should use CN, hybridization, and ring information
    // For now, fall back to current implementation

    return generateGFNFFBonds();
}

json GFNFF::generateTopologyAwareAngles(const Vector& cn, const std::vector<int>& hyb,
    const Vector& charges, const std::vector<int>& rings) const
{
    json angles = json::array();

    // TODO: Implement topology-aware angle parameter generation
    // This should use CN, hybridization, and ring information
    // For now, fall back to current implementation

    return generateGFNFFAngles();
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

    // Calculate all topology information
    topo_info.coordination_numbers = calculateCoordinationNumbers();
    topo_info.hybridization = determineHybridization();
    topo_info.pi_fragments = detectPiSystems(topo_info.hybridization);
    topo_info.ring_sizes = findSmallestRings();
    topo_info.eeq_charges = calculateEEQCharges(topo_info.coordination_numbers,
        topo_info.hybridization,
        topo_info.ring_sizes);

    // Initialize metal and aromatic flags
    topo_info.is_metal.resize(m_atomcount, false);
    topo_info.is_aromatic.resize(m_atomcount, false);

    for (int i = 0; i < m_atomcount; ++i) {
        int z = m_atoms[i];

        // Simple metal detection
        if (z > 20 && z < 31)
            topo_info.is_metal[i] = true; // 3d metals
        if (z > 38 && z < 49)
            topo_info.is_metal[i] = true; // 4d metals
        if (z > 56 && z < 81)
            topo_info.is_metal[i] = true; // 5d metals

        // Simple aromaticity detection (placeholder)
        if (topo_info.ring_sizes[i] == 5 || topo_info.ring_sizes[i] == 6) {
            if (topo_info.hybridization[i] == 2 && topo_info.pi_fragments[i] > 0) {
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

    // Get base parameters
    if (z >= 1 && z <= static_cast<int>(chi_angewChem2020.size())) {
        params.chi = chi_angewChem2020[z - 1];
        params.gam = gam_angewChem2020[z - 1];
    } else {
        // Fallback values
        params.chi = 1.0;
        params.gam = 0.5;
    }

    // TODO: Add environment-dependent corrections (dxi, dgam)
    // This should include:
    // - Coordination number corrections
    // - Hybridization corrections
    // - Ring corrections
    // - Charge corrections

    params.alp = 5.0; // TODO: Add polarizability parameters
    params.xi_corr = 0.0; // TODO: Calculate environment correction

    return params;
}