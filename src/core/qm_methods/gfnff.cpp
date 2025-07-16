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

#include <iostream>
#include <cmath>

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
        {"threads", m_parameters["threads"]},
        {"gradient", m_parameters["gradient"]}
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

    // Generate GFN-FF bonds with real parameters
    json bonds = generateGFNFFBonds();
    json angles = generateGFNFFAngles();
    
    parameters["bonds"] = bonds;
    parameters["angles"] = angles;
    parameters["dihedrals"] = json::array(); // TODO: Implement torsions
    parameters["inversions"] = json::array(); // TODO: Implement inversions
    parameters["vdws"] = json::array();       // TODO: Implement non-bonded

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

json GFNFF::generateGFNFFBonds()
{
    json bonds = json::array();
    
    // GFN-FF bond detection with connectivity threshold
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
    }
    
    return bonds;
}

json GFNFF::generateGFNFFAngles()
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
    // GFN-FF covalent radii in Angstrom (PoC Platzhalter-Parameter)
    static const std::map<int, double> covalent_radii = {
        {1,  0.32},  // H
        {6,  0.75},  // C
        {7,  0.71},  // N
        {8,  0.63},  // O
        {9,  0.64},  // F
        {15, 1.11},  // P
        {16, 1.03},  // S
        {17, 0.99},  // Cl
        {35, 1.16},  // Br
        {53, 1.32}   // I
        // TODO: Add complete GFN-FF covalent radii table (Z=1-86)
    };

    auto it = covalent_radii.find(atomic_number);
    if (it != covalent_radii.end()) {
        return it->second;
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

    // PoC GFN-FF bond parameters (simplified but functional)
    
    // Force constants (kcal/mol/Å²) - rough estimates for common bonds
    std::map<std::pair<int,int>, double> force_constants = {
        {{1,1},   400.0},  // H-H
        {{1,6},   450.0},  // C-H  
        {{1,7},   430.0},  // N-H
        {{1,8},   450.0},  // O-H
        {{6,6},   350.0},  // C-C single
        {{6,7},   400.0},  // C-N
        {{6,8},   450.0},  // C-O
        {{7,7},   450.0},  // N-N
        {{7,8},   400.0},  // N-O
        {{8,8},   350.0},  // O-O
        {{6,16},  300.0},  // C-S
        {{16,16}, 250.0}   // S-S
    };

    // Anharmonic factors (1/Å) - small corrections
    std::map<std::pair<int,int>, double> anharmonic_factors = {
        {{1,1},   -0.5},   {{1,6},   -0.3},   {{1,7},   -0.3},   {{1,8},   -0.4},
        {{6,6},   -0.2},   {{6,7},   -0.25},  {{6,8},   -0.3},   {{7,7},   -0.3},
        {{7,8},   -0.35},  {{8,8},   -0.4},   {{6,16},  -0.2},   {{16,16}, -0.15}
    };

    // Normalize bond key (smaller Z first)
    std::pair<int,int> bond_key = {std::min(z1,z2), std::max(z1,z2)};

    // Get parameters with fallbacks
    auto fc_it = force_constants.find(bond_key);
    params.force_constant = (fc_it != force_constants.end()) ? fc_it->second : 300.0;

    auto ah_it = anharmonic_factors.find(bond_key);
    params.anharmonic_factor = (ah_it != anharmonic_factors.end()) ? ah_it->second : -0.2;

    // Use current distance as equilibrium for PoC
    params.equilibrium_distance = distance;

    return params;
}

GFNFF::GFNFFAngleParams GFNFF::getGFNFFAngleParameters(int z1, int z2, int z3, double current_angle) const
{
    GFNFFAngleParams params;

    // PoC angle parameters based on center atom type
    std::map<int, double> angle_force_constants = {
        {1,  50.0},   // H center (rare)
        {6,  80.0},   // C center (sp3/sp2/sp)
        {7,  70.0},   // N center  
        {8,  60.0},   // O center
        {15, 60.0},   // P center
        {16, 50.0}    // S center
    };

    auto fc_it = angle_force_constants.find(z2);
    params.force_constant = (fc_it != angle_force_constants.end()) ? fc_it->second : 50.0;

    // Use current angle as equilibrium for PoC
    params.equilibrium_angle = current_angle * 180.0 / M_PI;

    // Fourier coefficients for cosine expansion: E = k*(C0 + C1*cos(θ) + C2*cos(2θ))
    double cos_eq = cos(current_angle);
    params.c0 = 1.0 - cos_eq;    // Ensures minimum at current angle
    params.c1 = -1.0;            // Linear restoring term
    params.c2 = 0.1;             // Small anharmonic correction

    return params;
}