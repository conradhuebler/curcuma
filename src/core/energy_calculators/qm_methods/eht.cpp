/*
 * <Extended Hückel Theory Implementation in Curcuma>
 * Copyright (C) 2023 Conrad Hübler <Conrad.Huebler@gmx.net>
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
 *
 * This file implements Extended Hückel Theory calculations using
 * Slater-type orbitals and the Wolfsberg-Helmholz approximation.
 */

#include "eht.h"

#include "interface/abstract_interface.h"
#include "src/core/curcuma_logger.h"
#include "src/core/global.h"
#include "src/core/molecule.h"
#include <fmt/format.h>

#include "ParallelEigenSolver.hpp"
#include "external/CxxThreadPool/include/CxxThreadPool.hpp"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <set>
#include <stdexcept>
#include <vector>

#include <Eigen/Dense>

#include "GTOIntegrals.hpp"
#include "STOIntegrals.hpp"
#include "eht_parameters.h"

#include "json.hpp"

EHT::EHT()
    : m_K(EHTParams::ParameterDatabase::getWolfsbergHelmholzConstant())
{
}

EHT::EHT(double K)
    : m_K(K)
{
}

Basisset EHT::MakeBasis()
{
    m_num_electrons = 0;
    std::vector<STO::Orbital> orbitals;

    // Iterate through all atoms in the molecule
    for (int i = 0; i < m_atoms.size(); ++i) {
        int atomic_number = m_atoms[i];

        // Check if parameters are available for this element
        if (!EHTParams::ParameterDatabase::hasElement(atomic_number)) {
            throw std::runtime_error("EHT parameters not available for element Z=" + std::to_string(atomic_number) + " (" + EHTParams::ParameterDatabase::getElementSymbol(atomic_number) + ")");
        }

        // Get all orbital parameters for this element
        EHTParams::ElementParameters element_params = EHTParams::ParameterDatabase::getElementParameters(atomic_number);

        // Create orbitals for each orbital type of this element
        for (const auto& orbital_pair : element_params) {
            STO::OrbitalType orbital_type = orbital_pair.first;
            EHTParams::OrbitalParameters params = orbital_pair.second;

            STO::Orbital orbital;
            orbital.type = orbital_type;
            orbital.x = m_geometry(i, 0); // x coordinate in Angstrom
            orbital.y = m_geometry(i, 1); // y coordinate in Angstrom
            orbital.z = m_geometry(i, 2); // z coordinate in Angstrom
            orbital.zeta = params.zeta; // Slater exponent
            orbital.VSIP = params.VSIP; // Valence State Ionization Potential
            orbital.atom = i; // Atom index

            orbitals.push_back(orbital);

            // Add electrons contributed by this orbital
            m_num_electrons += params.n_electrons;
        }
    }

    // Verbosity Level 2+: Basis set information
    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("EHT basis set constructed");
        CurcumaLogger::param("total_orbitals", static_cast<int>(orbitals.size()));
        CurcumaLogger::param("total_electrons", m_num_electrons);
    }

    // Detailed basis set summary - Level 3+
    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("Detailed basis set summary:");
        for (size_t i = 0; i < orbitals.size(); ++i) {
            const auto& orb = orbitals[i];
            std::string orbital_name;
            switch (orb.type) {
            case STO::S:
                orbital_name = "s";
                break;
            case STO::PX:
                orbital_name = "px";
                break;
            case STO::PY:
                orbital_name = "py";
                break;
            case STO::PZ:
                orbital_name = "pz";
                break;
            default:
                orbital_name = "unknown";
                break;
            }

            std::string orbital_info = fmt::format("{}:{}{}-{} (ζ={:.3f}, VSIP={:.2f} eV)",
                i + 1,
                EHTParams::ParameterDatabase::getElementSymbol(m_atoms[orb.atom]),
                orb.atom + 1,
                orbital_name,
                orb.zeta,
                orb.VSIP);
            CurcumaLogger::param(fmt::format("Orbital_{}", i + 1), orbital_info);
        }
    }

    return orbitals;
}

bool EHT::InitialiseMolecule()
{
    if (m_atoms.size() == 0) {
        CurcumaLogger::error("No atoms in molecule");
        return false;
    }

    // Validate that all atoms are supported
    if (!validateMolecule()) {
        return false;
    }

    m_num_electrons = 0;
    // Matrix dimensions will be set after basis set construction
    return true;
}

double EHT::Calculation(bool gradient)
{
    m_threads = 1; // EHT typically uses single thread for matrix operations
    // std::cout << m_geometry << std::endl << std::endl;
    //  Check for gradient request (not supported in EHT)
    if (gradient) {
        // std::cout << "Warning: EHT does not support analytical gradients." << std::endl;
    }

    // Validate molecule
    if (m_atoms.size() == 0) {
        CurcumaLogger::error("No molecule set for EHT calculation");
        return 0.0;
    }

    // Level 1+: Starting calculation
    if (CurcumaLogger::get_verbosity() >= 1) {
        CurcumaLogger::info("Starting Extended Hückel Theory calculation");
    }

    // Level 2+: Molecular parameters
    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::param("atoms", static_cast<int>(m_atoms.size()));
        CurcumaLogger::param("molecular_charge", m_charge);
        CurcumaLogger::param("spin_multiplicity", m_spin + 1);
    }

    try {
        // Step 1: Construct basis set
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info("Step 1: Constructing EHT basis set");
        }
        auto basisset = MakeBasis();

        // Initialize matrix dimensions
        int basis_size = static_cast<int>(basisset.size());
        m_mo = Matrix::Zero(basis_size, basis_size);
        m_energies = Vector::Zero(basis_size);

        // Step 2: Construct overlap matrix
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info("Step 2: Computing overlap matrix");
        }
        Matrix S = MakeOverlap(basisset);

        // Step 3: Construct Hamiltonian matrix
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info("Step 3: Constructing Hamiltonian matrix");
        }
        Matrix H = MakeH(S, basisset);

        // Step 4: Solve generalized eigenvalue problem
        if (CurcumaLogger::get_verbosity() >= 3) {
            CurcumaLogger::info("Step 4: Solving generalized eigenvalue problem HC = SCE");
            CurcumaLogger::param("matrix_size", fmt::format("{}x{}", basis_size, basis_size));
        }

        // Use parallel eigenvalue solver for efficiency
        ParallelEigenSolver solver(500, 128, 1e-10, false);
        solver.setThreadCount(m_threads);
        Eigen::MatrixXd mo_coefficients;

        bool success = solver.solve(S, H, m_energies, mo_coefficients, m_threads, false);

        if (!success) {
            CurcumaLogger::error("Eigenvalue solution failed!");
            return 0.0;
        }

        m_mo = mo_coefficients;

        // Step 5: Calculate electronic energy
        double electronic_energy = calculateElectronicEnergy() / eV2Eh;

        // Level 1+: Final energy result
        if (CurcumaLogger::get_verbosity() >= 1) {
            CurcumaLogger::energy_abs(electronic_energy, "EHT Electronic Energy");
        }

        // Level 2+: Orbital analysis
        if (CurcumaLogger::get_verbosity() >= 2) {
            printOrbitalAnalysisVerbose();
        }

        // Level 3+: Detailed orbital analysis (old style)
        if (CurcumaLogger::get_verbosity() >= 3) {
            printOrbitalAnalysis();
        }

        return electronic_energy;

    } catch (const std::exception& e) {
        CurcumaLogger::error(fmt::format("Error during EHT calculation: {}", e.what()));
        return 0.0;
    }
}

Matrix EHT::MakeOverlap(Basisset& basisset)
{
    Matrix S = Eigen::MatrixXd::Zero(basisset.size(), basisset.size());

    // Aktualisiere die Atompositionen in allen Orbitalen
    for (int i = 0; i < basisset.size(); ++i) {
        basisset[i].x = m_geometry(basisset[i].atom, 0) / 0.529177;
        basisset[i].y = m_geometry(basisset[i].atom, 1) / 0.529177;
        basisset[i].z = m_geometry(basisset[i].atom, 2) / 0.529177;
    }

    // Berechne die Überlappungsmatrix direkt
    for (int i = 0; i < basisset.size(); ++i) {
        for (int j = 0; j <= i; ++j) {
            // Die verbesserte calculateOverlap Funktion gibt bereits
            // normalisierte Werte zurück
            double overlap = STO::calculateOverlap(basisset[i], basisset[j]);
            S(i, j) = S(j, i) = overlap;

            // Debug-Ausgabe für kritische Werte
            // if (m_verbose) {
            //    std::cout << "Overlap between orbital " << i << " and " << j
            //              << ": " << overlap << std::endl;
            //}
        }
    }

    return S;
}

Matrix EHT::MakeH(const Matrix& S, const Basisset& basisset)
{
    Matrix H = Eigen::MatrixXd::Zero(basisset.size(), basisset.size());

    // Apply Wolfsberg-Helmholz approximation
    for (size_t i = 0; i < basisset.size(); ++i) {
        const STO::Orbital& orbital_i = basisset[i];

        for (size_t j = 0; j <= i; ++j) { // Only compute lower triangle due to symmetry
            const STO::Orbital& orbital_j = basisset[j];

            if (i == j) {
                // Diagonal elements: H_ii = VSIP_i
                H(i, j) = orbital_i.VSIP;
            } else {
                // Off-diagonal elements: H_ij = K * S_ij * (VSIP_i + VSIP_j) / 2
                double h_ij = m_K * S(i, j) * (orbital_i.VSIP + orbital_j.VSIP) / 2.0;
                H(i, j) = h_ij;
                H(j, i) = h_ij; // Enforce symmetry
            }
        }
    }

    // Level 3+: Hamiltonian construction details
    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("Hamiltonian matrix constructed using Wolfsberg-Helmholz approximation");
        CurcumaLogger::param("K_constant", fmt::format("{:.4f}", m_K));
    }

    return H;
}

bool EHT::validateMolecule() const
{
    for (int i = 0; i < m_atoms.size(); ++i) {
        int atomic_number = m_atoms[i];
        if (!EHTParams::ParameterDatabase::hasElement(atomic_number)) {
            CurcumaLogger::error(fmt::format("EHT parameters not available for element Z={} ({}) at atom {}",
                atomic_number,
                EHTParams::ParameterDatabase::getElementSymbol(atomic_number),
                i + 1));
            return false;
        }
    }
    return true;
}

double EHT::calculateElectronicEnergy() const
{
    if (m_energies.size() == 0 || m_num_electrons == 0) {
        return 0.0;
    }

    double energy = 0.0;
    int occupied_orbitals = m_num_electrons / 2; // Assuming closed shell

    // Sum up energies of occupied orbitals (each orbital holds 2 electrons)
    for (int i = 0; i < occupied_orbitals && i < m_energies.size(); ++i) {
        energy += 2.0 * m_energies(i); // Factor of 2 for electron pairs
    }

    // Handle odd number of electrons
    if (m_num_electrons % 2 == 1 && occupied_orbitals < m_energies.size()) {
        energy += m_energies(occupied_orbitals);
    }

    return energy;
}

double EHT::getHOMOLUMOGap() const
{
    if (m_energies.size() == 0 || m_num_electrons == 0) {
        return 0.0;
    }

    int homo_index = m_num_electrons / 2 - 1;
    int lumo_index = homo_index + 1;

    if (homo_index < 0 || lumo_index >= m_energies.size()) {
        return 0.0;
    }

    return m_energies(lumo_index) - m_energies(homo_index);
}

double EHT::getHOMOEnergy() const
{
    if (m_energies.size() == 0 || m_num_electrons == 0) {
        return 0.0;
    }

    int homo_index = m_num_electrons / 2 - 1;

    if (homo_index < 0 || homo_index >= m_energies.size()) {
        return 0.0;
    }

    return m_energies(homo_index);
}

double EHT::getLUMOEnergy() const
{
    if (m_energies.size() == 0 || m_num_electrons == 0) {
        return 0.0;
    }

    int lumo_index = m_num_electrons / 2;

    if (lumo_index >= m_energies.size()) {
        return 0.0;
    }

    return m_energies(lumo_index);
}

void EHT::printOrbitalAnalysis(int num_orbitals_around_gap) const
{
    if (m_energies.size() == 0 || m_num_electrons == 0) {
        CurcumaLogger::warn("No orbital energies available for analysis");
        return;
    }

    CurcumaLogger::info("=== Extended Hückel Theory Orbital Analysis ===");
    CurcumaLogger::param("total_electrons", m_num_electrons);
    CurcumaLogger::param("total_orbitals", static_cast<int>(m_energies.size()));
    CurcumaLogger::param("K_constant", fmt::format("{:.4f}", m_K));

    // Calculate electronic energy
    double electronic_energy = calculateElectronicEnergy();
    CurcumaLogger::param("electronic_energy", fmt::format("{:.4f} Eh", electronic_energy / eV2Eh));

    // HOMO-LUMO analysis
    int homo_index = m_num_electrons / 2 - 1;
    int lumo_index = homo_index + 1;

    if (homo_index >= 0 && lumo_index < m_energies.size()) {
        double homo_energy = m_energies(homo_index);
        double lumo_energy = m_energies(lumo_index);
        double gap = lumo_energy - homo_energy;

        CurcumaLogger::param("HOMO", fmt::format("{:.4f} eV (orbital {})", homo_energy, homo_index + 1));
        CurcumaLogger::param("LUMO", fmt::format("{:.4f} eV (orbital {})", lumo_energy, lumo_index + 1));
        CurcumaLogger::param("HOMO-LUMO_gap", fmt::format("{:.4f} eV", gap));
    }

    // Print orbitals around HOMO-LUMO gap
    CurcumaLogger::info("Orbitals around HOMO-LUMO gap:");

    int start_orbital = std::max(0, homo_index - num_orbitals_around_gap + 1);
    int end_orbital = std::min(lumo_index + num_orbitals_around_gap - 1,
        static_cast<int>(m_energies.size()) - 1);

    for (int i = start_orbital; i <= end_orbital; ++i) {
        std::string type = "";
        if (i == homo_index)
            type = " HOMO";
        else if (i == lumo_index)
            type = " LUMO";

        int occupation = 0;
        if (i < m_num_electrons / 2)
            occupation = 2;
        else if (i == m_num_electrons / 2 && m_num_electrons % 2 == 1)
            occupation = 1;

        std::string orbital_info = fmt::format("{:.4f} eV (occ={}){}",
            m_energies(i), occupation, type);
        CurcumaLogger::param(fmt::format("Orbital_{}", i + 1), orbital_info);
    }

    // Print lowest and highest energy orbitals
    CurcumaLogger::param("lowest_orbital", fmt::format("{:.4f} eV (orbital 1)", m_energies(0)));
    CurcumaLogger::param("highest_orbital", fmt::format("{:.4f} eV (orbital {})", m_energies(m_energies.size() - 1), m_energies.size()));

    CurcumaLogger::info("==========================================");
}

void EHT::printOrbitalAnalysisVerbose() const
{
    if (m_energies.size() == 0 || m_num_electrons == 0) {
        CurcumaLogger::warn("No orbital energies available for EHT analysis");
        return;
    }

    // Basic molecular orbital information
    CurcumaLogger::param("total_electrons", m_num_electrons);
    CurcumaLogger::param("total_orbitals", static_cast<int>(m_energies.size()));
    CurcumaLogger::param("K_constant", fmt::format("{:.3f}", m_K));

    // Calculate electronic energy
    double electronic_energy = calculateElectronicEnergy();
    CurcumaLogger::param("electronic_energy", fmt::format("{:.4f} Eh", electronic_energy / eV2Eh));

    // HOMO-LUMO analysis
    int homo_index = m_num_electrons / 2 - 1;
    int lumo_index = homo_index + 1;

    if (homo_index >= 0 && lumo_index < m_energies.size()) {
        double homo_energy = m_energies(homo_index);
        double lumo_energy = m_energies(lumo_index);
        double gap = lumo_energy - homo_energy;

        CurcumaLogger::param("HOMO", fmt::format("{:.4f} eV", homo_energy));
        CurcumaLogger::param("LUMO", fmt::format("{:.4f} eV", lumo_energy));
        CurcumaLogger::param("HOMO-LUMO_gap", fmt::format("{:.4f} eV ({:.2f} Eh)", gap, gap / 27.211));
    }

    // Orbital energy extremes
    CurcumaLogger::param("lowest_orbital", fmt::format("{:.4f} eV", m_energies(0)));
    CurcumaLogger::param("highest_orbital", fmt::format("{:.4f} eV", m_energies(m_energies.size() - 1)));
}
