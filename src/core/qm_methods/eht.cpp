/*
 * <Extendend Hückel Theory Implementation in Cucuma. >
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
 */

#include "interface/abstract_interface.h"
#include "src/core/global.h"
#include "src/core/molecule.h"

#include "external/CxxThreadPool/include/CxxThreadPool.hpp"

#include "ParallelEigenSolver.hpp"

#include <set>
#include <vector>

#include <Eigen/Dense>

#include "GTOIntegrals.hpp"
#include "STOIntegrals.hpp"
#include "basissetparser.hpp"

#include "json.hpp"

#include "eht.h"

#include <iomanip>

EHT::EHT()
{
}

Basisset EHT::MakeBasis()
{
    m_num_electrons = 0;
    std::vector<STO::Orbital> orbitals;
    for (int i = 0; i < m_mol.m_number_atoms; ++i) {
        if (m_mol.m_atoms[i] == 1) {
            m_num_electrons += 1;
            STO::Orbital o;
            o.type = STO::S;
            o.x = m_mol.m_geometry(i, 0);
            o.y = m_mol.m_geometry(i, 1);
            o.z = m_mol.m_geometry(i, 2);
            o.zeta = 1.2;
            o.VSIP = -0.5;
            o.atom = i;
            orbitals.push_back(o);
        } else if (m_mol.m_atoms[i] == 6) {
            m_num_electrons += 4;
            STO::Orbital o;
            o.type = STO::S;
            o.x = m_mol.m_geometry(i, 0);
            o.y = m_mol.m_geometry(i, 1);
            o.z = m_mol.m_geometry(i, 2);
            o.zeta = 1.625;
            o.VSIP = -0.7144;
            o.atom = i;
            orbitals.push_back(o);

            o.type = STO::PX;
            o.x = m_mol.m_geometry(i, 0);
            o.y = m_mol.m_geometry(i, 1);
            o.z = m_mol.m_geometry(i, 2);
            o.zeta = 1.625;
            o.VSIP = -0.3921;
            orbitals.push_back(o);

            o.type = STO::PY;
            o.x = m_mol.m_geometry(i, 0);
            o.y = m_mol.m_geometry(i, 1);
            o.z = m_mol.m_geometry(i, 2);
            o.zeta = 1.625;
            o.VSIP = -0.3921;
            orbitals.push_back(o);

            o.type = STO::PZ;
            o.x = m_mol.m_geometry(i, 0);
            o.y = m_mol.m_geometry(i, 1);
            o.z = m_mol.m_geometry(i, 2);
            o.zeta = 1.625;
            o.VSIP = -0.3921;
            orbitals.push_back(o);
        } else if (m_mol.m_atoms[i] == 8) {
            m_num_electrons += 6;
            STO::Orbital o;
            o.type = STO::S;
            o.x = m_mol.m_geometry(i, 0);
            o.y = m_mol.m_geometry(i, 1);
            o.z = m_mol.m_geometry(i, 2);
            o.zeta = 2.2399;
            o.VSIP = -1.1904;
            o.atom = i;
            orbitals.push_back(o);

            o.type = STO::PX;
            o.x = m_mol.m_geometry(i, 0);
            o.y = m_mol.m_geometry(i, 1);
            o.z = m_mol.m_geometry(i, 2);
            o.zeta = 2.0477;
            o.VSIP = -1.1904;
            o.atom = i;
            orbitals.push_back(o);

            o.type = STO::PY;
            o.x = m_mol.m_geometry(i, 0);
            o.y = m_mol.m_geometry(i, 1);
            o.z = m_mol.m_geometry(i, 2);
            o.zeta = 2.0477;
            o.VSIP = -1.1904;
            o.atom = i;
            orbitals.push_back(o);

            o.type = STO::PZ;
            o.x = m_mol.m_geometry(i, 0);
            o.y = m_mol.m_geometry(i, 1);
            o.z = m_mol.m_geometry(i, 2);
            o.zeta = 2.0477;
            o.VSIP = -1.1904;
            o.atom = i;
            orbitals.push_back(o);
        }
    }
    return orbitals;
}

bool EHT::InitialiseMolecule()
{
    m_num_electrons = 0;
    m_mo = Matrix::Zero(m_mol.m_number_atoms, m_mol.m_number_atoms);
    m_energies = Vector::Zero(m_mol.m_number_atoms);
    return true;
}

double EHT::Calculation(bool gradient, bool verbose)
{
    m_threads = 1;
    double energy = 0;
    m_verbose = verbose;
    if (gradient) {
        std::cout << "EHT does not support gradients." << std::endl;
    }

    if (m_mol.m_number_atoms == 0) {
        std::cout << "No molecule set." << std::endl;
        return 0;
    }

    // m_molecule.print_geom();
    auto basisset = MakeBasis();
    if (m_verbose) {
        std::cout << basisset.size() << std::endl;
        std::cout << "Making overlap matrix..." << std::endl;
    }
    Matrix S = MakeOverlap(basisset);
    /*
    if (m_verbose)
        std::cout << S << std::endl;
        */
    if (m_verbose)
        std::cout << "Making EHT hamiltonian eg Fock matrix." << std::endl;
    Matrix H = MakeH(S, basisset);
    /*if (m_verbose) {
        std::cout << std::endl
                  << std::endl;
        std::cout << H << std::endl;
    }*/
    std::cout << m_num_electrons << std::endl;
    if (m_verbose)
        std::cout << "Diagonalizing Fock matrix parallel..." << std::endl;

    ParallelEigenSolver solver(500, 128, 1e-10, true);
    solver.setThreadCount(m_threads);
    Eigen::MatrixXd mo;
    // solver.setDebug(true);
    bool success = solver.solve(S, H, m_energies, mo, m_threads, false);
    m_mo = mo;
    /*

    Eigen::SelfAdjointEigenSolver<Matrix> es_S(S);

    // Direkte Berechnung von S^(-1/2) aus den Eigenwerten/Eigenvektoren
    // S = U * D * U^T, also S^(-1/2) = U * D^(-1/2) * U^T
    const Matrix& U = es_S.eigenvectors();
    const Eigen::VectorXd& D = es_S.eigenvalues();

    // Berechnung der Diagonalmatrix D^(-1/2)
    Matrix D_1_2 = Matrix::Zero(S.rows(), S.cols());
    for (int i = 0; i < D.size(); ++i) {
        if (D(i) > 1e-10) { // Numerische Stabilität
            D_1_2(i, i) = 1.0 / std::sqrt(D(i));
        }
    }

    // Berechnung von S^(-1/2)
    Matrix S_1_2 = U * D_1_2 * U.transpose();

    // Berechnung der transformierten Fock-Matrix
    Matrix F = S_1_2.transpose() * H * S_1_2;

    // Diagonalisierung der transformierten Fock-Matrix
    Eigen::SelfAdjointEigenSolver<Matrix> es_F(F, Eigen::ComputeEigenvectors);

    // Extrahieren der Eigenwerte und Eigenvektoren
    m_energies = es_F.eigenvalues();
    m_mo = S_1_2 * es_F.eigenvectors(); // Rücktransformation der MOs
    */
    std::cout << m_energies << std::endl;

    if (m_verbose) {
        std::cout << "\n\nOrbitalenergien und Besetzung:\n";

        // Orbital mit niedrigster Energie (erstes Orbital)
        std::cout << "\nOrbital mit niedrigster Energie:\n";
        std::cout << "Orbital 1: " << m_energies(0) << " eV (Besetzung: 2)\n";

        // HOMO-LUMO Gap und Umgebung (5 Orbitale vor HOMO und 5 nach LUMO)
        int homo_index = m_num_electrons / 2 - 1;
        int lumo_index = homo_index + 1;
        double homo_lumo_gap = 0; // m_mo.eigenvalues()(lumo_index) - m_mo.eigenvalues()(homo_index);

        std::cout << "\nHOMO-LUMO Gap: " << homo_lumo_gap << " eV\n\n";

        // 5 Orbitale vor HOMO (oder weniger, falls nicht genug vorhanden)
        int start_homo = std::max(0, homo_index - 4);
        std::cout << "Orbitale um HOMO:\n";
        for (int i = start_homo; i <= homo_index; ++i) {
            std::cout << "Orbital " << i + 1 << ": " << m_energies(i)
                      << " eV (Besetzung: " << (i < m_num_electrons / 2 ? 2 : 0) << ")";
            if (i == homo_index)
                std::cout << " [HOMO]";
            std::cout << "\n";
        }

        // 5 Orbitale nach LUMO (oder weniger, falls nicht genug vorhanden)
        int end_lumo = std::min(lumo_index + 4, int(m_energies.size()) - 1);
        std::cout << "\nOrbitale um LUMO:\n";
        for (int i = lumo_index; i <= end_lumo; ++i) {
            std::cout << "Orbital " << i + 1 << ": " << m_energies(i)
                      << " eV (Besetzung: " << (i < m_num_electrons / 2 ? 2 : 0) << ")";
            if (i == lumo_index)
                std::cout << " [LUMO]";
            std::cout << "\n";
        }

        std::cout << "\n";
    }
    if (m_verbose)
        std::cout << "Total electronic energy = " << energy << " Eh." << std::endl;

    // solveFockMatrix(S, H, m_energies, m_mo);

    return energy;

    // std::cout << diag_F.eigenvalues().sum();
}

Matrix EHT::MakeOverlap(Basisset& basisset)
{
    Matrix S = Eigen::MatrixXd::Zero(basisset.size(), basisset.size());

    // Aktualisiere die Atompositionen in allen Orbitalen
    for (int i = 0; i < basisset.size(); ++i) {
        basisset[i].x = m_mol.m_geometry(basisset[i].atom, 0) / 0.529177;
        basisset[i].y = m_mol.m_geometry(basisset[i].atom, 1) / 0.529177;
        basisset[i].z = m_mol.m_geometry(basisset[i].atom, 2) / 0.529177;
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
    double K = 1.75;
    Matrix H = Eigen::MatrixXd::Zero(basisset.size(), basisset.size());
    for (int i = 0; i < basisset.size(); ++i) {
        STO::Orbital bi = basisset[i];
        for (int j = 0; j < basisset.size(); ++j) {
            STO::Orbital bj = basisset[j];
            if (i == j)
                H(i, j) = bi.VSIP;
            else {
                H(i, j) = K * S(i, j) * (bi.VSIP + bj.VSIP) / 2.0;
                H(j, i) = K * S(j, i) * (bi.VSIP + bj.VSIP) / 2.0;
            }
        }
    }
    return H;
}
