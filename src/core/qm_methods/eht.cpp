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
    if (m_verbose)
        std::cout << basisset.size() << std::endl;
    Matrix S = MakeOverlap(basisset);
    if (m_verbose)
        std::cout << S << std::endl;
    Matrix H = MakeH(S, basisset);
    if (m_verbose) {
        std::cout << std::endl
                  << std::endl;
        std::cout << H << std::endl;
    }
    Matrix S_1_2 = Matrix::Zero(basisset.size(), basisset.size());
    Eigen::JacobiSVD<Matrix> svd(S, Eigen::ComputeThinU | Eigen::ComputeThinV);

    Matrix IntS;
    IntS.setZero(basisset.size(), basisset.size());

    for (int i = 0; i < basisset.size(); ++i)
        IntS(i, i) = 1 / sqrt(svd.singularValues()(i));

    S_1_2 = svd.matrixU() * IntS * svd.matrixV().transpose();

    Eigen::SelfAdjointEigenSolver<Matrix> diag_F;
    Matrix F = S_1_2.transpose() * H * S_1_2;
    diag_F.compute(F);
    std::cout << std::endl;
    m_energies = diag_F.eigenvalues();
    m_mo = diag_F.eigenvectors();

    for (int i = 0; i < m_num_electrons / 2; ++i) {
        energy += diag_F.eigenvalues()(i) * 2;
        if (m_verbose)
            std::cout << diag_F.eigenvalues()(i) << " 2" << std::endl;
    }
    if (m_verbose)
        std::cout << "Total electronic energy = " << energy << " Eh." << std::endl;
    return energy;

    // std::cout << diag_F.eigenvalues().sum();
}

Matrix EHT::MakeOverlap(Basisset& basisset)
{
    Matrix S = Eigen::MatrixXd::Zero(basisset.size(), basisset.size());

    // Aktualisiere die Atompositionen in allen Orbitalen
    for (int i = 0; i < basisset.size(); ++i) {
        basisset[i].x = m_mol.m_geometry(basisset[i].atom, 0);
        basisset[i].y = m_mol.m_geometry(basisset[i].atom, 1);
        basisset[i].z = m_mol.m_geometry(basisset[i].atom, 2);
    }

    // Berechne die Überlappungsmatrix direkt
    for (int i = 0; i < basisset.size(); ++i) {
        for (int j = 0; j <= i; ++j) {
            // Die verbesserte calculateOverlap Funktion gibt bereits
            // normalisierte Werte zurück
            double overlap = STO::calculateOverlap(basisset[i], basisset[j]);
            S(i, j) = S(j, i) = overlap;

            // Debug-Ausgabe für kritische Werte
            if (m_verbose) {
                std::cout << "Overlap between orbital " << i << " and " << j
                          << ": " << overlap << std::endl;
            }
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
