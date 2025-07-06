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

QMDriver::QMDriver()
{
}

bool QMDriver::InitialiseMolecule()
{
    m_num_electrons = 0;
    m_mo = Matrix::Zero(m_mol.m_number_atoms, m_mol.m_number_atoms);
    m_energies = Vector::Zero(m_mol.m_number_atoms);
    return true;
}

Matrix QMDriver::MakeOverlap(Basisset& basisset)
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

Matrix QMDriver::MakeH(const Matrix& S, const Basisset& basisset)
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
