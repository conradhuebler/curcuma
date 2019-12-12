/*
 * <RMSD calculator for chemical structures.>
 * Copyright (C) 2019  Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "src/core/molecule.h"
#include "src/tools/geometry.h"

#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <vector>

#include "rmsd.h"

RMSDDriver::RMSDDriver(const Molecule &reference, const Molecule &target)
    : m_reference(reference)
    , m_target(target)
{

}


double RMSDDriver::CalculateRMSD()
{
    /* The rmsd kabsch algorithmn was adopted from here:
     * https://github.com/oleg-alexandrov/projects/blob/master/eigen/Kabsch.cpp
     * The specific git commit was
     * https://github.com/oleg-alexandrov/projects/blob/e7b1eb7a4d83d41af563c24859072e4ddd9b730b/eigen/Kabsch.cpp
     */
    double rmsd = 0;
    int factor = 1;

    Geometry reference = CenterMolecule(m_reference);
    Geometry target = CenterMolecule(m_target);

    Eigen::MatrixXd ref = reference.transpose();
    Eigen::MatrixXd tar = target.transpose();


    Eigen::MatrixXd Cov = ref*tar.transpose();
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(Cov, Eigen::ComputeThinU | Eigen::ComputeThinV);

     double d = (svd.matrixV() * svd.matrixU().transpose()).determinant();
           if (d > 0)
               d = factor*1.0;
           else
               d = factor*-1.0;
           Eigen::Matrix3d I = Eigen::Matrix3d::Identity(3, 3);
           I(2, 2) = d;



    const Eigen::Matrix3d R =  svd.matrixV() * I * svd.matrixU().transpose();

    std::cout << std::endl << R << std::endl << std::endl;

    Geometry rotated = tar.transpose()*R;

    for(int i = 0; i < rotated.rows(); ++i)
    {
        rmsd += (rotated(i, 0) - reference(i, 0))*(rotated(i, 0) - reference(i, 0)) +
                      (rotated(i, 1) - reference(i, 1))*(rotated(i, 1) - reference(i, 1)) +
                      (rotated(i, 2) - reference(i, 2))*(rotated(i, 2) - reference(i, 2));
    }
    rmsd /= double(rotated.rows());

    return sqrt(rmsd);
}

Geometry RMSDDriver::CenterMolecule(const Molecule &mol) const
{
    return GeometryTools::TranslateMolecule(mol, mol.Centroid(), Position{0, 0, 0});
}
