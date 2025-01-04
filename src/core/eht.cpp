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

#include "src/core/global.h"
#include "src/core/molecule.h"

#include "external/CxxThreadPool/include/CxxThreadPool.h"

#include <set>
#include <vector>

#include <Eigen/Dense>

#include "json.hpp"

#include "eht.h"

#include <iomanip>

EHT::EHT()
{
}

void EHT::CalculateEHT(bool gradient, bool verbose)
{
    m_verbose = verbose;
    if (gradient) {
        std::cout << "EHT does not support gradients." << std::endl;
    }

    if (m_molecule.AtomCount() == 0) {
        std::cout << "No molecule set." << std::endl;
        return;
    }

    start();
}

std::vector<STO_6G> EHT::MakeBasis()
{
    std::vector<STO_6G> basisset;
    for (int i = 0; i < m_molecule.Atoms().size(); ++i) {
        if (m_molecule.Atom(i).first == 1) {
            m_num_electrons += 1;
            STO_6G b;
            b.index = i;

            ehtSTO_6GHs t;
            b.alpha = t.alpha;
            b.coeff = t.coeff;
            b.sym = "s";
            b.e = -0.5;
            basisset.push_back(b);
            std::cout << "H" << std::endl;
        } else if (m_molecule.Atom(i).first == 6) {
            m_num_electrons += 4;

            STO_6G b;
            b.index = i;
            b.sym = "s";
            ehtSTO_6GCs t;
            b.alpha = t.alpha;
            b.coeff = t.coeff;
            b.e = -0.7144;

            basisset.push_back(b);
            ehtSTO_6GCp t2;
            b.alpha = t2.alpha;
            b.coeff = t2.coeff;
            b.x = 1;
            b.sym = "pz";
            b.e = -0.3921;

            basisset.push_back(b);

            b.x = 0;
            b.y = 1;
            b.sym = "px";
            basisset.push_back(b);

            b.y = 0;
            b.z = 1;
            b.sym = "py";
            basisset.push_back(b);
            std::cout << "C" << std::endl;
        } else if (m_molecule.Atom(i).first == 8) {
            m_num_electrons += 6;

            STO_6G b;
            b.index = i;
            b.sym = "s";
            ehtSTO_6GOs t;
            b.alpha = t.alpha;
            b.coeff = t.coeff;
            b.e = -1.1904;

            basisset.push_back(b);
            ehtSTO_6GOp t2;
            b.alpha = t2.alpha;
            b.coeff = t2.coeff;
            b.x = 1;
            b.sym = "pz";
            b.e = -0.5827;

            basisset.push_back(b);

            b.x = 0;
            b.y = 1;
            b.sym = "px";
            basisset.push_back(b);

            b.y = 0;
            b.z = 1;
            b.sym = "py";
            basisset.push_back(b);
            std::cout << "O" << std::endl;
        } else if (m_molecule.Atom(i).first == 7) {
            m_num_electrons += 5;

            STO_6G b;
            b.index = i;
            b.sym = "s";
            ehtSTO_6GNs t;
            b.alpha = t.alpha;
            b.coeff = t.coeff;
            b.e = -0.9404;

            basisset.push_back(b);
            ehtSTO_6GNp t2;
            b.alpha = t2.alpha;
            b.coeff = t2.coeff;
            b.x = 1;
            b.sym = "pz";
            b.e = -0.5110;

            basisset.push_back(b);

            b.x = 0;
            b.y = 1;
            b.sym = "px";
            basisset.push_back(b);

            b.y = 0;
            b.z = 1;
            b.sym = "py";
            basisset.push_back(b);
            std::cout << "O" << std::endl;
        }
    }
    return basisset;
}

void EHT::start()
{
    m_molecule.print_geom();
    auto basisset = MakeBasis();
    if (m_verbose)
        std::cout << basisset.size() << std::endl;
    Matrix S = MakeOverlap(basisset);
    // std::cout << S << std::endl;
    Matrix H = MakeH(S, basisset);
    if (m_verbose) {
        std::cout << std::endl
                  << std::endl;
        // std::cout << H << std::endl;
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

    double energy = 0;
    for (int i = 0; i < m_num_electrons / 2; ++i) {
        energy += diag_F.eigenvalues()(i) * 2;
        if (m_verbose)
            std::cout << diag_F.eigenvalues()(i) << " 2" << std::endl;
    }
    std::cout << "Total electronic energy = " << energy << " Eh." << std::endl;

    // std::cout << diag_F.eigenvalues().sum();
}

Matrix EHT::MakeOverlap(const std::vector<STO_6G>& basisset)
{
    Matrix S = Eigen::MatrixXd::Zero(basisset.size(), basisset.size());

    std::vector<double> norm(basisset.size(), 0);

    for (int i = 0; i < basisset.size(); ++i) {
        STO_6G bi = basisset[i];
        Vector pi = m_molecule.Atom(bi.index).second;
        STO_6G bj = basisset[i];
        Vector pj = m_molecule.Atom(bi.index).second;
        double xx = 0;
        if (bi.sym.compare("s") == 0 && bj.sym.compare("s") == 0) {
            for (int c1 = 0; c1 < 6; ++c1)
                for (int c2 = 0; c2 < 6; ++c2)
                    xx += ss(pi(0), pj(0), pi(1), pj(1), pi(2), pj(2), bi.alpha[c1], bj.alpha[c1], bi.coeff[c2], bj.coeff[c2]);
        } else if ((bi.sym.compare("px") == 0 && bj.sym.compare("px") == 0)) {
            for (int c1 = 0; c1 < 6; ++c1)
                for (int c2 = 0; c2 < 6; ++c2)
                    xx += p2(pi(0), pj(0), pi(1), pj(1), pi(2), pj(2), bi.alpha[c1], bj.alpha[c1], bi.coeff[c2], bj.coeff[c2]);
        } else if ((bi.sym.compare("py") == 0 && bj.sym.compare("py") == 0)) {
            for (int c1 = 0; c1 < 6; ++c1)
                for (int c2 = 0; c2 < 6; ++c2)
                    xx += p2(pi(1), pj(1), pi(0), pj(0), pi(2), pj(2), bi.alpha[c1], bj.alpha[c1], bi.coeff[c2], bj.coeff[c2]);
        } else if ((bi.sym.compare("pz") == 0 && bj.sym.compare("pz") == 0)) {
            for (int c1 = 0; c1 < 6; ++c1)
                for (int c2 = 0; c2 < 6; ++c2)
                    xx += p2(pi(2), pj(2), pi(1), pj(1), pi(0), pj(0), bi.alpha[c1], bj.alpha[c1], bi.coeff[c2], bj.coeff[c2]);
        }
        norm[i] = 1 / sqrt(xx);
    }

    for (int i = 0; i < basisset.size(); ++i) {
        STO_6G bi = basisset[i];
        Vector pi = m_molecule.Atom(bi.index).second;
        for (int j = i; j < basisset.size(); ++j) {
            STO_6G bj = basisset[j];
            Vector pj = m_molecule.Atom(bj.index).second;

            double xx = 0;

            if (bi.sym.compare("s") == 0 && bj.sym.compare("s") == 0) { /* s-s orbitals */
                for (int c1 = 0; c1 < 6; ++c1)
                    for (int c2 = 0; c2 < 6; ++c2)
                        xx += ss(pi(0), pj(0), pi(1), pj(1), pi(2), pj(2), bi.alpha[c1], bj.alpha[c1], bi.coeff[c2], bj.coeff[c2]);
            } else if ((bi.sym.compare("s") == 0 && bj.sym.compare("px") == 0) || (bi.sym.compare("px") == 0 && bj.sym.compare("s") == 0)) { /* s-px orbitals */
                if (bj.sym.compare("px") == 0) {
                    for (int c1 = 0; c1 < 6; ++c1)
                        for (int c2 = 0; c2 < 6; ++c2)
                            xx += sp(pi(0), pj(0), pi(1), pj(1), pi(2), pj(2), bi.alpha[c1], bj.alpha[c1], bi.coeff[c2], bj.coeff[c2]);
                } else if (bi.sym.compare("px") == 0) {
                    for (int c1 = 0; c1 < 6; ++c1)
                        for (int c2 = 0; c2 < 6; ++c2)
                            xx += sp(pj(0), pi(0), pj(1), pi(1), pj(2), pi(2), bj.alpha[c2], bi.alpha[c2], bj.coeff[c1], bi.coeff[c1]);
                }
            } else if ((bi.sym.compare("s") == 0 && bj.sym.compare("py") == 0) || (bi.sym.compare("py") == 0 && bj.sym.compare("s") == 0)) { /* s-py orbitals */
                if (bj.sym.compare("py") == 0) {
                    for (int c1 = 0; c1 < 6; ++c1)
                        for (int c2 = 0; c2 < 6; ++c2)
                            xx += sp(pi(0), pj(1), pi(1), pj(0), pi(2), pj(2), bi.alpha[c1], bj.alpha[c1], bi.coeff[c2], bj.coeff[c2]);
                } else if (bi.sym.compare("py") == 0) {
                    for (int c1 = 0; c1 < 6; ++c1)
                        for (int c2 = 0; c2 < 6; ++c2)
                            xx += sp(pj(1), pi(0), pj(0), pi(1), pj(2), pi(2), bj.alpha[c2], bi.alpha[c2], bj.coeff[c1], bi.coeff[c1]);
                }
            } else if ((bi.sym.compare("s") == 0 && bj.sym.compare("pz") == 0) || (bi.sym.compare("pz") == 0 && bj.sym.compare("s") == 0)) { /* s-pz orbitals */
                if (bj.sym.compare("pz") == 0) {
                    for (int c1 = 0; c1 < 6; ++c1)
                        for (int c2 = 0; c2 < 6; ++c2)
                            xx += sp(pi(0), pj(2), pi(1), pj(1), pi(2), pj(0), bi.alpha[c1], bj.alpha[c1], bi.coeff[c2], bj.coeff[c2]);
                } else if (bi.sym.compare("pz") == 0) {
                    for (int c1 = 0; c1 < 6; ++c1)
                        for (int c2 = 0; c2 < 6; ++c2)
                            xx += sp(pj(2), pi(0), pj(1), pi(1), pj(0), pi(2), bj.alpha[c2], bi.alpha[c2], bj.coeff[c1], bi.coeff[c1]);
                }
            } else if ((bi.sym.compare("px") == 0 && bj.sym.compare("px") == 0)) { /* px-px orbitals */
                for (int c1 = 0; c1 < 6; ++c1)
                    for (int c2 = 0; c2 < 6; ++c2)
                        xx += p2(pi(0), pj(0), pi(1), pj(1), pi(2), pj(2), bi.alpha[c1], bj.alpha[c1], bi.coeff[c2], bj.coeff[c2]);
            } else if ((bi.sym.compare("py") == 0 && bj.sym.compare("py") == 0)) { /* py-py orbitals */
                for (int c1 = 0; c1 < 6; ++c1)
                    for (int c2 = 0; c2 < 6; ++c2)
                        xx += p2(pi(1), pj(1), pi(0), pj(0), pi(2), pj(2), bi.alpha[c1], bj.alpha[c1], bi.coeff[c2], bj.coeff[c2]);
            } else if ((bi.sym.compare("pz") == 0 && bj.sym.compare("pz") == 0)) { /* pz-pz orbitals */
                for (int c1 = 0; c1 < 6; ++c1)
                    for (int c2 = 0; c2 < 6; ++c2)
                        xx += p2(pi(2), pj(2), pi(1), pj(1), pi(0), pj(0), bi.alpha[c1], bj.alpha[c1], bi.coeff[c2], bj.coeff[c2]);
            } else if ((bi.sym.compare("px") == 0 && bj.sym.compare("py") == 0) || (bi.sym.compare("py") == 0 && bj.sym.compare("px") == 0)) { /* px-py orbitals */
                if (bj.sym.compare("py") == 0) {
                    for (int c1 = 0; c1 < 6; ++c1)
                        for (int c2 = 0; c2 < 6; ++c2)
                            xx += pp(pi(0), pj(1), pi(1), pj(0), pi(2), pj(2), bi.alpha[c1], bj.alpha[c1], bi.coeff[c2], bj.coeff[c2]);
                } else if (bi.sym.compare("py") == 0) {
                    for (int c1 = 0; c1 < 6; ++c1)
                        for (int c2 = 0; c2 < 6; ++c2)
                            xx += pp(pi(1), pj(0), pi(0), pj(1), pi(2), pj(2), bi.alpha[c1], bj.alpha[c1], bi.coeff[c2], bj.coeff[c2]);
                }
            } else if ((bi.sym.compare("px") == 0 && bj.sym.compare("pz") == 0) || (bi.sym.compare("pz") == 0 && bj.sym.compare("px") == 0)) { /* px-pz-orbitals */
                if (bj.sym.compare("pz") == 0) {
                    for (int c1 = 0; c1 < 6; ++c1)
                        for (int c2 = 0; c2 < 6; ++c2)
                            xx += pp(pi(0), pj(2), pi(1), pj(1), pi(2), pj(0), bi.alpha[c1], bj.alpha[c1], bi.coeff[c2], bj.coeff[c2]);
                } else if (bi.sym.compare("pz") == 0) {
                    for (int c1 = 0; c1 < 6; ++c1)
                        for (int c2 = 0; c2 < 6; ++c2)
                            xx += pp(pi(2), pj(0), pi(1), pj(1), pi(0), pj(2), bi.alpha[c1], bj.alpha[c1], bi.coeff[c2], bj.coeff[c2]);
                }
            } else if ((bi.sym.compare("py") == 0 && bj.sym.compare("pz") == 0)) { /* py-pz-orbitals */

                for (int c1 = 0; c1 < 6; ++c1)
                    for (int c2 = 0; c2 < 6; ++c2)
                        xx += pp(pi(1), pj(2), pi(0), pj(1), pi(2), pj(0), bi.alpha[c1], bj.alpha[c1], bi.coeff[c2], bj.coeff[c2]);

            } else if (bi.sym.compare("pz") == 0 && bj.sym.compare("py") == 0) { /* pz-py-orbitals */
                for (int c1 = 0; c1 < 6; ++c1)
                    for (int c2 = 0; c2 < 6; ++c2)
                        xx += pp(pi(2), pj(1), pi(1), pj(0), pi(0), pj(2), bi.alpha[c1], bj.alpha[c1], bi.coeff[c2], bj.coeff[c2]);
            } else
                std::cout << "something is missing " << bi.sym << " " << bj.sym << std::endl;
            S(i, j) = xx * norm[i] * norm[j];
            S(j, i) = xx * norm[j] * norm[i];

            // S(j,i) += xx*norm[j]*norm[i]*0.5;

            // std::cout << i << " " << j << " "  << bi.sym << " " << bj.sym << " " << " "<< xx << " " <<xx*norm[i]*norm[j] <<std::endl;
        }
    }
    return S;
}

Matrix EHT::MakeH(const Matrix& S, const std::vector<STO_6G>& basisset)
{
    double K = 1.75;
    Matrix H = Eigen::MatrixXd::Zero(basisset.size(), basisset.size());
    for (int i = 0; i < basisset.size(); ++i) {
        STO_6G bi = basisset[i];
        Vector pi = m_molecule.Atom(bi.index).second;
        for (int j = 0; j < basisset.size(); ++j) {
            STO_6G bj = basisset[j];
            Vector pj = m_molecule.Atom(bj.index).second;
            if (i == j)
                H(i, j) = bi.e;
            else {
                H(i, j) = K * S(i, j) * (bi.e + bj.e) / 2.0;
                H(j, i) = K * S(j, i) * (bi.e + bj.e) / 2.0;
            }
        }
    }
    return H;
}

double EHT::ss(double x1, double x2, double y1, double y2, double z1, double z2, double alpha1, double alpha2, double c1, double c2)
{
    double ex = exp(-((alpha1 * alpha2) / (alpha1 + alpha2)) * ((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2)));
    return (c1 * c2 * (pow(2, 1.5) * pow(alpha1 * alpha2, 0.75) * pow(alpha1 + alpha2, -1.5) * ex));
}

double EHT::sp(double x1, double x2, double y1, double y2, double z1, double z2, double alpha1, double alpha2, double c1, double c2)
{
    double ex = exp(-((alpha1 * alpha2) / (alpha1 + alpha2)) * ((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2)));
    return ((c1 * c2 * (4 * pow(2, 0.5) * pow(alpha1, 1.75) * pow(alpha2, 1.25) * ex * (x1 - x2))) * (pow(alpha1 + alpha2, -2.5)));
}

double EHT::pp(double x1, double x2, double y1, double y2, double z1, double z2, double alpha1, double alpha2, double c1, double c2)
{
    double ex = exp(-((alpha1 * alpha2) / (alpha1 + alpha2)) * ((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2)));
    return (c1 * c2 * ((-8 * (pow(2, 0.5) * (pow(alpha1 * alpha2, 2.25)) * ex * (x1 - x2) * (y1 - y2))) * (pow(alpha1 + alpha2, -3.5))));
}

double EHT::p2(double x1, double x2, double y1, double y2, double z1, double z2, double alpha1, double alpha2, double c1, double c2)
{
    double ex = exp(-((alpha1 * alpha2) / (alpha1 + alpha2)) * ((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2)));
    return (c1 * c2 * ((4 * pow(2, 0.5) * pow(alpha1 * alpha2, 1.25) * ex * (alpha2 + alpha1 * (1 - 2 * alpha2 * (x1 - x2) * (x1 - x2))) * (pow(alpha1 + alpha2, -3.5)))));
}
