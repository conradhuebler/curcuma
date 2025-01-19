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

#pragma once
#include "src/core/global.h"

#include "external/CxxThreadPool/include/CxxThreadPool.h"

#include <set>
#include <vector>

#include <Eigen/Dense>

#include "json.hpp"

struct ehtSTO_6GHs {
    std::vector<double> alpha = { 0.3552322122E+02, 0.6513143725E+01, 0.1822142904E+01, 0.6259552659E+00, 0.2430767471E+00, 0.1001124280E+00 };
    std::vector<double> coeff = { 0.9163596281E-02, 0.4936149294E-01, 0.1685383049E+00, 0.3705627997E+00, 0.4164915298E+00, 0.1303340841E+00 };
};

struct ehtSTO_6GCs {
    std::vector<double> alpha = { 0.3049723950E+02, 0.6036199601E+01, 0.1876046337E+01, 0.7217826470E+00, 0.3134706954E+00, 0.1436865550E+00 };
    std::vector<double> coeff = { -0.1325278809E-01, -0.4699171014E-01, -0.3378537151E-01, 0.2502417861E+00, 0.5951172526E+00, 0.2407061763E+00 };
};

struct ehtSTO_6GCp {
    std::vector<double> alpha = { 0.3049723950E+02, 0.6036199601E+01, 0.1876046337E+01, 0.7217826470E+00, 0.3134706954E+00, 0.1436865550E+00 };
    std::vector<double> coeff = { 0.3759696623E-02, 0.3767936984E-01, 0.1738967435E+00, 0.4180364347E+00, 0.4258595477E+00, 0.1017082955E+00 };
};

struct ehtSTO_6GNs {
    std::vector<double> alpha = { 0.3919880787E+02, 0.7758467071E+01, 0.2411325783E+01, 0.9277239437E+00, 0.4029111410E+00, 0.1846836552E+00 };
    std::vector<double> coeff = { -0.1325278809E-01, -0.4699171014E-01, -0.3378537151E-01, 0.2502417861E+00, 0.5951172526E+00, 0.2407061763E+00 };
};

struct ehtSTO_6GNp {
    std::vector<double> alpha = { 0.3919880787E+02, 0.7758467071E+01, 0.2411325783E+01, 0.9277239437E+00, 0.4029111410E+00, 0.1846836552E+00 };
    std::vector<double> coeff = { 0.3759696623E-02, 0.3767936984E-01, 0.1738967435E+00, 0.4180364347E+00, 0.4258595477E+00, 0.1017082955E+00 };
};

struct ehtSTO_6GOs {
    std::vector<double> alpha = { 0.5218776196E+02, 0.1032932006E+02, 0.3210344977E+01, 0.1235135428E+01, 0.5364201581E+00, 0.2458806060E+00 };
    std::vector<double> coeff = { -0.1325278809E-01, -0.4699171014E-01, -0.3378537151E-01, 0.2502417861E+00, 0.5951172526E+00, 0.2407061763E+00 };
};

struct ehtSTO_6GOp {
    std::vector<double> alpha = { 0.5218776196E+02, 0.1032932006E+02, 0.3210344977E+01, 0.1235135428E+01, 0.5364201581E+00, 0.2458806060E+00 };
    std::vector<double> coeff = { 0.3759696623E-02, 0.3767936984E-01, 0.1738967435E+00, 0.4180364347E+00, 0.4258595477E+00, 0.1017082955E+00 };
};

struct STO_6G {
    std::vector<double> alpha, coeff;
    int x = 0, y = 0, z = 0;
    int index = 0;
    std::string sym = "";
    double e = 0;
};

class EHT {
public:
    EHT();

    void setMolecule(const Mol& mol) { m_mol = mol; }
    void Initialise();
    void Calculate(double* gradient = nullptr, bool verbose = false);

    Matrix MolecularOrbitals() const { return m_mo; }
    Vector Energies() const { return m_energies; }
    int NumElectrons() const { return m_num_electrons; }

private:
    std::vector<STO_6G> MakeBasis();

    Mol m_mol;
    /* Some integrals */
    /* s - s - sigma bond */
    double ss(double x1, double x2, double y1, double y2, double z1, double z2, double alpha1, double alpha2, double c1, double c2);

    /* s - p - sigma bond */
    double sp(double x1, double x2, double y1, double y2, double z1, double z2, double alpha1, double alpha2, double c1, double c2);

    /* p - p - pi bond */
    double pp(double x1, double x2, double y1, double y2, double z1, double z2, double alpha1, double alpha2, double c1, double c2);

    /* p - p - sigma bond */
    double p2(double x1, double x2, double y1, double y2, double z1, double z2, double alpha1, double alpha2, double c1, double c2);

    Matrix MakeOverlap(const std::vector<STO_6G>& basisset);
    Matrix MakeH(const Matrix& S, const std::vector<STO_6G>& basisset);

    int m_num_electrons = 0;

    Matrix m_mo;
    Vector m_energies;

    bool m_verbose = false;
};
