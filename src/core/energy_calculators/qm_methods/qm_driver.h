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

#include "external/CxxThreadPool/include/CxxThreadPool.hpp"

#include <set>
#include <vector>

#include <Eigen/Dense>

#include "GTOIntegrals.hpp"
#include "STOIntegrals.hpp"

#include "interface/abstract_interface.h"

#include "json.hpp"

typedef std::vector<STO::Orbital> Basisset;

class QMDriver : public QMInterface {
public:
    QMDriver();

    virtual bool InitialiseMolecule() override;
    virtual double Calculation(bool gradient = false, bool verbose = false) = 0;

    Matrix MolecularOrbitals() const { return m_mo; }
    Vector Energies() const { return m_energies; }
    int NumElectrons() const { return m_num_electrons; }

    double TotalEnergy() const
    {
        return m_total_energy;
    }

private:
    // std::vector<STO_6G> MakeBasis();
    // Basisset MakeBasis();

    virtual Matrix MakeOverlap(Basisset& basisset) = 0;
    virtual Matrix MakeH(const Matrix& S, const Basisset& basisset) = 0;

protected:
    double m_total_energy = 0;
    // Mol m_mol;
    Matrix m_H, m_S;
    int m_num_electrons = 0;
    int m_threads = 4;
    Matrix m_mo;
    Vector m_energies;

    bool m_verbose = false;
};
