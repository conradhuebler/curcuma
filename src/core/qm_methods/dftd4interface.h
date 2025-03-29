/*
 * < C++ XTB and tblite Interface >
 * Copyright (C) 2020 - 2024 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "src/tools/general.h"

#include "dftd_damping.h"
#include "dftd_dispersion.h"
#include "dftd_econv.h"
#include "dftd_geometry.h"

#include "interface/abstract_interface.h"

#include "src/core/global.h"

static json DFTD4Settings{
    { "d4_s6", 1.00 },
    { "d4_s8", 1.20065498 },
    { "d4_s10", 0 },
    { "d4_s9", 1 },
    { "d4_a1", 0.40085597 },
    { "d4_a2", 5.02928789 },
    { "d4_alp", 16 },
    { "d4_func", "pbe0" },
    { "d4_atm", true }
};

class DFTD4Interface : public QMInterface {
public:
    DFTD4Interface(const json& controller);
    DFTD4Interface();

    ~DFTD4Interface();

    bool InitialiseMolecule(const Mol& mol, double factor = 1);
    bool InitialiseMolecule(const std::vector<int>& atomtype);

    double Calculation(bool gradient = 0, bool verbose = false) override;
    void UpdateParameters(const json& controller);

    void clear();
    dftd4::dparam Parameter() const { return m_par; }

    inline void UpdateGeometry(const double* coord)
    {
        for (int i = 0; i < m_natoms; ++i) {
            m_mol.xyz(i, 0) = coord[3 * i + 0];
            m_mol.xyz(i, 1) = coord[3 * i + 1];
            m_mol.xyz(i, 2) = coord[3 * i + 2];
        }
    }

    inline void UpdateAtom(int index, double x, double y, double z, int element = -1)
    {
        m_mol.UpdateAtom(index, x, y, z, element);
    }

    void PrintParameter() const;

private:
    dftd4::dparam m_par;
    dftd4::TMolecule m_mol;
    int m_charge = 0;
    int m_natoms;
};
