/*
 * < C++ XTB and tblite Interface >
 * Copyright (C) 2020 - 2023 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "external/xtb/include/xtb.h"

#include "src/core/molecule.h"

static json XTBSettings{
    { "xtb_ac", 2 },
    { "xtb_maxiter", 100 },
    { "xtb_temp", 298 }
};

class UFF;

class XTBInterface {
public:
    XTBInterface(const json& xtbsettings = XTBSettings);
    ~XTBInterface();

    bool InitialiseMolecule(const Molecule& molecule);
    bool InitialiseMolecule(const Molecule* molecule);
    bool InitialiseMolecule(const int* attyp, const double* coord, const int natoms, const double charge, const int spin);

    bool UpdateMolecule(const Molecule& molecule);
    bool UpdateMolecule(const double* coord);

    /* int parameter
     * 66 = xtb GFN FF
     * 0 = xtb GFN 0
     * 1 = xtb GFN 1
     * 2 = xtb GFN 2
     * */
    double GFNCalculation(int parameter = 2, double* grad = 0);

    void clear();
    std::vector<double> Charges() const;
    std::vector<double> Dipole() const;

    std::vector<std::vector<double>> BondOrders() const;

private:
    int m_atomcount = 0;
    double m_thr = 1.0e-10;
    double* m_coord;
    int* m_attyp;
    int m_accuracy = 1, m_maxiter = 100;
    double m_temp = 298;
    xtb_TEnvironment m_env = NULL;
    xtb_TMolecule m_xtb_mol = NULL;
    xtb_TCalculator m_xtb_calc = NULL;
    xtb_TResults m_xtb_res = NULL;
    bool m_initialised = false;
    json m_xtbsettings;
};
