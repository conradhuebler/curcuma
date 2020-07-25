/*
 * <Docking tool for structures. >
 * Copyright (C) 2020 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "external/xtb/include/xtb.h"

#include "src/core/molecule.h"

class XTBInterface {
public:
    XTBInterface();
    ~XTBInterface();

    bool InitialiseMolecule(const Molecule& molecule);
    bool InitialiseMolecule(const int* attyp, const double* coord, const int natoms, const double charge);

    bool UpdateMolecule(const Molecule& molecule);
    bool UpdateMolecule(const double* coord);

    /* int parameter
     * 66 = xtb GFN FF
     * 0 = xtb GFN 0
     * 1 = xtb GFN 1
     * 2 = xtb GFN 2
     * */
    double GFNCalculation(int parameter = 2, double* grad = 0);

private:
#ifdef USE_XTB
    double m_thr = 1.0e-10;
    xtb_TEnvironment m_env;
    xtb_TMolecule m_mol;
    xtb_TCalculator m_calc;
    xtb_TResults m_res;
#endif
};
