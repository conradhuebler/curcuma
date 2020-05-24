/*
 * < C++ XTB Interface >
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

#include "external/xtb/include/xtb.h"

#include "src/core/global.h"
#include "src/core/molecule.h"

#include <assert.h>
#include <iostream>
#include <math.h>
#include <stdio.h>

#include "xtbinterface.h"

XTBInterface::XTBInterface()
{
}

double XTBInterface::GFNCalculation(const int* attyp, const double* coord, const int natoms, const double charge, int parameter, double* grad)
{
    double energy = 0;
#ifdef USE_XTB
    //double dipole[3];
    //double q[natoms];
    //double qp[6 * natoms];
    //double wbo[natoms * natoms];
    //char output;
    xtb_TEnvironment env;
    xtb_TMolecule mol;
    xtb_TCalculator calc;
    xtb_TResults res;

    env = xtb_newEnvironment();
    calc = xtb_newCalculator();
    res = xtb_newResults();
    mol = xtb_newMolecule(env, &natoms, attyp, coord, &charge, NULL, NULL, NULL);
    if (xtb_checkEnvironment(env)) {
        xtb_showEnvironment(env, NULL);
        return 1;
    }

    xtb_setVerbosity(env, XTB_VERBOSITY_MUTED);
    if (xtb_checkEnvironment(env)) {
        xtb_showEnvironment(env, NULL);
        return 2;
    }

    if (parameter == -1) {
        xtb_loadGFNFF(env, mol, calc, NULL);
        if (xtb_checkEnvironment(env)) {
            xtb_showEnvironment(env, NULL);
            return 3;
        }
    } else if (parameter == 0) {
        xtb_loadGFN0xTB(env, mol, calc, NULL);
        if (xtb_checkEnvironment(env)) {
            xtb_showEnvironment(env, NULL);
            return 3;
        }
    } else if (parameter == 1) {
        xtb_loadGFN1xTB(env, mol, calc, NULL);
        if (xtb_checkEnvironment(env)) {
            xtb_showEnvironment(env, NULL);
            return 3;
        }
    }

    else if (parameter == 2) {
        xtb_loadGFN2xTB(env, mol, calc, NULL);
        if (xtb_checkEnvironment(env)) {
            xtb_showEnvironment(env, NULL);
            return 3;
        }
    }

    xtb_singlepoint(env, mol, calc, res);
    if (xtb_checkEnvironment(env)) {
        xtb_showEnvironment(env, NULL);
        return 4;
    }

    xtb_getEnergy(env, res, &energy);
    if (grad != NULL)
        xtb_getGradient(env, res, grad);
#else
    throw("XTB is not included, sorry for that");
#endif
    return energy;
}

double XTBInterface::GFNCalculation(const Molecule& molecule, int parameter, double* grad)
{
    double energy = 0;
#ifdef USE_XTB

    int const natoms = molecule.AtomCount();
    double const charge = 0.0;

    int attyp[natoms];
    std::vector<int> atoms = molecule.Atoms();
    double coord[3 * natoms];

    for (int i = 0; i < atoms.size(); ++i) {
        std::pair<int, Position> atom = molecule.Atom(i);
        coord[3 * i + 0] = atom.second(0) / au;
        coord[3 * i + 1] = atom.second(1) / au;
        coord[3 * i + 2] = atom.second(2) / au;
        attyp[i] = atoms[i];
    }

    xtb_TEnvironment env;
    xtb_TMolecule mol;
    xtb_TCalculator calc;
    xtb_TResults res;

    env = xtb_newEnvironment();
    calc = xtb_newCalculator();
    res = xtb_newResults();
    mol = xtb_newMolecule(env, &natoms, attyp, coord, &charge, NULL, NULL, NULL);
    if (xtb_checkEnvironment(env)) {
        xtb_showEnvironment(env, NULL);
        return 1;
    }

    xtb_setVerbosity(env, XTB_VERBOSITY_MUTED);
    if (xtb_checkEnvironment(env)) {
        xtb_showEnvironment(env, NULL);
        return 2;
    }

    if (parameter == -1) {
        xtb_loadGFNFF(env, mol, calc, NULL);
        if (xtb_checkEnvironment(env)) {
            xtb_showEnvironment(env, NULL);
            return 3;
        }
    } else if (parameter == 0) {
        xtb_loadGFN0xTB(env, mol, calc, NULL);
        if (xtb_checkEnvironment(env)) {
            xtb_showEnvironment(env, NULL);
            return 3;
        }
    } else if (parameter == 1) {
        xtb_loadGFN1xTB(env, mol, calc, NULL);
        if (xtb_checkEnvironment(env)) {
            xtb_showEnvironment(env, NULL);
            return 3;
        }
    }

    else if (parameter == 2) {
        xtb_loadGFN2xTB(env, mol, calc, NULL);
        if (xtb_checkEnvironment(env)) {
            xtb_showEnvironment(env, NULL);
            return 3;
        }
    }
    xtb_singlepoint(env, mol, calc, res);
    if (xtb_checkEnvironment(env)) {
        xtb_showEnvironment(env, NULL);
        return 4;
    }

    xtb_getEnergy(env, res, &energy);
    if (grad != NULL)
        xtb_getGradient(env, res, grad);

#else
    throw("XTB is not included, sorry for that");
#endif
    return energy;
}
