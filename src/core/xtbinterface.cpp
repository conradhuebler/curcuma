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

double XTBInterface::GFN0Energy(const Molecule& molecule)
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

    xtb::PEEQ_options const opt = (xtb::PEEQ_options){ 0, 1, 2.0, 300, false, false, "none" };

    double dipole[3];
    double q[natoms];
    double qp[6 * natoms];
    double wbo[natoms * natoms];
    char output;
    double* grad = 0;
    int stat = xtb::GFN0_calculation(&natoms, attyp, &charge, NULL, coord, &opt, &output, &energy, grad);
#else
    throw("XTB is not included, sorry for that");
#endif
    return energy;
}

double XTBInterface::GFN1Energy(const Molecule& molecule)
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

    xtb::SCC_options const opt = (xtb::SCC_options){ 0, 0, 1.0, 300.0, true, true, true, 30, "none" };

    double dipole[3];
    double q[natoms];
    double qp[6 * natoms];
    double wbo[natoms * natoms];
    char output;
    int stat = xtb::GFN1_calculation(&natoms, attyp, &charge, NULL, coord, &opt, &output,
        &energy, NULL, dipole, q, NULL);
#else
    throw("XTB is not included, sorry for that");
#endif
    return energy;
}

double XTBInterface::GFN2Energy(const Molecule& molecule)
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

    xtb::SCC_options const opt = (xtb::SCC_options){ 0, 1, 2.0, 300.0, false, false, true, 30, "none" };

    double dipole[3];
    double q[natoms];
    double qp[6 * natoms];
    double wbo[natoms * natoms];
    char output;
    int stat = xtb::GFN2_calculation(&natoms, attyp, &charge, NULL, coord, &opt, &output,
        &energy, NULL, dipole, q, NULL, qp, wbo);
#else
    throw("XTB is not included, sorry for that");
#endif
    return energy;
}

double XTBInterface::GFN2Gradient(const int* attyp, const double* coord, const int natoms, const double charge, double* grad)
{
    double energy = 0;
#ifdef USE_XTB
    xtb::SCC_options const opt = (xtb::SCC_options){ 0, 1, 2.0, 300.0, true, true, true, 30, "none" };

    double dipole[3];
    double q[natoms];
    double qp[6 * natoms];
    double wbo[natoms * natoms];
    char output;
    int stat = xtb::GFN2_calculation(&natoms, attyp, &charge, NULL, coord, &opt, &output,
        &energy, grad, dipole, q, NULL, qp, wbo);
#else
    throw("XTB is not included, sorry for that");
#endif
    return energy;
}

double XTBInterface::GFN2Energy(const int* attyp, const double* coord, const int natoms, const double charge)
{
    double energy = 0;
#ifdef USE_XTB
    xtb::SCC_options const opt = (xtb::SCC_options){ 0, 1, 2.0, 300.0, false, false, true, 30, "none" };

    double dipole[3];
    double q[natoms];
    double qp[6 * natoms];
    double wbo[natoms * natoms];
    char output;
    double* grad = 0;
    int stat = xtb::GFN2_calculation(&natoms, attyp, &charge, NULL, coord, &opt, &output,
        &energy, grad, dipole, q, NULL, qp, wbo);
#else
    throw("XTB is not included, sorry for that");
#endif
    return energy;
}

double XTBInterface::GFN1Gradient(const int* attyp, const double* coord, const int natoms, const double charge, double* grad)
{
    double energy = 0;
#ifdef USE_XTB
    xtb::SCC_options const opt = (xtb::SCC_options){ 0, 1, 2.0, 300.0, true, false, true, 30, "none" };

    double dipole[3];
    double q[natoms];
    double qp[6 * natoms];
    double wbo[natoms * natoms];
    char output;
    int stat = xtb::GFN1_calculation(&natoms, attyp, &charge, NULL, coord, &opt, &output,
        &energy, grad, dipole, q, wbo);
#else
    throw("XTB is not included, sorry for that");
#endif
    return energy;
}

double XTBInterface::GFN1Energy(const int* attyp, const double* coord, const int natoms, const double charge)
{
    double energy = 0;
#ifdef USE_XTB
    xtb::SCC_options const opt = (xtb::SCC_options){ 0, 1, 2.0, 300.0, false, false, true, 30, "none" };

    double dipole[3];
    double q[natoms];
    double qp[6 * natoms];
    double wbo[natoms * natoms];
    char output;
    double* grad = 0;
    int stat = xtb::GFN1_calculation(&natoms, attyp, &charge, NULL, coord, &opt, &output,
        &energy, grad, dipole, q, wbo);
#else
    throw("XTB is not included, sorry for that");
#endif
    return energy;
}

double XTBInterface::GFN0Gradient(const int* attyp, const double* coord, const int natoms, const double charge, double* grad)
{
    double energy = 0;
#ifdef USE_XTB
    xtb::PEEQ_options const opt = (xtb::PEEQ_options){ 0, 1, 2.0, 300, true, false, "none" };

    double dipole[3];
    double q[natoms];
    double qp[6 * natoms];
    double wbo[natoms * natoms];
    char output;
    int stat = xtb::GFN0_calculation(&natoms, attyp, &charge, NULL, coord, &opt, &output, &energy, grad);
#else
    throw("XTB is not included, sorry for that");
#endif
    return energy;
}

double XTBInterface::GFN0Energy(const int* attyp, const double* coord, const int natoms, const double charge)
{
    double energy = 0;
#ifdef USE_XTB
    xtb::PEEQ_options const opt = (xtb::PEEQ_options){ 0, 1, 2.0, 300, false, false, "none" };

    double dipole[3];
    double q[natoms];
    double qp[6 * natoms];
    double wbo[natoms * natoms];
    char output;
    double* grad = 0;
    int stat = xtb::GFN0_calculation(&natoms, attyp, &charge, NULL, coord, &opt, &output, &energy, grad);
#else
    throw("XTB is not included, sorry for that");
#endif
    return energy;
}
