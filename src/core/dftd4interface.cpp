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

#include "dftd_damping.h"
#include "dftd_dispersion.h"
#include "dftd_econv.h"
#include "dftd_geometry.h"

#include "src/core/global.h"
#include "src/tools/general.h"

#include "src/core/molecule.h"

#include <iostream>
#include <math.h>
#include <stdio.h>

#include "dftd4interface.h"

DFTD4Interface::DFTD4Interface(const json& controller)
{
    json parameter = MergeJson(DFTD4Settings, controller);
    m_par.a1 = parameter["d4_a1"];
    m_par.a2 = parameter["d4_a2"];
    m_par.alp = parameter["d4_alp"];

    m_par.s6 = parameter["d4_s6"];
    m_par.s8 = parameter["d4_s8"];
    m_par.s10 = parameter["d4_s10"];

    m_par.s9 = parameter["d4_s9"];

    dftd4::d4par(parameter["d4_func"], m_par, parameter["d4_atm"]);
    PrintParameter();
}

DFTD4Interface::~DFTD4Interface()
{
    m_mol.FreeMemory();
}

void DFTD4Interface::PrintParameter() const
{
    std::cout << m_par.s8 << " " << m_par.s8 << " " << m_par.s9 << " " << m_par.s10 << " " << m_par.a1 << " " << m_par.a2 << " " << m_par.alp << std::endl;
}

void DFTD4Interface::UpdateParameters(const json& controller)
{
    json parameter = MergeJson(DFTD4Settings, controller);

    dftd4::d4par(parameter["d4_func"], m_par, parameter["d4_atm"]);
    m_par.a1 = parameter["d4_a1"];
    m_par.a2 = parameter["d4_a2"];
    m_par.alp = parameter["d4_alp"];

    m_par.s6 = parameter["d4_s6"];
    m_par.s8 = parameter["d4_s8"];
    m_par.s10 = parameter["d4_s10"];

    m_par.s9 = parameter["d4_s9"];
    PrintParameter();
}

bool DFTD4Interface::InitialiseMolecule(const Molecule& molecule, double factor)
{
    m_mol.GetMemory(molecule.AtomCount());
    for (int i = 0; i < molecule.AtomCount(); ++i) {
        std::pair<int, Position> atom = molecule.Atom(i);
        m_mol.UpdateAtom(i, atom.second(0) * factor, atom.second(1) * factor, atom.second(2) * factor, atom.first);
    }
    return true;
}

bool DFTD4Interface::InitialiseMolecule(const std::vector<int>& atomtype)
{
    m_mol.GetMemory(atomtype.size());
    for (int i = 0; i < atomtype.size(); ++i) {
        m_mol.UpdateAtom(i, 0, 0, 0, atomtype[i]);
    }
    return true;
}

double DFTD4Interface::DFTD4Calculation(double* grad)
{
    double energy = 0;
    dftd4::TCutoff cutoff;

    int info;
    /*
    int charge;
    std::cout << std::endl << std::endl;
    for (int i = 0; i < atoms.size(); ++i) {
        std::pair<int, Position> atom = m_molecule.Atom(i);
        m_mol.UpdateAtom(i, atom.second(0) / au, atom.second(1) / au, atom.second(2) / au);
    }*/
    // m_mol.Print();
    info = dftd4::get_dispersion(m_mol, m_charge, m_par, cutoff, energy, grad);
    // std::cout << energy << " " << info << std::endl;
    return energy;
}

void DFTD4Interface::clear()
{
}
