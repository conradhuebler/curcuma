/*
 * < C++ XTB and tblite Interface >
 * Copyright (C) 2020 - 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "interface/abstract_interface.h"

#include "src/core/global.h"
#include "src/tools/general.h"

#include <iostream>
#include <math.h>
#include <stdio.h>

#include "dftd4interface.h"

DFTD4Interface::DFTD4Interface(const ConfigManager& config)
    : m_config(config)
{
    // Claude Generated 2025: ConfigManager migration - Phase 3B
    m_par.a1 = m_config.get<double>("a1", 0.40085597);
    m_par.a2 = m_config.get<double>("a2", 5.02928789);
    m_par.alp = m_config.get<double>("alpha", 16.0);

    m_par.s6 = m_config.get<double>("s6", 1.00);
    m_par.s8 = m_config.get<double>("s8", 1.20065498);
    m_par.s10 = m_config.get<double>("s10", 0.0);

    m_par.s9 = m_config.get<double>("s9", 1.0);

    std::string functional = m_config.get<std::string>("functional", "pbe0");
    bool three_body = m_config.get<bool>("three_body", true);
    dftd4::d4par(functional, m_par, three_body);
    PrintParameter();
}

DFTD4Interface::DFTD4Interface()
    : m_config("dftd4", json{})
{
}

DFTD4Interface::~DFTD4Interface()
{
    m_mol.FreeMemory();
}

void DFTD4Interface::PrintParameter() const
{
    // std::cout << m_par.s8 << " " << m_par.s8 << " " << m_par.s9 << " " << m_par.s10 << " " << m_par.a1 << " " << m_par.a2 << " " << m_par.alp << std::endl;
}

// Claude Generated 2025: ConfigManager migration - Phase 2.2
// ConfigManager-based parameter update (modern version)
void DFTD4Interface::UpdateParameters(const ConfigManager& config)
{
    // Update functional and three-body settings first
    std::string functional = config.get<std::string>("functional", "pbe0");
    bool three_body = config.get<bool>("three_body", true);
    dftd4::d4par(functional, m_par, three_body);

    // Update damping parameters, keeping existing values as defaults
    m_par.a1 = config.get<double>("a1", m_par.a1);
    m_par.a2 = config.get<double>("a2", m_par.a2);
    m_par.alp = config.get<double>("alpha", m_par.alp);

    // Update scaling parameters
    m_par.s6 = config.get<double>("s6", m_par.s6);
    m_par.s8 = config.get<double>("s8", m_par.s8);
    m_par.s10 = config.get<double>("s10", m_par.s10);
    m_par.s9 = config.get<double>("s9", m_par.s9);

    PrintParameter();
}

// Backward compatibility: JSON version delegates to ConfigManager
void DFTD4Interface::UpdateParameters(const json& controller)
{
    ConfigManager param_config("dftd4", controller);
    UpdateParameters(param_config);
}

bool DFTD4Interface::InitialiseMolecule(const Mol& mol, double factor)
{
    m_mol.GetMemory(mol.m_number_atoms);
    for (int i = 0; i < mol.m_number_atoms; ++i) {
        int element = mol.m_atoms[i];
        m_mol.UpdateAtom(i, mol.m_geometry(i, 0) * factor, mol.m_geometry(i, 1) * factor, mol.m_geometry(i, 2), element);
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

double DFTD4Interface::Calculation(bool gradient)
{
    double energy = 0;
    dftd4::TCutoff cutoff;
    dftd4::TD4Model d4;
    // exit(0);
    dftd4::get_dispersion(m_mol, m_charge, d4, m_par, cutoff, energy, m_gradient.data());
    return energy;
}

void DFTD4Interface::clear()
{
}
