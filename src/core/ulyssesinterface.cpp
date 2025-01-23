/*
 * < C++ Ulysses Interface >
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "interface/abstract_interface.h"

#include <string>
#include <vector>

#include "src/core/interface/ulysses.h"

#include "ulyssesinterface.h"

UlyssesInterface::UlyssesInterface(const json& ulyssessettings)
{
    m_ulyssessettings = MergeJson(UlyssesSettings, ulyssessettings);
    m_Tele = m_ulyssessettings["Tele"];
    m_SCFmaxiter = m_ulyssessettings["SCFmaxiter"];
    m_solvent = m_ulyssessettings["ulysses_solvent"];
    m_method = m_ulyssessettings["method"];
    m_mult = m_ulyssessettings["mult"];
    m_ulysses = new UlyssesObject();

}

UlyssesInterface::~UlyssesInterface()
{
    delete m_ulysses;
}

bool UlyssesInterface::InitialiseMolecule()
{
    m_ulysses->setMethod(m_method);
    m_ulysses->setMolecule(m_geometry, m_atoms, m_charge, m_mult, "C1");
    std::cout << "Initialising Ulysses with method " <<m_method << " and SCFmaxiter " << m_SCFmaxiter << std::endl;

    return true;
}

bool UlyssesInterface::UpdateMolecule(const Geometry& geometry)
{
    m_geometry = geometry;
    m_ulysses->UpdateGeometry(geometry);
    return true;
}

double UlyssesInterface::Calculation(bool gradient, bool verbose)
{
    m_ulysses->setTele(m_Tele);
    m_ulysses->setMaxIter(m_SCFmaxiter);
    m_ulysses->Calculate(gradient, verbose);
    if (gradient) {
        m_gradient = m_ulysses->Gradient();
        // std::cout << m_gradient << std::endl;
    }
    return m_ulysses->Energy();
}