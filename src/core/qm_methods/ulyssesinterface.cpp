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

#include <string>
#include <vector>

#include "interface/abstract_interface.h"
#ifdef USE_ULYSSES
#include "interface/ulysses.h"
#endif

#include "ulyssesinterface.h"

UlyssesInterface::UlyssesInterface(const json& ulyssessettings)
{
    m_ulyssessettings = MergeJson(UlyssesSettings, ulyssessettings);
    m_Tele = m_ulyssessettings["Tele"];
    m_SCFmaxiter = m_ulyssessettings["SCFmaxiter"];
    m_solvent = m_ulyssessettings["solvent"];
    m_method = m_ulyssessettings["method"];
    m_mult = m_ulyssessettings["mult"];
    m_verbose = m_ulyssessettings["verbose"];
    if (std::find(m_solvents.begin(), m_solvents.end(), m_solvent) == m_solvents.end()) {
        std::cout << "Solvent " << m_solvent << " is not supported by Ulysses" << std::endl;
        m_solvent = "none";
    }
#ifdef USE_ULYSSES
    m_ulysses = new UlyssesObject();
#endif
}

UlyssesInterface::~UlyssesInterface()
{
#ifdef USE_ULYSSES
    delete m_ulysses;
#endif
}

bool UlyssesInterface::InitialiseMolecule()
{
#ifdef USE_ULYSSES
    m_ulysses->setMethod(m_method);
    m_ulysses->setMolecule(m_geometry, m_atoms, m_charge, m_mult, "C1");
    if (m_verbose > 0)
        std::cout << "Initialising Ulysses with method " << m_method << " and SCFmaxiter " << m_SCFmaxiter << std::endl;

    return true;
#else
    return false;
#endif
}

bool UlyssesInterface::UpdateMolecule(const Geometry& geometry)
{
#ifdef USE_ULYSSES
    m_geometry = geometry;
    m_ulysses->UpdateGeometry(geometry);
    return true;
#else
    return false;
#endif
}

double UlyssesInterface::Calculation(bool gradient, bool verbose)
{
#ifdef USE_ULYSSES
    m_ulysses->setTele(m_Tele);
    m_ulysses->setMaxIter(m_SCFmaxiter);
    m_ulysses->setSolvent(m_solvent);
    m_ulysses->Calculate(gradient, verbose);
    if (gradient) {
        m_gradient = m_ulysses->Gradient();
    }
    return m_ulysses->Energy();
#else
    return 0;
#endif
}

Vector UlyssesInterface::Charges() const
{
#ifdef USE_ULYSSES
    return m_ulysses->Charges();
#else
    return Vector{};
#endif
}

Vector UlyssesInterface::OrbitalEnergies() const
{
#ifdef USE_ULYSSES
    return m_ulysses->OrbitalEnergies();
#else
    return Vector{};
#endif
}