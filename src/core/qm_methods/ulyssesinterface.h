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

#pragma once

#include "interface/abstract_interface.h"

#include <string>
#include <vector>

static json UlyssesSettings{
    { "Tele", 300 },
    { "solvent", "none" },
    { "method", "GFN2" },
    { "SCFmaxiter", 100 },
    { "verbose", 0 },
    { "mult", 1 }

};
class UlyssesObject;

class UlyssesInterface : public QMInterface {
public:
    UlyssesInterface(const json& ulyssessettings = UlyssesSettings);
    ~UlyssesInterface();

    virtual bool InitialiseMolecule() override;
    virtual bool UpdateMolecule(const Geometry& geometry) override;

    double Calculation(bool gradient = false, bool verbose = false) override;
    virtual Vector Charges() const override;
    virtual Vector OrbitalEnergies() const override;

private:
    Mol m_mol;
    UlyssesObject* m_ulysses = nullptr;

    double m_Tele = 300;
    double m_SCFmaxiter = 100;
    bool m_verbose = 0;
    int m_mult = 1;

    StringList m_solvents = { "acetone", "acetone", "acetonitrile", "aniline", "benzaldehyde", "benzene", "dichloromethane", "chloroform", "carbon disulfide", "dioxane", "dmf", "dmso", "ethanol", "diethyl ether", "ethyl acetate", "furane", "hexadecane", "hexane", "methanol", "nitromethane", "octanol", "phenol", "thf", "toluene", "water", "octanol wet" };
    std::string m_solvent = "none";

    json m_ulyssessettings;
};
