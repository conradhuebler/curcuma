/*
 * <Extended Hückel Theory Parameters Module Implementation>
 * Copyright (C) 2023 Conrad Hübler <Conrad.Huebler@gmx.net>
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
 */

#include "eht_parameters.h"
#include <iostream>
#include <stdexcept>

namespace EHTParams {

// Static member initialization
std::map<int, ElementParameters> ParameterDatabase::s_parameters;
bool ParameterDatabase::s_initialized = false;

std::map<std::pair<int, STO::OrbitalType>, STO6GParameters> STO6GDatabase::s_sto6g_parameters;
bool STO6GDatabase::s_initialized = false;

void ParameterDatabase::initializeDatabase()
{
    if (s_initialized)
        return;

    // Hydrogen (Z=1) - Parameters from Hoffmann, R. J. Chem. Phys. 39, 1397 (1963)
    ElementParameters hydrogen;
    hydrogen[STO::S] = OrbitalParameters(1.20, -13.6, 1); // 1s: ζ=1.20, VSIP=-13.6 eV
    s_parameters[1] = hydrogen;

    // Carbon (Z=6) - Parameters optimized for organic molecules
    ElementParameters carbon;
    carbon[STO::S] = OrbitalParameters(1.625, -19.43, 2); // 2s: ζ=1.625, VSIP=-19.43 eV
    carbon[STO::PX] = OrbitalParameters(1.625, -10.64, 1); // 2px: ζ=1.625, VSIP=-10.64 eV
    carbon[STO::PY] = OrbitalParameters(1.625, -10.64, 1); // 2py: ζ=1.625, VSIP=-10.64 eV
    carbon[STO::PZ] = OrbitalParameters(1.625, -10.64, 0); // 2pz: ζ=1.625, VSIP=-10.64 eV
    s_parameters[6] = carbon;

    // Nitrogen (Z=7) - Parameters for sp3 hybridization
    ElementParameters nitrogen;
    nitrogen[STO::S] = OrbitalParameters(1.950, -25.56, 2); // 2s: ζ=1.950, VSIP=-25.56 eV
    nitrogen[STO::PX] = OrbitalParameters(1.950, -13.18, 1); // 2px: ζ=1.950, VSIP=-13.18 eV
    nitrogen[STO::PY] = OrbitalParameters(1.950, -13.18, 1); // 2py: ζ=1.950, VSIP=-13.18 eV
    nitrogen[STO::PZ] = OrbitalParameters(1.950, -13.18, 1); // 2pz: ζ=1.950, VSIP=-13.18 eV
    s_parameters[7] = nitrogen;

    // Oxygen (Z=8) - Parameters for divalent oxygen
    ElementParameters oxygen;
    oxygen[STO::S] = OrbitalParameters(2.275, -32.38, 2); // 2s: ζ=2.275, VSIP=-32.38 eV
    oxygen[STO::PX] = OrbitalParameters(2.275, -15.85, 2); // 2px: ζ=2.275, VSIP=-15.85 eV
    oxygen[STO::PY] = OrbitalParameters(2.275, -15.85, 1); // 2py: ζ=2.275, VSIP=-15.85 eV
    oxygen[STO::PZ] = OrbitalParameters(2.275, -15.85, 1); // 2pz: ζ=2.275, VSIP=-15.85 eV
    s_parameters[8] = oxygen;

    // Fluorine (Z=9) - Parameters for halogen chemistry
    ElementParameters fluorine;
    fluorine[STO::S] = OrbitalParameters(2.425, -40.17, 2); // 2s: ζ=2.425, VSIP=-40.17 eV
    fluorine[STO::PX] = OrbitalParameters(2.425, -18.65, 2); // 2px: ζ=2.425, VSIP=-18.65 eV
    fluorine[STO::PY] = OrbitalParameters(2.425, -18.65, 2); // 2py: ζ=2.425, VSIP=-18.65 eV
    fluorine[STO::PZ] = OrbitalParameters(2.425, -18.65, 1); // 2pz: ζ=2.425, VSIP=-18.65 eV
    s_parameters[9] = fluorine;

    s_initialized = true;
}

OrbitalParameters ParameterDatabase::getOrbitalParameters(int atomic_number, STO::OrbitalType orbital_type)
{
    initializeDatabase();

    auto element_it = s_parameters.find(atomic_number);
    if (element_it == s_parameters.end()) {
        throw std::runtime_error("EHT parameters not available for element Z=" + std::to_string(atomic_number));
    }

    auto orbital_it = element_it->second.find(orbital_type);
    if (orbital_it == element_it->second.end()) {
        throw std::runtime_error("EHT parameters not available for orbital type " + std::to_string(static_cast<int>(orbital_type)) + " of element Z=" + std::to_string(atomic_number));
    }

    return orbital_it->second;
}

ElementParameters ParameterDatabase::getElementParameters(int atomic_number)
{
    initializeDatabase();

    auto it = s_parameters.find(atomic_number);
    if (it == s_parameters.end()) {
        throw std::runtime_error("EHT parameters not available for element Z=" + std::to_string(atomic_number));
    }

    return it->second;
}

bool ParameterDatabase::hasElement(int atomic_number)
{
    initializeDatabase();
    return s_parameters.find(atomic_number) != s_parameters.end();
}

int ParameterDatabase::getValenceElectrons(int atomic_number)
{
    initializeDatabase();

    auto it = s_parameters.find(atomic_number);
    if (it == s_parameters.end()) {
        throw std::runtime_error("EHT parameters not available for element Z=" + std::to_string(atomic_number));
    }

    int total_electrons = 0;
    for (const auto& orbital_pair : it->second) {
        total_electrons += orbital_pair.second.n_electrons;
    }

    return total_electrons;
}

std::string ParameterDatabase::getElementSymbol(int atomic_number)
{
    switch (atomic_number) {
    case 1:
        return "H";
    case 6:
        return "C";
    case 7:
        return "N";
    case 8:
        return "O";
    case 9:
        return "F";
    default:
        return "X";
    }
}

void STO6GDatabase::initializeDatabase()
{
    if (s_initialized)
        return;

    // Hydrogen 1s STO-6G parameters
    s_sto6g_parameters[{ 1, STO::S }] = STO6GParameters(
        { 3.552322122E+01, 6.513143725E+00, 1.822142904E+00, 6.259552659E-01, 2.430767471E-01, 1.001124280E-01 },
        { 9.163596281E-03, 4.936149294E-02, 1.685383049E-01, 3.705627997E-01, 4.164915298E-01, 1.303340841E-01 });

    // Carbon 2s STO-6G parameters
    s_sto6g_parameters[{ 6, STO::S }] = STO6GParameters(
        { 3.049723950E+01, 6.036199601E+00, 1.876046337E+00, 7.217826470E-01, 3.134706954E-01, 1.436865550E-01 },
        { -1.325278809E-02, -4.699171014E-02, -3.378537151E-02, 2.502417861E-01, 5.951172526E-01, 2.407061763E-01 });

    // Carbon 2p STO-6G parameters (same for px, py, pz)
    STO6GParameters carbon_p(
        { 3.049723950E+01, 6.036199601E+00, 1.876046337E+00, 7.217826470E-01, 3.134706954E-01, 1.436865550E-01 },
        { 3.759696623E-03, 3.767936984E-02, 1.738967435E-01, 4.180364347E-01, 4.258595477E-01, 1.017082955E-01 });
    s_sto6g_parameters[{ 6, STO::PX }] = carbon_p;
    s_sto6g_parameters[{ 6, STO::PY }] = carbon_p;
    s_sto6g_parameters[{ 6, STO::PZ }] = carbon_p;

    // Nitrogen 2s STO-6G parameters
    s_sto6g_parameters[{ 7, STO::S }] = STO6GParameters(
        { 3.919880787E+01, 7.758467071E+00, 2.411325783E+00, 9.277239437E-01, 4.029111410E-01, 1.846836552E-01 },
        { -1.325278809E-02, -4.699171014E-02, -3.378537151E-02, 2.502417861E-01, 5.951172526E-01, 2.407061763E-01 });

    // Nitrogen 2p STO-6G parameters
    STO6GParameters nitrogen_p(
        { 3.919880787E+01, 7.758467071E+00, 2.411325783E+00, 9.277239437E-01, 4.029111410E-01, 1.846836552E-01 },
        { 3.759696623E-03, 3.767936984E-02, 1.738967435E-01, 4.180364347E-01, 4.258595477E-01, 1.017082955E-01 });
    s_sto6g_parameters[{ 7, STO::PX }] = nitrogen_p;
    s_sto6g_parameters[{ 7, STO::PY }] = nitrogen_p;
    s_sto6g_parameters[{ 7, STO::PZ }] = nitrogen_p;

    // Oxygen 2s STO-6G parameters
    s_sto6g_parameters[{ 8, STO::S }] = STO6GParameters(
        { 5.218776196E+01, 1.032932006E+01, 3.210344977E+00, 1.235135428E+00, 5.364201581E-01, 2.458806060E-01 },
        { -1.325278809E-02, -4.699171014E-02, -3.378537151E-02, 2.502417861E-01, 5.951172526E-01, 2.407061763E-01 });

    // Oxygen 2p STO-6G parameters
    STO6GParameters oxygen_p(
        { 5.218776196E+01, 1.032932006E+01, 3.210344977E+00, 1.235135428E+00, 5.364201581E-01, 2.458806060E-01 },
        { 3.759696623E-03, 3.767936984E-02, 1.738967435E-01, 4.180364347E-01, 4.258595477E-01, 1.017082955E-01 });
    s_sto6g_parameters[{ 8, STO::PX }] = oxygen_p;
    s_sto6g_parameters[{ 8, STO::PY }] = oxygen_p;
    s_sto6g_parameters[{ 8, STO::PZ }] = oxygen_p;

    s_initialized = true;
}

STO6GParameters STO6GDatabase::getParameters(int atomic_number, STO::OrbitalType orbital_type)
{
    initializeDatabase();

    auto key = std::make_pair(atomic_number, orbital_type);
    auto it = s_sto6g_parameters.find(key);

    if (it == s_sto6g_parameters.end()) {
        throw std::runtime_error("STO-6G parameters not available for element Z=" + std::to_string(atomic_number) + ", orbital type " + std::to_string(static_cast<int>(orbital_type)));
    }

    return it->second;
}

bool STO6GDatabase::hasParameters(int atomic_number, STO::OrbitalType orbital_type)
{
    initializeDatabase();
    auto key = std::make_pair(atomic_number, orbital_type);
    return s_sto6g_parameters.find(key) != s_sto6g_parameters.end();
}

} // namespace EHTParams