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

#include "s-dftd3.h"

#include "src/core/global.h"
#include "src/tools/general.h"

#include "interface/abstract_interface.h"

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <cstring>

#include "dftd3interface.h"

DFTD3Interface::DFTD3Interface(const ConfigManager& config)
    : m_config(config)
{
    // Claude Generated 2025: ConfigManager migration - Phase 3B
    m_d3_a1 = m_config.get<double>("a1", 0.0);
    m_d3_a2 = m_config.get<double>("a2", 0.0);
    m_d3_alp = m_config.get<double>("alpha", 14.0);

    m_d3_s6 = m_config.get<double>("s6", 1.0);
    m_d3_s8 = m_config.get<double>("s8", 0.0);
    m_d3_s9 = m_config.get<double>("s9", 0.0);
    m_bet = m_config.get<double>("beta", 8.0);

    m_atm = m_config.get<bool>("three_body", false);
    m_damping = m_config.get<std::string>("damping", "bj");
    m_functional = m_config.get<std::string>("functional", "pbe0");
    CreateParameter();
#ifdef USE_D3
    m_error = dftd3_new_error();
#endif
}

DFTD3Interface::DFTD3Interface()
    : m_config("dftd3", json{})
{
#ifdef USE_D3
    m_error = dftd3_new_error();
#endif
}

DFTD3Interface::~DFTD3Interface()
{
#ifdef USE_D3
    dftd3_delete_model(&m_disp);
    dftd3_delete_structure(&m_mol);
    dftd3_delete_error(&m_error);
    dftd3_delete_param(&m_param);

    delete m_coord;
    delete m_attyp;
#endif
}

void DFTD3Interface::PrintParameter() const
{
    std::cout << m_d3_s6 << " " << m_d3_s8 << " " << m_d3_s9 << " " << m_d3_a1 << " " << m_d3_a2 << " " << m_d3_alp << std::endl;
}

void DFTD3Interface::CreateParameter()
{
#ifdef USE_D3
    if (m_d3_a1 > 1e-8 || m_d3_a2 > 1e-8 || m_d3_s6 > 1e-8 || m_d3_s8 > 1e-8 || m_d3_s9 > 1e-8) {
        if (m_damping.compare("bj")) {
            m_param = dftd3_new_rational_damping(m_error, m_d3_s6, m_d3_s8, m_d3_s9, m_d3_a1, m_d3_a2, m_d3_alp);
        } else if (m_damping.compare("zero")) {
            m_param = dftd3_new_zero_damping(m_error, m_d3_s6, m_d3_s8, m_d3_s9, m_d3_a1, m_d3_a2, m_d3_alp);
        } else if (m_damping.compare("bjm")) {
            m_param = dftd3_new_mrational_damping(m_error, m_d3_s6, m_d3_s8, m_d3_s9, m_d3_a1, m_d3_a2, m_d3_alp);
        } else if (m_damping.compare("zerom")) {
            m_param = dftd3_new_mzero_damping(m_error, m_d3_s6, m_d3_s8, m_d3_s9, m_d3_a1, m_d3_a2, m_d3_alp, m_bet);
        } else if (m_damping.compare("op")) {
            m_param = dftd3_new_optimizedpower_damping(m_error, m_d3_s6, m_d3_s8, m_d3_s9, m_d3_a1, m_d3_a2, m_d3_alp, m_bet);
        }
    } else {
        char* cstr = new char[m_functional.length() + 1];
        strcpy(cstr, m_functional.c_str());
        if (m_damping.compare("bj") == 0) {
            m_param = dftd3_load_rational_damping(m_error, cstr, m_atm);
        } else if (m_damping.compare("zero") == 0) {
            m_param = dftd3_load_zero_damping(m_error, cstr, m_atm);
        } else if (m_damping.compare("bjm") == 0) {
            m_param = dftd3_load_mrational_damping(m_error, cstr, m_atm);
        } else if (m_damping.compare("zerom") == 0) {
            m_param = dftd3_load_mzero_damping(m_error, cstr, m_atm);
        } else if (m_damping.compare("op") == 0) {
            m_param = dftd3_load_optimizedpower_damping(m_error, cstr, m_atm);
        }
        delete[] cstr;
    }
#endif
}

// Claude Generated 2025: ConfigManager migration - Phase 2.1
// ConfigManager-based parameter update (modern version)
void DFTD3Interface::UpdateParameters(const ConfigManager& config)
{
    // Update parameters, keeping existing values as defaults
    m_d3_a1 = config.get<double>("a1", m_d3_a1);
    m_d3_a2 = config.get<double>("a2", m_d3_a2);
    m_d3_alp = config.get<double>("alpha", m_d3_alp);

    m_d3_s6 = config.get<double>("s6", m_d3_s6);
    m_d3_s8 = config.get<double>("s8", m_d3_s8);
    m_d3_s9 = config.get<double>("s9", m_d3_s9);

    m_bet = config.get<double>("beta", m_bet);
    m_atm = config.get<bool>("three_body", m_atm);

    // Reinitialize DFT-D3 parameter structure with new values
    CreateParameter();
}

// Backward compatibility: JSON version delegates to ConfigManager
void DFTD3Interface::UpdateParameters(const json& controller)
{
    ConfigManager param_config("dftd3", controller);
    UpdateParameters(param_config);
}

// Claude Generated 2025: ConfigManager migration - Phase 2.1
// D3-specific parameter update (modern version)
void DFTD3Interface::UpdateParametersD3(const ConfigManager& config)
{
    // Same implementation as UpdateParameters for D3
    UpdateParameters(config);
}

// Backward compatibility: JSON version delegates to ConfigManager
void DFTD3Interface::UpdateParametersD3(const json& controller)
{
    ConfigManager param_config("dftd3", controller);
    UpdateParametersD3(param_config);
}

bool DFTD3Interface::InitialiseMolecule(const std::vector<int>& atomtypes)
{
    m_natoms = atomtypes.size();
    m_attyp = new int[atomtypes.size()];
    m_coord = new double[3 * atomtypes.size()];
    for (int i = 0; i < atomtypes.size(); ++i) {
        m_coord[3 * i + 0] = (3 * i + 0) / au;
        m_coord[3 * i + 1] = (3 * i + 1) / au;
        m_coord[3 * i + 2] = (3 * i + 2) / au;
        m_attyp[i] = atomtypes[i];
    }
#ifdef USE_D3
    m_mol = dftd3_new_structure(m_error, atomtypes.size(), m_attyp, m_coord, NULL, NULL);
    m_disp = dftd3_new_d3_model(m_error, m_mol);
#endif

    return true;
}

void DFTD3Interface::UpdateAtom(int index, double x, double y, double z)
{
    m_coord[3 * index + 0] = x / au;
    m_coord[3 * index + 1] = y / au;
    m_coord[3 * index + 2] = z / au;
}

void DFTD3Interface::UpdateGeometry(const double* coord)
{
    // Copy all coordinates (assumed to be in Bohr already)
    if (m_coord && m_natoms > 0) {
        std::memcpy(m_coord, coord, 3 * m_natoms * sizeof(double));
    }
}

double DFTD3Interface::Calculation(bool gradient)
{
    double energy = 0;
    double sigma[9];

#ifdef USE_D3
    dftd3_update_structure(m_error, m_mol, m_coord, NULL);
    if (gradient)
        dftd3_get_dispersion(m_error, m_mol, m_disp, m_param, &energy, m_gradient.data(), sigma);
#endif

    return energy;
}

void DFTD3Interface::clear()
{
#ifdef USE_D3
    dftd3_delete_model(&m_disp);
    dftd3_delete_structure(&m_mol);
    dftd3_delete_error(&m_error);
#endif
}
