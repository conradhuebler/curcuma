/*
 * < C++ XTB and tblite Interface >
 * Copyright (C) 2020 - 2024 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "src/tools/general.h"
#ifdef USE_D3
#include "s-dftd3.h"
#endif

#include "src/core/molecule.h"

static json DFTD3Settings{
    { "d_s6", 0 },
    { "d_s8", 0 },
    { "d_s9", 0 },
    { "d_a1", 0 },
    { "d_a2", 0 },
    { "d_alp", 0 },
    { "d_bet", 8 },
    { "d_func", "pbe0" },
    { "d_atm", false },
    { "d_damping", "bj" }
};

class DFTD3Interface {
public:
    DFTD3Interface(const json& controller);
    DFTD3Interface();

    ~DFTD3Interface();

    bool InitialiseMolecule(const std::vector<int>& atomtypes);

    double DFTD3Calculation(double* grad = 0);
    void UpdateParameters(const json& controller);
    void UpdateParametersD3(const json& controller);

    void clear();

    void UpdateGeometry(const double* coord);

    void PrintParameter() const;
    void UpdateAtom(int index, double x, double y, double z);

    inline double ParameterA1() const { return m_d3_a1; }
    inline double ParameterA2() const { return m_d3_a2; }
    inline double ParameterAlp() const { return m_d3_alp; }
    inline double ParameterS6() const { return m_d3_s6; }
    inline double ParameterS8() const { return m_d3_s8; }
    inline double ParameterS9() const { return m_d3_s9; }

private:
    void CreateParameter();

    std::string m_damping;
    std::string m_functional;

    int m_charge = 0;
    double* m_coord;
    int* m_attyp;

    double m_d3_a1 = 1;
    double m_d3_a2 = 1;
    double m_d3_alp = 1;

    double m_d3_s6 = 1;
    double m_d3_s8 = 1;
    double m_d3_s9 = 1;

    double m_bet = 8;
    bool m_atm = false;

#ifdef USE_D3
    dftd3_error m_error;
    dftd3_structure m_mol;
    dftd3_model m_disp;
    dftd3_param m_param;
#endif
};
