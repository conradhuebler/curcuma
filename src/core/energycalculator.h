/*
 * < General Energy and Gradient Calculator >
 * Copyright (C) 2022 - 2024 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#ifdef USE_TBLITE
#include "src/core/tbliteinterface.h"
#endif

#ifdef USE_XTB
#include "src/core/xtbinterface.h"
#endif

#ifdef USE_D3
#include "src/core/dftd3interface.h"
#endif

#ifdef USE_D4
#include "src/core/dftd4interface.h"
#endif

#include "src/core/eht.h"
#include "src/core/eigen_uff.h"
#include "src/core/forcefield.h"
#include "src/core/qmdff.h"

#include <functional>

#include "json.hpp"
using json = nlohmann::json;

static json EnergyCalculatorJson{
    { "param_file", "none" }
};

class EnergyCalculator {
public:
    EnergyCalculator(const std::string& method, const json& controller);
    ~EnergyCalculator();

    void setMolecule(const Molecule& molecule);

    void updateGeometry(const double* coord);
    void updateGeometry(const std::vector<double>& geometry);

    void updateGeometry(const Matrix& geometry);
    void updateGeometry(const Eigen::VectorXd& geometry);

    void getGradient(double* coord);

    Matrix getGradient() const { return m_gradient; }

    Matrix Gradient() const;

    double CalculateEnergy(bool gradient = false, bool verbose = false);

    bool HasNan() const { return m_containsNaN; }

    bool Error() const { return m_error; }
#ifdef USE_TBLITE
    TBLiteInterface* getTBLiterInterface() const
    {
        return m_tblite;
    }
#endif

#ifdef USE_XTB
    XTBInterface* getXTBInterface() const
    {
        return m_xtb;
    }
#endif

#ifdef USE_D3
    DFTD3Interface* getD3Interface() const
    {
        return m_d3;
    }
#endif

#ifdef USE_D4
    DFTD4Interface* getD4Interface() const
    {
        return m_d4;
    }
#endif

    QMDFF* getQMDFFInterface() const
    {
        return m_qmdff;
    }

    inline void setParameter(const json& parameter)
    {
        m_parameter = parameter;
        // if (m_qmdff != NULL)
        //     m_qmdff->setParameter(parameter);
    }

    inline json Parameter() const { return m_parameter; }
    Vector Energies() const { return m_orbital_energies; }
    Vector OrbitalOccuptations() const { return m_orbital_occupation; }

    std::vector<double> Charges() const;
    Position Dipole() const;

    std::vector<std::vector<double>> BondOrders() const;
    int NumElectrons() const { return m_num_electrons; }

private:
    void InitialiseUFF();
    void CalculateUFF(bool gradient, bool verbose = false);

    void InitialiseTBlite();
    void CalculateTBlite(bool gradient, bool verbose = false);

    void InitialiseXTB();
    void CalculateXTB(bool gradient, bool verbose = false);

    void InitialiseD4();
    void CalculateD4(bool gradient, bool verbose = false);

    void InitialiseD3();
    void CalculateD3(bool gradient, bool verbose = false);

    void InitialiseQMDFF();
    void CalculateQMDFF(bool gradient, bool verbose = false);

    void InitialiseFF();
    void CalculateFF(bool gradient, bool verbose = false);

    json m_controller;

#ifdef USE_TBLITE
    TBLiteInterface* m_tblite = NULL;
#endif

#ifdef USE_XTB
    XTBInterface* m_xtb = NULL;
#endif
#ifdef USE_D3
    DFTD3Interface* m_d3 = NULL;
#endif
#ifdef USE_D4
    DFTD4Interface* m_d4 = NULL;
#endif

    eigenUFF* m_uff = NULL;
    QMDFF* m_qmdff = NULL;
    ForceField* m_forcefield = NULL;
    EHT* m_eht = NULL;
    StringList m_uff_methods = { "fuff" };
    StringList m_ff_methods = { "uff", "uff-d3", "qmdff" };
    StringList m_qmdff_method = { "fqmdff" };
    StringList m_tblite_methods = { "ipea1", "gfn1", "gfn2" };
    StringList m_xtb_methods = { "gfnff", "xtb-gfn1", "xtb-gfn2" };
    StringList m_d3_methods = { "d3" };
    StringList m_d4_methods = { "d4" };

    std::function<void(bool, bool)> m_ecengine;
    std::function<std::vector<double>()> m_charges;
    std::function<Position()> m_dipole;
    std::function<std::vector<std::vector<double>>()> m_bonds;
    json m_parameter;
    std::string m_method, m_param_file;
    Matrix m_geometry, m_gradient, m_molecular_orbitals;
    Vector m_orbital_energies, m_orbital_occupation;

    double m_energy;
    double *m_coord, *m_grad;

    int m_atoms, m_num_electrons = 0;
    int m_gfn = 2;
    int* m_atom_type;

    bool m_initialised = false;
    bool m_containsNaN = false;
    bool m_error = false;
    bool m_writeparam = false;
};
