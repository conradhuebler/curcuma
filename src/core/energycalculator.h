/*
 * < General Energy and Gradient Calculator >
 * Copyright (C) 2022 - 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
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

// #include "src/tools/general.h"

#ifdef USE_TBLITE
#include "src/core/qm_methods/tbliteinterface.h"
#endif

#ifdef USE_XTB
#include "src/core/qm_methods/xtbinterface.h"
#endif

#ifdef USE_D3
#include "src/core/qm_methods/dftd3interface.h"
#endif

#ifdef USE_D4
#include "src/core/qqm_methods/dftd4interface.h"
#endif

#include "src/core/qm_methods/eht.h"
#include "src/core/qm_methods/gfnff.h"
// #include "src/core/eigen_uff.h"
#include "src/core/forcefield.h"
// #include "src/core/qmdff.h"

#ifdef USE_ULYSSES
#include "src/core/qm_methods/ulyssesinterface.h"
#endif

#include <functional>

#include "json.hpp"
using json = nlohmann::json;

static json EnergyCalculatorJson{
    { "param_file", "none" },
    { "multi", 1 },
    { "method", "uff" },
    { "SCFmaxiter", 100 },
    { "Tele", 300 },
    { "solvent", "none" },
    { "verbose", 0 }

};

class EnergyCalculator {
public:
    EnergyCalculator(const std::string& method, const json& controller);
    // Claude Generated: Constructor with basename for parameter caching
    EnergyCalculator(const std::string& method, const json& controller, const std::string& basename);
    ~EnergyCalculator();

    void setMolecule(const Mol& mol);

    void updateGeometry(const double* coord);
    void updateGeometry(const std::vector<double>& geometry);

    void updateGeometry(const Matrix& geometry);
    void updateGeometry(const Eigen::VectorXd& geometry);

    Matrix Gradient() const { return m_gradient; }

    double CalculateEnergy(bool gradient = false, bool verbose = false);

    bool HasNan() const { return m_containsNaN; }

    bool Error() const { return m_error; }

    QMInterface* Interface() const { return m_qminterface; }

    inline void setParameter(const json& parameter)
    {
        m_parameter = parameter;
        // if (m_qmdff != NULL)
        //     m_qmdff->setParameter(parameter);
    }

    inline json Parameter() const { return m_parameter; }

    // Claude Generated: Set geometry filename for parameter caching
    void setGeometryFile(const std::string& filename) { m_geometry_file = filename; }
    void setBasename(const std::string& basename) { m_geometry_file = basename + ".xyz"; }
    Vector Energies() const { return m_orbital_energies; }
    Vector OrbitalOccuptations() const { return m_orbital_occupation; }

    Vector Charges() const;
    Position Dipole() const;

    std::vector<std::vector<double>> BondOrders() const;
    int NumElectrons() const { return m_num_electrons; }
    Eigen::MatrixXd NumGrad();

private:
    void initializeCommon(const json& controller);
    int SwitchMethod(const std::string& method);

    // Claude Generated: Centralized method dispatch to reduce code duplication
    void dispatchMethodAction(const std::string& action, const Mol* mol = nullptr);

    // Claude Generated: Helper function to check compilation flags
    bool isCompiled(const std::string& flag) const;

    void InitialiseUFF();
    void CalculateUFF(bool gradient, bool verbose = false);

    void InitialiseTBlite();
    void CalculateTBlite(bool gradient, bool verbose = false);

    void InitialiseXTB();
    void CalculateXTB(bool gradient, bool verbose = false);

    void InitialiseUlysses();
    void CalculateUlysses(bool gradient, bool verbose = false);

    void CalculateQMInterface(bool gradient, bool verbose = false);

    void InitialiseD4();
    void CalculateD4(bool gradient, bool verbose = false);

    void InitialiseD3();
    void CalculateD3(bool gradient, bool verbose = false);

    void InitialiseQMDFF();
    void CalculateQMDFF(bool gradient, bool verbose = false);

    void InitialiseFF();
    void CalculateFF(bool gradient, bool verbose = false);

    // Claude Generated: Auto-save and auto-load parameter functionality
    void saveAutoParameters(const Mol& mol, const json& parameters);
    bool tryLoadAutoParameters(const Mol& mol);

    json m_controller;

    QMInterface* m_qminterface = NULL;
    Mol m_mol;
    // eigenUFF* m_uff = NULL;// will be switch_method 8
    // QMDFF* m_qmdff = NULL;// will be switch_method 7
    ForceField* m_forcefield = NULL; // will be switch_method 0
    EHT* m_eht = NULL; // will be switch_method 6

    StringList m_uff_methods = { "fuff" };
    StringList m_ff_methods = { "uff", "uff-d3", "qmdff", "cgfnff" };
    StringList m_qmdff_method = { "fqmdff" };
    StringList m_tblite_methods = { "ipea1", "gfn1", "gfn2" };
    StringList m_xtb_methods = { "gfnff", "xtb-gfn1", "xtb-gfn2" };
    StringList m_ulysses_methods = { "ugfn2", "GFN2L", "pm3", "PM3PDDG", "MNDOPDDG", "PM3BP", "RM1", "AM1", "MNDO", "MNDOd", "pm6" };

    StringList m_d3_methods = { "d3" };
    StringList m_d4_methods = { "d4" };

    std::function<void(bool, bool)> m_ecengine;
    Vector m_charges;
    std::function<Position()> m_dipole;
    std::function<std::vector<std::vector<double>>()> m_bonds;
    json m_parameter;
    std::string m_method, m_param_file, m_geometry_file, m_basename;
    Matrix m_geometry, m_gradient, m_molecular_orbitals;
    Vector m_orbital_energies, m_orbital_occupation;
    Vector m_xtb_gradient;
    double m_energy, m_Tele = 300;
    double *m_coord, *m_grad;

    int m_atoms, m_num_electrons = 0;
    int m_mult = 1, m_SCFmaxiter = 100;
    std::string m_solvent = "none";

    int m_gfn = 2;
    int* m_atom_type;

    bool m_initialised = false;
    bool m_containsNaN = false;
    bool m_error = false;
    bool m_writeparam = false;
};
