/*
 * < General Energy and Gradient Calculator >
 * Copyright (C) 2022 - 2023 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#ifdef USE_TBLITE
#include "src/core/tbliteinterface.h"
#endif

#ifdef USE_XTB
#include "src/core/xtbinterface.h"
#endif

#include <functional>

#include "energycalculator.h"

EnergyCalculator::EnergyCalculator(const std::string& method, const json& controller)
    : m_method(method)
{
    m_charges = []() {
        return std::vector<double>{};
    };
    m_dipole = []() {
        return std::vector<double>{};
    };
    m_bonds = []() {
        return std::vector<std::vector<double>>{ {} };
    };

    if (std::find(m_uff_methods.begin(), m_uff_methods.end(), m_method) != m_uff_methods.end()) { // UFF energy calculator requested
        m_uff = new eigenUFF(controller);
        m_ecengine = [this](bool gradient, bool verbose) {
            this->CalculateUFF(gradient, verbose);
        };
    } else if (std::find(m_tblite_methods.begin(), m_tblite_methods.end(), m_method) != m_tblite_methods.end()) { // TBLite energy calculator requested
#ifdef USE_TBLITE
        m_tblite = new TBLiteInterface(controller);
        m_ecengine = [this](bool gradient, bool verbose) {
            this->CalculateTBlite(gradient, verbose);
            m_error = this->m_tblite->Error();
        };
        m_charges = [this]() {
            return this->m_tblite->Charges();
        };
        m_dipole = [this]() {
            return this->m_tblite->Dipole();
        };
        m_bonds = [this]() {
            return this->m_tblite->BondOrders();
        };
#else
        std::cout << "TBlite was not included ..." << std::endl;
        exit(1);
#endif

    } else if (std::find(m_xtb_methods.begin(), m_xtb_methods.end(), m_method) != m_xtb_methods.end()) { // XTB energy calculator requested
#ifdef USE_XTB
        m_xtb = new XTBInterface(controller);
        m_ecengine = [this](bool gradient, bool verbose) {
            this->CalculateXTB(gradient, verbose);
        };
        m_charges = [this]() {
            return this->m_xtb->Charges();
        };
        m_dipole = [this]() {
            return this->m_xtb->Dipole();
        };
        m_bonds = [this]() {
            return this->m_xtb->BondOrders();
        };
#else
        std::cout << "XTB was not included ..." << std::endl;
        exit(1);
#endif
    } else if (std::find(m_d3_methods.begin(), m_d3_methods.end(), m_method) != m_d3_methods.end()) { // Just D4 energy calculator requested
#ifdef USE_D3
        m_d3 = new DFTD3Interface(controller);
        m_ecengine = [this](bool gradient, bool verbose) {
            this->CalculateD3(gradient, verbose);
        };
#else
        std::cout << "D4 was not included ..." << std::endl;
        exit(1);
#endif
    } else if (std::find(m_d4_methods.begin(), m_d4_methods.end(), m_method) != m_d4_methods.end()) { // Just D4 energy calculator requested
#ifdef USE_D4
        m_d4 = new DFTD4Interface(controller);
        m_ecengine = [this](bool gradient, bool verbose) {
            this->CalculateD4(gradient, verbose);
        };
#else
        std::cout << "D4 was not included ..." << std::endl;
        exit(1);
#endif
    } else if (std::find(m_qmdff_method.begin(), m_qmdff_method.end(), m_method) != m_qmdff_method.end()) { // Just D4 energy calculator requested
        m_qmdff = new QMDFF(controller);
        if (m_parameter.size())
            m_qmdff->setParameter(m_parameter);
        m_ecengine = [this](bool gradient, bool verbose) {
            this->CalculateQMDFF(gradient, verbose);
        };

    } else if (std::find(m_ff_methods.begin(), m_ff_methods.end(), m_method) != m_ff_methods.end()) { // Just D4 energy calculator requested
        m_forcefield = new ForceField(controller);
        std::string file = controller["param"];
        nlohmann::json parameter;
        std::ifstream parameterfile(file);
        try {
            parameterfile >> parameter;
        } catch (nlohmann::json::type_error& e) {
        } catch (nlohmann::json::parse_error& e) {
        }
        m_forcefield->setParameter(parameter);
        m_ecengine = [this](bool gradient, bool verbose) {
            this->CalculateFF(gradient, verbose);
        };

    } else { // Fall back to UFF?
        m_uff = new eigenUFF(controller);
    }
}

EnergyCalculator::~EnergyCalculator()
{
    if (std::find(m_uff_methods.begin(), m_uff_methods.end(), m_method) != m_uff_methods.end()) { // UFF energy calculator requested
        delete m_uff;
    } else if (std::find(m_tblite_methods.begin(), m_tblite_methods.end(), m_method) != m_tblite_methods.end()) { // TBLite energy calculator requested
#ifdef USE_TBLITE
        delete m_tblite;
#endif

    } else if (std::find(m_xtb_methods.begin(), m_xtb_methods.end(), m_method) != m_xtb_methods.end()) { // XTB energy calculator requested
#ifdef USE_XTB
        delete m_xtb;
#endif
    } else if (std::find(m_d3_methods.begin(), m_d3_methods.end(), m_method) != m_d3_methods.end()) { // XTB energy calculator requested
#ifdef USE_D3
        delete m_d3;
#endif
    } else if (std::find(m_d4_methods.begin(), m_d4_methods.end(), m_method) != m_d4_methods.end()) { // XTB energy calculator requested
#ifdef USE_D4
        delete m_d4;
#endif
    } else if (std::find(m_qmdff_method.begin(), m_qmdff_method.end(), m_method) != m_qmdff_method.end()) { // Just D4 energy calculator requested
        delete m_qmdff;
    } else { // Fall back to UFF?
        delete m_uff;
    }
}

void EnergyCalculator::setMolecule(const Molecule& molecule)
{
    m_atoms = molecule.AtomCount();

    // m_atom_type[m_atoms];
    std::vector<int> atoms = molecule.Atoms();
    m_coord = new double[3 * m_atoms];
    m_grad = new double[3 * m_atoms];
    // std::vector<std::array<double, 3>> geom(m_atoms);
    m_gradient = Eigen::MatrixXd::Zero(m_atoms, 3);
    m_geometry = Eigen::MatrixXd::Zero(m_atoms, 3);
    for (int i = 0; i < m_atoms; ++i) {
        std::pair<int, Position> atom = molecule.Atom(i);
        m_geometry(i, 0) = atom.second(0);
        m_geometry(i, 1) = atom.second(1);
        m_geometry(i, 2) = atom.second(2);
        /*
        geom[i][0] = atom.second(0);
        geom[i][1] = atom.second(1);
        geom[i][2] = atom.second(2);*/
        // m_gradient.push_back({ 0, 0, 0 });
    }
    // m_geometry = geom;
    if (std::find(m_uff_methods.begin(), m_uff_methods.end(), m_method) != m_uff_methods.end()) { // UFF energy calculator requested
        m_uff->setMolecule(atoms, m_geometry);
        m_uff->Initialise(molecule.Bonds());
    } else if (std::find(m_tblite_methods.begin(), m_tblite_methods.end(), m_method) != m_tblite_methods.end()) { // TBLite energy calculator requested
#ifdef USE_TBLITE
        m_tblite->InitialiseMolecule(molecule);
        if (m_method.compare("gfn1") == 0)
            m_gfn = 1;
        else if (m_method.compare("gfn2") == 0)
            m_gfn = 2;
#endif

    } else if (std::find(m_xtb_methods.begin(), m_xtb_methods.end(), m_method) != m_xtb_methods.end()) { // XTB energy calculator requested
#ifdef USE_XTB
        m_xtb->InitialiseMolecule(molecule);
        if (m_method.compare("xtb-gfn1") == 0)
            m_gfn = 1;
        else if (m_method.compare("xtb-gfn2") == 0)
            m_gfn = 2;
        else if (m_method.compare("xtb-gfn0") == 0)
            m_gfn = 0;
        else if (m_method.compare("gfnff") == 0)
            m_gfn = 66;
#endif
    } else if (std::find(m_d3_methods.begin(), m_d3_methods.end(), m_method) != m_d3_methods.end()) { // D3 energy calculator requested
#ifdef USE_D3
        m_d3->InitialiseMolecule(molecule.Atoms());
#endif
    } else if (std::find(m_d4_methods.begin(), m_d4_methods.end(), m_method) != m_d4_methods.end()) { // D4 energy calculator requested
#ifdef USE_D4
        m_d4->InitialiseMolecule(molecule, 1 / au);
#endif
    } else if (std::find(m_qmdff_method.begin(), m_qmdff_method.end(), m_method) != m_qmdff_method.end()) { //
        m_qmdff->setMolecule(atoms, m_geometry);
        m_qmdff->Initialise();

    } else if (std::find(m_ff_methods.begin(), m_ff_methods.end(), m_method) != m_ff_methods.end()) { //
        /*
        m_forcefield->setMolecule(atoms, geom);
        m_forcefield->Initialise();
        */
    } else { // Fall back to UFF?
    }
    m_initialised = true;
}

void EnergyCalculator::updateGeometry(const Eigen::VectorXd& geometry)
{
    for (int i = 0; i < m_atoms; ++i) {
        m_coord[3 * i + 0] = geometry[3 * i + 0] / au;
        m_coord[3 * i + 1] = geometry[3 * i + 1] / au;
        m_coord[3 * i + 2] = geometry[3 * i + 2] / au;

        m_geometry(i, 0) = geometry[3 * i + 0];
        m_geometry(i, 1) = geometry[3 * i + 1];
        m_geometry(i, 2) = geometry[3 * i + 2];
    }
    // m_containsNaN = std::isnan(m_geometry[m_atoms - 1][0]);
}

void EnergyCalculator::updateGeometry(const double* coord)
{
    for (int i = 0; i < m_atoms; ++i) {
        m_coord[3 * i + 0] = coord[3 * i + 0] / au;
        m_coord[3 * i + 1] = coord[3 * i + 1] / au;
        m_coord[3 * i + 2] = coord[3 * i + 2] / au;

        m_geometry(i, 0) = coord[3 * i + 0];
        m_geometry(i, 1) = coord[3 * i + 1];
        m_geometry(i, 2) = coord[3 * i + 2];
    }
    // m_containsNaN = std::isnan(m_geometry[m_atoms - 1][0]);
}

void EnergyCalculator::updateGeometry(const std::vector<double>& geometry)
{
    for (int i = 0; i < m_atoms; ++i) {
        m_coord[3 * i + 0] = geometry[3 * i + 0] / au;
        m_coord[3 * i + 1] = geometry[3 * i + 1] / au;
        m_coord[3 * i + 2] = geometry[3 * i + 2] / au;

        m_geometry(i, 0) = geometry[3 * i + 0];
        m_geometry(i, 1) = geometry[3 * i + 1];
        m_geometry(i, 2) = geometry[3 * i + 2];
    }
    // m_containsNaN = std::isnan(m_geometry[m_atoms - 1][0]);
}
void EnergyCalculator::updateGeometry(const std::vector<std::array<double, 3>>& geometry)
{
    for (int i = 0; i < m_atoms; ++i) {
        m_coord[3 * i + 0] = geometry[i][0] / au;
        m_coord[3 * i + 1] = geometry[i][1] / au;
        m_coord[3 * i + 2] = geometry[i][2] / au;

        m_geometry(i, 0) = geometry[i][0];
        m_geometry(i, 1) = geometry[i][1];
        m_geometry(i, 2) = geometry[i][2];
    }
}

double EnergyCalculator::CalculateEnergy(bool gradient, bool verbose)
{
    m_ecengine(gradient, verbose);
    return m_energy;
}

void EnergyCalculator::CalculateUFF(bool gradient, bool verbose)
{
    m_uff->UpdateGeometry(m_geometry);
    m_energy = m_uff->Calculate(gradient, verbose);
    if (gradient) {
        m_gradient = m_uff->Gradient();
        // m_gradient = m_uff->NumGrad();
    }
}

void EnergyCalculator::CalculateTBlite(bool gradient, bool verbose)
{
#ifdef USE_TBLITE
    m_tblite->UpdateMolecule(m_coord);

    if (gradient) {
        m_energy = m_tblite->GFNCalculation(m_gfn, m_grad);
        for (int i = 0; i < m_atoms; ++i) {
            m_gradient(i, 0) = m_grad[3 * i + 0] * au;
            m_gradient(i, 1) = m_grad[3 * i + 1] * au;
            m_gradient(i, 2) = m_grad[3 * i + 2] * au;
        }
    } else
        m_energy = m_tblite->GFNCalculation(m_gfn);
#endif
}

void EnergyCalculator::CalculateXTB(bool gradient, bool verbose)
{
#ifdef USE_XTB
    m_xtb->UpdateMolecule(m_coord);

    if (gradient) {
        m_energy = m_xtb->GFNCalculation(m_gfn, m_grad);
        for (int i = 0; i < m_atoms; ++i) {
            m_gradient(i, 0) = m_grad[3 * i + 0] * au;
            m_gradient(i, 1) = m_grad[3 * i + 1] * au;
            m_gradient(i, 2) = m_grad[3 * i + 2] * au;
        }
    } else
        m_energy = m_xtb->GFNCalculation(m_gfn);
#endif
}

void EnergyCalculator::CalculateD3(bool gradient, bool verbose)
{
#ifdef USE_D3
    for (int i = 0; i < m_atoms; ++i) {
        m_d3->UpdateAtom(i, m_geometry(i, 0), m_geometry(i, 1), m_geometry(i, 2));
    }
    if (gradient) {
        m_energy = m_d3->DFTD3Calculation(m_grad);
        for (int i = 0; i < m_atoms; ++i) {
            m_gradient(i, 0) = m_grad[3 * i + 0] * au;
            m_gradient(i, 1) = m_grad[3 * i + 1] * au;
            m_gradient(i, 2) = m_grad[3 * i + 2] * au;
        }
    } else
        m_energy = m_d3->DFTD3Calculation();
#endif
}

void EnergyCalculator::CalculateD4(bool gradient, bool verbose)
{
#ifdef USE_D4
    for (int i = 0; i < m_atoms; ++i) {
        m_d4->UpdateAtom(i, m_geometry(i, 0) / au, m_geometry(i, 1) / au, m_geometry(i, 2) / au);
    }
    if (gradient) {
        m_energy = m_d4->DFTD4Calculation(m_grad);
        for (int i = 0; i < m_atoms; ++i) {
            m_gradient(i, 0) = m_grad[3 * i + 0] * au;
            m_gradient(i, 1) = m_grad[3 * i + 1] * au;
            m_gradient(i, 2) = m_grad[3 * i + 2] * au;
        }
    } else
        m_energy = m_d4->DFTD4Calculation();
#endif
}

void EnergyCalculator::CalculateQMDFF(bool gradient, bool verbose)
{
    m_qmdff->UpdateGeometry(m_geometry);
    m_energy = m_qmdff->Calculate(gradient, verbose);
    if (gradient) {
        m_gradient = m_qmdff->Gradient();
    }
}

void EnergyCalculator::CalculateFF(bool gradient, bool verbose)
{
    m_forcefield->UpdateGeometry(m_geometry);
    m_energy = m_forcefield->Calculate(gradient, verbose);
    if (gradient) {
        m_gradient = m_forcefield->Gradient();
    }
}

void EnergyCalculator::getGradient(double* gradient)
{
    for (int i = 0; i < m_atoms; ++i) {
        gradient[3 * i + 0] = m_gradient(i, 0);
        gradient[3 * i + 1] = m_gradient(i, 1);
        gradient[3 * i + 2] = m_gradient(i, 2);
    }
}

Matrix EnergyCalculator::Gradient() const
{
    return m_gradient;
}

std::vector<double> EnergyCalculator::Charges() const
{
    return m_charges();
}

std::vector<double> EnergyCalculator::Dipole() const
{
    return m_dipole();
}

std::vector<std::vector<double>> EnergyCalculator::BondOrders() const
{
    return m_bonds();
}
