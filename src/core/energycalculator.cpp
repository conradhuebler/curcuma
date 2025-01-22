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

#ifdef USE_TBLITE
#include "src/core/tbliteinterface.h"
#endif

#ifdef USE_XTB
#include "src/core/xtbinterface.h"
#endif

#ifndef _WIN32
#include <filesystem>
namespace fs = std::filesystem;
#endif

#include <filesystem>
#include <functional>

#include "forcefieldgenerator.h"

#include "energycalculator.h"

EnergyCalculator::EnergyCalculator(const std::string& method, const json& controller)
    : m_method(method)
{
    std::transform(m_method.begin(), m_method.end(), m_method.begin(), [](unsigned char c) { return std::tolower(c); });

    m_controller = controller;
    if (controller.contains("param_file")) {
        m_param_file = controller["param_file"];
    }

    if (controller.contains("write_param")) {
        m_writeparam = controller["write_param"];
    }

    m_bonds = []() {
        return std::vector<std::vector<double>>{ {} };
    };

    switch (SwitchMethod(m_method)) {
    case 8:
        /*
             m_qmdff = new QMDFF(controller);
             if (m_parameter.size())
                 m_qmdff->setParameter(m_parameter);
             m_ecengine = [this](bool gradient, bool verbose) {
                 this->CalculateQMDFF(gradient, verbose);
             };
             */
        break;
    case 7:
        /*
            m_uff = new eigenUFF(controller);
            m_ecengine = [this](bool gradient, bool verbose) {
                this->CalculateUFF(gradient, verbose);
            };*/
        break;

    case 6:
        delete m_eht;
        break;

    case 5:
        // m_d4 = new DFTD4Interface(controller);
        //  m_ecengine = [this](bool gradient, bool verbose) {
        //      this->CalculateD4(gradient, verbose);
        //  };
        break;

    case 4:
        m_qminterface = new DFTD3Interface(controller);
        m_ecengine = [this](bool gradient, bool verbose) {
            this->CalculateD3(gradient, verbose);
        };
        break;

    case 3:
        m_qminterface = new UlyssesInterface(controller);
        m_qminterface->setMethod(m_method);
        m_ecengine = [this](bool gradient, bool verbose) {
            this->CalculateUlysses(gradient, verbose);
        };
        break;

    case 2:
#ifdef USE_XTB
        m_qminterface = new XTBInterface(controller);
        m_qminterface->setMethod(m_method);
        m_ecengine = [this](bool gradient, bool verbose) {
            this->CalculateXTB(gradient, verbose);
        };
#else
        std::cout << "XTB was not included ..." << std::endl;
        exit(1);
#endif
        /*
                m_bonds = [this]() {
                    return this->m_qminterface->BondOrders();
                };
        */
        break;

    case 1:

        m_qminterface = new TBLiteInterface(controller);
        m_qminterface->setMethod(m_method);

        m_ecengine = [this](bool gradient, bool verbose) {
            this->CalculateTBlite(gradient, verbose);
            m_error = this->m_qminterface->Error();
        };
        /*
        m_bonds = [this]() {
            return this->m_tblite->BondOrders();
        };
      */
        break;

    case 0:
    default:
        m_forcefield = new ForceField(controller);
        m_ecengine = [this](bool gradient, bool verbose) {
            this->CalculateFF(gradient, verbose);
        };
        break;
    }

    if (std::find(m_uff_methods.begin(), m_uff_methods.end(), m_method) != m_uff_methods.end()) { // UFF energy calculator requested

    } else if (std::find(m_tblite_methods.begin(), m_tblite_methods.end(), m_method) != m_tblite_methods.end()) { // TBLite energy calculator requested
#ifdef USE_TBLITE

#else
        std::cout << "TBlite was not included ..." << std::endl;
        exit(1);
#endif

    } else if (std::find(m_xtb_methods.begin(), m_xtb_methods.end(), m_method) != m_xtb_methods.end()) { // XTB energy calculator requested
#ifdef USE_XTB

#else
        std::cout << "XTB was not included ..." << std::endl;
        exit(1);
#endif
    } else if (std::find(m_ulysses_methods.begin(), m_ulysses_methods.end(), m_method) != m_ulysses_methods.end()) { // XTB energy calculator requested
#ifdef USE_XTB

#else
        std::cout << "XTB was not included ..." << std::endl;
        exit(1);
#endif
    } else if (std::find(m_d3_methods.begin(), m_d3_methods.end(), m_method) != m_d3_methods.end()) { // Just D4 energy calculator requested
#ifdef USE_D3

#else
        std::cout << "D4 was not included ..." << std::endl;
        exit(1);
#endif
    } else if (std::find(m_d4_methods.begin(), m_d4_methods.end(), m_method) != m_d4_methods.end()) { // Just D4 energy calculator requested
#ifdef USE_D4

#else
        std::cout << "D4 was not included ..." << std::endl;
        exit(1);
#endif
    } else if (std::find(m_qmdff_method.begin(), m_qmdff_method.end(), m_method) != m_qmdff_method.end()) { // Just D4 energy calculator requested

    } else if (std::find(m_ff_methods.begin(), m_ff_methods.end(), m_method) != m_ff_methods.end()) { // Just D4 energy calculator requested

    } else if (m_method.compare("eht") == 0) {

    } else { // Fall back to UFF?
        //    m_uff = new eigenUFF(controller);
    }
}

EnergyCalculator::~EnergyCalculator()
{
    switch (SwitchMethod(m_method)) {
    case 8:
        // delete m_qmdff;
        break;
    case 7:
        //  delete m_uff;
        break;

    case 6:
        delete m_eht;
        break;

    case 5:
        //    delete m_d4;
        break;

    case 4:
        delete m_qminterface;
        break;

    case 3:
        delete m_qminterface;
        break;

    case 2:
        delete m_qminterface;
        break;

    case 1:
        delete m_qminterface;
        break;

    case 0:
    default:
        delete m_forcefield;
        break;
    }
    // delete[] m_coord;
    // delete[] m_grad;
}

void EnergyCalculator::setMolecule(const Mol& mol)
{
    m_mol = mol;
    m_atoms = mol.m_number_atoms;
    m_geometry = mol.m_geometry;
    m_gradient = Eigen::MatrixXd::Zero(m_atoms, 3);
    m_xtb_gradient = Eigen::VectorXd::Zero(3 * m_atoms);

    /*
    std::vector<int> atoms = molecule.Atoms();
    m_coord = new double[3 * m_atoms];
    m_grad = new double[3 * m_atoms];
    m_gradient = Eigen::MatrixXd::Zero(m_atoms, 3);
    m_geometry = molecule.getGeometry();
    */

    switch (SwitchMethod(m_method)) {
    case 8:
        // m_qmdff->setMolecule(m_atoms, m_geometry);
        // m_qmdff->Initialise();
        break;
    case 7:

        break;

    case 6:
        m_eht->setMolecule(mol);
        m_eht->Initialise();
        break;

    case 5:
#ifdef USE_D4
        static_cast<DFTD4Interface*>(m_qminterface)->InitialiseMolecule(mol, 1 / au);
#endif
        break;

    case 4:
        static_cast<DFTD3Interface*>(m_qminterface)->InitialiseMolecule(mol.m_atoms);
        break;

    case 3:
        m_qminterface->InitialiseMolecule(mol);
        break;

    case 2:
        m_qminterface->InitialiseMolecule(mol);
        m_qminterface->setMethod(m_method);
        break;

    case 1:
        m_qminterface->InitialiseMolecule(mol);
        m_qminterface->setMethod(m_method);

        break;

    case 0:
    default:
        if (m_parameter.size() == 0) {
            if (!std::filesystem::exists(m_param_file)) {
                ForceFieldGenerator ff(m_controller);
                ff.setMolecule(mol);
                ff.Generate();
                m_parameter = ff.getParameter();
                if (m_writeparam) {
                    std::ofstream parameterfile("ff_param.json");
                    parameterfile << m_parameter;
                }
            } else {
                std::ifstream parameterfile(m_param_file);
                try {
                    parameterfile >> m_parameter;
                } catch (nlohmann::json::type_error& e) {
                } catch (nlohmann::json::parse_error& e) {
                }
            }
        }
        m_forcefield->setAtomTypes(mol.m_atoms);

        m_forcefield->setParameter(m_parameter);
        break;
    }
    m_initialised = true;
}

int EnergyCalculator::SwitchMethod(const std::string& method)
{
    int switch_method = 0;
    if (std::find(m_uff_methods.begin(), m_uff_methods.end(), m_method) != m_uff_methods.end()) { // UFF energy calculator requested
        switch_method = 0;
        return switch_method;
    } else if (std::find(m_tblite_methods.begin(), m_tblite_methods.end(), m_method) != m_tblite_methods.end()) { // TBLite energy calculator requested
#ifdef USE_TBLITE
        switch_method = 1;
        return switch_method;
#endif
        switch_method = 3;
    } else if (std::find(m_xtb_methods.begin(), m_xtb_methods.end(), m_method) != m_xtb_methods.end()) { // XTB energy calculator requested
#ifdef USE_XTB
        switch_method = 2;
        return switch_method;
#endif
        switch_method = 3;
    } else if (std::find(m_d3_methods.begin(), m_d3_methods.end(), m_method) != m_d3_methods.end()) { // D3 energy calculator requested
#ifdef USE_D3
        switch_method = 4;
#endif
    } else if (std::find(m_d4_methods.begin(), m_d4_methods.end(), m_method) != m_d4_methods.end()) { // D4 energy calculator requested
#ifdef USE_D4
        switch_method = 5;
#endif
    } else if (std::find(m_qmdff_method.begin(), m_qmdff_method.end(), m_method) != m_qmdff_method.end()) { //
        switch_method = 7;

    } else if (std::find(m_ff_methods.begin(), m_ff_methods.end(), m_method) != m_ff_methods.end()) { //
        switch_method = 0;
    } else if (m_method.compare("eht") == 0) {
        switch_method = 6;
    } else if (std::find(m_ulysses_methods.begin(), m_ulysses_methods.end(), m_method) != m_ulysses_methods.end()) {
        switch_method = 3;
    } else {
        switch_method = 0;
    }
    return switch_method;
}
void EnergyCalculator::updateGeometry(const Eigen::VectorXd& geometry)
{
#pragma message("Eigen::VectorXd ....")
    for (int i = 0; i < m_atoms; ++i) {
        m_geometry(i, 0) = geometry[3 * i + 0];
        m_geometry(i, 1) = geometry[3 * i + 1];
        m_geometry(i, 2) = geometry[3 * i + 2];
    }

    // m_containsNaN = std::isnan(m_geometry[m_atoms - 1][0]);
}

void EnergyCalculator::updateGeometry(const double* coord)
{

    for (int i = 0; i < m_atoms; ++i) {
        m_geometry(i, 0) = coord[3 * i + 0];
        m_geometry(i, 1) = coord[3 * i + 1];
        m_geometry(i, 2) = coord[3 * i + 2];
    }
    // m_containsNaN = std::isnan(m_geometry[m_atoms - 1][0]);
}

void EnergyCalculator::updateGeometry(const std::vector<double>& geometry)
{
    for (int i = 0; i < m_atoms; ++i) {
        m_geometry(i, 0) = geometry[3 * i + 0];
        m_geometry(i, 1) = geometry[3 * i + 1];
        m_geometry(i, 2) = geometry[3 * i + 2];
    }
    // m_containsNaN = std::isnan(m_geometry[m_atoms - 1][0]);
}

void EnergyCalculator::updateGeometry(const Matrix& geometry)
{
    m_geometry = geometry;
}

double EnergyCalculator::CalculateEnergy(bool gradient, bool verbose)
{
    m_ecengine(gradient, verbose);
    return m_energy;
}

void EnergyCalculator::CalculateUFF(bool gradient, bool verbose)
{
    /*
    m_uff->UpdateGeometry(m_geometry);
    m_energy = m_uff->Calculate(gradient, verbose);
    if (gradient) {
        m_gradient = m_uff->Gradient();
        // m_gradient = m_uff->NumGrad();
    }
    */
}

void EnergyCalculator::CalculateTBlite(bool gradient, bool verbose)
{
#ifdef USE_TBLITE
    m_qminterface->UpdateMolecule(m_geometry);
    m_energy = m_qminterface->Calculation(gradient, verbose);
    if (gradient)
        m_gradient = m_qminterface->Gradient();

#endif
}

void EnergyCalculator::CalculateXTB(bool gradient, bool verbose)
{
#ifdef USE_XTB
    m_qminterface->UpdateMolecule(m_geometry);
    m_energy = m_qminterface->Calculation(gradient, verbose);
    if (gradient)
        m_gradient = m_qminterface->Gradient();

#endif
}

void EnergyCalculator::CalculateUlysses(bool gradient, bool verbose)
{

    m_qminterface->UpdateMolecule(m_geometry);
    m_energy = m_qminterface->Calculation(gradient, verbose);
    if (gradient) {
        m_gradient = m_qminterface->Gradient();
    }
}

void EnergyCalculator::CalculateD3(bool gradient, bool verbose)
{
#ifdef USE_D3
    for (int i = 0; i < m_atoms; ++i) {
        static_cast<DFTD3Interface*>(m_qminterface)->UpdateAtom(i, m_geometry(i, 0), m_geometry(i, 1), m_geometry(i, 2));
    }
    if (gradient) {
        m_energy = m_qminterface->Calculation(m_xtb_gradient.data());
        for (int i = 0; i < m_atoms; ++i) {
            m_gradient(i, 0) = m_xtb_gradient[3 * i + 0] * au;
            m_gradient(i, 1) = m_xtb_gradient[3 * i + 1] * au;
            m_gradient(i, 2) = m_xtb_gradient[3 * i + 2] * au;
        }
    } else
        m_energy = m_qminterface->Calculation(false, verbose);
#endif
}

void EnergyCalculator::CalculateD4(bool gradient, bool verbose)
{
#ifdef USE_D4
    for (int i = 0; i < m_atoms; ++i) {
        static_cast<DFTD4Interface*>(m_qminterface)->UpdateAtom(i, m_geometry(i, 0) / au, m_geometry(i, 1) / au, m_geometry(i, 2) / au);
    }
    if (gradient) {
        m_energy = m_qminterface->Calculation(m_grad);
        for (int i = 0; i < m_atoms; ++i) {
            m_gradient(i, 0) = m_grad[3 * i + 0] * au;
            m_gradient(i, 1) = m_grad[3 * i + 1] * au;
            m_gradient(i, 2) = m_grad[3 * i + 2] * au;
        }
    } else
        m_energy = m_qminterface->Calculation();
#endif
}

void EnergyCalculator::CalculateQMDFF(bool gradient, bool verbose)
{
    /*
    m_qmdff->UpdateGeometry(m_geometry);
    m_energy = m_qmdff->Calculate(gradient, verbose);
    if (gradient) {
        m_gradient = m_qmdff->Gradient();
    }*/
}

void EnergyCalculator::CalculateFF(bool gradient, bool verbose)
{
    m_forcefield->UpdateGeometry(m_geometry);
    m_energy = m_forcefield->Calculate(gradient, verbose);
    if (gradient) {
        m_gradient = m_forcefield->Gradient();
    }
}

Vector EnergyCalculator::Charges() const
{
    if (m_qminterface == nullptr)
        return Vector{};
    return m_qminterface->Charges();
}

Position EnergyCalculator::Dipole() const
{
    if (m_qminterface == nullptr)
        return Position{};
    return m_qminterface->Dipole();
}

std::vector<std::vector<double>> EnergyCalculator::BondOrders() const
{
    return m_bonds();
}
