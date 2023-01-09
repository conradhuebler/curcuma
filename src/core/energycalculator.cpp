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
    if (std::find(m_uff_methods.begin(), m_uff_methods.end(), m_method) != m_uff_methods.end()) { // UFF energy calculator requested
        m_uff = new UFF(controller);
        m_ecengine = [this](bool gradient, bool verbose) {
            this->CalculateUFF(gradient, verbose);
        };
    } else if (std::find(m_tblite_methods.begin(), m_tblite_methods.end(), m_method) != m_tblite_methods.end()) { // TBLite energy calculator requested
#ifdef USE_TBLITE
        m_tblite = new TBLiteInterface(controller);
        m_ecengine = [this](bool gradient, bool verbose) {
            this->CalculateTBlite(gradient, verbose);
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
    } else { // Fall back to UFF?
        m_uff = new UFF(controller);
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
    } else { // Fall back to UFF?
        delete m_uff;
    }
}

void EnergyCalculator::setMolecule(const Molecule& molecule)
{
    m_atoms = molecule.AtomCount();

    // m_atom_type[m_atoms];
    std::vector<int> atoms = molecule.Atoms();
    m_coord[3 * m_atoms];
    std::vector<std::array<double, 3>> geom(m_atoms);
    for (int i = 0; i < m_atoms; ++i) {
        std::pair<int, Position> atom = molecule.Atom(i);
        geom[i][0] = atom.second(0);
        geom[i][1] = atom.second(1);
        geom[i][2] = atom.second(2);
        m_gradient.push_back({ 0, 0, 0 });
    }
    m_geometry = geom;
    if (std::find(m_uff_methods.begin(), m_uff_methods.end(), m_method) != m_uff_methods.end()) { // UFF energy calculator requested
        m_uff->setMolecule(atoms, geom);
        m_uff->Initialise();
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
    } else { // Fall back to UFF?
    }
    m_initialised = true;
}

void EnergyCalculator::updateGeometry(const Eigen::VectorXd& geometry)
{
    for (int i = 0; i < m_atoms; ++i) {
        m_geometry[i][0] = geometry[3 * i + 0];
        m_geometry[i][1] = geometry[3 * i + 1];
        m_geometry[i][2] = geometry[3 * i + 2];
    }
}

void EnergyCalculator::updateGeometry(const double* coord)
{
    for (int i = 0; i < m_atoms; ++i) {
        m_geometry[i][0] = coord[3 * i + 0];
        m_geometry[i][1] = coord[3 * i + 1];
        m_geometry[i][2] = coord[3 * i + 2];
    }
}

void EnergyCalculator::updateGeometry(const std::vector<double>& geometry)
{
    for (int i = 0; i < m_atoms; ++i) {
        m_geometry[i][0] = geometry[3 * i + 0];
        m_geometry[i][1] = geometry[3 * i + 1];
        m_geometry[i][2] = geometry[3 * i + 2];
    }
}
void EnergyCalculator::updateGeometry(const std::vector<std::array<double, 3>>& geometry)
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
    double coord[3 * m_atoms];
    for (int i = 0; i < m_atoms; ++i) {
        coord[3 * i + 0] = m_geometry[i][0] / au;
        coord[3 * i + 1] = m_geometry[i][1] / au;
        coord[3 * i + 2] = m_geometry[i][2] / au;
    }
    m_tblite->UpdateMolecule(coord);

    if (gradient) {
        double grad[3 * m_atoms];
        m_energy = m_tblite->GFNCalculation(m_gfn, grad);
        for (int i = 0; i < m_atoms; ++i) {
            m_gradient[i][0] = grad[3 * i + 0] * au;
            m_gradient[i][1] = grad[3 * i + 1] * au;
            m_gradient[i][2] = grad[3 * i + 2] * au;
        }
    } else
        m_energy = m_tblite->GFNCalculation(m_gfn);
#endif
}

void EnergyCalculator::CalculateXTB(bool gradient, bool verbose)
{
#ifdef USE_XTB
    double coord[3 * m_atoms];
    for (int i = 0; i < m_atoms; ++i) {
        coord[3 * i + 0] = m_geometry[i][0] / au;
        coord[3 * i + 1] = m_geometry[i][1] / au;
        coord[3 * i + 2] = m_geometry[i][2] / au;
    }
    m_xtb->UpdateMolecule(coord);

    if (gradient) {
        double grad[3 * m_atoms];
        m_energy = m_xtb->GFNCalculation(m_gfn, grad);
        for (int i = 0; i < m_atoms; ++i) {
            m_gradient[i][0] = grad[3 * i + 0] * au;
            m_gradient[i][1] = grad[3 * i + 1] * au;
            m_gradient[i][2] = grad[3 * i + 2] * au;
        }
    } else
        m_energy = m_xtb->GFNCalculation(m_gfn);
#endif
}

void EnergyCalculator::CalculateD3(bool gradient, bool verbose)
{
#ifdef USE_D3
    for (int i = 0; i < m_atoms; ++i) {
        m_d3->UpdateAtom(i, m_geometry[i][0], m_geometry[i][1], m_geometry[i][2]);
    }
    if (gradient) {
        double grad[3 * m_atoms];
        m_energy = m_d3->DFTD3Calculation(grad);
        for (int i = 0; i < m_atoms; ++i) {
            m_gradient[i][0] = grad[3 * i + 0] * au;
            m_gradient[i][1] = grad[3 * i + 1] * au;
            m_gradient[i][2] = grad[3 * i + 2] * au;
        }
    } else
        m_energy = m_d3->DFTD3Calculation();
#endif
}

void EnergyCalculator::CalculateD4(bool gradient, bool verbose)
{
#ifdef USE_D4
    for (int i = 0; i < m_atoms; ++i) {
        m_d4->UpdateAtom(i, m_geometry[i][0] / au, m_geometry[i][1] / au, m_geometry[i][2] / au);
    }
    if (gradient) {
        double grad[3 * m_atoms];
        m_energy = m_d4->DFTD4Calculation(grad);
        for (int i = 0; i < m_atoms; ++i) {
            m_gradient[i][0] = grad[3 * i + 0] * au;
            m_gradient[i][1] = grad[3 * i + 1] * au;
            m_gradient[i][2] = grad[3 * i + 2] * au;
        }
    } else
        m_energy = m_d4->DFTD4Calculation();
#endif
}

void EnergyCalculator::getGradient(double* gradient)
{
    for (int i = 0; i < m_atoms; ++i) {
        gradient[3 * i + 0] = m_gradient[i][0];
        gradient[3 * i + 1] = m_gradient[i][1];
        gradient[3 * i + 2] = m_gradient[i][2];
    }
}
