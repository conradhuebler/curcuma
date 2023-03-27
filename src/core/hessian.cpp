/*
 * < General Calculator for the Hessian>
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
 *
 */

#include <Eigen/Dense>

#include "external/CxxThreadPool/include/CxxThreadPool.h"

#include "src/core/energycalculator.h"

#include "hessian.h"

HessianThread::HessianThread(const std::string& method, const json& controller, int i, int j, int xi, int xj)
    : m_method(method)
    , m_controller(controller)
    , m_i(i)
    , m_j(j)
    , m_xi(xi)
    , m_xj(xj)
{
    setAutoDelete(true);
}

HessianThread::~HessianThread()
{
}

void HessianThread::setMolecule(const Molecule& molecule)
{
    m_molecule = molecule;
    m_geom_ip_jp = m_molecule.Coords();
    m_geom_im_jp = m_molecule.Coords();
    m_geom_ip_jm = m_molecule.Coords();
    m_geom_im_jm = m_molecule.Coords();

    m_geom_ip_jp[m_i][m_xi] += m_d;
    m_geom_ip_jp[m_j][m_xj] += m_d;

    m_geom_im_jp[m_i][m_xi] -= m_d;
    m_geom_im_jp[m_j][m_xj] += m_d;

    m_geom_ip_jm[m_i][m_xi] += m_d;
    m_geom_ip_jm[m_j][m_xj] -= m_d;

    m_geom_im_jm[m_i][m_xi] -= m_d;
    m_geom_im_jm[m_j][m_xj] -= m_d;
}

void HessianThread::updateGeometry(const double* coord)
{
}

int HessianThread::execute()
{
    EnergyCalculator energy(m_method, m_controller);
    energy.setMolecule(m_molecule);
    double d2 = 1 / (4 * m_d * m_d);

    energy.updateGeometry(m_geom_ip_jp);
    double energy_ip_jp = energy.CalculateEnergy(false, false);

    energy.updateGeometry(m_geom_im_jp);
    double energy_im_jp = energy.CalculateEnergy(false, false);

    energy.updateGeometry(m_geom_ip_jm);
    double energy_ip_jm = energy.CalculateEnergy(false, false);

    energy.updateGeometry(m_geom_im_jm);
    double energy_im_jm = energy.CalculateEnergy(false, false);
    m_dd = d2 * (energy_ip_jp - energy_im_jp - energy_ip_jm + energy_im_jm);

    return 0;
}

Hessian::Hessian(const std::string& method, const json& controller, int threads)
    : m_method(method)
    , m_controller(controller)
    , m_threads(threads)
{
}

void Hessian::setMolecule(const Molecule& molecule)
{
    m_molecule = molecule;
}

void Hessian::CalculateHessian()
{
    m_hessian = Eigen::MatrixXd::Ones(3 * m_molecule.AtomCount(), 3 * m_molecule.AtomCount());
    auto geometry = m_molecule.Coords();
    CxxThreadPool* pool = new CxxThreadPool;
    pool->setActiveThreadCount(m_threads);
    for (int i = 0; i < m_molecule.AtomCount(); ++i) {
        for (int j = 0; j < m_molecule.AtomCount(); ++j) {
            double mass = 1 / sqrt(Elements::AtomicMass[m_molecule.Atoms()[i]] * Elements::AtomicMass[m_molecule.Atoms()[j]]);
            for (int xi = 0; xi < 3; ++xi)
                for (int xj = 0; xj < 3; ++xj) {
                    HessianThread* thread = new HessianThread(m_method, m_controller, i, j, xi, xj);
                    thread->setMolecule(m_molecule);
                    pool->addThread(thread);
                    m_hessian(3 * i + xi, 3 * j + xj) = mass;
                }
        }
    }
    pool->DynamicPool(2);
    pool->StartAndWait();
    for (auto* t : pool->Finished()) {
        HessianThread* thread = static_cast<HessianThread*>(t);
        m_hessian(3 * thread->I() + thread->XI(), 3 * thread->J() + thread->XJ()) *= thread->DD();
    }
    Eigen::SelfAdjointEigenSolver<Geometry> diag_I;
    diag_I.compute(m_hessian);
    std::cout << diag_I.eigenvalues() << std::endl;
    std::cout << diag_I.eigenvalues().cwiseSqrt() * 2720.7 /* /2.0/3.14*amu2au */ << std::endl;
    delete pool;
}
