/*
 * < General Calculator for the Hessian and thermodynamics>
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

HessianThread::HessianThread(const json& controller, int i, int j, int xi, int xj, bool fullnumerical)
    : m_controller(controller)
    , m_i(i)
    , m_j(j)
    , m_xi(xi)
    , m_xj(xj)
    , m_fullnumerical(fullnumerical)
{
    setAutoDelete(true);
    m_method = m_controller["method"];

    if (m_fullnumerical)
        m_schema = [this, i, j, xi, xj]() {
            this->Numerical();
        };
    else
        m_schema = [this]() {
            this->Seminumerical();
        };
}

HessianThread::~HessianThread()
{
}

void HessianThread::setMolecule(const Molecule& molecule)
{
    m_molecule = molecule;
}

int HessianThread::execute()
{
    m_schema();
    return 0;
}

void HessianThread::Numerical()
{
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

    EnergyCalculator energy(m_method, m_controller);
    energy.setMolecule(m_molecule);
    energy.setParameter(m_parameter);

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
}
void HessianThread::Seminumerical()
{
    m_geom_ip_jp = m_molecule.Coords();
    m_geom_im_jp = m_molecule.Coords();

    m_geom_ip_jp[m_i][m_xi] += m_d;
    // std::cout << m_controller << std::endl;

    EnergyCalculator energy(m_method, m_controller);
    energy.setMolecule(m_molecule);
    energy.setParameter(m_parameter);

    energy.updateGeometry(m_geom_ip_jp);
    energy.CalculateEnergy(true, false);
    Matrix gradientp = energy.Gradient();

    m_geom_im_jp[m_i][m_xi] -= m_d;

    energy.updateGeometry(m_geom_im_jp);
    energy.CalculateEnergy(true, false);
    Matrix gradientm = energy.Gradient();
    m_gradient = (gradientp - gradientm) / (2 * m_d) / au / au;
}

Hessian::Hessian(const std::string& method, const json& controller, bool silent)
    : CurcumaMethod(HessianJson, controller, silent)
    , m_method(method)
    , m_controller(controller)

{
    UpdateController(controller);
    m_threads = m_controller["threads"];
    /* Yeah, thats not really correct, but it works a bit */
    if (m_method.compare("gfnff") == 0) {
        m_scale_functions = [](double val) -> double {
            return val * 2720.57 - 0.0338928;
        };
    } else {
        m_scale_functions = [](double val) -> double {
            return val * 2720.57 - 0.0338928;
        };
    }
}

Hessian::Hessian(const json& controller, bool silent)
    : CurcumaMethod(HessianJson, controller, silent)
    , m_controller(controller)
{
    UpdateController(controller);
    /* Yeah, thats not really correct, but it works a bit */
    if (m_method.compare("gfnff") == 0) {
        m_scale_functions = [](double val) -> double {
            return val * 2720.57 - 0.0338928;
        };
    } else {
        m_scale_functions = [](double val) -> double {
            return val * 2720.57 - 0.0338928;
        };
    }
}

void Hessian::LoadControlJson()
{
    m_hess_calc = Json2KeyWord<bool>(m_defaults, "hess_calc");
    m_write_file = Json2KeyWord<std::string>(m_defaults, "hess_write_file");
    m_hess_read = Json2KeyWord<bool>(m_defaults, "hess_read");
    m_read_file = Json2KeyWord<std::string>(m_defaults, "hess_read_file");
    m_read_xyz = Json2KeyWord<std::string>(m_defaults, "hess_read_xyz");

    m_write_file = Json2KeyWord<std::string>(m_defaults, "hess_write_file");
    if (m_hess_read) {
        m_hess_calc = false;
    }

    m_freq_scale = Json2KeyWord<double>(m_defaults, "freq_scale");
    m_thermo = Json2KeyWord<double>(m_defaults, "thermo");
    m_freq_cutoff = Json2KeyWord<double>(m_defaults, "freq_cutoff");
    m_hess = Json2KeyWord<int>(m_defaults, "hess");
    m_method = Json2KeyWord<std::string>(m_defaults, "method");
}

void Hessian::setMolecule(const Molecule& molecule)
{
    m_molecule = molecule;
}

void Hessian::LoadMolecule(const std::string& file)
{
    m_molecule = Files::LoadFile(file);
    m_atom_count = m_molecule.AtomCount();
}

void Hessian::LoadHessian(const std::string& file)
{
    if (file.compare("hessian") == 0) {
        std::ifstream f(file);
        m_hessian = Matrix::Ones(3 * m_atom_count, 3 * m_atom_count);
        int row = 0, column = 0;
        for (std::string line; getline(f, line);) {
            if (line.compare("$hessian") == 0)
                continue;
            StringList entries = Tools::SplitString(line);
            for (int i = 0; i < entries.size(); ++i) {
                if (!Tools::isDouble(entries[i]))
                    continue;
                m_hessian(row, column) = std::stod(entries[i]) / au / au; // hessian files are some units ...
                column++;
                if (column == 3 * m_atom_count) {
                    column = 0;
                    row++;
                }
            }
        }
    } else {
        // will be json file
    }
}

void Hessian::start()
{
    if (m_hess_calc) {
        if (m_hess == 1) {
            CalculateHessianSemiNumerical();
        } else {
            CalculateHessianNumerical();
        }
    } else {
        LoadMolecule(m_read_xyz);
        LoadHessian(m_read_file);
    }

    m_frequencies = ConvertHessian(m_hessian);

    // PrintVibrationsPlain(eigenvalues);

    auto hessian2 = ProjectHessian(m_hessian);
    auto projected = ConvertHessian(hessian2);
    if (!m_silent)
        PrintVibrations(m_frequencies, projected);
}

Matrix Hessian::ProjectHessian(const Matrix& hessian)
{
    Matrix D = Eigen::MatrixXd::Random(3 * m_molecule.AtomCount(), 3 * m_molecule.AtomCount());
    Eigen::RowVector3d ex = { 1, 0, 0 };
    Eigen::RowVector3d ey = { 0, 1, 0 };
    Eigen::RowVector3d ez = { 0, 0, 1 };

    for (int i = 0; i < 3 * m_molecule.AtomCount(); ++i) {
        D(i, 0) = int(i % 3 == 0);
        D(i, 2) = int((i + 1) % 3 == 0);
        D(i, 1) = int((i + 2) % 3 == 0);
    }

    for (int i = 0; i < m_molecule.AtomCount(); ++i) {
        auto dx = ex.cross(m_molecule.Atom(i).second);
        auto dy = ey.cross(m_molecule.Atom(i).second);
        auto dz = ez.cross(m_molecule.Atom(i).second);

        D(3 * i, 3) = dx(0);
        D(3 * i + 1, 3) = dx(1);
        D(3 * i + 2, 3) = dx(2);

        D(3 * i, 4) = dy(0);
        D(3 * i + 1, 4) = dy(1);
        D(3 * i + 2, 4) = dy(2);

        D(3 * i, 5) = dz(0);
        D(3 * i + 1, 5) = dz(1);
        D(3 * i + 2, 5) = dz(2);
    }

    Eigen::MatrixXd XtX = D.transpose() * D;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(XtX);
    Eigen::MatrixXd S = es.operatorInverseSqrt();
    Matrix R = D * S;
    Matrix f = R.transpose() * hessian * R;
    for (int i = 0; i < 3 * m_molecule.AtomCount(); ++i) {
        for (int j = 0; j < 3 * m_molecule.AtomCount(); ++j) {
            if (i < 6 || j < 6)
                f(i, j) = 0;
        }
    }
    return f;
}

Vector Hessian::ConvertHessian(Matrix& hessian)
{
    Vector vector;

    for (int i = 0; i < m_molecule.AtomCount() * 3; i++) {
        for (int j = 0; j < m_molecule.AtomCount() * 3; j++) {
            double mass = 1 / sqrt(Elements::AtomicMass[m_molecule.Atoms()[i / 3]] * Elements::AtomicMass[m_molecule.Atoms()[j / 3]]);
            hessian(i, j) *= mass;
        }
    }

    Eigen::SelfAdjointEigenSolver<Geometry> diag_I;
    diag_I.compute(hessian);
    return diag_I.eigenvalues();
}

void Hessian::PrintVibrationsPlain(const Vector& eigenvalues)
{
    Vector eigval = eigenvalues.cwiseSqrt();
    std::cout << std::endl
              << " Frequencies: " << std::endl;

    for (int i = 0; i < m_molecule.AtomCount() * 3; ++i) {
        if (i % 6 == 0)
            std::cout << std::endl;
        std::cout << m_scale_functions(eigval(i)) << " ";
    }
    std::cout << std::endl;
}

void Hessian::PrintVibrations(Vector& eigenvalues, const Vector& projected)
{
    Vector eigval = eigenvalues.cwiseSqrt();
    std::cout << std::endl
              << " Frequencies: " << std::endl;

    for (int i = 0; i < m_molecule.AtomCount() * 3; ++i) {
        if (i % 6 == 0)
            std::cout << std::endl;
        if (projected(i) < 0)
            std::cout << m_scale_functions(sqrt(std::abs(eigenvalues(i)))) << "(i) "
                      << " ";
        else {
            if (projected(i) < 1e-10) {
                std::cout << projected(i) << "(*) ";

            } else {
                if (eigenvalues(i) < 0)
                    std::cout << m_scale_functions(sqrt(std::abs(eigenvalues(i)))) << "(*) "
                              << " ";
                else
                    std::cout << m_scale_functions(eigval(i)) << " ";
            }
        }
    }
    std::cout << std::endl;
}

void Hessian::CalculateHessianNumerical()
{
    m_hessian = Eigen::MatrixXd::Ones(3 * m_molecule.AtomCount(), 3 * m_molecule.AtomCount());

    CxxThreadPool* pool = new CxxThreadPool;
    pool->setActiveThreadCount(m_threads);
    if (m_silent)
        pool->setProgressBar(CxxThreadPool::ProgressBarType::None);
    else
        std::cout << "Starting Numerical Hessian Calculation" << std::endl;

    for (int i = 0; i < m_molecule.AtomCount(); ++i) {
        for (int j = 0; j < m_molecule.AtomCount(); ++j) {
            for (int xi = 0; xi < 3; ++xi)
                for (int xj = 0; xj < 3; ++xj) {
                    HessianThread* thread = new HessianThread(m_controller, i, j, xi, xj, true);
                    thread->setMolecule(m_molecule);
                    thread->setParameter(m_parameter);
                    pool->addThread(thread);
                }
        }
    }
    pool->DynamicPool(2);
    pool->StartAndWait();
    for (auto* t : pool->Finished()) {
        HessianThread* thread = static_cast<HessianThread*>(t);
        m_hessian(3 * thread->I() + thread->XI(), 3 * thread->J() + thread->XJ()) *= thread->DD();
    }

    delete pool;
}

void Hessian::CalculateHessianSemiNumerical()
{
    m_hessian = Eigen::MatrixXd::Ones(3 * m_molecule.AtomCount(), 3 * m_molecule.AtomCount());

    CxxThreadPool* pool = new CxxThreadPool;
    pool->setActiveThreadCount(m_threads);
    if (m_silent)
        pool->setProgressBar(CxxThreadPool::ProgressBarType::None);
    else
        std::cout << "Starting Seminumerical Hessian Calculation" << std::endl;

    for (int i = 0; i < m_molecule.AtomCount(); ++i) {
        for (int xi = 0; xi < 3; ++xi) {
            HessianThread* thread = new HessianThread(m_controller, i, 0, xi, 0, false);
            thread->setMolecule(m_molecule);
            thread->setParameter(m_parameter);
            pool->addThread(thread);
        }
    }
    pool->DynamicPool(2);
    pool->StartAndWait();

    for (auto* t : pool->Finished()) {
        HessianThread* thread = static_cast<HessianThread*>(t);
        int i = thread->I();
        int xi = thread->XI();
        Matrix gradient = thread->Gradient();
        for (int j = 0; j < gradient.rows(); ++j) {
            for (int k = 0; k < gradient.cols(); ++k) {
                m_hessian(3 * i + xi, 3 * j + k) = thread->Gradient()(j, k);
            }
        }
    }

    for (int i = 0; i < m_molecule.AtomCount(); ++i) {
        for (int j = 0; j < m_molecule.AtomCount(); ++j) {
            double mass = 1 / sqrt(Elements::AtomicMass[m_molecule.Atoms()[i]] * Elements::AtomicMass[m_molecule.Atoms()[j]]);
            for (int xi = 0; xi < 3; ++xi) {
                for (int xj = 0; xj < 3; ++xj) {
                    double value = (m_hessian(3 * i + xi, 3 * j + xj) + m_hessian(3 * j + xj, 3 * i + xi)) / 2.0;
                    m_hessian(3 * i + xi, 3 * j + xj) = value;
                    m_hessian(3 * j + xj, 3 * i + xi) = value;
                }
            }
        }
    }
    delete pool;
}
