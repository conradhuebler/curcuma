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

// #include <catch2/catch_test_macros.hpp>
// #include <catch2/generators/catch_generators_all.hpp>

// #include <finitediff.hpp>
// #include <spdlog/spdlog.h>

#include "external/CxxThreadPool/include/CxxThreadPool.h"

#include "src/core/energycalculator.h"

#include "hessian.h"

// using namespace fd;

HessianThread::HessianThread(const std::string& method, const json& controller, int i, int j, int xi, int xj, bool fullnumerical)
    : m_method(method)
    , m_controller(controller)
    , m_i(i)
    , m_j(j)
    , m_xi(xi)
    , m_xj(xj)
    , m_fullnumerical(fullnumerical)
{
    setAutoDelete(true);

    if (m_fullnumerical)
        m_schema = [this]() {
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
    m_d = 5e-3;

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
    m_d = 5e-3;
    m_geom_ip_jp = m_molecule.Coords();
    m_geom_im_jp = m_molecule.Coords();

    m_geom_ip_jp[m_i][m_xi] += m_d;

    EnergyCalculator energy(m_method, m_controller);
    energy.setMolecule(m_molecule);

    energy.updateGeometry(m_geom_ip_jp);
    double energy_ip_jp = energy.CalculateEnergy(true, false);
    Matrix gradientp = energy.Gradient();

    m_geom_im_jp[m_i][m_xi] -= m_d;

    energy.updateGeometry(m_geom_im_jp);
    double energy_im_jp = energy.CalculateEnergy(true, false);
    Matrix gradientm = energy.Gradient();
    // std::cout << gradientp << std::endl << std::endl << gradientm << std::endl;
    m_gradient = (gradientp - gradientm) / (2 * m_d) / au / au;
    // std::cout << m_gradient << std::endl << std::endl;
}

Hessian::Hessian(const std::string& method, const json& controller, int threads, bool silent)
    : CurcumaMethod(HessianJson, controller, silent)
    , m_method(method)
    , m_controller(controller)
    , m_threads(threads)
{
    UpdateController(controller);
}

Hessian::Hessian(const json& controller, bool silent)
    : CurcumaMethod(HessianJson, controller, silent)
    , m_controller(controller)
{
    UpdateController(controller);
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
    auto freqs = ConvertHessian(m_hessian);
}

void Hessian::CalculateHessian(bool fullnumerical)
{
    start();
    //   m_hessian = Eigen::MatrixXd::Ones(3 * m_molecule.AtomCount(), 3 * m_molecule.AtomCount());
    /*
        if (fullnumerical)
            CalculateHessianNumerical();
        else
            CalculateHessianSemiNumerical();
    */

    /*

        CalculateHessianNumerical();
        auto freqs = ConvertHessian(m_hessian);
        std::cout << " Frequencies " << std::endl << std::endl;
        std::cout << freqs.transpose() << std::endl;
        std::cout <<std::endl << std::endl;
    */
    // CalculateHessianSemiNumerical();
    // auto freqs = ConvertHessian(m_hessian);
    //  std::cout << freqs.transpose() << std::endl;
    //  std::cout <<std::endl << std::endl;

    /*

        FiniteDiffHess();
        freqs = ConvertHessian(m_hessian);
        std::cout << " Frequencies " << std::endl << std::endl;
        std::cout << freqs.transpose() << std::endl;
        std::cout <<std::endl << std::endl;
    */
}

Vector Hessian::ConvertHessian(Matrix& hessian)
{
    // std::cout << hessian << std::endl;
    {
        Eigen::SelfAdjointEigenSolver<Geometry> diag_I;
        diag_I.compute(hessian);
        // std::cout << diag_I.eigenvalues().cwiseSqrt() << std::endl;
    }
    Vector vector;
    constexpr double Eh_to_cm = 27.211386 * 8065.54;
    double Eh_to_Ha = 1.0 / 27.211386; // convert from Eh to Ha
    double Ha_to_eV = 27.211386; // convert from Ha to eV
    double amu_to_kg = 1.66054E-27; // convert from amu to kg
    double Angstrom_to_m = 1E-10; // convert from Angstrom to meters
    double c = 2.99792E10; // speed of light in cm/s
    double h = 6.62607E-34; // Planck's constant in J*s
    double eV_to_cm_inv = 8065.54; // convert from eV to cm^-1
    const double conv = 1.0e12; // conversion factor from Angstroms to cm
    auto amu = m_molecule.Atoms();

    for (int i = 0; i < m_molecule.AtomCount() * 3; i++) {
        for (int j = 0; j < m_molecule.AtomCount() * 3; j++) {
            double mass = 1 / sqrt(Elements::AtomicMass[m_molecule.Atoms()[i / 3]] * Elements::AtomicMass[m_molecule.Atoms()[j / 3]]);

            hessian(i, j) *= Ha_to_eV; // convert from Eh/A^2 to eV/A^2
            hessian(i, j) /= (Angstrom_to_m * Angstrom_to_m); // from A^2 to m^2
            hessian(i, j) *= mass * amu_to_kg; // sqrt(Elements::AtomicMass[amu[i/3]]*Elements::AtomicMass[amu[j/3]]); // convert from Ha/m^2 to Ha/(amu*m^2)
        }
    }

    Eigen::SelfAdjointEigenSolver<Geometry> diag_I;
    diag_I.compute(hessian);

    vector = diag_I.eigenvalues().cwiseSqrt();
    // std::cout << vector.transpose() << std::endl;
    /* this conversation factor has to be made pure unit based ... */
    vector *= 1.2796e+06; //(2.0 * M_PI * c) * eV_to_cm_inv*conv; // / sqrt(Elements::AtomicMass[amu[i/3]]);
    // std::cout << eV_to_cm_inv*conv/(2.0 * M_PI * c) << std::endl;
    // vector *= Ha_to_eV*eV_to_cm_inv*conv/(2.0 * M_PI * c);
    // std::cout << vector.transpose() << std::endl;

    std::cout << std::endl
              << " Frequencies: " << std::endl;

    for (int i = 0; i < m_molecule.AtomCount() * 3; ++i) {
        if (i % 6 == 0)
            std::cout << std::endl;
        if (diag_I.eigenvalues()(i) < 0)
            std::cout << sqrt(std::abs(diag_I.eigenvalues()(i))) * 1.2796e+06 << "(i)"
                      << " ";
        else
            std::cout << vector(i) << " ";
    }
    std::cout << std::endl;
    return vector;
}

void Hessian::CalculateHessianNumerical()
{
    m_hessian = Eigen::MatrixXd::Ones(3 * m_molecule.AtomCount(), 3 * m_molecule.AtomCount());

    std::cout << "Starting Numerical Hessian Calculation" << std::endl;

    CxxThreadPool* pool = new CxxThreadPool;
    pool->setActiveThreadCount(m_threads);
    for (int i = 0; i < m_molecule.AtomCount(); ++i) {
        for (int j = 0; j < m_molecule.AtomCount(); ++j) {
            for (int xi = 0; xi < 3; ++xi)
                for (int xj = 0; xj < 3; ++xj) {
                    HessianThread* thread = new HessianThread(m_method, m_controller, i, j, xi, xj, true);
                    thread->setMolecule(m_molecule);
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

    std::cout << "Starting Seminumerical Hessian Calculation" << std::endl;
    CxxThreadPool* pool = new CxxThreadPool;
    pool->setActiveThreadCount(m_threads);

    for (int i = 0; i < m_molecule.AtomCount(); ++i) {
        for (int xi = 0; xi < 3; ++xi) {
            HessianThread* thread = new HessianThread(m_method, m_controller, i, 0, xi, 0, false);
            thread->setMolecule(m_molecule);
            pool->addThread(thread);
        }
    }
    pool->DynamicPool(2);
    pool->StartAndWait();
    // std::cout <<std::endl << m_hessian << std::endl << std::endl;

    for (auto* t : pool->Finished()) {
        HessianThread* thread = static_cast<HessianThread*>(t);
        int i = thread->I();
        int xi = thread->XI();
        // std::cout << i << " " << xi << std::endl;
        // std::cout << thread->Gradient() << std::endl;
        Matrix gradient = thread->Gradient();
        for (int j = 0; j < gradient.rows(); ++j) {
            for (int k = 0; k < gradient.cols(); ++k) {
                //         std::cout << thread->Gradient()(j, k) <<  " ";
                m_hessian(3 * i + xi, 3 * j + k) = thread->Gradient()(j, k);
            }
        }
        //  std::cout <<std::endl << m_hessian << std::endl << std::endl;
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
    /*
    std::cout << std::endl
              << std::endl;
    std::cout << m_hessian << std::endl;
    */
    delete pool;
}

void Hessian::FiniteDiffHess()
{
    /*
    m_hessian = Eigen::MatrixXd::Ones(3 * m_molecule.AtomCount(), 3 * m_molecule.AtomCount());

      EnergyCalculator energy(m_method, m_controller);
      energy.setMolecule(m_molecule);

      const auto f = [&energy](const Eigen::VectorXd& x) -> double {
          energy.updateGeometry(x);
          return energy.CalculateEnergy(false, false);
      };

      Eigen::VectorXd x = Eigen::VectorXd::Zero(3*m_molecule.AtomCount());
    for(int i = 0; i < m_molecule.AtomCount(); ++i)
    {
        x(3*i + 0) = m_molecule.Atom(i).second[0];
        x(3*i + 1) = m_molecule.Atom(i).second[1];
        x(3*i + 2) = m_molecule.Atom(i).second[2];

    }
      Eigen::MatrixXd hess = Eigen::MatrixXd::Zero(3*m_molecule.AtomCount(), 3*m_molecule.AtomCount());

      Eigen::MatrixXd fhess;
      finite_hessian(x, f, fhess, EIGHTH);
      m_hessian = fhess;
      */
}
