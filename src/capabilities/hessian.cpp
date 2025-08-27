/*
 * < General Calculator for the Hessian Matrix and Vibrational Frequencies>
 * Copyright (C) 2019 - 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 * Claude Generated: Enhanced with configurable finite difference parameters and modern unit system
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
#include <iomanip>
#include <sstream>

#include "external/CxxThreadPool/include/CxxThreadPool.hpp"

#include "src/core/curcuma_logger.h"
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
    if (i == j && xi == i && xj == i && i == 0)
        m_schema = [this]() {
            this->Threaded();
        };
    else if (m_fullnumerical)
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
    // Update finite difference step size from controller
    m_d = Json2KeyWord<double>(m_controller, "finite_diff_step");

    // Numerical stability validation
    if (m_d <= 0.0) {
        CurcumaLogger::error("Finite difference step size must be positive, got: " + std::to_string(m_d));
        m_d = 5e-3; // Fallback to default
    }
    if (m_d < 1e-6) {
        CurcumaLogger::warn("Very small finite difference step (" + std::to_string(m_d) + " Bohr) may cause numerical instability");
    }
    if (m_d > 1e-1) {
        CurcumaLogger::warn("Large finite difference step (" + std::to_string(m_d) + " Bohr) may reduce accuracy");
    }
}

int HessianThread::execute()
{
    m_schema();
    return 0;
}

/**
 * @brief Full numerical Hessian calculation using four-point finite difference
 *
 * Educational Note: This implements the second derivative formula:
 * H_ij = ∂²E/∂x_i∂x_j ≈ [E(x_i+δ,x_j+δ) - E(x_i-δ,x_j+δ) - E(x_i+δ,x_j-δ) + E(x_i-δ,x_j-δ)] / (4δ²)
 *
 * This is the most accurate but computationally expensive method, requiring
 * 4 energy evaluations per Hessian element.
 */
void HessianThread::Numerical()
{
    m_geom_ip_jp = m_molecule.Coords();
    m_geom_im_jp = m_molecule.Coords();
    m_geom_ip_jm = m_molecule.Coords();
    m_geom_im_jm = m_molecule.Coords();

    m_geom_ip_jp(m_i, m_xi) += m_d;
    m_geom_ip_jp(m_j, m_xj) += m_d;

    m_geom_im_jp(m_i, m_xi) -= m_d;
    m_geom_im_jp(m_j, m_xj) += m_d;

    m_geom_ip_jm(m_i, m_xi) += m_d;
    m_geom_ip_jm(m_j, m_xj) -= m_d;

    m_geom_im_jm(m_i, m_xi) -= m_d;
    m_geom_im_jm(m_j, m_xj) -= m_d;

    EnergyCalculator energy(m_method, m_controller);
    energy.setVerbosity(0); // Silent mode for Hessian sub-calculations - set BEFORE initialization
    energy.setParameter(m_parameter);
    energy.setMolecule(m_molecule.getMolInfo());

    const double d2 = 1.0 / (4.0 * m_d * m_d); // Finite difference coefficient
    // Educational Note: Factor 1/(4δ²) comes from the four-point finite difference formula

    energy.updateGeometry(m_geom_ip_jp);
    double energy_ip_jp = energy.CalculateEnergy(false);

    energy.updateGeometry(m_geom_im_jp);
    double energy_im_jp = energy.CalculateEnergy(false);

    energy.updateGeometry(m_geom_ip_jm);
    double energy_ip_jm = energy.CalculateEnergy(false);

    energy.updateGeometry(m_geom_im_jm);
    double energy_im_jm = energy.CalculateEnergy(false);
    // Apply four-point finite difference formula for second derivative
    m_dd = d2 * (energy_ip_jp - energy_im_jp - energy_ip_jm + energy_im_jm);

    // Numerical stability check
    if (std::isnan(m_dd) || std::isinf(m_dd)) {
        CurcumaLogger::error("Numerical instability in Hessian element calculation: NaN or Inf detected");
        m_dd = 0.0; // Set to zero for safety
    }

    // Educational Note: This computes one element H_ij of the Hessian matrix
    // The result is in atomic units (Hartree/Bohr²)
}
/**
 * @brief Seminumerical Hessian calculation using gradient finite differences
 *
 * Educational Note: This implements the mixed derivative formula:
 * H_ij = ∂²E/∂x_i∂x_j ≈ [∇E(x_i+δ) - ∇E(x_i-δ)] / (2δ)
 *
 * More efficient than full numerical (2 gradient vs 4 energy evaluations per element)
 * but requires the energy method to provide accurate gradients.
 */
void HessianThread::Seminumerical()
{
    m_geom_ip_jp = m_molecule.Coords();
    m_geom_im_jp = m_molecule.Coords();

    m_geom_ip_jp(m_i, m_xi) += m_d;
    EnergyCalculator energy(m_method, m_controller);
    energy.setVerbosity(0); // Silent mode for Hessian sub-calculations - set BEFORE initialization
    energy.setParameter(m_parameter);
    energy.setMolecule(m_molecule.getMolInfo());

    energy.updateGeometry(m_geom_ip_jp);
    energy.CalculateEnergy(true);
    Matrix gradientp = energy.Gradient();

    m_geom_im_jp(m_i, m_xi) -= m_d;

    energy.updateGeometry(m_geom_im_jp);
    energy.CalculateEnergy(true);
    Matrix gradientm = energy.Gradient();
    m_gradient = (gradientp - gradientm) / (2.0 * m_d); // Central difference for gradient
}

void HessianThread::Threaded()
{
    m_hessian = Eigen::MatrixXd::Zero(3 * m_molecule.AtomCount(), 3 * m_molecule.AtomCount());

    EnergyCalculator energy(m_method, m_controller);
    energy.setVerbosity(0); // Silent mode for Hessian sub-calculations - set BEFORE initialization
    energy.setParameter(m_parameter);
    energy.setMolecule(m_molecule.getMolInfo());
    for (int i : m_atoms) {
        for (int xi = 0; xi < 3; ++xi) {

            m_geom_ip_jp = m_molecule.Coords();
            m_geom_im_jp = m_molecule.Coords();
            m_geom_ip_jp(i, xi) += m_d;

            energy.updateGeometry(m_geom_ip_jp);
            energy.CalculateEnergy(true);
            Matrix gradientp = energy.Gradient();

            m_geom_im_jp(i, xi) -= m_d;

            energy.updateGeometry(m_geom_im_jp);
            energy.CalculateEnergy(true);
            Matrix gradientm = energy.Gradient();
            m_gradient = (gradientp - gradientm) / (2.0 * m_d); // Central difference for gradient

            for (int j = 0; j < m_gradient.rows(); ++j) {
                for (int k = 0; k < m_gradient.cols(); ++k) {
                    m_hessian(3 * i + xi, 3 * j + k) = m_gradient(j, k);
                }
            }
        }
    }
}

Hessian::Hessian(const std::string& method, const json& controller, bool silent)
    : CurcumaMethod(HessianJson, controller, silent)
{
    UpdateController(controller);
    m_controller = MergeJson(m_defaults, controller);
    m_threads = m_controller["threads"];
    // Frequency scaling: Convert from √(Hartree/atomic_mass_unit) to cm⁻¹
    // Reference: Computational chemistry frequency scaling factors
    m_scale_functions = [](double eigenvalue_sqrt) -> double {
        // Convert from atomic units to wavenumbers (cm⁻¹)
        // Factor: 5140.4 cm⁻¹ * √(Eh/u) + empirical correction
        return eigenvalue_sqrt * CurcumaUnit::Energy::HARTREE_TO_WAVENUMBER / std::sqrt(CurcumaUnit::Constants::AMU_TO_AU) + 47.349;
    };
    m_method = method;
}

Hessian::Hessian(const json& controller, bool silent)
    : CurcumaMethod(HessianJson, controller, silent)
{
    UpdateController(controller);
    m_controller = MergeJson(m_defaults, controller);

    m_threads = m_controller["threads"];

    // Frequency scaling: Convert from √(Hartree/atomic_mass_unit) to cm⁻¹
    // Reference: Computational chemistry frequency scaling factors
    m_scale_functions = [](double eigenvalue_sqrt) -> double {
        // Convert from atomic units to wavenumbers (cm⁻¹)
        // Factor: 5140.4 cm⁻¹ * √(Eh/u) + empirical correction
        return eigenvalue_sqrt * CurcumaUnit::Energy::HARTREE_TO_WAVENUMBER / std::sqrt(CurcumaUnit::Constants::AMU_TO_AU) + 47.349;
    };
}
Hessian::~Hessian()
{
}
void Hessian::LoadControlJson()
{
    m_controller = MergeJson(m_defaults, m_controller);
    m_hess_calc = Json2KeyWord<bool>(m_controller, "hess_calc");
    m_write_file = Json2KeyWord<std::string>(m_controller, "hess_write_file");
    m_hess_read = Json2KeyWord<bool>(m_controller, "hess_read");
    m_read_file = Json2KeyWord<std::string>(m_controller, "hess_read_file");
    m_read_xyz = Json2KeyWord<std::string>(m_controller, "hess_read_xyz");

    m_write_file = Json2KeyWord<std::string>(m_controller, "hess_write_file");
    if (m_hess_read) {
        m_hess_calc = false;
    }

    m_freq_scale = Json2KeyWord<double>(m_controller, "freq_scale");
    m_thermo = Json2KeyWord<double>(m_controller, "thermo");
    m_freq_cutoff = Json2KeyWord<double>(m_controller, "freq_cutoff");
    m_finite_diff_step = Json2KeyWord<double>(m_controller, "finite_diff_step");
    m_verbosity = Json2KeyWord<int>(m_controller, "verbosity");
    m_hess = Json2KeyWord<int>(m_controller, "hess");
    m_method = Json2KeyWord<std::string>(m_controller, "method");
    m_threads = m_controller["threads"];

    // Initialize CurcumaLogger with verbosity level
    CurcumaLogger::set_verbosity(m_verbosity);
}

void Hessian::setMolecule(const Molecule& molecule)
{
    m_molecule = molecule;
    m_atoms_j.resize(m_molecule.AtomCount());
    m_atoms_i.resize(m_molecule.AtomCount());

    for (int i = 0; i < m_molecule.AtomCount(); ++i) {
        m_atoms_i[i] = i;
        m_atoms_j[i] = i;
    }
}

void Hessian::LoadMolecule(const std::string& file)
{
    m_molecule = Files::LoadFile(file);
    m_atom_count = m_molecule.AtomCount();
    m_atoms_j.resize(m_atom_count);
    m_atoms_i.resize(m_atom_count);

    for (int i = 0; i < m_atom_count; ++i) {
        m_atoms_i[i] = i;
        m_atoms_j[i] = i;
    }
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
                m_hessian(row, column) = std::stod(entries[i]) / CurcumaUnit::Legacy::au / CurcumaUnit::Legacy::au; // Convert from Å⁻² to Bohr⁻²
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
            CalculateHessianThreaded();
        } else if (m_hess == 2) {
            CalculateHessianNumerical();
        } else {
            CalculateHessianSemiNumerical();
        }
    } else {
        LoadMolecule(m_read_xyz);
        LoadHessian(m_read_file);
    }

    m_frequencies = ConvertHessian(m_hessian);

    // PrintVibrationsPlain(eigenvalues);

    auto hessian2 = ProjectHessian(m_hessian);
    auto projected = ConvertHessian(hessian2);
    if (m_verbosity >= 1)
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

/**
 * @brief Apply mass-weighting and diagonalize Hessian to get frequencies
 *
 * Educational Background:
 * The Hessian matrix H contains force constants in Cartesian coordinates.
 * To obtain vibrational frequencies, we need to:
 *
 * 1. Mass-weight the Hessian: H'_ij = H_ij / √(m_i * m_j)
 *    This accounts for the kinetic energy in the vibrational equation of motion
 *
 * 2. Diagonalize: H'φ_k = λ_kφ_k
 *    Eigenvalues λ_k are frequency², eigenvectors φ_k are normal modes
 *
 * 3. Convert to frequencies: ω_k = √|λ_k| (in atomic units)
 *    Negative eigenvalues indicate imaginary frequencies (unstable geometry)
 *
 * Physical interpretation:
 * - Eigenvalues > 0: Real vibrational frequencies
 * - Eigenvalues ≈ 0: Translation/rotation modes (should be ~6 for non-linear molecules)
 * - Eigenvalues < 0: Imaginary frequencies (saddle point, unstable geometry)
 *
 * @param hessian Hessian matrix to be mass-weighted and diagonalized
 * @return Vector of eigenvalues (frequency² in atomic units)
 */
Vector Hessian::ConvertHessian(Matrix& hessian)
{
    Vector vector;

    // Apply mass-weighting: H'_ij = H_ij / sqrt(m_i * m_j)
    for (int i = 0; i < m_molecule.AtomCount() * 3; i++) {
        for (int j = 0; j < m_molecule.AtomCount() * 3; j++) {
            double mass_factor = 1.0 / sqrt(Elements::AtomicMass[m_molecule.Atoms()[i / 3]] * Elements::AtomicMass[m_molecule.Atoms()[j / 3]]);
            hessian(i, j) *= mass_factor;
        }
    }

    // Numerical stability checks before diagonalization
    if (m_verbosity >= 3) {
        // Check for NaN or Inf in Hessian matrix
        bool hasNaN = false, hasInf = false;
        for (int i = 0; i < hessian.rows() && !hasNaN && !hasInf; ++i) {
            for (int j = 0; j < hessian.cols() && !hasNaN && !hasInf; ++j) {
                if (std::isnan(hessian(i, j)))
                    hasNaN = true;
                if (std::isinf(hessian(i, j)))
                    hasInf = true;
            }
        }
        if (hasNaN)
            CurcumaLogger::error("NaN values detected in Hessian matrix!");
        if (hasInf)
            CurcumaLogger::error("Infinite values detected in Hessian matrix!");
    }

    // Diagonalize mass-weighted Hessian to obtain eigenvalues (ω²) and eigenvectors (normal modes)
    Eigen::SelfAdjointEigenSolver<Geometry> diag_I;
    diag_I.compute(hessian);

    // Check diagonalization success
    if (diag_I.info() != Eigen::Success) {
        CurcumaLogger::error("Hessian diagonalization failed! Check matrix conditioning.");
        return Vector::Zero(hessian.rows());
    }

    // Educational Note: Eigenvalues are ω² in atomic units
    // To get frequencies in cm⁻¹: ν = (1/2π) * √(λ) * conversion_factor
    return diag_I.eigenvalues();
}

void Hessian::PrintVibrationsPlain(const Vector& eigenvalues)
{
    Vector eigval = eigenvalues.cwiseSqrt();
    CurcumaLogger::success("\nVibrational Frequencies (cm⁻¹):");

    std::ostringstream freq_output;
    for (int i = 0; i < m_molecule.AtomCount() * 3; ++i) {
        if (i % 6 == 0 && i > 0)
            freq_output << "\n";
        freq_output << std::fixed << std::setprecision(1) << m_scale_functions(eigval(i)) << " ";
    }
    CurcumaLogger::success(freq_output.str());
}

void Hessian::PrintVibrations(Vector& eigenvalues, const Vector& projected)
{
    Vector eigval = eigenvalues.cwiseSqrt();
    CurcumaLogger::success("\nVibrational Frequencies (cm⁻¹):");

    if (m_verbosity >= 2) {
        CurcumaLogger::info("Legend: (i)=imaginary, (*)=negative/zero eigenvalue, unmarked=normal mode");
    }

    std::ostringstream freq_output;
    int real_modes = 0, imaginary_modes = 0, zero_modes = 0;

    for (int i = 0; i < m_molecule.AtomCount() * 3; ++i) {
        if (i % 6 == 0 && i > 0)
            freq_output << "\n";

        if (projected(i) < 0) {
            freq_output << std::fixed << std::setprecision(1) << m_scale_functions(sqrt(std::abs(eigenvalues(i)))) << "(i) ";
            imaginary_modes++;
        } else {
            if (projected(i) < 1e-10) {
                freq_output << std::fixed << std::setprecision(6) << projected(i) << "(*) ";
                zero_modes++;
            } else {
                if (eigenvalues(i) < 0) {
                    freq_output << std::fixed << std::setprecision(1) << m_scale_functions(sqrt(std::abs(eigenvalues(i)))) << "(*) ";
                } else {
                    freq_output << std::fixed << std::setprecision(1) << m_scale_functions(eigval(i)) << " ";
                    real_modes++;
                }
            }
        }
    }
    CurcumaLogger::success(freq_output.str());

    if (m_verbosity >= 2) {
        CurcumaLogger::param("Real vibrational modes", std::to_string(real_modes));
        CurcumaLogger::param("Zero/translational modes", std::to_string(zero_modes));
        if (imaginary_modes > 0)
            CurcumaLogger::warn("Imaginary frequencies found: " + std::to_string(imaginary_modes) + " (indicates saddle point or incorrect geometry)");
    }
}

void Hessian::CalculateHessianThreaded()
{
    m_hessian = Eigen::MatrixXd::Zero(3 * m_molecule.AtomCount(), 3 * m_molecule.AtomCount());

    CxxThreadPool* pool = new CxxThreadPool;
    pool->setActiveThreadCount(m_threads);
    if (m_verbosity == 0) {
        pool->setProgressBar(CxxThreadPool::ProgressBarType::None);
    } else {
        pool->setProgressBar(CxxThreadPool::ProgressBarType::Continously); // Always disable for cleaner output
        CurcumaLogger::info("Starting Hessian calculation using method: " + m_method);
        CurcumaLogger::param("Calculation method", (m_hess == 1 ? "Threaded seminumerical" : (m_hess == 2 ? "Full numerical" : "Seminumerical")));
        CurcumaLogger::param("Number of threads", std::to_string(m_threads));
        CurcumaLogger::param("Finite difference step", std::to_string(m_finite_diff_step) + " Bohr (" + std::to_string(m_finite_diff_step * CurcumaUnit::Length::BOHR_TO_ANGSTROM) + " Å)");
    }

    std::vector<int> atoms;
    if (m_threads > m_molecule.AtomCount())
        m_threads = m_molecule.AtomCount();
    if (m_method.compare("gfnff") == 0) {
        m_threads = 1;
        CurcumaLogger::warn("GFN-FF enforces single thread approach for numerical stability");
    }
    std::vector<std::vector<int>> threads(m_threads);

    for (int i = 0; i < m_molecule.AtomCount(); ++i) {
        atoms.push_back(i);
        threads[i % m_threads].push_back(i);
    }

    for (int i = 0; i < threads.size(); ++i) {
        HessianThread* thread = new HessianThread(m_controller, 0, 0, 0, 0, true);
        thread->setMethod(m_method);
        thread->setMolecule(m_molecule);
        thread->setParameter(m_parameter);
        thread->setIndices(threads[i]);
        pool->addThread(thread);
    }

    pool->StaticPool();
    pool->StartAndWait();
    int i = 0;
    for (auto* t : pool->getFinishedThreads()) {
        HessianThread* thread = static_cast<HessianThread*>(t);
        m_hessian += thread->getHessian();
        ++i;
    }

    delete pool;

    // Symmetrize Hessian matrix and check for numerical errors (threaded method)
    double max_asymmetry = 0.0;
    int asymmetric_elements = 0;

    for (int i = 0; i < m_molecule.AtomCount(); ++i) {
        for (int j = 0; j < m_molecule.AtomCount(); ++j) {
            for (int xi = 0; xi < 3; ++xi) {
                for (int xj = 0; xj < 3; ++xj) {
                    int idx_ij = 3 * i + xi;
                    int idx_ji = 3 * j + xj;

                    double h_ij = m_hessian(idx_ij, idx_ji);
                    double h_ji = m_hessian(idx_ji, idx_ij);

                    // Check symmetry before averaging
                    double asymmetry = std::abs(h_ij - h_ji);
                    if (asymmetry > max_asymmetry)
                        max_asymmetry = asymmetry;
                    if (asymmetry > 1e-8)
                        asymmetric_elements++;

                    // Symmetrize by averaging
                    double symmetric_value = (h_ij + h_ji) / 2.0;
                    m_hessian(idx_ij, idx_ji) = symmetric_value;
                    m_hessian(idx_ji, idx_ij) = symmetric_value;
                }
            }
        }
    }

    // Report symmetry analysis
    if (m_verbosity >= 2) {
        CurcumaLogger::param("Maximum asymmetry", std::to_string(max_asymmetry) + " Eh/Bohr²");
        if (asymmetric_elements > 0) {
            CurcumaLogger::param("Asymmetric elements", std::to_string(asymmetric_elements));
        }
        if (max_asymmetry > 1e-6) {
            CurcumaLogger::warn("Significant Hessian asymmetry detected - may indicate numerical issues");
        }
    }
}

void Hessian::CalculateHessianNumerical()
{
    m_hessian = Eigen::MatrixXd::Ones(3 * m_molecule.AtomCount(), 3 * m_molecule.AtomCount());

    CxxThreadPool* pool = new CxxThreadPool;
    pool->setActiveThreadCount(m_threads);
    if (m_verbosity == 0) {
        pool->setProgressBar(CxxThreadPool::ProgressBarType::None);
    } else {
        pool->setProgressBar(CxxThreadPool::ProgressBarType::Continously);
        CurcumaLogger::info("Starting full numerical Hessian calculation");
        CurcumaLogger::param("Total energy evaluations", std::to_string(m_molecule.AtomCount() * m_molecule.AtomCount() * 9));
        CurcumaLogger::param("Finite difference step", std::to_string(m_finite_diff_step) + " Bohr");
    }

    for (/*const auto i : m_atoms_i */ int i = 0; i < m_molecule.AtomCount(); ++i) {
        for (/*const auto j : m_atoms_j */ int j = 0; j < m_molecule.AtomCount(); ++j) {
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
    for (auto* t : pool->getFinishedThreads()) {
        HessianThread* thread = static_cast<HessianThread*>(t);
        m_hessian(3 * thread->I() + thread->XI(), 3 * thread->J() + thread->XJ()) *= thread->DD();
    }

    delete pool;
}

void Hessian::CalculateHessianSemiNumerical()
{
    m_hessian = Eigen::MatrixXd::Zero(3 * m_molecule.AtomCount(), 3 * m_molecule.AtomCount());

    CxxThreadPool* pool = new CxxThreadPool;
    pool->setActiveThreadCount(m_threads);
    if (m_verbosity == 0) {
        pool->setProgressBar(CxxThreadPool::ProgressBarType::None);
    } else {
        pool->setProgressBar(CxxThreadPool::ProgressBarType::Continously);
        CurcumaLogger::info("Starting seminumerical Hessian calculation");
        CurcumaLogger::param("Total gradient evaluations", std::to_string(m_molecule.AtomCount() * 3 * 2));
        CurcumaLogger::param("Finite difference step", std::to_string(m_finite_diff_step) + " Bohr");
    }

    for (/*auto int i : m_atoms_i */ int i = 0; i < m_molecule.AtomCount(); ++i) {
        for (int xi = 0; xi < 3; ++xi) {
            HessianThread* thread = new HessianThread(m_controller, i, 0, xi, 0, false);
            thread->setMolecule(m_molecule);
            thread->setParameter(m_parameter);
            pool->addThread(thread);
        }
    }
    pool->DynamicPool(2);
    pool->StartAndWait();

    for (auto* t : pool->getFinishedThreads()) {
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

    // Symmetrize Hessian matrix (seminumerical method)
    double max_asymmetry = 0.0;

    for (int i = 0; i < m_molecule.AtomCount(); ++i) {
        for (int j = 0; j < m_molecule.AtomCount(); ++j) {
            for (int xi = 0; xi < 3; ++xi) {
                for (int xj = 0; xj < 3; ++xj) {
                    int idx_ij = 3 * i + xi;
                    int idx_ji = 3 * j + xj;

                    double h_ij = m_hessian(idx_ij, idx_ji);
                    double h_ji = m_hessian(idx_ji, idx_ij);

                    // Check symmetry
                    double asymmetry = std::abs(h_ij - h_ji);
                    if (asymmetry > max_asymmetry)
                        max_asymmetry = asymmetry;

                    // Symmetrize by averaging
                    double symmetric_value = (h_ij + h_ji) / 2.0;
                    m_hessian(idx_ij, idx_ji) = symmetric_value;
                    m_hessian(idx_ji, idx_ij) = symmetric_value;
                }
            }
        }
    }

    // Report symmetry for seminumerical method
    if (m_verbosity >= 2) {
        CurcumaLogger::param("Matrix symmetry check", "Max asymmetry: " + std::to_string(max_asymmetry) + " Eh/Bohr²");
    }
    delete pool;
}
