/*
 * <Native LBFGS Optimizer Implementation - Claude 3.5 Generated, Cleaned by Claude 4>
 * Copyright (C) 2025 Claude AI - Generated Code
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
 */

#include "lbfgs.h"
#include <cmath>
#include <iomanip>
#include <iostream>

// Claude Generated - Constructor with modern defaults
LBFGS::LBFGS(int memory_size)
    : m_memory_size(memory_size)
    , stepCount(0)
    , m_energy(0.0)
    , m_step_size(1.0) // More conservative default
{
    if (memory_size <= 0) {
        CurcumaLogger::warn("Invalid LBFGS memory size, using default value 10");
        m_memory_size = 10;
    }
}

// Claude Generated - Modern initialization with logging
void LBFGS::initialize(int natoms, const Vector& initial_coordinates)
{
    m_atoms = natoms;
    const int dof = 3 * m_atoms; // Degrees of freedom

    // Initialize state vectors
    m_gradient = Vector::Zero(dof);
    m_constraints = Vector::Ones(dof); // All coordinates free by default
    m_current_coordinates = initial_coordinates;
    stepCount = 0;
    m_error = false;

    // Initialize Hessian approximations
    m_hessian = Matrix::Identity(dof, dof);
    m_inverse_hessian = Matrix::Identity(dof, dof);

    // Clear optimization history
    m_step_history.clear();
    m_gradient_history.clear();
    m_rho_history.clear();
    m_diis_solutions.clear();
    m_diis_errors.clear();

    // Logging - Claude Generated
    if (m_verbosity >= 1) {
        CurcumaLogger::info_fmt("Native LBFGS optimizer initialized:");
        CurcumaLogger::param("Atoms", m_atoms);
        CurcumaLogger::param("Degrees of freedom", dof);
        CurcumaLogger::param("Memory size", m_memory_size);
        CurcumaLogger::param("Method", static_cast<int>(m_method));

        if (m_method == Method::DIIS) {
            CurcumaLogger::param("DIIS history", m_diis_hist);
            CurcumaLogger::param("DIIS start", m_diis_start);
        }

        if (m_method == Method::RFO) {
            CurcumaLogger::param("RFO lambda", m_lambda);
        }
    }
}

// Claude Generated - Modern step method with method switching logic
Vector LBFGS::step()
{
    Vector new_coordinates = m_current_coordinates;
    m_last_step_size = m_step_size;

    // Method selection with fallback logic
    std::string method_name;
    if (m_method == Method::LBFGS || (m_method == Method::DIIS && stepCount < m_diis_start)) {
        method_name = "L-BFGS";
        new_coordinates = LBFGSStep();
    } else if (m_method == Method::DIIS) {
        method_name = "DIIS";
        new_coordinates = DIISStep();
    } else if (m_method == Method::RFO) {
        method_name = "RFO";
        new_coordinates = RFOStep();
    }

    // Log method switching if applicable - Claude Generated
    if (m_method == Method::DIIS && stepCount == m_diis_start && m_verbosity >= 2) {
        CurcumaLogger::info_fmt("Switching from {} to {} after {} iterations", "L-BFGS", "DIIS", stepCount);
    }

    stepCount++;

    // Clean up DIIS history if needed
    if (m_diis_solutions.size() >= static_cast<size_t>(m_diis_hist)) {
        m_diis_solutions.erase(m_diis_solutions.begin());
        m_diis_errors.erase(m_diis_errors.begin());
    }

    return new_coordinates;
}
Eigen::MatrixXd LBFGS::update_inverse_hessian_approximation(Eigen::MatrixXd& B, const Eigen::VectorXd& deltaX, const Eigen::VectorXd& deltaG)
{
    double rho = 1.0 / deltaG.dot(deltaX);
    Eigen::VectorXd B_deltaX = B * deltaX;

    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(B.rows(), B.cols());

    B = (I - rho * deltaX * deltaG.transpose()) * B * (I - rho * deltaG * deltaX.transpose()) + rho * (deltaX * deltaX.transpose());

    return B;
}

void LBFGS::setHessian(const Matrix& hessian)
{
    m_hessian = hessian;

    Eigen::SelfAdjointEigenSolver<Matrix> solver(m_hessian);
    if (solver.info() != Eigen::Success) {
        std::cerr << "Eigenvalue decomposition failed!" << std::endl;
    }

    m_eigenvalues = solver.eigenvalues();
    m_eigenvectors = solver.eigenvectors();
    m_inverse_hessian = m_hessian.inverse();
}

double LBFGS::lineSearchRFO(const Vector& x, const Vector& p, double c1, double c2)
{
    double alpha = m_step_size;
    double alpha_prev = m_last_step_size;

    double f_x = m_energy;
    double gdotp = m_gradient.dot(p);
    int maxiter = 100;
    int iter = 0;
    while (std::abs(alpha - alpha_prev) > 1e-2) {
        Vector x_new = x + alpha * p;

        double f_x_new = getEnergy(x_new);

        if (f_x_new > f_x + c1 * alpha * gdotp || (alpha > 1e-10 && f_x_new >= getEnergy(x + alpha_prev * p))) {
            double alpha_temp = alpha;
            alpha = (alpha_prev + alpha) / 2.0;
            alpha_prev = alpha_temp;
        } else {
            getEnergyGradient(x_new);
            Vector grad_new = m_gradient;

            if (std::abs(grad_new.dot(p)) <= -c2 * gdotp) {
                return alpha;
            } else if (grad_new.dot(p) >= 0) {
                double alpha_temp = alpha;
                alpha = (alpha_prev + alpha) / 2.0;
                alpha_prev = alpha_temp;
            } else {
                alpha_prev = alpha;
                alpha *= 2.0;
            }
        }
        iter++;
    }
    return alpha;
}

double LBFGS::line_search_backtracking(const Vector& x, const Vector& p, double alpha_init, double rho, double c1)
{
    double alpha = alpha_init;
    double f_x = getEnergy(x);
    double gradient_dot_p = m_gradient.dot(p);

    while (getEnergy(x + alpha * p) > f_x + c1 * alpha * gradient_dot_p) {
        alpha *= rho;
    }
    return alpha;
}

Vector LBFGS::diisExtrapolation()
{
    int numErrors = m_diis_errors.size();
    Matrix B = Matrix::Zero(numErrors + 1, numErrors + 1);
    Matrix A = Matrix::Zero(numErrors, numErrors);
    Vector rhs = Vector::Zero(numErrors + 1);
    rhs(numErrors) = 1.0;

    for (int i = 0; i < numErrors; ++i) {
        for (int j = 0; j < numErrors; ++j) {
            B(i, j) = (m_inverse_hessian * m_diis_errors[i]).dot(m_inverse_hessian * m_diis_errors[j]);
            A(i, j) = B(i, j);
        }
        B(numErrors, i) = 1.0;
        B(i, numErrors) = 1.0;
    }
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A);
    double conditionNumber = svd.singularValues()(0) / svd.singularValues().tail(1)(0);

    if (conditionNumber > 1e10) {
        m_diis_solutions.erase(m_diis_solutions.begin());
        m_diis_errors.erase(m_diis_errors.begin());
        std::cout << m_diis_errors.size() << " " << conditionNumber << std::endl;
        int numErrors = m_diis_errors.size();
        Matrix B = Matrix::Zero(numErrors + 1, numErrors + 1);
        Vector rhs = Vector::Zero(numErrors + 1);
        rhs(numErrors) = -1.0;

        for (int i = 0; i < numErrors; ++i) {
            for (int j = 0; j < numErrors; ++j) {
                B(i, j) = (m_inverse_hessian * m_diis_errors[i]).dot(m_inverse_hessian * m_diis_errors[j]);
            }
            B(numErrors, i) = -1.0;
            B(i, numErrors) = -1.0;
        }
        Vector coeff = B.colPivHouseholderQr().solve(rhs);

        Vector diisSolution = Vector::Zero(m_diis_solutions[0].size());
        for (int i = 0; i < numErrors; ++i) {
            diisSolution += coeff[i] * m_diis_solutions[i];
        }

        return diisSolution;
    } else {

        Vector coeff = B.colPivHouseholderQr().solve(rhs);

        Vector diisSolution = Vector::Zero(m_diis_solutions[0].size());
        for (int i = 0; i < numErrors; ++i) {
            diisSolution += coeff[i] * m_diis_solutions[i];
        }

        return diisSolution;
    }
}
Vector LBFGS::DIISStep()
{
    double energy = m_energy;

    if (m_diis_solutions.size() == m_diis_hist) {
        m_diis_solutions.erase(m_diis_solutions.begin());
        m_diis_errors.erase(m_diis_errors.begin());
    }
    getEnergyGradient(m_current_coordinates);

    m_diis_solutions.push_back(m_current_coordinates);
    m_diis_errors.push_back(m_gradient);

    auto tmp = diisExtrapolation();
    m_current_coordinates = tmp;
    return m_current_coordinates;
}

Vector LBFGS::LBFGSStep()
{
    double energy = m_energy;
    if (stepCount == 0) {
        getEnergyGradient(m_current_coordinates);
        m_search_direction = -m_gradient;
    } else {
        Vector q = m_gradient;
        std::vector<double> alpha_list;

        for (int i = static_cast<int>(m_step_history.size()) - 1; i >= 0; --i) {
            double alpha = m_rho_history[i] * m_step_history[i].dot(q);
            alpha_list.push_back(alpha);
            q -= alpha * m_gradient_history[i];
        }

        Vector z = q;

        for (size_t i = 0; i < m_step_history.size(); ++i) {
            double beta = m_rho_history[i] * m_gradient_history[i].dot(z);
            z += m_step_history[i] * (alpha_list[m_step_history.size() - i - 1] - beta);
        }

        m_search_direction = -z;
    }
    // m_step_size = 2;
    // m_step_size = line_search_backtracking(x,p, m_step_size, 1e-1, 1e-1);
    // m_step_size = lineSearchRFO(x, p);

    Vector x_old = m_current_coordinates;
    Vector gradient_old = m_gradient;
    m_current_coordinates += m_step_size * m_search_direction;
    getEnergyGradient(m_current_coordinates);
    if ((energy - m_energy) < 0) {
        m_step_size *= 0.9;
        //    std::cout << m_step_size << std::endl;
    } else {
        //  m_step_size *= 1.01;
    }
    Vector s = m_current_coordinates - x_old;
    Vector y = m_gradient - gradient_old;

    double rho = 1.0 / y.dot(s);

    m_diis_solutions.push_back(m_current_coordinates);
    m_diis_errors.push_back(m_gradient);

    if (m_step_history.size() == static_cast<size_t>(m_memory_size)) {
        m_step_history.erase(m_step_history.begin());
        m_gradient_history.erase(m_gradient_history.begin());
        m_rho_history.erase(m_rho_history.begin());
    }

    m_step_history.push_back(s);
    m_gradient_history.push_back(y);
    m_rho_history.push_back(rho);

    return m_current_coordinates;
}

void LBFGS::updateHessian()
{
    Eigen::SelfAdjointEigenSolver<Matrix> solver(m_hessian);
    if (solver.info() != Eigen::Success) {
        std::cerr << "Eigenvalue decomposition failed!" << std::endl;
    }

    // Erhalte die Eigenwerte und Eigenvektoren
    m_eigenvalues = solver.eigenvalues();
    m_eigenvectors = solver.eigenvectors();

    // Bereite die RFO-Koeffizienten vor
    Vector tau(m_atoms);
    for (int i = 0; i < m_atoms; ++i) {
        tau(i) = 1.0 / (m_eigenvalues(i) + m_lambda); // Rationale Funktionsannäherung
    }
    // Aktualisiere das Hessian mit RFO
    Matrix B = Eigen::MatrixXd::Identity(3 * m_atoms, 3 * m_atoms);
    for (int i = 0; i < 3 * m_atoms; i++) {
        B += tau(i) * m_eigenvectors.col(i) * m_eigenvectors.col(i).transpose();
    }
    m_hessian = (m_hessian + B) / 2;
}

Vector LBFGS::RFOStep()
{
    // Wähle den minimalsten/niedrigsten (negativsten) Eigenwert & Eigenvektor
    int minIndex = 0;
    m_eigenvalues.minCoeff(&minIndex);
    //  std::cout << minIndex << "  " << m_eigenvalues.minCoeff(&minIndex) << std::endl;
    Vector direction = m_eigenvectors.col(minIndex);
    Vector gradient_mass_weighted = m_gradient;
    for (int i = 0; i < m_gradient.size(); ++i) {
        gradient_mass_weighted[i] /= std::sqrt(m_masses[i]);
        //  std::cout << std::sqrt(m_masses[i]) << " "  << std::endl;// 3 Dimensionen (x, y, z) pro Atom
    }

    Vector gradient_normal = m_eigenvectors.transpose() * gradient_mass_weighted;
    // std::cout << gradient_normal.transpose() << std::endl << m_gradient.transpose() << std::endl;
    //  Berechne den Schritt: nutze direction für EVF
    double lambda = m_eigenvalues[minIndex]; // kleinster Eigenwert
    double tau = 0.5; // Dämpfungsfaktor zur vorgeschlagenen Schrittgröße

    Vector p_rfo = -tau * (gradient_normal - lambda * direction); // RFO Schritte mit EVF

    // Maßnahme: Line-Search (falls notwendig)
    // if(std::abs(m_step_size  + 1)<1e-10)
    // m_step_size = line_search_backtracking(x,p_rfo, 1, 1e-4);
    //   m_step_size = 1e-4;
    //   std::cout << m_step_size << " ";
    m_step_size = -0.1; // lineSearchRFO(x, p_rfo);
    // m_step_size = 1;
    //    std::cout << m_step_size << " " <<std::endl;

    // Aktualisiere die aktuelle Position
    Vector x_new = m_current_coordinates + m_step_size * p_rfo;
    Vector gradient_prev = m_gradient;
    getEnergyGradient(x_new);
    // updateHessian();

    // Update des Hessians (SR1)
    Eigen::VectorXd s = x_new - m_current_coordinates;
    Eigen::VectorXd y = m_gradient - gradient_prev;
    Matrix H = sr1Update(s, y);
    setHessian(H);
    // std::cout << ( x - x_new).transpose() << std::endl;
    m_current_coordinates = x_new;
    return x_new;
}
Matrix LBFGS::sr1Update(const Eigen::VectorXd& s, const Eigen::VectorXd& y)
{
    Eigen::VectorXd diff = y - m_hessian * s;
    double denom = diff.dot(s);
    if (std::abs(denom) > 1e-10) {
        return m_hessian + m_lambda * (diff * diff.transpose()) / denom;
    }
    return m_hessian;
}
Vector LBFGS::getCurrentSolution() const
{
    return m_current_coordinates;
}

double LBFGS::getCurrentObjectiveValue() const
{
    return m_energy;
}

Vector LBFGS::getCurrentGradient() const
{
    return m_gradient;
}

double LBFGS::getEnergy(const Vector& x)
{
    double fx = 0.0;
    double charge = 0;
    m_interface->updateGeometry(x);
    if (m_interface->HasNan()) {
        m_error = true;
        return 0;
    }
    fx = m_interface->CalculateEnergy(false);
    m_energy = fx;
    return fx;
}

double LBFGS::getEnergyGradient(const Vector& x)
{
    double fx = 0.0;
    double charge = 0;
    m_interface->updateGeometry(x);
    if (m_interface->HasNan()) {
        m_error = true;
        return 0;
    }
    fx = m_interface->CalculateEnergy(true);
    auto gradient = m_interface->Gradient();
    m_error = std::isnan(fx);

    for (int i = 0; i < m_atoms; ++i) {
        m_gradient[3 * i + 0] = gradient(i, 0) * (m_constraints[i]);
        m_gradient[3 * i + 1] = gradient(i, 1) * (m_constraints[i]);
        m_gradient[3 * i + 2] = gradient(i, 2) * (m_constraints[i]);
    }
    m_energy = fx;
    return fx;
}
