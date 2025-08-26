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
    // 🧪 DIIS Extrapolation - Pulay (1980) Chem. Phys. Lett. 73, 393 - Claude Corrected
    int numErrors = m_diis_errors.size();
    Matrix B = Matrix::Zero(numErrors + 1, numErrors + 1);
    Vector rhs = Vector::Zero(numErrors + 1);
    rhs(numErrors) = 1.0;  // Constraint: Σ c_i = 1

    // Build B-matrix: B_ij = <e_i|e_j> (simple dot products, NOT Hessian-weighted)
    for (int i = 0; i < numErrors; ++i) {
        for (int j = 0; j < numErrors; ++j) {
            // CORRECTED: Use direct error vector dot products (Pulay's original method)
            B(i, j) = m_diis_errors[i].dot(m_diis_errors[j]);
        }
        // Lagrange multiplier constraints
        B(numErrors, i) = 1.0;
        B(i, numErrors) = 1.0;
    }

    // Check condition number for numerical stability
    Matrix B_errors = B.block(0, 0, numErrors, numErrors);  // Extract error-error block
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(B_errors);
    double conditionNumber = svd.singularValues()(0) / svd.singularValues().tail(1)(0);

    // If B-matrix is ill-conditioned, remove oldest error vector and rebuild
    if (conditionNumber > 1e12) {
        m_diis_solutions.erase(m_diis_solutions.begin());
        m_diis_errors.erase(m_diis_errors.begin());
        
        // Recursively call with reduced history
        if (m_diis_errors.size() >= 2) {
            return diisExtrapolation();
        } else {
            // Fall back to simple gradient step if insufficient history
            return m_current_coordinates - 0.1 * m_gradient;
        }
    }

    // Solve DIIS equations: B * c = rhs where c = [c_1, ..., c_n, λ]
    Vector coeff = B.colPivHouseholderQr().solve(rhs);

    // Construct DIIS solution: x_new = Σ c_i * x_i
    Vector diisSolution = Vector::Zero(m_diis_solutions[0].size());
    for (int i = 0; i < numErrors; ++i) {
        diisSolution += coeff[i] * m_diis_solutions[i];
    }

    return diisSolution;
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
        // 🧪 L-BFGS Two-Loop Recursion - Nocedal & Wright Algorithm 7.4 - Claude Corrected
        Vector q = m_gradient;  // Start with current gradient
        std::vector<double> alpha_list;
        alpha_list.reserve(m_step_history.size());

        // First loop: compute α_i and update q (backward through history)
        for (int i = static_cast<int>(m_step_history.size()) - 1; i >= 0; --i) {
            double alpha = m_rho_history[i] * m_step_history[i].dot(q);
            alpha_list.push_back(alpha);
            q -= alpha * m_gradient_history[i];
        }

        // Apply initial Hessian approximation H_0 = γI (simplified to identity)
        Vector z = q;  // z = H_0 * q, H_0 = I

        // Second loop: apply corrections (forward through history) 
        for (size_t i = 0; i < m_step_history.size(); ++i) {
            double beta = m_rho_history[i] * m_gradient_history[i].dot(z);
            // CORRECTED: proper indexing from end of alpha_list
            double alpha = alpha_list[m_step_history.size() - 1 - i];
            z += m_step_history[i] * (alpha - beta);
        }

        m_search_direction = -z;  // Search direction: -H * g
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
    Vector s = m_current_coordinates - x_old;  // Step difference
    Vector y = m_gradient - gradient_old;     // Gradient difference

    // 🧪 CURVATURE CONDITION CHECK - Claude Corrected
    // Ensure y·s > 0 for positive definite Hessian approximation (Nocedal & Wright)
    double sy = y.dot(s);
    if (sy > 1e-10) {  // Sufficient curvature condition
        double rho = 1.0 / sy;
        
        // Store LBFGS history with proper bounds checking
        if (m_step_history.size() == static_cast<size_t>(m_memory_size)) {
            m_step_history.erase(m_step_history.begin());
            m_gradient_history.erase(m_gradient_history.begin());
            m_rho_history.erase(m_rho_history.begin());
        }

        m_step_history.push_back(s);
        m_gradient_history.push_back(y);
        m_rho_history.push_back(rho);
    }
    // If curvature condition fails, skip this update (keeps previous Hessian info)

    // DIIS data (separate from LBFGS)
    m_diis_solutions.push_back(m_current_coordinates);
    m_diis_errors.push_back(m_gradient);

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
    // 🧪 Rational Function Optimization - Banerjee et al. (1985) J. Phys. Chem. 89, 52 - Claude Corrected
    
    // Mass-weight the gradient for proper RFO implementation
    Vector gradient_mass_weighted = Vector::Zero(m_gradient.size());
    for (int i = 0; i < m_gradient.size(); ++i) {
        int atom_index = i / 3;  // Map coordinate index to atom index
        gradient_mass_weighted[i] = m_gradient[i] / std::sqrt(m_masses[atom_index]);
    }
    
    // Mass-weight the Hessian (should now be properly initialized)
    Matrix hessian_mass_weighted = m_hessian;
    for (int i = 0; i < m_hessian.rows(); ++i) {
        int atom_i = i / 3;
        for (int j = 0; j < m_hessian.cols(); ++j) {
            int atom_j = j / 3;
            hessian_mass_weighted(i,j) /= std::sqrt(m_masses[atom_i] * m_masses[atom_j]);
        }
    }
    
    // Eigendecomposition of mass-weighted Hessian
    Eigen::SelfAdjointEigenSolver<Matrix> eigensolver(hessian_mass_weighted);
    if (eigensolver.info() != Eigen::Success) {
        // Fallback to simple gradient step if eigendecomposition fails
        m_current_coordinates -= 0.1 * m_gradient;
        getEnergyGradient(m_current_coordinates);
        return m_current_coordinates;
    }
    
    m_eigenvalues = eigensolver.eigenvalues();
    m_eigenvectors = eigensolver.eigenvectors();
    
    // Find lowest eigenvalue for RFO shift
    int minIndex = 0;
    double lambda_min = m_eigenvalues.minCoeff(&minIndex);
    
    // RFO eigenvalue shift using configurable parameter
    double lambda_shift = 0.0;
    if (lambda_min < 0.0) {
        lambda_shift = lambda_min - m_eigenvalue_shift;  // Shift beyond negative eigenvalue
    } else {
        lambda_shift = -m_eigenvalue_shift;  // Small negative shift for positive definite case
    }
    
    // Solve: (H + λI) * Δx = -g  using eigendecomposition
    // (H + λI) = Q * (Λ + λI) * Q^T
    Vector g_transformed = m_eigenvectors.transpose() * gradient_mass_weighted;
    Vector delta_x_transformed = Vector::Zero(g_transformed.size());
    
    for (int i = 0; i < g_transformed.size(); ++i) {
        double shifted_eigenvalue = m_eigenvalues[i] + lambda_shift;
        if (std::abs(shifted_eigenvalue) > 1e-12) {
            delta_x_transformed[i] = -g_transformed[i] / shifted_eigenvalue;
        }
    }
    
    // Transform back to original coordinates
    Vector delta_x_mass_weighted = m_eigenvectors * delta_x_transformed;
    
    // Convert from mass-weighted to Cartesian coordinates
    Vector delta_x = Vector::Zero(delta_x_mass_weighted.size());
    for (int i = 0; i < delta_x.size(); ++i) {
        int atom_index = i / 3;
        delta_x[i] = delta_x_mass_weighted[i] / std::sqrt(m_masses[atom_index]);
    }
    
    // Apply step with configurable trust radius control
    double step_norm = delta_x.norm();
    if (step_norm > m_trust_radius) {
        delta_x *= m_trust_radius / step_norm;
    }
    
    // Store old state for potential backtracking
    Vector x_old = m_current_coordinates;
    Vector gradient_old = m_gradient;
    double energy_old = m_energy;
    
    // Try the full step
    m_current_coordinates += delta_x;
    getEnergyGradient(m_current_coordinates);
    
    // Line search with configurable energy threshold
    if (m_energy > energy_old + m_energy_threshold) {
        // Backtrack with smaller step
        m_current_coordinates = x_old + 0.5 * delta_x;
        getEnergyGradient(m_current_coordinates);
        
        // If still bad, use even smaller step
        if (m_energy > energy_old + m_energy_threshold) {
            m_current_coordinates = x_old + 0.1 * delta_x;
            getEnergyGradient(m_current_coordinates);
            
            // If STILL bad, reject step entirely
            if (m_energy > energy_old + m_energy_threshold) {
                m_current_coordinates = x_old;  // Restore original coordinates
                m_gradient = gradient_old;
                m_energy = energy_old;
                // Reduce trust radius for next iteration
                m_trust_radius = std::max(m_trust_radius * 0.5, m_trust_radius_min);
            }
        }
    } else if (m_energy < energy_old - m_energy_threshold) {
        // Good step, increase trust radius
        m_trust_radius = std::min(m_trust_radius * 1.2, m_trust_radius_max);
    }
    
    // Update Hessian using SR1 update with actual step taken
    Vector s = m_current_coordinates - x_old;
    Vector y = m_gradient - gradient_old;
    if (s.norm() > 1e-12) {  // Only update if we actually moved
        Matrix H_new = sr1Update(s, y);
        setHessian(H_new);
    }
    
    return m_current_coordinates;
}
Matrix LBFGS::sr1Update(const Eigen::VectorXd& s, const Eigen::VectorXd& y)
{
    // 🧪 SR1 Hessian Update - Broyden (1970), Dennis & Moré (1977) - Claude Corrected
    // Formula: H_new = H + ((y - H*s) ⊗ (y - H*s)) / ((y - H*s)·s)
    Eigen::VectorXd diff = y - m_hessian * s;  // y - H*s
    double denom = diff.dot(s);                 // (y - H*s)·s
    
    // Skip update if denominator too small (would make Hessian singular)
    if (std::abs(denom) > 1e-12) {
        // CORRECTED: Remove m_lambda factor (belongs to RFO, not SR1)
        return m_hessian + (diff * diff.transpose()) / denom;
    }
    
    // Keep previous Hessian if update would be numerically unstable
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
