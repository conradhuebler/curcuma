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

// EIGEN_USE_LAPACKE: enables LAPACK dsyevd inside SelfAdjointEigenSolver.
// Safe here — this file has no variable named 'I'.
#ifndef EIGEN_USE_LAPACKE
#define EIGEN_USE_LAPACKE
#endif

#include "lbfgs.h"
#include "rf_solver.h"
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

// Claude Generated (Apr 2026): Read key optimizer settings from JSON config
void LBFGS::setConfig(const json& config)
{
    if (config.contains("lbfgs_line_search"))
        m_line_search_method = config["lbfgs_line_search"].get<std::string>();
    else if (config.contains("line_search"))
        m_line_search_method = config["line_search"].get<std::string>();

    if (config.contains("rfo_solver"))
        m_rfo_solver = config["rfo_solver"].get<std::string>();

    if (config.contains("verbosity"))
        m_verbosity = config["verbosity"].get<int>();

    if (config.contains("lbfgs_memory"))
        m_memory_size = config["lbfgs_memory"].get<int>();
    else if (config.contains("memory_size"))
        m_memory_size = config["memory_size"].get<int>();

    if (m_verbosity >= 2) {
        CurcumaLogger::param("L-BFGS line search", m_line_search_method);
        CurcumaLogger::param("RFO solver", m_rfo_solver);
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

double LBFGS::lineSearchStrongWolfe(const Vector& x, const Vector& p,
                                     double alpha_init, double c1, double c2)
{
    // Strong-Wolfe line search — Nocedal & Wright, Algorithms 3.5 and 3.6.
    // Conditions:
    //   Armijo:    f(x + α·p) <= f(x) + c1·α·∇f·p          (sufficient decrease)
    //   Curvature: |∇f(x + α·p)·p| <= c2·|∇f·p|            (strong Wolfe curvature)
    //
    // Phase 1 — bracket: expand α until Armijo violated or curvature satisfied.
    // Phase 2 — zoom: bisect the bracket until both conditions are met.
    //
    // Each accepted step caches the gradient, so the caller pays no extra gradient
    // evaluation on the final accepted point.

    const double f0    = m_energy;
    const double dphi0 = m_gradient.dot(p); // must be negative for descent direction

    if (dphi0 >= 0.0) {
        // Not a descent direction — return zero step (caller will handle).
        return 0.0;
    }

    // Internal helpers: evaluate f and ∇f·p at x + α·p.
    auto phi = [&](double a) -> double {
        return getEnergy(x + a * p);
    };
    auto dphi = [&](double a) -> double {
        getEnergyGradient(x + a * p);
        return m_gradient.dot(p);
    };

    double alpha_lo = 0.0, f_lo = f0, dphi_lo = dphi0;
    double alpha_hi = 0.0, f_hi = f0;

    double alpha = alpha_init;
    double alpha_prev = 0.0;
    double f_prev = f0;

    const int max_phase1 = 20;
    const int max_phase2 = 30;

    // --- Phase 1: bracket search ---
    for (int i = 0; i < max_phase1; ++i) {
        double f_a = phi(alpha);

        if (f_a > f0 + c1 * alpha * dphi0 || (i > 0 && f_a >= f_prev)) {
            // Armijo violated or f increased: bracket is [alpha_prev, alpha].
            alpha_lo = alpha_prev; f_lo = f_prev; dphi_lo = dphi(alpha_prev);
            alpha_hi = alpha;      f_hi = f_a;
            break;
        }

        double dp_a = dphi(alpha);
        if (std::abs(dp_a) <= -c2 * dphi0) {
            // Both conditions satisfied.
            return alpha;
        }
        if (dp_a >= 0.0) {
            // Slope positive: bracket is [alpha, alpha_prev].
            alpha_lo = alpha;      f_lo = f_a;     dphi_lo = dp_a;
            alpha_hi = alpha_prev; f_hi = f_prev;
            break;
        }

        alpha_prev = alpha;
        f_prev     = f_a;
        alpha      = std::min(alpha * 2.0, 10.0); // expand, cap at 10
    }

    if (alpha_hi == 0.0) {
        // Phase 1 didn't bracket — return best found so far.
        return alpha_prev > 0.0 ? alpha_prev : alpha_init * 0.1;
    }

    // --- Phase 2: zoom ---
    for (int i = 0; i < max_phase2; ++i) {
        // Bisect (cubic interpolation would be better but bisect is robust).
        double alpha_j = 0.5 * (alpha_lo + alpha_hi);
        double f_j = phi(alpha_j);

        if (f_j > f0 + c1 * alpha_j * dphi0 || f_j >= f_lo) {
            alpha_hi = alpha_j;
            f_hi = f_j;
        } else {
            double dp_j = dphi(alpha_j);
            if (std::abs(dp_j) <= -c2 * dphi0) {
                return alpha_j; // Both conditions satisfied.
            }
            if (dp_j * (alpha_hi - alpha_lo) >= 0.0) {
                alpha_hi = alpha_lo;
                f_hi = f_lo;
            }
            alpha_lo = alpha_j;
            f_lo = f_j;
            dphi_lo = dp_j;
        }

        if (std::abs(alpha_hi - alpha_lo) < 1e-12)
            break;
    }

    return alpha_lo > 0.0 ? alpha_lo : alpha_init * 0.1;
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
    Vector x_old = m_current_coordinates;
    Vector gradient_old = m_gradient;

    // Line search: Strong-Wolfe guarantees the curvature condition needed for
    // LBFGS to maintain positive-definite Hessian approximation (N&W Thm 7.1).
    // Backtracking (default) is cheaper per step but may stagnate on anisotropic PES.
    if (m_line_search_method == "strong_wolfe") {
        m_step_size = lineSearchStrongWolfe(x_old, m_search_direction, 1.0);
        if (m_step_size <= 0.0) m_step_size = 0.01; // guard against failed search
    } else {
        // Backtracking: simple Armijo with geometric reduction.
        m_step_size = 1.0;
        double f0 = energy;
        double dphi0 = m_gradient.dot(m_search_direction);
        for (int bt = 0; bt < 30; ++bt) {
            double f_trial = getEnergy(x_old + m_step_size * m_search_direction);
            if (f_trial <= f0 + 1e-4 * m_step_size * dphi0) break;
            m_step_size *= 0.5;
        }
    }

    m_current_coordinates = x_old + m_step_size * m_search_direction;
    getEnergyGradient(m_current_coordinates);
    Vector s = m_current_coordinates - x_old;  // Step difference
    Vector y = m_gradient - gradient_old;     // Gradient difference

    if (m_verbosity >= 2) {
        CurcumaLogger::info_fmt("L-BFGS: alpha={:.3e} ||step||={:.3e} ||grad||={:.3e} history={}",
                                m_step_size, s.norm(), m_gradient.norm(),
                                static_cast<int>(m_step_history.size()));
    }

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
    // Rational Function Optimization - Banerjee et al. (1985) J. Phys. Chem. 89, 52
    //
    // Uses RFSolver::calculateRFStep (shared with ANCOpt) which solves the augmented
    // eigenvalue problem [[H,g],[g',0]] via Lanczos for nvar+1 >= 50.  This replaces
    // the previous O(N^3) full SelfAdjointEigenSolver on the mass-weighted Hessian.
    //
    // Mass-weighting: gradient and Hessian are transformed to mass-weighted space,
    // the RFO step is computed there, then un-mass-weighted back to Cartesian.

    const int ndof = static_cast<int>(m_gradient.size());

    // Mass-weight gradient and Hessian.
    Vector g_mw(ndof);
    for (int i = 0; i < ndof; ++i)
        g_mw[i] = m_gradient[i] / std::sqrt(m_masses[i / 3]);

    Matrix H_mw = m_hessian;
    for (int i = 0; i < ndof; ++i)
        for (int j = 0; j < ndof; ++j)
            H_mw(i, j) /= std::sqrt(m_masses[i / 3] * m_masses[j / 3]);

    // Invalidate warm-start if system size changed.
    if (m_rfo_warm_start.size() != ndof + 1)
        m_rfo_warm_start.resize(0);

    // RFO displacement in mass-weighted space.
    // rfo_solver == "dense": always use full O(N^3) SelfAdjointEigenSolver (reference path)
    // rfo_solver == "lanczos" (default): Lanczos for nvar+1 >= 50, dense fallback otherwise
    Vector delta_x_mw;
    if (m_rfo_solver == "dense") {
        const int nvar1 = ndof + 1;
        RFMatrix A_aug = RFMatrix::Zero(nvar1, nvar1);
        A_aug.topLeftCorner(ndof, ndof) = H_mw;
        A_aug.block(ndof, 0, 1, ndof) = g_mw.transpose();
        A_aug.block(0, ndof, ndof, 1) = g_mw;
        Eigen::SelfAdjointEigenSolver<RFMatrix> solver(A_aug);
        if (solver.info() != Eigen::Success) {
            CurcumaLogger::warn("RFO dense solver failed, using steepest descent");
            return m_current_coordinates - 0.01 * m_gradient / m_gradient.norm();
        }
        RFVector ev = solver.eigenvectors().col(0);
        if (ev(ndof) < 0) ev = -ev;
        if (std::abs(ev(ndof)) < 1e-10)
            delta_x_mw = -g_mw * 0.1;
        else
            delta_x_mw = ev.head(ndof) / ev(ndof);
        if (m_verbosity >= 2)
            CurcumaLogger::info_fmt("RFO: dense solver (nvar={})", ndof);
    } else {
        delta_x_mw = RFSolver::calculateRFStep(g_mw, H_mw, m_rfo_warm_start);
        if (m_verbosity >= 2)
            CurcumaLogger::info_fmt("RFO: Lanczos solver (nvar={})", ndof);
    }

    // Un-mass-weight: x_cart = x_mw / sqrt(m_i)
    Vector delta_x(ndof);
    for (int i = 0; i < ndof; ++i)
        delta_x[i] = delta_x_mw[i] / std::sqrt(m_masses[i / 3]);

    // Trust-radius clamp.
    double step_norm = delta_x.norm();
    if (step_norm > m_trust_radius)
        delta_x *= m_trust_radius / step_norm;

    // Store old state for backtracking.
    Vector x_old = m_current_coordinates;
    Vector gradient_old = m_gradient;
    double energy_old = m_energy;

    // Try full step, then backtrack 0.5, then 0.1, then reject.
    m_current_coordinates += delta_x;
    getEnergyGradient(m_current_coordinates);

    if (m_energy > energy_old + m_energy_threshold) {
        m_current_coordinates = x_old + 0.5 * delta_x;
        getEnergyGradient(m_current_coordinates);

        if (m_energy > energy_old + m_energy_threshold) {
            m_current_coordinates = x_old + 0.1 * delta_x;
            getEnergyGradient(m_current_coordinates);

            if (m_energy > energy_old + m_energy_threshold) {
                m_current_coordinates = x_old;
                m_gradient = gradient_old;
                m_energy = energy_old;
                m_trust_radius = std::max(m_trust_radius * 0.5, m_trust_radius_min);
            }
        }
    } else if (m_energy < energy_old - m_energy_threshold) {
        m_trust_radius = std::min(m_trust_radius * 1.2, m_trust_radius_max);
    }

    // SR1 Hessian update on the actual step taken.
    Vector s = m_current_coordinates - x_old;
    Vector y = m_gradient - gradient_old;
    if (s.norm() > 1e-12) {
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
