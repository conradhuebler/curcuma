/*
 * <Native LBFGS Optimizer - Claude 3.5 Generated, Cleaned by Claude 4>
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

#pragma once

#include "src/core/curcuma_logger.h"
#include "src/core/energycalculator.h"

#include <Eigen/Dense>
#include <vector>

using Vector = Eigen::VectorXd;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;

/**
 * @brief Native LBFGS Optimizer with DIIS and RFO support - Claude Generated
 *
 * This is a comprehensive optimization class that supports:
 * - L-BFGS (Limited-memory Broyden-Fletcher-Goldfarb-Shanno)
 * - DIIS (Direct Inversion of Iterative Subspace)
 * - RFO (Rational Function Optimization)
 * - Various line search algorithms
 */
class LBFGS {
public:
    enum class Method {
        LBFGS = 1, // Standard L-BFGS
        DIIS = 2, // DIIS with LBFGS fallback
        RFO = 3 // Rational Function Optimization
    };

    explicit LBFGS(int memory_size = 10);
    ~LBFGS() = default;

    // Core optimization interface
    void initialize(int atoms, const Vector& initial_coordinates);
    Vector step();

    // Method-specific step functions
    Vector LBFGSStep();
    Vector DIISStep();
    Vector RFOStep();

    // Result accessors
    Vector getCurrentSolution() const;
    double getCurrentObjectiveValue() const;
    double getCurrentEnergy() const { return m_energy; }
    double Energy() const { return m_energy; } // Legacy compatibility - Claude Generated
    Vector getCurrentGradient() const;
    int getStepCount() const { return stepCount; }

    // Status checks
    bool hasError() const { return m_error; }
    bool isError() const { return m_error; } // Legacy compatibility - Claude Generated
    bool isConverged() { return false; } // Placeholder implementation

    // Configuration setters
    void setEnergyCalculator(EnergyCalculator* calculator) { m_interface = calculator; }
    void setConstraints(const Vector& constraints) { m_constraints = constraints; }
    void setOptimizationMethod(Method method) { m_method = method; }
    void setOptimMethod(int method)
    { // Legacy compatibility - Claude Generated
        if (method == 1)
            m_method = Method::LBFGS;
        else if (method == 2)
            m_method = Method::DIIS;
        else if (method == 3)
            m_method = Method::RFO;
        else
            m_method = Method::LBFGS;
    }
    void setLambda(double lambda) { m_lambda = lambda; }
    void setMasses(const std::vector<double>& masses) { m_masses = masses; }
    void setDIISParameters(int history_size, int start_iteration)
    {
        m_diis_hist = history_size;
        m_diis_start = start_iteration;
    }
    void setDIIS(int history_size, int start_iteration)
    { // Legacy compatibility - Claude Generated
        setDIISParameters(history_size, start_iteration);
    }
    void setHessian(const Matrix& hessian);
    void setVerbosity(int level) { m_verbosity = level; }

private:
    // Line search algorithms
    double lineSearchRFO(const Vector& x, const Vector& p, double c1 = 1e-4, double c2 = 0.9);
    double lineSearchBacktracking(const Vector& x, const Vector& p, double alpha_init,
        double rho = 0.5, double c1 = 1e-4);

    // DIIS implementation
    Vector diisExtrapolation();

    // Energy and gradient evaluation
    double getEnergy(const Vector& x);
    double getEnergyGradient(const Vector& x);

    // Hessian management
    void updateHessian();
    Matrix sr1Update(const Eigen::VectorXd& s, const Eigen::VectorXd& y);
    Eigen::MatrixXd update_inverse_hessian_approximation(Eigen::MatrixXd& B,
        const Eigen::VectorXd& deltaX,
        const Eigen::VectorXd& deltaG);
    double line_search_backtracking(const Vector& x, const Vector& p, double alpha_init,
        double rho, double c1 = 1e-4);

    // Core data members
    EnergyCalculator* m_interface = nullptr;

    // State variables
    bool m_error = false;
    int m_atoms = 0;
    int m_memory_size; // L-BFGS memory size (formerly 'm')
    int stepCount = 0;
    int m_verbosity = 2; // CurcumaLogger verbosity level

    // Method selection
    Method m_method = Method::LBFGS;

    // DIIS parameters
    int m_diis_hist = 10; // DIIS history size
    int m_diis_start = 10; // When to start DIIS

    // RFO parameters
    double m_lambda = 0.1; // RFO lambda parameter

    // Step size control
    double m_step_size = 1.0; // Current step size
    double m_last_step_size = 1.0;

    // Current optimization state
    double m_energy = 0.0;
    Vector m_current_coordinates;
    Vector m_gradient;
    Vector m_search_direction;
    Vector m_constraints; // Coordinate constraints (1=free, 0=fixed)

    // Hessian information
    Matrix m_hessian;
    Matrix m_eigenvectors;
    Matrix m_inverse_hessian;
    Vector m_eigenvalues;

    // Mass weighting (for RFO)
    std::vector<double> m_masses;

    // L-BFGS memory
    std::vector<Vector> m_step_history; // s_list
    std::vector<Vector> m_gradient_history; // y_list
    std::vector<double> m_rho_history; // rho_list

    // DIIS memory
    std::vector<Vector> m_diis_solutions; // Previous solutions
    std::vector<Vector> m_diis_errors; // Previous errors (gradients)
};
