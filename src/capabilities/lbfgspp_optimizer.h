/*
 * <LBFGSpp Optimizer Strategy - External Library Implementation>
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
 *
 */

#pragma once

#include "optimizer_driver.h"
#include <LBFGS.h>
#include <LBFGSB.h>

using namespace LBFGSpp;
using Eigen::VectorXd;

namespace Optimization {

/**
 * @brief LBFGSpp external library interface - Claude Generated
 * Adapter between LBFGSpp and Curcuma's energy calculator
 */
class LBFGSppObjectiveFunction {
public:
    LBFGSppObjectiveFunction(EnergyCalculator* calc, Molecule* mol, const std::vector<int>& constraints);

    // LBFGSpp required interface
    double operator()(const VectorXd& x, VectorXd& grad);

    // State accessors
    double getLastEnergy() const { return m_last_energy; }
    Vector getLastParameters() const { return m_last_parameters; }
    bool hasError() const { return m_error; }
    void clearError() { m_error = false; }

private:
    EnergyCalculator* m_energy_calculator;
    Molecule* m_molecule;
    std::vector<int> m_constraints;

    double m_last_energy = 0.0;
    Vector m_last_parameters;
    bool m_error = false;
    int m_atom_count = 0;
};

/**
 * @brief LBFGSpp Strategy Implementation - Claude Generated
 * Concrete optimizer using external LBFGSpp library
 */
class LBFGSppOptimizer : public OptimizerDriver {
private:
    // LBFGSpp specific parameters
    int m_lbfgs_m = 2000; // Memory parameter
    int m_lbfgs_past = 0; // History size for delta-based stopping
    double m_lbfgs_eps_abs = 1e-5; // Absolute stopping criterion
    double m_lbfgs_eps_rel = 1e-5; // Relative stopping criterion
    double m_lbfgs_delta = 0.0; // Delta for past function values
    int m_lbfgs_line_search = 3; // Line search algorithm
    int m_lbfgs_max_line_search = 20; // Maximum line search trials
    double m_lbfgs_min_step = 1e-20; // Minimum step length
    double m_lbfgs_max_step = 1e20; // Maximum step length
    double m_lbfgs_ftol = 1e-4; // Function tolerance for line search
    double m_lbfgs_wolfe = 0.9; // Wolfe condition parameter

    // LBFGSpp objects
    std::unique_ptr<LBFGSSolver<double>> m_solver;
    std::unique_ptr<LBFGSParam<double>> m_param;
    std::unique_ptr<LBFGSppObjectiveFunction> m_objective;

    // State tracking
    Vector m_current_coordinates;
    bool m_solver_converged = false;

protected:
    // OptimizerDriver interface implementation
    bool InitializeOptimizerInternal() override;
    Vector CalculateOptimizationStep(const Vector& current_coordinates,
        const Vector& gradient) override;
    bool CheckMethodSpecificConvergence() const override;
    void UpdateOptimizerState(const Vector& new_coordinates,
        const Vector& new_gradient,
        double new_energy) override;
    void FinalizeOptimizationInternal() override;

    // Helper methods
    void loadLBFGSParameters(const json& config);
    void configureLBFGSParam();

public:
    LBFGSppOptimizer();
    virtual ~LBFGSppOptimizer() = default;

    // OptimizerInterface implementation
    std::string getName() const override { return "LBFGSpp External Library"; }
    OptimizerType getType() const override { return OptimizerType::LBFGSPP; }
    bool supportsConstraints() const override { return true; }
    bool requiresHessian() const override { return false; }
    std::vector<std::string> getRequiredParameters() const override;

    json GetDefaultConfiguration() const override;

    // Configuration accessors
    void setMemoryParameter(int m) { m_lbfgs_m = m; }
    void setConvergenceParameters(double eps_abs, double eps_rel)
    {
        m_lbfgs_eps_abs = eps_abs;
        m_lbfgs_eps_rel = eps_rel;
    }
    void setLineSearchParameters(int algorithm, int max_trials, double ftol, double wolfe)
    {
        m_lbfgs_line_search = algorithm;
        m_lbfgs_max_line_search = max_trials;
        m_lbfgs_ftol = ftol;
        m_lbfgs_wolfe = wolfe;
    }
};

} // namespace Optimization