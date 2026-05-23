/*
 * <LBFGSpp Optimizer Strategy Implementation>
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

#include "lbfgspp_optimizer.h"
#include "src/core/citation_registry.h"
#include "src/tools/general.h"
#include <stdexcept>

namespace Optimization {

// Claude Generated - LBFGSpp objective function adapter
LBFGSppObjectiveFunction::LBFGSppObjectiveFunction(EnergyCalculator* calc, Molecule* mol,
    const std::vector<int>& constraints, bool use_numerical_gradient, double numerical_gradient_step)
    : m_energy_calculator(calc)
    , m_molecule(mol)
    , m_constraints(constraints)
    , m_use_numerical_gradient(use_numerical_gradient)
    , m_numerical_gradient_step(numerical_gradient_step)
{
    if (mol) {
        m_atom_count = mol->AtomCount();
    }

    // Initialize constraints if empty (all atoms mobile)
    if (m_constraints.empty() && mol) {
        m_constraints.assign(mol->AtomCount(), 1);
    }
}

double LBFGSppObjectiveFunction::operator()(const Vector& x, Vector& grad)
{
    m_error = false;

    try {
        // Update molecule geometry
        CoordinatesToMolecule(x, *m_molecule);

        // Check for NaN coordinates
        for (int i = 0; i < x.size(); ++i) {
            if (std::isnan(x[i]) || std::isinf(x[i])) {
                m_error = true;
                return 0.0;
            }
        }

        // Update geometry only (avoid full GFN-FF reinit via setMolecule)
        m_energy_calculator->updateGeometry(m_molecule->getGeometry());

        double energy;
        if (m_use_numerical_gradient) {
            // Use numerical gradient (for debugging)
            energy = m_energy_calculator->CalculateEnergy(false);
            if (std::isnan(energy) || std::isinf(energy)) {
                m_error = true;
                return 0.0;
            }
            Geometry gradient_geom = m_energy_calculator->NumGrad();
            grad = Vector::Map(gradient_geom.data(), gradient_geom.size());
        } else {
            // Use analytical gradient
            energy = m_energy_calculator->CalculateEnergy(true);
            if (std::isnan(energy) || std::isinf(energy)) {
                m_error = true;
                return 0.0;
            }
            Geometry gradient_geom = m_energy_calculator->Gradient();
            grad = Vector::Map(gradient_geom.data(), gradient_geom.size());
        }

        // Apply constraints to gradient
        for (int i = 0; i < m_atom_count && i < m_constraints.size(); ++i) {
            if (m_constraints[i] == 0) { // Fixed atom
                grad[3 * i] = 0.0;
                grad[3 * i + 1] = 0.0;
                grad[3 * i + 2] = 0.0;
            }
        }

        m_last_energy = energy;
        m_last_parameters = x;

        return energy;

    } catch (const std::exception& e) {
        CurcumaLogger::error_fmt("LBFGSpp objective function error: {}", e.what());
        m_error = true;
        return 0.0;
    }
}

// Claude Generated - LBFGSppOptimizer implementation
LBFGSppOptimizer::LBFGSppOptimizer()
    : OptimizerDriver()
{
    // Load default configuration
    m_configuration = GetDefaultConfiguration();
}

bool LBFGSppOptimizer::InitializeOptimizerInternal()
{
    try {
        // Load LBFGS-specific parameters from configuration
        loadLBFGSParameters(m_configuration);

        // Create LBFGSpp parameter object
        m_param = std::make_unique<LBFGSParam<double>>();
        configureLBFGSParam();

        // Create solver with LineSearchBacktracking (matching legacy CurcumaOpt)
        m_solver = std::make_unique<LBFGSSolver<double, LineSearchBacktracking>>(*m_param);

        // Create objective function
        // Claude Generated (Feb 21, 2026): Pass numerical gradient flags from context
        m_objective = std::make_unique<LBFGSppObjectiveFunction>(
            m_context.energy_calculator, &m_molecule, m_context.atom_constraints,
            m_context.use_numerical_gradient, m_context.numerical_gradient_step);

        // Initialize coordinate vector and LBFGSpp single-step workspace
        m_current_coordinates = MoleculeToCoordinates(m_molecule);
        double fx = m_current_energy;
        int already_converged = m_solver->InitializeSingleSteps(*m_objective, m_current_coordinates, fx);
        if (already_converged) {
            m_solver_converged = true;
        }

        CurcumaLogger::success("LBFGSpp optimizer initialized");

        CitationRegistry::cite("lbfgspp");
        CitationRegistry::cite("lbfgs", "lbfgspp");
        CurcumaLogger::param("Memory parameter (m)", m_lbfgs_m);
        CurcumaLogger::param("Absolute tolerance", fmt::format("{:.2e}", m_lbfgs_eps_abs));
        CurcumaLogger::param("Relative tolerance", fmt::format("{:.2e}", m_lbfgs_eps_rel));
        CurcumaLogger::param("Line search algorithm", m_lbfgs_line_search);

        // Size advisory for large systems (informational, verbosity >= 2)
        int dof = m_current_coordinates.size();
        if (dof > 3000) {
            double mb_history = static_cast<double>(dof) * m_lbfgs_m * 8.0 / (1024.0 * 1024.0);
            CurcumaLogger::info_fmt(
                "LBFGSpp: {} DOF ({} atoms) - m={} history vectors, ~{:.0f} MB memory per step. "
                "L-BFGS scales O(N*m) per step and handles large systems well.",
                dof, m_molecule.AtomCount(), m_lbfgs_m, mb_history);
        }

        return true;

    } catch (const std::exception& e) {
        CurcumaLogger::error_fmt("LBFGSpp initialization failed: {}", e.what());
        return false;
    }
}

Vector LBFGSppOptimizer::CalculateOptimizationStep(const Vector& current_coordinates,
    const Vector& gradient)
{
    // Update current coordinates
    m_current_coordinates = current_coordinates;
    Vector x = current_coordinates;
    double energy = m_current_energy;

    try {
        // Perform single LBFGS step — calls objective function internally (line search)
        m_solver->SingleStep(*m_objective, x, energy);

        if (m_objective->hasError()) {
            CurcumaLogger::warn("LBFGSpp objective function reported error");
            m_objective->clearError();
            return Vector::Zero(current_coordinates.size());
        }

        // Check if LBFGSpp considers itself converged
        m_solver_converged = m_solver->isConverged();

    } catch (const std::logic_error& e) {
        // LBFGSpp throws logic_error when line search step < min_step
        if (m_solver->Step() < 1e-8) {
            m_solver_converged = true;
        } else {
            CurcumaLogger::warn_fmt("LBFGSpp logic error: {}", e.what());
        }
        return Vector::Zero(current_coordinates.size());
    } catch (const std::runtime_error& e) {
        // LBFGSpp throws runtime_error when the line search step shrinks below min_step.
        // Treat as convergence: we cannot make further progress, output the current structure.
        const double gnorm = gradient.norm();
        if (gnorm < m_lbfgs_eps_abs * 100.0) {
            CurcumaLogger::info_fmt("LBFGSpp: line search step < min_step with ||g||={:.2e} — converged",
                                    gnorm);
        } else {
            CurcumaLogger::warn_fmt("LBFGSpp: line search failed (||g||={:.2e}), stopping: {}", gnorm, e.what());
        }
        m_solver_converged = true;
        return Vector::Zero(current_coordinates.size());
    }

    // SingleStep already evaluated energy+gradient at new coordinates —
    // signal to OptimizerDriver to skip redundant re-evaluation
    m_context.step_evaluated_energy = true;
    m_context.step_energy = energy;
    m_context.step_gradient = m_solver->final_grad();

    // Apply constraints to the pre-computed gradient
    if (m_context.use_constraints && !m_context.atom_constraints.empty()) {
        for (size_t i = 0; i < m_context.atom_constraints.size() && i < static_cast<size_t>(m_molecule.AtomCount()); ++i) {
            if (m_context.atom_constraints[i] == 0) {
                m_context.step_gradient[3 * i] = 0.0;
                m_context.step_gradient[3 * i + 1] = 0.0;
                m_context.step_gradient[3 * i + 2] = 0.0;
            }
        }
    }

    Vector step = x - current_coordinates;

    if (m_context.verbosity >= 2) {
        fmt::print("  LBFGSpp: alpha={:.3e} ||step||={:.3e} ||grad||={:.3e}{}\n",
                   m_solver->Step(), step.norm(), m_context.step_gradient.norm(),
                   m_solver_converged ? " [converged]" : "");
    }
    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info_fmt("LBFGSpp step norm: {:.6e}", step.norm());
    }

    return step;
}

bool LBFGSppOptimizer::CheckMethodSpecificConvergence() const
{
    // LBFGSpp handles its own convergence criteria
    return m_solver_converged;
}

void LBFGSppOptimizer::UpdateOptimizerState(const Vector& new_coordinates,
    const Vector& new_gradient,
    double new_energy)
{
    m_current_coordinates = new_coordinates;
    // LBFGSpp manages its own internal state
}

void LBFGSppOptimizer::FinalizeOptimizationInternal()
{
    // Cleanup LBFGSpp objects
    m_objective.reset();
    m_solver.reset();
    m_param.reset();

    CurcumaLogger::success("LBFGSpp optimization finalized");
}

json LBFGSppOptimizer::GetDefaultConfiguration() const
{
    json config = OptimizerInterfaceJson;

    // LBFGSpp-specific parameters - matching legacy CurcumaOpt defaults
    config["lbfgs_m"] = 2000;
    config["lbfgs_past"] = 0;
    config["lbfgs_eps_abs"] = 1e-5;
    config["lbfgs_eps_rel"] = 1e-5;
    config["lbfgs_delta"] = 0.0;
    config["lbfgs_line_search"] = 3; // LBFGS_LINESEARCH_DEFAULT
    config["lbfgs_max_line_search"] = 2;  // Legacy: LBFGS_ls_iter = 2
    config["lbfgs_min_step"] = 1e-4;      // Legacy: LBFGS_min_step = 1e-4 (not 1e-20!)
    config["lbfgs_max_step"] = 1e20;
    config["lbfgs_ftol"] = 1e-4;
    config["lbfgs_wolfe"] = 0.9;

    return config;
}

std::vector<std::string> LBFGSppOptimizer::getRequiredParameters() const
{
    return {
        "lbfgs_m", "lbfgs_eps_abs", "lbfgs_eps_rel",
        "lbfgs_line_search", "lbfgs_ftol", "lbfgs_wolfe"
    };
}

// Claude Generated - Helper methods
void LBFGSppOptimizer::loadLBFGSParameters(const json& config)
{
    if (config.contains("lbfgs_m"))
        m_lbfgs_m = config["lbfgs_m"];
    if (config.contains("lbfgs_past"))
        m_lbfgs_past = config["lbfgs_past"];
    if (config.contains("lbfgs_eps_abs"))
        m_lbfgs_eps_abs = config["lbfgs_eps_abs"];
    if (config.contains("lbfgs_eps_rel"))
        m_lbfgs_eps_rel = config["lbfgs_eps_rel"];
    if (config.contains("lbfgs_delta"))
        m_lbfgs_delta = config["lbfgs_delta"];
    if (config.contains("lbfgs_line_search"))
        m_lbfgs_line_search = config["lbfgs_line_search"];
    if (config.contains("lbfgs_max_line_search"))
        m_lbfgs_max_line_search = config["lbfgs_max_line_search"];
    if (config.contains("lbfgs_min_step"))
        m_lbfgs_min_step = config["lbfgs_min_step"];
    if (config.contains("lbfgs_max_step"))
        m_lbfgs_max_step = config["lbfgs_max_step"];
    if (config.contains("lbfgs_ftol"))
        m_lbfgs_ftol = config["lbfgs_ftol"];
    if (config.contains("lbfgs_wolfe"))
        m_lbfgs_wolfe = config["lbfgs_wolfe"];
}

void LBFGSppOptimizer::configureLBFGSParam()
{
    // Set LBFGSpp parameters
    m_param->m = m_lbfgs_m;
    m_param->past = m_lbfgs_past;
    m_param->epsilon = m_lbfgs_eps_abs;
    m_param->epsilon_rel = m_lbfgs_eps_rel;
    m_param->delta = m_lbfgs_delta;

    // Line search parameters
    m_param->linesearch = static_cast<LINE_SEARCH_TERMINATION_CONDITION>(m_lbfgs_line_search);
    m_param->max_linesearch = m_lbfgs_max_line_search;
    m_param->min_step = m_lbfgs_min_step;
    m_param->max_step = m_lbfgs_max_step;
    m_param->ftol = m_lbfgs_ftol;
    m_param->wolfe = m_lbfgs_wolfe;

    // Additional parameters
    m_param->max_iterations = m_context.max_iterations;
}

} // namespace Optimization