/*
 * <Optimizer Driver Base Class - Template Method Pattern>
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

#include "optimizer_interface.h"
#include "src/capabilities/rmsd.h"
#include <chrono>
#include <memory>

namespace Optimization {

/**
 * @brief Optimization context for shared state management - Claude Generated
 * Analog to QM system's parameter management, contains all shared data
 */
class OptimizationContext {
public:
    // Energy calculation
    EnergyCalculator* energy_calculator = nullptr;

    // Convergence criteria (using unit system from CLAUDE.md)
    double energy_threshold = 0.1; // kJ/mol - relative energy changes
    double rmsd_threshold = 0.01; // Å - geometric changes
    double gradient_threshold = 5e-4; // Eh/Bohr - gradient norm
    int max_iterations = 5000;
    int convergence_count = 11; // Bit field for convergence requirements

    // Performance settings
    int threads = 1;
    bool single_step_mode = false;
    double max_energy_rise = 100.0; // kJ/mol - maximum allowed energy increase

    // Output control
    bool write_trajectory = true;
    int verbosity = 1; // Default to minimal output (0=silent, 1=minimal, 2=scientific, 3=debug)
    bool print_output = true;
    std::string trajectory_filename;
    std::string output_basename;

    // Constraints and advanced features
    std::vector<int> atom_constraints; // 0 = fixed, 1 = mobile
    bool use_constraints = false;
    bool use_hessian = false;
    Matrix initial_hessian;

    // Molecular properties
    int charge = 0;
    int spin = 0;

    // RMSD calculation for convergence
    std::unique_ptr<RMSDDriver> rmsd_driver = nullptr;

    // Validation methods
    bool isValid() const;
    std::string getValidationErrors() const;

    // Factory method from JSON
    static OptimizationContext fromJson(const json& config, EnergyCalculator* calc);
};

/**
 * @brief Base class for optimizers using Template Method Pattern - Claude Generated
 * Analog to QMDriver - provides common functionality with customizable steps
 */
class OptimizerDriver : public OptimizerInterface {
protected:
    // Common data members (analog to QMDriver)
    Molecule m_molecule; // Current molecule
    std::vector<Molecule> m_trajectory; // Optimization trajectory
    std::vector<double> m_energy_trajectory; // Energy at each step
    Vector m_current_gradient; // Current gradient
    double m_current_energy = 0.0; // Current energy
    double m_initial_energy = 0.0; // Starting energy

    // Configuration and context
    OptimizationContext m_context;
    json m_configuration;

    // Convergence tracking
    int m_current_iteration = 0;
    bool m_converged = false;
    std::string m_convergence_reason;

    // Timing
    std::chrono::high_resolution_clock::time_point m_start_time;
    std::chrono::high_resolution_clock::time_point m_end_time;

    // Pure virtual methods for derived classes (Template Method Pattern)
    virtual bool InitializeOptimizerInternal() = 0;
    virtual Vector CalculateOptimizationStep(const Vector& current_coordinates,
        const Vector& gradient)
        = 0;
    virtual bool CheckMethodSpecificConvergence() const = 0;
    virtual void UpdateOptimizerState(const Vector& new_coordinates,
        const Vector& new_gradient,
        double new_energy)
        = 0;
    virtual void FinalizeOptimizationInternal() = 0;

    // Common helper methods
    bool evaluateEnergyAndGradient(const Vector& coordinates, double& energy, Vector& gradient);
    bool checkConvergence(double energy_change, double rmsd_change, double gradient_norm) const;
    void updateTrajectory(const Molecule& new_structure, double energy);
    void logOptimizationStep(int iteration, double energy, double energy_change,
        double rmsd_change, double gradient_norm) const;
    double calculateRMSD(const Molecule& mol1, const Molecule& mol2) const;

public:
    OptimizerDriver();
    virtual ~OptimizerDriver() = default;

    // OptimizerInterface implementation (Template Method)
    bool InitializeOptimization(const Molecule& molecule) override;
    bool InitializeOptimization(const Molecule* molecule) override;
    bool InitializeOptimization(const double* coordinates, int atom_count) override;

    bool UpdateGeometry(const Molecule& molecule) override;
    bool UpdateGeometry(const double* coordinates) override;

    // Main optimization method (Template Method implementation)
    OptimizationResult Optimize(bool write_trajectory = false, int verbosity = 1) final override;

    // Property accessors
    Vector GetCurrentGradient() const override { return m_current_gradient; }
    double GetCurrentEnergy() const override { return m_current_energy; }
    Matrix GetCurrentHessian() const override { return Matrix::Zero(0, 0); } // Default: no Hessian
    std::vector<Molecule> GetTrajectory() const override { return m_trajectory; }

    // Configuration management
    void LoadConfiguration(const json& config) override;
    json GetDefaultConfiguration() const override { return OptimizerInterfaceJson; }
    json GetCurrentConfiguration() const override { return m_configuration; }

    // Setters for context
    void setEnergyCalculator(EnergyCalculator* calculator);
    void setConstraints(const std::vector<int>& constraints);
    void setConvergenceCriteria(double energy_thresh, double rmsd_thresh, double grad_thresh);
    void setTrajectoryFile(const std::string& filename);
    void setBasename(const std::string& basename);
};

} // namespace Optimization