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
#include <functional>
#include <memory>

namespace Optimization {

// Coordinate conversion helpers used by multiple optimizer implementations
Vector MoleculeToCoordinates(const Molecule& mol);
void CoordinatesToMolecule(const Vector& coords, Molecule& mol);


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
    int convergence_count = 7; // Bit field: bits 0(energy)+1(RMSD)+2(gradient) all required

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

    // Step-includes-evaluation: when subclass step already evaluated energy/gradient
    // (e.g. LBFGSpp SingleStep calls objective internally), driver skips re-evaluation
    bool step_evaluated_energy = false;
    double step_energy = 0.0;
    Vector step_gradient;

    // Constraints and advanced features
    std::vector<int> atom_constraints; // 0 = fixed, 1 = mobile
    bool use_constraints = false;
    bool use_hessian = false;
    bool use_numerical_gradient = false; // Use numerical gradient instead of analytical (for debugging)
    double numerical_gradient_step = 1e-5; // Step size for numerical gradient (Bohr)
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
 * @brief Abstract base class for all optimizers — Template Method Pattern
 * Provides the optimization loop, convergence checking, and output table.
 * Subclasses implement only their step logic via CalculateOptimizationStep().
 * Claude Generated
 */
class OptimizerDriver {
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

    // Per-step notification callback (set by interactive simulation callers)
    // Signature: (iteration, current_molecule, current_energy) → continue?
    // Return false to request early stop.
    using StepCallback = std::function<bool(int, const Molecule&, double)>;
    StepCallback m_step_callback;

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
        double rmsd_change, double gradient_norm, double elapsed_time) const;
    double calculateRMSD(const Molecule& mol1, const Molecule& mol2) const;

public:
    OptimizerDriver();
    virtual ~OptimizerDriver() = default;

    /** @brief Register a callback invoked after every optimization step.
     *  The callback receives (iteration, current molecule, current energy).
     *  Return false to request early termination (e.g. stop button pressed). */
    void setStepCallback(StepCallback cb) { m_step_callback = std::move(cb); }

    // Initialization
    bool InitializeOptimization(const Molecule& molecule);
    bool InitializeOptimization(const Molecule* molecule);
    bool InitializeOptimization(const double* coordinates, int atom_count);

    // Geometry update
    bool UpdateGeometry(const Molecule& molecule);
    bool UpdateGeometry(const double* coordinates);

    // Main optimization loop — all subclasses use this, no overrides allowed
    OptimizationResult Optimize(bool write_trajectory = false, int verbosity = 1);

    // Property accessors
    Vector GetCurrentGradient() const { return m_current_gradient; }
    double GetCurrentEnergy() const { return m_current_energy; }
    Matrix GetCurrentHessian() const { return Matrix::Zero(0, 0); }
    std::vector<Molecule> GetTrajectory() const { return m_trajectory; }

    // Configuration management
    void LoadConfiguration(const json& config);
    json GetDefaultConfiguration() const { return OptimizerInterfaceJson; }
    json GetCurrentConfiguration() const { return m_configuration; }

    // Setters for context
    void setEnergyCalculator(EnergyCalculator* calculator);
    void setConstraints(const std::vector<int>& constraints);
    void setConvergenceCriteria(double energy_thresh, double rmsd_thresh, double grad_thresh);
    void setTrajectoryFile(const std::string& filename);
    void setBasename(const std::string& basename);

    // Method identification — implemented by each concrete optimizer
    virtual std::string getName() const = 0;
    virtual OptimizerType getType() const = 0;
    virtual bool supportsConstraints() const = 0;
    virtual bool requiresHessian() const = 0;
    virtual std::vector<std::string> getRequiredParameters() const = 0;
};

} // namespace Optimization