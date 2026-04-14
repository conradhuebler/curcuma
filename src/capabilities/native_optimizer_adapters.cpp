/*
 * <Native Optimizer Adapters Implementation>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (Apr 2026)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 */

#include "native_optimizer_adapters.h"

namespace Optimization {

// --- Base adapter ---

NativeOptimizerAdapter::NativeOptimizerAdapter(LBFGS::Method method)
    : OptimizerDriver()
    , m_native_method(method)
{
    m_configuration = OptimizerInterfaceJson;
}

bool NativeOptimizerAdapter::InitializeOptimizerInternal()
{
    try {
        int memory = m_configuration.value("lbfgs_memory", 10);
        m_lbfgs = std::make_unique<LBFGS>(memory);

        // Configure method
        m_lbfgs->setOptimizationMethod(m_native_method);
        m_lbfgs->setEnergyCalculator(m_context.energy_calculator);
        m_lbfgs->setVerbosity(0); // Driver handles output, not the LBFGS class

        // Translate constraints: OptimizerDriver uses per-atom vector<int> (0=fixed, 1=free)
        // LBFGS class uses per-coordinate Vector (0=fixed, 1=free)
        if (m_context.use_constraints && !m_context.atom_constraints.empty()) {
            int dof = 3 * m_molecule.AtomCount();
            Vector constraints = Vector::Ones(dof);
            for (size_t i = 0; i < m_context.atom_constraints.size()
                 && i < static_cast<size_t>(m_molecule.AtomCount()); ++i) {
                if (m_context.atom_constraints[i] == 0) {
                    constraints[3 * i] = 0.0;
                    constraints[3 * i + 1] = 0.0;
                    constraints[3 * i + 2] = 0.0;
                }
            }
            m_lbfgs->setConstraints(constraints);
        }

        // DIIS-specific parameters
        if (m_native_method == LBFGS::Method::DIIS) {
            int diis_hist = m_configuration.value("diis_history", 5);
            int diis_start = m_configuration.value("diis_start", 5);
            m_lbfgs->setDIISParameters(diis_hist, diis_start);
        }

        // RFO-specific parameters
        if (m_native_method == LBFGS::Method::RFO) {
            double trust_radius = m_configuration.value("trust_radius", 0.1);
            double energy_thresh = m_configuration.value("rfo_energy_threshold", 1e-4);
            double eigenvalue_shift = m_configuration.value("rfo_eigenvalue_shift", 1e-3);
            m_lbfgs->setRFOParameters(trust_radius, energy_thresh, eigenvalue_shift);
        }

        // Initialize with current coordinates
        Vector coords = MoleculeToCoordinates(m_molecule);
        m_lbfgs->initialize(m_molecule.AtomCount(), coords);

        return true;
    } catch (const std::exception& e) {
        CurcumaLogger::error_fmt("Native optimizer initialization failed: {}", e.what());
        return false;
    }
}

Vector NativeOptimizerAdapter::CalculateOptimizationStep(const Vector& current_coordinates,
    const Vector& gradient)
{
    try {
        // LBFGS::step() internally calls getEnergyGradient() and returns new coordinates
        Vector new_coords = m_lbfgs->step();

        if (m_lbfgs->hasError()) {
            CurcumaLogger::warn("Native optimizer reported error during step");
            return Vector::Zero(current_coordinates.size());
        }

        // Signal that energy was already evaluated inside step()
        m_context.step_evaluated_energy = true;
        m_context.step_energy = m_lbfgs->getCurrentEnergy();
        m_context.step_gradient = m_lbfgs->getCurrentGradient();

        // Apply constraints to the pre-computed gradient
        if (m_context.use_constraints && !m_context.atom_constraints.empty()) {
            for (size_t i = 0; i < m_context.atom_constraints.size()
                 && i < static_cast<size_t>(m_molecule.AtomCount()); ++i) {
                if (m_context.atom_constraints[i] == 0) {
                    m_context.step_gradient[3 * i] = 0.0;
                    m_context.step_gradient[3 * i + 1] = 0.0;
                    m_context.step_gradient[3 * i + 2] = 0.0;
                }
            }
        }

        // Return displacement, not absolute coordinates
        return new_coords - current_coordinates;

    } catch (const std::exception& e) {
        CurcumaLogger::error_fmt("Native optimizer step failed: {}", e.what());
        return Vector::Zero(current_coordinates.size());
    }
}

bool NativeOptimizerAdapter::CheckMethodSpecificConvergence() const
{
    // No method-specific convergence beyond what OptimizerDriver checks
    return true;
}

void NativeOptimizerAdapter::UpdateOptimizerState(const Vector& new_coordinates,
    const Vector& new_gradient,
    double new_energy)
{
    // LBFGS class manages its own internal state during step()
}

void NativeOptimizerAdapter::FinalizeOptimizationInternal()
{
    m_lbfgs.reset();
}

// --- Concrete adapters ---

NativeLBFGSAdapter::NativeLBFGSAdapter()
    : NativeOptimizerAdapter(LBFGS::Method::LBFGS)
{
}

NativeDIISAdapter::NativeDIISAdapter()
    : NativeOptimizerAdapter(LBFGS::Method::DIIS)
{
}

NativeRFOAdapter::NativeRFOAdapter()
    : NativeOptimizerAdapter(LBFGS::Method::RFO)
{
}

} // namespace Optimization
