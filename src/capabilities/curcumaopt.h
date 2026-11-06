/*
 * <Handling optimisation of structures. >
 * Copyright (C) 2020 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
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

// Only the PARAM block of module "opt" lives here; the registry generator scans this header.
// The optimisation itself is OptimizerDriver (optimizer_driver.h).
#include "src/core/parameter_macros.h"

// ============================================================================
// Parameter Registry Integration - Claude Generated (October 2025)
// ============================================================================
namespace {
    BEGIN_PARAMETER_DEFINITION(opt)
    MODULE_INFO("Geometry optimisation and single-point energies: optimiser choice, convergence and the energy method.", "Optimisation", {"opt", "sp"})

    // Basic
    // master's fuller help text, this branch's tier/bounds annotations.
    PARAM(method, String, "gfnff", "Energy/gradient method. Default gfnff: the fast general-purpose choice (native GFN-FF). gfn2 is the accurate one (native GFN2-xTB, ~100x slower); gfn1, eht, pm3, uff, qmdff also available.", "Basic", {},
        "tier=primary")
    PARAM(optimizer, String, "auto", "Optimizer algorithm: auto, lbfgspp, lbfgs, diis, rfo, ancopt", "Basic", {},
        "tier=primary")
    PARAM(threads, Int, 1, "Number of parallel threads", "Basic", {},
        "min=1")
    PARAM(charge, Int, 0, "Total molecular charge", "Basic", { "Charge" },
        "tier=primary")
    PARAM(spin, Int, 0, "Spin multiplicity", "Basic", { "Spin" },
        "tier=primary")
    PARAM(verbosity, Int, 1, "Output level: 0=silent, 1=table, 2=detailed, 3=debug", "Basic", { "verbose" })

    // Convergence
    PARAM(convergence_preset, String, "normal", "Convergence preset: loose, normal, tight, verytight (sets all thresholds)", "Convergence", { "conv_preset", "ConvPreset" })
    PARAM(energy_threshold, Double, 0.1, "Energy change threshold [kJ/mol]", "Convergence", { "d_e", "dE" })
    PARAM(rmsd_threshold, Double, 0.01, "RMSD change threshold [Angstrom]", "Convergence", { "d_rmsd", "dRMSD" })
    PARAM(gradient_threshold, Double, 5e-4, "Gradient norm threshold [Eh/Bohr]", "Convergence", { "grad_norm", "GradNorm" })
    PARAM(max_iterations, Int, 5000, "Maximum number of optimization steps", "Convergence", { "max_iter", "MaxIter" })
    PARAM(convergence_count, Int, 7, "Convergence bit field: 1=energy, 2=RMSD, 4=gradient, 8=optimizer-internal (7=all main criteria)", "Convergence", { "conv_count", "ConvCount" })
    PARAM(max_energy_rise, Double, 100.0, "Maximum allowed energy rise [kJ/mol] before abort", "Convergence", { "maxrise" })
    PARAM(stall_steps, Int, 20, "Stop as not converged, no progress, after this many consecutive steps in which the geometry moved less than stall_rmsd; 0 disables the check", "Convergence", {})
    PARAM(stall_rmsd, Double, 1e-6, "Per-step RMSD [Angstrom] below which a step counts as no movement for stall_steps", "Convergence", {})

    // Output
    PARAM(write_trajectory, Bool, true, "Write optimization trajectory to .trj.xyz", "Output", { "write_xyz", "writeXYZ" })

    // L-BFGS tuning (lbfgspp optimizer)
    PARAM(lbfgs_m, Int, 2000, "L-BFGS memory (number of stored steps)", "LBFGS", {})
    PARAM(lbfgs_line_search, Int, 3, "Line search type: 1=Armijo, 2=Wolfe, 3=StrongWolfe, 4=Backtracking", "LBFGS", { "lbfgs_lst" })
    PARAM(lbfgs_max_line_search, Int, 20, "Maximum line search iterations", "LBFGS", { "lbfgs_ls_iter" })
    PARAM(lbfgs_min_step, Double, 1e-4, "Minimum step size for line search", "LBFGS", {})
    PARAM(lbfgs_ftol, Double, 1e-4, "Armijo condition parameter (sufficient decrease)", "LBFGS", {})
    PARAM(lbfgs_wolfe, Double, 0.9, "Wolfe condition parameter (curvature)", "LBFGS", {})
    PARAM(lbfgs_eps_abs, Double, 1e-5, "L-BFGS gradient convergence epsilon (absolute)", "LBFGS", {})

    // Native optimizer parameters (lbfgs/diis/rfo)
    PARAM(diis_history, Int, 5, "DIIS: number of stored error vectors", "Native", { "diis_hist" })
    PARAM(diis_start, Int, 5, "DIIS: first iteration to apply extrapolation", "Native", {})
    PARAM(rfo_lambda, Double, 0.1, "RFO: initial trust radius", "Native", { "lambda" })
    PARAM(numgrad, Bool, false, "Use numerical gradient (finite differences, for debugging)", "Advanced", {})
    PARAM(numerical_gradient_step, Double, 1e-5, "Step size for numerical gradient [Bohr]", "Advanced", {})

    END_PARAMETER_DEFINITION
}
// ^^^^^^^^^^^^ PARAMETER DEFINITION BLOCK ^^^^^^^^^^^^
