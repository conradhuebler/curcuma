/*
 * <Shared SCF parameters of the quantum-chemical methods>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0.
 *
 * Claude Generated: registry module "scf" holds the SCF knobs that several methods
 * honour (native gfn1/gfn2, native NDDO pm3/am1/mndo/pm6, the external tblite, xtb
 * and Ulysses wrappers). Each method keeps its own default, so every default here is
 * a "use the method's own value" marker (0, -1 or "auto"). Knobs that only one solver
 * has (level shift, eigensolver, mixed precision, GPU, preconditioner) stay in module xtb.
 *
 * Resolution, highest first: the method's own scope (-xtb.scf_threshold, -tblite.max_iterations),
 * then this scope (-scf_threshold, -scf.scf_threshold), then the method default.
 * MethodFactory::applySharedScfScope() resolves aliases, drops the markers and maps the
 * values onto the external wrappers' own parameter names. See docs/SCF_MODES.md.
 */

#pragma once

#include "src/core/parameter_macros.h"

BEGIN_PARAMETER_DEFINITION(scf)
    PARAM(scf_threshold, Double, 0.0, "SCF convergence threshold; 0 uses the method's own: 1e-5 on max|dq_shell| for gfn1/gfn2, 1e-6 on the density change ||dP||_F for pm3/am1/mndo/pm6.", "SCF", {})
    PARAM(scf_max_iterations, Int, 0, "Most SCF iterations; 0 uses the method's own (gfn1/gfn2 150, pm3/am1/mndo/pm6 100, tblite/xtb/Ulysses their max_iterations). Raise it for slowly converging systems.", "SCF", {"max_iterations", "SCFmaxiter"})
    PARAM(scf_damping, Double, 0.0, "Fraction of the new charges (gfn1/gfn2, tblite) or density (pm3/am1/mndo/pm6) taken per step; 0 uses the method's own (0.4, NDDO 0.6). Lower values are more robust and slower.", "SCF", {})
    PARAM(scf_allow_unconverged, Bool, false, "Accept an SCF that reaches the iteration limit without converging (gfn1/gfn2, pm3/am1/mndo/pm6); by default this is an error. Non-finite energies are always an error.", "SCF", {})
    PARAM(scf_mode, String, "auto", "SCF strategy; 'auto' is broyden for gfn1/gfn2 and diis for pm3/am1/mndo/pm6. gfn1/gfn2 also take diis, plain, level-shift; NDDO takes plain (damping only).", "SCF", {})
    PARAM(scf_diis_start, Int, -1, "Damped iterations before Pulay DIIS starts; -1 uses the method's own (gfn1/gfn2 5, in scf_mode diis and level-shift; pm3/am1/mndo/pm6 3). The name carries the scf_ prefix because -diis_start belongs to the geometry optimiser (-sp, -opt).", "SCF", {"diis_start"})
    PARAM(scf_diis_subspace, Int, 0, "Number of Fock matrices kept in the DIIS history; 0 uses the method's own (gfn1/gfn2 6, pm3/am1/mndo/pm6 8).", "SCF", {"diis_subspace"})
END_PARAMETER_DEFINITION
