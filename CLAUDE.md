# CLAUDE.md - Curcuma Development Guide

## Overview

**Curcuma** is a molecular modelling and simulation toolkit for computational chemistry, force fields, and quantum calculations.

**Educational Purpose**: Teaching and research platform prioritizing pedagogical clarity over complex software engineering. Goal: learn, understand, and implement computational chemistry methods without getting lost in C++ abstractions.

## Very General Instructions for AI Coding
- Avoid flattery, compliments, or positive language. Be clear and concise. Do not use agreeable language to deceive.
- Do comprehensive verification before claiming completion
- Show me proof of completion, don’t just assert it
- Prioritize thoroughness over speed
- If I correct you, adapt your method for the rest of the task
- No completion claims until you can demonstrate zero remaining instances
- Dont use git -A to blindly add files

## AI-Generated Content and Validation Policy

### Status Labels — Definitions
These labels are used throughout CLAUDE.md files and documentation. Only the human operator may assign ✅ TESTED or ✅ APPROVED.

| Label | Meaning | Who sets it |
|-------|---------|-------------|
| 🤖 AI-generated | Code written by AI, not reviewed by human | AI |
| ⚙️ Machine-tested | Passes automated tests (CI, ctest) | AI |
| 👁️ Human-reviewed | Human has read and understood the code | Human only |
| ✅ TESTED | Human has run it on real problems and it behaves correctly | **Human only** |
| ✅ APPROVED | Human confirms correctness, ready for production | **Human only** |

**The AI must never write ✅ TESTED or ✅ APPROVED on its own work.**

### Conservative Self-Assessment Rules for AI
When documenting implemented features, the AI must apply these rules:

1. **Automated tests pass ≠ correct** — tests only cover what was anticipated. Unknown failure modes exist.
2. **Agreement with reference on test molecules ≠ general correctness** — the reference comparison is only as broad as the test set.
3. **No gaps visible ≠ no gaps exist** — absence of a known bug is not the same as correctness. Especially for AI-generated scientific code: the most dangerous bugs are those that produce plausible but wrong results.
4. **"Implemented" means the code compiles and runs** — it does not imply physical correctness, numerical stability across all inputs, or completeness relative to the reference method.
5. **When in doubt, add a caveat** — a caveat that turns out to be unnecessary is harmless. A missing caveat on wrong code causes user errors.

### Required Documentation for New AI-Generated Features
Every new method or capability added by AI must include in its CLAUDE.md:
- What was tested (which molecules, which conditions)
- What was **not** tested (system classes, edge cases, conditions)
- What is **not implemented** relative to the reference method
- A note that human production testing is pending until the human removes it

## General Instructions

- Each source code dir has a CLAUDE.md with basic information of the code and logic
- **Keep CLAUDE.md files FOCUSED and CONCISE** - ONE clear idea per bullet, max 1-2 lines
  - ❌ DON'T: Multi-paragraph explanations, code examples, historical details
  - ✅ DO: Brief statements with links to detailed docs if needed
  - ✅ DO: "✅ **Feature name** - Brief description" for completed items
- Remove completed/resolved items after 2-3 updates (move to git history)
- Tasks corresponding to code must be placed in the correct CLAUDE.md file
- Each CLAUDE.md has a variable part (short-term info, bugs) and preserved part (permanent knowledge)
- **Instructions blocks** contain operator-defined future tasks and visions for code development
- Only include information important for ALL subdirectories in main CLAUDE.md
- Preserve new knowledge from conversations but keep it brief
- Always suggest improvements to existing code
- **Keep entries concise and focused to save tokens**
- **Keep git commits concise and focused**
- **Rule of thumb**: If a CLAUDE.md section exceeds 20 lines, consider if it's better placed elsewhere
- Newly added features need a precise and short documentation under docs/, a link to the documentation from claude.md and a note in the readme
## Development Guidelines

### Code Organization
- Each `src/` subdirectory contains detailed CLAUDE.md documentation
- Variable sections updated regularly with short-term information
- Preserved sections contain permanent knowledge and patterns
- Instructions blocks contain operator-defined future tasks and visions

### Implementation Standards

#### Educational-First Design Principles
- **Core functionality visibility**: Always provide clear, direct access to the computational chemistry implementation
- **Minimal abstraction layers**: Avoid unnecessary templates, inheritance hierarchies, or design patterns that obscure the scientific content
- **Algorithm transparency**: The actual mathematical/physical implementation should be easily locatable and readable
- **Documentation focus**: Emphasize *what* the code does scientifically, not just *how* it's structured
- **Learning-oriented comments**: Include references to equations, papers, and theoretical background in code comments
- **Method implementation clarity**: Each computational method should have a clear entry point with minimal indirection
- **Accuracy:** 100 % with respect to referenz implementation for any scientific method

#### Code Organization for Learning
- **Flat over hierarchical**: Prefer simple, direct implementations over complex class hierarchies
- **Self-contained modules**: Each computational method should be understandable without deep knowledge of the entire system
- **Clear naming**: Function and variable names should reflect their scientific meaning
- **Minimal templates**: Only use templates when absolutely necessary; prefer explicit types for clarity
- **Direct implementations**: Avoid hiding core algorithms behind layers of abstractions

#### Standard Development Practices
- Mark new functions as "Claude Generated" for traceability
- Document new functions briefly (doxygen ready) with scientific context
- Document existing undocumented functions if appearing regularly (briefly and doxygen ready)
- Remove TODO Hashtags and text if done and approved
- Implement comprehensive error handling and logging
- Maintain backward compatibility where possible
- **Always check and consider instructions blocks** in relevant CLAUDE.md files before implementing
- Reformulate and clarify task and vision entries if not already marked as CLAUDE formatted
- In case of compiler warning for deprecated functions, replace the old function call with the new one
- Implement timing analysis for complex functions
- Keep track of significant improvements in AIChangelog.md, one line per fact
- **Complex Architecture Documentation**: Factory patterns, dispatchers, and multi-step workflows require comprehensive inline documentation following ARCHITECTURE_DOCUMENTATION.md standards
- **BMT output compatibility (MANDATORY)**: Every capability that writes output files MUST route them through `outputPath()` (CurcumaMethod subclasses) or `BMTUtils::outputPath()` (standalone handlers). Hardcoded CWD paths are not permitted. Verify with `-no_bmt` (legacy) and default BMT mode before merging.
- **No UTF symbols in terminal output**: Do not use Unicode box-drawing characters, emoji, arrows (->), checkmarks, or any non-ASCII symbols in fmt::print/std::cout output. Use plain ASCII only. Reason: breaks output in many terminal emulators, log files, and remote shells. CurcumaLogger colored output is exempt (uses ANSI codes, not Unicode).

#### Parameter Definition Standards (MANDATORY for new capabilities)
- **ALL new capabilities MUST use Parameter Registry System** - no static JSON configurations
- **Definition Location**: Define parameters in capability header using PARAM macros within `BEGIN_PARAMETER_DEFINITION(module)` block
- **Naming Convention**: Use **snake_case** exclusively (`max_iterations`, not `MaxIterations` or `maxIterations`)
- **Include Required**: Add `#include "src/core/parameter_macros.h"` to capability header
- **Help Text**: Provide comprehensive, user-facing descriptions for each parameter
- **Categories**: Group parameters logically (Basic, Algorithm, Output, Advanced)
- **Aliases**: Add old parameter names as aliases for backward compatibility during migration
- **Type Safety**: Use correct ParamType (String, Int, Double, Bool) matching C++ type
- **Constructor**: Use `ParameterRegistry::getInstance().getDefaultJson("module")` instead of static JSON
- **Build Verification**: Run `make GenerateParams` and check for validation warnings
- **Documentation**: See reference implementation in `src/capabilities/analysis.h`
- **Migration Guide**: Follow [docs/PARAMETER_MIGRATION_GUIDE.md](docs/PARAMETER_MIGRATION_GUIDE.md) for existing capabilities

#### Copyright and File Headers
- **Copyright ownership**: All copyright remains with Conrad Hübler as the project owner and AI instructor
- **Year updates**: Always update copyright year to current year when modifying files
- **Claude contributions**: Mark Claude-generated code sections but copyright stays with Conrad
- **Format**: `Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>`
- **AI acknowledgment**: Add Claude contribution notes in code comments, not copyright headers

#### Code Structure Guidelines
- **Main computational functions**: Should be easily findable and readable without deep C++ knowledge
- **Algorithm documentation**: Include mathematical formulations and literature references
- **Parameter transparency**: Make method parameters and their physical meaning obvious
- **Debugging accessibility**: Provide easy ways to inspect intermediate results and algorithm steps
- **Educational examples**: Include well-commented example usage in documentation

## Current Capabilities

### 1. Quantum Mechanical Methods

#### Native Implementations (Educational, No External Dependencies)

> ⚠️ **All native QM methods are AI-implemented and machine-tested only — not human production tested.**
> Results should be validated against external references (TBLite, Ulysses, XTB) before use in research.

- ⚠️ **Extended Hückel Theory (EHT)** - AI-implemented, machine-tested
- ⚠️ **GFN2-xTB (Native)** - AI-implemented, machine-tested; canonical `gfn2` backend; `-opt` works; Broyden SCF default (`-scf_mode diis|plain|level-shift`, `-scf_guess h0|eeq`); 11/12 sqm_reference @1e-8 vs tblite (only 231-atom `complex` open) — [docs/SQM_VALIDATION.md](docs/SQM_VALIDATION.md), [docs/SCF_MODES.md](docs/SCF_MODES.md)
  - **d-shell support (X-I1, June 2026)**: S/P/Cl/Si/… (main-group d) now compute via the cartesian→spherical transform, ≤1e-8 Eh vs tblite. CPU + **CUDA GPU** (validated on a GTX 1660: energy bit-identical to CPU, gradient ~1e-16) + **ROCm GPU** (B6 port, validated on a Radeon 890M/gfx1150: energy bit-identical to CPU at 8 dp, `-opt` tracks the CPU trajectory); **Vulkan d still falls back to CPU** (GLSL shaders pending). Transition metals enabled but **unvalidated**. See [docs/SQM_DSHELL_WP.md](docs/SQM_DSHELL_WP.md)
  - **Threading**: intra-molecule `-threads N` (setup 4×, gradient 3.6×); eigensolve capped at 8 threads (`CURCUMA_EIG_MAX_THREADS`) — [docs/SQM_THREADING.md](docs/SQM_THREADING.md)
  - **Integral setup** (Jul 2026): shell-pair-blocked overlap/multipole/gradient kernels — setup 209→91 ms, gfn2 total 1199→1083 ms on complex/231 — [docs/SQM_PERFORMANCE.md](docs/SQM_PERFORMANCE.md)
  - **Integral numerics caveat**: those kernels are algebraically exact but ~1 ulp off the old values (GCC FMA contraction); energies bit-identical, gradients ≤1.7e-14 Eh/Bohr
  - **Eigensolvers**: opt-in `-eigensolver native|purify|lobpcg`, `CURCUMA_EIG_TRED2=blocked` (MKL-free / GPU-portable) — [docs/SQM_EIGENSOLVE_GPU.md](docs/SQM_EIGENSOLVE_GPU.md)
  - **Large systems**: `-large_system_mode fragments|dc|sparse` scales SCF past ~1000 atoms — [docs/SQM_LARGE_SYSTEMS.md](docs/SQM_LARGE_SYSTEMS.md)
  - **SCF extrapolation**: `-scf_extrapolation aspc|gauss` cuts SCF iters in opt/MD (caffeine gfn2 215→90); experimental `xlbomd` = extended-Lagrangian MD — [docs/SQM_SCF_EXTRAPOLATION.md](docs/SQM_SCF_EXTRAPOLATION.md)
  - **GPU backends** (`-gpu cuda|rocm|vulkan|auto`): all three device-resident through `-opt`/`-md`; ROCm has FP32 mixed-precision ON by default (real win), Vulkan is opt-in only (eigensolve-bound, no speedup). Detail in [docs/SQM_GPU.md](docs/SQM_GPU.md) / [docs/SQM_ROCM.md](docs/SQM_ROCM.md) / [docs/SQM_VULKAN.md](docs/SQM_VULKAN.md); roadmap in [docs/SQM_GPU_ROADMAP.md](docs/SQM_GPU_ROADMAP.md)
  - **All GPU backends are runtime `dlopen` plugins** (`libcurcuma_{cuda,rocm,vulkan}.so`; CUDA Jul 2026, ROCm/Vulkan Sep 2026): the main binary has no backend symbol or `#ifdef`, CPU-only runs start in ~9 ms; `-gpu <backend>` loads the plugin lazily and warns + falls back to CPU when it is absent; `-methods` lists the plugins found. ROCm plugin build unverified (no SDK here). See [docs/GPU_PLUGIN_STARTUP.md](docs/GPU_PLUGIN_STARTUP.md)
- ⚠️ **GFN1-xTB (Native)** - AI-implemented, machine-tested; canonical `gfn1` backend; `-opt` works; Broyden SCF default; 10/12 sqm_reference @1e-8 vs tblite (He2 + `complex` remain) — [docs/SQM_WP2_gfn1_accuracy.md](docs/SQM_WP2_gfn1_accuracy.md)
- ⚠️ **PM3/AM1/MNDO (Native NDDO)** - AI-implemented, machine-tested; 21/21 tests vs Ulysses reference (< 4 µEh)
- ⚠️ **Native GFN-FF** - AI-implemented, machine-tested; see [docs/GFNFF_STATUS.md](docs/GFNFF_STATUS.md); S30L vs xtb 6.6.1 (Jul 2026): 29/30 within ~0.5 kcal/mol after F1/F2/CLI/F3/ipis + Jul 8 torsion fixes (MAD 433→1.20); residuals: 23 (.CHRG/harness quirk, curcuma judged correct), 30/AB (bond pibo) — see [docs/S30L_GFNNF_VALIDATION.md](docs/S30L_GFNNF_VALIDATION.md)

#### External Interfaces (Production Quality, Requires Compilation)
- **TBLite Interface** - Tight-binding DFT methods (GFN1, GFN2, iPEA1) + **Solvation** (CPCM, GB, ALPB)
- **XTB Interface** - Extended tight-binding methods (GFN-FF, GFN1, GFN2)
- **Ulysses Interface** - Semi-empirical methods (PM3, PM6, AM1, MNDO, RM1, etc.) + **Solvation** (GBSA)
- **Native GFN-FF** - Curcuma's own implementation (`gfnff`) - ✅ **IMPLEMENTED**

### 2. Force Field Methods
- **Universal Force Field (UFF)** - General-purpose molecular mechanics
- **GFN-FF** (`gfnff`) - ✅ **FULLY IMPLEMENTED** - See [docs/GFNFF_STATUS.md](docs/GFNFF_STATUS.md)
- ⚠️ **Coarse-grained beads** (`cg`, Sep 2026, AI/machine-tested) - LJ spheres/ellipsoids on the workspace engine, input via `-load_ff_json FILE` (`cg_default`, `cg_per_atom`, `pair_interactions`); analytic sphere gradient - see [docs/CLEANUP_2026_09.md](docs/CLEANUP_2026_09.md)
- **QMDFF** - Quantum Mechanically Derived Force Fields
- **Universal Parameter Caching** - Automatic save/load for all FF methods

### 3. Solvation Models (Implicit Solvent)
- ✅ **TBLite Solvation** - CPCM, GB (Generalized Born), ALPB for GFN methods. ⚠️ Pending in the `USE_TBLITE=OFF` dev build (see [docs/SOLVATION.md](docs/SOLVATION.md)) - unaffected: native ALPB/GBSA below.
- ✅ **Ulysses Solvation** - GBSA (Generalized Born + SA) for GFN/MNDO methods
- ⚠️ **Native GFN1/GFN2 ALPB + GBSA** (June 2026, AI/machine-tested) - self-consistent
  ALPB (`-xtb.solvent_model alpb`, P16 kernel) and GBSA (`-xtb.solvent_model gbsa`, Still kernel)
  in the native xTB SCF, matching tblite total ΔG (Born + CDS + shift; CM5 for gfn1) to
  ≤1e-8 Eh on the validation set (CPU + GPU); `-method gfn2 -xtb.solvent water -xtb.solvent_model gbsa`
  (legacy numeric codes 3/2 still accepted). CPCM native solvation still pending.
  See [docs/SQM_SOLVATION_WP.md](docs/SQM_SOLVATION_WP.md)
- ⚠️ **Native GFN-FF ALPB** (June 2026, AI/machine-tested) - self-consistent: the Born
  reaction field couples into the EEQ solve (`A_eeq += B`), so charges polarize in the solvent.
  `-method gfnff -gfnff.solvent water -gfnff.solvent_model alpb` matches **xtb 6.7.1** (`--gfnff
  --alpb`) to **≤1e-8 Eh** (7 mol × 4 solvents); analytic gradient FD-validated. GFN-FF has
  no separate GBSA (reference uses ALPB), so `-gfnff.solvent_model gbsa` maps to ALPB. See
  [docs/SQM_SOLVATION_WP.md](docs/SQM_SOLVATION_WP.md) WP5
- **25+ Solvents** - water, methanol, DMSO, acetone, benzene, etc.
- **Auto-Activation** - Specify `-solvent water` to enable
- **Documentation** - See [docs/SOLVATION.md](docs/SOLVATION.md) for details

### 4. Dispersion and Non-Covalent Corrections
- **DFT-D3** - Grimme's D3 dispersion correction
- **DFT-D4** - Next-generation D4 dispersion correction
- **H4 Correction** - Hydrogen bonding and halogen bonding corrections

### 5. Geometry Optimization
- **LBFGS Optimizer** - Limited-memory Broyden-Fletcher-Goldfarb-Shanno
- **Multiple Convergence Criteria** - Energy, gradient, RMSD-based
- **Constrained Optimization** - Distance, angle, and dihedral constraints

### 6. Conformational Analysis ✅ REFACTORED 2025
- **ConfSearch** - Systematic conformational searching (unified trajectory framework); supports **dual-method** runs (`-md_method` explore + pre-opt, `-opt_method` refine + rank; both fall back to `-method`) — see [docs/CONFSEARCH_DUAL_METHOD.md](docs/CONFSEARCH_DUAL_METHOD.md); **restartable** via `-restart` (self-contained checkpoint: bias pool + cumulative + seeds + schedule, written to CWD + BMT) — see [docs/CONFSEARCH_RESTART.md](docs/CONFSEARCH_RESTART.md); **registry-backed since Jul 2026** (67 PARAMs) so its flags are no longer auto-routed away, and every child computation (MD / 4 opt sites / 2 ConfScan passes) shares one `ChildConfig()` carrying charge, spin, gpu and the method sub-scopes; **RMSD-MTD bias speedup (Jul 2026)**: a rigorous Gaussian-cutoff screen (`-rmsd_mtd_screen`, default ON, physics-preserving) skips far hills before the Kabsch, plus an enforced pool cap (`-rmsd_mtd_max_gaussians`) — see [docs/CONFSEARCH_MTD_SCREEN.md](docs/CONFSEARCH_MTD_SCREEN.md)
- **ConfScan** - Conformational scanning along reaction coordinates
- **RMSD Analysis** - Structure comparison and alignment
- **Energy-based Filtering** - Automatic conformer ranking
- **Refactored Geometry Commands** - TrajectoryWriter for JSON format (Phase 5)

### 7. Molecular Dynamics
- **SimpleMD** - Basic molecular dynamics simulation
- ⚠️ **Temperature ramps / live T / thermal regions** (Jun 2026, AI/machine-tested) - `setTargetTemperature()` live setpoint; multi-stage `temp_ramp`/`temp_schedule` (`steps`/`reach` modes); per-atom-subset `temp_regions` (Berendsen/CSVR/Andersen; NH falls back to global). No-region path byte-identical to legacy. See [docs/TEMPERATURE_RAMP.md](docs/TEMPERATURE_RAMP.md)
- **NEB Docking** - Nudged elastic band for transition states
- **Trajectory Analysis** - Analysis of MD trajectories
- **PLUMED Metadynamics** - Enhanced sampling via PLUMED2 plugin (`-mtd` flag) — see [docs/PLUMED_HELP.md](docs/PLUMED_HELP.md)

### 8. Analysis Tools
- **✅ Parallel Analysis** - Frame-level parallelization with CxxThreadPool (3-8x speedup, January 2026)

### 9. Output Directory System
- **🤖 BMT (Basename.Method.Timestamp)** - Default output directory for all commands — see `src/tools/CLAUDE.md`
- **`-bak` flag** - Copy specified files from BMT directory back to CWD
- **`-no_bmt`** - Disable BMT, write output to CWD (legacy behavior)
- **✅ TrajectoryWriter** - Unified output system for Human/CSV/JSON/DAT formats
- **✅ Scattering Analysis** - P(q)/S(q) with logarithmic q-spacing and automatic gnuplot visualization (2026)
- **RMSD Calculations** - Root-mean-square deviation analysis
- **Persistent Diagram** - Topological data analysis
- **Hessian Analysis** - Second derivative calculations
- **Orbital Analysis** - Molecular orbital visualization and analysis

### 10. Core Computational Libraries
- ✅ **MNDO Integrals** - Dewar-Thiel multipole expansion for semi-empirical 2e⁻ integrals, see [docs/MNDO_INTEGRALS.md](docs/MNDO_INTEGRALS.md)

## Architecture

### Core Components

#### Energy Calculator (`src/core/energycalculator.h/cpp`)
**COMPLETELY REFACTORED (January 2025)** - New unified polymorphic architecture:

##### **New Architecture**
- **Polymorphic Design**: Single `ComputationalMethod` interface for all QM/MM methods
- **MethodFactory**: Priority-based method resolution with hierarchical fallbacks
- **Unified Interface**: `calculateEnergy()`, `getGradient()`, consistent API across all methods
- **Method Priority System**: `gfn2`/`gfn1` → Native xTB (AP3); `ipea1` → TBLite; `ugfn2` → Ulysses
- **Thread-Safe**: Full multi-threading support maintained
- **Universal Verbosity**: Consistent output control across all computational methods

##### **Method Resolution (New)**
```cpp
// New MethodFactory system (replaces old SwitchMethod)
std::unique_ptr<ComputationalMethod> method = 
    MethodFactory::createMethod("gfn2", config);
double energy = method->calculateEnergy();
```

##### **Supported Method Hierarchies** (AP3, April 2026)
- **gfn1/gfn2**: Native curcuma xTB (canonical); `ipea1`/`ugfn2` for other providers
- **xtb-gfn1/xtb-gfn2**: External GFN — TBLite (USE_TBLITE) → XTB binary (USE_XTB), like `xtb-gfnff`
- **tblite-gfn1/tblite-gfn2**: TBLite explicitly (forces that backend)
- **eht**: Native only (always available, no dependencies)
- **pm3**: Native only (H, C, N, O supported, no dependencies)
- **uff/qmdff**: ForceField wrapper with parameter generation
- **gfnff**: Native C++ GFN-FF (always available, ✅ **COMPLETE**)
- **xtb-gfnff**: Fortran/XTB GFN-FF — ExternalGFNFF (USE_GFNFF) → XTB (USE_XTB)

#### Force Field System (`src/core/forcefield.h/cpp`)
Modern force field engine with:
- **Threading** via the single `FFWorkspace` engine (the legacy `ForceFieldThread` engine was removed Sep 2026)
- **Universal parameter caching** - automatic save/load as JSON
- **Method-aware loading** - validates parameter compatibility
- **Multi-threading safety** - controllable caching for concurrent calculations

#### QM Interface (`src/core/qm_methods/interface/abstract_interface.h`)
Unified interface for all quantum mechanical methods:
```cpp
class QMInterface {
    virtual bool InitialiseMolecule() = 0;
    virtual double Calculation(bool gradient, bool verbose) = 0;
    virtual bool hasGradient() const = 0;
    virtual Vector Charges() const = 0;
    virtual Vector BondOrders() const = 0;
};
```

### File Organization

```
curcuma/
├── src/
│   ├── capabilities/          # High-level molecular modeling tasks
│   │   ├── confscan.cpp      # Conformational scanning
│   │   ├── confsearch.cpp    # Conformational searching  
│   │   ├── curcumaopt.cpp    # Geometry optimization
│   │   ├── simplemd.cpp      # Molecular dynamics
│   │   └── rmsd.cpp          # Structure analysis
│   ├── core/                 # Core computational engines
│   │   ├── energycalculator.cpp      # NEW: Unified polymorphic dispatcher
│   │   ├── molecule.cpp              # Molecular data structures
│   │   ├── curcuma_logger.cpp        # Universal logging system
│   │   ├── energy_calculators/       # NEW: All computational methods organized here
│   │   │   ├── computational_method.h     # Base interface for all methods
│   │   │   ├── method_factory.cpp         # Priority-based method creation
│   │   │   ├── qm_methods/                # QM method implementations & wrappers
│   │   │   │   ├── eht.cpp                # Extended Hückel Theory + verbosity
│   │   │   │   ├── xtbinterface.cpp       # XTB interface + verbosity
│   │   │   │   ├── tbliteinterface.cpp    # TBLite interface + verbosity
│   │   │   │   ├── ulyssesinterface.cpp   # Ulysses interface + verbosity
│   │   │   │   ├── gfnff_method.cpp         # ComputationalMethod wrapper
│   │   │   │   ├── orcainterface.cpp      # ORCA interface
│   │   │   │   ├── dftd3interface.cpp     # DFT-D3 dispersion corrections
│   │   │   │   ├── dftd4interface.cpp     # DFT-D4 dispersion corrections
│   │   │   │   ├── *_method.cpp           # Polymorphic method wrappers
│   │   │   │   └── interface/             # Abstract interfaces
│   │   │   └── ff_methods/                # Force field implementations
│   │   │       ├── forcefield.cpp         # Force field engine + verbosity  
│   │   │       ├── forcefieldgenerator.cpp # Parameter generation + verbosity
│   │   │       ├── ff_workspace*.cpp      # Partitioned energy/gradient engine (all FF methods)
│   │   │       ├── gfnff_method.cpp      # Native GFN-FF implementation (4329 lines)
│   │   │       ├── gfnff.h               # GFN-FF class interface
│   │   │       ├── gfnff_inversions.cpp   # GFN-FF inversion terms
│   │   │       ├── gfnff_torsions.cpp     # GFN-FF torsion terms
│   │   │       ├── qmdff.cpp              # QMDFF implementation
│   │   │       ├── eigen_uff.cpp          # UFF implementation
│   │   │       └── *_par.h                # Parameter databases
│   ├── tools/                # Utilities and file I/O
│   │   ├── formats.h         # File format handling (XYZ, MOL2, SDF)
│   │   └── geometry.h        # Geometric calculations
│   └── helpers/              # Development and testing tools
├── test_cases/               # Validation and benchmark molecules
├── external/                 # Third-party dependencies
└── CMakeLists.txt           # Build configuration
```

## Completed Developments (2026)

✅ **GFN-FF / GFN1 / GFN2 cleanup + speedup** (Sep 2026, AI/machine-tested) - one FF engine
(`FFWorkspace`; legacy `ForceField`/`ForceFieldThread` GFN-FF path and its per-step feed
removed), native xTB decoupled from `QMDriver`, duplicate parameter headers / old NDDO classes /
dead EEQ paths deleted (~20k lines), table-driven `MethodFactory` registry + one shared
sub-scope list, exact hot-path fixes (LAPACK scratch, HB-gradient index, EEQ views). Energies
identical to the last digit; GFN-FF SP 1.35-1.9x, GFN1 polymer 1.45x. Numbers, method and
open items in [docs/CLEANUP_2026_09.md](docs/CLEANUP_2026_09.md)

> Older 2025 work (parameter registry, polymorphic EnergyCalculator, native GFN2/GFN1/PM3, MNDO integrals, GFN-FF full implementation, scattering, analysis parallelization, dependency gating) is in `AIChangelog.md` + git history.

✅ **`-interaction` capability** (June 2026) - supramolecular interaction energy `E(AB)−E(A)−E(B)` for the S30L host-guest set; modes: S30L A/B/AB dir (+`.CHRG`), batch vs `reference_s30l` (MAD/RMSD), explicit `-fragA/-fragB`, single-AB auto-split
✅ **GFN-FF aromatic ring torsions fixed** (June 2026) - acyclic-only pi-sp3 rules were wrongly applied to ring torsions; gated on `!in_ring`; S30L host A now bit-identical to Fortran, validation 18/18 — see [docs/GFNFF_STATUS.md](docs/GFNFF_STATUS.md)
✅ **GFN-FF GPU HB-freeze resolved + per-frame gradient diagnostic** (June 2026) - the GPU HB-charge freeze is correct (matches CPU+Fortran); `test_gfnff_grad_traj` is the clean force metric (MD heat-exchange is not) — see [docs/GPU_GFNNF_DISCREPANCIES.md](docs/GPU_GFNNF_DISCREPANCIES.md)

## Build and Test Commands

**ALWAYS build and test in the `release/` directory.** This is the canonical build for regression comparisons. The `build/` directory is for development experiments only.

```bash
# Build — always use release/
cd release
make -j4

# Run all tests
ctest --output-on-failure

# Run specific test categories
ctest -R "cli_rmsd_" --output-on-failure      # RMSD CLI tests (6/6 passing)
ctest -R "cli_confscan_" --output-on-failure  # ConfScan CLI tests (7/7 passing)
ctest -R "cli_simplemd_" --output-on-failure  # SimpleMD CLI tests (7/7 passing)
ctest -R "cli_curcumaopt_" --output-on-failure # Opt CLI tests (6/6 passing)

# Run individual CLI test with verbose output
ctest -R "cli_rmsd_01" --verbose

# Legacy: Manual curcuma execution
./curcuma -sp input.xyz -method uff           # UFF single point
./curcuma -rmsd ref.xyz target.xyz            # RMSD calculation
./curcuma -opt input.xyz -method gfn2         # GFN2 optimization
```

**Test Status**: 26/26 CLI Tests passing (100%) ✅

## Project Management

- **Prioritized TODO List**: See [TODO.md](TODO.md)
- **Module Docs**: Each `src/` subdirectory has CLAUDE.md with specific tasks
- **GFN-FF Status**: See [docs/GFNFF_STATUS.md](docs/GFNFF_STATUS.md) for implementation details
- **Technical Debt**: Identified debt in GFN-FF / native xTB / QM interfaces / EnergyCalculator / RMSD-alignment API — [docs/TECHNICAL_DEBT.md](docs/TECHNICAL_DEBT.md)
- **rev-gfnff backlog**: places where curcuma's GFN-FF could be *better* than the reference, each one checked against an external reference (r2SCAN-3c / GFN2 / DLPNO / experiment) rather than against xtb — [docs/REV_GFNFF_TODO.md](docs/REV_GFNFF_TODO.md). Port fidelity stays the default; this collects the deliberate deviations and the open method limitations.
- **Benchmark Test-Set Retrieval** (Sep 2026): `scripts/fetch_testset.py` fetches Grimme-group sets (MOR41 auto, GMTKN55 auto via git, S30L manual/paywalled) on demand into the layout `scripts/mor41_validation.py`/`scripts/s30l_*.py` expect; `scripts/gmtkn55_compare.py` runs the same curcuma-vs-xtb per-structure check across GMTKN55's 54 subsets — gfn2 MAD 0.000, gfn1 MAD 0.047, gfnff MAD **2.117** kcal/mol (isolated-ion EEQ bugs Known Issue #8, nitro/N-H bond bugs Known Issue #12); skips open-shell for gfn1/gfn2 (Known Issues #9) — see [docs/GMTKN55_VALIDATION.md](docs/GMTKN55_VALIDATION.md). **Caches energies in `_run/energies.json` and reuses them unless `--recompute` is passed** — after a code change, drop the `<subset>/<name>|cur|<method>` keys first, or the "comparison" silently re-reports the old numbers; `scripts/testset_perf.py` benchmarks CPU/threading/GPU wall-clock on whatever set is fetched — see [docs/TESTSET_RETRIEVAL.md](docs/TESTSET_RETRIEVAL.md)
- **S30L-CI dataset** (Sep 2026, manually supplied — not in `fetch_testset.py`'s registry): 30 host-guest complexes, the counterion variant of S30L (charged guests 23-30 get an explicit counterion atom instead of a bare net charge). `test_cases/s30lci_test_set/` (gitignored per-structure dirs, tracked `reference_s30lci`/`README`), run via `scripts/s30lci_gfnff_compare.py`. Result: curcuma vs xtb MAD **0.386** kcal/mol over all 30 structures (was 1.92 before the pyrrole-Hückel fix, Known Issue #10, and the amide-H chi fix, Known Issue #11). Per-**term** against pprcht over all 90 fragments (A/B/AB): MAD **0.0007**, max 0.047 kcal/mol, **0/90 above 0.05** (was 15/90, max 6.91) after Known Issue #12.

## Workflow States
- **ADD**: Features to be added
- **WIP**: Currently being worked on
- **ADDED**: Basically implemented
- **TESTED**: Works (by operator feedback)
- **APPROVED**: Move to changelog, remove from CLAUDE.md

### Documentation Update Rules
- **Replace debugging details with architecture decisions** when issues are resolved
- **Remove unnecessary pointer addresses and crash investigation specifics**
- **Focus on architectural clarity** rather than technical debugging information
- **Document the "why" behind design decisions** for future reference
- **Eliminate redundant information** that doesn't add architectural value
- **Prioritize clean, maintainable documentation** over verbose troubleshooting history

## Git Best Practices
- **Only commit source files**: Use `git add <file>` for specific files, never `git add -A` without review
- **Review before committing**: Always check `git diff` and `git status` to avoid accidental commits
- **Build before commit**: Ensure `make -j4` succeeds and no compiler warnings/errors exist
- **Commit message format**: Start with action verb (Fix, Add, Improve, Refactor), follow with brief description
- **Include Co-Author info**: All commits include Claude contribution notes with proper attribution
- **Test artifacts stay local**: Build outputs and temporary test files are ignored by .gitignore
- **Branch names**: `feature/<topic>` for new capabilities, `fix/<topic>` for bug fixes; `<topic>` is 2-4 lowercase ASCII words in kebab-case naming the subject, no dates or issue numbers (e.g. `feature/gfnff-solvation`, `fix/bmt-dir-collision`)
- **Why the scheme matters**: `.github/workflows/ccpp.yml` builds `feature/**` and `fix/**` automatically and publishes each as its own `ci-<branch>` prerelease — a branch outside the scheme has to be added to the workflow by hand

## Standards

### Universal Logging System (`src/core/curcuma_logger.h/.cpp`)
**✅ FULLY IMPLEMENTED** across all computational methods

#### **Verbosity Levels**
- **Level 0**: **Silent Mode** - Zero output (critical for optimization/MD)
- **Level 1**: **Minimal Results** - Final energies, convergence status
- **Level 2**: **Scientific Analysis** - HOMO/LUMO, energy decomposition, molecular properties  
- **Level 3**: **Complete Debug** - Full orbital listings, timing, algorithm details

#### **Color Scheme & Functions**
- `CurcumaLogger::error()` - Always visible (red)
- `CurcumaLogger::warn()` - Level ≥1 (orange)
- `CurcumaLogger::success()` - Level ≥1 (green)
- `CurcumaLogger::result()` - Level ≥1 (white, neutral reporting of scientific results)
- `CurcumaLogger::info()` - Level ≥2 (default)
- `CurcumaLogger::param()` - Level ≥2 (blue, structured output)
- `CurcumaLogger::energy_abs()` - Energy output with units
- `CurcumaLogger::citation()` - Level ≥2 (green)

#### **Implementation Status**
- **QM Methods**: EHT, XTB, TBLite, Ulysses ✅
- **Force Fields**: ForceField, ForceFieldGenerator ✅  
- **Native Libraries**: XTB/TBLite verbosity synchronized ✅
- **Thread Safety**: Zero overhead at Level 0 ✅

### Unit System (`src/core/units.h`)
- **Centralized**: All constants in `CurcumaUnit` namespace with CODATA-2018 values
- **Internal**: Atomic units (Hartree, Bohr, atomic time) 
- **Output**: Auto-select user-friendly units (kJ/mol, Å, fs)
- **Educational**: Clear naming and comprehensive documentation
- **Migration**: Replace scattered constants with centralized functions

### JSON Controller System (`src/main.cpp`, CLI2Json)
- **Consistent Parameter Passing**: All methods use `controller["methodname"]` subdocuments
- **Examples**: `Hessian(controller["hessian"])`, `QMDFFFit(controller["qmdfffit"])`, `ModernOptimizer(..., controller["opt"])`
- **CLI Arguments**: Automatically split into controller subdocuments via `CLI2Json()`
- **Structure**: `controller[keyword][parameter]` - e.g. `controller["opt"]["verbosity"]`, `controller["hessian"]["MaxIter"]`
- **Global Parameters**: `verbosity`, `threads`, `method`, `gpu` are additionally duplicated at top-level
- **Flat-flag auto-routing (2026)**: Any registered PARAM is reachable by its flat name (`-cn_cutoff_bohr 5.5` routes to `controller["gfnff"]["cn_cutoff_bohr"]` because the registry records ownership). Same-name in the active command's module wins; truly ambiguous names (multiple owners, none matching) warn and stay in the command module. Dotted form `-<module>.<param>` always works for disambiguation. Unregistered/legacy flags stay in the command module (unchanged).
- **JSON round-trip (2026)**: `-export_run file.json` writes the resolved controller plus `_command`, `_input`, and full registry defaults for every touched module. `-import_config file.json` performs a recursive deep merge (CLI wins at every depth). Invoking `curcuma -import_config run.json` reads `_command`/`_input` from the JSON, so the file alone is enough to replay a run. See [docs/CLI_ROUND_TRIP.md](docs/CLI_ROUND_TRIP.md).

## Planned Development

### TRAJECTORY ANALYSIS CONSOLIDATION
**Status**: ✅ Phases 1-3 complete (Jan 2026) — see [docs/ANALYSIS_CONSOLIDATION_PLAN.md](docs/ANALYSIS_CONSOLIDATION_PLAN.md)
- ✅ TrajectoryWriter, analysis.cpp migration, TrajectoryStatistics extended
- ⏳ Phase 4: Migrate `trajectoryanalysis.cpp` + `rmsdtraj.cpp` (optional)
- ⏳ Phase 5-6: Cleanup geometry commands, ProgressTracker (optional)

---

### Breaking Changes (Test-Driven)
- **Molecule data structure refactoring**: Hybrid SOA/AOS design for better performance
  - **PHASE 1**: ✅ Comprehensive test suite with refactoring-specific validation
    - `src/core/test_molecule.cpp`: 15 test categories covering all functionality
    - `src/core/REFACTORING_ROADMAP.md`: Detailed phase-by-phase plan
    - Tests include current behavior AND validation for planned improvements
    - Specific tests for: XYZ parser unification, cache granularity, fragment O(1) lookup, type safety
  - **PHASE 2**: XYZ Comment Parser unification (eliminate 10 duplicate functions)
    - **CRITICAL**: Production comment formats must not break (ORCA, XTB, simple energy)
    - See `XYZ_COMMENT_FORMATS.md` for required format compatibility
  - **PHASE 3**: Granular cache system (replace single m_dirty flag)
  - **PHASE 4**: Fragment system O(1) lookups (replace std::map)
  - **PHASE 5**: Type-safe ElementType enum (replace int elements)
  - **PHASE 6**: Unified atom structure with zero-copy geometry access
  - **CRITICAL**: All existing functionality must remain API-compatible

## Known Issues

1. **GFN-FF Limitations**: See [docs/GFNFF_STATUS.md](docs/GFNFF_STATUS.md#known-limitations) for details (D4 dispersion, EEQ integration, metal parameters)

2. **GFN-FF S30L validation (Jul 2026, AI/machine-tested; updated Jul 8, 2026)**: see [docs/S30L_GFNNF_VALIDATION.md](docs/S30L_GFNNF_VALIDATION.md) — vs xtb 6.6.1: **29/30 within ~0.5 kcal/mol (MAD 433→1.20)**, only 23 and 30 outside — a later torsion-term fix (Jul 8, 2026) resolved the 27/28+7/8 residual this entry originally reported (was 28/30, MAD 1.38, four residual groups; superseded). Fixed: (a) charged complexes (23-30) EEQ qfrag=[0,0] → xtb-style both-assignment trial (F2); (b) `-charge -N` parsed as +N → negative-number CLI fix; (c) `AngleBending` 1/sinθ NaN guard (F1); (d) F2 cache bug — inline topo-cache load dropped qfrag (cached charged nfrag==2 re-runs neutralised Coulomb; now restores qfrag, write guard checks SUM); (e) F3 bond/Hückel/Coulomb — S CN=2 sp2→sp3 (11/12), halogen hyb + hoffdiag default 0 + π-system membership aligned to xtb (F/Cl in, Br/I out), `metal_type[86]` array had 83 entries (Kr/I/Xe missing → I read as TM → dgam ff=-0.9 not -0.07 → 15/16 Coulomb -14 kcal); array rewritten to xtb metal(86), 15/16 fixed; (f) F3d Hückel ipis — xtb subtracts the π-system charge from nelpi (`nelpi -= ipis`); curcuma didn't → charged hosts had nelpi too large → wrong pibo → bond off; now computes ipis (qheavy + neutralize-fragment + re-EEQ + dqa·1.1) and subtracts; fixed 25/26 (→0.0) and 30/B. Then (g) a torsion-term fix (Jul 8, 2026, see [docs/S30L_GFNNF_VALIDATION.md](docs/S30L_GFNNF_VALIDATION.md)) resolved the 27/28 (torsion deficit) and 7/8 residuals (cur ~80%→~100% of xtb). **Current residuals: only 23 (.CHRG/harness quirk, curcuma correct) and 30/AB (bond pibo, ipis=0).** Pre-existing from merge: `cli_curcumaopt_07_opt_multixyz` golden-value drift.

3. **Verbosity scoping across threads (residual)**: the global `CurcumaLogger` verbosity is now RAII-scoped via `CurcumaMethod` (ctor/dtor save-restore), but it is a shared static and cannot be cleanly scoped across `CxxThreadPool` workers — the pool-owning helpers (`PerformMolecularDynamics`/`PerformOptimisation`) and energy-method setup re-assert the level at their boundaries (kept deliberately); a `thread_local` verbosity (the only full fix) is out of scope. See [docs/CONFSEARCH_ROADMAP.md](docs/CONFSEARCH_ROADMAP.md) #1.

4. **ConfSearch Phase A-C**: efficiency/robustness features (RATTLE threshold, topo/Epot abort, seed funnel, opt→bias feedback, permutation-aware + adaptive MTD bias) — roadmap, open TODOs and experimental caveats in [docs/CONFSEARCH_ROADMAP.md](docs/CONFSEARCH_ROADMAP.md). Cross-run bias heating (shared-pool hills `W=k·counter` grow unbounded → `<T>` climbs run-by-run → NaN) is bounded by **defaults ON for ConfSearch**: `rmsd_mtd_freeze_inherited`+`temp_abort` (measured best: 0 blow-ups, best conformer yield; `rmsd_mtd_max_height` is opt-in for tighter T). The bare `-startT 500` run no longer blows up (roadmap TODO #4; intra-run wide-hill blow-up still open).

5. **Native GFN transition metals FIXED (MOR41, Jul 16, 2026, AI/machine-tested)**: three bugs, all in `xtb_native.cpp`/`STO_CGTO.hpp`. (a) `reference_occ`/`p_kcn`/`p_shpoly` are angular-momentum-indexed in tblite but were read by shell index; TMs order shells `[d,s,p]` so this scrambled them (main-group unaffected, ang==pos). GFN2: index by angular momentum. GFN1: same, but valence-aware (its valence+polarisation shells share an l — only H's two s-shells — so occupation goes to the first shell of each l; `p_kcn` stays shell-indexed per tblite gfn1). (b) STO-NG tables stopped at n=5, so 5d metals' 6s/6p shells used wrong Gaussians (~0.1–0.5 Eh); added tblite's dedicated 6s/6p STO-6G arrays. Further fixes closed the residual, all in the overlap/dispersion (electronic params were already exact): (c) **5p STO-4G transcription error** in `STO_CGTO.hpp` (pAlpha4 5p 3rd primitive) corrupted every 5p overlap (I + 4d/5d); (d) **D4 `r4/r2` table truncated at Z≤36** (placeholder 10.0) → heavy-element D4 under-binding; (e) **D3/D4 CN-Gaussian weight-normalization threshold** (`sum_weights>1e-10`) wrongly collapsed C6 for sparse-reference-CN metals → GFN1 D3 heavy-metal under-binding; (f) **missing D4 `sscale` entry** for reference-system Na (refsys=11) in `d4_corrections_data.cpp` → the 4th (high-CN) reference of Sc/Ti/V/Zr/Nb/Hf/Ta had C6 ~2.4× too large (gfn2 PR40/Ti over-binding). **Final: native GFN1 and GFN2 reproduce tblite for ALL 95 MOR41 structures to <1e-6 Eh (95/95 both methods; ~85/95 at the 1e-8 print floor).** MOR41 reaction MAD native-GFN2-vs-xtb **507→0.03 kcal/mol**. No main-group/3d/GFN-FF regression; all energy ctests pass. **Separate/open**: GFN-FF transition metals. The metal `btyp>=5` bond branch **is now implemented** (`53d6aeb`/`14fc648`); with the four-list neighbour port MOR41 GFN-FF per-structure MAD is **7.30** kcal/mol (max **37.80** = ED07, 39/95 within 1). Still open: `btyp=6` (eta) promotion unwired, TM-TM `mchar` attenuation omitted, and **ED07** (largest residual, entirely in the bond term, no eta ligand — unexplained). [docs/GFNFF_METAL_BOND_ANALYSIS.md](docs/GFNFF_METAL_BOND_ANALYSIS.md) is the **pre-fix** diagnostic (flagged superseded). Runner `scripts/mor41_validation.py`, results [docs/MOR41_VALIDATION.md](docs/MOR41_VALIDATION.md).

6. **GFN-FF FT-HMO π-occupation fixes + reference split (Jul 23, 2026, AI/machine-tested)**: two faithful ports in `huckel_solver.cpp` fixed the metal-coordinated aromatic rings — (a) **open-shell `occu` split for odd nelpi** (was closed-shell `nel/2` doubled → lost the odd electron + wrong biradical index; now `ihomoa=nel/2+1`/`ihomob=nel/2`, two `fermismear` passes; even nelpi bit-identical → COT/CB untouched) fixed the 5e Cp rings; (b) **`pisip>0.40` "wrong pi occupation" fallback** (`gfnff_ini.f90:1082`, xtb variant, NOT print-gated, redo `nelpi-1` at et=4000) fixed the 5-atom/7e N-heteroaromatics (ED21/PR16/ED16a: +180→~+2 kcal). Then (c) **carbene itag→FT-HMO** (the `calculatePiBondOrders` call site hard-coded an all-zeros itag → 2-coordinate carbene C mis-counted +1 π-e) and (d) **carbene angle θ0=145°** (`gfnff_ini.f90:1573`, was missing → 5-ring carbene kept θ0=109°) fixed ED16b (organic amidine) −20.9→**0.0000** kcal. Then (e) **halogen-bond B-atom topological filter** (`detectHalogenBondsNative`, `gfnff_method.cpp:9310`) — only excluded B directly bonded to X; Fortran `gfnff_ini.f90:872` needs `bpair(B,X)>3` (B must be A…B, not X-B). S/P/metal donors admitted B atoms 2-3 bonds away → ~11× too many X-bonds (PR34 114/−0.0201 vs ref 22/−0.0018 Eh). Now filters on `topo_info.bpair[X][B]<=3`; PR34 −10.5→+1.0, 31 structures improved. Then (f) **EEQ `gam`/`chi` heavy-element array corruption** (`gfnff_par.h`) — `gam_eeq` was placeholder garbage for **all Z=56-86** (W gam +0.064240 vs Fortran −0.003724) and `chi_eeq` wrong for Z=57-71 (La-Lu); replaced with the verbatim Fortran `*_angewChem2020` arrays. Wrong W hardness → W EEQ charge 0.326 (ref 0.351) → wrong Coulomb + metal-bond fqq. Fixed **every 5d-metal (W/Ir/Pt) at once: ED07 +8.5→+0.09, PR07 +6.7→+0.10, PR22 +2.9→+0.23, ED18/ED22/PR31/…** (14 improved). Then (g) **eta-aware X-bond bpair** (`detectHalogenBondsNative`) — the Fortran `bpair` (nbondmat/pairsbond) needs SYMMETRIC reachability and eta bonds are stored asymmetrically (metal lists the eta-C, eta-C omits the metal), so eta never bridges; curcuma's plain-BFS bpair DID bridge through the eta Ru-C, shortcutting X…B (ED33 P-Ru-C_eta=2 vs ref 5) and dropping the valid far X-bond. The filter now rebuilds the distance on an eta-free adjacency (metal↔itag=−1 edges removed; normal Ru-P/Ru-S kept, PR34 unchanged); ED33 +6.8→+0.45. Then (h) **SP3-specials torsion order** (`gfnff_torsions.cpp`) — the Fortran `gfnff_ini.f90:1746` SP3-specials block (sp3 group5-group5 N-N/P-P/N-P → nrot=3, phi0=60, f1=3.0, raw hyb) comes AFTER and overrides the pi-sp3 case (:1733); curcuma applied it BEFORE its (relocated-to-end) pi-sp3 override, so aminophosphine P-N torsions (N in a pi ring) were reset to phi0=180/f1=0.2 → half the torsion energy. Moved the SP3-specials override after the pi-sp3 block; ED30 −3.44→+0.00, PR30 −3.72→+0.08. Then (i) **bond `fcn` heavy-atom neighbour count** (`gfnff_method.cpp:~4754`, `getGFNFFBondParameters`) — the Fortran heavy-heavy bond weakener `fcn = 1/(1+0.007·nb(20,i)²)/(1+0.007·nb(20,j)²)` (`gfnff_ini.f90:1181-1183`) uses `topo%nb(20,i)`, the **bonded-neighbour COUNT** (slot 20 of the `nb` array holds the degree). curcuma called `countNeighborsWithin20Bohr()` — a literal 20-Bohr distance sphere (~40 atoms in a compact complex) — collapsing fcn to ~0.007 and nuking every **non-metal** heavy-heavy bond (P-P, P-S, S-S, …; metal bonds were unaffected — the metal branch `gfnff_ini.f90:1254-1259` already used the bonded degree via `topo.neighbor_lists`). PR27 has a spurious cis P…P bond (R=2.645 Å, both P on the Ru): its fc was −0.0002 vs reference −0.046, losing −18.5 kcal — the **entire** +17.7 kcal PR27 residual. Now uses `topo.neighbor_lists[atom].size()` (== `nb(20,i)`, consistent with the metal branch); PR27 +17.67→−0.76, P-P bond E −0.00015→−0.02961 Eh (ref −0.02963). Only PR27 changed (sole MOR41 non-metal heavy-heavy bond); 35/35 golden-value ctests byte-identical. Then (j) **FT-HMO π-occupation second-attempt temperature** (`huckel_solver.cpp`) — the "probably wrong pi occupation" fallback (`pisip>0.40` → redo with `nelpi−1`) was re-solving the reduced-electron Hückel Hamiltonian at **et=4000** (the xtb variant); the pprcht/gfnff reference (`gfnff_ini.f90:993`) deliberately uses **et=300** for the second attempt (colder → sharper occupation; only the FIRST solve uses 4000). Wrong et left a systematic ~0.005−0.008 piBO offset on the metal-coordinated N-heteroaromatic fragments, entirely in the bond term (both r0 via `pi_shift` and fc via `fpi` scale with piBO). Threaded an `et` arg through `solveAndBuildDensity` (default 4000, main solve unchanged); the redo passes 300.0. Fallback still always fires (xtb-like, NOT print-gated like pprcht's `if(pr2)` — else the energy would depend on the verbosity flag). piBO now bit-identical (<1e-6) to the reference on these fragments; PR25 +2.87→+0.10, ED25 +2.77→+0.03, PR16 +1.95→−0.00, ED16a +1.02→+0.04, ED21/PR21 +2.3/+2.5→+0.3/+0.5 (the 6 fallback-firing structures; zero regressions). Then (k) **CN neighbour-list cutoff too tight** (`gfnff.h` PARAM `cn_cutoff_bohr` + fallbacks) — the dynamic bond r0 = `(r0_base+cnfak·CN)·ff` needs the reference's CN, and the reference `cnthr` (`gfnff_param.f90:551`, accuracy=1) is **100 Bohr² = 10 Bohr**, but curcuma's default was **6.0 Bohr** — TIGHTER than the reference. Heavy/metal atoms have large covalent radii, so their erf-CN transition extends past 6 Bohr; the tight cutoff truncated real CN contributions. The FFWorkspace energy path (`calculateGFNFFCNWithNeighbors`) uses this cutoff, so its dynamic-r0 CN was wrong (PR23 Ir: CN(P) 3.491 vs correct 3.621 → bond +1.63 kcal). The legacy per-bond thread path uses the full O(N²) `calculateGFNFFCN` and was always correct — that mismatch (FF params right, workspace CN truncated) is how it surfaced. Raised the default to 10.0 (CN converged: 10==20==40); PR23 −4.083676→−4.086275 (analyzer −4.086276, 0.0003 kcal). Then (l) **BATM 1,4-pair test not eta-aware** (`calculateTopologyInfoOnce`, gfnff_method.cpp) — the bonded-ATM triple list is built from 1,4 pairs (`bpair==3`) + neighbours (`gfnff_ini.f90:708`); curcuma tested `bpair` against the plain-BFS `topo_distances`, which bridges through eta Ru-C bonds. Same asymmetric-eta issue as the X-bond bpair (k'): the Fortran `topo%bpair` (nbondmat/pairsbond) needs SYMMETRIC reachability so eta never bridges. The BFS shortcut generated hundreds of spurious triples (PR28 1085 vs reference 656) → BATM over-binding ~1 kcal on the eta complexes (PR26/PR28/PR27/ED33 — the whole remaining residual). Now uses the eta-free adjacency (metal↔itag=−1 edges removed) for the `bpair==3` test; triple counts match exactly (650/656/645/908), the 4 totals match the analyzer to **<0.002 kcal**. Then (m) **torsion periodicity for metal central bonds** (`gfnff_torsions.cpp`) — the Fortran nrot (`gfnff_ini.f90:1730-1744`) is keyed on the central BOND type and pi membership, not raw atom hyb: `nrot=1` default, `=3` if both central atoms hyb==3 (Me case), `=2` if `btyp(m)==2` (pi bond), `=3` if pi-sp3 (the sp2 atom has `piadr>0`). curcuma inferred nrot from hyb_j/hyb_k alone, so metal bonds (btyp=5) were mis-assigned: a TM with 3 neighbours is hyb=2 (`gfnff_ini2.f90:241`), so its M-C(sp3) torsion took the sp2-sp3 branch (nrot=3 + pi-sp3 f1=0.5, half barrier) and its M-C(sp2) torsion took sp2-sp2 (nrot=2) — but a metal atom has piadr=0 and a metal bond is btyp=5≠2, so the reference keeps nrot=1, f1=1.0. Left the torsion term ~0.4-0.56 kcal low (ED40a/PR41/ED14 — Ru/Ti/Ni). Gated the sp2-sp2→nrot=2 branch on `bond_type==2` and the sp2-sp3→nrot=3 periodicity + pi-sp3 f1=0.5 on the sp2 atom being pi; the 3 match the analyzer to <0.001 kcal. Then (n) **misaligned HB basicity/acidity arrays** (`gfnff_par.h`) — `hb_basicity`(xhbas)/`hb_acidity`(xhaci) had 17 zeros instead of 15 on the Ar-Ge (Z18-32) filler row, shifting every later entry by two, so As/Se/Br/Sb/Te/I read the wrong HB params (iodine got hbbas=3.5/hbaci=1.5 instead of 1.9/2.50). Surfaced on ED11 (Pd/I complex): residual entirely in the HB term, the C-H…I(iodide) bond; per-triple dump vs `abhgfnff_eg1` showed rdamp/qhoutl matched but bas/aci didn't → wrong iodine basicity/acidity. Removed the 2 extra zeros; ED11 −0.174→−0.001, PR35 (Se) +0.053→+0.014, organic HB unchanged. vs the **port reference pprcht/gfnff** (built standalone at `external/gfnff`, `-Dbuild_exe=ON`): per-structure MAD **7.27→…→0.023→0.007→0.004** (within-1 52→**95/95**, **all within 0.1 kcal**), **reaction-level MAD 0.008 kcal/mol** (max 0.07); 71/71 runnable gfnff ctests pass, 35/35 golden values unchanged; 3d/4d metals (Z≤55) unchanged. Fixes (i)-(n) are **per-structure** correctness, not reaction-level cancellation. **OPEN PORTING QUESTION — reference split**: the two GFN-FF impls themselves disagree on TMs — **pprcht/gfnff vs xtb 6.7.1 MAD 11.6, max 178 kcal** (PR10/PR04/ED10/…). curcuma now tracks pprcht (its port source) to MAD 1.4; the large "vs-xtb" residual is that divergence, not a curcuma bug. xtb changed TMs since the pprcht snapshot (cf. Moradi et al., JCC 2026, TM xTB extensions). **RESOLVED vs DLPNO-CCSD(T) (Table S1, 41 reactions): pprcht MAD 62.6, curcuma 63.1, xtb 71.7 kcal — pprcht (which curcuma tracks) is closer to the true QM than xtb (6/8 of the most-divergent reactions), so no reason to re-target xtb. But ALL GFN-FF impls are ~60-70 kcal from DLPNO: GFN-FF is a force field, fundamentally unsuitable for MOR41 reaction thermochemistry (a method limitation, not a port bug; GFN2 reaches ~12 MAD).** Pure-port residuals remaining (pprcht==xtb) are all small, all in the **metal-bond term**, and were investigated to the precision floor — **no clean single-parameter fix exists**: the `bond_params` array matches the Fortran `bond_angewChem2020` exactly (Z 1-86, 0 mismatches), and per bond the `fheavy`/`fpi`/`fqq`/`fcn`/`bstrength` factors all match. The residual is a sub-0.005 Å difference in the metal-bond equilibrium r0 = `gfnffrab(CN)` + metal shifts — a CN/geometry-dependent quantity, scattered across element pairs (PR23 Ir: Ir-C 0, Ir-P +0.005, Ir-Cl +0.002, Ir-I +0.001 Å) and below the analyzer's print precision, that accumulates to ED33/PR23/PR35 residuals. This was originally (wrongly) called "irreducible metal-bond r0 fine-precision"; the apples-to-apples egbond decomposition (instrumenting the Fortran energy-time dynamic r0, which differs from the printed setup r0) later proved the metal bonds are bit-identical at the energy level, and the real cause was the CN cutoff (k). After the fcn (i), FT-HMO et (j), CN-cutoff (k), BATM eta-aware (l), metal-torsion (m) and HB-array (n) fixes the port is at per-structure MAD **0.004**, reaction MAD **0.008** vs pprcht/gfnff, **all 95 within 0.1 kcal**; worst per-structure now PR30/PR31 ~0.07 kcal (residual fine-precision). Every MOR41 GFN-FF residual >0.1 kcal has been traced to a specific bug and fixed. Analyzer for per-bond/angle/torsion ground-truth: `external/gfnff/_build/gfnff <xyz>` (prints pibo/fqq/fc + angle + torsion tables with `pr=.true.`). See [docs/MOR41_VALIDATION.md](docs/MOR41_VALIDATION.md).

7. **Cross-platform (Wine vs. native Windows) GFN-FF determinism (Aug 2026, AI/machine-tested)**: same `curcuma.exe` gave a different bond term (MOR-testset PR38) run natively on Windows vs. under Wine on Linux. Root cause: `erf`/`acos`/`exp`/`log` are dynamically imported from `api-ms-win-crt-math-l1-1-0.dll`; Wine's reimplementation of that DLL rounds them differently in the last bit than native Windows' `ucrtbase`, which can flip GFN-FF/EEQ's hard-coded classification thresholds (CN cutoffs, 150°/170°/~350° angle thresholds) into a different discrete decision. Fixed with self-contained fdlibm-derived replacements behind `-DUSE_PORTABLE_MATH=ON` (off by default; on for the Windows nightly build). `round()`/`lround()` deliberately not vendored (IEEE 754 makes them exact regardless of input). See [docs/PORTABLE_ERF.md](docs/PORTABLE_ERF.md). **Not yet confirmed on a real Windows machine** (none available) — only that the identified call sites are now architecturally immune to CRT-DLL divergence.

8. **GFN-FF: two isolated-ion EEQ bugs FIXED (Sep 2026, AI/machine-tested, found+fixed via GMTKN55)**: (a) native GFN-FF returned exactly 0.0 Eh for any single free (charged) atom — the Coulomb/electrostatics term was built from a pairwise list, empty for N=1, so the per-atom diagonal EEQ self-energy (`chi*q + 0.5*gamma_AA*q^2`, nonzero even for one atom) was never added (xtb Mg2+: +3.8766 Eh vs curcuma 0.0 Eh — matched the 2432.6 kcal/mol GMTKN55 DIPCS10 outlier exactly). Fixed: `GFNFF::generateCoulombSelfEnergyNative()` (`gfnff_method.cpp`) derives the per-atom self-energy inputs straight from the topology, independent of the pair list (mirrors Fortran `gfnff_engrad.F90:1378-1389`'s unconditional per-atom self term); `FFWorkspace::setInteractionLists()` prefers it, falling back to the old pair-derived vectors only for UFF/QMDFF (no EEQ). (b) After (a), Mg2+/Na+ matched xtb exactly but Li+/Be+/Be2+ still showed a ~`-0.04*q^3` Eh offset: `EEQSolver::calculateDgam()` (`eeq_solver.cpp`) had the metal-type hardness correction (`imetal==1 → ff=-0.08`) nested inside `else if (Z > 10)`, but the Fortran reference (`gfnff_ini.f90:658-660`) applies it as an unconditional independent `if` — Li(Z=3)/Be(Z=4), the only metals with Z<=10, never reached the check while Na/Mg (Z>10) did. Fixed by moving the metal + noble-gas correction out of the `Z>10` branch. Result: 107/107 GMTKN55 single-atom structures now match xtb to <0.01 kcal/mol (was 0/107); GMTKN55 gfnff MAD **19.98 → 3.389** kcal/mol (RMSD 131.5→17.45, max 2432.6→500.0). Full `ctest` re-run clean (58/58 gfnff-labelled pass; the 21 full-suite failures are pre-existing/environmental — missing `release_tblite/` dump tree, documented `cli_curcumaopt_07` golden-value drift, 3 unit binaries stale from before this session — none touch GFN-FF/EEQ). Remaining GMTKN55 gfnff residual (MAD 3.389, unrelated to either fix) clusters into charged-HB (AHB21), SN2/proton-transfer TS (BH76) and non-classical-bonding (MB16-43/AL2X6/YBDE18/DC13) categories, same bug class as Known Issues #6, not individually root-caused. See [docs/GMTKN55_VALIDATION.md](docs/GMTKN55_VALIDATION.md).

9. **`-sp` has no working open-shell (UHF) path for native gfn1/gfn2 (Sep 2026, AI/machine-tested, found via GMTKN55)**: `-spin N` only sets `Molecule::m_spin`, which is read back in `src/main.cpp`/`src/capabilities/curcumaopt.cpp` for display/`-opt` but nowhere under `src/core/energy_calculators/` — the actual occupation channel `EnergyCalculator::m_mult` reads top-level `controller["multi"]` (`energycalculator.cpp:143`), and `-multi N` on the `-sp` command line is not force-routed to that top-level key the way `-charge`/`-spin` are (`main.cpp:693`), so it never reaches the SCF. Verified on GMTKN55 RC21/me (CH3 radical, GFN2): `-spin 0` and `-spin 1` give the *identical* energy (−3.56397536 Eh, silently closed-shell), vs xtb's true UHF doublet −3.56265832 Eh (0.83 kcal/mol off here, likely much larger for other radicals/ions). `scripts/gmtkn55_compare.py` works around this by skipping every structure with a nonzero `.UHF` for gfn1/gfn2 rather than reporting a wrong closed-shell number (gfnff unaffected — no explicit open-shell term). Affects roughly a third of GMTKN55 (RC21, RSE43, G21EA/G21IP, parts of BH76, …) and any other open-shell `-sp`/`-opt` use of native gfn1/gfn2. Not yet fixed.

10. **GFN-FF: missing "N/S is not a pi-atom" veto let a metal/ion-coordinated pyrrole N corrupt its ring's Hückel treatment — FIXED (Sep 2026, AI/machine-tested)**: found via the new S30L-CI dataset. `test_cases/s30lci_test_set/29` (a calix[4]pyrrole-type host binding an anion with an explicit Na+ counterion — S30L-CI adds explicit counterions to the classic S30L's charged complexes 23-30) showed curcuma GFN-FF vs xtb off by −55 kcal/mol for the AB complex, with A and B alone matching to <1e-6 Eh. Root-caused against the actual Fortran port source (`pprcht/gfnff`, built at `/home/conrad/src/curcuma/external/gfnff/_build/gfnff` — not present in this checkout's `external/`; add `pr=.true.` for the per-bond/Hückel debug tables used below), which is a **closer, more reliable ground truth than xtb 6.7.1** here (pprcht DE=+85.5 kcal vs xtb +95.5, only 10 kcal apart; curcuma was +40.2 kcal, 45 kcal off pprcht) — same "curcuma tracks pprcht, xtb has since diverged" pattern as Known Issue #6. Two false leads eliminated with hard evidence: (a) the 1-vs-2-fragment disagreement with xtb's diagnostic print is xtb-only noise (pprcht **also** reports 1 fragment/122 atoms); (b) an extra Na···C ionic contact (2.88 Å) curcuma perceived is already correctly dropped by the documented q-loop second pass (`calculateTopologyInfo()`, `gfnff_method.cpp:8606`) — the resulting 3 Na-metal bonds match pprcht to <1 kcal/mol. **Actual cause**: diffing all 129 per-bond energies against pprcht's printed bond table isolated the error to one 5-membered N-heteroaromatic ring (atoms 7-8-25-24-15). pprcht's own instrumented Hückel dump showed its pi-system for this ring has **4 atoms/4 electrons — the pyrrole N is excluded** (its 5th bond, to Na, means its lone pair is no longer free to donate); curcuma built a **5-atom/6-electron** system, wrongly including the N. Traced to `GFNFF::detectPiSystems()` (`gfnff_method.cpp:6305`): it was missing two unconditional pi-atom vetoes present in the Fortran source (`gfnff_ini.f90:918-919`, `if (at(i)==7 .and. topo%nb(20,i)>3) cycle` / `if (at(i)==16 .and. hyb(i)==5) cycle` — "NR3-X is not a pi" / "SO3 is not a pi"). Ported both verbatim (using the full, non-eta-filtered neighbour count `nb_full[i].size()`, matching Fortran's `topo%nb(20,i)`, added as a new parameter to `detectPiSystems()`). **Result**: system 29 AB −22.40293→**−22.33042 Eh** (pprcht −22.330686, now 0.16 kcal/mol off, was 45); S30L-CI set-wide curcuma-vs-xtb MAD **1.920→0.403** kcal/mol, system 29's own deviation −55.27→**−9.77** kcal/mol (residual matches the known pprcht-vs-xtb metal-bond-r0 fine-precision gap, not a curcuma bug). **Regression-checked against all 3 reference sets**: MOR41 gfnff — all 95 structures bit-identical before/after (no N/S atom in this set trips the new veto); GMTKN55 gfnff — MAD 3.389→3.462 kcal/mol (noise-level shift, same pre-existing AHB21/BH76/MB16-43/DC13 outlier clusters as Known Issue #8, no new outlier category); full `ctest` — same 3 pre-existing unrelated failures as before the fix (`test_orca_interface`, `xtb_cpscf`, `cli_curcumaopt_07_opt_multixyz` golden-value drift), confirmed bit-identical to the pre-fix binary. Since Na is `metal_type>0` but not a 3d/4d/5d transition metal, and MOR41/GMTKN55 don't exercise a plain closed-shell N-heterocycle coordinated by a light/ionic metal, this path was apparently untested before S30L-CI. Two small env-gated debug additions from this session, following the existing `CURCUMA_NBDIAG`/`CURCUMA_HYBDIFF` convention (zero-cost when unset): `CURCUMA_BONDDUMP=1` (`gfnff_method.cpp`) prints every perceived bond and per-bond dynamic r0/fc/alpha/fqq/CN; `CURCUMA_HUCKELDUMP=1` (`huckel_solver.cpp`) prints each pi-system's per-atom (Z, hyb, tag, electron count).

11. **GFN-FF: amide-hydrogen chi correction silently skipped whenever a molecule's pi atoms aren't a contiguous block starting at atom 1 — FIXED (Sep 2026, AI/machine-tested)**: found via a term-by-term audit of the S30L-CI set's smaller (<1 kcal/mol) curcuma-vs-xtb residuals, requested after Known Issue #10's fix. `test_cases/s30lci_test_set/17`'s guest B (a small symmetric cyclic bis-amide, 14 atoms) was off by −5.26 kcal/mol vs pprcht, isolated to the Coulomb/EEQ term alone (every other term — bond, angle, torsion, dispersion, HB, BATM — agreed to <0.002 kcal/mol). Root-caused by instrumenting the Fortran reference directly (a temporary `write(*,*)` in `goed_gfnff`'s caller, `gfnff_engrad.F90:307`, printing its converged charges — reverted after use): pprcht's amide N/H charges (N=−0.2796, H=+0.1900) were markedly less polarised than curcuma's (N=−0.3161, H=+0.2290). Traced to the Fortran source having **two different pi-membership arrays** feeding the *same* `amide()` function for two *different* corrections: `gfnff_ini.f90:651` (the dgam `ff=-0.16` branch) passes `piadr` — the deliberately-preserved buggy index-cutoff array (`piadr(i)!=0` behaves as `i<=npiall`, not "is atom i a pi atom", because of how the array is populated; GFN-FF's parameters are fit against this quirk, per Known Issue #10's sibling finding) — while `gfnff_ini.f90:673` (`chieeq(H) -= 0.02` via `amideH()`) passes `piadr2`, the *correct* atom-indexed pi-membership array. `EEQSolver::detectAmideNitrogens()` (`eeq_solver.cpp`) only ever computed the buggy variant, and `calculateFinalCharges()` fed that same result into **both** `calculateDgam()` (correct) **and** `detectAmideHydrogens()` (should have been the correct-membership variant). Whenever a molecule's pi atoms are not exactly original-atom-numbers 1..npiall — true for any simple cyclic amide with an sp3 CH2 interleaved between the pi atoms, as here — the buggy check wrongly returns "not amide" for real amide nitrogens, dropping the −0.02 chi correction on their H's and overpolarising the whole N-H-C(=O) unit. Fixed by adding an `exact_pi_membership` parameter to `detectAmideNitrogens()` (false = preserve the existing piadr-bug replication for the dgam branch, true = real `is_pi_atom[]` membership) and computing a second, correct-membership `is_amide` specifically for the `detectAmideHydrogens()` call. **Result**: system 17 B's charges now match pprcht's to 6 decimal places exactly; Coulomb energy −0.15515→**−0.14755 Eh** (pprcht −0.147550, from 5.26 kcal off to exact); S30L-CI MAD 0.403→**0.392** kcal/mol, system 17's own curcuma-vs-xtb deviation −0.29→**0.00**. Caveat hit during testing, not the bug itself: a stale per-structure `.topo.json` cache (format unchanged, so silently reused) masked the fix on the first rebuild — always delete `*.topo.json` under a test directory before re-measuring after a GFN-FF topology/charge change. **Regression-checked**: MOR41 gfnff all 95 structures bit-identical before/after; GMTKN55 gfnff MAD unchanged at 3.462 (identical to three decimals); full `ctest` — same result set as Known Issue #10 (`confscan_dtemplate` is flaky/pre-existing, reproduced identically with and without this fix; the other 3 failures are the same pre-existing ones). **Correction (Known Issue #12)**: the "GMTKN55 MAD unchanged at 3.462" claim here and in #10 was measured through `gmtkn55_compare.py`'s energy cache and is therefore not evidence either way.

12. **GFN-FF: four independent term-level bugs found by a systematic per-term diff against pprcht — FIXED (Sep 2026, AI/machine-tested)**: instead of chasing single outliers, every S30L-CI fragment (30 x A/B/AB = 90) was diffed against pprcht **term by term** (bond/angle/torsion/repulsion/Coulomb/dispersion/HB/XB/BATM, curcuma `-verbosity 2` vs the reference's own decomposition). 75/90 matched to <0.001 kcal; the other 15 fell into exactly three clusters, each one bug. A fourth was then found the same way in GMTKN55.
    - **(a) Nitro nitrogen lost one pi electron — up to 6.9 kcal/mol.** `HuckelSolver::countPiElectrons()` (`huckel_solver.cpp`) had nitrogen as an `else if` chain. The Fortran (`gfnff_ini.f90:917-920`) uses three **independent, overlapping** `if`s: a nitro N (Z=7, hyb=2, itag=1) matches *both* the itag rule and `hyb<=2` and therefore contributes **2** electrons — its own comment says so ("the itag=1 avoids an odd el number for the nitro group (its 4)"). curcuma gave it 1, so every polynitroaromatic entered the Hückel solve with an odd electron count and took the open-shell path, producing a qualitatively wrong pi-bond-order pattern. Only nitrogen has overlapping conditions; the `else if` form is correct for B/C/O/F/S/Cl. Clean check: **nitrobenzene −2.82982873 → −2.83730279 Eh**, matching pprcht (−2.83730275) *and* xtb (−2.83730277). S30L-CI 3/4/5 (benzofurazan guests): +6.91/+4.95/+4.92 → 0.0000 kcal.
    - **(b) Triple-bond torsion (sTors): the reference evaluates only its LAST entry, m times.** `gfnff_engrad.F90:494-503` loops `do i=1,m` but calls `sTors_eg(m, n, ...)`, and `sTors_eg` reads `topo%sTorsl(:,m)` — the array SIZE, not the loop index. Both pprcht **and xtb 6.7.1** therefore drop every entry but the last (and give exactly 0 when the last slot was never filled). Proven on S30L-CI 7/B (alkyne-linked C40H20 macrocycle): the 5 detected quartets match the reference exactly, curcuma's sum 6.125186e-4 Eh equals the observed difference, and `5 x E(entry 5)` = 0 equals the reference's contribution; on 7/A the same model reproduces the *opposite*-signed 8 x E(entry 8). curcuma keeps the correct summation as the default — `erefhalf` is a DLPNO-CCSD(T) diphenylacetylene value, not a fitted parameter — and reproduces the reference bit-for-bit under the new opt-in `-gfnff.storsion_reference_loop_bug true`. Systems 7/8: ±0.14…0.58 → 0.0006 kcal under the flag.
    - **(c) The N-outer torsion `fkl` halving used heuristics instead of the reference's pi array — 0.17 kcal/mol.** `gfnff_ini.f90:1693-1694` halves `fkl` when an outer torsion atom is a nitrogen with `piadr==0`. Crucially, **`gfnff_ini.f90:1016` replaces `piadr` with `itmp` right after the Hückel section** (`piadr = itmp`), so from there on `piadr` is a plain atom-indexed membership flag — *not* the index-cutoff array that Known Issue #11 documents for the earlier EEQ call sites. curcuma's `getGFNFFTorsionParameters()` instead re-derived "N in pi" from `pibo > 0.1` or "has an sp/sp2 neighbour", which reports a nitrogen kept OUT of the pi system but ringed by sp2 carbons (the Na-coordinated pyrrole N of S30L-CI 29) as "in pi" — 14 torsions at exactly twice the reference force constant. Fixed by materialising the reference array: `HuckelSolver::calculatePiBondOrders()` now also fills `TopologyInfo::pi_atoms_final` (1 for both ends of every bond inside a **solved** pi-system, mirroring the Fortran `itmp` loop), and the torsion rule tests it directly. Topology cache **version bumped 2 -> 3** so a v2 cache cannot silently supply the missing array. 29/AB +0.166 → +0.047 kcal (remainder is Coulomb).
    - **(d) The X-bond B-search cutoff was the wrong pair AND 2x too short — 0.099 kcal/mol on S30L-CI, and it was the real cause of the MOR41 "irreducible fine precision" residual.** `detectHalogenBondsNative()` pruned candidate B atoms at a hardcoded 10 Bohr on the **X-B** distance; the reference prunes on the **A-B** distance against `hbthr2` (`gfnff_ini2.f90:751-757`, 450 Bohr^2 = 21.2 Bohr at its accuracy 0.1). It also rejected any non-group-4 B with `qa > 0.05`, a filter the reference applies **only** to group 4 ("must be a (pi)base", `gfnff_ini.f90:874-876`) — that dropped the mildly positive thiophene sulfurs of host 11/12 from their own X-bond list. curcuma now derives the threshold from the same `hb_accuracy`/`hb_thr2_bohr2` PARAMs as the HB list: 140 → 890 triples on 11/A, XB −0.006759 → −0.006920 Eh (reference −0.006917).
    - **(e) Hydrogen's hybridization was folded onto sp3, firing an N-sp2 rule meant for sp3 nitrogen — 21.8 kcal/mol on one bond.** `getGFNFFBondParameters()` mapped "anything outside 1..3" to 3 before the `bsmat` lookup. For the lookup itself that fold is value-preserving (row/column 0 and 3 of `bsmat` are identical), which is why it survived so long — but the special cases keyed on the folded value then saw hydrogen (Fortran `hyb=0`) as sp3. An **sp2-N–H** bond therefore matched `hybi==3 && hybj==2 && N` and got `bstrength = 1.24*1.04 = 1.2896` instead of `bsmat(2,0) = 1.0792`, ratio 1.19497 — exactly the measured force-constant ratio. The same fold made the hypervalent branch (`hyb==5`) unreachable, and curcuma's `btyp` N-sp2 rule was additionally non-directional (the reference requires the **nitrogen** to be the sp3 partner, `gfnff_ini.f90:1108-1109`). Found on GMTKN55 `Amino20x4/ARG_xak`, where curcuma was 22 kcal/mol below **both** references (this is not a reference split — pprcht and xtb agree): the whole error sat in one bond, the guanidine imine N19-H21. Fixed by testing the special cases on the raw Fortran hybridizations. ARG_xak −5.89217389 → **−5.85709879** (pprcht −5.85709874).
    **Combined validation** (all measured with `.topo.json` caches cleared and, for GMTKN55, its energy cache dropped): S30L-CI per-term vs pprcht **0/90 structures above 0.05 kcal** (was 15/90, max 6.91), MAD 0.0007, max 0.047; S30L-CI vs xtb MAD 0.392 → **0.386**; MOR41 per-structure vs pprcht MAD **0.00429 → 0.00067**, max **0.070 → 0.012**, structures >0.01 kcal 9 → 1 — this **resolves the PR30/PR31 ~0.07 kcal residual that Known Issue #6 called "residual fine-precision"** (it was the X-bond cutoff, (d)); GMTKN55 gfnff vs xtb MAD **3.462 → 2.117** (RMSD 17.73 → 16.04), `Amino20x4` MAD 5.664 → 0.001 (max 26.9 → 0.016), remaining outlier clusters unchanged (AHB21/DC13/AL2X6/MB16-43/BH76 — Known Issue #8). `ctest`: the same 4 pre-existing failures as before (`confscan_dtemplate` flaky, `test_orca_interface`, `xtb_cpscf`, `cli_curcumaopt_07_opt_multixyz` golden-value drift), each verified to fail identically with the pre-fix binary; gfnff-labelled 55/56.
    **Method note**: the per-term sweep is the reusable part — it turns "system X is off by N kcal" into "the bond term of system X is off", which in every one of these four cases pointed straight at the responsible code path. Where the term alone was not enough, the reference was instrumented with a temporary `write(*,*)` (per-bond energies in `egbond`/`egbond_hb`, torsion factors, `piadr`/`sTorsl` contents), rebuilt via `make gfnff-lib && cd app && make`, then reverted with `git checkout --`.

13. **GFN-FF: the two-fragment charge placement followed dead reference code — FIXED (Sep 2026, AI/machine-tested)**: after #12, GMTKN55 gfnff still had 184/2458 structures above 1 kcal/mol. Rather than guess, the worst ~77 were **arbitrated** by running pprcht as a third opinion: 66 turned out to be genuine curcuma port errors (curcuma != pprcht, pprcht == xtb), 5 reference splits, 6 mixed — the opposite of MOR41/S30L, where the split dominates. Per-term fingerprints then separated them into a Coulomb family (all charged) and a bond/topology family.
    - **The Coulomb family, root cause.** For a charged molecule that falls into exactly two fragments, curcuma tried both placements of the net charge and kept the lower EEQ electrostatic energy — mirroring `gfnff_ini.f90:536-560`. But that reference block is **dead code**: it is gated on `sum(topo%qfrag(2:nfrag)) > 999` while `qfrag` is pre-initialised to `[ichrg, 0, ...]`, so the sum is 0 and it never runs, in pprcht *or* xtb (same finding as Known Issue #2's note about xtb's `.CHRG` path, but this is the no-charges-file branch). The effective reference rule is simply **"whole charge on fragment 0"**. Evidence the trial is not merely different but wrong: GMTKN55 `AHB21/21` is HCOO⁻···HF with the proton unambiguously on F (1.0 Å), so the −1 belongs on the formate; the trial put it on the two-atom HF fragment, giving that hydrogen a charge of **−0.52** and the Coulomb term −0.971 instead of −1.349 Eh — **237 kcal/mol on a 6-atom structure**. Fixed by making the reference rule the default; the trial stays reachable via `-gfnff.frag_charge_autodetect true`. `AHB21/21` −1.646456 → **−2.027889 Eh**, exactly pprcht. Set-wide: `AHB21` MAD 22.05→**0.606** (max 239→4.8), `CHB6` 7.73→**1.51** (max 112→14.5), `BH76` 6.90→**4.79**.
    - **A self-inflicted regression, caught by the same sweep.** The `nh_linear_fix` guard ported from `confsearch` two commits earlier assumed every genuine sp N-H is already caught by the structural rules preceding the GEODEP angle fallback. A full GMTKN55 scan (every structure, guard on vs off) showed it changes exactly **2 of 2462** — and both were genuine sp centres it wrongly demoted: `DIPCS10/n2h2_2+` (linear HN=NH²⁺, **+165.5 kcal**) and `NBPRC/nh-bh` (linear HN=BH, **+78.4**). Neither is reached by the structural rules (their partner is N resp. B, not a 1-coordinate C/N). Discriminator: a genuinely sp nitrogen sits in a **linear chain**, so its heavy partner is itself 2-coordinate, whereas the artefact's =N-H hangs off a 3-coordinate sp2 carbon — the guard now additionally requires the heavy partner to have ≥3 neighbours. `DIPCS10` MAD 8.28→**0.000**, `NBPRC` 3.84→**0.108**, artefact still guarded, MOR41+S30L-CI bit-identical. **This makes the guard differ from `confsearch`'s version; the ORIGIN comments in `gfnff_method.cpp`/`gfnff.h`/`docs/GFNFF_STATUS.md` say so and mark this branch as the newer one.** Worth backporting.
    - **Combined**: GMTKN55 gfnff MAD **2.117 → 1.445** (RMSD 16.04 → 12.83; >1 kcal 184→172, >20 kcal 69→53). MOR41 and S30L-CI bit-identical (the charged-two-fragment path is not exercised there — S30L-CI is neutral by construction). `ctest` unchanged.
    - **Still open** (all measured, none root-caused): a bond/topology-perception family — `DC13/c20bowl` (+564 kcal in the bond term alone), `AL2X6` bridged dimers, `ALK8` Li clusters, `HEAVY28`/`HEAVYSB11` heavy main-group hydrides — plus strained/hypervalent cases (oxiranes +24…+47, H2S2O7 +12, N-ylides −36) and `MB16-43`, the one cluster where pprcht and xtb genuinely disagree with each other. Separately, Known Issue #9 keeps 320 GMTKN55 structures out of the gfn1/gfn2 comparison entirely.
    - **Method note**: `scripts/gmtkn55_compare.py` caches energies in `_run/energies.json`; several "MAD unchanged" observations earlier in this work were cache artefacts. Drop the `<subset>/<name>|cur|<method>` keys (or pass `--recompute`) before believing any before/after number.

15. **GFN-FF: three metal/perception rules — main-group metals miscounted, TM-TM over-assigned, GEODEP sp2->sp3 missing — FIXED (Sep 2026, AI/machine-tested)**: the AL2X6 subset (bridged Al2X6 dimers, MAD 14.68 after #14) turned out to be three unrelated bugs, all found by the same per-term-then-per-bond drill.
    - **(a) The metal-neighbour count `nm` only knew transition metals.** `EEQSolver::calculateDxi()` counted metal neighbours with a hardcoded 3d/4d/5d range, but the reference tests `imetal(j) /= 0` (`gfnff_ini.f90:372`), and `imetal` is `param%metal(Z)` — which **includes the main-group metals** — demoted to 0 only for a low-coordinate element of group > 3 (`:273-274`). Aluminium was therefore invisible, so the polyvalent-halogen dxi rule took its `nm == 0` branch: the bridging chlorines of `al2cl6` got `-nn*0.021` instead of `+nn*0.05`, a 0.142 shift in `chieeq` that **inverted their topology charge** (−0.261 vs the reference's +0.158) and, through `fqq`, cost 52 kcal/mol in the bond term. `al2cl6` −60.07 → **+0.00 kcal**.
    - **(b) Any two bonded metals were classified TM-TM.** `classifyBondType()` carried an explicit TODO ("Simplified: If both are metals, assume TM-TM"); the reference requires `imetal==2` on both (`gfnff_ini.f90:1121`). The Al-Al bond of `al2me6` was promoted to btyp=7, which swaps `bstren(5)=1.00` for `bstren(7)=3.40`. (This classifier feeds torsion filtering; `getGFNFFBondParameters` already had the imetal==2 test, so the energy effect here was indirect.)
    - **(c) The GEODEP sp2->sp3 promotion was missing entirely.** `gfnff_ini.f90:1024-1032` re-classifies a group-4 atom with three neighbours that came out sp2 but sits outside every pi-system as **sp3** when it is markedly pyramidal (out-of-plane angle > 40 deg). It runs right after the Hückel section because it tests the post-Hückel `piadr` (our `pi_atoms_final`). Without it the bridging methyl carbons of `al2me6` stayed sp2, so their C-H bonds lost the X-sp3 r0 shift of −0.022 and took `bsmat(2,0)=1.0792` instead of `bsmat(3,0)=1.0000` — **8.6 kcal/mol on each of six bonds**. `al2me6` −46.77 → −3.67 kcal; its bond and angle terms now match the reference exactly.
    - Also corrected while in the area: the high-coordinate-metal torsion skip tested curcuma's broad `is_metal` flag where the reference tests `param%metal(Z) > 1` (transition metals only). No measured effect on any reference set — a faithfulness fix, not a bug fix. **Honest note**: it was briefly credited with the `al2me6` torsion residual on the strength of a diagnostic that had run against a stale binary (the build was failing at the time and `make` exited non-zero unnoticed). Always check the build's exit status, not just a grep for "error:".
    - **Result**: GMTKN55 gfnff MAD **1.244 → 0.966** — under 1 kcal/mol for the first time — RMSD 7.92 → 7.01, structures above 20 kcal 52 → 37. `AL2X6` MAD 14.68 → 3.34, `MB16-43` 13.66 → 10.77, `ALK8` drops out of the worst list. MOR41 per-structure vs pprcht unchanged (MAD 0.00067); S30L-CI bit-identical; `ctest` unchanged.
    - **Still open on al2me6** (−3.67 kcal): its torsion term is exactly 0.0 where the reference has 0.005849 Eh. 110 torsion quartets are generated (the same count as the reference) but every barrier comes out below the storage threshold. Not root-caused; bond, angle and all non-bonded terms of that molecule are exact.

16. **GFN-FF: two angle-rule ordering/duplication bugs and a missing imetal demotion — FIXED (Sep 2026, AI/machine-tested)**: taking the *unambiguous* GMTKN55 outliers first (curcuma != pprcht while pprcht == xtb), two families each showed a **constant** error per structural motif, which is the signature of a single mis-set parameter rather than an accumulation.
    - **(a) Oxiranes: the ring angle was overwritten by the generic ether rule — +25 kcal per 3-ring.** `oxirane`, `propyloxirane`, `dimethyloxirane` and `ISO34/P25` were all +24.3…+25.5 kcal off, `dioxirane` (two such centres) +47.4 — purely in the **angle** term. curcuma's oxygen block ("O with two neighbours -> theta0 = 104.5 deg") sat AFTER the ring-strain block, so a 3-ring ether lost its `theta0 = 82` again. The reference has the oxygen rules at `gfnff_ini.f90:1582-1598` and the ring rules at `:1635-1649` — the ring wins. Moved the oxygen block ahead of the ring block (it only reads `nh/nsi/nmet/npi` and the current angle, all computed further up). **All five structures now exact (0.00 kcal).**
    - **(b) PbH4: a duplicate metal block and a raw imetal — +27.2 kcal per PbH4 unit.** `pbh4`, `pbh4_hi`, `pbh4_hcl`, `pbh4_hbr`, `pbh4_teh2`, `pbh4_h2o` were each +27.2 kcal off and `pbh4_2` **exactly twice** that; `pbme3` −40.0, `pbh4_bih3` +36.0. Two independent causes, both the same missing rule — `gfnff_ini.f90:273-274` demotes a low-coordinate element of group > 3 back to a non-metal ("Sn, Pb, Bi, with small CN are better described as non-metals"):
      - `getGFNFFBondParameters()` read `metal_type[Z]` raw, so PbH4 took the METAL bond branch: `fqq` 1.1437 instead of 1.033 and the force constant −0.0523 instead of −0.064 (the bond half, +29.3 kcal).
      - The angle code carried **two** metal-centre blocks. The second is the documented one and applies the demotion; the first was an undocumented duplicate that did not, and it re-imposed `theta0 = 109.5` after the heavy-main-group rules had correctly produced `109.5 - nh*5 = 99.5`. For a perfectly tetrahedral PbH4 that makes theta0 equal the actual angle and the whole angle term collapses to **exactly zero** (vs 0.003308 Eh). Deleted the duplicate.
      **pbh4 +27.23 → −0.39, pbh4_2 +54.42 → −0.80, pbme3 −40.03 → −0.07, pbh4_bih3 +36.02 → −1.15.**
    - **Result**: GMTKN55 gfnff MAD **0.966 → 0.715**, RMSD 7.01 → 6.46, structures above 20 kcal 37 → 22 and above 5 kcal 81 → 53. `HEAVY28` MAD 8.68 → **0.407** (max 54.4 → 2.1), `HEAVYSB11` 4.27 → **0.031** (max 37.8 → 0.2), `ISO34` → 0.012, `FH51` → 0.060. MOR41 per-structure vs pprcht unchanged (MAD 0.00067, the demotion only touches `imetal==1` so transition metals are unaffected); S30L-CI bit-identical; `ctest` unchanged.
    - **Method note**: a *constant* offset repeated across structures that share one motif, and doubling when the motif appears twice, points at a single wrong parameter for that motif — worth checking before any per-structure analysis. Both families here were found that way from the arbitration table alone.

17. **GFN-FF: the q-loop's second pass re-detected fragments, silently dropping the EEQ constraints of every charged complex — FIXED (Sep 2026, AI/machine-tested)**: the anionic/cationic cluster (AHB21, BH76 SN2 transition states, G21EA, SIE4x4, CHB6 — 50-150 kcal/mol each) was one bug. All of them were **Coulomb-only** in the term fingerprint, with identical fragment counts and charge placement reported by both codes — yet the energy differed by 100+ kcal. Tracing the actual numbers the energy path uses (not the ones the EEQ diagnostic prints) showed the workspace running with `nfrag=1`, charges (−0.5, −0.5) instead of (−1, 0), and every charge-dependent quantity following (`alpeeq` computed from qa = −0.5, `dgam` likewise, hence the Coulomb self-energy). **Cause**: GFN-FF derives the topology twice — pass 1 with qa = 0, pass 2 with the pass-1 charges shrinking the bond radii. For a charged species the anionic atom's radius *grows*, and a contact can cross the bond threshold: GMTKN55 `G21EA/EA_25` is the dichlorine radical anion Cl₂⁻, whose 2.73 Å contact becomes a bond in pass 2 — **in both codes** (`#bonds : 1`). But the reference gates its entire fragment block on `if (topo%nfrag <= 1)` (`gfnff_ini.f90:467`), so pass 2 **keeps** pass 1's fragmentation and `qfrag`; curcuma re-detected them, so the two fragments merged into one and the per-fragment EEQ constraints vanished. Fixed by carrying pass 1's `nfrag`/`fraglist`/`qfrag` into pass 2 (three `m_frag_carry_*` members, cleared around the call).
    - **Result**: GMTKN55 gfnff MAD **0.715 → 0.460**, RMSD 6.46 → **3.84**, set maximum **151.5 → 105.5**, structures above 20 kcal 22 → 14. `BH76` MAD 4.79 → **0.080** (max 151.5 → 2.1), `SIE4x4` 2.60 → **0.000**, `CHB6` 1.51 → **0.000**, `G21EA` 2.65 → 0.653. `EA_25` −99.78 → **exact**, `hoch3fts` −150.09 → **exact**, `h2o2+_1.0` −59.86 → **exact**, `fch3fts` −151.45 → +1.18. MOR41 per-structure vs pprcht unchanged (MAD 0.00067); S30L-CI bit-identical; `ctest` unchanged.
    - **The physics, checked against an external reference** (this is port fidelity, not accuracy — see [docs/REV_GFNFF_TODO.md](docs/REV_GFNFF_TODO.md) #3): for Cl₂⁻ at that geometry, E(Cl₂⁻)−E(Cl⁻)−E(Cl) is **−41.5 kcal/mol at r²SCAN-3c** and −33.9 at GFN2 (experiment ~−30), while GFN-FF gives **−6.6 with nfrag=2** (the reference behaviour we now reproduce) and **−106.3 with nfrag=1**. Neither is right; the fragment count is the only lever and it is a discrete switch that brackets the truth. The bond term is essentially identical in both regimes (−11.5 vs −11.2 kcal) — the entire difference is the EEQ, which values the delocalisation from q=(−1,0) to (−0.5,−0.5) at −94.8 kcal/mol where the truth is ~−35, because its self-energy scales as q². Fixing that needs a different delocalisation term, i.e. a method change. We took the reference behaviour because it is the **smaller** error and it is self-consistent.
    - `BH76/fch3fts` was the last member and is now exact — see Known Issue #19.

19. **GFN-FF: a non-reference "charges must be within ±1" guard silenced the angle `fqq` for every charged fragment — FIXED (Sep 2026, AI/machine-tested)**: `BH76/fch3fts` (the F···CH3···F⁻ SN2 transition state) was +1.18 kcal/mol off, angle term only. Its angle *set* and its H-C-H force constants matched the reference exactly; only the three **F-C-H** constants came out 0.2631 against 0.224. Cause: curcuma wrapped the angle charge factor in `if (|qa_center| < 1.0 && |qa_i| < 1.0 && |qa_k| < 1.0)`, a guard the reference (`gfnff_ini.f90:1426-1430`) does not have. It fires exactly when an atom carries a full unit charge — which is what a charged fragment does: here the leaving fluoride is its own fragment with `qa = -1.000000` exactly, so all three F-C-H bends silently kept `fqq = 1.000` instead of 0.8521. Hand-check with the reference formula, `fqq = 1 - (q_C q_F + q_C q_H)·qfacBEN` with qfacBEN = −0.54, gives 0.852051, and 0.381·0.852051·0.853·0.810 = **0.2243** = the reference's 0.224. Guard removed. While there, the same function's metal branch was switched from curcuma's broad `is_metal` flag to `imetal` (`param%metal(Z)` with the low-coordinate group>3 demotion), matching the reference and the fixes of Known Issue #15 — no measured effect on any reference set.
    - **Result**: `fch3fts` +1.18 → **exact**. GMTKN55 gfnff MAD 0.447 → **0.446**, `BH76` MAD 0.080 → **0.057** (max 2.1). MOR41 per-structure vs pprcht unchanged; S30L-CI bit-identical; `ctest` unchanged.
    - **Checked against higher-level references** ([docs/REV_GFNFF_TODO.md](docs/REV_GFNFF_TODO.md) #6): an umbrella scan of that transition state gives 6.43 kcal/mol at r²SCAN-3c and 6.03 at GFN2 for a 10° distortion, against **2.13 for GFN-FF** — the angle term there is a factor of three too soft, and this fix softens the F-C-H constants further. So port fidelity again moves away from the physics (as in Known Issue #17), and it is still the right call: the accidental stiffening applied only to angles at a full-unit-charge fragment, so it was never a principled correction.

20. **GFN-FF: the hypervalent torsion correction was never implemented, multiplying every torsion at a hypervalent centre by 5 — FIXED (Sep 2026, AI/machine-tested)**: `ICONF/H2S2O7_*` and `PArel/h2s2o7*` are the same molecule, disulfuric acid, and were +12.6 kcal/mol off in the **torsion** term. A per-torsion diff showed **all twelve** torsions at exactly 5.00x the reference force constant, and instrumenting the reference's own `fctot` assembly narrowed it to a single factor: every one of `f1`, `f2`, `fqq` and `fkl` agreed to six digits, only `fij` differed — 0.059490 against 0.011898, a ratio of exactly 5. That is `gfnff_ini.f90:1811`, `if (btyp(m) .eq. 4) fij = fij*0.2d0`, the hypervalent central-bond correction. curcuma carried it as an explicit TODO ("Bond type detection not yet fully implemented … Deferred (low impact - rare in organic molecules)") — and it could not have worked anyway, because `classifyBondType()` folded hyb 5 onto 3 until Known Issue #15 removed that fold, so `btyp == 4` was unreachable. Implemented at the reference's position, last, after the alphaCO and amide scalings. **`ICONF/H2S2O7_1` and `PArel/h2s2o71` both exact.**
    - **Result**: GMTKN55 gfnff MAD **0.446 → 0.416**, `ICONF` MAD 1.78 → **0.398** (max 12.6 → 7.8), `PArel` 0.95 → **0.348** (max 14.0 → 2.6). MOR41 per-structure vs pprcht unchanged; S30L-CI bit-identical; `ctest` unchanged.
    - **Checked against higher-level references**: ICONF is a *conformer* set, so its published relative energies are the right yardstick. For the three H2S2O7 conformers, ΔE(2−1)/ΔE(3−1) in kcal/mol: published **0.55 / 3.55**, r²SCAN-3c 0.28 / 2.39, GFN2 1.23 / 6.53, GFN-FF after the fix 1.47 / **1.11**, before 1.41 / 0.65. The fix moves ΔE(3−1) towards the reference and ΔE(2−1) marginally away; GFN-FF is ~2 kcal wrong either way. So the correction matters for **absolute** energies (port fidelity) while conformer differences stay dominated by other errors — worth knowing before reading too much into a conformer-set MAD.
    - **Still open**: `BHDIV10/ts2` (+12.73, bond +9.49 and torsion +3.49) and `ISOL24/i13p` (−16.69, almost entirely bond) — both curcuma-vs-pprcht, xtb agrees with pprcht. `AL2X6/al2me6` (−26.3 vs xtb) is now **exact against pprcht**, i.e. a pure reference split in the repulsion term.

18. **GFN-FF: the q-loop's second pass never saw the charges it exists for — FIXED (Sep 2026, AI/machine-tested)**: `G21EA/EA_9` is the methylene anion CH2- (H-C-H = 99.9°), and curcuma had it 24.7 kcal/mol too low, split between the angle and Coulomb terms. Three defects in a row, all around the same rule — `gfnff_ini2.f90:250`, `if (topo%qa(i) < -0.4) then hyb=2; itag=0`, which un-tags a two-coordinate group-4 atom that the angle criterion had called a carbene:
    - **The q-loop gate was too narrow.** Curcuma skipped pass 2 whenever the charge-shrunk radii left the bond list unchanged, on the stated premise that "the charges only ever reach the rest of the model THROUGH the bond list". They do not: `topo%qa` enters the perception at exactly two places in `gfnff_ini2.f90` — the radius shrink at :122 and this carbene override at :250. Pass 1 runs with qa = 0, so the override can only ever fire in pass 2. Gate extended with the second channel.
    - **Pass 2's hybridization could not see pass 1's charges.** In the reference `topo%qa` is a member that survives the q-loop; curcuma rebuilds `topo_info` per pass, and its `topology_charges` are still empty when `determineHybridizationFortran` runs (the Phase-1 EEQ comes later in the same function). `m_bond_qa` already carries exactly those pass-1 charges for the radius shrink — the override now falls back to them.
    - **`calculateDxi()` built its own `TopologyInput` without `itag`.** The Sep 2026 itag plumbing (Known Issue #14) reached `eeq_topology_input`, `m_eeq_topo_cache` and the ipis re-solve, but not this fourth one, so the dxi carbene rule still re-derived the tag from the bond angle and kept `dxi = -0.15` after the override had cleared it — worth exactly `q(C) * 0.15 = 0.0427 Eh` in the Coulomb self-energy.
    **EA_9 −24.70 → exact.** GMTKN55 gfnff MAD **0.460 → 0.447**, `G21EA` MAD 0.653 → **0.131** (max 24.7 → 3.3). MOR41 per-structure vs pprcht unchanged; S30L-CI bit-identical; `ctest` unchanged.
    - **Checked against higher-level references** (see [docs/REV_GFNFF_TODO.md](docs/REV_GFNFF_TODO.md) #3b): a bending scan of CH2- puts the GFN2 and r²SCAN-3c minima together at **100°** (experiment ~102°), the reference-faithful GFN-FF at 130° and the carbene variant at 145°. So here port fidelity and physics point the same way — the fix also moves GFN-FF towards the truth — but a 30° residual error remains, which no parameter in GFN-FF can absorb.

14. **GFN-FF: the missing aryne rule made every carbon-cage rim atom a carbene — FIXED (Sep 2026, AI/machine-tested)**: `DC13/c20bowl` was GMTKN55's single worst structure (+500 kcal/mol, the global max since the set was first run) and its term fingerprint put +564 kcal in the **bond** term alone. Cause: the Hückel π-electron count came out **10 instead of 20** for the 20-carbon bowl, because every 2-coordinate rim carbon carried `itag=1` ("carbene", contributes no π electron). The reference tags them too — but then **clears the tags again**: `gfnff_ini2.f90:341-351` runs after `topo%nb = nbdum` and, for two BONDED carbons that are both tagged, resets both ("the very special situation of two carbene C bonded which is an arine"). Along a cage rim every tagged carbon has a tagged neighbour, so all 10 tags disappear. Ported with the reference's in-place semantics (both ends cleared as the scan proceeds), so chains of three or more adjacent carbene carbons resolve identically. Bond term +564 → **exact**; π electrons 10 → 20 = reference.
    - **Second half of the same bug.** With the tags now correct the structure was still −46.6 kcal off, all of it Coulomb (−0.1017 vs the reference's −0.000295 Eh — essentially zero, as it must be for a neutral homoatomic cage). `EEQSolver::calculateDxi()` did not read `itag` at all: it re-derived "is this a carbene" from the bond angle, which reproduces only the first half of `gfnff_ini2.f90:244-253` and misses both later corrections the reference applies to that array — the `qa < -0.4` override and the aryne rule above. Every rim carbon therefore kept a spurious `dxi = -0.15`, driving the EEQ charges ~18x too large. `TopologyInput` now carries the real `itag` and `calculateDxi()` tests it (falling back to the old heuristic when the caller supplies none). **c20bowl total: +500.0 → −6.29 kcal/mol.**
    - **Result**: GMTKN55 gfnff MAD **1.445 → 1.244**, RMSD 12.83 → **7.92**, and the set-wide maximum **500.0 → 151.5** — the first time the global max has moved since the set was introduced. `DC13` max 500.0 → 6.29 and it drops out of the worst-subset list. MOR41 and S30L-CI bit-identical; `ctest` unchanged.
    - **Residual on c20bowl** (−6.29 kcal, Coulomb only): charges still ~6x the reference's tiny ones. dxi is now 0 for every atom (pure carbon reaches no other dxi rule), so the remainder is in the CN or the Phase-1/Phase-2 coupling, not in the tagging. Not root-caused.


21. **GFN-FF: twelve porting defects found by arbitrating every remaining GMTKN55 outlier against pprcht — FIXED (Sep 2026, AI/machine-tested)**: the method from Known Issue #13 was applied to the whole residual list: run the port source **pprcht/gfnff** as a third opinion on each outlier, so that curcuma != pprcht == xtb marks a curcuma bug and curcuma == pprcht != xtb the known reference split (#6); then diff term by term, and where that is not enough, instrument the Fortran (`REFBOND` in `egbond`, `REFTORS`/`REFANGL` in `gfnff_ini.f90`, `REFHYB` in `gfnff_ini2.f90`) and compare the factors one at a time. Twelve defects came out, in five different parts of the code.
    - **(a) Boron was in the `picon` set — +12.7 kcal/mol.** `GFNFF::detectPiSystems()` admitted a non-sp/sp2 boron as a lone-pair donor, but the reference's gate is `nofs()` = {N,O,F,S,Cl} (`gfnff_ini2.f90:1397-1402`); boron is in `pilist` and therefore enters only through the `piat` branch, i.e. when it is itself sp or sp2. It has no lone pair to donate. In `BHDIV10/ts2` (a bicyclic B/N transition state) the sp3 boron linked two isolated C=C units into ONE 5-atom/4-electron Hückel system where the reference solves TWO 2-atom/2-electron ones. `ts2` +12.73 → **0.000**.
    - **(b) `detectPiSystems()` rewritten as a verbatim per-atom port — up to 31 kcal/mol.** curcuma classified pi membership per BOND; the reference decides per ATOM (`gfnff_ini.f90:312-336`): `piat = (hyb==1|2) and pilist(Z)`, `picon = (any neighbour is sp or sp2) and nofs(Z)`, then links every bonded pair of members with no bond-type test at all, and lets a member without any pi partner form its own one-atom system (skipped later by `npi < 2`). Two consequences were missing: the sp/sp2 neighbour that triggers `picon` may be of ANY element — **hydrogen included** — and membership does not require a partner. In `PX13/hf_2_ts` (the HF-dimer proton-transfer TS) both fluorines are picon precisely because the bridging hydrogens come out sp, so the reference builds a 2-atom F···F pi-system with piBO≈1.5e-4 and curcuma built none. `ISOL24/i13p` −16.69 → 0.000, `i6e` −11.31 → 0.000, `W4-11/s4-c2v` +31.49 → 0.000.
    - **(c) The topology cache did not invalidate the EEQ Cholesky cache — 22.5 kcal/mol, and the energy depended on whether a file existed.** `solveWithSchurCholesky` keys its factorization cache on geometry + CN only (`eeq_solver.cpp:2130-2146`), while the matrix also depends on dgam/alpeeq/dxi — the quantities Phase 1 produces. The design is sound only because Phase 1's local-solve branch sets `valid=false` before every Phase 2. A `.topo.json` cache hit SKIPS Phase 1, so that invalidation never happened, and in the q-loop (identical geometry and CN in both passes) pass 2 scored a false cache hit and solved with pass 1's factorization. `PX13/hf_2_ts` came out −0.20648 with the cache file present and −0.24236 without it. Same failure class as the Jul 2026 pending-buffer fix at `eeq_solver.cpp:2750`. **Not a port bug — a reproducibility bug**, and it silently contaminated earlier measurements, since 2462 `*.topo.json` files sit next to the GMTKN55 structures.
    - **(d) The pi r0 shift must REPLACE the local bond shift, not add to it — 19.2 kcal/mol.** `gfnff_ini.f90:1174` is `shift = gen%hueckelp*(gen%bzref-pibo(i))`, an assignment that discards the hypervalent/XH/F-F/X-sp3/X-sp specials set at :1143-1149; only the heavy-atom shifts at :1268 are added afterwards. curcuma computed the pi shift in a later step and ADDED it. The F-F contact of `hf_2_ts` got both the F-F rule (+0.22) and the pi shift: r0 3.0555 instead of 2.8355 Bohr. Two smaller faithfulness fixes in the same block: rules 1 and 2 there are independent `if`s that ASSIGN (curcuma had an else-if chain and a `+=`), so a hypervalent X-H bond kept `hyper_shift` where the reference overwrites it with `rabshifth`.
    - **(e) The `fxh` aldehyde rule must override the 3-ring rule — 21.3 kcal/mol.** Inside the C-H branch (`gfnff_ini.f90:1150-1157`) the two rules are independent `if`s and the aldehyde one comes second. curcuma had them as an else-if chain, so `W4-11/oxirene`'s ring carbons — a 3-ring AND a `ctype==1` carbon, since the ring oxygen is a pi atom — took fxh 1.05 instead of 0.95. That single factor is the ratio 1.10526 between curcuma's C-H force constant (−0.18250) and the reference's (−0.16512).
    - **(f) An invented sp3→sp2 promotion changed the bond type of conjugated ethers.** The reference assigns `btyp` from the RAW hyb array (`:1106-1121`); curcuma first promoted an sp3 atom sitting in a pi-system to sp2, citing a line that adjusts no hybridization at all. Oxirene's ring C-O bonds went from btyp=1 to btyp=2, giving their ring torsions periodicity 2 instead of 1 — and in a planar molecule every n=2 torsion sits exactly at its own minimum, so the whole torsion term collapsed to 0.0 against the reference's 0.001608 Eh. For nitrogen the promotion was invisible (an sp3 N on an sp2 partner reaches btyp=2 through the N-sp2 rule anyway), which is why it survived. While there, the pibo>0.1 promotion was widened from btyp==3 only to the reference's `bbtyp /= 3 .and. bbtyp < 5`.
    - **(g) The torsion pi scaling `f1 *= 0.55` must be applied once and LAST — 2.4 kcal/mol.** In the reference it is the closing line of the `if (pibo(m) > 0)` block, which sits after the ring/acyclic case AND after the SP3 specials, so it multiplies the FINAL f1. curcuma applied it up in the f2 section, where every later f1 assignment simply overwrote it — the ring cases (fr3/fr4/fr5/fr6, the terminal-atom 0.30) and the acyclic pi-sp3 override (0.5/0.2) all lost it, and two spot-fixes had been added to re-apply it for the CB7 and SP3-special branches only. `hf_2_ts`'s two 3-ring torsions had barrier 0.31958 instead of 0.17577, exactly the factor 1/0.55. Both spot-fixes removed.
    - **(h) `piadr` inside the torsion and angle setup is the POST-Hückel array — 0.56 kcal/mol, and it is what made (b) look like a regression.** `gfnff_ini.f90:1016` replaces `piadr` with `itmp` right after the Hückel section, and `itmp` marks only both ends of a bond inside a SOLVED pi-system. Everything after that line — the torsion pi-sp3 rule, `notpicon`, the f2 heavy-outer boost, the `amide`/`alphaCO`/`ctype` calls, and the angle `npi` counter — therefore reads the post-Hückel membership, which is strictly smaller than the pre-Hückel candidate list. curcuma read the candidate list (`pi_fragments`). The difference is not academic: a transition metal with three neighbours has hyb=2, so every amine/phosphine/ether donor on such a metal is a pre-Hückel candidate that the Hückel never solves. Fix (b) made the candidate list correctly larger and thereby exposed this: MOR41 `PR41` picked up 12 extra torsions at the pi-sp3 f1=0.5 instead of 1.0. Switched to `pi_atoms_final` throughout the torsion loop, the angle `npi`, and `ctype`.
    - **(i) The X-bond A-X list needs BOTH orderings of a homonuclear X-X bond — 7.0 kcal/mol.** The reference enumerates "every atom i, every neighbour ix of i that is an xatom" (`:843-846`), so a bond whose both ends are xatoms yields two ordered pairs. curcuma scanned the bond list with an if/else-if and kept only the first endpoint as the donor — and which one that was depended on the atom order in the input file. In `HAL59/BrBr_pyr` it picked the OUTER bromine, whose A-X···B angle is ~0 instead of ~180, so the X-bond term evaluated to exactly 0.0 against the reference's −0.011225 Eh. Affects every X-X bond: Br-Br, I-I, disulfide S-S, P-P.
    - **(j) The NR3 angle branch fired on an invented condition — 5.1 kcal/mol.** The reference's only test is `if (npi > 0)` (`:1518`), npi counting how many of the two OUTER atoms are pi atoms; curcuma additionally took the conjugated branch whenever an outer atom merely had hyb 1 or 2. In `WCPT18/ts4` the bridging hydrogen is two-coordinate and hence hyb=1, so the H-N-H angles of both nitrogens got theta0 = 113 deg instead of the saturated-pyramidal 104 deg (and f2 = 1−sumppi·0.7 instead of 0.40+nh·0.19). At the transition state those angles are opened to 152.8 deg, where 9 degrees of theta0 are worth 5.1 kcal/mol on a 7-atom molecule.
    - **(k) The q-loop must run twice, always.** The reference loop is `do while (qloop_count < 2 .and. gen%rqshrink > 1e-3)` with rqshrink = 0.23 fixed, i.e. exactly two passes for every molecule. curcuma gated pass 2 on "the bond list changed" plus, since Known Issue #18, the carbene override. The premise that those are the only channels was wrong a second time — see (l) — so the gate is gone. Measured afterwards: MOR41 and S30L-CI bit-identical, so the "calculateTopologyInfoOnce is not re-entrant" worry recorded at that gate does not show up on any reference structure.
    - **(l) The icase-2/3 neighbour lists are not filtered views of icase 1 — 2.7 kcal/mol.** `getnb` (`gfnff_ini2.f90:361-419`) applies the metal radius enlargement `fm` ONLY in icase 1; icase 2 (no high-coordination atoms) and icase 3 (no metals) keep fm = 1 and therefore use a STRICTER distance threshold for any pair involving a metal. curcuma derived all three lists from the single icase-1 bond list and only applied the hc_crit / metal filters, so `nb_hc` and `nb_nometal` could never lose a bond that `nbf` kept. Invisible until the charge shrink of pass 2 moves a metal pair between the two thresholds: in `AL2X6/al2f6` the Al-Al contact stays inside the enlarged icase-1 radius but falls outside the plain one, so the reference gets nbf = 5 / topo%nb = 4 and hence `nbdiff = 1`, which stops the group-3 rule (`nb20i > 4 .and. ati > 10 .and. nbdiff == 0 -> hyb = 5`) from firing. Aluminium is sp3 in the reference and hypervalent in curcuma — theta0 = 109.5 vs 90 deg and f2 = 1.0 vs 0.11 on every Al-centred angle. `al2f6` +2.85 → **0.000**. The two lists now run their own distance test.
    - **Combined result** (all measured with `*.topo.json` cleared and the GMTKN55 energy cache dropped): GMTKN55 gfnff vs xtb MAD **0.416 → 0.268**, RMSD 3.77 → **3.20**, set maximum 105.5 → **84.5**. Per subset: `BHDIV10` 1.182 → **0.043**, `PX13` 0.633 → **0.023**, `ISOL24` 0.881 → **0.300**, `W4-11` 0.716 → **0.197**, `CARBHB12` 0.445 → **0.000**, `HAL59` 0.312 → **0.000**, `WCPT18` 0.487 → 0.237, `AL2X6` 2.659 → 2.398, `MB16-43` 10.405 → 7.901, `BH76` 0.057 → 0.004, `NBPRC`/`YBDE18`/`IL16`/`RSE43`/`RC21` all lower; no subset rose by more than 0.004 (`BHROT27`, `WATER27`, noise). MOR41 per-structure vs pprcht **unchanged** at MAD 0.00068 / max 0.0123; S30L-CI vs xtb **unchanged** at MAD 0.386. `ctest`: the same four pre-existing failures (`confscan_dtemplate`, `test_orca_interface`, `xtb_cpscf`, `cli_curcumaopt_07_opt_multixyz`), and the two gfnff-labelled ones were re-run against a rebuild of the pre-session source and fail with **bit-identical numbers** there.
    - **What the residual now is.** The 14 largest non-`MB16-43` deviations were arbitrated explicitly: 13 are pprcht-vs-xtb reference splits with curcuma matching pprcht to <0.03 kcal (`AL2X6/al2me6` −26.3, `W4-11/b2h6` −16.8, `ISOL24/i11p` −8.4, `ICONF/N3P3H12_1` +7.8, `W4-11/s4-c2v` −6.3, `ISOL24/i23p` −5.8, `RSE43/P41` −4.9, the six `AHB21` 11-14, `BHPERI/13r_2` −4.1, `W4-11/hnnn` −4.1, `RSE43/P40` −4.1, both `C60ISO`), and one is a small genuine residual (`HEAVY28/pbh4_teh2` −0.39, not root-caused). So above ~0.4 kcal/mol and outside `MB16-43`, the remaining GMTKN55 gfnff deviation is the reference split, not a curcuma port error.
    - **Method note**: two of these (c and h) were only findable because a *different* fix changed the picture — (c) because a stale cache made a number move that should not have, (h) because (b) made a previously-correct-by-accident array wrong. Both are arguments for re-running all three reference sets after every single change rather than at the end.
