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
  - **Threading**: intra-molecule `-threads N` (setup 4×, gradient 3.6×) — [docs/SQM_THREADING.md](docs/SQM_THREADING.md)
  - **Eigensolvers**: opt-in `-eigensolver native|purify|lobpcg`, `CURCUMA_EIG_TRED2=blocked` (MKL-free / GPU-portable) — [docs/SQM_EIGENSOLVE_GPU.md](docs/SQM_EIGENSOLVE_GPU.md)
  - **Large systems**: `-large_system_mode fragments|dc|sparse` scales SCF past ~1000 atoms — [docs/SQM_LARGE_SYSTEMS.md](docs/SQM_LARGE_SYSTEMS.md)
  - **SCF extrapolation**: `-scf_extrapolation aspc|gauss` cuts SCF iters in opt/MD (caffeine gfn2 215→90); experimental `xlbomd` = extended-Lagrangian MD — [docs/SQM_SCF_EXTRAPOLATION.md](docs/SQM_SCF_EXTRAPOLATION.md)
  - **GPU backends** (`-gpu cuda|rocm|vulkan|auto`): all three device-resident through `-opt`/`-md`; ROCm has FP32 mixed-precision ON by default (real win), Vulkan is opt-in only (eigensolve-bound, no speedup). Detail in [docs/SQM_GPU.md](docs/SQM_GPU.md) / [docs/SQM_ROCM.md](docs/SQM_ROCM.md) / [docs/SQM_VULKAN.md](docs/SQM_VULKAN.md); roadmap in [docs/SQM_GPU_ROADMAP.md](docs/SQM_GPU_ROADMAP.md)
- ⚠️ **GFN1-xTB (Native)** - AI-implemented, machine-tested; canonical `gfn1` backend; `-opt` works; Broyden SCF default; 10/12 sqm_reference @1e-8 vs tblite (He2 + `complex` remain) — [docs/SQM_WP2_gfn1_accuracy.md](docs/SQM_WP2_gfn1_accuracy.md)
- ⚠️ **PM3/AM1/MNDO (Native NDDO)** - AI-implemented, machine-tested; 21/21 tests vs Ulysses reference (< 4 µEh)
- ⚠️ **Native GFN-FF** - AI-implemented, machine-tested; see [docs/GFNFF_STATUS.md](docs/GFNFF_STATUS.md); S30L vs xtb 6.6.1 (Jul 2026): 28/30 within ~1.5 kcal/mol after F1/F2/CLI/F3/ipis fixes (MAD 433→1.38); residuals: 23 (.CHRG quirk), 27/28+7/8 (torsion), 30/AB (bond) — see [docs/S30L_GFNNF_VALIDATION.md](docs/S30L_GFNNF_VALIDATION.md)

#### External Interfaces (Production Quality, Requires Compilation)
- **TBLite Interface** - Tight-binding DFT methods (GFN1, GFN2, iPEA1) + **Solvation** (CPCM, GB, ALPB)
- **XTB Interface** - Extended tight-binding methods (GFN-FF, GFN1, GFN2)
- **Ulysses Interface** - Semi-empirical methods (PM3, PM6, AM1, MNDO, RM1, etc.) + **Solvation** (GBSA)
- **Native GFN-FF** - Curcuma's own implementation (`gfnff`) - ✅ **IMPLEMENTED**

### 2. Force Field Methods
- **Universal Force Field (UFF)** - General-purpose molecular mechanics
- **GFN-FF** (`gfnff`) - ✅ **FULLY IMPLEMENTED** - See [docs/GFNFF_STATUS.md](docs/GFNFF_STATUS.md)
- **QMDFF** - Quantum Mechanically Derived Force Fields
- **Universal Parameter Caching** - Automatic save/load for all FF methods

### 3. Solvation Models (Implicit Solvent)
- ✅ **TBLite Solvation** - CPCM, GB (Generalized Born), ALPB for GFN methods
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
- **ConfSearch** - Systematic conformational searching (unified trajectory framework)
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
- **Multi-threading support** via `ForceFieldThread`
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
│   │   │       ├── forcefieldthread.cpp   # Multi-threading support
│   │   │       ├── gfnff_method.cpp      # Native GFN-FF implementation (4329 lines)
│   │   │       ├── gfnff.h               # GFN-FF class interface
│   │   │       ├── gfnff_advanced.cpp     # Advanced GFN-FF parameters
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

2. **GFN-FF S30L validation (Jul 2026, AI/machine-tested)**: see [docs/S30L_GFNNF_VALIDATION.md](docs/S30L_GFNNF_VALIDATION.md) — vs xtb 6.6.1: 28/30 within ~1.5 kcal/mol (MAD 433→1.38). Fixed: (a) charged complexes (23-30) EEQ qfrag=[0,0] → xtb-style both-assignment trial (F2); (b) `-charge -N` parsed as +N → negative-number CLI fix; (c) `AngleBending` 1/sinθ NaN guard (F1); (d) F2 cache bug — inline topo-cache load dropped qfrag (cached charged nfrag==2 re-runs neutralised Coulomb; now restores qfrag, write guard checks SUM); (e) F3 bond/Hückel/Coulomb — S CN=2 sp2→sp3 (11/12), halogen hyb + hoffdiag default 0 + π-system membership aligned to xtb (F/Cl in, Br/I out), `metal_type[86]` array had 83 entries (Kr/I/Xe missing → I read as TM → dgam ff=-0.9 not -0.07 → 15/16 Coulomb -14 kcal); array rewritten to xtb metal(86), 15/16 fixed; (f) F3d Hückel ipis — xtb subtracts the π-system charge from nelpi (`nelpi -= ipis`); curcuma didn't → charged hosts had nelpi too large → wrong pibo → bond off; now computes ipis (qheavy + neutralize-fragment + re-EEQ + dqa·1.1) and subtracts; fixed 25/26 (→0.0) and 30/B. Residuals: 23 (.CHRG quirk, curcuma correct), 27/28 + 7/8 (torsion deficit, cur ~80% of xtb), 30/AB (bond pibo, ipis=0). Pre-existing from merge: `cli_curcumaopt_07_opt_multixyz` golden-value drift.

3. **Verbosity scoping across threads (residual)**: the global `CurcumaLogger` verbosity is now RAII-scoped via `CurcumaMethod` (ctor/dtor save-restore), but it is a shared static and cannot be cleanly scoped across `CxxThreadPool` workers — the pool-owning helpers (`PerformMolecularDynamics`/`PerformOptimisation`) and energy-method setup re-assert the level at their boundaries (kept deliberately); a `thread_local` verbosity (the only full fix) is out of scope. See [docs/CONFSEARCH_ROADMAP.md](docs/CONFSEARCH_ROADMAP.md) #1.

4. **ConfSearch Phase A-C**: efficiency/robustness features (RATTLE threshold, topo/Epot abort, seed funnel, opt→bias feedback, permutation-aware + adaptive MTD bias) — roadmap, open TODOs and experimental caveats in [docs/CONFSEARCH_ROADMAP.md](docs/CONFSEARCH_ROADMAP.md). Cross-run bias heating (shared-pool hills `W=k·counter` grow unbounded → `<T>` climbs run-by-run → NaN) is bounded by **defaults ON for ConfSearch**: `rmsd_mtd_freeze_inherited`+`temp_abort` (measured best: 0 blow-ups, best conformer yield; `rmsd_mtd_max_height` is opt-in for tighter T). The bare `-startT 500` run no longer blows up (roadmap TODO #4; intra-run wide-hill blow-up still open).

5. **Native GFN transition metals (MOR41, Jul 16, 2026, AI/machine-tested)**: **GFN2 3d metals FIXED** — `reference_occ`/`p_kcn`/`p_shpoly` are angular-momentum-indexed in tblite but were read by shell index; transition metals order shells `[d,s,p]` so this scrambled them (main-group unaffected). After the gfn2-only fix (`xtb_native.cpp`), **3d metals (Cr/Mn/Fe/Co/Ni) are exact to tblite (1e-8)**; MOR41 reaction MAD native-GFN2-vs-xtb dropped 507→4.29 kcal/mol; no main-group regression. **Still open**: GFN2 4d/5d residual (Mo/Ru/Rh/Pd ~1e-3 Eh, W/Ir/Pt ~0.1–0.5 Eh — heavy-element band energy, not an indexing bug); **GFN1 TMs unfixed** (its valence+polarisation shells share angular momenta, e.g. H=2 s-shells, so the ang-remap collapses them — needs a gfn1-specific fix); GFN-FF TM deviations (force-field params). xtb-gfn1 is itself SCF-unstable on this set (PR40 → spurious −93.6 Eh). Runner `scripts/mor41_validation.py`, data `test_cases/MOR41-testset/reactions.dat`, results [docs/MOR41_VALIDATION.md](docs/MOR41_VALIDATION.md).

