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
- **Format**: `Copyright (C) 2019 - 2025 Conrad Hübler <Conrad.Huebler@gmx.net>`
- **AI acknowledgment**: Add Claude contribution notes in code comments, not copyright headers

#### Code Structure Guidelines
- **Main computational functions**: Should be easily findable and readable without deep C++ knowledge
- **Algorithm documentation**: Include mathematical formulations and literature references
- **Parameter transparency**: Make method parameters and their physical meaning obvious
- **Debugging accessibility**: Provide easy ways to inspect intermediate results and algorithm steps
- **Educational examples**: Include well-commented example usage in documentation

## Current Capabilities

### 1. Quantum Mechanical Methods
- **Extended Hückel Theory (EHT)** - Semi-empirical quantum chemistry
- **TBLite Interface** - Tight-binding DFT methods (GFN1, GFN2, iPEA1)
- **XTB Interface** - Extended tight-binding methods (GFN-FF, GFN1, GFN2)
- **Ulysses Interface** - Various semi-empirical methods (PM3, AM1, MNDO, etc.)
- **Native GFN-FF** - Curcuma's own implementation (`cgfnff`) - ✅ **IMPLEMENTED**

### 2. Force Field Methods
- **Universal Force Field (UFF)** - General-purpose molecular mechanics
- **GFN-FF** (`cgfnff`) - ✅ **FULLY IMPLEMENTED** - See [docs/GFNFF_STATUS.md](docs/GFNFF_STATUS.md)
- **QMDFF** - Quantum Mechanically Derived Force Fields
- **Universal Parameter Caching** - Automatic save/load for all FF methods

### 3. Dispersion and Non-Covalent Corrections
- **DFT-D3** - Grimme's D3 dispersion correction
- **DFT-D4** - Next-generation D4 dispersion correction  
- **H4 Correction** - Hydrogen bonding and halogen bonding corrections

### 4. Geometry Optimization
- **LBFGS Optimizer** - Limited-memory Broyden-Fletcher-Goldfarb-Shanno
- **Multiple Convergence Criteria** - Energy, gradient, RMSD-based
- **Constrained Optimization** - Distance, angle, and dihedral constraints

### 5. Conformational Analysis
- **ConfSearch** - Systematic conformational searching
- **ConfScan** - Conformational scanning along reaction coordinates
- **RMSD Analysis** - Structure comparison and alignment
- **Energy-based Filtering** - Automatic conformer ranking

### 6. Molecular Dynamics
- **SimpleMD** - Basic molecular dynamics simulation
- **NEB Docking** - Nudged elastic band for transition states
- **Trajectory Analysis** - Analysis of MD trajectories

### 7. Analysis Tools
- **RMSD Calculations** - Root-mean-square deviation analysis
- **Persistent Diagram** - Topological data analysis
- **Hessian Analysis** - Second derivative calculations
- **Orbital Analysis** - Molecular orbital visualization and analysis

## Architecture

### Core Components

#### Energy Calculator (`src/core/energycalculator.h/cpp`)
**COMPLETELY REFACTORED (January 2025)** - New unified polymorphic architecture:

##### **New Architecture**
- **Polymorphic Design**: Single `ComputationalMethod` interface for all QM/MM methods
- **MethodFactory**: Priority-based method resolution with hierarchical fallbacks
- **Unified Interface**: `calculateEnergy()`, `getGradient()`, consistent API across all methods
- **Method Priority System**: `gfn2` → TBLite > Ulysses > XTB (automatic fallback)
- **Thread-Safe**: Full multi-threading support maintained
- **Universal Verbosity**: Consistent output control across all computational methods

##### **Method Resolution (New)**
```cpp
// New MethodFactory system (replaces old SwitchMethod)
std::unique_ptr<ComputationalMethod> method = 
    MethodFactory::createMethod("gfn2", config);
double energy = method->calculateEnergy();
```

##### **Supported Method Hierarchies**
- **gfn2**: TBLite → Ulysses → XTB (automatic priority resolution)
- **gfn1**: TBLite → XTB → Ulysses  
- **uff**: ForceField wrapper with parameter generation
- **eht**: Extended Hückel Theory (native implementation)
- **cgfnff**: Native GFN-FF (✅ **COMPLETE**)

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

## Completed Developments (2025)

✅ **Platform-Independent External Dependency Discovery** - Phase 2a+2b complete: Plumed2 via find_library(), D4 via find_package(LAPACK), both with fallback support, portable across Linux/macOS/Windows
✅ **External Dependency Conditional Compilation** - Phase 1b complete: D3/D4 guards in QMDFF/UFF/ForceField, MethodFactory runtime checks, all 14 external libs properly gated
✅ **Parameter Registry System** - Macro-based parameter definitions, auto-help, type validation
✅ **ConfigManager Layer** - Type-safe parameter access, hierarchical dot notation
✅ **EnergyCalculator Refactoring** - Polymorphic interface, priority-based method resolution
✅ **Universal Verbosity System** - Consistent 4-level output across all methods (0-3)
✅ **Method Hierarchies** - `gfn2` auto-falls back: TBLite → Ulysses → XTB
✅ **Physical Architecture** - QM/MM methods organized under `src/core/energy_calculators/`
✅ **Topological Data Analysis** - dMatrix legacy functionality integrated as TDAEngine
✅ **Parameter Routing Fix** - Multi-module parameter hierarchies now work (json null-error fixed)
✅ **GFN-FF Implementation** - Complete and operational in ff_methods/ - See [docs/GFNFF_STATUS.md](docs/GFNFF_STATUS.md)
✅ **EEQ Phase 1 Full Implementation** (December 28, 2025) - Complete dxi calculation with pi-system detection, neighbor EN averaging, and environment-dependent corrections; 75% reduction in charge error (5.0× → 1.3×); fixes 4/6 GFN-FF energy terms automatically

## Build and Test Commands

```bash
# Build (use release/ for stable builds)
cd release
make -j4

# Or create new build directory
mkdir build && cd build
cmake .. && make -j4

# Run all tests (NEW: October 2025 - CLI Tests added)
ctest --output-on-failure

# Run specific test categories
ctest -R "cli_rmsd_" --output-on-failure      # RMSD CLI tests (5/6 passing)
ctest -R "cli_confscan_" --output-on-failure  # ConfScan CLI tests (6/7 passing)
ctest -R "cli_simplemd_" --output-on-failure  # SimpleMD CLI tests (1/7 passing - opt bug)
ctest -R "cli_curcumaopt_" --output-on-failure # Opt CLI tests (1/6 passing - JSON bug)

# Run individual CLI test with verbose output
ctest -R "cli_rmsd_01" --verbose

# Legacy: Manual curcuma execution
./curcuma -sp input.xyz -method uff           # UFF single point
./curcuma -rmsd ref.xyz target.xyz            # RMSD calculation
./curcuma -opt input.xyz -method gfn2         # GFN2 optimization (currently broken)
```

**Test Status** (Stand 2025-10-26):
- **26 CLI Tests** implementiert in `test_cases/cli/`
- **19/26 bestanden** (73%) - RMSD: 6/6 ✅, ConfScan: 7/7 ✅, Curcumaopt: 6/6 ✅, SimpleMD: 0/7 (JSON crash blocking)

## Project Management

- **Prioritized TODO List**: See [TODO.md](TODO.md) - 20 tasks extracted from documentation audit
- **Critical**: SimpleMD JSON null crash blocks 7 tests - parameter routing needs fix (see TODO.md:CRITICAL)
- **Test Blocking Issues**: Memory optimization for large systems
- **Module Docs**: Each `src/` subdirectory has CLAUDE.md with specific tasks
- **GFN-FF Status**: See [docs/GFNFF_STATUS.md](docs/GFNFF_STATUS.md) for implementation details

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
- **Global Parameters**: `verbosity`, `threads` are additionally duplicated at top-level

## Planned Development

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

1. **GFN-FF Large Molecule Coulomb Deviations**: While small systems match Fortran exactly, larger systems show non-proportional deviations (e.g., Triose +26 mEh, Complex +6.5 mEh). See [docs/GFNFF_STATUS.md](docs/GFNFF_STATUS.md).

2. **GFN-FF Limitations**: See [docs/GFNFF_STATUS.md](docs/GFNFF_STATUS.md#known-limitations) for details (D4 dispersion, EEQ integration, metal parameters)

3. **Unit migration**: Some legacy code still uses hardcoded constants instead of CurcumaUnit functions

## Recently Resolved ✅

- ✅ **EEQ Solver Numerical Stability** (Jan 29, 2026): Replaced PartialPivLU with ColPivHouseholderQR + iterative refinement for better numerical stability on large systems. Investigation confirmed triose deviation is parameter-related, not numerical.
- ✅ **GFN-FF Coulomb Exact Match** (Jan 29, 2026): Fixed dgam correction factors and enabled charge-corrected parameters (gameeq, alpeeq); Coulomb energy now matches Fortran reference within < 1 nEh (0.0001% error) for small systems.
- ✅ **EEQ Charge Fine-Tuning** (Jan 17, 2026): Aligned parameters with Fortran reference; 82% Coulomb energy error reduction (8.3% → 1.5%)
- ✅ **GFN-FF HB/XB Refinement** (Jan 17, 2026): Implemented context-dependent basicity/acidity overrides for non-covalent interactions
- ✅ **D4 Dispersion Weighting Fix** (Jan 17, 2026): Corrected to CN-only weighting for GFN-FF hybrid model; D4 set as default
- ✅ **JSON Null-Error & Parameter Routing** (October 2025): SimpleMD/curcumaopt fixed - 11 tests now pass
- ✅ **ForceField Inversion Bug** (October 2025): Vector bounds crash in UFF parameter generation fixed
- ✅ **CLI Test Infrastructure** (October 2025): 26 End-to-End validation tests with scientific accuracy
- ✅ **GFN-FF Architecture** (December 2025): Implementation complete and tests passing - [docs/GFNFF_STATUS.md](docs/GFNFF_STATUS.md)
