# CLAUDE.md - Curcuma Development Guide

## Overview

**Curcuma** is a molecular modelling and simulation toolkit for computational chemistry, force fields, and quantum calculations.

**Educational Purpose**: Teaching and research platform prioritizing pedagogical clarity over complex software engineering. Goal: learn, understand, and implement computational chemistry methods without getting lost in C++ abstractions.

## General Instructions

- Each source code dir has a CLAUDE.md with basic information of the code and logic
- **Keep CLAUDE.md files concise and focused** - avoid verbose descriptions and redundant information
- Tasks corresponding to code must be placed in the correct CLAUDE.md file
- Each CLAUDE.md has a variable part (short-term info, bugs) and preserved part (permanent knowledge)
- **Instructions blocks** contain operator-defined future tasks and visions for code development
- Only include information important for ALL subdirectories in main CLAUDE.md
- Preserve new knowledge from conversations but keep it brief
- Always suggest improvements to existing code

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
- **Native GFN-FF** - Curcuma's own implementation (`cgfnff`) - *WORK IN PROGRESS*

### 2. Force Field Methods
- **Universal Force Field (UFF)** - General-purpose molecular mechanics
- **GFN-FF** - Geometry/Frequency/Noncovalent Force Field (via XTB and native)
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
Central hub for all energy and gradient calculations. Routes method calls to appropriate interfaces:

```cpp
// Method routing via SwitchMethod()
case 9: GFNFF (cgfnff)           // Native GFN-FF implementation
case 6: EHT                      // Extended Hückel Theory  
case 4: DFT-D3                   // Dispersion corrections
case 3: Ulysses                  // Semi-empirical methods
case 2: XTB                      // Extended tight-binding
case 1: TBLite                   // Tight-binding DFT
case 0: ForceField               // UFF, QMDFF, etc.
```

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
│   │   ├── energycalculator.cpp  # Energy/gradient dispatcher
│   │   ├── forcefield.cpp        # Force field engine
│   │   ├── forcefieldthread.cpp  # Parallel FF calculations
│   │   ├── molecule.cpp          # Molecular data structures
│   │   └── qm_methods/           # Quantum chemistry interfaces
│   │       ├── gfnff.cpp         # Native GFN-FF implementation
│   │       ├── eht.cpp           # Extended Hückel Theory
│   │       ├── xtbinterface.cpp  # XTB interface
│   │       └── tbliteinterface.cpp # TBLite interface
│   ├── tools/                # Utilities and file I/O
│   │   ├── formats.h         # File format handling (XYZ, MOL2, SDF)
│   │   └── geometry.h        # Geometric calculations
│   └── helpers/              # Development and testing tools
├── test_cases/               # Validation and benchmark molecules
├── external/                 # Third-party dependencies
└── CMakeLists.txt           # Build configuration
```

## Recent Developments (July 2025)

### Native GFN-FF Implementation
- **Status**: WORK IN PROGRESS - Architecture complete, debugging needed
- **Method name**: `cgfnff` (curcuma-gfnff) vs `gfnff` (XTB-gfnff)
- **Integration**: Fully integrated in EnergyCalculator as case 9
- **Current issue**: JSON parameter generation bug (null values)

### Universal Parameter Caching 
- **Status**: COMPLETED and WORKING
- **Auto-naming**: `input.xyz` → `input.param.json`
- **Method validation**: Prevents loading incompatible parameters
- **Multi-threading safety**: `setParameterCaching(false)` for concurrent access
- **Performance**: 96% speedup demonstrated

## Build and Test Commands

```bash
# Build
mkdir build && cd build
cmake .. && make -j4

# Test UFF (working)
./curcuma -sp input.xyz -method uff

# Test native GFN-FF (work in progress)  
./curcuma -sp input.xyz -method cgfnff
```

## Standards

### Logging System (`src/core/curcuma_logger.h/.cpp`)
**Verbosity Levels**: 0 (silent) → 1 (results only) → 2 (normal) → 3 (debug)
**Color Scheme**: Red (errors), Orange (warnings), Green (success/citations), Blue (parameters), Magenta (debug)

**Core Functions**:
- `CurcumaLogger::error()` - Always visible (red)
- `CurcumaLogger::warn()` - Level ≥1 (orange) 
- `CurcumaLogger::success()` - Level ≥1 (green)
- `CurcumaLogger::info()` - Level ≥2 (default)
- `CurcumaLogger::param()` - Level ≥2 (blue)
- `CurcumaLogger::citation()` - Level ≥2 (green)

**Debug System** (`CURCUMA_DEBUG` cmake option):
- `CurcumaLogger::debug()` - Compile-time eliminated in release builds
- Auto-color detection with `CURCUMA_NO_COLOR` environment override

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

1. **cgfnff JSON bug**: Parameter generation creates null values
2. **Missing real GFN-FF parameters**: Currently uses placeholder values  
3. **Logging system**: Needs implementation of new standards
4. **Unit migration**: Need to replace hardcoded constants throughout codebase with CurcumaUnit functions