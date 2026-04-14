# CLAUDE.md - Capabilities Directory

## Overview

The capabilities directory contains high-level molecular modeling applications and computational chemistry tasks. These modules use the core computational engines to perform complex multi-step calculations and analyses.

## Structure

```
capabilities/
├── confscan.cpp/h         # Conformational scanning along reaction coordinates
├── confsearch.cpp/h       # Systematic conformational searching
├── curcumaopt.cpp/h       # Geometry optimization algorithms
├── simplemd.cpp/h         # Molecular dynamics simulation
├── rmsd.cpp/h            # Structure comparison and alignment
├── rmsdtraj.cpp/h        # Trajectory RMSD analysis
├── hessian.cpp/h         # Second derivative calculations
├── persistentdiagram.cpp/h # Topological data analysis
├── optimiser/            # Optimization algorithms
│   ├── lbfgs.cpp/h       # LBFGS optimization
│   └── LevMar*.h         # Levenberg-Marquardt variants
└── c_code/               # C interface code (Hungarian algorithm)
```

## Key Capabilities

### Conformational Analysis
- **ConfScan**: Systematic scanning of conformational space along defined coordinates
- **ConfSearch**: Automated conformational searching with energy filtering
- **RMSD Analysis**: Structure comparison, alignment, and trajectory analysis

### Optimization

> **All non-external optimizers below are 🤖 AI-generated. None are ✅ TESTED or ✅ APPROVED.**

- **CurcumaOpt**: Geometry optimization dispatcher — legacy system, human-tested
- **LBFGSpp**: Wrapper around external LBFGSpp library — external code, wrapper is AI-generated
- **ANCOPT** (`ancopt_optimizer.cpp/h`): 🤖 AI-generated port of XTB's AncOpt (Grimme)
  - ⚙️ Machine-tested: water, ethane with UFF converge
  - Not tested: QM gradients, large systems, transition states, linear molecules
  - Not implemented vs. reference: exact XTB ANC regeneration thresholds, micro-cycle logic
  - Known issue: crash in cleanup phase (signal 11, under investigation)
- **OptimizerDriver** (`optimizer_driver.cpp/h`): 🤖 AI-generated base class (Template Method)
- **OptimizerFactory** / **OptimizationDispatcher**: 🤖 AI-generated
- **Native L-BFGS / DIIS / RFO** (`optimisation/lbfgs.cpp`): 🤖 AI-generated — see `optimisation/CLAUDE.md`
- **Constrained optimization**: Fix atomic positions by setting gradient = 0

### Molecular Dynamics
- **SimpleMD**: Basic molecular dynamics with various thermostats
  - ✅ Coarse-graining support with automatic system detection
  - ✅ PBC wrapping for periodic boundary conditions
  - ✅ 10x timestep scaling for pure CG systems
  - ✅ VTF trajectory output for CG systems
  - ✅ Orientational dynamics infrastructure (prepared for Phase 6 ellipsoids)
- **NEB Docking**: Nudged elastic band for transition state searches
- **Trajectory Analysis**: Analysis tools for MD trajectories

### Advanced Analysis
- **Hessian**: Second derivative calculations for normal modes
- **Persistent Diagrams**: Topological data analysis for molecular structures
- **Enhanced TDA (dMatrix replacement)**: Complete topological data analysis with TDAEngine
- **Pairmapper**: Advanced structure matching algorithms

## Development Guidelines

### Interface Design
- All capabilities use the core EnergyCalculator for energy/gradient evaluations
- Consistent parameter handling through JSON-based configuration
- Progress reporting and result persistence for long calculations

### Performance Considerations
- Multi-threading support where applicable
- Memory-efficient handling of large trajectory data
- Automatic checkpointing for resumable calculations

## Instructions Block

**PRESERVED - DO NOT EDIT BY CLAUDE**

*Future development tasks and visions to be defined by operator/programmer*

## Variable Section

### Parameter System - ConfigManager Layer
✅ **Status**: Production-ready across 4+ modules (analysis, rmsd, simplemd, confscan)
- ConfigManager provides type-safe parameter access with hierarchical dot notation
- Multi-module parameter routing fixed (Oct 26, 2025)
- Migration ongoing: Replace Json2KeyWord calls with `config.get<T>("param")`
- Reference: `src/core/config_manager.h`, example: `analysis.cpp`

### Current Development
- Enhanced conformational search algorithms
- Improved trajectory analysis tools
- Better integration with quantum chemical methods
- ⚙️ **IMPLEMENTED** (AI): Parameter Registry System, dMatrix Integration, RMSD Code Restructuring
- ⚙️ **IMPLEMENTED** (AI, Oct 29, 2025): Help System Dynamic Generation
- ⚙️ **IMPLEMENTED** (AI, Oct 30, 2025): SimpleMD CG Integration Phase 1-4
- ⚙️ **IMPLEMENTED** (AI, Nov 2025): SimpleMD CG Integration Phase 5 - VTF trajectory output
- ⚙️ **IMPLEMENTED** (AI, Jan 2026): Analysis Output Refactoring - Registry-based handler architecture
- Pending: Unit system migration, RMSD Strategy pattern (Phase 3), CG Phase 6 (ellipsoidal extensions)

### New Analysis Output Architecture
- **New Analysis Output Architecture**: Handler-based system with registry pattern
- **File Naming Schema**: basename.general.csv, basename.NNN.type.csv, basename.type_statistics.csv
- **Extensible Design**: IAnalysisOutputHandler interface for new analysis types
- **Benefits**: Eliminates duplication, single point of change, automatic file generation

### Known Issues
- Memory optimization needed for large systems (>1000 atoms)
- ConfScan verbosity: Accept/Reject messages not visible at default level
- SimpleMD wall potential: Boundary logic and force calculation accuracy issues (TODO)

### Unported Features from Old CurcumaOpt (TODO)

The following features existed in `curcumaopt.cpp` but are **not yet ported** to the new `OptimizerDriver`/`OptimizationDispatcher` system:

1. **Parallel batch optimization** (`curcumaopt.cpp:276`) — `ProcessMolecules()` used `CxxThreadPool` (SPThread/OptThread) to optimize multiple molecules in parallel. The new system processes batches serially (loop over single-molecule calls). Implement `OptimizationBatchRunner` using CxxThreadPool.

2. **Hydrogen-only optimization** (`opt_h`, `curcumaopt.cpp:521`) — per-atom constraints derived from element type (element==1 → movable). Useful for optimizing H positions while heavy atoms are fixed. Needs integration into `OptimizerDriver` constraint system.

3. **Hessian after optimization** (`hessian=1/2`, `curcumaopt.cpp:237,323`) — computes normal modes using `Hessian` class after optimization, saves `hessian.json` and `scf.json`. Port by calling `Hessian` class inside `executeOptimization()` in `main.cpp` when `hessian` parameter is set.

4. **Molecular orbital diagram** (`mo_scheme`, `curcumaopt.cpp:404`) — `WriteMO()` generates TikZ LaTeX orbital diagram, `WriteMOAscii()` ASCII variant. Parameters: `mo_homo`, `mo_lumo`, `mo_scale`. Needs EnergyCalculator `Energies()` / `NumElectrons()` results plumbed through result struct.

5. **"stop" file interrupt** (`curcumaopt.cpp:665`) — checks for `./stop` file on disk during optimization loop to allow early termination. Simple to add to `OptimizerDriver::Optimize()` loop.

6. **`fusion` mode** (`curcumaopt.cpp:688`) — skips `Molecule::Check()` validity gate, needed for unusual bonding / fusion compounds. Add `bool fusion` flag to `OptimizationContext` and check in driver loop.

7. **Dipole moment output for GFN2** (`curcumaopt.cpp:395`) — printed via `EnergyCalculator::Dipole()` after SP/opt when method is gfn2. Not available in current `OptimizationResult` struct.

**Priority**: Items 1 (parallel batch) and 2 (opt_h constraints) are most commonly used. Items 3-7 are lower priority.

---

*This documentation covers all molecular modeling capabilities and applications*