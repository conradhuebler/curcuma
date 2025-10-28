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
- **CurcumaOpt**: Geometry optimization using various algorithms
- **LBFGS**: Limited-memory Broyden-Fletcher-Goldfarb-Shanno optimizer
- **Constrained optimization**: Fix atomic positions by setting gradient = 0

### Molecular Dynamics
- **SimpleMD**: Basic molecular dynamics with various thermostats
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
- ✅ **COMPLETED**: Parameter Registry System, dMatrix Integration, RMSD Code Restructuring
- ✅ **COMPLETED** (Oct 29, 2025): Help System Dynamic Generation - all modules use ParameterRegistry for help
- Pending: Unit system migration, RMSD Strategy pattern (Phase 3)

### Known Issues
- Memory optimization needed for large systems (>1000 atoms)
- ConfScan verbosity: Accept/Reject messages not visible at default level
- SimpleMD wall potential: Boundary logic and force calculation accuracy issues (TODO)

---

*This documentation covers all molecular modeling capabilities and applications*