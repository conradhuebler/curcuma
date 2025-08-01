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

### Current Development
- Enhanced conformational search algorithms
- Improved trajectory analysis tools
- Better integration with quantum chemical methods

### Known Issues
- Memory usage optimization needed for large systems
- Enhanced error handling for failed optimizations

---

*This documentation covers all molecular modeling capabilities and applications*