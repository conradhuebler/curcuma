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
- **ConfScan verbosity issue**: Accept/Reject messages not showing at default verbosity level - need to investigate default verbosity settings and ensure user feedback is always visible

### TODO - SimpleMD Wall Potential Fixes (January 2025)
**CRITICAL CORRECTNESS ISSUES IDENTIFIED:**

1. **Fix boundary logic error** (lines 814-829 in simplemd.cpp)
   - Current: `if (m_wall_x_min - m_wall_x_max < 1)` always triggers auto-config
   - Should: `if (m_wall_x_max - m_wall_x_min <= 0)` to detect invalid bounds
   - Impact: User-defined boundaries are ignored, always auto-configured

2. **Fix rectangular harmonic force calculation** (lines 2430-2434)
   - Current: Uses `std::abs(x - x_min)` which gives wrong force direction
   - Should: `(x - x_min)` without abs() for correct restoring force
   - Impact: Atoms pushed in wrong direction, simulation instability

3. **Fix log-Fermi force denominator** (lines 2318-2320, 2357-2359)
   - Current: Uses `(1 - exp_expr)` in denominator
   - Should: `(1 + exp_expr)` - correct derivative of log(1 + e^x)
   - Impact: Wrong force magnitudes, potential numerical instability

4. **Standardize parameter meanings**
   - Current: `wall_temp` means force constant k OR temperature kbT
   - Current: `wall_beta` unused in harmonic, critical for log-Fermi
   - Should: Consistent parameter interpretation across wall types
   - Add parameter validation and clear documentation

5. **Add numerical stability safeguards**
   - Check for distance = 0 in spherical calculations
   - Prevent exponential overflow in log-Fermi potentials
   - Add parameter bounds checking

6. **Add comprehensive wall potential tests**
   - Unit tests for force calculations vs analytical derivatives
   - Boundary condition tests (atoms at/beyond walls)
   - Parameter validation tests
   - Comparison with reference implementations

**Priority**: HIGH - These are correctness bugs that cause wrong physics
**Status**: IDENTIFIED - Analysis complete, fixes needed

---

*This documentation covers all molecular modeling capabilities and applications*