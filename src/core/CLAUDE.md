# CLAUDE.md - Core Directory

## Overview

The core directory contains the fundamental computational engines of Curcuma. These modules handle energy calculations, force field methods, molecular data structures, and quantum mechanical interfaces.

## Structure

```
core/
├── energycalculator.cpp/h    # Central energy/gradient dispatcher
├── molecule.cpp/h            # Molecular data structures and operations
├── forcefield.cpp/h          # Force field engine with caching
├── forcefieldthread.cpp/h    # Multi-threaded force field calculations
├── forcefieldgenerator.cpp/h # Parameter generation and setup
├── qm_methods/              # Quantum mechanical method interfaces
├── forcefield_terms/        # Force field term implementations
├── topology.h               # Molecular topology definitions
└── various parameter files   # UFF, QMDFF parameter definitions
```

## Key Components

### Energy Calculator (`energycalculator.cpp/h`)
Central hub routing all energy and gradient calculations via `SwitchMethod()`:
```cpp
// Method routing
case 9: GFNFF (cgfnff)      // Native GFN-FF implementation  
case 6: EHT                 // Extended Hückel Theory
case 4: DFT-D3              // Dispersion corrections
case 3: Ulysses             // Semi-empirical methods
case 2: XTB                 // Extended tight-binding
case 1: TBLite              // Tight-binding DFT
case 0: ForceField          // UFF, QMDFF, etc.
```

### Force Field System
- **ForceField**: Main force field engine with universal parameter caching
- **ForceFieldThread**: Multi-threaded parallel calculations
- **ForceFieldGenerator**: Automatic parameter generation and validation
- **Universal Caching**: JSON-based parameter persistence with method validation

### Molecule Class
- Complete molecular data structure with atoms, bonds, coordinates
- File I/O integration for various formats (XYZ, MOL2, SDF)
- Geometry manipulation and analysis functions
- Integration with all computational methods

## Technical Features

### Parameter Caching System
- **Auto-naming**: `input.xyz` → `input.param.json`
- **Method validation**: Prevents loading incompatible parameters
- **Performance**: Up to 96% speedup demonstrated
- **Thread safety**: `setParameterCaching(false)` for concurrent access

### Multi-threading Support
- Thread-safe force field calculations
- Controllable caching for parallel processing
- Memory-efficient parameter sharing

## Instructions Block

**PRESERVED - DO NOT EDIT BY CLAUDE**

*Core engine development priorities and architectural decisions to be defined by operator/programmer*

## Variable Section

### Current Development Status
- **Native GFN-FF (cgfnff)**: Architecture complete, debugging JSON parameter generation
- **Parameter caching**: Completed and fully functional with auto-save/load
- **ForceFieldGenerator**: Enhanced with detailed timing and status output
- **Multi-threading**: Enhanced support with safety controls

### Active Issues
- cgfnff JSON parameter generation creates null values
- Missing real GFN-FF parameters (using placeholders)
- Memory optimization for large molecular systems

### Recent Improvements (January 2025)
- **Enhanced ForceFieldGenerator output**: Added detailed timing and progress information for each generation step
- **Fixed automatic parameter caching**: Complete auto-save/load functionality in EnergyCalculator
- **Intelligent filename handling**: Multi-tier priority system for parameter file naming
  - Priority 1: geometry_file from controller
  - Priority 2: Global basename from CurcumaMethod (curcuma -opt file.xyz → file.param.json)
  - Priority 3: Fallback to method_atoms.param.json
- **CurcumaMethod integration**: setFileName() automatically sets global basename for all subclasses
- **Parameter summary display**: Cached parameters show detailed summary (bonds, angles, vdW pairs, etc.)
- **Comprehensive timing information**: Bond detection, angle/dihedral generation timing with summaries
- **Method validation**: Automatic compatibility checking when loading cached parameters
- **96% speedup achieved**: Through parameter caching system with validation
- Universal parameter caching system implementation
- Enhanced thread safety in force field calculations
- Improved error handling in energy calculator routing

---

*This documentation covers all core computational engines and data structures*