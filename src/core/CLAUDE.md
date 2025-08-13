# CLAUDE.md - Core Directory

## Overview

Core computational engines: energy calculations, force fields, molecular data, quantum interfaces.

## Key Components

### Energy Calculator (`energycalculator.cpp/h`)
Central dispatcher routing all calculations via `SwitchMethod()`:
```cpp
case 9: GFNFF (cgfnff)      // Native GFN-FF
case 6: EHT                 // Extended Hückel Theory
case 4: DFT-D3              // Dispersion corrections
case 3: Ulysses             // Semi-empirical methods
case 2: XTB                 // Extended tight-binding
case 1: TBLite              // Tight-binding DFT
case 0: ForceField          // UFF, QMDFF, etc.
```

### Molecule Class (`molecule.cpp/h`)
**Educational Focus**: Core molecular data structure - direct, minimal abstractions
- **Core Data**: Atoms, coordinates, bonds, charges, fragments
- **File I/O**: XYZ, MOL2, SDF formats integrated
- **Geometry Operations**: Distance, angle calculations, structure manipulation
- **Performance Note**: Large molecules (>1000 atoms) may need memory optimization

### Force Field System
- **ForceField**: Main engine with universal JSON parameter caching (96% speedup)
- **ForceFieldThread**: Multi-threaded calculations for large systems
- **Performance Critical**: Parameter caching essential for iterative calculations

## Instructions Block

**PRESERVED - DO NOT EDIT BY CLAUDE**

## Variable Section

### Active Issues
- cgfnff JSON parameter generation creates null values
- Missing real GFN-FF parameters (using placeholders)
- Memory optimization for large molecular systems (>1000 atoms)

### Performance Notes
- **Parameter caching**: 96% speedup for iterative calculations - critical for optimization/MD
- **Distance matrix caching**: 2-5x speedup for repeated distance/topology calculations
- **Thread safety**: Use `setParameterCaching(false)` for concurrent force field access
- **Large systems**: Optimized algorithms with bounds checking for molecules >1000 atoms

### Recent Status (January 2025)
- **Scientific fixes**: Dipole moments now use center of mass (physically correct)
- **Performance caching**: Distance matrices cached with automatic invalidation
- **Centralized unit system**: All conversions moved to `CurcumaUnit` namespace in `units.h`
- **Algorithm optimization**: Improved angle and distance calculations with bounds checking
- **Documentation**: Added scientific references and formula documentation for learning

### Unit System (`units.h`)
- **Centralized constants**: CODATA-2018 values in `CurcumaUnit` namespace
- **Educational focus**: Clear function names and comprehensive documentation
- **Consistency**: Eliminates scattered unit definitions across codebase
- **Backward compatibility**: Legacy aliases for existing code