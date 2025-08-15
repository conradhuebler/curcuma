# CLAUDE.md - Core Directory

## Overview

Core computational engines: energy calculations, force fields, molecular data, quantum interfaces.

## Key Components

### Energy Calculator (`energycalculator.cpp/h`)
**🚀 COMPLETELY REFACTORED (January 2025)** - Unified polymorphic architecture

#### **New Architecture**
- **Single Method Pointer**: `std::unique_ptr<ComputationalMethod> m_method`
- **Polymorphic Interface**: All QM/MM methods through same `calculateEnergy()` API
- **MethodFactory Integration**: Automatic method creation with priority fallbacks
- **Thread-Safe**: Enhanced concurrency support maintained
- **API Compatible**: All existing EnergyCalculator calls work unchanged

#### **Method Resolution (New)**
```cpp
// Old: SwitchMethod() with 10+ cases
// New: Single polymorphic method pointer
m_method = MethodFactory::createMethod(method_name, config);
return m_method->calculateEnergy(gradient);
```

#### **Benefits**
- **Eliminates SwitchMethod()**: No more giant switch statements
- **Automatic Fallbacks**: `gfn2` tries TBLite → Ulysses → XTB
- **Consistent API**: Same interface for all computational methods
- **Enhanced Error Handling**: Method-specific error reporting
- **Universal Verbosity**: Integrated CurcumaLogger support

### Molecule Class (`molecule.cpp/h`)
**Educational Focus**: Core molecular data structure - direct, minimal abstractions
- **Core Data**: Atoms, coordinates, bonds, charges, fragments
- **File I/O**: XYZ, MOL2, SDF formats integrated
- **Geometry Operations**: Distance, angle calculations, structure manipulation
- **Performance Note**: Large molecules (>1000 atoms) may need memory optimization

### Force Field System
- **ForceField**: Main engine with universal JSON parameter caching (96% speedup) + **CurcumaLogger verbosity**
- **ForceFieldGenerator**: Parameter generation with **progress tracking and timing**
- **ForceFieldThread**: Multi-threaded calculations for large systems
- **Universal Verbosity**: Energy decomposition, timing analysis, silent mode support
- **Performance Critical**: Parameter caching essential for iterative calculations

### New Directory Structure
```cpp
core/
├── energycalculator.cpp        # NEW: Unified polymorphic dispatcher
├── energy_calculators/         # NEW: Polymorphic method wrappers  
│   ├── computational_method.h  # Base interface for all methods
│   ├── method_factory.cpp      # Priority-based method creation
│   ├── qm_methods/            # QM method wrappers (EHT, XTB, TBLite, Ulysses)
│   └── ff_methods/            # Force field wrappers
├── curcuma_logger.cpp          # Universal logging system
├── forcefield.cpp              # + verbosity integration
├── forcefieldgenerator.cpp     # + progress tracking
└── qm_methods/                # Native implementations + verbosity
    ├── eht.cpp                # Extended Hückel + CurcumaLogger
    ├── xtbinterface.cpp       # XTB + synchronized verbosity
    ├── tbliteinterface.cpp    # TBLite + synchronized verbosity  
    └── ulyssesinterface.cpp   # Ulysses + CurcumaLogger
```

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

### Major Refactoring Completed (January 2025) ✅
- **🚀 EnergyCalculator Architecture**: Complete polymorphic refactoring - eliminates SwitchMethod  
- **🎯 Universal Verbosity System**: All QM/MM methods support 4-level verbosity (0-3)
- **🏗️ MethodFactory Integration**: Priority-based method resolution with automatic fallbacks
- **🔧 Enhanced Error Handling**: Method-specific error reporting with CurcumaLogger
- **📈 Performance Maintained**: Threading, caching, and optimization preserved
- **✅ API Compatibility**: All existing EnergyCalculator usage works unchanged

### Previous Improvements
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