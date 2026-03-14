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

### Configuration Management System
- **ConfigManager** (`config_manager.h/cpp`) - Claude Generated 2025
  - **Purpose**: Modern type-safe parameter access, eliminates Json2KeyWord boilerplate
  - **Architecture**: Wrapper around ParameterRegistry with hierarchical dot notation support
  - **API**: `config.get<T>("key")` with case-insensitive lookup and default value support
  - **Features**:
    - Automatic default merging
    - Hierarchical keys (`"topological.save_image"`)
    - Type safety
    - **Alias Resolution** (October 2025): Resolves aliases to canonical names via ParameterRegistry
    - **Case-Insensitive**: `-MaxTime`, `-maxtime` beide akzeptiert
  - **Status**: Production-ready, proof-of-concept in analysis.cpp (37 Json2KeyWord calls eliminated)
- **ParameterRegistry** (`parameter_registry.h/cpp`) - Claude Generated 2025
  - **Backend**: Stores all module parameters from build-time extraction
  - **Used By**: ConfigManager for default values and validation
  - **Alias Resolution** (October 2025): Case-insensitive alias lookup via `resolveAlias()`

### Force Field System
- **ForceField**: Main engine with universal JSON parameter caching (96% speedup) + **CurcumaLogger verbosity**
- **ForceFieldGenerator**: Parameter generation with **progress tracking and timing**
- **ForceFieldThread**: Multi-threaded calculations for large systems
- **Universal Verbosity**: Energy decomposition, timing analysis, silent mode support
- **Performance Critical**: Parameter caching essential for iterative calculations

### Physical Directory Structure (Completed Restructuring)
```cpp
core/
├── energycalculator.cpp        # NEW: Unified polymorphic dispatcher
├── molecule.cpp                # Core molecular data structures
├── curcuma_logger.cpp          # Universal logging system
├── config_manager.cpp/h        # NEW: Modern parameter access layer (Oct 2025)
├── parameter_registry.cpp/h    # Parameter registry backend (Oct 2025)
├── energy_calculators/         # NEW: All computational methods consolidated here
│   ├── computational_method.h      # Base interface for all methods
│   ├── method_factory.cpp          # Priority-based method creation
│   ├── qm_methods/                 # ALL QM methods (moved from src/core/qm_methods/)
│   │   ├── eht.cpp                 # Extended Hückel + CurcumaLogger
│   │   ├── xtbinterface.cpp        # XTB + synchronized verbosity
│   │   ├── tbliteinterface.cpp     # TBLite + synchronized verbosity
│   │   ├── ulyssesinterface.cpp    # Ulysses + CurcumaLogger
│   │   ├── gfnff.cpp               # Native GFN-FF (WIP)
│   │   ├── orcainterface.cpp       # ORCA interface
│   │   ├── dftd3interface.cpp      # DFT-D3 dispersion
│   │   ├── dftd4interface.cpp      # DFT-D4 dispersion
│   │   ├── *_method.cpp            # Polymorphic method wrappers
│   │   └── interface/              # Abstract interfaces
│   └── ff_methods/                 # ALL force field methods (moved from src/core/)
│       ├── forcefield.cpp          # Main FF engine + verbosity
│       ├── forcefieldgenerator.cpp # Parameter generation + progress tracking
│       ├── forcefieldthread.cpp    # Multi-threading support
│       ├── qmdff.cpp               # QMDFF implementation
│       ├── eigen_uff.cpp           # UFF implementation
│       └── *_par.h                 # Parameter databases (UFF, QMDFF)
└── (other core files...)           # topology.cpp, fileiterator.cpp, etc.
```

## Instructions Block

**PRESERVED - DO NOT EDIT BY CLAUDE**

## Variable Section

### Active Issues
- Memory optimization for large molecular systems (>1000 atoms)

### Performance Notes
- **Parameter caching**: 96% speedup for iterative calculations - critical for optimization/MD
- **Distance matrix caching**: 2-5x speedup for repeated distance/topology calculations
- **Thread safety**: Use `setParameterCaching(false)` for concurrent force field access
- **Large systems**: Optimized algorithms with bounds checking for molecules >1000 atoms

### Completed Developments ✅
- ✅ EnergyCalculator polymorphic refactoring (eliminates SwitchMethod)
- ✅ Universal Verbosity System (4-level output 0-3)
- ✅ MethodFactory with priority fallbacks (gfn2: TBLite → Ulysses → XTB)
- ✅ Parameter routing fix (Oct 26, 2025) - multi-module hierarchies now work
- ✅ ConfigManager type-safe parameter layer
- ✅ Physical architecture restructuring (energy_calculators/)

### Unit System (`units.h`)
- **Centralized constants**: CODATA-2018 values in `CurcumaUnit` namespace
- **Educational focus**: Clear function names and comprehensive documentation
- **Consistency**: Eliminates scattered unit definitions across codebase
- **Backward compatibility**: Legacy aliases for existing code