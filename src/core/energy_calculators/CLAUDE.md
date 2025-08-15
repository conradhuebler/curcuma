# CLAUDE.md - Energy Calculators Directory

## Overview

The energy_calculators directory contains the **new polymorphic architecture** for computational methods in Curcuma. This represents a complete refactoring of the previous EnergyCalculator system, implemented in January 2025.

## Architecture

### Core Design Principles

- **Single Interface**: All computational methods (QM and MM) implement the same `ComputationalMethod` interface
- **Polymorphism Over Templates**: Educational clarity through direct virtual function calls
- **Method Hierarchies**: Priority-based resolution with automatic fallbacks
- **API Compatibility**: All existing EnergyCalculator usage works unchanged
- **Thread Safety**: Enhanced concurrency support maintained

### Directory Structure

```
energy_calculators/
├── computational_method.h      # Base interface for all methods
├── method_factory.cpp/h        # Priority-based method creation
├── qm_methods/                 # Quantum method wrappers
│   ├── eht_method.cpp/h        # Extended Hückel Theory wrapper
│   ├── xtb_method.cpp/h        # XTB interface wrapper  
│   ├── tblite_method.cpp/h     # TBLite interface wrapper
│   └── ulysses_method.cpp/h    # Ulysses interface wrapper
└── ff_methods/                 # Force field wrappers
    └── forcefield_method.cpp/h # ForceField wrapper with threading
```

## Core Components

### ComputationalMethod Interface

**Base class for all computational methods** - provides unified API:

```cpp
class ComputationalMethod {
public:
    virtual ~ComputationalMethod() = default;
    virtual bool setMolecule(const Mol& mol) = 0;
    virtual double calculateEnergy(bool gradient = false) = 0;
    virtual Matrix getGradient() const = 0;
    virtual Vector getCharges() const = 0;
    virtual bool supportsGradients() const = 0;
    virtual std::string getMethodName() const = 0;
};
```

### MethodFactory System

**Priority-based method creation** with automatic fallbacks:

#### **Method Hierarchies**
```cpp
// gfn2 priority: TBLite → Ulysses → XTB
"gfn2" -> TBLiteMethod("gfn2") [preferred]
       -> UlyssesMethod("ugfn2") [fallback]  
       -> XTBMethod("gfn2") [last resort]

// gfn1 priority: TBLite → XTB → Ulysses
"gfn1" -> TBLiteMethod("gfn1") [preferred]
       -> XTBMethod("gfn1") [fallback]
       -> UlyssesMethod("ugfn1") [last resort]
```

#### **Method Creation**
```cpp
std::unique_ptr<ComputationalMethod> method = 
    MethodFactory::createMethod("gfn2", config);
// Automatically tries TBLite first, falls back to Ulysses, then XTB
```

## Method Implementations

### QM Method Wrappers

All QM methods are wrapped to provide the same interface while preserving their individual capabilities:

- **EHTMethod**: Wraps native Extended Hückel Theory implementation
- **XTBMethod**: Wraps XTB interface with synchronized verbosity
- **TBLiteMethod**: Wraps TBLite interface with context verbosity control  
- **UlyssesMethod**: Wraps Ulysses interface with SCF progress tracking

### Force Field Wrapper

**ForceFieldMethod** wraps the ForceField class while maintaining:
- **Multi-threading support** via ForceFieldThread
- **Parameter generation** integration with ForceFieldGenerator  
- **Universal caching** with 96% speedup for iterative calculations
- **Thread safety** controls for concurrent access

## Integration with EnergyCalculator

### Old vs New Architecture

#### **Old System (Pre-January 2025)**
```cpp
// Complex SwitchMethod with 10+ cases
switch(method_id) {
    case 0: /* ForceField code */ break;
    case 1: /* TBLite code */ break;  
    case 2: /* XTB code */ break;
    // ... 10+ more cases
}
```

#### **New System (January 2025)**
```cpp
// Single polymorphic method pointer
std::unique_ptr<ComputationalMethod> m_method;

// Simple method resolution
m_method = MethodFactory::createMethod(method_name, config);
return m_method->calculateEnergy(gradient);
```

### Benefits of New Architecture

1. **Eliminates Giant Switch**: No more 200+ line SwitchMethod function
2. **Automatic Fallbacks**: Method hierarchies with priority resolution  
3. **Enhanced Error Handling**: Method-specific error reporting
4. **Universal Verbosity**: Consistent CurcumaLogger integration
5. **Educational Clarity**: Direct polymorphic calls vs complex conditionals
6. **Maintainability**: Easy to add new methods without touching dispatcher

## Universal Verbosity Integration

All wrapped methods support the **4-level verbosity system**:

- **Level 0**: Silent mode (zero output for optimization/MD)
- **Level 1**: Final results only (energies, convergence)  
- **Level 2**: Scientific analysis (HOMO/LUMO, properties, energy decomposition)
- **Level 3**: Complete debug (full orbital listings, timing, algorithm details)

### Native Library Integration

- **XTB**: `XTB_VERBOSITY_MUTED/MINIMAL/FULL` synchronized with CurcumaLogger
- **TBLite**: Context verbosity (0/1/3) controlled by CurcumaLogger levels
- **Ulysses**: Direct CurcumaLogger integration for SCF progress
- **EHT**: Native CurcumaLogger implementation with orbital analysis

## Development Status

### ✅ Completed
- **Polymorphic Architecture**: Full implementation with all method wrappers
- **MethodFactory**: Priority-based resolution with hierarchical fallbacks
- **API Compatibility**: All existing EnergyCalculator usage preserved
- **Universal Verbosity**: Complete integration across all computational methods
- **Thread Safety**: Enhanced concurrency support maintained
- **Error Handling**: Comprehensive CurcumaLogger-based error reporting

### 🔧 Future Enhancements
- **Additional Methods**: Easy integration of new computational methods
- **Performance Optimization**: Method-specific performance tuning
- **Extended Hierarchies**: More complex fallback chains for specialized methods

---

*This directory represents the future of computational method integration in Curcuma - providing educational clarity through polymorphism while maintaining high performance and scientific accuracy.*