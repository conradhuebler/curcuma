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

## ConfigManager Integration (October 2025)

**Complete parameter management modernization** - All computational methods now use type-safe ConfigManager system.

### Parameter Flow Architecture

The ConfigManager system provides **end-to-end type-safe parameter passing** from user input through the polymorphic architecture to low-level interfaces:

```
User CLI/JSON
    ↓
Capability (e.g., Opt, SimpleMD)
    ↓ Creates ConfigManager("opt", controller)
EnergyCalculator(method, ConfigManager&)
    ↓ Delegates to ConfigManager constructor
    ↓ Exports JSON for backward compatibility
MethodFactory::createMethod(method, json)
    ↓ Creates method-specific wrapper
Method Wrapper (e.g., XTBMethod, TBLiteMethod)
    ↓ Creates ConfigManager("xtb", json)
QM/FF Interface (e.g., XTBInterface, ForceFieldGenerator)
    ↓ Type-safe access: config.get<int>("accuracy")
External Library (XTB, TBLite, ForceField engine)
```

### Implementation Details

#### **EnergyCalculator Constructor (Phase 3C - Oktober 2025)**

**Delegating Constructor Pattern** for backward compatibility:

```cpp
// Modern ConfigManager constructor
EnergyCalculator::EnergyCalculator(const std::string& method, const ConfigManager& config)
    : m_method_name(method)
{
    m_controller = config.exportConfig();  // Convert to JSON for compatibility
    createMethod(method, m_controller);
}

// Old JSON constructor delegates to ConfigManager version
EnergyCalculator::EnergyCalculator(const std::string& method, const json& controller)
    : EnergyCalculator(method, ConfigManager("energycalculator", controller))
{
}
```

**Benefits**:
- Single implementation point (no code duplication)
- All 12 capabilities work unchanged
- Seamless JSON ↔ ConfigManager conversion

#### **Method Wrapper Integration (Phase 3B)**

All method wrappers create ConfigManager instances for their interfaces:

```cpp
// Example: TBLiteMethod wrapper
TBLiteMethod::TBLiteMethod(const std::string& method_name, const json& config)
{
    ConfigManager tblite_config("tblite", config);  // Extract tblite-specific params
    m_tblite = std::make_unique<TBLiteInterface>(tblite_config);
}

// Example: XTBMethod wrapper
XTBMethod::XTBMethod(const std::string& method_name, const json& config)
{
    ConfigManager xtb_config("xtb", config);  // Extract xtb-specific params
    m_xtb = std::make_unique<XTBInterface>(xtb_config);
}
```

#### **QM/FF Interface Constructors (Phase 3B)**

All interfaces accept ConfigManager for type-safe parameter access:

```cpp
// Example: XTBInterface constructor
XTBInterface::XTBInterface(const ConfigManager& config)
    : m_config(config)
{
    // Type-safe parameter access with defaults from PARAM definitions
    m_accuracy = m_config.get<int>("accuracy", 2);
    m_SCFmaxiter = m_config.get<int>("max_iterations", 100);
    m_Tele = m_config.get<double>("electronic_temperature", 300.0);
    m_spin = m_config.get<double>("spin", 0.0);
}
```

### Migration Status (October 2025)

| Component | Status | Details |
|-----------|--------|---------|
| **Phase 3A**: Parameter Definitions | ✅ COMPLETE | 240 parameters across 12 modules with PARAM macros |
| **Phase 3B**: Interface Constructors | ✅ COMPLETE | All 8 QM/FF interfaces accept ConfigManager |
| **Phase 3B**: Method Wrappers | ✅ COMPLETE | All 8 wrappers create ConfigManager for interfaces |
| **Phase 3C**: EnergyCalculator Integration | ✅ COMPLETE | Delegating constructors, full backward compatibility |

**Remaining TODOs**:
- DFT-D3/D4 `UpdateParameters()` methods still use JSON (low priority, seldom used)
- Native GFN-FF parameter generation completeness (theoretical work needed)

### Key Benefits

1. **Type Safety**: `config.get<int>("accuracy")` catches type errors at compile time
2. **Single Source of Truth**: PARAM definitions in headers generate all defaults
3. **Automatic Validation**: Build-time parameter extraction via `make GenerateParams`
4. **Zero Runtime Overhead**: ConfigManager resolves to direct member variable access
5. **Educational Clarity**: Easy to trace parameter flow from CLI to library call

### Documentation Resources

For complete parameter flow examples and detailed architecture diagrams, see:
- **`ENERGY_SYSTEM_OVERVIEW.md`**: Complete parameter journey with step-by-step code examples
- **`qm_methods/QM_ARCHITECTURE.md`**: ConfigManager integration details and extension guidelines
- **`energy_modules_migration_guide.md`**: Phase-by-phase migration documentation with all 3 phases

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

### ✅ Completed (January-October 2025)
- **Polymorphic Architecture**: Full implementation with all method wrappers
- **MethodFactory**: Priority-based resolution with hierarchical fallbacks
- **API Compatibility**: All existing EnergyCalculator usage preserved
- **Universal Verbosity**: Complete integration across all computational methods
- **Thread Safety**: Enhanced concurrency support maintained
- **Error Handling**: Comprehensive CurcumaLogger-based error reporting
- **ConfigManager Integration**: Complete type-safe parameter system (Phases 3A-3C)
  - 240 parameters with PARAM macro definitions
  - All 8 QM/FF interfaces accept ConfigManager
  - Delegating constructors for backward compatibility
  - End-to-end parameter flow from CLI to external libraries

### 🟡 Remaining TODOs
- **DFT-D3/D4 UpdateParameters**: Migrate to ConfigManager (low priority, seldom used)
- **Native GFN-FF Parameters**: Complete parameter validation and real force field values

### 🔧 Future Enhancements
- **Additional Methods**: Easy integration of new computational methods
- **Performance Optimization**: Method-specific performance tuning
- **Extended Hierarchies**: More complex fallback chains for specialized methods

---

*This directory represents the future of computational method integration in Curcuma - providing educational clarity through polymorphism while maintaining high performance and scientific accuracy.*