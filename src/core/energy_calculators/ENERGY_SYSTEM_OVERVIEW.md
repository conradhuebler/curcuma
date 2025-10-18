# Energy System Overview - Complete Architecture Guide

## Introduction

This document provides a **comprehensive overview** of the Curcuma Energy System architecture, bridging the three main documentation sources:

1. **energy_modules_migration_guide.md** - ConfigManager migration process (Phases 1-3)
2. **CLAUDE.md** - Polymorphic architecture and design patterns
3. **QM_ARCHITECTURE.md** - Detailed quantum methods implementation

**Target Audience**: Developers adding new computational methods or understanding the parameter flow from CLI to external libraries.

---

## High-Level System Architecture

```
┌─────────────────────────────────────────────────────────────────────────┐
│                         USER INTERACTION LAYER                          │
│  CLI: curcuma -opt input.xyz -method gfn2 -max_iterations 200          │
│  JSON: { "method": "gfn2", "max_iterations": 200, "accuracy": 3 }      │
└──────────────────────────────────┬──────────────────────────────────────┘
                                   ↓
┌─────────────────────────────────────────────────────────────────────────┐
│                         CAPABILITY LAYER                                 │
│  Modules: CurcumaOpt, ConfScan, SimpleMD, Hessian, etc.                │
│  Role: Task orchestration, workflow management                          │
│  ├─ Create ConfigManager with module-specific defaults                  │
│  └─ Instantiate EnergyCalculator with method + config                   │
└──────────────────────────────────┬──────────────────────────────────────┘
                                   ↓
┌─────────────────────────────────────────────────────────────────────────┐
│                    ENERGYCALCULATOR DISPATCHER                          │
│  File: src/core/energycalculator.cpp                                   │
│  Role: Unified interface for all computational methods                  │
│  ├─ ConfigManager constructor (modern, preferred)                       │
│  ├─ JSON constructor (backward compatible, delegates to ConfigManager)  │
│  └─ Delegates to MethodFactory for method creation                      │
└──────────────────────────────────┬──────────────────────────────────────┘
                                   ↓
┌─────────────────────────────────────────────────────────────────────────┐
│                         METHOD FACTORY                                   │
│  File: src/core/energy_calculators/method_factory.cpp                  │
│  Role: Priority-based method resolution with automatic fallbacks        │
│  ├─ gfn2: TBLite → Ulysses → XTB (hierarchical priority)               │
│  ├─ gfn1: TBLite → XTB → Ulysses                                       │
│  ├─ ipea1: TBLite only                                                  │
│  ├─ uff/qmdff: ForceField methods                                       │
│  └─ eht/cgfnff: Native implementations                                  │
└──────────────────────────────────┬──────────────────────────────────────┘
                                   ↓
┌─────────────────────────────────────────────────────────────────────────┐
│                         METHOD WRAPPERS                                  │
│  Files: src/core/energy_calculators/qm_methods/*_method.cpp            │
│  Role: Polymorphic ComputationalMethod interface implementation         │
│  ├─ XTBMethod, TBLiteMethod, UlyssesMethod (QM)                        │
│  ├─ ForceFieldMethod (MM)                                               │
│  ├─ EHTMethod, GFNFFMethod (native)                                     │
│  └─ DispersionMethod (D3/D4)                                            │
│  Action: Creates ConfigManager("module", json) for interface            │
└──────────────────────────────────┬──────────────────────────────────────┘
                                   ↓
┌─────────────────────────────────────────────────────────────────────────┐
│                        QM/MM INTERFACES                                  │
│  Files: src/core/energy_calculators/qm_methods/*interface.cpp          │
│  Role: Type-safe parameter extraction and library interface             │
│  ├─ XTBInterface, TBLiteInterface, UlyssesInterface                     │
│  ├─ DFTD3Interface, DFTD4Interface, GFNFFInterface                      │
│  └─ ForceField, ForceFieldGenerator                                     │
│  Action: ConfigManager::get<T>("param") → member variables              │
└──────────────────────────────────┬──────────────────────────────────────┘
                                   ↓
┌─────────────────────────────────────────────────────────────────────────┐
│                      EXTERNAL LIBRARIES                                  │
│  XTB (libxtb), TBLite (libtblite), Ulysses, DFT-D3/D4                  │
│  Role: Actual computational chemistry calculations                      │
└─────────────────────────────────────────────────────────────────────────┘
```

---

## Complete Parameter Flow Example

### User Command
```bash
curcuma -opt input.xyz -method gfn2 -max_iterations 200 -accuracy 3 -convergence 1e-7
```

### Step-by-Step Transformation

#### **1. CLI Parsing → JSON Controller**
**File**: `src/main.cpp`

```json
controller = {
  "opt": {
    "method": "gfn2",
    "max_iterations": 200,
    "accuracy": 3,
    "convergence": 1e-7
  }
}
```

#### **2. Capability Instantiation**
**File**: `src/capabilities/curcumaopt.cpp`

```cpp
ConfigManager opt_config("opt", controller["opt"]);
// Merges user input with default parameters from ParameterRegistry

EnergyCalculator calculator("gfn2", opt_config);
```

#### **3. EnergyCalculator Constructor**
**File**: `src/core/energycalculator.cpp`

```cpp
EnergyCalculator::EnergyCalculator(const std::string& method, const ConfigManager& config)
    : m_method_name(method)
{
    // Export ConfigManager to JSON for backward compatibility with MethodFactory
    m_controller = config.exportConfig();

    // Create computational method using factory
    createMethod(method_name, m_controller);
}
```

**Result**: JSON with merged defaults
```json
{
  "method": "gfn2",
  "max_iterations": 200,      // User override
  "accuracy": 3,              // User override
  "convergence": 1e-7,        // User override
  "damping": 0.4,             // Default from ParameterRegistry
  "electronic_temperature": 300.0  // Default
}
```

#### **4. MethodFactory Resolution**
**File**: `src/core/energy_calculators/method_factory.cpp`

```cpp
std::unique_ptr<ComputationalMethod> MethodFactory::create(
    const std::string& method_name, const json& config)
{
    // method_name = "gfn2"
    // Check priority methods: gfn2 → TBLite > Ulysses > XTB

    // Priority 1: Try TBLite
    if (hasTBLite()) {
        return std::make_unique<TBLiteMethod>("gfn2", config);  ✅
    }

    // Priority 2: Try Ulysses (fallback if TBLite unavailable)
    if (hasUlysses()) {
        return std::make_unique<UlyssesMethod>("ugfn2", config);
    }

    // Priority 3: Try XTB (last resort)
    if (hasXTB()) {
        return std::make_unique<XTBMethod>("gfn2", config);
    }
}
```

#### **5. Method Wrapper Creation**
**File**: `src/core/energy_calculators/qm_methods/tblite_method.cpp`

```cpp
TBLiteMethod::TBLiteMethod(const std::string& method_name, const json& config)
    : m_method_name(method_name)
{
    // Create ConfigManager for TBLite-specific parameters
    ConfigManager tblite_config("tblite", config);

    // Instantiate TBLite interface with ConfigManager
    m_tblite = std::make_unique<TBLiteInterface>(tblite_config);

    // Set method type
    m_tblite->setMethod("gfn2");
}
```

#### **6. QM Interface Parameter Extraction**
**File**: `src/core/energy_calculators/qm_methods/tbliteinterface.cpp`

```cpp
TBLiteInterface::TBLiteInterface(const ConfigManager& config)
    : m_config(config)
{
    // Type-safe parameter extraction with defaults from PARAM definitions
    m_acc = m_config.get<int>("accuracy", 1);                      // → 3 (user)
    m_SCFmaxiter = m_config.get<int>("max_iterations", 100);       // → 200 (user)
    m_damping = m_config.get<double>("damping", 0.4);              // → 0.4 (default)
    m_Tele = m_config.get<double>("electronic_temperature", 300.0);// → 300.0 (default)
    m_Tele /= 315775.326864009;  // Convert Kelvin to atomic units

    // Initialize TBLite library context
    m_error = tblite_new_error();
    m_ctx = tblite_new_context();
    m_tblite_res = tblite_new_result();
}
```

#### **7. External Library Configuration**
**File**: `src/core/energy_calculators/qm_methods/tbliteinterface.cpp`

```cpp
double TBLiteInterface::Calculation(bool gradient)
{
    // Create TBLite calculator
    m_tblite_calc = tblite_new_gfn2_calculator(m_ctx, m_tblite_mol);

    // Apply user parameters to TBLite
    tblite_set_calculator_accuracy(m_ctx, m_tblite_calc, m_acc);           // 3
    tblite_set_calculator_max_iter(m_ctx, m_tblite_calc, m_SCFmaxiter);    // 200
    tblite_set_calculator_mixer_damping(m_ctx, m_tblite_calc, m_damping);  // 0.4
    tblite_set_calculator_temperature(m_ctx, m_tblite_calc, m_Tele);       // 300 K → AU

    // Perform SCF calculation
    tblite_get_singlepoint(m_ctx, m_tblite_mol, m_tblite_calc, m_tblite_res);
    tblite_get_result_energy(m_error, m_tblite_res, &energy);

    return energy;
}
```

---

## ConfigManager System Deep Dive

### Parameter Definition (PARAM Macros)

**Location**: Method header files (e.g., `xtbinterface.h`)

```cpp
#include "src/core/parameter_macros.h"
#include "src/core/config_manager.h"

BEGIN_PARAMETER_DEFINITION(xtb)
    PARAM(accuracy, Int, 2,
          "Accuracy level for XTB calculations (0-4).",
          "SCF",
          {"xtb_ac"})  // Alias for backward compatibility

    PARAM(max_iterations, Int, 100,
          "Maximum number of SCF iterations.",
          "SCF",
          {"SCFmaxiter"})

    PARAM(electronic_temperature, Double, 300.0,
          "Electronic temperature in Kelvin for Fermi smearing.",
          "SCF",
          {"Tele"})

    PARAM(spin, Double, 0.0,
          "Total spin of the system (0.0 = singlet).",
          "Molecular",
          {})
END_PARAMETER_DEFINITION
```

### Build-Time Parameter Extraction

```bash
# During build, CMake runs the parameter parser
make GenerateParams

# Parser scans all header files for PARAM blocks
curcuma_param_parser src/core/energy_calculators/qm_methods/*.h \
    → build/generated/parameter_registry.h

# Generated registry contains:
class ParameterRegistry {
    static const std::vector<ParameterDef> m_params;
    // Contains all 240+ parameters from 12 modules
};
```

### Runtime Parameter Access

```cpp
ConfigManager config("xtb", user_json);

// Type-safe access with default fallback
int accuracy = config.get<int>("accuracy", 2);           // Default: 2
int maxiter = config.get<int>("max_iterations", 100);    // Default: 100

// Case-insensitive alias resolution
int iter1 = config.get<int>("max_iterations");  // ✅
int iter2 = config.get<int>("MaxIterations");   // ✅ Case-insensitive
int iter3 = config.get<int>("SCFmaxiter");      // ✅ Alias resolution

// All three return the same value!
```

---

## Method Resolution & Fallback Hierarchy

### Priority-Based Resolution

The MethodFactory implements **hierarchical fallbacks** to ensure robustness:

```
┌─────────────────────────────────────────────┐
│   User requests: "gfn2"                     │
└──────────────────┬──────────────────────────┘
                   │
        ┌──────────┴──────────┐
        │  Priority Method?   │
        └──────────┬──────────┘
                   │ YES
        ┌──────────┴──────────────────┐
        │  Try providers in order:    │
        │  1. TBLite (preferred)      │
        │  2. Ulysses (fallback)      │
        │  3. XTB (last resort)       │
        └──────────┬──────────────────┘
                   │
    ┌──────────────┴───────────────┐
    │   Provider 1: TBLite?        │
    └──────────────┬───────────────┘
                   │
         ┌─────────┴─────────┐
         │  YES (compiled)   │  NO (not compiled)
         │                   │
         ↓                   ↓
    ┌────────────┐    ┌──────────────────┐
    │  SUCCESS   │    │  Try Provider 2  │
    │  Return    │    │    Ulysses?      │
    │  TBLite    │    └─────────┬────────┘
    └────────────┘              │
                     ┌──────────┴──────────┐
                     │  YES      │   NO    │
                     │           │         │
                     ↓           ↓         ↓
               ┌─────────┐  ┌─────────┐  ┌──────┐
               │ SUCCESS │  │Try XTB? │  │ FAIL │
               │ Return  │  │         │  │Throw │
               │ Ulysses │  └────┬────┘  │Error │
               └─────────┘       │        └──────┘
                             ┌───┴───┐
                             │SUCCESS│
                             │ XTB   │
                             └───────┘
```

### Explicit vs Priority Methods

**Priority Methods** (shared, with fallbacks):
- `gfn2`: TBLite → Ulysses → XTB
- `gfn1`: TBLite → XTB → Ulysses
- `ipea1`: TBLite only
- `gfnff`: External GFN-FF → XTB → Native cgfnff

**Explicit Methods** (single provider):
- `xtb-gfn1`, `xtb-gfn2`: XTB library explicitly
- `ugfn2`, `pm6`, `am1`: Ulysses methods explicitly
- `uff`, `uff-d3`, `qmdff`: Force field methods
- `eht`, `cgfnff`: Native implementations
- `d3`, `d4`: Dispersion corrections

---

## Polymorphic Architecture (ComputationalMethod Interface)

### Unified Interface

**File**: `src/core/energy_calculators/computational_method.h`

```cpp
class ComputationalMethod {
public:
    virtual ~ComputationalMethod() = default;

    // Molecule management
    virtual bool setMolecule(const Mol& mol) = 0;
    virtual bool updateGeometry(const Matrix& geometry) = 0;

    // Calculation
    virtual double calculateEnergy(bool gradient = false) = 0;
    virtual Matrix getGradient() const = 0;

    // Properties
    virtual Vector getCharges() const = 0;
    virtual Position getDipole() const = 0;
    virtual Vector getBondOrders() const = 0;
    virtual Vector getOrbitalEnergies() const = 0;
    virtual Vector getOrbitalOccupations() const = 0;

    // Metadata
    virtual std::string getMethodName() const = 0;
    virtual bool hasGradient() const = 0;
    virtual bool isThreadSafe() const = 0;

    // Configuration
    virtual void setParameters(const json& params) = 0;
    virtual json getParameters() const = 0;

    // Error handling
    virtual bool hasError() const = 0;
    virtual std::string getErrorMessage() const = 0;
    virtual void clearError() = 0;
};
```

### Benefits vs Old Switch-Based System

**Old System** (Pre-January 2025):
```cpp
// 200+ line SwitchMethod with duplicated logic
double EnergyCalculator::CalculateEnergy(bool gradient) {
    switch(m_method_id) {
        case 0:  // ForceField
            if (m_forcefield) {
                return m_forcefield->Calculate(gradient);
            }
            break;
        case 1:  // TBLite
            if (m_tblite) {
                return m_tblite->Calculation(gradient);
            }
            break;
        case 2:  // XTB
            if (m_xtb) {
                return m_xtb->Calculation(gradient);
            }
            break;
        // ... 10+ more cases
    }
}
```

**New System** (January 2025):
```cpp
// Single polymorphic call
double EnergyCalculator::CalculateEnergy(bool gradient) {
    if (!m_method) {
        handleMethodError("no method available");
        return 0.0;
    }
    return m_method->calculateEnergy(gradient);  // ✅ Polymorphic dispatch
}
```

**Advantages**:
1. ✅ No giant switch statements
2. ✅ Easy to add new methods (no dispatcher changes)
3. ✅ Consistent error handling
4. ✅ Educational clarity (direct polymorphism)
5. ✅ Reduced code duplication
6. ✅ Type-safe method access

---

## Verbosity System Integration

All computational methods support **4-level verbosity** via CurcumaLogger:

### Verbosity Levels

**Level 0: Silent** (optimization/MD)
```cpp
// ZERO output - critical for iterative calculations
// Only errors that terminate calculation
```

**Level 1: Minimal Results**
```cpp
CurcumaLogger::energy_abs(energy, "TBLite Energy");
// Output: TBLite Energy: -42.123456 Eh
```

**Level 2: Scientific Analysis**
```cpp
CurcumaLogger::param("HOMO", fmt::format("{:.4f} Eh", homo));
CurcumaLogger::param("LUMO", fmt::format("{:.4f} Eh", lumo));
CurcumaLogger::param("HOMO-LUMO_gap", fmt::format("{:.4f} Eh ({:.2f} eV)", gap, gap * 27.211));
// Output:
// HOMO: -0.2841 Eh
// LUMO: -0.1052 Eh
// HOMO-LUMO_gap: 0.1789 Eh (4.87 eV)
```

**Level 3: Complete Debug**
```cpp
for (int i = 0; i < orbital_energies.size(); ++i) {
    CurcumaLogger::param(fmt::format("Orbital_{}", i + 1),
        fmt::format("{:.6f} Eh (occ={:.3f})", orbital_energies(i), occ));
}
// Output: Full orbital listing, timing, algorithm details
```

### Native Library Verbosity Synchronization

**XTB Interface**:
```cpp
if (CurcumaLogger::get_verbosity() >= 3) {
    xtb_setVerbosity(m_env, XTB_VERBOSITY_FULL);
} else if (CurcumaLogger::get_verbosity() >= 2) {
    xtb_setVerbosity(m_env, XTB_VERBOSITY_MINIMAL);
} else {
    xtb_setVerbosity(m_env, XTB_VERBOSITY_MUTED);
}
```

**TBLite Interface**:
```cpp
if (CurcumaLogger::get_verbosity() >= 3) {
    tblite_set_context_verbosity(m_ctx, 3);  // Full debug
} else if (CurcumaLogger::get_verbosity() >= 2) {
    tblite_set_context_verbosity(m_ctx, 1);  // Minimal
} else {
    tblite_set_context_verbosity(m_ctx, 0);  // Silent
}
```

---

## Migration Status & Remaining TODOs

### ✅ Completed (October 2025)

**Phase 1: Parameter Definitions**
- ✅ 240 parameters across 12 modules
- ✅ XTB (4), TBLite (13), DFT-D3 (10), DFT-D4 (9)
- ✅ Ulysses (5), GFN-FF (3), ORCA (10)
- ✅ ForceFieldGenerator (36 parameters)

**Phase 2: Interface Constructors**
- ✅ 8 QM/MM interfaces migrated to ConfigManager constructors
- ✅ All interfaces use `config.get<T>()` type-safe access

**Phase 3: Method Wrappers**
- ✅ 8 method wrappers updated
- ✅ All create ConfigManager for interface instantiation

**Phase 3C: EnergyCalculator Integration**
- ✅ ConfigManager constructor overloads
- ✅ Delegating constructor pattern for backward compatibility
- ✅ All 12 capabilities work unchanged

### 🟡 Remaining TODOs

**High Priority**:
1. **DFT-D3/D4 UpdateParameters Migration**
   - Current: `void UpdateParameters(const json& controller)` ❌
   - Target: `void UpdateParameters(const ConfigManager& config)` ✅
   - Impact: 3 #pragma TODO markers in code
   - Effort: ~30 minutes

2. **Legacy Json2KeyWord Elimination**
   - Remaining: ~13 calls outside energy_calculators/
   - Target: Replace with `config.get<T>()`
   - Impact: Code consistency, reduced technical debt
   - Effort: ~1-2 hours

**Low Priority** (Functional, but incomplete):
3. **Native GFN-FF (cgfnff) Completion**
   - Status: 28 TODOs for advanced features
   - Impact: Full native implementation vs external library
   - Effort: Multi-week project (research + implementation)

4. **Ulysses Output Filtering**
   - Status: 1 TODO for stdout capture
   - Impact: Better verbosity control
   - Effort: ~1 hour

---

## Cross-References

### Related Documentation

1. **energy_modules_migration_guide.md**
   - Detailed migration steps (Phases 1-3C)
   - Before/After code examples
   - Verification procedures

2. **CLAUDE.md** (`src/core/energy_calculators/CLAUDE.md`)
   - Polymorphic architecture overview
   - Method hierarchies and fallbacks
   - Development status

3. **QM_ARCHITECTURE.md** (`qm_methods/QM_ARCHITECTURE.md`)
   - Detailed quantum methods implementation
   - Integral calculation modules
   - Extension guidelines with ConfigManager

4. **PARAMETER_SYSTEM.md** (`docs/PARAMETER_SYSTEM.md`)
   - Parameter Registry architecture
   - PARAM macro reference
   - Build-time extraction details

5. **PARAMETER_MIGRATION_GUIDE.md** (`docs/PARAMETER_MIGRATION_GUIDE.md`)
   - Step-by-step migration workflow
   - Common patterns and pitfalls
   - Validation procedures

### Quick Navigation

**Adding a new QM method?** → See QM_ARCHITECTURE.md "Extension Guidelines"

**Understanding parameter flow?** → See this document "Complete Parameter Flow Example"

**Migrating old code?** → See energy_modules_migration_guide.md

**ConfigManager API reference?** → See PARAMETER_SYSTEM.md

**Understanding method resolution?** → See this document "Method Resolution & Fallback Hierarchy"

---

## Summary

The Curcuma Energy System provides a **modern, type-safe, and educational** architecture for computational chemistry:

✅ **ConfigManager**: End-to-end type safety from CLI to external libraries
✅ **Polymorphism**: Clean, extensible method interface without switch statements
✅ **Priority Fallbacks**: Robust method resolution with hierarchical providers
✅ **Universal Verbosity**: Consistent 4-level output control across all methods
✅ **Build-Time Validation**: Parameter extraction catches errors before runtime
✅ **Backward Compatibility**: Old JSON code still works via delegation
✅ **Educational Clarity**: Easy to trace and understand parameter flow

**For new developers**: Start with this document's "Complete Parameter Flow Example" to understand how a single CLI parameter travels through the entire system.

**For contributors**: Follow the "Extension Guidelines" in QM_ARCHITECTURE.md to add new computational methods using the modern ConfigManager patterns.
