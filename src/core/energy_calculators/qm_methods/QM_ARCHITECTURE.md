# Quantum Methods Architecture Documentation

## Overview

This document describes the architecture and design patterns used in the quantum methods (QM) module of the Curcuma computational chemistry package. The module provides a flexible, extensible framework for implementing various quantum chemical methods.

## Core Design Principles

### 1. **Interface Segregation and Polymorphism**
- All QM methods inherit from a common `QMInterface` base class
- Polymorphic method dispatch allows uniform treatment of different QM methods
- Clear separation between interface definition and implementation

### 2. **Strategy Pattern Implementation**
- Different QM methods (EHT, XTB, TBLite, etc.) implement the same interface
- Methods can be swapped at runtime without changing client code
- Enables easy comparison and benchmarking of different methods

### 3. **Template Method Pattern**
- `QMDriver` provides common functionality for matrix-based QM methods
- Derived classes implement specific parts (overlap, Hamiltonian construction)
- Reduces code duplication while maintaining flexibility

## Architecture Layers

### Layer 1: Abstract Interface (`QMInterface`)
**File**: `interface/abstract_interface.h`

```cpp
class QMInterface {
public:
    // Molecule initialization with multiple overloads
    virtual bool InitialiseMolecule(const Mol& molecule);
    virtual bool InitialiseMolecule(const Mol* molecule);
    virtual bool InitialiseMolecule(const int* attyp, const double* coord, ...);
    
    // Geometry updates
    virtual bool UpdateMolecule(const Mol& molecule);
    virtual bool UpdateMolecule(const double* coord);
    
    // Core calculation method
    virtual double Calculation(bool gradient = false, bool verbose = false) = 0;
    
    // Property accessors
    virtual Vector Charges() const;
    virtual Vector Dipole() const;
    virtual Vector BondOrders() const;
    virtual Geometry Gradient() const;
};
```

**Key Features**:
- **Multiple initialization methods**: Supports Mol objects, raw arrays, pointers
- **Flexible geometry updates**: Coordinate updates without full reinitialization
- **Property interface**: Standardized access to molecular properties
- **Default JSON configuration**: Common settings for all QM methods

### Layer 2: Driver Base Class (`QMDriver`)
**File**: `qm_driver.h`

```cpp
class QMDriver : public QMInterface {
protected:
    // Common data members
    Mol m_mol;                    // Molecule object
    Matrix m_H, m_S;              // Hamiltonian and overlap matrices
    Matrix m_mo;                  // Molecular orbital coefficients
    Vector m_energies;            // Orbital energies
    int m_num_electrons;          // Total number of electrons
    int m_threads;                // Thread count for parallel operations

    // Pure virtual methods for derived classes
    virtual Matrix MakeOverlap(Basisset& basisset) = 0;
    virtual Matrix MakeH(const Matrix& S, const Basisset& basisset) = 0;
};
```

**Key Features**:
- **Matrix storage**: Manages overlap, Hamiltonian, and MO coefficient matrices
- **Threading support**: Configurable parallelization
- **Template methods**: Common workflow with customizable steps

### Layer 3: Concrete Implementations

#### Extended Hückel Theory (EHT)
**Files**: `eht.h`, `eht.cpp`, `eht_parameters.h`, `eht_parameters.cpp`

```cpp
class EHT : public QMDriver {
private:
    double m_K;  // Wolfsberg-Helmholz constant
    
    // Implementation-specific methods
    Basisset MakeBasis();
    Matrix MakeOverlap(Basisset& basisset) override;
    Matrix MakeH(const Matrix& S, const Basisset& basisset) override;
    
    // Analysis methods
    double getHOMOLUMOGap() const;
    void printOrbitalAnalysis(int num_orbitals_around_gap = 5) const;
};
```

**Features**:
- **Parameterized approach**: Uses literature VSIP values and Slater exponents
- **Modular parameters**: Separate parameter management in `EHT::ParameterDatabase`
- **Orbital analysis**: HOMO-LUMO gap calculation and orbital visualization

#### XTB Interface
**Files**: `xtbinterface.h`, `xtbinterface.cpp`

```cpp
class XTBInterface : public QMInterface {
private:
    // XTB-specific data
    std::vector<double> m_charges;
    std::vector<double> m_dipole;
    std::vector<double> m_bond_orders;
    
    // XTB method configuration
    int m_method;  // GFN-FF, GFN0, GFN1, GFN2
    double m_accuracy;
    int m_max_iterations;
};
```

**Features**:
- **External library integration**: Interface to XTB quantum chemistry package
- **Multiple method support**: GFN-FF, GFN0, GFN1, GFN2 methods
- **Property calculation**: Charges, dipoles, bond orders

## Mathematical Infrastructure

### Integral Calculation Modules

#### Slater-Type Orbitals (STO)
**File**: `STOIntegrals.hpp`

```cpp
namespace STO {
    enum OrbitalType { S, PX, PY, PZ, DXX, DYY, DZZ, DXY, DXZ, DYZ };
    
    struct Orbital {
        OrbitalType type;
        double x, y, z;      // Position
        double zeta;         // Slater exponent
        double VSIP;         // Valence State Ionization Potential
        int atom;            // Atom index
    };
    
    // Integral calculation functions
    double calculateOverlap(const Orbital& a, const Orbital& b);
    double calculateNormalization(const Orbital& orbital);
}
```

**Key Features**:
- **Complete STO implementation**: s, p, d orbital support
- **Analytical integrals**: Exact overlap calculations
- **Proper normalization**: Ensures orthogonality conditions

#### Gaussian-Type Orbitals (GTO)
**File**: `GTOIntegrals.hpp`

```cpp
namespace GTO {
    // Primitive Gaussian overlap integrals
    double overlapIntegral(double alpha_a, double alpha_b, 
                          const Vector3& Ra, const Vector3& Rb);
    
    // Gaussian product theorem
    double gaussianProduct(double alpha_a, double alpha_b, 
                          const Vector3& Ra, const Vector3& Rb);
}
```

### Parallel Computing Support
**File**: `ParallelEigenSolver.hpp`

```cpp
class ParallelEigenSolver {
public:
    // Constructor with convergence parameters
    ParallelEigenSolver(int max_iter, int block_size, double tolerance, bool use_divide_conquer);
    
    // Solve generalized eigenvalue problem HC = SCE
    bool solve(const Matrix& S, const Matrix& H, 
               Vector& eigenvalues, Matrix& eigenvectors, 
               int num_threads, bool debug = false);
    
    // Configuration methods
    void setThreadCount(int threads);
    void setBlockSize(int block_size);
    void setTolerance(double tolerance);
};
```

**Features**:
- **Divide-and-conquer algorithms**: Efficient for large matrices
- **Thread-based parallelization**: Configurable thread count
- **Block diagonalization**: Memory-efficient for very large systems

## Configuration Management

### ConfigManager & Parameter Registry System (2025)

**Modern Architecture**: All QM methods now use the **ConfigManager + Parameter Registry** system for type-safe, validated parameter access.

#### **Parameter Definition (in Method Headers)**
**Example from `xtbinterface.h`**:

```cpp
#include "src/core/parameter_macros.h"
#include "src/core/config_manager.h"

// Define parameters using PARAM macros
BEGIN_PARAMETER_DEFINITION(xtb)
    PARAM(accuracy, Int, 2, "Accuracy level for XTB calculations (0-4).", "SCF", {"xtb_ac"})
    PARAM(max_iterations, Int, 100, "Maximum number of SCF iterations.", "SCF", {"SCFmaxiter"})
    PARAM(electronic_temperature, Double, 300.0, "Electronic temperature in Kelvin.", "SCF", {"Tele"})
    PARAM(spin, Double, 0.0, "Total spin of the system (0.0 = singlet).", "Molecular", {})
END_PARAMETER_DEFINITION

class XTBInterface : public QMInterface {
public:
    XTBInterface(const ConfigManager& config);  // Modern constructor
    // ...
private:
    mutable ConfigManager m_config;
};
```

#### **Type-Safe Parameter Access (in Implementation)**
**Example from `xtbinterface.cpp`**:

```cpp
XTBInterface::XTBInterface(const ConfigManager& config)
    : m_config(config)
{
    // Type-safe parameter access with defaults
    m_accuracy = m_config.get<int>("accuracy", 2);
    m_SCFmaxiter = m_config.get<int>("max_iterations", 100);
    m_Tele = m_config.get<double>("electronic_temperature", 300.0);
    m_spin = m_config.get<double>("spin", 0.0);
}
```

#### **Build-Time Parameter Extraction**
Parameters are automatically extracted from header files during build:
```bash
make GenerateParams  # Scans PARAM blocks → generates parameter_registry.h
```

**Benefits**:
- **Type safety**: Compile-time type checking via `get<T>()`
- **Auto-generated help**: Parameters documented in code generate user help
- **Single source of truth**: No manual JSON/help synchronization
- **Validation**: Build-time duplicate detection and type validation
- **Aliases**: Backward compatibility via old parameter names
- **Case-insensitive**: `-MaxIter` and `-maxiter` both work

**See also**: `docs/PARAMETER_SYSTEM.md`, `docs/PARAMETER_MIGRATION_GUIDE.md`

## Data Flow and Lifecycle

### 1. **Initialization Phase**
```
Client Code → QMInterface::InitialiseMolecule() → Method-specific setup
    ↓
Molecule validation → Parameter loading → Memory allocation
```

### 2. **Calculation Phase**
```
QMInterface::Calculation() → Basis set construction → Matrix setup
    ↓
Overlap matrix → Hamiltonian matrix → Eigenvalue solution
    ↓
Property calculation → Results storage
```

### 3. **Property Access**
```
QMInterface::Charges() → Method-specific property getter
QMInterface::Gradient() → Gradient vector return
QMInterface::OrbitalEnergies() → Energy eigenvalues
```

## Error Handling Strategy

### 1. **Validation at Interface Level**
- Molecule structure validation
- Parameter range checking
- Memory allocation verification

### 2. **Exception-Safe Implementation**
- RAII principles for resource management
- Exception propagation with meaningful messages
- Graceful degradation for non-critical failures

### 3. **Method-Specific Error Handling**
- Library integration error checking (XTB, TBLite)
- Convergence failure detection
- Numerical stability monitoring

## Extension Guidelines

### Adding a New QM Method

1. **Inherit from appropriate base class**:
   - `QMInterface` for completely custom methods
   - `QMDriver` for matrix-based methods

2. **Implement required virtual methods**:
   ```cpp
   class NewMethod : public QMDriver {
   public:
       double Calculation(bool gradient, bool verbose) override;
   private:
       Matrix MakeOverlap(Basisset& basisset) override;
       Matrix MakeH(const Matrix& S, const Basisset& basisset) override;
   };
   ```

3. **Define method-specific parameters** (ConfigManager system):
   ```cpp
   // In newmethod.h
   #include "src/core/parameter_macros.h"
   #include "src/core/config_manager.h"

   BEGIN_PARAMETER_DEFINITION(newmethod)
       PARAM(max_iterations, Int, 100, "SCF iteration limit.", "SCF", {})
       PARAM(convergence, Double, 1e-6, "Energy convergence threshold.", "SCF", {})
       // Add all method-specific parameters here
   END_PARAMETER_DEFINITION

   class NewMethod : public QMDriver {
   public:
       NewMethod(const ConfigManager& config);  // Modern constructor
   private:
       mutable ConfigManager m_config;
   };
   ```

4. **Implement ConfigManager-based constructor**:
   ```cpp
   // In newmethod.cpp
   NewMethod::NewMethod(const ConfigManager& config)
       : m_config(config)
   {
       // Type-safe parameter access
       m_max_iter = m_config.get<int>("max_iterations", 100);
       m_conv_thr = m_config.get<double>("convergence", 1e-6);
   }
   ```

5. **Integrate with build system**:
   - Run `make GenerateParams` to extract parameter definitions
   - Update CMakeLists.txt
   - Add unit tests
   - Update documentation

**See also**: `docs/PARAMETER_MIGRATION_GUIDE.md` for complete migration workflow

### Best Practices

1. **Memory Management**:
   - Use RAII for all resources
   - Prefer smart pointers for dynamic allocation
   - Clear large matrices when not needed

2. **Performance Optimization**:
   - Use Eigen's optimized matrix operations
   - Implement threading where beneficial
   - Profile critical calculation paths

3. **Code Quality**:
   - Follow consistent naming conventions
   - Add comprehensive documentation
   - Write unit tests for all new functionality

## ConfigManager Integration (2025)

### Complete Parameter Flow Through the System

**The ConfigManager system provides end-to-end type-safe parameter management** from user input to low-level QM interfaces.

#### **Data Flow Chain**

```
User CLI/JSON Input
    ↓
Capability (e.g., Opt, ConfScan, SimpleMD)
    ↓  Creates ConfigManager with module-specific defaults
EnergyCalculator(method, ConfigManager&)
    ↓  Delegates to ConfigManager constructor
MethodFactory::create(method, json)
    ↓  Exports JSON for backward compatibility
Method Wrapper (e.g., XTBMethod, TBLiteMethod)
    ↓  Creates ConfigManager("xtb", json)
QM Interface (e.g., XTBInterface)
    ↓  Type-safe parameter access: config.get<int>("accuracy")
External Library (XTB, TBLite, Ulysses)
```

#### **Example: Complete Parameter Journey**

**1. User Input (CLI)**:
```bash
curcuma -opt input.xyz -method gfn2 -max_iterations 200 -accuracy 3
```

**2. Capability Layer** (`curcumaopt.cpp`):
```cpp
ConfigManager opt_config("opt", controller["opt"]);
EnergyCalculator calculator("gfn2", opt_config);
```

**3. EnergyCalculator** (`energycalculator.cpp`):
```cpp
EnergyCalculator::EnergyCalculator(const std::string& method, const ConfigManager& config)
    : m_method_name(method)
{
    m_controller = config.exportConfig();  // For MethodFactory compatibility
    createMethod(method_name, m_controller);
}
```

**4. MethodFactory** (`method_factory.cpp`):
```cpp
std::unique_ptr<ComputationalMethod> MethodFactory::create(
    const std::string& method_name, const json& config)
{
    // Priority resolution: gfn2 → TBLite > Ulysses > XTB
    return std::make_unique<TBLiteMethod>("gfn2", config);
}
```

**5. Method Wrapper** (`tblite_method.cpp`):
```cpp
TBLiteMethod::TBLiteMethod(const std::string& method_name, const json& config)
{
    ConfigManager tblite_config("tblite", config);  // Extract tblite-specific params
    m_tblite = std::make_unique<TBLiteInterface>(tblite_config);
}
```

**6. QM Interface** (`tbliteinterface.cpp`):
```cpp
TBLiteInterface::TBLiteInterface(const ConfigManager& config)
    : m_config(config)
{
    // Type-safe parameter access with defaults from PARAM definitions
    m_acc = m_config.get<int>("accuracy", 1);              // User: 3
    m_SCFmaxiter = m_config.get<int>("max_iterations", 100); // User: 200
    m_damping = m_config.get<double>("damping", 0.4);       // Default: 0.4
}
```

**7. External Library**:
```cpp
tblite_set_calculator_accuracy(m_ctx, m_tblite_calc, m_acc);         // 3
tblite_set_calculator_max_iter(m_ctx, m_tblite_calc, m_SCFmaxiter);  // 200
```

#### **Key Benefits of This Architecture**

1. **Type Safety Throughout**: Every layer uses typed access (`get<int>`, `get<double>`)
2. **Single Source of Truth**: PARAM definitions in headers generate all defaults
3. **Automatic Validation**: Build-time parameter extraction catches errors early
4. **Backward Compatibility**: Old JSON code still works via `exportConfig()`
5. **Educational Clarity**: Easy to trace parameter flow from CLI to library
6. **Zero Runtime Overhead**: ConfigManager resolves to direct member access

#### **Migration Status (October 2025)**

**✅ Completed**:
- ✅ Phase 1: Parameter definitions (240 parameters across 12 modules)
- ✅ Phase 2: Interface constructors (XTB, TBLite, Ulysses, DFT-D3/D4, GFN-FF, ForceFieldGenerator)
- ✅ Phase 3: Method wrappers (8 wrappers updated)
- ✅ Phase 3C: EnergyCalculator integration (delegating constructors)

**🟡 Remaining TODOs**:
- DFT-D3/D4 `UpdateParameters()` method signatures (still use JSON)
- Native GFN-FF (cgfnff) parameter generation completeness

**See Full Migration Details**: `energy_modules_migration_guide.md`

## Future Architecture Considerations

### 1. **Plugin Architecture**
- Dynamic loading of QM method implementations
- Runtime method discovery and registration
- Separation of core from method-specific code

### 2. **GPU Acceleration**
- CUDA/OpenCL integration for matrix operations
- GPU-accelerated integral calculation
- Heterogeneous computing support

### 3. **Distributed Computing**
- MPI support for large-scale calculations
- Cluster-aware job scheduling
- Network-transparent property access

This architecture provides a solid foundation for quantum chemical method implementation while maintaining flexibility for future extensions and optimizations.