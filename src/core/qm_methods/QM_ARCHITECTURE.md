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

### JSON-Based Configuration
**Location**: `interface/abstract_interface.h`

```cpp
static json QMInterfaceJson{
    { "threads", 1 },
    { "charge", 0 },
    { "muli", 1 },
    { "solver", "eigen" },
    { "method", "none" },
    { "basis", "sto-3g" },
    { "verbose", false },
    { "gradient", false },
    { "maxiter", 100 },
    { "scfconv", 1e-6 },
    // ... additional parameters
};
```

**Benefits**:
- **Consistent defaults**: All methods inherit common settings
- **Runtime configuration**: Parameters can be modified without recompilation
- **Method-specific overrides**: Each method can customize its configuration

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

3. **Add method-specific configuration**:
   - Extend JSON configuration
   - Add parameter validation
   - Document new parameters

4. **Integrate with build system**:
   - Update CMakeLists.txt
   - Add unit tests
   - Update documentation

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