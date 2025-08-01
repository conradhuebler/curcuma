# CLAUDE.md - Quantum Methods Directory

## Overview

The qm_methods directory contains the quantum mechanical method implementations and interfaces for Curcuma. This module provides a flexible, extensible framework for various quantum chemical methods including semi-empirical, tight-binding, and ab initio approaches.

## Structure

```
qm_methods/
├── interface/
│   ├── abstract_interface.h  # Base QMInterface class
│   ├── ulysses.cpp/h         # Ulysses semi-empirical interface
├── qm_driver.cpp/h           # Base driver for matrix-based QM methods
├── eht.cpp/h                 # Extended Hückel Theory implementation
├── eht_parameters.cpp/h      # EHT parameter database
├── gfnff.cpp/h               # Native GFN-FF implementation (cgfnff)
├── gfnff_advanced.cpp/h      # Advanced GFN-FF features
├── xtbinterface.cpp/h        # XTB method interface
├── tbliteinterface.cpp/h     # TBLite method interface
├── dftd3interface.cpp/h      # DFT-D3 dispersion corrections
├── dftd4interface.cpp/h      # DFT-D4 dispersion corrections
├── STOIntegrals.hpp          # Slater-type orbital integrals
├── GTOIntegrals.hpp          # Gaussian-type orbital integrals
├── ParallelEigenSolver.hpp   # Parallel matrix diagonalization
├── basissetparser.hpp        # Basis set file parsing
└── QM_ARCHITECTURE.md        # Detailed architecture documentation
```

## Architecture Overview

### Interface Layer (`QMInterface`)
Unified polymorphic interface for all quantum methods:
```cpp
class QMInterface {
    virtual bool InitialiseMolecule() = 0;
    virtual double Calculation(bool gradient, bool verbose) = 0;
    virtual Vector Charges() const = 0;
    virtual Vector BondOrders() const = 0;
    virtual Geometry Gradient() const = 0;
};
```

### Driver Layer (`QMDriver`)
Base class for matrix-based quantum methods providing:
- Common matrix storage (Hamiltonian, overlap, MO coefficients)
- Threading support and parallel computation
- Template method pattern for customizable calculation steps

## Method Implementations

### Native Methods
- **Extended Hückel Theory (EHT)**: Complete semi-empirical implementation
- **Native GFN-FF (cgfnff)**: Curcuma's own GFN-FF implementation (WORK IN PROGRESS)

### External Interfaces
- **XTB Interface**: Extended tight-binding methods (GFN-FF, GFN1, GFN2)
- **TBLite Interface**: Tight-binding DFT methods (GFN1, GFN2, iPEA1)
- **Ulysses Interface**: Various semi-empirical methods (PM3, AM1, MNDO)
- **DFT-D3/D4**: Dispersion correction interfaces

### Integral Support
- **STO Integrals**: Analytical Slater-type orbital overlap integrals
- **GTO Integrals**: Primitive Gaussian overlap integrals
- **Parallel Solver**: Optimized eigenvalue solver with threading

## EnergyCalculator Integration

Methods are routed via `SwitchMethod()` in energycalculator.cpp:
```cpp
case 9: GFNFF (cgfnff)      // Native GFN-FF - WORK IN PROGRESS
case 6: EHT                 // Extended Hückel Theory
case 4: DFT-D3              // Dispersion corrections
case 3: Ulysses             // Semi-empirical methods
case 2: XTB                 // Extended tight-binding
case 1: TBLite              // Tight-binding DFT
```

## Configuration System

JSON-based configuration with common defaults:
- Threading control
- SCF convergence parameters
- Method-specific settings
- Debugging and verbosity options

## Instructions Block

**PRESERVED - DO NOT EDIT BY CLAUDE**

*Quantum method development priorities, theoretical enhancements, and performance optimization goals to be defined by operator/programmer*

## Variable Section

### Current Development Status
- **Native GFN-FF (cgfnff)**: Architecture complete, JSON parameter generation debugging needed
- **EHT Implementation**: Fully functional with orbital analysis
- **All interfaces**: Working and integrated with EnergyCalculator

### Active Issues
- cgfnff JSON parameter generation creates null values causing crashes
- Missing real GFN-FF parameters (currently using placeholders)
- Memory optimization needed for large basis sets

### Recent Developments
- Complete QMInterface standardization across all methods
- Enhanced parallel eigenvalue solver implementation
- Improved basis set parsing and validation
- Better error handling in method initialization

### Performance Optimizations
- Threading support in matrix operations
- Efficient integral calculation algorithms
- Memory-optimized basis set handling

---

*This documentation covers all quantum mechanical methods and computational infrastructure. See QM_ARCHITECTURE.md for detailed technical specifications.*