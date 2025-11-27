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

## Verbosity System Standards

### Universal Verbosity Levels (All QM/MM Methods)

#### Level 0: Silent Mode
- **ABSOLUTELY NO OUTPUT** - Critical for optimization and MD
- Only internal calculations, zero console output
- Exception: Critical errors that terminate calculation

#### Level 1: Minimal Results
- Final energy and convergence status only
- Brief method identification
- Essential warnings (convergence failures)

#### Level 2: Properties and Analysis  
- **SCF/Iteration Progress**: Energy per cycle with iteration count
- **Orbital Properties**: HOMO, LUMO, band gap (formatted tables)
- **Molecular Properties**: Charges, bond orders, dipole moment
- **Method Parameters**: Key settings and convergence criteria

#### Level 3: Complete Analysis
- **Full Orbital Lists**: All energies in structured tables
- **Coefficients**: MO contributions, population analysis  
- **Debug Details**: SCF parameters, convergence mechanisms
- **Method Internals**: Basis sets, grid details, algorithm specifics

### Implementation Guidelines
```cpp
// Standard pattern for all QM methods
if (CurcumaLogger::get_verbosity() >= 1) {
    CurcumaLogger::energy_abs(final_energy, "SCF Energy");
}
if (CurcumaLogger::get_verbosity() >= 2) {
    // Orbital properties, molecular analysis
}
if (CurcumaLogger::get_verbosity() >= 3) {
    // Full coefficient matrices, debug output
}
```

## Variable Section

### Current Development Status ✅
- **✅ Universal Verbosity**: **COMPLETED** - All QM methods fully integrated with CurcumaLogger
- **✅ ConfigManager Integration**: **COMPLETED (Oktober 2025)** - All QM interfaces accept ConfigManager (Phases 3A-3C)
  - 240 parameters with PARAM macro definitions across 12 modules
  - All 8 QM/FF interfaces migrated to ConfigManager constructors
  - Type-safe parameter access replacing JSON-based configuration
  - End-to-end parameter flow from CLI to external libraries
- **✅ EHT Implementation**: Fully functional with orbital analysis and verbosity control
- **✅ XTB/TBLite Interfaces**: Native library verbosity synchronized with CurcumaLogger
- **✅ Ulysses Interface**: Complete CurcumaLogger integration with SCF progress
- **🔧 Native GFN-FF (cgfnff)**: Architecture complete, parameter debugging in progress

### Verbosity Integration Status ✅
- **✅ EHT**: Complete integration with `printOrbitalAnalysisVerbose()` for Level 2
- **✅ XTBInterface**: XTB native verbosity synchronized (XTB_VERBOSITY_MUTED/MINIMAL/FULL)
- **✅ TBLiteInterface**: TBLite context verbosity synchronized with CurcumaLogger levels
- **✅ UlyssesInterface**: Full SCF progress, orbital analysis, and molecular properties
- **✅ All Methods**: Silent mode (Level 0) for optimization/MD, debug mode (Level 3) available

### Active Issues

#### ConfigManager Migration
- **DFT-D3/D4 UpdateParameters**: 🟡 LOW PRIORITY - `UpdateParameters()` methods still accept JSON instead of ConfigManager
  - Affects: `dftd3interface.h/cpp`, `dftd4interface.cpp`, `dispersion_method.cpp`
  - Impact: Seldom used, low priority for migration
  - Solution: Add ConfigManager overloads with delegating JSON constructors (Pattern established in Phase 3B)

#### Native GFN-FF (cgfnff)
- **Parameter Generation**: JSON serialization creates null values for some parameters
- **Placeholder Parameters**: Missing real GFN-FF force field parameters from literature
- **Validation**: Incomplete parameter consistency checks with external GFN-FF reference

**GFN-FF Implementation Architecture** (for complete checklist see `../ff_methods/CLAUDE.md`):

- ✅ **Parameter Generation** (GFNFF class in gfnff.cpp)
  - CN-dependent radii, EEQ charges, topology detection, hybridization
  - Methods: generateTopologyAwareBonds(), generateGFNFFDispersionPairs(), etc.

- ✅ **Term Calculation** (ForceFieldThread in ../ff_methods/forcefieldthread.cpp)
  - Multi-threaded energy/gradient calculations
  - Methods: CalculateGFNFFBondContribution(), CalculateGFNFFDispersionContribution(), etc.

- ✅ **All 7 Terms Implemented**: Bond, Angle, Torsion, Inversion, Dispersion, Repulsion, Coulomb

**GFN-FF Components** (topology-aware parameter generation):
- ✅ CN-Berechnung (Coordination Numbers) - D3 method in `gfnff.cpp:calculateCoordinationNumbers()`
- ✅ CN-dependent radii - `r0_gfnff[86]` and `cnfak_gfnff[86]` in `gfnff.cpp`
- ✅ Row-dependent electronegativity - `p_enpoly[6][2]` for 6 periods
- ✅ Hybridization detection - Topology-based in `gfnff.cpp:determineHybridization()`
- ✅ Hybridization bond-strength matrix - `bsmat[4][4]` for mixed hybridizations
- ✅ EEQ charges - `gfnff.cpp:calculateEEQCharges()` with angewChem2020 parameters
- ✅ EEQ charge-dependent corrections (fqq) - Sigmoid function for bond strength
- ✅ Ring strain (ringf) - `findSmallestRings()`, fringbo=0.020 for 3-6 membered rings
- ✅ XH bond corrections (fxh) - Element-specific: BH(+10%), NH(+6%), OH(-7%)
- ✅ CN-dependent heavy atom (fcn) - Coordination-based bond weakening for Z>10

#### Other Issues
- **Memory Optimization**: Needed for large basis sets (>1000 atoms)
- **Ulysses D3H4X/D3H+ Corrections**: Corrections calculated internally but not extractable via getter methods - energies remain identical with/without corrections

### Recent Major Achievements (January-October 2025)
- **🎯 Universal Verbosity System**: All QM methods support consistent 4-level output control
- **🔧 ConfigManager Integration**: Complete type-safe parameter system (Phases 3A-3C)
  - 240 parameters with PARAM macro definitions
  - All 8 QM/FF interfaces accept ConfigManager constructors
  - Type-safe parameter access: `config.get<int>("accuracy")`
  - End-to-end parameter flow from CLI to external libraries
- **🔧 Native Library Integration**: XTB and TBLite verbosity controlled by CurcumaLogger
- **📊 Scientific Output**: HOMO/LUMO analysis, orbital properties, molecular analysis at Level 2
- **🚀 Performance**: Zero overhead silent mode for iterative calculations
- **🏗️ Enhanced Error Handling**: All methods use CurcumaLogger for consistent error reporting
- **✅ Complete Ulysses Integration**: All 27 semi-empirical methods (9 base × 3 correction modes) functional in Curcuma
- **🔧 MethodFactory Enhancement**: Fixed AM1/PM3 method recognition and universal method calculation support

### Performance Optimizations
- Threading support in matrix operations
- Efficient integral calculation algorithms
- Memory-optimized basis set handling
- **Silent Mode**: Zero-overhead Level 0 for iterative calculations

---

*This documentation covers all quantum mechanical methods and computational infrastructure. See QM_ARCHITECTURE.md for detailed technical specifications.*