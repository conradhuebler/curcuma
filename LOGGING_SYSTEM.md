# Curcuma Universal Logging System

## Overview

All computational modules in Curcuma follow a standardized verbosity system using `CurcumaLogger`. This ensures consistent output behavior across the entire codebase and supports different use cases from automated scripting to educational exploration.

## Universal Verbosity Levels

### **Level 0: Silent Mode**
- **Purpose**: Complete silence for automated workflows
- **Output**: Zero output (critical for optimization/MD where calculations are intermediate steps)
- **Use Cases**: 
  - Optimization loops where intermediate Hessian/energy calculations should be silent
  - Batch processing and scripting
  - Production calculations in automated workflows
- **Logger Functions**: Only `CurcumaLogger::error()` (always visible)

### **Level 1: Standard Results (DEFAULT)**
- **Purpose**: Essential results for typical users
- **Output**: 
  - Final results (energies, frequencies, etc.)
  - Method selection and basic settings
  - Convergence status
  - Critical warnings and information
- **Logger Functions**: `error()`, `warn()`, `success()`
- **Use Cases**: Normal interactive use, standard calculations

### **Level 2: Scientific Analysis (Educational)**
- **Purpose**: Educational and research use - show the science
- **Output**:
  - All Level 1 content
  - Scientific context and physical interpretation
  - Algorithm parameters with their physical meaning
  - Intermediate scientific results (HOMO/LUMO, energy components)
  - Method comparisons and analysis
  - Performance timing for optimization
- **Logger Functions**: `error()`, `warn()`, `success()`, `info()`, `param()`, `citation()`
- **Use Cases**: Learning, research, method development

### **Level 3: Full Debug**
- **Purpose**: Algorithm development and troubleshooting
- **Output**:
  - All Level 2 content
  - Detailed algorithmic progress
  - Intermediate matrices and vectors (if reasonable size)
  - Thread assignments and synchronization
  - Memory allocation details
  - Numerical stability checks
- **Logger Functions**: All logger functions available
- **Use Cases**: Debugging, development, algorithm validation

### **Level 4+: Deep Debug (Future Extension)**
- **Purpose**: Extreme detail for specialized debugging
- **Potential Output**:
  - Complete matrices even for large systems
  - Individual iteration details
  - Memory usage tracking
  - Profiling information
  - External library verbose output

## CurcumaLogger Functions by Verbosity Level

| Function | Level | Color | Purpose |
|----------|-------|-------|---------|
| `error()` | Always | Red | Critical errors (always visible) |
| `warn()` | ≥1 | Orange | Important warnings |
| `success()` | ≥1 | Green | Successful completion messages |
| `info()` | ≥2 | Default | Educational/scientific information |
| `param()` | ≥2 | Blue | Structured parameter output |
| `energy_abs()` | ≥2 | Default | Energy values with appropriate units |
| `citation()` | ≥2 | Green | Scientific references and citations |

## Implementation Guidelines

### For All Computational Modules

1. **Default Verbosity**: Always use Level 1 as default (`m_verbosity = 1`)
2. **Silent Mode**: Level 0 must produce zero output except errors
3. **Educational Content**: Level 2 includes physical interpretation and scientific context
4. **Method Information**: Always show method selection at Level 1
5. **Timing**: Performance information at Level 2+
6. **Thread Safety**: All logging must be thread-safe

### Module-Specific Content

#### **Energy Calculators (QM/MM Methods)**
- **Level 1**: Method, final energy, convergence
- **Level 2**: HOMO/LUMO, energy components, citations
- **Level 3**: SCF iterations, orbital details, matrix dimensions

#### **Optimization Methods**
- **Level 1**: Method, final energy, convergence criteria
- **Level 2**: Step sizes, gradient norms, constraint satisfaction
- **Level 3**: Individual optimization steps, line search details

#### **Hessian Calculations**
- **Level 1**: Method (numerical/seminumerical), threading, final frequencies
- **Level 2**: Finite difference parameters, mass-weighting, normal mode projection
- **Level 3**: Matrix symmetry, eigenvalue analysis, thread assignments

#### **Molecular Dynamics**
- **Level 1**: Method, final statistics, completion
- **Level 2**: Thermodynamic properties, trajectory analysis
- **Level 3**: Individual time steps, constraint violations

#### **Force Fields**
- **Level 1**: Method, parameter loading status, final energies
- **Level 2**: Energy decomposition, parameter sources, caching statistics
- **Level 3**: Individual force calculations, parameter generation details

## External Library Integration

When wrapping external libraries (XTB, TBLite, etc.), synchronize their verbosity:
- **Curcuma Level 0** → External library silent mode
- **Curcuma Level 1** → External library minimal output
- **Curcuma Level 2+** → External library detailed output

## JSON Parameter Integration

All modules should accept verbosity through JSON configuration:
```json
{
  "verbosity": 1,  // Global default
  "module_name": {
    "verbosity": 2  // Module-specific override
  }
}
```

## Educational Philosophy

The logging system serves Curcuma's educational mission:
- **Level 0**: Production use
- **Level 1**: Practical results
- **Level 2**: Learning and understanding
- **Level 3**: Deep algorithmic insight

Each level adds educational value while maintaining practical usability.

---

*This system ensures consistent behavior across all Curcuma modules while supporting both practical computation and educational exploration.*