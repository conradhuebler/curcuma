# CLAUDE.md - Capabilities Directory

## Overview

The capabilities directory contains high-level molecular modeling applications and computational chemistry tasks. These modules use the core computational engines to perform complex multi-step calculations and analyses.

## Structure

```
capabilities/
├── confscan.cpp/h         # Conformational scanning along reaction coordinates
├── confsearch.cpp/h       # Systematic conformational searching
├── curcumaopt.cpp/h       # Geometry optimization algorithms
├── simplemd.cpp/h         # Molecular dynamics simulation
├── rmsd.cpp/h            # Structure comparison and alignment
├── rmsdtraj.cpp/h        # Trajectory RMSD analysis
├── hessian.cpp/h         # Second derivative calculations
├── persistentdiagram.cpp/h # Topological data analysis
├── optimiser/            # Optimization algorithms
│   ├── lbfgs.cpp/h       # LBFGS optimization
│   └── LevMar*.h         # Levenberg-Marquardt variants
└── c_code/               # C interface code (Hungarian algorithm)
```

## Key Capabilities

### Conformational Analysis
- **ConfScan**: Systematic scanning of conformational space along defined coordinates
- **ConfSearch**: Automated conformational searching with energy filtering
- **RMSD Analysis**: Structure comparison, alignment, and trajectory analysis

### Optimization
- **CurcumaOpt**: Geometry optimization using various algorithms
- **LBFGS**: Limited-memory Broyden-Fletcher-Goldfarb-Shanno optimizer
- **Constrained optimization**: Fix atomic positions by setting gradient = 0

### Molecular Dynamics
- **SimpleMD**: Basic molecular dynamics with various thermostats
- **NEB Docking**: Nudged elastic band for transition state searches
- **Trajectory Analysis**: Analysis tools for MD trajectories

### Advanced Analysis
- **Hessian**: Second derivative calculations for normal modes
- **Persistent Diagrams**: Topological data analysis for molecular structures
- **Enhanced TDA (dMatrix replacement)**: Complete topological data analysis with TDAEngine
- **Pairmapper**: Advanced structure matching algorithms

## Development Guidelines

### Interface Design
- All capabilities use the core EnergyCalculator for energy/gradient evaluations
- Consistent parameter handling through JSON-based configuration
- Progress reporting and result persistence for long calculations

### Performance Considerations
- Multi-threading support where applicable
- Memory-efficient handling of large trajectory data
- Automatic checkpointing for resumable calculations

## Instructions Block

**PRESERVED - DO NOT EDIT BY CLAUDE**

*Future development tasks and visions to be defined by operator/programmer*

## Variable Section

### Parameter System Migration Status (October 2025)
**Status**: ✅ Infrastructure complete, ConfigManager layer implemented, ongoing capability migration
**Documentation**: [docs/PARAMETER_SYSTEM.md](../../docs/PARAMETER_SYSTEM.md) + [docs/PARAMETER_MIGRATION_GUIDE.md](../../docs/PARAMETER_MIGRATION_GUIDE.md)

**ConfigManager Implementation** (October 2025):
- ✅ **Core Implementation**: `src/core/config_manager.{h,cpp}` - Type-safe parameter access layer
- ✅ **Proof-of-Concept**: analysis.cpp - 37 Json2KeyWord calls eliminated, hierarchical dot notation working
- 🎯 **Goal**: Replace 395 Json2KeyWord calls across 17 files with elegant `config.get<T>()` API

**Migration Pattern**:
```cpp
// Old approach (legacy):
std::string format = Json2KeyWord<std::string>(m_config, "output_format");
bool flag = Json2KeyWord<bool>(m_config, "topological_save_image");

// New approach (ConfigManager):
ConfigManager m_config("module", controller);
std::string format = m_config.get<std::string>("output_format");
bool flag = m_config.get<bool>("topological.save_image");  // hierarchical!
```

**Migrated Modules**:
- ✅ **analysis** (25 parameters, 37 Json2KeyWord calls) - Complete with ConfigManager, production-ready reference
- ✅ **rmsd** (~15 parameters, 33 Json2KeyWord calls) - Complete with ConfigManager, enum-based method selection
- ✅ **simplemd** (48 parameters, 82 Json2KeyWord calls) - LARGEST migration, ~350 LOC boilerplate eliminated
- ✅ **confscan** (41 parameters, 59 Json2KeyWord calls) - Multi-Module architecture, imports RMSD parameters, sLX parsing simplified

**Pending Migration** (Priority Order):
- ⏳ **casino** (~14 parameters, 36 Json2KeyWord calls) - Monte Carlo simulation, simple structure
- ⏳ **opt** (~15 parameters, estimated 50+ calls) - Critical optimization module
- ⏳ **hessian** (~10 parameters, 14 calls) - Second derivatives

**Migration Requirements for New Capabilities**:
1. **Parameter Definitions**: Include `src/core/parameter_macros.h` in header
2. **PARAM Block**: Add `BEGIN_PARAMETER_DEFINITION(module)` block with PARAM entries
3. **ConfigManager Integration** (Modern Approach):
   - Include `src/core/config_manager.h` in header
   - Member variable: `ConfigManager m_config;` (replaces `json m_config`)
   - Constructor: `m_config("module", controller)` (no MergeJson needed!)
   - Parameter access: `m_config.get<T>("parameter")` or `m_config.get<T>("nested.parameter")`
   - Legacy JSON support: `m_config_legacy = m_config.exportConfig()` if needed for external APIs
4. Use **snake_case** for all parameter names (MANDATORY)
5. Remove static JSON configuration from .cpp files
6. Build: `make GenerateParams` before compilation
7. Verify: No validation warnings, help output correct
8. Test: Default params + custom CLI args + aliases work

### Current Development
- Enhanced conformational search algorithms
- Improved trajectory analysis tools
- Better integration with quantum chemical methods
- **✅ COMPLETED: Parameter Registry System (Phase 1 & 2)** - Production-ready infrastructure
  - Auto-generated help, type validation, alias resolution, JSON export/import
  - Build-time parameter extraction from PARAM macros
  - Reference implementation: analysis module (25 parameters)
- **✅ COMPLETED: dMatrix Integration** - Legacy -dMatrix functionality fully integrated into analysis.cpp
  - TDAEngine class provides all original features plus enhancements
  - Migration guide available in DMATRIX_MIGRATION.md
  - Research-grade topological data analysis maintained
- **Unit system migration**: Capabilities should migrate to centralized `CurcumaUnit` namespace
  - Replace hardcoded constants (2625.5, 627.5, etc.) with `CurcumaUnit::Energy::*` functions
  - Use centralized constants from `src/core/units.h` for consistency
- **RMSD Code Restructuring**: ✅ Phase 1 & 2 completed
  - **Phase 1**: Extracted CostMatrixCalculator and AssignmentSolver utilities
  - **Phase 2**: Refactored LoadControlJson() from monolithic (200+ lines) into 6 thematic methods
  - **Phase 3**: TODO - Strategy pattern for 7 alignment methods (incr, template, hybrid, subspace, inertia, molalign, dtemplate)

### Known Issues
- Memory usage optimization needed for large systems
- Enhanced error handling for failed optimizations
- **ConfScan verbosity issue**: Accept/Reject messages not showing at default verbosity level - need to investigate default verbosity settings and ensure user feedback is always visible

### ✅ COMPLETED - Universal Restart Validation System (October 2025)
**ROBUSTNESS ENHANCEMENT FOR ALL CAPABILITIES**

#### Implementation
**CurcumaMethod Base Class** - Universal validation infrastructure:
- `RestartValidationResult` struct for validation results
- `computeRestartChecksum()` - Checksumme über angegebene Felder (Standard C++ `<functional>`)
- `isValidDoubleString()` - Erkennt malformierte Zahlen wie "-na" statt "-nan"
- `validateRestartData()` - Zentrale Validierung mit:
  - Format-Version Check (mit Warning für Inkompatibilitäten)
  - Required Fields Validation
  - Checksummen-Validierung (erkennt Disk-Korruption, manuelle Edits)
  - String-Validierung für pipe-separated doubles

**Automatic Checksumming** in `TriggerWriteRestart()`:
- Fügt automatisch `format_version: "1.0"` hinzu
- Fügt `timestamp: <unix_time>` hinzu
- Berechnet automatisch Checksumme über alle pipe-separated Felder (geometry, velocities, etc.)
- Speichert `checksum_fields: [...]` (welche Felder gehasht wurden)

**SimpleMD Integration** (Proof-of-Concept):
```cpp
auto validation = validateRestartData(state,
    {"method", "geometry", "velocities"},  // required fields
    {"geometry", "velocities", "xi", "Q"}); // validate for malformed doubles

if (!validation.valid) {
    cerr << "[ERROR] " << validation.error_message << endl;
    return false;  // Start fresh simulation
}
```

#### Benefits
✅ **Keine Dependencies** - 100% Standard C++ (`<functional>`, `<sstream>`)
✅ **Automatisch** - Alle Capabilities profitieren via `TriggerWriteRestart()`
✅ **Früherkennung** - Validierung VOR Parsing verhindert `stod()` Crashes
✅ **Checksumme** - Erkennt Hardware-Fehler, Disk-Korruption, manuelle Edits
✅ **Zukunftssicher** - `format_version` erlaubt zukünftige Migrations-Logik
✅ **Wiederverwendbar** - Jede Capability kann `validateRestartData()` aufrufen

#### Testing
- ✅ Corrupted restart files with "-na" truncated NaN: **Detected & clean error**
- ✅ New restart files automatically include checksums and metadata
- ✅ Case-insensitive parameter matching via ConfigManager

#### Future Work
- Extend validation to other capabilities (opt, confscan, etc.)
- Add checksum verification UI feedback
- Implement automatic backup of corrupted files

---

### TODO - SimpleMD Wall Potential Fixes (January 2025)
**CRITICAL CORRECTNESS ISSUES IDENTIFIED:**

1. **Fix boundary logic error** (lines 814-829 in simplemd.cpp)
   - Current: `if (m_wall_x_min - m_wall_x_max < 1)` always triggers auto-config
   - Should: `if (m_wall_x_max - m_wall_x_min <= 0)` to detect invalid bounds
   - Impact: User-defined boundaries are ignored, always auto-configured

2. **Fix rectangular harmonic force calculation** (lines 2430-2434)
   - Current: Uses `std::abs(x - x_min)` which gives wrong force direction
   - Should: `(x - x_min)` without abs() for correct restoring force
   - Impact: Atoms pushed in wrong direction, simulation instability

3. **Fix log-Fermi force denominator** (lines 2318-2320, 2357-2359)
   - Current: Uses `(1 - exp_expr)` in denominator
   - Should: `(1 + exp_expr)` - correct derivative of log(1 + e^x)
   - Impact: Wrong force magnitudes, potential numerical instability

4. **Standardize parameter meanings**
   - Current: `wall_temp` means force constant k OR temperature kbT
   - Current: `wall_beta` unused in harmonic, critical for log-Fermi
   - Should: Consistent parameter interpretation across wall types
   - Add parameter validation and clear documentation

5. **Add numerical stability safeguards**
   - Check for distance = 0 in spherical calculations
   - Prevent exponential overflow in log-Fermi potentials
   - Add parameter bounds checking

6. **Add comprehensive wall potential tests**
   - Unit tests for force calculations vs analytical derivatives
   - Boundary condition tests (atoms at/beyond walls)
   - Parameter validation tests
   - Comparison with reference implementations

**Priority**: HIGH - These are correctness bugs that cause wrong physics
**Status**: IDENTIFIED - Analysis complete, fixes needed

---

*This documentation covers all molecular modeling capabilities and applications*