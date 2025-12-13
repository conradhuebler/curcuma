# CLAUDE.md - Force Field Methods Directory

## Overview

Force field implementation system with multi-threading support for UFF, QMDFF, and native GFN-FF.

## Architecture

### Core Separation of Concerns

**ForceFieldThread** (`forcefieldthread.cpp/h`):
- **Responsibility**: Energy and gradient calculations for all force field terms
- **Threading**: Multi-threaded via CxxThreadPool
- **Methods**: Calculate{Method}{Term}Contribution() pattern (e.g., CalculateGFNFFBondContribution)

**ForceField** (`forcefield.cpp/h`):
- **Responsibility**: Thread pool management, geometry updates, energy accumulation
- **Role**: Dispatcher that coordinates ForceFieldThread instances

**ForceFieldGenerator** (`forcefieldgenerator.cpp/h`):
- **Responsibility**: UFF parameter generation from atom types
- **Used by**: UFF method only

**GFNFF Class** (`gfnff_method.cpp/h` + `../qm_methods/gfnff.cpp/h`):
- **Responsibility**: GFN-FF parameter generation (topology-aware) + ConfigManager integration
- **Methods**: generateTopologyAwareBonds(), generateGFNFFDispersionPairs(), etc.
- **Output**: JSON parameter sets passed to ForceField

### GFN-FF Implementation Pattern

**Two-Phase Architecture:**

1. **Parameter Generation** (in GFNFF class):
   ```cpp
   // gfnff.cpp with ConfigManager migration
   ConfigManager config("gfnff", parameters);
   json params = config.exportConfig();
   params["bonds"] = generateTopologyAwareBonds(...);
   params["angles"] = generateTopologyAwareAngles(...);
   params["gfnff_dispersions"] = generateGFNFFDispersionPairs();
   // etc.
   m_forcefield->setParameter(params);
   ```

2. **Term Calculation** (in ForceFieldThread with parameter flags):
   ```cpp
   // forcefieldthread.cpp - Phase 2 parameter flag checks
   if (m_dispersion_enabled) {
       CalculateGFNFFDispersionContribution();  // Only if enabled
   }
   if (m_hbond_enabled) {
       CalculateGFNFFHydrogenBondContribution(); // Only if enabled
   }
   ```

### Adding New GFN-FF Terms - Checklist

To add a new GFN-FF energy term (e.g., "CrossTerm"), you MUST modify:

1. **Parameter Structures** (`forcefieldthread.h`):
   - Add new struct (e.g., `GFNFFCrossTerm { int i,j; double param1; }`)
   - Add member vector (e.g., `std::vector<GFNFFCrossTerm> m_gfnff_crossterms;`)

2. **Parameter Generation** (`gfnff.cpp`):
   - Add generation method (e.g., `generateGFNFFCrossTerms()`)
   - Call in `generateGFNFFParameters()` and add to JSON

3. **Term Calculation** (`forcefieldthread.cpp`):
   - Add calculation method (e.g., `CalculateGFNFFCrossTermContribution()`)
   - Call in `ForceFieldThread::execute()` under `if (m_method == 3)` block

4. **Energy Accumulation** (`forcefield.cpp`):
   - Add energy component member variable if needed
   - Collect from threads in `ForceField::Calculate()`

5. **Parameter Setter** (`forcefield.h/.cpp`):
   - Add setter method (e.g., `setGFNFFCrossTerms(const json&)`)
   - Call in `setParameter()` dispatcher

## Current Implementation Status

### ✅ Fully Implemented Terms
- Bond stretching (exponential potential)
- Angle bending (cosine + damping + fqq correction Phase 5A)
- Dihedral torsion (cosine series)
- Inversion (out-of-plane)
- Dispersion (D3/D4 Becke-Johnson damping - Phase 4 Pairwise)
- Repulsion (exponential r^-1.5 - Phase 4 Pairwise)
- Coulomb (EEQ + erf damping - Phase 4 Pairwise)

### ✅ Parameter Management (Phase 2 - December 2025)
- **ConfigManager Integration**: Type-safe parameter access with validation
- **Parameter Flags**: Selective term calculation (dispersion, hbond, repulsion, coulomb enabled/disable)
- **Legacy Code Removal**: CalculateGFNFFvdWContribution deprecated and removed
- **Test Coverage**: Parameter flag combinations test suite with 5 scenarios
  - Dispersion/hbond disabled tests
  - All non-bonded terms disabled test
  - Edge case (atoms at cutoff distance)
  - Metal-specific correction handling (Fe atom)

### 🟡 Lower Priority TODOs
- Phase 5B: Metal-specific fqq correction (2.5x factor)

## Performance

**Multi-threading Benchmarks** (water.xyz, 4 cores):
- 1 thread: 0.320s
- 4 threads: 0.120s
- Speedup: 2.67x ✅

## References

- ForceFieldThread implements formulas from Fortran `gfnff_engrad.F90`
- GFNFF parameter generation follows Spicher/Grimme J. Chem. Theory Comput. 2020
