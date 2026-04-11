# Curcuma Development TODO List

**Stand**: 2026-03-13
**Quelle**: Extrahiert aus allen CLAUDE.md-Dateien
**Test Status**: 26/26 CLI tests passing (100%) ✅

---

## 🔴 KRITISCH / HÖCHSTE PRIORITÄT

*No critical blockers remaining!* All SimpleMD issues resolved.

---

## 🟢 COARSE GRAINING DEVELOPMENT

### ✅ CG Integration - Phases 1-5 (COMPLETE - Oct/Nov 2025)
- Phase 1: Molecule helper functions (isCGSystem, hasMixedSystem, getCGAtoms)
- Phase 2: JSON parameter loading in forcefield.cpp
- Phase 3: VTF reader/writer in formats.h
- Phase 4: Testing & validation (unit tests + CLI test)
- Phase 5: SimpleMD CG integration (PBC, 10x timestep scaling, VTF output, orientational infrastructure)

### 🔵 CG Potentials - Phase 6: Ellipsoidal Extensions (OPTIONAL - LOWEST PRIORITY)
- **Status**: 🟡 PREPARED
- **Priority**: 🔵 LOW - After SimpleMD integration
- **Task**: Implement angle-dependent energy for ellipsoidal particles
- **Betroffene Dateien**: src/core/energy_calculators/ff_methods/cg_potentials.cpp, src/capabilities/casino.cpp
- **Aufwand**: ~4-5 h (after Phase 5 complete)
- **Current Status**: All infrastructure prepared (rotation matrices, ellipsoid detection, fallback energy)
- **Remaining**:
  - Complete `calculateEffectiveDistance()` for ellipsoid-ellipsoid interactions
  - Implement orientation-dependent potential in `calculateCGPairEnergy()`
  - Activate rotational moves in Casino
  - Test ellipsoidal shape calculations

---

## 🟡 TESTING & VALIDATION (test_cases/CLAUDE.md:296-303)

### Scientific Validation Enhancement
- **Status**: ⏳ PENDING
- **Task**: Erweitere wissenschaftliche Validierung für alle Tests (RMSD tolerances, energy convergence)
- **Dateien**: test_cases/cli/test_utils.sh, test_cases/cli/*/run_test.sh
- **Note**: All tests now passing (26/26), ready for enhanced validation

### Expected Failure Pattern für invalid_method Tests
- **Status**: ⏳ PENDING
- **Task**: Implementiere Pattern für Tests, die bewusst fehlschlagen sollen
- **Tests betroffen**:
  - test_cases/cli/curcumaopt/03_invalid_method/run_test.sh
  - test_cases/cli/rmsd/03_invalid_method/run_test.sh
  - test_cases/cli/confscan/03_invalid_rmsd_method/run_test.sh
- **Dateien**: test_cases/cli/test_utils.sh

### Performance Benchmarks & Regression Tests
- **Status**: ⏳ PLANNED
- **Task**: Füge Performance-Benchmarks hinzu für Regressions-Detection
- **Dateien**: test_cases/CMakeLists.txt, test_cases/cli/test_utils.sh

### test_molecule.cpp Extension
- **Status**: ⏳ PENDING
- **Task**: Erweitere für geplantes Molecule SOA/AOS Refactoring (Phase 2-6)
- **Dateien**: test_cases/test_molecule.cpp
- **Abhängigkeit**: Molecule Refactoring Phases müssen geplant sein

---

## 🟢 CORE DEVELOPMENT (src/core/)

### cgfnff Parameter Generation Bug
- **Status**: ❌ OPEN
- **Problem**: Parameter generation creates null JSON values
- **Betroffene Dateien**: src/core/energy_calculators/qm_methods/gfnff.cpp
- **Verweis**: CLAUDE.md - Known Issues, src/core/CLAUDE.md:104
- **Konsequenz**: Native GFN-FF nicht einsatzbereit

### Missing Real GFN-FF Parameters
- **Status**: ❌ OPEN
- **Problem**: Placeholders statt echter physikalischer Parameter
- **Betroffene Dateien**: src/core/energy_calculators/qm_methods/gfnff.cpp
- **Verweis**: src/core/CLAUDE.md:105
- **Abhängigkeit**: Benötigt theoretische Implementierung oder externe Daten

### Unit System Migration (CODATA-2018)
- **Status**: ⏳ IN PROGRESS
- **Task**: Replace hardcoded constants mit `CurcumaUnit` namespace functions
- **Betroffene Dateien**: Multiple legacy files mit hardcoded constants
- **Verweis**: src/core/CLAUDE.md:113, CLAUDE.md:296
- **Gewinn**: Centralized, documented, CODATA-2018 compliant constants

### Memory Optimization for Large Systems (>1000 atoms)
- **Status**: ⏳ PLANNED
- **Task**: Optimize Molecule data structure and distance matrix caching
- **Betroffene Dateien**: src/core/molecule.cpp/h, src/core/energy_calculators/ff_methods/forcefield.cpp
- **Verweis**: src/core/CLAUDE.md:106, CLAUDE.md - Performance Notes
- **Performance Impact**: Critical for large molecular systems

---

## 🔵 CAPABILITIES & ANALYSIS (src/capabilities/)

### ConfScan Verbosity Enhancement
- **Status**: ⏳ PENDING
- **Problem**: Accept/Reject messages not visible at default verbosity level
- **Task**: Adjust CurcumaLogger calls to be visible at level ≥1
- **Betroffene Dateien**: src/capabilities/confscan.cpp
- **Verweis**: src/capabilities/CLAUDE.md:84

### SimpleMD Wall Potential Physics
- **Status**: ⏳ PENDING (LOW PRIORITY)
- **Task**: Verify wall potential physics (boundary logic, force calculations)
- **Betroffene Dateien**: src/capabilities/simplemd.cpp/h
- **Verweis**: src/capabilities/CLAUDE.md:85
- **Note**: Tests now passing (7/7) - functionality works, physics validation pending

### RMSD Strategy Pattern - Phase 3
- **Status**: ⏳ PENDING
- **Task**: Complete Strategy pattern refactoring for RMSD module
- **Betroffene Dateien**: src/capabilities/rmsd.cpp/h
- **Verweis**: src/capabilities/CLAUDE.md:80
- **Phases**: Phase 1-2 ✅, Phase 3 ⏳

### Enhanced Conformational Search Algorithms
- **Status**: ⏳ PLANNED
- **Betroffene Dateien**: src/capabilities/confsearch.cpp/h
- **Verweis**: src/capabilities/CLAUDE.md:76

### Improved Trajectory Analysis Tools
- **Status**: ⏳ PLANNED
- **Betroffene Dateien**: src/capabilities/rmsdtraj.cpp/h
- **Verweis**: src/capabilities/CLAUDE.md:77

---

## 🟣 MOLECULE REFACTORING - Test-Driven Breaking Changes

**Status**: ✅ Phase 1 DONE, ⏳ Phase 2-6 PENDING
**Test-Driven Approach**: src/core/test_molecule.cpp (15 test categories)
**Critical Requirement**: All existing functionality must remain API-compatible

### Phase 2: XYZ Comment Parser Unification
- **Status**: ⏳ PENDING
- **Task**: Eliminate 10 duplicate XYZ parser functions
- **Betroffene Dateien**: src/core/molecule.cpp/h
- **Critical Constraint**: Production comment formats must NOT break (ORCA, XTB, simple energy)
- **Reference**: docs/XYZ_COMMENT_FORMATS.md
- **Verweis**: CLAUDE.md - Planned Development

### Phase 3: Granular Cache System
- **Status**: ⏳ PLANNED
- **Task**: Replace single `m_dirty` flag with fine-grained cache invalidation
- **Betroffene Dateien**: src/core/molecule.cpp/h
- **Performance Impact**: Selective recalculation instead of full cache invalidation

### Phase 4: Fragment System O(1) Lookups
- **Status**: ⏳ PLANNED
- **Task**: Replace std::map with optimized lookup structure
- **Betroffene Dateien**: src/core/molecule.cpp/h

### Phase 5: Type-Safe ElementType Enum
- **Status**: ⏳ PLANNED
- **Task**: Replace int element indices with ElementType enum
- **Betroffene Dateien**: src/core/molecule.cpp/h, src/core/units.h

### Phase 6: Unified Atom Structure with Zero-Copy Geometry
- **Status**: ⏳ PLANNED
- **Task**: Implement SOA/AOS hybrid design for geometry access
- **Betroffene Dateien**: src/core/molecule.cpp/h

---

## 📊 PRIORITY MATRIX

| Priority | Component | Count | Status |
|----------|-----------|-------|--------|
| 🔴 KRITISCH | None | 0 | ✅ ALL RESOLVED |
| 🟢 DONE (CG Phase 5) | SimpleMD CG Integration | 1 | ✅ COMPLETE |
| 🟢 DONE (CG Phases 1-4) | CG Core + VTF + Testing | 4 | ✅ COMPLETE |
| 🔵 LOW (CG Phase 6) | Ellipsoidal Extensions | 1 | 🟡 PREPARED |
| 🟡 TESTING | Scientific validation | 4 | ⏳ PENDING |
| 🟢 CORE | Parameter/Memory/Units | 4 | ⏳ PENDING |
| 🔵 CAPABILITIES | Confscan/RMSD/SimpleMD Physics | 4 | ⏳ PENDING |
| 🟣 REFACTORING | Molecule Phase 2-6 | 5 | ⏳ PLANNED |
| **TOTAL** | | **23** | **5 ✅ + 16 ⏳ + 2 🟡** |

---

## 🔧 TESTING STATUS BY MODULE

- **RMSD**: 6/6 ✅ (100%)
- **ConfScan**: 7/7 ✅ (100%)
- **CurcumaOpt**: 6/6 ✅ (100%)
- **SimpleMD**: 7/7 ✅ (100%) - **FIXED October 28, 2025**
- **Overall**: 26/26 ✅ (100%) 🎯

---

---

## 🔴 BUILD SYSTEM - CONDITIONAL COMPILATION FIXES (November 2025)

### Status: PARTIALLY FIXED - Build 2 Working, Others Need Deeper Fixes

**Session Date**: November 2025
**Task**: Test all 5 build configurations with different CMake options (USE_D3, USE_D4, USE_TBLITE, USE_ULYSSES, USE_XTB)

### Build Test Results

| Build | Config | Status | Notes |
|-------|--------|--------|-------|
| **Build 1** | Minimal (UFF, EHT only) | ❌ FAILED | `s-dftd3.h` not found - gfnff.cpp includes D3 unconditionally |
| **Build 2** | Standard (TBLite, Ulysses, D3) | ✅ **SUCCESS** | ✅ Verified working: UFF & GFN2 methods tested |
| **Build 3** | Full QM (+ XTB) | ❌ FAILED | Same D3 header chain issue |
| **Build 4** | D4 Dispersion | ❌ FAILED | D3 dependency blocks build |
| **Build 5** | TBLite only | ❌ FAILED | D3/D4 guard chain |

### Root Cause Analysis

**Problem**: Cascading `#ifdef` guards only protect headers, not implementations:
```
gfnff.cpp (always compiles)
  → #include "forcefield.h"
    → #include "forcefieldthread.h"
      → #ifdef USE_D3 #include "dftd3interface.h" #endif
        → #include "s-dftd3.h"  (UNGUARDED INCLUDE!)
```

When `USE_D3=OFF`, compiler skips the `#ifdef` but still tries to compile the file, causing `s-dftd3.h` not found error.

### Fixes Applied ✅

1. **`dftd3interface.h`** - Wrapped entire header with `#ifdef USE_D3...#endif`
2. **`dftd4interface.h`** - Wrapped entire header with `#ifdef USE_D4...#endif`
3. **`forcefieldthread.h`** - Wrapped `D3Thread` class definition with `#ifdef USE_D3...#endif`
4. **`forcefieldthread.cpp`** - Wrapped D3Thread implementation with `#ifdef USE_D3...#endif`

### Remaining Work (To Fix Builds 1, 3, 4, 5)

**Critical files still needing protection**:
- [ ] `gfnff.cpp` - Always compiled, needs conditional compilation or restructuring
- [ ] `forcefield.cpp` - D3Thread instantiation (line ~470) needs `#ifdef USE_D3` guard
- [ ] `gfnff.h` - Consider lazy-loading or factory pattern for D3/D4 dependencies
- [ ] Add H4Thread guards (same pattern as D3Thread)

### Build 2 (Standard) - PRODUCTION READY ✅

**Verified Methods:**
```bash
./curcuma -sp water.xyz -method uff    # ✅ Works
./curcuma -sp water.xyz -method gfn2   # ✅ Works (via TBLite)
```

**Available Methods in Build 2:**
- UFF (universal force field) - native
- EHT (extended Hückel theory) - native
- GFN2 (tight-binding DFT via TBLite) - recommended
- GFN1 (TBLite or Ulysses fallback)
- iPEA1 (TBLite)
- PM6, PM3, AM1, MNDO (Ulysses semi-empirical)

### Recommendations

**Short term**: Protect gfnff.cpp and forcefield.cpp with conditional compilation
**Medium term**: Refactor include chain - move D3/D4/H4 threads to separate file
**Long term**: Use CMake-level validation + CI/CD pipeline for all build configurations

### Test Infrastructure

Created `/home/conrad/src/curcuma/build_test/` with 5 isolated builds for regression testing.

See `docs/BUILD_SYSTEM_ROADMAP.md` for detailed implementation plan.

---

**Last Updated**: 2025-11-08 (Build System Testing)
**Next Review**: After implementing remaining Build 1/3/4/5 fixes or when starting CG Phase 5
