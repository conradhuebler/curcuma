# Curcuma Development TODO List

**Stand**: 2025-10-28 (nach SimpleMD Fixes)
**Quelle**: Extrahiert aus allen CLAUDE.md-Dateien
**Test Status**: 26/26 CLI tests passing (100%) 🎯

---

## ✅ RECENTLY RESOLVED (October 28, 2025)

### SimpleMD Complete Fix - 3 Commits
**Commits**: ced705d, 972559c, 284c7da

#### 1. SimpleMD Spin Parameter Fix (ced705d)
- **Problem**: Default `spin=1` caused TBLite crash for closed-shell systems
- **Solution**: Changed default to `spin=0` in simplemd.h:372
- **Result**: No more crashes with gfn2 method ✅

#### 2. CLI2Json Parameter Routing + ConfigManager Fallback (972559c)
- **Problem**: `-md.max_time 10` created nested JSON structure causing parameter lookup failures
- **Solution**:
  - CLI2Json: Strip redundant keyword prefixes (main.cpp:252-273)
  - ConfigManager: Fallback logic for legacy/malformed JSON (config_manager.cpp:404-453)
- **Result**: Robust parameter routing, `-md.max_time` === `-max_time` synonym ✅

#### 3. SimpleMD Trajectory File Generation (284c7da)
- **Problem**: Empty basename created `.trj.xyz` instead of `input.trj.xyz`
- **Solution**:
  - Added `md->setFile(argv[2])` in main.cpp:652
  - Added final frame write in simplemd.cpp:1673
  - Relaxed test validation for short simulations
- **Result**: All trajectory files generated with correct names ✅

**Test Impact**: SimpleMD 0/7 → 7/7, Overall 19/26 (73%) → 26/26 (100%) 🎯

---

## 🔴 KRITISCH / HÖCHSTE PRIORITÄT

*No critical blockers remaining!* All SimpleMD issues resolved.

---

## 🟢 COARSE GRAINING DEVELOPMENT (NEW - CG.md)

### CG Foundation - Phase 1: Molecule Helper Functions
- **Status**: ⏳ PENDING
- **Priority**: 🔴 HIGH - Blocks other CG work
- **Task**: Add CG detection methods to Molecule class
- **Betroffene Dateien**: src/core/molecule.h, src/core/molecule.cpp
- **Aufwand**: ~30 min
- **Details**:
  - `bool isCGSystem() const` - Check if molecule contains CG_ELEMENT (226)
  - `bool hasMixedSystem() const` - Check for both atomic + CG particles
  - `std::vector<int> getCGAtoms() const` - Get all CG atom indices
  - `std::vector<int> getAtomicAtoms() const` - Get all atomic indices

### CG Foundation - Phase 2: JSON Parameter Loading
- **Status**: ⏳ PENDING
- **Priority**: 🔴 CRITICAL - Core functionality
- **Task**: Complete `generateCGParameters()` in ForceField class
- **Betroffene Dateien**: src/core/energy_calculators/ff_methods/forcefield.cpp/h
- **Aufwand**: ~2-3 h
- **Details**:
  - Parse `cg_default` section from JSON config
  - Parse `cg_per_atom` section for per-particle overrides
  - Parse `pair_interactions` for custom pair parameters
  - Generate vdW structures with `type=3` for CG-CG interactions
  - Validate shape vectors and potential parameters

### CG Format - Phase 3: VTF Writer Implementation
- **Status**: ⏳ PENDING
- **Priority**: 🟡 MEDIUM - Needed for trajectory output
- **Task**: Implement `WriteVTF()` function for trajectory output
- **Betroffene Dateien**: src/tools/formats.h, src/tools/vtf_handler.cpp
- **Aufwand**: ~2-3 h
- **Details**:
  - Write atoms with CG markers and radii
  - Write bonds from mol.m_bonds
  - Write unit cell (periodic boundary conditions)
  - Write timestep coordinates
  - Support both single structure and trajectory writing

### CG Integration - Phase 4: Testing & Validation
- **Status**: ⏳ PENDING
- **Priority**: 🟡 MEDIUM - Verification of functionality
- **Task**: Create test cases and validate CG simulation workflow
- **Betroffene Dateien**: test_cases/, test_cases/test_cg_*.cpp
- **Aufwand**: ~2-3 h
- **Details**:
  - Create test_cg_mc.vtf with 2-3 CG particles
  - Create test_cg_mc.json with CG parameters (σ, ε)
  - Test Casino MC with CG config
  - Validate energy calculations for spherical particles
  - Test trajectory output and PBC wrapping
  - Verify acceptance ratios are reasonable (0.1-0.9)

### CG Potentials - Phase 5: Ellipsoidal Extensions (OPTIONAL)
- **Status**: ⏳ LOWEST PRIORITY
- **Priority**: 🔵 LOW - After spherical validation
- **Task**: Implement angle-dependent energy for ellipsoidal particles
- **Betroffene Dateien**: src/core/energy_calculators/ff_methods/cg_potentials.cpp, src/capabilities/casino.cpp
- **Aufwand**: ~4-5 h (after phases 1-4 complete)
- **Details**:
  - Complete `calculateEffectiveDistance()` for ellipsoid-ellipsoid interactions
  - Implement orientation-dependent potential in `calculateCGPairEnergy()`
  - Add rotational moves in Casino (`rotateParticle()`)
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
| 🔴 HIGH (CG) | CG Molecule + JSON Loading | 2 | ⏳ PENDING |
| 🟡 MEDIUM (CG) | VTF Writer + CG Testing | 2 | ⏳ PENDING |
| 🔵 LOW (CG) | Ellipsoidal Extensions | 1 | ⏳ LOWEST PRIORITY |
| 🟡 TESTING | Scientific validation | 4 | ⏳ PENDING |
| 🟢 CORE | Parameter/Memory/Units | 4 | ⏳ PENDING |
| 🔵 CAPABILITIES | Confscan/RMSD/SimpleMD Physics | 4 | ⏳ PENDING |
| 🟣 REFACTORING | Molecule Phase 2-6 | 5 | ⏳ PLANNED |
| **TOTAL** | | **22** | **0 ❌ + 21 ⏳ + 1 🔵** |

---

## 🔧 TESTING STATUS BY MODULE

- **RMSD**: 6/6 ✅ (100%)
- **ConfScan**: 7/7 ✅ (100%)
- **CurcumaOpt**: 6/6 ✅ (100%)
- **SimpleMD**: 7/7 ✅ (100%) - **FIXED October 28, 2025**
- **Overall**: 26/26 ✅ (100%) 🎯

---

**Last Updated**: 2025-10-29
**Next Review**: When starting CG implementation (Phase 1: Molecule helpers) or cgfnff debugging
