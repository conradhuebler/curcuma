# Curcuma Development TODO List

**Stand**: 2025-10-27 (nach Dokumentations-Audit)
**Quelle**: Extrahiert aus allen CLAUDE.md-Dateien
**Test Status**: 19/26 CLI tests passing (73%)

---

## 🔴 KRITISCH / HÖCHSTE PRIORITÄT

### SimpleMD JSON Null-Crash (URGENT)
- **Status**: ❌ NICHT GELÖST (Dokumentation widerspricht Realität!)
- **Problem**: Parameter-Routing in SimpleMD broken, crasht mit `json.exception.type_error.306`
- **Evidenz**: `test_cases/cli/simplemd/01_nve_short_run/stderr.log` zeigt Crash
- **Betroffene Dateien**: `src/capabilities/simplemd.cpp/h`, `src/main.cpp`
- **Verweis**: test_cases/CLAUDE.md:290, test_cases/CLAUDE.md:149
- **Aktion**: Parameter-Registry Migration für SimpleMD (analog zu gelöstem CurcumaOpt-Bug vom 26.10.2025)

### Dokumentations-Inkonsistenzen korrigieren
- **Status**: ⏳ PENDING
- **Test Status korrekt**: 19/26 (73%) ✅, aber SimpleMD falsch dokumentiert
- **Korrekturen erforderlich**:
  1. test_cases/CLAUDE.md:290 - "✅ FIXED SimpleMD JSON" → "⏳ PARTIALLY FIXED (curcumaopt works, simplemd crashes)"
  2. test_cases/CLAUDE.md:149 - "0/7 (no crashes, output file issue)" → "0/7 (crashes with JSON null error)"
- **Betroffene Dateien**: test_cases/CLAUDE.md

---

## 🟡 TESTING & VALIDATION (test_cases/CLAUDE.md:296-303)

### Scientific Validation Enhancement
- **Status**: ⏳ PENDING
- **Task**: Erweitere wissenschaftliche Validierung für alle Tests (RMSD tolerances, energy convergence)
- **Dateien**: test_cases/cli/test_utils.sh, test_cases/cli/*/run_test.sh
- **Abhängigkeit**: SimpleMD crash muss gelöst sein für vollständige Testabdeckung

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

### SimpleMD Wall Potential Issues
- **Status**: ⏳ PENDING
- **Problems**:
  1. Boundary logic accuracy
  2. Force calculation correctness
  3. Output file generation (0/7 tests generate trajectory)
- **Betroffene Dateien**: src/capabilities/simplemd.cpp/h
- **Verweis**: src/capabilities/CLAUDE.md:85
- **Test Coverage**: test_cases/cli/simplemd/* (0/7 passing)
- **Abhängigkeit**: SimpleMD JSON crash muss gelöst sein

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
| 🔴 KRITISCH | SimpleMD crash + Docs | 2 | ❌ BLOCKING |
| 🟡 TESTING | Scientific validation | 4 | ⏳ PENDING |
| 🟢 CORE | Parameter/Memory/Units | 4 | ⏳ PENDING |
| 🔵 CAPABILITIES | Confscan/SimpleMD/RMSD | 5 | ⏳ PENDING |
| 🟣 REFACTORING | Molecule Phase 2-6 | 5 | ⏳ PLANNED |
| **TOTAL** | | **20** | **2 ❌ + 13 ⏳ + 5 ⏳** |

---

## 🔧 TESTING STATUS BY MODULE

- **RMSD**: 6/6 ✅ (100%)
- **ConfScan**: 7/7 ✅ (100%)
- **CurcumaOpt**: 6/6 ✅ (100%)
- **SimpleMD**: 0/7 ❌ (JSON crash - BLOCKING)
- **Overall**: 19/26 ✅ (73%)

---

**Last Updated**: 2025-10-27
**Next Review**: After SimpleMD JSON crash resolved
