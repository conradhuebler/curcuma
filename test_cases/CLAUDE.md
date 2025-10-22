# CLAUDE.md - Test Cases Directory

## Overview

Das `test_cases/` Verzeichnis enthält alle Tests für Curcuma, organisiert in mehrere Testkategorien:
- **Unit Tests** (C++) - Direktes Testen von Core-Funktionalität
- **CLI Tests** (Bash) - End-to-End Tests der Command-Line Interface
- **Integration Tests** (C++) - Capability-spezifische Integrationstests

## Structure

```
test_cases/
├── cli/                          # ✅ CLI End-to-End Tests (Bash)
│   ├── curcumaopt/              # Optimization CLI tests (6 tests)
│   ├── rmsd/                    # RMSD CLI tests (6 tests)
│   ├── confscan/                # ConfScan CLI tests (7 tests)
│   ├── simplemd/                # SimpleMD CLI tests (7 tests)
│   ├── test_utils.sh            # Shared test utilities
│   ├── CMakeLists.txt           # CTest integration
│   ├── README.md                # CLI test documentation
│   ├── KNOWN_BUGS.md            # Documented bugs affecting tests
│   └── GOLDEN_REFERENCES.md     # Scientific reference values
├── AAAbGal.cpp                  # RMSD integration test (5 methods)
├── reorder/                     # RMSD reordering test
├── rmsd/                        # RMSD unit test
├── confscan.cpp                 # ConfScan integration test (10 scenarios)
├── test_energy_methods.cpp      # Energy method validation suite
├── test_molecule.cpp            # Molecule data structure tests
└── CMakeLists.txt               # Main test configuration
```

## CLI Tests (NEW - October 2025)

### Design Philosophy

**End-to-End Testing**: Tests verwenden curcuma genau wie ein User es verwenden würde - über die Command Line.

**Scientific Validation**: Tests validieren nicht nur Exit-Codes, sondern wissenschaftliche Korrektheit:
- **RMSD**: Numerische Werte mit Toleranz
- **Optimierung**: Energie-Konvergenz
- **ConfScan**: Anzahl akzeptierter Konformere
- **SimpleMD**: Trajektorien-Länge und Energie-Drift

### Test Organization Pattern

Jeder Test folgt dieser Struktur:
```
cli/capability/XX_test_name/
├── run_test.sh        # Haupttest-Skript
├── input.xyz          # Input-Molekül
├── expected_*         # Erwartete Ausgaben (optional)
├── stdout.log         # Curcuma stdout (generiert)
└── stderr.log         # Curcuma stderr (generiert)
```

### Test Script Pattern

```bash
#!/bin/bash
# Test: Descriptive Name
# Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated - Based on testing_plan_*.md

set -e
source "../../test_utils.sh"

run_test() {
    $CURCUMA -capability input.xyz > stdout.log 2> stderr.log
    assert_exit_code $? 0 "Capability should succeed"
}

validate_results() {
    # Scientific validation
    local value=$(extract_value_from_output stdout.log)
    assert_scientific_value "2.87214" "$value" "0.0001" "Scientific metric"
}

main() {
    test_header "Test Name"
    cleanup_before
    run_test && validate_results
    print_test_summary
    [ $TESTS_FAILED -eq 0 ] && exit 0 || exit 1
}

main "$@"
```

### test_utils.sh - Shared Utilities

Zentrale Bibliothek mit wiederverwendbaren Funktionen:

**Test Management**:
- `test_header()` - Formatierte Test-Überschrift
- `print_test_summary()` - Zusammenfassung am Ende
- `cleanup_test_artifacts()` - Clean up generierte Dateien

**Assertions (Exit Codes)**:
- `assert_exit_code()` - Vergleiche Exit-Code
- `assert_curcuma_success()` - Pragmatische Erfolgs-Prüfung (Datei-Existenz)

**Assertions (Datei-Inhalte)**:
- `assert_file_exists()` - Datei existiert
- `assert_string_in_file()` - String in Datei vorhanden
- `assert_file_not_empty()` - Datei nicht leer

**Scientific Validation**:
- `assert_scientific_value()` - Numerischer Vergleich mit Toleranz
- `compare_float()` - Floating-Point Arithmetik (awk-basiert)
- `extract_energy_from_xyz()` - Energie aus XYZ-Kommentar
- `extract_rmsd_from_output()` - RMSD aus Curcuma-Output
- `count_xyz_structures()` - Strukturen in Multi-XYZ Datei

**Golden References**:
- `store_golden_reference()` - Referenz speichern
- `load_golden_reference()` - Referenz laden

### CTest Integration

```cmake
# test_cases/cli/CMakeLists.txt
add_test(
    NAME cli_capability_XX_test_name
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/capability/XX_test_name/run_test.sh
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/capability/XX_test_name
)
set_tests_properties(cli_capability_XX_test_name PROPERTIES TIMEOUT 30)
```

**Ausführung**:
```bash
# Alle CLI Tests
ctest -R "cli_" --output-on-failure

# Spezifische Kategorie
ctest -R "cli_rmsd_" --output-on-failure

# Einzelner Test
ctest -R "cli_rmsd_01" --verbose
```

### Current Status (Stand 2025-10-19)

**Tests implementiert**: 26 Tests
**Tests bestanden**: 13/26 (50%)
**Kategorien**:
- RMSD: 5/6 ✅ (83%)
- ConfScan: 6/7 ✅ (86%)
- SimpleMD: 1/7 ❌ (14% - Optimierungs-Bug)
- curcumaopt: 1/6 ❌ (17% - JSON null-Fehler)

**Bekannte Probleme**:
- **curcumaopt JSON null-Fehler**: Betrifft 11/26 Tests (siehe KNOWN_BUGS.md)
- **SimpleMD**: 6 Tests fehlgeschlagen wegen Optimierungs-Bug
- **Invalid Method Tests**: 3 Tests erwarten Fehler-Codes (noch nicht implementiert)

**Golden References dokumentiert**:
- RMSD: 2.87214 Å (AAA-bGal, 90 Atome)
- Weitere Referenzen in GOLDEN_REFERENCES.md

## Unit Tests (C++)

### test_molecule.cpp
**Status**: ✅ Vollständiger Test-Suite (15 Kategorien)
**Zweck**: Vorbereitung für Molecule Refactoring (SOA/AOS Design)

Tests decken ab:
- XYZ Parser (10 duplicate functions zu vereinheitlichen)
- Fragment System (O(1) lookups geplant)
- Cache System (granulare Flags geplant)
- Geometry Access (zero-copy geplant)
- Type Safety (ElementType enum geplant)

### test_energy_methods.cpp
**Status**: ✅ Production-Ready
**Zweck**: Validierung aller Energy Methods gegen Referenz-Energien

Test-Molekül: `A.xyz` (117 Atome)

**Getestete Methoden**:
- QM: gfn2, gfn1, ipea1, pm6, pm3, eht
- Force Fields: uff, gfnff
- Providers: TBLite, Ulysses, XTB

**Referenz-Energien** (Beispiele):
```cpp
m_reference_energies["gfn2"] = -165.75590286;  // TBLite GFN2-xTB
m_reference_energies["uff"] = 1.25494377;      // UFF
m_reference_energies["eht"] = -190.28490972;   // Extended Hückel
```

**Toleranzen**:
- QM Methods: 1e-6 Eh (hochpräzise)
- Force Fields: 1e-5 Eh (Force Field Toleranz)

## Integration Tests

### AAAbGal.cpp - RMSD Multi-Method Test
**Zweck**: RMSD mit verschiedenen Alignment-Methoden testen

**Migriert (2025-10-19)**: `RMSDJson` → `ParameterRegistry::getInstance().getDefaultJson("rmsd")`

**Getestete Methoden**:
1. `AAAbGal_dtemplate()` - Dimer template
2. `AAAbGal_free()` - Free reordering
3. `AAAbGal_template()` - Template with nomunkres
4. `AAAbGal_subspace()` - Subspace search
5. `AAAbGal_incr()` - Incremental

**Referenz-RMSD**: 0.457061 (mit Reordering)

### confscan.cpp - ConfScan Multi-Scenario Test
**Zweck**: ConfScan mit verschiedenen Methoden testen

**Migriert (2025-10-19)**: `ConfScanJson` → `ParameterRegistry::getInstance().getDefaultJson("confscan")`

**Getestete Szenarien** (10):
1. `free()` - Free alignment
2. `subspace()` - Subspace search
3. `template_method()` - Template-based
4. `dtemplate()` - Dimer template
5. `molalign()` - External molalign
6. `sLX1()`, `sLX2()`, `sLX2Reset()` - sLX logic
7. `sLX20()`, `sLX20Reset()` - Extended sLX

## Development Guidelines

### Adding New CLI Tests

1. **Erstelle Test-Verzeichnis**:
   ```bash
   mkdir -p test_cases/cli/capability/XX_test_name
   cd test_cases/cli/capability/XX_test_name
   ```

2. **Kopiere Template**:
   ```bash
   cp ../../rmsd/01_default_rmsd/run_test.sh .
   ```

3. **Passe an**:
   - Test-Name und Beschreibung
   - Curcuma Command
   - Validierungs-Logik

4. **Registriere in CMakeLists.txt**:
   ```cmake
   add_test(NAME cli_capability_XX ...)
   ```

5. **Dokumentiere Referenzwerte** in GOLDEN_REFERENCES.md

### Adding New Unit Tests

1. **Erstelle `test_newfeature.cpp`**
2. **Folge test_energy_methods.cpp Pattern**:
   - Test class mit Setup
   - Reference values
   - Tolerance definition
   - Result reporting
3. **Registriere in CMakeLists.txt**:
   ```cmake
   add_executable(test_newfeature test_newfeature.cpp)
   target_link_libraries(test_newfeature curcuma_core ...)
   add_test(NAME test_newfeature COMMAND test_newfeature)
   ```

## Testing Best Practices

### CLI Tests
✅ **DO**:
- Teste reale User-Workflows
- Validiere wissenschaftliche Korrektheit
- Nutze Golden References mit Toleranzen
- Dokumentiere erwartete Fehler-Szenarien
- Cleanup generierte Dateien

❌ **DON'T**:
- Verlasse dich nur auf Exit-Codes
- Teste interne Implementation-Details
- Hardcode absolute Pfade
- Ignoriere wissenschaftliche Validierung

### Unit Tests
✅ **DO**:
- Teste einzelne Funktionen isoliert
- Nutze klare Referenz-Werte
- Dokumentiere Toleranzen wissenschaftlich
- Teste Edge-Cases

❌ **DON'T**:
- Teste mehrere Features in einem Test
- Verwende undokumentierte Magic Numbers
- Überspringe Error-Cases

## Known Issues (Stand 2025-10-19)

### CRITICAL: curcumaopt JSON null-Fehler
**Betrifft**: 11/26 CLI Tests (42%)
```
[ERROR] Optimization failed with lbfgspp:
[json.exception.type_error.306] cannot use value() with null
```

**Impact**:
- ❌ 5 curcumaopt Tests
- ❌ 6 simplemd Tests

**Root Cause**: Wahrscheinlich ConfigManager/ParameterRegistry Migration
**Location**: `src/capabilities/curcumaopt.cpp` oder `optimiser/`
**Priority**: HÖCHSTE - Blockiert 42% der Tests

### Minor Issues
- 3 "invalid_method" Tests: Expected failure scenarios noch nicht implementiert
- SimpleMD wall potential: Ausgabedateien nicht generiert (hängt mit opt-Bug zusammen)

## Instructions Block

**PRESERVED - DO NOT EDIT BY CLAUDE**

### Future Development
- [ ] Fix curcumaopt JSON null-Fehler (HÖCHSTE PRIORITÄT)
- [ ] Erweitere wissenschaftliche Validierung für alle Tests
- [ ] Implementiere "expected failure" Pattern für invalid_method Tests
- [ ] Füge Performance-Benchmarks hinzu (Regression Tests)
- [ ] Erweitere test_molecule.cpp für geplantes SOA/AOS Refactoring

### Testing Philosophy
- **Test-Driven Development**: Schreibe Tests vor Refactorings
- **Scientific Accuracy**: Validiere Physik, nicht nur Code
- **Regression Prevention**: Capture current behavior als Golden Reference
- **Documentation**: Tests sind lebende Dokumentation der Features

---

*Diese Dokumentation beschreibt die Test-Infrastruktur und Best Practices für Curcuma*
