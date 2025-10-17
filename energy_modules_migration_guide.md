# Migrations-Guide: EnergyCalculator, QM- und FF-Module

Dieses Dokument beschreibt die Migration der zentralen Energie-Berechnungs-Module auf das neue `ConfigManager`-System.

## Motivation

Die Kern-Module von Curcuma (`EnergyCalculator`, `XTBInterface`, `TBLiteInterface`, `ForceFieldMethod` etc.) verwenden derzeit eine Mischung aus `MergeJson`-Aufrufen und direkter Zuweisung von `json`-Objekten, um ihre Konfiguration zu verwalten. Dies führt zu inkonsistentem Code und erschwert die Nachverfolgung von Parametern.

Durch die Umstellung dieser fundamentalen Schichten auf den `ConfigManager` wird das gesamte System von der höchsten bis zur untersten Ebene vereinheitlicht.

## Allgemeine Strategie

Die Migration folgt einem klaren Muster:

1.  **Parameter-Definition**: Jedes Modul (z.B. `xtb`, `dftd3`) erhält einen `PARAMETER_DEFINITION_BLOCK` in seiner Header-Datei, der alle spezifischen Einstellungen enthält.
2.  **`ConfigManager` im Konstruktor**: Anstatt eines `json`-Objekts wird den Konstruktoren der Module eine `ConfigManager`-Instanz übergeben.
3.  **Typsicherer Zugriff**: Alle internen Zuweisungen (z.B. `m_accuracy = m_xtbsettings["xtb_ac"]`) werden durch `m_config.get<T>()`-Aufrufe ersetzt.
4.  **`EnergyCalculator` als Verteiler**: Der `EnergyCalculator` wird dafür verantwortlich sein, die korrekte, untergeordnete Konfiguration an die von ihm erstellten Methoden weiterzugeben.

## Phase 1: Parameter-Definitionen

Füge die entsprechenden Blöcke in die jeweiligen Header-Dateien ein.

**Beispiel für `xtbinterface.h`:**
```cpp
// In src/core/energy_calculators/qm_methods/xtbinterface.h

BEGIN_PARAMETER_DEFINITION(xtb)
    PARAM(accuracy, Double, 1.0, "Accuracy for xTB calculations.", "SCF", {"xtb_ac"})
    PARAM(max_iterations, Int, 200, "Max SCF iterations.", "SCF", {"SCFmaxiter"})
    PARAM(electronic_temperature, Double, 300.0, "Electronic temperature in Kelvin.", "SCF", {"Tele"})
    PARAM(spin, Double, 0.0, "Total spin of the system.", "SCF", {})
END_PARAMETER_DEFINITION
```

**Beispiel für `dftd3interface.h`:**
```cpp
// In src/core/energy_calculators/qm_methods/dftd3interface.h

BEGIN_PARAMETER_DEFINITION(dftd3)
    PARAM(version, Int, 3, "Version of the D3 correction (3 or 4).", "General", {"d_version"})
    PARAM(damping, Int, 2, "Damping function (0: zero, 1: bj, 2: z-damping).", "Damping", {"d_damping"})
    PARAM(s6, Double, 1.0, "Scaling factor for the C6 term.", "Scaling", {"d_s6"})
    PARAM(s8, Double, 0.78, "Scaling factor for the C8 term.", "Scaling", {"d_s8"})
    // ... weitere D3/D4 Parameter ...
END_PARAMETER_DEFINITION
```

## Phase 2: Refactoring der Implementierungen

### 2.1. `XTBInterface` (und ähnliche QM-Interfaces)

**Vorher (`xtbinterface.cpp`):**
```cpp
XTBInterface::XTBInterface(const json& xtbsettings)
{
    m_xtbsettings = MergeJson(XTBSettings, xtbsettings);
    m_accuracy = m_xtbsettings["xtb_ac"];
    m_SCFmaxiter = m_xtbsettings["SCFmaxiter"];
    // ...
}
```

**Nachher (`xtbinterface.cpp`):**
```cpp
// Der Konstruktor erhält nun einen ConfigManager
XTBInterface::XTBInterface(const ConfigManager& config)
    : m_config(config) // ConfigManager als Member speichern
{
    // Direkter, typsicherer Zugriff
    m_accuracy = m_config.get<double>("accuracy");
    m_SCFmaxiter = m_config.get<int>("max_iterations");
    // ...
}
```

### 2.2. `EnergyCalculator`

Der `EnergyCalculator` muss angepasst werden, um den `ConfigManager` an die von ihm erstellten Methoden weiterzugeben.

**Vorher (`energycalculator.cpp`):**
```cpp
bool EnergyCalculator::createMethod(const std::string& method_name, const json& config) {
    // ...
    m_method = MethodFactory::create(method_name, method_config);
    // ...
}
```

**Nachher (`energycalculator.cpp`):**
```cpp
bool EnergyCalculator::createMethod(const std::string& method_name, const ConfigManager& config) {
    // ...
    // Erstelle einen spezifischen ConfigManager für das Sub-Modul
    ConfigManager method_config(method_name, config.exportConfig()); 
    m_method = MethodFactory::create(method_name, method_config);
    // ...
}

// Die Factory-Funktion wird ebenfalls angepasst:
// std::unique_ptr<ComputationalMethod> MethodFactory::create(const std::string& name, const ConfigManager& config)
```

## 📊 PHASE 3A: Abgeschlossen ✅

### Status: Parameter-Definitionen erfolgreich implementiert

**Zeitstempel**: Januar 2025
**Parameter insgesamt**: 240 Parameter über 12 Module
**Build-Verifikation**: `make GenerateParams` ✅

### Implementierte Module

#### QM-Methoden (7 Module, 54 Parameter insgesamt)
| Modul | File | Parameter | Status |
|-------|------|-----------|--------|
| XTB | `xtbinterface.h` | 4 (accuracy, max_iterations, electronic_temperature, spin) | ✅ |
| TBLite | `tbliteinterface.h` | 13 (SCF, Solvation: CPCM, GB, ALPB) | ✅ |
| DFT-D3 | `dftd3interface.h` | 10 (C6/C8 Scaling, Damping, ATM) | ✅ |
| DFT-D4 | `dftd4interface.h` | 9 (D4-spezifisch, s10 Term) | ✅ |
| Ulysses | `ulyssesinterface.h` | 5 (Method, SCF, Multiplicity, Solvation) | ✅ |
| GFN-FF (extern) | `gfnffinterface.h` | 3 (Charge, Print, Solvent) | ✅ |
| ORCA | `orcainterface.h` | 10 (TB, SCF, Damping, Initial Guess, Solvation) | ✅ |

#### FF-Methoden (1 Modul, 36 Parameter)
| Modul | File | Parameter | Status |
|-------|------|-----------|--------|
| ForceFieldGenerator | `forcefieldgenerator.h` | 36 (Method, Dispersion, Scaling, H4, Force Constants) | ✅ |

### PHASE 3B: Konstruktor-Migration ✅ (ABGESCHLOSSEN)

**Status**: Januar 2025 - Erfolgreich abgeschlossen
**Umfang**: Alle 8 QM/FF-Interface-Konstruktoren + 8 Method-Wrapper + 3 Generator-Implementierungen

**Implementierte Komponenten**:
1. ✅ Interface-Konstruktoren (8 Dateien):
   - XTBInterface, TBLiteInterface, DFTD3Interface, DFTD4Interface
   - UlyssesInterface, GFNFFInterface, OrcaInterface

2. ✅ Method-Wrapper Updates (8 Dateien):
   - xtb_method.cpp, tblite_method.cpp, dispersion_method.cpp
   - ulysses_method.cpp, external_gfnff_method.cpp, eht_method.cpp
   - gfnff_method.cpp, forcefield_method.cpp

3. ✅ ForceFieldGenerator Migration:
   - forcefieldgenerator.h/.cpp - Constructor updated
   - forcefield_method.cpp - Updated FFGenerator creation
   - qmdfffit.cpp - Updated FFGenerator instantiation

**Verifizierung**: ✅ Build successful, Unit tests passing, End-to-end verification complete

---

### PHASE 3C: EnergyCalculator-Integration & End-to-End-Tests ✅ (ABGESCHLOSSEN)

**Status**: Oktober 2025 - Erfolgreich abgeschlossen

**Implementierte Komponenten**:
1. ✅ ConfigManager Constructor Overloads (2 neue Konstruktoren)
   - `EnergyCalculator(const std::string& method, const ConfigManager& config)`
   - `EnergyCalculator(const std::string& method, const ConfigManager& config, const std::string& basename)`

2. ✅ Delegating Constructor Pattern
   - Alte json-basierte Konstruktoren delegieren an ConfigManager-Versionen
   - Vollständige Rückwärtskompatibilität: Alle 12 Capabilities funktionieren unchanged
   - Einzelner Konvertierungspunkt verhindert Code-Duplikation

3. ✅ initializeCommonFromConfig Implementation
   - Neue Methode handles ConfigManager-basierte Initialization
   - Verwendet `config.exportConfig()` für interne JSON-Kompatibilität
   - Seamless Integration mit bestehender createMethod() Infrastruktur

4. ✅ Cleanup: Entfernte statische EnergyCalculatorJson-Deklaration
   - War nicht verwendet
   - Reduziert statisches Speicher-Overhead

**Verifizierung**: ✅ Build erfolgreich (100%), Backward Compatibility Test bestätigt (gfn2 Berechnung erfolgreich)

---

## Zusammenfassung der Vorteile

-   **Konsistenz**: Alle Berechnungsmethoden, von den high-level Capabilities bis zu den low-level QM-Engines, verwenden dieselbe Konfigurations-API.
-   **Entkopplung**: Die QM/FF-Interfaces müssen nichts mehr über JSON wissen. Sie arbeiten nur noch mit dem `ConfigManager`.
-   **Klarheit**: Die `...Settings.h`-Dateien mit den statischen JSON-Objekten werden überflüssig und durch die klareren `PARAM`-Blöcke in den Headern ersetzt.
-   **Wartbarkeit**: Das Hinzufügen einer neuen Einstellung zu z.B. `xtb` erfordert nur noch eine Zeile im `PARAM`-Block in `xtbinterface.h`, und der Rest des Systems (Hilfe, Default-Werte) wird automatisch aktualisiert.
