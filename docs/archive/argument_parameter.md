# Curcuma Argument and Parameter Handling: Inventory and Consolidation Plan

This document provides a comprehensive overview of the current state of command-line argument and parameter handling in `curcuma` and proposes a plan for a unified, maintainable, and user-friendly system.

## 1. Bestandsaufnahme (Inventory)

The current system for handling command-line arguments and parameters is functional but suffers from several inconsistencies that make it difficult to maintain, extend, and understand.

### 1.1. Argumenten-Verarbeitung (CLI Parsing)

- **Zentraler Parser**: Die Funktion `CLI2Json` in `src/main.cpp` ist der zentrale Einstiegspunkt. Sie konvertiert Kommandozeilen-Argumente (`argc`, `argv`) in ein verschachteltes `nlohmann::json` Objekt (`controller`).
- **Struktur**: Der erste Parameter (z.B. `-casino`) wird als Haupt-Schlüsselwort für das Submodul verwendet. Alle folgenden `-key value` Paare werden in ein untergeordnetes JSON-Objekt für dieses Schlüsselwort eingefügt.
- **Beispiel**: `curcuma -casino -steps 10000 -temp 300.0` wird zu:
  ```json
  {
    "casino": {
      "steps": 10000,
      "temperature": 300.0
    }
  }
  ```
- **Dot-Notation**: Die Hilfsfunktion `setNestedJsonValue` ermöglicht eine Punktnotation für verschachtelte Parameter (z.B. `-topological.save_image true`), was eine tiefere Strukturierung erlaubt.

### 1.2. Inkonsistenzen und Probleme

#### a) Uneinheitliche Groß-/Kleinschreibung (Argument Casing)

Die Benennung der JSON-Schlüssel für Parameter ist uneinheitlich. Es gibt eine Mischung aus `camelCase`, `PascalCase` und `snake_case`.

| Modul (`.cpp`) | Beispiele für Schlüssel                               | 
| :--------------- | :---------------------------------------------------- |
| `simplemd`       | `MaxTime`, `rmrottrans`, `noCOLVARfile`, `hmass`        |
| `casino`         | `steps`, `temperature`, `step_size`, `move_type`      |
| `curcumaopt`     | `SinglePoint`, `MaxIter`, `GradNorm`, `LBFGS_m`       |
| `confscan`       | `forceReorder`, `maxenergy`, `sLX`, `UseOrders`       |
| `docking`        | `Pos_X`, `Step_X`, `AutoPos`, `CentroidMaxDistance`   |
| `hessian`        | `hess_calc`, `hess_read_file`, `finite_diff_step`    |

#### b) Inkonsistenter Zugriff auf Submodul-Konfigurationen

Submodule greifen auf unterschiedliche Weise auf ihre Konfigurationen zu.

- **Direkter Merge**: `SimpleMD` und `Casino` erhalten das `controller`-Objekt und erwarten ihre Parameter direkt unter dem Haupt-Schlüsselwort (z.B. `controller["casino"]`).
- **Expliziter verschachtelter Zugriff**: `executeOptimization` in `main.cpp` greift explizit auf `controller["opt"]` zu.
- **Manuelle Umstrukturierung**: Der `-compare` Befehl in `main.cpp` ist ein extremer Sonderfall, bei dem der Inhalt von `controller["compare"]` manuell in ein neues Objekt unter den Schlüssel `rmsd` kopiert wird.
- **Workaround in `analysis.cpp`**: Das `UnifiedAnalysis`-Modul muss verschachtelte Parameter manuell aus `controller["analysis"]` in seine eigene Konfiguration `m_config` kopieren, da die `MergeJson`-Funktion dies nicht automatisch für verschachtelte Objekte tut.

#### c) Manuelles und fehleranfälliges Hilfe-System

Die Hilfe-Texte (`printHelp()`) sind in den meisten Modulen manuell als lange, formatierte String-Literale implementiert.

- **Wartungsaufwand**: Jede Änderung, jedes Hinzufügen oder Entfernen eines Parameters erfordert eine manuelle Anpassung des Hilfe-Strings.
- **Fehleranfälligkeit**: Es gibt keine Garantie, dass die Hilfe-Texte mit den tatsächlich im Code verwendeten Parametern (Standardwerte, Typen) übereinstimmen.
- **Inkonsistentes Format**: Das Layout und der Detailgrad der Hilfe-Texte variieren stark zwischen den Modulen.

#### d) Mangelnde Zentralisierung der Parameter-Definition

Parameter werden dezentral an verschiedenen Stellen definiert und verwendet:
1.  Als hartcodierte Strings in `Json2KeyWord<T>(...)`.
2.  Als Schlüssel in den `...Json`-Standardkonfigurationsobjekten.
3.  Als manuell geschriebene Texte in den `printHelp()`-Funktionen.

Diese Trennung macht es unmöglich, eine "Single Source of Truth" für einen Parameter zu haben.

---

## 2. Konsolidierungsplan (Version 4.1 - Finale Architektur)

Diese finale Version des Plans integriert das Feedback, eine vollständig C++-basierte Lösung ohne externe Skript-Abhängigkeiten zu schaffen und die Modulzugehörigkeit explizit zu machen.

### 2.1. Parameter-Definition direkt im Code mittels Makros

Die Parameter werden direkt in den `.h`-Dateien der jeweiligen Module in einem speziell markierten Block deklariert. Der Modulname wird dabei explizit übergeben.

**Regel: Alle Bezeichner (Parameter-Namen, Aliase, Kategorien) müssen auf Englisch sein, um die Konsistenz im Code zu gewährleisten.**

**Beispiel (`src/capabilities/casino.h`):**
```cpp
class Casino : public CurcumaMethod {
public:
    // ...

private:
    // vvvvvvvvvvvv PARAMETER DEFINITION BLOCK vvvvvvvvvvvv
    #define BEGIN_PARAMETER_DEFINITION(casino)

    PARAM(steps, Int, 10000, "Number of MC steps to perform.", "Basic", { "nstep1" })
    PARAM(temperature, Double, 300.0, "Simulation temperature in Kelvin.", "Basic", { "temp" })

    #define END_PARAMETER_DEFINITION
    // ^^^^^^^^^^^^ PARAMETER DEFINITION BLOCK ^^^^^^^^^^^^
};
```

**Beispiel für geteilten Parameter (`src/capabilities/simplemd.h`):**
```cpp
class SimpleMD : public CurcumaMethod {
private:
    #define BEGIN_PARAMETER_DEFINITION(simplemd)

    PARAM(max_time, Double, 1000.0, "Maximum simulation time in fs.", "Basic", { "MaxTime" })
    // Der Parser erkennt, dass "temperature" auch in anderen Modulen definiert ist
    PARAM(temperature, Double, 300.0, "Simulation temperature in Kelvin.", "Basic", { "temp" })

    #define END_PARAMETER_DEFINITION
};
```

### 2.2. C++ Präprozessor-Tool (`curcuma_param_parser`)

Ein kleines, eigenständiges C++-Kommandozeilen-Tool wird als Teil des Projekts erstellt.

-   **Aufgabe**: Das Tool parst alle relevanten `.h`-Dateien des Projekts.
-   **Logik**:
    1.  Es sucht nach `BEGIN_PARAMETER_DEFINITION(module_name)` und extrahiert den `module_name`.
    2.  Innerhalb des Blocks parst es die `PARAM(...)`-Makros und extrahiert alle Metadaten.
    3.  Es baut den hierarchischen Pfad (z.B. `casino.steps`).
    4.  **Automatische Erkennung**: Wenn das Tool denselben Parameternamen in mehreren Modulen findet, erkennt es diesen automatisch als geteilten Parameter.
-   **Output**: Das Tool generiert die finale C++ Header-Datei `generated_parameter_registry.h`.

### 2.3. Integration in den Build-Prozess (CMake)

Der Prozess wird vollständig in `CMakeLists.txt` automatisiert.

```cmake
# 1. Eigenes Executable für den Parser definieren und bauen
add_executable(curcuma_param_parser scripts/param_parser/main.cpp)

# 2. Alle relevanten Header-Dateien des Projekts für den Parser sammeln
file(GLOB_RECURSE ALL_PROJECT_HEADERS "src/*.h" "src/**/*.h")

# 3. Custom Command definieren, der den Parser ausführt, um die Registry zu generieren
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/generated/parameter_registry.h
    COMMAND curcuma_param_parser 
            --output ${CMAKE_CURRENT_BINARY_DIR}/generated/parameter_registry.h 
            --inputs ${ALL_PROJECT_HEADERS}
    DEPENDS curcuma_param_parser ${ALL_PROJECT_HEADERS}
    COMMENT "Generating C++ parameter registry from source headers"
)

# 4. Ein Target für den Generierungsschritt erstellen, um Abhängigkeiten zu managen
add_custom_target(GenerateParams DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/generated/parameter_registry.h)

# 5. Das Haupt-Target (curcuma_core) von der Generierung abhängig machen
add_dependencies(curcuma_core GenerateParams)

# 6. Das generierte Verzeichnis zu den Include-Pfaden hinzufügen
include_directories(${CMAKE_CURRENT_BINARY_DIR}/generated)
```

### 2.4. Finale Architektur und Vorteile

Das Endergebnis ist ein System, bei dem:

1.  **Parameter direkt beim Code definiert werden**, was die Übersicht und Wartung für Entwickler stark vereinfacht.
2.  Ein **C++-Tool die Konsistenz sicherstellt**, indem es eine einzige, autoritative Header-Datei generiert.
3.  **Keine externen Abhängigkeiten** wie Python oder zusätzliche Konfigurationsdateien für den Build-Prozess oder das Deployment nötig sind.
4.  Die **Hilfe-Texte und Standardwerte immer synchron** mit der Implementierung sind, da sie aus derselben Quelle generiert werden.
5.  **Komplexe Beziehungen** wie Vererbung und gemeinsame Nutzung zentral und automatisch verwaltet werden können.
