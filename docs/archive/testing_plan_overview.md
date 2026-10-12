# Test-Plan: Übergreifende Strategie für die Curcuma Test-Suite

## 1. Ziele und Motivation

Das primäre Ziel ist die Entwicklung einer robusten, automatisierten Test-Suite für `curcuma`. Diese Suite soll sicherstellen, dass das Programm wie erwartet funktioniert und zukünftige Änderungen oder Refactorings keine bestehenden Funktionalitäten beeinträchtigen (Regressionsschutz).

**Hauptziele:**
- **Validierung aller CLI-Argumente:** Jeder für den Benutzer verfügbare Parameter muss auf korrekte und fehlerhafte Eingaben getestet werden.
- **Korrektheit der Ergebnisse:** Sicherstellen, dass die wissenschaftlichen Berechnungen für definierte Testfälle korrekte Ergebnisse liefern.
- **Stabilität:** Gewährleisten, dass das Programm bei unerwarteten Eingaben nicht abstürzt, sondern kontrolliert mit einer verständlichen Fehlermeldung endet.
- **Automatisierung:** Die gesamte Test-Suite muss über einen einzigen Befehl ausgeführt und in den Continuous Integration (CI) Prozess integriert werden können.

## 2. Vorgeschlagene Werkzeuge

- **Test-Runner:** **CTest**. Als Teil von CMake ist CTest das natürliche Werkzeug, um Tests zu definieren, zu organisieren und auszuführen.
- **Test-Framework (für Unit-/Integration-Tests):** **GoogleTest (gtest)**. Ein mächtiges und weit verbreitetes Framework für C++, das eine reiche Auswahl an Assertions und Features bietet. Es lässt sich nahtlos in CMake/CTest integrieren.
- **CLI-Test-Skripting:** **Bash-Skripte**. Für reine End-to-End-CLI-Tests sind einfache, portable Shell-Skripte ideal. Sie können den `curcuma`-Binary ausführen, Exit-Codes prüfen und Ausgaben mit `grep` und `diff` analysieren.

## 3. Verzeichnisstruktur

Alle Tests sollen in einem neuen, dedizierten `tests/`-Verzeichnis im Root des Projekts angesiedelt werden.

```
curcuma/
├── src/
├── tests/
│   ├── CMakeLists.txt         # Haupt-CMake-Datei für alle Tests
│   ├── cli/                   # End-to-End CLI-Tests
│   │   ├── confscan/
│   │   │   ├── 01_basic_run/
│   │   │   │   ├── run_test.sh
│   │   │   │   └── input.xyz
│   │   │   └── 02_invalid_rmsd/
│   │   ├── simplemd/
│   │   └── ...
│   ├── core/                  # Unit-Tests für Kern-Klassen (z.B. ConfigManager)
│   │   ├── test_config_manager.cpp
│   │   └── ...
│   ├── capabilities/          # Unit-/Integration-Tests für Capabilities
│   │   ├── test_simplemd_thermostats.cpp
│   │   └── ...
│   └── test_data/             # "Golden"-Referenzdateien für Vergleiche
│       ├── confscan_basic_output.xyz
│       └── ...
└── ...
```

## 4. Arten von Tests

### Typ 1: Unit-Tests (Fokus: Isolation)
- **Ziel:** Testen einer einzelnen Klasse oder Funktion isoliert von ihren Abhängigkeiten.
- **Beispiel:** Testen der `Tools::String2DoubleVec`-Funktion mit verschiedenen Eingabe-Strings.
- **Implementierung:** In `tests/core/` oder `tests/capabilities/` unter Verwendung von GoogleTest.

### Typ 2: Integration-Tests (Fokus: Zusammenspiel)
- **Ziel:** Testen, wie mehrere Komponenten zusammenarbeiten.
- **Beispiel:** Testen, ob der `EnergyCalculator` bei der Methode `gfn2-d3` korrekt die `XTBInterface` und die `DFTD3Interface` instanziiert und aufruft.
- **Implementierung:** In `tests/capabilities/` unter Verwendung von GoogleTest, wobei mehrere "echte" Klassen zusammenarbeiten.

### Typ 3: End-to-End (CLI) Tests (Fokus: Benutzer-Sicht)
- **Dies ist der Kern der Anforderung zur Validierung der CLI-Argumente.**
- **Ziel:** Testen des gesamten Programms aus der Sicht des Benutzers. Der kompilierte `curcuma`-Binary wird mit einer Reihe von Kommandozeilen-Argumenten aufgerufen.
- **Implementierung:**
    1. Für jeden Testfall wird ein Verzeichnis unter `tests/cli/<capability>/` erstellt.
    2. Jedes Verzeichnis enthält ein `run_test.sh`-Skript.
    3. Das Skript führt `curcuma` mit spezifischen Argumenten aus und leitet `stdout` und `stderr` in Log-Dateien um.
    4. **Assertions im Skript:**
        - **Exit-Code prüfen:** `if [ $? -ne 0 ] ...`
        - **Ausgaben prüfen:** `grep "Error: Invalid value" stderr.log`
        - **Dateien vergleichen:** `diff output.xyz /path/to/test_data/golden_file.xyz`
    5. Das `CMakeLists.txt` in `tests/` registriert jedes `run_test.sh` als einen Testfall bei CTest mittels `add_test()`.

## 5. Integration in den Build-Prozess

- Das Haupt-`CMakeLists.txt` wird um `add_subdirectory(tests)` erweitert.
- Der `GenerateParams`-Schritt ist eine wichtige Abhängigkeit, da die Tests die korrekte Registrierung der Parameter prüfen.
- Der CI-Workflow (z.B. `.github/workflows/build.yml`) wird erweitert, um nach dem Build-Schritt den Befehl `ctest --output-on-failure` auszuführen. Ein Fehlschlagen eines Tests führt zum Fehlschlagen des gesamten CI-Laufs.

## 6. Roadmap

1.  **Setup:** Einrichten der Verzeichnisstruktur, Integration von GoogleTest in das CMake-Projekt.
2.  **Planung:** Erstellung der thematischen Test-Pläne für jede Capability (siehe separate Markdown-Dateien).
3.  **Implementierung (CLI):** Schrittweise Implementierung der in den Plänen definierten `run_test.sh`-Skripte, beginnend mit den wichtigsten Anwendungsfällen und Fehlerfällen für jede Capability.
4.  **Implementierung (Unit/Integration):** Parallel dazu können Unit-Tests für kritische Kernkomponenten (insbesondere `ConfigManager`) und komplexe Algorithmen (z.B. Thermostate) entwickelt werden.
5.  **CI-Integration:** Einbindung des `ctest`-Aufrufs in die GitHub-Actions.
