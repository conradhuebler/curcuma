# Test-Plan: Capability `simplemd`

## 1. Ziele

Dieser Plan umreißt die Teststrategie zur Validierung der `simplemd`-Capability. Die Tests sollen alle wichtigen CLI-Parameter abdecken, von der grundlegenden Simulations-Einrichtung über Thermostate und Constraints bis hin zu den Wand-Potenzialen. Das Ziel ist es, die Korrektheit und Robustheit der Molekulardynamik-Simulationen zu gewährleisten.

## 2. Test-Szenarien (CLI-Tests)

Die Tests werden im Verzeichnis `tests/cli/simplemd/` implementiert.

--- 

### Szenario 1: Kurze NVE-Simulation (Standard)

- **Beschreibung:** Eine sehr kurze MD-Simulation im NVE-Ensemble (kein Thermostat), um die grundlegende Funktionsweise des Integrators zu testen.
- **Verzeichnis:** `01_nve_short_run/`
- **Kommando:** `curcuma -i input.xyz -simplemd -simplemd.max_time 10 -simplemd.thermostat none`
- **Erwartungen:**
    - Exit-Code: `0`
    - `STDOUT`: Status-Updates für die Simulation werden ausgegeben.
    - **Artefakte:** Eine Trajektorien-Datei `input.trj.xyz` wird erstellt.
- **Prüfung im Skript:**
    1. Prüfe Exit-Code `0`.
    2. `grep "STARTING MD" stdout.log`
    3. `test -f input.trj.xyz`
    4. Vergleiche die erstellte Trajektorie mit einer "Golden-Standard"-Datei: `diff input.trj.xyz golden_nve.trj.xyz`.

--- 

### Szenario 2: NVT-Simulation mit Berendsen-Thermostat

- **Beschreibung:** Testet die Verwendung des Berendsen-Thermostats.
- **Verzeichnis:** `02_nvt_berendsen/`
- **Kommando:** `curcuma -i input.xyz -simplemd -simplemd.max_time 10 -simplemd.thermostat berendsen -simplemd.temperature 300`
- **Erwartungen:**
    - Exit-Code: `0`
    - `STDOUT`: Die Ausgabe sollte bestätigen, dass der Berendsen-Thermostat aktiv ist.
    - **Artefakte:** `input.trj.xyz`.
- **Prüfung im Skript:**
    1. Prüfe Exit-Code `0`.
    2. `grep "Using Berendsen thermostat" stdout.log`
    3. `test -f input.trj.xyz`
    4. Die resultierende Trajektorie wird sich von der NVE-Simulation unterscheiden. Ein `diff` gegen eine `golden_nvt_berendsen.trj.xyz` ist erforderlich.

--- 

### Szenario 3: Ungültiger Thermostat-Name

- **Beschreibung:** Stellt sicher, dass das Programm bei einem unbekannten Thermostat-Namen mit einem Fehler abbricht.
- **Verzeichnis:** `03_invalid_thermostat/`
- **Kommando:** `curcuma -i input.xyz -simplemd -simplemd.thermostat non_existent_thermo`
- **Erwartungen:**
    - Exit-Code: `> 0`
    - `STDERR`: Eine klare Fehlermeldung über den ungültigen Wert.
- **Prüfung im Skript:**
    1. Prüfe, ob der Exit-Code nicht `0` ist.
    2. `grep "Error: Invalid value for parameter 'thermostat'" stderr.log`

--- 

### Szenario 4: Simulation mit RATTLE-Constraints

- **Beschreibung:** Testet die Aktivierung von RATTLE für alle Bindungen.
- **Verzeichnis:** `04_rattle_all/`
- **Kommando:** `curcuma -i input.xyz -simplemd -simplemd.max_time 5 -simplemd.rattle 1`
- **Erwartungen:**
    - Exit-Code: `0`
    - `STDOUT`: Die Ausgabe sollte die Aktivierung von RATTLE bestätigen.
- **Prüfung im Skript:**
    1. Prüfe Exit-Code `0`.
    2. `grep "RATTLE constraints will be applied" stdout.log`
    3. Ein `diff` der Trajektorie gegen eine `golden_rattle.trj.xyz`.

--- 

### Szenario 5: Simulation mit sphärischem Wand-Potenzial

- **Beschreibung:** Testet die `wall`-Funktionalität.
- **Verzeichnis:** `05_spheric_wall/`
- **Kommando:** `curcuma -i input.xyz -simplemd -simplemd.max_time 10 -simplemd.wall_type spheric -simplemd.wall_radius 10.0`
- **Erwartungen:**
    - Exit-Code: `0`
    - `STDOUT`: Ausgabe bestätigt das sphärische Wand-Potenzial.
- **Prüfung im Skript:**
    1. Prüfe Exit-Code `0`.
    2. `grep "Using spheric wall potential" stdout.log`
    3. Analyse der Trajektorie: Ein Python-Skript könnte prüfen, ob alle Atome innerhalb des Radius von 10.0 Å bleiben.

--- 

### Szenario 6: Testen eines Alias (z.B. `T` für `temperature`)

- **Beschreibung:** Stellt die Abwärtskompatibilität für wichtige Aliase sicher.
- **Verzeichnis:** `06_alias_temperature/`
- **Kommando:** `curcuma -i input.xyz -simplemd -simplemd.max_time 10 -simplemd.thermostat berendsen -simplemd.T 300`
- **Erwartungen:**
    - Das Ergebnis muss exakt identisch zu Szenario 2 sein.
- **Prüfung im Skript:**
    1. Führe diesen Test und Test 2 aus.
    2. Vergleiche die `stdout.log`-Dateien und die Trajektorien. Sie müssen identisch sein.

--- 

### Szenario 7: Neustart aus einer Restart-Datei

- **Beschreibung:** Testet die `write_restart` und `restart_file` Funktionalität.
- **Verzeichnis:** `07_restart_simulation/`
- **Ablauf:**
    1. **Lauf 1:** `curcuma -i input.xyz -simplemd -simplemd.max_time 5 -simplemd.write_restart_frequency 5`
    2. **Lauf 2:** `curcuma -i input.xyz -simplemd -simplemd.max_time 10 -simplemd.restart_file input.restart.json`
- **Erwartungen:**
    - Beide Läufe haben Exit-Code `0`.
    - Lauf 1 erzeugt eine `input.restart.json`-Datei.
    - Lauf 2 liest diese Datei und setzt die Simulation fort. Die `stdout.log` von Lauf 2 sollte eine entsprechende Meldung enthalten.
- **Prüfung im Skript:**
    1. Führe Lauf 1 aus, prüfe Exit-Code und die Existenz der `.restart.json`-Datei.
    2. Führe Lauf 2 aus, prüfe Exit-Code.
    3. `grep "Restarting simulation from" stdout_lauf2.log`

## 3. Unit-Tests

- **Thermostat-Logik:** Wie im Refactoring-Plan zum Strategy Pattern beschrieben, könnten die einzelnen Thermostat-Klassen (`BerendsenThermostat` etc.) isoliert mit einem Mock-`SystemState` getestet werden, um ihre mathematische Korrektheit zu verifizieren.
- **RATTLE-Algorithmus:** Die Kernlogik des RATTLE-Algorithmus könnte als separate Funktion oder Klasse extrahiert und mit bekannten Molekül-Zuständen und Constraints getestet werden.
