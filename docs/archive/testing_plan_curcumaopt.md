# Test-Plan: Capability `curcumaopt`

## 1. Ziele

Dieser Plan beschreibt die Tests zur Validierung aller CLI-Argumente und der Kernfunktionalität der `curcumaopt`-Capability. Die Tests sollen sicherstellen, dass sowohl Geometrieoptimierungen als auch Single-Point-Berechnungen korrekt konfiguriert und ausgeführt werden können.

## 2. Test-Szenarien (CLI-Tests)

Alle Tests werden als End-to-End-Tests im Verzeichnis `tests/cli/curcumaopt/` implementiert. Jeder Testfall erhält ein eigenes Unterverzeichnis und ein `run_test.sh`-Skript.

--- 

### Szenario 1: Standard-Optimierung (UFF)

- **Beschreibung:** Eine einfache Geometrieoptimierung mit dem Standard-Force-Field (UFF).
- **Verzeichnis:** `01_default_uff_opt/`
- **Kommando:** `curcuma -i input.xyz -opt`
- **Erwartungen:**
    - Exit-Code: `0`
    - `STDOUT`: Enthält die Konvergenz-Informationen und die finale Energie.
    - **Artefakte:** Eine optimierte `input.opt.xyz`-Datei wird erstellt.
- **Prüfung im Skript:**
    1. Prüfe, ob der Exit-Code `0` ist.
    2. `grep "Optimization has converged" stdout.log`
    3. `test -f input.opt.xyz`
    4. `diff input.opt.xyz golden_uff_opt.xyz` (Vergleich mit einer Referenz-Struktur)

--- 

### Szenario 2: Single-Point-Berechnung mit GFN2-xTB

- **Beschreibung:** Eine einzelne Energieberechnung (Single Point) mit der GFN2-xTB-Methode.
- **Verzeichnis:** `02_gfn2_single_point/`
- **Kommando:** `curcuma -i input.xyz -opt -opt.method gfn2 -opt.single_point true`
- **Erwartungen:**
    - Exit-Code: `0`
    - `STDOUT`: Enthält die Energie-Informationen für einen einzelnen Punkt, **keine** Konvergenz-Schleife.
    - **Artefakte:** Keine `.opt.xyz`-Datei, da keine Optimierung stattfindet.
- **Prüfung im Skript:**
    1. Prüfe Exit-Code `0`.
    2. `grep "Total Energy" stdout.log`
    3. `! grep "Optimization has converged" stdout.log` (Stelle sicher, dass keine Optimierung stattfand).
    4. `! test -f input.opt.xyz` (Stelle sicher, dass keine Optimierungs-Datei geschrieben wurde).

--- 

### Szenario 3: Ungültiger Parameterwert

- **Beschreibung:** Testet die Reaktion des Programms auf eine ungültige Methode.
- **Verzeichnis:** `03_invalid_method/`
- **Kommando:** `curcuma -i input.xyz -opt -opt.method non_existent_method`
- **Erwartungen:**
    - Exit-Code: `> 0` (Fehler)
    - `STDERR`: Eine klare Fehlermeldung, dass die Methode nicht existiert.
- **Prüfung im Skript:**
    1. Prüfe, ob der Exit-Code nicht `0` ist.
    2. `grep -i "Error: Method 'non_existent_method' not recognized" stderr.log`

--- 

### Szenario 4: LBFGS-Parameter-Übergabe

- **Beschreibung:** Stellt sicher, dass die verschachtelten LBFGS-Parameter korrekt an den Optimierer übergeben werden.
- **Verzeichnis:** `04_lbfgs_params/`
- **Kommando:** `curcuma -i input.xyz -opt -opt.optimizer 0 -opt.lbfgs.m 10`
- **Erwartungen:**
    - Exit-Code: `0`
    - `STDOUT`: Die Ausgabe sollte anzeigen, dass die LBFGS-Parameter verwendet werden (erfordert eventuell eine Anpassung der Verbosity im Programm).
- **Prüfung im Skript:**
    1. Prüfe Exit-Code `0`.
    2. `grep "LBFGS corrections to store: 10" stdout.log` (Annahme: Programm gibt bei höherer Verbosity die Parameter aus).

--- 

### Szenario 5: Verwendung eines Alias

- **Beschreibung:** Testet, ob ein alter Alias (`SinglePoint`) noch funktioniert und zur neuen Nomenklatur (`single_point`) äquivalent ist.
- **Verzeichnis:** `05_alias_singlepoint/`
- **Kommando:** `curcuma -i input.xyz -opt -opt.method gfn2 -opt.SinglePoint true`
- **Erwartungen:**
    - Das Verhalten muss exakt identisch zu Szenario 2 sein.
- **Prüfung im Skript:**
    1. Führe diesen Test und Test 2 aus.
    2. Vergleiche die `stdout.log`-Dateien beider Läufe. Sie sollten (bis auf Laufzeit-Informationen) identisch sein.

--- 

### Szenario 6: Optimierung nur von Wasserstoffatomen

- **Beschreibung:** Testet die `optH`-Funktionalität.
- **Verzeichnis:** `06_opth_only/`
- **Kommando:** `curcuma -i input.xyz -opt -opt.optimize_h true`
- **Erwartungen:**
    - Exit-Code: `0`
    - Die resultierende Geometrie in `input.opt.xyz` sollte nur veränderte H-Positionen im Vergleich zur Eingabe aufweisen, während die Schweratome fixiert sind.
- **Prüfung im Skript:**
    1. Prüfe Exit-Code `0`.
    2. `test -f input.opt.xyz`
    3. Ein spezielles Skript (z.B. Python mit `numpy`) könnte die XYZ-Dateien laden und verifizieren, dass die Koordinaten der Nicht-Wasserstoff-Atome zwischen `input.xyz` und `input.opt.xyz` unverändert sind.

## 3. Unit-Tests

Zusätzlich zu den CLI-Tests könnten Unit-Tests für folgende Bereiche sinnvoll sein:
- **Konvergenz-Logik:** Die `ConvergenceChecker`-Klasse könnte isoliert getestet werden, um sicherzustellen, dass sie bei verschiedenen Eingabe-Serien von Energien und Gradienten korrekt `true` oder `false` zurückgibt.
