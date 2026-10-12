# Test-Plan: Capability `confscan`

## 1. Ziele

Dieser Test-Plan fokussiert auf die Validierung der `confscan`-Capability. Die Tests sollen die korrekte Funktionsweise des Konformer-Scans, die Filter-Mechanismen (RMSD, Energie) und die komplexe Parameter-Logik (z.B. für `sLX`, `rmsd_element`) sicherstellen.

## 2. Test-Szenarien (CLI-Tests)

Die Tests werden im Verzeichnis `tests/cli/confscan/` implementiert.

--- 

### Szenario 1: Standard-Scan

- **Beschreibung:** Ein einfacher Konformer-Scan mit Standard-RMSD- und Energie-Schwellenwerten.
- **Verzeichnis:** `01_default_scan/`
- **Kommando:** `curcuma -i conformers.xyz -confscan`
- **Erwartungen:**
    - Exit-Code: `0`
    - `STDOUT`: Zeigt den Fortschritt des Scans und eine Zusammenfassung der gefundenen einzigartigen Konformere.
    - **Artefakte:** `conformers.accepted.xyz`, `conformers.rejected.xyz`.
- **Prüfung im Skript:**
    1. Prüfe Exit-Code `0`.
    2. `grep "Found [0-9]* unique conformers" stdout.log`
    3. `test -f conformers.accepted.xyz`
    4. `diff conformers.accepted.xyz golden_default.accepted.xyz`

--- 

### Szenario 2: Dynamischer RMSD-Schwellenwert

- **Beschreibung:** Testet die `get_rmsd`-Funktionalität.
- **Verzeichnis:** `02_get_rmsd/`
- **Kommando:** `curcuma -i conformers.xyz -confscan -confscan.get_rmsd true`
- **Erwartungen:**
    - Exit-Code: `0`
    - `STDOUT`: Eine Meldung, die den dynamisch bestimmten RMSD-Wert anzeigt.
- **Prüfung im Skript:**
    1. Prüfe Exit-Code `0`.
    2. `grep "Dynamically determined RMSD threshold: [0-9.]*" stdout.log`

--- 

### Szenario 3: Ungültiger Wert für `rmsd_method`

- **Beschreibung:** Stellt sicher, dass das Programm bei einer ungültigen RMSD-Methode fehlschlägt.
- **Verzeichnis:** `03_invalid_rmsd_method/`
- **Kommando:** `curcuma -i conformers.xyz -confscan -confscan.rmsd_method non_existent`
- **Erwartungen:**
    - Exit-Code: `> 0`
    - `STDERR`: Eine Fehlermeldung, die den ungültigen Wert für `rmsd_method` beanstandet.
- **Prüfung im Skript:**
    1. Prüfe, ob der Exit-Code nicht `0` ist.
    2. `grep "Error: Invalid value for parameter 'rmsd_method'" stderr.log`

--- 

### Szenario 4: Test der `sLX`-Logik

- **Beschreibung:** Testet die spezielle `sLX`-Parsing-Logik, die mehrere Schwellenwerte gleichzeitig setzt.
- **Verzeichnis:** `04_slx_logic/`
- **Kommando:** `curcuma -i conformers.xyz -confscan -confscan.sLX "1.5,2.5"`
- **Erwartungen:**
    - Exit-Code: `0`
    - `STDOUT`: Die Ausgabe (bei erhöhter Verbosity) sollte zeigen, dass die Schwellenwerte für Energie, Inertia und Ripser entsprechend `1.5,2.5` gesetzt wurden.
- **Prüfung im Skript:**
    1. Prüfe Exit-Code `0`.
    2. `grep "Loose energy thresholds set by sLX: 1.5, 2.5" stdout.log` (Annahme: Programm gibt dies bei Verbosity > 0 aus).

--- 

### Szenario 5: Hybrid-RMSD mit spezifischen Elementen

- **Beschreibung:** Testet die `rmsd_element`-Option, die eine Liste von Elementen als String akzeptiert.
- **Verzeichnis:** `05_hybrid_rmsd_elements/`
- **Kommando:** `curcuma -i conformers.xyz -confscan -confscan.rmsd_method hybrid -confscan.rmsd_element "7,8"`
- **Erwartungen:**
    - Exit-Code: `0`
    - `STDOUT`: Eine Meldung bestätigt, dass der Hybrid-RMSD-Modus mit den Elementen N und O verwendet wird.
- **Prüfung im Skript:**
    1. Prüfe Exit-Code `0`.
    2. `grep "Using hybrid RMSD with element templates: 7, 8" stdout.log`

--- 

### Szenario 6: Nur Schweratome verwenden (`heavy_only`)

- **Beschreibung:** Testet die `heavy_only` (bzw. `heavy`) Option.
- **Verzeichnis:** `06_heavy_only/`
- **Kommando:** `curcuma -i conformers.xyz -confscan -confscan.heavy_only true`
- **Erwartungen:**
    - Exit-Code: `0`
    - Das Ergebnis (Anzahl der Konformere und die `accepted.xyz`-Datei) sollte sich von Szenario 1 unterscheiden.
- **Prüfung im Skript:**
    1. Prüfe Exit-Code `0`.
    2. `diff conformers.accepted.xyz golden_heavy_only.accepted.xyz`

--- 

### Szenario 7: Neustart eines Scans (`restart`)

- **Beschreibung:** Testet die `restart`-Funktionalität.
- **Verzeichnis:** `07_restart_scan/`
- **Ablauf:**
    1. Führe einen unvollständigen Scan auf einer großen Datei aus und beende ihn vorzeitig.
    2. **Kommando:** `curcuma -i large_conformers.xyz -confscan -confscan.restart true`
- **Erwartungen:**
    - Exit-Code: `0`
    - `STDOUT`: Eine Meldung wie "Restarting scan, found existing .confscan_progress file." sollte erscheinen.
- **Prüfung im Skript:**
    1. Erstelle eine Dummy-Fortschrittsdatei.
    2. Führe den `restart`-Befehl aus.
    3. `grep "Restarting scan" stdout.log`

## 3. Unit-Tests

- **Descriptor-Parsing:** Die komplexe Logik zur Verarbeitung der `sLX`, `sLE`, `sLI`, `sLH`-Parameter könnte in eine eigene, testbare Funktion oder Klasse ausgelagert werden. Man könnte sie mit verschiedenen Eingabe-Strings ("default", "1.0", "1.5,2.5") aufrufen und die resultierenden Schwellenwert-Vektoren prüfen.
- **Filter-Logik:** Die Kern-Entscheidungslogik ("akzeptiere ich dieses Konformer oder nicht?") könnte isoliert mit Mock-Konformeren und verschiedenen Schwellenwerten getestet werden.
