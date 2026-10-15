# Test-Plan: Capability `rmsd`

## 1. Ziele

Dieser Plan definiert die Tests für die `rmsd`-Capability. Die Tests sollen alle wichtigen Alignment-Methoden, die Atom-Auswahl-Mechanismen und die Parameter für den Kuhn-Munkres-Algorithmus validieren. Ziel ist es, die Korrektheit der RMSD-Berechnungen für verschiedene Anwendungsfälle sicherzustellen.

## 2. Test-Szenarien (CLI-Tests)

Die Tests werden im Verzeichnis `tests/cli/rmsd/` implementiert.

--- 

### Szenario 1: Standard-RMSD-Berechnung

- **Beschreibung:** Eine einfache RMSD-Berechnung zwischen zwei Molekülen mit Standardeinstellungen (inkl. Reordering).
- **Verzeichnis:** `01_default_rmsd/`
- **Kommando:** `curcuma -rmsd -ref ref.xyz -target target.xyz`
- **Erwartungen:**
    - Exit-Code: `0`
    - `STDOUT`: Gibt den berechneten RMSD-Wert aus.
- **Prüfung im Skript:**
    1. Prüfe Exit-Code `0`.
    2. `grep "RMSD after alignment: [0-9.]*" stdout.log`
    3. Vergleiche den RMSD-Wert mit einem erwarteten Wert: `grep "RMSD after alignment: 1.2345" stdout.log` (Annahme eines bekannten Ergebnisses).

--- 

### Szenario 2: Deaktiviertes Reordering

- **Beschreibung:** Testet die `no_reorder`-Option für Moleküle mit identischer Atomreihenfolge.
- **Verzeichnis:** `02_no_reorder/`
- **Kommando:** `curcuma -rmsd -ref ref.xyz -target target_same_order.xyz -rmsd.no_reorder true`
- **Erwartungen:**
    - Exit-Code: `0`
    - `STDOUT`: Sollte bestätigen, dass kein Reordering durchgeführt wurde. Der RMSD-Wert wird sich von dem mit Reordering unterscheiden.
- **Prüfung im Skript:**
    1. Prüfe Exit-Code `0`.
    2. `grep "Reordering disabled by user" stdout.log`

--- 

### Szenario 3: Ungültige Alignment-Methode

- **Beschreibung:** Stellt sicher, dass das Programm bei einer ungültigen `method` fehlschlägt.
- **Verzeichnis:** `03_invalid_method/`
- **Kommando:** `curcuma -rmsd -ref ref.xyz -target target.xyz -rmsd.method non_existent`
- **Erwartungen:**
    - Exit-Code: `> 0`
    - `STDERR`: Eine Fehlermeldung, die den ungültigen Wert für `method` beanstandet.
- **Prüfung im Skript:**
    1. Prüfe, ob der Exit-Code nicht `0` ist.
    2. `grep "Error: Invalid value for parameter 'method'" stderr.log`

--- 

### Szenario 4: Template-Alignment mit Element-Auswahl

- **Beschreibung:** Testet die `template`-Methode mit einer Auswahl von Elementen als String.
- **Verzeichnis:** `04_template_elements/`
- **Kommando:** `curcuma -rmsd -ref ref.xyz -target target.xyz -rmsd.method template -rmsd.element "6,7"`
- **Erwartungen:**
    - Exit-Code: `0`
    - `STDOUT`: Bestätigt die Verwendung der Template-Methode mit den Elementen C und N.
- **Prüfung im Skript:**
    1. Prüfe Exit-Code `0`.
    2. `grep "Using template alignment with elements: 6, 7" stdout.log`

--- 

### Szenario 5: Fragment-basierte RMSD-Berechnung

- **Beschreibung:** Testet die Berechnung des RMSD nur für spezifische Fragmente der Moleküle.
- **Verzeichnis:** `05_fragment_rmsd/`
- **Kommando:** `curcuma -rmsd -ref ref.xyz -target target.xyz -rmsd.fragment_reference 0 -rmsd.fragment_target 1`
- **Erwartungen:**
    - Exit-Code: `0`
    - `STDOUT`: Bestätigt, dass nur die angegebenen Fragmente verwendet wurden. Der RMSD-Wert wird sich deutlich von der Gesamt-Molekül-Berechnung unterscheiden.
- **Prüfung im Skript:**
    1. Prüfe Exit-Code `0`.
    2. `grep "Calculating RMSD between fragment 0 of reference and fragment 1 of target" stdout.log`
    3. Vergleiche den RMSD-Wert mit einem bekannten Referenzwert für diesen speziellen Fall.

--- 

### Szenario 6: Verwendung von Aliases (`reorder`)

- **Beschreibung:** Stellt sicher, dass der alte Alias `reorder` äquivalent zu `force_reorder` ist.
- **Verzeichnis:** `06_alias_reorder/`
- **Kommando:** `curcuma -rmsd -ref ref.xyz -target target.xyz -rmsd.reorder true`
- **Erwartungen:**
    - Das Verhalten sollte identisch zu einer Berechnung mit `-rmsd.force_reorder true` sein.
- **Prüfung im Skript:**
    1. Führe zwei Läufe aus, einen mit `reorder`, einen mit `force_reorder`.
    2. Vergleiche die `stdout.log`-Dateien. Sie sollten (bis auf Laufzeit) identisch sein.

## 3. Unit-Tests

- **Kuhn-Munkres-Algorithmus:** Die Implementierung des KM-Algorithmus (falls intern) ist ein perfekter Kandidat für Unit-Tests. Man kann ihn mit kleinen, bekannten Kosten-Matrizen füttern und prüfen, ob die korrekte Zuordnung mit der minimalen Summe gefunden wird.
- **Atom-Mapper/Fragmenter:** Die Logik, die Atom-Indizes basierend auf Fragment-Definitionen oder expliziten Listen (`reference_atoms`) extrahiert, könnte isoliert getestet werden, um sicherzustellen, dass sie immer die korrekten Atome aus einer `Molecule`-Struktur auswählt.
- **Alignment-Strategien:** Falls das Refactoring zum Strategy Pattern (wie im `refactoring_guide_strategy_pattern.md` vorgeschlagen) umgesetzt wird, kann jede Alignment-Strategie (`HungarianStrategy`, `InertiaStrategy` etc.) isoliert mit einfachen 2- oder 3-Atom-Systemen getestet werden, deren optimales Alignment analytisch bekannt ist.
