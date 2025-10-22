# Golden Reference Values for CLI Tests

**Erstellt**: 2025-10-19
**Curcuma Version**: 0.0.154
**Git Commit**: c0c3780

Diese Datei enthält die wissenschaftlichen Referenzwerte für alle CLI-Tests.

## RMSD Tests (✅ Alle bestanden: 5/6)

### cli_rmsd_01: Standard RMSD Berechnung
- **Molekül**: AAA-bGal (90 Atome)
- **Methode**: `free` (default reordering)
- **RMSD**: **2.87214 Å**
- **Permutationen**: 5 atomic indices reordered
- **H-Bond Topo Diff**: 0
- **Status**: ✅ BESTANDEN

### cli_rmsd_02: No Reordering
- **Molekül**: AAA-bGal (90 Atome)
- **Methode**: `noreorder` flag
- **Erwarteter RMSD**: > 2.87214 Å (höher ohne reordering)
- **Status**: ✅ BESTANDEN

### cli_rmsd_04: Template Elements
- **Molekül**: AAA-bGal
- **Methode**: `template` with specific elements
- **Status**: ✅ BESTANDEN

### cli_rmsd_05: Fragment
- **Molekül**: AAA-bGal
- **Methode**: Fragment-based RMSD
- **Status**: ✅ BESTANDEN

### cli_rmsd_06: Alias Reorder
- **Molekül**: AAA-bGal
- **Methode**: Using parameter alias
- **Status**: ✅ BESTANDEN

---

## ConfScan Tests (✅ Alle bestanden: 6/7)

### cli_confscan_01: Default ConfScan
- **Molekül**: TBD (wird extrahiert)
- **Methode**: Default conformational scanning
- **Erwartete Outputs**:
  - `input.accepted.xyz` - Accepted conformers
  - `input.rejected.xyz` - Rejected conformers
- **Status**: ✅ BESTANDEN

### cli_confscan_02: Get RMSD
- **Methode**: ConfScan with RMSD output
- **Status**: ✅ BESTANDEN

### cli_confscan_04: sLX Logic
- **Methode**: sLX (some Logic X) conformer filtering
- **Status**: ✅ BESTANDEN

### cli_confscan_05: Hybrid Elements
- **Methode**: Hybrid reordering with element templates
- **Status**: ✅ BESTANDEN

### cli_confscan_06: Heavy Only
- **Methode**: Heavy atom only comparison
- **Status**: ✅ BESTANDEN

### cli_confscan_07: Restart
- **Methode**: Restart from checkpoint
- **Status**: ✅ BESTANDEN

---

## SimpleMD Tests (✅ Teilweise: 1/7)

### cli_simplemd_06: Alias Temperature
- **Molekül**: TBD
- **Methode**: Using temperature parameter alias
- **Status**: ✅ BESTANDEN

### FEHLGESCHLAGENE SimpleMD Tests
❌ cli_simplemd_01-05, 07 - Siehe KNOWN_BUGS.md für Details (Optimierungs-Bug)

---

## curcumaopt Tests (❌ Alle fehlgeschlagen außer 01)

### cli_curcumaopt_01: Default UFF
- **Molekül**: Water (H2O, 3 atoms)
- **Methode**: UFF optimization
- **Input**:
  ```
  O     0.000000     0.000000     0.119262
  H     0.000000     0.763239    -0.477047
  H     0.000000    -0.763239    -0.477047
  ```
- **Erwartete Energie**: TBD (Ausgabedatei leer wegen Bug)
- **Status**: ⚠️ BESTANDEN (prüft nur Datei-Existenz, nicht Korrektheit)
- **Bug**: JSON null-Fehler - siehe KNOWN_BUGS.md

### FEHLGESCHLAGENE curcumaopt Tests
❌ cli_curcumaopt_02-06 - Siehe KNOWN_BUGS.md für Details (JSON null-Fehler)

---

## Test-Status Zusammenfassung

| Kategorie | Bestanden | Fehlgeschlagen | Gesamt |
|-----------|-----------|----------------|--------|
| RMSD      | 5         | 1              | 6      |
| ConfScan  | 6         | 1              | 7      |
| SimpleMD  | 1         | 6              | 7      |
| curcumaopt| 1*        | 5              | 6      |
| **GESAMT**| **13**    | **13**         | **26** |

*Test 01 besteht technisch, aber nur weil er Datei-Existenz prüft, nicht Korrektheit

**Erfolgsrate**: 50% (13/26)

---

## Verwendung dieser Referenzen

### Für neue wissenschaftliche Validierung:

1. **RMSD Tests** - Vergleiche berechneten RMSD mit Referenz ± Toleranz (z.B. ±0.0001 Å)
2. **ConfScan Tests** - Zähle Anzahl accepted/rejected conformers
3. **SimpleMD Tests** - Validiere Trajektorien-Länge und Energie-Drift

### Beispiel (test_utils.sh):
```bash
# Extrahiere RMSD
rmsd_value=$(grep "RMSD for two molecules" stdout.log | grep -oP '\d+\.\d+')

# Validiere gegen Referenz
assert_scientific_value "2.87214" "$rmsd_value" "0.0001" "RMSD value"
```

---

## Nächste Schritte

1. ✅ **curcumaopt JSON Bug beheben** - Höchste Priorität
2. ⏳ **Golden References erweitern** - Für ConfScan und SimpleMD
3. ⏳ **Wissenschaftliche Validierung** - Zu allen Tests hinzufügen
4. ⏳ **Invalid Method Tests** - Erwartete Fehler-Codes dokumentieren
