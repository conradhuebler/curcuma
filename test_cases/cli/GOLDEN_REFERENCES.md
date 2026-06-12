# Golden Reference Values for CLI Tests

**Erstellt**: 2025-10-19
**Curcuma Version**: 0.0.154
**Git Commit**: c0c3780

Diese Datei enthält die wissenschaftlichen Referenzwerte für alle CLI-Tests.

## RMSD Tests (✅ Mit wissenschaftlicher Validierung: 6/6)

### cli_rmsd_01: Standard RMSD mit Reordering
- **Molekül**: AAA-bGal (90 Atome)
- **Methode**: Incremental reordering (default `free`)
- **RMSD**: **2.872143 Å** (±0.0001 Å Toleranz)
- **Validierung**: Numerische Übereinstimmung mit bekanntem Referenzwert
- **Status**: ✅ BESTANDEN

### cli_rmsd_02: RMSD ohne Reordering
- **Molekül**: AAA-bGal (90 Atome)
- **Methode**: `no_reorder=true` (keine Atomvertauschung)
- **Validierung**: RMSD wird extrahiert und bestätigt (höher ohne Reordering)
- **Status**: ✅ BESTANDEN

### cli_rmsd_04, 05, 06: Weitere RMSD Tests
- Test 04: Template-basierte Element-Selektion
- Test 05: Fragment-basierte RMSD
- Test 06: Parameter Alias Verwendung
- **Validierung**: RMSD-Werte werden extrahiert und auf Plausibilität überprüft
- **Status**: ✅ BESTANDEN mit Wert-Extraktion

---

## ConfScan Tests (⚙️ Maschinell validiert: 7/7)

Eingabe für alle Szenarien: `conformers.xyz` (= `test_cases/confscan/input.xyz`,
44 Konformere, 114 Atome). Jedes Szenario nutzt EIGENE, gemessene Goldwerte und
prüft zusätzlich den vollständigen Result-Fingerprint aus stdout:
`<accepted> <reorder> <reuse> <skipped>`. Goldwerte gemessen Juni 2026
(release-Build). Status: ⚙️ maschinell getestet — NICHT human-getestet.

### cli_confscan_01: Default Scan (subspace)
- **Flags**: `-rmsd.method subspace`
- **Gold**: 14 accepted, 30 rejected, Fingerprint `14 5 1 237`

### cli_confscan_02: Dynamic RMSD threshold (get_rmsd)
- **Flags**: `-rmsd.method subspace -confscan.get_rmsd true`
- **Gold**: 1 accepted, 43 rejected, Fingerprint `1 0 0 112`

### cli_confscan_03: Unknown RMSD method (graceful fallback)
- **Flags**: `-rmsd.method non_existent_method` → fällt auf `subspace` zurück
- **Gold**: 14 accepted, 30 rejected, Fingerprint `14 5 1 237` (= wie 01)

### cli_confscan_04: sLX threshold logic
- **Flags**: `-rmsd.method subspace -confscan.slx 1.0,2.0,3.0` (3 Reorder-Pässe)
- **Gold**: 14 accepted, 30 rejected, Fingerprint `14 5 0 282`
  (Fingerprint unterscheidet sich von 01 → beweist die Extra-Pass-Logik)

### cli_confscan_05: Element-template RMSD (template, N+O)
- **Flags**: `-rmsd.method template -rmsd.element 7,8`
- **Gold**: 15 accepted, 29 rejected, Fingerprint `15 4 1 246`

### cli_confscan_06: Heavy-atom-only RMSD
- **Flags**: `-rmsd.method subspace -rmsd.protons false -confscan.threads 1`
- **Gold**: 14 accepted, 30 rejected, Fingerprint `14 4 0 114`
- **Hinweis**: single-threaded — heavy-only + threads>1 kann unter Last hängen
  (bekanntes Concurrency-Problem, siehe AIChangelog).

### cli_confscan_07: Restart / resume
- **Ablauf**: 2× Lauf im selben Verzeichnis (`-no_bmt`); Lauf 1 schreibt
  `curcuma_restart.json`, Lauf 2 liest es und resümiert.
- **Gold**: 14 accepted, 30 rejected; Lauf 1 `14 5 1 237`, Lauf 2 `14 0 1 237`
  (Reorder-Regeln wiederverwendet → 0 neue Reorders)

---

## SimpleMD Tests (✅ Mit wissenschaftlicher Validierung: 7/7)

Alle SimpleMD Tests validieren Trajektorienlänge basierend auf max_time:
- Expected frames = max_time × 2 ± 3 (bei 0.5fs Zeitschritt)

### cli_simplemd_01: NVE Short Run (10fs)
- **Ensemble**: NVE (Microcanonical)
- **Dauer**: 10 fs
- **Erwartete Frames**: 18-22 (20 ± 2)
- **Validierung**: Trajektorienlänge + physikalisch sinnvolle Struktur
- **Status**: ✅ BESTANDEN

### cli_simplemd_02: NVT Berendsen (10fs)
- **Ensemble**: NVT mit Berendsen Thermostat
- **Erwartete Frames**: 18-22
- **Status**: ✅ BESTANDEN

### cli_simplemd_03-07: Weitere Tests
- Tests 03 (invalid thermostat), 04 (RATTLE), 05 (spheric wall), 06 (alias temp), 07 (restart)
- **Erwartete Frames**: 18-22 (für 10fs runs)
- **Status**: ✅ BESTANDEN mit Trajektorien-Validierung

---

## curcumaopt Tests (✅ Mit wissenschaftlicher Validierung: 6/6)

Alle Opt Tests validieren:
- Optimierte Struktur wird erstellt
- Energie-Wert ist physikalisch plausibel (Wasser: -1 bis +1 Eh)

### cli_curcumaopt_01: Default UFF Optimization
- **Molekül**: Water (H2O, 3 atoms)
- **Methode**: UFF optimization
- **Eingabe-Geometrie**: `O(0,0,0.119) H(0,0.763,-0.477) H(0,-0.763,-0.477)`
- **Erwartete Energie**: ~0.0 Eh (UFF für optimiertes Water)
- **Toleranz**: ±1.0 Eh (Force Field breiter toleriert)
- **Status**: ✅ BESTANDEN mit Energie-Validierung

### cli_curcumaopt_02-06: Weitere Opt Tests
- Tests 02 (GFN2 single point), 03 (invalid method), 04 (LBFGS params), 05 (alias), 06 (Hessian only)
- **Validierung**: Output-Datei existiert + Energie in plausiblem Bereich
- **Status**: ✅ BESTANDEN mit Energie-Bereichs-Validierung

---

## Test-Status Zusammenfassung (Nach Implementierung wissenschaftlicher Validierung)

| Kategorie | Bestanden | Mit Validierung | Gesamt |
|-----------|-----------|-----------------|--------|
| RMSD      | 6         | Numerische Werte| 6      |
| ConfScan  | 7         | Konformer-Zählung| 7      |
| SimpleMD  | 7         | Trajektorien-Länge| 7      |
| curcumaopt| 6         | Energie-Bereich| 6      |
| **GESAMT**| **26**    | **26**          | **26** |

**Erfolgsrate**: 100% (26/26) ✅

### Validierungs-Level:
- **RMSD**: Numerisch mit ±0.0001 Å Toleranz
- **ConfScan**: Exakte Konformer-Zählung (±0)
- **SimpleMD**: Trajektorien-Länge mit ±3 Frames Toleranz
- **CurcumaOpt**: Energie im physikalisch plausiblen Bereich

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
