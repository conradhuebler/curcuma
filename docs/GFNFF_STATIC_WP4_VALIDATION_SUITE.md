# WP-S4: Cross-System Validation-Suite + adaptiver Re-Capture

**Kategorie**: Mittel
**Aufwand**: ~1.5 Tage (Validation-Skripte + adaptiver Re-Capture)
**Wirkung**: Sicherheitsnetz gegen wissenschaftlich falsche Static-Mode-Anwendung
**Abhängigkeiten**: WP-S1 (Static-Mode), WP-S2 (Diagnostics), WP-S3 (Cutoff)
**Status**: 🤖 Geplant
**Priorität**: erst nach WP-S1/S2/S3

---

## Motivation

Static-Mode (WP-S1) und Auto-Cutoff (WP-S3) sind beide **physikalisch
approximativ**. Ohne systematische Validation ist nicht klar, ab welcher
Geometrie-Auslenkung oder welchem System-Typ die Approximation versagt.

Ziel: Drei Liefer-Komponenten:
1. **Validation-Suite** — vier Test-Systeme × drei Static-Modi × zwei Cutoff-Modi
2. **Adaptiver Re-Capture** — Drift-Detector mit automatischem Refresh
3. **Forensik-Report-Generator** — extrahiert aus WP-S2-JSONL-Logs

---

## Komponente 1: Validation-Suite

### Test-Matrix

4 Systeme × {full, static_charges, static_cn, static_all} × {cutoff=0, cutoff=30}
= 32 Konfigurationen, jeweils 1000-Step NVE.

| System | Atome | Charakter | Erwartete Static-Mode-Eignung |
|--------|-------|-----------|--------------------------------|
| caffeine | 24 | rigid heterocyclic | ✅ alle Modi OK |
| triose | 66 | polare OH-reiche Kette | ⚠️ `static_charges` riskant (H-bond switching) |
| acetic_acid_dimer | 16 | dynamische H-Brücken | ❌ `static_all` versagt (HB rearrangement) |
| polymer (HDPE) | 472 | apolar, MD-Hauptanwendung | ✅ Bestes Static-Mode-Profil |

### Metriken pro Konfiguration

Aus `<run>.diag.jsonl` (WP-S2) extrahiert:

| Metrik | Schwellwert OK | Schwellwert ⚠️ | Schwellwert ❌ |
|--------|----------------|----------------|----------------|
| Total-E-Drift / 1000 Steps | <2 mEh | <10 mEh | ≥10 mEh |
| max\|Δq_i\| Step 0 vs 1000 | <0.01 e | <0.05 e | ≥0.05 e |
| max\|ΔCN_i\| Step 0 vs 1000 | <0.1 | <0.5 | ≥0.5 |
| RMSD vs full-mode-Lauf | <0.1 Å | <0.3 Å | ≥0.3 Å |
| HB count change | 0 | ±1 | ±2+ |

### Output

`docs/validation/static_mode_matrix.md` mit ampelfarbiger Tabelle —
32 Zellen, jeweils farbcodiert per Metrik-Aggregation.

### Implementierung

`scripts/static_mode_validation.py`:
- Iteriert Konfigurationen, ruft `curcuma -simplemd` mit den Flag-Kombis
- Parst pro Lauf die `.diag.jsonl`, extrahiert Metriken
- Schreibt Matrix-Markdown raus

CTest-Integration: `ctest -R static_validation` (nightly, nicht in default-Test-Set).

---

## Komponente 2: Adaptiver Re-Capture

Static-Mode mit hartem Einfrieren bricht für lange MDs zwangsläufig. Lösung:
periodischer Drift-Check + automatischer Re-Capture.

### Logik

```cpp
// gfnff_method.cpp::Calculation()
if (m_static_state_captured && m_static_drift_check_freq > 0
    && m_step_count % m_static_drift_check_freq == 0) {

    double geom_drift = (current_geometry - m_init_geometry).norm() / sqrt(N);

    if (geom_drift > m_static_drift_threshold) {
        CurcumaLogger::info(fmt::format(
            "Static-mode drift {:.3f} Bohr exceeds threshold {:.3f} — recapturing CN/charges",
            geom_drift, m_static_drift_threshold));
        m_static_state_captured = false;  // Forciert Re-Capture im nächsten Block
        m_recapture_count++;
    }
}
```

### Neue PARAMs

```cpp
PARAM(static_drift_check_freq, Int, 0,
      "Frequency (steps) to check geometry drift in static-mode. 0=disabled. "
      "When |geom - init_geom|_rms exceeds static_drift_threshold, "
      "CN/charges are recaptured. Recommended: 100 for production MD.",
      "Performance", {})

PARAM(static_drift_threshold, Double, 0.3,
      "RMS geometry displacement (Bohr) above which static-mode triggers "
      "re-capture of CN/charges. Default 0.3 Bohr ~ 0.16 Å.",
      "Performance", {})
```

### Performance-Erwartung

Polymer-MD mit `static_drift_check_freq=100`, `threshold=0.3`:
- Re-Capture vermutlich alle 200–500 Steps für typische 300 K NVT
- Amortisierte Kosten: ~55 ms/300 Steps ≈ 0.2 ms/Step zusätzlich
- Damit immer noch ~50 ms/Step Save → polymer 1000 NVE ≤ 82 s

---

## Komponente 3: Forensik-Report-Generator

`scripts/diag_report.py` — wandelt JSONL-Log in HTML-Report:
- Multi-Panel-Plot via matplotlib oder gnuplot
  - Total-Energy + Komponenten über Zeit
  - max|q_i| und max|CN_i| Drift
  - HB/XB count timeline
  - Gradient-norm-Histogramm
- Heuristik-Bewertung: "Run looks stable" / "Charge drift suspicious" / "Energy
  drift exceeds 10 mEh — check static-mode validity"
- Export: `<run>.report.html`

CLI:
```bash
python scripts/diag_report.py run.diag.jsonl --output run.report.html
```

---

## Tests

`test_cases/cli/static_validation/`:
1. `01_drift_recapture_polymer` — 500-step polymer mit `freq=50`, prüft mind. 1 Re-Capture
2. `02_validation_matrix_smoke` — 32 Konfigurationen × 100 Steps (smoke, nicht
   full 1000); Matrix-Tabelle wird generiert
3. `03_diag_report_html` — `diag_report.py` lauft ohne Crash, HTML enthält
   erwartete Sections

---

## Acceptance-Checkliste

- [ ] `scripts/static_mode_validation.py` läuft 32-Konfig-Matrix
- [ ] `docs/validation/static_mode_matrix.md` mit aktueller Ampel-Tabelle
- [ ] 2 neue PARAMs (`static_drift_check_freq`, `static_drift_threshold`)
- [ ] Re-Capture-Logik + CurcumaLogger-Info bei jedem Refresh
- [ ] Counter `m_recapture_count` im Final-Summary von simplemd
- [ ] `scripts/diag_report.py` HTML-Generator committed
- [ ] 3 neue CLI-Tests grün
- [ ] Nightly-CTest-Label `static_validation` registriert (excluded from default)
- [ ] `docs/GFNFF_STATUS.md` mit Verweis auf Validation-Matrix
