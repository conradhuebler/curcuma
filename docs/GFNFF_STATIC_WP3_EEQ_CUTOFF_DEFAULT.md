# WP-S3: EEQ-Coulomb-Cutoff — Default-Aktivierung mit Validation

**Kategorie**: Klein
**Aufwand**: ~0.5 Tag Code + ~1 Tag Validation
**Wirkung**: Coulomb-Matrix-Aufbau 17 ms → 5 ms (~12 ms/Step gespart), polymer 1000 NVE ~12 s schneller
**Abhängigkeiten**: keine — orthogonal zu WP-S1/S2
**Status**: 🤖 Geplant — Cutoff-Pfad existiert bereits, nur Default wechseln
**Risiko**: 🟠 hoch — verletzt Hellmann-Feynman gegenüber voller Coulomb-Energie

---

## Motivation

`eeq_distance_cutoff` ist heute auf `0.0` defaultet (matches Fortran `goed_gfnff`),
weil ein non-zero-Cutoff die Coulomb-Matrix sparsifiziert und damit Energie-
Erhaltung in MD verschlechtern kann (siehe Help-Text `eeq_solver.h:965-966`).

In Praxis ist der Effekt bei `cutoff ≥ 30 Bohr ≈ 16 Å` für neutral-polare
Systeme jedoch klein (Coulomb-Term fällt mit 1/r ab, aber abgeschirmt durch
EEQ-Gleichungsstruktur). Wenn ein Default-Cutoff für **nicht-ionische** Systeme
validierbar wäre, ergäbe das ~12 ms/Step polymer-Save.

**Konservativer Vorschlag**: Default bleibt `0.0`, aber neuer Auto-Modus
`eeq_distance_cutoff=auto` aktiviert `30.0` für Systeme mit `nfrag==1` und
`max|q|<0.5 e` aus Phase-1.

---

## Implementierungspunkte

### 1. PARAM-Erweiterung (`eeq_solver.h:965`)

Behalte Default `0.0`, aber dokumentiere Auto-Modus:

```cpp
PARAM(eeq_distance_cutoff, Double, 0.0,
      "Distance cutoff in Bohr for Coulomb matrix sparsification. "
      "0 = no cutoff (matches Fortran goed_gfnff, full energy/gradient consistency). "
      "Set explicitly to 30.0 for neutral molecular systems with weak Coulomb tails — "
      "saves ~12 ms/step for polymer-sized systems but violates strict Hellmann-Feynman "
      "vs. full Coulomb. INVALID for ionic systems, electrolytes, or charged fragments.",
      "Advanced", {})

PARAM(eeq_distance_cutoff_auto, Bool, false,
      "Auto-enable eeq_distance_cutoff=30 Bohr when Phase-1 yields nfrag==1 "
      "and max|q| < 0.5 e (heuristic for neutral molecular systems). "
      "Falls back to 0.0 (full Coulomb) for ionic or multi-fragment systems.",
      "Advanced", {})
```

### 2. Auto-Detection-Logik (`gfnff_method.cpp` nach Phase-1)

```cpp
double resolveEEQCutoff(double user_cutoff, bool auto_flag,
                       int nfrag, double max_abs_charge) const {
    if (auto_flag && user_cutoff == 0.0) {
        if (nfrag == 1 && max_abs_charge < 0.5) {
            CurcumaLogger::info(
                "EEQ cutoff auto-enabled: 30.0 Bohr (neutral mono-fragment system)");
            return 30.0;
        }
        CurcumaLogger::info(fmt::format(
            "EEQ cutoff auto: keeping 0.0 (nfrag={}, max|q|={:.3f} e — not neutral/single)",
            nfrag, max_abs_charge));
    }
    return user_cutoff;
}
```

Aufruf direkt nach Phase-1, vor Phase-2-Matrix-Aufbau.

### 3. Validation-Suite (`test_cases/validation/eeq_cutoff/`)

5 Systeme, jeweils mit `cutoff=0` vs `cutoff=30`:

| System | Erwartet OK | Anmerkung |
|--------|-------------|-----------|
| caffeine | ✅ ja | neutral, kompakt |
| triose | ✅ ja | neutral, polare OHs |
| acetic_acid_dimer | ✅ ja | H-bonded, neutral |
| polymer (HDPE-Modell) | ✅ ja | hauptanwendung |
| NaCl-cluster (4 Ionen) | ❌ nein — `cutoff=0` muss bleiben | Validiert Auto-Skip |

**Akzeptanz-Kriterien**:
- |ΔE_total| < 0.5 mEh für 4 neutrale Systeme
- |Δq_max| < 0.005 e/Atom für 4 neutrale Systeme
- Auto-Detection schaltet NaCl-cluster nicht ein

### 4. MD-Stabilität (`test_cases/cli/simplemd/12_eeq_cutoff_nve`)

100-Step NVE polymer mit `cutoff=30`:
- Energie-Drift < 1 mEh / 100 Steps
- max|Δq| zwischen Step 0 und 100 < 0.01 e

---

## Performance-Acceptance

| System | Baseline (cutoff=0) | Mit cutoff=30 | Save |
|--------|---------------------|----------------|------|
| polymer 1000 NVE ICX | 130 s | ≤ 118 s | ≥ 12 s, ≥ 9 % |
| triose 100 NVE ICX | — | — | proportional |

Kombiniert mit WP-S1 (Static-Mode): polymer ~63 s = **40 % schneller als xtb-Fortran**.

---

## Risiken & Mitigation

1. **Hellmann-Feynman-Verletzung im Gradient**: gradient von abgeschnittener
   Coulomb-Matrix ist nicht der exakte Energie-Gradient. Mitigation:
   Validation-Test `12_eeq_cutoff_nve` muss Energie-Erhaltung explizit prüfen.
2. **Ionische Systeme falsch detektiert**: `max|q|<0.5` ist Heuristik.
   Mitigation: Schwellwert konservativ; Doc empfiehlt explizites `cutoff=0`
   für Salzlösungen.
3. **Default nicht ändern**: WP setzt explizit nur Auto-Flag dazu, keine
   bestehende Test-Regression möglich.

---

## Acceptance-Checkliste

- [ ] PARAM `eeq_distance_cutoff_auto` registriert
- [ ] Auto-Detection-Routine implementiert + CurcumaLogger-Info
- [ ] 5-System-Validation-Suite committed (incl. NaCl-Cluster als Negativ-Test)
- [ ] CLI-Test `12_eeq_cutoff_nve` grün (Energie-Drift < 1 mEh / 100 Steps)
- [ ] Polymer-Benchmark erreicht ≥ 9 % Speedup
- [ ] `docs/GFNFF_STATUS.md` mit "Recommended for production MD: `cutoff_auto=true`"
- [ ] WP-S1-Kompatibilität: kombiniert getestet
