# WP-S1: Static-Mode — eingefrorene CN/Charges für MD

**Kategorie**: Klein–Mittel
**Aufwand**: ~1 Tag (~150 LoC)
**Wirkung**: Polymer-MD 1000 Steps: ICX 130 s → 75 s (≈42 % schneller, schlägt xtb-Fortran)
**Abhängigkeiten**: keine — orthogonal zu Cell-List (WP-S3)
**Status**: 🤖 Geplant
**Approximations-Profil**: liegt im Rauschen der GFN-FF-eigenen Toleranz (~mEh/Term);
für charge-transfer/Reaktionsdynamik ungeeignet

---

## Motivation

Pro MD-Schritt werden derzeit teure Schritte wiederholt, die sich bei kleinen
Geometrie-Fluktuationen (NVE/NVT um equilibrium) kaum ändern:

| Term | Kosten ICX | Drift bei MD um equilibrium |
|------|------------|------------------------------|
| CN-Berechnung (`erf`-Schleife) | ~5 ms | <1 % pro Step |
| `dcn` (CN-Gradient-Map) | ~10 ms | folgt CN |
| EEQ Phase-2 (PCG + dgam-Refinement) | ~30 ms | 0.001–0.01 e/Atom |
| D4 Gaussian-Weights + `dc6dcn` | ~10 ms | folgt CN |
| **Summe** | **~55 ms/Step** | |

Polymer-Baseline 130 s ICX → mit Static-Mode 75 s.

Annahme: bei `|Δr_i| < 0.1 Bohr` ändern sich CN/Charges um weniger, als die
GFN-FF-Methode selbst an absoluten Toleranzen mitbringt (~mEh/Term, je nach
System 0.01–0.5 e Charge-Unsicherheit gegenüber DFT). Static-Mode bewegt
sich also im methoden-immanenten Rauschen.

Klare Ausschluss-Bereiche bleiben: charge-transfer-Reaktionen, sp²/sp³-Übergänge,
große Konformations-Sprünge — da werden CN/Charge-Änderungen makroskopisch.

---

## Drei orthogonale Flags

```cpp
// gfnff.h, im BEGIN_PARAMETER_DEFINITION(gfnff)-Block
PARAM(static_charges, Bool, false,
      "Skip Phase-2 EEQ refinement after first call; reuse Phase-1 topology_charges. "
      "Saves ~30 ms/step. Accuracy: 0.001-0.01 e/atom drift around equilibrium. "
      "INVALID for charge-transfer reactions, ionic dynamics, or geometry "
      "deviations >0.5 Bohr from initial structure.",
      "Performance", {})

PARAM(static_cn, Bool, false,
      "Skip CN/dcn/D4-Gaussian-weight recomputation after first call; reuse init values. "
      "Saves ~25 ms/step (CN 5ms + dcn 10ms + D4 weights 10ms). "
      "INVALID for bond-bending across hybridization (sp2<->sp3) or large conformational changes.",
      "Performance", {})

PARAM(static_all, Bool, false,
      "Shorthand: enables static_charges AND static_cn. "
      "Effectively a classical fixed-charge MM mode with GFN-FF topology+parameters. "
      "Use only for stable production NVT/NPT in equilibrium regime.",
      "Performance", {})
```

`static_all=true` setzt intern beide Einzel-Flags. Aliases optional.

---

## Implementierungspunkte

### 1. Member-Cache in `GFNFF` (`gfnff.h`)

Bereits vorhanden: `m_last_cn`, `topology_charges` (in `GFNFFTopology`).
Neu zu ergänzen:

```cpp
// gfnff.h, im GFNFF-Klassenrumpf nach m_last_cn
Vector m_init_cn;                    // CN beim ersten Calculation()-Call
Matrix m_init_dcn;                   // CN-Gradient initial (für static_cn-Mode)
Vector m_init_d4_weights;            // D4 Gaussian-Weights je Atom×Reference
Vector m_init_d4_dc6dcn;             // dC6/dCN je Pair, initial
bool   m_static_state_captured = false;

// Flags (gespiegelt aus ConfigManager)
bool   m_static_charges = false;
bool   m_static_cn      = false;
```

### 2. Capture-Phase in `gfnff_method.cpp::Calculation()`

Vor dem üblichen CN-Recompute:

```cpp
if ((m_static_cn || m_static_charges) && !m_static_state_captured) {
    // Erst einmal regulär berechnen
    calculateCoordinationNumbers();          // füllt m_last_cn, m_dcn
    runEEQPhase1AndPhase2();                 // füllt m_eeq_charges
    precomputeD4WeightsAndDc6dcn();          // füllt D4-State
    m_init_cn       = m_last_cn;
    m_init_dcn      = m_dcn;
    m_init_d4_weights = m_d4_gw;             // exakte Member-Namen prüfen
    m_init_d4_dc6dcn  = m_d4_dc6dcn;
    m_static_state_captured = true;
}
```

### 3. Branch-Logik in `Calculation()`

```cpp
if (m_static_cn) {
    m_last_cn = m_init_cn;
    m_dcn     = m_init_dcn;                   // gradient-chain bleibt korrekt
    // D4 weights/dc6dcn werden aus init wieder eingespielt
} else {
    calculateCoordinationNumbers();
}

if (m_static_charges) {
    m_eeq_charges = topology_charges;         // Phase-1 reuse
    // Phase-2 PCG skip
} else {
    runEEQPhase2();
}
```

### 4. Gradient-Chain-Skips

- `forcefieldthread.cpp::CalculateGFNFFCoulombContribution` (TERM 1b): wenn
  `m_static_cn`, dann `dcn=0` → Term 1b liefert exakt null. **Kein Code-Pfad
  ändern** — `m_dcn` zeigt auf init-Werte (oder explizit `setZero()` für
  noch mehr Konsistenz).
- D4-Dispersion: wenn `m_static_cn`, dann `dC6/dCN`-Kette = 0, weil `dcn=0`.
  Energie-Beitrag bleibt korrekt, Gradient verliert chain-Term — bei kleinen
  Auslenkungen vernachlässigbar.

### 5. Init-Warnung und Drift-Monitor

```cpp
// gfnff_method.cpp ctor / Initialise(), nach Flag-Read
if (m_static_charges || m_static_cn) {
    CurcumaLogger::warn(fmt::format(
        "GFN-FF Static-Mode active: charges={}, cn={}. "
        "Energy gradient may deviate up to ~1 mEh/Bohr for large displacements. "
        "Reset advised every 500-1000 MD steps for production runs.",
        m_static_charges, m_static_cn));
}
```

Optionaler Drift-Monitor (separat steuerbar via `static_drift_check_freq`,
Default 0 = aus): vergleicht `(geometry - m_init_geometry).norm()` periodisch
gegen Schwellwert (0.5 Bohr) und warnt.

---

## Tests (`test_cases/cli/gfnff/`)

Fünf neue CLI-Tests:

1. `15_static_charges_caffeine` — single point, `static_charges=true`, vergleicht E gegen full
2. `16_static_cn_caffeine` — single point, `static_cn=true`
3. `17_static_all_caffeine` — single point, `static_all=true`
4. `18_static_all_md_polymer` — 100-step NVE, Energie-Drift <2 mEh/100 Steps
5. `19_static_all_md_acetic_dimer` — 100-step NVE, prüft H-Brücken bleiben stabil

Toleranzen:
- Single point E: |ΔE| < 0.5 mEh gegenüber full mode
- MD drift over 100 steps: < 2 mEh

---

## Performance-Acceptance

Polymer 1000 NVE-Steps ICX, Vergleich:

| Modus | Erwartet | Akzeptanz |
|-------|----------|-----------|
| full dynamic (baseline) | 130 s | — |
| `static_charges` allein | ≤ 105 s | mind. 18 % schneller |
| `static_cn` allein | ≤ 110 s | mind. 13 % schneller |
| `static_all` | ≤ 80 s | mind. 38 % schneller, **schlägt xtb-Fortran (~100 s)** |

---

## Out of Scope

- Adaptiver Re-Capture (z.B. bei Drift > Threshold automatisch refreshen): WP-S4
- GPU-Pfad-Anpassung: GPU recomputed CN/charges auf-device; Static-Mode initial
  nur CPU-relevant. GPU-Static-Mode kann eigenes Folge-WP werden.

---

## Acceptance-Checkliste

- [ ] 3 PARAMs in `gfnff.h` mit aussagekräftiger Help-Text-Warnung
- [ ] Member-Cache + Capture-Logik in `gfnff_method.cpp`
- [ ] Branch-Skips in `Calculation()` und (implizit via `m_dcn=0`) in
      `forcefieldthread.cpp` validiert
- [ ] Init-Warnung via CurcumaLogger
- [ ] 5 CLI-Tests grün
- [ ] Polymer-Benchmark erreicht Akzeptanzwerte
- [ ] `docs/GFNFF_STATUS.md` aktualisiert (neuer Static-Mode-Abschnitt mit
      Validitäts-Disclaimer)
