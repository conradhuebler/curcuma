# AP 8 — DIIS-Mischer und SCF-Robustheit

**Status:** Offen — kann parallel zu AP6/AP7 entwickelt werden
**Erstellt:** 2026-04-26
**Vorbedingung:** AP1–4 (lauffähiges SCF); AP6 empfohlen (bekannte Energiefehler vor DIIS-Tuning)

---

## Ziel

Der native SCF-Loop nutzt ausschließlich lineares Dämpfen (`damping = 0.4`). Das führt bei:
- **Stark polarisierten Systemen** (z.B. Aminosäuren mit Carboxylat)
- **Konjugierten π-Systemen** (z.B. Benzol, Naphthalin)  
- **Übergangsmetall-Komplexen** (nach AP7)

zu sehr langsamer Konvergenz (50+ Iterationen) oder gelegentlichem Nicht-Konvergieren.

`diis_accelerator.h` ist bereits im Codebase vorhanden aber nicht angebunden. AP8 drahtet DIIS ein und macht es robust.

---

## Aufgabenliste

### 8.1 — `diis_accelerator.h` prüfen und anpassen

```bash
grep -n "class DIIS\|apply\|push\|reset" src/core/energy_calculators/qm_methods/diis_accelerator.h
```

Erwartete Schnittstelle:
```cpp
class DIISAccelerator {
    void push(const Matrix& F, const Matrix& error);  // Fehlervektor = FPS - SPF
    Matrix extrapolate();                              // neue Fock-Matrix
    void reset();
};
```

Falls Schnittstelle nicht passt: minimal anpassen.

### 8.2 — DIIS in `xtb_scf.cpp` einbinden

Fehlervektor für DIIS in xTB: `e = FPS - SPF` (Commutator), Elemente im AO-Basis.

```cpp
// In SCF-Loop (xtb_scf.cpp):
if (use_diis && iter >= diis_start) {
    Matrix error = F * P * S - S * P * F;
    diis.push(F, error);
    if (iter >= diis_start + diis_vectors)
        F = diis.extrapolate();
}
```

DIIS-Start: ab Iteration 3 (erste Iterationen lineares Dämpfen als Initialisierung).
DIIS-Vektoren: 6–8 (Standardwert in xTB).

### 8.3 — Fallback auf lineares Dämpfen

Bei DIIS-Divergenz (Fehler wächst statt fällt) automatisch zurück zu linarem Dämpfen:
```cpp
if (diis_error > 10 * prev_error) {
    diis.reset();
    use_diis = false;
    CurcumaLogger::warn("DIIS diverged, switching to linear damping");
}
```

### 8.4 — Parameter via ConfigManager

```cpp
// In xtb_native.h / PARAM-Definitionen:
PARAM(use_diis, Bool, true, "Enable DIIS SCF accelerator")
PARAM(diis_vectors, Int, 6, "Number of DIIS vectors")
PARAM(diis_start, Int, 3, "SCF iteration to start DIIS")
```

### 8.5 — Validierung

**Konvergenz-Vergleich**:
| Molekül | Vorher (Iter) | Nachher (Iter) |
|---------|--------------|---------------|
| H₂O | — | — |
| C₆H₆ | — | — |
| Aspirin | — | — |

**Energiekonsistenz**: Energien nach DIIS müssen identisch zu linearem Dämpfen sein (innerhalb SCF-Toleranz).

---

## Akzeptanzkriterien

- [ ] Alle bestehenden Testmoleküle konvergieren weiterhin (Regression)
- [ ] C₆H₆ konvergiert mit < 20 Iterationen (vorher oft > 50)
- [ ] Energien nach DIIS stimmen mit linearem Dämpfen überein (< 1e-6 Eh Differenz)
- [ ] DIIS-Fallback auf lineares Dämpfen funktioniert

## Fortschritt

| Datum | Aufgabe | Status | Notizen |
|-------|---------|--------|---------|
| 2026-04-26 | 8.1 diis_accelerator.h | Offen | Datei existiert, Schnittstelle unklar |
| 2026-04-26 | 8.2 DIIS in xtb_scf.cpp | Offen | — |
| 2026-04-26 | 8.3 Fallback | Offen | — |
| 2026-04-26 | 8.4 Parameter | Offen | — |

## Schwierigkeiten / Blocker

- DIIS-Fehlervektor `FPS - SPF` muss im MO-Basis oder AO-Basis berechnet werden (TBLite nutzt MO-Basis) — konsistent mit `xtb_scf.cpp` wählen

---

## Referenzen

- Curcuma: `src/core/energy_calculators/qm_methods/diis_accelerator.h`
- Curcuma: `src/core/energy_calculators/qm_methods/xtb_scf.cpp` (SCF-Loop)
- TBLite: `external/tblite/src/tblite/scf/mixer.f90` (DIIS-Referenz)
- Literatur: Pulay 1980, CHEM. PHYS. LETTERS 73, 393
