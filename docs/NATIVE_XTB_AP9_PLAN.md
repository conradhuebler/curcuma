# AP 9 — DFT-D4 Dispersion und D4-Gradient

**Status:** Offen — nach AP6 (korrekte Energien ohne Dispersion Voraussetzung für sinnvollen D4-Vergleich)
**Erstellt:** 2026-04-26
**Vorbedingung:** AP6 abgeschlossen; `dftd4interface.cpp` vorhanden (externe D4-Bibliothek)

---

## Ziel

GFN2-xTB enthält D4-Dispersion als Pflichtterm. Im nativen Code ist dieser Beitrag als Stub implementiert:

```cpp
// xtb_native.cpp (aktuell):
m_E_dispersion = 0.0;  // TODO: DFT-D4 interface
```

Die D4-Bibliothek (`dftd4`) ist bereits in Curcuma integriert (`dftd4interface.cpp`). AP9 drahtet sie in die native GFN2/GFN1-Berechnung ein.

**Energetische Bedeutung** (Beispiel H₂O·H₂O Dimer):
- D4-Dispersionsenergie: ~−5 kJ/mol (−2 mEh)
- Für größere Moleküle (C₆H₆, Peptide): −10–50 mEh

---

## Aufgabenliste

### 9.1 — D4Interface in `xtb_native.cpp` einbinden

Die `dftd4interface.cpp` bietet eine C++-Schnittstelle zu `dftd4`. Einbindung in `Calculation()`:

```cpp
// xtb_native.cpp: nach SCF-Konvergenz
#ifdef USE_D4
    DFTD4Interface d4;
    d4.setMolecule(m_atoms, m_geometry);
    d4.setMethod("gfn2");  // GFN2-D4 Parameter
    m_E_dispersion = d4.CalculateEnergy();
    if (gradient) {
        Matrix d4_grad = d4.getGradient();
        m_gradient += d4_grad;  // Additiv zum analytischen Gradient
    }
#else
    m_E_dispersion = 0.0;
    if (CurcumaLogger::get_verbosity() >= 1)
        CurcumaLogger::warn("D4 dispersion not available (USE_D4=OFF)");
#endif
```

### 9.2 — GFN2-D4 vs. GFN1-D3

- **GFN2**: D4 Dispersion (Ladungs-abhängig, höhere Genauigkeit)
- **GFN1**: D3(BJ) Dispersion (Ladungs-unabhängig)

`dftd3interface.cpp` ist ebenfalls vorhanden → GFN1 bekommt D3, GFN2 bekommt D4.

```cpp
if (m_method == MethodType::GFN2) {
    // D4
} else if (m_method == MethodType::GFN1) {
    // D3(BJ)
}
```

### 9.3 — GFN2-D4 Parameter

Die D4-Parameter für GFN2 sind in TBLite definiert (`tblite/xtb/gfn2.f90:get_dispersion`):
```fortran
s6 = 1.0_wp, s8 = 2.7_wp, a1 = 0.52_wp, a2 = 5.0_wp
```

Diese müssen korrekt an `DFTD4Interface` übergeben werden.

### 9.4 — D4-Gradient-Integration

Der D4-Gradient (`d4.getGradient()`, nat×3, Eh/Å) wird additiv zum analytischen Gradient addiert:

```cpp
// In calculateGradient():
if (m_E_dispersion != 0.0) {
    m_gradient += m_d4_gradient;  // nach AP9.1 gesetzt
}
```

Oder: D4 liefert Gradient zusammen mit Energie in einem Aufruf.

### 9.5 — Validierung

**Energievergleich** (GFN2 mit D4, gegen TBLite):

| Molekül | E_native (kein D4) | E_native (mit D4) | E_TBLite | Verbesserung |
|---------|-------------------|-------------------|----------|--------------|
| H₂O | — | — | — | — |
| C₆H₆ | — | — | — | — |

**FD-Gradient** nach AP5b: muss weiterhin passen (D4-Gradient ist analytisch von d4-Bibliothek).

---

## Akzeptanzkriterien

- [ ] GFN2-Energie für C₆H₆ liegt nach D4-Integration näher an TBLite als vorher
- [ ] D4-Gradient konsistent mit FD-Gradient (< 5e-4 Eh/Å)
- [ ] `#ifndef USE_D4` Guards vorhanden: Code kompiliert auch ohne D4
- [ ] GFN1 mit D3: analog funktional

## Fortschritt

| Datum | Aufgabe | Status | Notizen |
|-------|---------|--------|---------|
| 2026-04-26 | 9.1 D4-Einbindung | Offen | dftd4interface.cpp vorhanden |
| 2026-04-26 | 9.2 GFN1-D3 | Offen | dftd3interface.cpp vorhanden |
| 2026-04-26 | 9.3 Parameter | Offen | aus TBLite gfn2.f90 |
| 2026-04-26 | 9.4 Gradient | Offen | — |
| 2026-04-26 | 9.5 Validierung | Offen | — |

## Schwierigkeiten / Blocker

- D4 ist Ladungs-abhängig: braucht Mulliken-Ladungen aus konvergiertem SCF — Reihenfolge in `Calculation()` muss stimmen (SCF zuerst, dann D4 mit Ladungen)
- Falls D4-Bibliothek nicht angebunden: Energie bleibt 0.0 und Warnung (akzeptabler Fallback)

---

## Referenzen

- Curcuma: `src/core/energy_calculators/qm_methods/dftd4interface.cpp/h`
- Curcuma: `src/core/energy_calculators/qm_methods/dftd3interface.cpp/h`
- TBLite: `external/tblite/src/tblite/xtb/gfn2.f90` (D4-Parameter: s6, s8, a1, a2)
- Literatur: Caldeweyher et al., J. Chem. Phys. 150, 154122 (2019) — D4
