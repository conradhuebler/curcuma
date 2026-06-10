# AP 6 — GFN2 Energiegenauigkeit: vat_extra-Debugging

**Status:** Offen — nach AP5b (vollständiger Gradient vor Energiedebugging sinnvoll)
**Erstellt:** 2026-04-26
**Vorbedingung:** AP5b abgeschlossen (FD-Gradient passt) — damit ist sicher, dass Fehler aus Energie, nicht Gradient kommt

---

## Ziel

Die native GFN2-Energie ist systematisch ~35–60 mEh gegenüber TBLite zu hoch.

**Bekannte Befunde:**
- `vat_extra`-Term in `addMultipolePotential()` weicht ~1.5e-4 Eh von TBLite ab
- Verdacht: `updown_to_magnet`-Funktion in der Quadrupol-Populationsanalyse (`xtb_scf.cpp`)
- TBLite-Dumps zeigen: H0-Matrix (nicht Fock-Matrix) ist Referenz für `dump.hamiltonian`
- Energiedekomposition (`getEnergyDecomposition()`) zeigt, welcher Term abweicht

---

## Diagnosestrategie

### 6.1 — Energiedekomposition gegen TBLite

TBLite-Referenz via `dump_tblite` oder `ipea1`-Methode (die TBLite aufruft):

```bash
./curcuma -sp H2O.xyz -method gfn2 -verbosity 2    # Native Dekomposition
./curcuma -sp H2O.xyz -method ipea1 -verbosity 2   # TBLite (als Referenz)
```

Erwartetes Ausgabeformat (via `getEnergyDecomposition()`):
```
Electronic:   -X.XXXX Eh
Coulomb ES2:  -X.XXXX Eh
ThirdOrder:   -X.XXXX Eh
Multipole:    -X.XXXX Eh
Repulsion:    +X.XXXX Eh
Total:        -X.XXXX Eh
```

Welcher Einzelterm stimmt nicht mit TBLite überein?

### 6.2 — `updown_to_magnet` in `xtb_scf.cpp`

TBLite trennt Quadrupol-Komponenten in „up-down" (Spherical-Harmonic-Basis) und konvertiert zu kartesischen Momenten via `updown_to_magnet`. Curcuma nutzt bereits kartesische Quadrupole — es ist unklar, ob diese Konvertierung korrekt in Curcumas kartesischer Darstellung ist.

**Test:** Vergleich von `m_wfn.qp_at` (kartesisch, 6 Komponenten, packed) gegen TBLite-Dump nach erster SCF-Iteration.

### 6.3 — `vat_extra` in `addMultipolePotential()`

```cpp
// xtb_multipole.cpp: addMultipolePotential()
// vat_extra = on-site dkernel + qkernel contribution to atomic potential
```

Vergleich der `pot.v_at`-Einträge nach `addMultipolePotential()` gegen TBLite-Dump:
- Nutze `dump_tblite_multipole` für TBLite-Referenz-Potenziale
- Identifiziere, ob Abweichung in `dkernel` oder `qkernel`

### 6.4 — Mulliken-Quadrupolmomente

GFN2 Mulliken-Quadrupole (`qp_at`, 6×nat) akkumulieren Fehler aus der Basis-Darstellung. Vergleich von `m_wfn.qp_at` nach konvergiertem SCF gegen TBLite-Dump.

---

## Aufgabenliste

| Sub-AP | Aufgabe | Datei | Aufwand |
|--------|---------|-------|---------|
| 6.1 | Energiedekomposition-Vergleich Native vs. TBLite (tabellarisch) | — | niedrig |
| 6.2 | `updown_to_magnet`-Analyse: Quadrupol-Koordinatenbasis | `xtb_scf.cpp` | mittel |
| 6.3 | `vat_extra` debuggen: dkernel/qkernel gegen TBLite-Dump | `xtb_multipole.cpp` | mittel |
| 6.4 | Mulliken-Quadrupole gegen TBLite vergleichen | `xtb_scf.cpp` | niedrig |
| 6.5 | Fix implementieren und FD-Test + Energie-Toleranzen prüfen | — | hoch |

---

## Akzeptanzkriterien

- [ ] GFN2-Energie für H₂O: |ΔE_native − E_TBLite| < 1e-3 Eh
- [ ] GFN2-Energie für H₂O, CH₄, NH₃, C₆H₆: alle < 1e-3 Eh
- [ ] `sqm_scf_*_gfn2` CTests (aktuell FAIL wegen Δε_max > 1e-4): bestehen nach Fix
- [ ] Kein Rückschritt bei GFN1-Energien

## Fortschritt

| Datum | Aufgabe | Status | Notizen |
|-------|---------|--------|---------|
| 2026-04-26 | 6.1 Dekomposition | Offen | — |
| 2026-04-26 | 6.2 updown_to_magnet | Offen | Hauptverdacht |
| 2026-04-26 | 6.3 vat_extra | Offen | ~1.5e-4 Eh bekannt |
| 2026-04-26 | 6.4 Mulliken-qp_at | Offen | — |

## Schwierigkeiten / Blocker

- Fehler kann in der Quadrupol-Populationsanalyse (Mulliken), in den Potentialen, oder in der Energieformel selbst liegen — erfordert systematischen SCF-Dump-Vergleich
- TBLite-Dumps (`dump_tblite_multipole`) müssen aktuell sein (USE_TBLITE nötig)

---

## Referenzen

- Memory: `project_gfn2_multipole_debug.md` — bekannte Befunde (vat_extra, updown_to_magnet)
- Curcuma: `src/core/energy_calculators/qm_methods/xtb_multipole.cpp` (vat_extra, addMultipolePotential)
- Curcuma: `src/core/energy_calculators/qm_methods/xtb_scf.cpp` (Mulliken, Quadrupole)
- TBLite: `external/tblite/src/tblite/scf/potential.f90` (v_at Berechnung)
- TBLite: `external/tblite/src/tblite/xtb/gfn2.f90` (on-site dkernel/qkernel Parameter)
- Dumps: `test_cases/sqm_reference/dump_tblite_multipole.cpp`
