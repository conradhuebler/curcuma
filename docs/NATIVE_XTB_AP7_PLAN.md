# AP 7 — d-Schalen-Unterstützung (S, P, Cl, Transition Metals)

**Status:** Offen — nach AP6 (stabile Energien für s/p-Elemente Voraussetzung)
**Erstellt:** 2026-04-26
**Vorbedingung:** AP5b + AP6 abgeschlossen (Gradient korrekt, Energie für H/C/N/O korrekt)

---

## Ziel

Alle d-Schalen-Beiträge sind in der nativen Implementierung mit `if (typeA < 0) continue` übersprungen. Das betrifft:
- **Elemente ab 2. Periode mit d-Orbitalen**: S (3d), P (3d), Cl (3d), Br (4d)  
- **Übergangsmetalle**: Cu, Zn, Fe, Ni, Co etc.
- **Praktische Bedeutung**: Viele biologische Systeme und anorganische Chemie

Nach AP7 können GFN2/GFN1-Rechnungen mit diesen Elementen korrekt durchgeführt werden.

---

## Betroffene Komponenten

| Komponente | Aktueller Stand | Änderung nötig |
|-----------|-----------------|----------------|
| `ao_to_type()` in `STO_CGTO.hpp` | gibt -1 für ang ≥ 2 | d-Typ-Mapping hinzufügen |
| `cgto_overlap` in `STO_CGTO.hpp` | d-Integrale fehlen | Obara-Saika für l=2 ergänzen |
| `cgto_multipole_ints` | `if (typeA < 0) continue` | d-Multipol-Integrale hinzufügen |
| `cgto_multipole_grad` (AP5b) | `if (typeA < 0) continue` | d-Gradientintegrale hinzufügen |
| `xtb_h0.cpp` | d-Polynome (`shpoly`) | `poly_d[]` Parameter aus TBLite |
| `gfn2_params.hpp` / `gfn1_params.hpp` | H,C,N,O vollständig | d-Elemente extrahieren |

---

## Aufgabenliste

### 7.1 — Parameter-Extraktion für d-Elemente

```bash
# scripts/extract_xtb_params.py erweitern:
# Alle Elemente mit d-Schalen (Z=11-18, 21-36, etc.)
python3 scripts/extract_xtb_params.py --elements all --method gfn2
```

`gfn2_params.hpp` muss vollständige Tabellen für alle Elemente bis Z=86 enthalten (TBLite deckt Z=1–86 ab).

### 7.2 — `cgto_overlap` für d-Funktionen

Obara-Saika-Integrale bis l=2 (d-Funktionen). TBLite-Referenz: `tblite/integral/overlap.f90`.

Die 5 kartesischen d-Orbitale: dxy, dxz, dyz, dx²-y², dz² (reeller sphärischer Harmonischer Basis oder kartesisch — konsistent mit TBLite wählen).

**Komplexität hoch**: Obara-Saika für l=2 hat viele Terme.

### 7.3 — d-Typ-Mapping in `ao_to_type()`

```cpp
// STO_CGTO.hpp: ao_to_type()
// Aktuell:
case 2: return -1;  // d — skip
// Nach AP7:
case 2: return 2;   // d — type 2
```

Alle Funktionen, die `-1` als Abbruchbedingung nutzen, müssen entsprechend angepasst werden.

### 7.4 — d-Multipol-Integrale (optional, abhängig von AP5b)

`cgto_multipole_ints` und `cgto_multipole_grad` für d-Schalen-AOs.
TBLite-Referenz: `tblite/integral/multipole.f90:multipole_cgto` für `ang=2`.

### 7.5 — H0-Parameter für d-Elemente

`shpoly`-Polynom-Parameter für d-Schalen aus TBLite extrahieren (analog zu s/p in `gfn2_params.hpp`).

### 7.6 — Validierung gegen TBLite

Testmoleküle mit d-Elementen:
- H₂S (Schwefel, Z=16)  
- HCl (Chlor, Z=17)
- ZnH₂ (Zink, Z=30) — wenn Parameter vorhanden

```bash
./curcuma -sp H2S.xyz -method gfn2
./curcuma -sp H2S.xyz -method ipea1   # TBLite-Referenz
```

---

## Akzeptanzkriterien

- [ ] H₂S mit GFN2: Energie endlich und konvergiert
- [ ] H₂S: |ΔE_native − E_TBLite| < 5e-3 Eh
- [ ] HCl: Energie endlich und konvergiert
- [ ] Keine Regressions für H₂O, CH₄, NH₃ (s/p-Elemente)
- [ ] Gradiententest für H₂S passiert (nach AP5b)

## Fortschritt

| Datum | Aufgabe | Status | Notizen |
|-------|---------|--------|---------|
| 2026-04-26 | 7.1 Parameter | Offen | — |
| 2026-04-26 | 7.2 cgto_overlap d | Offen | Obara-Saika l=2, hohe Komplexität |
| 2026-04-26 | 7.3 ao_to_type | Offen | einfach |
| 2026-04-26 | 7.4 d-Multipole | Offen | abhängig von AP5b-Code |
| 2026-04-26 | 7.5 H0 d-Poly | Offen | — |
| 2026-04-26 | 7.6 Validierung | Offen | — |

## Schwierigkeiten / Blocker

- Obara-Saika für l=2 ist aufwendig — alternativ: numerische Integration als Fallback mit analytischer Grenze für l≤1
- Übergangsmetalle (d-Blocks, Z=21–30) haben komplexere Elektronenstruktur → möglicherweise separate Validierung nötig

---

## Referenzen

- TBLite: `external/tblite/src/tblite/integral/overlap.f90` (Obara-Saika bis l=2)
- TBLite: `external/tblite/src/tblite/xtb/gfn2.f90` (d-Schalen-Parameter)
- Curcuma: `src/core/energy_calculators/qm_methods/STO_CGTO.hpp` (cgto_overlap, ao_to_type)
- Parameter: `scripts/extract_xtb_params.py` (erweitern für alle Elemente)
