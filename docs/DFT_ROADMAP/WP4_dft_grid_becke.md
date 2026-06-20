# WP4 — DFT-Quadratur-Grid (Euler-Maclaurin + Lebedev + Becke)

- **Status:** ⚠️ AI-generated / offen — ⚙️ machine-tested (nach Erfüllung)
- **Abhängigkeit:** WP3 (für Dichte zum Testen)
- **Validierung:** Elektronennummer + Gitterkonvergenz + Lebedev-Sphärenharmonischen-Test

## Ziel

Atomzentriertes Mehrzentren-Grid für Moleküle (xcDFTs einzentriger Fix); radial
Euler-Maclaurin, angular Lebedev-Laikov, Becke-Atompartition.

## Deliverables

- [ ] `dft_grid.hpp/.cpp`: `EulMac` + `Lebdev` portiert aus xcDFT `dft_grid.f`
      (zitiert, gekürzte Lebedev-Ordnungen 50/110/194/302), `read_grid`-Map
      (SG-0…SG-3 aus xcDFT `read_grid.f90`). **Neu:** Becke-Mehrzentren-Atompartition
      (Becke 1988) mit Atomradien-Skala; Grid-Punkte + Gewichte pro Atom.
- [ ] AO-Auswertung auf Grid (`AO_values_grid.f90`-Port) + kartesische AO-Gradienten
      (für WP8). `density`/`gradient_density` (xcDFT-Port).

## Erfolgskriterien

- [ ] `∫ρ·w ≈ N_e` (Elektronennummer, xcDFT `electron_number.f90`) ≤1e-3 auf He/H2O.
- [ ] Gitterkonvergenz: E_LDA stabil bei SG-1→SG-2→SG-3 (Δ <1e-4 Eh).
- [ ] Lebedev-Symmetrie (Gewichte/Integrale Sphärenfunktionen) ≤1e-10.

## Quellen & Lektüre

- xcDFT-Grid (Hauptlektüre, direkt portieren): `dft_grid.f:3` (`EulMac` —
  Punkte `R·(i/(N-i+1))²`, Gewichte `2R³(N+1)i⁵/(N-i+1)⁷`), `:35` (`Lebdev` —
  Laikov-Generator mit DATA-Blöcken für N=6…590; auf 50/110/194/302 kürzen);
  Block-Kommentar `:1-33` (Formeln + Zitate Zh. Vychisl. Mat. Mat. Fiz. 1975/76,
  Sibirsk. Mat. Zh. 1977 — in Doku übernehmen).
- `quadrature_grid.f90:33-46` (radial×angular), `read_grid.f90:21-43`
  (SG-0=23×170, SG-1=50×194, SG-2=75×302, SG-3=99×590 — Tabelle portieren).
- AO auf Grid (portieren): `AO_values_grid.f90:1` (Shell→`generate_shell.f90:1`
  kartesische Potenzen; `:75` Primitive `d·NormCoeff·exp(-αr²)`, `:94`
  Polynom `xA^ax·yA^ay·zA^az`, kartesische Gradienten `:78-90` — auch Basis für
  WP8); `generate_shell.f90:1`, `NormCoeff.f90:27-29`.
- Dichte/∇ρ (portieren): `density.f90:26-32` (`ρ=ΣP_μν AO_μ AO_ν`),
  `gradient_density.f90:1` (∇ρ, in xcDFT **nie aufgerufen** — hier aktivieren),
  `electron_number.f90:18` (`∫ρw` Sanity-Check, als Testkriterium nutzen).
- Becke-Partition (**neu**, nicht in xcDFT): Becke, JCP 88, 2547 (1988);
  Atomradien (CSD) in xcDFT `elements.f90:84-89` als Startpunkt für die
  Skalierungsfunktion f(r) = (3/2)r − (1/2)r³.
- Becke-Fix-Begründung (Doku): xcDFT-Grid ist einzentrish → nur Atome; curcuma
  braucht Mehrzentren für Moleküle.

## Validierungsergebnisse (nach Ausführung eintragen)

```
∫ρw vs N_e (He/H2O):     [ ]
Gitterkonvergenz ΔE:     [ ]
Lebedev-Sphärenharmonischen: [ ]
```