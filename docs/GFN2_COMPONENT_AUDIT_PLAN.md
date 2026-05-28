# GFN2 Component Audit Plan

Konkretes Vorgehen, um curcumas native GFN2 Komponente-für-Komponente gegen
tblite (Fortran) zu validieren. Master-Roadmap: [`GFN2_NATIVE_ROADMAP.md`](GFN2_NATIVE_ROADMAP.md).

> **Phase A+B+C (Coarse) implementiert, 2026-05-28** — `ctest -L gfn2_align`,
> 5/5 grün. Details siehe Abschnitt "Status" am Ende.

## Grundprinzip

Bei *fixierter tblite-Dichte* die Energie jeder Komponente parallel berechnen
und auf 1e-9 vergleichen. Wenn jede einzelne Komponente bei festem P stimmt,
ist die Gesamtsumme exakt. Vorteil gegenüber Total-Energy-Diff: kompensierende
Fehler werden sichtbar.

## Stand der bestehenden tblite-Patches

Bereits exponiert via C-API (per AP6b, in `external/tblite/include/tblite/result.h`):
- `tblite_get_result_dipole_integral`, `_quadrupole_integral`
- `tblite_get_result_atomic_dipole`, `_atomic_quadrupole`
- `tblite_get_result_potential_vat`, `_vdp`, `_vqp`
- `tblite_get_result_energy`, `_energies` (Total + per-atom — keine per-container Aufschlüsselung)

**Fehlt:** per-container Energien (E_repulsion, E_coulomb_shell, E_third_order,
E_multipole, E_dispersion, E_band/Tr(P·H0)).

## Phase A — Tblite-Patch für Energie-Container

Tblite ruft im Singlepoint-Pfad jeden Container's `get_energy` einzeln auf und
akkumuliert in `wfn%energy`. Die per-Container-Werte werden aktuell nicht
gespeichert.

**Patch-Skizze** (in `external/tblite/src/`, ähnlich zu den bestehenden
multipole-Patches):

1. **`xtb/singlepoint.f90`**: nach jedem `container%get_energy(...)`-Call den
   skalaren Beitrag (E_total minus E_total_vor_call, oder direkt aus dem Container)
   in einen neuen `results%e_components(:)` Array kopieren. Index-Mapping z.B.:
   ```
   e_components(1) = E_band  ! Tr(P·H0)
   e_components(2) = E_repulsion
   e_components(3) = E_coulomb_es2
   e_components(4) = E_third_order
   e_components(5) = E_multipole
   e_components(6) = E_dispersion
   ```
2. **`results.f90`**: `real(wp), allocatable :: e_components(:)` Feld hinzufügen.
3. **`api/result.f90`**: neuer Getter `tblite_get_result_energy_components`.
4. **`include/tblite/result.h`**: neue C-Signatur.

Damit hat die JSON-Dump-Pipeline die per-Container-Werte zum Vergleich.

**Aufwand:** ~30-60 Min (modellbar auf den existierenden multipole-Patch).

## Phase B — Dump-Tool erweitern

`test_cases/sqm_reference/dump_tblite_multipole.cpp` schreibt aktuell `gf2.*`
und `*_tblite`-Felder. Hinzufügen:

```json
"e_components_tblite": {
  "band":         <double>,
  "repulsion":    <double>,
  "coulomb_es2":  <double>,
  "third_order":  <double>,
  "multipole":    <double>,
  "dispersion":   <double>
}
```

## Phase C — Curcuma-Side Diff-Tool

Neues Tool `test_cases/sqm_reference/diag_curcuma_energy_components.cpp`,
modelliert auf `diag_curcuma_d4_potential`:

1. Lese Dump, hole tblite's `e_components_tblite` + `density` (P) + `atomic_charges` (q).
2. Konstruiere curcumas `XTB` Instanz, **injiziere tblite's konvergierte
   Dichte** (so dass curcuma das gleiche P benutzt, nicht selbst SCF rechnet),
   und ruf einmal jede `energy*()` Helper-Funktion auf:
   - `Tr(P · H0)` → band
   - `calcRepulsionEnergy()` → repulsion
   - `energyCoulombShell()` → coulomb_es2
   - `energyThirdOrder()` → third_order
   - `energyMultipole()` → multipole
   - `calcDispersionEnergy()` → dispersion
3. Per-Komponente max|diff| ausgeben + optionaler `--tol`/`--quiet` Modus
   (gleicher Stil wie diag_curcuma_d4_potential).

**Subtilität:** die curcuma-Funktionen lesen `m_wfn.P` und `m_wfn.q_at`/`q_sh`
intern. Wir müssen sie auf tblite's Werte setzen, bevor wir die energy*()
Helper aufrufen. Das vermeidet die Komplikation eines anderen SCF-Fixpunkts.

## Phase D — Ctest-Integration

Erweitere `test_cases/sqm_reference/CMakeLists.txt` analog zum bestehenden
`d4_diag`-Block:

```cmake
set(_E_COMP_MOLS_TOL
    "H2:1e-9"
    "H2O:1e-9"
    "NH3:1e-9"
    "CH4:1e-9"
    "triose:1e-7"  # 66 atoms × per-atom precision
)
foreach(_mt ${_E_COMP_MOLS_TOL})
    add_test(NAME ecomp_${_mol} COMMAND diag_curcuma_energy_components --quiet --tol ${_tol} ...)
    set_tests_properties(ecomp_${_mol} PROPERTIES LABELS "gfn2_align")
endforeach()
```

Run via `ctest -L gfn2_align`.

## Phase E — Lokalisieren und Fixen

Mit den Per-Komponenten-Werten:

| Diff-Pattern                                | Implikation |
|---------------------------------------------|-------------|
| Nur D4 weicht ab                            | C-Pfad-Audit fortsetzen (`d4_model%c6` Diff) |
| Repulsion weicht ab                         | Repulsion-Parameter / Damping prüfen |
| Coulomb-ES2 weicht ab                       | γ-Matrix-Konstruktion, Hubbard-U Tabelle |
| Third-order weicht ab                       | `p_hubbard_derivs` Tabelle, shell-resolved Logik |
| Multipole bei tblite-Dichte weicht ab       | (überraschend; AP6b sagte exakt) — amat_sd/dd/sq, mrad |
| Mehrere weichen ab mit kleinen Beträgen     | sammelt sich auf große Total-Diff (triose-Szenario) |

## Reihenfolge der Arbeit

1. **Phase A (tblite-Patch)** — einmaliger Aufwand.
2. **Phase B (Dump erweitern)** — neue Dumps generieren für H₂O/CH₄/NH₃/triose.
3. **Phase C (curcuma-Diff-Tool)** — Hauptarbeit, ~2-3 h.
4. **Phase D (ctest)** — schnell, sobald C steht.
5. **Phase E (Bug-Fixes)** — abhängig von Phase-C-Output.

## Open-Question-Box

- Wie injiziert man tblite's Dichte in curcuma's `XTB` Instanz, ohne den ganzen
  SCF-Pfad zu durchlaufen? Vermutlich: nach `InitialiseMolecule()` direkt
  `m_wfn.P = parseDensity(dump)`, dann `m_wfn.q_at = parseCharges(dump)` etc.,
  dann die Komponenten-Berechnungen einzeln aufrufen.
- Ist die curcuma-Datenstruktur (`Wavefunction`) kompatibel mit tblite's
  Dichte-Layout (orbital-ordering, sortierung)? Vermutlich ja, da die multipole-
  AO-Integrale schon bit-identisch übereinstimmen (per AP6b).
- Was, falls Phase-A-Patch nicht alle Container so feinkörnig auflöst? Plan B:
  Komponenten-Energien aus tblite via mehrere Singlepoint-Runs mit selektivem
  Container-Off-Schalten extrahieren.

## Status (2026-05-28)

**Phase A+B+C+D Coarse-Variante implementiert.**

- **Phase A (tblite-Patch)**: 5 neue C-API-Getter
  `tblite_get_result_energy_component_{halogen,repulsion,dispersion,interactions,electronic}`,
  jeweils nat-resolved. `results_type` um 5 `ecomp_*` Felder erweitert,
  `xtb/singlepoint.f90` kopiert die `exbond/erep/edisp/eint`-Arrays direkt
  nach `get_engrad`, `eelec` nach SCF-Loop. Patch-Stil identisch zum AP6b-Patch.
- **Phase B**: `dump_tblite_multipole.cpp` schreibt jetzt einen
  `e_components_tblite`-Block mit `{halogen, repulsion, dispersion,
  interactions, electronic}` plus `sums.{...}`. Container ohne
  Calculator-Zuweisung werden gracefully übersprungen (tblite_check_error +
  clear_error).
- **Phase C**: `XTB::evaluateComponentsAtFixedDensity(P, q_at, q_sh, dp_at, qp_at)`
  überspringt SCF, injiziert die externe Dichte/Moments, ruft jeden
  Energie-Helper einmal auf. `XTB::getReferenceShellOccupations()` exponiert
  `n0_sh` für die `q_sh = n0 - n_sh`-Rekonstruktion (kritisch: GFN2 hat
  *nicht-ganzzahlige* Referenz, z.B. N=`{1.5, 3.5}` — eine Spiegelung der
  Tabelle wäre fehleranfällig). `diag_curcuma_energy_components` reportet
  per-Container-Diff plus die **kombinierte** disp+elec-Diff (die einzige
  apples-to-apples-Zahl, weil tblite D4 SCF-coupled in `eelec` ablegt,
  während curcuma alles in `m_E_dispersion` hat).
- **Phase D**: `ctest -L gfn2_align`, 5 Moleküle, alle grün:
  | Molekül | combined \|disp+elec\| |
  |---------|----------------------|
  | H₂      | 1.50e-10 |
  | H₂O     | 6.69e-7 |
  | NH₃     | 1.42e-6 |
  | CH₄     | 1.79e-5 |
  | triose  | 3.21e-3 |

  Reproduziert die Total-Energy-Ist-Stand-Tabelle im Master exakt.
  Repulsion separat bit-identisch (≤ 1e-17).

**Folgerung**: Der gesamte triose-3.2-mEh-Residuum lebt in der D4+electronic-
Maschinerie, *nicht* in Repulsion. Das per-Container-Bild zeigt zusätzlich,
dass tblites `dispersion`-Container nur die pre-SCF-D4-Baseline (q=0)
einfängt; die SCC-gekoppelte D4-Energie wandert in `eelec`. Curcuma packt
alles in `m_E_dispersion`. Die Differenz pro Container ist also gross
(z.B. 7e-2 für triose), die Summe stimmt aber bis auf das 3.2-mEh-Residuum.

**Nächster Schritt** (P1 fortgeführt): per-Paar `d4_model%c6`-Export aus
tblite + curcuma `m_c6_flat_cache`-Diff, um die D4-Tiefe zu auditieren.
Oder Phase A2 (Fine): GFN2-Coulomb-Container intern in ES2/3rd/mp splitten,
um eine zweite Komponente auszuschließen.

## Verwandte Dokumente

- [`GFN2_NATIVE_ROADMAP.md`](GFN2_NATIVE_ROADMAP.md) — Master
- [`GFN2_D4_STATUS.md`](GFN2_D4_STATUS.md) — D4-Tiefe (eine der Komponenten)
- [`GFN2_SCF_STATUS.md`](GFN2_SCF_STATUS.md) — SCF-Stabilität
