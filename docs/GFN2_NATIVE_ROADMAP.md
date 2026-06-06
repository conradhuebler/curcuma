# GFN2-xTB native — Alignment-Roadmap

**Ziel:** Die native GFN2-Implementierung (`-method gfn2`, kanonischer Backend in
`release/`) bit-identisch zur Fortran-Referenz (tblite, in
`external/tblite` / `release_tblite/_deps/dftd4-src`). Sub-µEh Total-Energy
und sub-1e-5 Eh Orbital-Energie für alle Test-Moleküle vs. tblite.

**Strategie:** Per-Komponenten-Vergleich bei *fixierter Dichte*. Erst wenn jede
Komponente einzeln zu tblite passt, kann der SCF-Fixpunkt match'en und das
Gesamtenergie-Residuum verschwindet. Reine Total-Energy-Diffs (was wir bisher
hatten) verstecken kompensierende Fehler — die brauchen wir, um vom 3 mEh
Triose-Residuum runter auf µEh zu kommen.

## Energie-Komponenten

| Komponente | Curcuma-Datei | Tblite-Referenz | Curcuma-Status | Audit-Status |
|------------|---------------|------------------|----------------|--------------|
| **H0 / Tr(P·H0)** | `xtb_native.cpp:213` (`m_E_electronic`)<br>`xtb_h0.cpp::getHamiltonianH0` | `external/tblite/src/tblite/xtb/h0.f90` | Implementiert, EHT-Stil; `test_xtb_h0` validiert vs Dump | ✓ existiert |
| **Overlap (S)** | `xtb_native.cpp::buildBasis` / STO-CGTO | `external/tblite/src/tblite/integral/overlap.f90` | Implementiert; `test_xtb_overlap` validiert vs Dump | ✓ existiert |
| **Repulsion** | `xtb_native.cpp:263` (`m_E_repulsion`)<br>`calcRepulsionEnergy()` | `external/tblite/src/tblite/repulsion/effective.f90` | Implementiert, GFN2-effective | bit-identisch (≤ 1e-17, alle Test-Moleküle) ✓ |
| **Coulomb ES2 (γ)** | `xtb_native.cpp:208` (`m_E_coulomb_shell`)<br>`energyCoulombShell()` | `external/tblite/src/tblite/coulomb/shell.f90` | Implementiert, shell-resolved | im Electronic-Lump enthalten — Lump bit-identisch H₂/H₂O/NH₃ (siehe `ecomp_*`) |
| **Multipol-Elektrostatik** | `xtb_native.cpp:210` (`m_E_multipole`)<br>`energyMultipole()`, `xtb_multipole.cpp` | `external/tblite/src/tblite/coulomb/multipole.f90` | Implementiert, AP6b Multipol-Integrale ≤ 1e-15 vs Dump | bei *fixierter Dichte* exakt (siehe `diff_multipole_potential.py`) |
| **Third-order Coulomb** | `xtb_native.cpp:209` (`m_E_third_order`)<br>`energyThirdOrder()` | `external/tblite/src/tblite/coulomb/thirdorder.f90` | Implementiert, geht in `pot%vsh` | im Electronic-Lump enthalten — Lump bit-identisch H₂/H₂O/NH₃ |
| **D4-Dispersion** | `xtb_native.cpp:265` (`m_E_dispersion`)<br>`calcDispersionEnergy()`, dispersion module + `D4Evaluator::computeATM` | dftd4-src + tblite/disp/d4.f90 | Energie + analytischer Gradient + SCF-Potential + ∂E/∂q + **ATM-3-Body** | **GEFIXT (2026-05-29), sub-nEh vs tblite (von 3.2 mEh triose):** (1) `refcn`→`refcovcn` in der CN-Gauß-Gewichtung (2-Body C6 war 14% daneben auf C-C); (2) ATM-3-Body-Term implementiert (`computeATM`, Port von dftd4 `get_atm_dispersion`, q=0-C6, kein q-response; Gradient FD-validiert Maschinengenauigkeit). triose disp+elec 5.9e-9. Siehe [`GFN2_D4_STATUS.md`](GFN2_D4_STATUS.md) |
| **Halogen-Bond** | `xtb_native.cpp:264` (`m_E_halogen_bond`)<br>`calcHalogenBondEnergy()` | tblite halogen container | Implementiert; nur für Halogen-haltige Systeme aktiv | **ungemessen** (für H/C/N/O = 0) |
| **EEQ-Ladungen** (D4) | `dispersion/d4_charge_model.cpp` | dftd4-src/charge.f90 | Single-shot + analytisches ∂q/∂x | validiert via `test_d4_dedq` |
| **CN (D4-Pfad)** | `dispersion/d4_ncoord.cpp` | dftd4-src/ncoord.f90 | Exact cpp-d4 ncoord_d4 Port (EN-weighted) | bit-identisch zur Formel ✓ |
| **CN (GFN2-H0-Pfad)** | `gfn2.h:142` exp-logistisch | tblite/ncoord/gfn.f90 | Implementiert | **ungemessen** |
| **SCF-Iteration / DIIS** | `xtb_native.cpp:158-200` | tblite/scf/iterator.f90 | Damped warmup + DIIS ab Iter 5 | divergiert auf complex (231 Atome) — siehe [`GFN2_SCF_STATUS.md`](GFN2_SCF_STATUS.md) |

## Aktueller Ist-Stand (Total-Energy vs. tblite)

Werte sind das `disp+elec`-Residuum bei *fixierter tblite-Dichte*
(`diag_curcuma_energy_components`) — die apples-to-apples-Komponentenzahl.

| Molekül | vor refcovcn | + refcovcn (2-Body) | + ATM (3-Body) | Quelle |
|---------|--------------|---------------------|----------------|--------|
| H₂      | -3.3e-9   | -3.3e-9   | -3.3e-9   | exakt |
| H₂O     | -0.67 µEh | 3.0e-10   | **2.7e-10** | im Wesentlichen exakt |
| NH₃     | -1.42 µEh | 1.0e-9    | **1.4e-10** | im Wesentlichen exakt |
| CH₄     | -1.79e-5  | 5.9e-9    | **2.4e-11** | C-Pfad-Rest GEFIXT |
| triose  | -3.2 mEh  | -0.78 mEh | **5.9e-9**  | refcovcn + ATM → sub-nEh |
| complex | divergiert | divergiert | divergiert | SCF charge-sloshing, NICHT D4 |

**Update 2026-05-29**: GFN2-D4 jetzt sub-nEh vs tblite (von 3.2 mEh triose).
Zwei Bugs gefixt: (1) CN-Gauß-Gewichtung benutzte `refcn` statt `refcovcn`
(2-Body C6 14% daneben auf C-C); (2) ATM-3-Body-Term fehlte komplett
(`D4Evaluator::computeATM`, Port von dftd4 `get_atm_dispersion`, Energie +
Gradient FD-validiert auf Maschinengenauigkeit). Repulsion, Referenz-C6,
gewichteter 2-Body-C6, ATM — alle bit-identisch / sub-nEh. GFN-FF unberührt.

(Historisch) **Wichtig**: Triose's 3.2 mEh sind möglicherweise *nicht alleine D4*. Eine zweite
Komponente (Repulsion, Coulomb-γ, third-order, …) könnte mit einem ähnlichen
~5%-Effekt mitlaufen, aber durch die Total-Energy-Aggregation maskiert sein. Erst
Komponenten-für-Komponenten zeigt, ob D4 wirklich der dominante Rest ist.

**Phase A+B+C beantwortet das (2026-05-28):** Bei *injizierter tblite-Dichte*
reproduziert die kombinierte (disp+elec)-Diff den dokumentierten Total-Energy-
Rest aus der obigen Tabelle exakt — H₂ 1.5e-10, H₂O 0.67 µEh, NH₃ 1.42 µEh,
CH₄ 1.79e-5, triose 3.21 mEh. Repulsion stimmt zu ≤ 1e-17 (bit-identisch).
Damit ist gezeigt: **der triose-3.2-mEh-Rest steckt in der D4+electronic-
Maschinerie, nicht in Repulsion**. Weiterer Schritt: Phase A2 ("Fine") oder
direkter `d4_model%c6`-Diff (siehe `GFN2_D4_STATUS.md`), um zwischen
SCC-D4-baseline und ES2/3rd/multipole-Drift zu trennen.

## Test-Infrastruktur

Bereits vorhanden (`ctest -L sqm_*` und Module wie `xtb_*`):
- `sqm_reference_dump_*`: Generiert tblite-Dumps in `release_tblite/dumps/*_gfn2.json` (USE_TBLITE)
- `sqm_overlap_*`: `test_xtb_overlap` ≡ overlap S vs Dump
- `sqm_h0_*`: `test_xtb_h0` ≡ H0 vs Dump
- `sqm_coulomb_*`: `test_xtb_coulomb` ≡ Coulomb γ vs Dump (offen: was wird genau verglichen?)
- `sqm_scf_*`: `test_xtb_scf_snapshot` ≡ SCF-Konvergenz von tblite-Dichte
- `d4_diag_*`: `diag_curcuma_d4_potential` ≡ D4-SCF-Potential vs `vat_tblite - gf2.vat_extra`

**Zu ergänzen** (siehe [`GFN2_COMPONENT_AUDIT_PLAN.md`](GFN2_COMPONENT_AUDIT_PLAN.md)):
- ✅ `diag_curcuma_energy_components` (2026-05-28, Phase A+B+C) — Coarse-Variante: tblite
  liefert per-container `halogen / repulsion / dispersion / interactions / electronic-lump`
  via fünf neue `tblite_get_result_energy_component_*` C-APIs. Curcuma injiziert die
  tblite-Dichte (`XTB::evaluateComponentsAtFixedDensity`) und vergleicht jede Container-
  Energie. Repulsion bit-identisch; combined |disp+elec| reproduziert den Total-Diff
  exakt. ctest-label `gfn2_align`, 5/5 grün.
- offen (Fine, optional): tblite-internen GFN2-Coulomb-Container in ES2 / 3rd / mp
  aufschlüsseln — nur nötig falls die Phase-C-Sub-pieces-Auswertung nicht ausreicht,
  einen Sub-Term zu lokalisieren
- offen: `diag_curcuma_d4_c6` — per-Paar C₆(iat,jat) vs tblite (offene D4-Tiefe)

## Priorisierte Arbeit

**P1 — Triose vom 3 mEh runter (höchste Priorität):**
1. Tblite-Patch für per-Komponenten-Energie-Dump (siehe Component-Audit-Plan).
2. Component-Diff curcuma↔tblite auf CH₄ und triose. Zeigt, ob D4 wirklich der
   dominante Rest ist oder ob auch Repulsion/Coulomb-γ mitlaufen.
3. Wenn D4 dominant: tblite-Patch für `d4_model%c6`-Export, dann per-(iref,jref)
   Diff gegen `m_c6_flat_cache` — siehe [`GFN2_D4_STATUS.md`](GFN2_D4_STATUS.md).
4. Falls eine andere Komponente betroffen: einzelnen Bug fixen.

**P2 — Complex zur Konvergenz bringen** (orthogonal zur Genauigkeit):
- SCF-Stabilität, DIIS-Tuning. Siehe [`GFN2_SCF_STATUS.md`](GFN2_SCF_STATUS.md).

**P3 — Restliche ungemessene Komponenten** (Repulsion, Coulomb-γ, third-order):
- Nach P1 kann jeder einzeln auditiert werden. Aktuell vermutlich OK, aber
  unverifiziert.

## Verwandte Dokumente

- [`GFN2_D4_STATUS.md`](GFN2_D4_STATUS.md) — D4-Tiefe, Audit-Trail
- [`GFN2_COMPONENT_AUDIT_PLAN.md`](GFN2_COMPONENT_AUDIT_PLAN.md) — konkrete
  Patches und Diff-Strategie für jede Komponente
- [`GFN2_SCF_STATUS.md`](GFN2_SCF_STATUS.md) — Konvergenz-Probleme (complex)
- [`AP6b_multipole_F_discrepancy_WP.md`](AP6b_multipole_F_discrepancy_WP.md)
- [`PHASE3B4_MULLIKEN_RESPONSE_WP.md`](PHASE3B4_MULLIKEN_RESPONSE_WP.md) — Gradient
- [`PHASE3B5_MULTIPOLE_RESPONSE_WP.md`](PHASE3B5_MULTIPOLE_RESPONSE_WP.md) — Gradient

## Wo der Code lebt

- Native GFN2 driver: `src/core/energy_calculators/qm_methods/xtb_native.cpp`
- Sub-Module: `xtb_h0.cpp`, `xtb_multipole.cpp`, `xtb_gradient.cpp`,
  `xtb_response.cpp`, `diis_accelerator.h`
- Parameter: `parameters/gfn2_params.hpp`
- Dispersion: `src/core/energy_calculators/dispersion/`
- Test-Helpers: `test_cases/sqm_reference/`
- Tblite: `external/tblite/` (Fetched, kann gepatcht werden)
- Dftd4-src: `release_tblite/_deps/dftd4-src/` (über tblite FetchContent)
