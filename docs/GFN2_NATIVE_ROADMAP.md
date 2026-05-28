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
| **Repulsion** | `xtb_native.cpp:263` (`m_E_repulsion`)<br>`calcRepulsionEnergy()` | `external/tblite/src/tblite/repulsion/effective.f90` | Implementiert, GFN2-effective | **ungemessen** |
| **Coulomb ES2 (γ)** | `xtb_native.cpp:208` (`m_E_coulomb_shell`)<br>`energyCoulombShell()` | `external/tblite/src/tblite/coulomb/shell.f90` | Implementiert, shell-resolved | **ungemessen** |
| **Multipol-Elektrostatik** | `xtb_native.cpp:210` (`m_E_multipole`)<br>`energyMultipole()`, `xtb_multipole.cpp` | `external/tblite/src/tblite/coulomb/multipole.f90` | Implementiert, AP6b Multipol-Integrale ≤ 1e-15 vs Dump | bei *fixierter Dichte* exakt (siehe `diff_multipole_potential.py`) |
| **Third-order Coulomb** | `xtb_native.cpp:209` (`m_E_third_order`)<br>`energyThirdOrder()` | `external/tblite/src/tblite/coulomb/thirdorder.f90` | Implementiert, geht in `pot%vsh` | **ungemessen** |
| **D4-Dispersion** | `xtb_native.cpp:265` (`m_E_dispersion`)<br>`calcDispersionEnergy()`, dispersion module | dftd4-src + tblite/disp/d4.f90 | Energie + analytischer Gradient + SCF-Potential + ∂E/∂q | sub-µEh H₂O/NH₃; ~5% Rest auf C-Pfad (CH₄/triose) — siehe [`GFN2_D4_STATUS.md`](GFN2_D4_STATUS.md) |
| **Halogen-Bond** | `xtb_native.cpp:264` (`m_E_halogen_bond`)<br>`calcHalogenBondEnergy()` | tblite halogen container | Implementiert; nur für Halogen-haltige Systeme aktiv | **ungemessen** (für H/C/N/O = 0) |
| **EEQ-Ladungen** (D4) | `dispersion/d4_charge_model.cpp` | dftd4-src/charge.f90 | Single-shot + analytisches ∂q/∂x | validiert via `test_d4_dedq` |
| **CN (D4-Pfad)** | `dispersion/d4_ncoord.cpp` | dftd4-src/ncoord.f90 | Exact cpp-d4 ncoord_d4 Port (EN-weighted) | bit-identisch zur Formel ✓ |
| **CN (GFN2-H0-Pfad)** | `gfn2.h:142` exp-logistisch | tblite/ncoord/gfn.f90 | Implementiert | **ungemessen** |
| **SCF-Iteration / DIIS** | `xtb_native.cpp:158-200` | tblite/scf/iterator.f90 | Damped warmup + DIIS ab Iter 5 | divergiert auf complex (231 Atome) — siehe [`GFN2_SCF_STATUS.md`](GFN2_SCF_STATUS.md) |

## Aktueller Ist-Stand (Total-Energy vs. tblite)

| Molekül | ΔE | Quelle (vermutet) |
|---------|-----|-------------------|
| H₂      | -3.3e-9   | exakt |
| H₂O     | -0.67 µEh | sub-µEh (Multipol/D4 sauber) |
| NH₃     | -1.42 µEh | sub-µEh |
| CH₄     | -1.79e-5  | C-Pfad D4-Rest |
| triose  | -3.2 mEh  | 66 × C-Pfad-Rest, **ggf. + andere Komponente** |
| complex | divergiert | SCF charge-sloshing, NICHT D4 |

**Wichtig**: Triose's 3.2 mEh sind möglicherweise *nicht alleine D4*. Eine zweite
Komponente (Repulsion, Coulomb-γ, third-order, …) könnte mit einem ähnlichen
~5%-Effekt mitlaufen, aber durch die Total-Energy-Aggregation maskiert sein. Erst
Komponenten-für-Komponenten zeigt, ob D4 wirklich der dominante Rest ist.

## Test-Infrastruktur

Bereits vorhanden (`ctest -L sqm_*` und Module wie `xtb_*`):
- `sqm_reference_dump_*`: Generiert tblite-Dumps in `release_tblite/dumps/*_gfn2.json` (USE_TBLITE)
- `sqm_overlap_*`: `test_xtb_overlap` ≡ overlap S vs Dump
- `sqm_h0_*`: `test_xtb_h0` ≡ H0 vs Dump
- `sqm_coulomb_*`: `test_xtb_coulomb` ≡ Coulomb γ vs Dump (offen: was wird genau verglichen?)
- `sqm_scf_*`: `test_xtb_scf_snapshot` ≡ SCF-Konvergenz von tblite-Dichte
- `d4_diag_*`: `diag_curcuma_d4_potential` ≡ D4-SCF-Potential vs `vat_tblite - gf2.vat_extra`

**Zu ergänzen** (siehe [`GFN2_COMPONENT_AUDIT_PLAN.md`](GFN2_COMPONENT_AUDIT_PLAN.md)):
- `diag_curcuma_repulsion` — pairwise E_rep + per-atom Beitrag
- `diag_curcuma_coulomb_es2` — E_ES2 + per-shell q·γ·q'
- `diag_curcuma_third_order` — E_3 + per-shell q²·Γ
- `diag_curcuma_multipole_energy` — E_mp Aufschlüsselung (SD/DD/SQ)
- `diag_curcuma_d4_c6` — per-Paar C₆(iat,jat) vs tblite (offene D4-Tiefe)

Erfordert tblite-Patch für eine generische `tblite_get_energy_components(...)` C-API,
die jede Container-Energie einzeln zurückgibt.

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
