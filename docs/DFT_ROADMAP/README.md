# Native KS-DFT in Curcuma — Implementierungs-Roadmap

Didaktische, referenzvalidierte Kohn-Sham-DFT für curcuma, portiert/erweitert
aus **xcDFT (TCCM winter school 2019: DFT)** unter
`/home/conrad/src/claude_curcuma/xcDFT`. Status-Labels nach CLAUDE.md
(⚠️ AI-generated / ⚙️ machine-tested; kein selbst vergebenes ✅ TESTED/APPROVED).

## Method-Namen

**Keine Sammel-`-method dft`.** Jedes Funktional ist ein eigener Method-Name
(wie `gfn1`/`gfn2`/`pm3` in curcuma): `-method hf`, `-method lda`, `-method pbe`,
`-method b3lyp`. Der Functional-Wert wird im `MethodFactory` auf
`QMMethod(QMFunctional::..., config)` abgebildet (bis Sep 2026 `DFTMethod`/`DFTFunctional`).
Gemeinsame Parameter (Basis, Grid, SCF-Optionen) liegen im Parameter-Modul `qm`
(`-qm.basis`, `-qm.grid`, `-qm.scf_*`; das alte `-dft.*` wird weiter gelesen);
`functional` ist **kein** Parameter. HF-3c ist seit Sep 2026 nativ (`-method hf-3c`),
siehe [NATIVE_QM_IMPLEMENTATION.md](../NATIVE_QM_IMPLEMENTATION.md).

## Quellenangabe

- **Ursprung:** xcDFT / TCCM winter school 2019: DFT (Fortran, RKS, LDA+HF,
  Euler-Maclaurin×Lebedev-Grid, DIIS, Gauß-Basis). Direkt portierte Module
  tragen `// Ported from xcDFT (TCCM winter school 2019: DFT), src/<file>.f90`.
- **curcuma-native Erweiterungen** (nicht in xcDFT): Becke-Mehrzentren-Grid
  (xcDFT-Grid ist einzentrish → nur Atome), GGA/Hybrid (xcDFT-Stubs), ERI-Engine
  (xcDFT liest vorberechnete Integrale), analytischer Gradient, VWN5-Korrektur
  der fehlerhaften xcDFT-LDA-Korrelationsformel.
- **Zweite Referenz:** ORCA 6.1 (`/opt/orca_6_1`, via `ORCA_PATH=/opt/orca_6_1/orca`,
  bestehendes curcuma-ORCA-Interface) mit def2-SVP.

## WP-Reihenfolge & Abhängigkeitsgraph

```
WP0 (Gerüst/Setup)
 └─ WP1 (GTO-1e-Integrale)
     └─ WP2 (4-Zentren-ERI, McMurchie-Davidson)
         ├─ WP3 (HF-SCF, hartes ERI-Gate)
         │   └─ WP4 (DFT-Grid: Euler-Maclaurin + Lebedev + Becke)
         │       └─ WP5 (LDA: Slater-Dirac + VWN5)
         │           └─ WP6 (PBE-GGA, +∇ρ)
         │               └─ WP7 (B3LYP-Hybrid)
         │                   └─ WP8 (Analytischer Gradient)
         └────────────────────┴─ WP9 (Integration/CLI/Doku/CTest)  [benötigt WP1–WP8]
```

## Status-Tabelle

| WP | Titel | Status | Datei |
|----|-------|--------|-------|
| WP0 | Gerüst, Quellenangabe, Setup | ⚙️ machine-tested | `WP0_geruest_setup.md` |
| WP1 | GTO-1e-Integrale (S/T/V) | ⚙️ machine-tested (`ctest -L qm_1e` 10/10) | `WP1_gto_1e_integrale.md` |
| WP2 | 4-Zentren-ERI (McMurchie-Davidson) | ⚙️ machine-tested (`ctest -L qm_2e` 10/10) | `WP2_eri_mcmurchie_davidson.md` |
| WP3 | HF-SCF (rung 666) — hartes Gate | ⚙️ machine-tested (10/10 vs ORCA ≤4e-9) | `WP3_hf_scf_gate.md` |
| WP4 | DFT-Grid (Euler-Maclaurin + Lebedev + Becke) | ⚠️ offen | `WP4_dft_grid_becke.md` |
| WP5 | LDA (Slater-Dirac + VWN5) | ⚠️ offen | `WP5_lda_vwn5.md` |
| WP6 | PBE-GGA (+ ∇ρ) | ⚠️ offen | `WP6_pbe_gga.md` |
| WP7 | B3LYP-Hybrid (+ exakter Austausch) | ⚠️ offen | `WP7_b3lyp_hybrid.md` |
| WP8 | Analytischer Gradient | ⚙️ machine-tested für `hf`/`hf-3c` (Sep 2026, vs PySCF ≤2e-10 Eh/Bohr); DFT-Anteil offen | `WP8_analytischer_gradient.md` |
| WP9 | Integration, CLI, Parameter, Doku, CTest | ⚠️ offen | `WP9_integration_doku_ctest.md` |

Vor jedem WP: `make -j4` in `release/` grün ohne Warning. Nach jedem WP: `ctest`.
Volle Plandatei: `/home/conrad/.claude/plans/wir-wollen-dft-in-clever-gizmo.md`.