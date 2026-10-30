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
| WP4 | DFT-Grid (Treutler-Ahlrichs/Euler-Maclaurin + Lebedev + Becke/Stratmann, `qm_grid.*`) | ⚙️ machine-tested (Okt 2026, `ctest -L qm_grid` 19/19: Punkte/Gewichte vs PySCF ≤8e-13, ∫ρ vs PySCF ≤2e-13, N_el-Tabelle H2O..Benzol); nur H-Ne, geschlossene Schale, E_xc-Konvergenz offen (WP5) | `WP4_dft_grid_becke.md` |
| WP5 | LDA (Slater-Dirac + VWN5/VWN-RPA, `qm_xc.*`, `-method lda`) | ⚙️ machine-tested (Okt 2026, `ctest -L qm_lda` 20/20: Kernels vs libxc ≤6e-14 rel., E vs PySCF gleiches Gitter ≤5e-12 Eh, 12 Fälle inkl. RI; vs ORCA 6.1 LSD ≤5.4e-8 Eh, Benzol 4e-7, Gitter nicht identisch); Gradient siehe WP8, geschlossene Schale, H-Ne | `WP5_lda_vwn5.md` |
| WP6 | PBE-GGA (PBE-X + PBE-C mit PW92, `qm_xc.*`, `-method pbe`) | ⚙️ machine-tested (Okt 2026, `ctest -L qm_pbe` 21/21: Kerne vs libxc ≤5.1e-15 skaliert, V_xc vs FD ≤7.8e-9, E vs PySCF gleiches Gitter ≤3.9e-12 Eh, 11 Fälle inkl. RI; vs ORCA 6.1 mit libxc-μ ≤7e-8 Eh, Benzol 1e-6 (ORCA-Winkelgitter), ORCA-Standard-μ 0.21952 gibt 4e-6..5.5e-5); Gradient siehe WP8, geschlossene Schale, H-Ne | `WP6_pbe_gga.md` |
| WP7 | B3LYP-Hybrid (0.08 Slater + 0.72 B88 + 0.20 HF + 0.19 VWN + 0.81 LYP, `qm_xc.*`, `-method b3lyp`, VWN per `-qm.vwn_variant`, Default vwn5 = ORCA B3LYP) | ⚙️ machine-tested (Okt 2026, `ctest -L qm_b3lyp` 29/29: Kerne vs libxc ≤3.2e-14 skaliert, E vs PySCF gleiches Gitter ≤3.0e-11 Eh, 19 Fälle beide VWN-Varianten inkl. RI-J/RI-JK; vs ORCA 6.1 B3LYP und B3LYP/G ≤5.3e-8 Eh, Benzol 6.2e-7 (ORCA-Winkelgitter; curcuma 250x590 3.6e-9)); Gradient siehe WP8, geschlossene Schale, H-Ne | `WP7_b3lyp_hybrid.md` |
| WP8 | Analytischer Gradient (HF; KS: a_x-gewichtetes K, XC-Basisteil mit AO-Zweitableitungen, Becke-Gitterantwort `-qm.grid_response`, Default true) | ⚙️ machine-tested: `hf`/`hf-3c` (Sep 2026, vs PySCF ≤2e-10 Eh/Bohr); `lda`/`pbe`/`b3lyp` (Okt 2026, `ctest -L qm_dft_grad` 55/55: vs PySCF gleiches Gitter mit/ohne Gitterantwort ≤3.6e-11 Eh/Bohr, 27 Fälle inkl. RI-J/RI-JK, beide VWN; vs FD der eigenen Energie ≤3.4e-10; vs ORCA 6.1 fein 6e-8..1.3e-6 (Gitter); `-opt` erreicht das PySCF-Minimum); geschlossene Schale, H-Ne; GPU siehe unten | `WP8_analytischer_gradient.md` |
| WP9 | Integration, CLI, Parameter, Doku, CTest | ⚠️ offen | `WP9_integration_doku_ctest.md` |
| WP10 | RI-J / RI-JK fuer HF (`-qm.ri j\|jk`, def2/J, def2/JK, CPU + CUDA) inkl. analytischem Gradienten (`-opt`) | ⚙️ machine-tested (Okt 2026, Energie vs PySCF ≤1e-10 Eh, vs ORCA 6.1 ≤2e-9 Eh, 5 Molekuele; Gradient vs PySCF ≤6.1e-11 Eh/Bohr, vs ORCA ≤8.8e-9 Eh/Bohr); RI-C offen | `WP10_korrelation_ri_mpn_ccsd.md` |
| WP11 | kanonisches RMP2/SCS/SOS (Schritt 1+3, CPU + CUDA, `-method mp2`) und sein analytischer Gradient (Z-Vektor, frozen core; CPU, 2e-Teil auch CUDA) | ⚙️ machine-tested (Okt 2026, Energie vs PySCF ≤1.8e-9 Eh, vs ORCA ≤2.0e-9 Eh, 16 Faelle; Gradient vs PySCF-FD ≤1.4e-9 Eh/Bohr, 16 Faelle, vs ORCA ≤7.7e-8, `-opt` erreicht das PySCF-Minimum); RI-MP2 eigene Zeile | `WP10_korrelation_ri_mpn_ccsd.md` |
| WP11/2 | RI-MP2/SCS/SOS (`-method ri-mp2`, Hilfsbasis def2-universal-JKFIT/-JFIT statt def2-SVP/C: Methodenvariante; CPU + CUDA) | ⚙️ machine-tested (Okt 2026, vs PySCF mit denselben Hilfsbasisdateien ≤1.7e-9 Eh, vs ORCA 6.1 AuxC def2/JK ≤1.8e-9 Eh, 16 Faelle + JFIT/RI-JK-SCF); def2-SVP/C und Gradient offen | `WP10_korrelation_ri_mpn_ccsd.md` |
| WP11 (Meta-GGA) | TPSS und SCAN (`-method tpss\|scan`, τ-abhängige Kerne mit Dual-Zahl-Ableitungen in `qm_xc_kernels.h`, V_xc mit ½ v_τ ∇φ·∇φ) | Stufe 1 (Energien, CPU) ⚙️ machine-tested (30. Okt 2026, `ctest -L qm_mgga` 50/50: Kerne vs libxc ≤7.7e-11 skaliert auf 1760 Punkten, E vs PySCF gleiches Gitter ≤1.4e-11 Eh, 28 Fälle inkl. RI; vs ORCA 6.1 LibXC(TPSS) fein ≤8.5e-8 Eh, Benzol 2.4e-6; SCAN stark gitterempfindlich: sg1-sg3 bis 3e-3 Eh, 1e-6 Eh erst bei ~400x590); Stufe 2 (Gradient, CPU) ⚙️ machine-tested (30. Okt 2026, `ctest -L qm_mgga_grad`: vs PySCF gleiches Gitter ≤9.7e-12 Eh/Bohr mit/ohne Gitterantwort, FD ≤2e-10 mitwanderndes Gitter, `-opt` TPSS); Stufe 3 (GPU) offen, geschlossene Schale, H-Ne | `WP11_meta_gga.md` |
| WP12-15 | RCCSD, lokales PNO-MP2 (WP15), Bewertung Schwarz-Ansaetze (WP12 gestrichen) | ⚠️ nur Plan | `WP10_korrelation_ri_mpn_ccsd.md` |

**KS-DFT auf der GPU** (30. Okt 2026, `-gpu cuda`, `-qm.gpu_xc`): V_xc/E_xc, die KS-SCF-Schleife (geräteresident) und der KS-Gradient (XC-Teil mit Gitterantwort, 2e-Teil mit kx) laufen auf dem Gerät; Formelkerne einmal in `qm_xc_kernels.h` für CPU und GPU. ⚙️ machine-tested nur auf einer RTX 5090 (`ctest -L qm_gpu_dft` 34/34: V_xc/E_xc vs CPU gleiche Dichte ≤4.4e-15 relativ, SCF-Energie ≤1e-10 Eh, Gradient gleiche Dichte ≤4.5e-12 Eh/Bohr, `-opt` erreicht das CPU-Minimum). Messungen, nicht Getestetes: [QM_GPU_ROADMAP.md Abschnitt 4](../QM_GPU_ROADMAP.md#4-kohn-sham-dft-on-the-gpu-oct-2026).

**AO-Screening pro Gitterblock und Gewichtsableitung ohne Atomgrenze** (30. Okt 2026, `-qm.xc_ao_threshold`, Default 1e-12, 0 = alter Pfad bitgleich): Punkte jedes Atoms in kompakte Blöcke sortiert, pro Block nur die Schalen mit Reichweite (`AOEvaluator::shellRadii`, eine Formel für CPU und GPU); Energie/Gradient gegen Schwelle 0 ≤6.4e-12 Eh / ≤1.4e-12 Eh/Bohr (H2O bis C40H82), V_xc pro Aufbau C40H82 CPU 23.2 → 0.92 s, GPU 3.0 → 0.49 s, nahezu linear ab C20; Becke-Gewichtsableitung auf dem Gerät ohne 700-Atom-Grenze (gleiche Bits). 🤖 AI-generated, ⚙️ machine-tested; Tabellen: [QM_GPU_ROADMAP.md Abschnitt 4.1](../QM_GPU_ROADMAP.md#41-ao-screening-per-block-and-weight-derivative-without-atom-limit-oct-30-2026).

Vor jedem WP: `make -j4` in `release/` grün ohne Warning. Nach jedem WP: `ctest`.
Volle Plandatei: `/home/conrad/.claude/plans/wir-wollen-dft-in-clever-gizmo.md`.