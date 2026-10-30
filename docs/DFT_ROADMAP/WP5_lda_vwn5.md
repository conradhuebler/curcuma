# WP5 — LDA (Slater-Dirac + VWN)

- **Status:** 🤖 AI-generated, ⚙️ machine-tested (`ctest -L qm_lda` 20/20, Okt 2026); menschlicher Produktionstest ausstehend
- **Abhängigkeit:** WP4 (Gitter)
- **Code:** `src/core/energy_calculators/qm_methods/qm_xc.{h,cpp}` (Kernels, V_xc), KS-Teile in `qm_scf.cpp`
  (`composeFock`, `scfEnergy`, `buildJKForSCF`, `xcFor`, `ensureXCGrid`), Test `test_cases/qm_lda/`,
  Referenz `scripts/qm_lda_pyscf_reference.py` (PySCF 2.14, libxc 7.0.0)

## Was gebaut ist

- `-method lda` = geschlossenschaliges RKS: F = h + J[P] + V_xc, E = Tr(Ph) + 1/2 Tr(PJ) + E_xc + E_nn; kein K-Build.
- Austausch Dirac/Slater, Korrelation VWN (nur paramagnetischer Fit, zeta = 0), Potential v = d(rho eps)/d rho analytisch.
- `-qm.vwn_variant vwn5` (Default) | `vwn_rpa` (Alias `vwn3`); die Variante steht in der Ergebniszeile
  (`LDA SVWN5 (Slater exchange + VWN5 correlation)/def2-SVP grid sg1`).
- V_mn = sum_g w_g v_xc(rho_g) phi_m phi_n blockweise (128 Punkte, eine GEMM pro Block), OpenMP über eine feste,
  nur von nbf abhängige Chunk-Zahl: Ergebnis unabhängig von `-qm.threads` (1 vs 4 Threads bitgleich).
- `-qm.xc_density_threshold` (Default 1e-12): Punkte mit rho darunter fehlen in E_xc und V_xc.
- Gitter + `AOEvaluator` einmal pro Geometrie (`ensureXCGrid`), verworfen in `buildOneElectronIntegrals`,
  `UpdateMolecule`, `resetForNewMolecule`; V_xc-Cache an die Dichte gebunden (bitweiser Vergleich).
- J-Pfade: gespeicherter Tensor nur `buildCoulomb`; RI (`-qm.ri j|jk`) nur J aus dem Fit (jk = J mit JKFIT-Hilfsbasis);
  Device-Tensor und integral-direct (CPU/GPU) berechnen K mit und verwerfen es (kein J-only-Direct-Kernel).
  Device-SCF-Schleife (`-qm.gpu_scf`) läuft für LDA nicht (kein V_xc auf dem Device), J/K dürfen vom Device kommen.
- Gradient: seit WP8 DFT (Okt 2026) analytisch, `-opt` funktioniert ([WP8](WP8_analytischer_gradient.md)).
- WP6/WP7-Schnittstelle: `XCFunctional::exact_exchange` (a_x, Fock -a_x/2 K, Energie -a_x/4 Tr(PK)), `needsGradient()`;
  seit WP6 ist das Funktional eine gewichtete `XCKernel`-Liste, der GGA-Zweig von `buildVxc` existiert ([WP6](WP6_pbe_gga.md)).

## Welche VWN wer benutzt (geprüft Okt 2026)

| Name | Bedeutung für zeta = 0 |
|---|---|
| curcuma `vwn5`, libxc `LDA_C_VWN` (id 7), PySCF `'LDA,VWN'`/`SVWN`/`VWN5`, ORCA `LSD`/`VWN5`, ORCA `B3LYP` | Monte-Carlo-Fit (VWN V) |
| curcuma `vwn_rpa`, libxc `LDA_C_VWN_RPA` (id 8), PySCF `VWN3`, ORCA `VWN3`, ORCA `B3LYP/G`, PySCF/libxc `B3LYP` (id 402) | RPA-Fit |
| libxc `LDA_C_VWN_1` ... `_4` | geschlossenschalig = VWN5 (Falle: `LDA_C_VWN_3` ist nicht RPA) |
| PySCF `LDA` allein | nur Slater, **ohne** Korrelation |

ORCA-Ausgabe kontrolliert (`Correlation Functional .... VWN-5` / `VWN-3`, B3LYP `LDAOpt VWN-5`, B3LYP/G `VWN-3`).
ORCA VWN3 = RPA-Fit belegt: H2O, feinstes Gitter, curcuma `vwn_rpa` - ORCA `VWN3` = 5.5e-8 Eh, derselbe Rest wie bei VWN5.
xcDFTs LDA-Korrelation (Gewicht im Nenner) ist kein Dichtefunktional und wurde nicht portiert.

## Validierung (Okt 2026, def2-SVP sphärisch, SCF 1e-10 / dE 1e-12 / |g| 1e-9)

Kernels gegen libxc 7.0.0 (`eval_xc`, 61 Dichten 1e-12 ... 1e2): Slater eps/v <= 4.1e-16 relativ, VWN5 <= 6.1e-14,
VWN-RPA <= 1.3e-14; v gegen FD von rho eps <= 3.4e-11.

Gegen PySCF RKS mit identischem Gitter (Treutler x Lebedev, Becke, Treutler-Radienanpassung, `prune=None`,
`small_rho_cutoff=0`):

| Fall | E curcuma (Eh) | dE vs PySCF | dE_xc | dN_grid |
|---|---|---|---|---|
| H2O sg1 | -75.7865448343 | 5.7e-14 | -2.1e-10 | 0 |
| CH4 sg1 | -40.0680638858 | 3.6e-14 | -2.9e-10 | 4e-14 |
| NH3 sg1 | -56.0221647353 | 1.7e-13 | 3.2e-10 | 3e-14 |
| HF sg1 | -99.7015950889 | -8.5e-14 | -2.6e-10 | 1e-14 |
| LiH sg1 | -7.9095485372 | 4.3e-14 | -2.6e-10 | -2e-15 |
| Ne sg1 | -128.0530990613 | 5.7e-14 | 1.3e-10 | -4e-15 |
| Benzol sg1 | -229.9291566382 | -5.1e-13 | 1.4e-10 | 2e-13 |
| H2O sg3 | -75.7865464797 | 2.8e-14 | -2.1e-10 | 5e-15 |
| H2O VWN-RPA sg1 | -75.9822532279 | 7.1e-14 | 2.7e-10 | 0 |
| NH3 VWN-RPA sg2 | -56.2141489235 | 3.0e-13 | -5.9e-10 | 1e-14 |
| H2O RI-J (JFIT) sg1 | -75.7866483143 | 1.3e-13 | -2.6e-10 | -1e-14 |
| Benzol RI-JK-Hilfsbasis sg1 | -229.9292414002 | -4.3e-12 | 1.4e-10 | 2e-13 |

dE_xc ~1e-10: E_xc ist erster Ordnung in der Dichte, die beiden SCFs enden an verschiedenen Dichten innerhalb ihrer
Schwellen; die Gesamtenergie ist stationär.

Gegen ORCA 6.1 (`! LSD def2-SVP NoRI VeryTightSCF`, Gitter nicht identisch, ORCA prunt; „fein" =
`%method IntAcc 8 AngularGrid 7 GridPruning 0`, IntAcc 12 ändert H2O/HF/Benzol um <= 1.1e-9):

| Molekül | ORCA DefGrid1 | DefGrid2 | DefGrid3 | ORCA fein | curcuma 250x1454 | curcuma - ORCA fein |
|---|---|---|---|---|---|---|
| H2O | -75.786513805 | -75.786544587 | -75.786546415 | -75.786546540 | -75.786546485 | 5.4e-8 |
| CH4 | -40.068064816 | -40.068064280 | -40.068064247 | -40.068064313 | -40.068064304 | 8.9e-9 |
| NH3 | -56.022157748 | -56.022172968 | -56.022172169 | -56.022172271 | -56.022172276 | -5.3e-9 |
| HF | -99.701539093 | -99.701594168 | -99.701594815 | -99.701595114 | -99.701595080 | 3.4e-8 |
| LiH | -7.909563182 | -7.909546887 | -7.909549490 | -7.909549565 | -7.909549573 | -8.1e-9 |
| Ne | -128.053097953 | -128.053098819 | -128.053099107 | -128.053099101 | -128.053099101 | -1.4e-10 |
| Benzol | -229.929159268 | -229.929123707 | -229.929123199 | -229.929124202 | -229.929123805 | 4.0e-7 |

Rest <= 5.4e-8 Eh für die kleinen Moleküle, 4.0e-7 für Benzol. Nicht belegt, woher er kommt: ORCAs feinstes hier
benutztes Winkelgitter ist Lebedev-770 (H eine Stufe kleiner), curcuma ändert sich zwischen 590 und 974 Winkelpunkten
für Benzol um 4e-7. Als Gitterunterschied vermutet, nicht bewiesen. WP6-Nachtrag: curcuma 250x590 = -229.929124204,
2e-9 von ORCA fein; ORCA gibt H nur Lebedev-590, AngularGrid 8/9 lehnt ORCA ab. Konsistent mit ORCAs Winkelauflösung.

Gitterkonvergenz (curcuma; je Stufe dE_xc / dE_total in Eh als Differenz zu 250x1454):

| Molekül | sg0 | sg1 | sg2 | sg3 | 150x974 | E_xc 250x1454 |
|---|---|---|---|---|---|---|
| H2O | -7.8e-5 / -7.4e-5 | 1.1e-6 / 1.7e-6 | -7.3e-7 / -5.4e-7 | 6.7e-9 / 5.6e-9 | 1.0e-9 / 1.3e-9 | -8.8251333457 |
| CH4 | -1.1e-4 / -9.3e-5 | -1.8e-6 / 4.2e-7 | -1.6e-7 / -1.8e-7 | 2.7e-9 / -3.7e-9 | 2.0e-10 / 1.0e-10 | -6.4557100861 |
| NH3 | -3.8e-5 / -3.9e-5 | 5.4e-6 / 7.5e-6 | -8.9e-8 / 1.1e-7 | -9.7e-9 / -9.0e-9 | -1.0e-10 / -3.0e-10 | -7.4784543158 |
| HF | -1.7e-4 / -1.8e-4 | -5.5e-7 / -8.9e-9 | 1.4e-7 / 8.1e-8 | 1.0e-8 / 8.6e-9 | -6.0e-10 / -7.0e-10 | -10.1901934200 |
| LiH | -2.0e-4 / -9.4e-5 | 2.6e-6 / 1.0e-6 | -2.4e-7 / -1.1e-7 | 1.5e-8 / 8.5e-9 | -2.0e-10 / -2.0e-10 | -2.0564802497 |
| Ne | -6.2e-4 / -5.8e-4 | 4.2e-8 / 4.0e-8 | 0 / 0 | 0 / 0 | 0 / 0 | -11.7553660966 |
| Benzol | -3.2e-4 / -3.1e-4 | -1.6e-5 / -3.3e-5 | 1.3e-6 / -3.4e-6 | -3.9e-7 / -4.0e-7 | -8.2e-9 / 1.3e-8 | -32.5420331248 |

Invarianz (Gesamtenergie, |dE| in Eh): umgekehrte Atomreihenfolge und Translation <= 7e-14 (H2O, NH3). Rotation
(Lebedev-Gitter raumfest, 5 zufällige Rotationen): H2O sg1 3.4e-6, sg2 8.3e-7, sg3 9.0e-9, 150x974 1.5e-9;
NH3 sg1 1.2e-5, sg3 2.3e-8; Benzol sg1 3.7e-5, sg2 6.7e-6, sg3 6.5e-7, 150x974 5.4e-8.
`UpdateMolecule` A -> B = frisches B auf 3e-14 (Gitter-/V_xc-Cache folgt der Geometrie). Integral-direct = Tensor
auf 1e-13; `-gpu cuda` (RTX 5090, mit und ohne direct) = CPU auf 10 Stellen (Benzol, einmal geprüft, kein ctest).
Laufzeit Benzol sg1: 10 V_xc-Builds 1.75 s (1 Thread) / 0.23 s (8 Threads).

## Nicht getestet / nicht implementiert

- Offene Schalen, geladene Moleküle, Z > 10 (Basis), andere Basen als def2-SVP, kartesische d nur im Gittertest (WP4).
- V_xc auf der GPU; Gitter-Pruning und Screening der AOs pro Block.
- Kein J-only-Kernel für integral-direct (K wird berechnet und verworfen); RI mixed precision mit LDA läuft über den
  vollen J/K-Build (nicht gemessen).
- ORCA-Rest (<= 5.4e-8 Eh) nicht aufgeklärt, Benzol 4e-7 nur plausibel (Winkelgitter, s. o.); xcDFT-Vergleich nicht möglich (Quellen fehlen in dieser Sitzung).
- Große Systeme: Speicher der Chunk-Matrizen (bis 128 MB) und Laufzeit nur bis Benzol gemessen.

## Quellen

W. Kohn, L. J. Sham, Phys. Rev. 140, A1133 (1965); P. A. M. Dirac, Proc. Cambridge Philos. Soc. 26, 376 (1930);
J. C. Slater, Phys. Rev. 81, 385 (1951); S. H. Vosko, L. Wilk, M. Nusair, Can. J. Phys. 58, 1200 (1980);
libxc `lda_c_vwn.c` (Parameter); Gitter: siehe WP4. Zitiert über `qmcite::lda` (`qm_citations.h`).
