# WP7 — B3LYP-Hybrid (+ exakter Austausch, rung 4)

- **Status:** 🤖 AI-generated, ⚙️ machine-tested (`ctest -L qm_b3lyp` 29/29, Okt 2026); menschlicher Produktionstest ausstehend
- **Abhängigkeit:** WP6 (GGA-Zweig von `buildVxc`), WP5 (VWN), WP2/WP10 (K-Build)
- **Code:** `src/core/energy_calculators/qm_methods/qm_xc.{h,cpp}` (Kerne `b88Exchange`, `lypCorrelation`, Tabelle
  `b3lypFunctional`), Auswahl in `qm_engine.cpp`; Test `test_cases/qm_b3lyp/`, Referenz
  `scripts/qm_b3lyp_pyscf_reference.py` (PySCF 2.14, libxc 7.0.0)

## Was gebaut ist

- `-method b3lyp` = geschlossenschaliges RKS, E_xc = 0.08 Slater + 0.72 B88 + 0.20 HF + 0.19 VWN + 0.81 LYP
  (Stephens et al. 1994, Becke 1993); a_x = 0.2 in `XCFunctional::exact_exchange`, sonst nur eine Tabellenzeile.
- F = h + J + V_xc - (a_x/2) K, E = Tr(Ph) + 1/2 Tr(PJ) - (a_x/4) Tr(PK) + E_xc + E_nn (spin-summierte Konvention von
  `composeFock`/`scfEnergy` in `qm_scf.cpp`, seit WP5 vorbereitet, für WP7 nachgerechnet, nicht geändert).
- B88 (Becke 1988 Gl. 8) und LYP (Lee-Yang-Parr 1988 in der Laplace-freien Form von Miehlich et al. 1989 Gl. 2, wie
  libxc `GGA_C_LYP`): ε, v_ρ, v_σ analytisch, geschlossenschalige Reduktion und Ableitungen in `qm_xc.cpp`.
- **VWN-Variante** über `-qm.vwn_variant` (wie `lda`), Default `vwn5`; die Ergebniszeile nennt sie immer
  (`B3LYP-VWN5 (...)` bzw. `B3LYP-VWN-RPA (...)`). Kein eigener Methodenname `b3lyp-g` (MethodFactory kennt keine
  Varianten-Namen für QM; ein Parameter genügt und ist dort, wo `lda` ihn schon hat).

| curcuma | entspricht |
|---|---|
| `-method b3lyp` (Default `vwn5`) | ORCA `B3LYP` (Ausgabe `LDAOpt VWN-5`), PySCF `B3LYP5`, libxc `HYB_GGA_XC_B3LYP5` (475) |
| `-method b3lyp -qm.vwn_variant vwn_rpa` | Gaussian `B3LYP`, libxc/PySCF `B3LYP` = `HYB_GGA_XC_B3LYP` (402), ORCA `B3LYP/G` (`VWN-3`) |

Begründung Default VWN5: gleicher Default wie `-method lda`, gleicher Wert wie ORCA (zweite Referenz des Projekts);
Gaussian-Nutzer müssen `vwn_rpa` setzen. Unterschied der Varianten: 0.037 Eh (H2O) bis 0.15 Eh (Benzol).
PySCF-Falle: der Alias `B3LYP` hat zwischen PySCF-Versionen seine Bedeutung gewechselt; das Referenzskript nutzt die
libxc-IDs.

- K-Pfade: gespeicherter Tensor, integral-direct (CPU, GPU), Device-Tensor: wie HF (`buildJK`), V_xc auf der CPU.
  `-qm.ri jk`: J und K gefittet (wie HF). `-qm.ri j`: J gefittet, K exakt (PySCF `only_dfj=True`, ORCA „RIJONX“),
  Ausgabe `RI-J (exact K)`; das ist eine gültige Variante, keine stille Näherung. Device-SCF-Schleife läuft für KS nicht.
- Gradient: seit WP8 DFT (Okt 2026) analytisch, `-opt` funktioniert ([WP8](WP8_analytischer_gradient.md)).
- PBE0 o. ä. wäre eine weitere Funktion wie `b3lypFunctional` (0.75 PBE-X + PBE-C, a_x 0.25); nicht angelegt, da
  nicht validiert (ORCA bräuchte wieder libxc-μ, s. WP6).

## Schwellen und Sonderfälle

- libxc 7.0.0 (über PySCF gemessen): `GGA_X_B88` und `GGA_C_LYP` verwerfen ρ < 1e-15; σ wird auf
  (1e-15)^(8/3) = 1e-40 angehoben (aus der dens_threshold abgeleitet, wie bei PBE gemessen; direkt messen geht für LYP
  nicht, dessen σ-Term ~exp(-c ρ^(-1/3)) bei kleinem ρ verschwindet).
- curcuma: `-qm.xc_density_threshold` (Default 1e-12) plus fester GGA-Boden 1e-20 (WP6). B88 leitet nach y = x² ab, nicht
  nach x = √y: σ = 0 ist regulär (g'(0) = 1); x/√(1+x²) für große x als 1/√(1+1/y). LYP: exp(-c r) wird unterhalb
  ρ ≈ 4e-11 zu 0, dann bleibt nur der σ-freie Term; r⁵ bleibt oberhalb 1e-20 endlich.
- NaN/Inf-Scan (ctest): ρ 1e-20..1e2 × σ 0..1e20, 10858 Kern-Auswertungen, 0 nicht endlich.
- σ = 0 ohne libxc-Boden gegen libxc: B88 7.2e-11 relativ bei ρ = 1e-12 (Wirkung von libxcs Boden), LYP 1.6e-16.

## Validierung (Okt 2026, def2-SVP sphärisch, SCF 1e-10 / dE 1e-12 / |g| 1e-9)

Kerne gegen libxc (`eval_xc`, 472 Paare: ρ 1e-12..1e2 × σ ∈ {0, 1e-30..1e10} ∪ s ∈ {0.01..1e3} ∪ σ = 8e12; skaliert =
|Δ| / max(|Wert|, |σ-freier Teil|), für v_σ max(|Wert|, ρ|ε₀|/σ); Definition in `check_qm_b3lyp.py`):

| Kern | ε | v_ρ | v_σ | v_ρ vs FD | v_σ vs FD |
|---|---|---|---|---|---|
| B88 (`GGA_X_B88`) | 6.7e-16 | 3.2e-14 (ρ 1e-9, σ 1e10) | 1.5e-15 | 9.5e-12 | 1.0e-12 |
| LYP (`GGA_C_LYP`) | 2.6e-15 | 6.9e-15 | 1.1e-15 | 1.0e-9 | 2.4e-13 |

V_xc gegen FD von E_xc in der Dichtematrix (4-Punkt, h = 1e-4): H2O 2.4e-10, Benzol 1.5e-9 relativ; 1 vs 4 Threads 0.

Gegen PySCF RKS mit identischem Gitter (Treutler × Lebedev, Becke, Treutler-Radien, `prune=None`,
`small_rho_cutoff=0`); E_HFX = -(a_x/4) Tr(PK) gegen PySCFs `get_k` auf dessen Dichte:

| Fall | E curcuma (Eh) | dE vs PySCF | dE_xc | dE_HFX | Iter. |
|---|---|---|---|---|---|
| H2O VWN5 sg1 | -76.3134099519 | 7.5e-12 | -1.2e-10 | -3.3e-11 | 12 |
| CH4 VWN5 sg1 | -40.4514930294 | 1.9e-11 | 2.5e-10 | 5.1e-11 | 12 |
| NH3 VWN5 sg1 | -56.4692896643 | 1.8e-11 | -2.0e-10 | -5.4e-11 | 13 |
| HF VWN5 sg1 | -100.3154002576 | 6.0e-12 | 1.5e-10 | 3.6e-11 | 11 |
| LiH VWN5 sg1 | -8.0655391033 | 3.2e-12 | -2.1e-12 | -1.5e-12 | 12 |
| Ne VWN5 sg1 | -128.7647367078 | 7.1e-12 | -2.4e-10 | -6.2e-11 | 6 |
| Benzol VWN5 sg1 | -231.9293723520 | 3.0e-11 | -3.9e-10 | -9.9e-11 | 14 |
| H2O VWN-RPA sg1 | -76.3506035075 | 7.5e-12 | -1.3e-10 | -3.4e-11 | 12 |
| CH4 VWN-RPA sg1 | -40.4874822265 | 1.8e-11 | 2.7e-10 | 5.5e-11 | 12 |
| NH3 VWN-RPA sg1 | -56.5057775542 | 1.8e-11 | -2.0e-10 | -5.2e-11 | 13 |
| HF VWN-RPA sg1 | -100.3531561831 | 6.1e-12 | 9.7e-11 | 2.2e-11 | 11 |
| LiH VWN-RPA sg1 | -8.0794797662 | 3.0e-12 | -1.0e-11 | -3.6e-12 | 12 |
| Ne VWN-RPA sg1 | -128.8030841596 | 7.0e-12 | -2.2e-10 | -5.6e-11 | 6 |
| Benzol VWN-RPA sg1 | -232.0837055532 | 3.0e-11 | -4.1e-10 | -1.0e-10 | 14 |
| H2O VWN5 sg3 | -76.3134120218 | 6.1e-12 | -1.2e-10 | -3.2e-11 | 12 |
| NH3 VWN-RPA sg2 | -56.5057835198 | 1.4e-11 | -2.0e-10 | -5.1e-11 | 13 |
| H2O VWN5 RI-J (JFIT, K exakt) | -76.3135216073 | 7.6e-12 | -3.3e-10 | -8.2e-11 | 12 |
| Benzol VWN5 RI-JK (JKFIT) | -231.9294218501 | 2.9e-11 | -2.6e-10 | -7.1e-11 | 14 |
| H2O VWN-RPA RI-JK (JKFIT) | -76.3506150727 | 7.5e-12 | -1.3e-10 | -3.4e-11 | 12 |

dE ≤ 3.0e-11 Eh; die Komponenten (erster Ordnung in der Dichte) bis 4e-10, die Gesamtenergie ist stationär (wie WP5/WP6).
Warum der Rest größer ist als bei PBE (≤ 3.9e-12), wurde nicht untersucht.

Gegen ORCA 6.1 (`! B3LYP def2-SVP NoRI VeryTightSCF` bzw. `B3LYP/G`; „fein" = `%method IntAcc 8 AngularGrid 7
GridPruning 0`). ORCA-Ausgabe geprüft: `Exchange B88, XAlpha 0.666667, XBeta 0.004200`, `Correlation LYP`,
`LDAOpt VWN-5` (B3LYP) bzw. `VWN-3` (B3LYP/G), `ScalHFX 0.200000, ScalDFX 0.720000, ScalDFC 0.810000, ScalLDAC 1.000000`
(ORCA zählt die GGA-Anteile relativ zum LDA-Teil: 0.08 = 1 - 0.2 - 0.72, 0.19 = 1 - 0.81); kein RIJCOSX, keine RI.
Keine abweichende Konstante gefunden (anders als PBE-μ in WP6); LYP-a..d druckt ORCA nicht.

| Molekül | ORCA DefGrid1 | DefGrid2 | DefGrid3 | ORCA fein | curcuma 250x1454 | curcuma - ORCA fein |
|---|---|---|---|---|---|---|
| **B3LYP (VWN5)** | | | | | | |
| H2O | -76.313383137 | -76.313412534 | -76.313411442 | -76.313412080 | -76.313412026 | 5.3e-8 |
| CH4 | -40.451494966 | -40.451493297 | -40.451493124 | -40.451493188 | -40.451493184 | 4.4e-9 |
| NH3 | -56.469284825 | -56.469296234 | -56.469295699 | -56.469295664 | -56.469295669 | -4.8e-9 |
| HF | -100.315368837 | -100.315401765 | -100.315399601 | -100.315400712 | -100.315400724 | -1.2e-8 |
| LiH | -8.065549618 | -8.065537083 | -8.065539354 | -8.065539299 | -8.065539306 | -6.8e-9 |
| Ne | -128.764740196 | -128.764736487 | -128.764736248 | -128.764736223 | -128.764736223 | -1.4e-10 |
| Benzol | -231.929372869 | -231.929390817 | -231.929382852 | -231.929383875 | -231.929384495 | -6.2e-7 (250x590: 3.6e-9) |
| **B3LYP/G (VWN-RPA)** | | | | | | |
| H2O | -76.350576486 | -76.350606092 | -76.350605009 | -76.350605647 | -76.350605593 | 5.3e-8 |
| CH4 | -40.487484203 | -40.487482522 | -40.487482348 | -40.487482412 | -40.487482408 | 4.5e-9 |
| NH3 | -56.505772641 | -56.505784212 | -56.505783679 | -56.505783644 | -56.505783649 | -4.7e-9 |
| HF | -100.353124585 | -100.353157688 | -100.353155530 | -100.353156641 | -100.353156653 | -1.2e-8 |
| LiH | -8.079490403 | -8.079477719 | -8.079480011 | -8.079479957 | -8.079479963 | -6.6e-9 |
| Ne | -128.803087647 | -128.803083939 | -128.803083700 | -128.803083675 | -128.803083675 | -1.7e-10 |
| Benzol | -232.083705937 | -232.083723747 | -232.083715696 | -232.083716721 | -232.083717338 | -6.2e-7 (250x590: 3.8e-9) |

- Kleine Moleküle ≤ 5.3e-8 Eh, für beide Varianten gleich und mit demselben Muster wie LDA (H2O 5.4e-8) und PBE (7e-8):
  der Rest hängt am Gitter, nicht am Funktional (vermutet, ORCAs Gitter nicht nachgebaut).
- Benzol 6.2e-7; mit curcuma 250x590 (ORCAs H-Winkelgitter ist Lebedev-590) 3.6e-9 / 3.8e-9, wie WP5/WP6.
- SCF-Iterationen curcuma 6-14 (Ziel < 80), ORCA 11 (H2O).

Gitterkonvergenz (curcuma VWN5; je Stufe dE_xc / dE_total in Eh als Differenz zu 250x1454):

| Molekül | sg0 | sg1 | sg2 | sg3 | 150x974 | E_xc 250x1454 |
|---|---|---|---|---|---|---|
| H2O | -2.5e-4 / -2.0e-4 | 2.0e-6 / 2.1e-6 | -2.8e-7 / -1.6e-7 | 9.0e-9 / 4.6e-9 | -5.1e-9 / -4.0e-9 | -7.5674171471 |
| CH4 | 2.4e-5 / -1.5e-5 | -1.8e-6 / 1.5e-7 | -2.4e-8 / -4.6e-8 | -2.2e-9 / -6.1e-9 | -1.0e-10 / -1.0e-10 | -5.5465358732 |
| NH3 | 8.9e-5 / 4.3e-5 | 4.2e-6 / 6.0e-6 | -3.2e-8 / 1.3e-7 | -4.8e-9 / -5.2e-9 | 1.6e-9 / 5.0e-10 | -6.4209283983 |
| HF | -4.8e-5 / -1.3e-4 | 2.8e-7 / 4.7e-7 | -2.5e-8 / -6.0e-8 | 1.0e-8 / 5.7e-9 | 1.2e-9 / 1.0e-9 | -8.7360966590 |
| LiH | -1.2e-4 / -4.7e-5 | 1.5e-6 / 2.0e-7 | -2.0e-7 / -1.1e-7 | 3.1e-9 / 8.1e-9 | -5.9e-9 / 3.0e-10 | -1.7990947575 |
| Ne | -7.3e-4 / -5.8e-4 | -4.6e-7 / -4.8e-7 | 0 / 1.0e-10 | 0 / 0 | 0 / 0 | -10.0703769021 |
| Benzol | 4.4e-4 / 2.1e-4 | 2.3e-5 / 1.2e-5 | -5.9e-6 / -6.8e-6 | 1.1e-6 / 6.2e-7 | -3.6e-7 / -1.9e-7 | -27.9773430712 |

Invarianz (Gesamtenergie, |dE| in Eh; Rotation: Maximum über 5 zufällige Rotationen, Lebedev-Gitter raumfest):

| Molekül | Gitter | umgekehrte Reihenfolge | Translation | Rotation |
|---|---|---|---|---|
| H2O | sg1 / sg2 / sg3 / 150x974 | ≤ 1.4e-14 | ≤ 7.1e-14 | 4.5e-6 / 4.0e-7 / 1.6e-8 / 7.0e-9 |
| NH3 | sg1 / sg2 / sg3 / 150x974 | ≤ 1.7e-13 | ≤ 1.3e-13 | 1.4e-5 / 4.7e-7 / 1.4e-8 / 8.6e-10 |
| Benzol | sg1 / sg2 / sg3 / 150x974 | ≤ 9.7e-13 | ≤ 2.3e-13 | 2.0e-5 / 1.3e-5 / 9.4e-7 / 2.0e-7 |

`UpdateMolecule` A → B = frisches B (2.8e-14); integral-direct = Tensor (H2O 4.3e-14, Benzol 6e-13); `-gpu cuda`
(RTX 5090, Device-Tensor und `-qm.scf_direct on` mit FP32-Anfangsphase) = CPU auf 8 ausgegebenen Stellen (Benzol, einmal
geprüft, kein ctest).

## Nicht getestet / nicht implementiert

- Offene Schalen, geladene Moleküle, Z > 10, andere Basen als def2-SVP, kartesische d mit B3LYP.
- andere Hybride (PBE0), range-separated Hybride, DFT-D (`b3lyp` ist reines B3LYP ohne Dispersion).
- V_xc auf der GPU, `-qm.gpu_scf` für KS, `-qm.ri_mixed_precision` und `-qm.ri_fit_radius` mit B3LYP (laufen über den
  HF-Pfad von `buildJK`, nicht gemessen), RI-JK `batched`.
- ORCA mit identischem Gitter nicht reproduziert; Rest ≤ 5.3e-8 Eh (klein) bzw. 6.2e-7 (Benzol) nur plausibel erklärt.
- LYP-Konstanten in ORCA nicht ausgegeben, nur über die Energie-Übereinstimmung indirekt geprüft.
- Große Systeme: nur bis Benzol gemessen.

## Quellen

A. D. Becke, Phys. Rev. A 38, 3098 (1988), DOI 10.1103/PhysRevA.38.3098; C. Lee, W. Yang, R. G. Parr, Phys. Rev. B 37,
785 (1988), DOI 10.1103/PhysRevB.37.785; B. Miehlich, A. Savin, H. Stoll, H. Preuss, Chem. Phys. Lett. 157, 200
(1989), DOI 10.1016/0009-2614(89)87234-3; A. D. Becke, J. Chem. Phys. 98, 5648 (1993), DOI 10.1063/1.464913;
P. J. Stephens, F. J. Devlin, C. F. Chabalowski, M. J. Frisch, J. Phys. Chem. 98, 11623 (1994),
DOI 10.1021/j100096a001 (alle DOIs gegen Crossref geprüft); libxc `gga_x_b88.c`, `gga_c_lyp.c`, `hyb_gga_xc_b3lyp.c`
(Konstanten). Zitiert über `qmcite::kohnSham` (`qm_citations.h`).
