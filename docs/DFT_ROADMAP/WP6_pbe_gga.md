# WP6 — PBE-GGA (+ ∇ρ, rung 2)

- **Status:** 🤖 AI-generated, ⚙️ machine-tested (`ctest -L qm_pbe` 21/21, Okt 2026); menschlicher Produktionstest ausstehend
- **Abhängigkeit:** WP5 (KS-SCF, V_xc-Gerüst), WP4 (Gitter, AO-Ableitungen)
- **Code:** `src/core/energy_calculators/qm_methods/qm_xc.{h,cpp}` (Kerne `pbeExchange`, `pbeCorrelation`,
  `pw92Correlation`, GGA-Zweig in `buildVxc`), Test `test_cases/qm_pbe/`, Referenz `scripts/qm_pbe_pyscf_reference.py`
  (PySCF 2.14, libxc 7.0.0)

## Was gebaut ist

- `-method pbe` = geschlossenschaliges RKS mit PBE-Austausch + PBE-Korrelation (libxc `GGA_X_PBE,GGA_C_PBE`, PySCF `PBE`),
  a_x = 0, kein K-Build; SCF, Caches, J-Pfade (gespeichert, direct, RI) wie WP5.
- PBE 1996: F_x(s) = 1 + κ - κ/(1 + μs²/κ), κ = 0.804, μ = 0.2195149727645171 (= βπ²/3); Korrelation PW92 + H(r_s, t),
  β = 0.06672455060314922, γ = (1 - ln 2)/π²; LDA-Teil PW92 mit libxc-`LDA_C_PW_MOD`-Konstanten (A = 0.0310907), nicht VWN.
- ε, v_ρ = ∂(ρε)/∂ρ, v_σ = ∂(ρε)/∂σ analytisch (σ = |∇ρ|², Gesamtdichte); Herleitung in `qm_xc.cpp`.
- V_mn = Σ_g w_g [v_ρ φ_m φ_n + 2 v_σ ∇ρ·∇(φ_m φ_n)], gebaut als X = φᵀ[(w v_ρ/2) φ + 2 w v_σ (∇ρ·∇φ)], V = X + Xᵀ
  (symmetrisch per Konstruktion); AO-Werte und -Ableitungen (deriv 1) pro 128-Punkt-Block, feste Chunk-Zahl wie LDA:
  1 vs 4 Threads bitgleich.
- WP7-Schnittstelle: `XCFunctional::terms` = gewichtete Liste von `XCKernel` (Slater, VWN5, VWN-RPA, PW92, PBE-X, PBE-C)
  plus `exact_exchange`; B3LYP braucht nur B88- und LYP-Kerne und die Liste 0.08 Slater + 0.72 B88 + 0.19 VWN-RPA + 0.81 LYP,
  a_x = 0.2. `needsGradient()` schaltet den GGA-Pfad ein, sobald ein GGA-Kern mit Gewicht ≠ 0 dabei ist.
- Gradient: seit WP8 DFT (Okt 2026) analytisch, `-opt` funktioniert ([WP8](WP8_analytischer_gradient.md)).

## Schwellen (wie hier gehalten)

- libxc 7.0.0 (gemessen über PySCF): PBE-X verwirft ρ < 1e-15, PBE-C ρ < 1e-12; σ wird auf 1e-40 (X) bzw. 1e-32 (C) angehoben.
- curcuma: `-qm.xc_density_threshold` (Default 1e-12) gilt für alle Terme; σ wird nicht angehoben (σ geht nur als σ ein, nie
  als √σ, σ = 0 ist regulär). Zusätzlich fester Boden 1e-20 nur für GGA-Kerne: bei ρ ≈ 1e-27 rundet ln(1 + 1/Q1) auf 0, A = ∞,
  NaN. Benzol sg1: Schwelle 0 / 1e-15 / 1e-12 ändern E um ≤ 4e-13 Eh.
- σ = 0 ohne libxc-Boden gegen libxc: X 5.7e-11, C 6.4e-3 relativ bei ρ = 1e-12 (Wirkung von libxcs Boden, kein Fehler).
- Große t (s → ∞): ε_c = ε_c^PW92 + H → 0 durch Auslöschung; Wert dort nur auf ~1e-16·|ε_c^PW92| genau (libxc ebenso,
  libxc liefert dort teils v_σ < 0 oder 0, curcuma positive Werte). Der Kernel-Test vergleicht deshalb mit skaliertem Fehler.

## Validierung (Okt 2026, def2-SVP sphärisch, SCF 1e-10 / dE 1e-12 / |g| 1e-9)

Kerne gegen libxc (`eval_xc`, 464 Paare: ρ 1e-12..1e2 × σ ∈ {0, 1e-30..1e10} ∪ s ∈ {0.01..1e3}; skaliert = |Δ|/max(|Wert|,
Größe der auslöschenden Terme), Definition in `check_qm_pbe.py`):

| Kern | ε | v_ρ | v_σ | v_ρ vs FD | v_σ vs FD |
|---|---|---|---|---|---|
| PBE-X (`GGA_X_PBE`) | 6.6e-16 | 8.6e-16 | 7.4e-16 (rein rel. 1.4e-15) | 1.0e-12 | 3.0e-13 |
| PBE-C (`GGA_C_PBE`) | 5.1e-15 | 3.7e-15 | 9.6e-16 | 9.7e-8 | 3.7e-11 |
| PW92 (`LDA_C_PW_MOD`) | 5.9e-16 | 5.7e-16 | – | 9.7e-8 | – |

FD-Rest 1e-7 bei ρ = 1e-12: ln(1 + x) statt log1p (wie libxc) verliert dort ~1e-10 relativ, die Differenz verstärkt es.
V_xc gegen FD von E_xc in der Dichtematrix (4-Punkt, h = 1e-4): H2O 1.6e-10, Benzol 7.8e-9 relativ; V 1 vs 4 Threads 0.

Gegen PySCF RKS mit identischem Gitter (Treutler x Lebedev, Becke, Treutler-Radien, `prune=None`, `small_rho_cutoff=0`):

| Fall | E curcuma (Eh) | dE vs PySCF | dE_xc | dN_grid |
|---|---|---|---|---|
| H2O sg1 | -76.2630324241 | 9.9e-14 | -2.7e-10 | 4e-15 |
| CH4 sg1 | -40.4143231313 | 1.6e-13 | -4.4e-10 | 6e-14 |
| NH3 sg1 | -56.4309554971 | 1.3e-13 | -4.4e-10 | 2e-14 |
| HF sg1 | -100.2538446404 | 5.7e-14 | -3.4e-10 | 7e-15 |
| LiH sg1 | -8.0378769530 | 2.7e-14 | 3.6e-10 | 9e-16 |
| Ne sg1 | -128.6878678656 | -5.7e-14 | -7.8e-11 | -2e-15 |
| Benzol sg1 | -231.7711867597 | -1.3e-12 | 5.7e-11 | 9e-14 |
| H2O sg3 | -76.2630349660 | 1.4e-13 | -2.7e-10 | -1e-14 |
| NH3 sg2 | -56.4309626598 | 4.3e-14 | -3.9e-10 | 3e-14 |
| H2O RI-J (JFIT) sg1 | -76.2631417480 | 1.4e-13 | -2.6e-10 | -9e-15 |
| Benzol RI-JK-Hilfsbasis sg1 | -231.7712747119 | -3.9e-12 | 2.6e-11 | 1e-13 |

dE_xc ~1e-10: wie WP5 (E_xc erster Ordnung in der Dichte, Gesamtenergie stationär).

Gegen ORCA 6.1 (`! PBE def2-SVP NoRI VeryTightSCF`; „fein" = `%method IntAcc 8 AngularGrid 7 GridPruning 0`).
ORCA-Ausgabe: `Exchange PBE, XKappa 0.804000, XMuePBE 0.219520`, `Correlation PBE, CBetaPBE 0.066725`,
`LDAOpt PW91-LDA` (= PW92). ORCAs Standard-μ ist 0.21952, libxc/PySCF/curcuma 0.2195149727645171; „fein, libxc-μ/β" setzt
`XMuePBE 0.2195149727645171 CBetaPBE 0.06672455060314922` (Ausgabe dann 0.219515):

| Molekül | ORCA DefGrid1 | DefGrid2 | DefGrid3 | ORCA fein | ORCA fein, libxc-μ/β | curcuma 250x1454 | - ORCA fein | - ORCA fein libxc-μ/β |
|---|---|---|---|---|---|---|---|---|
| H2O | -76.263016736 | -76.263046470 | -76.263049296 | -76.263049151 | -76.263035052 | -76.263034982 | 1.4e-5 | 7.0e-8 |
| CH4 | -40.414336122 | -40.414334587 | -40.414334284 | -40.414334243 | -40.414323677 | -40.414323682 | 1.1e-5 | -5.0e-9 |
| NH3 | -56.430964947 | -56.430975864 | -56.430975000 | -56.430975042 | -56.430962853 | -56.430962867 | 1.2e-5 | -1.5e-8 |
| HF | -100.253811741 | -100.253860329 | -100.253861260 | -100.253861743 | -100.253845565 | -100.253845568 | 1.6e-5 | -3.9e-9 |
| LiH | -8.037889227 | -8.037877511 | -8.037881531 | -8.037881555 | -8.037877511 | -8.037877518 | 4.0e-6 | -7.1e-9 |
| Ne | -128.687894340 | -128.687889292 | -128.687885506 | -128.687885436 | -128.687867042 | -128.687867043 | 1.8e-5 | -1.1e-10 |
| Benzol | -231.771231028 | -231.771232802 | -231.771234175 | -231.771233440 | -231.771179532 | -231.771178521 | 5.5e-5 | 1.0e-6 |

- 4e-6..5.5e-5 Eh gegen ORCAs Standard-PBE sind das μ (Erstordnungs-Abschätzung H2O mit PySCF-Dichte: μ = 0.21952 statt
  libxc-μ gibt -1.5e-5 Eh, gemessen 1.4e-5). Mit libxc-μ/β in ORCA bleiben ≤ 7e-8 Eh (kleine Moleküle).
- Benzol-Rest 1e-6 (und WP5-LDA 4e-7): ORCAs feinstes Winkelgitter ist Lebedev-770, H eine Stufe kleiner (590);
  AngularGrid 8/9 lehnt ORCA ab. curcuma 250x590 = -231.771179507 (ORCA libxc-μ/β: -231.771179532, Δ 2.5e-8), 250x770
  -231.771178609, 250x1202 -231.771178564, 250x2030 -231.771178509, 300x2702 -231.771178504. LDA: curcuma 250x590
  -229.929124204 gegen ORCA LSD fein -229.929124202 (Δ 2e-9). Der Rest ist damit mit ORCAs Winkelauflösung konsistent;
  ORCAs Gitter wurde nicht nachgebaut, also nicht bewiesen.
- ORCA IntAcc 12 statt 8 ändert ≤ 9e-9 Eh.

Gitterkonvergenz (curcuma; je Stufe dE_xc / dE_total in Eh als Differenz zu 250x1454):

| Molekül | sg0 | sg1 | sg2 | sg3 | 150x974 | E_xc 250x1454 |
|---|---|---|---|---|---|---|
| H2O | -4.8e-4 / -6.0e-4 | 2.2e-6 / 2.6e-6 | -2.7e-7 / -7.1e-8 | 2.6e-8 / 1.6e-8 | -5.3e-9 / -1.1e-8 | -9.3133758406 |
| CH4 | 2.3e-4 / 2.6e-4 | -1.9e-6 / 5.5e-7 | -7.1e-8 / -9.7e-8 | 3.8e-9 / -5.8e-9 | -1.0e-10 / -5.2e-10 | -6.8209653692 |
| NH3 | 3.2e-4 / 3.5e-4 | 4.7e-6 / 7.4e-6 | -1.0e-8 / 2.1e-7 | -1.1e-9 / -2.6e-9 | 7.1e-9 / 5.1e-9 | -7.9009654053 |
| HF | 4.3e-4 / 3.7e-4 | 6.4e-7 / 9.3e-7 | -6.3e-8 / -1.2e-7 | 2.8e-8 / 1.9e-8 | 1.4e-9 / 3.4e-9 | -10.7548991416 |
| LiH | -2.5e-4 / -1.0e-4 | 1.3e-6 / 5.6e-7 | -2.9e-7 / -1.4e-7 | 1.8e-8 / 6.0e-9 | -1.3e-8 / -2.2e-10 | -2.1937145142 |
| Ne | -1.7e-3 / -1.2e-3 | -3.5e-8 / -8.2e-7 | 1.6e-9 / 3.0e-9 | -2.5e-11 / -4.6e-12 | -1.1e-13 / -8.5e-14 | -12.4062138478 |
| Benzol | 1.7e-3 / 1.8e-3 | 2.7e-5 / -8.2e-6 | 2.5e-6 / -1.6e-6 | -1.6e-6 / -9.8e-7 | 3.1e-7 / 1.6e-7 | -34.4470306083 |

sg0 ist für PBE 3-10x schlechter als für LDA (bis 1.8e-3 Eh); sg1→sg2 < 1e-4 Eh (Ziel erfüllt).

Invarianz (Gesamtenergie, |dE| in Eh; Rotation: Maximum über 5 zufällige Rotationen, das Lebedev-Gitter bleibt raumfest):

| Molekül | Gitter | umgekehrte Reihenfolge | Translation | Rotation |
|---|---|---|---|---|
| H2O | sg1 / sg2 / sg3 / 150x974 | ≤ 7e-14 | ≤ 7e-14 | 9.9e-6 / 2.8e-7 / 3.2e-8 / 1.7e-8 |
| NH3 | sg1 / sg2 / sg3 / 150x974 | ≤ 2e-13 | ≤ 1.3e-13 | 1.1e-5 / 1.6e-6 / 2.8e-8 / 8.7e-9 |
| Benzol | sg1 / sg2 / sg3 / 150x974 | ≤ 5e-13 | ≤ 3e-13 | 1.9e-5 / 4.9e-6 / 1.5e-6 / 2.3e-7 |

Rotationsschwankung PBE gegenüber LDA (WP5): H2O sg1 9.9e-6 statt 3.4e-6, Benzol sg1 1.9e-5 statt 3.7e-5, sg3 1.5e-6 statt 6.5e-7.
`UpdateMolecule` A → B = frisches B (Differenz 0); integral-direct = Tensor auf 6e-14; `-gpu cuda` (RTX 5090, Tensor und
`-qm.scf_direct on`) = CPU auf 8 ausgegebenen Stellen (Benzol, einmal geprüft, kein ctest).

## Nicht getestet / nicht implementiert

- Offene Schalen (ζ ≠ 0, φ ≠ 1), geladene Moleküle, Z > 10, andere Basen als def2-SVP, kartesische d mit PBE.
- meta-GGA; V_xc auf der GPU (`-gpu cuda` liefert nur J), `-qm.gpu_scf` läuft für KS nicht.
- Kein AO-Screening pro Block, kein Gitter-Pruning. Laufzeit Benzol sg1, 10 V_xc-Builds: PBE 2.46 s (1 Thread) / 0.33 s (8),
  LDA 1.75 / 0.24 s.
- ORCA mit identischem Gitter nicht reproduziert; Benzol-Rest 1e-6 nur plausibel erklärt (Winkelgitter).
- Große Systeme: nur bis Benzol gemessen.

## Quellen

J. P. Perdew, K. Burke, M. Ernzerhof, Phys. Rev. Lett. 77, 3865 (1996), DOI 10.1103/PhysRevLett.77.3865;
J. P. Perdew, Y. Wang, Phys. Rev. B 45, 13244 (1992), DOI 10.1103/PhysRevB.45.13244 (DOIs gegen Crossref geprüft);
libxc `gga_x_pbe.c`, `gga_c_pbe.c`, `lda_c_pw.c` (Konstanten). Zitiert über `qmcite::kohnSham` (`qm_citations.h`).
