# WP11 — Meta-GGA TPSS und SCAN (+ τ, rung 3)

- **Status Stufe 1 (Energien, CPU):** 🤖 AI-generated, ⚙️ machine-tested (`ctest -L qm_mgga` 50/50, 30. Okt 2026);
  **Stufe 2 (analytischer Gradient, CPU):** 🤖 AI-generated, ⚙️ machine-tested (`ctest -L qm_mgga_grad`, 30. Okt 2026,
  Abschnitt „Stufe 2“); **Stufe 3 (CUDA: V_xc, KS-SCF-Loop, XC-Gradient):** 🤖 AI-generated, ⚙️ machine-tested
  (`ctest -L qm_gpu_dft`, 30. Okt 2026, nur RTX 5090, Abschnitt „Stufe 3“); menschlicher Produktionstest ausstehend.
- Name: "WP11" steht in der README-Tabelle auch für das kanonische RMP2 (historisch); diese Datei ist der Meta-GGA-Teil.
- **Code:** Kerne `qm_xc_kernels.h` (Meta-GGA-Block: `evaluateKernelMGGA`, `evaluatePointMGGA`, Dual-Zahlen `Dual3`),
  Tabelle `kMetaGGATable` und `buildVxcMeta` in `qm_xc.cpp`, `QMFunctional::TPSS/SCAN` in `qm_engine.*`, Methodennamen in
  `method_factory.cpp`, Zitate `tpss`/`scan` (`qm_citations.h`, `citation_database.cpp`, DOIs gegen Crossref geprüft).
  Test `test_cases/qm_mgga/`, Referenz `scripts/qm_mgga_pyscf_reference.py` (PySCF 2.14, libxc 7.0.0).

## Was gebaut ist

- `-method tpss` (libxc `MGGA_X_TPSS,MGGA_C_TPSS`, ORCA `TPSS`/`LibXC(TPSS)`), `-method scan` (`MGGA_X_SCAN,MGGA_C_SCAN`,
  ORCA 6.1 nur `LibXC(SCAN)`, kein `! SCAN`-Schlüsselwort); geschlossene Schale, a_x = 0, H-Ne.
- τ = ½ Σ_i |∇ψ_i|² über beide Spins = ½ Σ_mn P_mn ∇φ_m·∇φ_n (Konvention von libxc unpolarisiert und PySCF `eval_rho`);
  pro Block aus AO-Ableitungen (deriv 1): τ_q = ½ Σ_k (∂_kφ P ∂_kφᵀ)_qq.
- V_mn = Σ_g w_g [v_ρ φ_mφ_n + 2 v_σ ∇ρ·∇(φ_mφ_n) + ½ v_τ ∇φ_m·∇φ_n] (∂τ/∂P_mn = ½ ∇φ_m·∇φ_n); eigene Funktion
  `buildVxcMeta` (mit und ohne AO-Screening, Chunks unabhängig von der Threadzahl). LDA/GGA/Hybrid-Pfade unberührt.
- Formeln: TPSS-PRL 91, 146401 Gl. 3-14 (Austausch Gl. 5, 7, 8, 10; Korrelation Gl. 11-14 mit C(0,0) = 0.53, ξ = 0
  bei geschlossener Schale, PBE-β konstant); SCAN nach libxc `mgga_x_scan.mpl`/`mgga_c_scan.mpl` (Gleichungsnummern der
  SCAN-PRL nicht geprüft, Paper hier nicht zugänglich). Kommentare mit Gleichungen in `qm_xc_kernels.h`.
- Ableitungen: Vorwärts-AD mit Dual-Zahlen (Wert + ∂/∂ρ, ∂/∂σ, ∂/∂τ), eine Quelle für ε und v. Grund: TPSS-C schachtelt
  drei PBE-Auswertungen, max() und min(), SCAN hat stückweise Interpolation. Der vorhandene analytische PBE-C-Kern geht
  unverändert ein (in einen Dual-Wert gehoben). Kein Näherungsverfahren: Ableitung des ausgewerteten Ausdrucks auf Rundung.
- Funktionaltabelle: eine Zeile pro Funktional (Name, Label, a_x, Terme); TPSSh wäre `{0.9 TPSS_X, TPSS_C, a_x 0.1}`,
  r2SCAN braucht eigene Kerne. Beides **nicht** implementiert.

## Schwellen (libxc 7.0.0, Quellen `functionals.c`, `work_mgga_inc.c`)

- libxc: dens_threshold 1e-15, σ ≥ (1e-15^(4/3))² = 1e-40, τ ≥ 1e-20, σ ≤ 8ρτ (τ ≥ τ_W), ζ-Schwelle DBL_EPSILON.
- curcuma (`evaluateKernelMGGA`): τ-Boden 1e-20 und σ ≤ 8ρτ wie libxc (Ableitungen am geklemmten Punkt, wie libxc);
  σ-Boden nicht (σ = 0 exakt behandelt: TPSS-Wurzel mit σ ausgeklammert, SCAN g_x = 1 für a1/p^(1/4) > 700);
  Dichteboden 1e-20 (GGA-Regel) plus `-qm.xc_density_threshold`. Für P positiv semidefinit gilt τ ≥ τ_W exakt.
- TPSS-C ist bei z = τ_W/τ = 1 nicht differenzierbar (min(z, 1), Gl. 11): bei α = 0 entscheidet die Rundung die Seite
  (curcuma und libxc können verschiedene Seiten nehmen; v-Werte dort im Kerneltest ausgenommen, ε verglichen).

## Validierung (30. Okt 2026, def2-SVP sphärisch, SCF 1e-10 / dE 1e-12 / |g| 1e-9)

Kerne gegen libxc (1760 Punkte: ρ 1e-12..1e2, s 0..30, α 0..1e3 inkl. 0.999/1/1.001 und beide Seiten der vier
SCAN-Abschneidepunkte, 110 Punkte mit τ < τ_W; skalierter Fehler |Δ|/max(|Wert|, Größe der auslöschenden Terme), Definition
in `check_qm_mgga.py`; FD: 4-Punkt, Schritt so, dass α sich um ≤ 1e-5·max(α, 0.01) bewegt, bester von h/10h/100h):

| Kern | ε | v_ρ | v_σ | v_τ | v_ρ / v_σ / v_τ vs FD |
|---|---|---|---|---|---|
| TPSS-X | 8.2e-16 | 1.5e-15 | 5.0e-16 | 1.5e-16 | 2.8e-9 / 2.5e-9 / 3.0e-9 |
| TPSS-C | 1.5e-11 | 1.8e-11 | 1.6e-11 | 1.5e-11 | 2.9e-7 / 5.7e-8 / 3.4e-12 |
| SCAN-X | 6.1e-14 | 7.7e-11 | 5.2e-11 | 5.1e-11 | 2.4e-9 / 2.1e-9 / 2.1e-9 |
| SCAN-C | 3.1e-14 | 2.0e-11 | 2.0e-11 | 2.0e-11 | 1.8e-6 / 1.2e-8 / 1.3e-9 |

FD-Grenzen: ρ ≥ 1e-10, τ ≥ 1e-18, α ≥ 0.01 (darunter Rundungsrauschen der Auslöschung bei großem t bzw. τ-Boden);
SCAN-C v_ρ 1.8e-6 bei ρ = 1e-10, α = 1 ist rauschbegrenzt (libxc-Abweichung dort 2e-11). V_xc gegen FD von E_xc in P
(H2O, Benzol, mit/ohne Screening): ≤ 1.2e-9 relativ, 1 vs 4 Threads bitgleich, V symmetrisch.

Energien gegen PySCF RKS auf identischem Gitter (`ctest -L qm_mgga`, 28 Fälle):

| Fall | TPSS E | ΔE vs PySCF | SCAN E | ΔE vs PySCF |
|---|---|---|---|---|
| H2O sg1 | -76.3511459352 | +3.4e-12 | -76.3220299701 | +7.1e-14 |
| CH4 sg1 | -40.4937398357 | +3.1e-12 | -40.4582874096 | +1.7e-13 |
| NH3 sg1 | -56.5101967983 | +3.3e-12 | -56.4768651814 | -1.4e-14 |
| HF sg1 | -100.3545987083 | +3.5e-12 | -100.3270220873 | -1.4e-14 |
| LiH sg1 | -8.0750496408 | +1.8e-12 | -8.0599160892 | +2.8e-14 |
| Ne sg1 | -128.8031176989 | +3.7e-12 | -128.7724924467 | 0.0 |
| Benzol sg1 | -232.1294666864 | +1.4e-11 | -231.9711879008 | -1.2e-12 |
| H2O sg2 / sg3 | -76.3511540822 / -76.3511541768 | +3.2e-12 / +3.3e-12 | -76.3190848093 / -76.3199663398 | +8.5e-14 / +5.7e-14 |
| NH3 sg2 / sg3 | -56.5102045096 / -56.5102047832 | +3.3e-12 / +3.4e-12 | -56.4762269166 / -56.4755776690 | +2.7e-13 / +2.6e-13 |
| Benzol sg3 | -232.1294454034 | +1.3e-11 | -231.9671744125 | -1.2e-12 |
| H2O RI-J sg1 | -76.3512551408 | +3.4e-12 | -76.3221351844 | +1.1e-13 |
| Benzol RI-J (JK-aux) sg1 | -232.1295572578 | +1.0e-11 | -231.9712785910 | -3.0e-12 |

E_xc ≤ 5.3e-10 Eh (Dichte erster Ordnung), N_Gitter ≤ 1.8e-13. TPSS liegt systematisch +2..4e-12 pro Molekül über PySCF;
Ursache nicht bestimmt (Kandidat, nicht geprüft: libxcs ζ-Schwelle in φ(ζ=1) des spinpolarisierten PBE-Terms).
SCF: 6-16 Iterationen (TPSS 6-14, SCAN 7-16), keine Oszillation beobachtet, keine Schwelle gelockert.

Gegen ORCA 6.1 (`! LibXC(TPSS)|TPSS|LibXC(SCAN) def2-SVP NoRI VeryTightSCF DefGrid1/2/3`, „fein" = DefGrid3 +
`IntAcc 8 AngularGrid 7 GridPruning 0`), Differenz ORCA − curcuma (curcuma 300x1454, Benzol TPSS 250x1454):

| Molekül | LibXC(TPSS) DG1 / DG2 / DG3 / fein | TPSS nativ fein | LibXC(SCAN) DG1 / DG2 / DG3 / fein |
|---|---|---|---|
| H2O | +4.5e-5 / +8.0e-6 / -1.7e-7 / -8.5e-8 | +7.9e-7 | -3.5e-3 / +1.3e-3 / -1.3e-3 / -2.4e-5 |
| CH4 | +2.2e-5 / +2.1e-6 / +2.0e-8 / +4.6e-8 | +7.6e-7 | -4.3e-4 / -9.6e-4 / +5.8e-4 / -7.1e-5 |
| NH3 | +3.1e-5 / +4.8e-8 / +2.3e-7 / +5.4e-9 | +8.0e-7 | -1.6e-3 / -2.1e-4 / -2.7e-4 / -1.2e-6 |
| HF | +3.9e-5 / +2.5e-6 / +7.7e-7 / +7.3e-8 | +1.0e-6 | -4.3e-3 / +2.7e-3 / -1.8e-3 / -9.3e-5 |
| LiH | -1.1e-5 / +3.0e-6 / -1.1e-6 / +1.0e-8 | +2.3e-7 | +6.6e-5 / +1.1e-6 / +3.9e-5 / -6.3e-7 |
| Ne | +6.5e-5 / -7.8e-6 / -1.8e-6 / +3.1e-10 | +1.0e-6 | -4.5e-3 / -2.1e-3 / +4.5e-3 / -9.7e-5 |
| Benzol | +8.0e-5 / +2.9e-5 / -3.8e-6 / -2.4e-6 | +8.2e-7 | -2.6e-3 / -5.4e-3 / +3.0e-3 / -3.2e-4 |

- ORCA LibXC(TPSS) druckt die Konstanten (b 0.40, c 1.59096, e 1.537, κ 0.804, μ 0.21951, β 0.06672455, d 2.8,
  C0 0.53/0.87/0.50/2.26) = curcuma; fein ≤ 8.5e-8 Eh (kleine Moleküle), Benzol 2.4e-6 (ORCAs Winkelgitter, wie WP6).
- ORCAs natives TPSS liegt 2.2e-7..3.3e-6 Eh über ORCAs eigenem LibXC(TPSS); Konstanten werden nicht gedruckt,
  `CBetaPBE` ändert nichts. Die Differenz ist zwischen zwei ORCA-Implementierungen; Ursache nicht bestimmt.
- SCAN: ORCA-Gitter selbst nicht konvergiert (DG1-DG3 springen um mEh); nur Größenordnung vergleichbar.

## Gitterkonvergenz (curcuma; ΔE gegen das feinste Gitter der Zeile, Eh)

| Funktional / Molekül | sg1 | sg2 | sg3 | 99x974 | 150x974 | 250x1454 | Referenz |
|---|---|---|---|---|---|---|---|
| TPSS H2O | +8.2e-6 | +1.0e-7 | +4.9e-9 | +1.3e-8 | +1.4e-8 | 1e-14 | 300x1454 |
| TPSS CH4 | +2.9e-6 | +3.2e-8 | +2.0e-8 | -1.2e-9 | -1.2e-9 | 3e-14 | 300x1454 |
| TPSS Benzol | -2.4e-5 | +9.5e-6 | -2.4e-6 | +1.4e-7 | +1.3e-7 | 0 | 250x1454 |
| SCAN H2O | -2.1e-3 | +8.8e-4 | -2.7e-6 | -3.5e-6 | -4.5e-5 | +8.9e-7 | 300x1454 |
| SCAN CH4 | -5.4e-4 | -6.3e-4 | +1.8e-4 | +1.8e-4 | +6.1e-8 | -9.6e-6 | 300x1454 |
| SCAN Benzol | -2.9e-3 | -3.2e-3 | +1.1e-3 | +1.1e-3 | +4.4e-5 | -3.9e-5 | 300x1454 |

Radial (x590) und Winkel (300x) getrennt, gegen 500x590 bzw. 300x2030:

| | 99 | 150 | 200 | 250 | 300 | 400 | Winkel 302 / 590 / 974 / 1454 |
|---|---|---|---|---|---|---|---|
| SCAN H2O | -3.6e-6 | -4.5e-5 | -2.4e-6 | +6.0e-7 | -2.9e-7 | +1.4e-7 | -1.6e-6 / +1.1e-6 / +2.7e-7 / +4.3e-7 |
| SCAN CH4 | +1.8e-4 | -5.6e-7 | -2.6e-5 | -1.0e-5 | -3.6e-7 | -2.8e-7 | -6.2e-6 / +8.5e-8 / -3.2e-8 / -3.4e-8 |
| SCAN Benzol | +1.1e-3 | +4.1e-5 | -8.9e-5 | -4.2e-5 | -2.9e-6 | -7.4e-7 | (nicht gerechnet) |
| TPSS Benzol | +5.9e-8 | +4.8e-10 | -4e-13 | -6e-13 | 1e-13 | -1e-13 | +9.6e-6 / -2.4e-6 / +1.6e-7 / +3.8e-8 |

- **TPSS:** radial ab 150 Punkten konvergiert; Winkel bestimmt den Fehler: 1e-6 Eh braucht 974 Winkelpunkte für Benzol
  (H2O/CH4 schon sg2/sg3). Empfehlung: `-qm.grid 99x974`.
- **SCAN:** radiale Konvergenz langsam und nicht monoton (bekannte SCAN-Eigenschaft der α-Interpolation); dieselbe
  Schwankung zeigt PySCF auf identischem Gitter (CH4 x590: PySCF 150/250/300 = -40.4577517125/-40.4577611678/
  -40.4577515148, curcuma gleich auf 2e-9). 1e-6 Eh: 300x590 für H2O/CH4, etwa 400x590 für Benzol; sg1-sg3 liegen bis
  3e-3 Eh daneben. Kein Default geändert (sg1); `QMEngine` warnt für TPSS unter 974 Winkelpunkten, für SCAN unter 400x590.

Invarianz (Atomreihenfolge umgekehrt, Translation: ≤ 6.8e-13 Eh auf allen Gittern) und Rotation (Spannweite über die
Ausgangslage und 5 zufällige Drehungen; das Lebedev-Gitter dreht nicht mit):

| | sg1 | sg3 | 150x974 | 250x1454 |
|---|---|---|---|---|
| TPSS H2O / NH3 / Benzol | 1.5e-5 / 1.1e-5 / 5.5e-5 | 4.4e-8 / 3.2e-8 / 4.3e-6 | 1.5e-8 / 2.0e-8 / 4.0e-7 | 6.8e-9 / 4.2e-9 / 7.9e-8 |
| SCAN H2O / NH3 / Benzol | 2.9e-5 / 3.5e-5 / 2.5e-4 | 1.7e-6 / 2.4e-6 / 1.0e-5 | 5.3e-7 / 4.8e-7 / 1.4e-6 | 4.8e-7 / 2.6e-7 / 6.9e-7 |

Weitere Checks (ctest): UpdateMolecule A → B = frisches B (≤ 2.8e-14 Eh); 1 vs 4 Threads Energie bitgleich;
AO-Screening 1e-12 vs 0: ≤ 4.3e-14 Eh (H2O); `-gpu cuda` (Stand Stufe 1, seit Stufe 3 nicht mehr): XC-Teil als gezählter
CPU-Fallback `qm_xc`, Energie = CPU (Δ 0 auf 8 Stellen), `-gpu_strict true` endet mit Exit 3; Gradientanfrage und `-opt`
enden mit Fehler (`hasGradient()` false). LDA/PBE/B3LYP/HF: `-dump_gradient` gegen das Binary vor der Änderung (16b2b122)
60/60 CPU-Läufe bitgleich (H2O, NH3, Benzol; 1/4 Threads; Screening an/aus; RI-JK); GPU-Läufe innerhalb des
Lauf-zu-Lauf-Rauschens der Referenz (Energie bitgleich, Gradient ≤ 1e-13 wie Referenz gegen Referenz).

## Stufe 2: analytischer Gradient (CPU, 30. Okt 2026)

- **Code:** `buildXCGradientMeta` in `qm_xc.cpp` (eigene Funktion, LDA/GGA-Pfade unverändert), Aufruf über
  `qmxc::buildXCGradient`; `QMEngine::hasGradient()` true für alle KS-Funktionale; Testtreiber `test_qm_dft_grad`
  (`--functional tpss|scan`), Tests `test_cases/qm_mgga_grad/`, Referenz `scripts/qm_dft_grad_pyscf_reference.py --set mgga`.
- **Formel** (Basisteil, Johnson/Gill/Pople 1993 um τ erweitert): ∂τ/∂A_k = −Σ_{m auf A} Σ_j ∂_k∂_jφ_m (P ∂_jφ)_m,
  also im WP8-Ausdruck zusätzlich Σ_j (w v_τ/2) ∂_k∂_jφ_m (P ∂_jφ)_m. Es reichen AO-Zweitableitungen (deriv 2), deriv 3
  ist nicht nötig (τ enthält erste Ableitungen, seine Kernableitung eine mehr).
- **Gitterantwort:** Punktteil über dieselbe Translationsidentität wie WP8 (der τ-Term steckt in der Punktsumme);
  die Becke-Gewichtsableitung braucht nur e_p = ρ ε_xc(ρ, σ, τ) (`edens`) und ist unverändert wiederverwendet.
  Geprüft durch FD mit mitwanderndem Gitter und PySCF `grid_response=True` (Tabellen unten).
- Ableitungen am geklemmten Punkt (τ-Boden, σ ≤ 8ρτ) wie V_xc: wo eine Klemme greift, ist der Gradient nicht die
  Ableitung der geklemmten Energie (für P positiv semidefinit nur auf Rundungsniveau).
- `-gpu cuda` (Stand Stufe 2, seit Stufe 3 auf dem Device): XC-Gradient als gezählter CPU-Fallback `qm_xc`, 2e-Teil auf
  dem Device; H2O TPSS/SCAN Gradient gegen `-gpu none` ≤ 3.6e-14 Eh/Å, `-gpu_strict true` endete mit Exit 3.

**(1) gegen PySCF 2.14** (`nuc_grad_method()`, identisches Gitter, def2-SVP, max |g − g_PySCF| Eh/Bohr, ohne / mit
Gitterantwort; ctest `qm_mgga_grad_pyscf_*`, 17 Fälle, Gate 1e-9):

| Fall | TPSS | SCAN | Größe der Gitterantwort TPSS / SCAN |
|---|---|---|---|
| H2O / NH3 / CH4 / HF / Benzol sg1 | ≤ 9.7e-12 / ≤ 9.7e-12 | ≤ 7.3e-12 / ≤ 7.6e-12 | 1.2e-5..9.5e-5 / 2.6e-4..9.9e-2 |
| H2O sg3, H2O RI-J sg1, NH3 150x974 | ≤ 8.6e-12 / ≤ 8.5e-12 | ≤ 6.6e-12 / ≤ 6.7e-12 | 6.3e-8..1.1e-5 / 5.0e-3..3.9e-2 |
| Benzol TPSS RI-JK sg1 | 9.7e-12 / 9.7e-12 | – | 5.5e-5 |
| H2O / NH3 250x1454 (ad hoc, nicht im ctest) | 5.4e-12 / – | 1.9e-12 / 6.2e-12 | |

**(2) FD der eigenen Energie** (4-Punkt, Schritt h in Bohr, in Å verschoben, alle Atome, max |g − FD| Eh/Bohr):

| Fall | mitwanderndes Gitter vs mit Antwort | eingefrorenes Gitter vs ohne Antwort |
|---|---|---|
| H2O TPSS sg1 | 2.7e-11 (h 2e-3) | 3.1e-7 (h 2e-3), 1.0e-6 (5e-4), 2.7e-9 (1e-4), 2.1e-9 (2e-5) |
| H2O SCAN sg1 | 7.2e-10 (2e-3), 2.0e-10 (5e-4), 9.6e-10 (1e-4) | 1.4e-5 (2e-3), 5.2e-8 (5e-4), 1.3e-9 (1e-4), 5.3e-9 (2e-5) |
| NH3 TPSS sg1 | 8.0e-11 (2e-3), 3.0e-10 (5e-4) | 1.0e-7 (2e-3), 3.0e-9 (1e-4), 2.0e-9 (5e-5) |
| H2O TPSS / SCAN 250x1454 | 1.9e-10 / 1.3e-10 (5e-4) | – / 1.1e-9 (1e-4) |

- Eingefrorenes Gitter, großes h: abweichend ist jeweils nur die Komponente des Schweratoms entlang der Symmetrieachse
  (O_z, N_z); sie fällt mit kleinerem h auf das FD-Rauschen (≤ 5e-9), nicht monoton (TPSS 5e-4 schlechter als 2e-3).
  Deutung (nicht bewiesen): Punkte des eigenen Atoms auf der Achse nahe am Kern, über die der Kern bei raumfestem
  Gitter wandert; der Meta-GGA-Integrand ist dort bei Schritten ≥ Abstand der innersten Radialschalen nicht glatt
  genug für den 4-Punkt-Stern. LDA/GGA zeigten das nicht (WP8). Die PySCF-Werte ohne Antwort (= Ableitung bei
  raumfestem Gitter) stimmen unabhängig davon auf ≤ 6e-12. ctest: mitwandernd h 5e-4 (Gate 2e-9), eingefroren
  h 1e-4 (Gate 1e-8).

**(3) Gitterstreuung des Gradienten** (H2O/NH3, max |g − g(500x1454)| Eh/Bohr, mit / ohne Gitterantwort; die
Referenz teilt das 1454er Winkelgitter, Winkelfehler daher getrennt):

| | sg1 | sg3 | 99x974 | 150x974 | 250x1454 | 300x1454 | 400x590 |
|---|---|---|---|---|---|---|---|
| TPSS H2O mit | 7.4e-6 | 2.4e-7 | 2.6e-7 | 2.8e-7 | 4e-13 | 5e-13 | 2.5e-7 |
| TPSS NH3 mit | 3.8e-5 | 4.4e-7 | 5.8e-8 | 5.9e-8 | 5e-13 | 5e-13 | 4.6e-7 |
| SCAN H2O mit / ohne | 5.8e-4 / 4.0e-2 | 1.1e-4 / 5.0e-3 | 8.3e-5 / 5.0e-3 | 3.9e-5 / 1.4e-2 | 9.3e-6 / 2.2e-3 | 6.2e-7 / 2.1e-4 | 2.3e-5 / 5.0e-5 |
| SCAN NH3 mit / ohne | 6.6e-4 / 2.3e-2 | 1.1e-4 / 1.0e-2 | 1.0e-4 / 1.0e-2 | 3.5e-5 / 5.2e-3 | 4.5e-6 / 3.2e-4 | 1.0e-6 / 8.0e-5 | 1.8e-5 / 1.2e-4 |

Winkel 300x1454 gegen 300x2030 (mit Antwort): TPSS H2O/NH3 1.5e-8/1.9e-8, SCAN 1.0e-5/3.7e-6. Rotation (Gitter dreht
nicht mit, max Δg zurückgedreht): TPSS sg1 4.0e-5/7.7e-5, 150x974 4.2e-7/1.7e-7, 250x1454 5.0e-8/4.6e-8; SCAN sg1
1.8e-4/1.2e-4, 150x974 1.4e-5/1.2e-5, 250x1454 1.2e-5/6.4e-6 (H2O/NH3). Translation: Gradient unverändert ≤ 3.3e-13,
Σ_A g_A mit Antwort ≤ 1e-12.

- **TPSS:** Gradient auf 250x1454 bis ~5e-8 Eh/Bohr gitterkonvergiert, sg1 bis 8e-5.
- **SCAN:** ohne Gitterantwort unbrauchbar (Fehler bis 4e-2 Eh/Bohr auf sg1, 2e-3 noch bei 250x1454); mit Antwort
  bleibt eine Streuung von ~1e-5 Eh/Bohr selbst bei 250x1454/300x1454 (Rotation, Winkelgitter). `-qm.grid_response
  false` mit SCAN ist daher nicht zu empfehlen; für Optimierungen feine Gitter (≥ 300x1454) und lockere
  Gradientenkriterien erwarten.

**(4) `-opt`** (`-opt.convergence_preset tight`, sg1, Referenz scipy-BFGS auf PySCF-Energie + Gradient mit Antwort):

| Start | ΔE (curcuma-Min − PySCF-Min) | max Δd_ij (Å) | max\|g\| am Ende | curcuma am PySCF-Minimum ΔE / max\|g\| |
|---|---|---|---|---|
| H2O_distorted TPSS | 3.5e-11 | 5.5e-6 | 2.6e-6 | 3.3e-12 / 1.6e-8 |
| NH3 (`qm_1e/NH3.xyz`) TPSS | 7.3e-11 | 8.7e-6 | 3.5e-6 | 3.0e-12 / 5.3e-8 |
| formaldehyde_distorted TPSS | 1.2e-11 | 2.5e-6 | 1.3e-6 | 5.3e-12 / 2.0e-7 (= PySCF-Restgradient) |

CLI ohne tight (Default-Kriterien, 4 Threads): H2O_distorted TPSS 3.9 s, E 9.3e-8 über dem PySCF-Minimum, max |g| am
Ende 4.3e-6; NH3 TPSS 5.7 s, ΔE 1.1e-9, 6.0e-6; SCAN H2O/NH3 laufen durch (5.8/4.1 s, Endgradient 1.9e-6/9.8e-7),
für SCAN ohne PySCF-Geometrievergleich.

**(5) Threads, UpdateMolecule, Regression:** 1 vs 4 Threads Gradient bitgleich (NH3 TPSS/SCAN, Benzol TPSS);
UpdateMolecule gegen neue Instanz max Δg 1.6e-12 (TPSS) / 1.2e-12 (SCAN). LDA/PBE/B3LYP/HF: `-dump_gradient` gegen das
Binary von 3ff3ae37 60/60 CPU-Läufe bitgleich (H2O_distorted, NH3, Benzol; 1/4 Threads; beide VWN-Varianten,
`grid_response false`, RI-JK, Screening an/aus); nach dem vollen Neubau 5/5 Stichproben (H2O) erneut bitgleich.
`ctest -L qm_mgga_grad` 38/38.

**(6) gegen ORCA 6.1** (`! LibXC(TPSS)|TPSS def2-SVP NoRI VeryTightSCF DefGrid3 EnGrad`, `IntAcc 8 AngularGrid 7
GridPruning 0`; curcuma 250x1454 mit Antwort; max |g_ORCA − g| Eh/Bohr):

| | LibXC(TPSS) | TPSS nativ | ΔE LibXC(TPSS) |
|---|---|---|---|
| H2O | 2.1e-6 | 2.1e-6 | -8.6e-8 |
| NH3 | 5.6e-7 | 5.8e-7 | +5.3e-9 |

ORCA IntAcc 10 statt 8 ändert ORCAs H2O-Gradient um 1.7e-7; curcuma 250x770 gegen 250x1454 2.7e-7. Die 2e-6 bei H2O
(O_z) sind damit nicht durch die Winkelpunktzahl erklärt; Größenordnung wie WP8 (PBE H2O 1.3e-6 gegen ORCA). Ursache
nicht bestimmt (ORCAs Gitter, Partition und Gewichtsableitung nicht nachgebaut); SCAN nicht gegen ORCA verglichen
(ORCAs SCAN-Energie selbst schon gitterunkonvergiert, Stufe 1).

## Stufe 3: Meta-GGA auf CUDA (30. Okt 2026)

- **Code:** `cuda/qm_gpu_xc.cu` (`densityXCMetaKernel`, `metaDensity`, `metaVxcProduct`, `combineTauKernel`,
  `gradTerms<true>`), Flag `QMGpuXCSetup::mgga` (`qm_gpu_xc_backend.h`, angehängt; Plugin-Eintrag umbenannt zu
  `curcuma_cuda_create_qm_xc_mgga`), `QMEngine::ensureGpuXC` (`qm_scf.cpp`: der Meta-GGA-Fallback ist entfernt).
  Punktfunktional = `qmxc::kern::evaluatePointMGGA` aus `qm_xc_kernels.h` (dieselbe Quelle wie die CPU, keine Kopie).
- **Ablauf pro Block** (AO-Screening-Gruppen wie LDA/GGA): ein strided-batched DGEMM bildet [φP, ∂_xφP, ∂_yφP, ∂_zφP],
  τ_q = ½ Σ_k Σ_m ∂_kφ_qm (∂_kφP)_qm; V_xc: W += φᵀY + Σ_k (∂_kφ)ᵀ diag(w v_τ/4) ∂_kφ (vier DGEMMs in fester Folge),
  V = W + Wᵀ (= Term von `buildVxcMeta`); Gradient: T_k += (w v_τ/2) Σ_j ∂_k∂_jφ ∘ (∂_jφP) (AO-Ableitung 2), Punktteil
  und Becke-Gewichtsableitung (e_p = ρ ε_xc) unverändert. Keine FP-Atomics im XC-Teil, Blockgröße und Wellen nicht
  vom freien Speicher abhängig (wie LDA/GGA). Der resident KS-SCF-Loop nimmt V_xc über `vxcDevice` unverändert.
- Verworfen: die vier W-Produkte als ein strided-batched DGEMM plus Summe war auf der RTX 5090 3x langsamer
  (Benzol 150x974: 9.2 s statt 3.1 s für 13 Aufbauten). Ein Laufzeit-`if` im Gradientenkern änderte die LDA/GGA-Bits
  (1.7e-15 Eh/Å, andere FMA-Kontraktion); daher Template-Parameter, `<false>` ist der alte Code.
- Fallback-Zähler: tpss/scan erzeugen keinen `qm_xc`-Eintrag mehr; `-gpu_strict true` läuft durch (ctest
  `qm_gpu_dft_scf_*_tpss*`/`*_scan*` mit `-gpu_strict true`). Gezählter Fallback bleiben die LDA/GGA-Gründe (kein
  Plugin, Schalen über f, > 32 Primitive, zu wenig Speicher: ctest `qm_gpu_dft_fallback_h2o_scan_memory`, Gerätefehler).

GPU gegen CPU (RTX 5090, def2-SVP; „gleiche Dichte" = `CURCUMA_QM_GPU_XC_CHECK=1`, jeder Aufbau auch auf der CPU;
SCF/Gradient = zwei Läufe, `-qm.scf_threshold 1e-10`, 4 CPU-Threads, `check_qm_gpu_dft.py`):

| Molekül (n, Punkte sg1) | Funktional | Gitter | max rel\|dV\| | rel\|dE_xc\| | \|dE\| SCF | Iter. CPU/GPU | XC-Gradient gleiche Dichte (Basis/Punkte/Gewichte) | \|dg\| Lauf |
|---|---|---|---:|---:|---:|---:|---:|---:|
| H2O (24, 29100) | tpss / scan | sg1 | 5.5e-16 / 4.0e-16 | 7.9e-16 / 9.5e-16 | 0 / 0 | 13/13, 13/13 | ≤ 1.5e-14 / ≤ 7.7e-14 | 7.0e-14 / 6.5e-14 |
| H2O | tpss / scan | 150x974 | 5.5e-16 / 9.4e-16 | 2.3e-15 / 4.0e-15 | 0 / 0 | 13/13, 13/13 | ≤ 3.8e-13 / ≤ 4.1e-13 | 3.8e-13 / 3.9e-13 |
| NH3 (29, 38800) | tpss / scan | sg1 | 7.9e-16 / 9.4e-16 | 6.6e-16 / 8.9e-16 | 0 / 0 | 13/13, 13/13 | ≤ 9.8e-15 / ≤ 6.2e-14 | 1.2e-14 / 4.7e-14 |
| NH3 | tpss / scan | 150x974 | 1.4e-15 / 9.3e-16 | 2.4e-15 / 2.2e-15 | 0 / 0 | 13/13, 13/13 | ≤ 2.1e-13 / ≤ 3.1e-13 | 2.4e-13 / 3.0e-13 |
| Benzol (114, 116400) | tpss / scan | sg1 | 1.1e-15 / 9.1e-16 | 8.0e-16 / 8.4e-16 | 0 / 0 | 14/14, 15/15 | ≤ 9.0e-14 / ≤ 8.4e-14 | 1.5e-11 / 1.0e-11 |
| Benzol | tpss / scan | 150x974 | 1.5e-15 / 1.7e-15 | 6.7e-15 / 7.4e-15 | 0 / 0 | 14/13, 14/14 | ≤ 3.0e-13 / ≤ 1.4e-13 | 1.3e-11 / 9.3e-12 |
| Koffein (246, 232800) | tpss | sg1 | 6.9e-16 | 9.3e-16 | 1e-10 | 30/32 | ≤ 1.8e-13 | 5.2e-11 |

\|dE\| auf den 10 gedruckten Stellen, \|dg\| in Eh/Bohr; der 2e-Teil (Q-G7, Atomics) trägt die 1e-11 der Läufe
(gleiche Dichte ≤ 1.3e-11 Benzol, 5.3e-11 Koffein). Koffein TPSS 1e-10: 32 statt 30 Iterationen bei gleicher Energie
(langsam konvergierende SCF, Rundung verschiebt das Ende); bei 1e-8 (Zeittabelle) 26/26. Koffein auf 150x974 nur
in der Zeittabelle (Energie gleich), nicht gegen gleiche Dichte geprüft (CPU-Referenz ~4 s pro Aufbau x 30).

Zeiten (CPU 8 Threads `-gpu none`, GPU `-gpu cuda -qm.gpu_scf on -gpu_strict true`, `-qm.scf_threshold 1e-8`, Wand in
s; V_xc pro Aufbau = letzter Aufbau, verbosity 3; Gradient = XC-Basis- und Punktteil, CPU plus Gewichtsableitung):

| System | Funktional | Gitter | \|dE\| | Iter. CPU/GPU | SCF-Loop CPU / GPU | V_xc pro Aufbau CPU / GPU | Wand CPU / GPU |
|---|---|---|---:|---:|---:|---:|---:|
| Benzol | tpss | sg1 | 0 | 12/12 | 1207 / 464 ms | 66.7 / 33.6 ms | 1.59 / 1.07 |
| Benzol | scan | sg1 | 0 | 13/13 | 1292 / 506 ms | 65.9 / 33.9 ms | 1.68 / 1.01 |
| Benzol | tpss | 150x974 | 0 | 12/12 | 11.2 / 5.7 s | 896 / 468 ms | 12.6 / 6.9 |
| Benzol | scan | 150x974 | 0 | 12/12 | 11.2 / 5.7 s | 897 / 472 ms | 12.6 / 7.0 |
| Koffein | tpss | sg1 | 1e-10 | 26/26 | 17.3 / 5.6 s | 310 / 180 ms | 18.3 / 7.1 |
| Koffein | scan | sg1 | 0 | 24/24 | 16.9 / 5.1 s | 308 / 181 ms | 17.8 / 6.7 |
| Koffein | tpss | 150x974 | 0 | 26/26 | 111.8 / 56.0 s | 3982 / 2120 ms | 117.0 / 60.1 |

| Gradient (sg1) | XC-Gradient CPU (+ Gewichtsabl.) / GPU | Wand CPU / GPU | max\|g(GPU) − g(CPU)\| Eh/Bohr |
|---|---:|---:|---:|
| Benzol tpss / scan | 75.1 (+11.3) / 33.0 ms, 74.9 (+10.7) / 33.3 ms | 4.7 / 1.2 s, 4.9 / 1.3 s | 1.3e-11 / 1.0e-11 |
| Koffein tpss / scan | 343 (+83) / 154 ms, 342 (+85) / 154 ms | 87.7 / 9.2 s, 87.1 / 8.8 s | 3.5e-10 / 2.8e-10 |

- Der XC-Teil ist auf dieser Karte (FP64 schwach) nur 1.7-2x schneller als 8 CPU-Threads; Stufen pro Aufbau Benzol
  150x974 (`CURCUMA_QM_GPU_XC_PROFILE=1`): AO 37, X = [φ, ∂φ]P 152, Dichte+Funktional 64, W 237 ms (PBE gleiches Gitter
  37 / 49 / 26 / 65 ms): die viermal so vielen FP64-GEMMs dominieren. Der Gewinn der Läufe kommt vom 2e-Teil (J/K,
  2e-Gradient) auf dem Device. Koffein-Gradient 1e-8: Differenz der Läufe aus der SCF-Konvergenz.
- `-opt` Benzol tpss (sg1, Default-Kriterien, 8 Threads): CPU 35.1 s, GPU 5.8 s; gleiche Endenergie (E_elec
  -434.2210594551, E_nn 202.0898086630 Eh auf 10 Stellen), Koordinaten auf den 6 geschriebenen Stellen gleich (nur das
  Vorzeichen von Nullen verschieden). ctest `qm_gpu_dft_opt_h2o_tpss` (tight): ΔE_total 1.0e-10 Eh, max\|Δr_ij\| 6.7e-7 Å.
- Regression: CPU `-dump_gradient` gegen das Binary von 7842f654 (git worktree, gleiche CMake-Optionen inkl. MKL) 51/51
  bitgleich (H2O_distorted, NH3, Benzol x hf, lda beide VWN, pbe, pbe ohne Antwort, pbe ohne Screening, b3lyp beide VWN,
  b3lyp RI-JK, tpss, scan, tpss ohne Screening, scan 150x974; 1/4 Threads). GPU (`-qm.gpu_grad off -qm.gpu_grad1e off`)
  lda/pbe/b3lyp 15/18 bitgleich, die 3 RI-JK-Fälle Energie gleich, \|dg\| ≤ 4.0e-14 Eh/Å; das Referenz-Binary selbst gibt
  bei RI-JK in 3 Läufen 3 verschiedene md5 (Lauf-zu-Lauf-Rauschen des RI-Pfads).
- ctest `-L qm_gpu_dft` um 22 Fälle (`*_tpss*`, `*_scan*`, nur mit `USE_CUDA` registriert): V_xc gleiche Dichte (H2O,
  NH3, Benzol, H2O SCAN 150x974, Benzol ohne Screening), SCF im Device-Loop (Benzol tpss/scan, NH3 SCAN 150x974; RI-JK
  Host-Loop), Gradienten (H2O, Formaldehyd, Benzol, ohne Antwort, 64er-Blöcke, 1024er-Gruppen, ohne Screening), `-opt`
  H2O, Speicher-Fallback. Iterationstoleranz die bestehende (±1): auf H2O/NH3/Benzol sind die Zahlen gleich bis auf
  Benzol TPSS 150x974 (CPU 14, GPU 13). Der opt-Modus von `check_qm_gpu_dft.py` vergleicht seitdem E_elec + E_nn statt E_elec
  (TPSS-Minima 8e-7 Å auseinander gaben 1.2e-6 Eh in E_elec, 1e-10 Eh in der Gesamtenergie). `qm_mgga_gpu_H2O_*`
  erwarten jetzt den Device-Pfad (kein `qm_xc`-Fallback, `-gpu_strict true` Exit 0, Energie = CPU ≤ 1e-8 Eh) und
  melden ohne CUDA-Plugin Exit 77 (skip; geprüft mit einer Binary-Kopie ohne Plugin).
- `ctest -L "qm_mgga|qm_mgga_grad|qm_lda|qm_pbe|qm_b3lyp|qm_dft_grad|qm_gpu_dft|qm_grid|qm_hf3c"` (-j4): 320/322 im
  ersten Lauf, die zwei Fehlschläge waren `qm_mgga_gpu_H2O_tpss/scan` mit der alten Fallback-Erwartung; nach deren
  Anpassung 2/2 und der volle Satz erneut 322/322 (275 s).

## Nicht getestet / nicht implementiert

- Nicht implementiert: offene Schale, Z > 10, r2SCAN/revTPSS/TPSSh/Hybride, DFT-D, Hessian; auf dem Device: τ aus den
  besetzten Orbitalen (∂_kφ C_occ statt ∂_kφ P, für Koffein ~5x weniger FLOPs im X-Schritt; das Device-Objekt bekommt
  nur P), FP32-Phasen, ROCm/Vulkan.
- Stufe 3 nicht getestet: andere GPUs als die RTX 5090, Moleküle über Koffein, Stratmann-Partition auf dem Device,
  `-qm.cartesian_d`, Multi-GPU, `-md`, SCAN-Optimierung auf dem Device (nur Gradienten), Koffein 150x974 gegen gleiche
  Dichte.
- Nicht getestet: Moleküle größer als Benzol, geladene Systeme, `-qm.cartesian_d`, Stratmann-Partition und
  `radii_adjust becke|none` mit Meta-GGA, integral-direct J mit Meta-GGA (Pfad wie GGA, nicht gemessen), `-md`;
  Gradient: SCAN-Optimierung gegen eine Referenzgeometrie, ORCA-Gradientenvergleich (nicht gerechnet), Wirkung der
  Klemmen (τ-Boden, σ ≤ 8ρτ) auf den Gradienten außer über die Testmoleküle.
- SCAN-Gleichungsnummern nicht am Paper geprüft; Formeln folgen libxc 7.0.0 (über libxc-Werte validiert, nicht über das
  Paper).
- Menschlicher Produktionstest ausstehend.
