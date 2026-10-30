# WP11 — Meta-GGA TPSS und SCAN (+ τ, rung 3)

- **Status Stufe 1 (Energien, CPU):** 🤖 AI-generated, ⚙️ machine-tested (`ctest -L qm_mgga` 50/50, 30. Okt 2026);
  menschlicher Produktionstest ausstehend. **Stufe 2 (Gradient) und Stufe 3 (GPU) offen** (TODO.md).
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
AO-Screening 1e-12 vs 0: ≤ 4.3e-14 Eh (H2O); `-gpu cuda`: XC-Teil als gezählter CPU-Fallback `qm_xc` („meta-GGA ... has no
device XC kernel yet"), Energie = CPU (Δ 0 auf 8 Stellen), `-gpu_strict true` endet mit Exit 3; Gradientanfrage und `-opt`
enden mit Fehler (`hasGradient()` false). LDA/PBE/B3LYP/HF: `-dump_gradient` gegen das Binary vor der Änderung (16b2b122)
60/60 CPU-Läufe bitgleich (H2O, NH3, Benzol; 1/4 Threads; Screening an/aus; RI-JK); GPU-Läufe innerhalb des
Lauf-zu-Lauf-Rauschens der Referenz (Energie bitgleich, Gradient ≤ 1e-13 wie Referenz gegen Referenz).

## Nicht getestet / nicht implementiert

- Nicht implementiert: Gradient (Stufe 2), GPU-Kern (Stufe 3), offene Schale, Z > 10, r2SCAN/revTPSS/TPSSh/Hybride, DFT-D.
- Nicht getestet: Moleküle größer als Benzol, geladene Systeme, `-qm.cartesian_d`, Stratmann-Partition und
  `radii_adjust becke|none` mit Meta-GGA, integral-direct J mit Meta-GGA (Pfad wie GGA, nicht gemessen), `-md`.
- SCAN-Gleichungsnummern nicht am Paper geprüft; Formeln folgen libxc 7.0.0 (über libxc-Werte validiert, nicht über das
  Paper).
- Menschlicher Produktionstest ausstehend.
