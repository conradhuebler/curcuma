# WP8 — Analytischer Gradient (HF Sep 2026, Kohn-Sham Okt 2026)

- **Status:** 🤖 AI-generated, ⚙️ machine-tested: `hf`/`hf-3c` (Sep 2026, `ctest -L qm_grad`), `lda`/`pbe`/`b3lyp`
  (Okt 2026, `ctest -L qm_dft_grad` 55/55); menschlicher Produktionstest ausstehend
- **Code:** `QMEngine::computeGradient()` (`qm_engine.cpp`), `QMEngine::xcGradient()` (`qm_scf.cpp`),
  `qmxc::buildXCGradient` (`qm_xc.cpp`), `qmgrid::gridWeightGradient` und AO-Zweitableitungen `AOEvaluator` deriv 2
  (`qm_grid.cpp`), `qmint::gradientTwoElectronScaled` (`qm_integrals.cpp`, deklariert in `qm_ri.h`),
  `RIJK::gradientDensities(P, kx, ...)` (`qm_ri.cpp`)
- **Test:** `test_cases/qm_dft_grad/` (Treiber `test_qm_dft_grad.cpp`, Prüfskript ohne numpy), Referenz
  `scripts/qm_dft_grad_pyscf_reference.py` (PySCF 2.14, scipy 1.18.1, venv `~/.venvs/qmref`)

## Was gebaut ist (Kohn-Sham, geschlossene Schale)

dE/dR_A = Σ P dH/dR − Σ W dS/dR + ½ Σ D (mn|ls)^A + dE_nn/dR + dE_xc/dR, mit

- W = 2 Σ_i ε_i C C^T aus den **KS**-Orbitalenergien (= ½ P F_KS P inkl. V_xc und −(a_x/2) K): der Pulay-Term
  kommt unverändert aus `gradientOneElectron` (Q-G8-Device-Pfad gilt auch für KS).
- D = P_mn P_ls − (a_x/4)(P_ml P_ns + P_ms P_nl): `gradientTwoElectronScaled(kx)`, kx = 0 (LDA/PBE: nur J),
  0.2 (B3LYP); kx = 1 ist bitweise der HF-Ausdruck. RI: `gradientDensities` mit G = c P − 2 kx Γ, M2 = ½ c cᵀ − kx M
  (kx = 0: Γ entfällt); RI-J addiert kx · exakten Austauschgradienten.
- XC-Basisteil (Johnson/Gill/Pople 1993): dE_xc/dA_k = −2 Σ_g w_g Σ_{m auf A} [∂_k φ_m (Z P)_m
  + Σ_j v_j ∂_k∂_j φ_m (φ P)_m], Z = v_ρ φ + Σ_j v_j ∂_j φ, v_j = 2 v_σ ∂_j ρ (LDA: v_j = 0). Zweite
  AO-Ableitungen neu in `AOEvaluator::evaluate(deriv = 2)`.
- **Gitterantwort** (`-qm.grid_response true`, Default): (i) Punkte wandern mit ihrem Atom B: Σ_p w_p ∂_r e(r_p)
  auf B, berechnet als −Σ aller Basisanteile des Punktes (Translationsidentität); (ii) Ableitung der Becke-Gewichte
  w_p = W_p P_B/Σ_C P_C für A ≠ B bei festem Punkt (Becke 1988, Größenanpassung ν = μ + a(1 − μ²), Stratmann-Polynom
  ebenso), d w/d R_B = −Σ_{A≠B} d w/d R_A (Johnson/Gill/Pople 1993; Gill/Johnson/Pople/Frisch, CPL 197, 499 (1992)).
  „P_C ohne einen Faktor“ über Präfix/Suffix-Produkte, keine Division durch eine eventuell verschwindende Zellfunktion.
- Ohne Gitterantwort (`false`) = PySCF-Default `grid_response=False`: nur der Basisteil.
- Beziehung: die gemeldete Energie wird mit dem **mitwandernden** Gitter berechnet. Gradient MIT Gitterantwort =
  exakte Ableitung dieser Energie (= FD mit mitwanderndem Gitter). Gradient OHNE Gitterantwort = exakte Ableitung der
  Energie bei raumfest eingefrorenem Gitter (= FD mit `QMEngine::freezeXCGrid`, nur Validierungs-Hook). Der Unterschied
  ist der Gitterfehler des Gradienten (sg1 ~1e-5 Eh/Bohr, sg3 ~3e-7).
- Determinismus: XC-Teile mit fester Chunk-Zahl; die Reduktion von `gradientOneElectron`, des 2e-Treibers (auch HF, MP2)
  und von `RIJK::threeIndexGradient`/`metricGradient` geändert von „pro Thread + omp critical“ auf feste
  Round-Robin-Chunks: 1 vs 4 Threads bitweise gleich (exakter Pfad, RI-J). HF-Gradienten ändern sich dadurch nur in der
  Summationsreihenfolge (~1e-14).
- Kosten (Benzol B3LYP/def2-SVP sg1, 8 Threads): 2e 2.9 s, 1e 62 ms, XC 65 ms, Gewichtsableitung +11 ms.
- Zitate (`qmcite::ksGradient`): Johnson/Gill/Pople 1993 (10.1063/1.464906), mit Gitterantwort zusätzlich
  Gill/Johnson/Pople/Frisch 1992 (10.1016/0009-2614(92)85807-M) und Becke 1988; DOIs gegen Crossref geprüft
  (Okt 30 2026). Der Auftrag nannte „Gill/Johnson/Pople/Frisch 1993“; Crossref findet die Arbeit mit diesen vier Autoren
  in CPL 1992.

## Default-Entscheidung `-qm.grid_response true`

| Kriterium (sg1) | mit Gitterantwort | ohne |
|---|---|---|
| Analytisch vs FD der gemeldeten Energie (H2O) | ≤ 3.6e-11 | 1.1e-5..1.5e-5 |
| Σ_A g_A (H2O / NH3) | 1e-14 | 1.7e-5 / 1.5e-5 |
| `-opt` PBE/B3LYP H2O, H2CO (sg1, Endgradient mit Antwort gemessen) | 1.7e-6..4.2e-6, ΔE ≤ 1.1e-7 | 2.5e-5..2.4e-4, ΔE 8e-7..5.6e-6 |
| Laufzeit dieser vier `-opt` | 6.8 / 10.8 / 36.3 / 36.7 s | 26.2 / 6.3 / 56.4 / 56.9 s |
| Zusatzzeit (Benzol) | +11 ms von ~3 s | 0 |

Mit Antwort ist der Gradient die Ableitung der Energie, die der Optimierer sieht; ohne bleibt ein Gitterfehler, der
`-opt` neben dem Minimum stehen lässt. PySCF-Default-kompatibel per `-qm.grid_response false`.

## Validierung (Okt 2026, def2-SVP sphärisch, SCF 1e-10 / dE 1e-12 / |g| 1e-9, Eh/Bohr)

**(1) FD der eigenen Energie** (4-Punkt, h = 2e-3 Bohr in Å verschoben, alle Atome; max |g − FD|):

| Fall | FD mitwanderndes Gitter vs mit Antwort | FD eingefrorenes Gitter vs ohne Antwort | FD mitwandernd vs ohne Antwort |
|---|---|---|---|
| H2O LDA (VWN5) | 3.2e-11 | 2.2e-11 | 1.5e-5 |
| H2O PBE | 2.2e-11 | 3.4e-10 | 1.2e-5 |
| H2O B3LYP | 3.6e-11 | 1.2e-10 | 1.1e-5 |
| H2O B3LYP RI-J / RI-JK | 3.6e-11 / 3.7e-11 | ctest | – |
| H2O PBE RI-J / RI-JK, LDA RI-J | 5.2e-11 / 3.7e-11 / 3.7e-11 | ctest | – |
| NH3 LDA / B3LYP VWN-RPA | 6.9e-11 / 5.7e-11 | ctest | – |

**(2) gegen PySCF 2.14** (`nuc_grad_method()`, `grid_response` False/True, identisches Gitter wie WP5-WP7, 27 Fälle):
max |g − g_PySCF| ohne Antwort ≤ 3.6e-11, mit Antwort ≤ 2.6e-11 (HF B3LYP); LDA/PBE ≤ 8.9e-12; inkl. RI-J (J gefittet,
K exakt), RI-JK, Benzol, sg3, beide VWN-Varianten. Größe der Gitterantwort selbst: 1e-5..6.4e-5 (sg1), 3.2e-7 (sg3).

| Molekül | LDA | PBE | B3LYP | B3LYP VWN-RPA | (ohne / mit Antwort) |
|---|---|---|---|---|---|
| H2O | 2.0e-13 / 1.9e-13 | 7.4e-12 / 7.4e-12 | 2.1e-11 / 1.5e-11 | 2.1e-11 / 1.5e-11 | |
| NH3 | 2.4e-12 / 2.4e-12 | 5.0e-12 / 5.0e-12 | 2.0e-11 / 7.8e-12 | 2.0e-11 / 8.0e-12 | |
| CH4 | 4.2e-12 / 4.2e-12 | 1.7e-12 / 1.7e-12 | 1.4e-11 / 5.3e-13 | 1.8e-11 / 4.3e-12 | |
| HF | 4.0e-12 / 4.0e-12 | 7.4e-12 / 7.4e-12 | 3.6e-11 / 2.6e-11 | 3.5e-11 / 2.5e-11 | |
| Benzol | 8.9e-12 / 8.9e-12 | 7.9e-12 / 7.9e-12 | 1.3e-11 / 9.5e-12 | 1.5e-11 / 1.1e-11 | |

RI: H2O B3LYP RI-J 2.1e-11/1.5e-11, RI-JK 2.2e-11/1.5e-11; NH3 PBE RI-JK 8.5e-12, LDA RI-J 4.1e-12; Benzol B3LYP-RPA
RI-JK 1.2e-11/9.6e-12; H2O B3LYP sg3 1.8e-11/1.4e-11.

**(3) gegen ORCA 6.1 grob** (`! LSD|PBE|B3LYP def2-SVP NoRI VeryTightSCF EnGrad`, `IntAcc 8 AngularGrid 7
GridPruning 0`, PBE mit libxc-μ/β wie WP6; curcuma 250x1454 mit Antwort; ORCA-Ausgabe nennt keine
Gitterantwort-Einstellung, „Difference to translation invariance“ 0):

| | LSD | PBE | B3LYP |
|---|---|---|---|
| H2O | 4.4e-7 | 1.3e-6 | 7.8e-7 |
| NH3 | 1.8e-7 | 3.5e-7 | 1.7e-7 |
| CH4 | 6.0e-8 | 1.0e-7 | 6.1e-8 |
| HF | 4.1e-7 | 4.7e-7 | 5.0e-7 |

Einordnung: curcuma-Gitterkonvergenz H2O (gegen 300x2030): PBE 250x590 5.9e-7, 250x770 7.4e-8, 250x1454 3.1e-8;
ORCAs feinstes Gitter ist Lebedev-770 (H: 590, OptM3-Radialgitter, mBecke-Gewichte), der Rest ist also von der Größe
von ORCAs Gitterfehler; mit/ohne curcuma-Antwort bei 250x1454 gleich. Nicht bewiesen (ORCA-Gitter nicht nachgebaut).

**(4) Translation, Rotation** (B3LYP, Gradient der gedrehten Geometrie zurückgedreht; Gitter dreht nicht mit):

| | Σ g (mit / ohne Antwort) | ΔE Rotation | max Δg Rotation (mit / ohne) | Drehmoment (mit) |
|---|---|---|---|---|
| H2O sg1 / sg2 / sg3 | 1e-14 / 1.8e-5, 2.3e-6, 6.4e-7 | 2.9e-6 / 2.2e-7 / 5.4e-9 | 2.2e-5/2.3e-5, 2.5e-6/2.3e-6, 2.8e-7/5.5e-7 | ≤ 1.8e-5 |
| NH3 sg1 / sg2 / sg3 | 4e-14 / 1.5e-5, 3.2e-6, 9.2e-8 | 5.7e-6 / 2.3e-7 / 2.9e-9 | 6.4e-5/1.7e-5, 1.2e-5/1.6e-6, 1.1e-7/2.0e-7 | 2.1e-5 / 1.4e-5 / 1.7e-7 |
| Benzol sg1 / sg2 / sg3 | 6e-14 / 6e-14 (Symmetrie) | 5.2e-5 / 9.3e-7 / 6.2e-7 | 1.0e-4/4.0e-5, 2.2e-5/2.3e-5, 5.5e-6/1.6e-6 | 1.5e-4 / 1.4e-4 / 2.9e-5 (gedreht) |

Translation: Gradient unverändert ≤ 3.8e-13. Drehmoment ≠ 0 ist echt (die Energie hängt von der Orientierung ab).

**(5) `-opt`** (`-opt.convergence_preset tight`, sg1, Start `qm_grad/{H2O,formaldehyde}_distorted.xyz`; Referenz
scipy-BFGS auf PySCF-Energie + Gradient mit Antwort, gtol 1e-7):

| | ΔE (curcuma-Min − PySCF-Min) | max Δd_ij (Å) | max\|g\| am Ende | curcuma am PySCF-Minimum ΔE / max\|g\| |
|---|---|---|---|---|
| H2O LDA / PBE / B3LYP | 1.9e-11 / 1.1e-7 / 3.2e-9 | 5.7e-6 / 3.9e-5 / 1.8e-6 | 2.6e-6 / 2.6e-6 / 1.7e-6 | ≤ 1e-11 / ≤ 8.5e-8 |
| H2CO LDA / PBE / B3LYP | 5.5e-11 / 2.4e-11 / 8.6e-11 | ≤ 9.5e-6 | 2.0e-6..4.2e-6 | ≤ 2.4e-11 / ≤ 2.0e-7 |

H2O PBE: das curcuma-Minimum liegt in einer anderen Orientierung (Gitter nicht rotationsinvariant), dort ist die
Energie 1.1e-7 höher; die PySCF-Geometrie ist auch für curcuma stationär (8.5e-8). Gate daher 1e-6 Eh plus Test am
PySCF-Minimum.

**(6) Threads, UpdateMolecule:** 1 vs 4 Threads bitweise gleich (LDA/PBE/B3LYP NH3, B3LYP Benzol, RI-J); RI-JK nicht
bitweise (die RI-JK-SCF-Energie selbst weicht schon um 1e-14 ab, vorbestehend). UpdateMolecule vs neue Instanz
max Δg 3.1e-13 (B3LYP), 3.3e-13 (PBE).

Energien unverändert: H2O HF -75.9575137984; `ctest -L "qm_1e|qm_2e|qm_grad|qm_hf3c|qm_grid|qm_lda|qm_pbe|qm_b3lyp|qm_ri|qm_mp2|qm_dft_grad"`
340/342 (die zwei `qm_opt_hf3c_*` brauchen numpy im System-Python).

## Nicht getestet / nicht implementiert

- Nicht implementiert: offene Schalen, Z > 10, GPU-KS-Gradient (2e-Teil für KS immer CPU; das Device-2e-Kernel hat die
  RHF-Gewichte), Hessian, Dispersion, Gitter-Pruning, SSF-Screening. Meta-GGA-Gradient: WP11 Stufe 2
  (`WP11_meta_gga.md`).
- Nicht getestet: `-md` mit KS, Stratmann-Partition und `radii_adjust becke|none` im Gradienten (Code vorhanden, nur
  über die gemeinsame Formel), `euler_maclaurin`-Radialgitter, kartesische d (`-qm.cartesian_d`), integral-direct und
  Device-1e-Pfad mit KS, Moleküle über Benzol, geladene Systeme, Gitter-Rotationsinvarianz über 3 Orientierungen hinaus.
- `-qm.xc_density_threshold` macht die Energie an der Schwelle unstetig (Punkte fallen weg); im Gradienten nicht
  sichtbar (≤ 1e-10), aber nicht systematisch untersucht.
