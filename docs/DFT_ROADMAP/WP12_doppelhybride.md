# WP12 — Doppelhybride (B2PLYP, B2GP-PLYP, DSD-*, revDSD-PBEP86)

- **Status Stufe 1 (Energien, CPU):** 🤖 AI-generated, ⚙️ machine-tested (`ctest -L qm_dh`, 31. Okt 2026); menschlicher
  Produktionstest ausstehend. **Stufe 2 (analytischer Gradient, CPU):** 🤖 AI-generated, ⚙️ machine-tested (`ctest -L
  qm_dh_grad`, 31. Okt 2026), menschlicher Produktionstest ausstehend. **Stufe 3 (CUDA):** 🤖 AI-generated, ⚙️ machine-tested
  nur auf einer RTX 5090 (`ctest -L qm_gpu_dft`, Fälle `*dh*`, 31. Okt 2026), menschlicher Produktionstest ausstehend; Abschnitt
  „Stufe 3“ unten. **Stufe 4 (D4, CPU):** 🤖 AI-generated, ⚙️ machine-tested (`ctest -L qm_d4`, 31. Okt 2026), menschlicher
  Produktionstest ausstehend; Abschnitt „Stufe 4“ unten.
- **Verhaltensänderung 31. Okt 2026** (Betreiberentscheidung, Abschnitt „Entscheidungen vom 2026-10-31“): Defaults
  `-qm.p86_lda vwn5` (vorher `pz81`) und frozen core für alle Doppelhybride (vorher aus); KS- und PT2-Zweielektronen-
  ableitung in einem Durchlauf; MP2-Gradient über Threadzahlen bitgleich. Energien von `dsd-pbep86*`/`revdsd-pbep86*`
  ändern sich um 1.3e-3..5.8e-3 Eh (P86), die aller Doppelhybride durch frozen core; ältere Zahlen unten nennen ihre
  Einstellungen.
- Name: „WP12“ stand in der README-Tabelle für die gestrichene Bewertung der Schwarz-Ansätze; diese Datei ist das neue WP12.
- **Code:** `qm_methods/double_hybrid_method.{h,cpp}` (Tabelle `doubleHybridTable()`, Klasse `DoubleHybridMethod`, Parameter
  `-qm.dispersion`, `-qm.dh_pt2`, `-qm.p86_lda`), Kerne `P86_C`, `P86VWN_C`, `B95_C` in `qm_xc_kernels.h`,
  `QMFunctional::DoubleHybrid` + `QMEngine::setXCFunctional` (`qm_engine.*`), PT2 = `computeRIMP2`/`computeMP2` (nur die
  HF-Abfrage erweitert), D3(BJ) = `D3ParameterGenerator` (der von hf-3c), Zitate `qmcite::doubleHybrid`.
  Tests `test_cases/qm_dh/`, Referenzen `scripts/qm_dh_pyscf_reference.py` (PySCF 2.14, libxc 7.0.0),
  `scripts/qm_dh_orca_compare.py` (ORCA 6.1.0).
- **Code Stufe 2:** `QMEngine::computeDoubleHybridGradient` (`qm_dh_grad.cpp`: RI-PT2-Zwischengrößen, KS-Orbitalhessian
  `dhG4`, Z-Vektor, Kräfte), `qmxc::buildFxc`/`buildXCResponseGradient` (`qm_xc_response.cpp`), Jet-Kerne `qm_xc_fxc.h` (seit Stufe 3 `qm_xc_fxc_kernels.h`, Host und Device),
  `qmint::gradientTwoElectronDHPT2` (`qm_integrals.cpp`, deklariert in `qm_ri.h`), `RIJK::fitTransposeRows`/`auxNormalisation`,
  D3(BJ)-Gradient wie hf-3c in `DoubleHybridMethod::calculateEnergy`, Zitat `dh_gradient` (Neese/Schwabe/Grimme 2007).
  Tests `test_cases/qm_dh_grad/`, Referenzen `scripts/qm_dh_grad_pyscf_reference.py`, `scripts/qm_dh_grad_orca_compare.py`,
  `scripts/qm_dh_opt_orca_compare.py`.

## Was gebaut ist

E = E_KS[a_x HF + Σ w_k E_xc,k] + c_os E_os + c_ss E_ss + E_disp: geschlossenschalige Hybrid-KS-SCF (Fock wie WP7), dann die
MP2-Paarsummen von `qm_mp2.h` mit den KS-Orbitalen und -Eigenwerten (keine Einfachanregungen, wie in allen DH-Definitionen),
getrennt nach Spin (os/ss), skaliert. PT2 standardmäßig RI (`-qm.ri_mp2_aux`, Default def2-universal-JKFIT = def2/JK; kein
def2-SVP/C vorhanden: Methodenvariante, Ausgabe sagt es), `-qm.dh_pt2 canonical` = exakte Integrale. Frozen core über
`-qm.frozen_core` (seit 31. Okt 2026 Default an für Doppelhybride, wie ORCA; mp2 und die anderen bleiben aus). Seit Stufe 2 mit analytischem Gradient
(`hasGradient()` true, seit Stufe 4 auch mit D4; `-opt` läuft); `-gpu cuda` (seit Stufe 3): Energien aller Doppelhybride ganz auf dem Device
(J/K, V_xc/E_xc mit P86/B95, resident KS-SCF-Loop, RI-MP2-Summe), kein gezählter Fallback, `-gpu_strict true` läuft durch; im
Gradienten KS-Teil, F_xc des Z-Vektors, XC-Antwortgradient und PT2-1e/2e/RI-3-Index-Ableitungen auf dem Device, die
PT2-Zwischengrößen (Amplituden, Lagrangian, Z-Vektor-Algebra) auf dem Host als gezählter Fallback `qm_dh_grad_pt2`.

## Parameter und Quellen

| Methode | a_x | DFT-X | DFT-C | c_os | c_ss | Dispersion | Quelle der Zahlen |
|---|---|---|---|---|---|---|---|
| `b2plyp` | 0.53 | 0.47 B88 | 0.73 LYP | 0.27 | 0.27 | – | Grimme 2006; = ORCA 6.1, Psi4 |
| `b2plyp-d3bj` | 0.53 | 0.47 B88 | 0.73 LYP | 0.27 | 0.27 | s6 0.64, a1 0.3065, s8 0.9147, a2 5.0570 | D3(BJ): Goerigk/Grimme 2011 (s-dftd3-Tabelle) |
| `b2gp-plyp` | 0.65 | 0.35 B88 | 0.64 LYP | 0.36 | 0.36 | – | Karton et al. 2008; = ORCA, Psi4 |
| `b2gp-plyp-d3bj` | 0.65 | 0.35 B88 | 0.64 LYP | 0.36 | 0.36 | s6 0.56, a1 0, s8 0.2597, a2 6.3332 | Goerigk/Grimme 2011 |
| `dsd-blyp` | 0.75 | 0.25 B88 | 0.53 LYP | 0.46 | 0.60 | – | Kozuch/Martin 2013 „ohne Dispersion“; = ORCA `DSD-BLYP`, Psi4 |
| `dsd-blyp-d3bj` | 0.71 | 0.29 B88 | 0.54 LYP | 0.47 | 0.40 | s6 0.57, a1 0, s8 0, a2 5.4 | nur Psi4 (eine Quelle) |
| `dsd-pbep86` | 0.72 | 0.28 PBE | 0.44 P86 | 0.51 | 0.36 | – | KM2013; = ORCA, Psi4 |
| `dsd-pbep86-d3bj` | 0.69 | 0.31 PBE | 0.44 P86 | 0.52 | 0.22 | s6 0.48, a1 0, s8 0, a2 5.6 | Psi4; D3 auch s-dftd3 |
| `dsd-pbeb95-d3bj` | 0.66 | 0.34 PBE | 0.55 B95 | 0.46 | 0.09 | s6 0.61, a1 0, s8 0, a2 6.2 | Psi4; D3 auch s-dftd3 |
| `revdsd-pbep86` | 0.69 | 0.31 PBE | 0.4224 P86 | 0.5935 | 0.0566 | – | Mehta/Santra/Martin 2023 Tab. 1 (PZ81-Zeile); = ORCA `revDSD-PBEP86/2021` |
| `revdsd-pbep86-d4` | wie oben | | | | | D4(EEQ)-ATM s6 0.5917, s8 0, a1 0.371, a2 4.201, s9 1 | Mehta et al. 2023 Tab. 1 (PZ81-Zeile, gelesen 31. Okt 2026, arXiv 2211.03105); seit Stufe 4 |
| `b2plyp-d4` | 0.53 | 0.47 B88 | 0.73 LYP | 0.27 | 0.27 | D4(EEQ)-ATM s6 0.64, s8 1.16888646, a1 0.44154604, a2 4.73114642, s9 1 | dftd4 4.3.1 `parameters.toml` (Fit Caldeweyher 2019); = cpp-d4, = ORCA 6.1 `B2PLYP D4` |
| `b2gp-plyp-d4` | 0.65 | 0.35 B88 | 0.64 LYP | 0.36 | 0.36 | D4(EEQ)-ATM s6 0.56, s8 0.94633372, a1 0.42907301, a2 5.18802602, s9 1 | dftd4 4.3.1 (Caldeweyher 2019); = cpp-d4, = ORCA 6.1 `B2GP-PLYP D4` |

- Kozuch/Martin 2013 (JCC 34, 2327) und Santra et al. 2019 (JPCA 123, 5129) sind nicht gelesen (nicht zugänglich); die
  DSD-Zahlen stammen aus ORCA 6.1 (gedruckte ScalHFX/ScalDFX/ScalDFC/„Scaling for aa/bb, ab pairs“) und Psi4
  `dh_functionals.py`, die revDSD-Zahlen aus Mehta et al. (arXiv 2211.03105, Can. J. Chem. 101, 656 (2023)).
- **Weggelassen:** `dsd-pbeb95` ohne Dispersion (ORCA 0.70 HF / 0.30 PBE / 0.52 B95, Psi4 0.69 / 0.31 / 0.54, beide 0.48 / 0.22,
  beide zitieren KM2013); `revdsd-pbep86-d3bj` (s6 0.4339 bei Mehta et al., 0.4377 in s-dftd3; c_os/c_ss nur eine Quelle).
- Widersprüche gefunden: die D4-Parameter von revDSD-PBEP86 in der dftd4-Tabelle (s6 0.5132, a1 0.44, a2 3.60, DOI des
  2019er Papers) weichen von Mehta et al. ab (gehören zu den XC/PT2-Zahlen von 2019, nicht zu dieser Zeile); ORCA 6.1
  `revDSD-PBEP86-D4/2021` rechnet a2 = 4.2014 statt 4.201 (Stufe 4); die s-dftd3-Einträge `dsdblyp` (s6 0.50, s8 0.2130, a2 6.0519) gehören zum DSD-BLYP von 2010, nicht zu KM2013.
- **P86-Lokalteil:** libxc `GGA_C_P86` (Perdew 1986, Gaussian, Q-Chem, Psi4) setzt den Gradiententerm auf PZ81, ORCA (kein
  PZ81) auf VWN5 (Ausgabe `LDAOpt VWN-5`). Default seit 31. Okt 2026 `-qm.p86_lda vwn5` (vorher `pz81`); `pz81` = libxc/PySCF.
  Unterschied siehe Tabelle unten (1.3e-3..5.8e-3 Eh). **Widerspruch:** die revDSD-Parameter wurden mit PZ81 gefittet (Mehta
  et al.); der Default nutzt sie wie ORCAs `revDSD-PBEP86/2021` auf VWN5 (Mehta et al. geben für VWN5 eigene: 0.4200 / 0.6090 /
  0.0476, hier nicht eingebaut); `revdsd-pbep86*` warnt bei verbosity >= 1, `-qm.p86_lda pz81` rechnet die gefittete Definition.
  Für DSD-PBEP86 (KM2013, Paper nicht gelesen) ist nicht geprüft, auf welchem Lokalteil gefittet wurde.
- ORCA-Falle: `%mp2 PS` ist opposite spin, `PT` same spin (Ausgabe „aa/bb pairs“ = PT).
- D3(BJ): zweikörperig (s9 = 0), wie die DSD-Fits (Psi4 „d3bj2b“) und ORCAs D3BJ (nur E6/E8).

## Neue Kerne (`qm_xc_kernels.h`)

- `P86_C`: PZ81 (Perdew/Zunger 1981, libxc `LDA_C_PZ`, nicht `_MOD`) + e_H = σ ρ^(-4/3) C(r_s) e^(-Φ), Φ = f̃ C_∞/C |∇ρ| ρ^(-7/6)
  (Perdew 1986 Gl. 6, 8, 9; f̃ = 1.745·0.11 wie libxc `GGA_C_P86`, nicht `_FT`); Ableitungen von Hand, v_σ endlich bei σ = 0.
  `P86VWN_C`: derselbe Term auf VWN5 (libxc `GGA_C_P86VWN`). Auf dem Host `noinline` (sonst änderte GCCs Inlining die
  B3LYP-Gradientenbits).
- `B95_C`: libxc `mgga_c_bc95.mpl` (Becke 1996), Stoll-Zerlegung mit PW92 (`LDA_C_PW_MOD`, ζ = 0 und ζ = 1), g_par =
  (t_s − x_s²/8)/(K_F (1 + c_ss x_s²)²), g_perp = 1/(1 + 2 c_ab x_s²), c_ss = 0.038, c_ab = 0.0031; Ableitungen mit den
  Dual-Zahlen von WP11, τ-Boden und σ ≤ 8ρτ wie die Meta-GGAs.

Gegen libxc (`ctest qm_dh_kernel`; 66 (ρ, σ)-Punkte ρ 1e-12..1e2 × s 0..30, B95 zusätzlich α 0..1e3 und 110 Punkte τ < τ_W;
skalierter Fehler wie WP11):

| Kern | ε | v_ρ | v_σ | v_τ | v vs FD (ρ ≥ 1e-8) |
|---|---|---|---|---|---|
| P86 (PZ81) | 6.8e-16 | 1.6e-15 | 5.6e-15 | – | ≤ 2.4e-12 |
| P86 (VWN5) | 5.4e-14 | 4.8e-14 | 5.6e-15 | – | ≤ 6.8e-12 |
| B95 | 6.4e-16 | 4.7e-15 | 1.2e-15 | 7.9e-16 | ≤ 1.6e-10 (v_ρ bei ρ = 1e-10: 6.6e-6, Stencil-Rauschen) |

## Validierung (31. Okt 2026, def2-SVP sphärisch, SCF 1e-10 / dE 1e-12 / |g| 1e-9)

**(1) gegen PySCF 2.14** (kein DH-Funktional in PySCF/libxc: Komponenten; KS mit `xc = 'a_x*HF + w*GGA_X_..., c*GGA_C_...'`
auf identischem Gitter, dann RI-PT2 mit derselben JKFIT-Datei, `cholesky_eri`, bzw. `mp.MP2(mf)`; max |Δ| in Eh über die Fälle,
ctest `qm_dh_energy_*`, `qm_dh_fc_*`, `qm_dh_canonical_*`, Gates E_KS 1e-8, E_xc 1e-7, E_os/E_ss 1e-9):

| Funktional (Moleküle) | Fälle | E_KS | E_xc | E_os | E_ss | SCF-Iter. |
|---|---|---|---|---|---|---|
| b2plyp (H2O, CH4, NH3, HF, LiH, Ne, Benzol; sg1) | 7 | 1.5e-11 | 1.2e-10 | 1.2e-11 | 3.2e-12 | 6-14 |
| b2gp-plyp | 7 | 8.2e-12 | 1.1e-10 | 6.4e-12 | 5.1e-12 | 7-14 |
| dsd-blyp | 7 | 7.0e-12 | 5.5e-11 | 7.0e-12 | 3.2e-12 | 7-15 |
| dsd-pbep86 (PZ81) | 7 | 1.6e-12 | 6.6e-11 | 7.1e-12 | 2.8e-12 | 6-15 |
| revdsd-pbep86 | 7 | 2.5e-12 | 9.0e-11 | 7.9e-12 | 4.5e-12 | 6-14 |
| dsd-pbeb95-d3bj (ohne LiH) | 6 | 2.3e-13 | 9.9e-11 | 4.4e-12 | 4.2e-12 | 6-14 |
| dsd-blyp-d3bj, dsd-pbep86-d3bj, dsd-pbep86 VWN5 (H2O, NH3, Benzol) | 9 | 8.0e-12 | 9.0e-11 | 1.1e-11 | 7.7e-12 | 14-15 |
| b2plyp / dsd-pbeb95-d3bj H2O sg3 | 2 | 3.5e-12 | 7.5e-11 | 3.9e-13 | 1.1e-12 | 13-14 |
| alle 52 mit frozen core | 52 | 1.5e-11 | 1.2e-10 | 1.1e-11 | 7.6e-12 | |
| kanonisches PT2 (H2O, Benzol; 18) | 18 | 1.5e-11 | 7.5e-11 | 1.1e-11 | 8.0e-12 | |

Beispiele E_KS (sg1): H2O b2plyp −76.2190522793, Benzol b2plyp −231.5589833705, Benzol dsd-pbeb95-d3bj −231.2524112476 Eh.

**(2) gegen ORCA 6.1** (`scripts/qm_dh_orca_compare.py`: ORCA-Schlüsselwort, NoRI-SCF, `%mp2 RI true`, AuxC def2/JK, frozen core,
DefGrid3 + IntAcc 8 AngularGrid 7 GridPruning 0, XMuePBE = libxc-μ; Overrides der Skalierung/D3-Parameter für die -d3bj-Zeilen
und dsd-pbeb95-d3bj; curcuma 250x1454, `-qm.frozen_core true`, P86 mit `-qm.p86_lda vwn5`; Differenz curcuma − ORCA in Eh):

| Funktional | H2O | NH3 | CH4 | HF | Ne | Benzol | ΔE_PT2 max | ΔE_D3 max |
|---|---|---|---|---|---|---|---|---|
| b2plyp | 3.3e-8 | -1.7e-10 | 2.3e-9 | -1.1e-8 | -2.5e-10 | -3.9e-7 | 1.5e-8 | – |
| b2gp-plyp | 2.5e-8 | -2.7e-9 | 6.6e-10 | -8.6e-9 | -1.6e-10 | -2.8e-7 | 1.6e-8 | – |
| dsd-blyp | 1.9e-8 | -3.0e-8 | -1.6e-9 | -6.8e-9 | -1.3e-10 | -1.9e-7 | 2.8e-8 | – |
| dsd-pbep86 (VWN5) | 2.1e-8 | -3.1e-8 | -5.1e-9 | 4.9e-9 | 2.1e-10 | 3.2e-7 | 2.4e-8 | – |
| revdsd-pbep86 (VWN5) | 2.3e-8 | -2.3e-8 | -4.9e-9 | 3.5e-9 | -5.7e-11 | 3.0e-7 | 1.6e-8 | – |
| dsd-pbeb95-d3bj | 4.0e-8 | -2.4e-8 | -9.4e-9 | 1.5e-8 | 9.0e-12 | **1.9e-5** | 6.5e-7 (Benzol) | 3.8e-10 |
| b2plyp-d3bj / b2gp-plyp-d3bj / dsd-blyp-d3bj / dsd-pbep86-d3bj (H2O, NH3, Benzol) | ≤ 3.3e-8 | ≤ 2.3e-8 | | | | ≤ 3.9e-7 | 2.0e-8 | 5.6e-10 |

- Kleine Moleküle ≤ 4e-8 Eh, Benzol 2-4e-7: dasselbe Muster wie WP5-WP7 (Gitter, ORCAs Winkelgitter). E_D3 ≤ 5.6e-10 Eh (ORCA
  druckt 9 Stellen); ctest `qm_dh_d3_*` (Gate 1e-9) liest die ORCA-Werte aus `reference_orca.json`.
- **B95 Benzol 1.9e-5:** Winkelgitter. curcuma E_KS dsd-pbeb95-d3bj Benzol (250 radial, frozen core) 590: −231.252343550,
  770: −231.252322207, 974 (150 radial): −231.252316563, 1202: −231.252324536, 1454: −231.252324039, 2030: −231.252322370
  (400x1454 = 250x1454 auf 12 Stellen); ORCA fein −231.2523435468 = curcuma 250x590 auf 3e-9. B95 braucht ≥ 770 Winkelpunkte
  und streut dann noch um 8e-6 Eh (Benzol); kleine Moleküle unauffällig. Ursache der Winkelempfindlichkeit nicht untersucht.

**(3) Weitere Checks (ctest):** 1 vs 4 Threads Gesamtenergie bitgleich (b2plyp, dsd-pbeb95-d3bj); UpdateMolecule A → B = frisch
(≤ 2.8e-14 Eh); Gradientanfrage = Fehler, `hasGradient()` false; `curcuma -sp -method b2plyp` druckt E, `-opt` endet mit „no
analytic gradient“; `revdsd-pbep86-d4` endet mit „D4 ... not available“; `-qm.dispersion d3bj` bei `b2plyp` = Fehler mit Hinweis
auf `b2plyp-d3bj` (seit Stufe 4 rechnet `revdsd-pbep86-d4`, `qm_dh_cli_d4` prüft das); `-gpu cuda` (RTX 5090): gezählter `qm_xc`-Fallback, Energie = CPU (Δ 0 auf 12 Stellen), `-gpu_strict true`
Exit 3. `ctest -L qm_dh` 78/78; `ctest -L "qm_dh|qm_mgga|qm_mgga_grad|qm_lda|qm_pbe|qm_b3lyp|qm_dft_grad|qm_gpu_dft|qm_grid|qm_hf3c|qm_mp2|qm_lmp2|qm_ccsd|qm_ri|qm_cite"`
(-j4) 573/573 (372 s); nach der letzten Umstellung (Tabelle ohne std::string, Zitierfunktion außer Linie) `qm_dh|qm_cite` 86/86.
Laufzeit Benzol b2plyp sg1, 8 Threads: 1.6 s Wand (10 SCF-Iterationen, V_xc 0.53 s, RI-PT2 61 ms inkl. eigenem RI-Setup).

**(4) Regression:** gegen das Binary von 602f727f (vor der Änderung kopiert, gleiche CMake-Optionen), CPU (`-gpu none`),
`-dump_gradient`-md5 + Energiezeilen, H2O_distorted/NH3/Benzol, 1 und 4 Threads: hf, lda (beide VWN), pbe (+ ohne Gitterantwort,
ohne Screening), b3lyp (beide VWN, RI-JK), tpss, scan, mp2/scs-mp2 (+ fc; Gradient), ri-mp2 (+ fc), sos-ri-mp2, ccsd, lmp2, hf-3c:
112/112 bitgleich (erster Lauf: b3lyp-Gradienten 18/18 im letzten Bit anders, Ursache GCC-Inlining des neuen P86-Falls in
`evaluateKernel`; mit `noinline` behoben; nach der letzten Codeänderung erneut 112/112).

## Stufe 2: analytischer Gradient (CPU) — Gleichungen

Geschlossene Schale, reelle kanonische KS-Orbitale C, ε (Fock des Hybrids f = h + J − (a_x/2) K + V_xc). I, J
eingefroren, i, j aktiv besetzt, a, b virtuell, k, l alle besetzt, p, q alle. Herleitung wie der MP2-Gradient
(`qm_mp2_grad.cpp`, Kopfblock); was sich durch die KS-Referenz ändert:

1. **Energie.** E = E_KS + E2 + E_disp, E2 = Σ_ijab (ia|jb) T̃_ijab, t = (ia|jb)/D, D = ε_i + ε_j − ε_a − ε_b,
   T̃ = (c_os + c_ss) t − c_ss t^ba. E_KS ist variationell in C (sein Gradient ist der WP8/WP11-KS-Gradient mit
   W_KS = 2 Σ_k ε_k C_k C_kᵀ und Gitterantwort). E2 ist es nicht: Orbitalantwort über Z-Vektor in der **KS**-Basis.
2. **Unrelaxierte Dichte und Lagrangian** wie MP2: P_ij = −2 Σ T̃_ikab t_jkab, P_ab = 2 Σ T̃_ijac t_ijbc,
   L1_pi = 4 Σ_jab T̃_ijab (pa|jb), L2_pa = 4 Σ_ijb T̃_ijab (ip|jb). RI: (pq|rs) ≈ Σ_R B^R_pq B^R_rs,
   B = L⁻¹ (R|pq) (Cholesky J = L Lᵀ), Γ^R_ia = Σ_jb T̃_ijab B^R_jb, L1_pi = 4 Σ_aR B^R_pa Γ^R_ia,
   L2_pa = 4 Σ_iR B^R_ip Γ^R_ia.
3. **KS-Orbitalhessian.** Die Ableitung von f bei Orbitaldrehung U enthält statt der HF-Matrix
   A_pq,rs = 4(pq|rs) − a_x[(pr|qs) + (ps|qr)] + 4 (pq|f_xc|rs), also
   G4(X) = Cᵀ [4 J(X̂) − 2 a_x K(X̂) + 4 F_xc[X̂]] C, X̂ = C X Cᵀ, mit dem XC-Kern (zweite Funktionalableitung,
   geschlossene Schale, Singulett) F_xc[D]_mn = ∫ w [δv_ρ φ_mφ_n + 2 δv_σ ∇ρ·∇(φ_mφ_n) + 2 v_σ ∇ρ_D·∇(φ_mφ_n)
   + ½ δv_τ ∇φ_m·∇φ_n], ρ_D = Σ D φφ, σ_D = 2 ∇ρ·∇ρ_D, τ_D = ½ Σ D ∇φ·∇φ, δv_x = Σ_y ∂²e/∂x∂y · y_D (x, y ∈ ρ, σ, τ).
   Damit unverändert gegenüber MP2: L_ak = L1_ak[k aktiv] − L2_ka + G4(P + W_fz)_ak, w_Ij = L1_Ij/(ε_j − ε_I),
   (ε_a − ε_k) z_ak + G4(Z)_ak = −L_ak (PCG), P2 = P + w/2 + z/2, Ws wie MP2 mit diesem G4.
4. **Kraft.** dE2/dX = Σ P2 h^X − Σ (−C Ws Cᵀ) S^X + ½ Σ D2 (mn|ls)^X + [XC-Term] + [RI-Term], mit
   D2 = P2_mn P_ls + P_mn P2_ls − (a_x/4)(P2_ml P_ns + P_ml P2_ns + P2_ms P_nl + P_ms P2_nl) (+ 4 Γ^sym kanonisch).
   XC-Term = Σ P2_mn ∂V_xc,mn/∂X bei festen AO-Dichten = d/dt [dE_xc[P + t P2]/dX]_{t=0}: die Richtungsableitung
   des KS-XC-Gradienten (Basis-, Punkt- und Becke-Gewichtsteil) in Richtung P2. Sie braucht nur die **zweite**
   Funktionalableitung (Potentiale v_x an P, δv_x in Richtung P2, Energiedichte der Gewichtsableitung
   δe = v_ρ ρ2 + v_σ σ2 + v_τ τ2); eine dritte Ableitung tritt im Z-Vektor-Weg nicht auf (sie käme erst bei
   Kontraktionen mit zwei Antwortdichten, z. B. TDDFT-Gradient oder Hessian). Ohne `-qm.grid_response` nur der Basisteil.
   RI-Term = Σ_P Σ_mn G^P_mn (mn|P)^X − Σ M2_PQ (P|Q)^X, G = 2(H + Hᵀ), H^P = C_o c^P C_vᵀ, c = J⁻¹-Koeffizienten von Γ
   (D Fᵀ Γ im Rohbasis-Aux), M2 = 2 sym(Σ_ia d^P_ia c^Q_ia), d = J⁻¹ (P|ia) (`RIJK::threeIndexGradient`,
   `metricGradient`). Kanonisch: 4 Γ^sym im 2e-Treiber wie MP2.
5. **Dispersion:** D3(BJ)-Gradient des hf-3c-Codes (Paar + CN-Kette).
6. **Kerne f_xc:** zweite Ableitungen per Vorwärts-AD zweiter Ordnung (Jet: Wert, 3 Gradienten, 6 Hessian-Einträge)
   auf eigenen, nur-Energie-Templates der Kerne B88, LYP, PBE-X, PZ81, VWN5, P86-Gradiententerm, PW92, B95
   (`qm_xc_fxc.h`, seit Stufe 3 `qm_xc_fxc_kernels.h`); Begründung: von Hand gesetzte zweite Ableitungen von LYP/P86/B95 sind fehleranfällig, AD ist exakt
   bis Rundung. Die Produktionskerne (Energie, V_xc) bleiben unberührt; die Jet-Ableitungen erster Ordnung werden gegen
   sie geprüft (ctest), die zweiten gegen libxc (PySCF `eval_xc(deriv=2)`) und gegen FD der Produktionspotentiale.

## Stufe 2: Validierung (31. Okt 2026, def2-SVP sphärisch, SCF 1e-10 / dE 1e-12 / |g| 1e-9, Eh/Bohr)

**Kerne** (`ctest qm_dh_grad_kernel`, 105 Punkte ρ 1e-6..1e2 × s 0.01..10 × α 0.1..5): Jet-Ableitungen 1. Ordnung gegen die
Produktionspotentiale ≤ 4.6e-15 relativ; 2. Ableitungen gegen libxc 7.0.0 (`eval_xc` deriv 2: v2rho2, v2rhosigma,
v2sigma2, bei B95 auch v2rhotau, v2sigmatau, v2tau2) B88 4.1e-14, LYP 1.9e-14, PBE-X 5.3e-15, P86 (PZ81) 6.1e-14,
P86 (VWN5) 6.7e-14, B95 2.5e-15 (relativ, Skala max(|ref|, 1e-3 max|Zeile|)); gegen FD der Produktionspotentiale ≤ 2.3e-12
nach Abzug des 10-fachen Stencil-Rundungsrauschens. **F_xc und XC-Richtungsgradient** gegen FD (h 1e-4) von V_xc[P + tD] bzw.
des XC-Gradienten (`qm_dh_grad_fxc_*`, H2O, relativ): b2plyp 7.4e-11 / Basis 1.0e-9 / Punkte 2.6e-9 / Gewichts-Energiedichte
1.7e-10; dsd-pbep86 (PZ81) 1.8e-10 / 1.2e-9 / 1.7e-9 / 1.8e-10; dsd-pbeb95-d3bj 2.9e-10 / 4.0e-9 / 1.7e-9 / 2.4e-10.

**(1) FD der eigenen Gesamtenergie** (4-Punkt, h = 1e-3 Bohr in Å, alle Atome; mitwanderndes Gitter vs Gradient mit
Gitterantwort, eingefrorenes Gitter (`freezeXCGrid`) vs ohne; H2O, NH3, HF, CH4 × RI/kanonisch × mit/ohne frozen core,
sg1; max |g − FD| über die 4 Moleküle; 416 Läufe):

| Funktional | mitw. RI | mitw. RI+fc | mitw. kan. | mitw. kan.+fc | eingefr. RI | eingefr. RI+fc | eingefr. kan. | eingefr. kan.+fc | Z-Iter. |
|---|---|---|---|---|---|---|---|---|---|
| b2plyp | 1.3e-10 | 1.3e-10 | 1.3e-10 | 1.3e-10 | 1.3e-10 | 1.3e-10 | 1.3e-10 | 1.3e-10 | 7-9 |
| b2plyp-d3bj | 1.3e-10 | 1.3e-10 | 1.3e-10 | 1.3e-10 | 1.3e-10 | 1.3e-10 | 1.3e-10 | 1.3e-10 | 7-9 |
| b2gp-plyp | 1.3e-10 | 1.3e-10 | 1.2e-10 | 1.3e-10 | 1.2e-10 | 1.2e-10 | 1.2e-10 | 1.2e-10 | 8-9 |
| b2gp-plyp-d3bj | 1.3e-10 | 1.3e-10 | 1.3e-10 | 1.3e-10 | 1.2e-10 | 1.1e-10 | 1.2e-10 | 1.2e-10 | 8-9 |
| dsd-blyp | 1.1e-10 | 1.2e-10 | 1.1e-10 | 1.2e-10 | 1.6e-10 | 1.6e-10 | 1.6e-10 | 1.6e-10 | 9 |
| dsd-blyp-d3bj | 1.3e-10 | 1.3e-10 | 1.3e-10 | 1.3e-10 | 1.5e-10 | 1.6e-10 | 1.5e-10 | 1.6e-10 | 8-9 |
| dsd-pbeb95-d3bj | 1.2e-10 | 1.1e-10 | 1.2e-10 | 1.2e-10 | 1.2e-10 | 1.2e-10 | 1.2e-10 | 1.2e-10 | 8-10 |
| dsd-pbep86 (VWN5) | 1.6e-10 | 1.7e-10 | 1.7e-10 | 1.7e-10 | 1.6e-10 | 1.6e-10 | 1.6e-10 | 1.6e-10 | 8-9 |
| dsd-pbep86-d3bj (VWN5) | 1.3e-10 | 1.2e-10 | 1.3e-10 | 1.2e-10 | 1.6e-10 | 1.6e-10 | 1.6e-10 | 1.6e-10 | 8-10 |
| revdsd-pbep86 (VWN5) | 1.5e-10 | 1.5e-10 | 1.4e-10 | 1.5e-10 | 1.3e-10 | 1.3e-10 | 1.3e-10 | 1.3e-10 | 8-10 |
| dsd-pbep86 (PZ81, Default bis 31. Okt 2026) | 3.3e-5 | 3.3e-5 | 3.3e-5 | 3.3e-5 | 3.3e-5 | 3.3e-5 | 3.3e-5 | 3.3e-5 | 8-9 |
| dsd-pbep86-d3bj (PZ81) | 3.7e-5 | 3.7e-5 | 3.7e-5 | 3.7e-5 | 5.4e-5 | 5.4e-5 | 5.4e-5 | 5.4e-5 | 8-10 |
| revdsd-pbep86 (PZ81) | 4.3e-5 | 4.3e-5 | 4.3e-5 | 4.3e-5 | 4.7e-5 | 4.7e-5 | 4.7e-5 | 4.7e-5 | 8-10 |

- **PZ81-Naht:** PZ81 (Perdew/Zunger 1981, libxc `LDA_C_PZ`, der Lokalteil von `GGA_C_P86`) ist bei r_s = 1 unstetig
  (γ/(1 + β1 + β2) = −0.059632 gegen B + D = −0.0596, Δε ≈ 3e-5 Eh; libxc `LDA_C_PZ_MOD` korrigiert das, ist aber nicht die
  P86-Definition). Die DH-Energie springt deshalb, wenn ein Gitterpunkt die Fläche ρ = 3/(4π) überquert: Scan H2O
  dsd-pbep86, H-Atom in z, eingefrorenes Gitter, Schritt 1e-4 Bohr: dritte Differenzen 1-3e-9 Eh an Einzelstellen (Sprung-
  muster), mit VWN5-P86 ≤ 9e-13 überall. Die FD-Abweichung der PZ81-Zeilen (h 1e-3: 3.3e-5, h 2e-4: 1.9e-5) ist dieser Sprung,
  kein Gradientenfehler: dieselbe Zahl gegen die FD der PySCF-Energie (3.31e-5, unten), gleiche Größe für RI/kanonisch/fc und
  mitwanderndes/eingefrorenes Gitter, VWN5-P86 ≤ 1.7e-10. Der analytische Gradient ist die Ableitung des glatten Teils.
  Folge für `-opt`: H2CO mit dsd-pbep86-d3bj (PZ81) bleibt bei |g| ≈ 3.1e-5 stehen (Liniensuche, Exit 1, Geometrie ≤ 1e-9 Eh
  vom Minimum), mit `-qm.p86_lda vwn5` konvergiert es (48 Schritte, |g| 1.4e-6); H2O und NH3 mit PZ81 konvergieren. Deshalb
  seit 31. Okt 2026 der VWN5-Default (Nachweis unten).

**(2) gegen PySCF 2.14** (kein DH-Gradient in PySCF/libxc: 4-Punkt-FD, h 1e-3 Bohr, der PySCF-DH-Energie aus Komponenten wie
Stufe 1, gleiches Gitter, PySCF baut es pro Geometrie neu = mitwanderndes Gitter; curcuma mit Gitterantwort; ctest
`qm_dh_grad_pyscf_*`): H2O b2plyp RI 9.6e-11, RI+fc 1.3e-10, kanonisch 1.5e-10; NH3 dsd-blyp RI 5.9e-10; H2O dsd-pbep86
(VWN5) RI 1.3e-10; HF revdsd-pbep86 (VWN5) RI+fc 1.1e-10; H2O dsd-pbeb95-d3bj RI 4.4e-10 (D3 abgezogen); CH4 b2gp-plyp
kanonisch+fc 1.4e-10; Energien ≤ 9.1e-12 Eh. PZ81: H2O dsd-pbep86 3.31e-5, HF revdsd-pbep86 2.1e-6 (PySCFs Energie hat
dieselben Sprünge; nicht als Gate).

**(3) gegen ORCA 6.1** (`scripts/qm_dh_grad_orca_compare.py`: `! <DH> def2-SVP NoRI VeryTightSCF DefGrid3 EnGrad`, AuxC def2/JK,
`%mp2 RI true`, frozen core, IntAcc 8 AngularGrid 7 GridPruning 0, XMuePBE; curcuma 250x1454, `-qm.frozen_core true`,
P86 auf VWN5, mit Gitterantwort; ob ORCA eine Gitterantwort rechnet, steht nicht in der Ausgabe):

| | b2plyp | dsd-blyp | dsd-pbep86 | b2plyp-d3bj | revdsd-pbep86 |
|---|---|---|---|---|---|
| H2O max \|Δg\| / ΔE | 3.2e-7 / 3.3e-8 | 1.3e-6 / 1.9e-8 | 1.5e-6 / 2.1e-8 | 3.2e-7 / 3.3e-8 | – |
| NH3 max \|Δg\| / ΔE | 1.8e-6 / −1.7e-10 | 3.5e-7 / −3.0e-8 | 7.3e-7 / −3.1e-8 | – | 6.4e-7 / −2.3e-8 |

Größenordnung wie der WP8-KS-Gradient gegen ORCA (4e-7..1.3e-6, Gitter). dsd-pbeb95-d3bj nicht verglichen: ORCA 6.1 bricht ab
(„Gradient is not yet implemented for the PWPB95, PW6B95 and DSD-PBEB95 funtionals“).

**(4) `-opt`** (ctest `qm_dh_opt_*`, sg1, `-opt.convergence_preset tight`): H2O_distorted / NH3 / H2CO_distorted mit b2plyp
und dsd-pbep86-d3bj, H2O mit dsd-pbeb95-d3bj: Endgradient 5.4e-7..3.9e-6 (H2CO dsd-pbep86-d3bj nur mit VWN5-P86, s. o.).
Gegen ORCA-`TightOpt` mit den Einstellungen von (3) (`scripts/qm_dh_opt_orca_compare.py`, curcuma 250x1454 + fc; max |Δd_ij|):
b2plyp H2O 5.2e-6 Å, NH3 5.8e-5 Å (flaches Inversionspotential), H2CO 8.7e-6 Å; Energien gleich auf 6 Stellen (die
curcuma-Ausgabe der xyz-Datei hat 6). dsd-pbep86 (VWN5-P86 wie ORCA) H2O 2.1e-6 Å, H2CO 6.9e-6 Å.

**(5) Weitere Checks:** Σ_A g_A mit Gitterantwort ≤ 4.1e-14 (ohne 5-9e-6, Gitterfehler); verschobene Geometrie max |Δg| ≤ 2.8e-13;
gedrehte Geometrie (Gitter dreht nicht mit) max |Δg| 1.3e-5..2.8e-5 bei sg1 (Gitterfehler, wie WP8); `updateGeometry` gegen
frisch ≤ 1.9e-12; 1 vs 4 Threads: RI-PT2 bitgleich (NH3 b2plyp, Benzol dsd-pbep86-d3bj), kanonisches PT2 damals nur auf 1e-15
(Benzol b2plyp 8.3e-16, NH3 dsd-pbeb95-d3bj 1.3e-16): `mp2GradientHost` (mit `-method mp2` geteilt) teilte die erste
Viertel-Transformation nach Threadzahl auf. Seit 31. Okt 2026 feste Blöcke, 1/4/8 Threads bitgleich (nächster Abschnitt).
`-gpu cuda` (RTX 5090): gezählter `qm_dh_grad`-Fallback, Gradient = CPU bis 4e-14 Eh/Å, `-gpu_strict true` Exit 3.
`ctest -L qm_dh_grad` 66/66 (201 s, -j4).

**(6) Regression:** gegen das Binary von a5596a8c (vor der Änderung kopiert, gleicher Build), CPU, `-dump_gradient`-md5 und
Energiezeilen (`-verbosity 2`), H2O_distorted/NH3/Benzol, 1 und 4 Threads: hf, lda (beide VWN), pbe (+ ohne Gitterantwort,
ohne Screening), b3lyp (beide VWN, RI-JK), tpss, scan, mp2/scs-mp2 (+ fc, Gradient), ri-mp2 (+ fc), sos-ri-mp2, ccsd, lmp2,
hf-3c, Doppelhybrid-Energien (b2plyp + fc + kanonisch, dsd-pbep86, dsd-pbep86-d3bj, dsd-pbeb95-d3bj, revdsd-pbep86 VWN5,
b2gp-plyp-d3bj): 160/160 bitgleich. Nach dem letzten Build (Zitat `dh_gradient`) erneut 160/160. `ctest -L "qm_dh|qm_dh_grad|qm_mgga|qm_mgga_grad|qm_lda|qm_pbe|qm_b3lyp|qm_dft_grad|qm_gpu_dft|qm_grid|qm_hf3c|qm_mp2|qm_cite"` (-j4) 552/552 (453 s).

**Kosten** (Benzol b2plyp sg1, 8 Threads, ruhige Maschine, einmal gemessen, Stand a06631c0): Energie 1.6 s, Energie + Gradient
RI-PT2 8.8 s, kanonisch 9.5 s; davon Z-Vektor (9 Iterationen, je J/K + F_xc) 1.0 s, PT2-Zweielektronenableitung 2.9 s (eigener
Durchlauf neben dem KS-2e-Teil). Seit 31. Okt 2026 ein gemeinsamer Durchlauf (Tabelle im nächsten Abschnitt).

## Entscheidungen vom 2026-10-31 (Betreiber) und ihre Umsetzung

**(a) `-qm.p86_lda` Default `vwn5` (vorher `pz81`).** Grund: PZ81 (libxc `LDA_C_PZ`) ist bei r_s = 1 unstetig (Sprung
~1e-9 Eh pro Gitterpunkt, der die Fläche ρ = 3/(4π) überquert, s. Stufe 2 (1)); `-opt` H2CO dsd-pbep86-d3bj blieb stehen.
VWN5 ist das, was ORCA 6.1 rechnet. Nachweis (`formaldehyde_distorted`, `-qm.scf_threshold 1e-9 -opt.convergence_preset
tight`, 4 Threads): Binary a06631c0 mit seinen Defaults (PZ81, ohne fc) Exit 1 nach 17 Schritten, |g| 3.1e-5; neues Binary mit
seinen Defaults (VWN5, fc) Exit 0 nach 49 Schritten, |g| 9e-6 (Normzeile der Ausgabe); neuer ctest
`qm_dh_opt_formaldehyde_distorted_dsd-pbep86-d3bj_default` (keine P86-/fc-Flags, Endgradient 4.4e-6 Eh/Bohr). Mit dem
normalen Preset konvergierten beide Binaries (12 Schritte). **Widerspruch, offen benannt:** die revDSD-PBEP86-Parameter sind
auf PZ81 gefittet (Mehta et al. 2023); `revdsd-pbep86*` rechnet jetzt per Default ORCAs Variante (dieselben Zahlen auf
VWN5) und warnt (verbosity >= 1); `-qm.p86_lda pz81` gibt die gefittete Definition, Mehtas eigener VWN5-Satz ist nicht
eingebaut (TODO.md, Entscheidungen). Unterschied VWN5 − PZ81 (`test_qm_dh_grad grad`, def2-SVP, sg1, RI-PT2, frozen core,
Gitterantwort, Eh bzw. Eh/Bohr; H2CO = `qm_grad/formaldehyde_hf3c_min_reference.xyz`):

| Funktional | Molekül | E(PZ81) | E(VWN5) − E(PZ81) | max \|Δg\| | max \|g\| (PZ81) |
|---|---|---|---|---|---|
| dsd-pbep86 | H2O | −76.2119985940 | −1.976e-3 | 2.96e-4 | 8.4e-2 |
| dsd-pbep86 | NH3 | −56.3847124298 | −1.338e-3 | 3.41e-4 | 5.0e-2 |
| dsd-pbep86 | H2CO | −114.1895884246 | −2.805e-3 | 4.72e-4 | 1.3e-2 |
| dsd-pbep86 | Benzol | −231.5592803541 | −5.795e-3 | 2.32e-4 | 9.3e-3 |
| dsd-pbep86-d3bj | H2O | −76.2063166108 | −1.976e-3 | 2.97e-4 | 8.5e-2 |
| dsd-pbep86-d3bj | NH3 | −56.3804751086 | −1.338e-3 | 3.41e-4 | 4.8e-2 |
| dsd-pbep86-d3bj | H2CO | −114.1821968626 | −2.806e-3 | 4.72e-4 | 1.2e-2 |
| dsd-pbep86-d3bj | Benzol | −231.5421769606 | −5.798e-3 | 2.32e-4 | 9.8e-3 |
| revdsd-pbep86 | H2O | −76.2026674540 | −1.897e-3 | 2.85e-4 | 8.6e-2 |
| revdsd-pbep86 | NH3 | −56.3783196306 | −1.285e-3 | 3.27e-4 | 4.7e-2 |
| revdsd-pbep86 | H2CO | −114.1767461156 | −2.694e-3 | 4.53e-4 | 1.1e-2 |
| revdsd-pbep86 | Benzol | −231.5182280189 | −5.568e-3 | 2.23e-4 | 1.1e-2 |

Die Gradientendifferenz (2-5e-4 Eh/Bohr) verschiebt Gleichgewichtsgeometrien; nicht gemessen, um wie viel.

**(b) frozen core Default an, nur für Doppelhybride** (wie ORCA). `-qm.frozen_core` ist jetzt ein String-Parameter mit der
Marke `auto` (wie die `scf`-Marken): Doppelhybride lesen `auto` als true, mp2/scs-mp2/sos-mp2/ri-mp2/scs-ri-mp2/sos-ri-mp2/
ccsd/lmp2 als false (Verhalten unverändert); explizit `true`/`false` gilt für jede Methode. `-export_run` schreibt `auto`, ein
Replay behält also den Methoden-Default (geprüft: H2O b2plyp, Replay mit 1 eingefrorenen Orbital, Energie gleich auf 12 Stellen). Eingefroren pro Atom (`qmint::frozenCoreOrbitals`, ladungsunabhängig): H-He 0, Li-Ne
1 (1s), Na-Ar 5, K-Kr 9; im Engine-Umfang (H-Ne) also 1 pro Atom ab Li. Gegen ORCA mit frozen core verglichen nur für H, C, N, O,
F, Ne (Abschnitte (2), (3)); ob ORCA bei Li, Be, B gleich zählt, ist nicht geprüft. ctest `qm_dh_cli_defaults` (H2O
dsd-pbep86: 1 eingefroren und P86 auf VWN5; `-qm.frozen_core false`: 0; mp2: 0; mp2 `-qm.frozen_core true`: 1).

**Tests und Referenzen:** die PySCF-/libxc-Vergleiche (`qm_dh_energy_*`, `qm_dh_fc_*`, `qm_dh_canonical_*`,
`qm_dh_grad_pyscf_*`, `qm_dh_grad_fd_*`) laufen über die Testtreiber, die `p86_lda` und `frozen_core` immer explizit setzen
(Treiber-Default = PySCF-Konvention: pz81, ohne fc; `_vwn5`/`@vwn5`, `_fc`/`--fc` wie bisher); `reference.json` bleibt
unverändert gültig. Die ORCA-Vergleiche setzten schon `-qm.frozen_core true` und VWN5 (entspricht jetzt den Defaults). Die
`qm_dh_opt_*`-Fälle setzen beide Flags jetzt explizit auf der curcuma-Kommandozeile (vorher erbten sie die alten Defaults).

**(c) KS- und PT2-Zweielektronenableitung in einem Durchlauf.** `computeGradient(defer = true)` lässt den CPU-KS-2e-Teil weg,
`qmint::gradientTwoElectronDHPT2(..., include_ks = true)` addiert P_mn P_ls − (a_x/4)(P_ml P_ns + P_ms P_nl) zu den
PT2-Quartettgewichten (separabel P2 × P, kanonisch zusätzlich 4 Γ^sym): ein Durchlauf über die Integralableitungen
(mn|ls)^X für RI- und kanonischen Pfad. **Nicht zusammengelegt:** der RI-PT2-3-Index- und Metrikteil ((mn|P)^X, (P|Q)^X):
andere Integrale, und die KS-SCF eines Doppelhybrids ist exakt (`-qm.ri j|jk` ist ein Fehler), es gibt also keinen
KS-3-Index-Durchlauf, mit dem er zusammengehen könnte; der 1e-Teil bleibt ebenfalls getrennt (60 ms). Läuft der KS-2e-Teil auf
dem Device (`-gpu cuda`, Q-G7), bleibt er dort und der PT2-Teil getrennt auf der CPU (wie vorher); geprüft b2plyp RTX 5090
gegen CPU: H2O_distorted 1.4e-13, Benzol 3.6e-11 Eh/Bohr (a06631c0 ebenfalls 3.6e-11: Quartett-Screen des Device-2e-Teils). Zeiten (Benzol b2plyp sg1,
def2-SVP, `-qm.frozen_core false`, 8 Threads, CPU, ruhige Maschine, je 3 Läufe, Median, Wandzeit des ganzen `-sp` mit
`-dump_gradient`):

| PT2 | vorher (a06631c0) Wand | nachher Wand | vorher KS-2e + PT2-2e | nachher gemeinsamer Durchlauf |
|---|---|---|---|---|
| RI | 8.43 s (8.34 / 8.43 / 8.64) | 5.61 s (5.58 / 5.61 / 5.67) | 2785 + 2928 ms | 2885 ms |
| kanonisch | 9.04 s (8.95 / 9.04 / 9.14) | 6.17 s (6.16 / 6.17 / 6.24) | 2818 + 2869 ms | 2855 ms |

**(d) MP2-Gradient über Threadzahlen bitgleich** (Korrektur des Betreibers an der ersten Fassung des Auftrags, die MP2 unberührt
lassen sollte): `mp2GradientHost` (`qm_mp2_grad.cpp`) schnitt die erste Viertel-Transformation in Zeilenblöcke, deren Höhe aus
der Threadzahl folgte; jetzt feste 8192 Zeilen. Die MP2- und kanonischen DH-Gradienten weichen deshalb im letzten Bit von
a06631c0 ab (Verhaltensänderung, kein bitgleicher MP2-Gradient mehr zu a06631c0). `mp2EnergyHost` (`qm_mp2.cpp`) hat dieselbe
threadabhängige Aufteilung und blieb unverändert (Auftrag: Energien bitgleich); MP2-Energien waren 1/4/8 Threads auf H2O_distorted
und NH3 bitgleich, größere Systeme nicht geprüft.

**Regression** (Binary a06631c0 = `release/` vor der Änderung kopiert, gleiche CMake-Optionen, gegen das neue; CPU `-gpu none`,
`-verbosity 2` Energiezeilen + md5 bzw. Zahlen des `-dump_gradient`; H2O_distorted, NH3, Benzol (ohne mp2/scs-mp2/ccsd/kanonisch),
1 und 4 Threads; Doppelhybride mit den alten Defaults explizit `-qm.frozen_core false -qm.p86_lda pz81`):
hf, lda (2), pbe (3), b3lyp (3 inkl. RI-JK), tpss, scan, hf-3c Gradienten und ri-mp2 (+fc), sos-ri-mp2, ccsd, lmp2, mp2,
scs-mp2, mp2 fc Energien: 112/112 bitgleich (md5 und Energiezeilen). mp2, mp2 fc, scs-mp2 Gradient (12 Läufe): Energiezeilen
gleich, Gradient max |Δg| 5.0e-15 Eh/Bohr; Doppelhybrid-Gradienten (b2plyp RI/RI+fc/kanonisch, dsd-pbep86, dsd-pbep86-d3bj,
dsd-pbeb95-d3bj, revdsd-pbep86 VWN5, b2gp-plyp-d3bj kanonisch; 44 Läufe): Energiezeilen gleich, max |Δg| 5.7e-14 Eh/Bohr (Gate
1e-12). **Threads** (neues Binary, 1/4/8, md5 + Energiezeilen): mp2, scs-mp2, mp2 fc, b2plyp RI/kanonisch, dsd-pbep86-d3bj
RI/kanonisch (H2O_distorted, NH3; RI-DH auch Benzol) und ri-mp2-/mp2-Energien: 23/23 bitgleich; das alte Binary: 11/23 (alle 12
Abweichungen kanonische Gradienten, 1 Thread gegen 4 = 8). `qm_dh_grad_threads_*` prüft jetzt 1/4/8 Threads auf Bitgleichheit,
auch kanonisch.

## Stufe 3: CUDA (31. Okt 2026)

- **Code:** keine neue Kernformel. P86 (PZ81/VWN5) und B95 laufen über `evaluatePoint`/`evaluatePointMGGA` aus
  `qm_xc_kernels.h` durch den GGA- bzw. Meta-GGA-Pfad von `cuda/qm_gpu_xc.cu`; entfernt ist nur die CPU-only-Sperre in
  `QMEngine::ensureGpuXC` (`qm_scf.cpp`) samt `isCpuOnlyKernel`. Gradient: `qm_dh_grad_gpu.cpp` (Device-Dispatch; in
  `qm_dh_grad.cpp` nur Haken), `QMGpuXC::fxc`/`responseGradient` mit `densityPairKernel`, `respRowKernel`, `respColKernel`
  (`cuda/qm_gpu_xc.cu`), Jets in `qm_xc_fxc_kernels.h` (aus `qm_xc_fxc.h` herausgelöst, `QM_HD`, Host-Arithmetik
  unverändert), Setup-Flag `QMGpuXCSetup::fxc` (angehängt, nur für Doppelhybride gesetzt), Plugin-Eintrag umbenannt zu
  `curcuma_cuda_create_qm_xc_fxc`.
- **Aufwand/Nutzen pro Baustein** (Auftrag (a)-(d)): (a) KS-Anteil inkl. XC-Gradient: fällt mit P86/B95 auf dem Device
  von selbst ab (vorhandene Q-G7/Q-G8- und XC-Gradientenpfade). (b) F_xc im Z-Vektor: die Jets waren reiner C++-Code
  ohne Container, also device-tauglich; ein Kernel mit zwei Dichten pro Block plus die V_xc-Montage des vorhandenen Pfads.
  (c) XC-Richtungsableitung: Basis- und Punktteil als eigener T_k-Kern, der Gewichtsteil über die vorhandene
  Becke-Ableitung mit der Richtungs-Energiedichte. (d) PT2-2e: ein GPU-MP2-Gradient existiert (`cuda/qm_gpu_mp2_grad.cu`),
  hat aber a_x = 1 fest; statt eines neuen Kerns die Polarisation des Q-G7-Kerns (kx = a_x): separabler Teil =
  1/2 [G(P + P2) − G(P − P2)] (zwei Läufe), kanonisches Γ über den MP2-Kern mit P = P2 = 0, RI-3-Index über
  `cuda/qm_gpu_ri_grad.cu`, PT2-1e über Q-G8. **Nicht portiert:** die PT2-Zwischengrößen (RI: B für alle MO-Paare,
  (ia|jb), Amplituden, P_ij/P_ab, Γ, L1/L2, Rücktransformation; kanonisch `mp2GradientHost`) und die Z-Vektor-Algebra:
  gezählter Fallback `qm_dh_grad_pt2`. Der RI-Metrikteil (P|Q)^X bleibt Host-Code wie im RI-JK-Gradienten (nicht gezählt).
- **Fallbacks danach:** Energie: keiner (`-gpu_strict true` läuft durch, ctest `qm_gpu_dft_scf_*dsd*`, `qm_dh_gpu_H2O`).
  Gradient: genau `qm_dh_grad_pt2`; Device-Fehler je Teil: `qm_dh_xc` (F_xc/Antwortgradient), `qm_grad1e`, `qm_grad`,
  `qm_mp2_grad`, `qm_ri_grad`. `-gpu_strict true` mit Gradient endet mit Exit 3 (ctest
  `qm_gpu_dft_fallback_dh_h2o_b2plyp_grad_strict`).

GPU gegen CPU auf identischem Gitter (RTX 5090, def2-SVP, `check_qm_gpu_dft.py`, 83 Fälle; Tabellen und Funktionale:
[QM_GPU_ROADMAP.md 4.3](../QM_GPU_ROADMAP.md#43-double-hybrids-on-the-device-oct-31-2026-wp12-stage-3)): V_xc gleiche Dichte
max rel ≤ 2.1e-15, E_xc ≤ 6.4e-15 relativ; SCF H2O/NH3/Benzol (sg1: alle 10 Funktionale, 150x974: 4) E_KS gleich auf
10 Stellen, Gesamtenergie (12 Stellen) ≤ 1.0e-12 Eh, Iterationen gleich bis auf Benzol 150x974 (±1); Koffein E_KS gleich,
Gesamt ≤ 1.5e-11 Eh, 2-3 Iterationen auseinander (langsam konvergierende SCF, wie WP11 Koffein TPSS 30/32). Gradient
gleiche Dichte (H2O, NH3, Benzol, Koffein; sg1, 150x974): F_xc rel ≤ 3.8e-15, XC-Antwortgradient ≤ 1.5e-16, PT2-2e ≤ 7.9e-12,
RI-3-Index ≤ 2.6e-15 Eh/Bohr; ganze Läufe GPU gegen CPU ≤ 8.4e-12 (Benzol), Koffein 1.6e-10 Eh/Bohr (Ende der SCF).
Iterationstoleranz der neuen ctests: die bestehende ±1 (`--tol-it` unverändert): auf H2O/NH3/Benzol sg1 sind die Zahlen
in allen 30 Fällen gleich, die ctests nutzen nur diese Moleküle und sg1.

**Zeiten** (CPU 8 Threads gegen RTX 5090, sg1, `-qm.scf_threshold 1e-8`, Defaults fc + VWN5-P86, je 3 Läufe, Median;
Tabellen in [QM_GPU_ROADMAP.md 4.3](../QM_GPU_ROADMAP.md#43-double-hybrids-on-the-device-oct-31-2026-wp12-stage-3)): Energie
Benzol dsd-pbep86-d3bj 1.53 / 0.76 s (SCF-Loop 1091 / 197 ms, V_xc pro Aufbau 34.7 / 12.0 ms), dsd-pbeb95-d3bj 1.93 / 1.04 s;
Koffein 12.65 / 3.98 s bzw. 16.40 / 6.71 s (Iterationen bei 1e-8 gleich). Energie + Gradient Benzol 5.83 / 1.74 s
(dsd-pbep86-d3bj), 5.80 / 1.71 s (b2plyp), 6.97 / 2.39 s (dsd-pbeb95-d3bj); Koffein 121.0 / 14.6, 120.1 / 14.4, 161.0 / 21.1 s.
`-opt` Benzol: dsd-pbep86-d3bj 69.5 / 15.7 s, b2plyp 80.5 / 18.0 s, Endenergien gleich auf 12 Stellen. **Befund Koffein:** der
CPU-Z-Vektor bricht nach 100 Iterationen mit max|R| 3.4e-10 (dsd-pbep86-d3bj) bzw. 7.6e-9 (dsd-pbeb95-d3bj) ab (Warnung,
akzeptiert unter 1e-6), auf der GPU (exakte J/K aus dem Device-Tensor) konvergiert er in 18 Iterationen; vermutete Ursache
das Screening-Rauschen der integral-direkten Host-J/K (246 Funktionen, Tensor 29 GB). Ein Teil
des Koffein-Faktors (8.3x) stammt daher; Benzol 3.4x. Ursache gefunden und behoben: Abschnitt „Z-Vektor mit
integral-direkten J/K" unten.

**Regression** (Binary d0f21497 = `release/` vor der Änderung kopiert, mit seinem Plugin): CPU `-gpu none`, md5 der
`-dump_gradient`-Dateien und Energiezeilen, H2O_distorted/NH3/Benzol, 1 und 4 Threads, Methoden wie Stufe 2 plus
Doppelhybrid-Gradienten (b2plyp RI/ohne fc/kanonisch, dsd-pbep86 (VWN5, PZ81), dsd-pbep86-d3bj, dsd-pbeb95-d3bj,
revdsd-pbep86, b2gp-plyp-d3bj kanonisch, dsd-blyp ohne Gitterantwort) und -Energien: 174/180 bitgleich. **Abweichung:** die
6 Fälle `dsd-pbeb95-d3bj`-Gradient (Energiezeilen gleich, 1 = 4 Threads) weichen im letzten Bit ab: max |Δg| 8.3e-16
(H2O), 8.4e-16 (NH3), 2.0e-14 Eh/Å (Benzol). Ursache nicht gefunden. Geprüft: Objektcode von `qm_xc.cpp`, `qm_grid.cpp`,
`qm_xc_response.cpp` gleich dem von HEAD; die vom Linker genommenen schwachen Jet-Kopien `b95`/`kernelJet` in
`qm_dh_grad.cpp.o` bytegleich; `sizeof(QMEngine)` unverändert; ein Nachbau von HEAD im Scratch-Verzeichnis reproduzierte
schon b2plyp nicht bitgleich (als Referenz unbrauchbar). Zwischenstand: mit dem Device-Code direkt in `qm_dh_grad.cpp`
dieselbe B95-Abweichung, mit `functionalCovered` außer Linie (eine Jet-Quelle) wichen alle Doppelhybrid-Gradienten im
letzten Bit ab (verworfen). GPU (`-gpu cuda -qm.gpu_scf on`), lda/pbe/b3lyp (beide VWN)/tpss/scan/pbe ohne Antwort ×
H2O_distorted/NH3/Benzol: ohne Device-1e/2e-Gradient (keine Atomics) 21/21 bitgleich; mit (Atomics) Energiezeilen
21/21 gleich, max |Δg| ≤ 1.0e-13 Eh/Bohr, das Referenz-Binary gegen sich selbst 7.9e-14 (Benzol scan).

## Stufe 4: D4 (CPU, 31. Okt 2026)

**Befund vorab (was es gab):** nativ vorhanden war das D4 der GFN2-Implementierung (`dispersion/`: Referenzdaten,
Gauß-Gewichte und Ladungs-Zeta pro Referenzzustand `D4ParameterGenerator::buildAtomRefW`/`contractC6Gfn2`,
Zweikörperterm `D4Evaluator::computeEnergyAndGradient`, ATM `computeATM`, kovalente EN-gewichtete CN `d4_ncoord`, EEQ mit
analytischem dq/dR `D4ChargeModel`). Es fehlte für DFT-D4: (a) die EEQ-Variante der Referenzladungen (dftd4 `clsq`/`clsh`;
curcuma hatte nur die GFN2-Variante `refq`/`refh`), (b) dftd4s Ladungsmodell EEQ(2019) (`D4ChargeModel` rechnete den
angewChem2020-Satz von GFN-FF mit der GFN-FF-CN, Kappung 4.4), (c) funktionalspezifische Parameter. Extern: `external/cpp-d4`
(C++-Port, nur mit `USE_D4`, in `release/` aus; nicht benutzt), kein dftd4-Binary; Referenz: dftd4 4.3.1 von PyPI in einem
venv im Scratchpad (`pip install dftd4==4.3.1 tomli`), Quellen dftd4 v4.3.1 und multicharge v0.6.0 (GitHub-Tarballs) für
die Daten; ORCA 6.1 (enthält dftd4 3.4.0 laut Ausgabe).

**Code:** `D4ParameterGenerator::setReferenceChargeVariant(Eeq)` lädt `clsq`/`clsh` (Default `Gfn2` = bisheriges Verhalten),
`D4ChargeModel::setModel(Eeq2019)` = multicharge `new_eeq2019_model` (chi, eta, kcnchi, rad; erf-CN mit kcn 7.5 über
D3-Kovalenzradien × 4/3, 25 Bohr, log-Kappung bei 8; adjungierter dq/dR-Gradient), Default `GfnffAngew2020` unverändert;
Daten generiert in `dispersion/d4_eeq_tables.h` (`scripts/d4_export_eeq_tables.py`). `curcuma::dispersion::D4DftDispersion`
(`dispersion/d4_dft.*`) ruft genau die GFN2-Bausteine auf: EEQ-Ladungen → `setTopologyCharges` → `GenerateParameters`
(kovalente CN an) → `computeEnergyAndGradient` (per-reference charge, dE/dCN, dE/dq) → `computeATM` (q = 0, 40 Bohr) →
`addD4CovalentCNGradient` (30 Bohr) → `addChargeResponseGradient`. Parameter: `d4DftParameters()` (eine Tabelle mit Quelle
pro Satz), Doppelhybride über die neuen Zeilen `b2plyp-d4`, `b2gp-plyp-d4` und `revdsd-pbep86-d4` (Parameter nur in der
D4-Tabelle) bzw. `-qm.dispersion d4` (für `b2plyp`, `b2gp-plyp`, `revdsd-pbep86`, deren `-d4`-Zeile dieselben XC/PT2-Zahlen
hat), KS über `pbe-d4`, `b3lyp-d4`, `tpss-d4`, `scan-d4` oder `-qm.dispersion d4` (`QMMethod`). Native GFN2 ruft
`D4DftDispersion` nicht auf; sein Code wurde nur um die zwei Schalter (Defaults = alt) erweitert.

**Parameter (D4(EEQ)-ATM, s9 1, alpha 16, a2 in Bohr):**

| Schlüssel | s6 | s8 | a1 | a2 | Quelle |
|---|---|---|---|---|---|
| pbe | 1.0 | 0.95948085 | 0.38574991 | 4.80688534 | dftd4 4.3.1 `parameters.toml` (Fit Caldeweyher 2019, DOI 10.1063/1.5090222); = cpp-d4; ORCA `PBE D4` gleich |
| b3lyp | 1.0 | 2.02929367 | 0.40868035 | 4.53807137 | wie oben; ORCA `B3LYP D4` gleich |
| tpss | 1.0 | 1.76596355 | 0.42822303 | 4.54257102 | wie oben; ORCA `TPSS D4` gleich |
| scan | 1.0 | 1.46126056 | 0.62930855 | 6.31284039 | wie oben; ORCA `SCANfunc D4` gleich |
| b2plyp | 0.64 | 1.16888646 | 0.44154604 | 4.73114642 | wie oben; ORCA `B2PLYP D4` gleich |
| b2gp-plyp | 0.56 | 0.94633372 | 0.42907301 | 5.18802602 | wie oben; ORCA `B2GP-PLYP D4` gleich |
| revdsd-pbep86 | 0.5917 | 0 | 0.371 | 4.201 | Mehta et al. 2023 Tab. 1, PZ81-Zeile (arXiv 2211.03105 gelesen; „s8 = 0 throughout; cATM = 1.0 for D4“) |

- „ORCA gleich“: ORCAs eingebauter Satz gibt dieselbe D4-Energie wie der explizit gesetzte (`%method D4S6 ... end`) auf allen 9
  gedruckten Stellen, 5 Moleküle. **revDSD-PBEP86-D4/2021 in ORCA weicht ab:** ein Fit an ORCAs Energien (5 Moleküle) ergibt
  s6 0.5917, a1 0.371, a2 **4.2014** (Residuen ≤ 5e-10 Eh, bei a2 4.2013/4.2015 bis 6e-7 Eh); Unterschied zu 4.201: H2O 6.6e-8,
  Benzol 2.4e-6 Eh. curcuma folgt dem Paper; die ORCA-Vergleiche setzen `D4A2 4.201`.
- **Nicht eingebaut:** D4 für `dsd-blyp`, `dsd-pbep86`, `dsd-pbeb95` (dftd4 hat Sätze mit s6 0.54/0.47/0.54, die zu keinem der
  DSD-XC/PT2-Sätze dieser Tabelle passen, KM2013 D3BJ-s6 0.57/0.48/0.61; welche DSD-Definition dftd4 meint, ist nicht geklärt),
  `hf`/`lda` (keine Anforderung), Mehtas VWN5-Satz von revDSD-PBEP86-D4 (0.4200/0.6090/0.0476, s6 0.5850, a1 0.370, a2 4.183).
- B3LYP: die dftd4-Parameter sind für „B3LYP“ gefittet; ob mit VWN5 oder VWN-RPA, ist nicht geprüft; `b3lyp-d4` nimmt sie für
  beide `-qm.vwn_variant`.

**Validierung** (`ctest -L qm_d4`, Referenzen `test_cases/qm_dh/reference_d4.json` = `scripts/qm_d4_dftd4_reference.py`,
dftd4 4.3.1 Python, Modell d4, Default-Cutoffs; `reference_orca_d4.json` = `scripts/qm_d4_orca_compare.py`):

| Größe (35 Fälle: H2O, NH3, CH4, H2CO, Benzol × 7 Sätze) | max \|Δ\| | Gate |
|---|---|---|
| E_D4 gegen dftd4 4.3.1 | 8.7e-18 Eh | 1e-12 |
| Zweikörper / ATM gegen dftd4 (`get_pairwise_dispersion`) | 1.6e-17 / 8.5e-21 Eh | 1e-12 |
| EEQ-Ladungen / CN gegen dftd4 | 3.4e-15 / 4.4e-16 | 1e-10 |
| Gradient gegen dftd4 | 8.7e-19 Eh/Bohr | 1e-11 |
| Gradient gegen 4-Punkt-FD der eigenen D4-Energie (h 5e-4 Bohr) | 3.9e-14 Eh/Bohr | 1e-10 |
| E_D4 gegen ORCA 6.1 (9 Dezimalen) | 4.9e-10 Eh | 1e-9 |
| D4-Gradient gegen ORCA 6.1 (`DISPERSION GRADIENT`, 9 Dezimalen) | 5.0e-10 Eh/Bohr | 1e-9 |

- Geladenes System (nicht als ctest): H3O+ (Geometrie von Hand), pbe-Satz, Ladung +1 gegen dftd4 `charge=1`: ΔE 2.7e-19 Eh,
  Ladungen 8.3e-16, Gradient 9.8e-20, gegen FD 5.3e-16 Eh/Bohr.
- **Methode mit D4 = Methode ohne + D4-Term** (`qm_dh_d4_*`, `qm_ks_d4_*`; def2-SVP sg1, VWN5-P86, fc): ΔE ≤ 6.1e-15 Eh, Gradient
  ≤ 2.0e-17 Eh/Bohr (H2O revdsd-pbep86-d4 und b2plyp-d4, NH3 b2gp-plyp-d4, H2CO revdsd-pbep86-d4; H2O pbe/b3lyp/scan, NH3
  tpss); Gesamtgradient gegen 4-Punkt-FD der Gesamtenergie (mitwanderndes Gitter, h 1e-3 Bohr): H2O revdsd-pbep86-d4 9.0e-11,
  b2plyp-d4 8.8e-11 Eh/Bohr.
- **revDSD-PBEP86-D4 gegen ORCA 6.1** (Einstellungen wie Stufe 1/2: def2-SVP, NoRI, AuxC def2/JK, RI-PT2, frozen core,
  DefGrid3 + IntAcc 8 AngularGrid 7 GridPruning 0, XMuePBE, VWN5-P86 auf beiden Seiten, `D4A2 4.201`; curcuma 250x1454):

| Molekül | ΔE gesamt | ΔE_D4 | max \|Δg\| (Eh/Bohr) |
|---|---|---|---|
| H2O | +2.3e-8 | +1.4e-10 | 1.5e-6 |
| NH3 | −2.3e-8 | −3.2e-10 | 6.4e-7 |
| H2CO | +4.8e-9 | −1.3e-10 | 5.1e-7 |

  Größenordnung wie revdsd-pbep86 ohne D4 in Stufe 1/2 (Gitter); der D4-Anteil liegt an ORCAs Druckgenauigkeit.
- **`-opt`** (`qm_dh_cli_opt_d4`): H2O_distorted mit revdsd-pbep86-d4, `-opt.convergence_preset tight`: Exit 0 nach 10 Schritten,
  letzte gedruckte Gradientennorm 7e-6. Gegen ORCA-`TightOpt` (revDSD-PBEP86-D4/2021 mit `D4A2 4.201`, Einstellungen wie oben, curcuma
  250x1454 + fc, `scripts/qm_dh_opt_orca_compare.py`-Funktionen): max |Δd_ij| H2O 5.3e-6 Å, H2CO 8.3e-6 Å, Energien gleich auf den 6
  gedruckten Stellen der curcuma-xyz-Datei.
- **GPU** (`-gpu cuda`, RTX 5090): D4 läuft auf der CPU und zählt keinen Fallback; revdsd-pbep86-d4 H2O Energie GPU = CPU auf 12
  Stellen, `-gpu_strict true` Energie Exit 0; Gradient GPU gegen CPU 2.0e-14 Eh/Å (einziger Fallback der bekannte
  `qm_dh_grad_pt2`); pbe-d4 Gradient GPU gegen CPU 8.3e-14 Eh/Å, `-gpu_strict true` mit Gradient Exit 0.
- **CLI:** `qm_dh_cli_d4` (revdsd-pbep86-d4 rechnet und druckt E(disp D4); `-qm.dispersion d4` bei `dsd-blyp` = Fehler mit
  Liste der D4-Namen), `qm_dh_cli_ks_d4` (pbe-d4 = pbe + E(disp D4) und = `pbe -qm.dispersion d4`).
- **Regression:** HEAD (76e5cbe3) und HEAD + nur diese Änderung, beide im selben Scratch-Build-Verzeichnis gebaut (Worktree, gleiche
  CMake-Optionen wie `release/`, `curcuma` ohne neu gebautes CUDA-Plugin gelinkt; die Arbeitskopie enthielt parallel Änderungen
  anderer Aufträge an `qm_dh_grad*` und `method_catalog*`, die hier nicht enthalten sind), CPU `-gpu none`, md5 der
  `-dump_gradient`-Datei + Energiezeilen (`-verbosity 2`), H2O_distorted/NH3/Benzol, 1 und 4 Threads: gfn1, gfn2, gfnff, hf, lda,
  pbe, b3lyp (+ RI-JK), tpss, scan, hf-3c, mp2 (Gradient), ri-mp2, ccsd, lmp2, b2plyp (Gradient; Energie ohne fc),
  dsd-pbep86-d3bj, revdsd-pbep86, dsd-pbeb95-d3bj (Gradienten): 58/58 (4 Threads) und 58/58 (1 Thread) bitgleich. **Befund
  dazu:** dasselbe HEAD im Scratch-Verzeichnis gebaut ist gegen das `release/`-Binary NICHT bitgleich (alle Gradienten
  ≤ 7e-13 Eh/Å, gfn1/gfn2 bis 1.4e-6 Eh/Å bei gleichen 8-stelligen Energien; jedes Binary für sich reproduzierbar), wie
  schon in Stufe 3 beobachtet; Ursache nicht untersucht (Kandidaten: im Worktree fehlendes `external/ulysses-main`,
  `-DCURCUMA_SOURCE_DIR`). Deshalb der Vergleich innerhalb eines Build-Verzeichnisses. ctest auf dem `release/`-Build mit
  dieser Änderung (und den parallelen Änderungen an `qm_dh_grad*`): `-L "qm_dh|qm_dh_grad|qm_mgga|qm_mgga_grad|qm_lda|qm_pbe|
  qm_b3lyp|qm_dft_grad|qm_gpu_dft|qm_grid|qm_hf3c|qm_mp2|qm_cite|qm_ri|qm_lmp2|qm_ccsd"` (-j6) 703/703, `-R
  "sqm_val|d4|xtb_gradient|gradient_unit_contract"` ohne die neuen qm_d4-Fälle 82/82, `-L qm_d4` 45/45 (dazu `qm_dh_cli_d4`). Energie eines
  Doppelhybrids mit D4 1 = 4 Threads bitgleich (Benzol revdsd-pbep86-d4).

**Befunde an den vorhandenen GFN2-Daten** (nicht geändert, GFN2 bleibt bitgleich; `d4_export_eeq_tables.py` druckt den
Vergleich mit dftd4 4.3.1 `reference.inc`): `d4_refq`/`d4_refh_charges` haben 10 statt 14 Dezimalen, Sc (Z = 21) fehlt die
4. Referenzladung (refq 1.0291, refh −0.3430), `d4_refcovcn` hat für Sr/Ba einen 4. Eintrag (15.121/14.677, dftd4 0). Auswirkung
auf GFN2 nicht gemessen. Die DFT-D4-Paarliste endet bei 50 Bohr (`PAIR_EVAL_CUTOFF_BOHR`), dftd4 bei 60: Paare zwischen 50 und
60 Bohr fehlen (keine Testmoleküle dieser Größe).

**Nicht getestet:** Elemente außer H, C, N, O (die DFT-Engine kann H-Ne; F, Ne, Li-B nicht verglichen), Moleküle größer als
Benzol, offene Schale, Ladung nur H3O+ (pbe-Satz, ohne SCF), `-md`, GPU außer RTX 5090, Ferngrenzen der Cutoffs (CN 30 Bohr,
ATM 40 Bohr, Paarliste 50 Bohr), Threads > 1 im D4-Teil (läuft seriell).

### Z-Vektor mit integral-direkten J/K (Befund und Behebung, 31. Okt 2026)

**Ursache:** `qmint::DirectJK::build` überspringt ein Quartett, wenn Q_AB Q_CD max|X| (über die sechs berührten
Schalenblöcke) unter `eri_screening` (1e-12) liegt, eine ABSOLUTE Schranke. Für die SCF-Dichte ist das gewollt, für den
linearen Operator G4(X) des Z-Vektors nicht: J/K(aX) ≠ a J/K(X). Die CG-Suchrichtungen schrumpfen mit dem Residuum
(max|p| ~ 1e-9 ab Iteration 8-15), dann fällt der größte Teil der Quartette weg, der Operator wird unsymmetrisch und
inhomogen, die CG stagniert bei max|R| ~3e-10 (bzw. 5e-9 bei dsd-pbeb95-d3bj). Koffein läuft integral-direkt, weil der
Tensor über `eri_max_memory_mb` liegt; Benzol, H2O, NH3 laufen sonst mit gespeichertem Tensor (Schwarz-Screening
dichteunabhängig, linear). Die GPU nutzt den Device-Tensor, darum dort 18 Iterationen.

**Messung** (Diagnose-Binary = HEAD 76e5cbe3 + temporärer Prüfblock `CURCUMA_DH_ZVEC_CHECK`, nicht eingecheckt; X =
erster Z-Vektor-Schätzer rhs/diag, dsd-pbep86-d3bj, sg1, def2-SVP, 8 Threads): relativer Fehler max|J(aX)/a − J(X)| / max|J(X)|.

| a | Benzol direkt J / K | Koffein direkt J / K | F_xc (beide) | Benzol Tensor J | Koffein GPU-Tensor J |
|---|---|---|---|---|---|
| 1e-2 | 2.5e-7 / 5.9e-8 | 9.5e-8 / 3.6e-8 | ≤ 2.2e-15 | 1.3e-14 | 3.8e-15 |
| 1e-4 | 1.7e-5 / 1.1e-5 | 1.0e-5 / 3.4e-6 | ≤ 2.2e-15 | 1.4e-14 | 3.4e-15 |
| 1e-6 | 1.9e-3 / 1.3e-3 | 8.6e-4 / 3.2e-4 | ≤ 2.2e-15 | 1.1e-14 | 3.6e-15 |
| 1e-8 | 1.4e-1 / 9.3e-2 | 6.6e-2 / 3.1e-2 | ≤ 2.2e-15 | 1.1e-14 | 3.8e-15 |
| 1e-10 | 4.6 / 0.50 | 1.9 / 0.48 | ≤ 2.2e-15 | 1.3e-14 | 3.4e-15 |

Ganzer Operator M = diag + G4 (applyM): Symmetrie |<y|Mx> − <x|My>| / |<y|Mx>| und Homogenität bei a = 1, 1e-6, 1e-9:
Benzol direkt 4.2e-9, 9.5e-3, 5.4 (Homogenität 3.0e-2, 20); Koffein 6.7e-9, 1.7e-2, 9.6 (1.0e-2, 8.8); Benzol Tensor und
Koffein GPU ≤ 1.3e-13. Bei a = 1 gegen DirectJK ohne Screening (thr 0): Benzol dJ 1.8e-11, dK 7.4e-12; Koffein 3.4e-11 /
1.6e-11; GPU-Tensor 1.5e-13 / 1.0e-13. Wiederholbarkeit: dieselbe direkte Anwendung zweimal ist mit 1 Thread bitgleich, mit
4/8 Threads nicht (dJ 2.6e-16 .. 5.9e-16; Thread-Teilsummen in `omp critical` in zufälliger Reihenfolge, TODO.md);
Tensorpfad bitgleich, 1 = 8 Threads. **Nicht die Ursache:** F_xc (linear bis 2e-15, inkl. `xc_ao_threshold`),
Preconditioner, Abbruchkriterium und CG-Rundung (dieselbe CG mit exakten J/K läuft auf der GPU und nach der Behebung auf
der CPU mit identischer Residuenfolge in 18 Iterationen), inkrementeller Fock-Aufbau (die Z-Vektor-Aufbauten sind volle
Aufbauten, `m_direct_thr_now` = −1). `-qm.eri_prim_density false` hilft nicht (Benzol direkt 100 Iterationen);
`-qm.eri_screening 1e-15` oder 0 behebt es (Benzol 11 Iterationen mit beiden, Koffein 18 mit 1e-15), verlangsamt aber auch die SCF.

Iterationen (HEAD, CPU, 8 Threads, dsd-pbep86-d3bj): H2O direkt 9 (Tensor 9), NH3 direkt 35 (Tensor 10), Benzol direkt
100 (Tensor 11), Koffein 100; mit `-qm.eri_screening 1e-10` auch H2O 100, NH3 100. Der Effekt wächst mit dem Anteil
kleiner Schwarz-Produkte, also mit der Systemgröße und der Schwelle.

**Behebung** (nur `QMEngine::dhG4`, `qm_dh_grad.cpp`): vor `buildJKExact` wird X mit einer Zweierpotenz s auf max|X| in
[0.5, 1) skaliert, J und K danach mit 1/s. Das ist in Gleitkomma exakt (gespeicherter Tensor: dieselben Bits) und macht
die direkte Schranke relativ zu max|X| (für Host-DirectJK und Device Q-G3 gleich). Nachher Benzol direkt Symmetrie
1.3e-10 / 5.5e-11 / 7.0e-11, Homogenität 2.1e-11 / 2.6e-11; Koffein 2.5e-10 / 3.3e-10 / 2.2e-10, Homogenität 8.3e-11 /
7.3e-11 (Restfehler = relative Schranke 1e-12 auf Quartettebene).

**Koffein danach** (eigene Builds im Worktree: HEAD gegen HEAD + Fix, gleiche CMake-Optionen; CPU 8 Threads, Defaults fc +
VWN5-P86, sg1; die Maschine war durch einen parallelen Build/Lauf belastet (Last ~35), Zeiten nur paarweise vergleichbar;
Gradienten in Eh/Bohr, max über alle Komponenten):

| Funktional | Iter. HEAD / Fix | max\|R\| HEAD / Fix | Z-Vektor-Zeit HEAD / Fix | gesamt HEAD / Fix | HEAD − GPU | Fix − GPU | HEAD − Fix |
|---|---|---|---|---|---|---|---|
| dsd-pbep86-d3bj | 100 / 18 | 3.3e-10 / 6.6e-11 | 117.6 / 25.9 s | 276.8 / 174.0 s | 2.4e-10 | 9.2e-11 | 1.5e-10 |
| dsd-pbeb95-d3bj | 100 / 18 | 4.6e-9 / 3.9e-11 | 84.8 / 23.2 s | 189.5 / 126.7 s | 1.4e-9 | 8.3e-11 | 1.4e-9 |
| b2plyp | 100 / 17 | 3.1e-10 / 5.9e-11 | 41.2 / 15.4 s | 141.0 / 115.0 s | 1.2e-10 | 5.5e-11 | 1.2e-10 |

Die Residuenfolge des Fix ist die der GPU (18/18/17 Iterationen, gleiche letzte Werte). Mit `-qm.eri_screening 1e-15`
(SCF und Z-Vektor): HEAD gegen Fix 1.5e-13 (dsd-pbep86-d3bj) bzw. 1.4e-13 (dsd-pbeb95-d3bj), GPU gegen 1e-15 je 2.2e-11;
Fix mit Default-Screening gegen Fix mit 1e-15: 7.6e-11 bzw. 6.7e-11. Der verbleibende Abstand Fix − GPU (≤ 9.2e-11) kommt
also aus der Default-Schwelle der direkten KS-SCF, nicht aus dem Z-Vektor. GPU-Läufe HEAD gegen Fix: ≤ 9.7e-14 Eh/Å
(Device-Atomics; HEAD gegen sich selbst 1.4e-13). **Hinweis:** dsd-pbeb95-d3bj ändert sich um 1.4e-9 Eh/Bohr, mehr als
die vorher gesetzte Grenze 1e-9; das ist der Fehler des alten, nicht konvergierten Z-Vektors (gegen GPU und 1e-15).

**Regression** (Worktree-Builds HEAD gegen HEAD + Fix, CPU `-gpu none`, md5 der `-dump_gradient`-Dateien und der
Energiezeilen, H2O_distorted/NH3 mit 1 und 4 Threads, Benzol 4 Threads; hf, pbe, b3lyp, tpss, scan, hf-3c, lda, mp2
(nicht Benzol), b2plyp, dsd-pbep86, dsd-pbep86-d3bj, dsd-pbeb95-d3bj, revdsd-pbep86, b2gp-plyp-d3bj, b2plyp kanonisch):
73/73 bitgleich (die DH-Läufe dieser Moleküle nutzen den gespeicherten Tensor; die Zweierpotenz-Skalierung ändert dort kein Bit).
Kleine Moleküle direkt erzwungen (`-qm.scf_direct on`): Benzol HEAD 100 Iterationen, Fix 11, Gradient gegen direkt ohne
Screening 7.1e-11 (HEAD) bzw. 4.4e-12 (Fix); NH3 35 → 10, 3.2e-11 → 3.0e-14 Eh/Bohr. ctest `qm_dh_grad_direct_*` (H2O
dsd-pbep86-d3bj VWN5, NH3 b2plyp; erzwungen direkt mit `eri_screening` 1e-10 gegen den Tensor; Grenzen: Residuum ≤ 1e-10,
höchstens 2x + 2 Iterationen des Tensorlaufs, Gradient ≤ 1e-10): HEAD 100 Iterationen, 5.0e-9 / 9.2e-8, ΔG 4.0e-10 /
2.9e-9 (FAIL); Fix 9 Iterationen, 1.2e-11 / 1.9e-11, 2.8e-14 / 1.8e-12 (OK).

Nicht geprüft: MP2-Gradient (`mp2G4`, nutzt nach direkter SCF den Tensor, unverändert), `-gpu` mit direkten Device-J/K
(Q-G3, `-qm.gpu_direct on`) im Z-Vektor, andere Basen, Moleküle größer als Koffein, `-opt` mit direkter SCF.

## Nicht getestet / nicht implementiert

- Nicht implementiert: PT2-Gradient mit RI-J/RI-JK-SCF (`-qm.ri j|jk`: Fehler), PT2-Zwischengrößen des Gradienten auf dem
  Device (gezählter Fallback `qm_dh_grad_pt2`, Stufe 3), ROCm/Vulkan, Hessian, D4 für die DSD-Zeilen (Stufe 4), VV10/NL, offene Schale, Z > 10, `dsd-pbeb95` ohne
  Dispersion, `revdsd-pbep86-d3bj`, PWPB95, ωB2PLYP u. a. range-separierte DH, SCS-Varianten mit eigenen Einfachanregungen.
- Stufe 3 nicht getestet: andere GPUs als die RTX 5090, Moleküle über Koffein, `-md`, kanonischer PT2-Gradient über Benzol
  auf dem Device (n⁴-Γ wird hochgeladen), Stratmann-Partition, kartesische d, `-qm.xc_ao_threshold 0` mit den neuen
  Kernen, Multi-GPU; `-opt` nur Benzol (dsd-pbep86-d3bj, b2plyp) und ctest H2O dsd-pbep86-d3bj.
- Nicht getestet: andere Basen als def2-SVP, def2-SVP/C (nicht vorhanden), `-qm.ri j|jk` im SCF der DH (läuft über den
  B3LYP-Pfad, nicht gemessen), integral-direct SCF + kanonisches PT2 großer Systeme, Moleküle größer als Benzol, geladene Systeme;
  Gradient: `-md`, Moleküle größer als Benzol (Gradient nur bis H2CO validiert, Benzol nur Threads), integral-direct SCF,
  kartesische d, sg2/sg3-Gitter, Stratmann-Partition, `-qm.xc_ao_threshold 0` (die Tests laufen mit dem Default-Screening),
  Speicherbedarf des RI-Gradienten (n²·naux für B aller MO-Paare) bei großen Systemen.
- Stufe 2: die Jet-Kerne in `qm_xc_fxc_kernels.h` (bis Stufe 2 `qm_xc_fxc.h`) sind eine zweite Niederschrift der Kernformeln (nur für f_xc); eine Änderung an einem
  Produktionskern muss dort nachgezogen werden (`qm_dh_grad_kernel` vergleicht die ersten Ableitungen).
- LiH mit dsd-pbeb95-d3bj auf sg1: SCF konvergiert nicht (curcuma und PySCF stagnieren bei |g| ≈ 6e-5, E_KS beider −8.0158776);
  auf 99x590 in 10 Iterationen konvergiert. Kein ctest-Fall, Ursache nicht untersucht.
- Die Parameter der DSD-Zeilen stützen sich auf Programme (ORCA, Psi4), nicht auf das gelesene Paper.
- Menschlicher Produktionstest ausstehend.

## Quellen

S. Grimme, J. Chem. Phys. 124, 034108 (2006), DOI 10.1063/1.2148954; A. Karton et al., J. Phys. Chem. A 112, 12868 (2008),
DOI 10.1021/jp801805p; S. Kozuch, D. Gruzman, J. M. L. Martin, J. Phys. Chem. C 114, 20801 (2010), DOI 10.1021/jp1070852;
S. Kozuch, J. M. L. Martin, Phys. Chem. Chem. Phys. 13, 20104 (2011), DOI 10.1039/c1cp22592h; dies., J. Comput. Chem. 34, 2327
(2013), DOI 10.1002/jcc.23391; G. Santra, N. Sylvetsky, J. M. L. Martin, J. Phys. Chem. A 123, 5129 (2019),
DOI 10.1021/acs.jpca.9b03157; N. Mehta, G. Santra, J. M. L. Martin, Can. J. Chem. 101, 656 (2023), DOI 10.1139/cjc-2022-0277;
L. Goerigk, S. Grimme, Phys. Chem. Chem. Phys. 13, 6670 (2011), DOI 10.1039/c0cp02984j; J. P. Perdew, Phys. Rev. B 33, 8822
(1986), DOI 10.1103/PhysRevB.33.8822; J. P. Perdew, A. Zunger, Phys. Rev. B 23, 5048 (1981), DOI 10.1103/PhysRevB.23.5048;
A. D. Becke, J. Chem. Phys. 104, 1040 (1996), DOI 10.1063/1.470829 (alle DOIs gegen Crossref geprüft, 31. Okt 2026);
libxc 7.0.0 `gga_c_p86.mpl`, `gga_c_p86vwn.mpl`, `lda_c_pz.mpl`, `mgga_c_bc95.mpl` (Formeln, Konstanten); Psi4
`dh_functionals.py`, simple-dftd3 und dftd4 `parameters.toml` (Stand Okt 2026).
Stufe 4: E. Caldeweyher, C. Bannwarth, S. Grimme, J. Chem. Phys. 147, 034112 (2017), DOI 10.1063/1.4993215; E. Caldeweyher
et al., J. Chem. Phys. 150, 154122 (2019), DOI 10.1063/1.5090222; E. Caldeweyher et al., Phys. Chem. Chem. Phys. 22, 8499 (2020),
DOI 10.1039/D0CP00502A (alle drei und Mehta et al. 2023, Santra et al. 2019 gegen Crossref geprüft, 31. Okt 2026); dftd4 v4.3.1
(`src/dftd4/reference.inc`, `model/d4.f90`, `disp.f90`, `cutoff.f90`, `assets/parameters.toml`) und multicharge v0.6.0
(`param.f90`, `param/eeq2019.f90`, `model/eeq.f90`), GitHub-Tarballs vom 31. Okt 2026.
