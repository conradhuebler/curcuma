# WP12 — Doppelhybride (B2PLYP, B2GP-PLYP, DSD-*, revDSD-PBEP86)

- **Status Stufe 1 (Energien, CPU):** 🤖 AI-generated, ⚙️ machine-tested (`ctest -L qm_dh`, 31. Okt 2026); menschlicher
  Produktionstest ausstehend. **Stufe 2 (Gradient)** und **Stufe 3 (GPU-Kerne, D4)** offen (TODO.md).
- Name: „WP12“ stand in der README-Tabelle für die gestrichene Bewertung der Schwarz-Ansätze; diese Datei ist das neue WP12.
- **Code:** `qm_methods/double_hybrid_method.{h,cpp}` (Tabelle `doubleHybridTable()`, Klasse `DoubleHybridMethod`, Parameter
  `-qm.dispersion`, `-qm.dh_pt2`, `-qm.p86_lda`), Kerne `P86_C`, `P86VWN_C`, `B95_C` in `qm_xc_kernels.h`,
  `QMFunctional::DoubleHybrid` + `QMEngine::setXCFunctional` (`qm_engine.*`), PT2 = `computeRIMP2`/`computeMP2` (nur die
  HF-Abfrage erweitert), D3(BJ) = `D3ParameterGenerator` (der von hf-3c), Zitate `qmcite::doubleHybrid`.
  Tests `test_cases/qm_dh/`, Referenzen `scripts/qm_dh_pyscf_reference.py` (PySCF 2.14, libxc 7.0.0),
  `scripts/qm_dh_orca_compare.py` (ORCA 6.1.0).

## Was gebaut ist

E = E_KS[a_x HF + Σ w_k E_xc,k] + c_os E_os + c_ss E_ss + E_disp: geschlossenschalige Hybrid-KS-SCF (Fock wie WP7), dann die
MP2-Paarsummen von `qm_mp2.h` mit den KS-Orbitalen und -Eigenwerten (keine Einfachanregungen, wie in allen DH-Definitionen),
getrennt nach Spin (os/ss), skaliert. PT2 standardmäßig RI (`-qm.ri_mp2_aux`, Default def2-universal-JKFIT = def2/JK; kein
def2-SVP/C vorhanden: Methodenvariante, Ausgabe sagt es), `-qm.dh_pt2 canonical` = exakte Integrale. Frozen core über
`-qm.frozen_core` (Default false wie mp2; ORCA rechnet DH mit frozen core). Nur Energien: `hasGradient()` false, Gradient
und `-opt` enden mit Fehler; `-gpu cuda`: J/K und RI-MP2-Summe auf dem Device (vorhandene Pfade), XC mit P86/B95 als gezählter
`qm_xc`-Fallback auf der CPU, `-gpu_strict true` endet mit Exit 3.

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
| `revdsd-pbep86-d4` | wie oben | | | | | D4 s6 0.5917, a1 0.371, a2 4.201 | Mehta et al. 2023; **D4 nicht implementiert** (Fehler) |

- Kozuch/Martin 2013 (JCC 34, 2327) und Santra et al. 2019 (JPCA 123, 5129) sind nicht gelesen (nicht zugänglich); die
  DSD-Zahlen stammen aus ORCA 6.1 (gedruckte ScalHFX/ScalDFX/ScalDFC/„Scaling for aa/bb, ab pairs“) und Psi4
  `dh_functionals.py`, die revDSD-Zahlen aus Mehta et al. (arXiv 2211.03105, Can. J. Chem. 101, 656 (2023)).
- **Weggelassen:** `dsd-pbeb95` ohne Dispersion (ORCA 0.70 HF / 0.30 PBE / 0.52 B95, Psi4 0.69 / 0.31 / 0.54, beide 0.48 / 0.22,
  beide zitieren KM2013); `revdsd-pbep86-d3bj` (s6 0.4339 bei Mehta et al., 0.4377 in s-dftd3; c_os/c_ss nur eine Quelle).
- Widersprüche gefunden: die D4-Parameter von revDSD-PBEP86 in der dftd4-Tabelle (s6 0.5132, a1 0.44, a2 3.60) weichen von
  Mehta et al. ab; die s-dftd3-Einträge `dsdblyp` (s6 0.50, s8 0.2130, a2 6.0519) gehören zum DSD-BLYP von 2010, nicht zu KM2013.
- **P86-Lokalteil:** libxc `GGA_C_P86` (Perdew 1986, Gaussian, Q-Chem, Psi4) setzt den Gradiententerm auf PZ81, ORCA (kein
  PZ81) auf VWN5 (Ausgabe `LDAOpt VWN-5`). Default `-qm.p86_lda pz81`; `vwn5` reproduziert ORCA. Unterschied dsd-pbep86/def2-SVP:
  2.0e-3 Eh (H2O) bis 5.8e-3 Eh (Benzol). Die revDSD-Parameter wurden mit PZ81 gefittet (Mehta et al.); ORCAs
  `revDSD-PBEP86/2021` nutzt dieselben Zahlen auf VWN5 (Mehta et al. geben für VWN5 eigene: 0.4200 / 0.6090 / 0.0476).
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
auf `b2plyp-d3bj`; `-gpu cuda` (RTX 5090): gezählter `qm_xc`-Fallback, Energie = CPU (Δ 0 auf 12 Stellen), `-gpu_strict true`
Exit 3. `ctest -L qm_dh` 78/78; `ctest -L "qm_dh|qm_mgga|qm_mgga_grad|qm_lda|qm_pbe|qm_b3lyp|qm_dft_grad|qm_gpu_dft|qm_grid|qm_hf3c|qm_mp2|qm_lmp2|qm_ccsd|qm_ri|qm_cite"`
(-j4) 573/573 (372 s); nach der letzten Umstellung (Tabelle ohne std::string, Zitierfunktion außer Linie) `qm_dh|qm_cite` 86/86.
Laufzeit Benzol b2plyp sg1, 8 Threads: 1.6 s Wand (10 SCF-Iterationen, V_xc 0.53 s, RI-PT2 61 ms inkl. eigenem RI-Setup).

**(4) Regression:** gegen das Binary von 602f727f (vor der Änderung kopiert, gleiche CMake-Optionen), CPU (`-gpu none`),
`-dump_gradient`-md5 + Energiezeilen, H2O_distorted/NH3/Benzol, 1 und 4 Threads: hf, lda (beide VWN), pbe (+ ohne Gitterantwort,
ohne Screening), b3lyp (beide VWN, RI-JK), tpss, scan, mp2/scs-mp2 (+ fc; Gradient), ri-mp2 (+ fc), sos-ri-mp2, ccsd, lmp2, hf-3c:
112/112 bitgleich (erster Lauf: b3lyp-Gradienten 18/18 im letzten Bit anders, Ursache GCC-Inlining des neuen P86-Falls in
`evaluateKernel`; mit `noinline` behoben; nach der letzten Codeänderung erneut 112/112).

## Nicht getestet / nicht implementiert

- Nicht implementiert: Gradient (Stufe 2), P86/B95 auf dem Device (Stufe 3), D4, VV10/NL, offene Schale, Z > 10, `dsd-pbeb95` ohne
  Dispersion, `revdsd-pbep86-d3bj`, PWPB95, ωB2PLYP u. a. range-separierte DH, SCS-Varianten mit eigenen Einfachanregungen.
- Nicht getestet: andere Basen als def2-SVP, def2-SVP/C (nicht vorhanden), `-qm.ri j|jk` im SCF der DH (läuft über den
  B3LYP-Pfad, nicht gemessen), integral-direct SCF + kanonisches PT2 großer Systeme, Moleküle größer als Benzol, geladene Systeme,
  `-md`/Hessian (kein Gradient).
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
