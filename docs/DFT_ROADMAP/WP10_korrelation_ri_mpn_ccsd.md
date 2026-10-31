# WP10-WP14 -- Dichtefitting (RI), MPn und CCSD, plus Bewertung von Schwarz-Ansaetzen

- **Status:** WP10 RI-J/RI-JK (HF, CPU + CUDA), WP11 Schritt 1+3 (kanonisches RMP2, SCS/SOS, CPU + CUDA) und
  WP11 Schritt 2 (RI-MP2 mit def2-universal-JKFIT/-JFIT als Korrelations-Hilfsbasis, CPU + CUDA) sind
  implementiert, 🤖 AI-generated, ⚙️ machine-tested (Okt 2026, Abschnitte WP10, WP11 und RI-MP2 unten), ebenso der
  analytische Gradient fuer RI-J/RI-JK-HF (Abschnitt "WP10 -- Gradient") und der analytische Gradient des kanonischen
  MP2/SCS/SOS (CPU, 2e-Ableitungsteil auch CUDA; Abschnitt "WP11 -- MP2-Gradient"); RI-MP2 mit def2-SVP/C und der
  RI-MP2-Gradient sind Plan; WP13 (kanonisches RCCSD, CPU) ist implementiert (Abschnitt WP13); fuer WP15 (und die Metrikfrage aus WP14) liegt eine Python-Machbarkeitsmessung
  vor, und das lokale PNO-RI-MP2 ist als Methode `lmp2` implementiert (Abschnitt WP15, nur Energie, CPU).
- **Abhaengigkeit:** WP3 (RHF, ERI-Engine, `QMEngine`), fuer die GPU-Teile die Stufen Q-G1/Q-G2 aus
  [../QM_GPU_ROADMAP.md](../QM_GPU_ROADMAP.md).
- **Referenzen auf dieser Maschine:** PySCF (venv `~/.venvs/qmref`), ORCA 6.1 (`~/Downloads/orca_6_1_0_...`).
  Fuer Energien und Gradienten kann eine unabhaengige zweite Pruefung gegen ORCA sinnvoll bleiben; fuer Zeitmessungen nicht (Betreibervorgabe 2026-10-26: andere CPU, nicht vergleichbar).
- Reihenfolge und Abhaengigkeiten: WP11 (kanonisches MP2, O(N^5)) zuerst, dann WP10 (RI) und RI-MP2, danach
  WP13 (RCCSD); WP12 (MP3-MP5) ist **gestrichen** (Betreiber, Okt 2026: gemeint war MP2). WP14 ist die Auswertung
  der Schwarz-Frage.
- **Laufzeitvergleich mit ORCA 6.1 und PySCF 2.14** (HF, RI-JK, RI-MP2, LMP2/DLPNO-MP2, CCSD; 1 und 24 Threads,
  GPU; Okt 2026, 🤖 AI-generated, ⚙️ machine-tested, eine Maschine): [../QM_TIMING_ORCA_PYSCF_2026-10.md](../QM_TIMING_ORCA_PYSCF_2026-10.md);
  Nachmessung mit den aktuellen Voreinstellungen (25. Okt 2026, ORCA/PySCF am selben Tag neu) im Abschnitt "Re-measurement after the optimisation work".

## Warum diese Reihenfolge

Alles hier braucht die 4-Index-Integrale in der MO-Basis. Die native Engine hat den AO-Tensor
(gespeichert oder direkt) und einen RHF-Referenzzustand. Der Engpass wechselt dabei von der
Integralberechnung zu Matrix-Multiplikationen (n^5 fuer die Transformation, n^5 bis n^8 fuer die
Korrelation). Das ist GEMM-Arbeit und passt auf die GPU-Stufen; es ist aber **neu zu messen**, nicht aus den
xTB-Erfahrungen zu uebernehmen.

## WP10 -- Dichtefitting (RI-J, RI-JK; RI-C offen)

**Stand Okt 2026: RI-J und RI-JK fuer closed-shell HF, H-Ne, CPU und CUDA. 🤖 AI-generated, ⚙️ machine-tested;
menschlicher Produktionstest steht aus.** Code: `src/core/energy_calculators/qm_methods/qm_ri.{h,cpp}` (Integrale,
Metrik, Fit, Parser, Engine-Anbindung), `cuda/qm_gpu_ri.cu` + `qm_gpu_ri_backend.h` (GPU), Tests `test_cases/qm_ri/`.

- **Schalter** (Parameter-Registry, Modul `qm`): `-qm.ri off|j|jk` (Standard off; der exakte Pfad bleibt die Referenz
  aller bestehenden Tests), `-qm.ri_aux auto|def2/J|def2/JK|<Datei>` (auto: j -> def2-universal-JFIT, jk ->
  def2-universal-JKFIT), `-qm.ri_metric cholesky|eigen`, `-qm.ri_metric_cutoff` (1e-10), `-qm.ri_max_memory_mb`
  (16000), `-qm.gpu_ri auto|on|off`, `-qm.gpu_ri_fit off|on` (Metrik-Loesung auf dem Device, Abschnitt "Skalierung RI-JK");
  `-qm.ri_storage`, `-qm.gpu_ri_stream`, `-qm.ri_mixed_precision`, `-qm.ri_fit_radius` (Abschnitt "RI-JK ressourcensparend").
- **Verfahren:** RI-J: `c = J^-1 g`, `g_P = (P|ls) D_ls`, `J_mn = (mn|P) c_P`; Austausch exakt (gespeicherter Tensor
  oder integral-direkt wie ohne RI, entspricht ORCA `RIJONX`). RI-JK (Weigend 2002): `B = L^-1 (P|mn)` mit
  `J = L L^T`, `K_mn = sum_Q (B^Q D B^Q)_mn` ueber `D = V S V^T` (Eigenvektoren von D, mit sqrt|lambda| skaliert),
  einmal gespeichert (naux x n(n+1)/2, gepacktes unteres Dreieck je Zeile, Abschnitt "Gepackter 3-Index-Tensor"). Integrale: McMurchie-Davidson ueber Schalenpaare; eine Hilfsfunktion ist das Paar
  aus ihr und einer s-Funktion mit Exponent 0 am selben Zentrum. Hilfsfunktionen kugelfoermig (Nullraum des
  Laplace-Operators je L, gleiche Spanne wie die reellen Kugelflaechenfunktionen), mit `-qm.cartesian_d true`
  kartesisch wie die Orbitalbasis. **Max L der Hilfsbasen fuer H-Ne: 4 (g)** (JFIT Li-Ne, JKFIT B-Ne); mit der
  d-Orbitalbasis braucht `(mn|P)` L bis 8. Der vorhandene ERI-Code (GTO::Orbital, s-d) wird nicht benutzt.
- **Metrik:** Cholesky ist Standard (stabil, kein Abschneiden, gleiche Wahl wie PySCF); faellt auf die
  Eigenzerlegung zurueck, wenn J numerisch nicht positiv definit ist. Gemessene Eigenwerte der
  Coulomb-normierten Metrik ((P|P) = 1): H2O def2/J 2.0e-4..8.8 (cond 4.4e4), H2O def2/JK 1.0e-6..10.9 (1.1e7),
  Benzol def2/JK 2.7e-8..32 (1.2e9), Anthracen def2/JK 1.0e-8..49 (4.7e9), kartesische JK-Hilfsbasis H2O 1.5e-10..32
  (2.1e11). Cholesky und `eigen` (Cutoff 1e-10, nichts verworfen) geben fuer Benzol RI-JK dasselbe E_elec auf 10
  Stellen; `eigen` mit Cutoff 1e-6 verwirft 11 von 558 Funktionen und aendert E um -3.4e-6 Eh. Setup Benzol: Cholesky
  24 ms, `eigen` 190 ms.
- **Hilfsbasen (Herkunft, Lizenz):** `qm_methods/def2-universal-JFIT.turbomole` und `def2-universal-JKFIT.turbomole`,
  H-Ne-Bloecke unveraendert aus dem Basis Set Exchange (Version 0.12, Turbomole-Format, Kopf "Data from Turbomole 7.3",
  Kopf in den Dateien erhalten; byte-gleich mit dem Download geprueft). Der BSE-Kopf nennt keine Lizenz; BSE-Daten
  werden ueblicherweise mit Zitat der Originalarbeit genutzt (aus dem Gedaechtnis, nicht nachgeschlagen: Weigend,
  PCCP 8, 1057 (2006) fuer def2/J; Weigend, PCCP 4, 4285 (2002) fuer def2/JK) -- **Lizenzfrage fuer den Betreiber
  offen**. Der lokale JKFIT-Download bricht
  nach Pd ab (H-Ne vollstaendig). Dieselben Daten wie ORCA `def2/J`/`def2/JK` (Energien unten auf 2e-9 Eh).
  **def2-SVP/C (RI-C) fehlt**: RI-MP2 benutzt deshalb def2-universal-JKFIT (Standard) oder -JFIT als
  Korrelations-Hilfsbasis (Betreiberentscheidung Okt 2026; Abschnitt RI-MP2 unten).
- **Nicht implementiert:** RI-C mit def2-SVP/C, hf-3c mit RI (ignoriert mit Warnung: HF-3c ist mit exakten Integralen in MINIX
  definiert, `qm_ri_hf3c_ignored`), offene Schalen, Elemente ab Na, device-residenter SCF-Loop mit RI (RI laeuft im Host-Loop). Integral-direktes RI
  (batched) und lokaler Austausch-Fit seit Okt 2026, Abschnitt "RI-JK ressourcensparend".

**Validierung (alle Zahlen gemessen, HF/def2-SVP, `-qm.scf_threshold 1e-10`):**

- **Integrale vs PySCF** (H2O, kartesische Orbital- und Hilfsbasis, `-qm.eri_screening 0`, Dump ueber
  `CURCUMA_QM_RI_DUMP=<prefix>`): (P|Q) Coulomb-normiert max. Abweichung 2.7e-11 (def2/J) / 2.9e-11 (def2/JK),
  (mn|P) normiert 2.5e-11 / 2.6e-11, gefitteter Tensor sum_PQ (mn|P) J^-1 (Q|ls) 3.5e-11 (enthaelt f und g).
- **Gesamtenergien** (Eh; PySCF mit denselben Hilfsbasisdateien, ORCA 6.1 `RIJONX def2/J` bzw. `RIJK def2/JK`,
  `VeryTightSCF`):

  | Molekuel | Modus | curcuma | -PySCF | -ORCA | RI-Fehler vs exakt (Eh, Eh/Atom) |
  |---|---|---|---|---|---|
  | H2O (verzerrt) | RI-J | -75.9535706759 | +1e-10 | +7.1e-10 | -1.26e-4, -4.2e-5 |
  | H2O | RI-JK | -75.9533909924 | 0 | +6.3e-10 | +5.34e-5, +1.8e-5 |
  | NH3 | RI-J | -56.1402393646 | 0 | -1.3e-9 | -1.42e-4, -3.5e-5 |
  | NH3 | RI-JK | -56.1400758658 | 0 | -1.3e-9 | +2.17e-5, +5.4e-6 |
  | Formaldehyd (verzerrt) | RI-J | -113.7704540851 | +1e-10 | -1.6e-9 | -1.24e-4, -3.1e-5 |
  | Formaldehyd | RI-JK | -113.7702284010 | 0 | -1.6e-9 | +1.02e-4, +2.5e-5 |
  | Benzol | RI-J | -230.5362806561 | 0 | -2.3e-10 | -3.03e-4, -2.5e-5 |
  | Benzol | RI-JK | -230.5358608794 | +1e-10 | -3.9e-10 | +1.16e-4, +9.7e-6 |
  | Naphthalin | RI-J | -383.0669931516 | 0 | -1.9e-9 | -4.18e-4, -2.3e-5 |
  | Naphthalin | RI-JK | -383.0663624198 | 0 | -1.4e-9 | +2.13e-4, +1.2e-5 |

  Kartesisch (H2O, PySCF `cart=True`): RI-J und RI-JK gleich auf 1e-10. GPU vs CPU (Benzol, erster J/K-Aufbau):
  max|dJ| 2.8e-14, max|dK| 8.3e-14, E_elec gleich auf 10 Stellen. Getestet nur auf diesen Molekuelen (Elemente H, C,
  N, O; Anthracen nur Zeiten). Nicht getestet: Li, Be, B, F, Ne, He (Hilfsbasisdaten geladen, aber kein Molekuel
  gerechnet), geladene Systeme, lineare Abhaengigkeiten in groesseren Hilfsbasen, Systeme ueber 24 Atome, andere
  Orbitalbasen als def2-SVP.

**Zeiten** (Intel Core Ultra 9 285K, 24 Kerne, RTX 5090, `-qm.scf_threshold 1e-8`, ms je Fock-Aufbau = J/K-Zeit / Iterationen;
"exakt" CPU = gespeicherter Tensor bei Benzol (1.3 GB), integral-direkt bei Naphthalin/Anthracen (Tensor ueber
`eri_max_memory_mb`); "exakt" GPU = Tensor auf dem Device (Q-G2; 8.4 bzw. 29 GB), Anthracen mit Device-SCF-Loop):

| System (n, naux JK) | Pfad | Setup | ms/Fock | Wand | max RSS |
|---|---|---|---|---|---|
| Benzol (114, 558) | exakt 1 Thr | ERI 1935 ms | 133 | 4.2 s | 1323 MB |
| | RI-J 1 Thr | RI 172 ms + ERI 1933 ms | 136 | 4.5 s | 1360 MB |
| | RI-JK 1 Thr | 391 ms | 47 | 1.3 s | 97 MB |
| | exakt 16 Thr | ERI 285 ms | 51 | 1.25 s | 1323 MB |
| | RI-JK 16 Thr | 56 ms | 4.5 | 0.27 s | 97 MB |
| | exakt GPU | ERI 3668 ms (inkl. CUDA-Init) | 0.9 | 4.1 s | 634 MB |
| | RI-JK GPU | 114 ms + Upload 80 ms | 12.5 | 0.64 s | 643 MB |
| Naphthalin (180, 894) | exakt 1 Thr (direkt) | - | 5601 | 131 s | 45 MB |
| | RI-JK 1 Thr | 1908 ms | 219 | 7.0 s | 281 MB |
| | exakt 16 Thr (direkt) | - | 406 | 10.0 s | 53 MB |
| | RI-J 16 Thr | 89 ms | 411 | 10.2 s | 199 MB |
| | RI-JK 16 Thr | 235 ms | 25 | 1.3 s | 282 MB |
| | exakt GPU | ERI 400 ms | 6.7 | 1.3 s | 517 MB |
| | RI-JK GPU | 236 ms + Upload 80 ms | 17.7 | 1.3 s | 841 MB |
| Anthracen (246, 1230) | exakt 1 Thr (direkt) | - | 14423 | 365 s | 59 MB |
| | RI-J 1 Thr | 1407 ms | 13480 | 357 s | 431 MB |
| | RI-JK 1 Thr | 5764 ms | 769 | 25.6 s | 655 MB |
| | exakt 16 Thr (direkt) | - | 1045 | 27.8 s | 74 MB |
| | RI-J 16 Thr | 193 ms | 1190 | 33.1 s | 447 MB |
| | RI-JK 16 Thr | 699 ms | 87 | 4.2 s | 655 MB |
| | exakt GPU | ERI 839 ms | 28.0 | 2.6 s | 756 MB |
| | RI-JK GPU | 716 ms + Upload 105 ms | 29.6 | 3.0 s | 1192 MB |

- **Ehrlich:** RI-J spart fuer HF nichts (der exakte Austausch rechnet J mit; nur Vorstufe fuer reine Funktionale).
  Auf der GPU ist RI-JK nicht schneller als der exakte Pfad, solange der 4-Index-Tensor in den Device-Speicher
  passt (Benzol 12.5 gegen 0.9 ms je Fock-Aufbau); dort zaehlt der Speicher (29 GB gegen 0.6 GB bei Anthracen). Auf
  dem Device sind die Kontraktionen bei Anthracen etwa 3x schneller als 16 CPU-Threads, bei Benzol langsamer
  (Latenz, Transfers je Aufbau; die Host-Zerlegung von D kostet bei Anthracen nur ~4 ms je Aufbau). Setup RI-JK
  Anthracen bei 8 Threads: 3-Index-Integrale 363 ms, Vorwaertssubstitution B = L^-1 T 894 ms.
- BLAS: die schweren RI-Schritte (Vorwaertssubstitution, g = T vec(D), J, K) sind eigene OpenMP-Schleifen. Der Agent
  begruendete das mit der einthreadigen Referenz-netlib-BLAS (Setup Naphthalin mit Eigen-GEMM 8.5 s statt 0.1 s). Diese
  Begruendung gilt fuer das Release-Binary nicht, siehe Korrektur im Abschnitt "BLAS (Korrektur)" unten; die Schleifen
  bleiben, sind aber nicht durch langsames BLAS erzwungen. **Ueberholt fuer Metrik-Loesung und J/K-Aufbau, siehe
  "Skalierung RI-JK" unten** (dort Level-3-Aufrufe je Thread).

### Skalierung RI-JK (Okt 2026)

**🤖 AI-generated, ⚙️ machine-tested; menschlicher Produktionstest steht aus.** Profil vorher (heptacen, n=510,
naux=2574, 16 Threads): "factor/B" 16.3 s, davon Cholesky der Metrik 0.1 s, der Rest die Vorwaertssubstitution
B = L^-1 T als skalare Schleife (~100 GFLOPS); auf der GPU war das die Haelfte der Wandzeit. Geaendert:

- **Metrik-Loesung** (`lowerSolveRows`, ebenso `upperSolveRowsLT`, `leftMultiplyRows` in `qm_ri.cpp`): je Thread ein
  Spaltenblock (256 Spalten) des zeilenweisen Tensors, darauf EIN Level-3-Aufruf (Eigen `solveInPlace<OnTheRight>`,
  mit `USE_BLAS` dtrsm der einfaedigen MKL), das OpenMP darueber ist unseres. Nur die Summationsreihenfolge aendert sich.
- **Optional auf dem Device:** `-qm.gpu_ri_fit off|on|auto` (Standard off; `cuda/qm_gpu_ri_fit.cu`, Eintrag
  `curcuma_cuda_qm_ri_lower_solve`): cublasDtrsm auf gestreamten Spaltenbloecken (2 Streams, gepinnte Puffer), B kommt
  auf den Host zurueck (Gradient, RI-MP2, CPU-Fallback brauchen ihn). Fehler vor dem ersten Schreiben: Host-Loesung,
  danach: Setup ohne Device wiederholt; beides gezaehlt als `qm_ri_fit`. Standard off, weil die RTX 5090 (FP64 1/64)
  nur knapp schneller ist als der neue Host-Pfad und unter n ~ 450 langsamer (Tabelle).
- **J/K-Aufbau CPU** (`RIJK::build`): ein Durchlauf ueber B je Aufbau statt drei: je Thread und Block von nq Zeilen
  Y = M^T V (dgemm), g_q = sum_k s_k v_k . X^q_k, J += g_q B^q, K +-= Y Y^T (dsyrk); gleiche Zerlegung wie auf dem
  Device. Vorher skalare Dreifachschleife fuer K plus zwei Speicherdurchlaeufe fuer J.
- **Orbitalfaktor statt Eigenzerlegung:** `runSCF` merkt sich sqrt(2) C_occ zu der Dichte, die es daraus baut;
  `buildJK` nimmt ihn, wenn es genau diese Matrix bekommt (bitweiser Vergleich), sonst wie bisher `factorDensity`
  (SAD-Start, Daempfung). CPU und GPU. Spart die n^3-Zerlegung (heptacen ~21 ms je Aufbau).
- **1e-Integrale** (`qm_integrals.cpp` S, T, V): Zeilen parallel (vorher seriell), V mit den Kernen innen in derselben
  Reihenfolge, also bitgleich (Naphthalin exakt HF E_elec auf 10 Stellen gleich, `ctest -L qm_1e` gruen);
  heptacen 3.5 s -> 0.28 s. Team-Groesse aus `-threads` (`ScopedBlasThreads` in `buildOneElectronIntegrals`).

Zeiten (s Wand, je eine Messung, 16 Threads, Core Ultra 9 285K, RTX 5090, Standard-SCF-Schwelle, `-qm.ri jk`; GPU mit
`-gpu cuda -qm.gpu_ri true -qm.gpu_scf true`):

| System (n, naux) | CPU vorher | CPU nachher | GPU vorher | GPU nachher | GPU nachher + `gpu_ri_fit on` |
|---|---|---|---|---|---|
| Benzol (114, 558) | 0.21 | 0.12 | 0.61 | 0.54 | 0.47 |
| Naphthalin (180, 894) | 0.90 | 0.46 | 0.98 | 0.74 | 0.89 |
| Anthracen (246, 1230) | 2.93 | 1.38 | 2.02 | 1.31 | 1.49 |
| Pentacen (378, 1902) | 16.7 | 6.55 | 9.16 | 5.16 | 5.18 |
| C20H42 (490, 2256) | 32.6 | 11.3 | 22.6 | 7.44 | 7.26 |
| Heptacen (510, 2574) | 65.7 | 23.0 | 31.9 | 13.8 | 13.4 |
| C40H82 (970, 4476) | nicht gemessen | 145.5 | nicht gemessen | - | 131.8 (Kontraktionen CPU, s.u.) |

- factor/B (Host / Device): Pentacen 3.41 -> 0.66 / 0.82 s, C20H42 12.9 -> 1.55 / 1.41 s, Heptacen 16.6 -> 2.14 /
  1.77 s, C40H82 23.9 (Host) / 15.3 s (Device). Device-Werte enthalten die CUDA-Initialisierung (~0.15 s), die sonst
  im ersten J/K-Aufbau anfaellt.
- ms je J/K-Aufbau CPU: Pentacen 381 -> 159, C20H42 728 -> 373, Heptacen 1264 -> 508, C40H82 5060; GPU: 107 -> 100,
  182 -> 165, 245 -> 226. Auf dem Device ist K (GEMM + SYRK, ~3 naux n^2 n_occ Flops, ~2e11 bei Heptacen) bei etwa der
  Haelfte der FP64-Spitzenleistung der Karte; dort gibt es ohne weniger Flops nichts zu holen. CPU: 8 gegen 16 Threads
  nur 173 gegen 159 ms (Pentacen), 24 Threads 175 ms (P-/E-Kerne, statische Verteilung fuer Determinismus).
- **Energien:** E_elec (10 Stellen, `CURCUMA_QM_SCF_TRACE=1`) fuer alle sechs Systeme und alle Pfade gleich dem
  Stand vorher, Ausnahme C20H42 GPU 1e-10 Eh; Iterationszahlen gleich. C40H82: Host- und Device-Loesung gleiche
  E_elec, aber 19 gegen 18 Iterationen (Konvergenztest auf dP, nur Summationsreihenfolge verschieden).
- **Speicher (Grenze fuer n):** Tensor 8 naux n^2 Byte, naux ~ 4.6 n (Alkane) bis 5.05 n (Acene), also ~37-40 n^3.
  Device (32 GB, 32.6 GB frei gemeldet): Tensor bis ~32.5 GB, also n <= ~930 (Acene) bis ~955 (Alkane); C40H82
  (n = 970, 33.7 GB) passt nicht, die Kontraktionen fallen gezaehlt auf die CPU (`qm_ri`). Host: Standard
  `-qm.ri_max_memory_mb 16000` begrenzt auf n <= ~735-755; darueber muss der Wert erhoeht werden (C40H82 mit 40000).
  Der RI-JK-Gradient braucht zusaetzlich naux n^2 fuer G. (Stand vor der gepackten Speicherung, naechster Abschnitt.)
- **Nicht getestet:** `-qm.gpu_ri_fit` Fehlerpfade (Device-Fehler vor/nach dem ersten Schreiben), ROCm/Vulkan (kein
  Eintrag, faellt gezaehlt zurueck), RI-JK mit `-qm.ri_metric eigen` bei grossen Systemen (`leftMultiplyRows` nur ueber
  `ctest -L qm_ri`), RI-Gradient und RI-MP2 ueber n = 180 hinaus, andere Thread-Zahlen als 8/16/24, andere CPUs/GPUs,
  wiederholte Messungen (alle Zeiten Einzelmessungen).

### Gepackter 3-Index-Tensor (Okt 2026)

**🤖 AI-generated, ⚙️ machine-tested; menschlicher Produktionstest steht aus.** (P|mn) und B^Q_mn sind symmetrisch in
(m, n); jede Zeile von `RIJK::tensor()` haelt jetzt nur das untere Dreieck, Spalte `pairIndex(m, k) = m(m+1)/2 + k`
(m >= k), n(n+1)/2 statt n^2 Spalten. Das ist die einzige Speicherform (kein Schalter, Begruendung unten).

- **Setup:** schreibt nur das Dreieck; im Diagonalblock A == B bleibt der Wert mit m >= n, also derselbe, den die volle
  Form (letzter Schreibzugriff) fuer beide Eintraege behielt. Metrik-Loesung (Host und `-qm.gpu_ri_fit`) wirkt nur auf
  den Hilfsindex und laeuft unveraendert auf den halb so vielen Spalten. `CURCUMA_QM_RI_DUMP` entpackt auf n^2.
- **Kontraktionen:** g_Q = sum_p T_Q,p w_p mit w = P_mm bzw. P_mk + P_km (`packContractionWeights`); J gepackt
  akkumuliert und am Ende entpackt (`unpackSymmetric`). K (CPU, `RIJK::build`): je Zeile wird das Dreieck in einen
  Puffer kopiert und als symmetrische Matrix mit V multipliziert (dsymm), dann SYRK wie bisher. Device
  (`cuda/qm_gpu_ri.cu`): gepackt hochgeladen, g/J als dgemv auf den gepackten Spalten, fuer K entpackt ein Kernel
  (`unpackPairs`) je Block von Zeilen in einen Puffer, danach GEMM + SYRK wie bisher. Plugin-Eintrag umbenannt in
  `curcuma_cuda_create_qm_ri_packed` (ein altes Plugin wird gezaehlt abgelehnt statt den Tensor falsch zu lesen).
- **Gradient:** G (naux x n(n+1)/2) ebenfalls gepackt, P darin symmetrisiert (bei exakt symmetrischem P ohne Wirkung);
  `weightBlock` liest gepackt, die Device-Gewichte entstehen weiter auf dem Host. **RI-MP2:** je Zeile entpackt, dann die
  bisherige Transformation (gleiche Zahlen).

Zeiten (s, 16 Threads, Core Ultra 9 285K, RTX 5090, `-qm.ri jk`, Standard-SCF-Schwelle; Median aus 3 Laeufen ausser
C40H82 = 1 Lauf; "vorher" = Binary von Commit 777d9708; GPU mit `-gpu cuda -qm.gpu_ri true -qm.gpu_scf true`):

| System (n, naux) | Tensor MB | Setup ms CPU | ms/J-K CPU | Wand CPU | Upload ms GPU | ms/J-K GPU | Wand GPU |
|---|---|---|---|---|---|---|---|
| Benzol (114, 558) | 58 -> 29 | 44 -> 36 | 2.5 -> 1.5 | 0.14 -> 0.12 | 81 -> 66 | 13.4 -> 13.0 | 0.56 -> 0.45 |
| Naphthalin (180, 894) | 232 -> 116 | 147 -> 112 | 12.0 -> 6.7 | 0.47 -> 0.34 | 80 -> 68 | 17.6 -> 16.7 | 0.75 -> 0.69 |
| Anthracen (246, 1230) | 596 -> 299 | 380 -> 284 | 38.8 -> 21.2 | 1.40 -> 0.95 | 99 -> 80 | 26.4 -> 26.1 | 1.33 -> 1.21 |
| Pentacen (378, 1902) | 2174 -> 1090 | 1425 -> 1023 | 161 -> 121 | 6.60 -> 5.17 | 212 -> 135 | 99.3 -> 97.5 | 5.15 -> 4.61 |
| C20H42 (490, 2256) | 4333 -> 2171 | 2775 -> 1777 | 372 -> 289 | 11.3 -> 8.8 | 364 -> 213 | 165 -> 160 | 7.48 -> 6.25 |
| Heptacen (510, 2574) | 5356 -> 2683 | 3851 -> 2589 | 510 -> 413 | 23.2 -> 18.7 | 444 -> 250 | 226 -> 228 | 13.8 -> 12.4 |
| C40H82 (970, 4476) | 33692 -> 16863 | 33320 -> 18180 | 5048 -> 4152 | 147.6 -> 123.9 | - -> 1235 | 5082 (CPU) -> 1717 | 146.2 -> 62.1 |

- **C40H82 passt jetzt auf die Karte:** vorher fiel die Kontraktion gezaehlt auf die CPU (`qm_ri`, 33.7 GB Tensor),
  jetzt 16.9 GB auf dem Device, Wand 146 -> 62 s (mit `-qm.gpu_ri_fit on` 58.3 s, factor/B 12.1 -> 8.1 s). Beide mit
  `-qm.ri_max_memory_mb 40000`: der Host-Standard 16000 reicht knapp nicht (16863 MB), die Meldung nennt die Groesse.
- **Grenzen fuer n** (8 naux n(n+1)/2 Byte ~ 18.4 n^3 Alkane, 20.2 n^3 Acene): Host-Standard 16000 MB n <= ~925 (Acene)
  bis ~955 (Alkane), vorher ~735-755; Device (32 GB, abzueglich ~1.1 GB Puffer) n <= ~1155 bis ~1195, vorher ~930-955.
  Geschaetzt aus der Formel, nur C40H82 gemessen. Der RI-JK-Gradient braucht zusaetzlich naux n(n+1)/2 fuer G.
- **Energien:** E_elec (10 Stellen) aller sieben Systeme, CPU und GPU, gleich dem Binary vorher, Ausnahme C20H42 GPU:
  vorher -2163.0352272371, jetzt ...372 (= CPU-Wert; 1e-10, war schon vorher die GPU/CPU-Differenz).
  Iterationen gleich bis auf C40H82 (CPU 19 -> 21, GPU 18; Konvergenztest auf dP, nur Summationsreihenfolge
  verschieden, E_elec gleich).
- **Kein Schalter:** nirgends langsamer als 5 % (schlechtester Fall Heptacen GPU ms/J-K +1 %, Wand -10 %); Energien
  innerhalb 1e-10 Eh. Die CPU gewinnt durch den halben Speicherstrom im K-Aufbau (dsymm auf dem Dreieck), die GPU vor
  allem durch Setup und Upload; der K-Aufbau auf dem Device bleibt FP64-rechenbegrenzt, der Entpack-Kernel kostet
  dort ~1-2 %.
- **Gradient und RI-MP2 gegen das Binary vorher** (`-qm.scf_threshold 1e-10`, 16 Threads): RI-MP2 (eigenes Setup und
  RI-JK-Tensor der SCF wiederverwendet), Naphthalin und Anthracen, CPU und GPU: E(MP2 corr) auf 12 Stellen gleich.
  RI-J-Gradient: max. Abweichung 1.4e-12 Eh/A. RI-JK-Gradient: 8.9e-11 (Naphthalin), 1.2e-10 (Anthracen), 2.7e-11 Eh/A
  (Benzol); die Energie gleich auf 1e-12. Eingegrenzt an Benzol: der alte Host-Pfad war der Ausreisser -- alt mit
  Device-Loesung (`-qm.gpu_ri_fit on`), neu Host, neu Device und `-qm.ri_metric eigen` liegen untereinander innerhalb
  3e-13 bis 4e-12, alt Host liegt 2.7e-11 daneben; Threadzahl (1/8/16) und SCF-Schwelle (1e-10/1e-11) aendern den
  Gradienten nur um <= 1.5e-12. Vermutete Ursache (nicht nachgewiesen): in der vollen Form wurden (m,n) und (n,m) in
  verschiedenen Spaltenbloecken geloest, B war nicht bitweise symmetrisch, und der Gradient (zweite Metrik-Loesung) ist
  darauf empfindlicher als die Energie. `ctest -L qm_ri` (Gradient vs PySCF, 5e-10 Eh/Bohr) gruen;
  `CURCUMA_QM_RI_DUMP` (H2O) bitgleich zum Binary vorher (md5).
- **Test angepasst:** `qm_ri_gradient_error` nutzte `-qm.ri_max_memory_mb 0.8`, mit dem halben Tensor (0.27 MB) passt
  der Gradient dort jetzt hinein; Grenze auf 0.4 gesenkt (geprueft: mit 0.8 kein Fehler mehr, mit 0.4 die Meldung).
- **Nicht getestet:** keine Messung von RI-Gradient/RI-MP2-Zeiten, `-qm.ri_metric
  eigen` nur ueber `ctest -L qm_ri`, ROCm/Vulkan (keine Implementierung), andere GPUs/CPUs, Device-Speichergrenze
  nicht ausgereizt (kein System zwischen n = 970 und ~1190 gerechnet).

### RI-JK ressourcensparend (Okt 2026)

**🤖 AI-generated, ⚙️ machine-tested; menschlicher Produktionstest steht aus.** Code: `qm_ri.cpp` Abschnitte 4a
(batched, FP32), 4a' (lokaler Fit), 5 (Engine-Anbindung), `cuda/qm_gpu_ri_lowmem.cu` + `qm_gpu_ri_lowmem_backend.h`
(Streaming, FP32 auf dem Device, batched Austausch auf dem Device), `qm_scf.cpp` (Konvergenztest bei gemischter
Genauigkeit). Tests `test_cases/qm_ri/` (`qm_ri_batched_*`, `qm_ri_auto_batched_*`, `qm_ri_mixed_*`, `qm_ri_local_*`,
`qm_ri_grad_batched_*`, `qm_gpu_ri_stream_*`, `qm_gpu_ri_mixed_*`, `qm_gpu_ri_batched_*`).

**Entwurf (vor der Umsetzung festgelegt), was exakt ist und was nicht:**
- B = L^-1 T mit vollem unteren Cholesky-Faktor L: jede Zeile von B haengt von allen Zeilen von T ab. Eine
  Vorwaertssubstitution ueber Hilfsfunktions-Batches braucht alle frueheren B-Zeilen; **reines Batching von B gibt
  keinen O(Batch)-Speicher**. Verworfen.
- Exakt und ohne B: **J** aus g_P = T_P . w(P), c = J^-1 g, J_mn = sum_P c_P T_P,mn (zwei Durchlaeufe ueber neu
  gerechnete T-Batches). **K** ueber halbtransformierte Integrale: Y^P = T^P V (V = sqrt(2) C_occ bzw. Dichtefaktor),
  X = L^-1 Y auf dem Hilfsindex, K = sum_Q X^Q X^Q^T. Dieselbe Groesse wie sum_Q (B^Q V)(B^Q V)^T, nur andere
  Summationsreihenfolge. Speicher naux n r statt naux n(n+1)/2 (Faktor ~n/(2 n_occ): Heptacen 2.6, C40H82 3.0); die
  Orbitalspalten sind unabhaengig, also teilt der Modus Y in Orbital-Batches auf, wenn es nicht unter das Limit passt
  (max(2, Batches) Integral-Durchlaeufe je Aufbau). Kosten: 2 Integral-Durchlaeufe + Y + Dreiecksloesung
  (naux^2 n r) je Aufbau.

**Schalter** (Modul `qm`, alle mit Hilfetext in der Registry):
- `-qm.ri_storage stored|batched|auto` (Standard auto): stored wie bisher; auto nimmt stored, wenn der Tensor unter
  `-qm.ri_max_memory_mb` passt (dann bitgleich zum bisherigen Verhalten), sonst batched statt des bisherigen Abbruchs.
  `-qm.ri_batch_mb` (0 = min(256 MB, Limit/8)) Groesse eines Zeilen-Batches.
- `-qm.gpu_ri_stream auto|on|off` (Standard auto), `-qm.gpu_ri_chunk_mb`: passt der Tensor nicht auf die Karte (oder
  `-qm.gpu_max_memory_mb`), bleibt er im Host-RAM (cudaHostRegister) und wird je Aufbau in Zeilen-Chunks gestreamt
  (ein Kopier- und ein Rechen-Stream, zwei Puffer, Ereignisse ordnen die Wiederverwendung); mode jk in EINEM
  Durchlauf (g, J, K je Chunk). Vorher fiel ein zu grosser Tensor gezaehlt auf die CPU.
- `-qm.ri_mixed_precision true|false|auto` (Standard false), `-qm.ri_fp32_switch` (1e-3): Semantik wie GFN2
  `-scf_mixed_precision`. B als FP32 (Setup spaltenweise, kein FP64-Spitzenwert). Stufe 1: FP32-GEMM/SYRK (CPU Eigen
  float, Device SGEMM/SSYRK), Teilsummen je Batch in FP64 addiert, g und J in FP64. Stufe 2 (ab dP < Schalter oder
  drei Iterationen ohne Halbierung): derselbe FP32-Tensor mit FP64-Arithmetik plus Korrektur dJ/dK = exakt -
  FP32-Tensor aus einem exakten batched FP64-Aufbau. Konvergenz wird nur auf einem exakten Aufbau akzeptiert (erreicht
  die korrigierte SCF dP < Schwelle, ist der naechste Aufbau exakt), der Schlussaufbau ist exakt. auto = FP32 nur,
  wenn der FP64-Tensor nicht ins Limit passt (Host, oder `-qm.gpu_max_memory_mb`), der FP32-Tensor samt Puffern aber
  schon.
- `-qm.ri_fit_radius <Bohr>` (0 = aus): lokaler Austausch-Fit, **eine Naeherung** (unten).

**Warum Stufe 2 und nicht sofort exakt** (gemessen, Benzol): FP32-Arithmetik hat einen Rauschboden bei
||dP|| ~ 2-4e-5; der Rundungsfehler des FP32-Tensors verschiebt den Fixpunkt um ~1e-6 Eh und die Dichte um ~1e-4
(Frobenius). Direkt auf exakte batched Aufbauten zu wechseln brauchte 12 teure Aufbauten (Benzol); mit der Korrektur
sind es 3 (Wechsel, Pruefung bei Konvergenz, Schlussaufbau). Die Korrektur haengt von P nur ueber die Rundung von B ab
(relativ 6e-8). DIIS-Eintraege aus Stufe 1 werden beim Wechsel verworfen (mit ihnen 6 Iterationen mehr).

**Energien** (`-qm.scf_threshold` Standard, E_elec auf 10 Stellen, Einzelmessungen):
- stored (Standard) gegen das Binary vom Ausgangscommit d38a9daa (gleicher Compiler, gleiche MKL): alle "SCF iter"-
  Zeilen von Benzol, Naphthalin, Anthracen, Pentacen, C20H42, Heptacen identisch (md5), CPU und GPU; Rohtensor
  `CURCUMA_QM_RI_DUMP` (H2O) md5-gleich.
- batched, gemischt (CPU und GPU, resident und gestreamt), gestreamt FP64: E_elec gleich stored auf 10 Stellen fuer alle
  sechs Systeme; Ausnahme C20H42 auf der GPU 1e-10 (die bekannte GPU/CPU-Differenz). ctest bei
  `-qm.scf_threshold 1e-10`: batched (auch geteilt in viele Batches, eigen-Metrik, RI-J), gemischt, lokal mit
  R = 1e6: |dE| = 0 auf 10 Stellen gegen stored, <= 1e-10 gegen PySCF.
- RI-Gradient nach batched SCF (eigener gespeicherter Tensor): Benzol 2.6e-11 Eh/Bohr gegen PySCF.

**Zeiten und Speicher, mittelgrosse Systeme** (s Wand / ms je J/K-Aufbau / max RSS MB, je eine Messung, 16 Threads,
Core Ultra 9 285K, RTX 5090, Standard-SCF-Schwelle, `-qm.ri jk`; Iterationen in Klammern, wenn sie von stored
abweichen; GPU mit `-gpu cuda -qm.gpu_ri true -qm.gpu_scf true`):

| System (n, naux) | stored CPU | batched CPU | gemischt CPU | stored GPU | gestreamt GPU (`on`) | gemischt GPU | batched GPU |
|---|---|---|---|---|---|---|---|
| Benzol (114, 558) | 0.11 / 1.6 / 149 | 0.49 / 38 / 200 | 0.23 / 9.4 / 183 (13) | 0.46 / 16 | 0.56 / 17 | 0.58 / 22 (13) | 0.82 / 51 |
| Naphthalin (180, 894) | 0.31 / 6.9 / 252 | 2.3 / 155 / 347 | 0.79 / 33 / 356 (16) | 0.67 / 19 | 0.79 / 22 | 1.15 / 41 (17) | 2.7 / 171 |
| Anthracen (246, 1230) | 0.83 / 22 / 461 | 6.3 / 378 / 635 | 1.9 / 79 / 677 (17) | 1.11 / 29 | 1.29 / 34 | 2.2 / 82 (17) | 6.7 / 391 |
| Pentacen (378, 1902) | 3.6 / 126 / 1429 | 22.0 / 1257 / 1281 | 7.5 / 289 / 1479 (19) | 3.5 / 106 | 3.7 / 113 | 6.4 / 219 (20) | 22.4 / 1269 |
| C20H42 (490, 2256) | 6.8 / 299 / 2515 | 28.9 / 1987 / 2098 | 12.0 / 598 / 2356 (14) | 5.2 / 176 | 5.4 / 177 | 9.8 / 422 (15) | 27.4 / 1870 |
| Heptacen (510, 2574) | 11.2 / 415 / 3020 | 55.6 / 3007 / 2421 | 20.4 / 729 / 2963 (21) | 8.2 / 244 | 8.5 / 247 | 14.5 / 463 (21) | 51.8 / 2786 |

- batched: 3-index-Integrale sind der groesste Posten (Heptacen 2.0 von 3.0 s je Aufbau, zwei Durchlaeufe; der
  Integralkern ist die lesbare McMurchie-Davidson-Form, nicht die getunte des exakten Pfads). Unter dem 16-GB-Limit
  passt Y hier in einen Orbital-Batch, deshalb spart batched bei diesen Groessen nur 0-20 % RSS; der Gewinn ist, dass
  es ueberhaupt unter jedem Limit laeuft (Y teilt sich auf).
- gemischt CPU: die FP32-Tensor-Aufbaeute kosten im Mittel 320 ms (Heptacen, beide Stufen) gegen 415 ms FP64, C60H122
  17 s gegen 24.7 s (Faktor 1.3-1.45, nicht 2; Ursache nicht untersucht, ein Teil sind die Kopien der gepackten Zeilen und die FP64-Akkumulation von g und J); dazu
  3 exakte batched Aufbauten. Unter diesen Groessen ist gemischt deshalb langsamer als stored und spart RSS nur, wenn
  das Limit Y klein haelt; Energie gleich stored auf 10 Stellen (alle sechs Systeme, CPU und GPU).
- gestreamt: bis 7 % langsamer als resident (Heptacen 247 gegen 244 ms); gemischt GPU: Device-J/K in FP32 98 ms je
  Durchlauf (Heptacen), die Wand bestimmen die 3 exakten Aufbauten auf der CPU (Integrale).

**Grosse Systeme: C60H122** (n = 1450, naux = 6696, FP64-Tensor 56.4 GB; vorher mit Standardlimit ein Abbruch und auf der
Karte nicht moeglich; E_elec stored CPU = -8176.9793892180):

| Lauf | Wand s | ms je J/K | Iter. | max RSS | Device-Peak | E_elec - stored |
|---|---|---|---|---|---|---|
| CPU Standard (auto -> batched, 16000 MB) | 1096 | 76006 | 13 | 17.0 GB | - | +1e-10 |
| CPU stored, `-qm.ri_max_memory_mb 60000` | 462 | 24716 | 13 | 57.3 GB | - | 0 |
| CPU gemischt, 32000 MB | 839 | 44210 (FP32-Aufbauten 17 s) | 15 | 33.1 GB | - | +1e-10 |
| CPU lokal R = 8 Bohr (Standardlimit) | 522 | 32538 | 14 | 20.7 GB | - | +1.38e-4 (7.6e-7/Atom) |
| GPU Standard (batched + Device-Austausch) | 872 | 59979 | 13 | 17.7 GB | 5.6 GB | -1e-10 |
| GPU gestreamt FP64, 60000 MB | 234 | 8092 | 13 | 57.2 GB | 3.3 GB | 0 |
| GPU gemischt (FP32 resident 28.2 GB), 32000 MB | 586 | 26776 | 16 | 33.4 GB | 32.1 GB | 0 |

- Gestreamt: 789 GB in 14 Durchlaeufen hochgeladen, 112 s Device-Zeit; die Karte ist FP64-rechenbegrenzt (K ~1e13 Flops
  je Aufbau), die Kopie verdeckt. Setup 85 s (Loesung B = L^-1 T auf dem Host 70 s).
- Gemischt GPU: vor der Speicherkorrektur (Y-Budget ohne den FP32-Tensor) 405 s bei 48.9 GB RSS ueber dem Limit; jetzt
  586 s bei 33.4 GB (mehr Orbital-Batches in den 3 exakten Aufbauten, 21 statt 6 Integral-Durchlaeufe).
- Das Limit `-qm.ri_max_memory_mb` deckt die RI-Felder; RSS liegt ~1 GB darueber (Fock, Dichte, MO, Eigenloeser).
  Lokal: Y_k haelt sich nicht an das Limit (1450 x 1444 x 241 x 8 B ~ 4 GB, nicht geprueft).

**Lokaler Fit (`-qm.ri_fit_radius`, Naeherung)**: Fehler E(R) - E(global) in Eh / Eh je Atom (mittlere Domaene in
Hilfsfunktionen), Standardschwelle; der Fehler ist immer positiv (K_lokal <= K_global, Schur-Komplement):

| System (Atome, naux) | R = 4 | 6 | 8 | 10 | 12 |
|---|---|---|---|---|---|
| C10H22 (32, 1146) | 2.6e-4 / 8.1e-6 (682) | 3.6e-5 / 1.1e-6 (862) | 5.5e-6 / 1.7e-7 (949) | 7.8e-7 / 2.5e-8 (1022) | 2.1e-7 / 6.5e-9 (1047) |
| Anthracen (24, 1230) | 1.2e-4 / 4.9e-6 (953) | 1.2e-5 / 4.9e-7 (1081) | 1.0e-6 / 4.4e-8 (1154) | 1.6e-7 / 6.7e-9 (1191) | 3.1e-9 / 1.3e-10 (1213) |
| (H2O)20 (60, 2260) | 2.2e-3 / 3.6e-5 (274) | 3.5e-4 / 5.8e-6 (420) | 5.0e-5 / 8.4e-7 (668) | 6.9e-6 / 1.1e-7 (950) | 2.0e-6 / 3.4e-8 (1217) |
| C20H42 (62, 2256) | 5.5e-4 / 8.8e-6 (836) | 1.2e-4 / 1.9e-6 (1090) | 3.6e-5 / 5.9e-7 (1236) | 1.4e-5 / 2.2e-7 (1381) | 5.0e-6 / 8.1e-8 (1455) |
| Pentacen (36, 1902) | 6.3e-4 / 1.8e-5 (1261) | 1.4e-4 / 3.8e-6 (1471) | 1.9e-5 / 5.4e-7 (1608) | 2.3e-6 / 6.5e-8 (1696) | 2.5e-7 / 7.0e-9 (1756) |
| C30H62 (92, 3366) | 8.2e-4 / 8.9e-6 (895) | 1.8e-4 / 2.0e-6 (1169) | 6.1e-5 / 6.6e-7 (1337) | 2.2e-5 / 2.4e-7 (1510) | 7.7e-6 / 8.4e-8 (1603) |

- Verfahren: P = L L^T (pivotisierter Cholesky, keine Lokalisierungsiterationen), Domaene je Spalte = Hilfsfunktionen auf
  Atomen innerhalb R der Kernatome (Pivot-Atom plus die Atome mit den groessten sum l_mk^2 bis 99.9 %); fuer eine Dichte
  sind robuster (Dunlap) und symmetrischer lokaler Fit gleich. Domaenen bleiben nach zwei gleichen Aufbauten fest (mit
  wandernden Domaenen konvergierte Benzol bei R = 4 nicht in 100 Iterationen); die Energie haengt daher etwas vom
  SCF-Pfad ab. Eine erste Variante (Kernatome mit >= 2 % Anteil) hatte 5-30x groessere Fehler (Pentacen R = 12:
  3.2e-4 Eh/Atom), verworfen. R = 1e6 gibt den globalen Fit (ctest `qm_ri_local_exact_benzene_jk`).
- Nutzen nur im batched-Betrieb: Y_k braucht sum_k n |F_k| statt naux n n_occ (C60H122 R = 8: Domaene 22 % von naux) und
  die lokalen Loesungen sind billiger; C60H122 32.5 statt 76 s je Aufbau. Die Integrale laufen weiter voll. J bleibt
  global. RI-Gradient mit R > 0: Fehler (`qm_ri_local_gradient_error`).

**C80H162** (n = 1930, naux = 8916, FP64-Tensor 132.9 GB): `-gpu cuda -qm.ri_max_memory_mb 140000`, Tensor gepinnt im
Host-RAM, je Aufbau gestreamt: 694 s Wand (Setup 237 s, davon Loesung B = L^-1 T auf dem Host 208 s), 13 Iterationen,
27.0 s je J/K-Aufbau (1861 GB hochgeladen, 375 s Device-Zeit), max RSS 133.2 GB, Device-Peak 3.5 GB,
E_elec = -11508.2009585823 (keine Vergleichsrechnung; kein Rueckfall). Vorher war die Device-Grenze n ~1190 und ohne
Erhoehung des Limits ein Abbruch. Andere Modi fuer C80H162 nicht gerechnet (batched CPU geschaetzt ~3x C60H122 je Aufbau).

**Nicht getestet / nicht implementiert:** ROCm/Vulkan (kein Eintrag, gezaehlter Rueckfall); batched mit
`-qm.ri_metric eigen` nur Formaldehyd; Streaming bei mode j nur ein Handlauf (Benzol, E_elec gleich stored auf 10 Stellen);
gemischte Genauigkeit mit eigen-Metrik (Device-Austausch dann CPU); RI-MP2/LMP2 nach batched SCF rechnen ihre eigene
gespeicherte Einrichtung (nur Code-Pfad, kein Test); Gradient mit lokalem Fit (Fehler); batched/lokal auf dem Device
nur die Austausch-Loesung (Integrale, Y, J auf der CPU); Wiederholungsmessungen (alle Zeiten Einzelmessungen auf einer
Maschine mit Desktop-Last ~1-2).

## WP10 -- Gradient (RI-J, RI-JK)

**Stand Okt 2026: analytischer Kerngradient der RI-HF-Energie, CPU und CUDA. 🤖 AI-generated, ⚙️ machine-tested;
menschlicher Produktionstest steht aus.** Code: `qm_ri.cpp` Abschnitt 6 (Herleitung als Kommentar),
`QMEngine::riGradientTwoElectron`, `cuda/qm_gpu_ri_grad.cu` + `qm_gpu_ri_grad_backend.h`; Tests `test_cases/qm_ri/`
(`qm_ri_grad_*`, `qm_ri_opt_*`, `qm_gpu_ri_grad_*`), Referenzen `reference_grad.json`
(`scripts/ri_gradient_reference.py`, `scripts/ri_opt_reference.py`).

- **Verfahren:** dE/dx = sum_P sum_mn G^P_mn (mn|P)^x - sum_PQ M2_PQ (P|Q)^x mit G^P = c_P P - 2 Gam^P,
  M2 = 1/2 c c^T - M (RI-J: Gam = M = 0); c = J^-1 g, Gam = 1/4 J^-1 (P (Q|..) P), M = 1/4 J^-1 Tr((R|..)P(S|..)P) J^-1
  (Weigend/Haeser/Patzelt/Ahlrichs 1998, Weigend 2002; dieselben Ausdruecke wie PySCF `auxbasis_response`). Die
  Metrik-Loesung laeuft auf der Dichtefaktor-Form V (F^T W) V^T (r^2 statt n^2 Spalten). Ableitungsintegrale:
  angehobene Hermite-Tabellen auf beiden Orbitalschalen, Hilfszentrum per Translationsinvarianz. 1e-Teil, W-Term
  (mit den RI-SCF-Orbitalenergien) und E_nn wie im exakten Pfad. RI-J: Austausch-Gradient exakt ueber
  `qmint::gradientTwoElectronExchange` (gleicher Quartett-Kern, nur Austauschdichte; der exakte Pfad ist bitgleich,
  Benzol 1 Thread md5 des `-dump_gradient` vorher/nachher gleich).
- **Schalter:** `-qm.gpu_ri_grad auto|on|off` (3-Index-Teil auf dem Device, Fallback gezaehlt als `qm_ri_grad`).
- **Validierung** (HF/def2-SVP, `-qm.scf_threshold 1e-10`, Eh/Bohr):
  - vs PySCF 2.14 (`density_fit(JKFIT)` bzw. `density_fit(JFIT, only_dfj=True)`, `nuc_grad_method()`, dieselben
    Hilfsbasisdateien, `conv_tol 1e-14`, `conv_tol_grad 1e-11`): max. Abweichung H2O 4.6e-11 (j, jk), NH3 2.5e-11,
    Formaldehyd 6.1e-11, Benzol 2.5e-11/2.7e-11, Naphthalin 2.5e-11/2.6e-11, H2O kartesisch 4.6e-11/4.5e-11;
    Formaldehyd jk mit `-qm.ri_metric eigen` (nichts verworfen) gleich. **Falle:** mit PySCF `conv_tol_grad 1e-8`
    lag die Referenz selbst 1e-9..2e-9 daneben (PySCF-FD bestaetigte unseren Wert erst mit enger Konvergenz).
  - vs ORCA 6.1 (`! HF def2-SVP RIJK def2/JK | RIJONX def2/J ExtremeSCF EnGrad`): H2O 6.9e-9 (j) / 7.5e-9 (jk),
    Formaldehyd 8.8e-9 / 8.3e-9, Benzol 4.1e-9 / 4.8e-9; ORCA weicht von PySCF genauso ab (eigene Schwellen).
  - Finite Differenz der eigenen RI-Energie (zentral, h = 1e-4 Bohr, `dump_qm_gradient --ri --fd`): H2O 2.7e-9 (j),
    2.4e-9 (jk), NH3 2.1e-9 / 1.6e-9, H2O kartesisch 3.0e-9 / 2.6e-9 (FD-Abbruchfehler-Niveau).
  - Summe der Komponenten <= 2.3e-12 Eh/Bohr (Naphthalin 2.3e-12, Formaldehyd 9e-13, sonst <= 2e-13).
  - `-opt -qm.ri jk` (tight): H2O E = -75.9612828916, Formaldehyd -113.7796614858, gleich dem PySCF-RI-JK-Minimum
    (scipy BFGS gtol 1e-7) auf 1e-12 Eh, Abstaende auf 5e-7 A (geomeTRIC nicht installiert).
  - GPU vs CPU (3-Index-Teil, `CURCUMA_QM_GPU_RI_GRAD_CHECK=1`): Benzol 1.4e-13, Anthracen 3.7e-13 (max |g| ~20).
- **Zeiten** (Intel Core Ultra 9 285K, 16 Threads, RTX 5090; `-qm.scf_threshold 1e-8`, 2e-Gradient in ms):

  | System (n) | Pfad | 1e | 2e gesamt | davon Dichten / 3-Index / Metrik | Wand -sp -gradient |
  |---|---|---|---|---|---|
  | Benzol (114) | exakt CPU | 33 | 1564 | - | 2.8 s |
  | | RI-J CPU | 34 | 1596 | 0.5 / 26 / 0.7, exakter Austausch 1568 | 2.9 s |
  | | RI-JK CPU | 35 | 56 | 17 / 37 / 1.1 | 0.32 s |
  | | exakt GPU | 3 | 173 | - | 0.76 s |
  | | RI-JK GPU | 3 | 64 | 21 / 42 (Host 14, Upload 5, Kernel 22) / 1.3 | 0.58 s |
  | Naphthalin (180) | exakt CPU | 116 | 10124 | - | 20.0 s |
  | | RI-JK CPU | 118 | 204 | 66 / 135 / 2.5 | 1.35 s |
  | | exakt GPU | 4 | 789 | - | 1.94 s |
  | | RI-JK GPU | 4 | 196 | 65 / 125 (Host 46, Upload 13, Kernel 65) / 4.8 | 1.27 s |
  | Anthracen (246) | exakt CPU | 340 | 43350 | - | 73.9 s |
  | | RI-JK CPU | 331 | 596 | 218 / 371 / 5.3 | 4.3 s |
  | | exakt GPU | 6 | 1960 | - | 4.6 s |
  | | RI-JK GPU | 6 | 470 | 195 / 266 (Host 100, Upload 28, Kernel 135) / 7.1 | 2.7 s |

- **Ehrlich:** RI-J spart beim HF-Gradienten nichts (der exakte Austausch-Gradient dominiert). Der Device-Kern fuer
  den 3-Index-Teil ist nur so schnell wie 16 CPU-Threads (Benzol 42 gegen 37 ms, Anthracen 266 gegen 371 ms); der
  Kernel selbst braucht bei Anthracen 135 ms, Host-Vorbereitung der Gewichtsbloecke (321 MB) und Upload 128 ms. Die
  Plugin-Bibliothek wird fuer sm_75 gebaut; auf der RTX 5090 (sm_120) wird beim ersten Lauf PTX JIT-kompiliert
  (~100 ms, danach im CUDA-Cache). Dichten und Metrikteil laufen nur auf dem Host.
- **Getestet:** H2O, NH3, Formaldehyd, Benzol (ctest gegen PySCF), Naphthalin (manuell gegen PySCF), Anthracen (nur
  Zeiten, GPU gegen CPU); kartesisch nur H2O; eigen-Metrik nur Formaldehyd jk, nichts verworfen.
- **Nicht getestet:** Li, Be, B, F, Ne, He; geladene Systeme; eigen-Metrik mit verworfenen Eigenwerten (dann ist der
  Gradient nicht die exakte Ableitung); Systeme ueber 24 Atome; ROCm/Vulkan (kein Device-Pfad); `-md` mit RI.
- **Nicht implementiert:** RI-MP2-Gradient (eigenes Paket, Fehler bleibt), Device-Pfad fuer Dichten und Metrikteil,
  RI-J-Austausch-Gradient auf dem Device (der exakte Q-G7-Kern rechnet J und K zusammen; mit RI-J laeuft der
  Austauschteil auf der CPU).

## WP11 -- MO-Transformation und MP2

**Stand Okt 2026: Schritt 1 (CPU und GPU) und Schritt 3 umgesetzt; 🤖 AI-generated, ⚙️ machine-tested, menschlicher
Produktionstest steht aus.** Schritt 2 (RI-MP2): eigener Abschnitt unten.

- **Code:** `qm_mp2.{h,cpp}` (Parameter, `qmint::mp2EnergyHost`, `QMEngine::computeMP2`), `mp2_method.{h,cpp}`
  (Methoden `mp2`, `scs-mp2`, `sos-mp2`), `cuda/qm_gpu_mp2.cu` + Naht `qm_gpu_mp2_backend.h` (Eintrag
  `curcuma_cuda_qm_mp2_energy`). Parameter: `-qm.frozen_core` (Standard aus; H-He 0, Li-Ne 1s, Na-Ar [He]2s2p ...),
  `-qm.mp2_batch_mb` (2000), `-qm.mp2_max_memory_mb` (32000), `-qm.gpu_mp2 auto|on|off`. `-sp` gibt E(HF)+E_corr der
  gewaehlten Variante, `-verbosity 2` zusaetzlich E_OS, E_SS, SCS (1.2/0.333), SOS (1.3) und die Zeitaufteilung.
- **Algorithmus:** vier Ein-Index-Transformationen in Batches ueber i (Schritt 1 als hohe schmale GEMM, OpenMP ueber
  Zeilenbloecke; Schritte 2-4 je i; Energie je i, weil (ia|jb) und (ib|ja) im selben i-Block liegen). Zusatzspeicher
  batch*(n^3 + o v^2) auf der CPU, batch*(n^3 + n^2 v) auf der GPU, nie der o^2 v^2-Tensor.
- **Tensorquelle:** gespeicherter Host-Tensor der SCF; sonst Download des Device-Tensors; sonst Neuaufbau mit
  `buildERI`, wenn n^4*8 Byte <= `-qm.mp2_max_memory_mb` (bei integral-direkter SCF), darueber Abbruch mit Fehler.
  **Kein integral-direktes MP2**, also Grenze n^4 im Hostspeicher. GPU: gespeicherter Device-Tensor (Q-G1/Q-G2,
  `QMGpuJK::deviceERI`), sonst Host-Tensor hochgeladen; Fehler -> CPU, gezaehlt als `qm_mp2`.
- **Gradient:** seit Okt 2026 analytisch (mp2, scs-mp2, sos-mp2, mit und ohne frozen core), eigener Abschnitt
  "WP11 -- MP2-Gradient" unten; der fruehere Fehlertest `qm_mp2_opt_error` ist durch `qm_mp2_opt_*` ersetzt.
- **Validierung** (def2-SVP sphaerisch, `-qm.scf_threshold 1e-10`): `ctest -L qm_mp2` (30 Tests) gegen PySCF 2.14
  `mp.MP2` (RHF conv_tol 1e-12, `scripts/mp2_reference.py` -> `test_cases/qm_mp2/reference.json`) auf H2O, NH3, CH4,
  HF, H2O/Formaldehyd verzerrt, Formaldehyd, Benzol, je voll und frozen core: alle E_OS, E_SS, E_corr, SCS, SOS
  innerhalb 1e-8, groesste Abweichung 1.8e-9 Eh (Benzol voll). ORCA 6.1 `! RHF MP2|SCS-MP2 def2-SVP ExtremeSCF
  NoFrozenCore|FrozenCore` (druckt 9 Nachkommastellen) auf denselben 16 Faellen: E_corr und E_SCS innerhalb 2.0e-9 Eh
  (Tabelle unten). Mit `VeryTightSCF` lag ORCA selbst bis 7.6e-8 von PySCF entfernt (Formaldehyd verzerrt);
  ExtremeSCF ist fuer den Vergleich noetig. Naphthalin (180 BF, nicht im ctest): 5.0e-10 (voll), 5.1e-10 (fc) gegen
  PySCF. GPU gegen CPU: E_corr auf allen 12 gedruckten Stellen gleich (Benzol, Naphthalin, Formaldehyd).
- **Zeiten** (RTX 5090, Hostspeicher 256 GB; CPU-BLAS: laut Korrektur im Abschnitt "BLAS (Korrektur)" bindet das Binary
  MKL sequential, einfaedig, nicht netlib). MP2-Anteil ohne SCF, Transformation / Energie:

  | System (n, o, v) | CPU 1 Thread | CPU 16 Threads | GPU (gespeicherter Device-Tensor) |
  |---|---|---|---|
  | Benzol voll (114, 21, 93) | 267 / 6.1 ms (Nachmessung, siehe Korrektur) | 68 / 1.0 ms (Nachmessung) | 64 / 0.3 ms |
  | Benzol fc (114, 15, 93) | 2458 / 3.6 ms (ueberholt) | 656 / 0.5 ms (ueberholt) | 63 / 0.2 ms |
  | Naphthalin voll (180, 34, 146) | 57171 / 41 ms (ueberholt, nicht nachgemessen) | 551 / 7.1 ms (Nachmessung) | 199 / 0.7 ms |
  | Naphthalin fc (180, 24, 146) | nicht gemessen | 7364 / 3.0 ms (ueberholt) | 135 / 0.4 ms |

  **Korrektur (Okt 2026, Nachmessung nach dem Merge, `-method mp2 -gpu none -qm.scf_threshold 1e-10`, Benzol, n=1):**
  Transformation / Energie 267 / 6.1 ms (1 Thread), 117 / 2.0 ms (4 Threads), 68 / 1.0 ms (16 Threads); Naphthalin
  (`-qm.scf_direct off`, 16 Threads) 551 / 7.1 ms. Die mit "ueberholt" markierten CPU-Werte der Tabelle sind 10 bis 20 mal
  groesser als diese Nachmessung und waren nicht reproduzierbar; Ursache nicht geklaert (Vermutung, nicht geprueft:
  Messung unter Last durch einen gleichzeitigen Build). Die GPU-Werte wurden nicht neu gemessen.

  CPU-Naphthalin mit `-qm.scf_direct off` (Tensor gespeichert); mit Standard (direkte SCF) kommt der Neuaufbau des
  Tensors dazu (8398 MB, 1572 ms). GPU mit `-qm.gpu_direct on`: Neuaufbau auf dem Host + Upload (599 ms) vor der
  Transformation. Device-Speicher Naphthalin (gerechnet, nicht gemessen): Tensor 8398 MB + Batch 34*84.5 MB.
- **Nicht getestet:** offene Schale (nicht implementiert), Elemente ausser H, C, N, O, F (Li-Ne-frozen-core-Regel nur auf
  C, N, O, F gegen PySCF/ORCA geprueft), andere Basen als def2-SVP, kartesisch 6d, Systeme ueber 180 BF, ROCm/Vulkan
  (keine Implementierung, `-gpu rocm` faellt gezaehlt auf die CPU), geladene Molekuele.
- **Nicht implementiert gegenueber ORCA/PySCF:** integral-direktes MP2, MP2-Dichte als Ausgabe (die relaxierte
  Dichte existiert nur intern im Gradienten), frei waehlbare SCS-Faktoren, Ausnutzung der Permutationssymmetrie in der Transformation (Faktor ~2 offen).

Vergleich (Eh; curcuma `-qm.scf_threshold 1e-10 -qm.threads 16`; ORCA ExtremeSCF, 9 gedruckte Stellen; "SCS - ORCA"
gegen `! SCS-MP2`, das zusammen mit E_corr die OS/SS-Aufteilung prueft):

| Molekuel | E_corr voll | - PySCF | - ORCA | SCS - ORCA | E_corr fc | - PySCF | - ORCA | SCS - ORCA |
|---|---|---|---|---|---|---|---|---|
| H2O | -0.1999842841 | 6.0e-10 | -5.2e-11 | 1.3e-10 | -0.1973855431 | 6.3e-10 | -1.4e-10 | 7.5e-10 |
| NH3 | -0.1936363035 | -9.3e-11 | 5.3e-10 | 4.2e-10 | -0.1908720455 | -9.3e-11 | 4.8e-10 | 8.4e-10 |
| CH4 | -0.1647349435 | 1.5e-10 | 4.6e-10 | 1.1e-10 | -0.1612389364 | 1.4e-10 | 6.0e-10 | 5.4e-10 |
| HF | -0.2015058020 | -4.1e-10 | -4.9e-11 | 1.9e-10 | -0.1994674300 | -4.1e-10 | 3.1e-11 | 1.4e-10 |
| H2O_distorted | -0.2000156496 | -3.1e-10 | 3.5e-10 | 6.5e-10 | -0.1974146809 | -3.2e-10 | 9.9e-11 | 3.3e-10 |
| formaldehyde_distorted | -0.3259997291 | 3.0e-10 | 9.0e-10 | 1.1e-09 | -0.3208916624 | 2.9e-10 | 6.3e-10 | 1.1e-09 |
| formaldehyde | -0.3214541495 | -1.7e-09 | 5.4e-10 | 6.1e-10 | -0.3162181368 | -1.7e-09 | 1.2e-09 | 6.1e-10 |
| benzene | -0.7970301248 | -1.8e-09 | 1.2e-09 | 1.1e-09 | -0.7794624731 | -1.8e-09 | 1.9e-09 | 2.0e-09 |

Die Restabweichungen (1e-10 bis 2e-9) liegen in der Groessenordnung der SCF-Konvergenz der Referenzen (PySCF
conv_tol_grad = sqrt(conv_tol) = 1e-6) und der Schwarz-Schwelle 1e-12; nicht weiter zerlegt.

Urspruenglicher Plan:

- **Schritt 1 (kanonisch, ohne RI):** halbtransformierter Tensor (mn|ia), dann (ia|jb). O(n^5), Speicher
  O(o v n^2). Closed-shell RMP2 mit eingefrorenem Rumpf (`-qm.frozen_core`). Validierung gegen PySCF `mp.MP2`
  und ORCA `MP2` auf den vorhandenen Testmolekuelen (H2O, NH3, CH4, HF, Formaldehyd, Benzol zur Skalierung).
- **Schritt 2 (RI-MP2):** (ia|P)-Tensor aus WP10, Energie ueber `B_ia^P`. Vergleich gegen den kanonischen Wert
  zeigt den RI-Fehler; gegen ORCA `RI-MP2` mit `def2-SVP/C` als Test.
- **Schritt 3 (SCS-/SOS-MP2):** nur Skalierungsfaktoren auf die Spinkomponenten; billig, aber erst nach der
  Spin-Aufteilung der Energie.
- **GPU:** Transformation und Energiesumme sind GEMM; der device-residente Tensor aus Q-G1/Q-G2 ist der
  Startpunkt. Zuerst die CPU-Version fertigstellen und validieren.
- **Offen:** Gradient (Z-Vektor, relaxierte Dichte) ist ein eigenes Arbeitspaket; ohne ihn ist MP2 nur
  Einzelpunkt, kein `-opt`.

## WP11 -- MP2-Gradient

**Stand Okt 2026: analytischer Kerngradient des kanonischen closed-shell MP2 (mp2, scs-mp2, sos-mp2; voller und
eingefrorener Rumpf), CPU; der 2e-Ableitungsteil auch auf CUDA. 🤖 AI-generated, ⚙️ machine-tested, menschlicher
Produktionstest steht aus.**

- **Code:** `qm_mp2_grad.cpp` (Herleitung im Kommentarblock am Dateianfang; `qmint::mp2GradientHost`,
  `QMEngine::computeMP2Gradient`, `QMEngine::mp2G4`, `qmint::mp2ActiveToCartesian`), `qmint::gradientTwoElectronMP2`
  in `qm_integrals.cpp` (deklariert in `qm_mp2.h`; teilt `gradientTwoElectronDriver` mit dem RHF-Gradienten),
  `cuda/qm_gpu_mp2_grad.cu` + Naht `qm_gpu_mp2_grad_backend.h` (Eintrag `curcuma_cuda_create_qm_mp2_grad2e`, Kopie
  des Q-G7-Kerns mit anderem Quartettgewicht), `mp2_method.cpp` (Gradient in Eh/Angstrom). `-sp -gradient`,
  `-opt -method mp2` und alles, was `getGradient()` nutzt, laufen damit.
- **Verfahren:** Lagrange-/Z-Vektor-Formulierung (Handy-Schaefer) mit relaxierter Dichte, energiegewichteter Dichte und
  nicht-separabler 2-Teilchen-Dichte Gamma = sum T~ C C C C im AO-Raum (n^4 auf dem Host); Frozen-Core ueber die
  Rotationen Rumpf-aktiv (w_Ij = L1_Ij/(e_j - e_I)) und Rumpf-virtuell (im Z-Vektor). Amplituden und
  Halbtransformationen in Batches ueber i wie die Energie (`-qm.mp2_batch_mb`), GEMMs mit `qmint::gemmABt`.
  Z-Vektor: vorkonditioniertes CG mit Diagonale e_a - e_k, J/K vom exakten Pfad (gespeicherter Tensor bzw. Device);
  nach integral-direkter SCF mit dem fuer MP2 neu gebauten Tensor (der direkte Builder hat einen Schwellen-Rauschboden,
  Naphthalin blieb damit nach 100 Iterationen bei 1.7e-9 stehen). Abbruch bei max|r| <= 1e-10, Fehler ueber 1e-6.
- **Validierung** (def2-SVP sphaerisch, `-qm.scf_threshold 1e-10`): `ctest -L qm_mp2` (`qm_mp2_grad_*`,
  `qm_mp2_opt_*`, `qm_gpu_mp2_grad_*`, 29 neue Tests). Referenz `test_cases/qm_mp2/reference_grad.json` aus
  `scripts/mp2_gradient_reference.py` (PySCF 2.14): **PySCFs analytischer MP2-Gradient
  (`mp.MP2(mf).nuc_grad_method()`) weicht von der finiten Differenz der eigenen PySCF-MP2-Energie um bis zu 1.1e-7
  Eh/Bohr ab** (H2O; unabhaengig von der CPHF-Toleranz 1e-9 oder 1e-13; Ursache nicht untersucht). Deshalb ist die
  enge Referenz die Vierpunkt-Differenz der PySCF-Energie (h = 1e-3 Bohr, RHF conv_tol 1e-13).
  curcuma gegen diese FD: H2O, NH3, CH4, HF, H2O/Formaldehyd verzerrt, Formaldehyd je voll und fc <= 6.5e-10,
  Benzol 1.3e-9 (voll) / 1.4e-9 (fc) Eh/Bohr; gegen PySCF analytisch <= 1.1e-7; Energie <= 1.1e-11 Eh.
  Gegen die zentrale Differenz der eigenen curcuma-Energie (h = 1e-4 Bohr): H2O verzerrt 2.8e-9 (voll und fc), NH3
  1.6e-9 (voll und fc), SCS-MP2 Formaldehyd verzerrt 3.7e-9, SCS-MP2 H2O fc 2.6e-9, SOS-MP2 NH3 fc 1.6e-9.
  Summe ueber Atome <= 1e-12. ORCA 6.1 (`! MP2 def2-SVP NoRI ExtremeSCF EnGrad NoAutoStart [No]FrozenCore`, .engrad)
  gegen curcuma: H2O verzerrt 2.6e-8 / 2.0e-8 (fc), Formaldehyd verzerrt 2.4e-8 / 2.5e-8, Benzol 7.7e-8 / 7.3e-8,
  Naphthalin (voll, nicht im ctest) 7.6e-8 Eh/Bohr (Energie 5.0e-10 Eh); ORCA liegt um denselben Betrag neben der PySCF-FD (2.0e-8 bis 7.6e-8), die Abweichung ist also ORCAs
  eigene Genauigkeit (vermutlich Z-Vektor-Toleranz, nicht geprueft). `-opt -method mp2` (Preset tight) von den
  verzerrten Startstrukturen: H2O E - E_min(PySCF) = 1.6e-13 Eh, max |d r_ij| = 1.7e-7 A; Formaldehyd 4.7e-11 Eh,
  5.7e-6 A (PySCF-Minimum: scipy BFGS gtol 1e-7, `scripts/mp2_opt_reference.py`). GPU (2e-Teil) gegen CPU-2e-Teil
  (`qm_gpu_mp2_grad_*`): Formaldehyd verzerrt 1.1e-14 / 2.5e-14 (fc), Benzol 1.4e-13 / 1.3e-13 Eh/Bohr.
- **Zeiten** (ms, RTX 5090, CPU-BLAS = MKL sequential (siehe Korrektur), `curcuma -sp X -method mp2 -gradient -qm.scf_threshold 1e-10
  -qm.threads N -gpu none|cuda -verbosity 2`, Zeile "MP2 gradient timing"; Maschine sonst unbelastet; GPU = 2e-Ableitungsteil, 1e-Teil, SCF und Z-Vektor-J/K auf dem Device,
  Amplituden/Ruecktransformation auf der CPU mit 16 Threads):

  | System (n, o, v), voll | SCF | Amplituden | Rueck-transf. | Z-Vektor (It.) | Dichten | 2e-Abl. | 1e-Abl. | Wand |
  |---|---|---|---|---|---|---|---|---|
  | Benzol (114, 21, 93), 1 Thread | 2418 | 435 | 1179 | 1474 (10) | 930 | 21100 | 421 | 30.3 s |
  | Benzol, 16 Threads | 937 | 120 | 334 | 606 (10) | 158 | 1605 | 34 | 4.3 s |
  | Benzol, GPU | 42 | 113 | 312 | 13 (10) | 109 | 396 | 18 | 1.8 s |
  | Naphthalin (180, 34, 146), 1 Thread | 158196 + 9779 Tensor | 4064 | 10635 | 11282 (14) | 6220 | 136051 | 1564 | 346 s |
  | Naphthalin, 16 Threads | 11032 + 1585 Tensor | 1051 | 2872 | 4775 (14) | 1025 | 11372 | 122 | 34.9 s |
  | Naphthalin, GPU | 236 (+ Download 8.4 GB) | 1024 | 2887 | 118 (14) | 743 | 2639 | 39 | 10.3 s |

  Naphthalin CPU: integral-direkte SCF (Tensor 8398 MB > `-qm.eri_max_memory_mb` 4000), Tensor fuer MP2 neu gebaut.
  GPU-2e-Teil Naphthalin: 2PDM nach kartesisch 682 ms (Host), Upload 780 ms, Kerne 1157 ms, 10.4 GB Device-Speicher;
  Benzol 78 / 122 / 192 ms. GPU gegen CPU, ganzer Gradient Naphthalin: 2.6e-12 Eh/Bohr. Speicher Host (gerechnet):
  Tensor n^4 + 2PDM n^4 (Naphthalin 2 x 8.4 GB) plus o v o v- und o n o v-Felder.
- **Threadzahl (31. Okt 2026, Verhaltensaenderung):** `mp2GradientHost` schnitt die erste Viertel-Transformation in
  threadzahlabhaengige Zeilenbloecke, der Gradient war 1 gegen 4/8 Threads nicht bitgleich (md5, H2O_distorted, NH3); jetzt
  feste 8192-Zeilen-Bloecke: 1/4/8 Threads bitgleich, gegen a06631c0 max |dg| 5.0e-15 Eh/Bohr, Energien bitgleich
  (`mp2EnergyHost` unveraendert, hat dieselbe Aufteilung; dort 1/4/8 Threads auf H2O_distorted und NH3 bitgleich, groessere
  Systeme nicht geprueft). Messungen: [WP12](WP12_doppelhybride.md), Abschnitt "Entscheidungen vom 2026-10-31".
- **Nicht getestet:** offene Schale (nicht implementiert), Elemente ausser H, C, N, O, F, andere Basen als def2-SVP,
  kartesisch 6d, geladene Molekuele, Systeme ueber 180 BF, `-md` und `-hessian` mit mp2 (laufen ueber `getGradient()`,
  nicht geprueft), ROCm/Vulkan (keine Implementierung, `-gpu rocm` faellt gezaehlt auf die CPU), entartete
  Rumpf-/Valenz-Orbitalenergien (w_Ij teilt durch e_j - e_I).
- **Nicht implementiert:** RI-MP2-Gradient (Fehler wie zuvor), MP2-Gradient mit RI-J/RI-JK-SCF (Fehler),
  integral-direkter Gradient (n^4 Tensor und n^4 2-Teilchen-Dichte auf dem Host; Grenze `-qm.mp2_max_memory_mb` je
  Feld), GPU fuer Amplituden, Ruecktransformation und Z-Vektor-Algebra (nur die J/K des Z-Vektors laufen ueber den
  vorhandenen Device-Pfad, wenn die SCF auf dem Device lief), Dichte-/Dipolausgabe.

## WP11 Schritt 2 -- RI-MP2

**Stand Okt 2026: closed-shell RI-MP2 (CPU und CUDA), 🤖 AI-generated, ⚙️ machine-tested; menschlicher
Produktionstest steht aus.** Code: `qm_rimp2.cpp` (`QMEngine::computeRIMP2`, `qmint::riMp2EnergyHost`,
`qmint::gemmABt`), `RIJK::fitRows` in `qm_ri.cpp`, `cuda/qm_gpu_rimp2.cu` + Naht `qm_gpu_rimp2_backend.h`
(Eintrag `curcuma_cuda_qm_rimp2_energy`), Methoden in `mp2_method.cpp`. Tests `test_cases/qm_mp2/` (`qm_rimp2_*`),
Referenzen `scripts/rimp2_reference.py` -> `test_cases/qm_mp2/reference_ri.json`.

- **Hilfsbasis (Methodenvariante):** curcuma hat keine def2-SVP/C-Datei. Nach Betreiberentscheidung dient
  def2-universal-JKFIT (Standard) bzw. -JFIT als Korrelations-Hilfsbasis (`-qm.ri_mp2_aux`, auch def2/JK, def2/J
  oder ein Dateipfad). Das ist **nicht** RI-MP2/def2-SVP/C anderer Programme; die Ausgabe sagt es ("RI-MP2 (aux
  def2-universal-JKFIT; ... method variant, not RI-MP2/def2-SVP/C)").
- **Methoden:** `ri-mp2`, `scs-ri-mp2`, `sos-ri-mp2` (gleiche OS/SS-Aufteilung und Faktoren wie kanonisch),
  `-qm.frozen_core` wie kanonisch. Die SCF bleibt exakt (`-qm.ri off`); mit `-qm.ri jk` ist E(HF) die RI-JK-Energie,
  die Ausgabe markiert das (`E(HF) = ... (RI-JK def2-universal-JKFIT SCF)`), und bei gleicher Hilfsbasisdatei wird
  der RI-Tensor der SCF wiederverwendet (kein zweites Setup).
- **Algorithmus:** eigenes RIJK-Setup im Modus j (rohe (Q|mn) + Metrikfaktor, Cholesky Standard, `eigen` ueber
  `-qm.ri_metric`), dann (Q|ia) = Co^T T^Q Cv je Hilfsfunktion, Metrikfaktor auf die o v Spalten (`fitRows`; billiger
  als auf n^2 Spalten), Transposition zu B[i][a][P]. Energie CPU: Paare i >= j, je Paar eine v x v-Kontraktion
  K = B_i B_j^T, Faktor 2 fuer i != j, Paarenergien in fester Reihenfolge summiert (unabhaengig von der Threadzahl).
  GPU: B einmal hochgeladen, je i-Block eine cuBLAS-DGEMM (FP64), Paarenergie-Kernel wie kanonisch (ohne
  i <-> j-Symmetrie). Kein Batching auf der CPU noetig (Zusatzspeicher je Thread v^2).
- **Speicher:** Setup-Tensor naux n(n+1)/2 (gepackt seit Okt 2026; Grenze `-qm.ri_max_memory_mb`, Fehler mit Groessenangabe), B o v naux
  zweimal (vor/nach Transposition, gleiche Grenze). Anthracen: Tensor 596 MB (gepackt 299 MB), B 92 MB.
- **GPU:** `-qm.gpu_mp2 auto|on|off` gilt fuer den kanonischen und den RI-Pfad (kein eigener Schalter); Fehler ->
  CPU, gezaehlt als `qm_ri_mp2`. Nur die Energie (GEMM + Summe) laeuft auf dem Device; Integrale, Metrik und
  B-Transformation bleiben auf dem Host.
- **Gradient:** keiner; `-opt`/`-gradient` brechen mit "ri-mp2 has no analytic gradient" ab (`qm_rimp2_opt_error`).
- **BLAS (gemessen):** das Release-Build linkt die Referenz-netlib-libblas (einfaedig); Eigen-GEMM laeuft darueber
  mit 2.7 GFLOPS (199 x 1230 x 9353). Eigens interne GEMM (ohne `EIGEN_USE_BLAS`) erreicht 58-60 GFLOPS auf einem
  Thread, ist aber nicht pro Datei abschaltbar (`EIGEN_USE_BLAS` gilt global ueber PCH/`curcuma_eigen_config.h`,
  gemischte Instanziierungen waeren ein ODR-Verstoss). Deshalb `qmint::gemmABt` (4x4-Registerbloecke, Skalarprodukte
  entlang des zusammenhaengenden Index): 29 GFLOPS auf 1 Thread, 369 GFLOPS mit 16 Threads in der Paarschleife
  (netlib in derselben Schleife: 5.6 bzw. 78 GFLOPS). BLAS-Auswahl des Builds nicht geaendert.

**Validierung** (HF/def2-SVP sphaerisch, `-qm.scf_threshold 1e-10`, curcuma 16 Threads). PySCF 2.14: Hilfsbasis aus
curcumas eigenen `.turbomole`-Dateien gelesen (gleiche Daten), manuelle Auswertung ueber `df.incore.cholesky_eri`,
gegen `mp.dfmp2.DFMP2` <= 3e-15 Eh; gleiche naux (z.B. Benzol 558). ORCA 6.1: `! RHF RI-MP2 def2-SVP ExtremeSCF
NoFrozenCore|FrozenCore` + `%basis AuxC "def2/JK" end` (Ausgabe: "AuxC ... def2/JK", 113 Funktionen fuer H2O, keine
J/K-Hilfsbasis in der SCF, also exakte SCF), SCS ueber `%mp2 DoSCS true end` (dann ist die gedruckte
RI-MP2-Korrelationsenergie die SCS-Energie). RI-Fehler = RI-MP2 minus kanonisches MP2 (WP11, gleiche SCF).

| Molekuel | | E_corr RI-MP2 (JKFIT) | - PySCF | - ORCA | SCS - ORCA | RI - kanonisch (Eh) | je Atom |
|---|---|---|---|---|---|---|---|
| H2O | voll | -0.2000132566 | 5.9e-10 | 4.0e-10 | 4.9e-11 | -2.90e-05 | -9.7e-06 |
| H2O | fc | -0.1974146621 | 6.3e-10 | -9.1e-11 | 5.2e-11 | -2.91e-05 | -9.7e-06 |
| NH3 | voll | -0.1935490206 | -9.4e-11 | 3.5e-10 | 1.8e-10 | 8.73e-05 | 2.2e-05 |
| NH3 | fc | -0.1907850152 | -9.4e-11 | 8.3e-10 | 1.8e-10 | 8.70e-05 | 2.2e-05 |
| CH4 | voll | -0.1645809798 | 1.5e-10 | 1.9e-10 | -6.6e-11 | 1.54e-04 | 3.1e-05 |
| CH4 | fc | -0.1610852210 | 1.4e-10 | -2.2e-11 | 6.9e-10 | 1.54e-04 | 3.1e-05 |
| HF | voll | -0.2015397788 | -4.1e-10 | 2.0e-10 | -1.3e-10 | -3.40e-05 | -1.7e-05 |
| HF | fc | -0.1995015462 | -4.1e-10 | -1.9e-10 | 2.0e-10 | -3.41e-05 | -1.7e-05 |
| H2O verzerrt | voll | -0.2000448210 | -3.1e-10 | -1.6e-11 | 4.3e-10 | -2.92e-05 | -9.7e-06 |
| H2O verzerrt | fc | -0.1974439973 | -3.2e-10 | 7.0e-10 | 5.8e-10 | -2.93e-05 | -9.8e-06 |
| Formaldehyd verzerrt | voll | -0.3259877876 | 3.0e-10 | 3.5e-10 | 1.1e-09 | 1.19e-05 | 3.0e-06 |
| Formaldehyd verzerrt | fc | -0.3208798889 | 2.9e-10 | 1.1e-09 | 9.5e-10 | 1.18e-05 | 2.9e-06 |
| Formaldehyd | voll | -0.3214440762 | -1.7e-09 | 8.0e-10 | 5.4e-10 | 1.01e-05 | 2.5e-06 |
| Formaldehyd | fc | -0.3162082526 | -1.7e-09 | 4.3e-10 | 9.0e-10 | 9.88e-06 | 2.5e-06 |
| Benzol | voll | -0.7969025894 | -2.5e-10 | 1.6e-09 | 1.2e-09 | 1.28e-04 | 1.1e-05 |
| Benzol | fc | -0.7793356368 | -2.5e-10 | 1.2e-09 | 1.2e-09 | 1.27e-04 | 1.1e-05 |
| Naphthalin | voll | -1.3240815899 | 6.3e-10 | - | - | 1.46e-04 | 8.1e-06 |
| Naphthalin | fc | -1.2945641517 | 6.8e-10 | - | - | 1.44e-04 | 8.0e-06 |
| Anthracen | voll | -1.8526760336 (`-qm.scf_threshold 1e-8`) | - | - | - | 1.61e-04 | 6.7e-06 |

- JFIT (H2O, Benzol, voll und fc): - PySCF <= 5.8e-10, - ORCA (`AuxC "def2/J"`) <= 1.8e-9, SCS - ORCA <= 1.6e-9.
- RI-JK-SCF + RI-MP2 (`-qm.ri jk`, Tensor wiederverwendet; H2O, Formaldehyd, Benzol, voll und fc): E(HF) - PySCF
  <= 4.6e-12, - ORCA (`! RIJK def2/JK`) <= 7.7e-10; E_corr - PySCF <= 1.7e-9, - ORCA <= 1.7e-9, SCS - ORCA <= 1.3e-9.
- `ctest -L qm_mp2`: 58/58 (davon 28 neue RI-MP2-Tests, 3 davon auf der GPU); groesste Abweichung zu PySCF in den Tests
  1.7e-9 Eh (Formaldehyd). GPU gegen CPU: E_corr auf allen 12 Stellen gleich (Benzol voll/fc, Formaldehyd fc;
  Naphthalin/Anthracen 1.3e-11 bzw. 1.1e-11 Unterschied, andere Summationsreihenfolge).
- Der RI-Fehler (3e-6 bis 3e-5 Eh je Atom) ist eine Messung, kein Testkriterium. Mit def2-SVP/C waere er nach der
  Literatur kleiner (aus dem Gedaechtnis, nicht gemessen; keine def2-SVP/C-Datei vorhanden).

**Zeiten** (Intel Core Ultra 9 285K, RTX 5090, MKL-sequential-BLAS (urspruenglich als netlib angegeben, siehe Korrektur), `-qm.scf_threshold 1e-8`, volle Korrelation; ms des
MP2-Teils: kanonisch Transformation / Energie, RI Setup / B-Transformation / Energie; Wand und max RSS fuer den ganzen
Lauf inkl. SCF; kanonisch CPU mit `-qm.scf_direct off` (gespeicherter Tensor), GPU mit dem Device-Tensor der SCF; RI
"exakt" = exakte SCF (Standard `scf_direct auto`, Naphthalin/Anthracen integral-direkt), "RI-JK" = `-qm.ri jk`):

| System (n, o, v, naux) | Pfad | 1 Thread | 16 Threads | GPU |
|---|---|---|---|---|
| Benzol (114, 21, 93, 558) | kanonisch | 271 / 6.1 ms; 4.3 s, 1.6 GB | 68 / 1.1 ms; 1.25 s, 1.8 GB | 63 / 0.3 ms; 0.63 s |
| | RI, SCF exakt | 236 / 58 / 74 ms; 4.3 s, 1.4 GB | 35 / 8.4 / 5.1 ms; 1.21 s | 35 / 9.3 / 55 ms; 0.66 s |
| | RI, SCF RI-JK | 0 / 26 / 74 ms; 1.31 s, 138 MB | 0 / 5.3 / 5.2 ms; 0.23 s | 0 / 6.8 / 3.9 ms; 0.62 s |
| Naphthalin (180, 34, 146, 894) | kanonisch | 2584 / 30 ms; 27.8 s, 9.8 GB | 561 / 6.9 ms; 8.9 s, 10.4 GB | 199 / 0.7 ms; 1.30 s |
| | RI, SCF exakt | 867 / 309 / 771 ms; 131 s, 372 MB | 114 / 33 / 51 ms; 9.9 s | 115 / 34 / 77 ms; 1.33 s |
| | RI, SCF RI-JK | 0 / 147 / 772 ms; 7.7 s, 369 MB | 0 / 18 / 51 ms; 1.08 s | 0 / 20 / 30 ms; 1.19 s |
| Anthracen (246, 47, 199, 1230) | kanonisch | nicht gerechnet | 3685 / 35 ms; 34.2 s, 31.3 GB | 1125 / 2.1 ms; 3.6 s |
| | RI, SCF exakt | 2002 / 1301 / 3828 ms; 371 s, 852 MB | 278 / 148 / 327 ms; 27.7 s | 275 / 142 / 135 ms; 3.0 s |
| | RI, SCF RI-JK | 0 / 518 / 3809 ms; 28.9 s, 839 MB | 0 / 54 / 288 ms; 3.6 s | 0 / 56 / 135 ms; 2.4 s |

- **Ehrlich:** der RI-Korrelationsteil ist auf der CPU (16 Threads) 3x (Naphthalin) bis 5x (Anthracen) schneller als
  die kanonische Transformation, der Speicher faellt von n^4 (31 GB Anthracen) auf unter 1 GB. Die Wandzeit bestimmt
  mit exakter SCF aber die SCF (integral-direkt: Naphthalin 1 Thread 131 s gegen 28 s kanonisch mit gespeichertem
  Tensor); erst mit `-qm.ri jk` wird der ganze Lauf schnell (Anthracen 3.6 s gegen 34 s). Auf der GPU ist die
  RI-Energie bei Benzol/Naphthalin mit exakter SCF langsamer als 16 CPU-Threads (55 bzw. 77 ms): der erste
  cuBLAS-Aufruf im Prozess kostet ~50 ms (mit RI-JK-SCF, cuBLAS schon geladen, 3.9 ms); bei Anthracen laeuft die GEMM
  mit ~1.7 TFLOPS (126 ms). Der kanonische GPU-Pfad bleibt bei diesen Groessen insgesamt gleich schnell, weil die
  Device-SCF den Tensor ohnehin haelt.
- Die kanonischen CPU-Zeiten hier (Benzol 68 ms bei 16 Threads) sind viel kleiner als die der WP11-Tabelle oben
  (913 ms); bestaetigt durch Nachmessung des Leiters (siehe Korrektur in der WP11-Zeitentabelle).
- **Nicht getestet:** Elemente ausser H, C, N, O, F (Li-Ne-Hilfsbasisdaten geladen, kein Molekuel gerechnet),
  geladene Molekuele, kartesische Basis (`-qm.cartesian_d true`), andere Orbitalbasen als def2-SVP, Systeme ueber
  246 BF, `-qm.ri j` + RI-MP2 nur auf H2O (JFIT, Tensor der SCF wiederverwendet: E(HF) und E_corr gegen PySCF mit
  DF-Coulomb + exaktem Austausch 4e-12 bzw. 8.0e-11 Eh, von Hand, nicht im ctest), Metrik `eigen` mit verworfenen Funktionen, `-hessian`/`-md` mit ri-mp2, ROCm/Vulkan (keine
  Implementierung), Speichergrenzen-Fehler (nur Codepfad).
- **Nicht implementiert gegenueber ORCA/PySCF:** def2-SVP/C (und jede andere -/C-Basis als Datei), RI-MP2-Gradient,
  offene Schalen, integral-direkte 3-Index-Integrale (Setup-Tensor naux n(n+1)/2 voll im Speicher), Batching der
  B-Transformation fuer sehr grosse Systeme, i <-> j-Symmetrie auf der GPU.

## WP12 -- MP3 bis MP5: gestrichen

Betreiberentscheidung Okt 2026: gemeint war MP2, und zwar die O(N^5)-Variante (kanonisches MP2 mit
4-Index-Transformation, danach RI-MP2, das formal ebenfalls N^5 skaliert, aber mit kleinerem Vorfaktor). MP3-MP5
sind nicht Teil des Plans; der fruehere Abschnitt (MP5 mit O(n^8), Pruefung ueber den vollen Determinantenraum)
steht in der Git-Historie dieser Datei.

## WP13 -- Closed-shell CCSD (RCCSD)

**Stand 18.10.2026: kanonisches RCCSD in der Engine (CPU), Methode `ccsd`** (🤖 AI-generated, ⚙️ machine-tested;
menschlicher Produktionstest steht aus). Code: `qm_methods/qm_ccsd.{h,cpp}` (Integrale, Residuum, Loeser,
`QMEngine::computeCCSD`), `ccsd_method.{h,cpp}` (Wrapper), Tests `test_cases/qm_ccsd/` (`ctest -L qm_ccsd`), Referenz
`scripts/ccsd_reference.py`. Vorlage war der Python-Prototyp `scripts/local_ccsd/ccsd_canonical.py`
([../LOCAL_CCSD_PROTOTYPE_2026-10.md](../LOCAL_CCSD_PROTOTYPE_2026-10.md)); jede Kontraktion traegt im C++ denselben
Index-String wie dort.

- **Verfahren:** spinadaptierte Amplitudengleichungen nach Hirata et al., JCP 120, 2581 (2004), Gl. 35-45 (Form von
  PySCF `rintermediates.py`), als Residuum R = N - e t; Jacobi-Schritt + Pulay-DIIS (8 Vektoren), Start t1 = 0,
  t2 = MP2, also Iteration 0 = MP2-Energie. Kanonische Fock-Bloecke aus den Orbitalenergien (f_ov = 0).
  Ladder als eine (o^2 x v^2)(v^2 x v^2)-GEMM, t1-Anteil der Ladder als o^3 v^3-Term (kein zweites v^4-Feld).
- **Integrale:** volle Vierfach-Transformation des AO-Tensors (Quelle wie MP2: gespeicherter SCF-Tensor, Download vom
  Device oder Neuaufbau bis `-qm.mp2_max_memory_mb`) in den aktiven MO-Raum, dann Bloecke oooo, ovoo, ovov, oovv, ovvo,
  ovvv, vvvv. Speicher n^3 m + n^2 m^2 waehrend der Transformation, danach v^4 + o v^3 + Amplituden; Schaetzung gegen
  `-qm.cc_max_memory_mb` (32000) vor dem Start, darueber Fehler mit der benoetigten Groesse.
- **Kontraktionen:** `qmint::ccsdContract("kcld,ilcd->ki", ...)`: Permutation (nur wenn noetig, gekachelt) + eine
  BLAS-GEMM mit Transpositions-Flags; OpenMP ueber einen Block je Thread (MKL sequentiell, siehe "BLAS (Korrektur)").
- **Schalter** (Modul `qm`): `-qm.cc_max_iter` (100), `-qm.cc_threshold` (1e-10, |dE|), `-qm.cc_residual_threshold`
  (1e-8, max|R|), `-qm.cc_diis` (8), `-qm.cc_max_memory_mb`; `-qm.frozen_core` wie bei MP2. Nicht konvergiert = Fehler.
- **Ausgabe:** E(HF), E(MP2 corr) (Iteration 0), E(CCSD corr), Iterationen mit letztem dE und max|R|, T1-Diagnostik
  (sqrt(sum t_ia^2 / 2o), dieselbe Konvention wie PySCF und ORCA, an H2O geprueft), Gesamtenergie.
- **Fehler statt Ergebnis:** Gradient (`-opt`, `-gradient`: "ccsd has no analytic gradient"), `-qm.ri j|jk` (ein
  RI-HF-Referenzzustand mit exakten CC-Integralen waere eine Mischmethode), ungerade Elektronenzahl, Speicherlimit.
  `-gpu cuda`: keine Device-Implementierung, gezaehlter Fallback `qm_ccsd` (SCF kann auf dem Device laufen).

**Validierung** (def2-SVP, `-qm.scf_threshold 1e-10`, je ein Lauf; PySCF 2.14 RHF conv_tol 1e-12 + `cc.CCSD`
conv_tol 1e-12; ORCA 6.1 `! RHF CCSD def2-SVP ExtremeSCF [No]FrozenCore`, `%mdci STol 1e-10`, Werte aus
`job.property.txt`). curcuma - Referenz fuer E_corr(CCSD), Eh:

| Molekuel | n | o / v (voll) | E_corr voll | - PySCF | - ORCA | E_corr fc | - PySCF | - ORCA | it |
|---|---|---|---|---|---|---|---|---|---|
| H2O | 24 | 5 / 19 | -0.209360921262 | 2.2e-10 | 5.0e-10 | -0.207022884026 | 2.0e-10 | 4.8e-10 | 12 |
| NH3 | 29 | 5 / 24 | -0.209721196637 | -3.3e-10 | 3.2e-11 | -0.207197008337 | -3.5e-10 | 3.0e-11 | 16 |
| CH4 | 34 | 5 / 29 | -0.187388088172 | -2.4e-11 | 1.8e-10 | -0.184068822300 | -2.4e-11 | 2.2e-10 | 13 |
| Formaldehyd | 38 | 8 / 30 | -0.336554568128 | 4.4e-11 | 7.9e-10 | -0.331728862546 | 4.3e-11 | 8.2e-10 | 16 |
| Ethan | 58 | 9 / 49 | -0.344209177048 | 2.7e-11 | 3.7e-10 | -0.338032497971 | 2.7e-11 | 4.4e-10 | 15 |
| Wasserdimer | 48 | 10 / 38 | -0.427867429254 | -2.7e-11 | 5.9e-10 | -0.423398893268 | -8.8e-11 | 5.5e-10 | 16 / 15 |
| Benzol | 114 | 21 / 93 | -0.836086277674 | 5.5e-11 | 1.3e-09 | -0.819682777012 | 5.6e-11 | 1.5e-09 | 15 |

- Maximal 3.5e-10 Eh gegen PySCF, 1.5e-9 Eh gegen ORCA (Ziel 1e-8); E(HF) gegen PySCF <= 7.3e-12 Eh, Iteration 0
  gegen PySCF-MP2 <= 2.1e-10 Eh; T1-Diagnostik gegen PySCF <= 5e-9 (Druckgenauigkeit 1e-8). Der Rest gegen PySCF
  (2-3e-10 bei H2O/NH3) folgt aus `-qm.cc_threshold 1e-10` (Prototyp mit 1e-11: 8e-12); der Rest gegen ORCA folgt
  ORCAs RHF (E(HF) curcuma - ORCA bis 1.2e-9).
- Weitere Tests: integral-direktes SCF (Tensor fuer CCSD neu aufgebaut), 1 Thread, die drei Fehlerpfade. ORCA lief nur
  hier; eine zweite Maschine ist fuer Zeiten nicht noetig (Betreibervorgabe 2026-10-26).
- **Zeiten** (CPU, `-qm.threads 16`, CCSD-Teil ohne SCF, Transformation / Iterationen; je ein Lauf, **unter Last**: ein
  paralleler Job belegte etwa 16 der 24 Kerne, die Zahlen sind daher obere Schranken und die Thread-Skalierung ist
  nicht aussagekraeftig):

  | System (n, o, v) | voll | fc |
  |---|---|---|
  | H2O (24, 5, 19) | 5 / 23 ms | 4 / 13 ms |
  | NH3 (29, 5, 24) | 10 / 48 ms | 7 / 29 ms |
  | CH4 (34, 5, 29) | 16 / 75 ms | 16 / 46 ms |
  | Formaldehyd (38, 8, 30) | 24 / 271 ms | 21 / 140 ms |
  | Wasserdimer (48, 10, 38) | 65 / 939 ms | 61 / 569 ms |
  | Ethan (58, 9, 49) | 131 / 1387 ms | 107 / 814 ms |
  | Benzol (114, 21, 93) | 1.6 / 40.6 s | 1.3 / 22.3 s |

  Benzol fc mit 1 / 4 / 16 Threads: Iterationen 42.3 / 24.8 / 22.3 s (unter derselben Last). Schaetzung (nicht
  gemessen): die Iteration hat bei Benzol fc ~1e11 Flops; ein unbelasteter Lauf mit 16 Threads sollte deutlich unter
  der Messung liegen. Speicher Benzol voll geschaetzt 3.6 GB.
- **Nicht getestet:** offene Schale (nicht implementiert), Elemente ausser H, C, N, O, andere Basen als def2-SVP,
  kartesisch 6d, Systeme ueber 114 Basisfunktionen, geladene Molekuele, Thread-Skalierung auf freier Maschine,
  schlecht konvergierende Faelle (grosse T1, fast entartete HOMO/LUMO).
- **Nicht implementiert gegenueber PySCF/ORCA:** (T), Gradient, Dichten/Eigenschaften, RI- oder Cholesky-Integrale,
  AO-getriebene Ladder (kein v^4 im Speicher), Device-Pfad, Symmetrie, offene Schalen, EOM, lokale Varianten.
- **Folgearbeit:** (T)-Korrektur (Standard, Ziel der DLPNO-CCSD(T)-Referenz im Vault); Ladder auf dem Device.

## WP15 -- Lokales MP2 mit PNO-Truncation in ueberlappenden Domaenen (Idee des Betreibers, Okt 2026)

**Stand 18.10.2026: Python-Prototyp gemessen (17.10.), C++-Methode `lmp2` in der Engine (Abschnitt "Engine-Implementierung" unten)** (🤖 AI-generated, ⚙️ machine-tested;
[../LOCAL_MP2_FEASIBILITY_2026-10.md](../LOCAL_MP2_FEASIBILITY_2026-10.md), Skripte `scripts/local_mp2/`). Alkane bis
C20H42, (H2O)4-16, Benzol bis Pentacen, def2-SVP, RI mit JKFIT, Referenz volles RI-MP2:
- Paarenergien fallen mit R^-6 (Acene R^-5.7), Zahl signifikanter Paare und PNO-Gesamtzahl linear; BP(0.02)-Domaenen
  verlieren 1-2.5 %, BP(0.02) + 7 Bohr und T_CutPNO 1e-8 geben 3e-6 (Wasser) bis 3e-5 Eh/Atom (Pentacen), ORCA
  NormalPNO 2-3e-5 Eh/Atom.
- Schwarz-Vorkonditionierer (Stufe 5): halbiert die Iterationen (7-9 -> 3-5), kostet aber das 3-30fache; kein Gewinn.
- Empfehlung des Dokuments: C++-LMP2 erst, wenn MP2 fuer > 100 Atome Ziel ist und die SCF dort skaliert.

### Engine-Implementierung: Methode `lmp2` (18.10.2026)

**Stand: lokales PNO-RI-MP2 in C++, closed shell, nur Energie. 🤖 AI-generated, ⚙️ machine-tested; menschlicher
Produktionstest steht aus.** Code: `qm_lmp2.{h,cpp}` (Verfahren, Gleichungen, Literatur im Kopf), `lmp2_method.{h,cpp}`,
`QMEngine::correlationRI` (aus `qm_rimp2.cpp` herausgezogen, RI-MP2 und LMP2 teilen den RI-Tensor). Tests
`test_cases/qm_lmp2/` (`ctest -L qm_lmp2`), Referenzskript `scripts/lmp2_reference.py`.

- **Aufruf:** `-method lmp2` (Alias `local-ri-mp2`), beliebige SCF (exakt oder `-qm.ri j|jk`), `-qm.frozen_core`,
  Hilfsbasis `-qm.ri_mp2_aux` (Standard def2-universal-JKFIT, wie RI-MP2). Schalter (Registry, Modul `qm`):
  `lmp2_localizer pm|boys`, `lmp2_bp_threshold` 0.02, `lmp2_domain_radius` 7 Bohr, `lmp2_pno_cut` 1e-8,
  `lmp2_pno_correction` true, `lmp2_pair_cut` 15 Bohr, `lmp2_distant_pairs dipole|drop`, `lmp2_fock_cut` 1e-5,
  `lmp2_threshold` 1e-9, `lmp2_residual_threshold` 1e-7, `lmp2_max_iter`, `lmp2_diis`, `lmp2_max_memory_mb`,
  `lmp2_check` (rechnet zusaetzlich kanonisches RI-MP2 mit derselben Hilfsbasis und druckt die Differenz).
- **Verfahren wie der Prototyp** (`lmp2_common.py`): PM (Mulliken-Populationen, Jacobi-Sweeps ab Cholesky-Orbitalen,
  disjunkte Paare je Runde, daher threadzahlunabhaengig) oder Boys; BP(0.02) + Radius; PAO-Paarraeume (Eigenwertschnitt
  1e-6); `B_{i,r}^P` (lokales i, PAO r) per Halbtransformation des RI-Tensors (`fitRows` wie RI-MP2), nur fuer die PAOs
  der erweiterten Domaene von i; K_ij in der PAO-Basis, PNOs aus semikanonischer Paardichte; Residual mit
  PNO-Ueberlappungen, diagonaler Vorkonditionierer + DIIS; kein o^2 v^2-Tensor.
- **Korrekturen:** (a) wie der Prototyp die semikanonische PNO-Abschneidekorrektur je Paar; (b) **neu gegenueber dem
  Prototyp** (dort wurden alle Paare gehalten): Paare mit Schwerpunktabstand > `lmp2_pair_cut` bekommen die
  Dipol-Dipol-Abschaetzung (Pinski 2015) ueber die PAO-Raeume von [i] und [j], Austausch vernachlaessigt.
- **Abweichung von der Vorgabe:** PM mit Mulliken- statt Meta-Loewdin-Populationen (PySCF-Standard im Prototyp);
  Folge: andere LMOs, Wasserfehler 1.5x hoeher als im Prototyp (unten).
- **DIIS-Falle (gefunden und behoben):** die unskalierte DIIS-Matrix (Fehlerprodukte ~1e-20 neben dem Rand -1) liess
  die Iteration bei max|R| ~1e-9 stehen; Fehlerblock jetzt auf seine groesste Diagonale normiert (14-16 statt > 50 Iterationen).

**Exaktheit (alle Paare, vollstaendige PAO-Raeume, T_CutPNO 0, alle F_ik):** LMP2 = curcumas RI-MP2 (gleiche Hilfsbasis)
auf <= 2e-12 Eh: Butan fc/voll, Wasserdimer (Boys), Benzol fc (PM und Boys), Wasserdimer mit RI-JK-SCF und mit JFIT,
Essigsaeuredimer fc mit 1 Thread. Damit sind Lokalisierung, PAO-Raeume, B-Faktoren und Residual geprueft.

**Standardschwellen gegen kanonisches RI-MP2** (curcuma, gleiche Hilfsbasis, RI-JK-SCF, def2-SVP, fc; Fehler =
E(LMP2) - E(RI-MP2) in Eh/Atom; "alle Paare" = `-qm.lmp2_pair_cut 1e6 -qm.lmp2_fock_cut 0`, das Setting des Prototyps
Q3 "BP(0.02)+r7, 1e-8"; je ein Lauf):

| System | Atome | Standard | alle Paare | Prototyp | PNOs/Paar | Paare gehalten |
|---|---|---|---|---|---|---|
| C4H10 | 14 | 5.96e-6 | 5.96e-6 | - | 42.5 | 91/91 |
| C8H18 | 26 | 1.21e-5 | 1.18e-5 | 1.2e-5 | 35.2 | 306/325 |
| C12H26 | 38 | 1.65e-5 | 1.59e-5 | 1.6e-5 | 33.0 | 538/703 |
| C16H34 | 50 | 1.90e-5 | 1.81e-5 | 1.8e-5 | 32.1 | 770/1225 |
| C20H42 | 62 | 2.05e-5 | 1.95e-5 | 1.9e-5 | 31.7 | 1002/1891 |
| C30H62 | 92 | 2.26e-5 | 2.15e-5 | - | 31.1 | 1582/4186 |
| C40H82 | 122 | 2.36e-5 | 2.24e-5 | - | 30.8 | 2162/7381 |
| Benzol | 12 | 1.57e-5 | 1.57e-5 | - | 46.3 | 120/120 |
| Naphthalin | 18 | 2.11e-5 | 2.11e-5 | 2.1e-5 | 46.9 | 300/300 |
| Anthracen | 24 | 2.59e-5 | 2.58e-5 | 2.4e-5 | 43.9 | 553/561 |
| Tetracen | 30 | 5.26e-5 | 2.90e-5 | 2.9e-5 | 43.1 | 812/903 |
| Pentacen | 36 | 8.60e-5 | 3.15e-5 | 3.1e-5 | 42.9 | 1069/1326 |
| (H2O)4 | 12 | 4.51e-6 | 4.51e-6 | - | 15.4 | 136/136 |
| (H2O)8 | 24 | 3.68e-6 | 3.68e-6 | 2.5e-6 | 11.1 | 527/528 |
| (H2O)12 | 36 | 4.83e-6 | 4.78e-6 | 2.9e-6 | 9.8 | 1035/1176 |
| (H2O)16 | 48 | 5.03e-6 | 4.97e-6 | 3.1e-6 | 8.9 | 1556/2080 |
| (H2O)20 | 60 | 5.78e-6 | 5.69e-6 | - | 8.4 | 2028/3240 |

- Alkane und Acene treffen den Prototyp (<= 8 %), Wasser liegt 1.5x darueber (Mulliken-PM, s.o.). Fehler/Atom der
  Alkane waechst noch langsam (C40 2.4e-5), wie im Prototyp.
- `lmp2_fock_cut` 1e-5 aendert nichts (<= 4e-10 Eh bei C20H42, (H2O)20, Pentacen gegen 0).
- Paarschnitt 15 Bohr: bei Alkanen und Wasser +0.1e-5 bzw. +0.1e-6 Eh/Atom; die Dipolabschaetzung gewinnt 61-63 % der
  fehlenden Paarenergie zurueck (C20, C40), bei Wasser mehr als die Differenz (Fehler 5.8e-6 statt 6.5e-6 mit `drop`).
  **Acene: 15 Bohr reicht nicht** (Pentacen 8.6e-5, Dipol holt 6 % der fehlenden 2e-3 Eh; die pi-Paare bei 15 Bohr sind
  nicht dipolartig); `-qm.lmp2_pair_cut 20` gibt 3.9e-5, `25` = alle Paare (3.15e-5). Fuer pi-Systeme 25 Bohr setzen.

**Zeit** (24 Threads, Korrelationsteil ohne SCF; RI-MP2 = `lmp2_check` im selben Prozess, gleicher RI-Tensor; je ein
Lauf, Maschine sonst ruhig):

| System | Basis | SCF (RI-JK) s | LMP2 s | RI-MP2 s | LMP2-Speicher MB |
|---|---|---|---|---|---|
| C8H18 | 202 | 0.40 | 0.25 | 0.05 | 193 |
| C12H26 | 298 | 1.5 | 0.66 | 0.36 | 399 |
| C16H34 | 394 | 4.1 | 1.24 | 2.3 | 640 |
| C20H42 | 490 | 9.3 | 2.07 | 13.2 | 917 |
| C30H62 | 730 | 42 | 5.8 | 153 | 1765 |
| C40H82 | 970 | 150 | 12.4 | 696 | 2834 |
| C60H122 | 1450 | 701 | 36.4 | nicht gerechnet | 5638 |
| C80H162 | 1930 | 2715 | 92.9 | nicht gerechnet | 9329 |
| Pentacen | 378 | 5.0 | 3.9 | 2.3 | 1190 |
| (H2O)12 | 288 | 1.6 | 0.58 | 0.51 | 219 |
| (H2O)20 | 480 | 11.4 | 2.0 | 18.0 | 674 |

- **Break-even** gegen curcumas RI-MP2: Alkane zwischen C12 und C16 (300-400 Basisfunktionen), Wasser bei (H2O)12,
  Acene bis Pentacen nicht erreicht. Die SCF dominiert jede Rechnung (C60: 701 s SCF, 36 s LMP2; C80: 2715 s / 93 s).
- Groesstes gerechnetes System: C80H162 (242 Atome, 1930 Funktionen, 8916 Hilfsfunktionen; 4482 von 29161 Paaren
  iteriert, 30.5 PNOs/Paar, 11 Iterationen; RI-JK-Tensor ~133 GB). Kanonisches RI-MP2 wurde ab C60 nicht gerechnet
  (C40 schon 696 s, Skalierung N^5). E_corr(C60) bzw. E_corr(C80) weichen von der linearen Extrapolation aus C20-C40
  bzw. C30-C60 um 9e-9 bzw. 3e-9 Eh ab (Groessenkonsistenz, kein Referenzwert); der Fehler gegen RI-MP2 ist dort unbekannt.
- **Grenzen:** die SCF bzw. der gespeicherte RI-JK-Tensor (naux x n(n+1)/2: C60 56 GB, C100 ~260 GB > RAM) begrenzt,
  nicht LMP2. Im LMP2 waechst der Paarschritt (K_ij in der PAO-Basis ueber alle naux Hilfsfunktionen) wie N^2 und die
  B-Halbtransformation wie N^4 (Paarschritt C40/C60/C80: 6.6/15/28 s, B: 4.0/17/52 s); lineare Skalierung braeuchte lokale Fit-Domaenen (WP14-Ergebnis:
  robust lokal, R >= 8-10 Bohr) und eine integraldirekte 3-Index-Erzeugung.

**ORCA 6.1 DLPNO-MP2** (NormalPNO, `RIJK def2/JK`, `%basis AuxC "def2/JK"`, fc; nur Plausibilitaet, andere Domaenen und
Lokalisierung): kanonisches RI-MP2 ORCA vs curcuma 2e-9 (C12H26), 8.5e-8 ((H2O)16), 2.5e-7 Eh (Anthracen; SCF-Schwelle
1e-8 in curcuma); Fehler/Atom curcuma LMP2 / ORCA DLPNO: C12H26 1.65e-5 / 2.55e-5, (H2O)16 5.0e-6 / 1.39e-5, Anthracen
2.59e-5 / 3.27e-5; PNOs/Paar 33.0 / 26.3, 8.9 / 9.3, 43.9 / 49.3.

- **Getestet:** H, C, O; def2-SVP; Alkane C4-C80 (Fehler nur bis C40), Acene 1-5 Ringe, (H2O)4-20, Wasser-/Essigsaeuredimer; fc und voll;
  exakte und RI-JK-SCF; JKFIT und JFIT; 1 und 24 Threads (Exaktheitstest).
- **Nicht getestet:** andere Elemente als H, C, O; andere Basen (TZVP, diffuse:
  groessere Domaenen); geladene Systeme; kleine HOMO-LUMO-Luecken jenseits Pentacen; 3D-kompakte Systeme > 60 Atome;
  Thread-Reproduzierbarkeit bei Standardschwellen nur Butan (ctest); Laufzeiten auf anderer Hardware.
- **Nicht implementiert gegenueber DLPNO-MP2 (ORCA):** PNO-Domaenen (TCutDO), Vorselektion ueber semikanonische
  Paarenergien (nur Abstand), offene Schalen, Gradient, SCS, Device-Pfad. Lokale Fit-Domaenen und integraldirekte
  3-Index-Integrale: seit 18.10.2026 (`-qm.lmp2_fit local`, naechster Abschnitt).

### Lokales Fitting und integraldirekte 3-Index-Integrale (`-qm.lmp2_fit local`, 18.10.2026)

**Stand: 🤖 AI-generated, ⚙️ machine-tested; menschlicher Produktionstest steht aus.** Ziel: kein globaler
3-Index-Tensor (Speicher ~N^3, B-Transformation ~N^4) mehr im LMP2. Code: `qm_lmp2.cpp` Abschnitt 3b
(`localThreeIndex`, `FitFactor`), `qmri::ThreeIndexDirect` in `qm_ri.{h,cpp}` (Abschnitt 4b: dieselben Helfer
`makePair`/`pairBlock`/`transformAxis3`/OutMaps und derselbe Schwarz-Ausdruck wie `RIJK::setup`, das selbst unveraendert ist).

- **Verfahren:** Fit-Domaene F_i = Hilfsfunktionen auf allen Atomen im Abstand <= `lmp2_fit_radius` (Bohr) von einem
  BP-Atom des Orbitals i. Pro Orbital werden die rohen Integrale X_i(P, r) = (P|i r) fuer P in G_i = Vereinigung der
  F_j aller behaltenen Partner j und r = PAOs der erweiterten Domaene integraldirekt erzeugt: je Hilfsschale c auf
  Atom K fuer alle Orbitale mit K in G_i die Bloecke (c|AB) ueber die AO-Traeger der Orbitale (|L| > `lmp2_ao_cut`)
  und ihrer PAOs, dann L^T (eine GEMM fuer alle diese Orbitale) und C_PAO; keine (P|mn) werden gespeichert.
  Paarintegrale: `lmp2_fit_domain orbital` (Standard) robust nach Dunlap mit den Orbitaldomaenen,
  (ia|jb) = d_ia^T (P|jb)_{F_i} + (ia|Q)_{F_j} d_jb - d_ia^T J_{F_i F_j} d_jb, d_i = J_{F_i}^-1 X_i(F_i); umgesetzt als
  d_i^T Z_j(F_i) + X_i(F_j)^T d_j mit Z_j = X_j - J(G_j, F_j) d_j (einmal pro Orbital). `pair`: symmetrischer Fit in
  F_ij = F_i u F_j (Cholesky je verschiedener Domaene). Metrik: Cholesky des Coulomb-normierten Teilblocks
  (eigen mit `-qm.ri_metric eigen` oder wenn Cholesky scheitert). Die 2-Index-Metrik wird voll gespeichert (naux^2).
- **Abweichung von der Vorgabe:** Standard ist die robuste Orbitaldomaene statt F_ij: gleicher Fehler ab 6 Bohr (Tabelle),
  eine Faktorisierung pro Orbital statt pro Paardomaene, Paarschritt 1.5x (C30-C80) bis 2.7x (PEO8) schneller. Die
  Integralmenge haengt an G_i (Partnerdomaenen), nicht an F_i, und ist fuer beide Varianten gleich.
- **SCF:** das lokale LMP2 ruft `QMEngine::correlationRI` nicht auf; als Referenz geht jede SCF (exakt
  integraldirekt `-qm.ri off`, auch `-gpu cuda -qm.gpu_direct on`, oder RI-JK).

**(1) Exaktheit** (unendlicher Fit-Radius 1e6, `lmp2_ao_cut 0`, vollstaendige Raeume, alle Paare, T_CutPNO 0, alle
F_ik) gegen curcumas RI-MP2 (gleiche Hilfsbasis): Butan fc (pair und orbital) und Wasserdimer (Boys; RI-JK-SCF)
<= 1e-13, Essigsaeuredimer fc 6e-14, Benzol fc mit JFIT 2.6e-12 Eh; die Integralzahl ist dort gleich der von
`RIJK::setup` (Butan 235880 Bloecke).

**(2) Fehler des lokalen Fittings allein** (alles andere auf Standard; E(lokal) - E(global) in Eh/Atom, RI-JK-SCF,
def2-SVP, fc, je ein Lauf; globale Referenz = die Laeufe der Tabelle oben, mit dem aktuellen Binary fuer C20 bitgleich
reproduziert). `lmp2_ao_cut 0` (reiner Fit-Fehler), `orbital`; `pair` stimmt ab R = 6 Bohr auf 2 Stellen ueberein
(R = 4: bis 14 % anders):

| System | R = 4 | R = 6 | R = 8 | R = 10 | R = 12 |
|---|---|---|---|---|---|
| (H2O)4 | 4.3e-7 | 7.7e-8 | -1.7e-8 | 2.9e-9 | 1.8e-9 |
| (H2O)8 | 4.9e-7 | 1.4e-7 | 2.3e-8 | 4.9e-9 | 3.0e-9 |
| (H2O)12 | 8.4e-7 | 2.3e-7 | 3.2e-8 | 1.0e-8 | 2.8e-9 |
| (H2O)16 | 8.0e-7 | 2.2e-7 | 4.0e-8 | 1.2e-8 | 4.5e-9 |
| (H2O)20 | 7.4e-7 | 1.9e-7 | 3.4e-8 | 1.5e-8 | 7.1e-9 |
| Benzol | 6.7e-6 | 3.9e-7 | 0 | 0 | 0 |
| Naphthalin | 8.1e-6 | 1.3e-6 | 2.1e-7 | 1.5e-8 | 0 |
| Anthracen | 9.2e-6 | 1.7e-6 | 3.9e-7 | 1.4e-7 | 5.0e-8 |
| Tetracen | 9.5e-6 | 1.9e-6 | 4.8e-7 | 2.0e-7 | 8.8e-8 |
| Pentacen | 9.9e-6 | 2.0e-6 | 5.3e-7 | 2.3e-7 | 1.1e-7 |
| C8H18 | 7.3e-6 | 7.4e-7 | 2.5e-7 | 8.3e-8 | 4.0e-8 |
| C12H26 | 7.8e-6 | 8.7e-7 | 3.3e-7 | 1.3e-7 | 7.0e-8 |
| C20H42 | 8.3e-6 | 9.6e-7 | 4.0e-7 | 1.7e-7 | 9.4e-8 |
| C30H62 | 8.5e-6 | 1.0e-6 | 4.3e-7 | 1.8e-7 | 1.1e-7 |
| C40H82 | 8.6e-6 | 1.0e-6 | 4.4e-7 | 1.9e-7 | 1.1e-7 |
| PEO8 (C16H34O9, `test_cases/structures/medium/peo8-chain.h-capped`, `lmp2_ao_cut` 1e-6) | - | - | 9.9e-7 | 3.1e-7 | 8.1e-8 |

- "0" = Fit-Domaene umfasst schon alle Hilfsfunktionen. Ziel <= 1e-6 Eh/Atom: **R = 8 Bohr (Standard) haelt es in allen
  Systemen, PEO8 nur knapp (9.9e-7)**; R = 6 verfehlt es bei Acenen (bis 2e-6) und knapp bei langen Alkanen. Der Fehler
  ist positiv (Korrelation fehlt) und pro Atom in Alkanen ab C30 konstant (Groessenkonsistenz des Fehlers).
- `lmp2_ao_cut` (bei R = 8, `pair`; Fehler gegen global): 1e-7 und 1e-6 gleich dem reinen Fit-Fehler (<= 1e-9 Eh/Atom
  Unterschied), 1e-5 aendert ihn um bis -2e-7 (Alkane, teils Kompensation), 1e-4 gibt -3e-5..-5e-5 (unbrauchbar).
  Die Laufzeit haengt kaum am Schnitt (C30: 2.80 / 2.75 s fuer 1e-6 / 1e-5), weil die Schalenvereinigung ueber alle
  Orbitale mit K in G_i die Arbeit bestimmt. Standard 1e-6.

**(3) Zeit und Speicher** (24 Threads, Maschine sonst ruhig, je ein Lauf; "lokal" = Standard-Lokalfit, orbital, R = 8,
ao_cut 1e-6; "Integrale" = integraldirekte Erzeugung von X_i inkl. beider Transformationen bzw. global die
B-Halbtransformation; Spitzen-RSS des Prozesses aus /proc VmHWM, schliesst die SCF ein):

| System | Basis / aux | SCF | SCF s | LMP2 global s | LMP2 lokal s | Integrale global (B) / lokal s | Paarschritt lokal s | Spitzen-RSS global / lokal GB |
|---|---|---|---|---|---|---|---|---|
| C20H42 | 490 / 2256 | RI-JK | 9.3 | 2.09 | 2.62 | 0.24 / 0.95 | 1.17 | 3.5 / 3.8 |
| C30H62 | 730 / 3366 | RI-JK | 43 | 5.82 | 5.88 | 1.41 / 2.83 | 2.10 | 9.1 / 9.2 |
| C40H82 | 970 / 4476 | RI-JK | 150 | 12.45 | 10.32 | 3.97 / 5.43 | 3.06 | 19.9 / 19.1 |
| C60H122 | 1450 / 6696 | GPU direkt (lokal), RI-JK (global, alte Messung) | 472 / 701 | 36.4 | 21.24 | 17 / 11.65 | 5.12 | ~56 GB Tensor / 5.6 |
| C80H162 | 1930 / 8916 | GPU direkt (lokal), RI-JK (global, alte Messung) | 1521 / 2715 | 92.9 | 38.53 | 52.3 / 18.10 | 7.66 | ~133 GB Tensor / 7.1 |
| C80H162, lokal mit RI-JK-SCF (8 Threads, nur Fehler) | 1930 / 8916 | RI-JK | 2991 | - | - | - | - | 132.7 (SCF-Tensor) |

- Global C60/C80: Zeiten aus der Tabelle oben (gleiches Binary bis auf die Zeitmessung des Paarschritts; die LMP2-Gesamtzeit
  ist davon nicht betroffen). C20-C40 global mit dem aktuellen Binary neu gemessen. "GPU direkt" = exakte integraldirekte
  SCF auf der RTX 5090 (`-qm.ri off -gpu cuda -qm.gpu_direct on -qm.gpu_scf on`, SCF-Code unveraendert); sie braucht kein
  RI und haelt den Prozess bei 5.6 / 7.1 GB, ist aber nicht billiger: C40 194 s (RI-JK CPU 150 s), C60 472 s, C80 1521 s
  (42 s pro J/K-Aufbau, ~N^2.7).
- **Skalierungsexponenten** (log-log-Steigung zwischen zwei Messpunkten ueber die Zahl der Basisfunktionen, keine
  Ausgleichsrechnung): lokale Integrale C40->C80 1.74, C60->C80 1.54 (globale B-Transformation C40->C80 3.7); lokaler
  Paarschritt C40->C80 1.32, C60->C80 1.40; LMP2 lokal gesamt C60->C80 2.1, davon die Lokalisierung (unveraendert,
  1.6 -> 7.8 s) der am schnellsten wachsende Teil; LMP2 ohne Lokalisierung C60->C80 1.58. Linear ist die Integralstufe
  damit noch nicht: G_i (Vereinigung der Partner-Fitdomaenen, ~23 Bohr um das Orbital) und die AO-Traeger saettigen
  erst ab etwa C60 (mittlere X-Zeilen 1636 / 1889 / 1974 fuer C20 / C40 / C60).
- Fit-Fehler bei C80 (RI-JK-SCF, lokal gegen global, alles andere Standard): 4.6e-7 Eh/Atom (C40: 4.4e-7).
- Speicher LMP2 selbst (Ausgabe "memory"): lokal C80 6.9 GB (X_i, d_i, Z_i, Metrik 0.6 GB, PNO-Daten), global 9.3 GB
  B-Faktoren zusaetzlich zum 133-GB-Tensor.

**Groesstes gerechnetes System:** C80H162 (242 Atome, 1930 Basisfunktionen, 8916 Hilfsfunktionen) mit GPU-SCF: SCF
1521 s, LMP2 38.5 s, Prozess-Spitze 7.1 GB. **C120H242 (362 Atome) nicht fertig:** die exakte GPU-direkt-SCF lief 3.9 h
ohne Ergebnis und wurde abgebrochen (Ursache nicht untersucht; erwartet waren ~1.2 h aus der C40-C80-Steigung); RI-JK
ist dort aus Speichergruenden (Tensor ~450 GB) nicht moeglich. C160 deshalb nicht versucht. Die SCF, nicht LMP2, ist die
Grenze.
**Nachtrag 19.10.2026 -- die SCF-Grenze ist verschoben:** die exakte GPU-direkt-SCF braucht nach der Ueberarbeitung
(sortierte Paarlisten statt O(npair^2)-Rechteck, zweistufiger Hermite-Kern fuer die grossen Klassen, Atom-HF-Startdichte
`-qm.scf_guess atom`, lockere Schwelle in der Anfangsphase, dichtegewichtetes Primitiv-Screening; alles energieneutral auf
10 Stellen) C40 31.6 s, C60 72 s, C80 130 s (vorher 1381 s) und C120H242 (2890 Funktionen) 306 s, je 13 Iterationen,
Wandzeit ~n^2.1. Messungen und Grenzen: [QM_GPU_ROADMAP.md, Abschnitt 2 "Large systems"](../QM_GPU_ROADMAP.md). LMP2 auf
C120 wurde in dieser Sitzung nicht gerechnet. Der alte Code oszilliert auf C120 (||dP|| 4e2..6e3 ueber 20 Iterationen, ~170 s pro Fock-Aufbau,
nach ~1 h abgebrochen): Ursache der fruehere 3.9-h-Rechnung war die alte Startdichte zusammen mit den langsamen Aufbauten.

- **Getestet:** H, C, O; def2-SVP, JKFIT (JFIT nur Exaktheit); Systeme der Tabelle (2) plus PEO8; RI-JK- und exakte
  (GPU-direkt) SCF; 1 und 4 Threads (ctest, bitgleich auf 12 Stellen); ctest `qm_lmp2_*` (24 Tests).
- **Nicht getestet:** andere Elemente als H, C, O; andere Basen; geladene und 3D-kompakte Systeme > 60 Atome; der
  Fit-Fehler ueber C40 nur fuer C80 (RI-JK-SCF) gegen global; das Thread-Verhalten nur im ctest (1 vs 4 Threads); die
  `pair`-Variante ueber C80 hinaus.
- **Nicht implementiert:** lineare Skalierung der 2-Index-Metrik (voll naux^2 gespeichert, Aufbau N^2; C160 ~2.6 GB),
  Lokalisierung (Jacobi-Sweeps ~N^3) und Paar-Vorauswahl (o^2 Abstandstest) bleiben wie vorher; keine PNO-Domaenen
  fuer den Fit (ORCA fittet im PNO-Raum mit Paardomaenen); keine Device-Variante der 3-Index-Erzeugung.

### Skalierung der Lokalisierung (`-qm.lmp2_localizer_algo fast|classic`, 19.10.2026)

**Stand: 🤖 AI-generated, ⚙️ machine-tested; menschlicher Produktionstest steht aus.** Code: `qm_lmp2.cpp` Abschnitt 2
(`pmLocalizeFast`, `productParallel`), Schalter `lmp2_localizer_algo` (Standard `fast`, `classic` = alter Weg) und
`lmp2_localizer_screen` (1e-8). Zeitaufteilung bei `-verbosity 2` in der Zeile "LMP2 localisation".

- **Profil vorher** (C40H82, n = 970, 24 Threads, ms): sphaerische Transformation der Dipolmatrizen (dichtes Q^T M Q,
  n^3) 214, PM-Sweeps 120, Schwerpunkte (o Matrix-Vektor-Produkte mit n x n, n^2 o) 54, Populationsmatrizen 28,
  Dipolintegrale 23, Cholesky-Start 17, U 13. Die Sweeps drehen natm Matrizen o x o (O(natm o^3) je Sweep, 111 MB bei
  C80, Zeilenzugriffe mit Stride); die Zahl der Sweeps waechst nicht (Alkane 6, Wasser 5, Acene 6-27).
- **Aenderungen:** (a) PM-Sweeps ohne Populationsmatrizen: Q^A_st je Paar aus den Spalten s, t von L und SL = S L, nach
  einer Drehung nur die zwei Orbitale neu bewertet; gleiche Zielfunktion, gleiche Paarformel, gleicher Round-Robin-Plan,
  gleiches Konvergenzkriterium; O(o^2 n) je Sweep, parallel ueber die disjunkten Paare einer Runde (threadzahlunabhaengig);
  (b) Abschirmung: Atomtraeger je Orbital (Norm der Koeffizienten > screen), Paar nur ueber gemeinsame Atome und nur,
  wenn die Schranke des Drehwinkels > screen ist; (c) sparse sphaerische Transformation; (d) Schwerpunkte und
  Dipol-Abschaetzung (distant) aus einem threadparallelen Produkt Dip L (die BLAS ist sequentiell), das die
  distant-Rechnung wiederverwendet. Boys bleibt auf dem Matrixweg (3 Matrizen; eine paarparallele Variante war bis
  Pentacen langsamer und wurde verworfen).
- **Ergebnisgleichheit** (gleiche SCF, beide Wege im selben Prozess, `CURCUMA_LMP2_LOC_COMPARE`): |dE_corr| <= 2.2e-14 Eh,
  |dF| <= 7.1e-13 auf C8, C20, C40, C60, C80, Acenen 1-5, (H2O)4-20; Boys (Benzol, C20, Pentacen) <= 1.4e-14 Eh. Mit
  screen 1e-6 statt 1e-8: dE bis 3.3e-11 Eh und F um bis 1.7e-10 kleiner (C80), daher Standard 1e-8.
- **Zeit Lokalisierung** (ms, 24 Threads, je ein Lauf; vorher = Binary vor der Aenderung, C60/C80 aus den Laeufen vom
  18.10.; C40, Pentacen und (H2O)20 vorher bei Last ~24 durch einen fremden Job, nachher Last < 1 bzw. C80 Last ~17):

| System | n | o | vorher | classic (neue Infrastruktur) | fast | Sweeps | Anteil an LMP2 (fast) |
|---|---|---|---|---|---|---|---|
| C20H42 | 490 | 61 | 73 | 41 | 25 | 6 | 0.9 % |
| C40H82 | 970 | 121 | 468 | 210 | 96 | 6 | 0.9 % |
| C60H122 | 1450 | 181 | 1618 | 657 | 213 | 6 | 0.9 % |
| C80H162 | 1930 | 241 | 7763 | 5593 | 386 | 6 | 1.2 % |
| Pentacen | 378 | 51 | 87 | 67 | 31 | 19 | 0.6 % |
| (H2O)20 | 480 | 80 | 81 | 51 | 28 | 5 | 0.7 % |

- **Exponenten** (gegen n): vorher C40 -> C80 4.1, C60 -> C80 5.5; fast C40 -> C80 2.0 (Sweeps 1.6, Dipolintegrale 1.5,
  Cholesky-Start 3.0). Bei C80: Dipolintegrale 117, Cholesky-Start 141, Sweeps 102, Schwerpunkte 12, U 14 ms.
- **Abschirmung greift erst ab C60:** Traeger bei 1e-8 umfassen bei C40 113 von 122 Atomen, bei C80 146 von 242; 66 % der
  Paare werden bei C80 noch ausgewertet. Die Sweeps sind daher bis C80 praktisch O(o^2 n), nicht linear.
- **Rest-N^3:** Cholesky-Start (D = C C^T und L^T S L, je ein sequentieller BLAS-Aufruf, n^2 o) ist bei C80 der groesste
  Teil (37 %); nicht geaendert, weil eine andere Rundung die Pivotwahl bei gleichwertigen Atomen aendern kann.
- **Naechstgroesste Teile von LMP2 bei C80** (GPU-SCF, ms): integraldirektes (P|i r) 18963, Paare/PNOs 7872,
  Iterationen 1923, Ueberlappungen 873, distant 58 (vorher 245); nicht angefasst.
- **Nicht getestet:** Systeme > C80, kompakte 3D-Systeme > 60 Atome, andere Elemente und Basen, Wirkung der Abschirmung
  dort, wo sie viele Paare ueberspringt (erst jenseits C80 zu erwarten); Thread-Unabhaengigkeit von `fast` nur durch
  Konstruktion und ctest (1 Thread Essigsaeuredimer/Butan), nicht systematisch.

**Abgrenzung:** kein Nachbau von DLPNO-CCSD(T) (zu gross, viele Naeherungsschichten), sondern ein lokales MP2 nach
dem Muster Pulay-Saebo-LMP2 + PNO. Die ueberlappenden Domaenen sind genau die Schwarz-artigen Teilgebiete aus WP14,
aber auf der Ebene des **Wellenfunktionsraums**, nicht der Metrikinversion -- dort liegt der Skalierungshebel.
Ehrlicher Rahmen: das ist **Neuimplementierung bekannter Methoden** (LMP2, PNO-LMP2, DLPNO-MP2), keine neue Methode;
die Analogie "Domaene = Schwarz-Teilgebiet" ist eine Lesart, kein Befund. Aus der Literatur erinnert, nicht
nachgeschlagen.

Stufen, jede mit eigenem Pruefkriterium gegen kanonisches MP2 (WP11):

1. **Lokalisierung** der besetzten Orbitale (Pipek-Mezey oder Boys; Energie des SCF invariant, das ist der Test).
2. **Paarvorauswahl:** weit entfernte Paare (i,j) nur ueber ein Abstands-/Dipolabschaetzung der Paarenergie; Schwelle
   einstellbar, Fehler gemessen.
3. **Domaenen:** je Paar die Vereinigung atomzentrierter Bereiche (PAOs oder projizierte AOs) um die beteiligten
   lokalen Orbitale; **Ueberlappung zwischen Paaren** ist beabsichtigt.
4. **PNO-Kompression:** je Paar die MP2-Paardichte diagonalisieren (semikanonische Amplituden), Natuerliche Orbitale mit
   Besetzung unter Schwelle (`T_CutPNO`) verwerfen. Schwelle einstellbar; Korrelationsenergie-Rueckgewinn gegen
   kanonisch ist die Messgroesse.
5. **Iteration** der Residualgleichungen im lokalen Raum (Fock-Matrix dort nicht diagonal), Konvergenz und Iterationszahl.

**Neuheit (Okt 2026, flache Suche):** zwei Standard-Websuchen fanden keine Arbeit, die einen additiven Schwarz-Vorkonditionierer
auf die gekoppelten LMP2/PNO-Amplitudengleichungen anwendet; das beweist keine Neuheit (zwei Suchen, abweichende Fachsprache
moeglich: "Domaene", "Fragment", "Divide-and-Conquer"). Bekannt ist: die PNO-Gleichungen werden iterativ geloest, mit einem
diagonalen Vorkonditionierer aus Pseudo-Virtualenergien. Die Kopplung zwischen Paaren laeuft ueber die Nichtdiagonalelemente
der Fock-Matrix im lokalen Orbitalraum; nur dort koennte ein Schwarz-artiger Vorkonditionierer ansetzen. Erwartung (Hypothese):
die diagonale Variante konvergiert schon in wenigen Iterationen, der Gewinn waere also klein; wie bei EEQ zaehlt der Vergleich
gegen den bestehenden Vorkonditionierer. Fuer eine Neuheitsaussage braucht es eine gruendliche Literaturrecherche.

**Kriterien (vorab):** Rueckgewinn der kanonischen MP2-Korrelationsenergie ueber mindestens drei Systemklassen
(Alkane, Aromaten, Wasserstoffbruecken-Cluster) und je mehrere Groessen; Fehler pro Schwelle als Kurve, nicht ein Punkt.
Die Schwellen sind Parameter, keine Konstanten.

**Wichtige Randbedingungen:**
- Auf kleinen Systemen (Benzol, 114 Funktionen) bringt lokale Korrelation **nichts**; der Break-even liegt erfahrungsgemaess
  bei einigen zehn schweren Atomen (nicht gemessen). Testsaetze deshalb Alkanketten und Wassercluster, solange das
  kanonische MP2 noch als Referenz rechenbar ist (Speichergrenze der Maschine pruefen).
- Referenzen: kanonisch aus WP11 und PySCF; **ORCA `DLPNO-MP2`** ist als externer Vergleich verfuegbar. Der Vergleich mit
  ORCA gilt nur qualitativ (andere Domaenen und Schwellen), PySCF/kanonisch ist die exakte Referenz.
- Abhaengigkeit: WP11 und WP10 (RI-Integrale fuer die Paar-Integrale `(ia|jb)` in der Domaene).
- Nicht im Umfang: Gradient, offene Schalen, Triples, CCSD.

## WP14 -- Bewertung: "Schwarz-Attention" fuer Matrixoperationen in der QM-Engine

**Begriffsklaerung vorab.** "Schwarz" meint hier die **Schwarzsche Alternierungsmethode** (Gebietszerlegung mit
Ueberlappung, der Vorkonditionierer-Gedanke aus dem Paper von Koehler und Rheinbach, in den Vault-Notizen
*Overlapping Schwarz Attention* und *curcuma Schwarz-Attention*). Das ist **nicht** die Cauchy-Schwarz-Ungleichung
der Integralabschaetzung `|(ab|cd)| <= sqrt((ab|ab)(cd|cd))`, die in `qm_integrals.cpp` schon laeuft
(`-qm.eri_screening`). Der Namensanklang ist Zufall; beide Dinge nicht vermischen.

**Was bisher im Vault belegt ist** (Notiz *curcuma Schwarz-Attention*, Stand 10.10.): Fuer die EEQ-Matrix von
GFN-FF -- langreichweitig, Inverse lokal -- schlaegt die **ungelernte** lokale Inverse (`localinv`, ueberlappende
Atomnachbarschaften, exakt invertierte Bloecke) den Jacobi-Vorkonditionierer deutlich (PCG-Iterationen 68 -> 15 auf
7320 Atomen), und eine gelernte Variante hat gegen sie keinen erkennbaren Vorteil; das Paper selbst zeigt nur 1D-Poisson
mit n=1. Daraus folgt die Arbeitshypothese fuer die QM-Engine: **zuerst der ungelernte Schwarz-artige Loeser, gelernt
nur bei nachgewiesener Luecke.**

Die folgende Tabelle ist eine **Hypothesenliste aus dem Quelltext und der Literatur, nichts davon gemessen.**

| Matrix-Operation in der Engine | Struktur | Passt der Schwarz-Gedanke? | Bewertung |
|---|---|---|---|
| **RI-Metrik** `J_PQ`, `J^-1` bzw. `J^-1/2` (WP10, WP11) | Coulomb-Matrix, langreichweitig; ob die Inverse im Hilfsbasissatz lokal ist, ist **nicht gemessen** (bei EEQ ja, dort aber eine Punktladungs-Matrix mit einer Funktion je Atom; ein Hilfsbasissatz hat viele, stark ueberlappende Funktionen je Atom) | ueberlappende atomzentrierte Fit-Domaenen mit lokaler exakter Inverse (lokales Dichtefitting, aus der Literatur erinnert, nicht nachgeschlagen) | **Nebenkandidat, kleiner Hebel** (Kostenrechnung unten). Erst Python-Messung der Inverse-Abfallstruktur |
| **Fock-Diagonalisierung** (O(n^3), `solveFock`) | dicht, nicht lokal im Spektrum | Gebietszerlegung entspricht Divide-and-Conquer / lokalisierten Orbitalen; fuer xTB existiert `-large_system_mode dc` | Fuer die ab-initio-Engine erst sinnvoll bei vielen hundert Atomen; Genauigkeitsverlust gegenueber exaktem Eigenproblem ist **die** Kosten-Nutzen-Frage. Niedrige Prioritaet |
| **SCF-Konvergenz** (DIIS/Level-Shift) | iterative Fixpunktsuche | Vorkonditionierung des SCF-Schrittes (dielektrische Antwort, siehe Vault *SCF-Vorkonditionierung und dielektrische Antwort*, bei GFN2 schon umgesetzt) | Uebertragbarkeit auf HF plausibel, aber HF hat anders strukturierte Antwort; nur nach Messung der Iterationszahlen auf grossen Systemen |
| **J/K-Aufbau** | Kontraktion mit 4-Index-Tensor | Cauchy-Schwarz-Screening (existiert); Schwarz-Zerlegung hilft hier nicht, die Arbeit ist die Integralauswertung | Kein Ansatzpunkt fuer Gebietszerlegung; Screening ist bereits der richtige Mechanismus |
| **Amplitudengleichungen** (WP13) | nicht-lineare Gleichungen, Vorkonditionierung durch Orbitalenergie-Nenner | Lokale Korrelation (PNO/DLPNO) *ist* eine Gebietszerlegung mit Ueberlappung, aber mit physikalisch begruendeten Domaenen | Gelerntes `M^-1` fuer die Amplitudenstartwerte ist reine Spekulation; DLPNO ist der etablierte Weg und ausdruecklich ausserhalb des Umfangs |
| **MO-Transformation** | GEMM-Folge | keine lokale Struktur noetig | Kein Ansatzpunkt |

**Folgerungen:**

1. **Korrektur (Okt 2026, auf Nachfrage des Betreibers):** die erste Fassung nannte RI-J den aussichtsreichsten
   Kandidaten. Das war zu stark. Kostenrechnung (n Basisfunktionen, N_aux Hilfsfunktionen, grobe Ordnungen):

   | Groesse | RI-J-Kontraktion n^2 N_aux (je SCF-Iteration) | Metrik-Anwendung N_aux^2 (je Iteration) | Metrik-Zerlegung N_aux^3 (einmalig) |
   |---|---|---|---|
   | n=114, N_aux~500 (Benzol/def2-SVP) | 6.5e6 | 2.5e5 | 1.2e8 |
   | n=1000, N_aux~5000 | 5.0e9 | 2.5e7 | 1.2e11 |
   | n=5000, N_aux~25000 | 6.2e11 | 6.2e8 | 1.6e13 (dichte Inverse 5 GB) |

   (Annahme N_aux ~ 5n; Rechnung hier im Gespraech mit Python-Taschenrechner, nicht an der Engine gemessen.)
   Fuer HF mit RI-J ist die Metrik also **nicht** der Engpass: sie macht ein Hundertstel der Kontraktion aus, die
   einmalige Zerlegung wird erst ab einigen tausend Basisfunktionen und wegen des Speichers relevant. Der Hebel bei
   RI-J ist die Dreizentren-Integralformung samt Screening, nicht die Metrik.
   Im **RI-MP2** ist die Metrik-Anwendung `J^-1/2` auf `(ia|Q)` mit o v N_aux^2 spuerbarer, aber immer noch
   kleiner als die Energiekontraktion o^2 v^2 N_aux (Benzol/def2-SVP: 4.9e8 gegen 1.9e9, bei o=100, v=900,
   N_aux=5000: 2.2e12 gegen 4.0e13). Auch dort ist der N^5-Hauptterm eine GEMM-Kontraktion, kein Loeser-Problem.
   Der Betreiber hatte die Vermutung, dass Schwarz-artige Ansaetze eher bei korrelierten Methoden greifen; die
   Rechnung stuetzt das **teilweise**: die Metrik spielt dort eine groessere Rolle als bei RI-J, der wirkliche
   Skalierungshebel der korrelierten Methoden ist aber die lokale Korrelation (Paar- und Virtuellraum-Truncation,
   PNO/DLPNO), also eine Gebietszerlegung des Wellenfunktionsraums und nicht der Metrikinversion. Das ist
   ausdruecklich ausserhalb des Umfangs (siehe WP13).
2. Ein **gelerntes** Modul ist fuer QM-Matrizen nicht begruendet: Trainingsdaten waeren exakt und kostenlos (wie bei
   EEQ), aber die EEQ-Messung zeigt, dass die exakte lokale Inverse billig ist und kaum Raum laesst. Das gilt hier
   bis zum Gegenbeweis ebenso.
3. Ein ML-Anteil fuer die **Energie** (Lesart B der Vault-Notiz, Delta-Lernen) ist ein anderes Vorhaben und gehoert
   nicht in diesen Plan.

**Messung Okt 2026** ([../LOCAL_MP2_FEASIBILITY_2026-10.md](../LOCAL_MP2_FEASIBILITY_2026-10.md), Q4): die Inverse
der JKFIT-Metrik ist nicht lokal (bei 30-50 Bohr nur 1-1.5 Groessenordnungen unter J); robuste lokale Fits erfuellen den
ersten Falsifikator fuer RI-J (<= 6e-8 Eh/Atom bei 8 Bohr) und RI-MP2 in Alkanen/Wasser, verfehlen ihn fuer RI-MP2 in
Acenen ab Anthracen (1.3-2.2e-6); lokales Cholesky schneller erst ab naux ~7000-14000, dort dauert das globale 0.3-2 s.

**Falsifikatoren (vor der Messung festgelegt):**

- RI mit lokaler Fit-Domaene: Gesamtenergiefehler gegen globales RI kleiner als 1e-6 Eh pro Atom bei Domaenenradius
  nicht groesser als 8 Bohr, sonst keine Verwendung. Ueber mindestens drei Systeme und drei Groessen, nicht einzeln.
- Skalierung: lokale Variante muss ab einer messbaren Systemgroesse schneller sein als das globale Cholesky der
  Metrik, sonst bleibt es beim globalen Verfahren.
- Kein Befund aus einer einzelnen Konfiguration (Paper und Vorarbeit liefern n=1).

## Offene Entscheidungen fuer den Betreiber

0. Betreibervorschlag Okt 2026: den WP14-Entwurf erst nach der Python-Messung der Metrik-Inversen von einem starken
   Modell (Opus oder Fable) durchdenken lassen.

1. ~~MP5 woertlich?~~ Erledigt: MP2 (O(N^5)), WP12 gestrichen.
2. Hilfsbasissaetze: def2/J und def2/JK liegen vor (BSE, H-Ne, Okt 2026); offen ist die Lizenzfrage. RI-MP2 nutzt
   nach Betreiberentscheidung def2-universal-JKFIT/-JFIT als Korrelations-Hilfsbasis; ob def2-SVP/C nachgeladen
   werden soll, ist offen.
3. (T)-Korrektur und Gradienten in den Plan aufnehmen?
4. ~~Vault-Fuehrung der Schwarz-Bewertung~~ Erledigt: Pfad 5 in der Vault-Notiz *curcuma Schwarz-Attention*.

## BLAS (Korrektur, Okt 2026)

- **Befund:** mehrere Abschnitte oben (RI-MP2, RI-Gradient, MP2-Gradient) nennen die einthreadige Referenz-netlib-BLAS
  als Bremse des Release-Builds. Das ist fuer das Binary falsch. `LD_DEBUG=bindings ./curcuma ...` zeigt, dass `dgemm_`,
  `dgemv_`, `dtrsm_` aus `libmkl_intel_lp64.so.3` kommen; `libblas.so.3` (netlib) wird nur als Abhaengigkeit von
  `liblapack.so.3` geladen. Das Binary linkt `libmkl_intel_lp64`, `libmkl_sequential`, `libmkl_core` (Threading-Schicht
  **sequential**, also einfaedig), obwohl `USE_MKL=OFF` im Cache steht (CMake findet MKL ueber `find_package(BLAS)`).
- **Gemessen** (Eigenstaendiges Testprogramm, Core Ultra 9 285K, n=1 je Zahl): `dgemm` 2000^3: netlib 4.4 GFLOPS, MKL
  sequential 76.7, MKL gnu_thread 8 Threads 566, 16 Threads 751. Eigen-Produkt 199 x 1230 x 9353 (`EIGEN_USE_BLAS`):
  netlib allein 4.2 GFLOPS, MKL sequential 68.3, Eigen eigene GEMM 23.9 (1 Thread) und 188 (OpenMP, 16 Threads).
  Die 2.7 bis 4 GFLOPS der Agenten stammen aus einem separat gegen `-lblas` gelinkten Test, nicht aus dem Binary.
- **Folge:** `qmint::gemmABt` (29 GFLOPS je Thread) ist pro Thread langsamer als MKL sequential (68), gewinnt nur durch die
  Parallelisierung ueber Threads. Ob MKL sequential mit OpenMP-Schleifen aussenherum oder MKL mit Threads schneller ist,
  ist **nicht in den Produktivpfaden gemessen**. Offen: MKL-Threading-Schicht fuer grosse dichte Operationen (Eigenloeser,
  Eigen-GEMM im Host-SCF) bei n > 1000.
- **Test mit gethreadeter MKL (`-DUSE_MKL=ON`, Build `build_mkl/`, nur CPU, Okt 2026, n=1 bis 3 je Zahl, 16 Threads):** das
  Projekt unterstuetzt `USE_MKL=ON` (explizit `mkl_gnu_thread`; `find_package(BLAS)` loest sonst auf die sequentielle MKL
  auf, CMake-Kommentar WP1 Juni 2026); `release/` ist ohne diese Option konfiguriert. Ergebnis: HF, MP2, RI-MP2 keine
  Aenderung (Benzol HF 1.01 gegen 1.03 s, Naphthalin HF 7.50 gegen 7.65 s, MP2-Transformation Naphthalin 552 gegen 541 ms),
  RI-Setup Naphthalin 148 gegen 238 ms (eine Messung, nicht wiederholt). GFN2 (`-method gfn2`, dichter Eigenloeser)
  Polystyrol20 (325 Atome): 1.30 s gegen 0.68 s (3 Laeufe je, Streuung unter 2%), GFN1 1.29 gegen 0.60 s; bei 173 Atomen
  (jl22-bgal) 0.25 gegen 0.26 s (5% langsamer). 48 Energien (gfn2 und gfn1, 24 Strukturen bis 120 Atome) beider Binaries
  stimmen auf alle 8 gedruckten Dezimalstellen ueberein. Nicht gemessen: Systeme ueber 325 Atome, MD und Mehrmolekuel-Pools
  (Oversubscription), `-gpu`-Pfade (der Testbuild hat kein CUDA-Plugin), ConfScan, vollstaendiges `ctest`.
