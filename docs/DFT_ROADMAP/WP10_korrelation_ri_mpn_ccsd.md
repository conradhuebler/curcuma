# WP10-WP14 -- Dichtefitting (RI), MPn und CCSD, plus Bewertung von Schwarz-Ansaetzen

- **Status:** WP10 RI-J/RI-JK (HF, CPU + CUDA), WP11 Schritt 1+3 (kanonisches RMP2, SCS/SOS, CPU + CUDA) und
  WP11 Schritt 2 (RI-MP2 mit def2-universal-JKFIT/-JFIT als Korrelations-Hilfsbasis, CPU + CUDA) sind
  implementiert, 🤖 AI-generated, ⚙️ machine-tested (Okt 2026, Abschnitte WP10, WP11 und RI-MP2 unten), ebenso der
  analytische Gradient fuer RI-J/RI-JK-HF (Abschnitt "WP10 -- Gradient") und der analytische Gradient des kanonischen
  MP2/SCS/SOS (CPU, 2e-Ableitungsteil auch CUDA; Abschnitt "WP11 -- MP2-Gradient"); RI-MP2 mit def2-SVP/C, der
  RI-MP2-Gradient und WP13-WP15 sind Plan.
- **Abhaengigkeit:** WP3 (RHF, ERI-Engine, `QMEngine`), fuer die GPU-Teile die Stufen Q-G1/Q-G2 aus
  [../QM_GPU_ROADMAP.md](../QM_GPU_ROADMAP.md).
- **Referenzen auf dieser Maschine:** PySCF (venv `~/.venvs/qmref`), ORCA 6.1 (`~/Downloads/orca_6_1_0_...`).
  Jeder Test muss spaeter auf einer anderen Maschine noch einmal gegen ORCA laufen (Betreibervorgabe).
- Reihenfolge und Abhaengigkeiten: WP11 (kanonisches MP2, O(N^5)) zuerst, dann WP10 (RI) und RI-MP2, danach
  WP13 (RCCSD); WP12 (MP3-MP5) ist **gestrichen** (Betreiber, Okt 2026: gemeint war MP2). WP14 ist die Auswertung
  der Schwarz-Frage.

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
  (16000), `-qm.gpu_ri auto|on|off`.
- **Verfahren:** RI-J: `c = J^-1 g`, `g_P = (P|ls) D_ls`, `J_mn = (mn|P) c_P`; Austausch exakt (gespeicherter Tensor
  oder integral-direkt wie ohne RI, entspricht ORCA `RIJONX`). RI-JK (Weigend 2002): `B = L^-1 (P|mn)` mit
  `J = L L^T`, `K_mn = sum_Q (B^Q D B^Q)_mn` ueber `D = V S V^T` (Eigenvektoren von D, mit sqrt|lambda| skaliert),
  einmal gespeichert (naux x n^2). Integrale: McMurchie-Davidson ueber Schalenpaare; eine Hilfsfunktion ist das Paar
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
  definiert, `qm_ri_hf3c_ignored`), offene Schalen, Elemente ab Na, integral-direktes RI, lokales Fitten,
  Speicherkompression des 3-Index-Tensors (voll n^2, nicht n(n+1)/2), device-residenter SCF-Loop mit RI (RI laeuft
  im Host-Loop).

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
- Das System-BLAS dieser Maschine ist die Referenz-BLAS (netlib, ein Thread); die schweren RI-Schritte
  (Vorwaertssubstitution, g = T vec(D), J, K) sind deshalb eigene OpenMP-Schleifen. Mit Eigen-GEMM auf netlib-BLAS
  dauerte das Setup fuer Naphthalin 8.5 s statt 0.1 s.

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
- **Zeiten** (RTX 5090, Hostspeicher 256 GB; CPU-BLAS dieser Maschine ist die **Referenz-netlib-libblas**,
  einfaedig, das bremst die CPU-Zahlen; MKL/OpenBLAS nicht gemessen). MP2-Anteil ohne SCF, Transformation / Energie:

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
- **Zeiten** (ms, RTX 5090, CPU-BLAS = Referenz-netlib, `curcuma -sp X -method mp2 -gradient -qm.scf_threshold 1e-10
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
- **Speicher:** Setup-Tensor naux n^2 (Grenze `-qm.ri_max_memory_mb`, Fehler mit Groessenangabe), B o v naux
  zweimal (vor/nach Transposition, gleiche Grenze). Anthracen: Tensor 596 MB, B 92 MB.
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

**Zeiten** (Intel Core Ultra 9 285K, RTX 5090, netlib-BLAS, `-qm.scf_threshold 1e-8`, volle Korrelation; ms des
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
  offene Schalen, integral-direkte 3-Index-Integrale (Setup-Tensor naux n^2 voll im Speicher), Batching der
  B-Transformation fuer sehr grosse Systeme, i <-> j-Symmetrie auf der GPU.

## WP12 -- MP3 bis MP5: gestrichen

Betreiberentscheidung Okt 2026: gemeint war MP2, und zwar die O(N^5)-Variante (kanonisches MP2 mit
4-Index-Transformation, danach RI-MP2, das formal ebenfalls N^5 skaliert, aber mit kleinerem Vorfaktor). MP3-MP5
sind nicht Teil des Plans; der fruehere Abschnitt (MP5 mit O(n^8), Pruefung ueber den vollen Determinantenraum)
steht in der Git-Historie dieser Datei.

## WP13 -- Closed-shell CCSD (RCCSD)

- Spinadaptierte Amplitudengleichungen (T1, T2), kanonische RHF-Referenz, DIIS auf den Amplituden,
  MO-Integrale aus WP11. Speicher O(v^4) fuer das `(ab|cd)`-Integral im Standardschema; Skalierung O(o^2 v^4).
- **Validierung:** gegen PySCF `cc.CCSD` und ORCA `CCSD` (`! CCSD def2-SVP`, eingefrorener Rumpf wie in
  der Referenz), Energien auf 1e-8 Eh; Konvergenzverhalten (Iterationszahl) als zweite Kontrolle.
- **(T)-Korrektur:** nicht angefragt; wird als sinnvolle Folgearbeit gefuehrt (CCSD(T) ist der Standard und die
  vorhandene Benchmark-Referenz DLPNO-CCSD(T) im Vault zielt darauf).
- **Nicht im Umfang:** DLPNO/Domaenen-Truncation, offene Schalen (UHF/ROHF-CC), Gradient und EOM.
- **Skalierung:** praktisch bis etwa 200 Basisfunktionen (Benzol/def2-SVP: 114 Funktionen, v ~ 93) im
  Speicher; die (ab|cd)-Schleife ist der GPU-Kandidat.

## WP15 -- Lokales MP2 mit PNO-Truncation in ueberlappenden Domaenen (Idee des Betreibers, Okt 2026)

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
