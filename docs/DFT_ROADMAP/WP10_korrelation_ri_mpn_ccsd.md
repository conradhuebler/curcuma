# WP10-WP14 -- Dichtefitting (RI), MPn und CCSD, plus Bewertung von Schwarz-Ansaetzen

- **Status:** WP10 RI-J/RI-JK (HF, CPU + CUDA) implementiert, 🤖 AI-generated, ⚙️ machine-tested (Okt 2026, Abschnitt
  WP10 unten); RI-C, WP11-WP15 sind Plan.
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
  **def2-SVP/C (RI-C) fehlt**: RI-MP2 ist nicht implementiert und ist nach dem Merge mit WP11 zu planen.
- **Nicht implementiert:** RI-C/RI-MP2, Gradient mit RI (`-opt`/`-gradient` mit `-qm.ri j|jk` brechen mit Fehler ab,
  `ctest qm_ri_gradient_error`), hf-3c mit RI (ignoriert mit Warnung: HF-3c ist mit exakten Integralen in MINIX
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

## WP11 -- MO-Transformation und MP2

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
2. Hilfsbasissaetze: def2/J und def2/JK liegen vor (BSE, H-Ne, Okt 2026); offen sind die Lizenzfrage und
   def2-SVP/C fuer RI-MP2.
3. (T)-Korrektur und Gradienten in den Plan aufnehmen?
4. ~~Vault-Fuehrung der Schwarz-Bewertung~~ Erledigt: Pfad 5 in der Vault-Notiz *curcuma Schwarz-Attention*.
