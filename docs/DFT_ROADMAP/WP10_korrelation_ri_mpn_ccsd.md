# WP10-WP14 -- Dichtefitting (RI), MPn und CCSD, plus Bewertung von Schwarz-Ansaetzen

- **Status:** 🤖 AI-generated Plan. Umgesetzt (Okt 2026, AI-generated, machine-tested): nur WP11 Schritt 1 (kanonisches
  RMP2, CPU + CUDA) und Schritt 3 (SCS/SOS), siehe WP11. Alles andere ist nicht implementiert.
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

## WP10 -- Dichtefitting (RI-J, RI-JK, RI-C)

**Ziel:** 3-Index-Integrale `(mn|P)` mit Hilfsbasis `P`, Metrik `J_PQ = (P|Q)`, Fitkoeffizienten
`c_P = sum_Q J^-1_PQ (Q|ls) P_ls`. Speicher O(n^2 N_aux) statt O(n^4).

- **Nutzen:** RI-J (Coulomb) und RI-JK fuer HF/DFT; **RI-C** (Korrelations-Hilfsbasis) macht MP2 und
  CCSD ueberhaupt erst fuer mittelgrosse Systeme moeglich. Auf der GPU sind es GEMMs.
- **Preis, ehrlich:** die Energie aendert sich um etwa 1e-5 Eh pro Atom (hilfsbasisabhaengig). Das ist eine
  **Methodenvariante** (`-qm.ri on|off`, Standard aus), kein Ersatz des exakten Pfads; der exakte Pfad bleibt die
  Referenz fuer alle bestehenden Tests.
- **Fehlt im Repo:** Hilfsbasissaetze. Vorhanden sind nur `def2-SVP`, `def2-SV(P)`, `MINIX` (Orbitalbasen).
  Benoetigt: `def2/J`, `def2/JK` und `def2-SVP/C` (Weigend et al.), im selben Dateiformat. Beschaffung und
  Lizenz/Quelle (EMSL/Basis-Set-Exchange) klaeren, nichts aus dem Gedaechtnis eintippen.
- **Metrikinverse** (Kostenanteil klein, siehe WP14): `J_PQ` ist SPD, aber schlecht konditioniert. Optionen: Cholesky mit Pivotierung,
  Eigenzerlegung mit Abschneiden kleiner Eigenwerte (`J^-1/2`), oder Cholesky-Zerlegung der AO-ERI
  (pivotiert, keine Hilfsbasis, Fehler ueber Schwelle einstellbar -- naeher an der Referenz, siehe WP14).
- **Validierung:** gegen PySCF `density_fit()` mit identischer Hilfsbasis (Erwartung: Uebereinstimmung nahe
  Rundungsfehler, weil dasselbe Verfahren) und gegen ORCA `RIJK`/`RI-MP2`; zusaetzlich der Fehler gegen den exakten Pfad
  als Messgroesse, nicht als Testkriterium.
- **Nicht im Umfang:** COSX/RIJCOSX (ORCA-spezifisch), Gradient (spaeter, eigener Schritt).

## WP11 -- MO-Transformation und MP2

**Stand Okt 2026: Schritt 1 (CPU und GPU) und Schritt 3 umgesetzt; 🤖 AI-generated, ⚙️ machine-tested, menschlicher
Produktionstest steht aus.** Schritt 2 (RI-MP2) offen (wartet auf WP10).

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
- **Gradient:** keiner. Eine Gradientenanforderung ist ein Fehler ("mp2 has no analytic gradient"): `-opt` und
  `-sp -gradient` brechen damit ab (geprueft, `qm_mp2_opt_error`); `-hessian` meldet ihn fuer jede Auslenkung, laeuft
  aber weiter und schreibt keine Frequenzen (Verhalten von `hessian.cpp`, nicht geaendert); `-md` nicht getestet.
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
  | Benzol voll (114, 21, 93) | 3517 / 6.0 ms | 913 / 1.0 ms | 64 / 0.3 ms |
  | Benzol fc (114, 15, 93) | 2458 / 3.6 ms | 656 / 0.5 ms | 63 / 0.2 ms |
  | Naphthalin voll (180, 34, 146) | 57171 / 41 ms | 10196 / 6.8 ms | 199 / 0.7 ms |
  | Naphthalin fc (180, 24, 146) | nicht gemessen | 7364 / 3.0 ms | 135 / 0.4 ms |

  CPU-Naphthalin mit `-qm.scf_direct off` (Tensor gespeichert); mit Standard (direkte SCF) kommt der Neuaufbau des
  Tensors dazu (8398 MB, 1572 ms). GPU mit `-qm.gpu_direct on`: Neuaufbau auf dem Host + Upload (599 ms) vor der
  Transformation. Device-Speicher Naphthalin (gerechnet, nicht gemessen): Tensor 8398 MB + Batch 34*84.5 MB.
- **Nicht getestet:** offene Schale (nicht implementiert), Elemente ausser H, C, N, O, F (Li-Ne-frozen-core-Regel nur auf
  C, N, O, F gegen PySCF/ORCA geprueft), andere Basen als def2-SVP, kartesisch 6d, Systeme ueber 180 BF, ROCm/Vulkan
  (keine Implementierung, `-gpu rocm` faellt gezaehlt auf die CPU), geladene Molekuele.
- **Nicht implementiert gegenueber ORCA/PySCF:** RI-MP2, integral-direktes MP2, MP2-Gradient und -Dichte,
  frei waehlbare SCS-Faktoren, Ausnutzung der Permutationssymmetrie in der Transformation (Faktor ~2 offen).

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
2. Hilfsbasissaetze beschaffen (Quelle/Lizenz) -- Voraussetzung fuer WP10/WP11-RI.
3. (T)-Korrektur und Gradienten in den Plan aufnehmen?
4. ~~Vault-Fuehrung der Schwarz-Bewertung~~ Erledigt: Pfad 5 in der Vault-Notiz *curcuma Schwarz-Attention*.
