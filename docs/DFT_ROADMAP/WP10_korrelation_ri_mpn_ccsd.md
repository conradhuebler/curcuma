# WP10-WP14 -- Dichtefitting (RI), MPn und CCSD, plus Bewertung von Schwarz-Ansaetzen

- **Status:** 🤖 AI-generated Plan, **nichts davon ist implementiert oder getestet.** Stand Okt 2026.
- **Abhaengigkeit:** WP3 (RHF, ERI-Engine, `QMEngine`), fuer die GPU-Teile die Stufen Q-G1/Q-G2 aus
  [../QM_GPU_ROADMAP.md](../QM_GPU_ROADMAP.md).
- **Referenzen auf dieser Maschine:** PySCF (venv `~/.venvs/qmref`), ORCA 6.1 (`~/Downloads/orca_6_1_0_...`).
  Jeder Test muss spaeter auf einer anderen Maschine noch einmal gegen ORCA laufen (Betreibervorgabe).
- Reihenfolge und Abhaengigkeiten: WP10 (RI) -> WP11 (MO-Transformation + MP2) -> WP12 (MP3-MP5) und
  WP13 (RCCSD); WP14 ist die Auswertung der Schwarz-Frage und haengt an WP10.

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
- **Metrikinverse:** `J_PQ` ist SPD, aber schlecht konditioniert. Optionen: Cholesky mit Pivotierung,
  Eigenzerlegung mit Abschneiden kleiner Eigenwerte (`J^-1/2`), oder Cholesky-Zerlegung der AO-ERI
  (pivotiert, keine Hilfsbasis, Fehler ueber Schwelle einstellbar -- naeher an der Referenz, siehe WP14).
- **Validierung:** gegen PySCF `density_fit()` mit identischer Hilfsbasis (Erwartung: Uebereinstimmung nahe
  Rundungsfehler, weil dasselbe Verfahren) und gegen ORCA `RIJK`/`RI-MP2`; zusaetzlich der Fehler gegen den exakten Pfad
  als Messgroesse, nicht als Testkriterium.
- **Nicht im Umfang:** COSX/RIJCOSX (ORCA-spezifisch), Gradient (spaeter, eigener Schritt).

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

## WP12 -- MP3, MP4, MP5

**Praezisierung noetig:** gemeint ist hier die Moeller-Plesset-Reihe in dritter bis fuenfter Ordnung (MP5
woertlich genommen). Falls eigentlich MP2 gemeint war, ist das mit WP11 erledigt.

- **MP3:** zusaetzlich zu MP2 zweite Ordnung der Amplituden, Skalierung O(n^6) -- im Wesentlichen ein
  CCD-Schritt mit festem T2. Moderat, gut machbar.
- **MP4 (SDQ, dann mit T):** SDQ O(n^6); die Tripel (T) O(n^7) und deutlich mehr Terme.
- **MP5:** O(n^8), sehr viele Terme (Singles bis Quadruples); handgeschriebene Formeln sind fehleranfaellig.
  Realistischer Weg: ueber ein generisches Tensorkontraktionsgeruest (Diagramme/Einstein-Summen) oder Ableitung
  aus der CC-Maschinerie (WP13), dann Vergleich gegen eine unabhaengige Referenz.
- **Ehrliche Einordnung:** die MPn-Reihe divergiert oder oszilliert fuer viele Systeme; MP3-MP5 sind in der
  Praxis kaum noch im Einsatz (CC ist robuster). Der Wert ist hier **didaktisch** (Konvergenzverhalten der Reihe
  zeigen) und als Test der Tensor-Infrastruktur.
- **Referenz:** nicht geprueft, welche der Programme auf dieser Maschine MP3-MP5 anbieten (PySCF: MP2 sicher,
  hoehere Ordnungen vermutlich nicht; ORCA: vor dem Einsatz in der Handbuchversion 6.1 nachlesen). Eine
  unabhaengige Pruefung von MP3-MP5 ist in kleinen Systemen (H2, HeH+, LiH in minimaler Basis) moeglich, indem die
  Stoerungsreihe im vollen Determinantenraum (volle CI-Matrix, PySCF `fci`) mit Hamiltonian-Vektor-Produkten
  ausgewertet wird -- aufwendig, aber unabhaengig von der Diagrammformel-Implementierung.

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
| **RI-Metrik** `J_PQ` (WP10) und `J^-1` | Coulomb-Matrix, langreichweitig, Inverse lokal -- **dieselbe Struktur wie die EEQ-Matrix** | Ueberlappende atomzentrierte Fit-Domaenen mit lokaler exakter Inverse = lokales Dichtefitting (aus der Literatur erinnert, nicht nachgeschlagen: atompaar-/domaenenlokales Dichtefitting, wie es lokale Korrelationsmethoden nutzen) | **Aussichtsreichster Kandidat.** Wirkt als Loeser/Vorkonditionierer oder als Naeherung der Fitkoeffizienten. Neu zu messen: Energiefehler gegen globales DF, Konditionszahl, Skalierung |
| **Fock-Diagonalisierung** (O(n^3), `solveFock`) | dicht, nicht lokal im Spektrum | Gebietszerlegung entspricht Divide-and-Conquer / lokalisierten Orbitalen; fuer xTB existiert `-large_system_mode dc` | Fuer die ab-initio-Engine erst sinnvoll bei vielen hundert Atomen; Genauigkeitsverlust gegenueber exaktem Eigenproblem ist **die** Kosten-Nutzen-Frage. Niedrige Prioritaet |
| **SCF-Konvergenz** (DIIS/Level-Shift) | iterative Fixpunktsuche | Vorkonditionierung des SCF-Schrittes (dielektrische Antwort, siehe Vault *SCF-Vorkonditionierung und dielektrische Antwort*, bei GFN2 schon umgesetzt) | Uebertragbarkeit auf HF plausibel, aber HF hat anders strukturierte Antwort; nur nach Messung der Iterationszahlen auf grossen Systemen |
| **J/K-Aufbau** | Kontraktion mit 4-Index-Tensor | Cauchy-Schwarz-Screening (existiert); Schwarz-Zerlegung hilft hier nicht, die Arbeit ist die Integralauswertung | Kein Ansatzpunkt fuer Gebietszerlegung; Screening ist bereits der richtige Mechanismus |
| **Amplitudengleichungen** (WP13) | nicht-lineare Gleichungen, Vorkonditionierung durch Orbitalenergie-Nenner | Lokale Korrelation (PNO/DLPNO) *ist* eine Gebietszerlegung mit Ueberlappung, aber mit physikalisch begruendeten Domaenen | Gelerntes `M^-1` fuer die Amplitudenstartwerte ist reine Spekulation; DLPNO ist der etablierte Weg und ausdruecklich ausserhalb des Umfangs |
| **MO-Transformation** | GEMM-Folge | keine lokale Struktur noetig | Kein Ansatzpunkt |

**Folgerungen:**

1. Die einzige Stelle mit derselben mathematischen Struktur wie der belegte Fall ist die **RI-Metrik**. WP14
   beginnt daher mit WP10 und einer Python-Rekonstruktion (wie in der EEQ-Erkundung: erst Python, dann C++).
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

1. **MP5 woertlich?** (siehe WP12; MP2/MP3 sind die ueblichen Ziele.)
2. Hilfsbasissaetze beschaffen (Quelle/Lizenz) -- Voraussetzung fuer WP10/WP11-RI.
3. (T)-Korrektur und Gradienten in den Plan aufnehmen?
4. Soll der Schwarz-Teil (WP14) im Vault als eigener Pfad in *curcuma Schwarz-Attention* gefuehrt werden
   (Vault-Regeln, Pfad 5: QM-Matrizen) oder nur hier?
