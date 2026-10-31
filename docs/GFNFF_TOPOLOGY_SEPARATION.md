# GFN-FF Topologie-Separation — Analyse

**Status**: 🤖 AI-generiert, ⚙️ maschinell gemessen (Juli 2026)
**Frage**: Lässt sich der topologische Teil von GFN-FF separieren, als Basis für
(a) einen molekülübergreifenden Delta-ML-Korrektursatz und (b) ein schnelles GFN-FF?
**Antwort**: Ja — die Separation existiert bereits strukturell, und der eingefrorene
Block ist **beweisbar eine reine Funktion des Molekülgraphen**, nicht der Ankergeometrie.

Diese Analyse korrigiert eine frühere pessimistische Einschätzung, die den
Topologie-Anker als blockierendes Hindernis für ein generisches Modell einstufte.
Die Messungen widerlegen das.

---

## 1. Kernbefund: die Topologie ist graph-abgeleitet, nicht anker-abgeleitet

### Code-Beleg

`EEQSolver::calculateTopologyChargesMultiRHS` (`eeq_solver.cpp:2884-2896`) baut die
Phase-1-Coulomb-Matrix aus **topologischen** Distanzen:

```cpp
Matrix topo_dist;
if (topology.has_value()) {
    // Dijkstra topological distances (sparse graph, O(N*E*logN))
    topo_dist = computeTopologicalDistancesSparse(atoms, *topology);
    ...
    double r = topo_dist(i, j);   // Topological distance in Bohr
```

Der geometrische Pfad (`eeq_solver.cpp:2917ff`) ist ausdrücklich nur **Fallback**
für Nicht-GFN-FF-Nutzer (D4). Referenz: `gfnff_ini.f90:431-461` (Floyd-Warshall
über die Bindungsmatrix, skaliert mit `gen%rfgoed1`).

Konsequenz: `topology_charges` (qa), und damit alles daraus Abgeleitete
(`dxi`, `dgam`, `dalpha`, `alpeeq`), hängt **nicht** von den realen Kernkoordinaten ab,
sondern nur von Elementliste + Bindungsgraph.

### Messung — bit-identisch über Konformere

Testsystem: `test_cases/confs/conf.xyz`, 90 Atome, 353 Konformere.
Verglichen wurden Frame 0, Frame 176 und Frame 352 (maximal auseinanderliegend),
jeweils mit frisch gelöschtem `.topo.json` (kalter Lauf, volle Neuberechnung).

| Größe | max abs. Differenz Frame 0 vs 352 |
|-------|-----------------------------------|
| `topology_charges` (Phase-1 EEQ) | **0.000e+00** (exakt gleich) |
| `dxi` | **0.000e+00** |
| `dgam` | **0.000e+00** |
| `alpeeq`, `eeq_chi`, `eeq_gam`, `eeq_alp`, `eeq_cnf` | 0.000e+00 |
| `hybridization`, `pi_fragments`, `ring_sizes`, `rings` | 0.000e+00 |
| `itag`, `is_metal`, `is_aromatic`, `is_amide_h`, `fraglist`, `qfrag` | 0.000e+00 |
| `coordination_numbers` | 8.475e-03 |

Die einzige abweichende Größe ist CN — und CN liegt **per Design** in
`GFNFFDynamicState`, nicht in `GFNFFTopology` (`gfnff.h:407`), wird also ohnehin
pro Geometrie neu gerechnet.

Diskrete Topologie ebenfalls über alle drei Konformere identisch:
93 Bindungen, 27/93 pi-Bindungsordnungen, 927 BATM-Triples.
Nur die HB-Liste variiert (2177/1902/1950) — die wird pro Schritt neu detektiert
und ist kein Teil des eingefrorenen Blocks.

### Der Cache-Fingerprint ist bereits ein Graph-Schlüssel

`GFNFF::computeTopologyFingerprint()` (`gfnff_method.cpp:2555`) hasht
`N | Elementliste | Bindungsliste` — **keine Koordinaten**.
Alle drei Konformere ergeben denselben Fingerprint `7738153530184087840`.

Damit ist eine Bibliothek vorberechneter Topologien, adressiert über den
Molekülgraphen, vom bestehenden Design bereits getragen.

---

## 2. Was das für Delta-ML bedeutet

### Der Anker ist für Konformere nicht nötig

`DeltaMLMethod` (`deltaml_method.cpp:170-186`) initialisiert die Basismethode aus
einer gespeicherten Ankergeometrie, um "Topologie-Rauschen" aus Delta zu halten
(begründet in `deltaml_method.h:80-85`). Die Messung zeigt: für Konformere
desselben Moleküls **gibt es dieses Rauschen nicht** — die Topologie ist identisch,
mit oder ohne Anker.

Der Anker schützt real gegen genau ein Ereignis: eine Änderung des **Bindungsgraphen**
(Bindung gebrochen/geknüpft). Das ist diskret und selten, nicht kontinuierlich.

### Ein generisches Modell ist damit wohldefiniert

Für einen molekülübergreifenden Satz bekommt jedes Molekül seine eigene Topologie.
Das ist kein Anker-Artefakt, sondern eine wohldefinierte Funktion seines Graphen.
Delta = E_ref − E_GFNFF ist folglich pro Molekül eindeutig, solange der detektierte
Bindungsgraph der chemisch richtige ist.

Das ursprünglich als blockierend eingestufte Hindernis entfällt damit.

### Das verbleibende Risiko ist ein anderes

Nicht der Anker, sondern **Bindungsgraph-Diskontinuität**: liegt ein Kontakt nahe
am Detektions-Cutoff (`rthr = 1.25`, `gfnff_ini2.f90:361-419`), kippt der Graph, und
Delta springt unstetig. Das ist testbar und sollte vor einem generischen Modell
quantifiziert werden — es ersetzt die Anker-Frage als offenen Punkt.

Nicht entkräftet bleiben die übrigen Einwände gegen ein fertiges Universalmodell:
Elementzahl-Baseline wird zum Hauptterm, Modellgröße (M > 12000 hart gesperrt,
Normalengleichungen M x M), und der unbelegte Nutzen molekülübergreifender Modelle.
Ein **Startmodell zum Weiterlernen** für `-activelearn` bleibt die realistischere Form.

---

## 3. Was das für ein schnelles GFN-FF bedeutet

### Die Kostenverteilung ist extrem schief

Gemessen mit `-verbosity 2` (Timing-Tabelle aus `printGFNFFEnergyReport`):

| System | Atome | Topologie | Parametergen. | Setup gesamt | **Energie-Call** | Verhältnis |
|--------|-------|-----------|---------------|--------------|------------------|------------|
| caffeine | 24 | 1.86 ms | 3.00 ms | 4.9 ms | **0.13 ms** | **38x** |
| complex | 231 | 4.06 ms | 11.00 ms | 15.1 ms | **1.64 ms** | **9x** |

Für kleine Organik — also genau das Regime eines generischen Trainingssatzes
(QM9-artig, 10-30 Atome) — sind rund **97 % der Kosten Setup, 3 % Energie**.

### Der bestehende Cache greift nur zu einem Teil

Kalt vs. warm (`.topo.json`-Treffer), caffeine:

| | topo | param | Energie |
|---|------|-------|---------|
| kalt | 1.2 ms | 2.0 ms | 0.08 ms |
| warm (Cache-Treffer) | 0.8 ms | 1.0 ms | 0.07 ms |

Der Cache überspringt nur Phase-1 EEQ (~1.4 ms von 3.2 ms). Bindungsdetektion,
Hückel, Ringperzeption, BATM-Liste und Dispersions-Paare laufen weiterhin.
Selbst warm bleibt das Verhältnis Setup:Energie bei ~26x.

### Hebel

Ein vollständig eingefrorener Topologie-Pfad — der gesamte `GFNFFTopology`-Block
**plus** der daraus abgeleitete Parametersatz (Bonds/Angles/Torsions/Inversions,
`generateGFNFFParameterSet`) — würde die Kosten pro Struktur für kleine Organik
um schätzungsweise **20-40x** senken. Das ist deutlich mehr als die 15 %, die
`gfnff-fast` heute in MD erreicht, weil `gfnff-fast` an einer anderen Stelle ansetzt
(pro-Schritt-EEQ/CN einfrieren, nicht Setup überspringen).

Anwendungsfall ist nicht MD, sondern **Datensatzerzeugung**: viele Strukturen,
jede mit eigenem Setup. Genau der Delta-ML-Trainingsloop.

---

## 4. Bereits vorhandene Bausteine

| Baustein | Ort | Was es tut |
|----------|-----|------------|
| `GFNFFTopology` / `GFNFFDynamicState` | `gfnff.h:337` / `:407` | Struktureller Split statisch/dynamisch, seit März 2026 |
| `topology_mode = "constant"` | `gfnff_method.cpp:899` | Topologie nie neu rechnen |
| `.topo.json`-Cache + Fingerprint | `gfnff_method.cpp:2555, 9813` | Graph-adressierte Persistenz von Phase-1 |
| `static_charges` / `static_cn` / `static_all` | `docs/GFNFF_STATIC_WP1_FROZEN_STATE.md` | Pro-Schritt-EEQ/CN einfrieren |
| `-method gfnff-fast` | `docs/GFNFF_FAST_WP.md` | Preset aus den beiden static-Flags |
| `-forced_bonds` (externe Bindungsliste) | `gfnff_method.cpp:965` | Bindungsgraph von außen vorgeben (aus polymerbuild) |

Der letzte Punkt ist bemerkenswert: es gibt bereits einen Weg, den Bindungsgraphen
zu **injizieren** statt zu detektieren. Damit wäre die Diskontinuitäts-Frage aus
Abschnitt 2 im Prinzip adressierbar — ein einmal festgelegter Graph kippt nicht mehr.

---

## 5. Vorgeschlagene nächste Schritte

Nach Nutzen pro Aufwand:

1. **Diskontinuität quantifizieren.** Über einen MD- oder Konformer-Satz zählen,
   wie oft der Bindungsgraph-Fingerprint wechselt, und die Delta-Sprunghöhe an
   solchen Wechseln messen. Das ist die offene Frage, die den Anker ersetzt.
   Billig: Fingerprint pro Frame loggen.

2. **Setup-Freeze über Strukturen** (`topology_mode=constant` + Parametersatz-Reuse
   über den Fingerprint als Schlüssel). Direkter 20-40x-Hebel für die
   Datensatzerzeugung. Genauigkeit ist per Konstruktion exakt, solange der
   Fingerprint gleich bleibt — kein Approximationsprofil nötig, anders als bei
   `gfnff-fast`.

3. **Anker-Anforderung in Delta-ML lockern.** Der Anker kann auf seine tatsächliche
   Funktion reduziert werden: Fingerprint-Vergleich statt Geometrie-Reinitialisierung.
   Passt Fingerprint, ist kein Anker nötig; passt er nicht, ist die Warnung
   berechtigt und schärfer als heute.

---

## 6. Grenzen dieser Analyse

- **Ein Molekül** für die Konformer-Stabilität (90 Atome, C/H/N/O, neutral).
  Nicht geprüft: geladene Systeme, Metalle, Systeme nahe am Bindungs-Cutoff.
- **Keine MD-Frames** getestet, nur optimierte Konformere. MD-Strukturen sind
  stärker verzerrt; die Bindungsgraph-Stabilität dort ist offen (siehe Schritt 1).
- Der Fingerprint hasht die Bindungsliste über **Atomindizes**, ist also
  reihenfolgeabhängig und kein kanonischer Graph-Hash. Für eine molekülübergreifende
  Topologie-Bibliothek wäre eine Kanonisierung nötig — dasselbe Molekül in anderer
  Atomreihenfolge ergibt heute einen anderen Schlüssel.
- Timing auf einer Maschine, 1 Thread, Sequential-Backend. Die Verhältnisse sind
  robust, die Absolutwerte nicht.
- Die 20-40x-Schätzung in Abschnitt 3 ist aus der Kostenaufteilung extrapoliert,
  **nicht gemessen** — es existiert kein Pfad, der den Parametersatz über
  Strukturen hinweg wiederverwendet.

Menschliche Produktionsprüfung steht aus.
