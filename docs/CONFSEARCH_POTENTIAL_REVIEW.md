# ConfSearch gegen GOAT: Potential-Review der WEKLQ-Kampagne (Stand 2026-10-01)

> 🤖 AI-generierte Auswertung, nur lesend (keine Läufe gestartet, kein src geändert). Jede Energie
> unten ist an der Quelle aus den xyz-Dateien nachgerechnet, der Tautomer-Check der Bestwerte ist
> jeweils gegen `WEKLQ/input.xyz` gerechnet. Ausnahmen sind ausdrücklich als „Journal, nicht
> nachgerechnet“ markiert. Bezug: GOAT-Referenz **−161.6578433336 Eh** (Kommentarzeile von
> `WEKLQ/goat.xyz`), kJ/mol = (E − Ref) × 2625.5. Negativ heißt: tiefer als GOAT.
> Human production testing pending.

## 0. Kurzbriefing (Entscheidungsvorlage)

**Empfehlung:** Vor jeder neuen Zugmenge sollte der Keim-Zufall repariert werden. Seit dem
Race-Fix `17a103e4` (2.9.) hängen die MD-Startgeschwindigkeiten nur noch von `-seed` und der
Startgeometrie ab. In den GOAT-Test-Armen G1–G3 starteten **57–65 % aller MDs von genau dem Keim
der Vorwiederholung**, bei gleicher Temperatur also mit **bitgleichen Anfangsgeschwindigkeiten**.
Die Änderung ist eine Codezeile und kostet keine Rechenzeit. Danach sollte der billigste Weg
repliziert werden: das Hybrid-Protokoll H1 mit zwei weiteren Seeds. Für einen Anspruch „GOAT
geschlagen“ fehlt noch eine eigene GOAT-Rechnung mit protokolliertem Budget.

**Begriffe:**
- **Keim**: Struktur, von der eine MD-Wiederholung startet.
- **Wiederholung (`repeat`)**: erneuter MD-Durchgang derselben Temperaturstufe mit den bis dahin besten Keimen.
- **beschränkt / frei**: optimiert mit bzw. ohne die Zwangsbedingungen `-hold_polar_h` auf den polaren X–H-Bindungen. Beschränkte Energien sind obere Schranken.
- **Release**: der Schlussschritt `-release_polar_h_rank`. Er optimiert die Top 100 frei nach und prüft dabei die Topologie.
- **Tautomer**: Ein Proton ist an ein anderes Schweratom gewandert, es ist also ein anderes Molekül.
- **core-h**: Laufzeit × Threads. Das ist eine Obergrenze, denn die G-Arme liefen überbucht.
- **Race-Fix**: Commit `17a103e4`. Vorher kam die Vielfalt zwischen den Walkern aus einem Thread-Datenrennen, nicht aus `-seed`.

| Rang | Hebel | erwarteter Nutzen | Kosten | billigstes Prüfexperiment |
|---|---|---|---|---|
| 1 | Zufall je Wiederholung variieren (`seed ^ hash(cycle, rep)`) | **Vermutung**. Indiz: vor dem Race-Fix lagen beide repeat-5-Läufe bei ≤ −11 kJ/mol (n = 2), danach alle fünf bei ≥ −2.2 (n = 5, p ≈ 0.10, konfundiert) | 1 Codezeile, 0 Rechenzeit | Diversitätstest: G3-Keim mit G3-Pool, 10 MDs mit gleichem gegen 10 mit variiertem Zufall, danach Quench; ~20–40 core-h |
| 2 | Hybrid H1 (gfnff-MD, gfn2-Relax/Rang) als Standardweg | H1 kam auf −7.8 kJ/mol bei 406 core-h, das sind **0.45×** eines G-Arms (n = 1). Pro Treffer 406 gegen 1690–2720 core-h. **Vermutung** | 2 Seeds × ~400 core-h | 2 H1-Replikate; vorher Release nachträglich auf die alten Hybride (0829: +0.6 beschränkt), einige core-h |
| 3 | Mehrere kurze Läufe statt eines langen (Portfolio) | 1 von 4 repeat-2-Läufen unter GOAT bei ~100 core-h (n = 4, vor dem Fix). **Vermutung** | 4 × 100 core-h | 4 Seeds mit repeat 2 und `-explore_md_time 1950` (physikalisch gleich lang wie früher 1000) |
| 4 | Neuheitsregel: tiefen Nachbarn auch bei gleichem H-Brücken-Muster behalten | **belegt n = 1**: G1 verwarf in Zyklus 3 einen gültigen Vorschlag bei **−2.42** (12.9 kJ/mol unter dem laufenden Bestwert). 1 Fall unter 4090 Vorschlägen aus 6 Läufen | trivial | Regel ab z. B. > 5 kJ/mol Gewinn ohne Musterbedingung; retrospektiver Replay ist gemacht (Abschnitt B) |
| 5 | Kosten Schluss-Dedup und Bias-Pool | reine Ersparnis: Schlussschritt 10–14 % der Wanduhr (G), **37 % (77 h)** im seed_rank-20-Lauf | mittel | Schluss-ConfScan auf dem G1-Ensemble mit Energie-Vorfilter timen; Top 100 müssen identisch bleiben |

**Die drei wichtigsten Unsicherheiten:**
1. **Budget und Vollständigkeit von GOAT sind unbekannt.** In `WEKLQ/` liegt keine ORCA-Ausgabe. Laut Vault war der GOAT-Lauf nicht abgeschlossen. Belegt ist deshalb nur „tiefer als die beste bekannte GOAT-Struktur“, nicht „schneller oder besser als GOAT bei gleichem Budget“.
2. **Die Stichproben sind klein und die Seed-Streuung ist groß** (~10 kJ/mol zwischen identisch konfigurierten Armen). Jeder Gruppenvergleich unten hat n ≤ 5. Der Vorher/Nachher-Unterschied zum Race-Fix fällt mit dem ConfGen-Umbau vom 4./5.9. zusammen.
3. **Die MD-Uhr seit dem master-Merge vom 30.9.** Früher lief eine nominale `-explore_md_time 2000` physikalisch 3903 fs. Der laufende Lauf (`input.confsearch.20260930_235013`, Seed 1) hat bei gleichem Parameter also die **halbe physikalische MD-Länge**. Dazu ist COM-Entfernung aus (`-remove_com_mode 0 -remove_com_motion 0`). Er ist deshalb **kein Replikat** von G1–G3.

---

## A. Alle Läufe und Arme (nachgerechnet)

Methode: `scratchpad/review/scan.py` (Session 8056c0d5) liest alle 81 BMT-Verzeichnisse unter
`WEKLQ/`. Gewertet werden alle optimierten Dateien mit gfn2-Skala. Nicht gewertet werden
`s1_*` (MD-Schnappschüsse), `failed`, `dropped`, `topo_rejected`, unoptimierte `proposals.xyz`
und `s0_initial.*.xyz`. Die tiefste Struktur mit korrekter H-Zuordnung gilt als Bestwert. Der
Check vergleicht Kovalenzbindungen (Faktor 1.25), den H-Partner und die Schweratombindungen
gegen `input.xyz`; die GOAT-Topologie ist mit `input.xyz` identisch. „Tautomer über dem
Bestwert“ zählt die tieferen, aber ungültigen Strukturen desselben Laufs. Das Budget ist
Wanduhr (erste bis letzte Datei) × Threads (aus der fish-History bzw. `goat_test.sh`) und
`n_opt` (Frames in `s2_relax` + `proposals.opt` + `s5_refine*.opt`).
**Alle Läufe liegen vor dem 30.9.; die physikalische MD-Länge war 1.95 × nominal.**

| Lauf (Verzeichnis) | Datum | Methode MD/Rang | Kernparameter | Wand h × Thr = core-h | n_opt | Bestwert gültig kJ/mol | Skala | Tautomere darüber |
|---|---|---|---|---|---|---|---|---|
| Juli/Anf. Aug. (~25 Läufe, `20260726_*` … `20260805_*`) | 26.7.–5.8. | gfnff/gfn2 bzw. gfn2/gfn2, meist 1–2 Thr | repeat 1–3, ohne hold_polar_h | 0.3–55 h | 11–63 594 | Median ≈ +58, bester **+26.9** (`0804_233419`) | frei | ja, in 13 Läufen war der rohe Bestwert ein Tautomer (bis **−17.5** in `0803_164833`/`0804_095937`) |
| `0806_172652` | 6.8. | gfnff/gfn2 (Relax gfnff) | repeat 3, 2 Thr | 64 × 2 = 127 | 14 231 | +28.5 | frei | 3 |
| `0809_090038` | 9.8. | dito + surface_feedback | repeat 3 | 27 × 2 = 54 | 8 396 | +28.6 | frei | 2 |
| `0810_095046` | 10.8. | dito | repeat 5 | 257 × 2 = 514 | 24 888 | +38.4 | frei | 12 (roh +3.9) |
| `0822_194926` | 22.8. | gfn2/gfn2, ohne hold_polar_h | repeat 5, 10 Thr | 25 × 10 = 253 | 5 958 | **+2.9** | frei | 6 (Top 6 des Pools Tautomere, roh −10.3) |
| `0823_211047` **Rekord** | 23.8. | gfn2/gfn2, hold_polar_h | repeat 5, 600→300 K, 10 Thr | 109 × 10 = 1088 | 10 677 | **−19.77** | frei¹ | 0 |
| `0824_233218` | 24.8. | gfn2/gfn2 | 600→500, time 1000, repeat 10, 6 Thr | 24 × 6 = 145 | 2 512 | +16.2 | beschr. | 0 |
| `0825_235049` | 25.8. | gfn2/gfn2 | 600→400, time 1000, repeat 2 | 12 × 6 = 73 | 1 308 | +12.4 | beschr. | 0 |
| `0826_004137` | 26.8. | gfnff/gfn2 (Relax gfn2) | repeat 5, 600→450, 2 Thr | 178 × 2 = 356 | 5 981 | +15.1 | beschr. | **11 (Top 11 des Pools Tautomere, dienten 12 Wdh. als Keime)** |
| `0827_082258` | 27.8. | gfn2/gfn2, ohne ConfGen | repeat 2, time 1000, 600→300, 4 Thr | 23 × 4 = **93** | 1 690 | **−1.29** | beschr. | 0 |
| `0827_082327` | 27.8. | dito + ConfGen | dito | 34 × 4 = 134 | 2 184 | +6.4 | beschr. | 3 |
| `0828_154027` seed_bias_penalty 20 | 28.8. | gfn2/gfn2 | repeat 5, 10 Thr (per SIGHUP abgebrochen) | 124 × 10 = 1236 | 9 433 | **−11.27**² | frei¹ | 50 |
| `0828_180536` | 28.8. | gfn2/gfn2 + ConfGen | repeat 2, time 1000 | 27 × 4 = 107 | 1 626 | +11.4 | beschr. | 0 |
| `hybrid_explore_refine` `0829_231850` | 29.8. | gfnff/gfn2, refine_md | repeat 2, time 1000, 4 Thr | 37 × 4 = 147 | 2 735 | **+0.6** | beschr. | 0 |
| `hybrid_v2` `0901_003512` | 1.9. | gfnff-Explore, gfn2-Refine-MD | repeat 2 | 43 × 4 = 171 | 3 468 | +8.8 | beschr. | 0 |
| Golden `0902_205547` | 2.9. | gfn2/gfn2 | repeat 5, 600→450, 2000 fs, 14 Thr | 50 × 14 = 706 | 6 430 | +4.41 (frei +2.44, Journal, nicht nachgerechnet) | beschr. | 0 |
| arm0/1/2/3 (550 K) | 4./5.9. | gfn2/gfn2 | repeat 3, 550 K, 5 Thr | je ~7 × 5 = 36 | 654–666 | +38.3 / +41.5 / +38.2 / +28.5 | beschr. | 0 / 2 / 0 / 0 |
| **G1** s101 | 5.9. | gfn2/gfn2 | Golden-Protokoll, 8 Thr | 111 × 8 = 891 | 6 638 | **−2.19** (verworfen: −2.42, s. B4) | frei | 1 (Release-Ausgabe) |
| **G2** s102 | 5.9. | dito | dito | 119 × 8 = 950 | 7 322 | +4.97 | frei | 1 |
| **G3** s103 | 5.9. | dito | dito | 104 × 8 = 831 | 6 589 | **−1.12** | frei | 0 |
| **H1** s101 | 5.9. | gfnff/gfn2, refine_md | repeat 2, 600→300, 1000 fs, 6 Thr | 68 × 6 = **406** | 3 816 | **−7.77** | frei | 5 in der Release-Ausgabe (bis −32.1) |
| gfnff-Trichter `0907_223811` | 7.9. | gfnff/gfnff | seed_rank 20, 4 Thr | 90 × 4 = 359 (+30 gfn2) | 14 974 | +39.9 (Journal; die gfn2-Nachoptimierungen lagen im gelöschten Scratchpad 2f925b59, **nicht nachrechenbar**) | frei | 62 von 2011 (Journal) |
| `0916_104429` seed_rank 20 | 16.9. | gfn2/gfn2 | Golden-Protokoll mit seed_rank 20, seed 42, 10 Thr (Threads angenommen, Kommando in der History bei 10:25) | 207 × 10 = 2067 | 13 749 | **+11.09** | frei | 1 (roh −5.7, vom Release-Guard korrekt verworfen) |
| laufend `0930_235013` | 30.9. | gfn2/gfn2 | G-Protokoll, Seed 1, 10 Thr, COM aus, **neue MD-Uhr** | läuft (Zyklus 1) | – | – | – | – |

¹ Vor dem 4.9. liefen die ConfGen-Optimierungen ohne Zwangsbedingungen. Der Rekord (−19.77) und
der Strafwert (−11.27) sind freie Vorschlagsenergien mit gültiger Topologie. Der Rekord war der
frei nachoptimierte beschränkte Elter bei −15.73.
² Das Journal/Projekte nennt −8.9 (beschränkt). Nachgerechnet liegt im selben Lauf ein
**gültiger** Vorschlag `nci_from_1_d1` bei −11.27 (frei). Er wurde fünfmal in Zyklus 4 r5 und 5 r1
gebaut und nie ins Ensemble übernommen. Damit war auch dieser Lauf unter GOAT.

**Verdichtung nach Gruppe** (Bestwert gültig, kJ/mol; beschränkte Werte sind obere Schranken):

| Gruppe | n | Median | min | max | unter GOAT | core-h Median | core-h pro Treffer |
|---|---|---|---|---|---|---|---|
| gfn2/gfn2 lang (≥ 20 h) | 12 | +3.7 | −19.8 | +16.2 | 5/12 | 768 | 1700 |
| davon repeat 5 **nach** Race-Fix (Golden, G1–3, sr20) | 5 | +4.4 | −2.2 | +11.1 | 2/5 | 891 | 2722 |
| davon repeat 5 **vor** Race-Fix (Rekord, Strafe) | 2 | −15.5 | −19.8 | −11.3 | 2/2 | 1162 | 1162 |
| gfn2 repeat 2, time 1000 nominal | 4 | +8.9 | −1.3 | +12.4 | 1/4 | 100 | 407 |
| Hybrid gfnff-MD / gfn2-Relax (0826, 0829, v2, H1) | 4 | +4.7 | −7.8 | +15.1 | 1/4 | 264 | 1080 |
| Hybrid mit gfnff-Relax (Aug.) | 3 | +28.6 | +28.5 | +38.4 | 0/3 | 127 | – |
| 550-K-Arme | 4 | +38.3 | +28.5 | +41.5 | 0/4 | 36 | – |

Zeitpunkt des Endbestwerts (beschränkte Skala, `traj.py`): Golden nach 10.6 von 43.2 h
(**75 % der Suchzeit ohne Fortschritt**), G3 nach 41 von 94 h (56 %), sr20 nach 24 von 129 h
(**81 %**), H1 nach 40 von 67 h (40 %), Rekord nach 40 von 105 h (61 %). G1 und G2 fanden ihren
Bestwert erst spät (85 bzw. 91 h, 14 bzw. 13 % ohne Fortschritt).

---

## B. Misserfolge und Sackgassen, jeweils mit der belegenden Messung

1. **gfnff als Rangfläche.** Pearson gfnff gegen gfn2-SP beträgt +0.10 (n = 1077). Die 8 gfn2-besten Strukturen stehen auf der gfnff-Liste auf Rang 273–1631 von 1949. Innerhalb eines Zyklus liegt r bei −0.32/−0.46 (n = 13/14). Quelle: Journal 7.–11.9. und DUAL_METHOD; hier nicht neu gerechnet. Die Aussage ist gut belegt.
2. **Reine gfnff-Suche.** Von 2011 gfnff-Strukturen landete frei auf gfn2 keine unter +25 und keine unter GOAT; die beste lag bei +39.9 (Journal). Die Rohdaten sind nicht mehr vorhanden. Bestätigt wird das Bild von den drei nachgerechneten August-Läufen mit gfnff-Relax (+28.5/+28.6/+38.4). Dazu wachsen die Bias-Kosten quadratisch: 3.6 s pro MD-Schritt bei 23 400 Hügeln und 451 Permutationen, 43 GB RSS, weil jeder Aufruf den Pool vollständig kopiert.
3. **Rg-Bias/Rg-Flutung verworfen.** Die Rg-Verteilung der Ablagen wandert nicht (Median 5.46–5.51 Å in allen vier Armen). Die Bestwerte der Arme (+38.3/+41.5/+38.2/+28.5) liegen innerhalb der Streuung der beiden Kontrollen (10 kJ/mol). Jeder Arm hat n = 1.
4. **Neuheitsregel verwirft einen tieferen Nachbarn (neu gefunden).** G1, `cycle03_T500K_r5`: `nci_from_2_d1_unconverged` bei **−2.42 kJ/mol** (beschränkt, unkonvergiert, also obere Schranke) mit gültiger Topologie. Der nächste Nachbar im Templat-Ensemble lag 0.85 Å entfernt bei +13.3. Der Vorschlag steht nur in `proposals.opt`, nicht in `.new`, und das Log meldet den Zyklusbesten bei −161.653506. Laut `confgen.cpp:3539-3556` muss ein tieferer Nachbar zusätzlich ≥ 1 abweichende H-Brücke haben (`new_pattern_min 1`). Wahrscheinlicher Grund der Ablehnung ist deshalb ein gleiches Muster (Hamming 0); nicht direkt geloggt. Der Lauf fand −2.19 erst einen Zyklus später. Häufigkeit: **1 Fall bei 4090 Vorschlägen aus 6 Läufen nach dem 2.9.** (Skript `lost.py`). Im Strafe-Lauf (vor dem Fix) wurde ebenso ein gültiger Vorschlag bei −11.27 nie übernommen. Dort lagen allerdings Skalen gemischt vor (frei gegen beschränkt), eine saubere Attribution ist nicht möglich.
5. **seed_rank 20 (gfn2).** Ergebnis +11.09 bei 2067 core-h und 13 749 Optimierungen, gegen G1–G3 mit −2.2/−1.1/+5.0 bei 830–950 core-h. n = 1 gegen n = 3: doppeltes Budget, schlechter als jeder G-Arm. Nicht im Journal (Lauf nach dem 11.9.). Dazu **77 h Schluss-Dedup**, 37 % der Wanduhr (Zeitstempel-Lücke 21.9. 20:04 bis 25.9. 01:25; ein Log fehlt, die Zuordnung ist plausibel, aber nicht belegt).
6. **seed_bias_penalty 20.** Ergebnis −11.3 (nachgerechnet) gegen −19.8 im Rekordlauf bei 1.5× Optimierungen. n = 1 gegen 1. Einordnung „Breite statt Tiefe“: das Signal liegt knapp über der Seed-Streuung, ist aber nicht belastbar.
7. **relax_polish_rank war ein Messfehler** (Restraint-Energie). Belegt mit Kontrolle, n = 6: beschränkt nachoptimiert 0.00 kJ/mol. Zurückgenommen.
8. **„Release als Hebel für die Suche“ ist widerlegt (neu gemessen).** Die freie Nachoptimierung verschiebt fast konstant: Median 2.5–3.2 kJ/mol, P90 4.8–6.2, max 8–17. Die Rangfolge bleibt dabei erhalten: Spearman beschränkt gegen frei **0.94–0.97**. Der freie Beste war in 4 von 5 Läufen der beschränkte Rang 1, einmal Rang 2. Höchstens 1 der freien Top 10 lag außerhalb der beschränkten Top 10 (n = 5 Läufe × 100, `rel.py`). Konsequenz: Die Sub-GOAT-Werte von G1/G3 bestehen fast ganz aus diesem Offset. Beschränkt standen sie bei +0.2/+0.3; der **Suchvorsprung ist ≈ 0**.
9. **Tautomere.** In freien Läufen sowie im Hybrid 0826 vor dem ConfGen-Fix waren die tiefsten Pool-Strukturen Tautomere. Beispiel 0826: H107 ging von O57 nach N25, das ist ein Zwitterion. Die Top 11 waren betroffen und dienten 12 Wiederholungen als Keime; gültig blieb am Ende nur +15.1. Im aktuellen Code ist der Pool sauber: H1 hat Tautomere nur in der verworfenen Release-Ausgabe. Die Release-Tautomerrate hängt aber an der Explore-Methode: H1 38 %, G 2–10 %. Freie gfn2-Tautomere liegen bis ~25 kJ/mol unter dem richtigen Molekül (H1: −32.1 gegen −7.8).
10. **Zugsätze (Route, NCI, Klammer, Kompaktierung, Crossover).** Alle Vorschlagstypen liegen auf gfn2 bei einem Median von +80…120 kJ/mol. Die Kompaktierung scheitert an Clashes (0 von 45). Quelle: Journal 4./5.9., n = 1 pro Variante. Einzelne Treffer gibt es dennoch: H1 −4.5/−4.6 (Journal).
11. **Stagnations-Abbruch wäre schädlich.** Lange trockene Strecken gingen späteren Sprüngen voraus: G2 hatte 8 trockene Wiederholungen, danach +23 → +18 → +11; der Strafe-Lauf 9 trockene, danach +6 → −7 (nachgerechnet per Wiederholung). Eine Abbruchregel mit ≤ 8 hätte beide Sprünge gekostet.

---

## C. Wo liegt Vorsprung vor GOAT?

**Die GOAT-Referenz, an der Quelle geprüft.** Sie besteht aus **einer** Struktur mit
−161.6578433336 Eh (`goat.xyz`). Ihre Topologie ist identisch mit `input.xyz`. Curcumas
gfn2-Einzelpunkt weicht 0.02 kJ/mol ab (Journal, nicht neu gerechnet). **Budget, Laufzeit und
Ensemble von GOAT fehlen im Verzeichnis.** Laut Vault war der GOAT-Lauf nicht abgeschlossen.

**Belegt** (Topologie geprüft, gleiche gfn2-Skala, freie Energien sind mit GOAT direkt
vergleichbar):
- Die gfn2-Fläche hat mindestens eine gültige Mulde **19.8 kJ/mol unter GOAT** (Rekord, n = 1). curcuma kann sie finden.
- **6 Läufe** erreichten gültige Strukturen unter GOAT: Rekord −19.8, Strafe −11.3, H1 −7.8, G1 −2.2 (−2.4 verworfen), 0827a −1.3 (beschränkt), G3 −1.1. Das sind 6 von 20 gfn2-gerankten Läufen.
- Der Topologie-Schutz im Schlussschritt funktioniert: 0 Tautomere unter den gemeldeten Endwerten von G1–G3/H1/sr20, obwohl 2–38 % der Release-Ausgaben Tautomere waren.

**Nur Vermutung:**
- „Sicher schlagen“: Das G-Protokoll nach dem Fix liegt bei 2/5 (einschließlich sr20) bzw. 2/4 (ohne sr20). Bei unabhängigen Läufen mit p ≈ 0.4–0.5 ergäben drei parallele Läufe 78–88 % Chance auf mindestens einen Treffer. n = 5 ist dafür zu klein.
- „Hybrid ist billiger“: H1 lag bei 406 core-h pro Treffer (n = 1). Die übrigen Hybride verfehlten, liefen aber auf älterem Code (0826 mit Tautomer-Keimen) oder hatten keinen Release-Schritt (0829 +0.6 beschränkt; mit dem gemessenen Median-Offset von 2.5 kJ/mol **vermutlich** unter GOAT, nicht gerechnet).
- „Vor dem Race-Fix war die Suche tiefer“: 2/2 gegen 0/5 unter −2.2 kJ/mol. Konfundiert mit dem ConfGen-Umbau (Zwangsbedingungen, Route, Neuheitsregel) und der Leiter bis 300 K. Beide Rekordwerte fielen aber schon bei 450 K. Golden (Race-Fix ja, ConfGen alt) kam auf +2.4 frei, was eher auf den Zufall/Race-Fix deutet als auf ConfGen. n = 1.

**Wo wir hinken:**
- **Reproduzierbarkeit.** Identische Konfiguration streut ~10 kJ/mol zwischen den Seeds (G: −2.2/−1.1/+5.0; 550-K-Kontrollen: 10 kJ/mol).
- **Suchvorsprung ≈ 0 auf der beschränkten Skala** (B8): G1/G3 endeten dort bei +0.2/+0.3.
- **Budget.** 830–950 core-h pro G-Arm, davon 10–14 % im Schluss-Dedup. Ob GOAT billiger ist, lässt sich nicht beurteilen (Unsicherheit 1).
- **Die GOAT-Mulde selbst wurde nie gefunden** (kleinster RMSD 1.93 Å; eine GOAT-Brücke kommt im Ensemble nie vor). Das Ensemble ist bei GOATs Mulde also unvollständig, bei tieferen Mulden nicht.

**Widerlegungsfrage zum Gesamtfazit** („curcuma kann GOAT unterbieten, aber nicht zuverlässig“):
Widerlegt wäre es, wenn GOAT bei protokolliertem, vergleichbarem Budget selbst ≤ −2 kJ/mol
erreicht (dann kein Vorsprung) oder wenn drei weitere G-Seeds alle unter GOAT landen (dann
zuverlässig). **Keine der beiden Messungen liegt vor.**

---

## D. Hebel im Detail (Rangfolge wie im Briefing)

### 1. Zufall je Wiederholung variieren
- **Befund (Code):** `SimpleMD::seedRandomEngine()` (`simplemd.cpp:1271-1283`) seedet mit `m_seed XOR hash(Startgeometrie)`. ConfSearch gibt jedem Kind dasselbe `-seed`. Gleicher Keim bedeutet deshalb gleiche Normalverteilungs-Ziehungen; bei anderer Temperatur werden sie nur mit √T skaliert. Der stochastische Thermostat hängt am selben Generator.
- **Befund (Log, n = 3 G-Arme × 190 Keime):** 108/117/123 Keime (57/62/65 %) waren identisch mit einem Keim der Vorwiederholung, davon 93/94/105 bei gleicher Temperatur. In H1 waren es 55/130 (42 %). In G3 blieb Keim 1 zehn Wiederholungen in Folge dieselbe Struktur (Zyklus 3 r1 bis Zyklus 4 r5), in genau der Phase ohne Fortschritt.
- **Befund (Geometrie):** Walker 0 aus demselben Keim, G3 Zyklus 3: RMSD zwischen aufeinanderfolgenden Wiederholungen 0.7–1.5 Å bei 210–390 fs, 3.5–5.9 Å ab ~800 fs (4 Paare). Die Trajektorien laufen also nur über den gewachsenen Bias-Pool auseinander. Der Vergleich mit dem Rekordlauf (1.3 Å schon bei 40 fs) ist schwach: Es ist dort nicht geprüft, ob der Keim derselbe war.
- **Nutzen: Vermutung**, keine Zahl belegbar. Obere Indizgrenze: repeat-5 vor dem Fix Median −15.5 gegen +4.4 danach, konfundiert.
- **Kosten:** eine Zeile, z. B. `m_seed ^ hash(cycle, rep)` beim Kind-Aufruf. Die Reproduzierbarkeit bleibt erhalten.
- **Billigstes Experiment:** G3-Keim `cycle02_T550K_r5_w0_i6535_t60` mit dem G3-Pool (Restart-Stand Zyklus 3). 10 MDs bei 500 K mit identischem gegen 10 mit variiertem Zufall, 2000 fs, jeweils Quench der Schnappschüsse. Ablesen: Zahl verschiedener Minima innerhalb 20 kJ/mol und Bestwert; ~20–40 core-h. **Widerlegt, wenn** der identische Zufall gleich viele verschiedene Minima liefert. Danach folgt erst der teure A/B (3 Seeds × G-Protokoll).

### 2. Hybrid H1 als Standardweg prüfen
- **Nutzen:** −7.8 kJ/mol bei 406 core-h und 3816 Optimierungen, gegen G mit 830–950 core-h und 6600–7300 Optimierungen. n = 1, **Vermutung**.
- **Risiko:** 38 % Tautomere im Release. Der Guard fängt sie, sie verlieren aber den mittleren Release-Gewinn von ~2.5 kJ/mol. Seit dem master-Merge ist die gfnff-MD eine andere: C6 wird pro Schritt aufgefrischt, und die MD-Uhr ist korrigiert (`-explore_md_time 1000` war früher 1950 fs physikalisch).
- **Billigstes Experiment:** (a) Release nachträglich (Top 100 frei plus Topologie-Check) auf `hybrid_explore_refine_20260830` (+0.6 beschränkt), `hybrid_v2` (+8.8) und die vier gfn2-repeat-2-Läufe; einige core-h. Kippt 0829 unter GOAT, steht der Hybrid bei 2/4. (b) Zwei H1-Replikate mit `-explore_md_time 1950`, ~800 core-h. **Widerlegt, wenn** beide über GOAT enden.

### 3. Portfolio kurzer Läufe
- **Nutzen:** 0827a erreichte −1.3 (beschränkt) bei **93 core-h** (1 von 4 kurzen Läufen). Vier kurze Läufe kosten etwa die Hälfte eines G-Arms. **Vermutung**, und alle vier liefen vor dem Race-Fix.
- **Gegenargument (liegt vor):** Die Rekordkette brauchte die Wiederholungen r3–r5 derselben Stufe (CONFSEARCH_LEVERS). Kurze Läufe schneiden genau diesen Mechanismus ab.
- **Experiment:** 4 Seeds mit repeat 2 und `-explore_md_time 1950`, je 4 Threads. Vergleich gegen die G-Arme pro core-h.

### 4. Neuheitsregel für tiefe Nachbarn
- **Nutzen belegt, n = 1:** −2.42 statt −2.19, einen Zyklus früher. Ein früherer, tieferer Keim ist dabei nicht eingerechnet. Häufigkeit 1 auf 4090 Vorschläge.
- **Kosten:** trivial. Vorschlag: `gain > 5 kJ/mol` akzeptiert ohne Musterbedingung, darunter gilt die bestehende Regel.
- **Experiment:** Unit-Test mit genau diesem Vorschlag und seinem Templat. Retrospektiv liegt kein weiterer Fall vor.

### 5. Budget-Lecks
- **Schluss-Dedup:** 12.6/13.9/10.1 h in G1–G3 (10–14 %), 77 h in sr20 (37 %). Es ist quadratisch in der Ensemblegröße, bei H1 mit 401 Konformeren nur 1 h.
- **Bias-Pool:** Vollkopie pro Auswertung (`shared_bias_pool.cpp:16-20`), dazu Permutationen ohne Deckel. Für gfn2-Läufe ist der Anteil nicht gemessen, für gfnff dominiert er (3.6 s/Schritt).
- **Nutzen:** 10–37 % Wanduhr, frei für zusätzliche Wiederholungen. Das Ergebnis ändert sich nur, wenn Hebel 1 greift.
- **Experiment:** Schluss-ConfScan auf `g1_s101/.../input.cumulative.opt.xyz` mit Energie-Vorfilter (z. B. 50 kJ/mol) timen; Top 100 müssen gleich bleiben.

### Nicht weiterverfolgen (Messung liegt vor)
gfnff als Rangfläche oder Selektor; reine gfnff-Suche; Rg-Bias und Rg-Flutung; seed_rank 20;
relax_polish; Release als Rangverbesserung (B8); Stagnations-Abbruch mit ≤ 8 Wiederholungen
(B11); Kompaktierungszug.

---

## Anhang: Messmethode und Fallen

- **Skripte** (Session-Scratchpad, temporär): `…/8056c0d5-…/scratchpad/review/`
  - `scan.py`: alle Läufe, Bestwert und Topologie
  - `top.py`: Top-N pro Lauf mit Herkunft
  - `lost.py`: verworfene Vorschläge unter dem laufenden Bestwert
  - `traj.py`: laufender Bestwert pro Wiederholung
  - `rel.py`: beschränkte gegen freie Rangfolge
  - `same_seed.py`: Walker-0-Trajektorien bei gleichem Keim
  - `hdiff.py`: welches H gewandert ist
- **Fallen berücksichtigt:**
  - Energieskalen: frei und beschränkt stehen getrennt in der Spalte „Skala“.
  - Budget: core-h und n_opt sind angegeben. Die G-Arme liefen überbucht (30 Threads auf 16 Kerne), ihre core-h sind also überschätzt.
  - Provenienz: aus Namen gelesen, nicht per RMSD.
  - Tautomere: jeder Bestwert ist geprüft.
- **Nicht nachgerechnet:** Golden frei (+2.44), gfnff-Trichter (+39.9, 1949/62), Korrelationen von DUAL_METHOD/Journal, GOAT-Einzelpunkt 0.02 kJ/mol.
- **Threads des Laufs `0916_104429`** sind aus der History bei 10:25 übernommen (20 Walker-Verzeichnisse bestätigen seed_rank 20); das exakte Startkommando fehlt.

## E. Nachtrag 2026-10-02: Messungen zu den Hebeln (AI-erzeugt, maschinell ausgewertet, nicht von Hand geprüft)

Alle Zahlen aus Ausgabedateien nachgerechnet; gfnff-Tests auf dem WEKLQ-Peptid (107 Atome), Rohdaten unter `/m2scratch/conrad/Projekte/2026/confsearch/seedvary_{pilot,paired}/`.

- **Hebel 1 (`-seed_vary`, neu, Default aus):** Mechanismus belegt, Nutzen nicht. Ohne den Schalter wiederholen sich 21-55 % der MD-Starts bitgleich (gleicher Seed + gleiche Startgeometrie, `simplemd.cpp` `seedRandomEngine()`), mit Schalter 0 %. Gepaarter Test (gfnff, 4 Temperaturstufen, je 3 Wiederholungen, ConfGen aus): tiefste Energie vary minus fixed, n = 7 Paare, Median +2.7 kJ/mol, SD der Paardifferenzen 19.8, Vorzeichentest 5:2, p = 0.45 -- nicht unterscheidbar, Teststärke gering (für 10 kJ/mol Effekt bräuchte es grob 30 Paare, Überschlag). Nur gfnff, ein Molekül, nicht gfn2. Der Pilot mit 1 Zyklus (n = 3) war ebenfalls nicht aussagekräftig; das Maß "Zahl verschiedener Minima" taugt nicht (fast jede relaxierte Struktur zählt als verschieden).
- **Seed-Streuung ist größer als jeder getestete Schalter:** `fixed`-Arm, tiefste Energie je Lauf 0.9 bis 61.6 kJ/mol über dem Besten (n = 8, Median 26.9). Ein einzelner Lauf ist zu großen Teilen ein Zufallstreffer; das stützt Hebel 3 (mehrere kürzere Läufe) mehr als Hebel 1. Das ist eine Vermutung auf n = 8 und einem Molekül.
- **Hebel 5 (Schluss-Dedup):** an Logs bestätigt für große Pools (G1-G3: 9.7-11.6 % der Wanduhr, n = 3 Läufe; sr20 37 % nur aus Dateizeitstempeln, Obergrenze). Bei kleinen Pools < 2 %. Dominant ist überall die Optimierung (46-91 %).
- **Release-Schritt (`-release_polar_h_rank`) funktioniert wie gedacht:** Im Lauf `20260930_235013` (gfn2, `-hold_polar_h`) gewannen 78 von 100 freigegebenen Strukturen im Median 3.33 kJ/mol (max 8.75). Bei 22 von 100 wanderte das Proton beim Lösen der Restraints; diese behalten ihre restringierte Geometrie und stehen nicht als Tautomere im Ergebnis (`cumulative.opt.accepted.xyz`: 1002 Strukturen, 0 gewanderte H). Die Rohdatei `s7_release.gfn2.opt.xyz` enthält diese 22 Tautomere (tiefster -161.660638 Eh, 16.2 kJ/mol unter der Suche) und ist KEIN Ergebnis. Rang 100 reicht: die 100 tiefsten Strukturen umspannen 48.6 kJ/mol, der größte Release-Gewinn war 8.75.
- **Ergebnis des Laufs `20260930_235013`** (gfn2, 600-450 K, je 5 Wiederholungen, `-seed 1`, COM-Entfernung aus, Laufzeit ca. 29.5 h, 10 Threads): tiefste Struktur im akzeptierten Ensemble nach Release -161.656039 Eh (0 gewanderte H), das sind -92.9 kJ/mol gegenüber der ersten Struktur (-161.620669 Eh) und 4.1 kJ/mol unter dem Wert vor dem Release (-161.654465 Eh). Ihr Vorsprung ist isoliert: im akzeptierten Ensemble liegt nur 1 Struktur innerhalb 10 kJ/mol des Besten. Der Bestwert wurde einmal gefunden (Zyklus 4, Wiederholung 2), davor Plateau seit Zyklus 2, Wiederholung 2 (n = 1 Lauf). Kein GOAT-Vergleich: es fehlt eine GOAT-Referenz mit protokolliertem Budget auf gleicher Skala.
- **Beobachtungen zu Zeit und Budget:** Relaxation einer Wiederholung ca. 1-1.5 h (155-207 Strukturen à 107 Atome gfn2); die kürzeren MD-Schritte seit dem MD-Uhr-Fix (`ef462fcf`) liefern pro Wiederholung gut 10-mal mehr Schnappschüsse als früher (207 gegenüber 20 in Zyklus 1).
- **Hebel 3 (k kurze Läufe gegen einen langen): NICHT beantwortet.** Versuch: gfnff, 4 Temperaturstufen, kurze Läufe `-repeat 3` (fixed_s21-28, n = 8, alle fertig: Wanduhr 27-62 min, 1.4-3.1 CPU-h, tiefste Energie rel. zum Besten Median 26.0, SD 18.7 kJ/mol) gegen lange Läufe `-repeat 6` (Seeds 31-34). Alle 4 langen Läufe wurden am 150-min-Limit abgebrochen (7.5 CPU-h, noch in der letzten Temperaturstufe 450 K, Wiederholung 1-4 von 6; tiefste Energie Median 30.6, SD 9.4, n = 4); die zweite Charge (Seeds 35-38) wurde gestoppt. Befund: **Die Laufzeit wächst überproportional mit den Wiederholungen** (lang/kurz mindestens 4.5x statt 2x, wegen des wachsenden Pools; Last der Maschine 23-35 gegenüber ca. 10 bei den kurzen, Wanduhrzeiten daher nicht 1:1 vergleichbar). Minima aus je 2 kurzen Läufen: Median 20.2 (alle 28 Paare, nicht unabhängig) bzw. 15.8 (4 disjunkte Paare) bei ca. 3.4 CPU-h gegen 30.6 bei 7.5 CPU-h der unfertigen langen -- Vermutung (n = 4 gegen 8, unfertig, lastverfälscht), kein Beleg. Widerlegung wäre: fertige lange Läufe mit gleichem oder tieferem Median bei vergleichbarer Last; diese Messung liegt nicht vor (0 fertige lange Läufe).
  - **Neu-Design (nicht gestartet):** Limit 6-7 h je Lauf (Schätzung 3.5-5 h für `-repeat 6`), bei Last > 20 höchstens 2 parallel, erst Pilot mit 2 langen Läufen zur Zeitmessung; oder billiger `-repeat 4` gegen 2x `-repeat 2` (unter ca. 2 h); gepaarte Seeds, n >= 8 je Arm, Budget als "CPU-h bis zum selben Punkt des Temperaturplans" definieren.
- **Hebel 3, billiges Design (3.10.2026, gfnff, gepaart n = 8, alle 24 Läufe fertig, 0 gewanderte H):** `-repeat 4` einzeln (L, Seeds 41-48, 93-130 min, mittel 4.68 CPU-h) gegen Minimum aus 2 Läufen `-repeat 2` (S, Seeds 51-66, mittel 0.85 CPU-h je Lauf). Budget nicht gleich: 2 S kosten im Median 0.41 des L-Budgets. Paardifferenz (S-Minimum minus L, kJ/mol): Median -2.8, Mittel -3.5, SD 14.9, Vorzeichen 4:4 -- nicht unterscheidbar, aber S braucht nur ca. 40 % der CPU-h. Einzelläufe rel. zum Besten: L Median 19.8 (SD 6.9, n = 8), S Median 25.3 (SD 15.0, n = 16): L ist zuverlässiger, S hat die größere Streuung; das globale Minimum kam aus einem S-Lauf. Explorativ (nachträglich gruppiert, n = 4): Best-of-4-S (je 2 Paare, ca. 3.4 CPU-h) lag in 4 von 4 Gruppen vor Best-of-2-L (ca. 9.4 CPU-h). Aussage: **bei gfnff, einem Molekül und Maschinenlast 9-21 (Wanduhren nicht 1:1) liefert die kurze Variante pro CPU-h mindestens gleich viel, vermutlich mehr; ein Effekt bei gleichem Budget ist bei n = 8 nicht belegt.** Für gfn2 nicht geprüft.
- **gfn2-Portfolio 4 kurze Läufe gegen 1 langen (5.-6.10.2026, AI-erzeugt, maschinell ausgewertet):** 4 unabhängige Läufe `-repeat 2` (Seeds 101-104, je 3 Threads, 26-30 h Wanduhr, ca. 340 Thread-h in Summe), sonst Parameter des Referenzlaufs `20260930_235013` (`-repeat 5`, 10 Threads, 29.5 h, ca. 295 Thread-h; Budgets also etwa gleich, nicht ein Bruchteil). Endergebnis je Lauf aus `cumulative.opt.accepted.xyz` nach Release, jeweils 0 gewanderte H: s102 -161.656119 Eh (n = 390; **-0.21 kJ/mol gegenüber dem Referenzwert -161.656039 Eh**, also gleichauf), s101 -161.652252 (n = 442; +9.9), s103 -161.649927 (n = 605; +16.1), s104 -161.649759 (n = 606; +16.5). Median +13.0, Mittel +10.6, Spanne 16.7 kJ/mol. Aussage: **das Beste aus 4 kurzen Läufen erreicht bei etwa gleichem Budget den langen Lauf, ein einzelner kurzer Lauf liegt im Mittel 10 kJ/mol darüber.** n = 4 gegen 1: ob der lange Lauf bei anderen Seeds ebenfalls streut, ist unbekannt (gfnff: Seed-Streuung SD ca. 20 kJ/mol), die Aussage ist eine Vermutung. Auffällig und nicht untersucht: je Lauf 13000-19000 `[ERROR]`-Zeilen "native SCF did not converge within 151 iterations (max|dq| ~ 4)" in den Logs (vermutlich unphysische Geometrien bei Optimierungen und ConfGen-Vorschlägen; im ConfGen-Batch von Wiederholung 1 "5-10 failed" bei 19-22 Strukturen, Text im Log abgeschnitten); die Läufe liefen trotzdem durch. Ob der lange Lauf diese Ausfälle ebenfalls hatte, ist unbekannt (Ausgabe ging an das Terminal).
