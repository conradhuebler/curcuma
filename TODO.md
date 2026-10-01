# Curcuma Development TODO List

**Stand**: 2026-10-01. Die frühere Fassung (1380 Zeilen, mit Messprotokollen und erledigten Untersuchungen) steht
unverändert in [docs/TODO_ARCHIVE_2026-10.md](docs/TODO_ARCHIVE_2026-10.md); die Detailmessungen zu den Einträgen unten
sind dort unter der jeweiligen Überschrift zu finden.

Regeln für diese Datei: ein Eintrag hat höchstens drei Zeilen (Was, warum es offen ist, wo die Details stehen). Erledigtes
gehört in `AIChangelog.md` bzw. das Archiv, Fehlerberichte in [docs/KNOWN_ISSUES_ARCHIVE.md](docs/KNOWN_ISSUES_ARCHIVE.md).
"Geprüft" heißt: am angegebenen Datum gegen den Quelltext gelesen, nicht neu gemessen. Dateien und Zeilen gelten für diesen
Stand.

## 1. Korrektheit (offen)

- **GFN-FF-Energie hängt von der Atomnummerierung ab** (911 von 2557 Strukturen, bis 237 kcal/mol bei geladenen Mehrfragment-Systemen;
  in pprcht und curcuma). Referenzeigenschaft, kein Portierungsfehler. [docs/REV_GFNFF_TODO.md](docs/REV_GFNFF_TODO.md) #12, Known Issue #34.
- **ROCm-GFN-FF ohne Coulomb-Term** als Voreinstellung seit `ab6e3f5e` (F-1, Umgehung `-gfnff.gpu_coulomb_implicit false`); GFN1/GFN2 auf ROCm
  ignorieren `-scf_mixed_precision false` (G2-13); Energiesummen nehmen wave32 an (F-17). Keine ROCm-Hardware, nur dokumentiert:
  [docs/MULTI_GPU_GAPS.md](docs/MULTI_GPU_GAPS.md).
- **UFF- und QMDFF-Gradient falsch** (Eintrag vom 2026-09-27: UFF-Winkel ohne den Faktor -sin(theta); Fix vorbereitet, Betreiber: "erstmal nicht").
  Nicht erneut geprüft; Details im Archiv unter "UFF- und QMDFF-Gradient falsch".
- **`-seed` wirkt nicht auf die Anfangsgeschwindigkeiten der MD.** Geprüft 2026-10-01: `SimpleMD::InitVelocities()` (`simplemd.cpp:1141`)
  zieht aus `static std::default_random_engine generator` (Z. 1143) mit festem Standardseed. Ein Fix ändert jede MD-Startgeschwindigkeit.
- **`shouldUpdateHBXB()`: RMSD-Formel ohne Faktor `sqrt(natoms)`** (`gfnff_method.cpp:3058`, geprüft 2026-10-01). Bei 7320 Atomen
  ~13,6 Å mittlere Verschiebung bis zum Auslösen. Aus der Fortran-Referenz übernommen (`gfnff_ini2.f90:717`); vor einer Änderung gegen
  MOR41/GMTKN55 prüfen und das Abweichen von der Referenz entscheiden.
- **D4-ATM-C6 auf den Tripeln bleibt auf dem Setup-Wert eingefroren.** Ohne Wirkung, solange `dispersion_atm` aus ist (Voreinstellung).
  Ein Paar, das erst später in 60 Bohr eintritt, wurde nicht durch einen Test erzeugt. [docs/GFNFF_PAIR_LIST_REFRESH.md](docs/GFNFF_PAIR_LIST_REFRESH.md).
- **Optimierer**: (a) der eigene L-BFGS bleibt bei Koffein/GFN-FF ab Schritt 33 stehen (die Stillstandserkennung `-opt.stall_steps` beendet
  den Lauf inzwischen, Ursache nur per Codelesung); (b) `mixture2`: kein Abstieg an einer N-H···O=C-Brücke, nicht eingegrenzt;
  (c) der `auto`-Optimierer kann bei hochsymmetrischen Systemen mit nur totalsymmetrischer Mode abbrechen (Td-CH4, gfn1; `-opt.optimizer lbfgs` geht).
- **GFN-FF-NVE-Drift auf `complex`** (231 Atome): 6,1e-3 Eh über 1 ps bei dt 0,5 fs, unabhängig vom EEQ-Löser. Nicht untersucht.
- **SIGSEGV am Ursprung nicht gefunden** (2026-09-24; ASan, compute-sanitizer, gdb ohne Treffer). Eintrag im Archiv.
- **`DC13/c20bowl`**: Coulomb -6,29 kcal/mol gegen die Referenz, Ladungen ~6x zu groß, nicht eingegrenzt (Known Issue #14).
- **Reproduzierbarkeit**: GPU läuft nicht bitgleich (`atomicAdd`-Reihenfolge); CPU nicht bitgleich über verschiedene Threadzahlen
  (`FFWorkspace` partitioniert nach Threadzahl). Betreiber 2026-09-27: später. Erster Schritt: Repulsion als Sammel-Kernel hinter
  einem Schalter, messen auf A4500 und H200. Archiv: "CPU/GPU-Trajektorien-Divergenz" und "Bitgleichheit".
- **Verbosity** ist ein geteilter Static und lässt sich über `CxxThreadPool`-Worker nicht sauber scopen (Known Issue #3).
- **ConfSearch**: Wide-Hill-MTD-Aufblähung innerhalb eines Laufs; `threads=1` im Legacy-Modus mit Fortschrittsbalken (SIGSEGV nach der
  Optimierung, Eintrag von 2025, nicht erneut geprüft). [docs/CONFSEARCH_ROADMAP.md](docs/CONFSEARCH_ROADMAP.md).
- **Windows**: `-DUSE_PORTABLE_MATH=ON` nicht auf einer echten Windows-Maschine bestätigt ([docs/PORTABLE_ERF.md](docs/PORTABLE_ERF.md)).
- **Ladungsplatzierung** (Known Issue #31) nicht getestet für periodische Systeme, GPU-Laufzeit, viele Fragmente in Kontakt und lange MD mit
  `frag_charge_s_max` > 1.

## 2. Leistung

- **GFN2 auf der GPU, Eigenlöser = 73 % des Laufs** (polymer_2x, 1x A4500). Offen:
  (a) `-scf_pseudo_diag` als Standard einschalten (Betreiberentscheidung, Messungen sprechen dafür);
  (b) `-scf_pseudo_diag_fp64` auf einer GPU mit vollem FP64 (H200) messen;
  (c) CPU-Pfad für `-scf_pseudo_diag` in `XTB::solveEigen` (Schätzung 20 bis 25 %, nicht implementiert, kein Gewinn unter ~1500 Basisfunktionen);
  (d) Teilspektrum: der residente Loop übergibt `n_eig = 0`, Grund unbekannt, schließt sich derzeit mit der verteilten FP64-Lösung aus;
  (e) Schwellen `gpu_eigensolver_min_nao` (4000) und `gpu_density_min_nao` nachmessen;
  (f) `purify`/`lobpcg` gibt es nur auf der CPU, auf der GPU ist es der einzige Hebel, der die Skalierung ändert;
  (g) Setup und Integrale (32 s) nie im Detail profiliert.
  **`scf_fp32_stall_patience` nicht ändern**, bis die Fixpunktfrage (Schritt 1 weicht um 2,7e-6 Eh ab) geklärt ist.
  Zahlen und Begründungen: Archiv, Abschnitt "GPU-SCF GFN1/GFN2"; [docs/GFN2_GPU_COST_PLAN.md](docs/GFN2_GPU_COST_PLAN.md).
- **GPU-Gradient GFN1/GFN2**: Punkte 1 bis 6 umgesetzt; offen: der H0-Paar-Kernel (454 ms) ist einzelkartig, und `densityPatternDistributed`
  kopiert die C-Spaltenscheiben bei jedem Aufruf neu. [docs/SQM_PERFORMANCE.md](docs/SQM_PERFORMANCE.md).
- **GFN-FF-MD, EEQ-Löser = 69 % des Schritts** (polymer_2x, 758,7 von 1095 ms, 16 CPU-Threads). Ob das nahe am Erreichbaren liegt, ist unbelegt.
  Billiger erster Schritt: `-verbosity 2` und die Iterationszahl der projizierten PCG lesen.
- **GFN-FF Static-Mode-Initiative** (Mai 2026): die PARAMs `static_charges`/`static_cn` (WP-S1) sind im Code (`gfnff.h:323`), und der
  Changelog nennt WP-S1 bis S3 als umgesetzt. WP-S4 (Validierungssuite, 32-Konfigurationen-Matrix) nicht geprüft.
  [docs/GFNFF_STATIC_WP4_VALIDATION_SUITE.md](docs/GFNFF_STATIC_WP4_VALIDATION_SUITE.md).
- **Speicher für Systeme > 1000 Atome** (Molecule-Datenstruktur, Distanzmatrix-Cache): geplant, nicht begonnen.

## 3. Tests und Validierung

- `cli_simplemd_08/09` (Essigsäuredimer, CSVR, dt 1 fs) gegen das gemergte Binary mit der MD-Uhr-Korrektur (`ef462fcf`) erneut laufen lassen (Known Issue #32).
- `test_cg_potentials` wird gebaut, hat aber keinen ctest-Eintrag; der CG-Beispielaufruf mit VTF-Eingabe endet mit "Failed to initialize ForceField
  engine" (`docs/archive/CG.md`, geprüft 2026-10-01).
- `molecule_comprehensive` verweist auf das Target `test_molecule`, dessen Quelle `src/core/test_molecule.cpp` nicht im Baum liegt.
- Wissenschaftliche Validierung der CLI-Tests ausbauen (RMSD-Toleranzen, Energiekonvergenz); Muster für absichtlich fehlschlagende
  Tests (`03_invalid_method` in `curcumaopt`, `rmsd`, `confscan`); Performance-Benchmarks für Regressionserkennung.
- ConfScan: Accept/Reject-Meldungen bei Standard-Verbosity nicht sichtbar (Eintrag von 2025, nicht erneut geprüft).

## 4. Betreiber-Prüfung offen (🤖/⚙️, nur Sie vergeben ✅)

- CODATA-2018-Einheiten (`191cebe3`; Rest 0,0036 kcal/mol gegen pprcht auf polymer_2x nicht verfolgt), Known Issue #35.
- `eeq_refactor_eps_bohr` Voreinstellung 0 (Faktor-Cache aus) und `nonbonded_skin_bohr` 2,0; Ladungsplatzierung `ensemble` (#31);
  Stillstandserkennung `-opt.stall_steps` 20; `-scf_guess fragments`.

## 5. Funktionen und Refactoring (niedrige Priorität)

- **CG Phase 6**: winkelabhängige Ellipsoid-Energie (`calculateEffectiveDistance`, `calculateCGPairEnergy`, Rotationszüge in Casino). Phase 1 bis 5 sind umgesetzt.
- SimpleMD: Physik des Wandpotentials prüfen; RMSD-Strategy-Pattern Phase 3; erweiterte ConfSearch-Algorithmen; bessere Trajektorienanalyse;
  ConfSearch mit GPU und mehreren Threads (nur noch relevant, solange der GPU-Geräte-Pool inaktiv ist, Korrektur 2026-09-28).
- **Molecule-Refactoring** (Phasen 2 bis 6): XYZ-Kommentar-Parser vereinheitlichen, granulare Caches, O(1)-Fragmentzugriff, `ElementType`-Enum,
  SOA/AOS. Plan: `src/core/REFACTORING_ROADMAP.md`, Formate: `src/core/XYZ_COMMENT_FORMATS.md`.
- **Native QM** (Stand November 2025, nicht geprüft): GFN2-Parametererweiterung, PM3-Elementumfang (F u. a.), Validierungsmoleküle, Doku
  "Wann GFN2, GFN1 oder PM3". Der Elementumfang von PM3 im Code ist nicht verifiziert.
- **Build-System** (Stand November 2025, nicht geprüft): bedingte Kompilierung, fünf Build-Varianten, nur der Standardbuild war damals fehlerfrei.
  Seit 2026-09-30 baut `make` in `release/` wieder mit Exit 0; die übrigen Varianten vor einer Wiederaufnahme neu messen.

## 6. Aus der alten Liste entfernt (überholt oder erledigt)

| Eintrag | Grund |
|---|---|
| "cgfnff Parameter Generation Bug", "Missing Real GFN-FF Parameters" | GFN-FF ist vollständig, die Dateien und das Problem existieren nicht mehr |
| "`gpu_strict` fehlt" | umgesetzt, `src/core/gpu_fallback.h`, [docs/GPU_TUNING.md](docs/GPU_TUNING.md) Abschnitt 1 |
| CG Phase 1 bis 5 | umgesetzt |
| EEQ-Warmstart über die q-Loop-Durchgänge | umgesetzt und nachgemessen 2026-09-18 |
| `make` in `release/` scheitert an CUDA-Unittests | behoben 2026-09-30 |
| Dispersionspaarliste, Repulsionspaarliste | behoben, Known Issue #33 |
| Optimierer: Liniensuchfehlschlag als Konvergenz, Abbruch ohne Struktur, SCF-Unsinn als Ergebnis | behoben 2026-09-27 |
| Hessian: SCF-Schwelle erreicht die Worker nicht | behoben 2026-09-29 |
| Cholesky-Faktor-Cache in der MD nicht exakt | Cache standardmäßig aus (`eeq_refactor_eps_bohr` 0), Betreiber-Prüfung offen (Abschnitt 4) |
