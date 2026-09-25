# Curcuma Development TODO List

**Stand**: 2026-03-13
**Quelle**: Extrahiert aus allen CLAUDE.md-Dateien
**Test Status**: 26/26 CLI tests passing (100%) ✅

---

## 🔴 KRITISCH / HÖCHSTE PRIORITÄT

*No critical blockers remaining!* All SimpleMD issues resolved.

---

## 🟢 COARSE GRAINING DEVELOPMENT

### ✅ CG Integration - Phases 1-5 (COMPLETE - Oct/Nov 2025)
- Phase 1: Molecule helper functions (isCGSystem, hasMixedSystem, getCGAtoms)
- Phase 2: JSON parameter loading in forcefield.cpp
- Phase 3: VTF reader/writer in formats.h
- Phase 4: Testing & validation (unit tests + CLI test)
- Phase 5: SimpleMD CG integration (PBC, 10x timestep scaling, VTF output, orientational infrastructure)

### 🔵 CG Potentials - Phase 6: Ellipsoidal Extensions (OPTIONAL - LOWEST PRIORITY)
- **Status**: 🟡 PREPARED
- **Priority**: 🔵 LOW - After SimpleMD integration
- **Task**: Implement angle-dependent energy for ellipsoidal particles
- **Betroffene Dateien**: src/core/energy_calculators/ff_methods/cg_potentials.cpp, src/capabilities/casino.cpp
- **Aufwand**: ~4-5 h (after Phase 5 complete)
- **Current Status**: All infrastructure prepared (rotation matrices, ellipsoid detection, fallback energy)
- **Remaining**:
  - Complete `calculateEffectiveDistance()` for ellipsoid-ellipsoid interactions
  - Implement orientation-dependent potential in `calculateCGPairEnergy()`
  - Activate rotational moves in Casino
  - Test ellipsoidal shape calculations

---

## 🟡 TESTING & VALIDATION (test_cases/CLAUDE.md:296-303)

### Scientific Validation Enhancement
- **Status**: ⏳ PENDING
- **Task**: Erweitere wissenschaftliche Validierung für alle Tests (RMSD tolerances, energy convergence)
- **Dateien**: test_cases/cli/test_utils.sh, test_cases/cli/*/run_test.sh
- **Note**: All tests now passing (26/26), ready for enhanced validation

### Expected Failure Pattern für invalid_method Tests
- **Status**: ⏳ PENDING
- **Task**: Implementiere Pattern für Tests, die bewusst fehlschlagen sollen
- **Tests betroffen**:
  - test_cases/cli/curcumaopt/03_invalid_method/run_test.sh
  - test_cases/cli/rmsd/03_invalid_method/run_test.sh
  - test_cases/cli/confscan/03_invalid_rmsd_method/run_test.sh
- **Dateien**: test_cases/cli/test_utils.sh

### Performance Benchmarks & Regression Tests
- **Status**: ⏳ PLANNED
- **Task**: Füge Performance-Benchmarks hinzu für Regressions-Detection
- **Dateien**: test_cases/CMakeLists.txt, test_cases/cli/test_utils.sh

### test_molecule.cpp Extension
- **Status**: ⏳ PENDING
- **Task**: Erweitere für geplantes Molecule SOA/AOS Refactoring (Phase 2-6)
- **Dateien**: test_cases/test_molecule.cpp
- **Abhängigkeit**: Molecule Refactoring Phases müssen geplant sein

---

## 🟢 CORE DEVELOPMENT (src/core/)

### GFN-FF Static-Mode Performance Initiative (4-WP-Bundle, 2026-05)
- **Status**: 🤖 GEPLANT — basiert auf Analyse "C++ schneller als xtb"
- **Ziel**: Polymer 1000 NVE ICX 130 s → ~63 s (40 % schneller als xtb-Fortran)
- **WP-S1**: [Frozen CN/Charges](docs/GFNFF_STATIC_WP1_FROZEN_STATE.md) — 3 PARAMs, ~150 LoC, ~1 Tag
- **WP-S2**: [MD-Diagnostics JSONL](docs/GFNFF_STATIC_WP2_MD_DIAGNOSTICS.md) — Streaming-Dump, ~1 Tag
- **WP-S3**: [EEQ-Cutoff Auto-Default](docs/GFNFF_STATIC_WP3_EEQ_CUTOFF_DEFAULT.md) — 12 ms/Step save, ~0.5+1 Tag
- **WP-S4**: [Validation-Suite + Re-Capture](docs/GFNFF_STATIC_WP4_VALIDATION_SUITE.md) — 32-Konfig-Matrix, ~1.5 Tage
- **Reihenfolge**: S1 → S2 (parallel zu S3) → S4
- **Risiko**: S1/S3 sind physikalische Approximationen — S4 ist die Sicherheits-Klammer

### cgfnff Parameter Generation Bug
- **Status**: ❌ OPEN
- **Problem**: Parameter generation creates null JSON values
- **Betroffene Dateien**: src/core/energy_calculators/qm_methods/gfnff.cpp
- **Verweis**: CLAUDE.md - Known Issues, src/core/CLAUDE.md:104
- **Konsequenz**: Native GFN-FF nicht einsatzbereit

### Missing Real GFN-FF Parameters
- **Status**: ❌ OPEN
- **Problem**: Placeholders statt echter physikalischer Parameter
- **Betroffene Dateien**: src/core/energy_calculators/qm_methods/gfnff.cpp
- **Verweis**: src/core/CLAUDE.md:105
- **Abhängigkeit**: Benötigt theoretische Implementierung oder externe Daten

### Unit System Migration (CODATA-2018)
- **Status**: ⏳ IN PROGRESS
- **Task**: Replace hardcoded constants mit `CurcumaUnit` namespace functions
- **Betroffene Dateien**: Multiple legacy files mit hardcoded constants
- **Verweis**: src/core/CLAUDE.md:113, CLAUDE.md:296
- **Gewinn**: Centralized, documented, CODATA-2018 compliant constants

### GFN-FF EEQ warm start across the q-loop passes — DONE, entry was stale (re-measured Sep 18, 2026)
- The warm start is in place and works: polymer_2x with the topology cache deleted, 7320 atoms /
  1502 fragments, pass 1 needs **75 PCG iterations / 877 ms**, pass 2 **0 iterations / 26 ms** —
  it starts from the pass-1 solution and is converged immediately. The "32 + 71 iterations,
  ~2.6 s" this entry used to quote no longer happens.
- What is left of the EEQ cost is the FIRST solve (75 iterations, 877 ms of the GFN-FF setup).
  Reducing that needs a better preconditioner or a better initial guess, not a warm start.
- Reproduce: `rm <name>.topo.json && CURCUMA_GFNFF_PROFILE=1 curcuma -sp <xyz> -method gfnff
  -verbosity 2` and read the two `projected PCG converged in N iterations` lines.

### Memory Optimization for Large Systems (>1000 atoms)
- **Status**: ⏳ PLANNED
- **Task**: Optimize Molecule data structure and distance matrix caching
- **Betroffene Dateien**: src/core/molecule.cpp/h, src/core/energy_calculators/ff_methods/forcefield.cpp
- **Verweis**: src/core/CLAUDE.md:106, CLAUDE.md - Performance Notes
- **Performance Impact**: Critical for large molecular systems

### GPU-SCF GFN1/GFN2: der Eigenloeser ist 73 % des Laufs (2026-09)
- **Status**: ⏳ OFFEN, Zerlegung gemessen, kein Punkt umgesetzt
- **Messung** (polymer_2x, 7320 Atome, nao 15444, GFN2, 1x RTX A4500, `-sp -gradient` 254 s,
  `CURCUMA_GPU_PROFILE=1`, Geraete festgenagelt):

  | Phase | s | Anteil |
  |---|---:|---:|
  | **Eigenloeser gesamt** | **174.2** | **73 %** |
  | Dichte + Populationen | 39.8 | 17 % |
  | Setup + Integrale | 32.0 | 13 % |
  | Gradient (nach der Optimierung) | 6.3 | 2.6 % |

  Einzeln: `eig FP32 syevd` 82.7 s / 11 Aufrufe, `eig FP64 syevd` 49.3 s / 1,
  `eig FP64 reduce` 21.4, `eig FP64 back-transform` 10.5, `eig FP32 copy+reduce` 6.2,
  `eig FP32 back-transform` 4.1.

**1. Weniger FP64-Iterationen — GEMESSEN, gemischtes Bild, Fixpunkt-Frage offen**
- Eine FP64-Iteration kostet 49.3 + 21.4 + 10.5 = **81.2 s**, eine FP32-Iteration **8.5 s**
  (Faktor 9.6). Der Einzelpunkt braucht 11 FP32 + 1 FP64.
- **Im Einzelpunkt (kalter SCF) ist `scf_fp32_stall_patience` wirkungslos.** Sweep 1/2/3/5
  liefert exakt dieselbe Aufteilung (11 FP32 + 1 FP64, `stall 0` in allen vier) und dieselbe
  Zeit auf die Sekunde. Der Waechter feuert dort nicht — FP32 macht echte Arbeit bis zur
  Konvergenz. `-scf_fp32_threshold 1e-6` (enger als die FP32-Rauschgrenze) kostet **+38 s**
  (3 zusaetzliche FP32-Iterationen, dann greift der Waechter doch): Beleg fuer die Rauschgrenze,
  keine brauchbare Schraube.
- **Im MD-Schritt (extrapolierter Start) feuert der Waechter, und `patience=1` spart 33 %**
  (325.0 -> 216.1 s, 1 FP64-Runde weniger, 5 FP32-Iterationen weniger). Mechanismus: nach
  `-scf_extrapolation aspc` startet der SCF bereits unterhalb der FP32-Rauschgrenze, FP32
  kann dort nichts mehr leisten, und `patience` bestimmt nur, wie schnell das bemerkt wird.
  `mixed_precision=false` (reines FP64) ist dagegen katastrophal: 716.4 s, 28 FP64-Runden.
- **Offen, nicht kleinreden**: Schritt 0 (kalt, Waechter feuert nicht) ist bit-identisch
  zwischen Vorgabe und `patience=1` — die Geometrie fuer Schritt 1 ist also gleich. Schritt 1
  selbst unterscheidet sich um 2.7e-6 Eh. Ob das derselbe Fixpunkt in einem groesseren
  Konvergenzfenster ist oder ein anderer, ist **noch nicht geklaert** (Test bei
  `-scf_threshold 1e-7` lief methodisch falsch — fehlendes `-dump_frequency 1`, `Final Energy`
  existiert im MD-Log nicht, nur die 6-Nachkommastellen-Tabelle — und wird wiederholt).
  **Bis zur Klaerung `patience` nicht als Vorgabe aendern.**
- **Nur Consumer-Karten.** Auf vollwertigem FP64 ist gemischte Genauigkeit standardmaessig aus
  (`xtb_gpu_method.cpp:703`), dort existiert der Effekt nicht.

**2. Teildiagonalisierung ist implementiert, wird aber nicht benutzt**
- `eigensolveResidentFock` kann das Teilspektrum (`xtb_gpu_context.cu:2943-2950`, AP1,
  `cusolverDnDsyevdx/Ssyevdx`, `il=1..n_eig`). Die Dichte braucht nur die besetzten Spalten,
  `syevd` rechnet alle 15444 Eigenpaare.
- **Der residente Loop uebergibt `n_eig = 0`** (`:4994`), also volles Spektrum. Warum, ist
  **unbekannt** — Versaeumnis oder verworfener Versuch; erst nachsehen, dann messen.
- **Vorbehalt aus dem Code**: die verteilte FP64-Loesung ist auf `!partial` gegattert
  (`:2959`) — Teilspektrum und Mehrkartenloesung schliessen sich derzeit aus. Fuer eine Karte
  waere es ein Gewinn, fuer mehrere muss man waehlen.

**3. Schwellen der vorhandenen Verteilung nachmessen**
- Eigenloeser verteilt: **1.95x** auf 4 PCIe-Karten gemessen. Dichte verteilt: **3.4x** gemessen
  (39758/12 -> 11558/13 Aufrufe).
- Die Gates `gpu_eigensolver_min_nao` (4000) und `gpu_density_min_nao` sind gesetzt, aber
  **nie nachgemessen** worden. Wo genau sich Verteilen lohnt, ist offen.
- Auf einer H200 ist NVLink der einzige Unterschied, der hier nicht simulierbar ist — der
  Datenweg zwischen den Karten ist die gemessene Bremse.

**4. Struktureller Hebel: gar nicht diagonalisieren**
- `purify` und `lobpcg` existieren in curcuma, aber **nur auf der CPU** — im CUDA-Pfad kein
  einziger Treffer (geprueft). Dichtematrix-Purifikation (McWeeny) kommt ohne Diagonalisierung
  aus und nutzt Duennbesetzung.
- **Der einzige Punkt dieser Liste, der die SKALIERUNG aendert** statt der Konstanten, und
  entsprechend der aufwendigste.

**5. Setup und Integrale (32 s) sind nie im Detail profiliert worden**
- Im Einzelpunkt einmalig, in einer **MD pro Schritt**. `integrals: overlap + H0` 11.4 s,
  `setup` 20.6 s — was darin steckt, ist unbekannt.

**6. GFN-FF ist ein voellig anderer Fall**
- Dort sind **93 % des Schritts der EEQ-Loeser** (1345 von 1474 ms auf der GPU, 758.7 von
  1095 ms auf 16 CPU-Kernen). Keiner der Punkte 1-5 greift. Siehe den eigenen TODO-Eintrag.

### GPU-Gradient GFN1/GFN2 — Punkte 1-6 UMGESETZT (2026-09)
- **Status**: alle sechs Punkte implementiert und gemessen. Offen geblieben: der H0-Paar-Kernel
  (454 ms) ist noch einzelkartig, und `densityPatternDistributed` kopiert die C-Spaltenscheiben
  bei jedem Aufruf neu auf die Helferkarten (nur die Musterindizes sind zwischengespeichert).
- **Ergebnis** (polymer_2x, 1x A4500, selbst nachgemessen): `-sp -gradient` **287.1 -> 256 s**,
  Gradientenblock **18747 -> 6499 ms**, W-DGEMM **12748 -> 3309 ms**, `finalize: download P and C`
  **15332 -> 0.0 ms**, Geraetespeicher im Gradienten **12337 -> 8143 MiB**. Energie unveraendert
  (`-11799.19965134 Eh`), Gradient gegen den Referenzbau **1.4e-14** bei Standard-SCF-Schwelle.
  ctest: gpu 200/200, gpu_gradient 24/24, sqm 335/335.
- **Messung vorher** (polymer_2x, 7320 Atome, nao 15444, 1x RTX A4500): `-sp` 249 s,
  `-sp -gradient` **287.1 s** -> der Gradient kostet **38 s**, nicht die 96 s, die eine
  fruehere Fassung dieses Eintrags aus MD-Schritt minus Einzelpunkt gebildet hatte.
  Davon: Geraet **13.6 s** (94 % davon EINE DGEMM, s. Punkt 3), Host-Gradientenarbeit
  ~5.2 s, Host-Rueckholung **19.1 s** (Punkt 1). MD-Schritt 332.4 s; die restlichen 45 s
  sind eine zusaetzliche FP64-Eigenloesung, kein Gradient (docs/SQM_PERFORMANCE.md).
- **Vollstaendige Analyse mit Zeilennummern**: [docs/SQM_PERFORMANCE.md](docs/SQM_PERFORMANCE.md)
  "The GPU gradient at 7320 atoms: where the 96 s go"
- **Reihenfolge ist wichtig** — erst messen, dann bauen; Schritt 3 zuletzt:

1. **ERLEDIGT** — **`-gradient` schaltet einen Sparpfad ab, den der Gradient nicht braucht — GEMESSEN 19.1 s,
   und zwar in JEDEM MD-Schritt** (Download P/C 15.3 s, Potential 1.86 s, Energien 1.67 s,
   Bandenergie 0.28 s; ohne Gradient sind alle vier exakt 0.0 ms). `xtb_native.cpp:1519-1523` faellt bei
   `gradient == true` in `finalize()` und laedt P und C herunter, obwohl
   `xtb_gradient.cpp:862-868` leere Matrizen mit `pc_resident=true` uebergibt. Host-Rueckfall
   ist abgesichert (`xtb_native.cpp:1678-1681`). **~15 Zeilen, eine Datei, Gradient-Mathematik
   unveraendert.** Gate: Gradient muss bit-identisch bleiben.
2. **ERLEDIGT**: `computeGradient` hat jetzt 8 `profMark`-Aufrufe (env-gated, ohne Kosten wenn
   `CURCUMA_GPU_PROFILE` nicht gesetzt ist). Damit sind die Punkte 3-5 beziffert statt geschaetzt.
3. **ERLEDIGT** — **die `W`-DGEMM ist 12.75 s = 94 % des Geraetegradienten.**
   `W = C_occ*diag(2 eps)*C_occ^T` (`:5306`) wird dicht ueber nao^2 gebaut, gelesen wird es nur
   an den 6.1 % gespeicherten Paaren (`:1218-1219`). `k_density_sp` (`:1687`) ist dieselbe
   SDDMM und kann `W` mit Gewicht `2*eps` auf dem Muster bauen — ~16x weniger Flops und ~1.9 GB
   weniger. `ensureDenseDensity` (`:5280`) kostet dagegen **0.0 ms** (P ist nach dem residenten
   finalize schon dicht), die urspruengliche Vermutung dazu war falsch.
4. **ERLEDIGT** (gethreadet ueber `parallelStripes`, Abweichung 2.1e-14) — **Zwei `nat^2`-Schleifen auf einem Host-Kern**: GFN2-Multipol-Wechselwirkungsgradient
   (`xtb_gradient.cpp:894-968`) und CN-Kettenregel (`:970-1000`), 2.68e7 Paare, kein Cutoff,
   nicht gethreadet — waehrend Host-Abschnitt 2b via `parallelStripes` (`:231-234`) sehr wohl
   threadet. Die CN-Kettenregel hat zudem keinen Cutoff, wo die Energie bei 25 Bohr
   abschneidet (`xtb_gpu_context.cu:832`) — inkonsistent **und** langsam.
5. **ERLEDIGT** (Coulomb 233.8 -> 27.4 ms, Faktor 8.5; Repulsions-Cutoff geprueft und BEWUSST
   NICHT gemacht, weil `calcRepulsionEnergy` (`xtb_h0.cpp:369-400`) selbst keinen hat) — **`gexp == 2.0` hart am einzigen Aufrufer** (`:5362`), aber vier FP64 `pow` je Schalenpaar
   in `k_grad_coulomb` (`:1330`); `k_grad_repulsion` (`:1102`) drei `pow` je Atompaar ohne
   Cutoff. Beide in der `for j<i`-Form mit `atomicAdd` — Gather-Umbau steht bereits in
   [docs/SQM_GPU_ROADMAP.md](docs/SQM_GPU_ROADMAP.md):39-46.
6. **ERLEDIGT fuer W** (3308.7 -> 855.2 ms auf 4 Karten, 1-vs-4-GPU-Gradient 2.13e-10 bei
   `-scf_threshold 1e-9`; der H0-Paar-Kernel mit 454 ms ist weiterhin einzelkartig und waere der
   naechste Kandidat) — **Gradient ueber GPUs verteilen**: `k_grad_h0_pulay_sp` ist eine reine Reduktion
   ueber den Paarbereich, Ausgabe nur `grad` (3*nat) + `dEdcn` (nat) = **176 KB** — billiger zu
   verteilen als der Eigenloeser. Skelett existiert in `densityPatternDistributed` (`:3118-3258`).
- **Warnung vor dem naheliegenden Ansatz**: die schalenpaar-blockierte Form des
  Multipol-Gradienten wurde auf der CPU gemessen und **als langsamer verworfen**
  (docs/SQM_PERFORMANCE.md "Blocking the multipole GRADIENT: tried, measured, reverted").
- **Validierung**: `ctest -L gpu_gradient` (`sqm_cuda_gradient_*`, tol 1e-7),
  `scripts/gradient_compare.py`, `gradient_unit_contract`.

### GFN-FF-MD: der EEQ-Loeser ist 69 % des Schritts — Stand explorativ, nicht optimiert (2026-09)
- **Status**: ⏳ OFFEN, Kostenanteil gemessen, Optimalitaet **unbekannt**
- **Messung** (polymer_2x, 7320 Atome, ~1500 Fragmente, CPU 16 Threads, MD dt 1 fs,
  `-md_diagnostics_timing`): Schritt 1095 ms, davon **`eeq_solve` 758.7 ms**. Naechstgroesster
  Posten `d4_gw` 65.7 ms. Das ist Pro-Schritt-Arbeit (Ladungen haengen an der Geometrie),
  **kein** Setup.
- **Was bereits da ist** (nachgelesen, nicht nachgemessen): die Voreinstellung `cholesky`
  routet fuer dieses System automatisch auf projiziertes PCG (`eeq_solver.cpp:1583-1584`,
  `ppcg_auto` bei `nfrag >= 1 && natoms >= 500`), und dieser Pfad hat einen Warmstart aus den
  Ladungen des Vorschritts (`m_ppcg_last_q`, in `solveWithProjectedPCG`). Die 758.7 ms sind
  also **mit** Warmstart.
- **Was NICHT gezeigt ist**: dass dieser Stand nahe am Erreichbaren liegt. Die Auswahl ist eine
  Heuristik mit handgesetzten Schwellen (`eeq_pcg_expected_iters` 30, `eeq_pcg_nfrag_threshold` 4,
  `pcg_large_threshold` 500), und die einzige dokumentierte Skalierungszahl ist
  „polymer/1410: 44 -> 16 ms" aus der PARAM-Beschreibung (`eeq_solver.h:1099`) — ein 5x
  kleineres System. Ob 758.7 ms bei 7320 Atomen gut oder schlecht sind, ist unbelegt.
- **Naechster Schritt, billig**: ein Lauf mit `-verbosity 2` gibt die Iterationszahl aus
  (`[EEQ] projected PCG converged in N iterations (|Pr|=..., nfrag=...)`). Wenige Iterationen
  = der Warmstart traegt und die Kosten liegen im Matrixprodukt; viele = die Projektion auf
  ~1500 Nebenbedingungen dominiert und ist der Hebel.
- **Bereits gemessen und erledigt**: `-eeq_solver.solve_method pcg` ist hier **78x langsamer**
  (59.4 s gegen 761 ms je Loesung), weil einfaches PCG `nfrag+1` Loesungen braucht — im Code
  bei `:1600-1603` vorhergesagt. `ppcg` explizit zu setzen aendert nichts (identischer Pfad,
  bit-identische Energie). Kein Flag-Gewinn abzuholen.
- **Relevanz fuer Hardware**: solange 69 % des Schritts in einem nebenbedingungs-behafteten
  linearen Loeser stecken, entscheidet dessen Implementierung, ob eine schnellere GPU bei
  GFN-FF etwas bringt — nicht die Kraftfeld-Kernel.

### `make` in release/ scheitert: fuenf CUDA-Unittests ohne `USE_CUDA` (2026-09)
- **Status**: ⏳ OFFEN, vorbestehend (nicht von der Gradienten-Instrumentierung verursacht)
- **Symptom**: `cd release && make -j8` endet mit **Exit 2**. Die Hauptziele bauen
  (`curcuma_cuda` 7 %, `curcuma_core` 62 %, `curcuma` 63 %); es scheitern nur
  `test_xtb_cuda_{cn,eeq,gamma,gradient,h0,multipole,overlap,qat}` mit
  „`gpu` in Namensbereich `curcuma::xtb` bezeichnet keinen Typ".
- **Ursache**: `namespace gpu` in `cuda/xtb_gpu_context.h:28` steht hinter `#ifdef USE_CUDA`
  (`:20`); diesen Testzielen fehlt das Define. Konfigurationsfehler in
  `test_cases/sqm_reference/CMakeLists.txt`, nicht im Quelltext.
- **Anzahl korrigiert (21.9.2026)**: hier stand „fuenf". Das war aus einem ABGEBROCHENEN
  `make -j8`-Protokoll gezaehlt — der Lauf haelt an, sobald eine Zielgruppe scheitert, und
  welche Ziele ueberhaupt versucht werden, schwankt. Ein vollstaendiger Lauf zeigt acht.
- **Warum es lange unbemerkt blieb**: `ctest -L gpu` meldet 200/200 und `ctest -R sqm` 335/335,
  weil diese Binaries dort nicht registriert sind bzw. nie gebaut werden. Ein gruener
  `ctest` beweist hier also **nicht**, dass `make` durchlaeuft — den Exit-Status separat pruefen
  (dieselbe Falle wie Known Issue #15).

### GFN-FF: Dispersionspaarliste hat denselben Architekturfehler wie die Repulsion (2026-09)
- **Status**: 🤖 BEHOBEN (uncommitted, zur Pruefung), siehe docs/GFNFF_PAIR_LIST_REFRESH.md. Offen: GPU-C6-Refresh nur in Gradientenaufrufen; Eintritt eines Paares von >60 Bohr nicht durch einen Test erzeugt.
- **Beim Review gefunden, vom Fix NICHT erfasst**: die ATM-Dreikoerperterme des D4-Terms
  (`t.C6_ij`/`t.C6_ik`/`t.C6_jk` auf jedem `ATMTriple`) teilen denselben Architekturfehler wie
  der Zweikoerperterm — einmalig bei `generateDispersionPairsNative()` berechnet, nirgends im
  Fix aktualisiert. Die Tripel-**Zugehoerigkeit** ist bindungsbasiert und korrekt unveraenderlich;
  die **C6-Werte** darauf nicht. Gemessen an einem kleinen Molekuel (Koffein, GFN-FF `-sp`):
  ATM-Term -7.5e-9 Eh gegen -1.8e-2 Eh Zweikoerperterm (~4e-7) — vermutlich meist vernachlaessigbar,
  aber nicht allgemein geprueft (grosse, dicht gepackte Systeme koennten anders liegen). Nicht
  weiterverfolgt angesichts der Groessenordnung und der Kosten einer erneuten vollen
  MOR41/GMTKN55/ctest-Validierung.
- Nach `ba0319dd` (Repulsions-Paarliste periodisch neu aufgebaut): die Dispersionspaarliste teilt
  denselben Aufbau — einmalig bei `InitialiseMolecule()`, nie neu aufgebaut. Bestaetigt durch
  den bestehenden `updateHBXBIfNeeded()`-Verbose-Log, der explizit `"Dispersion pairs" ... "(static)"`
  ausgibt.
- **Task**: analog zu `GFNFF::updateNonbondedRepulsionIfNeeded()` einen Rebuild-und-Re-Upload-Pfad
  fuer die Dispersion bauen, eigene Kostenmessung (Dispersion hat andere Paarzahlen/Cutoffs als
  Repulsion, die 3 ms/Schritt von dort uebertragen sich nicht automatisch).
- **Regressionsmassstab bereits etabliert**: `scripts/refset_regression.py` gegen MOR41/GMTKN55,
  muss wie bei der Repulsion bit-identisch bleiben.

### GFN-FF/xTB: `shouldUpdateHBXB()`s RMSD-Formel ist bei grossen Systemen praktisch wirkungslos (2026-09)
- **Status**: ⏳ OFFEN, quantifiziert, nicht behoben
- **Messung**: `rmsd = sqrt(sum_sq_diff) / natoms` (`gfnff_method.cpp:2630-2658`) statt der
  korrekten Pro-Atom-RMSD `sqrt(sum_sq_diff / natoms)` — es fehlt ein Faktor `sqrt(natoms)`.
  Bei `natoms=7320` braucht es **~13.6 Angstroem** mittlere Verschiebung pro Atom, bis der
  Standard-Schwellwert (`hb_update_rmsd_bohr` 0.3 Bohr, „pro Atom" gemeint) ueberhaupt feuert —
  praktisch nie in einer realen MD.
- **Herkunft**: treue Portierung der Fortran-Referenz (`gfnff_ini2.f90:717`), der Fehler steckt
  vermutlich auch dort — nicht curcuma-eigen, aber unbehoben.
- **Folge**: HB/XB-Paare werden bei grossen Systemen ebenso selten neu klassifiziert wie vorher
  die Repulsion es war — nur ohne die katastrophale Konsequenz, weil HB/XB energetisch schwaecher
  ist. Nicht als Ursache des polymer_2x-Absturzes bestaetigt, aber derselbe Fehlerklasse.
- **Task**: korrekte Formel (`sqrt(sum_sq_diff / natoms)`), gegen MOR41/GMTKN55 pruefen — HB/XB-
  Paarzahlen duerfen sich fuer kleine Systeme nicht aendern (dort ist der Faktor `sqrt(N)` klein
  genug, dass der Unterschied meist unter der Schwelle bleibt, aber nicht garantiert unter allen
  MOR41/GMTKN55-Strukturen).

### GFN-FF/EEQ: CPU/GPU-Trajektorien-Divergenz — Ursache direkt gezeigt, kein Fix (2026-09)
- **Status**: ⏳ OFFEN, Mechanismus bestaetigt, kein Loesungsweg umgesetzt
- **Befund**: ein 10-ps-GFN-FF-MD-Vergleich CPU vs. GPU auf `polymer_2x` (7320 Atome, nfrag=1500)
  divergiert reproduzierbar (nahezu deckungsgleich bis ~400 fs, ab ~1200 fs vollstaendig
  entkoppelt). Ursache **direkt gezeigt, nicht nur vermutet**: bei `nfrag=1500` waehlt der
  CPU-Pfad automatisch den iterativen PCG-Loeser fuer die EEQ-Gleichung
  (`eeq_ppcg_min_atoms=500`, `eeq_ppcg_min_nfrag=1`), der GPU-Pfad loest stets dicht/exakt
  (Schur-Cholesky) — zwei verschiedene, je fuer sich konvergierte numerische Verfahren fuer
  dieselbe Gleichung. Kontrollexperiment auf einem 8-Wasser-Cluster (24 Atome, nfrag=8),
  **identisches CPU-Binary**, nur `-eeq_solver.solve_method cholesky` gegen `ppcg`
  unterschiedlich: |dEpot| > 1e-5 Eh bei 0.18 ps, > 1e-3 Eh bei 0.94 ps, > 1e-2 Eh bei 1.52 ps —
  dasselbe Muster wie beim 300-fach groesseren CPU/GPU-Fall. Details:
  `Labor/curcuma MD-Stabilität großer Systeme.md` (Vault, Eintrag 23.9.2026).
- **Blockierter Loesungsversuch — SIGSEGV auf CUDA**: `-eeq_rocm_cpu_fragment_threshold N`
  (trotz Namens plattformuebergreifend, auf ROCm bereits validiert) sollte den GPU-Pfad auf den
  CPU-Loeser zwingen und damit gleiche Physik erzwingen. Auf CUDA stuerzt das reproduzierbar
  ab (unabhaengig von `-threads`): `Program received signal SIGSEGV` in
  `Eigen::internal::call_dense_assignment_loop` <- `EEQSolver::calculateFinalCharges` <-
  `GFNFF::prepareCNAndEEQ` <- `GFNFFGpuMethodImpl<GFNFFCudaBackend>::calculateEnergy` — reiner
  Eigen-Host-Code, kein CUDA-Kernel.
- **Root Cause**: `gfnff.h:1128` und `gfnff_method.cpp:1335` dokumentieren explizit „Safe for
  GPU path where CUDA corrupts heap metadata" — Eigen-Heap-Allokationen sind in einem
  CUDA-aktiven Prozess nicht sicher, deshalb existieren an mehreren Stellen memcpy-in-
  vorallozierte-Puffer-Umwege. `EEQSolver::calculateFinalCharges` ist eine grosse, generische
  Funktion (auch vom reinen CPU-Build genutzt), die frei alloziert — sie wurde nie dafuer
  gehaertet, weil sie vorher nie aus einem CUDA-Prozess heraus aufgerufen wurde (die
  `eeq_rocm_cpu_fragment_threshold`-Route war bisher nur auf ROCm/HIP validiert, das dieses
  Heap-Problem offenbar nicht hat).
- **Task**: (a) `EEQSolver::calculateFinalCharges` alloc-frei machen fuer den CUDA-Fall (invasiv,
  geteilter Code, betrifft auch den reinen CPU-Pfad) — dann erneut auf CUDA testen; (b)
  alternativ ein bezahlbarer dritter Loeser, der auf CPU und GPU identisch reproduzierbar ist,
  auch bei nfrag~1500; (c) alternativ akzeptieren und Trajektorien ab ~ps-Zeitskala nur noch
  ensemble-/statistisch statt punktweise vergleichen. Keine Entscheidung getroffen.
- **Nicht versucht, weil unbezahlbar**: CPU auf den GPU-exakten (dichten) Loeser zwingen
  (`-eeq_solver.eeq_ppcg_min_nfrag` sehr hoch) — bei nfrag=1500 kostet die exakte Loesung laut
  eigener Parameterbeschreibung „+nfrag zusaetzliche Faktorisierungen"; >30 Minuten reichten
  nicht fuer einen einzigen MD-Schritt.
- **SIGSEGV bestaetigt als Race Condition, kein billiger Fix (2026-09, selbe Sitzung)**:
  `MALLOC_ARENA_MAX=1` (Standard-Gegenmassnahme fuer CUDA/glibc-Heap-Interaktionen) getestet,
  nicht nur vermutet — 3 Wiederholungen je Konfiguration: 0/3 erfolgreich MIT der Variable, 1/3
  erfolgreich OHNE. Kein Unterschied, beide nichtdeterministisch. Bestaetigt exakt die
  Charakterisierung in `docs/TECHNICAL_DEBT.md` F-Q9/D-26/D-46 („Root cause uninvestigated") —
  kein neuer Fund, aber eine konkrete Absturzrate fuer diesen Aufruf. Task (a) oben bleibt der
  einzige echte Loesungsweg fuer diesen Zweig; kein Kurzschluss ueber Umgebungsvariablen.
- **Toleranz-Straffung getestet und verworfen (2026-09, selbe Sitzung)**: die naheliegende Idee,
  eine straffere `eeq_ppcg_tol` wuerde die Divergenz verzoegern, ist widerlegt. `eeq_ppcg_tol`
  1e-6/1e-9/1e-12/1e-14 auf dem 8-Wasser-Kontrollsystem liefern **bit-identische** Trajektorien —
  PCG konvergiert dort schon bei der lockersten Anforderung in einer Iteration auf
  Maschinengenauigkeit (|Pr|=2.66e-16, per `-verbosity 3` bestaetigt). PCG ist hier nicht „zu
  ungenau"; die Abweichung zu Cholesky sitzt in einem winzigen algebraischen/Rundungsunterschied
  zwischen den Loesungswegen selbst, nicht im PCG-Konvergenzgrad — ein Toleranz-Regler kann diese
  Klasse Divergenz also nicht beheben, wenn beide Seiten bereits nahe Maschinengenauigkeit loesen.
- **GPU-eigener PCG-Pfad (WP7-C) getestet und verworfen — nicht an Korrektheit, an Performance
  (2026-09, selbe Sitzung)**: `-gfnff.solve_method pcg` engagiert sauber den geraeteresidenten
  PCG-Loeser (`cuda/eeq_solver_gpu.cu`, WP7-C Block-Jacobi, umgeht den SIGSEGV komplett — kein
  CUDA-Prozess ruft CPU-EEQ-Code auf). Smoke-Test (3 MD-Schritte, `polymer_2x`, nfrag=1502)
  lief fehlerfrei, aber **~9 Minuten fuer 3 Schritte**; die geplante 2-ps-Zielmessung (2000
  Schritte) wurde nach 39,5 Minuten ohne einen einzigen fertigen Print-Zyklus abgebrochen
  (~3 Min/Schritt hochgerechnet, ~100 h fuer 2 ps). Prozess dabei aktiv (113% CPU, GPU-Auslastung
  96%), nicht haengend. Gegenprobe auf dem 8-Wasser-Cluster (nfrag=8): 100 Schritte in 1.45 s —
  **dasselbe nfrag-Skalierungsmuster wie der exakte CPU-Loeser** (Known-Issue-Eintrag oben,
  „+nfrag zusaetzliche Faktorisierungen"). WP7-C ist also fuer nfrag~1500 technisch korrekt,
  aber praktisch unbezahlbar — kein Kurzschluss fuer dieses System.
- **WP7-C-Ursache gefunden (2026-09-24, neues Teilprojekt „EEQ-Loeser-Benchmark")**: nicht die
  Fragmentzahl ist das Problem, sondern die Groesse des **groessten** Fragments.
  `buildBlockJacobiFactors` (`cuda/eeq_solver_gpu.cu:1441-1516`) baut je Fragment eine explizite
  dichte Inverse (`cusolverDnDpotrf`+`cusolverDnDpotri`, je 2 `cudaStreamSynchronize` seriell in
  einer Host-Schleife). Fuer `polymer_2x`s 1410-Atom-Polymerfragment ist das eine dichte
  1410x1410-Matrixinversion — vergleichbar teuer wie der exakte Loeser selbst. Kontrolliert
  nachgewiesen: bei **konstanter** Fragmentzahl (500) skaliert die Zeit allein mit der Groesse
  des groessten Fragments (3 Atome: 5.8s / 500 Atome: 11.4s / 1410 Atome: 52.4s, fuer 2
  MD-Schritte) — reine Fragmentzahl-Skalierung (uniform klein) ist dagegen mild (250->1000
  Fragmente: 2.4->20.9s). Die vorhandene Sicherung `GPU_BLOCK_JACOBI_MAX_NF = 2048`
  (`eeq_solver_gpu.cu:1434`) greift nicht (1410 < 2048). **Naheliegender Fix, nicht umgesetzt**:
  Schwellenwert deutlich senken oder bei grossen Fragmenten auf den guenstigeren
  Mittelwert-Projektions-Vorkonditionierer zurueckfallen, den die CPU-Variante
  `solveWithProjectedPCG` bereits benutzt. Details, Messwerte, Reproduktionskommandos:
  [[curcuma EEQ-Löser-Benchmark]] (Vault, Projekt + Labor).
- **Vollstaendige Matrix (2026-09-24, selbe Sitzung, 6 Gruppen, ~50 Systeme) verfeinert den
  obigen Befund: ZWEI unabhaengige Kostentreiber, nicht einer.** Ein reiner
  Fragmentgroessen-Schwellenwert (Fix oben) reicht NICHT fuer `polymer_2x`-artige Systeme.
  - **G1 (reines Wasser, kein Ausreisser)**: GPU-PCG waechst schon bei uniform KLEINEN
    Fragmenten ueberlinear (~nfrag^1.8) und wird bei nfrag~2440 unbezahlbar (>150s) — rein durch
    die FragmentANZAHL, kein grosses Fragment beteiligt.
  - **G3 (feste Polymer-Atomzahl ~2820 + 1500 Wasser, Aufspaltung 1/2/4/8/16 Ketten, nfrag~1500)**:
    GPU-PCG timeoutet (>120s) bei JEDER Aufspaltung — die hohe Grund-Fragmentzahl allein
    dominiert schon so stark, dass die Fragmentgroessen-Variation keinen Unterschied mehr macht.
    CPU-ppcg und GPU-chol werden dagegen mit mehr/kleineren Fragmenten deutlich schneller
    (>120s bei 1 Riesenfragment -> 7.6-10.0s bei 16 kleinen Ketten) — ein einzelnes sehr grosses
    Fragment ist auch fuer die "normalen" Loeser schlecht, nicht nur fuer WP7-C.
  - **G4 (fixe 2x1410-Ketten, Wasser 0->4500)**: scharfer Schwelleneffekt — bei water=0 (nfrag=2)
    ist GPU-PCG mit 24.75s noch im normalen Bereich; sobald 500 Wassermolekuele dazukommen
    (nfrag=502), springt es sofort auf >120s.
  - **G5 (gleiche Atomzahl ~7320, extreme Fragmentzahl-Kontraste, DIREKTER Vergleich wie
    angefordert)**: 1 Riesenfragment (nfrag=1) -> CPU-ppcg UND GPU-chol timeouten EBENFALLS
    (>150s, s. separater Befund unten); polymer_2x-artig (nfrag=1502) -> CPU-ppcg 33.4s,
    GPU-chol 33.2s, GPU-PCG >150s; 2440 reines Wasser -> CPU-ppcg 6.3s (schnellster Fall!),
    GPU-chol 16.4s, GPU-PCG >150s. Bei GLEICHER Atomzahl variiert CPU-ppcg allein durch
    Fragmentierung um Faktor >24x.
  - **G6 (`polymer_2x` verdoppelt, 4x1410+3000 Wasser, 14640 Atome)**: CPU-ppcg 77-85s, GPU-chol
    98-106s (beide schon spuerbar langsam), GPU-PCG >180s.
  - **Konsequenz fuer einen Fix**: (a) Fragmentgroessen-Schwellenwert senken UND (b) den
    seriellen Zwei-Sync-pro-Fragment-Aufbau (`buildBlockJacobiFactors`) umbauen (z. B. batched
    cuSOLVER statt serieller Host-Schleife, ungeprueft ob cuSOLVER das unterstuetzt) — nur (a)
    reicht nicht, weil die Fragmentanzahl-Komponente bei `polymer_2x`-typischem nfrag~1500
    unabhaengig von der Fragmentgroesse schon dominiert.
  - **Eigenstaendiger Nebenbefund, NICHT Teil der Solver-Frage**: ein einzelnes
    7320-Atom-Riesenfragment (G5) laesst auch die normalerweise schnellen Pfade (CPU-ppcg,
    GPU-chol) timeouten. Log-Pruefung: haengt bereits bei "Move structure to the origin...",
    also VOR jeder EEQ-Solver-Aktivitaet. Eine O(N^3)-Floyd-Warshall-Ursache in
    `EEQSolver::computeTopologicalDistancesSparse` (`eeq_solver.cpp:3140`) wurde geprueft und
    beim Codelesen verworfen — die dortige Dijkstra-Variante hat einen fruehen Abbruch bei 12 A
    topologischem Abstand, sollte fuer eine Kettentopologie schnell sein. **Ursache offen**,
    eigene Untersuchung noetig, nicht mit WP7-C vermischen.
  - **Betriebslehre**: ein per Inline-Timeout automatisch in den Hintergrund verschobener
    Prozess (kein explizites `run_in_background:true`) ueberlebte einen Turn-/Kontextwechsel
    NICHT zuverlaessig (Sweep brach unbemerkt mitten in einer Gruppe ab). Mit explizitem
    `run_in_background:true` lief er durch. Fuer kuenftige lange Hintergrundlaeufe: immer
    explizit `run_in_background:true` verwenden.
  - **Verweis**: [[curcuma EEQ-Löser-Benchmark]] (Vault, Projekt + Labor, vollstaendige Tabellen).
- **Fix (a) umgesetzt, gebaut, real getestet — hilft, loest aber NICHT `polymer_2x` (2026-09-24,
  selbe Sitzung, dirty, nicht committet)**: zwei neue PARAMs (`eeq_solver.h`)
  `gpu_block_jacobi_max_frag_atoms` (Default 300) und `gpu_block_jacobi_max_nfrag` (Default 400),
  durchgereicht `buildBlockJacobiFactors` -> `solveWithDeviceRHSAndGPUPCG` ->
  `GFNFFGpuMethodImpl` (liest per `gfnff_cfg.value(...)`, gleiches Muster wie
  `max_pcg_iterations`). ROCm-Stub-Signatur mitgezogen (`EEQSolverHip` ist Alias auf
  `EEQSolverGPU`), unkompiliert (kein SDK hier).
  - **Validiert (klein, funktioniert wie geplant)**: 8-Wasser-Cluster (nfrag=8), 40 Schritte in
    1.16s, Block-Jacobi bleibt korrekt aktiv (`EEQ GPU PCG: block-Jacobi preconditioner active`
    in jedem Schritt) — kein Kollateralschaden.
  - **KORREKTUR der eigenen Einschaetzung, durch echten Test auf `polymer_2x` gefunden**: bei
    nfrag=1502 (reales UND synthetisches `polymer_2x`-Aequivalent getestet) bleibt WP7-C
    weiterhin >120s ohne Abschluss — **keine messbare Verbesserung**. Root Cause im Code
    gefunden: `solveWithDeviceRHSAndGPUPCG` (`eeq_solver_gpu.cu:1772-1785`) laeuft UNABHAENGIG
    vom Block-Jacobi-Status eine Schleife `for (f=0; f<nfrag; ++f) { runSinglePCG(...); }` —
    **(nfrag+1) vollstaendig separate PCG-Solves**, bei nfrag=1502 also 1503 einzelne Laeufe mit
    je bis zu 200 Iterationen und je einem O(N^2)-Matvec. Das ist die GPU-Entsprechung der
    "+nfrag zusaetzliche Solves", die der CPU-seitige projizierte Ansatz
    (`solveWithProjectedPCG`, EIN Solve statt nfrag+1) gezielt vermeidet — laut eigener
    PARAM-Beschreibung der CPU-Variante, nur nie auf WP7-C uebertragen erkannt. Der gefixte
    Block-Jacobi-Aufbau war ein ZUSAETZLICHER, nicht der einzige Kostentreiber.
  - **Der Fix bleibt sinnvoll**, nur enger als gedacht: echte, nie-schlechtere Verbesserung fuer
    grosse Einzelfragmente bei moderater Fragmentzahl (bis ~400 Fragmente), NICHT fuer
    `polymer_2x`-Groessenordnungen. Ein Fix dafuer braucht einen GPU-Port des projizierten
    Ansatzes (ein Solve statt nfrag+1) — eine neue Implementierung, kein Parameter-Tuning mehr,
    deutlich groesserer Umfang als heute umgesetzt.
  - **Methodenlehre**: ein Fix ist erst geprueft, wenn er auf dem Fall getestet wurde, der ihn
    ausgeloest hat (hier: `polymer_2x`) — nicht nur auf dem Fall, den die eigene Diagnose als
    Ursache identifiziert hatte (hier: Block-Jacobi-Aufbau bei moderatem nfrag). "Fix compiliert
    und verhaelt sich wie geplant" haette fast zu einer stillschweigend falschen
    Erfolgsmeldung gefuehrt.
  - **Verweis**: [[curcuma EEQ-Löser-Benchmark]] (Vault, Projekt + Labor, volles Testprotokoll).
- **Methodenfalle bei diesem Benchmark**: `-sp` (Einzelpunkt) verfehlt den Block-Jacobi-Aufbau
  komplett — der laeuft laut Code nur bei einem "PCG refactor", nicht beim ersten Solve. Ein
  `-sp`-Sweep bis nfrag=1000 zeigte faelschlich "kein Problem" (alle ~15-20ms). Erst ein
  2-Schritte-MD-Lauf deckte den Effekt auf. Fuer jeden kuenftigen Solver-Benchmark: mindestens
  2 MD-Schritte, nie nur `-sp`.
- **Zwischenstand nach drei geprueften Hebeln**: fuer `polymer_2x` (nfrag=1502) sind jetzt alle
  drei naheliegenden „gleicher Algorithmus auf beiden Seiten"-Wege empirisch gescheitert — zwei
  an Korrektheit/Stabilitaet (SIGSEGV, Toleranz wirkungslos), einer an Performance (GPU-PCG).
  Verbleibend realistisch: (a) SIGSEGV am Ursprung fixen — einzige Option ohne zusaetzliche
  Performance-Baustelle, da nur die ohnehin-PCG-nutzende CPU betroffen waere; (b) WP7-C fuer
  grosse Fragmentzahlen performant machen (eigenes, ungeschaetztes Arbeitspaket, z. B. warum
  1502 Fragmente ~2000x langsamer sind als 8 — nicht root-caused); (c) akzeptieren und nur noch
  ensemble-/statistisch statt punktweise vergleichen. Keine Entscheidung getroffen.
- **Verweis**: [[Governance-Regeln für KI-Coding-Agenten an wissenschaftlichem Code]] (Vault) —
  generalisierbare Lehre „zwei numerisch verschiedene, je fuer sich korrekte Loeser lassen
  dieselbe Physik chaotisch auseinanderlaufen" (inkl. der Toleranz-Korrektur im selben
  Abschnitt).
- **GPU-Port des projizierten Ansatzes (WP7-E) umgesetzt, gebaut, real auf `polymer_2x` getestet
  — loest den nfrag+1-Solves-Engpass (2026-09-24, AI-implemented, machine-tested)**: Option (b)
  aus obigem Zwischenstand umgesetzt, kein Parameter-Tuning — ein neuer GPU-Loeser
  `EEQSolverGPU::solveWithDeviceRHSAndGPUProjectedPCG` (`cuda/eeq_solver_gpu.cu`/`.h`), der
  direkt `EEQSolver::solveWithProjectedPCG` (CPU, `eeq_solver.cpp`) portiert: EIN PCG-Solve im
  Tangentialraum der Constraints (Projektion `v -= Mittelwert(v)` pro Fragment via 5 neuen
  Kernen `k_frag_sum/scatter_add/project_delta/feasibility_delta/precond_correct`) statt WP7-Cs
  `nfrag+1` unabhaengiger Solves + Block-Jacobi-Aufbau. Nur diagonaler Jacobi-Praekonditionierer
  (wie die CPU-Referenz) — der teure Block-Jacobi-Aufbau (WP7-D) entfaellt komplett, nicht nur
  fuer grosse Fragmente. Eingehaengt ueber die bestehenden CPU-PARAMs `eeq_ppcg_min_nfrag`/
  `eeq_ppcg_min_atoms` (Auto-Modus bevorzugt WP7-E ab denselben Schwellen wie die CPU — WP7-C
  bleibt nur bei explizitem `-gfnff.solve_method pcg` erreichbar); ROCm-Stub (`return false`,
  wie die anderen unportierten WP7-Pfade).
  - **Korrektheit**: `-sp` auf `polymer_2x_gfnff_opt.xyz` (7320 Atome), `solve_method cholesky`
    (WP7-A, exakt) vs. `solve_method ppcg` (WP7-E): **-917.33663023 Eh beide, auf 8 Nachkomma-
    stellen identisch**. Wasser8-Cluster (nfrag=8, N=24, unter der Auto-Schwelle) bleibt
    unveraendert auf WP7-A, keine Regression.
  - **Performance, auf dem tatsaechlich ausloesenden Fall getestet** (nicht nur dem, auf den die
    eigene Diagnose zeigte — Lehre aus dem vorherigen Fix-Versuch beherzigt): 8-Schritte-GFN-FF-
    MD auf `polymer_2x_gfnff_opt.xyz` (N=7320, nfrag=1502, `-gpu cuda -gfnff.solve_method auto`,
    dt=0.5 fs) **fertig nach 23 s**, jeder Schritt zeigt `path=WP7-E GPU-Schur (projected-pcg)`.
    WP7-C blieb auf demselben System zuvor bei 120 s UND 300 s Timeout ohne jeden Fortschritt
    haengen (s. Eintrag oben). Energien ueber die 8 Schritte physikalisch plausibel
    (-917.3 -> -914.1 Eh, glatte Gleichgewichtseinstellung), Ladungen unauffaellig
    (|q|max ~0.7-0.72, Mittelwert 0).
  - **Damit vorlaeufig entschieden**: Ensemblegleichheit (Option c) war nicht noetig — echte
    Loeser-Reproduzierbarkeit (WP7-A cholesky == WP7-E ppcg auf 8 Nachkommastellen) war
    erreichbar, wie in der urspruenglichen Anweisung gefordert ("Wir akzeptieren
    Ensemblegleichheit erst, wenn wir keine andere Chance haben").
  - **Noch offen / nicht getestet**: laengere MD-Laeufe (Energieerhaltung ueber >8 Schritte),
    ROCm-Portierung, Vergleich CPU-ppcg-Timing vs. GPU-WP7-E-Timing auf derselben
    Groessenordnung (nur WP7-C-vs-WP7-E GPU-intern gemessen).
  - **Verweis**: [[curcuma EEQ-Löser-Benchmark]] (Vault, Projekt + Labor, wird nachgezogen).
- **G6-Doppelsystem (14640 Atome, nfrag=3004) mit WP7-E getestet — konvergiert korrekt, aber die
  Gesamtlaufzeit ist fast vollstaendig CPU-Setup, nicht der EEQ-Solve (2026-09-24)**:
  `g6_polymer4x.xyz` (4x1410-Atom-Ketten + 3000 Wasser), 8-Schritte-GFN-FF-MD,
  `-gfnff.solve_method auto`. `path=WP7-E GPU-Schur (projected-pcg)` in allen 8 Schritten,
  Energien/Ladungen plausibel, Gesamtlaufzeit **118 s**. Die eingebaute Timing-Aufschluesselung
  zeigt: der eigentliche EEQ-Solve kostet nur **~260 ms/Schritt** (~2 s ueber alle 8 Schritte) —
  WP7-E ist bei dieser Groesse kein Engpass. Die restlichen ~116 s sind **einmaliges CPU-Setup**
  (`One-time setup topo=49076.7 ms param=103560.0 ms`, identisch in jedem Schritt-Report
  ausgegeben, also wirklich einmalig, nicht kumulativ pro Schritt): das synthetische
  3000-Wasser-Sub-System erzeugt **1.028.354 Wasserstoffbruecken-Kandidaten** in der
  Topologie-Erkennung — derselbe, in `docs/GFNFF_PERFORMANCE_LEVERS.md` ("Lever 1")
  dokumentierte, vom EEQ-Solver unabhaengige Engpass. **Folge**: der naive Vergleich mit den
  alten Matrix-Zahlen (CPU-ppcg 77-85 s, GPU-chol 98-106 s, GPU-PCG/WP7-C >180 s Timeout) ist
  irrefuehrend, weil diese vermutlich denselben solverunabhaengigen Setup-Anteil enthalten —
  ohne deren eigene Aufschluesselung laesst sich der reine Solver-Unterschied nicht sauber
  herausrechnen. Sauber belegt ist nur: WP7-E selbst bleibt auch bei doppelter Systemgroesse
  (gegenueber `polymer_2x`) trivial schnell; der offene Punkt fuer diese Systemklasse (viele
  kleine Fragmente/grosse Wasserboxen) ist die CPU-Setup-Zeit, ein separates, bereits bekanntes
  Problem.

- **HB-/XB-Fix (Lever 1) umgesetzt, gebaut, korrekt — loest aber NICHT den G6-Setup-Engpass;
  eine Kopplung an den Bond-Term erst gefunden und dann korrekt eingegrenzt (2026-09-24,
  AI-implemented, machine-tested)**: Auftrag „HB-Fix und XB-Fix umsetzen" (energiebasiertes
  Pruning statt distanzbasiert, s. Briefing oben). Umgesetzt: zwei neue GFNFF-Methoden
  `estimateHBStrengthCase1`/`estimateXBStrength` (`gfnff_method.cpp`), die die EXAKTE
  Energie-Formel aus dem Energie-Kernel (`ff_workspace_gfnff.cpp calcHydrogenBonds`/
  `calcHalogenBonds`) vorzeitig — vor der Struct-Allokation — auswerten und Kandidaten unterhalb
  eines neuen Schwellenwerts (`hb_min_pair_energy_eh`/`xb_min_pair_energy_eh`, Default 1e-9 Eh)
  verwerfen. Die vier gemeinsam genutzten Damping-Grundfunktionen (`ws_damping_*`,
  `ws_charge_scaling`) wurden aus `ff_workspace_gfnff.cpp`s anonymem Namespace in den geteilten
  Header `gfnff_par.h` verschoben, damit Detektion und Energie-Kernel garantiert dieselbe
  Implementierung nutzen (kein Drift-Risiko).
  - **Echter Bug gefunden und korrigiert, BEVOR er unbemerkt geblieben waere**: ein erster
    Versuch wendete das Pruning auf ALLE HB-Faelle an (case 1 „unbound" UND case 2/3/4
    „donor-gebunden"). Der HB-Term selbst blieb dabei korrekt (Energieaenderung <1e-8 Eh auf
    `triose.xyz`), aber die GESAMTENERGIE verschob sich um **0,18 kcal/mol** (2,87e-4 Eh) — fast
    vollstaendig im BOND-Term. Ursache: jedes case-2/3/4-Kandidat mit N/O-Akzeptor speist
    `bond_hb_data`/`hb_cn_H` (`ff_workspace_gfnff.cpp computeHBCoordinationNumbers`), eine rein
    GEOMETRISCHE erf-basierte Zaehlgroesse, die den donor-H-BOND-Term reskaliert (`egbond_hb`,
    `VBOND_SCALE`) — und diese Groesse korreliert NICHT mit |E_HB| (ein Akzeptor mit schwacher
    Basizitaet kann geometrisch nah sein und voll zur Zaehlung beitragen, aber wenig zur Energie).
    Gefunden nur, weil die VOLLE Energie-Dekomposition verglichen wurde, nicht nur der HB-Term
    selbst — eine reine "HB-Term unveraendert"-Pruefung haette den Bug durchgelassen.
  - **Fix**: Pruning NUR fuer case 1 (unbound A...H...B) und XB — beide strukturell sicher, da
    ein case-1-H per Konstruktion an keinen der beiden flankierenden Atome gebunden ist und daher
    nie den `bond_hb_data`-Lookup-Key treffen kann (kein XB-Analogon zu `hb_cn_H` gefunden).
    case 2/3/4 werden NIE mehr geprueft, unabhaengig vom Schwellenwert. Die dadurch unsicher
    gewordene Methode `estimateHBStrengthCase2or4` wurde vollstaendig entfernt (Definition +
    Deklaration), nicht nur deaktiviert.
  - **Korrektheit nach dem Fix, real getestet**: `triose.xyz`, volle Energie-Dekomposition,
    `hb_min_pair_energy_eh 0` (alt) vs. Default (neu) — **bitidentisch** in JEDER Komponente
    (Bond -9.6737006656 Eh beide, H-bonds -0.0122278523 Eh beide, case-1/2-Zaehlungen 1495/107
    beide identisch — bei diesem kleinen Molekuel greift die case-1-Schwelle gar nicht). `ctest -L
    gfnff`: **77/78 bestanden**, einziger Fehlschlag der bereits bekannte vorbestehende
    `cli_curcumaopt_07_opt_multixyz` (golden-value drift, nicht mit diesem Fix zusammenhaengend).
  - **Wirkung auf `g6_polymer4x.xyz`**: HB-Kandidaten **1.028.354 -> 205.284 (5x weniger)**,
    Energien der 8-Schritte-MD gegenueber dem reinen-WP7-E-Lauf um **~1,1e-5 Eh** verschoben
    (relativ ~5e-9, konsistent mit dem Verwerfen echt vernachlaessigbarer case-1-Kandidaten).
    **Aber**: Gesamtlaufzeit **124 s vs. 118 s Baseline — keine messbare Verbesserung**, weil
    (siehe Eintrag oben) die HB-Detektion selbst nur ~1,4 s der ~150 s Setup-Zeit ausmacht; der
    eigentliche G6-Engpass liegt woanders (Topologie-Erkennung + sonstige Parametergenerierung).
    Der Fix reduziert aber die PRO-SCHRITT-Kosten des HB-Energie-Kernels (5x weniger Eintraege
    in `calcHydrogenBonds` pro Aufruf) — bei laengeren MD-Laeufen (viele Schritte statt 8) sollte
    das kumulativ sichtbar werden, wurde in dieser Sitzung aber nicht gemessen.
  - **Methodenlehre**: dieselbe Regel wie beim WP7-E-Fix — ein Fix ist erst geprueft, wenn er auf
    der VOLLEN Rechnung getestet wurde, nicht nur auf dem Term, den er direkt aendert. Ein „der
    HB-Term aendert sich kaum" waere hier eine stillschweigend falsche Erfolgsmeldung gewesen.
  - **Noch offen**: der eigentliche G6-Setup-Engpass (Topologie-Erkennung, ~49 s, und
    sonstige Parametergenerierung, ~100 s abzueglich HB) ist nicht root-caused. Lever 2 aus dem
    Briefing (die verbleibenden O(N²)-Schleifen cell-listen — `nb_hc`/`nb_nometal`,
    BATM-Scan, Bond-BFS, = Phase 5 des Multi-GPU-Plans) ist der naheliegende naechste Schritt.

- **G6-Root-Cause AUFGEKLAERT: der dominante Teil (48s/98.5s Topologie, ueber beide q-Loop-
  Paesse) ist ein Testdaten-Artefakt von `gen_system.py`, kein GFN-FF-Bug (2026-09-25,
  verifiziert)**: `CURCUMA_GFNFF_PROFILE=1` (summiert ueber BEIDE q-Loop-Paesse — der normale
  verbosity-2-Report ueberschreibt Pass 1 mit Pass 2 und zeigt nur die Haelfte) deckte auf:
  „Hueckel pi bond orders" allein kostet **48.275 s von 98.540 s Gesamttopologie (49%)** —
  mehr als die urspruenglich vermuteten O(N²)-Nachbarlisten (nb_hc+nb_nometal, 8.8s, 9%)
  zusammen. Root-Cause-Recherche (per `CURCUMA_HUCKELDUMP=1`): die "4 π-Systeme", die der
  Report zeigt, sind KEINE kleinen Fragmente — es sind die 4 KOMPLETTEN synthetischen
  Polymerketten (~2820 Atome je System, 77% des ganzen Molekuels!), weil `GFNFF::
  detectPiSystems` (`gfnff_method.cpp:6866-6982`) **2816 von 2819 gesaettigten
  Alkyl-Kohlenstoffen faelschlich als sp2 (hyb=2) statt sp3 klassifiziert** fand. Die 48s
  kommen exakt aus `HuckelSolver::solveAndBuildDensity` (`huckel_solver.cpp:445`): eine
  DICHTE O(N³)-Eigenzerlegung (`Eigen::SelfAdjointEigenSolver`) mit ndim≈2820, viermal, mal
  zwei q-Loop-Paesse = 8 dichte ~2820×2820-Diagonalisierungen — rechnerisch exakt konsistent
  mit den gemessenen 48s. **Entscheidende Verifikation auf dem ECHTEN Molekuel**
  (`polymer_2x_gfnff_opt.xyz`, 7320 Atome [korrigiert 2026-09-25, stand hier faelschlich 1410 = ein Polymerstrang], reale Geometrie statt synthetisch generiert):
  derselbe `CURCUMA_HUCKELDUMP`/`CURCUMA_GFNFF_PROFILE`-Lauf findet **0 von 0 π-Systemen**,
  Hueckel-Kosten **89.1 ms** (statt 48.275.000 ms — Faktor 540000x weniger). **Damit ist
  geklaert: das ist KEIN GFN-FF-Korrektheitsbug, der bei echten Molekuelen (lange
  Alkylketten, Lipide, reale Polymere) auftreten wuerde — es ist ein Artefakt der Art, wie
  `gen_system.py` (Scratchpad-Skript dieser Sitzung, nicht Teil des Repos) seine synthetischen
  Ketten baut (vermutlich unrealistische Bindungswinkel/-laengen, die die Hybridisierungs-
  Heuristik taeuschen).** Kein dringender Fix noetig; die 48s aus dem G6-Test sind nicht
  repraesentativ fuer echte Systeme und sollten nicht als Grundlage fuer weitere
  Performance-Entscheidungen dienen. **Konsequenz fuer Lever 2**: bleibt der legitime,
  allgemeingueltige naechste Schritt (nb_hc/nb_nometal cell-listen, 8.5s auf G6, aber
  distanzbasiert und damit unabhaengig vom Testdaten-Artefakt — gilt genauso auf echten
  Systemen). Der groessere, weniger cell-list-geeignete Posten dahinter (`bpair`/
  `topo_distances`, dichte N×N-Matrizen, ~1,7 GB bei N=14640) braucht stattdessen eine
  Sparse-Umstellung (Aufwand M-L, mittleres Risiko — mehrere Call-Sites mit bereits frueher
  gefixten Bugs, s. Known Issues #6(i)/#21(l)/#24(a)).

- **Lever 2 umgesetzt: `nb_hc`/`nb_nometal` per SpatialCellList, 17x schneller, korrekt
  (2026-09-25, AI-implemented, machine-tested)**: `pair_bonded_no_fm` (der distanzbasierte,
  elementpaar- und ladungsabhaengige Bindungstest hinter icase 2/3) bleibt UNVERAENDERT — nur
  die Kandidatengenerierung wird von O(N²) auf eine `SpatialCellList` (bereits fuer HB/XB/
  Repulsion im Einsatz) umgestellt, mit einem strengen, im Code selbst berechneten oberen
  Abstands-Bound (Scan ueber alle 86 Elementpaare bei grosszuegiger CN-Klammer, `ff`-Faktor
  NICHT als <=1 angenommen, da `rab_p` negative Eintraege hat, plus 3,5-Bohr-Marge fuer den
  ladungsabhaengigen `qshift`-Term, der bei negativer Ladung `rco` VERGROESSERT statt
  verkleinert). Kandidatenreihenfolge innerhalb einer Zeile ist beim Cell-List-Pfad nicht mehr
  aufsteigend nach Atomindex (anders als der alte O(N²)-Scan) — deshalb `std::sort` nach dem
  Aufbau, um Bit-Identitaet zu garantieren statt sie anzunehmen. Gated auf
  `nb_cell_list_min_atoms` (Default 800, dasselbe PARAM wie bei HB/XB), kleine Systeme bleiben
  auf dem alten O(N²)-Pfad.
  - **Korrektheit, real getestet**: `triose.xyz` (66 Atome, unter der Schwelle, alter Pfad)
    bitidentisch zum committeten Stand. `polymer_2x_gfnff_opt.xyz` (7320 Atome [korrigiert 2026-09-25, stand hier faelschlich 1410], ueber der
    Schwelle, NEUER Cell-List-Pfad aktiv) — **volle Energie-Dekomposition bitidentisch** zum
    Vor-Lever-2-Lauf in JEDER Komponente (Bond -820,6570687703 Eh, H-bonds -4,5521241259 Eh,
    case-1/2-Zaehlungen 10370/43703, alle identisch). `ctest -L gfnff`: **77/78 bestanden**,
    einziger Fehlschlag der bekannte vorbestehende `cli_curcumaopt_07_opt_multixyz`.
  - **Performance, auf G6 gemessen** (`CURCUMA_GFNFF_PROFILE=1`, ueber beide q-Loop-Paesse):
    `nb_hc list` **4309,2 ms → 253,6 ms**, `nb_nometal list` **4153,6 ms → 247,5 ms** — Faktor
    **~17x** auf beiden, ~8 s gespart, exakt wie aus der Lever-2-Recherche erwartet. Da das
    Kriterium rein distanzbasiert ist (nicht an das Hueckel-Testdaten-Artefakt gekoppelt),
    gilt der Gewinn genauso auf echten grossen Systemen.
  - **Noch offen** (zum Zeitpunkt dieses Eintrags): `bpair`/`topo_distances` — erledigt im
    naechsten Eintrag.

- **`bpair`/`topo_distances` sparse umgestellt — bit-identisch, ~1,75 GB weniger Peak-RSS
  bei 14640 Atomen (2026-09-25, AI-implemented, machine-tested)**: beide Tabellen (und die
  dichte `inL`-Hilfsmatrix in `computeBpairNbondmat`) sind jetzt eine `SparseTopoTable`
  (`gfnff.h`): pro Atom nur die Partner bis 3 (bpair) bzw. 5 (BFS) Bindungen, Zeilen nach j
  sortiert, alles andere liest den Fernwert (5 bzw. 999). Algorithmen unveraendert (nbondmat
  Level 1 + zwei pairsbond-Runden, tiefenbegrenzter BFS); die neuen pairsbond-Tags werden pro
  Zeile gesammelt und danach seriell in beide Zeilen gemischt. BATM-Scan laeuft nur noch ueber
  gespeicherte Paare `j<i` (gleiche Reihenfolge wie der dichte Scan), XB-Filter und H-H-Repulsion
  lesen per `get()`.
  - **Korrektheit**: temporaerer Eintrag-fuer-Eintrag-Vergleich dicht vs. sparse (vor dem
    Commit entfernt) — **0 Abweichungen** in beiden Tabellen und beiden q-Loop-Paessen auf
    triose, MOR41 PR26/PR27/PR28/PR34/ED33 (eta-Komplexe, XB) und `polymer_2x` (7320).
    `scripts/refset_regression.py` gegen ein Referenzbinary aus einem Worktree auf `526fcb8a`:
    MOR41 **95/95**, GMTKN55 gfnff **2462/2462** ohne jede Abweichung. `polymer_2x` und G6: alle
    14 Energiekomponenten auf 10 Nachkommastellen identisch, G6 22504 BATM-Tripel beide.
    `ctest -L gfnff` 77/78 (bekannter `cli_curcumaopt_07`; mit Alt- und Neubinary identische
    Frame-Energien und Geometrien). S30L-CI nicht geprueft (Strukturen lokal nicht vorhanden).
  - **Messung** (`-sp -threads 16`, `CURCUMA_GFNFF_PROFILE=1`, beide q-Loop-Paesse): Phase
    „topo distances + BATM list" G6 (14640) **1083,8 → 11,6 ms**, `polymer_2x` (7320)
    **309,3 → 7,0 ms**. Peak-RSS G6 **14,49 → 12,74 GB**, `polymer_2x` 3,63 → 3,44 GB. Wandzeit
    G6 60,5 → 59,9 s — die Phase war nur ~1 s davon; der Gewinn ist hauptsaechlich Speicher.
  - **Naechste dichte Posten** (nicht umgesetzt): `pi_bond_orders` (Dreiecksarray N(N+1)/2
    doubles, ~860 MB bei N=14640, 44 Lesestellen in 3 Dateien) und die dichte topologische
    Abstandsmatrix der EEQ Phase 1 (`computeTopologicalDistancesSparse`, `eeq_solver.cpp:3140`,
    liefert trotz des Namens eine dichte N×N-`Matrix` plus ein N²-float-Arbeitsfeld).
  - **Offene Portierungsfrage, gefunden beim Gegenlesen der Referenz, NICHT geaendert**: die
    H-H-Repulsionsfaktoren `hh13rep`/`hh14rep` liest die Referenz aus `topo%bpair`
    (`gfnff_ini.f90:755-756`), curcuma aus dem BFS-`topo_distances`. Beide stimmen ueberein,
    solange die Nachbarliste symmetrisch ist; ueber eine asymmetrisch gespeicherte eta-Bindung
    (Metall listet das eta-C, das C nicht das Metall) ist der BFS-Abstand richtungsabhaengig und
    kann 2/3 liefern, wo bpair 5 liest (z. B. H am Cp-Ring vs. Hydrid am Metall). Betrifft nur
    eta-Komplexe mit solchen H-H-Paaren; ob MOR41 einen Fall enthaelt, ist nicht gemessen.
    Umstellung auf `bpair` waere ein eigener Commit mit MOR41-Arbitrierung gegen pprcht.
  - **Gemessen per Opt-in `-gfnff.hh_repulsion_bpair true` (2026-09-25, AI-implemented,
    machine-tested; Default unveraendert BFS)**:
    - **Die Asymmetrie ist nicht auf eta beschraenkt**: GMTKN55 aendert sich an 5 von 2462
      Strukturen, alle mit Hauptgruppenmetall (Al/Li/Mg); MOR41 an 0 von 95.
    - **Gegen xtb 6.7.1** (`~/Downloads/xtb-dist`, kein pprcht verfuegbar — kein gfortran auf
      dieser Maschine), curcuma − xtb in kcal/mol, BFS → bpair: `AL2X6/al2me5` +0,0664 →
      **+0,0004**, `ALK8/li2_ch4` +0,0311 → **+0,0001**, `MB16-43/11` +0,0048 → **−0,0001**.
      `MB16-43/23` −83,633 → −83,641 und `/37` −13,519 → −13,522 liegen in der bekannten
      pprcht-vs-xtb-Spaltung und sind so nicht beurteilbar. al2me5 war die groesste verbleibende
      curcuma-vs-pprcht-Abweichung (Known Issue #25: 0,066) — gleiche Groesse wie die Verschiebung,
      Vorzeichen gegen pprcht aber ungemessen.
    - **Der BFS-Default ist atomreihenfolgeabhaengig**: der BFS startet beim kleineren Index, und
      ueber eine einseitig gespeicherte Bindung erreicht er den Partner nur in einer Richtung.
      Idealisiertes CpFe(CO)2H: Hydrid zuerst −3,05533075, Hydrid zuletzt −3,05561154 Eh
      (0,18 kcal/mol, nur Nonbond-Repulsion); mit bpair beide −3,05561154. Das Molekuel taugt nur
      als Schalterkontrolle — curcuma liegt dort 150–177 kcal/mol neben xtb, und xtb selbst gibt
      fuer die zwei Reihenfolgen 23,6 kcal/mol verschiedene Energien.
    - **Gegen pprcht/gfnff** (`external/gfnff` @ `0491df2f`, `-Dbuild_exe=ON`, gfortran 16.2.1;
      Provenienz: `AHB21/21` −2,027888583 = dokumentiert −2,027889), curcuma − pprcht in kcal/mol,
      BFS → bpair: al2me5 +0,0664 → **+0,0005**, li2_ch4 +0,0311 → **+0,0001**, MB16-43/23
      +0,0082 → **+0,0008**, MB16-43/11 +0,0048 → **−0,0001**, CpFe(CO)2H Hydrid zuerst
      +0,1763 → **+0,0001** (Hydrid zuletzt +0,0001 beide; pprcht selbst ist
      reihenfolgeunabhaengig, curcuma trifft es — die 150–177 kcal/mol oben waren reine
      xtb-Abweichung). MB16-43/37 +0,0001 → −0,0023: dort wird der Repulsionsterm mit bpair
      **exakt** (−6e-8 Eh gegen +3,8e-6 mit BFS); die Gesamtuebereinstimmung mit BFS war
      Fehlerkompensation mit einer Dispersionsabweichung von −3,7e-6 Eh, die jetzt sichtbar ist
      (alle anderen Terme ≤5e-5 kcal/mol). Diese Dispersionsdifferenz ist ein eigener, offener
      Kleinstbefund, nicht Folge des Schalters.
    - **Erledigt (2026-09-25)**: `hh_repulsion_bpair` ist Default. Die dabei sichtbar gewordene
      MB16-43/37-Dispersionsdifferenz war curcumas eigener ATM-Dreikoerperterm (in der Referenz
      nicht vorhanden) — jetzt aus (`dispersion_atm`), MOR41+GMTKN55 vs pprcht MAD 0,00014 → 0,00007,
      Strukturen >0,001 kcal/mol 58 → 8. Ein anschliessender Permutationstest fand den
      Amid-H-Reihenfolgefehler der Referenz (curcuma korrekt, Opt-in `amideh_acidity_order_bug`)
      und 911/2557 Strukturen, deren Energie in pprcht UND curcuma von der Atomnummerierung
      abhaengt — Details: CLAUDE.md Known Issue #32, `docs/REV_GFNFF_TODO.md` #11-#13.
    - **Offen**: Reihenfolgeabhaengigkeit der Referenz (REV_GFNFF_TODO #12: Ladungsplatzierung,
      Inversionsterm, Coulomb bei neutralen Molekuelen, Bindungsterm MB16-43); verbleibende
      Portierungsreste vs pprcht >0,01 kcal/mol: MB16-43/04 +0,0445, HEAVYSB11/pbme3 −0,0389,
      MB16-43/43 +0,0137, MB16-43/01 +0,0131.
- **`cli_curcumaopt_07_opt_multixyz` repariert (2026-09-25)**: der Test verglich den Multi-XYZ-Pfad
  mit einer Golden-Datei vom Juni 2026, die seitdem mit jeder GFN-FF-Korrektur gedriftet war
  (Frame 02/12 lagen 0,15–0,17 Eh ueber dem echten Minimum, weitere um bis zu 1,9e-5). Einzel- und
  Multi-XYZ-Optimierung stimmen mit dem aktuellen Binary fuer alle 17 Frames auf <=1e-6 Eh ueberein
  — der Pfad war nie falsch. Der Test optimiert jetzt jeden Frame selbst einzeln (+6 s) und
  vergleicht dagegen, `golden_energies.txt` ist entfernt. Dadurch laeuft erstmals auch der zweite
  Durchlauf mit `-threads 4` (wurde nach dem ersten Fehlschlag nie erreicht): 40/40. Negativkontrolle
  (eine Referenz um 1e-4 Eh verfaelscht) schlaegt an. `ctest -L gfnff` jetzt **78/78**.

### SIGSEGV am Ursprung untersucht (Auftrag „fix den SIGSEGV am Ursprung") — nicht gefunden, Werkzeuge sind blind dafuer (2026-09-24)
- **Status**: ⏳ OFFEN. Root Cause NICHT gefunden trotz gruendlicher Untersuchung mit ASan,
  compute-sanitizer und gdb. Kein Fix umgesetzt.
- **Eigener ASan-Build von Grund auf gebaut** (`build_asan/`, existierte vorher nicht), passend
  zu `release/`'s CUDA/GFN-FF/BLAS-Konfiguration. **Dabei einen echten, aber unabhaengigen
  Werkzeugkettenfehler gefunden und aussortiert**: Eigens AVX-512/AVX2-`pstore`-Intrinsics
  crashen unter dieser GCC-16.1.1+ASan-Kombination deterministisch bei JEDER vektorisierten
  `Eigen::Vector`-Zuweisung — reproduzierbar auf einem 24-Atom-System, reinem `-sp -gpu cuda`,
  ganz ohne MD oder Schwellenwert-Flag. Bestaetigt abwesend im normalen `release/`-Binary,
  bestaetigt weg bei `-O0` und mit `-DEIGEN_DONT_VECTORIZE`. Kein curcuma-Bug, nicht
  weiterverfolgt — aber **jeder ASan-Befund aus diesem Baum braucht ab sofort eine Gegenprobe
  gegen das normale Binary**, bevor er geglaubt wird.
- **Mit deaktivierter Vektorisierung (ASan nutzbar) reproduziert sich der ECHTE Bug weder unter
  ASan noch unter `compute-sanitizer --tool memcheck`**: 5/5 bzw. 4/5-abstuerzend-aber-0-Fehler
  im selben Testlauf, der das normale Binary 4/5-mal crashen laesst. Das ist ein **aussagekraeftiges
  Negativergebnis**: beide Werkzeuge sind genau fuer diese Fehlerklassen gebaut (ASan fuer
  Host-Heap-Missbrauch, compute-sanitizer fuer CUDA-Host/Device-Missbrauch) und finden nichts,
  obwohl der Absturz im selben Lauf mit aehnlicher Rate auftritt.
- **Ausgeschlossen, jeweils mit Gegenprobe**: CPU-seitige Threading-Race in der O(N^2)
  EEQ-Matrix-Fuellung (std::thread-parallel ab 64 Atomen) — crasht identisch bei `-threads 1`;
  GPU-Kernel-seitiger Fehler — `info threads` am Absturzpunkt zeigt genau EINEN Thread; CLI-
  Parsing des langen `-eeq_rocm_cpu_fragment_threshold`-Flag-Namens — crasht identisch, wenn
  derselbe Wert per `-import_config`-JSON statt CLI gesetzt wird; Groessen-Mismatch im
  wiederverwendeten `m_phase2_A`-Puffer — per Breakpoint bestaetigt, dass `ensurePhase2Buffers`
  bei jedem Aufruf die korrekten, gecachten `natoms`/`nfrag` sieht. **Vor allem: ein reiner
  CPU-Lauf derselben MD (`-threads 8`, kein `-gpu`) crasht nie (0/3)** — der Fehler braucht
  zwingend einen aktiven CUDA-Kontext.
  Kernel-Log (`dmesg`/Xid): kein NVIDIA-Treiberfehler zu irgendeinem Absturz protokolliert
  (uneindeutig — Xid deckt GPU-Kernel-/Reset-Fehler ab, nicht notwendigerweise Userspace-
  CUDA-Bibliotheks-Heap-Interaktionen).
- **Echter Fund, keine reine Ausschlussliste**: gdb am `this`-Zeiger des allerersten
  `EEQSolver::ensurePhase2Buffers`-Aufrufs nach dem Absturz-ausloesenden CUDA-Sync
  (`FFWorkspaceGPU::finalizeCNForCPU`s `cudaStreamSynchronize`) zeigt echte Korruption — kein
  Debugger-Artefakt, denn ein nachfolgendes `print m_phase2_buf_natoms` ueber denselben Zeiger
  schlaegt mit „Cannot access memory" fehl. Der korrupte Zeigerwert dekodiert zu plausiblen
  ASCII-Fragmenten; das war NICHT reproduzierbar als festes Muster (aendert sich von Lauf zu
  Lauf) und die naheliegende Spur (ein langer CLI-Flag-Name) hielt der Gegenprobe oben nicht
  stand — als zufaelligen Byteinhalt IRGENDEINES nahen Heap-/Rodata-Strings behandeln, nicht
  als geloeste Spur.
- **Wo das F-Q9/D-26/D-46 zuruecklaesst** (`docs/TECHNICAL_DEBT.md`): „Root cause
  uninvestigated" stimmt nicht mehr — es wurde mit den Standardwerkzeugen untersucht, bis
  bestaetigt war, dass diese Werkzeuge es nicht sehen, und der Ausloeser ist eingegrenzt auf
  „CPU beruehrt einen grossen heap-residenten Puffer kurz nach einem `cudaStreamSynchronize`,
  im selben Prozess wie ein aktiver CUDA-Kontext" (deckt sich direkt mit der bestehenden
  `gfnff.h:1128`-Charakterisierung „CUDA corrupts heap metadata"). Genauer geht es ohne
  NVIDIA-interne Treiber-Werkzeuge nicht, die hier nicht verfuegbar sind. **Kein Fix
  umgesetzt** — das bereits vorhandene Umgehungsmuster (vorallozierte Puffer, memcpy statt
  Eigen-Zuweisung, keine frische Heap-Aktivitaet nahe einem CUDA-Sync) ist die einzige bekannte
  Gegenmassnahme; es auf `EEQSolver`s Phase-2-Pfad auszuweiten ist der naechste konkrete,
  abgegrenzte Schritt, kein Ursprungsfix.
- **Verweis**: `docs/TECHNICAL_DEBT.md` F-Q9 (dortiger datierter Eintrag mit vollen Details),
  `Labor/curcuma MD-Stabilität großer Systeme.md` (Vault) fuer das vollstaendige Sitzungsprotokoll
  mit jedem verworfenen Hebel und dem exakten Kommando dazu.
- **Aufraeumen**: `build_asan/` (4,1 GB, nicht versioniert) steht noch, fuer den Fall, dass die
  Untersuchung fortgesetzt wird. Kann geloescht werden, wenn nicht.

### `gpu_strict` fehlt — stiller CPU-Rueckfall (2026-09)
- **Status**: ⏳ PLANNED
- **Problem**: `xtb_gpu_context.cu:3875` warnt nur, wenn eine Rechnung nicht auf die Karte
  passt, und rechnet dann auf der CPU weiter. Bei Laufzeitmessungen auf fremder Hardware
  kostet das Stunden, bevor es auffaellt.
- **Task**: PARAM `gpu_strict` (Bool, default false) — GPU-Fehler/OOM als harter Fehler.
  Im Multi-GPU-Plan vorgesehen, nie implementiert.
- **Verweis**: [docs/MULTI_GPU.md](docs/MULTI_GPU.md):88, Testprotokoll H200

---

## 🔵 CAPABILITIES & ANALYSIS (src/capabilities/)

### ConfScan Verbosity Enhancement
- **Status**: ⏳ PENDING
- **Problem**: Accept/Reject messages not visible at default verbosity level
- **Task**: Adjust CurcumaLogger calls to be visible at level ≥1
- **Betroffene Dateien**: src/capabilities/confscan.cpp
- **Verweis**: src/capabilities/CLAUDE.md:84

### SimpleMD Wall Potential Physics
- **Status**: ⏳ PENDING (LOW PRIORITY)
- **Task**: Verify wall potential physics (boundary logic, force calculations)
- **Betroffene Dateien**: src/capabilities/simplemd.cpp/h
- **Verweis**: src/capabilities/CLAUDE.md:85
- **Note**: Tests now passing (7/7) - functionality works, physics validation pending

### RMSD Strategy Pattern - Phase 3
- **Status**: ⏳ PENDING
- **Task**: Complete Strategy pattern refactoring for RMSD module
- **Betroffene Dateien**: src/capabilities/rmsd.cpp/h
- **Verweis**: src/capabilities/CLAUDE.md:80
- **Phases**: Phase 1-2 ✅, Phase 3 ⏳

### ConfSearch: CxxThreadPool Crash mit threads=1 + Progress-Bar
- **Status**: ⏳ PENDING
- **Problem**: CxxThreadPool im Legacy-Modus crashed bei threads=1 (SIGSEGV nach Optimization). Discrete-Progress-Bar funktioniert nur im Legacy-Modus (ParallelLoop ruft updateStatus()), nicht im Worker-Pool-Modus.
- **Task**: 
  1. Crash in CxxThreadPool Legacy-Modus bei threads=1 debuggen und fixen
  2. Alternative: Progress-Bar im Worker-Pool-Modus aktivieren (updateStatus() während wait() poll)
  3. Danach: Serial-Path in PerformOptimisation() entfernen, unified CxxThreadPool für alle threads
- **Betroffene Dateien**: src/capabilities/confsearch.cpp, external/CxxThreadPool/include/CxxThreadPool.hpp
- **Verweis**: confsearch-branch, May 2026
- **Workaround**: Aktuell Serial-Path für threads <= 1, Parallel-Path (Legacy-Modus) für threads > 1

### ConfSearch: GPU + Multi-Threading (Future)
- **Status**: ⏳ PLANNED
- **Problem**: Bei threads > 1 konkurrieren mehrere MD-Instanzen um die GPU. Aktuell wird GPU deaktiviert wenn threads > 1.
- **Task**: 
  1. Implementiere GPU-Lock (Mutex) oder Queue, sodass nur 1 Thread gleichzeitig die GPU nutzt
  2. Alternative: Ein Thread bekommt GPU-CUDA, andere nutzen CPU-Fallback
  3. Dann kann GPU auch bei threads > 1 aktiviert werden
- **Betroffene Dateien**: src/capabilities/confsearch.cpp, src/core/energycalculator.cpp
- **Verweis**: confsearch-branch, May 2026
- **Workaround**: Aktuell wird GPU auf "none" gesetzt wenn threads > 1

### Enhanced Conformational Search Algorithms
- **Status**: ⏳ PLANNED
- **Betroffene Dateien**: src/capabilities/confsearch.cpp/h
- **Verweis**: src/capabilities/CLAUDE.md:76

### Improved Trajectory Analysis Tools
- **Status**: ⏳ PLANNED
- **Betroffene Dateien**: src/capabilities/rmsdtraj.cpp/h
- **Verweis**: src/capabilities/CLAUDE.md:77

---

## 🟣 MOLECULE REFACTORING - Test-Driven Breaking Changes

**Status**: ✅ Phase 1 DONE, ⏳ Phase 2-6 PENDING
**Test-Driven Approach**: src/core/test_molecule.cpp (15 test categories)
**Critical Requirement**: All existing functionality must remain API-compatible

### Phase 2: XYZ Comment Parser Unification
- **Status**: ⏳ PENDING
- **Task**: Eliminate 10 duplicate XYZ parser functions
- **Betroffene Dateien**: src/core/molecule.cpp/h
- **Critical Constraint**: Production comment formats must NOT break (ORCA, XTB, simple energy)
- **Reference**: docs/XYZ_COMMENT_FORMATS.md
- **Verweis**: CLAUDE.md - Planned Development

### Phase 3: Granular Cache System
- **Status**: ⏳ PLANNED
- **Task**: Replace single `m_dirty` flag with fine-grained cache invalidation
- **Betroffene Dateien**: src/core/molecule.cpp/h
- **Performance Impact**: Selective recalculation instead of full cache invalidation

### Phase 4: Fragment System O(1) Lookups
- **Status**: ⏳ PLANNED
- **Task**: Replace std::map with optimized lookup structure
- **Betroffene Dateien**: src/core/molecule.cpp/h

### Phase 5: Type-Safe ElementType Enum
- **Status**: ⏳ PLANNED
- **Task**: Replace int element indices with ElementType enum
- **Betroffene Dateien**: src/core/molecule.cpp/h, src/core/units.h

### Phase 6: Unified Atom Structure with Zero-Copy Geometry
- **Status**: ⏳ PLANNED
- **Task**: Implement SOA/AOS hybrid design for geometry access
- **Betroffene Dateien**: src/core/molecule.cpp/h

---

## 📊 PRIORITY MATRIX

| Priority | Component | Count | Status |
|----------|-----------|-------|--------|
| 🔴 KRITISCH | None | 0 | ✅ ALL RESOLVED |
| 🟢 DONE (CG Phase 5) | SimpleMD CG Integration | 1 | ✅ COMPLETE |
| 🟢 DONE (CG Phases 1-4) | CG Core + VTF + Testing | 4 | ✅ COMPLETE |
| 🔵 LOW (CG Phase 6) | Ellipsoidal Extensions | 1 | 🟡 PREPARED |
| 🟡 TESTING | Scientific validation | 4 | ⏳ PENDING |
| 🟢 CORE | Parameter/Memory/Units | 4 | ⏳ PENDING |
| 🔵 CAPABILITIES | Confscan/RMSD/SimpleMD Physics | 4 | ⏳ PENDING |
| 🟣 REFACTORING | Molecule Phase 2-6 | 5 | ⏳ PLANNED |
| **TOTAL** | | **23** | **5 ✅ + 16 ⏳ + 2 🟡** |

---

## 🔧 TESTING STATUS BY MODULE

- **RMSD**: 6/6 ✅ (100%)
- **ConfScan**: 7/7 ✅ (100%)
- **CurcumaOpt**: 6/6 ✅ (100%)
- **SimpleMD**: 7/7 ✅ (100%) - **FIXED October 28, 2025**
- **Overall**: 26/26 ✅ (100%) 🎯

---

**Last Updated**: 2025-11-14
**Next Review**: When starting CG implementation (Phase 1: Molecule helpers) or cgfnff debugging

## Native QM Methods (GFN2, GFN1, PM3) - November 2025

### ✅ Completed (Phase 1-4)

- [x] GFN2-xTB native implementation structure
- [x] GFN2 core algorithms (CN, Hamiltonian, SCF, energies)
- [x] GFN2 analytical gradients (Electronic, Repulsion, Coulomb, CN) - ✅ **NEW Feb 2026**
- [x] GFN1-xTB native implementation with halogen bond correction
- [x] PM3 NDDO implementation (H, C, N, O)
- [x] Integration into MethodFactory with priority fallbacks
- [x] Educational documentation (NATIVE_QM_IMPLEMENTATION_STATUS.md)
- [x] **Parameter loader infrastructure** (`gfn2_params_loader.h/cpp`) with real TBLite parameters for all 86 elements - ✅ **NEW Feb 2026**
- [x] **Ulysses methods documentation** - Complete guide for 27 semi-empirical methods (AM1, MNDO, PM6, etc.)
- [x] **Overlap Derivatives** analytical implementation in `STOIntegrals.hpp` - ✅ **NEW Feb 2026**

### 🔧 TODO: Parameter Expansion (Medium Priority)

**Files**: `gfn2_params_loader.cpp`, `gfn2_params_loader.h`
**Status**: ✅ Infrastructure and 86-element DB complete, ⏳ Extension needed

**GFN2 Real Parameters from TBLite** (Foundation ✅, Extension needed):
- [x] ✅ Parameter loader class structure (`ParameterDatabase`)
- [x] ✅ Shell-resolved parameter structures (`ShellParams`, `ElementParams`, `PairParams`)
- [x] ✅ Hardcoded real parameters for all 86 elements (basic set) - ✅ **UPDATED Feb 2026**
- [x] ✅ Element-pair specific Hamiltonian scaling (C-H, C-C, C-N, C-O, N-H, O-H)
- [x] ✅ Exact shell-Hubbard corrections from TBLite (SHELL_HUBBARD_CORR) - ✅ **NEW March 2026**
- [x] ✅ D4 dispersion integration (USE_D4 conditional) - ✅ **NEW March 2026**
- [ ] ⏳ Full TOML parser implementation (currently: stub)
- [x] ✅ Complete periodic table coverage (86 elements) - ✅ **UPDATED Feb 2026**
- [ ] ⏳ Extract polynomial corrections poly(r) for all pairs
- [ ] ⏳ Extract complete gamma-AB Coulomb kernel parameters
- [ ] ⏳ Extract full AES2 multipole parameters (quadrupoles)
- [ ] ⏳ Validate against TBLite reference energies (<1% error target)

**GFN1 Real Parameters from TBLite**:
- [ ] Create `gfn1_params_loader.h/cpp` (analogous structure)
- [ ] Extract simplified parameter set (no ES3)
- [ ] Extract halogen bond parameters for F, Cl, Br, I, At
- [ ] Validate against TBLite GFN1 reference

**Implementation Notes**:
- ✅ Foundation in place: See `gfn2_params_loader.cpp:41-318`
- ✅ Educational transparency: Comments explain each parameter's physical meaning
- ⏳ Next step: Implement `parseSimpleTOML()` or use external TOML library (toml11)
- ⏳ Alternative: Continue hardcoding parameters from TBLite source for remaining elements

### 🎯 TODO: PM3 Element Extension (Medium Priority)

**Critical Missing Elements**:
- [ ] F (Fluorine) - Important for pharmaceuticals
- [ ] Cl (Chlorine) - Common in organic chemistry
- [ ] S (Sulfur) - Biochemistry (cysteine, methionine)
- [ ] P (Phosphorus) - DNA, ATP, phosphates

**Parameter Source**: MOPAC parameter database
- URL: http://openmopac.net/manual/parameters.html
- Format: U_ss, U_pp, beta_s, beta_p, zeta_s, zeta_p, alpha, Gaussian terms

**Extended Elements** (Lower Priority):
- [ ] Br, I (heavier halogens)
- [ ] Si, Se (semiconductors, proteins)
- [ ] Transition metals (Fe, Cu, Zn for catalysis)

### ⚡ TODO: Performance Improvements (Optional)

**Analytical Gradients** (Speedup: 10-20x for optimizations):
- [x] Implement Hellmann-Feynman theorem derivatives - ✅ **DONE Feb 2026**
- [x] GFN2: dH/dR, dS/dR analytical formulas (STOIntegrals) - ✅ **DONE Feb 2026**
- [ ] GFN1: Simplified derivative terms
- [ ] PM3: NDDO gradient formulas from MOPAC
- [ ] Benchmark: H2O optimization (numerical vs analytical)

**SCF Acceleration**:
- [ ] DIIS (Direct Inversion in Iterative Subspace)
- [ ] Level shifting for difficult convergence
- [ ] Adaptive damping parameters

### 🧪 TODO: Validation and Testing

**Test Molecules** (Small):
- [ ] H2O - HOMO/LUMO, dipole moment
- [ ] CH4 - Symmetry, C-H bonds
- [ ] NH3 - Lone pair, pyramidal geometry
- [ ] H2CO - Carbonyl, planarity

**Comparison Targets**:
- [x] GFN2 native vs TBLite (energy error < 1% with real params) - D4 + Shell-Hubbard completed March 2026
- [ ] GFN1 native vs TBLite (energy error < 1%)
- [ ] PM3 native vs MOPAC (energy error < 5%)

**Properties to Validate**:
- [ ] Total energy (Hartree)
- [ ] HOMO/LUMO gap (eV)
- [ ] Mulliken charges
- [ ] Dipole moment (Debye)
- [ ] Gradient accuracy (force = -gradient)

### 📋 TODO: Integration Tasks

**D3/D4 Dispersion**:
- [x] GFN2 D4 integration (March 2026) - calculateDispersionEnergy() implemented
- [ ] Fix `dftd3interface.h/cpp` issues
- [ ] Fix `dftd4interface.h/cpp` issues
- [ ] Connect GFN1 to D3 (replace stub)
- [ ] Connect GFN2 to D4 (replace stub)
- [ ] Validate dispersion energies for Ar2, benzene dimer

**Heat of Formation** (PM3-specific):
- [ ] Implement ΔH_f calculation from atomization energy
- [ ] Add experimental reference data for validation
- [ ] Compare with MOPAC heats of formation

### 📚 Documentation TODO

- [ ] Add example usage to CLAUDE.md for each method
- [ ] Create tutorial: "When to use GFN2 vs GFN1 vs PM3"
- [ ] Document parameter extraction workflow
- [ ] Add benchmark comparison tables (TBLite, MOPAC, Ulysses)

### 🚫 NOT TODO (Use Existing Interfaces)

**These methods already available via Ulysses - no native implementation needed**:
- ❌ AM1 - Available: `method = "am1"` (Ulysses)
- ❌ MNDO - Available: `method = "mndo"` (Ulysses)
- ❌ PM6 - Available: `method = "pm6"` (Ulysses)
- ❌ RM1 - Available: `method = "rm1"` (Ulysses)
- ❌ PM3PDDG - Available: `method = "pm3pddg"` (Ulysses)

**Reason**: Ulysses interface provides production-quality implementations with full validation. Native implementations only needed for educational purposes or when external dependencies unavailable.

### 📊 Status Summary

**Implementation**: ✅ 100% Complete (3 methods, ~2767 lines)
**Integration**: ✅ 100% Complete (MethodFactory, CMakeLists.txt)
**Parameters**: ⚠️ 30% Complete (approximations work, real params TODO)
**Validation**: ⏸️ 0% Complete (needs test molecules)
**Performance**: ⚠️ 50% Complete (numerical gradients slow, analytical TODO)

**Next Recommended Action**: Extract GFN2 parameters from TBLite TOML for production accuracy.

---

*Last Updated: November 2025*
*See: docs/NATIVE_QM_IMPLEMENTATION_STATUS.md for full details*

---

## 🔴 BUILD SYSTEM - CONDITIONAL COMPILATION FIXES (November 2025)

### Status: PARTIALLY FIXED - Build 2 Working, Others Need Deeper Fixes

**Session Date**: November 2025
**Task**: Test all 5 build configurations with different CMake options (USE_D3, USE_D4, USE_TBLITE, USE_ULYSSES, USE_XTB)

### Build Test Results

| Build | Config | Status | Notes |
|-------|--------|--------|-------|
| **Build 1** | Minimal (UFF, EHT only) | ❌ FAILED | `s-dftd3.h` not found - gfnff.cpp includes D3 unconditionally |
| **Build 2** | Standard (TBLite, Ulysses, D3) | ✅ **SUCCESS** | ✅ Verified working: UFF & GFN2 methods tested |
| **Build 3** | Full QM (+ XTB) | ❌ FAILED | Same D3 header chain issue |
| **Build 4** | D4 Dispersion | ❌ FAILED | D3 dependency blocks build |
| **Build 5** | TBLite only | ❌ FAILED | D3/D4 guard chain |

### Root Cause Analysis

**Problem**: Cascading `#ifdef` guards only protect headers, not implementations:
```
gfnff.cpp (always compiles)
  → #include "forcefield.h"
    → #include "forcefieldthread.h"
      → #ifdef USE_D3 #include "dftd3interface.h" #endif
        → #include "s-dftd3.h"  (UNGUARDED INCLUDE!)
```

When `USE_D3=OFF`, compiler skips the `#ifdef` but still tries to compile the file, causing `s-dftd3.h` not found error.

### Fixes Applied ✅

1. **`dftd3interface.h`** - Wrapped entire header with `#ifdef USE_D3...#endif`
2. **`dftd4interface.h`** - Wrapped entire header with `#ifdef USE_D4...#endif`
3. **`forcefieldthread.h`** - Wrapped `D3Thread` class definition with `#ifdef USE_D3...#endif`
4. **`forcefieldthread.cpp`** - Wrapped D3Thread implementation with `#ifdef USE_D3...#endif`

### Remaining Work (To Fix Builds 1, 3, 4, 5)

**Critical files still needing protection**:
- [ ] `gfnff.cpp` - Always compiled, needs conditional compilation or restructuring
- [ ] `forcefield.cpp` - D3Thread instantiation (line ~470) needs `#ifdef USE_D3` guard
- [ ] `gfnff.h` - Consider lazy-loading or factory pattern for D3/D4 dependencies
- [ ] Add H4Thread guards (same pattern as D3Thread)

### Build 2 (Standard) - PRODUCTION READY ✅

**Verified Methods:**
```bash
./curcuma -sp water.xyz -method uff    # ✅ Works
./curcuma -sp water.xyz -method gfn2   # ✅ Works (via TBLite)
```

**Available Methods in Build 2:**
- UFF (universal force field) - native
- EHT (extended Hückel theory) - native
- GFN2 (tight-binding DFT via TBLite) - recommended
- GFN1 (TBLite or Ulysses fallback)
- iPEA1 (TBLite)
- PM6, PM3, AM1, MNDO (Ulysses semi-empirical)

### Recommendations

**Short term**: Protect gfnff.cpp and forcefield.cpp with conditional compilation
**Medium term**: Refactor include chain - move D3/D4/H4 threads to separate file
**Long term**: Use CMake-level validation + CI/CD pipeline for all build configurations

### Test Infrastructure

Created `/home/conrad/src/curcuma/build_test/` with 5 isolated builds for regression testing.

See `docs/BUILD_SYSTEM_ROADMAP.md` for detailed implementation plan.

---

**Last Updated**: 2025-11-08 (Build System Testing)
**Next Review**: After implementing remaining Build 1/3/4/5 fixes or when starting CG Phase 5
