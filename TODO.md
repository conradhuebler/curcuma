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
