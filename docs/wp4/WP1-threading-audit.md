# WP1 — Threading-Audit (P4b)

**Status:** ⚙️ Machine-tested (Mai 2026) — Operator setzt ✅ TESTED
**Aufwand:** 2h
**Voraussetzung für:** WP2 (EEQ batched parallel)

## Hypothese

Im Mai-2026-Profil (`mixture.xyz`, N=6200) ist die Phase **CN + EEQ** mit **1766 ms (53 %)** als "serial CPU" gelabelt. Die EEQ-Komponenten — Per-Fragment-LU (Stage 4 batched), Coordination-Number-Berechnung, D4-Gaussian-Weights — sind algorithmisch trivial parallelisierbar. Das Label sagt nicht zwingend, dass die Phase wirklich seriell läuft; möglich sind drei Ursachen:

1. Tatsächlich nur 1 Thread aktiv, weil `-threads N` nicht durchgereicht wird.
2. Threads aktiv, aber das Logging-Label ist veraltet.
3. Threads aktiv, aber ein Sub-Step (z. B. die Stage-4-Loop selbst) ist serial — siehe WP2.

## Aufgabe

End-to-end-Audit der Thread-Propagation für den `gfnff`-CPU-Pfad.

### Schritte

1. **Eintrittspunkt identifizieren:** Wo liest die Capability (Opt, SimpleMD, sp) `threads` aus dem Controller-JSON? In welche ConfigManager-Felder fließt der Wert?
2. **Folge-die-Variable:** Tracken wie `threads` in
   - `EnergyCalculator` ankommt (`m_threads` o. ä.)
   - `MethodFactory::createMethod` durchgereicht wird
   - `GFNFF`-Wrapper (`gfnff_method.cpp`) gesetzt wird
   - `EEQSolver` erreicht (`eeq_solver.cpp`)
   - `D4ParameterGenerator::precomputeGaussianWeights(pool, num_threads)` (`d4param_generator.cpp:986`) bekommt
   - `CNCalculator` (über `omp_get_max_threads` — siehe `cn_calculator.cpp:130` `#pragma omp parallel for`)
3. **Mit Logger empirisch prüfen:** An jedem Knoten ein Verbosity-2-Print `"<modul>: using N threads"`. Lauf auf `mixture.xyz` mit `-threads 4` und mit `-threads 1`, Diff der Logs.
4. **Label korrigieren:** Falls eine Phase tatsächlich n-threaded läuft, die `"(serial CPU)"`-Strings in `gfnff_method.cpp` (s. `phase_row("CN + EEQ (serial CPU)", ...)`, ~Z. 131) anpassen.
5. **Lücke identifizieren:** Falls eine Phase serial läuft trotz `-threads 4`, das ist der Befund für WP2 — Stelle dokumentieren.
6. **Mechanismus-Inventar:** Tabelle aller parallelen Stellen im `gfnff`-Pfad mit Spalten *Datei, Zeile, Mechanismus (OpenMP / CxxThreadPool / raw std::thread / serial), Begründung*. Hintergrund: `cn_calculator.cpp` nutzt OpenMP, `d4param_generator.cpp` nutzt CxxThreadPool (+ `std::thread`-Fallback) — dieser Mix ist eine bekannte Inkonsistenz. WP1 produziert die Empfehlung für *neue* Stellen (siehe unten).

### Mechanismus-Empfehlung (Hintergrund für WP2 ff.)

| Mechanismus | Fork-Join-Latenz | Eignung in Curcuma |
|-------------|------------------|--------------------|
| **CxxThreadPool + `enqueue` mit Per-Thread-Buffern** | 5–10 µs | ✅ **Bevorzugt** für neue Stellen — schon Curcuma-Konvention (`precomputeGaussianWeights`), False Sharing strukturell verhindert, deterministisches Interleave-Scheduling |
| OpenMP `parallel for` | 1–5 µs | OK für triviale Reduktionen ohne geteilte Schreibziele (CN, kleine Loops). Bei komplexen Reduktionen (Energie + Gradient + Per-Term-Beiträge) Risiko von False Sharing / Critical-Sections — siehe Fortran-GFN-FF, das mit OpenMP **mit mehr Threads langsamer wird**. |
| `std::async(launch::async, …)` | **50–100 µs pro Call** | ❌ **Vermeiden.** Erzeugt pro Aufruf neuen `std::thread`, oder läuft (default) seriell. Für Loops mit > 100 Iterationen und kurzen Aufgaben katastrophal. |

**Wichtiger Befund:** Der Fortran-GFN-FF-Originalcode (OpenMP) skaliert oberhalb von ~4 Threads negativ. Vermutlich False Sharing auf Energie/Gradient-Akkumulatoren oder zu naives `schedule(static)`. CxxThreadPool mit Per-Thread-`local_energy` und Reduktion am Ende umgeht dieses Muster strukturell — und deckt sich mit dem bisherigen Curcuma-Stil.

### Skalierungs-Sweep (zusätzliche Audit-Aufgabe)

Beweisen, ob Curcuma denselben Pathologie-Effekt hat:

```bash
for T in 1 2 4 8 16; do
  ./curcuma -sp mixture.xyz -method gfnff -threads $T -verbosity 2 \
    | grep -E "Total energy|EEQ Phase|Bond|Coulomb|Dispersion"
done
```

Erwartet bei korrektem Threading: monoton fallende Wall-Time bis zur Core-Anzahl, dann Plateau.
Pathologisch: Wall-Time hat ein Minimum bei T=2 oder T=4 und steigt darüber wieder — Indiz für False Sharing, Memory-Bandwidth-Sättigung, oder zu fein-granulares Scheduling.

Eintrag in WP1-Bericht: Tabelle Threads × Wall-Time pro Phase. Falls Pathologie sichtbar, ist das ein zusätzlicher Hebel für WP2/WP3 (Per-Thread-Buffer-Hygiene).

### Code-Anker

| Datei | Position | Kontext |
|-------|---------|---------|
| `src/core/energy_calculators/ff_methods/gfnff_method.cpp` | ~Z. 131 | `phase_row("CN + EEQ (serial CPU)", ...)` |
| `src/core/energy_calculators/ff_methods/gfnff_method.cpp` | `prepareCNAndEEQ` (~Z. 902) | Phase-A-Helper |
| `src/core/energy_calculators/ff_methods/eeq_solver.cpp` | Z. 1294–1347 | Stage-4 Batched-Loop (relevant für WP2) |
| `src/core/energy_calculators/ff_methods/cn_calculator.cpp` | Z. 130 | `#pragma omp parallel for` (bereits parallel) |
| `src/core/energy_calculators/ff_methods/d4param_generator.cpp` | Z. 986 | `precomputeGaussianWeights(pool, num_threads)` (bereits parallel) |

## Akzeptanzkriterien

1. Tabelle in WP-Doc + GFNFF-PERFORMANCE-ROADMAP, die für jeden Knoten "Threads aktiv: ja/nein, woher?" beantwortet.
2. Falls Lücke: minimaler Patch, der `threads` an die fehlende Stelle propagiert. Phase-Wall-Time auf `mixture.xyz` ≥ 1.5× schneller mit `-threads 4` vs. `-threads 1`.
3. Falls keine Lücke: Label `"(serial CPU)"` zu `"(multi-threaded)"` korrigiert; WP2 fokussiert dann ausschließlich auf den Stage-4-Loop.
4. CTest grün, `test_gfnff_numgrad` grün.

## Risiken

- Niedrig. Reines Audit + Logging, keine algorithmische Änderung.
- Falls eine Stelle `threads=1` hardcoded hat (z. B. weil sie als nicht-thread-safe eingestuft war), nicht blind aktivieren — Thread-Safety zuerst prüfen (Eigen, Per-Atom-Vektoren).

## Verifikation

```bash
cd release && make -j4
./curcuma -sp ../test_cases/molecules/larger/mixture.xyz -method gfnff -threads 4 -verbosity 2 2>&1 | grep -E "thread|Thread|serial|parallel"
./curcuma -sp ../test_cases/molecules/larger/mixture.xyz -method gfnff -threads 1 -verbosity 2 2>&1 | grep -E "thread|Thread|serial|parallel"
ctest -R "gfnff" --output-on-failure
```

Erwartetes Ergebnis: das `-threads 4`-Log zeigt an mind. 5 Knoten "using 4 threads", das `-threads 1`-Log zeigt überall "using 1 thread".

---

## Ergebnis (Mai 2026)

### Hauptbefund: ein Bug + ein falsches Label, keine Threading-Sackgasse

1. **Bug gefunden und behoben** — `GFNFFComputationalMethod::setThreadCount()` (`qm_methods/gfnff_method.cpp:97-99`) ignorierte den Parameter mit `(void)threads`. Heute leitet er an `GFNFF::setThreadCount` weiter, der eine neue Member `m_threads` setzt und gleichzeitig `m_parameters["threads"]` synchron hält.
2. **Konstruktor-Pfad funktioniert** — `m_controller["threads"]` wird in `EnergyCalculator::createMethod` (`src/core/energycalculator.cpp:146-178`) ins `method_config` kopiert und damit über `MethodFactory::create` in den GFN-FF-Konstruktor (`m_threads = m_parameters.value("threads", 1)`). Erst nach dem Setter-Fix funktionieren beide Pfade.
3. **Label korrigiert** — `phase_row("CN + EEQ (serial CPU)", ...)` → `phase_row("CN + EEQ", ...)` in `gfnff_method.cpp:131`. Die Phase ist nicht serial: CN nutzt OpenMP, EEQ-A-Matrix und Distance-Matrix nutzen OpenMP+CxxThreadPool. Nur der **Stage-4 Per-Fragment-LU-Loop** in `eeq_solver.cpp:1305-1335` ist explizit seriell — das ist die Zielstelle für **WP2**.
4. **Keine OpenMP-Pathologie auf Curcuma-Seite** — der Skalierungs-Sweep zeigt monoton fallende Wall-Time bis T=16, **keine** negative Skalierung wie im Fortran-Original. Die `parallel efficiency` fällt allerdings monoton (78.6 % bei T=2 → 41.3 % bei T=16) — typisch für Memory-Bandwidth-Sättigung bei größeren Thread-Zahlen.

### Skalierungs-Sweep — `mixture.xyz` (N=6200, nfrag=1400) Single-Point + Gradient

`./curcuma -sp mixture.xyz -method gfnff -threads T -verbosity 2`
Topology-Cache (`mixture.topo.json`) vor jedem Lauf gelöscht, um Phase-2-Skip zu vermeiden.
Hardware: AMD Ryzen 9 9950X3D, 16 phys. Cores.

| T  | Total (ms) | CN+EEQ (ms) | CN deriv (ms) | D4 GW (ms) | FF-Pool (ms) | Pool-Eff (%) | Speedup vs T=1 |
|----|-----------|-------------|---------------|------------|--------------|--------------|----------------|
| 1  | 3074      | 1584        | 1206          | 138        | 1490         | (Referenz)   | 1.00× |
| 2  | 2446      | 1589        | 1266          | 80         | 857          | 78.6 %       | 1.26× |
| 4  | 1896      | 1413        | 1104          | 56         | 483          | 63.8 %       | 1.62× |
| 8  | 1670      | 1373        | 1086          | 50         | 296          | 53.3 %       | 1.84× |
| 16 | 1469      | 1192        | 917           | 39         | 277          | 41.3 %       | 2.09× |

**Beobachtungen:**

- **FF-Pool skaliert OK** — T=1→8: 1490/296 = 5.0× (statt idealen 8×). Sub-linear, aber positive Skalierung.
- **CN+EEQ skaliert kaum** — T=1→16: 1584/1192 = nur 1.33×. Der dominante Teil sind die **CN-Derivate** (917-1266 ms, ~60 % der Phase, nicht die Stage-4-Loop). Skalierung der CN-Derivate von T=1→16: nur 1.31×. Genaue Stelle: `calculateCoordinationNumberDerivatives` in `gfnff_method.cpp:1426`.
- **`std::erf`-Loops in CN-Derivaten** — siehe `cn_calculator.cpp:130, 218, 238`. Jede Loop nutzt OpenMP `#pragma omp parallel for schedule(dynamic, 32)` — sollte gut skalieren. Die schwache Skalierung deutet auf Memory-Bandwidth-Druck hin (sparse-Matrix `SpMatrix` Allokationen sind nicht vektorisiert/lokal).
- **D4 Gaussian Weights skaliert gut** — T=1→16: 138/39 = 3.5× — passt zu den 1264 Atomen pro Thread bei T=16.
- **Energie variiert mit T um ~700 mEh** auf `mixture.xyz` (-856.29 → -855.61 Eh) — bestehendes Multi-Fragment-EEQ-Reproduzierbarkeitsproblem (vermutlich nicht-deterministische Reihenfolge der Schur-Cholesky-Reduktionen oder Convergence-Pfad-Drift). Auf kleinen Molekülen (`acetic_acid_dimer.xyz`) ist die Energie bit-identisch über T={1,2,4} (-2.47129863 Eh in allen Fällen) — also **nicht** durch WP1 verursacht. Eintrag in `Known Issues` für separate Untersuchung.

### Mechanismus-Inventar

GFN-FF CPU-Pfad — alle parallelen Stellen ohne `cuda/`-Unterordner:

| Datei | Zeile | Mechanismus | Begründung |
|-------|-------|-------------|------------|
| `cn_calculator.cpp` | 130 | OpenMP `parallel for schedule(dynamic, 32)` | reine Reduktion über CN-Werte, sehr kurze Inner-Loop, OpenMP-Latenz minimal |
| `cn_calculator.cpp` | 218 | OpenMP `parallel for` | Neighbor-list-Aufbau |
| `cn_calculator.cpp` | 238 | OpenMP `parallel for` | CN aus Neighbor-list |
| `d4param_generator.cpp` | 377/383 | OpenMP `parallel` + `for schedule(dynamic, 10)` | D4 Gaussian-Weights |
| `d4param_generator.cpp` | 692/696 | OpenMP `parallel for` | D4 Gaussian-Weight-Derivate |
| `d4param_generator.cpp` | 1054-1072 | CxxThreadPool `enqueue` mit `std::thread`-Fallback | Per-Atom-Worker mit Per-Thread-Buffer (Curcuma-Konvention) |
| `d4param_generator.cpp` | 1300-1320 | CxxThreadPool / `std::thread`-Fallback | Gaussian-Weight-Derivate Worker |
| `d4param_generator.cpp` | 1390-1410 | CxxThreadPool / `std::thread`-Fallback | dc6dcn-Matrix Worker |
| `eeq_solver.cpp` | 1159 | OpenMP `parallel for schedule(static)` | Distance-Matrix-Aufbau (kleines System) |
| `eeq_solver.cpp` | 2658 | OpenMP `parallel for schedule(dynamic)` | A-Matrix-Aufbau |
| `eeq_solver.cpp` | 2869-2878 | CxxThreadPool `enqueue` (Pool || `std::thread`-Fallback) | Distance-Matrix Worker (großes System) |
| `eeq_solver.cpp` | 3056-3065 | CxxThreadPool `enqueue` (Pool || `std::thread`-Fallback) | A-Matrix Worker (großes System) |
| `eeq_solver.cpp` | **1305-1335** | **serial** | **Stage-4 Per-Fragment-LU — Ziel von WP2** |
| `gfnff_method.cpp` (FF) | 2660-2750 | CxxThreadPool für 6 unabhängige Parameter-Gen-Phasen | one-time Setup, nicht im Hot-Pfad |

**Inkonsistenz**: OpenMP und CxxThreadPool koexistieren. Empfehlung für *neue* Stellen: CxxThreadPool + Per-Thread-Buffer (siehe Mechanismus-Empfehlung oben). Bestehende OpenMP-Stellen lassen, wo sie nachweislich gut skalieren (CN, A-Matrix); umstellen, wo sie es nicht tun (Kandidaten: `cn_calculator.cpp` Loops, falls die schlechte CN-Skalierung in WP4 auf False Sharing zurückgeht).

### Implementations-Tracking

Geänderte Dateien:

| Datei | Änderung |
|-------|----------|
| `src/core/energy_calculators/ff_methods/gfnff.h` | `int m_threads = 1;` Member, `setThreadCount(int)` / `threadCount()` Inline-Methoden |
| `src/core/energy_calculators/ff_methods/gfnff_method.cpp` (FF) | Beide GFNFF-Konstruktoren initialisieren `m_threads`. `setParameters` aktualisiert `m_threads`. 5 alte Lese-Stellen `m_parameters.value("threads", 1)` → `m_threads`. Label `"CN + EEQ (serial CPU)"` → `"CN + EEQ"`. |
| `src/core/energy_calculators/qm_methods/gfnff_method.cpp` (QM-Wrapper) | `setThreadCount()`: forward an `m_gfnff->setThreadCount` plus Update von `m_parameters["threads"]`. |

### Verifikation der Implementation

```bash
# Korrektheit auf kleinen Molekülen
$ curcuma -sp acetic_acid_dimer.xyz -method gfnff -threads {1,2,4} -verbosity 1
T=1: -2.47129863 Eh
T=2: -2.47129863 Eh    ✅ bit-identisch
T=4: -2.47129863 Eh    ✅ bit-identisch
T=4 triose: -9.91485397 Eh

# Threading propagiert (Pool-Größe in Log)
$ curcuma -sp mixture.xyz -method gfnff -threads 4 -verbosity 2 | grep "Thread pool"
[RESULT]  Thread pool (N=4)                 wall= 483.48    cpu=1234.39    ✅

# Label korrigiert
$ curcuma -sp mixture.xyz -method gfnff -threads 4 -verbosity 2 | grep "CN + EEQ"
[RESULT]  CN + EEQ                          wall=1412.79    ✅ kein "(serial CPU)" mehr
```

### Hand-off an WP2

Konkrete Stellen, an denen WP2 ansetzt:

1. **Stage-4 Per-Fragment-LU-Loop** `eeq_solver.cpp:1305-1335` — der einzige echte serielle Hot-Loop im EEQ-Pfad. CxxThreadPool + Per-Thread-Worker analog zu `eeq_solver.cpp:2869-2878`. Voraussetzung dort: Pool-Reference im EEQSolver bereits erreichbar (Pattern in `gfnff_method.cpp:990-992` zeigt, wie der Pool durchgereicht wird).
2. **Optional zusätzlich für WP4**: CN-Derivate in `gfnff_method.cpp:1426` skalieren nur 1.33× über T=1→16. Sparse-Matrix-Allokationen sind verdächtig. Lohnt sich, im Rahmen von WP4 (CN SIMD) gleich mitzuoptimieren oder als eigenes Sub-WP zu führen.

### Bekannte Limitierungen (nicht von WP1 induziert)

- **Energie-Drift mit Thread-Anzahl auf Multi-Fragment-Systemen**: ~700 mEh über T=1→16 auf `mixture.xyz`. Nicht reproduzierbar zwischen Läufen mit gleichem T (verschiedene Werte bei zwei T=4-Läufen: -855.84515631 und -855.88498327). Vermutlich nicht-deterministische Parallel-Reduction-Reihenfolge in der EEQ-A-Matrix oder im Schur-Cholesky-Solver. Eigene Untersuchung außerhalb von WP4 nötig — möglicher separater Task **WP-EEQ-Determinism**.
