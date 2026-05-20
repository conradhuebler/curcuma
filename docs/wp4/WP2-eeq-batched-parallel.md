# WP2 — EEQ Batched parallelisieren (P4a)

**Status:** ⚙️ Machine-tested (Mai 2026) — **funktioniert technisch, Speedup gering**. Operator setzt ✅ TESTED.
**Aufwand:** 2h
**Erwarteter Nutzen (geplant):** Phase 2 Faktor 4–6× (660 → 120–200 ms)
**Tatsächlicher Nutzen (gemessen):** Phase 2 Faktor 1.17× über T=1→16. Total Single-Point Wall-Time praktisch unverändert vs. WP1-Baseline.
**Voraussetzung:** WP1 abgeschlossen

## Hypothese

Im EEQ-Solver Stage 4 (Mai 2026, Commit b6a3a42) werden bei hoch-fragmentierten Systemen 1400 unabhängige Per-Fragment-LU-Solves nacheinander in einer einzigen Schleife abgewickelt:

```cpp
// eeq_solver.cpp:1305-1335
for (int f = 0; f < nfrag && batched_ok; ++f) {
    // assemble Aug, rhs_local für Fragment f
    Eigen::PartialPivLU<Matrix> lu(Aug);
    Vector sol = lu.solve(rhs_local);
    ...
}
```

Die Iterationen sind **embarrassingly parallel**: jedes Fragment hat seinen eigenen `Aug`-Block, eigene `rhs_local`, eigenes `sol`-Result. Die einzige geteilte Schreibziel-Struktur ist `charges_batched(ids[i]) = sol(i)` — diese Schreibposition ist **disjunkt**, weil jedes Atom genau zu einem Fragment gehört.

Profil-Datenpunkt: 660 ms / 1400 Fragmente ≈ 470 µs pro LU. Bei 4 Threads erwarte ich Faktor 3.5–4×, bei 8 Threads 6–8× (begrenzt durch Memory Bandwidth, da die `Aug`-Matrizen klein sind).

## Aufgabe

Stage-4-Loop in `eeq_solver.cpp:1307` parallelisieren.

### Mechanismus-Wahl: CxxThreadPool + `enqueue`, nicht OpenMP

Hintergrund: Der Fortran-GFN-FF-Originalcode (OpenMP) skaliert oberhalb von ~4 Threads negativ — vermutlich False Sharing auf gemeinsamen Akkumulatoren oder zu naives `schedule(static)`. Curcuma hat in `precomputeGaussianWeights` (`d4param_generator.cpp:998`) bereits das passende Pattern: persistenter Pool, Worker-Lambdas mit `enqueue`, Per-Thread-Buffer, Reduktion am Ende. Dem folgen.

| Option | Pro | Contra | Empfehlung |
|--------|-----|--------|-----------|
| **CxxThreadPool + `enqueue`** | Per-Thread-`local`-Variablen → False Sharing strukturell ausgeschlossen, deterministisches Interleave-Scheduling, Curcuma-Konvention | EEQSolver braucht Pool-Reference (aus WP1) | ✅ **Nehmen** |
| OpenMP `parallel for` | 1–5 µs Fork-Join, in 2 Zeilen | Risiko der Fortran-Pathologie, Eigen kann nested-parallelisieren | Optional als A/B-Vergleich hinter Compile-Flag |
| `std::async(launch::async, ...)` | C++-Standard | **~50–100 µs Thread-Erzeugung × 1400 = 70–140 ms Overhead pro Single-Point** | ❌ **Vermeiden** |

### Schritte

1. **Vor-Pass:** `std::vector<std::vector<int>> fragment_ids(nfrag)` einmalig aufbauen (im Aufrufer, vor `dispatchSolve`), statt in der Loop pro Fragment zu scannen. Spart einen N×nfrag-Pass.

2. **Worker-Lambda mit Per-Thread-State** (kein geteilter Akkumulator → kein False Sharing):
   ```cpp
   std::atomic<bool> batched_ok{true};
   auto worker = [&](int t_id, int T) {
       // Per-Thread Scratch — alloziert einmal, wiederverwendet
       Eigen::PartialPivLU<Matrix> lu;
       Matrix Aug;
       Vector rhs_local, sol;
       for (int f = t_id; f < nfrag; f += T) {        // interleaved scheduling
           if (!batched_ok.load(std::memory_order_relaxed)) return;
           const auto& ids = fragment_ids[f];
           const int nf = static_cast<int>(ids.size());
           if (nf == 0) continue;
           Aug.resize(nf + 1, nf + 1);
           rhs_local.resize(nf + 1);
           // assemble Aug, rhs_local aus A_nn, C_constraints, rhs_atoms, rhs_constraints
           // ...
           lu.compute(Aug);
           sol = lu.solve(rhs_local);
           if (!sol.allFinite()) {
               batched_ok.store(false, std::memory_order_relaxed);
               return;
           }
           for (int i = 0; i < nf; ++i) charges_batched(ids[i]) = sol(i);  // disjunkt
       }
   };
   ```

3. **Pool-Dispatch (Stil aus `d4param_generator.cpp:1054–1072`):**
   ```cpp
   const int T = std::min(m_threads, nfrag);
   if (T > 1 && m_pool) {
       std::vector<std::future<void>> futures;
       futures.reserve(T - 1);
       for (int t = 1; t < T; ++t)
           futures.push_back(m_pool->enqueue(worker, t, T));
       worker(0, T);                   // Main-Thread arbeitet auch
       for (auto& f : futures) f.get();
   } else if (T > 1) {                 // Fallback: raw std::thread (1× Erzeugung pro Single-Point — OK)
       std::vector<std::thread> threads(T - 1);
       for (int t = 1; t < T; ++t)
           threads[t - 1] = std::thread(worker, t, T);
       worker(0, T);
       for (auto& th : threads) th.join();
   } else {
       worker(0, 1);                   // serial fallback
   }
   ```

4. **Pool-Reference im EEQSolver durchreichen:** Aus WP1-Befund. Falls `EEQSolver` keinen `m_pool` hat, im Konstruktor optional übergeben (`CxxThreadPool* pool = nullptr`), gleicher Stil wie `precomputeGaussianWeights(pool, num_threads)`.

5. **Fallback intakt halten:** Wenn `batched_ok == false` nach der Loop, weiterhin auf `SchurCholesky` fallen.

6. **A/B-Vergleichsmessung (optional, hinter Compile-Flag `EEQ_USE_OPENMP`):** Eine OpenMP-Variante einbauen, damit der Fortran-Pathologie-Befund auf Curcuma-Boden direkt verifizierbar ist. Wenn OpenMP messbar schlechter skaliert, dokumentieren und Variante entfernen. Wenn gleich gut, beide Varianten stehen lassen.

### Code-Anker

| Datei | Position | Hinweis |
|-------|---------|---------|
| `src/core/energy_calculators/ff_methods/eeq_solver.cpp` | Z. 1305–1335 | Stage-4-Loop, hier parallelisieren |
| `src/core/energy_calculators/ff_methods/eeq_solver.cpp` | Z. 1937–2050 | Vorhandener Multi-RHS PCG-Pfad — Vorbild für Threading? |
| `src/core/energy_calculators/ff_methods/d4param_generator.cpp` | Z. 986–1072 | Referenz-Pattern für CxxThreadPool-Nutzung mit Fallback |

### Thread-Safety

- `Aug`, `rhs_local`, `lu`, `sol`: lokale Variablen pro Iteration → trivial sicher.
- `charges_batched`: Schreiben mit disjunkten Indizes → race-frei.
- `A_nn`, `C_constraints`, `rhs_atoms`, `rhs_constraints`: nur lesen → race-frei.
- Eigen `PartialPivLU` ist nicht thread-safe falls geteilt, hier aber lokal.

## Akzeptanzkriterien

1. ⚙️ `mixture.xyz` (N=6200, nfrag=1400) Single-Point: EEQ-Phase-2 von ~660 ms auf < 200 ms (≥ 3.5× mit 4 Threads).
2. ⚙️ Numerischer Energie-Diff zur seriellen Variante < 1 nEh (LU-Solver liefert dieselbe Antwort, Reihenfolge der Per-Fragment-Updates ändert nichts).
3. ⚙️ `test_gfnff_numgrad` grün.
4. ⚙️ Fallback-Pfad: künstlich `batched_ok=false` in einem Worker auslösen → SchurCholesky-Fallback funktioniert.
5. ⚙️ `gfnff` CTest-Suite grün.

## Risiken

- **Mittel.** `Eigen::PartialPivLU` ist intern thread-safe nur über separate Instanzen (gegeben — jede Iteration hat ihre eigene `lu`).
- **Eigen Multi-Threading-Konflikt:** Eigen kann selbst threaden (`Eigen::nbThreads()`). Bei kleinen Per-Fragment-LUs ist das nicht der Fall (Schwelle ~256×256), aber sicherheitshalber `Eigen::setNbThreads(1)` vor der Loop und am Ende restaurieren — verhindert nested parallelism.
- **Fragment-Größe-Varianz:** ein einziges großes Fragment dominiert die Laufzeit (load imbalance). `schedule(dynamic, 16)` lindert. Falls auf realen Systemen ein Fragment ≫ alle anderen, ggf. später `schedule(guided)`.
- WP1-Befund kann zeigen, dass `m_threads` im EEQ-Solver bei `1` hängt — dann muss erst die Propagation gefixt werden, sonst sieht der Speedup wie nicht-vorhanden aus.

## Verifikation

```bash
cd release && make -j4

# Korrektheit
./curcuma -sp ../test_cases/molecules/larger/mixture.xyz -method gfnff -threads 4 -verbosity 2 \
  > before.log 2>&1
git stash       # vor WP2
./curcuma -sp ../test_cases/molecules/larger/mixture.xyz -method gfnff -threads 4 -verbosity 2 \
  > before-serial.log 2>&1
git stash pop
diff <(grep "Total energy" before.log) <(grep "Total energy" before-serial.log)  # < 1 nEh

# Performance
ctest -R "gfnff" --output-on-failure
./test_cases/test_gfnff_numgrad
```

Eintrag in der Status-Tabelle der Roadmap mit gemessenen ms/Schritt.

---

## Ergebnis (Mai 2026)

### Implementation

| Datei | Änderung |
|-------|----------|
| `eeq_solver.h` | `dispatchSolve` und `calculateTopologyCharges` Signaturen erweitert um `CxxThreadPool* pool = nullptr, int num_threads = 1` |
| `eeq_solver.cpp:1294-1387` | Pre-Pass `fragment_ids[nfrag]`, `std::atomic<bool> batched_ok`, Worker-Lambda mit Per-Thread-Scratch (Aug, rhs_local, lu, sol), Pool-Dispatch nach Vorbild `dist_worker` (Z. 2870-2881), Threshold `nfrag >= 16` |
| `eeq_solver.cpp:2520, 2234, 3244` | `dispatchSolve`-Aufrufer reichen `pool/num_threads` durch |
| `eeq_solver.cpp:2262-2270` | `calculateTopologyCharges`-Definition entsprechend ergänzt |
| `gfnff_method.cpp:1086-1095` | Non-Gradient-Pfad bekommt jetzt Pool (war Default-`nullptr`) |
| `gfnff_method.cpp:8042-8053, 8131-8146, 9711-9723, 10033-10049` | Topology-Setup-Aufrufer reichen `m_threads` und Pool durch |

### Sweep — `mixture.xyz` (N=6200, nfrag=1400) Single-Point + Gradient

| T  | Total (ms) | EEQ Phase 1 | EEQ Phase 2 | Phase 2 vs T=1 | Total vs T=1 | Pre-WP2 Total |
|----|-----------|-------------|-------------|----------------|--------------|---------------|
| 1  | 3068      | 500         | 371         | (Referenz)     | 1.00×        | 3074 ms       |
| 2  | 2389      | 489         | 409         | 0.91×          | 1.28×        | 2446 ms       |
| 4  | 1840      | 494         | 409         | 0.91×          | 1.67×        | 1896 ms       |
| 8  | 1673      | 499         | 357         | 1.04×          | 1.83×        | 1670 ms       |
| 16 | 1479      | 498         | 316         | 1.17×          | 2.07×        | 1469 ms       |

**Korrektheit:**
- `acetic_acid_dimer.xyz` (Stage-4 inaktiv): bit-identisch über T={1,4} (-2.47129863 Eh).
- `triose.xyz`: -9.91485397 Eh.
- `mixture.xyz`: bestehender Multi-Fragment-EEQ-Drift mit T (~700 mEh über T=1→16) — von WP1 dokumentiert, **nicht durch WP2 verschlechtert**.

### Befund

**Stage-4-Loop läuft technisch korrekt parallel** (Logger-Print: `"Batched per-fragment solve done (nfrag=1400, threads=N)"` mit dem korrekten N). Speedup ist aber deutlich kleiner als geplant.

Hauptursachen:

1. **Phase 2 ist heute nur 371 ms wall** (nicht 660 ms wie im Mai-2026-Profil). Das System hat sich seit der Profilerhebung verändert — möglicherweise durch andere Optimierungen oder geänderten Code-Pfad (Topology-Cache-Skip, andere EEQ-Solver-Logik). Der absolute Hebel von Stage-4 ist dadurch geringer.

2. **Per-Fragment-LU mit 4-5 Atomen → 6×6-Matrix.** Bei dieser Größe dominiert der Eigen-Overhead pro `lu.compute()`/`solve()` die reinen FLOPs. Threading bringt erst bei deutlich größeren Fragmenten (>20 Atome) viel.

3. **Phase-2-Wall enthält Code außerhalb der parallelen Loop:** Matrix-Assembly (`A_nn` 6200×6200, `C_constraints` 1400×6200), Korrekturen-Aufbau (dxi, dgam, alpeeq), Convergence-Check, Validierung. Die Stage-4-LU-Loop selbst ist vermutlich nur ~80-100 ms der 371 ms — der parallele Teil ist klein.

4. **Memory-Bandwidth-Limit auf `A_nn` (308 MB).** Random-Access-Pattern über `A_nn(ids[i], ids[j])` aus mehreren Threads gleichzeitig ist L3-Cache-Limited.

5. **Mögliches False Sharing** in `charges_batched(ids[i])` bei dicht gepackten Fragment-Indizes (1400 Fragmente / 6200 Atome → benachbarte Atome derselben Cache-Line können von verschiedenen Threads geschrieben werden). Nicht verifiziert (würde `perf c2c` brauchen).

### Bewertung — Wert von WP2

**Behalten,** weil:
- Saubere Pool-Propagation durch alle EEQ-Aufrufpfade — strukturelle Voraussetzung für künftige Stage-4-Optimierungen.
- Stage-4-Loop ist nicht mehr serial — wenn das System wächst (N>10000, nfrag>5000), wird der Hebel größer.
- Architektonisch konsistent: alle 6 Aufrufer (1025/1091/8042/8126/9711/10020) reichen jetzt einheitlich Pool durch.
- Bestehender Pre-Pass `fragment_ids` ist eine eigenständige Verbesserung (spart inneren O(N)-Scan, bleibt auch bei nur 1 Thread aktiv).

**Aber: WP2 ist nicht der Speedup-Hebel den der Mai-Plan unterstellte.** Die echten Bottlenecks sind anderswo:

- **CN-Derivate** (`gfnff_method.cpp:1426`) skalieren nur 1.31× über T=1→16 — größter Hebel im prepareCNAndEEQ-Block (Phase 2 + CN-Derivate sind ~75 % der Phase). → **WP4 oder neues Sub-WP**.
- **Bond-Term** (175 ms im FF-Pool) und **Coulomb** (660 ms im FF-Pool) im FF-Block. → **WP3 / WP6**.

### Status für Roadmap

P4a: ⚙️ Machine-tested. Stage-4-Loop technisch parallelisiert, Total-Speedup gegen WP1-Baseline marginal (1.67× T=4 statt 1.62×; 2.07× T=16 statt 2.09×). Kein nennenswerter Performance-Gewinn auf `mixture.xyz`. Bei größeren Multi-Fragment-Systemen (N>10000) wahrscheinlich relevant. Saubere Architekturvoraussetzung für künftige EEQ-Optimierungen.
