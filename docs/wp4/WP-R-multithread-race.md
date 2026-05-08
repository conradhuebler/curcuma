# WP-R — Multi-Thread Race auf Single-Fragment-Energie (triose T=4)

**Status:** 🆕 Vorgeschlagen
**Aufwand:** 4–8h Diagnose + Fix
**Erwarteter Nutzen:** Reproducibility — vorerst kein Performance-Effekt, aber Voraussetzung für saubere Regression-Tests, MD-Stabilität und Optimizer-Konvergenz auf großen Systemen
**Voraussetzung:** keine

## Symptom

`triose.xyz` (22 Atome, **1 Fragment**, also kein Multi-Fragment-EEQ-Pfad) liefert bei `-threads 4` nicht-deterministische Energien zwischen aufeinanderfolgenden Läufen:

```
Run 1 T=4: -9.91485397 Eh   (häufigster Wert, matches T=1)
Run 2 T=4: -9.91485397 Eh
Run 3 T=4: -9.91337978 Eh   ← +1.47 mEh
Run 4 T=4: -9.91485397 Eh
Run 5 T=4: -9.91463123 Eh   ← +0.22 mEh
```

Drei diskrete stabile Endzustände, **kein** kontinuierlicher Round-Off-Drift. Bei `-threads 1` immer `-9.91485397 Eh` (deterministisch).

**Verifiziert pre-WP-G:** 5 Läufe vor der WP-G-Implementation zeigen exakt dieselbe Verteilung. Die Race ist **pre-existing**, nicht durch WP-G eingeführt. Auch nicht durch WP4 oder andere kürzliche WPs — das ist ein bestehender latenter Bug.

**Verifiziert bei kleineren Systemen:** `acetic_acid_dimer.xyz` (16 Atome, 2 Fragmente) ist bei T={1,4} bit-identisch (`-2.47129863 Eh`). Die Race triggert also erst ab einer bestimmten Komplexität.

## Abgrenzung zu anderem T-Drift

`mixture.xyz` (N=6200, nfrag=1400) zeigt T-Drift in der Größenordnung 700 mEh. Das ist ein **anderes Phänomen** — Multi-Fragment-EEQ-Reproduzierbarkeit, dokumentiert in WP1. Der EEQ-Solver hat parallel-non-deterministische Reduktion bei stark fragmentierten Systemen.

WP-R fokussiert auf **Single-Fragment**-Drift, der diesen Pfad nicht erklärt.

## Hypothesen

Drei stabile Endzustände bei kleinen mEh-Drifts deuten auf **algorithmische Inkonsistenz** hin, die je nach Thread-Scheduling zu leicht unterschiedlichen Konvergenz-Pfaden in einem iterativen Solver führt. Kandidaten:

### H1: Iterativer EEQ-Solver (PCG/SchurCholesky)

Die EEQ-Charge-Berechnung in `eeq_solver.cpp` hat einen iterativen Refinement-Pfad (Phase 2). Wenn die Konvergenz threading-empfindlich ist (z. B. unterschiedliche Reduktion-Reihenfolge in der A-Matrix-Konstruktion → leicht unterschiedliche Eigenwerte → leicht unterschiedliche LU/Cholesky-Pivot-Wahl), kann der finale Charges-Vektor minimal anders sein.

**Test:** Logge `m_eeq_charges` nach Phase 2 für jeden T=4-Lauf. Wenn die Charges driften, ist EEQ die Quelle.

### H2: Hückel-Solver (Pi-bond-orders)

`HuckelSolver` (Jan 2026) macht einen iterativen SCF auf einer Sub-Hamiltonian-Matrix. SCF-Konvergenz kann threading-empfindlich sein (Eigenvektor-Pivoting, Linear-Algebra-Pivot-Wahl).

**Test:** `huckel_solver.cpp` Konvergenz-Iterations-Count und finalen Eigenwert pro Run loggen. Falls Iterations-Count zwischen Runs schwankt → SCF konvergiert zu anderen Lösungen.

### H3: HB/XB-Detection ist atom-order-abhängig

`computeHBCoordinationNumbers` und das HB/XB-Setup könnten eine atom-order-abhängige Selektion machen (z. B. "nimm den nächsten H-Donor" → bei equal distance wählt der Code den ersten in der Iteration). Bei T=4 wird die Reihenfolge der Atom-Iteration durch Pool-Scheduling beeinflusst.

**Test:** `m_hb_pairs.size()` und Liste der HB-Atom-Tripel pro Run vergleichen.

### H4: Pool-Scheduling der Parameter-Generierung

`generateGFNFFParameterSet` läuft selbst parallel (T>1) und nutzt `CxxThreadPool::enqueue` für mehrere Sub-Phasen (Bonds, Angles, Coulomb, Dispersion etc., siehe `gfnff_method.cpp:2660`). Wenn diese Sub-Phasen sich auf gemeinsamen Member-Variablen schreiben (oder ein Sub-Phase auf das Output einer anderen wartet), kann die Reihenfolge variieren.

**Test:** `n_threads_for_param_gen` auf 1 setzen, T=4 nur für den Energie-Pfad → wenn deterministisch, ist Param-Gen die Quelle.

### H5: Topology-Caching (`mixture.topo.json`-Pfad)

Der erste Lauf nach Topo-Cache-Reset baut die Topologie neu. Beim zweiten Lauf wird sie aus dem Cache geladen und Phase-2-EEQ wird ggf. übersprungen (`skip_phase2`). Bei triose ohne Topo-Cache-Reset könnten Cache-vs-Neu-Berechnung-Pfade subtle unterschiedliche Charges erzeugen.

**Test:** Mit explizitem Topo-Cache-Reset zwischen jedem Lauf (das wurde aber schon im Test verifiziert: Drift trotzdem). Eher nicht die Quelle.

### H6: Eigen-internal Threading

`Eigen::nbThreads()` ist nicht explizit auf 1 gesetzt. Wenn Eigen intern threaded und mehrere FF-Methoden gleichzeitig Eigen-Operationen ausführen (über CxxThreadPool-Worker), kann das Eigen-Thread-Pool selbst Race-anfällig sein.

**Test:** `Eigen::setNbThreads(1)` global → wenn deterministisch, ist Eigen die Quelle.

## Diagnose-Schritte (nach Aufwand sortiert)

### Phase 1 — Quick Wins (1h)

1. **Eigen-Threading off:** `Eigen::setNbThreads(1)` global vor erstem FF-Call. Run 5× triose T=4. Falls deterministisch → H6 bestätigt.

2. **Topo-Cache eliminieren:** `cache_topology=false` für triose. Run 5× T=4. Falls Drift weg → H5 bestätigt (unwahrscheinlich, schon getestet).

3. **Param-Gen seriell:** `generateGFNFFParameterSet`-Pfad mit `thread_count=1` aber Energy-Pfad mit T=4. Run 5× T=4. Falls deterministisch → H4 bestätigt.

### Phase 2 — Component-Logging (2-3h)

4. **EEQ-Charges loggen:** Nach `prepareCNAndEEQ` finalen `m_eeq_charges`-Vektor per Atom drucken. 5 Läufe T=4 → Diff. Wenn Charges driften → H1.

5. **Hückel-Iteration-Count:** SCF-Iterationen pro Run loggen. Falls verschiedene Counts → H2.

6. **HB/XB-Listen:** Größe und Inhalt der HB-Pair-Liste pro Run drucken. Falls verschiedene → H3.

### Phase 3 — Fix (1-4h, abhängig von Hypothese)

Abhängig vom Befund:

- **H1 (EEQ):** Reduktion-Reihenfolge im A-Matrix-Aufbau deterministisch machen — z. B. immer `i+=T`-Stride mit fester Thread-Range, oder `Eigen::setNbThreads(1)` lokal um Eigen-Solver-Aufrufe herum.
- **H2 (Hückel):** SCF konvergenz-Toleranz erhöhen oder Initial-Guess deterministisch machen.
- **H3 (HB/XB):** Selektions-Tie-Breaker explizit machen (`if (dist == best_dist && atom_id < best_atom)`).
- **H4 (Param-Gen):** Param-Generation auf seriellen Pfad fixieren oder Reduktions-Reihenfolge garantieren.
- **H5/H6:** Konfiguration anpassen.

## Kritische Dateien

| Datei | Vermuteter Bezug |
|-------|-------------------|
| `src/core/energy_calculators/ff_methods/eeq_solver.cpp` | H1 — A-Matrix-Aufbau parallel, Schur/Cholesky |
| `src/core/energy_calculators/ff_methods/huckel_solver.cpp` | H2 — SCF |
| `src/core/energy_calculators/ff_methods/gfnff_method.cpp` | H3 — `computeHBCoordinationNumbers`, HB/XB-Detection |
| `src/core/energy_calculators/ff_methods/gfnff_method.cpp:2660` | H4 — `generateGFNFFParameterSet` parallel sub-phases |

**Bestehende Helper:**
- WP1-Skalierungs-Sweep-Methodik mit `verbosity 2` Phase-Logging
- `CurcumaLogger::result` für deterministisches Output-Logging

## Verification

```bash
cd release && make -j4

# Reproduzierbarkeitstest — 10 Läufe T=4 müssen identische Energie liefern
for i in $(seq 1 10); do
  rm -f ../test_cases/molecules/larger/triose.topo.json
  ./curcuma -sp ../test_cases/molecules/larger/triose.xyz -method gfnff -threads 4 -verbosity 1 \
    | grep "Final Energy"
done | sort -u | wc -l
# Erwartung: 1 (also alle Läufe identisch)

# Skalierungs-Konsistenz: T=1 vs T=4 vs T=8 — alle bit-identisch
for T in 1 2 4 8 16; do
  rm -f ../test_cases/molecules/larger/triose.topo.json
  ./curcuma -sp ../test_cases/molecules/larger/triose.xyz -method gfnff -threads $T -verbosity 1 \
    | grep "Final Energy"
done | sort -u | wc -l
# Erwartung: 1

# Cross-System-Test
for sys in acetic_acid_dimer triose CH3OCH3; do
  for T in 1 4; do
    rm -f ../test_cases/molecules/larger/${sys}.topo.json
    echo -n "$sys T=$T: "
    ./curcuma -sp ../test_cases/molecules/larger/${sys}.xyz -method gfnff -threads $T -verbosity 1 \
      | grep "Final Energy"
  done
done
# Erwartung: pro System sind T=1 und T=4 bit-identisch
```

## Akzeptanzkriterien

1. ⚙️ `triose.xyz` 10 Läufe T=4 mit Topo-Cache-Reset → genau **1 unique Energie-Wert**.
2. ⚙️ T=1 vs T=4 vs T=8 vs T=16 für triose, CH3OCH3, caffeine → bit-identisch.
3. ⚙️ `mixture.xyz` Multi-Fragment-Drift bleibt unverändert (das ist ein **anderes** Problem, nicht WP-R-Scope).
4. ⚙️ Performance-Regression < 5 % auf `mixture.xyz` T=4 (falls Fix Eigen-Threading deaktiviert oder serielle Reduktion erzwingt).
5. ⚙️ Eintrag in Roadmap-Status für P4i auf `⚙️ Machine-tested` mit Hypothesen-Befund + Fix.
6. Operator setzt ✅ TESTED.

## Risks

- **Mittel.** Race-Condition-Bugs sind notorisch schwer zu reproduzieren — die Diagnose-Phase könnte länger dauern als geschätzt.
- **Performance-Trade-off:** Falls H1 oder H6 bestätigt wird, muss möglicherweise serieller Code-Pfad eingeführt werden, der Performance kostet. Akzeptabel, weil Reproduzierbarkeit > Speedup.
- **Nicht alle Drift-Quellen sind reine Bugs:** Floating-Point-Reduktion-Reihenfolge ist mathematisch nicht-assoziativ. Eine geringfügige Drift im µEh-Bereich ist akzeptabel — aber 1.47 mEh ist eindeutig zu groß.

## Out of Scope

- `mixture.xyz` Multi-Fragment-EEQ-Drift mit T (700 mEh) — eigenes Phänomen, eigenes WP wert.
- GPU-Reproduzierbarkeit — separate Testumgebung nötig.
- MD-Long-Run-Stabilität (1000+ Schritte) — Folgewirkung der Race auf Trajektorien-Drift, separat zu untersuchen.

## Beziehung zu anderen WPs

- **WP1 dokumentierte den `mixture.xyz` Multi-Fragment-Drift** als bestehende Limitierung — WP-R ist verwandt, fokussiert aber auf Single-Fragment.
- **WP4 + WP-G:** beide haben in ihrer Doku den triose T=4-Drift als pre-existing notiert — WP-R ist die Folgeanalyse.
- **Vorbedingung für künftige Performance-WPs:** ohne deterministische Energien sind Speedup-Vergleiche unsicher (z. B. WP-G zeigt 3-7% Total-Speedup, davon ist ein Teil vermutlich Run-Varianz statt echter Gewinn).
