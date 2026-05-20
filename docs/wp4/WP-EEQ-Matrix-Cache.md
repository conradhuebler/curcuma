# WP-EEQ-Matrix-Cache — Zweite Cache-Stufe für die EEQ-A_nn-Off-Diagonale

**Status:** 🆕 Vorgeschlagen (Mai 2026)
**Aufwand:** ~1-2 Tage
**Erwarteter Nutzen:** Polymer N=1410: **0.2-0.5 s / 100 MD-Steps** zusätzlich zu WP-EEQ-Cache (Cholesky-Cache).
**Voraussetzung:** WP-EEQ-Cache (gefixt) bereits aktiv.

## Hypothese

WP-EEQ-Cache (Mai 2026, gefixt) cached die **Cholesky-Faktorisierung** von A_nn über MD-Steps. Sehr gut für SchurCholesky-Solver. Aber für polymer N=1410 wird der **Schur-Cholesky-Pfad mit eps=0.05 nur in 25% der Calls reused** — die anderen 75% bauen die A_nn-Matrix neu UND faktorisieren neu.

Das **Neu-Bauen** der A_nn besteht aus:
1. Off-Diagonale: `A_ij = erf(γ_ij · r_ij) / r_ij` — O(N²) erf-Aufrufe pro Step → ~2-5 ms für polymer
2. Diagonale: `A_ii = gam_corrected(i) + 2/sqrt(π) · α_corrected(i) + ∑_j A_ij` — O(N) pro Step
3. Verschiedene Korrekturen (dxi, dgam, alpeeq) — O(N), schon gecacht

Wenn die Geometrie pro Step nur leicht ändert, sind die meisten `A_ij`-Werte praktisch konstant — sie ändern sich nur, wenn die Distanz `r_ij` sich relevant ändert. Bei Polymer mit dt=1 fs ändert sich vielleicht 5-10% der Paar-Distanzen um >0.05 Bohr. Die restlichen 90-95% der A_ij-Einträge könnten **inkrementell** geupdated werden statt vollständig neu gebaut.

Wenn der LLT-Cache nicht greift (Hit-Rate 25%), ist das primär die Cholesky-Faktorisierung (22 ms). Aber der **Matrix-Aufbau** vor der Faktorisierung kostet weitere 2-5 ms — und der ist immer mit O(N²) erf-Calls verbunden.

## Aufgabe

### 1. Inkrementeller A_nn-Off-Diagonal-Cache

```cpp
// eeq_solver.h
struct EEQMatrixCache {
    Matrix A_offdiag_last;     ///< Last computed A_ij off-diagonal (NxN)
    Matrix last_geometry;      ///< Geometry at last full rebuild
    Vector last_cn;
    double rebuild_eps_bohr = 0.20;  ///< if max_dr > this, full rebuild
    double update_eps_bohr  = 0.02;  ///< if pair-r changes by > this, update that pair only
    bool valid = false;
    void reset() { valid = false; }
};
mutable EEQMatrixCache m_matrix_cache;
```

In `calculateFinalCharges` vor dem Coulomb-Matrix-Build:

```cpp
const double max_dr = (geometry_bohr - m_matrix_cache.last_geometry).rowwise().norm().maxCoeff();
const double max_dcn = (cn - m_matrix_cache.last_cn).cwiseAbs().maxCoeff();

if (!m_matrix_cache.valid || max_dr > m_matrix_cache.rebuild_eps_bohr || max_dcn > 0.05) {
    // Full O(N²) rebuild (current behavior)
    rebuildAnnOffDiagonal(geometry_bohr, cn, dgam, A);
    m_matrix_cache.A_offdiag_last = ...;  // store result
    m_matrix_cache.last_geometry  = geometry_bohr;
    m_matrix_cache.last_cn        = cn;
    m_matrix_cache.valid          = true;
} else {
    // Incremental update: only recompute A_ij for atom-pairs where |Δr_ij| > update_eps
    // Atoms that moved more than update_eps_bohr — find them, mark all their pairs dirty
    std::vector<int> dirty_atoms;
    for (int i = 0; i < natoms; ++i) {
        double dr = (geometry_bohr.row(i) - m_matrix_cache.last_geometry.row(i)).norm();
        if (dr > m_matrix_cache.update_eps_bohr) dirty_atoms.push_back(i);
    }
    // Re-compute only A_ij where i OR j is in dirty_atoms
    A = m_matrix_cache.A_offdiag_last;  // start from last
    for (int i : dirty_atoms) {
        for (int j = 0; j < natoms; ++j) {
            if (i == j) continue;
            double r = (geometry_bohr.row(i) - geometry_bohr.row(j)).norm();
            double gamma_ij = ...;
            A(i, j) = std::erf(gamma_ij * r) / r;
            A(j, i) = A(i, j);  // symmetric
        }
    }
    // Update cache for next step (do NOT update last_geometry — only on full rebuild)
}
```

**Wichtig**: `last_geometry` und `A_offdiag_last` werden **nur beim Full-Rebuild** aktualisiert, nicht beim inkrementellen Update. Sonst akkumulieren Fehler.

### 2. Validierung gegen Pre-Cache

Long-Run-MD (1000 Steps) vergleichen mit Cache disabled. Energie-Drift sollte < 1 mEh sein (Konsistenz mit WP-EEQ-Cache-Risiko-Diskussion).

### 3. PARAM-Exposition

```cpp
// eeq_solver.h PARAMETER_DEFINITION
PARAM(eeq_matrix_rebuild_eps_bohr, Double, 0.20,
      "WP-EEQ-Matrix-Cache: max atom displacement before full O(N²) rebuild of EEQ "
      "Coulomb matrix. 0.0 = always rebuild. Default 0.20 covers typical MD steps.",
      "Algorithm", {})
PARAM(eeq_matrix_update_eps_bohr, Double, 0.02,
      "WP-EEQ-Matrix-Cache: per-atom displacement threshold for incremental matrix update. "
      "Pairs involving atoms moved by less than this keep cached A_ij values.",
      "Algorithm", {})
```

### 4. Threading

Der inkrementelle Update-Loop ist embarrassingly parallel über `dirty_atoms`. CxxThreadPool-Pattern analog `dist_worker` in eeq_solver.cpp:2990.

### Code-Anker

| Datei | Position | Aufgabe |
|-------|----------|---------|
| `src/core/energy_calculators/ff_methods/eeq_solver.h` | nach `EEQCholeskyCache` (~L911) | `EEQMatrixCache` struct |
| `src/core/energy_calculators/ff_methods/eeq_solver.cpp` | calculateFinalCharges, vor Coulomb-Build (~L3030+) | Cache-Check + Inkrementelles Update |
| `src/core/energy_calculators/ff_methods/eeq_solver.h` | PARAMETER_DEFINITION block (~L1014) | Neue PARAMs |
| `src/core/energy_calculators/ff_methods/gfnff_method.cpp` | `getCachedTopology` (~L879) | `m_eeq_solver->invalidateMatrixCache()` zusätzlich zu Cholesky |

## Akzeptanzkriterien

1. ⚙️ Polymer 100-Step MD: Wall-Time **≤ 9.4 s** (Cholesky-Cache aktiv, Matrix-Cache aktiv).
2. ⚙️ `eeq_matrix_rebuild_eps_bohr = 0`: bit-identisch zu pre-WP (alle Pfade Full-Rebuild).
3. ⚙️ Energie-Drift über 1000-Step-MD: < 1 mEh.
4. ⚙️ Single-Point caffeine/triose: Energie bit-identisch (erster Call ist immer Full-Rebuild).
5. ⚙️ `ctest -R cli_simplemd` 7/7 pass.

## Risiken

1. **Hellmann-Feynman-Konsistenz** — wenn A_ij für ein Paar stale ist, aber das Paar gerade dann gradient-relevant ist, gibt es einen Gradient-Fehler. Lösung: `update_eps_bohr` deutlich kleiner als `refactor_eps_bohr` setzen (0.02 vs 0.05).
2. **Speicher** — `A_offdiag_last` ist N² × 8 B = 15 MB für N=1410. Akzeptabel.
3. **Komplexität** — der gesamte Cache-Stack (Cholesky + Matrix + Pending) wird unübersichtlich. Educational-First-Constraint sagt: dokumentieren, einfache Helper-Strukturen, keine Templates.
4. **dirty_atoms-Ineffizienz** — falls 50% der Atome dirty sind, ist der inkrementelle Pfad genauso teuer wie Full-Rebuild. Threshold gut wählen, ggf. Fallback zu Full-Rebuild wenn |dirty| > 0.3·N.

## Verifikation

```bash
cd release && make -j4

# Korrektheit erst
./curcuma -sp test_cases/molecules/larger/caffeine.xyz -method gfnff -v 2
# Energie sollte bit-identisch sein

# Long-Run-Drift
/bin/rm -f test_cases/molecules/larger/polymer.topo.json
./release/curcuma -md test_cases/molecules/larger/polymer.xyz -method gfnff \
   -print 1 -maxtime 1000 -threads 1 -md.md_diagnostics true > md1000.log
# Energy column drift over 1000 fs < 1 mEh

# Performance
time ./release/curcuma -md test_cases/molecules/larger/polymer.xyz -method gfnff \
   -print 1 -maxtime 1e2 -threads 1
# Erwartet: ≤ 9.4 s

# Falsifizierungstest (Cache off)
./release/curcuma -md polymer.xyz -method gfnff -maxtime 1e2 -threads 1 \
   -eeq_solver.eeq_matrix_rebuild_eps_bohr 0
# Sollte identisch zu pre-WP-Verhalten sein
```

## Beziehung zu anderen WPs

- **WP-EEQ-Cache (Mai 2026)**: orthogonal. Cholesky-Cache spart Faktorisierung; dieser WP spart Matrix-Aufbau. Beide profitieren von gleicher Geometrie-Stabilität.
- **WP-FF-DistMatrix-Sharing**: wenn `srab` zentral geteilt ist, kann der Matrix-Cache auch `srab(triIdx(i,j))` als Distanz-Quelle nehmen statt eigene Berechnung.
- **WP-FF-Packed-Triangular**: `A_offdiag_last` könnte auch packed gespeichert werden — −7 MB.
- **WP-S1 (static_charges)**: alternative, radikalere Strategie — komplettes Skippen von EEQ-Phase-2 in stabilen MD-Trajektorien. Wenn aktiv: dieser WP irrelevant.

XTB-Vorbild gibt es hier nicht — XTB faktorisiert tatsächlich jeden Step neu. Curcuma geht hier über XTB hinaus.
