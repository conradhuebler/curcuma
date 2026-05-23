# WP-FF-Packed-Triangular — Distanzen/EEQ-Matrix in packed lower-triangular Format

**Status:** 🆕 Vorgeschlagen (Mai 2026)
**Aufwand:** ~1 Tag (folgt auf WP-FF-DistMatrix-Sharing)
**Erwarteter Nutzen:** Polymer N=1410: **−7-8 MB Speicher + 5-10% schnellere Distanz-intensive Loops** durch halbierte Cache-Pressure.
**Hebel:** mittel — additiv zu WP-FF-DistMatrix-Sharing.

## Hypothese

XTB speichert symmetrische N×N-Distanz-Matrizen `sqrab`, `srab`, `eeqtmp` immer **packed triangular** (`N(N+1)/2` Einträge statt N²). Die Indexierung `ij = i*(i-1)/2 + j` mit `i > j` ist konstanter Overhead, dafür halbiert sich der Speicher und damit die L2/L3-Cache-Pressure.

Curcuma's `EEQSolver::m_phase2_distances` ist heute eine **full N×N Eigen::Matrix**. Für N=1410 sind das 1.99M doubles = 15.2 MB. Packed wären es 7.6 MB. Die Distanz-Matrix wird in vielen Funktionen sequentiell gelesen — der zweite Triangle hat per Symmetrie identische Daten und ist nur durch Streaming-Cost teuer.

## Aufgabe

### 1. `PackedTri<T>` Helper-Typ

Educational-First-Constraint: Templates vermeiden, lieber dünner Wrapper:

```cpp
// src/core/energy_calculators/ff_methods/packed_tri.h
#pragma once
#include <Eigen/Dense>

/// Lower-triangular packed buffer for symmetric N×N matrices.
/// Storage: row-major lower triangle, N(N+1)/2 elements.
/// Index: at(i, j) with i >= j returns element at row i, col j.
/// Diagonal access via at(i, i).
/// Mirrors Fortran XTB convention: rab(N(N+1)/2), ij = i*(i-1)/2 + j (1-based)
///                                                = i*(i-1)/2 + j (0-based, i > j)
struct PackedTriDistance {
    Eigen::VectorXd data;
    int natoms = 0;

    void resize(int N) { natoms = N; data.resize(N*(N+1)/2); data.setZero(); }
    int  size() const { return natoms; }

    inline int idx(int i, int j) const noexcept {
        // i > j or i == j; caller is responsible for ordering or use the safe variant
        return i*(i+1)/2 + j;
    }
    inline int idxSafe(int i, int j) const noexcept {
        if (j > i) std::swap(i, j);
        return i*(i+1)/2 + j;
    }
    inline double& at(int i, int j)        noexcept { return data[idxSafe(i, j)]; }
    inline double  at(int i, int j)  const noexcept { return data[idxSafe(i, j)]; }
    inline double& operator()(int i, int j)        noexcept { return at(i, j); }
    inline double  operator()(int i, int j)  const noexcept { return at(i, j); }
};
```

Note: XTB uses `i*(i-1)/2 + j` (1-based). 0-based equivalent is `i*(i+1)/2 + j` for `i >= j`. Diagonal `i==j` maps to `i*(i+1)/2 + i = i(i+3)/2`. Off-diagonal `i > j`: `i*(i+1)/2 + j`.

### 2. EEQ Phase-2 verwendet packed

`eeq_solver.cpp:2977-3043` — `m_phase2_distances` als `PackedTriDistance` umstellen:

```cpp
// Header
mutable PackedTriDistance m_phase2_distances_pkg;

// calculateFinalCharges
m_phase2_distances_pkg.resize(natoms);  // already-correct-size noop after first call
// fill (parallelized via dist_worker as before, but writing packed)
m_phase2_distances_pkg(i, j) = r;       // i > j naturally in loop
```

Alle Lookups später (Coulomb-Matrix-Build, dispatchSolve-Vorbereitung) auf `m_phase2_distances_pkg(i, j)` umstellen.

### 3. Coulomb-Matrix-Build profitiert

Die EEQ-A_nn-Off-Diagonale wird aus `r_ij` per `A_ij = erf(γ·r)/r` gebaut. Mit packed distances + symmetrischer Loop:

```cpp
// vorher (eeq_solver.cpp im Coulomb-build-Block)
for (int i = 0; i < natoms; ++i)
    for (int j = 0; j < natoms; ++j) {
        if (i == j) continue;
        double r = m_phase2_distances(i, j);   // full matrix read
        A(i, j) = erf(gamma_ij * r) / r;
    }

// nachher
for (int i = 0; i < natoms; ++i)
    for (int j = 0; j < i; ++j) {
        double r = m_phase2_distances_pkg(i, j);  // packed read
        double aij = erf(gamma_ij * r) / r;
        A(i, j) = aij;
        A(j, i) = aij;                           // symmetric write
    }
```

Speicher gleich (A ist hier dense, weil Cholesky es braucht), aber Lookup-Pattern besser. Marginal-Win.

### 4. Optional: A_nn auch packed (advanced)

Eigen unterstützt `SelfAdjointView<>` für symmetrische Matrizen + `LLT` darauf. Würde A_nn-Storage von N² auf N(N+1)/2 halbieren (7.6 MB Ersparnis). Aber: Eigen's LLT auf SelfAdjointView ist möglicherweise nicht so optimiert wie auf full Matrix. **Vor Implementierung benchmarken**.

### Code-Anker

| Datei | Position | Aufgabe |
|-------|----------|---------|
| `src/core/energy_calculators/ff_methods/packed_tri.h` (NEU) | — | `PackedTriDistance` Helper |
| `src/core/energy_calculators/ff_methods/eeq_solver.h` | L905 | `m_phase2_distances` → `PackedTriDistance` |
| `src/core/energy_calculators/ff_methods/eeq_solver.cpp` | L2977-3043 + alle `m_phase2_distances`-Reader | Migration |
| `external/gfnff/src/gfnff_engrad.F90:168, 180-194` | Referenz | XTB Storage-Pattern |

## Akzeptanzkriterien

1. ⚙️ Polymer N=1410: Phase-2-Distanz-Buffer 15 MB → 7.6 MB.
2. ⚙️ Energie und Charges bit-identisch (gleiche Float-Operationen, andere Lese-Reihenfolge nicht).
3. ⚙️ Wall-time: messbarer kleiner Win 100-300 ms / 100 MD-Steps (cache-locality).
4. ⚙️ `ctest -R "eeq"` grün (NaCl-cluster, polymer).

## Risiken

1. **Index-Bugs** — packed-Indexierung ist fehleranfällig. Helper-Struct mit `idxSafe` minimiert das. Unit-Test des Helpers vorab schreiben.
2. **Diagonal-Sonderfall** — i==j-Lookup ist semantisch 0 für Distanzen, aber im packed Layout existiert ein Slot. Sollte explizit auf 0 gesetzt werden (oder über `if (i!=j)` Branches gehandhabt).
3. **Eigen-LLT-Kompatibilität** — falls A_nn auf packed/SelfAdjointView umgestellt wird, muss `Eigen::LLT<Eigen::SelfAdjointView<...>>` getestet werden. Falls Performance schlechter: Variante 4 zurücknehmen.

## Verifikation

```bash
cd release && make -j4

# Korrektheit
ctest -R "eeq\|gfnff" --output-on-failure | tail -5

# Memory check (massif oder /proc/PID/status)
./release/curcuma -md test_cases/molecules/larger/polymer.xyz -method gfnff \
   -print 1 -maxtime 1e2 &
PID=$!
sleep 5
grep VmPeak /proc/$PID/status
wait $PID
# Erwartet: ~7 MB weniger als vor dem Patch
```

## Beziehung zu anderen WPs

- **Setzt WP-FF-DistMatrix-Sharing voraus** — Sinn macht es nur, wenn `sqrab/srab` global geteilt werden. Davor ist die packed-Konvertierung ein isolierter EEQ-internal Refactor mit kleinem Win.
- **WP-EEQ-Matrix-Cache** profitiert: zweite Cache-Stufe (A_nn-Off-Diag-Werte) kann auch packed gehalten werden.
