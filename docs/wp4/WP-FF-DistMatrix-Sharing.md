# WP-FF-DistMatrix-Sharing — Distanz-Matrix einmal pro Step, alle Terme teilen sie

**Status:** 🆕 Vorgeschlagen (Mai 2026)
**Aufwand:** ~1-2 Tage
**Erwarteter Nutzen:** Polymer N=1410 single-thread: **0.5-1.0 s / 100 MD-Steps** durch Eliminierung redundanter `sqrt`-Berechnungen.
**Hebel:** HOCH — die größte verbleibende Architektur-Lücke gegenüber der XTB-Fortran-Referenz.

## Hypothese

Die XTB-Referenz (`external/gfnff/src/gfnff_engrad.F90:175-195`) berechnet als **erste Aktion** in `gfnff_eg`:

```fortran
do i = 1,n
  ij = i*(i-1)/2
  do j = 1,i-1
    k = ij+j
    sqrab(k) = (xyz(1,i)-xyz(1,j))**2 + (xyz(2,i)-xyz(2,j))**2 + (xyz(3,i)-xyz(3,j))**2
    srab(k)  = sqrt(sqrab(k))
  end do
end do
```

Die Arrays `sqrab(N(N+1)/2)` und `srab(N(N+1)/2)` werden danach **von ALLEN Termen** (Repulsion, Coulomb, Dispersion, HB, XB, BATM, EEQ-Phase-2) per `ij = i*(i-1)/2 + j` indiziert konsumiert — kein einziger `sqrt` in den Term-Subroutinen.

Curcuma hingegen hat heute:
- **95 `sqrt`/`.norm()`** Aufrufe in `forcefieldthread.cpp` (verteilt über `CalculateGFNFFBondedRepulsion`, `…NonbondedRepulsion`, `…Coulomb`, `…HydrogenBond`, `…BATM`, …)
- **53 `sqrt`/`.norm()`** in `eeq_solver.cpp` (Phase-2-Distanzmatrix wird unabhängig gebaut)
- Plus CN/dcn in `gfnff_method.cpp` mit eigener Pair-Schleife

Jeder Term iteriert über seine eigene Paar-Liste, berechnet die Distanzen neu. Für Polymer N=1410 single-thread sind das geschätzt 0.5-1.5 ms pro Term × ~10 Terme = 5-15 ms pro MD-Step nur an Distanz-Berechnung.

## Aufgabe

### 1. Shared Distanz-Provider in `gfnff_method.cpp`

Neue Methode `GFNFF::computeSharedDistances()`, aufgerufen einmal pro `Calculation()` direkt nach `prepareCNAndEEQ`:

```cpp
// gfnff_method.h
private:
    mutable Eigen::VectorXd m_shared_sqrab;  ///< Packed triangular: sqrab[i*(i-1)/2 + j] = (r_ij)²
    mutable Eigen::VectorXd m_shared_srab;   ///< Packed triangular: srab[i*(i-1)/2 + j] = r_ij
    mutable bool            m_shared_dist_valid = false;

public:
    /// Compute or refresh shared packed distance arrays. Called once per energy/gradient call.
    /// Pattern matches XTB gfnff_eg lines 175-195.
    void computeSharedDistances() const;
    const Eigen::VectorXd& sharedSqrab() const { return m_shared_sqrab; }
    const Eigen::VectorXd& sharedSrab()  const { return m_shared_srab; }

    /// Helper: packed triangular index for atom pair (i, j), i > j.
    /// Mirrors XTB convention; valid for i != j only.
    static inline int triIdx(int i, int j) {
        if (j > i) std::swap(i, j);
        return i*(i-1)/2 + j;
    }
```

Die Berechnung kann threaded laufen (existierender CxxThreadPool) — analog `dist_worker` in `eeq_solver.cpp:2990-2999`.

### 2. ForceFieldThread konsumiert die Arrays

`ForceField::AutoRanges` (oder ein neuer Setter) übergibt `srab`/`sqrab` per `const std::shared_ptr<const Eigen::VectorXd>` an jede `ForceFieldThread`-Instanz. In den `CalculateGFNFF*Contribution`-Methoden:

```cpp
// vorher
double r = (m_geometry.row(pair.i) - m_geometry.row(pair.j)).norm();

// nachher
double r = (*m_shared_srab)(GFNFF::triIdx(pair.i, pair.j));
```

Migration pro Term:
- `CalculateGFNFFBondedRepulsionContribution` (~6000 Paare polymer)
- `CalculateGFNFFNonbondedRepulsionContribution` (~hunderte-tausende Paare)
- `CalculateGFNFFCoulombContribution` (alle nicht-1,2/1,3-Paare)
- `CalculateGFNFFDispersionContribution` (D4-Paarliste)
- `CalculateGFNFFHydrogenBondContribution` (A-H...B-Tripel: 3 Distanzen pro Tupel)
- `CalculateGFNFFHalogenBondContribution` (analog HB)
- `CalculateGFNFFBatmContribution` (Triplet: 3 Distanzen)

Bond/Angle/Dihedral/Torsion/Inversion sind topologisch — ihre Paare sind kleine Listen (~5000 polymer), `srab`-Lookup hilft, aber der Win ist klein. Niedrige Priorität für diese.

### 3. EEQSolver konsumiert die Arrays

`calculateFinalCharges` baut aktuell sein eigenes `m_phase2_distances(N,N)` (eeq_solver.cpp:2977-3043). Stattdessen:

```cpp
// Drop the local distance computation. Use externally provided packed array.
const Eigen::VectorXd* srab = m_external_srab;  // set via new setter from GFNFF
```

Phase 2 verwendet derzeit `distances(i,j)` ~10 Mal pro Coulomb-Matrix-Eintrag — Lookup-Code muss auf `srab(triIdx(i,j))` umgestellt werden. Diagonalfall (i==j): trivial 0.

Alternative: behalten von `m_phase2_distances` als unpacked-View, **aber gefüllt aus externem srab** statt frisch berechnet — geringerer Refactoring-Aufwand, weniger Speicher-Win, gleiche Compute-Ersparnis.

### 4. CN-Pfad

`calculateCoordinationNumberDerivatives` (gfnff_method.cpp:5676+) iteriert ebenfalls über Paare. Wenn `sqrab` vorliegt: direkt nachschauen statt `(geom(i)-geom(j)).squaredNorm()`. Schon WP-D Stage A/D haben CN/dCN optimiert; diese Stufe ergänzt sich.

### Code-Anker

| Datei | Position | Aufgabe |
|-------|----------|---------|
| `src/core/energy_calculators/ff_methods/gfnff.h` | neue Members + API | `m_shared_sqrab/srab`, `computeSharedDistances`, `triIdx` |
| `src/core/energy_calculators/ff_methods/gfnff_method.cpp` | nach prepareCNAndEEQ in `Calculation()` ~L1611 | `computeSharedDistances()` Aufruf |
| `src/core/energy_calculators/ff_methods/forcefield.h/.cpp` | Setter + Storage | `srab/sqrab`-Reference an Threads |
| `src/core/energy_calculators/ff_methods/forcefieldthread.{h,cpp}` | Inner-Loop-Replacements | `m_geom.row(i)-row(j).norm()` → `(*m_srab)(triIdx(i,j))` |
| `src/core/energy_calculators/ff_methods/eeq_solver.cpp` | calculateFinalCharges:2977-3043 | Distanz-Loop ersetzen |
| `external/gfnff/src/gfnff_engrad.F90:175-195` | Referenz | XTB-Vorbild |

## Akzeptanzkriterien

1. ⚙️ Polymer N=1410, single-thread: `-md -maxtime 1e2` wall-time **≤ 9.0 s** (war 9.7 s).
2. ⚙️ Single-Point caffeine/triose/monosaccharide: Energie bit-identisch oder Δ < 10 nEh (Float-Order ändert sich nicht, weil dieselben Summen).
3. ⚙️ `test_gfnff_numgrad` grün (Gradient unverändert).
4. ⚙️ `ctest -R cli_simplemd` 7/7 pass.
5. ⚙️ Speicher-Footprint pro Step: `m_phase2_distances` (15 MB) entfällt — Win-Win mit Speicher.

## Risiken

1. **Indexierungs-Bugs** — die Fortran 1-basierte Formel `i*(i-1)/2 + j` ist in C++ 0-basiert identisch, aber für i==j ist die packed-Form undefiniert. Lookup-Helper muss `i != j` voraussetzen oder explizit Diagonal-Branch haben.
2. **Wide-Refactor-Risiko** — 95 sqrt/norm-Stellen in einer Datei. Vorsicht-pro-Term Migration empfohlen, nicht alles auf einmal.
3. **Multi-Thread-Konsistenz** — `m_shared_srab` wird einmal pro Energie-Aufruf gefüllt, dann RO geteilt. Keine Schreib-Konflikte solange der Sync-Point (`computeSharedDistances` Aufruf) vor den Term-Threads liegt.
4. **Packed-Index in Inner-Loop** — `i*(i-1)/2` Multiplikation pro Lookup. Bei polymer mit ~10⁵ Lookups/Step ≈ 200 µs zusätzlich. Vernachlässigbar.

## Verifikation

```bash
cd /home/conrad/src/curcuma/release && make -j4

# Korrektheit
./curcuma -sp test_cases/molecules/larger/caffeine.xyz -method gfnff
# Energie sollte bit-identisch zu pre-WP sein

# Performance
/bin/rm -f test_cases/molecules/larger/polymer.topo.json
time ./release/curcuma -md test_cases/molecules/larger/polymer.xyz -method gfnff \
     -print 1 -maxtime 1e2 -threads 1
# Erwartet: ≤ 9.0 s

# Gradient-Regression
ctest -R test_gfnff_numgrad --output-on-failure
ctest -R cli_simplemd     --output-on-failure
```

## Beziehung zu anderen WPs

- **WP-EEQ-Cache** (Mai 2026, gefixt): orthogonal — Cache spart Cholesky-Faktorisierung, dieser WP spart Distanz-Berechnung.
- **WP-FF-Packed-Triangular**: Ergänzung — wenn `sqrab/srab` bereits packed sind, kann Phase-2-EEQ-Matrix auch packed geführt werden (-50% memory).
- **WP-D Stages A-D**: haben CN/dCN inner-loop schon optimiert. Dieser WP ergänzt sie um shared-distance-arrays statt jede Funktion eine eigene Matrix bauen zu lassen.
- **WP-FF-SoA-pair-loops**: SoA-Pattern profitiert von shared-srab, weil Pass-A nur noch index-lookup statt sqrt macht.

XTB-Fortran zeigt: **ein einziges Distanz-Array, viele Konsumenten, alle Term-Subroutinen distanz-frei**. Das ist das Architektur-Ziel.
