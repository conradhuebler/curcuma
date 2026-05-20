# WP-G — `m_gradient`-Layout-Optimierung (P4h)

**Status:** ⚙️ Machine-tested (Mai 2026) — **funktioniert technisch, Speedup gering**. Operator setzt ✅ TESTED.
**Aufwand:** ~3h Implementation
**Erwarteter Nutzen (geplant):** Bond cpu-time 170 → ~50 ms (3×), Total Wall 1840 → ~600-700 ms
**Tatsächlicher Nutzen (gemessen):** Total Wall **3-7% schneller** vs. post-WP4 (mixture.xyz). Bond-cpu-Wall im Noise. Hypothese aus WP3 (Cache-Miss-dominiert) hat sich nicht bestätigt.

## Hypothese

WP3 hat aufgedeckt, dass das Eliminieren der Per-Bond Heap-Allokation **keinen** messbaren Speedup brachte. Bond cpu-time bleibt bei ~170 ms / 4800 Bonds = ~35 µs/Bond für eine ~50-FLOP-Auswertung. Die echte Hot-Stelle ist mit hoher Wahrscheinlichkeit:

```cpp
m_gradient.row(bond.i) += g.transpose();   // forcefieldthread.cpp:1001
m_gradient.row(bond.j) -= g.transpose();
```

Begründung:

- `m_gradient` ist `Matrix` = `Eigen::MatrixXd` (definition: `forcefieldthread.h:489`).
- Eigen-Default-Storage-Order ist **ColumnMajor**. Für eine N×3-Matrix mit N=6200 ist Stride zwischen Zeilen-Elementen `N * sizeof(double) = 6200 × 8 = 49.6 KB`.
- L1-Cache typischerweise 32 KB — d. h. die x/y/z-Werte eines Atoms liegen in **drei verschiedenen Cache-Lines, jeweils >32 KB auseinander**. Jede `row(i) += ...`-Operation triggert 3 Cache-Misses.
- Bei 4800 Bonds × 2 Atome × 3 Cache-Misses × ~50 ns = **~140 ms cpu-time** — passt fast exakt zur gemessenen Bond-Wall.

Andere FF-Terme haben dasselbe Problem: Coulomb (mehrere Tausend Pair-Schreibops), Dispersion (~2000 Pair-Ops), Repulsion. Alle nutzen `m_gradient.row(i) +=` mit demselben Cache-Pattern. Eine Layout-Änderung würde **alle** profitieren.

## Aufgabe

### Variante A: RowMajor (preferred)

Datentyp ändern:
```cpp
// Vorher (forcefieldthread.h:489):
Matrix m_gradient;  // = Eigen::MatrixXd, ColumnMajor

// Nachher:
Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> m_gradient;
```

Damit liegen x/y/z eines Atoms **kontigous** im Speicher (24 Bytes), eine `row(i) +=`-Operation triggert höchstens **eine** Cache-Line.

### Variante B: Per-Thread Stack-Akkumulator

Innerhalb der Inner-Loop in einem `std::array<double, 3>` o. ä. akkumulieren, am Loop-Ende einmal in `m_gradient` schreiben. Reduziert Schreibvorgänge bei wiederholtem Zugriff auf dasselbe Atom (z. B. Bond-Loop für ein zentrales Atom mit vielen Nachbarn). Greift weniger als A bei großen, sparsen Bond-Listen, kann aber zusätzlich helfen.

### Empfehlung

**Variante A zuerst** — Layout-Änderung greift überall ohne Code-Logik anfassen zu müssen. Falls B noch zusätzlich nötig: später anhängen.

### Schritte

1. **Typ ändern** in `forcefieldthread.h:489` und `forcefield.cpp` (wenn `m_gradient` dort referenziert wird) auf `Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>`.
2. **Aufrufer prüfen:** wo wird `m_gradient` als `Matrix&`/`MatrixXd&` typisiert übergeben? Eigen ist strict mit Storage-Order — eine `Matrix&` (ColumnMajor) kann nicht auf eine RowMajor-Matrix referenzieren. Mögliche Fixes: typedef `using GradientMatrix = ...`, oder die Funktionssignaturen auf den Konkret-Typ umstellen.
3. **Reduktion in `ForceField::Calculate()`:** Wenn dort die Per-Thread-`m_gradient` aufsummiert wird, muss die Aggregat-Variable denselben Storage-Order haben.
4. **Initialisierung:** `m_gradient = Matrix::Zero(natoms, 3)` ggf. anpassen.
5. **Rückgabewert** `getGradient()`: wenn Curcuma-API `Matrix` (ColumnMajor) erwartet, am Ende kopieren oder mit `Eigen::Map`-View arbeiten. Wahrscheinlich: explizit konvertieren mit `Matrix grad = m_gradient_rowmajor;` für die externe API-Grenze.
6. **Compile** mit `-O3` und prüfen, dass Eigen die Inner-Loop wirklich vektorisiert (objdump auf `CalculateGFNFFBondContribution` muss `vmovupd`/`vfmadd*` zeigen).

## Critical Files

| Datei | Änderung |
|-------|----------|
| `src/core/energy_calculators/ff_methods/forcefieldthread.h:489` | `m_gradient`-Datentyp auf `Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>` umstellen |
| `src/core/energy_calculators/ff_methods/forcefieldthread.cpp` | Lokale Temp-Variablen mit gleichem Typ; `setZero(natoms, 3)`-Aufrufe ggf. anpassen |
| `src/core/energy_calculators/ff_methods/forcefield.cpp` | Reduktion der Per-Thread-Gradients; Setter/Getter ggf. mit Konvertierung |
| `src/core/energy_calculators/ff_methods/forcefield.h` | Falls `m_gradient` auch dort als Member existiert |

**Bestehende Helper, die bleiben:**
- `geom().row(...)` und alle anderen Lese-Zugriffe — unverändert.
- Inner-Loop-Mathematik — unverändert.

## Verification

```bash
cd release && make -j4

# Korrektheit — bit-identisch (Storage-Order ändert nicht die Werte)
./curcuma -sp ../test_cases/molecules/larger/acetic_acid_dimer.xyz -method gfnff -threads 4 -verbosity 1
# Erwartung: -2.47129863 Eh
./curcuma -sp ../test_cases/molecules/larger/triose.xyz -method gfnff -threads 4 -verbosity 1
# Erwartung: -9.91485397 Eh

# Numerischer Gradient
./test_cases/test_gfnff_numgrad
# Erwartung: < 1e-5 Hartree/Bohr

# Performance — alle FF-Terme
rm -f ../test_cases/molecules/larger/mixture.topo.json
./curcuma -sp ../test_cases/molecules/larger/mixture.xyz -method gfnff -threads 4 -verbosity 2 \
  | grep -E "Bond  |Angle  |Dihedral|Dispersion|Coulomb|Repulsion|H-bonds|Force Field|Total energy"
```

## Akzeptanzkriterien

1. ⚙️ `acetic_acid_dimer.xyz` und `triose.xyz` Energie bit-identisch zu pre-WP-G.
2. ⚙️ Numerischer Gradient `test_gfnff_numgrad` < 1e-5 Hartree/Bohr.
3. ⚙️ Auf `mixture.xyz` mindestens **einer** der dominanten FF-Terme (Bond, Coulomb, Dispersion) zeigt cpu-time-Reduktion ≥ 1.8×. Erwartung: Bond 170 → ~50 ms, Coulomb 660 → ~250 ms.
4. ⚙️ Total Wall auf `mixture.xyz` (T=4) ≤ 1500 ms (vs. 1840 ms heute) — d. h. **300+ ms gespart** durch eine reine Layout-Änderung.
5. ⚙️ Eintrag in `docs/GFNFF-PERFORMANCE-ROADMAP.md` Status-Tabelle für P4h auf `⚙️ Machine-tested` mit Vorher/Nachher-Tabelle.
6. Operator setzt ✅ TESTED.

## Risiken

- **Mittel.** Storage-Order-Mismatch zwischen Member-Typ und Funktionssignaturen kann zu Compile-Fehlern führen. Eigen ist strict — `Matrix& m_gradient` (ColumnMajor) kann nicht auf eine RowMajor-Variable referenziert werden.
- **API-Kompatibilität:** `getGradient()` extern wird vermutlich `Matrix` zurückgeben (ColumnMajor-Typ). Conversion-Cost am API-Rand ist gering (eine Kopie pro Energy-Call), aber muss explizit gemacht werden.
- **Numerik:** Storage-Order ändert nicht die Werte, aber **Reihenfolge der Reduktion** über Threads könnte sich minimal verschieben → letzte 2-3 Dezimalstellen-Drift möglich auf großen Systemen. Auf kleinen Molekülen bit-identisch.
- **Eigen-Vektorisierung:** mit `Dynamic, 3, RowMajor` und einer 3-Komponenten-Zeile vektorisiert AVX2/AVX-512 nicht voll (Padding zu 4 Komponenten könnte helfen, dann müssen sich die Zugriffe aber konsistent verhalten). Falls Performance-Gewinn ausbleibt, als Variante C eine `Matrix<double, Dynamic, 4, RowMajor>` mit ungenutzter 4. Spalte — dann kann Eigen vmm-vektorisieren.

## Out of Scope

- Bond-Term-Mathematik selbst (war WP3-Scope, dort kein Hebel gefunden).
- Coulomb-Cutoff (WP6) — orthogonale Optimierung, kann zusätzlich greifen.
- CN-Derivate (separates WP).

## Beziehung zu WP3

WP3 hat das Bond-Loop-Heap-Allocation-Problem strukturell aufgeräumt, aber den Speedup nicht erreicht. WP-G adressiert die echte Ursache (Cache-Miss-Lawine durch ColumnMajor-`m_gradient`). Beide WPs sind komplementär: WP3-Refactor ist Voraussetzung dafür, dass die Inner-Loop sauber ist; WP-G bringt den Performance-Gewinn.

## Vorbedingung

- WP1 abgeschlossen (Threading klar)
- WP3 abgeschlossen (Bond-Loop ist Stack-Vector statt Heap-Matrix — sauberer Code-Ausgangspunkt für die Layout-Migration)

---

## Ergebnis (Mai 2026)

### Implementation

**Typ-Alias** in `forcefieldthread.h`:
```cpp
using GeoGradMatrix = Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>;
```

**Migrierte Member (~30 Variablen):**
- `forcefieldthread.h`: `m_geometry`, `m_gradient`, 10 Component-Gradients (`m_gradient_bond/angle/...`); Per-Component-Getter Return-Typ auf `const GeoGradMatrix&`
- `forcefield.h`: `m_geometry`, `m_gradient`, 3 CN-Korrektionen (`m_bond/disp/coulomb_cn_correction`)
- `ff_workspace.h`: `m_geometry`, `m_result_gradient`, `m_grad_before_cn`, 10 Result-Gradients, `FFAccumulator::gradient` + 10 Komponenten; Getter Return-Typen
- `gfnff.h`: `m_geometry`, `m_geometry_bohr`, `m_gradient`, `GeometryChangeDetector::m_last_geometry`
- `cuda/ff_workspace_gpu.h/.cu`: `setGeometry`-Signatur

**Lambdas + Function-Pointers** angepasst:
- `runWithGradCapture` Lambda (forcefieldthread.cpp:122) — Parameter `GeoGradMatrix&`
- `FFAccumulator::reset` Lambda — Parameter `GeoGradMatrix&`
- `initGradientComponents` Lambda — Parameter `GeoGradMatrix&`
- `sumComponentGradient` (forcefield.cpp:2990) — Function-Pointer Return-Typ `const GeoGradMatrix&`; interne Aggregation als `GeoGradMatrix`, Konversion zu `Matrix` am API-Rand

**Externe API stable:** `ForceField::Gradient()`, `getGradient()` etc. geben weiter `Matrix` (ColumnMajor) zurück. Eigen konvertiert implizit beim by-value Return. Konversionskosten 50 KB pro Energy-Call — vernachlässigbar.

**ALPB-Boundary:** `m_solvation->addGradient` erwartet ColumnMajor `Matrix&`. Konversion in `gfnff_method.cpp:1471` mit lokalem Temp-Buffer. ALPB nicht im Hot-Path.

### Korrektheit

T=1 (deterministisch): ✅ bit-identisch zu pre-WP-G
- `acetic_acid_dimer.xyz` T=1: -2.47129863 Eh
- `triose.xyz` T=1: -9.91485397 Eh
- `mixture.xyz` T=1: -856.28814787 Eh

T=4: triose zeigt Multi-Thread-Drift zwischen 3 Werten (-9.91485397 / -9.91337978 / -9.91463123). **Verified pre-WP-G hat denselben Drift** (5 Läufe pre-WP-G zeigen dieselbe Verteilung) — bestehende Race-Condition, **nicht durch WP-G eingeführt**. Eigenes WP für die Race-Untersuchung nötig.

### Performance — `mixture.xyz` (N=6200)

| T | CN-deriv post-WPG | post-WP4 | Δ | Total Wall post-WPG | post-WP4 | Δ |
|---|-------------------|----------|---|---------------------|----------|---|
| 1 | 426 ms | 446 ms | -4% | 2099 ms | 2269 ms | **-7%** |
| 4 | 232 ms | 243 ms | -5% | 967 ms | 1012 ms | **-4%** |
| 8 | 174 ms | 165 ms | +5% | 708 ms | 730 ms | **-3%** |

WP-G bringt **3-7% Total-Wall-Speedup** — deutlich unter den erwarteten 40-50%. Bond-cpu-time bleibt im Run-zu-Run-Noise.

### Warum die Hypothese nur teilweise zutrifft

1. **Eigen Expression-Templates verstecken den Storage-Order:** `m_gradient.row(i) +=` wird intern auf Element-Zugriffe abgebildet — der Compiler optimiert beide Layouts zu ähnlichen Maschinencode-Mustern.
2. **N×3 ColumnMajor mit 6200 Atomen:** Stride zwar 49.6 KB, aber moderner Hardware-Prefetcher cached die drei x/y/z Cache-Lines effektiv beim sequentiellen Atom-Iteration-Pattern.
3. **Hot-Pfad ist nicht eindeutig identifizierbar:** Bond + Coulomb + Dispersion + Repulsion verteilen den Anteil, kein Term ist >25% der FF-Pool-Phase.
4. **Multi-Thread-Reduktion war schon kontiguös:** `m_gradient += t->Gradient()` mit Storage-Order-Crossing kostet 50 KB Konversion — Eigen optimiert das fast gleich gut wie ohne Konversion.

### Bewertung

**Behalten,** weil:
- Architektur-Konsistenz mit GPU-Pfad (`ff_workspace_gpu.h:81` Doku sagt schon "row-major")
- Saubere Datenstruktur — `N×3` mit fester Spaltenzahl 3 ist semantisch korrekter als `Dynamic×Dynamic`
- Voraussetzung für künftige SIMD-Inner-Loops, falls direkte Vektorisierung implementiert wird
- Kein Performance-Verlust, kleine Verbesserung (3-7%) ohne Risiko

**Aber:** Die WP3-Hypothese "Cache-Miss-Lawine in `m_gradient.row()`" hat sich **nicht bestätigt**. Der echte Bond-Hotspot bleibt unaufgeklärt. Weitere Optimierung würde tiefere Profile-Analyse benötigen (perf c2c für False Sharing, perf stat für IPC/Cache-Miss-Rate).

### Status für Roadmap

P4h: ⚙️ Machine-tested. Layout-Migration sauber durchgeführt, Korrektheit auf T=1 bit-identisch, T=4-Multi-Thread-Drift pre-existiert. Total-Wall mixture.xyz T=4: 1012 → 967 ms (4% schneller). **Kein dramatischer Speedup**, aber strukturell wertvoll und ohne Regression.
