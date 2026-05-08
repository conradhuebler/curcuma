# WP-G — `m_gradient`-Layout-Optimierung (P4h)

**Status:** 🆕 Vorgeschlagen
**Aufwand:** 4–8h (mehrere Aufrufer betroffen)
**Erwarteter Nutzen:** **Größter Hebel im FF-Pool laut WP3-Befund** — Bond, Coulomb, Dispersion, Repulsion, ATM, H-Bonds profitieren alle. Schätzung Bond cpu-time 170 → ~50 ms (3×), Coulomb 660 → ~250 ms, Dispersion 168 → ~60 ms.

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
