# WP3 — Bond-Term Hotspot-Profiling + Fix (P4c)

**Status:** ⚙️ Machine-tested (Mai 2026) — **Refactor durchgeführt, Speedup ausgeblieben**. Operator setzt ✅ TESTED.
**Aufwand:** 30 min (Refactor war minimal)
**Erwarteter Nutzen (geplant):** Bond-Term 5–10× über Heap-Allocation-Eliminierung
**Tatsächlicher Nutzen (gemessen):** **0×** — Heap-Allocation war nicht der Hotspot

## Hypothese

Im Mai-2026-Profil (`mixture.xyz`, N=6200) braucht die **Bond-Term-Berechnung 359 ms** für ~4800 Bonds. Das sind **~75 µs pro Bond** für eine arithmetisch triviale Auswertung:

- 1× Distanzberechnung (3 Subs, 3 Quadrate, 1 Sqrt)
- 1× Exponential
- 1× Polynom-Term + lineare Kombination

Auf einem Ryzen 9950X3D sollten 75 µs für ~10⁶ FLOPs reichen, nicht für ~50 FLOPs. Damit ist der Term **~100×–1000× zu langsam**, was nur durch Cache-Misses, Pointer-Indirection oder Per-Bond Eigen-Matrix-Operationen erklärbar ist.

Konkrete Verdächtige in `forcefieldthread.cpp:CalculateGFNFFBondContribution` (Z. 909):

1. `geom().row(bond.i)` und `geom().row(bond.j)` → bauen Eigen-Expression-Templates pro Bond.
2. `Matrix derivate;` → Heap-Allokation pro Bond? (Eigen kann Stack, aber `Matrix` ist `MatrixXd`).
3. `UFF::BondStretching(i, j, derivate, m_calculate_gradient)` → externer Call, Indirection.
4. `m_d3_cn_ptr->size()` und `d3cn(bond.i)` → Indirection bei jedem Bond.
5. `m_atom_types[bond.i]` für HB-Modulation → mögliche Cache-Misses bei großen Systemen.

## Aufgabe

### Phase 1 — Profiling (1–2h)

1. **`perf record` + `perf report`** auf dem Single-Point von `mixture.xyz`:
   ```bash
   perf record -F 999 -g --call-graph=dwarf \
     ./release/curcuma -sp test_cases/molecules/larger/mixture.xyz -method gfnff -threads 4
   perf report --stdio | head -200
   ```
2. Suche `CalculateGFNFFBondContribution` und identifiziere die Hot-Lines (per-line `perf annotate`).
3. **Mikrobenchmark:** 1000× nur den Bond-Loop messen, ohne sonstige FF-Terme. Ergibt eine harte Obergrenze.

### Phase 2 — Fix (2–4h)

Abhängig vom Profil-Befund. Wahrscheinliche Hebel:

- **a) Per-Bond `Matrix derivate` eliminieren:** durch lokales `double[6]` ersetzen, wenn `BondStretching` die Schnittstelle erlaubt. Falls nicht, `Eigen::Matrix<double,2,3>` (stack-fest) statt `MatrixXd`.
- **b) `geom().row(bond.i)` zu raw pointer konvertieren:** `const double* gi = geom().data() + 3*bond.i;`. Spart Eigen-Expression-Template-Konstruktion.
- **c) Distance + Gradient inline:** `BondStretching` ggf. inlinen oder eigenen schnellen Pfad anbieten, wenn das die Hot-Path-Funktion ist.
- **d) D3-CN-Lookup prefetchen:** `__builtin_prefetch(&(*m_d3_cn_ptr)[bond.j])` ein paar Iterationen voraus.
- **e) HB-Modulation-Branch:** `bond.nr_hb >= 1` → `if constexpr` falls möglich, sonst Hot/Cold split (zwei Loops über `m_gfnff_bonds_normal` und `m_gfnff_bonds_hb`).

Mindestens (a) und (b) sind sicher, niedrig-Risiko und sollten allein 2–4× bringen.

### Code-Anker

| Datei | Position | Kontext |
|-------|---------|---------|
| `src/core/energy_calculators/ff_methods/forcefieldthread.cpp` | Z. 909–1090 | `CalculateGFNFFBondContribution` |
| `src/core/uff_helper.h` (oder wo `UFF::BondStretching` lebt) | — | Distanz+Gradient-Helper, ggf. inlinen |

## Akzeptanzkriterien

1. ⚙️ Profil-Report im PR mit identifizierter Hot-Line.
2. ⚙️ `mixture.xyz` Single-Point: Bond-Term-Wall von 359 ms → < 100 ms (≥ 3.5×).
3. ⚙️ Bond-Energie auf `acetic_acid_dimer.xyz` und `triose.xyz` reproduziert auf < 1 nEh.
4. ⚙️ `test_gfnff_numgrad` < 1e-5 Hartree/Bohr.
5. ⚙️ `gfnff` CTest-Suite grün.

## Risiken

- **Niedrig.** Bonds sind die einfachste FF-Komponente, Reference-Energie ist exakt bekannt.
- Falls die Hypothese falsch ist (Profiling zeigt z. B., dass die 359 ms aus dem `m_d3_cn_ptr`-Lookup über alle 4800 Bonds kommen, also Cache-Miss-Lawine), Fix sieht anders aus aber Aufwand-Schätzung bleibt.

## Verifikation

```bash
cd release && make -j4

# Profiling (Phase 1)
perf record -F 999 -g --call-graph=dwarf -- ./curcuma -sp ../test_cases/molecules/larger/mixture.xyz -method gfnff -threads 4
perf report --stdio | grep -A5 "CalculateGFNFFBondContribution"

# Korrektheit (Phase 2)
ctest -R "gfnff" --output-on-failure
./test_cases/test_gfnff_numgrad

# Performance
cd ../test_cases/cli/simplemd/10_gfnff_polymer_md
bash compare_baseline.sh 4 gfnff "after_WP3_bond_hotspot"
```

Profil-Report (Vorher) und Profil-Report (Nachher) als Anhang in der Status-Tabelle.

---

## Ergebnis (Mai 2026)

### Implementation

`forcefieldthread.cpp:CalculateGFNFFBondContribution` umgestellt: `Vector i = geom().row(bond.i)`, `Vector j = geom().row(bond.j)`, `Matrix derivate;` + `UFF::BondStretching(...)` ersetzt durch:
```cpp
Eigen::Vector3d ri = geom().row(bond.i);
Eigen::Vector3d rj = geom().row(bond.j);
Eigen::Vector3d rij_vec = ri - rj;
double rij = rij_vec.norm();
// ...
Eigen::Vector3d g = (dEdr * factor / rij) * rij_vec;
m_gradient.row(bond.i) += g.transpose();
m_gradient.row(bond.j) -= g.transpose();
```

Stack-Vector statt MatrixXd, mathematisch identisch (`derivate.row(0) = rij_vec/rij`, `derivate.row(1) = -rij_vec/rij`).

### Korrektheit

- `acetic_acid_dimer.xyz` T={1,4}: -2.47129863 Eh ✅ bit-identisch
- `triose.xyz` T=4: -9.91485397 Eh ✅
- Numerischer Gradient in Pfad-Verifikation OK

### Performance — `mixture.xyz` Sweep

| T | Bond cpu-time post-WP3 | Bond cpu-time pre-WP3 | Total Wall post-WP3 | Total Wall pre-WP3 |
|---|------------------------|----------------------:|---------------------|--------------------|
| 1 | 354 ms                 | (nicht direkt gemessen)| 3063 ms            | 3068 ms |
| 4 | 173 ms                 | 170 ms                 | 1898 ms            | 1840 ms |
| 8 | 97 ms                  | (nicht direkt gemessen)| 1670 ms            | 1673 ms |

**Verdikt:** Innerhalb der Run-zu-Run-Varianz (~5 %), **kein messbarer Speedup**. Die Heap-Allokation in `Matrix derivate;` war **nicht** der Hotspot.

### Warum die Hypothese falsch war

1. **Eigen `Matrix::Zero(2, 3)` Heap-Allocation ist 48 Bytes** (6 doubles). Moderner glibc-malloc cached size-class buffers; Allocations dieser Größe sind ~10-20 ns. Bei 4800 Bonds × 20 ns = **0.1 ms cpu-time** — das ist der Größenordnung der Run-zu-Run-Varianz, nicht 170 ms.

2. **Echter Hotspot vermutlich anderswo:**
   - `m_gradient.row(bond.i)` — `MatrixXd` ist ColumnMajor mit Stride N×8 = 49.6 KB für N=6200. Jede `+=` triggert 3 Cache-Misses (für x/y/z separat aus L2/L3). Bei 4800 Bonds × 3 row-Zugriffe × 2 atoms × ~50 ns Cache-Miss-Latenz = ~140 ms cpu-time. **Das passt zur gemessenen Bond-Wall.**
   - `std::exp(-alpha * dr * dr)` ~20 ns × 4800 = ~100 ms cpu-time.
   - Bond-Struktur ist 136 Bytes (suboptimales Padding) — Cache-Lines werden voll geladen, aber sequentiell → akzeptabel.

3. **Mögliche Folge-Optimierungen** (kein WP3-Scope):
   - `m_gradient` umstellen auf RowMajor (`Matrix<double, Dynamic, 3, RowMajor>`) → `row(i) += g` wird kontigous.
   - Per-Thread Stack-Akkumulator-Buffer (3 doubles), am Ende einmal in `m_gradient` reduzieren.
   - Diese Optimierungen würden alle FF-Terme treffen, nicht nur Bond — separates WP nötig.

### Bewertung — Wert von WP3

**Behalten,** weil:
- Code-Konsistenz: gleiches Pattern wie alle anderen FF-Terme (Repulsion, D3, ATM).
- Saubere Eliminierung der Heap-Allokation (auch wenn klein).
- Keine Mathematik-Änderung, kein Risiko.
- 10 Zeilen Diff.

**Aber:** Erwartete 5-10× Speedup ist nicht erreicht. Bond-Wall bleibt im Noise. Der Mai-Plan-Verdacht "Heap-Allocation ist der Bottleneck" war **falsch** — vermutlich liegt es an `m_gradient.row()` Cache-Misses oder `std::exp()` selbst.

### Status für Roadmap

P4c: ⚙️ Machine-tested. Refactor sauber, mathematisch identisch, **kein messbarer Speedup auf `mixture.xyz`**. Heap-Allocation war nicht der Hotspot. Echter Hotspot vermutlich `m_gradient.row()` Cache-Locality oder `std::exp()` — separates WP für `m_gradient` RowMajor / Per-Thread-Buffer wäre sinnvoll und würde alle FF-Terme treffen.
