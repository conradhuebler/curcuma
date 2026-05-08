# WP5 — D4 Gaussian-Weights vektorisieren (P4f)

**Status:** 🆕 Vorgeschlagen
**Aufwand:** 3–4h
**Erwarteter Nutzen:** D4-Weights ~2× (auf Top von schon vorhandener Threading-Parallelisierung)

## Hypothese

In `d4param_generator.cpp:986–1072` (`precomputeGaussianWeights`) läuft pro Atom eine kleine Inner-Loop über `MAX_REF=7` Referenzzustände:

```cpp
for (int ref = 0; ref < nref && ref < MAX_REF; ++ref) {
    double cni_ref = m_refcn[elem_i][ref];
    double diff_cn = cni - cni_ref;
    weights[ref] = std::exp(-wf * diff_cn * diff_cn);  // skalar
    sum_weights += weights[ref];
}
```

Mai-2026-Profil (`mixture.xyz`, 6200 Atome, MAX_REF=7): 144 ms / (6200 × 7) ≈ **3.3 µs pro Atom-Reference-Paar**, was wieder grotesk hoch ist (sollte < 100 ns sein). Vermutlich:

1. Per-Atom `std::vector<double> weights(nref, 0.0)` Heap-Allokation.
2. `m_refcn[elem_i][ref]`: 2× Vector-of-Vector-Indirection.
3. `std::exp` skalar (kein SIMD über 7 Werte aber 4× SIMD über 4 Atome wäre ein Gewinn).

Existierende Parallelisierung (Mar 2026) ist über **Atome (outer loop)** verteilt — jeder Thread macht seine eigene Atomliste. SIMD-Vektorisierung kann **innerhalb** der Inner-Loop wirken (Eigen `Array<double, 7, 1>`) **oder** **über die Outer-Loop** (4 Atome simultan).

## Aufgabe

### Schritte

1. **`m_refcn` flatten:** Statt `std::vector<std::vector<double>>` ein `Eigen::Matrix<double, Dynamic, MAX_REF>` mit Padding (MAX_ELEM × MAX_REF = 118 × 7 ≈ 800 doubles, ~6 KB). Cache-fest, eine Indirection.
2. **Per-Atom Inner-Loop als Eigen-Array:**
   ```cpp
   Eigen::Array<double, MAX_REF, 1> diff = cni - m_refcn_flat.row(elem_i).transpose().head<MAX_REF>();
   Eigen::Array<double, MAX_REF, 1> w = (-wf * diff.square()).exp();
   double sum = w.head(nref).sum();
   if (sum > 1e-10) w.head(nref) /= sum;
   else { w.setZero(); w(0) = 1.0; }
   ```
   Eigen kann die `exp`-Vektorisierung für `Array<double,4,1>` und `Array<double,8,1>` automatisch (AVX2/AVX-512). Für 7 Elemente nutzt es 4-er- oder 8-er-Vektoren mit Maskierung — Compile-Output verifizieren.
3. **Heap-Allokation eliminieren:** `m_gaussian_weights[i]` als `Eigen::Matrix<double, MAX_REF, 1>` statt `std::vector<double>`. Überall wo `m_gaussian_weights[i][ref]` gelesen wird, an die neue Datenstruktur anpassen.
4. **`computeGaussianWeightDerivatives` analog:** wenn dieselbe Per-Atom-Loop dort existiert, gleiches Pattern.

### Code-Anker

| Datei | Position | Kontext |
|-------|---------|---------|
| `src/core/energy_calculators/ff_methods/d4param_generator.cpp` | Z. 986–1072 | `precomputeGaussianWeights` |
| `src/core/energy_calculators/ff_methods/d4param_generator.cpp` | Z. 1026–1036 | Inner Hot-Loop mit `std::exp` |
| `src/core/energy_calculators/ff_methods/d4param_generator.h` | — | `m_refcn`, `m_gaussian_weights` Datentypen |
| `src/core/energy_calculators/ff_methods/cuda/` | — | GPU-Variante (referenz, nicht ändern) |

### Welche Stellen hängen am Datentyp `m_gaussian_weights`?

Vor Patch:
```bash
grep -rn "m_gaussian_weights" src/core/energy_calculators/ff_methods/ | grep -v cuda
```

Wenn ≤ 5 Lese-Stellen: Refactor in einem Schritt. Falls > 5: Adapter-Pattern (Accessor-Methode behalten, intern umstellen).

## Akzeptanzkriterien

1. ⚙️ `mixture.xyz` Single-Point: D4-Weights-Wall (Phase im Profile) von ~144 ms → < 80 ms.
2. ⚙️ Gaussian-Weights `weights[i][ref]` numerisch identisch zu pre-WP5 auf 1e-12 absolut.
3. ⚙️ Caffeine, triose D4-Energie reproduziert auf 1 nEh.
4. ⚙️ `gfnff` CTest, `test_gfnff_numgrad` grün.
5. ⚙️ GPU-Variante in `cuda/` unverändert und weiterhin funktional (Phase-6-Pipeline läuft).

## Risiken

- **Niedrig–Mittel.** D4-Weights gehen direkt in C6-Werte und damit in Dispersionsenergie und Gradient. Ein Bit-Drift hier verändert alle Dispersionsenergien.
- **Datenstruktur-Refactor:** `std::vector<double>` → `Eigen::Matrix` — alle Konsumenten müssen mitziehen. Fehleranfällig wenn versteckte `[]`-Zugriffe existieren.
- **Bound-checks:** `MAX_REF` ist eine Compile-Konstante (gut), aber `nref` variiert pro Element. Sicherstellen dass Read-Bereich `[0..nref)` bleibt.

## Verifikation

```bash
cd release && make -j4

# Numerische Identität
./curcuma -sp ../test_cases/molecules/larger/caffeine.xyz -method gfnff -verbosity 3 \
  > caffeine_after.log 2>&1
git stash
./curcuma -sp ../test_cases/molecules/larger/caffeine.xyz -method gfnff -verbosity 3 \
  > caffeine_before.log 2>&1
git stash pop
diff <(grep "Dispersion" caffeine_before.log) <(grep "Dispersion" caffeine_after.log)  # < 1 nEh

# Gradient
./test_cases/test_gfnff_numgrad
ctest -R "gfnff" --output-on-failure
```

Inner-Loop SIMD-Verifikation: `objdump -d release/curcuma | grep -A30 precomputeGaussianWeights | grep -E "vmulpd|vexp|vfmadd"` muss Vektor-Ops zeigen.
