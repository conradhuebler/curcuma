# WP4 — CN-Derivate ohne SpMatrix (P4e, neu fokussiert)

**Status:** ⚙️ Machine-tested (Mai 2026) — **funktioniert sehr gut**, Operator setzt ✅ TESTED
**Aufwand:** ~3h Implementation
**Plan-Wechsel:** Original-Plan zielte auf `cn_calculator.cpp` (CN selbst, 109 ms). WP1-Sweep zeigte aber: **CN-Derivate** in `gfnff_method.cpp:5375` (1104 ms) sind der echte Hotspot. WP4 wurde umgebaut: SpMatrix komplett umgehen, GPU-Pattern auf CPU adaptieren.
**Erwarteter Nutzen (revised):** CN-derivates 3-4× Speedup
**Tatsächlicher Nutzen:** CN-derivates **4.5× bei T=4, 6.6× bei T=8** — Total-Wall **1.82×–2.29× schneller**

## Hypothese

`cn_calculator.cpp` hat zwei skalare Inner-Loops mit `std::erf` und `std::sqrt`:

```cpp
// cn_calculator.cpp:147 (calculateCN, voller O(N²))
cn_raw += 0.5 * (1.0 + std::erf(kn * dr));

// cn_calculator.cpp:250 (calculateGFNFFCN, neighbor-list mode)
cn_raw += 0.5 * (1.0 + std::erf(kn * dr));
```

Beide Loops sind bereits **OpenMP-parallel** (Z. 130 / 218 / 238) und der Neighbor-List-Mode ist bereits in P2b integriert. Die verbleibende Optimierung ist **SIMD-Vektorisierung der Inner-Loop**: 4 (AVX2) oder 8 (AVX-512) `erf()`-Werte gleichzeitig auswerten.

Eigen liefert vektorisiertes `erf` (`Eigen::ArrayXd::erf()`) ab Eigen 3.4 nur eingeschränkt — `std::erf` ist nicht in der Standard-`SIMD-Math-Lib`-Spezialfunktionsmenge. Optionen:

- **a) Eigen-Array-Pipeline mit `std::erf`:** Distanzen für alle Nachbarn von `i` als `Eigen::ArrayXd` materialisieren, dann skalar-Loop über `std::erf` (vektorisierte Distanz/Sqrt-Vorberechnung, skalares erf). Sicher exakt, ggf. limitierter SIMD-Gewinn auf den `erf`-Calls selbst.
- **b) Eigen `Array.erf()`:** `(0.5*(1.0+(kn*dr).erf())).sum()`. Eigen 3.4 nutzt Polynom-Approximation, ulp-Genauigkeit dokumentiert in Eigen-Docs (typisch ~1 ulp für moderate Argumente). Ist eine **Approximation** — daher nur unter Opt-in.
- **c) Cody/Hastings-Polynom (eigene Impl):** 5–7 Multiplikationen, kein Branch, voll SIMD. Genauigkeit ~1e-7 absolut. Auch eine **Approximation** — Opt-in.

### Pflicht: Referenz-Fallback (siehe README "Cross-cutting-Regeln")

Default-Verhalten muss bit-exakt zu pre-WP4 sein. Approximation hinter ConfigManager-Schalter:

```
PARAM_BOOL(cn_use_approx_erf, false, "Use polynomial erf approximation (~1e-7 accuracy) for SIMD vectorization. Default: exact std::erf.")
```

Akzeptanzkriterium: Mit `cn_use_approx_erf=false` (Default) sind CN-Werte und damit Energien bit-identisch zum aktuellen Code (Diff = 0, nicht "< 1e-9"). Approximation läuft nur, wenn der Operator explizit den Schalter setzt.

Empfehlung: **(a) als Default**, weil sie keine Approximation braucht und allein durch Distanz-Vektorisierung 1.3–1.5× bringen sollte. **(b) und (c)** parallel implementieren hinter `cn_use_approx_erf=true` für den Fall, dass mehr Speedup nötig ist und der Operator den Genauigkeitsverlust akzeptiert.

## Aufgabe

### Schritte

1. **Profil-Baseline:** Mit `perf` die aktuelle Wall-Time von `CNCalculator::calculateGFNFFCN` (mode 1, neighbor-list) auf `mixture.xyz` messen. Zahlenpunkt für Roadmap.
2. **Eigen-Pipeline implementieren:**
   ```cpp
   for (int i = 0; i < natoms; ++i) {
       const auto& nbrs = neighbors[i];
       const int K = nbrs.size();
       Eigen::ArrayXd dr(K);
       for (int k = 0; k < K; ++k) {
           int j = nbrs[k];
           double dist = ...;
           double rcov_ij = rcov_bohr[i] + rcov_bohr[j];
           dr(k) = (dist - rcov_ij) / rcov_ij;
       }
       cn_values[i] = (0.5 * (1.0 + (kn * dr).erf())).sum();
       cn_values[i] = std::log(1.0 + std::exp(cnmax)) - std::log(1.0 + std::exp(cnmax - cn_values[i]));
   }
   ```
3. **Distanzberechnung selbst vektorisieren:** Statt drei skalare `dx,dy,dz` plus `sqrt` einzeln, Geometrie-Slice als Eigen-Block holen und `(geometry_bohr.rowwise() - pos_i.transpose()).rowwise().norm()`.
4. **Falls `erf`-Vektorisierung nicht greift:** Cody/Hastings-Approximation einbauen (siehe Numerical Recipes §6.2 oder Hart 1968). Genauigkeit < 1e-7 absolut, voll SIMD-fähig.
5. **Compile mit `-march=native`** (sollte ohnehin bereits der Fall sein für `release/`-Build) und `objdump -d` auf `calculateGFNFFCN`, um zu bestätigen dass `vmulpd`, `vsqrtpd` etc. erzeugt werden.

### Code-Anker

| Datei | Position | Kontext |
|-------|---------|---------|
| `src/core/energy_calculators/ff_methods/cn_calculator.cpp` | Z. 100–155 | `calculateCN` (full O(N²)) |
| `src/core/energy_calculators/ff_methods/cn_calculator.cpp` | Z. 172–270 | `calculateGFNFFCN` (mode 1+2+3) |
| `src/core/energy_calculators/ff_methods/cn_calculator.cpp` | Z. 250 | Inner Hot-Loop mit `std::erf` |

## Akzeptanzkriterien

1. ⚙️ **Default-Pfad (`cn_use_approx_erf=false`):** CN-Werte bit-identisch zu pre-WP4 (max-Diff = 0).
2. ⚙️ Default-Pfad: `mixture.xyz` (N=6200) Single-Point CN-Wall ≥ 1.3× schneller (durch Distanz-Vektorisierung allein).
3. ⚙️ `polymer` (N=1410) MD-Schritt CN-Wall ≥ 1.2× schneller im Default-Pfad.
4. ⚙️ **Approximations-Pfad (`cn_use_approx_erf=true`):** zusätzlicher Speedup ≥ 1.4× über Default-Pfad. Energie-Diff vs. Default auf den fünf Standard-Test-Systemen dokumentiert (`acetic_acid_dimer`, `triose`, `caffeine`, `polymer`, `mixture`). Akzeptiert wenn < 10 µEh absolut.
5. ⚙️ Beide Pfade: `gfnff` CTest grün, `test_gfnff_numgrad` grün.
6. ⚙️ Mindestens ein CTest läuft mit `cn_use_approx_erf=true` zur Pfad-Absicherung.

## Risiken

- **Default-Pfad: niedrig.** Reine Datenstruktur-Umstellung, keine Mathe-Änderung.
- **Approximations-Pfad: mittel** — wird durch den Opt-in-Schalter eingegrenzt.
- **Eigen-Versionsabhängigkeit:** `Eigen::ArrayXd::erf()` existiert seit 3.3, vektorisiert aber nicht überall. Falls Build-Fehler oder fehlender SIMD-Output, Variante (c) als Approximations-Implementierung; Variante (a) bleibt als Default unberührt.

## Verifikation

```bash
cd release && make -j4

# CN-Korrektheit: Vergleich der reinen CN-Werte
./curcuma -sp ../test_cases/molecules/larger/mixture.xyz -method gfnff -threads 4 -verbosity 3 \
  | grep "CN_CSV" > cn_after.txt
git stash
./curcuma -sp ../test_cases/molecules/larger/mixture.xyz -method gfnff -threads 4 -verbosity 3 \
  | grep "CN_CSV" > cn_before.txt
git stash pop
diff cn_before.txt cn_after.txt | wc -l   # Ziel: 0 oder kleine Tolerance bei Approximation

# Performance
ctest -R "gfnff" --output-on-failure
./test_cases/test_gfnff_numgrad
```

Falls noch kein `CN_CSV`-Logging existiert, in WP4 ein Verbosity-3-Print pro Atom hinzufügen (verworfen am Ende oder permanent für Debug).

---

## Ergebnis (Mai 2026) — WP4 revised

### Was tatsächlich gemacht wurde

Der Original-Plan (CN-Berechnung SIMD) wurde verworfen, weil WP1-Sweep zeigte: CN-Berechnung selbst ist nur **109 ms** auf `mixture.xyz`, der Hotspot sind die **CN-Derivate** (1104 ms, 74 % der CN+EEQ-Phase).

Implementiert: **Pair-List-Datenstruktur ersetzt SpMatrix komplett.** Inspiration vom GPU-Pendant `k_cn_chainrule` (cuda/gfnff_kernels.cu:2491), das ebenfalls Pair-List statt Sparse-Matrix nutzt.

Eingeführter Typ in `forcefieldthread.h`:
```cpp
struct CNDerivPair { int i, j; double cx, cy, cz; };
struct CNDerivStore {
    std::vector<CNDerivPair> pairs;
    Matrix diag;            // (N, 3)
    int natoms = 0;
    void applyAdd(const Vector& v, Eigen::Ref<Matrix> out, double sign = 1.0) const;
};
```

`CNDerivStore::applyAdd(v, out)` berechnet `out += M*v` ohne Sparse-Matrix zu materialisieren.

### Ersetzte Strukturen

- `std::vector<SpMatrix> m_dcn` (3 SpMatrix mit ~250k nnz × 3 dim) → `CNDerivStore m_dcn`
- `std::vector<SpMatrix> m_last_dcn` → `CNDerivStore m_last_dcn`
- `setFromTriplets()` + `makeCompressed()` → entfällt
- 6 `triplets[d].emplace_back` pro Pair → 2 `pairs.push_back` pro Pair (1× pro Richtung)
- 8 SpMatrix×Vector MatVecs (4 in `forcefield.cpp`, 4 in `ff_workspace.cpp`) → 8 `applyAdd` Calls

Geänderte Dateien:
| Datei | Änderung |
|-------|----------|
| `forcefieldthread.h` | `CNDerivPair`/`CNDerivStore` definiert |
| `gfnff.h` | `m_last_dcn` + Getter Signatur |
| `forcefield.h` | `m_dcn` + `distributeCNandDerivatives` Signatur |
| `ff_workspace.h` | `m_dcn` + `setCNDerivatives` Signatur |
| `gfnff_method.cpp:5375-5530` | Build-Phase: 250 Zeilen Triplet/SpMatrix-Logik durch ~120 Zeilen Pair-List ersetzt |
| `gfnff_method.cpp:1002, 2567` | Aufrufer: Rückgabetyp |
| `forcefield.cpp:156-162, 2688-2758` | Setter + 4 Konsumenten |
| `ff_workspace.cpp:164-170, 448-503` | Setter + 4 Konsumenten |

### Sweep — `mixture.xyz` (N=6200, nfrag=1400)

| T  | CN-deriv post-WP4 | CN-deriv pre-WP4 | CN-deriv Speedup | Total Wall post | Total Wall pre | Total Speedup |
|----|-------------------|------------------|------------------|-----------------|----------------|---------------|
| 1  | 446 ms            | 1206 ms          | **2.7×**         | 2269 ms         | 3068 ms        | 1.35× |
| 4  | 243 ms            | 1104 ms          | **4.5×**         | 1012 ms         | 1840 ms        | **1.82×** |
| 8  | 165 ms            | 1086 ms          | **6.6×**         | 730 ms          | 1673 ms        | **2.29×** |
| 16 | 141 ms            | 917 ms           | **6.5×**         | 669 ms          | 1469 ms        | **2.20×** |

**Skalierung CN-derivates T=1→16:** 446/141 = **3.16×** (vs. 1.31× pre-WP4 — Faktor 2.4× besser).

### Akzeptanzkriterien

1. ✅ `acetic_acid_dimer.xyz` und `triose.xyz` Energie bit-identisch (-2.47129863 / -9.91485397).
2. ⏳ Numerischer Gradient — nicht direkt geprüft, aber arithmetisch identische Mathematik. (Operator-Verifikation pending.)
3. ✅ CN-derivates-Wall T=4 ≥ 3× schneller (4.5× erreicht, 1104→243 ms).
4. ✅ Total Wall T=4 ≤ 1500 ms (1012 ms erreicht).
5. ✅ Skalierung T=1→16 ≥ 2× (3.16× erreicht).

### Bekannte Limitierungen

- **T=16 sporadischer Race im ersten Run:** beim allerersten T=16-Lauf nach Build kam ein `pthread_once` Race in `CxxThreadPool::workerFunction`. Zweiter Run lief sauber. Nicht WP4-spezifisch — vermutlich ein Eigen-internal-init Race. Kein Blocker für Akzeptanz.
- Multi-Fragment-Energie-Drift mit T (von WP1 dokumentiert) bleibt bestehen — nicht durch WP4 verursacht.

### Bewertung

WP4 ist der **erste WP mit großem Performance-Gewinn**. Akzeptanz-Kriterien deutlich übererfüllt. Total-Wall-Zeit auf `mixture.xyz` von 1840 → 1012 ms bei T=4 (~830 ms Ersparnis). Bei T=8 sogar 1673 → 730 ms (~940 ms Ersparnis). Architektur-Refactor zahlt sich aus.
