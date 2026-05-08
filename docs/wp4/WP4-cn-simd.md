# WP4 — CN-Berechnung SIMD-vektorisieren (P4e)

**Status:** 🆕 Vorgeschlagen
**Aufwand:** 4–6h
**Erwarteter Nutzen:** CN-Term ~2× (sekundärer Effekt: D4 Gaussian Weights profitieren auch, weil sie auf CN-Werten aufbauen)

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
