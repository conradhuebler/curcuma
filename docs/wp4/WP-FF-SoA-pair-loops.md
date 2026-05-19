# WP-FF-SoA — Pass-A/B/C SoA-Pattern für FF-Pair-Loops

**Status:** ⚙️ Implementiert (Mai 2026) — neutral für Repulsion, TODO siehe unten
**Aufwand:** ~1 Tag (4 Loops × ~2 h each)
**Erwarteter Nutzen:** Polymer N=1410: −1 bis −3 ms/step über FF-total durch Memory-Traffic-Reduktion.
**Voraussetzung:** WP-D Stage B hat `fast_exp_neg_sq_block` + SoA-Pattern etabliert.

## Messergebnisse (Mai 2026, Polymer N=1410, 4 Threads)

| Build | ms/step | Anmerkung |
|-------|---------|-----------|
| FAST_EXP=ON (SoA aktiv) | 57.9 ms | Bonded + Nonbonded Repulsion mit Pass-A/B/C |
| FAST_EXP=OFF (alter Loop) | 57.2 ms | Original-Code, bit-identisch |
| WP-D-Stage-D Baseline (1d382b7) | 44.2 ms | Referenz vor WP-EEQ-Cache |

**Befund**: SoA-Pattern für Repulsion bringt keinen messbaren Gewinn (~±0.7 ms, Rauschen).
Die Regression 44.2→57 ms stammt aus dem WP-EEQ-Cache-Commit (7e50d27), nicht aus WP-FF-SoA.

**Ursache**: Für die Repulsions-Pair-Listen (klein bis mittel) überwiegt der 3-Pass-Overhead
den exp()-Batch-Gewinn. WP-D Stage B gewann durch SoA, weil der dcn-Loop echte O(N²)
Pairs ohne Cutoff hatte. Repulsion-Listen sind deutlich kleiner (cutoff-gefiltert).

## TODOs — Was noch geprüft werden muss

- [ ] **WP-EEQ-Cache Regression untersuchen**: 44.2 → 57 ms kommt aus Commit 7e50d27.
      EEQ-Cholesky-Cache-Overhead messen (Invalidierungslogik? Cache-miss-Rate? RMSD-check-Kosten?).
      `md_diagnostics_timing` auf Polymer laufen lassen, EEQ-Phase isoliert messen.
- [ ] **Repulsion-Listengröße messen**: Wie viele Paare pro Thread bei Polymer N=1410?
      Falls < 500 Paare/Thread: SoA-Overhead > Gewinn — dann `#ifdef`-Block entfernen.
- [ ] **Mixture N=6200 messen**: Akzeptanzkriterium war "≥1.2× schneller (T=4, mixture)".
      Bei großen Systemen könnte das SoA-Pattern doch greifen. Noch nicht gemessen.
- [ ] **BATM/HBond vorher mit `md_diagnostics_timing` profilieren** bevor SoA dort eingebaut wird.
      BATM hat kein exp() (std::pow statt exp), HBond exp() nur in verschachteltem Nachbarloop.

## Hypothese

WP-D Stage B (Mai 2026) zeigte: für GFN-FF-Inner-Loops mit Distance-Cutoff dominiert **Per-Paar-Memory-Traffic** über die mathematische Operation. Der Pass-A/B/C-Split mit SoA-Scratch in Eigen-Matrix gewann ~1.86 ms in `dcn` step 3 — primär durch Layout-Verbesserung, nicht durch SIMD-exp.

Dieselbe Loop-Form findet sich in vier weiteren FF-Termen in `forcefieldthread.cpp`:

| Term | Loop-Pattern | exp/erf-Anteil | Hebel-Schätzung |
|------|--------------|----------------|-----------------|
| **Repulsion (bonded + nonbonded)** | Pair-Loop mit `exp(-α·r)` | 2 exp pro Paar | 1-2 ms |
| **HBond Gradient-Distribution** | 3er-Tupel-Loop (A-H...B) mit Damping | 1-2 exp pro Tupel | 0.5-1 ms |
| **BATM** | Triplet-Loop mit Tang-Toennies | 3+ exp pro Triplet | 0.5-1 ms |
| **Angle-Damping** | gfnffdampa pro Angle | 1 exp pro Angle | 0.3-0.6 ms |

Jeweils nicht groß — aber **fünf zusammen** sind ~3-5 ms und der Code wird einheitlicher. Pattern ist mechanisch:

```cpp
// Pass A: collect surviving pairs into SoA scratch
Eigen::Matrix<double, Dynamic, K_COLS, RowMajor> buf_scatter(N_max, K_COLS);
std::vector<double> buf_x_sq(N_max), buf_exp(N_max);
int K = 0;
for (const auto& pair : pair_list) {
    // distance check (or cell-list)
    if (skip) continue;
    buf_x_sq[K] = ...;
    buf_scatter(K, 0..K_COLS-1) = ...precomputed factors...;
    ++K;
}

// Pass B: SIMD exp
curcuma::gfnff::fast_exp_neg_sq_block(buf_x_sq.data(), buf_exp.data(), K);

// Pass C: scatter
for (int k = 0; k < K; ++k) {
    double e = buf_exp[k];
    // scatter to local.diag, local.gradient, local.energy
}
```

## Aufgabe

### 1. Repulsion (bonded + nonbonded)

`forcefieldthread.cpp:CalculateGFNFFBondedRepulsionContribution` und
`forcefieldthread.cpp:CalculateGFNFFNonbondedRepulsionContribution`:

```cpp
// Aktueller Inner-Loop
for (const auto& rep : m_gfnff_repulsions_*) {
    double r = ...;
    double exp_val = std::exp(-rep.alpha * r);  // Skalar exp
    double e = ... * exp_val * ...;
    // gradient
}
```

Restrukturieren in 3-Pass-SoA. Helper `fast_exp_block(arg_buffer, out, K)` (Neuvariante von `fast_exp_neg_sq_block`, akzeptiert beliebiges `arg` statt nur `-x²`).

Erweiterung `fast_exp.h`:
```cpp
void fast_exp_block(const double* x, double* out, std::size_t n) noexcept;  // exp(x[k])
void fast_exp_neg_block(const double* x, double* out, std::size_t n) noexcept;  // exp(-x[k])
```

### 2. HBond Gradient-Distribution

`forcefieldthread.cpp:CalculateGFNFFHydrogenBondContribution` enthält eine 3er-Tupel-Loop (A−H...B). Innerhalb der Loop sind drei separate Damping-Funktionen mit eigenen `exp`-Calls.

Pass-A sammelt für jedes Tupel die drei Argumente; Pass-B führt drei `fast_exp_block`-Calls aus (oder einen kombinierten); Pass-C scattert.

### 3. BATM

`forcefieldthread.cpp:CalculateGFNFFBatmContribution` (Body-Above-Three-Mention). Triplet-Loop mit 3-Atom-Geometrie pro Eintrag. Tang-Toennies-Damping enthält `exp`-Calls.

Triplet-Loops sind komplexer als Pair-Loops; lohnen sich nur wenn die `exp`-Anteil signifikant ist. Vorher messen via objdump/perf bevor refactored.

### 4. Angle-Damping (`gfnffdampa`)

In `forcefieldthread.cpp` werden im Angle-Term per Angle 2 `gfnffdampa`-Calls gemacht (jeweils mit `exp`). Für polymer ~5000 Angles × 2 = 10k exp-Calls pro Step → klein (~0.3 ms), aber das ist die einfachste Stelle zum Einbauen.

### Code-Anker

| Datei | Funktion | Zeilen-Bereich (approx) |
|-------|----------|-------------------------|
| `src/core/energy_calculators/ff_methods/forcefieldthread.cpp` | `CalculateGFNFFBondedRepulsionContribution` | Repulsion-Pair-Loop |
| `src/core/energy_calculators/ff_methods/forcefieldthread.cpp` | `CalculateGFNFFNonbondedRepulsionContribution` | Repulsion-Pair-Loop |
| `src/core/energy_calculators/ff_methods/forcefieldthread.cpp` | `CalculateGFNFFHydrogenBondContribution` | A-H-B-Triplet-Loop |
| `src/core/energy_calculators/ff_methods/forcefieldthread.cpp` | `CalculateGFNFFBatmContribution` | Triplet-Loop |
| `src/core/energy_calculators/ff_methods/forcefieldthread.cpp` | `CalculateGFNFFAngleContribution` (inner `gfnffdampa`) | per Angle |
| `src/core/energy_calculators/ff_methods/fast_exp.h` | erweitert um `fast_exp_block`, `fast_exp_neg_block` | neue Interfaces |

## Akzeptanzkriterien

Pro Term-Patch separat verifizieren (nicht alle auf einmal — sonst lokalisiert man Drift nicht):

1. ⚙️ Default `USE_GFNFF_FAST_EXP=OFF`: bit-identisch zu pre-WP.
2. ⚙️ FAST_EXP=ON: per-term Single-Point-Drift ≤ 1 µEh auf caffeine/triose/polymer.
3. ⚙️ `test_gfnff_numgrad` grün (term-spezifisch).
4. ⚙️ Per-Term-Wall-Clock-Test (vor und nach Patch):
   - Repulsion total: ≥ 1.2× schneller (T=4, mixture)
   - HBond: ≥ 1.1× schneller
   - Angle-Damping: messbar (mind. Noise nicht überschritten)

## Risiken

1. **Code-Duplizierung** — fünf separate Loop-Restrukturen ist viel mechanischer Code. Risiko: Inkonsistenzen zwischen den Patches. Lösung: **template/inline helper** `make_simd_loop(...)` der die Pass-A/B/C-Pipeline abstrahiert, plus Loop-spezifische Lambdas. Educational-First-Constraint sagt aber: Templates vermeiden, lieber 5× direkter Code.

2. **Inner-Loop-Kompexität** — manche Loops (HBond, BATM) haben mehr State pro Iteration als der einfache dcn-Loop. Pass-A-Scratch braucht u.U. 6-8 Spalten statt 4. Eigen-Matrix-Layout-Wahl wichtig (RowMajor mit ≥1 Cache-Line pro row optimal).

3. **Kleine Loops keine Wins** — Angle-Damping mit nur 10k exp-Calls/Step ist möglicherweise SCHLECHTER nach Restrukturierung (Pass-A-Overhead > exp-Gewinn). Pre-Messung mit `md_diagnostics_timing` pro Term machen, **bevor** der Patch geschrieben wird. WP-D Stage B v1 war ein Beispiel dafür.

4. **fast_exp-Genauigkeit für Repulsion** — Repulsion-`exp(-α·r)` mit `α ∈ [1, 8] Bohr⁻¹` und r ∈ [0.5, 15] Bohr gibt Argumente bis ~120. Domäne ähnlich zu dcn aber Sign-Flip. `fast_exp_neg_block` muss diese Domäne abdecken. Test entsprechend erweitern.

5. **MD-Stabilität bei kombiniertem Patch** — wenn alle 4 Terme zusammen patched werden: kumulativer Float-Order-Drift kann MD-Trajektorie messbar verschieben. Long-Run-Test (1000+ steps) Pflicht.

## Verifikation

```bash
cd release && cmake -DUSE_GFNFF_FAST_EXP=ON -DUSE_MARCH_NATIVE=ON ..
make -j4

# Pro Term: numgrad + single-point gegenüber pre-patch
./test_cases/test_gfnff_numgrad

# Per-Term Phase-Profile (welche Phase hat sich verändert)
cd test_cases/cli/simplemd/10_gfnff_polymer_md
rm -f input.diag.jsonl
../../release/curcuma -md input.xyz -method gfnff -maxtime 300 -threads 4 \
    -md.dump_frequency 10 -md.md_diagnostics true -md.md_diagnostics_timing true
# input.diag.jsonl enthält per-Term-Times (siehe ff_total Breakdown)

python3 -c "
import json, statistics
r = [json.loads(l) for l in open('input.diag.jsonl') if json.loads(l).get('step',0) >= 50]
# Pro-Term-Phase-Felder (depends on what md_diagnostics_timing exposes)
keys = ['repulsion', 'hbond', 'angle', 'bond', 'ff_total']
for k in keys:
    vals = [x['timing_ms'].get(k, 0) for x in r]
    if vals: print(f'{k:>12}: med={statistics.median(vals):.2f} ms')"

# 1000-Step Wall mit Cache disabled (to isolate FF-Term-effect)
bash compare_baseline.sh 4 gfnff "WP-FF-SoA"
```

## Beziehung zu WP-D Stage B

WP-D Stage B etablierte:
1. `fast_exp_neg_sq_block` AVX2-Helper.
2. Pass-A/B/C-SoA-Pattern mit Eigen-Matrix + counter K.
3. Erkenntnis: Memory-Traffic > exp() für die typischen GFN-FF-Inner-Loops.

Dieser WP wendet die gleiche Methodik mechanisch auf die FF-Term-Pair/Triplet-Loops an. Erwartete Wins pro Term sind klein (~1 ms), aber **kumulativ und konsistent** mit der bestehenden Architektur. Macht das `fast_exp.h`-Helper-Module zu einem etablierten Curcuma-Building-Block statt einer Einmal-Optimierung für dcn.

Kein dramatischer Hebel — eher Code-Konsolidierung mit kleinem Performance-Side-Effect.
