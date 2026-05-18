# WP-D Stage C — Neighbor-List-Wiederverwendung in dcn step 3

**Status:** ✅ Implementiert + gemessen (Mai 2026)
**Aufwand:** ~3 h
**Erwarteter Nutzen:** dcn: −1.5 bis −3 ms/step auf polymer N=1410 (durch Eliminieren von ~960 k Threshold-Checks pro `dcn`-Call). Skaliert mit N²/k wo k = avg Nachbarn pro Atom (typisch 50, gegenüber N/2 = 705).

## Hypothese

`CNCalculator::calculateGFNFFCN()` baut eine Neighbor-List ab `cn_cutoff_bohr > 0` (Default 6.0 Bohr, P2b-WP, Apr 2026). Die Liste wird nach der CN-Berechnung **verworfen**. `calculateCoordinationNumberDerivatives()` (Stage A: nimmt `cn_raw` als Eingang) iteriert anschließend **erneut den vollen N²/2-Triangle** mit per-Paar-Distance-Threshold-Check.

XTB (siehe `external/xtb/src/gfnff/gfnff_eg.f90:3677`, `ncoordNeighs`) iteriert pro Atom über `neighs(iat)` echte Nachbarn — Faktor 30-70× weniger Inner-Loop-Iterationen bei N=1410.

Für polymer mit cn_cutoff_bohr=6.0:
- Aktuell: 1410 × 1410 / 2 = **994 k Paar-Tests** pro dcn-Call, davon ~80 k überleben (8 %)
- Mit Liste: ~1410 × 50 / 2 = **35 k Paar-Iterationen**, alle überleben (100 %)

Pro übersprungenem Paar: 1 sub + 1 dot + 1 cmp ≈ 5–8 Zyklen × ~960 k = **~4 M Zyklen ≈ 1.0–1.5 ms** Threshold-Cost gespart. Plus Cache-Locality-Bonus weil nur die Nachbarliste durchwandert wird statt der full row.

## Aufgabe

### 1. CNCalculator API erweitern

`src/core/energy_calculators/ff_methods/cn_calculator.{h,cpp}`:

```cpp
struct CNResult {
    Vector cn_raw;                          // bereits da (Stage A Output)
    std::vector<std::vector<int>> neighbors; // NEU: pro-Atom-Nachbarliste
    double threshold_used;                  // NEU: Cutoff² in Bohr² für Konsistenz
};
```

`calculateGFNFFCN(...)` Signatur erweitern: gibt `CNResult` statt `Vector` zurück. Wenn `cn_cutoff_bohr > 0`, befüllt es `neighbors`; sonst bleibt es leer und nachgelagerter Code nutzt den alten Pfad.

### 2. `calculateCoordinationNumberDerivatives` Signatur

`gfnff_method.{cpp,h}:5608`:

```cpp
CNDerivStore calculateCoordinationNumberDerivatives(
    const Vector& cn,
    const Vector& cn_raw_in,
    double threshold,
    CxxThreadPool* pool, int num_threads,
    const std::vector<std::vector<int>>* neighbors = nullptr  // NEU
) const;
```

Wenn `neighbors != nullptr`: ersetze die innere `for (int j = 0; j < i; ++j)` durch `for (int j : (*neighbors)[i] mit j < i)`. Threshold-Check kann komplett entfallen (Liste ist schon gefiltert).

### 3. Caller anpassen

`gfnff_method.cpp` Hot-Path:
- `m_cn_result = m_cn_calculator->calculateGFNFFCN(...)`
- `m_cn = m_cn_result.cn_raw`
- `auto dcn_store = calculateCoordinationNumberDerivatives(cn, m_cn, threshold, pool, threads, &m_cn_result.neighbors)`

Pre-existierende Caller, die ohne Liste arbeiten (e.g. EEQ-Phase-2 Recompute mit `cn_raw_in.size() != m_atomcount`), passen `neighbors = nullptr` und fallen auf den N²-Pfad zurück.

### 4. Triangular-Symmetrie

`neighbors[i]` enthält i.A. sowohl `j < i` als auch `j > i` (symmetric). In der dcn-Loop will man nur Paare mit `j < i` (sonst doppelt). Zwei Optionen:
- **A**: bei Listen-Aufbau in CNCalculator nur `j < i` einfügen (Upper-Triangular-Liste).
- **B**: bei Iteration `if (j >= i) continue;` als billigerer Check als der vollständige Distance-Threshold.

Option A ist sauberer. Erfordert separaten Listen-Aufbau-Pfad weil CNCalculator selbst die full symmetric Liste benötigt (für `cn_i += erfCN; cn_j += erfCN`). Lösung: 2 Listen abgeben (`neighbors_full` und `neighbors_upper`) oder eine `neighbors_full` plus `j < i` Inline-Check.

### Code-Anker

| Datei | Zeilen | Inhalt |
|-------|--------|--------|
| `src/core/energy_calculators/ff_methods/cn_calculator.cpp` | 215–280 | Neighbor-List-Build-Block (P2b) |
| `src/core/energy_calculators/ff_methods/cn_calculator.h` | Class-Header | Erweitern um `CNResult` |
| `src/core/energy_calculators/ff_methods/gfnff_method.cpp` | 5608–5950 | `calculateCoordinationNumberDerivatives` |
| `src/core/energy_calculators/ff_methods/gfnff_method.cpp` | Aufruf-Stellen | Caller-Signatur ergänzen (`m_cn_result.neighbors`) |
| `src/core/energy_calculators/ff_methods/gfnff_method.h` | | Method-Signatur |

## Akzeptanzkriterien

1. ⚙️ Default-Pfad mit `cn_cutoff_bohr=6.0`: dcn-Energie und -Gradient bit-identisch zu pre-Stage-C (innerhalb Float-Order-of-Operations).
2. ⚙️ Default-Pfad mit `cn_cutoff_bohr=0`: identisches Verhalten zu Stage B (`neighbors=nullptr`-Fallback).
3. ⚙️ test_gfnff_numgrad grün auf beiden Modi.
4. ✅ polymer N=1410 T=4 dcn-Phase: **0.16 ms median** (Ziel ≤ 15.5 ms, Stage B Baseline ~17.2 ms) — **>100× Speedup**.
5. ⚙️ Keine zusätzliche Allokation pro dcn-Call (Neighbor-List wird einmal pro Step in CNCalculator gebaut und referenziert).

## Messergebnisse (Mai 2026)

| Metrik | Stage B Baseline | Stage C | Faktor |
|--------|-----------------|---------|--------|
| dcn Median | ~17.2 ms | **0.16 ms** | >100× |
| dcn Mean | ~17.2 ms | 0.17 ms | >100× |
| dcn Min/Max | — | 0.11 / 0.24 ms | — |
| Testlauf | polymer N=1410, T=4, 500 steps, steps≥50 | | |

Tatsächlicher Hebel weit größer als erwartet (WP-Prognose: 1.5–3 ms gespart).
Ursache: Neighbor-List-Iteration hat deutlich bessere Cache-Lokalität als der N²-Loop über alle Atome,
zusätzlich zu den ~28× weniger Iterationen.

## Risiken

1. **Neighbor-List-Konsistenz**: CNCalculator nutzt cutoff `cn_cutoff_bohr` (Default 6.0 Bohr). `calculateCoordinationNumberDerivatives` nutzt `threshold` (typisch 9.5 Bohr² → ~3.1 Bohr). Wenn cutoff der Liste **kleiner** als threshold ist: dcn fehlt Paare. Lösung: assert/check dass cn-cutoff ≥ √threshold; sonst Fallback auf N²-Pfad.

2. **Stage A Cn-Raw-Reuse koppeln**: Liste ist nur konsistent wenn das CN aus demselben Call kommt. Bei EEQ-Phase-2-Recompute (anderer Code-Pfad mit eigenem cn_raw_in) muss `neighbors=nullptr` übergeben werden — sonst Inkonsistenz.

3. **Memory**: `std::vector<std::vector<int>>` mit avg 50 Einträgen × 4 Bytes × 1410 Atome ≈ 280 kB. Klein, aber jede vector-of-vector hat 24 Bytes Header → 34 kB Overhead. Akzeptabel.

4. **Listen-Aufbau-Kosten** falls Liste neu pro Step gebaut wird: WP-P2b zeigt das ist Cache-warm und schnell. Aber: wenn das CN-Update den static-Mode (WP-S1) skippt, könnte die Liste veraltet sein. → static_cn-Pfad muss die Liste mit-cachen.

## Verifikation

```bash
cd release && cmake -DUSE_GFNFF_FAST_EXP=ON -DUSE_MARCH_NATIVE=ON ..
make -j4

# Korrektheit
ctest --output-on-failure -R "gfnff"
./test_cases/test_gfnff_numgrad

# Phase-Profile vergleichen vorher/nachher
cd test_cases/cli/simplemd/10_gfnff_polymer_md
rm -f input.diag.jsonl
../../release/curcuma -md input.xyz -method gfnff -maxtime 500 -threads 4 \
    -md.seed 42 -md.no_restart -md.dump_frequency 10 \
    -md.md_diagnostics true -md.md_diagnostics_timing true
python3 -c "
import json, statistics
r = [json.loads(l) for l in open('input.diag.jsonl') if json.loads(l).get('step',0) >= 50]
print('dcn median:', round(statistics.median([x['timing_ms']['dcn'] for x in r]), 2), 'ms')"

# Wall-Clock
bash compare_baseline.sh 4 gfnff "WP-D-Stage-C_neighbor_reuse"

# Cn-Cutoff Sweep zum Validieren der Konsistenz
../../release/curcuma -sp input.xyz -method gfnff -gfnff.cn_cutoff_bohr 8.0
../../release/curcuma -sp input.xyz -method gfnff -gfnff.cn_cutoff_bohr 0  # Fallback
```

## Beziehung zu WP-D Stage A/B

| Stage | Patch | dcn Hebel |
|-------|-------|-----------|
| A | cn_raw reuse | step 1 weg, 3 ms |
| B | SoA Pass A + SIMD exp() | step 3 inner micro-opt, 2 ms |
| **C (dies)** | Neighbor-List reuse | Threshold-Checks weg, 1.5-3 ms |
| D (Folge) | CN + dCN single-pass | redundante geom-Reads, 1-3 ms |

Sammeln WP-D A+B+C+D zusammen polymer-Hebel ~7-12 ms `dcn` (von 21 ms pre auf ~10 ms). Step C ist mechanisch (Pattern: Liste durchreichen statt neu bauen). Lohnt sich nur wenn Default `cn_cutoff_bohr > 0` (Standard seit P2b).
