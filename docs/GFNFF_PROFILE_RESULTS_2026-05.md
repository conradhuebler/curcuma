# GFN-FF MD Per-Step Profile (May 2026, WP-P1)

**Status**: 🤖 AI-generated initial dataset — ⚙️ machine-tested, single-system,
human production validation pending.

WP-P1 instrumentation (PARAM `md_diagnostics_timing=true`) liefert per-Step-
Wall-Clock-Breakdown ins `<basename>.diag.jsonl`. Erste Mess-Daten von
polymer-NVE-MD zeigen, dass unser bisheriges mentales Kostenmodell der
Hot-Loop teilweise falsch war.

---

## Mess-Setup

| Item | Wert |
|------|------|
| System | `polymer.xyz` (N=1410, nfrag=1) |
| MD | NVE, dt=1.0 fs, max_time=20 fs, dump_frequency=5 |
| Threads | 4 |
| Compiler | g++ (release/) |
| Hardware | AMD Ryzen 9 9950X3D, NVIDIA RTX 5080 (CPU-Pfad) |
| Datenquelle | `polymer.diag.jsonl`, Record 2 (Step 5, nicht der Init-Step) |

Reproduzierbar via:

```bash
echo '{"simplemd":{"method":"gfnff","max_time":20,"dt":1.0,
  "dump_frequency":5,"md_diagnostics":true,"md_diagnostics_timing":true,
  "thermostat":"none","seed":42}}' > cfg.json
./curcuma -md polymer.xyz -import_config cfg.json -threads 4
```

---

## Per-Step Profile (Step 5 = 5 fs)

```
step_total       91.109 ms   <-- gesamter step()
|-- integrator   88.503 ms   <-- Velocity-Verlet inkl. Energy-Call
|   '-- ff_total 88.482 ms   <-- FF-Berechnung (Energy+Gradient)
|       '-- prep_total 68.145 ms   <-- prepareCNAndEEQ (Phase 1+2)
|           |-- dcn         23.322 ms (#1 Bottleneck)
|           |-- eeq_solve   22.224 ms (#2)
|           |-- cn          11.442 ms (#3)
|           |-- d4_gw        3.421 ms (#5)
|           |-- eeq_topo     1.291 ms
|           |-- charge_dist  0.009 ms
|           '-- cnf          0.001 ms
'-- Rest          2.6 ms     <-- AverageQuantities + Restart-Trigger + Loop-Overhead

FF-Terme (Energie/Gradient ohne Prep) = ff_total - prep_total = 20.3 ms (#4)
```

---

## Verifizierte / falsifizierte Hypothesen

### ✓ Bestätigt

1. **WP-S1 `static_charges` ≈ 30 ms/Step Ersparnis** — realer Phase-2-EEQ-Anteil
   ist 22.2 ms. Knapp unter der ursprünglichen Schätzung, aber im Bereich.

2. **WP-S1 `static_cn` Wirkung deutlich höher als angenommen** — geplant ~25 ms
   (cn+dcn+D4-weights). Realer Anteil: cn 11.4 + dcn 23.3 + d4_gw 3.4 =
   **38.2 ms**. `dcn`-Recompute dominiert, war bisher als ~10 ms eingeschätzt.

3. **WP-S1 `static_all` Sparpotential 60.3 ms/Step (66 % von step_total)**
   passt zum gemessenen 3.2× Speedup im Quick-Test (polymer NVT 100 fs:
   10.76 s → 3.34 s).

### ✗ Widerlegt / korrigiert

4. **Roadmap-P4c "75 µs/Bond" trifft hier nicht zu**.
   - Erwartung aus Mai-2026-mixture-Profil: 75 µs × ~3000 polymer-Bonds = ~225 ms
     Bond-Anteil.
   - Realität: gesamte FF-Term-Berechnung (Bond + Angle + Dihedral + Coulomb +
     Dispersion + HB/XB + Repulsion) = **20.3 ms**.
   - Interpretation: Mai-Zahl war wahrscheinlich Memory-Pressure-getrieben bei
     N=6200 (Cache-Spillover) und gilt nicht für N=1410. **WP-P3-Microbench**
     muss das isolieren.

5. **Velocity-Verlet ist praktisch frei** — `integrator − ff_total ≈ 0.02 ms`.
   Optimierung des Verlet-Loops würde nichts bringen.

6. **`hbxb_update`-Trennung nicht möglich** — Wert ist 0 im Output, weil die
   HB/XB-Re-Detection in `Calculation()` verschachtelt ist und durch
   `eeq_topo` (1.3 ms) bereits mit-erfasst wird. Kein eigener Optimierungs-
   punkt.

### ⊘ Nicht ablesbar / offene Fragen

7. **GPU-Profil**: nicht gemessen in dieser Session (CPU-only Smoke-Test).
   Per `setRecordKernelTimings(true)` aktivierbar; `gpu`-Block im JSONL
   wäre dann gefüllt.

8. **Single-Term-Static-Effekt** (statt `static_all`): vorhanden in
   früheren Tests, aber Anteil unter Mess-Rauschen (~0.1 s/100 Steps).
   Erklärung mit obigen Zahlen: `static_cn` allein bringt 38.2 ms/Step ≈
   3.8 s/100 Steps — sollte messbar sein, war es aber nicht im Quick-Test.
   Mögliche Ursachen:
   - Capture-Phase im ersten Step kostet doppelt
   - `dcn`-Skip nur greift, wenn `m_static_state_captured` schon true ist
   - Topology-Reset triggert öfter als gedacht (alle 1-2 Steps?)
   → WP-P2-Sweep mit threads={1,4,8} + Capture-Step-Tracking würde klären.

---

## Konsequenzen für Folge-Optimierungen

### Neue Bottleneck-Reihenfolge (Polymer NVE, CPU)

| Rank | Phase | Zeit | Anteil | Optimierungs-Hebel |
|------|-------|------|--------|---------------------|
| 1 | `dcn` (CN-Derivative pair list) | 23.3 ms | 26 % | SIMD-Vektorisierung, Pair-Cache (siehe unten) |
| 2 | `eeq_solve` (Phase 2 EEQ + Cholesky) | 22.2 ms | 24 % | Bereits WP-S1 (`static_charges`) und WP-S3 (`eeq_distance_cutoff_auto`); EEQ-Pool-Skalierung |
| 3 | `cn` (Coordination Numbers) | 11.4 ms | 13 % | SIMD `erf()` (P4e), OpenMP-Skalierung-Audit |
| 4 | FF-Terme zusammen | 20.3 ms | 22 % | Aktuell genauer aufschlüsseln (Bond vs Angle vs Coulomb) — WP-P3-Microbench |
| 5 | `d4_gw` (D4 Gaussian weights) | 3.4 ms | 4 % | P4f SIMD-Vektorisierung |

### Konkrete Aktions-Vorschläge

- **Neuer WP-Kandidat**: **dcn-Vektorisierung** in `calculateCoordinationNumberDerivatives()`.
  - 26 % des step_total — höher als der dokumentierte Bond-Hotspot.
  - Pair-list ist O(N²)-erf()-Loop — gleiche SIMD-Lücke wie CN selbst.
  - Erwartet: 2-3× Speedup → 8-12 ms gespart.

- **WP-P3 Mikrobench (Bond-Term)** wird durch diese Daten **weniger dringend**.
  - Bond-Anteil ist im Polymer-Profil ein Teil der 20.3 ms FF-Terme, nicht
    die 225 ms wie aus Mai-Profil extrapoliert.
  - Microbench bleibt als saubere Klärung sinnvoll, ist aber kein
    high-priority-Target mehr.

- **WP-P2 Cross-System-Sweep** wird **wichtiger**: bestätigt oder widerlegt,
  ob die polymer-Befunde auf mixture (N=6200) skalieren oder ob das Mai-2026-
  Profil eine andere Bottleneck-Topographie zeigt.

- **WP-S1 `static_all`-Default-Empfehlung**: für equilibrium NVT-Production-
  MD auf neutralen organischen Systemen ist `static_all=true` jetzt durch
  Profil-Daten gerechtfertigt (3.2× speedup, 66 % saubere Phase-Skip-Decke).

- **Folge-WP für `dcn`**: separater Aufwand, ~1 Tag. Skizze:
  1. Profilieren der `calculateCoordinationNumberDerivatives()`-Inner-Loop
  2. Pair-list ist bereits da (`m_cn_deriv_pair_list`?), aber erf()-Calls skalar
  3. Eigen-Array-Vektorisierung der `dr/r * dcn_factor`-Operation
  4. Test gegen finite-difference (`test_gfnff_gradients`)

---

## Anwendung des WP-P1-Hooks (für Nachnutzer)

```bash
# Profil-Daten generieren
echo '{"simplemd":{"method":"gfnff","max_time":100,"dt":1.0,
       "md_diagnostics":true,"md_diagnostics_timing":true}}' > cfg.json
./curcuma -md mol.xyz -import_config cfg.json

# Schnellanalyse
python3 scripts/analyse_diag.py mol.diag.jsonl
# Oder direkt: jq für Phasen-Aggregation
python3 -c "
import json
recs = [json.loads(l) for l in open('mol.diag.jsonl')]
# Skip Record 0 (Init), aggregiere Rest
ts = [r['timing_ms'] for r in recs[1:]]
avg = {k: sum(t[k] for t in ts)/len(ts) for k in ts[0] if isinstance(ts[0][k], (int,float))}
for k, v in sorted(avg.items(), key=lambda kv: -kv[1]):
    print(f'{k:14s} {v:>8.2f} ms')
"
```

---

## Caveats

- Nur **ein System** (polymer N=1410), **ein Zustand** (NVE, dt=1.0 fs),
  **ein Step** (Step 5). Über mehrere Steps mitteln und Cross-System-Sweep
  liefert WP-P2.
- **NVE statt NVT**: bei csvr-Thermostat würde der Thermostat-Hook hinzukommen
  (~0.1-0.5 ms), für reine Profile-Analyse aber irrelevant.
- **Capture-Phase im ersten Step** ist nicht separat ausgewiesen — Record 0
  (Step 0) hat mit Init-Cost andere Werte.
- **Hardware-spezifisch**: AMD Zen 4 mit 4 Threads. Auf anderen Maschinen
  können Skalierungen anders aussehen.

---

## Links

- WP-P1 Plan: [GFNFF_PROFILE_WP1_TIMING_INSTRUMENTATION.md](GFNFF_PROFILE_WP1_TIMING_INSTRUMENTATION.md)
- WP-P2 Plan (Cross-System-Sweep): [GFNFF_PROFILE_WP2_BENCHMARK_SWEEP.md](GFNFF_PROFILE_WP2_BENCHMARK_SWEEP.md)
- WP-P3 Plan (Bond-Mikrobench): [GFNFF_PROFILE_WP3_BOND_HOTSPOT.md](GFNFF_PROFILE_WP3_BOND_HOTSPOT.md)
- Performance-Roadmap (Mai 2026): [GFNFF-PERFORMANCE-ROADMAP.md](GFNFF-PERFORMANCE-ROADMAP.md)
- WP-S1 Static-Mode: [GFNFF_STATIC_WP1_FROZEN_STATE.md](GFNFF_STATIC_WP1_FROZEN_STATE.md)

---

## WP-P2 Cross-System Sweep (May 2026)

**Setup**: `scripts/wp_profile_sweep.py`, MD csvr 200 fs / dt=1 fs /
dump=10, seed=42, threads=4, AMD Ryzen 9 9950X3D, release/.
Mittelwert über Records 1–20 pro Konfig (Record 0 = Init-Step geskippt).

**Geskippte Configs**:
- `mixture` (N=6200): MD instabil unter dt=1 fs/csvr — skip per User-Hinweis.
  Mai-2026-Mixture-SP-Profil bleibt der einzige Datenpunkt für N>>1500.
- `gpu`-Pfad: USE_CUDA=OFF auf der Test-Maschine (Standard-release/-Build).

### caffeine (N=24, CPU)

| Modus | step_total | dcn | eeq_solve | cn | d4_gw | ff_total | speedup |
|-------|------------|-----|-----------|-----|------|----------|---------|
| baseline    | 0.38 | 0.01 | 0.08 | 0.01 | 0.01 | 0.27 | 1.00× |
| static_all  | 0.15 | 0.00 | 0.00 | 0.00 | 0.00 | 0.08 | 2.55× |
| cutoff_auto | 0.29 | 0.01 | 0.06 | 0.01 | 0.01 | 0.20 | 1.29× |
| both        | 0.15 | 0.00 | 0.00 | 0.00 | 0.00 | 0.07 | 2.49× |

### polymer (N=1410, CPU)

| Modus | step_total | dcn | eeq_solve | cn | d4_gw | ff_total | speedup |
|-------|------------|-----|-----------|-----|------|----------|---------|
| baseline    | 102.09 | 20.96 | 21.63 | 10.66 | 3.54 | 99.26 | 1.00× |
| static_all  |  28.58 |  0.00 |  0.00 |  0.00 | 0.00 | 25.84 | **3.57×** |
| cutoff_auto |  99.02 | 21.03 | 20.24 | 10.69 | 3.46 | 96.27 | 1.03× |
| both        |  28.70 |  0.00 |  0.00 |  0.00 | 0.00 | 25.91 | 3.56× |

---

## Auswertung gegen die Bottleneck-Lücken (siehe oben "Nicht ablesbar")

### ✓ Bestätigt

1. **WP-P1-Polymer-Snapshot reproduziert** über 20 Records:
   `dcn` 21.0 ms ≈ Snapshot 23.3, `eeq_solve` 21.6 ≈ 22.2, `cn` 10.7 ≈ 11.4.
   Single-Step-Wert war nicht zufällig.

2. **`static_all` ist der Workhorse**:
   - Caffeine 2.5× (klein, FF-Terme dominieren auch im baseline)
   - Polymer 3.6× (große System-Klasse, Prep dominiert)
   - Bestätigt die geometrische Skala-Abhängigkeit: je größer N, desto
     mehr lohnt sich Static-Mode.

3. **Bond-Term-Anteil bei Polymer**:
   `ff_total` (99.3 ms) − `total`-Prep (≈ 56 ms aus dcn+eeq+cn+d4_gw+...)
   ≈ **43 ms** für alle FF-Energie-Terme (Bond+Angle+Dihedral+Coulomb+
   Disp+HB+Repulsion) zusammen. Bei ~3000 polymer-Bonds wäre der
   Mai-2026-Hypothese-Wert von 75 µs/Bond = 225 ms — wir messen **43 ms total
   für alle Bonded+Nonbonded-Terme**. P4c-Hypothese ist für polymer-Größe
   eindeutig **falsifiziert**. Mai-Befund stammt sehr wahrscheinlich aus
   N=6200-Memory-Pressure und gilt nicht für mittelgroße Systeme.

### ✗ Falsifiziert / korrigiert

4. **`cutoff_auto` allein bringt fast nichts**:
   - Caffeine: 1.29× (durch N<200-Gate in eeq_solver.cpp nicht aktiv;
     Speedup vermutlich Rausch- oder Topology-Cache-Effekt)
   - Polymer: **1.03×** (effektiv null trotz aktivem Cutoff bei N>200)
   - Erklärung: bei polymer-Größe (1410 Atome) ist `eeq_solve` nur
     21.6 ms, der Cutoff-Anteil davon noch kleiner. Erst bei mixture-Größen
     (N=6200, Mai-Profil 660 ms eeq_solve) lohnt sich Cutoff-Sparsifizierung.
   - **Konsequenz**: WP-S3 (`eeq_distance_cutoff_auto`) sollte als
     "großes-System-Feature" dokumentiert werden — kein Default-Win für
     mittelgroße Polymere.

5. **`static_all` + `cutoff_auto` = `static_all`**:
   - Polymer `both` 28.70 ms vs `static_all` 28.58 ms (Differenz im Rausch).
   - Logisch: bei aktivem Static-Mode wird Phase 2 EEQ skippt, der Cutoff
     hat keinen Angriffspunkt mehr.

### ⊘ Offen (nicht aus diesem Sweep ablesbar)

6. **GPU-Pfad**: nicht gemessen (USE_CUDA=OFF auf Test-Maschine).
7. **WP-S3 + GPU**: per Konstruktion CPU-only (`eeq_solver.cpp:1080`
   greift im CPU-EEQSolver-Pfad). GPU hat eigene `k_coulomb`-Cutoff-Logik.
8. **`static_cn` allein vs `static_charges` allein**: nicht im Sweep
   ausgeführt. Mit dem 0%-Speedup-Befund aus dem ursprünglichen Quick-Test
   (28 ms Polymer-Speedup-Lücke) bleibt ungeklärt, ob der erste Step
   die Capture-Phase teuer macht und so den Mean verzerrt.
9. **Threading-Audit**: nicht im Default-Sweep, aber via `--thread-sweep`
   ausführbar.

---

## Neue Optimierungs-Empfehlungen (Stand Mai 2026)

| Rang | Maßnahme | Hebel | Aufwand |
|------|----------|-------|---------|
| 1 | **`dcn`-Vektorisierung** (`calculateCoordinationNumberDerivatives()`) | Polymer 21 ms Phase, eigenes ~2-3× erwartet → 7-14 ms Save | ~1 Tag |
| 2 | **`cn`-SIMD-Loop** (`std::erf()` vektorisieren, P4e) | Polymer 10.7 ms Phase, ~2× erwartet → 5 ms Save | ~0.5 Tag |
| 3 | **GPU-Pfad messen + optimieren** (USE_CUDA-Build aktivieren) | Unbekannt — Roadmap-Mai sagt 7.8× GPU vs CPU bei polymer; Static-Mode-Effekt auf GPU nicht verifiziert | ~1 Tag inkl. Mess-Run |
| 4 | **`eeq_solve` parallelisieren bei großen Systemen** | Bei polymer 21.6 ms ≈ 80% serial; bei mixture entscheidend (660 ms aus Mai-Profil) | ~1.5 Tage |
| 5 | **Bond-Term-Mikrobench WP-P3** | Niedrig priorisiert: 43 ms aller FF-Terme zusammen bei polymer, nicht Bond allein | ~0.5 Tag |

WP-S3 (`eeq_distance_cutoff_auto`) bleibt sinnvoll für N>4000-Systeme,
aber **nicht** als Default für mittelgroße MD-Production-Läufe.

---

## Reproduzierbarkeit

```bash
cd /home/conrad/src/curcuma
python3 scripts/wp_profile_sweep.py \
    --curcuma release/curcuma \
    --systems caffeine,polymer --paths cpu \
    --max-time 200 --timeout 900 \
    --output /tmp/sweep_cpu_full
cat /tmp/sweep_cpu_full/wp_profile_summary.md
```

Erweitert zu `caffeine,polymer,mixture` falls mixture-Stabilität
gelöst (siehe `docs/GFNFF_STATIC_WP4_VALIDATION_SUITE.md` für
Re-Capture-Idee).
