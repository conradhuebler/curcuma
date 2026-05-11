# WP-P1: GFN-FF Per-Step Timing Instrumentation (CPU + GPU)

**Kategorie**: Klein–Mittel
**Aufwand**: ~1 Tag (~120 LoC)
**Wirkung**: Per-Step-Breakdown der MD-Hot-Loop in `<basename>.diag.jsonl`; macht Bottleneck-Analyse messbar statt heuristisch
**Abhängigkeiten**: WP-S2 (Diagnostics-JSONL); CPU-seitige `PrepTiming` existiert bereits
**Status**: ✅ Implementiert (May 2026) — ⚙️ machine-tested, single-system smoke; first
results in [GFNFF_PROFILE_RESULTS_2026-05.md](GFNFF_PROFILE_RESULTS_2026-05.md)
**Priorität**: Voraussetzung für WP-P2 (Benchmark-Sweep) und WP-P3 (Bond-Microbench)

---

## Motivation

WP-S1/S2/S3 wurden auf Heuristik-Basis gebaut. Der "schnelle Lauf"-Test
(polymer NVT 100 fs) zeigte: nur `static_all + cutoff_auto` liefert Speedup
(3.2×); einzelne Flags hatten 0 % messbaren Effekt. Das mentale Kostenmodell
ist offenbar teilweise falsch. Voraussetzung für saubere Optimierung: harte
Zahlen pro Phase pro Step.

`PrepTiming`-Struct existiert bereits in `gfnff_method.cpp`
(cn/eeq_topo/cnf/dcn/d4_gw/eeq_solve/charge_dist mit ms-Präzision, gedruckt
bei Verbosity ≥ 2). GPU-Stream-Timings ebenfalls (`setRecordKernelTimings`
mit Per-Stream-CUDA-Events). Beide werden bisher nur als One-Shot-Log
ausgegeben, nicht in die MD-JSONL-Trajektorie geschrieben.

---

## Output-Erweiterung

In jedem WP-S2-Snapshot (`<basename>.diag.jsonl`) zusätzliches Feld:

```json
{
  "step": 50, "time_fs": 50.0,
  "energy": { ... },
  "charges": [ ... ],
  "timing_ms": {
    "prep_cn":        2.1,
    "prep_eeq_topo":  0.4,
    "prep_dcn":       8.7,
    "prep_d4_gw":     9.3,
    "prep_eeq_solve": 28.4,
    "ff_total":       42.6,
    "integrator":     5.2,
    "hbxb_update":    0.0,
    "write_geo":      0.0,
    "step_total":     97.8,
    "gpu": {
      "ms_sA":  6.1,
      "ms_sB":  3.4,
      "ms_sC":  4.2,
      "ms_p2": 12.7,
      "total": 26.4
    }
  }
}
```

GPU-Block nur bei aktivem GPU-Pfad. Default des neuen PARAM ist `false`,
weil Hook-Overhead nicht-null (~1-2 µs pro `chrono::high_resolution_clock`-Call).

---

## Implementierungspunkte

### 1. Neuer PARAM in `src/capabilities/simplemd.h`

Im `BEGIN_PARAMETER_DEFINITION(simplemd)`-Block neben `md_diagnostics`:

```cpp
PARAM(md_diagnostics_timing, Bool, false,
      "Add a 'timing_ms' block to each <basename>.diag.jsonl record (per-phase wall-clock breakdown: CN/EEQ/dcn/D4-weights/FF/integrator/HBXB/I-O). For GPU runs adds a 'gpu' sub-block with per-stream CUDA-event timings. Requires md_diagnostics=true. ~1-2 us overhead per hook.", "Output", {})
```

### 2. CPU: `PrepTiming`-Exposure auf `GFNFF`

`GFNFF` cached die letzte `PrepTiming` als Public-Member, vergleichbar zu
`m_eeq_cutoff_auto_active`:

```cpp
// gfnff.h, neben WP-S1/S3-Diagnostic-Members
mutable PrepTiming m_last_prep_timing{};  ///< WP-P1: cached after each prepareCNAndEEQ()

const PrepTiming& getLastPrepTiming() const { return m_last_prep_timing; }
```

In `prepareCNAndEEQ()` (ab `gfnff_method.cpp:957`): am Funktionsende
unbedingt `m_last_prep_timing = local_timing;`, auch wenn `out_timing == nullptr`.

### 3. CPU: Forwarder auf `ComputationalMethod` + `EnergyCalculator`

`computational_method.h` — zwei neue virtuelle Methoden:

```cpp
/// WP-P1: latest CPU phase breakdown (FF method only; default empty json)
virtual json getLastPrepTiming() const { return {}; }

/// WP-P1: latest GPU stream timings (GFN-FF-GPU only; default empty json)
virtual json getStreamTimings() const { return {}; }
```

`gfnff_method.h/.cpp` (Override → `m_gfnff->getLastPrepTiming()` als JSON).
`gfnff_gpu_method.h/.cpp` (Override → `m_gpu_workspace->getLastStreamTimings()`).

`EnergyCalculator` bekommt Forwarder `LastPrepTiming()` und `StreamTimings()`.

### 4. GPU: `FFWorkspaceGPU::getLastStreamTimings()` API

`cuda/ff_workspace_gpu.cu` hat schon `m_record_kernel_timings` und
elapsed-Time-Berechnung in `launchChargeDependentAndFinish()`. Neuer Member:

```cpp
struct StreamTimings { double ms_sA, ms_sB, ms_sC, ms_p2, total; };

private:
StreamTimings m_last_stream_timings{};  ///< WP-P1: cached after each calculate()

public:
const StreamTimings& getLastStreamTimings() const { return m_last_stream_timings; }
```

Ausgabepfad-Block (`cuda/ff_workspace_gpu.cu:2480-2498`): statt nur Log
schreiben, auch `m_last_stream_timings` setzen.

### 5. SimpleMD: Sub-Phase-Wall-Clocks in `step()`

Im `step()`-Body (`simplemd.cpp` ab Z. 1720) mit `std::chrono::high_resolution_clock`:

```cpp
auto t0 = std::chrono::high_resolution_clock::now();
// ... Integrator() ...
auto t_after_integrator = std::chrono::high_resolution_clock::now();
double integrator_ms = std::chrono::duration<double, std::milli>(t_after_integrator - t0).count();

// ... Energy()/CalculateEnergy() ...
auto t_after_energy = ...;
double ff_total_ms = ...;

// ... HBXB-Update läuft in Energy() — Sub-Timing kommt aus PrepTiming.eeq_topo + dedicated tracker ...

// ... WriteGeometry() ...
double write_geo_ms = ...;
```

Im JSONL-Snapshot-Block (`simplemd.cpp:1785-1800`) sammelt SimpleMD die
Wall-Clocks und übergibt sie an den Writer:

```cpp
if (m_md_diagnostics && m_diag_writer) {
    json timing;
    if (m_md_diagnostics_timing) {
        timing = m_interface->LastPrepTiming();
        timing["ff_total"]   = ff_total_ms;
        timing["integrator"] = integrator_ms;
        timing["hbxb_update"]= hbxb_update_ms;
        timing["write_geo"]  = write_geo_ms;
        timing["step_total"] = step_total_ms;
        json gpu_t = m_interface->StreamTimings();
        if (!gpu_t.empty()) timing["gpu"] = gpu_t;
    }
    m_diag_writer->writeSnapshot(..., timing /* optional argument */);
}
```

### 6. `MDDiagnosticsWriter` erweitern

Neue Overload mit optionalem Timing-Block:

```cpp
void writeSnapshot(int step, double time_fs,
                   const json& energy_decomp,
                   const Vector& charges,
                   const Vector& cn,
                   const Matrix& gradient,
                   int hb_count, int xb_count,
                   const json& timing = json{});
```

Backward-compatible: default-Argument `json{}` wird nur geschrieben, wenn
non-empty.

---

## Kritische Files

| Datei | Änderung | LoC |
|-------|----------|-----|
| `src/capabilities/simplemd.h` | 1 PARAM + 1 Member | ~6 |
| `src/capabilities/simplemd.cpp` | Sub-phase Wall-Clocks + JSONL-Block-Build + Param-Read | ~40 |
| `src/capabilities/md_diagnostics.h/.cpp` | optionaler timing-Parameter in writeSnapshot | ~10 |
| `src/core/energy_calculators/computational_method.h` | 2 virtuals mit default empty json | ~8 |
| `src/core/energy_calculators/qm_methods/gfnff_method.h/.cpp` | Override CPU | ~12 |
| `src/core/energy_calculators/qm_methods/gfnff_gpu_method.h/.cpp` | Override GPU | ~12 |
| `src/core/energy_calculators/ff_methods/gfnff.h/.cpp` | `m_last_prep_timing`-Cache + Getter | ~10 |
| `src/core/energy_calculators/ff_methods/cuda/ff_workspace_gpu.h/.cu` | `StreamTimings`-Struct + Cache + Getter | ~20 |
| `src/core/energycalculator.h/.cpp` | 2 Forwarder | ~12 |

Build-Verify: `cd release && make GenerateParams && make -j4`. Erwartet:
469 PARAMs (+1 ggü WP-S3).

---

## Verifikation

### a) Build & PARAM-Registrierung
```bash
cd /home/conrad/src/curcuma/release
MKL_ROOT=/opt/intel/oneapi/mkl/latest make GenerateParams && make -j4
grep -A 3 '"md_diagnostics_timing"' generated/parameter_registry.h
```

### b) Smoke-Test CPU
```bash
echo '{"simplemd": {"method":"gfnff","max_time":20,"dt":1.0,
"dump_frequency":5,"md_diagnostics":true,"md_diagnostics_timing":true}}' \
   > /tmp/wpP1.json
./curcuma -md ../test_cases/molecules/larger/caffeine.xyz -import_config /tmp/wpP1.json
python3 -c "import json; rec=json.loads(open('caffeine.diag.jsonl').readline()); \
            print(json.dumps(rec.get('timing_ms', {}), indent=2))"
```
Akzeptanz: `timing_ms`-Block mit ≥ 8 numerischen Feldern, alle ≥ 0.
`step_total ≥ ff_total + integrator + write_geo`.

### c) Smoke-Test GPU (USE_CUDA)
```bash
./curcuma -md polymer.xyz -method gfnff -gpu cuda -maxtime 20 -dt 1.0 \
  -dump 5 -md_diagnostics true -md_diagnostics_timing true
python3 -c "import json; rec=json.loads(open('polymer.diag.jsonl').readline()); \
            print(json.dumps(rec['timing_ms'].get('gpu', {}), indent=2))"
```
Akzeptanz: `gpu`-Subblock mit `ms_sA/sB/sC/p2/total`, alle ≥ 0.

### d) Default-Off / kein Codepfad-Wechsel
```bash
./curcuma -md caffeine.xyz -method gfnff -maxtime 20 -dt 1.0
ctest --output-on-failure
```
Akzeptanz: ohne `md_diagnostics_timing=true` keine JSONL-Erweiterung,
Default-Pfad unverändert. Regression bleibt grün.

---

## Out of Scope (für WP-P2/P3)

- Benchmark-Sweep über mehrere Systeme/Modi → WP-P2
- Bond-Mikrobench im Microsecond-Bereich → WP-P3
- Aggregation/Reporting der gesammelten Daten → WP-P2

---

## Risiken

| Risiko | Mitigation |
|--------|-----------|
| Timing-Hooks brechen async GPU-Pipeline | Wall-Clock-Hooks nur CPU-seitig; GPU nutzt die bereits stabile `setRecordKernelTimings`-Infrastruktur |
| Overhead bei `md_diagnostics_timing=true` macht das Ergebnis ungenau | ~1-2 µs/Hook × 8 Hooks = ~15 µs/Step, vernachlässigbar gegenüber 50-100 ms/Step polymer. Doku-Hinweis im PARAM. |
| Bestehende Tests brechen | Default `false` → kein Codepfad-Wechsel |
