# WP-P2: GFN-FF Cross-System Benchmark-Sweep + Profile Report

**Kategorie**: Klein–Mittel (Skript-lastig)
**Aufwand**: ~1 Tag (~250 LoC inkl. Skript + Doku)
**Wirkung**: 24-Konfig Mess-Matrix mit Pro-Phase-Breakdown; klärt welche WP-S-Flags wirklich Speedup bringen und wo die wirklichen Bottlenecks sitzen
**Abhängigkeiten**: WP-P1 (Timing-Hooks im JSONL); WP-S1/S2/S3 (Static-Mode, Diagnostics, cutoff_auto)
**Status**: 🤖 Geplant
**Priorität**: nach WP-P1, vor WP-P3 (Microbench-Targets ergeben sich aus den Sweep-Ergebnissen)

---

## Motivation

Bisher haben wir genau zwei Performance-Messpunkte: polymer NVT 100 fs
(10.76 s baseline → 3.34 s static_all+cutoff_auto, 3.2×) und caffeine NVT
500 fs (0.17 s → 0.07 s, 2.3×). Einzelne Flags zeigten 0% Speedup, was nicht
zur Theorie passt. Ohne strukturierten Sweep über (System × Pfad × Mode)
bleibt die WP-S-Wirksamkeitsfrage qualitativ.

WP-P1 liefert die per-Phase-Wall-Clocks; WP-P2 ist der Treiber, der sie
über eine 24-Konfig-Matrix erhebt und in einem Ergebnis-Bericht
aggregiert.

---

## Sweep-Matrix

3 Systeme × 2 Pfade × 4 Modi = **24 Konfigurationen**.

| System | N | Charakter |
|--------|---|-----------|
| caffeine | 24 | klein, rigid, isolated |
| polymer | 1410 | mittelgroß, hauptanwendung, lange Kette |
| mixture | 6200 | sehr groß, nfrag=1400 (Test verfügbar?) |

| Pfad | Auslöser |
|------|----------|
| CPU | `-method gfnff -threads 4` |
| GPU | `-method gfnff -gpu cuda -threads 4` (nur falls USE_CUDA) |

| Modus | Flags |
|-------|-------|
| `baseline` | (keine) |
| `static_all` | `-static_all true` |
| `cutoff_auto` | `-eeq_distance_cutoff_auto true` |
| `both` | `-static_all true -eeq_distance_cutoff_auto true` |

Jeder Lauf: 200 fs NVT (csvr, dt=1.0 fs), `dump_frequency=10`,
`md_diagnostics=true`, `md_diagnostics_timing=true` (WP-P1), `seed=42`.

200 fs / 10 = 20 Records pro Lauf; insgesamt ~480 JSONL-Records über alle
Configs. Pro Polymer-Lauf ~22 s, pro Mixture-Lauf evtl. ~10 min — Sweep
sollte unter 1 h durchlaufen.

---

## Implementierungspunkte

### 1. Sweep-Skript `scripts/wp_profile_sweep.py`

Python 3 ohne externe Dependencies (`json`, `subprocess`, `argparse`,
`statistics`). Struktur:

```python
import json, subprocess, argparse, statistics, csv, time
from pathlib import Path

SYSTEMS = {
    "caffeine": "test_cases/molecules/larger/caffeine.xyz",
    "polymer":  "test_cases/molecules/larger/polymer.xyz",
    "mixture":  "test_cases/molecules/larger/mixture.xyz",  # optional, check upfront
}
PATHS = ["cpu", "gpu"]
MODES = {
    "baseline":    {},
    "static_all":  {"static_all": True},
    "cutoff_auto": {"eeq_distance_cutoff_auto": True},
    "both":        {"static_all": True, "eeq_distance_cutoff_auto": True},
}

def run_one(system_path, path, mode_flags, out_dir, max_time=200, dt=1.0):
    cfg = {
        "gfnff": dict(mode_flags),
        "simplemd": {
            "method": "gfnff", "max_time": max_time, "dt": dt,
            "thermostat": "csvr", "dump_frequency": 10,
            "md_diagnostics": True, "md_diagnostics_timing": True,
            "seed": 42, "threads": 4,
        },
        "global": {"method": "gfnff", "threads": 4, "gpu": (path == "gpu")},
    }
    cfg_path = out_dir / "cfg.json"
    cfg_path.write_text(json.dumps(cfg))
    t0 = time.monotonic()
    cp = subprocess.run([CURCUMA, "-md", system_path,
                         "-import_config", str(cfg_path)],
                        cwd=out_dir, capture_output=True, text=True)
    wall = time.monotonic() - t0
    return wall, cp.returncode

def parse_jsonl(jsonl_path):
    records = []
    for line in open(jsonl_path):
        records.append(json.loads(line))
    return records

def aggregate(records):
    if not records: return {}
    fields = list(records[1].get("timing_ms", {}).keys())  # skip step 0 (no full step)
    out = {}
    for f in fields:
        vals = [r["timing_ms"][f] for r in records[1:] if f in r["timing_ms"]]
        if isinstance(vals[0], (int, float)):
            out[f] = statistics.mean(vals)
    return out

def main():
    # iterate, collect, write CSV + summary markdown
    ...
```

Outputs:
- `<out_dir>/wp_profile_results.csv` — alle Records: `system, path, mode, step, time_fs, step_total, prep_cn, prep_eeq_solve, ff_total, integrator, write_geo, ms_p2_gpu, ...`
- `<out_dir>/wp_profile_summary.md` — aggregierte Tabelle:

```markdown
## System: polymer (N=1410)

### CPU path
| Modus | step_total | prep_cn | prep_eeq | ff_total | integrator | speedup |
|-------|-----------|---------|----------|----------|-----------|---------|
| baseline    | 105.4 | 5.1 | 28.7 | 65.3 | 5.2 | 1.0× |
| static_all  |  31.8 | 0.0 |  0.0 | 26.4 | 5.2 | 3.3× |
| cutoff_auto |  94.1 | 5.1 | 18.5 | 65.4 | 5.2 | 1.1× |
| both        |  28.6 | 0.0 |  0.0 | 23.2 | 5.2 | 3.7× |
```

(Werte fiktiv; echte Werte kommen aus dem Sweep.)

### 2. Ergebnis-Bericht `docs/GFNFF_PROFILE_RESULTS_2026-05.md`

Statisches Markdown mit den gesammelten Tabellen, kommentiert:

- 24-Konfig-Übersicht
- Pro-Phase-Heatmap (CPU vs GPU, baseline vs Modi)
- Identifizierte Top-5-Bottlenecks mit Größenordnung
- Antworten auf die 10 Bottleneck-Lücken (siehe Kommentar im Sweep-Skript)
- Empfehlungen für Folge-WPs (welche P4*-Maßnahmen aus
  `GFNFF-PERFORMANCE-ROADMAP.md` versprechen wirklich Speedup)

### 3. Roadmap-Update

`docs/GFNFF-PERFORMANCE-ROADMAP.md` erhält am Ende einen neuen Block:

```markdown
## Mai 2026 — MD-Loop-Profil (WP-P2)

Quelle: `docs/GFNFF_PROFILE_RESULTS_2026-05.md`, 24-Konfig Sweep über
caffeine/polymer/mixture, CPU+GPU, 4 Modi.

Hauptergebnisse (polymer NVT 200 fs, 4 threads, CSVR):
- ...
```

### 4. CI-Integration (optional)

Nightly-Test (`add_test` mit Label `nightly`), der den Sweep für eine
reduzierte Matrix laufen lässt (nur polymer + caffeine, nur 50 fs):

```cmake
add_test(NAME wp_profile_sweep_nightly
    COMMAND python3 ${CMAKE_SOURCE_DIR}/scripts/wp_profile_sweep.py
            --systems caffeine,polymer --max-time 50 --output ${CMAKE_BINARY_DIR}/profile_nightly
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
set_tests_properties(wp_profile_sweep_nightly PROPERTIES
    TIMEOUT 600 LABELS "nightly;profile")
```

Default-CTest greift den nicht (ohne `-L nightly`).

---

## Kritische Files

| Datei | Änderung | LoC |
|-------|----------|-----|
| `scripts/wp_profile_sweep.py` (neu) | Sweep + CSV + Markdown-Aggregation | ~180 |
| `docs/GFNFF_PROFILE_RESULTS_2026-05.md` (neu) | Ergebnisbericht | ~100 |
| `docs/GFNFF-PERFORMANCE-ROADMAP.md` | neuer "Mai 2026 MD-Loop-Profil"-Block | ~30 |
| `test_cases/cli/CMakeLists.txt` | optionales Nightly-Test-Eintrag | ~4 |

Kein C++ — reiner Python + Markdown + CMake. Build entfällt.

---

## Verifikation

### a) Sweep-Smoke-Test (kurze Matrix)
```bash
python3 scripts/wp_profile_sweep.py --systems caffeine --paths cpu --max-time 20 \
        --output /tmp/sweep_smoke
ls /tmp/sweep_smoke/wp_profile_*
cat /tmp/sweep_smoke/wp_profile_summary.md
```
Akzeptanz: 4 Configs ausgeführt (4 Modi × 1 System × 1 Pfad), Summary
enthält Tabelle mit ≥ 4 Zeilen, Speedup-Spalte vorhanden.

### b) Voller Sweep
```bash
python3 scripts/wp_profile_sweep.py --systems caffeine,polymer --paths cpu,gpu \
        --max-time 200 --output /tmp/sweep_full
```
Akzeptanz: 16 Configs (2 Systeme × 2 Pfade × 4 Modi). Bei vorhandenem
mixture.xyz erweitert auf 24 Configs.

### c) Optional Mixture
```bash
[ -f test_cases/molecules/larger/mixture.xyz ] && \
    python3 scripts/wp_profile_sweep.py --systems mixture --paths cpu --max-time 200 \
        --output /tmp/sweep_mixture
```
Akzeptanz: läuft, falls Mixture verfügbar. Skript skippt cleanly mit Warning,
falls nicht.

### d) CSV-Konsistenz
```bash
python3 -c "import csv; rows=list(csv.DictReader(open('/tmp/sweep_full/wp_profile_results.csv'))); \
            print('rows:', len(rows), 'cols:', list(rows[0].keys()))"
```
Akzeptanz: ≥ 4 × 20 = 80 Records bei 200 fs/dump=10/4 Configs.

### e) Roadmap-Update
```bash
grep -c "Mai 2026 — MD-Loop-Profil" docs/GFNFF-PERFORMANCE-ROADMAP.md
```
Akzeptanz: 1 (Block existiert).

---

## Out of Scope (für WP-P3 und Folge-WPs)

- Bond-Term-Mikrobench im µs-Bereich → WP-P3
- Aktuelle Bond-Optimierung (P4c-Fix) → eigenes WP nach WP-P3-Mikrobench-Ergebnis
- SoA-Refactoring der CPU-Geometrie → eigenes großes WP, hängt von WP-P2-Ergebnissen ab
- Vergleich gegen xtb-Fortran als externer Maßstab → eigenes "Cross-Implementation-Bench"-WP

---

## Risiken

| Risiko | Mitigation |
|--------|-----------|
| Mixture-System nicht im Repo | Skript skippt mit Warning, Sweep reduziert auf 16 Configs (caffeine+polymer) |
| Polymer-NVT explodiert während Sweep | timeout pro Run (per `subprocess.run(timeout=300)`), Skip + Logging |
| GPU-Pfad nicht verfügbar (USE_CUDA=OFF) | Pfad-Auswahl in Sweep-Args; CPU-Only-Sweep liefert 12 statt 24 Configs |
| Sweep total > 1 h | `--max-time 50` als schnellerer Smoke-Modus; voller Sweep nur Nightly |
| Statistische Variation einzelner Runs hoch | Sweep nimmt Mittelwert über alle Records (außer Step 0), nicht nur Total-Wall |
