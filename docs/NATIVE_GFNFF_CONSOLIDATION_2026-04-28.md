# Native-GFN-FF Branch-Konsolidierung — 28. April 2026

## Kontext

Die zwei Linien `native-gfnff` (lokal) und `origin/native-gfnff` (remote) waren um 12 / 17 Commits divergent (gemeinsamer Vorfahre `e73246e`). Lokal hatte zusätzlich Cherry-Picks aus `feature/gfnff-gpu-optimizations` (Apr 24-25), darunter den HB-Re-Detection-Fix (`a0ab407`) und EEQ-Robustheits-Patches (`ada6925`, `7389e2e`). Origin hatte im selben Zeitfenster eine eigenständige Performance-Reihe (P1a–P6, G1a–G3, G2a) plus den CPU-GFN-FF analytical-vs-numerical-Gradient-Test (`7c36da2`) und PBC minimum-image (`f8c6aa0`).

Auslöser: Frage nach Triage der GPU-Performance-Commits gegen den GFN-FF-Gradienten ("welche Commits brechen den Gradienten, welche sind sinnvoll").

## Was passiert ist

### 1. Stash-Aufräumen
`stash@{0}` enthielt einen abgebrochenen `e37db8b`-Cherry-Pick-Versuch mit un-aufgelösten Konflikt-Markern. Inhalt nach `/tmp/native-gfnff-stash0.patch` gesichert, dann verworfen.

### 2. Merge `origin/native-gfnff` → lokal `native-gfnff` (`f4d0c5f`)
Drei Konfliktdateien manuell aufgelöst (alle wertvollen Blöcke beider Seiten erhalten):

- `src/core/energy_calculators/qm_methods/gfnff_gpu_method.cpp` — HEAD's HB-Re-Detection-Block (`rebuildBondHBData()`, `updateHBAlphaPairs()`, `updateBondHBMetadata()`) erhalten **plus** origin's PBC-Unit-Cell-Upload (`f8c6aa0`) und D4 Gaussian-Weights-Skip-Threshold (`d75048e`, P1a) automatisch gemergt.
- `src/capabilities/simplemd.cpp` — HEAD's Stepwise-API (`step()` / `prepareRun()` / `finalizeRun()`-Refactor, `30704c5`) erhalten; origin's verbessertes DOF-Logging (`InitConstrainedBonds()`, "1-2, 1-3" Aufschlüsselung) automatisch gemergt.
- `test_cases/CMakeLists.txt` — trivial (Leerzeile).

### 3. Repair-Commit (`9701792`)
Der Merge legte zwei latente Defekte aus den Apr-25 Cherry-Picks frei:

**`eeq_solver.h`** — der `ada6925`-Cherry-Pick hatte **nur die `.cpp`-Änderungen, nicht die `.h`-Änderungen** aus `58892a4`+`87c1a8e`+`e37db8b` übernommen. Folge: alle Member-Variablen, Helfer und PARAM-Definitionen fehlten im Header (Build-Fehler `m_allow_unconverged undefined`, etc.). Repariert: vollständigen Header von `e37db8b` übernommen, dann die aggressiven Defaults daraus zurück auf konservative Werte gestellt:

| PARAM | nach `e37db8b` | konservativ wiederhergestellt |
|---|---|---|
| `solve_method` default | `pcg` | `schur_cholesky` |
| `max_pcg_iterations` | 100 | 200 |
| `pcg_large_system_iterations` | 100 | 5000 |
| `m_selected_method` initial | `PCG` | `SchurCholesky` |

**`eeq_solver.cpp`** in `dispatchSolve` — der `7389e2e`-Cherry-Pick (`e37db8b`) wurde mit Konflikt-Markern committed. Nach Auflösung blieben zwei strukturelle Spuren:

- duplizierter `pcg_viable`-Tuning-Block (`if (pcg_viable) { ... pcg_tol_factor ... }` doppelt aufgeführt mit zweiter `bool pcg_viable = true;` Re-Deklaration und doppeltem Gershgorin-Check) — entfernt.
- `schur_ok` in `if (schur_ok)` benutzt ohne Deklaration — `bool schur_ok = true;` plus degenerate-Schur-Komplement-Guard ergänzt (matches stash content).

**`gfnff_method.cpp`** im HB/XB-Detection-Block — der HB-Cache-Code (`m_last_hbonds = params.hbonds`, etc.) lag durch ein Auto-Merge-Artefakt **außerhalb** des `if (m_parameters.value("hbond", true))`-Blocks (mit verwaister `}` und doppeltem `t_hbxb`-Timing). Block korrekt eingerückt, doppelte Zeile entfernt.

## Verifikation

| Item | Status |
|---|---|
| `make -C release -j4` | ✅ exit 0, 100 % built |
| `cli_curcumaopt_*` | ✅ 6/6 |
| `cli_rmsd_*` | ✅ 6/6 |
| `cli_confscan_*` | ✅ 7/7 |
| `cli_gfnff_02_caffeine_gfnff_energy_components` | ✅ |
| `cli_gfnff_gpu_01_gfnff_gpu_singlepoint` | ✅ |
| **`cli_gfnff_01_complex_gfnff_singlepoint`** | ❌ **3.5 mEh Drift** (Toleranz 5 µEh, 700× Überschreitung) |

## Befund: 3.5 mEh Regression in `cli_gfnff_01_complex`

Komponenten-Aufschlüsselung gegen Fortran-XTB-Referenz für `complex.xyz` (231 Atome):

| Komponente | Calc (Eh) | Reference (Eh) | Diff |
|---|---|---|---|
| Bond | -37.025273 | -37.025276 | +3 µEh ✓ |
| Angle | +0.056189 | +0.056189 | ✓ |
| **Torsion** | +0.060880 | +0.067179 | **-6.3 mEh** |
| Coulomb | -0.935935 | -0.935935 | ✓ |
| Dispersion | -0.254176 | -0.254177 | ✓ |
| **HBond** | -0.026207 | -0.022713 | **-3.5 mEh** |
| Inversion | +0.006299 | (nicht aufgelistet) | — |
| BATM | -0.019361 | (nicht aufgelistet) | — |
| **Total** | **-37.244142** | **-37.240652** | **-3.5 mEh** |

**Hauptverdächtige:**

1. **HBond** (-3.5 mEh): passt **exakt zur Total-Diff** in Vorzeichen + Größe. HEAD-eigener `a0ab407`-HB-Re-Detection-Fix ändert HB-Pair-Listen während `calculateEnergy()`. Konsistent mit `GPU_Gradient_Debug_Analysis_2026-04-10.md`-Befund: HB-Atome (Indices 11/15) waren bereits Quelle der GPU/CPU-Diskrepanzen.
2. **Torsion** (-6.3 mEh): teilweise kompensierender Drift. Ursache offen — origin-P-Reihe (P1a D4 cache, P2b CN neighbor list, struct-shrinks) oder Stepwise-API-Refactor.

Reference fasst Repulsion/Inversion/BATM/ATM offenbar anders zusammen als die aktuelle JSON-Decomposition (Reference: explizites `REPULSION=0.893442` Eh; Curcuma JSON: kein separater `Repulsion`-Schlüssel, aber Inversion+BATM nicht in der Reference-Liste). Die Rest-Differenz (~0 nach Berücksichtigung dieser Re-Gruppierung) ist innerhalb numerischer Toleranz.

## Was noch offen ist

- **HB-Pair-Anzahl-Diagnose**: vergleiche `m_last_hbonds.size()` mit/ohne `a0ab407`-Re-Detection auf `complex.xyz`, um zu prüfen ob die Re-Detection in dieser Geometrie tatsächlich Pair-Listen ändert.
- **`feature/gfnff-gpu-optimizations`-CPU-Anteile** (`SpatialCellList` aus `4b7eb24`): rückportierbarer CPU-Speedup, noch nicht extrahiert.
- **GPU-Commit-Audit** (`4b7eb24`, `a06841c`, `1aa018d`): separater Vergleich auf GPU-Pfad-Diskrepanz an HB-Atomen 11/15, vor einem Merge dieser Commits.

## Status der GPU-Branch-Triage

Bezogen auf die ursprüngliche Frage "welche Commits aus `feature/gfnff-gpu-optimizations` sind sinnvoll":

| Commit | Inhalt | Status |
|---|---|---|
| `58892a4` EEQ robustness + AccuracyProfile | CPU-only | **drin** (via Merge-Repair) |
| `87c1a8e` Phase-2 EEQ skip wenn Geometrie unverändert | CPU-only | **drin** (via Merge-Repair) |
| `3610463` Timing-Output | CPU-only | drin (`81a251e` + Merge) |
| `e37db8b` HB/XB skip + Phase-2 cache + D4 cutoff | CPU-only | **drin** (via Merge-Repair) |
| `4b7eb24` SpatialCellList + GPU dlogdcn offload | gemischt | offen — CPU-Anteil rückportierbar |
| `a06841c` GPU HB-CN offload + CN pair list | reine CUDA | offen — Audit ausstehend |
| `1aa018d` Konfigurierbare GPU block size | reine CUDA | offen — Audit ausstehend |
