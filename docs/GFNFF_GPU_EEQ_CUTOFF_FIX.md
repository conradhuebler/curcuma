# GFN-FF EEQ Coulomb-Cutoff — Polymer-Energiediskrepanz und Hellmann-Feynman-Konsistenz

*Claude Generated, Mai 2026 — fix für Polymer (N=1410) GPU/CPU Coulomb-Mismatch.*

## Symptom

Auf einem Polymer-System (N=1410, `test_cases/molecules/larger/polymer.xyz`):

- GPU `gfnff -gpu cuda`: **-203.566 Eh**
- CPU `gfnff`: **-202.588 Eh**
- Δ ≈ **0.978 Eh**, exakt im Coulomb-Term lokalisiert (GPU `-8.286` vs CPU `-7.308`)

Kleine Moleküle (H₂O, Caffeine N=24, Triose N=66, Complex N≈30) zeigten CPU=GPU bit-genau. Der Bug skaliert mit der Atomzahl.

## Ursache

In `eeq_solver.cpp::calculateFinalCharges` (Zeilen 2641-2643) wurde für Systeme mit `natoms > 200` ein 30-Bohr-Cutoff auf die EEQ-Coulomb-Matrix angewandt:

```cpp
double eeq_dist_cutoff = m_config.get<double>("eeq_distance_cutoff", 30.0);
const double EEQ_CUTOFF_SQ = (natoms > 200 && eeq_dist_cutoff > 0.0)
    ? eeq_dist_cutoff * eeq_dist_cutoff : 0.0;
```

Der Cutoff war eine PCG-Conditioning-Optimierung. Mit dem inzwischen aktiven Cholesky-/LDLT-Solver ist er numerisch nicht nötig.

Die GPU-Pfade (`solve()`, `solveWithDeviceRHS()`, `solveWithDeviceRHSAndGPUSchur()`) hingegen riefen `k_eeq_build_matrix` immer mit `cutoff_sq = 0.0` auf. Daraus resultierte:

| N    | CPU EEQ A-Matrix     | GPU EEQ A-Matrix | Konsequenz |
|------|----------------------|------------------|------------|
| ≤200 | full (cutoff inaktiv)| full             | CPU = GPU  |
| >200 | truncated (30 Bohr)  | full             | Δ ≈ 1 Eh   |

## Fortran-Referenz

`external/gfnff/src/gfnff_engrad.F90`, Subroutine `goed_gfnff` (Zeilen 1274-1391):

- Zeilen 1319-1331: A-Matrix wird **ohne jeden Cutoff** gebaut (alle N(N-1)/2 Paare).
- Zeilen 1379-1389: Energie summiert über alle Paare.
- Solver: `dsytrf`/`dsytrs` (LDLT mit Pivoting).

Fortran kennt kein `eeq_distance_cutoff`. Die Referenz-Energie ist die ohne Cutoff: für das Polymer **-203.566 Eh**.

## Hellmann-Feynman-Argument

Wird ein Cutoff nur in einem Teil der Pipeline angewandt, sind die EEQ-Charges nicht stationär bezüglich der berechneten Coulomb-Energie:

- EEQ-Lösung minimiert: `E_eeq(q) = ½ q^T A_truncated q + b^T q` unter `Σq_i = q_total`
- Energie wird berechnet mit: `E_coul,full(q) = ½ q^T A_full q + b^T q`
- Daher: `∂E_coul,full/∂q = A_full q + b = (A_full − A_truncated) q − λ·𝟙 ≠ 0`

Der korrekte Gradient lautet:

```
dE/dx = (∂E/∂x)|_q  +  (∂E/∂q) · (∂q/∂x)
                       └────── ≠ 0 wenn cutoff inkonsistent ──────┘
```

Wird der zweite Term ignoriert (so wie es jeder Code tut, der HF voraussetzt), driftet die MD-Energie. Mit `eeq_distance_cutoff = 0` ist `A_full = A_truncated`, der zweite Term verschwindet, HF ist erfüllt.

## Fix (commits ab Mai 2026)

1. **GPU-Solver-Signaturen** akzeptieren `cutoff_sq` (`eeq_solver_gpu.h/.cu`):
   - `solveWithDeviceRHS(..., double cutoff_sq, bool force_refactor)`
   - `solveWithDeviceRHSAndGPUSchur(..., double cutoff_sq, bool force_refactor)`
   - `solve()` hatte den Parameter bereits.

2. **GPU-Caller** liest `eeq_distance_cutoff` aus der Config und reicht ihn weiter (`gfnff_gpu_method.cpp:493`):
   ```cpp
   const double eeq_cutoff_sq = (N > 200 && m_eeq_distance_cutoff > 0.0)
                                  ? m_eeq_distance_cutoff * m_eeq_distance_cutoff
                                  : 0.0;
   ```
   Logik spiegelt das CPU-Pendant.

3. **EEQSolver-Default** `eeq_distance_cutoff: 30.0 → 0.0` (`eeq_solver.h:875`):
   ```cpp
   PARAM(eeq_distance_cutoff, Double, 0.0,
         "Distance cutoff in Bohr ... (0 = no cutoff, matches Fortran goed_gfnff). "
         "Non-zero values violate Hellmann-Feynman vs. the full Coulomb energy and "
         "degrade MD energy conservation.", "Advanced", {})
   ```

4. **GPU-Caller-Default** ebenfalls auf `0.0` (`gfnff_gpu_method.cpp:48`).

## Verifikation

| Molekül   | Atoms | CPU (Eh)         | GPU (Eh)         | Δ        |
|-----------|-------|------------------|------------------|----------|
| H₂O       | 3     | -0.32728217      | -0.32728217      | 0        |
| Caffeine  | 24    | -4.67273711      | -4.67273706      | 5e-8     |
| Triose    | 66    | -9.91485397      | -9.91485455      | 6e-7     |
| Polymer   | 1410  | -203.56552633    | -203.56553523    | 9e-6     |

Polymer GPU↔CPU jetzt konsistent (vorher 0.978 Eh Diff). Reststreuung im µEh-Bereich ist Floating-Point-Rauschen (atomicAdd-Reihenfolge in CN-Pair-Kernel, akzeptabel).

MD-Energieerhaltung auf Polymer (operator-bestätigt Mai 2026): "sieht gut aus" — der vorher beobachtete übermäßige Energie-Austausch mit dem Thermostat ist verschwunden.

## Performance-Optimierung mit Cutoff: Anforderungen für künftige Experimente

Der Cutoff bleibt als CLI-Parameter verfügbar (`-eeq_distance_cutoff <Bohr>`) für Experimente bei sehr großen Systemen, wo der `O(N²)`-Aufbau der EEQ-Matrix oder der Coulomb-Pair-List spürbar wird.

**Wenn ein Cutoff > 0 verwendet wird, MUSS er an allen vier Sites konsistent angewandt werden, sonst bricht der Gradient:**

1. **EEQ A-Matrix** (gegeben):
   - CPU: `eeq_solver.cpp:2641-2706` (in `calculateFinalCharges`)
   - GPU: `k_eeq_build_matrix` (gfnff_kernels.cu / eeq_solver_gpu.cu) via `cutoff_sq` Argument

2. **Coulomb-Pair-List** (Filter beim Erzeugen):
   - `gfnff_method.cpp:8220` (`c.r_cut = 100.0`) — Cutoff hier auf min(100, eeq_distance_cutoff) absenken, sonst kennen die Coulomb-Kernel die in der EEQ ignorierten Paare nicht und summieren sie trotzdem auf.

3. **CPU-Coulomb-Energie + Gradient**:
   - `forcefieldthread.cpp:CalculateGFNFFCoulombContribution` — filtert pro Paar via `coul.r_cut`, also automatisch konsistent wenn Punkt 2 gemacht ist.

4. **GPU-Coulomb-Kernel**:
   - `k_coulomb` (gfnff_kernels.cu:334-379) — filtert pro Paar via `r_cut[tid]`, ebenfalls automatisch konsistent über Punkt 2.

**Validierung jedes Cutoff-Experiments:**

- Energie-Match CPU↔GPU innerhalb 1e-5 Eh (für jede getestete N).
- Gradient-Match CPU↔GPU innerhalb 1e-6 Eh/Bohr (`test_gfnff_gradients`).
- MD-Energie-Drift in NVE über ≥10 ps unter 1e-4 Eh/ps; ohne sichtbares Aufschaukeln im NVT-Thermostat-Austausch.
- Vergleich gegen XTB-Referenz auf einer Auswahl von Systemen — die Referenz hat **keinen** Cutoff, jeder Cutoff produziert eine andere Energiekurve. Daher Cutoff-Energien nur untereinander vergleichen, nicht gegen XTB.

**Self-Energie und Self-Interaction-Terme** (`goed_gfnff` Zeilen 1387-1388) sind atomweise Beiträge ohne Paardistanz — vom Cutoff unberührt.

## Verwandte Dokumente

- `docs/GFNFF_STATUS.md` — GFN-FF Gesamt-Status
- `docs/GFNFF-PERFORMANCE-ROADMAP.md` — Performance-Optimierungsplan
- `docs/GPU_WP5_FULL_GPU_EEQ_PIPELINE.md` — WP5-A GPU-Schur-Komplement-Pfad

## Geänderte Dateien

- `src/core/energy_calculators/ff_methods/eeq_solver.h` (Default 30.0 → 0.0)
- `src/core/energy_calculators/ff_methods/cuda/eeq_solver_gpu.h` (Signatur)
- `src/core/energy_calculators/ff_methods/cuda/eeq_solver_gpu.cu` (cutoff_sq weiter­reichen)
- `src/core/energy_calculators/qm_methods/gfnff_gpu_method.h` (Member, Default 0.0)
- `src/core/energy_calculators/qm_methods/gfnff_gpu_method.cpp` (Cutoff-Berechnung, drei Aufrufstellen)
