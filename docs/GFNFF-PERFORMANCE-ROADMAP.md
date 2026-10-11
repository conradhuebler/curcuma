# GFN-FF Performance Roadmap

Konsolidiert aus: GFNFF-CPU.md (CPU-Bottleneck-Analyse) und GFNFF-GPU.md (GPU-Analyse, März-April 2026).

---

## Ziel

Messbare Speedups für MD-Simulationen (>200 Atome) und Batch-Optimierungen ohne physikalische Korrektheit zu riskieren.
Referenz-Validierung: numerischer vs. analytischer Gradient muss nach jeder Änderung bestehen.

---

## Priorisierte Maßnahmen

### Stufe 1 — Kurzfristig (jeweils 2-4h, hohes Risiko/Nutzen-Verhältnis)

| ID  | Beschreibung | Datei | Aufwand | Nutzen |
|-----|-------------|-------|---------|--------|
| P1a | **D4 Gaussian-Gewichte Threshold-Cache** — `if (max_cn_change < 0.01) skip` | `gfnff_method.cpp` | 2-3h | ~50% D4-Term (MD) |
| P1b | **Round-Robin Thread-Split** — ersetze linearen Strip in AutoRanges | `forcefield.cpp:1545-1640` | 1-2h | 10-20% Gesamt (4+ Threads) |
| P1c | **GFNFFDispersion Struct schrumpfen** — 5 Legacy-Felder (`C8`, `s6`, `s8`, `a1`, `a2`) auslagern; 88 → 48 Bytes | `gfnff_parameters.h:83-97` | 3-4h | 20-30% Dispersion-Term |

### Stufe 2 — Mittelfristig (je ~1 Tag, hoher Nutzen bei größeren Systemen)

| ID  | Beschreibung | Datei | Aufwand | Nutzen |
|-----|-------------|-------|---------|--------|
| P2a | **N×N Distanzmatrix eliminieren** — CN-Cutoff-Loop statt vollständiger Matrix; `TopologyInfo.distance_matrix` → Optional | `gfnff_method.cpp:651-693` | 1 Tag | −90% Speicher bei N>500 |
| P2b | **CN Neighbor-List mit Cutoff 6 Bohr** — ersetzt O(N²)-`erf()`-Loop | `gfnff_method.cpp:669-693` | 1 Tag | 5-10× CN-Term N>200 |

### Stufe 3 — Langfristig (je 4-8h, moderater Nutzen)

| ID  | Beschreibung | Datei | Aufwand | Nutzen |
|-----|-------------|-------|---------|--------|
| P3a | **GFNFFCoulomb: `q_i`/`q_j` aus Struct** — direkt `charges[i]` nutzen; 144 → 128 Bytes | `gfnff_parameters.h:118-131` | 4-6h | ~10% Coulomb-Term |
| P3b | **`std::set` → Bitset Bonded-Ausschluss** — obere Dreiecksbitmatrix, O(1)-Lookup | `forcefieldthread.h:609` | 2-3h | 5-10% NB-Terme |

---

## GPU-spezifische Maßnahmen

### Stufe G1 — Kurzfristig (je 2-4h)

| ID  | Beschreibung | Datei | Aufwand | Nutzen |
|-----|-------------|-------|---------|--------|
| G1a | **Dispersion-Paare nach `idx_i` sortieren** — vor Upload, O(N log N) einmalig; verbessert L2-Locality | `ff_workspace_gpu.cu` (Upload-Phase) | 2-3h | 5-15% Dispersion |
| G1b | **`__ldg()` für SoA-Pair-Parameter** — `__ldg(&C6[tid])`, `__ldg(&r4r2ij[tid])` etc. | `gfnff_kernels.cu` (alle Pair-Kernel) | 2h | 3-8% L1-Entlastung |

### Stufe G2 — Mittelfristig

| ID  | Beschreibung | Datei | Aufwand | Nutzen |
|-----|-------------|-------|---------|--------|
| G2a | **`k_dihedrals` / `k_angles` Register-Druck reduzieren** — große Kernel aufteilen oder Hilfsvariablen in shared mem auslagern | `gfnff_kernels.cu:564-1037` | 1-2 Tage | +30% Occupancy |
| G2b | **Nsight Compute Profiling** — Register-Spilling und Speicherbandbreite empirisch messen, bevor G2a begonnen wird | — | 2-4h | Validierung der Hypothesen |
| G2c | **Konsistenter EEQ + Coulomb Distance-Cutoff** (CPU **und** GPU, Energie **und** Gradient) — opt-in via `eeq_distance_cutoff > 0`, derzeit Default 0 (matcht Fortran). **Achtung**: Cutoff weicht von Fortran-Referenz ab und verändert das Energie-Surface, daher nur als Performance-Experiment für sehr große Systeme; siehe `docs/GFNFF_GPU_EEQ_CUTOFF_FIX.md` für die vier Pflicht-Sites (EEQ A-Matrix, Coulomb-Pair-List, CPU-Coulomb, GPU-Coulomb). Wenn auch nur **eine** Site den Cutoff anders behandelt → Hellmann-Feynman-Verletzung → MD-Drift. | `eeq_solver.cpp`, `gfnff_method.cpp:8220`, `forcefieldthread.cpp` (Coulomb), `eeq_solver_gpu.cu` + `gfnff_kernels.cu` (`k_coulomb`) | 1 Tag | Moderate Coulomb-Speedup bei N>5000; nicht der Bottleneck (das bleibt EEQ-Cholesky O(N³/6)). Niedrige Priorität bis Cholesky selbst optimiert ist. |

---

## Nicht empfohlen

- OpenMP/SIMD-Pragmas CPU: CxxThreadPool + Eigen-SIMD reichen. Kein messbarer Mehrwert.
- Coulomb/Dispersion-Cutoff einseitig kürzen: Inkonsistenz zwischen Energie und EEQ-Lösung verletzt Hellmann-Feynman → MD driftet (verifiziert Mai 2026, Polymer N=1410). Nur über die **konsistente** Maßnahme G2c oben anpacken — und dann zusätzlich gegen XTB validieren, weil der gesamte Cutoff-Pfad nicht-Fortran ist.
- CUTLASS / Tensor Cores: GFN-FF ist kein Matrix-Problem.
- Texture Memory: Legacy-API; L2-Cache reicht für aktuelle Kernel-Balance.
- Kooperative CUDA-Gruppen: Kein Gewinn gegenüber bestehendem Warp-Shuffle-Reduktionscode.

---

## Verifikation nach jeder Änderung

```bash
cd release && make -j4
ctest -R "gfnff" --output-on-failure
./test_cases/test_gfnff_numgrad    # Gradient-Regression (CPU)
# GPU:
./curcuma -sp test_cases/molecules/larger/acetic_acid_dimer.xyz -method gfnff -gpu cuda
```

Performance-Messung (Polymer, 1410 Atome, 1000 Schritte):
```bash
cd test_cases/cli/simplemd/10_gfnff_polymer_md
bash compare_baseline.sh 4 gfnff     "after_P1a_d4_threshold"
bash compare_baseline.sh 4 gfnff-gpu "after_G1a_pair_sorting"
```

---

## Abhängigkeiten zwischen Maßnahmen

```
P2a (Distanzmatrix) → P2b (CN Neighbor-List): P2b vereinfacht sich nach P2a
P1c (Disp-Struct)   → G1a (Pair-Sortierung): G1a ist einfacher wenn Struct kleiner
G2b (Profiling)     → G2a (Kernel-Split): Profiling zuerst, um richtige Kernel zu priorisieren
```

---

## Status

| ID | Status | Ergebnis |
|----|--------|---------|
| P1a | ⚙️ Machine-tested | CPU: 1.01× (156.3→157.6 ms/step je nach Threshold), GPU: 1.03× — Cache bei 1410 Atomen selten aktiv. Threshold-Benchmark: 0.0=156.3ms, 0.01=157.6ms, 0.05=156.8ms, 0.1=155.2ms — alles innerhalb Messungenauigkeit. D4-Term ist kleiner Anteil der Gesamtzeit. |
| P1b | ⚙️ Machine-tested | CPU: 1.04× (157.7→151.1 ms/step) — Round-Robin für Dispersion/Coulomb/Repulsion/D4. Effekt minimal bei 1410 Atomen: Cache-Lines bereits ausreichend groß. Erwarteter Nutzen tritt erst bei >4 Threads oder unbalancierter Last auf. |
| P1c | ⚙️ Machine-tested | CPU: 1.00× (157.9 vs 157.7 ms/step) — Struct 88→48 Bytes, kein messbarer Effekt bei 1410 Atomen. Cache-Limitierung tritt erst bei größeren Systemen auf. |
| P2a | ⚙️ Machine-tested | CPU: 1.05× (157.7→150.5 ms/step) — distance_matrix aus Per-Step-Pfad entfernt, squared_dist_matrix eliminiert (dead weight), CNCalculator statt D3-Style CN in updateDynamicState(), countNeighborsWithin20Bohr on-the-fly, HuckelSolver erhält geometry_bohr statt N×N-Matrix. Speicherersparnis: ~32MB bei 1410 Atomen. Performance-Effekt minimal bei 1410 Atomen (Cache-Limitierung nicht erreicht). |
| P2b | ⚙️ Machine-tested | CPU: 1.08× (157.7→146.6 ms/step, 6 Bohr cutoff) — Neighbor-List mit konfigurierbarem Cutoff. Referenzmodus (voll O(N²)): 159.6 ms/step. Parameter: `-gfnff.cn_cutoff_bohr 6` (default), `-gfnff.cn_cutoff_bohr 0 -gfnff.cn_accuracy 0` (Referenz). CPU-GPU-Inkonsistenz (30 vs 40 Bohr) dokumentiert. |
| P3a | ⚙️ Machine-tested | CPU: 0.99× (159.8 vs 157.7 ms/step) — GFNFFCoulomb 120→56 Bytes, redundante Felder entfernt (q_i/j, chi_i/j, gam_i/j, alp_i/j). Per-Atom-Vektoren statt Per-Pair-Daten. Kein messbarer Effekt bei 1410 Atomen. |
| P3b | ⚙️ Machine-tested | CPU: 0.99× (159.8 vs 157.7 ms/step) — m_bonded_pairs Dead Code entfernt, std::set→vector<bool> Bitset in generateRepulsionPairsNative(). Cache-Effekt minimal bei 1410 Atomen. |
| G1a | ⚙️ Machine-tested | CPU: 0.99× (159.8 vs 157.7 ms/step) — Dispersion-Paare nach idx_i sortiert vor GPU-Upload. Kein CPU-Effekt. GPU: 1.50× (20.0 vs 30.0 ms/step) — deutlich bessere L2-Lokalität durch Sortierung. |
| G1b | ⚙️ Machine-tested | CPU: N/A — __ldg() nur für GPU relevant. GPU: kombiniert mit G1a-Messung (1.50×). __ldg() für alle Paar-SoA-Parameter in k_dispersion, k_repulsion, k_coulomb, k_coulomb_self, k_coulomb_postprocess. |
| G1c | ⚙️ Machine-tested | CPU: N/A. GPU: 1.58× (19.0 vs 30.0 ms/step) — Coulomb + Repulsion-Paare ebenfalls nach idx_i sortiert. Kein zusätzlicher Effekt bei 1410 Atomen (Sortierung bringt nur bei Dispersion Gewinn wegen engerem Cutoff). |
| G2a | ⚙️ Machine-tested | CPU: 1.06× (148.8 vs 157.7 ms/step) — kein CPU-Effekt. GPU: 1.58× (19.0 vs 30.0 ms/step) — Register-Spilling reduziert durch `__launch_bounds__(256, 2)` für k_angles, k_dihedrals, k_inversions, k_storsions, k_hbonds. k_angles/storsions/inversions: 0 Spills (von 88-264 Bytes Stack). k_dihedrals: 88 Bytes Stack (von 432). k_hbonds: 312 Bytes Stack (von 632), noch 548/780 Bytes Spill — benötigt weitere Optimierung. |
| G2b | ⏭️ Übersprungen | Nsight Compute nicht installiert. Register-Spilling stattdessen via ptxas -v analysiert (siehe G2a-Ergebnisse). |
| G2c | 🟡 Nicht aktiv | Default `eeq_distance_cutoff = 0` (matcht Fortran) seit Mai 2026. Ein Cutoff-Experiment muss alle vier Sites konsistent abdecken — sonst HF-Verletzung und MD-Drift (verifiziert auf Polymer N=1410: 30-Bohr-Cutoff einseitig in EEQ-Matrix führte zu 0.978 Eh Coulomb-Mismatch GPU/CPU + spürbarem Thermostat-Energieaustausch). Siehe `docs/GFNFF_GPU_EEQ_CUTOFF_FIX.md`. |

*Nur der Operator darf Status auf ✅ TESTED setzen.*

---

## WP-P1 MD Per-Step-Profil (Mai 2026)

Per-Step-Wall-Clock-Breakdown via `md_diagnostics_timing=true` Hook
(WP-P1). Quelle: [GFNFF_PROFILE_RESULTS_2026-05.md](GFNFF_PROFILE_RESULTS_2026-05.md).

**Polymer NVE (N=1410, dt=1 fs, 4 Threads), Step 5**:

| Phase | Zeit | Anteil step_total | Hebel |
|-------|------|--------------------|-------|
| `dcn` (CN-Derivative pair list) | 23.3 ms | 26 % | **Neu erkannter #1 Bottleneck** — bisher unterschätzt |
| `eeq_solve` (Phase 2) | 22.2 ms | 24 % | WP-S1 `static_charges` + WP-S3 `cutoff_auto` |
| `cn` (Coordination Numbers) | 11.4 ms | 13 % | P4e SIMD, OpenMP-Audit |
| FF-Terme zusammen | 20.3 ms | 22 % | Bond / Angle / Coulomb / Disp; nur ~7 ms Bond-Anteil |
| `d4_gw` (D4 Gaussian weights) | 3.4 ms | 4 % | P4f SIMD |

**Korrektur zur Mai-2026-Mixture-Hypothese**:
- P4c "75 µs/Bond" (mixture-Profil) reproduziert sich **nicht** auf Polymer.
  Gesamte FF-Term-Berechnung ist 20.3 ms, nicht ~225 ms wie aus 75 µs ×
  ~3000 Bonds extrapoliert. Mixture-Befund war wahrscheinlich Memory-
  Pressure-getrieben bei N=6200 — WP-P3-Mikrobench würde das klären.

**Konsequenzen**:
- **WP-S1 `static_all`** ist nun durch Profil-Daten gerechtfertigt:
  `dcn + eeq_solve + cn + d4_gw = 60.3 ms = 66 %` der step_total werden
  durch Static-Mode wegrationalisiert. Passt zum gemessenen 3.2×-Speedup
  im Quick-Test.
- **Neuer high-priority WP-Kandidat**: `dcn`-Vektorisierung (SIMD oder
  Pair-Cache-Optimierung) für die `calculateCoordinationNumberDerivatives()`-
  Inner-Loop.
- **WP-P3 Bond-Microbench** wird sekundär; WP-P2 Cross-System-Sweep wird
  priorisiert, um Polymer-vs-Mixture-Diskrepanz zu klären.

---

## WP-P2 Cross-System Sweep (Mai 2026)

Treiber: `scripts/wp_profile_sweep.py`, MD csvr 200 fs/dt=1 fs, 4 Threads.
Vollständiger Datenpunkt: caffeine+polymer × {baseline, static_all,
cutoff_auto, both}. Mixture-MD ist bei dt=1 fs instabil und wurde
skippt; GPU-Pfad nicht gemessen (USE_CUDA=OFF in der Standard-release/).

| System | Modus | step_total | dcn | eeq_solve | speedup |
|--------|-------|------------|-----|-----------|---------|
| caffeine | baseline    | 0.38 ms |  0.01 |  0.08 | 1.00× |
| caffeine | static_all  | 0.15 ms |  0.00 |  0.00 | **2.55×** |
| caffeine | cutoff_auto | 0.29 ms |  0.01 |  0.06 | 1.29× |
| polymer  | baseline    | 102.1 ms | 21.0 | 21.6 | 1.00× |
| polymer  | static_all  |  28.6 ms |  0.0 |  0.0 | **3.57×** |
| polymer  | cutoff_auto |  99.0 ms | 21.0 | 20.2 | 1.03× |

**Hauptbefunde**:
- **WP-P1-Polymer-Snapshot reproduziert über 20 Records** (dcn ≈ 21, eeq_solve ≈ 21.6).
- **`static_all` ist der Workhorse**: polymer 3.6×, caffeine 2.5× — skala-abhängig.
- **`cutoff_auto` allein ist nutzlos bei mittelgroßen Systemen** (polymer 1.03×).
  Lohnt sich erst bei N>4000 (mixture-Klasse) wenn Phase-2-EEQ dominiert.
- **Bond-Hotspot-Hypothese definitiv falsifiziert** für polymer-Klasse:
  alle FF-Energie-Terme zusammen = 43 ms ≪ 75 µs × 3000 Bonds = 225 ms-Hypothese.
  Mai-2026-mixture-Befund war Memory-Pressure-getrieben.

**Konsequenzen** (Details in [GFNFF_PROFILE_RESULTS_2026-05.md](GFNFF_PROFILE_RESULTS_2026-05.md)):
- Empfehlung: **`dcn`-Vektorisierung** (~21 ms Phase, ~2× Hebel) als
  high-priority Folge-WP.
- `cn`-SIMD (~11 ms Phase, ~2× Hebel) als Sekundär-Folge.
- WP-S3 (`eeq_distance_cutoff_auto`) sollte als "großes-System-Feature"
  dokumentiert werden, nicht als Default-Win.
- WP-P3 Bond-Mikrobench bleibt niedrig priorisiert.

---

## WP-D Stage A — cn_raw-Reuse (Mai 2026)

Erste Umsetzung der dcn-Optimierung. `CNCalculator::calculateGFNFFCN` gibt
optional `cn_raw` zurück; `calculateCoordinationNumberDerivatives` skippt
seine eigene N²-erf-Schleife, wenn der gecachete Wert übergeben wird.

| Phase | Pre-WP-D | Stage A | Delta |
|-------|----------|---------|-------|
| polymer `dcn` (baseline) | 21.0 ms | 18.0 ms | −14 % |
| polymer `step_total` | 102.1 ms | 100.9 ms | −1.2 % |

**Kleiner als erwartet**: Step 1 war seit März 2026 OpenMP-parallel; bei
4 Threads ist die Wiederholung schon günstig. **Bit-identische Energien**
(polymer SP −203.51758419 Eh pre/post). gfnff_numgrad-Tests bleiben grün.

**Stage B (Pair-Cache) deferred**: würde nur ~1–2 ms zusätzlich sparen
für ~120 LoC. Effort/Reward schlecht. **Echter Folge-WP-Kandidat**: SIMD-
Vektorisierung der `std::exp()`-Loop in dcn step 3 (~10 ms Hebel via
SLEEF / Polynomial-Approximation).

Details: [GFNFF_PROFILE_RESULTS_2026-05.md](GFNFF_PROFILE_RESULTS_2026-05.md) WP-D-Sektion.

---

## WP-D Stage B — SIMD `exp(-x²)` in dcn step 3 (Mai 2026)

Hand-rolled AVX2-Remez-Polynom-Kernel (`fast_exp.h/.cpp`) plus Pass-A/B/C
Loop-Restruktur in `calculateCoordinationNumberDerivatives()`. Opt-in via
`USE_GFNFF_FAST_EXP=OFF`-Default — Default-Build bleibt bit-identisch zu
pre-WP. Polymer A/B auf gleichem Hash (`beacabe`) und System-Load:

| Phase | OFF Reference | ON SoA-v2 | Delta |
|-------|---------------|-----------|-------|
| polymer `dcn` median | 19.39 ms | 17.21 ms | −2.18 ms (−11 %) |
| polymer 1000-Step Wall | 105.7 ms/step | 104.2 ms/step | −1.5 ms (−1.4 %) |

**Kleiner als die 10-ms-Hebel-Erwartung des Stage-A-Eintrags**, weil
`exp()` mit ~0.2 ms Gesamt-Anteil nicht der `dcn`-Bottleneck war —
dominant sind Memory-Reads/Writes pro überlebendem (i,j)-Paar.

Erkenntnis: v1 (push_back in 7 std::vectors) machte `dcn` **schlechter**
(21.2 ms). v2 (Eigen-Matrix SoA + counter `K`, fixed-size pre-alloc)
fängt den Pass-A-Overhead ab. Layout-Effekt > exp-Vektorisierung selbst.

**`fast_exp_neg_sq_block`** ist als wiederverwendbarer Helper angelegt,
nicht Teil dieses WP, aber gelistet als Roadmap-Anker:
`dcn` step 1, `dcn` step 2, CN erf-Loop, D4-Gaussian-Weights, Angle-Damping.

Details: [docs/wp4/WP-D-Stage-B-dcn-simd-exp.md](wp4/WP-D-Stage-B-dcn-simd-exp.md).

| Datum | Variante | dcn mean | Wall 1000-step |
|-------|----------|----------|----------------|
| Roadmap | Pre-WP-D | 21.0 ms | 102.1 ms |
| Roadmap | Stage A (cn_raw reuse) | 18.0 ms | 100.9 ms |
| beacabe 2026-05 | Stage B v1 push_back | 21.2 ms | 105.0 ms |
| beacabe 2026-05 | Stage B v2 SoA (ON) | 17.9 ms | 104.2 ms |
| beacabe 2026-05 | Reference (OFF) | 19.8 ms | 105.7 ms |

---

## Kumulative Bilanz WP1–WP5 (Mai 2026)

`mixture.xyz` (N=6200, nfrag=1400), Single-Point + Gradient, AMD Ryzen 9 9950X3D, Topo-Cache vor jedem Lauf gelöscht:

| T | Pre-WP1 (Total Wall) | Post-WP4 | Post-WP-G | Post-WP5 | **Cumulative Speedup** |
|---|----------------------|----------|-----------|----------|------------------------|
| 1 | 3074 ms | 2269 ms | 2099 ms | **2078 ms** | **1.48×** |
| 4 | 1840 ms | 1012 ms | 967 ms | **933 ms** | **1.97×** |
| 8 | 1673 ms | 730 ms  | 708 ms  | **712 ms**  | **2.35×** |
| 16| 1469 ms | 669 ms  | n.g.    | n.g.    | **~2.4×** |

**Single-Thread (T=1) gewann 996 ms (32 %)** ohne Threading-Maßnahme — kommt aus
- WP4 (CN-deriv Pair-Liste, 1206→426 ms bei T=1)
- WP-G (Layout, ~70 ms)
- WP5 (Heap-Reduktion, ~21 ms)

**Multi-Thread (T=4, Hauptarbeitsmodus): fast halbiert** (1840→933 ms = 49 % Ersparnis).

WP1, WP2, WP3 lieferten keinen direkten Speedup, aber:
- WP1 fixte den `setThreadCount`-Bug (`(void)threads` → ignorierte das CLI-Argument), Voraussetzung für korrekte Threading-Skalierung
- WP2 räumte die EEQ-Pool-Propagation strukturell auf, Voraussetzung für künftige Multi-Fragment-Optimierungen
- WP3 vereinheitlichte den Bond-Inner-Loop mit dem Repulsion/D3/ATM-Pattern (Stack-Vector statt Heap-Matrix)

**Open issues (nicht performance):**
- WP-R: Multi-Thread-Race auf triose T=4 (1.47 mEh Drift, pre-existing, eigenes WP)
- WP-V: Gradient-Validierung CPU vs GPU vs Fortran (MD-heat-bath-Drift, pre-existing, eigenes WP)
- WP6: Coulomb-Cutoff (blockiert auf G2c)

---

## Ausgangsmessung — 2026-04-11 (git 1c27621)

**System:** AMD Ryzen 9 9950X3D, NVIDIA GeForce RTX 5080, 4 Threads
**Testfall:** Polymer, 1410 Atome, 1000 MD-Schritte (GFN-FF, CSVR, Seed 42)

| Methode | Wall-Zeit | ms/Schritt | Faktor vs. GPU |
|---------|-----------|-----------|----------------|
| `gfnff-gpu` (`-method gfnff -gpu cuda`) | 30.0 s | 30.0 ms | 1× |
| `gfnff` (nativ CPU, 4 Threads) | 157.7 s | 157.7 ms | 5.3× langsamer |
| `xtb-gfnff` (Fortran/XTB, 4 Threads) | 422.1 s | 422.1 ms | 14.1× langsamer |

**Bewertung:**
- Native CPU-Implementierung ist bereits **2.7× schneller** als xtb-gfnff.
- GPU ist **5.3× schneller** als native CPU — Ziel der CPU-Optimierungen ist es, diesen Abstand auf ≤3× zu verkleinern.
- Bei 1410 Atomen dominieren O(N²)-Terme (Dispersion, Coulomb, CN): P1c, P2a, P2b haben hier den größten Hebel.
- Rohdaten: `test_cases/cli/simplemd/10_gfnff_polymer_md/baseline.dat`

---

## WP4 — CPU-Pfad Multi-Fragment-Profil (Mai 2026)

Nachfolger zur EEQ-Multi-Fragment-Optimierung (Commit `b6a3a42`, Mai 2026). Profil-Datenpunkt: `mixture.xyz` (N=6200, nfrag=1400), single-point + gradient, CPU-Pfad, 4 Threads.

### Profil

| Phase | Zeit | Anteil |
|-------|------|--------|
| **Total energy call** | **3341 ms** | 100% |
| CN + EEQ (serial CPU) | 1766 ms | 53% |
|   davon EEQ Phase 2 (Stage 4 batched) | ~660 ms | 20% |
|   davon Param/CN/D4-Gaussian | ~1100 ms | 33% |
| Force Field Energy (multi-threaded) | 1574 ms | 47% |
|   Coulomb | 640 ms | 19% |
|   Bond | 359 ms | 11% |
|   H-bonds | 205 ms | 6% |
|   Dispersion | 185 ms | 6% |
|   Repulsion (nonbond) | 60 ms | 2% |

GPU-Vergleich für dasselbe System: 1.5 s (2.2× schneller als CPU). Force-Field-Energie auf GPU: ~25 ms (63× schneller als CPU). Pro OPT-Schritt zusätzlich ~17 s für `RMSDDriver::start()` (Kabsch+Topologie-Rebuild) — bereits gefixt in Commit `b6a3a42` via O(N) Coordinate-Diff für N≥500.

### Geplante Maßnahmen

#### Stufe 4a — Quick Wins (je 1-3h, hohes Risiko/Nutzen-Verhältnis)

| ID  | Beschreibung | Datei | Aufwand | Geschätzter Nutzen |
|-----|-------------|-------|---------|---------------------|
| P4a | **EEQ Batched parallelisieren** — 1400 unabhängige Fragment-LU via `CxxThreadPool` oder OpenMP. Embarrassingly parallel. | `eeq_solver.cpp` (Stage-4-Loop, Mai 2026 Commit) | 1-2h | Phase-2 4-6× (660→120-200 ms) |
| P4b | **Threading-Audit aller "serial CPU"-Phasen** — `CN + EEQ (serial CPU)` ist explizit als seriell gelabelt; Default-Thread-Anzahl im OPT-Pfad prüfen, sicherstellen dass `-threads N` durchgereicht wird. | `gfnff_method.cpp` Phase-2-Loop, capability-Layer | 2h | 4-8× insgesamt wenn Threads <2 |
| P4c | **Bond-Term Hotspot-Profiling** — 359 ms / 4800 Bonds = 75 µs/Bond ist ~100× zu langsam für eine simple `exp((r-r0)·α)` Auswertung. Vermutung: Cache-Misses durch Pointer-Indirection oder einzelne Eigen-Matrix-Ops pro Bond statt Batch. `perf record` zur Lokalisierung. | `forcefieldthread.cpp:CalculateGFNFFBondContribution` | 1-2h Profiling, dann 2-4h Fix | Bond 5-10× wenn Hypothese stimmt |

#### Stufe 4b — Mittelfristig (je ~1 Tag, hoher Nutzen bei N>2000)

| ID  | Beschreibung | Datei | Aufwand | Geschätzter Nutzen |
|-----|-------------|-------|---------|---------------------|
| P4d | **Coulomb Cell-List / Verlet-List** — räumliche Lokalisierung statt voller O(N²)-Iteration. `erf(γ·r)/r` < 1e-7 jenseits ~12 Bohr; bei mixture sind die meisten Paare > 12 Bohr entfernt. Alternative: einfache Cutoff-Liste (ohne Cell-List) reicht für mittelgroße Systeme. **Achtung**: Konsistenz mit EEQ-Matrix nötig (siehe G2c-Warnung), sonst Hellmann-Feynman-Verletzung. | `forcefieldthread.cpp:CalculateGFNFFCoulombContribution` | 1 Tag | Coulomb 5-20× bei N>2000 |
| P4e | **CN-Berechnung SIMD-vektorisieren** — `std::erf()` per-Atom-Loop ersetzen durch Eigen array op auf Distance-Matrix-Slice. Funktioniert direkt mit Eigen-AVX2/AVX-512 Vector-Math. | `gfnff_method.cpp:calculateCN` (oder `cn_calculator.cpp`) | 4-6h | CN ~2× |
| P4f | **D4 Gaussian-Weights vektorisieren** — 144 ms / 6200 Atome / 7 Refs = 3.3 µs pro Atom-Reference-Paar. Inner Loop (MAX_REF=7) als Eigen array op, exp() vektorisiert. | `d4param_generator.cpp:precomputeGaussianWeights` | 3-4h | D4-Weights ~2× |

#### Stufe 4c — Langfristig (Architektur-Eingriff)

| ID  | Beschreibung | Datei | Aufwand | Geschätzter Nutzen |
|-----|-------------|-------|---------|---------------------|
| P4g | **Cell-List statt Pair-Listen für alle Non-Bonded-Terme** — gemeinsame räumliche Datenstruktur für Coulomb, Dispersion, Repulsion, H-bonds. Komplett O(N) statt O(N²). | mehrere Dateien, Architektur-Refactor | 1-2 Wochen | 10-30× bei N>5000 |

### Realistisches CPU-Ziel

Mit Stufe 4a+4b umgesetzt: Single-point CPU von 3.3 s → ~0.6-0.8 s (mixture.xyz). Faktor vs. GPU sinkt von 2.2× auf ~1.3-1.5×. Bei kleineren Systemen (N<1000) wird der Abstand komplett verschwinden.

### Abhängigkeiten

```
P4b (Threading-Audit)  →  P4a (EEQ batched parallel): wenn Threads schon korrekt
                                                       durchgereicht werden, P4a
                                                       ist trivial
P4c (Bond-Profiling)   →  ggf. weitere FF-Term-Optimierungen via gleichem Muster
P4d (Coulomb Cutoff)   →  G2c (konsistenter EEQ-Cutoff) muss zuerst aktiviert
                          und gegen XTB validiert sein
```

### Status

| ID | Status | Ergebnis |
|----|--------|---------|
| P4a | ⚙️ Machine-tested | Stage-4-Loop in `eeq_solver.cpp:1294-1387` mit CxxThreadPool parallelisiert (Pool an alle 6 EEQ-Aufrufer durchgereicht). `mixture.xyz` (N=6200, nfrag=1400) Sweep: Total 3068→1479 ms (T=1→16, 2.07×) — **kein nennenswerter Gewinn vs. WP1-Baseline** (war 3074→1469 ms = 2.09×). Phase 2 nur 1.17× über T=1→16. Hebel war im Mai-Plan überschätzt: die 660-ms-Profilzahl gilt heute nicht mehr (Phase 2 ist nur 371 ms wall, davon LU-Loop ~80-100 ms). Per-Fragment-LU mit 4-5 Atomen ist zu klein für effizientes Threading. Strukturell wertvoll als Voraussetzung. Vollständiger Bericht: `docs/wp4/WP2-eeq-batched-parallel.md`. |
| P4b | ⚙️ Machine-tested | **Bug gefunden + behoben** in `qm_methods/gfnff_method.cpp:97-99` (`setThreadCount` ignorierte den Parameter). Neue `m_threads`-Member in `GFNFF`. Label `"CN + EEQ (serial CPU)"` → `"CN + EEQ"` korrigiert. Skalierungs-Sweep auf `mixture.xyz` T={1,2,4,8,16}: Total 3074→1469 ms (2.09×), FF-Pool skaliert 5.0× von T=1→8 (parallel eff 78.6→41.3 %). **Keine** OpenMP-Pathologie wie im Fortran-Original. Schwachpunkt: CN-Derivate skalieren nur 1.31× — Kandidat für WP4 oder eigenes Sub-WP. Stage-4 Per-Fragment-LU in `eeq_solver.cpp:1305-1335` als einzige echt serielle Stelle bestätigt — direktes WP2-Ziel. Vollständiger Bericht: `docs/wp4/WP1-threading-audit.md`. |
| P4c | ⚙️ Machine-tested | `Matrix derivate;` Heap-Allocation in `forcefieldthread.cpp:934` durch Stack-`Eigen::Vector3d` ersetzt (analog Repulsion/D3/ATM-Pattern). Korrektheit bit-identisch (`acetic_acid_dimer` -2.47129863). **Kein messbarer Speedup** auf `mixture.xyz` (Bond cpu-time T=4: 173 ms post-WP3 vs. 170 ms pre-WP3, im Noise). Heap-Allocation war **nicht** der Hotspot — 48-byte malloc bei moderner glibc ~20 ns × 4800 Bonds = 0.1 ms cpu-time. Echter Hotspot vermutlich `m_gradient.row()` Cache-Misses (ColumnMajor MatrixXd Stride N×8 = 49.6 KB) oder `std::exp()` selbst. Vollständiger Bericht: `docs/wp4/WP3-bond-hotspot.md`. Folge-WP: `m_gradient` auf RowMajor oder Per-Thread-Buffer (würde alle FF-Terme treffen). |
| P4d | 🆕 Vorgeschlagen | Blockiert durch G2c |
| P4e | ⚙️ Machine-tested | **Größter WP-Gewinn bisher.** Original-Plan (CN-SIMD) verworfen, weil WP1 zeigte: CN-Berechnung selbst nur 109 ms, Hotspot waren CN-**Derivate** (1104 ms). Implementiert: `std::vector<SpMatrix> m_dcn` ersetzt durch `CNDerivStore` (Pair-List + (N,3)-Diag). Pattern entlehnt aus GPU-`k_cn_chainrule`. Triplet-Allokation + `setFromTriplets` + `makeCompressed` entfällt komplett. `mixture.xyz` Sweep: CN-deriv T=1→16: 446→141 ms (3.16× Skalierung, vs. 1.31× pre-WP4). T=4 CN-deriv 1104→243 ms (**4.5× schneller**). T=4 Total-Wall 1840→1012 ms (**1.82× schneller**); T=8 Total 1673→730 ms (**2.29×**). Korrektheit bit-identisch auf kleinen Molekülen. Vollständiger Bericht: `docs/wp4/WP4-cn-simd.md`. |
| P4f | ⚙️ Machine-tested | **D4 Gaussian-Weights Inner-Loop-Refactor.** `precomputeGaussianWeights` + `computeGaussianWeightDerivatives` in `d4param_generator.cpp` umgestellt auf Stack-`Eigen::Array<double, MAX_REF=7, 1>` Pipeline — eliminiert ~37k Heap-Allocations pro Energy-Step. Bit-identisch zu pre-WP5. **D4-GW-Wall im Noise** (55.5 ms post-WP5 vs ~56 ms pre-WP5 bei T=4 — 3 Runs für Stabilität). objdump zeigt: libm `exp@plt` wird weiterhin skalar aufgerufen — Eigen `Array<double, 7, 1>::exp()` vektorisiert nicht ohne MKL-VML / vexp-Intrinsic. Total Wall T=4: ~967 → ~931 ms (~4% schneller, vermutlich Heap-Reduktions-Effekt). Vollständiger Bericht: `docs/wp4/WP5-d4-gaussian-simd.md`. |
| P4g | 🆕 Vorgeschlagen | Größter Refactor; nicht vor P4d |
| P4h | ⚙️ Machine-tested | **`m_gradient` + `m_geometry` Layout-Migration** auf RowMajor (`GeoGradMatrix = Matrix<double, Dynamic, 3, RowMajor>`). ~30 Member umgestellt in `forcefieldthread.h/.cpp`, `forcefield.h/.cpp`, `ff_workspace.h/.cpp`, `gfnff.h`, `cuda/ff_workspace_gpu.h/.cu`. Externe API (`getGradient()`) bleibt `Matrix` (ColumnMajor) — Eigen konvertiert am Boundary. ALPB Solvation-Boundary ebenfalls mit Konversion. **Korrektheit T=1 bit-identisch** zu pre-WP-G. T=4 Multi-Thread-Drift auf triose pre-existiert (5 Läufe pre-WP-G zeigen dieselbe Verteilung — eigenständiger Race-Bug, nicht WP-G). **Performance-Gewinn nur 3-7%** auf `mixture.xyz` (Total Wall T=4: 1012 → 967 ms), Bond-cpu-Wall im Run-zu-Run-Noise. Die WP3-Hypothese "Cache-Miss-Lawine in `m_gradient.row()`" hat sich **nicht bestätigt** — moderner Hardware-Prefetcher cached die ColumnMajor-3-Cache-Lines effektiv. Strukturell wertvoll (Konsistenz mit GPU-RowMajor-Erwartung), aber kein dramatischer Speedup. Vollständiger Bericht: `docs/wp4/WP-G-gradient-layout.md`. |
| P4i | 🆕 Vorgeschlagen | **WP-R Multi-Thread-Race-Bug-Hunt.** `triose.xyz` (22 Atome, 1 Fragment) zeigt bei T=4 nicht-deterministische Energie zwischen Läufen — drei stabile Endzustände mit bis zu 1.47 mEh Drift. T=1 deterministisch. Verifiziert pre-WP-G und post-WP-G — bestehende Race, nicht durch jüngere WPs eingeführt. 1.47 mEh ist KEIN Round-Off → algorithmische Inkonsistenz in iterativem Solver (Kandidaten: EEQ-Phase 2, Hückel-SCF, HB/XB-Detection, Param-Gen-Pool, Eigen-internal-Threading). Vorbedingung für saubere Regression-Tests und MD-Stabilität auf großen Systemen. Vollständiger Plan: `docs/wp4/WP-R-multithread-race.md`. |
| P4j | 🆕 Vorgeschlagen | **WP-V Gradient-Validierung CPU vs GPU vs Fortran.** Bestehende, dokumentierte Diskrepanz (siehe `docs/GPU_GFNNF_DISCREPANCIES.md` Mar–Apr 2026): CPU/GPU/XTB-Fortran-Gradienten weichen voneinander ab; Energie ist auf < 1 μEh konsistent, der Drift sitzt im Gradient. MD heat-bath-Exchange auf polymer-System driftet ~0.2 Eh (CPU vs Fortran), GPU näher an Fortran. **Pre-existing**, nicht durch WP4-Reihe eingeführt. Strategie: numerische FD-Referenz aus XTB pro Term, dann Term-Bisect, dann Fix. Hot-Kandidaten: HB-Gradient-Distribution, CN chain-rule Reduktions-Reihenfolge, atomicAdd vs sequential summation, Repulsion (CLAUDE.md sagt "⚠️ Partial"). Vorbedingung: WP-R sollte zuerst aufgeklärt sein. Vollständiger Plan: `docs/wp4/WP-V-gradient-validation.md`. |

## d4_gw im CPU-MD: C6-Refresh einmal, Blockpartition (Okt 2026)

`d4_gw` (Gaussgewichte, C6-Refresh, dC6/dCN) kostete bei polymer_2x (7320 Atome) 103 ms je Schritt bei 16 Threads. Drei exakte Änderungen: (1) der zweite C6-Refresh am Ende von `GFNFF::prepareCNAndEEQ`
entfällt, wenn der Gradientenpfad im selben Aufruf Gewichte und Halbkontraktion für dasselbe `m_last_cn` neu gebaut (`updateCNValuesForGradient` gibt jetzt `true` zurück) und `refreshDispersionC6()` die C6 wirklich
geschrieben hat (gibt `bool` zurück); beim ersten Aufruf (Setup-Geometrie nicht veraltet), bei Energy-only-Aufrufen und wenn der `d4_cn_cache_threshold` die Neuberechnung übersprungen hat, läuft er wie bisher. (2) Zusammenhängende Paarblöcke statt Stride
in `refreshDispersionC6` und in der Paarschleife von `computeDC6DCN` (jede Zelle hat einen Schreiber, Werte unabhängig von der Partition). (3) `m_dc6dcn` wird in parallelen Spaltenblöcken per memset genullt statt seriell (429 MB). Maschinengeprüft:
md5 der `-dump_gradient`-Datei alt/neu identisch in 52 Vergleichen (6 Bulk-Systeme inkl. polymer_2x bei 1 und 16 Threads, 31 GMTKN55, 15 MOR41), `md_diagnostics`-Datei (Ladungen, CN, Energien) identisch bei 50 Schritten (16 Threads) und 4 Schritten (1 Thread).
Zeit polymer_2x MD 20 fs, abwechselnd alt/neu, n=3, Median: 16 Threads `d4_gw` 102.7 -> 58.1 ms, Schritt 382.8 -> 329.3 ms; 1 Thread 165.3 -> 111.1 ms, Schritt 773.1 -> 726.6 ms. Nicht geprüft: GPU/ROCm, periodisch,
`d4_cn_cache_threshold` > 0 im MD, `-gfnff.dispersion_c6_update false`. Offen: Paar-Arrays statt N x N (B), ein Durchlauf für C6 und dc6 (D), C6 erst bei der Energie (E).

## d4_gw: dC6/dCN als Paar-Arrays statt N x N (Kandidat B, Okt 2026)

`D4ParameterGenerator::m_dc6dcn` war eine dichte N x N-Matrix (429 MB bei 7320 Atomen), die jeder Gradientenschritt nullte, obwohl nur die Zellen (i,j)/(j,i) der Paarliste gefüllt und gelesen wurden.
Jetzt `D4DC6DCNPairs` (`gfnff_parameters.h`): sortierte Paarliste (i<j) mit `dc6_di[p]`, `dc6_dj[p]`; `computeDC6DCN` rechnet je Paar mit demselben Kern und derselben Summationsreihenfolge. Leser: `FFWorkspace::calcDispersion`/`calcD4Dispersion`
und `D4Evaluator::pairEnergyAndGradient` über `lookup()` (Positionshinweis, sonst Binärsuche; ein Paar außerhalb der Liste hatte Zelle 0 und trägt nichts bei). Der frühere Voll-N²-Pfad ohne Host-Paarliste (nur GPU, `gpu_disp_pairs_on_device`) hatte keinen Host-Leser und rechnet nichts mehr.
Maschinengeprüft: md5 der `-dump_gradient`-Datei alt/neu identisch in 54 Vergleichen (8 Bulk/Molekül-Fälle inkl. polymer_2x bei 1 und 16 Threads, `-gfnff.dispersion_cutoff_bohr 30` (gefilterte Liste), caffeine `gfn2`; 31 GMTKN55, 15 MOR41), `md_diagnostics`-Datei identisch bei 50 Schritten (16 Threads) und 4 Schritten (1 Thread).
Zeit polymer_2x MD 20 fs, abwechselnd alt/neu, n=3, Median, Last 3-6: 16 Threads `d4_gw` 61.3 -> 33.7 ms, Schritt 374.6 -> 348.1 ms; 1 Thread `d4_gw` 119.9 -> 79.5 ms, Schritt 759.1 -> 724.9 ms (Schritt-Zeiten bei der Last unsicher).
GPU (`-gpu cuda`): Host-`computeDC6DCN` läuft nur beim Setup und beim Neuaufbau der Paarliste, sein Ergebnis liest die GPU nicht (dC6/dCN auf dem Gerät, `k_dc6dcn_per_pair`); `updateCNValuesForGradient`, `refreshDispersionC6`, `prepareCNAndEEQ` laufen im GPU-Pfad nicht (gdb-Aufrufzählung, caffeine SP und MD 5 fs).
Vor A+C / HEAD / B je n=3 auf der GPU: Unterschiede innerhalb der GPU-Streuung (polymer_2x |dE| <= 3.4e-13 Eh bei Streuung 1.1e-13, |dg| <= 3.8e-16 Eh/A; caffeine, hydronium+ mit `frag_charge_model reference` bitgleich; MD 20 fs |dE| <= 6.8e-13 bei Streuung 6.8e-13); GPU-CPU unverändert, GPU-Schritt 38.0-38.9 ms in allen drei. Nicht geprüft: ROCm/Vulkan, `gpu_disp_pairs_on_device`, periodisch, Pfad mit Paarlisten-Neuaufbau auf der GPU.

## GPU-MD: CN-Paarliste als Verlet-Liste, Restart-Zustand nur bei `-rescue` (K1, K6, Okt 2026)

**K1** Ursache: `checkDisplacement(0.5)` vergleicht mit der Topologie-Referenz, die im GPU-Pfad nach dem Setup nie neu gesetzt wird; ab der ersten Auslenkung > 0,5 Bohr (polymer_2x: Schritt 4) war die Prüfung dauerhaft wahr und `generateCNPairListOnGPU()` (O(N^2), 4,75 ms) lief jeden Schritt (CUPTI: 397 Neubauten in 400 Schritten).
Jetzt hat die Liste eine eigene Referenzgeometrie (gesetzt bei jedem Neubau, `FFWorkspaceGPU::cnPairListStale()`), hält Paare bis f*rs + 2s (s = `-gfnff.gpu_cn_pair_skin_bohr`, 1,0) und wird vor `computeCN` neu gebaut, sobald ein Atom sich weiter als D = s + rs_min/2 (f - 1 - 8/7,5) bewegt hat (f = 2,5, H-H: D = 1,35 Bohr).
Begründung: ein fehlendes Paar hat dann r >= rs (1 + 8/7,5), sein erf-CN-Beitrag ist <= 0,5 erfc(8) = 5,6e-30, die Ableitung <= 4,2e-28/Bohr. Die Topologieprüfung selbst ist unverändert; erste Liste wie bisher nach `computeCN` (Single Point unverändert). ROCm: unverändert (Trait `has_cn_pair_skin = false`). `-gfnff.gpu_cn_pair_skin_bohr -1` stellt den alten Zustand her.
Gemessen (polymer_2x, RTX 5080, 16 Threads): 94 Builds in 1000 Schritten; `k_cn_compute_pairs` 9 -> 17 us, `k_cn_chainrule` 12 -> 19 us; ff_total Median 24,0 -> 19,3 ms (n=3 x 198 Schritte). Ergebnisse alt/neu innerhalb der GPU-Streuung (Known Issue #36): SP n=3 |dE| <= 2,3e-13 Eh (Streuung 2,3e-13), |dg| <= 3,6e-16; MD 100 fs Epot gleich auf 1e-6 Eh, Ladungen |dq| <= 8,5e-14 (Streuung 6,6e-14); 1000 fs NVE/CSVR alt/neu im Rahmen der Streuung, Etot-Drift NVE -0,0274 Eh in beiden.
**K6** `SimpleMD::step()` legte jeden Dump `WriteRestartInformation()` in `m_run_states` ab (640 kB bei polymer_2x), gelesen nur im `-rescue`-Zweig; jetzt nur bei `-rescue`, und dort wird der Zustand aus `TriggerWriteRestart()` (gibt ihn jetzt zurück) übernommen statt ein zweites Mal gebaut. Ausgabedateien alt/neu (CPU, 8 und 20 Dumps, mit/ohne `-rescue`) md5-gleich bis auf Zeitstempel/Seed.
Zeit (Wanduhr 400 gegen 100 fs, abwechselnd, n=3): Standard-Dump 38,4 -> 34,0 ms/Schritt, `-dump_frequency 1` 51,6 -> 41,9 ms. Nicht geprüft: ROCm (nicht kompiliert), periodische Systeme, `gpu_cn_pair_cutoff_factor` != 2,5, Optimierung mit Energie-only-Aufrufen.

## GPU-MD: EEQ-Löser und Dispersion parallel zu dc6dcn (K3, K2, Okt 2026)

**Abhängigkeiten (Code gelesen).** `k_dispersion` liest Koordinaten, die Dispersions-SoA (`idx`, `r4r2ij`, `r0_sq`, `r_cut`, `zetac6` aus den Phase-1-Topologieladungen beim Paarlistenbau) und C6/dC6dCN aus `k_dc6dcn_per_pair` (Gaussgewichte der CN); keine EEQ-Ladung. Die frühere Notiz "needs dc6dcn from the post-EEQ charges" war falsch. Der Gerät-EEQ (WP7-E ppcg) liest Koordinaten, `rhs_atoms` (`k_build_eeq_rhs`, Phase 1) und alpha/gam (synchron beim Topologiewechsel); er schreibt nur eigene Puffer. Gemeinsame Schreibziele von Dispersion und Coulomb/Postprocess sind `d_grad`, `d_dEdcn`, `d_energies` (alle `atomicAdd`).
**Ursache der Serialisierung.** Der ppcg wartete über ein Ereignis auf dem Legacy-Stream auf jede vorherige Arbeit aller Streams (auch dc6dcn und Phase-1-Arbeitsströme); `k_dispersion` wurde erst nach der Rückkehr des Lösers auf dem Host eingereiht.
**Umsetzung.** K3: `FFWorkspaceGPU::recordEEQInputsReady()` nach Phase 1, vor dc6dcn; `EEQSolverGPU::setInputReadyEvent()` ersetzt im ppcg die Legacy-Sperre (einmalig je Schritt, andere Löserpfade unverändert); der Löser-Stream bekommt die höchste Priorität (`setStreamPriorityHigh`). K2: `launchDispersionGradientEarly()` startet `k_dispersion` nach dc6dcn auf dem Pairwise-Stream mit 1 Block je SM (Grid-Stride-Schleife; volles Grid unverändert eine Paarberechnung je Thread); Phase 2 lässt den Haupt-Stream vorab auf dessen Ende warten (vor jeder Graph-Aufnahme). Schalter `-gfnff.gpu_eeq_overlap`, `-gfnff.gpu_dispersion_overlap` (true; beide false = alter Ablauf, CUPTI und Wanduhr gleich alt). ROCm: Trait `has_stream_overlap = false`.
**Ohne Priorität und mit vollem Grid** brachte die Überlappung nichts (CUPTI-Median 19,5 ms wie alt): rechengebundene FP64-Kernel füllen das Gerät, die kleinen CG-Kernel warteten dahinter. Grid der Dispersion (CUPTI-Median, n=1 je Wert): voll 18,4, 0,75/SM 16,9, 1/SM 15,9, 1,5/SM 16,1 ms (dc6dcn volles Grid; dc6dcn zusätzlich begrenzt war nicht besser).
**Gemessen** (polymer_2x, RTX 5080, 16 Threads). CUPTI (n=1 Lauf, 24-28 Normalschritte): Schritt-Median alt 19,59, beide aus 19,63, K3 18,79, K3+K2 16,18 ms; Gerät belegt 18,81 / 18,80 / 18,23 / 15,65. CG-Iteration 0,33 -> 0,48 (K3) -> 0,55 ms (K3+K2), weil dc6dcn bzw. Dispersion gleichzeitig laufen; die Schranke -7 ms wird daher nicht erreicht. Wanduhr 400 gegen 100 fs, abwechselnd, n=3, Last 1,2-2,9: alt 34,11, aus 34,18, K3 33,50, K3+K2 30,91 ms/Schritt. K4 (Coulomb auf eigenem Stream) entfällt: Coulomb beginnt 0,1 ms nach dem CG-Ende, die Dispersion ist dann fertig.
**Ergebnisse** innerhalb der GPU-Streuung (Known Issue #36): SP n=3 polymer_2x |dE| alt-alt 2,3e-13, alt-neu 3,4e-13 Eh, |dg| 3,7e-16 / 4,0e-16; caffeine, hydronium-water2 gleich; GPU gegen CPU |dg| 3,2e-14 / 1,5e-13 / 5,3e-14 wie alt. MD 100 fs n=3: Epot gleich auf 1e-6, |dq| alt-alt 5,8e-14, alt-neu 7,5e-14. 400 fs NVE n=2: |dq| 3,6e-11 / 4,1e-11, Etot-Drift -0,0548 Eh in allen vier Läufen. CPU bitgleich (8 SP, MD 8 Schritte). compute-sanitizer racecheck/synccheck/memcheck auf caffeine (5 Schritte, ppcg und cholesky) ohne Befund; Wettläufe über Streams im globalen Speicher prüft keines dieser Werkzeuge.
**Nicht geprüft**: ROCm, Mehr-GPU-Split, andere GPUs und Systemgrößen, Energie-only-Aufrufe mit ppcg, Schrittfolgen mit Löser-Fallback (Jacobi-Wiederholung) unter Überlappung.

## GPU-MD: Host-Neubauten der Paarlisten exakt beschleunigt (Okt 2026)

**Messung vorher** (polymer_2x, `-gpu cuda`, 16 Threads, 1000 fs, Zeitmarken im Hostcode, Instrumentierung wieder entfernt). Alle drei Neubauten laufen am Schrittanfang seriell auf dem Host, das Gerät wartet:
| Neubau | Auslöser | Anzahl in 1000 Schritten | Dauer je Neubau (Median) | im Mittel je Schritt |
|---|---|---|---|---|
| HB/XB-Neuerkennung | erzwungen alle 10 Schritte (`hb_update_force_every`), RMSD 0,3 Bohr | 100 | 39,8 ms (Erkennung 25,5, CPU-Arbeitsbereich 4,8, Hochladen 4,4, `rebuildBondHBData` 3,5) | 3,97 ms |
| Repulsions-Paarliste | Verschiebung > Skin/2 = 1 Bohr | 126 (alle 7-9) | 19,0 ms Host + 5,8 ms Hochladen | 3,14 ms |
| D4-Paarliste (9,4 M Paare) | Verschiebung > 5 Bohr (Skin 10) | 12 (alle 64-99) | 584 ms Host + 125 ms Hochladen | 8,52 ms |
Ursache D4: im MD steht OpenMP prozessweit auf 1 Thread (CxxThreadPool), CN (85 ms) und Paarschleife (450 ms) liefen seriell.
**Änderungen** (alle exakt: gleiche Listen, gleiche Reihenfolge, gleiche Werte):
- `GFNFF::updateDispersionPairsIfNeeded()`: Neubau mit dem Thread-Budget der Instanz (`ScopedBlasThreads`, Pool an den Generator); `D4ParameterGenerator::GenerateDispersionPairsNative()` zweiphasig (Zählen, Schreiben an Zeilen-Offsets) statt Zeilenvektoren plus serieller Verkettung, zeta je Atom und r4r2/R0^2 je Elementpaar aus Tabellen mit identischen Ausdrücken, dC6/dCN-Paarliste im selben Durchlauf.
- `detectHydrogenBondsNative()`: erster Bindungspartner per Tabelle, Bindungstest über Partnerlisten, AB-Paare parallel über Atombereiche auf einer gefilterten Kopie der Zellenliste (`SpatialCellList::subset()`, `forEachPairOfAtom()`), Fall-1-Suche nur über HB-Wasserstoffe, Einträge verschoben statt kopiert; `rebuildBondHBData()` mit Hash-Map.
- `generateRepulsionPairsNative()`: Kandidatensuche parallel über Atombereiche, Bindungstest über Partnerlisten.
- CUDA: `FFWorkspaceGPU::updateDispersion()/updateRepulsion()` laden über einen gepinnten 2 x 32-MB-Puffer, gefüllt von `-threads` Host-Threads, asynchron in Stücken (`-gfnff.gpu_staged_list_upload`, Standard true; false = alter Weg). Debug-Kopie der HB/XB-Liste nur bei Verbosität >= 3.
**Messung nachher** (gleiche Auslöseschritte): HB 12,8 ms, Repulsion 8,4 + 2,8 ms, D4 76,7 + 48,7 ms; im Mittel 15,6 -> 4,3 ms je Schritt. Wanduhr (1000 gegen 250 fs, abwechselnd, n=3, Last 1,4-2,7): 33,6/33,8/33,7 -> 21,7/21,7/21,7 ms je Schritt. `-md_diagnostics_timing`, dump 1, 1000 Schritte: step_total Mittel 40,3 -> 29,0, Median 23,8 -> 23,8 ms.
**Prüfung**: CPU md5 gleich (SP: 6 Bulk inkl. polymer_2x 1/16 Threads, 31 GMTKN55 inkl. geladen/offenschalig, 15 MOR41, 5 GFN2; MD 8 Schritte 16 Threads Standard und mit Neubau in jedem Schritt, 4 Threads). GPU: Gerätelisten (Hash) alt gleich neu an 44 Neubaupunkten (200 fs) und 31 erzwungenen Neubauten; Host gleich Gerät an 33 D4- und 56 Repulsions-Neubauten; SP/MD-Streuung alt-neu wie alt-alt; 400 fs NVE Etot-Drift -0,0548 Eh unverändert.
**Nicht geprüft**: ROCm (nicht kompiliert, Trait `has_staged_list_upload = false`), Mehr-GPU-Split, periodische Systeme, andere Systemgrößen.
