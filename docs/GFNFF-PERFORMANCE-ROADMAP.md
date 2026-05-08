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
| P4f | 🆕 Vorgeschlagen | — |
| P4g | 🆕 Vorgeschlagen | Größter Refactor; nicht vor P4d |
| P4h | 🆕 Vorgeschlagen | **`m_gradient` Layout-Optimierung** (RowMajor oder Per-Thread-Buffer). Nach WP3 hinzugefügt: Heap-Allocation war nicht der Bond-Hotspot, vermutete echte Ursache ist Cache-Miss-Lawine durch ColumnMajor-`m_gradient` mit Stride 49.6 KB für N=6200. Eine Layout-Änderung würde **alle** FF-Terme treffen (Bond, Coulomb, Dispersion, Repulsion, ATM, H-Bonds). Erwartung: Bond cpu-time 170 → ~50 ms, Coulomb 660 → ~250 ms; Total Wall mixture.xyz T=4 von 1840 → ~1500 ms (≥300 ms Ersparnis). Aufwand 4-8h, mehrere Aufrufer betroffen. Vollständiger Plan: `docs/wp4/WP-G-gradient-layout.md`. |
