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

---

## Nicht empfohlen

- OpenMP/SIMD-Pragmas CPU: CxxThreadPool + Eigen-SIMD reichen. Kein messbarer Mehrwert.
- Cutoff-Reduktion: Physikalisch riskant ohne Referenz-Validierung.
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
| G1a | ⏳ offen | — |
| G1b | ⏳ offen | — |
| G2a | ⏳ offen | — |
| G2b | ⏳ offen | — |

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
