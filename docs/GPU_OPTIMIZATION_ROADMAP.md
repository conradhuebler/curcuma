# GFN-FF GPU-Optimierungs-Roadmap — RTX 5080

**Stand**: Mai 2026 | **Hardware**: AMD Ryzen 9 9950X3D + RTX 5080 (Blackwell SM_120, 64 MB L2)  
**Gemessene Baseline** (1000 MD-Schritte, Polymer N=1410):

| Methode | ms/Schritt | Bemerkung |
|---------|-----------|-----------|
| gfnff-gpu (aktuell, 85ca3cf clean) | **~29 ms** | nach `m_cn_pairs_generated`-Bugfix |
| gfnff-gpu (mit G2a-Experiment) | ~28.7 ms | Gaussian-Weights vor Phase 1, deaktiviert |
| gfnff-gpu (historisch ccefa22) | ~19 ms | Referenz-Baseline (Ziel wiederherstellen) |
| gfnff CPU (4 Threads) | 148.8 ms | — |
| xtb-gfnff | 422.1 ms | — |

Siehe [WP1](GPU_WP1_L2_PERSISTENCE_LAUNCH_BOUNDS.md) für Analyse der Ursachen und Experimente.

GPU aktuell **7.8× schneller als CPU**. Details: [GFNFF-GPU.md](GFNFF-GPU.md)

**Hinweis**: Die ptxas-Register-Daten in GFNFF-GPU.md wurden auf **sm_75 (Turing)** gemessen — für SM_120 (Blackwell) gelten andere Register-Budgets. WP1 kalibriert gezielt für SM_120.

**Ziel**: 5–10× GPU-Speedup vs. Baseline (≤6 ms/Schritt, N=1410)

---

## Überblick: Warum ist die GPU aktuell untergenutzt?

```
Pro Schritt (aktuell):
  GPU: k_cn_compute [O(N²)]
  GPU→CPU: finalizeCNForCPU()  ← Sync-Punkt 1: GPU wartet auf CPU
  CPU: prepareEEQParametersForGPU() [O(N) loop]
  CPU: eeq_gpu.solve() [Matrix-Build + Cholesky O(N³)]
  CPU→GPU: setEEQCharges() H2D
  GPU: launchChargeDependentAndFinish() [7 Kernel-Launches, kein Graph]
  GPU→CPU: download energy + gradient ← Sync-Punkt 2
```

GROMACS-Vergleich: keine CPU-Sync-Punkte, kein O(N³) im Hot-Path, >80% GPU-Auslastung.

---

## Arbeitspakete

| WP | Titel | Aufwand | Wirkung | Abhängig von |
|----|-------|---------|---------|--------------|
| [WP1](GPU_WP1_L2_PERSISTENCE_LAUNCH_BOUNDS.md) | L2-Persistenz + Blackwell Launch-Bounds | ~4h | Mittel–Hoch | — |
| [WP2](GPU_WP2_EEQ_RHS_KERNEL.md) | k_build_eeq_rhs — GPU-seitiger RHS | ✅ **Implementiert** | Infrastruktur (kein Timing-Gewinn) | — |
| [WP3](GPU_WP3_PAIRLIST_CN_KERNEL.md) | Pair-List CN — O(N²)→O(N·k) | ~8h | **Hoch — nächster Schritt** | — |
| [WP4](GPU_WP4_PHASE2_CUDA_GRAPH.md) | Phase-2-CUDA-Graph | ~6h | Mittel | WP2 ✅ |
| [WP5](GPU_WP5_FULL_GPU_EEQ_PIPELINE.md) | Vollständig GPU-residenter EEQ | ~10 Tage | **Transformativ** | WP2 ✅ + WP3 + WP4 |
| [WP6](GPU_WP6_TILED_CN_SHARED_MEMORY.md) | Tiled k_cn_compute (Shared Memory) | ~7 Tage | Hoch (N≥3000) | WP3 optional |

---

## Empfohlene Reihenfolge

```
WP2 ✅ (EEQ RHS Kernel — Infrastruktur für WP4/WP5)
   ↓
WP3 ← NÄCHSTER SCHRITT (CN Pair-List, O(N²)→O(N·k), ~3–5 ms Gewinn)
   ↓
WP1 (L2-Persistenz + Launch-Bounds, unabhängig, Quick Win)
   ↓
WP4 (Phase-2-Graph, baut auf WP2 ✅ auf)
   ↓
WP5 (vollständig GPU, baut auf WP2 ✅ + WP3 + WP4 auf)
   ↓
WP6 (Shared Memory CN, für N≥3000)
```

---

## Bottleneck-Analyse

### Aktuelle Pipeline (per Schritt, N=1410, RTX 5080)

| Schritt | Zeit (geschätzt) | Bottleneck | WP |
|---------|-----------------|------------|----|
| GPU Phase 1 (CUDA Graph) | ~3–5 ms | Kernel-Compute (O(N²) erf) | WP3, WP6 |
| `finalizeCNForCPU` Sync | ~0.1 ms | Stream-Sync | WP5 |
| `prepareEEQParametersForGPU` | ~0 ms | **Entfällt** (WP2 ✅) | — |
| `eeq_gpu.solveWithDeviceRHS()` fresh | ~10–15 ms | **O(N³/6) Cholesky ← Haupt-Bottleneck** | WP5 |
| `eeq_gpu.solveWithDeviceRHS()` lazy | ~1–2 ms | O(N²) potrs | WP4 |
| `setEEQCharges` H2D | ~0.05 ms | 11 KB Upload | WP5 |
| GPU Phase 2 (7 Kernel-Launches) | ~3–5 ms | Kernel-Compute | WP4 |
| D2H energy+gradient | ~0.1 ms | 11 KB Download | WP5 |

**Dominierender Faktor**: EEQ Cholesky O(N³/6) (~10–15 ms fresh solve, ~50% des Schritts).  
WP2 hat die O(N) CPU-Vorbereitung (~0.2 ms) eliminiert — im Rauschen.  
Realer Gewinn kommt von **WP3** (O(N²) CN-Kernel → O(N·k)) und **WP5** (Cholesky aus hot path).

### Nach WP3 (nächster Schritt)

| Schritt | Zeit vorher | Zeit nachher |
|---------|------------|-------------|
| `k_cn_compute` (N=1410) | ~3–5 ms | ~0.3–0.6 ms (5–10×) |
| **Gesamtschritt** | ~29 ms | ~25–27 ms |

CN ist derzeit O(N²) mit erf(): 2 Mio. Aufrufe pro Schritt für N=1410, davon ~80% außerhalb Cutoff. Pair-List reduziert auf ~400K effektive Paare.

### Nach allen WPs

| Schritt | Zeit (erwartet) | Verbesserung |
|---------|----------------|-------------|
| Gesamt GPU-Schritt | ~3–5 ms | Alle CPU-Syncs weg |
| EEQ (lazy, GPU-seitig) | ~0.5–1 ms | Kein H2D/D2H |
| CN (Pair-List oder Tiling) | ~0.2–0.5 ms | 5–10× schneller |
| Phase 2 (CUDA Graph) | ~2–3 ms | Zero-Overhead Launch |

---

## Schlüssel-Dateien

| Datei | Relevant für |
|-------|-------------|
| `src/core/energy_calculators/qm_methods/gfnff_gpu_method.cpp` | Hot-Path Orchestration (WP2, WP5) |
| `src/core/energy_calculators/ff_methods/cuda/ff_workspace_gpu.cu` | Graph-Capture, CN, L2 (WP1, WP3, WP4) |
| `src/core/energy_calculators/ff_methods/cuda/gfnff_kernels.cu` | Kernel-Implementierungen (WP2, WP3, WP6) |
| `src/core/energy_calculators/ff_methods/cuda/gfnff_kernels.cuh` | Kernel-Deklarationen, launch_bounds (WP1) |
| `src/core/energy_calculators/ff_methods/cuda/eeq_solver_gpu.cu` | EEQ Cholesky + Schur (WP5) |
| `src/core/energy_calculators/ff_methods/gfnff.h` | EEQGPUParams-Erweiterung (WP2) |
| `src/core/energy_calculators/ff_methods/gfnff_method.cpp` | prepareEEQParametersForGPU (WP2) |
