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
| [WP2](GPU_WP2_EEQ_RHS_KERNEL.md) | k_build_eeq_rhs — GPU-seitiger RHS | ~8h | **Hoch** | — |
| [WP3](GPU_WP3_PAIRLIST_CN_KERNEL.md) | Pair-List CN — O(N²)→O(N·k) | ~8h | **Hoch** | — |
| [WP4](GPU_WP4_PHASE2_CUDA_GRAPH.md) | Phase-2-CUDA-Graph | ~6h | Mittel | WP2 empfohlen |
| [WP5](GPU_WP5_FULL_GPU_EEQ_PIPELINE.md) | Vollständig GPU-residenter EEQ | ~10 Tage | **Transformativ** | WP2+WP3+WP4 |
| [WP6](GPU_WP6_TILED_CN_SHARED_MEMORY.md) | Tiled k_cn_compute (Shared Memory) | ~7 Tage | Hoch (N≥3000) | WP3 optional |

---

## Empfohlene Reihenfolge

```
WP1 (Quick Win, unabhängig)
   ↓
WP3 (CN Pair-List, unabhängig)
   ↓
WP2 (EEQ RHS Kernel, unabhängig)
   ↓
WP4 (Phase-2-Graph, baut auf WP2 auf)
   ↓
WP5 (vollständig GPU, baut auf WP2+WP3+WP4 auf)
   ↓
WP6 (Shared Memory CN, für große N)
```

---

## Bottleneck-Analyse

### Aktuelle Pipeline (per Schritt, N=1410, RTX 5080)

| Schritt | Zeit (geschätzt) | Bottleneck |
|---------|-----------------|------------|
| GPU Phase 1 (CUDA Graph) | ~3–5 ms | Kernel-Compute |
| `finalizeCNForCPU` Sync | ~0.1 ms | Stream-Sync |
| `prepareEEQParametersForGPU` | ~0.2 ms | CPU O(N) |
| `eeq_gpu.solve()` (fresh) | ~10–15 ms | O(N³/6) Cholesky |
| `eeq_gpu.solve()` (lazy, RMSD>0) | ~1–2 ms | O(N²) potrs |
| `setEEQCharges` H2D | ~0.05 ms | 11 KB Upload |
| GPU Phase 2 (7 Kernel-Launches) | ~3–5 ms | Kernel-Compute |
| D2H energy+gradient | ~0.1 ms | 11 KB Download |

**Dominierender Faktor**: EEQ Cholesky (~10–15 ms für fresh solve). WP2+WP5 adressieren dies durch GPU-seitige RHS-Konstruktion und vollständig asynchrone EEQ-Pipeline.

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
