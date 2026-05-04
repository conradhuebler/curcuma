# WP3: Pair-List-basierter CN-Kernel — O(N²) → O(N·k) Reduktion

**Kategorie**: Mittlerer Aufwand  
**Aufwand**: ~6–8 Stunden  
**Wirkung**: Hoch — k_cn_compute wird jeden Schritt aufgerufen; 5–10× weniger Arbeit für N=1410  
**Abhängigkeiten**: Keine (unabhängig von WP1/WP2)  
**Status**: 🤖 Geplant

---

## Kontext

`k_cn_compute` weist jeden Atom einen Thread zu, der alle anderen N-1 Atome durchläuft und den erf()-basierten CN-Beitrag berechnet. Für N=1410 sind das ~2 Millionen erf()-Auswertungen pro Schritt, von denen ~80% den Distanz-Cutoff überschreiten und nichts beitragen.

Die Infrastruktur für GPU-Pair-Lists existiert bereits (`k_generate_cn_pairs_count/write` in `gfnff_kernels.cuh`), wird aber aktuell nur für `k_cn_chainrule` verwendet.

**Plan**: 
1. CN-Pair-List einmalig bei Topologie-Build auf GPU generieren (mit CUB-Prefix-Sum statt atomicCounter)
2. Neuer Kernel `k_cn_compute_pairs`: 1 Thread/Paar, atomicAdd auf `d_cn_raw[i/j]`
3. `k_dlogdcn` unverändert danach

---

## Schritt 1: CUB-basierte Pair-List-Generierung

### Problem mit dem aktuellen Ansatz
`k_generate_cn_pairs_count` zählt gültige Paare pro Atom mit atomicAdd auf einem globalen Counter. Das führt zu Contention bei vielen Threads.

### CUB-Lösung (contentionfrei)

**Dateien**: `ff_workspace_gpu.cu`, `CMakeLists.txt`

```cpp
// CUB ist Header-only und im CUDA-Toolkit enthalten:
#include <cub/cub.cuh>

// Neue Methode in FFWorkspaceGPU: buildCNPairList(stream)
void FFWorkspaceGPU::buildCNPairListGPU(cudaStream_t stream)
{
    auto& impl = *m_impl;
    const int N = impl.natoms;

    // Phase 1: Anzahl gültiger Paare pro Atom zählen (1 Thread/Atom)
    DeviceBuffer<int> d_counts(N + 1);  // +1 für Prefix-Sum
    k_generate_cn_pairs_count<<<gridFor(N), 256, 0, stream>>>(
        N, impl.coords.d_x.ptr, impl.coords.d_y.ptr, impl.coords.d_z.ptr,
        impl.atom_types.ptr, impl.cn_cutoff_factor,
        d_counts.ptr  // d_count[i] = Anzahl Paare für Atom i
    );

    // Phase 2: Exclusive Prefix Sum via CUB (O(N), contentionfrei)
    size_t temp_bytes = 0;
    cub::DeviceScan::ExclusiveSum(nullptr, temp_bytes, d_counts.ptr, d_counts.ptr, N + 1, stream);
    DeviceBuffer<char> d_temp(temp_bytes);
    cub::DeviceScan::ExclusiveSum(d_temp.ptr, temp_bytes, d_counts.ptr, d_counts.ptr, N + 1, stream);
    // Nach ExclusiveSum: d_counts[i] = start-Index für Atom i, d_counts[N] = total_pairs

    // Total-Paar-Anzahl aus Device runterladen (Pinned-Buffer für Async)
    int total_pairs = 0;
    cudaMemcpyAsync(&total_pairs, d_counts.ptr + N, sizeof(int),
                    cudaMemcpyDeviceToHost, stream);
    cudaStreamSynchronize(stream);  // Einmalig pro Topologie-Build — akzeptabel

    // Buffers (re-)allozieren wenn nötig
    if (total_pairs > impl.cn_pairs.capacity) {
        impl.cn_pairs.d_idx_i.reallocate(total_pairs);
        impl.cn_pairs.d_idx_j.reallocate(total_pairs);
        impl.cn_pairs.d_rcov_sum.reallocate(total_pairs);
        impl.cn_pairs.capacity = total_pairs;
    }
    impl.cn_pairs.n = total_pairs;

    // Phase 3: Paare in vorberechnete Positionen schreiben (kein atomicAdd nötig)
    k_generate_cn_pairs_write_indexed<<<gridFor(N), 256, 0, stream>>>(
        N, impl.coords.d_x.ptr, impl.coords.d_y.ptr, impl.coords.d_z.ptr,
        impl.atom_types.ptr, impl.cn_cutoff_factor,
        d_counts.ptr,                    // Prefix-Sum-Offsets
        impl.cn_pairs.d_idx_i.ptr,
        impl.cn_pairs.d_idx_j.ptr,
        impl.cn_pairs.d_rcov_sum.ptr
    );
}
```

### Neues k_generate_cn_pairs_write_indexed Kernel

Der bestehende `k_generate_cn_pairs_write` benutzt einen atomicAdd-Counter. Die neue Version nutzt den CUB-Prefix-Sum-Offset:

```cpp
// gfnff_kernels.cuh:
/// Write CN pair list using precomputed per-atom offsets (contention-free).
/// d_offsets[i] = start index for atom i (from CUB ExclusiveSum on d_counts).
__global__ GFNFF_KERNEL_BOUNDS_LIGHT void k_generate_cn_pairs_write_indexed(
    int N,
    const double* __restrict__ cx, const double* __restrict__ cy, const double* __restrict__ cz,
    const int*    __restrict__ atom_types,
    double        cutoff_factor,
    const int*    __restrict__ d_offsets,  ///< [N] start index per atom (CUB prefix sum)
    int*          idx_i, int* idx_j, double* rcov_sum
);

// gfnff_kernels.cu — Implementierung:
__global__ GFNFF_KERNEL_BOUNDS_LIGHT void k_generate_cn_pairs_write_indexed(...)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;
    int zi = atom_types[i];
    double xi = cx[i], yi = cy[i], zi_coord = cz[i];
    double rcov_i = d_rcov_d3[zi - 1];  // constant memory

    int write_pos = d_offsets[i];  // Startposition für Atom i
    for (int j = i + 1; j < N; ++j) {
        double dx = xi - cx[j];
        double dy = yi - cy[j];
        double dz = zi_coord - cz[j];
        double r2 = dx*dx + dy*dy + dz*dz;
        double rcov_ij = rcov_i + d_rcov_d3[atom_types[j] - 1];
        double thresh = cutoff_factor * rcov_ij;
        if (r2 < thresh * thresh) {
            idx_i[write_pos] = i;
            idx_j[write_pos] = j;
            rcov_sum[write_pos] = rcov_ij;
            ++write_pos;
        }
    }
}
```

**Hinweis**: `k_generate_cn_pairs_count` muss analog (nur Zählen, kein Schreiben) angepasst werden, um die Counts pro Atom (nicht global) zurückzugeben.

---

## Schritt 2: Neuer k_cn_compute_pairs Kernel

**Deklaration in `gfnff_kernels.cuh`**:
```cpp
/// Pair-list-based CN computation: 1 thread per (i,j) pair, atomicAdd to d_cn_raw.
/// Replaces O(N^2) k_cn_compute for cached topology (pair list reused across steps).
/// Reference: gfnff_cn.f90:66-126
__global__ GFNFF_KERNEL_BOUNDS_LIGHT void k_cn_compute_pairs(
    int n_pairs,
    const int*    __restrict__ idx_i,     ///< [n_pairs] atom i indices
    const int*    __restrict__ idx_j,     ///< [n_pairs] atom j indices
    const double* __restrict__ rcov_sum,  ///< [n_pairs] sum of covalent radii
    const double* __restrict__ cx, const double* __restrict__ cy, const double* __restrict__ cz,
    double* cn_raw,                        ///< [N] output: erf-sum (atomicAdd, zero beforehand)
    double  kn                             ///< decay constant (-7.5)
);
```

**Implementierung in `gfnff_kernels.cu`**:
```cpp
__global__ GFNFF_KERNEL_BOUNDS_LIGHT void k_cn_compute_pairs(
    int n_pairs,
    const int* __restrict__ idx_i, const int* __restrict__ idx_j,
    const double* __restrict__ rcov_sum,
    const double* __restrict__ cx, const double* __restrict__ cy, const double* __restrict__ cz,
    double* cn_raw, double kn)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= n_pairs) return;

    int i = idx_i[tid], j = idx_j[tid];
    double dx = cx[i]-cx[j], dy = cy[i]-cy[j], dz = cz[i]-cz[j];
    double r = __dsqrt_rn(dx*dx + dy*dy + dz*dz);

    // erf-basierter CN-Beitrag: 0.5 * erfc(kn * (r/rcov_sum - 1))
    // erfc(x) = 1 - erf(x); hier: 0.5*(1 + erf(-kn*(r/rcov-1))) = 0.5*erfc(kn*(r/rcov-1))
    double arg = kn * (r / rcov_sum[tid] - 1.0);
    double contrib = 0.5 * erfc(arg);  // CUDA hat keine native erfc für double, aber erf():
    // contrib = 0.5 * (1.0 - erf(arg));

    // Symmetrisch: Beitrag für beide Atome i und j
    atomicAdd(&cn_raw[i], contrib);
    atomicAdd(&cn_raw[j], contrib);
}
```

**Wichtig**: Vor dem Kernel muss `d_cn_raw` auf 0 gesetzt werden (`k_zero_double<<<gridFor(N), 256>>>`).

---

## Schritt 3: Integration in ff_workspace_gpu.cu

### Struct-Erweiterung (ff_workspace_gpu.h):

```cpp
struct CNPairList {
    DeviceBuffer<int>    d_idx_i;
    DeviceBuffer<int>    d_idx_j;
    DeviceBuffer<double> d_rcov_sum;
    int n = 0;           ///< Aktuelle Paarzahl
    int capacity = 0;    ///< Allozierte Kapazität
    bool valid = false;  ///< Initialisiert?
};
CNPairList cn_pairs;  ///< In FFWorkspaceGPUImpl hinzufügen
```

### Anpassung in computeCNOnGPU():

```cpp
void FFWorkspaceGPU::computeCNOnGPU(cudaStream_t stream, double kn, double cnmax, double threshold_sq)
{
    auto& impl = *m_impl;
    const int N = impl.natoms;

    // WP3: Pair-List-Pfad (wenn verfügbar)
    if (impl.cn_pairs.valid && impl.cn_pairs.n > 0) {
        // cn_raw nullen
        auto cfg_n = getLaunchConfig(N);
        k_zero_double<<<cfg_n.gridSize, cfg_n.blockSize, 0, stream>>>(
            impl.cn.d_cn_raw.ptr, N);

        // Pair-List-CN
        auto cfg = getLaunchConfig(impl.cn_pairs.n);
        k_cn_compute_pairs<<<cfg.gridSize, cfg.blockSize, 0, stream>>>(
            impl.cn_pairs.n,
            impl.cn_pairs.d_idx_i.ptr, impl.cn_pairs.d_idx_j.ptr,
            impl.cn_pairs.d_rcov_sum.ptr,
            impl.coords.d_x.ptr, impl.coords.d_y.ptr, impl.coords.d_z.ptr,
            impl.cn.d_cn_raw.ptr, kn
        );
    } else {
        // Fallback: Original O(N²) Kernel
        auto cfg = getLaunchConfig(N);
        k_cn_compute<<<cfg.gridSize, cfg.blockSize, 0, stream>>>(
            N, ..., kn, cnmax, threshold_sq);
    }

    // log-Transformation unverändert:
    {
        auto cfg = getLaunchConfig(N);
        k_dlogdcn_and_logcn<<<cfg.gridSize, cfg.blockSize, 0, stream>>>(...);
    }
}
```

### buildCNPairListGPU() aufrufen

In `uploadTopologyToGPU()` bzw. nach `updateReferenceGeometry()`:

```cpp
// Nach dem Upload aller SoA-Strukturen bei Init/Topologie-Rebuild:
m_gpu_workspace->buildCNPairListGPU(stream);
```

---

## Abschätzung der Paaranzahl

Für GFN-FF-Moleküle mit typischer Dichte (1 Å³/Atom):
- Covalent-Radii-Cutoff: ca. 2× Summe der kovalenten Radii ≈ 2×(1.5 Bohr × 2) = 6 Bohr = 3.2 Å
- Mit cutoff_factor ≈ 2.0: Cutoff ≈ 6 Å
- Durchschnittlich ~30 Nachbarn/Atom für Polymer

→ N=1410 Atome × 30 Paare/2 (symmetrisch) ≈ 21.000 Paare.  
Statt 2.000.000 Iterationen → **~95× weniger Arbeit**.

---

## Verifikation

```bash
# Build:
cd release && make -j4

# CN-Korrektheit:
./curcuma -sp ../test_cases/molecules/larger/polymer.xyz -method gfnff -gpu cuda -verbosity 3
# → CN-Werte müssen identisch zu CPU-Pfad sein (in der Ausgabe vergleichen)

# Energie-Regression:
./curcuma -sp ../test_cases/molecules/larger/polymer.xyz -method gfnff -gpu cuda
./curcuma -sp ../test_cases/molecules/larger/polymer.xyz -method gfnff
# → Differenz < 1 nEh

# CTest:
ctest -R "cli_gfnff\|cli_simplemd" --output-on-failure

# Timing (k_cn_compute sollte deutlich kürzer sein):
./curcuma -md ../test_cases/cli/simplemd/10_gfnff_polymer_md/input.xyz -method gfnff -gpu cuda -steps 50 -verbosity 1
```

---

## Sonderfälle

- **n_pairs = 0** (z.B. Einzelatom): Fallback auf `k_cn_compute` (Guard `if (n > 0)`)
- **Topologie-Rebuild**: Pair-List wird neu generiert; `cn_pairs.valid = false` während Rebuild
- **PBC**: `k_generate_cn_pairs_write_indexed` muss Minimum-Image-Convention berücksichtigen — analog zu bestehenden `applyMIC()`-Aufrufen in anderen Kerneln
