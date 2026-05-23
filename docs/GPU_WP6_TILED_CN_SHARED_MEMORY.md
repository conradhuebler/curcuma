# WP6: Tiled k_cn_compute mit Shared Memory

**Kategorie**: Großer Aufwand  
**Aufwand**: ~5–8 Tage  
**Wirkung**: Hoch für k_cn_compute (5–20× schneller), besonders relevant wenn N > 5000  
**Abhängigkeiten**: WP3 kann alternativ implementiert werden; WP6 ersetzt WP3 für große N  
**Status**: 🤖 Geplant — niedrigere Priorität als WP1–WP5

---

## Kontext

WP3 (Pair-List CN) reduziert die Arbeit von O(N²) auf O(N·k) durch Eliminierung von Paaren außerhalb des Cutoffs. WP6 adressiert ein anderes Problem: **Speicherbandbreite**. Auch mit Pair-List muss jeder Thread Koordinaten für zwei Atome aus Global Memory lesen. Für N=1410 und k=30 Nachbarn: ~42.000 Paare × 6 doubles (2 Atome × 3 Koordinaten) × 8 Byte = 2 MB Bandbreite pro CN-Schritt.

Auf Blackwell SM_120 hat jede SM **128 KB L1/Shared Memory**. Für N=1410 passen alle Koordinaten (3×1410×8 = 34 KB) in eine einzige SM. Mit Tiling kann jede Koordinate einmal geladen und von allen Threads verwendet werden.

Das entspricht dem klassischen **N-Body-Tiling** aus der GPU-Astrophysik (vgl. Nyland et al., GPU Gems 3, Kapitel 31).

---

## Algorithmus: Tiled O(N²) CN-Kernel

```
Grid: ceil(N/TILE_I) × ceil(N/TILE_J) Blöcke  (N×N-Zerlegung in TILE×TILE-Kacheln)
Jeder Block: TILE_I × TILE_J Threads
Shared Memory: TILE_I × 3 doubles (Atom-i-Koordinaten) + TILE_J × 3 doubles (Atom-j-Koordinaten)
```

**Warum nicht Pair-List + Shared Memory?** Die Pair-List ermöglicht sparse Zugriff auf Koordinaten — die Zugriffsmuster sind nicht mehr tile-freundlich. Das Tiling funktioniert am besten mit dem dichten O(N²) Algorithmus, weil dann die Zugriffe vorhersehbar und blockweise sind.

---

## Implementierung

### Parameter

```cpp
// gfnff_kernels.cuh:
#define CN_TILE_SIZE 64  // 64×64 = 4096 Threads/Block für große N
                          // 64×3×8 × 2 = 3 KB Shared Memory (weit unter 128 KB)
```

**Warum 64?** 
- Shared Memory: 2 × 64 × 3 × 8 = 3072 Bytes → passt locker in 128 KB
- 64² = 4096 Threads/Block → zu groß. Besser: jeder Thread bearbeitet eine (i,j)-Kachel.
- Realistisch: TILE=64 mit TILE×1 Threads (64 Threads/Block, jeder loopt über 64 j-Atome)

**Alternative**: Klassisches GPU-N-Body-Tiling mit square blocks:
- Block: TILE×1 Threads, jeder Thread = ein i-Atom
- Shared Memory: Tile der j-Koordinaten (TILE×3 doubles = 1.5 KB)
- Jeder Thread loopt über alle j in der Kachel → TILE j-Atome pro Iteration

```cpp
// Optimales Layout für Blackwell:
#define CN_TILE 128  // 128 j-Atome pro Kachel
// Shared Memory: 128 × 3 × 8 = 3 KB (trivial)
// Grid: N × ceil(N/TILE) Blöcke (jeder Block = 1 i-Atom, TILE j-Atome)
// Threads/Block: 128 (1 Thread/j-Atom in der Kachel)
```

### Kernel-Deklaration in gfnff_kernels.cuh

```cpp
/// Tiled CN computation using shared memory coordinate cache.
/// Each block processes TILE j-atoms for one or more i-atoms.
/// Reduces global memory bandwidth by TILE× for j-coordinate reads.
/// Reference: GPU Gems 3, Chapter 31 (N-body tiling); gfnff_cn.f90:66-126
/// Requires: __shared__ cache for j-coordinates (3*TILE*8 bytes per block)
__global__ void k_cn_compute_tiled(
    int natoms,
    const double* __restrict__ cx,          ///< [N] SoA x-coordinates (Bohr)
    const double* __restrict__ cy,          ///< [N] SoA y-coordinates (Bohr)
    const double* __restrict__ cz,          ///< [N] SoA z-coordinates (Bohr)
    const int*    __restrict__ atom_types,  ///< [N] 1-based atomic numbers
    double*       cn_raw,                   ///< [N] output (atomicAdd)
    double        kn,                       ///< decay constant (-7.5)
    double        cnmax,
    double        threshold_sq
);
```

### Kernel-Implementierung in gfnff_kernels.cu

```cpp
#define CN_TILE 128

__global__ void k_cn_compute_tiled(
    int N,
    const double* __restrict__ cx, const double* __restrict__ cy, const double* __restrict__ cz,
    const int* __restrict__ atom_types,
    double* cn_raw,
    double kn, double cnmax, double threshold_sq)
{
    // Shared memory: j-Atom-Koordinaten für diesen Tile
    __shared__ double s_cx[CN_TILE], s_cy[CN_TILE], s_cz[CN_TILE];
    __shared__ int    s_zt[CN_TILE];  // Atom-Typen für rcov-Lookup

    int i = blockIdx.x;  // 1 Block = 1 i-Atom (oder mehrere mit Stride)
    if (i >= N) return;

    double xi = cx[i], yi = cy[i], zi = cz[i];
    int    zi_type = atom_types[i];
    double rcov_i = d_rcov_d3[zi_type - 1];  // constant memory

    double cn_sum = 0.0;  // Lokale Akkumulation für Atom i

    // Über alle j-Tiles iterieren
    int num_tiles = (N + CN_TILE - 1) / CN_TILE;
    for (int tile = 0; tile < num_tiles; ++tile) {
        // j-Koordinaten kollaborativ in Shared Memory laden (1 Thread lädt 1 j-Atom)
        int j_global = tile * CN_TILE + threadIdx.x;
        if (j_global < N) {
            s_cx[threadIdx.x] = cx[j_global];
            s_cy[threadIdx.x] = cy[j_global];
            s_cz[threadIdx.x] = cz[j_global];
            s_zt[threadIdx.x] = atom_types[j_global];
        }
        __syncthreads();

        // Jeder Thread (repräsentiert Atom i) berechnet Beiträge aller j im Tile
        int tile_end = min(CN_TILE, N - tile * CN_TILE);
        for (int t = 0; t < tile_end; ++t) {
            int j = tile * CN_TILE + t;
            if (j == i) continue;  // Diagonale überspringen

            double dx = xi - s_cx[t];
            double dy = yi - s_cy[t];
            double dz = zi - s_cz[t];
            double r2 = dx*dx + dy*dy + dz*dz;

            // Cutoff-Check (threshold_sq ist der quadratische Cutoff)
            if (r2 > threshold_sq) continue;

            double r = __dsqrt_rn(r2);
            double rcov_ij = rcov_i + d_rcov_d3[s_zt[t] - 1];
            double arg = kn * (r / rcov_ij - 1.0);
            cn_sum += 0.5 * (1.0 - erf(arg));  // erfc via 1-erf
        }
        __syncthreads();
    }

    // Atomares Schreiben: Atom i hat seinen CN-Beitrag vollständig berechnet
    // (kein atomicAdd nötig — jedes i hat seinen eigenen Block)
    cn_raw[i] = cn_sum;  // Kein atomicAdd! (i ist eindeutig pro Block)
}
```

**Wichtiger Unterschied zu k_cn_compute_pairs (WP3)**: Beim Tiled-Kernel schreibt jeder Block exklusiv auf `cn_raw[i]` — **kein atomicAdd** nötig. Das ist effizienter als der WP3-Ansatz (der atomicAdd auf `cn_raw[i]` und `cn_raw[j]` verwendet).

---

## Grid-Konfiguration

```cpp
// In computeCNOnGPU() — Tiled-Variante:
void FFWorkspaceGPU::computeCNOnGPU_tiled(cudaStream_t stream, ...)
{
    const int N = impl.natoms;
    // 1 Block pro i-Atom, CN_TILE Threads pro Block
    dim3 grid(N);
    dim3 block(CN_TILE);
    size_t shared_bytes = CN_TILE * (3 * sizeof(double) + sizeof(int));

    // cn_raw nullen ist nicht nötig (direktes Schreiben, kein atomicAdd)
    k_cn_compute_tiled<<<grid, block, shared_bytes, stream>>>(
        N, impl.coords.d_x.ptr, impl.coords.d_y.ptr, impl.coords.d_z.ptr,
        impl.atom_types.ptr, impl.cn.d_cn_raw.ptr,
        kn, cnmax, threshold_sq
    );
    // k_dlogdcn unverändert danach
}
```

---

## WP3 vs WP6: Wann was?

| | WP3 (Pair-List) | WP6 (Tiling) |
|---|---|---|
| **Algorithmus** | Sparse O(N·k) | Dense O(N²), aber mit Datenwiederverwertung |
| **atomicAdd** | Ja (auf cn_raw[i] + cn_raw[j]) | Nein (exklusiver Schreib-Zugriff) |
| **Shared Memory** | Nein | Ja (~1.5 KB/Block) |
| **Cutoff-Effizienz** | ~95% nutzlose Paare entfernt | Cutoff-Check im Kernel (aber Koordinaten schon im SM) |
| **Ideal für** | Kleine N (≤2000) | Große N (≥3000) |
| **Graph-Kompatibilität** | Gut | Gut |
| **Code-Komplexität** | Mittel | Mittel |

**Empfehlung**: WP3 für aktuelle Molekülgrößen (N≤2000), WP6 für zukünftige große Systeme (N≥3000, Polymere, Proteine).

**Hybrid-Ansatz**: In `computeCNOnGPU()` nach N wählen:
```cpp
if (N > 3000) {
    computeCNOnGPU_tiled(stream, ...);   // WP6: Shared Memory
} else if (impl.cn_pairs.valid) {
    computeCNOnGPU_pairlist(stream, ...); // WP3: Pair-List
} else {
    computeCNOnGPU_original(stream, ...); // Original O(N²) Fallback
}
```

---

## Blackwell Thread Block Clusters (SM_120, fortgeschritten)

Blackwell SM_120 unterstützt **Thread Block Clusters** bis zu 8 Blöcke. Mehrere `k_cn_compute_tiled`-Blöcke können koordiniert werden, um ihre Shared-Memory-Tiles via **Distributed Shared Memory** auszutauschen — effektiv ein 8×-größeres Tile ohne L2-Traffic.

```cpp
// Für N=1410 mit CN_TILE=128: 1410/128 = 11 Tiles
// Mit Cluster-Größe 8: 8 Blöcke teilen ihre j-Tiles direkt
// Effektive Tile-Größe: 8 × 128 = 1024 j-Atome → 1 Tile deckt 73% von N=1410 ab
// → Fast alle j-Koordinaten aus Distributed Shared Memory statt L2

// CUDA 12.0+ API:
cudaLaunchConfig_t config = {};
cudaLaunchAttribute attr[1];
attr[0].id = cudaLaunchAttributeClusterDimension;
attr[0].val.clusterDim = {8, 1, 1};  // 8 Blöcke pro Cluster
config.attrs = attr;
config.numAttrs = 1;
config.gridDim = N;
config.blockDim = CN_TILE;
config.sharedMemBytes = shared_bytes;
cudaLaunchKernelEx(&config, k_cn_compute_tiled_clustered, ...);
```

Diese Optimierung ist Blackwell-spezifisch und sollte **nach** der grundlegenden Tiling-Implementierung als separate Phase implementiert werden.

---

## Verifikation

```bash
# Build:
cd release && make -j4

# Korrektheit (CN-Werte identisch zu Original):
./curcuma -sp ../test_cases/molecules/larger/polymer.xyz -method gfnff -gpu cuda -verbosity 3
# → CN-Werte in Ausgabe müssen bit-identisch sein (gleicher erf-Aufruf)

# Energie-Regression:
./curcuma -sp ../test_cases/molecules/larger/complex.xyz -method gfnff -gpu cuda
./curcuma -sp ../test_cases/molecules/larger/complex.xyz -method gfnff

# Skalierung testen (N-Abhängigkeit des Speedups):
for mol in H2O.xyz CH4.xyz benzene.xyz polymer.xyz; do
    echo "=== $mol ==="
    time ./curcuma -sp ../test_cases/molecules/larger/$mol -method gfnff -gpu cuda -verbosity 0
done

# nvprof/nsys zur Validierung der Shared-Memory-Hit-Rate:
nsys profile ./curcuma -md polymer.xyz -method gfnff -gpu cuda -steps 20
# → Ziel: k_cn_compute_tiled zeigt hohe L1/Shared-Memory-Hit-Rate in nsys
```
