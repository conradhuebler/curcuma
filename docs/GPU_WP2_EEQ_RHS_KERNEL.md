# WP2: k_build_eeq_rhs — GPU-seitiger EEQ-RHS-Aufbau

**Kategorie**: Mittlerer Aufwand  
**Aufwand**: ~6–8 Stunden  
**Wirkung**: Hoch — eliminiert `finalizeCNForCPU`-Sync + O(N) CPU-Loop aus dem kritischen Pfad  
**Abhängigkeiten**: WP1 empfohlen, aber unabhängig implementierbar  
**Status**: ✅ Implementiert (Mai 2026) — kein messbarer Timing-Gewinn, siehe Diagnose unten

---

## Kontext

Jeder GFN-FF-Schritt führt aktuell folgende Sequenz aus:

```
GPU: k_cn_compute → d_cn_final
GPU→CPU: finalizeCNForCPU() [Sync-Punkt 1]
CPU: prepareCNAndEEQ()  → CN-Derivate (O(N))
CPU: prepareEEQParametersForGPU(cn) → rhs_atoms[i] = chi_corrected[i] + cnf[i]*sqrt(max(cn[i],0))
CPU→GPU: eeq_gpu.solve() → H2D alpha_corrected, gam_corrected, rhs_atoms
```

Das `rhs_atoms`-Array ist der einzige CN-abhängige Schritt auf CPU. `chi_corrected[i]` und `cnf[i]` sind **topologie-konstant** (berechnet einmalig aus `topo.eeq_chi`, `topo.dxi`, `topo.eeq_cnf`). Nur der `sqrt(cn[i])`-Term ändert sich jeden Schritt.

Ziel: `k_build_eeq_rhs<<<ceil(N/256), 256>>>(d_cn_final, d_chi_corrected, d_cnf, d_rhs, N)` direkt nach `k_cn_compute` auf demselben Stream starten → kein CPU-Sync nötig.

---

## Neue Struct: EEQTopologyParams (GPU-Upload, einmalig)

Topologie-konstante EEQ-Parameter als GPU-Buffer in `ff_workspace_gpu.h`:

```cpp
// In FFWorkspaceGPUImpl (ff_workspace_gpu.h, Impl-Struct):
struct EEQTopologyBuffers {
    DeviceBuffer<double> d_alpha_corrected; ///< [N] alpha^2 pro Atom (topologiekonstant)
    DeviceBuffer<double> d_gam_corrected;   ///< [N] gamma + dgam (topologiekonstant)
    DeviceBuffer<double> d_chi_corrected;   ///< [N] -chi + dxi + amide-Korrektur (topologiekonstant)
    DeviceBuffer<double> d_cnf;             ///< [N] cnf_eeq pro Atom (topologiekonstant)
    DeviceBuffer<double> d_rhs_atoms;       ///< [N] RHS-Vektor (per-step, GPU-beschreibbar)
    DeviceBuffer<double> d_rhs_constraints; ///< [nfrag] Ziellladungen pro Fragment
    int                  nfrag = 1;
    bool                 valid = false;     ///< Initialisiert nach erstem Topologie-Upload
};
EEQTopologyBuffers eeq_topo; ///< hinzufügen in FFWorkspaceGPUImpl
```

---

## Neue Methode: uploadEEQTopologyParams()

Wird nach jeder Topologie-Berechnung aufgerufen (d.h. einmal bei Init und bei Topologie-Rebuild).

**Deklaration in `ff_workspace_gpu.h`**:
```cpp
/// Upload topology-constant EEQ parameters to GPU once per topology build.
/// After this, k_build_eeq_rhs can construct rhs_atoms[i] entirely on GPU.
void uploadEEQTopologyParams(const GFNFF::EEQGPUParams& params, cudaStream_t stream = 0);
```

**Implementierung in `ff_workspace_gpu.cu`**:
```cpp
void FFWorkspaceGPU::uploadEEQTopologyParams(const GFNFF::EEQGPUParams& params, cudaStream_t stream)
{
    auto& impl = *m_impl;
    const int N = impl.natoms;
    auto& et = impl.eeq_topo;

    // alpha_corrected und gam_corrected: konstant, kein sqrt(cn) drin
    et.d_alpha_corrected.uploadAsync(params.alpha_corrected.data(), N, stream);
    et.d_gam_corrected.uploadAsync(params.gam_corrected.data(), N, stream);

    // chi_corrected: rhs_atoms ohne den cnf*sqrt(cn)-Term
    // rhs_atoms[i] = chi_corrected[i] + cnf[i]*sqrt(cn[i])
    // → chi_corrected[i] = rhs_atoms[i] - cnf[i]*sqrt(cn[i])
    // ABER: wir brauchen chi_corrected und cnf getrennt für den GPU-Kernel.
    // Also: neue Methode getChiCorrected() und getCNF() in GFNFF erzeugen
    // ODER: EEQGPUParams um zwei neue Felder erweitern.
    et.d_chi_corrected.uploadAsync(params.chi_corrected_static.data(), N, stream);
    et.d_cnf.uploadAsync(params.cnf.data(), N, stream);

    // rhs_constraints (nfrag doubles) — klein, einfach hochladen
    et.d_rhs_constraints.uploadAsync(params.rhs_constraints.data(), params.nfrag, stream);
    et.nfrag = params.nfrag;
    et.valid = true;
}
```

### Erweiterung von EEQGPUParams in gfnff.h

```cpp
struct EEQGPUParams {
    std::vector<double> alpha_corrected;      // [N] — wie bisher
    std::vector<double> gam_corrected;        // [N] — wie bisher
    std::vector<double> rhs_atoms;            // [N] — wie bisher (CPU-Pfad)
    std::vector<double> rhs_constraints;      // [nfrag]
    std::vector<int>    fraglist;             // [N]
    int nfrag = 1;
    // NEU: Aufspaltung von rhs_atoms für GPU-Kernel
    std::vector<double> chi_corrected_static; // [N] = -chi + dxi + amide_corr (topologiekonstant)
    std::vector<double> cnf;                  // [N] cnf_eeq (topologiekonstant)
};
```

**Anpassung in `prepareEEQParametersForGPU()` in `gfnff_method.cpp`**:

```cpp
// Zeile 1029-1037: Aufspaltung der rhs_atoms-Berechnung
double chi_corrected = -chi_base + dxi_i;
if (z == 1 && i < static_cast<int>(topo.is_amide_h.size()) && topo.is_amide_h[i]) {
    chi_corrected -= 0.02;
}
// NEU: getrennte Speicherung
params.chi_corrected_static[i] = chi_corrected;  // topologiekonstant
params.cnf[i] = cnf_i;                            // topologiekonstant
// Alt bleibt für CPU-Pfad:
params.rhs_atoms[i] = chi_corrected + cnf_i * std::sqrt(std::max(cn_i, 0.0));
```

---

## Neuer Kernel: k_build_eeq_rhs

**Deklaration in `gfnff_kernels.cuh`**:
```cpp
/// Construct EEQ RHS vector entirely on GPU: rhs[i] = chi_corrected[i] + cnf[i]*sqrt(max(cn[i],0))
/// Run directly after k_cn_compute on same stream — eliminates finalizeCNForCPU() sync.
/// Reference: gfnff_method.cpp:prepareEEQParametersForGPU:1015-1037
__global__ GFNFF_KERNEL_BOUNDS_LIGHT void k_build_eeq_rhs(
    int N,
    const double* __restrict__ d_cn,          ///< [N] GPU CN (d_cn_final nach k_cn_compute)
    const double* __restrict__ d_chi_corr,    ///< [N] -chi+dxi+amide_corr (topologiekonstant)
    const double* __restrict__ d_cnf,         ///< [N] cnf_eeq (topologiekonstant)
    double*       __restrict__ d_rhs          ///< [N] output RHS
);
```

**Implementierung in `gfnff_kernels.cu`**:
```cpp
__global__ GFNFF_KERNEL_BOUNDS_LIGHT void k_build_eeq_rhs(
    int N,
    const double* __restrict__ d_cn,
    const double* __restrict__ d_chi_corr,
    const double* __restrict__ d_cnf,
    double*       __restrict__ d_rhs)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;
    double cn_i = d_cn[i];
    double sqrt_cn = (cn_i > 0.0) ? __dsqrt_rn(cn_i) : 0.0;  // fused sqrt, kein Branch-Overhead
    d_rhs[i] = d_chi_corr[i] + d_cnf[i] * sqrt_cn;
}
```

---

## Integration in den Per-Schritt-Pfad

### Anpassung in `ff_workspace_gpu.cu`: `prepareAndLaunchChargeIndependent()`

Nach `k_cn_compute` und `k_dlogdcn` direkt anfügen:

```cpp
// WP2: EEQ RHS auf GPU aufbauen (nach k_cn_compute, vor EEQ solve)
if (impl.eeq_topo.valid) {
    auto cfg = getLaunchConfig(impl.natoms);
    k_build_eeq_rhs<<<cfg.gridSize, cfg.blockSize, 0, stream>>>(
        impl.natoms,
        impl.cn.d_cn_final.ptr,          // d_cn nach log-Transformation
        impl.eeq_topo.d_chi_corrected.ptr,
        impl.eeq_topo.d_cnf.ptr,
        impl.eeq_topo.d_rhs_atoms.ptr    // Ausgabe: device-seitiger RHS
    );
    CUDA_CHECK_ASYNC(stream);
}
```

### Neue Methode: `getDeviceRHSPtr()` in `ff_workspace_gpu.h`

```cpp
/// Return device pointer to EEQ RHS vector (valid after k_build_eeq_rhs ran)
const double* getDeviceRHSPtr() const;
```

### Anpassung in `eeq_solver_gpu.cu`: `solve()` akzeptiert Device-Pointer

Neue Überladung (Rückwärtskompatibilität bleibt):

```cpp
/// Overload: rhs_atoms already on GPU (from k_build_eeq_rhs).
/// Eliminates H2D upload of rhs_atoms per step.
bool solveWithDeviceRHS(
    int natoms, int nfrag,
    const double* cx, const double* cy, const double* cz,  // device ptrs
    const double* d_alpha_corrected,     // device ptr (uploadEEQTopologyParams)
    const double* d_gam_corrected,       // device ptr
    const double* d_rhs_atoms,          // device ptr (k_build_eeq_rhs output)
    const double* d_rhs_constraints,    // device ptr
    const std::vector<int>& fraglist,   // host (nur für Schur-Komplement)
    double* out_z1,                     // host output (für CPU-Schur)
    double* out_Z2,
    bool force_refactor = true
);
```

### Anpassung in `gfnff_gpu_method.cpp`

In `calculateEnergy()`, Zeilen 407–504 ersetzen:

```cpp
// VORHER (Zeilen 407-504):
m_gpu_workspace->finalizeCNForCPU(m_gpu_cn_final);          // GPU→CPU Sync
m_gfnff->prepareCNAndEEQ(gradient, true, &m_gpu_cn_final, true);
GFNFF::EEQGPUParams eeq_params = m_gfnff->prepareEEQParametersForGPU(cn);
// ... eeq_gpu.solve() mit host-seitigem rhs_atoms

// NACHHER (WP2):
// k_build_eeq_rhs läuft bereits auf main stream (gestartet in prepareAndLaunchChargeIndependent)
// prepareCNAndEEQ nur noch für CN-Derivate (skip_eeq=true, kein rhs_atoms mehr nötig):
if (gradient) {
    m_gpu_workspace->finalizeCNForCPU(m_gpu_cn_final);       // Noch nötig für CNF-Derivate
    m_gfnff->prepareCNAndEEQ(gradient, true, &m_gpu_cn_final, /*skip_eeq=*/true);
}
// EEQ solve mit device-seitiger RHS:
bool eeq_ok = m_eeq_gpu->solveWithDeviceRHS(
    N, 1,
    m_gpu_workspace->getDeviceXPtr(),
    m_gpu_workspace->getDeviceYPtr(),
    m_gpu_workspace->getDeviceZPtr(),
    m_gpu_workspace->getDeviceAlphaPtr(),   // via uploadEEQTopologyParams
    m_gpu_workspace->getDeviceGamPtr(),
    m_gpu_workspace->getDeviceRHSPtr(),     // k_build_eeq_rhs output
    m_gpu_workspace->getDeviceRHSConstraintsPtr(),
    eeq_fraglist,
    m_eeq_z1.data(), m_eeq_Z2.data(),
    force_refactor);
```

---

## Topologie-Update-Handling

`uploadEEQTopologyParams()` muss nach jedem Topologie-Rebuild aufgerufen werden:

```cpp
// In gfnff_gpu_method.cpp, nach consumeFullTopologyUpdate():
if (m_gfnff->consumeFullTopologyUpdate()) {
    m_gpu_workspace->updateReferenceGeometry();
    m_gpu_workspace->invalidateGraph();
    m_eeq_has_ref_geom = false;

    // WP2: Topologie-Parameter neu hochladen
    const Vector& dummy_cn = m_gfnff->getLastCN();
    GFNFF::EEQGPUParams topo_params = m_gfnff->prepareEEQParametersForGPU(dummy_cn);
    m_gpu_workspace->uploadEEQTopologyParams(topo_params);
}
```

Und einmalig bei `setMolecule()` nach `InitialiseMolecule()`.

---

## Diagnose: Warum kein Timing-Gewinn?

WP2 eliminiert O(N) CPU-Arbeit (~0.2 ms) und 3×N H2D-Uploads (~0.1 ms). Der dominierende Bottleneck ist der EEQ-Cholesky (`eeq_gpu.solve()`), der ~10–15 ms kostet (fresh) bzw. ~1–2 ms (lazy). Ein 0.3 ms Gewinn auf einem 29 ms Schritt ist im Rauschen.

**Was WP2 tatsächlich liefert**: Infrastruktur für WP5 (vollständig GPU-residenter EEQ). `d_rhs_atoms` liegt nun als stabiler Device-Pointer vor — WP4 (Phase-2-Graph) kann ihn direkt in den Graphen einbetten ohne H2D-Dependency.

**Was WP2 NICHT implementiert hat**: `finalizeCNForCPU()` bleibt bedingungslos (sollte laut Spec auf `if (gradient)` beschränkt werden). Erfordert Restructurierung der Fallback-Pfade — kleiner Gewinn (~0.1 ms), vorerst zurückgestellt.

## Implementierte Änderungen (Mai 2026)

| Datei | Änderung |
|-------|---------|
| `gfnff.h` | `EEQGPUParams` um `chi_corrected_static`, `cnf` erweitert |
| `gfnff_method.cpp` | `prepareEEQParametersForGPU()` füllt neue Felder |
| `gfnff_kernels.cuh` | `GFNFF_KERNEL_BOUNDS_LIGHT`-Makro + `k_build_eeq_rhs` Deklaration |
| `gfnff_kernels.cu` | `k_build_eeq_rhs` Implementierung |
| `ff_workspace_gpu.h` | `uploadEEQTopologyParams()`, `isEEQTopoValid()`, 4 Device-Ptr-Getter |
| `ff_workspace_gpu.cu` | `EEQTopologyBuffers` in Impl; `uploadEEQTopologyParams()`; Kernel-Launch in `prepareAndLaunchChargeIndependent()` |
| `eeq_solver_gpu.h` | `solveWithDeviceRHS()` Überladung deklariert |
| `eeq_solver_gpu.cu` | `solveWithDeviceRHS()` implementiert (D2D für col0, H2D nur für Constraint-Spalten) |
| `gfnff_gpu_method.h` | `m_eeq_fraglist`, `m_eeq_rhs_constraints`, `m_eeq_nfrag` gecacht |
| `gfnff_gpu_method.cpp` | Topo-Upload bei init + rebuild; `solve()` → `solveWithDeviceRHS()` im Haupt-Pfad |

## Verifikation

```bash
# Build:
cd release && make -j4

# Energie-Regression (GPU vs CPU muss <1 nEh bleiben):
./curcuma -sp ../test_cases/molecules/larger/polymer.xyz -method gfnff -gpu cuda
./curcuma -sp ../test_cases/molecules/larger/polymer.xyz -method gfnff

# Gradient-Test:
./curcuma -sp ../test_cases/molecules/larger/polymer.xyz -method gfnff -gpu cuda -gradient

# CTest:
ctest -R "cli_gfnff\|cli_simplemd" --output-on-failure

# Timing-Vergleich (vorher/nachher):
./curcuma -md ../test_cases/molecules/larger/polymer.xyz -method gfnff -gpu cuda -steps 100 -verbosity 1
# Erwartung: "EEQ" Timing reduziert, "CPU EEQ sync" verschwindet aus Ausgabe
```

---

## Nächste Schritte nach WP2

WP4 (Phase-2-Graph) wird nach WP2 einfacher, da `d_rhs_atoms` in fester Device-Adresse liegt.
