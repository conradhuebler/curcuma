# WP1: L2-Cache-Persistenz + Blackwell Launch-Bounds-Kalibrierung

**Kategorie**: Quick Win  
**Aufwand**: ~4 Stunden  
**Wirkung**: Mittel–Hoch (15–30% Reduktion L2-Miss-Rate, Occupancy-Steigerung auf SM_120)  
**Abhängigkeiten**: keine  
**Status**: ❌ BLOCKIERT — beide Implementierungsversuche führten zu Energie=0 (Mai 2026)

---

## ❌ Implementierungsversuche (Mai 2026) — beide gescheitert

### Teilaufgabe A: L2-Persistenzfenster (`cudaStreamSetAttribute`)

**Implementiert**: `setGeometry()` in `ff_workspace_gpu.cu` nach den drei `cudaMemcpyAsync`-Aufrufen.

**Ergebnis**: Energie=0, Gradient=0 nach dem ersten Schritt.

**Ursache**: Die drei SoA-Arrays `d_x`, `d_y`, `d_z` werden einzeln via `cudaMalloc` allokiert. Ihr kombinerter Adress-Span übersteigt potenziell das gerätespezifische Persisting-L2-Limit. `cudaStreamSetAttribute` gibt dann `cudaErrorInvalidValue` zurück, was den Stream in einen Fehlerzustand versetzt und alle nachfolgenden Kernel-Starts lautlos unterdrückt.

Versuchte Gegenmaßnahmen (alle gescheitert):
- `num_bytes` auf `l2_bytes` geclampt → Energie weiterhin 0
- `cudaGetLastError()` nach fehlgeschlagenem `cudaStreamSetAttribute` → Energie weiterhin 0

**Voraussetzung für korrekten Fix**: `d_x`, `d_y`, `d_z` in einer einzigen `cudaMalloc(3*N*sizeof(double))`-Allokation mit Offset-Pointern zusammenlegen. Dann ist `num_bytes = 3*N*8 ≈ 34 KB` garantiert kleiner als das L2-Limit.

---

### Teilaufgabe B: Launch-Bounds `GFNFF_KERNEL_BOUNDS_LIGHT` (`__launch_bounds__(256, 4)`)

**Implementiert**: Kernel-Deklarationen in `gfnff_kernels.cuh` von `GFNFF_KERNEL_BOUNDS` auf `GFNFF_KERNEL_BOUNDS_LIGHT` umgestellt für alle als „leicht" klassifizierten Kernel.

**Ergebnis**: Energie=0, Gradient=0.

**Mechanismus**: nvcc propagiert `__launch_bounds__` aus Forward-Deklarationen (`.cuh`) auf Kernel-Definitionen (`.cu`), auch wenn die Definition selbst keine Bounds-Annotation hat. Eine Änderung in der Deklaration ist damit funktional äquivalent zur Änderung der Definition.

**Ursache unklar**: `__launch_bounds__(256, 4)` reduziert das maximale Register-Budget von 128 auf 64 Register/Thread. Mindestens ein Kernel überschreitet dieses Budget auf SM_120, was zu stillem Fehlverhalten führt (kein Launch-Error, nur falsche Ergebnisse: 0). Welcher Kernel betroffen ist, wurde noch nicht isoliert — das erfordert schrittweises Umstellen einzelner Kernels mit Test nach jedem Schritt.

**Gesichert funktionierend** (bereits vor WP1 auf LIGHT):
- `k_build_eeq_rhs`, `k_cn_compute_pairs`, `k_logcn` — diese drei WP2/WP3-Kernels laufen mit `GFNFF_KERNEL_BOUNDS_LIGHT` korrekt.

**Voraussetzung für erneuten Versuch**:
1. `ptxas -v` im Build aktivieren, um echte Register-Zahlen auf SM_120 zu messen
2. Kernels einzeln auf LIGHT umstellen und nach jedem Schritt testen
3. Nur Kernels mit bestätigten ≤64 Registern auf SM_120 auf LIGHT setzen

---

## Durchgeführte Analysen und Experimente (Mai 2026)

### Bug-Fix: `m_cn_pairs_generated` (behoben)

**Problem**: `generateCNPairListOnGPU()` enthält 2× `cudaStreamSynchronize` und wurde **jeden MD-Schritt** aufgerufen, weil `m_cn_pairs_generated = true` nach dem GPU-Aufruf fehlte. Die CPU-Version `generateCNPairList()` setzt das Flag korrekt, wird aber im GPU-Pfad nie aufgerufen.

**Fix** (`gfnff_gpu_method.cpp`):
```cpp
if (!m_cn_pairs_generated) {
    m_gpu_workspace->generateCNPairListOnGPU();
    m_cn_pairs_generated = true;  // war fehlend
}
```

**Einsparung**: ~4–5 ms/Schritt (gemessen).

---

### Experiment G1a: Dispersion-Pair-Sorting (ausprobiert, verworfen)

**Idee**: Dispersionspaar-Liste nach `idx_i` sortieren → koaleszierte Koordinaten-Lesezugriffe für `k_dispersion`.

**Ergebnis**: Kein messbarer Vorteil. Code auskommentiert in `ff_workspace_gpu.cu` (Konstruktor).

---

### Experiment G2a: Gaussian-Weights vor Phase 1 (ausprobiert, derzeit deaktiviert)

**Idee**: `computeGaussianWeightsOnGPU(use_cn_final=true)` vor `prepareAndLaunchChargeIndependent()` aufrufen, sodass `k_dispersion` bereits in Phase 1 (anstatt Phase 2) laufen kann. Phase 1 läuft parallel zur CPU-EEQ → `k_dispersion` überschneidet sich mit dem Cholesky-Solve.

**Messungen** (N=1410, 1000 Schritte):
- Gaussian-Weights-Experiment aktiv: **~25 ms/Schritt** (ca. 1 ms schneller als aktuell)
- Deaktiviert (aktueller Stand): **~26 ms/Schritt**
- Differenz: **~1 Sekunde über 1000 Schritte**

**Warum deaktiviert**: Auf Anfrage zurückgesetzt, um einen sauberen Ausgangspunkt zu haben. Kann reaktiviert werden — die Logik ist korrekt (`event_upload`-Fence in `prepareAndLaunchChargeIndependent` stellt sicher, dass `sA` auf `dc6dcn` wartet).

**Code-Änderungen (wenn reaktiviert)**:
- `gfnff_gpu_method.cpp`: `computeGaussianWeightsOnGPU(use_cn_final=true)` **vor** `prepareAndLaunchChargeIndependent()`
- `ff_workspace_gpu.cu`: `if (m_dispersion_enabled)` statt `if (m_dispersion_enabled && !gradient)` in Phase 1; Phase-2-Dispersion entfernen

---

### CSVR-Thermostat: Analyse `m_atom_temp`

**Problem (identifiziert, dead code)**: In `SimpleMD::CSVR()` (simplemd.cpp) wird in der Atom-Schleife `m_atom_temp[i].push_back(...)` aufgerufen — für N=1410 sind das 1410 Heap-Operationen pro MD-Schritt. `m_atom_temp` wird **nirgends gelesen**. Nach 10.000 Schritten: ~113 MB verschwendeter Speicher.

**Hypothese (widerlegt)**: Der Dead-Code könnte ~13 ms/Schritt kosten.

**Messung**: Thermostat `none` vs. `csvr` zeigt **keine signifikante Zeitdifferenz** → `m_atom_temp` ist **nicht** der Flaschenhals. Der Dead Code bleibt vorerst erhalten (Aufräumen separat).

**Beibehaltene Fixes in `CSVR()`**:
- Guard gegen Division durch null: `if (m_Ekin <= 0.0) return;`
- `m_seed++` kommentiert: Die `static`-Generatoren nutzen `m_seed` nicht → das Inkrement hatte nie einen Effekt.

---

### Aktuelle Performance-Zahlen (Mai 2026, N=1410, RTX 5080, aus baseline.dat)

| Konfiguration | ms/Schritt | Bemerkung |
|---------------|-----------|-----------|
| Historische Referenz (ccefa22) | 19.0 ms | Ziel wiederherstellen |
| 85ca3cf, clean state | 29.2 ms | aktueller Stand ohne Experimente |
| G2a: k_dispersion Phase 1 | 28.9 ms | +0.3 ms Einsparung |
| G2a + event_cn_ready | 28.7 ms | +0.5 ms Einsparung gesamt |
| Ziel (nach WP1–WP3) | ≤15 ms | geschätzt |

**Hauptflaschenhals** (unverändert): CPU-Cholesky EEQ ~10–15 ms/Schritt (dominiert Hot-Path). WP2+WP5 adressieren dies.

**Regression gegenüber Referenz**: ~10 ms/Schritt — Ursache unklar, nicht vollständig auf bekannte Experimente zurückführbar.

---

---

## Kontext

Der RTX 5080 nutzt die Blackwell-Architektur (SM_120) mit 64 MB L2-Cache und 256 KB Registerfile/SM. Die aktuellen `__launch_bounds__(512, 2)` wurden auf Turing (sm_75, 64 KB Registerfile) optimiert. Auf Blackwell sind für register-leichte Kernel (`k_repulsion`, `k_bonds`, `k_coulomb`) höhere Occupancy-Werte möglich. Zusätzlich werden Koordinaten (34 KB für N=1410) von allen 22 Kerneln gelesen, aber nicht explizit im L2 gehalten.

---

## Teilaufgabe A: L2-Cache-Persistenz für Koordinaten

### Dateien
- `src/core/energy_calculators/ff_methods/cuda/ff_workspace_gpu.cu`
- `src/core/energy_calculators/ff_methods/cuda/ff_workspace_gpu.h`

### Implementierung

**Wo einfügen**: In `FFWorkspaceGPU::setGeometry()` direkt nach dem `cudaMemcpyAsync` der Koordinaten.

```cpp
// Nach: cudaMemcpyAsync(impl.coords.d_x.ptr, ...) etc.

// WP1: L2-Persistenzfenster für Koordinaten setzen (Blackwell: 64 MB L2)
// Koordinaten (3*N*8 = 34 KB für N=1410) von allen 22 Kerneln gelesen.
// Ziel: Koordinaten bleiben zwischen Kernel-Wellen im L2.
#if CUDART_VERSION >= 11000
{
    // Alle drei SoA-Koordinaten-Arrays als persisting markieren
    auto set_l2_persist = [&](const double* ptr, size_t bytes) {
        cudaStreamAttrValue attr = {};
        attr.accessPolicyWindow.base_ptr   = const_cast<double*>(ptr);
        attr.accessPolicyWindow.num_bytes  = bytes;
        attr.accessPolicyWindow.hitRatio   = 1.0f;
        attr.accessPolicyWindow.hitProp    = cudaAccessPropertyPersisting;
        attr.accessPolicyWindow.missProp   = cudaAccessPropertyStreaming;
        cudaStreamSetAttribute(stream, cudaStreamAttributeAccessPolicyWindow, &attr);
    };
    size_t coord_bytes = static_cast<size_t>(N) * sizeof(double);
    set_l2_persist(impl.coords.d_x.ptr, coord_bytes);
    set_l2_persist(impl.coords.d_y.ptr, coord_bytes);
    set_l2_persist(impl.coords.d_z.ptr, coord_bytes);
}
#endif
```

**Zurücksetzen**: In `FFWorkspaceGPU::downloadResults()` oder am Ende des Schritts (optional, da dasselbe Fenster in jedem Schritt überschrieben wird).

**Voraussetzung prüfen**: `cudaDeviceGetAttribute(&l2_size, cudaDevAttrL2CacheSize, dev)`. Nur aktivieren wenn L2 ≥ 4 * 3 * N * 8 Bytes (Faktor 4 für etwas Puffer). Für N=1410 sind das 164 KB, problemlos auf RTX 5080 mit 64 MB.

---

## Teilaufgabe B: Launch-Bounds für SM_120 kalibrieren

### Dateien
- `src/core/energy_calculators/ff_methods/cuda/gfnff_kernels.cuh`
- `src/core/energy_calculators/ff_methods/cuda/ff_workspace_gpu.cu` (getLaunchConfig)

### Strategie

Auf Turing/Ampere war `__launch_bounds__(512, 2)` der Kompromiss zwischen Registerdruck und Occupancy. Auf Blackwell mit 256 KB Registerfile/SM sind 4 Blöcke/SM realistisch für leichte Kernel.

**Schritt 1**: ptxas-Output für SM_120 prüfen (nach Build mit `-arch=sm_120`):
```bash
cd release && make -j4 2>&1 | grep "ptxas info.*k_repulsion\|k_bonds\|k_coulomb\|k_dispersion"
```

**Schritt 2**: Kernel-Kategorien nach Registerverbrauch aufteilen.

Aktuelle Messungen (sm_75):
| Kernel | Register | Stack | Status |
|--------|----------|-------|--------|
| k_dispersion | 43 | 0 | leicht |
| k_repulsion | 40 | 0 | leicht |
| k_bonds | ~50 | 0 | leicht |
| k_coulomb | ~55 | 0 | leicht |
| k_coulomb_self | ~30 | 0 | sehr leicht |
| k_angles | 82 | 40 B | mittel |
| k_storsions | 106 | 40 B | schwer |
| k_inversions | 116 | 40 B | schwer |
| k_dihedrals | 128 | 88 B | sehr schwer |
| k_hbonds | 128 | 312+548 B | maximum |

**Schritt 3**: Zwei Makros definieren statt einem:

```cpp
// gfnff_kernels.cuh — WP1: Blackwell-spezifische launch_bounds

// Leichte Kernel: ≤64 Register — 4 Blöcke/SM möglich auf SM_120
#define GFNFF_KERNEL_BOUNDS_LIGHT  __launch_bounds__(256, 4)

// Schwere Kernel: >100 Register — 2 Blöcke/SM wie bisher
#define GFNFF_KERNEL_BOUNDS_HEAVY  __launch_bounds__(256, 2)

// Rückwärtskompatibilität (behält altes Verhalten für unveränderte Kernel)
#define GFNFF_KERNEL_BOUNDS  GFNFF_KERNEL_BOUNDS_HEAVY
```

**Schritt 4**: Makros in Kernel-Signaturen anpassen:

```cpp
// Leichte Kernel — GFNFF_KERNEL_BOUNDS_LIGHT:
__global__ GFNFF_KERNEL_BOUNDS_LIGHT void k_dispersion(...)
__global__ GFNFF_KERNEL_BOUNDS_LIGHT void k_repulsion(...)
__global__ GFNFF_KERNEL_BOUNDS_LIGHT void k_repulsion_mixed(...)
__global__ GFNFF_KERNEL_BOUNDS_LIGHT void k_bonds(...)
__global__ GFNFF_KERNEL_BOUNDS_LIGHT void k_coulomb(...)
__global__ GFNFF_KERNEL_BOUNDS_LIGHT void k_coulomb_self(...)
__global__ GFNFF_KERNEL_BOUNDS_LIGHT void k_coulomb_postprocess(...)
__global__ GFNFF_KERNEL_BOUNDS_LIGHT void k_subtract_qtmp(...)
__global__ GFNFF_KERNEL_BOUNDS_LIGHT void k_cn_chainrule(...)
__global__ GFNFF_KERNEL_BOUNDS_LIGHT void k_dlogdcn(...)
__global__ GFNFF_KERNEL_BOUNDS_LIGHT void k_zero_double(...)

// Schwere Kernel — GFNFF_KERNEL_BOUNDS_HEAVY (wie bisher):
__global__ GFNFF_KERNEL_BOUNDS_HEAVY void k_angles(...)
__global__ GFNFF_KERNEL_BOUNDS_HEAVY void k_dihedrals(...)
__global__ GFNFF_KERNEL_BOUNDS_HEAVY void k_inversions(...)
__global__ GFNFF_KERNEL_BOUNDS_HEAVY void k_storsions(...)
__global__ GFNFF_KERNEL_BOUNDS_HEAVY void k_hbonds(...)
__global__ GFNFF_KERNEL_BOUNDS_HEAVY void k_batm(...)
__global__ GFNFF_KERNEL_BOUNDS_HEAVY void k_atm(...)
```

**Schritt 5**: `getLaunchConfig()` in `ff_workspace_gpu.cu` anpassen. Aktuell wählt es Blockgröße 256 für n < 16384. Da SM_120 4 Blöcke/SM erlaubt, ist 256 mit `__launch_bounds__(256, 4)` optimal. Keine Änderung an der Logik nötig — nur das Makro steuert die SM-Belegung.

### Verifikation

```bash
# Build und ptxas-Output prüfen:
cd release && make -j4 2>&1 | grep -A2 "ptxas info.*k_dispersion\|k_bonds"
# Erwartung: k_dispersion mit (256,4): occupancy ~50% statt ~25%

# Regression:
ctest --output-on-failure
```

---

## Erwartetes Ergebnis

- `k_dispersion`, `k_repulsion`, `k_bonds`, `k_coulomb`: ~2× höhere SM-Occupancy auf SM_120
- L2-Hit-Rate für Koordinaten: +15–30% über alle 22 Kernel
- Keine Änderung an wissenschaftlicher Korrektheit (rein mechanische GPU-Optimierung)
