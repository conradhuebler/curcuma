# WP7: EEQ-Löser für große Systeme — GPU und CPU, auswählbar

**Kategorie**: Mittel–Groß  
**Aufwand**: ~15–25 Tage (alle drei Varianten)  
**Wirkung**: Ermöglicht GFN-FF für Systeme >10k Atome; schließt nfrag>1 GPU-Lücke  
**Abhängigkeiten**: WP5 (GPU-EEQ-Pipeline), WP2 ✅  
**Status**: 🤖 Geplant — Infrastruktur teilweise vorhanden

---

## Motivation

Der aktuelle EEQ-Löser verwendet Cholesky-Faktorisierung: O(N³/6). Für N=1410
kostet das ~10–15 ms/Schritt und dominiert die gesamte Rechenzeit. Für N=10k wäre
das ~5 Sekunden pro Schritt — nicht für MD geeignet.

Gleichzeitig fehlt für **nfrag>1** (Komplexe, Dimere, Polymere) der vollständig
GPU-seitige Schur-Komplement-Schritt: `solveWithDeviceRHS()` lädt z1+Z₂ auf die
CPU herunter (~N×nfrag×8 Byte) und macht den Schur-Schritt CPU-seitig.

Ziel: **drei EEQ-Löserstategien**, jeweils für GPU und CPU, über einen einzigen
Parameter `eeq_solver` auswählbar.

---

## Überblick: Drei Strategien

| Strategie | Physik | GPU-Komplexität | CPU-Komplexität | Ideal für |
|-----------|--------|----------------|----------------|-----------|
| **A: Vollständiges Schur (nfrag-allgemein)** | exakt | O(N³) Cholesky + O(nfrag²·N) Reduktionen | O(nfrag²·N) Reduktionen + O(nfrag³) | N<5k, nfrag 1–20 |
| **B: Gebatchtes per-Fragment Cholesky** | approximativ (no cross-frag) | O(Σ N_f³) gebatcht | O(Σ N_f³) parallelisiert | weit getrennte Fragmente, nfrag groß |
| **C: PCG-Löser mit Schur-Projektion** | exakt | O(k·N²) GEMV auf GPU | O(k·N²) parallelisiert | N>1k, MD/opt (warm start) |

---

## Strategie A: Allgemeines Schur-Komplement (nfrag-allgemein)

### Physik

Behält den vollen N×N-Cholesky — cross-fragment Coulomb-Terme sind korrekt enthalten.
Erweitert WP5-A (nfrag=1 Skalar) auf beliebiges nfrag.

### Mathematik

Nach `cusolverDnDpotrs` liegt auf Device:
```
d_rhs[ 0   .. N-1 ]  = z₁       (Atom-RHS-Lösung)
d_rhs[ N   .. 2N-1]  = Z₂[:,0]  (Fragment-0-Constraintlösung)
d_rhs[ 2N  .. 3N-1]  = Z₂[:,1]  (Fragment-1-Constraintlösung)
...
d_rhs[ k·N .. (k+1)N-1 ] = Z₂[:,k-1]
```

Benötigt werden:
```
Cz1[f]  = Σᵢ∈frag_f  z₁[i]             → nfrag Summen
S[f,g]  = Σᵢ∈frag_f  Z₂[i,g]           → nfrag² Summen (Schur-Matrix)
```

Dann CPU: löse `S·λ = Cz1 − Q_frag` (nfrag×nfrag-System, trivial für nfrag≤20).  
Dann GPU: `q[i] = z₁[i] − Σ_f Z₂[i,f]·λ[f]` (ein Kernel).

Download: **(2·nfrag + nfrag²) doubles** — für nfrag=2: 48 Byte statt ~5 KB.

### GPU-Implementierung

**Neuer Kernel** `k_eeq_reduce_fragment_sums`:
```cuda
// Für jedes Fragment f: reduziere z₁ und Z₂[:,f] über Atome i∈frag_f
// Input:  d_rhs (N×(nfrag+1)), d_frag_atom_map, d_frag_atom_offsets
// Output: d_Cz1[nfrag], d_S[nfrag*nfrag]
__global__ void k_eeq_reduce_fragment_sums(
    int N, int nfrag,
    const double* d_rhs,           // [N*(nfrag+1)], column-major
    const int*    d_frag_atom_map, // [N] globale Atomindizes sortiert nach Fragment
    const int*    d_frag_offsets,  // [nfrag+1] Startindex pro Fragment
    double* d_Cz1,                 // [nfrag] output
    double* d_S);                  // [nfrag*nfrag] output
```

**Neuer Kernel** `k_eeq_schur_general`:
```cuda
// q[i] = z₁[i] − Σ_f Z₂[i,f] · λ[f]
__global__ void k_eeq_schur_general(
    int N, int nfrag,
    const double* d_rhs,    // [N*(nfrag+1)] z₁ + Z₂
    const double* d_lambda, // [nfrag] Lagrange-Multiplikatoren (von CPU)
    double* d_charges);     // [N] output
```

**Neuer EEQSolverGPU-Aufruf**:
```cpp
bool EEQSolverGPU::solveWithDeviceRHSAndGPUSchurGeneral(
    int natoms, int nfrag,
    const double* cx, const double* cy, const double* cz,
    const double* d_alpha, const double* d_gam,
    const double* d_rhs_atoms,
    const std::vector<int>& fraglist,
    const std::vector<double>& rhs_constraints,
    double cutoff_sq,
    bool force_refactor);
// Gibt charges in getDeviceChargesPtr() zurück (d_rhs[0..N-1])
```

### CPU-Pendant (parallelisiert)

Der bestehende `EEQSolver::dispatchSolve()` + `solveWithSchurCholesky()` behandelt
beliebiges nfrag korrekt. Parallelisierbar über OpenMP für die Reduktionsschritte:

```cpp
// In EEQSolver::solveWithSchurCholesky():
// Step 2: Löse A·z₁ = b  (bereits Eigen, kann mit Eigen::BDCSVD::compute parallelisiert werden)
// Step 3: Löse A·Z₂ = Cᵀ  → nfrag unabhängige Rückwärtssubstitutionen
// → per OpenMP thread je Spalte parallelisieren (nfrag ≤ 20: trivial)
Matrix Z2(natoms, nfrag);
#pragma omp parallel for num_threads(std::min(nfrag, omp_get_max_threads()))
for (int f = 0; f < nfrag; ++f) {
    Z2.col(f) = llt.solve(C.row(f).transpose());
}
```

**Aufwand GPU**: ~1–2 Tage (2 neue Kernels + `solveWithDeviceRHSAndGPUSchurGeneral`)  
**Aufwand CPU**: ~4h (OpenMP in `solveWithSchurCholesky`)

---

## Strategie B: Gebatchtes per-Fragment Cholesky

### Physik

**Approximation**: Die N×N-Matrix wird in nfrag unabhängige N_f×N_f-Blöcke aufgeteilt.
Cross-Fragment Coulomb-Terme werden ignoriert.

**Gültigkeitsbereich**: Fragmente mit Abstand ≥ ~8 Å (≥ 15 Bohr):
```
J_cross = erf(γ·r)/r  für r=15 Bohr ≈ 0.067 Eh/e²
A_diag  ≈ η + √(2/π)/√α ≈ 0.5 Eh/e²
Relativer Fehler: ~13% der Diagonale → Ladungsfehler typisch <5%
```

Für eng gepackte Komplexe (π-Systeme, Kristalle) **nicht verwenden**.

### Wann sinnvoll?

- Polymer-Modelle mit unabhängigen Monomeren (Abstand > 8 Å)
- Gas-Phase-Komplexe mit schwacher Wechselwirkung
- Voruntersuchungen mit N > 5k wo exakter Solver zu langsam

### GPU-Implementierung

`solveWithDeviceRHSAndGPUSchurBatched()` ist bereits vollständig implementiert
(`eeq_solver_gpu.cu:930`). Kernels `k_eeq_build_fragment_matrices` und
`k_eeq_gather_rhs_fragments` existieren in `gfnff_kernels.cu`.

**Status**: Deaktiviert in `gfnff_gpu_method.cpp:580` — muss nur aktiviert und
hinter dem `eeq_solver`-Parameter geschaltet werden:

```cpp
// gfnff_gpu_method.cpp — Dispatch nach Strategie:
if (m_eeq_solver_strategy == "batched" && m_eeq_gpu->isFragmentTopoValid()) {
    eeq_ok = m_eeq_gpu->solveWithDeviceRHSAndGPUSchurBatched(...);
} else if (m_eeq_nfrag == 1) {
    eeq_ok = m_eeq_gpu->solveWithDeviceRHSAndGPUSchur(...);  // WP5-A
} else {
    eeq_ok = m_eeq_gpu->solveWithDeviceRHSAndGPUSchurGeneral(...);  // WP7-A
}
```

### CPU-Pendant (parallelisiert)

N unabhängige Cholesky-Faktorisierungen — ideal für OpenMP:

```cpp
// In EEQSolver (neue Methode calcBatchedFragmentCharges):
#pragma omp parallel for schedule(dynamic)
for (int f = 0; f < nfrag; ++f) {
    // Baue N_f × N_f Matrix für Fragment f
    Matrix A_f = buildFragmentMatrix(f, ...);
    // Löse EEQ für Fragment f (Schur exakt auf kleiner Matrix)
    charges_f = solveFragmentEEQ(A_f, rhs_f, Q_f);
}
```

Speedup: nfrag × O(N_f³) statt O(N³) — für N=10k, nfrag=100, N_f=100:
- Batched: 100 × (100³/6) ≈ 17M FLOP
- Full: (10000³/6) ≈ 167G FLOP  
→ ~10000× schneller (Cholesky-Schritt allein)

**Aufwand GPU**: ~2 Tage (Dispatch + Test)  
**Aufwand CPU**: ~2 Tage (neue Methode + OpenMP)

---

## Strategie C: PCG-Löser mit Fragment-Projektion

### Physik

Exakt (kein Approximationsfehler). Löst `A·q = b` iterativ in O(k·N²) wobei k die
Iterationsanzahl ist. Mit Warm-Start (q aus Vorschritt) typisch k=5–20 für MD.

Fragment-Nebenbedingungen werden als Projektion nach jedem CG-Schritt erzwungen:
```
nach jeder CG-Iteration: für jedes Fragment f:
  delta = (Σᵢ∈f qᵢ − Q_f) / N_f
  qᵢ -= delta  ∀ i ∈ f
```
Diese Projektion erhält die CG-Konvergenz und erzwingt Ladungserhaltung exakt.

### GPU-Implementierung

Jede CG-Iteration besteht aus:
1. `cublasDgemv`: `A·d` — Matrix-Vektor-Produkt, O(N²), GPU-nativ
2. `cublasDdot`: Skalarprodukte für α, β — trivial
3. `k_pcg_update`: Vektorupdates q, r, d — elementweise
4. `k_pcg_project_fragments`: Fragment-Projektion — O(N) Reduktion + Korrektur

```cuda
// k_pcg_project_fragments: erzwingt Σᵢ∈f qᵢ = Q_f für alle f
__global__ void k_pcg_project_fragments(
    int N, int nfrag,
    double* d_q,                   // [N] Ladungsvektor (in-place)
    const int* d_frag_atom_map,    // [N] sortierte Atomindizes
    const int* d_frag_offsets,     // [nfrag+1]
    const double* d_Q_frag,        // [nfrag] Zielladungen
    const int* d_frag_sizes);      // [nfrag]
```

**Warm-Start**: `d_q` vom Vorschritt wird beibehalten (D2D — kein Transfer).

**Jacobi-Vorkonditionierer**: `M = diag(A)` — ist bereits als `d_gam + sqrt(2/π)/sqrt(d_alpha)`,
topology-konstant nach WP2. Ein Kernel lädt aus `eeq_topo.d_gam` und `eeq_topo.d_alpha`.

**Auto-Benchmark wie CPU-PCG**: Beim ersten Aufruf Benchmark zwischen Cholesky und PCG;
für große N (>500) PCG bevorzugen.

### CPU-Pendant (parallelisiert)

Der CPU-PCG-Löser `EEQSolver::solveWithPCG()` ist bereits implementiert
(`eeq_solver.cpp`, `EEQSolveMethod::PCG`). Parallelisierung der Matrix-Vektor-
Multiplikation `A·d`:

```cpp
// In solveWithPCG() — Matrix-Vektor-Produkt parallelisieren:
// Aktuell: serielle Eigen-Multiplikation
Vector Ad = A * d;  // → Eigen::EIGEN_VECTORIZE mit OpenMP Threads

// Besser: explizites OpenMP
#pragma omp parallel for schedule(static)
for (int i = 0; i < natoms; ++i) {
    Ad[i] = A.row(i).dot(d);  // Vektorzugriff optimiert
}
```

Alternativ: Eigen nutzt intern BLAS-Parallelisierung wenn Eigen::setNbThreads() gesetzt.

**Fragment-Projektion** in CPU-PCG (nach jedem CG-Schritt):
```cpp
for (int f = 0; f < nfrag; ++f) {
    double sum = 0.0;
    for (int i : frag_atoms[f]) sum += q[i];
    double delta = (sum - Q_frag[f]) / frag_sizes[f];
    for (int i : frag_atoms[f]) q[i] -= delta;
}
```

**Aufwand GPU**: ~5–7 Tage (cublasDgemv-Integration + Projektions-Kernel + Warm-Start)  
**Aufwand CPU**: ~2–3 Tage (Fragment-Projektion + OpenMP MATVEC)

---

## Parameter-Interface

Alle Strategien über einen einzigen Parameter `eeq_solver` steuerbar.

### Parameterdefinition (eeq_solver.h)

```cpp
PARAM(eeq_solver_strategy, String, "auto",
      "EEQ linear solver strategy:\n"
      "  auto          — Cholesky für N<500, PCG für N≥500 (Benchmark beim ersten Schritt)\n"
      "  cholesky      — Voller N×N Cholesky + allgemeines Schur-Komplement (exakt, O(N³))\n"
      "  batched       — Per-Fragment Cholesky (approximativ, ignoriert cross-fragment Coulomb)\n"
      "  pcg           — Preconditioned Conjugate Gradient mit Fragment-Projektion (exakt, O(k·N²))\n"
      "  cholesky_gpu  — Wie 'cholesky', aber erzwingt GPU-Pfad auch wenn nfrag>1\n"
      "  batched_gpu   — Wie 'batched', GPU-gebatcht (erfordert weit getrennte Fragmente)\n"
      "  pcg_gpu       — Wie 'pcg', cublasDgemv auf GPU\n",
      "Algorithm", {})

PARAM(eeq_pcg_max_iter, Int, 200,
      "PCG: maximale Iterationen (wird für große N automatisch skaliert)", "Algorithm", {})

PARAM(eeq_pcg_tol, Double, 1e-10,
      "PCG: Konvergenztoleranz auf |r|", "Algorithm", {})

PARAM(eeq_batched_warn_threshold, Double, 15.0,
      "Batched-Strategie: Warnung wenn Fragment-Abstand < diesem Wert (Bohr). "
      "0 = keine Warnung.", "Advanced", {})
```

### CLI-Beispiele

```bash
# Standard (auto): Cholesky für kleine, PCG für große Systeme
./curcuma -md protein.xyz -method gfnff -gpu cuda

# Erzwinge PCG für großes Polymer (N=10k)
./curcuma -md polymer_10k.xyz -method gfnff -gpu cuda -eeq_solver_strategy pcg_gpu

# Batched für weit getrennte Fragmentsysteme (schnell, approximativ)
./curcuma -sp cluster_100frag.xyz -method gfnff -eeq_solver_strategy batched

# Debug: Vergleich exakt vs. approximativ
./curcuma -sp complex.xyz -method gfnff -eeq_solver_strategy cholesky  # Referenz
./curcuma -sp complex.xyz -method gfnff -eeq_solver_strategy batched   # Approximation
```

---

## Dispatch-Logik (gfnff_gpu_method.cpp)

```cpp
// Strategie-Auswahl im EEQ-Solve-Block:
std::string strategy = m_eeq_solver_strategy;

// Auto-Auswahl falls nicht explizit gesetzt
if (strategy == "auto") {
    if (N >= m_eeq_pcg_threshold)
        strategy = use_device_rhs ? "pcg_gpu" : "pcg";
    else
        strategy = use_device_rhs ? "cholesky_gpu" : "cholesky";
}

// Dispatch
if (strategy == "cholesky_gpu" || strategy == "cholesky") {
    if (m_eeq_nfrag == 1) {
        eeq_ok = m_eeq_gpu->solveWithDeviceRHSAndGPUSchur(...);      // WP5-A
    } else {
        eeq_ok = m_eeq_gpu->solveWithDeviceRHSAndGPUSchurGeneral(...);//WP7-A
    }
} else if (strategy == "batched_gpu" || strategy == "batched") {
    if (m_eeq_gpu->isFragmentTopoValid()) {
        warnIfFragmentsTooClose();  // Warnung wenn Abstand < threshold
        eeq_ok = m_eeq_gpu->solveWithDeviceRHSAndGPUSchurBatched(...);//WP6
    } else {
        // Fallback auf cholesky wenn Topo nicht initialisiert
        eeq_ok = m_eeq_gpu->solveWithDeviceRHSAndGPUSchurGeneral(...);
    }
} else if (strategy == "pcg_gpu" || strategy == "pcg") {
    eeq_ok = m_eeq_gpu->solveWithDeviceRHSAndGPUPCG(...);            // WP7-C
}

// Letzter Fallback: CPU EEQ
if (!eeq_ok) {
    // bestehender CPU-Pfad bleibt immer verfügbar
}
```

---

## Skalierungsanalyse für N=10k

| Strategie | Cholesky-Zeit (geschätzt) | Annahmen |
|-----------|--------------------------|---------|
| A: Voll Cholesky | ~500 ms/Schritt | O(10000³/6), RTX 5080 ~8 TFLOP |
| B: Batched (nfrag=100, N_f=100) | ~0.05 ms/Schritt | 100×O(100³/6), parallelisiert |
| C: PCG (k=20 Iterationen) | ~20 ms/Schritt | 20×O(10000²), cublasDgemv |

Für N=10k ist Strategie A nicht praktikabel. Strategie B setzt physikalische
Separierung voraus. Strategie C (PCG) ist der allgemeine Weg für große Systeme.

---

## Implementierungsreihenfolge

```
WP7-A: solveWithDeviceRHSAndGPUSchurGeneral (nfrag>1, korrekte Physik)
  → Aufwand: ~2 Tage
  → Schließt die nfrag>1 GPU-Lücke für complex.xyz

WP7-B: Batched aktivieren + Fragment-Abstandswarnung + CPU-OpenMP
  → Aufwand: ~3 Tage
  → Bereits implementiert, nur Aktivierung + Validierung

WP7-C: PCG auf GPU (cublasDgemv + Fragment-Projektion)
  → Aufwand: ~7 Tage
  → Ermöglicht N>5k MD

Parameter-Interface (eeq_solver_strategy)
  → Aufwand: ~1 Tag (PARAM-Makro + Dispatch-Logik)
  → Parallel zu WP7-A implementierbar
```

---

## Validierung

Für jede neue Strategie muss gelten:
1. Energie = CPU-Referenz-EEQ ± 1 µEh (Strategie A, C)
2. Energie ≈ CPU-Referenz ± 1 mEh (Strategie B, nur für gut separierte Fragmente)
3. Gradient-Konsistenz: numerischer Gradient ≈ analytischer Gradient (±1e-5)
4. MD-Energieerhaltung: Drift < 1 mEh/ns (alle Strategien)

```bash
# Validierungsbefehle:
./curcuma -sp complex.xyz -method gfnff -eeq_solver_strategy cholesky   # Referenz
./curcuma -sp complex.xyz -method gfnff -eeq_solver_strategy cholesky_gpu
./curcuma -sp complex.xyz -method gfnff -eeq_solver_strategy batched -eeq_batched_warn_threshold 0

# MD-Energieerhaltung:
./curcuma -md polymer.xyz -method gfnff -gpu cuda -eeq_solver_strategy pcg_gpu -steps 1000
# → Energiedrift in Ausgabe prüfen
```

---

## Dateien

| Datei | Änderungen |
|-------|-----------|
| `eeq_solver_gpu.cu` | `solveWithDeviceRHSAndGPUSchurGeneral()`, `solveWithDeviceRHSAndGPUPCG()` |
| `eeq_solver_gpu.h` | Neue Methoden-Deklarationen |
| `gfnff_kernels.cu` | `k_eeq_reduce_fragment_sums`, `k_eeq_schur_general`, `k_pcg_*` |
| `gfnff_kernels.cuh` | Kernel-Deklarationen |
| `gfnff_gpu_method.cpp` | Dispatch-Logik, `m_eeq_solver_strategy` Member |
| `eeq_solver.h` | PARAM `eeq_solver_strategy` + `eeq_pcg_max_iter` etc. |
| `eeq_solver.cpp` | Fragment-Projektion in PCG, OpenMP MATVEC |

---

## WP7-D: Block-Jacobi-Präkonditionierer + kontaktbewusste Auswahl (Jun 2026) ✅

**Status**: ⚙️ AI-implementiert, maschinell getestet (CPU + CUDA)

Zwei Schwächen im Viel-Fragment-Regime (Solvensboxen, Wirt-Gast, vdW-Komplexe)
geschlossen:

### 1. GPU-PCG: Block-Jacobi statt diagonalem Jacobi

Der GPU-PCG (Strategie C) nutzte nur einen **diagonalen** Jacobi-Präkonditionierer
(`d_pcg_M_inv = 1/A[i,i]`) → langsame Konvergenz bei vielen Fragmenten. Portiert
wurde der **per-Fragment Block-Jacobi** der CPU (`EEQSolver::buildBlockJacobi`):
`M⁻¹ = blockdiag(A_ff⁻¹)`. Senkt die PCG-Iterationen von ~30–100 auf ~2–5, **exakt**.

- `buildBlockJacobiFactors()` (`eeq_solver_gpu.cu`): baut die per-Fragment-Blöcke
  (`k_eeq_build_fragment_matrices` → `d_A_blocks`), faktorisiert (`cusolverDnDpotrf`)
  und invertiert (`cusolverDnDpotri`) je Block, symmetrisiert (`k_eeq_symmetrize_blocks`).
  Nicht-SPD oder zu großer Block (N_f > 2048) ⇒ Rückfall auf diagonalen Jacobi.
- Apply (`k_eeq_block_jacobi_apply`, 1 Block/Fragment): `z = blockdiag(A_ff⁻¹)·r` als
  symmetrisches GEMV, Gather/Scatter über `frag_atom_map`; r/z bleiben in globaler
  Atomreihenfolge (wie CPU `BlockJacobiPC::apply`). Nur beim Refaktorisieren gebaut,
  über MD-Schritte amortisiert. Wiederverwendet die WP6-Fragmenttopologie.
- Verifiziert: GPU-PCG+Block-Jacobi == GPU-SchurCholesky bit-identisch (≤1e-8,
  27-Wasser-Cluster, `-opt`, Schritt-für-Schritt). `m_pcg_block_jacobi_valid` wird bei
  Topologiewechsel und im (d_A_blocks überschreibenden) Batched-Pfad invalidiert.

### 2. Kontaktbewusste Auswahl (CPU): exakt statt Batched bei Kontakt

Der CPU-Dispatcher wählte Batched (Strategie B, **lässt cross-fragment Coulomb still
weg**) allein anhand der Fragmentdichte (`nfrag/N`), unabhängig davon, ob Fragmente
in Kontakt sind. Neu: `m_contact_min_dist` (minimaler Inter-Fragment-Atomabstand,
aus den bereits gebauten gepackten Distanzen + `fraglist`, O(N²) auf der Distanzphase).
Bei Kontakt (`< eeq_batched_min_distance`, Default 15 Bohr) oder unbekannt wählt
Auto/PCG den **exakten** Löser (SchurCholesky bzw. PCG+Block-Jacobi). Neuer Schalter
`eeq_contact_prefer_exact` (Default `true`); `false` stellt die alte dichte-basierte
Batched-Auswahl wieder her. Batched bleibt als explizite Option (`-solve_method
batched`) für gut getrennte Fragmente erhalten.

- Verifiziert (27-Wasser-Cluster, in Kontakt, 3.19 Bohr): Default-Auto wählt jetzt
  „exact solver", E = −8.77203498 Eh (== SchurCholesky == PCG+Block-Jacobi); Batched
  weicht um 3.58 mEh ab (cross-fragment Coulomb fehlt).

**GPU-Dispatch:** Auf der GPU wählt `Auto` ohnehin nie Batched (nur PCG bzw.
SchurCholesky, beide exakt); Batched läuft nur bei explizitem `-gfnff.solve_method
batched` (mit bestehender Kontakt-Warnung). Daher genügt GPU-seitig der bessere
PCG-Präkonditionierer; keine Dispatch-Änderung nötig.

### Offen / nicht umgesetzt

- **FMM/Treecode** für das cross-fragment Coulomb-Matvec (O(N log N) statt O(N²)) —
  großer Aufwand, Nische; das PCG-Matvec bleibt dicht O(N²).
- **ROCm-Spiegelung** des Block-Jacobi-Applys (`rocm/eeq_solver_hip`) — naheliegender
  Hipify-Folgeschritt, hier (CUDA-Scope) nicht umgesetzt.

### Geänderte Dateien (WP7-D)

| Datei | Änderungen |
|-------|-----------|
| `gfnff_kernels.cu/.cuh` | `k_eeq_block_jacobi_apply`, `k_eeq_symmetrize_blocks` |
| `eeq_solver_gpu.cu` | `buildBlockJacobiFactors()`, `applyPrecondPCG()`, `m_pcg_block_jacobi_valid`, PCG-Verdrahtung |
| `eeq_solver.h` | PARAM `eeq_contact_prefer_exact`, Member `m_contact_min_dist` |
| `eeq_solver.cpp` | Kontaktmetrik in `calculateFinalCharges`, kontaktbewusste Batched-Auswahl in `dispatchSolve` |
