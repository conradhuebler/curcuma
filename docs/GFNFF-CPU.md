# GFN-FF CPU Optimierungsanalyse

## Kontext

Die native GFN-FF CPU-Implementierung ist funktional korrekt, hat aber strukturelle Engpässe die bei größeren Molekülen (>200 Atome) und MD-Simulationen spürbar werden. Diese Analyse identifiziert konkrete, messbare Optimierungen — geordnet nach Aufwand/Nutzen.

**Analysierte Dateien:**
- `src/core/energy_calculators/ff_methods/forcefieldthread.cpp` — innere Schleifen
- `src/core/energy_calculators/ff_methods/forcefieldthread.h` — Structs, Datenhaltung
- `src/core/energy_calculators/ff_methods/forcefield.cpp` — Thread-Split, AutoRanges
- `src/core/energy_calculators/ff_methods/gfnff_parameters.h` — Pair-Structs
- `src/core/energy_calculators/ff_methods/gfnff_method.cpp` — Topologie, EEQ, D4

---

## Befund: Konkrete Bottlenecks

### 1. `std::set<pair<int,int>>` für Bonded-Pair-Ausschluss
**Datei:** `forcefieldthread.h:609`, `forcefieldthread.cpp:106-115`
**Problem:** Für jeden nicht-gebundenen Paar-Term (Dispersion, Coulomb, Repulsion) wird geprüft ob das Paar gebunden ist — via `m_bonded_pairs.count({i,j})`, das ist O(log N_bonds) pro Lookup mit Heap-Allokation pro Knoten im Baum.
**Lösung:** Ersetzen durch flache obere Dreiecksbitmatrix `std::vector<bool>` (N×(N+1)/2 Bits), Index-Formel `min(i,j)*N + max(i,j)`. Lookup O(1), cache-freundlich.
- Für N=100: 1250 Bytes statt ~6 KB (std::set + tree-Overhead)

**Aufwand:** 2-3h. **Nutzen:** 5-10% für Coulomb/Repulsion-Terme.

---

### 2. N×N Distanzmatrix — O(N²) Speicher und Redundanz
**Datei:** `gfnff_method.cpp:651-663`
**Problem:** Jede Geometrieänderung berechnet `topo.distance_matrix` vollständig (N×N doubles). Größe: 80 KB bei N=100, 8 MB bei N=1000, **800 MB bei N=3000**. Die Matrix wird für CN-Berechnung (Zeile 686) und EEQ verwendet, aber CN braucht nur Distanzen zu kurzen Nachbarn (r < 6 Bohr).
**Lösung:**
1. CN-Berechnung: Nur Paare mit r < CN_cutoff (≈6 Bohr) — O(N×avg_neighbors) statt O(N²)
2. EEQ Phase 2: Distanzen on-the-fly in der Coulomb-Paar-Schleife (bereits vorhanden in `generateGFNFFCoulombPairs`)
3. Distanzmatrix aus `TopologyInfo` entfernen oder auf Optional<> reduzieren

**Aufwand:** 1 Tag. **Nutzen:** Speicher −90% bei N>500; kritisch für MD mit >500-Atom-Systemen.

---

### 3. `GFNFFDispersion` Struct: Legacy-Felder verschwenden Cache
**Datei:** `gfnff_parameters.h:83-97`, heißer Loop `forcefieldthread.cpp:1893`
**Problem:** Struct ist 88 Bytes, enthält aber 5 Legacy-Felder (`C8`, `s6`, `s8`, `a1`, `a2`) die in GFN-FF nicht gelesen werden (nur D3-Kompatibilität). Der innerste Dispersion-Loop lädt 88 Bytes/Paar, braucht aber effektiv nur ~48 Bytes (`i`, `j`, `C6`, `r4r2ij`, `r0_squared`, `r_cut`, `zetac6`). Bei 500K Paaren: 44 MB Cache-Durchsatz statt 24 MB.
**Lösung:** Legacy-Felder in separaten `D3DispersionExtra`-Vektor auslagern (nur für UFF-D3 verwendet), `GFNFFDispersion` auf 48 Bytes schrumpfen.

**Aufwand:** 3-4h. **Nutzen:** 20-30% Speedup im Dispersion-Term (dominanter Term für größere Systeme).

---

### 4. Thread-Lastbalancierung: Linearer Strip führt zu Imbalance
**Datei:** `forcefield.cpp:1545-1640` (AutoRanges)
**Problem:** Alle Terme verwenden linearen Strip `[i*N/T, (i+1)*N/T)`. Aber Dispersion hat ~3000 Paare, Coulomb ~5000, HB nur ~50. Wenn Thread 0 alle HB-Paare bekommt (schnell fertig), wartet er auf Thread 3 mit Dispersion.
**Lösung:** Round-Robin-Verteilung für dominante Terme:
```cpp
for (int j = i; j < pairs.size(); j += thread_count)
    thread->addGFNFFDispersion(pairs[j]);
```
Das verteilt gleichmäßig wenn Pair-Rechenzeiten variieren.

**Aufwand:** 1-2h. **Nutzen:** 10-20% Gesamt-Speedup bei 4+ Threads.

---

### 5. CN-Berechnung: O(N²) `erf()`-Aufrufe jedes MD-Schritt
**Datei:** `gfnff_method.cpp:669-693`
**Problem:** CN-Loop iteriert N×N Paare mit `std::erf()` pro Paar. Bei N=200: 40.000 teure Transcendental-Aufrufe pro MD-Schritt. `erf()` ist ~20× langsamer als einfache Arithmetik.
**Lösung (empfohlen):** Neighbor-List mit Cutoff r_CN ≈ 6 Bohr. Nur Nachbarn innerhalb des Cutoffs tragen zu CN bei — praktisch alle anderen Terme sind exp(-x) → 0. Reduziert von O(N²) auf O(N × 10).

**Aufwand:** 1 Tag. **Nutzen:** 5-10× Speedup im CN-Term bei N>200; bei MD dominiert dieser Term ohne Topologie-Caching.

---

### 6. D4 Gaussian-Gewichte: Unnötige Neuberechnung bei MD
**Datei:** D4ParameterGenerator::precomputeGaussianWeights, aufgerufen aus gfnff_method.cpp
**Problem:** Bei jedem Energie-Aufruf werden `gw = exp(-4*(CN-CN_ref)²)` für N×7 Referenzzustände neu berechnet. In MD-Simulationen ändern sich CNs pro Schritt um <0.01 — Gewichte bleiben praktisch konstant.
**Lösung:** Threshold-Caching: `if (max_cn_change < 0.01) return;` — nutzt vorhandene CN-Tracking-Infrastruktur aus `needsFullTopologyUpdate()`.

**Aufwand:** 2-3h. **Nutzen:** ~50% Zeitersparnis im D4-Term für MD.

---

### 7. `GFNFFCoulomb`: Dynamische Ladungen unnötig im Struct gespeichert
**Datei:** `gfnff_parameters.h:118-131`
**Problem:** 144-Byte-Struct speichert `q_i`, `q_j` pro Paar (16 Bytes), obwohl `m_eeq_charges_ptr` bereits als Zeiger auf den Ladungsvektor existiert. Im Loop könnte man `charges[coulomb.i]` direkt lesen statt das Struct zu updaten.
**Lösung:** `q_i`, `q_j` aus Struct entfernen; im Coulomb-Loop `m_eeq_charges_ptr->operator()(coulomb.i)` verwenden (Zeiger existiert bereits, `forcefieldthread.h:217`). Struct schrumpft von 144 auf 128 Bytes.

**Aufwand:** 4-6h. **Nutzen:** ~10% Coulomb-Term, kein EEQ-Update-Overhead mehr.

---

## Priorisierte Roadmap

| # | Optimierung | Aufwand | Nutzen | Zielgruppe |
|---|-------------|---------|--------|------------|
| A1 | D4 Gaussian-Gewichte Threshold-Cache | 2-3h | 50% D4-Term (MD) | MD >100 Atome |
| A2 | Round-Robin Thread-Split (AutoRanges) | 1-2h | 10-20% Gesamt | Multi-Core |
| B1 | N×N Distanzmatrix eliminieren | 1 Tag | −90% Speicher bei N>500 | MD, Optimierung |
| B2 | GFNFFDispersion Struct verkleinern | 3-4h | 20-30% Dispersion-Term | Alle Systeme |
| C1 | CN Neighbor-List-Cutoff | 1 Tag | 5-10× CN-Term N>200 | Große Moleküle |
| C2 | GFNFFCoulomb q_i/q_j entfernen | 4-6h | 10% Coulomb-Term | Alle Systeme |
| C3 | std::set → Bitset Bonded-Ausschluss | 2-3h | 5-10% NB-Terme | Kleine Moleküle |

**Nicht empfohlen (geringer Effekt):**
- OpenMP/SIMD-Pragmas: CxxThreadPool + Eigen-SIMD reichen. Kein messbarer Mehrwert.
- Cutoff-Reduktion: Physikalisch riskant ohne Referenz-Validierung.
- JSON-Roundtrip-Eliminierung: Nativer `GFNFFParameterSet`-Pfad existiert bereits (`forcefield.cpp:449`).

## Verifikation nach jeder Änderung

```bash
cd release && make -j4
ctest -R "gfnff" --output-on-failure
./test_cases/test_gfnff_numgrad               # Gradient-Regression
```
MD-Performance: SimpleMD-Laufzeit mit 200-Atom-System vor/nach messen.

---

# GFN-FF Periodic Boundary Conditions (CPU + GPU)

## Context

The native GFN-FF implementation (`gfnff`) currently computes all pair distances as simple Euclidean distances. For periodic systems (crystals, periodic MD boxes) the minimum-image convention (MIC) must be applied so each atom pair uses the nearest periodic image. The Fortran GFN-FF reference does NOT support PBC — this is a Curcuma extension. SimpleMD already wraps coordinates post-step, but the energy/gradient calculation ignores periodicity entirely.

**What already exists (no reimplementation needed):**
- `src/tools/pbc_utils.h` — complete MIC library: `applyMinimumImage(r, cell, cell_inv)`, `calculateDistancePBC()`, triclinic support
- `Molecule::hasPBC()`, `getUnitCell()`, `getUnitCellInverse()` — already set by SimpleMD and XYZ readers
- `SimpleMD::wrapAllAtomsIntoCentralCell()` — wraps after each step (already working)

**What is missing:**
- ForceField/ForceFieldThread have no lattice data — all distances are Euclidean
- GPU kernels have no lattice upload or MIC helper
- GFNFF class does not extract the lattice from the Molecule and pass it down

## Architecture Overview

```
Molecule (has_pbc, unit_cell)
  ↓ GFNFF::InitialiseMolecule() — extract lattice
GFNFF (m_unit_cell, m_has_pbc)
  ↓ after setParameter()
ForceField::setUnitCell(cell, has_pbc)
  ↓ AutoRanges() distributes to threads
ForceFieldThread (m_unit_cell, m_unit_cell_inv)
  ↓ in all pair distance loops
PBCUtils::applyMinimumImage(r_vec, cell, cell_inv)

GPU path:
GFNFFGPUMethod → upload d_lattice[9] + d_lattice_inv[9]
All 22 kernels → __device__ applyMIC(dx,dy,dz, lattice, inv)
```

## Implementation Plan

### Step 1: ForceField / ForceFieldThread (CPU)

**File: `src/core/energy_calculators/ff_methods/forcefield.h`**
- Add members: `Eigen::Matrix3d m_unit_cell`, `Eigen::Matrix3d m_unit_cell_inv`, `bool m_has_pbc = false`
- Add method: `void setUnitCell(const Eigen::Matrix3d& cell, bool has_pbc)`

**File: `src/core/energy_calculators/ff_methods/forcefield.cpp`**
- Implement `setUnitCell()`: store cell + compute inverse (`cell.inverse()`), set flag
- In `AutoRanges()`: pass cell, inv, flag to each ForceFieldThread via `thread->setUnitCell()`

**File: `src/core/energy_calculators/ff_methods/forcefieldthread.h`**
- Add members: `Eigen::Matrix3d m_unit_cell`, `Eigen::Matrix3d m_unit_cell_inv`, `bool m_has_pbc = false`
- Add method: `void setUnitCell(const Eigen::Matrix3d& cell, const Eigen::Matrix3d& inv, bool has_pbc)`
- Add `#include "src/tools/pbc_utils.h"`

**File: `src/core/energy_calculators/ff_methods/forcefieldthread.cpp`**

Apply MIC to pair vector `r = pos_j - pos_i` before taking `.norm()` in these locations:

| Term | Approx. location | Change |
|------|-----------------|--------|
| Dispersion | ~line 1899 | `if(m_has_pbc) r_vec = PBCUtils::applyMinimumImage(r_vec, m_unit_cell, m_unit_cell_inv);` |
| Coulomb | ~line 2189 | same pattern |
| Bonded repulsion | ~line 2081 | same |
| Non-bonded repulsion | ~line 2189 | same |
| Hydrogen bonds (r_AH, r_HB, r_AB) | ~lines 2374-2428 | apply MIC to each pair vector |
| Bonds | ~line 637 | same (bonded, but needed for images across cell face) |
| Angles | distance pairs | same |
| Torsions | distance pairs | same |

**Pattern for each location:**
```cpp
Eigen::Vector3d r_vec = pos_j - pos_i;
if (m_has_pbc)
    r_vec = PBCUtils::applyMinimumImage(r_vec, m_unit_cell, m_unit_cell_inv);
double r = r_vec.norm() * m_au;
```

### Step 2: GFNFF class — extract and propagate lattice

**File: `src/core/energy_calculators/ff_methods/gfnff_method.cpp`**

In `GFNFF::InitialiseMolecule(const Mol& mol)` (or wherever Molecule is first processed):
```cpp
m_has_pbc = mol.hasPBC();
if (m_has_pbc)
    m_unit_cell = mol.getUnitCell();  // Angstrom 3×3
```

In `GFNFF::initializeForceField()` (after `m_forcefield->setParameter()`):
```cpp
if (m_has_pbc)
    m_forcefield->setUnitCell(m_unit_cell, true);
```

In `GFNFF::UpdateMolecule(const Matrix& geometry)`:
```cpp
if (m_has_pbc)
    m_forcefield->setUnitCell(m_unit_cell, true);  // cell doesn't change, but threads may be recreated
```

**File: `src/core/energy_calculators/ff_methods/gfnff.h`**
- Add members: `Eigen::Matrix3d m_unit_cell`, `bool m_has_pbc = false`

**Note on Fortran reference**: The external GFN-FF (`external/gfnff/src/`) has NO active PBC implementation. Only dead parameter arrays (`pbcq`, `pbch`, `pbcs` in `dftd4param.f90`) are defined but never used. All distances are pure Euclidean. PBC in Curcuma is therefore a new extension beyond the Fortran reference — no reference values exist for validation; translational invariance and energy conservation in MD serve as the primary validation criteria.

**Note on units**: PBCUtils works in Angstrom; ForceFieldThread uses Bohr internally (multiplies by `m_au` = Bohr_to_Angstrom⁻¹). The unit cell from `Molecule::getUnitCell()` is in Angstrom. Pair vectors from `geom()` are in Bohr. Therefore: store cell in Bohr in ForceFieldThread (convert once in `setUnitCell()`).

### Step 3: GPU — lattice upload + MIC device helper

**File: `src/core/energy_calculators/ff_methods/cuda/ff_workspace_gpu.h`**
- Add method: `void setUnitCell(const Eigen::Matrix3d& cell_bohr, bool has_pbc)`
- Add members: `bool m_has_pbc = false`
- Add pinned/device buffers: `double m_h_lattice[9]`, `double m_h_lattice_inv[9]`

**File: `src/core/energy_calculators/ff_methods/cuda/ff_workspace_gpu.cu`**
- Add device constants or kernel parameters: `__constant__ double d_lattice[9]`, `__constant__ double d_lattice_inv[9]`, `__constant__ bool d_has_pbc`
- Implement `setUnitCell()`: compute inverse, copy to constant memory via `cudaMemcpyToSymbol`
- Add `__device__` helper (in kernel file, before all kernels):
```cuda
__device__ inline void applyMIC(double& dx, double& dy, double& dz) {
    if (!d_has_pbc) return;
    // frac = L_inv * r
    double fx = d_lattice_inv[0]*dx + d_lattice_inv[3]*dy + d_lattice_inv[6]*dz;
    double fy = d_lattice_inv[1]*dx + d_lattice_inv[4]*dy + d_lattice_inv[7]*dz;
    double fz = d_lattice_inv[2]*dx + d_lattice_inv[5]*dy + d_lattice_inv[8]*dz;
    // wrap to [-0.5, 0.5)
    fx -= round(fx); fy -= round(fy); fz -= round(fz);
    // back to Cartesian
    dx = d_lattice[0]*fx + d_lattice[3]*fy + d_lattice[6]*fz;
    dy = d_lattice[1]*fx + d_lattice[4]*fy + d_lattice[7]*fz;
    dz = d_lattice[2]*fx + d_lattice[5]*fy + d_lattice[8]*fz;
}
```
- In all 22 kernels, after computing `dx = cx[j] - cx[i]` etc., add: `applyMIC(dx, dy, dz);`
- Gradient sign: MIC affects position, not gradient direction — the gradient uses the same (MIC-corrected) displacement vector, so gradient is automatically correct.

**File: `src/core/energy_calculators/qm_methods/gfnff_gpu_method.cpp`**
- In `GFNFFGPUMethod::calculateEnergy()`, before launching kernels:
```cpp
if (m_gfnff->hasPBC())
    m_gpu_workspace->setUnitCell(m_gfnff->getUnitCellBohr(), true);
```

### Step 4: Topology generation for periodic systems

For bonded terms (bonds, angles, torsions, inversions): topology is generated from connectivity; no change needed — MIC in the energy calculation handles image bonds correctly.

For non-bonded pair lists (dispersion, coulomb, repulsion): currently all N×(N-1)/2 pairs. For periodic systems with a cutoff smaller than L/2 this is fine. For very small boxes (cutoff > L/2), explicit image enumeration would be needed — **out of scope for this implementation**; document as known limitation.

### Step 5: Input — reading the unit cell

XYZ files with PBC need the lattice encoded. Two common formats:
1. Comment line: `Lattice="a11 a12 a13 a21 ..."` (ASE/OVITO convention)
2. POSCAR / CIF

For now: **check if the Molecule reader already parses `Lattice=` from XYZ comments** and use that. If not, add minimal support in `src/tools/formats.h`.

Check: `grep -n "Lattice\|lattice\|unit_cell" src/tools/formats.h src/core/molecule.cpp`

### Step 6: CLI flag

Add `-pbc` flag that reads a lattice from CLI or activates PBC from the XYZ comment. Or simply: if the XYZ file has `Lattice=` in the comment, PBC is automatically active (zero-config for users who supply ASE-format XYZ).

## Critical Files

| File | Change |
|------|--------|
| `src/core/energy_calculators/ff_methods/forcefield.h` | Add unit cell members + setter |
| `src/core/energy_calculators/ff_methods/forcefield.cpp` | Implement setter, distribute to threads |
| `src/core/energy_calculators/ff_methods/forcefieldthread.h` | Add unit cell members, include pbc_utils |
| `src/core/energy_calculators/ff_methods/forcefieldthread.cpp` | Apply MIC in ~10 distance loops |
| `src/core/energy_calculators/ff_methods/gfnff.h` | Add unit cell members |
| `src/core/energy_calculators/ff_methods/gfnff_method.cpp` | Extract lattice from Molecule, propagate |
| `src/core/energy_calculators/ff_methods/cuda/ff_workspace_gpu.h` | Lattice setter, host buffers |
| `src/core/energy_calculators/ff_methods/cuda/ff_workspace_gpu.cu` | Constant memory, applyMIC(), 22 kernels |
| `src/core/energy_calculators/qm_methods/gfnff_gpu_method.cpp` | Upload lattice before GPU launch |
| `src/tools/pbc_utils.h` | Reuse as-is (no change) |
| `src/core/molecule.h` | Reuse as-is (no change) |

## Reuse (no reimplementation)

- `PBCUtils::applyMinimumImage(r_vec, cell, cell_inv)` — `src/tools/pbc_utils.h`
- `Molecule::getUnitCell()`, `hasPBC()`, `getUnitCellInverse()` — `src/core/molecule.h`

## Known Limitations (document, do not implement)

- Image enumeration for boxes smaller than cutoff radius
- Pressure coupling / variable cell (SimpleMD)
- Topology regeneration when box changes during MD (currently uses cached topology)
- No PBC for XB parameter generation (halogen bonds across cell faces)

## Verification

1. **Unit test**: Create a water dimer in a box of 10×10×10 Å. Compute energy with and without PBC — should match for large box (no image contribution). Then shrink box to 5 Å — energy should change as image pairs interact.
2. **Translational invariance**: Shift all atoms by one lattice vector — energy and gradient must be unchanged.
3. **Gradient consistency**: Numerical gradient (finite difference with PBC) must match analytical gradient with PBC.
4. **SimpleMD run**: Run 100 steps of water box — energy should be conserved, no atoms escape the box.
5. **GPU vs CPU**: Energy and gradient from GPU path must match CPU path for periodic system.
6. Build: `cd release && make -j4`; Tests: `ctest --output-on-failure`
