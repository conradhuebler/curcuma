# GFN-FF CPU Optimierungsanalyse

Updated 2026-04-12. All CPU optimizations (P1a–P3b) are machine-tested.

---

## Completed Optimizations

| ID | Optimization | Status | Result |
|----|-------------|--------|--------|
| P1a | D4 Gaussian weights threshold cache | ⚙️ Machine-tested | CPU: 1.01× (156.3→157.6 ms/step) — Cache rarely active at 1410 atoms. D4 term is small fraction of total time. |
| P1b | Round-robin thread distribution | ⚙️ Machine-tested | CPU: 1.04× (157.7→151.1 ms/step) — Effect minimal at 1410 atoms: cache lines already large enough. Expected benefit at >4 threads or unbalanced loads. |
| P1c | Shrink GFNFFDispersion struct (88→48 bytes) | ⚙️ Machine-tested | CPU: 1.00× (157.9 vs 157.7 ms/step) — No measurable effect at 1410 atoms. Cache limit not reached. |
| P2a | Eliminate N×N distance matrix | ⚙️ Machine-tested | CPU: 1.05× (157.7→150.5 ms/step) — distance_matrix removed from per-step path, squared_dist_matrix eliminated (dead weight). CNCalculator replaces D3-style CN. Memory savings: ~32 MB at 1410 atoms. |
| P2b | CN neighbor-list with configurable cutoff | ⚙️ Machine-tested | CPU: 1.08× (157.7→146.6 ms/step, 6 Bohr cutoff) — Reference mode (full O(N²)): 159.6 ms/step. Parameter: `-gfnff.cn_cutoff_bohr 6` (default). |
| P3a | Shrink GFNFFCoulomb struct (120→56 bytes) | ⚙️ Machine-tested | CPU: 0.99× (159.8 vs 157.7 ms/step) — Redundant fields removed (q_i/j, chi_i/j, gam_i/j, alp_i/j). Per-atom vectors instead of per-pair data. |
| P3b | `std::set` → bitset bonded exclusion | ⚙️ Machine-tested | CPU: 0.99× (159.8 vs 157.7 ms/step) — m_bonded_pairs dead code removed, `std::set`→`vector<bool>` bitset in `generateRepulsionPairsNative()`. |

**Summary:** At 1410 atoms, CPU optimizations are within measurement noise (±2%). Benefits (cache locality, memory savings) will appear at larger system sizes (N>2000) where O(N²) terms dominate cache capacity.

---

## Performance Summary (2026-04-12)

**System:** AMD Ryzen 9 9950X3D, NVIDIA GeForce RTX 5080, 4 Threads
**Testfall:** Polymer, 1410 Atome, 1000 MD-Schritte

| Methode | ms/Schritt | Faktor vs. Baseline |
|---------|-----------|---------------------|
| gfnff (CPU, nach P1a–P3b) | 148.8 | 1.06× |
| gfnff (CPU, Baseline) | 157.7 | 1.00× |
| gfnff-gpu (nach G1a+b+G1c+G2a) | 19.0 | GPU 7.8× faster than CPU |

---

## Not Recommended

- OpenMP/SIMD pragmas: CxxThreadPool + Eigen-SIMD sufficient, no measurable benefit
- Cutoff reduction: Physically risky without reference validation
- CUTLASS / Tensor Cores: GFN-FF is not a matrix problem
- JSON roundtrip elimination: Native `GFNFFParameterSet` path already exists (`forcefield.cpp:449`)

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
