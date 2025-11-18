# Multi-CV Metadynamics Examples

**Claude Generated (November 2025)**

This document provides practical examples for using the Multi-CV Metadynamics (CV-MTD) framework in Curcuma SimpleMD.

## Table of Contents

1. [Basic 1D Metadynamics](#basic-1d-metadynamics)
2. [2D Metadynamics](#2d-metadynamics)
3. [3D Metadynamics](#3d-metadynamics)
4. [Available CV Types](#available-cv-types)
5. [Parameter Reference](#parameter-reference)
6. [Output Files](#output-files)

---

## Basic 1D Metadynamics

**Scenario**: Enhance sampling along a distance between two atoms (e.g., bond breaking/forming).

### Input Structure

Create a simple molecule (e.g., `ethane.xyz`):
```
8
Ethane molecule
C   0.0   0.0   0.0
C   1.5   0.0   0.0
H  -0.5   0.9   0.0
H  -0.5  -0.9   0.0
H  -0.5   0.0   0.9
H   2.0   0.9   0.0
H   2.0  -0.9   0.0
H   2.0   0.0   0.9
```

### Command Line

```bash
./curcuma -md ethane.xyz \
  -method uff \
  -temperature 300 \
  -max_time 50000 \
  -time_step 1.0 \
  -cv_mtd true \
  -cv1_type distance \
  -cv1_atoms "0,1" \
  -cv1_sigma 0.05 \
  -cv_mtd_height 1.0 \
  -cv_mtd_pace 500 \
  -cv_mtd_well_tempered true \
  -cv_mtd_bias_factor 10.0 \
  -verbosity 2
```

### Explanation

- **`-cv_mtd true`**: Enable multi-CV metadynamics
- **`-cv1_type distance`**: First CV is a distance
- **`-cv1_atoms "0,1"`**: Between atoms 0 and 1 (C-C bond)
- **`-cv1_sigma 0.05`**: Gaussian width 0.05 Å (narrow for precise sampling)
- **`-cv_mtd_height 1.0`**: Initial Gaussian height 1.0 kJ/mol
- **`-cv_mtd_pace 500`**: Deposit Gaussian every 500 MD steps
- **`-cv_mtd_well_tempered true`**: Use well-tempered MTD (converges)
- **`-cv_mtd_bias_factor 10.0`**: Bias factor γ = 10 (controls convergence speed)

### Expected Output

- **`COLVAR`**: CV values and bias energy at each timestep
  ```
  # time(fs)  CV1_distance  bias_energy(kJ/mol)
  0.0         1.52          0.0
  1.0         1.51          0.0
  ...
  ```

- **`HILLS`**: Deposited Gaussians for restarting
  ```
  # time(fs)  CV1  height  sigma1
  500.0       1.53  1.000  0.050
  1000.0      1.49  0.950  0.050
  ...
  ```

---

## 2D Metadynamics

**Scenario**: Sample conformational space of alanine dipeptide using φ and ψ dihedral angles (Ramachandran plot).

### Input Structure

Download alanine dipeptide structure or use:
```
22
Alanine dipeptide
C  -2.5   0.0   0.0
...
(full structure omitted for brevity)
```

### Command Line

```bash
./curcuma -md alanine_dipeptide.xyz \
  -method gfn2 \
  -temperature 300 \
  -max_time 100000 \
  -time_step 0.5 \
  -cv_mtd true \
  -cv1_type dihedral \
  -cv1_atoms "4,6,8,14" \
  -cv1_sigma 0.3 \
  -cv1_periodic true \
  -cv1_period 360.0 \
  -cv2_type dihedral \
  -cv2_atoms "6,8,14,16" \
  -cv2_sigma 0.3 \
  -cv2_periodic true \
  -cv2_period 360.0 \
  -cv_mtd_height 0.5 \
  -cv_mtd_pace 1000 \
  -cv_mtd_well_tempered true \
  -cv_mtd_bias_factor 15.0 \
  -verbosity 2
```

### Explanation

- **Two CVs**: φ (phi) and ψ (psi) dihedral angles
- **`cv1_atoms "4,6,8,14"`**: Atoms defining φ dihedral (C-N-Cα-C)
- **`cv2_atoms "6,8,14,16"`**: Atoms defining ψ dihedral (N-Cα-C-N)
- **`cv1_periodic true`**: Dihedrals are periodic (-180° to +180°)
- **`cv1_period 360.0`**: Full period is 360 degrees
- **Lower height (0.5 kJ/mol)**: 2D requires smaller Gaussians (curse of dimensionality)
- **Slower pace (1000 steps)**: More time to explore before adding bias

### Post-Processing

Create 2D free energy surface (FES):
```bash
# Using plumed sum_hills (if PLUMED is installed)
plumed sum_hills --hills HILLS --outfile fes.dat --mintozero

# Or manually bin COLVAR data to create histogram
```

Expected FES shows canonical Ramachandran regions (α-helix, β-sheet, PPII).

---

## 3D Metadynamics

**Scenario**: Sample complex molecular rearrangement using 3 distance CVs.

### Command Line

```bash
./curcuma -md complex_molecule.xyz \
  -method uff \
  -temperature 350 \
  -max_time 200000 \
  -time_step 1.0 \
  -cv_mtd true \
  -cv1_type distance \
  -cv1_atoms "0,5" \
  -cv1_sigma 0.1 \
  -cv2_type distance \
  -cv2_atoms "5,10" \
  -cv2_sigma 0.1 \
  -cv3_type distance \
  -cv3_atoms "10,15" \
  -cv3_sigma 0.1 \
  -cv_mtd_height 0.3 \
  -cv_mtd_pace 2000 \
  -cv_mtd_well_tempered true \
  -cv_mtd_bias_factor 20.0 \
  -verbosity 1
```

### Important Notes for 3D

⚠️ **Curse of Dimensionality**: Number of Gaussians required scales as (L/σ)^d, where d=3.
- Use **smaller heights** (0.1-0.5 kJ/mol)
- Use **slower pace** (>1000 steps)
- Use **larger bias_factor** (15-25) for well-tempered MTD
- Consider **grid acceleration** for >500 Gaussians (future feature)

---

## Available CV Types

### 1. Distance

**Description**: Distance between two atoms
**Atoms required**: 2
**Example**: `-cv1_type distance -cv1_atoms "0,5"`
**Typical σ**: 0.05-0.1 Å
**Units**: Angstroms (Å)

### 2. Angle

**Description**: Valence angle between three atoms
**Atoms required**: 3
**Example**: `-cv1_type angle -cv1_atoms "0,5,10"`
**Typical σ**: 5-10 degrees
**Units**: Degrees (°)

### 3. Dihedral

**Description**: Torsion angle between four atoms
**Atoms required**: 4
**Example**: `-cv1_type dihedral -cv1_atoms "0,5,10,15"`
**Typical σ**: 10-20 degrees
**Periodic**: Usually yes (`-cv1_periodic true -cv1_period 360.0`)
**Units**: Degrees (°)

### 4. Gyration Radius

**Description**: Radius of gyration (compactness measure)
**Atoms required**: Arbitrary (list of atoms to include)
**Example**: `-cv1_type gyration -cv1_atoms "0,1,2,3,4,5"`
**Typical σ**: 0.1-0.2 Å
**Units**: Angstroms (Å)

### 5. Coordination Number

**Description**: Smooth coordination number between two groups
**Atoms required**: Two groups (A and B) or self-coordination
**Example**: `-cv1_type coordination -cv1_atoms "0,1,2" -cv1_atoms_B "10,11,12"`
**Typical σ**: 0.5-1.0
**Units**: Dimensionless

**Special parameters**:
- `-cv1_cutoff 3.0`: Distance cutoff (Å)
- `-cv1_n 6`: Numerator exponent
- `-cv1_m 12`: Denominator exponent

---

## Parameter Reference

### CV-MTD Global Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `cv_mtd` | Bool | false | Enable multi-CV metadynamics |
| `cv_mtd_pace` | Int | 500 | Deposit Gaussian every N steps |
| `cv_mtd_height` | Double | 1.0 | Initial Gaussian height (kJ/mol) |
| `cv_mtd_well_tempered` | Bool | true | Use well-tempered MTD |
| `cv_mtd_bias_factor` | Double | 10.0 | Well-tempered bias factor γ |
| `cv_mtd_write_stride` | Int | 1 | Write COLVAR/HILLS every N steps |

### CV Definition Parameters (1-3 CVs)

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `cv1_type` | String | "none" | CV type: distance/angle/dihedral/gyration/coordination |
| `cv1_atoms` | String | "" | Comma-separated atom indices (e.g., "0,5") |
| `cv1_sigma` | Double | 0.1 | Gaussian width for CV1 |
| `cv1_periodic` | Bool | false | Enable periodicity (for dihedrals) |
| `cv1_period` | Double | 360.0 | Period if periodic (degrees) |

*Repeat for `cv2_*` and `cv3_*` for additional CVs.*

### Grid Acceleration Parameters (November 2025 - Claude Generated)

Grid acceleration provides **O(1)** bias evaluation instead of O(N_gaussians) for large simulations.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `cv_mtd_grid_auto` | Bool | true | Auto-enable grid when Gaussian count reaches threshold |
| `cv_mtd_grid_threshold` | Int | 500 | Gaussian count to trigger grid acceleration |
| `cv_mtd_grid_bins` | Int | 100 | Grid bins per CV dimension |

**When to Use Grid Acceleration**:
- **Always recommended**: Grid auto-enables at 500 Gaussians by default
- **1D simulations**: Negligible overhead (~1 KB memory)
- **2D simulations**: Moderate memory (~80 KB for 100×100 grid)
- **3D simulations**: Significant memory (~8 MB for 100×100×100 grid)

**Performance Gains**:
- Without grid: 1000 Gaussians = 1000× slower per MD step
- With grid: O(1) per MD step regardless of Gaussian count
- Break-even: ~100-500 Gaussians (dimensionality dependent)

**Example (disable auto-grid)**:
```bash
./curcuma -md system.xyz \
  -cv_mtd true \
  -cv1_type distance -cv1_atoms "0,5" \
  -cv_mtd_grid_auto false  # Disable auto, use direct summation always
```

**Example (manual grid threshold)**:
```bash
./curcuma -md system.xyz \
  -cv_mtd true \
  -cv1_type distance -cv1_atoms "0,5" \
  -cv_mtd_grid_threshold 200  # Enable grid at 200 Gaussians (earlier)
  -cv_mtd_grid_bins 150       # Higher resolution grid
```

### Choosing Parameters

**Gaussian Height (`cv_mtd_height`)**:
- **1D**: 0.5-2.0 kJ/mol (aggressive sampling)
- **2D**: 0.3-1.0 kJ/mol (moderate)
- **3D**: 0.1-0.5 kJ/mol (conservative)

**Gaussian Width (`cvN_sigma`)**:
- Should cover ~10-20% of CV range
- Too small: Many Gaussians required (slow convergence)
- Too large: Poor resolution, bias "leaks" across barriers

**Deposition Pace (`cv_mtd_pace`)**:
- Fast (100-500 steps): Aggressive exploration
- Medium (500-1000 steps): Balanced
- Slow (1000-5000 steps): Conservative, better accuracy

**Bias Factor (`cv_mtd_bias_factor`)**:
- **Standard MTD**: Set `cv_mtd_well_tempered false` (γ → ∞)
- **Well-tempered MTD**: γ = 5-50
  - Lower γ (5-10): Fast convergence, moderate accuracy
  - Higher γ (15-50): Slower, better free energy accuracy
- Relates to fictitious temperature: T + ΔT, where ΔT = T(γ-1)

---

## Output Files

### COLVAR

**Format**: Space-separated columns
**Contents**: Timestep, CV values, bias energy

```
# time(fs)  CV1  CV2  bias_energy(kJ/mol)
0.0         1.52  120.5  0.0
1.0         1.51  121.0  0.023
2.0         1.50  119.8  0.045
...
```

**Usage**: Visualize CV evolution, check sampling quality

### HILLS

**Format**: Space-separated columns
**Contents**: Deposition time, CV center, height, widths

```
# time(fs)  CV1    CV2     height  sigma1  sigma2
500.0       1.53   120.3   1.000   0.050   10.0
1000.0      1.49   118.7   0.950   0.050   10.0
1500.0      1.55   122.1   0.903   0.050   10.0
...
```

**Usage**: Restart simulations, reconstruct FES, debugging

### Trajectory Files

- **`basename.xyz`**: Full MD trajectory
- **`basename.mtd.xyz`**: Structures at Gaussian deposition times (legacy RMSD-MTD)

---

## Tips and Best Practices

### 1. Start with 1D

Always test new systems with 1D metadynamics first:
- Verify CV behaves as expected
- Tune parameters (height, sigma, pace)
- Check for numerical stability

### 2. Validate with Known Systems

Test on systems with known free energy profiles:
- Alanine dipeptide (Ramachandran plot)
- Butane torsion (3-fold barrier)
- Simple bond stretching

### 3. Monitor Convergence

Check COLVAR file:
- CV should visit all relevant regions
- Bias energy should plateau (well-tempered MTD)
- Oscillations indicate insufficient sampling

### 4. Restart Simulations

To continue a simulation:
```bash
# Save HILLS file from previous run
cp HILLS HILLS_previous

# Restart (not yet implemented - future feature)
./curcuma -md ... -cv_mtd_restart HILLS_previous
```

### 5. Post-Processing

Reconstruct free energy surface:
```python
import numpy as np
import matplotlib.pyplot as plt

# Load HILLS file
hills = np.loadtxt('HILLS', skiprows=1)
time = hills[:,0]
cv1 = hills[:,1]
cv2 = hills[:,2]  # If 2D
height = hills[:,3]

# Create 2D histogram
H, xedges, yedges = np.histogram2d(cv1, cv2, bins=50)

# Plot
plt.contourf(H.T, extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
plt.xlabel('CV1')
plt.ylabel('CV2')
plt.colorbar(label='Free Energy (kJ/mol)')
plt.show()
```

---

## Troubleshooting

### Problem: Simulation crashes with "CV gradient NaN"

**Cause**: Atoms in collinear or overlapping configuration
**Solution**:
- Add initial equilibration without bias
- Check CV atom indices are correct
- Use softer force field (GFN-FF instead of UFF)

### Problem: Bias energy grows without bound

**Cause**: `cv_mtd_well_tempered false` (standard MTD)
**Solution**: Enable well-tempered: `-cv_mtd_well_tempered true`

### Problem: Poor sampling, CV stuck in one region

**Cause**: Barriers too high or Gaussians too small
**Solution**:
- Increase height: `-cv_mtd_height 2.0`
- Decrease pace: `-cv_mtd_pace 200`
- Increase temperature: `-temperature 400`

### Problem: Compilation error "CVFactory::create not found"

**Cause**: Missing include or incorrect CV type string
**Solution**: Check `cv1_type` spelling (lowercase: distance, angle, dihedral, gyration, coordination)

---

## References

1. **Laio & Parrinello (2002)**: Original metadynamics paper
   DOI: 10.1073/pnas.202427399

2. **Barducci et al. (2008)**: Well-tempered metadynamics
   DOI: 10.1103/PhysRevLett.100.020603

3. **Bonomi et al. (2009)**: PLUMED plugin (architecture inspiration)
   DOI: 10.1016/j.cpc.2009.05.011

4. **Iannuzzi et al. (2003)**: Coordination number CV
   DOI: 10.1103/PhysRevLett.90.238302

5. **Blondel & Karplus (1996)**: Dihedral angle gradients
   DOI: 10.1002/(SICI)1096-987X(199604)17:5/6<689::AID-JCC10>3.0.CO;2-T

---

**For more information, see**:
- `docs/COLLECTIVE_VARIABLES.md` - Mathematical formulations and theory
- `docs/MULTI_CV_METADYNAMICS.md` - Algorithm details and implementation
- `curcuma -md -help` - Full parameter list

**Report issues**: https://github.com/conradhuebler/curcuma/issues
