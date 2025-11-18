# Multi-Dimensional Metadynamics with Multiple Collective Variables

**Author:** Claude (Anthropic) & Conrad Hübler
**Date:** November 2025
**Status:** Implementation Guide

---

## Table of Contents

1. [Motivation](#motivation)
2. [Theoretical Background](#theoretical-background)
3. [Mathematical Formulation](#mathematical-formulation)
4. [Implementation Strategy](#implementation-strategy)
5. [Curse of Dimensionality](#curse-of-dimensionality)
6. [Parameter Selection Guidelines](#parameter-selection-guidelines)
7. [Output Format](#output-format)
8. [References](#references)

---

## 1. Motivation

### Why Multiple CVs?

**Problem:** A single collective variable (CV) often cannot describe complex molecular processes.

**Examples where multiple CVs are essential:**

1. **Protein-ligand binding:**
   - CV1: Distance (ligand-protein COM)
   - CV2: Radius of gyration (ligand compactness)
   - CV3: Number of contacts (coordination)

2. **Conformational transitions:**
   - CV1: Dihedral φ (backbone)
   - CV2: Dihedral ψ (backbone)
   - CV3: RMSD to native structure

3. **Chemical reactions:**
   - CV1: Bond distance (reactant)
   - CV2: Bond distance (product)
   - CV3: Coordination number (catalyst)

**Key Insight:** Multi-CV metadynamics explores a **multi-dimensional free energy landscape** F(s₁, s₂, ..., s_d) where d is the number of CVs.

---

## 2. Theoretical Background

### 2.1 Free Energy Landscape

For a system described by collective variables **s** = (s₁, s₂, ..., s_d), the free energy is:

```
F(s) = -k_B T * ln(P(s))
```

where:
- P(s) is the probability density in CV space
- k_B is Boltzmann's constant
- T is temperature

### 2.2 Metadynamics Bias Potential

In metadynamics, we add a time-dependent bias potential V(s, t) that fills the free energy landscape:

```
V(s, t) = Σ_i w_i * exp(-Σ_α (s_α - s_α,i)² / (2σ_α²))
```

where:
- **s** = (s₁, s₂, ..., s_d): Current CV values
- **s_i** = (s₁,i, s₂,i, ..., s_d,i): CV values when i-th Gaussian was deposited
- w_i: Height of i-th Gaussian
- σ_α: Width of Gaussian along CV α
- α runs from 1 to d (number of CVs)

**Key Point:** This is a **product of independent Gaussians** along each CV dimension:

```
V(s, t) = Σ_i w_i * ∏_α exp(-(s_α - s_α,i)² / (2σ_α²))
       = Σ_i w_i * exp(-(s₁ - s₁,i)²/(2σ₁²)) *
                  exp(-(s₂ - s₂,i)²/(2σ₂²)) *
                  ... *
                  exp(-(s_d - s_d,i)²/(2σ_d²))
```

### 2.3 Well-Tempered Metadynamics (Multi-CV)

For well-tempered MTD, the Gaussian height decreases over time:

```
w_i = w_0 * exp(-V(s_i, t_i) / (k_B * ΔT))
```

where:
- w_0: Initial Gaussian height
- ΔT: Bias temperature (controls tempering)
- V(s_i, t_i): Current bias at position s_i when depositing Gaussian i

**Convergence Property:** As t → ∞, the bias converges to:

```
lim_{t→∞} V(s, t) = -(T + ΔT)/ΔT * F(s) + C
```

This allows direct reconstruction of the free energy surface!

---

## 3. Mathematical Formulation

### 3.1 Bias Force Calculation

The bias force on atom j is computed via the **multidimensional chain rule**:

```
F_j^bias = -∇_r_j V(s, t)
        = -Σ_α (∂V/∂s_α) * (∂s_α/∂r_j)
```

**Step-by-step derivation:**

1. **Compute ∂V/∂s_α** (derivative of bias w.r.t. CV α):
   ```
   ∂V/∂s_α = Σ_i w_i * exp(-Σ_β (s_β - s_β,i)²/(2σ_β²)) *
                      [-(s_α - s_α,i) / σ_α²]
   ```

2. **Compute ∂s_α/∂r_j** (gradient of CV α w.r.t. atom j position):
   - This comes from the CV implementation (e.g., Distance CV, Angle CV, etc.)
   - See `docs/COLLECTIVE_VARIABLES.md` for individual CV gradients

3. **Sum over all CVs**:
   ```
   F_j^bias = -Σ_α [∂V/∂s_α * ∂s_α/∂r_j]
   ```

**Computational Complexity:** O(N_gaussians * N_cvs * N_atoms)

### 3.2 Efficient Computation

**Optimization 1: Factor out common exponentials**

Since `exp(-Σ_β (s_β - s_β,i)²/(2σ_β²)) = ∏_β exp(-(s_β - s_β,i)²/(2σ_β²))`:

```cpp
// Precompute for each Gaussian i
double total_exp = 1.0;
for (int alpha = 0; alpha < N_cvs; ++alpha) {
    double ds = s_alpha - s_alpha_i;
    total_exp *= exp(-ds * ds / (2 * sigma_alpha * sigma_alpha));
}

// Then for each CV alpha
for (int alpha = 0; alpha < N_cvs; ++alpha) {
    double dV_ds = -w_i * total_exp * (s_alpha - s_alpha_i) / (sigma_alpha * sigma_alpha);
    // ... apply chain rule with CV gradient
}
```

**Optimization 2: Early termination**

If any Gaussian is far from current position in ANY CV dimension:
```cpp
bool skip_gaussian = false;
for (int alpha = 0; alpha < N_cvs; ++alpha) {
    if (abs(s_alpha - s_alpha_i) > 6.0 * sigma_alpha) {
        skip_gaussian = true;  // exp(-36) ≈ 10^-16 (negligible)
        break;
    }
}
```

---

## 4. Implementation Strategy

### 4.1 Data Structure Changes

**Current (single CV - RMSD only):**
```cpp
struct BiasStructure {
    Geometry geometry;        // Single RMSD CV
    double rmsd_reference;    // CV value
    double time, energy, factor;
};
```

**New (multiple CVs):**
```cpp
struct BiasStructure {
    std::vector<double> cv_values;  // s = (s₁, s₂, ..., s_d)
    double time, energy, factor;
    int counter;
};

class BiasThread {
private:
    std::vector<CVPtr> m_cvs;  // Multiple CVs
    std::vector<double> m_sigmas;  // σ_α for each CV
    std::vector<BiasStructure> m_biased_structures;
};
```

### 4.2 Bias Calculation Algorithm

**Pseudocode:**
```
function ComputeBias(current_molecule):
    // Step 1: Evaluate all CVs
    for alpha = 1 to N_cvs:
        s_alpha = cvs[alpha].calculate(current_molecule)

    // Step 2: Compute bias energy
    V_bias = 0
    for i = 1 to N_gaussians:
        // Compute product of Gaussians
        total_exp = 1.0
        for alpha = 1 to N_cvs:
            ds = s_alpha - s_alpha_i
            total_exp *= exp(-ds^2 / (2 * sigma_alpha^2))

        // Add contribution
        V_bias += w_i * total_exp

    // Step 3: Compute gradients
    gradient = zeros(N_atoms, 3)
    for alpha = 1 to N_cvs:
        dV_ds_alpha = 0
        for i = 1 to N_gaussians:
            // Same total_exp as above
            dV_ds_alpha += -w_i * total_exp * (s_alpha - s_alpha_i) / sigma_alpha^2

        // Chain rule: multiply by CV gradient
        cv_gradient = cvs[alpha].gradient(current_molecule)
        gradient += dV_ds_alpha * cv_gradient

    return V_bias, gradient
```

---

## 5. Curse of Dimensionality

### 5.1 The Problem

**Key Challenge:** The number of Gaussians needed to fill CV space grows **exponentially** with dimensionality d:

```
N_gaussians ≈ (L / σ)^d
```

where:
- L: Range of each CV
- σ: Gaussian width
- d: Number of CVs

**Example (catastrophic scaling):**
- 1D: L=10 Å, σ=0.5 Å → N ≈ 20 Gaussians
- 2D: N ≈ 400 Gaussians
- 3D: N ≈ 8,000 Gaussians
- 4D: N ≈ 160,000 Gaussians (INFEASIBLE!)

### 5.2 Practical Limits

**RECOMMENDATION:**
- **1 CV:** Always fine
- **2 CVs:** Good, most common in literature
- **3 CVs:** Challenging, requires careful selection
- **4+ CVs:** Usually impractical, consider alternative methods

**Alternative for >3 CVs:**
- Path collective variables (s, z)
- Committor analysis
- Transition path sampling
- Harmonic restraints (not full MTD)

### 5.3 Mitigating Strategies

1. **Well-tempered MTD:** Reduces over-filling (recommended)
2. **Grid acceleration:** Interpolate bias on grid (>100x faster)
3. **Variable bias width:** Larger σ in less important regions
4. **Bias exchange MTD:** Multiple replicas with different CVs
5. **Careful CV selection:** Choose minimal set that describes process

---

## 6. Parameter Selection Guidelines

### 6.1 Gaussian Widths (σ_α)

**Rule of Thumb:**
```
σ_α ≈ (1/3 to 1/2) * ΔCVΑ
```

where ΔCVΑ is the expected fluctuation of CV α.

**How to estimate ΔCVΑ:**
1. Run short unbiased MD (100-1000 steps)
2. Compute standard deviation of CV α
3. Set σ_α = 0.5 * std(CV_α)

**Example (alanine dipeptide, φ and ψ dihedrals):**
- ΔCVΦ ≈ 30° → σ_φ = 15° = 0.26 rad
- ΔCVΨ ≈ 30° → σ_ψ = 15° = 0.26 rad

### 6.2 Gaussian Heights (w)

**Well-tempered MTD (RECOMMENDED):**
```
w_0 = k_B * T / (20 to 100)
```

**Standard MTD:**
```
w = 0.5 to 2.0 kJ/mol  (depends on barrier height)
```

**ΔT (bias temperature, WT-MTD only):**
```
ΔT = (5 to 10) * T
```

Higher ΔT → Faster convergence but less accurate
Lower ΔT → Slower convergence but more accurate

### 6.3 Deposition Rate

**How often to add new Gaussian:**
```
Pace = 100 to 500 MD steps
```

**Trade-off:**
- Fast (50 steps): Rapid exploration, may overshoot barriers
- Slow (1000 steps): Accurate, but slow convergence

### 6.4 Simulation Length

**Convergence criterion:**
Run until the free energy difference between consecutive time windows is below threshold:

```
ΔF_converged < k_B * T  (≈ 0.6 kcal/mol at 300 K)
```

**Typical times (2D metadynamics):**
- Simple system (alanine dipeptide): 10-50 ns
- Complex system (protein folding): 100 ns - 1 μs

---

## 7. Output Format

### 7.1 COLVAR File (Multiple CVs)

**Format (space-separated):**
```
# time  s1  s2  ...  s_d  V_bias  N_gaussians
0.0     1.5  120.0  ...  -30.0   0.0     0
0.5     1.6  118.0  ...  -28.0   0.1     1
1.0     1.7  115.0  ...  -25.0   0.2     1
...
```

**Columns:**
1. Time (ps)
2-[d+1]. CV values (s₁, s₂, ..., s_d)
[d+2]. Bias potential V(s, t)
[d+3]. Number of Gaussians deposited so far

### 7.2 HILLS File (Gaussian Parameters)

**Format:**
```
# time  s1_center  s2_center  ...  sigma1  sigma2  ...  height  bias_factor
0.5     1.6        118.0      ...  0.2     15.0    ...  0.5     1.0
1.0     1.7        115.0      ...  0.2     15.0    ...  0.49    0.98
...
```

### 7.3 Free Energy Surface (Post-Processing)

**1D FES (if only 1 CV):**
```
# s1  F(s1)  error
-180  5.2    0.3
-170  4.8    0.2
...
```

**2D FES (if 2 CVs):**
```
# s1   s2    F(s1,s2)  error
-180  -180  10.5      0.5
-180  -170  9.8       0.4
...
```

**File formats:**
- `.dat`: ASCII, human-readable
- `.dx`: OpenDX format (VMD compatible)
- `.cube`: Gaussian Cube format (Gaussian, Chimera)

---

## 8. References

### Core Metadynamics Papers

[1] **Laio, A. & Parrinello, M.** (2002). *Escaping free-energy minima.*
    Proc. Natl. Acad. Sci. USA **99**, 12562-12566.
    DOI: [10.1073/pnas.202427399](https://doi.org/10.1073/pnas.202427399)
    **→ Original metadynamics paper**

[2] **Barducci, A., Bussi, G. & Parrinello, M.** (2008). *Well-tempered metadynamics: A smoothly converging and tunable free-energy method.*
    Phys. Rev. Lett. **100**, 020603.
    DOI: [10.1103/PhysRevLett.100.020603](https://doi.org/10.1103/PhysRevLett.100.020603)
    **→ Well-tempered MTD (RECOMMENDED method)**

[3] **Valsson, O., Tiwary, P. & Parrinello, M.** (2016). *Enhancing important fluctuations: Rare events and metadynamics from a conceptual viewpoint.*
    Annu. Rev. Phys. Chem. **67**, 159-184.
    DOI: [10.1146/annurev-physchem-040215-112229](https://doi.org/10.1146/annurev-physchem-040215-112229)
    **→ Comprehensive review of modern MTD variants**

### Multi-CV Metadynamics

[4] **Bussi, G., Laio, A. & Parrinello, M.** (2006). *Equilibrium free energies from nonequilibrium metadynamics.*
    Phys. Rev. Lett. **96**, 090601.
    DOI: [10.1103/PhysRevLett.96.090601](https://doi.org/10.1103/PhysRevLett.96.090601)
    **→ Multi-dimensional MTD convergence**

[5] **Branduardi, D., Gervasio, F. L. & Parrinello, M.** (2007). *From A to B in free energy space.*
    J. Chem. Phys. **126**, 054103.
    DOI: [10.1063/1.2432340](https://doi.org/10.1063/1.2432340)
    **→ Path collective variables (alternative to multi-CV)**

### Applications

[6] **Ensing, B. et al.** (2006). *Metadynamics as a tool for exploring free energy landscapes of chemical reactions.*
    Acc. Chem. Res. **39**, 73-81.
    DOI: [10.1021/ar040198i](https://doi.org/10.1021/ar040198i)
    **→ Chemical reactions with 2-3 CVs**

[7] **Stirling, A., Iannuzzi, M., Laio, A. & Parrinello, M.** (2004). *Azobenzene at coinage metal surfaces: Role of dispersive van der Waals interactions.*
    ChemPhysChem **5**, 1558-1568.
    DOI: [10.1002/cphc.200400063](https://doi.org/10.1002/cphc.200400063)
    **→ Multiple CVs for surface chemistry**

### Software

[8] **Bonomi, M. et al.** (2009). *PLUMED: A portable plugin for free-energy calculations with molecular dynamics.*
    Comp. Phys. Comm. **180**, 1961-1972.
    DOI: [10.1016/j.cpc.2009.05.011](https://doi.org/10.1016/j.cpc.2009.05.011)
    **→ PLUMED implementation (reference for multi-CV)**

[9] **Tribello, G. A. et al.** (2014). *PLUMED 2: New feathers for an old bird.*
    Comp. Phys. Comm. **185**, 604-613.
    DOI: [10.1016/j.cpc.2013.09.018](https://doi.org/10.1016/j.cpc.2013.09.018)
    **→ PLUMED 2 architecture**

---

## Appendix A: Example - 2D Metadynamics for Alanine Dipeptide

**System:** Alanine dipeptide (small protein model)

**CVs:**
- CV1: φ (backbone dihedral C-N-CA-C)
- CV2: ψ (backbone dihedral N-CA-C-N)

**Parameters:**
- σ_φ = 0.35 rad (≈ 20°)
- σ_ψ = 0.35 rad (≈ 20°)
- w_0 = 1.2 kJ/mol
- ΔT = 3000 K (WT-MTD)
- Pace = 500 steps (1 ps)
- Total time = 50 ns

**Result:**
Free energy surface F(φ, ψ) shows:
- **C7eq** (φ ≈ -80°, ψ ≈ +80°): Global minimum
- **C7ax** (φ ≈ +60°, ψ ≈ -60°): Metastable state
- **α_R** (φ ≈ -60°, ψ ≈ -40°): Helical region

**Barrier heights:**
- C7eq → C7ax: ~10 kJ/mol
- C7eq → α_R: ~15 kJ/mol

**Reference:** Laio & Parrinello (2002), Fig. 3

---

**End of Document**
