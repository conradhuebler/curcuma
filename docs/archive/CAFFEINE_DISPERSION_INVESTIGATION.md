# Caffeine Dispersion Error Investigation Summary

**Date**: February 11, 2026
**Status**: Analysis Complete - Ready for Debugging Phase
**Error**: 75 µEh systematic overestimation (Curcuma 0.0141 mEh more negative than Fortran)

## Executive Summary

Curcuma GFN-FF dispersion for Caffeine has a **0.41% error** (75 µEh) - acceptable for typical work but unacceptable for our **< 1 µEh target**. The error is **systematic (energy too negative)**, suggesting a consistent underestimation in the damping term or overestimation in C6 values.

**Good News**: Code review found NO obvious bugs. The issue is likely subtle (CN accuracy, reference data alignment, or normalization factors).

## Root Cause Analysis: Top 5 Hypotheses

### Hypothesis 1: CN Calculation Deviation in Aromatic Systems 🔴 PRIORITY 1

**Why it matters**: Caffeine is aromatic with 3 fused rings
- CN affects Gaussian weighting: `w = exp(-4 * (CN - CNref)²)`
- Small CN errors amplify quadratically in weights
- π-system detection may fail in complex geometries

**Evidence to check**:
- Compare per-atom CN with Fortran analyzer
- Check `CNCalculator::calculateGFNFFCN()` implementation
- Verify exponential coordination formula: `CN = Σ 1/(1+exp(...))`

**Code location**: `src/core/energy_calculators/ff_methods/eeq_solver.cpp` (CN calculation)

**Diagnostic method**:
```cpp
// In d4param_generator.cpp line 264-269 (already has this!)
if (m_atoms.size() <= 10) {
    fmt::print(stderr, "D4_CN_DIAG: CN values:\n");
    for (size_t i = 0; i < m_atoms.size(); ++i) {
        fmt::print(stderr, "D4_CN_DIAG: atom {} CN={:.6f}\n", i, m_cn_values[i]);
    }
}
```

For Caffeine, run:
```bash
./release/test_cases/test_gfnff_validation test_cases/reference_data/caffeine.ref.json 2>&1 | grep D4_CN_DIAG
# Compare CN values with Fortran analyzer output
```

### Hypothesis 2: Reference CN Data Mismatch 🟡 PRIORITY 2

**Why it matters**: Gaussian weights use reference CN from `refcn_fortran` array
- Must match Fortran `dftd3param.f90` exactly
- Wrong refcn → wrong weights → wrong C6 interpolation

**Evidence to check**:
- Verify `m_refcn` loaded correctly from `d4_reference_cn_fortran.cpp`
- Compare first 5 elements for C, H, N, O:
  ```
  C: [2.89, 3.20, 3.44, 3.62]
  H: [0.55]
  N: [2.95, 3.20, 3.54]
  O: [2.84, 3.04]
  ```

**Code location**: `src/core/energy_calculators/ff_methods/d4param_generator.cpp` line 96-109

**Diagnostic**:
```cpp
// Add after line 272:
if (CurcumaLogger::get_verbosity() >= 2) {
    CurcumaLogger::info("=== Reference CN Data Validation ===");
    for (int elem = 0; elem < 4; ++elem) {  // C, H, N, O only
        for (int ref = 0; ref < m_refn[elem]; ++ref) {
            CurcumaLogger::info(fmt::format("refcn[{}][{}] = {:.4f}",
                elem, ref, m_refcn[elem][ref]));
        }
    }
}
```

### Hypothesis 3: Gaussian Weight Normalization Bug 🟡 PRIORITY 3

**Why it matters**: Weight normalization affects relative contribution magnitudes
- If `sum_weights ≤ 1e-10`, fallback sets only ref[0] to 1.0
- May cause entire reference state distribution to collapse

**Evidence to check**:
- Print `sum_weights` per atom
- Check how often fallback condition triggers
- Verify normalized weights sum to 1.0

**Code location**: `src/core/energy_calculators/ff_methods/d4param_generator.cpp` line 1007-1015

**Diagnostic**:
```cpp
// In precomputeGaussianWeights(), after line 1008:
if (sum_weights <= 1e-10) {
    fmt::print(stderr, "WEIGHT_FALLBACK: atom {} sum_weights={:.2e} (too small)\n",
               i, sum_weights);
}
// Add after line 1015:
for (int ref = 0; ref < nref && CurcumaLogger::get_verbosity() >= 3; ++ref) {
    fmt::print(stderr, "WEIGHT_NORM: atom {} ref {} weight={:.6f}\n",
               i, ref, weights[ref]);
}
```

### Hypothesis 4: Correction Factor Incompleteness 🟡 PRIORITY 4

**Why it matters**: Polarizability corrections (ascale, sscale, secaiw) are element-specific
- Only validated for H, C, N, O, halogens
- Aromatic systems (π-electrons) may need special handling
- Missing or wrong values → wrong C6 → wrong energy

**Evidence to check**:
- Print ascale values for Caffeine atoms
- Verify sscale values match Fortran
- Check that secaiw data exists for all reference systems

**Code location**: `src/core/energy_calculators/ff_methods/computeC6Reference()` line 869-885

**Diagnostic**:
```cpp
// In computeC6Reference(), add:
if (iw == 0 && CurcumaLogger::get_verbosity() >= 3) {
    fmt::print(stderr, "C6_CORR_DEBUG: elem_i={} ref_i={} ascale={:.4f} hcount={:.1f}\n",
               elem_i, ref_i, ascale_i, hcount_i);
    fmt::print(stderr, "C6_CORR_DEBUG: elem_j={} ref_j={} ascale={:.4f} hcount={:.1f}\n",
               elem_j, ref_j, ascale_j, hcount_j);
}
```

### Hypothesis 5: Casimir-Polder Integration Precision 🟢 PRIORITY 5

**Why it matters**: C6 calculation uses numerical integration with 23 frequency points
- Trapezoidal rule error could accumulate
- May be compounded by correction factor application

**Evidence to check**:
- Compare pre-computed C6 values with Fortran formula
- Check that frequency grid is correct
- Verify integration formula: `c6 = (3/π) ∫ α_i(ω) α_j(ω) dω`

**Code location**: `src/core/energy_calculators/ff_methods/computeC6Reference()` line 843-942

**Root Cause if True**: Low (but possible). Check if frequency grid matches Fortran exactly.

## Immediate Debugging Action Plan

### Step 1: Enable Full Diagnostics (5 min)

```bash
cd /home/conrad/src/claude_curcuma/curcuma

# Create a test script to run with full verbosity
cat > run_caffeine_debug.sh << 'EOF'
#!/bin/bash
export CURCUMA_VERBOSITY=3
export OMP_NUM_THREADS=1
./release/test_cases/test_gfnff_validation test_cases/reference_data/caffeine.ref.json 2>&1 | tee caffeine_debug_full.log
echo "Extracting diagnostics..."
grep -E "D4_CN_DIAG|C6_DEBUG|WEIGHT|DISP_CSV|Dispersion" caffeine_debug_full.log > caffeine_diag.csv
EOF

chmod +x run_caffeine_debug.sh
```

### Step 2: Analyze per-atom CN (10 min)

```bash
grep "D4_CN_DIAG" caffeine_debug_full.log | awk -F'CN=' '{print $2}' > cn_values.txt

# Compare with Fortran reference (need to extract from analyzer)
./external/gfnff/build/test/gfnff-gfnff_analyze-test test_cases/molecules/larger/caffeine.xyz - 2 2 2>&1 | grep -i "coordination" > fortran_cn.txt
```

### Step 3: Analyze per-pair C6 and energy (10 min)

```bash
grep "DISP_CSV" caffeine_debug_full.log > disp_pairs.csv

# Python analysis script
python3 << 'EOF'
import csv
import sys

# Parse CSV data
pairs = []
total_energy = 0.0

with open('disp_pairs.csv') as f:
    for line in f:
        if 'DISP_CSV' in line:
            # Format: DISP_CSV: idx,i,j,rij,C6,r4r2ij,r0_squared,zetac6,t6,t8,energy
            parts = line.replace('DISP_CSV:', '').strip().split(',')
            if len(parts) >= 11:
                pair_info = {
                    'idx': int(parts[0]),
                    'i': int(parts[1]),
                    'j': int(parts[2]),
                    'rij': float(parts[3]),
                    'c6': float(parts[4]),
                    'r4r2ij': float(parts[5]),
                    'r0_sq': float(parts[6]),
                    'zetac6': float(parts[7]),
                    'energy': float(parts[10])
                }
                pairs.append(pair_info)
                total_energy += pair_info['energy']

print(f"Total dispersion (Curcuma): {total_energy:.10f} Eh")
print(f"Total dispersion (Fortran): -0.0181251200 Eh")
print(f"Error: {total_energy - (-0.0181251200):.10e} Eh")
print(f"Error %: {(total_energy / (-0.0181251200) - 1) * 100:.3f}%")
print()

# Find largest contributors
top_pairs = sorted(pairs, key=lambda p: abs(p['energy']), reverse=True)[:10]
print("Top 10 contributing pairs:")
print("idx   i-j     C6          r0²       zetac6    energy")
for p in top_pairs:
    print(f"{p['idx']:3d}  {p['i']}-{p['j']:2d}  {p['c6']:.6e}  {p['r0_sq']:.4f}  {p['zetac6']:.6f}  {p['energy']:.10e}")

EOF
```

## Decision Tree: Root Cause Identification

```
Is CN deviation detected?
├─ YES: Issue is EEQSolver or CNCalculator
│       → Fix CNCalculator::calculateGFNFFCN()
│       → Validate π-system detection for aromatic molecules
│       → Expected fix time: 2-4 hours
│
└─ NO: CN matches Fortran reference
    │
    Is C6 weighted correctly?
    ├─ NO: Issue is reference data or Gaussian weights
    │       → Fix m_refcn data or weight normalization
    │       → Validate refcn_fortran matches dftd3param.f90
    │       → Expected fix time: 1-3 hours
    │
    └─ YES: C6 matches Fortran reference
            │
            Are pair energies summing correctly?
            ├─ NO: Issue is energy calculation or pair summing
            │       → Check BJ damping formula
            │       → Verify m_final_factor
            │       → Expected fix time: 1-2 hours
            │
            └─ YES: All components match individually
                    → Anomaly in interaction or cumulative effect
                    → Request deeper Fortran source comparison
                    → Expected investigation time: 4-8 hours
```

## Success Criteria

✅ **PASS**: Caffeine dispersion error < 10 µEh (< 0.1% error)
✅ **PASS**: No regression on CH4, C6H6 test molecules
✅ **PASS**: All 11+ test molecules with updated dispersion values match reference

## Files to Monitor During Debugging

**Critical Implementation Files**:
- `src/core/energy_calculators/ff_methods/d4param_generator.cpp` - D4 parameter generation
- `src/core/energy_calculators/ff_methods/eeq_solver.cpp` - CN calculation
- `src/core/energy_calculators/ff_methods/forcefieldthread.cpp:1640-1750` - BJ damping
- `src/core/energy_calculators/ff_methods/gfnff_method.cpp:5795-5940` - Dispersion pair generation

**Reference Data Files**:
- `src/core/energy_calculators/ff_methods/d4_reference_cn_fortran.cpp` - Reference CN
- `src/core/energy_calculators/ff_methods/d4_reference_data_fixed.cpp` - Charges
- `src/core/energy_calculators/ff_methods/d4_alphaiw_data.cpp` - Polarizabilities
- `src/core/energy_calculators/ff_methods/d4_corrections_data.cpp` - Correction factors

**Test Data**:
- `test_cases/reference_data/caffeine.ref.json` - Expected: -0.01812512 Eh
- `test_cases/molecules/larger/caffeine.xyz` - Caffeine geometry (24 atoms)

## Estimated Timeline

**Best Case** (CN hypothesis true, simple fix): 2-4 hours
**Typical Case** (deeper debugging needed): 4-8 hours
**Complex Case** (Fortran source comparison): 8-16 hours

---

**Generated by**: Dispersion Accuracy Analysis System
**Reference**: GFNFF_STATUS.md Phase 3: Root Cause Analysis & Fix Strategy
