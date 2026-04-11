# GFN-FF Dispersion Debugging Plan: Caffeine 75 µEh Error

**Status**: Investigation in progress (Feb 11, 2026)
**Target**: Reduce Caffeine dispersion error from **75 µEh** to **< 10 µEh**
**Current Accuracy**: CH4 4 nEh ✅ | Caffeine 75 µEh ❌ | C6H6 0 nEh ✅

## Problem Statement

Curcuma GFN-FF dispersion for Caffeine (24 atoms):
- **Curcuma**: -0.01804996 Eh
- **Fortran Ref**: -0.01812512 Eh
- **Error**: +75 µEh (+0.41%)
- **Status**: Unacceptable (75× over < 1 µEh target)

The error is **systematic overestimation** (energy too negative). This is critical because:
1. Affects geometry optimization convergence
2. Impacts intermolecular interactions in MD simulations
3. Represents 0.41% error in a key energy component for medium-sized systems

## Root Cause Candidates (Priority Order)

### Hypothesis 1: CN Calculation Deviation in Complex π-Systems 🔴 HIGHEST PRIORITY
**Evidence**: Caffeine has extended aromatic system (3 fused rings)
- CN is crucial for D4 weighting: `weight[ref] = exp(-wf * (CN - CN_ref)²)`
- Small CN deviations amplify in Gaussian weighting
- π-system detection (EEQSolver) may fail in fused ring geometry

**Test**: Compare per-atom CN values with Fortran reference
```bash
# Phase 3A: CN Validation
# Add verbosity >= 2 output in d4param_generator.cpp:precomputeGaussianWeights()
# Print for each atom: "CN[i]={} vs refcn[0]={}" for neutral state
```

**If True**: Fix in EEQSolver coordination number calculation
**If False**: Move to Hypothesis 2

### Hypothesis 2: Reference State Weighting Formula Mismatch 🟡 HIGH PRIORITY
**Evidence**: D4 uses CN-only (not CN+charge like full D4)
- Formula: `weight = exp(-4.0 * (CN_atom - CN_ref)²)` per `gfnff_gdisp0.f90:405`
- Current implementation matches this (line 999-1003 in d4param_generator.cpp)
- But reference CN values (refcn) may differ from Fortran

**Test**: Verify refcn array matches Fortran dftd3param.f90
```cpp
// Check: m_refcn values for each element
// Should match Fortran dftd3param.f90 line ~1150
for (int elem = 0; elem < 5; ++elem) {  // Check C, H, N, O
    for (int ref = 0; ref < m_refn[elem]; ++ref) {
        printf("refcn[%d][%d] = %.6f\n", elem, ref, m_refcn[elem][ref]);
    }
}
```

**If True**: Fix reference CN data import
**If False**: Move to Hypothesis 3

### Hypothesis 3: C6 Casimir-Polder Integration Precision 🟡 MEDIUM PRIORITY
**Evidence**: Casimir-Polder uses numerical integration (trapezoidal rule)
- Frequency grid has 23 points
- Each point multiply alphaiw values with corrections
- Correction factors (ascale, sscale, secaiw) may have precision issues

**Test**: Compare pre-computed C6 values with Fortran
```cpp
// Add diagnostic output in precomputeC6ReferenceMatrix()
// For first 10 element-pair combinations:
printf("C6[C-C ref0-ref0] = %.10e (should match Fortran)\n",
       m_c6_flat_cache[c6FlatIndex(5, 5, 0, 0)]);
```

**If True**: Increase integration precision or verify frequency grid
**If False**: Move to Hypothesis 4

### Hypothesis 4: Correction Factor Table Incompleteness 🟡 MEDIUM PRIORITY
**Evidence**: Hydrogen count (d4_refh) was previously 0 (fixed Feb 8)
- Correction formula: `alpha = ascale * (alphaiw - hcount * sscale * secaiw)`
- Only validated for C, H, N, O, halogens
- May be incomplete for aromatic system treatment

**Test**: Audit correction factors for Caffeine elements
```cpp
// In computeC6Reference(), log:
// For each atom pair: ascale_i, ascale_j, hcount_i, hcount_j
// For Caffeine (C, H, N, O), all should be ≤ 1.0 for ascale
```

**If True**: Update correction factor tables from Fortran
**If False**: Move to Hypothesis 5

### Hypothesis 5: Damping Formula Parameter Deviation 🟢 LOWER PRIORITY
**Evidence**: GFN-FF modified BJ damping (fixed Jan 25)
- Parameters: a1=0.58, a2=4.80 (hardcoded in forcefieldthread.cpp:1653)
- R0 formula: `r0² = (a1*sqrt(r4r2ij) + a2)²`
- If r4r2ij incorrect → R0 incorrect → damping incorrect

**Test**: Extract and compare R0 values for largest dispersion pairs
```cpp
// Add debug output for top 5 pairs by |energy|
if (index < 5) {
    printf("Pair %d-%d: r0_sq=%.6f, r4r2ij=%.6f, a1=0.58 a2=4.80\n",
           disp.i, disp.j, disp.r0_squared, disp.r4r2ij);
}
```

**If True**: Verify damping parameters or r4r2ij calculation
**If False**: Anomaly in pair-specific calculation

## Systematic Debugging Workflow

### Phase 1: Enable Diagnostics (Est. 30 min)

**Step 1.1**: Modify forcefieldthread.cpp to enable dispersion logging for all molecules
```cpp
// Line 1655: Change from
bool disp_diag = (m_gfnff_dispersions.size() > 0 && m_gfnff_dispersions.size() <= 50);

// To:
bool disp_diag = CurcumaLogger::get_verbosity() >= 3;
```

**Step 1.2**: Add CN diagnostic output in d4param_generator.cpp
```cpp
// In precomputeGaussianWeights(), after line 1001:
if (CurcumaLogger::get_verbosity() >= 3 && i < 5) {  // First 5 atoms
    CurcumaLogger::info(fmt::format("CN diag[{}]: cn={:.4f} vs refcn_neutral={:.4f} weight[0]={:.4f}",
        i, cni, cni_ref, weights[0]));
}
```

**Step 1.3**: Create test runner script
```bash
#!/bin/bash
export OMP_NUM_THREADS=1
./release/test_cases/test_gfnff_validation test_cases/reference_data/caffeine.ref.json 2>&1 | tee caffeine_debug.log
```

### Phase 2: Collect Data (Est. 15 min)

Run with verbosity=3:
```bash
cd release
export CURCUMA_VERBOSITY=3
../release/test_cases/test_gfnff_validation ../test_cases/reference_data/caffeine.ref.json 2>&1 | tee /tmp/caffeine_diag.log
grep "DISP_CSV\|CN diag" /tmp/caffeine_diag.log > /tmp/dispersion_data.csv
```

**Expected Output Format**:
```
DISP_CSV: idx,i,j,rij,C6,r4r2ij,r0_squared,zetac6,t6,t8,energy
DISP_CSV: 0,0,1,1.891234,1.234567e-05,0.123456,23.456789,0.950000,0.012345,0.000123,-0.000012345
CN diag[0]: cn=3.5234 vs refcn_neutral=2.8900 weight[0]=0.8234
```

### Phase 3: Analysis (Est. 30 min)

**Step 3.1**: Per-atom CN analysis
```python
#!/usr/bin/env python3
import re
cn_data = {}
with open('/tmp/dispersion_data.csv') as f:
    for line in f:
        if 'CN diag' in line:
            match = re.search(r'CN diag\[(\d+)\]: cn=([\d.]+) vs refcn_neutral=([\d.]+)', line)
            if match:
                atom, cn, refcn = int(match.group(1)), float(match.group(2)), float(match.group(3))
                cn_data[atom] = (cn, refcn, abs(cn - refcn))

# Check if any atom has CN deviation > 0.01
errors = [(atom, cn, refcn, dev) for atom, (cn, refcn, dev) in cn_data.items() if dev > 0.01]
if errors:
    print("❌ HYPOTHESIS 1 CONFIRMED: CN deviation detected")
    for atom, cn, refcn, dev in sorted(errors, key=lambda x: x[3], reverse=True):
        print(f"  Atom {atom}: CN={cn:.4f} vs ref={refcn:.4f} (deviation {dev:.4f})")
else:
    print("✅ HYPOTHESIS 1 REJECTED: CN values match (deviation < 0.01)")
```

**Step 3.2**: Per-pair C6 analysis
```python
# Parse DISP_CSV data
import csv
pairs = []
with open('/tmp/dispersion_data.csv') as f:
    for line in f:
        if 'DISP_CSV:' in line:
            parts = line.replace('DISP_CSV:', '').split(',')
            if len(parts) >= 11:
                pairs.append({
                    'idx': int(parts[0]),
                    'i': int(parts[1]),
                    'j': int(parts[2]),
                    'rij': float(parts[3]),
                    'c6': float(parts[4]),
                    'r4r2ij': float(parts[5]),
                    'r0_sq': float(parts[6]),
                    'zetac6': float(parts[7]),
                    'energy': float(parts[10])
                })

# Sum all energies
total_disp = sum(p['energy'] for p in pairs)
print(f"Total dispersion from Curcuma: {total_disp:.10f} Eh")
print(f"Expected from Fortran:         -0.0181251200 Eh")
print(f"Error: {total_disp - (-0.0181251200):.10e} Eh ({(total_disp / (-0.0181251200) - 1) * 100:.2f}%)")

# Find largest contributors
top_pairs = sorted(pairs, key=lambda p: abs(p['energy']), reverse=True)[:10]
print("\nTop 10 dispersion pairs by magnitude:")
for p in top_pairs:
    print(f"  Pair {p['i']}-{p['j']}: C6={p['c6']:.4e}, E={p['energy']:.10e} Eh")
```

### Phase 4: Root Cause Identification (Est. 20 min)

1. **If CN deviations found**: Problem is EEQSolver or coordination number definition
   - Check EEQSolver::calculateCoordinationNumbers()
   - Compare exponential formula with Fortran: `CN = Σ 1/(1+exp(...))`

2. **If CN matches but C6 deviates**: Problem is Casimir-Polder integration
   - Verify alphaiw data matches Fortran (269 states)
   - Check frequency grid (23 points)
   - Verify correction factors (ascale, sscale, secaiw)

3. **If both CN and C6 match but total energy wrong**: Problem is pair summing or damping
   - Check BJ damping parameters (a1=0.58, a2=4.80)
   - Verify r4r2ij calculation (3*sqrt(z4r2_i)*sqrt(z4r2_j))
   - Check zetac6 calculation (charge weighting - should be ~1.0)

### Phase 5: Implement Fix (Est. 60 min)

Once root cause identified, implement fix in:
- **EEQSolver** if CN issue
- **d4param_generator.cpp** if C6/weighting issue
- **forcefieldthread.cpp** if damping/pair summing issue

Then re-run test and verify:
```bash
./release/test_cases/test_gfnff_validation test_cases/reference_data/caffeine.ref.json | grep "Dispersion"
# Expected: "Dispersion" line should show error < 10 µEh
```

## Files to Monitor

**Critical Files**:
- `src/core/energy_calculators/ff_methods/d4param_generator.cpp` (CN weighting, C6 calculation)
- `src/core/energy_calculators/ff_methods/forcefieldthread.cpp:1640-1750` (BJ damping, pair sum)
- `src/core/energy_calculators/ff_methods/eeq_solver.cpp` (CN calculation)

**Reference Data**:
- `src/core/energy_calculators/ff_methods/d4_reference_cn_fortran.cpp` (reference CN values)
- `src/core/energy_calculators/ff_methods/d4_reference_data_fixed.cpp` (reference charges)
- `src/core/energy_calculators/ff_methods/d4_alphaiw_data.cpp` (polarizabilities, 269 states)

**Test Data**:
- `test_cases/reference_data/caffeine.ref.json` (reference dispersion: -0.01812512 Eh)
- `test_cases/molecules/larger/caffeine.xyz` (test geometry, 24 atoms)

## Expected Completion

**Ideal Case** (CN hypothesis true): 2-4 hours (diagnosis + fix)
**Typical Case** (deeper investigation needed): 4-8 hours
**Complex Case** (architectural issue): 8-16 hours

## Success Criteria

✅ **PASS**: Caffeine dispersion error < 10 µEh (< 0.1% error)
✅ **PASS**: No regression on CH4 (< 10 nEh) or C6H6 (< 10 nEh)
✅ **PASS**: All 11+ test molecules < 100 µEh total error

---

*This plan follows the structure from GFNFF_STATUS.md Phase 3: Root Cause Analysis & Fix Strategy*
