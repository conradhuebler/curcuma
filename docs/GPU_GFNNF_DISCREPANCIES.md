# GPU GFN-FF Gradient Discrepancies

**Date**: 2026-03-27 (original) / **2026-06-19 (RESOLVED)**
**Status**: ✅ RESOLVED — the GPU GFN-FF gradient is correct; the HB-charge freeze is correct.

## RESOLUTION (June 2026)

The long-open question "is the GPU HB-charge freeze a bug?" is settled: **no, the freeze is
correct.** The GPU GFN-FF gradient matches the CPU at every geometry; the earlier
GPU-analytical-vs-GPU-FD discrepancy (4.2e-3) and the CPU-vs-GPU MD heat-exchange differences
were **measurement artifacts**, not force bugs.

Evidence (this machine, GTX 1660, current code):
- **Static gradient, per-geometry**: a new per-frame diagnostic
  (`test_cases/test_gfnff_grad_traj.cpp`) evaluates the gradient of each backend at every
  frame of a fixed trajectory — the clean metric (no integrator, no chaos). On a 51-frame
  acetic-acid-dimer trajectory: **GPU gradient == CPU gradient, mean 8.2e-10, max 7e-9** per
  frame. The GPU force is bit-correct at every geometry.
- **The HB charges are FROZEN by design.** The CPU HB term uses `hb.q_H/q_A/q_B` stored at
  topology build (`forcefieldthread.cpp:2525+`), and so does the Fortran reference. So the
  GPU's frozen HB charges already match the reference. A device-resident gather that refreshes
  the HB charges to the live per-step EEQ charges was implemented and tested — it makes MD
  **diverge**, not converge:

  | acetic dimer, 1000 fs | "Exchange with heat bath" | Δ vs CPU |
  |---|---|---|
  | CPU (reference)       | 0.00709227 Eh | — |
  | GPU **frozen** (default) | 0.00708362 Eh | **8.7e-6** ✓ |
  | GPU resident (opt-in) | 0.00610579 Eh | 5.1e-4 ✗ |

  The April-2026 "GPU analytical vs GPU FD = 4.2e-3" was the charge-re-solving FD artifact
  flagged as a caveat in that section: the FD perturbation re-solved the EEQ charges while the
  analytical gradient (correctly) held the HB charges frozen, so the two disagreed — that is a
  property of the finite-difference reference, not a gradient error.

- **MD "Exchange with heat bath" is NOT a valid cross-backend force diagnostic.** It conflates
  force differences with integrator/thermostat/velocity-init differences and chaotic
  amplification, and the GPU's own run-to-run `atomicAdd` nondeterminism is ~5e-5 Eh on a
  floppy 1000 fs trajectory — larger than the CPU-vs-GPU-frozen difference. Use the per-frame
  gradient diagnostic instead.

**Code outcome**: the device-resident gather (`k_gather_hb_charges`,
`FFWorkspaceGPU::refreshHBChargesFromDevice()`) is kept **opt-in, default OFF** — enable with
`CURCUMA_GFNFF_GPU_RESIDENT_HBQ=1` for experimentation only; the correct (frozen) behavior is
the default. The per-frame harness (`test_gfnff_grad_traj`) is the recommended force metric.

The historical investigation below is retained for context.

---

## Summary (historical, 2026-03-27)

GPU GFN-FF gradient validation tests show discrepancies between CPU and GPU implementations. After fixing the dEdcn snapshot bug and dlogdcn formula, 16/19 tests pass but gradient errors remain for larger molecules.

## Test Results

### Energy Comparison (Complex Molecule, 231 atoms)

| Method | Single Point Energy | Diff from XTB |
|--------|---------------------|----------------|
| Native GFN-FF (GPU) | -37.24064797 Eh | +0.72 μEh |
| XTB GFN-FF (Fortran) | -37.24064869 Eh | reference |
| CPU GFN-FF | -37.24064797 Eh | +0.72 μEh |

Energy agreement is excellent (< 1 μEh).

### MD "Exchange with Heat Bath" Comparison

| Method | Exchange Value | Diff |
|--------|----------------|------|
| Native GFN-FF (GPU) | -7.70874 Eh | +0.204 Eh |
| XTB GFN-FF (Fortran) | -7.91275 Eh | reference |

This ~0.2 Eh difference in MD energy exchange suggests gradient accumulation errors during dynamics.

### Gradient Validation Results

| Test | dEdcn | Grad Max Diff | Grad Norm Ratio | Status |
|------|-------|----------------|-----------------|--------|
| H2O | PASS | <1e-4 | ~1.0 | ✅ PASS |
| HCN | PASS | <1e-4 | ~1.0 | ✅ PASS |
| HCl | PASS | <1e-4 | ~1.0 | ✅ PASS |
| HH | PASS | <1e-4 | ~1.0 | ✅ PASS |
| O3 | PASS | <1e-4 | ~1.0 | ✅ PASS |
| OH | PASS | <1e-4 | ~1.0 | ✅ PASS |
| CH3OH | PASS | <1e-4 | ~1.0 | ✅ PASS |
| CH4 | PASS | 6e-5 | 0.976 | ❌ FAIL (norm) |
| CH3OCH3 | PASS | 4e-4 | 1.017 | ❌ FAIL |
| C6H6 | PASS | 2e-4 | 1.051 | ❌ FAIL |
| acetic_acid_dimer | PASS | 3.7e-3 | 1.005 | ❌ FAIL |
| complex (231 atoms) | PASS | 7.6e-3 | 1.008 | ❌ FAIL |
| polymer (1410 atoms) | PASS | 1e-3 | 1.01 | ❌ FAIL |

**Pattern**: Gradient errors scale with molecule size, suggesting accumulation of small errors or missing terms.

## Fixes Applied

### 1. dEdcn Snapshot Bug (ff_workspace_gpu.cu)

**Issue**: `dEdcnTotal()` returned zeros because dEdcn snapshot was only populated at verbosity ≥ 3.

**Fix**: Always populate and download dEdcn snapshot regardless of verbosity.

```cpp
// Before (broken):
if (need_snapshots) {
    cudaMemcpyAsync(impl.d_dEdcn_snapshot.ptr, impl.d_dEdcn.ptr, ...);
}

// After (fixed):
cudaMemcpyAsync(impl.d_dEdcn_snapshot.ptr, impl.d_dEdcn.ptr, ...);
```

**Result**: dEdcn values now match between CPU and GPU (< 1e-10 error).

### 2. dlogdcn Formula (ggfnff_method.cpp)

**Issue**: GPU dlogdcn formula was mathematically incorrect, giving wrong values (even negative for cn > 4.4).

**Original (wrong)**:
```cpp
double expval = std::exp(lncnmax - m_gpu_cn_final[i]);
dlogdcn[i] = (expval - 1.0) / expval;
```

**Fixed (matches Fortran gfnff_cn.f90)**:
```cpp
dlogdcn[i] = exp_cnmax / (exp_cnmax + std::exp(m_gpu_cn_final[i]));
```

**Result**: dlogdcn values now match Fortran reference, but gradient errors persist.

### 3. File Path Fixes

Three reference JSON files had incorrect geometry paths:
- `HCN.ref.json`: Fixed `../external/gfnff/test/HCN.xyz` → `test_cases/molecules/trimers/HCN.xyz`
- `methane.ref.json`: Fixed `../external/gfnff/test/methane.xyz` → `test_cases/molecules/larger/CH4.xyz`
- `acetic_acid_dimer.ref.json`: Fixed `external/gfnff/test/acetic_acid_dimer.xyz` → `test_cases/molecules/larger/acetic_acid_dimer.xyz`

## Remaining Issues

### CN Chain-Rule Gradient Computation

The GPU kernel `k_cn_chainrule` implements the CN chain-rule gradient contribution. After analysis:

**CPU (gfnff_method.cpp)**:
```cpp
// Sparse matrix dcn built with dlogdcn-weighted derivatives
// g(i) = -(dlogdcn_i * dEdcn(i) + dlogdcn_j * dEdcn(j)) * comp
```

**GPU (gfnff_kernels.cu)**:
```cpp
double fac = dSdr / rij * (dEdcn[i] * dlogdcn[i] + dEdcn[j] * dlogdcn[j]);
add_grad(grad, i, fac * dx, fac * dy, fac * dz);
add_grad(grad, j, -fac * dx, -fac * dy, -fac * dz);
```

**Fortran Reference (gfnff_cn.f90 lines 112-115)**:
```fortran
dlogCN(:,j,j) = dlogdcnj*rij + dlogCN(:,j,j)
dlogCN(:,i,j) = -dlogdcnj*rij
dlogCN(:,j,i) = dlogdcni*rij
dlogCN(:,i,i) = -dlogdcni*rij + dlogCN(:,i,i)
```

The formulas appear mathematically equivalent, yet gradient discrepancies remain. Possible causes:

1. **Different CN pair enumeration** - CPU iterates j < i, GPU uses stored pair list
2. **Floating-point accumulation order** - GPU parallel reduction may differ from CPU sequential
3. **Other gradient terms** - Bonds, dispersion, Coulomb chain-rule terms may have subtle differences

### Investigation Needed

1. **CN pair list verification** - Compare CPU and GPU pair lists for identical molecules
2. **Gradient term isolation** - Test each gradient term independently
3. **Numerical precision** - Check if double precision accumulations differ
4. **Stream synchronization** - Verify all GPU streams complete before gradient download

## Files Modified

1. `src/core/energy_calculators/ff_methods/cuda/ff_workspace_gpu.cu` - dEdcn snapshot fix
2. `src/core/energy_calculators/qm_methods/ggfnff_method.cpp` - dlogdcn formula fix
3. `test_cases/reference_data/HCN.ref.json` - file path fix
4. `test_cases/reference_data/methane.ref.json` - file path fix
5. `test_cases/reference_data/acetic_acid_dimer.ref.json` - file path fix

## April 2026: Isolated HB Gradient Test Results

### Code Analysis: HB Formulas Identical

Systematic comparison of all HB gradient code (Cases 1-4: eg1, eg2new, eg2_rnr, eg3) across CPU (`forcefieldthread.cpp`), GPU (`gfnff_kernels.cu`), and Fortran (`gfnff_engrad.F90`) confirmed:

- **All formulas are mathematically identical** — damping, out-of-line, neighbor, angle/torsion terms
- Neither CPU, GPU, nor Fortran compute ∂E_HB/∂q·∂q/∂x (charges treated as constant)
- HB-CN bond modulation chain rule is correctly implemented on GPU (`k_bonds` + `k_hb_alpha_chainrule`)

### Isolated HB Gradient Test (`test_hb_gradient_isolation`)

New test comparing CPU analytical, GPU analytical, and finite-difference gradients (acetic acid dimer, all terms enabled):

| Comparison | Max Diff (Eh/Bohr) | Tolerance | Result |
|-----------|---------------------|-----------|--------|
| CPU anal vs CPU FD | 5.7e-9 | 1e-4 | ✅ PASS |
| GPU anal vs GPU FD | 4.2e-3 | 1e-4 | ❌ FAIL |
| CPU anal vs GPU anal | 3.7e-3 | 1e-6 | ❌ FAIL |
| CPU FD vs GPU FD | 2.0e-3 | 1e-6 | ❌ FAIL |

**Key finding**: CPU analytical gradient is perfect (5.7e-9 vs numerical). GPU has real gradient errors (4.2e-3).

### Error Pattern

- **Affected atoms**: Donor O (atom 2), acceptor O (atom 5), neighboring C atoms — all H-bond participants
- **Unaffected atoms**: Methyl H atoms, H-bonding H atoms — perfect GPU gradients
- **Implication**: Bug is in gradient distribution to A/B (donor/acceptor) atoms in GPU HB kernel, not in the formula itself

### Fixes Applied (April 2026)

1. **Term flag forwarding to GPU** (`ff_workspace_gpu.cu`, `gfnff_method.cpp`): GFNFFParameterSet flags now populated from m_parameters and read in FFWorkspaceGPU constructor
2. **GPU method flag forwarding** (`gfnff_gpu_method.cpp`): Added setDispersionEnabled/setHBondEnabled/etc. calls
3. **ConfigManager sub-object stripping**: Discovered that `config_manager.cpp:72` drops JSON sub-objects — flags must be at top-level in controller JSON

### Test Infrastructure

- `test_cases/test_hb_gradient_isolation.cpp` — Standalone test executable
- Added to `test_cases/CMakeLists.txt` with CTest registration
- Supports `--all-terms` flag and custom molecule path

## April 2026 Update: GPU Gradient May Be More Correct Than CPU

**Date**: 2026-04-29  
**Status**: Revised interpretation — investigation ongoing

### New Evidence

MD simulations (polymer system) and heat-exchange analysis comparing all three implementations:

| Method | Heat exchange (MD) | Agreement with XTB |
|--------|--------------------|--------------------|
| Native GFN-FF CPU | small deviations from XTB | worse |
| Native GFN-FF GPU (`-gpu cuda`) | closer to XTB | better |
| XTB GFN-FF (Fortran reference) | reference | — |

The GPU implementation produces MD trajectories that agree better with XTB GFN-FF than the CPU implementation, as measured by heat-bath exchange values. The CPU path runs without instability but shows systematic deviations that are still open for production use.

### Revised Interpretation

The CPU vs GPU gradient discrepancies documented above (e.g. 3.7e-3 Eh/Bohr on acetic acid dimer) were previously assumed to indicate GPU errors. The new MD evidence suggests the opposite may be true:

- **The CPU gradient may contain errors** not present in the GPU path — possibly in HB gradient distribution, CN chain-rule accumulation order, or atomicAdd vs sequential summation.
- The GPU's `atomicAdd`-based accumulation may accidentally reproduce the Fortran summation order more faithfully than the CPU sequential loops.
- The isolated HB test result (CPU anal vs FD: 5.7e-9 ✅, GPU anal vs FD: 4.2e-3 ❌) is not overturned, but needs to be re-evaluated: the FD reference on the CPU may itself contain the same CPU error.

### What This Means for the Roadmap

1. **Do not treat GPU gradient as the bug source** until a careful per-term comparison against the Fortran reference is done (not against the CPU).
2. **Priority**: Compare GPU and CPU gradients separately against XTB numerical gradients (not against each other).
3. **The GPU path is currently the more reliable production path** for MD simulations.

### Still Unknown

- Which specific term(s) cause the CPU gradient to diverge from XTB
- Whether the GPU advantage holds for all molecule classes (only polymer/complex tested)
- Whether the GPU `atomicAdd` ordering is deterministic across runs

## Next Steps

1. **Direct comparison against XTB numerical gradient** — Run `xtb --grad` on test molecules, compare both native CPU and GPU against this reference, not against each other.
2. **Per-term CPU gradient audit** — Isolate which gradient term diverges in the CPU path relative to XTB (HB most likely candidate based on April isolation tests).
3. **Confirm GPU stability on more systems** — Repeat MD heat-exchange test on acetic acid dimer and caffeine.
4. **HB SoA data verification** — Confirm GPU H-bond pair/parameter arrays match CPU data (still relevant to understand the source of divergence).