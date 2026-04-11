# GPU GFN-FF Gradient Discrepancies

**Date**: 2026-03-27
**Status**: Investigation in progress

## Summary

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

## Next Steps

1. **Per-term GPU gradient isolation** — Add diagnostic output to `k_hbonds` kernel separating damping/out-of-line/neighbor gradient contributions per atom
2. **Compare GPU vs CPU per-atom HB gradient accumulation** — Identify which specific gradient sub-term diverges on donor/acceptor atoms
3. **atomicAdd accumulation order** — Test if Kahan summation or deterministic reduction eliminates the error (4.2e-3 is too large for simple FP rounding)
4. **HB SoA data verification** — Confirm GPU H-bond pair/parameter arrays match CPU data