# Next Session: Dispersion Test Reference Values Integration

**Status**: Test suite fully compiled and operational ✅

## What's Done

1. ✅ Created `test_cases/test_dispersion.cpp` (860 lines, 36 tests)
2. ✅ Updated `test_cases/CMakeLists.txt` with test configuration
3. ✅ Created `test_cases/DISPERSION_REFERENCE_GUIDE.md` (400 lines)
4. ✅ Fixed D3 interface linker error (implemented `UpdateGeometry()`)
5. ✅ Compiled test_dispersion executable successfully

## What Needs To Be Done (Next Session)

### Phase 1: Generate Reference Values with XTB

**Command template:**
```bash
xtb test_cases/validation/{molecule}.xyz --gfn 2 --d3 --func {functional}
# Extract: "dispersion energy: X.XXXXXXXXX Eh"

xtb test_cases/validation/{molecule}.xyz --gfn 2 --d4 --func {functional}
# Extract: "dispersion energy: X.XXXXXXXXX Eh"
```

**Test Matrix Summary:**
- **D3**: 18 tests (6 molecules × 3 scenarios: pbe0/bj, b3lyp/bj, pbe0/bj+atm)
- **D4**: 18 tests (6 molecules × 3 scenarios: pbe0, b3lyp, pbe0 no-atm)
- **Molecules**: HH.xyz, benzene.xyz, ethene.xyz, butane.xyz, HCl.xyz, OH.xyz

### Phase 2: Update test_dispersion.cpp with Reference Values

**Location of test cases:**
- D3 references: Lines ~180-280 in setupReferences()
- D4 references: Lines ~282-400 in setupReferences()

**Replace all `TODO_REFERENCE` (value 0.0) with actual values:**

Example:
```cpp
// OLD:
    TODO_REFERENCE,  // TODO: Fill from XTB 6.6.1

// NEW:
    -0.012345678,    // XTB 6.6.1: xtb benzene.xyz --gfn 2 --d3 --func pbe0
```

### Phase 3: Run Tests and Validate

```bash
cd release
cmake .. -DUSE_D3=ON -DUSE_D4=OFF  # Or -DUSE_D4=ON if available
make test_dispersion

# Run tests:
./test_dispersion --verbose
# Output: dispersion_test_results.csv

# CTest integration:
ctest -R "test_dispersion" --output-on-failure
```

## Test Reference Data Structure

### D3 Tests (Organized by molecule/scenario)
```
D3-1 to D3-6:   6 molecules × PBE0/BJ (no ATM)
D3-7 to D3-10:  4 molecules × B3LYP/BJ (no ATM)
D3-11 to D3-16: 6 molecules × PBE0/BJ + ATM
D3-17 to D3-18: 2 molecules × PBE0/Zero damping
```

### D4 Tests
```
D4-1 to D4-10:  6 molecules × (PBE0, B3LYP) with 3-body
D4-11 to D4-16: 6 molecules × PBE0 without 3-body
D4-17 to D4-18: 2 molecules × Custom parameters (s8=1.5, s6=0.8)
```

## Expected Value Ranges (for sanity checking)

| Molecule | Size | D3 Range | D4 Range |
|----------|------|----------|----------|
| H2 | 2 | ~-5e-6 to -1e-5 Eh | ~-1e-5 to -2e-5 Eh |
| HCl, OH | 2 | ~-1e-4 to -3e-4 Eh | ~-2e-4 to -5e-4 Eh |
| Ethene | 6 | ~-3e-3 to -5e-3 Eh | ~-4e-3 to -7e-3 Eh |
| Benzene | 12 | ~-1e-2 to -2e-2 Eh | ~-1.5e-2 to -2.5e-2 Eh |
| Butane | 14 | ~-5e-3 to -1e-2 Eh | ~-7e-3 to -1.2e-2 Eh |

**Note**: D4 typically ~20-30% more negative than D3

## Important Files

1. **Test suite**: `/home/conrad/src/claude_curcuma/curcuma/test_cases/test_dispersion.cpp`
   - Lines 180-400: setupReferences() with TODO markers
   - D3 tests: Lines ~180-270
   - D4 tests: Lines ~272-400

2. **Executable**: `/home/conrad/src/claude_curcuma/curcuma/release/test_dispersion`

3. **Reference guide**: `/home/conrad/src/claude_curcuma/curcuma/test_cases/DISPERSION_REFERENCE_GUIDE.md`
   - XTB/ORCA command examples
   - Expected value ranges
   - Verification methods

## Recent Fixes Applied

### D3 Interface UpdateGeometry Implementation
- **File**: `src/core/energy_calculators/qm_methods/dftd3interface.{h,cpp}`
- **Changes**:
  - Added `int m_natoms` member to header (line 94)
  - Set `m_natoms = atomtypes.size()` in `InitialiseMolecule()` (line 159)
  - Implemented `UpdateGeometry(const double* coord)` (lines 183-189)
  - Added `#include <cstring>` for memcpy

## Build Command

```bash
cd release
cmake .. -DUSE_D3=ON -DUSE_D4=OFF  # Configure with D3
make test_dispersion               # Build executable
./test_dispersion --verbose        # Run tests
```

## Expected Output

```
DFT-D3/D4 Dispersion Energy Test Suite
================================================================================
Total test cases: 36

[ OK ] D3_HH_pbe0_bj      Energy: -0.000012345 Eh (ref: -0.000012345) Err: 0.000000e+00
[ OK ] D3_Benzene_pbe0_bj Energy: -0.012345678 Eh (ref: -0.012345678) Err: 0.000000e+00
...
================================================================================
TEST SUMMARY
================================================================================
Total tests:   36
Passed:        36
Failed:        0
Skipped:       0 (method not compiled)
Success rate:  100.0%

CSV report saved to: dispersion_test_results.csv
```

## Timeline Estimate

- **Phase 1** (Generate XTB values): ~30-60 min (parallel for multiple molecules)
- **Phase 2** (Update test_dispersion.cpp): ~30 min (manual entry of 36 values)
- **Phase 3** (Validate and test): ~15 min (run tests, verify pass/fail)
- **Total**: ~2 hours

## Quick Checklist

- [ ] Run XTB for all 6 molecules with D3/D4 parameters
- [ ] Extract dispersion energies from XTB output
- [ ] Update all 36 TODO_REFERENCE values in test_dispersion.cpp
- [ ] Rebuild: `make test_dispersion`
- [ ] Run tests: `./test_dispersion --verbose`
- [ ] Verify: All 36 tests pass (0 failures, 0 skipped)
- [ ] Check CSV output: `dispersion_test_results.csv`
- [ ] Commit changes to git with proper message

## Notes

- XTB 6.6.1 is the reference implementation (document if using different version)
- Coordinates in test_dispersion.cpp are in Angstrom; conversion to Bohr handled automatically
- D4 tests may skip if `-DUSE_D4=OFF` was used in cmake (that's OK)
- Test tolerances: 1e-6 Eh (D3) and 1e-5 Eh (D4) - adjust if needed
- CSV report useful for tracking which tests pass/fail and error magnitudes

---

**Created**: 2025-12-15
**Status**: Ready for reference value integration in next session
