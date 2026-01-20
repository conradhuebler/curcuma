# GFN-FF CLI Tests Documentation

## Overview

The GFN-FF CLI tests validate the native GFN-FF implementation in Curcuma against reference values from XTB.

## Test Cases

### 01_complex_gfnff_singlepoint
Tests GFN-FF single point calculation on a large molecule (231 atoms) to validate:
- Energy accuracy against XTB reference
- Computational performance
- Gradient norm consistency

### 02_caffeine_gfnff_energy_components
Tests GFN-FF energy component breakdown and charges for caffeine (24 atoms):
- Validates all energy components (bond, angle, torsion, repulsion, electrostatic, dispersion, HB)
- Validates EEQ charges against XTB reference (24 charge values)
- Identifies discrepancies in bond energy calculation
- Expected to fail initially due to known ~0.6 Eh bond energy error
- Provides baseline measurement for debugging GFN-FF implementation

## Reference Values (XTB 6.6.1)

### Complex Molecule (complex.xyz)
- **Total Energy**: -37.241047037682 Eh
- **Gradient Norm**: 0.875832449077 Eh/a0
- **Calculation Time**: ~0.1 seconds

### Caffeine Molecule (caffeine.xyz)
- **Total Energy**: -4.672736999448 Eh
- **Bond Energy**: -4.798987726686 Eh
- **Angle Energy**: 0.018051428064 Eh
- **Torsion Energy**: 0.000898394346 Eh
- **Repulsion Energy**: 0.300917862858 Eh
- **Electrostatic Energy**: -0.174362982362 Eh
- **Dispersion Energy**: -0.018125118259 Eh
- **HB Energy**: -0.000000034030 Eh
- **Charges**: 24 reference values in reference_charges.txt

## Validation Criteria

### Energy Accuracy
- Tolerance: 1e-6 Eh (absolute difference)
- Validates energetic correctness of the implementation

### Performance
- Target: <1 second computation time for the test molecule
- Monitors scalability and efficiency

## Execution

To run all GFN-FF tests:
```bash
ctest -R "cli_gfnff_" --output-on-failure
```

To run a specific test:
```bash
ctest -R "cli_gfnff_01" --verbose
```

## Test Structure

Each test follows the standard CLI test pattern:
```
01_complex_gfnff_singlepoint/
├── run_test.sh        # Main test script
├── complex.xyz        # Test molecule (231 atoms)
├── stdout.log         # Generated output
└── stderr.log         # Generated errors

02_caffeine_gfnff_energy_components/
├── run_test.sh        # Main test script
├── caffeine.xyz       # Test molecule (24 atoms)
├── reference_charges.txt  # XTB reference charges (24 values)
├── stdout.log         # Generated output
└── stderr.log         # Generated errors
```

## Development Guidelines

1. Add new test cases for different molecule sizes and complexities
2. Update reference values when XTB or Curcuma implementations change
3. Monitor performance regressions through timing metrics
4. Document any deviations from expected behavior

## Related Files

- `test_cases/molecules/larger/complex.xyz` - Source molecule
- `gfnff_goal.md` - Implementation targets and requirements
- `test_cases/cli/test_utils.sh` - Shared test utilities