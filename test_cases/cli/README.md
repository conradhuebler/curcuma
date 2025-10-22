# CLI Test Suite - Curcuma Argument Validation

**Created:** October 2025
**Author:** Claude (based on testing_plan_*.md specifications)
**Copyright:** (C) 2025 Conrad Hübler

## Overview

Comprehensive Bash-based end-to-end test suite validating all CLI arguments for Curcuma's main capabilities. These tests ensure correct behavior for valid inputs and proper error handling for invalid inputs.

## Test Coverage

### 📊 Statistics
- **Total Tests:** 26 scenarios
- **Capabilities:** 4 (curcumaopt, rmsd, confscan, simplemd)
- **Test Types:** Success paths, error handling, parameter validation, alias compatibility

### 🎯 Test Scenarios

#### curcumaopt (6 tests)
1. ✅ Default UFF optimization
2. ✅ GFN2 single point calculation
3. ❌ Invalid method error handling
4. ✅ LBFGS parameter passing
5. ✅ Alias backward compatibility (SinglePoint)
6. ✅ Hydrogen-only optimization

#### rmsd (6 tests)
1. ✅ Standard RMSD calculation
2. ✅ No reordering mode
3. ❌ Invalid alignment method error
4. ✅ Template with element selection
5. ✅ Fragment-based RMSD
6. ✅ Alias compatibility (reorder)

#### confscan (7 tests)
1. ✅ Default conformer scan
2. ✅ Dynamic RMSD threshold
3. ❌ Invalid RMSD method error
4. ✅ sLX threshold logic
5. ✅ Hybrid RMSD with elements
6. ✅ Heavy atoms only mode
7. ✅ Restart functionality

#### simplemd (7 tests)
1. ✅ Short NVE simulation
2. ✅ NVT with Berendsen thermostat
3. ❌ Invalid thermostat error
4. ✅ RATTLE constraints
5. ✅ Spherical wall potential
6. ✅ Temperature alias (T)
7. ✅ Restart simulation

## Directory Structure

```
cli/
├── README.md                   # This file
├── test_utils.sh               # Shared helper functions
├── template_test.sh            # Template for new tests
├── CMakeLists.txt              # CTest integration
├── curcumaopt/                 # Optimization tests
│   ├── 01_default_uff_opt/
│   │   ├── run_test.sh         # Test script
│   │   └── input.xyz           # Test molecule
│   └── ...
├── rmsd/                       # RMSD tests
├── confscan/                   # Conformer scanning tests
└── simplemd/                   # Molecular dynamics tests
```

## Running Tests

### All CLI Tests
```bash
cd build
ctest -R cli_ --output-on-failure
```

### Specific Capability
```bash
ctest -R cli_curcumaopt --output-on-failure
ctest -R cli_rmsd --output-on-failure
ctest -R cli_confscan --output-on-failure
ctest -R cli_simplemd --output-on-failure
```

### Individual Test
```bash
ctest -R cli_curcumaopt_01 --verbose
```

### Direct Execution
```bash
cd test_cases/cli/curcumaopt/01_default_uff_opt
./run_test.sh
```

## Test Utilities (`test_utils.sh`)

### Available Assertions
- `assert_exit_code <actual> <expected> <msg>` - Verify exit code
- `assert_file_exists <path> <msg>` - Check file existence
- `assert_file_not_exists <path> <msg>` - Verify file absence
- `assert_string_in_file <pattern> <file> <msg>` - Grep validation
- `assert_string_not_in_file <pattern> <file> <msg>` - Negative grep
- `assert_numeric_match <expected> <actual> <tolerance> <msg>` - Float comparison

### Helper Functions
- `extract_energy_from_xyz <file>` - Parse energy from XYZ comment
- `test_header <name>` - Print formatted test header
- `print_test_summary` - Display pass/fail statistics
- `cleanup_test_artifacts` - Remove common output files

## Creating New Tests

1. **Copy Template:**
   ```bash
   cp cli/template_test.sh cli/capability/XX_test_name/run_test.sh
   ```

2. **Modify Test:**
   - Set `TEST_NAME`
   - Implement `run_test()` function
   - Implement `validate_results()` function
   - Make executable: `chmod +x run_test.sh`

3. **Add to CMake:**
   Edit `cli/CMakeLists.txt` and add:
   ```cmake
   add_test(
       NAME cli_capability_XX_test_name
       COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/capability/XX_test_name/run_test.sh
       WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/capability/XX_test_name
   )
   ```

4. **Test Locally:**
   ```bash
   ./run_test.sh  # Direct execution
   cd build && ctest -R cli_capability_XX --verbose  # Via CTest
   ```

## Design Principles

### ✅ Best Practices
- **Small molecules:** Fast execution (water, ethane, methanol)
- **Short runs:** MD simulations limited to 5-10 fs
- **Isolated tests:** Each test in own directory with own input files
- **Clear output:** Color-coded pass/fail messages
- **Error validation:** Explicit tests for expected failures

### ❌ Anti-Patterns
- Don't use large molecules (>100 atoms) unless necessary
- Don't run long simulations in tests
- Don't share input files between tests (use copies)
- Don't rely on specific energy values (tolerances needed)

## Integration with CI/CD

Tests are automatically run via:
```yaml
# .github/workflows/build.yml
- name: Run tests
  run: |
    cd build
    ctest --output-on-failure
```

All CLI tests must pass for CI to succeed.

## Troubleshooting

### Test Fails with "curcuma binary not found"
- Ensure you built curcuma: `cd build && make`
- Check that binary exists at `build/curcuma` or `release/curcuma`
- Tests automatically search both locations

### Test Times Out
- Default timeout: 30 seconds
- Increase in `cli/CMakeLists.txt`: `set_tests_properties(...  PROPERTIES TIMEOUT 60)`

### Numerical Mismatch
- Energy/RMSD values may vary slightly between systems
- Adjust tolerance in `assert_numeric_match()`
- Or make validation informational (non-critical)

## References

- **Test Plans:** `testing_plan_*.md` files in project root
- **Overview:** `testing_plan_overview.md`
- **Parameter System:** `docs/PARAMETER_SYSTEM.md`

## Future Work

- [ ] Add golden reference comparisons (binary diff)
- [ ] Unit tests for ConfigManager, Thermostat classes
- [ ] Performance benchmarking tests
- [ ] Cross-platform compatibility tests (Linux/Mac)
