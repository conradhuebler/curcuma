# Parameter Registry Migration Guide

**For**: Curcuma developers migrating capabilities from static JSON to Parameter Registry
**Difficulty**: Easy (~30 minutes per module)
**Prerequisites**: Basic C++ knowledge, understanding of JSON

## Quick Start Checklist

- [ ] Identify all parameters in static JSON
- [ ] Add PARAM block to capability header
- [ ] Update constructor to use registry
- [ ] Remove static JSON configuration
- [ ] Build and verify
- [ ] Test functionality

## Step-by-Step Migration

### Example Module: Casino (Monte Carlo Simulation)

#### Step 1: Analyze Existing Code

**Current implementation** (`src/capabilities/casino.h` + `.cpp`):

```cpp
// casino.h
class Casino : public CurcumaMethod {
public:
    Casino(const json& controller, bool silent);
    // ...
};

// casino.cpp
static const json CasinoJson = {
    { "steps", 10000 },
    { "temperature", 300.0 },
    { "step_size", 0.1 },
    { "move_type", "random" },
    { "save_trajectory", false },
    { "output_interval", 100 }
};

Casino::Casino(const json& controller, bool silent)
    : CurcumaMethod(CasinoJson, controller, silent)
{
    m_config = MergeJson(CasinoJson, controller);
}
```

#### Step 2: Include Parameter Macros

Add to top of `casino.h`:

```cpp
#include "src/core/parameter_macros.h"
```

#### Step 3: Create PARAM Block

Add to private section of class in `casino.h`:

```cpp
class Casino : public CurcumaMethod {
public:
    Casino(const json& controller, bool silent);
    // ... rest of interface ...

private:
    // vvvvvvvvvvvv PARAMETER DEFINITION BLOCK vvvvvvvvvvvv
    BEGIN_PARAMETER_DEFINITION(casino)

    PARAM(steps, Int, 10000,
          "Number of Monte Carlo steps to perform",
          "Basic",
          {"nstep", "n_steps"})

    PARAM(temperature, Double, 300.0,
          "Simulation temperature in Kelvin",
          "Basic",
          {"temp", "T"})

    PARAM(step_size, Double, 0.1,
          "Maximum displacement per MC move in Angstrom",
          "Algorithm",
          {"delta", "displacement"})

    PARAM(move_type, String, "random",
          "Type of MC move: random|systematic|adaptive",
          "Algorithm",
          {})

    PARAM(save_trajectory, Bool, false,
          "Save trajectory to .trj.xyz file",
          "Output",
          {"save_trj"})

    PARAM(output_interval, Int, 100,
          "Write structure every N steps",
          "Output",
          {"write_interval"})

    END_PARAMETER_DEFINITION
    // ^^^^^^^^^^^^ PARAMETER DEFINITION BLOCK ^^^^^^^^^^^^

    json m_config;
    // ... other member variables ...
};
```

**IMPORTANT**:
- Use **snake_case** for parameter names (`step_size`, not `StepSize`)
- Match types exactly: `Int`, `Double`, `String`, `Bool`
- Provide clear, user-facing help text
- Add old parameter names as aliases for backward compatibility

#### Step 4: Update Constructor

Modify `casino.cpp`:

```cpp
// REMOVE static JSON:
// static const json CasinoJson = { ... };  // DELETE THIS

// UPDATE constructor:
#include "src/core/parameter_registry.h"  // Add this include

Casino::Casino(const json& controller, bool silent)
    : CurcumaMethod(ParameterRegistry::getInstance().getDefaultJson("casino"),
                    controller, silent)
{
    UpdateController(controller);
    m_config = MergeJson(ParameterRegistry::getInstance().getDefaultJson("casino"),
                        controller);

    // Rest of constructor unchanged...
}
```

#### Step 5: Build and Test

```bash
cd release

# Rebuild parser (if modified)
make curcuma_param_parser

# Generate new registry
make GenerateParams

# Should output:
# Successfully generated .../parameter_registry.h with X parameter definitions

# Check for errors
cat generated/parameter_registry.h | grep casino
# Verify all 6 parameters present

# Full rebuild
make curcuma

# Test for validation errors
./curcuma -help 2>&1 | grep -i "duplicate\|error"
# Should be clean (no output)

# Test module still works
./curcuma -casino input.xyz -steps 5000 -temperature 350
# Should run with custom parameters
```

#### Step 6: Verify Functionality

**Test default values**:
```bash
./curcuma -casino input.xyz
# Should use: steps=10000, temperature=300.0, etc.
```

**Test aliases**:
```bash
./curcuma -casino input.xyz -temp 350 -nstep 5000
# Aliases should resolve to canonical names
```

**Test help**:
```bash
./curcuma -help-module casino  # Once implemented
# Should show all 6 parameters with descriptions
```

## Common Patterns

### Pattern 1: Simple Value Parameters

**Before**:
```cpp
static const json MyJson = {
    { "max_iterations", 100 },
    { "threshold", 1e-6 }
};
```

**After**:
```cpp
BEGIN_PARAMETER_DEFINITION(my_module)

PARAM(max_iterations, Int, 100,
      "Maximum optimization iterations",
      "Algorithm",
      {"max_iter", "MaxIter"})

PARAM(threshold, Double, 1e-6,
      "Convergence threshold",
      "Algorithm",
      {"conv_thresh"})

END_PARAMETER_DEFINITION
```

### Pattern 2: Boolean Flags

**Before**:
```cpp
{ "verbose", true },
{ "save_output", false }
```

**After**:
```cpp
PARAM(verbose, Bool, true,
      "Enable detailed logging",
      "Output",
      {})

PARAM(save_output, Bool, false,
      "Write results to file",
      "Output",
      {"write_results"})
```

### Pattern 3: String Choices

**Before**:
```cpp
{ "method", "lbfgs" },
{ "format", "xyz" }
```

**After**:
```cpp
PARAM(method, String, "lbfgs",
      "Optimization method: lbfgs|cg|sd",
      "Algorithm",
      {})

PARAM(format, String, "xyz",
      "Output format: xyz|mol2|pdb",
      "Output",
      {"output_format"})
```

### Pattern 4: Nested JSON (Flatten)

**Before**:
```cpp
{
    "optimization", {
        {"max_iter", 100},
        {"threshold", 1e-6}
    }
}
```

**After** (flatten with prefix):
```cpp
PARAM(optimization_max_iter, Int, 100,
      "Maximum optimization iterations",
      "Optimization",
      {"opt_max_iter"})

PARAM(optimization_threshold, Double, 1e-6,
      "Optimization convergence threshold",
      "Optimization",
      {"opt_thresh"})
```

**Alternative** (CLI still supports dot notation):
```cpp
PARAM(max_iter, Int, 100, "...", "Optimization", {})
PARAM(threshold, Double, 1e-6, "...", "Optimization", {})

// CLI: -optimization.max_iter 200 (dot notation works automatically)
```

## Migration Checklist Per Module

### Pre-Migration
- [ ] Read existing static JSON
- [ ] Identify all parameter names, types, defaults
- [ ] Note any special handling or nested structures
- [ ] Check for parameter aliases in CLI parsing code

### During Migration
- [ ] Add `#include "src/core/parameter_macros.h"`
- [ ] Create PARAM block with all parameters
- [ ] Use snake_case for parameter names
- [ ] Provide comprehensive help text
- [ ] Add old names as aliases
- [ ] Update constructor to use registry
- [ ] Remove static JSON declaration

### Post-Migration
- [ ] Build succeeds without errors
- [ ] Generated registry contains all parameters
- [ ] No validation warnings (`./curcuma -help`)
- [ ] Module runs with default parameters
- [ ] Module runs with custom CLI arguments
- [ ] Aliases work correctly
- [ ] Help text displays properly (once implemented)

## Handling Special Cases

### Shared Parameters

If parameter appears in multiple modules (e.g., `temperature` in casino + simplemd):

**Best Practice**:
- Use **identical** type, default, and help text
- Parser auto-detects shared parameters
- Validation warns if types differ

```cpp
// casino.h
PARAM(temperature, Double, 300.0,
      "Simulation temperature in Kelvin",
      "Basic",
      {"temp", "T"})

// simplemd.h (IDENTICAL definition)
PARAM(temperature, Double, 300.0,
      "Simulation temperature in Kelvin",
      "Basic",
      {"temp", "T"})
```

### Complex Nested Configurations

For deeply nested JSON like `topological.*` in analysis:

**Approach 1**: Flatten with prefix
```cpp
PARAM(topological_save_image, Bool, false, "...", "Topology", {})
PARAM(topological_colormap, String, "hot", "...", "Topology", {})
```

**Approach 2**: Keep flat in PARAM, reconstruct in code
```cpp
// Constructor
json topological = {
    {"save_image", Json2KeyWord<bool>(m_config, "topological_save_image")},
    {"colormap", Json2KeyWord<std::string>(m_config, "topological_colormap")}
};
```

### Optional Parameters

For parameters that may or may not be set:

```cpp
PARAM(output_file, String, "",
      "Optional output file path (empty = no file output)",
      "Output",
      {"out"})

// In code:
std::string out_file = Json2KeyWord<std::string>(m_config, "output_file");
if (!out_file.empty()) {
    // Write to file
}
```

## Testing Strategy

### Unit Test Template

```bash
#!/bin/bash
# test_casino_migration.sh

MODULE="casino"

echo "Testing $MODULE migration..."

# 1. Build
echo "[1/5] Building..."
cd release
make curcuma_param_parser > /dev/null 2>&1
make GenerateParams > /dev/null 2>&1
make curcuma > /dev/null 2>&1 || {
    echo "❌ Build failed"
    exit 1
}
echo "✅ Build succeeded"

# 2. Check validation
echo "[2/5] Validation..."
./curcuma -help 2>&1 | grep -i "duplicate.*$MODULE" && {
    echo "❌ Duplicate parameters detected"
    exit 1
}
echo "✅ No validation errors"

# 3. Test defaults
echo "[3/5] Default parameters..."
./curcuma -$MODULE input.xyz > /dev/null 2>&1 || {
    echo "❌ Failed with defaults"
    exit 1
}
echo "✅ Defaults work"

# 4. Test custom params
echo "[4/5] Custom parameters..."
./curcuma -$MODULE input.xyz -steps 500 -temperature 350 > /dev/null 2>&1 || {
    echo "❌ Failed with custom params"
    exit 1
}
echo "✅ Custom params work"

# 5. Test aliases
echo "[5/5] Aliases..."
./curcuma -$MODULE input.xyz -nstep 500 -temp 350 > /dev/null 2>&1 || {
    echo "❌ Aliases don't work"
    exit 1
}
echo "✅ Aliases work"

echo "
🎉 All tests passed for $MODULE!"
```

## Troubleshooting

### Build Error: "error: 'ParameterRegistry' has not been declared"
**Fix**: Add `#include "src/core/parameter_registry.h"` to `.cpp` file

### Runtime: Parameters not found
**Fix**: Check `initialize_generated_registry()` is called in `main.cpp`

### Validation: "Duplicate parameter 'foo'"
**Fix**: Check for duplicate PARAM entries or alias conflicts in PARAM block

### Parser: "Malformed PARAM around line X"
**Fix**: Verify PARAM syntax: `PARAM(name, Type, value, "help", "category", {})`
- Check all commas present
- Verify balanced quotes
- Ensure parentheses closed

### Functionality broken after migration
**Debug**: Compare old JSON with registry output:
```bash
# Export current defaults
./curcuma -export-config casino

# Compare with old CasinoJson
# Values should match exactly
```

## Migration Priority

Recommended order (easiest → most complex):

1. ✅ **analysis** (Done - 25 params) - Reference implementation
2. **casino** (~14 params) - Simple, no nesting
3. **hessian** (~10 params) - Small, straightforward
4. **simplemd** (~20 params) - Has shared `temperature` parameter
5. **opt** (~15 params) - Critical module, moderate complexity
6. **confscan** (~12 params) - Moderate complexity
7. **rmsd** (~15 params) - Complex, many nested options

## Getting Help

- Technical details: [PARAMETER_SYSTEM.md](PARAMETER_SYSTEM.md)
- Reference implementation: `src/capabilities/analysis.h`
- Parser code: `scripts/param_parser/main.cpp`
- Registry implementation: `src/core/parameter_registry.{h,cpp}`

## Contributing

After migrating a module:
1. Test thoroughly (use checklist above)
2. Update migration status in this doc
3. Add example to `PARAMETER_SYSTEM.md` if pattern is novel
4. Commit with message: `✅ Migrate <module> to Parameter Registry`
