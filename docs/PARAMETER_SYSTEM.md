# Curcuma Parameter Registry System

**Status**: Production-ready (Phase 1 & 2 completed January 2025)
**Author**: Claude Code (supervised by Conrad Hübler)

## Overview

The Parameter Registry System provides automated, type-safe parameter management for Curcuma capabilities. Parameters are defined directly in capability headers using macros, extracted at build time, and compiled into a unified registry with auto-generated help, validation, and JSON export/import.

### Key Benefits

- **Single Source of Truth**: Parameters defined once in source code
- **Auto-Generated Help**: No manual help text synchronization
- **Type Safety**: Compile-time type checking for all parameters
- **Validation**: Automatic duplicate detection and type consistency checks
- **Alias Support**: Multiple names for the same parameter
- **JSON I/O**: Export defaults, import custom configurations

## Architecture

```
┌─────────────────────────────────────────────────────────────┐
│ 1. DEFINITION (Developer writes)                            │
│    src/capabilities/analysis.h:                             │
│                                                             │
│    BEGIN_PARAMETER_DEFINITION(analysis)                    │
│    PARAM(output_format, String, "human", "...", "Output", {})│
│    END_PARAMETER_DEFINITION                                │
└─────────────────┬───────────────────────────────────────────┘
                  │
                  ▼
┌─────────────────────────────────────────────────────────────┐
│ 2. BUILD TIME (CMake invokes)                              │
│    curcuma_param_parser scans all *.h files                │
│    Extracts PARAM(...) definitions                         │
│    Generates: release/generated/parameter_registry.h       │
└─────────────────┬───────────────────────────────────────────┘
                  │
                  ▼
┌─────────────────────────────────────────────────────────────┐
│ 3. RUNTIME (main.cpp initializes)                          │
│    initialize_generated_registry()                         │
│    validateRegistry()  // Checks duplicates, types         │
│    ParameterRegistry::getInstance()  // Singleton access   │
└─────────────────┬───────────────────────────────────────────┘
                  │
                  ▼
┌─────────────────────────────────────────────────────────────┐
│ 4. USAGE (Capabilities query)                              │
│    - getDefaultJson("analysis") → JSON with all defaults   │
│    - printHelp("analysis") → Formatted help output         │
│    - resolveAlias("analysis", "format") → "output_format"  │
│    - findDefinition("analysis", "output_format") → metadata│
└─────────────────────────────────────────────────────────────┘
```

## Developer Guide

### Adding Parameters to a New Capability

**Step 1**: Include the macro header in your capability header:

```cpp
// src/capabilities/my_capability.h
#include "src/core/parameter_macros.h"
#include "curcumamethod.h"

class MyCapability : public CurcumaMethod {
public:
    MyCapability(const json& controller, bool silent);
    void start() override;

private:
    // vvvvvvvvvvvv PARAMETER DEFINITION BLOCK vvvvvvvvvvvv
    BEGIN_PARAMETER_DEFINITION(my_capability)

    PARAM(max_iterations, Int, 100,
          "Maximum number of optimization iterations",
          "Algorithm",
          {"max_iter", "MaxIter"})

    PARAM(convergence_threshold, Double, 1e-6,
          "Energy convergence threshold in Hartree",
          "Algorithm",
          {"conv_thresh"})

    PARAM(output_file, String, "",
          "Optional output file path",
          "Output",
          {"out", "output"})

    PARAM(verbose, Bool, true,
          "Enable detailed logging",
          "Output",
          {})

    END_PARAMETER_DEFINITION
    // ^^^^^^^^^^^^ PARAMETER DEFINITION BLOCK ^^^^^^^^^^^^

    json m_config;
};
```

**Step 2**: Update constructor to use registry instead of static JSON:

```cpp
// BEFORE (old way):
static const json MyCapabilityJson = {
    {"max_iterations", 100},
    {"convergence_threshold", 1e-6},
    // ...
};

MyCapability::MyCapability(const json& controller, bool silent)
    : CurcumaMethod(MyCapabilityJson, controller, silent)
{
    m_config = MergeJson(MyCapabilityJson, controller);
}

// AFTER (new way):
MyCapability::MyCapability(const json& controller, bool silent)
    : CurcumaMethod(ParameterRegistry::getInstance().getDefaultJson("my_capability"),
                    controller, silent)
{
    UpdateController(controller);
    m_config = MergeJson(ParameterRegistry::getInstance().getDefaultJson("my_capability"),
                        controller);
}
```

**Step 3**: Build and verify:

```bash
cd release
make curcuma_param_parser  # Rebuild parser if needed
make GenerateParams        # Extract parameters → generated/parameter_registry.h
make curcuma              # Compile with new registry

# Verify no validation errors
./curcuma -help | grep -i "duplicate\|error"
```

### PARAM Macro Syntax

```cpp
PARAM(name, type, default_value, help_text, category, aliases)
```

**Arguments**:

| Argument | Type | Description | Example |
|----------|------|-------------|---------|
| `name` | identifier | Canonical parameter name (snake_case) | `max_iterations` |
| `type` | ParamType | `String`, `Int`, `Double`, `Bool` | `Int` |
| `default_value` | value | Default value matching type | `100` |
| `help_text` | string | User-facing description | `"Maximum iterations"` |
| `category` | string | Grouping for help display | `"Algorithm"` |
| `aliases` | list | Alternative names: `{}` or `{"alias1", ...}` | `{"max_iter"}` |

**Type Requirements**:
- **String**: Use quoted strings: `"default_value"`
- **Int**: Plain integers: `42`, `-10`, `0`
- **Double**: Floating point: `1.5`, `1e-6`, `0.0`
- **Bool**: `true` or `false` (lowercase)

### Parameter Categories

Standard categories for consistent help organization:

- **Basic**: Essential, frequently-used parameters
- **Algorithm**: Algorithm-specific settings
- **Output**: Output format, file paths, verbosity
- **Advanced**: Expert-level options
- **Topology**: Topological analysis (analysis module)
- **Trajectory**: Trajectory-specific (analysis module)

## Build System Integration

### CMake Dependency Chain

```cmake
# 1. Build the parser tool
add_executable(curcuma_param_parser scripts/param_parser/main.cpp)

# 2. Collect all project headers
file(GLOB_RECURSE ALL_PROJECT_HEADERS "src/*.h" "src/**/*.h")

# 3. Generate registry header (depends on parser + all headers)
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/generated/parameter_registry.h
    COMMAND curcuma_param_parser
            --output ${CMAKE_CURRENT_BINARY_DIR}/generated/parameter_registry.h
            --inputs ${ALL_PROJECT_HEADERS}
    DEPENDS curcuma_param_parser ${ALL_PROJECT_HEADERS}
    COMMENT "Generating C++ parameter registry from source headers"
)

# 4. Create target for generation step
add_custom_target(GenerateParams
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/generated/parameter_registry.h)

# 5. Core library depends on generated registry
add_dependencies(curcuma_core GenerateParams)

# 6. Include generated directory
include_directories(${CMAKE_CURRENT_BINARY_DIR}/generated)
```

**Incremental Builds**: Changing a PARAM block triggers:
1. Regeneration of `parameter_registry.h`
2. Recompilation of `parameter_registry.cpp`
3. Recompilation of `main.cpp` (includes generated header)
4. Relinking of `curcuma` executable

## Runtime API

### ParameterRegistry Singleton

```cpp
auto& registry = ParameterRegistry::getInstance();
```

### Core Methods

#### Get Default JSON Configuration

```cpp
json defaults = registry.getDefaultJson("analysis");
// Returns:
// {
//   "output_format": "human",
//   "max_iterations": 100,
//   ...
// }
```

#### Print Auto-Generated Help

```cpp
registry.printHelp("analysis");
// Outputs formatted help grouped by category:
// Parameters for module: analysis
// ============================================================
//
// [Output]
//   -output_format <string> (default: human)
//       Output format: human|json|csv
//       Aliases: format
//
// [Algorithm]
//   -max_iterations <int> (default: 100)
//       Maximum number of iterations
//       Aliases: max_iter, MaxIter
```

#### Validate Registry

```cpp
bool valid = registry.validateRegistry();
// Checks:
// - Duplicate parameter names within modules
// - Alias conflicts
// - Type consistency for shared parameters across modules
// Outputs errors to stderr if invalid
```

#### Resolve Alias

```cpp
std::string canonical = registry.resolveAlias("analysis", "format");
// Returns: "output_format"
```

#### Find Parameter Definition

```cpp
const ParameterDefinition* def = registry.findDefinition("analysis", "output_format");
if (def) {
    std::cout << "Type: " << (int)def->type << std::endl;
    std::cout << "Default: " << std::any_cast<std::string>(def->defaultValue) << std::endl;
}
```

## JSON Export/Import

### Export Default Configuration

```bash
./curcuma -export-config analysis > analysis_config.json
```

Generates:
```json
{
  "fragments": true,
  "metrics": "gyration,rout,end2end",
  "output_file": "",
  "output_format": "human",
  "properties": "all",
  "ripser": false,
  "statistics": "none",
  "topological_save_image": false,
  "verbose": true,
  "window": 10
}
```

### Import Custom Configuration

```bash
# Edit analysis_config.json (modify parameters)
./curcuma -analysis input.xyz -import-config analysis_config.json
```

### List All Modules

```bash
./curcuma -list-modules
# Output:
# Available modules:
#   analysis (25 parameters)
#   casino (14 parameters)  # After migration
#   ...
```

## Advanced Features

### Shared Parameters

When the same parameter name appears in multiple modules:

```cpp
// casino.h
PARAM(temperature, Double, 300.0, "Simulation temperature in K", "Basic", {"temp"})

// simplemd.h
PARAM(temperature, Double, 300.0, "Simulation temperature in K", "Basic", {"temp"})
```

**Automatic Detection**: Parser recognizes shared parameter
**Validation**: Warns if types differ between modules
**Benefit**: Consistent meaning across capabilities

### Multi-line PARAM Definitions

Parser supports spreading parameters across multiple lines:

```cpp
PARAM(topological_save_persistence_image,
      Bool,
      false,
      "Save persistence images (.PI files) for topological analysis",
      "Topology",
      {})
```

### Nested Parameters (Dot Notation in CLI)

While PARAM definitions are flat, CLI supports nested notation:

```bash
./curcuma -analysis file.xyz -topological.save_image true -topological.colormap viridis
```

Internally stored as:
```json
{
  "topological": {
    "save_image": true,
    "colormap": "viridis"
  }
}
```

## Migration Status

| Module | Status | Parameter Count | Notes |
|--------|--------|----------------|-------|
| analysis | ✅ Complete | 25 | Proof-of-concept, tested |
| casino | ⏳ Pending | ~14 | Monte Carlo simulation |
| simplemd | ⏳ Pending | ~20 | Molecular dynamics |
| opt | ⏳ Pending | ~15 | Geometry optimization |
| confscan | ⏳ Pending | ~12 | Conformational scanning |
| hessian | ⏳ Pending | ~10 | Second derivatives |
| rmsd | ⏳ Pending | ~15 | Structure comparison |

## Troubleshooting

### Build Errors

**Error**: `generated/parameter_registry.h: No such file or directory`
**Solution**: Run `make GenerateParams` before main build

**Error**: `Duplicate parameter 'foo' in module 'bar'`
**Solution**: Check PARAM block for duplicate names or alias conflicts

### Parser Warnings

**Warning**: `Malformed PARAM in file.h around line X`
**Cause**: Syntax error in PARAM macro (missing comma, unbalanced quotes)
**Fix**: Check PARAM syntax matches: `PARAM(name, Type, value, "help", "cat", {})`

### Runtime Issues

**Issue**: Parameter not found at runtime
**Check**: Was `initialize_generated_registry()` called in main.cpp?
**Verify**: Registry contains module: `./curcuma -list-modules`

## Best Practices

1. **snake_case Names**: Use `max_iterations`, not `MaxIterations` or `maxIterations`
2. **Comprehensive Help**: Write clear, user-facing descriptions
3. **Meaningful Categories**: Group related parameters logically
4. **Preserve Aliases**: Keep old parameter names as aliases during migration
5. **Test After Changes**: Always verify with `./curcuma -help | grep warning`

## Future Enhancements

- [ ] Auto-generate documentation website from registry
- [ ] Parameter dependency checking (if X set, Y required)
- [ ] Range validation (min/max values)
- [ ] Enum types for string parameters
- [ ] Parameter deprecation warnings

## References

- Implementation: `src/core/parameter_registry.{h,cpp}`
- Parser: `scripts/param_parser/main.cpp`
- Macros: `src/core/parameter_macros.h`
- Example: `src/capabilities/analysis.h`
- Migration Guide: [PARAMETER_MIGRATION_GUIDE.md](PARAMETER_MIGRATION_GUIDE.md)
