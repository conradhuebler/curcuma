# Parameter Flow Architecture - Complete Guide

## Overview

This document describes the complete parameter flow from CLI arguments through to capability implementation. Understanding this architecture is essential for:
- Adding new capabilities with proper parameter handling
- Fixing parameter-related bugs
- Implementing multi-module capabilities (like ConfScan + RMSD)
- Debugging parameter routing issues

## Parameter Flow Diagram

```
┌─────────────────────────────────────────────────────────────────┐
│ User Input: ./curcuma -confscan conf.xyz -rmsd.method subspace  │
└────────────────────────────┬────────────────────────────────────┘
                             │
                             ▼
┌─────────────────────────────────────────────────────────────────┐
│ LAYER 1: CLI2Json (src/main.cpp lines 184-342)                  │
│ Purpose: Convert command-line arguments to nested JSON structure │
│ Input:  argv array with -flag value pairs                        │
│ Output: controller dict with module-specific subdocuments        │
└────────────────────────────┬────────────────────────────────────┘
                             │
   ┌─────────────────────────┼─────────────────────────────┐
   │                         │                             │
   ▼                         ▼                             ▼
┌─────────────────┐   ┌──────────────────┐   ┌──────────────────┐
│ Parse Arguments │   │ setNestedJsonValue│   │ Extract Modules  │
│ -flag value     │   │ Create nested     │   │ "rmsd.method" →  │
│                 │   │ structures        │   │ contr["rmsd"][]  │
└────────┬────────┘   └────────┬─────────┘   └────────┬─────────┘
         │                     │                      │
         └─────────────────────┼──────────────────────┘
                               │
                               ▼
                    ┌────────────────────────┐
                    │ controller dict ready  │
                    │ {                      │
                    │   "confscan": {...},   │
                    │   "rmsd": {...},       │
                    │   "verbosity": 3       │
                    │ }                      │
                    └────────────┬───────────┘
                                 │
┌────────────────────────────────┴──────────────────────────────┐
│ LAYER 2: execute* Function (src/main.cpp)                     │
│ Purpose: Route to appropriate capability                       │
│ Examples: executeConfScan, executeRMSD, executeOptimization   │
│ 💡 Key: Just pass controller to capability - NO merging needed!│
└────────────────────────────────┬──────────────────────────────┘
                                 │
                                 ▼
┌────────────────────────────────────────────────────────────┐
│ LAYER 3: Capability Constructor (src/capabilities/*.cpp)   │
│ Example: ConfScan::ConfScan(const json& controller)        │
│                                                            │
│ 1. Calls CurcumaMethod() with:                           │
│    - defaults from ParameterRegistry                      │
│    - controller (user parameters)                         │
│    - silent flag                                          │
│                                                            │
│ 2. Initializes ConfigManager as Multi-Module:            │
│    m_config({modules...}, controller)                    │
│                                                            │
│ 3. ConfigManager automatically:                           │
│    - Loads defaults for each module                       │
│    - Merges user parameters from controller              │
│    - Handles alias resolution                            │
│    - Provides type-safe parameter access                 │
└────────────┬───────────────────────────────────────────────┘
             │
             ▼
┌──────────────────────────────────────────────────┐
│ LAYER 4: ConfigManager (src/core/config_manager)│
│                                                 │
│ Multi-Module Initialization:                   │
│ ConfigManager m_config({"confscan", "rmsd"},   │
│                        controller);             │
│                                                 │
│ For each module:                               │
│ 1. Load defaults: ParameterRegistry::getDefault│
│ 2. Extract user: controller[module_name]       │
│ 3. Merge: defaults + user parameters           │
│ 4. Resolve aliases & validate types            │
│ 5. Store in m_module_configs[module]           │
│                                                 │
│ Type-safe access:                              │
│ m_config.get<std::string>("rmsd.method")      │
└──────────┬───────────────────────────────────────┘
           │
           ▼
┌─────────────────────────────────────┐
│ Capability Implementation           │
│ - Loads parameters via ConfigManager│
│ - Executes algorithm with parameters│
│ - Returns results                   │
└─────────────────────────────────────┘
```

## Detailed Layer Descriptions

### Layer 1: CLI2Json Function (src/main.cpp lines 184-342)

**Purpose**: Transform command-line arguments into structured JSON

**Example Input**:
```bash
./curcuma -confscan conf.xyz -rmsd.method subspace -verbose 3 -threads 4
```

**Parsing Process**:

1. **Extract command keyword**: `-confscan` → `keyword = "confscan"`

2. **Parse arguments into flat dictionary**:
   ```cpp
   key["rmsd"]["method"] = "subspace"   // setNestedJsonValue handles "rmsd.method"
   key["verbose"] = 3
   key["threads"] = 4
   ```

3. **Extract module parameters** (lines 278-306):
   - **NEW (Oct 2025)**: Handles BOTH formats:
     - Flat dotted keys: `"rmsd.method"` (for backward compatibility)
     - Nested structures: `key["rmsd"] = {"method": "subspace"}` (created by setNestedJsonValue)

   ```cpp
   // Check if value is an object AND key is not the current command
   bool is_nested_module = param_value.is_object() && param_name != keyword;

   if (is_nested_module) {
       module_params[param_name] = param_value;  // Extract to top level
       keys_to_remove.push_back(param_name);
   }
   ```

4. **Build controller structure**:
   ```json
   {
     "confscan": { /* parsed confscan params */ },
     "rmsd": { "method": "subspace" },
     "verbose": 3,
     "threads": 4
   }
   ```

**Key Fix (October 26, 2025)**:
- **Before**: Extraction logic looked for `param_name.find('.')` but setNestedJsonValue creates nested structures
- **After**: Now checks both `is_flat_dotted` AND `is_nested_module` conditions
- **Impact**: Multi-module parameters (`-rmsd.method`) now correctly routed to `controller["rmsd"]["method"]`

---

### Layer 2: execute* Functions (src/main.cpp)

**Purpose**: Route to appropriate capability based on command

**Example Functions**:
- `executeConfScan()` - lines 473-510
- `executeRMSD()` - lines 381-395
- `executeOptimization()` - lines 413-471
- `executeSimpleMD()` - lines 549-578

**Architecture Change (October 26, 2025)**:
```cpp
// OLD (redundant manual merging):
json confscan_config = {};
confscan_config["confscan"] = ParameterRegistry::getInstance().getDefaultJson("confscan");
confscan_config["rmsd"] = ParameterRegistry::getInstance().getDefaultJson("rmsd");
// Manual merge of controller parameters...
auto* scan = new ConfScan(confscan_config);

// NEW (simplified, lets ConfigManager handle it):
auto* scan = new ConfScan(controller);  // Pass directly!
```

**Why This Works**:
1. `controller` already has correct structure from CLI2Json:
   - `controller["confscan"]` - confscan params
   - `controller["rmsd"]` - rmsd params

2. ConfScan constructor uses ConfigManager which:
   - Automatically loads defaults from ParameterRegistry
   - Automatically merges user params
   - Single source of truth for defaults
   - No redundant merging!

**Benefits**:
- ✅ Cleaner, simpler code (removed 22 lines of boilerplate)
- ✅ Single responsibility per layer
- ✅ Less chance of bugs from duplicate merging
- ✅ Future-proof: ConfigManager improvements automatically benefit all capabilities

---

### Layer 3: Capability Constructor

**Example: ConfScan Constructor**

```cpp
ConfScan::ConfScan(const json& controller, bool silent)
    : CurcumaMethod(ParameterRegistry::getInstance().getDefaultJson("confscan"),
                    controller, silent),
      m_config(std::vector<std::string>{"confscan", "rmsd"}, controller)
{
    LoadControlJson();
}
```

**Three Important Parts**:

1. **Base Class Initialization (CurcumaMethod)**:
   - Stores defaults from ParameterRegistry
   - Stores raw controller
   - Useful for restart validation, error handling

2. **Multi-Module ConfigManager**:
   - Declared as: `ConfigManager m_config({"confscan", "rmsd"}, controller)`
   - Handles BOTH modules automatically
   - Loads defaults for each, merges user params
   - Provides hierarchical access: `m_config.get<T>("rmsd.method")`

3. **LoadControlJson()**:
   - Reads parameters from ConfigManager
   - Initializes internal data structures
   - Validates parameter values

**Multi-Module ConfigManager Details**:
```cpp
// Constructor that accepts multiple modules
ConfigManager(const std::vector<std::string>& modules, const json& controller)

// For each module:
// 1. Load module defaults: ParameterRegistry::getInstance().getDefaultJson(module_name)
// 2. Extract user params: controller[module_name] if present
// 3. Merge: defaults + user_params (user params take precedence)
// 4. Validate: check aliases, types, required fields
// 5. Store: m_module_configs[module_name] = final_config

// Provides access:
m_config.get<std::string>("confscan.method")    // Get confscan module param
m_config.get<std::string>("rmsd.method")        // Get rmsd module param
m_config.get<int>("threads")                     // Get global param
```

---

### Layer 4: ConfigManager (src/core/config_manager.h/cpp)

**Purpose**: Type-safe, hierarchical parameter access with automatic default merging

**Two Usage Patterns**:

#### Pattern 1: Single Module
```cpp
ConfigManager config("rmsd", controller);
std::string method = config.get<std::string>("method");
```

#### Pattern 2: Multi-Module (for composite capabilities)
```cpp
ConfigManager config({"confscan", "rmsd"}, controller);
std::string rmsd_method = config.get<std::string>("rmsd.method");
std::string confscan_method = config.get<std::string>("confscan.method");
```

**Automatic Default Loading**:
```cpp
// Constructor automatically:
// 1. Loads defaults from ParameterRegistry for EACH module
// 2. Extracts user params from controller
// 3. Merges (user params override defaults)
// 4. Resolves aliases via ParameterRegistry
// 5. Validates types

// Result: Every parameter has a value (default if user didn't specify)
```

**Parameter Access**:
```cpp
// Type-safe access with defaults
T value = config.get<T>("key");
T value = config.get<T>("key", default_value);  // With custom default

// Hierarchical access with dot notation
config.get<std::string>("module.parameter")
config.get<bool>("nested.module.parameter")

// Case-insensitive (via alias resolution)
config.get<int>("MaxIterations")  // Finds "max_iterations"
config.get<int>("max_iterations")  // Both work
```

---

## Common Parameter Routing Issues and Fixes

### Issue 1: Module parameters not reaching capability

**Symptom**: `-rmsd.method subspace` loads default `incr` instead

**Root Cause**: CLI2Json not extracting nested module parameters

**Fix (October 26, 2025)**:
- Updated parameter extraction logic in CLI2Json (lines 285-306)
- Now checks `is_nested_module` in addition to `is_flat_dotted`
- Properly extracts `controller["rmsd"] = {"method": "subspace"}`

**Check**:
```cpp
// Verify parameter made it to controller
if (controller.contains("rmsd") && controller["rmsd"].contains("method")) {
    std::string method = controller["rmsd"]["method"].get<std::string>();
    // method should be "subspace", not "incr"
}
```

### Issue 2: Default parameters missing in capability

**Symptom**: Capability fails because required parameter is null/missing

**Root Cause**: ConfigManager not loading defaults from ParameterRegistry

**Check**:
1. Verify parameter is defined in ParameterRegistry:
   ```bash
   ./curcuma -help | grep -i "your_parameter"
   ```

2. Verify ConfigManager is initialized:
   ```cpp
   ConfigManager m_config("capability_name", controller);
   // Constructor automatically loads defaults
   ```

3. Access parameter correctly:
   ```cpp
   auto value = m_config.get<T>("parameter_name");
   // Will have default if user didn't specify
   ```

### Issue 3: Parameter merging happening twice (redundant)

**Symptom**: Code is slow, excessive logging, potential bugs from double-merge

**Root Cause**: execute* function manually merges, then ConfigManager merges again

**Fix (October 26, 2025)**:
- Removed manual merging from executeConfScan
- Let ConfigManager handle it automatically
- Result: cleaner code, single source of truth

**Pattern**:
```cpp
// ❌ WRONG (redundant):
json config = ParameterRegistry::getInstance().getDefaultJson("module");
for (auto& [k, v] : controller["module"].items()) {
    config[k] = v;
}
capability = new Capability(config);

// ✅ CORRECT:
capability = new Capability(controller);
// Capability constructor uses ConfigManager to merge automatically
```

---

## Best Practices

### For Adding New Capabilities

1. **Define Parameters**:
   - Add PARAM macros in capability header
   - Use snake_case naming
   - Include help text, category, aliases

2. **Constructor**:
   ```cpp
   YourCapability::YourCapability(const json& controller, bool silent)
       : CurcumaMethod(ParameterRegistry::getInstance().getDefaultJson("yourmodule"),
                       controller, silent),
         m_config("yourmodule", controller)  // Single module
   {
       LoadControlJson();  // Initialize from ConfigManager
   }
   ```

3. **Access Parameters**:
   ```cpp
   std::string method = m_config.get<std::string>("method");
   int iterations = m_config.get<int>("max_iterations", 100);  // With default
   ```

4. **For Multi-Module** (like ConfScan + RMSD):
   ```cpp
   m_config(std::vector<std::string>{"confscan", "rmsd"}, controller);

   // Then access both modules:
   auto rmsd_method = m_config.get<std::string>("rmsd.method");
   auto confscan_method = m_config.get<std::string>("confscan.method");
   ```

### For CLI Arguments

- Use dot notation for module-specific parameters:
  ```bash
  ./curcuma -confscan file.xyz -rmsd.method subspace -opt.max_iter 1000
  ```

- Global parameters (no module prefix):
  ```bash
  -verbose 3 -threads 4 -verbosity 2
  ```

- Aliases work automatically:
  ```bash
  -RMSDmethod subspace  # Finds "rmsd.method" via alias resolution
  -max_iterations 100   # Works with -MaxIterations too
  ```

### For Debugging Parameter Issues

1. **Check CLI2Json routing**:
   ```cpp
   // Add temporary debug in executeXxx:
   std::cerr << "controller keys: ";
   for (auto& [k, v] : controller.items()) {
       std::cerr << k << " ";
   }
   std::cerr << std::endl;
   ```

2. **Check ConfigManager merging**:
   ```cpp
   // Define DEBUG_CONFIG_MANAGER to see merge results
   // In config_manager.cpp around line 313-315
   ```

3. **Verify ParameterRegistry defaults**:
   ```bash
   ./curcuma -help | grep -i parameter_name
   ```

4. **Check capability receives parameters**:
   ```cpp
   // In LoadControlJson():
   auto value = m_config.get<T>("parameter");
   std::cerr << "Parameter value: " << value << std::endl;
   ```

---

## Architecture Summary

| Layer | Responsibility | Files | Key Classes |
|-------|---|---|---|
| **CLI2Json** | Parse args → controller dict | src/main.cpp | setNestedJsonValue() |
| **execute*** | Route to capability | src/main.cpp | executeConfScan(), etc. |
| **Capability** | Initialize with ConfigManager | src/capabilities/*.cpp | ConfScan, RMSD, etc. |
| **ConfigManager** | Merge defaults + load params | src/core/config_manager.* | ConfigManager |
| **ParameterRegistry** | Store param definitions | src/core/parameter_registry.* | ParameterRegistry |

**Golden Rule**:
> Each layer has ONE clear responsibility. Don't do work in a layer that should happen elsewhere. Trust the layers below you.

- ✅ execute* functions: Just pass controller (don't merge)
- ✅ Capability constructors: Initialize ConfigManager (let it merge)
- ✅ ConfigManager: Load defaults + merge (provide clean API)
- ✅ ParameterRegistry: Store definitions (single source of truth)

---

## Migration Checklist

When refactoring parameter handling in existing capabilities:

- [ ] Identify all PARAM definitions and move to header
- [ ] Create ParameterRegistry entries via build system
- [ ] Update capability constructor:
  - [ ] Remove `m_defaults` static JSON
  - [ ] Add ConfigManager initialization
  - [ ] Specify modules (single or multi)
- [ ] Update parameter access:
  - [ ] Replace `Json2KeyWord()` calls with `m_config.get<T>()`
  - [ ] Use hierarchical dot notation for clarity
- [ ] Remove manual default merging in execute* functions
- [ ] Test with various parameter combinations:
  - [ ] No parameters (all defaults)
  - [ ] Single parameter override
  - [ ] Multiple parameter overrides
  - [ ] Aliases work correctly
- [ ] Run help output to verify:
  - [ ] Parameter descriptions appear
  - [ ] Types are correct
  - [ ] Aliases are listed

---

**Document Version**: 1.0
**Last Updated**: 2025-10-26
**Author**: Claude (AI Code Assistant)
**Related Files**: PARAMETER_SYSTEM.md, PARAMETER_MIGRATION_GUIDE.md
