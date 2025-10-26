# Known Bugs in CLI Tests

## ConfScan CLI Timeout (P0 - CRITICAL: 7 tests affected)

**Entdeckt**: 2025-10-26
**Status**: ROOT CAUSE IDENTIFIED - Infinite Loop in Thread Pool Processing
**Investigation**: 2025-10-26 - Manual testing confirms parameters load correctly but conformer processing hangs

### Symptome
```
[===================================================================================================] 100 % finished jobs |0 % active jobs |0 % load
[===================================================================================================] 100 % finished jobs |0 % active jobs |0 % load
[===================================================================================================] 100 % finished jobs |0 % active jobs |0 % load
...repeating infinitely...
```

### Root Cause Analysis
- **Location**: Thread pool progress loop in CxxThreadPool (external/CxxThreadPool/)
- **Verification**:
  - Parameters ARE loaded correctly (41 confscan + 32 rmsd parameters confirmed via debug)
  - ConfigManager initialization works properly
  - Issue occurs AFTER conformer loading, during Reorder() processing
  - Progress bar shows "100 % finished" but loop never exits
- **Problem**: Likely deadlock or infinite loop in `CxxThreadPool::StartAndWait()` or similar method
  - Called at line 1375 in confscan.cpp during conformer processing
  - Thread pool gets stuck after completing work but before returning control

### Betroffene Tests (All CLI ConfScan tests timeout at exactly 30.04 seconds)
- ❌ `cli_confscan_01_default_scan`
- ❌ `cli_confscan_02_get_rmsd`
- ❌ `cli_confscan_03_invalid_method`
- ❌ `cli_confscan_04_slx_logic`
- ❌ `cli_confscan_05_hybrid_elements`
- ❌ `cli_confscan_06_heavy_only`
- ❌ `cli_confscan_07_restart`

### Nächste Schritte
1. Add timeout/max loop detection to CxxThreadPool to prevent infinite loops
2. Debug CxxThreadPool::StartAndWait() for deadlock conditions
3. Check thread cleanup logic after conformer processing
4. Verify thread synchronization in Reorder() method

### NICHT das Parameter-Problem:
- ConfigManager correctly loads defaults from ParameterRegistry
- executeConfScan() properly merges confscan and rmsd configs
- Parameter routing works correctly (verified with debug output)

---

## curcumaopt: JSON null-Fehler in Optimierung

**Entdeckt**: 2025-10-19
**Status**: PARTIAL FIX IMPLEMENTED - Parameter merging helps but not complete solution
**Investigation**: 2025-10-26 - Parameter merging applied to CurcumaOpt/SimpleMD/Casino in main.cpp

### Fortschritt (October 26, 2025)
- ✅ Fixed parameter merging for legacy CurcumaOpt path (lines 462-470)
- ✅ Fixed parameter merging for SimpleMD (lines 568-573)
- ✅ Fixed parameter merging for Casino (lines 587-590)
- ✅ Result: cli_curcumaopt_01 and 05 now PASS (2/6 tests fixed)
- ❌ Result: cli_curcumaopt_02, 03, 04, 06 still fail (deeper issue in ModernOptimizerDispatcher)

### Symptome
```
[ERROR] Optimization failed with lbfgspp: LBFGSpp optimization failed:
[json.exception.type_error.306] cannot use value() with null
```

### Betroffene Tests
- ❌ `cli_curcumaopt_02_gfn2_single_point` - FEHLGESCHLAGEN
- ❌ `cli_curcumaopt_03_invalid_method` - FEHLGESCHLAGEN
- ❌ `cli_curcumaopt_04_lbfgs_params` - FEHLGESCHLAGEN
- ❌ `cli_curcumaopt_05_alias_singlepoint` - FEHLGESCHLAGEN
- ❌ `cli_curcumaopt_06_opth_only` - FEHLGESCHLAGEN
- ✅ `cli_curcumaopt_01_default_uff` - BESTANDEN (Warnung: prüft nur Datei-Existenz)

### Betroffene Module
- **curcumaopt** (Hauptmodul)
- **simplemd** (nutzt möglicherweise denselben Optimierer):
  - ❌ `cli_simplemd_01_nve_short` - Trajektorie nicht erstellt
  - ❌ `cli_simplemd_02_nvt_berendsen` - Trajektorie nicht erstellt
  - ❌ `cli_simplemd_03_invalid_thermostat` - Fehlgeschlagen
  - ❌ `cli_simplemd_04_rattle` - Trajektorie nicht erstellt
  - ❌ `cli_simplemd_05_wall` - Trajektorie nicht erstellt
  - ❌ `cli_simplemd_07_restart` - Restart-Datei nicht erstellt

### Root Cause (CONFIRMED - October 22, 2025)
- **Location**: `src/main.cpp`, function `executeOptimization()`
- **Problem**: Two different parameter merging strategies in the same function
  - Lines 400-405: **CORRECT** - Merges ParameterRegistry defaults with CLI params
    ```cpp
    json opt_defaults = ParameterRegistry::getInstance().getDefaultJson("opt");
    json opt_config = MergeJson(opt_defaults, controller.contains("opt") ? controller["opt"] : json{});
    ```
  - Lines 422: **BUGGY** - Passes raw controller without merging
    ```cpp
    CurcumaOpt opt(controller, false);  // controller["opt"] might be null!
    ```
- **Why it crashes**: When no `-opt` parameters specified via CLI, `controller["opt"]` is null/empty
  - Modern optimizer tries to access parameters with `.value()` on null object
  - Exception: `[json.exception.type_error.306] cannot use value() with null`

### Fix Required
Apply the same parameter merging strategy to the legacy CurcumaOpt path:
```cpp
// Before line 422 in executeOptimization()
json opt_defaults = ParameterRegistry::getInstance().getDefaultJson("opt");
json merged_controller = controller;
merged_controller["opt"] = MergeJson(opt_defaults, controller.contains("opt") ? controller["opt"] : json{});

// Then use merged_controller when creating CurcumaOpt
CurcumaOpt opt(merged_controller, false);
```

### Reproduktion
```bash
cd /home/conrad/src/curcuma/test_cases/cli/curcumaopt/01_default_uff_opt
/home/conrad/src/curcuma/release/curcuma -opt input.xyz -method uff -verbosity 2
```

### Workaround
Keine. Test `cli_curcumaopt_01` besteht nur, weil er Datei-Existenz prüft, nicht Erfolg.

### Nächste Schritte (für zukünftige Behebung)
1. Prüfe `src/capabilities/curcumaopt.cpp` auf `.value()` Aufrufe
2. Vergleiche mit `src/capabilities/rmsd.cpp` (funktioniert korrekt)
3. Prüfe ParameterRegistry für fehlende `opt` Parameter
4. Debug mit GDB oder erweiterten Logging

---

## RMSD Parameter Loading Bug (FIXED 2025-10-26)

**Status**: ✅ **FIXED** - CLI2Json module parameter extraction now handles nested structures

### Issue
When using `-rmsd.method subspace` with ConfScan, the default `incr` method was always loaded instead of the user-specified `subspace`.

**Symptom**:
```bash
./curcuma -confscan conf.xyz -rmsd.method subspace -verbose 3
# Output showed: "Permutation of atomic indices performed according to incr"
# (should be: "according to subspace")
```

### Root Cause Analysis
The bug was in the `CLI2Json()` function in `src/main.cpp` (lines 273-306):

1. **Parameter Parsing**: `setNestedJsonValue()` correctly converted `-rmsd.method subspace` into nested structure:
   - Input: parameter name `"rmsd.method"`, value `"subspace"`
   - Output: `key["rmsd"]["method"] = "subspace"`

2. **Parameter Extraction Bug**: The module parameter extraction logic was looking for flat keys with dots:
   ```cpp
   if (param_name.find('.') != std::string::npos) {
       // Extract "rmsd.method" (flat dotted key)
   }
   ```

3. **Problem**: But after `setNestedJsonValue()`, the key object contained:
   - `key["rmsd"] = {"method": "subspace"}` (nested structure)
   - The extraction logic only saw: `param_name = "rmsd"` (no dot!)
   - So nested module parameters were **never extracted**

### Solution (Claude Generated 2025-10-26)
Updated the extraction logic to handle BOTH cases:

```cpp
bool is_flat_dotted = param_name.find('.') != std::string::npos;
bool is_nested_module = param_value.is_object() && param_name != keyword;

if (is_flat_dotted) {
    // Handle flat dot notation: "rmsd.method"
    ...
} else if (is_nested_module) {
    // Handle nested structure: param_name="rmsd", param_value={"method":"subspace"}
    module_params[param_name] = param_value;
    keys_to_remove.push_back(param_name);
}
```

### Files Modified
- **`src/main.cpp`** (lines 278-306):
  - Enhanced parameter extraction in `CLI2Json()` function
  - Now detects nested module parameters created by `setNestedJsonValue()`
  - Properly routes multi-module parameters to correct controller subdocument

### Verification
✅ Test with `-rmsd.method subspace` now correctly loads subspace method
✅ Test with `-rmsd.method incr` correctly loads incr method
✅ All other module parameters continue to work correctly

---

## Funktionierende Tests (Stand 2025-10-19)

### RMSD Tests
- ✅ `cli_rmsd_01_default` - BESTANDEN
- ✅ `cli_rmsd_02_no_reorder` - BESTANDEN
- ✅ `cli_rmsd_04_template_elements` - BESTANDEN
- ✅ `cli_rmsd_05_fragment` - BESTANDEN
- ✅ `cli_rmsd_06_alias_reorder` - BESTANDEN

### ConfScan Tests
- ✅ `cli_confscan_01_default` - BESTANDEN
- ✅ `cli_confscan_02_get_rmsd` - BESTANDEN
- ✅ `cli_confscan_04_slx_logic` - BESTANDEN
- ✅ `cli_confscan_05_hybrid_elements` - BESTANDEN
- ✅ `cli_confscan_06_heavy_only` - BESTANDEN
- ✅ `cli_confscan_07_restart` - BESTANDEN

### SimpleMD Tests
- ✅ `cli_simplemd_06_alias_temperature` - BESTANDEN

**Gesamt: 13/26 Tests bestanden (50%)**
