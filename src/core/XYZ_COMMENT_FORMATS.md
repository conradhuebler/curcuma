# XYZ Comment Line Formats - Production Requirements

## Overview
This document catalogs the specific XYZ comment line formats used in production workflows. **These formats MUST NOT break during refactoring.**

## Critical Production Formats

### 1. ORCA Format
```
"Coordinates from ORCA-job input E -674.785305664016"
```
- **Pattern**: `"Coordinates from ORCA-job input E " + <energy>`
- **Energy**: `-674.785305664016` Hartree
- **Usage**: Output from ORCA quantum chemistry calculations
- **Parser**: Should extract energy value for `mol.setEnergy()`

### 2. XTB Format  
```
" energy: -22.066967268618 gnorm: 0.057841939709 xtb: 6.6.0 (8843059)"
```
- **Pattern**: `" energy: " + <energy> + " gnorm: " + <gradient_norm> + " xtb: " + <version>`
- **Energy**: `-22.066967268618` Hartree
- **Gradient Norm**: `0.057841939709` (convergence indicator)
- **Version**: `6.6.0 (8843059)` (XTB version and commit)
- **Usage**: Output from XTB tight-binding calculations
- **Parser**: Should extract energy and optionally gradient norm

### 3. Simple Energy Format
```
"-3323.538813022354"
```
- **Pattern**: Just a floating-point number (positive or negative)
- **Energy**: `-3323.538813022354` Hartree
- **Usage**: Common in optimization trajectories, minimal output
- **Parser**: Direct energy extraction with `std::stod()`

### 4. Empty Comment
```
""
```
- **Pattern**: Empty string or whitespace only
- **Usage**: Default/minimal XYZ files
- **Parser**: Should not crash, no energy extraction

### 5. Multi-Field Format  
```
"Step 42 E: -123.456 RMS: 0.001 MAX: 0.005"
```
- **Pattern**: Multiple fields with keywords
- **Energy**: `-123.456` (after "E:")
- **RMS**: `0.001` (RMS gradient)
- **MAX**: `0.005` (Max gradient component)
- **Usage**: Optimization step information
- **Parser**: Should extract energy after "E:" keyword

## Parser Requirements

### Backward Compatibility
- All existing comment formats MUST continue to work
- No changes to molecule energy values after parsing
- No crashes on any comment format (graceful fallback)

### Energy Extraction Priority
1. **ORCA**: Extract energy after "E " 
2. **XTB**: Extract energy after "energy: "
3. **Simple**: Parse entire string as double
4. **Multi-Field**: Extract energy after "E:"
5. **Empty**: No extraction, keep existing energy
6. **Unknown**: No extraction, store comment as-is

### Format Detection Logic
```cpp
CommentFormat detectFormat(const std::string& comment) {
    if (comment.empty() || isWhitespaceOnly(comment)) 
        return EMPTY;
    if (comment.find("Coordinates from ORCA-job input E") != std::string::npos)
        return ORCA_FORMAT;
    if (comment.find(" energy:") != std::string::npos && comment.find("gnorm:") != std::string::npos)
        return XTB_FORMAT;  
    if (comment.find("E:") != std::string::npos)
        return MULTI_FIELD;
    if (isSimpleNumber(comment))
        return SIMPLE_ENERGY;
    return UNKNOWN;
}
```

## Test Validation

### Required Test Cases
Each format must be tested with:
- **Parsing Success**: No crashes or exceptions
- **Energy Extraction**: Correct energy value extracted (where applicable)
- **Robustness**: Variations in whitespace, case, etc.

### Edge Cases
- Leading/trailing whitespace
- Scientific notation: `"-1.234567e-03"`
- Very large/small numbers
- Invalid number formats (should not crash)
- Unicode/special characters (graceful fallback)

## Legacy Parser Functions
The current implementation uses 10 separate functions:
- `setXYZComment_0()` through `setXYZComment_8()`
- `setXYZComment_10()` (note: `_9()` is missing!)

These will be unified into a single parser with format detection, but **all functionality must be preserved**.

## Migration Strategy
1. **Document current behavior**: Run tests to establish baseline
2. **Identify format patterns**: Map each `setXYZComment_X()` to format type
3. **Implement unified parser**: Single entry point with format detection
4. **Validate**: All production formats work identically
5. **Remove legacy code**: Only after validation passes

---

**CRITICAL**: Any change that breaks these production formats will cause workflow failures. Test thoroughly before deployment.