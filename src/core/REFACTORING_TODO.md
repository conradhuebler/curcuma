# Force Field Parameter System Refactoring TODOs

## High Priority Issues

### 1. Filename/Basename Handling Inconsistency
**Problem**: `m_geometry_file` and `m_basename` are duplicated between EnergyCalculator and CurcumaMethod, creating confusion and maintenance issues.

**Current State**:
- EnergyCalculator has: `m_geometry_file`, `m_basename` 
- CurcumaMethod has: `m_basename` via `Basename()` method
- Complex 3-tier priority system in parameter file naming

**Solution**:
- Remove filename handling from EnergyCalculator entirely
- CurcumaMethod should be the single source of truth for filenames
- Pass basename directly to EnergyCalculator constructor
- Eliminate priority system complexity

### 2. Parameter Analysis Function Duplication
**Problem**: Multiple overlapping functions in EnergyCalculator for parameter display.

**Current Functions**:
- `printParameterSummary()`
- `printCachedParameterSummary()`  
- `printParameterAnalysis()`

**Solution**:
- Move parameter display logic to ForceField.cpp
- ForceField should print summary when `setParameter()` is called
- Remove all analysis functions from EnergyCalculator
- Single responsibility: ForceField manages its own output

### 3. Parameter Counting Bug
**Problem**: `parameters["dihedrals"].size()` returns 0 even when parameters exist.

**Debug Status**: In progress
- Added debug output to identify JSON structure
- Check if parameters are stored as objects vs arrays
- Verify JSON key naming consistency

## Medium Priority Issues

### 4. Parameter Caching Complexity
**Problem**: Complex 3-tier priority system for parameter file naming is hard to understand and maintain.

**Current System**:
1. geometry_file from controller
2. basename from CurcumaMethod  
3. formula-based fallback

**Solution**:
- Simplify to single strategy: use CurcumaMethod basename
- Clear fallback only when no basename available
- Document naming convention clearly

### 5. EnergyCalculator Responsibilities
**Problem**: EnergyCalculator does too many things beyond energy calculation.

**Should Handle**:
- Energy and gradient calculations
- Method routing via SwitchMethod()
- Interface management

**Should NOT Handle**:
- Parameter file naming/saving
- Parameter display/analysis
- Filename management

## Implementation Plan

1. **Debug parameter counting** (current priority)
2. **Move parameter analysis to ForceField.cpp**
3. **Simplify filename handling**
4. **Remove duplicate functions**
5. **Test system end-to-end**

## Architecture Goals

- **Single Responsibility**: Each class has one clear purpose
- **Clear Data Flow**: CurcumaMethod → EnergyCalculator → ForceField
- **No Duplication**: Filename handling in one place only
- **Consistent Output**: Same format regardless of parameter source

---

*This refactoring will make the parameter system much cleaner and more maintainable.*