# Curcuma Build System - Implementation Roadmap

**Version**: 1.0
**Date**: November 2025
**Status**: Partially Fixed - Build 2 Working, Roadmap for Remaining 4 Builds
**Priority**: HIGH - Enables flexible build configurations for different use cases

---

## Executive Summary

The Curcuma build system has 5 primary configurations controlled by CMake options:
- `USE_TBLITE`, `USE_ULYSSES`, `USE_XTB`, `USE_D3`, `USE_D4`

**Current Status:**
- ✅ **Build 2 (Standard)** - FULLY WORKING
- ❌ **Builds 1, 3, 4, 5** - BLOCKED by conditional compilation issues

**Root Cause:** Header guards (`#ifdef`) protect includes but not implementations. When a library is disabled, source files still try to use classes/functions that don't exist.

**Solution:** Systematic protection of all dependent code with conditional compilation guards.

---

## Current Build Configurations

### Build 1: Minimal (EHT + UFF Only)
```cmake
-DUSE_TBLITE=OFF
-DUSE_ULYSSES=OFF
-DUSE_D3=OFF
-DUSE_XTB=OFF
-DUSE_D4=OFF
```
**Available Methods**: UFF, EHT (native only)
**Status**: ❌ BLOCKED - gfnff.cpp depends on D3

### Build 2: Standard (Recommended) ✅
```cmake
-DUSE_TBLITE=ON
-DUSE_ULYSSES=ON
-DUSE_D3=ON
-DUSE_XTB=OFF
-DUSE_D4=OFF
```
**Available Methods**: UFF, EHT, GFN2 (TBLite), GFN1, iPEA1, PM3, PM6, AM1, MNDO
**Status**: ✅ WORKING - TESTED AND VERIFIED

### Build 3: Full QM (Maximum Methods)
```cmake
-DUSE_TBLITE=ON
-DUSE_ULYSSES=ON
-DUSE_D3=ON
-DUSE_XTB=ON
-DUSE_D4=OFF
```
**Available Methods**: All from Build 2 + XTB (GFN-FF, GFN1, GFN2)
**Status**: ❌ BLOCKED - Same D3/D4 header chain issue

### Build 4: With D4 (Dispersion)
```cmake
-DUSE_TBLITE=ON
-DUSE_ULYSSES=ON
-DUSE_D3=ON
-DUSE_D4=ON
-DUSE_XTB=OFF
```
**Available Methods**: Build 2 methods + D4 dispersion
**Status**: ❌ BLOCKED - D3/D4 dependency chain

### Build 5: TBLite Only (Single Library)
```cmake
-DUSE_TBLITE=ON
-DUSE_ULYSSES=OFF
-DUSE_D3=OFF
-DUSE_XTB=OFF
-DUSE_D4=OFF
```
**Available Methods**: GFN2, GFN1, iPEA1 (TBLite only)
**Status**: ❌ BLOCKED - D3 included indirectly

---

## Problem Analysis

### Issue: Cascading Include Chain

```
gfnff.cpp (ALWAYS COMPILED)
  ↓
#include "forcefield.h"
  ↓
#include "forcefieldthread.h"
  ↓ [Line 28-30 in forcefieldthread.h]
#ifdef USE_D3
  #include "dftd3interface.h"    ← Header is conditional
#endif
  ↓
#include "s-dftd3.h"             ← INNER include NOT conditional!
```

When `USE_D3=OFF`:
1. CMake skips dftd3interface.cpp compilation ✓
2. But forcefieldthread.h still includes dftd3interface.h ✓
3. And dftd3interface.h always includes "s-dftd3.h" ✗
4. The file doesn't exist when D3 is disabled → **COMPILATION ERROR**

### Files Involved

**Protected (Fixes Applied):**
- ✅ `dftd3interface.h` - Full header wrapped with `#ifdef USE_D3`
- ✅ `dftd4interface.h` - Full header wrapped with `#ifdef USE_D4`
- ✅ `forcefieldthread.h` - D3Thread class wrapped with `#ifdef USE_D3`
- ✅ `forcefieldthread.cpp` - D3Thread implementation wrapped with `#ifdef USE_D3`

**Still Unprotected (Causing Failures):**
- ❌ `gfnff.cpp` - Always compiled, indirectly includes D3
- ❌ `gfnff.h` - Includes forcefield.h unconditionally
- ❌ `forcefield.cpp` - Instantiates D3Thread without guards (line ~470)
- ❌ `H4Thread` - No guards at all (similar to D3Thread)

---

## Fix Strategy

### Phase 1: Immediate Fixes (Get Builds 1, 3, 4, 5 Working)

#### 1.1 Protect gfnff.cpp Compilation

**File**: `src/core/energy_calculators/qm_methods/gfnff.cpp`

Option A: Conditional Compilation
```cpp
#ifndef USE_D3
#ifndef USE_D4
// If neither D3 nor D4, this file should not include the problematic headers
#define GFNFF_STANDALONE 1
#endif
#endif
```

Option B (Recommended): Restructure Includes
```cpp
// Only include what's needed
#include "gfnff.h"
// Don't include forcefield.h directly - use forward declarations
class ForceField;  // Forward declare
class D3Thread;    // Forward declare
```

**Priority**: 🔴 CRITICAL - Blocks 4/5 builds

#### 1.2 Protect D3Thread/H4Thread Instantiation

**File**: `src/core/energy_calculators/ff_methods/forcefield.cpp`
**Location**: Line ~466-486

```cpp
#ifdef USE_D3
int d3 = false;  // m_parameters["d3"];
if (d3) {
    if (free_threads > 1)
        free_threads--;
    D3Thread* thread = new D3Thread(m_threads - 1, free_threads);
    thread->setParamater(m_parameters);
    thread->Initialise(m_atom_types);
    m_threadpool->addThread(thread);
    m_stored_threads.push_back(thread);
}
#endif

#ifdef USE_D4
// Similar pattern for D4 (currently all zeros)
#endif
```

**Priority**: 🔴 CRITICAL - Runtime failures even if headers compile

#### 1.3 Add H4Thread Guards

**File**: `src/core/energy_calculators/ff_methods/forcefieldthread.h/cpp`

Same pattern as D3Thread:
```cpp
#ifdef USE_HBONDS4
class H4Thread : public ForceFieldThread { ... };
#endif
```

**Note**: Need to define `USE_HBONDS4` CMake option or use existing flag

**Priority**: 🟡 HIGH - Prevents potential future issues

### Phase 2: Architecture Refactoring (Clean Design)

#### 2.1 Decouple Dispersion from Core ForceField

**Current Problem**: `forcefield.h` pulls in dispersion-related classes even when not used

**Solution**: Create separate header for optional features
```
src/core/energy_calculators/ff_methods/dispersion_threads.h
  - Contains D3Thread, H4Thread classes
  - Only included when USE_D3 || USE_D4 || USE_HBONDS4
  - Reduces pollution of forcefield.h
```

**Benefit**: Cleaner separation of concerns

#### 2.2 Factory Pattern for Thread Creation

**Instead of**: Direct instantiation in forcefield.cpp
```cpp
D3Thread* thread = new D3Thread(...);  // Hard-coded dependency
```

**Use**: Factory pattern
```cpp
class DispersionThreadFactory {
    static std::unique_ptr<ForceFieldThread> createD3Thread(...) {
        #ifdef USE_D3
        return std::make_unique<D3Thread>(...);
        #else
        return nullptr;  // or throw exception
        #endif
    }
};

// Usage:
auto thread = DispersionThreadFactory::createD3Thread(...);
if (thread) m_threadpool->addThread(thread.release());
```

**Benefit**: Single point of control for optional features

#### 2.3 CMake-Level Validation

**Add to CMakeLists.txt**:
```cmake
# Warn about incompatible option combinations
if(USE_D4 AND NOT USE_D3)
    message(WARNING "D4 requires D3 - enabling USE_D3=ON")
    set(USE_D3 ON CACHE BOOL "" FORCE)
endif()

# Document what methods are available
message(STATUS "Available QM Methods:")
if(USE_TBLITE)
    message(STATUS "  - TBLite: GFN2, GFN1, iPEA1")
endif()
if(USE_ULYSSES)
    message(STATUS "  - Ulysses: PM3, PM6, AM1, MNDO")
endif()
# ... etc
```

**Benefit**: Prevents invalid configurations early

### Phase 3: Testing Infrastructure

#### 3.1 Automated Build Testing

**Location**: `/home/conrad/src/curcuma/build_test/`

**Script**: `test_all_builds.sh`
```bash
#!/bin/bash
BUILDS=(
    "build_1_minimal:OFF:OFF:OFF:OFF:OFF"
    "build_2_standard:ON:ON:OFF:ON:OFF"
    "build_3_full_qm:ON:ON:ON:ON:OFF"
    "build_4_d4:ON:ON:OFF:ON:ON"
    "build_5_tblite_only:ON:OFF:OFF:OFF:OFF"
)

for config in "${BUILDS[@]}"; do
    IFS=':' read -r name tblite ulysses xtb d3 d4 <<< "$config"
    echo "Building $name..."
    cd "build_$name"
    cmake .. -DUSE_TBLITE=$tblite -DUSE_ULYSSES=$ulysses \
             -DUSE_XTB=$xtb -DUSE_D3=$d3 -DUSE_D4=$d4
    make -j4 && echo "✅ $name SUCCESS" || echo "❌ $name FAILED"
done
```

#### 3.2 CI/CD Integration

**For GitHub Actions**:
```yaml
name: Build Matrix Test
on: [push, pull_request]
jobs:
  build:
    strategy:
      matrix:
        config:
          - { name: minimal, tblite: 'OFF', ulysses: 'OFF', d3: 'OFF', d4: 'OFF', xtb: 'OFF' }
          - { name: standard, tblite: 'ON', ulysses: 'ON', d3: 'ON', d4: 'OFF', xtb: 'OFF' }
          - { name: full-qm, tblite: 'ON', ulysses: 'ON', d3: 'ON', d4: 'OFF', xtb: 'ON' }
          # ... more configs
    steps:
      - uses: actions/checkout@v3
      - run: |
          mkdir build_${{ matrix.config.name }}
          cd build_${{ matrix.config.name }}
          cmake .. -DUSE_TBLITE=${{ matrix.config.tblite }} ...
          make -j4
```

**Benefit**: Catch regressions immediately

---

## Implementation Timeline

### Immediate (This Week)
- [ ] Phase 1.1: Protect gfnff.cpp compilation
- [ ] Phase 1.2: Guard D3Thread instantiation in forcefield.cpp
- [ ] Phase 1.3: Add H4Thread guards
- [ ] Test all 5 builds - should all compile successfully

### Short-term (Next 1-2 Weeks)
- [ ] Phase 2.1: Create dispersion_threads.h
- [ ] Phase 2.2: Implement factory pattern
- [ ] Phase 2.3: Add CMake validation
- [ ] Document build options in README

### Medium-term (Next Month)
- [ ] Phase 3.1: Automated local build test script
- [ ] Phase 3.2: CI/CD integration (GitHub Actions)
- [ ] Performance benchmarking across configurations

---

## Verification Checklist

Once fixes are implemented, verify:

### Compilation Verification
- [ ] Build 1 compiles without D3/D4 errors
- [ ] Build 2 compiles with standard options
- [ ] Build 3 compiles with XTB enabled
- [ ] Build 4 compiles with D4 enabled
- [ ] Build 5 compiles with TBLite only

### Functional Verification
```bash
# Test each build
cd build_2_standard
./curcuma -sp water.xyz -method uff     # Should work
./curcuma -sp water.xyz -method gfn2    # Should work
./curcuma -rmsd ref.xyz target.xyz      # Should work
```

### Regression Testing
```bash
ctest --output-on-failure  # All 26 tests should pass in each build
```

---

## References

### Related Documentation
- `TODO.md` - High-level task list
- `CLAUDE.md` - Project overview
- `src/core/CLAUDE.md` - Core subsystem documentation

### Build Configuration Details
```cmake
# Main CMakeLists.txt controls these options
option(USE_TBLITE "Enable TBLite tight-binding DFT")
option(USE_ULYSSES "Enable Ulysses semi-empirical")
option(USE_XTB "Enable XTB extended tight-binding")
option(USE_D3 "Enable DFT-D3 dispersion")
option(USE_D4 "Enable DFT-D4 dispersion")
```

### Method Availability Matrix

| Method | Build 1 | Build 2 | Build 3 | Build 4 | Build 5 |
|--------|---------|---------|---------|---------|---------|
| UFF | ✅ | ✅ | ✅ | ✅ | ✅ |
| EHT | ✅ | ✅ | ✅ | ✅ | ✅ |
| GFN2 (TBLite) | ❌ | ✅ | ✅ | ✅ | ✅ |
| GFN1 (TBLite) | ❌ | ✅ | ✅ | ✅ | ✅ |
| iPEA1 | ❌ | ✅ | ✅ | ✅ | ✅ |
| GFN-FF (XTB) | ❌ | ❌ | ✅ | ❌ | ❌ |
| PM3/PM6/AM1/MNDO | ❌ | ✅ | ✅ | ✅ | ❌ |
| D3 Correction | ❌ | ✅ | ✅ | ✅ | ❌ |
| D4 Correction | ❌ | ❌ | ❌ | ✅ | ❌ |

---

## Success Metrics

✅ **Achieved**:
- Build 2 (Standard) fully functional
- UFF and GFN2 methods verified working
- Root cause identified and partially fixed

🎯 **Upcoming**:
- All 5 builds compile cleanly
- All 26 CLI tests pass on each build configuration
- CI/CD pipeline validates builds on every commit
- Clean separation of optional dependencies

---

**Last Updated**: November 2025
**Status**: Roadmap Document - Ready for Implementation
