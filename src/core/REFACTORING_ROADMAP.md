# Molecule Class Refactoring Roadmap

## Overview
This document outlines the planned refactoring of the Molecule class with specific test validations for each phase. All changes are designed to maintain API compatibility while improving performance and code quality.

## Pre-Refactoring Baseline
**CRITICAL**: Run comprehensive test suite to establish baseline:
```bash
cd build
make test_molecule
./test_molecule
```
Expected output: All tests pass, establishing current behavior as reference.

## Phase 1: XYZ Comment Parser Unification (HIGH PRIORITY)

### Current Problem
- 10 separate `setXYZComment_X()` functions (X = 0,1,2,3,4,5,6,7,8,10)
- Code duplication and maintenance nightmare
- Missing `setXYZComment_9()` function (potential bug)

### CRITICAL: Production Comment Formats (MUST NOT BREAK!)
```
1. ORCA:   "Coordinates from ORCA-job input E -674.785305664016"
2. XTB:    " energy: -22.066967268618 gnorm: 0.057841939709 xtb: 6.6.0 (8843059)"
3. Simple: "-3323.538813022354"
4. Empty:  ""
5. Multi:  "Step 42 E: -123.456 RMS: 0.001 MAX: 0.005"
```

### Refactoring Plan
```cpp
// Replace current approach with Strategy Pattern
class XYZCommentParser {
    enum class CommentFormat {
        ORCA_FORMAT,        // "Coordinates from ORCA-job input E <energy>"
        XTB_FORMAT,         // " energy: <E> gnorm: <grad> xtb: <version>"
        SIMPLE_ENERGY,      // "<energy>" (just a number)
        MULTI_FIELD,        // "Step X E: Y RMS: Z"
        EMPTY,              // "" or whitespace only
        UNKNOWN             // Fallback - store as-is
    };
    
    static CommentFormat detectFormat(const std::string& comment);
    static bool parseComment(const std::string& comment, Molecule& mol);
    
    // Specific parsers for each format
    static bool parseORCA(const std::string& comment, Molecule& mol);
    static bool parseXTB(const std::string& comment, Molecule& mol);
    static bool parseSimpleEnergy(const std::string& comment, Molecule& mol);
};
```

### Test Validation
- `test_xyz_comment_parsing()` validates ALL production formats
- **CRITICAL**: ORCA, XTB, and simple energy formats must work unchanged
- Energy extraction should work for parseable formats
- Unknown formats should not crash (graceful fallback)
- Performance should be equal or better

## Phase 2: Granular Cache System (HIGH PRIORITY)

### Current Problem
- Single `m_dirty` flag invalidates ALL caches
- Changing charges invalidates geometry-based caches (unnecessary)
- No fine-grained cache control

### Refactoring Plan
```cpp
enum class CacheType {
    GEOMETRY = 0,     // Distance matrix, topology
    CONNECTIVITY = 1, // Fragments, bonds  
    PROPERTIES = 2,   // Charges, dipole
    ANALYSIS = 3      // Persistent diagrams, etc.
};

class CacheManager {
    mutable std::bitset<4> m_cache_valid;
    mutable std::array<std::any, 4> m_cache_storage;
    
public:
    void invalidate(CacheType type);
    void invalidateAll();
    template<typename T> bool isValid(CacheType type) const;
    template<typename T> T& getCache(CacheType type) const;
};
```

### Test Validation
- `test_cache_granularity()` should show improved behavior
- Charge changes should NOT invalidate geometry caches
- Performance should improve for mixed operations

## Phase 3: Fragment System Optimization (MEDIUM PRIORITY)

### Current Problem
- `std::map<int, int> m_fragment_assignment` → O(log n) lookups
- Fragment lookups inefficient for large molecules

### Refactoring Plan
```cpp
struct FragmentInfo {
    std::vector<std::vector<int>> fragments;        // Fragment → atoms
    std::vector<int> atom_to_fragment;             // Atom → fragment (O(1)!)
    std::vector<double> fragment_masses;           // Cached masses
    mutable bool valid = false;
};
```

### Test Validation
- `test_fragment_system_performance()` should show speed improvement
- All fragment-based operations must remain correct
- Memory usage should be similar or better

## Phase 4: Type-Safe Elements (MEDIUM PRIORITY)

### Current Problem
- `std::vector<int> m_atoms` → no type safety
- Magic numbers for elements (6 = Carbon, 8 = Oxygen)
- No compile-time validation

### Refactoring Plan
```cpp
enum class ElementType : uint8_t {
    H = 1, He = 2, Li = 3, Be = 4, B = 5, C = 6, N = 7, O = 8, F = 9,
    Ne = 10, Na = 11, Mg = 12, Al = 13, Si = 14, P = 15, S = 16, Cl = 17,
    // ... all elements
    
    // Utility methods
    constexpr double atomicMass() const;
    constexpr double covalentRadius() const;
    constexpr std::string_view symbol() const;
};

// Replace: std::vector<int> m_atoms
// With:    std::vector<ElementType> m_atoms
```

### Test Validation
- `test_element_type_safety()` should work with new enum
- All existing element-based operations must work
- Performance should be equal (uint8_t vs int might be faster)

## Phase 5: Unified Atom Structure (LOW PRIORITY - BREAKING)

### Current Problem
- Separate containers for different atom properties
- Cache-unfriendly memory layout
- No zero-copy geometry access

### Refactoring Plan
```cpp
// Hybrid approach: Fast access + Rich metadata
class Molecule {
private:
    // Performance-critical: Separate arrays for hot paths
    mutable Geometry m_geometry;           // Zero-copy matrix access
    std::vector<ElementType> m_elements;   // Type-safe elements
    std::vector<double> m_charges;         // Partial charges
    
    // Rich metadata: Separate structure for less frequent access
    struct AtomMetadata {
        double mass = 0.0;
        int fragment_id = -1;
        std::vector<int> bonded_to;
    };
    std::vector<AtomMetadata> m_metadata;
    
public:
    // Zero-copy access (educational focus)
    const Geometry& coordinates() const { return m_geometry; }
    Geometry& coordinates() { invalidateCache(CacheType::GEOMETRY); return m_geometry; }
    
    // Type-safe element access
    ElementType element(int i) const { return m_elements[i]; }
    double charge(int i) const { return m_charges[i]; }
};
```

### Test Validation
- ALL existing tests must pass without modification
- `test_geometry_matrix_access()` should show zero-copy access
- Performance tests should show improvement

## Regression Testing Strategy

### Continuous Validation
For each refactoring phase:
1. **Before**: Run `./test_molecule` → All tests pass
2. **Implement**: Make changes according to plan
3. **After**: Run `./test_molecule` → All tests still pass
4. **Performance**: Compare benchmark results

### Performance Benchmarks
Key metrics to track:
- Distance matrix calculation time
- Fragment detection time
- Memory usage for large molecules (>1000 atoms)
- Cache hit rates for repeated operations

### API Compatibility
- All public methods must remain unchanged
- Return types and parameters must be identical
- Existing user code should compile without modification

## Timeline Estimate
- **Phase 1** (XYZ Parser): 1-2 days
- **Phase 2** (Cache System): 2-3 days  
- **Phase 3** (Fragments): 1-2 days
- **Phase 4** (Type Safety): 2-3 days
- **Phase 5** (Unified Structure): 3-5 days
- **Total**: ~2 weeks with thorough testing

## Success Criteria
✅ All tests pass throughout refactoring process
✅ API remains backward compatible
✅ Performance improves or remains equal
✅ Code quality and maintainability improved
✅ Educational clarity enhanced (algorithm visibility)

---
*This roadmap ensures systematic, test-driven refactoring with minimal risk of regression.*