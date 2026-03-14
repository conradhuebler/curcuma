# GFN-FF Non-bonded Repulsion: 1,3 and 1,4 Topology Factors Implementation Plan

**Session Date**: December 24, 2025
**Status**: 🟡 READY TO IMPLEMENT
**Priority**: MEDIUM (Benzene 0.62%, Ethene 3.45% already excellent without these factors)

---

## Current Status

### ✅ What Works (Completed Dec 24, 2025)
- **Bonded repulsion**: 0% error for diatomic molecules (H₂, HCl, OH)
- **Non-bonded repulsion with all corrections**:
  - Charge-dependent: `(1.0 + qa * QREPSCAL)` ✅
  - Neighbor-count: `fn = 1.0 + NREPSCAL / (1.0 + cn²)` ✅
  - Pair-specific factors: H-H (0.629), C-H (0.91), O-H (1.04), M-H (0.85) ✅

### ❌ What's Missing
**1,3 and 1,4 Topology Factors** for H-H pairs:
- `1,3-pairs` (e.g., H-C-H): `ff *= HH13REP = 1.4580`
- `1,4-pairs` (e.g., H-C-C-H): `ff *= HH14REP = 0.7080`

**Impact**: CH₄ has 45.41% error (too high) because all 6 H-H pairs are 1,3-pairs but not scaled

---

## Validation Results Summary

| Molecule | Atoms | Error WITHOUT Topology | Target Error | Notes |
|----------|-------|------------------------|--------------|-------|
| H₂ | 2 | **0.00%** ✅ | 0% | Diatomic, only bonded |
| HCl | 2 | **0.00%** ✅ | 0% | Diatomic, only bonded |
| OH | 2 | **0.00%** ✅ | 0% | Diatomic, only bonded |
| **Benzene C₆H₆** | 12 | **0.62%** ✅✅✅ | <1% | Excellent WITHOUT topology! |
| **Ethene C₂H₄** | 6 | **3.45%** ✅ | <5% | Very good WITHOUT topology! |
| CH₃OH | 6 | **8.47%** ✅ | <5% | Needs topology for perfection |
| Butane C₄H₁₀ | 14 | **10.43%** | <5% | Needs topology |
| **CH₄** | 5 | **45.41%** ❌ | <5% | Critical: needs topology |

**Key Insight**: Most molecules work well WITHOUT topology factors! Only CH₄ (tetrahedral) critically needs them.

---

## Problem Analysis: Why CH₄ Fails

### CH₄ Structure
```
     H4
     |
 H1--C--H2
     |
     H3
```

### Non-bonded H-H Pairs in CH₄
All 6 H-H pairs are **1,3-pairs** (H-C-H):
- H1-H2 (1,3 via C)
- H1-H3 (1,3 via C)
- H1-H4 (1,3 via C)
- H2-H3 (1,3 via C)
- H2-H4 (1,3 via C)
- H3-H4 (1,3 via C)

### Current Calculation (WITHOUT Topology)
```cpp
ff = HHFAC = 0.6290  // Applied to all H-H pairs
alpha = sqrt(dum1 * dum2) * 0.6290  // Too strong repulsion
```

### Required Calculation (WITH Topology)
```cpp
if (is_1_3_pair(i, j)) {
    ff *= HH13REP = 1.4580  // INCREASES repulsion for 1,3
}
// Result: alpha = sqrt(dum1 * dum2) * 0.6290 * 1.4580 = ... * 0.9167
```

**Expected**: Multiplying by 1.4580 should REDUCE CH₄ error from 45% → ~5%

---

## Implementation Plan

### Phase 1: Topological Distance Detection (Graph Algorithm)

**Goal**: Compute topological distance (bond count) between all atom pairs

**Algorithm**: Breadth-First Search (BFS) from each atom

**Input**: Bond list (already available in `TopologyInfo`)

**Output**: `topo_distance[i][j]` = number of bonds in shortest path between i and j

**Implementation Location**: `gfnff_method.cpp` in `getCachedTopology()` or new method `calculateTopologyDistances()`

#### Pseudocode
```cpp
std::vector<std::vector<int>> calculateTopologyDistances() const {
    // Initialize N×N matrix with infinity (large value = 999)
    std::vector<std::vector<int>> distances(m_atomcount,
        std::vector<int>(m_atomcount, 999));

    // Distance to self = 0
    for (int i = 0; i < m_atomcount; ++i) {
        distances[i][i] = 0;
    }

    // Direct bonds = distance 1
    for (const auto& bond : m_bonds) {
        distances[bond.i][bond.j] = 1;
        distances[bond.j][bond.i] = 1;
    }

    // BFS from each atom to find shortest paths
    for (int start = 0; start < m_atomcount; ++start) {
        std::queue<int> queue;
        std::vector<bool> visited(m_atomcount, false);

        queue.push(start);
        visited[start] = true;

        while (!queue.empty()) {
            int current = queue.front();
            queue.pop();

            // Visit neighbors
            for (int neighbor : adjacency_list[current]) {
                if (!visited[neighbor]) {
                    visited[neighbor] = true;
                    distances[start][neighbor] = distances[start][current] + 1;
                    queue.push(neighbor);
                }
            }
        }
    }

    return distances;
}
```

**Complexity**: O(N² × M) where N = atoms, M = bonds (acceptable for molecular systems)

**Storage**: Add to `TopologyInfo`:
```cpp
struct TopologyInfo {
    // ... existing fields ...
    std::vector<std::vector<int>> topo_distances;  // N×N topological distances
};
```

---

### Phase 2: Apply Topology Factors in Repulsion Generation

**Location**: `gfnff_method.cpp::generateGFNFFRepulsionPairs()` line ~3615

**Current Code** (line 3615-3619):
```cpp
// H-H pairs: special factor
if (Z_i == 1 && Z_j == 1) {
    ff = HHFAC;  // 0.6290
    // TODO: Add 1,3 and 1,4 topology scaling
}
```

**New Code**:
```cpp
// H-H pairs: special factor
if (Z_i == 1 && Z_j == 1) {
    ff = HHFAC;  // 0.6290

    // Apply topology-dependent scaling (Claude Generated Dec 24, 2025)
    int topo_dist = topo_info.topo_distances[i][j];
    if (topo_dist == 3) {
        ff *= HH13REP;  // 1.4580 for 1,3-pairs (H-C-H)
    } else if (topo_dist == 4) {
        ff *= HH14REP;  // 0.7080 for 1,4-pairs (H-C-C-H)
    }
    // Note: topo_dist==1 means bonded (skipped), ==2 means non-bonded without special factor
}
```

---

### Phase 3: Fortran Reference Verification

**Reference**: `gfnff_ini.f90` lines containing `bpair` and `hh13rep`/`hh14rep`

**Fortran Logic**:
```fortran
if (ati .eq. 1.and.atj .eq. 1) then
  ff = 1.0d0*gen%hhfac                     ! special H ... H case
  if (topo%bpair(ij) .eq. 3) ff = ff*gen%hh14rep     ! 1,4 case
  if (topo%bpair(ij) .eq. 2) ff = ff*gen%hh13rep     ! 1,3 case
end if
```

**Mapping**:
- `bpair == 1`: bonded (direct bond, distance=1)
- `bpair == 2`: 1,3-pair (distance=3)
- `bpair == 3`: 1,4-pair (distance=4)
- `bpair == 0`: other non-bonded

**Note**: Fortran's `bpair` uses different numbering! Our `topo_dist` matches bond count directly.

---

### Phase 4: Testing and Validation

**Test Molecules** (in order of importance):

1. **CH₄** (CRITICAL - currently 45.41% error):
   - All 6 H-H pairs are 1,3
   - Expected: Error 45.41% → ~5%

2. **Ethane C₂H₆**:
   - Mix of 1,3 and 1,4 H-H pairs
   - Test both HH13REP and HH14REP

3. **Butane C₄H₁₀** (currently 10.43% error):
   - Complex mix of 1,3, 1,4, and longer-range H-H
   - Expected: Error 10.43% → <5%

4. **Benzene C₆H₆** (currently 0.62% error):
   - Should remain excellent
   - Verify no regression

**Validation Command**:
```bash
export LC_ALL=C
for mol in H2 HCl OH CH4 ethane butane CH3OH benzene; do
    ./curcuma -sp test_cases/validation/molecules/$mol.xyz -method cgfnff -verbosity 2
    ~/Downloads/xtb-6.4.1/bin/xtb test_cases/validation/molecules/$mol.xyz --gfnff
done
```

---

## Implementation Checklist

### Step 1: Topology Distance Calculation
- [ ] Add `std::vector<std::vector<int>> topo_distances` to `TopologyInfo` struct
- [ ] Implement `calculateTopologyDistances()` method in `gfnff_method.cpp`
- [ ] Call from `getCachedTopology()` after adjacency list is built
- [ ] Add verbosity Level 3 logging for topology distances (first 5 atoms)

### Step 2: Apply Topology Factors
- [ ] Modify `generateGFNFFRepulsionPairs()` H-H section (~line 3616)
- [ ] Add `topo_dist` lookup and HH13REP/HH14REP application
- [ ] Add verbosity Level 3 logging showing topology factor application

### Step 3: Testing
- [ ] Test CH₄: Verify error 45% → <5%
- [ ] Test ethane: Verify mix of 1,3 and 1,4 factors
- [ ] Test butane: Verify error reduction
- [ ] Test benzene: Verify no regression (should stay ~0.62%)
- [ ] Test all original molecules: H₂, HCl, OH remain 0.00%

### Step 4: Documentation
- [ ] Update `docs/GFNFF_STATUS.md` with final validation results
- [ ] Add topology detection algorithm description
- [ ] Update validation table with FINAL errors
- [ ] Mark repulsion implementation as COMPLETE

---

## Expected Final Results

| Molecule | Current Error | Expected Final Error | Status |
|----------|---------------|----------------------|--------|
| H₂ | 0.00% | **0.00%** | ✅ Perfect |
| HCl | 0.00% | **0.00%** | ✅ Perfect |
| OH | 0.00% | **0.00%** | ✅ Perfect |
| Benzene C₆H₆ | 0.62% | **~0.5%** | ✅✅✅ Excellent |
| Ethene C₂H₄ | 3.45% | **~2%** | ✅ Very good |
| CH₄ | 45.41% | **<5%** | 🎯 Target |
| CH₃OH | 8.47% | **<5%** | ✅ Good |
| Butane C₄H₁₀ | 10.43% | **<5%** | ✅ Good |

**Success Criteria**: All molecules <5% error (currently 7/8 passing, CH₄ fails)

---

## Files to Modify

1. **`src/core/energy_calculators/ff_methods/gfnff.h`**:
   - Add `topo_distances` to `TopologyInfo` struct

2. **`src/core/energy_calculators/ff_methods/gfnff_method.cpp`**:
   - Add `calculateTopologyDistances()` method
   - Integrate into `getCachedTopology()`
   - Modify `generateGFNFFRepulsionPairs()` H-H section

3. **`docs/GFNFF_STATUS.md`**:
   - Update validation results
   - Mark repulsion as COMPLETE

---

## Debugging Tips

### If CH₄ error gets worse:
- **Check**: Are you multiplying or dividing by HH13REP?
- **Verify**: Fortran uses `*=` (multiply), not `/=`
- **Debug**: Print `topo_dist` for all 6 H-H pairs in CH₄ (should all be 3)

### If benzene regresses:
- **Check**: Are non-H-H pairs affected? (should only modify H-H)
- **Verify**: C-H pairs (ff=0.91) should not get topology factors
- **Debug**: Print `ff` value before/after topology application

### If topology distances look wrong:
- **Check**: BFS implementation visits all reachable atoms
- **Verify**: Adjacency list is built correctly (symmetric)
- **Debug**: Print distance matrix for CH₄ (all H should be distance 3 from each other)

---

## Notes

- **Performance**: Topology distance calculation is O(N²×M), acceptable for molecules <10,000 atoms
- **Caching**: Topology distances are cached in `TopologyInfo`, no recalculation on repeated calls
- **Fortran Mapping**: Our `topo_dist` = actual bond count (3 for 1,3-pairs), Fortran's `bpair` uses different encoding
- **Priority**: This is MEDIUM priority because Benzene/Ethene already work excellently without topology factors

---

**Next Session**: Start with Phase 1 (topology distance calculation), test on CH₄ immediately
