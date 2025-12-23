# Trajectory Analysis Framework Consolidation Plan

## Overview

This document outlines the plan to consolidate and refactor fragmented trajectory analysis code across 4+ modules (`analysis.cpp`, `trajectoryanalysis.cpp`, `rmsdtraj.cpp`, and geometry commands in `main.cpp`) into a unified framework.

**Status**: ✅ Plan documented (December 2025)
**Implementation**: Pending (6 iterations, ~18 hours estimated)
**Priority**: MEDIUM (maintenance + extensibility)
**ROI**: Very high - Eliminates long-term maintenance burden, improves consistency

---

## Problem Statement

### Current State: Fragmentation
- **4 different output implementations**: `outputHumanToStream()`, `outputCSVToStream()`, `exportTimeSeries()`, inline code
- **3 different statistics engines**: `TrajectoryStatistics`, `TimeSeriesStats`, `Tools::mean/stdev`
- **3 different progress reporters**: 10-frame steps, 100-frame steps, 5% percentage tracking
- **~800 lines of duplicate code** across modules

### Impact
- **High maintenance burden**: Changes to output format require updating 4+ locations
- **Inconsistent user experience**: Same data shown differently across modules
- **Poor testability**: Output logic only tested in specific contexts
- **Difficult to extend**: New output formats (HDF5, Parquet) require modifying multiple modules

---

## Solution Architecture

### 3 New Utility Modules

#### 1. TrajectoryWriter (`src/tools/trajectory_writer.h/cpp`)
Unified output system for all trajectory data formats.

**Formats Supported**:
- Human-readable tables (fixed-width, with cumulative statistics)
- CSV (row-per-frame with moving window statistics)
- JSON (structured metadata)
- DAT (legacy RMSD format)

**Key Features**:
- Standardized JSON input schema across all modules
- Consistent notation: `<X>` = cumulative mean, `σ(X)` = std dev, `<X>_w` = moving avg
- Simple API: `writeHumanTable()`, `writeCSV()`, `writeJSON()`, `writeDAT()`

#### 2. TrajectoryStatistics Extension (`src/capabilities/trajectory_statistics.h`)
Enhance existing `TrajectoryStatistics` class with:
- Min/Max/Median calculations
- Advanced features: autocorrelation, equilibration time, convergence checking
- Export methods: `exportStatistics()`, `exportAllStatistics()`

**Strategy**: Extend existing class rather than create new one

#### 3. ProgressTracker (`src/tools/progress_tracker.h/cpp`) [OPTIONAL]
Unified progress reporting with automatic timing information.

---

## Implementation Phases

### Phase 1: TrajectoryWriter Foundation (PRIORITY: HIGH)
- **Files**: Create `src/tools/trajectory_writer.h/cpp`
- **Effort**: ~6 hours
- **Impact**: Eliminates 80% of output-code duplication

**Steps**:
1. Extract `outputHumanToStream()` from `analysis.cpp`
2. Extract `outputCSVToStream()` from `analysis.cpp`
3. Implement `writeJSON()` wrapper
4. Implement `writeDAT()` for RMSD compatibility
5. Add unit tests

### Phase 2: Migration - Analysis Module (PRIORITY: HIGH)
- **Files**: Update `src/capabilities/analysis.cpp`
- **Effort**: ~2 hours
- **Changes**:
  - Replace private `outputHumanToStream()`/`outputCSVToStream()` calls
  - Remove private output methods (~184 lines)

### Phase 3: TrajectoryStatistics Extension (PRIORITY: MEDIUM)
- **Files**: Extend `src/capabilities/trajectory_statistics.h/cpp`
- **Effort**: ~3 hours
- **Changes**:
  - Add `getMin()`, `getMax()`, `getMedian()`
  - Add advanced stats (optional: autocorrelation, equilibration)
  - Add export methods for Writer integration

### Phase 4: Migration - Analysis + RMSD (PRIORITY: MEDIUM)
- **Files**: Update `src/capabilities/trajectoryanalysis.cpp`, `src/capabilities/rmsdtraj.cpp`
- **Effort**: ~4 hours
- **Changes**:
  - Replace custom CSV export with TrajectoryWriter
  - Replace `TimeSeriesStats` and `Tools::mean/stdev` with TrajectoryStatistics
  - Remove custom output code (~60 lines total)

### Phase 5: Cleanup - Geometry Commands (PRIORITY: LOW)
- **Files**: Update `src/main.cpp` geometry handlers
- **Effort**: ~1 hour
- **Changes**: Replace inline output code with TrajectoryWriter calls

### Phase 6: ProgressTracker Integration [OPTIONAL] (PRIORITY: LOW)
- **Files**: Create `src/tools/progress_tracker.h/cpp`
- **Effort**: ~2 hours
- **Changes**: Unified progress reporting across all analysis modules

---

## Benefits

### Quantitative
- **Code reduction**: ~800 lines eliminated
- **Test coverage**: Output logic tested once instead of 4x
- **Maintenance**: Single point of change for output formats

### Qualitative
- **Consistency**: All trajectory commands use identical notation and formats
- **Extensibility**: New formats added in one location
- **Clarity**: Clear separation of concerns (Writer vs. Statistics vs. Progress)
- **Testability**: Can test output independently of analysis logic

---

## Risks & Mitigations

### Risk: Breaking existing output format
**Mitigation**: All output formats remain identical; internal refactoring only

### Risk: Performance regression
**Mitigation**: Profile before/after; TrajectoryWriter is pure output (single pass through data)

### Risk: JSON schema incompatibility
**Mitigation**: Standardized schema with backwards-compatible export methods

---

## Success Criteria

- [ ] All trajectory analysis modules use TrajectoryWriter
- [ ] TrajectoryStatistics used uniformly across codebase
- [ ] Output formats unchanged (bit-for-bit identical)
- [ ] All existing tests pass
- [ ] Code reduction > 600 lines
- [ ] New output formats can be added without touching existing modules

---

## Timeline Estimate

**Total**: 6 iterations, ~18 hours
**Recommended pace**: 1-2 iterations per week
**Critical path**: Phases 1 → 2 (must complete before Phase 4)

**Iteration Breakdown**:
- Iteration 1: TrajectoryWriter Foundation (6 hrs) - HIGH PRIORITY
- Iteration 2: Analysis.cpp Migration (2 hrs) - HIGH PRIORITY
- Iteration 3: TrajectoryStatistics Extension (3 hrs) - MEDIUM PRIORITY
- Iteration 4: Trajectoryanalysis + RMSD Migration (4 hrs) - MEDIUM PRIORITY
- Iteration 5: Geometry Cleanup (1 hr) - LOW PRIORITY
- Iteration 6: ProgressTracker [OPTIONAL] (2 hrs) - LOW PRIORITY

---

## Related Files in Codebase

- **Reference implementation**: `src/capabilities/analysis.cpp` - Study `outputHumanToStream()` and `outputCSVToStream()` patterns
- **Existing statistics**: `src/capabilities/trajectory_statistics.h` - Base class to extend with Min/Max/Median
- **Module to migrate #1**: `src/capabilities/trajectoryanalysis.cpp` - Uses custom `TimeSeriesStats` and CSV export
- **Module to migrate #2**: `src/capabilities/rmsdtraj.cpp` - Uses inline DAT output and `Tools::mean/stdev`
- **Geometry commands**: `src/main.cpp` - `executeBond()`, `executeAngle()`, `executeTorsion()` functions (already supports trajectories)

---

## Questions & Decisions

### Q: Why extend TrajectoryStatistics instead of creating TimeSeriesAnalyzer?
**A**: The existing TrajectoryStatistics class already has:
- Complete per-metric tracking with maps
- Moving average window support
- Well-tested implementation

Extending is less risky than reimplementing.

### Q: Should we deprecate TimeSeriesStats?
**A**: Yes, with a DEPRECATED comment pointing to TrajectoryStatistics

### Q: What about progress reporting in geometry commands?
**A**: Optional for Phase 6; geometry commands are fast enough without it

---

## Approval Status

- **Plan Author**: Claude (December 2025)
- **Owner**: Conrad Hübler
- **Status**: Ready for Phase 1 implementation
