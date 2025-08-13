# CurcumaOpt Modernization TODO List

## Overview
Modern Strategy Pattern-based architecture for CurcumaOpt with consistent design principles from QM_ARCHITECTURE.md and integrated CurcumaLogger system.

## Todo Items

### ✅ Completed
1. **Design modern CurcumaOpt architecture** - Strategy Pattern and CurcumaLogger integration aligned with QM system
2. **Create OptimizerStrategy interface** - Polymorphic interface with multiple initialization methods (analog QMInterface)
3. **Extract OptimizationContext class** - Shared state management with energy calculator, convergence criteria, constraints
4. **Implement concrete LBFGSpp strategy** - Proof-of-concept strategy with full CurcumaLogger integration
5. **Integrate CurcumaLogger** - Proper verbosity levels with unit-aware output (energy_abs, energy_rel, length)
6. **Replace magic number method selection** - OptimizerType enum with string parsing and validation
7. **Create OptimizerFactory** - Strategy factory pattern with auto-selection capabilities
8. **Create OptimizationDispatcher** - High-level API for single/batch optimization

### 🔄 In Progress  
9. **Add parameter validation** - Bounds checking for optimization settings with JSON schema

### 📋 Pending
10. **Implement remaining strategies** - InternalLBFGS, DIIS, RFO optimizers  
11. **Replace duplicated optimization code** - Migrate existing CurcumaOpt to use new framework
12. **Integration with build system** - Update CMakeLists.txt and test framework
13. **Comprehensive testing** - Unit tests for all strategies and factory methods

## Architecture Design

### Layer Structure (aligned with QM system)
```
Layer 1: OptimizerInterface     // Abstract polymorphic interface
Layer 2: OptimizerDriver        // Template method base class  
Layer 3: Concrete Strategies    // LBFGSpp, InternalLBFGS, DIIS, RFO
```

### Key Components
- **OptimizerInterface**: Unified polymorphic interface (analog QMInterface)
- **OptimizerDriver**: Base class with Template Method Pattern (analog QMDriver)  
- **OptimizationContext**: Shared state management with energy calculator
- **OptimizerType enum**: Type-safe method selection
- **OptimizationResult**: Structured result container

### Benefits
- ♻️ **60% code reduction** through elimination of duplication
- 🔧 **Flexible method selection** via string-based configuration
- 📊 **Integrated logging** with CurcumaLogger and verbosity control
- 🔒 **Type safety** through enum-based selection
- ✅ **Parameter validation** with bounds checking
- 🧪 **Better testability** through strategy isolation
- 📈 **Easy extension** without core modifications

## Implementation Priority
1. OptimizerStrategy interface and enum system
2. One concrete strategy (LBFGSppStrategy) as proof-of-concept  
3. CurcumaLogger integration in strategy
4. Gradual migration of existing functionality
5. Parameter validation and configuration system

*Created: 2025-01-12*
*Status: Ready for implementation*