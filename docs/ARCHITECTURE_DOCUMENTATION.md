# Architecture Documentation Standards

## Purpose

This document defines comprehensive documentation standards for complex architectural patterns in Curcuma, particularly factory patterns, dispatchers, and multi-step workflows. The goal is to ensure that elegant but opaque code remains maintainable and debuggable.

## When to Apply These Standards

Apply these documentation standards when implementing:
- Factory patterns (e.g., `MethodFactory`, `OptimizerFactory`)
- Dispatcher patterns (e.g., `ModernOptimizerDispatcher`)
- Complex multi-step workflows (e.g., `optimizeStructure()` chains)
- Any code where the execution path involves multiple files and indirect calls

## Documentation Templates

### 1. Architectural Decision Record (ADR) for Complex Functions

For any complex dispatcher or factory function, include a comprehensive ADR-style comment:

```cpp
/*
 * ARCHITECTURAL DECISION RECORD: [System Name] 
 * 
 * CONTEXT: [What problem does this solve?]
 * - List the challenges this architecture addresses
 * - Explain why simple approaches are insufficient
 * 
 * DECISION: [What architecture was chosen?]
 * - Factory pattern with priority-based selection
 * - Dispatcher with automatic fallback mechanisms
 * - Template/Strategy pattern implementation
 * 
 * IMPLEMENTATION CHAIN: [Complete execution path]
 * 1. Entry Point: main.cpp:664 → ModernOptimizerDispatcher::optimizeStructure()
 * 2. Factory Call: optimizer_factory.cpp:45 → selectBestAlgorithm() 
 * 3. Method Creation: lbfgs_optimizer.cpp:34 → LBFGSOptimizer constructor
 * 4. Execution: lbfgs_optimizer.cpp:123 → optimize() main loop
 * 
 * RUNTIME BEHAVIOR: [How selection/dispatch works]
 * - "lbfgs" → LBFGSOptimizer class → src/capabilities/optimisation/lbfgs.cpp:optimize()
 * - "fire" → FIREOptimizer class → src/capabilities/optimisation/fire_optimizer.cpp:minimize()
 * - "auto" → Heuristic selection based on system size and constraints
 * 
 * DEBUGGING ENTRY POINTS: [Where to look when things go wrong]
 * - Set verbosity ≥2 to see algorithm selection logic
 * - Check OptimizerFactory::selectBestAlgorithm() for selection criteria
 * - Each optimizer reports convergence via CurcumaLogger at verbosity ≥1
 * - Common failures: Check energy calculator initialization in optimizer constructors
 */
```

### 2. Factory Class Documentation Standard

Every factory class must include this comprehensive header:

```cpp
/**
 * @file optimizer_factory.cpp
 * @brief Factory for creating optimization algorithm instances
 * 
 * SUPPORTED ALGORITHMS:
 * - "lbfgs": Limited-memory BFGS → src/capabilities/optimisation/lbfgs.cpp
 * - "fire": Fast Inertial Relaxation Engine → src/capabilities/optimisation/fire_optimizer.cpp
 * - "cg": Conjugate Gradient → src/capabilities/optimisation/cg_optimizer.cpp
 * - "auto": Automatic selection based on system characteristics
 * 
 * SELECTION HIERARCHY:
 * 1. Explicit user choice via controller["opt_algorithm"]
 * 2. Constraint-based selection (constrained systems prefer FIRE)
 * 3. Size-based selection (>1000 atoms prefer L-BFGS)
 * 4. Default fallback to L-BFGS
 * 
 * INTEGRATION POINTS:
 * - main.cpp:664: Primary entry point from CLI arguments
 * - src/capabilities/curcumaopt.cpp:123: Legacy integration point
 * - Python bindings: python/optimization_wrapper.cpp:45
 * 
 * ERROR HANDLING:
 * - Unknown algorithms → throws OptimizerCreationException with suggested alternatives
 * - Missing dependencies → automatic fallback to available algorithms with warning
 * - All errors logged via CurcumaLogger with file:line information for debugging
 * 
 * PERFORMANCE CONSIDERATIONS:
 * - Factory instances are lightweight - no caching needed
 * - Algorithm objects are created on-demand
 * - Thread-safe for concurrent factory calls
 */
```

### 3. Call Chain Documentation

For complex workflows, document the complete execution path:

```cpp
/*
 * CALL CHAIN DOCUMENTATION: Optimization Workflow
 * 
 * USER INPUT: curcuma -opt input.xyz -opt_algorithm lbfgs -gradient_threshold 1e-4
 * 
 * EXECUTION FLOW:
 * main.cpp:445 → CLI2Json() converts arguments to controller["opt_algorithm"] = "lbfgs"
 * main.cpp:664 → if(controller["modernopt"]) branch determines new vs legacy system
 * main.cpp:665 → ModernOptimizerDispatcher::optimizeStructure(molecule, method, calc, params)
 *   ↳ optimizer_factory.cpp:67 → OptimizerFactory::create("lbfgs", controller)
 *   ↳ optimizer_factory.cpp:89 → validateAlgorithmAvailability("lbfgs") 
 *   ↳ lbfgs_optimizer.cpp:34 → LBFGSOptimizer(controller, energy_calculator)
 *   ↳ lbfgs_optimizer.cpp:45 → validateParameters(controller) checks convergence criteria
 * main.cpp:666 → optimizer->optimize(molecule) starts optimization loop
 *   ↳ lbfgs_optimizer.cpp:123 → optimize() main iteration loop
 *   ↳ lbfgs_optimizer.cpp:145 → computeSearchDirection() using L-BFGS history
 *   ↳ lbfgs_optimizer.cpp:189 → lineSearch() calls EnergyCalculator repeatedly
 *   ↳ src/core/energycalculator.cpp:324 → CalculateEnergy(gradient=true) for each step
 *   ↳ lbfgs_optimizer.cpp:234 → checkConvergence() evaluates stopping criteria
 * main.cpp:667 → result.converged check and final structure output
 * 
 * RESULT FLOW:
 * LBFGSOptimizer::optimize() → OptimizationResult{converged, energy, iterations, time}
 * main.cpp:668 → Success logging via CurcumaLogger
 * main.cpp:674 → Write optimized structure to .opt.xyz file
 * 
 * ERROR PATHS:
 * Any step failure → OptimizerException → main.cpp:670 catch block → error logging
 * Energy calculation failure → EnergyCalculator error handling → optimizer abort
 * Convergence failure → Warning logged, partial result returned
 */
```

### 4. Debugging-Oriented Comments

Include specific comments for common debugging scenarios:

```cpp
// DEBUGGING NOTE: If optimization fails silently, check these locations in order:
// 1. EnergyCalculator initialization: Set controller["verbosity"] ≥ 2 to see method creation
// 2. Gradient validation: Enable lbfgs_optimizer.cpp:156 validateGradient() checks
// 3. Line search failures: Check lbfgs_optimizer.cpp:189 lineSearch() return codes
// 4. Convergence criteria: Verify controller["gradient_threshold"] is reasonable (1e-4 to 1e-6)
// 5. Step size issues: Monitor lbfgs_optimizer.cpp:203 step_size scaling
// TROUBLESHOOTING GUIDE: See docs/optimization/DEBUGGING.md for step-by-step diagnosis
```

### 5. Cross-Reference Documentation

Establish bidirectional references between related architectures:

```cpp
// RELATED ARCHITECTURE: This factory follows the same pattern as EnergyCalculator method creation
// PRIMARY REFERENCE: src/core/energy_calculators/method_factory.cpp:67 for implementation pattern
// PARALLEL USAGE: main.cpp:664 (optimization) mirrors main.cpp:445 (energy method selection)
// SIMILAR PATTERNS: See MethodFactory::create() for priority-based method resolution
// 
// INTEGRATION NOTE: Changes to this factory should be coordinated with:
// - EnergyCalculator method selection logic (src/core/energycalculator.cpp:109)
// - Legacy optimization system (src/capabilities/curcumaopt.cpp:123)
// - Python bindings optimization wrapper (python/optimization_wrapper.cpp:45)
```

## Implementation Examples

### Example 1: ModernOptimizerDispatcher Documentation

This is how `main.cpp:664` should be documented:

```cpp
/*
 * ARCHITECTURAL DECISION RECORD: Modern Optimization Dispatcher
 * 
 * CONTEXT: Multiple optimization algorithms with different capabilities and performance
 * - L-BFGS: Best for unconstrained problems, memory efficient for large systems
 * - FIRE: Superior for constrained problems and surface relaxation  
 * - Conjugate Gradient: Good balance for medium-sized systems
 * - Legacy system: Backward compatibility for existing workflows
 * 
 * DECISION: Factory pattern with automatic algorithm selection
 * - Priority-based selection: user choice > system characteristics > defaults
 * - Graceful fallback: if preferred algorithm unavailable, select next best
 * - Unified interface: same OptimizationResult for all algorithms
 * 
 * IMPLEMENTATION CHAIN:
 * 1. main.cpp:664 → ModernOptimizerDispatcher::optimizeStructure() 
 * 2. src/capabilities/optimisation/optimizer_factory.cpp:45 → selectBestAlgorithm()
 * 3. src/capabilities/optimisation/lbfgs_optimizer.cpp:34 → LBFGSOptimizer()
 * 4. src/capabilities/optimisation/lbfgs_optimizer.cpp:123 → optimize()
 * 
 * RUNTIME BEHAVIOR:
 * - controller["opt_algorithm"] = "lbfgs" → Direct LBFGSOptimizer creation
 * - controller["opt_algorithm"] = "auto" → Heuristic: >1000 atoms→LBFGS, constraints→FIRE
 * - Invalid algorithm → Warning + fallback to L-BFGS with user notification
 * 
 * DEBUGGING ENTRY POINTS:
 * - Verbosity ≥ 2: Shows algorithm selection reasoning and parameters
 * - OptimizerFactory::selectBestAlgorithm() logs selection criteria
 * - Each optimizer reports detailed convergence info at verbosity ≥ 1
 */
if (controller["modernopt"].get<bool>()) {
    // Use modern optimization system with automatic algorithm selection
    auto result = ModernOptimization::ModernOptimizerDispatcher::optimizeStructure(
        &molecule, method, &energy_calc, controller["opt"]);
        
    if (result.success) {
        CurcumaLogger::success("Optimization completed successfully!");
        // ... result handling
    }
} else {
    // FALLBACK: Legacy optimization system for backward compatibility
    // See: src/capabilities/curcumaopt.cpp for implementation
    CurcumaOpt opt(controller, false);
    opt.start();
}
```

### Example 2: Factory Class Implementation

```cpp
/**
 * @file method_factory.cpp  
 * @brief Factory for creating computational method instances
 * 
 * SUPPORTED METHODS:
 * - "gfn2": GFN2-xTB → TBLite > Ulysses > XTB (priority order)
 * - "gfn1": GFN1-xTB → TBLite > XTB > Ulysses
 * - "uff": Universal Force Field → src/core/energy_calculators/ff_methods/forcefield_method.cpp
 * - "eht": Extended Hückel Theory → src/core/energy_calculators/qm_methods/eht_method.cpp
 * 
 * PRIORITY RESOLUTION:
 * 1. Check TBLite availability and compatibility
 * 2. Fallback to Ulysses interface if TBLite unavailable  
 * 3. Final fallback to XTB interface
 * 4. Throw MethodCreationException if all interfaces fail
 * 
 * INTEGRATION POINTS:
 * - src/core/energycalculator.cpp:109: Primary usage in EnergyCalculator constructor
 * - src/capabilities/hessian.cpp:113: Silent method creation for numerical derivatives
 * - src/capabilities/curcumaopt.cpp:145: Method creation for optimization
 * 
 * ERROR HANDLING:
 * - Missing libraries → automatic fallback with CurcumaLogger::warn()
 * - Invalid method names → MethodCreationException with suggested alternatives
 * - Interface initialization failure → detailed error with troubleshooting steps
 */
std::unique_ptr<ComputationalMethod> MethodFactory::create(
    const std::string& method_name, const json& config) {
    
    // DEBUGGING NOTE: Add verbosity ≥ 3 here to see method resolution process
    CurcumaLogger::info("Creating computational method: " + method_name);
    
    // Method resolution with priority-based fallback...
}
```

## Mandatory Standards for New Code

### 1. All Factory Classes Must Include:
- Complete ADR documentation explaining the architecture decision
- Exhaustive list of supported types/algorithms with file paths
- Clear integration points with file:line references  
- Comprehensive error handling documentation

### 2. All Dispatcher Functions Must Include:
- Complete call chain documentation from user input to result
- Runtime behavior explanation for each supported path
- Explicit debugging entry points with verbosity recommendations
- Cross-references to related architectural patterns

### 3. Complex Workflows Must Include:
- Step-by-step execution flow with file:line references
- Error path documentation with recovery mechanisms
- Performance considerations and optimization notes
- Integration testing and validation procedures

## Enforcement

These standards are **mandatory** for all new complex architectural implementations. Code reviews must verify:
- ADR documentation completeness
- Call chain accuracy and file:line correctness
- Debugging information usefulness
- Cross-reference bidirectionality

Existing complex code should be upgraded to these standards during maintenance cycles or when debugging issues arise.