# optimisation notes (former CLAUDE.md, as of 2026-10-01)
Source: the complete previous `src/capabilities/optimisation/CLAUDE.md` (179 lines, commit dad35624), moved unchanged.
Date: 2026-10-01. Statements are as old as their date; the maintained short version is [src/capabilities/optimisation/CLAUDE.md](../../src/capabilities/optimisation/CLAUDE.md).

Statements below that the code contradicts (checked 2026-10-01):
- Line numbers: `LBFGSStep()` is at `lbfgs.cpp:398`, `DIISStep()` 276 (with `diisExtrapolation()` 227), `RFOStep()` 511, `sr1Update()` 615.
- DIIS "SVD-based coefficient solution": the coefficients are solved with `colPivHouseholderQr`; `JacobiSVD` only gives the condition number (`lbfgs.cpp:248-266`).
- SR1 formula "(s·(y-Hs))/(s·(y-Hs)) * outer_product(y-Hs, y-Hs)": the code computes H + (y-Hs)(y-Hs)^T / ((y-Hs)·s) (`lbfgs.cpp:619-625`).
- "RFO-specific line search with curvature conditions": `lineSearchRFO()` exists but has no caller, nor has `line_search_backtracking()`; RFO scales its step 1 / 0.5 / 0.1 or rejects it (`lbfgs.cpp:582-603`).
- "Modern Optimizer System (`modern_optimizer_simple.h/.cpp`)", its auto-selection, `optimizeWithLBFGSpp()`, `optimizeWithInternal()`, legacy bridge: the file is in no CMake source list and no file includes it. The live dispatcher is `OptimizerFactory`/`OptimizationDispatcher` in `../optimizer_factory.*`.
- "LBFGSpp start geometry evaluated once ... `LBFGSppObjectiveFunction::prime()`": true, but it lives in `../lbfgspp_optimizer.h:64` / `.cpp:176`, not in `modern_optimizer_simple`.
- CLI block: `-optimizer internal` is not a recognised name (falls back to lbfgspp with a warning, `optimizer_interface.cpp:98`); `auto` always selects LBFGSpp (`optimizer_factory.cpp:62-66`), not "L-BFGS above 50 atoms"; without `-optimizer` the default is `auto`, there is no fallback to legacy CurcumaOpt (`main.cpp:2004`).
- "Automatic Fallback: Legacy system für unbekannte Optimizer": unknown names fall back to LBFGSpp.
- "Educational Output" example: the code prints "Using native Curcuma optimizer: <name>" for every name except `auto` (`main.cpp:2006-2008`), lbfgspp included.
- "Phase 5: Strong-Wolfe line search ... optional via `m_line_search_method`" and "Verbosity 2: per-step alpha ...": implemented in `lbfgs.cpp`, but `NativeOptimizerAdapter` never calls `setConfig`/`setLineSearch`/`setRFOSolver` and sets verbosity 0 (`native_optimizer_adapters.cpp:35`); only the unbuilt `modern_optimizer_simple.cpp` called `setConfig`.
- "EIGEN_USE_LAPACKE: Enabled per-file": only when CMake found LAPACKE (`HAVE_LAPACKE`, `CMakeLists.txt:1478-1487`).
- File organization tree omits `rf_solver.h/.cpp`.

---

# CLAUDE.md - Optimization Algorithms Directory

## Overview

This directory contains the **core optimization algorithms** for molecular geometry optimization. Focus is on **educational clarity** and **scientific transparency** rather than software engineering abstractions.

## 🧪 **SCIENTIFIC CONTENT** (Core Algorithms)

### Native LBFGS Implementation (`lbfgs.h/.cpp`)
**Mathematical algorithms - the actual computational chemistry:**

#### L-BFGS Algorithm (`LBFGSStep()` lines 263-320)
- **Two-loop recursion** for limited-memory BFGS
- **Scientific core**: Quasi-Newton method using gradient history
- **Literature**: Nocedal & Wright "Numerical Optimization", Chapter 7
- **Key equations**: H₀⁻¹ approximation, α/β coefficients, search direction computation

#### DIIS Algorithm (`DIISStep()` lines 245-261) 
- **Direct Inversion of Iterative Subspace** - acceleration method
- **Scientific core**: Linear combination of previous solutions to minimize residual
- **Literature**: Pulay, P. (1980) Chem. Phys. Lett. 73, 393
- **Key equations**: B-matrix construction, SVD-based coefficient solution

#### RFO Algorithm (`RFOStep()` lines 347-391)
- **Rational Function Optimization** - eigenvector following
- **Scientific core**: Newton-Raphson in eigenvector space of Hessian
- **Literature**: Banerjee et al. (1985) J. Phys. Chem. 89, 52
- **Key equations**: Eigenvalue shifting, mass-weighted coordinates, τ coefficients

#### SR1 Hessian Update (`sr1Update()` lines 390-400)
- **Symmetric Rank-1 update** formula
- **Scientific core**: (s·(y-Hs))/(s·(y-Hs)) * outer_product(y-Hs, y-Hs)
- **Literature**: Broyden (1970), Dennis & Moré (1977)

### Line Search Methods
- **Backtracking line search** with Armijo condition
- **RFO-specific line search** with curvature conditions
- **Step size control** and convergence monitoring

## 🏗️ **ORGANIZATIONAL CODE** (Software Engineering)

### Modern Optimizer System (`modern_optimizer_simple.h/.cpp`)
**Pure software engineering - no computational chemistry:**

#### Strategy Pattern Dispatcher
- **OptimizerType enum**: Type-safe method selection (no magic numbers)
- **parseOptimizerType()**: String-to-enum conversion with error handling
- **Auto-selection logic**: Size-based algorithm choice (>200 atoms → LBFGSpp, etc.)

#### Integration Layer  
- **optimizeWithLBFGSpp()**: Wrapper for external LBFGSpp library
- **LBFGSpp start geometry evaluated once** (Sep 29, 2026, 🤖): `InitializeSingleSteps()` reuses the driver's initial energy/gradient via `LBFGSppObjectiveFunction::prime()` (one fewer SCF per optimisation)
- **optimizeWithInternal()**: Wrapper for internal GPTLBFGS
- **Legacy compatibility**: Bridge to existing CurcumaOpt interface

#### Logging & Results
- **CurcumaLogger integration**: Verbosity-controlled output with colors
- **SimpleOptimizationResult**: Standardized result structure
- **Unit-aware output**: Energy in kJ/mol, distances in Å, time in seconds

## File Organization

```
optimisation/
├── lbfgs.h/.cpp           # 🧪 SCIENCE: Core optimization algorithms (L-BFGS, DIIS, RFO)
├── modern_optimizer_simple.h/.cpp  # 🏗️ ORGANIZATION: Strategy pattern dispatcher
└── CLAUDE.md              # This documentation
```

## Key Principles

### Educational-First Design
- **Algorithm visibility**: Core mathematical methods directly accessible
- **Minimal abstraction**: No unnecessary templates or inheritance hiding the science
- **Clear separation**: Scientific algorithms vs. software organization clearly marked
- **Literature references**: Papers and equations documented in code comments

### Scientific Transparency
- **Direct implementation**: L-BFGS two-loop recursion visible and understandable
- **Mathematical naming**: Variables reflect scientific meaning (gradient, hessian, eigenvalues)
- **Parameter meaning**: Physical/mathematical significance clear from context
- **Debugging access**: Easy inspection of intermediate algorithmic steps

## Instructions Block

**PRESERVED - DO NOT EDIT BY CLAUDE**

*Future tasks and visions defined by operator/programmer*

## 🚀 **EINHEITLICHES OPTIMIZER-SYSTEM** (Claude Generated)

### **Kommandozeilen-Interface:**
```bash
# 🧪 CURCUMA'S EIGENE OPTIMIERER:
./curcuma -opt input.xyz -optimizer native_lbfgs   # Native L-BFGS two-loop recursion
./curcuma -opt input.xyz -optimizer diis           # Native DIIS acceleration  
./curcuma -opt input.xyz -optimizer rfo            # Native RFO eigenvector following
./curcuma -opt input.xyz -optimizer auto           # Auto-selection (>50 Atome → L-BFGS)

# 🏗️ EXTERNE BIBLIOTHEKEN:
./curcuma -opt input.xyz -optimizer lbfgspp        # External LBFGSpp library
./curcuma -opt input.xyz -optimizer internal       # Legacy internal LBFGS

# Standard (ohne -optimizer): Fallback zu legacy system
./curcuma -opt input.xyz                           # → Legacy CurcumaOpt
```

### **Einheitliche Systematik - FESTGEHALTEN:**
- **`-method`**: **WAS** berechnet wird (Potential: uff, gfnff, eht, etc.) 
- **`-optimizer`**: **WIE** optimiert wird (native_lbfgs, diis, rfo, etc.)
- **Einheitlich für alle Module**: ConfScan, RMSD, etc. können das gleiche System nutzen

### **Educational Output:**
```
🧪 Using native Curcuma optimizer: diis
🧪 SCIENCE: Running native DIIS acceleration method  
Algorithm: Direct Inversion of Iterative Subspace
Citation: Pulay, P. (1980) Chem. Phys. Lett. 73, 393
DIIS history: 10
DIIS start iteration: 5

Iteration 5 (L-BFGS): E = -154.123456 Eh, |grad| = 1.23e-04
Iteration 10 (DIIS): E = -154.124567 Eh, |grad| = 5.67e-05
🧪 Native DIIS converged!

🎉 Optimization successful!
Method: native_diis
Iterations: 12  
Time: 2.345 seconds
Optimized structure saved to: input.opt.xyz
```

## Variable Section

### AI Implementation Status (Apr 2026)

> **All native optimizers (L-BFGS, DIIS, RFO, SR1) are 🤖 AI-generated.**
> None have been ✅ TESTED or ✅ APPROVED by the human operator.

| Method | Status | Tested on | Not tested |
|--------|--------|-----------|------------|
| Native L-BFGS | 🤖 AI-generated, ⚙️ compiles | water, ethane (UFF) | large systems, QM methods, constrained opt |
| Native DIIS | 🤖 AI-generated, ⚙️ compiles | small UFF molecules | convergence stability, near-degenerate cases |
| Native RFO | 🤖 AI-generated, ⚙️ compiles | small UFF molecules | saddle point searches, transition states |
| SR1 update | 🤖 AI-generated, ⚙️ compiles | indirectly via DIIS | standalone correctness vs. reference |

**Known gaps vs. reference implementations:**
- **L-BFGS stalls permanently** once the backtracking search is exhausted (takes the ~1e-9 step, no history reset): caffeine/GFN-FF frozen from step 33 at |g| 1.6e-3, lbfgspp converges in 38 (Sep 28, 2026, TODO.md). The driver's stall detection now ends such runs; the optimizer itself is not fixed.
- No independent numerical gradient check for any native method
- Step size control not validated against reference optimizer (XTB, Gaussian)
- DIIS extrapolation: no guarantee of convergence on difficult PES
- RFO: mode following untested

**Human production testing pending — do not use for production calculations.**

### Large-System Optimizations (Apr 2026, all 🤖 AI-generated, ⚙️ compiles, not ✅ TESTED)

**`rf_solver.h/.cpp`** (Phase 6 — shared Rational Function solver):
- `RFSolver::lanczosLowestEigenpair()` — Lanczos for lowest eigenpair of augmented RF matrix
- `RFSolver::calculateRFStep()` — dispatches Lanczos (nvar≥49) or dense fallback
- Used by ANCOpt and native RFO (replaces per-optimizer duplicate code)
- Enables iterative RF step scaling for large systems without O(N³) eigendecomp

**`lbfgs.cpp` enhancements**:
- Phase 5: Strong-Wolfe line search (`lineSearchStrongWolfe`, N&W Alg 3.5/3.6) — optional via `m_line_search_method = "strong_wolfe"`
- Phase 6: `RFOStep()` replaced O(N³) `SelfAdjointEigenSolver` with `RFSolver::calculateRFStep()`
- Verbosity 2: per-step alpha, step norm, gradient norm, history size

**EIGEN_USE_LAPACKE**: Enabled per-file in `rf_solver.cpp` and `lbfgs.cpp` (safe — no local variable `I`). Enables LAPACK D&C backend for `SelfAdjointEigenSolver` in fallback path.

### 🎯 **EINHEITLICHES DESIGN FESTGEHALTEN:**
1. **`-optimizer` Parameter**: Bestimmt Optimierungsalgorithmus (nicht `-method`)
2. **Native First**: Curcuma's eigene Implementierungen bevorzugt  
3. **Educational Output**: Literaturzitate und wissenschaftliche Parameter sichtbar
4. **Automatic Fallback**: Legacy system für unbekannte Optimizer
5. **Erweiterbar**: System kann einfach um weitere Optimizer erweitert werden

---

**Educational Focus**: This directory demonstrates computational chemistry optimization algorithms with minimal software engineering overhead, prioritizing learning and understanding over complex abstractions.
